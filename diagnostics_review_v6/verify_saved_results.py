"""Replay inputs, balanced predictive draws, MC precision and baseline warnings.

Only reads completed inference. Verification output is separate from both the
experiment and its generated report; no inference array is altered.
"""
import argparse,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from models import get_benchmark_a_lvl,get_benchmark_b_water_table
from forward_mfj import ActiveSurveyGeometry,PassiveSurveyGeometry,FullWavefieldMFJOperator
from fusion import ActivePassiveSpectralFuser,BroadbandFusionConfig
from smvi_inversion import SMVIConfig,PickingFreeSpectralLikelihood
from ablation_experiments import profile_metrics
from run_benchmark import code_hashes


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--benchmark',required=True,choices=['benchmark_a','benchmark_b'])
    parser.add_argument('--results-dir',default='results_v6');parser.add_argument('--output-dir',default='results_v6_validation');args=parser.parse_args()
    source=ROOT/args.results_dir;out=ROOT/args.output_dir;out.mkdir(parents=True,exist_ok=True)
    report=json.loads((source/'benchmark_verification_summary.json').read_text());manifest=json.loads((source/'run_manifest.json').read_text())
    assert report['status']==manifest['status']=='completed'
    assert report['source_sha256']==manifest['source_sha256']==code_hashes()
    for name,sha in report['source_sha256'].items():
        assert hashlib.sha256((source/'source_snapshot'/name).read_bytes()).hexdigest()==sha
    code=args.benchmark
    model,space=(get_benchmark_a_lvl if code=='benchmark_a' else get_benchmark_b_water_table)()
    with np.load(source/(code+'_observations.npz')) as data:obs={k:data[k] for k in data.files}
    with np.load(source/(code+'_inference.npz')) as data:res={k:data[k] for k in data.files}
    cfg=SMVIConfig(**report[code]['config'])
    op=FullWavefieldMFJOperator(ActiveSurveyGeometry(),PassiveSurveyGeometry(),obs['freqs_active'],obs['freqs_passive'],obs['c_grid'],report['quadrature_density'])
    fuser=ActivePassiveSpectralFuser(op,BroadbandFusionConfig(norm_power=1.35))
    replay=fuser.generate_fused_spectrum_for_model(model,True,np.random.default_rng(cfg.random_seed))
    for key in ['u_active','u_passive','d_fused','w_active','w_passive']:
        np.testing.assert_allclose(replay[key],obs[key],rtol=1e-10,atol=1e-20 if key.startswith('u_') else 1e-10)
    likelihood=PickingFreeSpectralLikelihood(fuser,space,obs['d_fused'],cfg)
    np.testing.assert_allclose(res['posterior_samples_phys'].mean(axis=0),res['posterior_mean_phys'])
    _,j,n,w,mean_spectrum=likelihood.evaluate_physical_model(res['posterior_mean_phys'])
    np.testing.assert_allclose(mean_spectrum,res['d_syn_posterior_mean'],atol=1e-10,rtol=1e-10)
    stored=report[code]['mean_model_metrics']
    np.testing.assert_allclose([j,1-n,w],[stored['misfit'],stored['ncc_similarity'],stored['wasserstein_w2_mps']],rtol=1e-9,atol=1e-9)
    labels=res['predictive_comp_labels'];counts=np.bincount(labels,minlength=cfg.n_components)
    np.testing.assert_array_equal(counts,res['predictive_component_counts'])
    assert np.all(counts==counts[0])
    assert cfg.predictive_draws<=len(labels)<=cfg.predictive_max_draws
    banks={k:[] for k in range(cfg.n_components)}
    for index,(z,k) in enumerate(zip(res['predictive_samples_z'],labels)):
        banks[k].append(likelihood.evaluate_physical_model(space.latent_to_physical(z))[-1])
        if index%64==0:print(code,'replayed predictive draw',index+1,flush=True)
    banks={k:np.array(images) for k,images in banks.items()};weights=res['weights']
    predictive=sum(weights[k]*banks[k].mean(axis=0) for k in banks)
    variance_mean=sum(weights[k]**2*banks[k].var(axis=0,ddof=1)/len(banks[k]) for k in banks)
    response_variance=np.maximum(sum(weights[k]*np.mean(banks[k]**2,axis=0) for k in banks)-predictive**2,0.)
    for a,b in [(predictive,res['d_posterior_predictive_mean']),
                (np.sqrt(variance_mean),res['d_posterior_predictive_mc_se']),
                (np.sqrt(response_variance),res['d_posterior_predictive_std'])]:
        np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-9)
    pm=report[code]['predictive_metrics'];j,n,w=likelihood.spectral_metrics(predictive)
    np.testing.assert_allclose([j,1-n,w],[pm['misfit'],pm['ncc_similarity'],pm['wasserstein_w2_mps']],atol=1e-9,rtol=1e-9)
    rng=np.random.default_rng(cfg.random_seed+200003)
    for h in pm['precision_history']:
        per=h['draws']//cfg.n_components;arrays={k:bank[:per] for k,bank in banks.items()}
        image=sum(weights[k]*arrays[k].mean(axis=0) for k in banks)
        half=sum(weights[k]*arrays[k][:per//2].mean(axis=0) for k in banks)
        variance=sum(weights[k]**2*arrays[k].var(axis=0,ddof=1)/per for k in banks)
        boot=[]
        for _ in range(cfg.predictive_bootstrap_repeats):
            im=sum(weights[k]*arrays[k][rng.integers(0,per,per)].mean(axis=0) for k in banks)
            boot.append(likelihood.spectral_metrics(im)[0])
        se=float(np.std(boot,ddof=1));im_se=float(np.sqrt(np.mean(variance)))
        delta=abs(likelihood.spectral_metrics(half)[0]-likelihood.spectral_metrics(image)[0])
        np.testing.assert_allclose([se,im_se,delta],[h['misfit_mc_se'],h['mean_spectrum_mc_se_rms'],h['nested_half_misfit_difference']],atol=1e-10,rtol=1e-9)
        assert h['precision_target_met']==bool(se<=cfg.predictive_mc_tolerance and im_se<=cfg.predictive_mc_tolerance and delta<=2*cfg.predictive_mc_tolerance)
    assert pm['precision_target_met']==pm['precision_history'][-1]['precision_target_met']
    baseline=report[code]['baseline_comparison']
    prior=space.latent_to_physical(np.zeros(space.dim));np.testing.assert_allclose(prior,res['prior_midpoint_phys'])
    _,jp,np_,wp,_=likelihood.evaluate_physical_model(prior)
    for key,value in [('misfit',jp),('ncc_similarity',1-np_),('wasserstein_w2_mps',wp)]:
        assert np.isclose(value,baseline['prior_midpoint'][key])
    prior_profile=profile_metrics(prior[None,:],space);actual=profile_metrics(res['posterior_samples_phys'],space)
    expected=[]
    if stored['misfit']>jp:expected.append('mean_model_fit_worse_than_prior_midpoint')
    if pm['misfit']>jp:expected.append('predictive_point_estimate_worse_than_prior_midpoint')
    if not pm['precision_target_met']:expected.append('predictive_monte_carlo_precision_target_not_met')
    if abs(pm['misfit']-jp)<=2*pm['misfit_mc_se']:expected.append('predictive_prior_difference_not_resolved_by_mc_precision')
    if actual['vs_profile_relative_rmse_pct']>prior_profile['vs_profile_relative_rmse_pct']:expected.append('posterior_profile_RMSE_worse_than_informative_prior_midpoint')
    assert expected==baseline['warnings']
    assert len(res['history_elbo'])==cfg.n_iterations
    record={'benchmark':code,'status':'saved_arrays_metrics_and_mc_precision_reproduced','actual_iterations':cfg.n_iterations,
        'balanced_predictive_draws':len(labels),'component_counts':counts.tolist(),'predictive_metrics':pm,
        'baseline_comparison':baseline,'posterior_convergence_validated':False,
        'source_sha256':code_hashes(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/(code+'_saved_verification.json')).write_text(json.dumps(record,indent=2,allow_nan=False))
    print(code,'all stored values and diagnostics reproduced',flush=True)

if __name__=='__main__':main()
