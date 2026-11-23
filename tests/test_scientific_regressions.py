"""Regression checks for mathematical validity, not expected geological outcomes."""
from types import SimpleNamespace
import ast
from pathlib import Path
import numpy as np
import pytest
from scipy.stats import multivariate_normal
from models import get_benchmark_a_lvl, get_benchmark_b_water_table
from smvi_inversion import (SMVIConfig, SteinMixtureVariationalInference,
    AffineInvariantEnsembleMCMC, compute_spectral_misfit_numba, discrete_wasserstein2)
from forward_mfj import (ActiveSurveyGeometry, PassiveSurveyGeometry, FullWavefieldMFJOperator,
                        _psv_surface_compliance_single_wk)

@pytest.mark.parametrize('getmodel',[get_benchmark_a_lvl,get_benchmark_b_water_table])
def test_conditional_transform_and_jacobian(getmodel):
    _,space=getmodel();rng=np.random.default_rng(3)
    for z in rng.normal(0,2,(8,space.dim)):
        m=space.latent_to_physical(z)
        assert np.all(m[4:8]>=1.45*m[:4])
        np.testing.assert_allclose(space.physical_to_latent(m),z,atol=1e-10)
        eps=1e-5;eye=np.eye(space.dim)*eps
        jac=np.column_stack([(space.latent_to_physical(z+e)-space.latent_to_physical(z-e))/(2*eps) for e in eye])
        np.testing.assert_allclose(space.latent_log_abs_det_jacobian(z),np.linalg.slogdet(jac)[1],rtol=1e-8)

def test_discrete_wasserstein_translation_nonuniform_grid():
    c=np.array([100.,110.,130.,200.]);p=np.array([0.,1.,0.,0.]);q=np.array([0.,0.,1.,0.])
    assert discrete_wasserstein2(p,q,c)==pytest.approx(400.)
    _,_,w=compute_spectral_misfit_numba(p[None],q[None],c,np.ones(1))
    assert w==pytest.approx(20.,rel=1e-9)
    assert discrete_wasserstein2(p,p,c)==0

def test_frequency_mask_uses_hz():
    c=np.linspace(100,1500,40);f=np.array([6.5,7.,11.,13.])
    obs=np.zeros((4,40));obs[:,2]=1;syn=obs.copy();syn[:,-1]=1
    for i in range(4):
        weights=np.zeros(4);weights[i]=1
        _,loss,_=compute_spectral_misfit_numba(obs,syn,c,weights,freqs=f)
        assert (loss==0) == (f[i]<=12)

def test_full_covariance_density_score_and_sampling_agree():
    space=SimpleNamespace(dim=2,latent_to_physical=lambda z:z)
    cfg=SMVIConfig(n_components=1)
    e=SteinMixtureVariationalInference(SimpleNamespace(param_space=space,config=cfg))
    e.means_z[:]=[.2,-.3];e.log_stds_z[:]=np.log([.4,.7]);e.offdiag_z[0,1,0]=.5
    L=e._cholesky(0);cov=L@L.T;z=np.array([.6,.1])
    lp,score=e._compute_mixture_log_prob_and_score(z)
    assert lp==pytest.approx(multivariate_normal.logpdf(z,e.means_z[0],cov))
    np.testing.assert_allclose(score,-np.linalg.solve(cov,z-e.means_z[0]))
    draws,_,_=e.sample_posterior(30000)
    np.testing.assert_allclose(np.cov(draws.T),cov,atol=.015)

def test_initialization_does_not_read_truth():
    class NoTruth:
        dim=2
        @property
        def true_vector(self):raise AssertionError('truth leakage')
    likelihood=SimpleNamespace(param_space=NoTruth(),config=SMVIConfig(n_components=3),
                               log_posterior_latent=lambda z:(-.5*z@z,0,0,0))
    engine=SteinMixtureVariationalInference(likelihood);engine._initialize_mixture_components()
    assert len(np.unique(engine.means_z,axis=0))==3

def test_mcmc_normal_has_no_clipping_atoms():
    class Normal:
        param_space=SimpleNamespace(dim=1,latent_to_physical=lambda z:z)
        def log_posterior_latent(self,z):return -.5*((z[0]-3.4)/.8)**2,0,0,0
    r=AffineInvariantEnsembleMCMC(Normal(),n_walkers=12,n_steps=2500,burn_in=500,seed=77).run()
    z=r['samples_z'][:,0]
    assert abs(z.mean()-3.4)<.08
    assert abs(z.std()-.8)<.08
    assert not np.any(z==3.5)
    assert np.mean(z>3.5)>.3

def test_geometry_preserves_actual_distances():
    p=PassiveSurveyGeometry();xy=p.station_xy
    raw=np.array([np.linalg.norm(xy[j]-xy[i]) for i in range(25) for j in range(i+1,25)])
    assert p.interstation_offsets.min()==pytest.approx(raw.min())
    assert p.pair_redundancy_weights.sum()==300
    coords=np.array([[0.,0.],[10.,0.],[0.,10.]])
    custom=PassiveSurveyGeometry(n_stations=3,coordinates=coords,use_option4_2d_nested=False)
    np.testing.assert_allclose(custom.interstation_offsets,[10.,np.sqrt(200.)])

def test_homogeneous_compliance_matches_analytic_rayleigh_function():
    h=np.array([0.]);vs=np.array([300.]);vp=np.array([600.]);rho=np.array([2000.]);q=np.array([100.])
    omega=2*np.pi*8+.1j
    for c in [180.,280.,450.,900.]:
        k=2*np.pi*8/c
        vs_c=vs[0]*(1-.5j/q[0]);vp_c=vp[0]*(1-.5j/q[0]);mu=rho[0]*vs_c**2
        pp=np.sqrt(k*k-omega**2/vp_c**2);ss=np.sqrt(k*k-omega**2/vs_c**2)
        rayleigh=(2*k*k-omega**2/vs_c**2)**2-4*k*k*pp*ss
        expected=-pp*(omega**2/vs_c**2)/(mu*rayleigh)
        got=_psv_surface_compliance_single_wk(omega,k,h,vs,vp,rho,q,q)
        np.testing.assert_allclose(got,expected,rtol=1e-8,atol=1e-20)

def small_operator(density=1.):
    return FullWavefieldMFJOperator(ActiveSurveyGeometry(n_channels=12),
        PassiveSurveyGeometry(n_stations=8,dx_station=12.,use_option4_2d_nested=False),
        np.array([8.,20.]),np.array([2.,8.]),np.linspace(100,1200,50),density)

def test_imaged_spectrum_depends_on_compliance(monkeypatch):
    import forward_mfj as fm
    m,_=get_benchmark_a_lvl();op=small_operator()
    u,_=op.synthesize_frequency_gather(m)
    before=op.compute_mfj_spectrum(u)
    def fake(freq,k,eta,*args):return np.full((len(freq),len(k)),1+2j)
    monkeypatch.setattr(fm,'_psv_surface_compliance_grid_numba',fake)
    u,_=op.synthesize_frequency_gather(m)
    after=op.compute_mfj_spectrum(u)
    assert np.max(abs(before-after))>.05

def test_noise_snr_is_measured_before_imaging():
    from ablation_experiments import add_complex_noise
    rng=np.random.default_rng(2);signal=rng.normal(size=(9,20))+1j*rng.normal(size=(9,20))
    observed,measured=add_complex_noise(signal,2.,rng)
    assert measured==pytest.approx(2.)
    assert 20*np.log10(np.linalg.norm(signal)/np.linalg.norm(observed-signal))==pytest.approx(2.)

def test_synthetic_conventional_transforms_are_actual_kernels():
    from acquisition_design import transform_spectrum
    op=small_operator();freq=op.freqs_active;r=op.r_act
    u=np.exp(2j*np.pi*freq[:,None]*r[None,:]/300.)
    spectra=[transform_spectrum(op,u,'active',kind) for kind in ['mfj','plane_wave','standing_bessel']]
    assert all(np.all(np.isfinite(s)) for s in spectra)
    assert np.max(abs(spectra[0]-spectra[2]))>.01

def test_single_gaussian_optimizer_learns_correlated_target():
    mu=np.array([.6,-.4]);cov=np.array([[.49,.175],[.175,.3125]]);inv=np.linalg.inv(cov)
    class GaussianTarget:
        param_space=SimpleNamespace(dim=2,latent_to_physical=lambda z:z)
        config=SMVIConfig(n_components=1,n_iterations=220,n_mc_samples_per_comp=32,
                          lr_mean=.035,lr_chol=.018,n_posterior_draws=2000,random_seed=73)
        eval_count=0
        def log_posterior_latent(self,z):
            dz=z-mu;return -.5*dz@inv@dz,0.,0.,0.
        def score_and_post_latent(self,z):
            self.eval_count+=1
            lp,*rest=self.log_posterior_latent(z)
            return -inv@(z-mu),lp,*rest
        def evaluate_physical_model(self,m):
            return 0.,0.,0.,0.,np.ones((2,2))
    result=SteinMixtureVariationalInference(GaussianTarget()).run(verbose=False)
    np.testing.assert_allclose(result['means_z'][0],mu,atol=.13)
    L=result['cholesky_z'][0]
    np.testing.assert_allclose(L@L.T,cov,atol=.10)

def test_layer_partition_does_not_change_homogeneous_compliance():
    omega=2*np.pi*8+.2j;k=.2
    single=(np.array([0.]),np.array([300.]),np.array([600.]),np.array([2000.]),np.array([50.]),np.array([100.]))
    layered=(np.array([4.,6.,0.]),np.full(3,300.),np.full(3,600.),np.full(3,2000.),np.full(3,50.),np.full(3,100.))
    a=_psv_surface_compliance_single_wk(omega,k,*single)
    b=_psv_surface_compliance_single_wk(omega,k,*layered)
    np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-20)

def test_selected_quadrature_converges_on_active_and_passive_offsets():
    from run_benchmark import quadrature_check
    m,_=get_benchmark_a_lvl()
    result=quadrature_check(m,small_operator(density=2.))
    assert result['status']=='passed_5pct_selected_active_and_passive_offsets'

def test_render_refresh_rejects_changed_numerical_sources(monkeypatch):
    import json,sys
    import replot_saved_results as render
    path=Path('results_v4/benchmark_verification_summary.json')
    if not path.exists():pytest.skip('completed integration run required')
    report=json.loads(path.read_text())
    changed=dict(report.get('render_source_sha256',report['source_sha256']))
    changed['forward_mfj.py']='modified-numerical-source'
    monkeypatch.setattr(render,'code_hashes',lambda:changed)
    monkeypatch.setattr(sys,'argv',['replot_saved_results.py','--results-dir','results_v4'])
    with pytest.raises(ValueError,match='Numerical source changed: forward_mfj.py'):
        render.main()

def test_manuscript_numbers_come_from_completed_report():
    import json
    from generate_pdf_manuscript import build_manuscript
    path=Path('results_v4/benchmark_verification_summary.json')
    if not path.exists():pytest.skip('completed integration run required')
    report=json.loads(path.read_text())
    root=Path('.').resolve();text=build_manuscript(report,root/'results_v4',root)
    for key in ['benchmark_a','benchmark_b']:
        m=report[key]['mean_model_metrics']
        assert f"{m['ncc_similarity']:.5f}" in text
        assert f"{m['wasserstein_w2_mps']:.3f}" in text
    assert report['ablation_studies']['mcmc']['convergence_status'] in text
    assert 'bukan uji kalibrasi interval' in text
