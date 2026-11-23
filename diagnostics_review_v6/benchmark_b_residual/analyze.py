"""Decompose the frozen v6 Benchmark B residual without refitting."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from forward_mfj import ActiveSurveyGeometry, PassiveSurveyGeometry, FullWavefieldMFJOperator
from fusion import ActivePassiveSpectralFuser, BroadbandFusionConfig
from models import get_benchmark_b_water_table
from run_benchmark import code_hashes
from smvi_inversion import PickingFreeSpectralLikelihood, SMVIConfig, compute_spectral_misfit_numba

source = ROOT / "results_v6"
report = json.loads((source / "benchmark_verification_summary.json").read_text())
assert report["status"] == "completed" and code_hashes() == report["source_sha256"]
code = "benchmark_b"
inputs = {}
for suffix in ("observations", "inference"):
    path = source / f"{code}_{suffix}.npz"
    artifact_hashes = json.loads((ROOT / "results_v6_validation/verified_artifact_sha256.json").read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == artifact_hashes[path.name]
    with np.load(path) as data:
        inputs[suffix] = dict(data)
obs, fit = inputs["observations"], inputs["inference"]
cfg = SMVIConfig(**report[code]["config"])
_, space = get_benchmark_b_water_table()
op = FullWavefieldMFJOperator(ActiveSurveyGeometry(), PassiveSurveyGeometry(n_radial_bins=0),
    obs["freqs_active"], obs["freqs_passive"], obs["c_grid"], report["quadrature_density"])
fuser = ActivePassiveSpectralFuser(op, BroadbandFusionConfig(norm_power=1.35))
channels = ("fused", "active_bb", "passive_bb")
targets = {ch: PickingFreeSpectralLikelihood(fuser, space, obs["d_"+ch], cfg, channel=ch) for ch in channels}
f = fuser.freqs_broadband
c = fuser.c_grid

def spectra(physical):
    result = fuser.generate_fused_spectrum_for_model(space.vector_to_model(physical))
    return np.stack([result["d_"+ch] for ch in channels])

prior = spectra(fit["prior_midpoint_phys"])
centers = np.array([spectra(space.latent_to_physical(z)) for z in fit["means_z"]])
weights = fit["weights"]
sums = np.zeros_like(centers)
counts = np.zeros(len(weights), dtype=int)
draw_losses = [[] for _ in weights]
for i, (z, k) in enumerate(zip(fit["predictive_samples_z"], fit["predictive_comp_labels"])):
    fields = spectra(space.latent_to_physical(z))
    sums[k] += fields
    counts[k] += 1
    draw_losses[k].append(targets["fused"].spectral_metrics(fields[0])[0])
    if (i+1) % 128 == 0:
        print(f"Replayed {i+1} draws", flush=True)
means = sums / counts[:,None,None,None]
pred = np.tensordot(weights, means, axes=1)
np.testing.assert_allclose(pred[0], fit["d_posterior_predictive_mean"], atol=1e-10)
np.testing.assert_allclose(targets["fused"].spectral_metrics(prior[0])[0],
    report[code]["baseline_comparison"]["prior_midpoint"]["misfit"], atol=1e-10)

def row_terms(o, s, active_rows):
    ncc = np.zeros(len(f))
    transport = ncc.copy()
    velocity = np.zeros((len(f),3))
    for i in np.flatnonzero(active_rows):
        limit = len(c) if f[i] > cfg.passive_cutoff_hz else np.searchsorted(c, cfg.passive_velocity_max, side="right")
        a, b = o[i,:limit], s[i,:limit]
        _, loss, w2 = compute_spectral_misfit_numba(a[None,:], b[None,:], c[:limit], np.ones(1),
            wasserstein_scale_mps=cfg.wasserstein_scale_mps)
        ncc[i] = (1-cfg.wasserstein_weight)*loss
        transport[i] = cfg.wasserstein_weight*(w2/cfg.wasserstein_scale_mps)**2
        an, bn = np.linalg.norm(a), np.linalg.norm(b)
        assert an > 0 and bn > 0
        delta = .5*(a/an-b/bn)**2*(1-cfg.wasserstein_weight)
        for j, region in enumerate((c[:limit]<=680, (c[:limit]>680)&(c[:limit]<=1200), c[:limit]>1200)):
            velocity[i,j] = delta[region].sum()
    np.testing.assert_allclose(velocity.sum(axis=1),ncc,atol=1e-11)
    return ncc,transport,velocity

bands = {"1–<8 Hz": f<8, "8–12 Hz": (f>=8)&(f<=12), ">12–30 Hz": (f>12)&(f<=30),
         ">30–42 Hz": (f>30)&(f<=42), ">42–55 Hz": f>42}
result = {"source_sha256":report["source_sha256"],"channels":{}, "components":[]}
for ch_idx,ch in enumerate(channels):
    target=targets[ch]
    w=target.freq_weights/target.freq_weights.sum()
    n0,t0,v0=row_terms(target.d_obs,prior[ch_idx],w>0)
    n1,t1,v1=row_terms(target.d_obs,pred[ch_idx],w>0)
    p0,p1=w@(n0+t0),w@(n1+t1)
    np.testing.assert_allclose([p0,p1],[target.spectral_metrics(x[ch_idx])[0] for x in (prior,pred)],atol=1e-10)
    result["channels"][ch]={"prior":p0,"prediction":p1,"delta_ncc":w@(n1-n0),"delta_w2":w@(t1-t0),
        "frequency_bands":[{"band":label,"weight":w[mask].sum(),"delta_ncc":w[mask]@(n1-n0)[mask],
        "delta_w2":w[mask]@(t1-t0)[mask],"delta_total":w[mask]@(n1+t1-n0-t0)[mask],
        "delta_ncc_by_velocity":(w[mask]@(v1-v0)[mask]).tolist()} for label,mask in bands.items() if np.any(w[mask]>0)],
        "velocity_bins_mps":["≤680",">680–1200",">1200"],"delta_ncc_by_velocity":(w@(v1-v0)).tolist()}
target=targets["fused"]
for k in range(len(weights)):
    dropped=(pred-weights[k]*means[k])/(1-weights[k])
    result["components"].append({"component":k+1,"weight":weights[k],"draws":int(counts[k]),
        "center_J":target.spectral_metrics(centers[k,0])[0],
        "conditional_mean_spectrum_J":target.spectral_metrics(means[k,0])[0],
        "mean_draw_J":float(np.mean(draw_losses[k])),
        "leave_one_component_out_J":target.spectral_metrics(dropped[0])[0]})
result["scope"]="In-sample diagnostics of the frozen v6 fit. Channel spectra use their own observations and native support, without channel-specific refits. Velocity contributions partition NCC only; Wasserstein is reported by frequency."
np.savez_compressed(HERE/"spectra.npz", frequency=f, velocity=c, prior=prior, prediction=pred,
    observation=np.stack([obs["d_"+ch] for ch in channels]), component_centers=centers, component_means=means)
(HERE/"residual_summary.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
print("Fused v6 result reproduced and decomposition saved.",flush=True)
