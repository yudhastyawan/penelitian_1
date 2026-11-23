"""Measure the formerly failing deterministic controls under the repaired code."""
import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from types import SimpleNamespace
from smvi_inversion import discrete_wasserstein2,compute_spectral_misfit_numba,PickingFreeSpectralLikelihood,SMVIConfig
from fusion import ActivePassiveSpectralFuser
from run_benchmark import benchmark_frequency_grids,code_hashes

scales=[]
for high in [598.,2450.]:
    c=np.arange(130.,high+.01,4.)
    p=np.exp(-.5*((c-250)/12)**2)[None,:];q=np.exp(-.5*((c-300)/12)**2)[None,:]
    j,n,w=compute_spectral_misfit_numba(p,q,c,np.ones(1),freqs=np.array([20.]),wasserstein_scale_mps=200.)
    scales.append({'maximum_velocity':high,'J':j,'W2':w,'W2_contribution':j-.65*n})
np.testing.assert_allclose(scales[0]['W2_contribution'],scales[1]['W2_contribution'],rtol=1e-8)
fa,fp=benchmark_frequency_grids(True)
fuser=ActivePassiveSpectralFuser(SimpleNamespace(freqs_active=fa,freqs_passive=fp,c_grid=np.arange(130.,2450.01,4.)))
obs=np.ones((len(fuser.freqs_broadband),len(fuser.c_grid)))
full=PickingFreeSpectralLikelihood(fuser,SimpleNamespace(),obs,SMVIConfig())
window=PickingFreeSpectralLikelihood(fuser,SimpleNamespace(),obs,SMVIConfig(),velocity_max=600.)
assert full.wasserstein_scale_mps==window.wasserstein_scale_mps
passive=PickingFreeSpectralLikelihood(fuser,SimpleNamespace(),obs,SMVIConfig(),channel='passive_bb')
assert passive.freq_weights.sum()==11.
assert passive.freq_weights[fuser.freqs_broadband==12][0]==.25
try:discrete_wasserstein2(np.array([np.nan,1.]),np.ones(2),np.array([0.,1.]))
except ValueError:nan_rejected=True
else:raise AssertionError('nonfinite transport input was accepted')
result={'status':'passed','transport_nonfinite_rejected':nan_rejected,'fixed_scale_window_checks':scales,
        'W2_contribution_ratio':scales[0]['W2_contribution']/scales[1]['W2_contribution'],
        'full_parent_scale_mps':full.wasserstein_scale_mps,'restricted_window_scale_mps':window.wasserstein_scale_mps,
        'passive_band_integral_hz':float(passive.freq_weights.sum()),'passive_12hz_weight':float(passive.freq_weights[fuser.freqs_broadband==12][0]),
        'source_sha256':code_hashes(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'results_v6_validation/followup_controls.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
