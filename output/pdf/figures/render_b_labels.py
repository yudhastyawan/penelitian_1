"""Render revised Benchmark B labels from the completed v6 run."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
)
from fusion import ActivePassiveSpectralFuser, BroadbandFusionConfig
from models import get_benchmark_b_water_table
from run_benchmark import code_hashes
from visualize import (
    plot_fc_spectrum_comparison,
    plot_posterior_profiles_and_marginals,
    plot_waveforms_and_fusion,
)

source = ROOT / "results_v6"
report = json.loads((source / "benchmark_verification_summary.json").read_text())
assert report["status"] == "completed" and report["source_sha256"] == code_hashes()
hashes = json.loads(
    (ROOT / "results_v6_validation/verified_artifact_sha256.json").read_text()
)
data = []
for suffix in ["observations", "inference"]:
    path = source / f"benchmark_b_{suffix}.npz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == hashes[path.name]
    with np.load(path) as d:
        data.append(dict(d))
obs, fit = data
truth, space = get_benchmark_b_water_table()
model = replace(
    truth, name="Benchmark B: Kontras Vp pada batas muka air tanah sintetik"
)
op = FullWavefieldMFJOperator(
    ActiveSurveyGeometry(),
    PassiveSurveyGeometry(n_radial_bins=0),
    obs["freqs_active"],
    obs["freqs_passive"],
    obs["c_grid"],
    report["quadrature_density"],
)
fuser = ActivePassiveSpectralFuser(op, BroadbandFusionConfig(norm_power=1.35))
rng = np.random.default_rng(report["benchmark_b"]["config"]["random_seed"])
replay = fuser.generate_fused_spectrum_for_model(truth, True, rng)
for key in ("u_active", "u_passive", "d_fused", "w_active", "w_passive"):
    np.testing.assert_allclose(
        replay[key], obs[key], rtol=1e-10, atol=1e-20 if key.startswith("u_") else 1e-10
    )
active = op.synthesize_time_gather(truth, "active", True, rng)
passive = op.synthesize_time_gather(truth, "passive", True, rng)
plot_waveforms_and_fusion(
    model, op, fuser, obs, *active, *passive, str(HERE / "benchmark_b_gelombang_v6.png")
)
plot_fc_spectrum_comparison(
    model, fuser, obs["d_fused"], fit, str(HERE / "benchmark_b_spektrum_v6.png")
)
plot_posterior_profiles_and_marginals(
    model, space, fit, str(HERE / "benchmark_b_profil_v6.png")
)
print("Verified v6 data replay and rendered three Benchmark B panels.")
