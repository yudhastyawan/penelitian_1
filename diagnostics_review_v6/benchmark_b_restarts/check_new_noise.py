"""Evaluate fixed Benchmark B predictions on eight new noise realizations."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
)
from fusion import ActivePassiveSpectralFuser, BroadbandFusionConfig
from models import get_benchmark_b_water_table
from run_benchmark import code_hashes
from smvi_inversion import PickingFreeSpectralLikelihood, SMVIConfig

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prepare", action="store_true")
args = parser.parse_args()
source = ROOT / "results_v6"
report = json.loads((source / "benchmark_verification_summary.json").read_text())
assert code_hashes() == report["source_sha256"]
cfg = SMVIConfig(**report["benchmark_b"]["config"])
truth, space = get_benchmark_b_water_table()
with np.load(source / "benchmark_b_observations.npz") as d:
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(n_radial_bins=0),
        d["freqs_active"],
        d["freqs_passive"],
        d["c_grid"],
        report["quadrature_density"],
    )
fuser = ActivePassiveSpectralFuser(op, BroadbandFusionConfig(norm_power=1.35))
seeds = np.arange(3101, 3109)
data_path = HERE / "noise_evaluation_inputs.npz"
if args.prepare:
    spectra = np.array(
        [
            fuser.generate_fused_spectrum_for_model(
                truth, True, np.random.default_rng(int(s))
            )["d_fused"]
            for s in seeds
        ]
    )
    np.savez_compressed(
        data_path, seeds=seeds, fused=spectra, frequencies=fuser.freqs_broadband
    )
    print("Saved eight noise realizations for evaluation only.")
    sys.exit(0)
with np.load(data_path) as data:
    np.testing.assert_array_equal(data["seeds"], seeds)
    np.testing.assert_array_equal(data["frequencies"], fuser.freqs_broadband)
    spectra = data["fused"]
with np.load(source / "benchmark_b_inference.npz") as d:
    prior = space.vector_to_model(d["prior_midpoint_phys"])
    predictions = {2027: d["d_posterior_predictive_mean"]}
for fit_seed in (2028, 2029, 2030):
    path = HERE / f"seed_{fit_seed}/inference.npz"
    with np.load(path) as d:
        predictions[fit_seed] = d["d_posterior_predictive_mean"]
prior_spec = fuser.generate_fused_spectrum_for_model(prior)["d_fused"]
rows = []
for seed, o in zip(seeds, spectra):
    lik = PickingFreeSpectralLikelihood(fuser, space, o, cfg)
    j0 = lik.spectral_metrics(prior_spec)[0]
    rows.append(
        {
            "noise_seed": int(seed),
            "J_prior": j0,
            "by_fit_seed": {
                s: {"J": j, "delta_from_prior": j - j0}
                for s, p in predictions.items()
                for j in [lik.spectral_metrics(p)[0]]
            },
        }
    )
out = {
    "scope": "Eight new noise realizations of the same synthetic truth and forward model. Predictions held fixed at their fit to the original v6 observations. No channel-specific refits.",
    "source_sha256": report["source_sha256"],
    "results": rows,
}
(HERE / "new_noise_evaluation.json").write_text(
    json.dumps(out, indent=2, allow_nan=False) + "\n"
)
for s in predictions:
    d = np.array([r["by_fit_seed"][s]["delta_from_prior"] for r in rows])
    print(
        s, "median delta", float(np.median(d)), "wins", int(np.sum(d < 0)), "of", len(d)
    )
