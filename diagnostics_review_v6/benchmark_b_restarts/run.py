"""Restart the frozen v6 Benchmark B fit on its original observation arrays."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from dataclasses import replace, asdict

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
from smvi_inversion import (
    PickingFreeSpectralLikelihood,
    SMVIConfig,
    SteinMixtureVariationalInference,
)
from run_benchmark import code_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    source = ROOT / "results_v6"
    report = json.loads((source / "benchmark_verification_summary.json").read_text())
    assert report["status"] == "completed" and report["source_sha256"] == code_hashes()
    artifacts = json.loads(
        (ROOT / "results_v6_validation/verified_artifact_sha256.json").read_text()
    )
    input_file = source / "benchmark_b_observations.npz"
    observation_hash = hashlib.sha256(input_file.read_bytes()).hexdigest()
    assert observation_hash == artifacts[input_file.name]
    b = report["benchmark_b"]
    assert args.seed != b["config"]["random_seed"]
    cfg = replace(SMVIConfig(**b["config"]), random_seed=args.seed)
    out = HERE / f"seed_{args.seed}"
    out.mkdir(exist_ok=False)
    (out / "config.json").write_text(json.dumps(asdict(cfg), indent=2))
    with np.load(input_file) as data:
        obs = dict(data)
    _, space = get_benchmark_b_water_table()
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(n_radial_bins=0),
        obs["freqs_active"],
        obs["freqs_passive"],
        obs["c_grid"],
        report["quadrature_density"],
    )
    fuser = ActivePassiveSpectralFuser(op, BroadbandFusionConfig(norm_power=1.35))
    np.testing.assert_array_equal(fuser.freqs_broadband, obs["frequencies"])
    likelihood = PickingFreeSpectralLikelihood(fuser, space, obs["d_fused"], cfg)
    print(f"Benchmark B, seed {args.seed}, iterations {cfg.n_iterations}", flush=True)
    result = SteinMixtureVariationalInference(likelihood, cfg).run(verbose=True)
    np.testing.assert_equal(len(result["history_elbo"]), cfg.n_iterations)
    np.savez_compressed(
        out / "inference.npz",
        **{k: v for k, v in result.items() if isinstance(v, np.ndarray)},
    )
    print("Saved inference. Replaying channel predictions for diagnostics.", flush=True)

    prior_model = space.latent_to_physical(np.zeros(space.dim))
    prior = fuser.generate_fused_spectrum_for_model(space.vector_to_model(prior_model))
    prior_j = likelihood.spectral_metrics(prior["d_fused"])[0]
    np.testing.assert_allclose(
        prior_j, b["baseline_comparison"]["prior_midpoint"]["misfit"], rtol=1e-10
    )
    channel_means = {
        ch: np.zeros_like(obs[f"d_{ch}"]) for ch in ("active_bb", "passive_bb")
    }
    counts = result["predictive_component_counts"]
    for z, k in zip(result["predictive_samples_z"], result["predictive_comp_labels"]):
        predicted = fuser.generate_fused_spectrum_for_model(
            space.vector_to_model(space.latent_to_physical(z))
        )
        for ch, mean in channel_means.items():
            mean += result["weights"][k] / counts[k] * predicted[f"d_{ch}"]
    band = (fuser.freqs_broadband > 30) & (fuser.freqs_broadband <= 42)
    channels = {}
    for ch, mean in channel_means.items():
        target = PickingFreeSpectralLikelihood(
            fuser, space, obs[f"d_{ch}"], cfg, channel=ch
        )
        j0, j1 = (
            float(target.spectral_metrics(s)[0]) for s in (prior[f"d_{ch}"], mean)
        )
        m = {"prior": j0, "prediction": j1, "delta": j1 - j0}
        if ch == "active_bb":
            previous = target.freq_weights.copy()
            target.freq_weights = previous * band
            change = (
                target.spectral_metrics(mean)[0]
                - target.spectral_metrics(prior[f"d_{ch}"])[0]
            )
            m["band_gt30_le42_contribution"] = float(
                change * previous[band].sum() / previous.sum()
            )
        channels[ch] = m
    vp_idx = space.param_names.index("Vp_3")
    vp = result["posterior_samples_phys"][:, vp_idx]
    pred = result["predictive_metrics"]
    summary = {
        "status": "completed",
        "seed": args.seed,
        "observation_sha256": observation_hash,
        "source_sha256": report["source_sha256"],
        "optimization_diagnostics": result["optimization_diagnostics"],
        "iterations": cfg.n_iterations,
        "posterior_predictive": pred,
        "prior_J": float(prior_j),
        "mean_model_metrics": result["mean_model_metrics"],
        "best_component_metrics": result["best_component_metrics"],
        "channel_diagnostics": channels,
        "Vp_3": {
            "mean": float(vp.mean()),
            "quantiles_95": np.percentile(vp, [2.5, 97.5]).tolist(),
            "truth": float(space.true_vector[vp_idx]),
        },
    }
    assert code_hashes() == report["source_sha256"]
    (out / "summary.json.tmp").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (out / "summary.json.tmp").replace(out / "summary.json")
    print(f"Completed seed {args.seed}: J={pred['misfit']:.6f}", flush=True)


if __name__ == "__main__":
    main()
