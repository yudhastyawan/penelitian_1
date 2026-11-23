"""Check active-band misfit sensitivity on nested f–c grids, for fixed models."""

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
    compute_trapezoidal_weights,
)
from fusion import normalize_frequency_slices_numba
from models import get_benchmark_b_water_table
from run_benchmark import code_hashes
from smvi_inversion import compute_spectral_misfit_numba

report = json.loads(
    (ROOT / "results_v6/benchmark_verification_summary.json").read_text()
)
assert code_hashes() == report["source_sha256"]
cfg = report["benchmark_b"]["config"]
truth, space = get_benchmark_b_water_table()
with np.load(ROOT / "results_v6/benchmark_b_inference.npz") as data:
    centers = space.latent_to_physical(data["means_z"])
    models = {
        "prior": data["prior_midpoint_phys"],
        "model_mean": data["posterior_mean_phys"],
        **{f"center_{k}": m for k, m in enumerate(centers, 1)},
    }
f = np.arange(30.0, 42.01, 0.25)
c = np.arange(130.0, 2450.01, 2.0)
op = FullWavefieldMFJOperator(
    ActiveSurveyGeometry(),
    PassiveSurveyGeometry(n_radial_bins=0),
    f,
    np.array([8.0, 12.0]),
    c,
    report["quadrature_density"],
)


def clean(model):
    gather, _ = op.synthesize_frequency_gather(model, "active", False)
    return op.compute_mfj_spectrum(gather, "active")


obs = clean(truth)
result = []
for label, physical in models.items():
    syn = clean(space.vector_to_model(physical))
    for df, dc in [(2.0, 4.0), (1.0, 4.0), (0.5, 4.0), (0.25, 4.0), (0.25, 2.0)]:
        fi = int(df / 0.25)
        ci = int(dc / 2.0)
        o, s = [
            normalize_frequency_slices_numba(a[::fi, ::ci].copy(), 1.35)
            for a in (obs, syn)
        ]
        j, n, w = compute_spectral_misfit_numba(
            o,
            s,
            c[::ci].copy(),
            compute_trapezoidal_weights(f[::fi]),
            wasserstein_weight=cfg["wasserstein_weight"],
            wasserstein_scale_mps=cfg["wasserstein_scale_mps"],
        )
        result.append(
            {
                "model": label,
                "df_hz": df,
                "dc_mps": dc,
                "J": j,
                "NCC": 1 - n,
                "W2_mps": w,
            }
        )
out = {
    "scope": "Clean synthetic active spectra, fixed models and geometry, closed 30–42 Hz band. Dense synthesis subsampled at exact native nodes, no interpolation. No refit or field data.",
    "source_sha256": report["source_sha256"],
    "results": result,
}
(HERE / "resolution_30_42Hz.json").write_text(
    json.dumps(out, indent=2, allow_nan=False) + "\n"
)
for label in models:
    print(label, [round(r["J"], 5) for r in result if r["model"] == label])
