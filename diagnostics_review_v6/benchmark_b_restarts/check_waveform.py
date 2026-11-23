"""Compare Benchmark B band-limited active traces with a separate propagator."""

import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import j0
from scipy.integrate import trapezoid
from pyprop8 import LayeredStructureModel
from pyprop8._propagators import compute_H_matrices, sourceVector

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
)
from models import get_benchmark_b_water_table
from run_benchmark import code_hashes

report = json.loads(
    (ROOT / "results_v6/benchmark_verification_summary.json").read_text()
)
assert code_hashes() == report["source_sha256"]
truth, _ = get_benchmark_b_water_table()
h, vs, vp, rho, qs, qp = truth.to_flat_arrays()
earth = LayeredStructureModel(
    [(h[i] if i < len(h) - 1 else np.inf, vp[i], vs[i], rho[i]) for i in range(len(h))]
)
earth.mu = rho * (vs * (1 - 0.5j / qs)) ** 2
earth.sigma = rho * (vp * (1 - 0.5j / qp)) ** 2
dz, sigma, mu, rho2, isrc, irec, *_ = earth.with_interfaces(1e-6, 0.0)
ratio = earth.sigma[0] / earth.mu[0]
static = (1 - (0.5 * ratio - 1) / (ratio - 1)) / earth.mu[0]
freq = np.arange(30.0, 42.01, 0.5)
offsets = np.array([4.0, 20.0, 60.0, 114.0])
op = FullWavefieldMFJOperator(
    ActiveSurveyGeometry(),
    PassiveSurveyGeometry(n_radial_bins=0),
    freq,
    np.array([8.0, 12.0]),
    np.arange(130.0, 2450.01, 4.0),
    report["quadrature_density"],
)
own = op._green_gather(truth, freq, offsets)


def independent(dk, cutoff):
    k = np.arange(1, int(np.ceil(cutoff / dk)) + 1) * dk
    g = []
    for f in freq:
        H, _ = compute_H_matrices(
            k, 2 * np.pi * f * (1 + 0.005j), dz, sigma, mu, rho2, isrc, irec
        )
        src = sourceVector(
            np.zeros((3, 3)), np.array([0.0, 0.0, 1.0]), k, sigma[isrc], mu[isrc]
        )[0]
        integrand = (k * (H @ src).value[:, 0, 2] - static)[:, None] * j0(
            k[:, None] * offsets
        )
        integrand = np.vstack([np.full(len(offsets), -static), integrand])
        g.append(trapezoid(integrand, np.r_[0.0, k], axis=0) + static / offsets)
    return np.array(g)


ref = independent(0.001, 20.0)
refine = independent(0.0005, 25.0)
x = (freq / op.active_geom.f_center) ** 2
wavelet = x * np.exp(-x)
bins = np.rint(freq / 0.5).astype(int)


def to_time(g):
    spec = np.zeros((257, len(offsets)), complex)
    spec[bins] = np.conj(g) * wavelet[:, None]
    return np.fft.irfft(spec, n=512, axis=0)


def relative(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


record = {
    "scope": "Independent PyProp8 layered propagator and linear k quadrature with shared physical moduli, static-tail identity, source and geometry. Noise-free; inverse FFT limited to closed 30–42 Hz band.",
    "source_sha256": report["source_sha256"],
    "frequencies_hz": freq.tolist(),
    "offsets_m": offsets.tolist(),
    "reference_refinement_relative_l2": relative(ref, refine),
    "complex_spectrum_relative_l2": relative(own, refine),
    "bandlimited_waveform_relative_l2_by_offset": [
        relative(a, b) for a, b in zip(to_time(own).T, to_time(refine).T)
    ],
}
(HERE / "independent_waveform_30_42Hz.json").write_text(
    json.dumps(record, indent=2, allow_nan=False) + "\n"
)
print(
    json.dumps(
        {
            k: v
            for k, v in record.items()
            if k not in ("source_sha256", "frequencies_hz")
        },
        indent=2,
    )
)
