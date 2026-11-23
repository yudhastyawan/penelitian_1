"""Independent checks and fine-grid diagnostics for the reassessment.

Writes diagnostic evidence, never a research manuscript. The independent P-SV
reference uses PyProp8 1.1.5's propagator/source-vector implementation and the
same constant-Q moduli. It checks compliance, not real ambient interferometry.
"""

from pathlib import Path
import json
import numpy as np
from models import get_benchmark_a_lvl, get_benchmark_b_water_table
from forward_mfj import _psv_surface_compliance_single_wk


def independent_compliance_check(factory):
    from pyprop8 import LayeredStructureModel
    from pyprop8._propagators import compute_H_matrices, sourceVector

    model, _ = factory()
    h, vs, vp, rho, qs, qp = model.to_flat_arrays()
    earth = LayeredStructureModel(
        [
            (h[i] if i < len(h) - 1 else np.inf, vp[i], vs[i], rho[i])
            for i in range(len(h))
        ]
    )
    earth.mu = rho * (vs * (1 - 0.5j / qs)) ** 2
    earth.sigma = rho * (vp * (1 - 0.5j / qp)) ** 2
    # Positive source depth is required by PyProp8; approach the surface limit.
    dz, sigma, mu, rho2, isrc, irec, *_ = earth.with_interfaces(1e-6, 0.0)
    k = np.array([0.02, 0.05, 0.1, 0.2, 0.4, 0.8])
    frequencies = [2.0, 10.0, 40.0]
    errors = []
    for f in frequencies:
        w = 2 * np.pi * f * (1 + 0.005j)
        H, _ = compute_H_matrices(k, w, dz, sigma, mu, rho2, isrc, irec)
        src, _ = sourceVector(
            np.zeros((3, 3)), np.array([0.0, 0.0, 1.0]), k, sigma[isrc], mu[isrc]
        )
        ref = (H @ src).value[:, 0, 2]
        own = np.array(
            [_psv_surface_compliance_single_wk(w, x, h, vs, vp, rho, qs, qp) for x in k]
        )
        errors.extend((abs(ref - own) / abs(ref)).tolist())
    worst = float(max(errors))
    return {
        "status": "passed" if worst < 1e-4 else "failed",
        "maximum_relative_error": worst,
        "frequencies_hz": frequencies,
        "wavenumbers_per_m": k.tolist(),
        "reference": "PyProp8 layered propagator; surface source depth limit 1e-6 m",
        "scope": "Layered frequency-wavenumber compliance; not full waveform/field validation",
    }


def main():
    from run_benchmark import code_hashes

    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="results_v6_validation")
    args = parser.parse_args()
    out = Path(args.output_dir)
    if (out / "independent_physics_validation.json").exists():
        parser.error(
            "Validation record already exists; use a new --output-dir to preserve prior evidence"
        )
    out.mkdir(parents=True, exist_ok=True)
    checks = {
        factory()[0].code: independent_compliance_check(factory)
        for factory in [get_benchmark_a_lvl, get_benchmark_b_water_table]
    }
    report = {
        "status": (
            "passed"
            if all(v["status"] == "passed" for v in checks.values())
            else "failed"
        ),
        "independent_compliance": checks,
        "source_sha256": code_hashes(),
    }
    (out / "independent_physics_validation.json").write_text(
        json.dumps(report, indent=2)
    )
    print(json.dumps(report["independent_compliance"], indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
