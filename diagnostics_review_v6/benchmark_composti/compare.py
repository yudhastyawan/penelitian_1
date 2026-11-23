"""Compare v6 with the COMPOSTI Levin solver using the active A/B geometry."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXTERNAL = ROOT / "tmp/composti_external"
UPSTREAM = EXTERNAL / "upstream"


def run_composti(exchange):
    import reflectivityCPP

    with np.load(exchange) as data:
        layers = np.asfortranarray(data["layers"])
        frequency = data["frequency"]
        offsets = data["offsets"]
        wavelet = data["wavelet"].astype(np.complex128)
    output = np.zeros((len(frequency), len(offsets)), dtype=np.complex128, order="F")
    # Include direct arrival and multiples, exclude frequency taper, return displacement.
    options = np.array([1, 1, 0, 0], dtype=np.int32)
    reflectivityCPP.compute_displ_Levin(layers, frequency, wavelet, offsets, options, output)
    if not np.isfinite(output).all() or not np.any(output):
        raise ValueError("COMPOSTI produced invalid or empty output")
    np.save(exchange.with_suffix(".npy"), output)


def native_without_shift(model, operator, frequency, offsets, nlog, cutoff_multiplier=4,
                         matched_slowness=False):
    """Evaluate the same v6 compliance with real omega and a configurable quadrature."""
    from scipy.special import j0
    from forward_mfj import _psv_surface_compliance_grid_numba, compute_trapezoidal_weights

    h, vs, vp, rho, qs, qp = model.to_flat_arrays()
    mu = rho[0] * (vs[0] * (1 - 0.5j / qs[0])) ** 2
    ratio = (vp[0] * (1 - 0.5j / qp[0]) / (vs[0] * (1 - 0.5j / qs[0]))) ** 2
    static = (1 - (0.5 * ratio - 1) / (ratio - 1)) / mu
    result = []
    for f in frequency:
        omega = 2 * np.pi * f
        if matched_slowness:
            k_min = omega * 1e-3 / vp.max()
            k_taper = omega * 1.2 / vs.min()
            k_max = 5 * k_taper
            k = np.r_[np.geomspace(k_min, k_taper, nlog)[:-1],
                      np.linspace(k_taper, k_max, nlog)]
        else:
            k_split = omega / 5000
            k_max = cutoff_multiplier * max(
                3 * omega / max(operator.c_grid[0] / 2, 40), 8 / offsets.min()
            )
            k = np.r_[np.linspace(1e-10, k_split, 33)[:-1],
                      np.geomspace(k_split, k_max, nlog)]
        compliance = _psv_surface_compliance_grid_numba(
            np.array([f]), k, 0.0, h, vs, vp, rho, qs, qp
        )[0]
        integrand = k * compliance
        if matched_slowness:
            q = np.clip((k - k_taper) / (k_max - k_taper), 0, 1)
            integrand *= 0.5 * (1 + np.cos(np.pi * q))
            correction = 0.0
        else:
            integrand -= static
            correction = static / offsets
        result.append(
            (integrand * compute_trapezoidal_weights(k)) @ j0(k[:, None] * offsets[None, :])
            + correction
        )
    return np.array(result) / (2 * np.pi)


def time_gather(spectrum, frequency):
    nt, dt = 1024, 1 / 512
    bins = np.rint(frequency * nt * dt).astype(int)
    np.testing.assert_allclose(bins / (nt * dt), frequency)
    positive = np.zeros((nt // 2 + 1, spectrum.shape[1]), complex)
    positive[bins] = np.conj(spectrum)
    return np.arange(nt) * dt, np.fft.irfft(positive, n=nt, axis=0)


def relative_l2(estimate, reference, axis=None):
    return np.linalg.norm(estimate - reference, axis=axis) / np.linalg.norm(reference, axis=axis)


def compare_case(label, model):
    from forward_mfj import ActiveSurveyGeometry, FullWavefieldMFJOperator, PassiveSurveyGeometry

    with np.load(ROOT / f"results_v6/benchmark_{label}_observations.npz") as stored:
        frequency = stored["freqs_active"]
        offsets = stored["active_offsets"]
        c_grid = stored["c_grid"]
    operator = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(), PassiveSurveyGeometry(n_stations=4, use_option4_2d_nested=False),
        frequency, np.array([1.0, 2.0]), c_grid, quadrature_density=2.0,
    )
    np.testing.assert_array_equal(operator.r_act, offsets)
    dense = np.arange(frequency.min(), frequency.max() + 0.25, 0.5)
    selected = np.searchsorted(dense, frequency)
    np.testing.assert_array_equal(dense[selected], frequency)
    x = (dense / operator.active_geom.f_center) ** 2
    wavelet = x * np.exp(-x)
    h, vs, vp, rho, qs, qp = model.to_flat_arrays()
    exchange = EXTERNAL / f"benchmark_{label}_input.npz"
    np.savez(exchange, layers=np.column_stack([vp, qp, vs, qs, rho, h]),
             frequency=dense, offsets=offsets, wavelet=wavelet)
    env = {**os.environ, "PYTHONPATH": str(UPSTREAM / "src/reflectivityCPP")}
    start = perf_counter()
    subprocess.run(
        [str(EXTERNAL / "venv/bin/python"), str(__file__), "--worker", str(exchange)],
        env=env, check=True,
    )
    external_sec = perf_counter() - start
    # COMPOSTI follows the positive-frequency convention of irfft. Convert to the
    # native exp(-i omega t) wavefield convention before the common MFJ transform.
    external = np.conj(np.load(exchange.with_suffix(".npy")))
    native = operator._green_gather(model, dense, offsets, cache=False)
    native *= wavelet[:, None] / (2 * np.pi)
    same_grid = native_without_shift(model, operator, dense, offsets, 2048, cutoff_multiplier=1)
    q_only = native_without_shift(model, operator, dense, offsets, 8192)
    matched = native_without_shift(model, operator, dense, offsets, 8192, matched_slowness=True)
    refinement = native_without_shift(model, operator, dense, offsets, 4096, matched_slowness=True)
    wavefields = {
        "composti": external,
        "v6": native,
        "v6_no_shift_same_grid": same_grid * wavelet[:, None],
        "v6_no_shift_refined_grid": q_only * wavelet[:, None],
        "v6_q_and_slowness_matched": matched * wavelet[:, None],
    }
    images = {
        key: operator.compute_mfj_spectrum(value[selected])
        for key, value in wavefields.items()
    }
    traces = {key: time_gather(value, dense)[1] for key, value in wavefields.items()}
    time = time_gather(external, dense)[0]
    band = (dense >= 30) & (dense <= 42)
    trace_offsets = [4, 20, 60, 114]
    indices = np.searchsorted(offsets, trace_offsets)
    np.testing.assert_array_equal(offsets[indices], trace_offsets)
    report = {"model": model.name, "composti_wall_sec": external_sec, "trace_offsets_m": trace_offsets,
              "native_frequency_hz": frequency.tolist(), "trace_frequency_hz": dense.tolist(),
              "offsets_m": offsets.tolist(),
              "matched_quadrature_refinement_relative_l2": float(relative_l2(
                  refinement * wavelet[:, None], wavefields["v6_q_and_slowness_matched"]
              ))}
    for key in wavefields:
        if key == "composti":
            continue
        ref_img, img = images["composti"], images[key]
        report[key] = {
            "native_complex_relative_l2": float(relative_l2(wavefields[key][selected], external[selected])),
            "30_42_hz_complex_relative_l2": float(relative_l2(wavefields[key][band], external[band])),
            "trace_relative_l2_by_offset": relative_l2(traces[key], traces["composti"], axis=0).tolist(),
            "mfj_relative_l2": float(relative_l2(img, ref_img)),
            "mfj_cosine_similarity": float(np.vdot(img, ref_img) / (np.linalg.norm(img) * np.linalg.norm(ref_img))),
        }
    np.savez_compressed(
        HERE / f"benchmark_{label}_comparison.npz",
        frequencies_hz=dense, image_frequencies_hz=frequency, offsets_m=offsets,
        time_s=time, phase_velocity_mps=c_grid,
        **{f"u_{key}": value for key, value in wavefields.items()},
        **{f"image_{key}": value for key, value in images.items()},
        **{f"trace_{key}": value for key, value in traces.items()},
    )
    make_figure(label, frequency, c_grid, offsets, time, traces, images, indices)
    print(label.upper(), {key: {
        metric: round(value, 6) for metric, value in report[key].items()
        if isinstance(value, float)
    } for key in wavefields if key != "composti"}, flush=True)
    return report


def make_figure(label, frequency, c_grid, offsets, time, traces, images, indices):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(14, 6.6), layout="constrained")
    fig.suptitle(f"Benchmark {label.upper()} | COMPOSTI dan forward v6 pada kanal aktif", fontsize=15)
    for ax, i in zip(axes[0], indices):
        scale = np.max(np.abs(traces["composti"][:, i]))
        for key, legend, style in [("composti", "COMPOSTI", "-"), ("v6", "v6", "--"),
                                   ("v6_q_and_slowness_matched", "v6 diselaraskan", ":")]:
            ax.plot(time, traces[key][:, i] / scale, style, lw=1.2, label=legend)
        ax.set(title=f"Offset {offsets[i]:.0f} m", xlabel="Waktu (s)", xlim=(0, 0.95))
        ax.grid(alpha=0.2)
    axes[0, 0].set_ylabel("Perpindahan relatif")
    axes[0, -1].legend(fontsize=8)
    labels = ["COMPOSTI", "v6", "v6 diselaraskan", "v6 − COMPOSTI"]
    fields = [images["composti"], images["v6"], images["v6_q_and_slowness_matched"],
              images["v6"] - images["composti"]]
    for i, (ax, title, img) in enumerate(zip(axes[1], labels, fields)):
        is_delta = i == 3
        mesh = ax.pcolormesh(frequency, c_grid, img.T, shading="auto",
                             cmap="RdBu_r" if is_delta else "magma",
                             vmin=-1 if is_delta else 0, vmax=1)
        ax.set(title=title, xlabel="Frekuensi (Hz)")
        fig.colorbar(mesh, ax=ax, shrink=0.78)
    axes[1, 0].set_ylabel("Kecepatan fase (m/s)")
    fig.savefig(HERE / f"benchmark_{label}_composti.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    if args.worker is not None:
        run_composti(args.worker)
        return
    sys.path.insert(0, str(ROOT))
    from models import get_benchmark_a_lvl, get_benchmark_b_water_table

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=UPSTREAM, text=True).strip()
    report = {
        "composti_commit": commit,
        "composti_api": "reflectivityCPP.compute_displ_Levin",
        "scope": "Clean active vertical surface displacement, source 25 Hz, same native grid and offsets.",
        "macos_port": "Clang serial build, POSIX jn in place of std::cyl_bessel_j (orders 0 and 1).",
        "v6_scaling": "Divide by 2*pi for standard point-force normalization.",
        "imaginary_frequency_shift_v6_ratio": 0.005,
        "imaginary_frequency_shift_controls_ratio": 0.0,
        "benchmarks": {label: compare_case(label, factory()[0]) for label, factory in
                       [("a", get_benchmark_a_lvl), ("b", get_benchmark_b_water_table)]},
    }
    (HERE / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
