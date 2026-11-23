"""Summarize three restarts against the original v6 fit and its prior."""

import json
from pathlib import Path

import numpy as np

out = Path(__file__).resolve().parent
root = out.parents[1]
source = root / "results_v6"
report = json.loads((source / "benchmark_verification_summary.json").read_text())
b = report["benchmark_b"]
residual = json.loads(
    (out.parent / "benchmark_b_residual/residual_summary.json").read_text()
)
names = list(b["parameters"])
with np.load(source / "benchmark_b_observations.npz") as a:
    lower, upper = a["prior_lower"], a["prior_upper"]
assert np.all(lower[4:8] >= 1.45 * upper[:4]), "Conditional Vp bounds are active"
indices = {name: names.index(name) for name in ("Vp_2", "Vp_3", "h_1")}
baseline = float(b["baseline_comparison"]["prior_midpoint"]["misfit"])
print("Prior J", baseline)
rows = []
params = []
prior_intervals = {
    name: (lower[i] + np.array([0.025, 0.975]) * (upper[i] - lower[i]))
    for name, i in indices.items()
}
for seed in (2027, 2028, 2029, 2030):
    if seed == 2027:
        sample_path = source / "benchmark_b_inference.npz"
        p = b["predictive_metrics"]
        active = residual["channels"]["active_bb"]
        band = next(
            x["delta_total"]
            for x in active["frequency_bands"]
            if x["band"] == ">30–42 Hz"
        )
        active_delta = active["prediction"] - active["prior"]
        passive = residual["channels"]["passive_bb"]
        passive_delta = passive["prediction"] - passive["prior"]
    else:
        sample_path = out / f"seed_{seed}/inference.npz"
        r = json.loads(sample_path.with_name("summary.json").read_text())
        assert (
            r["status"] == "completed"
            and r["iterations"] == b["config"]["n_iterations"]
        )
        assert r["source_sha256"] == report["source_sha256"]
        np.testing.assert_allclose(r["prior_J"], baseline, rtol=1e-10)
        p = r["posterior_predictive"]
        active = r["channel_diagnostics"]["active_bb"]
        band = active["band_gt30_le42_contribution"]
        active_delta = active["delta"]
        passive_delta = r["channel_diagnostics"]["passive_bb"]["delta"]
    rows.append(
        (
            seed,
            p["misfit"],
            p["misfit_mc_se"],
            p["draws"],
            p["precision_target_met"],
            active_delta,
            band,
            passive_delta,
        )
    )
    with np.load(sample_path) as arrays:
        s = arrays["posterior_samples_phys"]
        v = [
            (
                name,
                float(s[:, i].mean()),
                *np.percentile(s[:, i], [2.5, 97.5]),
                b["parameters"][name]["true"],
            )
            for name, i in indices.items()
        ]
        params.extend((seed, *x) for x in v)

lines = [
    "# Restart Benchmark B pada observasi v6",
    "",
    "Tiga restart memakai seed 2028, 2029, dan 2030. Seed 2027 adalah run asli. Setiap restart menggunakan konfigurasi SMVI, 60 iterasi, geometri, dan observasi asli yang sama. Seed mengubah inisialisasi, sampling optimasi, dan sampling prediksi.",
    "",
    f"Misfit prior spektrum gabungan **{baseline:.5f}**. Δ positif berarti misfit prediksi di atas prior. Pita aktif memakai $30<f\leq42$ Hz dengan bobot v6.",
    "",
    "| Seed | J gabungan | Galat baku MC | Draw | Target presisi | Δ aktif | Kontribusi pita aktif | Δ pasif |",
    "|--:|--:|--:|--:|:--|--:|--:|--:|",
]
for seed, j, se, draw, precision, act, band, pas in rows:
    lines.append(
        f"| {seed} | {j:.5f} | {se:.5f} | {draw} | {'terpenuhi' if precision else 'belum'} | {act:+.5f} | {band:+.5f} | {pas:+.5f} |"
    )
lines += [
    "",
    "## Parameter target",
    "",
    "| Seed | Parameter | Rata-rata | Kuantil 95% campuran | Kuantil 95% prior | Acuan |",
    "|--:|:--|--:|--:|--:|--:|",
]
for seed, name, mean, lo, hi, truth in params:
    plo, phi = prior_intervals[name]
    lines.append(
        f"| {seed} | {name} | {mean:.2f} | [{lo:.2f}, {hi:.2f}] | [{plo:.2f}, {phi:.2f}] | {truth:.2f} |"
    )
lines += [
    "",
    "Satuan Vp dalam m/s dan h dalam m. Prior ketiga parameter berbentuk seragam pada batas yang disimpan. Kuantil campuran mengukur sebaran parameter bersyarat pada satu observasi dan satu hasil pengoptimasian.",
    "",
    "## Pemeriksaan operator pada pita 30–42 Hz",
    "",
    "Uji grid memakai spektrum bersih dari model yang tetap dan data acuan sintetik. Pembandingan gelombang memakai model berlapis Benchmark B, empat offset, serta pita 30–42 Hz. Lihat hasil JSON masing-masing untuk grid, frekuensi, dan galat.",
    "",
    "## Cakupan",
    "",
    "Perbandingan terhadap prior memakai observasi v6 yang sama. Validasi data baru, ablasi dengan refit kanal terpisah, dan pengukuran cakupan kuantil pada model acuan berulang merupakan eksperimen tersendiri.",
]
(out / "RINGKASAN.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines[:14]))
