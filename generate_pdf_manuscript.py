"""Generate the current manuscript from the completed experiment JSON.

Markdown is generated once; Pandoc creates LaTeX from that same text, then
pdflatex compiles it. Scalar experiment results are never hardcoded here.
"""

import argparse, json, shutil, subprocess
from pathlib import Path

NAME = "NASKAH_LENGKAP_PENELITIAN_PHD_MFJ_SMVI"


def build_manuscript(report, results_dir, root):
    mode = report["mode"]
    cfg = report["config"]
    lines = [
        "# Inversi Spektrum Aktif-Pasif MFJ-SMVI: Implementasi dan Verifikasi Numerik",
        "",
        f"Versi: `{report['version']}`. Mode eksperimen: `{mode}`.",
        f"Waktu mulai UTC: {report['run_started_utc']}. Durasi total: {report['total_runtime_sec']:.1f} detik.",
        "",
        "## Abstrak",
        "",
        "Kerangka ini menggabungkan sintesis fungsi Green P-SV berlapis, pencitraan MFJ, "
        "fusi aktif-pasif, dan inferensi campuran Gaussian berkorelasi. Revisi ini menghapus "
        "cabang modal heuristik, perubahan manual sampel pembanding, serta angka hasil yang "
        "ditentukan sebelumnya. Misfit memakai NCC dan transport Wasserstein kuadratik diskrit. "
        "Perbandingan dilakukan melalui inversi terpisah, noise ditambahkan pada gather sebelum "
        "imaging, dan seluruh sampel serta konfigurasi disimpan.",
        "",
        "**Status bukti:** eksperimen sintetik yang berbagi solver antara generasi dan inversi. "
        "Keberhasilan menjalankan alur tidak membuktikan akurasi lapangan, konvergensi posterior, "
        "atau kalibrasi interval kredibel. Hasil anggaran cepat tidak dipakai untuk menyatakan "
        "keunggulan metodologi.",
        "",
        "## 1. Masalah penelitian dan asumsi",
        "",
        "Benchmark A menguji lapisan lunak terjebak pada kedalaman 4-10 m dengan Vs 150 m/s. "
        "Benchmark B menguji batas jenuh pada kedalaman 3.5 m dengan lonjakan Vp 550 menjadi "
        "1650 m/s. Model terdiri dari tiga lapisan dan half-space; densitas dan Q dianggap "
        "diketahui. Rasio Poisson diturunkan dari Vp dan Vs, bukan parameter independen.",
        "",
        "## 2. Metode yang diimplementasikan",
        "",
        (root / "MATHEMATICAL_FORMULATION.md")
        .read_text()
        .split("\n", 2)[2]
        .replace("\n## ", "\n### "),
        "",
        "## 3. Konfigurasi eksperimen",
        "",
        f"Komponen campuran: {cfg['n_components']}; iterasi: {cfg['n_iterations']}; "
        f"sampel Monte Carlo per komponen/iterasi: {cfg['n_mc_samples_per_comp']}; "
        f"sampel posterior: {cfg['n_posterior_draws']}. Sigma likelihood: {cfg['sigma_likelihood']}. "
        f"Densitas kuadratur: {report['quadrature_density']}.",
        "",
        "Frekuensi sintesis aktif: "
        + ", ".join(f"{x:.2f}" for x in report["frequency_grids"]["active"])
        + " Hz.",
        "",
        "Frekuensi sintesis pasif: "
        + ", ".join(f"{x:.2f}" for x in report["frequency_grids"]["passive"])
        + " Hz.",
        "",
        "Prior, koordinat, offset aktual, spektrum input, model sebenarnya dan sampel tersimpan "
        "di berkas NPZ. SHA256 kode dicatat dalam run manifest. Model sebenarnya dipakai "
        "pada evaluasi pemulihan parameter, bukan inisialisasi optimizer.",
        "",
        "Referensi NCC untuk setiap transformasi adalah hasil tanpa noise dari transformasi "
        "dan geometri yang sama. Nilai tersebut mengukur ketahanan terhadap noise, "
        "bukan peringkat ketepatan fisika antar-transformasi.",
        "",
        "## 4. Hasil benchmark",
        "",
    ]
    for key in ["benchmark_a", "benchmark_b"]:
        b = report[key]
        m = b["mean_model_metrics"]
        lines += [
            f"### {b['benchmark_name']}",
            "",
            f"NCC pada model posterior mean: {m['ncc_similarity']:.5f}; "
            f"W2: {m['wasserstein_w2_mps']:.3f} m/s; misfit komposit: {m['misfit']:.5f}. "
            f"RMSE relatif profil Vs: {b['vs_profile_relative_rmse_pct']:.3f}%; "
            f"MAPE profil Vs: {b['vs_profile_mape_pct']:.3f}%.",
            "",
            f"Coverage empiris parameter: {b['empirical_parameter_coverage']}. "
            "Coverage ini bukan uji kalibrasi interval pada dataset berulang.",
            "",
            "| Parameter | Nilai sejati | Mean | SD | CI 95% | Error mean (%) |",
            "|:--|--:|--:|--:|:--|--:|",
        ]
        for name, p in b["parameters"].items():
            lines.append(
                f"| {name.replace('_', ' ')} | {p['true']:.2f} | {p['posterior_mean']:.2f} | {p['posterior_std']:.2f} | [{p['ci95_low']:.2f}, {p['ci95_high']:.2f}] | {p['rel_error_pct']:.2f} |"
            )
        lines += [
            "",
            f"Evaluasi forward: {b['forward_evaluations']}; durasi inferensi: {b['elapsed_sec']:.1f} detik.",
            "",
        ]
    ab = report["ablation_studies"]
    lines += [
        "## 5. Eksperimen pembanding",
        "",
        "| Metode | RMSE relatif profil Vs (%) | MAPE profil Vs (%) | Evaluasi forward |",
        "|:--|--:|--:|--:|",
    ]
    for name in [
        "active_only",
        "passive_only",
        "joint",
        "rayleigh_only",
        "single_gaussian",
    ]:
        d = ab[name]
        lines.append(
            f"| {name.replace('_',' ')} | {d['vs_profile_relative_rmse_pct']:.3f} | {d['vs_profile_mape_pct']:.3f} | {d['forward_evaluations']} |"
        )
    mc = ab["mcmc"]
    lines += [
        "",
        "Rayleigh-only menggunakan Benchmark B; metode lain pada tabel memakai Benchmark A. "
        "Anggaran inferensi harus diperiksa sebelum membandingkan efisiensi.",
        "",
        f"MCMC: `{mc['convergence_status']}`. Acceptance rate: {mc['acceptance_rate']:.3f}. "
        f"Minimum ESS estimasi: {min(mc['ess']):.1f}. "
        f"Mean Vs2: {mc['vs2_mean']:.3f} m/s; SD: {mc['vs2_std']:.3f} m/s. "
        "MCMC tidak disebut ground truth selama konvergensi dan eksplorasi moda belum dibuktikan.",
        "",
        f"Picked-curve: `{ab['picked_curve']['status']}`. "
        "Pembanding ini menghasilkan estimasi titik melalui least squares, bukan posterior.",
        "",
        "## 6. Stress-test pada domain rekaman",
        "",
        "| Kasus | SNR aktif terukur (dB) | SNR pasif terukur (dB) | RMSE relatif Vs (%) |",
        "|:--|--:|--:|--:|",
    ]
    for name, d in report["waveform_noise_stress_test"].items():
        if isinstance(d, dict):
            lines.append(
                f"| {name.replace('_',' ')} | {d['measured_active_snr_db']:.3f} | {d['measured_passive_snr_db']:.3f} | {d['vs_profile_relative_rmse_pct']:.3f} |"
            )
    q = report["quadrature_check"]
    lines += [
        "",
        "Sumber lemah dimodelkan dengan pelemahan frekuensi dan offset sebelum transformasi. "
        "Background mengacu sumber kuat awal. Kurva fusi merupakan keluaran fungsi fusi "
        "aktual; tidak ada crossover frekuensi yang dipaksakan.",
        "",
        "## 7. Verifikasi dan keterbatasan",
        "",
        "Regresi mencakup compliance analitik, kesetaraan lapisan homogen, Jacobian prior, "
        "W2 diskrit, densitas dan sampling berkorelasi, optimasi target Gaussian, "
        "inisialisasi tanpa model sejati, MCMC tanpa clipping, serta SNR terukur.",
        "",
        f"Perubahan relatif L2 gather pada penggandaan densitas kuadratur: "
        f"{q['relative_complex_gather_l2_change']:.6f}. Status: {q['status'].replace('_', ' ')}. "
        "Pemeriksaan ini mencakup offset aktif dan pasif terpilih, bukan seluruh "
        "model dalam prior. Detail perubahan maksimum per titik tersedia di JSON.",
        "",
        "Rekaman waktu pada lampiran adalah realisasi representatif pada grid FFT rapat; "
        "inversi mencocokkan gather frekuensi tersimpan, bukan rekaman waktu tersebut. ",
        "",
        "Generasi dan inversi memakai solver yang sama. Noise pasif masih berupa model "
        "fenomenologis. Kuadratur terhingga memerlukan pemeriksaan sensitivitas terhadap "
        "jumlah titik dan batas integrasi. SMVI adalah optimizer eksperimental; kalibrasi "
        "CI, multi-start, pengulangan noise, dan sensitivitas prior belum diselesaikan. "
        "Tidak ada klaim speedup, global MAP, atau keunggulan universal geometri tertentu.",
        "",
        "## 8. Pekerjaan riset berikutnya",
        "",
        "Prioritas lanjutan adalah konvergensi kuadratur pada seluruh offset, perbandingan "
        "waveform dengan solver independen, kalibrasi posterior pada dataset berulang, "
        "MCMC dari beberapa awal independen, dan validasi lapangan terhadap logging serta "
        "piezometer. Sasaran publikasi dalam proposal bergantung pada bukti tersebut.",
        "",
        "## Referensi metode",
        "",
        "- Kennett (1983), Seismic Wave Propagation in Stratified Media: rekursi medium berlapis.",
        "- Goodman dan Weare (2010), Ensemble samplers with affine invariance: stretch-move MCMC.",
        "- [POT, dokumentasi transport 1D](https://pythonot.github.io/gen_modules/ot.lp.html#ot.lp.wasserstein_1d).",
        "- [emcee, diagnostik autokorelasi](https://emcee.readthedocs.io/en/stable/tutorials/autocorr/).",
        "",
        "Lampiran berikut memakai halaman lanskap agar label gambar lebih mudah dibaca.",
        "",
    ]
    figures = [
        (
            "fig0_field_acquisition_design_comparison.png",
            "Geometri dan diagnostik array terukur.",
        ),
        (
            "fig0b_dispersion_image_comparison.png",
            "Kernel transformasi sebenarnya pada empat geometri.",
        ),
        (
            "benchmark_a_fig1_waveforms_and_fusion.png",
            "Benchmark A: waveform, spektrum, dan bobot fusi.",
        ),
        (
            "benchmark_a_fig2_fc_spectrum_comparison.png",
            "Benchmark A: spektrum observasi dan posterior mean.",
        ),
        (
            "benchmark_a_fig3_posterior_profiles_and_marginals.png",
            "Benchmark A: profil serta ketidakpastian eksperimental.",
        ),
        (
            "benchmark_b_fig1_waveforms_and_fusion.png",
            "Benchmark B: waveform, spektrum, dan bobot fusi.",
        ),
        (
            "benchmark_b_fig2_fc_spectrum_comparison.png",
            "Benchmark B: spektrum observasi dan posterior mean.",
        ),
        (
            "benchmark_b_fig3_posterior_profiles_and_marginals.png",
            "Benchmark B: profil serta ketidakpastian eksperimental.",
        ),
        (
            "fig4_ablation_and_mcmc.png",
            "Ablasi independen dan MCMC dengan status diagnostik.",
        ),
        (
            "fig5_waveform_noise_stress.png",
            "Noise pada gather dan bobot fusi yang dihitung ulang.",
        ),
    ]
    for name, caption in figures:
        path = results_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        rel = path.relative_to(root) if path.is_relative_to(root) else path
        lines += [
            "```{=latex}",
            r"\clearpage",
            r"\begin{landscape}",
            "```",
            "",
            f"![{caption}]({rel})" + "{height=85%}",
            "",
            "```{=latex}",
            r"\end{landscape}",
            "```",
            "",
        ]
    return "\n".join(lines) + "\n"


def build_diagnostic_report(report, editorial_notes=""):
    """Regenerate critical findings from structured data; notes live separately."""
    from research_quality import assess_research_readiness

    checks = assess_research_readiness(report)
    lines = [
        "# Laporan diagnostik — bukan hasil penelitian tervalidasi",
        "",
        f"Versi: {report['version']}. Mode: {report['mode']}.",
        "",
        "Status penerbitan: " + ("lolos" if checks["publishable"] else "ditahan"),
        "",
        "Pemeriksaan yang belum lolos: " + ", ".join(checks["blocking_checks"]),
        "",
    ]
    for key in ["benchmark_a", "benchmark_b"]:
        b = report[key]
        pred = b.get("predictive_metrics")
        lines += [
            f"## {b['benchmark_name']}",
            "",
            f"Metrik pada model mean: `{json.dumps(b['mean_model_metrics'])}`.",
            f"Metrik prediksi rata-rata sampel: `{json.dumps(pred)}`.",
            "",
            "Konvergensi posterior belum dinyatakan. Angka ini adalah diagnosis numerik.",
            "",
            "### Perbandingan terhadap prior dan peringatan",
            "",
        ]
        baseline = b.get("baseline_comparison")
        if baseline:
            lines += [
                f"Model tengah prior: `{json.dumps(baseline['prior_midpoint'])}`.",
                f"Evaluasi model sejati (validasi saja): `{json.dumps(baseline['truth_model_evaluation_only'])}`.",
                "",
                "**Peringatan hasil:** "
                + (
                    ", ".join(baseline["warnings"])
                    if baseline["warnings"]
                    else "Tidak ada flag pembanding; ini tidak membuktikan validasi posterior."
                ),
                "",
            ]
        else:
            lines += [
                "**Perbandingan prior belum tersedia dalam data terstruktur; belum diperiksa.**",
                "",
            ]
        if pred:
            lines += [
                f"Galat Monte Carlo misfit: {pred.get('misfit_mc_se','belum dihitung')}. "
                f"Jumlah draw: {pred.get('draws','tidak tersedia')}. "
                f"Target presisi terpenuhi: {pred.get('precision_target_met','belum diperiksa')}.",
                "Galat Monte Carlo adalah ketelitian estimasi dari campuran yang dipasang, bukan kalibrasi posterior.",
                "",
            ]
    if editorial_notes.strip():
        lines += ["## Catatan editorial terpisah", "", editorial_notes.strip(), ""]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results_v6")
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument(
        "--diagnostic-report",
        action="store_true",
        help="Write a clearly marked diagnostic Markdown inside the results directory",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    results = (root / args.results_dir).resolve()
    report = json.loads((results / "benchmark_verification_summary.json").read_text())
    if report.get("status") != "completed":
        raise ValueError("completed experiment required")
    manifest = json.loads((results / "run_manifest.json").read_text())
    if manifest.get("status") != "completed" or manifest.get(
        "source_sha256"
    ) != report.get("source_sha256"):
        raise ValueError("Run manifest is incomplete or inconsistent with the report")
    from run_benchmark import code_hashes

    if report.get("render_source_sha256", report["source_sha256"]) != code_hashes():
        raise ValueError(
            "source hashes differ from completed run; rerun before publishing"
        )
    notes_path = results / "DIAGNOSTIC_EDITORIAL_NOTES.md"
    notes = notes_path.read_text() if notes_path.exists() else ""
    if args.diagnostic_report:
        path = results / "DIAGNOSTIC_VERIFICATION_REPORT.md"
        path.write_text(build_diagnostic_report(report, notes))
        print(path)
        return
    from research_quality import require_research_readiness

    require_research_readiness(report)
    md = build_manuscript(report, results, root)
    (root / (NAME + ".md")).write_text(md)
    (results / "DIAGNOSTIC_VERIFICATION_REPORT.md").write_text(
        build_diagnostic_report(report, notes)
    )
    if args.no_pdf:
        return
    pandoc = shutil.which("pandoc")
    latex = shutil.which("pdflatex")
    if not pandoc or not latex:
        raise RuntimeError("pandoc and pdflatex must be on PATH")
    header = root / "tmp" / "pdfs" / "manuscript_header.tex"
    header.parent.mkdir(parents=True, exist_ok=True)
    header.write_text("\\usepackage{pdflscape}\n")
    subprocess.run(
        [
            pandoc,
            str(root / (NAME + ".md")),
            "--standalone",
            "--toc",
            "--variable",
            "papersize=a4",
            "--variable",
            "geometry:margin=22mm",
            "--variable",
            "fontsize=10pt",
            "--variable",
            "colorlinks=true",
            "--include-in-header",
            str(header),
            "--resource-path",
            str(root),
            "-o",
            str(root / (NAME + ".tex")),
        ],
        check=True,
        cwd=root,
    )
    # Old class/language auxiliaries are incompatible with a new Pandoc document.
    for suffix in [".aux", ".toc", ".out"]:
        (root / (NAME + suffix)).unlink(missing_ok=True)
    for _ in range(3):
        process = subprocess.run(
            [latex, "-interaction=nonstopmode", "-halt-on-error", NAME + ".tex"],
            cwd=root,
            capture_output=True,
            text=True,
        )
        if process.returncode:
            raise RuntimeError(process.stdout[-4000:])
    print(
        "Generated Markdown, LaTeX and PDF from completed run:", root / (NAME + ".pdf")
    )


if __name__ == "__main__":
    main()
