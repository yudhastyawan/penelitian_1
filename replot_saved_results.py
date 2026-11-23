"""Refresh figures from saved inference arrays without rerunning optimization.

Reject changed numerical code. Keep the original experiment hashes immutable,
and record separate rendering hashes for provenance.
"""

import argparse, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from run_benchmark import code_hashes
from models import get_benchmark_a_lvl, get_benchmark_b_water_table
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
)
from fusion import ActivePassiveSpectralFuser, BroadbandFusionConfig
from visualize import (
    plot_waveforms_and_fusion,
    plot_fc_spectrum_comparison,
    plot_posterior_profiles_and_marginals,
)

RENDER_FILES = {"visualize.py", "generate_pdf_manuscript.py", "replot_saved_results.py"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results_v6")
    parser.add_argument(
        "--document-only",
        action="store_true",
        help="Refresh document provenance only; require unchanged visualization source",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    out = (root / args.results_dir).resolve()
    path = out / "benchmark_verification_summary.json"
    report = json.loads(path.read_text())
    if report.get("status") != "completed":
        raise ValueError("Completed experiment required")
    original = report["source_sha256"]
    current = code_hashes()
    for name, sha in original.items():
        if name not in RENDER_FILES and current.get(name) != sha:
            raise ValueError(f"Numerical source changed: {name}. Rerun the experiment.")
    if set(current) - set(original) - RENDER_FILES:
        raise ValueError(
            "Unrecognized new source files; inspect provenance before rendering"
        )
    if args.document_only:
        previous = report.get("render_source_sha256", original)
        if previous.get("visualize.py") != current.get("visualize.py"):
            raise ValueError("Visualization source changed; refresh figures first")
    else:
        for code, getmodel in [
            ("benchmark_a", get_benchmark_a_lvl),
            ("benchmark_b", get_benchmark_b_water_table),
        ]:
            model, space = getmodel()
            with np.load(out / (code + "_observations.npz")) as data:
                observed = {k: data[k] for k in data.files}
            with np.load(out / (code + "_inference.npz")) as data:
                result = {k: data[k] for k in data.files}
            operator = FullWavefieldMFJOperator(
                ActiveSurveyGeometry(),
                PassiveSurveyGeometry(n_radial_bins=0),
                observed["freqs_active"],
                observed["freqs_passive"],
                observed["c_grid"],
                report["quadrature_density"],
            )
            fuser = ActivePassiveSpectralFuser(
                operator, BroadbandFusionConfig(norm_power=1.35)
            )
            # Recreate the original random-stream position before representative time gathers.
            rng = np.random.default_rng(report[code]["config"]["random_seed"])
            replay = fuser.generate_fused_spectrum_for_model(model, True, rng)
            for key in ["u_active", "u_passive", "d_fused"]:
                np.testing.assert_allclose(
                    replay[key], observed[key], rtol=1e-12, atol=1e-30
                )
            ta, ra, sa = operator.synthesize_time_gather(model, "active", True, rng)
            tp, rp, sp = operator.synthesize_time_gather(model, "passive", True, rng)
            prefix = str(out / code)
            plot_waveforms_and_fusion(
                model,
                operator,
                fuser,
                observed,
                ta,
                ra,
                sa,
                tp,
                rp,
                sp,
                prefix + "_fig1_waveforms_and_fusion.png",
            )
            plot_fc_spectrum_comparison(
                model,
                fuser,
                observed["d_fused"],
                result,
                prefix + "_fig2_fc_spectrum_comparison.png",
            )
            plot_posterior_profiles_and_marginals(
                model,
                space,
                result,
                prefix + "_fig3_posterior_profiles_and_marginals.png",
            )
            print("Reproduced observations and refreshed figures:", code, flush=True)
    revision = {
        "rendered_utc": datetime.now(timezone.utc).isoformat(),
        "figures_replotted": not args.document_only,
        "note": "Only rendering sources changed; numerical modules verified against original experiment hashes. Metrics and posterior samples unchanged. Representative time gathers are separate fine-grid realizations, not the sparse frequency data fitted by inversion.",
    }
    report["render_source_sha256"] = current
    report["render_revision"] = revision
    path.write_text(json.dumps(report, indent=2, allow_nan=False))
    manifest_path = out / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["render_source_sha256"] = current
    manifest["render_revision"] = revision
    manifest_path.write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
