"""Reproducible synthetic experiment driver. --quick runs all paths at reduced budget.

A completed run is not evidence of posterior calibration or field readiness.
All figures, scalar metrics, samples, input data and source hashes are persisted.
"""

from __future__ import annotations
import importlib.metadata
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib, json, os, platform, sys, time
from pathlib import Path
import numpy as np
from models import get_benchmark_a_lvl, get_benchmark_b_water_table
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
    NUMBA_AVAILABLE,
    DISBA_AVAILABLE,
)
from fusion import ActivePassiveSpectralFuser, BroadbandFusionConfig
from smvi_inversion import (
    SMVIConfig,
    PickingFreeSpectralLikelihood,
    SteinMixtureVariationalInference,
)
from acquisition_design import (
    plot_acquisition_design_comparison,
    plot_dispersion_image_acquisition_comparison,
)
from ablation_experiments import (
    run_phd_ablation_suite,
    run_noise_and_weak_source_stress_test,
    save_inference,
    profile_metrics,
)
from visualize import (
    plot_waveforms_and_fusion,
    plot_fc_spectrum_comparison,
    plot_posterior_profiles_and_marginals,
)

VERSION = "v6_followup_fixes"


def code_hashes():
    root = Path(__file__).resolve().parent
    paths = list(root.glob("*.py")) + list((root / "tests").glob("*.py"))
    paths += [
        root / name
        for name in ["requirements.txt", "requirements-lock.txt", "environment.yml"]
    ]
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(paths)
    }


def model_summary(model, space, result):
    samples = result["posterior_samples_phys"]
    lo, hi = np.percentile(samples, [2.5, 97.5], axis=0)
    truth = space.true_vector
    mean = samples.mean(axis=0)
    std = samples.std(axis=0)
    params = {
        name: {
            "true": float(truth[i]),
            "posterior_mean": float(mean[i]),
            "posterior_std": float(std[i]),
            "best_component": float(result["best_component_phys"][i]),
            "ci95_low": float(lo[i]),
            "ci95_high": float(hi[i]),
            "rel_error_pct": float(100 * abs(mean[i] - truth[i]) / abs(truth[i])),
            "covered_by_95ci": bool(lo[i] <= truth[i] <= hi[i]),
        }
        for i, name in enumerate(space.param_names)
    }
    covered = int(np.sum((truth >= lo) & (truth <= hi)))
    return {
        "benchmark_name": model.name,
        "config": result["config"],
        "parameters": params,
        "mean_model_metrics": result["mean_model_metrics"],
        "best_component_metrics": result["best_component_metrics"],
        "predictive_metrics": result.get("predictive_metrics"),
        "optimization_diagnostics": result["optimization_diagnostics"],
        "baseline_comparison": result.get("baseline_comparison"),
        "empirical_parameter_coverage": f"{covered}/{space.dim}",
        "coverage_caveat": "Coverage across parameters of one synthetic model is not repeated-data calibration.",
        "forward_evaluations": int(result["n_forward_evals"][0]),
        "elapsed_sec": float(result["elapsed_sec"][0]),
        **profile_metrics(samples, space),
    }


def run_single_benchmark(
    true_model,
    param_space,
    c_grid,
    output_dir,
    seed=2026,
    config=None,
    freqs_active=None,
    freqs_passive=None,
    quadrature_density=2.0,
):
    config = replace(config or SMVIConfig(), random_seed=seed)
    if config.wasserstein_scale_mps is None:
        config = replace(
            config, wasserstein_scale_mps=float(0.15 * (c_grid[-1] - c_grid[0]))
        )
    fa = np.linspace(8, 56, 33) if freqs_active is None else freqs_active
    fp = np.linspace(0.5, 12, 24) if freqs_passive is None else freqs_passive
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(n_radial_bins=0),
        fa,
        fp,
        c_grid,
        quadrature_density,
    )
    fuser = ActivePassiveSpectralFuser(op, BroadbandFusionConfig(norm_power=1.35))
    rng = np.random.default_rng(seed)
    obs = fuser.generate_fused_spectrum_for_model(true_model, add_noise=True, rng=rng)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    prefix = output / true_model.code
    np.savez_compressed(
        str(prefix) + "_observations.npz",
        **obs,
        frequencies=fuser.freqs_broadband,
        freqs_active=fa,
        freqs_passive=fp,
        c_grid=c_grid,
        active_offsets=op.r_act,
        passive_offsets=op.r_pas,
        station_coordinates=op.passive_geom.station_xy,
        prior_lower=param_space.lower_bounds,
        prior_upper=param_space.upper_bounds,
        true_model=param_space.true_vector,
    )
    print("[Figures] " + true_model.code, flush=True)
    ta, ra, sa = op.synthesize_time_gather(true_model, "active", True, rng)
    tp, rp, sp = op.synthesize_time_gather(true_model, "passive", True, rng)
    plot_waveforms_and_fusion(
        true_model,
        op,
        fuser,
        obs,
        ta,
        ra,
        sa,
        tp,
        rp,
        sp,
        str(prefix) + "_fig1_waveforms_and_fusion.png",
    )
    likelihood = PickingFreeSpectralLikelihood(
        fuser, param_space, obs["d_fused"], config
    )
    print(f"[Benchmark] {true_model.code}", flush=True)
    res = SteinMixtureVariationalInference(likelihood, config).run(verbose=True)
    res["config"] = asdict(config)
    prior = param_space.latent_to_physical(np.zeros(param_space.dim))
    _, jp, np_, wp, _ = likelihood.evaluate_physical_model(prior)
    _, jt, nt, wt, _ = likelihood.evaluate_physical_model(param_space.true_vector)
    prior_profile = profile_metrics(prior[None, :], param_space)
    posterior_profile = profile_metrics(res["posterior_samples_phys"], param_space)
    warnings = []
    if res["mean_model_metrics"]["misfit"] > jp:
        warnings.append("mean_model_fit_worse_than_prior_midpoint")
    predictive = res.get("predictive_metrics")
    if predictive:
        if predictive["misfit"] > jp:
            warnings.append("predictive_point_estimate_worse_than_prior_midpoint")
        if not predictive["precision_target_met"]:
            warnings.append("predictive_monte_carlo_precision_target_not_met")
        if abs(predictive["misfit"] - jp) <= 2 * predictive["misfit_mc_se"]:
            warnings.append("predictive_prior_difference_not_resolved_by_mc_precision")
    if (
        posterior_profile["vs_profile_relative_rmse_pct"]
        > prior_profile["vs_profile_relative_rmse_pct"]
    ):
        warnings.append("posterior_profile_RMSE_worse_than_informative_prior_midpoint")
    res["prior_midpoint_phys"] = prior
    res["baseline_comparison"] = {
        "prior_midpoint": {
            "misfit": float(jp),
            "ncc_similarity": float(1 - np_),
            "wasserstein_w2_mps": float(wp),
            **prior_profile,
        },
        "truth_model_evaluation_only": {
            "misfit": float(jt),
            "ncc_similarity": float(1 - nt),
            "wasserstein_w2_mps": float(wt),
        },
        "warnings": warnings,
        "validation_forward_evaluations": 2,
        "scope": "A prior midpoint is a point baseline, not an independently calibrated posterior.",
    }
    save_inference(str(prefix) + "_inference.npz", res)
    plot_fc_spectrum_comparison(
        true_model,
        fuser,
        obs["d_fused"],
        res,
        str(prefix) + "_fig2_fc_spectrum_comparison.png",
    )
    plot_posterior_profiles_and_marginals(
        true_model,
        param_space,
        res,
        str(prefix) + "_fig3_posterior_profiles_and_marginals.png",
    )
    return model_summary(true_model, param_space, res), fuser, obs, res


def quadrature_check(model, operator):
    base = operator.quadrature_density
    checks = {}
    for name, f, r in [
        ("active", np.array([8.0, 20.0, 40.0]), np.array([4.0, 20.0, 60.0, 114.0])),
        (
            "passive",
            np.array([1.0, 4.0, 8.0, 12.0]),
            np.array([10.92, 100.0, 300.0, 562.63]),
        ),
    ]:
        operator.quadrature_density = base
        a = operator._green_gather(model, f, r)
        operator.quadrature_density = 2 * base
        b = operator._green_gather(model, f, r)
        relative = abs(a - b) / np.maximum(abs(b), 1e-30)
        checks[name] = {
            "frequencies_hz": f.tolist(),
            "offsets_m": r.tolist(),
            "relative_l2_change": float(
                np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30)
            ),
            "maximum_relative_point_change": float(np.max(relative)),
        }
    operator.quadrature_density = base
    return {
        "densities": [base, 2 * base],
        "checks": checks,
        "relative_complex_gather_l2_change": max(
            v["relative_l2_change"] for v in checks.values()
        ),
        "status": (
            "passed_5pct_selected_active_and_passive_offsets"
            if max(v["maximum_relative_point_change"] for v in checks.values()) < 0.05
            else "requires_quadrature_refinement"
        ),
        "caveat": "Selected frequencies/offsets only; this does not establish convergence of every model in the prior.",
    }


def benchmark_frequency_grids(quick=False):
    # The overlap is synthesized at identical native nodes in both channels.
    fa = np.unique(
        np.r_[
            np.arange(8.0, 55.01, 2.0 if quick else 1.0),
            np.arange(8.0, 12.01, 0.5),
            55.0,
        ]
    )
    fp = np.arange(0.5, 12.01, 0.5)
    return fa, fp


def _benchmark_worker(*args):
    """Independent numerical process; return small summaries, not Green caches."""
    summary, fuser, _, _ = run_single_benchmark(*args)
    return summary, quadrature_check(args[0], fuser.operator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Reduced-budget end-to-end verification, not production inference",
    )
    parser.add_argument("--output-dir", default="results_v6")
    parser.add_argument(
        "--benchmarks-only",
        action="store_true",
        help="Save two benchmark inversions; do not claim ablations were run",
    )
    parser.add_argument(
        "--parallel-benchmarks",
        action="store_true",
        help="Run the two independent benchmark processes concurrently (benchmarks-only)",
    )
    parser.add_argument("--predictive-draws", type=int, default=128)
    parser.add_argument("--predictive-max-draws", type=int, default=512)
    parser.add_argument("--predictive-bootstrap-repeats", type=int, default=32)
    parser.add_argument("--predictive-mc-tolerance", type=float, default=0.01)
    parser.add_argument("--mc-samples", type=int, default=4)
    parser.add_argument("--components", type=int, default=4)
    parser.add_argument("--refinement-steps", type=int, default=25)
    parser.add_argument(
        "--gradient-backend",
        choices=["autodiff", "finite_difference"],
        default="autodiff",
    )
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    if args.parallel_benchmarks and not args.benchmarks_only:
        parser.error("--parallel-benchmarks currently requires --benchmarks-only")
    root = Path(__file__).resolve().parent
    output = Path(args.output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(
            "Output directory must be empty; existing experiments are preserved. Use a new --output-dir."
        )
    output.mkdir(parents=True, exist_ok=True)
    quick = args.quick
    cfg = SMVIConfig(
        n_components=args.components,
        n_iterations=args.iterations or (8 if quick else 120),
        n_mc_samples_per_comp=args.mc_samples,
        n_posterior_draws=1200,
        random_seed=args.seed,
        gradient_backend=args.gradient_backend,
        initialization_refinement_steps=args.refinement_steps,
        predictive_draws=args.predictive_draws,
        predictive_max_draws=args.predictive_max_draws,
        predictive_bootstrap_repeats=args.predictive_bootstrap_repeats,
        predictive_mc_tolerance=args.predictive_mc_tolerance,
    )
    fa, fp = benchmark_frequency_grids(quick)
    density = 2.0
    start = time.time()
    manifest = {
        "version": VERSION,
        "status": "running",
        "run_started_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "quick_verification" if quick else "research_experiment",
        "parallel_benchmarks": args.parallel_benchmarks,
        "config": asdict(cfg),
        "frequency_grids": {"active": fa.tolist(), "passive": fp.tolist()},
        "quadrature_density": density,
        "source_sha256": code_hashes(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numba": NUMBA_AVAILABLE,
            "disba": DISBA_AVAILABLE,
            "packages": {
                name: importlib.metadata.version(name)
                for name in [
                    "numpy",
                    "scipy",
                    "numba",
                    "matplotlib",
                    "disba",
                    "torch",
                    "pyprop8",
                ]
            },
        },
        "limitations": [
            "Synthetic generation and inversion share the Green solver.",
            "Passive EGF noise is phenomenological; no ambient-source interferometry simulation.",
            "SMVI and its intervals are experimental and not calibrated.",
            "No global MAP, independent field validation, or speedup claim is made.",
        ],
    }
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    try:
        geometry = plot_acquisition_design_comparison(
            str(output / "fig0_field_acquisition_design_comparison.png")
        )
        ma, pa = get_benchmark_a_lvl()
        mb, pb = get_benchmark_b_water_table()
        ca = np.arange(120.0, 1250.01, 4.0 if quick else 2.0)
        cb = np.arange(130.0, 2450.01, 4.0 if quick else 2.0)
        images = plot_dispersion_image_acquisition_comparison(
            ma,
            ca,
            str(output / "fig0b_dispersion_image_comparison.png"),
            fa,
            fp,
            density,
        )
        if args.parallel_benchmarks:
            from concurrent.futures import ProcessPoolExecutor
            import multiprocessing

            with ProcessPoolExecutor(
                max_workers=2, mp_context=multiprocessing.get_context("spawn")
            ) as pool:
                af = pool.submit(
                    _benchmark_worker,
                    ma,
                    pa,
                    ca,
                    output,
                    args.seed,
                    cfg,
                    fa,
                    fp,
                    density,
                )
                bf = pool.submit(
                    _benchmark_worker,
                    mb,
                    pb,
                    cb,
                    output,
                    args.seed + 1,
                    cfg,
                    fa,
                    fp,
                    density,
                )
                a, qa = af.result()
                b, qb = bf.result()
        else:
            a, fua, oa, ra = run_single_benchmark(
                ma, pa, ca, output, args.seed, cfg, fa, fp, density
            )
            b, fub, ob, rb = run_single_benchmark(
                mb, pb, cb, output, args.seed + 1, cfg, fa, fp, density
            )
            qa = quadrature_check(ma, fua.operator)
            qb = quadrature_check(mb, fub.operator)
        ablations = stress = {
            "status": "not_run",
            "reason": "benchmarks-only diagnostic run",
        }
        if not args.benchmarks_only:
            abl_cfg = replace(
                cfg,
                n_iterations=8 if quick else 60,
                n_components=3 if quick else 4,
                initialization_refinement_steps=10,
            )
            ablations = run_phd_ablation_suite(
                ma,
                pa,
                fua,
                oa["d_fused"],
                ra,
                mb,
                pb,
                fub,
                ob["d_fused"],
                rb,
                str(output / "fig4_ablation_and_mcmc.png"),
                abl_cfg,
                mcmc_steps=80 if quick else 1000,
                burn_in=25 if quick else 300,
                obs_fusion_a=oa,
            )
            stress = run_noise_and_weak_source_stress_test(
                ma,
                pa,
                fua,
                oa["d_fused"],
                ra,
                str(output / "fig5_waveform_noise_stress.png"),
                abl_cfg,
            )
        report = {
            **manifest,
            "status": "completed",
            "total_runtime_sec": time.time() - start,
            "benchmark_a": a,
            "benchmark_b": b,
            "field_acquisition_evaluation": geometry,
            "dispersion_image_comparison": images,
            "ablation_studies": ablations,
            "waveform_noise_stress_test": stress,
            "quadrature_check": qa,
            "quadrature_check_b": qb,
        }
        from research_quality import assess_research_readiness

        report["research_readiness"] = assess_research_readiness(report)
        if code_hashes() != manifest["source_sha256"]:
            raise RuntimeError(
                "Source changed during run; do not publish mixed-version output"
            )
        temp = output / "benchmark_verification_summary.json.tmp"
        temp.write_text(json.dumps(report, indent=2, allow_nan=False))
        temp.replace(output / "benchmark_verification_summary.json")
        manifest.update(
            status="completed", total_runtime_sec=report["total_runtime_sec"]
        )
        print(f'[Done] {report["total_runtime_sec"]:.1f} seconds; {output}', flush=True)
    except KeyboardInterrupt:
        manifest.update(status="interrupted", error="Interrupted before completion")
        raise
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
