"""Controlled inversions and waveform-level noise experiments.

Every numerical comparison comes from its own forward/likelihood calculation.
Posterior samples are never shifted to enforce an expected scientific outcome.
"""

from dataclasses import replace, asdict
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from models import evaluate_depth_profiles
from smvi_inversion import (
    SMVIConfig,
    PickingFreeSpectralLikelihood,
    SteinMixtureVariationalInference,
    AffineInvariantEnsembleMCMC,
)


def profile_metrics(samples, space):
    z = np.linspace(0, 32, 250)
    actual = evaluate_depth_profiles(space.true_vector[None, :], space, z)["vs_mean"]
    estimated = evaluate_depth_profiles(samples, space, z)["vs_mean"]
    relative = (estimated - actual) / actual
    return {
        "vs_profile_relative_rmse_pct": float(100 * np.sqrt(np.mean(relative**2))),
        "vs_profile_mape_pct": float(100 * np.mean(abs(relative))),
    }


def fit_case(
    fuser,
    space,
    observed,
    config,
    frequency_weights=None,
    velocity_max=None,
    channel="fused",
):
    lik = PickingFreeSpectralLikelihood(
        fuser, space, observed, config, velocity_max, channel
    )
    if frequency_weights is not None:
        lik.set_frequency_weights(frequency_weights)
    result = SteinMixtureVariationalInference(lik, config).run(verbose=False)
    result["config"] = {
        **asdict(config),
        "wasserstein_scale_mps": lik.wasserstein_scale_mps,
    }
    return result


def save_inference(path, result):
    arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
    np.savez_compressed(path, **arrays)
    if "config" in result:
        import json

        Path(str(path) + ".json").write_text(json.dumps(result["config"], indent=2))


def summary_case(result, space):
    samples = result["posterior_samples_phys"]
    return {
        **profile_metrics(samples, space),
        "config": result.get("config"),
        "posterior_mean": samples.mean(axis=0).tolist(),
        "posterior_std": samples.std(axis=0).tolist(),
        "mean_model_metrics": result["mean_model_metrics"],
        "forward_evaluations": int(result["n_forward_evals"][0]),
        "elapsed_sec": float(result["elapsed_sec"][0]),
    }


def run_phd_ablation_suite(
    model_a,
    space_a,
    fuser_a,
    d_obs_a,
    smvi_res_a,
    model_b,
    space_b,
    fuser_b,
    d_obs_b,
    smvi_res_b,
    output_fig_path,
    config=None,
    mcmc_steps=1000,
    burn_in=300,
    obs_fusion_a=None,
):
    config = config or SMVIConfig(n_components=4, n_iterations=30)
    output = Path(output_fig_path).parent
    print("[Ablation] Active-only and passive-only inversions", flush=True)
    # Use each channel's observed spectrum, not a frequency-masked fused overlap.
    if obs_fusion_a is None:
        raise ValueError(
            "Observed active/passive spectra are required for controlled ablations"
        )
    f = fuser_a.freqs_broadband
    cfg = replace(config, random_seed=3001)
    active = fit_case(
        fuser_a, space_a, obs_fusion_a["d_active_bb"], cfg, f >= 8, channel="active_bb"
    )
    passive = fit_case(
        fuser_a,
        space_a,
        obs_fusion_a["d_passive_bb"],
        replace(cfg, random_seed=3002),
        f <= 12,
        channel="passive_bb",
    )
    print(
        "[Ablation] Rayleigh-window inversion and optimized single Gaussian", flush=True
    )
    ray = fit_case(
        fuser_b, space_b, d_obs_b, replace(cfg, random_seed=3003), velocity_max=600.0
    )
    advi = fit_case(
        fuser_a,
        space_a,
        d_obs_a,
        replace(cfg, n_components=1, repulsion_strength=0.0, random_seed=3004),
    )
    print(
        "[Ablation] Independent-start reference MCMC (convergence is diagnosed)",
        flush=True,
    )
    mlik = PickingFreeSpectralLikelihood(fuser_a, space_a, d_obs_a, cfg)
    mcmc = AffineInvariantEnsembleMCMC(
        mlik, n_walkers=24, n_steps=mcmc_steps, burn_in=burn_in, seed=4001
    ).run()
    for name, result in [
        ("active_only", active),
        ("passive_only", passive),
        ("rayleigh_only", ray),
        ("single_gaussian", advi),
        ("mcmc", mcmc),
    ]:
        save_inference(output / (name + "_samples.npz"), result)
    # Conventional picked-curve baseline is an independent constrained optimizer.
    picked = {
        "status": "not_available",
        "reason": "disba fundamental solver not installed",
    }
    try:
        from disba import PhaseDispersion
        from scipy.optimize import least_squares

        mask = fuser_a.c_grid <= 600
        curve = fuser_a.c_grid[mask][np.argmax(d_obs_a[:, mask], axis=1)]
        periods = 1 / f[::-1]
        curve_evaluations = [0]

        def residual(z):
            curve_evaluations[0] += 1
            mod = space_a.vector_to_model(space_a.latent_to_physical(z))
            try:
                res = PhaseDispersion(
                    mod.thickness / 1000, mod.vp / 1000, mod.vs / 1000, mod.rho / 1000
                )(periods, mode=0, wave="rayleigh")
                if len(res.period) < 3:
                    return np.full(len(f), 1000.0)
                pred = np.interp(f, (1 / res.period)[::-1], res.velocity[::-1] * 1000)
                return (pred - curve) / 20
            except Exception:
                return np.full(len(f), 1000.0)

        opt = least_squares(
            residual,
            np.zeros(space_a.dim),
            max_nfev=35 if config.n_iterations < 20 else 150,
        )
        picked = {
            "status": "optimized_point_estimate",
            "optimizer_success": bool(opt.success),
            "model": space_a.latent_to_physical(opt.x).tolist(),
            "picked_curve_rmse_mps": float(
                np.sqrt(np.mean((residual(opt.x) * 20) ** 2))
            ),
            "posterior_available": False,
            "curve_evaluations": curve_evaluations[0],
            "nfev": int(opt.nfev),
            "seed_initialization": "prior midpoint",
        }
        np.savez_compressed(
            output / "picked_curve_fit.npz",
            frequencies=f,
            picked_velocity=curve,
            model=space_a.latent_to_physical(opt.x),
        )
    except ImportError:
        pass
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    z = np.linspace(0, 32, 250)
    for label, result in [
        ("Active only", active),
        ("Passive only", passive),
        ("Joint", smvi_res_a),
    ]:
        prof = evaluate_depth_profiles(result["posterior_samples_phys"], space_a, z)
        axes[0].plot(prof["vs_mean"], z, label=label)
    true = evaluate_depth_profiles(space_a.true_vector[None, :], space_a, z)
    axes[0].plot(true["vs_mean"], z, "k--", label="True")
    axes[0].invert_yaxis()
    axes[0].set(xlabel="Vs (m/s)", ylabel="Depth (m)", title="Independent inversions")
    axes[0].legend(fontsize=8)
    axes[1].hist(
        ray["posterior_samples_phys"][:, 5],
        bins=25,
        density=True,
        alpha=0.5,
        label="Rayleigh window",
    )
    axes[1].hist(
        smvi_res_b["posterior_samples_phys"][:, 5],
        bins=25,
        density=True,
        alpha=0.5,
        label="Full window",
    )
    axes[1].set(xlabel="Vp2 (m/s)", title="Estimated uncertainty, not calibrated")
    axes[1].legend(fontsize=8)
    for name, result in [
        ("SMVI", smvi_res_a),
        ("Single Gaussian", advi),
        ("MCMC", mcmc),
    ]:
        key = "samples_phys" if name == "MCMC" else "posterior_samples_phys"
        axes[2].hist(
            result[key][:, 1], bins=25, density=True, histtype="step", label=name
        )
    axes[2].set(xlabel="Vs2 (m/s)", title="MCMC: " + mcmc["convergence_status"])
    axes[2].legend(fontsize=8)
    fig.savefig(output_fig_path, dpi=160)
    plt.close(fig)
    return {
        "active_only": summary_case(active, space_a),
        "passive_only": summary_case(passive, space_a),
        "joint": summary_case(smvi_res_a, space_a),
        "rayleigh_only": summary_case(ray, space_b),
        "single_gaussian": summary_case(advi, space_a),
        "picked_curve": picked,
        "mcmc": {
            "n_steps": mcmc_steps,
            "burn_in": burn_in,
            "n_walkers": 24,
            "random_seed": 4001,
            "convergence_status": mcmc["convergence_status"],
            "converged": mcmc["converged"],
            "acceptance_rate": float(mcmc["acceptance_rate"][0]),
            "autocorrelation_time_steps": mcmc["autocorrelation_time_steps"].tolist(),
            "ess": mcmc["ess"].tolist(),
            "forward_evaluations": int(mcmc["n_evals"][0]),
            "elapsed_sec": float(mcmc["elapsed_sec"][0]),
            "vs2_mean": float(mcmc["samples_phys"][:, 1].mean()),
            "vs2_std": float(mcmc["samples_phys"][:, 1].std()),
        },
        "figure_path": str(output_fig_path),
    }


def add_complex_noise(signal, snr_db, rng):
    """Calibrate additive noise power exactly on the supplied complex gather."""
    noise = rng.normal(size=signal.shape) + 1j * rng.normal(size=signal.shape)
    noise *= np.linalg.norm(signal) / (np.linalg.norm(noise) * 10 ** (snr_db / 20))
    measured = 20 * np.log10(np.linalg.norm(signal) / np.linalg.norm(noise))
    return signal + noise, float(measured)


def stress_observation(
    fuser, model, weak_source=False, snr_db=18.0, coupling_sigma=0.35, seed=2030
):
    op = fuser.operator
    rng = np.random.default_rng(seed)
    ua, _ = op.synthesize_frequency_gather(model, "active", False)
    up, _ = op.synthesize_frequency_gather(model, "passive", False)
    original_rms = np.linalg.norm(ua)
    if weak_source:
        ua *= 1 / (1 + np.exp(-(op.freqs_active[:, None] - 15) / 2.0))
        ua *= np.exp(-op.r_act[None, :] / 26.0)
    # Apply receiver coupling before imaging, then calibrate additive noise.
    ua *= rng.lognormal(-0.5 * coupling_sigma**2, coupling_sigma, (1, len(op.r_act)))
    ua *= np.exp(1j * rng.normal(0, 0.08, (1, len(op.r_act))))
    if weak_source:
        # Fixed background relative to the original strong source; no invented spectral gap.
        norm = np.linalg.norm(ua)
        effective_snr = snr_db + 20 * np.log10(max(norm, 1e-30) / original_rms)
    else:
        effective_snr = snr_db
    ua, sa = add_complex_noise(ua, effective_snr, rng)
    up, sp = add_complex_noise(up, snr_db, rng)
    fused = fuser.fuse_spectra(
        op.compute_mfj_spectrum(ua, "active"), op.compute_mfj_spectrum(up, "passive")
    )
    fused.update(
        measured_active_snr_db=sa,
        measured_passive_snr_db=sp,
        coupling_lognormal_sigma=coupling_sigma,
        weak_source=weak_source,
    )
    return fused


def run_noise_and_weak_source_stress_test(
    model_a, space_a, fuser_a, d_obs_ideal, smvi_res_ideal, output_fig_path, config=None
):
    from acquisition_design import build_four_candidate_field_designs, make_operator
    from fusion import ActivePassiveSpectralFuser

    config = config or SMVIConfig(n_components=4, n_iterations=30)
    base = fuser_a.operator
    op1 = make_operator(
        build_four_candidate_field_designs()[0],
        base.freqs_active,
        base.freqs_passive,
        base.c_grid,
        base.quadrature_density,
    )
    fuser1 = ActivePassiveSpectralFuser(op1, fuser_a.config)
    cases = [
        ("weak_source_option1", fuser1, True, 18.0, 0.35),
        ("weak_source_option4", fuser_a, True, 18.0, 0.35),
        ("high_noise_option4", fuser_a, False, 2.0, 0.65),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    report = {}
    for i, (name, fuser, weak, snr, coupling) in enumerate(cases):
        print("[Stress] " + name, flush=True)
        obs = stress_observation(fuser, model_a, weak, snr, coupling, seed=5010)
        res = fit_case(
            fuser, space_a, obs["d_fused"], replace(config, random_seed=5010)
        )
        save_inference(Path(output_fig_path).parent / (name + "_samples.npz"), res)
        np.savez_compressed(
            Path(output_fig_path).parent / (name + "_observation.npz"),
            spectrum=obs["d_fused"],
            active_weight=obs["w_active"],
            passive_weight=obs["w_passive"],
        )
        axes[0, i].pcolormesh(
            fuser.freqs_broadband,
            fuser.c_grid,
            obs["d_fused"].T,
            shading="auto",
            cmap="turbo",
        )
        axes[0, i].set(title=name, xlabel="Frequency (Hz)", ylabel="Velocity (m/s)")
        axes[1, i].plot(
            fuser.freqs_broadband, obs["w_active"], label="Actual active weight"
        )
        axes[1, i].plot(
            fuser.freqs_broadband, obs["w_passive"], label="Actual passive weight"
        )
        axes[1, i].set(xlabel="Frequency (Hz)", ylabel="Fusion weight")
        axes[1, i].legend(fontsize=7)
        report[name] = {
            **summary_case(res, space_a),
            "measured_active_snr_db": obs["measured_active_snr_db"],
            "measured_passive_snr_db": obs["measured_passive_snr_db"],
            "coupling_lognormal_sigma": coupling,
            "weak_source": weak,
            "data_seed": 5010,
            "nominal_snr_db": snr,
        }
    fig.savefig(output_fig_path, dpi=160)
    plt.close(fig)
    report["figure_path"] = str(output_fig_path)
    report["limitations"] = (
        "Same forward solver for generation/inversion; no field validation or preselected winner."
    )
    return report
