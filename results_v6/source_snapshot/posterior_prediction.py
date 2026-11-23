"""Balanced mixture predictions with explicit, bounded Monte Carlo precision.

Response spread and the standard error of its estimated mean are separate.
Bootstrap uncertainty here is numerical MC error conditional on fitted q; it
is not posterior calibration or measurement-noise uncertainty.
"""

import numpy as np


def evaluate_mixture_prediction(likelihood, means, cholesky, weights, config):
    weights = np.asarray(weights, dtype=float)
    if (
        not np.all(np.isfinite(weights))
        or np.any(weights < 0)
        or not np.isclose(weights.sum(), 1.0)
    ):
        raise ValueError("mixture weights must be finite probabilities")
    active = np.flatnonzero(weights > 0)
    k = len(active)
    dim = means.shape[1]
    minimum = int(config.predictive_draws)
    maximum = int(config.predictive_max_draws)
    repeats = int(config.predictive_bootstrap_repeats)
    tol = float(config.predictive_mc_tolerance)
    if (
        k < 1
        or minimum < 2 * k
        or maximum < minimum
        or repeats < 8
        or not np.isfinite(tol)
        or tol <= 0
    ):
        raise ValueError(
            "invalid predictive draw budget, MC tolerance or bootstrap count"
        )
    per = int(np.ceil(minimum / k))
    cap = maximum // k
    if per > cap:
        raise ValueError(
            "predictive budget cannot cover the requested balanced draw count"
        )
    rng = np.random.default_rng(config.random_seed + 100003)
    bootstrap_rng = np.random.default_rng(config.random_seed + 200003)
    banks = {i: [] for i in active}
    latents = {i: [] for i in active}
    history = []
    while True:
        for i in active:
            additional = per - len(banks[i])
            zs = means[i] + rng.normal(size=(additional, dim)) @ cholesky[i].T
            for z in zs:
                model = likelihood.param_space.latent_to_physical(z)
                image = np.asarray(likelihood.evaluate_physical_model(model)[-1])
                if not np.all(np.isfinite(image)):
                    raise FloatingPointError("nonfinite predictive spectrum")
                banks[i].append(image)
                latents[i].append(z)
        arrays = {i: np.asarray(banks[i]) for i in active}
        mean = sum(weights[i] * arrays[i].mean(axis=0) for i in active)
        variance_mean = sum(
            weights[i] ** 2 * arrays[i].var(axis=0, ddof=1) / per for i in active
        )
        mc_se = np.sqrt(variance_mean)
        response_var = np.maximum(
            sum(weights[i] * np.mean(arrays[i] ** 2, axis=0) for i in active) - mean**2,
            0.0,
        )
        jt, jn, jw = likelihood.spectral_metrics(mean)
        half = sum(weights[i] * arrays[i][: per // 2].mean(axis=0) for i in active)
        half_delta = abs(float(likelihood.spectral_metrics(half)[0]) - jt)
        boot = []
        for _ in range(repeats):
            pred = sum(
                weights[i] * arrays[i][bootstrap_rng.integers(0, per, per)].mean(axis=0)
                for i in active
            )
            boot.append(likelihood.spectral_metrics(pred)[0])
        metric_se = float(np.std(boot, ddof=1))
        image_se_rms = float(np.sqrt(np.mean(variance_mean)))
        stable = bool(
            metric_se <= tol and image_se_rms <= tol and half_delta <= 2 * tol
        )
        history.append(
            {
                "draws": int(k * per),
                "misfit": float(jt),
                "misfit_mc_se": metric_se,
                "mean_spectrum_mc_se_rms": image_se_rms,
                "nested_half_misfit_difference": float(half_delta),
                "precision_target_met": stable,
            }
        )
        if stable or per == cap:
            break
        per = min(2 * per, cap)
    labels = np.concatenate([np.full(per, i, dtype=int) for i in active])
    counts = np.bincount(labels, minlength=len(weights))
    return {
        "d_posterior_predictive_mean": mean,
        "d_posterior_predictive_std": np.sqrt(response_var),
        "d_posterior_predictive_mc_se": mc_se,
        "predictive_samples_z": np.concatenate([latents[i] for i in active]),
        "predictive_comp_labels": labels,
        "predictive_component_counts": counts,
        "predictive_metrics": {
            "misfit": float(jt),
            "ncc_similarity": float(1 - jn),
            "wasserstein_w2_mps": float(jw),
            "draws": int(len(labels)),
            "component_draw_counts": counts.tolist(),
            "misfit_mc_se": metric_se,
            "misfit_mc_resampling_interval95": np.percentile(
                boot, [2.5, 97.5]
            ).tolist(),
            "mean_spectrum_mc_se_rms": image_se_rms,
            "precision_target": tol,
            "precision_target_met": stable,
            "bootstrap_repeats": repeats,
            "sampling": "equal draws per component, component means combined with mixture weights",
            "scope": "Monte Carlo precision conditional on fitted mixture; not posterior calibration",
            "precision_history": history,
        },
    }
