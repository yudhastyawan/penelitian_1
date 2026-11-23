"""Separate successful execution from scientific acceptance/publication."""


def assess_research_readiness(report):
    failures = []
    if report.get("status") != "completed":
        failures.append("experiment_incomplete")
    if report.get("mode") != "research_experiment":
        failures.append("diagnostic_or_quick_run")
    for key in ["benchmark_a", "benchmark_b"]:
        if (
            not report.get(key, {})
            .get("optimization_diagnostics", {})
            .get("converged", False)
        ):
            failures.append(key + "_optimization_not_validated")
    for key in ["benchmark_a", "benchmark_b"]:
        if not (report.get(key, {}).get("predictive_metrics") or {}).get(
            "precision_target_met", False
        ):
            failures.append(key + "_predictive_mc_precision_not_validated")
    for key in ["quadrature_check", "quadrature_check_b"]:
        if (
            report.get(key, {}).get("status")
            != "passed_5pct_selected_active_and_passive_offsets"
        ):
            failures.append(key + "_failed_or_missing")
    mcmc = report.get("ablation_studies", {}).get("mcmc", {})
    if mcmc.get("convergence_status") != "validated_independent_chains":
        failures.append("independent_posterior_reference_missing")
    for name in [
        "frequency_velocity_resolution",
        "repeated_data_calibration",
        "independent_waveform_validation",
    ]:
        if report.get("validation", {}).get(name, {}).get("status") != "passed":
            failures.append(name + "_missing_or_failed")
    return {
        "publishable": not failures,
        "blocking_checks": failures,
        "meaning": "Computation completion is distinct from scientific validation. No automatic override.",
    }


def require_research_readiness(report):
    assessment = assess_research_readiness(report)
    if not assessment["publishable"]:
        raise ValueError(
            "Research manuscript blocked: " + ", ".join(assessment["blocking_checks"])
        )
    return assessment
