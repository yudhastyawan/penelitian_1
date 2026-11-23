"""Regression tests for the five independently reproduced follow-up findings."""

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from smvi_inversion import (
    SMVIConfig,
    PickingFreeSpectralLikelihood,
    compute_spectral_misfit_numba,
)


def fake_fuser():
    f = np.r_[np.arange(1.0, 12.01, 0.5), np.arange(14.0, 54.01, 2.0), 55.0]
    c = np.arange(130.0, 2450.01, 4.0)
    return SimpleNamespace(
        freqs_broadband=f,
        c_grid=c,
        operator=SimpleNamespace(freqs_active=f[f >= 8], freqs_passive=f[f <= 12]),
    )


def test_transport_rejects_nonfinite_masses_without_hanging():
    # Isolate a regression of the former infinite loop and enforce a timeout.
    program = """
import numpy as np
from smvi_inversion import discrete_wasserstein2
q=np.ones(2);c=np.array([0.,1.])
discrete_wasserstein2(q,q,c)
for p in [np.array([np.nan,1.]),np.array([np.inf,1.]),np.array([-1.,2.]),np.zeros(2)]:
    try: discrete_wasserstein2(p,q,c)
    except ValueError: continue
    raise AssertionError('invalid masses were accepted')
print('rejected invalid masses')
"""
    r = subprocess.run(
        [sys.executable, "-c", program],
        capture_output=True,
        text=True,
        timeout=15,
        cwd=Path(__file__).resolve().parents[1],
    )
    assert r.returncode == 0, r.stderr
    assert "rejected invalid masses" in r.stdout


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0])
def test_likelihood_rejects_invalid_observations_before_forward(bad):
    f = fake_fuser()
    obs = np.ones((len(f.freqs_broadband), len(f.c_grid)))
    obs[0, 0] = bad
    with pytest.raises(ValueError, match="observations must be finite"):
        PickingFreeSpectralLikelihood(f, SimpleNamespace(), obs, SMVIConfig())


def test_fixed_wasserstein_scale_is_invariant_to_empty_velocity_extension():
    metrics = []
    for high in [598.0, 2450.0]:
        c = np.arange(130.0, high + 0.01, 4.0)
        a = np.exp(-0.5 * ((c - 250) / 12) ** 2)[None, :]
        b = np.exp(-0.5 * ((c - 300) / 12) ** 2)[None, :]
        metrics.append(
            compute_spectral_misfit_numba(
                a, b, c, np.ones(1), freqs=np.array([20.0]), wasserstein_scale_mps=200.0
            )
        )
    np.testing.assert_allclose(metrics[0], metrics[1], rtol=1e-8, atol=1e-8)
    f = fake_fuser()
    obs = np.ones((len(f.freqs_broadband), len(f.c_grid)))
    full = PickingFreeSpectralLikelihood(f, SimpleNamespace(), obs, SMVIConfig())
    restricted = PickingFreeSpectralLikelihood(
        f, SimpleNamespace(), obs, SMVIConfig(), velocity_max=600.0
    )
    assert full.wasserstein_scale_mps == restricted.wasserstein_scale_mps == 348.0


def test_channel_frequency_integral_and_boundary_weights():
    f = fake_fuser()
    obs = np.ones((len(f.freqs_broadband), len(f.c_grid)))
    lik = PickingFreeSpectralLikelihood(
        f, SimpleNamespace(), obs, SMVIConfig(), channel="passive_bb"
    )
    assert lik.freq_weights.sum() == pytest.approx(11.0)
    assert lik.freq_weights[f.freqs_broadband == 12][0] == pytest.approx(0.25)
    joint = PickingFreeSpectralLikelihood(f, SimpleNamespace(), obs, SMVIConfig())
    joint.set_frequency_weights(f.freqs_broadband <= 12)
    np.testing.assert_allclose(joint.freq_weights, lik.freq_weights)
    joint.set_frequency_weights((f.freqs_broadband <= 3) | (f.freqs_broadband >= 14))
    assert joint.freq_weights.sum() == pytest.approx(
        2.0 + 41.0
    )  # Two retained bands, not the excluded gap.


def prediction_target():
    class Target:
        param_space = SimpleNamespace(latent_to_physical=lambda z: z)

        def evaluate_physical_model(self, m):
            return (
                0.0,
                0.0,
                0.0,
                0.0,
                np.array([[m[0], 2 * m[0]], [1 + m[0], 2 - m[0]]]),
            )

        def spectral_metrics(self, s):
            return float((s.mean() - 1) ** 2), 0.0, float(abs(s[0, 0]))

    return Target()


def test_balanced_prediction_uses_mixture_weights_and_separates_mc_error():
    from posterior_prediction import evaluate_mixture_prediction

    means = np.array([[-1.0], [2.0]])
    L = np.array([[[0.3]], [[0.5]]])
    w = np.array([0.25, 0.75])
    cfg = SMVIConfig(
        predictive_draws=64,
        predictive_max_draws=128,
        predictive_bootstrap_repeats=16,
        predictive_mc_tolerance=0.2,
    )
    result = evaluate_mixture_prediction(prediction_target(), means, L, w, cfg)
    labels = result["predictive_comp_labels"]
    z = result["predictive_samples_z"][:, 0]
    count = result["predictive_component_counts"]
    assert count[0] == count[1]
    mu = sum(w[k] * z[labels == k].mean() for k in range(2))
    se = np.sqrt(
        sum(w[k] ** 2 * z[labels == k].var(ddof=1) / count[k] for k in range(2))
    )
    assert result["d_posterior_predictive_mean"][0, 0] == pytest.approx(mu)
    assert abs(mu - 1.25) < 0.2
    assert result["d_posterior_predictive_mc_se"][0, 0] == pytest.approx(se)
    assert result["d_posterior_predictive_std"][0, 0] > 10 * se
    assert result["predictive_metrics"]["precision_target_met"]


def test_predictive_precision_failure_is_reported_at_the_budget_limit():
    from posterior_prediction import evaluate_mixture_prediction

    cfg = SMVIConfig(
        predictive_draws=32,
        predictive_max_draws=64,
        predictive_bootstrap_repeats=8,
        predictive_mc_tolerance=1e-15,
    )
    result = evaluate_mixture_prediction(
        prediction_target(),
        np.array([[-1.0], [2.0]]),
        np.array([[[0.3]], [[0.5]]]),
        np.array([0.5, 0.5]),
        cfg,
    )
    m = result["predictive_metrics"]
    assert m["draws"] == 64 and not m["precision_target_met"]
    assert [x["draws"] for x in m["precision_history"]] == [32, 64]


def test_report_regeneration_preserves_structured_warnings_and_separate_notes(
    tmp_path, monkeypatch
):
    import generate_pdf_manuscript as generator
    from run_benchmark import code_hashes

    b = {
        "benchmark_name": "test",
        "mean_model_metrics": {"misfit": 0.4},
        "predictive_metrics": {
            "misfit": 0.5,
            "misfit_mc_se": 0.02,
            "draws": 128,
            "precision_target_met": False,
        },
        "baseline_comparison": {
            "prior_midpoint": {"misfit": 0.3},
            "truth_model_evaluation_only": {"misfit": 0.2},
            "warnings": ["predictive_point_estimate_worse_than_prior_midpoint"],
        },
    }
    hashes = code_hashes()
    report = {
        "status": "completed",
        "version": "test",
        "mode": "quick_verification",
        "source_sha256": hashes,
        "benchmark_a": b,
        "benchmark_b": b,
    }
    (tmp_path / "benchmark_verification_summary.json").write_text(json.dumps(report))
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"status": "completed", "source_sha256": hashes})
    )
    notes = tmp_path / "DIAGNOSTIC_EDITORIAL_NOTES.md"
    notes.write_text("catatan terpisah yang harus bertahan")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_pdf_manuscript.py",
            "--results-dir",
            str(tmp_path),
            "--diagnostic-report",
        ],
    )
    generator.main()
    first = (tmp_path / "DIAGNOSTIC_VERIFICATION_REPORT.md").read_text()
    generator.main()
    second = (tmp_path / "DIAGNOSTIC_VERIFICATION_REPORT.md").read_text()
    assert first == second
    assert "predictive_point_estimate_worse_than_prior_midpoint" in second
    assert "Model tengah prior:" in second and "0.02" in second
    assert notes.read_text() in second


def test_likelihood_rejects_complex_spectra_instead_of_discarding_imaginary_data():
    f = fake_fuser()
    obs = np.ones((len(f.freqs_broadband), len(f.c_grid)), dtype=complex)
    with pytest.raises(ValueError, match="real-valued"):
        PickingFreeSpectralLikelihood(f, SimpleNamespace(), obs, SMVIConfig())


def test_likelihood_rejects_dc_frequency_before_green_evaluation():
    f = fake_fuser()
    f.freqs_broadband[0] = 0.0
    obs = np.ones((len(f.freqs_broadband), len(f.c_grid)))
    with pytest.raises(ValueError, match="frequency grid"):
        PickingFreeSpectralLikelihood(f, SimpleNamespace(), obs, SMVIConfig())
