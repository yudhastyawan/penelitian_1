"""Checks against independent physics, known targets and numerical derivatives."""

from types import SimpleNamespace
import numpy as np
import pytest
from scipy.special import logsumexp
from models import get_benchmark_a_lvl, get_benchmark_b_water_table
from forward_mfj import (
    ActiveSurveyGeometry,
    PassiveSurveyGeometry,
    FullWavefieldMFJOperator,
)
from fusion import ActivePassiveSpectralFuser
from smvi_inversion import (
    SMVIConfig,
    PickingFreeSpectralLikelihood,
    SteinMixtureVariationalInference,
)


@pytest.mark.parametrize("factory", [get_benchmark_a_lvl, get_benchmark_b_water_table])
def test_layered_compliance_against_independent_pyprop8(factory):
    pytest.importorskip("pyprop8")
    from validate_reassessment import independent_compliance_check

    result = independent_compliance_check(factory)
    assert result["status"] == "passed", result


def test_fusion_uses_native_slices_without_inventing_frequencies():
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(),
        np.array([8.0, 12.0, 20.0, 30.0]),
        np.array([1.0, 3.0, 8.0, 10.0, 12.0]),
        np.arange(120.0, 1000.0, 20.0),
    )
    f = ActivePassiveSpectralFuser(op)
    np.testing.assert_array_equal(f.freqs_broadband, [1.0, 3.0, 8.0, 12.0, 20.0, 30.0])
    a = np.arange(4 * len(op.c_grid)).reshape(4, -1)
    copied = f._interp_spectrum_to_broadband(a, op.freqs_active)
    assert (
        10.0 not in f.freqs_broadband
    )  # A passive-only interior slice must not break the transition.
    np.testing.assert_array_equal(copied[4], a[2])
    fused = f.fuse_spectra(a, np.flip(a[:1].repeat(5, axis=0), axis=1))
    assert fused["w_active"][2] == 0.0
    assert fused["w_active"][3] == pytest.approx(1.0)
    assert fused["w_active"][-1] == pytest.approx(1.0)


@pytest.mark.parametrize(
    "factory,channel,velocity_max",
    [
        (get_benchmark_a_lvl, "fused", None),
        (get_benchmark_b_water_table, "fused", None),
        (get_benchmark_a_lvl, "active_bb", None),
        (get_benchmark_b_water_table, "passive_bb", 600.0),
    ],
)
def test_autodiff_matches_numba_spectra_target_and_finite_difference(
    factory, channel, velocity_max
):
    torch = pytest.importorskip("torch")
    from autodiff_likelihood import AutodiffSpectralLikelihood

    model, space = factory()
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(),
        np.array([8.0, 10.0, 16.0, 32.0]),
        np.array([1.0, 4.0, 8.0, 10.0, 12.0]),
        np.arange(120.0, 1801.0, 8.0),
        2.0,
    )
    f = ActivePassiveSpectralFuser(op)
    obs = f.generate_fused_spectrum_for_model(model, True, np.random.default_rng(5))[
        "d_" + channel
    ]
    cfg = SMVIConfig()
    base = PickingFreeSpectralLikelihood(f, space, obs, cfg, velocity_max, channel)
    ad = AutodiffSpectralLikelihood(f, space, obs, cfg, velocity_max, channel)
    z = np.random.default_rng(17).normal(0, 0.2, space.dim)
    image, lp, *_ = ad.spectrum_and_target(torch.tensor(z))
    syn = f.generate_fused_spectrum_for_model(
        space.vector_to_model(space.latent_to_physical(z))
    )["d_" + channel]
    np.testing.assert_allclose(image.detach().numpy(), syn, atol=2e-8, rtol=2e-8)
    assert float(lp.detach()) == pytest.approx(
        base.log_posterior_latent(z)[0], abs=2e-7
    )
    grad = ad.score_and_post_latent(z)[0]
    fd = base.score_and_post_latent(z, fd_eps=1e-5)[0]
    assert np.linalg.norm(grad - fd) / max(np.linalg.norm(fd), 1e-10) < 3e-4


def test_wiggles_include_aperture_ends_and_do_not_overlap():
    from visualize import select_radial_traces

    r = PassiveSurveyGeometry().interstation_offsets
    indices, width = select_radial_traces(r, 24)
    assert indices[0] == 0 and indices[-1] == len(r) - 1
    assert np.all(np.diff(r[indices]) > 2 * width)


def test_stein_kernel_acts_on_covariances_and_has_correct_normalization():
    engine = SteinMixtureVariationalInference(
        SimpleNamespace(
            param_space=SimpleNamespace(dim=1),
            config=SMVIConfig(n_components=2, repulsion_strength=0.3),
        )
    )
    theta = np.array(
        [[0.0, -0.7], [0.0, 0.7]]
    )  # Same means, different covariance parameters.
    force = np.array([[1.0, 2.0], [3.0, 4.0]])
    diff = theta[:, None] - theta[None, :]
    ds = (diff**2).sum(axis=-1)
    h = ds[0, 1] / np.log(3)
    k = np.exp(-ds / h)
    expected = (k.T @ force - 0.3 * 2 * np.sum(k[:, :, None] * diff, axis=0) / h) / 2
    np.testing.assert_allclose(engine._stein_kernel_update(theta, force), expected)
    assert abs(engine._stein_kernel_update(theta, np.zeros_like(theta))[0, 1]) > 0


def test_multiple_components_recover_both_known_modes():
    class Target:
        param_space = SimpleNamespace(dim=1, latent_to_physical=lambda z: z)
        config = SMVIConfig(
            n_components=4,
            n_iterations=350,
            n_mc_samples_per_comp=24,
            lr_mean=0.025,
            lr_chol=0.012,
            repulsion_strength=0.05,
            n_posterior_draws=10000,
            random_seed=72,
        )
        eval_count = 0

        def log_posterior_latent(self, z):
            return (
                logsumexp(-0.5 * ((z[0] - np.array([-2.5, 2.5])) / 0.4) ** 2),
                0.0,
                0.0,
                0.0,
            )

        def score_and_post_latent(self, z):
            self.eval_count += 1
            logs = -0.5 * ((z[0] - np.array([-2.5, 2.5])) / 0.4) ** 2
            grad = np.exp(logs - logsumexp(logs)) @ (
                -(z[0] - np.array([-2.5, 2.5])) / 0.16
            )
            return np.array([grad]), *self.log_posterior_latent(z)

        def evaluate_physical_model(self, m):
            return 0.0, 0.0, 0.0, 0.0, np.ones((2, 2))

    result = SteinMixtureVariationalInference(Target()).run(False)
    samples = result["posterior_samples_z"][:, 0]
    np.testing.assert_array_equal(result["weights"], np.full(4, 0.25))
    assert abs(samples.mean()) < 0.15
    assert abs(samples.std() - np.sqrt(0.16 + 2.5**2)) < 0.15
    assert 0.45 < np.mean(samples < 0) < 0.55
    assert np.all(abs(result["means_z"].ravel()) > 2.0)


def test_completed_quick_run_cannot_be_published():
    from research_quality import require_research_readiness

    with pytest.raises(ValueError, match="diagnostic_or_quick_run"):
        require_research_readiness(
            {"status": "completed", "mode": "quick_verification"}
        )


def test_independent_waveform_adapter_runs_without_hidden_failure():
    pytest.importorskip("pyprop8")
    from forward_mfj import PyProp8ForwardAdapter

    model, _ = get_benchmark_a_lvl()
    gather = PyProp8ForwardAdapter.compute_seismograms(
        model, np.array([20.0, 40.0]), 64, 0.02, 5.0
    )
    assert gather.shape == (64, 2)
    assert np.all(np.isfinite(gather))


def test_fft_convention_gives_expected_homogeneous_rayleigh_arrivals():
    from scipy.signal import hilbert
    from models import LayeredEarthModel

    model = LayeredEarthModel(
        "halfspace",
        "halfspace",
        np.array([0.0]),
        np.array([300.0]),
        np.array([600.0]),
        np.array([2000.0]),
        np.array([1000.0]),
        np.array([2000.0]),
    )
    op = FullWavefieldMFJOperator(
        ActiveSurveyGeometry(),
        PassiveSurveyGeometry(),
        np.array([10.0, 20.0, 40.0]),
        np.array([2.0, 4.0]),
        np.arange(150.0, 800.0, 2.0),
        2.0,
    )
    t, r, u = op.synthesize_time_gather(model, "active", False)
    indices = [np.argmin(abs(r - x)) for x in [40.0, 60.0, 80.0, 100.0, 114.0]]
    roots = np.roots([1.0, -8.0, 20.0, -12.0])
    root = next(x.real for x in roots if abs(x.imag) < 1e-10 and 0 < x.real < 1)
    expected = 0.045 + r[indices] / (300 * np.sqrt(root))
    mask = (t > 0.08) & (t < 1.0)
    arrivals = t[mask][np.argmax(abs(hilbert(u[mask][:, indices], axis=0)), axis=0)]
    np.testing.assert_allclose(arrivals, expected, atol=0.04, rtol=0.0)


def test_benchmark_overlap_is_synthesized_at_matching_native_frequencies():
    from run_benchmark import benchmark_frequency_grids

    for quick in [False, True]:
        fa, fp = benchmark_frequency_grids(quick)
        expected = np.arange(8.0, 12.01, 0.5)
        np.testing.assert_array_equal(fa[(fa >= 8) & (fa <= 12)], expected)
        np.testing.assert_array_equal(fp[(fp >= 8) & (fp <= 12)], expected)
        op = FullWavefieldMFJOperator(
            ActiveSurveyGeometry(),
            PassiveSurveyGeometry(),
            fa,
            fp,
            np.arange(120.0, 1000.0, 20.0),
        )
        f = ActivePassiveSpectralFuser(op)
        # Identical shape proxies isolate the transition policy from noise.
        shape = np.linspace(0.0, 1.0, len(op.c_grid))
        fused = f.fuse_spectra(
            np.tile(shape, (len(fa), 1)), np.tile(shape, (len(fp), 1))
        )
        mask = (f.freqs_broadband >= 8) & (f.freqs_broadband <= 12)
        weight = fused["w_active"][mask]
        assert np.all(np.diff(weight) >= -1e-12)
        assert weight[0] == 0.0 and weight[-1] == pytest.approx(1.0)
