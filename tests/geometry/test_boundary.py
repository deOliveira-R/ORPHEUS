r"""Unit tests for tensor-decomposed boundary conditions.

The BCs in :mod:`orpheus.geometry.boundary` are rank-1 (vacuum,
specular, white, periodic, albedo) or rank-N (mixed) realisations of
the tensor decomposition

.. math::

    R = \sum_\alpha G_\alpha \otimes A_\alpha

acting on the outgoing angular flux at a boundary face. These tests
pin the algebraic semantics of each primitive: vacuum is zero,
specular is index-permutation by a chosen axis, white is the
cosine-weighted average of the outgoing hemisphere broadcast over the
incoming hemisphere, periodic is identity (the spatial pushforward
lives outside the protocol), albedo is scalar multiplication, and
mixed is linear in the components.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBC,
    MixedBC,
    PeriodicBC,
    ResolvedBC,
    SpecularBC,
    VacuumBC,
    WhiteBC,
)
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    ProductQuadrature,
)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — protocol semantics on synthetic data
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_vacuum_bc_returns_zero() -> None:
    """``VacuumBC.apply_to_incoming`` returns zeros exactly."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(0).standard_normal((quad.N, 3))
    bc = VacuumBC()

    psi_in = bc.apply_to_incoming(psi_out, quad)

    assert psi_in.shape == psi_out.shape
    np.testing.assert_array_equal(psi_in, 0.0)


@pytest.mark.foundation
def test_vacuum_bc_is_resolved_bc() -> None:
    """The Protocol is runtime-checkable: VacuumBC qualifies."""
    assert isinstance(VacuumBC(), ResolvedBC)


@pytest.mark.foundation
def test_specular_bc_indexes_through_reflection_partner() -> None:
    """``SpecularBC.apply_to_incoming(psi)[n] == psi[reflection_index[n]]``."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.arange(quad.N * 2, dtype=float).reshape(quad.N, 2)
    bc = SpecularBC(axis="x", albedo=1.0)
    ref = quad.reflection_index("x")

    psi_in = bc.apply_to_incoming(psi_out, quad)

    # psi_in[n] should equal psi_out[ref[n]]
    np.testing.assert_array_equal(psi_in, psi_out[ref])


@pytest.mark.foundation
def test_specular_bc_with_partial_albedo() -> None:
    """SpecularBC scales by ``albedo``."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.array([[1.0], [2.0], [3.0], [4.0]])
    bc = SpecularBC(axis="x", albedo=0.5)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    # ref[n] = N - 1 - n: psi_out reversed.
    expected = 0.5 * psi_out[::-1]
    np.testing.assert_array_equal(psi_in, expected)


@pytest.mark.foundation
def test_specular_bc_axis_y_on_lebedev() -> None:
    """Lebedev y-reflection partner: ``apply(axis='y')`` matches index."""
    quad = LebedevSphere.create(order=9)
    psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
    bc = SpecularBC(axis="y", albedo=1.0)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    np.testing.assert_array_equal(psi_in, psi_out[quad.reflection_index("y")])


@pytest.mark.foundation
def test_white_bc_returns_cosine_weighted_average() -> None:
    """White BC: incoming = ``Σ w·μ·ψ_out / Σ w·μ`` over outgoing hemisphere."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    rng = np.random.default_rng(2)
    psi_out = rng.standard_normal((quad.N, 1))
    bc = WhiteBC(axis="x", outward_sign=+1, albedo=1.0)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    mu = quad.mu_x
    w = quad.weights
    # Hand-compute the expected cosine-weighted average over outgoing
    # hemisphere (μ > 0 for outward_sign = +1).
    out_mask = mu > 0.0
    cos_w = np.where(out_mask, w * mu, 0.0)
    expected_avg = (cos_w[:, None] * psi_out).sum(axis=0) / cos_w.sum()
    expected = np.broadcast_to(expected_avg, psi_out.shape)

    np.testing.assert_allclose(psi_in, expected, rtol=1e-15, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_is_angle_independent() -> None:
    """White BC produces the same value at every ordinate index."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(3).standard_normal((quad.N, 2))
    bc = WhiteBC(axis="x", outward_sign=+1, albedo=0.7)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    # All rows identical.
    for n in range(1, quad.N):
        np.testing.assert_allclose(psi_in[n], psi_in[0], rtol=1e-15, atol=0.0)


@pytest.mark.foundation
def test_white_bc_albedo_scales_linearly() -> None:
    """Two white BCs at α=1 and α=0.4 differ by exactly the ratio."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(4).standard_normal((quad.N, 2))
    bc1 = WhiteBC(axis="x", outward_sign=+1, albedo=1.0)
    bc2 = WhiteBC(axis="x", outward_sign=+1, albedo=0.4)

    psi1 = bc1.apply_to_incoming(psi_out, quad)
    psi2 = bc2.apply_to_incoming(psi_out, quad)

    np.testing.assert_allclose(psi2, 0.4 * psi1, rtol=1e-14, atol=0.0)


@pytest.mark.foundation
def test_white_bc_4_point_quadrature_hand_computed() -> None:
    """White BC on 4-point GL: explicit hand calculation."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    # 4-point GL: μ = ±0.339981..., ±0.861136..., w = 0.652145..., 0.347855...
    mu_pos = quad.mu_x[quad.mu_x > 0]
    w_pos = quad.weights[quad.mu_x > 0]
    # Set ψ_out: 1.0 on ordinates with μ > 0 (outgoing for outward_sign=+1),
    # 7.0 on incoming (irrelevant for the average).
    psi_out = np.where(quad.mu_x[:, None] > 0, 1.0, 7.0)

    bc = WhiteBC(axis="x", outward_sign=+1, albedo=1.0)
    psi_in = bc.apply_to_incoming(psi_out, quad)

    # Average should be (Σ_outgoing w·μ·1) / (Σ_outgoing w·μ) = 1.
    np.testing.assert_allclose(psi_in, 1.0, rtol=1e-14, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_axis_z_on_product_quadrature() -> None:
    """White BC on z-axis routes through ``mu_z``."""
    quad = ProductQuadrature.create(n_mu=8, n_phi=4)
    psi_out = np.where(quad.mu_z[:, None] > 0, 2.0, 9.0)
    bc = WhiteBC(axis="z", outward_sign=+1, albedo=1.0)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    np.testing.assert_allclose(psi_in, 2.0, rtol=1e-13, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_z_axis_unsupported_on_1d_quadrature() -> None:
    """1-D GL has no ``mu_z`` — z-axis raises ValueError."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.zeros((quad.N, 1))
    bc = WhiteBC(axis="z", outward_sign=+1, albedo=1.0)

    with pytest.raises(ValueError, match="mu_z"):
        bc.apply_to_incoming(psi_out, quad)


@pytest.mark.foundation
def test_periodic_bc_returns_input_unchanged() -> None:
    """PeriodicBC is identity on the angular axis (smoke test).

    The spatial pushforward is the caller's responsibility; this
    primitive only certifies that the angular structure is identity,
    which is a no-op return of the (caller-supplied partner-face
    outgoing) buffer.
    """
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(5).standard_normal((quad.N, 3))
    bc = PeriodicBC()

    psi_in = bc.apply_to_incoming(psi_out, quad)

    np.testing.assert_array_equal(psi_in, psi_out)
    # Returned array is a copy (independent of caller's buffer).
    psi_in[0, 0] = 1e9
    assert psi_out[0, 0] != 1e9


@pytest.mark.foundation
def test_albedo_bc_scales_outgoing() -> None:
    """``AlbedoBC(α).apply(ψ_out) == α·ψ_out``."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(6).standard_normal((quad.N, 2))
    bc = AlbedoBC(albedo=0.5)

    psi_in = bc.apply_to_incoming(psi_out, quad)

    np.testing.assert_array_equal(psi_in, 0.5 * psi_out)


@pytest.mark.foundation
def test_albedo_bc_zero_is_vacuum_equivalent() -> None:
    """``AlbedoBC(0)`` is equivalent to ``VacuumBC`` value-wise."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(7).standard_normal((quad.N, 2))

    np.testing.assert_array_equal(
        AlbedoBC(albedo=0.0).apply_to_incoming(psi_out, quad),
        VacuumBC().apply_to_incoming(psi_out, quad),
    )


@pytest.mark.foundation
def test_mixed_bc_is_linear_combination() -> None:
    """``MixedBC([(c1, A), (c2, B)]).apply == c1·A.apply + c2·B.apply``."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(8).standard_normal((quad.N, 2))
    spec = SpecularBC(axis="x", albedo=1.0)
    white = WhiteBC(axis="x", outward_sign=+1, albedo=1.0)

    bc = MixedBC([(0.3, spec), (0.7, white)])

    psi_mixed = bc.apply_to_incoming(psi_out, quad)
    psi_expected = (
        0.3 * spec.apply_to_incoming(psi_out, quad)
        + 0.7 * white.apply_to_incoming(psi_out, quad)
    )

    np.testing.assert_allclose(psi_mixed, psi_expected, rtol=1e-14, atol=0.0)


@pytest.mark.foundation
def test_mixed_bc_empty_is_vacuum() -> None:
    """``MixedBC([])`` returns zero (the empty sum)."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(9).standard_normal((quad.N, 1))

    bc = MixedBC([])
    psi_in = bc.apply_to_incoming(psi_out, quad)

    np.testing.assert_array_equal(psi_in, 0.0)


@pytest.mark.foundation
def test_mixed_bc_singleton_matches_primitive() -> None:
    """``MixedBC([(1.0, A)])`` is equivalent to A directly."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(10).standard_normal((quad.N, 2))
    primitive = SpecularBC(axis="x", albedo=0.8)
    bc = MixedBC([(1.0, primitive)])

    np.testing.assert_array_equal(
        bc.apply_to_incoming(psi_out, quad),
        primitive.apply_to_incoming(psi_out, quad),
    )


@pytest.mark.foundation
def test_mixed_bc_coefficients_chain() -> None:
    """``MixedBC`` coefficients distribute through linear primitives.

    ``MixedBC([(c, MixedBC([(d, X)]))])`` = ``MixedBC([(c·d, X)])`` via
    the apply method's linearity.
    """
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(11).standard_normal((quad.N, 1))
    primitive = SpecularBC(axis="x", albedo=1.0)

    inner = MixedBC([(0.5, primitive)])
    outer = MixedBC([(0.4, inner)])
    flat = MixedBC([(0.5 * 0.4, primitive)])

    np.testing.assert_allclose(
        outer.apply_to_incoming(psi_out, quad),
        flat.apply_to_incoming(psi_out, quad),
        rtol=1e-15,
        atol=0.0,
    )


@pytest.mark.foundation
def test_all_primitives_are_resolved_bc() -> None:
    """The runtime-checkable Protocol accepts every primitive."""
    instances: list[ResolvedBC] = [
        VacuumBC(),
        SpecularBC(axis="x", albedo=1.0),
        WhiteBC(axis="x", outward_sign=+1, albedo=1.0),
        PeriodicBC(),
        AlbedoBC(albedo=0.5),
        MixedBC([]),
        MixedBC([(0.3, SpecularBC(axis="x")), (0.7, WhiteBC(axis="x"))]),
    ]
    for bc in instances:
        assert isinstance(bc, ResolvedBC), f"{type(bc).__name__} not ResolvedBC"
