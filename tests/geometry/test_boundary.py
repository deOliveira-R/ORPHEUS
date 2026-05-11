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
    AlbedoBoundaryOperator,
    MixedBoundaryOperator,
    PeriodicBoundaryOperator,
    BoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    ProductQuadrature,
)


# ─────────────────────────────────────────────────────────────────────
# Wave 10 helper: migrate 2-arg ``bc.apply(psi, quad)`` legacy call
# sites to the realizer path. See ``transient-giggling-cake.md`` Wave
# 10 for the scope-narrowing decision (the C10.2 "drop 2-arg overload
# + delete _BoundBoundaryOperator" steps are deferred — see the
# follow-up issue filed at end of Wave 10).
# ─────────────────────────────────────────────────────────────────────


def _realize_for_sn(bc, quad):
    """Realize a BC via the SN realizer with a minimal method space.

    Used by tests that previously called ``bc.apply(psi, quad)`` directly.
    Returns a 1-arg :class:`LinearOperator` whose ``apply(psi)`` matches
    the legacy 2-arg ``bc.apply(psi, quad)`` semantics for all rank-1
    BCs **except** vacuum (which requires per-face inflow indices; see
    the vacuum-equivalent special cases below that retain the legacy
    2-arg form intentionally).

    Realizer-path output is bit-equivalent (or ``nulp <= 4`` at worst)
    to the legacy ``bc.apply(psi, quad)`` for non-vacuum BCs. Verified
    by Wave 5's ``test_sn_boundary_realizer.py`` equivalence tests.
    """
    realizer = SNBoundaryRealizer()
    method_space = SNMethodSpace.minimal(quad)
    return realizer.realize(bc, method_space)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — protocol semantics on synthetic data
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_vacuum_bc_returns_zero() -> None:
    """``VacuumBoundaryOperator.apply`` returns zeros exactly."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(0).standard_normal((quad.N, 3))
    bc = VacuumBoundaryOperator()

    # Wave 10 deviation: legacy 2-arg form retained — vacuum's
    # standalone apply tests the zeros-all contract per Wave 7
    # (Option a); realizer path is tested separately in
    # test_bc_equivalence_snapshot.py::TestVacuumLebedev17Snapshot::
    # test_realizer_zeroes_only_inflow_per_section_16A5.
    psi_in = bc.apply(psi_out, quad)

    assert psi_in.shape == psi_out.shape
    np.testing.assert_array_equal(psi_in, 0.0)


@pytest.mark.foundation
def test_vacuum_bc_is_resolved_bc() -> None:
    """The Protocol is runtime-checkable: VacuumBoundaryOperator qualifies."""
    assert isinstance(VacuumBoundaryOperator(), BoundaryOperator)


@pytest.mark.foundation
def test_specular_bc_indexes_through_reflection_partner() -> None:
    """``SpecularBoundaryOperator.apply(psi)[n] == psi[reflection_index[n]]``."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.arange(quad.N * 2, dtype=float).reshape(quad.N, 2)
    bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
    ref = quad.reflection_index("x")

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # psi_in[n] should equal psi_out[ref[n]]
    np.testing.assert_array_equal(psi_in, psi_out[ref])


@pytest.mark.foundation
def test_specular_bc_with_partial_albedo() -> None:
    """SpecularBoundaryOperator scales by ``albedo``."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.array([[1.0], [2.0], [3.0], [4.0]])
    bc = SpecularBoundaryOperator(axis="x", albedo=0.5)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # ref[n] = N - 1 - n: psi_out reversed.
    expected = 0.5 * psi_out[::-1]
    np.testing.assert_array_equal(psi_in, expected)


@pytest.mark.foundation
def test_specular_bc_axis_y_on_lebedev() -> None:
    """Lebedev y-reflection partner: ``apply(axis='y')`` matches index."""
    quad = LebedevSphere.create(order=9)
    psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
    bc = SpecularBoundaryOperator(axis="y", albedo=1.0)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, psi_out[quad.reflection_index("y")])


@pytest.mark.foundation
def test_white_bc_returns_cosine_weighted_average() -> None:
    """White BC: incoming = ``Σ w·μ·ψ_out / Σ w·μ`` over outgoing hemisphere."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    rng = np.random.default_rng(2)
    psi_out = rng.standard_normal((quad.N, 1))
    bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

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
    bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=0.7)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # All rows identical.
    for n in range(1, quad.N):
        np.testing.assert_allclose(psi_in[n], psi_in[0], rtol=1e-15, atol=0.0)


@pytest.mark.foundation
def test_white_bc_albedo_scales_linearly() -> None:
    """Two white BCs at α=1 and α=0.4 differ by exactly the ratio."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(4).standard_normal((quad.N, 2))
    bc1 = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
    bc2 = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=0.4)

    psi1 = _realize_for_sn(bc1, quad).apply(psi_out)
    psi2 = _realize_for_sn(bc2, quad).apply(psi_out)

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

    bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # Average should be (Σ_outgoing w·μ·1) / (Σ_outgoing w·μ) = 1.
    np.testing.assert_allclose(psi_in, 1.0, rtol=1e-14, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_axis_z_on_product_quadrature() -> None:
    """White BC on z-axis routes through ``mu_z``."""
    quad = ProductQuadrature.create(n_mu=8, n_phi=4)
    psi_out = np.where(quad.mu_z[:, None] > 0, 2.0, 9.0)
    bc = WhiteBoundaryOperator(axis="z", outward_sign=+1, albedo=1.0)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_allclose(psi_in, 2.0, rtol=1e-13, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_z_axis_unsupported_on_1d_quadrature() -> None:
    """1-D GL has no ``mu_z`` — z-axis raises ValueError."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.zeros((quad.N, 1))
    bc = WhiteBoundaryOperator(axis="z", outward_sign=+1, albedo=1.0)

    # Realizer constructs AngularAverageOperator.from_quadrature(quad,
    # axis="z", ...) eagerly; the ValueError fires inside the realize
    # call (before .apply runs). The ``pytest.raises`` context wraps
    # both, so the test discriminates on the ValueError regardless of
    # which call site (realize vs. apply) raises.
    with pytest.raises(ValueError, match="mu_z"):
        _realize_for_sn(bc, quad).apply(psi_out)


@pytest.mark.foundation
def test_periodic_bc_returns_input_unchanged() -> None:
    """PeriodicBoundaryOperator is identity on the angular axis (smoke test).

    The spatial pushforward is the caller's responsibility; this
    primitive only certifies that the angular structure is identity,
    which is a no-op return of the (caller-supplied partner-face
    outgoing) buffer.
    """
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(5).standard_normal((quad.N, 3))
    bc = PeriodicBoundaryOperator()

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, psi_out)
    # Returned array is a copy (independent of caller's buffer).
    psi_in[0, 0] = 1e9
    assert psi_out[0, 0] != 1e9


@pytest.mark.foundation
def test_albedo_bc_scales_outgoing() -> None:
    """``AlbedoBoundaryOperator(α).apply(ψ_out) == α·ψ_out``."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(6).standard_normal((quad.N, 2))
    bc = AlbedoBoundaryOperator(albedo=0.5)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, 0.5 * psi_out)


@pytest.mark.foundation
def test_albedo_bc_zero_is_vacuum_equivalent() -> None:
    """``AlbedoBoundaryOperator(0)`` is equivalent to ``VacuumBoundaryOperator`` value-wise."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(7).standard_normal((quad.N, 2))

    # Wave 10 deviation: legacy 2-arg form retained — vacuum's
    # standalone apply tests the zeros-all contract per Wave 7
    # (Option a). The realizer path returns IncomingOrdinateMaskTensor
    # which zeroes ONLY inflow rows (§16A.5); legacy returns all-zero.
    # Equivalence-on-inflow is tested in
    # test_bc_equivalence_snapshot.py.
    np.testing.assert_array_equal(
        AlbedoBoundaryOperator(albedo=0.0).apply(psi_out, quad),
        VacuumBoundaryOperator().apply(psi_out, quad),
    )


@pytest.mark.foundation
def test_mixed_bc_is_linear_combination() -> None:
    """``MixedBoundaryOperator([(c1, A), (c2, B)]).apply == c1·A.apply + c2·B.apply``."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    psi_out = np.random.default_rng(8).standard_normal((quad.N, 2))
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)

    bc = MixedBoundaryOperator([(0.3, spec), (0.7, white)])

    psi_mixed = _realize_for_sn(bc, quad).apply(psi_out)
    psi_expected = (
        0.3 * _realize_for_sn(spec, quad).apply(psi_out)
        + 0.7 * _realize_for_sn(white, quad).apply(psi_out)
    )

    np.testing.assert_allclose(psi_mixed, psi_expected, rtol=1e-14, atol=0.0)


@pytest.mark.foundation
def test_mixed_bc_empty_is_vacuum() -> None:
    """``MixedBoundaryOperator([])`` returns zero (the empty sum)."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(9).standard_normal((quad.N, 1))

    bc = MixedBoundaryOperator([])
    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, 0.0)


@pytest.mark.foundation
def test_mixed_bc_singleton_matches_primitive() -> None:
    """``MixedBoundaryOperator([(1.0, A)])`` is equivalent to A directly."""
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(10).standard_normal((quad.N, 2))
    primitive = SpecularBoundaryOperator(axis="x", albedo=0.8)
    bc = MixedBoundaryOperator([(1.0, primitive)])

    np.testing.assert_array_equal(
        _realize_for_sn(bc, quad).apply(psi_out),
        _realize_for_sn(primitive, quad).apply(psi_out),
    )


@pytest.mark.foundation
def test_mixed_bc_coefficients_chain() -> None:
    """``MixedBoundaryOperator`` coefficients distribute through linear primitives.

    ``MixedBoundaryOperator([(c, MixedBoundaryOperator([(d, X)]))])`` = ``MixedBoundaryOperator([(c·d, X)])`` via
    the apply method's linearity.
    """
    quad = GaussLegendre1D.create(n_ordinates=4)
    psi_out = np.random.default_rng(11).standard_normal((quad.N, 1))
    primitive = SpecularBoundaryOperator(axis="x", albedo=1.0)

    inner = MixedBoundaryOperator([(0.5, primitive)])
    outer = MixedBoundaryOperator([(0.4, inner)])
    flat = MixedBoundaryOperator([(0.5 * 0.4, primitive)])

    np.testing.assert_allclose(
        _realize_for_sn(outer, quad).apply(psi_out),
        _realize_for_sn(flat, quad).apply(psi_out),
        rtol=1e-15,
        atol=0.0,
    )


@pytest.mark.foundation
def test_all_primitives_are_resolved_bc() -> None:
    """The runtime-checkable Protocol accepts every primitive."""
    instances: list[BoundaryOperator] = [
        VacuumBoundaryOperator(),
        SpecularBoundaryOperator(axis="x", albedo=1.0),
        WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0),
        PeriodicBoundaryOperator(),
        AlbedoBoundaryOperator(albedo=0.5),
        MixedBoundaryOperator([]),
        MixedBoundaryOperator([(0.3, SpecularBoundaryOperator(axis="x")), (0.7, WhiteBoundaryOperator(axis="x"))]),
    ]
    for bc in instances:
        assert isinstance(bc, BoundaryOperator), f"{type(bc).__name__} not BoundaryOperator"


# ═══════════════════════════════════════════════════════════════════════
# Issue 9.6 — registry + composition + transpose
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_registry_contains_all_primitives() -> None:
    """All six concrete BC subtypes self-register under their key.

    Wave 7 merged the legacy ``BoundaryOperator`` ABC into
    :class:`BoundaryTraceLaw`, so the registry now also holds test
    stubs from :mod:`tests.geometry.test_boundary_trace_law` and the
    Wave-7 additions (``prescribed_inflow``). Assert ``>=`` rather
    than equality.
    """
    expected_keys = {
        "vacuum", "reflective", "white", "periodic", "albedo", "mixed",
    }
    assert expected_keys <= set(BoundaryOperator.registry.keys())


@pytest.mark.foundation
def test_registry_create_returns_concrete_instance() -> None:
    bc = BoundaryOperator.create("vacuum")
    assert isinstance(bc, VacuumBoundaryOperator)
    bc = BoundaryOperator.create("reflective", axis="x", albedo=1.0)
    assert isinstance(bc, SpecularBoundaryOperator)
    assert bc.axis == "x"
    assert bc.albedo == 1.0


@pytest.mark.foundation
def test_registry_create_unknown_key_raises() -> None:
    with pytest.raises(KeyError):
        BoundaryOperator.create("nonsense")


@pytest.mark.foundation
def test_specular_capabilities_include_transpose() -> None:
    """SpecularBoundaryOperator advertises apply_transpose."""
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    assert "apply_transpose" in spec.capabilities


@pytest.mark.foundation
def test_specular_apply_transpose_reciprocity_unweighted() -> None:
    r"""``<B(psi_out), phi_in> == <psi_out, B^T(phi_in)>``.

    Tests the Euclidean inner-product reciprocity identity for clean
    axis specular reflection. Index permutations under axis
    reflection are involutions, so the transpose acts as the same
    permutation.
    """
    quad = GaussLegendre1D.create(n_ordinates=8)
    rng = np.random.default_rng(42)
    psi_out = rng.standard_normal((quad.N, 2))
    phi_in = rng.standard_normal((quad.N, 2))

    spec = SpecularBoundaryOperator(axis="x", albedo=0.7)
    realized = _realize_for_sn(spec, quad)
    Bpsi = realized.apply(psi_out)
    BTphi = realized.apply_transpose(phi_in)

    lhs = float(np.sum(Bpsi * phi_in))
    rhs = float(np.sum(psi_out * BTphi))
    assert np.isclose(lhs, rhs, rtol=1e-13)


@pytest.mark.foundation
def test_specular_self_inverse_identity() -> None:
    r"""``B(B(x)) == albedo^2 * x`` for a clean axis reflection."""
    quad = GaussLegendre1D.create(n_ordinates=8)
    rng = np.random.default_rng(7)
    x = rng.standard_normal((quad.N, 2))

    spec = SpecularBoundaryOperator(axis="x", albedo=0.7)
    realized = _realize_for_sn(spec, quad)
    once = realized.apply(x)
    twice = realized.apply(once)

    expected = (0.7 ** 2) * x
    np.testing.assert_allclose(twice, expected, rtol=1e-13, atol=1e-14)


@pytest.mark.foundation
def test_operator_sum_of_bcs_matches_mixed_baseline() -> None:
    r"""``0.7 * Specular + 0.3 * White`` matches Mixed BC baseline.

    With Issue 9.6, BCs participate in operator algebra: scalar
    multiplication and operator addition produce composable
    objects whose ``apply(psi_out, quad)`` reproduces the explicit
    weighted sum that ``MixedBoundaryOperator`` already provided.
    """
    quad = GaussLegendre1D.create(n_ordinates=8)
    rng = np.random.default_rng(99)
    psi_out = rng.standard_normal((quad.N, 2))

    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)

    # Build the composition on the REALIZED 1-arg LinearOperators so
    # the Wave-0 OperatorSum algebra sees them as plain
    # ``LinearOperator``s rather than legacy 2-arg BCs. The realized
    # ``MixedBoundaryOperator`` path is recursive — the realizer
    # decomposes the components and reassembles via ``OperatorSum``
    # internally — so the two sides are structurally identical and
    # remain bit-comparable.
    composed = 0.7 * _realize_for_sn(spec, quad) + 0.3 * _realize_for_sn(white, quad)
    baseline = MixedBoundaryOperator([(0.7, spec), (0.3, white)])

    np.testing.assert_allclose(
        composed.apply(psi_out),
        _realize_for_sn(baseline, quad).apply(psi_out),
        rtol=1e-13,
        atol=1e-14,
    )


@pytest.mark.foundation
def test_boundary_operator_keys_match_class_attribute() -> None:
    """The ``key`` ClassVar must match the registry insertion key."""
    assert VacuumBoundaryOperator.key == "vacuum"
    assert SpecularBoundaryOperator.key == "reflective"
    assert WhiteBoundaryOperator.key == "white"
    assert PeriodicBoundaryOperator.key == "periodic"
    assert AlbedoBoundaryOperator.key == "albedo"
    assert MixedBoundaryOperator.key == "mixed"
