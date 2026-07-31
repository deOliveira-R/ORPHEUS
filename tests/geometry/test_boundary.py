r"""Unit tests for tensor-decomposed boundary conditions.

The BCs in :mod:`orpheus.geometry.boundary` are rank-1 (vacuum,
specular, white, periodic, albedo) realisations of the tensor
decomposition

.. math::

    R = \sum_\alpha G_\alpha \otimes A_\alpha

acting on the outgoing angular flux at a boundary face. These tests
pin the algebraic semantics of each primitive: vacuum is zero,
specular is index-permutation by a chosen axis, white is the
cosine-weighted average of the outgoing hemisphere broadcast over the
incoming hemisphere, periodic is identity (the spatial pushforward
lives outside the protocol), and albedo is scalar multiplication.
Rank-N compositions are tested through the Wave-0 operator algebra
acting on realised primitives (the Wave-11 replacement for the
removed ``MixedBoundaryOperator`` composer); see
``tests/geometry/test_law_composition.py`` for the tree-walker tests.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    PeriodicBoundary,
    BoundaryTraceLaw,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature


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
    BCs **except** vacuum (which requires per-face inflow indices —
    use :func:`_realize_vacuum_for_face_right` for that).

    Realizer-path output is bit-equivalent (or ``nulp <= 4`` at worst)
    to the legacy ``bc.apply(psi, quad)`` for non-vacuum BCs. Verified
    by Wave 5's ``test_sn_boundary_realizer.py`` equivalence tests.
    """
    realizer = SNBoundaryRealizer()
    method_space = SNMethodSpace.minimal(quad)
    return realizer.realize(bc, method_space)


def _realize_vacuum_for_face_right(bc, quad):
    """Realize a vacuum BC for the right face (outward normal +x).

    Issue #176 / C176.2: the bare :class:`VacuumInflow.apply`
    returns ``np.zeros_like(psi)`` (legacy zeros-all, per Wave 7
    Option a). The §16A.5-correct path returns
    :class:`IncomingOrdinateMaskTensor` that zeros ONLY the inflow
    rows. This helper supplies the per-face inflow indices a vacuum
    realizer needs.

    For a 1-D-style "right" face (outward normal :math:`+\\hat x`),
    inflow ordinates are those with :math:`\\mu_x < -\\epsilon`
    (:math:`\\Omega \\cdot \\hat n_{\\text{right}} = +\\mu_x < 0`).
    """
    inflow_indices = np.flatnonzero(quad.mu_x < -1e-12)
    method_space = SNMethodSpace(
        quadrature=quad, face="xmax", inflow_indices=inflow_indices,
    )
    return SNBoundaryRealizer().realize(bc, method_space)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — protocol semantics on synthetic data
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_vacuum_bc_realizer_zeros_only_inflow_per_section_16A5() -> None:
    """Vacuum BC via the realizer zeros ONLY the inflow rows.

    Issue #176 / C176.2: migrated from the legacy zeros-all
    ``bc.apply(psi, quad)`` contract to the §16A.5-correct
    inflow-only-mask semantics. The realizer's vacuum branch
    returns an :class:`IncomingOrdinateMaskTensor` that zeros
    inflow rows (where :math:`\\Omega \\cdot \\hat n < 0`) and
    passes outflow rows through unchanged. This is the
    production contract used by every SN sweep after Issue #188
    (curvilinear) was wired through the realizer.

    The bare :class:`VacuumInflow.apply(psi, quad)`
    direct call still returns the legacy zeros-all output as a
    backward-compat fallback (documented in the BC docstring).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(0).standard_normal((quad.N, 3))
    bc = VacuumInflow()

    psi_in = _realize_vacuum_for_face_right(bc, quad).apply(psi_out)

    assert psi_in.shape == psi_out.shape
    # Inflow rows (mu_x < 0 at the right face) zeroed.
    inflow = np.flatnonzero(quad.mu_x < -1e-12)
    np.testing.assert_array_equal(psi_in[inflow], 0.0)
    # Outflow rows passed through unchanged.
    outflow = np.flatnonzero(quad.mu_x > 1e-12)
    np.testing.assert_array_equal(psi_in[outflow], psi_out[outflow])


@pytest.mark.foundation
def test_vacuum_bc_is_resolved_bc() -> None:
    """The Protocol is runtime-checkable: VacuumInflow qualifies."""
    assert isinstance(VacuumInflow(), BoundaryTraceLaw)


@pytest.mark.foundation
def test_specular_bc_indexes_through_reflection_partner() -> None:
    """``ReflectiveBoundary.apply(psi)[n] == psi[reflection_index[n]]``."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.arange(quad.N * 2, dtype=float).reshape(quad.N, 2)
    bc = ReflectiveBoundary(axis="x", albedo=1.0)
    ref = quad.reflection_index("x")

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # psi_in[n] should equal psi_out[ref[n]]
    np.testing.assert_array_equal(psi_in, psi_out[ref])


@pytest.mark.foundation
def test_specular_bc_with_partial_albedo() -> None:
    """ReflectiveBoundary scales by ``albedo``."""
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_out = np.array([[1.0], [2.0], [3.0], [4.0]])
    bc = ReflectiveBoundary(axis="x", albedo=0.5)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # ref[n] = N - 1 - n: psi_out reversed.
    expected = 0.5 * psi_out[::-1]
    np.testing.assert_array_equal(psi_in, expected)


@pytest.mark.foundation
def test_specular_bc_axis_y_on_lebedev() -> None:
    """Lebedev y-reflection partner: ``apply(axis='y')`` matches index."""
    quad = Quadrature.lebedev(order=9)
    psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
    bc = ReflectiveBoundary(axis="y", albedo=1.0)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, psi_out[quad.reflection_index("y")])


@pytest.mark.foundation
def test_white_bc_returns_cosine_weighted_average() -> None:
    """White BC: incoming = ``Σ w·μ·ψ_out / Σ w·μ`` over outgoing hemisphere."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(2)
    psi_out = rng.standard_normal((quad.N, 1))
    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

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
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(3).standard_normal((quad.N, 2))
    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.7)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # All rows identical.
    for n in range(1, quad.N):
        np.testing.assert_allclose(psi_in[n], psi_in[0], rtol=1e-15, atol=0.0)


@pytest.mark.foundation
def test_white_bc_albedo_scales_linearly() -> None:
    """Two white BCs at α=1 and α=0.4 differ by exactly the ratio."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(4).standard_normal((quad.N, 2))
    bc1 = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    bc2 = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.4)

    psi1 = _realize_for_sn(bc1, quad).apply(psi_out)
    psi2 = _realize_for_sn(bc2, quad).apply(psi_out)

    np.testing.assert_allclose(psi2, 0.4 * psi1, rtol=1e-14, atol=0.0)


# ─────────────────────────────────────────────────────────────────────
# Textbook Gauss-Legendre tables (Abramowitz & Stegun, Table 25.4) —
# the STRUCTURALLY INDEPENDENT reference constants for the two
# hand-computed white-BC gates below. Typed in from the published
# table, NOT read off ``quad.weights`` / ``quad.mu_*``: the threat
# model for those gates is the operator's own ``w·|Ω·n̂|`` weight
# formula, so a reference re-derived from the code's arrays would be a
# re-transcription (vv-principles L11 — procedural, not structural,
# independence). Each gate asserts the quadrature MATCHES this table
# first, so a change in ``Quadrature`` reds loudly instead of silently
# moving the reference.
# ─────────────────────────────────────────────────────────────────────

# 4-point GL on [-1, 1]: (|node|, weight), inner pair then outer pair.
_GL4_TABLE = (
    (0.3399810435848563, 0.6521451548625461),
    (0.8611363115940526, 0.3478548451374538),
)
# 8-point GL on [-1, 1]: positive half only (the rule is symmetric).
_GL8_POSITIVE_TABLE = (
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)


@pytest.mark.foundation
def test_white_bc_4_point_quadrature_hand_computed() -> None:
    r"""White BC on 4-point GL: an explicit hand calculation, against a
    ψ_out that is NOT constant over the outgoing hemisphere.

    B0.3 REPAIR — this gate was **blind** (measured). It previously fed
    a hemisphere-CONSTANT ψ_out (1.0 outgoing / 7.0 incoming) and
    asserted the output was that same constant. A normalised weighted
    average of a constant IS that constant **for any weights**, so the
    measured functional's invariance group contained the entire
    ``w·|Ω·n̂|`` weight formula (``vv-principles`` Mode 12): the test
    PASSED with the ``|Ω·n̂|`` factor dropped, with the wrong
    normaliser, and with the outgoing-hemisphere mask removed. It also
    contained no hand-computed number at all — the expected ``1.0``
    fell out of normalisation alone, so the docstring over-claimed.

    The repair breaks the invariance by giving the two outgoing
    ordinates DIFFERENT values, which makes the answer a genuine
    cosine-weighted mean:

    .. math::

        \psi_- \;=\;
        \frac{w_1\mu_1\,\psi_1 \;+\; w_2\mu_2\,\psi_2}
             {w_1\mu_1 \;+\; w_2\mu_2}

    with :math:`(\mu_1, w_1) = (0.33998104…,\,0.65214515…)`,
    :math:`(\mu_2, w_2) = (0.86113631…,\,0.34785485…)` from the
    published GL-4 table, :math:`\psi_1 = 1`, :math:`\psi_2 = 4`
    :math:`\Rightarrow \psi_- = 2.723973656470134`.

    The incoming ordinates carry 7.0 so the outgoing-hemisphere mask is
    also constrained (they must NOT enter the average).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    (mu_1, w_1), (mu_2, w_2) = _GL4_TABLE

    # Precondition: the quadrature IS the published GL-4 rule. If this
    # ever moves, the gate fails loudly rather than silently rebasing
    # its reference on the new arrays.
    np.testing.assert_allclose(
        np.sort(np.abs(quad.mu_x)), [mu_1, mu_1, mu_2, mu_2],
        rtol=1e-15, atol=0.0,
        err_msg="GL-4 nodes no longer match the published table",
    )
    np.testing.assert_allclose(
        np.sort(quad.weights), [w_2, w_2, w_1, w_1],
        rtol=1e-15, atol=0.0,
        err_msg="GL-4 weights no longer match the published table",
    )

    psi_1, psi_2, psi_incoming = 1.0, 4.0, 7.0
    psi_out = np.full((quad.N, 1), psi_incoming)
    psi_out[np.isclose(quad.mu_x, +mu_1)] = psi_1
    psi_out[np.isclose(quad.mu_x, +mu_2)] = psi_2

    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    expected = (w_1 * mu_1 * psi_1 + w_2 * mu_2 * psi_2) / (
        w_1 * mu_1 + w_2 * mu_2
    )
    assert abs(expected - 2.723973656470134) < 1e-14  # the docstring's number
    np.testing.assert_allclose(
        psi_in, expected, rtol=1e-14, atol=0.0,
        err_msg=(
            "white BC is not the |Ω·n̂|-weighted outgoing-hemisphere "
            "mean on GL-4"
        ),
    )


@pytest.mark.foundation
def test_white_bc_axis_z_on_product_quadrature() -> None:
    r"""White BC on the z-axis routes through ``mu_z`` — and returns the
    genuine cosine-weighted mean over the outgoing hemisphere.

    B0.3 REPAIR — the same Mode-12 blindness as the GL-4 gate above:
    the old ψ_out was hemisphere-constant (2.0 outgoing / 9.0 incoming)
    and the assertion was ``psi_in == 2.0``, which holds for any weight
    formula. Measured: it PASSED with the ``|Ω·n̂|`` factor dropped.

    The repair sets :math:`\psi_+ = \mu_z` on the outgoing hemisphere,
    so the answer is a true weighted mean of a NON-constant field.
    Because the product rule is (GL-8 in :math:`\mu`) ⊗ (uniform in
    :math:`\varphi`) and :math:`\psi` is :math:`\varphi`-independent,
    the azimuthal weights cancel and the reference reduces to the
    published GL-8 half-table:

    .. math::

        \psi_- \;=\;
        \frac{\sum_{\mu_k>0} w_k\,\mu_k^2}{\sum_{\mu_k>0} w_k\,\mu_k}
        \;=\; \frac{1/3}{0.50576403\ldots} \;=\; 0.659068878836894 .

    Two references, deliberately of different pillars:

    * the DISCRETE hand calculation above (tight, ``rtol=1e-13``);
    * the CONTINUUM closed form, which is the analytic identity the
      audit found the suite lacked —
      :math:`\bigl(\int_{2\pi}|\Omega\cdot\hat n|\mu\,d\Omega\bigr) /
      \bigl(\int_{2\pi}|\Omega\cdot\hat n|\,d\Omega\bigr)
      = (2\pi/3)/\pi = 2/3`.
      The numerator is EXACT under GL-8 (:math:`\mu^2` is degree 2),
      the denominator is not (:math:`|\mu|` has a kink at 0), so the
      discrete answer sits 1.14 % below 2/3 — a genuine, explained
      quadrature gap, gated at ``rtol=2e-2``. That is loose by design
      and still catches every white mutation (the dropped-cosine
      variant reads 0.5058, i.e. 24 % low).
    """
    quad = Quadrature.product(n_mu=8, n_phi=4)
    mu_k = np.array([m for m, _ in _GL8_POSITIVE_TABLE])
    w_k = np.array([w for _, w in _GL8_POSITIVE_TABLE])

    # Precondition: the product rule's polar nodes ARE the published
    # GL-8 nodes, and each carries the 2π-scaled GL-8 weight.
    np.testing.assert_allclose(
        np.unique(np.round(quad.mu_z, 12)),
        np.concatenate([-mu_k[::-1], mu_k]),
        rtol=1e-11, atol=0.0,
        err_msg="product(8, 4) polar nodes are no longer GL-8",
    )

    psi_incoming = 9.0
    psi_out = np.where(
        quad.mu_z[:, None] > 0, quad.mu_z[:, None], psi_incoming,
    )
    bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    # (1) discrete hand calculation off the published GL-8 half-table.
    expected = float((w_k * mu_k**2).sum() / (w_k * mu_k).sum())
    assert abs(expected - 0.659068878836894) < 1e-14  # docstring's number
    np.testing.assert_allclose(
        psi_in, expected, rtol=1e-13, atol=0.0,
        err_msg=(
            "white BC is not the |Ω·n̂|-weighted outgoing-hemisphere "
            "mean on product(8, 4)"
        ),
    )

    # (2) closed-form continuum anchor: ∫|Ω·n̂|μ dΩ / ∫|Ω·n̂| dΩ = 2/3.
    np.testing.assert_allclose(
        psi_in, 2.0 / 3.0, rtol=2e-2, atol=0.0,
        err_msg=(
            "white BC departs from the closed-form hemispheric "
            "cosine-weighted mean 2/3 by more than the GL-8 half-range "
            "quadrature gap (1.14 %)"
        ),
    )


@pytest.mark.foundation
def test_white_bc_z_axis_unsupported_on_1d_quadrature() -> None:
    r"""Slab GL has zero z-cosines — z-axis white BC fails the
    downstream "no outgoing ordinates" guard.

    R-1 Phase A detour-C: the unified :class:`Quadrature` class
    always exposes ``mu_z`` as a :func:`@property` view (returning
    zeros for slab measures); the failure mode that the OLD
    ``GaussLegendre1D`` adapter triggered via :exc:`AttributeError`
    on missing ``mu_z`` now surfaces one layer downstream — the
    :class:`AngularAverageOperator.from_quadrature` realizer sees
    that NO ordinate has outward :math:`\mu_z > 0` and raises
    ValueError with a "no outgoing ordinates" message. Same
    defensive behavior (slab z-axis is structurally degenerate);
    different layer of catch.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_out = np.zeros((quad.N, 1))
    bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)

    with pytest.raises(ValueError, match="no outgoing ordinates"):
        _realize_for_sn(bc, quad).apply(psi_out)


@pytest.mark.foundation
def test_periodic_bc_returns_input_unchanged() -> None:
    """PeriodicBoundary is identity on the angular axis (smoke test).

    The spatial pushforward is the caller's responsibility; this
    primitive only certifies that the angular structure is identity,
    which is a no-op return of the (caller-supplied partner-face
    outgoing) buffer.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(5).standard_normal((quad.N, 3))
    bc = PeriodicBoundary()

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, psi_out)
    # Returned array is a copy (independent of caller's buffer).
    psi_in[0, 0] = 1e9
    assert psi_out[0, 0] != 1e9


@pytest.mark.foundation
def test_albedo_bc_scales_outgoing() -> None:
    """``AlbedoBoundary(α).apply(ψ_out) == α·ψ_out``."""
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_out = np.random.default_rng(6).standard_normal((quad.N, 2))
    bc = AlbedoBoundary(albedo=0.5)

    psi_in = _realize_for_sn(bc, quad).apply(psi_out)

    np.testing.assert_array_equal(psi_in, 0.5 * psi_out)


@pytest.mark.foundation
def test_albedo_zero_and_vacuum_agree_on_inflow_rows() -> None:
    """``AlbedoBoundary(0)`` and vacuum agree on inflow rows.

    Issue #176 / C176.2: migrated from the legacy zeros-all
    equivalence (both BCs returned ``np.zeros_like(psi)`` via the
    2-arg ``apply(psi, quad)`` direct call) to the §16A.5-correct
    inflow-rows equivalence. The realizer's ``AlbedoBoundary(0)``
    branch returns :class:`ZeroOperator` (zeros all rows); the
    vacuum branch returns :class:`IncomingOrdinateMaskTensor`
    (zeros only inflow rows). The two agree on the inflow rows
    (both zero them) — that's the production-relevant
    equivalence for SN sweeps that read inflow only.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_out = np.random.default_rng(7).standard_normal((quad.N, 2))

    albedo_zero = _realize_for_sn(
        AlbedoBoundary(albedo=0.0), quad,
    ).apply(psi_out)
    vacuum_realized = _realize_vacuum_for_face_right(
        VacuumInflow(), quad,
    ).apply(psi_out)

    inflow = np.flatnonzero(quad.mu_x < -1e-12)
    np.testing.assert_array_equal(albedo_zero[inflow], 0.0)
    np.testing.assert_array_equal(vacuum_realized[inflow], 0.0)
    np.testing.assert_array_equal(
        albedo_zero[inflow], vacuum_realized[inflow],
    )


@pytest.mark.foundation
def test_wave0_sum_of_realized_bcs_acts_as_weighted_sum() -> None:
    r"""Wave-11 replacement for ``MixedBoundaryOperator.apply`` linearity.

    The rank-N composer was removed in Wave 11; rank-N compositions are
    now expressed by realising each leaf BC against the method space
    and composing the realised primitives with the Wave-0 algebra
    dunders. The composed operator's ``apply(psi)`` equals the explicit
    weighted sum of the leaves' realised ``apply(psi)`` outputs.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(8).standard_normal((quad.N, 2))
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

    spec_realized = _realize_for_sn(spec, quad)
    white_realized = _realize_for_sn(white, quad)
    composed = 0.3 * spec_realized + 0.7 * white_realized

    psi_composed = composed.apply(psi_out)
    psi_expected = (
        0.3 * spec_realized.apply(psi_out)
        + 0.7 * white_realized.apply(psi_out)
    )

    np.testing.assert_allclose(psi_composed, psi_expected, rtol=1e-14, atol=0.0)


@pytest.mark.foundation
def test_all_primitives_are_resolved_bc() -> None:
    """The runtime-checkable Protocol accepts every primitive."""
    instances: list[BoundaryTraceLaw] = [
        VacuumInflow(),
        ReflectiveBoundary(axis="x", albedo=1.0),
        WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0),
        PeriodicBoundary(),
        AlbedoBoundary(albedo=0.5),
        ZeroFluxBoundary(),
    ]
    for bc in instances:
        assert isinstance(bc, BoundaryTraceLaw), f"{type(bc).__name__} not BoundaryTraceLaw"


# ═══════════════════════════════════════════════════════════════════════
# Issue 9.6 — registry + composition + transpose
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_registry_contains_all_primitives() -> None:
    """All five rank-1 concrete BC subtypes self-register under their key.

    Wave 7 merged the legacy ``BoundaryOperator`` ABC into
    :class:`BoundaryTraceLaw` (Wave O O.4a.1-β retired the deprecated
    alias), so the registry now also holds test
    stubs from :mod:`tests.geometry.test_boundary_trace_law` and the
    Wave-7 additions (``prescribed_inflow``). Wave 11 removed the
    ``"mixed"`` key alongside the ``MixedBoundaryOperator`` class —
    the registry MUST no longer hold it.
    """
    expected_keys = {
        "vacuum", "reflective", "white", "periodic", "albedo",
        "zero_flux",
    }
    registry_keys = set(BoundaryTraceLaw.registry.keys())
    assert expected_keys <= registry_keys
    assert "mixed" not in registry_keys, (
        "Wave 11 removed MixedBoundaryOperator; the 'mixed' key "
        "MUST no longer be in BoundaryTraceLaw.registry."
    )


@pytest.mark.foundation
def test_registry_create_returns_concrete_instance() -> None:
    bc = BoundaryTraceLaw.create("vacuum")
    assert isinstance(bc, VacuumInflow)
    bc = BoundaryTraceLaw.create("reflective", axis="x", albedo=1.0)
    assert isinstance(bc, ReflectiveBoundary)
    assert bc.axis == "x"
    assert bc.albedo == 1.0
    bc = BoundaryTraceLaw.create("zero_flux")
    assert isinstance(bc, ZeroFluxBoundary)


@pytest.mark.foundation
def test_sn_realizer_refuses_zero_flux() -> None:
    """SN REFUSES the zero-flux Dirichlet law (#290 P3 / ruling 3).

    φ_Γ = 0 is the albedo-family member 𝒜 = −1 in the partial-current
    basis; a NEGATIVE angular inflow is unrepresentable in a
    non-negative ψ, so the SN realizer must refuse rather than realize
    something else. The error redirects to the physical zero-incoming
    law (VacuumInflow) and to the diffusion realizer.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    with pytest.raises(BoundaryError, match="VacuumInflow"):
        _realize_for_sn(ZeroFluxBoundary(), quad)


@pytest.mark.foundation
def test_registry_create_unknown_key_raises() -> None:
    with pytest.raises(KeyError):
        BoundaryTraceLaw.create("nonsense")


@pytest.mark.foundation
def test_specular_realized_op_advertises_apply_transpose() -> None:
    """The realised specular op is adjointable.

    Issue #186 (B3 + β2): the structural predicates are an operator-tree
    concept; descriptors do not carry them. The realised
    :class:`PermutationOperator` (α=1 fast path) carries the
    ``is_adjointable`` truth that consumers (sensitivity adjoints) inspect.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    realized = _realize_for_sn(spec, quad)
    assert realized.is_adjointable


@pytest.mark.foundation
def test_specular_apply_transpose_reciprocity_unweighted() -> None:
    r"""``<B(psi_out), phi_in> == <psi_out, B^T(phi_in)>``.

    Tests the Euclidean inner-product reciprocity identity for clean
    axis specular reflection. Index permutations under axis
    reflection are involutions, so the transpose acts as the same
    permutation.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(42)
    psi_out = rng.standard_normal((quad.N, 2))
    phi_in = rng.standard_normal((quad.N, 2))

    spec = ReflectiveBoundary(axis="x", albedo=0.7)
    realized = _realize_for_sn(spec, quad)
    Bpsi = realized.apply(psi_out)
    BTphi = realized.apply_transpose(phi_in)

    lhs = float(np.sum(Bpsi * phi_in))
    rhs = float(np.sum(psi_out * BTphi))
    assert np.isclose(lhs, rhs, rtol=1e-13)


@pytest.mark.foundation
def test_specular_self_inverse_identity() -> None:
    r"""``B(B(x)) == albedo^2 * x`` for a clean axis reflection."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(7)
    x = rng.standard_normal((quad.N, 2))

    spec = ReflectiveBoundary(axis="x", albedo=0.7)
    realized = _realize_for_sn(spec, quad)
    once = realized.apply(x)
    twice = realized.apply(once)

    expected = (0.7 ** 2) * x
    np.testing.assert_allclose(twice, expected, rtol=1e-13, atol=1e-14)


@pytest.mark.foundation
def test_operator_sum_of_bcs_acts_as_weighted_sum() -> None:
    r"""``0.7 * Specular + 0.3 * White`` realises the explicit weighted sum.

    With Issue 9.6, realised BCs participate in operator algebra:
    scalar multiplication and operator addition produce composable
    objects whose ``apply(psi)`` reproduces the explicit weighted sum
    of the leaves' ``apply`` outputs. Wave 11 removed the
    ``MixedBoundaryOperator`` composer; this test pins the algebraic
    identity that replaces the deleted "matches MixedBC baseline"
    contract, comparing the composed operator's ``apply`` against the
    pointwise weighted sum (the structurally-independent reference is
    the linearity of ``LinearOperator.apply`` itself, not any other
    composer).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(99)
    psi_out = rng.standard_normal((quad.N, 2))

    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

    spec_realized = _realize_for_sn(spec, quad)
    white_realized = _realize_for_sn(white, quad)
    composed = 0.7 * spec_realized + 0.3 * white_realized
    expected = 0.7 * spec_realized.apply(psi_out) + 0.3 * white_realized.apply(psi_out)

    np.testing.assert_allclose(
        composed.apply(psi_out), expected, rtol=1e-13, atol=1e-14,
    )


@pytest.mark.foundation
def test_boundary_operator_keys_match_class_attribute() -> None:
    """The ``key`` ClassVar must match the registry insertion key.

    Wave 11 removed ``MixedBoundaryOperator``; only the rank-1
    concretes are pinned here.
    """
    assert VacuumInflow.key == "vacuum"
    assert ReflectiveBoundary.key == "reflective"
    assert WhiteBoundary.key == "white"
    assert PeriodicBoundary.key == "periodic"
    assert AlbedoBoundary.key == "albedo"
