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
    IsotropicReturn,
    PeriodicBoundary,
    BoundaryTraceLaw,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.numerics.spaces.angular_trace_space import build_omega_dot_n
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature
from tests._harness.references import mirror_partner_indices
from tests.sn._test_helpers import face_method_space


# ─────────────────────────────────────────────────────────────────────
# Wave 10 helper: migrate 2-arg ``bc.apply(psi, quad)`` legacy call
# sites to the realizer path. See ``transient-giggling-cake.md`` Wave
# 10 for the scope-narrowing decision (the C10.2 "drop 2-arg overload
# + delete _BoundBoundaryOperator" steps are deferred — see the
# follow-up issue filed at end of Wave 10).
# ─────────────────────────────────────────────────────────────────────


def _realize_narrowed_for_face_right(bc, quad):
    r"""Realize a NARROWED BC for the right face (outward normal :math:`+\hat x`).

    Inflow ordinates there are those with :math:`\mu_x < -\epsilon`
    (:math:`\Omega\cdot\hat n_{\rm right} = +\mu_x < 0`); outflow the mirror.

    RENAMED at campaign phase **B3.2** from ``_realize_vacuum_for_face_right``.
    Pre-B3.2 only vacuum needed per-face data (for its inflow mask); since a
    law's DOMAIN is now :math:`\Gamma_+`, **reflective needs a face too** —
    the specular table must be re-indexed into outflow-local coordinates, and
    a faceless space has no :math:`\Gamma_+` to index into.
    """
    return SNBoundaryRealizer().realize(bc, face_method_space(quad, face="xmax"))


def _half_traces(quad):
    r"""The right face's ``(inflow, outflow)`` index sets, hand-derived.

    Written out rather than read back off the method space so the tests below
    keep an independent statement of WHICH ordinates are which — the sign
    convention is exactly what a boundary test exists to pin.
    """
    return (
        np.flatnonzero(quad.mu_x < -1e-12),
        np.flatnonzero(quad.mu_x > +1e-12),
    )


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — protocol semantics on synthetic data
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_vacuum_bc_realizes_to_the_zero_map_onto_gamma_minus() -> None:
    r"""Vacuum BC via the realizer is the zero map :math:`\Gamma_+ \to \Gamma_-`.

    Issue #176 / C176.2 migrated this from the legacy zeros-all
    ``bc.apply(psi, quad)`` contract to the §16A.5 inflow-only mask.
    Campaign phase **B3.2** completes the arc: with the law's DOMAIN narrowed
    to :math:`\Gamma_+`, there are no outflow rows to preserve, and vacuum's
    whole content (:math:`R = 0`) is the honest zero map between the two
    half-traces.

    The retired half of the old assertion — "outflow rows pass through
    unchanged" — is not deleted but SUPERSEDED: those rows are no longer in
    the operator's domain, so the claim has no referent. What replaces it is
    the codomain shape, asserted below.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.random.default_rng(0).standard_normal((quad.N, 3))
    bc = VacuumInflow()
    inflow, outflow = _half_traces(quad)

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_out[outflow])

    assert psi_in.shape == (inflow.size, 3), (
        f"vacuum emitted {psi_in.shape}; the narrowed codomain is Γ₋, i.e. "
        f"{(inflow.size, 3)} — NOT the full face {psi_out.shape}."
    )
    np.testing.assert_array_equal(psi_in, 0.0)


@pytest.mark.foundation
def test_vacuum_bc_is_resolved_bc() -> None:
    """The Protocol is runtime-checkable: VacuumInflow qualifies."""
    assert isinstance(VacuumInflow(), BoundaryTraceLaw)


@pytest.mark.foundation
def test_specular_bc_indexes_through_reflection_partner() -> None:
    r"""``B(γ₊ψ)[j] == ψ[mirror_partner[inflow[j]]]``.

    RE-POSED at **B3.2**: the pre-B3.2 claim was the full-face gather
    ``psi_in == psi_out[mirror_partner]``. Narrowed, row ``j`` of the image
    is the mirror of the ``j``-th INFLOW ordinate — which is the same gather
    restricted to :math:`\Gamma_-`, and is written here from the independent
    geometric partner reference and the inflow set alone (no ``to_local``,
    no production pairing derivation), so neither a remap error nor a
    pairing error can cancel against the reference.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    psi_out = np.arange(quad.N * 2, dtype=float).reshape(quad.N, 2)
    bc = ReflectiveBoundary(axis="x", albedo=1.0)
    ref = mirror_partner_indices(quad, "x")
    inflow, outflow = _half_traces(quad)

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_out[outflow])

    np.testing.assert_array_equal(psi_in, psi_out[ref][inflow])


@pytest.mark.foundation
def test_specular_bc_with_partial_albedo() -> None:
    r"""ReflectiveBoundary scales by ``albedo``, on the narrowed domain.

    ``gauss_legendre(4)`` at ``xmax``: ``ref[n] = N-1-n``, inflow = ``[0, 1]``
    (μ<0), outflow = ``[2, 3]``. So the image is ``0.5 * ψ[[3, 2]]`` — note
    the REVERSAL, which is what makes the slab the discriminating fixture for
    the reflective narrowing's ``to_local`` remap (a naive ``arange`` would
    give ``0.5 * ψ[[2, 3]]``). Hand-written as a literal below so the leg
    states the expected values rather than re-deriving them.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_out = np.array([[1.0], [2.0], [3.0], [4.0]])
    bc = ReflectiveBoundary(axis="x", albedo=0.5)
    inflow, outflow = _half_traces(quad)
    np.testing.assert_array_equal(inflow, [0, 1])
    np.testing.assert_array_equal(outflow, [2, 3])

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_out[outflow])

    np.testing.assert_array_equal(psi_in, np.array([[2.0], [1.5]]))


@pytest.mark.foundation
def test_specular_bc_axis_y_on_lebedev() -> None:
    """Lebedev y-reflection partner: ``apply(axis='y')`` matches index.

    B3.2: realized on the ``ymax`` face (outward normal +ŷ), so the domain is
    that face's :math:`\\Gamma_+` and the reference is the y-reflection gather
    restricted to its :math:`\\Gamma_-`.
    """
    quad = Quadrature.lebedev(order=9)
    psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
    bc = ReflectiveBoundary(axis="y", albedo=1.0)
    space = face_method_space(
        quad, face="ymax", faces=("xmin", "xmax", "ymin", "ymax"),
    )
    # Ω·n_ymax = +μ_y, so inflow ⟺ μ_y < 0 — stated independently of the
    # helper, because a face↔axis mix-up is precisely this test's threat.
    np.testing.assert_array_equal(
        space.inflow_indices, np.flatnonzero(quad.mu_y < -1e-12),
    )

    psi_in = SNBoundaryRealizer().realize(bc, space).apply(
        psi_out[space.outflow_indices]
    )

    np.testing.assert_array_equal(
        psi_in, psi_out[mirror_partner_indices(quad, "y")][space.inflow_indices],
    )


@pytest.mark.foundation
def test_white_bc_returns_cosine_weighted_average() -> None:
    r"""White BC: :math:`\gamma_-\psi = \sum_{\Gamma_+} w\mu\psi / \sum_{\Gamma_+} w\mu`.

    **B3.4a** narrowed the law to :math:`\Gamma_+ \to \Gamma_-`, so the probe
    lives on :math:`\Gamma_+` and the image on :math:`\Gamma_-`. The reference
    stays hand-computed from ``quad.mu_x`` / ``quad.weights`` — but it is now
    summed over the RESTRICTED index set rather than over a masked full-``N``
    array, which is the whole point: the mask was a second outflow classifier.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    inflow, outflow = _half_traces(quad)
    rng = np.random.default_rng(2)
    psi_out = rng.standard_normal((outflow.size, 1))
    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_out)

    cos_w = quad.weights[outflow] * quad.mu_x[outflow]
    expected_avg = (cos_w[:, None] * psi_out).sum(axis=0) / cos_w.sum()
    expected = np.broadcast_to(expected_avg, (inflow.size, 1))

    assert psi_in.shape == (inflow.size, 1)
    np.testing.assert_allclose(psi_in, expected, rtol=1e-15, atol=1e-15)


@pytest.mark.foundation
def test_white_bc_is_angle_independent() -> None:
    """White BC produces the same value at every INFLOW ordinate."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    inflow, outflow = _half_traces(quad)
    psi_out = np.random.default_rng(3).standard_normal((outflow.size, 2))
    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.7)

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_out)

    assert psi_in.shape[0] == inflow.size
    # All rows identical — Lambertian emission has no angular structure.
    for n in range(1, inflow.size):
        np.testing.assert_allclose(psi_in[n], psi_in[0], rtol=1e-15, atol=0.0)


@pytest.mark.foundation
def test_white_bc_albedo_scales_linearly() -> None:
    """Two white BCs at α=1 and α=0.4 differ by exactly the ratio."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    _inflow, outflow = _half_traces(quad)
    psi_out = np.random.default_rng(4).standard_normal((outflow.size, 2))
    bc1 = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    bc2 = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.4)

    psi1 = _realize_narrowed_for_face_right(bc1, quad).apply(psi_out)
    psi2 = _realize_narrowed_for_face_right(bc2, quad).apply(psi_out)

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

    **B3.4a re-pose.** The old body also fed the INCOMING ordinates 7.0, to
    constrain the outgoing-hemisphere mask (they must not enter the average).
    That mask is gone: the law's domain is :math:`\Gamma_+`, so an inflow
    ordinate cannot enter the average because it cannot be handed to the
    operator at all. The constraint is preserved — strictly strengthened — as
    the final leg: a full-face input is REFUSED. An erasure became an absence,
    and an absence is not silently removable the way a mask was.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    inflow, outflow = _half_traces(quad)
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

    psi_1, psi_2 = 1.0, 4.0
    mu_out = quad.mu_x[outflow]
    psi_out = np.empty((outflow.size, 1))
    psi_out[np.isclose(mu_out, +mu_1)] = psi_1
    psi_out[np.isclose(mu_out, +mu_2)] = psi_2
    # Activation guard, OUTSIDE the claim: the two outgoing ordinates must
    # genuinely carry DIFFERENT values, or the weighted mean is invariant
    # under the whole w·|Ω·n̂| formula and the gate is the blind one B0.3
    # repaired (vv Mode 12).
    assert psi_out.min() != psi_out.max()

    bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    op = _realize_narrowed_for_face_right(bc, quad)
    psi_in = op.apply(psi_out)

    expected = (w_1 * mu_1 * psi_1 + w_2 * mu_2 * psi_2) / (
        w_1 * mu_1 + w_2 * mu_2
    )
    assert abs(expected - 2.723973656470134) < 1e-14  # the docstring's number
    assert psi_in.shape == (inflow.size, 1)
    np.testing.assert_allclose(
        psi_in, expected, rtol=1e-14, atol=0.0,
        err_msg=(
            "white BC is not the |Ω·n̂|-weighted outgoing-hemisphere "
            "mean on GL-4"
        ),
    )

    # The re-posed mask constraint: an inflow ordinate cannot enter the
    # average because the whole face is not in the domain.
    with pytest.raises(ValueError, match=r"expected \|Γ₊\|"):
        op.apply(np.full((quad.N, 1), 7.0))


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

    # B3.4a — the probe lives on Γ₊ of the z-max face. The old body wrote 9.0
    # on the incoming ordinates to constrain the mask; those ordinates are no
    # longer in the domain, so the constraint is now structural.
    from orpheus.numerics.spaces.angular_trace_space import (
        TANGENTIAL_EPS, build_omega_dot_n,
    )
    odn_z = build_omega_dot_n(quad, ("zmax",))[0]
    outflow_z = np.flatnonzero(odn_z > +TANGENTIAL_EPS)
    inflow_z = np.flatnonzero(odn_z < -TANGENTIAL_EPS)
    psi_out = quad.mu_z[outflow_z][:, None]
    bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)

    psi_in = SNBoundaryRealizer().realize(
        bc, face_method_space(quad, face="zmax", faces=("zmin", "zmax")),
    ).apply(psi_out)
    assert psi_in.shape == (inflow_z.size, 1)

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
    r"""Slab GL has zero z-cosines — a z-axis white BC is a RANK MISMATCH
    between the face and the cubature, and is refused as one.

    R-1 Phase A detour-C: the unified :class:`Quadrature` class always exposes
    ``mu_z`` as a :func:`@property` view (returning zeros for slab measures),
    so the OLD ``GaussLegendre1D`` adapter's :exc:`AttributeError` on a missing
    ``mu_z`` cannot fire; the defect has to be caught by VALUE.

    **B3.4a moved the catch one layer further in, and sharpened it.** The
    operator no longer runs its own outflow classification, so it no longer
    reaches a "no outgoing ordinates" state of its own; it asks
    :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n` —
    the single face-name → signed-projection primitive — which refuses a face
    whose normal-axis cosines are IDENTICALLY zero, naming the mismatch between
    the face layout and the angular cubature rather than reporting the
    downstream symptom. Same defect, one fewer place that can classify it, and
    a message that says what is actually wrong.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)

    # Re-posed at G6.3 step 3b: ``AngularAverageOperator.from_quadrature``
    # is retired, and the refusal was never its own — this docstring already
    # credits "the single face-name → signed-projection primitive". Asserted
    # directly against that primitive, which is what the realizer calls.
    with pytest.raises(ValueError, match="requires genuine mu_z"):
        build_omega_dot_n(quad, ("zmax",))


@pytest.mark.foundation
def test_periodic_bc_returns_input_unchanged() -> None:
    r"""PeriodicBoundary is the identity on the ordinate index — and that is
    EARNED by the partner identification, not assumed.

    RE-POSED at campaign phase **B3.4c**, the last of the narrowings (vacuum /
    reflective at B3.2, white / prescribed inflow at B3.4a, albedo at B3.4b).
    Two things about this test were false before it:

    * it realized on a FACELESS method space, and its docstring said "the
      spatial pushforward is the caller's responsibility". No caller ever had
      that responsibility, and no mechanism existed for one to discharge it —
      the composition fed every law its OWN face's :math:`\Gamma_+`, so
      periodic returned a face's own outflow as its inflow. Since B3.4c the
      channel is real (:meth:`~orpheus.geometry.boundary.PairedDeck.domain_face` names the partner and
      ``SNBoundaryOperator._face_domains`` supplies it), and a law that cannot
      name its installation face is refused rather than realized;
    * the probe was full-face, so "identity" was a statement about an
      endomorphism. Narrowed, it is a statement about
      :math:`\Gamma_+(f') \to \Gamma_-(f)` — a map between two DIFFERENT index
      sets that happens to be the identity on the global ordinate index
      because the faces' outward normals are opposite. The assertion below
      pins that those two sets really are different-but-equal: disjoint from
      :math:`\Gamma_+(f)`, equal to each other.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    bc = PeriodicBoundary(axis="x")
    inflow, outflow = _half_traces(quad)          # of the RIGHT face (xmax)
    partner_outflow = np.flatnonzero(quad.mu_x < -1e-12)   # Γ₊ of xmin

    # The identification that makes the identity body correct, stated here
    # independently of the realizer's own guard: xmax's Γ₋ IS xmin's Γ₊, and
    # is NOT xmax's own Γ₊ (so the probe below discriminates).
    np.testing.assert_array_equal(inflow, partner_outflow)
    assert not set(inflow) & set(outflow)

    psi_partner_out = np.random.default_rng(5).standard_normal(
        (partner_outflow.size, 3)
    )
    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_partner_out)

    assert psi_in.shape == (inflow.size, 3)
    np.testing.assert_array_equal(psi_in, psi_partner_out)
    # ⚠ The "returned array is a copy" leg that stood here RETIRED at B3.4c
    # with the leaf that guaranteed it. ``PeriodicWrapOperator``'s body was
    # ``x.copy()``; the successor is a bare
    # :class:`~orpheus.numerics.operator.IdentityOperator`, which returns its
    # argument BY REFERENCE — deliberately, because it is algebra-closed
    # (``I.inverse() is I``). So the retired leaf's claimed "project-wide
    # convention that op.apply(psi) may be mutated freely" was never universal:
    # ``IdentityOperator`` had always violated it.
    #
    # The contract that MATTERS is at the composition, and it is pinned there
    # rather than dropped —
    # ``test_sn_boundary_operator.py::TestPeriodicReadsThePartnerFace::
    # test_the_output_never_aliases_the_input``. It holds for a structural
    # reason: the trace restriction is fancy indexing (which copies) and the
    # image is scattered into a freshly-zeroed sink, so no production caller
    # ever hands the leaf a buffer it still owns.


@pytest.mark.foundation
def test_albedo_bc_scales_the_flux_its_CLOSURE_pairs() -> None:
    r"""``AlbedoBoundary(α, SpecularReturn(a))`` returns ``α`` of the MIRROR
    partner's flux — not of the array-position partner's.

    **RE-POSED at B3.4b.** This asserted ``psi_in == 0.5 * psi_out`` on the
    whole face, which was the pre-narrowing endomorphism. That statement is
    now unspellable: a bare albedo is refused, and the completed law is typed
    :math:`\Gamma_+ \to \Gamma_-`, so "α times the outgoing flux" has to say
    WHICH outgoing direction feeds each incoming one.

    The slab is the DISCRIMINATING fixture and that is the whole point of
    keeping this row here rather than folding it into the ≡ theorems: on
    ``gauss_legendre`` the mirror REVERSES order (``perm[inflow] = [3, 2]`` →
    local ``[1, 0]``), so the old positional reading and the specular reading
    genuinely disagree. On ``product(2,4)`` / ``level_symmetric(6)`` they
    coincide and this row would prove nothing (MEASURED 2026-08-01).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_full = np.random.default_rng(6).standard_normal((quad.N, 2))
    inflow, outflow = _half_traces(quad)
    bc = AlbedoBoundary(0.5, SpecularReturn(axis="x"))

    psi_in = _realize_narrowed_for_face_right(bc, quad).apply(psi_full[outflow])

    perm = mirror_partner_indices(quad, "x")
    np.testing.assert_array_equal(psi_in, 0.5 * psi_full[perm[inflow]])
    # …and the retired positional reading is DIFFERENT here, which is what
    # makes the assertion above discriminating rather than a coincidence.
    assert not np.array_equal(psi_in, 0.5 * psi_full[outflow]), (
        "on this quadrature the mirror pairing and the array-position pairing "
        "coincide, so this row cannot tell them apart — pick a slab or "
        "lebedev fixture."
    )


@pytest.mark.foundation
def test_albedo_zero_and_vacuum_are_THE_SAME_narrowed_zero_map() -> None:
    r"""``AlbedoBoundary(0, closure)`` and vacuum are the same object.

    **RE-POSED at B3.4b, and the claim STRENGTHENED.** Before the narrowing
    these two lived on different domains (albedo(0) was a full-face
    endomorphism), so the honest statement was the weaker "they agree on the
    inflow rows". B3.4b's fold routes every α=0 kernel through
    ``_narrowed_zero_operator`` — the same body vacuum uses — so the two are
    now the SAME map :math:`\Gamma_+ \to \Gamma_-`, and that is what this
    asserts: a perfectly absorbing surface IS a vacuum.

    The non-zero probe matters (``vv-principles`` Mode 12): ``|Γ₊| == |Γ₋|``
    on every quadrature in the tree, so shape alone cannot distinguish a
    genuine two-space zero map from an endomorphic ``0.0 * x`` echo.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    psi_full = np.random.default_rng(7).standard_normal((quad.N, 2))
    inflow, outflow = _half_traces(quad)
    probe = psi_full[outflow]

    vacuum_image = _realize_narrowed_for_face_right(VacuumInflow(), quad).apply(probe)
    for closure in (
        SpecularReturn(axis="x"),
        IsotropicReturn(axis="x", outward_sign=+1),
    ):
        albedo_image = _realize_narrowed_for_face_right(
            AlbedoBoundary(0.0, closure), quad,
        ).apply(probe)
        assert albedo_image.shape == (inflow.size, 2)
        np.testing.assert_array_equal(albedo_image, 0.0)
        np.testing.assert_array_equal(albedo_image, vacuum_image)
    # The probe is NOT zero — otherwise every row above is 0 == 0.
    assert np.count_nonzero(probe)


# The two mixed-law gates below carried ``xfail(strict=True)`` from B3.2 until
# **B3.4c**, on the shared reason "a NARROWED law cannot be summed with an
# UN-NARROWED one — white is still a full-face endomorphism". That debt is now
# discharged: B3.4a narrowed white, B3.4b albedo, B3.4c periodic, so every law
# is typed Γ₊ → Γ₋ and ``α·spec + β·white`` has matching factor shapes. The
# markers are DELETED (an xfail that XPASSes is a failure, by design — the
# marker set is the todo list, and this item is done).
#
# ⚠ The markers had also gone MISATTRIBUTED, which is why this is worth a note
# rather than a silent deletion (vv Mode 8, class 4). From B3.4a onward they no
# longer reddened for the documented reason: white WAS narrowed, and what
# actually raised was the faceless ``_realize_for_sn`` helper handing it a
# method space with no Γ₊ — a setup error, not the shape mismatch the reason
# names. An xfail hides ANY failure, so the rows read as committed coverage of
# "narrowing incomplete" while asserting nothing about it. Routing both
# operands through the same face-carrying space is what makes the claim
# stateable AND what makes the pass mean what it says.


@pytest.mark.foundation
def test_wave0_sum_of_realized_bcs_acts_as_weighted_sum() -> None:
    r"""Wave-11 replacement for ``MixedBoundaryOperator.apply`` linearity.

    The rank-N composer was removed in Wave 11; rank-N compositions are
    now expressed by realising each leaf BC against the method space
    and composing the realised primitives with the Wave-0 algebra
    dunders. The composed operator's ``apply(psi)`` equals the explicit
    weighted sum of the leaves' realised ``apply(psi)`` outputs.

    RE-POSED at **B3.2** onto the narrowed domain :math:`\Gamma_+`, which is
    where a mixed BC's summands will BOTH live once B3.4 lands. The claim is
    verbatim what it was; only the vector it is asserted on moved.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

    spec_realized = _realize_narrowed_for_face_right(spec, quad)
    white_realized = _realize_narrowed_for_face_right(white, quad)
    _, outflow = _half_traces(quad)
    psi_out = np.random.default_rng(8).standard_normal((outflow.size, 2))
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
    # The refusal fires BEFORE any method-space read (it is a state-cone
    # refusal, not a missing-data one), so a face-carrying space proves the
    # point more strongly than a faceless one: the law is unrealizable even
    # when everything needed to realize it is present.
    with pytest.raises(BoundaryError, match="VacuumInflow"):
        _realize_narrowed_for_face_right(ZeroFluxBoundary(), quad)


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
    realized = _realize_narrowed_for_face_right(spec, quad)
    assert realized.is_adjointable


@pytest.mark.foundation
def test_specular_apply_transpose_reciprocity_unweighted() -> None:
    r"""``⟨B γ₊ψ, φ⟩_{Γ₋} == ⟨γ₊ψ, Bᵀφ⟩_{Γ₊}`` — the RECTANGULAR reciprocity.

    RE-POSED at **B3.2**. The pre-B3.2 leg paired two full-face vectors, which
    was the right statement for a square endomorphism. A narrowed law is a
    rectangular :math:`\Gamma_+ \to \Gamma_-` map, so its Euclidean adjoint
    pairs a :math:`\Gamma_+` vector with a :math:`\Gamma_-` one — the identity
    is unchanged in content, but each side now lives in its own space, and
    getting THAT wrong is the new failure mode this leg guards.

    Note the deliberate choice of INDEPENDENT random vectors on the two
    half-traces (rather than two slices of one full-face draw): a shared draw
    would let a domain/codomain confusion partially cancel.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(42)
    inflow, outflow = _half_traces(quad)
    psi_plus = rng.standard_normal((outflow.size, 2))
    phi_minus = rng.standard_normal((inflow.size, 2))

    spec = ReflectiveBoundary(axis="x", albedo=0.7)
    realized = _realize_narrowed_for_face_right(spec, quad)
    Bpsi = realized.apply(psi_plus)
    BTphi = realized.apply_transpose(phi_minus)
    assert Bpsi.shape == phi_minus.shape, "B does not land on Γ₋"
    assert BTphi.shape == psi_plus.shape, "Bᵀ does not land on Γ₊"

    lhs = float(np.sum(Bpsi * phi_minus))
    rhs = float(np.sum(psi_plus * BTphi))
    assert np.isclose(lhs, rhs, rtol=1e-13)
    # Non-vacuity: an all-zero pairing would satisfy the identity trivially.
    assert abs(lhs) > 1e-6, f"reciprocity pairing is ~0 ({lhs:.3e}) — vacuous"


@pytest.mark.foundation
def test_specular_narrowed_law_composed_with_its_transpose_is_alpha_squared() -> None:
    r"""``Bᵀ(B(x)) == albedo² · x`` on :math:`\Gamma_+`.

    RE-POSED at **B3.2** from ``B(B(x)) == albedo² · x``. The old spelling is
    no longer a well-formed statement: a :math:`\Gamma_+ \to \Gamma_-` map
    cannot be composed with itself, and it silently "worked" pre-B3.2 only
    because the law was a square endomorphism of the whole face.

    The mathematical content — a specular mirror is an INVOLUTION, so applying
    it twice returns you where you started up to :math:`\alpha^2` — survives
    intact in the honest composition ``ιᵀ ∘ ι``: the narrowed ``B`` is
    ``α · (a bijection Γ₊ → Γ₋)``, so ``Bᵀ B = α² I`` on :math:`\Gamma_+`.
    That is the same claim, correctly typed, and it is STRONGER than the old
    one: it additionally pins that the mirror is a bijection rather than
    merely an involution on a larger space.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(7)
    _, outflow = _half_traces(quad)
    x = rng.standard_normal((outflow.size, 2))

    spec = ReflectiveBoundary(axis="x", albedo=0.7)
    realized = _realize_narrowed_for_face_right(spec, quad)
    round_trip = realized.apply_transpose(realized.apply(x))

    np.testing.assert_allclose(
        round_trip, (0.7 ** 2) * x, rtol=1e-13, atol=1e-14,
    )


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

    RE-POSED at **B3.2** onto :math:`\Gamma_+` — see the sibling above.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    rng = np.random.default_rng(99)
    _, outflow = _half_traces(quad)
    psi_out = rng.standard_normal((outflow.size, 2))

    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)

    spec_realized = _realize_narrowed_for_face_right(spec, quad)
    white_realized = _realize_narrowed_for_face_right(white, quad)
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
