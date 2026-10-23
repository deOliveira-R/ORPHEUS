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
``tests/gates/geometry/test_law_composition.py`` for the tree-walker tests.
"""

from __future__ import annotations

import math
from fractions import Fraction

from tests._harness.float_bounds import mp_fraction

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    PeriodicBoundary,
    BoundaryTraceLaw,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.numerics.spaces.angular_trace_space import build_omega_dot_n
from orpheus.numerics.quadrature import Quadrature


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — protocol semantics on synthetic data
# ═══════════════════════════════════════════════════════════════════════




@pytest.mark.foundation
# ─────────────────────────────────────────────────────────────────────
# RETIRED at wave 2 S2b (retirement-audit D.14 demotion check): the 18 rows below were the law classes' SEMANTICS pinned through
# ``SNBoundaryRealizer``, which the carve deletes; a law class has no realization of its own to test any more. Each claim has a
# new home, in ``tests/gates/sn/operators/test_sn_boundary_realizer.py`` unless named:
#   vacuum is the zero map onto Gamma_-           -> TestRealizeVacuum
#   specular indexes through the partner, y axis, partial weight -> TestRealizeSpecular (hand-computed rows, bit-exact)
#   white: cosine-weighted average, angle-independent, weight linear, 4-point hand computation, z axis -> TestRealizeWhite
#   periodic returns its input unchanged          -> TestRealizeTranslation
#   albedo scales the flux its closure pairs; albedo zero == vacuum (the same zero map) -> TestRealizeVacuum.test_a_wall_of_zero_weight_is_the_same_zero_map
#   the sum of realized walls acts as the weighted sum -> TestRealizeTheWeightedWalls.test_a_mixed_wall_is_the_sum_of_the_images_in_the_canonical_order
#   zero flux is refused                          -> tests/gates/sn/operators/test_b3_domain_narrowing.py::TestEveryLawsDomain::test_a_negative_return_is_refused_on_the_state_cone_axis
#   the specular op advertises apply_transpose; Bt B = alpha^2 -> TestTheRealizedFaceIsMETRICCorrect and TestRealizeSpecular.test_the_scaled_mirror_composed_with_its_transpose_is_alpha_squared
# The rows that stay below pin the law CLASSES (registry, resolved tags, the published GL4 table), until their classes retire at WF2.4.


def test_vacuum_bc_is_resolved_bc() -> None:
    """The Protocol is runtime-checkable: VacuumInflow qualifies."""
    assert isinstance(VacuumInflow(), BoundaryTraceLaw)














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


def _weighted_mean(g: list[Fraction], psi: list[Fraction]) -> Fraction:
    return sum(gi * pi for gi, pi in zip(g, psi)) / sum(g)


#: One unit in the 16th decimal place: the typed entries carry 16 correct
#: decimals (truncated or rounded), the table's OWN precision contract.
#: `[M]` 2026-10-01 against mpmath at 40 digits: the six GL-4/GL-8 entries
#: sit 0.05 to 0.68 units from the exact values, so they are NOT uniformly
#: rounded to 16 decimals (half a unit would red 0.3478548451374538 at 0.57).
_TABLE_UNIT = Fraction(1, 10**16)




@pytest.mark.foundation
def test_the_gl4_rule_is_the_published_table() -> None:
    r"""The PRODUCER claim the white-BC row no longer welds in: the shipped GL-4 rule is the A&S table.

    Two contracts, each stated rather than picked:

    * the table's own: every printed entry carries 16 correct decimals, so it
      sits within one unit of its 16th decimal (:math:`10^{-16}`) of the
      exact node or weight — checked against mpmath at 40 digits,
      independently of production;
    * the producer's: the shipped rule is CORRECTLY ROUNDED (gate (a)), so
      each float is within :math:`\tfrac12\,\mathrm{ulp}` of the exact value.

    Together: :math:`|\text{shipped} - \text{published}| \le
    \tfrac12\mathrm{ulp} + 10^{-16}` per entry. This row is RED on a
    Golub–Welsch rule that is not correctly rounded (8 ULP on :math:`w_2`
    on macOS 27.0.1) and is the producer pin gate (a) already carries; it is
    kept here only so the boundary module still states which table its
    docstring numbers come from. (It may equally be retired into gate (a).)
    """
    import mpmath

    with mpmath.workdps(40):
        nodes, weights = mpmath.gauss_quadrature(4, "legendre")
        exact = sorted(
            ((mp_fraction(abs(x)), mp_fraction(w)) for x, w in zip(nodes, weights) if x > 0),
            key=lambda t: t[0],
        )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    shipped_mu = sorted({abs(float(m)) for m in quad.mu_x})
    shipped_w = {abs(float(m)): float(w) for m, w in zip(quad.mu_x, quad.weights)}
    for (mu_t, w_t), (mu_x, w_x), mu_s in zip(_GL4_TABLE, exact, shipped_mu):
        for printed, true, shipped in ((mu_t, mu_x, mu_s), (w_t, w_x, shipped_w[mu_s])):
            if abs(Fraction(str(printed)) - true) > _TABLE_UNIT:
                pytest.fail(f"the typed table entry {printed} does not carry 16 correct decimals")
            contract = Fraction(math.ulp(shipped)) / 2 + _TABLE_UNIT
            if abs(Fraction(shipped) - Fraction(str(printed))) > contract:
                pytest.fail(
                    f"shipped GL-4 entry {shipped!r} is {float(abs(Fraction(shipped) - Fraction(str(printed)))):.3e} "
                    f"from the published {printed}: outside 1/2 ulp + 1e-16 = {float(contract):.3e}"
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
def test_all_primitives_are_resolved_bc() -> None:
    """The runtime-checkable Protocol accepts every primitive."""
    instances: list[BoundaryTraceLaw] = [
        VacuumInflow(),
        ReflectiveBoundary(axis="x"),
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
    stubs from :mod:`tests.gates.geometry.test_boundary_trace_law` and the
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
    bc = BoundaryTraceLaw.create("reflective", axis="x")
    assert isinstance(bc, ReflectiveBoundary)
    assert bc.axis == "x"
    bc = BoundaryTraceLaw.create("zero_flux")
    assert isinstance(bc, ZeroFluxBoundary)




@pytest.mark.foundation
def test_registry_create_unknown_key_raises() -> None:
    with pytest.raises(KeyError):
        BoundaryTraceLaw.create("nonsense")










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
