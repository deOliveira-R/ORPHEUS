r"""Universal ``assert_*`` invariant tests (Wave 7 / C7.6).

Each concrete :class:`BoundaryTraceLaw` overrides the relevant
``assert_*`` invariants from the §16A.12 universal-invariants
catalog (Grand Report v3). This file pins the invariant contract
PER CONCRETE BC × PER APPLICABLE INVARIANT.

V&V tags
========

* ``@pytest.mark.l1`` — these tests verify a numerical contract
  (the invariant holds on physically valid inputs and raises on
  physically invalid ones), not a software invariant. They pair
  with the matching ERR-NNN error in
  ``docs/theory/verification/error_catalog.rst``.
* ``@pytest.mark.catches("ERR-NNN")`` decorators wire each negative
  test to its catalogue entry per the ERR-040..ERR-047 mapping
  shipped in Wave 3.

Invariant → ERR-NNN map
=======================

* :meth:`ReflectiveBoundary.assert_is_involutive` → ERR-044
  (``ReflectionNotInvolutiveError``)
* :meth:`ReflectiveBoundary.assert_geometry_map_measure_preserving`
  → ERR-042 (``BoundaryGeometryMapNotMeasurePreservingError``) — the
  DIRECT ``w·|μ|`` pushforward check (#52; the pre-#52 delegation to
  the involution check left the unequal-weight-pairing hole open)
* :meth:`ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow`
  → ERR-045 (``ReflectionDidNotMapInflowToOutflowError``) — the
  third independent reflection-table invariant (#52)
* :meth:`BoundaryTraceLaw.assert_source_is_placeable` →
  ERR-047 (``BoundarySourceNotOnIncomingTraceError``) — the
  universal certification that a source-carrying law has a NAMEABLE
  :math:`\\Gamma_-` to deliver into (#52; real body on the ABC;
  RE-POSED at P6 from ``assert_source_lives_on_incoming_trace``,
  which certified the source's VALUES through a declinable probe)
* :meth:`WhiteBoundary.assert_response_positive_if_declared` →
  ERR-043 (``BoundaryResponseNotPositiveError``)
* :meth:`WhiteBoundary.assert_submarkov` → ERR-046
  (``SubmarkovViolationError``)
* :meth:`AlbedoBoundary.assert_response_positive_if_declared` →
  ERR-043
* :meth:`AlbedoBoundary.assert_submarkov` → ERR-046

Production wiring (#52): every invariant above fires at REALIZE time
through :meth:`BoundaryTraceLaw.assert_realizable` (called by
:meth:`SNBoundaryRealizer.realize` at entry), so the negative legs
below include realize-time catchers — the invariants are production
guards, not test-only helpers. The ERR-041 realizer-seam guard
(vacuum trace orientation) lives with the realizer's own tests in
``tests/gates/sn/operators/test_sn_boundary_realizer.py``.

Plan: ``.claude/plans/transient-giggling-cake.md`` Wave 7 / C7.6;
#52 (the four unbuilt BC-layer invariants, 2026-07-12).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectiveBoundary,
    SelfPairedDeck,
    SubmarkovViolationError,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.geometry.boundary._errors import ReflectionNotInvolutiveError
from orpheus.geometry.transformation import Permutation
from orpheus.geometry.faces import AXIS_NAMES
from orpheus.numerics.quadrature import Quadrature
# RETIRED at wave 2 S2b (retirement-audit D.14 demotion check). Nine tests that ran a law through ``SNBoundaryRealizer`` (deleted by
# the carve) are gone; the 41 rows that stay pin the law CLASSES until they retire at WF2.4. Where each deleted claim lives:
#   the three ``test_raises_at_realize_time`` rows (ERR-044 involution, ERR-042 measure preservation, ERR-045 inflow-to-outflow on the
#     poisoned tables) and ``test_raises_at_realize_time_on_the_albedo_route`` (the certification on the second carrier) ->
#     ``tests/gates/numerics/test_ordinate_pairing.py`` (the certification is numerics', keyed on the quadrature and the motion; the
#     poison fixtures ``_TableQuad`` and ``_BROKEN`` stay here until WF2.4, WF2.2 rider (a));
#   ``TestSourceIsPlaceableInvariant.test_raises_at_realize_time`` (ERR-047, a source that cannot be sampled on the face) ->
#     ``tests/gates/transport/test_boundary_source_from_specs.py`` (the ERR-047 marker migrated there at P3: the SN source realizer
#     refuses a space that is not the inflow half-trace);
#   ``TestPrescribedInflowRealizesTheZeroMap`` (4) and ``TestSNRealizerPrescribedInflowDispatch`` (3): a pure inflow realizes to the
#     zero map, linear, source-blind, with a refused faceless space -> ``test_operator_block_role.py::TestBoundaryLeaves::
#     test_prescribed_inflow_realizes_the_same_object_vacuum_does``, ``test_capability_survival.py::
#     TestRealizedPrescribedInflowCapabilities``, ``test_declared_law_is_linear.py`` and ``test_sn_boundary_realizer.py::TestRealizeVacuum``;
#   ``TestFactorAdjointabilityMatchesTheRealizedOperator`` (3 of 4): the law factors' adjointability against the realized operator's
#     is a question about a stored string kind; the typed operator DERIVES ``is_adjointable`` from its inner operator
#     (``RealizedBoundary.is_adjointable``), gated by ``test_capability_survival.py::test_boundary_law_wrapper_delegates_predicates``.


# ─────────────────────────────────────────────────────────────────────
# ReflectiveBoundary.assert_is_involutive (ERR-044)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-044")
class TestReflectiveInvolutionInvariant:
    """:meth:`ReflectiveBoundary.assert_is_involutive` correctness."""

    def test_passes_for_default_axis_reflection_lebedev(self) -> None:
        """Lebedev-17's x-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="x")
        quad = Quadrature.lebedev(17)
        # No raise.
        bc.assert_is_involutive(quad)

    def test_passes_for_level_symmetric_y_axis(self) -> None:
        """Level-symmetric S6 y-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="y")
        quad = Quadrature.level_symmetric(sn_order=6)
        bc.assert_is_involutive(quad)

    def test_passes_for_gauss_legendre_x_axis(self) -> None:
        """1-D Gauss-Legendre x-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="x")
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        bc.assert_is_involutive(quad)

    def test_raises_for_non_involutive_perm(self) -> None:
        """An injected non-involutive pairing raises.

        Construct a fake quadrature whose ``ordinate_permutation`` (the
        certification's read since G6.3 step 7d) returns a 1-cycle
        rotation (NOT its own inverse). The invariant must catch this —
        a rotation is a genuine bijection, so it passes the
        :class:`Permutation` construction and only THIS check sees it.
        """
        quad = Quadrature.lebedev(17)
        N = quad.n_ordinates

        # A fake pairing that is a rotation (i → (i+1) % N) — not an
        # involution for N ≠ 2.
        class _FakeQuad:
            n_ordinates = quad.n_ordinates

            @staticmethod
            def ordinate_permutation(motion) -> Permutation:  # noqa: ARG004
                return Permutation(np.roll(np.arange(N), 1))

        bc = ReflectiveBoundary(axis="x")
        with pytest.raises(ReflectionNotInvolutiveError):
            bc.assert_is_involutive(_FakeQuad())  # type: ignore[arg-type]


class _TableQuad:
    """Quadrature stand-in: REAL nodes/weights, injectable specular
    pairing.

    The negative legs below need pairings that are wrong in exactly ONE
    invariant while every other datum stays production-real — the
    surgical way to prove the three pairing invariants are
    independent (the ERR-045 catalog lesson: "checking only one or two
    leaves a hole"). Served through ``ordinate_permutation`` — the
    certification's read since G6.3 step 7d — so every mutant must be a
    genuine bijection (:class:`Permutation` refuses anything else at
    construction), which the three mutant families all are: each
    violates its ONE invariant while remaining a permutation.
    """

    def __init__(self, base: Quadrature, table: np.ndarray) -> None:
        self.n_ordinates = base.n_ordinates
        self.weights = base.weights
        self._base = base
        self._table = np.asarray(table)

    @property
    def measure(self):
        # The surrogate mutates the PERMUTATION only; its nodes are the real
        # rule's, so the measure passes straight through. Needed because the
        # reference embeds via ``quad.measure.support.ambient_points(...)`` — a duck-typed
        # stand-in owes every member of the contract it stands in for, not
        # only the ones the last reader happened to call (``plan-authoring``
        # §6b: a surrogate is a call site no symbol grep can return).
        return self._base.measure

    def axis_cosines(self, axis_index: int) -> np.ndarray:
        return self._base.axis_cosines(axis_index)

    def mean_axis_cosine(self, axis_index: int) -> np.ndarray:
        return self._base.mean_axis_cosine(axis_index)

    def ordinate_permutation(self, motion) -> Permutation:  # noqa: ARG002
        return Permutation(self._table)


@pytest.mark.l1
@pytest.mark.catches("ERR-042")
class TestReflectiveMeasurePreservingInvariant:
    """:meth:`ReflectiveBoundary.assert_geometry_map_measure_preserving`
    checks the ``w·|μ|`` pushforward DIRECTLY (#52).

    The pre-#52 body delegated to the involution check on the claim
    that weight symmetry is "implied by construction" — leaving open
    the exact ERR-042 hole: an involutive table pairing ordinates from
    DIFFERENT weight classes. The neighbor-pair mutant below is that
    hole made flesh; its raise is the invariant's tooth, and its clean
    pass through the involution check is the independence proof.
    """

    def test_passes_for_real_quadratures(self) -> None:
        for quad in (
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.product(n_mu=8, n_phi=8),
        ):
            for axis in ("x", "y", "z"):
                ReflectiveBoundary(axis=axis).assert_geometry_map_measure_preserving(quad)

    def test_raises_for_involutive_weight_class_mispairing(self) -> None:
        """Neighbor-pair table on GL-8: a perfect involution that
        pairs ordinates with UNEQUAL weights and |μ| — the measure
        check must red where the involution check stays green."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.array([1, 0, 3, 2, 5, 4, 7, 6]))
        bc = ReflectiveBoundary(axis="x")

        # Independence: the mutant sails through the ERR-044 check…
        bc.assert_is_involutive(mutant)  # type: ignore[arg-type]

        # …and reds the ERR-042 check.
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
            bc.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]



# ─────────────────────────────────────────────────────────────────────
# ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow (ERR-045)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-045")
class TestReflectiveInflowToOutflowInvariant:
    """:meth:`ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow`
    (#52) — the THIRD independent reflection-table invariant.

    The identity table is the canonical ERR-045 mutant: every ordinate
    self-maps, so it is trivially involutive (ERR-044 green) and
    trivially measure-preserving (ERR-042 green — ``m[id] ≡ m``), yet
    every non-tangential inflow ordinate maps to ITSELF instead of its
    outflow partner — a length-1 self-loop in the sweep dependency
    graph. Only this check can see it.
    """

    def test_passes_for_real_quadratures(self) -> None:
        for quad in (
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.product(n_mu=8, n_phi=8),
        ):
            for axis in ("x", "y", "z"):
                ReflectiveBoundary(axis=axis).assert_reflection_maps_inflow_to_outflow(quad)

    def test_identity_table_raises_while_sibling_invariants_pass(self) -> None:
        """The independence triad on ONE mutant: identity passes
        involution AND measure, reds ONLY inflow→outflow."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.arange(quad.n_ordinates))
        bc = ReflectiveBoundary(axis="x")

        bc.assert_is_involutive(mutant)  # type: ignore[arg-type]
        bc.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]

        with pytest.raises(ReflectionDidNotMapInflowToOutflowError):
            bc.assert_reflection_maps_inflow_to_outflow(mutant)  # type: ignore[arg-type]

    def test_tangential_self_map_is_exempt(self) -> None:
        """A 1-D quadrature's y-axis table is the identity over
        ALL-tangential ordinates (``mu_y ≡ 0``) — the geometrically
        correct image. The exemption must keep it green (an
        over-strict check would red every degenerate axis)."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        pi = quad.ordinate_permutation(
            SelfPairedDeck.mirror(axis="y", dimension=3).motion
        )
        assert pi is not None
        np.testing.assert_array_equal(pi.indices, np.arange(quad.n_ordinates))
        ReflectiveBoundary(axis="y").assert_reflection_maps_inflow_to_outflow(quad)



# ─────────────────────────────────────────────────────────────────────
# WhiteBoundary.assert_response_positive_if_declared (ERR-043)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-043")
class TestWhiteResponsePositiveInvariant:
    """:meth:`WhiteBoundary.assert_response_positive_if_declared`."""

    def test_passes_for_zero_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.0).assert_response_positive_if_declared()

    def test_passes_for_unit_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0).assert_response_positive_if_declared()

    def test_passes_for_partial_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.5).assert_response_positive_if_declared()

    def test_raises_for_negative_albedo(self) -> None:
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=-0.1)
        with pytest.raises(BoundaryResponseNotPositiveError):
            bc.assert_response_positive_if_declared()


# ─────────────────────────────────────────────────────────────────────
# WhiteBoundary.assert_submarkov (ERR-046)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-046")
class TestWhiteSubmarkovInvariant:
    """:meth:`WhiteBoundary.assert_submarkov`."""

    def test_passes_for_zero_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.0).assert_submarkov()

    def test_passes_for_unit_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0).assert_submarkov()

    def test_passes_for_partial_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.5).assert_submarkov()

    def test_raises_for_supercritical_albedo(self) -> None:
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.2)
        with pytest.raises(SubmarkovViolationError):
            bc.assert_submarkov()


# ─────────────────────────────────────────────────────────────────────
# AlbedoBoundary.assert_response_positive_if_declared (ERR-043)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-043")
class TestAlbedoResponsePositiveInvariant:
    """:meth:`AlbedoBoundary.assert_response_positive_if_declared`."""

    def test_passes_for_zero_albedo(self) -> None:
        AlbedoBoundary(albedo=0.0).assert_response_positive_if_declared()

    def test_passes_for_unit_albedo(self) -> None:
        AlbedoBoundary(albedo=1.0).assert_response_positive_if_declared()

    def test_passes_for_partial_albedo(self) -> None:
        AlbedoBoundary(albedo=0.5).assert_response_positive_if_declared()

    def test_raises_for_negative_albedo(self) -> None:
        with pytest.raises(BoundaryResponseNotPositiveError):
            AlbedoBoundary(albedo=-0.1).assert_response_positive_if_declared()


# ─────────────────────────────────────────────────────────────────────
# AlbedoBoundary.assert_submarkov (ERR-046)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-046")
class TestAlbedoSubmarkovInvariant:
    """:meth:`AlbedoBoundary.assert_submarkov`."""

    def test_passes_for_zero_albedo(self) -> None:
        AlbedoBoundary(albedo=0.0).assert_submarkov()

    def test_passes_for_unit_albedo(self) -> None:
        AlbedoBoundary(albedo=1.0).assert_submarkov()

    def test_passes_for_partial_albedo(self) -> None:
        AlbedoBoundary(albedo=0.5).assert_submarkov()

    def test_raises_for_supercritical_albedo(self) -> None:
        with pytest.raises(SubmarkovViolationError):
            AlbedoBoundary(albedo=1.2).assert_submarkov()


# ReflectiveBoundary's albedo bounds (ERR-043, ERR-046): the class
# ``TestReflectiveAlbedoBounds`` stood here from 2026-09-30 until the
# reflective cleanup (2026-10-01) removed the mirror's albedo, and with it
# the two checks the class pinned. Its successor is
# ``tests/gates/geometry/test_reflective_is_a_mirror.py`` (the albedo is not
# a parameter); ERR-043 and ERR-046 keep their catchers above (white and
# albedo).


# ─────────────────────────────────────────────────────────────────────
# Wave-7 sweep-cycle flag (§15A.2)
# ─────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────
# PrescribedInflow (Wave 7 / C7.5)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
class TestPrescribedInflowConstruction:
    """:class:`PrescribedInflow` construction + registry + sweep-cycle
    classification."""

    def test_default_source_is_no_source(self) -> None:
        from orpheus.geometry.boundary import NoSource, PrescribedInflow

        bc = PrescribedInflow()
        assert isinstance(bc.source, NoSource)

    def test_constant_inflow_source_attaches(self) -> None:
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        bc = PrescribedInflow(source=ConstantInflowSource(value=2.5))
        assert isinstance(bc.source, ConstantInflowSource)
        assert bc.source.value == 2.5

    def test_kind_is_prescribed_inflow(self) -> None:
        from orpheus.geometry.boundary import PrescribedInflow

        assert PrescribedInflow().kind == "prescribed_inflow"

    def test_registered_under_prescribed_inflow_key(self) -> None:
        from orpheus.geometry.boundary import (
            BoundaryTraceLaw,
            PrescribedInflow,
        )

        assert (
            BoundaryTraceLaw.registry["prescribed_inflow"] is PrescribedInflow
        )

    def test_create_via_registry(self) -> None:
        from orpheus.geometry.boundary import (
            BoundaryTraceLaw,
            PrescribedInflow,
        )

        bc = BoundaryTraceLaw.create("prescribed_inflow")
        assert isinstance(bc, PrescribedInflow)




# ``TestIncomingSourceOperator`` lived here until **P3** (2026-08-05): three
# rows exercising the retired ``IncomingSourceOperator`` standalone. Where each
# claim went:
#
# * ``test_apply_ignores_input`` — the rank-0 contract. RETIRED WITH ITS
#   SUBJECT; the successor is not rank-0, it is zero, and
#   ``TestPrescribedInflowRealizesTheZeroMap`` above states that directly.
# * ``test_predicates_are_apply_only`` — INVERTED (the zero map advertises a
#   transpose) and re-posed in
#   ``tests/gates/sn/operators/test_capability_survival.py::TestRealizedPrescribedInflowCapabilities``.
# * ``test_apply_fills_the_codomain_not_the_input_shape`` — the ``vv`` Mode-12
#   row, and the only one that needed a NEW home rather than a rewrite. Its
#   argument was that ``|Γ₊| == |Γ₋|`` on every reachable face, so only a
#   hand-built operator with unequal ends can tell "emits the codomain" from
#   "echoes the input". That hazard did not retire with the operator — it is
#   exactly what ``ZeroOperator``'s two space hooks exist to prevent, and
#   ``_narrowed_zero_operator``'s own docstring calls relying on the
#   endomorphic ``0.0 * x`` echo "wrong in principle and merely lucky in
#   practice". It moved to ``tests/gates/numerics/test_zero_operator_spaces.py``,
#   which is where the unequal-size construction is reachable.




# ─────────────────────────────────────────────────────────────────────
# BoundaryTraceLaw.assert_source_is_placeable (ERR-047)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-047")
class TestSourceIsPlaceableInvariant:
    """:meth:`BoundaryTraceLaw.assert_source_is_placeable`
    (#52, RE-POSED at P6) — the real universal body replacing the no-op default.

    ⭐ **What this class asserts CHANGED at P6, and the change is the point.**
    It used to certify the source's VALUES: an :class:`InflowSourceSpec` filled
    whatever *shape* it was handed, so the delivered-:math:`q` guarantee needed
    the realizer's inflow mask, and the body probed the source at ``(N,)`` to
    see whether :math:`q \\not\\equiv 0` before demanding one.

    That probe was **declinable** — `[M]` a source returning zeros at the probe
    shape and ``7.0`` at the delivery shape skipped the certification and
    delivered anyway. It is gone. A spec now receives :math:`\\Gamma_-(f)`
    itself and returns its shape, so :math:`q \\in \\Gamma_-` holds **by
    construction** and there is no value-level claim left to make.

    What survives is STRUCTURAL and is asserted unconditionally: a law carrying
    a source needs a face that can NAME its inflow set. Note the consequence
    the rows below pin — the check can no longer look at the source's values,
    because in the arm where it would matter (``inflow_indices is None``) there
    is no space to evaluate against. So it discriminates on whether a source is
    :class:`NoSource`, not on what a source currently contains."""

    def test_nonzero_source_without_inflow_set_raises(self) -> None:
        """Law-level negative: the exact catalog catcher —
        ``ConstantInflowSource(2.0)`` on a mixed inflow/outflow
        quadrature, no outflow mask available."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        with pytest.raises(BoundarySourceNotOnIncomingTraceError):
            law.assert_source_is_placeable(None)

    def test_nonzero_source_with_inflow_set_passes(self) -> None:
        """Positive: with the face's inflow indices supplied, the
        realizer masks and the invariant certifies."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        law.assert_source_is_placeable(np.flatnonzero(quad.mu_x < 0))

    def test_homogeneous_laws_certify_masklessly(self) -> None:
        """Positive sweep: every homogeneous law's ``NoSource`` is
        :math:`q \\equiv 0` — trivially on :math:`\\Gamma_-`, with or
        without a mask. The universal body must not demand trace data
        it does not need."""
        for law in (
            VacuumInflow(),
            ReflectiveBoundary(axis="x"),
            WhiteBoundary(axis="x", outward_sign=+1),
            AlbedoBoundary(albedo=0.5),
            PeriodicBoundary(),
        ):
            law.assert_source_is_placeable(None)

    def test_a_zero_VALUED_constant_is_still_a_source_and_still_needs_a_trace(
        self,
    ) -> None:
        """⭐ RE-POSED at P6 — this row asserted the OPPOSITE until 2026-08-06.

        It read ``test_zero_valued_constant_certifies_masklessly``:
        ``ConstantInflowSource(0.0)`` has zero support, so it certified
        without an inflow set, exactly as :class:`NoSource` does.

        That was the presence probe speaking, and the probe is what made the
        guard declinable. The check can no longer ask "is this source's value
        zero?" — in the arm where the answer would matter there is no
        :math:`\\Gamma_-` to evaluate the source against, which is precisely
        what is being complained about. So the discriminator is the source's
        IDENTITY: :class:`NoSource` means *no delivery*, and anything else
        means *a delivery whose destination must be nameable*, whatever it
        currently holds.

        ⚠ This is a deliberate behaviour change and it is strictly safer: the
        old arm let a source that is zero **today** past a check it would fail
        the moment its value changed, and value-at-check-time is not
        value-at-delivery-time — which was ERR-047's whole mechanism.

        The escape hatch for a genuinely sourceless prescribed law is
        unchanged and typed: ``PrescribedInflow()`` defaults to
        :class:`NoSource` and takes the first arm."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            NoSource,
            PrescribedInflow,
        )

        with pytest.raises(BoundarySourceNotOnIncomingTraceError):
            PrescribedInflow(
                source=ConstantInflowSource(value=0.0)
            ).assert_source_is_placeable(None)

        # …and the sourceless spelling still certifies masklessly, which is
        # what keeps the row above a claim about DELIVERY rather than a blanket
        # refusal of the prescribed family.
        assert PrescribedInflow().assert_source_is_placeable(None) is None
        assert isinstance(PrescribedInflow().source, NoSource)


    # ``test_delivered_q_has_no_row_off_the_incoming_trace`` and
    # ``test_codomain_size_is_the_construction_contract`` lived here until
    # **P3** (2026-08-05). Both asked the REALIZED OPERATOR to deliver ``q``
    # and then checked where the rows landed. The operator no longer delivers
    # ``q`` — that was the P3 defect — so the postcondition moved with the
    # delivery:
    #
    # * the delivered-``q``-lands-on-``Γ₋`` leg is now
    #   ``tests/gates/transport/test_boundary_source_from_specs.py::
    #   TestTheRecipeLandsOnTheInflowSlotsOnly``, which carries the
    #   ``catches("ERR-047")`` marker with it. Its
    #   ``test_every_other_row_of_that_face_is_zero`` row is the same claim
    #   against the channel that now delivers.
    # * ``test_codomain_size_is_the_construction_contract`` guarded
    #   ``n_inflow`` — a row count the operator stored ALONGSIDE the space it
    #   describes. P3 removed the count; the space is now the only source, and
    #   ``_narrowed_zero_operator`` runs ``checked_space_extent`` against the
    #   half-trace index arrays it sizes the hooks from, so a space that
    #   contradicts the count is refused at construction rather than asserted
    #   about afterwards.
    #
    # The LAW-tier rows above keep this class's ERR-047 marker earned: they
    # exercise ``assert_source_is_placeable`` itself, which is the
    # invariant the catalog entry names and which P3 did not touch.


# ─────────────────────────────────────────────────────────────────────
# B3.4b — ONE certification, TWO carriers of the specular pairing
# ─────────────────────────────────────────────────────────────────────


def _non_involutive_but_otherwise_valid(
    quadrature: Quadrature, *, axis: str,
) -> np.ndarray:
    r"""A table that breaks ONLY the involution (ERR-044).

    Why not the obvious ``np.roll(arange(N), 1)``: an n-cycle rotation breaks
    all three invariants at once (it pairs unequal weight classes, so ERR-042
    fires FIRST in the aggregate and the row proves nothing about ERR-044's
    independence). Worse, **no pure 3-cycle can ever be the ERR-044-only
    mutant**: a cycle :math:`a \to b \to c \to a` needs
    :math:`\mathrm{sign}(b) = -\mathrm{sign}(a)`,
    :math:`\mathrm{sign}(c) = -\mathrm{sign}(b) = \mathrm{sign}(a)` and then
    :math:`\mathrm{sign}(a) = -\mathrm{sign}(c) = -\mathrm{sign}(a)` — a
    contradiction, so every odd cycle breaks ERR-045 too.

    The construction that works is :math:`\pi \circ \sigma`, where :math:`\pi`
    is the TRUE mirror and :math:`\sigma` swaps two ordinates of the SAME sign
    and the SAME measure :math:`w|\mu_a|`:

    * **measure** — :math:`m[\pi[\sigma[i]]] = m[\sigma[i]] = m[i]`, since
      :math:`\sigma` stays inside a measure class and :math:`\pi` preserves it;
    * **sign** — :math:`\mathrm{sign}(\pi[\sigma[i]]) = -\mathrm{sign}(i)`,
      since :math:`\sigma` preserves sign and :math:`\pi` flips it;
    * **involution** — BROKEN: for the swapped pair :math:`\{a, b\}`,
      :math:`T[T[a]] = \pi[\sigma[\pi[b]]] = \pi[\pi[b]] = b \neq a`
      (:math:`\sigma` does not touch :math:`\pi[b]`, which has the opposite
      sign).
    """
    true_pi = quadrature.ordinate_permutation(
        SelfPairedDeck.mirror(axis=axis, dimension=3).motion
    )
    if true_pi is None:  # pragma: no cover — fixtures are mirror-closed
        raise AssertionError(f"{quadrature!r} has no {axis}-mirror pairing")
    perm = true_pi.indices
    # ``mean_axis_cosine``: the w|mu| cosine measure is a FLUX question, so a
    # suppressed axis contributes zero rather than having no answer. Mirrors
    # the production spelling in ``geometry/boundary/_specular.py``.
    mu = quadrature.mean_axis_cosine(AXIS_NAMES.index(axis))
    measure = quadrature.weights * np.abs(mu)
    # Group same-sign ordinates by measure; the first class with ≥2 members
    # supplies the swap. Rounded because the class is a float coincidence of
    # a symmetric node set, not an exact-compare question.
    classes: dict[tuple[int, float], list[int]] = {}
    for n in range(int(quadrature.n_ordinates)):
        classes.setdefault(
            (int(np.sign(mu[n])), round(float(measure[n]), 12)), []
        ).append(n)
    for members in classes.values():
        if len(members) >= 2:
            a, b = members[0], members[1]
            sigma = np.arange(int(quadrature.n_ordinates))
            sigma[a], sigma[b] = b, a
            return perm[sigma]
    raise AssertionError(
        f"{quadrature!r} has no same-sign measure class with two members, so "
        f"the ERR-044-only mutant is not constructible here — pick a "
        f"quadrature with a degenerate measure class."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-042", "ERR-044", "ERR-045")
class TestSpecularPairingCertifiedOnBothCarriers:
    r"""The three table invariants fire for BOTH laws that stand on the pairing.

    Campaign phase **B3.4b** moved the specular pairing's certification out of
    :class:`ReflectiveBoundary`'s methods and into
    :mod:`~orpheus.geometry.boundary._specular`, because the 2026-08-01 ruling
    gave :class:`AlbedoBoundary` a :class:`SpecularReturn` closure that realizes
    through the SAME specular pairing — with the
    pairing in :math:`R` instead of :math:`G`.

    That module's central claim is *"one certification, two laws"*, and until
    this class it was a **docstring claim with no test**: MEASURED
    2026-08-01, deleting the
    ``isinstance(self.reemission, SpecularReturn)`` clause from
    ``AlbedoBoundary.assert_realizable`` left the entire B3.4b gate file
    (``test_reemission_closure.py``, 178 tests) GREEN. A wrong table would have
    been caught on the reflective route and **silently realized** on the albedo
    one — exactly the twin-path failure this campaign exists to remove.

    Structure: each broken table is wrong in exactly ONE invariant (the
    independence the ERR-045 lesson demands), and each row asserts BOTH
    carriers raise the SAME error type with their OWN ``law=`` attribution.
    The attribution leg is what makes this more than a duplicate of the three
    classes above — it proves the shared module reports the caller, not a
    hard-coded ``"reflective"``.
    """

    #: ``(id, quadrature_factory, table_factory, error_type)``. Each table is
    #: production-real except for the ONE property it breaks.
    #:
    #: The QUADRATURE is per-row, not shared, and that is load-bearing: the
    #: ERR-044-only mutant needs a **degenerate measure class** (two same-sign
    #: ordinates carrying equal ``w·|μ_x|``) to swap inside, and MEASURED,
    #: ``gauss_legendre(8)`` has none — its four per-side measures are all
    #: distinct. ``level_symmetric(6)`` supplies the degeneracy.
    _BROKEN = [
        # Neighbour-pair on GL-8: a perfect involution pairing UNEQUAL weight
        # classes — ERR-042 only (the involution and sign checks stay green).
        ("weight_class_mispairing",
         lambda: Quadrature.gauss_legendre(n_ordinates=8),
         lambda q: np.array([1, 0, 3, 2, 5, 4, 7, 6]),
         BoundaryGeometryMapNotMeasurePreservingError),
        # The identity table: every ordinate self-maps — trivially involutive
        # and trivially measure-preserving, but every inflow ordinate maps to
        # ITSELF instead of its outflow partner. ERR-045 only.
        ("self_map_identity",
         lambda: Quadrature.gauss_legendre(n_ordinates=8),
         lambda q: np.arange(q.n_ordinates),
         ReflectionDidNotMapInflowToOutflowError),
        # The true mirror composed with a swap of two SAME-SIGN, SAME-MEASURE
        # ordinates: measure-preserving and sign-crossing, but NOT its own
        # inverse. ERR-044 only. See :func:`_non_involutive_but_otherwise_valid`
        # for why the obvious n-cycle will not do.
        ("mirror_after_same_class_swap",
         lambda: Quadrature.level_symmetric(sn_order=6),
         lambda q: _non_involutive_but_otherwise_valid(q, axis="x"),
         ReflectionNotInvolutiveError),
    ]

    #: The two carriers, with the tier each holds the pairing in.
    def _carriers(self):
        from orpheus.geometry.boundary import SpecularReturn

        return [
            ("reflective", ReflectiveBoundary(axis="x"), "G"),
            ("albedo", AlbedoBoundary(0.5, SpecularReturn(axis="x")), "R"),
        ]

    @pytest.mark.parametrize(
        "table_id,build_quad,build_table,error_type", _BROKEN,
        ids=[row[0] for row in _BROKEN],
    )
    def test_both_carriers_reject_the_same_broken_table(
        self, table_id, build_quad, build_table, error_type,
    ) -> None:
        """A table broken in ONE invariant reddens BOTH laws' certification."""
        quad = build_quad()
        mutant = _TableQuad(quad, build_table(quad))
        seen = {}
        for law_id, law, tier in self._carriers():
            with pytest.raises(error_type) as exc:
                law.assert_realizable(mutant)  # type: ignore[arg-type]
            seen[law_id] = exc.value.law
        assert seen == {"reflective": "reflective", "albedo": "albedo"}, (
            f"{table_id}: each carrier must attribute the failure to ITSELF; "
            f"got {seen}. A hard-coded law_key in "
            f"orpheus.geometry.boundary._specular would blame the wrong law "
            f"and send a reader to the wrong file."
        )

    def test_the_real_table_certifies_on_both_carriers(self) -> None:
        """POSITIVE control (``vv-principles`` anti-#11).

        Without it the negatives above validate only the RAISING behaviour,
        not the invariant claim: a certification that raised on every table
        would pass all three rows.
        """
        for quad in (
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
        ):
            for law_id, law, tier in self._carriers():
                law.assert_realizable(quad)

    def test_the_albedo_carrier_leaves_the_G_HOOK_a_no_op(self) -> None:
        r"""The TIER claim, not just the check-fires claim.

        ``AlbedoBoundary``'s :math:`G` **is** the identity map and **does**
        preserve the measure, so the base template's polymorphic
        ``assert_geometry_map_measure_preserving`` hook is correctly a no-op
        here even on a table that is catastrophically broken — the pairing is
        certified as a law-specific extension because it lives in :math:`R`.

        This is the gate that reds if someone "closes the hole" by overriding
        the :math:`G` hook instead of extending ``assert_realizable``: that
        would re-assert the very conflation the 2026-08-01 ruling forbids
        (a wall's constitutive return posing as a symmetry of the domain).
        """
        from orpheus.geometry.boundary import SpecularReturn

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.array([1, 0, 3, 2, 5, 4, 7, 6]))
        law = AlbedoBoundary(0.5, SpecularReturn(axis="x"))
        # The G hook is silent — G is the identity deck element, which fixes
        # no geometry.
        law.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]
        # …while the aggregate still reds, via the R-tier extension.
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
            law.assert_realizable(mutant)  # type: ignore[arg-type]

    def test_a_diffuse_closure_does_NOT_fire_the_pairing_checks(self) -> None:
        """SCOPE control: only a SPECULAR closure stands on the mirror pairing.

        An ``IsotropicReturn`` closure realizes through the Lambertian average
        and never touches ``ordinate_permutation``, so a poisoned pairing must
        leave it untouched. Without this leg, an ``assert_realizable`` that
        fired the pairing checks unconditionally would pass every row above
        while refusing perfectly good diffuse laws.
        """
        from orpheus.geometry.boundary import IsotropicReturn

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.roll(np.arange(quad.n_ordinates), 1))
        AlbedoBoundary(
            0.5, IsotropicReturn(axis="x", outward_sign=+1),
        ).assert_realizable(mutant)  # type: ignore[arg-type]
        # …and the bare law likewise (no closure, no pairing).
        AlbedoBoundary(0.5).assert_realizable(mutant)  # type: ignore[arg-type]



# ─────────────────────────────────────────────────────────────────────
# B3.4c — the factor tier's adjointability claim must match what the
#         realizer actually builds
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestFactorAdjointabilityMatchesTheRealizedOperator:
    r"""``G.is_adjointable and R.is_adjointable == realized.is_adjointable``.

    Two independent encodings of one claim, compared — the ERR-041 discipline
    applied to adjointability. The factor tier DECLARES whether the map it
    names exposes an honest transpose; the realized operator ANSWERS from its
    own composition tree. Nothing compared them before **B3.4c**, and they had
    silently drifted apart: ``SpatialWrap`` (the wrap deck's axis-lettered
    predecessor, retired onto
    :class:`~orpheus.geometry.boundary.PairedDeck`) declared ``False`` while the
    operator realizing it answered ``True``, so a consumer got opposite
    answers depending on which one it asked. (The declaration was reporting an
    implementation gap — #183, periodic's unbuilt partner channel — in a slot
    whose contract is a property of the MAP. B3.4c built the channel, the
    declaration became true, and this gate is what keeps the two tiers from
    parting again.)

    The conjunction is the right shape because the realized operator is
    :math:`R \circ G`: a composition is Euclidean-adjointable exactly when both
    factors are. It bites in both directions — a factor that over-claims
    (``True`` with no realized transpose) and one that under-claims (the
    B3.4c-era ``SpatialWrap``) are both red.
    """

    #: ``(law_id, law)`` over every law with a realizable SN arm. The albedo
    #: rows carry both closures because since B3.4b an albedo face answers by
    #: its CLOSURE, not its class, so a single row would pin only one arm.
    #: (That comment continued "a specular closure is adjointable and a diffuse
    #: one is not" until 2026-08-04; G6.3 step 3 factored the diffuse
    #: realization and both are now adjointable. Both rows STAY — the
    #: closure-dispatch mechanism is what they pin, not the answers.)
    def _laws(self):
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            IsotropicReturn,
            SpecularReturn,
        )

        return [
            ("vacuum", VacuumInflow()),
            ("reflective", ReflectiveBoundary(axis="x")),
            ("white", WhiteBoundary(axis="x", outward_sign=-1, albedo=1.0)),
            ("albedo_specular",
             AlbedoBoundary(0.5, SpecularReturn(axis="x"))),
            ("albedo_isotropic",
             AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=-1))),
            ("periodic", PeriodicBoundary(axis="x")),
            # ⭐ P3 — prescribed inflow JOINED this set, and the set is now
            # every shipped law with no documented exception. It stood outside
            # until P3 because its realized affine operator declined a
            # transpose while both its factors declare adjointable; realizing
            # the LINEAR factor (the zero map, 0ᵀ = 0) closed that gap. See
            # ``test_prescribed_inflow_joined_the_conjunction``.
            ("prescribed_inflow",
             PrescribedInflow(source=ConstantInflowSource(value=2.5))),
        ]



    def test_every_shipped_kernel_is_now_adjointable(self) -> None:
        """Records the state change that emptied the ``False`` column.

        Kept as its own row so the one-valued set is an ASSERTED fact rather
        than an accident nobody noticed. If a future law arrives that genuinely
        cannot expose a transpose, this reddens and the gate above regains a
        natural counterexample.
        """
        answers = {
            law_id: (
                law.geometry_map.is_adjointable
                and law.response_kernel.is_adjointable
            )
            for law_id, law in self._laws()
        }
        assert set(answers.values()) == {True}, (
            f"a shipped law now declares non-adjointable: "
            f"{ {k: v for k, v in answers.items() if not v} }"
        )

