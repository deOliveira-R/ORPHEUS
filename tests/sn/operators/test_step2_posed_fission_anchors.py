r"""Consumers-campaign step 2 (ii) — PRE-carve anchors for the ONE posed ``F``.

Campaign: ``.claude/plans/cs4c_binding_design.md`` §27.2 (R-cc6 (ii)).
Verification plan: ``scratch/_consumers/test_architect_step2.md`` §3.
Sibling module (sub-steps (i) and (iii)):
``tests/sn/architecture/test_step2_terminal_object_anchors.py``.

The claim, in the domain's terms
================================

The eigenvalue Problem's terminal object is the pencil :math:`(A, F)` **on one
space**.  Today the forward k-problem and its adjoint build two different
:math:`F` objects on two different spaces, which is why the adjoint
*re-derives* fission instead of *viewing* it: ``F_adjoint`` cannot be ``F.H``
when ``F`` does not exist as one object.

``[M]`` 2026-09-12 at ``b0fd3e7e`` (``scratch/_consumers/probes2/p2_counts2.py``):

==============  =========================================  ================================
route            object                                     space
==============  =========================================  ================================
forward          ``IsotropicFission`` (``solver.py:1544``)   ``bulk_space`` ``(ng, *spatial)``
adjoint SEEDLESS ``FissionOperator``  (``solver.py:2772``)   ``FullFieldSpace``
adjoint CARRYING ``OperatorProduct``  (``solver.py:2821``)   ``CoupledSpace``
==============  =========================================  ================================

⛔⛔ The brief's *"adopt the adjoint path's posed F spelling forward"* therefore
names THREE spellings, not one: ``F_posed = stack @ restrict_bulk`` is the
CARRYING spelling only, and a gate written as *"the forward F is the adjoint's
F_posed"* is a false red on every slab and every 2-D Cartesian row.

⭐⭐ And the bit-identity is reachable exactly ONE way
====================================================

``[M]`` ``probes2/p3_F_equivalence.py`` / ``p4_spaces_and_seeds.py``, 200
seeds each:

* ``FissionOperator(...).isotropic_energy.apply(φ)`` vs today's forward
  ``IsotropicFission.from_material_xs(...).apply(φ)`` — **``array_equal``
  200/200, ``max|Δ| = 0.000000e+00``**;
* ``FissionOperator(...).apply(ψ)`` on the composite vs today's
  ``AngularSourceSink.from_isotropic(IsotropicFission.apply(∫ψ dΩ))`` lift —
  **``array_equal`` 0/200**, ``max|Δ| = 2.775558e-17``, ``max rel =
  2.214908e-16`` (≤ 1 nulp, draw-STABLE).

⟹ *"the hub owns ONE ``F = FissionOperator``; the forward k-outer reads its
DERIVED ``isotropic_energy``"* is bit-identical; *"the forward applies ``F`` on
the composite"* costs 1 nulp on every eigen solve.  Nine of the fourteen DD
regression cases already drift 1–11 ULP (the DriftWarning delta table in the
plan's §5.2), so that second route's cost is not absorbable as noise — it
would move the drift SET.  Plan open ruling **O-1**.

Marks
=====

``foundation`` — software/architecture invariants of the fission binding; no
theory ``:label:``, hence no ``verifies(...)``.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.geometry.coord import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver, _adjoint_posing_parts, _as_sn_mesh
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_transfer import IsotropicFission
from orpheus.transport.source_sinks import AngularSourceSink

pytestmark = pytest.mark.foundation


def _slab_hub() -> SNMesh:
    """Fuel | moderator 2G slab, GL-8, reflective | vacuum — SEEDLESS.

    The seedless arm is the one the brief's ``stack @ restrict_bulk``
    description does NOT cover, so it is the arm every row here runs on
    unless it says otherwise.
    """
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 9),
        mat_ids=np.array([0, 0, 0, 0, 1, 1, 1, 1]),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return _as_sn_mesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=8), materials,
        scattering_order=0,
    )


def _sphere_hub() -> SNMesh:
    """The CARRYING arm — a sphere, where the adjoint's ``F`` is the lift."""
    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.01, 2.0, 9),
        mat_ids=np.array([0, 0, 0, 0, 1, 1, 1, 1]),
        coord=CoordSystem.SPHERICAL,
    )
    return _as_sn_mesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=8), materials,
        scattering_order=0,
    )


def _random_composite(hub: SNMesh, seed: int) -> "tuple[FullField, AngularFlux]":
    """The composite AND its interior, narrowed.

    ``FullField.interior`` is typed by the ``BulkField`` protocol, which
    carries no ``integrate_angular``; the angular reduction is a fact about
    :class:`AngularFlux` specifically.  Returning both is the narrowing
    (``coding-elegance`` #19 — the principled spelling, not a cast).
    """
    rng = np.random.default_rng(seed)
    interior = AngularFlux(
        values=rng.random(hub.angular_trial_space.shape),
        space=hub.angular_trial_space,
    )
    return FullField(
        interior=interior,
        boundary=AngularBoundaryFlux.zeros(hub.angular_trace),
    ), interior


# ═══════════════════════════════════════════════════════════════════════
# RECORD — today's two (three) F spellings.  DELETE at sub-step (ii).
# ═══════════════════════════════════════════════════════════════════════


class TestRecordTwoFissionObjects:
    """The forward and the adjoint do not share an ``F``."""

    def test_record_the_forward_and_adjoint_F_are_different_types_on_different_spaces(
        self,
    ) -> None:
        """RECORD — ⛔ DELETE at sub-step (ii); do not repair.

        This row is what makes the API change LOUD.  Its ruled successor is
        :meth:`TestRuledOneFissionPerProblem.test_ruled_one_F_mint_per_problem`,
        and a ``strict`` xfail alone would stay silent if the carve landed one
        ``F`` with the wrong semantics (``vv`` Mode-8, fourth class).
        """
        hub = _slab_hub()
        mat_xs = hub.material_xs_field()
        forward = IsotropicFission.from_material_xs(
            mat_xs, space=mat_xs.mesh.bulk_space,
        )
        _implicit, _gain, adjoint_F, _template = _adjoint_posing_parts(hub)

        assert type(forward) is IsotropicFission
        assert type(adjoint_F) is FissionOperator
        assert forward.domain is hub.bulk_space
        assert adjoint_F.domain is hub.full_field_space
        assert forward.domain != adjoint_F.domain, (
            "the forward and adjoint F now share a domain — if this reds, "
            "sub-step (ii) has landed and this RECORD row must be DELETED."
        )

    def test_record_the_carrying_adjoint_is_a_THIRD_spelling(self) -> None:
        r"""RECORD — ``stack @ restrict_bulk`` is the CARRYING spelling only.

        ``[M]`` on a sphere the adjoint's ``F`` is an ``OperatorProduct`` on
        the ``CoupledSpace``; on the slab above it is a plain
        ``FissionOperator`` on the ``FullFieldSpace``.  Any step-2 gate that
        says "the adjoint's ``F_posed``" without naming the arm is a false
        red on 2 of 3 charts.
        """
        from orpheus.numerics.operator import OperatorProduct

        hub = _sphere_hub()
        _implicit, _gain, adjoint_F, _template = _adjoint_posing_parts(hub)
        assert isinstance(adjoint_F, OperatorProduct), (
            f"the carrying adjoint's F is {type(adjoint_F).__name__}, not the "
            f"rectangular lift the campaign's §27.1 table describes."
        )
        assert adjoint_F.domain is not hub.full_field_space

    def test_record_two_F_mints_per_problem(self) -> None:
        """RECORD — one hub, a forward solver and an adjoint posing: TWO mints.

        ``[M]`` ``probes2`` — ``IsotropicFission.from_material_xs`` fires once
        (``solver.py:1544``, ``SNSolver.__init__``) and
        ``FissionOperator.from_solver_data`` once (``solver.py:2772``,
        ``_adjoint_posing_parts``).  ⛔ DELETE at sub-step (ii).
        """
        counts = _count_fission_mints(_slab_hub())
        assert counts == {"energy": 1, "angular": 1}, (
            f"expected one mint of each binding per Problem, saw {counts}."
        )


def _count_fission_mints(hub: SNMesh) -> dict[str, int]:
    """Mint census over ONE hub: build the forward solver, then pose the adjoint.

    Patched on the CLASS (both factories are ``classmethod``s), so every
    module binding resolves through the spy — ``lessons`` L46e's rebinding
    rule with the receiver that admits only one binding site.
    """
    counter: collections.Counter[str] = collections.Counter()
    originals = {}
    for cls, name, key in (
        (IsotropicFission, "from_material_xs", "energy"),
        (FissionOperator, "from_solver_data", "angular"),
    ):
        originals[(cls, name)] = getattr(cls, name).__func__

        def make(original=originals[(cls, name)], key=key):  # noqa: B008
            def spy(owner, *args, **kwargs):  # noqa: ANN001, ANN002, ANN003
                counter[key] += 1
                return original(owner, *args, **kwargs)

            return classmethod(spy)

        setattr(cls, name, make())
    try:
        SNSolver(hub)
        _adjoint_posing_parts(hub)
    finally:
        for (cls, name), original in originals.items():
            setattr(cls, name, classmethod(original))
    if not counter:
        raise RuntimeError(
            "the fission-mint census counted nothing — the instrument is "
            "dead and its zero carries no information (vv anti-#17).",
        )
    return dict(counter)


# ═══════════════════════════════════════════════════════════════════════
# THEOREM — the walls sub-step (ii) must not break
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheDerivedEnergyBindingIsBitIdentical:
    r"""⭐⭐ The bit-identity wall for the RECOMMENDED route.

    ``FissionOperator`` is an ``AngularLift[IsotropicFission]``
    (``fission.py:214``) and DERIVES its energy binding from the ONE datum it
    carries.  ``[M]`` that derived binding's ``apply`` is ``array_equal`` to
    today's forward mint on **200/200 seeds** (``max|Δ| = 0.000000e+00``), and
    the derived object is cached (``F.isotropic_energy is F.isotropic_energy``
    → ``True``).

    ⟹ *the hub owns ONE ``F``; the forward k-outer reads ``F.isotropic_energy``*
    is a re-homing with NO arithmetic in it.  This row is the wall that says
    so, and it is green before AND after the carve.
    """

    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 5, 8, 13, 21])
    def test_law_derived_binding_equals_todays_forward_mint(
        self, seed: int,
    ) -> None:
        hub = _slab_hub()
        mat_xs = hub.material_xs_field()
        forward = IsotropicFission.from_material_xs(
            mat_xs, space=mat_xs.mesh.bulk_space,
        )
        angular = FissionOperator.from_solver_data(
            mat_xs=mat_xs, space=hub.full_field_space,
        )
        phi = np.random.default_rng(seed).random(hub.bulk_space.shape)
        np.testing.assert_array_equal(
            np.asarray(forward.apply(phi)),
            np.asarray(angular.isotropic_energy.apply(phi)),
            err_msg=(
                "the derived energy binding is no longer bit-identical to the "
                "forward mint — sub-step (ii)'s bit-identity claim rests on "
                "exactly this (200/200 seeds at b0fd3e7e); a tolerance here "
                "would admit the very drift the claim denies."
            ),
        )

    def test_law_the_derived_binding_carries_the_hubs_bulk_space(self) -> None:
        r"""The ends agree by ``==``, and ⚠ NOT by ``is``.

        ``[M]`` ``F.isotropic_energy.domain == hub.bulk_space`` is **True**
        while ``is`` is **False** (the lift derives a scalar sub-space off its
        codomain rather than reading the hub's instance).  ``lessons`` L77a /
        L72b: an ends gate spelled ``is`` here is a false red, so this row
        asserts ``==`` and records the identity fact instead of asserting it —
        if a later step interns the space, nothing here reds.
        """
        hub = _slab_hub()
        angular = FissionOperator.from_solver_data(
            mat_xs=hub.material_xs_field(), space=hub.full_field_space,
        )
        assert angular.isotropic_energy.domain == hub.bulk_space
        assert angular.isotropic_energy.codomain == hub.bulk_space
        assert angular.domain is hub.full_field_space
        assert angular.codomain is hub.full_field_space


class TestLawTheCompositeRouteIsOneNulpAwayNotBitIdentical:
    r"""⛔ The route the brief's words describe, PRICED.

    ``[M]`` 200 seeds: the composite route
    (``FissionOperator.apply`` on a ``FullField``) and today's forward route
    (``IsotropicFission`` + ``AngularSourceSink.from_isotropic``) are
    ``array_equal`` on **0/200**, with ``max|Δ| = 2.775558e-17`` and
    ``max rel = 2.214908e-16`` — draw-STABLE at ≤ 1 nulp.  The gap is the
    ``/W`` ordering: the lift normalises before the broadcast, the angular
    binding after (``vv`` §bit-identity criterion 3 — an IEEE re-association,
    not an algebraic change).

    This row is the PRICE tag on open ruling O-1.  It is a measured band, not
    a contract about which route production takes: it stays green either way,
    and its value is that the number is committed rather than living in a
    memo (``lessons`` L39 — a measured number in a comment is not a gate).
    """

    @pytest.mark.parametrize("seed", [1000, 1001, 1002, 1007])
    def test_law_the_two_routes_agree_to_one_nulp(self, seed: int) -> None:
        hub = _slab_hub()
        mat_xs = hub.material_xs_field()
        forward = IsotropicFission.from_material_xs(
            mat_xs, space=mat_xs.mesh.bulk_space,
        )
        angular = FissionOperator.from_solver_data(
            mat_xs=mat_xs, space=hub.full_field_space,
        )
        psi, interior = _random_composite(hub, seed)

        scalar = interior.integrate_angular()
        lifted = AngularSourceSink.from_isotropic(
            np.asarray(forward.apply(np.asarray(getattr(scalar, "values", scalar)))),
            hub,
        )
        composite = angular.apply(psi)

        np.testing.assert_array_almost_equal_nulp(
            np.asarray(lifted.values),
            np.asarray(composite.interior.values),
            nulp=1,
        )

    def test_law_fission_sources_no_trace(self) -> None:
        r"""ACTIVATION + THEOREM — the composite ``F`` emits ZERO on the trace.

        Fission is a bulk emission; the boundary block of ``F·ψ`` is exactly
        ``0.0`` (``[M]`` ``max|boundary| = 0.0``).  Shipped as its own row so
        the nulp band above cannot be read as covering the trace, and with a
        non-vacuity leg on the bulk so a fixture that produced nothing at all
        would red rather than pass.
        """
        hub = _slab_hub()
        angular = FissionOperator.from_solver_data(
            mat_xs=hub.material_xs_field(), space=hub.full_field_space,
        )
        out = angular.apply(_random_composite(hub, 42)[0])
        bulk = np.asarray(out.interior.values)
        assert float(np.max(np.abs(bulk))) > 0.0, (
            "non-vacuity: F emitted nothing on the bulk, so the trace-zero "
            "claim below is 0 == 0."
        )
        np.testing.assert_array_equal(
            np.asarray(out.boundary.values),
            np.zeros_like(np.asarray(out.boundary.values)),
            err_msg="the composite fission binding emitted a boundary source.",
        )


# ═══════════════════════════════════════════════════════════════════════
# XFAIL(strict) — the ruled post-carve behaviour.  RED today.
# ═══════════════════════════════════════════════════════════════════════


class TestRuledOneFissionPerProblem:
    """R-cc6 (ii) — ONE ``F``, minted once, read by both faces."""

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "R-cc6 (ii) NOT LANDED — the forward k-outer mints its own "
            "IsotropicFission (solver.py:1544) while the adjoint mints a "
            "FissionOperator (solver.py:2772): two objects, two spaces, so "
            "F_adjoint cannot be F.H. WHEN THIS XPASSES: the Problem owns one "
            "F — delete this marker AND "
            "TestRecordTwoFissionObjects.test_record_two_F_mints_per_problem."
        ),
    )
    def test_ruled_one_F_mint_per_problem(self) -> None:
        """Counted, not named.

        The plan does not rule what the hub's member is CALLED (open ruling
        O-4), so this row asserts the COUNT — exactly one fission mint across
        a forward-solver construction and an adjoint posing on ONE hub —
        rather than a guessed attribute name (``plan-authoring`` §1: a name in
        an assertion is a guess wearing a contract).
        """
        counts = _count_fission_mints(_slab_hub())
        total = sum(counts.values())
        assert total == 1, (
            f"one Problem minted {total} fission bindings ({counts}); the "
            f"pencil's F is one object or the adjoint is a re-derivation."
        )
