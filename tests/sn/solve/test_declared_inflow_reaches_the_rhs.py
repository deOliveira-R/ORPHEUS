r"""A DECLARED prescribed inflow reaches the solve's right-hand side.

Phase **P2′** of `.claude/plans/archive/affine_boundary_source_channel.md`. The affine
boundary law is :math:`\gamma_-\psi = L\,\gamma_+\psi + q`, and the composite
source :math:`q = q_{\rm bulk} \oplus q_\partial` has always been the channel
:math:`q_\partial` belongs in. What was missing is the step *before* it: a
:class:`~orpheus.geometry.boundary.PrescribedInflow` **declared as a boundary
law** never reached that channel. It realized into an affine operator carrying
no :attr:`BlockRole.BOUNDARY` stamp, which nothing consumed — so the
declaration was **silently inert**.

⭐ **Why these gates exist in this shape** (user ruling, 2026-08-05):

    *Tests must route through the machinery that a user would exercise without
    bypassing code functionality. Or else it's not testing the path the users
    go through.*

The pre-existing non-vacuum coverage supplies :math:`q_\partial` by hand
(``AngularBoundarySourceSink.prescribed_inflow`` straight into a composite
source), deliberately bypassing the law tier — see
``docs/theory/verification/sn.rst``'s §4.6 MMS, which says so explicitly. That
verifies the CHANNEL and is silent about every step a user actually takes. These
gates start from the **declaration** and assert it arrives.

⚠ **The keystone is** :meth:`TestTheDeclarationIsNotInert.test_a_declared_inflow_changes_the_rhs`.
Every other row here would pass against a ``from_mesh_laws`` that read the laws
correctly while nothing called it — which was exactly the pre-P2′ state, one
level down. That row is the one that fails if the wiring is removed.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import (
    ConstantInflowSource,
    NoSource,
    PrescribedInflow,
    ReflectiveBoundary,
    VacuumInflow,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import _build_fixed_source_rhs
from orpheus.transport.source_sinks import AngularBoundarySourceSink
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l1

_VALUE = 2.5


def _slab(ng: int = 2, n_ord: int = 8, nx: int = 4) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _declare(sn: SNMesh, face: str, law) -> SNMesh:
    """Install a law the way the mesh's own resolve body does."""
    sn.bc[face] = sn.realize_boundary_law(law, face)
    return sn


def _bulk(sn: SNMesh) -> np.ndarray:
    return np.zeros((sn.quad.N, sn.ng, *sn.spatial_shape))


class TestOnlyPrescribedInflowContributes:
    r""":math:`q = 0` for every law whose content is the linear factor."""

    @pytest.mark.parametrize(
        "law",
        [VacuumInflow(), ReflectiveBoundary(axis="x", albedo=1.0),
         PrescribedInflow(source=NoSource())],
        ids=["vacuum", "reflective", "prescribed(NoSource)"],
    )
    def test_a_sourceless_law_contributes_nothing(self, law) -> None:
        """Vacuum, reflective, and a zero-source prescribed all give q = 0.

        The third is the discriminating one: it IS the prescribed family, so a
        family-only test would wrongly expect a contribution. What decides is
        the source's VALUE.
        """
        sn = _declare(_slab(), "xmin", law)
        assert AngularBoundarySourceSink.from_mesh_laws(sn).linf == 0.0

    def test_a_declared_inflow_lands_on_that_faces_inflow_rows_only(self) -> None:
        r""":math:`q \in \Gamma_-(f)` for the declaring face, zero elsewhere."""
        sn = _declare(
            _slab(), "xmin",
            PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
        )
        q = AngularBoundarySourceSink.from_mesh_laws(sn)
        inflow = sn.angular_trace.inflow_indices_for_face("xmin")
        view = np.asarray(q.face_view("xmin"))

        np.testing.assert_array_equal(view[inflow], _VALUE)
        np.testing.assert_array_equal(np.delete(view, inflow, axis=0), 0.0)
        # The OTHER face still declares vacuum and must stay untouched — the
        # per-face mapping is not broadcast.
        np.testing.assert_array_equal(np.asarray(q.face_view("xmax")), 0.0)

    def test_two_faces_can_declare_independently(self) -> None:
        """Different faces, different sources — neither overwrites the other."""
        sn = _slab()
        _declare(sn, "xmin", PrescribedInflow(ConstantInflowSource(value=1.0)))
        _declare(sn, "xmax", PrescribedInflow(ConstantInflowSource(value=3.0)))
        q = AngularBoundarySourceSink.from_mesh_laws(sn)
        for face, value in (("xmin", 1.0), ("xmax", 3.0)):
            rows = sn.angular_trace.inflow_indices_for_face(face)
            np.testing.assert_array_equal(
                np.asarray(q.face_view(face))[rows], value
            )


class TestTheDeclarationIsNotInert:
    r"""⭐ The keystone — the declaration reaches the RHS the solvers consume."""

    def test_a_declared_inflow_changes_the_rhs(self) -> None:
        r"""⭐ The gate that fails if P2′'s wiring is removed.

        ``_build_fixed_source_rhs`` is the single construction point for the
        composite ``q = q_bulk ⊕ q_∂`` that BOTH inner paths (source iteration
        and Krylov) consume. Before P2′ its bulk-array arm hard-coded
        ``zeros_on`` — so a user could declare an inflow and get vacuum, with
        no error anywhere.

        Both legs are necessary: the vacuum leg pins that the arm is not simply
        always non-zero, and the declared leg pins the arrival.
        """
        vacuum = _slab()
        rhs_vacuum = _build_fixed_source_rhs(_bulk(vacuum), vacuum)
        assert rhs_vacuum.boundary.linf == 0.0

        declared = _declare(
            _slab(), "xmin",
            PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
        )
        rhs_declared = _build_fixed_source_rhs(_bulk(declared), declared)
        assert rhs_declared.boundary.linf == pytest.approx(_VALUE), (
            "a DECLARED prescribed inflow did not reach the composite RHS — "
            "the declaration is inert, which is the whole defect P2′ closes"
        )

    def test_the_bulk_leaf_is_untouched_by_the_boundary_wiring(self) -> None:
        """The two leaves are independent — q_∂ arriving must not perturb q_bulk.

        Compared against the SAME bulk array through a vacuum mesh, so the only
        difference between the two runs is the declaration.
        """
        rng = np.random.default_rng(20260805)
        vacuum = _slab()
        bulk = rng.standard_normal((vacuum.quad.N, vacuum.ng, *vacuum.spatial_shape))
        declared = _declare(
            _slab(), "xmin",
            PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
        )
        np.testing.assert_array_equal(
            np.asarray(_build_fixed_source_rhs(bulk, vacuum).interior.values),
            np.asarray(_build_fixed_source_rhs(bulk, declared).interior.values),
        )


class TestTheSourceCannotBeSpecifiedTwice:
    r"""Two answers to one question is refused, not silently resolved."""

    def test_declaring_AND_supplying_a_composite_boundary_raises(self) -> None:
        """Adding would double-count; overriding would make one input a no-op.

        This is reachable: a caller may legitimately pass a composite source
        (the pre-P2′ way to get an inflow) against a mesh that now also
        declares one.
        """
        sn = _declare(
            _slab(), "xmin",
            PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
        )
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        inflow_rows = sn.angular_trace.inflow_indices_for_face("xmin")
        slot = np.zeros(
            AngularBoundarySourceSink.zeros(sn.angular_trace).face_view("xmin").shape
        )
        slot[inflow_rows] = 1.0
        composite = TimedFullField(
            interior=AngularSourceSink.from_isotropic(
                np.zeros((sn.ng, *sn.spatial_shape)), sn
            ),
            boundary=AngularBoundarySourceSink.prescribed_inflow(
                sn, {"xmin": slot}
            ),
            _history=(),
            history_depth=2,
        )
        with pytest.raises(ValueError, match="specified"):
            _build_fixed_source_rhs(composite, sn)

    def test_a_composite_boundary_is_still_accepted_without_a_declaration(
        self,
    ) -> None:
        """Positive leg: the pre-existing composite path is NOT broken.

        Without it the row above would pass against an implementation that
        rejected every composite boundary source — which would break the §4.6
        MMS and every other direct-supply consumer.
        """
        sn = _slab()  # all vacuum — nothing declared
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        inflow_rows = sn.angular_trace.inflow_indices_for_face("xmin")
        slot = np.zeros(
            AngularBoundarySourceSink.zeros(sn.angular_trace).face_view("xmin").shape
        )
        slot[inflow_rows] = 1.0
        composite = TimedFullField(
            interior=AngularSourceSink.from_isotropic(
                np.zeros((sn.ng, *sn.spatial_shape)), sn
            ),
            boundary=AngularBoundarySourceSink.prescribed_inflow(
                sn, {"xmin": slot}
            ),
            _history=(),
            history_depth=2,
        )
        assert _build_fixed_source_rhs(composite, sn).boundary.linf == 1.0


# ─────────────────────────────────────────────────────────────────────
# P3 — ⭐ SINGLE DELIVERY, asserted on the converged answer
# ─────────────────────────────────────────────────────────────────────


def _het_slab(n_ord: int = 8) -> SNMesh:
    r"""Heterogeneous, 2G, asymmetric ``SigS``, scattering-active.

    NOT ``placeholder_materials``: ``SigS ≡ 0`` there, so a delivered inflow
    never feeds a scattering source and the SI solve collapses to a single
    sweep. The delivery-COUNT claim survives that degeneracy, but a claim read
    off the converged flux would be config-blind — and ``c ≈ 0.90–0.96`` here
    means the inflow actually propagates through the scattering source.
    """
    from orpheus.derivations.common.xs_library import get_mixture

    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=1.0),
                 Region(mat_id=1, outer_thickness_cm=2.0)),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=6), RegionMesh(n_cells=6)),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=n_ord),
        {0: get_mixture("A", "2g"), 1: get_mixture("D", "2g")},
    )


def _composite(sn: SNMesh, *, boundary_value: float = 0.0,
               bulk_level: float = 1.0):
    r"""ONE bulk spelling for every leg of every row below — deliberately.

    ⚠ ``_build_fixed_source_rhs``'s two arms take DIFFERENT bulk sources: the
    bulk-array arm a per-ordinate ``(N, ng, *spatial)`` array, the composite arm
    an :class:`AngularSourceSink`. Feeding one to each leg of a two-channel
    comparison produces a difference in the BULK that reads as a difference in
    the boundary channel (measured: ``φ[0] = 3.083`` vs ``2.480``, which cost a
    probe). Always routing through this helper makes that trap unspellable.
    """
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.timed_full_field import TimedFullField

    slot = np.zeros(
        AngularBoundarySourceSink.zeros(sn.angular_trace).face_view("xmin").shape
    )
    slot[sn.angular_trace.inflow_indices_for_face("xmin")] = boundary_value
    return TimedFullField(
        interior=AngularSourceSink.from_isotropic(
            np.full((sn.ng, *sn.spatial_shape), bulk_level), sn,
        ),
        boundary=AngularBoundarySourceSink.prescribed_inflow(sn, {"xmin": slot}),
        _history=(),
        history_depth=2,
    )


def _solve(sn: SNMesh, inner: str, **kw):
    from orpheus.sn.solver import (
        SNSolver, _solve_fixed_source_krylov, _solve_fixed_source_si,
    )

    q = _build_fixed_source_rhs(_composite(sn, **kw), sn)
    solver = SNSolver(
        sn, inner_solver=inner, scattering_order=0,
        max_inner=2000, inner_tol=1e-13,
    )
    driver = (
        _solve_fixed_source_si if inner == "source_iteration"
        else _solve_fixed_source_krylov
    )
    return driver(solver, sn, q, 0.0, 2000, 1e-13)


def _gamma_minus(sol, sn: SNMesh, face: str) -> np.ndarray:
    """The converged INFLOW trace of ``face`` — the boundary law's own output."""
    rows = sn.angular_trace.inflow_indices_for_face(face)
    return np.asarray(sol.angular_flux.boundary.face_view(face))[rows]


@pytest.mark.catches("ERR-075")
@pytest.mark.verifies("bc-single-delivery")
@pytest.mark.parametrize("inner", ["source_iteration", "krylov"])
def test_the_declared_boundary_law_holds_on_the_answer(inner: str) -> None:
    r"""⭐⭐ SINGLE DELIVERY: :math:`\gamma_-\psi|_f = q_f` exactly, at convergence.

    **The row the whole campaign is for.** The affine boundary law is
    ``γ₋ψ = L γ₊ψ + q`` and prescribed inflow is the ``L = 0`` case, so the
    converged inflow trace IS the declared source. No reference solver, no
    discretization dependence, no tolerance: the boundary condition is a
    DEFINITION, and this reads it off the answer.

    `[M]` at ``ef4c3537`` (pre-P3), on this fixture, printed at 12 dp:

    ===========================================  ============
    configuration                                ``γ₋(xmin)``
    ===========================================  ============
    double delivery (``B`` affine AND ``q_∂``)   ``5.0``
    single delivery (either channel alone)       ``2.5``
    channel lost                                 ``0.0``
    ===========================================  ============

    Three configurations, three exactly-distinguishable readings, one
    assertion. ``assert_array_equal`` rather than ``allclose`` is MEASURED, not
    assumed: the double reading is ``2.5 + 2.5`` (exact in IEEE) and the single
    reading is a copy. Should a future moment-resolved trace make it inexact,
    descend to ``assert_array_almost_equal_nulp`` and say why — do NOT reach for
    ``rtol``, which would blur the 1× / 2× distinction this row exists to make.

    ⭐ **Why ``B(0) == 0`` at the leaf cannot replace this row.** That gate sees
    a non-linear ``B``, and the realistic regressions here are LINEAR: setting
    ``L := IdentityOperator`` for prescribed inflow is perfectly linear and
    perfectly invisible to it, while this row reddens immediately. Conversely a
    superposition/slope gate (``φ`` affine in ``q``) is a provable non-catcher —
    a doubled delivery is ``q → 2q``, which is *still exactly affine in* ``q``,
    so the slope test passes with the wrong slope (``vv`` Mode 12).

    ⭐ **The ``krylov`` parameter is load-bearing.** Pre-P3 it did not merely
    give a wrong number, it RAISED: an affine ``A(x) = A_lin(x) − c`` breaks
    GMRES's Arnoldi relation ``A V_k = V_{k+1} H_k``, so the residual scipy
    tracks is meaningless and ``_certify_within_group_exit`` catches the
    divergence (``ConvergenceCertificateError``, ``‖Aψ − q‖/‖q‖ = 1.718``).
    Declared prescribed inflow × Krylov was UNUSABLE — before P2′ as well as
    after it — and only #189 (the law is not a registered ``BC`` kind) kept that
    out of every production driver.
    """
    sn = _declare(
        _het_slab(), "xmin",
        PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
    )
    sol = _solve(sn, inner)

    got = _gamma_minus(sol, sn, "xmin")

    # Leg 1 — the DELIVERY COUNT, which is what this row is for. ``atol`` is
    # 1e-9 against readings of 0.0 / 2.5 / 5.0: eleven orders of margin, so the
    # three cases stay exactly distinguishable while the message stays readable.
    np.testing.assert_allclose(
        got, _VALUE, rtol=0.0, atol=1e-9,
        err_msg=(
            f"[{inner}] γ₋(xmin) = {got.min():.12f}..{got.max():.12f}, "
            f"expected {_VALUE}. {_VALUE * 2} means the inflow is delivered "
            f"TWICE (an affine operator is back in the B block AND q_∂ carries "
            f"it); 0.0 means the channel is lost."
        ),
    )

    # Leg 2 — the EXACTNESS claim, which differs by path and is MEASURED.
    if inner == "source_iteration":
        # SI writes the source into the inflow slot and sweeps from it, so the
        # converged trace is a COPY of q — bit-exact, and asserted as such.
        np.testing.assert_array_equal(got, _VALUE)
    else:
        # ⚠ Krylov SOLVES the trace rows instead of copying a seed into them, so
        # they carry the ITERATION RESIDUAL: `[M]` 2.500000000000008, i.e. 18 ULP
        # at 2.5 (7.99e-15 absolute) on this fixture.
        #
        # ⛔ The deviation scales with ``inner_tol``, NOT with machine epsilon.
        # `[M]` across fixtures at tol=1e-13: 1, 7, 18, 23 ULP (GL-8 and GL-16,
        # vacuum and reflective partner faces) — but at ``inner_tol = 1e-10`` the
        # same fixture reads 1.2e-11, i.e. **27 580 ULP**. So ``nulp=64`` is a
        # budget for THIS module's ``inner_tol = 1e-13`` (pinned in ``_solve``),
        # not a universal floating-point claim. Loosening the solver tolerance
        # here without raising the budget is a legitimate red, and the right
        # response is to re-derive the budget from the new tolerance rather than
        # to inflate this number.
        #
        # Kept in ULP rather than as an ``rtol`` even so: a ULP count reads as a
        # floating-point claim that someone must justify, whereas an ``rtol``
        # invites being scaled off ``inner_tol`` until it spans the 1×/2× gap
        # this row exists to separate.
        np.testing.assert_array_almost_equal_nulp(
            got, np.full(got.shape, _VALUE), nulp=64,
        )
    # The other face declares VACUUM, so its inflow trace is zero. Without this
    # leg the row is green against a q broadcast to every face.
    np.testing.assert_array_equal(
        _gamma_minus(sol, sn, "xmax"), 0.0,
        err_msg=(
            f"[{inner}] γ₋(xmax) is non-zero on a face declaring vacuum — the "
            f"declared source is reaching faces that did not declare it."
        ),
    )


@pytest.mark.parametrize("inner", ["source_iteration", "krylov"])
def test_an_undeclared_face_has_a_zero_inflow_trace(inner: str) -> None:
    r"""PC⁻, the negative control for the row above.

    An all-vacuum mesh under the identical bulk source must read
    :math:`\gamma_-\psi = 0` on both faces. Without this, the keystone would
    pass against an implementation that wrote ``_VALUE`` onto every inflow trace
    unconditionally — and its own second leg could not catch that, because it
    reads the same mesh that declared the source.
    """
    sn = _het_slab()  # nothing declared
    sol = _solve(sn, inner)
    for face in ("xmin", "xmax"):
        np.testing.assert_array_equal(
            _gamma_minus(sol, sn, face), 0.0,
            err_msg=f"[{inner}] γ₋({face}) is non-zero with no law declared",
        )


def test_the_two_user_paths_reach_the_same_fixed_point() -> None:
    r"""Gate B: DECLARING the law and SUPPLYING ``q_∂`` by hand agree.

    The two routes a user can take to the same physics:

    * declare ``PrescribedInflow(ConstantInflowSource(2.5))`` on the mesh and
      let :meth:`from_mesh_laws` assemble ``q_∂`` (the law tier);
    * install vacuum and hand :meth:`prescribed_inflow` the same values
      directly (the channel tier — what the §4.6 MMS does).

    ⭐ **BIT-IDENTICAL, and that is the post-carve claim.** `[M]`
    ``|φ_declared − φ_supplied|_inf = 0.000000e+00``, ``array_equal = True``, on
    both inner solvers at ``‖φ‖ = 6.838442``. After P3 the two routes are not two
    computations that agree — they are ONE float program: both assemble the same
    ``q_∂`` into the same composite RHS, and the realized ``B`` contributes zero
    on both. Nothing is left to differ.

    ⛔ **Do NOT weaken this to a tolerance.** Pre-P3 the same comparison read
    ``1.998e-13`` with ``array_equal = False`` — because *then* the declared route
    delivered through ``B`` while the supplied route delivered through ``q_ext``,
    two genuinely different computations reaching one fixed point. That number
    belongs to the DEFECT, not to the fix, and a gate built on it would be blind
    to the very thing this row exists for: a ``2.9e-14`` relative divergence
    sails through any ``rtol`` scaled off ``inner_tol``. If this row ever needs a
    tolerance again, the two channels have stopped being one program — which is
    the finding, not the noise.

    ⚠ Distinguish this from the keystone's Krylov leg, which is NOT exact. That
    leg compares ``γ₋ψ`` against the *literal* ``2.5``, and on Krylov the trace
    rows are SOLVED rather than copied, so they carry the iteration residual.
    Three different numbers, three different claims — conflating them is how a
    correct gate gets loosened.
    """
    declared = _declare(
        _het_slab(), "xmin",
        PrescribedInflow(source=ConstantInflowSource(value=_VALUE)),
    )
    supplied = _het_slab()  # vacuum everywhere; q_∂ comes in the composite

    phi_declared = np.asarray(
        _solve(declared, "source_iteration").scalar_flux.values
    )
    phi_supplied = np.asarray(
        _solve(supplied, "source_iteration",
               boundary_value=_VALUE).scalar_flux.values
    )
    scale = float(np.max(np.abs(phi_supplied)))
    np.testing.assert_array_equal(
        phi_declared, phi_supplied,
        err_msg=(
            "the law tier and the channel tier no longer produce the same float "
            "program. Post-P3 they must be BIT-identical: both assemble the same "
            "q_∂ into the same composite RHS and the realized B contributes zero "
            "to both. A non-zero difference here means a second delivery path "
            "has reappeared — do NOT relax this to a tolerance (see docstring)."
        ),
    )
    # Activation: the comparison is only meaningful if the inflow MOVED the
    # answer. Against the no-inflow control it must differ substantially.
    phi_vacuum = np.asarray(
        _solve(_het_slab(), "source_iteration").scalar_flux.values
    )
    assert np.max(np.abs(phi_declared - phi_vacuum)) > 1e-3 * max(scale, 1.0), (
        "the inflow barely changed the flux, so this row would pass even if "
        "both legs ignored q_∂ entirely"
    )
