r"""A DECLARED prescribed inflow reaches the solve's right-hand side.

Phase **P2′** of `.claude/plans/affine_boundary_source_channel.md`. The affine
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
            AngularBoundarySourceSink.zeros_on(sn).face_view("xmin").shape
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
            AngularBoundarySourceSink.zeros_on(sn).face_view("xmin").shape
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
