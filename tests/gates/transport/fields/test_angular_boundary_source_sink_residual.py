r"""Foundation tests for the B.3 boundary role leaves
:class:`~orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink` and
:class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`.

These complete the ``{Boundary}×{Source,Residual}`` cells of the field
role grid (siblings of
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`), giving
the :class:`~orpheus.transport.fields._bases.AngularBoundaryField` storage base
its 2nd and 3rd concrete instances. All three boundary leaves are empty
leaves — storage / validation / algebra / per-face access / factories
are inherited from ``AngularBoundaryField`` — so most of their machinery is
already pinned by ``test_boundary_flux.py``. This module adds:

* construction of the two NEW leaves (space-primary zeros /
  ``from_face_arrays``), and
* the **load-bearing cross-class invariant** unique to the boundary
  family: all three leaves share the SAME ``AngularTraceSpace`` (``angular_trace_space(mesh, quadrature, scheme)``),
  so the space gate would PASS — it is the **class-identity** gate that
  rejects ``AngularBoundarySourceSink ± AngularBoundaryFlux`` etc. This is the boundary
  analogue of the bulk "same units / same space ≠ same meaning"
  discipline.

``foundation`` tests — software invariants, not physics claims.

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3.
* ``orpheus/numerics/field.py`` — Layer-1 class-identity gate.
* ``tests/gates/transport/fields/test_boundary_flux.py`` — the AngularBoundaryFlux
  contract these leaves inherit.

Layer DAG Phase 4 (2026-11-01): the file builds material meshes and spells each space it reads from the raw
constructors of orpheus.transport.phase_spaces (the SN problem's and the diffusion mesh's cached mints are those
constructors on the problem's own ingredients, gated in
tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py and
tests/gates/diffusion/test_scalar_trace_is_the_raw_constructor.py).
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import angular_trace_space
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields._bases import AngularBoundaryField
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.residuals import AngularBoundaryResidual
from orpheus.transport.fields import AngularBoundarySourceSink

from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from tests.gates.transport._phase_space_cases import posing_quadrature, posing_scheme
from tests.gates.transport._phase_space_cases import slab_mesh


pytestmark = [pytest.mark.foundation]

# The two NEW boundary leaves (parametrized over for shared contract).
NEW_BOUNDARY_LEAVES = [AngularBoundarySourceSink, AngularBoundaryResidual]


# ── Fixtures ─────────────────────────────────────────────────────────


def _quad8_mesh(nx: int = 5, ng: int = 2) -> MaterialMesh:
    """GL(8) sibling of ``_slab_mesh`` — the boundary-trace CONTENT (and
    shape) differs, so the carrier mints an UNEQUAL trace space (F2)."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 5, ng: int = 2) -> MaterialMesh:
    mesh = Mesher(StructuredGeometry.sphere((0.0, 1.0), (0,), outer=0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def _cartesian_2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> MaterialMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
    )
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# Construction + inheritance (both new leaves share AngularBoundaryFlux's
# inherited machinery; spot-check it is wired on the new classes).
# ════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("Leaf", NEW_BOUNDARY_LEAVES)
class TestConstructionInherited:
    def test_inherits_field_and_boundary_field(self, Leaf) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = Leaf.zeros(trace)
        assert isinstance(bf, Field)
        assert isinstance(bf, AngularBoundaryField)
        assert isinstance(bf, Leaf)
        assert bf.space is trace

    # ``test_zeros_on_uses_mesh_trace`` RETIRED at CS4b S5: it pinned the
    # sugar factory's space SOURCE ("zeros_on reads angular_trace_space(mesh, quadrature, scheme)"),
    # and the space-primary spelling makes that identity true by
    # construction at the call site — the carrier-mint claims live on the
    # trace-space gates. The layout claim survives in
    # test_inherits_field_and_boundary_field above.

    def test_sphere_layout_only_xmax(self, Leaf) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = Leaf.zeros(trace)
        assert set(bf.layout.faces) == {"xmax"}

    def test_from_face_arrays_slab(self, Leaf) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        bf = Leaf.from_face_arrays_on(trace, {"xmin": np.full((N, 2), 1.0), "xmax": np.full((N, 2), 2.0)},
        )
        assert isinstance(bf, Leaf)
        np.testing.assert_array_equal(bf.face_view("xmin"), 1.0)
        np.testing.assert_array_equal(bf.face_view("xmax"), 2.0)

    def test_from_face_arrays_2d(self, Leaf) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        bf = Leaf.from_face_arrays_on(trace, {
            "xmin": np.full((N, 2, 2), 1.0),
            "xmax": np.full((N, 2, 2), 2.0),
            "ymin": np.full((N, 2, 3), 3.0),
            "ymax": np.full((N, 2, 3), 4.0),
        })
        np.testing.assert_array_equal(bf.face_view("ymin"), 3.0)

    def test_post_init_requires_trace_space(self, Leaf) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        from orpheus.numerics.space import FunctionSpace
        plain = FunctionSpace(name="sn_boundary_flat", shape=(trace.shape[0],))
        with pytest.raises(TypeError, match="AngularTraceSpace"):
            Leaf(values=np.zeros(trace.shape[0]), space=plain)


@pytest.mark.parametrize("Leaf", NEW_BOUNDARY_LEAVES)
class TestAlgebraClosedWithinClass:
    def test_add_sub_closed(self, Leaf) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        a = Leaf.zeros(trace)
        b = Leaf.zeros(trace)
        a.values[:] = 1.0
        b.values[:] = 2.0
        s = a + b
        d = b - a
        assert isinstance(s, Leaf) and isinstance(d, Leaf)
        np.testing.assert_array_equal(s.values, 3.0)
        np.testing.assert_array_equal(d.values, 1.0)

    def test_scalar_mul_div_neg(self, Leaf) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        a = Leaf.zeros(trace)
        a.values[:] = 4.0
        np.testing.assert_allclose((2.5 * a).values, 10.0)
        np.testing.assert_allclose((a / 4.0).values, 1.0)
        np.testing.assert_array_equal((-a).values, -4.0)
        assert isinstance(-a, Leaf)

    def test_cross_carrier_discrimination_is_trace_content(self, Leaf) -> None:
        """CS4b S3 (F2): twin carriers mix; a different quadrature's trace
        content refuses."""
        a = Leaf.zeros(angular_trace_space(slab_mesh(nx=5), posing_quadrature(slab_mesh(nx=5)), posing_scheme(slab_mesh(nx=5))))
        b = Leaf.zeros(angular_trace_space(slab_mesh(nx=5), posing_quadrature(slab_mesh(nx=5)), posing_scheme(slab_mesh(nx=5))))  # distinct instance, same content
        _ = a + b  # twin content — legal since the F2 re-key
        c = Leaf.zeros(angular_trace_space(_quad8_mesh(), Quadrature.gauss_legendre(n_ordinates=8), DiamondDifference()))
        with pytest.raises(ValueError, match="equal space"):
            _ = a + c

    def test_frozen(self, Leaf) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = Leaf.zeros(trace)
        with pytest.raises(FrozenInstanceError):
            bf.mesh = m  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# The load-bearing boundary invariant: cross-class arithmetic RAISES
# even though all three boundary leaves share the SAME AngularTraceSpace.
#
# This is structurally DIFFERENT from the bulk families, where
# cross-class operands also have different spaces/shapes (so either gate
# would catch them). Here the space gate would PASS — only the
# class-identity gate discriminates AngularBoundaryFlux / AngularBoundarySourceSink /
# AngularBoundaryResidual.
# ════════════════════════════════════════════════════════════════════


class TestCrossClassRejectionSharedSpace:
    def test_all_three_share_the_same_trace_space(self) -> None:
        """Pre-condition for the invariant: the three leaves are built on
        the IDENTICAL ``angular_trace_space(mesh, quadrature, scheme)`` object — so ``space == space`` and
        the space gate alone would NOT reject cross-class arithmetic."""
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        flux = AngularBoundaryFlux.zeros(trace)
        src = AngularBoundarySourceSink.zeros(trace)
        res = AngularBoundaryResidual.zeros(trace)
        assert flux.space is src.space is res.space is trace
        # The space gate would pass (equal spaces) — proving it is the
        # CLASS gate that must do the rejection below.
        assert flux.space == src.space == res.space

    @pytest.mark.parametrize(
        ("A", "B"),
        [
            (AngularBoundarySourceSink, AngularBoundaryFlux),
            (AngularBoundaryResidual, AngularBoundaryFlux),
            (AngularBoundarySourceSink, AngularBoundaryResidual),
        ],
    )
    def test_cross_class_add_sub_raises(self, A, B) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        a = A.zeros(trace)
        b = B.zeros(trace)
        # Same shape, same space, distinct class → class gate rejects.
        assert a.values.shape == b.values.shape
        assert a.space == b.space
        with pytest.raises(TypeError, match="same-class"):
            _ = a + b  # type: ignore[operator]
        with pytest.raises(TypeError, match="same-class"):
            _ = b - a  # type: ignore[operator]

    def test_cross_class_inner_product_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        src = AngularBoundarySourceSink.zeros(trace)
        res = AngularBoundaryResidual.zeros(trace)
        with pytest.raises(TypeError, match="same-class"):
            src.inner_product(res)  # type: ignore[arg-type]


# ════════════════════════════════════════════════════════════════════
# The prescribed-inflow generator (the ergonomic q-source constructor).
#
# Source-role-only: writes ONLY the inflow ordinate slots of the named
# faces (outflow slots of a prescribed inflow are physically meaningless
# — the sweep determines outflow), leaving everything else zero. This is
# the ergonomic specialisation of the general ``from_face_arrays`` that
# every prescribed-inflow consumer (non-vacuum MMS, splitting probe)
# previously hand-rolled as zero-allocate + ``face_view[inflow] = …``.
# ════════════════════════════════════════════════════════════════════


class TestPrescribedInflowGenerator:
    @pytest.mark.catches("ERR-047")
    def test_inflow_written_outflow_masked_absent_face_vacuum(self) -> None:
        m = slab_mesh(nx=5, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        mu = posing_quadrature(m).mu_x
        # Full (N, ng) slot with EVERY ordinate non-zero on purpose: the
        # generator MUST mask the outflow ordinates back to zero.
        vals = np.outer(mu + 2.0, np.array([1.0, 0.5]))  # (N, 2), all non-zero
        bss = AngularBoundarySourceSink.prescribed_inflow_on(trace, {"xmin": vals})
        assert isinstance(bss, AngularBoundarySourceSink)

        inflow = trace.inflow_indices_for_face("xmin")
        outflow = np.setdiff1d(np.arange(N), inflow)
        view = bss.face_view("xmin")
        # inflow ordinates carry the prescribed values ...
        np.testing.assert_array_equal(view[inflow, :], vals[inflow, :])
        # ... outflow ordinates are masked to zero (illegal-state prevention) ...
        np.testing.assert_array_equal(view[outflow, :], 0.0)
        # ... and a face not named is vacuum.
        np.testing.assert_array_equal(bss.face_view("xmax"), 0.0)

    @pytest.mark.catches("ERR-047")
    def test_sphere_single_face(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        vals = np.full((N, 2), 3.0)
        bss = AngularBoundarySourceSink.prescribed_inflow_on(trace, {"xmax": vals})
        inflow = trace.inflow_indices_for_face("xmax")
        outflow = np.setdiff1d(np.arange(N), inflow)
        view = bss.face_view("xmax")
        np.testing.assert_array_equal(view[inflow, :], 3.0)
        np.testing.assert_array_equal(view[outflow, :], 0.0)

    def test_unknown_face_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        with pytest.raises(ValueError, match="not a"):
            AngularBoundarySourceSink.prescribed_inflow_on(trace, {"nope": np.zeros((N, 2))})

    def test_shape_mismatch_raises(self) -> None:
        m = slab_mesh(nx=5, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        with pytest.raises(ValueError, match="shape"):
            AngularBoundarySourceSink.prescribed_inflow_on(trace, {"xmin": np.zeros((N, 1))})

    def test_generator_is_source_role_only(self) -> None:
        """The prescribed-inflow generator is the SOURCE-role leaf's verb;
        the flux / residual leaves do not carry it (a flux trace is swept,
        a residual is differenced — neither is a prescribed source)."""
        for verb in ("prescribed_inflow", "prescribed_inflow_on"):
            assert hasattr(AngularBoundarySourceSink, verb)
            assert not hasattr(AngularBoundaryFlux, verb)
            assert not hasattr(AngularBoundaryResidual, verb)

    def test_a_scalar_inflow_onto_a_moment_slot_seeds_the_average(self) -> None:
        """A face-AVERAGE ``(N, ng, *transverse)`` array prescribed onto a linear-discontinuous trace, whose 2-D slot is
        ``(N, ng, *transverse, 2)`` (#251), lands in the AVERAGE moment of the inflow rows; the transverse slope moment
        and every outflow row stay exactly zero (layer DAG Phase 4 spec §2 file 4: no row reached this branch of
        ``prescribed_inflow_on`` before; first red, the battery arm ``inflow_scalar_onto_last_moment``)."""
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), LinearDiscontinuous())
        N = posing_quadrature(m).n_ordinates
        slot = trace.layout.faces["xmin"].shape
        assert slot == (N, 2, 2, 2), f"precondition: an LD 2-D face slot carries the transverse moment, got {slot}"
        rng = np.random.default_rng(5)
        vals = rng.uniform(0.5, 1.5, size=(N, 2, 2))  # (N, ng, ny): every entry non-zero
        bss = AngularBoundarySourceSink.prescribed_inflow_on(trace, {"xmin": vals})
        view = bss.face_view("xmin")
        inflow = trace.inflow_indices_for_face("xmin")
        outflow = np.setdiff1d(np.arange(N), inflow)
        assert inflow.size > 0 and outflow.size > 0
        np.testing.assert_array_equal(view[inflow, ..., 0], vals[inflow])
        np.testing.assert_array_equal(view[inflow, ..., 1], 0.0)
        np.testing.assert_array_equal(view[outflow], 0.0)
