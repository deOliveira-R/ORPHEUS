r"""L0/L2 — pure-Field :class:`AngularBoundaryFlux` algebra, FaceLayout slice views,
and flat-buffer round-trips.

Pins the post-D-G contract for
:class:`orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`:

* Inherits :class:`~orpheus.numerics.field.Field`; no hand-coded
  dunders.
* Flat backing buffer + :class:`FaceLayout` produces per-face slice
  views with NO copies (memory-shared).
* Mesh-binding rejection (cross-mesh arithmetic raises even when
  spaces match).
* Round-trip equivalence: face-view writes ↔ flat-buffer reads.
* Frozen contract: attribute assignment raises ``FrozenInstanceError``.
* Geometry-conditional constructors (slab / sphere / 2-D Cartesian).
* ``from_face_arrays`` builder pattern (the post-D-G replacement for
  pre-D-G per-attribute mutation).

References
----------

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`
  §3.1 (the specification for this module).
* ``vv-principles`` "Bit-identity vs principled-equivalence" — this
  test is the bit-identity gate at the FIELD-ALGEBRA layer.

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
from orpheus.transport.spatial import DiamondDifference
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from tests.gates.transport._phase_space_cases import posing_quadrature, posing_scheme
from tests.gates.transport._phase_space_cases import slab_mesh


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures (mirrored from tests/gates/sn/test_typed_fields.py)
# ───────────────────────────────────────────────────────────────────────


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


# ───────────────────────────────────────────────────────────────────────
# A. Field-ABC inheritance and algebra
# ───────────────────────────────────────────────────────────────────────


class TestFieldAlgebraInherited:
    def test_inherits_field(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert isinstance(bf, Field)

    def test_flux_add_and_sub_are_plain_vector_ops(self) -> None:
        """Flux lives in V (campaign 1 CS3): bf + bf and bf − bf are the
        plain vector ops returning fresh SAME-typed fields, originals
        unchanged (until 2026-08-19 the #208 affine gate raised on + and
        minted an AngularBoundaryDisplacement on −)."""
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf1 = AngularBoundaryFlux.zeros(trace)
        bf2 = AngularBoundaryFlux.zeros(trace)
        bf1.values[:] = 3.0
        bf2.values[:] = 2.0
        s = bf1 + bf2
        if type(s) is not AngularBoundaryFlux:
            raise AssertionError("bf + bf left the leaf type")
        np.testing.assert_array_equal(s.values, 5.0)
        out = bf1 - bf2
        if out is bf1 or out is bf2:
            raise AssertionError("subtraction returned an operand, not a copy")
        if type(out) is not AngularBoundaryFlux:
            raise AssertionError("bf − bf left the leaf type")
        np.testing.assert_array_equal(out.values, 1.0)
        np.testing.assert_array_equal(bf1.values, 3.0)  # originals unchanged
        np.testing.assert_array_equal(bf2.values, 2.0)

    def test_sub_sphere(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf1 = AngularBoundaryFlux.zeros(trace)
        bf2 = AngularBoundaryFlux.zeros(trace)
        bf1.values[:] = 5.0
        bf2.values[:] = 2.0
        out = bf1 - bf2
        np.testing.assert_array_equal(out.values, 3.0)

    def test_scalar_mul_2d(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        bf.values[:] = 3.0
        out = bf * 2.5
        np.testing.assert_allclose(out.values, 7.5)

    def test_scalar_div_propagates(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        bf.values[:] = 10.0
        out = bf / 4.0
        np.testing.assert_allclose(out.values, 2.5)

    def test_neg_2d(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        bf.values[:] = 1.5
        out = -bf
        np.testing.assert_array_equal(out.values, -1.5)

    def test_displacement_distributive_property_2d(self) -> None:
        """L1: distributivity (d1 + d2) * c == d1*c + d2*c in V — since the
        CS3 cone carve the differences are ordinary same-typed boundary
        fluxes, so this is now a law of the ONE flux vector space."""
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        base = AngularBoundaryFlux.zeros(trace)
        bf1 = AngularBoundaryFlux.zeros(trace)
        bf2 = AngularBoundaryFlux.zeros(trace)
        rng = np.random.default_rng(0)
        bf1.values[:] = rng.standard_normal(bf1.values.shape)
        bf2.values[:] = rng.standard_normal(bf2.values.shape)
        d1 = bf1 - base  # same-class signed difference (V — campaign 1 CS3)
        d2 = bf2 - base
        left = (d1 + d2) * 1.7
        right = d1 * 1.7 + d2 * 1.7
        # Distributivity holds in real arithmetic but accumulates ULP drift
        # in IEEE-754 (FP non-associativity). Per vv-principles, bound is
        # (reduction depth) × ULP — for this 2x2-elem combine, ~1e-14.
        np.testing.assert_allclose(left.values, right.values, rtol=1e-13)

    def test_inner_product_with_self_is_l2_squared(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        bf.values[:] = 1.5
        # inner_product(x, x) == l2(x)² under the space metric. The AngularTraceSpace
        # now carries the |Ω·n|·w partial-current weights (Phase 4.1 G-adjoint),
        # so this is the weighted norm, NOT the Euclidean N·x² — assert the
        # genuine (metric-robust) identity against the space-induced l2.
        np.testing.assert_allclose(bf.inner_product(bf), bf.l2 ** 2)


# ───────────────────────────────────────────────────────────────────────
# B. FaceLayout — slice views, no copies, all geometries
# ───────────────────────────────────────────────────────────────────────


class TestFaceLayoutSliceViews:
    def test_slab_layout_has_xmin_and_xmax(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert set(bf.layout.faces) == {"xmin", "xmax"}

    def test_sphere_layout_has_only_xmax(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert set(bf.layout.faces) == {"xmax"}

    def test_cartesian_2d_layout_has_four_faces(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert set(bf.layout.faces) == {"xmin", "xmax", "ymin", "ymax"}

    def test_face_view_is_memory_shared_slab(self) -> None:
        """face_view returns a view, NOT a copy — writes propagate."""
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        xmin_view = bf.face_view("xmin")
        # Write through the view; flat buffer reflects the change.
        xmin_view[...] = 7.0
        flat_start = bf.layout.faces["xmin"].offset
        flat_end = flat_start + bf.layout.faces["xmin"].flat_size
        np.testing.assert_array_equal(bf.values[flat_start:flat_end], 7.0)

    def test_face_view_is_memory_shared_2d_xmin(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        xmin_view = bf.face_view("xmin")
        N = posing_quadrature(m).n_ordinates
        assert xmin_view.shape == (N, m.ng, m.spatial_shape[1])
        xmin_view[0, 0, 0] = 9.0
        flat_idx = bf.layout.faces["xmin"].offset
        assert bf.values[flat_idx] == 9.0

    def test_face_shapes_match_per_geometry(self) -> None:
        # 1-D slab.
        m = slab_mesh(nx=4, ng=2)
        N = posing_quadrature(m).n_ordinates
        bf = AngularBoundaryFlux.zeros(angular_trace_space(m, posing_quadrature(m), posing_scheme(m)))
        assert bf.face_view("xmin").shape == (N, 2)
        assert bf.face_view("xmax").shape == (N, 2)

        # 1-D sphere.
        m = _sphere_mesh(nx=4, ng=2)
        N = posing_quadrature(m).n_ordinates
        bf = AngularBoundaryFlux.zeros(angular_trace_space(m, posing_quadrature(m), posing_scheme(m)))
        assert bf.face_view("xmax").shape == (N, 2)

        # 2-D Cartesian.
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        N = posing_quadrature(m).n_ordinates
        bf = AngularBoundaryFlux.zeros(angular_trace_space(m, posing_quadrature(m), posing_scheme(m)))
        assert bf.face_view("xmin").shape == (N, 2, 2)
        assert bf.face_view("xmax").shape == (N, 2, 2)
        assert bf.face_view("ymin").shape == (N, 2, 3)
        assert bf.face_view("ymax").shape == (N, 2, 3)

    def test_face_views_dict_keys_match_layout(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        views = bf.face_views
        assert set(views) == set(bf.layout.faces)

    def test_unknown_face_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        with pytest.raises(KeyError, match="ymin"):
            bf.face_view("ymin")

    def test_total_size_consistent_with_face_sizes_2d(self) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        total = sum(slot.flat_size for slot in bf.layout.faces.values())
        assert bf.layout.total_size == total
        assert bf.values.shape == (total,)

    def test_no_interior_cells_in_layout_2d(self) -> None:
        """CRITICAL: 2-D layout contains ONLY face slots — interior
        wavefront cache lives on SweepScratch, NOT on AngularBoundaryFlux.

        Catches Rank 5 failure mode per the verification memo.
        """
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        bf = AngularBoundaryFlux.zeros(trace)
        # Post-D-G face-only total:
        face_only = 2 * N * 2 * 2 + 2 * N * 2 * 3  # x-faces + y-faces
        assert bf.values.size == face_only
        # Pre-D-G conflated total would be:
        legacy = N * 2 * (3 + 1) * 2 + N * 2 * 3 * (2 + 1)  # (nx+1)*ny + nx*(ny+1)
        assert bf.values.size < legacy


# ───────────────────────────────────────────────────────────────────────
# C. Mesh-binding rejection (preserved from legacy)
# ───────────────────────────────────────────────────────────────────────


class TestMeshBindingRejection:
    def test_cross_carrier_discrimination_is_trace_content(self) -> None:
        """CS4b S3 (F2): twin carriers mint EQUAL trace spaces (the
        content digest agrees) and boundary fluxes mix across them; a
        carrier with a different quadrature mints an UNEQUAL trace and
        both binary ops refuse."""
        m1 = slab_mesh(nx=5)
        trace_m1 = angular_trace_space(m1, posing_quadrature(m1), posing_scheme(m1))
        m2 = slab_mesh(nx=5)   # different instance, same content
        bf1 = AngularBoundaryFlux.zeros(trace_m1)
        bf2 = AngularBoundaryFlux.zeros(angular_trace_space(m2, posing_quadrature(m2), posing_scheme(m2)))
        _ = bf1 - bf2  # twin content — legal since the F2 re-key
        _ = bf1 + bf2
        bf3 = AngularBoundaryFlux.zeros(angular_trace_space(_quad8_mesh(), Quadrature.gauss_legendre(n_ordinates=8), DiamondDifference()))
        with pytest.raises(ValueError, match="equal space"):
            bf1 - bf3
        with pytest.raises(ValueError, match="equal space"):
            bf1 + bf3

    def test_wrong_type_add_rejected(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        with pytest.raises(TypeError, match="same-class"):
            bf + 42  # type: ignore[operator]


# ───────────────────────────────────────────────────────────────────────
# D. Flat-buffer round-trip (the new invariant)
# ───────────────────────────────────────────────────────────────────────


class TestFlatBufferRoundTrip:
    def test_face_write_reflected_in_flat_buffer(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        bf.face_view("xmin")[...] = 1.0
        bf.face_view("xmax")[...] = 2.0
        xmin_off = bf.layout.faces["xmin"].offset
        xmin_size = bf.layout.faces["xmin"].flat_size
        xmax_off = bf.layout.faces["xmax"].offset
        xmax_size = bf.layout.faces["xmax"].flat_size
        np.testing.assert_array_equal(bf.values[xmin_off:xmin_off + xmin_size], 1.0)
        np.testing.assert_array_equal(bf.values[xmax_off:xmax_off + xmax_size], 2.0)

    def test_flat_write_reflected_in_face(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        # Write something into the flat buffer; verify face view sees it.
        bf.values[0] = 8.0
        # xmin is the first face → element 0 maps to xmin face's [0, 0, 0].
        assert bf.face_view("xmin")[0, 0, 0] == 8.0

    def test_round_trip_arithmetic_preserves_face_values(self) -> None:
        """After ``out = bf1 − bf2`` every face of the resulting (same-typed)
        difference equals the per-face difference — bit-identical at the
        face-value level, the difference is face-structured."""
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        rng = np.random.default_rng(42)
        bf1 = AngularBoundaryFlux.zeros(trace)
        bf2 = AngularBoundaryFlux.zeros(trace)
        bf1.values[:] = rng.standard_normal(bf1.values.shape)
        bf2.values[:] = rng.standard_normal(bf2.values.shape)
        out = bf1 - bf2  # same-class signed difference — shares the trace layout
        for name in bf1.layout.faces:
            np.testing.assert_array_equal(
                out.face_view(name),
                bf1.face_view(name) - bf2.face_view(name),
            )


# ───────────────────────────────────────────────────────────────────────
# E. Construction factories
# ───────────────────────────────────────────────────────────────────────


class TestConstruction:
    def test_zeros_slab(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert isinstance(bf, AngularBoundaryFlux)
        np.testing.assert_array_equal(bf.values, 0.0)
        assert bf.space is trace

    def test_zeros_sphere(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert "xmin" not in bf.layout.faces
        np.testing.assert_array_equal(bf.values, 0.0)

    def test_zeros_2d(self) -> None:
        m = _cartesian_2d_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        assert set(bf.layout.faces) == {"xmin", "xmax", "ymin", "ymax"}

    def test_post_init_validates_shape(self) -> None:
        """Constructing with values.shape != (layout.total_size,) raises."""
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        with pytest.raises(ValueError, match="values.shape"):
            AngularBoundaryFlux(
                values=np.zeros(trace.shape[0] + 1),
                space=trace,
            )

    def test_post_init_requires_trace_space(self) -> None:
        """A.5: the space MUST be an AngularTraceSpace carrying a FaceLayout —
        a layout-less plain FunctionSpace cannot back an AngularBoundaryFlux
        (illegal-states-unrepresentable; ``face_view`` needs the layout)."""
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        from orpheus.numerics.space import FunctionSpace
        plain = FunctionSpace(name="sn_boundary_flat", shape=(trace.shape[0],))
        with pytest.raises(TypeError, match="AngularTraceSpace"):
            AngularBoundaryFlux(values=np.zeros(trace.shape[0]), space=plain)

    def test_from_face_arrays_slab(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        xmin_arr = np.full((N, 2), 1.0)
        xmax_arr = np.full((N, 2), 2.0)
        bf = AngularBoundaryFlux.from_face_arrays_on(trace, {"xmin": xmin_arr, "xmax": xmax_arr})
        np.testing.assert_array_equal(bf.face_view("xmin"), 1.0)
        np.testing.assert_array_equal(bf.face_view("xmax"), 2.0)

    def test_from_face_arrays_2d(self) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        face_arrays = {
            "xmin": np.full((N, 2, 2), 1.0),
            "xmax": np.full((N, 2, 2), 2.0),
            "ymin": np.full((N, 2, 3), 3.0),
            "ymax": np.full((N, 2, 3), 4.0),
        }
        bf = AngularBoundaryFlux.from_face_arrays_on(trace, face_arrays)
        np.testing.assert_array_equal(bf.face_view("xmin"), 1.0)
        np.testing.assert_array_equal(bf.face_view("ymax"), 4.0)

    def test_from_face_arrays_missing_face_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        # Slab needs both xmin and xmax; omit xmin.
        with pytest.raises(ValueError, match="do not match"):
            AngularBoundaryFlux.from_face_arrays_on(trace, {"xmax": np.zeros((N, 2))})

    def test_from_face_arrays_extra_face_raises(self) -> None:
        m = _sphere_mesh()
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        # Sphere only has xmax; passing xmin should raise.
        with pytest.raises(ValueError, match="do not match"):
            AngularBoundaryFlux.from_face_arrays_on(trace, {"xmax": np.zeros((N, 2)), "xmin": np.zeros((N, 2))},
            )

    def test_from_face_arrays_wrong_shape_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        N = posing_quadrature(m).n_ordinates
        with pytest.raises(ValueError, match="shape"):
            AngularBoundaryFlux.from_face_arrays_on(trace, {
                "xmin": np.zeros((N, 3)),   # ng=3 instead of 2
                "xmax": np.zeros((N, 2)),
            })


# ───────────────────────────────────────────────────────────────────────
# F. Frozen-instance gate
# ───────────────────────────────────────────────────────────────────────


class TestFrozenInstance:
    def test_assign_values_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        with pytest.raises(FrozenInstanceError):
            bf.values = np.zeros(bf.values.shape)  # type: ignore[misc]

    def test_assign_layout_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        with pytest.raises(FrozenInstanceError):
            bf.layout = bf.layout  # type: ignore[misc]

    def test_assign_mesh_raises(self) -> None:
        m = slab_mesh(nx=5)
        trace = angular_trace_space(m, posing_quadrature(m), posing_scheme(m))
        bf = AngularBoundaryFlux.zeros(trace)
        with pytest.raises(FrozenInstanceError):
            bf.mesh = m  # type: ignore[misc]
