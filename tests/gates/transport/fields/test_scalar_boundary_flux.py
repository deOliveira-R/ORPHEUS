"""Intrinsic-property gates for the scalar trace substrate (#290 P2).

Every math-bearing type ships a test of its DEFINING laws (project
standard): here the :class:`ScalarTraceSpace` space laws, the
:class:`ScalarBoundaryFlux` P1 partial-current dictionary

    J = J⁺ − J⁻          (net outward current)
    φ_Γ = 2(J⁺ + J⁻)     (P1 boundary scalar flux)

the albedo-family identities ``J⁻ = 𝒜·J⁺`` (vacuum 𝒜=0 / reflective
𝒜=1 / zero-flux Dirichlet 𝒜=−1), the flat-buffer/view storage
contract, the field-algebra guards inherited from ``BoundaryField``, and
the widened ``FullField``/``FullFieldSpace`` composite admission.

The carrier is the :class:`~orpheus.diffusion.augmented_mesh.MaterialMesh`
phase space (#290 P7a): a boundary trace is method behavior, so the
scalar trace family binds to the diffusion method-mesh (a bare
``MaterialMesh`` owns no trace — the negative gate below).

Convention contract: ``.claude/plans/diffusion_crosswalk.md``.
Foundation tier — software invariants on the new types; no
``verifies(...)``.

Layer DAG Phase 4 (2026-11-01): the file builds material meshes and spells each space it reads from the raw
constructors of orpheus.transport.phase_spaces (the SN problem's and the diffusion mesh's cached mints are those
constructors on the problem's own ingredients, gated in
tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py and
tests/gates/diffusion/test_scalar_trace_is_the_raw_constructor.py).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import scalar_trace_space
from tests.gates.transport._phase_space_cases import scalar_full_field_space
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces import FullFieldSpace, ScalarTraceSpace
from orpheus.transport.fields import ScalarBoundaryFlux, ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_mesh import MaterialMesh

pytestmark = [pytest.mark.foundation]


def _slab_mesh(nx: int = 4, width: float = 10.0) -> MaterialMesh:
    mesh1d = Mesher(StructuredGeometry.from_homogeneous(width, MirrorDeck())).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    return MaterialMesh(mesh1d, {0: get_mixture("A", "2g")})


def _one_group_mesh(nx: int = 4, width: float = 10.0) -> MaterialMesh:
    """1-group sibling — the scalar-trace CONTENT (and shape) differs, so
    the carrier mints an UNEQUAL trace space (F2)."""
    mesh1d = Mesher(StructuredGeometry.from_homogeneous(width, MirrorDeck())).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    return MaterialMesh(mesh1d, {0: get_mixture("A", "1g")})


# ─────────────────────────────────────────────────────────────────────
# ScalarTraceSpace laws
# ─────────────────────────────────────────────────────────────────────


class TestScalarTraceSpace:
    def test_slab_dim_faces_and_area_metric(self):
        """1-D slab: 2 faces × 2 components × ng DOFs; unit-area metric."""
        ts = scalar_trace_space(_slab_mesh())
        np.testing.assert_array_equal(ts.shape, (2 * 2 * 2,))
        if ts.face_names != ("xmin", "xmax"):
            pytest.fail(f"slab faces should be (xmin, xmax); got {ts.face_names}")
        if ts.slot_shape("xmax") != (2, 2):
            pytest.fail(f"slot shape should be (2, ng); got {ts.slot_shape('xmax')}")
        # Slab boundary areas are exactly 1 → the metric is exactly ones.
        np.testing.assert_array_equal(
            ts.inner_product_weights, np.ones(ts.shape)
        )

    def test_sphere_has_single_outer_face_with_surface_area(self):
        """The pole is NOT a face; the metric is the outer surface 4πR²."""
        r_outer = 5.0
        geometry = StructuredGeometry.sphere((0.0, r_outer), (0,), outer=0)
        mesh1d = Mesher(geometry).partition(CellsByCount.uniform_width(3)).mesh
        # A curved face has no mirror deck: the perfect mirror is the specular wall tag (rule R2).
        mm = MaterialMesh(mesh1d, {0: get_mixture("A", "2g")})
        ts = scalar_trace_space(mm)
        if ts.face_names != ("xmax",):
            pytest.fail(f"sphere trace must have only xmax; got {ts.face_names}")
        # Same-source claim: the metric IS mesh.areas[-1] broadcast.
        np.testing.assert_array_equal(
            ts.inner_product_weights,
            np.full(ts.shape, mm.areas[-1]),
        )

    def test_identity_is_name_and_shape(self):
        """Two same-shape scalar traces compare equal (metadata excluded)."""
        a = scalar_trace_space(_slab_mesh(nx=4))
        b = scalar_trace_space(_slab_mesh(nx=7))  # same boundary, different interior
        if a != b:
            pytest.fail("(name, shape)-identical scalar traces must be equal")
        if a == FunctionSpace(name="scalar_flux", shape=a.shape):
            pytest.fail("different-name space must not compare equal")

    def test_for_faces_rejects_bad_areas(self):
        with pytest.raises(ValueError, match="face_areas"):
            ScalarTraceSpace.for_faces(
                [("xmin", ()), ("xmax", ())], 2, {"xmax": 1.0},
            )
        with pytest.raises(ValueError, match="positive"):
            ScalarTraceSpace.for_faces(
                [("xmin", ()), ("xmax", ())], 2, {"xmin": 0.0, "xmax": 1.0},
            )


# ─────────────────────────────────────────────────────────────────────
# ScalarBoundaryFlux — the P1 dictionary and the storage contract
# ─────────────────────────────────────────────────────────────────────


class TestPartialCurrentLaws:
    def test_p1_dictionary_identities(self):
        """J = J⁺ − J⁻ and φ_Γ = 2(J⁺ + J⁻), exactly (same float ops)."""
        mm = _slab_mesh()
        j_plus = np.array([[3.0, 1.0]])   # placed on OUTFLOW_ROW below
        j_minus = np.array([[1.0, 0.25]])
        pc = ScalarBoundaryFlux.from_face_arrays_on(scalar_trace_space(mm),
            {
                "xmin": np.zeros((2, 2)),
                "xmax": np.vstack([j_plus, j_minus]),
            },
        )
        np.testing.assert_array_equal(
            pc.net_current("xmax"), (j_plus - j_minus)[0]
        )
        np.testing.assert_array_equal(
            pc.p1_boundary_scalar_flux("xmax"), 2.0 * (j_plus + j_minus)[0]
        )

    def test_albedo_family_identities(self):
        """vacuum 𝒜=0 ⟹ φ_Γ = 2J; reflective 𝒜=1 ⟹ J = 0; zero-flux
        𝒜=−1 ⟹ φ_Γ = 0 — the three named members of J⁻ = 𝒜·J⁺."""
        mm = _slab_mesh()
        j_plus = np.array([2.0, 0.5])

        def with_albedo(albedo: float) -> ScalarBoundaryFlux:
            slot = np.vstack([j_plus, albedo * j_plus])
            return ScalarBoundaryFlux.from_face_arrays_on(scalar_trace_space(mm), {"xmin": np.zeros((2, 2)), "xmax": slot},
            )

        vacuum = with_albedo(0.0)
        np.testing.assert_array_equal(vacuum.inflow_view("xmax"), np.zeros(2))
        np.testing.assert_array_equal(
            vacuum.p1_boundary_scalar_flux("xmax"),
            2.0 * vacuum.net_current("xmax"),
        )
        reflective = with_albedo(1.0)
        np.testing.assert_array_equal(
            reflective.net_current("xmax"), np.zeros(2)
        )
        zero_flux = with_albedo(-1.0)
        np.testing.assert_array_equal(
            zero_flux.p1_boundary_scalar_flux("xmax"), np.zeros(2)
        )

    def test_views_share_backing_buffer(self):
        """outflow/inflow views are zero-copy into the flat buffer."""
        mm = _slab_mesh()
        pc = ScalarBoundaryFlux.zeros(scalar_trace_space(mm))
        pc.outflow_view("xmin")[:] = 7.0
        np.testing.assert_array_equal(
            pc.face_view("xmin")[ScalarTraceSpace.OUTFLOW_ROW],
            np.full(2, 7.0),
        )
        if float(pc.values.sum()) != 14.0:
            pytest.fail("view write did not land in the backing buffer")


class TestPartialCurrentGuards:
    def test_requires_scalar_trace_space(self):
        """An angular/bare space is rejected with the family message."""
        mm = _slab_mesh()
        with pytest.raises(TypeError, match="ScalarTraceSpace"):
            ScalarBoundaryFlux(
                values=np.zeros(8),
                space=FunctionSpace(name="scalar_trace", shape=(8,)),
            )

    def test_space_content_guard(self):
        """CS4b S3 (F2): twin carriers mint EQUAL scalar traces and mix; a
        different group structure mints an UNEQUAL trace and refuses."""
        a = ScalarBoundaryFlux.zeros(scalar_trace_space(_slab_mesh()))
        b = ScalarBoundaryFlux.zeros(scalar_trace_space(_slab_mesh()))
        _ = a - b  # twin content — legal since the F2 re-key
        c = ScalarBoundaryFlux.zeros(scalar_trace_space(_one_group_mesh()))
        with pytest.raises(ValueError, match="equal space"):
            _ = a - c

    def test_bare_material_mesh_is_refused(self):
        """A boundary trace is method BEHAVIOR (#290 P7a): the bare
        MaterialMesh data carrier owns no scalar trace, so a trace
        field cannot be built on one — the family diagnosis points at
        the promotion."""
        mesh1d = Mesher(StructuredGeometry.from_homogeneous(10.0, MirrorDeck())).partition(
            CellsByCount.uniform_width(4),
        ).mesh
        mm = MaterialMesh(mesh1d, {0: get_mixture("A", "2g")})
        # Post-S5 the family diagnosis lives on the space hook (the surface
        # space_on and the keyed factories ride).
        with pytest.raises(ValueError, match="no scalar trace.*MaterialMesh"):
            ScalarBoundaryFlux._face_space_of(mm)

    def test_vector_algebra(self):
        """The V algebra on the scalar trace (campaign 1 CS3): state − state
        and state + state are plain same-typed vector ops (until 2026-08-19
        the #208 affine discipline minted a ScalarBoundaryDisplacement on −
        and refused +)."""
        mm = _slab_mesh()
        a = ScalarBoundaryFlux.zeros(scalar_trace_space(mm))
        b = ScalarBoundaryFlux.from_face_arrays_on(scalar_trace_space(mm), {"xmin": np.ones((2, 2)), "xmax": np.ones((2, 2))},
        )
        step = b - a
        if type(step) is not ScalarBoundaryFlux:
            pytest.fail(f"state − state left the leaf type: {type(step).__name__}")
        np.testing.assert_array_equal(step.values, b.values - a.values)
        moved = a + step
        if type(moved) is not ScalarBoundaryFlux:
            pytest.fail("the update step left the leaf type")
        np.testing.assert_array_equal(moved.values, b.values)
        s = a + b
        if type(s) is not ScalarBoundaryFlux:
            pytest.fail("state + state left the leaf type")
        np.testing.assert_array_equal(s.values, a.values + b.values)

    def test_from_face_arrays_rejects_wrong_faces(self):
        mm = _slab_mesh()
        with pytest.raises(ValueError, match="face_arrays keys"):
            ScalarBoundaryFlux.from_face_arrays_on(scalar_trace_space(mm), {"xmax": np.zeros((2, 2))})


# ─────────────────────────────────────────────────────────────────────
# The widened composite — FullField / FullFieldSpace admission
# ─────────────────────────────────────────────────────────────────────


class TestScalarComposite:
    def test_scalar_composite_is_an_element_of_the_carrier_composite_space(self):
        mm = _slab_mesh()
        full = FullField(
            interior=ScalarFlux.zeros(mm.bulk_space), boundary=ScalarBoundaryFlux.zeros(scalar_trace_space(mm)),
        )
        # CS4b S4 (F4): the composite's derived space content-equals the
        # carrier's cached mint (== not is — the wrapper is per-instance;
        # the MEMBERS are the carrier's cached objects).
        if full.space != scalar_full_field_space(mm):
            pytest.fail("composite space must content-equal the carrier mint")
        if full.space.interior_space is not scalar_full_field_space(mm).interior_space:
            pytest.fail("interior member must BE the carrier's cached mint")

    def test_twin_mesh_composite_blocks_mix(self):
        """CS4b S4 (F4 — the F2 doctrine's composite leg): two separately
        constructed but content-identical carriers mint content-equal
        spaces, so a composite may be assembled across them — the pre-S4
        ``interior.mesh is boundary.mesh`` refusal was mesh-IDENTITY
        doctrine, retired tree-wide at S3. Cross-slot coherence is not a
        construction question (no carrier reference exists there); the
        operator/solution admission seams hold the reference and refuse
        content mismatches (their own witnesses discriminate volumes /
        quadrature / family)."""
        a, b = _slab_mesh(), _slab_mesh()
        full = FullField(
            interior=ScalarFlux.zeros(a.bulk_space),
            boundary=ScalarBoundaryFlux.zeros(scalar_trace_space(b)),
        )
        if full.space != scalar_full_field_space(a):
            pytest.fail("twin-carrier composite must content-equal either mint")

    def test_direct_sum_space_metric_and_inner_product(self):
        """from_blocks admits the scalar trace; the direct-sum inner
        product is the bulk + trace block sum (hand value)."""
        mm = _slab_mesh()
        bulk = ScalarFlux.zeros(mm.bulk_space)
        pc = ScalarBoundaryFlux.zeros(scalar_trace_space(mm))
        pc.outflow_view("xmax")[:] = [3.0, 1.0]
        pc.inflow_view("xmax")[:] = [1.0, 0.25]
        full = FullField(interior=bulk, boundary=pc)
        ffs = FullFieldSpace.from_blocks(bulk.space, scalar_trace_space(mm))
        np.testing.assert_array_equal(
            ffs.shape, (bulk.values.size + pc.values.size,)
        )
        # Slab metric is exactly 1 ⟹ ⟨x,x⟩ = Σ values² = 9+1+1+0.0625.
        if ffs.inner_product(full, full) != 11.0625:
            pytest.fail(
                f"direct-sum inner product wrong: {ffs.inner_product(full, full)}"
            )
        # Metric application is the identity on the unit-area slab trace.
        weighted = ffs.apply_metric(full)
        np.testing.assert_array_equal(weighted.boundary.values, pc.values)
