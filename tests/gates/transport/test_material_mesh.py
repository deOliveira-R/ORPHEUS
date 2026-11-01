r"""MaterialMesh — the method-agnostic mesh + materials data carrier.

:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` is the
"mesh + materials" middle type between geometry ``Mesh1D`` (material
*ids*, no cross sections) and ``SNProblem`` (mesh + materials + quadrature +
sweep machinery).  These tests pin its **intrinsic data contract** (the
``coding-elegance`` standard: a type ships a test of its defining
invariants); the SN side of the data/behaviour split (``SNProblem`` is-a ``MaterialMesh``,
its data block bit-identical to a standalone carrier, the ``from_material_mesh`` promotion) is
``tests/gates/sn/mesh/test_problem_is_a_material_mesh.py`` since layer DAG Phase 4:

* the carrier holds mesh + materials and exposes the method-agnostic data
  accessors (``ng`` / ``volumes`` / ``volume_measure`` / ``ndim`` /
  ``spatial_shape`` / ``mat_map`` / ``material_xs_field``);
* ``ng`` consistency is enforced at construction;
* a materials declaration covers every referenced id.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.data.materials import InconsistentMaterialsError
from orpheus.transport.mesh import (
    MaterialMesh,
    MaterialXSField,
)
from orpheus.geometry.deck import MirrorDeck
from orpheus.numerics.axis import Axis, BasisKind, Declared

pytestmark = pytest.mark.foundation


def _mix(sig_t, *, ng):
    """A minimal balanced-free synthetic Mixture with a given ng."""
    sig_t = np.asarray(sig_t, dtype=float)
    return Mixture(
        SigC=0.5 * sig_t, SigL=np.zeros(ng), SigF=np.zeros(ng),
        SigP=np.zeros(ng), SigT=sig_t,
        SigS=[csr_matrix(np.diag(0.5 * sig_t))], Sig2=[csr_matrix((ng, ng))],
        chi=np.zeros(ng), eg=None,
    )


def _two_material_mesh(ng=2):
    mats = {0: _mix([1.0, 1.5][:ng], ng=ng), 1: _mix([1.2, 1.8][:ng], ng=ng)}
    geometry = StructuredGeometry.slab(
        (0.0, 2.0, 4.0, 5.0), (0, 1, 0), left=MirrorDeck(), right=MirrorDeck(),
    )
    mesh = Mesher(geometry).partition((
        CellsByCount.uniform_width(2),
        CellsByCount.uniform_width(2),
        CellsByCount.uniform_width(1),
    )).mesh
    return mesh, mats


# ── Carrier data contract ─────────────────────────────────────────────

def test_material_mesh_holds_mesh_and_materials():
    mesh, mats = _two_material_mesh()
    mm = MaterialMesh(mesh, mats)
    # Un-weld arc (R20/R21): the attribute is the parsed stage-1
    # DECLARATION, not an alias of the caller's dict — same entries by
    # identity, read-only mapping surface (the old aliasing let callers
    # mutate the mesh's materials after construction).
    from orpheus.data.materials import Materials
    assert isinstance(mm.materials, Materials)
    assert all(mm.materials[k] is mats[k] for k in mats)
    assert mm.ng == 2
    assert mm.ndim == 1
    assert mm.spatial_shape == (5,)
    np.testing.assert_array_equal(mm.mat_map, mesh.mat_ids)
    np.testing.assert_array_equal(mm.volumes, mesh.volumes)


def test_material_mesh_volume_measure_matches_mesh():
    mesh, mats = _two_material_mesh()
    mm = MaterialMesh(mesh, mats)
    mu = mm.volume_measure
    np.testing.assert_array_equal(mu.weights, mesh.volume_measure.weights)


def test_material_mesh_builds_xs_field():
    mesh, mats = _two_material_mesh()
    mm = MaterialMesh(mesh, mats)
    field = mm.mat_xs
    assert isinstance(field, MaterialXSField)
    # per-cell SigT view follows the mat_map (cell 2 is material 1).
    sig_t = field.total_cross_section  # (ng, nx)
    assert sig_t.shape == (2, 5)
    np.testing.assert_allclose(sig_t[:, 2], mats[1].SigT)
    np.testing.assert_allclose(sig_t[:, 0], mats[0].SigT)


# ── ng-consistency invariant (the defining law) ───────────────────────

def test_inconsistent_ng_raises_at_construction():
    geometry = StructuredGeometry.slab(
        (0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=MirrorDeck(),
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh
    mats = {0: _mix([1.0, 1.5], ng=2), 1: _mix([1.0], ng=1)}
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        MaterialMesh(mesh, mats)


def test_missing_material_id_raises():
    geometry = StructuredGeometry.slab(
        (0.0, 1.0, 2.0), (0, 7), left=MirrorDeck(), right=MirrorDeck(),
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh
    with pytest.raises(ValueError, match="references material ids"):
        MaterialMesh(mesh, {0: _mix([1.0, 1.0], ng=2)})


# ── Materials declaration vs assignment ──────────────────────────────
# Until the CS4c coda (2026-09-08) this section pinned the defining
# invariants of ``MaterialMesh.from_materials`` — the mesh-less one-cell
# carrier of the infinite homogeneous medium. The factory retired with the
# coda (the homogeneous problem poses on its own space and builds no
# carrier); the ONE live guard those rows exercised survives here with its
# witness on a genuine one-cell carrier.


def test_a_materials_declaration_must_cover_every_referenced_id():
    """A cell referencing material id 0 with a declaration lacking key 0
    fails loud at construction (parse, don't validate downstream) — the
    ``Materials.restrict`` guard, witnessed on a genuine one-cell carrier."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, MirrorDeck())).partition(
        CellsByCount.uniform_width(1),
    ).mesh
    with pytest.raises(ValueError, match="references material ids"):
        MaterialMesh(mesh, {3: _mix([1.0, 1.0], ng=2)})


def test_the_transport_mesh_package_re_exports_no_mesh_name() -> None:
    """``orpheus.transport.mesh`` (the material mesh) re-exports none of the mesh names that moved to :mod:`orpheus.mesh`: a re-export
    would let a consumer keep the old spelling and hide the layer from the linter. Moved 2026-10-26 from
    ``tests/gates/mesh/test_module_layout.py::test_no_old_package_re_exports_a_mesh_name[orpheus.transport.mesh]`` (layer DAG ruling 27:
    as an import statement it was a mesh gate importing transport)."""
    import orpheus.transport.mesh as transport_mesh

    names = ("AxisMesh", "RadialAxisMesh", "Axis1D", "AxisCoord", "FaceLabel", "AXIS_NAMES")
    leaked = [name for name in names if hasattr(transport_mesh, name)]
    assert not leaked, f"orpheus.transport.mesh still exports {leaked}; import them from orpheus.mesh"


# ── Moved from tests/gates/numerics/test_space_of_axes.py by layer DAG Phase 4 (2026-11-01): the rows whose subject
# is the material mesh's ``bulk_space`` mint (a transport claim; the numerics file keeps ``of_axes``'s laws) ──


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def test_bulk_space_on_a_MESHED_carrier_is_the_honest_scalar_bulk() -> None:
    """B10 — ``(ng, *spatial)`` with cell-volume weights; and DISTINCT
    from the method mesh's angular composite.

    ``bulk_space`` is inherited by ``SNProblem``/``DiffusionMesh``, so the
    uniform formula must be honest on EVERY member (the seed of CS2's
    single scalar-bulk mint). The second half (the scalar bulk differs from
    the diffusion composite ``bulk (+) scalar trace``, built here from the
    raw constructors since layer DAG Phase 4) is what makes D7's
    chain-ordering claim non-vacuous.

    # CS1.5 re-point: ``bulk_space`` moves to Medium; this gate re-points
    # with it.
    """
    from orpheus.data.macro_xs.synthetic import get_mixture
    from tests.gates.transport._phase_space_cases import scalar_full_field_space
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import CellEdges, Mesher
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    mix = get_mixture("A", "2g")
    # Non-uniform edges, so the volumes are NOT all-ones and the weights
    # survive canonicalization (a uniform unit mesh would collapse to the
    # counting spelling and this gate would assert nothing).
    mesh1d = Mesher(StructuredGeometry.from_homogeneous(6.0, MirrorDeck())).partition(
        CellEdges(np.array([0.0, 1.0, 3.0, 6.0])),
    ).mesh
    carrier = MaterialMesh(mesh1d, {0: mix})
    space = carrier.bulk_space
    _require(space.shape == (2, 3), f"scalar bulk shape {space.shape} != (2, 3)")
    axes = space.axes
    assert axes is not None
    spatial_weights = axes[1].weights
    assert spatial_weights is not None
    _require(
        bool(np.array_equal(spatial_weights, carrier.volumes)),
        "the spatial factor measure must BE the cell volumes",
    )

    composite = scalar_full_field_space(carrier)
    _require(
        carrier.bulk_space != composite,
        "the scalar bulk and the composite carrier collapsed — D7's "
        "chain-ordering claim would be vacuous",
    )


def test_the_spatial_axis_is_minted_through_the_carriers_own_measure() -> None:
    """G6a (CS5) — the mesh generates the spatial measure; the measure
    mints the axis; the axis is TODAY'S axis (same identity, same digest)
    plus its generator.

    Independent literal anchor (verification plan Q4): the volumes and
    centres asserted below are HAND-DERIVED from the edge list — not read
    back from the mesh — so the mesh→measure→axis chain has a pin that
    is structurally independent of every array it threads.
    """
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.mesh import CellEdges, Mesher
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    mix = get_mixture("A", "2g")
    # Cartesian slab, edges 0|1|3|6 ⟹ BY HAND: widths (=volumes) 1,2,3;
    # centres 0.5, 2.0, 4.5. Non-uniform, so canonicalization keeps them.
    carrier = MaterialMesh(
        Mesher(StructuredGeometry.from_homogeneous(6.0, MirrorDeck())).partition(
            CellEdges(np.array([0.0, 1.0, 3.0, 6.0])),
        ).mesh,
        {0: mix},
    )
    space = carrier.bulk_space
    sp_ax = space.axis("spatial")
    g = sp_ax.origin.measure  # the rank-1 spatial axis carries its measure
    npt = np.testing
    npt.assert_array_equal(g.weights, np.array([1.0, 2.0, 3.0]),
                           err_msg="volumes, hand-derived")
    npt.assert_array_equal(g.nodes, np.array([0.5, 2.0, 4.5]),
                           err_msg="centres, hand-derived")
    # The atoms are content (S7, ruling 71): the axis keeps the centres and
    # the line they lie on, so the volumes declared with no atoms are another
    # axis (until S7 the two were one).
    npt.assert_array_equal(sp_ax.nodes, np.array([0.5, 2.0, 4.5]), err_msg="centres, hand-derived")
    _require(sp_ax.support == g.support, "the support is the measure's")
    literal = Axis(
        "spatial", carrier.spatial_shape, weights=carrier.volumes,
        origin=Declared(BasisKind.NODAL),
    )
    _require(not (sp_ax == literal), "an atomless literal must not pass for the sampled axis")
    # G8-spatial — the section law on the rank-1 arm (the law's domain):
    _require(Axis.of_measure(g, sp_ax.label) == sp_ax, "spatial mint: not a section")


def test_the_rank_d_spatial_axis_is_generator_less_BY_CONTRACT() -> None:
    """G6b (CS5) — the CS2 rank-d seam, pinned as a CONTRACT.

    A ``DiscreteMeasure`` is a flat atom list, so ``Axis.of_measure`` mints
    rank 1; a rank-d spatial axis (shape ``(nx, ny, nz)``) has no rank-d
    measure→axis pairing yet (posed since layer DAG Phase 4 on a 2-D
    ``MaterialMesh``, d = 2 being rank-d; until then a d = 3 ``SNProblem``,
    which a transport gate may not build). [M] minting it flat would change every
    d≥2 space name (``spatial(12,)#3712…`` vs ``spatial(3, 4)#1dcb…`` —
    verification plan R1), so the rank-d arm deliberately stays
    generator-less. ⛔ The day CS2 mints the rank-d pairing, THIS row
    must be inverted DELIBERATELY — it is the seam's witness, not a
    permanent truth.
    """
    from orpheus.mesh import Mesh2D
    from tests.gates.transport._materials import placeholder_materials

    carrier = MaterialMesh(
        Mesh2D(
            edges_x=np.linspace(0.0, 1.0, 3), edges_y=np.linspace(0.0, 2.0, 4), mat_map=np.zeros((2, 3), dtype=int),
            face_roles={face: MirrorDeck() for face in ("xmin", "xmax", "ymin", "ymax")},
        ),
        placeholder_materials(ng=2),
    )
    sp_ax = carrier.bulk_space.axis("spatial")
    _require(sp_ax.shape == (2, 3), f"fixture rank moved: {sp_ax.shape}")
    _require(
        isinstance(sp_ax.origin, Declared) and sp_ax.nodes is None,
        "the rank-d arm must stay Declared until CS2 mints the "
        "rank-d measure→axis pairing (inverting this is a deliberate act)",
    )


def test_bulk_space_is_cached_and_content_stable() -> None:
    """B11 — one carrier, ONE instance (``is``); equal carriers, EQUAL
    spaces (the derived name is content, not identity).

    # CS1.5 re-point: ``bulk_space`` moves to Medium; this gate re-points
    # with it.
    """
    from orpheus.data.macro_xs.synthetic import get_mixture
    from tests.gates.transport._carrier_helpers import unit_cell_carrier

    mix = get_mixture("A", "2g")
    carrier = unit_cell_carrier({0: mix})
    _require(carrier.bulk_space is carrier.bulk_space, "must be cached")
    twin = unit_cell_carrier({0: mix})
    _require(
        carrier.bulk_space == twin.bulk_space
        and hash(carrier.bulk_space) == hash(twin.bulk_space),
        "equal carriers must mint equal spaces",
    )
