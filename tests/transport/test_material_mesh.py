r"""MaterialMesh — the method-agnostic mesh + materials data carrier.

:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` is the
"mesh + materials" middle type between geometry ``Mesh1D`` (material
*ids*, no cross sections) and ``SNMesh`` (mesh + materials + quadrature +
sweep machinery).  These tests pin its **intrinsic data contract** (the
``coding-elegance`` standard: a type ships a test of its defining
invariants) and the **data/behavior split** with ``SNMesh(MaterialMesh)``:

* the carrier holds mesh + materials and exposes the method-agnostic data
  accessors (``ng`` / ``volumes`` / ``volume_measure`` / ``ndim`` /
  ``spatial_shape`` / ``mat_map`` / ``material_xs_field``);
* ``ng`` consistency is enforced at construction;
* ``SNMesh`` **is-a** ``MaterialMesh`` (Liskov) — every carrier accessor
  works on an ``SNMesh``, and its data block is bit-identical to a
  standalone ``MaterialMesh`` built from the same inputs;
* ``SNMesh.from_material_mesh`` promotes a carrier to a solvable phase
  space (the homogenization re-solve path).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.mesh import (
    InconsistentMaterialsError,
    MaterialMesh,
    MaterialXSField,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh

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
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, 6),
        mat_ids=np.array([0, 0, 1, 1, 0]),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
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
    field = mm.material_xs_field()
    assert isinstance(field, MaterialXSField)
    # per-cell SigT view follows the mat_map (cell 2 is material 1).
    sig_t = field.total_cross_section  # (ng, nx)
    assert sig_t.shape == (2, 5)
    np.testing.assert_allclose(sig_t[:, 2], mats[1].SigT)
    np.testing.assert_allclose(sig_t[:, 0], mats[0].SigT)


# ── ng-consistency invariant (the defining law) ───────────────────────

def test_inconsistent_ng_raises_at_construction():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 3), mat_ids=np.array([0, 1]),
        coord=CoordSystem.CARTESIAN,
    )
    mats = {0: _mix([1.0, 1.5], ng=2), 1: _mix([1.0], ng=1)}
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        MaterialMesh(mesh, mats)


def test_missing_material_id_raises():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 3), mat_ids=np.array([0, 7]),
        coord=CoordSystem.CARTESIAN,
    )
    with pytest.raises(ValueError, match="references material ids"):
        MaterialMesh(mesh, {0: _mix([1.0, 1.0], ng=2)})


# ── Liskov: SNMesh IS-A MaterialMesh ──────────────────────────────────

def test_snmesh_is_a_material_mesh():
    assert issubclass(SNMesh, MaterialMesh)


def test_snmesh_data_block_bit_identical_to_standalone_carrier():
    """An SNMesh's inherited data block matches a standalone MaterialMesh
    built from the same mesh + materials (the split is bit-identical)."""
    mesh, mats = _two_material_mesh()
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    snm = SNMesh(mesh, quad, mats)
    mm = MaterialMesh(mesh, mats)
    assert snm.ng == mm.ng
    assert snm.ndim == mm.ndim
    assert snm.spatial_shape == mm.spatial_shape
    np.testing.assert_array_equal(snm.mat_map, mm.mat_map)
    np.testing.assert_array_equal(snm.volumes, mm.volumes)
    np.testing.assert_array_equal(
        snm.volume_measure.weights, mm.volume_measure.weights,
    )
    # Every carrier accessor is callable on the SNMesh (substitutability).
    assert isinstance(snm.material_xs_field(), MaterialXSField)


# ── from_material_mesh promotion (the data/behavior join) ─────────────

def test_from_material_mesh_round_trips_to_solvable_snmesh():
    mesh, mats = _two_material_mesh()
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mm = MaterialMesh(mesh, mats)
    snm = SNMesh.from_material_mesh(mm, quad)
    assert isinstance(snm, SNMesh)
    # carries the carrier's data verbatim …
    np.testing.assert_array_equal(snm.mat_map, mm.mat_map)
    np.testing.assert_array_equal(snm.volumes, mm.volumes)
    assert snm.materials is mm.materials
    # … plus the SN method layer (quadrature + streaming stencil).
    assert snm.quad is quad
    assert snm.reduced is not None


# ── Meshless single-cell carrier (from_materials) ─────────────────────
# The degenerate carrier an infinite homogeneous medium lives on: the
# phase space the k∞ solver assembles (C − K_iso − F/k) on with streaming
# dropped.  These pin its defining invariants (the coding-elegance
# standard: a type ships a test of its intrinsic contract).


def test_from_materials_is_meshless_single_cell():
    """One Cartesian cell, one region (id 0), no legacy mesh adapter."""
    mm = MaterialMesh.from_materials({0: _mix([1.0, 1.5], ng=2)})
    assert mm.mesh is None              # meshless
    assert mm.ndim == 1
    assert mm.spatial_shape == (1,)
    assert mm.nx == 1
    np.testing.assert_array_equal(mm.mat_map, [0])   # single region, id 0
    assert mm.ng == 2                   # derived from materials (no twin)


def test_from_materials_unit_volume_measure():
    """The lone cell has unit volume — volume_measure weights it 1.0.

    k∞ is a production/absorption ratio (volume-invariant); the unit
    choice keeps reaction rates equal to the bare ⟨Σ,φ⟩ contraction.
    """
    mu = MaterialMesh.from_materials({0: _mix([1.0, 1.5], ng=2)}).volume_measure
    np.testing.assert_array_equal(mu.weights, [1.0])
    assert mu.weights.shape == (1,)


def test_from_materials_resolves_xs_field():
    """material_xs_field() resolves on the meshless carrier and exposes the
    single material's cross sections in the lone cell."""
    mix = _mix([1.0, 1.5], ng=2)
    field = MaterialMesh.from_materials({0: mix}).material_xs_field()
    assert isinstance(field, MaterialXSField)
    sig_t = field.total_cross_section   # (ng, nx) = (2, 1)
    assert sig_t.shape == (2, 1)
    np.testing.assert_allclose(sig_t[:, 0], mix.SigT)


def test_from_materials_requires_material_id_zero():
    """The single cell carries material id 0, so a materials dict lacking
    key 0 fails loud at construction (parse, don't validate downstream)."""
    with pytest.raises(ValueError, match="references material ids"):
        MaterialMesh.from_materials({3: _mix([1.0, 1.0], ng=2)})
