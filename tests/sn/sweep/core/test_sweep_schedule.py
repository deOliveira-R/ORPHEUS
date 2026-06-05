r"""Phase 3 sub-step 3a: the polymorphic :class:`SweepSchedule` (Jacobi / G-S).

Foundation tests for the mesh-time octant-group schedule the SI resolvent
consumes (its uniform sweep-and-reflect loop). The schedule is the ONLY thing
that differs between Jacobi and Gauss-Seidel, so these pin its structure:

* Jacobi = ONE group, ALL octants, NO inter-group reflect.
* G-S = one group per in-plane octant (in quadrature sweep order); each group
  reflects the reflective faces its octants OUTFLOW through. Non-reflective
  outgoing faces (vacuum / white) and grazing axes contribute no reflect.
* INVARIANT: both schedules sweep the SAME ordinates (the schedule changes the
  splitting, never which ordinates are swept) — the fixed point is unchanged.

vv L11 discriminating: the slab cases hardcode the face↔direction wiring from
first principles (a −x-streaming octant exits through ``xmin``), catching a
sign/face swap in the schedule's outgoing-face map.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    Mesh1D,
    Mesh2D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep_graph import OctantLabel
from orpheus.sn.sweep_schedule import SweepSchedule
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _slab(bcs: tuple, nx: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=bcs,
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _box(*, nx: int = 4, ny: int = 4, ng: int = 1, **bc_kwargs) -> SNMesh:
    edges = np.linspace(0.0, 2.0, nx + 1)
    mesh = Mesh2D(
        edges_x=edges,
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        **bc_kwargs,
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _all_indices(schedule: SweepSchedule) -> list[int]:
    return sorted(
        idx
        for group in schedule.groups
        for sweep in group.sweeps
        for idx in sweep.indices
    )


def _reflective_faces(sn: SNMesh) -> set[str]:
    return {
        f for f in sn.trace.layout.faces if getattr(sn, f"bc_{f}") == "reflective"
    }


def _expected_outgoing(label: OctantLabel) -> set[str]:
    """First-principles outgoing faces for an octant label (independent of the
    schedule module's internal map — vv L11 structural independence)."""
    faces: set[str] = set()
    if label.sign_x > 0:
        faces.add("xmax")
    elif label.sign_x < 0:
        faces.add("xmin")
    if label.sign_y > 0:
        faces.add("ymax")
    elif label.sign_y < 0:
        faces.add("ymin")
    return faces


# ─── Jacobi structure ────────────────────────────────────────────────

def test_jacobi_is_one_group_no_reflect_slab():
    sn = _slab((BC.reflective, BC.reflective))
    sched = SweepSchedule.jacobi(sn)
    assert sched.kind == "jacobi"
    assert len(sched.groups) == 1
    assert sched.groups[0].reflect_faces == ()
    assert _all_indices(sched) == list(range(sn.quad.N))


def test_jacobi_is_one_group_no_reflect_box():
    sn = _box(
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    sched = SweepSchedule.jacobi(sn)
    assert len(sched.groups) == 1
    assert sched.groups[0].reflect_faces == ()
    assert _all_indices(sched) == list(range(sn.quad.N))


# ─── G-S slab (single-axis 2-cycle) ──────────────────────────────────

def test_gs_slab_reflective_specular_faces():
    """Discriminating (first principles): two groups in lexicographic order;
    the −x octant exits through ``xmin``, the +x octant through ``xmax``."""
    sn = _slab((BC.reflective, BC.reflective))
    sched = SweepSchedule.gauss_seidel(sn)
    assert sched.kind == "gauss_seidel"
    # one octant per group; lexicographic order puts (-1,0) before (+1,0).
    labels = [g.sweeps[0].label for g in sched.groups]
    assert labels == [OctantLabel(-1, 0), OctantLabel(+1, 0)]
    reflect = {g.sweeps[0].label: g.reflect_faces for g in sched.groups}
    assert reflect[OctantLabel(-1, 0)] == ("xmin",)
    assert reflect[OctantLabel(+1, 0)] == ("xmax",)


def test_gs_slab_vacuum_reflective_only_reflective_face_reflects():
    sn = _slab((BC.vacuum, BC.reflective))
    refl = _reflective_faces(sn)
    sched = SweepSchedule.gauss_seidel(sn)
    for g in sched.groups:
        label = g.sweeps[0].label
        expected = tuple(_expected_outgoing(label) & refl)
        assert set(g.reflect_faces) == set(expected)
    # exactly one octant reflects (the one exiting through the lone reflective
    # face), the other reflects nothing (it exits through vacuum).
    n_reflecting = sum(1 for g in sched.groups if g.reflect_faces)
    assert n_reflecting == 1


def test_gs_slab_full_vacuum_no_reflect_anywhere():
    sn = _slab((BC.vacuum, BC.vacuum))
    sched = SweepSchedule.gauss_seidel(sn)
    assert all(g.reflect_faces == () for g in sched.groups)


# ─── G-S 2-D box (per-axis 2-cycles) ─────────────────────────────────

def test_gs_box_all_reflective_each_octant_reflects_its_outgoing_faces():
    sn = _box(
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    sched = SweepSchedule.gauss_seidel(sn)
    labels = [g.sweeps[0].label for g in sched.groups]
    # one group per DISTINCT in-plane octant label (out-of-plane sign_z merged
    # into a single group, so a face's outflow is complete before it reflects).
    assert len(set(labels)) == len(labels)
    for g in sched.groups:
        # every sweep in a group shares the group's in-plane label.
        assert all(s.label == g.sweeps[0].label for s in g.sweeps)
        # all four faces reflective ⇒ every outgoing face reflects (1 face for
        # an axis-aligned octant, 2 for a diagonal — quadrature-dependent).
        assert set(g.reflect_faces) == _expected_outgoing(g.sweeps[0].label)
    # collectively, every reflective box face is reflected by the octant(s) that
    # outflow through it — the full reflective coupling is scheduled, none dropped.
    reflected = {f for g in sched.groups for f in g.reflect_faces}
    assert reflected == {"xmin", "xmax", "ymin", "ymax"}


def test_gs_box_half_reflective_no_reflect_on_vacuum_axis():
    # x reflective (both faces), y vacuum (both faces): only x faces reflect.
    sn = _box(
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.vacuum, bc_ymax=BC.vacuum,
    )
    sched = SweepSchedule.gauss_seidel(sn)
    for g in sched.groups:
        assert "ymin" not in g.reflect_faces
        assert "ymax" not in g.reflect_faces
        label = g.sweeps[0].label
        if label.sign_x > 0:
            assert g.reflect_faces == ("xmax",)
        elif label.sign_x < 0:
            assert g.reflect_faces == ("xmin",)


# ─── the fixed-point invariant: same ordinates, different splitting ──

def test_jacobi_and_gs_sweep_identical_ordinates():
    sn = _box(
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    full = list(range(sn.quad.N))
    assert _all_indices(SweepSchedule.jacobi(sn)) == full
    assert _all_indices(SweepSchedule.gauss_seidel(sn)) == full
