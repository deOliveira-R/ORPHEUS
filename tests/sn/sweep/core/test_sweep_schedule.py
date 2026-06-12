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


def _box(*, nx: int = 4, ny: int = 4, ng: int = 1, quad=None, **bc_kwargs) -> SNMesh:
    edges = np.linspace(0.0, 2.0, nx + 1)
    mesh = Mesh2D(
        edges_x=edges,
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        **bc_kwargs,
    )
    if quad is None:
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
        f for f in sn.trace.layout.faces if sn.bc[f] == "reflective"
    }


# Hand-transcribed (axis, sign) → face table — deliberately NOT derived from
# the production ``AXIS_NAMES`` f-string (mirror-not-import: production
# derives, the test hand-lists, so a drift in either is observable —
# vv L11 structural independence).
_FACE_OF = {
    (0, +1): "xmax", (0, -1): "xmin",
    (1, +1): "ymax", (1, -1): "ymin",
    (2, +1): "zmax", (2, -1): "zmin",
}


def _expected_outgoing(label: OctantLabel) -> set[str]:
    """First-principles outgoing faces for an octant label (independent of the
    schedule module's derivation — vv L11 structural independence)."""
    return {_FACE_OF[(a, s)] for a, s in enumerate(label.signs) if s != 0}


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
    the −x octant exits through ``xmin``, the +x octant through ``xmax``.

    HISTORY (C3.6): the d=1 labels are honest 1-tuples ``(±1,)`` — the
    schedule projects to the mesh's ``ndim``.  Before C3.6 they carried a
    phantom zero-padded second sign ``(±1, 0)`` that the walk re-truncated.
    """
    sn = _slab((BC.reflective, BC.reflective))
    sched = SweepSchedule.gauss_seidel(sn)
    assert sched.kind == "gauss_seidel"
    # one octant per group; lexicographic order puts (-1,) before (+1,).
    labels = [g.sweeps[0].label for g in sched.groups]
    assert labels == [OctantLabel((-1,)), OctantLabel((+1,))]
    reflect = {g.sweeps[0].label: g.reflect_faces for g in sched.groups}
    assert reflect[OctantLabel((-1,))] == ("xmin",)
    assert reflect[OctantLabel((+1,))] == ("xmax",)


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
        if label.signs[0] > 0:
            assert g.reflect_faces == ("xmax",)
        elif label.signs[0] < 0:
            assert g.reflect_faces == ("xmin",)


# ─── diagonal cubature: the shared-face correctness pin (regression) ──

@pytest.mark.catches("ERR-056")
def test_gs_diagonal_quadrature_shared_face_assigned_to_last_group_only():
    r"""⚠ Correctness pin (the Lebedev shared-face bug, 3c.2 / ERR-056).

    On a DIAGONAL / spherical cubature (``lebedev``) each octant outflows TWO
    faces and each reflective face is shared by ≥2 in-plane octant groups
    (e.g. ``xmax`` ← every +x octant: ``(+1,0)``, ``(+1,+1)``, ``(+1,−1)``).
    A reflective face MUST be assigned to EXACTLY ONE group — the LAST (in
    sweep order) that outflows it — so the inter-group reflect fires only when
    the face's outflow is COMPLETE.  Assigning to a non-last (premature) group
    would absorb the not-yet-swept octants' SEED value (the wavefront is rebuilt
    + seeded each solve) and reflect garbage → a WRONG fixed point (this drove
    the lebedev all-reflective flux to 3.4 instead of the flat 5.88)."""
    quad = Quadrature.lebedev(order=17)
    sn = _box(
        quad=quad,
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    sched = SweepSchedule.gauss_seidel(sn)
    order = {g.sweeps[0].label: i for i, g in enumerate(sched.groups)}
    # every reflective face appears in EXACTLY ONE group (no double-assignment).
    counts: dict[str, int] = {}
    for g in sched.groups:
        for f in g.reflect_faces:
            counts[f] = counts.get(f, 0) + 1
    assert set(counts) == {"xmin", "xmax", "ymin", "ymax"}
    assert all(c == 1 for c in counts.values()), f"shared face double-assigned: {counts}"
    # the assigned group must OUTFLOW the face AND be the LAST one to do so.
    for f in counts:
        outflowing = [
            g.sweeps[0].label for g in sched.groups
            if f in _expected_outgoing(g.sweeps[0].label)
        ]
        assigned = next(
            g.sweeps[0].label for g in sched.groups if f in g.reflect_faces
        )
        assert f in _expected_outgoing(assigned)
        assert order[assigned] == max(order[label] for label in outflowing), (
            f"{f} not assigned to the LAST group outflowing it"
        )


# ─── the fixed-point invariant: same ordinates, different splitting ──

def test_jacobi_and_gs_sweep_identical_ordinates():
    sn = _box(
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    full = list(range(sn.quad.N))
    assert _all_indices(SweepSchedule.jacobi(sn)) == full
    assert _all_indices(SweepSchedule.gauss_seidel(sn)) == full
