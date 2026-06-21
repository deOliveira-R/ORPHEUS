r"""C3.6 synthetic d≠2 schedule pins — the honest in-plane projection.

The schedule (:func:`~orpheus.sn.sweep_schedule._octant_sweep`) is the SOLE
in-plane projection site: it keeps the mesh's first ``ndim`` quadrature
signs and zero-pads missing axes.  Before C3.6 it truncated EVERY label to a
2-tuple and ``_outgoing_faces`` hand-listed x/y faces — at d=3 a +z octant
would NEVER shed ``zmax``, so a reflective z face got NO Gauss-Seidel
reflect group and the G-S splitting converged to a silently WRONG fixed
point (the ERR-056 class: shed-EXISTENCE here, where the original ERR-056
was shed-ORDER on a shared face).  These pins close that hole BY
CONSTRUCTION before a 3-axis mesh exists.

Scope boundary (vv Mode 7 — declare what the ansatz nulls): the duck-typed
pins prove the schedule's STRUCTURE (labels, face sets, reflect-group
assignment) cheaply at the four reads the schedule makes. Since C5.5
(#225) a d=3 ``SNMesh`` IS constructible (mesh-less ``from_axes``):
``test_gs_d3_schedule_from_real_mesh`` below pins the order-INVARIANT
facts live, and the d=3 flux-VALUE claim lands in
``tests/sn/solve/test_d3_admission.py`` (the Mode-9 G-S≡Jacobi box).
The schedule layer is mesh-LIGHT by design: ``_octant_sweep`` /
``_outgoing_faces`` are pure on quadrature partition entries / labels, and
``SweepSchedule.gauss_seidel`` reads only ``ndim``, ``quad.octants``,
``trace.layout.faces``, and the ``bc_{face}`` attributes — all forgeable
with :class:`types.SimpleNamespace` (the duck-types live HERE, not as a
production injection seam — one consumer, Pattern 6 defers the abstraction).
"""
from __future__ import annotations

import itertools
from types import SimpleNamespace

import pytest

from orpheus.sn.sweep_graph import OctantLabel
from orpheus.sn.sweep_schedule import (
    SweepSchedule,
    _octant_sweep,
    _outgoing_faces,
)

from tests.sn.sweep.core.test_sweep_schedule import _expected_outgoing

pytestmark = [pytest.mark.foundation]


def _entry(signs: tuple[int, ...], indices: tuple[int, ...] = (0,)):
    """Forge a quadrature octant-partition entry (label + indices is all the
    schedule reads; the partition's measure is untouched)."""
    return SimpleNamespace(label=signs, indices=indices)


# ─── Tier 1: the pure projection helpers, no mesh ────────────────────

def test_octant_sweep_keeps_all_ndim_signs_at_d3():
    """The in-plane projection keeps the mesh's first ndim signs — at d=3
    nothing is dropped (pre-C3.6 this truncated to ``(+1, -1)``)."""
    sweep = _octant_sweep(_entry((+1, -1, +1)), ndim=3)
    assert sweep.label == OctantLabel((+1, -1, +1))


def test_octant_sweep_projects_out_of_plane_signs_at_d2():
    """A 2-D mesh under an S² cubature drops sign_z — the pre-C3.6
    behaviour, preserved exactly (the d=2 bit-identity half of the pin)."""
    sweep = _octant_sweep(_entry((+1, -1, +1)), ndim=2)
    assert sweep.label == OctantLabel((+1, -1))


def test_octant_sweep_zero_pads_short_labels():
    """A quadrature with fewer signs than the mesh has axes (a slab
    quadrature over a multi-D mesh) zero-pads: sign 0 = no streaming on
    that axis."""
    sweep = _octant_sweep(_entry((+1,)), ndim=2)
    assert sweep.label == OctantLabel((+1, 0))


@pytest.mark.catches("ERR-056")
def test_outgoing_faces_include_z_at_d3():
    """A +z-streaming octant SHEDS zmax — the shed-EXISTENCE hole the
    hand-listed x/y face table left open at d=3."""
    assert _outgoing_faces(OctantLabel((+1, -1, +1))) == (
        "xmax", "ymin", "zmax",
    )
    assert _outgoing_faces(OctantLabel((-1, 0, -1))) == ("xmin", "zmin")
    # d=1 honest 1-tuple labels derive a 1-face set.
    assert _outgoing_faces(OctantLabel((-1,))) == ("xmin",)


def test_outgoing_faces_match_first_principles_table_all_octants():
    """Every d ∈ {1,2,3} sign combination against the hand-transcribed
    (axis, sign) → face table (vv L11: production derives, the test
    hand-lists — a x↔y↔z swap in either is observable)."""
    for ndim in (1, 2, 3):
        for signs in itertools.product((-1, 0, +1), repeat=ndim):
            label = OctantLabel(signs)
            assert set(_outgoing_faces(label)) == _expected_outgoing(label)


# ─── Tier 2: the real Gauss-Seidel assignment loop at d=3 ────────────

_D3_FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")


def _fake_mesh_3d(**bcs: str) -> SimpleNamespace:
    """Duck-typed 3-axis mesh: exactly the four reads the schedule makes.

    ``bcs`` is keyed by plain face name (``xmin="reflective"``); unnamed
    faces default to vacuum.  Post-C4 (#220) the consumer contract is the
    face-name-keyed ``bc`` dict, mirrored here.  The 8 sign-octants stand in for an S²
    cubature's octant partition (one ordinate each, lexicographic order —
    the assignment assertions are order-relative, matching how
    ``gauss_seidel`` consumes the entry order).
    """
    octants = tuple(
        _entry(signs, indices=(i,))
        for i, signs in enumerate(itertools.product((-1, +1), repeat=3))
    )
    return SimpleNamespace(
        ndim=3,
        quad=SimpleNamespace(octants=octants),
        trace=SimpleNamespace(layout=SimpleNamespace(faces=_D3_FACES)),
        bc={face: bcs.get(face, "vacuum") for face in _D3_FACES},
    )


@pytest.mark.catches("ERR-056")
def test_gs_d3_mixed_bc_reflects_only_reflective_axes():
    """x reflective / y vacuum / z reflective: every group reflects exactly
    its outgoing ∩ reflective faces — the MIXED-BC config makes an
    axis-name swap observable (an all-reflective box would mask it: the
    (+1,+1,+1) octant must reflect xmax and zmax but NEVER ymax)."""
    fake = _fake_mesh_3d(
        xmin="reflective", xmax="reflective",
        zmin="reflective", zmax="reflective",
    )
    sched = SweepSchedule.gauss_seidel(fake)
    assert sched.kind == "gauss_seidel"
    assert len(sched.groups) == 8  # one per sign-octant, none merged
    reflected = {f for g in sched.groups for f in g.reflect_faces}
    assert reflected == {"xmin", "xmax", "zmin", "zmax"}
    for g in sched.groups:
        label = g.sweeps[0].label
        assert label.ndim == 3
        # a group only ever reflects faces it actually outflows through.
        assert set(g.reflect_faces) <= _expected_outgoing(label)
        assert not {"ymin", "ymax"} & set(g.reflect_faces)


@pytest.mark.catches("ERR-056")
def test_gs_d3_every_reflective_face_assigned_to_last_outflowing_group():
    """The shared-face discipline at d=3: each reflective face (shared by
    4 of the 8 octants) is assigned to EXACTLY ONE group — the LAST in
    sweep order that outflows it (premature reflection would absorb
    not-yet-swept octants' seed values: the original ERR-056)."""
    fake = _fake_mesh_3d(**{face: "reflective" for face in _D3_FACES})
    sched = SweepSchedule.gauss_seidel(fake)
    counts: dict[str, int] = {}
    for g in sched.groups:
        for f in g.reflect_faces:
            counts[f] = counts.get(f, 0) + 1
    assert set(counts) == set(_D3_FACES)
    assert all(c == 1 for c in counts.values()), counts
    order = {g.sweeps[0].label: i for i, g in enumerate(sched.groups)}
    for face in _D3_FACES:
        outflowing = [
            g.sweeps[0].label for g in sched.groups
            if face in _expected_outgoing(g.sweeps[0].label)
        ]
        assert len(outflowing) == 4  # shared by every same-sign octant
        assigned = next(
            g.sweeps[0].label for g in sched.groups
            if face in g.reflect_faces
        )
        assert order[assigned] == max(order[lab] for lab in outflowing)


def test_jacobi_d3_one_group_all_octants():
    """Jacobi at d=3: ONE group, all 8 octants, no reflect — the schedule
    split is dimension-blind."""
    fake = _fake_mesh_3d()
    sched = SweepSchedule.jacobi(fake)
    assert len(sched.groups) == 1
    assert sched.groups[0].reflect_faces == ()
    labels = {s.label for s in sched.groups[0].sweeps}
    assert labels == {
        OctantLabel(signs)
        for signs in itertools.product((-1, +1), repeat=3)
    }


# ─── C5.5 (#225): the real mesh presents the duck-typed contract ─────


@pytest.mark.catches("ERR-056")
def test_gs_d3_schedule_from_real_mesh():
    """A real ``from_axes`` 3-axis mesh feeds ``gauss_seidel`` the same
    interface the synthetic pins above assumed — and the order-invariant
    schedule facts hold on it: 8 unmerged sign-octant groups, the
    reflected set is exactly the reflective faces, each shared face is
    assigned to EXACTLY ONE group (the ERR-056 discipline), and a group
    only reflects faces it outflows. (Order-RELATIVE facts — which
    group is "last" — depend on the cubature's octant enumeration and
    are pinned synthetically above.)
    """
    import numpy as np
    from orpheus.derivations.common.xs_library import make_mixture
    from orpheus.geometry import BC
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.axis import AxisMesh
    from orpheus.sn.geometry import SNMesh

    refl, vac = BC("reflective"), BC("vacuum")
    mix = make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([0.5]),
        sig_f=np.array([0.0]), nu=np.array([0.0]),
        chi=np.zeros(1), sig_s=np.array([[0.5]]),  # non-fissile ⇒ null χ (S10a guard)
    )
    mesh = SNMesh.from_axes(
        (
            AxisMesh(edges=np.linspace(0.0, 1.0, 3), bc_low=refl, bc_high=refl),
            AxisMesh(edges=np.linspace(0.0, 1.0, 4), bc_low=vac, bc_high=vac),
            AxisMesh(edges=np.linspace(0.0, 1.0, 5), bc_low=refl, bc_high=refl),
        ),
        Quadrature.level_symmetric(sn_order=4),
        {0: mix},
    )
    sched = SweepSchedule.gauss_seidel(mesh)
    assert sched.kind == "gauss_seidel"
    assert len(sched.groups) == 8
    counts: dict[str, int] = {}
    for g in sched.groups:
        label = g.sweeps[0].label
        assert label.ndim == 3
        assert set(g.reflect_faces) <= _expected_outgoing(label)
        for f in g.reflect_faces:
            counts[f] = counts.get(f, 0) + 1
    assert set(counts) == {"xmin", "xmax", "zmin", "zmax"}
    assert all(c == 1 for c in counts.values()), counts
