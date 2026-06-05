r"""Foundation — :class:`WavefrontFlux` (interior face-flux cochain) Phase 1.

Pins the minted contract for
:class:`orpheus.transport.fields.wavefront_flux.WavefrontFlux` and its space
:class:`orpheus.numerics.spaces.interior_face_space.InteriorFaceSpace`
(Wave O #205/#208, plan ``wavefront_flux_foundation.md`` Phase 1):

* **Units / class identity** — ``UNITS == ANGULAR_FLUX_UNITS``; same-class
  algebra works; cross-class arithmetic (vs ``BoundaryFlux``) raises.
* **Field + zero-copy views** — ``face(axis)`` shares memory with ``values``;
  writes propagate; per-axis shapes ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)``.
* **Biproduct laws** (``C¹ = C¹_int ⊕ C¹_∂``): ``ι* ∘ ι_* = id`` on the
  boundary chain ("absorption = identity"; project-after-inject = seed then
  absorb), and ``π_int ∘ ι_∂ = 0`` (the injection leaves the
  strictly-interior faces untouched).
* **The typed ι_*/ι* round-trip pin** (the cross-domain-attacker's Frame-1
  first test) — bit-for-bit vs the raw seed/absorb, with a NO-transpose
  discriminator and a negative control so the pin is not vacuous (L11).
* **Axis-parametric API** (§3a′ 3D-readiness) — the interior ``FaceLayout``
  builds for ``axes=(0,)``, ``(0,1)``, ``(0,1,2)`` through ONE code path.
* **Factories / frozen / mesh-binding** guards.

All ``@pytest.mark.foundation`` (software-invariant; no theory ``:label:``).
Contracts that wrap the SAME buffers use ``np.array_equal`` (bit-for-bit), per
``vv-principles`` Bit-identity — this is a type-only mint.
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.interior_face_space import InteriorFaceSpace
from orpheus.numerics.units import ANGULAR_FLUX_UNITS
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.fields.wavefront_flux import WavefrontFlux

from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 6, ng: int = 2, bc: str = "reflective") -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC(bc), bc_right=BC(bc),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 6, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cart2d_mesh(nx: int = 6, ny: int = 5, ng: int = 2, bc: str = "reflective") -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC(bc), bc_xmax=BC(bc), bc_ymin=BC(bc), bc_ymax=BC(bc),
    )
    quad = Quadrature.level_symmetric(4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _rng():
    return np.random.default_rng(20260604)


# ───────────────────────────────────────────────────────────────────────
# Units + class identity
# ───────────────────────────────────────────────────────────────────────


def test_units_are_angular_flux():
    assert WavefrontFlux.UNITS == ANGULAR_FLUX_UNITS


def test_is_a_field():
    sn = _cart2d_mesh()
    assert isinstance(WavefrontFlux.zeros_on(sn), Field)


def test_same_class_algebra_works():
    sn = _cart2d_mesh()
    a = WavefrontFlux.zeros_on(sn)
    a.values[:] = 1.0
    b = WavefrontFlux.zeros_on(sn)
    b.values[:] = 2.0
    assert np.array_equal((a + b).values, np.full_like(a.values, 3.0))
    assert np.array_equal((3.0 * a).values, np.full_like(a.values, 3.0))


def test_cross_class_arithmetic_raises():
    sn = _cart2d_mesh()
    wf = WavefrontFlux.zeros_on(sn)
    bf = BoundaryFlux.zeros_on(sn)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = wf + bf  # type: ignore[operator]


# ───────────────────────────────────────────────────────────────────────
# Field + zero-copy views
# ───────────────────────────────────────────────────────────────────────


def test_face_shapes_2d():
    sn = _cart2d_mesh(nx=6, ny=5)
    wf = WavefrontFlux.zeros_on(sn)
    N, ng = sn.quad.N, sn.ng
    assert wf.face(0).shape == (N, ng, 7, 5)   # x-normal: nx+1
    assert wf.face(1).shape == (N, ng, 6, 6)   # y-normal: ny+1
    assert wf.axes == (0, 1)


def test_face_shapes_1d():
    sn = _slab_mesh(nx=8)
    wf = WavefrontFlux.zeros_on(sn)
    N, ng = sn.quad.N, sn.ng
    assert wf.face(0).shape == (N, ng, 9)      # x-normal only
    assert wf.axes == (0,)


def test_face_is_zero_copy_view():
    sn = _cart2d_mesh()
    wf = WavefrontFlux.zeros_on(sn)
    for a in wf.axes:
        assert np.shares_memory(wf.face(a), wf.values)
    # writes through the view propagate to the backing buffer
    wf.face(0)[0, 0, 0, 0] = 7.0
    assert wf.face(0)[0, 0, 0, 0] == 7.0
    assert (wf.values == 7.0).sum() == 1


def test_face_missing_axis_raises():
    sn = _slab_mesh()           # 1-D: only axis 0
    wf = WavefrontFlux.zeros_on(sn)
    with pytest.raises(KeyError):
        wf.face(1)


# ───────────────────────────────────────────────────────────────────────
# Factories
# ───────────────────────────────────────────────────────────────────────


def test_zeros_on_is_all_zero_and_sized():
    sn = _cart2d_mesh(nx=6, ny=5)
    wf = WavefrontFlux.zeros_on(sn)
    N, ng = sn.quad.N, sn.ng
    expected = N * ng * (7 * 5) + N * ng * (6 * 6)   # x-faces + y-faces flat
    assert wf.values.shape == (expected,)
    assert not wf.values.any()


def test_from_mesh_round_trips_buffer():
    sn = _cart2d_mesh()
    wf0 = WavefrontFlux.zeros_on(sn)
    buf = _rng().standard_normal(wf0.values.shape)
    wf = WavefrontFlux.from_mesh(buf, sn)
    assert wf.values is buf
    assert isinstance(wf.space, InteriorFaceSpace)


# ───────────────────────────────────────────────────────────────────────
# Biproduct laws  C¹ = C¹_int ⊕ C¹_∂
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("mesh_fn", [_cart2d_mesh, _slab_mesh, _sphere_mesh])
def test_absorption_is_identity(mesh_fn):
    r"""ι* ∘ ι_* = id on the boundary chain (seed then absorb, no walk)."""
    sn = mesh_fn()
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf)                       # ι_*
    out = BoundaryFlux.zeros_on(sn)
    wf.absorb(out)                    # ι*
    assert np.array_equal(out.values, bf.values)


def test_pi_int_after_injection_is_zero_2d():
    r"""π_int ∘ ι_∂ = 0 — injecting the boundary leaves the STRICTLY-interior
    faces (positions 1..n-1 along each axis) untouched at zero."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape) + 1.0  # nonzero
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf)
    # x-normal interior (drop the two domain edges i=0, i=nx) is all zero
    assert not wf.face(0)[:, :, 1:-1, :].any()
    # y-normal interior (drop j=0, j=ny) is all zero
    assert not wf.face(1)[:, :, :, 1:-1].any()


# ───────────────────────────────────────────────────────────────────────
# The typed ι_*/ι* round-trip pin (attacker Frame-1 first test)
# ───────────────────────────────────────────────────────────────────────


def test_seed_matches_raw_no_transpose_2d():
    r"""The seeded domain-edge slots equal the raw edge-slice assignment
    bit-for-bit — pins that the interior FaceLayout edge ordering AGREES with
    the BoundaryFlux face_view ordering (NO transpose; attacker risk #2)."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf)
    nx, ny = sn.nx, sn.ny
    assert np.array_equal(wf.face(0)[:, :, 0, :], bf.face_view("xmin"))
    assert np.array_equal(wf.face(0)[:, :, nx, :], bf.face_view("xmax"))
    assert np.array_equal(wf.face(1)[:, :, :, 0], bf.face_view("ymin"))
    assert np.array_equal(wf.face(1)[:, :, :, ny], bf.face_view("ymax"))


def test_edge_view_reads_seeded_trace_2d():
    r"""edge_view(face) is the zero-copy domain-edge slot — equals the seeded
    boundary face (the read counterpart of seed/absorb), via the SAME
    _edge_slot mapping (no duplicated face→slot literals at the call site)."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf)
    for face in ("xmin", "xmax", "ymin", "ymax"):
        assert np.array_equal(wf.edge_view(face), bf.face_view(face))
        assert np.shares_memory(wf.edge_view(face), wf.values)


def test_roundtrip_matches_raw_seed_absorb_2d():
    r"""The typed seed/absorb reproduces the raw numpy seed/absorb the legacy
    sweep does, bit-for-bit (both edges of the carve agree)."""
    sn = _cart2d_mesh(nx=6, ny=5)
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)

    # raw
    px = np.zeros((N, ng, nx + 1, ny)); py = np.zeros((N, ng, nx, ny + 1))
    px[:, :, 0, :] = bf.face_view("xmin"); px[:, :, nx, :] = bf.face_view("xmax")
    py[:, :, :, 0] = bf.face_view("ymin"); py[:, :, :, ny] = bf.face_view("ymax")
    raw_out = BoundaryFlux.zeros_on(sn)
    raw_out.face_view("xmin")[:] = px[:, :, 0, :]
    raw_out.face_view("xmax")[:] = px[:, :, nx, :]
    raw_out.face_view("ymin")[:] = py[:, :, :, 0]
    raw_out.face_view("ymax")[:] = py[:, :, :, ny]

    # typed
    wf = WavefrontFlux.zeros_on(sn); wf.seed(bf)
    typed_out = BoundaryFlux.zeros_on(sn); wf.absorb(typed_out)

    assert np.array_equal(typed_out.values, raw_out.values)


def test_negative_control_transposed_seed_breaks_identity():
    r"""L11 negative control: a deliberately transposed seed (xmin↔xmax) MUST
    break the round-trip identity — proves the round-trip pin is not vacuous."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf = WavefrontFlux.zeros_on(sn)
    # WRONG: swap xmin/xmax into each other's edge slots
    wf.face(0)[:, :, 0, :] = bf.face_view("xmax")
    wf.face(0)[:, :, sn.nx, :] = bf.face_view("xmin")
    wf.face(1)[:, :, :, 0] = bf.face_view("ymin")
    wf.face(1)[:, :, :, sn.ny] = bf.face_view("ymax")
    out = BoundaryFlux.zeros_on(sn); wf.absorb(out)
    assert not np.array_equal(out.values, bf.values)


# ───────────────────────────────────────────────────────────────────────
# Face-restricted seed / absorb (Phase 3 sub-step 3b — G-S per-group trace)
# ───────────────────────────────────────────────────────────────────────


def test_seed_faces_subset_touches_only_those_edges():
    r"""seed(bf, faces=("xmin",)) injects ONLY the xmin domain-edge slot; the
    other domain edges stay at their prior (zero) value — the exact-restriction
    the G-S re-seed relies on."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape) + 1.0  # nonzero
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf, faces=("xmin",))
    nx, ny = sn.nx, sn.ny
    assert np.array_equal(wf.face(0)[:, :, 0, :], bf.face_view("xmin"))
    # every other domain edge untouched (still zero)
    assert not wf.face(0)[:, :, nx, :].any()    # xmax
    assert not wf.face(1)[:, :, :, 0].any()     # ymin
    assert not wf.face(1)[:, :, :, ny].any()    # ymax


def test_absorb_faces_subset_writes_only_those_faces():
    r"""absorb(out, faces=("xmax",)) writes ONLY the xmax outflow back; the
    other faces of ``out`` stay zero — the per-group outflow writeback the G-S
    schedule does before reflecting that face."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf = WavefrontFlux.zeros_on(sn)
    wf.seed(bf)                                  # whole-trace seed
    out = BoundaryFlux.zeros_on(sn)
    wf.absorb(out, faces=("xmax",))
    assert np.array_equal(out.face_view("xmax"), bf.face_view("xmax"))
    for other in ("xmin", "ymin", "ymax"):
        assert not out.face_view(other).any()


def test_seed_faces_none_equals_explicit_all():
    r"""faces=None ≡ explicitly listing every boundary face (the default IS the
    whole trace)."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf_none = WavefrontFlux.zeros_on(sn); wf_none.seed(bf)
    wf_all = WavefrontFlux.zeros_on(sn)
    wf_all.seed(bf, faces=tuple(bf.layout.faces))
    assert np.array_equal(wf_none.values, wf_all.values)


def test_single_face_seeds_partition_whole_trace():
    r"""Seeding each face individually reconstructs the whole-trace seed — the
    faces partition the interior edge with no overlap (block-diagonal exact
    restriction; mirrors the SNBoundaryOperator single-face partition pin)."""
    sn = _cart2d_mesh(nx=6, ny=5)
    bf = BoundaryFlux.zeros_on(sn)
    bf.values[:] = _rng().standard_normal(bf.values.shape)
    wf_whole = WavefrontFlux.zeros_on(sn); wf_whole.seed(bf)
    wf_piece = WavefrontFlux.zeros_on(sn)
    for face in bf.layout.faces:
        wf_piece.seed(bf, faces=(face,))
    assert np.array_equal(wf_piece.values, wf_whole.values)


def test_seed_unknown_face_raises():
    sn = _cart2d_mesh()
    bf = BoundaryFlux.zeros_on(sn)
    wf = WavefrontFlux.zeros_on(sn)
    with pytest.raises(ValueError, match="bogus"):
        wf.seed(bf, faces=("bogus",))


def test_absorb_unknown_face_raises():
    sn = _cart2d_mesh()
    bf = BoundaryFlux.zeros_on(sn)
    wf = WavefrontFlux.zeros_on(sn)
    with pytest.raises(ValueError, match="bogus"):
        wf.absorb(bf, faces=("bogus",))


# ───────────────────────────────────────────────────────────────────────
# Axis-parametric API (§3a′ 3D-readiness): one path for 1-D / 2-D / 3-D
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "dims",
    [(10,), (6, 5), (4, 5, 3)],
    ids=["1d", "2d", "3d"],
)
def test_interior_layout_axis_parametric(dims):
    r"""The interior FaceLayout builds for 1-D, 2-D AND 3-D through ONE path
    (proving 3-D is a parameter, not a fork). No mesh / solver required."""
    N, ng = 8, 2
    layout = InteriorFaceSpace.interior_layout(N, ng, dims)
    # one face-normal field per axis
    assert len(layout.faces) == len(dims)
    # per-axis shape: dims[b]+1 along the face's own axis, dims[b] elsewhere
    names = ("x", "y", "z")
    for a in range(len(dims)):
        expected = (N, ng) + tuple(
            dims[b] + 1 if b == a else dims[b] for b in range(len(dims))
        )
        assert layout.faces[names[a]].shape == expected
    # flat total = Σ per-axis sizes; no gaps/overlaps (FaceLayout validates)
    assert layout.total_size == sum(s.flat_size for s in layout.faces.values())


def test_wavefront_flux_runs_1d_and_2d_one_path():
    r"""WavefrontFlux.zeros_on works for BOTH a 1-D and a 2-D mesh through the
    same constructor (axes=(0,) and (0,1))."""
    for sn, n_axes in [(_slab_mesh(), 1), (_cart2d_mesh(), 2)]:
        wf = WavefrontFlux.zeros_on(sn)
        assert len(wf.axes) == n_axes
        assert wf.values.shape == (wf.space.layout.total_size,)


# ───────────────────────────────────────────────────────────────────────
# Frozen + mesh-binding guards
# ───────────────────────────────────────────────────────────────────────


def test_frozen_attribute_rebind_raises():
    sn = _cart2d_mesh()
    wf = WavefrontFlux.zeros_on(sn)
    with pytest.raises(FrozenInstanceError):
        wf.values = np.zeros(wf.values.shape)  # type: ignore[misc]


def test_cross_mesh_arithmetic_raises():
    a = WavefrontFlux.zeros_on(_cart2d_mesh())
    b = WavefrontFlux.zeros_on(_cart2d_mesh())   # distinct SNMesh instance
    with pytest.raises(ValueError, match="distinct SNMesh"):
        _ = a + b


# ───────────────────────────────────────────────────────────────────────
# InteriorFaceSpace identity
# ───────────────────────────────────────────────────────────────────────


def test_interior_face_space_identity():
    sn = _cart2d_mesh()
    s1 = InteriorFaceSpace.from_mesh(sn)
    s2 = InteriorFaceSpace.from_mesh(sn)
    assert s1 == s2                       # (name, shape) identity
    assert hash(s1) == hash(s2)
    assert s1.name == "sn_wavefront"
    assert s1.face_names == ("x", "y")
