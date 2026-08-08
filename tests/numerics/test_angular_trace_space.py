r"""Tests for :mod:`orpheus.numerics.spaces.angular_trace_space`.

The unified :class:`AngularTraceSpace` (#205 / #201) — ONE whole-boundary
trace function space carrying the signed :math:`\Omega\cdot\hat n` per
face. Inflow / outflow are *selectors* over it, not separate types.

L0 tests cover structure / identity / dtype / shape / selectors /
error raising. L1 tests cover the directional-predicate correctness
:math:`\mathrm{sign}(\Omega\cdot\hat n_f)` against hand-computed
expectations. The eps-gap guard pins the principled tangential
tolerance against every shipped quadrature.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace, TANGENTIAL_EPS
from orpheus.numerics.quadrature import Quadrature


# ─────────────────────────────────────────────────────────────────────
# Helpers — build the (mesh, quadrature, layout) triple the unified
# AngularTraceSpace consumes. Layouts mirror SNMesh.boundary_face_layout:
# slab xmin/xmax, curvilinear xmax-only, 2-D xmin/xmax/ymin/ymax.
# ─────────────────────────────────────────────────────────────────────


def _mesh1d(coord: CoordSystem, n: int = 4) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, n + 1),
        mat_ids=np.zeros(n, dtype=int),
        coord=coord,
    )


def _mesh2d(coord: CoordSystem, nx: int = 3, ny: int = 3) -> Mesh2D:
    return Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=coord,
    )


def _slab_layout(N: int, ng: int = 1) -> FaceLayout:
    return FaceLayout.from_named_shapes([("xmin", (N, ng)), ("xmax", (N, ng))])


def _curvilinear_layout(N: int, ng: int = 1) -> FaceLayout:
    return FaceLayout.from_named_shapes([("xmax", (N, ng))])


def _cartesian2d_layout(N: int, ng: int, nx: int, ny: int) -> FaceLayout:
    return FaceLayout.from_named_shapes([
        ("xmin", (N, ng, ny)), ("xmax", (N, ng, ny)),
        ("ymin", (N, ng, nx)), ("ymax", (N, ng, nx)),
    ])


def _trace_2d(quad, nx: int = 3, ny: int = 3, ng: int = 1) -> AngularTraceSpace:
    return AngularTraceSpace.from_quadrature_and_layout(
        quad, _cartesian2d_layout(quad.N, ng, nx, ny),
    )


def _trace_1d(coord: CoordSystem, quad, ng: int = 1) -> AngularTraceSpace:
    if coord is CoordSystem.CARTESIAN:
        layout = _slab_layout(quad.N, ng)
    else:
        layout = _curvilinear_layout(quad.N, ng)
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


# ─────────────────────────────────────────────────────────────────────
# L0 — structure, identity, frozen
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
def test_trace_space_constructible_and_frozen():
    """L0: AngularTraceSpace builds via factory; mutation raises; is a FunctionSpace."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad)
    assert isinstance(space, AngularTraceSpace)
    assert isinstance(space, FunctionSpace)
    with pytest.raises(dataclasses.FrozenInstanceError):
        space.name = "evil"  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        space.shape = (1,)  # type: ignore[misc]


@pytest.mark.l0
def test_shape_is_whole_boundary_flat_total_size():
    """L0: shape == (layout.total_size,) — the whole boundary, flat."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad, nx=3, ny=3, ng=2)
    # 4 faces: 2 x-faces (N*ng*ny) + 2 y-faces (N*ng*nx).
    expected = 2 * quad.N * 2 * 3 + 2 * quad.N * 2 * 3
    assert space.shape == (expected,)
    assert space.layout.total_size == expected


@pytest.mark.l0
def test_face_names_track_layout_order():
    """L0: face_names is the layout's ordered faces (omega_dot_n row order)."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad)
    assert space.face_names == ("xmin", "xmax", "ymin", "ymax")
    assert space.omega_dot_n.shape == (4, quad.N)


@pytest.mark.l0
def test_identity_independent_of_leaf_data():
    """L0: two trace spaces of the same (name, shape) compare equal —
    layout / omega_dot_n are ``compare=False`` leaf-data (FunctionSpace
    identity convention)."""
    quad = Quadrature.lebedev(11)
    a = _trace_2d(quad)
    b = _trace_2d(quad)
    assert a == b
    assert hash(a) == hash(b)
    np.testing.assert_array_equal(a.omega_dot_n, b.omega_dot_n)


# ─────────────────────────────────────────────────────────────────────
# L1 — directional-predicate correctness
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("inflow-mask-discrete")
def test_2d_cartesian_lebedev_inflow_selectors_per_face():
    """L1: per-face inflow indices match the hand-computed mu signs."""
    quad = Quadrature.lebedev(11)
    eps = TANGENTIAL_EPS
    space = _trace_2d(quad)
    # xmin: outward -x → inflow iff Ω·n = -mu_x < 0 iff mu_x > eps.
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("xmin"), np.flatnonzero(quad.mu_x > eps),
    )
    # xmax: outward +x → inflow iff mu_x < -eps.
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("xmax"), np.flatnonzero(quad.mu_x < -eps),
    )
    # ymin / ymax along the y-axis.
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("ymin"), np.flatnonzero(quad.mu_y > eps),
    )
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("ymax"), np.flatnonzero(quad.mu_y < -eps),
    )


@pytest.mark.l1
def test_outflow_selectors_are_sign_flipped_inflow():
    """L1: outflow at a face is the +sign half (complement of inflow on
    a non-tangential quadrature)."""
    quad = Quadrature.lebedev(11)
    eps = TANGENTIAL_EPS
    space = _trace_2d(quad)
    np.testing.assert_array_equal(
        space.outflow_indices_for_face("xmin"), np.flatnonzero(quad.mu_x < -eps),
    )
    np.testing.assert_array_equal(
        space.outflow_indices_for_face("xmax"), np.flatnonzero(quad.mu_x > eps),
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-040")
@pytest.mark.verifies("ordinate-partition-inflow-outflow")
def test_axis_aligned_ordinates_excluded_from_both_selectors():
    """L1: pure-axis ordinates (e.g. (0,0,±1)) are tangential to a
    perpendicular face — in NEITHER selector."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad)
    z_axis = np.where(
        (np.abs(quad.mu_x) < 1e-12) & (np.abs(quad.mu_y) < 1e-12)
        & (np.abs(quad.mu_z) > 0.5)
    )[0]
    assert z_axis.size >= 2, "expected (0,0,±1) ordinates in Lebedev 11"
    inflow = space.inflow_indices_for_face("xmin")
    outflow = space.outflow_indices_for_face("xmin")
    assert not np.intersect1d(inflow, z_axis).size
    assert not np.intersect1d(outflow, z_axis).size


@pytest.mark.l0
@pytest.mark.verifies("inflow-mask-discrete")
def test_1d_slab_inflow_bit_identical_to_mu_signs():
    """L0: 1-D slab xmin / xmax inflow matches the legacy left / right
    masks bit-for-bit (xmin↔old 'left', xmax↔old 'right')."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(CoordSystem.CARTESIAN, quad)
    assert space.face_names == ("xmin", "xmax")
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("xmin"),
        np.flatnonzero(quad.mu_x > TANGENTIAL_EPS),
    )
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("xmax"),
        np.flatnonzero(quad.mu_x < -TANGENTIAL_EPS),
    )


@pytest.mark.l1
@pytest.mark.verifies("ordinate-partition-inflow-outflow")
def test_inflow_xor_outflow_complementary_for_gl_1d():
    """L0: GL ordinates are strictly in (-1, 1) → no tangentials; inflow
    and outflow partition every ordinate per face."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(CoordSystem.CARTESIAN, quad)
    for face in space.face_names:
        inflow = set(space.inflow_indices_for_face(face).tolist())
        outflow = set(space.outflow_indices_for_face(face).tolist())
        assert inflow.isdisjoint(outflow)
        assert inflow | outflow == set(range(quad.N))


# ─────────────────────────────────────────────────────────────────────
# L0 — curvilinear: ONE boundary (the face-naming reconciliation)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL])
def test_curvilinear_has_only_outer_xmax_face(coord):
    """L0: a solid sphere / cylinder has exactly ONE boundary face
    (``xmax`` = outer radius). The pole r=0 is the angular closure's
    regularity condition, NOT a BC face — so there is no inner face
    (this is the latent-bug fix: the legacy trace fabricated a phantom
    ``left`` mask for r=0)."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(coord, quad)
    assert space.face_names == ("xmax",)
    # xmax inflow on curvilinear == the same predicate as slab's xmax.
    np.testing.assert_array_equal(
        space.inflow_indices_for_face("xmax"),
        np.flatnonzero(quad.mu_x < -TANGENTIAL_EPS),
    )
    with pytest.raises(ValueError, match="Unknown face"):
        space.inflow_indices_for_face("xmin")


@pytest.mark.l1
def test_curvilinear_xmax_matches_cartesian_xmax():
    """L1: the outer-radius predicate is identical to slab xmax — the
    per-face inflow predicate is a property of mu_x and the outward
    normal, both shared across 1-D coord systems."""
    quad = Quadrature.gauss_legendre(8)
    cart = _trace_1d(CoordSystem.CARTESIAN, quad)
    sph = _trace_1d(CoordSystem.SPHERICAL, quad)
    cyl = _trace_1d(CoordSystem.CYLINDRICAL, quad)
    np.testing.assert_array_equal(
        cart.inflow_indices_for_face("xmax"), sph.inflow_indices_for_face("xmax"),
    )
    np.testing.assert_array_equal(
        sph.inflow_indices_for_face("xmax"), cyl.inflow_indices_for_face("xmax"),
    )


# ─────────────────────────────────────────────────────────────────────
# L0 — error contracts
# ─────────────────────────────────────────────────────────────────────


# C5.3 (#225): the former ``test_2d_cylindrical_raises`` retired WITH the
# gate it pinned — AngularTraceSpace is geometry-blind (it never sees a mesh), so
# the 2-D-cylindrical refusal lives where the geometry enters: SNMesh
# construction (``axes_from_legacy_mesh`` raises NotImplementedError for
# 2-D non-Cartesian; pinned in
# tests/sn/primitives/test_axis_native_construction.py).


@pytest.mark.l0
def test_unknown_face_in_layout_raises():
    """L0: a layout naming a face absent from the normal table raises."""
    quad = Quadrature.gauss_legendre(8)
    # C5.3 (#225): "zmin" is a KNOWN face (the normal table derives
    # from AXIS_NAMES and carries all six axis-aligned faces) — the
    # unknown-face probe needs a name outside the axis world.
    bad = FaceLayout.from_named_shapes([("wmin", (quad.N, 1))])
    with pytest.raises(ValueError, match="Unknown face name"):
        AngularTraceSpace.from_quadrature_and_layout(quad, bad)


@pytest.mark.l0
def test_selector_on_unknown_face_raises():
    """L0: selecting a face not in the layout raises ValueError."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad)
    with pytest.raises(ValueError, match="Unknown face"):
        space.inflow_indices_for_face("bogus")


@pytest.mark.l0
def test_omega_dot_n_dtype_is_float():
    """L0: the signed projection table is float (not int / bool)."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad)
    assert space.omega_dot_n.dtype == np.float64


# ─────────────────────────────────────────────────────────────────────
# Foundation — principled tangential tolerance (the eps-gap guard)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_eps_is_four_machine_eps():
    """The tangential tolerance is the principled 4·machine_eps, NOT a
    hand-tuned magic number."""
    assert TANGENTIAL_EPS == 4.0 * np.finfo(np.float64).eps


@pytest.mark.foundation
def test_eps_sits_in_the_round_off_to_genuine_gap():
    r"""Across every shipped quadrature, the tangential eps must sit
    strictly between the round-off cluster (nominally-tangential cosines,
    which are EXACTLY 0.0) and the smallest GENUINE direction cosine.
    This is the durable form of ``eps_probe.py``: a future quadrature
    cannot silently violate the gap that makes the inflow/outflow masks
    bit-identical to the legacy ``1e-15`` / ``1e-12`` tolerances.
    """
    quads = [Quadrature.gauss_legendre(n) for n in (2, 3, 4, 5, 7, 8, 16, 32, 64)]
    for order in (3, 5, 7, 11, 17, 29, 53):
        try:
            quads.append(Quadrature.lebedev(order))
        except Exception:  # noqa: BLE001 — order may be unavailable
            pass
    # The family that MOTIVATED #325 — product rules carried the
    # trig-round-off tangential cosines this eps existed to absorb.
    # Post-E3 (roots-of-unity nodes) their tangential set is EXACTLY
    # 0.0, so they must sit in the same empty-band regime as GL and
    # Lebedev; enumerating them here is what makes this gate a claim
    # about EVERY shipped rule rather than the two families that were
    # never at risk. Odd n_phi exercises the quarter-point zeros.
    quads += [
        Quadrature.product(4, 8),
        Quadrature.product(6, 12),
        Quadrature.product(4, 5),
        Quadrature.folded_product(4, 8),
        Quadrature.folded_product(4, 16),
    ]
    quads += [Quadrature.level_symmetric(n) for n in (2, 4, 8, 12)]

    min_genuine = np.inf
    for q in quads:
        for axis in ("mu_x", "mu_y", "mu_z"):
            mu = getattr(q, axis, None)
            if mu is None:
                continue
            amu = np.abs(np.asarray(mu, dtype=float))
            # Nominally-tangential cosines are exactly 0.0 (quadrature
            # symmetry) — never in the (0, eps] round-off band.
            assert not np.any((amu > 0.0) & (amu <= TANGENTIAL_EPS)), (
                f"{type(q).__name__} {axis} has a cosine in the round-off "
                f"band (0, eps] — the exactly-zero assumption is violated"
            )
            genuine = amu[amu > TANGENTIAL_EPS]
            if genuine.size:
                min_genuine = min(min_genuine, float(genuine.min()))

    # eps must be far below the smallest genuine cosine (empirically
    # 2.44e-2, ~13 orders above eps).
    assert TANGENTIAL_EPS < min_genuine
    assert min_genuine > 1e-3, (
        f"smallest genuine |cosine| = {min_genuine:.3e} is suspiciously "
        f"small — investigate before trusting the gap"
    )


# ─────────────────────────────────────────────────────────────────────
# Foundation — the partial-current boundary metric G_s = |Ω·n_f| ⊙ w_n
# (Wave O / O.2b, #208).  Phase 4 populates inner_product_weights here;
# it is the metric under which the BoundaryOperator Hilbert adjoints are
# correct.  These pin the data shape/values; the discriminating proof
# that it makes the adjoint correct is the §2 reciprocity gate (4.3).
# ─────────────────────────────────────────────────────────────────────


def _reference_trace_metric(space: AngularTraceSpace, quad) -> np.ndarray:
    """Independent flat ``|Ω·n_f|·w_n`` build (per-face, broadcast over all
    trailing axes) — the cross-check ground for the production metric."""
    layout = space.layout
    ref = np.zeros((layout.total_size,), dtype=float)
    for f_idx, face in enumerate(layout.faces):
        slot = layout.faces[face]
        face_w = np.abs(space.omega_dot_n[f_idx]) * np.asarray(quad.weights)
        broadcast = np.broadcast_to(
            face_w.reshape((face_w.shape[0],) + (1,) * (len(slot.shape) - 1)),
            slot.shape,
        )
        ref[slot.offset : slot.offset + slot.flat_size] = broadcast.reshape(-1)
    return ref


@pytest.mark.foundation
@pytest.mark.parametrize(
    "coord",
    [CoordSystem.CARTESIAN, CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL],
)
def test_trace_metric_is_cosine_weighted_1d(coord):
    """O.2b: the 1-D trace metric is the partial-current weight
    ``|Ω·n_f|·w_n`` (NOT Euclidean ``None``) — bit-exact against an
    independent per-face build, on slab (2 faces) + curvilinear (1 face)."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(coord, quad, ng=2)
    assert space.inner_product_weights is not None
    np.testing.assert_array_equal(
        space.inner_product_weights, _reference_trace_metric(space, quad),
    )


@pytest.mark.foundation
def test_trace_metric_is_cosine_weighted_2d():
    """O.2b: the angular weight broadcasts across the group AND the
    along-edge spatial axis of a 2-D face slot ``(N, ng, n_face_cells)``."""
    quad = Quadrature.lebedev(11)
    space = _trace_2d(quad, nx=3, ny=3, ng=2)
    assert space.inner_product_weights is not None
    np.testing.assert_array_equal(
        space.inner_product_weights, _reference_trace_metric(space, quad),
    )


@pytest.mark.foundation
def test_trace_inner_product_is_weighted_not_euclidean():
    """O.2b: ``inner_product`` computes ``Σ |Ω·n|·w_n · a·b`` — equals the
    explicit weighted sum AND differs from the bare Euclidean ``Σ a·b``
    (the metric is load-bearing, not a no-op)."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(CoordSystem.CARTESIAN, quad, ng=2)
    rng = np.random.default_rng(2026)
    a = rng.standard_normal(space.shape)
    b = rng.standard_normal(space.shape)
    w = space.inner_product_weights
    assert space.inner_product(a, b) == pytest.approx(float(np.sum(w * a * b)))
    assert space.inner_product(a, b) != pytest.approx(float(np.sum(a * b)))


@pytest.mark.foundation
def test_trace_metric_group_independent():
    """O.2b: the metric is purely angular — every group column of a face
    slot equals ``|Ω·n_f|·w_n`` (no energy dependence)."""
    quad = Quadrature.gauss_legendre(8)
    space = _trace_1d(CoordSystem.CARTESIAN, quad, ng=3)
    slot = space.layout.faces["xmax"]
    w_face = slot.slice_view(space.inner_product_weights)  # (N, ng)
    f_idx = space.face_names.index("xmax")
    expected_col = np.abs(space.omega_dot_n[f_idx]) * np.asarray(quad.weights)
    for g in range(w_face.shape[1]):
        np.testing.assert_array_equal(w_face[:, g], expected_col)


# ─────────────────────────────────────────────────────────────────────
# C5.3 (#225) — z faces in the AXIS_NAMES-derived normal table
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_face_normals_z_derived_from_axis_names():
    """The normal parse carries all six axis-aligned faces (axis, sign).

    C5-G9: derived from AXIS_NAMES, not hand-listed — a table that
    silently lacked zmin/zmax was the pre-C5.3 d=3 blocker.

    Migrated at **B3.4c**: the module-private ``_FACE_NORMALS`` dict this
    asserted against was retired when the ``min``/``max`` ↔ sign convention
    collapsed onto the one bijection in
    :mod:`orpheus.numerics.face_layout`. The CLAIM is unchanged and is
    asserted here against the primitive this module now reads (the bijection's
    own laws — round-trip, involution — live in ``tests/numerics/
    test_face_layout.py``); what moved is where the answer comes from.
    """
    from orpheus.numerics.face_layout import face_normal

    expected = {
        "xmin": (0, -1), "xmax": (0, +1),
        "ymin": (1, -1), "ymax": (1, +1),
        "zmin": (2, -1), "zmax": (2, +1),
    }
    np.testing.assert_equal(
        {face: face_normal(face) for face in expected}, expected,
    )


@pytest.mark.foundation
def test_omega_dot_n_zmax_is_plus_mu_z():
    """C5-G10: the zmax row is exactly +mu_z, zmin exactly −mu_z."""
    quad = Quadrature.level_symmetric(sn_order=4)
    layout = FaceLayout.from_named_shapes([
        ("zmin", (quad.N, 1)), ("zmax", (quad.N, 1)),
    ])
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    np.testing.assert_array_equal(trace.omega_dot_n[0], -np.asarray(quad.mu_z))
    np.testing.assert_array_equal(trace.omega_dot_n[1], +np.asarray(quad.mu_z))


@pytest.mark.foundation
def test_omega_dot_n_all_six_faces_distinct_axes():
    """C5-G11: each face row maps to ITS axis cosine (Mode-2 swap detector).

    A diagonal cubature (level-symmetric) has genuinely distinct
    mu_x/mu_y/mu_z arrays, so a zmin→mu_x (the pre-C4 reflection bug
    class) or mu_y↔mu_z swap is an array mismatch, not a coincidence.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    faces = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
    layout = FaceLayout.from_named_shapes(
        [(f, (quad.N, 1)) for f in faces],
    )
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    mu = [np.asarray(quad.mu_x), np.asarray(quad.mu_y), np.asarray(quad.mu_z)]
    # Non-vacuity: the three cosine arrays must genuinely differ.
    if np.array_equal(mu[0], mu[1]) or np.array_equal(mu[1], mu[2]):
        pytest.fail("cubature does not distinguish the axis cosines")
    for row, face in enumerate(faces):
        axis, sign = {"x": 0, "y": 1, "z": 2}[face[0]], (-1 if face.endswith("min") else +1)
        np.testing.assert_array_equal(trace.omega_dot_n[row], sign * mu[axis])


@pytest.mark.foundation
def test_z_face_on_two_cosine_quadrature_fails_loud():
    """A z-face against a quadrature with NO genuine third cosine RAISES.

    C5.5 fail-loud guard (the C5.3 review carry, re-fixed after the
    first guard proved dead code): the per-axis cosines are PROPERTIES
    that zero-pad past the cubature's intrinsic dimensionality (1-D
    Gauss-Legendre carries ``mu_z == zeros(N)``, never an absent
    attribute), so the guard discriminates on VALUE — all-zero cosines
    on a layout-named normal axis are a rank-mismatch, not a legitimate
    all-tangential face. Without the raise, every ordinate at that face
    silently classifies as neither inflow nor outflow.
    """
    quad = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(np.asarray(quad.mu_z), 0.0)  # non-vacuity
    layout = FaceLayout.from_named_shapes([("zmin", (quad.N, 1))])
    with pytest.raises(ValueError, match="rank-mismatch"):
        AngularTraceSpace.from_quadrature_and_layout(quad, layout)
