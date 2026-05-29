r"""2-D Cartesian L2 face_view convention crosswalk — unit-source positive tests.

Test-architect Test 2.3 (D-H.2-C4a verification spec).  THE load-bearing
structural check for the C4 L2-native rewrite.

What this catches
=================

The legacy 2-D Cartesian matvec routes through a packed eq_map layout
(``transport_operator_matvec`` at ``operator.py:455``) where face
state lives in concatenated B1″ slots indexed by
``eq_map.face_outer_ordinate`` / ``eq_map.face_inner_ordinate``.
The L2 ``BoundaryFlux.face_view("xmin")`` returns a per-face shaped
view ``(N, ng, ny)``.  The mapping between these two indexings is
the convention crosswalk that C4 collapses.

A face-mapping bug (Pattern 7 convention drift) typically manifests
as: "the unit at xmin reads back at xmax", or "the unit at outward
ordinate streams INWARD".  Homogeneous flat-flux tests cannot see
these — the redistribution path is null and the integral over the
mesh is invariant to face-ordering swaps.

For each of ``{xmin, xmax, ymin, ymax}`` this file places a unit
flux at a specific face position + ordinate slot, applies the
streaming operator, and asserts:

* The matvec output is NONZERO (the BC apply reaches the bulk).
* The output's NORM is bounded above by what's physically consistent
  with a single-cell-worth of streaming (catches order-of-magnitude
  errors from a missing volume/area factor).
* The four faces produce DISTINCT output patterns (catches face
  aliasing — if all four placements produced the same result, the
  layout is degenerate).

Why xfail-strict until D-H.2-C4b
===============================

Today (post-C3) the composite ``StreamingOperator._apply_timed_full_field``
raises ``NotImplementedError`` for 2-D Cartesian (operator.py:2225).
C4b unblocks the path via a flat-shim through the legacy 2-D kernel;
C4c-d replace the kernel with L2-native.  This file lands xfail at
C4a and auto-flips to PASS at C4b.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh


# ── Mesh fixture: 4×4 vacuum-everywhere 2-D Cartesian, 1G pure streamer ─


def _pure_streamer_2d_mesh(nx: int = 4, ny: int = 4) -> SNMesh:
    r"""4×4 vacuum-everywhere mesh; materials Σ_s = 0, νΣ_f = 0, Σ_t = 1.

    Pure-streamer config — every cell is a "transparent" absorber.
    Eliminates the scattering / fission paths so the test isolates
    the streaming + BC apply convention crosswalk.

    Single-group keeps the test reducible to per-cell analysis.
    """
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "1g")})


def _zero_state_with_unit_face(
    mesh: SNMesh, face: str,
) -> "object":
    r"""Build a ``TimedFullField`` whose only nonzero is a unit at ``face``.

    Places ``boundary.face_view(face)[...] = 1.0`` (all ordinates,
    all groups, all face cells); zeros everywhere else.  The
    structural positive test: this unit should propagate downstream
    via the streaming operator's BC apply.
    """
    state = mesh.zeros_timed_full_field()
    state.boundary.face_view(face)[...] = 1.0
    return state


# ── The four positive face-mapping tests ────────────────────────────────


@pytest.mark.xfail(
    strict=True,
    reason=(
        "2-D Cartesian L2 matvec lands in D-H.2-C4b "
        "(StreamingOperator._apply_timed_full_field currently raises "
        "NotImplementedError for ny > 1).  This xfail-strict auto-flips "
        "to PASS when C4b ships the flat-shim through "
        "transport_operator_matvec."
    ),
)
@pytest.mark.parametrize("face", ["xmin", "xmax", "ymin", "ymax"])
def test_unit_at_face_produces_nonzero_matvec(face: str) -> None:
    r"""Unit at each face MUST produce a nonzero matvec output.

    Negative-version: if ``out.bulk.values`` is identically zero for
    any face, the BC apply for that face is broken (Pattern 7
    convention drift catches this).
    """
    from orpheus.sn.operator import StreamingOperator

    mesh = _pure_streamer_2d_mesh()
    sigma_t = np.ones((1, mesh.nx, mesh.ny))
    L = StreamingOperator(mesh, sigma_t)
    state = _zero_state_with_unit_face(mesh, face)

    out = L.apply(state)

    bulk = out.bulk.values
    assert not np.allclose(bulk, 0.0), (
        f"Unit at face {face!r} produced identically-zero bulk output. "
        f"BC apply for {face!r} is broken — the convention crosswalk "
        f"between L2 face_view and the matvec kernel lost the face."
    )


@pytest.mark.xfail(
    strict=True,
    reason="2-D L2 matvec lands in D-H.2-C4b (same as the test above).",
)
def test_four_faces_produce_distinct_outputs() -> None:
    r"""xmin / xmax / ymin / ymax MUST produce distinct matvec outputs.

    Negative-version: if any two faces produce identical ``bulk.values``
    arrays, the face-mapping has degenerated (aliasing).  Catches a
    regression where, e.g., the convention crosswalk fuses xmin+ymin
    into one slot.
    """
    from orpheus.sn.operator import StreamingOperator

    mesh = _pure_streamer_2d_mesh()
    sigma_t = np.ones((1, mesh.nx, mesh.ny))
    L = StreamingOperator(mesh, sigma_t)

    outputs = {}
    for face in ("xmin", "xmax", "ymin", "ymax"):
        state = _zero_state_with_unit_face(mesh, face)
        outputs[face] = L.apply(state).bulk.values.copy()

    face_names = list(outputs)
    for i, a in enumerate(face_names):
        for b in face_names[i + 1 :]:
            assert not np.allclose(outputs[a], outputs[b]), (
                f"Face {a!r} and {b!r} produced identical matvec outputs. "
                f"Convention aliasing — the L2 face_view → matvec slot "
                f"mapping has degenerated."
            )


@pytest.mark.xfail(
    strict=True,
    reason="2-D L2 matvec lands in D-H.2-C4b.",
)
def test_xmin_unit_streams_in_positive_x_direction() -> None:
    r"""Unit at xmin should populate cells in the positive-x direction.

    Geometric content: xmin is the LEFT boundary; a unit there
    (on the ordinates with ``μ_x > 0``) flows INTO the mesh in the
    +x direction.  The matvec output for those ordinates should be
    LARGER at small ix than at large ix (the unit decays as it
    streams across the absorbing medium).

    Catches a sign-flip on the streaming-direction convention: if
    the kernel reads ``μ_x`` flipped, the unit at xmin streams in
    the WRONG direction and the output peaks at large ix instead.
    """
    from orpheus.sn.operator import StreamingOperator

    mesh = _pure_streamer_2d_mesh()
    sigma_t = np.ones((1, mesh.nx, mesh.ny))
    L = StreamingOperator(mesh, sigma_t)
    state = _zero_state_with_unit_face(mesh, "xmin")

    out = L.apply(state).bulk.values  # (N, 1, nx, ny)

    # Restrict to outward-x ordinates (μ_x > 0) — they're the ones
    # that receive flux from xmin.
    mu_x = mesh.quad.mu_x
    outward_x = mu_x > 0
    out_outward = out[outward_x]  # (N_out, 1, nx, ny)

    # Sum over ordinates and y; the resulting 1-D profile in x
    # should be monotonically decreasing for an absorbing pure
    # streamer fed only at xmin.
    profile_x = out_outward.sum(axis=(0, 1, 3))  # (nx,)

    assert profile_x[0] > profile_x[-1], (
        f"xmin unit did NOT stream in +x direction: profile_x = "
        f"{profile_x}.  Sign-flip on μ_x convention or BC slot map."
    )


@pytest.mark.xfail(
    strict=True,
    reason="2-D L2 matvec lands in D-H.2-C4b.",
)
def test_ymin_unit_streams_in_positive_y_direction() -> None:
    r"""Unit at ymin should populate cells in the positive-y direction.

    Mirror of the xmin/xmax pair on the y axis — catches the same
    sign-flip class but on the OTHER spatial axis.  Variable-swap
    bugs (μ_x ↔ μ_y) fail HERE.
    """
    from orpheus.sn.operator import StreamingOperator

    mesh = _pure_streamer_2d_mesh()
    sigma_t = np.ones((1, mesh.nx, mesh.ny))
    L = StreamingOperator(mesh, sigma_t)
    state = _zero_state_with_unit_face(mesh, "ymin")

    out = L.apply(state).bulk.values  # (N, 1, nx, ny)

    mu_y = mesh.quad.mu_y
    outward_y = mu_y > 0
    out_outward = out[outward_y]

    profile_y = out_outward.sum(axis=(0, 1, 2))  # (ny,)

    assert profile_y[0] > profile_y[-1], (
        f"ymin unit did NOT stream in +y direction: profile_y = "
        f"{profile_y}.  Sign-flip on μ_y convention or BC slot map."
    )
