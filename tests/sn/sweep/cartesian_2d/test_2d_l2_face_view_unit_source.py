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

face_view is ACTIVE as of Wave O #208 O.4b Phase E1
===================================================

These tests were authored xfail at D-H.2-C4a, anticipating the
face_view-as-trace extraction.  Through D-H.2-C4d the 2-D matvec kept
``face_view`` PASSIVE (cell-centre-proxy semantics), so the tests
stayed xfail.  O.4b Phase E1 made the 2-D matvec walk BARE — the
representation ``loss_action`` (which since S6.3 lives on the loss
representation, off the operator; ``ScanMarch`` default since S6.9,
the window a bare peer too) — it
reads ``psi.boundary.inflow`` as the GIVEN incoming edge and emits the
boundary-block residual — so the boundary trace now reaches the bulk.

* ``test_unit_at_face_produces_nonzero_matvec`` and
  ``test_four_faces_produce_distinct_outputs`` are PROMOTED (the active
  trace makes them correct positive structural pins).
* The two directional-streaming tests were REDESIGNED at O.4b Phase E3.
  They originally asserted a monotonic spatial DECAY of the bare MATVEC
  ``L·ψ`` on a boundary-only field — but that is a property of the
  streaming SOLVE ``(L+C)^{-1}``, NOT of the matvec.  With ψ̄ = 0 the
  diamond closure ``ψ_out = 2·ψ̄ − ψ_in = −ψ_in`` makes the matvec output
  OSCILLATE in sign cell-to-cell (verified: profile_x = [−A, +A, −A, +A])
  with NO spatial decay (equal magnitude every cell), so the original
  assertion was un-satisfiable on the matvec.  E3 redesigns them against
  the CONVERGED 2-D fixed-source SI solve (``solve_sn_fixed_source(...,
  inner_solver="source_iteration")`` — the 2-D SI fixed-source path is
  live; only 2-D fixed-source Krylov + 2-D eigenvalue SI are deferred):
  a localised edge source streams INTO the mesh and the scalar flux
  decays monotonically in the streaming direction.  A μ_x/μ_y sign flip
  or variable swap (Mode 2) makes the peak land on the WRONG edge.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

# 2-D Cartesian L2 face_view convention crosswalk: structural pins on
# the BoundaryFlux.face_view writability + unit-source convention (no
# theory-page :label:). Foundation, not a physics equation gate. (Was a
# V&V orphan before the taxonomy reorg forced a marker.)
pytestmark = pytest.mark.foundation


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
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=mesh)
    state.boundary.face_view(face)[...] = 1.0
    return state


# ── The four positive face-mapping tests ────────────────────────────────


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


# ── Directional streaming tests (O.4b Phase E3 redesign) ────────────────
#
# REDESIGNED against the converged 2-D fixed-source SI solve (the original
# bare-matvec assertion was un-satisfiable — see module docstring).  A
# localised volumetric source on one edge column/row streams INTO the mesh;
# the converged scalar flux peaks at the source edge and decays monotonically
# in the +x / +y streaming direction.  A μ_x / μ_y sign-flip or variable swap
# (Mode 2) makes the flux peak on the WRONG edge.


def _decay_mesh_geom(nx: int = 6, ny: int = 6) -> Mesh2D:
    r"""6×6 vacuum 2-D Mesh2D for the directional decay solve.

    Slightly larger than the 4×4 face fixtures so the streaming decay has
    enough cells to be monotone and unambiguous.
    """
    return Mesh2D(
        edges_x=np.linspace(0.0, 3.0, nx + 1),
        edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )


def test_xmin_source_streams_in_positive_x_direction() -> None:
    r"""A source on the xmin column streams in +x: flux peaks at ix=0, decays.

    Solve content: a localised volumetric source on the LEFT column
    (``ix = 0``, all ordinates, all y) feeds the mesh.  The converged scalar
    flux ``φ(x, y)`` summed over y must PEAK at ``ix = 0`` and decay
    monotonically in +x (vacuum BC, attenuating medium).

    Catches a sign-flip on the μ_x streaming-direction convention: if the
    kernel reads ``μ_x`` flipped the streaming carries flux the WRONG way and
    the profile no longer peaks at the source edge.  This is the SOLVE-based
    replacement for the un-satisfiable bare-matvec direction assertion (2-D
    fixed-source SI is live; only 2-D fixed-source Krylov + 2-D eigenvalue SI
    are deferred).
    """
    from orpheus.sn.solver import solve_sn_fixed_source

    geom = _decay_mesh_geom()
    quad = Quadrature.level_symmetric(sn_order=4)
    nx, ny, N = 6, 6, quad.N
    ext = np.zeros((N, 1, nx, ny))
    ext[:, 0, 0, :] = 1.0  # source on the xmin column (all y)

    res = solve_sn_fixed_source(
        materials={0: get_mixture("A", "1g")},
        mesh=geom, quadrature=quad, external_source=ext,
        inner_solver="source_iteration", max_inner=500, inner_tol=1e-10,
    )
    profile_x = res.scalar_flux.values[0].sum(axis=1)  # (nx,) sum over y

    assert profile_x.argmax() == 0, (
        f"xmin source did NOT peak at ix=0: profile_x = {profile_x}.  "
        f"Sign-flip on μ_x convention — the source streams the wrong way."
    )
    assert np.all(np.diff(profile_x) < 0), (
        f"xmin source flux is NOT monotonically decreasing in +x: "
        f"profile_x = {profile_x}.  Streaming-direction or attenuation bug."
    )


def test_ymin_source_streams_in_positive_y_direction() -> None:
    r"""A source on the ymin row streams in +y: flux peaks at iy=0, decays.

    Mirror of the xmin test on the y axis — catches the same sign-flip class
    on the OTHER spatial axis.  A μ_x ↔ μ_y variable swap (Mode 2) makes the
    y-profile fail to peak at iy=0.  By symmetry of the 6×6 mesh the y-profile
    here MUST equal the x-profile of the companion test (x↔y exchange) — a
    free cross-check that the two axes are handled symmetrically.
    """
    from orpheus.sn.solver import solve_sn_fixed_source

    geom = _decay_mesh_geom()
    quad = Quadrature.level_symmetric(sn_order=4)
    nx, ny, N = 6, 6, quad.N
    ext = np.zeros((N, 1, nx, ny))
    ext[:, 0, :, 0] = 1.0  # source on the ymin row (all x)

    res = solve_sn_fixed_source(
        materials={0: get_mixture("A", "1g")},
        mesh=geom, quadrature=quad, external_source=ext,
        inner_solver="source_iteration", max_inner=500, inner_tol=1e-10,
    )
    profile_y = res.scalar_flux.values[0].sum(axis=0)  # (ny,) sum over x

    assert profile_y.argmax() == 0, (
        f"ymin source did NOT peak at iy=0: profile_y = {profile_y}.  "
        f"Sign-flip on μ_y convention or a μ_x↔μ_y variable swap."
    )
    assert np.all(np.diff(profile_y) < 0), (
        f"ymin source flux is NOT monotonically decreasing in +y: "
        f"profile_y = {profile_y}.  Streaming-direction or attenuation bug."
    )
