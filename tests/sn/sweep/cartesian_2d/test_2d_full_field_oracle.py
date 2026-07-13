r"""End-to-end ``windowed ≡ full-field`` 2-D equivalence — the storage-B oracle.

Phase 5b replaced the full per-axis interior face cochain with a rolling
``_MovingFrontier`` window for BOTH the 2-D sweep (`_sweep_jacobi`) and
its matvec twin (`MovingFrontierWindow.loss_action`). The retained
full-field paths — `FullFieldWavefront.sweep` + `.loss_action` (the d-generic
spine; since S6.4(d) full-cochain interior KERNELS on the shared
``_OctantWalk`` frame) — are the VERIFICATION ORACLES: the fuller view (every
interior face retained during the walk) that the moving frontier is
cross-checked against END-TO-END.

This is STRONGER than `test_sweep_graph_window_equivalence.py` (which pins only
the per-octant WALK at the graph level): it pins the whole sweep
/ matvec ORCHESTRATOR — boundary plumbing included — bit-for-bit. Both share
the cell kernel
(`DiamondDifference.cell_kernel_batch`/`residual_kernel_batch`) AND (since
S6.4) the octant frame, so the MATH cannot drift; only the interior storage
policy differs, which is exactly what this pins.  (History: through S6.4(c)
the oracle frames carried the cochain as a typed ``WavefrontFlux`` with
whole-trace ι_*/ι* seed/absorb; the (d) fold moved the boundary algebra onto
the shared frame, and the type was RETIRED at (f) — the cochain concept lives
on as ``_MovingFrontier`` (the front) + ``_octant_face_cochain`` (the
history).)

``foundation`` — a software invariant (window ≡ full-field), no theory label.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.loss_representation import FullFieldWavefront, MovingFrontierWindow
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation

# (nx, ny, level, ng, bc)
CASES = [
    (8, 8, 4, 2, "reflective"),
    (8, 8, 4, 2, "vacuum"),
    (12, 7, 6, 2, "reflective"),
    (5, 9, 4, 4, "reflective"),
]


def _build_mesh(nx, ny, lvl, ng, bc):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC(bc), bc_xmax=BC(bc), bc_ymin=BC(bc), bc_ymax=BC(bc),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(lvl), {0: get_mixture("A", f"{ng}g")})


def _random_sig_t(rng, ng, nx, ny):
    """Heterogeneous, strictly positive Σ_t (Cardinal Rule ≥2G + het)."""
    return rng.uniform(0.3, 3.0, size=(ng, nx, ny))


def _seed_random_inflow(rng, boundary):
    """Fill every boundary face slab with a random inflow trace (so the
    orchestrator's boundary handling is genuinely exercised)."""
    for face in boundary.layout.faces:
        fv = boundary.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_sweep_window_equals_full_field_end_to_end(nx, ny, lvl, ng, bc):
    """MovingFrontierWindow.sweep (windowed) ≡ FullFieldWavefront.sweep
    (oracle), bit-for-bit: angular flux, scalar flux, AND the post-sweep
    boundary trace.

    S6.4(d): the oracle leg drives the representation's ``sweep`` (the former
    free function ``_sweep_full_field`` dissolved into the shared
    ``_OctantWalk`` frame × the full-cochain kernel — retirement = test
    migration)."""
    rng = np.random.default_rng(abs(hash((nx, ny, lvl, ng, bc))) % (2**32))
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    N = sn_mesh.quad.N
    sig_t = _random_sig_t(rng, ng, nx, ny)
    Q = rng.uniform(0.0, 2.0, size=(N, ng, nx, ny))

    bf_win = AngularBoundaryFlux.zeros_on(sn_mesh)
    _seed_random_inflow(rng, bf_win)
    bf_full = AngularBoundaryFlux.from_mesh(bf_win.values.copy(), sn_mesh)

    ang_w, scal_w = MovingFrontierWindow(sn_mesh).sweep(Q, sig_t, bf_win)
    ang_f, scal_f = FullFieldWavefront(sn_mesh).sweep(Q, sig_t, bf_full)

    np.testing.assert_array_equal(ang_w, ang_f, err_msg="angular flux")
    np.testing.assert_array_equal(scal_w, scal_f, err_msg="scalar flux")
    np.testing.assert_array_equal(
        bf_win.values, bf_full.values, err_msg="post-sweep boundary trace",
    )


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_matvec_window_equals_full_field_end_to_end(nx, ny, lvl, ng, bc):
    """MovingFrontierWindow.loss_action (windowed) ≡ FullFieldWavefront.loss_action
    (oracle), bit-for-bit: bulk (L+C)ψ AND boundary-block residual.

    S6.3: both representations' ``loss_action`` return ``(L+C)ψ`` (the operator
    does the ``−C`` glue in :meth:`apply`); the window≡full comparison is
    unchanged — two storage policies over ONE walk.

    S6.4(a): SUPERSEDES the retired A2D-1 source-hash pin
    (``TestT4dApply2DCartesianSourceHashPin``) — this ``assert_array_equal``
    output oracle is the tripwire on the 2-D matvec body, now the shared
    ``_OctantWalk`` apply frame + the window's interior kernel.  A source-hash
    on a SHARED walk trips on every legitimate refactor with no behavior
    signal; any actual drift in the relocated body fails THIS test against the
    structurally-distinct full-field oracle (the NON-SQUARE cases catch an
    x↔y swap).  Doubles as the gate-memo §5 relocation-identity pin for every
    S6.4 sub-step that touches the octant frame."""
    rng = np.random.default_rng((abs(hash((nx, ny, lvl, ng, bc))) % (2**32)) ^ 9)
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    N = sn_mesh.quad.N
    sig_t = _random_sig_t(rng, ng, nx, ny)

    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
    state.interior.values[...] = rng.uniform(-1.0, 1.0, size=state.interior.values.shape)
    _seed_random_inflow(rng, state.boundary)

    out_win = MovingFrontierWindow(sn_mesh).loss_action(sig_t, state)
    out_full = FullFieldWavefront(sn_mesh).loss_action(sig_t, state)

    np.testing.assert_array_equal(
        out_win.interior.values, out_full.interior.values, err_msg="bulk residual",
    )
    np.testing.assert_array_equal(
        out_win.boundary.values, out_full.boundary.values,
        err_msg="boundary-block residual",
    )
