r"""Pin the joint-batch behavior of :func:`ordinate_scan` for SLAB.

Issue #196 PR-INDEX-1.  The slab 1-D sweep, after the principled
index migration (``(N, ng, nx, ny)`` layout), is expected to invoke
:func:`~orpheus.sn.spatial.scan.ordinate_scan` **exactly twice per
sweep** for SLAB: once per chain direction, joint-batched over all
ordinates in that direction.

This is the load-bearing performance win:

* PR-INDEX-1 *before*: one scan call per ordinate ⇒ ``N/2`` calls per
  direction, ``N`` calls per sweep (e.g. N=16 ⇒ 16 calls).
* PR-INDEX-1 *after* (slab): one scan call per chain direction ⇒ 2
  calls per sweep, regardless of ``N`` or ``ng``.

CURVILINEAR (sphere/cylinder) is NOT joint-batched — the
Morel-Montry angular thread couples ordinates sequentially within a
μ-level (Hébert §3.9.4 Eqs. 3.437/3.439).  That restructuring is
deferred research-level work (plan §7).  Curvilinear keeps its
per-ordinate scan call count; we pin that here too as a safety
sentinel.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn import sweep as sweep_module
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.sources import PerOrdinateSource
from orpheus.sn.sweep import transport_sweep
from tests.sn._test_helpers import placeholder_materials


def _slab_setup(N: int = 4, nx: int = 8, ng: int = 2):
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    sn_mesh = SNMesh(mesh, GaussLegendre1D.create(n_ordinates=N), placeholder_materials(ng=ng))
    # Issue #196 PR-INDEX-5: Q principled (ng, nx, ny).
    # Issue #197 PR-TYPED-4: strict typed source.
    Q = PerOrdinateSource.from_isotropic(np.full((ng, nx, 1), 1.0), sn_mesh)
    sig_t = np.full((ng, nx, 1), 0.5)  # (ng, nx, ny) — PR-INDEX-3
    return Q, sig_t, sn_mesh


def _sphere_setup(N: int = 4, nx: int = 8, ng: int = 2):
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    sn_mesh = SNMesh(mesh, GaussLegendre1D.create(n_ordinates=N), placeholder_materials(ng=ng))
    # Issue #196 PR-INDEX-5: Q principled (ng, nx, ny).
    # Issue #197 PR-TYPED-4: strict typed source.
    Q = PerOrdinateSource.from_isotropic(np.full((ng, nx, 1), 1.0), sn_mesh)
    sig_t = np.full((ng, nx, 1), 0.5)  # (ng, nx, ny) — PR-INDEX-3
    return Q, sig_t, sn_mesh


def _count_ordinate_scan_calls(Q, sig_t, sn_mesh, boundary_flux):
    """Run one sweep and return the number of ``ordinate_scan`` calls."""
    calls = []
    original = sweep_module.ordinate_scan

    def tracking_scan(a, b, psi_0):
        calls.append((a.shape, b.shape, np.asarray(psi_0).shape))
        return original(a, b, psi_0)

    sweep_module.ordinate_scan = tracking_scan
    try:
        transport_sweep(Q, sig_t, sn_mesh, boundary_flux)
    finally:
        sweep_module.ordinate_scan = original
    return calls


@pytest.mark.l0
def test_slab_joint_batch_one_scan_per_chain_direction() -> None:
    """SLAB invokes ordinate_scan exactly twice per sweep — one per chain."""
    Q, sig_t, sn_mesh = _slab_setup(N=4, nx=8, ng=2)
    calls = _count_ordinate_scan_calls(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())
    assert len(calls) == 2, (
        f"Slab sweep should invoke ordinate_scan exactly 2 times "
        f"(one per chain direction); got {len(calls)}: {calls}"
    )


@pytest.mark.l0
def test_slab_joint_batch_independent_of_N() -> None:
    """SLAB scan count is invariant in N (always 2)."""
    for N in (2, 4, 8, 16):
        Q, sig_t, sn_mesh = _slab_setup(N=N, nx=8, ng=2)
        calls = _count_ordinate_scan_calls(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())
        assert len(calls) == 2, (
            f"Slab N={N}: expected 2 scan calls, got {len(calls)}"
        )


@pytest.mark.l0
def test_slab_joint_batch_independent_of_ng() -> None:
    """SLAB scan count is invariant in ng (always 2)."""
    for ng in (1, 2, 4):
        Q, sig_t, sn_mesh = _slab_setup(N=8, nx=10, ng=ng)
        calls = _count_ordinate_scan_calls(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())
        assert len(calls) == 2, (
            f"Slab ng={ng}: expected 2 scan calls, got {len(calls)}"
        )


@pytest.mark.l0
def test_slab_joint_batch_call_shapes() -> None:
    """SLAB scan inputs carry the joint-batch shape ``(nx, K, ng)``.

    Pins the principled layout — leading axis is the cell/scan axis,
    middle axis is the ordinate batch (``K = N/2`` for symmetric GL),
    trailing axis is the group batch.
    """
    Q, sig_t, sn_mesh = _slab_setup(N=8, nx=10, ng=3)
    calls = _count_ordinate_scan_calls(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())
    assert len(calls) == 2
    nx, K, ng = 10, 4, 3  # GL-8 symmetric: 4 ordinates per direction
    for a_shape, b_shape, psi0_shape in calls:
        assert a_shape == (nx, K, ng), (
            f"a shape expected (nx, K, ng) = {(nx, K, ng)}; got {a_shape}"
        )
        assert b_shape == (nx, K, ng), (
            f"b shape expected (nx, K, ng) = {(nx, K, ng)}; got {b_shape}"
        )
        assert psi0_shape == (K, ng), (
            f"psi_0 shape expected (K, ng) = {(K, ng)}; got {psi0_shape}"
        )


@pytest.mark.l0
def test_sphere_per_ordinate_scan_safety_sentinel() -> None:
    """SPHERE invokes ordinate_scan once per non-degenerate ordinate.

    The M-M angular thread (Hébert §3.9.4) couples ordinates
    sequentially within a μ-level.  This sentinel pins the
    per-ordinate dispatch — a future joint-batch optimisation for
    curvilinear is deferred research-level work and would intentionally
    change this count.  Until then, drift in this count means an
    unintended structural change.
    """
    N = 4
    Q, sig_t, sn_mesh = _sphere_setup(N=N, nx=8, ng=2)
    calls = _count_ordinate_scan_calls(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())
    # Sphere has one level, all N ordinates are non-degenerate for
    # this test configuration.  Each ordinate gets exactly one scan
    # call.
    assert len(calls) == N, (
        f"Sphere N={N}: expected {N} scan calls (one per ordinate); "
        f"got {len(calls)}"
    )
