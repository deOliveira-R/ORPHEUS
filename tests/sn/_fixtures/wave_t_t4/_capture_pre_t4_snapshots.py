"""Wave T step T.4a — pre-T.4 snapshot capture script.

Captures the numerical output of ``StreamingOperator.apply``,
``(L+C).apply``, and ``(L+C).solve`` on deterministic fixtures BEFORE
the T.4 lift rewires ``StreamingOperator.apply`` into the
``M_spatial + M_angular_redist`` decomposition.  The captured arrays
are loaded back by the L4-1..L4-7 regression tests in
``tests/sn/test_streaming_operator.py::TestPreT4RegressionSnapshot``
and ``tests/sn/test_invertible_operator.py`` (extended).

Per the T.4 verification spec
(``.claude/agent-memory/test-architect/wave_t_t4_streaming_verification_spec.md``)
§3 substep T.4a: pre-T.4 numerical regression baseline.  Without these
snapshots the L4-1..L4-7 tests have no pre-T.4 reference to compare
against; the principled-equivalence three-criteria gate requires
structurally-independent backing (k_∞ closed-form, L1 MMS) AT VALUE
LEVEL — these snapshots provide bit-identity per dispatch arm so any
reduction-order drift surfaces as a controlled nulp delta rather than
silent contamination.

Usage
-----

    .venv/bin/python tests/sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py

Writes ``pre_t4_snapshots.npz`` AND ``pre_t4_walltime.json`` to the
same directory.  Re-run only when the upstream test fixtures (XS
libraries, mesh dimensions, seed) change AND the change is
intentional — otherwise the snapshot is stale.

Determinism
-----------

* All input ψ use ``np.random.default_rng(<fixed_seed>)``.
* All meshes use fixed ``nx`` / ``ny`` / edges.
* All materials use programmatically-built ``make_mixture`` with
  hardcoded P0+P1 arrays (no library lookup that could drift).
* All quadratures use fixed order; Gauss-Legendre N=8 for 1-D,
  level-symmetric S6 for 2-D.

Output schema (see Wave-T T.4 verification spec §3 for the full
table; this script lands a focused subset covering the load-bearing
arms).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

# Ensure the orpheus package on the path when running standalone.
_REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(_REPO_ROOT))

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import Mesh1D, Mesh2D
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.transport.fields.boundary_flux import BoundaryFlux

OUT_NPZ = Path(__file__).parent / "pre_t4_snapshots.npz"
OUT_JSON = Path(__file__).parent / "pre_t4_walltime.json"

# Seed scheme — one base seed plus deterministic offsets per arm.
SEED_BASE = 20260531  # Wave T T.4a capture date


# ─────────────────────────────────────────────────────────────────────
# Material builders — hardcoded so the snapshots are independent of
# the xs_library state.  Mixture A (2G) parameters from xs_library
# region A 2G entry (asymmetric SigS; non-zero P1).
# ─────────────────────────────────────────────────────────────────────


def _mix_1g() -> "make_mixture":
    """Single-group P0 mixture — `numerical-bug-signatures` Mode 0 fixture."""
    return make_mixture(
        sig_t=np.array([1.0]),
        sig_c=np.array([0.2]),
        sig_f=np.array([0.1]),
        nu=np.array([2.5]),
        chi=np.array([1.0]),
        sig_s=np.array([[0.7]]),
    )


def _mix_2g_p1_asymmetric() -> "make_mixture":
    """2G P1 asymmetric SigS — load-bearing for ERR-002 family detection
    (see ``numerical-bug-signatures`` Signature 3 — transpose drift).

    Mirrors the T.3 fixture's asymmetric P0 + non-zero P1 + (n,2n)
    pattern, which exercises the M-M angular redistribution's group-
    axis coupling via ``total_xs`` and the per-group cell balance.
    """
    p0 = np.array([[0.38, 0.10], [0.05, 0.90]])  # asymmetric ⇒ ERR-002 detector
    p1 = np.array([[0.02, 0.01], [0.00, 0.04]])  # non-zero P1 anisotropy
    mix = make_mixture(
        sig_t=np.array([0.5, 1.0]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.5, 2.5]),
        chi=np.array([1.0, 0.0]),
        sig_s=p0,
    )
    mix.SigS = [csr_matrix(p0), csr_matrix(p1)]
    return mix


# ─────────────────────────────────────────────────────────────────────
# Mesh builders — slab / sphere / cylinder / 2-D Cartesian, with the
# BC variants the spec mandates.
# ─────────────────────────────────────────────────────────────────────


def _slab_mesh(
    *, ng: int, bc_left: BC, bc_right: BC, nx: int = 20, N: int = 8,
) -> SNMesh:
    mix = _mix_1g() if ng == 1 else _mix_2g_p1_asymmetric()
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        bc_left=bc_left,
        bc_right=bc_right,
    )
    quad = Quadrature.gauss_legendre(N)
    return SNMesh(mesh, quad, {0: mix})


def _sphere_mesh(*, ng: int, nx: int = 20, N: int = 8) -> SNMesh:
    mix = _mix_1g() if ng == 1 else _mix_2g_p1_asymmetric()
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),  # pole (r=0); structural, NOT physical
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(N)
    return SNMesh(mesh, quad, {0: mix})


def _cylinder_mesh(*, ng: int, nx: int = 20, sn_order: int = 4) -> SNMesh:
    """Cylinder requires a level-structured quadrature (μ-levels)."""
    mix = _mix_1g() if ng == 1 else _mix_2g_p1_asymmetric()
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),  # pole
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=sn_order)
    return SNMesh(mesh, quad, {0: mix})


def _cart2d_mesh(
    *, ng: int, bc_kind: str, nx: int = 6, ny: int = 6, sn_order: int = 4,
) -> SNMesh:
    """2-D Cartesian — exercises the ``MovingFrontierWindow.loss_action`` path (the matvec walk that since S6.3 lives on the loss representation, off the operator)."""
    mix = _mix_1g() if ng == 1 else _mix_2g_p1_asymmetric()
    bc = BC("reflective") if bc_kind == "specular" else BC("vacuum")
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=bc, bc_xmax=bc, bc_ymin=bc, bc_ymax=bc,
    )
    quad = Quadrature.level_symmetric(sn_order=sn_order)
    return SNMesh(mesh, quad, {0: mix})


# ─────────────────────────────────────────────────────────────────────
# ψ + q fixture constructors (deterministic).
# ─────────────────────────────────────────────────────────────────────


def _make_state(sn_mesh: SNMesh, *, seed: int) -> TimedFullField:
    """Build a TimedFullField with a deterministic random bulk ψ."""
    rng = np.random.default_rng(seed)
    N = sn_mesh.quad.N
    ng = sn_mesh.ng
    nx, ny = sn_mesh.nx, sn_mesh.ny
    psi_values = rng.uniform(0.05, 1.0, size=(N, ng, nx, ny))
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    from dataclasses import replace
    return replace(
        state,
        bulk=replace(state.bulk, values=psi_values),
    )


def _make_sigma_t(sn_mesh: SNMesh) -> np.ndarray:
    """Per-group per-cell σ_t (the StreamingOperator constructor input).

    Layout ``(ng, nx, ny)`` under Issue #196 PR-INDEX-3.  Reads from
    each cell's material σ_t via ``sn_mesh.mat_map`` (the (nx, ny)
    material-id array reshaped from Mesh1D.mat_ids / Mesh2D.mat_map).
    """
    ng = sn_mesh.ng
    nx, ny = sn_mesh.nx, sn_mesh.ny
    sig_t = np.empty((ng, nx, ny), dtype=float)
    mat_map = sn_mesh.mat_map  # (nx, ny)
    for ix in range(nx):
        for iy in range(ny):
            mid = int(mat_map[ix, iy])
            mat = sn_mesh.materials[mid]
            for g in range(ng):
                sig_t[g, ix, iy] = float(mat.SigT[g])
    return sig_t


def _build_L_C(sn_mesh: SNMesh) -> tuple[StreamingOperator, CollisionOperator]:
    """Build the leaf L (StreamingOperator) and C (CollisionOperator)."""
    sigma_t = _make_sigma_t(sn_mesh)
    L = StreamingOperator(sn_mesh, sigma_t)
    C = CollisionOperator(sn_mesh, sigma_t)
    return L, C


# ─────────────────────────────────────────────────────────────────────
# Capture
# ─────────────────────────────────────────────────────────────────────


def _capture_apply(name: str, sn_mesh: SNMesh, *, seed: int,
                   snapshots: dict[str, np.ndarray]) -> None:
    """Capture L.apply(ψ) bulk + boundary for one geometry arm."""
    L, _ = _build_L_C(sn_mesh)
    state = _make_state(sn_mesh, seed=seed)
    out = L.apply(state)
    snapshots[f"{name}_apply_bulk"] = out.bulk.values.copy()
    snapshots[f"{name}_apply_boundary"] = out.boundary.values.copy()
    snapshots[f"seed_psi_{name}"] = state.bulk.values.copy()


def _capture_LpC_apply(name: str, sn_mesh: SNMesh, *, seed: int,
                       snapshots: dict[str, np.ndarray]) -> None:
    """Capture (L+C).apply(ψ) — verifies algebra-decomposition invariant.

    After T.4c, ``M_spatial.apply + M_angular_redist.apply ==
    (L+C).apply`` value-equal to this snapshot (with σ_t·ψ subtraction
    unwinding at the apply boundary per Q3 decision γ).
    """
    L, C = _build_L_C(sn_mesh)
    LpC = L + C  # InvertibleOperator
    state = _make_state(sn_mesh, seed=seed)
    out = LpC.apply(state)
    snapshots[f"{name}_LpC_apply_bulk"] = out.bulk.values.copy()
    snapshots[f"{name}_LpC_apply_boundary"] = out.boundary.values.copy()


def _capture_LpC_solve(name: str, sn_mesh: SNMesh, *, seed: int,
                       snapshots: dict[str, np.ndarray]) -> None:
    """Capture (L+C).solve(q) — verifies .solve path is UNTOUCHED by T.4.

    Per Q5 the .solve path is OUT OF SCOPE for T.4.  This snapshot
    pins L4-6 to detect accidental sweep-path perturbation.
    """
    L, C = _build_L_C(sn_mesh)
    LpC = L + C  # InvertibleOperator
    q_state = _make_state(sn_mesh, seed=seed)
    out = LpC.solve(q_state)
    snapshots[f"{name}_LpC_solve_bulk"] = out.bulk.values.copy()
    snapshots[f"{name}_LpC_solve_boundary"] = out.boundary.values.copy()
    snapshots[f"seed_q_{name}"] = q_state.bulk.values.copy()


def _capture_perf_baseline(snapshots: dict[str, np.ndarray]) -> dict:
    """Capture pre-T.4 walltime on slab 2G P1 vacuum, 200 (L+C).apply iterations.

    Per spec §6 R8 (perf regression gate).  Baseline matched by L5-1.
    """
    import platform

    sn_mesh = _slab_mesh(ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
                         nx=40, N=8)
    L, C = _build_L_C(sn_mesh)
    LpC = L + C
    state = _make_state(sn_mesh, seed=SEED_BASE + 100)

    # Warmup
    for _ in range(5):
        LpC.apply(state)

    # Measure — 200 iterations gives mean + std.  Per spec ≥ 1000 is
    # nominal; 200 keeps the snapshot-script wallclock < 30 s.
    n_iters = 200
    samples = np.empty(n_iters, dtype=float)
    for i in range(n_iters):
        t0 = time.perf_counter()
        LpC.apply(state)
        samples[i] = time.perf_counter() - t0

    return {
        "fixture": "slab_2g_p1_asymmetric_vacuum_nx40_S8GL_apply",
        "n_iterations": n_iters,
        "walltime_seconds_median": float(np.median(samples)),
        "walltime_seconds_p95": float(np.percentile(samples, 95)),
        "walltime_seconds_mean": float(np.mean(samples)),
        "walltime_seconds_stdev": float(np.std(samples)),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "rng_seed": SEED_BASE + 100,
        "capture_date": "2026-05-31",
    }


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────


def main() -> None:
    snapshots: dict[str, np.ndarray] = {}

    print("Capturing 1-D arms (slab, sphere, cylinder) ...")
    # ── Slab arms ─────────────────────────────────────────────────
    _capture_apply(
        "slab_1g_vacuum",
        _slab_mesh(ng=1, bc_left=BC("vacuum"), bc_right=BC("vacuum")),
        seed=SEED_BASE + 1, snapshots=snapshots,
    )
    _capture_apply(
        "slab_2g_vacuum",
        _slab_mesh(ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum")),
        seed=SEED_BASE + 2, snapshots=snapshots,
    )
    # Spec mandated "white BC" to catch BC convention drift; the SN
    # `BOUNDARY_OPERATOR_REGISTRY` (orpheus/sn/geometry.py:155) only
    # registers ``vacuum`` and ``reflective`` today (the other 5
    # realizer-handled kinds — white / periodic / albedo /
    # prescribed_inflow / mixed — await SN-sweep-side wiring per the
    # registry docstring).  Reflective is the closest BC variant we
    # can land here; it exercises the D-B+1 tensor-network specular
    # path (the only existing production TP instance) and is
    # structurally distinct from vacuum (permutation vs incoming-mask).
    _capture_apply(
        "slab_2g_reflective",
        _slab_mesh(ng=2, bc_left=BC("reflective"), bc_right=BC("vacuum")),
        seed=SEED_BASE + 3, snapshots=snapshots,
    )

    # ── Sphere arms ───────────────────────────────────────────────
    _capture_apply(
        "sphere_1g",
        _sphere_mesh(ng=1),
        seed=SEED_BASE + 11, snapshots=snapshots,
    )
    _capture_apply(
        "sphere_2g",
        _sphere_mesh(ng=2),
        seed=SEED_BASE + 12, snapshots=snapshots,
    )

    # ── Cylinder arms ─────────────────────────────────────────────
    _capture_apply(
        "cyl_1g",
        _cylinder_mesh(ng=1),
        seed=SEED_BASE + 21, snapshots=snapshots,
    )
    _capture_apply(
        "cyl_2g",
        _cylinder_mesh(ng=2),
        seed=SEED_BASE + 22, snapshots=snapshots,
    )

    print("Capturing 2-D Cartesian arms (untouched by T.4 per Q1) ...")
    # ── 2-D Cartesian arms (defensive pin via A2D-1; Q1 = hybrid) ─
    _capture_apply(
        "cart2d_1g_specular",
        _cart2d_mesh(ng=1, bc_kind="specular"),
        seed=SEED_BASE + 31, snapshots=snapshots,
    )
    _capture_apply(
        "cart2d_1g_vacuum",
        _cart2d_mesh(ng=1, bc_kind="vacuum"),
        seed=SEED_BASE + 32, snapshots=snapshots,
    )
    _capture_apply(
        "cart2d_2g_specular",
        _cart2d_mesh(ng=2, bc_kind="specular"),
        seed=SEED_BASE + 33, snapshots=snapshots,
    )

    print("Capturing (L+C).apply algebra invariants ...")
    # ── (L+C).apply for 4 geometries — verifies algebra-decomposition ─
    _capture_LpC_apply(
        "slab_2g_vacuum",
        _slab_mesh(ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum")),
        seed=SEED_BASE + 41, snapshots=snapshots,
    )
    _capture_LpC_apply(
        "sphere_2g",
        _sphere_mesh(ng=2),
        seed=SEED_BASE + 42, snapshots=snapshots,
    )
    _capture_LpC_apply(
        "cyl_2g",
        _cylinder_mesh(ng=2),
        seed=SEED_BASE + 43, snapshots=snapshots,
    )
    _capture_LpC_apply(
        "cart2d_2g_specular",
        _cart2d_mesh(ng=2, bc_kind="specular"),
        seed=SEED_BASE + 44, snapshots=snapshots,
    )

    print("Capturing (L+C).solve Q5-out-of-scope contract ...")
    # ── (L+C).solve for 3 geometries — verifies .solve UNTOUCHED ──
    _capture_LpC_solve(
        "slab_2g_vacuum",
        _slab_mesh(ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum")),
        seed=SEED_BASE + 51, snapshots=snapshots,
    )
    _capture_LpC_solve(
        "sphere_2g",
        _sphere_mesh(ng=2),
        seed=SEED_BASE + 52, snapshots=snapshots,
    )
    _capture_LpC_solve(
        "cyl_2g",
        _cylinder_mesh(ng=2),
        seed=SEED_BASE + 53, snapshots=snapshots,
    )

    # ── Persist ────────────────────────────────────────────────────────
    np.savez(OUT_NPZ, **snapshots)
    print(f"wrote {OUT_NPZ}")
    for k, v in snapshots.items():
        print(f"  {k:50s} shape={v.shape} dtype={v.dtype}")

    print("\nCapturing perf baseline (200 iterations) ...")
    perf = _capture_perf_baseline(snapshots)
    OUT_JSON.write_text(json.dumps(perf, indent=2))
    print(f"wrote {OUT_JSON}")
    print(f"  median walltime: {perf['walltime_seconds_median']*1000:.3f} ms")
    print(f"  p95    walltime: {perf['walltime_seconds_p95']*1000:.3f} ms")


if __name__ == "__main__":
    main()
