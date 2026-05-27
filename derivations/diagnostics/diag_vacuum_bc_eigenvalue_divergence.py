"""Diagnostic: SN vacuum-BC eigenvalue divergence — minimal reproducer.

Companion to:
    .claude/agent-memory/numerics-investigator/vacuum_bc_eigenvalue_divergence.md

Failing test:
    tests/sn/test_boundary_conditions.py::TestSNBCSweepBehavior::
    test_vacuum_keff_lower_than_reflective

Pre-existing on base refactor/sn-operator-algebra@62994ad.

The three probes below discriminate the surviving hypotheses:
  P1 — inner_solver="krylov": discriminates SI-specific bug from
        operator-construction bug
  P2 — fixed-source vacuum: isolates within-group SI from outer power
        iteration
  P3 — scipy.sparse.linalg.eigs spectrum: discriminates operator-side bug
        from power-iteration-interaction bug
"""
from __future__ import annotations
import sys
import pathlib
import warnings

# Force the script to load the COLOCATED orpheus, not whichever orpheus
# the venv's `pip install -e .` happened to bind to. (`pytest` uses
# <rootdir> auto-discovery via pyproject.toml so this isn't an issue
# for tests, but standalone diagnostics need the explicit prepend.)
_THIS = pathlib.Path(__file__).resolve()
_REPO = _THIS.parents[2]  # derivations/diagnostics/<this> -> repo
sys.path.insert(0, str(_REPO))

import numpy as np


def _build_2eg_slab(boundary, n_cells=20, length_cm=2.0):
    from orpheus.derivations.reference_values import get
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    geom_mesh = Mesh1D.from_geometry(
        StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=length_cm),),
            bcs=(boundary, boundary),
        ),
        region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.gauss_legendre(4)
    return materials, geom_mesh, quad


def probe_si_vacuum():
    """Reproduce the failure with SI inner solver."""
    from orpheus.geometry import BC
    from orpheus.sn.solver import solve_sn
    materials, mesh, quad = _build_2eg_slab(BC.vacuum)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = solve_sn(materials, mesh, quad)
    print(f"[P_SI] vacuum keff           = {result.keff:.10f}  (expected < 1.875)")
    print(f"[P_SI] terminal psi max      = "
          f"{float(np.max(np.abs(result.angular_flux.values))):.3e}")
    print(f"[P_SI] terminal phi max      = "
          f"{float(np.max(np.abs(result.scalar_flux.values))):.3e}")
    print(f"[P_SI] n_outer               = {result.history.n_outer}")
    print(f"[P_SI] keff_history last 5   = "
          f"{[f'{k:.6e}' for k in result.history.keff_history[-5:]]}")


def probe_krylov_vacuum():
    """Probe 1: SI-specific vs operator-construction bug."""
    from orpheus.geometry import BC
    from orpheus.sn.solver import solve_sn
    materials, mesh, quad = _build_2eg_slab(BC.vacuum)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = solve_sn(materials, mesh, quad, inner_solver="krylov")
    print(f"[P_Kr] krylov vacuum keff    = {result.keff:.10f}  (expected < 1.875)")
    print(f"[P_Kr] terminal psi max      = "
          f"{float(np.max(np.abs(result.angular_flux.values))):.3e}")
    print(f"[P_Kr] terminal phi max      = "
          f"{float(np.max(np.abs(result.scalar_flux.values))):.3e}")
    print(f"[P_Kr] n_outer               = {result.history.n_outer}")
    print(f"[P_Kr] keff_history last 5   = "
          f"{[f'{k:.6e}' for k in result.history.keff_history[-5:]]}")


def probe_reflective_ref():
    """Baseline: reflective should give k_inf = 1.875."""
    from orpheus.geometry import BC
    from orpheus.sn.solver import solve_sn
    materials, mesh, quad = _build_2eg_slab(BC.reflective)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = solve_sn(materials, mesh, quad)
    print(f"[P_Re] reflective keff       = {result.keff:.10f}  (expected = 1.875)")
    print(f"[P_Re] terminal psi max      = "
          f"{float(np.max(np.abs(result.angular_flux.values))):.3e}")
    print(f"[P_Re] n_outer               = {result.history.n_outer}")


if __name__ == "__main__":
    print("=" * 70)
    print("Probe baseline: reflective BC (should converge to k_inf=1.875)")
    print("=" * 70)
    probe_reflective_ref()
    print()
    print("=" * 70)
    print("Probe SI: vacuum BC + SI inner solver (REPRODUCES FAILURE)")
    print("=" * 70)
    probe_si_vacuum()
    print()
    print("=" * 70)
    print("Probe Krylov: vacuum BC + Krylov inner solver (discriminator P1)")
    print("=" * 70)
    probe_krylov_vacuum()
