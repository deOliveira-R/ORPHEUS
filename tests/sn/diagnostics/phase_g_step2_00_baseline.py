"""Diagnostic: Step 1 didn't perturb the Phase F empirical baseline.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

Hypothesis to verify: tip dda6f28 (after Step 1) still reproduces the
Phase F closeout numbers on sphere_2g_3reg n=40:
  k_eff(SI) = 1.38069560
  k_eff(Kr) = 1.38464040
  k_diff   ≈ 0.286 %
  sf[0]/sf[1] (SI) ≈ 0.778
  sf[0]/sf[1] (Kr) ≈ 1.029

If any of these moves materially, Step 1 broke something and Step 2
shouldn't proceed.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn


def _build_sphere_2g_3reg(n_cells: int = 40):
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective,),
    )
    n_per_region = (n_cells // 4, n_cells // 2, n_cells // 4)
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
    )
    return dict(
        materials={0: fuel, 1: mod},
        mesh=mesh,
        quadrature=GaussLegendre1D.create(n_ordinates=8),
        scattering_order=0,
    )


@pytest.mark.slow
def test_phase_f_baseline_si_vs_krylov_n40():
    kwargs = _build_sphere_2g_3reg(40)

    result_si = solve_sn(inner_solver="source_iteration", **kwargs)
    result_kr = solve_sn(inner_solver="krylov", **kwargs)

    keff_si = float(result_si.keff)
    keff_kr = float(result_kr.keff)

    sf_si = result_si.scalar_flux[:, 0, 0]  # (nx, ny=1, ng) → (nx,) at g=0
    sf_kr = result_kr.scalar_flux[:, 0, 0]
    r01_si = float(sf_si[0] / sf_si[1])
    r01_kr = float(sf_kr[0] / sf_kr[1])
    k_diff_pct = 100.0 * (keff_kr - keff_si) / keff_si

    print(f"\nSI:  k_eff={keff_si:.8f}  sf[0]/sf[1]={r01_si:.4f}")
    print(f"Kr:  k_eff={keff_kr:.8f}  sf[0]/sf[1]={r01_kr:.4f}")
    print(f"k_diff = {k_diff_pct:.3f} %")

    # Phase F closeout numbers
    assert abs(keff_si - 1.38069560) < 1e-5, (
        f"k_eff(SI) drifted from Phase F baseline 1.38069560 → {keff_si}"
    )
    assert abs(keff_kr - 1.38464040) < 1e-5, (
        f"k_eff(Kr) drifted from Phase F baseline 1.38464040 → {keff_kr}"
    )
    assert abs(r01_si - 0.7776) < 1e-3, (
        f"r01(SI) drifted from Phase F baseline 0.7776 → {r01_si}"
    )
    assert abs(r01_kr - 1.0288) < 5e-3, (
        f"r01(Kr) drifted from Phase F baseline 1.0288 → {r01_kr}"
    )


if __name__ == "__main__":
    test_phase_f_baseline_si_vs_krylov_n40()
