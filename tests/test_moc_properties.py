"""Unit tests for MOC solver properties.

Tests structural properties of the Method of Characteristics solution:
- Particle balance: production / absorption = keff
- Scalar flux positivity
- from_annular geometry factory produces correct material layout
"""

import numpy as np
import pytest

from derivations import get
from method_of_characteristics import MoCGeometry, solve_moc


def test_particle_balance():
    """For homogeneous fill, production / absorption = keff."""
    case = get("moc_cyl1D_1eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix, 1: mix, 2: mix}
    geom = MoCGeometry.default_pwr()
    result = solve_moc(materials, geom, max_outer=200)

    flux = result.scalar_flux  # (n, n, ng)
    vol = geom.volume  # (n, n)

    production = np.sum(flux[:, :, 0] * mix.SigP[0] * vol)
    absorption = np.sum(flux[:, :, 0] * (mix.SigC[0] + mix.SigF[0]) * vol)

    k_balance = production / absorption
    np.testing.assert_allclose(
        k_balance, result.keff, rtol=1e-4,
        err_msg=f"Particle balance: {k_balance:.6f} ≠ keff={result.keff:.6f}",
    )


def test_flux_positivity():
    """Scalar flux must be positive everywhere."""
    case = get("moc_cyl1D_1eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix, 1: mix, 2: mix}
    geom = MoCGeometry.default_pwr()
    result = solve_moc(materials, geom, max_outer=200)

    assert np.all(result.scalar_flux > 0), (
        f"Non-positive flux: min={result.scalar_flux.min():.6e}"
    )


def test_from_annular_material_layout():
    """from_annular must assign correct materials based on distance from center."""
    geom = MoCGeometry.from_annular(
        radii=[0.3, 0.5, 0.8],
        mat_ids=[2, 1, 0],
        pitch=2.0,
        n_cells=20,
    )
    center = 1.0  # pitch / 2
    delta = 0.1   # pitch / n_cells

    # Check a few cells
    for ix in range(20):
        for iy in range(20):
            cx = (ix + 0.5) * delta
            cy = (iy + 0.5) * delta
            r = np.sqrt((cx - center)**2 + (cy - center)**2)
            if r <= 0.3:
                expected = 2
            elif r <= 0.5:
                expected = 1
            else:
                expected = 0
            assert geom.mat_map[ix, iy] == expected, (
                f"Cell ({ix},{iy}): r={r:.3f}, expected mat={expected}, "
                f"got {geom.mat_map[ix, iy]}"
            )


def test_from_annular_homogeneous():
    """from_annular with single material should run and match analytical k."""
    case = get("moc_cyl1D_1eg_1rg")
    mix = next(iter(case.materials.values()))
    geom = MoCGeometry.from_annular(
        radii=[1.0],
        mat_ids=[0],
        pitch=2.0,
        n_cells=10,
    )
    result = solve_moc({0: mix}, geom, max_outer=200)
    assert abs(result.keff - case.k_inf) < 1e-3, (
        f"from_annular homogeneous: keff={result.keff:.6f} vs {case.k_inf:.6f}"
    )
