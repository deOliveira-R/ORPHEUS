"""Unit tests for Monte Carlo solver properties.

Tests structural properties of the MC geometry and solver:
- ConcentricPinCell material_id_at correctness
- SlabPinCell material_id_at correctness (backward compatibility)
- 1G homogeneous gives deterministic result (σ=0)
"""

import numpy as np
import pytest

from monte_carlo import (
    ConcentricPinCell, SlabPinCell, MCParams, solve_monte_carlo,
)
from derivations import get


# ── Geometry unit tests ──────────────────────────────────────────────

def test_concentric_pin_cell_materials():
    """ConcentricPinCell must return correct material for known positions."""
    geom = ConcentricPinCell(
        radii=[0.5, 0.8, 1.2],
        mat_ids=[2, 1, 0],
        pitch=3.0,
    )
    center = 1.5

    # At center → material 2 (innermost)
    assert geom.material_id_at(center, center) == 2

    # Just inside first boundary
    assert geom.material_id_at(center + 0.4, center) == 2

    # Between first and second boundary
    assert geom.material_id_at(center + 0.6, center) == 1

    # Outside all boundaries
    assert geom.material_id_at(center + 1.0, center) == 0


def test_slab_pin_cell_backward_compat():
    """SlabPinCell.default_pwr must match the original hardcoded geometry."""
    geom = SlabPinCell.default_pwr(pitch=3.6)

    # Original MATLAB regions:
    # fuel: 0.9 < x < 2.7, clad: 0.7-0.9 or 2.7-2.9, cool: rest
    assert geom.material_id_at(1.8, 1.0) == 2   # fuel
    assert geom.material_id_at(0.8, 1.0) == 1   # clad
    assert geom.material_id_at(2.8, 1.0) == 1   # clad
    assert geom.material_id_at(0.3, 1.0) == 0   # cool
    assert geom.material_id_at(3.2, 1.0) == 0   # cool


def test_concentric_4_regions():
    """ConcentricPinCell with 4 regions returns all 4 materials."""
    geom = ConcentricPinCell(
        radii=[0.3, 0.5, 0.7, 1.5],
        mat_ids=[3, 2, 1, 0],
        pitch=3.0,
    )
    center = 1.5

    assert geom.material_id_at(center, center) == 3
    assert geom.material_id_at(center + 0.4, center) == 2
    assert geom.material_id_at(center + 0.6, center) == 1
    assert geom.material_id_at(center + 1.0, center) == 0


# ── MC solver properties ─────────────────────────────────────────────

def test_1g_homogeneous_deterministic():
    """1G homogeneous MC gives σ=0 (every neutron sees the same XS)."""
    case = get("mc_cyl1D_1eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}

    # Use slab geometry with single material everywhere
    geom = SlabPinCell(boundaries=[], mat_ids=[0], pitch=3.6)
    params = MCParams(
        n_neutrons=100, n_inactive=20, n_active=100,
        seed=42, geometry=geom,
    )
    result = solve_monte_carlo(materials, params)

    assert result.sigma < 1e-10, (
        f"1G homogeneous MC should have σ≈0, got σ={result.sigma:.6e}"
    )
