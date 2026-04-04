"""Unit tests for collision probability matrix properties.

Tests algebraic properties that the CP matrices must satisfy,
independent of any eigenvalue solver:
- Row sums = 1 (neutron conservation)
- Reciprocity: Σ_t V P_ij = Σ_t V P_ji
- Non-negativity: P ≥ 0
- 1-region limit reproduces homogeneous k
"""

import numpy as np
import pytest

from geometry import CoordSystem, Mesh1D, Zone, mesh1d_from_zones
from collision_probability import (
    CPParams,
    _compute_slab_cp_group, _compute_cp_group,
    _build_ki_tables, _chord_half_lengths,
)
from derivations._xs_library import get_mixture, get_xs


# ── Fixtures ──────────────────────────────────────────────────────────

def _slab_pinf_1g():
    """Build P_inf for a 1G 2-region slab."""
    xs_a = get_xs("A", "1g")
    xs_b = get_xs("B", "1g")
    mesh = mesh1d_from_zones([
        Zone(outer_edge=0.5, mat_id=0, n_cells=1),
        Zone(outer_edge=1.0, mat_id=1, n_cells=1),
    ], coord=CoordSystem.CARTESIAN)
    sig_t_g = np.array([xs_a["sig_t"][0], xs_b["sig_t"][0]])
    return _compute_slab_cp_group(sig_t_g, mesh), sig_t_g, mesh.widths


def _cyl_pinf_1g():
    """Build P_inf for a 1G 2-region cylinder."""
    r_fuel, r_cell = 0.5, 1.0
    mesh = mesh1d_from_zones([
        Zone(outer_edge=r_fuel, mat_id=0, n_cells=1),
        Zone(outer_edge=r_cell, mat_id=1, n_cells=1),
    ], coord=CoordSystem.CYLINDRICAL)
    xs_a = get_xs("A", "1g")
    xs_b = get_xs("B", "1g")
    sig_t_g = np.array([xs_a["sig_t"][0], xs_b["sig_t"][0]])
    ki_x, _, ki4_v = _build_ki_tables(20000, 50.0)

    radii = mesh.edges[1:]
    gl_pts, gl_wts = np.polynomial.legendre.leggauss(64)
    breakpoints = mesh.edges
    y_all, w_all = [], []
    for seg in range(len(breakpoints) - 1):
        a, b = breakpoints[seg], breakpoints[seg + 1]
        y_all.append(0.5 * (b - a) * gl_pts + 0.5 * (b + a))
        w_all.append(0.5 * (b - a) * gl_wts)
    y_pts = np.concatenate(y_all)
    y_wts = np.concatenate(w_all)
    chords = _chord_half_lengths(radii, y_pts)

    P_inf = _compute_cp_group(sig_t_g, mesh, chords, y_pts, y_wts, ki_x, ki4_v)
    return P_inf, sig_t_g, mesh.volumes


# ── Row sums (neutron conservation) ──────────────────────────────────

def test_slab_row_sums():
    """P_inf row sums must equal 1 (every neutron collides somewhere)."""
    P_inf, _, _ = _slab_pinf_1g()
    row_sums = P_inf.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)


def test_cylinder_row_sums():
    """P_inf row sums must equal 1 for cylindrical geometry."""
    P_inf, _, _ = _cyl_pinf_1g()
    row_sums = P_inf.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)


# ── Reciprocity ──────────────────────────────────────────────────────

def test_slab_reciprocity():
    """Σ_t[i] V[i] P[i,j] = Σ_t[j] V[j] P[j,i] for all i≠j."""
    P_inf, sig_t, V = _slab_pinf_1g()
    N = P_inf.shape[0]
    for i in range(N):
        for j in range(i + 1, N):
            lhs = sig_t[i] * V[i] * P_inf[i, j]
            rhs = sig_t[j] * V[j] * P_inf[j, i]
            assert abs(lhs - rhs) < 1e-10, (
                f"Reciprocity violated: Σ_t[{i}]V[{i}]P[{i},{j}]={lhs:.6e} "
                f"≠ Σ_t[{j}]V[{j}]P[{j},{i}]={rhs:.6e}"
            )


def test_cylinder_reciprocity():
    """Reciprocity for cylindrical geometry."""
    P_inf, sig_t, V = _cyl_pinf_1g()
    N = P_inf.shape[0]
    for i in range(N):
        for j in range(i + 1, N):
            lhs = sig_t[i] * V[i] * P_inf[i, j]
            rhs = sig_t[j] * V[j] * P_inf[j, i]
            assert abs(lhs - rhs) < 1e-10, (
                f"Reciprocity violated at ({i},{j}): {lhs:.6e} ≠ {rhs:.6e}"
            )


# ── Non-negativity ───────────────────────────────────────────────────

def test_slab_non_negativity():
    """All collision probabilities must be non-negative."""
    P_inf, _, _ = _slab_pinf_1g()
    assert np.all(P_inf >= 0), f"Negative P_inf entry: min={P_inf.min():.6e}"


def test_cylinder_non_negativity():
    """All collision probabilities must be non-negative."""
    P_inf, _, _ = _cyl_pinf_1g()
    assert np.all(P_inf >= 0), f"Negative P_inf entry: min={P_inf.min():.6e}"


# ── 1-region homogeneous limit ───────────────────────────────────────

def test_slab_1region_homogeneous_limit():
    """1-region CP slab must give P=1 (all neutrons collide in the only region)."""
    mesh = mesh1d_from_zones([
        Zone(outer_edge=0.5, mat_id=0, n_cells=1),
    ], coord=CoordSystem.CARTESIAN)
    sig_t_g = np.array([1.0])
    P_inf = _compute_slab_cp_group(sig_t_g, mesh)
    np.testing.assert_allclose(P_inf[0, 0], 1.0, atol=1e-10)


def test_cylinder_1region_homogeneous_limit():
    """1-region CP cylinder must give P=1."""
    r_cell = 1.0
    mesh = mesh1d_from_zones([
        Zone(outer_edge=r_cell, mat_id=0, n_cells=1),
    ], coord=CoordSystem.CYLINDRICAL)
    sig_t_g = np.array([1.0])
    ki_x, _, ki4_v = _build_ki_tables(20000, 50.0)

    radii = mesh.edges[1:]
    gl_pts, gl_wts = np.polynomial.legendre.leggauss(64)
    y_pts = 0.5 * r_cell * gl_pts + 0.5 * r_cell
    y_wts = 0.5 * r_cell * gl_wts
    chords = _chord_half_lengths(radii, y_pts)

    P_inf = _compute_cp_group(sig_t_g, mesh, chords, y_pts, y_wts, ki_x, ki4_v)
    np.testing.assert_allclose(P_inf[0, 0], 1.0, atol=1e-10)
