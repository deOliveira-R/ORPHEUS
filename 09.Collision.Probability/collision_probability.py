"""Collision probability (CP) method for neutron transport.

Solves the multi-group neutron transport equation using the collision
probability method with white boundary condition for an infinite lattice.

Two geometries are supported:

* **Slab** — 1D half-cell with E₃ exponential-integral kernel.
  Requires a :class:`~geometry.mesh.Mesh1D` with
  ``coord = CoordSystem.CARTESIAN``.
* **Concentric cylinders** — Wigner-Seitz cell with Ki₃/Ki₄
  Bickley-Naylor kernel and numerical y-quadrature.
  Requires a :class:`~geometry.mesh.Mesh1D` with
  ``coord = CoordSystem.CYLINDRICAL``.

Both share the same power iteration: once the P_inf matrices are built,
the eigenvalue solve is geometry-independent.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.special import expn

from data.macro_xs.mixture import Mixture
from data.macro_xs.cell_xs import CellXS, assemble_cell_xs
from geometry import CoordSystem, Mesh1D
from numerics.eigenvalue import power_iteration


# ═══════════════════════════════════════════════════════════════════════
# Parameters and result container
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class CPParams:
    """Solver parameters for the collision probability method."""

    max_outer: int = 500
    keff_tol: float = 1e-6
    flux_tol: float = 1e-5
    n_ki_table: int = 20000
    ki_max: float = 50.0
    n_quad_y: int = 64


@dataclass
class CPResult:
    """Results of a collision probability calculation."""

    keff: float
    keff_history: list[float]
    flux: np.ndarray
    flux_fuel: np.ndarray
    flux_clad: np.ndarray
    flux_cool: np.ndarray
    geometry: Mesh1D
    eg: np.ndarray
    elapsed_seconds: float


# ═══════════════════════════════════════════════════════════════════════
# CP solver class (satisfies EigenvalueSolver protocol)
# ═══════════════════════════════════════════════════════════════════════

class CPSolver:
    """Geometry-independent CP eigenvalue solver.

    Once the infinite-lattice CP matrices P_inf are built (by the
    geometry-specific kernel), the eigenvalue iteration is identical
    for slab and cylindrical geometries:

        φ_g = P_inf_g^T · (V · Q_g) / (Σ_t · V)

    where V is the cell volume array (thicknesses for slab, areas for
    cylinder) and Q is the total source (fission + scattering + n,2n).
    """

    def __init__(
        self,
        P_inf: np.ndarray,
        xs: CellXS,
        volumes: np.ndarray,
        mat_ids: np.ndarray,
        materials: dict[int, Mixture],
        keff_tol: float = 1e-6,
        flux_tol: float = 1e-5,
    ) -> None:
        self.P_inf = P_inf        # (N, N, ng)
        self.xs = xs
        self.volumes = volumes    # (N,)
        self.mat_ids = mat_ids
        self.ng = xs.sig_t.shape[1]
        self.N = xs.sig_t.shape[0]
        self.keff_tol = keff_tol
        self.flux_tol = flux_tol

        # Cache scattering and (n,2n) matrices per material
        self._scat_mats = {mid: materials[mid].SigS[0] for mid in materials}
        self._n2n_mats = {mid: materials[mid].Sig2 for mid in materials}

    def initial_flux_distribution(self) -> np.ndarray:
        return np.ones((self.N, self.ng))

    def compute_fission_source(
        self, flux_distribution: np.ndarray, keff: float,
    ) -> np.ndarray:
        fission_rate = np.sum(self.xs.sig_p * flux_distribution, axis=1)
        return self.xs.chi * fission_rate[:, np.newaxis] / keff

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        N, ng = self.N, self.ng

        # Total source = fission + scattering + (n,2n)
        Q = fission_source.copy()
        for k in range(N):
            mid = self.mat_ids[k]
            Q[k, :] += self._scat_mats[mid].T @ flux_distribution[k, :]
            Q[k, :] += 2.0 * (self._n2n_mats[mid].T @ flux_distribution[k, :])

        # Apply CP matrices: φ_g = P_inf_g^T · (V · Q_g) / (Σ_t · V)
        phi = np.empty((N, ng))
        for g in range(ng):
            source = self.volumes * Q[:, g]
            phi[:, g] = self.P_inf[:, :, g].T @ source

        denom = self.xs.sig_t * self.volumes[:, np.newaxis]
        pos = denom > 0
        phi[pos] = phi[pos] / denom[pos]
        phi[~pos] = 0.0

        # Numerical conditioning (prevent overflow in subsequent iterations)
        phi *= 1.0 / np.max(phi)

        return phi

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        v = self.volumes[:, np.newaxis]
        production = np.sum(self.xs.sig_p * flux_distribution * v)
        absorption = np.sum(self.xs.sig_a * flux_distribution * v)
        return float(production / absorption)

    def converged(
        self, keff: float, keff_old: float,
        flux_distribution: np.ndarray, flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        if iteration <= 2:
            return False
        delta_k = abs(keff - keff_old)
        delta_phi = np.max(np.abs(flux_distribution - flux_old)) / \
            max(np.max(np.abs(flux_distribution)), 1e-30)

        if iteration <= 5 or iteration % 10 == 0:
            print(f"    iter {iteration:4d}  keff = {keff:.6f}  "
                  f"dk = {delta_k:.2e}  dphi = {delta_phi:.2e}")

        if delta_k < self.keff_tol and delta_phi < self.flux_tol:
            print(f"    iter {iteration:4d}  keff = {keff:.6f}  Converged.")
            return True
        return False


# ═══════════════════════════════════════════════════════════════════════
# Shared post-processing
# ═══════════════════════════════════════════════════════════════════════

def _volume_averaged_fluxes(
    phi: np.ndarray,
    volumes: np.ndarray,
    mat_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Volume-averaged flux per material (fuel=2, clad=1, cool=0)."""
    ng = phi.shape[1]
    vol_fuel = volumes[mat_ids == 2].sum()
    vol_clad = volumes[mat_ids == 1].sum()
    vol_cool = volumes[mat_ids == 0].sum()

    flux_fuel = np.zeros(ng)
    flux_clad = np.zeros(ng)
    flux_cool = np.zeros(ng)

    for k in range(len(mat_ids)):
        v = volumes[k]
        mid = mat_ids[k]
        if mid == 2:
            flux_fuel += phi[k, :] * v / vol_fuel
        elif mid == 1:
            flux_clad += phi[k, :] * v / vol_clad
        else:
            flux_cool += phi[k, :] * v / vol_cool

    return flux_fuel, flux_clad, flux_cool


# ═══════════════════════════════════════════════════════════════════════
# Slab E₃ kernel
# ═══════════════════════════════════════════════════════════════════════

def _e3(x):
    """Vectorised E_3(x) = integral_0^1 mu exp(-x/mu) dmu."""
    return expn(3, np.maximum(x, 0.0))


def _compute_slab_cp_group(
    sig_t_g: np.ndarray,
    mesh: Mesh1D,
) -> np.ndarray:
    """Within-cell CP matrix for one energy group in slab geometry.

    Uses the E₃ second-difference formula.  The half-cell has a reflective
    centre boundary and a white (isotropic re-entry) outer boundary.

    Returns the infinite-lattice CP matrix P_inf (N, N).
    """
    N = mesh.N
    t = mesh.widths
    tau = sig_t_g * t

    bnd_pos = np.zeros(N + 1)
    for k in range(N):
        bnd_pos[k + 1] = bnd_pos[k] + tau[k]

    rcp = np.zeros((N, N))

    for i in range(N):
        sti = sig_t_g[i]
        tau_i = tau[i]
        if sti <= 0:
            continue

        self_same = 0.5 * sti * (2 * t[i] - (2.0 / sti) * (0.5 - _e3(tau_i)))
        rcp[i, i] += self_same

        for j in range(N):
            tau_j = tau[j]

            if j > i:
                gap_d = bnd_pos[j] - bnd_pos[i + 1]
            elif j < i:
                gap_d = bnd_pos[i] - bnd_pos[j + 1]
            else:
                gap_d = None

            if gap_d is not None:
                gap_d = max(gap_d, 0.0)
                dd = (_e3(gap_d) - _e3(gap_d + tau_i)
                      - _e3(gap_d + tau_j) + _e3(gap_d + tau_i + tau_j))
            else:
                dd = 0.0

            gap_c = bnd_pos[i] + bnd_pos[j]
            dc = (_e3(gap_c) - _e3(gap_c + tau_i)
                  - _e3(gap_c + tau_j) + _e3(gap_c + tau_i + tau_j))

            rcp[i, j] += 0.5 * (dd + dc)

    P_cell = np.zeros((N, N))
    for i in range(N):
        if sig_t_g[i] * t[i] > 0:
            P_cell[i, :] = rcp[i, :] / (sig_t_g[i] * t[i])

    P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)
    P_in = sig_t_g * t * P_out
    P_inout = max(1.0 - P_in.sum(), 0.0)

    P_inf = P_cell.copy()
    if P_inout < 1.0:
        P_inf += np.outer(P_out, P_in) / (1.0 - P_inout)

    return P_inf


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical Ki₃/Ki₄ kernel
# ═══════════════════════════════════════════════════════════════════════

def _build_ki_tables(
    n_pts: int = 20000,
    x_max: float = 50.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tabulate Ki_3 and Ki_4 on a uniform grid.

    Ki_3(x) = int_0^{pi/2} exp(-x/sin t) sin t dt
    Ki_4(x) = int_x^inf Ki_3(t) dt

    Returns (x_grid, ki3_vals, ki4_vals).
    """
    x_grid = np.linspace(0, x_max, n_pts)
    ki3_vals = np.empty(n_pts)
    ki3_vals[0] = 1.0

    for i in range(1, n_pts):
        ki3_vals[i], _ = quad(
            lambda t, xx=x_grid[i]: np.exp(-xx / np.sin(t)) * np.sin(t),
            0, np.pi / 2,
        )

    dx = x_grid[1] - x_grid[0]
    ki4_vals = np.cumsum(ki3_vals[::-1])[::-1] * dx
    ki4_vals[-1] = 0.0

    return x_grid, ki3_vals, ki4_vals


def _ki4_lookup(x, x_grid, ki4_vals):
    """Vectorised Ki_4 lookup."""
    return np.interp(x, x_grid, ki4_vals, right=0.0)


def _chord_half_lengths(radii, y_pts):
    """Half-chord lengths l_k(y) for each annular region.  Shape (N, n_y)."""
    N = len(radii)
    chords = np.zeros((N, len(y_pts)))
    r_inner = np.zeros(N)
    r_inner[1:] = radii[:-1]
    y2 = y_pts**2

    for k in range(N):
        r_out, r_in = radii[k], r_inner[k]
        outer = np.sqrt(np.maximum(r_out**2 - y2, 0.0))
        if r_in > 0:
            inner = np.sqrt(np.maximum(r_in**2 - y2, 0.0))
            mask = y_pts < r_in
            chords[k, mask] = outer[mask] - inner[mask]
        mask_p = (y_pts >= r_in) & (y_pts < r_out)
        chords[k, mask_p] = outer[mask_p]

    return chords


def _compute_cp_group(
    sig_t_g: np.ndarray,
    mesh: Mesh1D,
    chords: np.ndarray,
    y_pts: np.ndarray,
    y_wts: np.ndarray,
    ki_x: np.ndarray,
    ki4_v: np.ndarray,
) -> np.ndarray:
    """Within-cell CP matrix for one energy group in cylindrical geometry.

    Uses the Ki₄ second-difference formula with numerical y-quadrature.
    Returns the infinite-lattice CP matrix P_inf (N, N) after applying
    the white boundary condition.
    """
    N = mesh.N
    V = mesh.volumes
    n_y = len(y_pts)

    tau = sig_t_g[:, None] * chords

    bnd_pos = np.zeros((N + 1, n_y))
    for k in range(N):
        bnd_pos[k + 1, :] = bnd_pos[k, :] + tau[k, :]

    ki4_0 = _ki4_lookup(np.zeros(n_y), ki_x, ki4_v)

    rcp = np.zeros((N, N))

    for i in range(N):
        tau_i = tau[i, :]
        sti = sig_t_g[i]
        if sti == 0:
            continue

        self_same = 2.0 * chords[i, :] - (2.0 / sti) * (
            ki4_0 - _ki4_lookup(tau_i, ki_x, ki4_v)
        )
        rcp[i, i] += 2.0 * sti * np.dot(y_wts, self_same)

        for j in range(N):
            tau_j = tau[j, :]

            if j > i:
                gap_d = bnd_pos[j, :] - bnd_pos[i + 1, :]
            elif j < i:
                gap_d = bnd_pos[i, :] - bnd_pos[j + 1, :]
            else:
                gap_d = None

            if gap_d is not None:
                gap_d = np.maximum(gap_d, 0.0)
                dd = (_ki4_lookup(gap_d, ki_x, ki4_v)
                      - _ki4_lookup(gap_d + tau_i, ki_x, ki4_v)
                      - _ki4_lookup(gap_d + tau_j, ki_x, ki4_v)
                      + _ki4_lookup(gap_d + tau_i + tau_j, ki_x, ki4_v))
            else:
                dd = np.zeros(n_y)

            gap_c = bnd_pos[i, :] + bnd_pos[j, :]
            dc = (_ki4_lookup(gap_c, ki_x, ki4_v)
                  - _ki4_lookup(gap_c + tau_i, ki_x, ki4_v)
                  - _ki4_lookup(gap_c + tau_j, ki_x, ki4_v)
                  + _ki4_lookup(gap_c + tau_i + tau_j, ki_x, ki4_v))

            rcp[i, j] += 2.0 * np.dot(y_wts, dd + dc)

    P_cell = np.zeros((N, N))
    for i in range(N):
        if sig_t_g[i] * V[i] > 0:
            P_cell[i, :] = rcp[i, :] / (sig_t_g[i] * V[i])

    P_out = np.maximum(1.0 - P_cell.sum(axis=1), 0.0)
    S_cell = mesh.surfaces[-1]  # 2*pi*r_cell for cylindrical
    P_in = sig_t_g * V * P_out / S_cell
    P_inout = max(1.0 - P_in.sum(), 0.0)

    P_inf = P_cell.copy()
    if P_inout < 1.0:
        P_inf += np.outer(P_out, P_in) / (1.0 - P_inout)

    return P_inf


# ═══════════════════════════════════════════════════════════════════════
# Public API: thin wrappers
# ═══════════════════════════════════════════════════════════════════════

def solve_cp_slab(
    materials: dict[int, Mixture],
    mesh: Mesh1D | None = None,
    max_outer: int = 500,
    keff_tol: float = 1e-6,
    flux_tol: float = 1e-5,
) -> CPResult:
    """Solve the slab collision probability eigenvalue problem.

    Parameters
    ----------
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID.
    mesh : Mesh1D, optional
        Cartesian 1-D mesh (half-cell).  Defaults to a standard PWR
        slab half-cell via :func:`geometry.factories.pwr_slab_half_cell`.
    """
    t_start = time.perf_counter()

    if mesh is None:
        from geometry import pwr_slab_half_cell
        mesh = pwr_slab_half_cell()

    if mesh.coord != CoordSystem.CARTESIAN:
        raise ValueError(
            f"CP slab solver requires CARTESIAN mesh, got {mesh.coord}"
        )

    _any_mat = next(iter(materials.values()))
    eg = _any_mat.eg
    ng = _any_mat.ng
    N = mesh.N

    xs = assemble_cell_xs(materials, mesh.mat_ids)

    # Build CP matrices using slab E₃ kernel
    print(f"  Computing slab CP matrices for {ng} groups, {N} regions ...")
    P_inf = np.empty((N, N, ng))
    for g in range(ng):
        P_inf[:, :, g] = _compute_slab_cp_group(xs.sig_t[:, g], mesh)
        if (g + 1) % 100 == 0:
            print(f"    group {g + 1}/{ng}")
    print("  CP matrices done.")

    # Eigenvalue solve
    print("  Starting power iteration ...")
    solver = CPSolver(P_inf, xs, mesh.volumes, mesh.mat_ids, materials,
                      keff_tol=keff_tol, flux_tol=flux_tol)
    keff, keff_history, phi = power_iteration(solver, max_iter=max_outer)

    flux_fuel, flux_clad, flux_cool = _volume_averaged_fluxes(
        phi, mesh.volumes, mesh.mat_ids)

    elapsed = time.perf_counter() - t_start
    print(f"  Elapsed: {elapsed:.1f}s")

    return CPResult(
        keff=keff, keff_history=keff_history, flux=phi,
        flux_fuel=flux_fuel, flux_clad=flux_clad, flux_cool=flux_cool,
        geometry=mesh, eg=eg, elapsed_seconds=elapsed,
    )


def solve_cp_concentric(
    materials: dict[int, Mixture],
    mesh: Mesh1D | None = None,
    params: CPParams | None = None,
) -> CPResult:
    """Solve the cylindrical collision probability eigenvalue problem.

    Parameters
    ----------
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID.
    mesh : Mesh1D, optional
        Cylindrical 1-D mesh (Wigner-Seitz equivalent cell).
        Defaults to a standard PWR pin via
        :func:`geometry.factories.pwr_pin_equivalent`.
    params : CPParams, optional
        Solver parameters (Ki table size, quadrature order, etc.).
    """
    t_start = time.perf_counter()

    if mesh is None:
        from geometry import pwr_pin_equivalent
        mesh = pwr_pin_equivalent()
    if params is None:
        params = CPParams()

    if mesh.coord != CoordSystem.CYLINDRICAL:
        raise ValueError(
            f"CP concentric solver requires CYLINDRICAL mesh, got {mesh.coord}"
        )

    _any_mat = next(iter(materials.values()))
    eg = _any_mat.eg
    ng = _any_mat.ng
    N = mesh.N

    print("  Building Ki3/Ki4 lookup tables ...")
    ki_x, _, ki4_v = _build_ki_tables(params.n_ki_table, params.ki_max)

    xs = assemble_cell_xs(materials, mesh.mat_ids)

    # y-quadrature (composite Gauss-Legendre over annular breakpoints)
    radii = mesh.edges[1:]  # outer radius of each annulus
    breakpoints = mesh.edges  # includes 0 at start
    gl_pts, gl_wts = np.polynomial.legendre.leggauss(params.n_quad_y)
    y_all, w_all = [], []
    for seg in range(len(breakpoints) - 1):
        a, b = breakpoints[seg], breakpoints[seg + 1]
        y_all.append(0.5 * (b - a) * gl_pts + 0.5 * (b + a))
        w_all.append(0.5 * (b - a) * gl_wts)
    y_pts = np.concatenate(y_all)
    y_wts = np.concatenate(w_all)

    chords = _chord_half_lengths(radii, y_pts)

    # Build CP matrices using cylindrical Ki₃/Ki₄ kernel
    print(f"  Computing CP matrices for {ng} groups, {N} regions ...")
    P_inf = np.empty((N, N, ng))
    for g in range(ng):
        P_inf[:, :, g] = _compute_cp_group(
            xs.sig_t[:, g], mesh, chords, y_pts, y_wts, ki_x, ki4_v,
        )
        if (g + 1) % 100 == 0:
            print(f"    group {g + 1}/{ng}")
    print("  CP matrices done.")

    # Eigenvalue solve
    print("  Starting power iteration ...")
    solver = CPSolver(P_inf, xs, mesh.volumes, mesh.mat_ids, materials,
                      keff_tol=params.keff_tol, flux_tol=params.flux_tol)
    keff, keff_history, phi = power_iteration(solver, max_iter=params.max_outer)

    flux_fuel, flux_clad, flux_cool = _volume_averaged_fluxes(
        phi, mesh.volumes, mesh.mat_ids)

    elapsed = time.perf_counter() - t_start
    print(f"  Elapsed: {elapsed:.1f}s")

    return CPResult(
        keff=keff, keff_history=keff_history, flux=phi,
        flux_fuel=flux_fuel, flux_clad=flux_clad, flux_cool=flux_cool,
        geometry=mesh, eg=eg, elapsed_seconds=elapsed,
    )
