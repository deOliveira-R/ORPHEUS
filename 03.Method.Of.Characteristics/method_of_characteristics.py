"""2D Method of Characteristics (MoC) neutron transport solver for a PWR pin cell.

Solves the multi-group transport equation on a 2D Cartesian mesh with
reflective boundary conditions using 8 characteristic ray directions
(E, NE, N, NW, W, SW, S, SE) and direct transport sweeps.

Port of MATLAB ``MoC_PWR.m`` + ``fly.m`` + ``flyFrom.m``.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from data.macro_xs.cell_xs import assemble_cell_xs
from data.macro_xs.mixture import Mixture
from numerics.eigenvalue import power_iteration


# ---------------------------------------------------------------------------
# Direction enumeration (matching MATLAB d.EAST=1..d.SOUTH_EAST=8)
# ---------------------------------------------------------------------------

class Dir(IntEnum):
    EAST = 0
    NORTH_EAST = 1
    NORTH = 2
    NORTH_WEST = 3
    WEST = 4
    SOUTH_WEST = 5
    SOUTH = 6
    SOUTH_EAST = 7

N_RAYS = 8
_SQRT2 = np.sqrt(2.0)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MoCGeometry:
    """2D Cartesian mesh for a PWR pin cell (full quarter-symmetry)."""

    n_cells: int
    delta: float  # mesh step (cm)
    mat_map: np.ndarray  # (n_cells, n_cells) int — 0=cool, 1=clad, 2=fuel
    volume: np.ndarray  # (n_cells, n_cells) cell volumes (cm^2)

    @classmethod
    def default_pwr(cls) -> MoCGeometry:
        """Standard 10x10 cell mesh: 5 fuel + 1 clad + 4 coolant columns."""
        n = 10
        delta = 0.2

        vol = np.full((n, n), delta**2)
        vol[0, :] /= 2
        vol[-1, :] /= 2
        vol[:, 0] /= 2
        vol[:, -1] /= 2

        # Each row is identical: fuel(5) + clad(1) + cool(4)
        row = np.array([2, 2, 2, 2, 2, 1, 0, 0, 0, 0], dtype=int)
        mat = np.tile(row, (n, 1)).T  # shape (10, 10)

        return cls(n_cells=n, delta=delta, mat_map=mat, volume=vol)


@dataclass
class MoCResult:
    """Results of a Method of Characteristics calculation."""

    keff: float
    keff_history: list[float]
    flux_fuel: np.ndarray     # (ng,) volume-averaged scalar flux in fuel
    flux_clad: np.ndarray     # (ng,) volume-averaged scalar flux in clad
    flux_cool: np.ndarray     # (ng,) volume-averaged scalar flux in coolant
    scalar_flux: np.ndarray   # (n_cells, n_cells, ng) scalar flux in each cell
    geometry: MoCGeometry
    eg: np.ndarray            # (ng+1,) energy group boundaries
    elapsed_seconds: float


# ---------------------------------------------------------------------------
# Core transport functions
# ---------------------------------------------------------------------------

def _fly(fi_a: np.ndarray, sig: np.ndarray, s: float, q: np.ndarray) -> np.ndarray:
    """Transport a ray across a half-node.

    fiB = fiA * exp(-sig*s) + q * (1 - exp(-sig*s)) / sig

    Port of MATLAB ``fly.m``.
    """
    tau = sig * s
    exp_tau = np.exp(-tau)
    return fi_a * exp_tau + q * (1.0 - exp_tau) / sig


def _fly_from(
    fi: np.ndarray,
    sig_t: np.ndarray,
    q: np.ndarray,
    ix: int, iy: int,
    direction: Dir,
    n: int,
    delta: float,
) -> None:
    """Propagate angular flux from node (ix,iy) to neighbor along direction.

    Modifies fi in-place. Port of MATLAB ``flyFrom.m``.

    Parameters
    ----------
    fi : (n_cells, n_cells, ng, N_RAYS) angular flux array (modified in-place).
    sig_t : (n_cells, n_cells, ng) total cross section.
    q : (n_cells, n_cells, ng) isotropic source per ray.
    """
    d = direction

    # Step offsets and path length for each direction
    if d == Dir.EAST:
        if ix == n - 1:
            return
        dx, dy, s = 1, 0, delta / 2
    elif d == Dir.NORTH_EAST:
        if ix == n - 1 or iy == 0:
            return
        dx, dy, s = 1, -1, delta / _SQRT2
    elif d == Dir.NORTH:
        if iy == 0:
            return
        dx, dy, s = 0, -1, delta / 2
    elif d == Dir.NORTH_WEST:
        if ix == 0 or iy == 0:
            return
        dx, dy, s = -1, -1, delta / _SQRT2
    elif d == Dir.WEST:
        if ix == 0:
            return
        dx, dy, s = -1, 0, delta / 2
    elif d == Dir.SOUTH_WEST:
        if ix == 0 or iy == n - 1:
            return
        dx, dy, s = -1, 1, delta / _SQRT2
    elif d == Dir.SOUTH:
        if iy == n - 1:
            return
        dx, dy, s = 0, 1, delta / 2
    elif d == Dir.SOUTH_EAST:
        if ix == n - 1 or iy == n - 1:
            return
        dx, dy, s = 1, 1, delta / _SQRT2
    else:
        return

    jx, jy = ix + dx, iy + dy
    di = int(d)

    # Fly from node center to boundary, then from boundary to next node center
    fi_mid = _fly(fi[ix, iy, :, di], sig_t[ix, iy, :], s, q[ix, iy, :])
    fi[jx, jy, :, di] = _fly(fi_mid, sig_t[jx, jy, :], s, q[jx, jy, :])


# ---------------------------------------------------------------------------
# MoC eigenvalue solver (EigenvalueSolver protocol)
# ---------------------------------------------------------------------------

class MoCSolver:
    """2D Method of Characteristics solver satisfying the EigenvalueSolver protocol."""

    def __init__(
        self,
        materials: dict[int, Mixture],
        geom: MoCGeometry,
        keff_tol: float = 1e-6,
        flux_tol: float = 1e-5,
    ) -> None:
        self.geom = geom
        self.materials = materials
        self.keff_tol = keff_tol
        self.flux_tol = flux_tol

        n = geom.n_cells
        ng = next(iter(materials.values())).ng
        self.n = n
        self.ng = ng

        # Per-cell cross sections via assemble_cell_xs (flattened, then reshaped)
        xs = assemble_cell_xs(materials, geom.mat_map)
        self.sig_t = xs.sig_t.reshape(n, n, ng)
        self.sig_a = xs.sig_a.reshape(n, n, ng)
        self.sig_p = xs.sig_p.reshape(n, n, ng)
        self.chi = xs.chi.reshape(n, n, ng)

        # Sparse scattering and (n,2n) matrices per material (only 3 unique)
        self.sig_s0_mat: dict[int, np.ndarray] = {}
        self.sig2_mat: dict[int, np.ndarray] = {}
        for mat_id, mix in materials.items():
            self.sig_s0_mat[mat_id] = mix.SigS[0]
            self.sig2_mat[mat_id] = mix.Sig2

    def initial_flux_distribution(self) -> np.ndarray:
        """Return ones array with shape (n, n, ng, N_RAYS)."""
        return np.ones((self.n, self.n, self.ng, N_RAYS))

    def compute_fission_source(
        self,
        flux_distribution: np.ndarray,
        keff: float,
    ) -> np.ndarray:
        """Build the fission source: chi * (SigP . phi) / keff.

        Also stashes the production rate for the per-iteration flux
        normalization applied at the end of solve_fixed_source.

        Returns shape (n, n, ng) — the fission component of the isotropic source.
        """
        scalar_flux = flux_distribution.sum(axis=3)  # (n, n, ng)
        n = self.n
        geom = self.geom
        fission_source = np.empty((n, n, self.ng))

        # Accumulate production rate for normalization (same flux as source)
        p_rate = 0.0
        for iy in range(n):
            for ix in range(n):
                v = geom.volume[ix, iy]
                FI = scalar_flux[ix, iy, :]
                mat_id = geom.mat_map[ix, iy]
                sig2_cs = np.array(self.sig2_mat[mat_id].sum(axis=1)).ravel()
                p_rate += (self.sig_p[ix, iy, :] + 2 * sig2_cs) @ FI * v
                fission_source[ix, iy, :] = (
                    self.chi[ix, iy, :] * (self.sig_p[ix, iy, :] @ FI) / keff
                )
        self._p_rate_for_norm = p_rate
        return fission_source

    def solve_fixed_source(
        self,
        fission_source: np.ndarray,
        flux_distribution: np.ndarray,
    ) -> np.ndarray:
        """Ray-tracing sweep over 8 directions with scattering + (n,2n) sources.

        Assembles the total isotropic source q = (fission + 2*Sig2^T + SigS0^T) * phi / N_RAYS,
        then performs the characteristic sweeps with reflective BCs.

        Returns updated angular flux (n, n, ng, N_RAYS).
        """
        n = self.n
        geom = self.geom
        fi = flux_distribution  # modified in-place by _fly_from
        scalar_flux = fi.sum(axis=3)  # (n, n, ng)

        # Isotropic source per ray: (fission + scattering + n2n) / N_RAYS
        q = np.empty((n, n, self.ng))
        for iy in range(n):
            for ix in range(n):
                FI = scalar_flux[ix, iy, :]
                mat_id = geom.mat_map[ix, iy]
                q[ix, iy, :] = (
                    fission_source[ix, iy, :]
                    + 2.0 * (self.sig2_mat[mat_id].T @ FI)
                    + self.sig_s0_mat[mat_id].T @ FI
                ) / N_RAYS

        # --- Ray sweeps with reflective boundary conditions ---

        # Horizontal: EAST then reflect, WEST then reflect
        for iy in range(n):
            for ix in range(n):
                _fly_from(fi, self.sig_t, q, ix, iy, Dir.EAST, n, geom.delta)
            fi[n - 1, iy, :, Dir.WEST] = fi[n - 1, iy, :, Dir.EAST]
            for ix in range(n - 1, -1, -1):
                _fly_from(fi, self.sig_t, q, ix, iy, Dir.WEST, n, geom.delta)
            fi[0, iy, :, Dir.EAST] = fi[0, iy, :, Dir.WEST]

        # Vertical: SOUTH then reflect, NORTH then reflect
        for ix in range(n):
            for iy in range(n):
                _fly_from(fi, self.sig_t, q, ix, iy, Dir.SOUTH, n, geom.delta)
            fi[ix, n - 1, :, Dir.NORTH] = fi[ix, n - 1, :, Dir.SOUTH]
            for iy in range(n - 1, -1, -1):
                _fly_from(fi, self.sig_t, q, ix, iy, Dir.NORTH, n, geom.delta)
            fi[ix, 0, :, Dir.SOUTH] = fi[ix, 0, :, Dir.NORTH]

        # Main diagonal SE, reflect at corner, NW back, reflect at corner
        ixx, iyy = 0, 0
        while ixx < n - 1:
            _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_EAST, n, geom.delta)
            ixx += 1
            iyy += 1
        # Reflect at bottom-right corner
        fi[n - 1, n - 1, :, Dir.NORTH_WEST] = fi[n - 1, n - 1, :, Dir.SOUTH_EAST]
        fi[n - 1, n - 1, :, Dir.NORTH_EAST] = fi[n - 1, n - 1, :, Dir.SOUTH_EAST]
        fi[n - 1, n - 1, :, Dir.SOUTH_WEST] = fi[n - 1, n - 1, :, Dir.SOUTH_EAST]
        ixx, iyy = n - 1, n - 1
        while ixx > 0:
            _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_WEST, n, geom.delta)
            ixx -= 1
            iyy -= 1
        # Reflect at top-left corner
        fi[0, 0, :, Dir.SOUTH_EAST] = fi[0, 0, :, Dir.NORTH_WEST]
        fi[0, 0, :, Dir.NORTH_EAST] = fi[0, 0, :, Dir.NORTH_WEST]
        fi[0, 0, :, Dir.SOUTH_WEST] = fi[0, 0, :, Dir.NORTH_WEST]

        # Anti-diagonal SW, reflect at corner, NE back, reflect at corner
        ixx, iyy = n - 1, 0
        while ixx > 0:
            _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_WEST, n, geom.delta)
            ixx -= 1
            iyy += 1
        # Reflect at bottom-left corner
        fi[0, n - 1, :, Dir.NORTH_EAST] = fi[0, n - 1, :, Dir.SOUTH_WEST]
        fi[0, n - 1, :, Dir.NORTH_WEST] = fi[0, n - 1, :, Dir.SOUTH_WEST]
        fi[0, n - 1, :, Dir.SOUTH_EAST] = fi[0, n - 1, :, Dir.SOUTH_WEST]
        ixx, iyy = 0, n - 1
        while ixx < n - 1:
            _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_EAST, n, geom.delta)
            ixx += 1
            iyy -= 1
        # Reflect at top-right corner
        fi[n - 1, 0, :, Dir.SOUTH_WEST] = fi[n - 1, 0, :, Dir.NORTH_EAST]
        fi[n - 1, 0, :, Dir.SOUTH_EAST] = fi[n - 1, 0, :, Dir.NORTH_EAST]
        fi[n - 1, 0, :, Dir.NORTH_WEST] = fi[n - 1, 0, :, Dir.NORTH_EAST]

        # Off-diagonal SE rays, reflect on right boundary
        for start_ix in range(n - 2, 0, -1):
            ixx, iyy = start_ix, 0
            while ixx < n - 1:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_EAST, n, geom.delta)
                ixx += 1
                iyy += 1
            fi[ixx, iyy, :, Dir.SOUTH_WEST] = fi[ixx, iyy, :, Dir.SOUTH_EAST]

        # Off-diagonal SE rays, reflect on bottom boundary
        for start_iy in range(1, n - 1):
            ixx, iyy = 0, start_iy
            while iyy < n - 1:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_EAST, n, geom.delta)
                ixx += 1
                iyy += 1
            fi[ixx, iyy, :, Dir.NORTH_EAST] = fi[ixx, iyy, :, Dir.SOUTH_EAST]

        # Off-diagonal SW rays, reflect on left boundary
        for start_ix in range(1, n - 1):
            ixx, iyy = start_ix, 0
            while ixx > 0:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_WEST, n, geom.delta)
                ixx -= 1
                iyy += 1
            fi[ixx, iyy, :, Dir.SOUTH_EAST] = fi[ixx, iyy, :, Dir.SOUTH_WEST]

        # Off-diagonal SW rays, reflect on bottom boundary
        for start_iy in range(1, n - 1):
            ixx, iyy = n - 1, start_iy
            while iyy < n - 1:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.SOUTH_WEST, n, geom.delta)
                ixx -= 1
                iyy += 1
            fi[ixx, iyy, :, Dir.NORTH_WEST] = fi[ixx, iyy, :, Dir.SOUTH_WEST]

        # Off-diagonal NW rays, reflect on left boundary
        for start_ix in range(1, n - 1):
            ixx, iyy = start_ix, n - 1
            while ixx > 0:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_WEST, n, geom.delta)
                ixx -= 1
                iyy -= 1
            fi[ixx, iyy, :, Dir.NORTH_EAST] = fi[ixx, iyy, :, Dir.NORTH_WEST]

        # Off-diagonal NW rays, reflect on top boundary
        for start_iy in range(n - 2, 0, -1):
            ixx, iyy = n - 1, start_iy
            while iyy > 0:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_WEST, n, geom.delta)
                ixx -= 1
                iyy -= 1
            fi[ixx, iyy, :, Dir.SOUTH_WEST] = fi[ixx, iyy, :, Dir.NORTH_WEST]

        # Off-diagonal NE rays, reflect on top boundary
        for start_iy in range(1, n - 1):
            ixx, iyy = 0, start_iy
            while iyy > 0:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_EAST, n, geom.delta)
                ixx += 1
                iyy -= 1
            fi[ixx, iyy, :, Dir.SOUTH_EAST] = fi[ixx, iyy, :, Dir.NORTH_EAST]

        # Off-diagonal NE rays, reflect on right boundary
        for start_ix in range(1, n - 1):
            ixx, iyy = start_ix, n - 1
            while ixx < n - 1:
                _fly_from(fi, self.sig_t, q, ixx, iyy, Dir.NORTH_EAST, n, geom.delta)
                ixx += 1
                iyy -= 1
            fi[ixx, iyy, :, Dir.NORTH_WEST] = fi[ixx, iyy, :, Dir.NORTH_EAST]

        # Normalize flux so total production rate = 100 n/s (prevents overflow)
        fi *= 100.0 / self._p_rate_for_norm

        return fi

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        """Compute k_eff = production / absorption from scalar flux."""
        n = self.n
        geom = self.geom
        scalar_flux = flux_distribution.sum(axis=3)  # (n, n, ng)

        p_rate = 0.0
        a_rate = 0.0
        for iy in range(n):
            for ix in range(n):
                v = geom.volume[ix, iy]
                FI = scalar_flux[ix, iy, :]
                mat_id = geom.mat_map[ix, iy]
                sig2_cs = np.array(self.sig2_mat[mat_id].sum(axis=1)).ravel()
                p_rate += (self.sig_p[ix, iy, :] + 2 * sig2_cs) @ FI * v
                a_rate += self.sig_a[ix, iy, :] @ FI * v

        keff = p_rate / a_rate
        print(f"  keff = {keff:9.5f}")
        return keff

    def converged(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: np.ndarray,
        flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        """Check convergence of eigenvalue and flux."""
        if iteration <= 2:
            return False
        keff_change = abs(keff - keff_old)
        flux_change = np.linalg.norm(flux_distribution - flux_old) / \
            max(np.linalg.norm(flux_distribution), 1e-30)
        return keff_change < self.keff_tol and flux_change < self.flux_tol


# ---------------------------------------------------------------------------
# Public wrapper
# ---------------------------------------------------------------------------

def solve_moc(
    materials: dict[int, Mixture],
    geom: MoCGeometry | None = None,
    max_outer: int = 200,
) -> MoCResult:
    """Run the 2D Method of Characteristics transport calculation.

    Parameters
    ----------
    materials : dict mapping material ID (0=cool, 1=clad, 2=fuel) to Mixture.
    geom : MoCGeometry (default: standard 10x10 PWR mesh).
    max_outer : int — maximum number of outer (power) iterations.
    """
    t_start = time.perf_counter()

    if geom is None:
        geom = MoCGeometry.default_pwr()

    eg = materials[2].eg
    ng = materials[2].ng

    solver = MoCSolver(materials, geom)
    keff, keff_history, fi = power_iteration(solver, max_iter=max_outer)

    # --- Post-processing: volume-averaged spectra ---
    n = geom.n_cells
    scalar_flux = fi.sum(axis=3)  # (n, n, ng)

    vol_fuel = geom.volume[geom.mat_map == 2].sum()
    vol_clad = geom.volume[geom.mat_map == 1].sum()
    vol_cool = geom.volume[geom.mat_map == 0].sum()

    flux_fuel = np.zeros(ng)
    flux_clad = np.zeros(ng)
    flux_cool = np.zeros(ng)

    for iy in range(n):
        for ix in range(n):
            v = geom.volume[ix, iy]
            FI = scalar_flux[ix, iy, :]
            mat_id = geom.mat_map[ix, iy]
            if mat_id == 2:
                flux_fuel += FI * v / vol_fuel
            elif mat_id == 1:
                flux_clad += FI * v / vol_clad
            else:
                flux_cool += FI * v / vol_cool

    elapsed = time.perf_counter() - t_start
    print(f"  Elapsed: {elapsed:.1f}s")

    return MoCResult(
        keff=keff_history[-1],
        keff_history=keff_history,
        flux_fuel=flux_fuel,
        flux_clad=flux_clad,
        flux_cool=flux_cool,
        scalar_flux=scalar_flux,
        geometry=geom,
        eg=eg,
        elapsed_seconds=elapsed,
    )
