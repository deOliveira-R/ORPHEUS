"""1D Discrete Ordinates (SN) neutron transport solver.

Solves the multi-group neutron transport equation on a 1D slab mesh
with reflective boundary conditions using Gauss-Legendre angular
quadrature and diamond-difference spatial discretization.

The 1D transport equation for direction mu_n:
    mu_n * dpsi/dx + Sigma_t * psi = Q / 2

where the /2 is the 1D isotropic source normalization (integral over
[-1,1] of (1/2) dmu = 1).
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.sparse import issparse

from data.macro_xs.mixture import Mixture


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GaussLegendreQuadrature:
    """Gauss-Legendre angular quadrature on [-1, 1].

    Weights sum to 2.0 (the measure of [-1, 1]).
    Points are symmetric about mu=0; the partner of index i is N-1-i.
    """

    mu: np.ndarray       # (N,) quadrature points
    weights: np.ndarray  # (N,) weights summing to 2.0
    N: int

    @classmethod
    def gauss_legendre(cls, n_ordinates: int) -> GaussLegendreQuadrature:
        """Build N-point Gauss-Legendre quadrature.

        Parameters
        ----------
        n_ordinates : int
            Number of quadrature points (must be even for SN).
        """
        mu, w = np.polynomial.legendre.leggauss(n_ordinates)
        return cls(mu=mu, weights=w, N=n_ordinates)


@dataclass
class Slab1DGeometry:
    """1D slab geometry for SN transport.

    A sequence of cells with specified widths and material IDs.
    Boundary conditions are reflective on both sides (infinite lattice).
    """

    cell_widths: np.ndarray   # (N_cells,) cell widths in cm
    mat_ids: np.ndarray       # (N_cells,) material ID per cell
    N: int                    # number of cells

    @classmethod
    def from_benchmark(
        cls,
        n_fuel: int,
        n_mod: int,
        t_fuel: float,
        t_mod: float,
    ) -> Slab1DGeometry:
        """Build a half-cell slab [fuel | moderator] for benchmark comparison.

        Reflective BCs on both sides give the infinite-lattice eigenvalue.
        The half-cell is [fuel(t_fuel) | mod(t_mod)], matching the CP
        benchmark geometry.

        Parameters
        ----------
        n_fuel : int — number of cells in fuel region
        n_mod  : int — number of cells in moderator region
        t_fuel : float — total fuel thickness (cm)
        t_mod  : float — total moderator thickness (cm)
        """
        dx_fuel = t_fuel / n_fuel
        dx_mod = t_mod / n_mod
        widths = np.concatenate([
            np.full(n_fuel, dx_fuel),
            np.full(n_mod, dx_mod),
        ])
        mat_ids = np.concatenate([
            np.full(n_fuel, 2, dtype=int),   # fuel = mat 2
            np.full(n_mod, 0, dtype=int),    # mod  = mat 0
        ])
        return cls(cell_widths=widths, mat_ids=mat_ids, N=n_fuel + n_mod)

    @classmethod
    def homogeneous(
        cls,
        n_cells: int,
        total_width: float,
        mat_id: int = 2,
    ) -> Slab1DGeometry:
        """Build a homogeneous slab of a single material."""
        dx = total_width / n_cells
        return cls(
            cell_widths=np.full(n_cells, dx),
            mat_ids=np.full(n_cells, mat_id, dtype=int),
            N=n_cells,
        )


@dataclass
class SN1DResult:
    """Results of a 1D SN calculation."""

    keff: float
    keff_history: list[float]
    flux: np.ndarray          # (N_cells, ng) scalar flux per cell
    geometry: Slab1DGeometry
    eg: np.ndarray            # (ng+1,) energy group boundaries
    elapsed_seconds: float


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

def solve_sn_1d(
    materials: dict[int, Mixture],
    geom: Slab1DGeometry,
    quad: GaussLegendreQuadrature | None = None,
    max_outer: int = 500,
    keff_tol: float = 1e-7,
    flux_tol: float = 1e-6,
    max_inner: int = 200,
    inner_tol: float = 1e-8,
) -> SN1DResult:
    """Solve the 1D multi-group SN eigenvalue problem.

    Uses source iteration with transport sweeps (diamond differencing).
    Reflective BCs on both sides (infinite lattice).

    Parameters
    ----------
    materials : dict mapping material ID to Mixture.
    geom : Slab1DGeometry — cell widths and material assignment.
    quad : GaussLegendreQuadrature (default: S16, 16 ordinates).
    max_outer : int — maximum power iterations.
    keff_tol : float — convergence tolerance on keff.
    flux_tol : float — convergence tolerance on scalar flux.
    max_inner : int — maximum scattering source iterations per outer.
    inner_tol : float — convergence tolerance for inner iterations.
    """
    t_start = time.perf_counter()

    if quad is None:
        quad = GaussLegendreQuadrature.gauss_legendre(16)

    _any_mat = next(iter(materials.values()))
    eg = _any_mat.eg
    ng = _any_mat.ng
    nc = geom.N

    # --- Pre-compute per-cell cross sections ---
    sig_t = np.empty((nc, ng))      # total XS
    sig_s0 = np.empty((nc, ng, ng)) # P0 scattering matrix (dense)
    nu_sig_f = np.empty((nc, ng))   # production XS
    chi = np.empty((nc, ng))        # fission spectrum
    sig_a = np.empty((nc, ng))      # absorption (for keff balance)

    for i in range(nc):
        m = materials[geom.mat_ids[i]]
        sig_t[i, :] = m.SigT
        nu_sig_f[i, :] = m.SigP
        chi[i, :] = m.chi

        # Dense scattering matrix
        S0 = m.SigS[0]
        sig_s0[i, :, :] = S0.toarray() if issparse(S0) else np.asarray(S0)

        # Absorption = total - scattering out
        sig_a[i, :] = m.SigT - np.asarray(m.SigS[0].sum(axis=1)).ravel()

    dx = geom.cell_widths

    # Pre-compute sweep coefficients per cell and per ordinate (positive mu only)
    # For mu > 0 (half of the ordinates, stored in the upper half of leggauss output):
    n_half = quad.N // 2
    mu_pos = quad.mu[n_half:]          # (n_half,) positive mu values
    w_pos = quad.weights[n_half:]      # (n_half,) corresponding weights

    # Coefficients for diamond difference:
    #   psi_out = alpha * psi_in + beta * Q
    #   where alpha = (2*mu - dx*sig_t) / (2*mu + dx*sig_t)   [NO — this is step]
    #   Diamond: psi_out = (2*mu*psi_in + dx*Q/2) / (2*mu + dx*sig_t)
    # Pre-compute denominator and numerator coefficients
    # denom[n, i, g] = 2*mu[n] + dx[i]*sig_t[i, g]
    # For each (ordinate, cell): psi_out = (2*mu*psi_in + dx*Q/2) / denom
    two_mu_pos = 2.0 * mu_pos                   # (n_half,)
    half_dx = 0.5 * dx                          # (nc,)
    # denom[n, i, g] = two_mu_pos[n] + dx[i]*sig_t[i,g]
    # Shape broadcast: (n_half, 1, 1) + (1, nc, ng)
    denom = two_mu_pos[:, None, None] + dx[None, :, None] * sig_t[None, :, :]  # (n_half, nc, ng)
    source_coeff = half_dx[None, :, None] / denom   # (n_half, nc, ng) -- multiply by Q to get source term
    stream_coeff = two_mu_pos[:, None, None] / denom  # (n_half, nc, ng) -- multiply by psi_in

    # --- Power iteration ---
    phi = np.ones((nc, ng))
    keff = 1.0
    keff_history: list[float] = []

    # Persistent boundary fluxes (reflective BC state)
    # psi_bc[n, g]: outgoing angular flux at the right boundary for positive direction n
    psi_bc_right = np.zeros((n_half, ng))
    psi_bc_left = np.zeros((n_half, ng))

    for n_outer in range(1, max_outer + 1):
        phi_old = phi.copy()
        keff_old = keff

        # Fission source
        fission_rate = np.sum(nu_sig_f * phi, axis=1)  # (nc,)
        Q_f = chi * fission_rate[:, np.newaxis] / keff  # (nc, ng)

        # --- Inner iterations (scattering + transport sweep) ---
        for n_inner in range(max_inner):
            phi_prev = phi.copy()

            # Total source = fission + scattering
            Q_s = np.einsum('ijk,ij->ik', sig_s0, phi)  # (nc, ng)
            Q = Q_f + Q_s  # (nc, ng)

            # Transport sweep with reflective BCs
            phi = _transport_sweep_fast(
                Q, stream_coeff, source_coeff, w_pos, n_half,
                psi_bc_right, psi_bc_left, nc, ng,
            )

            # Check inner convergence
            norm = np.linalg.norm(phi)
            if norm > 0:
                inner_res = np.linalg.norm(phi - phi_prev) / norm
                if inner_res < inner_tol:
                    break

        # Update keff
        production = np.sum(nu_sig_f * phi * dx[:, np.newaxis])
        absorption = np.sum(sig_a * phi * dx[:, np.newaxis])
        keff = production / absorption

        keff_history.append(keff)

        # Check convergence
        keff_change = abs(keff - keff_old)
        flux_change = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)

        if keff_change < keff_tol and flux_change < flux_tol and n_outer > 2:
            break

    elapsed = time.perf_counter() - t_start

    return SN1DResult(
        keff=keff,
        keff_history=keff_history,
        flux=phi,
        geometry=geom,
        eg=eg,
        elapsed_seconds=elapsed,
    )


def _transport_sweep_fast(
    Q: np.ndarray,
    stream_coeff: np.ndarray,
    source_coeff: np.ndarray,
    w_pos: np.ndarray,
    n_half: int,
    psi_bc_right: np.ndarray,
    psi_bc_left: np.ndarray,
    nc: int,
    ng: int,
) -> np.ndarray:
    """Vectorized diamond-difference transport sweep.

    Exploits GL symmetry: mu[N-1-n] = -mu[n], w[N-1-n] = w[n].
    Only sweeps positive directions; negative directions are the
    reverse sweep with the same |mu| and weight.

    The diamond-difference recurrence psi_out[i] = a[i]*psi_in[i] + b[i]*Q[i]
    is solved via cumulative products to avoid Python cell loops.

    Parameters
    ----------
    Q : (nc, ng) total source
    stream_coeff : (n_half, nc, ng) = 2*mu / (2*mu + dx*sig_t)
    source_coeff : (n_half, nc, ng) = (dx/2) / (2*mu + dx*sig_t)
    w_pos : (n_half,) weights for positive directions
    psi_bc_right, psi_bc_left : (n_half, ng) boundary flux storage (updated in-place)

    Returns new scalar flux phi (nc, ng).
    """
    phi_new = np.zeros((nc, ng))

    # Source terms b[n, i, g] = source_coeff * Q  (precompute once)
    bQ = source_coeff * Q[np.newaxis, :, :]  # (n_half, nc, ng)

    for n in range(n_half):
        w = w_pos[n]
        a = stream_coeff[n]   # (nc, ng)
        s = bQ[n]             # (nc, ng)

        # --- Forward sweep (mu > 0): left to right ---
        psi_fwd = _solve_recurrence_fwd(a, s, psi_bc_left[n])  # (nc, ng) cell-avg
        psi_bc_right[n, :] = _outgoing_fwd(a, s, psi_bc_left[n], nc)

        phi_new += w * psi_fwd

        # --- Backward sweep (mu < 0): right to left ---
        # Reverse arrays, sweep, reverse back
        a_rev = a[::-1]
        s_rev = s[::-1]
        psi_bwd_rev = _solve_recurrence_fwd(a_rev, s_rev, psi_bc_right[n])
        psi_bc_left[n, :] = _outgoing_fwd(a_rev, s_rev, psi_bc_right[n], nc)

        phi_new += w * psi_bwd_rev[::-1]

    return phi_new


def _solve_recurrence_fwd(
    a: np.ndarray, s: np.ndarray, psi0: np.ndarray,
) -> np.ndarray:
    """Solve diamond-difference recurrence and return cell-average fluxes.

    Recurrence: psi_out[i] = a[i]*psi_in[i] + s[i]
                psi_in[0] = psi0
                psi_in[i+1] = psi_out[i]
                psi_avg[i] = 0.5*(psi_in[i] + psi_out[i])

    Uses cumulative products for vectorization:
        psi_in[i] = (prod_{k=0}^{i-1} a[k]) * psi0
                   + sum_{j=0}^{i-1} (prod_{k=j+1}^{i-1} a[k]) * s[j]

    Parameters
    ----------
    a : (nc, ng) streaming coefficient
    s : (nc, ng) source coefficient (already multiplied by Q)
    psi0 : (ng,) incoming boundary flux

    Returns (nc, ng) cell-average angular flux.
    """
    nc, ng = a.shape

    # Cumulative product of a: cp[i] = prod_{k=0}^{i} a[k]
    cp = np.cumprod(a, axis=0)  # (nc, ng)

    # psi_in[0] = psi0
    # psi_in[i] = cp[i-1] * psi0 + sum_{j=0}^{i-1} (cp[i-1]/cp[j]) * s[j]
    #           = cp[i-1] * (psi0 + sum_{j=0}^{i-1} s[j]/cp[j])
    # Note: cp[j] can underflow for thick cells; use safe division.

    # s_over_cp[j] = s[j] / cp[j]
    s_over_cp = s / cp  # (nc, ng)

    # cumsum of s_over_cp: cs[i] = sum_{j=0}^{i} s[j]/cp[j]
    cs = np.cumsum(s_over_cp, axis=0)  # (nc, ng)

    # psi_in[0] = psi0
    # psi_in[i] = cp[i-1] * (psi0 + cs[i-1])  for i >= 1
    psi_in = np.empty((nc, ng))
    psi_in[0] = psi0
    if nc > 1:
        psi_in[1:] = cp[:-1] * (psi0[np.newaxis, :] + cs[:-1])

    # psi_out[i] = a[i]*psi_in[i] + s[i]
    psi_out = a * psi_in + s

    # Cell-average
    return 0.5 * (psi_in + psi_out)


def _outgoing_fwd(
    a: np.ndarray, s: np.ndarray, psi0: np.ndarray, nc: int,
) -> np.ndarray:
    """Compute the outgoing flux at the end of a forward sweep."""
    cp = np.cumprod(a, axis=0)
    s_over_cp = s / cp
    cs = np.cumsum(s_over_cp, axis=0)
    return cp[-1] * (psi0 + cs[-1])
