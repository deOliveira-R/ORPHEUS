"""Direct transport operator for Krylov inner solves.

Provides the explicit operator T: ψ → T·ψ via finite differences, used
by the ``bicgstab`` inner solver path in :class:`SNSolver`.

Three geometries are supported:

* **Cartesian 2D** — ``T = μ_x ∂/∂x + μ_y ∂/∂y + Σ_t``
* **Spherical 1D** — ``T = μ (A ∂/∂r)/V + (α ∂/∂μ)/V + Σ_t``
* **Cylindrical 1D** — per-level azimuthal redistribution

The sweep-based solver (source iteration) inverts T implicitly via
diamond-difference sweeps.  This module forms T explicitly so that
scipy's Krylov solvers (BiCGSTAB, GMRES) can solve  T·ψ = b  directly.

.. warning:: Operator-sweep inconsistency

   The operator T formed here does **not** use the same diamond-difference
   (DD) or weighted-diamond-difference (WDD) closure as the sweep paths
   in :mod:`orpheus.sn.sweep`. Instead it uses:

   * **Cartesian**: upwind cell-center finite differences for the
     streaming gradient — first-order accurate and consistent with DD
     on **uniform** meshes, but first-order inconsistent on
     non-uniform meshes (divides by the local ``dx[ix]`` rather than
     the cell-center distance ``(dx[ix]+dx[ix±1])/2``).
   * **Curvilinear**: arithmetic averages for spatial face fluxes and
     τ-weighted interpolation for angular face fluxes, rather than the
     DD closure ``ψ_out = 2·ψ_avg − ψ_in`` and the WDD closure
     ``ψ_angle_out = (ψ − (1−τ)·ψ_angle_in)/τ`` used by the sweeps.

   **Why not use DD face fluxes directly?** The DD/WDD closure turns
   the operator into a triangular system whose condition number grows
   exponentially with optical thickness. Forward-substitution (the
   sweep) is the natural solver for such systems; applying them as a
   matvec inside unpreconditioned BiCGSTAB produces catastrophically
   growing face fluxes in the Krylov search directions, leading to
   overflow. A DD-consistent Krylov solve would require a sweep-based
   preconditioner (DSA/TSA), which is a future enhancement.

   **Practical implication.** On uniform meshes, the BiCGSTAB and
   source-iteration paths converge to the same physics in the fine-mesh
   limit (both are first-order consistent with the transport equation).
   On non-uniform meshes, the BiCGSTAB Cartesian path has an additional
   O(h) inconsistency in the gradient stencil. On heterogeneous problems
   at coarse meshes, the two paths can differ at O(h²) in the
   eigenvalue. For verification work, **always use source iteration**
   (the default).

   See GitHub issues #96 (Cartesian) and #97 (curvilinear) for the
   full audit trail. The inconsistency was surfaced during the
   ERR-025 diagnosis (Phase 2.1b of the verification campaign).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse.linalg import LinearOperator

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperatorMixin,
)

from .quadrature import AngularQuadrature

if TYPE_CHECKING:
    from .geometry import SNMesh

__all__ = [
    "EquationMap",
    "build_equation_map",
    "build_equation_map_spherical",
    "build_equation_map_cylindrical",
    "solution_to_angular_flux",
    "solution_to_angular_flux_spherical",
    "solution_to_angular_flux_cylindrical",
    "angular_flux_to_scalar",
    "transport_operator_matvec",
    "transport_operator_matvec_spherical",
    "transport_operator_matvec_cylindrical",
    "build_transport_linear_operator",
    "build_transport_linear_operator_spherical",
    "build_transport_linear_operator_cylindrical",
    "build_rhs",
    "build_rhs_spherical",
    "build_rhs_cylindrical",
    "SNStreamingOperator",
]


# ═══════════════════════════════════════════════════════════════════════
# Equation map: which (ordinate, cell) pairs are unknowns
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class EquationMap:
    """Mapping between 1D solution vector and 4D angular flux."""

    n_eq: int               # number of angular unknowns (per group)
    n_unknowns: int         # n_eq * ng (total scalar unknowns)
    ordinate: np.ndarray    # (n_eq,) ordinate index for each unknown
    ix: np.ndarray          # (n_eq,) x-cell index
    iy: np.ndarray          # (n_eq,) y-cell index


def build_equation_map(
    nx: int, ny: int, quad: AngularQuadrature, ng: int,
) -> EquationMap:
    """Identify which (ordinate, cell) combos are unknowns.

    Filter: mu_z >= 0 (upper hemisphere), and NOT incoming at
    reflective boundaries (those are determined by reflection).
    """
    mu_x, mu_y = quad.mu_x, quad.mu_y
    # Need mu_z — for LebedevSphere it's available, for GL1D it's 0
    mu_z = getattr(quad, 'mu_z', np.zeros(quad.N))

    ords, ixs, iys = [], [], []
    for iy in range(ny):
        for ix in range(nx):
            for n in range(quad.N):
                if mu_z[n] < -1e-15:
                    continue
                if ix == 0 and mu_x[n] > 1e-15:
                    continue
                if ix == nx - 1 and mu_x[n] < -1e-15:
                    continue
                if iy == 0 and mu_y[n] > 1e-15:
                    continue
                if iy == ny - 1 and mu_y[n] < -1e-15:
                    continue
                ords.append(n)
                ixs.append(ix)
                iys.append(iy)

    n_eq = len(ords)
    return EquationMap(
        n_eq=n_eq,
        n_unknowns=n_eq * ng,
        ordinate=np.array(ords, dtype=int),
        ix=np.array(ixs, dtype=int),
        iy=np.array(iys, dtype=int),
    )


# ═══════════════════════════════════════════════════════════════════════
# Solution ↔ angular flux conversion
# ═══════════════════════════════════════════════════════════════════════

def solution_to_angular_flux(
    solution: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    nx: int, ny: int, ng: int,
) -> np.ndarray:
    """Convert 1D solution vector to 4D angular flux (ng, N, nx, ny).

    Applies z-reflection and reflective BCs to fill the full array
    from the reduced set of unknowns.
    """
    mu_x, mu_y = quad.mu_x, quad.mu_y
    mu_z = getattr(quad, 'mu_z', np.zeros(quad.N))
    ref_x = quad.reflection_index("x")
    ref_y = quad.reflection_index("y")
    # z-reflection: need ref_z for Lebedev
    ref_z = getattr(quad, '_ref_z', np.arange(quad.N))

    fi = np.zeros((ng, quad.N, nx, ny))

    # Scatter solution into fi
    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    for k in range(eq_map.n_eq):
        fi[:, eq_map.ordinate[k], eq_map.ix[k], eq_map.iy[k]] = flux[:, k]

    # Z-reflection
    for n in range(quad.N):
        if mu_z[n] < -1e-15:
            fi[:, n, :, :] = fi[:, ref_z[n], :, :]

    # X reflective BCs
    for n in range(quad.N):
        if mu_x[n] > 1e-15:
            fi[:, n, 0, :] = fi[:, ref_x[n], 0, :]
        if mu_x[n] < -1e-15:
            fi[:, n, -1, :] = fi[:, ref_x[n], -1, :]

    # Y reflective BCs
    for n in range(quad.N):
        if mu_y[n] > 1e-15:
            fi[:, n, :, 0] = fi[:, ref_y[n], :, 0]
        if mu_y[n] < -1e-15:
            fi[:, n, :, -1] = fi[:, ref_y[n], :, -1]

    return fi


def angular_flux_to_scalar(
    fi: np.ndarray, quad: AngularQuadrature, nx: int, ny: int, ng: int,
) -> np.ndarray:
    """Integrate angular flux to scalar flux: φ = Σ w_n ψ_n."""
    sf = np.zeros((nx, ny, ng))
    for iy in range(ny):
        for ix in range(nx):
            sf[ix, iy, :] = np.sum(
                fi[:, :, ix, iy] * quad.weights[None, :], axis=1,
            )
    return sf


# ═══════════════════════════════════════════════════════════════════════
# Finite-difference gradients (diamond scheme with reflective BCs)
# ═══════════════════════════════════════════════════════════════════════

def _compute_gradients(
    fi: np.ndarray,
    n: int, ix: int, iy: int,
    quad: AngularQuadrature,
    nx: int, ny: int,
    dx: np.ndarray, dy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Upwind cell-center gradients with reflective BCs.

    Returns (dfi/dx, dfi/dy), each shape (ng,).

    The gradient between adjacent cell centers is divided by the
    **cell-center distance** ``(dx[i] + dx[j]) / 2`` rather than
    the local cell width ``dx[i]``. On a uniform mesh these are
    identical; on a non-uniform mesh the cell-center distance is
    the correct denominator for a first-order consistent FD stencil.
    The original MATLAB code (Mikityuk, PSI 2015) used a scalar
    ``g.delta`` (uniform mesh only), so this distinction did not arise.
    """
    ref_x = quad.reflection_index("x")
    ref_y = quad.reflection_index("y")
    mu_x, mu_y = quad.mu_x, quad.mu_y

    # X gradient
    if mu_x[n] > 1e-15:
        if ix == 0:
            dfix = fi[:, ref_x[n], ix, iy] - fi[:, ref_x[n], ix + 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix + 1])
        else:
            dfix = fi[:, n, ix, iy] - fi[:, n, ix - 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix - 1])
    elif mu_x[n] < -1e-15:
        if ix == nx - 1:
            dfix = fi[:, ref_x[n], ix - 1, iy] - fi[:, ref_x[n], ix, iy]
            hx = 0.5 * (dx[ix - 1] + dx[ix])
        else:
            dfix = fi[:, n, ix + 1, iy] - fi[:, n, ix, iy]
            hx = 0.5 * (dx[ix + 1] + dx[ix])
    else:
        dfix = np.zeros(fi.shape[0])
        hx = 1.0

    # Y gradient
    if mu_y[n] > 1e-15:
        if iy == 0:
            dfiy = fi[:, ref_y[n], ix, iy] - fi[:, ref_y[n], ix, iy + 1]
            hy = 0.5 * (dy[iy] + dy[iy + 1])
        else:
            dfiy = fi[:, n, ix, iy] - fi[:, n, ix, iy - 1]
            hy = 0.5 * (dy[iy] + dy[iy - 1])
    elif mu_y[n] < -1e-15:
        if iy == ny - 1:
            dfiy = fi[:, ref_y[n], ix, iy - 1] - fi[:, ref_y[n], ix, iy]
            hy = 0.5 * (dy[iy - 1] + dy[iy])
        else:
            dfiy = fi[:, n, ix, iy + 1] - fi[:, n, ix, iy]
            hy = 0.5 * (dy[iy + 1] + dy[iy])
    else:
        dfiy = np.zeros(fi.shape[0])
        hy = 1.0

    return dfix / hx, dfiy / hy


# ═══════════════════════════════════════════════════════════════════════
# Transport operator  T: ψ → μ·∇ψ + Σ_t·ψ
# ═══════════════════════════════════════════════════════════════════════

def transport_operator_matvec(
    solution: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ny: int, ng: int,
    dx: np.ndarray, dy: np.ndarray,
) -> np.ndarray:
    """Apply the streaming + collision operator T·ψ.

    Parameters
    ----------
    solution : (n_unknowns,) flattened angular flux vector.
    sig_t : (nx, ny, ng) total cross section.

    Returns
    -------
    np.ndarray
        Shape ``(n_unknowns,)``. T applied to the angular flux.
    """
    fi = solution_to_angular_flux(solution, eq_map, quad, nx, ny, ng)

    lhs = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        n, ix, iy = eq_map.ordinate[k], eq_map.ix[k], eq_map.iy[k]
        dfidx, dfidy = _compute_gradients(fi, n, ix, iy, quad, nx, ny, dx, dy)
        lhs[:, k] = (
            quad.mu_x[n] * dfidx
            + quad.mu_y[n] * dfidy
            + sig_t[ix, iy, :] * fi[:, n, ix, iy]
        )

    return lhs.ravel(order='F')


def build_transport_linear_operator(
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ny: int, ng: int,
    dx: np.ndarray, dy: np.ndarray,
) -> LinearOperator:
    """Build scipy LinearOperator for T = μ·∇ + Σ_t."""
    def matvec(x):
        return transport_operator_matvec(
            x, eq_map, quad, sig_t, nx, ny, ng, dx, dy,
        )

    n = eq_map.n_unknowns
    return LinearOperator((n, n), matvec=matvec, dtype=float)


# ═══════════════════════════════════════════════════════════════════════
# RHS construction (fission + scattering + n2n, per ordinate)
# ═══════════════════════════════════════════════════════════════════════

def build_rhs(
    fission_source: np.ndarray,
    scalar_flux: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_s: dict[int, list[np.ndarray]],
    sig2: dict[int, np.ndarray],
    mat_map: np.ndarray,
    nx: int, ny: int, ng: int,
    scattering_order: int = 0,
    angular_flux: np.ndarray | None = None,
) -> np.ndarray:
    """Build the RHS source vector for T·ψ = b.

    All isotropic sources are divided by sum(weights) — the angular
    normalization for the discrete angular flux equation.

    For Pn scattering (scattering_order > 0), the scattering source
    is per-ordinate using Legendre moments of the angular flux::

        qS(n) = Σ_l (2l+1) · Σ_s^l^T @ [Σ_m fiL_lm · Y_lm(n)] / sum_w

    Parameters
    ----------
    fission_source : (nx, ny, ng) — already divided by sum(w) by the caller.
    scalar_flux : (nx, ny, ng) — current scalar flux for scattering.
    sig_s : dict[mat_id → list of (ng, ng)] Legendre scattering matrices.
    sig2 : dict[mat_id → (ng, ng)] (n,2n) matrices.
    scattering_order : Legendre order L (0 = P0 isotropic).
    angular_flux : (ng, N, nx, ny) angular flux for computing Legendre
        moments. Required if scattering_order > 0.

    Returns
    -------
    np.ndarray
        Shape ``(n_unknowns,)``. The RHS source vector.
    """
    sum_w = float(quad.weights.sum())
    L = scattering_order
    mu_z = getattr(quad, 'mu_z', np.zeros(quad.N))

    # Precompute Legendre moments if anisotropic scattering
    fiL = None
    Y = None
    if L > 0 and angular_flux is not None:
        Y = quad.spherical_harmonics(L)  # (N, L+1, 2L+1)
        w = quad.weights
        fiL = np.zeros((nx, ny, ng, L + 1, 2 * L + 1))
        for l in range(L + 1):
            for m in range(-l, l + 1):
                for n in range(quad.N):
                    fiL[:, :, :, l, l + m] += (
                        w[n] * angular_flux[:, n, :, :].T * Y[n, l, l + m]
                    )

    rhs = np.zeros((ng, eq_map.n_eq))
    eq_idx = 0
    for iy in range(ny):
        for ix in range(nx):
            mid = int(mat_map[ix, iy])
            phi_cell = scalar_flux[ix, iy, :]

            # Fission (already normalized by caller)
            qF = fission_source[ix, iy, :]

            # (n,2n) — isotropic
            q2 = 2.0 * (sig2[mid].T @ phi_cell) / sum_w

            for n in range(quad.N):
                if mu_z[n] < -1e-15:
                    continue
                if ix == 0 and quad.mu_x[n] > 1e-15:
                    continue
                if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                    continue
                if iy == 0 and quad.mu_y[n] > 1e-15:
                    continue
                if iy == ny - 1 and quad.mu_y[n] < -1e-15:
                    continue

                # Scattering source (Pn expansion)
                qS = np.zeros(ng)
                for l in range(L + 1):
                    if l == 0:
                        # P0: isotropic, use scalar flux
                        qS += sig_s[mid][0].T @ phi_cell / sum_w
                    elif fiL is not None:
                        # P1+: anisotropic, use Legendre moments
                        SUM = np.zeros(ng)
                        for m in range(-l, l + 1):
                            SUM += fiL[ix, iy, :, l, l + m] * Y[n, l, l + m]
                        qS += (2 * l + 1) * (sig_s[mid][l].T @ SUM) / sum_w

                rhs[:, eq_idx] = qF + q2 + qS
                eq_idx += 1

    return rhs.ravel(order='F')


# ═══════════════════════════════════════════════════════════════════════
# Spherical 1D operator: T = μ(A∂/∂r)/V + (α∂/∂μ)/V + Σ_t
# ═══════════════════════════════════════════════════════════════════════

def build_equation_map_spherical(
    nx: int, quad: AngularQuadrature, ng: int,
) -> EquationMap:
    """Equation map for spherical 1D: all (ordinate, cell) pairs except
    incoming directions at the outer reflective boundary."""
    mu_x = quad.mu_x
    ords, ixs, iys = [], [], []
    for ix in range(nx):
        for n in range(quad.N):
            # Skip incoming at outer boundary (reflective BC)
            if ix == nx - 1 and mu_x[n] < -1e-15:
                continue
            ords.append(n)
            ixs.append(ix)
            iys.append(0)  # always iy=0 for 1D

    n_eq = len(ords)
    return EquationMap(
        n_eq=n_eq,
        n_unknowns=n_eq * ng,
        ordinate=np.array(ords, dtype=int),
        ix=np.array(ixs, dtype=int),
        iy=np.array(iys, dtype=int),
    )


def solution_to_angular_flux_spherical(
    solution: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    nx: int, ng: int,
) -> np.ndarray:
    """Convert 1D solution vector to angular flux array (ng, N, nx, 1).

    Applies reflective BC at the outer boundary.
    """
    ref_x = quad.reflection_index("x")
    fi = np.zeros((ng, quad.N, nx, 1))

    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    for k in range(eq_map.n_eq):
        fi[:, eq_map.ordinate[k], eq_map.ix[k], 0] = flux[:, k]

    # Reflective BC at outer boundary: incoming (μ<0) = reflected partner
    for n in range(quad.N):
        if quad.mu_x[n] < -1e-15:
            fi[:, n, -1, 0] = fi[:, ref_x[n], -1, 0]

    return fi


def transport_operator_matvec_spherical(
    solution: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ng: int,
    face_areas: np.ndarray,
    volumes: np.ndarray,
    alpha_half: np.ndarray,
    redist_dAw: np.ndarray,
    tau_mm: np.ndarray,
) -> np.ndarray:
    r"""Apply the spherical transport operator T·ψ.

    .. math::

        (T\psi)_{n,i} = \frac{\mu_n}{V_i}
          \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
        + \frac{\Delta A_i}{w_n V_i}
          \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
        + \Sigma_t \psi_{n,i}

    The :math:`\Delta A / w` geometry factor (``redist_dAw``, precomputed
    in :class:`SNMesh`) ensures per-ordinate flat-flux consistency
    (Bailey et al. 2009).

    Face fluxes are approximated by arithmetic averages of cell-centre values.
    """
    fi = solution_to_angular_flux_spherical(solution, eq_map, quad, nx, ng)
    ref_x = quad.reflection_index("x")
    A = face_areas       # (nx+1,)
    V = volumes[:, 0]    # (nx,)
    dAw = redist_dAw     # (nx, N) precomputed ΔA_i/w_n
    alpha = alpha_half   # (N+1,) non-negative dome
    N = quad.N
    mu = quad.mu_x

    lhs = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        psi_ni = fi[:, n, i, 0]

        # ── Spatial streaming: μ (A ∂ψ/∂r) / V ──────────────────────
        if i < nx - 1:
            psi_right = 0.5 * (fi[:, n, i, 0] + fi[:, n, i + 1, 0])
        else:
            if mu[n] > 1e-15:
                psi_right = fi[:, n, i, 0]
            else:
                psi_right = fi[:, ref_x[n], i, 0]

        if i > 0:
            psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
        else:
            psi_left = 0.0

        streaming = mu[n] * (A[i + 1] * psi_right - A[i] * psi_left) / V[i]

        # ── Angular redistribution: (ΔA/w) (α ∂ψ/∂μ) / V ──────────
        # Angular face flux uses M-M weighted interpolation (τ).
        dA_w = dAw[i, n]  # precomputed geometry factor
        tau_n = tau_mm[n]

        if n < N - 1:
            psi_angle_right = tau_n * fi[:, n + 1, i, 0] + (1.0 - tau_n) * fi[:, n, i, 0]
        else:
            psi_angle_right = fi[:, n, i, 0]

        if n > 0:
            psi_angle_left = tau_mm[n - 1] * fi[:, n, i, 0] + (1.0 - tau_mm[n - 1]) * fi[:, n - 1, i, 0]
        else:
            psi_angle_left = fi[:, n, i, 0]

        redistribution = dA_w * (alpha[n + 1] * psi_angle_right
                                 - alpha[n] * psi_angle_left) / V[i]

        # ── Collision ────────────────────────────────────────────────
        collision = sig_t[i, 0, :] * psi_ni

        lhs[:, k] = streaming + redistribution + collision

    return lhs.ravel(order='F')


def build_transport_linear_operator_spherical(
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ng: int,
    face_areas: np.ndarray,
    volumes: np.ndarray,
    alpha_half: np.ndarray,
    redist_dAw: np.ndarray,
    tau_mm: np.ndarray,
) -> LinearOperator:
    """Build scipy LinearOperator for spherical T."""
    def matvec(x):
        return transport_operator_matvec_spherical(
            x, eq_map, quad, sig_t, nx, ng,
            face_areas, volumes, alpha_half, redist_dAw, tau_mm,
        )

    n = eq_map.n_unknowns
    return LinearOperator((n, n), matvec=matvec, dtype=float)


def build_rhs_spherical(
    fission_source: np.ndarray,
    scalar_flux: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_s: dict[int, list[np.ndarray]],
    sig2: dict[int, np.ndarray],
    mat_map: np.ndarray,
    nx: int, ng: int,
    scattering_order: int = 0,
    angular_flux: np.ndarray | None = None,
) -> np.ndarray:
    """Build the RHS source vector for spherical T·ψ = b.

    Same structure as Cartesian ``build_rhs`` but with spherical
    equation map (no y-direction, no z-reflection filtering).
    """
    sum_w = float(quad.weights.sum())
    L = scattering_order

    rhs = np.zeros((ng, eq_map.n_eq))
    eq_idx = 0
    for ix in range(nx):
        mid = int(mat_map[ix, 0])
        phi_cell = scalar_flux[ix, 0, :]

        qF = fission_source[ix, 0, :]
        q2 = 2.0 * (sig2[mid].T @ phi_cell) / sum_w

        for n in range(quad.N):
            if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                continue

            qS = sig_s[mid][0].T @ phi_cell / sum_w
            rhs[:, eq_idx] = qF + q2 + qS
            eq_idx += 1

    return rhs.ravel(order='F')


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical 1D operator: T = η(A∂/∂r)/V + (ΔA/w)(α∂/∂φ)/V + Σ_t
# ═══════════════════════════════════════════════════════════════════════

# Equation map and solution-to-flux reuse the spherical versions
# (both are 1D with reflective BC at the outer boundary).
build_equation_map_cylindrical = build_equation_map_spherical
solution_to_angular_flux_cylindrical = solution_to_angular_flux_spherical


def transport_operator_matvec_cylindrical(
    solution: np.ndarray,
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ng: int,
    face_areas: np.ndarray,
    volumes: np.ndarray,
    alpha_per_level: list[np.ndarray],
    redist_dAw_per_level: list[np.ndarray],
    tau_mm_per_level: list[np.ndarray],
) -> np.ndarray:
    r"""Apply the cylindrical transport operator T·ψ.

    Per-level azimuthal redistribution with geometry-weighted
    :math:`\Delta A / w` factor and Morel–Montry angular closure.
    """
    fi = solution_to_angular_flux_cylindrical(solution, eq_map, quad, nx, ng)
    ref_x = quad.reflection_index("x")
    A = face_areas       # (nx+1,)
    V = volumes[:, 0]    # (nx,)
    N = quad.N
    mu = quad.mu_x

    # Build reverse map: global ordinate → (level, local index)
    ord_to_level = np.empty(N, dtype=int)
    ord_to_local = np.empty(N, dtype=int)
    for p, level_idx in enumerate(quad.level_indices):
        for m_local, n in enumerate(level_idx):
            ord_to_level[n] = p
            ord_to_local[n] = m_local

    lhs = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        psi_ni = fi[:, n, i, 0]

        p = ord_to_level[n]
        m_local = ord_to_local[n]
        alpha = alpha_per_level[p]
        dAw = redist_dAw_per_level[p]
        tau_level = tau_mm_per_level[p]
        level_idx = quad.level_indices[p]
        M = len(level_idx)

        # ── Spatial streaming: η (A ∂ψ/∂r) / V ─────────────────────
        if i < nx - 1:
            psi_right = 0.5 * (fi[:, n, i, 0] + fi[:, n, i + 1, 0])
        else:
            if mu[n] > 1e-15:
                psi_right = fi[:, n, i, 0]
            else:
                psi_right = fi[:, ref_x[n], i, 0]

        if i > 0:
            psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
        else:
            psi_left = 0.0

        streaming = mu[n] * (A[i + 1] * psi_right - A[i] * psi_left) / V[i]

        # ── Angular redistribution: (ΔA/w)(α ∂ψ/∂φ) / V ───────────
        dA_w = dAw[i, m_local]
        tau_m = tau_level[m_local]

        if m_local < M - 1:
            n_next = level_idx[m_local + 1]
            psi_angle_right = tau_m * fi[:, n_next, i, 0] + (1.0 - tau_m) * fi[:, n, i, 0]
        else:
            psi_angle_right = fi[:, n, i, 0]

        if m_local > 0:
            n_prev = level_idx[m_local - 1]
            tau_prev = tau_level[m_local - 1]
            psi_angle_left = tau_prev * fi[:, n, i, 0] + (1.0 - tau_prev) * fi[:, n_prev, i, 0]
        else:
            psi_angle_left = fi[:, n, i, 0]

        redistribution = dA_w * (alpha[m_local + 1] * psi_angle_right
                                 - alpha[m_local] * psi_angle_left) / V[i]

        # ── Collision ────────────────────────────────────────────────
        collision = sig_t[i, 0, :] * psi_ni

        lhs[:, k] = streaming + redistribution + collision

    return lhs.ravel(order='F')


def build_transport_linear_operator_cylindrical(
    eq_map: EquationMap,
    quad: AngularQuadrature,
    sig_t: np.ndarray,
    nx: int, ng: int,
    face_areas: np.ndarray,
    volumes: np.ndarray,
    alpha_per_level: list[np.ndarray],
    redist_dAw_per_level: list[np.ndarray],
    tau_mm_per_level: list[np.ndarray],
) -> LinearOperator:
    """Build scipy LinearOperator for cylindrical T."""
    def matvec(x):
        return transport_operator_matvec_cylindrical(
            x, eq_map, quad, sig_t, nx, ng,
            face_areas, volumes,
            alpha_per_level, redist_dAw_per_level, tau_mm_per_level,
        )

    n = eq_map.n_unknowns
    return LinearOperator((n, n), matvec=matvec, dtype=float)


# RHS builder reuses the spherical version (same 1D isotropic structure).
build_rhs_cylindrical = build_rhs_spherical


# ═══════════════════════════════════════════════════════════════════════
# SNStreamingOperator — unified LinearOperator for L = Ω·∇ + Σ_t
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class SNStreamingOperator(LinearOperatorMixin):
    r"""Streaming-collision operator :math:`L = \Omega\cdot\nabla + \Sigma_t`
    as a :class:`~orpheus.numerics.operator.LinearOperator`.

    Wave D Round 3 capstone (Issue #160).  Implements the Wave A
    :class:`~orpheus.numerics.operator.LinearOperator` Protocol with
    three capabilities:

    * **apply** — matrix-free forward action :math:`L\,\psi`.
      Reuses the symmetric closure math from
      :func:`transport_operator_matvec` (Cartesian upwind FD) /
      :func:`transport_operator_matvec_spherical` (arithmetic face
      averages + τ-weighted symmetric M-M angular interpolation) /
      :func:`transport_operator_matvec_cylindrical` (per-level
      azimuthal redistribution with M-M closure).  The math is
      **extracted verbatim** from those functions; this class is a
      thin LinearOperator-Protocol wrapper.

    * **solve** — :math:`L^{-1}\,q` via the Wave D Round 2 unified
      sweep (:func:`orpheus.sn.sweep.transport_sweep`).  This uses
      the WDD asymmetric closure of
      :class:`~orpheus.sn.spatial.diamond.DiamondDifference` (the
      sweep's existing math, ERR-026 affected for curvilinear).

    * **apply_transpose** — adjoint action :math:`L^*\,\psi`.
      Constructed from the explicit transpose of the dense matrix
      assembled by probing :meth:`apply` with each unit basis vector
      (one-time cost cached in :attr:`_dense_matrix`).  This
      construction is exact by linear algebra: every linear operator
      has a transpose, and the transpose of the assembled matrix
      *is* the operator's transpose.  Reciprocity
      :math:`\langle L\psi, \varphi\rangle = \langle \psi, L^*\varphi\rangle`
      holds to round-off by construction (see
      ``tests/sn/test_snstreamingoperator.py`` for the gating tests).

    Why ``apply`` and ``solve`` use **different** closures (by design)
    -----------------------------------------------------------------

    The historical sweep (:func:`transport_sweep`) and the historical
    finite-difference operator (the ``transport_operator_matvec_*``
    functions in this module) were built at different times for
    different consumers (source iteration vs BiCGSTAB) and ship
    **two distinct discretisations** of the same continuous operator.

    * **apply** carries the **symmetric** discretisation: upwind
      cell-center FD on Cartesian; arithmetic face averages with
      τ-symmetric Morel-Montry angular interpolation on curvilinear.
      This is the closure that makes the BiCGSTAB Krylov path agree
      with analytical references (per ERR-026 evidence table — see
      :ref:`sn-streaming-operator`).

    * **solve** carries the **WDD asymmetric** closure
      :math:`\psi_{n+1/2} = (\overline\psi - (1-\tau)\,\psi_{n-1/2})
      /\tau`.  This is the historical SN sweep's closure (the
      forward-substitution-friendly upper-triangular form that lets
      a sweep run in :math:`O(N\cdot N_{\rm cells})` work).

    On uniform meshes both closures converge to the same physics in
    the fine-mesh limit; on coarse meshes they differ at
    :math:`O(h)` (Cartesian) or have a closure-bias-driven
    self-consistent fixed point on curvilinear (ERR-026, deferred to
    Wave E).  The Wave E Issue 15 reconciliation is **Krylov-on-apply
    with solve as preconditioner** — the Krylov outer iteration uses
    the symmetric closure (correct discretisation) while the sweep
    is invoked as a preconditioner only, so its closure bias does
    not poison the converged solution.

    Capability set
    --------------

    ``frozenset({"apply", "solve", "apply_transpose"})`` — the operator
    is a full citizen of the Wave A operator algebra: it composes with
    :class:`~orpheus.sn.scattering.ScatteringOperator` and
    :class:`~orpheus.sn.fission.FissionOperator` under
    :math:`(L - S - F)`; the composition's capability set falls out
    of the Wave A closure laws (see :ref:`operator-algebra`).

    Vector layout
    -------------

    :meth:`apply` and :meth:`apply_transpose` operate on the **packed
    1-D solution vector** (shape ``(n_unknowns,)``) used by the
    BiCGSTAB FD operator path: an :class:`EquationMap` selects
    which ``(ordinate, cell)`` combinations are unknowns (the rest
    are determined by reflective BCs and z-hemisphere reflection),
    and the vector is laid out group-major in Fortran order.  This
    is the natural input shape for
    :func:`scipy.sparse.linalg.bicgstab`.

    :meth:`solve` operates on **structured arrays**: source ``Q``
    shape ``(nx, ny, ng)`` plus persistent boundary-flux dict
    ``psi_bc`` and optional anisotropic source ``Q_aniso`` shape
    ``(N, nx, ny, ng)``.  It returns a ``(angular_flux, scalar_flux)``
    tuple matching :func:`transport_sweep`'s contract.  The shape
    mismatch between the packed-vector ``apply`` and the
    structured-array ``solve`` reflects the historical layouts of
    the two consumers; Wave E will normalise these via a single
    Krylov-on-apply path.

    Attributes
    ----------
    sn_mesh : SNMesh
        The augmented geometry carrying quadrature, materials, BCs,
        cell-update strategy, and (for curvilinear) the precomputed
        connection coefficients ``alpha_half``, ``redist_dAw``,
        ``tau_mm`` (or per-level analogues for cylindrical).
    sig_t : np.ndarray
        Total cross-section, shape ``(nx, ny, ng)``.  Held as a
        separate attribute (not derived from ``sn_mesh``) because
        the existing solver passes it around explicitly.
    capabilities : frozenset[str]
        ``{"apply", "solve", "apply_transpose"}``.

    Notes
    -----
    The bit-identical regression contract for the SN reshape
    campaign holds because :class:`SNStreamingOperator` is an
    **additive** code path: existing solver paths in
    :mod:`orpheus.sn.solver` continue to use the legacy
    ``transport_operator_matvec_*`` and ``transport_sweep`` APIs
    directly; nothing in the existing solver paths changes when
    this class is added.  Wave E Issue 15 wires the solver to the
    operator algebra; that is where the campaign's closure
    reconciliation actually happens.
    """

    sn_mesh: "SNMesh"
    sig_t: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        )
    )

    # Lazy caches (constructed on first call).
    _eq_map: EquationMap | None = field(default=None, init=False, repr=False)
    _dense_matrix: np.ndarray | None = field(
        default=None, init=False, repr=False,
    )
    _dense_matrix_T: np.ndarray | None = field(
        default=None, init=False, repr=False,
    )

    # ── EquationMap dispatch ──────────────────────────────────────────

    def _ensure_eq_map(self) -> EquationMap:
        """Lazily build the geometry-appropriate :class:`EquationMap`.

        The same dispatch logic the legacy BiCGSTAB paths use
        (:meth:`SNSolver._solve_bicgstab_*`).  Built once and cached
        on first :meth:`apply` / :meth:`apply_transpose` call.
        """
        if self._eq_map is None:
            nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.sig_t.shape[2]
            quad = self.sn_mesh.quad
            curv = getattr(self.sn_mesh, "curvature", None)
            if curv == "spherical":
                self._eq_map = build_equation_map_spherical(nx, quad, ng)
            elif curv == "cylindrical":
                self._eq_map = build_equation_map_cylindrical(nx, quad, ng)
            else:
                self._eq_map = build_equation_map(nx, ny, quad, ng)
        return self._eq_map

    @property
    def n_unknowns(self) -> int:
        """Total scalar unknowns ``n_eq * ng`` for the packed vector."""
        return self._ensure_eq_map().n_unknowns

    # ── apply: forward action L·ψ ─────────────────────────────────────

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Forward action :math:`L\,\psi` on the packed solution vector.

        Dispatches by ``self.sn_mesh.curvature`` to the existing
        finite-difference matvec routines:

        * Cartesian (curvature ``None``):
          :func:`transport_operator_matvec`.
        * Spherical (curvature ``"spherical"``):
          :func:`transport_operator_matvec_spherical`.
        * Cylindrical (curvature ``"cylindrical"``):
          :func:`transport_operator_matvec_cylindrical`.

        The math is **extracted verbatim** from those functions; this
        method is a thin protocol-conforming wrapper.

        Parameters
        ----------
        psi : np.ndarray
            Packed solution vector, shape ``(n_unknowns,)`` where
            ``n_unknowns = eq_map.n_eq * ng`` is determined by
            :meth:`_ensure_eq_map` for this geometry.

        Returns
        -------
        np.ndarray
            ``L\,\psi`` as a packed vector, same shape as ``psi``.
        """
        eq_map = self._ensure_eq_map()
        sn_mesh = self.sn_mesh
        nx, ny, ng = sn_mesh.nx, sn_mesh.ny, self.sig_t.shape[2]
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)

        if curv == "spherical":
            # Use ``self.reduced.{face_areas, alpha_half, redist_dAw,
            # tau_mm}`` directly (Wave D R1 canonical accessor) to
            # avoid the deprecated property's DeprecationWarning.
            reduced = sn_mesh.reduced
            return transport_operator_matvec_spherical(
                psi, eq_map, quad, self.sig_t,
                nx, ng,
                reduced.face_areas,
                sn_mesh.volumes,
                reduced.alpha_half,
                reduced.redist_dAw,
                reduced.tau_mm,
            )
        if curv == "cylindrical":
            reduced = sn_mesh.reduced
            return transport_operator_matvec_cylindrical(
                psi, eq_map, quad, self.sig_t,
                nx, ng,
                reduced.face_areas,
                sn_mesh.volumes,
                reduced.alpha_per_level,
                reduced.redist_dAw_per_level,
                reduced.tau_mm_per_level,
            )
        # Cartesian (1-D slab or 2-D)
        return transport_operator_matvec(
            psi, eq_map, quad, self.sig_t,
            nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
        )

    # ── solve: L^{-1}·q via the unified sweep ─────────────────────────

    def solve(
        self,
        Q: np.ndarray,
        psi_bc: dict | None = None,
        Q_aniso: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Inverse action :math:`L^{-1}\,q` via the Wave D Round 2 sweep.

        Delegates to :func:`orpheus.sn.sweep.transport_sweep` with
        ``self.sn_mesh.cell_update`` (defaulting to
        :class:`DiamondDifference`).  Bit-identical to a direct
        :func:`transport_sweep` call on the same arguments.

        **The closure used by :meth:`solve` is asymmetric (WDD)**, which
        differs from the symmetric closure :meth:`apply` carries.  See
        the class docstring "Why ``apply`` and ``solve`` use different
        closures (by design)" for the rationale.

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source density, shape ``(nx, ny, ng)``.
        psi_bc : dict or None
            Persistent boundary-flux dict storing ψ on each face
            between outer iterations.  If ``None``, a fresh empty
            dict is supplied; the caller cannot then carry state
            between sweeps.
        Q_aniso : np.ndarray or None
            Per-ordinate anisotropic source, shape ``(N, nx, ny, ng)``,
            for P1+ scattering.  ``None`` for isotropic-only (P0).

        Returns
        -------
        tuple
            ``(angular_flux, scalar_flux)`` matching the
            :func:`transport_sweep` contract:

            * ``angular_flux`` shape ``(N, nx, ny, ng)``,
            * ``scalar_flux`` shape ``(nx, ny, ng)``.
        """
        from .sweep import transport_sweep
        if psi_bc is None:
            psi_bc = {}
        return transport_sweep(Q, self.sig_t, self.sn_mesh, psi_bc, Q_aniso)

    # ── apply_transpose: adjoint action L*·ψ via dense transpose ──────

    def _ensure_dense_matrix(self) -> np.ndarray:
        """Assemble the dense matrix of :meth:`apply` by probing.

        Calls :meth:`apply` with each of the ``n_unknowns`` unit basis
        vectors, stacks the outputs as columns, and returns the
        resulting ``(n_unknowns, n_unknowns)`` matrix.  Cost is
        :math:`O(n_{\rm unknowns}^2)` time and space.

        The dense assembly is a one-time cost on the first
        :meth:`apply_transpose` call; the matrix is cached on
        ``self._dense_matrix`` and its transpose on
        ``self._dense_matrix_T``.

        For the small reciprocity test problems (slab GL N=4,
        nx=4 / sphere GL N=4, nx=4 — ``n_unknowns ~ 30-150``) the
        cost is negligible.  For production-scale problems
        (``n_unknowns ~ 10^4-10^6``) the dense path is unsuitable;
        Wave E will ship an :math:`O(n)` analytic-adjoint matvec.
        """
        if self._dense_matrix is None:
            n = self.n_unknowns
            mat = np.empty((n, n), dtype=float)
            basis = np.zeros(n, dtype=float)
            for j in range(n):
                basis[j] = 1.0
                mat[:, j] = self.apply(basis)
                basis[j] = 0.0
            self._dense_matrix = mat
            self._dense_matrix_T = mat.T.copy()
        return self._dense_matrix

    def apply_transpose(self, psi: np.ndarray) -> np.ndarray:
        r"""Adjoint action :math:`L^*\,\psi` on the packed vector.

        Implemented via the explicit transpose of the dense matrix
        assembled by probing :meth:`apply`.  Reciprocity
        :math:`\langle L\,\psi,\,\varphi\rangle = \langle\psi,\,
        L^*\,\varphi\rangle` holds to round-off by construction:
        every linear operator has a transpose, and the transpose of
        the assembled dense matrix *is* the operator's transpose.

        The reciprocity test in
        ``tests/sn/test_snstreamingoperator.py`` verifies this
        identity holds across slab, spherical, and cylindrical
        geometries for synthetic ``(ψ, φ)`` pairs.  A failure of
        that test would indicate :meth:`apply` is non-linear or
        the dense-assembly probe code is wrong — both are
        catastrophic operator-correctness failures.

        Parameters
        ----------
        psi : np.ndarray
            Packed solution vector, shape ``(n_unknowns,)``.

        Returns
        -------
        np.ndarray
            :math:`L^*\,\psi` as a packed vector.
        """
        self._ensure_dense_matrix()
        assert self._dense_matrix_T is not None
        return self._dense_matrix_T @ psi

