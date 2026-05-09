"""Direct transport operator for Krylov inner solves.

Provides the explicit symmetric-closure transport operator
:math:`L = \\Omega\\cdot\\nabla + \\Sigma_t` as
:class:`SNStreamingOperator`, plus the packed-vector layout helpers
(:class:`EquationMap`, :func:`solution_to_angular_flux*`,
:func:`transport_operator_matvec*`) that the operator's
:meth:`apply` path uses internally.

Three geometries are supported:

* **Cartesian 2D** — ``L = μ_x ∂/∂x + μ_y ∂/∂y + Σ_t``
* **Spherical 1D** — ``L = μ (A ∂/∂r)/V + (α ∂/∂μ)/V + Σ_t``
* **Cylindrical 1D** — per-level azimuthal redistribution

The sweep-based solver (source iteration) inverts :math:`L` implicitly
via diamond-difference sweeps.  This module forms :math:`L` explicitly
so that scipy's Krylov solvers (GMRES) can solve :math:`L\\,\\psi = b`
directly via :class:`SNStreamingOperator.apply`.  Wave E Round 2
retired the standalone ``build_rhs*`` and
``build_transport_linear_operator*`` helpers along with the
``angular_flux_to_scalar`` aggregator: the Krylov consumer in
:mod:`orpheus.sn.solver` now wraps :meth:`SNStreamingOperator.apply`
as a ``scipy.sparse.linalg.LinearOperator`` directly, with the sweep
as preconditioner.

.. note:: Symmetric-closure invariant

   The operator :math:`L` formed here uses the **symmetric** closure
   that makes the Krylov path agree with analytical references:

   * **Cartesian**: upwind cell-center finite differences for the
     streaming gradient — first-order accurate and consistent with DD
     on **uniform** meshes, but first-order inconsistent on
     non-uniform meshes (divides by the local ``dx[ix]`` rather than
     the cell-center distance ``(dx[ix]+dx[ix±1])/2``).
   * **Curvilinear**: arithmetic averages for spatial face fluxes and
     τ-weighted interpolation for angular face fluxes — the
     symmetric form, distinct from the WDD asymmetric closure
     :math:`\\psi_{\\rm out} = (\\overline{\\psi} - (1 - \\tau)\\,
     \\psi_{\\rm in})/\\tau` used by the sweeps.

   On uniform meshes the symmetric-closure operator and the WDD
   sweep converge to the same physics in the fine-mesh limit; on
   curvilinear, the sweep's WDD closure has the ERR-026 closure-bias-
   driven self-consistent fixed point that is now bypassed by the
   Krylov-on-:meth:`apply` path (Wave E Round 2).

.. note:: Boundary-condition handling — Wave E Round 3

   Wave E Round 3 extended :func:`solution_to_angular_flux*` to consume
   the :class:`~orpheus.geometry.boundary.BoundaryOperator` instances on the
   :class:`~orpheus.sn.geometry.SNMesh` (``bc_xmin``, ``bc_xmax``,
   ``bc_ymin``, ``bc_ymax``).  Each function fills incoming-at-boundary
   slots via ``bc.apply_to_incoming(outgoing, quad)`` (Wave B Issue 7
   tensor-decomposed BC algebra).  Bit-identity to the pre-Round 3
   reflective-only fill is preserved for :class:`SpecularBoundaryOperator` (the
   default ``BC.reflective`` factory), since
   ``SpecularBoundaryOperator(axis="x").apply_to_incoming(out, quad) == out[ref_x]``;
   :class:`VacuumBoundaryOperator` returns zero, which is the correct vacuum
   incoming flux.  This closes ERR-026 for the curvilinear
   ``solve_sn_fixed_source`` MMS path: the FD operator is now
   BC-faithful for vacuum / reflective / white / albedo / mixed BCs
   uniformly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperatorMixin,
)

from .quadrature import AngularQuadrature

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryOperator

    from .geometry import SNMesh

__all__ = [
    "EquationMap",
    "build_equation_map",
    "build_equation_map_spherical",
    "build_equation_map_cylindrical",
    "solution_to_angular_flux",
    "solution_to_angular_flux_spherical",
    "solution_to_angular_flux_cylindrical",
    "transport_operator_matvec",
    "transport_operator_matvec_spherical",
    "transport_operator_matvec_cylindrical",
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
    *,
    bc_xmin: "BoundaryOperator | None" = None,
    bc_xmax: "BoundaryOperator | None" = None,
    bc_ymin: "BoundaryOperator | None" = None,
    bc_ymax: "BoundaryOperator | None" = None,
) -> np.ndarray:
    """Convert 1D solution vector to 4D angular flux (ng, N, nx, ny).

    Applies z-hemisphere reflection (Lebedev) and BC-resolved fills at
    the four boundaries of the 2-D Cartesian domain.

    Wave E Round 3 (ERR-026 closure): the four ``bc_*`` keyword
    arguments accept :class:`~orpheus.geometry.boundary.BoundaryOperator`
    instances built by :class:`~orpheus.sn.geometry.SNMesh` from the
    mesh's :class:`~orpheus.geometry.mesh.BC` declarations.  Each BC's
    ``apply_to_incoming(outgoing, quad)`` method maps the boundary's
    outgoing angular flux to the incoming flux per the BC's tensor
    decomposition (vacuum → 0; specular → ``out[ref]``; white,
    albedo, mixed → their respective combinations).  When all four
    are ``None`` (legacy callers) the behaviour falls back to specular
    reflection on every face — bit-identical to the pre-Round 3
    hard-coded reflective fill that the BiCGSTAB FD path relied on.
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    if bc_xmin is None:
        bc_xmin = SpecularBoundaryOperator(axis="x", albedo=1.0)
    if bc_xmax is None:
        bc_xmax = SpecularBoundaryOperator(axis="x", albedo=1.0)
    if bc_ymin is None:
        bc_ymin = SpecularBoundaryOperator(axis="y", albedo=1.0)
    if bc_ymax is None:
        bc_ymax = SpecularBoundaryOperator(axis="y", albedo=1.0)

    mu_x, mu_y = quad.mu_x, quad.mu_y
    mu_z = getattr(quad, 'mu_z', np.zeros(quad.N))
    # z-reflection: need ref_z for Lebedev
    ref_z = getattr(quad, '_ref_z', np.arange(quad.N))

    fi = np.zeros((ng, quad.N, nx, ny))

    # Scatter solution into fi
    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    for k in range(eq_map.n_eq):
        fi[:, eq_map.ordinate[k], eq_map.ix[k], eq_map.iy[k]] = flux[:, k]

    # Z-reflection (Lebedev / 3-D quadratures): the eq_map only carries
    # mu_z >= 0; the lower hemisphere is filled by reflection across z.
    for n in range(quad.N):
        if mu_z[n] < -1e-15:
            fi[:, n, :, :] = fi[:, ref_z[n], :, :]

    # ── X-axis BC fills ───────────────────────────────────────────────
    # The eq_map skips incoming-at-boundary slots; we fill them here per
    # the BC.  Layout: ``apply_to_incoming`` consumes a full
    # ``(N, ng_axis_2)`` array (here ``(N, ng)``) and returns the same
    # shape; we slice the incoming entries out for the boundary fill.
    # For each y-row independently to avoid spurious coupling.
    for iy in range(ny):
        # Left face (x = xmin): outgoing = mu_x < 0, incoming = mu_x > 0.
        outgoing_xmin = fi[:, :, 0, iy].T   # (N, ng)
        incoming_xmin = bc_xmin.apply_to_incoming(outgoing_xmin, quad)
        for n in range(quad.N):
            if mu_x[n] > 1e-15:
                fi[:, n, 0, iy] = incoming_xmin[n]

        # Right face (x = xmax): outgoing = mu_x > 0, incoming = mu_x < 0.
        outgoing_xmax = fi[:, :, -1, iy].T  # (N, ng)
        incoming_xmax = bc_xmax.apply_to_incoming(outgoing_xmax, quad)
        for n in range(quad.N):
            if mu_x[n] < -1e-15:
                fi[:, n, -1, iy] = incoming_xmax[n]

    # ── Y-axis BC fills ───────────────────────────────────────────────
    for ix in range(nx):
        # Bottom face (y = ymin): outgoing = mu_y < 0, incoming = mu_y > 0.
        outgoing_ymin = fi[:, :, ix, 0].T   # (N, ng)
        incoming_ymin = bc_ymin.apply_to_incoming(outgoing_ymin, quad)
        for n in range(quad.N):
            if mu_y[n] > 1e-15:
                fi[:, n, ix, 0] = incoming_ymin[n]

        # Top face (y = ymax): outgoing = mu_y > 0, incoming = mu_y < 0.
        outgoing_ymax = fi[:, :, ix, -1].T  # (N, ng)
        incoming_ymax = bc_ymax.apply_to_incoming(outgoing_ymax, quad)
        for n in range(quad.N):
            if mu_y[n] < -1e-15:
                fi[:, n, ix, -1] = incoming_ymax[n]

    return fi


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
    *,
    bc_xmin: "BoundaryOperator | None" = None,
    bc_xmax: "BoundaryOperator | None" = None,
    bc_ymin: "BoundaryOperator | None" = None,
    bc_ymax: "BoundaryOperator | None" = None,
) -> np.ndarray:
    """Apply the streaming + collision operator T·ψ.

    Parameters
    ----------
    solution : (n_unknowns,) flattened angular flux vector.
    sig_t : (nx, ny, ng) total cross section.
    bc_xmin, bc_xmax, bc_ymin, bc_ymax :
        Wave E Round 3 — :class:`~orpheus.geometry.boundary.BoundaryOperator`
        instances threaded through to :func:`solution_to_angular_flux`
        for BC-faithful boundary fills.  ``None`` = legacy reflective
        fallback (bit-identical to pre-Round 3 hard-coded behaviour).

    Returns
    -------
    np.ndarray
        Shape ``(n_unknowns,)``. T applied to the angular flux.
    """
    fi = solution_to_angular_flux(
        solution, eq_map, quad, nx, ny, ng,
        bc_xmin=bc_xmin, bc_xmax=bc_xmax,
        bc_ymin=bc_ymin, bc_ymax=bc_ymax,
    )

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


# ═══════════════════════════════════════════════════════════════════════
# Spherical 1D operator: L = μ(A∂/∂r)/V + (α∂/∂μ)/V + Σ_t
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
    *,
    bc_outer: "BoundaryOperator | None" = None,
) -> np.ndarray:
    """Convert 1D solution vector to angular flux array (ng, N, nx, 1).

    Applies the outer-boundary BC (r = R) via the
    :class:`~orpheus.geometry.boundary.BoundaryOperator` ``bc_outer``.  The
    inner boundary at r = 0 is intrinsically symmetric (the spherical
    pole has zero face area; the matvec sets ``psi_left = 0`` there
    by construction), so no fill is needed at i = 0.

    Wave E Round 3 (ERR-026 closure): when ``bc_outer`` is ``None``
    the function falls back to specular reflection at the outer face,
    bit-identical to the pre-Round 3 hard-coded reflective fill that
    the BiCGSTAB FD path relied on.  Production callers should pass
    ``sn_mesh.bc_right`` so that vacuum / white / albedo / mixed BCs
    are honoured uniformly.
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    if bc_outer is None:
        bc_outer = SpecularBoundaryOperator(axis="x", albedo=1.0)

    fi = np.zeros((ng, quad.N, nx, 1))

    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    for k in range(eq_map.n_eq):
        fi[:, eq_map.ordinate[k], eq_map.ix[k], 0] = flux[:, k]

    # ── Outer-boundary BC fill (r = R) ────────────────────────────────
    # Build the outgoing flux at i = nx-1 indexed by all ordinates
    # (incoming-ordinate slots are still zero from np.zeros above) and
    # delegate to the BC's tensor-decomposed action.  For SpecularBoundaryOperator
    # this returns ``outgoing[ref_x[n]]`` — bit-identical to the
    # pre-Round 3 ``fi[:, n, -1, 0] = fi[:, ref_x[n], -1, 0]`` fill.
    # For VacuumBoundaryOperator it returns zeros (correct vacuum incoming).
    outgoing = fi[:, :, -1, 0].T   # (N, ng)
    incoming = bc_outer.apply_to_incoming(outgoing, quad)
    for n in range(quad.N):
        if quad.mu_x[n] < -1e-15:
            fi[:, n, -1, 0] = incoming[n]

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
    *,
    bc_outer: "BoundaryOperator | None" = None,
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

    Face fluxes are approximated by arithmetic averages of cell-centre
    values.  At the outer face (i = nx-1) the cell-centre extrapolation
    ``psi_right = fi[:, n, i, 0]`` is used uniformly: for outgoing
    ordinates (μ > 0) this is the symmetric closure's cell-center
    extrapolation; for incoming ordinates the slot was BC-filled by
    :func:`solution_to_angular_flux_spherical` (Wave E Round 3) with
    ``bc_outer.apply_to_incoming(...)``, so the same read accesses the
    BC-faithful value.  Unknown ordinates at i = nx-1 (per the eq_map)
    have μ ≥ 0; the inward-direction read at i = nx-2's outer face
    (``psi_right = 0.5*(fi[:, n, nx-2, 0] + fi[:, n, nx-1, 0])``) is
    where the BC fill enters the matvec.

    Parameters
    ----------
    bc_outer :
        :class:`~orpheus.geometry.boundary.BoundaryOperator` for the outer
        face (r = R).  ``None`` = legacy reflective fallback (bit-
        identical to the pre-Round 3 hard-coded reflective fill).
    """
    fi = solution_to_angular_flux_spherical(
        solution, eq_map, quad, nx, ng, bc_outer=bc_outer,
    )
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
        # At the outer face (i = nx-1), valid unknowns always carry
        # μ ≥ 0 by the eq_map skip rule.  The if/else branch is
        # preserved for full bit-identity to the pre-Round 3 reflective
        # path; the ``else`` (μ < 0) is dead code in normal operation
        # but the ``|μ| ≤ 1e-15`` corner case (near-zero ordinates from
        # non-GL quadratures) still falls through to it and would read
        # the partner's BC-filled value.
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


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical 1D operator: L = η(A∂/∂r)/V + (ΔA/w)(α∂/∂φ)/V + Σ_t
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
    *,
    bc_outer: "BoundaryOperator | None" = None,
) -> np.ndarray:
    r"""Apply the cylindrical transport operator T·ψ.

    Per-level azimuthal redistribution with geometry-weighted
    :math:`\Delta A / w` factor and Morel–Montry angular closure.

    Parameters
    ----------
    bc_outer :
        :class:`~orpheus.geometry.boundary.BoundaryOperator` for the outer
        face (r = R).  ``None`` = legacy reflective fallback (bit-
        identical to the pre-Round 3 hard-coded reflective fill).
        Wave E Round 3 (ERR-026 closure).
    """
    fi = solution_to_angular_flux_cylindrical(
        solution, eq_map, quad, nx, ng, bc_outer=bc_outer,
    )
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
                bc_outer=sn_mesh.bc_right,
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
                bc_outer=sn_mesh.bc_right,
            )
        # Cartesian (1-D slab or 2-D)
        return transport_operator_matvec(
            psi, eq_map, quad, self.sig_t,
            nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
            bc_xmin=sn_mesh.bc_xmin,
            bc_xmax=sn_mesh.bc_xmax,
            bc_ymin=sn_mesh.bc_ymin,
            bc_ymax=sn_mesh.bc_ymax,
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

