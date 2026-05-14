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
   slots via ``bc.apply(outgoing)`` (Wave B Issue 7
   tensor-decomposed BC algebra; Wave 9 migrated the call sites to the
   uniform 1-arg contract — the quadrature is bound on the shim at
   :meth:`SNMesh._resolve_one` time).  Bit-identity to the pre-Round 3
   reflective-only fill is preserved for :class:`SpecularBoundaryOperator` (the
   default ``BC.reflective`` factory), since
   ``SpecularBoundaryOperator(axis="x").apply(out) == out[ref_x]``;
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
    from .spatial.pole_angular_closure import PoleAngularClosure

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

    # ── Issue #168 Phase C: precomputed (cell, ordinate) → unknown index ──
    # Sparse inverse lookup for the sweep-frame matvec. Set lazily on
    # first call to ``unknowns_at_cell_for_mask`` so existing callers
    # are not forced to build it when they do not need it.
    _cell_ord_to_k: np.ndarray | None = field(default=None, repr=False)
    """``(nx, N) int`` with ``-1`` for absent (ordinate, cell) pairs.

    Populated lazily by :meth:`unknowns_at_cell_for_mask`.
    """

    def unknowns_at_cell_for_mask(
        self,
        cell_idx: int,
        ordinate_mask: np.ndarray,
    ) -> np.ndarray:
        r"""Return unknown indices ``k`` at cell ``cell_idx`` for ordinates in mask.

        Issue #168 Phase C — sweep-frame matvec helper. Built from
        :attr:`ordinate` and :attr:`ix` as a precomputed inverse
        lookup table (O(N · nx) memory, lazily allocated) so the
        per-cell call is an O(mask_count) array gather rather than
        the O(n_eq) linear scan the legacy per-equation pattern used.

        The first call materialises a ``(nx_inferred, N) int`` table
        where ``table[i, n] == k`` if unknown ``k`` is at cell ``i``
        ordinate ``n``, and ``-1`` otherwise. Subsequent calls reuse
        the table; no per-call allocation cost beyond the boolean
        mask gather.

        Parameters
        ----------
        cell_idx : int
            Cell index ``i`` in the spatial axis.
        ordinate_mask : np.ndarray
            Boolean array, shape ``(N,)``, ``True`` for ordinates to
            include.

        Returns
        -------
        np.ndarray
            ``int`` array of unknown indices (subset of ``[0, n_eq)``)
            ordered by the mask's ``True`` positions in ascending
            ordinate index. Returns an empty array if no unknown at
            ``cell_idx`` matches the mask.

        Raises
        ------
        ValueError
            If ``ordinate_mask`` is not 1-D or its length disagrees
            with the quadrature's ordinate count inferred from
            :attr:`ordinate`.
        """
        ordinate_mask = np.asarray(ordinate_mask, dtype=bool)
        if ordinate_mask.ndim != 1:
            raise ValueError(
                f"ordinate_mask must be 1-D; got shape "
                f"{ordinate_mask.shape}"
            )
        if self._cell_ord_to_k is None:
            # Infer (nx, N) from the eq_map content. Use the maximum
            # index plus one so the table covers every cell / ordinate
            # the eq_map references.
            nx = int(self.ix.max()) + 1 if self.n_eq > 0 else 0
            n_ord = int(self.ordinate.max()) + 1 if self.n_eq > 0 else 0
            n_ord = max(n_ord, ordinate_mask.size)
            table = -np.ones((nx, n_ord), dtype=np.intp)
            table[self.ix, self.ordinate] = np.arange(self.n_eq)
            self._cell_ord_to_k = table
        table = self._cell_ord_to_k
        if ordinate_mask.size != table.shape[1]:
            raise ValueError(
                f"ordinate_mask length {ordinate_mask.size} does not "
                f"match inferred N={table.shape[1]}"
            )
        ks = table[cell_idx][ordinate_mask]
        return ks[ks >= 0]


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
    ``apply(outgoing)`` method (Wave 9 1-arg contract — the quadrature
    is bound on the SNMesh-side shim) maps the boundary's outgoing
    angular flux to the incoming flux per the BC's tensor decomposition
    (vacuum → 0; specular → ``out[ref]``; white, albedo, mixed → their
    respective combinations).  When all four
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
    # the BC.  Layout: ``apply`` consumes a full
    # ``(N, ng_axis_2)`` array (here ``(N, ng)``) and returns the same
    # shape; we slice the incoming entries out for the boundary fill.
    # For each y-row independently to avoid spurious coupling.
    for iy in range(ny):
        # Left face (x = xmin): outgoing = mu_x < 0, incoming = mu_x > 0.
        outgoing_xmin = fi[:, :, 0, iy].T   # (N, ng)
        incoming_xmin = bc_xmin.apply(outgoing_xmin)
        for n in range(quad.N):
            if mu_x[n] > 1e-15:
                fi[:, n, 0, iy] = incoming_xmin[n]

        # Right face (x = xmax): outgoing = mu_x > 0, incoming = mu_x < 0.
        outgoing_xmax = fi[:, :, -1, iy].T  # (N, ng)
        incoming_xmax = bc_xmax.apply(outgoing_xmax)
        for n in range(quad.N):
            if mu_x[n] < -1e-15:
                fi[:, n, -1, iy] = incoming_xmax[n]

    # ── Y-axis BC fills ───────────────────────────────────────────────
    for ix in range(nx):
        # Bottom face (y = ymin): outgoing = mu_y < 0, incoming = mu_y > 0.
        outgoing_ymin = fi[:, :, ix, 0].T   # (N, ng)
        incoming_ymin = bc_ymin.apply(outgoing_ymin)
        for n in range(quad.N):
            if mu_y[n] > 1e-15:
                fi[:, n, ix, 0] = incoming_ymin[n]

        # Top face (y = ymax): outgoing = mu_y > 0, incoming = mu_y < 0.
        outgoing_ymax = fi[:, :, ix, -1].T  # (N, ng)
        incoming_ymax = bc_ymax.apply(outgoing_ymax)
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
) -> np.ndarray:
    r"""Convert 1D solution vector to **pure cell-centre** angular flux.

    Issue #168 Phase C — simplified to return only the cell-centre
    array. The companion ``boundary_face_flux`` output of the Phase A
    signature is RETIRED: under the sweep-frame matvec architecture
    the BC trace law is applied **inside the matvec body** on the
    WDD-propagated outflow face values, not on the cell-centre slots
    of this array. The §16A.3 contract that the inflow trace equals
    the BC operator applied to the outflow trace is now honoured by
    construction.

    Returns ``fi`` shape ``(ng, N, nx, 1)`` — cell-centre angular
    flux, DOFs scattered from ``solution`` at their
    ``(ordinate, cell)`` slots, plus an analytical extension at
    ``i = nx-1`` for inward ordinates: those slots receive the
    reflected-partner's cell-centre value. The eq_map intentionally
    excludes inward-at-boundary slots from the unknown set (the
    sweep determines them via the BC); the analytical extension is
    a smooth-flat-flux-preserving choice for the WDD recurrence in
    the sweep-frame matvec.

    History
    -------

    The Phase A signature returned ``(fi, boundary_face_flux)`` where
    the second array carried per-(group, ordinate, side) BC-resolved
    face fluxes; the matvec read the inward outer-face entries from
    this companion array. Phase C retires that two-storage design
    entirely. The §16A.3 contract requires the BC trace law to consume
    the boundary face TRACE (the WDD-propagated outflow face value),
    not interior cell-centre approximations. The sweep-frame matvec
    achieves this by computing the outflow face flux via the WDD
    diamond as it propagates outward, then calling
    ``bc_outer.apply`` ONCE at the boundary edge to produce the
    inflow face flux for the inward direction. There is no
    pre-staged BC face-value array.

    Parameters
    ----------
    solution
        Packed (ng·n_eq,) vector laid out group-major Fortran order.
    eq_map
        :class:`EquationMap` describing which (ordinate, cell) slots
        are unknowns.
    quad
        Angular quadrature (used only for ``quad.N``).
    nx, ng
        Spatial cell count + group count.

    Returns
    -------
    np.ndarray
        ``fi`` shape ``(ng, N, nx, 1)`` — cell-centre angular flux.
    """
    fi = np.zeros((ng, quad.N, nx, 1))
    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    for k in range(eq_map.n_eq):
        fi[:, eq_map.ordinate[k], eq_map.ix[k], 0] = flux[:, k]
    # Analytical extension for inward ordinates at the outer
    # boundary cell. The eq_map skipped these slots (the BC fixes
    # them); we set the cell-centre to the reflected-partner's
    # cell-centre value so the WDD recurrence on flat ψ preserves
    # the per-ordinate flat-flux invariant. For reflective BC this
    # is the exact value the inward sweep would consume; for
    # vacuum / other BCs it is a smooth extension that lets the
    # WDD diamond closure stay consistent with the BC trace law
    # applied at the face.
    ref_x = quad.reflection_index("x")
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
    *,
    sn_mesh: "SNMesh | None" = None,
    bc_outer: "BoundaryOperator | None" = None,
    pole_angular_closure: "PoleAngularClosure | None" = None,
) -> np.ndarray:
    r"""Apply the spherical transport operator T·ψ.

    .. math::

        (T\psi)_{n,i} = \frac{\mu_n}{V_i}
          \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
        + \frac{\Delta A_i}{w_n V_i}
          \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
        + \Sigma_t \psi_{n,i}

    The :math:`\Delta A / w` geometry factor (``redist_dAw``, precomputed
    in :class:`SNMesh`) carries the geometry-redistribution coefficient
    per Hébert (2009) §3.9.4 Eq. 3.428.

    Phase C — sweep-frame architecture (Issue #168)
    ------------------------------------------------

    The matvec is **one sweep iteration semantically**. For each
    bulk direction (outward μ ≥ 0 then inward μ < 0) the WDD diamond
    closure ``ψ_face_out = 2·ψ_cell - ψ_face_in`` propagates the
    face flux cell-by-cell along the direction's DAG order. The BC
    trace law owns the boundary edge: ``bc_outer.apply(outflow_face)``
    is called ONCE per matvec, consuming the WDD-propagated outflow
    face values (per the §16A.3 contract) and producing the inflow
    face values the inward sweep consumes.

    Vectorisation: the per-cell update operates on the WHOLE ordinate
    subset at once via ``outgoing_mask = quad.mu_x > +eps`` and
    ``incoming_mask = quad.mu_x < -eps``. No ``for n in range(quad.N)``
    loop appears in the sweep body — only the per-cell loop over the
    DAG-ordered cell sequence.

    Three Issue #168 defects close by construction in this rewrite:

    * Defect 1 (outer-face cell-centre truncation): WDD propagation
      gives the outflow face flux at full second order; no algebraic
      extrapolation through cell centres.
    * Defect 2 (BC-fill cell-centre contamination): the BC apply
      operates on face values directly; cell-centre storage at
      ``fi[..., -1, 0]`` is not consumed for inward ordinates.
    * Defect 3 (Bailey ΔA/w pole closure inconsistency): the
      ``PoleAngularClosure`` strategy evaluates the angular
      redistribution term consistently with the WDD spatial closure,
      so Krylov-on-apply preconditioned by the sweep converges
      cleanly (subject to the empirical Gate 1.1 outcome — see
      :file:`tests/sn/test_phase_c_gates.py`).

    Phase A's :class:`BoundaryFaceFlux` Protocol retires entirely.

    Parameters
    ----------
    sn_mesh :
        Optional :class:`SNMesh` for the new
        :meth:`SNMesh.iter_cells_by_direction` API. When ``None``,
        the function falls back to plain ``range(nx)`` iteration
        for the outward sweep and ``range(nx-1, -1, -1)`` for the
        inward sweep (the equivalent direct sweep order for
        Cartesian-style cell numbering).
    bc_outer :
        :class:`~orpheus.geometry.boundary.BoundaryOperator` for the outer
        face (r = R).  ``None`` = legacy reflective fallback.
    pole_angular_closure :
        :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
        strategy for evaluating the angular redistribution term.
        ``None`` defaults to
        :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
        (Phase B default — see :class:`SNMesh` for the rationale).
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    from .spatial.pole_angular_closure import (
        LegacyTauSymmetricInterpolation,
    )

    if bc_outer is None:
        # Match the pre-Phase-C default: specular reflection at the
        # outer face. SpecularBoundaryOperator is the deprecated alias
        # for ReflectiveBoundary; we pass through realize() to obtain
        # the same 1-arg PermutationOperator wired by SNMesh.
        from .boundary_realizer import SNBoundaryRealizer, SNMethodSpace
        spec_law = SpecularBoundaryOperator(axis="x", albedo=1.0)
        bc_outer = SNBoundaryRealizer().realize(
            spec_law, SNMethodSpace.minimal(quad),
        )
    if pole_angular_closure is None:
        pole_angular_closure = LegacyTauSymmetricInterpolation()

    fi = solution_to_angular_flux_spherical(
        solution, eq_map, quad, nx, ng,
    )
    A = face_areas       # (nx+1,)
    V = volumes[:, 0]    # (nx,)
    N = quad.N
    eps = 1e-15

    # Vectorised ordinate masks (per-user hard constraint — no
    # ``for n in range(quad.N)`` loops in Phase C-touched code).
    outgoing_mask = quad.mu_x > +eps
    incoming_mask = quad.mu_x < -eps
    mu_out = quad.mu_x[outgoing_mask]
    mu_in = quad.mu_x[incoming_mask]
    n_out = int(mu_out.size)
    n_in = int(mu_in.size)

    # ── Phase D Carlson coupled-pole sweep context ───────────────────
    # The M-M angular closure consumes this context (when supplied) to
    # build the canonical Hébert §3.9.4 Eqs. (3.432)-(3.435) inward
    # μ = −1 sweep seed for the recurrence's ``psi_half_left``.
    # Legacy / Bailey closures ignore it (Protocol-shape compatibility).
    #
    # ``bc_outer_value``: angular flux at the outer face at μ = −1,
    # derived from input ψ via the realized BC operator.  Applying the
    # BC to the cell-centred ψ at the outer cell (a first-order proxy
    # for the outflow face flux) gives the inward face flux per
    # ordinate; we extract the value at the most-inward ordinate
    # (μ closest to −1).  Linear in ψ — preserves operator linearity.
    #
    # For reflective BC + flat ψ = C: gives C ✓ (PermutationOperator
    # swaps outgoing↔incoming, both = C).
    # For vacuum BC: gives 0 ✓ (IncomingOrdinateMaskTensor zeroes
    # inflow ordinates).
    from .spatial.psi_half_angle_seed import CarlsonSweepContext
    if n_in > 0:
        # Apply BC realization to cell-centred outer-cell ψ (shape (N, ng))
        # to obtain the outer-face inflow trace.  Read the value at the
        # most-inward ordinate (μ closest to −1) per group.
        outer_inflow_estimate = bc_outer.apply(fi[:, :, -1, 0].T)  # (N, ng)
        # Most-inward ordinate index in the incoming-mask subset:
        # ordinate with smallest μ_x.
        most_inward_global_idx = int(np.argmin(quad.mu_x))
        bc_outer_value = outer_inflow_estimate[most_inward_global_idx, :]  # (ng,)
    else:
        bc_outer_value = np.zeros((ng,))

    # σ_t per (g, i): the matvec receives ``sig_t`` shape ``(nx, ny=1, ng)``.
    # Reshape to ``(ng, nx)`` for the Carlson sweep.
    sigma_t_gx = sig_t[:, 0, :].T  # (ng, nx)

    # Radial widths for the inward sweep's denominator.
    dr = sn_mesh.dx if sn_mesh is not None else np.diff(np.arange(nx + 1))

    carlson_ctx = CarlsonSweepContext(
        sigma_t=sigma_t_gx,
        dr=dr,
        mu_quad=quad.mu_x.copy(),
        weights=quad.weights.copy(),
        bc_outer_value=bc_outer_value,
    )

    # ── Phase B angular redistribution (precompute per (g, n, i)) ──
    redist_full = pole_angular_closure(
        fi[..., 0], alpha_half, redist_dAw, tau_mm, V,
        level_indices=None,  # spherical = single level
        carlson_context=carlson_ctx,
    )

    # Allocate output and the boundary-face buffer for the BC apply.
    lhs = np.empty((ng, eq_map.n_eq))
    outflow_at_boundary = np.zeros((ng, N))

    # Direction-keyed cell-visit sequences. Outward sweep iterates
    # i = 0 → nx-1; the WDD recurrence walks ψ_face_in → ψ_face_out
    # along the way. At the pole face (i=0) ψ_face = 0 by symmetry —
    # also multiplied by A[0] = 0 in the streaming term.
    if sn_mesh is not None:
        outward_visits = list(sn_mesh.iter_cells_by_direction(+1))
        inward_visits = list(sn_mesh.iter_cells_by_direction(-1))
    else:
        from .spatial.cell_update import CellVisit
        outward_visits = [
            CellVisit(cell_idx=i, streaming_terms=None,
                      face_area_downstream=0.0)
            for i in range(nx)
        ]
        inward_visits = [
            CellVisit(cell_idx=i, streaming_terms=None,
                      face_area_downstream=0.0)
            for i in range(nx - 1, -1, -1)
        ]

    # ── Phase 1: outgoing ordinates (μ > 0), i = 0 → nx-1 ─────────
    # Pole-face initial condition (Lewis-Miller §4.5; preserves the
    # per-ordinate flat-flux invariant on a homogeneous reflective
    # sphere): ψ_face_in at the pole face is the cell-centre value
    # ψ_cell[0]. This is the standard SN "Carlson starting direction"
    # initialisation for the apply matvec direction — for a true flat
    # field this gives ψ_face_out = 2·ψ_cell - ψ_cell = ψ_cell at
    # every cell, so streaming + redistribution sum to zero per
    # ordinate on flat ψ (the canonical curvilinear test).
    # The pole-face streaming contribution is multiplied by A[0]=0
    # so this choice does not introduce a spurious source there;
    # it only sets the recurrence's anchor for the cell-by-cell WDD
    # propagation across the interior.
    if n_out > 0:
        # ψ_face_in initialised at cell-centre value of the pole cell.
        if len(outward_visits) > 0:
            i0 = outward_visits[0].cell_idx
            psi_face_in = fi[:, outgoing_mask, i0, 0].copy()
        else:
            psi_face_in = np.zeros((ng, n_out))
        for visit in outward_visits:
            i = visit.cell_idx
            # (ng, n_out) cell-centre flux at outgoing ordinates.
            psi_cell = fi[:, outgoing_mask, i, 0]
            # WDD diamond closure.
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_out[None, :]
                * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
                / V[i]
            )
            redistribution = redist_full[:, outgoing_mask, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, outgoing_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out
        # The last cell's psi_face_out is the boundary outflow face
        # for the outgoing direction.
        outflow_at_boundary[:, outgoing_mask] = psi_face_out

    # ── BC trace law at boundary edge (Issue #168 Phase C fix) ────
    # The realised BC operator consumes the full (N, ng) outflow face
    # vector and returns the (N, ng) inflow face vector with the
    # inflow ordinates populated. The shim wraps a 1-arg
    # LinearOperator; its `apply` expects the leading axis to be the
    # ordinate axis (axis=0 of PermutationOperator /
    # IncomingOrdinateMaskTensor).
    inflow_full = bc_outer.apply(outflow_at_boundary.T)  # (N, ng)

    # ── Phase 2: incoming ordinates (μ < 0), i = nx-1 → 0 ─────────
    if n_in > 0:
        psi_face_in = inflow_full[incoming_mask, :].T  # (ng, n_in)
        for visit in inward_visits:
            i = visit.cell_idx
            psi_cell = fi[:, incoming_mask, i, 0]
            # WDD diamond closure (face_in is at A[i+1], face_out
            # — downstream — at A[i]).
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_in[None, :]
                * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
                / V[i]
            )
            redistribution = redist_full[:, incoming_mask, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, incoming_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out

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
    sn_mesh: "SNMesh | None" = None,
    bc_outer: "BoundaryOperator | None" = None,
    pole_angular_closure: "PoleAngularClosure | None" = None,
) -> np.ndarray:
    r"""Apply the cylindrical transport operator T·ψ.

    Per-level azimuthal redistribution with geometry-weighted
    :math:`\Delta A / w` factor and Morel–Montry angular closure.

    Phase C — sweep-frame architecture (Issue #168)
    ------------------------------------------------

    The cylindrical analogue of the sweep-frame spherical matvec.
    The cylindrical case has per-:math:`\mu`-level structure: each
    level has its own azimuthal ordinate subset with its own
    +η / -η partition. The sweep runs **per level**, with each
    level's outward (η > 0) ordinates iterating cells 0 → nx-1,
    then the BC trace law applies to the level's full outflow face
    vector, then the level's inward (η < 0) ordinates iterate cells
    nx-1 → 0.

    Pure-azimuthal degenerate ordinates (|η| < 1e-15) carry no
    spatial flow; their streaming term is zero by construction
    (mu_x ≈ 0) but they still contribute through the angular
    redistribution and collision. The matvec includes them
    explicitly via a degenerate mask handled in a third per-level
    pass.

    Three Issue #168 defects close by construction — see the
    spherical docstring for the full narrative. Phase A's
    :class:`BoundaryFaceFlux` Protocol retires entirely.

    Parameters
    ----------
    sn_mesh :
        Optional :class:`SNMesh` for the new
        :meth:`SNMesh.iter_cells_by_direction` API. When ``None``,
        the function uses ``range(nx)`` / ``range(nx-1, -1, -1)``.
    bc_outer :
        :class:`~orpheus.geometry.boundary.BoundaryOperator` for the outer
        face (r = R).  ``None`` = legacy reflective fallback.
    pole_angular_closure :
        :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
        strategy for evaluating the per-level azimuthal redistribution.
        ``None`` defaults to
        :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`.
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    from .spatial.pole_angular_closure import (
        LegacyTauSymmetricInterpolation,
    )

    if bc_outer is None:
        from .boundary_realizer import SNBoundaryRealizer, SNMethodSpace
        spec_law = SpecularBoundaryOperator(axis="x", albedo=1.0)
        bc_outer = SNBoundaryRealizer().realize(
            spec_law, SNMethodSpace.minimal(quad),
        )
    if pole_angular_closure is None:
        pole_angular_closure = LegacyTauSymmetricInterpolation()

    fi = solution_to_angular_flux_cylindrical(
        solution, eq_map, quad, nx, ng,
    )
    A = face_areas       # (nx+1,)
    V = volumes[:, 0]    # (nx,)
    N = quad.N
    mu = quad.mu_x
    eps = 1e-15

    # ── Phase D Carlson coupled-pole sweep context (per μ-level) ────
    # Cylindrical analog of the spherical Carlson context (see the
    # spherical matvec for the architectural rationale).  Each
    # μ-level has its own azimuthal subset of ordinates; the Carlson
    # context is built per-level so the M-M recurrence's seed is the
    # canonical Hébert form for each level's auxiliary inward
    # direction.  Per the Phase D diagnosis memo §7 observation 3,
    # cylindrical Gate 1.1 already passes empirically (per-level
    # α-domes telescope) — applying the Carlson seed per-level aligns
    # the architecture with the canonical form without regressing
    # the empirical pass.
    from .spatial.psi_half_angle_seed import CarlsonSweepContext

    # σ_t per (g, i).
    sigma_t_gx = sig_t[:, 0, :].T  # (ng, nx)
    dr = sn_mesh.dx if sn_mesh is not None else np.diff(np.arange(nx + 1))

    # Pre-compute outer-face inflow estimate (full ordinate vector)
    # via the BC realization on cell-centred ψ at the outer cell.
    # Shape ``(N, ng)``.
    outer_inflow_estimate = bc_outer.apply(fi[:, :, -1, 0].T)

    carlson_ctx_per_level: list[CarlsonSweepContext] = []
    for level_idx in quad.level_indices:
        level_idx_arr = np.asarray(level_idx)
        mu_level = mu[level_idx_arr]
        weights_level = quad.weights[level_idx_arr]
        # Most-inward ordinate within the level (smallest μ).
        # Within-level index of the most-inward ordinate.
        within_idx_most_inward = int(np.argmin(mu_level))
        global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
        bc_outer_value_level = outer_inflow_estimate[global_idx_most_inward, :]  # (ng,)
        carlson_ctx_per_level.append(
            CarlsonSweepContext(
                sigma_t=sigma_t_gx,
                dr=dr,
                mu_quad=mu_level.copy(),
                weights=weights_level.copy(),
                bc_outer_value=bc_outer_value_level,
            )
        )

    # ── Phase B per-level angular redistribution ──────────────────
    redist_full = pole_angular_closure(
        fi[..., 0],
        alpha_per_level,
        redist_dAw_per_level,
        tau_mm_per_level,
        V,
        level_indices=quad.level_indices,
        carlson_context=carlson_ctx_per_level,
    )

    lhs = np.empty((ng, eq_map.n_eq))
    outflow_at_boundary = np.zeros((ng, N))

    level_indices = quad.level_indices

    # ── Per-level outward sweep + BC apply + inward sweep ─────────
    # Phase 1: outward sweep across every level, populating the
    # boundary outflow face vector.
    for level_p, level_idx in enumerate(level_indices):
        level_idx = np.asarray(level_idx)
        eta_level = mu[level_idx]
        out_within = eta_level > +eps
        if not np.any(out_within):
            continue
        # Global ordinate indices for the outward subset of this
        # level.
        global_out = level_idx[out_within]
        mu_out = mu[global_out]
        n_out = int(global_out.size)
        # Use a representative outward ordinate for the cell ordering.
        if sn_mesh is not None:
            visits = list(
                sn_mesh.iter_cells_by_direction(+1, mu_level_idx=level_p)
            )
        else:
            from .spatial.cell_update import CellVisit
            visits = [
                CellVisit(cell_idx=i, streaming_terms=None,
                          face_area_downstream=None)
                for i in range(nx)
            ]

        # Pole-face initial condition: cell-centre value (Lewis-
        # Miller §4.5; preserves flat-flux invariant — see the
        # spherical matvec for the full rationale).
        if len(visits) > 0:
            i0 = visits[0].cell_idx
            psi_face_in = fi[:, global_out, i0, 0].copy()
        else:
            psi_face_in = np.zeros((ng, n_out))
        for visit in visits:
            i = visit.cell_idx
            psi_cell = fi[:, global_out, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_out[None, :]
                * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
                / V[i]
            )
            global_out_mask = np.zeros(N, dtype=bool)
            global_out_mask[global_out] = True
            redistribution = redist_full[:, global_out, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, global_out_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out
        # Boundary outflow face for this level's outward ordinates.
        outflow_at_boundary[:, global_out] = psi_face_out

    # ── BC trace law at boundary edge (single call across all levels) ──
    inflow_full = bc_outer.apply(outflow_at_boundary.T)  # (N, ng)

    # Phase 2: inward sweep across every level, consuming the BC's
    # inflow face values.
    for level_p, level_idx in enumerate(level_indices):
        level_idx = np.asarray(level_idx)
        eta_level = mu[level_idx]
        in_within = eta_level < -eps
        if not np.any(in_within):
            continue
        global_in = level_idx[in_within]
        mu_in = mu[global_in]
        n_in = int(global_in.size)
        if sn_mesh is not None:
            visits = list(
                sn_mesh.iter_cells_by_direction(-1, mu_level_idx=level_p)
            )
        else:
            from .spatial.cell_update import CellVisit
            visits = [
                CellVisit(cell_idx=i, streaming_terms=None,
                          face_area_downstream=None)
                for i in range(nx - 1, -1, -1)
            ]

        psi_face_in = inflow_full[global_in, :].T  # (ng, n_in)
        for visit in visits:
            i = visit.cell_idx
            psi_cell = fi[:, global_in, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_in[None, :]
                * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
                / V[i]
            )
            global_in_mask = np.zeros(N, dtype=bool)
            global_in_mask[global_in] = True
            redistribution = redist_full[:, global_in, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            ks = eq_map.unknowns_at_cell_for_mask(i, global_in_mask)
            if ks.size > 0:
                lhs[:, ks] = streaming + redistribution + collision
            psi_face_in = psi_face_out

    # Phase 3 (cylindrical pure-azimuthal degenerate ordinates): |η| < eps.
    # Spatial flow is zero (μ_x ≈ 0). The streaming term is zero by
    # construction, so only the redistribution + collision contribute.
    degenerate_mask = np.abs(mu) < eps
    if np.any(degenerate_mask):
        global_deg = np.where(degenerate_mask)[0]
        for i in range(nx):
            psi_cell = fi[:, global_deg, i, 0]
            redistribution = redist_full[:, global_deg, i]
            collision = sig_t[i, 0, :, None] * psi_cell
            # Streaming is identically zero for |η| < eps.
            ks = eq_map.unknowns_at_cell_for_mask(i, degenerate_mask)
            if ks.size > 0:
                lhs[:, ks] = redistribution + collision

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
                sn_mesh=sn_mesh,
                bc_outer=sn_mesh.bc_right,
                pole_angular_closure=sn_mesh.pole_angular_closure,
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
                sn_mesh=sn_mesh,
                bc_outer=sn_mesh.bc_right,
                pole_angular_closure=sn_mesh.pole_angular_closure,
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

