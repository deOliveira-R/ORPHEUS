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

    from .boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from .sources import IsotropicSource, PerOrdinateSource
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
    "StreamingOperator",
    "CollisionOperator",
]


# ═══════════════════════════════════════════════════════════════════════
# Equation map: which (ordinate, cell) pairs are unknowns
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class EquationMap:
    """Mapping between the packed 1-D FD-matvec solution vector and the
    4-D angular flux ``(N, ng, nx, ny)``.

    The packed solution is ``(n_unknowns,)`` where
    ``n_unknowns = n_eq * ng``.  Flat indexing under Fortran order:
    ``flat_idx = g + ng * k`` for group ``g`` and equation ``k``.  This
    is an **optimisation choice of the FD-matvec sparse-matrix CSR row
    order** — it is **orthogonal** to the user-visible field convention
    ``(N, ng, nx, ny)`` (which lives at the
    :func:`solution_to_angular_flux` decoder output and propagates
    end-to-end through every public API per
    :ref:`theory-sn-index-convention`).

    The packed vector is purely an internal contract between
    :func:`build_equation_map` (CSR row-order producer) and the
    ``transport_operator_matvec_*`` family (CSR row-order consumer);
    it never crosses a public boundary.  The Fortran-flat ``(g, k)``
    layout was retained at PR-INDEX-7 §B.4 because flipping it
    multiplies sparse-matrix blast radius without changing user-visible
    behaviour.

    Attributes
    ----------
    n_eq : int
        Number of angular unknowns per group.  Equals ``N * nx * ny``
        for Cartesian geometries; less for curvilinear geometries where
        inward-at-outer-boundary ordinate-cell slots are BC-resolved
        rather than solved.
    n_unknowns : int
        Total scalar unknowns ``= n_eq * ng``.
    ordinate, ix, iy : ndarray of shape ``(n_eq,)``
        Per-equation indices into the ``(N, ng, nx, ny)`` angular flux.
        ``solution_to_angular_flux*`` scatters: ``fi[ordinate[k], g,
        ix[k], iy[k]] = sol[g + ng * k]``.
    """

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
    """Convert 1D solution vector to 4D angular flux **(N, ng, nx, ny)**.

    Applies z-hemisphere reflection (Lebedev) and BC-resolved fills at
    the four boundaries of the 2-D Cartesian domain.

    Issue #196 PR-INDEX-7 — output layout is principled
    ``(N, ng, nx, ny)`` (n leading, g second, then mesh axes) per
    :ref:`theory-sn-index-convention`. Was the FD-matvec internal
    layout ``(ng, N, nx, ny)`` through PR-INDEX-6; the flip closes the
    PR-INDEX-4 §9.1 deferral. The packed-vector traversal in
    :func:`build_equation_map` is UNCHANGED — only the decoded 4-D
    array's axis order flips.

    Wave E Round 3 (ERR-026 closure): the four ``bc_*`` keyword
    arguments accept :class:`~orpheus.geometry.boundary.BoundaryOperator`
    instances built by :class:`~orpheus.sn.geometry.SNMesh` from the
    mesh's :class:`~orpheus.geometry.mesh.BC` declarations.  Each BC's
    ``apply(outgoing)`` method (Wave 9 1-arg contract — the quadrature
    is bound on the SNMesh-side shim) maps the boundary's outgoing
    angular flux to the incoming flux per the BC's tensor decomposition
    (vacuum → 0; specular → ``out[ref]``; white, albedo, mixed → their
    respective combinations). ``apply`` consumes ``(N, ng)`` — under
    PR-INDEX-7 the slice ``fi[:, :, 0, iy]`` IS ``(N, ng)`` natively
    (no transpose). When all four are ``None`` (legacy callers) the
    behaviour falls back to specular reflection on every face —
    bit-identical to the pre-Round 3 hard-coded reflective fill that
    the BiCGSTAB FD path relied on.
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

    # PR-INDEX-7: principled (N, ng, nx, ny) allocation.
    fi = np.zeros((quad.N, ng, nx, ny))

    # Scatter solution into fi: the packed vector reshape stays
    # ``(ng, n_eq, order='F')`` per the EquationMap traversal which
    # PR-INDEX-7 preserves (PR-INDEX-4 §9.1 deferral content). The
    # destination's axis 0 is the ordinate; the ``:`` slice covers
    # the ng axis matching ``flux[:, k]``'s shape ``(ng,)``.
    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    # Vectorised scatter: ``fi[ord_arr, :, ix_arr, iy_arr] = flux.T``
    # writes the (n_eq, ng) packed slab into the (N, ng, nx, ny) view at
    # the eq_map's unknown slots in one BLAS-friendly assignment.  Replaces
    # the per-k Python loop (L16 perf hit on Krylov-driven test suites).
    fi[eq_map.ordinate, :, eq_map.ix, eq_map.iy] = flux.T

    # Z-reflection (Lebedev / 3-D quadratures): the eq_map only carries
    # mu_z >= 0; the lower hemisphere is filled by reflection across z.
    for n in range(quad.N):
        if mu_z[n] < -1e-15:
            fi[n, :, :, :] = fi[ref_z[n], :, :, :]

    # ── X-axis BC fills ───────────────────────────────────────────────
    # The eq_map skips incoming-at-boundary slots; we fill them here per
    # the BC.  Layout: ``apply`` consumes a full ``(N, ng)`` array and
    # returns the same shape; under PR-INDEX-7 the slice
    # ``fi[:, :, 0, iy]`` IS ``(N, ng)`` natively (no transpose).
    # For each y-row independently to avoid spurious coupling.
    for iy in range(ny):
        # Left face (x = xmin): outgoing = mu_x < 0, incoming = mu_x > 0.
        outgoing_xmin = fi[:, :, 0, iy]   # (N, ng) natively under PR-INDEX-7
        incoming_xmin = bc_xmin.apply(outgoing_xmin)
        for n in range(quad.N):
            if mu_x[n] > 1e-15:
                fi[n, :, 0, iy] = incoming_xmin[n]

        # Right face (x = xmax): outgoing = mu_x > 0, incoming = mu_x < 0.
        outgoing_xmax = fi[:, :, -1, iy]  # (N, ng)
        incoming_xmax = bc_xmax.apply(outgoing_xmax)
        for n in range(quad.N):
            if mu_x[n] < -1e-15:
                fi[n, :, -1, iy] = incoming_xmax[n]

    # ── Y-axis BC fills ───────────────────────────────────────────────
    for ix in range(nx):
        # Bottom face (y = ymin): outgoing = mu_y < 0, incoming = mu_y > 0.
        outgoing_ymin = fi[:, :, ix, 0]   # (N, ng)
        incoming_ymin = bc_ymin.apply(outgoing_ymin)
        for n in range(quad.N):
            if mu_y[n] > 1e-15:
                fi[n, :, ix, 0] = incoming_ymin[n]

        # Top face (y = ymax): outgoing = mu_y > 0, incoming = mu_y < 0.
        outgoing_ymax = fi[:, :, ix, -1]  # (N, ng)
        incoming_ymax = bc_ymax.apply(outgoing_ymax)
        for n in range(quad.N):
            if mu_y[n] < -1e-15:
                fi[n, :, ix, -1] = incoming_ymax[n]

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

    Issue #196 PR-INDEX-7: consumes ``fi`` in principled
    ``(N, ng, nx, ny)`` layout. Per-(n, ix, iy) slice
    ``fi[n, :, ix, iy]`` is shape ``(ng,)``.

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
            dfix = fi[ref_x[n], :, ix, iy] - fi[ref_x[n], :, ix + 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix + 1])
        else:
            dfix = fi[n, :, ix, iy] - fi[n, :, ix - 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix - 1])
    elif mu_x[n] < -1e-15:
        if ix == nx - 1:
            dfix = fi[ref_x[n], :, ix - 1, iy] - fi[ref_x[n], :, ix, iy]
            hx = 0.5 * (dx[ix - 1] + dx[ix])
        else:
            dfix = fi[n, :, ix + 1, iy] - fi[n, :, ix, iy]
            hx = 0.5 * (dx[ix + 1] + dx[ix])
    else:
        # PR-INDEX-7: ng is fi.shape[1] (axis 1 = groups under principled layout).
        dfix = np.zeros(fi.shape[1])
        hx = 1.0

    # Y gradient
    if mu_y[n] > 1e-15:
        if iy == 0:
            dfiy = fi[ref_y[n], :, ix, iy] - fi[ref_y[n], :, ix, iy + 1]
            hy = 0.5 * (dy[iy] + dy[iy + 1])
        else:
            dfiy = fi[n, :, ix, iy] - fi[n, :, ix, iy - 1]
            hy = 0.5 * (dy[iy] + dy[iy - 1])
    elif mu_y[n] < -1e-15:
        if iy == ny - 1:
            dfiy = fi[ref_y[n], :, ix, iy - 1] - fi[ref_y[n], :, ix, iy]
            hy = 0.5 * (dy[iy - 1] + dy[iy])
        else:
            dfiy = fi[n, :, ix, iy + 1] - fi[n, :, ix, iy]
            hy = 0.5 * (dy[iy + 1] + dy[iy])
    else:
        dfiy = np.zeros(fi.shape[1])
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
    sig_t : (ng, nx, ny) total cross section (Issue #196 PR-INDEX-3).
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
        # PR-INDEX-7: fi is (N, ng, nx, ny); per-(n, ix, iy) slice is (ng,).
        lhs[:, k] = (
            quad.mu_x[n] * dfidx
            + quad.mu_y[n] * dfidy
            + sig_t[:, ix, iy] * fi[n, :, ix, iy]
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
        ``fi`` shape **``(N, ng, nx, 1)``** — cell-centre angular
        flux in principled axis order (Issue #196 PR-INDEX-7; was
        ``(ng, N, nx, 1)`` through PR-INDEX-6). The packed-vector
        traversal in :func:`build_equation_map_spherical` is
        unchanged; only the decoded 4-D array's axis order flips.
    """
    # PR-INDEX-7: principled (N, ng, nx, 1) allocation.
    fi = np.zeros((quad.N, ng, nx, 1))
    flux = solution.reshape(ng, eq_map.n_eq, order='F')
    # Vectorised scatter (L16 perf): n_eq Python iterations → one fancy-
    # index assignment.  ``fi[ord_arr, :, ix_arr, 0]`` has shape
    # ``(n_eq, ng)``; ``flux.T`` matches.
    fi[eq_map.ordinate, :, eq_map.ix, 0] = flux.T
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
            fi[n, :, -1, 0] = fi[ref_x[n], :, -1, 0]
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
        Optional :class:`SNMesh` for the
        :meth:`SNMesh.dag_walk` API. When ``None``,
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
        :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
        — the Phase D canonical closure that matches the sweep route
        AND exposes ``compute_psi_half_per_level`` for the upcoming
        unified matvec (Issue #197 PR-TYPED-6c). In practice the
        fallback is unreachable: all production and test callers
        pass ``sn_mesh.pole_angular_closure`` explicitly, and
        :class:`SNMesh` itself defaults to ``MorelMontryAngularSweep``.
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    from .spatial.pole_angular_closure import (
        MorelMontryAngularSweep,
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
        pole_angular_closure = MorelMontryAngularSweep()

    fi = solution_to_angular_flux_spherical(
        solution, eq_map, quad, nx, ng,
    )
    # PR-INDEX-7: ``fi`` is principled (N, ng, nx, 1) — the public
    # return of ``solution_to_angular_flux_spherical``. The matvec
    # body algebra below (per-(mask, i) slices, ``pole_angular_closure``
    # input shape, ``outflow_at_boundary`` buffer) is structured around
    # the ``(ng, N, ...)`` ordering. We compute one named-intermediate
    # transpose view here — the matvec-internal layout — and consume
    # that view through the rest of the body. NO data copy: numpy
    # transpose returns a view. The PUBLIC interface flip stays at
    # ``solution_to_angular_flux_spherical``; the matvec's INTERNAL
    # contract with ``pole_angular_closure`` (and its existing
    # ``(ng, N, nx)`` consumers in spatial/) is preserved per the same
    # §B.4 rationale that preserved the EquationMap traversal order.
    psi_g_first = fi.transpose(1, 0, 2, 3)  # (ng, N, nx, 1) view
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
        # Apply BC realization to cell-centred outer-cell ψ.  Under
        # PR-INDEX-7 the slice ``fi[:, :, -1, 0]`` IS ``(N, ng)``
        # natively (no transpose).
        outer_inflow_estimate = bc_outer.apply(fi[:, :, -1, 0])  # (N, ng)
        # Most-inward ordinate index in the incoming-mask subset:
        # ordinate with smallest μ_x.
        most_inward_global_idx = int(np.argmin(quad.mu_x))
        bc_outer_value = outer_inflow_estimate[most_inward_global_idx, :]  # (ng,)
    else:
        bc_outer_value = np.zeros((ng,))

    # σ_t per (g, i): the matvec receives ``sig_t`` shape ``(ng, nx, 1)``
    # (Issue #196 PR-INDEX-3).  Slice drops the degenerate ``ny`` axis.
    sigma_t_gx = sig_t[:, :, 0]  # (ng, nx)

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
    # ``pole_angular_closure`` consumes ``psi_cells`` in (ng, N, nx)
    # layout — we feed the matvec-internal view.
    redist_full = pole_angular_closure(
        psi_g_first[..., 0], alpha_half, redist_dAw, tau_mm, V,
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
        outward_visits = list(sn_mesh.dag_walk(direction_sign=+1))
        inward_visits = list(sn_mesh.dag_walk(direction_sign=-1))
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
        # PR-INDEX-7: use the matvec-internal ``psi_g_first`` view so
        # the (ng, n_out) algebra below stays bit-exact.
        if len(outward_visits) > 0:
            i0 = outward_visits[0].cell_idx
            psi_face_in = psi_g_first[:, outgoing_mask, i0, 0].copy()
        else:
            psi_face_in = np.zeros((ng, n_out))
        for visit in outward_visits:
            i = visit.cell_idx
            # (ng, n_out) cell-centre flux at outgoing ordinates.
            psi_cell = psi_g_first[:, outgoing_mask, i, 0]
            # WDD diamond closure.
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_out[None, :]
                * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
                / V[i]
            )
            redistribution = redist_full[:, outgoing_mask, i]
            collision = sig_t[:, i, 0, None] * psi_cell
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
            psi_cell = psi_g_first[:, incoming_mask, i, 0]
            # WDD diamond closure (face_in is at A[i+1], face_out
            # — downstream — at A[i]).
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_in[None, :]
                * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
                / V[i]
            )
            redistribution = redist_full[:, incoming_mask, i]
            collision = sig_t[:, i, 0, None] * psi_cell
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
        Optional :class:`SNMesh` for the
        :meth:`SNMesh.dag_walk` API. When ``None``,
        the function uses ``range(nx)`` / ``range(nx-1, -1, -1)``.
    bc_outer :
        :class:`~orpheus.geometry.boundary.BoundaryOperator` for the outer
        face (r = R).  ``None`` = legacy reflective fallback.
    pole_angular_closure :
        :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
        strategy for evaluating the per-level azimuthal redistribution.
        ``None`` defaults to
        :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
        — the Phase D canonical closure that matches the sweep route
        AND exposes ``compute_psi_half_per_level`` for the upcoming
        unified matvec (Issue #197 PR-TYPED-6c). In practice the
        fallback is unreachable: all production and test callers
        pass ``sn_mesh.pole_angular_closure`` explicitly, and
        :class:`SNMesh` itself defaults to ``MorelMontryAngularSweep``.
    """
    from orpheus.geometry.boundary import SpecularBoundaryOperator

    from .spatial.pole_angular_closure import (
        MorelMontryAngularSweep,
    )

    if bc_outer is None:
        from .boundary_realizer import SNBoundaryRealizer, SNMethodSpace
        spec_law = SpecularBoundaryOperator(axis="x", albedo=1.0)
        bc_outer = SNBoundaryRealizer().realize(
            spec_law, SNMethodSpace.minimal(quad),
        )
    if pole_angular_closure is None:
        pole_angular_closure = MorelMontryAngularSweep()

    fi = solution_to_angular_flux_cylindrical(
        solution, eq_map, quad, nx, ng,
    )
    # PR-INDEX-7: ``fi`` is principled (N, ng, nx, 1). The matvec
    # body algebra below uses (ng, n_mask) per-cell slices; we hold
    # one named-intermediate transpose view to preserve the body's
    # internal ordering. NO data copy — numpy transpose returns a
    # view. Same §B.4 rationale as the spherical matvec.
    psi_g_first = fi.transpose(1, 0, 2, 3)  # (ng, N, nx, 1) view
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

    # σ_t per (g, i): receives ``sig_t`` shape ``(ng, nx, 1)``
    # (Issue #196 PR-INDEX-3) — slice drops the degenerate ``ny`` axis.
    sigma_t_gx = sig_t[:, :, 0]  # (ng, nx)
    dr = sn_mesh.dx if sn_mesh is not None else np.diff(np.arange(nx + 1))

    # Pre-compute outer-face inflow estimate (full ordinate vector)
    # via the BC realization on cell-centred ψ at the outer cell.
    # Under PR-INDEX-7 ``fi[:, :, -1, 0]`` IS ``(N, ng)`` natively.
    outer_inflow_estimate = bc_outer.apply(fi[:, :, -1, 0])  # (N, ng)

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
    # ``pole_angular_closure`` consumes ``psi_cells`` in (ng, N, nx)
    # layout — feed the matvec-internal transpose view.
    redist_full = pole_angular_closure(
        psi_g_first[..., 0],
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
                sn_mesh.dag_walk(
                    direction_sign=+1, mu_level_idx=level_p,
                )
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
        # PR-INDEX-7: use ``psi_g_first`` view for the (ng, n_out)
        # algebra.
        if len(visits) > 0:
            i0 = visits[0].cell_idx
            psi_face_in = psi_g_first[:, global_out, i0, 0].copy()
        else:
            psi_face_in = np.zeros((ng, n_out))
        for visit in visits:
            i = visit.cell_idx
            psi_cell = psi_g_first[:, global_out, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_out[None, :]
                * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
                / V[i]
            )
            global_out_mask = np.zeros(N, dtype=bool)
            global_out_mask[global_out] = True
            redistribution = redist_full[:, global_out, i]
            collision = sig_t[:, i, 0, None] * psi_cell
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
                sn_mesh.dag_walk(
                    direction_sign=-1, mu_level_idx=level_p,
                )
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
            psi_cell = psi_g_first[:, global_in, i, 0]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            streaming = (
                mu_in[None, :]
                * (A[i + 1] * psi_face_in - A[i] * psi_face_out)
                / V[i]
            )
            global_in_mask = np.zeros(N, dtype=bool)
            global_in_mask[global_in] = True
            redistribution = redist_full[:, global_in, i]
            collision = sig_t[:, i, 0, None] * psi_cell
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
            psi_cell = psi_g_first[:, global_deg, i, 0]
            redistribution = redist_full[:, global_deg, i]
            collision = sig_t[:, i, 0, None] * psi_cell
            # Streaming is identically zero for |η| < eps.
            ks = eq_map.unknowns_at_cell_for_mask(i, degenerate_mask)
            if ks.size > 0:
                lhs[:, ks] = redistribution + collision

    return lhs.ravel(order='F')


# ═══════════════════════════════════════════════════════════════════════
# transport_operator_matvec_unified — Issue #197 PR-TYPED-6c Step 2+
#
# Single geometry-agnostic SN transport operator apply M(ψ; σ_t) =
# (L + C)·ψ in canonical (N, ng, nx, ny) 4-D layout. Replaces the three
# legacy helpers (transport_operator_matvec for Cartesian,
# transport_operator_matvec_spherical, transport_operator_matvec_cylindrical)
# whose per-cell math was a Pattern-2 twin of the sweep route's
# DiamondDifference.update / .residual algebra.
#
# The unified body delegates the per-cell algebra to
# cell_balance_for_streaming (PR-TYPED-6a foundation primitive),
# vectorised over the ordinate-mask axis (n_mask=N for matvec), and
# consumes MorelMontryAngularSweep.compute_psi_half_per_level via the
# _MMHalfGrid typed accessor (PR-TYPED-6b / Step 1.5) for the M-M
# upstream half-angle flux.
#
# Step-by-step rollout:
#   Step 2 (this commit): sphere implementation verified bit-exact
#       against transport_operator_matvec_spherical at ULP scale.
#       Cylinder + Cartesian raise NotImplementedError (Steps 3+4).
#   Step 3: extend to cylindrical (per-level outward / BC / inward /
#       degenerate structure) and verify against
#       transport_operator_matvec_cylindrical.
#   Step 4: extend to Cartesian (WDD via cell_balance_for_streaming;
#       differs from legacy FD _compute_gradients by an O(h) order-
#       of-accuracy delta — characterized, NOT bit-exact).
#   Step 5: wire StreamingOperator.apply to call the unified function;
#       legacy helpers retire in Step 7.
# ═══════════════════════════════════════════════════════════════════════


def transport_operator_matvec_unified(
    psi_view: np.ndarray,            # (N, ng, nx, ny) canonical layout
    sn_mesh: "SNMesh",
    sigma_t: np.ndarray,             # (ng, nx, ny)
    *,
    bc_outer: "BoundaryOperator | None" = None,
    pole_angular_closure: "PoleAngularClosure | None" = None,
) -> np.ndarray:                     # (N, ng, nx, ny) canonical layout
    r"""Unified geometry-agnostic SN transport operator apply.

    Computes :math:`M(\psi; \Sigma_t) = (L + C)\,\psi` over all
    geometries (sphere, cylinder, slab) via the shared per-cell
    algebra in :func:`cell_balance_for_streaming` and (for curvilinear)
    :meth:`MorelMontryAngularSweep.compute_psi_half_per_level`.

    Issue #197 PR-TYPED-6c — replaces the three legacy geometry-keyed
    helpers (Pattern 2: single source of truth for the discretisation
    algebra). The body is geometry-agnostic by data: the per-cell
    algebra reads cell geometry from :class:`StreamingTerms` (carried
    on each :class:`CellVisit` yielded by :meth:`SNMesh.dag_walk`) and
    the per-ordinate angular closure from
    :class:`~orpheus.sn.spatial.pole_angular_closure._MMHalfGrid`.

    Parameters
    ----------
    psi_view :
        Angular flux in canonical layout ``(N, ng, nx, ny)``. The
        matvec body uses a zero-copy transpose view to ``(ng, N, nx, ny)``
        for the per-cell ``(ng, n_mask)`` algebra, but the public
        contract is canonical.
    sn_mesh :
        :class:`SNMesh` carrying the per-cell geometry, quadrature,
        boundary realizations, and angular closure strategy. The
        function reads ``sn_mesh.quad``, ``sn_mesh.reduced``,
        ``sn_mesh.volumes``, ``sn_mesh.dx``, ``sn_mesh.bc_right``,
        ``sn_mesh.pole_angular_closure``.
    sigma_t :
        Per-group per-cell total cross section, shape ``(ng, nx, ny)``.
        Issue #196 PR-INDEX-3 — group-leading.
    bc_outer :
        Outer-face boundary operator. ``None`` (default) reads
        ``sn_mesh.bc_right``.
    pole_angular_closure :
        Angular closure strategy. ``None`` (default) reads
        ``sn_mesh.pole_angular_closure``.

    Returns
    -------
    np.ndarray
        Operator action ``M(ψ; σ_t) = (L + C)·ψ`` in canonical
        ``(N, ng, nx, ny)`` layout. Same shape as ``psi_view``.

    Raises
    ------
    NotImplementedError
        Cartesian curvature is Step 4.  Sphere + cylinder share a single
        per-level body — sphere is the M_p=N, n_levels=1 case of the
        cylinder algebra (the curvature read is a single data-only
        normalisation at the boundary).
    """
    from .spatial.cell_balance import cell_balance_for_streaming
    from .spatial.pole_angular_closure import MorelMontryAngularSweep
    from .spatial.psi_half_angle_seed import CarlsonSweepContext

    quad = sn_mesh.quad
    N = quad.N
    ng = psi_view.shape[1]
    nx = sn_mesh.nx
    ny = sn_mesh.ny
    eps = 1e-15
    curvature_raw = getattr(sn_mesh, "curvature", None)
    # ``SNMesh.curvature`` is ``None`` for Cartesian meshes.  Normalize.
    curvature = curvature_raw if curvature_raw is not None else "cartesian"

    if curvature not in ("spherical", "cylindrical", "cartesian"):
        raise ValueError(f"Unknown curvature: {curvature!r}")
    if curvature == "cartesian" and ny > 1:
        raise NotImplementedError(
            "transport_operator_matvec_unified 2-D Cartesian is not yet "
            "wired through dag_walk; only 1-D slab (ny=1) is implemented. "
            "2-D anti-diagonal wavefront sweeps remain on the legacy path."
        )

    if bc_outer is None:
        bc_outer = sn_mesh.bc_right
    bc_inner = sn_mesh.bc_left if curvature == "cartesian" else None
    if pole_angular_closure is None:
        pole_angular_closure = sn_mesh.pole_angular_closure
    if pole_angular_closure is None and curvature != "cartesian":
        pole_angular_closure = MorelMontryAngularSweep()

    # ── Per-level data sourced from the angular closure strategy ─────
    # Sphere, cylinder, slab all flow through the same per-level body.
    # Slab and sphere are the 1-level case (level_indices = (arange(N),),
    # M_p = N); cylinder iterates over its μ-levels.  The strategy owns
    # the per-level partition and the M-M coefficients (α, τ, ΔA/w, and
    # the derived c_in / c_out); the matvec body reads them via the
    # canonical tuple-of-arrays interface (PR-TYPED-6.5 Phase 3a.3).
    # ``IdentityAngularClosure`` (Cartesian) ships neutral values
    # (α = 0, τ = 1, ΔA/w = 0, c_in = c_out = 0) so the per-cell algebra
    # collapses to the slab form 2|μ|·1 + Σ_t·V (Step 2.5 principle —
    # ``cell_balance_for_streaming`` is geometry-blind by data).
    # ``sn_mesh.areas`` returns the geometry's face-area array
    # (Cartesian: ones; cylinder: 2πr; sphere: 4πr²) — Phase 1
    # canonicalised this on ``SNMesh`` directly.
    mu_x = quad.mu_x
    has_angular_closure = (curvature != "cartesian")
    level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
    alpha_per_level: tuple[np.ndarray, ...] = pole_angular_closure._alpha_per_level
    redist_dAw_per_level: tuple[np.ndarray, ...] = pole_angular_closure._dAw_per_level
    tau_mm_per_level: tuple[np.ndarray, ...] = pole_angular_closure._tau_per_level
    c_in_per_level: tuple[np.ndarray, ...] = pole_angular_closure._c_in_per_level
    c_out_per_level: tuple[np.ndarray, ...] = pole_angular_closure._c_out_per_level
    A = sn_mesh.areas                                                # (nx+1,)

    # ── Internal view: (ng, N, nx, ny) — group-leading for the
    # (ng, n_mask) per-cell algebra cell_balance_for_streaming consumes.
    # Zero-copy transpose; public interface stays canonical.
    psi_g_first = psi_view.transpose(1, 0, 2, 3)                     # (ng, N, nx, ny)
    out_g_first = np.zeros((ng, N, nx, ny))

    V = sn_mesh.volumes[:, 0]                                        # (nx,)
    sigma_t_gx = sigma_t[:, :, 0]                                    # (ng, nx)
    dr = sn_mesh.dx

    # ── Phase 1 spatial-upstream seed at the inner boundary ──────────
    # The predicate is structural, not curvature-keyed: ``bc_inner is
    # None`` means the inner edge is a pole (sphere / cylinder solid at
    # r=0 — Lewis-Miller §4.5: r=0 is the geometric pole, not a real
    # face, so the cell-centre proxy preserves the per-ordinate
    # flat-flux invariant). ``bc_inner is not None`` means there IS a
    # real inner boundary (slab at x=0, future hollow sphere /
    # annulus at r_inner), so apply ``bc_inner`` to the cell-centre ψ
    # at i=0 (symmetric to the ``bc_outer.apply`` on the outer side).
    # PR-TYPED-6.5 Phase 3a.1 retired the prior ``curvature ==
    # "cartesian"`` dispatch at this site.
    if bc_inner is None:
        pole_face_seed = psi_view[:, :, 0, 0].copy()                 # (N, ng)
    else:
        pole_face_seed = bc_inner.apply(psi_view[:, :, 0, 0])        # (N, ng)

    # ── Per-level Carlson contexts + _MMHalfGrid (curvilinear only) ──
    # Slab has no M-M angular closure (no curvature → ΔA/w = 0); skip
    # the entire half-grid computation and pass psi_angular_upstream=None
    # in the per-cell call.  cell_balance_for_streaming handles this
    # branch (slab's redistribution and angular-upstream terms vanish
    # naturally — no algebraic shortcut, just neutral data).
    half_grid_per_level: list | None = None
    if has_angular_closure:
        outer_inflow_estimate = bc_outer.apply(psi_view[:, :, -1, 0])    # (N, ng)
        carlson_ctx_per_level: list[CarlsonSweepContext] = []
        for level_idx in level_indices:
            level_idx_arr = np.asarray(level_idx)
            mu_level = mu_x[level_idx_arr]
            weights_level = quad.weights[level_idx_arr]
            within_idx_most_inward = int(np.argmin(mu_level))
            global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
            bc_outer_value_level = outer_inflow_estimate[global_idx_most_inward, :]
            carlson_ctx_per_level.append(
                CarlsonSweepContext(
                    sigma_t=sigma_t_gx,
                    dr=dr,
                    mu_quad=mu_level.copy(),
                    weights=weights_level.copy(),
                    bc_outer_value=bc_outer_value_level,
                )
            )

        # _MMHalfGrid (Pattern 4 typed accessor) per level.
        # Shape ``(ng, M_p, nx)`` per level via ``.upstream_per_ordinate``.
        psi_cells_g_first = psi_g_first[..., 0]                          # (ng, N, nx)
        half_grid_per_level = []
        for p, level_idx in enumerate(level_indices):
            level_idx_arr = np.asarray(level_idx)
            psi_level = psi_cells_g_first[:, level_idx_arr, :]           # (ng, M_p, nx)
            hg = pole_angular_closure.compute_psi_half_per_level(
                psi_level, tau_mm_per_level[p],
                carlson_context=carlson_ctx_per_level[p],
            )
            half_grid_per_level.append(hg)

    # Per-level M-M closure constants (c_in, c_out within-level) are
    # precomputed at strategy construction (see ``MorelMontryAngularSweep``
    # / ``IdentityAngularClosure``).  The matvec reads them via the
    # ``_c_in_per_level`` / ``_c_out_per_level`` tuples bound above.
    # For ordinate m at level p (curvilinear): α_in = α_half[m],
    # α_out = α_half[m+1], c_out_m = α_out / τ_m,
    # c_in_m = (1 − τ_m)/τ_m · α_out + α_in.  Identity ships zeros.

    outflow_at_boundary = np.zeros((ng, N))

    # ── Phase 1: outward sweep per level (μ_x > 0) ───────────────
    for p, level_idx in enumerate(level_indices):
        level_idx_arr = np.asarray(level_idx)
        mu_level = mu_x[level_idx_arr]
        outgoing_within = mu_level > +eps
        if not np.any(outgoing_within):
            continue
        global_out = level_idx_arr[outgoing_within]
        mu_out = mu_x[global_out]
        n_out_p = int(global_out.size)

        within_out_positions = np.where(outgoing_within)[0]
        c_in_sub = c_in_per_level[p][within_out_positions]
        c_out_sub = c_out_per_level[p][within_out_positions]
        redist_dAw_p = redist_dAw_per_level[p]                       # (nx, M_p)
        upstream_p_all = (
            half_grid_per_level[p].upstream_per_ordinate            # (ng, M_p, nx)
            if half_grid_per_level is not None else None
        )

        cell_indices_outward = list(
            sn_mesh.dag_walk_cell_indices(direction_sign=+1, mu_level_idx=p)
        )
        if not cell_indices_outward:
            continue
        # Spatial-upstream face seed at the inner boundary (i = i0).
        # Curvilinear: cell-centre proxy at the pole (no BC at r=0).
        # Slab: bc_xmin-applied cell-centre at x=0.
        psi_face_in = pole_face_seed[global_out, :].T                # (ng, n_out_p)

        for i in cell_indices_outward:
            psi_cell = psi_g_first[:, global_out, i, 0]
            # PR-TYPED-6.5 Phase 2.11: cell_balance_for_streaming no
            # longer accepts M-M-specific args (dA_w, c_in, c_out,
            # psi_angular_upstream).  Compute the closure's per-cell
            # contribution inline; Phase 3a delegates this to
            # ``closure.cell_contribution(...)``.
            dA_w_sub = redist_dAw_p[i, within_out_positions]    # (n_out_p,)
            angular_denom_term = dA_w_sub * c_out_sub            # (n_out_p,)
            if upstream_p_all is None:
                angular_numer_upstream = np.zeros((ng, n_out_p))
            else:
                angular_numer_upstream = (
                    dA_w_sub[None, :] * c_in_sub[None, :]
                    * upstream_p_all[:, within_out_positions, i]
                )                                                # (ng, n_out_p)

            denom, numer_upstream = cell_balance_for_streaming(
                abs_mu=mu_out,
                A_downstream=A[i + 1],
                A_total=A[i] + A[i + 1],
                total_xs=sigma_t_gx[:, i],
                volume=V[i],
                psi_face_in=psi_face_in,
                angular_denom_term=angular_denom_term,
                angular_numer_upstream=angular_numer_upstream,
            )
            m_full = (denom * psi_cell - numer_upstream) / V[i]
            out_g_first[:, global_out, i, 0] = m_full

            # WDD face propagation: ψ_face_out = 2·ψ̄ − ψ_face_in.
            psi_face_in = 2.0 * psi_cell - psi_face_in
        outflow_at_boundary[:, global_out] = psi_face_in

    # ── BC trace at outer face ───────────────────────────────
    # Apply ``bc_outer`` to the WDD-propagated outward outflow at the
    # outer face (r=R for curvilinear, x=L for slab).  The outward
    # sweep above populates ``outflow_at_boundary[:, global_out]`` with
    # the WDD-propagated face flux at the end of each level's outward
    # leg; that IS the outer-face trace.  The matvec stays linear in
    # ψ_view; the single-pass form does NOT chain inward outflow back
    # into outward (that's the SI sweep's job, not the matvec's).
    # PR-TYPED-6.5 Phase 3a.2 retired the Cartesian cell-centre-as-face
    # proxy here — free correctness gain for slab (the previous code
    # read ``psi_view[..., -1, 0]`` which is the cell-CENTRE at i=nx-1,
    # not the face flux at x=L).
    inflow_full = bc_outer.apply(outflow_at_boundary.T)              # (N, ng)

    # ── Phase 2: inward sweep per level (μ_x < 0) ────────────────
    for p, level_idx in enumerate(level_indices):
        level_idx_arr = np.asarray(level_idx)
        mu_level = mu_x[level_idx_arr]
        incoming_within = mu_level < -eps
        if not np.any(incoming_within):
            continue
        global_in = level_idx_arr[incoming_within]
        mu_in = mu_x[global_in]
        abs_mu_in = -mu_in
        n_in_p = int(global_in.size)

        within_in_positions = np.where(incoming_within)[0]
        c_in_sub = c_in_per_level[p][within_in_positions]
        c_out_sub = c_out_per_level[p][within_in_positions]
        redist_dAw_p = redist_dAw_per_level[p]
        upstream_p_all = (
            half_grid_per_level[p].upstream_per_ordinate
            if half_grid_per_level is not None else None
        )

        cell_indices_inward = list(
            sn_mesh.dag_walk_cell_indices(direction_sign=-1, mu_level_idx=p)
        )
        if not cell_indices_inward:
            continue
        psi_face_in = inflow_full[global_in, :].T                    # (ng, n_in_p)

        for i in cell_indices_inward:
            psi_cell = psi_g_first[:, global_in, i, 0]
            dA_w_sub = redist_dAw_p[i, within_in_positions]      # (n_in_p,)
            angular_denom_term = dA_w_sub * c_out_sub             # (n_in_p,)
            if upstream_p_all is None:
                angular_numer_upstream = np.zeros((ng, n_in_p))
            else:
                angular_numer_upstream = (
                    dA_w_sub[None, :] * c_in_sub[None, :]
                    * upstream_p_all[:, within_in_positions, i]
                )                                                 # (ng, n_in_p)

            denom, numer_upstream = cell_balance_for_streaming(
                abs_mu=abs_mu_in,
                A_downstream=A[i],
                A_total=A[i] + A[i + 1],
                total_xs=sigma_t_gx[:, i],
                volume=V[i],
                psi_face_in=psi_face_in,
                angular_denom_term=angular_denom_term,
                angular_numer_upstream=angular_numer_upstream,
            )
            m_full = (denom * psi_cell - numer_upstream) / V[i]
            out_g_first[:, global_in, i, 0] = m_full

            psi_face_in = 2.0 * psi_cell - psi_face_in

    # ── Phase 3: degenerate ordinates (|μ_x| < eps), no radial flow ──
    # Sphere quadratures (GL) have no exact zeros — this branch is a
    # no-op for sphere by construction (degenerate_mask is empty).
    degenerate_mask = np.abs(mu_x) < eps
    if np.any(degenerate_mask):
        global_deg = np.where(degenerate_mask)[0]
        # Map each degenerate global ordinate to its (level, within-level) position.
        deg_level: list[int] = []
        deg_within: list[int] = []
        for n_global in global_deg:
            for p, lvl in enumerate(level_indices):
                lvl_arr = np.asarray(lvl)
                pos = np.where(lvl_arr == n_global)[0]
                if pos.size > 0:
                    deg_level.append(p)
                    deg_within.append(int(pos[0]))
                    break
        n_deg = global_deg.size
        for i in range(nx):
            psi_upstream_collected = np.empty((ng, n_deg))
            dA_w_collected = np.empty(n_deg)
            c_in_collected = np.empty(n_deg)
            c_out_collected = np.empty(n_deg)
            for col_idx in range(n_deg):
                p = deg_level[col_idx]
                m = deg_within[col_idx]
                psi_upstream_collected[:, col_idx] = (
                    half_grid_per_level[p].upstream_per_ordinate[:, m, i]
                )
                dA_w_collected[col_idx] = redist_dAw_per_level[p][i, m]
                c_in_collected[col_idx] = c_in_per_level[p][m]
                c_out_collected[col_idx] = c_out_per_level[p][m]

            psi_cell = psi_g_first[:, global_deg, i, 0]              # (ng, n_deg)
            angular_denom_term = dA_w_collected * c_out_collected     # (n_deg,)
            angular_numer_upstream = (
                dA_w_collected[None, :] * c_in_collected[None, :]
                * psi_upstream_collected
            )                                                         # (ng, n_deg)
            denom, numer_upstream = cell_balance_for_streaming(
                abs_mu=np.abs(mu_x[global_deg]),
                A_downstream=0.0,
                A_total=A[i] + A[i + 1],
                total_xs=sigma_t_gx[:, i],
                volume=V[i],
                psi_face_in=np.zeros((ng, n_deg)),
                angular_denom_term=angular_denom_term,
                angular_numer_upstream=angular_numer_upstream,
            )
            m_full = (denom * psi_cell - numer_upstream) / V[i]
            out_g_first[:, global_deg, i, 0] = m_full

    return out_g_first.transpose(1, 0, 2, 3)                         # (N, ng, nx, ny)


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

    :meth:`solve` operates on **typed source carriers** in the
    principled storage layout (see :ref:`theory-sn-index-convention`):
    :class:`~orpheus.sn.sources.IsotropicSource` shape ``(ng, nx, ny)``
    plus persistent :class:`~orpheus.sn.boundary_flux.BoundaryFlux` and
    optional :class:`~orpheus.sn.sources.PerOrdinateSource` shape
    ``(N, ng, nx, ny)``.  It returns a ``(angular_flux, scalar_flux)``
    tuple matching :func:`transport_sweep`'s contract.  The shape
    mismatch between the packed-vector ``apply`` (FD-matvec internal
    Fortran-flat layout, deferred to PR-INDEX-7) and the
    principled-storage ``solve`` reflects the historical layouts of
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
        Total cross-section, shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-3).
        Held as a separate attribute (not derived from ``sn_mesh``) because
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
            # PR-INDEX-3: ``sig_t`` layout is ``(ng, nx, ny)`` so the
            # group count lives at axis 0.
            nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.sig_t.shape[0]
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
        # PR-INDEX-3: ``sig_t`` shape is ``(ng, nx, ny)``.
        nx, ny, ng = sn_mesh.nx, sn_mesh.ny, self.sig_t.shape[0]
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
        iso_source: "IsotropicSource",
        boundary_flux: "BoundaryFlux | None" = None,
        aniso_source: "PerOrdinateSource | None" = None,
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

        Issue #197 PR-TYPED-4 — strict typed inputs only.  ``iso_source``
        MUST be an :class:`IsotropicSource`; ``aniso_source`` MUST be a
        :class:`PerOrdinateSource` or ``None``.  The legacy keyword-only
        ``Q_aniso`` alias is GONE.

        Parameters
        ----------
        iso_source : IsotropicSource
            Isotropic source density, shape ``(ng, nx, ny)`` (Issue
            #196 PR-INDEX-5 — principled).
        boundary_flux : BoundaryFlux or None
            Persistent boundary state (Issue #197 PR-TYPED-2 — the
            typed replacement for the legacy ``psi_bc: dict``).  If
            ``None``, a fresh zero-initialised
            :class:`~orpheus.sn.boundary_flux.BoundaryFlux` is supplied;
            the caller cannot then carry state between sweeps.
        aniso_source : PerOrdinateSource or None
            Per-ordinate anisotropic source, shape ``(N, ng, nx, ny)``,
            for P1+ scattering.  ``None`` for isotropic-only (P0).

        Returns
        -------
        tuple
            ``(angular_flux, scalar_flux)`` matching the
            :func:`transport_sweep` contract:

            * ``angular_flux`` shape ``(N, ng, nx, ny)``,
            * ``scalar_flux`` shape ``(ng, nx, ny)``.
        """
        from .sweep import transport_sweep
        if boundary_flux is None:
            boundary_flux = self.sn_mesh.zeros_boundary_flux()
        # ``self.sig_t`` and ``transport_sweep``'s ``sig_t`` parameter
        # are both in the principled ``(ng, nx, ny)`` layout under
        # PR-INDEX-3 — no bridge.
        return transport_sweep(
            iso_source, self.sig_t, self.sn_mesh, boundary_flux,
            aniso_source=aniso_source,
        )

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


# ═══════════════════════════════════════════════════════════════════════
# Phase G four-operator algebra — :class:`StreamingOperator` and
# :class:`CollisionOperator` leaves (Issue #196 Step 3+4.b.i, Resolution A).
# ═══════════════════════════════════════════════════════════════════════
#
# The Grand Report v3 §6 four-operator unification splits the historical
# bundled ``L = Ω·∇ + Σ_t`` into two independently composable leaves:
#
#     L  →  StreamingOperator        (streaming + angular redistribution)
#     C  →  CollisionOperator        (σ·, diagonal in pos/group/direction)
#
# so the within-group operator becomes pure algebra
#
#     A_wg = L + C - S.foldable_part()
#
# (an :class:`~orpheus.numerics.operator.OperatorSum`; no wrapper class).
#
# Resolution A — subtractive definition of L
# ------------------------------------------
#
# An earlier attempt (reverted commit ``ad37ca0``) implemented
# ``StreamingOperator.apply`` by calling the matvec primitive with
# ``Σ_t = 0`` and treating that as "pure streaming". This is
# mathematically WRONG for curvilinear: the matvec is **rational in
# σ_t** (not affine) through Hébert §3.9.4 Eq. 3.434's Carlson
# coupled-pole seed denominator ``Δr·σ_t + 2``. Setting ``σ_t = 0``
# degenerates the seed and produces a different operator from
# :math:`L = \Omega\cdot\nabla + {\rm redistribution}`. The empirical
# signature was a 3-13% relative error on random ψ in the composition
# test (see ``derivations/diagnostics/diag_LC_decomposition_sn.py``
# evidence memo).
#
# **Resolution A** defines L subtractively:
#
# .. math::
#
#     L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
#                                \sigma_t \odot \psi_{\rm packed}
#
# i.e. call the matvec primitive at the user's σ_t (not zero), then
# subtract the cell-collision term :math:`\sigma_t \odot \psi` at the
# packed-vector level. Combined with
# :math:`C.{\rm apply}(\psi) = \sigma_t \odot \psi`, this gives
#
# .. math::
#
#     (L + C).{\rm apply}(\psi) \;=\; M(\psi;\;\sigma_t)
#
# **bit-exact, by construction** — the decomposition is an algebraic
# identity, not an approximation. Verified at ``rel_residual = 0.0``
# across slab/sphere/cylinder × 3 random seeds in
# :file:`tests/sn/test_streaming_operator_decomposition.py`.
#
# This requires both L and C to carry ``σ_t`` at constructor time. The
# discrete L's σ_t parametrisation is a property of Hébert's M-M
# angular closure (the Carlson seed uses σ_t · φ_0 / W as the inward
# equivalent source at μ = −1); it is intrinsic to the discretisation,
# not a defect. The CONTINUOUS L (:math:`\Omega\cdot\nabla\psi +
# (1-\mu^2)/r\,\partial\psi/\partial\mu`) is σ_t-independent, but the
# DISCRETE curvilinear streaming-plus-redistribution operator inherits
# σ_t through its closure. Cartesian L has no Carlson seed; for
# Cartesian the subtraction removes σ_t exactly, so Cartesian L is
# σ_t-independent at the apply level.
#
# This substep (Issue #196 Phase G Step 3+4.b.i) lands the LEAVES; the
# fusion of ``(L + C - S.foldable_part()).solve(q)`` through the
# within-group sweep is substep 3+4.b.ii.


@dataclass
class StreamingOperator(LinearOperatorMixin):
    r"""Pure streaming + angular-redistribution operator :math:`L` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "L" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Carries the spatial
    streaming math plus the curvilinear angular redistribution — the
    cell-collision term lives in the separate :class:`CollisionOperator`
    leaf. The split lets the within-group operator be pure algebra
    (``L + C - S.foldable_part()``); no ``WithinGroupOperator`` wrapper.

    Resolution A — subtractive ``apply``
    ------------------------------------

    .. math::

        L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                                  \sigma_t \odot \psi_{\rm packed}

    The matvec primitive is called at the user's σ_t (constructor-stored),
    then the cell-collision term :math:`\sigma_t \odot \psi` is
    subtracted at the packed-vector level. Combined with
    :math:`C.{\rm apply}(\psi) = \sigma_t \odot \psi`, this gives
    :math:`(L + C).{\rm apply}(\psi) = M(\psi;\;\sigma_t)` bit-exact.

    Why σ_t is a CONSTRUCTOR parameter (Pattern 4)
    ----------------------------------------------

    The discrete curvilinear matvec is RATIONAL in σ_t through the
    Carlson coupled-pole seed (Hébert §3.9.4 Eqs. 3.432-3.435):

    .. math::

        \bar\phi_i \;=\; \frac{\Delta r_i\,\bar Q_i + 2\,\bar\phi_{i+1/2}}
                              {\Delta r_i\,\Sigma_t(i) + 2}

    where :math:`\bar Q_i = \Sigma_t(i)\,\phi_0(i) / W` at ℓ=0. The
    discrete L (with M-M angular closure) is therefore intrinsically
    σ_t-coupled — analogous to the DD coefficient :math:`\alpha_{DD}
    (\sigma_t\,\Delta x)` in characteristic-line methods, or the
    exponential characteristic :math:`\exp(-\sigma_t\,s)` in MoC/CP.
    Discrete streaming operators routinely carry σ_t through their
    closure choices.

    The CONTINUOUS L (:math:`\Omega\cdot\nabla\psi + (1-\mu^2)/r\,
    \partial\psi/\partial\mu`) is σ_t-independent. The DISCRETE L's
    σ_t parametrisation is a property of the discretisation. Pattern 4
    (illegal states unrepresentable) — ``StreamingOperator(sn_mesh)``
    without σ_t is not a meaningful object for curvilinear SN; the
    constructor refusing to construct without σ_t encodes that fact.

    Capability set
    --------------

    ``frozenset({CAP_APPLY})`` — pure streaming alone is **not
    invertible** (the streaming operator is rank-deficient without a
    collision term to make the within-group cell balance non-singular).
    The ``solve`` capability appears at the
    :class:`~orpheus.numerics.operator.OperatorSum` level: ``(L + C
    - S_foldable).solve(q)`` routes to the within-group sweep via the
    fusion hook substep 3+4.b.ii adds. No ``apply_transpose`` yet —
    the analytic reverse-direction sweep is Step 6 work.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry carrying quadrature, BCs, pole closure,
        and (for curvilinear) the precomputed connection coefficients.
        ``StreamingOperator`` reads ``sn_mesh.bc_*`` directly — no
        ``boundary`` constructor parameter.
    sigma_t : np.ndarray
        Total cross-section, shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-3).
        Carried at constructor time per Resolution A's subtractive
        definition.

    Notes
    -----
    This is Phase G Step 3+4.b.i — an **additive** introduction. The
    legacy :class:`SNStreamingOperator` (which bundles ``L + C`` under
    a single ``Σ_t`` constructor with three capabilities) remains in
    place; substep 3+4.c retires it after
    :class:`~orpheus.sn.solver.SNSolver` is wired to consume the leaf
    algebra.
    """

    sn_mesh: "SNMesh"
    sigma_t: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    # Lazy EquationMap cache (geometry-dispatched, same logic as
    # SNStreamingOperator._ensure_eq_map).
    _eq_map: EquationMap | None = field(default=None, init=False, repr=False)

    def _ensure_eq_map(self, ng: int) -> EquationMap:
        """Lazily build the geometry-appropriate :class:`EquationMap`."""
        if self._eq_map is None:
            nx, ny = self.sn_mesh.nx, self.sn_mesh.ny
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
        """Total packed scalar unknowns inferred from ``sigma_t.shape[0]``.

        PR-INDEX-3: ``sigma_t`` layout is ``(ng, nx, ny)`` — group at
        axis 0.
        """
        ng = int(self.sigma_t.shape[0])
        return self._ensure_eq_map(ng=ng).n_unknowns

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Subtractive forward action :math:`L\,\psi = M(\psi;\sigma_t)
        - \sigma_t \odot \psi`.

        Issue #197 PR-TYPED-6c Step 5 — routes through the single
        geometry-agnostic :func:`transport_operator_matvec_unified`
        for 1-D slab, sphere, and cylinder.  The three legacy
        per-geometry helpers retire fully at Step 7 once
        :class:`SNStreamingOperator` is also retired in favour of the
        ``(L + C)`` operator-algebra path.

        2-D Cartesian (``ny > 1``) remains on the legacy
        :func:`transport_operator_matvec` (FD via
        :func:`_compute_gradients`) — anti-diagonal wavefront sweeps
        need a different iteration structure than
        :meth:`SNMesh.dag_walk` and lie outside the PR-TYPED-6c scope.

        By Resolution A: :math:`(L + C).{\rm apply}(\psi) =
        M(\psi;\sigma_t)` bit-exact (the ``+`` is operator addition;
        the per-cell σ_t·ψ cancellation lives at the algebra layer,
        not inside this leaf).

        Parameters
        ----------
        psi : np.ndarray
            Packed solution vector, shape ``(n_unknowns,)``.

        Returns
        -------
        np.ndarray
            ``L·ψ`` as a packed vector, same shape as ``psi``.
        """
        sn_mesh = self.sn_mesh
        # PR-INDEX-3: ``sigma_t`` shape is ``(ng, nx, ny)``.
        ng = int(self.sigma_t.shape[0])
        eq_map = self._ensure_eq_map(ng=ng)
        if eq_map.n_unknowns != psi.size:
            raise ValueError(
                f"StreamingOperator.apply: packed psi size {psi.size} "
                f"does not match eq_map.n_unknowns {eq_map.n_unknowns} "
                f"(ng={ng})."
            )
        nx, ny = sn_mesh.nx, sn_mesh.ny
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)

        # 2-D Cartesian remains on the legacy FD path (anti-diagonal
        # wavefront, not dag_walk).  PR-TYPED-7 will absorb it.
        if curv is None and ny > 1:
            m_full = transport_operator_matvec(
                psi, eq_map, quad, self.sigma_t,
                nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
                bc_xmin=sn_mesh.bc_xmin,
                bc_xmax=sn_mesh.bc_xmax,
                bc_ymin=sn_mesh.bc_ymin,
                bc_ymax=sn_mesh.bc_ymax,
            )
        else:
            # 1-D slab / sphere / cylinder — unified matvec.
            if curv is None:
                psi_view = solution_to_angular_flux(
                    psi, eq_map, quad, nx, ny, ng,
                    bc_xmin=sn_mesh.bc_xmin,
                    bc_xmax=sn_mesh.bc_xmax,
                    bc_ymin=sn_mesh.bc_ymin,
                    bc_ymax=sn_mesh.bc_ymax,
                )
            else:
                psi_view = solution_to_angular_flux_spherical(
                    psi, eq_map, quad, nx, ng,
                )
            m_view = transport_operator_matvec_unified(
                psi_view, sn_mesh, self.sigma_t,
            )
            # Re-pack at the unknown slots — vectorised advanced indexing.
            # ``m_view[ord_arr, :, ix_arr, iy_arr]`` returns ``(n_eq, ng)``;
            # transpose to ``(ng, n_eq)`` then ravel(order='F') matches the
            # ``packed[g + ng·k]`` layout the legacy helpers emit.
            m_full = m_view[
                eq_map.ordinate, :, eq_map.ix, eq_map.iy,
            ].T.ravel(order="F")

        # Subtract σ_t ⊙ ψ at the packed-vector level (Resolution A).
        # Layout matches CollisionOperator's gather: ``(ng, nx, ny)``
        # advanced-indexed by ``(ix, iy)`` returns ``(ng, n_eq)``
        # directly under PR-INDEX-3; Fortran ravel gives (n_unknowns,)
        # packed vector aligned with `psi`.
        sigma_packed = self.sigma_t[
            :, eq_map.ix, eq_map.iy
        ].ravel(order='F')
        return m_full - sigma_packed * psi


@dataclass
class CollisionOperator(LinearOperatorMixin):
    r"""Pure collision operator :math:`C = \sigma\cdot` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "C" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Diagonal in position,
    group, and direction:

    .. math::

        (C\psi)_{n,i,j,g} \;=\; \sigma(i,j,g)\,\psi_{n,i,j,g}

    The supplied ``sigma`` is per-cell per-group only — the operator
    broadcasts identically across every ordinate ``n`` because the
    collision cross-section is direction-independent.

    Convention-agnostic σ
    ---------------------

    The same operator class supports the full total cross-section
    :math:`\sigma_t` AND the within-group removal cross-section
    :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`. The operator
    does NOT carry an interpretation flag — both quantities are
    ``(ng, nx, ny)`` arrays applied as :math:`\sigma\cdot\psi` (see
    :ref:`theory-sn-index-convention` for the canonical layout). The
    fusion hook (substep 3+4.b.ii) builds a lazy :math:`\sigma_r`
    :class:`CollisionOperator` from
    :meth:`~orpheus.sn.scattering.ScatteringOperator.foldable_sigma` at
    the moment it detects ``L + C - S.foldable_part()`` in an
    :class:`~orpheus.numerics.operator.OperatorSum`. (Pattern 4 —
    illegal states unrepresentable via the type, not a runtime flag
    the operator never inspects.)

    Capability set
    --------------

    ``frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`` —
    collision is the simplest sort of operator. ``solve(q) = q / σ``
    is element-wise division; ``apply_transpose == apply`` because the
    operator is self-adjoint (diagonal in every basis the SN method
    operates in). All three capabilities are analytic.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — used for the packed-vector
        :class:`EquationMap` dispatch (same as
        :class:`StreamingOperator`).
    sigma : np.ndarray
        Per-cell per-group cross-section, shape ``(ng, nx, ny)``
        (Issue #196 PR-INDEX-3). May be σ_t (full collision) or σ_r
        (removal — within-group self-scatter folded). The operator's
        action is identical.
    """

    sn_mesh: "SNMesh"
    sigma: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        )
    )

    # Lazy EquationMap cache — same dispatch as StreamingOperator.
    _eq_map: EquationMap | None = field(default=None, init=False, repr=False)

    def _ensure_eq_map(self, ng: int) -> EquationMap:
        """Lazily build the geometry-appropriate :class:`EquationMap`."""
        if self._eq_map is None:
            nx, ny = self.sn_mesh.nx, self.sn_mesh.ny
            quad = self.sn_mesh.quad
            curv = getattr(self.sn_mesh, "curvature", None)
            if curv == "spherical":
                self._eq_map = build_equation_map_spherical(nx, quad, ng)
            elif curv == "cylindrical":
                self._eq_map = build_equation_map_cylindrical(nx, quad, ng)
            else:
                self._eq_map = build_equation_map(nx, ny, quad, ng)
        return self._eq_map

    def _sigma_at_unknowns(self, eq_map: EquationMap, ng: int) -> np.ndarray:
        r"""Gather ``σ`` at each packed-unknown's ``(g, ix, iy)`` slot.

        Builds the packed-vector-shaped coefficient array
        ``(σ[g, ix[k], iy[k]])`` for ``k`` ordered as the eq_map's
        Fortran-order ``(g, k)`` flatten — exactly what
        ``psi.reshape(ng, n_eq, order='F')`` produces on the input
        side. Returns shape ``(n_unknowns,)``; element-wise multiply
        with the input ``psi``.

        PR-INDEX-3: ``sigma`` is principled ``(ng, nx, ny)`` — advanced
        indexing on ``(ix, iy)`` returns ``(ng, n_eq)`` directly, no
        transpose required.
        """
        sigma_per_eq = self.sigma[:, eq_map.ix, eq_map.iy]    # (ng, n_eq)
        return sigma_per_eq.ravel(order='F')                  # (n_unknowns,)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Forward action :math:`C\,\psi = \sigma\cdot\psi`."""
        ng = int(self.sigma.shape[0])
        eq_map = self._ensure_eq_map(ng)
        if eq_map.n_unknowns != psi.size:
            raise ValueError(
                f"CollisionOperator.apply: packed psi size {psi.size} "
                f"does not match eq_map.n_unknowns {eq_map.n_unknowns} "
                f"(ng={ng})."
            )
        sigma_packed = self._sigma_at_unknowns(eq_map, ng)
        return sigma_packed * psi

    def solve(self, q: np.ndarray) -> np.ndarray:
        r"""Inverse action :math:`C^{-1}\,q = q/\sigma` element-wise.

        Trivially invertible: collision is diagonal, so the inverse is
        per-slot reciprocal scaling. Returns NaN / Inf at slots where
        ``σ == 0`` per the IEEE-754 division contract — consumers
        constructing :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`
        must guarantee positivity (the operator does not check).
        """
        ng = int(self.sigma.shape[0])
        eq_map = self._ensure_eq_map(ng)
        if eq_map.n_unknowns != q.size:
            raise ValueError(
                f"CollisionOperator.solve: packed q size {q.size} "
                f"does not match eq_map.n_unknowns {eq_map.n_unknowns} "
                f"(ng={ng})."
            )
        sigma_packed = self._sigma_at_unknowns(eq_map, ng)
        return q / sigma_packed

    def apply_transpose(self, psi: np.ndarray) -> np.ndarray:
        r"""Adjoint action :math:`C^*\,\psi = \sigma\cdot\psi`.

        Equal to :meth:`apply` — collision is self-adjoint (diagonal
        operator). Returned bit-equal to ``apply(psi)``.
        """
        return self.apply(psi)
