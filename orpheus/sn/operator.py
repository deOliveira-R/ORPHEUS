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
    OperatorSum,
)

from orpheus.numerics.quadrature import Quadrature

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryOperator

    from .angular_flux import AngularFlux
    from .boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from .sources import IsotropicSource, PerOrdinateSource
    from .spatial.pole_angular_closure import PoleAngularClosure

__all__ = [
    "EquationMap",
    "build_equation_map",
    "build_equation_map_spherical",
    "build_equation_map_cylindrical",  # alias of _spherical
    "build_equation_map_with_traces",
    "solution_to_angular_flux",
    "solution_to_angular_flux_spherical",
    "solution_to_angular_flux_cylindrical",  # alias of _spherical
    "solution_to_angular_flux_with_traces",
    "pack_with_traces",
    "transport_operator_matvec",  # Cartesian 2-D FD; retires in PR-TYPED-7
    "transport_operator_matvec_unified",
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
    ``n_unknowns = (n_eq + n_face_outer + n_face_inner) * ng``.  Flat
    indexing under Fortran order: ``flat_idx = g + ng * k`` for group
    ``g`` and equation ``k``.  This is an **optimisation choice of
    the FD-matvec sparse-matrix CSR row order** — it is **orthogonal**
    to the user-visible field convention ``(N, ng, nx, ny)`` (which
    lives at the :func:`solution_to_angular_flux` decoder output and
    propagates end-to-end through every public API per
    :ref:`theory-sn-index-convention`).

    The packed vector is purely an internal contract between the
    ``build_equation_map_*`` family (CSR row-order producer) and the
    ``transport_operator_matvec_*`` family (CSR row-order consumer);
    it never crosses a public boundary.  The Fortran-flat ``(g, k)``
    layout was retained at PR-INDEX-7 §B.4 because flipping it
    multiplies sparse-matrix blast radius without changing user-visible
    behaviour.

    PR-TYPED-6.5 Phase 3b — face state in the packed vector (B1''
    architecture).  The matvec previously consumed the cell-CENTRE at
    ``i=nx-1`` as a proxy for the FACE flux at ``r=R`` when seeding
    the Carlson coupled-pole sweep (and the cell-centre at ``i=0`` for
    the slab inner BC).  That cell-centre-as-face proxy produced an
    ``O(h)`` discretisation gap between the matvec path and the
    sweep route (cylinder twin-path divergence ``rel ≈ 4e-3`` at
    ``nx=40``).  B1'' adds the boundary face flux as genuine unknowns
    in the packed vector: ``psi_face_outer`` carries the outward-
    ordinate outflow at the outer face (``r=R`` / ``x=L``);
    ``psi_face_inner`` carries the inward-ordinate outflow at the
    inner face (``x=0`` slab; future hollow / annulus inner edge).

    Attributes
    ----------
    n_eq : int
        Number of cell-centre angular unknowns per group.  Equals
        ``N * nx * ny`` for Cartesian geometries; equals
        ``N * nx * ny`` (no compression) under the B1'' layout for
        curvilinear too — boundary-inward cell-centres are GENUINE
        unknowns (no BC-resolved-slot skipping), because the BC trace
        acts on the FACE state, not the cell-centre.  Legacy curvilinear
        ``build_equation_map_spherical`` / ``_cylindrical`` still
        produce the compressed layout (``n_face_outer=0``,
        ``n_face_inner=0``) and remain the contract the legacy
        :class:`SNStreamingOperator` consumes; the face-aware factory
        :func:`build_equation_map_with_traces` produces the B1''
        layout for :class:`StreamingOperator`.
    n_unknowns : int
        Total scalar unknowns
        ``= (n_eq + n_face_outer + n_face_inner) * ng``.
    ordinate, ix, iy : ndarray of shape ``(n_eq,)``
        Per-equation indices into the ``(N, ng, nx, ny)`` angular flux.
        ``solution_to_angular_flux*`` scatters: ``fi[ordinate[k], g,
        ix[k], iy[k]] = sol[g + ng * k]``.
    n_face_outer, n_face_inner : int
        Number of face-flux unknowns per group at the outer / inner
        face.  ``n_face_outer = N_outward`` (count of outward ordinates
        at the outer face) when the outer face is a real boundary;
        ``n_face_inner = N_inward`` when the inner face is a real
        boundary (slab at ``x=0``; future hollow sphere / annulus).
        ``= 0`` for pole edges (sphere / cylinder solid at ``r=0``) and
        for the legacy cell-only layout.
    face_outer_ordinate, face_inner_ordinate : ndarray or None
        Per-face-equation ordinate indices.  ``face_outer_ordinate[k]``
        is the ordinate index at face slot ``k`` (consecutive starting
        from ``n_eq``); the corresponding spatial position is ``ix =
        nx - 1`` for the outer face / ``ix = 0`` for the inner face.
        ``None`` when the corresponding ``n_face_*`` is ``0``.
    """

    n_eq: int               # number of cell-centre angular unknowns (per group)
    n_unknowns: int         # (n_eq + n_face_outer + n_face_inner) * ng
    ordinate: np.ndarray    # (n_eq,) ordinate index for each cell unknown
    ix: np.ndarray          # (n_eq,) x-cell index
    iy: np.ndarray          # (n_eq,) y-cell index

    # ── PR-TYPED-6.5 Phase 3b: face-state slots in the packed vector ──
    # ``n_face_outer * ng`` scalars follow the cell-centre block at
    # packed offset ``n_eq * ng``; ``n_face_inner * ng`` scalars follow
    # at ``(n_eq + n_face_outer) * ng``.  ``face_*_ordinate`` arrays
    # have length ``n_face_*``.  All defaults preserve the legacy
    # cell-only contract.
    n_face_outer: int = 0
    n_face_inner: int = 0
    face_outer_ordinate: np.ndarray | None = None
    face_inner_ordinate: np.ndarray | None = None

    # R-1 Step 4: ``unknowns_at_cell_for_mask`` + ``_cell_ord_to_k``
    # retired (zero production callers — see R-1 Step 4 retirement
    # audit).  Originally an Issue #168 Phase C optimization for the
    # sweep-frame matvec; under the typed :class:`AngularFlux` /
    # :class:`BoundaryFlux` architecture the per-cell ordinate lookup
    # is direct array indexing on ``psi.values[:, :, i, j]`` — no
    # precomputed inverse table needed.


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


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical aliases — same compressed eq_map + decoder as spherical
# ═══════════════════════════════════════════════════════════════════════
#
# Cylindrical and spherical share the ``("left", "right")`` 1-D face
# structure and the "skip incoming-at-outer-boundary" eq_map convention.
# The aliases let consumers dispatch on ``sn_mesh.curvature`` without
# duplicating the function bodies.  PR-TYPED-6c Step 7 retired the
# per-geometry matvec helpers (``transport_operator_matvec_spherical``
# + ``transport_operator_matvec_cylindrical``); the eq_map + decoder
# aliases survive because they support the LEGACY packed-vector layout
# consumed by :class:`SNStreamingOperator` (the legacy bundle) and by
# ``solve_sn``'s Krylov inner.  Both retire together when ``solve_sn``
# migrates to the B1''-aware (L+C) algebra (PR-TYPED-7 follow-up).

build_equation_map_cylindrical = build_equation_map_spherical
solution_to_angular_flux_cylindrical = solution_to_angular_flux_spherical


# ═══════════════════════════════════════════════════════════════════════
# PR-TYPED-6.5 Phase 3b — B1'' face-aware packed-vector primitives
# ═══════════════════════════════════════════════════════════════════════
#
# The legacy ``build_equation_map_spherical`` / ``_cylindrical`` /
# ``solution_to_angular_flux_spherical`` family compresses inward-at-
# outer cell-centre slots out of the packed vector and synthesises
# them at decode time via an analytical extension (reflected-partner
# cell-centre).  The matvec then reads the cell-CENTRE at ``i = nx − 1``
# as a proxy for the FACE flux at ``r = R`` for the Carlson seed —
# producing an ``O(h)`` cell-centre-as-face proxy error.
#
# B1'' (face state in packed vector) replaces both pieces:
#
# * Every cell-centre slot is a GENUINE unknown (no compression at
#   the inward-at-outer cells; the cell-balance equation determines
#   them just like every other cell-centre).
# * The outward-ordinate outflow flux at the outer face is a NEW
#   unknown ``psi_face_outer`` (shape ``(N_outward, ng)``); the matvec
#   body reads this DIRECTLY for the Carlson seed (no cell-centre
#   proxy).
# * For slab + future hollow / annulus, the inward-ordinate outflow
#   flux at the inner face is a NEW unknown ``psi_face_inner`` (shape
#   ``(N_inward, ng)``); the matvec body reads this for the
#   ``bc_inner.apply(...)`` slab seed.
#
# The face-residual equations the matvec returns (m_face_outer,
# m_face_inner) drive GMRES to align the stored face state with the
# WDD-propagated face state at convergence.

def build_equation_map_with_traces(
    nx: int,
    quad: AngularQuadrature,
    ng: int,
    *,
    has_inner_bc: bool,
) -> EquationMap:
    r"""Build a face-aware :class:`EquationMap` for the B1'' packed layout.

    Cell-centre block: every ``(ordinate, cell)`` pair is an unknown
    (``n_eq = N · nx``) — no compression at the inward-at-outer cells
    that the legacy curvilinear maps skipped.  The Bc trace acts on
    the FACE state (``psi_face_outer``) directly, so the inward-at-
    outer cell-centres are determined by the inward sweep's WDD
    propagation just like every other inward cell-centre.

    Outer face block: ``n_face_outer = N_outward`` outward-ordinate
    slots (one per outward ordinate ``n`` with ``μ_x[n] > 0``)
    immediately following the cell-centre block.

    Inner face block (slab / hollow / annulus only):
    ``n_face_inner = N_inward`` inward-ordinate slots (one per inward
    ordinate ``n`` with ``μ_x[n] < 0``) following the outer face block.

    Parameters
    ----------
    nx :
        Number of spatial cells.
    quad :
        Angular quadrature (consumed for ``N`` and ``mu_x``).
    ng :
        Number of energy groups (sets ``n_unknowns = (n_eq +
        n_face_outer + n_face_inner) · ng``).
    has_inner_bc :
        ``True`` for slab and future hollow sphere / annulus (the
        inner edge is a real boundary requiring face state).
        ``False`` for sphere / cylinder solid at ``r = 0`` (the inner
        edge is a pole — Lewis-Miller §4.5 — and carries no face
        state; the cell-centre proxy at ``i = 0`` is the canonical
        pole IC).
    """
    mu_x = quad.mu_x
    N = quad.N

    # Cell-centre block: every (ordinate, cell) is an unknown.
    cell_ords = np.repeat(np.arange(N, dtype=int), nx)
    cell_ixs = np.tile(np.arange(nx, dtype=int), N)
    cell_iys = np.zeros_like(cell_ixs)
    n_eq = N * nx

    # Outer face: all outward ordinates (μ_x > 0).
    outward_mask = mu_x > 1e-15
    face_outer_ord = np.where(outward_mask)[0].astype(int)
    n_face_outer = int(face_outer_ord.size)

    # Inner face: all inward ordinates (μ_x < 0), only when the
    # inner edge is a real boundary (slab today).
    if has_inner_bc:
        inward_mask = mu_x < -1e-15
        face_inner_ord = np.where(inward_mask)[0].astype(int)
        n_face_inner = int(face_inner_ord.size)
    else:
        face_inner_ord = None
        n_face_inner = 0

    n_unknowns = (n_eq + n_face_outer + n_face_inner) * ng
    return EquationMap(
        n_eq=n_eq,
        n_unknowns=n_unknowns,
        ordinate=cell_ords,
        ix=cell_ixs,
        iy=cell_iys,
        n_face_outer=n_face_outer,
        n_face_inner=n_face_inner,
        face_outer_ordinate=face_outer_ord,
        face_inner_ordinate=face_inner_ord,
    )


def solution_to_angular_flux_with_traces(
    solution: np.ndarray,
    eq_map: EquationMap,
    nx: int,
    ng: int,
    *,
    N: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    r"""Decode the B1'' packed vector into cell-centre + face arrays.

    Inverse of :func:`pack_with_traces`; pure reshape, no analytical
    extensions or BC fills.

    Parameters
    ----------
    solution :
        Packed vector, shape ``(eq_map.n_unknowns,)``.
    eq_map :
        Face-aware :class:`EquationMap` built by
        :func:`build_equation_map_with_traces`.
    nx, ng :
        Spatial cell count + group count.
    N :
        Quadrature ordinate count.

    Returns
    -------
    psi_cell : ndarray, shape ``(N, ng, nx, 1)``
        Cell-centre angular flux scattered from the cell-centre
        block of ``solution`` per the eq_map's ``(ordinate, ix)``.
    psi_face_outer : ndarray, shape ``(n_face_outer, ng)`` or empty
        Outward-ordinate outflow flux at the outer face, decoded
        from the face-outer block immediately following the
        cell-centre block.  Empty array when ``n_face_outer == 0``.
    psi_face_inner : ndarray, shape ``(n_face_inner, ng)`` or None
        Inward-ordinate outflow flux at the inner face (slab /
        hollow / annulus).  ``None`` when ``n_face_inner == 0``
        (sphere / cylinder solid at the pole).
    """
    n_cell_scalars = eq_map.n_eq * ng
    cell_block = solution[:n_cell_scalars].reshape(ng, eq_map.n_eq, order="F")

    # Scatter cell block into (N, ng, nx, 1).
    psi_cell = np.zeros((N, ng, nx, 1))
    psi_cell[eq_map.ordinate, :, eq_map.ix, 0] = cell_block.T

    # Outer face block — shape (n_face_outer, ng).
    if eq_map.n_face_outer > 0:
        face_outer_start = n_cell_scalars
        face_outer_end = face_outer_start + eq_map.n_face_outer * ng
        psi_face_outer = solution[face_outer_start:face_outer_end].reshape(
            ng, eq_map.n_face_outer, order="F",
        ).T  # (n_face_outer, ng)
    else:
        psi_face_outer = np.zeros((0, ng))

    # Inner face block — shape (n_face_inner, ng) or None.
    if eq_map.n_face_inner > 0:
        face_inner_start = (
            n_cell_scalars + eq_map.n_face_outer * ng
        )
        face_inner_end = face_inner_start + eq_map.n_face_inner * ng
        psi_face_inner = solution[face_inner_start:face_inner_end].reshape(
            ng, eq_map.n_face_inner, order="F",
        ).T  # (n_face_inner, ng)
    else:
        psi_face_inner = None

    return psi_cell, psi_face_outer, psi_face_inner


def pack_with_traces(
    m_cell: np.ndarray,
    m_face_outer: np.ndarray,
    m_face_inner: np.ndarray | None,
    eq_map: EquationMap,
) -> np.ndarray:
    r"""Encode (cell, outer-face, inner-face) residuals into a B1''
    packed vector aligned with :func:`solution_to_angular_flux_with_traces`.

    Parameters
    ----------
    m_cell :
        Cell-centre residual, shape ``(N, ng, nx, 1)`` matching
        ``solution_to_angular_flux_with_traces``'s ``psi_cell``.
    m_face_outer :
        Outer-face residual, shape ``(n_face_outer, ng)``.
    m_face_inner :
        Inner-face residual, shape ``(n_face_inner, ng)``, or
        ``None`` for pole geometries.
    eq_map :
        Face-aware :class:`EquationMap`.

    Returns
    -------
    np.ndarray
        Packed vector, shape ``(eq_map.n_unknowns,)``.
    """
    parts: list[np.ndarray] = []

    # Cell-centre block — gather at eq_map's (ordinate, ix) slots,
    # transpose to (ng, n_eq), Fortran-ravel for the packed layout.
    cell_packed = m_cell[
        eq_map.ordinate, :, eq_map.ix, 0,
    ].T.ravel(order="F")
    parts.append(cell_packed)

    # Outer face block — (n_face_outer, ng) → (ng, n_face_outer) → F-ravel.
    if eq_map.n_face_outer > 0:
        parts.append(m_face_outer.T.ravel(order="F"))

    # Inner face block — same pattern.
    if eq_map.n_face_inner > 0:
        assert m_face_inner is not None, (
            "pack_with_traces: m_face_inner must be supplied when "
            "eq_map.n_face_inner > 0"
        )
        parts.append(m_face_inner.T.ravel(order="F"))

    return np.concatenate(parts) if len(parts) > 1 else parts[0]


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
    psi: "AngularFlux",
    sigma_t: np.ndarray,             # (ng, nx, ny)
    *,
    bc_outer: "BoundaryOperator | None" = None,
    pole_angular_closure: "PoleAngularClosure | None" = None,
) -> "AngularFlux":
    r"""Unified geometry-agnostic SN transport operator apply (path-forward).

    Computes :math:`M(\psi; \Sigma_t) = (L + C)\,\psi` over all
    1-D geometries (sphere, cylinder, slab) via the shared per-cell
    algebra in :func:`cell_balance_for_streaming` and (for curvilinear)
    :class:`MorelMontryAngularSweep`.

    R-1 Step 4 Step G0 (2026-05-22) — native-shape signature.  Inputs
    and outputs are typed :class:`AngularFlux` carrying the cell flux
    ``(N, ng, nx, ny)`` AND the boundary face state
    (:class:`BoundaryFlux`) at the outer + slab-inner faces in their
    natural ``(N, ng)`` layout.  The legacy packed-face-slot signature
    (``psi_face_outer: (n_face_outer, ng)`` + ``face_outer_ordinate``
    int slot map) retired with this rewrite — the inflow / outflow
    ordinate masks are derived from the quadrature direction signs
    (``quad.mu_x > 0`` for outflow at the outer face, etc.), NOT from
    a precomputed slot map.  Companion design notes in
    ``.claude/plans/r1_step4_g_convention_crosswalk.md`` Axis 3 + 6.

    The body is geometry-agnostic by data: the per-cell algebra reads
    cell geometry from ``sn_mesh.volumes`` + ``sn_mesh.areas``, and
    the per-ordinate angular closure from
    :class:`~orpheus.sn.spatial.pole_angular_closure._MMHalfGrid`.

    Issue #197 PR-TYPED-6c — replaced the three legacy geometry-keyed
    helpers (Pattern 2: single source of truth for the discretisation
    algebra).  R-1 Step 4 Step G0 — replaced the packed-face-slot
    I/O with native typed-AngularFlux I/O (Pattern 2 again: the eq_map
    slot maps disappear; ``BoundaryFlux`` IS the face-state carrier).

    Parameters
    ----------
    psi :
        Typed :class:`AngularFlux` carrying:

        * ``psi.values`` — angular flux in canonical layout
          ``(N, ng, nx, ny)``.
        * ``psi.boundary.xmax_face`` — outer-face flux ``(N, ng)``.
          Outflow positions (``quad.mu_x > 0``) carry the unknown
          face state the matvec accumulates a residual for; inflow
          positions are read by ``bc_outer.apply`` to produce the
          inward sweep's seed.
        * ``psi.boundary.xmin_face`` — inner-face flux ``(N, ng)``
          (slab only — sphere / cylinder have no inner face).
          Outflow positions (``quad.mu_x < 0``) carry the slab's
          inner-face unknown; inflow positions are read by
          ``bc_inner.apply``.

        The mesh is read as ``psi.mesh`` (the typed-flux carries it).
    sigma_t :
        Per-group per-cell total cross section, shape ``(ng, nx, ny)``.
        Issue #196 PR-INDEX-3 — group-leading.
    bc_outer :
        Outer-face boundary operator. ``None`` (default) reads
        ``psi.mesh.bc_right``.
    pole_angular_closure :
        Angular closure strategy. ``None`` (default) reads
        ``psi.mesh.pole_angular_closure``.

    Returns
    -------
    AngularFlux
        ``M(ψ; σ_t) = (L + C)·ψ`` in typed form:

        * ``.values`` — cell-centre matvec result ``(N, ng, nx, ny)``.
        * ``.boundary.xmax_face`` — outer-face residual ``(N, ng)``.
          At outflow positions (``quad.mu_x > 0``): ``WDD-propagated
          face − psi.boundary.xmax_face`` (driven to zero by GMRES at
          convergence).  At inflow positions: identically zero
          (no equation; the inflow value is determined by the BC,
          not an unknown).
        * ``.boundary.xmin_face`` — inner-face residual ``(N, ng)``
          for slab (zero at inflow positions; ``WDD-propagated −
          stored`` at outflow positions).  For sphere / cylinder,
          ``xmin_face`` is ``None`` per :class:`BoundaryFlux.zeros`
          dispatch on the mesh's curvature.

    Raises
    ------
    NotImplementedError
        2-D Cartesian (``ny > 1``) still routes through the legacy
        :func:`transport_operator_matvec` (FD via
        :func:`_compute_gradients`); anti-diagonal wavefront sweeps
        absorb in Phase A.
    """
    from .angular_flux import AngularFlux
    from .boundary_flux import BoundaryFlux
    from .spatial.cell_balance import cell_balance_for_streaming
    from .spatial.pole_angular_closure import MorelMontryAngularSweep

    if not isinstance(psi, AngularFlux):
        raise TypeError(
            f"transport_operator_matvec_unified expects an AngularFlux; "
            f"got {type(psi).__name__}.  R-1 Step 4 Step G0 retired the "
            f"packed-face-slot signature."
        )

    sn_mesh = psi.mesh
    psi_view = psi.values                                            # (N, ng, nx, ny)
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

    # ── Geometry data sourced from the natural owners ────────────────
    # Sphere, cylinder, slab all flow through the same body.  Slab and
    # sphere are the 1-level case (``level_indices = (arange(N),)``);
    # cylinder iterates over its μ-levels.  ``closure.cell_contribution``
    # reads the M-M coefficients (α, τ, ΔA/w, c_in, c_out) directly
    # from the closure object; ``IdentityAngularClosure`` returns zeros
    # so the per-cell algebra collapses to the slab form
    # ``2|μ|·1 + Σ_t·V`` (Step 2.5 principle — ``cell_balance_for_streaming``
    # is geometry-blind by data).  ``sn_mesh.areas`` returns the
    # geometry's face-area array (Cartesian: ones; cylinder: 2πr;
    # sphere: 4πr²) — Phase 1 canonicalised this on SNMesh directly.
    mu_x = quad.mu_x
    level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
    A = sn_mesh.areas                                                # (nx+1,)

    # ── Internal view: (ng, N, nx, ny) — group-leading for the
    # (ng, n_mask) per-cell algebra cell_balance_for_streaming consumes.
    # Zero-copy transpose; public interface stays canonical.
    psi_g_first = psi_view.transpose(1, 0, 2, 3)                     # (ng, N, nx, ny)
    out_g_first = np.zeros((ng, N, nx, ny))

    V = sn_mesh.volumes[:, 0]                                        # (nx,)
    sigma_t_gx = sigma_t[:, :, 0]                                    # (ng, nx)

    # ── Boundary face state — read directly from psi.boundary ────────
    # R-1 Step 4 Step G0 — boundary face state arrives in native
    # ``(N, ng)`` layout via :class:`BoundaryFlux`.  Slab has both
    # ``xmin_face`` and ``xmax_face``; sphere / cylinder have
    # ``xmax_face`` only (the inner edge is a pole, not a real face).
    boundary = psi.boundary
    face_outer = boundary.xmax_face                                  # (N, ng)
    face_inner = boundary.xmin_face                                  # (N, ng) for slab; None otherwise

    # ── Phase 1 spatial-upstream seed at the inner boundary ──────────
    # The predicate is structural, not curvature-keyed: ``bc_inner is
    # None`` means the inner edge is a pole (sphere / cylinder solid at
    # r=0 — Lewis-Miller §4.5: r=0 is the geometric pole, not a real
    # face, so the cell-centre proxy preserves the per-ordinate
    # flat-flux invariant). ``bc_inner is not None`` means there IS a
    # real inner boundary (slab at x=0, future hollow sphere /
    # annulus at r_inner), so apply ``bc_inner`` to the FACE flux
    # ``face_inner`` directly.
    if bc_inner is None:
        pole_face_seed = psi_view[:, :, 0, 0].copy()                 # (N, ng)
    elif face_inner is not None:
        pole_face_seed = bc_inner.apply(face_inner)                  # (N, ng)
    else:
        raise ValueError(
            "Slab geometry requires psi.boundary.xmin_face to be "
            "populated (R-1 Step 4 Step G0 retired the cell-centre "
            "proxy fallback inside the matvec; legacy SN consumers "
            "must build a BoundaryFlux carrying the cell-centre proxy "
            "as face state at their call site)."
        )

    # ── Pre-compute the angular-closure state ───────────────────────
    # B1'' (PR-TYPED-6.5 Phase 3b) — the Carlson coupled-pole seed
    # consumes the OUTER FACE flux directly, not the cell-CENTRE at
    # ``i = nx-1``.  The Hébert §3.9.4 Eqs. 3.432-3.435 recurrence
    # demands the FACE trace; the pre-Phase-3b cell-centre-as-face
    # proxy produced an ``O(h)`` discretisation gap that drove the
    # cylinder twin-path divergence (``rel ≈ 4e-3`` at ``nx = 40``).
    outer_inflow_estimate = bc_outer.apply(face_outer)               # (N, ng)

    # ``closure.precompute_psi_state`` returns a ``tuple[_MMHalfGrid, ...]``
    # for curvilinear (one element per μ-level; sphere has one) and
    # ``None`` for Cartesian (``IdentityAngularClosure``).  The matvec
    # body then reads the per-cell angular contribution via
    # ``closure.cell_contribution(psi_state, i, p, within_positions)``
    # — single primitive across both strategies (Pattern 2).
    psi_state = pole_angular_closure.precompute_psi_state(
        psi_view, sigma_t=sigma_t_gx,
        bc_outer_inflow_estimate=outer_inflow_estimate,
    )

    # Per-level M-M closure constants (c_in, c_out within-level) are
    # precomputed at strategy construction (see ``MorelMontryAngularSweep``
    # / ``IdentityAngularClosure``); read via ``closure.cell_contribution``.

    # ── Directional sweep primitive ─────────────────────────────────
    # The outward and inward sweeps are the SAME procedure applied with
    # opposite direction signs.  PR-TYPED-6.5 Phase 4 — extracted as a
    # nested closure so the two calls read as one operator applied with
    # ``direction_sign = +1`` / ``-1`` (Cardinal Rule 2 — single source
    # of truth for per-cell streaming + WDD face propagation; the
    # angular contribution comes from ``closure.cell_contribution``,
    # the spatial+collision balance from ``cell_balance_for_streaming``).
    #
    # Side-effect: writes to ``out_g_first[:, dir_ords, i, 0]`` at every
    # visited (cell, ordinate) and returns the accumulated WDD-
    # propagated face flux at the end of each level's sweep — the
    # outer-face accumulator for ``direction_sign=+1`` (outflow at
    # ``i=nx-1``) and the inner-face accumulator for
    # ``direction_sign=-1`` (outflow at ``i=0``).
    def _sweep_direction(
        direction_sign: int,
        psi_face_in_init: np.ndarray,                                # (N, ng)
    ) -> np.ndarray:                                                 # (ng, N)
        outflow_at_end = np.zeros((ng, N))
        for p, level_idx in enumerate(level_indices):
            level_idx_arr = np.asarray(level_idx)
            mu_level = mu_x[level_idx_arr]
            within_mask = (
                mu_level > +eps if direction_sign > 0
                else mu_level < -eps
            )
            if not np.any(within_mask):
                continue
            global_dir = level_idx_arr[within_mask]
            abs_mu = np.abs(mu_x[global_dir])
            within_positions = np.where(within_mask)[0]

            cell_indices = list(
                sn_mesh.dag_walk_cell_indices(
                    direction_sign=direction_sign, mu_level_idx=p,
                )
            )
            if not cell_indices:
                continue
            psi_face_in = psi_face_in_init[global_dir, :].T          # (ng, n_dir_p)

            for i in cell_indices:
                psi_cell = psi_g_first[:, global_dir, i, 0]
                # Strategy-supplied per-cell angular contribution.
                # M-M reads from its precomputed half-grid; Identity
                # returns zeros (Cartesian has no angular redistribution).
                angular_denom_term, angular_numer_upstream = (
                    pole_angular_closure.cell_contribution(
                        psi_state, i, p, within_positions,
                    )
                )
                # WDD downstream face area: ``A[i+1]`` for outward,
                # ``A[i]`` for inward — the face the streaming term
                # propagates ψ onto.
                A_downstream = A[i + 1] if direction_sign > 0 else A[i]
                denom, numer_upstream = cell_balance_for_streaming(
                    abs_mu=abs_mu,
                    A_downstream=A_downstream,
                    A_total=A[i] + A[i + 1],
                    total_xs=sigma_t_gx[:, i],
                    volume=V[i],
                    psi_face_in=psi_face_in,
                    angular_denom_term=angular_denom_term,
                    angular_numer_upstream=angular_numer_upstream,
                )
                m_full = (denom * psi_cell - numer_upstream) / V[i]
                out_g_first[:, global_dir, i, 0] = m_full

                # WDD face propagation: ψ_face_out = 2·ψ̄ − ψ_face_in.
                psi_face_in = 2.0 * psi_cell - psi_face_in
            # At loop end ``psi_face_in`` carries the WDD-propagated
            # face flux at the streaming-direction's terminal face for
            # each ordinate in this level.
            outflow_at_end[:, global_dir] = psi_face_in
        return outflow_at_end

    # Outward sweep (μ_x > 0): seeds at the inner edge with
    # ``pole_face_seed`` (cell-centre proxy at the pole, or
    # ``bc_inner.apply`` at slab/hollow inner BC); accumulates the WDD-
    # propagated outflow at the outer face (``i = nx-1``).
    outflow_at_boundary = _sweep_direction(
        direction_sign=+1, psi_face_in_init=pole_face_seed,
    )

    # ── BC trace at outer face ───────────────────────────────
    # Apply ``bc_outer`` to the WDD-propagated outward outflow at the
    # outer face (r=R for curvilinear, x=L for slab).  The single-pass
    # matvec does NOT chain inward outflow back into outward (that's
    # the SI sweep's job, not the matvec's).
    inflow_full = bc_outer.apply(outflow_at_boundary.T)              # (N, ng)

    # Inward sweep (μ_x < 0): seeds at the outer edge with the bc-
    # trace inflow; accumulates the WDD-propagated outflow at the inner
    # face (``i = 0``).  For slab + future hollow / annulus this
    # populates ``m_face_inner``.
    outflow_at_inner = _sweep_direction(
        direction_sign=-1, psi_face_in_init=inflow_full,
    )

    # ── Phase 3: degenerate ordinates (|μ_x| < eps), no radial flow ──
    # Sphere quadratures (GL) have no exact zeros — this branch is a
    # no-op for sphere by construction (degenerate_mask is empty).
    # Cylinder quadratures (LevelSymmetricSN) can carry degenerate
    # ordinates at certain orders; slab carries none either way.
    #
    # Degenerate ordinates have no streaming direction — A_downstream = 0
    # and ψ_face_in = 0 throughout.  The cell-balance reduces to
    # collision + angular-redistribution alone.  Per-cell angular
    # contribution still comes from ``closure.cell_contribution``, but
    # ordinates here are SCATTERED across levels (rather than within
    # a single level as in the directional sweep), so we accumulate
    # the contributions per (level, within-level) before the per-cell
    # balance call.
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
        abs_mu_deg = np.abs(mu_x[global_deg])
        zero_face = np.zeros((ng, n_deg))
        for i in range(nx):
            # Collect the per-ordinate angular contribution by routing
            # each degenerate ordinate through ``closure.cell_contribution``
            # at its own (level, within-level) position.  Reads the SAME
            # closure data as the directional sweep does (M-M's half-grid +
            # c_in/c_out; Identity's zeros) — single source of truth.
            angular_denom_term = np.empty(n_deg)
            angular_numer_upstream = np.empty((ng, n_deg))
            for col_idx in range(n_deg):
                denom_one, numer_one = pole_angular_closure.cell_contribution(
                    psi_state, i, deg_level[col_idx],
                    np.array([deg_within[col_idx]]),
                )
                angular_denom_term[col_idx] = denom_one[0]
                angular_numer_upstream[:, col_idx] = numer_one[:, 0]

            psi_cell = psi_g_first[:, global_deg, i, 0]              # (ng, n_deg)
            denom, numer_upstream = cell_balance_for_streaming(
                abs_mu=abs_mu_deg,
                A_downstream=0.0,
                A_total=A[i] + A[i + 1],
                total_xs=sigma_t_gx[:, i],
                volume=V[i],
                psi_face_in=zero_face,
                angular_denom_term=angular_denom_term,
                angular_numer_upstream=angular_numer_upstream,
            )
            m_full = (denom * psi_cell - numer_upstream) / V[i]
            out_g_first[:, global_deg, i, 0] = m_full

    m_cell = out_g_first.transpose(1, 0, 2, 3)                       # (N, ng, nx, ny)

    # ── Face residuals (R-1 Step 4 Step G0 native-shape) ─────────────
    # ``outflow_at_boundary`` holds the WDD-propagated outflow at the
    # outer face (``i = nx-1``) accumulated by the outward sweep;
    # ``outflow_at_inner`` holds the WDD-propagated outflow at the inner
    # face (``i = 0``) accumulated by the inward sweep.  Each is shape
    # ``(ng, N)`` with non-zero columns only at the within-direction
    # ordinates (``mu_x > 0`` for outward; ``mu_x < 0`` for inward).
    #
    # The face residual is ``WDD-propagated − stored`` at outflow
    # positions (GMRES drives it to zero, aligning the stored face with
    # the WDD-propagated face at convergence).  Inflow positions carry
    # no equation (the inflow value comes from the BC, not an unknown);
    # their residual stays at zero (the BoundaryFlux default).
    # Outflow / inflow masks are derived from the quadrature direction
    # sign — NO precomputed ``face_outer_ordinate`` slot map (R-1 Step 4
    # Step G0 — removed legacy eq_map slot dispatch).
    #
    # Boundary arrays are sized to the INPUT flux's ng (``psi_view.shape[1]``)
    # NOT ``sn_mesh.ng`` — bare-ndarray callers may carry a sig_t whose
    # ng differs from the mesh's materials ng (the dimensional sin
    # acknowledged in #205).  Typed callers have ng-consistency by
    # construction; the input-ng allocation is correct for both.
    m_boundary = BoundaryFlux(mesh=sn_mesh)
    m_boundary.xmax_face = np.zeros((N, ng))
    if face_inner is not None:
        m_boundary.xmin_face = np.zeros((N, ng))

    outer_outflow_mask = mu_x > +eps                                 # (N,) bool
    if np.any(outer_outflow_mask):
        m_boundary.xmax_face[outer_outflow_mask, :] = (
            outflow_at_boundary[:, outer_outflow_mask].T
            - face_outer[outer_outflow_mask, :]
        )

    if face_inner is not None:
        inner_outflow_mask = mu_x < -eps                             # (N,) bool
        if np.any(inner_outflow_mask):
            m_boundary.xmin_face[inner_outflow_mask, :] = (
                outflow_at_inner[:, inner_outflow_mask].T
                - face_inner[inner_outflow_mask, :]
            )

    return AngularFlux(m_cell, sn_mesh, boundary=m_boundary)


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

        Issue #197 PR-TYPED-6c Step 7 — 1-D slab / sphere / cylinder
        all route through the single geometry-agnostic
        :func:`transport_operator_matvec_unified` (no face state — the
        LEGACY proxy fallback preserves the SN bundle's pre-B1''
        packed format and the cell-centre-as-face Carlson seed used by
        the production source-iteration paths).  2-D Cartesian
        (``ny > 1``) keeps the legacy :func:`transport_operator_matvec`
        FD path — anti-diagonal wavefront sweeps need a different
        iteration structure than :meth:`SNMesh.dag_walk` and absorb
        into the algebra in PR-TYPED-7.

        Side-effect of the rewire on cylinder: the legacy
        :func:`transport_operator_matvec_cylindrical` carried ERR-049
        (per-ordinate routing bug surfaced during PR-TYPED-6c Step 3);
        the unified body's per-cell algebra is structurally immune.
        SNStreamingOperator-via-unified therefore agrees with the
        unified-via-StreamingOperator path on cylinder up to the B1''
        face-state difference (the SN path still uses the cell-centre
        proxy at the Carlson seed — B1'' lives on the new
        :class:`StreamingOperator` leaf).  Slab's FD-vs-WDD difference
        is structural (FD 1st-order ≠ WDD 2nd-order) and persists as
        an order-of-accuracy delta documented in
        ``test_unified_matvec_slab.py``.

        Parameters
        ----------
        psi : np.ndarray
            Packed solution vector, shape ``(n_unknowns,)``.

        Returns
        -------
        np.ndarray
            ``L·ψ`` as a packed vector, same shape as ``psi``.
        """
        eq_map = self._ensure_eq_map()
        sn_mesh = self.sn_mesh
        # PR-INDEX-3: ``sig_t`` shape is ``(ng, nx, ny)``.
        nx, ny, ng = sn_mesh.nx, sn_mesh.ny, self.sig_t.shape[0]
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)

        # 2-D Cartesian remains on the legacy FD path (anti-diagonal
        # wavefront, not dag_walk).  PR-TYPED-7 will absorb it.
        if curv is None and ny > 1:
            return transport_operator_matvec(
                psi, eq_map, quad, self.sig_t,
                nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
                bc_xmin=sn_mesh.bc_xmin,
                bc_xmax=sn_mesh.bc_xmax,
                bc_ymin=sn_mesh.bc_ymin,
                bc_ymax=sn_mesh.bc_ymax,
            )

        # 1-D slab / sphere / cylinder — unified matvec via the LEGACY
        # cell-only packed format.  The decoder reconstructs ψ_view
        # with the analytical extension at inward-at-outer cells
        # (curvilinear); the unified matvec body uses the cell-centre
        # proxy at the Carlson seed (the bug B1'' fixes — but the
        # SNStreamingOperator path STAYS on the proxy through Step 7
        # to preserve the legacy bundle's packed-format contract for
        # ``solve_sn`` consumers).
        if curv in ("spherical", "cylindrical"):
            psi_view = solution_to_angular_flux_spherical(
                psi, eq_map, quad, nx, ng,
            )
        else:
            psi_view = solution_to_angular_flux(
                psi, eq_map, quad, nx, ny, ng,
                bc_xmin=sn_mesh.bc_xmin,
                bc_xmax=sn_mesh.bc_xmax,
                bc_ymin=sn_mesh.bc_ymin,
                bc_ymax=sn_mesh.bc_ymax,
            )
        # R-1 Step 4 Step G0 — wrap the legacy cell-only packed view into
        # a typed AngularFlux for the path-forward matvec.  The legacy
        # SN bundle's pre-B1'' contract is "cell-centre as face proxy";
        # build a synthetic BoundaryFlux carrying that proxy so the
        # path-forward matvec can do the BC-trace step.  This
        # synthesise-then-call shim retires with SNStreamingOperator at
        # G3f — the legacy bundle's contract is the only thing keeping
        # the proxy alive.
        from .angular_flux import AngularFlux
        from .boundary_flux import BoundaryFlux
        boundary_legacy = BoundaryFlux.zeros(sn_mesh)
        boundary_legacy.xmax_face[:, :] = psi_view[:, :, -1, 0]  # (N, ng)
        if boundary_legacy.xmin_face is not None:
            boundary_legacy.xmin_face[:, :] = psi_view[:, :, 0, 0]
        psi_typed = AngularFlux(psi_view, sn_mesh, boundary=boundary_legacy)
        result = transport_operator_matvec_unified(psi_typed, self.sig_t)
        m_view = result.values
        # Re-pack at the legacy eq_map's unknown slots.
        return m_view[
            eq_map.ordinate, :, eq_map.ix, eq_map.iy,
        ].T.ravel(order="F")

    # ── solve: L^{-1}·q via the unified sweep ─────────────────────────

    def solve(
        self,
        source: "PerOrdinateSource",
        boundary_flux: "BoundaryFlux | None" = None,
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

        R-1 Step 4 A1 — single per-ordinate source parameter.  The
        legacy ``(iso_source, aniso_source)`` parameter pair is GONE.
        External iso scalar sources project to per-ordinate at
        construction via
        :meth:`~orpheus.sn.sources.PerOrdinateSource.from_isotropic`;
        scattering- and fission-generated sources are already per-
        ordinate by producer contract
        (:meth:`~orpheus.sn.scattering.ScatteringOperator.apply`,
        :meth:`~orpheus.sn.fission.FissionOperator.apply`).

        Parameters
        ----------
        source : PerOrdinateSource
            Per-ordinate volumetric source, shape ``(N, ng, nx, ny)``.
            Per-ordinate magnitude (the producer has applied any
            required ``/W`` projection).
        boundary_flux : BoundaryFlux or None
            Persistent boundary state (Issue #197 PR-TYPED-2 — the
            typed replacement for the legacy ``psi_bc: dict``).  If
            ``None``, a fresh zero-initialised
            :class:`~orpheus.sn.boundary_flux.BoundaryFlux` is supplied;
            the caller cannot then carry state between sweeps.

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
            source, self.sig_t, self.sn_mesh, boundary_flux,
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
        """Lazily build the geometry-appropriate :class:`EquationMap`.

        PR-TYPED-6.5 Phase 3b — 1-D paths build a face-aware eq_map
        via :func:`build_equation_map_with_traces` so the packed vector
        carries ``ψ_face_outer`` (curvilinear + slab) and
        ``ψ_face_inner`` (slab + future hollow / annulus) alongside
        the cell-centre block.  2-D Cartesian (``ny > 1``) retains the
        legacy :func:`build_equation_map` since it routes through the
        unchanged 2-D FD path.
        """
        if self._eq_map is None:
            nx, ny = self.sn_mesh.nx, self.sn_mesh.ny
            quad = self.sn_mesh.quad
            curv = getattr(self.sn_mesh, "curvature", None)
            if curv is None and ny > 1:
                # 2-D Cartesian — legacy FD path
                self._eq_map = build_equation_map(nx, ny, quad, ng)
            else:
                # 1-D — B1'' face-aware (PR-TYPED-6.5 Phase 3b).
                # ``has_inner_bc`` distinguishes slab (real inner BC
                # at ``x=0``) from sphere / cylinder solid (pole at
                # ``r=0``, no face state).
                has_inner_bc = (curv is None)
                self._eq_map = build_equation_map_with_traces(
                    nx, quad, ng, has_inner_bc=has_inner_bc,
                )
        return self._eq_map

    @property
    def n_unknowns(self) -> int:
        """Total packed scalar unknowns inferred from ``sigma_t.shape[0]``.

        PR-INDEX-3: ``sigma_t`` layout is ``(ng, nx, ny)`` — group at
        axis 0.
        """
        ng = int(self.sigma_t.shape[0])
        return self._ensure_eq_map(ng=ng).n_unknowns

    def apply(
        self, psi: "np.ndarray | AngularFlux",
    ) -> "np.ndarray | AngularFlux":
        r"""Subtractive forward action :math:`L\,\psi = M(\psi;\sigma_t)
        - \sigma_t \odot \psi`.

        Issue #197 PR-TYPED-6c Step 5 — routes through the single
        geometry-agnostic :func:`transport_operator_matvec_unified`
        for 1-D slab, sphere, and cylinder.  The three legacy
        per-geometry helpers retire fully at Step 7 once
        :class:`SNStreamingOperator` is also retired in favour of the
        ``(L + C)`` operator-algebra path.

        PR-TYPED-6.5 Phase 3b — 1-D paths consume the B1'' face-aware
        packed vector (cell-centre block + outer-face block + inner-
        face block for slab) so the matvec body reads the actual face
        flux from the packed state instead of synthesising it from the
        cell-centre proxy at ``i = nx-1`` (curvilinear) or ``i = 0``
        (slab).  The cylinder twin-path divergence (rel ≈ 4e-3 at
        ``nx = 40``) gets fixed at its architectural source.

        R-1 Step 3d — typed :class:`AngularFlux` overload.  The Carlson
        seed and inner-BC seed are pulled from ``psi.boundary``
        directly (no flat-vector round-trip); the matvec returns the
        cell action and the face residual; the result is a typed
        :class:`AngularFlux` whose ``.boundary`` carries the B1'' face
        residual at the outer face (and, for slab, the inner face).
        This is the ONLY operator that emits a non-zero face residual;
        :class:`CollisionOperator`, :class:`~orpheus.sn.scattering.ScatteringOperator`,
        and :class:`~orpheus.sn.fission.FissionOperator` all leave the
        output ``.boundary`` at the auto-allocated zero.

        2-D Cartesian (``ny > 1``) remains on the legacy
        :func:`transport_operator_matvec` (FD via
        :func:`_compute_gradients`) — anti-diagonal wavefront sweeps
        need a different iteration structure than
        :meth:`SNMesh.dag_walk` and lie outside the PR-TYPED-6c scope.
        The R-1 typed overload routes 2-D through a flat round-trip
        via :meth:`AngularFlux.to_flat_with_traces` /
        :meth:`AngularFlux.from_flat_with_traces` for parity with 1-D
        (the underlying compute remains the legacy FD; Phase A absorbs
        the 2-D path).

        By Resolution A: :math:`(L + C).{\rm apply}(\psi) =
        M(\psi;\sigma_t)` bit-exact (the ``+`` is operator addition;
        the per-cell σ_t·ψ cancellation lives at the algebra layer,
        not inside this leaf).

        Parameters
        ----------
        psi : np.ndarray or AngularFlux
            * Bare ``np.ndarray`` — packed 1-D vector,
              shape ``(n_unknowns,)``.
            * :class:`~orpheus.sn.angular_flux.AngularFlux` — typed
              flux with embedded :class:`BoundaryFlux` carrying the
              face state at outer (and slab-inner) face.

        Returns
        -------
        np.ndarray or AngularFlux
            ``L·ψ`` in the matching type.  When :class:`AngularFlux`,
            ``.boundary.xmax_face`` carries the outer-face residual
            and (for slab) ``.boundary.xmin_face`` carries the
            inner-face residual.
        """
        from .angular_flux import AngularFlux
        from .boundary_flux import BoundaryFlux
        if isinstance(psi, AngularFlux):
            return self._apply_typed(psi)
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
            # Subtract σ_t ⊙ ψ at the legacy packed slots.
            sigma_packed = self.sigma_t[
                :, eq_map.ix, eq_map.iy
            ].ravel(order='F')
            return m_full - sigma_packed * psi

        # 1-D slab / sphere / cylinder — B1'' face-aware unified matvec.
        # R-1 Step 4 Step G0 — the path-forward matvec consumes a typed
        # AngularFlux + native BoundaryFlux.  The bare-ndarray packed
        # path stays as an adapter that decodes packed → typed, calls
        # the path-forward matvec, encodes typed → packed.  This shim
        # retires at G3d when the _with_traces family + EquationMap
        # retire (the legacy bare-ndarray caller chain dies with them).
        from .angular_flux import AngularFlux
        from .boundary_flux import BoundaryFlux
        psi_cell, psi_face_outer, psi_face_inner = (
            solution_to_angular_flux_with_traces(
                psi, eq_map, nx, ng, N=quad.N,
            )
        )
        # Scatter the packed face sub-arrays into native (N, ng) BoundaryFlux
        # slots so the path-forward matvec can read them directly.  Slots
        # not in the eq_map's face-ordinate maps stay at zero (the
        # correct "no equation here" value).  Allocate by the INPUT ng
        # (``sigma_t.shape[0]``) rather than ``sn_mesh.ng`` — bare-ndarray
        # callers may carry a sig_t whose ng differs from the mesh's
        # materials ng (the dimensional sin per #205).
        boundary_in = BoundaryFlux(mesh=sn_mesh)
        boundary_in.xmax_face = np.zeros((quad.N, ng))
        if curv is None:  # slab has a real inner face
            boundary_in.xmin_face = np.zeros((quad.N, ng))
        if eq_map.n_face_outer > 0:
            boundary_in.xmax_face[eq_map.face_outer_ordinate, :] = (
                psi_face_outer
            )
        if eq_map.n_face_inner > 0 and boundary_in.xmin_face is not None:
            boundary_in.xmin_face[eq_map.face_inner_ordinate, :] = (
                psi_face_inner
            )
        psi_typed = AngularFlux(psi_cell, sn_mesh, boundary=boundary_in)
        result = transport_operator_matvec_unified(psi_typed, self.sigma_t)
        m_cell = result.values
        # Gather face residuals back to packed sub-arrays at the eq_map
        # slot positions (the legacy packed layout's contract).
        m_face_outer = (
            result.boundary.xmax_face[eq_map.face_outer_ordinate, :]
            if eq_map.n_face_outer > 0
            else None
        )
        m_face_inner = (
            result.boundary.xmin_face[eq_map.face_inner_ordinate, :]
            if eq_map.n_face_inner > 0 and result.boundary.xmin_face is not None
            else None
        )
        # Pack matvec output into the B1'' packed layout.
        m_full = pack_with_traces(m_cell, m_face_outer, m_face_inner, eq_map)

        # Subtract σ_t ⊙ ψ at the CELL-CENTRE slots ONLY.  Face slots
        # carry no volumetric collision (the cell-balance σ·ψ term is
        # a CELL quantity; the face residual ``m_face_*`` is a TRACE
        # equation).  Resolution A's subtractive form thus zeros out
        # at face slots: ``L·ψ_face = M·ψ_face − 0 = m_face_*``.
        n_cell_scalars = eq_map.n_eq * ng
        sigma_packed_cell = self.sigma_t[
            :, eq_map.ix, eq_map.iy
        ].ravel(order='F')
        m_full[:n_cell_scalars] = (
            m_full[:n_cell_scalars] - sigma_packed_cell * psi[:n_cell_scalars]
        )
        return m_full

    def _apply_typed(self, psi: "AngularFlux") -> "AngularFlux":
        r"""Typed-AngularFlux body of :meth:`apply` (R-1 Step 3d).

        Routes natively through :func:`transport_operator_matvec_unified`
        for 1-D (no flat round-trip), and routes 2-D Cartesian through
        a flat round-trip via :meth:`AngularFlux.to_flat_with_traces`
        for parity (the underlying compute kernel for 2-D is the legacy
        :func:`transport_operator_matvec`; Phase A absorbs it).

        Returns an :class:`AngularFlux` whose ``.boundary`` carries the
        face residual (the only operator in the L/C/S/F set that emits
        a non-zero face residual).
        """
        from .angular_flux import AngularFlux
        from .boundary_flux import BoundaryFlux
        sn_mesh = self.sn_mesh
        if sn_mesh is not psi.mesh:
            raise ValueError(
                "StreamingOperator.apply(AngularFlux): operator and flux "
                "must share the SAME SNMesh instance (mesh-identity "
                "invariant)."
            )
        ng = int(self.sigma_t.shape[0])
        nx, ny = sn_mesh.nx, sn_mesh.ny
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)

        # 2-D Cartesian — flat round-trip onto the legacy compute.
        if curv is None and ny > 1:
            flat_in = psi.to_flat_with_traces()
            flat_out = self.apply(flat_in)
            return AngularFlux.from_flat_with_traces(flat_out, sn_mesh)

        # 1-D slab / sphere / cylinder — typed-native compute.  R-1 Step 4
        # Step G0 — the path-forward matvec consumes ``psi`` (AngularFlux)
        # natively.  No eq_map gather/scatter; ``psi.boundary.xmax_face``
        # carries the face state in (N, ng) layout, inflow vs outflow
        # masks derived from the quadrature direction sign inside the
        # matvec body.  ``_ensure_eq_map`` retires alongside this gather
        # block at G3b — the typed leaf operator carries no eq_map state.
        result = transport_operator_matvec_unified(psi, self.sigma_t)

        # Subtract σ_t ⊙ ψ at the CELL-CENTRE only — face slots carry
        # no volumetric collision (the cell-balance σ·ψ term is a CELL
        # quantity; the face residual is a TRACE equation).
        cell_values = result.values - self.sigma_t[None, :, :, :] * psi.values

        # Face residual is the matvec's boundary output as-is — the
        # path-forward matvec writes the residual at outflow positions
        # only (inflow stays zero per ``BoundaryFlux.zeros``).
        return AngularFlux(cell_values, sn_mesh, boundary=result.boundary)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`L + X`.

        When ``X`` is a :class:`CollisionOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        carrying the algebraic identity :math:`(L + C)^{-1} \approx
        \text{WDD sweep}`.  Otherwise falls through to the generic
        :class:`OperatorSum` via the mixin.
        """
        if isinstance(other, CollisionOperator):
            return InvertibleOperator(self, other)
        return super().__add__(other)


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
        """Lazily build the geometry-appropriate :class:`EquationMap`.

        Matches :meth:`StreamingOperator._ensure_eq_map` so the (L+C)
        algebra composes on the same packed format.  1-D uses the
        B1'' face-aware layout (PR-TYPED-6.5 Phase 3b); 2-D Cartesian
        retains the legacy FD layout.
        """
        if self._eq_map is None:
            nx, ny = self.sn_mesh.nx, self.sn_mesh.ny
            quad = self.sn_mesh.quad
            curv = getattr(self.sn_mesh, "curvature", None)
            if curv is None and ny > 1:
                self._eq_map = build_equation_map(nx, ny, quad, ng)
            else:
                has_inner_bc = (curv is None)
                self._eq_map = build_equation_map_with_traces(
                    nx, quad, ng, has_inner_bc=has_inner_bc,
                )
        return self._eq_map

    def _sigma_at_unknowns(self, eq_map: EquationMap, ng: int) -> np.ndarray:
        r"""Gather ``σ`` at each packed CELL-CENTRE unknown's
        ``(g, ix, iy)`` slot.

        Returns shape ``(n_eq · ng,)`` — covers the cell-centre block
        of the packed vector only.  Face slots (B1'') carry no
        volumetric collision, so this gather does NOT include them.
        Callers element-wise multiply the cell-centre block of
        ``psi`` by this; face slots are handled separately (zero
        contribution to ``apply``).

        PR-INDEX-3: ``sigma`` is principled ``(ng, nx, ny)`` — advanced
        indexing on ``(ix, iy)`` returns ``(ng, n_eq)`` directly, no
        transpose required.
        """
        sigma_per_eq = self.sigma[:, eq_map.ix, eq_map.iy]    # (ng, n_eq)
        return sigma_per_eq.ravel(order='F')                  # (n_eq · ng,)

    def apply(
        self, psi: "np.ndarray | AngularFlux",
    ) -> "np.ndarray | AngularFlux":
        r"""Forward action :math:`C\,\psi = \sigma\cdot\psi`.

        PR-TYPED-6.5 Phase 3b — face slots receive zero contribution
        (no volumetric collision on the boundary face trace); the
        cell-centre block carries the elementwise ``σ ⊙ ψ_cell``.

        R-1 Step 3c — typed :class:`AngularFlux` overload.  The diagonal
        action :math:`\sigma\cdot\psi` is per-cell per-group, broadcast
        across every ordinate.  Result's ``.boundary`` is the
        auto-allocated zero :class:`BoundaryFlux` — collision has no
        face-trace contribution.

        Parameters
        ----------
        psi : np.ndarray or AngularFlux
            * Bare ``np.ndarray`` — packed 1-D vector matching
              :attr:`_eq_map.n_unknowns`.  Returns bare ``np.ndarray``.
            * :class:`~orpheus.sn.angular_flux.AngularFlux` shape
              ``(N, ng, nx, ny)``.  Returns :class:`AngularFlux` with
              ``.values = σ ⊙ ψ.values`` and zero ``.boundary``.
        """
        from .angular_flux import AngularFlux
        if isinstance(psi, AngularFlux):
            # Typed branch — broadcast σ across ordinate axis.
            # σ is (ng, nx, ny); ψ.values is (N, ng, nx, ny);
            # the product broadcasts σ to every ordinate (collision
            # is direction-independent).  The output's ``.boundary``
            # is auto-allocated to zeros (R-1 Step 2), which is the
            # correct algebraic value — C has no face contribution.
            return AngularFlux(
                self.sigma[None, :, :, :] * psi.values, psi.mesh,
            )
        ng = int(self.sigma.shape[0])
        eq_map = self._ensure_eq_map(ng)
        if eq_map.n_unknowns != psi.size:
            raise ValueError(
                f"CollisionOperator.apply: packed psi size {psi.size} "
                f"does not match eq_map.n_unknowns {eq_map.n_unknowns} "
                f"(ng={ng})."
            )
        n_cell_scalars = eq_map.n_eq * ng
        out = np.zeros_like(psi)
        sigma_packed = self._sigma_at_unknowns(eq_map, ng)
        out[:n_cell_scalars] = sigma_packed * psi[:n_cell_scalars]
        return out

    def solve(
        self, q: "np.ndarray | AngularFlux",
    ) -> "np.ndarray | AngularFlux":
        r"""Inverse action :math:`C^{-1}\,q = q/\sigma` element-wise on
        the cell-centre block.

        Trivially invertible on the cell-centre block: collision is
        diagonal, so the inverse is per-slot reciprocal scaling.
        Returns NaN / Inf at slots where ``σ == 0`` per the IEEE-754
        division contract — consumers constructing
        :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}` must
        guarantee positivity (the operator does not check).

        PR-TYPED-6.5 Phase 3b — the face block has no σ (rank-deficient
        for pure ``C`` on the face).  Face slots of ``q`` pass through
        unchanged (identity-on-face) — this is the formal pseudoinverse
        for the rank-deficient face block and the right behaviour when
        ``(L + C).solve`` consumes a packed ``q`` with zero face entries
        (the standard within-group source).  Callers that genuinely
        need ``C^{-1}`` on a non-zero face block (an unusual scenario)
        must handle the face slots themselves.

        R-1 Step F — typed :class:`AngularFlux` overload.  The
        inverse-collision action is per-cell per-group, broadcast
        across every ordinate.  Result's ``.boundary`` is the
        auto-allocated zero :class:`BoundaryFlux` (collision is rank-
        deficient on the face block; pseudoinverse leaves it zero).
        """
        from .angular_flux import AngularFlux
        if isinstance(q, AngularFlux):
            # Typed branch — broadcast 1/σ across the ordinate axis.
            # σ is (ng, nx, ny); q.values is (N, ng, nx, ny).  The
            # result's ``.boundary`` is auto-allocated to zeros — the
            # pseudoinverse on the rank-deficient face block.
            return AngularFlux(
                q.values / self.sigma[None, :, :, :], q.mesh,
            )
        ng = int(self.sigma.shape[0])
        eq_map = self._ensure_eq_map(ng)
        if eq_map.n_unknowns != q.size:
            raise ValueError(
                f"CollisionOperator.solve: packed q size {q.size} "
                f"does not match eq_map.n_unknowns {eq_map.n_unknowns} "
                f"(ng={ng})."
            )
        n_cell_scalars = eq_map.n_eq * ng
        out = q.copy()
        sigma_packed = self._sigma_at_unknowns(eq_map, ng)
        out[:n_cell_scalars] = q[:n_cell_scalars] / sigma_packed
        # Face slots (out[n_cell_scalars:]) pass through unchanged.
        return out

    def apply_transpose(
        self, psi: "np.ndarray | AngularFlux",
    ) -> "np.ndarray | AngularFlux":
        r"""Adjoint action :math:`C^*\,\psi = \sigma\cdot\psi`.

        Equal to :meth:`apply` — collision is self-adjoint (diagonal
        operator). Returned bit-equal to ``apply(psi)``.

        R-1 Step F — typed :class:`AngularFlux` overload propagates
        through :meth:`apply`.
        """
        return self.apply(psi)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`C + X`.

        When ``X`` is a :class:`StreamingOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        with the streaming operator placed first (the canonical
        ``L + C`` ordering for the algebraic identity).  Otherwise
        falls through to the generic :class:`OperatorSum` via the
        mixin.
        """
        if isinstance(other, StreamingOperator):
            return InvertibleOperator(other, self)
        return super().__add__(other)


# ─────────────────────────────────────────────────────────────────────────
# InvertibleOperator — sweep-invertible composite (L + C)
# ─────────────────────────────────────────────────────────────────────────


class InvertibleOperator(OperatorSum):
    r"""Sweep-invertible composite :math:`L + C` carrying ``.solve`` = WDD sweep.

    R-1 Step C (2026-05-19) — the SN-specific algebraic identity

    .. math::

        (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
        \text{WDD sweep}

    has no generic ``(A+B)^{-1}`` formula — :class:`OperatorSum` by
    itself cannot ``solve``.  The WDD sweep IS the inverse algorithm
    for this specific composite — that's the algebraic foundation of
    the entire SN method (Lewis & Miller §3.2; Adams & Larsen 2002
    §III).  :class:`InvertibleOperator` is the specialisation that
    carries the identity at the type level: it inherits the
    :class:`OperatorSum` ``apply`` (the sum of the operand actions)
    and adds ``solve`` via :func:`~orpheus.sn.sweep.transport_sweep`.

    Construction
    ============

    Two equivalent paths:

    * **Operator algebra dispatch** — ``L + C`` where ``L`` is a
      :class:`StreamingOperator` and ``C`` is a
      :class:`CollisionOperator` returns this class automatically (see
      :meth:`StreamingOperator.__add__` and
      :meth:`CollisionOperator.__add__`).  The composite reads as math.
    * **Explicit construction** — ``InvertibleOperator(L, C)``.  Useful
      when composing variants such as
      ``InvertibleOperator(L_leaf, CollisionOperator(σ_r))`` where
      ``σ_r = σ_t - Σ_{s,0}^{g→g}`` is the removal cross-section that
      lets one fold the within-group self-scatter into the diagonal
      collision term (Adams & Larsen 2002 §III; tracked by issue
      `#200 <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_).

    The two paths produce structurally identical objects — the choice
    only changes the call-site readability.

    Capability set
    ==============

    ``frozenset({CAP_APPLY, CAP_SOLVE})`` — adds ``solve`` to the
    parent :class:`OperatorSum`'s ``apply``-only set.
    ``apply_transpose`` is reserved for Phase H (adjoint propagation
    through the composite).

    The ``.solve`` API
    ==================

    The ``rhs`` parameter is a typed
    :class:`~orpheus.sn.angular_flux.AngularFlux` carrying:

    * ``rhs.values`` — per-ordinate source ``(N, ng, nx, ny)``.  This
      is treated as the per-ordinate anisotropic source
      :math:`Q^{\rm aniso}` that the sweep consumes (the isotropic
      source is zero).
    * ``rhs.boundary`` — face source / BC inflow trace.  Typically
      zero for volumetric SI/Krylov sources (which carry no face
      contribution); the persistent reflective-BC state lives on the
      :class:`SNMesh` and is handled inside the sweep.
    * ``rhs(1)`` (the lag-1 frame, when ``rhs.history_depth >= 2``) —
      the previous iterate that GENERATED this source.  Threaded as
      :pydata:`initial_guess` to :func:`transport_sweep` so the
      curvilinear sweep can read the Carlson coupled-pole seed from
      it.  ``None`` (cold start) → the sweep falls back to its
      in-iteration-source default.

    Parameters
    ----------
    streaming : StreamingOperator
        :math:`L = \Omega\cdot\nabla + \text{angular redistribution}`.
        Resolution A subtractive form:
        ``L.apply(ψ) = M(ψ; σ_t) - σ_t ⊙ ψ_cell``.
    diagonal : CollisionOperator
        :math:`C = \sigma\cdot`.  Its ``.sigma`` attribute is the
        per-cell per-group coefficient used by the sweep (canonically
        ``σ_t``; can be ``σ_r`` for the foldable variant).

    Notes
    -----
    The validation ``σ > 0`` everywhere guards against the
    ``σ_r < 0`` case that can arise when within-group self-scatter
    exceeds total cross-section (rare; not physically meaningful but
    mathematically possible for ill-conditioned multi-group sets).
    The sweep would emit NaN at those cells — surfacing the
    inconsistency at construction is friendlier.
    """

    def __init__(
        self,
        streaming: "StreamingOperator",
        diagonal: "CollisionOperator",
    ) -> None:
        if not isinstance(streaming, StreamingOperator):
            raise TypeError(
                f"InvertibleOperator: 'streaming' must be a "
                f"StreamingOperator; got {type(streaming).__name__}."
            )
        if not isinstance(diagonal, CollisionOperator):
            raise TypeError(
                f"InvertibleOperator: 'diagonal' must be a "
                f"CollisionOperator; got {type(diagonal).__name__}."
            )
        if streaming.sn_mesh is not diagonal.sn_mesh:
            raise ValueError(
                "InvertibleOperator: streaming and diagonal operators "
                "must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if not np.all(diagonal.sigma > 0):
            min_sigma = float(np.min(diagonal.sigma))
            raise ValueError(
                f"InvertibleOperator: diagonal coefficient must be "
                f"strictly positive everywhere for the WDD sweep to be "
                f"well-defined; got min(sigma) = {min_sigma:.3e}.  If "
                f"sigma_r = sigma_t - Sigma_(s,0)^(g->g) is dipping "
                f"negative, the multi-group cross-section set is "
                f"physically inconsistent."
            )
        super().__init__(streaming, diagonal)
        # OperatorSum.__init__ set capabilities = {CAP_APPLY, ...};
        # we add CAP_SOLVE because this composite IS sweep-invertible.
        self.capabilities = self.capabilities | frozenset({CAP_SOLVE})

    # ── Convenience accessors ─────────────────────────────────────────

    @property
    def streaming(self) -> "StreamingOperator":
        """The streaming operand (alias for ``self.a``)."""
        return self.a  # type: ignore[return-value]

    @property
    def diagonal(self) -> "CollisionOperator":
        """The diagonal-collision operand (alias for ``self.b``)."""
        return self.b  # type: ignore[return-value]

    @property
    def sn_mesh(self) -> "SNMesh":
        """The shared :class:`SNMesh` (validated mesh-identity at init)."""
        return self.streaming.sn_mesh

    @property
    def sigma(self) -> np.ndarray:
        r"""The diagonal coefficient used by ``solve`` (σ_t or σ_r)."""
        return self.diagonal.sigma

    # ── solve: WDD sweep ─────────────────────────────────────────────

    def solve(
        self,
        rhs: "AngularFlux",
        *,
        initial_guess: "AngularFlux | None" = None,
    ) -> "AngularFlux":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` via the WDD sweep.

        The cell-balance equation
        :math:`(\Omega\cdot\nabla + \sigma)\,\psi = Q` is integrated
        cell-by-cell in inflow-to-outflow order; the angular closure
        (Cartesian → identity, curvilinear → Morel-Montry) is bound
        on the mesh.

        Parameters
        ----------
        rhs : AngularFlux
            * ``rhs.values`` — per-ordinate source :math:`Q^{\rm aniso}`.
            * ``rhs.boundary`` — BC inflow trace (typically zero for
              SI/Krylov volumetric sources; falls back to ``rhs.boundary``
              when no ``initial_guess`` is supplied).
        initial_guess : AngularFlux or None, optional
            Previous iterate :math:`\psi^{(k-1)}` for the curvilinear
            Carlson coupled-pole seed (M-M reads it as the level's
            angular flux to derive :math:`\bar Q` at :math:`\mu = -1`).
            Its ``.boundary`` also seeds the reflective-BC partner-flux
            trace when present (the previous outflow IS the partner-
            ordinate inflow).  ``None`` (default) → cold start: M-M's
            Carlson seed degenerates to the zero-input result (same
            as ``ZeroSeed`` ablation); BC trace falls back to
            ``rhs.boundary``.

            Explicit kwarg (post-Phase-1.2) — the seed is no longer
            piggy-backed on ``rhs(1)`` (the AngularFlux lag-1 history
            is reserved for time-derivative tracking, an unrelated
            concern).  Outer iterators (:class:`SourceIteration`,
            :class:`KEigenvalue`) pass the previous iterate
            explicitly; GMRES residual calls pass ``None`` so the
            preconditioner doesn't silently route through stateful
            seed state (closes the R-1 Step D silent-fallback bug
            class).

        Returns
        -------
        AngularFlux
            The angular flux satisfying :math:`(L + C)\,\psi =
            \text{rhs}`, with the sweep's outflow face state in
            ``.boundary`` and ``history_depth`` inherited from
            ``rhs.history_depth``.
        """
        from .angular_flux import AngularFlux
        from .sources import PerOrdinateSource
        from .sweep import transport_sweep

        if not isinstance(rhs, AngularFlux):
            raise TypeError(
                f"InvertibleOperator.solve: 'rhs' must be AngularFlux; "
                f"got {type(rhs).__name__}.  The typed-flux contract is "
                f"the only supported surface in R-1; bare-ndarray rhs "
                f"goes through the legacy SNSolver paths."
            )
        if rhs.mesh is not self.sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve: rhs and operator must share "
                "the same SNMesh instance (mesh-identity invariant)."
            )

        sn_mesh = self.sn_mesh
        # R-1 Step 4 A1 — single per-ordinate source convention.
        # ``rhs.values`` IS per-ordinate density by producer contract
        # (``ScatteringOperator.apply`` typed branch normalises by
        # ``/sum_w`` at the producer boundary; ``q_ext`` external
        # sources project via ``PerOrdinateSource.from_isotropic``).
        # The ``transport_sweep`` accepts a single
        # :class:`PerOrdinateSource` and does NOT apply ``/W``
        # internally; both sides of the producer-consumer contract
        # are per-ordinate density.
        source = PerOrdinateSource(rhs.values, sn_mesh)

        # Boundary-state plumbing for reflective / partner-flux BCs.
        # The sweep mutates ``boundary_flux`` in place (write-through);
        # the persistent partner-flux trace is the CRITICAL state for
        # reflective BCs to converge to the right fixed point.  Two
        # paths into the buffer, in priority order:
        #
        # 1. ``initial_guess.boundary`` — the previous iterate's outflow
        #    trace (which IS the partner flux for reflective BCs).
        # 2. ``rhs.boundary`` — fallback when no ``initial_guess``
        #    (typically zero for volumetric sources; non-zero for
        #    explicit BC-driven solves).
        #
        # Why this matters: for SI on reflective BCs, the matvec
        # iteration produces an outflow that becomes the next inflow.
        # Without threading the partner flux, the sweep sees zero
        # inflow each call → effectively VACUUM, and the fixed point
        # shifts away from the true reflective answer.
        boundary_buf = sn_mesh.zeros_boundary_flux()
        if initial_guess is not None:
            _copy_boundary_face_state(initial_guess.boundary, boundary_buf)
        else:
            _copy_boundary_face_state(rhs.boundary, boundary_buf)

        angular, _scalar = transport_sweep(
            source,
            self.sigma,
            sn_mesh,
            boundary_buf,
            initial_guess=initial_guess,
        )

        return AngularFlux(
            angular,
            sn_mesh,
            boundary=boundary_buf,
            history_depth=rhs.history_depth,
        )


def _copy_boundary_face_state(
    src: "BoundaryFlux", dst: "BoundaryFlux",
) -> None:
    """Copy non-None face buffers from ``src`` into ``dst`` in place.

    Both buffers must share the same :class:`SNMesh` (no shape
    validation needed — same mesh implies same buffer sizes).  Used by
    :meth:`InvertibleOperator.solve` to seed the sweep's mutable
    write-through buffer from the immutable rhs.boundary state.
    """
    if src.xmin_face is not None and dst.xmin_face is not None:
        dst.xmin_face[...] = src.xmin_face
    if src.xmax_face is not None and dst.xmax_face is not None:
        dst.xmax_face[...] = src.xmax_face
    if src.xmin_xmax_buf is not None and dst.xmin_xmax_buf is not None:
        dst.xmin_xmax_buf[...] = src.xmin_xmax_buf
    if src.ymin_ymax_buf is not None and dst.ymin_ymax_buf is not None:
        dst.ymin_ymax_buf[...] = src.ymin_ymax_buf
