r"""Unified S\ :sub:`N` transport sweep parameterized by a cell-update strategy.

Round 2 of Wave D of the SN reshape campaign (Issue #161).  Rewrites
the historical 4-path sweep dispatch (slab GL fast path / 2-D wavefront
/ spherical / cylindrical, each with its own inlined per-cell algebra)
into a unified algorithm that branches on a single boolean
(:attr:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.requires_upstream_angular_state`)
exposed by :class:`~orpheus.sn.geometry.SNMesh.reduced` (the
:class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
the mesh consumes from Wave B Issue #6 / Wave D Round 1):

* **Cartesian** (``requires_upstream_angular_state is False``)
  — 1-D + 2-D Cartesian.  Preserves the historical 1-D cumprod fast
  path (algebraic identity for Diamond Difference, gated by GL1D +
  isotropic) inside the unified Cartesian sweep; falls back to the
  2-D wavefront diagonal scheduling otherwise.
* **Curvilinear** (``requires_upstream_angular_state is True``)
  — spherical + cylindrical via the same μ-marching algorithm,
  parameterized by ``sn_mesh.cell_update.update(streaming_terms,
  total_xs, source, upstream_state)`` from the Wave C
  :class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol.  The
  cylindrical per-level loop structure is preserved.

The bit-identical contract — non-negotiable
===========================================

When ``sn_mesh.cell_update`` is :class:`DiamondDifference` (the
default), the 11 frozen regression snapshots at
``tests/sn/regression/snapshots/`` MUST stay ``np.array_equal``-bit-
identical to the pre-rewrite baseline.  The contract is preserved
because:

* The 1-D cumprod fast path is reused **verbatim** as a vectorised
  optimisation that is algebraically equivalent to a per-cell DD
  call sequence — the algebra is preserved at the operation level
  (same operands in the same order), so the per-cell average flux
  matches the per-cell DD output to 1 ULP.
* The 2-D wavefront path is reused **verbatim** with the inlined DD
  math for now.  Wave C-extension's introduction of LD / EC / Step
  strategies will require parameterising 2-D wavefront via
  ``cell_update.update()`` while preserving anti-diagonal
  vectorisation; that work is out of scope here per the
  algebra-of-record discipline (no behaviour change beyond unified
  dispatch for Wave D).
* The curvilinear sweeps invoke ``sn_mesh.cell_update.update(...)``
  per cell.  :class:`DiamondDifference` is a bit-identical
  extraction of the inlined sweep math (Wave C verified via
  hand-calc tests using ``np.array_equal``); the per-cell
  orchestration here matches the pre-rewrite construction order
  for ``streaming_terms``, ``total_xs``, ``source``, and
  ``upstream_state`` so 1-ULP equality holds.

ERR-026 stays open through Wave D
=================================

The curvilinear-sweep one-directional WDD closure
:math:`\psi_{n+1/2} = (\overline{\psi} - (1-\tau_{mm})\,
\psi_{n-1/2})/\tau_{mm}` is reproduced **bit-identically** by
:class:`DiamondDifference` (extracted from the original sweep).
The Wave-E migration of :func:`SNSolver.solve_sn_fixed_source` to
Krylov-on-``apply`` is what closes ERR-026; the unified sweep
preserves the bug bit-identically as the dispatch consolidation
keeps backward compatibility intact.

References
==========

* Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
  *A piecewise linear finite element discretization of the diffusion
  equation for arbitrary polyhedral grids.*  JCP 227, 3738–3757.
  Eq. 50 (sphere/cylinder dome recursion); Eq. 74 (Morel–Montry
  closure).
* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §4.5 (curvilinear DD); §5.3 (DD, weighted-DD,
  Step, LD); §6.4 (sweep ordering).

See also
========

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the
  Wave C Protocol the curvilinear sweeps dispatch through.
* :class:`~orpheus.sn.spatial.diamond.DiamondDifference` — the
  default cell-update strategy; bit-identical extraction of the
  inlined sweep math.
* :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
  — the geometry-layer primitive whose
  ``requires_upstream_angular_state`` boolean drives the dispatch.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry import CoordSystem

from .spatial.cell_update import UpstreamState
from .spatial.diamond import DiamondDifference

if TYPE_CHECKING:
    from .geometry import SNMesh


# ═══════════════════════════════════════════════════════════════════════
# Top-level dispatch
# ═══════════════════════════════════════════════════════════════════════

def transport_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.cell_update``
    (defaults to :class:`DiamondDifference`).

    Dispatch is on
    ``sn_mesh.reduced.requires_upstream_angular_state``: ``False`` for
    slab / 2-D Cartesian (no angular redistribution between successive
    half-angles), ``True`` for spherical / cylindrical (μ-marching with
    angular redistribution).  For 2-D Cartesian where
    ``sn_mesh.reduced is None`` the dispatch falls through to
    Cartesian.

    Parameters
    ----------
    Q : (nx, ny, ng) isotropic source density.
    sig_t : (nx, ny, ng) total macroscopic cross section.
    sn_mesh : SNMesh — augmented geometry with precomputed stencil,
        resolved boundary conditions, and a ``cell_update`` strategy.
    psi_bc : mutable dict storing persistent boundary fluxes
        between outer iterations.
    Q_aniso : (N, nx, ny, ng) per-ordinate anisotropic source (P1+
        scattering).  ``None`` for isotropic-only (P0).

    Returns
    -------
    angular_flux : (N, nx, ny, ng) angular flux per ordinate.
    scalar_flux : (nx, ny, ng) = Σ_n w_n ψ_n.
    """
    reduced = sn_mesh.reduced
    if reduced is not None and reduced.requires_upstream_angular_state:
        return _curvilinear_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
    return _cartesian_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


# ═══════════════════════════════════════════════════════════════════════
# Cartesian dispatch (1-D cumprod fast path or 2-D wavefront)
# ═══════════════════════════════════════════════════════════════════════

def _cartesian_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian dispatch: 1-D cumprod fast path or 2-D wavefront.

    The 1-D cumprod fast path is gated by three preconditions:

    * ``sn_mesh.cell_update`` is :class:`DiamondDifference` — the fast
      path is an algebraic identity that holds **only** for DD
      (Lewis & Miller §5.3; the recurrence
      :math:`\\psi_{\\rm out} = a\\psi_{\\rm in} + b Q` is DD-specific).
    * Quadrature is 1-D Gauss–Legendre (``ny == 1`` and all ``mu_y``
      vanish to machine precision).
    * Source is isotropic (``Q_aniso is None``).

    All three preconditions are needed for the cumprod identity to
    hold; if any fails, the 2-D wavefront path runs (which handles 1-D
    cases as a special case of 2-D).
    """
    quad = sn_mesh.quad
    ny = sn_mesh.ny

    is_gl_1d = (
        ny == 1
        and np.all(np.abs(quad.mu_y) < 1e-15)
        and Q_aniso is None
    )
    cell_update_is_dd = isinstance(sn_mesh.cell_update, DiamondDifference)

    if is_gl_1d and cell_update_is_dd:
        return _sweep_1d_cumprod(Q, sig_t, sn_mesh, psi_bc)
    return _sweep_2d_wavefront(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


# ═══════════════════════════════════════════════════════════════════════
# Curvilinear dispatch (spherical or cylindrical)
# ═══════════════════════════════════════════════════════════════════════

def _curvilinear_sweep(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Curvilinear dispatch: spherical or cylindrical μ-marching.

    Both branches dispatch per-cell to ``sn_mesh.cell_update.update(
    streaming_terms, total_xs, source, upstream_state)``.  When
    ``cell_update`` is :class:`DiamondDifference` (the default), the
    output is bit-identical to the historical inlined sweep math
    (Wave C verified).  Wave C-extension will swap in LD / EC / Step
    strategies without rewriting these orchestration loops.
    """
    if sn_mesh.reduced.coord is CoordSystem.SPHERICAL:
        return _sweep_1d_spherical(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
    # ``requires_upstream_angular_state == True`` AND ``coord !=
    # SPHERICAL`` ⇒ cylindrical.
    return _sweep_1d_cylindrical(Q, sig_t, sn_mesh, psi_bc, Q_aniso)


# ═══════════════════════════════════════════════════════════════════════
# 1D cumprod path (fast, for DD + GL quadrature on slab + isotropic)
# ═══════════════════════════════════════════════════════════════════════

def _sweep_1d_cumprod(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """1D sweep using cumulative products for the DD recurrence.

    Uses the precomputed streaming stencil from SNMesh:
    streaming_x[n, i] = 2|μ_x[n]| / dx[i].

    This path is an **algebraic identity** for Diamond Difference: the
    cumprod-of-``a`` + cumsum-of-``s/cumprod_a`` evaluates the per-cell
    DD recurrence over a whole row in vectorised numpy ops.  The
    operation order matches :func:`_solve_recurrence` /
    :func:`_outgoing` so the per-cell DD call sequence (Wave C
    :class:`DiamondDifference`) and this fast path produce identical
    bits.  The 1-D cumprod path is reused verbatim from the
    pre-Wave-D sweep at the operation level — preserving 1 ULP
    equality is the bit-identity contract for the campaign.
    """
    dx = sn_mesh.dx
    nx = len(dx)
    ng = Q.shape[2]
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights
    ref_x = quad.reflection_index("x")

    # Squeeze out the ny=1 dimension for 1D arrays
    Q_1d = Q[:, 0, :]        # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]  # (nx, ng)

    # Precompute DD face-flux recurrence coefficients for positive directions.
    # The diamond-difference face-flux recurrence is
    #     ψ_out = a·ψ_in + b·(Q/W)
    # where (see ``orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence``
    # for the symbolic derivation of Eq. dd-recurrence)
    #     a = (2μ − Δx·Σ_t) / (2μ + Δx·Σ_t)
    #     b = 2·Δx / (2μ + Δx·Σ_t)
    # and W = Σ w_n is the quadrature weight sum, needed because the
    # isotropic scalar source Q produced by :meth:`SNSolver._add_scattering_source`
    # is in scalar-flux units — the per-ordinate transport equation sees
    # Q/W as its right-hand side. (See ``_sweep_2d_wavefront``, which
    # applies the same ``weight_norm = 1/W`` factor via ``Q_scaled``.)
    n_half = N // 2
    mu_pos = np.abs(quad.mu_x[N // 2:])  # positive half
    w_pos = weights[N // 2:]
    weight_norm = 1.0 / weights.sum()

    denom = 2.0 * mu_pos[:, None, None] + dx[None, :, None] * sig_t_1d[None, :, :]
    stream_coeff = (
        2.0 * mu_pos[:, None, None] - dx[None, :, None] * sig_t_1d[None, :, :]
    ) / denom
    source_coeff = (2.0 * weight_norm) * dx[None, :, None] / denom

    # Initialize boundary fluxes
    if "bc_1d" not in psi_bc:
        psi_bc["bc_1d"] = {
            "left": np.zeros((n_half, ng)),
            "right": np.zeros((n_half, ng)),
        }
    bc = psi_bc["bc_1d"]

    angular_flux = np.zeros((N, nx, 1, ng))
    phi = np.zeros((nx, ng))
    bQ = source_coeff * Q_1d[None, :, :]

    # Tensor-decomposed BCs (R = Σ_α G_α ⊗ A_α). The factory in
    # ``SNMesh.BOUNDARY_OPERATOR_REGISTRY`` resolved the geometry-declared BC into a
    # :class:`BoundaryOperator` whose ``apply`` we call below.
    #
    # Reassemble the full per-ordinate outgoing-flux face arrays from
    # the persistent ``bc`` storage (which stores positive-half values
    # using a *post-reflection* indexing convention: ``bc["left"][n]``
    # is the outgoing flux at the partner ordinate
    # :math:`N - 1 - (n_{\text{half}} + n) = n_{\text{half}} - 1 - n`).
    # Calling ``apply`` then lets specular reflection reach
    # back to its partner via ``reflection_index("x")``; vacuum returns
    # zeros. Bit-equality with the previous string-kind dispatch is
    # preserved because this is an algebraic re-expression of the same
    # operation: ``SpecularBoundaryOperator(axis="x").apply(psi, quad)[k]``
    # is ``psi[ref_x[k]]``, which for GL is ``psi[N-1-k]``.
    bc_left_obj = sn_mesh.bc_left
    bc_right_obj = sn_mesh.bc_right

    # Per-ordinate outgoing-flux buffers at each face. The cumprod
    # convention writes positive-half values into ``bc["left"]`` /
    # ``bc["right"]`` indexed by the *positive* sweep counter ``n``;
    # we mirror that into the full-N face buffers used by the BC
    # protocol, taking care to lay each value at its outgoing-side
    # index so :meth:`SpecularBoundaryOperator.apply` retrieves it via
    # ``reflection_index("x")`` at the matching incoming ordinate.
    psi_face_left_out = np.zeros((N, ng))
    psi_face_right_out = np.zeros((N, ng))
    for n in range(n_half):
        # ``bc["left"][n]`` was written by the previous-iteration
        # backward sweep — it is the outgoing flux at the *left*
        # face for ordinate ``n_half - 1 - n`` (negative side).
        psi_face_left_out[n_half - 1 - n] = bc["left"][n]
        # ``bc["right"][n]`` was written by the previous-iteration
        # forward sweep — it is the outgoing flux at the *right* face
        # for ordinate ``n_half + n`` (positive side).
        psi_face_right_out[n_half + n] = bc["right"][n]

    psi_face_left_in = bc_left_obj.apply(psi_face_left_out, quad)

    for n in range(n_half):
        a = stream_coeff[n]  # (nx, ng)
        s = bQ[n]            # (nx, ng)

        # Forward sweep (positive direction, left → right):
        # incoming at the *left* face for ordinate ``n_half + n``.
        psi0_left = psi_face_left_in[n_half + n]
        psi_fwd = _solve_recurrence(a, s, psi0_left)
        bc["right"][n, :] = _outgoing(a, s, psi0_left)
        # Update the right-face outgoing buffer in place: the backward
        # sweep at this ``n`` consumes the just-written outgoing via
        # the BC operator (matches the original in-iteration coupling
        # where right-reflective backward read ``bc["right"][n]``
        # immediately after the forward sweep wrote it).
        psi_face_right_out[n_half + n] = bc["right"][n]
        psi_face_right_in = bc_right_obj.apply(
            psi_face_right_out, quad,
        )
        phi += w_pos[n] * psi_fwd
        angular_flux[n_half + n, :, 0, :] = psi_fwd

        # Backward sweep (negative direction via reversal, right → left):
        # incoming at the *right* face for ordinate ``n_half - 1 - n``.
        psi0_right = psi_face_right_in[n_half - 1 - n]
        psi_bwd = _solve_recurrence(a[::-1], s[::-1], psi0_right)
        bc["left"][n, :] = _outgoing(a[::-1], s[::-1], psi0_right)
        phi += w_pos[n] * psi_bwd[::-1]
        angular_flux[n_half - 1 - n, :, 0, :] = psi_bwd[::-1]

    return angular_flux, phi[:, None, :]  # restore ny=1 dim


def _solve_recurrence(
    a: np.ndarray, s: np.ndarray, psi0: np.ndarray,
) -> np.ndarray:
    """Solve DD recurrence via cumulative products. Returns cell-average flux."""
    nc = a.shape[0]
    cp = np.cumprod(a, axis=0)
    cs = np.cumsum(s / cp, axis=0)

    psi_in = np.empty_like(a)
    psi_in[0] = psi0
    if nc > 1:
        psi_in[1:] = cp[:-1] * (psi0[None, :] + cs[:-1])

    psi_out = a * psi_in + s
    return 0.5 * (psi_in + psi_out)


def _outgoing(
    a: np.ndarray, s: np.ndarray, psi0: np.ndarray,
) -> np.ndarray:
    """Outgoing flux at the end of a forward sweep."""
    cp = np.cumprod(a, axis=0)
    cs = np.cumsum(s / cp, axis=0)
    return cp[-1] * (psi0 + cs[-1])


# ═══════════════════════════════════════════════════════════════════════
# 1D spherical path (cell-by-cell with angular redistribution)
# ═══════════════════════════════════════════════════════════════════════

def _sweep_1d_spherical(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Spherical 1-D sweep with geometry-weighted angular redistribution.

    Processes ordinates sequentially from most negative :math:`\mu` to
    most positive, dispatching per-cell to
    ``sn_mesh.cell_update.update(...)`` for the closure algebra.

    The balance equation includes a geometry factor
    :math:`\Delta A_i / w_n` on the redistribution term
    (Bailey et al. 2009), ensuring per-ordinate flat-flux consistency:

    .. math::

        \psi_{n,i} = \frac{S_i V_i + |\mu_n|(A_{\rm in}+A_{\rm out})
            \psi^s_{\rm in} + \frac{\Delta A_i}{w_n}
            (\alpha_{n+\frac12}+\alpha_{n-\frac12})\psi_{n-\frac12}}
            {2|\mu_n| A^s_{\rm out}
            + 2\frac{\Delta A_i}{w_n}\alpha_{n+\frac12} + \Sigma_t V_i}

    The :math:`\alpha` coefficients are computed as
    :math:`\alpha_{n+1/2} = \alpha_{n-1/2} - w_n \mu_n` and form a
    non-negative dome for μ-sorted ordinates.

    The cell-update math itself lives in
    ``sn_mesh.cell_update.update(...)`` — this orchestrator only
    handles the sweep ordering, boundary conditions, and source
    construction.  When ``cell_update`` is :class:`DiamondDifference`
    (the default), the output is bit-identical to the pre-Wave-D
    inlined math (Wave C verified).
    """
    nx = sn_mesh.nx
    ng = Q.shape[2]
    quad = sn_mesh.quad
    N = quad.N
    mu = quad.mu_x
    weights = quad.weights

    Q_1d = Q[:, 0, :]          # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]  # (nx, ng)

    # Read connection coefficients via the canonical ReducedStreamingOperator
    # accessors (NOT the deprecated SNMesh.alpha_half / .redist_dAw / .tau_mm
    # properties, which emit DeprecationWarning).  The arrays here are the
    # same objects the deprecated properties would return — bit-identical
    # numerical values.
    reduced = sn_mesh.reduced
    A = reduced.face_areas       # (nx+1,) surface areas at cell faces
    V = sn_mesh.volumes[:, 0]    # (nx,) cell volumes
    cell_update = sn_mesh.cell_update

    # Boundary condition at the outer face (r = R) is a tensor-
    # decomposed :class:`~orpheus.geometry.boundary.BoundaryOperator`. The
    # spherical sweep stores the outgoing flux per ordinate in
    # ``bc_outer`` (indexed by *outgoing* — positive μ — ordinate) and
    # reads the incoming flux for negative μ via
    # ``apply``. For ``VacuumBoundaryOperator`` this returns zeros; for
    # ``SpecularBoundaryOperator(axis="x")`` it returns ``bc_outer[ref[n]]``,
    # which is bit-identical to the previous ``bc_outer[ref[n]].copy()``
    # call site. The incoming buffer is recomputed each iteration so
    # the loop's own outgoing updates feed back through the BC operator.
    bc_outer_obj = sn_mesh.bc_right

    # Persistent boundary flux at the outer face (per ordinate, indexed
    # by the *outgoing* — positive-μ — ordinate). Negative-μ entries
    # remain zero throughout the sweep; they are placeholders so the
    # array can be passed whole into ``apply``.
    if "bc_sph" not in psi_bc:
        psi_bc["bc_sph"] = np.zeros((N, ng))
    bc_outer = psi_bc["bc_sph"]

    # Angular "face flux" between successive ordinates: ψ_{n-1/2, i}
    # Shape (nx, ng). Initialised to zero for the first ordinate (α_{1/2}=0).
    psi_angle = np.zeros((nx, ng))

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    # Isotropic source → angular source density by dividing by sum(w)
    # Then multiply by cell volume for the balance equation
    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm  # (nx, ng)

    # Per-ordinate anisotropic source (MMS external / P1+), if present
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_1d = Q_aniso[:, :, 0, :] * weight_norm  # (N, nx, ng)

    for n in range(N):
        mu_n = mu[n]
        w_n = weights[n]

        # Per-ordinate volumetric source
        QV = QV_iso
        if has_aniso:
            QV = QV_iso + Q_aniso_1d[n] * V[:, None]

        # Set up incoming spatial flux at the start of the sweep.
        # For inward sweeps (μ < 0), boundary → centre: read incoming
        # from the outer-face BC.  For outward (μ ≥ 0), at r = 0 we
        # have A[0] = 4π(0)² = 0, so no spatial incoming flux.
        if mu_n < 0:
            # ``apply`` is bit-identical to the previous
            # ``bc_outer[ref[n]].copy()`` indexing for SpecularBoundaryOperator and
            # zeros for VacuumBoundaryOperator.
            psi_in_full = bc_outer_obj.apply(bc_outer, quad)
            psi_spatial_in = psi_in_full[n]
        else:
            psi_spatial_in = np.zeros(ng)

        # Iterate cells in topological-sort order for this ordinate.
        # The CellVisit packet pre-resolves the sweep direction —
        # face_area_downstream is the outgoing face (outer for
        # outward, inner for inward); spatial_upstream below is the
        # flux flowing INTO the cell.  No sign-of-μ branching here.
        for visit in sn_mesh.iter_cell_visits(ordinate_idx=n):
            i = visit.cell_idx
            upstream = UpstreamState(
                spatial_upstream=psi_spatial_in,
                angular_upstream=psi_angle[i],
            )
            result = cell_update.update(
                visit=visit,
                total_xs=sig_t_1d[i],
                source=QV[i],
                upstream_state=upstream,
            )

            psi = result.cell_average_flux
            psi_angle[i] = result.outgoing_angular_state

            angular_flux[n, i, 0, :] = psi
            scalar_flux[i] += w_n * psi

            # ``outgoing_spatial_flux`` is always populated for the
            # spherical (non-degenerate) curvilinear branch.
            psi_spatial_in = result.outgoing_spatial_flux

        # Store outgoing flux at outer boundary for reflective BC,
        # only on outward sweeps — this is the last visit's outgoing
        # face flux on cell nx-1.
        if mu_n >= 0:
            bc_outer[n] = psi_spatial_in

    return angular_flux, scalar_flux[:, None, :]  # restore ny=1 dim


# ═══════════════════════════════════════════════════════════════════════
# 1D cylindrical path (per-level azimuthal redistribution)
# ═══════════════════════════════════════════════════════════════════════

def _sweep_1d_cylindrical(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Cylindrical 1-D sweep with geometry-weighted azimuthal redistribution.

    For each μ-level *p*, processes azimuthal ordinates sequentially
    from most-inward (:math:`\eta = -\sin\theta`) to most-outward
    (:math:`\eta = +\sin\theta`), dispatching per-cell to
    ``sn_mesh.cell_update.update(...)`` for the closure algebra.

    The balance equation includes a geometry factor
    :math:`\Delta A_i / w_m` on the redistribution term
    (Bailey et al. 2009), ensuring per-ordinate flat-flux consistency:

    .. math::

        \psi_{m,i} = \frac{S_i V_i + |\eta_m|(A_{\rm in}+A_{\rm out})
            \psi^s_{\rm in} + \frac{\Delta A_i}{w_m}
            (\alpha_{m+\frac12}+\alpha_{m-\frac12})\psi_{m-\frac12}}
            {2|\eta_m| A^s_{\rm out}
            + 2\frac{\Delta A_i}{w_m}\alpha_{m+\frac12} + \Sigma_t V_i}

    The :math:`\alpha` coefficients are computed from the radial
    direction cosine :math:`\eta` (Bailey et al. Eq. 50) and form a
    non-negative dome, so the denominator is unconditionally positive.

    The cell-update math itself lives in
    ``sn_mesh.cell_update.update(...)``.  The pure-azimuthal degenerate
    case (``abs_mu < 1e-15``) is signalled by the strategy via
    ``CellResult.outgoing_spatial_flux is None`` — the orchestrator
    skips the face-flux update accordingly.
    """
    nx = sn_mesh.nx
    ng = Q.shape[2]
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights

    Q_1d = Q[:, 0, :]          # (nx, ng)
    sig_t_1d = sig_t[:, 0, :]  # (nx, ng)

    reduced = sn_mesh.reduced
    A = reduced.face_areas       # (nx+1,) = 2πr at edges
    V = sn_mesh.volumes[:, 0]    # (nx,) cell volumes
    cell_update = sn_mesh.cell_update

    # Boundary condition at the outer face (r = R) is a tensor-
    # decomposed :class:`~orpheus.geometry.boundary.BoundaryOperator`. The
    # cylindrical sweep stores per-ordinate outgoing flux in
    # ``bc_outer`` (indexed by global ordinate) and reads incoming
    # via ``apply``. For ``VacuumBoundaryOperator`` returns zeros; for
    # ``SpecularBoundaryOperator(axis="x")`` returns ``bc_outer[ref[n]]`` — bit-
    # identical to the previous ``bc_outer[ref[n]].copy()`` call site.
    bc_outer_obj = sn_mesh.bc_right

    # Persistent boundary flux at the outer face (per ordinate)
    if "bc_cyl" not in psi_bc:
        psi_bc["bc_cyl"] = np.zeros((N, ng))
    bc_outer = psi_bc["bc_cyl"]

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    # Isotropic source → angular source density
    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm  # (nx, ng)

    # Per-ordinate anisotropic source (MMS external / P1+), if present
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_1d = Q_aniso[:, :, 0, :] * weight_norm  # (N, nx, ng)

    # Process each μ-level independently
    for p, level_idx in enumerate(quad.level_indices):
        M = len(level_idx)

        # Azimuthal "face flux" between successive ordinates on this level.
        # Initialised to zero: α_{1/2} = 0 so the product α·ψ vanishes.
        psi_angle = np.zeros((nx, ng))

        for m_local in range(M):
            n = level_idx[m_local]  # global ordinate index
            eta_n = quad.mu_x[n]    # radial direction cosine
            abs_eta = abs(eta_n)
            w_n = weights[n]

            # Per-ordinate volumetric source
            QV = QV_iso
            if has_aniso:
                QV = QV_iso + Q_aniso_1d[n] * V[:, None]

            # Set up incoming spatial flux at the start of the sweep
            # for this ordinate.  Inward (η < 0): read from outer-face
            # BC.  Outward (η > 0): zero at r = 0.  Degenerate
            # (|η| < 1e-15): unused by the strategy; pass zeros.
            if eta_n < 0:
                psi_in_full = bc_outer_obj.apply(bc_outer, quad)
                psi_spatial_in = psi_in_full[n]
            else:
                psi_spatial_in = np.zeros(ng)

            # Iterate cells in DAG-topological order for this ordinate.
            # iter_cell_visits resolves the cylindrical
            # ``direction_idx = m_local`` + ``mu_level_idx = p`` to
            # the correct global η, populates the geometric
            # StreamingTerms (with the level-resolved global ordinate
            # for abs_mu), and yields the sweep-direction-resolved
            # face_area_downstream for each visit.
            for visit in sn_mesh.iter_cell_visits(
                ordinate_idx=m_local, mu_level_idx=p,
            ):
                i = visit.cell_idx
                upstream = UpstreamState(
                    spatial_upstream=psi_spatial_in,
                    angular_upstream=psi_angle[i],
                )
                result = cell_update.update(
                    visit=visit,
                    total_xs=sig_t_1d[i],
                    source=QV[i],
                    upstream_state=upstream,
                )

                psi = result.cell_average_flux
                psi_angle[i] = result.outgoing_angular_state

                angular_flux[n, i, 0, :] = psi
                scalar_flux[i] += w_n * psi

                # For the cylindrical pure-azimuthal degenerate case,
                # ``result.outgoing_spatial_flux`` is None — no
                # face-flux update.  For non-degenerate sweeps it
                # carries the next cell's incoming spatial flux.
                if result.outgoing_spatial_flux is not None:
                    psi_spatial_in = result.outgoing_spatial_flux

            # Store outgoing at outer boundary for reflective BC —
            # only on outward (non-degenerate) sweeps.  Mirrors the
            # original ``else`` branch (eta_n >= 0 and not degenerate).
            # The inward branch returns above without writing
            # ``bc_outer``; the degenerate case has no spatial flow.
            if eta_n >= 0 and abs_eta >= 1e-15:
                bc_outer[n] = psi_spatial_in

    return angular_flux, scalar_flux[:, None, :]  # restore ny=1 dim


# ═══════════════════════════════════════════════════════════════════════
# 2D wavefront path (vectorized along anti-diagonals; inlined DD math)
# ═══════════════════════════════════════════════════════════════════════
#
# The 2-D wavefront sweep retains its inlined DD math.  Parameterising
# anti-diagonal-vectorised cell updates by the
# :class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol while
# preserving 1-ULP equality is the open architectural problem for
# Wave C-extension's LD / EC / Step rollout — the per-cell call shape
# would need to accept ``(n_diag, ng)`` slice arguments rather than
# the contract's ``(ng,)`` per-cell shape.  Wave D's gating contract
# is bit-identity for ``cell_update = DiamondDifference``; the inlined
# DD math here meets that contract by construction.

def _sweep_2d_wavefront(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: SNMesh,
    psi_bc: dict,
    Q_aniso: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """2D sweep using wavefront parallelism along anti-diagonals.

    Uses the precomputed streaming stencil from SNMesh:
    streaming_x[n, i] = 2|μ_x[n]| / dx[i],
    streaming_y[n, j] = 2|μ_y[n]| / dy[j].

    Boundary conditions are read from ``sn_mesh.bc_xmin`` etc.

    The per-cell DD math is inlined here (rather than dispatched via
    ``sn_mesh.cell_update.update(...)``) because anti-diagonal
    vectorisation operates on ``(n_diag, ng)`` cell slices, a shape
    the per-cell ``CellUpdate`` Protocol does not currently accept.
    Wave C-extension will resolve this when LD / EC / Step strategies
    land — the Protocol's vectorisation contract is the open design
    point.  For Wave D, the inlined DD math is bit-identical to a
    per-cell :class:`DiamondDifference` call sequence by construction.
    """
    dx = sn_mesh.dx
    dy = sn_mesh.dy
    nx, ny, ng = Q.shape
    quad = sn_mesh.quad
    N = quad.N
    mu_x = quad.mu_x
    mu_y = quad.mu_y
    weights = quad.weights
    ref_x = quad.reflection_index("x")
    ref_y = quad.reflection_index("y")

    angular_flux = np.zeros((N, nx, ny, ng))
    scalar_flux = np.zeros((nx, ny, ng))

    # Persistent boundary flux arrays for reflective BCs
    if "bc_2d_x" not in psi_bc:
        psi_bc["bc_2d_x"] = np.zeros((N, nx + 1, ny, ng))
        psi_bc["bc_2d_y"] = np.zeros((N, nx, ny + 1, ng))

    psi_x = psi_bc["bc_2d_x"]  # (N, nx+1, ny, ng) face fluxes in x
    psi_y = psi_bc["bc_2d_y"]  # (N, nx, ny+1, ng) face fluxes in y

    weight_norm = 1.0 / weights.sum()

    # Precompute diagonal indices per sweep direction (4 directions).
    _diag_cache: dict[tuple[int, int], tuple] = {}
    for sx in (-1, 1):
        for sy in (-1, 1):
            ix_arr = np.arange(nx) if sx >= 0 else np.arange(nx - 1, -1, -1)
            iy_arr = np.arange(ny) if sy >= 0 else np.arange(ny - 1, -1, -1)
            diags = []
            for k in range(nx + ny - 1):
                i_start = max(0, k - ny + 1)
                i_end = min(nx - 1, k)
                local_i = np.arange(i_start, i_end + 1)
                local_j = k - local_i
                diags.append((ix_arr[local_i], iy_arr[local_j]))
            _diag_cache[(sx, sy)] = (
                0 if sx >= 0 else 1,   # ix_in
                1 if sx >= 0 else 0,   # ix_out
                0 if sy >= 0 else 1,   # iy_in
                1 if sy >= 0 else 0,   # iy_out
                diags,
            )

    # Precompute scaled source (avoids recomputing per diagonal)
    Q_scaled = Q * weight_norm
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_scaled = Q_aniso * weight_norm  # (N, nx, ny, ng)

    # Precomputed streaming stencil
    str_x = sn_mesh.streaming_x  # (N_ord, nx)
    str_y = sn_mesh.streaming_y  # (N_ord, ny)

    for n in range(N):
        mx = mu_x[n]
        my = mu_y[n]
        w = weights[n]

        # Per-ordinate source: isotropic + anisotropic (if present)
        Q_n = Q_scaled
        if has_aniso:
            Q_n = Q_scaled + Q_aniso_scaled[n]  # (nx, ny, ng)

        if abs(mx) < 1e-15 and abs(my) < 1e-15:
            # Pure z-directed ordinate: no streaming in x or y.
            psi_avg = Q_n / sig_t  # (nx, ny, ng)
            angular_flux[n, :, :, :] = psi_avg
            scalar_flux += w * psi_avg
            continue

        # Look up precomputed diagonal indices for this sweep direction
        key = (1 if mx >= 0 else -1, 1 if my >= 0 else -1)
        ix_in, ix_out, iy_in, iy_out, diags = _diag_cache[key]

        # Apply boundary conditions at incoming faces via the tensor-
        # decomposed :class:`BoundaryOperator` Protocol on each face. For
        # ``VacuumBoundaryOperator`` the result is zero (the buffers stay
        # zero-initialised and the slice assignment is a no-op
        # write); for ``SpecularBoundaryOperator(axis=...)`` the result is
        # ``psi_face[ref_axis[n]]`` — bit-identical to the previous
        # in-place ``psi_x[n, 0] = psi_x[ref_x[n], 0]`` copy.
        if mx >= 0:
            psi_x[n, 0, :, :] = sn_mesh.bc_xmin.apply(
                psi_x[:, 0, :, :], quad,
            )[n]
        else:
            psi_x[n, nx, :, :] = sn_mesh.bc_xmax.apply(
                psi_x[:, nx, :, :], quad,
            )[n]

        if my >= 0:
            psi_y[n, :, 0, :] = sn_mesh.bc_ymin.apply(
                psi_y[:, :, 0, :], quad,
            )[n]
        else:
            psi_y[n, :, ny, :] = sn_mesh.bc_ymax.apply(
                psi_y[:, :, ny, :], quad,
            )[n]

        # Precomputed streaming for this ordinate
        str_x_n = str_x[n]  # (nx,)
        str_y_n = str_y[n]  # (ny,)

        for ii, jj in diags:
            # Gather incoming face fluxes
            psi_in_x = psi_x[n, ii + ix_in, jj, :]   # (n_diag, ng)
            psi_in_y = psi_y[n, ii, jj + iy_in, :]    # (n_diag, ng)

            # Precomputed streaming coefficients for these cells
            sx_ii = str_x_n[ii, None]  # (n_diag, 1)
            sy_jj = str_y_n[jj, None]  # (n_diag, 1)

            # Diamond-difference equation using precomputed stencil
            denom = sig_t[ii, jj, :] + sx_ii + sy_jj

            psi_avg = (
                Q_n[ii, jj, :]
                + sx_ii * psi_in_x
                + sy_jj * psi_in_y
            ) / denom

            # Store outgoing face fluxes for next diagonal
            psi_x[n, ii + ix_out, jj, :] = 2.0 * psi_avg - psi_in_x
            psi_y[n, ii, jj + iy_out, :] = 2.0 * psi_avg - psi_in_y

            # Accumulate angular and scalar flux
            angular_flux[n, ii, jj, :] = psi_avg
            scalar_flux[ii, jj, :] += w * psi_avg

    return angular_flux, scalar_flux
