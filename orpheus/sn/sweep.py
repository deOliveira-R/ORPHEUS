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

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 §3.9.4
  (pp. 141-144) — primary source for the curvilinear S\ :sub:`N`
  cell-balance + DD difference relations the sweep implements.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169.  Auxiliary justification for the M-M clamp.
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
from .spatial.psi_half_angle_seed import carlson_inward_sweep_from_source
from .sweep_graph import OctantLabel

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

    psi_face_left_in = bc_left_obj.apply(psi_face_left_out)

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
        psi_face_right_in = bc_right_obj.apply(psi_face_right_out)
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

    # Isotropic source → angular source density by dividing by sum(w)
    # Then multiply by cell volume for the balance equation
    weight_norm = 1.0 / weights.sum()
    QV_iso = Q_1d * V[:, None] * weight_norm  # (nx, ng)

    # ── Phase G Step 2 Path C (Issue #196) two surgical SI-vs-apply
    # convention fixes —— DOMINANT defect (pole-face WDD IC) + sub-
    # dominant (Carlson seed source).  See
    # ``.claude/agent-memory/numerics-investigator/
    # issue_196_phase_g_step2_minimal_reproducer.md`` for the
    # algebraic walkthrough.  Without these the SI Picard fixed point
    # disagrees with the apply-matvec fixed point by ~22% L0 on the
    # homogeneous reflective sphere streaming-equilibrium test.
    #
    # Both fixes read previous-Picard-iter state cached in ``psi_bc``:
    #
    # - ``psi_bc["psi_pole"]``: (N, ng) cell-centre ψ at i=0 from the
    #   previous outer iteration; carries the Lewis-Miller §4.5
    #   Carlson starting-direction pole-face IC for μ ≥ 0.
    # - ``psi_bc["phi_0_prev"]``: (nx, ng) scalar flux φ_0 from the
    #   previous outer iteration; carries the apply-matvec's
    #   ``Q̄ = 0.5 · Σ_t · φ_0`` Carlson seed source convention.
    #
    # On the FIRST outer iteration neither key exists — both fixes
    # degrade gracefully to the legacy Phase F values (zero pole IC;
    # ``Q̄ = 0.5 · Q_1d``).  Subsequent iterations use the apply-matvec
    # convention; Picard converges geometrically to the apply-matvec
    # fixed point.

    # ── Phase F (Issue #168) Carlson coupled-pole seed for ψ_{1/2,i,g} ──
    # Mirror-image of the apply-matvec Phase D fix
    # (operator.py::transport_operator_matvec_spherical).  Replaces the
    # legacy ``psi_angle = np.zeros((nx, ng))`` Carlson-zero seed —
    # diagnosed as the wrong term initialisation (ERR-026 manifestation
    # #6, twin of the apply-path bug Phase D closed).
    #
    # Per Hébert §3.9.4 Eqs. (3.432)-(3.435), the canonical seed at
    # μ = −1 is the cell-centred output of the inward DD sweep:
    #
    #     φ̄_i = (Δr_i · Q̄_i + 2 · φ̄_{i+1/2})
    #            / (Δr_i · Σ_t,i + 2)
    #     φ̄_{i-1/2} = 2 · φ̄_i - φ̄_{i+1/2}
    #
    # ── Phase G Step 2 Path C Carlson seed source fix (Issue #196) ──
    # Apply-matvec convention (psi_half_angle_seed.py:569):
    # ``Q̄ = 0.5 · Σ_t · φ_0`` where φ_0 is built from the *input* ψ.
    # The SI's previous-iter scalar flux φ_0 plays the same role for
    # the Picard recursion.  At the within-group fixed point
    # ``Σ_t · φ_0 = Q_within`` exactly, so the two conventions agree
    # there; off the FP they differ.  Driving the Carlson seed from
    # ``Σ_t · φ_0_prev`` makes SI Picard consume the same seed as the
    # apply-matvec's matvec, so the two converge to the same fixed
    # point.  ERR-048 manifestation #2 (sub-dominant; ~3% rel
    # residual after the pole-face IC fix alone).
    #
    # ``bc_outer_value``: outer-face inflow at μ = −1.  Apply the BC
    # operator to the persistent outflow buffer ``bc_outer`` (carries
    # the previous outward sweep's outgoing flux per ordinate) and
    # extract the value at the most-inward ordinate.  Linear in the
    # current sweep state, BC-realised — vacuum → 0, reflective → the
    # mirrored outflow at μ = +1.
    most_inward_idx = int(np.argmin(mu))
    inflow_full = bc_outer_obj.apply(bc_outer)            # (N, ng)
    bc_outer_value = inflow_full[most_inward_idx, :]       # (ng,)

    # Build ``Q̄[g, i]`` per the apply-matvec convention when a previous
    # scalar flux is cached; fall back to the Phase F legacy form
    # ``Q_1d / Σw`` on the first outer iteration (cold start).
    #
    # Phase G Step 2 (Issue #196): the previous form ``0.5 · …`` was
    # the sphere-specific encoding of ``1/Σw`` (Σw = 2 for GL on
    # [−1, 1]).  Replaced with the geometry-general ``/ weights.sum()``
    # — bit-identical for sphere; correct for any 3D quadrature.
    # See :mod:`orpheus.sn.spatial.psi_half_angle_seed` for the full
    # Pomraning / Hébert §3.9.4 derivation.
    sigma_t_gx = sig_t_1d.T                                # (ng, nx)
    dr = sn_mesh.dx
    Sw = weights.sum()
    if "phi_0_prev" in psi_bc:
        phi_0_prev = psi_bc["phi_0_prev"]                  # (nx, ng)
        Q_bar = sigma_t_gx * phi_0_prev.T / Sw             # (ng, nx)
    else:
        Q_bar = Q_1d.T / Sw                                # (ng, nx) cold start

    phi_aux = carlson_inward_sweep_from_source(
        Q_bar=Q_bar,
        sigma_t=sigma_t_gx,
        dr=dr,
        bc_outer_value=bc_outer_value,
    )  # (ng, nx)

    # Angular "face flux" between successive ordinates: ψ_{n-1/2, i}.
    # Shape (nx, ng).  Carlson coupled-pole seed for the M-M
    # recurrence — replaces Phase B's hardcoded zero.  Transpose
    # ``(ng, nx) -> (nx, ng)`` to match the sweep's per-cell layout.
    psi_angle = phi_aux.T.copy()                           # (nx, ng)

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

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
        # from the outer-face BC.  For outward (μ ≥ 0), the pole-face
        # IC at i=0 follows Lewis-Miller §4.5 Carlson starting-
        # direction convention.
        if mu_n < 0:
            # ``apply`` is bit-identical to the previous
            # ``bc_outer[ref[n]].copy()`` indexing for SpecularBoundaryOperator and
            # zeros for VacuumBoundaryOperator.
            psi_in_full = bc_outer_obj.apply(bc_outer)
            psi_spatial_in = psi_in_full[n]
        else:
            # ── Phase G Step 2 Path C (Issue #196) pole-face IC fix ──
            # Lewis-Miller §4.5 Carlson starting-direction convention,
            # mirroring operator.py::transport_operator_matvec_spherical
            # lines 781-786.  At i=0 the face area A[0] = 4π·0² = 0, so
            # the streaming-IN term |μ|·A[0]·ψ_face_in is identically
            # zero regardless of ψ_face_in.  But the WDD diamond
            # recurrence ``ψ_face_out = 2·ψ_cell - ψ_face_in``
            # propagates the choice downstream: setting
            # ψ_face_in = ψ_cell[i=0] gives ψ_face_out = ψ_cell on a
            # flat ψ field, preserving Pomraning isotropy at r=0.
            # Legacy ψ_face_in = 0 yields ψ_face_out = 2·ψ_cell which
            # destroys flat-field invariance.  See
            # ``.claude/agent-memory/numerics-investigator/
            # issue_196_phase_g_step2_minimal_reproducer.md``.
            # Cold-start: pre-iter-1 there is no ψ_pole cache; fall
            # back to the legacy zero (matches Phase F behaviour).
            if "psi_pole" in psi_bc:
                psi_spatial_in = psi_bc["psi_pole"][n]      # (ng,)
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

    # ── Phase G Step 2 Path C (Issue #196) cache previous-iter state
    # for the next sweep's Carlson seed source + pole-face IC fixes.
    # Both keys must be supplied by the SAME outer Picard iteration so
    # they consistently descend from one converged sweep result.
    psi_bc["psi_pole"] = angular_flux[:, 0, 0, :].copy()      # (N, ng)
    psi_bc["phi_0_prev"] = scalar_flux.copy()                 # (nx, ng)

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

    # ── Phase F (Issue #168) per-level Carlson coupled-pole seed ────
    # Cylindrical analog of the spherical sweep fix.  Each μ-level's
    # azimuthal recurrence gets its own Hébert §3.9.4 seed.  Per the
    # Phase F diagnostic memo §4.2, cylindrical Gate 1.1 was passing
    # empirically (per-level α-domes telescope back to zero at the
    # level's edges so the wrong zero seed self-cancels) — but
    # Cardinal Rule 2 (architecture) requires the structural fix
    # alongside the spherical one.
    #
    # ── Phase G Step 2 Path C (Issue #196) two-fix backport ──
    # Mirrors the spherical sweep: previous-iter ψ at i=0 (pole-face
    # IC for outward sweeps) + previous-iter φ_0 (Carlson seed source)
    # both cached in ``psi_bc``.  Cold-start: keys absent → fall back
    # to Phase F legacy values.
    # Phase G Step 2 (Issue #196): geometry-general ``1/Σw_full``
    # normalisation replaces the sphere-only literal ``0.5``.  ``phi_0_prev``
    # is the FULL scalar flux (sum over all ordinates), so the full
    # quadrature weight sum ``Sw_full = weights.sum()`` is the matching
    # normalisation.  Sphere: ``Sw_full = 2``, bit-identical to the old
    # form; cylinder ProductQuadrature: ``Sw_full = 4π`` and the old
    # ``0.5`` was wrong by factor ``2π`` (ERR-048 manifestation #3).
    sigma_t_gx = sig_t_1d.T              # (ng, nx)
    dr = sn_mesh.dx
    Sw_full = weights.sum()
    if "phi_0_prev_cyl" in psi_bc:
        phi_0_prev = psi_bc["phi_0_prev_cyl"]            # (nx, ng)
        Q_bar_iso = sigma_t_gx * phi_0_prev.T / Sw_full  # (ng, nx)
    else:
        Q_bar_iso = Q_1d.T / Sw_full     # (ng, nx) — L=0, P_0(−1)=1 cold start
    mu_full = quad.mu_x                  # radial direction cosines
    # Outer-face inflow trace via BC operator: bc_outer carries the
    # outflow per ordinate from prior outward sweeps.  Apply once
    # before the per-level loop (linear, commutes with per-level
    # extraction).
    inflow_full_cyl = bc_outer_obj.apply(bc_outer)  # (N, ng)

    # Process each μ-level independently
    for p, level_idx in enumerate(quad.level_indices):
        M = len(level_idx)
        level_idx_arr = np.asarray(level_idx)
        # Most-inward (smallest η) ordinate within this level.
        within_idx_most_inward = int(np.argmin(mu_full[level_idx_arr]))
        global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
        bc_outer_value_level = inflow_full_cyl[global_idx_most_inward, :]  # (ng,)

        phi_aux_level = carlson_inward_sweep_from_source(
            Q_bar=Q_bar_iso,
            sigma_t=sigma_t_gx,
            dr=dr,
            bc_outer_value=bc_outer_value_level,
        )  # (ng, nx)

        # Azimuthal "face flux" between successive ordinates on this level.
        # Carlson coupled-pole seed (Phase F) replaces the legacy zero.
        psi_angle = phi_aux_level.T.copy()  # (nx, ng)

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
            # BC.  Outward (η > 0): pole-face IC at i=0 per Lewis-Miller
            # §4.5 Carlson starting-direction convention.  Degenerate
            # (|η| < 1e-15): unused by the strategy; pass zeros.
            if eta_n < 0:
                psi_in_full = bc_outer_obj.apply(bc_outer)
                psi_spatial_in = psi_in_full[n]
            else:
                # ── Phase G Step 2 Path C pole-face IC fix (cylindrical) ──
                # Mirror of the spherical fix above; see that branch for
                # the algebraic justification.  At i=0 (axis) the inner
                # face area A[0] = 2π·0 = 0, so the streaming-IN term is
                # identically zero — but the WDD recurrence propagates
                # ψ_face_in downstream.  ψ_face_in = ψ_cell[i=0]
                # preserves Pomraning isotropy at the axis.
                if "psi_pole_cyl" in psi_bc:
                    psi_spatial_in = psi_bc["psi_pole_cyl"][n]   # (ng,)
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

    # ── Phase G Step 2 Path C (Issue #196) cache previous-iter state
    # for the next sweep's Carlson seed source + pole-face IC fixes.
    psi_bc["psi_pole_cyl"] = angular_flux[:, 0, 0, :].copy()  # (N, ng)
    psi_bc["phi_0_prev_cyl"] = scalar_flux.copy()             # (nx, ng)

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
    r"""2-D wavefront sweep — per-octant batched (Wave 2 / C2.6).

    Iterates over angular **octants** (lexicographic order from
    :attr:`AngularQuadrature.octants`). For each in-plane octant
    :math:`\sigma = (\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y)`:

    1. **BC apply once** on the octant-incoming face(s) — replaces
       the per-ordinate ``bc.apply(...)[n]`` calls of the legacy
       implementation (saves ``N`` invocations per sweep).
    2. **Dispatch to ``SNMesh.sweep_graphs[OctantLabel(σ)]``** —
       the per-octant ``SweepDependencyGraph`` precomputed once at
       mesh construction (Wave 2 / C2.4).
    3. The graph's ``apply`` walks topological levels (anti-diagonals)
       and dispatches each level to
       ``sn_mesh.cell_update.update_batch(slice_args)`` — vectorised
       over ``(N_oct, n_diag, ng)`` simultaneously.

    The smoking gun ``for n in range(N)`` is gone: the outer loop is
    now ``for octant in quad.octants`` (4-8 iterations, structural),
    and within each octant the work is one BC apply + one
    ``SweepDependencyGraph.apply`` call. The ordinate axis is
    INTERNAL to every numpy operation.

    Bit-identity to legacy
    ----------------------

    For LS-family quadratures whose ordinate ordering is
    octant-grouped in lexicographic order (``LevelSymmetricSN``,
    ``ProductQuadrature``), this implementation is bit-identical to
    the legacy per-ordinate loop on every snapshot — verified by
    ``tests/sn/test_2d_octant_sweep_equivalence.py`` (the C2.5
    TESTS-FIRST harness). The argument:

    * BC apply for octant σ reads partners (in the x-reflected
      octant). For LS4, lex order matches legacy n-order at the
      partner-state granularity, so the same iteration's value is
      observed at the same point.
    * Within an octant, the per-ordinate cell sweeps are independent
      (different rows of ``psi_x`` / ``psi_y``); batching is
      bit-identical to any per-ordinate sequencing of the same set.

    For quadratures whose ordering is NOT lex (e.g. Lebedev), the
    converged answer is the same Gauss-Seidel fixed point but the
    iter-to-iter values may differ. Per ``vv-principles`` Bit-identity
    vs principled-equivalence: principled at every step (octant-batched
    BC apply IS the §15A.2 form), structurally-independent grounding
    (closed-form L1 anchor + MMS regression suite both pass),
    convergence to the same answer (verified empirically on regression
    snapshots).
    """
    nx, ny, ng = Q.shape
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights

    angular_flux = np.zeros((N, nx, ny, ng))
    scalar_flux = np.zeros((nx, ny, ng))

    # Persistent boundary-flux buffers for reflective BCs (mutated in
    # place across sweep calls — partner reads need the previous
    # iteration's outgoing-face writes).
    if "bc_2d_x" not in psi_bc:
        psi_bc["bc_2d_x"] = np.zeros((N, nx + 1, ny, ng))
        psi_bc["bc_2d_y"] = np.zeros((N, nx, ny + 1, ng))
    psi_x = psi_bc["bc_2d_x"]   # (N, nx+1, ny, ng)
    psi_y = psi_bc["bc_2d_y"]   # (N, nx, ny+1, ng)

    # Source pre-scale once (was inside the ordinate loop in legacy
    # for a no-op O(1) load — kept here for clarity).
    weight_norm = 1.0 / weights.sum()
    Q_scaled = Q * weight_norm
    has_aniso = Q_aniso is not None
    if has_aniso:
        Q_aniso_scaled = Q_aniso * weight_norm  # (N, nx, ny, ng)

    str_x = sn_mesh.streaming_x   # (N, nx)
    str_y = sn_mesh.streaming_y   # (N, ny)
    cell_update = sn_mesh.cell_update

    for octant in quad.octants:
        label_tuple = octant.label   # e.g. (-1, +1, +1) or (-1,) for 1-D
        oct_idx = octant.indices     # (N_oct,) int into N
        N_oct = oct_idx.size
        # In-plane signs — drop sign_z if the label is 3-D.
        sx = label_tuple[0] if len(label_tuple) >= 1 else +1
        sy = label_tuple[1] if len(label_tuple) >= 2 else 0
        # Pure-z degenerate octant: no in-plane streaming. The
        # angular flux is the volumetric balance ψ = Q_n / Σ_t and
        # the scalar flux gets a weighted contribution.
        if sx == 0 and sy == 0:
            Q_pure_z = Q_scaled[None, :, :, :]    # (1, nx, ny, ng)
            if has_aniso:
                Q_pure_z = Q_pure_z + Q_aniso_scaled[oct_idx]
            psi_avg_pure_z = Q_pure_z / sig_t     # broadcasts to (N_oct, nx, ny, ng)
            angular_flux[oct_idx] = psi_avg_pure_z
            scalar_flux += np.einsum(
                "nijg,n->ijg", psi_avg_pure_z, weights[oct_idx],
            )
            continue

        # Effective in-plane sign for sweep-graph lookup. Match
        # legacy's ``key = (1 if mx >= 0 else -1, ...)`` mapping:
        # ordinates with ``mx == 0`` are treated as ``+1`` (the BC
        # apply uses xmin, the streaming coefficient is zero, and
        # the WDD result is identical regardless of sign choice).
        sx_eff = +1 if sx == 0 else sx
        sy_eff = +1 if sy == 0 else sy
        sweep_graph = sn_mesh.sweep_graphs[OctantLabel(sx_eff, sy_eff)]

        # ── BC apply once per octant ───────────────────────────────
        #
        # The octant's incoming-x face is index 0 if sx >= 0 else nx.
        # Boundary operators take the FULL (N, ny, ng) face buffer
        # (so reflective BCs can read partner-octant rows) and return
        # an updated full buffer; we scatter only this octant's rows
        # back into psi_x / psi_y.
        if sx_eff >= 0:
            full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, 0, :, :])
            psi_x[oct_idx, 0, :, :] = full_face_x[oct_idx]
        else:
            full_face_x = sn_mesh.bc_xmax.apply(psi_x[:, nx, :, :])
            psi_x[oct_idx, nx, :, :] = full_face_x[oct_idx]

        if sy_eff >= 0:
            full_face_y = sn_mesh.bc_ymin.apply(psi_y[:, :, 0, :])
            psi_y[oct_idx, :, 0, :] = full_face_y[oct_idx]
        else:
            full_face_y = sn_mesh.bc_ymax.apply(psi_y[:, :, ny, :])
            psi_y[oct_idx, :, ny, :] = full_face_y[oct_idx]

        # ── Per-octant buffers for the graph apply ────────────────
        #
        # Fancy indexing creates copies (not views), so we extract
        # mutable per-octant buffers, run the graph, and scatter back.
        # The cost (~2 × N_oct × nx × ny × ng × 8 bytes per octant)
        # is negligible vs the loop savings.
        psi_x_oct = psi_x[oct_idx].copy()    # (N_oct, nx+1, ny, ng)
        psi_y_oct = psi_y[oct_idx].copy()    # (N_oct, nx, ny+1, ng)
        Q_octant = Q_scaled[None, :, :, :]   # (1, nx, ny, ng) broadcastable
        if has_aniso:
            Q_octant = Q_octant + Q_aniso_scaled[oct_idx]   # (N_oct, nx, ny, ng)
        angular_flux_oct = np.zeros((N_oct, nx, ny, ng))

        sweep_graph.apply(
            cell_update=cell_update,
            psi_x_octant=psi_x_oct,
            psi_y_octant=psi_y_oct,
            Q_octant=Q_octant,
            sig_t=sig_t,
            str_x_octant=str_x[oct_idx],
            str_y_octant=str_y[oct_idx],
            weights_octant=weights[oct_idx],
            angular_flux_octant=angular_flux_oct,
            scalar_flux_buf=scalar_flux,
        )

        # Scatter back the persistent BC buffers + angular flux.
        psi_x[oct_idx] = psi_x_oct
        psi_y[oct_idx] = psi_y_oct
        angular_flux[oct_idx] = angular_flux_oct

    return angular_flux, scalar_flux
