r"""Production (numpy) d-generic UBLD primitive + the shared d=1 closed form.

Branch-2 of sub-step **D5b-S1** (#240 / #38 / #37) — the numpy production
counterpart to the SymPy algebra-of-record in
:mod:`orpheus.derivations.discrete.sn.ld_ubld`.  Two layers, one source of
truth:

1. **The d-generic dense primitive** (:func:`assemble_ubld` /
   :func:`per_cell_solve`).  Assembles the :math:`2^d \times 2^d` Kronecker
   cell system ``A = G + F_out + Σ_t·M`` and ``R = M·S + Σ_axes F_in`` for
   ``d = len(hs)`` axes, BATCHED over a stack of cells / anti-diagonal slots,
   and solves it with a batched :func:`numpy.linalg.solve`.  This MIRRORS the
   symbolic assembler one-for-one (the 1-D Legendre-basis factors
   ``M_1d = diag(h, θh)``, ``G_1d = |μ|[[0,0],[-2,0]]``,
   ``F_out_1d = |μ|[[1,1],[1,1]]``, ``F_in`` trace ``B(-1) = [1,-1]·|μ|``); NO
   :math:`2^d \times 2^d` entry is hand-transcribed.  It is the CANONICAL
   source for BOTH ``d = 1`` (today, the production slab LD reduces from it)
   and ``d ≥ 2`` (S2 will wire the bilinear cell-batch kernel onto it — NOT
   done here).  In production, ``d = 1`` does NOT route through this dense
   solve (that would be the L16 per-cell-solve regression); the dense primitive
   is the d-generic *reference* the closed-form fast path is proven equal to.

2. **The shared d=1 closed form** (:func:`d1_closed_form`).  The analytic
   Schur complement of the primitive's ``d = 1`` ``2×2`` — VECTORIZED over the
   cell / ordinate / group stack (no per-cell dense solve), so the production
   ``d = 1`` SCAN stays on the fast path (L16).  This is the ONE place the slab
   LD ``2×2`` SCAN/per-cell algebra lives: the two ×V production views
   (:meth:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous._schur_terms`
   — ×V per-cell ``update``/``residual``;
   :meth:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.affine_scan_coefficients`
   — ×V ``CumprodScan`` coefficients) derive their coefficients from this
   helper, applying their ×V scaling at the call site.  The DAG-wavefront kernel
   pair (``cell_kernel_batch`` / ``residual_kernel_batch``) does NOT use the
   closed form — a matvec is a forward APPLY, intrinsically moment-valued, so it
   rides the d-generic dense primitive (1) for every d (#240 D5b-S3, the unified
   moment matvec); at d=1 the dense ``2×2`` IS this closed form's Schur, proven
   ``==`` (symbolically by
   :mod:`orpheus.derivations.discrete.sn.ld_ubld`, in code by
   ``tests/transport/spatial/test_ld_ubld_primitive.py``).

The scale-free invariants
=========================

The d=1 Schur is naturally parameterised by the **÷V streaming-over-volume**
``g = |μ| A_down / V`` and the local :math:`\Sigma_t`.  Two dimensionless
quantities carry the entire closure and are SCALE-INVARIANT (identical in the
×V and ÷V conventions):

.. math::

   k = \frac{g/\theta}{g/\theta + \Sigma_t}, \qquad w = \frac{1}{1 + k}.

``w`` is the cell-average blend weight (:math:`\bar\psi = (1-w)\psi_{\rm in} +
w\,\psi_{\rm out}`); ``k`` is the slope-elimination ratio.  Every view's
coefficients are an algebraic function of ``(g, Σ_t, k, w)`` times a power of
the cell volume ``V`` (the ×V vs ÷V choice).  See :class:`D1ClosedForm`.

References
==========

* :mod:`orpheus.derivations.discrete.sn.ld_ubld` — the SymPy algebra-of-record
  this module is the numpy Branch-2 of (the symbolic ``simplify(diff) == 0``
  proofs of the d=1 reduction + the three-view equality).
* :mod:`orpheus.transport.spatial.linear_discontinuous` — the production LD scheme
  whose two ×V 1-D views (per-cell + scan) single-source through
  :func:`d1_closed_form`; its DAG kernel rides the d-generic dense primitive.
* Larsen & Morel (1989). JCP 83(1):212-236 — the slab-LD moment system.
* Maginot, Ragusa & Morel (2016). NSE 185(1):17-42 — the d-D Cartesian UBLD.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ═══════════════════════════════════════════════════════════════════════
# The 1-D LD factor operators (Legendre moment basis {1, P₁} on width h)
# ═══════════════════════════════════════════════════════════════════════
#
# The numpy mirror of the symbolic factors in
# ``orpheus.derivations.discrete.sn.ld_ubld``.  ``mass_1d`` scales with the
# cell width ``h``; the streaming / face factors are dimensionless ``×|μ|``.

#: dimensionless 1-D streaming/stiffness ``∫B_i ∂_x B_j`` (per unit ``|μ|``).
_GRAD_1D = np.array([[0.0, 0.0], [-2.0, 0.0]])
#: dimensionless 1-D downstream-face ``outer(B(+1), B(+1))`` (per unit ``|μ|``).
_FOUT_1D = np.array([[1.0, 1.0], [1.0, 1.0]])
#: dimensionless 1-D upstream-face test-weighting ``B(-1)`` (per unit ``|μ|``).
_FIN_TRACE = np.array([1.0, -1.0])


def octant_moment_frame_signs(
    octant_signs: tuple[int, ...], spatial_basis_per_axis: int,
) -> "np.ndarray | None":
    r"""Sweep-frame ⇄ global-frame sign vector for the ``2^d`` moment axis.

    The per-cell LD moment vector is produced and consumed by the cell kernel
    in the per-ordinate SWEEP frame: each axis ``a`` is oriented so the
    "downstream" face is at the local coordinate ``+1``.  For an ordinate
    sweeping in the NEGATIVE global direction on axis ``a`` (``octant_signs[a]
    == -1``) the sweep coordinate is the REVERSE of the global coordinate, so
    the slope (``P₁``) moment on that axis is sign-FLIPPED relative to the
    global-frame slope (#240 D5b-S3 — the diffusion-limit root cause).

    The iterate / scattering source ``Σ_s·φ̂`` lives in the GLOBAL frame (the
    angular reduction ``φ̂ = Σ_n w_n ψ̂_n`` sums slopes across ordinates of BOTH
    sweep directions — they must share one frame or they partially cancel).
    This returns the ``2^d``-length sign vector that maps between the two
    frames in the tensor-Legendre Kronecker layout ``[o_0, …, o_{d-1}]``
    (``o_a ∈ {0,1}`` = P₀/P₁ on axis ``a``, axis 0 outer):

    .. math::

        \mathrm{sign}[o_0,\dots,o_{d-1}] = \prod_{a=0}^{d-1}
            (\mathrm{octant\_signs}[a])^{o_a}.

    The average moment (all ``o_a = 0``) is sign-invariant; a per-axis slope
    flips once if that axis sweeps backward; the d=2 cross moment ``x̂y`` flips
    when an ODD number of its active axes sweep backward.  The map is an
    INVOLUTION (its own inverse), so the SAME vector converts global→sweep on
    the source/probe INPUT and sweep→global on the moment/residual OUTPUT.

    Single-moment closures (DD/Step, ``spatial_basis_per_axis == 1``) have a
    sign-invariant average-only moment, so this returns ``None`` — the single
    source of the "no frame map" convention every consumer keys on (``_reframe``
    passes the array through untouched, so DD/Step stay byte-identical; the
    backward-compat invariant).  A zero octant sign (the degenerate no-streaming
    pure-transverse set) leaves the slope on that axis un-flipped (``0^0 = 1``
    for the average; a slope on a non-streaming axis is a degenerate ordinate
    the wavefront short-circuits — never reached here).
    """
    if spatial_basis_per_axis == 1:
        return None
    d = len(octant_signs)
    size = spatial_basis_per_axis ** d
    signs = np.ones(size)
    for flat in range(size):
        factor = 1
        for a in range(d):
            o_a = (flat // (spatial_basis_per_axis ** (d - 1 - a))) % spatial_basis_per_axis
            if o_a and octant_signs[a] < 0:
                factor = -factor
        signs[flat] = factor
    return signs


def mass_1d(h: np.ndarray, theta: float) -> np.ndarray:
    r"""Batched 1-D LD mass ``diag(h, θh)`` — shape ``(..., 2, 2)``.

    The numpy mirror of :func:`orpheus.derivations.discrete.sn.ld_ubld.mass_1d`.
    ``h`` is broadcast over a leading batch (cells / anti-diagonal slots).
    """
    h = np.asarray(h, dtype=np.float64)
    out = np.zeros(h.shape + (2, 2), dtype=np.float64)
    out[..., 0, 0] = h
    out[..., 1, 1] = theta * h
    return out


def _batched_kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    r"""Kronecker product over the trailing two axes, batched over the lead.

    For ``a`` of shape ``(..., p, q)`` and ``b`` of shape ``(..., r, s)`` (same
    leading batch), returns ``(..., p·r, q·s)`` — the per-batch Kronecker
    product.  The numpy batched analogue of :func:`sympy.kronecker_product`.
    """
    *batch, p, q = a.shape
    r, s = b.shape[-2:]
    # (..., p, 1, q, 1) * (..., 1, r, 1, s) -> (..., p, r, q, s) -> (..., p·r, q·s)
    prod = a[..., :, None, :, None] * b[..., None, :, None, :]
    return prod.reshape(*batch, p * r, q * s)


# ═══════════════════════════════════════════════════════════════════════
# The d-generic UBLD assembler (Kronecker products — single source)
# ═══════════════════════════════════════════════════════════════════════


def assemble_ubld(
    hs: list[np.ndarray],
    mus: list[np.ndarray],
    sig_t: np.ndarray,
    theta: float,
) -> dict:
    r"""Assemble the d-generic UBLD cell matrices via batched Kronecker products.

    The numpy Branch-2 mirror of
    :func:`orpheus.derivations.discrete.sn.ld_ubld.assemble_ubld`.  Builds the
    mass ``M``, streaming ``G``, downstream-face ``F_out``, and
    ``A = G + F_out + Σ_t·M`` for ``d = len(hs)`` axes, BATCHED over a leading
    stack (cells / anti-diagonal slots / ordinates / groups), so a single call
    assembles the whole level's :math:`2^d \times 2^d` systems.

    The streaming and downstream-face terms put the gradient / outflow on the
    ACTIVE axis and the **mass** on every transverse axis (the volume /
    transverse-face integral factorisation), exactly as the symbolic
    assembler — this is the SINGLE site of the UBLD algebra; ``d = 1, 2, 3``
    all flow through it.

    Parameters
    ----------
    hs :
        Per-axis cell widths, ``d`` arrays each broadcastable to the batch.
    mus :
        Per-axis :math:`|\mu_a|` (already sweep-pre-resolved, non-negative),
        ``d`` arrays each broadcastable to the batch.
    sig_t :
        Local :math:`\Sigma_t`, broadcastable to the batch.
    theta :
        The slope-moment weight :math:`\theta` (SN-exact ``1/3``).

    Returns
    -------
    dict
        ``{"A", "M", "G", "F_out", "size", "d"}`` — ``A`` / ``M`` / ``G`` /
        ``F_out`` are ``(..., 2^d, 2^d)`` batched matrices.
    """
    d = len(hs)
    if len(mus) != d:
        raise ValueError(f"hs and mus must agree on d; got {d} and {len(mus)}.")
    masses = [mass_1d(np.asarray(h, dtype=np.float64), theta) for h in hs]

    # Mass = Kronecker product of all 1-D masses.
    M = masses[0]
    for k in range(1, d):
        M = _batched_kron(M, masses[k])

    size = 2**d
    batch_shape = M.shape[:-2]
    G = np.zeros(batch_shape + (size, size), dtype=np.float64)
    F_out = np.zeros(batch_shape + (size, size), dtype=np.float64)
    for a in range(d):
        mu_a = np.asarray(mus[a], dtype=np.float64)
        # The active axis carries the (dimensionless) gradient / outflow factor,
        # broadcast across the batch; the transverse axes carry the mass.
        grad_factors = [
            np.broadcast_to(_GRAD_1D, batch_shape + (2, 2)) if k == a else masses[k]
            for k in range(d)
        ]
        face_factors = [
            np.broadcast_to(_FOUT_1D, batch_shape + (2, 2)) if k == a else masses[k]
            for k in range(d)
        ]
        Ga = grad_factors[0]
        Foa = face_factors[0]
        for k in range(1, d):
            Ga = _batched_kron(Ga, grad_factors[k])
            Foa = _batched_kron(Foa, face_factors[k])
        G += mu_a[..., None, None] * Ga
        F_out += mu_a[..., None, None] * Foa

    A = G + F_out + np.asarray(sig_t, dtype=np.float64)[..., None, None] * M
    return {"A": A, "M": M, "G": G, "F_out": F_out, "size": size, "d": d}


def assemble_inflow_axis(
    hs: list[np.ndarray],
    mus: list[np.ndarray],
    axis: int,
    upstream_face_moments: np.ndarray,
    theta: float,
) -> np.ndarray:
    r"""Inflow (upwind) RHS contribution from one upstream face — batched.

    The numpy mirror of
    :func:`orpheus.derivations.discrete.sn.ld_ubld.assemble_inflow_axis`.  The
    incoming face normal to ``axis`` is the upstream neighbour's outflow trace
    — a ``2^{d-1}``-moment transverse object (``upstream_face_moments``, shape
    ``(..., 2^{d-1})`` in the transverse-axes Kronecker order) — weighted into
    the active-axis test functions ``B(-1) = [1,-1]`` and the **mass** on the
    transverse axes, times ``|μ_axis|``.  Returns ``(..., 2^d)``.

    Supports ``axis ∈ {0, d-1}`` (the boundary slots — all the d≤2 cases need);
    the general interior-axis interleave is deferred (S2, the d≥3 cell kernel).
    """
    d = len(hs)
    masses = [mass_1d(np.asarray(h, dtype=np.float64), theta) for h in hs]
    # Transverse mass over the (2^{d-1}) upstream-face moment vector.
    transverse_mass = None
    for k in range(d):
        if k == axis:
            continue
        transverse_mass = (
            masses[k]
            if transverse_mass is None
            else _batched_kron(transverse_mass, masses[k])
        )
    face = np.asarray(upstream_face_moments, dtype=np.float64)
    weighted_face = (
        face
        if transverse_mass is None
        else np.einsum("...ij,...j->...i", transverse_mass, face)
    )
    mu_axis = np.asarray(mus[axis], dtype=np.float64)
    trace = _FIN_TRACE  # B(-1) = [1, -1]
    batch_shape = weighted_face.shape[:-1]
    if axis == 0:
        # kron(trace, weighted_face): outer = active axis, inner = transverse.
        out = trace[:, None] * weighted_face[..., None, :]
        out = out.reshape(*batch_shape, trace.size * weighted_face.shape[-1])
    elif axis == d - 1:
        out = weighted_face[..., :, None] * trace[None, :]
        out = out.reshape(*batch_shape, weighted_face.shape[-1] * trace.size)
    else:
        raise NotImplementedError(
            "assemble_inflow_axis supports axis in {0, d-1}; the general "
            f"interior-axis interleave (d={d}, axis={axis}) is deferred (S2)."
        )
    return mu_axis[..., None] * out


def per_cell_solve(assembled: dict, rhs: np.ndarray) -> np.ndarray:
    r"""Solve the batched ``A ψ⃗ = R`` for the ``2^d``-moment cell unknowns.

    ``rhs`` has shape ``(..., 2^d)`` matching the assembled ``A`` batch; returns
    ``(..., 2^d)``.  Uses a batched :func:`numpy.linalg.solve` (the d-generic
    dense reference — NOT the production d=1 fast path).
    """
    A = assembled["A"]
    return np.linalg.solve(A, rhs[..., None])[..., 0]


# ═══════════════════════════════════════════════════════════════════════
# The shared d=1 closed form — the ONE place the slab LD 2×2 algebra lives
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class D1ClosedForm:
    r"""The vectorized analytic Schur of the d=1 UBLD ``2×2`` (the fast path).

    The single source of the slab-LD ``2×2`` algebra.  Built by
    :func:`d1_closed_form` from the ÷V streaming-over-volume ``g = |μ|A_down/V``
    and the local :math:`\Sigma_t`; carries the SCALE-FREE invariants
    (``g_over_theta``, ``d2``, ``k``, ``w``) plus the ÷V Schur denominator
    ``eff_denom``.  The three production views derive their coefficients from
    these by applying the ×V / ÷V scaling at the call site:

    * **÷V (DAG kernel)** — ``eff_denom`` is the ÷V Schur ``S``; the RHS folds
      ``g·(1 + g_over_theta/d2)·ψ_in`` and the source ``Q`` (per-volume); the
      outflow rides ``w``.  See :meth:`kernel_rhs`.
    * **×V (per-cell / scan)** — multiply by the cell volume ``V``:
      ``S_×V = V·eff_denom``; the scan transmission ``a`` and the per-cell
      Schur intermediates follow.  See :meth:`schur_xV` / :meth:`scan_xV`.

    All arrays broadcast over the ``(N, ng, nx)`` (or ``(N_oct, ng, n_diag)``)
    stack — there is NO per-cell Python loop and NO per-cell dense solve, so
    the production d=1 path stays on the fast closed-form path (L16).
    """

    g: np.ndarray             #: ``|μ|·A_down/V`` — ÷V streaming-over-volume
    g_over_theta: np.ndarray  #: ``g/θ`` — the slope-row streaming
    d2: np.ndarray            #: ``g/θ + Σ_t`` — ÷V slope denominator (D₂/V)
    k: np.ndarray             #: ``(g/θ)/d2`` — slope-elimination ratio
    w: np.ndarray             #: ``1/(1+k)`` — cell-average blend weight
    eff_denom: np.ndarray     #: ``(g+Σ_t) + g·k`` — the ÷V Schur diagonal S/V
    theta: float              #: the slope-moment weight θ (SN-exact 1/3)

    def kernel_rhs(self, Q_cells: np.ndarray, psi_in: np.ndarray) -> np.ndarray:
        r"""The ÷V (DAG-kernel) RHS: ``Q + g·ψ_in + g·(g/θ)·ψ_in/d2``.

        ``Q_cells`` is the per-volume average-moment source (flat, ``Q̂ = 0``);
        ``psi_in`` the upstream face flux.  The ÷V Schur balance is then
        ``eff_denom · ψ̄ = kernel_rhs``.
        """
        return Q_cells + self.g * psi_in + self.g * self.g_over_theta * psi_in / self.d2

    def _xV(self, V: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        r"""The ×V streaming ``m = |μ|A_down = g·V`` and Schur diagonal ``S = V·eff_denom``.

        The shared ÷V→×V scaling both ×V views (:meth:`schur_xV`, :meth:`scan_xV`)
        build on — single-sourced so the convention lives in ONE place (it goes
        live the moment ÷V→×V changes, e.g. a curvilinear ``V_eff`` or the D6
        ``w(Σ)`` blend).
        """
        return self.g * V, V * self.eff_denom

    def _slope_fold(
        self, V: np.ndarray | float, s_hat: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        r"""The shared slope-row fold — the ONE site of the ``s_hat`` algebra.

        Single-sources the slope-moment elimination both ×V slope consumers ride
        — the per-cell Schur (:meth:`schur_xV`, the matvec/per-cell update) AND
        the moment-aware SCAN (:meth:`scan_slope_face_source` / :meth:`scan_reconstruct`,
        the 1-D production sweep, #240 D5b-S3 OWED-2).  With ``|μ|A_down = g·V``
        and the ×V slope denominator :math:`D_2' = \theta V d_2 = \Sigma_t\theta V
        + |\mu|A_{\rm down}` (the extra ``θ`` vs the ÷V ``d2``):

        * ``mu_Adown = |μ|A_down`` — the ×V streaming (slope reconstruction),
        * ``d2p = D₂'`` — the ×V slope denominator,
        * ``eff_source_shift = s_hat·θ·|μ|A_down/D₂'`` — the amount the slope
          source pulls OUT of the average-row effective source
          (``eff_source = s_bar − eff_source_shift``),
        * ``slope_source = θ·s_hat`` — the slope-row RHS source (drives ``ψ̂``).
        """
        mu_Adown = self.g * V
        d2p = self.theta * V * self.d2
        eff_source_shift = s_hat * mu_Adown * self.theta / d2p
        slope_source = self.theta * s_hat
        return mu_Adown, d2p, eff_source_shift, slope_source

    def schur_xV(
        self,
        V: np.ndarray | float,
        s_bar: np.ndarray,
        s_hat: np.ndarray,
        psi_in: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        r"""The ×V per-cell Schur intermediates (the ``_LDCellTerms`` fields).

        Reproduces the production
        :class:`~orpheus.transport.spatial.linear_discontinuous._LDCellTerms` fields in
        the ×V contract (``source`` = ``Q·V`` moments).  With ``|μ|A_down = g·V``
        the ×V streaming and the ×V slope denominator

        .. math::

           D_2' = \theta\,V\,d_2 = \Sigma_t\,\theta\,V + |\mu|\,A_{\rm down}

        (note the extra ``θ`` relative to the ÷V ``d2``: the production
        ``D₂'`` carries the slope-moment weight, ``d2`` is the ÷V slope denom):

        * ``S = V·eff_denom`` — the ×V Schur diagonal,
        * ``eff_source = s_bar − s_hat·θ·|μ|A_down/D₂'`` — the source folded
          through the slope row,
        * ``eff_numer = |μ|A_down·ψ_in·(D₂' + |μ|A_down)/D₂'`` — the inflow
          folded through the slope row,
        * ``slope_source = θ·s_hat`` — the slope-row RHS source,
        * ``mu_Adown = |μ|A_down`` — the ×V streaming (for the slope reconstr.),
        * ``D₂'`` — the ×V slope denominator.

        The slope-row fold is single-sourced through :meth:`_slope_fold` (shared
        with the moment-aware scan, #240 D5b-S3 OWED-2 — Pattern 2: ONE site for
        the ``s_hat`` algebra).  Returns ``(S, eff_source, eff_numer, slope_source,
        mu_Adown, D₂')`` so the production ``_schur_terms`` builds its bundle from
        this ONE site.
        """
        S = V * self.eff_denom                       # ×V Schur diagonal
        mu_Adown, d2p, eff_source_shift, slope_source = self._slope_fold(V, s_hat)
        eff_source = s_bar - eff_source_shift
        eff_numer = mu_Adown * psi_in * (d2p + mu_Adown) / d2p
        return S, eff_source, eff_numer, slope_source, mu_Adown, d2p

    def scan_xV(
        self, V: np.ndarray | float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""The ×V scan coefficients ``(a, inverse_denom, w)`` (transmission form).

        Reproduces
        :meth:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.affine_scan_coefficients`:
        with ``m = |μ|A_down = g·V``, ``S_×V = V·eff_denom``,

        .. math::

           \mathrm{inverse\_denom} = 1/S_{\times V}, \qquad
           a = m\,(1+k)^2 / S_{\times V} - k, \qquad w = 1/(1+k).

        ``a`` is source-independent (the affine transmission ``ψ_out = a·ψ_in +
        b``); ``w`` is the cell-average blend weight (scale-free).  These three
        coefficients are SOURCE-INDEPENDENT (one ``Σ_t``-epoch — the cache holds
        them); the SLOPE source ``s_hat`` enters the source-dependent ``b`` and
        the per-cell reconstruction via :meth:`scan_moment_terms` /
        :meth:`scan_reconstruct` (#240 D5b-S3 OWED-2).
        """
        m, S = self._xV(V)                          # |μ|·A_down, ×V Schur diagonal
        inverse_denom = 1.0 / S
        a = m * (1.0 + self.k) ** 2 * inverse_denom - self.k
        return a, inverse_denom, self.w

    def scan_slope_face_source(
        self, V: np.ndarray | float, s_hat: np.ndarray,
        inverse_denom: np.ndarray, w: np.ndarray | float,
    ) -> np.ndarray:
        r"""The SLOPE source's contribution to the face-chain affine source ``b``.

        The 1-D moment SCAN propagates the downstream FACE
        :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` along the chain (#240 D5b-S3
        OWED-2).  With a slope source ``s_hat`` (the scattering-slope source
        :math:`\Sigma_s\hat\phi` × V) the face-chain ``b`` gains a slope term on
        top of the flat-source emission ``source_emission(s_bar, inv, w)``:

        .. math::

           b = \underbrace{s_{\rm bar}\,\tfrac{\mathrm{inv}}{w}}
                          _{\text{flat-source (average) emission}}
             + \underbrace{\tfrac{\theta\,s_{\rm hat}}{D_2'}
               - \tfrac{\theta\,|\mu|A_{\rm down}\,s_{\rm hat}}{D_2'}
                 \cdot \tfrac{\mathrm{inv}}{w}}_{\text{this method (slope)}}.

        The first slope term ``θ s_hat / D₂'`` is the slope's direct lift onto the
        outgoing face (``ψ_out = ψ̄ + ψ̂``); the second is the slope source's pull
        on the average-row effective source (``−eff_source_shift`` run through the
        emission ``inv/w``).  The flat-source ``b`` stays the existing
        :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`
        call (DD/Step never reach here — ``s_hat`` is the LD-only slope moment),
        so this is a pure ADDITION, never a re-baseline of the average path.  The
        fold is single-sourced through :meth:`_slope_fold`.
        """
        _mu_Adown, d2p, eff_source_shift, slope_source = self._slope_fold(V, s_hat)
        return slope_source / d2p - eff_source_shift * inverse_denom / w

    def scan_reconstruct(
        self,
        V: np.ndarray | float,
        s_bar: np.ndarray,
        s_hat: np.ndarray,
        psi_in: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Per-cell ``(ψ̄, ψ̂)`` from the chained upstream face — the scan twin.

        Once the moment SCAN has the chained upstream face ``psi_in`` per cell
        (from :meth:`scan_slope_face_source` + the Blelloch scan), the
        cell-average ``ψ̄`` and slope ``ψ̂`` moments are the SAME Schur
        reconstruction the per-cell :meth:`schur_xV` does — re-used here so the
        moment iterate ``φ̂ = Σ_n w_n ψ̂_n`` the scan emits is byte-consistent with
        the DAG (#240 D5b-S3 OWED-2; both ride :meth:`_slope_fold`).  ``s_bar`` /
        ``s_hat`` are the ×V source moments (``Q̄·V`` / ``Q̂·V``); ``psi_in`` the
        per-cell upstream face.
        """
        S, eff_source, eff_numer, slope_source, mu_Adown, d2p = self.schur_xV(
            V, s_bar, s_hat, psi_in,
        )
        psi_bar = (eff_source + eff_numer) / S
        psi_hat = (mu_Adown * psi_bar + slope_source - mu_Adown * psi_in) / d2p
        return psi_bar, psi_hat


def d1_closed_form(
    g: np.ndarray | float, sig_t: np.ndarray | float, theta: float
) -> D1ClosedForm:
    r"""Build the shared d=1 LD closed form from ``g = |μ|A_down/V`` and ``Σ_t``.

    The ÷V streaming-over-volume ``g`` and the local :math:`\Sigma_t` fully
    determine the scale-free closure (``k``, ``w``) and the ÷V Schur diagonal
    ``eff_denom``.  This is the analytic Schur complement of the primitive's
    ``d = 1`` ``2×2`` — vectorized (no dense solve), so it is the production
    fast path (L16).  Proven ``== `` the dense :func:`per_cell_solve` at
    ``d = 1`` by ``tests/transport/spatial/test_ld_ubld_primitive.py``.
    """
    g = np.asarray(g, dtype=np.float64)
    sig_t = np.asarray(sig_t, dtype=np.float64)
    g_over_theta = g / theta
    d2 = g_over_theta + sig_t
    k = g_over_theta / d2
    w = 1.0 / (1.0 + k)
    eff_denom = (g + sig_t) + g * k
    return D1ClosedForm(
        g=g, g_over_theta=g_over_theta, d2=d2, k=k, w=w, eff_denom=eff_denom,
        theta=theta,
    )


__all__ = [
    "D1ClosedForm",
    "assemble_inflow_axis",
    "assemble_ubld",
    "d1_closed_form",
    "mass_1d",
    "octant_moment_frame_signs",
    "per_cell_solve",
]
