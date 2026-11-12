r"""Symbolic d-generic UBLD — the tensor-product (bi/tri)linear LD cell system.

Branch-1 (SymPy) algebra-of-record for the dimension-generic **Unlumped
BiLinear Discontinuous** (UBLD) per-cell Galerkin system of the SN
equations on a Cartesian cell.  This is sub-step **D5b-S1 Branch 1** of
issues #240 / #38 / #37 — the canonical symbolic algebra + the proofs
that de-risk the downstream numpy production primitive (Branch 2, a
SEPARATE follow-on).  It SHIPS NOTHING to production.

State **1A** (closed-form SymPy, ``algebra-of-record`` skill): the per-cell
system is a finite-dimensional dense linear solve in elementary symbolic
entries, so every claim closes as ``sympy.simplify(diff) == 0``.

The design fork (read the literature memo first)
================================================

The multi-D analog of LD on a *Cartesian* cell is **NOT** the simplex-P1
``{1, x, y}`` object (``1+d`` moments).  Adams 2001 proved simplex-LD
**fails** the thick diffusion limit on quadrilaterals, while the BILINEAR
DG-P1 (UBLD) — basis ``{1, x, y, xy}`` (``2^d`` moments) — **passes**.
The ``xy`` cross moment is diffusion-limit-load-bearing.  This module
builds the tensor-product bilinear/trilinear object, ``2^d`` moments per
(group, ordinate).  See
``.claude/agent-memory/literature-researcher/multi_d_ld_closure.md``
(MRM-2016 Eqs. 1-12; Adams-2001; BLA-1992).

The Kronecker build (Pattern 5 / Pattern 2 — single source of truth)
====================================================================

The streaming operator ``Ω·∇ = Σ_a μ_a ∂_a`` is a SUM over axes; each
term differentiates one axis and is the identity (in the L2 / mass sense)
on the others.  The tensor-product basis separates, so the three Galerkin
matrices factor as **Kronecker products of the verified 1-D LD factor
operators**:

* **Mass** ``M = M₁ ⊗ M₂ ⊗ …`` — the σ_t collision term.  Each 1-D factor
  is ``M_1d = diag(h, θh)`` (the Legendre-basis 1-D LD mass).
* **Streaming** ``G = Σ_a μ_a · (M₁ ⊗ … ⊗ G_1d ⊗ … ⊗ M_d)`` — the
  gradient/stiffness on the active axis ``a``, **mass on every transverse
  axis** (because the volume integral ``∫ B_i (∂_a B_j)`` factors as the
  axis-a stiffness times the transverse L2 inner products).  The 1-D
  factor is ``G_1d = |μ| · [[0, 0], [-2, 0]]`` (dimensionless ``∫ B_i
  ∂_x B_j`` in the Legendre basis on ``[-1,1]``-mapped width ``h``).
* **Downstream-face** ``F_out = Σ_a μ_a · (M₁ ⊗ … ⊗ F_out_1d ⊗ … ⊗ M_d)``
  — the DG-upwind outflow trace on the active axis, mass on the
  transverse axes (the transverse face integral).  The 1-D factor is
  ``F_out_1d = |μ| · outer(B(+1), B(+1)) = |μ| · [[1,1],[1,1]]`` (the
  downstream node ``B(+1) = [1, 1]`` in the moment basis).
* **Upstream-face** (inflow → RHS) ``F_in`` weights the upwind neighbour's
  outflow trace into the test functions ``B_i(-1) = [1, -1]`` on the
  active axis, mass on the transverse axes (the ``2^{d-1}``-moment face
  object: average + transverse slopes).

The assembled per-cell system (MRM-2016 Eq. 12) is

.. math::

   A_{\rm cell}\,\vec\psi = \vec R, \qquad
   A_{\rm cell} = G + F_{\rm out} + \Sigma_t M, \qquad
   \vec R = M\,\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces},

a ``2^d × 2^d`` dense non-symmetric solve.  ``d=1`` is the
"Kronecker-with-one-factor" identity — it reduces EXACTLY to the
production slab 2×2 (``orpheus.transport.spatial.linear_discontinuous``), the
``xy`` coupling falls out of the algebra for ``d≥2``, and NO 4×4/8×8
entry is hand-transcribed.

Verification functions
======================

* :func:`derive_d1_reduction_to_production` — Oracle (i): the assembled
  ``d=1`` ``A`` (the code binding of :math:`A_{\rm cell}`) / ``R`` / Schur ``S`` / slope ``ψ̂`` / ``D₂'`` reduce
  symbolically to the production ``_LDCellTerms`` / ``_schur_terms``
  closed forms.
* :func:`derive_d1_kernel_view_equals` — Oracle (i), ÷V view: the
  production ``_kernel_terms`` (``g, g/θ, D₂, eff_denom, w``) flat-source
  form equals the ``d=1`` reduction.
* :func:`derive_d1_scan_view_equals` — Oracle (i), ×V view: the production
  ``affine_scan_coefficients`` (``m, t, p, D₂, k, S, a, w``) form equals
  the ``d=1`` reduction.  Together these two prove "single-source the
  math" is achievable for Branch 2.
* :func:`derive_d2_exact_on_bilinear` — Oracle (ii): the ``d=2`` UBLD
  recovers ANY bilinear flux ``ψ = a + bx + cy + dxy`` exactly (the
  multi-D analog of the 1-D "exact on linear-in-x" oracle); the ``xy``
  coupling is exercised (``d ≠ 0``).
* :func:`derive_d1_transpose_equals_At_Minv` — the #310 C2 R1b keystone:
  the LD cell VJP is ``Aᵀ M⁻¹`` (mass-inverse FIRST — the order
  discriminant is nonzero), with the inflow-fold and face-trace pullback
  rows derived from the exact symbolic Jacobians.
* :func:`derive_octant_frame_sign_is_involution` — the moment-frame sign
  map ``D_s`` is an involution and frame conjugation commutes with
  transpose (the reverse walk's reframe-in/reframe-out law).
* :func:`derive_d1_integrated_slope_source_is_the_galerkin_moment` — the
  data side (#707 Q-D): the integrated moment source is the cell mass on
  the Legendre coefficients, ``M_cell·q_coef`` (slab ``V·diag(1, θ)``), equal
  to the production scan's ×V pre-scale times the closure's ×θ; the V-only
  reading is the negative.
* :func:`derive_d1_coordinate_pairing_is_the_integral` — in that chart the
  coordinate pairing of a state with stored data is ``∫ψq`` exactly.
* :func:`derive_curvilinear_ld_cell_gram` — forward-looking: the sphere and
  cylinder LD cell Gram is full and position-dependent.
* :func:`derive_d3_assembles` — structural readiness: the ``d=3``
  trilinear assembler produces an 8×8 system with the ``θ³`` triple-cross
  diagonal weight.

References
==========

* Maginot, Ragusa & Morel (2016). *Non-negative Methods for Bilinear
  Discontinuous Differencing of the S_N Equations on Quadrilaterals.* NSE
  185(1):17-42.  §2 Eqs. (1)-(12) — the 2-D Cartesian UBLD weak form.
* Adams (2001). *Discontinuous Finite Element Transport Solutions in
  Thick Diffusive Problems.* NSE 137(3):298-333 — simplex-LD FAILS /
  bilinear PASSES the thick diffusion limit on quadrilaterals.
* Börgers, Larsen & Adams (1992). *The asymptotic diffusion limit of a
  linear discontinuous discretization of a two-dimensional linear
  transport equation.* JCP 98(2):285-300.
* Larsen & Morel (1989). JCP 83(1):212-236 — the 1-D slab-LD moment
  system (the d=1 reduction oracle, θ=1/3).

See also
========

* :mod:`orpheus.transport.spatial.linear_discontinuous` — the production 1-D LD
  this module's d=1 reduction is proven equal to.
* :mod:`orpheus.derivations.discrete.sn.balance` — the sibling discrete-SN
  symbolic derivation (DD / WDD / cumprod).
"""

from __future__ import annotations

import sympy as sp

# ═══════════════════════════════════════════════════════════════════════
# The 1-D LD factor operators (Legendre moment basis {1, P₁} on width h)
# ═══════════════════════════════════════════════════════════════════════
#
# These four factors are the verified building blocks (reverse-engineered
# from the production slab 2×2 and confirmed below by the d=1 reduction
# oracle).  The d-generic assembler composes them via Kronecker products.

THETA = sp.Symbol("theta", positive=True)
r"""The slope-moment weight :math:`\theta` (LM-1989 Eq. 4.3b).  The numeric
SN-exact value is :math:`1/3`; kept symbolic so the algebra is general."""


def mass_1d(h: sp.Expr, theta: sp.Expr = THETA) -> sp.Matrix:
    r"""1-D LD mass matrix in the Legendre basis: ``diag(h, θh)``.

    The diagonal of ``∫_cell B_i B_j`` for the L2-orthogonal moment basis
    ``{P₀=1, P₁}``: ``⟨P₀,P₀⟩ = h``, ``⟨P₁,P₁⟩ = θh`` (the ``θ=1/3``
    normalisation IS this 1-D mass diagonal ratio).
    """
    return sp.diag(h, theta * h)


def grad_1d() -> sp.Matrix:
    r"""1-D streaming gradient/stiffness ``∫ B_i ∂_x B_j`` (per unit ``|μ|``).

    Dimensionless ``[[0, 0], [-2, 0]]`` in the Legendre basis: only
    ``∫ P₁ ∂_x P₀`` survives the IBP combination as the volumetric
    coupling, and the ``-2`` is the mapped derivative ``∂_x P₁ = 2/h``
    times the basis inner product (the ``h`` cancels against the face/mass
    scaling — see the d=1 reduction proof).
    """
    return sp.Matrix([[0, 0], [-2, 0]])


def fout_1d() -> sp.Matrix:
    r"""1-D downstream-face matrix ``outer(B(+1), B(+1))`` (per unit ``|μ|``).

    The DG-upwind OUTFLOW trace at the downstream node: ``B(+1) = [1, 1]``
    in the moment basis (``P₀(+1)=1``, ``P₁(+1)=1``), so the implicit
    self-coupling face term is the rank-1 ``[[1,1],[1,1]]``.
    """
    return sp.Matrix([[1, 1], [1, 1]])


def fin_trace_weight() -> sp.Matrix:
    r"""1-D upstream-face test-weighting ``B(-1) = [1, -1]`` (per unit ``|μ|``).

    The INFLOW (upwind) face weights the upstream neighbour's outflow
    trace into the test functions evaluated at the upstream node
    ``B(-1) = [1, -1]`` (``P₀(-1)=1``, ``P₁(-1)=-1``).  Goes to the RHS
    ``F_in`` (upwind / discontinuous — no continuity enforced).
    """
    return sp.Matrix([1, -1])


# ═══════════════════════════════════════════════════════════════════════
# The d-generic UBLD assembler (Kronecker products — single source)
# ═══════════════════════════════════════════════════════════════════════


def _kron(*factors: sp.Matrix) -> sp.Matrix:
    """The Kronecker product of ``factors`` in axis order (axis 0 outermost)."""
    out = factors[0]
    for factor in factors[1:]:
        out = sp.Matrix(sp.kronecker_product(out, factor))
    return out


def assemble_ubld(
    hs: list[sp.Expr],
    mus: list[sp.Expr],
    sig_t: sp.Expr,
    theta: sp.Expr = THETA,
) -> dict:
    r"""Assemble the d-generic UBLD per-cell matrices via Kronecker products.

    Builds ``M`` (mass), ``G`` (streaming), ``F_out`` (downstream face),
    and ``A = G + F_out + Σ_t M`` for ``d = len(hs)`` axes.  Axis order is
    the Kronecker order ``hs[0]`` outer … ``hs[-1]`` inner; the moment
    vector index ``Σ_a 2^{d-1-a}·o_a`` (``o_a ∈ {0,1}`` the moment order
    on axis ``a``) — e.g. d=2 moment order ``[ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]``.

    The streaming and downstream-face terms put the gradient / outflow on
    the ACTIVE axis and the **mass** on every transverse axis (the
    volume / transverse-face integral factorisation).  This is the SINGLE
    site of the UBLD algebra — d=1, d=2, d=3 all flow through it.

    Returns
    -------
    dict
        ``{"A", "M", "G", "F_out", "size", "d"}`` — the assembled
        ``2^d × 2^d`` SymPy matrices and the dimension metadata.
    """
    d = len(hs)
    if len(mus) != d:
        raise ValueError(f"hs and mus must agree on d; got {d} and {len(mus)}.")
    masses = [mass_1d(h, theta) for h in hs]

    # Mass = Kronecker product of all 1-D masses.
    M = _kron(*masses)

    size = 2**d
    G = sp.zeros(size, size)
    F_out = sp.zeros(size, size)
    for a in range(d):
        # Active axis a carries the gradient / outflow; transverse axes the mass.
        G += _kron(*[grad_1d() if k == a else masses[k] for k in range(d)]) * mus[a]
        F_out += _kron(*[fout_1d() if k == a else masses[k] for k in range(d)]) * mus[a]

    A = G + F_out + sig_t * M
    return {"A": A, "M": M, "G": G, "F_out": F_out, "size": size, "d": d}


def inflow_scatter_axis(
    hs: list[sp.Expr],
    mus: list[sp.Expr],
    axis: int,
    theta: sp.Expr = THETA,
) -> sp.Matrix:
    r"""The inflow lift of one upstream face: ``2^{d-1}`` face moments to the cell RHS.

    The incoming face normal to ``axis`` carries the upstream neighbour's
    outflow trace, a ``2^{d-1}``-moment transverse object ordered as the
    Kronecker moment vector of the transverse axes (in increasing axis
    order).  It is weighted into the test functions at the upstream node,
    ``B(-1) = [1, -1]``, on ``axis``, and by the 1-D mass on every
    transverse axis, times ``|μ_axis|``: the face integral
    ``|μ_a| ∮ B_i ψ_upwind`` factorises exactly as the volume terms of
    :func:`assemble_ubld` do.  Returns the ``2^d × 2^{d-1}`` matrix.

    The ``B(-1)`` column sits in the ``axis`` slot of the Kronecker product,
    so every axis, interior ones included (``0 < axis < d-1`` at ``d ≥ 3``),
    is the same one-line product.
    """
    d = len(hs)
    masses = [mass_1d(h, theta) for h in hs]
    factors = [fin_trace_weight() if k == axis else masses[k] for k in range(d)]
    # ``Matrix * scalar`` (scalar on the right) stays Matrix-typed; the
    # ``scalar * Matrix`` order returns an Expr in the SymPy stubs.
    return _kron(*factors) * mus[axis]


def assemble_inflow_axis(
    hs: list[sp.Expr],
    mus: list[sp.Expr],
    axis: int,
    upstream_face_moments: sp.Matrix,
    theta: sp.Expr = THETA,
) -> sp.Matrix:
    r"""Inflow (upwind) RHS contribution from one upstream face.

    :func:`inflow_scatter_axis` applied to the upstream face moments: a
    ``2^d``-vector RHS contribution.
    """
    return inflow_scatter_axis(hs, mus, axis, theta) * upstream_face_moments


def per_cell_solve(assembled: dict, rhs: sp.Matrix) -> sp.Matrix:
    r"""Solve ``A ψ⃗ = R`` for the ``2^d``-moment cell unknown vector."""
    return assembled["A"].LUsolve(rhs)


def outflow_trace_axis(d: int, axis: int) -> sp.Matrix:
    r"""The outflow trace on the downstream face normal to ``axis``.

    The cell function evaluated at the downstream node, ``B(+1) = [1, 1]``,
    on ``axis``, with every transverse moment kept: the ``2^{d-1} × 2^d``
    matrix from the cell moments to the face moments, in the transverse
    Kronecker order :func:`inflow_scatter_axis` consumes.  For ``d = 1`` it is
    the LM-1989 Eq. 4.3c closure ``ψ_out = ψ̄ + ψ̂``.
    """
    factors = [sp.Matrix([[1, 1]]) if k == axis else sp.eye(2) for k in range(d)]
    return _kron(*factors)


def downstream_face_trace(psi_moments: sp.Matrix, d: int) -> sp.Expr:
    r"""The d=1 outgoing face flux ``B(+1)·ψ⃗ = ψ̄ + ψ̂``, a scalar.

    For ``d ≥ 2`` the face is a ``2^{d-1}``-moment object:
    :func:`outflow_trace_axis`.
    """
    if d != 1:
        raise ValueError(
            f"downstream_face_trace is the d=1 scalar trace; got d={d}. "
            "Use outflow_trace_axis for the 2^(d-1)-moment face."
        )
    return (outflow_trace_axis(1, 0) * psi_moments)[0]


# ═══════════════════════════════════════════════════════════════════════
# Oracle (i) — the d=1 reduction to the production slab LD
# ═══════════════════════════════════════════════════════════════════════


def derive_d1_reduction_to_production() -> dict:
    r"""V_d1 — the assembled d=1 UBLD reduces to the production slab 2×2.

    Proves: the Kronecker-with-one-factor assembly reproduces the
    production ``A`` / ``R`` / Schur ``S`` / slope ``ψ̂`` / ``D₂'`` closed
    forms (``orpheus.transport.spatial.linear_discontinuous._LDCellTerms`` /
    ``_schur_terms``) symbolically.
    """
    mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
    Qbar, Qhat, psi_in = sp.symbols("Qbar Qhat psi_in", real=True)
    theta = THETA

    asm = assemble_ubld([h], [mu], sig_t, theta)
    A = asm["A"]
    # Inflow on the single axis 0: upstream face is a scalar trace ψ_in;
    # its "moment" vector is the 1-D scalar [ψ_in] -> use axis-0 inflow.
    R = asm["M"] * sp.Matrix([Qbar, Qhat]) + mu * fin_trace_weight() * psi_in

    # Production natural mass-weighted 2×2 (module docstring lines 52-64).
    A_prod = sp.Matrix([[sig_t * h + mu, mu], [-mu, sig_t * theta * h + mu]])
    R_prod = sp.Matrix([Qbar * h + mu * psi_in, theta * Qhat * h - mu * psi_in])
    diff_A = sp.simplify(A - A_prod)
    diff_R = sp.simplify(R - R_prod)

    # Per-cell solve + outgoing-face closure ψ_out = ψ̄ + ψ̂ = B(+1)·ψ⃗.
    psi = per_cell_solve(asm, R)
    psi_bar = sp.simplify(psi[0])
    psi_hat = sp.simplify(psi[1])
    psi_out = sp.simplify(psi_bar + psi_hat)
    diff_face = sp.simplify(psi_out - downstream_face_trace(psi, d=1))

    # Production Schur closed forms (linear_discontinuous.py lines 332-335,
    # 254-259).  s_bar = Qbar·h (= source[0]), s_hat = Qhat·h (= source[1]).
    s_bar, s_hat = Qbar * h, Qhat * h
    D2p = sig_t * theta * h + mu
    S = (mu**2 + (sig_t * h + mu) * D2p) / D2p
    eff_source = s_bar - s_hat * mu * theta / D2p
    eff_numer = mu * psi_in * (D2p + mu) / D2p
    psi_bar_prod = (eff_source + eff_numer) / S
    # slope reconstruction (_LDCellTerms.slope): (μψ̄ + θ·s_hat − μψ_in)/D₂'
    psi_hat_prod = (mu * psi_bar_prod + theta * s_hat - mu * psi_in) / D2p

    diff_bar = sp.simplify(psi_bar - psi_bar_prod)
    diff_hat = sp.simplify(psi_hat - psi_hat_prod)

    passed = all(
        x == 0 if not hasattr(x, "is_zero_matrix") else x.is_zero_matrix
        for x in (diff_A, diff_R)
    ) and all(x == 0 for x in (diff_face, diff_bar, diff_hat))

    return {
        "name": "V_d1: assembled d=1 UBLD == production slab LD 2×2 + Schur S/ψ̂/D₂'",
        "diff_A": diff_A,
        "diff_R": diff_R,
        "diff_face": diff_face,
        "diff_psi_bar": diff_bar,
        "diff_psi_hat": diff_hat,
        "S": sp.simplify(S),
        "D2_prime": D2p,
        "pass": passed,
    }


def derive_d1_kernel_view_equals() -> dict:
    r"""V_d1_kernel — the production ÷V ``_kernel_terms`` form equals d=1.

    Proves: the ÷V kernel coefficients (``g = |μ|/Δ``, ``g_over_theta``,
    ``d2``, ``eff_denom``, ``w = 1/(1+g_over_theta/d2)``) and the generic
    base reconstruction (``ψ_out = (ψ̄ − (1−w)ψ_in)/w``) reproduce the d=1
    assembled ``ψ̄`` / ``ψ_out`` at the FLAT source slice (``Q̂ = 0``).
    Half of the "single-source the math" proof for Branch 2.
    """
    mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
    Qbar, psi_in = sp.symbols("Qbar psi_in", real=True)
    theta = THETA

    asm = assemble_ubld([h], [mu], sig_t, theta)
    R = asm["M"] * sp.Matrix([Qbar, 0]) + mu * fin_trace_weight() * psi_in
    psi = per_cell_solve(asm, R)
    psi_bar = sp.simplify(psi[0])
    psi_out = sp.simplify(psi[0] + psi[1])

    # production ÷V kernel (_kernel_terms, lines 443-453); Q is per-volume Qbar.
    g = mu / h
    g_over_theta = g / theta
    d2 = g_over_theta + sig_t
    eff_denom = (g + sig_t) + g * g_over_theta / d2
    rhs_k = Qbar + g * psi_in + g * g_over_theta * psi_in / d2
    w = 1 / (1 + g_over_theta / d2)
    psi_avg_k = sp.simplify(rhs_k / eff_denom)
    # outgoing_face_from_average: ψ_out = (ψ̄ − (1−w)ψ_in)/w
    psi_out_k = sp.simplify((psi_avg_k - (1 - w) * psi_in) / w)

    diff_bar = sp.simplify(psi_avg_k - psi_bar)
    diff_out = sp.simplify(psi_out_k - psi_out)
    return {
        "name": "V_d1_kernel: production ÷V _kernel_terms == d=1 UBLD (flat Q̂=0)",
        "diff_psi_bar": diff_bar,
        "diff_psi_out": diff_out,
        "pass": diff_bar == 0 and diff_out == 0,
    }


def derive_d1_scan_view_equals() -> dict:
    r"""V_d1_scan — the production ×V ``affine_scan_coefficients`` form equals d=1.

    Proves: the ×V scan coefficients (``m = |μ|face_area_downstream``, ``t = Σ_t V``,
    ``p = m/θ``, ``D₂ = t+p``, ``k = p/D₂``, ``S = (t+m)+m·p/D₂``,
    ``a = m(1+k)²/S − k``, ``w = 1/(1+k)``) with the generic
    ``source_emission`` / ``cell_average`` reconstruction reproduce the
    d=1 assembled ``ψ̄`` / ``ψ_out`` (flat ``Q̂ = 0``).  The transmission
    ``a`` is source-independent (the affine recurrence ``ψ_out = aψ_in +
    b``).  The other half of "single-source the math" for Branch 2.
    """
    mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
    Qbar, psi_in = sp.symbols("Qbar psi_in", real=True)
    theta = THETA

    asm = assemble_ubld([h], [mu], sig_t, theta)
    R = asm["M"] * sp.Matrix([Qbar, 0]) + mu * fin_trace_weight() * psi_in
    psi = per_cell_solve(asm, R)
    psi_bar = sp.simplify(psi[0])
    psi_out = sp.simplify(psi[0] + psi[1])

    # production ×V scan (affine_scan_coefficients, lines 564-571); slab face_area_downstream=1, V=h.
    m = mu
    t = sig_t * h
    p = m / theta
    D2 = t + p
    k = p / D2
    S = (t + m) + m * p / D2
    inv = 1 / S
    a = m * (1 + k) ** 2 / S - k
    w = 1 / (1 + k)
    QV = Qbar * h
    # source_emission(QV, inv, w) = QV·inv/w; ψ_out = a·ψ_in + b
    b = QV * inv / w
    psi_out_s = a * psi_in + b
    # cell_average(ψ_in, ψ_out, w) = (1−w)ψ_in + w·ψ_out
    psi_avg_s = (1 - w) * psi_in + w * psi_out_s

    diff_bar = sp.simplify(psi_avg_s - psi_bar)
    diff_out = sp.simplify(psi_out_s - psi_out)
    a_source_independent = Qbar not in a.free_symbols
    return {
        "name": "V_d1_scan: production ×V affine_scan_coefficients == d=1 UBLD (flat Q̂=0)",
        "diff_psi_bar": diff_bar,
        "diff_psi_out": diff_out,
        "a_source_independent": a_source_independent,
        "pass": diff_bar == 0 and diff_out == 0 and a_source_independent,
    }


# ═══════════════════════════════════════════════════════════════════════
# Oracle (ii) — exact on a bilinear flux (the xy coupling is exercised)
# ═══════════════════════════════════════════════════════════════════════


def derive_d2_exact_on_bilinear() -> dict:
    r"""V_d2_bilinear — the d=2 UBLD recovers ψ = a + bx + cy + dxy exactly.

    The multi-D analog of the 1-D "exact on linear-in-x" oracle.  Feeds
    the manufactured cell traces (upstream-x and upstream-y DG-exact face
    moments) and source moments for a bilinear ψ, and proves the solved
    ``2^2 = 4`` tensor-Legendre moments equal the EXACT projections of ψ.
    The ``xy`` cross moment is exercised (``d ≠ 0`` symbolic).
    """
    theta = sp.Rational(1, 3)  # SN-exact value; closed-form integration
    a, b, c, d = sp.symbols("a b c d", real=True)
    mu_x, mu_y, sig_t = sp.symbols("mu_x mu_y Sigma_t", positive=True)
    hx, hy = sp.symbols("h_x h_y", positive=True)
    x0, y0 = sp.symbols("x0 y0", real=True)
    xx, yy = sp.symbols("x y", real=True)

    xmx, xmy = x0 + hx / 2, y0 + hy / 2  # cell centres

    # EXACT tensor-Legendre moments of ψ = a + bx + cy + dxy on the cell,
    # basis {1, L₁} with L₁(x)=2(x−xc)/h, L2-normalised (so the moment IS
    # the coefficient).  Moment order [bar, ŷ, x̂, x̂y] (kron x-outer,y-inner).
    pbar = a + b * xmx + c * xmy + d * xmx * xmy
    phat_x = b * hx / 2 + d * (hx / 2) * xmy
    phat_y = c * hy / 2 + d * xmx * (hy / 2)
    phat_xy = d * (hx / 2) * (hy / 2)
    psi_moments_exact = sp.Matrix([pbar, phat_y, phat_x, phat_xy])

    # Manufactured source Q = Ω·∇ψ + Σ_t ψ.
    Qfield = (
        mu_x * (b + d * yy)
        + mu_y * (c + d * xx)
        + sig_t * (a + b * xx + c * yy + d * xx * yy)
    )

    def moment(expr, ox, oy):
        Lx = 2 * (xx - xmx) / hx
        Ly = 2 * (yy - xmy) / hy
        basis = (Lx**ox) * (Ly**oy)
        norm = sp.Rational(1, 1)
        if ox == 1:
            norm *= sp.Rational(1, 3)
        if oy == 1:
            norm *= sp.Rational(1, 3)
        integ = sp.integrate(
            sp.integrate(expr * basis, (xx, x0, x0 + hx)), (yy, y0, y0 + hy)
        ) / (hx * hy)
        return sp.simplify(integ / norm)

    Sbar = moment(Qfield, 0, 0)
    Shx = moment(Qfield, 1, 0)
    Shy = moment(Qfield, 0, 1)
    Shxy = moment(Qfield, 1, 1)
    S_moments = sp.Matrix([Sbar, Shy, Shx, Shxy])  # [bar, ŷ, x̂, x̂y]

    asm = assemble_ubld([hx, hy], [mu_x, mu_y], sig_t, theta)

    # DG-exact upstream-x face trace = ψ(x0, y); its (bar_y, slope_y) moments.
    def y_face_moments(expr_y):
        Ly = 2 * (yy - xmy) / hy
        m0 = sp.integrate(expr_y, (yy, y0, y0 + hy)) / hy
        m1 = sp.integrate(expr_y * Ly, (yy, y0, y0 + hy)) / hy / sp.Rational(1, 3)
        return sp.Matrix([m0, m1])

    def x_face_moments(expr_x):
        Lx = 2 * (xx - xmx) / hx
        m0 = sp.integrate(expr_x, (xx, x0, x0 + hx)) / hx
        m1 = sp.integrate(expr_x * Lx, (xx, x0, x0 + hx)) / hx / sp.Rational(1, 3)
        return sp.Matrix([m0, m1])

    trace_x = a + b * x0 + c * yy + d * x0 * yy  # ψ(x0, y)
    trace_y = a + b * xx + c * y0 + d * xx * y0  # ψ(x, y0)
    Fin_x = assemble_inflow_axis([hx, hy], [mu_x, mu_y], 0, y_face_moments(trace_x), theta)
    Fin_y = assemble_inflow_axis([hx, hy], [mu_x, mu_y], 1, x_face_moments(trace_y), theta)

    R = asm["M"] * S_moments + Fin_x + Fin_y
    psi_solved = per_cell_solve(asm, R)
    diff = sp.simplify(psi_solved - psi_moments_exact)
    return {
        "name": "V_d2_bilinear: d=2 UBLD recovers a+bx+cy+dxy exactly (xy coupling exercised)",
        "diff": diff,
        "pass": diff.is_zero_matrix,
    }


def derive_d3_assembles() -> dict:
    r"""V_d3 — the d=3 trilinear UBLD assembles to an 8×8 with θ³ cross weight.

    Structural readiness for d=3 (no full bilinear-exactness oracle here —
    that is Branch-2 numerical territory).  Confirms the d-generic
    assembler produces the 8×8 trilinear system and the ``xyz`` triple-cross
    moment carries the ``Σ_t h_x h_y h_z θ³`` diagonal weight.
    """
    mu_x, mu_y, mu_z, sig_t = sp.symbols("mu_x mu_y mu_z Sigma_t", positive=True)
    hx, hy, hz = sp.symbols("h_x h_y h_z", positive=True)
    theta = THETA
    asm = assemble_ubld([hx, hy, hz], [mu_x, mu_y, mu_z], sig_t, theta)
    A = asm["A"]
    size_ok = A.shape == (8, 8)
    # The xyz moment is index 7 (all axes order-1); its collision diagonal weight.
    xyz_collision = sp.simplify(asm["M"][7, 7] * sig_t / (hx * hy * hz))
    weight_ok = sp.simplify(xyz_collision - sig_t * theta**3) == 0
    return {
        "name": "V_d3: d=3 trilinear UBLD assembles 8×8 with θ³ triple-cross weight",
        "size": A.shape,
        "xyz_collision_weight": sp.simplify(asm["M"][7, 7] / (hx * hy * hz)),
        "pass": size_ok and weight_ok,
    }


def derive_d1_transpose_equals_At_Minv() -> dict:
    r"""V_d1_T — the LD cell VJP IS ``Aᵀ M⁻¹`` (mass-inverse FIRST) + the pullbacks.

    The production ÷V residual relation
    (:meth:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.residual_kernel_batch`)
    is ``res = M⁻¹(A ψ⃗ − R)`` with ``R = M·S⃗ + μ·B(−1)·ψ_in`` and the
    outgoing face the downstream trace ``ψ_out = B(+1)·ψ⃗``.  This oracle
    builds the d=1 forward map SYMBOLICALLY, takes its exact Jacobians, and
    proves the reverse-mode pullbacks are what the #310 C2 kernel spells
    (`.claude/plans/archive/residue_verification_spec.md` §3.1 — the R1b keystone):

    * ``∂res/∂ψ⃗ = M⁻¹A``  ⟹  the ψ⃗-cotangent matrix is
      ``(M⁻¹A)ᵀ = Aᵀ M⁻¹`` (``M`` diagonal ⟹ ``M⁻ᵀ = M⁻¹``) — the mass
      inverse applies BEFORE ``Aᵀ`` (NEW-algebra (a); the
      M-R1b-MASSORDER mutation target);
    * the ORDER DISCRIMINANT: ``Aᵀ M⁻¹ ≠ M⁻¹ Aᵀ`` at generic symbols —
      had the two commuted, the mass-order gate would be toothless;
    * ``∂res/∂ψ_in = −μ·M⁻¹·B(−1)``  ⟹  the upstream-face pullback is
      ``−μ·B(−1)ᵀ·M⁻¹·r̄`` (the inflow-fold transpose);
    * ``∂ψ_out/∂ψ⃗ = B(+1) = [1, 1]``  ⟹  the face cotangent BROADCASTS
      equally into both moment slots (the trace transpose).
    """
    mu, h, sig_t = sp.symbols("mu h Sigma_t", positive=True)
    theta = THETA
    asm = assemble_ubld([h], [mu], sig_t, theta)
    A, M = asm["A"], asm["M"]
    Minv = M.inv()

    psi_syms = sp.Matrix(sp.symbols("psibar psihat", real=True))
    psi_in = sp.Symbol("psi_in", real=True)
    Sbar, Shat = sp.symbols("Sbar Shat", real=True)

    # The forward map, built honestly from the assembled record.
    R = M * sp.Matrix([Sbar, Shat]) + fin_trace_weight() * psi_in * mu
    res = Minv * (A * psi_syms - R)

    # (i) ψ⃗-pullback: jacobian(res, ψ⃗)ᵀ == Aᵀ·M⁻¹ (the kernel's order).
    jac_psi = res.jacobian(psi_syms)
    vjp_matrix = A.T * Minv
    diff_mass_order = sp.simplify(jac_psi.T - vjp_matrix)

    # (ii) the order discriminant — the wrong composition differs.
    discriminant = sp.simplify(vjp_matrix - Minv * A.T)

    # (iii) ψ_in-pullback: ∂res/∂ψ_in transposed == −μ·B(−1)ᵀ·M⁻¹.
    jac_in = res.diff(psi_in)                                # (2, 1)
    vjp_in_row = -(fin_trace_weight().T * Minv) * mu         # (1, 2)
    diff_inflow = sp.simplify(jac_in.T - vjp_in_row)

    # (iv) face-trace pullback: the ψ_out Jacobian row is B(+1) = [1, 1],
    # so the face cotangent broadcasts into BOTH moment slots.
    face = downstream_face_trace(psi_syms, d=1)
    face_jac = sp.Matrix([[sp.diff(face, s) for s in psi_syms]])
    diff_face = sp.simplify(face_jac.T - sp.Matrix([1, 1]))

    passed = (
        diff_mass_order.is_zero_matrix
        and not discriminant.is_zero_matrix
        and diff_inflow.is_zero_matrix
        and diff_face.is_zero_matrix
    )
    return {
        "name": "V_d1_T: LD cell VJP == Aᵀ·M⁻¹ (mass-inverse first) + pullback rows",
        "diff_mass_order": diff_mass_order,
        "order_discriminant_nonzero": not discriminant.is_zero_matrix,
        "diff_inflow_pullback": diff_inflow,
        "diff_face_pullback": diff_face,
        "pass": passed,
    }


def derive_interval_legendre_gram(n: int = 5) -> dict:
    r"""V_gram — the interval Legendre basis: Gram, endpoint values, parity, and θ.

    The algebra of record of
    :class:`orpheus.numerics.basis.interval_legendre_basis.IntervalLegendreBasis`
    (the moment-basis carve, #423). For :math:`P_k` on :math:`\xi \in [-1, 1]`
    (SymPy's ``legendre``, :math:`P_k(1) = 1`), against the normalised measure
    :math:`\mathrm d\xi/2`, proves for :math:`k, l < n`:

    * the Gram :math:`\int_{-1}^{1} P_k P_l\,\mathrm d\xi/2 = \delta_{kl}/(2k+1)`
      (DLMF Eq. 18.3, Table 18.3.1, with the measure halved);
    * the endpoint values :math:`P_k(\pm 1) = (\pm 1)^k`;
    * the parity :math:`P_k(-\xi) = (-1)^k P_k(\xi)`;
    * the slope-moment weight :math:`\theta = 1/3` (LM-1989 Eq. 4.3b) is the
      degree-1 Gram, and the degree-0 Gram is 1 (why the d = 1 closed forms
      are parameterised by θ alone).
    """
    xi = sp.Symbol("xi", real=True)
    P = [sp.legendre(k, xi) for k in range(n)]
    gram = sp.Matrix(n, n, lambda k, l: sp.integrate(P[k] * P[l], (xi, -1, 1)) / 2)
    gram_ok = sp.simplify(gram - sp.diag(*[sp.Rational(1, 2 * k + 1) for k in range(n)])).is_zero_matrix
    endpoints_ok = all(P[k].subs(xi, 1) == 1 and P[k].subs(xi, -1) == (-1) ** k for k in range(n))
    parity_ok = all(sp.expand(P[k].subs(xi, -xi) - (-1) ** k * P[k]) == 0 for k in range(n))
    theta_ok = gram[0, 0] == 1 and gram[1, 1] == sp.Rational(1, 3)
    return {
        "name": "V_gram: interval Legendre Gram, endpoint values, parity, theta",
        "gram": gram,
        "gram_ok": gram_ok,
        "endpoints_ok": endpoints_ok,
        "parity_ok": parity_ok,
        "theta_ok": theta_ok,
        "pass": bool(gram_ok and endpoints_ok and parity_ok and theta_ok),
    }


def derive_octant_frame_sign_is_involution() -> dict:
    r"""V_frame_T — the octant moment-frame sign map: involution, transpose-preserved.

    The sweep⇄global frame map is the diagonal ``D_s`` with
    ``D_s[o⃗, o⃗] = ∏_a signs[a]^{o_a}`` (the symbolic mirror of
    :func:`orpheus.transport.spatial._ubld.octant_moment_frame_signs`).
    Proves, at d=1 (backward sweep) and d=2 (backward-x; backward-both):

    * ``D_s·D_s = I`` — the map is an INVOLUTION, so ONE vector converts
      global→sweep on the input AND sweep→global on the output;
    * ``(D_s·A·D_s)ᵀ = D_s·Aᵀ·D_s`` — frame conjugation COMMUTES with
      transpose, so the reverse walk reframes the residual cotangent IN
      and the ψ⃗-cotangent OUT with the SAME vector the forward uses (the
      adjoint of ``x ↦ D·K(D·x)`` is ``y ↦ D·Kᵀ(D·y)``) — NEW-algebra
      (c); a one-sided or forgotten ``D_s`` breaks the identity
      (M-R1b-FRAMESIGN).
    """
    checks: dict[str, bool] = {}
    all_pass = True
    for d, octant_signs in ((1, (-1,)), (2, (-1, 1)), (2, (-1, -1))):
        hs = [sp.Symbol(f"h_{a}", positive=True) for a in range(d)]
        mus = [sp.Symbol(f"mu_{a}", positive=True) for a in range(d)]
        sig_t = sp.Symbol("Sigma_t", positive=True)
        asm = assemble_ubld(hs, mus, sig_t, THETA)
        size = asm["size"]
        # The ∏ signs[a]^{o_a} diagonal in the Kronecker moment order —
        # the same index arithmetic as the numpy octant_moment_frame_signs.
        signs = []
        for flat in range(size):
            factor = 1
            for a in range(d):
                o_a = (flat >> (d - 1 - a)) & 1
                if o_a and octant_signs[a] < 0:
                    factor = -factor
            signs.append(factor)
        D = sp.diag(*signs)
        involution = sp.simplify(D * D - sp.eye(size))
        conjugation = sp.simplify((D * asm["A"] * D).T - D * asm["A"].T * D)
        label = f"d={d},signs={octant_signs}"
        ok = involution.is_zero_matrix and conjugation.is_zero_matrix
        checks[label] = ok
        all_pass = all_pass and ok
    return {
        "name": "V_frame_T: D_s involution + frame conjugation commutes with transpose",
        "checks": checks,
        "pass": all_pass,
    }


# ═══════════════════════════════════════════════════════════════════════
# The data side (#707 Q-D): the integrated chart of the moment rows
# ═══════════════════════════════════════════════════════════════════════
#
# The data side (sources, residuals, responses) is stored INTEGRATED
# (#707 R6) and dual on the spatial factor only, the in-cell moment axis
# included (R8, premise confirmed by the user 2026-11-11). These three
# derivations fix what "integrated" means for a moment row: the Galerkin
# moment q_int,i = ∫_cell q b_i dV = (M_cell·q_coef)_i, never V·q_coef.


def derive_d1_integrated_slope_source_is_the_galerkin_moment() -> dict:
    r"""V_qint — the integrated moment source is the cell mass on the Legendre coefficients.

    For a cell source :math:`q = \bar Q + \hat Q P_1(\xi)` in the trial span
    (the Legendre coefficient chart, :math:`P_1(+1) = 1`), proves:

    * :math:`\int_{\rm cell} q\,b_i\,\mathrm dx = (M_{\rm cell}\,q_{\rm coef})_i`
      with :math:`M_{\rm cell} =` :func:`mass_1d` ``(h)``
      :math:`= h\,\mathrm{diag}(1, \theta)`, :math:`\theta = 1/3` (slab), and
      on the d=2 cell :math:`M_{\rm cell} = \mathrm{mass}_x \otimes \mathrm{mass}_y`,
      the rows weighted :math:`V\theta^{|o|}` (the bilinear row :math:`V\theta^2`);
    * the slope row of the upwind weak form, solved for :math:`\hat\psi`, equals
      the production ``D1ClosedForm.scan_reconstruct`` slope line
      :math:`\hat\psi = (|\mu|\bar\psi + \theta s_{\rm hat} - |\mu|\psi_{\rm in}) / D_2'`
      with :math:`s_{\rm hat} = \hat Q V` (retyped from
      ``orpheus/transport/spatial/_ubld.py`` ``_slope_fold``: the ×V pre-scale
      of the scan plus the closure's ×θ IS the mass);
    * NEGATIVE: the V-only reading (slope RHS :math:`= V\hat Q`, θ dropped)
      leaves the residual :math:`-2\hat Q h/(\Sigma_t h + 3|\mu|)`, nonzero.
    """
    xi = sp.Symbol("xi", real=True)
    mu, h, hx, hy, sig_t = sp.symbols("mu h h_x h_y Sigma_t", positive=True)
    Qbar, Qhat, psi_bar, psi_hat, psi_in = sp.symbols("Qbar Qhat psibar psihat psi_in", real=True)
    theta = sp.Rational(1, 3)
    P = [sp.Integer(1), xi]  # the trial basis b_i = P_i(ξ); dx = h dξ/2

    # (1) slab: ∫ q b_i dx == mass_1d(h)·q_coef
    q = Qbar + Qhat * xi
    q_int = sp.Matrix([sp.integrate(q * b * h / 2, (xi, -1, 1)) for b in P])
    diff_mass = sp.simplify(q_int - mass_1d(h, theta) * sp.Matrix([Qbar, Qhat]))

    # (1') d=2: ∫∫ q b_i b_j dxdy == (mass_x ⊗ mass_y)·q_coef, Kronecker order (x outer)
    eta = sp.Symbol("eta", real=True)
    c = sp.symbols("c0:4", real=True)
    basis2 = [sp.Integer(1), eta, xi, xi * eta]  # [ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]
    q2 = sum(ci * bi for ci, bi in zip(c, basis2))
    q2_int = sp.Matrix([
        sp.integrate(sp.integrate(q2 * b * (hx / 2) * (hy / 2), (xi, -1, 1)), (eta, -1, 1))
        for b in basis2
    ])
    M2 = _kron(mass_1d(hx, theta), mass_1d(hy, theta))
    diff_mass_d2 = sp.simplify(q2_int - M2 * sp.Matrix(c))
    V2 = hx * hy
    row_weights_d2 = [sp.simplify(M2[i, i] / V2) for i in range(4)]
    weights_ok = row_weights_d2 == [1, theta, theta, theta**2]

    # (2) the slope row of the weak form: test P₁, IBP, upwind ψ_in, ψ_out = ψ̄ + ψ̂
    psi = psi_bar + psi_hat * xi
    psi_out = psi.subs(xi, 1)
    slope_row = (
        mu * (psi_out * P[1].subs(xi, 1) - psi_in * P[1].subs(xi, -1))
        - mu * sp.integrate(psi * sp.diff(P[1], xi) * (2 / h) * (h / 2), (xi, -1, 1))
        + sig_t * sp.integrate(psi * P[1] * h / 2, (xi, -1, 1))
        - q_int[1]
    )
    psi_hat_weak = sp.solve(slope_row, psi_hat)[0]
    D2p = sig_t * theta * h + mu
    s_hat = Qhat * h  # the scan's ×V pre-scale
    psi_hat_prod = (mu * psi_bar + theta * s_hat - mu * psi_in) / D2p
    diff_slope_row = sp.simplify(psi_hat_weak - psi_hat_prod)

    # (3) negative: the V-only reading drops the closure's θ
    psi_hat_V = (mu * psi_bar + s_hat - mu * psi_in) / D2p
    residual_V_only = sp.simplify(psi_hat_weak - psi_hat_V)
    expected_residual = -2 * Qhat * h / (sig_t * h + 3 * mu)
    negative_ok = residual_V_only != 0 and sp.simplify(residual_V_only - expected_residual) == 0

    passed = bool(
        diff_mass.is_zero_matrix and diff_mass_d2.is_zero_matrix and weights_ok
        and diff_slope_row == 0 and negative_ok
    )
    return {
        "name": "V_qint: integrated moment source == M_cell·q_coef (slab V·diag(1,θ); d=2 V·θ^|o|)",
        "diff_mass_slab": diff_mass,
        "diff_mass_d2": diff_mass_d2,
        "row_weights_d2_over_V": row_weights_d2,
        "diff_slope_row_vs_scan_reconstruct": diff_slope_row,
        "residual_V_only": residual_V_only,
        "pass": passed,
    }


def derive_d1_coordinate_pairing_is_the_integral() -> dict:
    r"""V_pair — in the V·θ chart the coordinate pairing IS the integral.

    For a trial-span state :math:`\psi = \bar\psi + \hat\psi P_1` and ANY source
    :math:`q` (here :math:`q = e^{a x}`, not in the trial span), stored
    integrated as :math:`q_{\rm int, i} = \int q\,b_i\,\mathrm dx`, proves

    .. math:: \psi_{\rm coef}\cdot q_{\rm int} = \int_{\rm cell}\psi\,q\,\mathrm dx

    exactly, and for a trial-span :math:`q` that :math:`q_{\rm int} = M q_{\rm coef}`
    so the pairing is :math:`\psi_{\rm coef}^{\mathsf T} M q_{\rm coef}`. NEGATIVE:
    the V chart (:math:`V q_{\rm coef}`) misses by
    :math:`(1-\theta) h\,\hat\psi\hat Q = 2h\hat\psi\hat Q/3`. This is the
    discriminating question, not the transpose reciprocity
    :math:`\langle A^{-1}q, r\rangle = \langle q, A^{-\mathsf T} r\rangle`, which
    holds for every invertible matrix in every chart.
    """
    xi = sp.Symbol("xi", real=True)
    h, a = sp.symbols("h a", positive=True)
    x0, psi_bar, psi_hat, Qbar, Qhat = sp.symbols("x0 psibar psihat Qbar Qhat", real=True)
    theta = sp.Rational(1, 3)
    P = [sp.Integer(1), xi]
    psi = psi_bar + psi_hat * xi
    x = x0 + xi * h / 2

    # any source: q = exp(a x) on the cell
    q = sp.exp(a * x)
    q_int = sp.Matrix([sp.integrate(q * b * h / 2, (xi, -1, 1)) for b in P])
    pairing = (sp.Matrix([psi_bar, psi_hat]).T * q_int)[0]
    integral = sp.integrate(psi * q * h / 2, (xi, -1, 1))
    diff_any_source = sp.simplify(pairing - integral)

    # trial-span source: q_int = M q_coef; V-chart negative
    qs = Qbar + Qhat * xi
    integral_s = sp.integrate(psi * qs * h / 2, (xi, -1, 1))
    pair_mass = (sp.Matrix([psi_bar, psi_hat]).T * mass_1d(h, theta) * sp.Matrix([Qbar, Qhat]))[0]
    pair_V = h * (psi_bar * Qbar + psi_hat * Qhat)
    diff_mass_chart = sp.simplify(pair_mass - integral_s)
    miss_V_chart = sp.simplify(pair_V - integral_s)
    negative_ok = sp.simplify(miss_V_chart - 2 * h * psi_hat * Qhat / 3) == 0

    passed = bool(diff_any_source == 0 and diff_mass_chart == 0 and negative_ok)
    return {
        "name": "V_pair: ψ_coef·q_int == ∫ψq (V·θ chart); V chart misses 2hψ̂Q̂/3",
        "diff_any_source": diff_any_source,
        "diff_mass_chart": diff_mass_chart,
        "miss_V_chart": miss_V_chart,
        "pass": passed,
    }


def derive_curvilinear_ld_cell_gram() -> dict:
    r"""V_gram_curv — the LD cell Gram under r²dr and r dr (FORWARD-LOOKING).

    LD refuses curvilinear meshes today (``LinearDiscontinuous.supports_curvilinear
    = False``, Issue #158); this pins what the integrated chart must become there.
    On the cell :math:`[r_c - h/2, r_c + h/2]`, :math:`r = r_c + \xi h/2`,
    :math:`\varepsilon = h/r_c`, the Gram of :math:`\{1, \xi\}` is FULL and
    position-dependent, so the data-side map is the cell mass matrix, of which
    the slab ``V·diag(1, θ)`` is the flat-measure value:

    * sphere (:math:`r^2\mathrm dr`): :math:`V = h(h^2 + 12 r_c^2)/12`,
      :math:`M_{01}/V = 4\varepsilon/(\varepsilon^2 + 12)`,
      :math:`M_{11}/V = (3\varepsilon^2 + 20)/(5(\varepsilon^2 + 12))`;
      pole cell :math:`r_c = h/2`: :math:`M_{01}/V = 1/2`, :math:`M_{11}/V = 2/5`;
    * cylinder (:math:`r\,\mathrm dr`, per radian): :math:`V = h r_c`,
      :math:`M_{01}/V = \varepsilon/6`, :math:`M_{11}/V = 1/3`;
      pole cell: :math:`M_{01}/V = 1/3`, :math:`M_{11}/V = 1/3`;
    * the flat limit :math:`\varepsilon \to 0` recovers :math:`M/V = \mathrm{diag}(1, \theta)`
      in both geometries (:func:`mass_1d` ``/h``).
    """
    xi = sp.Symbol("xi", real=True)
    h, rc, eps = sp.symbols("h r_c epsilon", positive=True)
    theta = sp.Rational(1, 3)
    flat = mass_1d(sp.Integer(1), theta)
    expected = {
        "sphere": (2, h * (h**2 + 12 * rc**2) / 12, 4 * eps / (eps**2 + 12),
                   (3 * eps**2 + 20) / (5 * (eps**2 + 12)), sp.Rational(1, 2), sp.Rational(2, 5)),
        "cylinder": (1, h * rc, eps / 6, theta, theta, theta),
    }
    out: dict = {}
    passed = True
    for name, (power, V_exp, m01_exp, m11_exp, m01_pole, m11_pole) in expected.items():
        r = rc + xi * h / 2
        measure = r**power * h / 2
        M = sp.Matrix(2, 2, lambda i, j: sp.integrate(xi ** (i + j) * measure, (xi, -1, 1)))
        V = sp.simplify(M[0, 0])
        m01 = sp.simplify((M[0, 1] / V).subs(h, eps * rc))
        m11 = sp.simplify((M[1, 1] / V).subs(h, eps * rc))
        pole = (sp.simplify((M[0, 1] / V).subs(rc, h / 2)), sp.simplify((M[1, 1] / V).subs(rc, h / 2)))
        flat_limit = sp.Matrix([[1, sp.limit(m01, eps, 0)], [sp.limit(m01, eps, 0), sp.limit(m11, eps, 0)]])
        ok = (
            sp.simplify(V - V_exp) == 0
            and sp.simplify(m01 - m01_exp) == 0
            and sp.simplify(m11 - m11_exp) == 0
            and pole == (m01_pole, m11_pole)
            and (flat_limit - flat).is_zero_matrix
            and m01 != 0  # the Gram is NOT diagonal: a scalar θ cannot carry it
        )
        out[name] = {"V": V, "M01_over_V": m01, "M11_over_V": m11, "pole": pole, "ok": ok}
        passed = passed and ok
    return {
        "name": "V_gram_curv: LD cell Gram under r²dr / r dr (full, position-dependent; forward-looking)",
        "sphere": out["sphere"],
        "cylinder": out["cylinder"],
        "pass": bool(passed),
    }


# ═══════════════════════════════════════════════════════════════════════
# CLI — print all derivations
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    for fn in (
        derive_d1_reduction_to_production,
        derive_d1_kernel_view_equals,
        derive_d1_scan_view_equals,
        derive_d2_exact_on_bilinear,
        derive_d3_assembles,
        derive_d1_transpose_equals_At_Minv,
        derive_octant_frame_sign_is_involution,
        derive_d1_integrated_slope_source_is_the_galerkin_moment,
        derive_d1_coordinate_pairing_is_the_integral,
        derive_curvilinear_ld_cell_gram,
    ):
        result = fn()
        status = "PASS" if result["pass"] else "FAIL"
        print(f"[{status}] {result['name']}")
        for key, val in result.items():
            if key in ("name", "pass"):
                continue
            print(f"    {key} = {val}")
        print()
