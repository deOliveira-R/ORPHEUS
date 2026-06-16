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

   A\,\vec\psi = \vec R, \qquad
   A = G + F_{\rm out} + \Sigma_t M, \qquad
   \vec R = M\,\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces},

a ``2^d × 2^d`` dense non-symmetric solve.  ``d=1`` is the
"Kronecker-with-one-factor" identity — it reduces EXACTLY to the
production slab 2×2 (``orpheus.sn.spatial.linear_discontinuous``), the
``xy`` coupling falls out of the algebra for ``d≥2``, and NO 4×4/8×8
entry is hand-transcribed.

Verification functions
======================

* :func:`derive_d1_reduction_to_production` — Oracle (i): the assembled
  ``d=1`` ``A`` / ``R`` / Schur ``S`` / slope ``ψ̂`` / ``D₂'`` reduce
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

* :mod:`orpheus.sn.spatial.linear_discontinuous` — the production 1-D LD
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
    M = masses[0]
    for k in range(1, d):
        M = sp.Matrix(sp.kronecker_product(M, masses[k]))

    size = 2**d
    G = sp.zeros(size, size)
    F_out = sp.zeros(size, size)
    for a in range(d):
        # Active axis a carries the gradient / outflow; transverse axes the mass.
        grad_factors = [grad_1d() if k == a else masses[k] for k in range(d)]
        face_factors = [fout_1d() if k == a else masses[k] for k in range(d)]
        Ga = grad_factors[0]
        Foa = face_factors[0]
        for k in range(1, d):
            Ga = sp.Matrix(sp.kronecker_product(Ga, grad_factors[k]))
            Foa = sp.Matrix(sp.kronecker_product(Foa, face_factors[k]))
        G += mus[a] * Ga
        F_out += mus[a] * Foa

    A = G + F_out + sig_t * M
    return {"A": A, "M": M, "G": G, "F_out": F_out, "size": size, "d": d}


def assemble_inflow_axis(
    hs: list[sp.Expr],
    mus: list[sp.Expr],
    axis: int,
    upstream_face_moments: sp.Matrix,
    theta: sp.Expr = THETA,
) -> sp.Matrix:
    r"""Inflow (upwind) RHS contribution from one upstream face.

    The incoming face normal to ``axis`` is the upstream neighbour's
    outflow trace — a ``2^{d-1}``-moment transverse object
    (``upstream_face_moments``, ordered as the transverse-axes Kronecker
    moment vector).  It is weighted into the test functions ``B(-1) =
    [1,-1]`` on ``axis`` and the **mass** on the transverse axes, times
    ``|μ_axis|``.  Returns a ``2^d``-vector RHS contribution.
    """
    d = len(hs)
    masses = [mass_1d(h, theta) for h in hs]
    # Transverse mass applied to the (2^{d-1}) upstream face moments.
    transverse_mass = None
    for k in range(d):
        if k == axis:
            continue
        transverse_mass = (
            masses[k]
            if transverse_mass is None
            else sp.Matrix(sp.kronecker_product(transverse_mass, masses[k]))
        )
    weighted_face = (
        upstream_face_moments
        if transverse_mass is None
        else transverse_mass * upstream_face_moments
    )
    # Re-insert the active-axis test weighting B(-1)=[1,-1] by Kronecker'ing
    # it into the active slot, scaled by |μ_axis| (the streaming face flux is
    # |μ_axis| · ∮ B_i ψ_upwind — the same |μ_a| factor F_out carries).
    # ``weighted_face`` already lives in the transverse-axes Kronecker order,
    # so the active axis sits at the boundary of the product (axis 0 -> outer,
    # axis d-1 -> inner).  The oracles use d=2 (axis ∈ {0, 1}); the general
    # interior-axis interleave (d≥3, 0 < axis < d-1) is Branch-2 territory.
    mu_axis = mus[axis]
    # ``Matrix * scalar`` (scalar on the right) stays Matrix-typed; the
    # ``scalar * Matrix`` order returns an Expr in the SymPy stubs.
    if axis == 0:
        return sp.Matrix(sp.kronecker_product(fin_trace_weight(), weighted_face)) * mu_axis
    if axis == d - 1:
        return sp.Matrix(sp.kronecker_product(weighted_face, fin_trace_weight())) * mu_axis
    raise NotImplementedError(
        "assemble_inflow_axis supports axis in {0, d-1} (all the oracles need); "
        f"the general interior-axis interleave (d={d}, axis={axis}) is deferred "
        "to Branch 2."
    )


def per_cell_solve(assembled: dict, rhs: sp.Matrix) -> sp.Matrix:
    r"""Solve ``A ψ⃗ = R`` for the ``2^d``-moment cell unknown vector."""
    return assembled["A"].LUsolve(rhs)


def downstream_face_trace(psi_moments: sp.Matrix, d: int) -> sp.Expr:
    r"""Outgoing face flux = trace of the cell function at the downstream node.

    For d=1 this is ``B(+1)·ψ⃗ = ψ̄ + ψ̂`` (the LM-1989 Eq. 4.3c closure).
    For ``d≥2`` the downstream face (normal to the sweep axis) is a
    ``2^{d-1}``-moment transverse object — the face-cochain widening that
    Branch 2 owns; here we return the d=1 scalar trace (the oracle uses d=1).
    """
    if d == 1:
        return (sp.Matrix([1, 1]).T * psi_moments)[0]
    raise NotImplementedError(
        "downstream_face_trace returns the d=1 scalar trace; the d>=2 face is "
        "a 2^{d-1}-moment transverse object (Branch 2 / face-cochain widening)."
    )


# ═══════════════════════════════════════════════════════════════════════
# Oracle (i) — the d=1 reduction to the production slab LD
# ═══════════════════════════════════════════════════════════════════════


def derive_d1_reduction_to_production() -> dict:
    r"""V_d1 — the assembled d=1 UBLD reduces to the production slab 2×2.

    Proves: the Kronecker-with-one-factor assembly reproduces the
    production ``A`` / ``R`` / Schur ``S`` / slope ``ψ̂`` / ``D₂'`` closed
    forms (``orpheus.sn.spatial.linear_discontinuous._LDCellTerms`` /
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

    Proves: the ×V scan coefficients (``m = |μ|A_down``, ``t = Σ_t V``,
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

    # production ×V scan (affine_scan_coefficients, lines 564-571); slab A_down=1, V=h.
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
    ):
        result = fn()
        status = "PASS" if result["pass"] else "FAIL"
        print(f"[{status}] {result['name']}")
        for key, val in result.items():
            if key in ("name", "pass"):
                continue
            print(f"    {key} = {val}")
        print()
