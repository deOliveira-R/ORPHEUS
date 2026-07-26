r"""Adjoint-weighted (eigenvalue-consistent) collapse rules — the algebra of record.

The P6 (#281) taxonomy derivation: which per-channel collapse rule makes the
coarse multiplication factor **first-order stationary** in the flux shapes —
the eigenvalue-consistent generalization of the forward (flux-weighted)
homogenization/condensation that already ships
(:meth:`orpheus.sn.solution.Solution.homogenize` /
:meth:`~orpheus.sn.solution.Solution.condense`).

Every rule implemented by the production collapse MUST be derivable here
first (the algebra-of-record discipline): the SymPy theorems below are the
source of truth for the theory section in
:doc:`/theory/foundations/frame` and for the ``adjoint=`` arm of the
collapse machinery.

The worth functional and the two error sources
==============================================

For the generalized eigenproblem :math:`A\varphi = \tfrac{1}{k}F\varphi`
with left (adjoint) eigenvector :math:`\varphi^*`, a perturbation
:math:`(\delta A, \delta F)` shifts the eigenvalue by (T0, the keystone)

.. math::

    \delta\mu \;=\;
    \frac{\langle \varphi^*, (\delta A - \mu\,\delta F)\,\varphi\rangle}
         {\langle \varphi^*, F \varphi\rangle}
    \;+\; \mathcal O(\delta^2),
    \qquad \mu = \tfrac{1}{k} .

Replacing the fine per-cell cross sections by region-collapsed constants
**on the same fine mesh** is such a perturbation — the *XS-collapse worth*.
A collapse rule that zeroes each region's worth term therefore kills the
first-order eigenvalue error of the collapse itself.  The remaining coarse
re-solve error (the coarse *discretization* of streaming) is shared by every
weighting and is NOT touched by these rules — the honest scope of the C2
gate (`p6_adjoint_verification_spec.md` §4).

The taxonomy (spatial axis; theorems T1–T5)
===========================================

With :math:`w^{\pm}_{i,g} \equiv \varphi^*_{i,g}\varphi_{i,g}` the *bilinear
pair weight* on fine cell :math:`i`, group :math:`g` (volume :math:`V_i`):

* **T1 — vector channels** (:math:`\Sigma_t,\Sigma_c,\Sigma_L,\Sigma_f`):
  :math:`\Sigma_{R,g} = \sum_i V_i w^{\pm}_{i,g}\Sigma_{i,g} \big/
  \sum_i V_i w^{\pm}_{i,g}` — the product weight :math:`\varphi^*\odot\varphi`,
  NOT a bare :math:`\varphi\to\varphi^*` swap.
* **T2 — matrix channels** (:math:`\Sigma_{s,\ell}[g'\!\to\!g]`,
  :math:`\Sigma_{2n}`): the **per-pair** rule — weight
  :math:`\varphi^*_{i,g}\varphi_{i,g'}` (sink adjoint × source flux) per
  :math:`(g',g)` entry.  The source-product weight
  :math:`(\varphi^*\varphi)_{g'}` (what the existing per-(cell, group)
  plumbing would broadcast) does NOT zero the worth.
* **T3 — the fission dyad** (:math:`\chi\otimes\nu\Sigma_f`): a per-pair
  collapse of the rank-1 dyad is not rank-1, but a ``Mixture`` stores the
  factors.  The **mixed-fold factored rule** is exact for the *total*
  region worth (all the eigenvalue needs): with the cell emission
  importances :math:`\iota_i = \sum_g \varphi^*_{i,g}\chi_{i,g}` (fine) and
  :math:`\tilde\iota_i = \sum_g \varphi^*_{i,g}\chi_{R,g}` (collapsed),

  .. math::

      \nu\Sigma_{f,R,g'} \;=\;
      \frac{\sum_i V_i\,\iota_i\,\nu\Sigma_{f,i,g'}\,\varphi_{i,g'}}
           {\sum_i V_i\,\tilde\iota_i\,\varphi_{i,g'}}

  zeroes the total fission worth for **any** simplex :math:`\chi_R`
  (theorem); the canonical :math:`\chi_R` is the adjoint-weighted-emission
  convex average (weights :math:`q_i = V_i\,\iota_i\,p_i`,
  :math:`p_i=\sum_{g'}\nu\Sigma_{f,i,g'}\varphi_{i,g'}`), which stays a
  simplex and degenerates to the production-weighted forward rule.
* **T4 — the balance trade-off**: worth-exact rules BREAK the definitional
  total-XS balance :math:`\Sigma_t = \Sigma_c+\Sigma_L+\Sigma_f+
  \mathrm{rowsum}(\Sigma_{s0})+\mathrm{rowsum}(\Sigma_{2n})` (the per-pair
  matrix rowsums do not reproduce the T1-collapsed total).  Conversely,
  restoring balance by defining :math:`\Sigma_t := \text{sum of parts}`
  re-introduces a first-order worth residual.  Worth-exactness and balance
  are mutually exclusive under :math:`\varphi^*\neq` const — the classical
  property of bilinear-weighted constants (reactivity-preserving, not
  rate-preserving).  At :math:`\varphi^*=` const both hold (the forward
  degenerate).
* **T5 — the forward rule's first-order error**: the flux-weighted
  (forward) collapse leaves a generically nonzero first-order worth — the
  theoretical basis of the C2 comparative gate (forward gap first-order,
  adjoint-weighted gap second-order).

Angular honest-scope note (P0 pairing).  The collision worth pairs the
ANGULAR fluxes, :math:`\sum_n w_n \psi^*_n\psi_n`; the scalar product
:math:`\varphi^*\varphi` is its P0 truncation.  It is EXACT for the
isotropic operators (scattering source, fission) and truncated for the
:math:`\Sigma_t` collision term — the classical practice for
bilinear-weighted constants.  The angular-product weight is a documented
seam, not implemented in P6.

The energy axis (T6): condensation is pure projection
=====================================================

Energy condensation has NO streaming carve (streaming does not couple
groups): the coarse pencil can be an EXACT left-diagonal rescaling of the
Petrov-Galerkin projection of the fine pencil (test
:math:`\varphi^*\cdot\mathbf 1_G`, trial :math:`\varphi\cdot\mathbf 1_{G'}`)
— so with the *consistent-denominator* convention the condensed pencil
reproduces the fine :math:`k` EXACTLY when condensed with the true spectra
(T6a), and the error under a perturbed spectrum is second-order for the
bilinear rules vs first-order for the forward rules (T6b — the C2 basis on
the energy axis).  The exactness selects the denominators: every channel of
row :math:`G` divides by the same pairing so that ONE coarse flux vector
satisfies all rows.  On the energy axis the block pair weight is separable
(:math:`\sum_{g\in G,g'\in G'}\varphi^*_g\varphi_{g'} =
\Phi^*_G\,\Phi_{G'}`), which is what makes the exact convention available —
and it makes the condensed fission dyad exactly rank-1 (the factored form
survives condensation, unlike the spatial axis).  The two vector-rule
conventions (diagonal-pair bilinear vs separable-denominator) and the
coarse-flux normalization they imply are adjudicated by T6 and reconciled
against the classical literature (`.claude/plans/p6_literature_memo.md`).

Verification discipline
=======================

Every theorem below verifies its identity with exact (rational/symbolic)
arithmetic and raises :class:`DerivationError` on failure — no bare
``assert`` (this module is not pytest-collected, so bare asserts would be
stripped under the canonical ``python -O``; vv-principles Mode 8).  The
pinning suite is ``tests/derivations/test_homogenization_rules.py``.

Usage::

    .venv/bin/python -m orpheus.derivations.common.homogenization  # narrate all
"""

from __future__ import annotations

from typing import Callable

import sympy as sp

__all__ = [
    "DerivationError",
    "derive_first_order_k_shift",
    "derive_vector_channel_rule",
    "derive_matrix_channel_rule",
    "derive_fission_factored_rule",
    "derive_balance_tradeoff",
    "derive_forward_rule_first_order_error",
    "derive_energy_condensation_exactness",
    "run_all",
]


class DerivationError(AssertionError):
    """A symbolic identity this module claims failed to verify."""


def _require(condition: bool, message: str) -> None:
    """Raise :class:`DerivationError` unless ``condition`` (Mode-8-proof assert)."""
    if not condition:
        raise DerivationError(message)


def _is_zero(expr: sp.Expr, message: str) -> None:
    """Verify ``expr`` simplifies to exactly zero (raises otherwise)."""
    simplified = sp.simplify(sp.together(expr))
    _require(simplified == 0, f"{message}: residual = {simplified}")


def _is_nonzero(expr: sp.Expr, message: str) -> None:
    """Verify ``expr`` is NOT identically zero (raises if it vanishes)."""
    simplified = sp.simplify(sp.together(expr))
    _require(simplified != 0, f"{message}: unexpectedly vanished")


# ═══════════════════════════════════════════════════════════════════════
# Common symbol grids — 2 fine cells in one region, 2 groups.
#
# The worth identities are per-region sums of per-cell terms, so 2 cells ×
# 2 groups with fully generic symbols is the general case for every
# identity proved below (each theorem is linear in the per-cell channel
# symbols; no cancellation depends on the counts).
# ═══════════════════════════════════════════════════════════════════════

_N, _NG = 2, 2
_CELLS = range(_N)
_GROUPS = range(_NG)

V = [sp.Symbol(f"V_{i}", positive=True) for i in _CELLS]
phi = [[sp.Symbol(f"phi_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
phis = [[sp.Symbol(f"phis_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]


def _pair_weight(i: int, g: int) -> sp.Expr:
    """The bilinear pair weight w± = φ*·φ at (cell i, group g)."""
    return phis[i][g] * phi[i][g]


def _generic_rationals() -> dict[sp.Symbol, sp.Rational]:
    """A generic exact-rational instantiation for counterexample evaluation.

    Values are deliberately asymmetric (no two equal, no proportionality
    between the φ and φ* rows) so an expression that vanishes on the
    instance vanishes for structural reasons, and an expression claimed
    nonzero cannot vanish by accidental symmetry.
    """
    subs: dict[sp.Symbol, sp.Rational] = {}
    primes = iter([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43])
    for i in _CELLS:
        subs[V[i]] = sp.Rational(next(primes), 4)
        for g in _GROUPS:
            subs[phi[i][g]] = sp.Rational(next(primes), 6)
            subs[phis[i][g]] = sp.Rational(next(primes), 5)
    return subs


# ═══════════════════════════════════════════════════════════════════════
# T0 — the keystone: first-order eigenvalue shift of a generalized pencil
# ═══════════════════════════════════════════════════════════════════════

def derive_first_order_k_shift() -> None:
    r"""T0: verify :math:`\delta\mu = \langle x^*,(\delta A-\mu\,\delta F)x\rangle/\langle x^*,Fx\rangle`.

    The operator proof is three lines — differentiate
    :math:`(A(\varepsilon)-\mu(\varepsilon)F(\varepsilon))x(\varepsilon)=0`,
    left-multiply by the unperturbed left eigenvector :math:`x^*`, and the
    :math:`\delta x` term dies because :math:`x^{*\mathsf T}(A-\mu F)=0`.
    Here the identity is verified end-to-end with exact arithmetic: an
    asymmetric rational pencil, generic rational perturbations, the exact
    perturbed eigenvalue's :math:`\varepsilon`-series, and the formula —
    coefficient by coefficient.

    This is what connects the T1–T3 worth-zeroing rules to
    :math:`\delta k = \mathcal O(\delta^2)`: zero the region worths and the
    first-order term of the collapse error is identically zero.
    """
    print("=" * 68)
    print("T0. First-order eigenvalue shift of a generalized pencil")
    print("=" * 68)

    eps = sp.Symbol("varepsilon")

    # Asymmetric pencil with rational spectrum: A0 = F0 @ M, eig(pencil) = eig(M).
    M = sp.Matrix([[2, 1], [0, 3]])
    F0 = sp.Matrix([[3, 1], [1, 2]])       # symmetric positive definite, det=5
    A0 = F0 * M
    mu0 = sp.Integer(3)                     # dominant eigenvalue of M

    x = sp.Matrix([1, 1])                   # right: (M − 3I)x = 0
    y = sp.Matrix([0, 1])                   # left of M: yᵀ(M − 3I) = 0
    xs = (F0.T).solve(y)                    # pencil left eigenvector: x*ᵀF0 = yᵀ
    _is_zero((xs.T * (A0 - mu0 * F0)).norm() ** 2, "x* is not a left null vector")
    _is_zero(((A0 - mu0 * F0) * x).norm() ** 2, "x is not a right null vector")

    dA = sp.Matrix([[sp.Rational(1, 3), sp.Rational(2, 7)],
                    [sp.Rational(5, 11), sp.Rational(3, 13)]])
    dF = sp.Matrix([[sp.Rational(2, 5), sp.Rational(1, 9)],
                    [sp.Rational(4, 7), sp.Rational(1, 6)]])

    # Exact perturbed eigenvalue: the char-poly root continuous with μ0.
    mu = sp.Symbol("mu")
    charpoly = ((A0 + eps * dA) - mu * (F0 + eps * dF)).det()
    roots = sp.solve(sp.Eq(charpoly, 0), mu)
    branch = [r for r in roots if sp.simplify(r.subs(eps, 0) - mu0) == 0]
    _require(len(branch) == 1, "could not isolate the μ0-continuous branch")
    series = sp.series(branch[0], eps, 0, 2).removeO()

    predicted = (xs.T * (dA - mu0 * dF) * x)[0, 0] / (xs.T * F0 * x)[0, 0]
    measured = sp.simplify(series.coeff(eps, 1))
    _is_zero(measured - predicted, "first-order shift formula mismatch")
    print(f"  δμ (series)  = {measured}")
    print(f"  δμ (formula) = {sp.simplify(predicted)}")
    print("  ✓ δμ = ⟨x*,(δA − μδF)x⟩ / ⟨x*,Fx⟩ — exact to O(ε²)\n")


# ═══════════════════════════════════════════════════════════════════════
# T1 — vector channels: the bilinear pair-weight rule
# ═══════════════════════════════════════════════════════════════════════

def _vector_worth(sigma: list[list[sp.Symbol]], rule: Callable[[int], sp.Expr]) -> list[sp.Expr]:
    """Per-group P0 reaction worth W_g = Σ_i V_i φ*_{i,g}(Σ_{R,g} − Σ_{i,g})φ_{i,g}."""
    return [
        sum(V[i] * phis[i][g] * (rule(g) - sigma[i][g]) * phi[i][g] for i in _CELLS)
        for g in _GROUPS
    ]


def derive_vector_channel_rule() -> None:
    r"""T1: the pair weight :math:`\varphi^*\!\odot\varphi` zeroes the vector-channel worth.

    Solving :math:`W_g = 0` for :math:`\Sigma_{R,g}` yields exactly the
    bilinear :math:`\langle\varphi^*,\Sigma\varphi\rangle_R /
    \langle\varphi^*,\varphi\rangle_R` — and the forward (flux-weighted)
    rule leaves :math:`W_g \neq 0` on a generic instance (its residual is
    first-order in the :math:`\varphi^*` tilt; T5 consolidates).

    Production consequence: the sigma-frame test weight under ``adjoint=``
    is the PRODUCT array ``phi_star * phi`` — the drafted "(φ → φ*)"
    phrasings are corrected by B4 to this bilinear form.
    """
    print("=" * 68)
    print("T1. Vector channels — the bilinear pair weight φ*⊙φ")
    print("=" * 68)

    sigma = [[sp.Symbol(f"sigma_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]

    def bilinear(g: int) -> sp.Expr:
        num = sum(V[i] * _pair_weight(i, g) * sigma[i][g] for i in _CELLS)
        den = sum(V[i] * _pair_weight(i, g) for i in _CELLS)
        return num / den

    for g, worth in enumerate(_vector_worth(sigma, bilinear)):
        _is_zero(worth, f"T1: bilinear rule leaves worth in group {g}")
    print("  ✓ Σ_{R,g} = ⟨φ*,Σφ⟩_R / ⟨φ*,φ⟩_R zeroes W_g identically (all g)")

    # Uniqueness: solving W_g = 0 returns the same rule.
    sig_R = sp.Symbol("sigma_R", positive=True)
    worth_g0 = sum(
        V[i] * phis[i][0] * (sig_R - sigma[i][0]) * phi[i][0] for i in _CELLS
    )
    solved = sp.solve(sp.Eq(worth_g0, 0), sig_R)
    _require(len(solved) == 1, "T1: worth-zeroing rule is not unique")
    _is_zero(solved[0] - bilinear(0), "T1: solved rule differs from the bilinear")
    print("  ✓ W_g = 0 has the bilinear as its UNIQUE solution")

    def forward(g: int) -> sp.Expr:
        num = sum(V[i] * phi[i][g] * sigma[i][g] for i in _CELLS)
        den = sum(V[i] * phi[i][g] for i in _CELLS)
        return num / den

    subs = _generic_rationals()
    for i in _CELLS:
        for g in _GROUPS:
            subs[sigma[i][g]] = sp.Rational(7 * i + 3 * g + 2, 9)
    residual = _vector_worth(sigma, forward)[0].subs(subs)
    _is_nonzero(residual, "T1: forward rule accidentally zeroed the worth")
    print(f"  ✓ forward (flux-weighted) rule leaves W_0 = {sp.nsimplify(residual)} ≠ 0\n")


# ═══════════════════════════════════════════════════════════════════════
# T2 — matrix channels: the per-pair sink×source rule
# ═══════════════════════════════════════════════════════════════════════

def derive_matrix_channel_rule() -> None:
    r"""T2: matrix channels need the per-pair weight :math:`\varphi^*_g\,\varphi_{g'}`.

    The scattering worth decomposes per :math:`(g',g)` entry;
    each term zeroes iff

    .. math::

        \Sigma_{s,R}[g'\!\to\!g] \;=\;
        \frac{\sum_i V_i\,\varphi^*_{i,g}\,\Sigma_{s,i}[g'\!\to\!g]\,\varphi_{i,g'}}
             {\sum_i V_i\,\varphi^*_{i,g}\,\varphi_{i,g'}} .

    The source-product rule — weight :math:`(\varphi^*\varphi)_{g'}`, which
    is what the existing per-(cell, group) leading-aligned broadcast would
    produce — does NOT zero the off-diagonal worth (generic counterexample).
    Production consequence: the ``adjoint=`` arm needs per-pair plumbing for
    matrix channels (B2's structural work), not a weight-array swap.
    """
    print("=" * 68)
    print("T2. Matrix channels — per-pair sink×source weight φ*_g·φ_g'")
    print("=" * 68)

    sig_s = [
        [[sp.Symbol(f"sigs_{i}{gf}{gt}", positive=True) for gt in _GROUPS]
         for gf in _GROUPS]
        for i in _CELLS
    ]

    def per_pair(gf: int, gt: int) -> sp.Expr:
        num = sum(V[i] * phis[i][gt] * sig_s[i][gf][gt] * phi[i][gf] for i in _CELLS)
        den = sum(V[i] * phis[i][gt] * phi[i][gf] for i in _CELLS)
        return num / den

    def source_product(gf: int, gt: int) -> sp.Expr:
        num = sum(V[i] * _pair_weight(i, gf) * sig_s[i][gf][gt] for i in _CELLS)
        den = sum(V[i] * _pair_weight(i, gf) for i in _CELLS)
        return num / den

    def worth(rule: Callable[[int, int], sp.Expr], gf: int, gt: int) -> sp.Expr:
        return sum(
            V[i] * phis[i][gt] * (rule(gf, gt) - sig_s[i][gf][gt]) * phi[i][gf]
            for i in _CELLS
        )

    for gf in _GROUPS:
        for gt in _GROUPS:
            _is_zero(worth(per_pair, gf, gt), f"T2: per-pair leaves worth at ({gf},{gt})")
    print("  ✓ per-pair rule zeroes the scattering worth for every (g',g)")

    subs = _generic_rationals()
    for i in _CELLS:
        for gf in _GROUPS:
            for gt in _GROUPS:
                subs[sig_s[i][gf][gt]] = sp.Rational(11 * i + 5 * gf + 2 * gt + 1, 8)
    off_diag = worth(source_product, 0, 1).subs(subs)
    _is_nonzero(off_diag, "T2: source-product rule accidentally zeroed (0→1)")
    print(f"  ✓ source-product rule leaves W[0→1] = {sp.nsimplify(off_diag)} ≠ 0")

    # Degeneracy: on the diagonal (g'=g) — and for flat φ* — both coincide.
    for g in _GROUPS:
        _is_zero(
            per_pair(g, g) - source_product(g, g),
            "T2: diagonal per-pair ≠ source-product",
        )
    flat = {phis[i][g]: sp.Symbol("c", positive=True) for i in _CELLS for g in _GROUPS}

    def fwd_source(gf: int, gt: int) -> sp.Expr:
        num = sum(V[i] * phi[i][gf] * sig_s[i][gf][gt] for i in _CELLS)
        den = sum(V[i] * phi[i][gf] for i in _CELLS)
        return num / den

    _is_zero(
        (per_pair(0, 1) - fwd_source(0, 1)).subs(flat),
        "T2: flat-φ* degenerate is not the forward source-weighted rule",
    )
    print("  ✓ degenerates: diagonal ≡ source-product; flat φ* ≡ forward rule\n")


# ═══════════════════════════════════════════════════════════════════════
# T3 — the fission dyad: mixed-fold factored rule, χ_R free (then canonical)
# ═══════════════════════════════════════════════════════════════════════

def derive_fission_factored_rule() -> None:
    r"""T3: the factored fission rule that zeroes the TOTAL region worth exactly.

    The per-pair collapse of :math:`\chi\otimes\nu\Sigma_f` is not rank-1,
    but the eigenvalue only needs the TOTAL fission worth to vanish (T0
    contracts everything against one scalar).  Theorem: with

    .. math::

        \iota_i = \sum_g \varphi^*_{i,g}\chi_{i,g}, \qquad
        \tilde\iota_i = \sum_g \varphi^*_{i,g}\chi_{R,g}, \qquad
        \nu\Sigma_{f,R,g'} =
        \frac{\sum_i V_i\,\iota_i\,\nu\Sigma_{f,i,g'}\,\varphi_{i,g'}}
             {\sum_i V_i\,\tilde\iota_i\,\varphi_{i,g'}},

    the total worth :math:`W_F = \sum_{i,g,g'} V_i \varphi^*_{i,g}
    (\chi_{R,g}\nu\Sigma_{f,R,g'} - \chi_{i,g}\nu\Sigma_{f,i,g'})
    \varphi_{i,g'}` vanishes **identically for any** :math:`\chi_R`
    (verified with free simplex-constrained symbols).  The canonical
    :math:`\chi_R` — the adjoint-weighted-emission convex average with
    weights :math:`q_i = V_i\,\iota_i\,p_i` — is a simplex by construction,
    and the whole rule degenerates at flat :math:`\varphi^*` to today's
    forward pair (flux-weighted :math:`\nu\Sigma_f`, production-weighted
    :math:`\chi`).
    """
    print("=" * 68)
    print("T3. Fission dyad — mixed-fold νΣf rule; χ_R free, canonical convex")
    print("=" * 68)

    nsf = [[sp.Symbol(f"nsf_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
    chi = [[sp.Symbol(f"chi_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
    # Free collapsed spectrum, simplex-constrained: χ_R = (t, 1−t).
    t = sp.Symbol("t", positive=True)
    chi_R = [t, 1 - t]

    iota = [sum(phis[i][g] * chi[i][g] for g in _GROUPS) for i in _CELLS]
    iota_R = [sum(phis[i][g] * chi_R[g] for g in _GROUPS) for i in _CELLS]
    p = [sum(nsf[i][gp] * phi[i][gp] for gp in _GROUPS) for i in _CELLS]

    def nsf_R(gp: int) -> sp.Expr:
        num = sum(V[i] * iota[i] * nsf[i][gp] * phi[i][gp] for i in _CELLS)
        den = sum(V[i] * iota_R[i] * phi[i][gp] for i in _CELLS)
        return num / den

    total_worth = sum(
        V[i] * phis[i][g] * (chi_R[g] * nsf_R(gp) - chi[i][g] * nsf[i][gp]) * phi[i][gp]
        for i in _CELLS for g in _GROUPS for gp in _GROUPS
    )
    _is_zero(total_worth, "T3: mixed-fold rule leaves total fission worth")
    print("  ✓ W_F ≡ 0 for ANY simplex χ_R (t left free) — the χ slot is exactness-free")

    # Canonical χ_R: adjoint-weighted-emission convex average (q_i = V_i ι_i p_i).
    q = [V[i] * iota[i] * p[i] for i in _CELLS]
    chi_canonical = [
        sum(q[i] * chi[i][g] for i in _CELLS) / sum(q[i] for i in _CELLS)
        for g in _GROUPS
    ]
    # Simplex: needs Σ_g χ_{i,g} = 1 per cell — impose and check.
    simplex = {chi[i][1]: 1 - chi[i][0] for i in _CELLS}
    _is_zero(
        (chi_canonical[0] + chi_canonical[1]).subs(simplex) - 1,
        "T3: canonical χ_R does not sum to 1",
    )
    print("  ✓ canonical χ_R (q_i = V_i·ι_i·p_i convex average) is a simplex")

    # Degeneracy at flat φ* = c: νΣf_R → flux-weighted; χ_R → production-weighted.
    c = sp.Symbol("c", positive=True)
    flat = {phis[i][g]: c for i in _CELLS for g in _GROUPS}
    fwd_nsf = (
        sum(V[i] * nsf[i][0] * phi[i][0] for i in _CELLS)
        / sum(V[i] * phi[i][0] for i in _CELLS)
    )
    _is_zero(
        (nsf_R(0) - fwd_nsf).subs(simplex).subs({t: sp.Rational(1, 3)}).subs(flat),
        "T3: flat-φ* νΣf_R is not the forward flux-weighted rule",
    )
    fwd_chi0 = (
        sum(V[i] * p[i] * chi[i][0] for i in _CELLS)
        / sum(V[i] * p[i] for i in _CELLS)
    )
    _is_zero(
        (chi_canonical[0] - fwd_chi0).subs(simplex).subs(flat),
        "T3: flat-φ* χ_R is not the forward production-weighted rule",
    )
    print("  ✓ flat-φ* degenerates: νΣf_R → flux-weighted, χ_R → production-weighted\n")


# ═══════════════════════════════════════════════════════════════════════
# T4 — the balance trade-off
# ═══════════════════════════════════════════════════════════════════════

def derive_balance_tradeoff() -> None:
    r"""T4: worth-exactness and the total-XS balance are mutually exclusive.

    Fine materials satisfy :math:`\Sigma_{t} = \Sigma_a +
    \mathrm{rowsum}(\Sigma_{s})` per source group (collapse the removal
    channels into one :math:`\Sigma_a` symbol WLOG).  Collapsing
    :math:`\Sigma_t` and :math:`\Sigma_a` with T1 and :math:`\Sigma_s` with
    T2 leaves the balance residual

    .. math::

        B_{g'} = \Sigma_{t,R,g'} - \Sigma_{a,R,g'} -
                 \sum_g \Sigma_{s,R}[g'\!\to\!g]

    generically NONZERO (the per-pair rowsum re-weights each sink term by a
    different denominator).  Conversely defining :math:`\Sigma_{t,R} :=`
    sum-of-parts restores balance but re-introduces exactly that mismatch
    as a first-order collision worth.  At flat :math:`\varphi^*` the
    residual vanishes identically — the forward collapse enjoys both
    properties, the adjoint-weighted one must choose (and P6 chooses
    worth-exactness: the entire point of eigenvalue-consistent constants;
    the imbalance is documented and gate-pinned, never asserted away).
    """
    print("=" * 68)
    print("T4. Balance trade-off — worth-exact ⟺ balance broken (φ* non-flat)")
    print("=" * 68)

    sig_a = [[sp.Symbol(f"siga_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
    sig_s = [
        [[sp.Symbol(f"sigs_{i}{gf}{gt}", positive=True) for gt in _GROUPS]
         for gf in _GROUPS]
        for i in _CELLS
    ]
    # Fine balance: Σt is DERIVED (the identity holds cell-by-cell).
    sig_t = [
        [sig_a[i][g] + sum(sig_s[i][g][gt] for gt in _GROUPS) for g in _GROUPS]
        for i in _CELLS
    ]

    def t1(channel: list[list[sp.Expr]], g: int) -> sp.Expr:
        num = sum(V[i] * _pair_weight(i, g) * channel[i][g] for i in _CELLS)
        den = sum(V[i] * _pair_weight(i, g) for i in _CELLS)
        return num / den

    def t2(gf: int, gt: int) -> sp.Expr:
        num = sum(V[i] * phis[i][gt] * sig_s[i][gf][gt] * phi[i][gf] for i in _CELLS)
        den = sum(V[i] * phis[i][gt] * phi[i][gf] for i in _CELLS)
        return num / den

    balance_residual = [
        t1(sig_t, gf) - t1(sig_a, gf) - sum(t2(gf, gt) for gt in _GROUPS)
        for gf in _GROUPS
    ]

    subs = _generic_rationals()
    for i in _CELLS:
        for gf in _GROUPS:
            subs[sig_a[i][gf]] = sp.Rational(3 * i + 2 * gf + 1, 7)
            for gt in _GROUPS:
                subs[sig_s[i][gf][gt]] = sp.Rational(13 * i + 3 * gf + 5 * gt + 2, 12)
    residual_value = balance_residual[0].subs(subs)
    _is_nonzero(residual_value, "T4: balance residual vanished generically")
    print(f"  ✓ worth-exact collapse: balance residual B_0 = {sp.nsimplify(residual_value)} ≠ 0")

    flat = {phis[i][g]: sp.Symbol("c", positive=True) for i in _CELLS for g in _GROUPS}
    for gf in _GROUPS:
        _is_zero(
            balance_residual[gf].subs(flat),
            f"T4: flat-φ* balance residual nonzero in group {gf}",
        )
    print("  ✓ flat-φ* degenerate: balance restored identically (the forward property)")

    # The converse: Σt_R := sum-of-parts re-introduces the SAME mismatch as
    # a first-order collision worth (the worth of the difference).
    sum_of_parts = [
        t1(sig_a, gf) + sum(t2(gf, gt) for gt in _GROUPS) for gf in _GROUPS
    ]
    collision_worth = sum(
        V[i] * phis[i][gf] * (sum_of_parts[gf] - sig_t[i][gf]) * phi[i][gf]
        for i in _CELLS for gf in _GROUPS
    )
    mismatch_worth = sum(
        V[i] * phis[i][gf] * (sum_of_parts[gf] - t1(sig_t, gf)) * phi[i][gf]
        for i in _CELLS for gf in _GROUPS
    )
    _is_zero(
        collision_worth - mismatch_worth,
        "T4: converse decomposition failed",
    )
    _is_nonzero(collision_worth.subs(subs), "T4: sum-of-parts worth vanished")
    print("  ✓ converse: Σt_R := sum-of-parts ⇒ first-order collision worth = the mismatch\n")


# ═══════════════════════════════════════════════════════════════════════
# T5 — the forward rule's first-order error (the C2 basis)
# ═══════════════════════════════════════════════════════════════════════

def derive_forward_rule_first_order_error() -> None:
    r"""T5: the forward (flux-weighted) collapse leaves a first-order total worth.

    Consolidates the T1/T2 counterexamples into the single statement the C2
    gate rests on: with every channel collapsed by the forward rules, the
    TOTAL region worth :math:`\langle\varphi^*,(\delta A)\varphi\rangle_R`
    is generically nonzero and scales linearly with the adjoint tilt
    (verified by series in the tilt parameter) — while the T1+T2 rules
    leave exactly zero at every tilt.  Forward gap: first-order.
    Adjoint-weighted gap: the first-order term is identically zero.
    """
    print("=" * 68)
    print("T5. Forward rule — first-order in the adjoint tilt (the C2 basis)")
    print("=" * 68)

    sigma = [[sp.Symbol(f"sigma_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
    eps = sp.Symbol("varepsilon", positive=True)

    # φ* = 1 + ε·tilt (generic per-(i,g) tilts): ε parametrizes non-self-adjointness.
    tilt = [[sp.Rational(2 * i + 3 * g + 1, 5) for g in _GROUPS] for i in _CELLS]
    tilted = {phis[i][g]: 1 + eps * tilt[i][g] for i in _CELLS for g in _GROUPS}

    def forward(g: int) -> sp.Expr:
        num = sum(V[i] * phi[i][g] * sigma[i][g] for i in _CELLS)
        den = sum(V[i] * phi[i][g] for i in _CELLS)
        return num / den

    worth = sum(
        V[i] * phis[i][g] * (forward(g) - sigma[i][g]) * phi[i][g]
        for i in _CELLS for g in _GROUPS
    ).subs(tilted)

    subs = _generic_rationals()
    subs = {k: v for k, v in subs.items() if not str(k).startswith("phis")}
    for i in _CELLS:
        for g in _GROUPS:
            subs[sigma[i][g]] = sp.Rational(5 * i + 4 * g + 3, 11)
    series = sp.series(worth.subs(subs), eps, 0, 2).removeO()

    _is_zero(series.coeff(eps, 0), "T5: forward worth has an O(1) term (impossible)")
    _is_nonzero(series.coeff(eps, 1), "T5: forward worth first-order term vanished")
    print(f"  ✓ forward collapse worth = ({sp.nsimplify(series.coeff(eps, 1))})·ε + O(ε²)")
    print("  ✓ T1/T2 rules: worth ≡ 0 at every ε (proved in T1/T2)\n")


# ═══════════════════════════════════════════════════════════════════════
# T6 — the energy axis: condensation as exact projection
# ═══════════════════════════════════════════════════════════════════════

def derive_energy_condensation_exactness() -> None:
    r"""T6: consistent-denominator bilinear condensation reproduces k EXACTLY.

    0-D ∞-medium pencil, 4 fine groups → 2 coarse (blocks {0,1}, {2,3}),
    generic rational data (asymmetric :math:`\Sigma_s` with upscatter,
    :math:`\chi\not\parallel\nu\Sigma_f`).  With the fine solution pair
    :math:`\varphi = A^{-1}\chi` (up to scale), :math:`\varphi^* =
    A^{-\mathsf T}\nu\Sigma_f`, and :math:`k = \nu\Sigma_f^{\mathsf T}
    A^{-1}\chi` (rank-1 fission), condense with the CONSISTENT-DENOMINATOR
    convention — every channel of coarse row :math:`G` divides by the same
    :math:`\Phi^*_G` fold so the coarse system is a left-diagonal rescaling
    of the Petrov-Galerkin projection of the fine pencil (test
    :math:`\varphi^*\mathbf 1_G`, trial :math:`\varphi\mathbf 1_{G'}`).
    With the PLAIN condensed flux :math:`\Phi_G = \sum_{g\in G}\varphi_g`
    as the coarse carrier, the row-consistent constants are bilinear
    numerators over the row-uniform denominator
    :math:`\Phi^*_G \times (\text{trial fold})`:

    .. math::

        \Sigma_{t,C,G}
          = \frac{\sum_{g\in G}\varphi^*_g\Sigma_{t,g}\varphi_g}
                 {\Phi^*_G\,\Phi_G},
        \qquad
        \Sigma_{s,C}[G'\!\to\!G]
          = \frac{\sum_{g'\in G'}\sum_{g\in G}
                  \varphi^*_g\,\Sigma_s[g'\!\to\!g]\,\varphi_{g'}}
                 {\Phi^*_G\,\Phi_{G'}} .

    (The classical DIAGONAL-pair vector bilinear
    :math:`\langle\varphi^*\Sigma\varphi\rangle_G/\langle\varphi^*\varphi\rangle_G`
    is row-consistent with a DIFFERENT coarse carrier — the
    bilinear-weighted flux :math:`\tilde\Phi_G =
    \langle\varphi^*\varphi\rangle_G/\Phi^*_G`; each convention pair is
    exact with ITS OWN carrier, and mixing them is not (T6c).  Which pair
    is the classical prescription — and which carrier the re-solved coarse
    flux should MEAN — is the literature-reconciliation + user-checkpoint
    item.)  Then:

    * **T6a** — the condensed pencil's dominant eigenvalue equals the fine
      :math:`k` EXACTLY (rational arithmetic, zero error), and the
      condensed fission dyad is exactly rank-1 (the factored
      :math:`(\chi_C, \nu\Sigma_{f,C})` survives condensation — unlike the
      spatial axis, T3).
    * **T6b** — condensing with a PERTURBED spectrum pair
      (:math:`\varphi(1+\varepsilon h)`) leaves
      :math:`k_C(\varepsilon) - k = \mathcal O(\varepsilon^2)` for the
      bilinear convention vs :math:`\mathcal O(\varepsilon)` for the
      forward (flux-only) rules — the C2 comparative basis, exact on the
      energy axis.

    The DIAGONAL-pair vector bilinear (T1's spatial form) does NOT sit in a
    consistent row convention on the energy axis (T6c shows its k residual
    is nonzero at the true spectra): the energy-axis vector rule is the
    row-consistent one.  This is the convention the literature memo must
    reconcile (the classical "bilinear-weighted constants come with a
    bilinear-weighted coarse flux" subtlety).
    """
    print("=" * 68)
    print("T6. Energy condensation — exact projection (consistent denominators)")
    print("=" * 68)

    ng, blocks = 4, [(0, 1), (2, 3)]
    # Generic rational fine data: asymmetric scattering incl. upscatter.
    sig_t = [sp.Rational(x, 10) for x in (9, 8, 11, 12)]
    sig_s = sp.Matrix([
        [sp.Rational(3, 10), sp.Rational(1, 10), sp.Rational(1, 20), sp.Rational(1, 30)],
        [sp.Rational(1, 25), sp.Rational(4, 10), sp.Rational(1, 10), sp.Rational(1, 15)],
        [sp.Rational(1, 40), sp.Rational(1, 12), sp.Rational(5, 10), sp.Rational(1, 10)],
        [sp.Rational(1, 50), sp.Rational(1, 30), sp.Rational(1, 9), sp.Rational(6, 10)],
    ])  # [g_from, g_to]
    nsf = sp.Matrix([sp.Rational(x, 20) for x in (7, 5, 9, 4)])
    chi = sp.Matrix([sp.Rational(5, 10), sp.Rational(3, 10), sp.Rational(1, 10), sp.Rational(1, 10)])

    A = sp.diag(*sig_t) - sig_s.T          # row = sink group g
    phi_f = A.solve(chi)                    # fine flux (k-scaled out of the dyad)
    phis_f = A.T.solve(nsf)                 # fine adjoint
    k_fine = (nsf.T * phi_f)[0, 0]          # rank-1: k = νΣfᵀ A⁻¹ χ
    _require(k_fine > 0, "T6: fine k not positive")
    print(f"  fine k = {k_fine}")

    def condense(phi_v: sp.Matrix, phis_v: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
        """Row-consistent bilinear condensation → (Σt_C, Σs_C, νΣf_C, χ_C)."""
        Phi = [sum(phi_v[g] for g in blk) for blk in blocks]
        Phis = [sum(phis_v[g] for g in blk) for blk in blocks]
        sig_t_C = sp.Matrix([
            sum(phis_v[g] * sig_t[g] * phi_v[g] for g in blocks[G])
            / (Phis[G] * Phi[G])
            for G in range(2)
        ])
        sig_s_C = sp.zeros(2, 2)
        for Gf in range(2):
            for Gt in range(2):
                num = sum(
                    phis_v[g] * sig_s[gp, g] * phi_v[gp]
                    for gp in blocks[Gf] for g in blocks[Gt]
                )
                sig_s_C[Gf, Gt] = num / (Phis[Gt] * Phi[Gf])
        nsf_C = sp.Matrix([
            sum(nsf[gp] * phi_v[gp] for gp in blocks[Gf]) / Phi[Gf]
            for Gf in range(2)
        ])
        chi_raw = [
            sum(phis_v[g] * chi[g] for g in blocks[Gt]) / Phis[Gt]
            for Gt in range(2)
        ]
        s = chi_raw[0] + chi_raw[1]
        chi_C = sp.Matrix([chi_raw[0] / s, chi_raw[1] / s])   # simplex WLOG
        nsf_C = nsf_C * s                                     # rank-1 rescale freedom
        return sig_t_C, sig_s_C, nsf_C, chi_C

    def coarse_k(parts: tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]) -> sp.Expr:
        sig_t_C, sig_s_C, nsf_C, chi_C = parts
        A_C = sp.diag(*sig_t_C) - sig_s_C.T
        return (nsf_C.T * A_C.solve(chi_C))[0, 0]

    # T6a — exactness at the true spectra (+ rank-1 fission survives).
    k_C = coarse_k(condense(phi_f, phis_f))
    _is_zero(k_C - k_fine, "T6a: bilinear condensation is not exact at the true spectra")
    print("  ✓ T6a: k_C == k_fine EXACTLY (rational identity; rank-1 fission survives)")

    # T6c — the DIAGONAL-pair vector rule breaks row consistency on this axis.
    def condense_diag_sig_t(phi_v: sp.Matrix, phis_v: sp.Matrix) -> sp.Matrix:
        return sp.Matrix([
            sum(phis_v[g] * sig_t[g] * phi_v[g] for g in blocks[G])
            / sum(phis_v[g] * phi_v[g] for g in blocks[G])
            for G in range(2)
        ])

    parts = condense(phi_f, phis_f)
    diag_parts = (condense_diag_sig_t(phi_f, phis_f), parts[1], parts[2], parts[3])
    k_diag = coarse_k(diag_parts)
    _is_nonzero(k_diag - k_fine, "T6c: diagonal-pair Σt accidentally exact")
    print(f"  ✓ T6c: diagonal-pair Σt convention: k residual = {sp.nsimplify(k_diag - k_fine)} ≠ 0")

    # T6b — perturbed spectra: bilinear O(ε²) vs forward O(ε).
    eps = sp.Symbol("varepsilon")
    h = sp.Matrix([sp.Rational(1, 3), -sp.Rational(2, 5), sp.Rational(1, 7), -sp.Rational(1, 4)])
    hs = sp.Matrix([-sp.Rational(1, 6), sp.Rational(1, 5), -sp.Rational(2, 7), sp.Rational(1, 9)])
    phi_pert = sp.Matrix([phi_f[g] * (1 + eps * h[g]) for g in range(ng)])
    phis_pert = sp.Matrix([phis_f[g] * (1 + eps * hs[g]) for g in range(ng)])

    k_pert = coarse_k(condense(phi_pert, phis_pert))
    gap = sp.series(sp.together(k_pert - k_fine), eps, 0, 2).removeO()
    _is_zero(gap.coeff(eps, 0), "T6b: bilinear gap has an O(1) term")
    _is_zero(gap.coeff(eps, 1), "T6b: bilinear gap has a FIRST-order term")
    print("  ✓ T6b: bilinear condensation gap = O(ε²) — first-order term identically zero")

    def condense_forward(phi_v: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
        """Today's forward rules: flux-average source axis, marginalize sink + χ."""
        Phi = [sum(phi_v[g] for g in blk) for blk in blocks]
        sig_t_C = sp.Matrix([
            sum(sig_t[g] * phi_v[g] for g in blocks[G]) / Phi[G] for G in range(2)
        ])
        sig_s_C = sp.zeros(2, 2)
        for Gf in range(2):
            for Gt in range(2):
                num = sum(sig_s[gp, g] * phi_v[gp] for gp in blocks[Gf] for g in blocks[Gt])
                sig_s_C[Gf, Gt] = num / Phi[Gf]
        nsf_C = sp.Matrix([
            sum(nsf[gp] * phi_v[gp] for gp in blocks[Gf]) / Phi[Gf] for Gf in range(2)
        ])
        chi_C = sp.Matrix([sum(chi[g] for g in blocks[Gt]) for Gt in range(2)])
        return sig_t_C, sig_s_C, nsf_C, chi_C

    k_fwd_pert = coarse_k(condense_forward(phi_pert))
    gap_fwd = sp.series(sp.together(k_fwd_pert - k_fine), eps, 0, 2).removeO()
    # Forward condensation is ALSO exact at the true spectrum (ε=0) — the
    # classical property — but its gap is first-order in the spectrum error.
    _is_zero(gap_fwd.coeff(eps, 0), "T6b: forward gap has an O(1) term")
    _is_nonzero(gap_fwd.coeff(eps, 1), "T6b: forward gap first-order term vanished")
    print(f"  ✓ T6b: forward condensation gap = ({sp.nsimplify(gap_fwd.coeff(eps, 1))})·ε + O(ε²)\n")


# ═══════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════

def run_all() -> None:
    """Run every derivation (raises :class:`DerivationError` on any failure)."""
    derive_first_order_k_shift()
    derive_vector_channel_rule()
    derive_matrix_channel_rule()
    derive_fission_factored_rule()
    derive_balance_tradeoff()
    derive_forward_rule_first_order_error()
    derive_energy_condensation_exactness()
    print("ALL adjoint-weighted collapse theorems verified.")


if __name__ == "__main__":
    run_all()
