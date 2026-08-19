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

Angular pairing (T1b — USER-RULED implemented).  The collision worth
pairs the ANGULAR fluxes: the exact :math:`\Sigma_t` rule weights by
:math:`\rho_{i,g} = \sum_n w_n \psi^*_{i,g,n}\psi_{i,g,n}` (T1b — unique,
proved), of which the scalar :math:`\varphi^*\varphi` is the P0/isotropic
truncation (they coincide identically on isotropic shapes — also proved).
The user ruled (P6 open, option 2) that production implements the exact
angular rule for :math:`\Sigma_t` — both solutions carry :math:`\psi` —
rather than the classical scalar practice.  The scalar pair remains EXACT
for the isotropic operators (the P0 scattering source and fission), and
the moment-resolved refinement for anisotropic scattering orders
(:math:`\Sigma_{s,\ell}` pairing the :math:`\ell`-moments; :math:`\ell=0`
= T2's scalar pair; Parseval makes :math:`\rho` the all-moment sum) stays
a documented seam until an anisotropic-collapse consumer exists.

Programmatic documentation (USER-RULED).  Every production rule has a
grid-parameterized BUILDER (what the theorems verify) and a Sum-form
DISPLAY equation (:func:`collapse_rules` — what the B4 doc generator
renders into the theory pages); each theorem proves the two spellings
equal at the concrete grids, so the documented math and the verified math
cannot drift.

The energy axis (T6): condensation is pure projection
=====================================================

Energy condensation has NO streaming carve (streaming does not couple
groups): the coarse pencil is an EXACT left-diagonal rescaling of the
Petrov-Galerkin projection of the fine pencil (test
:math:`\varphi^*\cdot\mathbf 1_G`, trial :math:`\varphi\cdot\mathbf 1_{G'}`)
under a carrier-CONSISTENT convention — and the convention is **Bell &
Glasstone Ch. 6** (reconciled 2026-07-26, memo Source E): plain condensed
flux carrier (B&G (6.125)), flux-weighted-average adjoint carrier
:math:`\Psi^{\dagger}_G = \langle\varphi^*\varphi\rangle_G/\Phi_G`
((6.126)–(6.128)), diagonal-bilinear vector constants ((6.135)), per-block
sink×source matrix constants ((6.136)), factored fission with the simplex
rescale.  T6a proves the condensed :math:`k` exact at the true spectra;
T6r proves the row-scaling freedom (why :math:`k` alone cannot pin the
constants — the classical normalization does); T6c proves carrier
consistency is load-bearing (mixing folds loses exactness); T6b proves
the second-order-vs-first-order spectrum-error gap (B&G (6.90)) — the C2
basis, exact on the energy axis.

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
    "vector_bilinear_rule",
    "angular_sigma_t_rule",
    "matrix_per_pair_rule",
    "fission_nsf_mixed_fold_rule",
    "fission_chi_canonical_rule",
    "collapse_rules",
    "derive_first_order_k_shift",
    "derive_vector_channel_rule",
    "derive_angular_sigma_t_rule",
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
# The rule builders — the single source of the collapse-rule STRUCTURE.
#
# Each production rule exists in exactly two spellings, both here:
#
# * a grid-parameterized BUILDER (finite sums over explicit index sets) —
#   what the theorems verify and what any production cross-check calls;
# * a Sum-form DISPLAY equation over IndexedBase (``collapse_rules()``) —
#   what the doc generator renders into the theory pages (the user-ruled
#   programmatic algebra-of-record pipeline).
#
# The two spellings are welded by PROOF, not by convention: every theorem
# verifies ``display.subs(N, 2).doit() ≡ builder(concrete grids)`` before
# using the builder, so a drift between what the docs show and what the
# code implements cannot survive the pin suite.
# ═══════════════════════════════════════════════════════════════════════

def vector_bilinear_rule(cells, V_, phis_, phi_, sigma_, g: int) -> sp.Expr:
    r"""T1: :math:`\Sigma_{R,g} = \sum_i V_i\varphi^*_{i,g}\Sigma_{i,g}\varphi_{i,g}
    \big/ \sum_i V_i\varphi^*_{i,g}\varphi_{i,g}` — the pair weight φ*⊙φ."""
    num = sum(V_[i] * phis_[i][g] * sigma_[i][g] * phi_[i][g] for i in cells)
    den = sum(V_[i] * phis_[i][g] * phi_[i][g] for i in cells)
    return num / den


def angular_sigma_t_rule(cells, ordinates, V_, w_, psis_, psi_, sigma_, g: int) -> sp.Expr:
    r"""T1b: the EXACT collision rule — weight :math:`\rho_{i,g} = \sum_n
    w_n\psi^*_{i,g,n}\psi_{i,g,n}` (the angular pairing the collision worth
    actually carries; the scalar φ*⊙φ is its P0 truncation)."""
    rho = {
        i: sum(w_[n] * psis_[i][g][n] * psi_[i][g][n] for n in ordinates)
        for i in cells
    }
    num = sum(V_[i] * rho[i] * sigma_[i][g] for i in cells)
    den = sum(V_[i] * rho[i] for i in cells)
    return num / den


def matrix_per_pair_rule(cells, V_, phis_, phi_, sig_s_, gf: int, gt: int) -> sp.Expr:
    r"""T2: :math:`\Sigma_{s,R}[g'\!\to\!g] = \sum_i V_i\varphi^*_{i,g}
    \Sigma_{s,i}[g'\!\to\!g]\varphi_{i,g'} \big/ \sum_i V_i\varphi^*_{i,g}
    \varphi_{i,g'}` — sink adjoint × source flux, per pair."""
    num = sum(V_[i] * phis_[i][gt] * sig_s_[i][gf][gt] * phi_[i][gf] for i in cells)
    den = sum(V_[i] * phis_[i][gt] * phi_[i][gf] for i in cells)
    return num / den


def fission_nsf_mixed_fold_rule(
    cells, groups, V_, phis_, phi_, nsf_, chi_, chi_R_, gp: int,
) -> sp.Expr:
    r"""T3: the mixed-fold :math:`\nu\Sigma_{f,R,g'}` — numerator folded by the
    FINE emission importance :math:`\iota_i = \sum_g\varphi^*_{i,g}\chi_{i,g}`,
    denominator by the COLLAPSED :math:`\tilde\iota_i =
    \sum_g\varphi^*_{i,g}\chi_{R,g}` (exact total-worth for any simplex χ_R)."""
    iota = {i: sum(phis_[i][g] * chi_[i][g] for g in groups) for i in cells}
    iota_R = {i: sum(phis_[i][g] * chi_R_[g] for g in groups) for i in cells}
    num = sum(V_[i] * iota[i] * nsf_[i][gp] * phi_[i][gp] for i in cells)
    den = sum(V_[i] * iota_R[i] * phi_[i][gp] for i in cells)
    return num / den


def fission_chi_canonical_rule(cells, groups, V_, phis_, phi_, nsf_, chi_, g: int) -> sp.Expr:
    r"""T3: the canonical :math:`\chi_{R,g}` — the adjoint-weighted-emission
    convex average, weights :math:`q_i = V_i\,\iota_i\,p_i` (a simplex by
    construction; production-weighted at flat φ*)."""
    iota = {i: sum(phis_[i][gg] * chi_[i][gg] for gg in groups) for i in cells}
    p = {i: sum(nsf_[i][gp] * phi_[i][gp] for gp in groups) for i in cells}
    q = {i: V_[i] * iota[i] * p[i] for i in cells}
    return sum(q[i] * chi_[i][g] for i in cells) / sum(q[i] for i in cells)


# — the Sum-form display equations (the doc generator's input) —

_i = sp.Symbol("i", integer=True)
_n_idx = sp.Symbol("n", integer=True)
_g_sym = sp.Symbol("g", integer=True)
_gp_sym = sp.Symbol("g'", integer=True)
_NR = sp.Symbol("N_R", integer=True, positive=True)
_NW = sp.Symbol("N_Omega", integer=True, positive=True)

_Vb = sp.IndexedBase("V")
_wb = sp.IndexedBase("w")
_phib = sp.IndexedBase(r"\varphi")
_phisb = sp.IndexedBase(r"\varphi^{*}")
_psib = sp.IndexedBase(r"\psi")
_psisb = sp.IndexedBase(r"\psi^{*}")
_Sigb = sp.IndexedBase(r"\Sigma")
_Sigsb = sp.IndexedBase(r"\Sigma_{s}")
_nsfb = sp.IndexedBase(r"\nu\Sigma_{f}")
_chib = sp.IndexedBase(r"\chi")
_iotab = sp.IndexedBase(r"\iota")
_iotaRb = sp.IndexedBase(r"\tilde{\iota}")
_qb = sp.IndexedBase("q")


def collapse_rules() -> dict[str, sp.Eq]:
    r"""The display (Sum-form) collapse rules, keyed by channel class.

    These are what the theory pages RENDER (via the B4 generator).  Every
    entry is proof-welded to its builder inside the corresponding theorem
    (``display ≡ builder`` at the concrete grids), so the documented math
    and the verified math cannot drift apart.
    """
    S = sp.Sum
    vector = sp.Eq(
        sp.Symbol(r"\Sigma_{R,g}"),
        S(_Vb[_i] * _phisb[_i, _g_sym] * _Sigb[_i, _g_sym] * _phib[_i, _g_sym], (_i, 1, _NR))
        / S(_Vb[_i] * _phisb[_i, _g_sym] * _phib[_i, _g_sym], (_i, 1, _NR)),
    )
    rho = S(_wb[_n_idx] * _psisb[_i, _g_sym, _n_idx] * _psib[_i, _g_sym, _n_idx], (_n_idx, 1, _NW))
    angular = sp.Eq(
        sp.Symbol(r"\Sigma_{t,R,g}"),
        S(_Vb[_i] * rho * _Sigb[_i, _g_sym], (_i, 1, _NR))
        / S(_Vb[_i] * rho, (_i, 1, _NR)),
    )
    matrix = sp.Eq(
        sp.Symbol(r"\Sigma_{s,R}[g'\to g]"),
        S(_Vb[_i] * _phisb[_i, _g_sym] * _Sigsb[_i, _gp_sym, _g_sym] * _phib[_i, _gp_sym], (_i, 1, _NR))
        / S(_Vb[_i] * _phisb[_i, _g_sym] * _phib[_i, _gp_sym], (_i, 1, _NR)),
    )
    nsf_eq = sp.Eq(
        sp.Symbol(r"\nu\Sigma_{f,R,g'}"),
        S(_Vb[_i] * _iotab[_i] * _nsfb[_i, _gp_sym] * _phib[_i, _gp_sym], (_i, 1, _NR))
        / S(_Vb[_i] * _iotaRb[_i] * _phib[_i, _gp_sym], (_i, 1, _NR)),
    )
    chi_eq = sp.Eq(
        sp.Symbol(r"\chi_{R,g}"),
        S(_qb[_i] * _chib[_i, _g_sym], (_i, 1, _NR)) / S(_qb[_i], (_i, 1, _NR)),
    )
    return {
        "vector": vector,
        "angular_sigma_t": angular,
        "matrix_per_pair": matrix,
        "fission_nsf_mixed_fold": nsf_eq,
        "fission_chi_canonical": chi_eq,
    }


def _display_matches_builder(key: str, builder_expr: sp.Expr, subs_map: dict) -> None:
    """Prove ``collapse_rules()[key]`` expands (N_R→2, indices→grids) to the builder."""
    display = collapse_rules()[key].rhs
    expanded = display.subs({_NR: _N, _NW: _N}).doit().subs(subs_map)
    _is_zero(expanded - builder_expr, f"display '{key}' drifted from its builder")


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
        return vector_bilinear_rule(_CELLS, V, phis, phi, sigma, g)

    # Weld the display (Sum-form) spelling to the builder before using it.
    _display_matches_builder(
        "vector", bilinear(0),
        {
            **{_Vb[i + 1]: V[i] for i in _CELLS},
            **{_phisb[i + 1, _g_sym]: phis[i][0] for i in _CELLS},
            **{_phib[i + 1, _g_sym]: phi[i][0] for i in _CELLS},
            **{_Sigb[i + 1, _g_sym]: sigma[i][0] for i in _CELLS},
        },
    )
    print("  ✓ display Sum-form ≡ builder (the doc-generator weld)")

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
# T1b — the collision channel: the EXACT angular pairing (user-ruled)
# ═══════════════════════════════════════════════════════════════════════

def derive_angular_sigma_t_rule() -> None:
    r"""T1b: the collision worth pairs ANGULARLY — the exact Σt rule weights by
    :math:`\rho_{i,g} = \sum_n w_n\psi^*_{i,g,n}\psi_{i,g,n}`.

    The collision term of the pencil acts on the full angular flux, so its
    worth is :math:`\sum_{i,g,n} V_i w_n \psi^*_{i,g,n}(\Sigma_{R,g} -
    \Sigma_{i,g})\psi_{i,g,n}` — zeroed uniquely by the angular-product
    weight :math:`\rho` (proved), NOT by the scalar pair :math:`\varphi^*
    \varphi` (its P0 truncation; generic anisotropic counterexample).  On
    isotropic angular shapes (:math:`\psi = \varphi/W`, :math:`\psi^* =
    \varphi^*/W`) the two rules coincide identically — the classical
    scalar prescription is exactly the isotropic limit.

    USER RULING (P6 open, option 2): production implements THIS rule for
    :math:`\Sigma_t` — both solutions carry :math:`\psi` — rather than the
    P0 truncation.  The moment-resolved refinement for the anisotropic
    scattering orders (:math:`\Sigma_{s,\ell}` pairing the
    :math:`\ell`-moments :math:`\varphi^*_{\ell m}\varphi_{\ell m}`;
    :math:`\ell = 0` is exactly T2's scalar pair, and Parseval makes the
    :math:`\rho`-weight the all-moment sum) stays a documented seam until
    an anisotropic-collapse consumer exists.
    """
    print("=" * 68)
    print("T1b. Collision channel — the exact angular pairing ρ = Σ_n w ψ*ψ")
    print("=" * 68)

    sigma = [[sp.Symbol(f"sigma_{i}{g}", positive=True) for g in _GROUPS] for i in _CELLS]
    w = [sp.Symbol(f"w_{n}", positive=True) for n in range(_N)]
    psi = [
        [[sp.Symbol(f"psi_{i}{g}{n}", positive=True) for n in range(_N)] for g in _GROUPS]
        for i in _CELLS
    ]
    psis = [
        [[sp.Symbol(f"psis_{i}{g}{n}", positive=True) for n in range(_N)] for g in _GROUPS]
        for i in _CELLS
    ]
    ordinates = range(_N)

    def angular(g: int) -> sp.Expr:
        return angular_sigma_t_rule(_CELLS, ordinates, V, w, psis, psi, sigma, g)

    _display_matches_builder(
        "angular_sigma_t", angular(0),
        {
            **{_Vb[i + 1]: V[i] for i in _CELLS},
            **{_wb[n + 1]: w[n] for n in ordinates},
            **{_psisb[i + 1, _g_sym, n + 1]: psis[i][0][n] for i in _CELLS for n in ordinates},
            **{_psib[i + 1, _g_sym, n + 1]: psi[i][0][n] for i in _CELLS for n in ordinates},
            **{_Sigb[i + 1, _g_sym]: sigma[i][0] for i in _CELLS},
        },
    )
    print("  ✓ display Sum-form ≡ builder (the doc-generator weld)")

    def worth(rule: Callable[[int], sp.Expr], g: int) -> sp.Expr:
        return sum(
            V[i] * w[n] * psis[i][g][n] * (rule(g) - sigma[i][g]) * psi[i][g][n]
            for i in _CELLS for n in ordinates
        )

    for g in _GROUPS:
        _is_zero(worth(angular, g), f"T1b: angular rule leaves collision worth in group {g}")
    print("  ✓ ρ-weighted rule zeroes the ANGULAR collision worth identically (all g)")

    # Uniqueness.
    sig_R = sp.Symbol("sigma_R", positive=True)
    w0 = sum(
        V[i] * w[n] * psis[i][0][n] * (sig_R - sigma[i][0]) * psi[i][0][n]
        for i in _CELLS for n in ordinates
    )
    solved = sp.solve(sp.Eq(w0, 0), sig_R)
    _require(len(solved) == 1, "T1b: worth-zeroing rule is not unique")
    _is_zero(solved[0] - angular(0), "T1b: solved rule differs from the ρ-weighted builder")
    print("  ✓ the angular worth-zeroing rule is UNIQUE (= the ρ-weighted builder)")

    # Isotropic degeneracy: ψ = φ/W, ψ* = φ*/W ⇒ angular rule ≡ scalar bilinear.
    W_tot = sum(w[n] for n in ordinates)
    iso = {}
    for i in _CELLS:
        for g in _GROUPS:
            for n in ordinates:
                iso[psi[i][g][n]] = phi[i][g] / W_tot
                iso[psis[i][g][n]] = phis[i][g] / W_tot
    scalar = vector_bilinear_rule(_CELLS, V, phis, phi, sigma, 0)
    _is_zero(
        sp.simplify(angular(0).subs(iso) - scalar),
        "T1b: isotropic degenerate is not the scalar bilinear rule",
    )
    print("  ✓ isotropic ψ, ψ* ⇒ angular rule ≡ scalar φ*⊙φ rule (the P0 limit)")

    # The scalar rule does NOT zero the angular worth on anisotropic shapes.
    subs = _generic_rationals()
    primes = iter([47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131])
    for i in _CELLS:
        for g in _GROUPS:
            subs[sigma[i][g]] = sp.Rational(7 * i + 3 * g + 2, 9)
            for n in ordinates:
                subs[psi[i][g][n]] = sp.Rational(next(primes), 8)
                subs[psis[i][g][n]] = sp.Rational(next(primes), 11)
    subs[w[0]], subs[w[1]] = sp.Rational(1, 2), sp.Rational(1, 2)
    # The scalar rule evaluated with the CONSISTENT scalar moments of ψ, ψ*.
    phi_of_psi = {
        phi[i][g]: sum(w[n] * psi[i][g][n] for n in ordinates) for i in _CELLS for g in _GROUPS
    }
    phis_of_psi = {
        phis[i][g]: sum(w[n] * psis[i][g][n] for n in ordinates) for i in _CELLS for g in _GROUPS
    }
    scalar_rule_val = scalar.subs(phi_of_psi).subs(phis_of_psi)
    residual = worth(lambda g: scalar_rule_val, 0).subs(subs)
    _is_nonzero(residual, "T1b: scalar rule accidentally zeroed the angular worth")
    print(f"  ✓ scalar φ*⊙φ rule leaves angular worth = {sp.nsimplify(residual)} ≠ 0 (anisotropic)\n")


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
        return matrix_per_pair_rule(_CELLS, V, phis, phi, sig_s, gf, gt)

    _display_matches_builder(
        "matrix_per_pair", per_pair(0, 1),
        {
            **{_Vb[i + 1]: V[i] for i in _CELLS},
            **{_phisb[i + 1, _g_sym]: phis[i][1] for i in _CELLS},
            **{_phib[i + 1, _gp_sym]: phi[i][0] for i in _CELLS},
            **{_Sigsb[i + 1, _gp_sym, _g_sym]: sig_s[i][0][1] for i in _CELLS},
        },
    )
    print("  ✓ display Sum-form ≡ builder (the doc-generator weld)")

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
        return fission_nsf_mixed_fold_rule(
            _CELLS, _GROUPS, V, phis, phi, nsf, chi, chi_R, gp,
        )

    _display_matches_builder(
        "fission_nsf_mixed_fold", nsf_R(0),
        {
            **{_Vb[i + 1]: V[i] for i in _CELLS},
            **{_iotab[i + 1]: iota[i] for i in _CELLS},
            **{_iotaRb[i + 1]: iota_R[i] for i in _CELLS},
            **{_nsfb[i + 1, _gp_sym]: nsf[i][0] for i in _CELLS},
            **{_phib[i + 1, _gp_sym]: phi[i][0] for i in _CELLS},
        },
    )
    print("  ✓ display Sum-form ≡ builder (νΣf mixed-fold; the doc-generator weld)")

    total_worth = sum(
        V[i] * phis[i][g] * (chi_R[g] * nsf_R(gp) - chi[i][g] * nsf[i][gp]) * phi[i][gp]
        for i in _CELLS for g in _GROUPS for gp in _GROUPS
    )
    _is_zero(total_worth, "T3: mixed-fold rule leaves total fission worth")
    print("  ✓ W_F ≡ 0 for ANY simplex χ_R (t left free) — the χ slot is exactness-free")

    # Canonical χ_R: adjoint-weighted-emission convex average (q_i = V_i ι_i p_i).
    q = [V[i] * iota[i] * p[i] for i in _CELLS]
    chi_canonical = [
        fission_chi_canonical_rule(_CELLS, _GROUPS, V, phis, phi, nsf, chi, g)
        for g in _GROUPS
    ]
    _display_matches_builder(
        "fission_chi_canonical", chi_canonical[0],
        {
            **{_qb[i + 1]: q[i] for i in _CELLS},
            **{_chib[i + 1, _g_sym]: chi[i][0] for i in _CELLS},
        },
    )
    print("  ✓ display Sum-form ≡ builder (canonical χ; the doc-generator weld)")
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

    # The CANONICAL builders — never a private re-spelling (the enforcer's
    # weld discipline: a theorem constraining a rule must consume the same
    # single source the production code mirrors, else it silently pins a
    # stale copy the day the rule is refined).
    def t1(channel, g: int) -> sp.Expr:
        return vector_bilinear_rule(_CELLS, V, phis, phi, channel, g)

    def t2(gf: int, gt: int) -> sp.Expr:
        return matrix_per_pair_rule(_CELLS, V, phis, phi, sig_s, gf, gt)

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
    r"""T6: bilinear condensation in the B&G convention reproduces k EXACTLY.

    0-D ∞-medium pencil, 4 fine groups → 2 coarse (blocks {0,1}, {2,3}),
    generic rational data (asymmetric :math:`\Sigma_s` with upscatter,
    :math:`\chi\not\parallel\nu\Sigma_f`).  Fine pair :math:`\varphi =
    A^{-1}\chi`, :math:`\varphi^* = A^{-\mathsf T}\nu\Sigma_f`, :math:`k =
    \nu\Sigma_f^{\mathsf T} A^{-1}\chi` (rank-1 fission).

    **The convention is Bell & Glasstone Ch. 6** (reconciled 2026-07-26;
    `.claude/plans/archive/p6_literature_memo.md` Source E): coarse FLUX carrier =
    the plain condensed flux :math:`\Phi_G = \sum_{g\in G}\varphi_g`
    (B&G (6.125)); coarse ADJOINT carrier = the flux-weighted average
    :math:`\Psi^{\dagger}_G = \langle\varphi^*\varphi\rangle_G/\Phi_G`
    (B&G (6.126)–(6.128)); constants = the bilinear group constants
    (B&G (6.135)/(6.136)), spelled discrete and unnormalized as

    .. math::

        \Sigma_{C,G}
          = \frac{\langle\varphi^*\Sigma\varphi\rangle_G}
                 {\langle\varphi^*\varphi\rangle_G},
        \qquad
        \Sigma_{s,C}[G'\!\to\!G]
          = \frac{\sum_{g'\in G'}\sum_{g\in G}
                  \varphi^*_g\,\Sigma_s[g'\!\to\!g]\,\varphi_{g'}}
                 {\Phi_{G'}\;\Psi^{\dagger}_G},

    fission FACTORED (B&G (4.38)+(6.136), separable kernel):
    :math:`\nu\Sigma_{f,C,G'} = \langle\nu\Sigma_f\varphi\rangle_{G'}/
    \Phi_{G'}` (flux-weighted) and :math:`\chi^{\dagger}_G =
    \langle\chi\varphi^*\rangle_G/\Psi^{\dagger}_G`, with the rank-1
    rescale :math:`s = \sum_G\chi^{\dagger}_G` moved into
    :math:`\nu\Sigma_{f,C}` so :math:`\chi_C` is a simplex (the Mixture
    law; :math:`k`-invariant).  Every channel of coarse row :math:`G` then
    divides by the same fold :math:`\Psi^{\dagger}_G`, so the coarse system
    is a left-diagonal rescaling of the Petrov-Galerkin projection of the
    fine pencil (test :math:`\varphi^*\mathbf 1_G`, trial
    :math:`\varphi\mathbf 1_{G'}`) evaluated at the plain carrier:

    * **T6a** — the condensed pencil's :math:`k` equals the fine :math:`k`
      EXACTLY (rational identity; the factored fission survives
      condensation — unlike the spatial axis, T3).
    * **T6r** — row-scaling freedom: rescaling every channel of a coarse
      row by a common factor preserves the pencil's spectrum, so
      :math:`k`-exactness alone CANNOT pin the constants' values — an
      entire family of row-consistent conventions is exact.  What pins
      them is the classical carrier normalization (B&G (6.125)/(6.126));
      the Mode-12 lens, one level up.
    * **T6c** — MIXING carriers (the diagonal bilinear vector rule paired
      with a plain-sum adjoint fold in the matrix denominators) is NOT
      row-consistent and loses exactness — carrier consistency is
      load-bearing, not cosmetic.
    * **T6b** — condensing with a PERTURBED spectrum pair leaves
      :math:`k_C(\varepsilon) - k = \mathcal O(\varepsilon^2)` for the
      bilinear convention (B&G (6.90): the flux-only error term is *first*
      order) vs :math:`\mathcal O(\varepsilon)` for the forward rules —
      the C2 comparative basis, exact on the energy axis.
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
        """B&G-convention bilinear condensation → (Σt_C, Σs_C, νΣf_C, χ_C)."""
        Phi = [sum(phi_v[g] for g in blk) for blk in blocks]
        pair = [sum(phis_v[g] * phi_v[g] for g in blk) for blk in blocks]
        psi_dag = [pair[G] / Phi[G] for G in range(2)]        # Ψ†_G (B&G 6.126-6.128)
        sig_t_C = sp.Matrix([
            sum(phis_v[g] * sig_t[g] * phi_v[g] for g in blocks[G]) / pair[G]
            for G in range(2)                                  # B&G (6.135)
        ])
        sig_s_C = sp.zeros(2, 2)
        for Gf in range(2):
            for Gt in range(2):
                num = sum(
                    phis_v[g] * sig_s[gp, g] * phi_v[gp]
                    for gp in blocks[Gf] for g in blocks[Gt]
                )
                sig_s_C[Gf, Gt] = num / (Phi[Gf] * psi_dag[Gt])  # B&G (6.136)
        nsf_C = sp.Matrix([
            sum(nsf[gp] * phi_v[gp] for gp in blocks[Gf]) / Phi[Gf]
            for Gf in range(2)                                 # flux-weighted (B&G 4.38+6.136)
        ])
        chi_dag = [
            sum(phis_v[g] * chi[g] for g in blocks[Gt]) / psi_dag[Gt]
            for Gt in range(2)                                 # adjoint-contracted emission
        ]
        s = chi_dag[0] + chi_dag[1]
        chi_C = sp.Matrix([chi_dag[0] / s, chi_dag[1] / s])   # simplex (Mixture law)
        nsf_C = nsf_C * s                                     # rank-1 rescale freedom
        return sig_t_C, sig_s_C, nsf_C, chi_C

    def coarse_k(parts: tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]) -> sp.Expr:
        sig_t_C, sig_s_C, nsf_C, chi_C = parts
        A_C = sp.diag(*sig_t_C) - sig_s_C.T
        return (nsf_C.T * A_C.solve(chi_C))[0, 0]

    # T6a — exactness at the true spectra (+ rank-1 fission survives).
    parts = condense(phi_f, phis_f)
    k_C = coarse_k(parts)
    _is_zero(k_C - k_fine, "T6a: B&G bilinear condensation is not exact at the true spectra")
    print("  ✓ T6a: k_C == k_fine EXACTLY (B&G convention; rank-1 fission survives)")

    # T6r — row-scaling freedom: a common per-row rescale of every channel
    # preserves the pencil spectrum, so k cannot pin the constants' values;
    # the classical carrier normalization does.
    c_row = sp.Rational(7, 3)
    sig_t_r = sp.Matrix([parts[0][0] * c_row, parts[0][1]])
    sig_s_r = parts[1].copy()
    for Gf in range(2):
        sig_s_r[Gf, 0] = sig_s_r[Gf, 0] * c_row               # sink-row 0 rescale
    chi_r = sp.Matrix([parts[3][0] * c_row, parts[3][1]])     # (simplex waived here)
    k_r = coarse_k((sig_t_r, sig_s_r, parts[2], chi_r))
    _is_zero(k_r - k_fine, "T6r: row rescale moved k (it must not)")
    print("  ✓ T6r: a common per-row rescale leaves k EXACT — k cannot pin the constants")

    # T6c — MIXED carriers lose exactness: B&G diagonal-bilinear Σt paired
    # with matrix denominators folded by the PLAIN adjoint sum (Φ*_G) rather
    # than the flux-weighted-average carrier Ψ†_G.
    Phi_blocks = [sum(phi_f[g] for g in blk) for blk in blocks]
    Phis_plain = [sum(phis_f[g] for g in blk) for blk in blocks]
    sig_s_mixed = sp.zeros(2, 2)
    for Gf in range(2):
        for Gt in range(2):
            num = sum(
                phis_f[g] * sig_s[gp, g] * phi_f[gp]
                for gp in blocks[Gf] for g in blocks[Gt]
            )
            sig_s_mixed[Gf, Gt] = num / (Phi_blocks[Gf] * Phis_plain[Gt])
    k_mixed = coarse_k((parts[0], sig_s_mixed, parts[2], parts[3]))
    _is_nonzero(k_mixed - k_fine, "T6c: mixed-carrier convention accidentally exact")
    print(f"  ✓ T6c: mixed carriers (plain-sum adjoint fold): k residual = "
          f"{sp.nsimplify(k_mixed - k_fine)} ≠ 0")

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
    derive_angular_sigma_t_rule()
    derive_matrix_channel_rule()
    derive_fission_factored_rule()
    derive_balance_tradeoff()
    derive_forward_rule_first_order_error()
    derive_energy_condensation_exactness()
    print("ALL adjoint-weighted collapse theorems verified.")


if __name__ == "__main__":
    run_all()
