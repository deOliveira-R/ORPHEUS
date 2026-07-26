r"""Consistent DSA for slab weighted-diamond SN — the four-step derivation of record.

The algebra of record for Phase 3a of the DSA campaign (issue #2): executes
Larsen's four-step procedure [Larsen1982a]_ **symbolically** on the slab
weighted-diamond (WD) discrete-ordinates equations and proves, rather than
transcribes, the consistent low-order (synthetic diffusion) operator:

* the interior tridiagonal row in the **edge** scalar-flux corrections
  :math:`f_0` (Larsen Eq. (27)) with the coefficient set (23a–f),
* the cell-average update relations (28a/28b),
* the Marshak / reflecting boundary rows (38)/(39) closed via the one-sided
  relations (25)/(26),
* and the object-level identity behind the recipe: the four-step output IS
  the two-moment (:math:`\ell \le 1`) restriction of the promoted-minus-
  original discrete system followed by a Schur elimination of everything
  but the edge :math:`f_0` — "consistent" = *derived by moment reduction*,
  never *discretize-the-reduced* (the frame-analysis Q2 verdict,
  ``.claude/plans/dsa_rp_frame_analysis.md``).

Derivation discipline
---------------------

**Quadrature moments stay symbolic.** The S4-symmetric quadrature is carried
as four positive symbols :math:`(\mu_a, \mu_b, \omega_a, \omega_b)` with the
:math:`\pm` pairing imposed structurally. Angular sums reduce to polynomials
in these symbols, and the classic constants emerge as *properties*, not
transcribed literals — each at its true mechanism:

* the :math:`1/3` and :math:`2/3` in the P1 balance moment are the LEGENDRE
  RECURSION coefficients of :math:`\mu^2 = (2P_2 + P_0)/3` — polynomial
  identities, quadrature-independent (:func:`prove_p_recursion_lemma`);
* the :math:`1/3` in :math:`D = 1/[3(\sigma_T-\sigma_{S1})]` and the
  :math:`3` in the closure coupling :math:`(3/2)\rho` are the QUADRATURE
  moment :math:`W_2 = \sum_m \mu_m^2\omega_m` (exact 1/3 for any rule
  integrating :math:`\mu^2` exactly under :math:`W_0 = 1` — Gauss-Legendre
  from :math:`N \ge 2`);
* the Marshak half-range coefficients are :math:`\gamma_N =
  \sum_{\mu_m>0}\mu_m\omega_m` and :math:`3W_2^+` (NOT the continuum 1/4
  and 1/2; substituting continuum values at small :math:`N` is a quiet
  consistency break — [AdamsLarsen2002]_ §III.B item M2-4).

This dodges the two recorded transcription hazards outright (the ⅓→½ OCR
trap and the small-N γ substitution; literature memo §6.1.3).

**The slot assignment is the method.** Step 3 promotes ONLY the explicit
:math:`\ell \le 1` moments; the :math:`\phi_2` terms and the closure
residual functionals :math:`L_0[\gamma\psi]`, :math:`L_0[\beta\psi]` are
LAGGED — carried here as opaque symbols precisely so the promotion cannot
touch them (expanding them in moment coordinates would silently promote
part of a lagged term — a derivation bug this module's first draft caught).
Their exactness as rewrites of the raw angular reductions is proven
separately (:func:`derive_moment_equations`), for EVERY symmetric
quadrature, via the annihilation identities.

**Normalization.** This module works in Larsen's convention
:math:`W_0 = \sum_m \omega_m = 1` (:math:`\phi` = angular *average*),
imposed only where a printed target assumes it. The ORPHEUS slab
quadrature carries raw Gauss-Legendre weights :math:`\sum_m w_m = 2`
(``orpheus/numerics/quadrature/rules_1d.py``); the map
:math:`\omega_m = w_m/2` is applied ONCE at the numeric boundary of this
module (the production ties in ``tests/derivations/test_dsa_rules.py``)
with the weight sums asserted numerically there.

**Transcription reference.** [Larsen1982a]_ is the sole transcription
source for the target forms; [Alcouffe1977]_ is a cross-check ONLY — its
printed Eqs. (17) and (23) carry sign errata (memo §1.5), so nothing here
transcribes from them.

Verified by ``tests/derivations/test_dsa_rules.py`` (foundation level).
Every check raises :class:`DerivationError` — never a bare ``assert`` — so
the proofs hold under ``python -O`` (vv-principles Mode 8).

References
----------
.. [Larsen1982a] E. W. Larsen, "Unconditionally Stable Diffusion-Synthetic
   Acceleration Methods for the Slab Geometry Discrete Ordinates Equations.
   Part I: Theory," *Nucl. Sci. Eng.* **82**, 47–63 (1982). Bare equation
   numbers in this module refer to this paper.
.. [McCoyLarsen1982] D. R. McCoy and E. W. Larsen, Part II: Numerical
   Results, *Nucl. Sci. Eng.* **82**, 64–70 (1982).
.. [Alcouffe1977] R. E. Alcouffe, "Diffusion Synthetic Acceleration Methods
   for the Diamond-Differenced Discrete-Ordinates Equations," *Nucl. Sci.
   Eng.* **64**, 344–355 (1977).
.. [AdamsLarsen2002] M. L. Adams and E. W. Larsen, "Fast iterative methods
   for discrete-ordinates particle transport calculations," *Prog. Nucl.
   Energy* **40**, 3–159 (2002).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import sympy as sp


class DerivationError(RuntimeError):
    """A symbolic identity this module claims to prove failed to hold."""


def _require(condition: bool, message: str) -> None:
    """Mode-8-proof check: raises (never a bare ``assert``)."""
    if not condition:
        raise DerivationError(message)


def _is_zero(expr: sp.Expr) -> bool:
    """Whether ``expr`` simplifies to exactly zero."""
    return sp.simplify(sp.expand(expr)) == 0


def _rational_zero(expr: sp.Expr) -> bool:
    """Zero test for (possibly large) RATIONAL expressions.

    ``cancel`` + numerator-expansion — exact and far cheaper than
    ``simplify`` on multivariate rational functions.
    """
    num, _den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.expand(num) == 0


# ═══════════════════════════════════════════════════════════════════════
# The symbolic S4-symmetric slab quadrature
# ═══════════════════════════════════════════════════════════════════════

#: Positive polar cosines and weights of the symmetric rule. The ordinate
#: list is (+μ_a, +μ_b, −μ_b, −μ_a) with weights (ω_a, ω_b, ω_b, ω_a) —
#: the ± pairing is imposed STRUCTURALLY, so odd moments vanish
#: identically without any quadrature identity being assumed.
MU_A, MU_B = sp.symbols("mu_a mu_b", positive=True)
W_A, W_B = sp.symbols("omega_a omega_b", positive=True)

ORDINATES: tuple[tuple[sp.Expr, sp.Expr], ...] = (
    (MU_A, W_A),
    (MU_B, W_B),
    (-MU_B, W_B),
    (-MU_A, W_A),
)
N_ORD = len(ORDINATES)


def moment_sum(integrand) -> sp.Expr:
    r""":math:`\sum_m \omega_m\,\mathrm{integrand}(\mu_m)` over the rule."""
    return sp.expand(sum(w * integrand(mu) for mu, w in ORDINATES))


def half_range_sum(integrand, *, positive: bool) -> sp.Expr:
    r"""The same reduction restricted to :math:`\mu_m > 0` (or :math:`< 0`)."""
    return sp.expand(
        sum(
            w * integrand(mu)
            for mu, w in ORDINATES
            if (mu.is_positive if positive else mu.is_negative)
        )
    )


#: The named quadrature moments — the derivation's ONLY interface to the
#: rule. W1 = W3 = 0 hold structurally (± pairing); the rest stay symbolic.
W0 = moment_sum(lambda mu: 1)
W2 = moment_sum(lambda mu: mu**2)
GAMMA_HALF = half_range_sum(lambda mu: mu, positive=True)  # γ_N
W2_HALF = half_range_sum(lambda mu: mu**2, positive=True)  # W2⁺

#: Substitution realizing Larsen's normalization W0 = 1 (solved for ω_b).
#: Applied ONLY at printed-target comparisons — every structural identity
#: in this module holds without it.
LARSEN_NORMALIZATION = {W_B: sp.Rational(1, 2) - W_A}

#: Larsen's full printed convention: W0 = 1 AND the μ²-exactness W2 = 1/3
#: (what Gauss-Legendre delivers from N ≥ 2). The two constraints are
#: LINEAR in the weights, so solving for (ω_a, ω_b) keeps the whole
#: comparison RATIONAL in (μ_a, μ_b) — no radicals, and full theorem
#: strength on the entire constraint variety (μ_a, μ_b stay free).
PRINTED_CONVENTION = {
    W_A: (sp.Rational(1, 6) - MU_B**2 / 2) / (MU_A**2 - MU_B**2),
    W_B: (MU_A**2 / 2 - sp.Rational(1, 6)) / (MU_A**2 - MU_B**2),
}


def _to_printed_convention(expr: sp.Expr) -> sp.Expr:
    """Push an expression into Larsen's printed convention (W0=1, W2=1/3)."""
    return sp.cancel(sp.together(sp.expand(expr).subs(PRINTED_CONVENTION)))


# ═══════════════════════════════════════════════════════════════════════
# The WD closure weights and their annihilation identities (Larsen (14))
# ═══════════════════════════════════════════════════════════════════════

# Generic-cell data symbols (translation-generic derivation: cell "i" and
# its right neighbor carry independent symbols).
H_I, H_IP1 = sp.symbols("h_i h_ip1", positive=True)
SIG_T_I, SIG_T_IP1 = sp.symbols("sigma_T_i sigma_T_ip1", positive=True)
SIG_S0_I, SIG_S0_IP1 = sp.symbols("sigma_S0_i sigma_S0_ip1", nonnegative=True)
SIG_S1_I, SIG_S1_IP1 = sp.symbols("sigma_S1_i sigma_S1_ip1", real=True)
#: Fixed-weight WD closure magnitude a (α_m = a·sgn(μ_m); DD: a → 0).
A_WD_I, A_WD_IP1 = sp.symbols("a_i a_ip1", nonnegative=True)


def wd_alpha(a_cell: sp.Expr, mu: sp.Expr) -> sp.Expr:
    r""":math:`\alpha_{mi} = a_i\,\mathrm{sgn}(\mu_m)` — a Larsen (12a–c)
    member. DD is :math:`a_i = 0`; step is :math:`a_i = 1`."""
    return a_cell * sp.sign(mu)


def rho_of(a_cell: sp.Expr) -> sp.Expr:
    r""":math:`\rho_i = L_1\alpha_i = \sum_m \mu_m\alpha_{mi}\omega_m` (14a)."""
    return moment_sum(lambda mu: mu * wd_alpha(a_cell, mu))


def derive_closure_weight_identities() -> None:
    r"""Prove the annihilation set that closes step 2 (Larsen (14)–(15)).

    With :math:`\gamma_m = \alpha_m - \mu_m\rho/W_2` and :math:`\beta_m =
    \mu_m\alpha_m - \rho/W_0`,

    .. math::

        L_0\gamma = L_1\gamma = 0, \qquad L_0\beta = L_1\beta = 0

    hold for EVERY symmetric quadrature — Larsen's (14b) "3" is
    :math:`1/W_2` in disguise, recovered verbatim at :math:`W_2 = 1/3`.
    This is the property behind the exactness of the explicit/opaque
    split in :func:`derive_moment_equations`.
    """
    a = A_WD_I
    rho = rho_of(a)
    _require(
        _is_zero(rho - 2 * a * (MU_A * W_A + MU_B * W_B)),
        "ρ = L1α must reduce to the half-range sum 2a·Σ_{μ>0}μω",
    )

    def gamma(mu: sp.Expr) -> sp.Expr:
        return wd_alpha(a, mu) - mu * rho / W2

    def beta(mu: sp.Expr) -> sp.Expr:
        return mu * wd_alpha(a, mu) - rho / W0

    checks = {
        "L0_gamma": moment_sum(gamma),
        "L1_gamma": moment_sum(lambda mu: mu * gamma(mu)),
        "L0_beta": moment_sum(beta),
        "L1_beta": moment_sum(lambda mu: mu * beta(mu)),
    }
    for name, expr in checks.items():
        _require(_is_zero(expr), f"annihilation identity {name} must vanish")


def prove_p_recursion_lemma() -> None:
    r"""The P1-balance ⅓ mechanism: :math:`L_1[\mu v] = \tfrac{2}{3}L_2 v
    + \tfrac{1}{3}L_0 v` for ANY per-ordinate vector :math:`v`.

    Pure polynomial identity (:math:`\mu^2 = (2P_2(\mu) + P_0(\mu))/3`)
    threaded through the moment definitions — quadrature-INDEPENDENT,
    unlike the :math:`W_2` mechanism. Keeping the two mechanisms distinct
    is the honest spelling of "the ⅓" (memo §6.1.3 watch item).
    """
    v = sp.symbols(f"v0:{N_ORD}")
    l1_mu_v = sum(mu * mu * w * v[m] for m, (mu, w) in enumerate(ORDINATES))
    l2 = sum(
        sp.legendre(2, mu) * w * v[m] for m, (mu, w) in enumerate(ORDINATES)
    )
    l0 = sum(w * v[m] for m, (mu, w) in enumerate(ORDINATES))
    _require(
        _is_zero(l1_mu_v - (sp.Rational(2, 3) * l2 + sp.Rational(1, 3) * l0)),
        "L1[μv] = (2/3)L2[v] + (1/3)L0[v] must hold identically",
    )


# ═══════════════════════════════════════════════════════════════════════
# Steps 1–2 — the four moment equations, derived and identity-verified
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class CellSymbols:
    r"""The moment-level symbols of one cell.

    Explicit unknowns: :math:`(\phi_0, \phi_1)` at ``left`` / ``right``
    edges and the ``bar`` average, plus the lagged second moments
    ``phi2_left`` / ``phi2_right``. Opaque lagged closure functionals:
    ``g_left`` / ``g_right`` (:math:`L_0[\gamma\psi]` at each edge) and
    ``b_left`` / ``b_right`` (:math:`L_0[\beta\psi]`). Scattering-source
    iterate moments ``phi0_lag`` / ``phi1_lag`` (:math:`\phi_{ni}^l`) and
    external-source moments ``s0`` / ``s1``.
    """

    tag: str

    def __post_init__(self) -> None:
        s = lambda name: sp.Symbol(f"{name}_{self.tag}")  # noqa: E731
        for loc in ("left", "right", "bar"):
            object.__setattr__(self, loc, (s(f"phi0_{loc}"), s(f"phi1_{loc}")))
        object.__setattr__(self, "phi2_left", s("phi2_left"))
        object.__setattr__(self, "phi2_right", s("phi2_right"))
        object.__setattr__(self, "g_left", s("Lgamma_left"))
        object.__setattr__(self, "g_right", s("Lgamma_right"))
        object.__setattr__(self, "b_left", s("Lbeta_left"))
        object.__setattr__(self, "b_right", s("Lbeta_right"))
        object.__setattr__(self, "phi0_lag", s("phi0lag"))
        object.__setattr__(self, "phi1_lag", s("phi1lag"))
        object.__setattr__(self, "s0", s("s0"))
        object.__setattr__(self, "s1", s("s1"))


def _explicit_moment_equations(
    cell_h: sp.Expr,
    sig_t: sp.Expr,
    sig_s0: sp.Expr,
    sig_s1: sp.Expr,
    a_wd: sp.Expr,
    c: CellSymbols,
) -> dict[str, sp.Expr]:
    r"""The four moment residuals (16a–d) in explicit/opaque form.

    Derived once by :func:`derive_moment_equations` (which PROVES this
    form equals the raw angular reductions); consumed by the correction
    step. Larsen's forms are recovered at :math:`W_0 = 1, W_2 = 1/3`.
    """
    rho = rho_of(a_wd)
    return {
        # (16a) — exact in ℓ ≤ 1:
        "P0_balance": (
            (c.right[1] - c.left[1]) / cell_h
            + sig_t * c.bar[0]
            - sig_s0 * W0 * c.phi0_lag
            - c.s0
        ),
        # (16b) — the Legendre-recursion 2/3 and 1/3; W2 on the P1 source:
        "P1_balance": (
            (
                sp.Rational(2, 3) * (c.phi2_right - c.phi2_left)
                + sp.Rational(1, 3) * (c.right[0] - c.left[0])
            )
            / cell_h
            + sig_t * c.bar[1]
            - 3 * sig_s1 * W2 * c.phi1_lag
            - c.s1
        ),
        # (16c) — the ρ/(2W2) closure coupling + the opaque γ functional:
        "P0_closure": (
            c.bar[0]
            - (c.right[0] + c.left[0]) / 2
            - rho / (2 * W2) * (c.right[1] - c.left[1])
            - (c.g_right - c.g_left) / 2
        ),
        # (16d) — the ρ/(2W0) coupling + the opaque β functional:
        "P1_closure": (
            c.bar[1]
            - (c.right[1] + c.left[1]) / 2
            - rho / (2 * W0) * (c.right[0] - c.left[0])
            - (c.b_right - c.b_left) / 2
        ),
    }


def derive_moment_equations() -> None:
    r"""Steps 1–2: derive (16a–d) and PROVE the explicit/opaque form exact.

    Builds the per-ordinate balance (10a) and WD closure (10b) on raw
    edge/average flux symbols, takes the :math:`L_0`/:math:`L_1`
    reductions computationally, and verifies they equal
    :func:`_explicit_moment_equations` with every symbol REALIZED from
    the same raw fluxes (:math:`\phi_{n,loc} = L_n\psi_{loc}`,
    :math:`g_{loc} = L_0[\gamma\psi_{loc}]`,
    :math:`b_{loc} = L_0[\beta\psi_{loc}]`) — for the SYMBOLIC
    quadrature, i.e. for every symmetric rule. Nothing in (16) is
    transcribed; the split into explicit moments + opaque lagged
    functionals is exact by construction.
    """
    psi_l = sp.symbols(f"psiL0:{N_ORD}")
    psi_r = sp.symbols(f"psiR0:{N_ORD}")
    psi_b = sp.symbols(f"psiB0:{N_ORD}")
    c = CellSymbols("drv")
    cell = (H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I)
    a = A_WD_I
    rho = rho_of(a)

    def l_n(vec, n: int) -> sp.Expr:
        return sum(
            sp.legendre(n, ORDINATES[m][0]) * ORDINATES[m][1] * vec[m]
            for m in range(N_ORD)
        )

    def l_gamma(vec) -> sp.Expr:
        return sum(
            (wd_alpha(a, ORDINATES[m][0]) - ORDINATES[m][0] * rho / W2)
            * ORDINATES[m][1]
            * vec[m]
            for m in range(N_ORD)
        )

    def l_beta(vec) -> sp.Expr:
        return sum(
            (ORDINATES[m][0] * wd_alpha(a, ORDINATES[m][0]) - rho / W0)
            * ORDINATES[m][1]
            * vec[m]
            for m in range(N_ORD)
        )

    realization = {
        c.left[0]: l_n(psi_l, 0), c.left[1]: l_n(psi_l, 1),
        c.right[0]: l_n(psi_r, 0), c.right[1]: l_n(psi_r, 1),
        c.bar[0]: l_n(psi_b, 0), c.bar[1]: l_n(psi_b, 1),
        c.phi2_left: l_n(psi_l, 2), c.phi2_right: l_n(psi_r, 2),
        c.g_left: l_gamma(psi_l), c.g_right: l_gamma(psi_r),
        c.b_left: l_beta(psi_l), c.b_right: l_beta(psi_r),
    }

    def balance(m: int) -> sp.Expr:
        mu, _ = ORDINATES[m]
        return (
            mu / H_I * (psi_r[m] - psi_l[m])
            + SIG_T_I * psi_b[m]
            - SIG_S0_I * c.phi0_lag
            - 3 * SIG_S1_I * mu * c.phi1_lag
        )

    def closure(m: int) -> sp.Expr:
        mu, _ = ORDINATES[m]
        alpha = wd_alpha(A_WD_I, mu)
        return (
            psi_b[m]
            - (1 + alpha) / 2 * psi_r[m]
            - (1 - alpha) / 2 * psi_l[m]
        )

    raw = {
        "P0_balance": l_n([balance(m) for m in range(N_ORD)], 0) - c.s0,
        "P1_balance": l_n([balance(m) for m in range(N_ORD)], 1) - c.s1,
        "P0_closure": l_n([closure(m) for m in range(N_ORD)], 0),
        "P1_closure": l_n([closure(m) for m in range(N_ORD)], 1),
    }
    explicit = _explicit_moment_equations(*cell, c)
    for name in raw:
        realized = explicit[name].subs(realization, simultaneous=True)
        _require(
            _is_zero(raw[name] - realized),
            f"(16) [{name}]: the explicit/opaque form must equal the raw "
            f"angular reduction for the symbolic quadrature",
        )


def _self_check_steps12() -> None:
    derive_closure_weight_identities()
    prove_p_recursion_lemma()
    derive_moment_equations()


# ═══════════════════════════════════════════════════════════════════════
# Steps 3–4 — promotion, subtraction, elimination
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class CorrectionSymbols:
    r"""The correction unknowns of one cell.

    :math:`(f_0, f_1)` at ``left`` / ``right`` edges and the ``bar``
    average, plus the sweep displacement moments ``d0`` / ``d1``
    (:math:`d_n = \phi_{ni}^{l+1/2} - \phi_{ni}^{l}` — the raw material
    of the DSA residual source; Larsen's :math:`g` sources are
    cross-section-weighted multiples of these).
    """

    tag: str

    def __post_init__(self) -> None:
        s = lambda name: sp.Symbol(f"{name}_{self.tag}")  # noqa: E731
        for loc in ("left", "right", "bar"):
            object.__setattr__(self, loc, (s(f"f0_{loc}"), s(f"f1_{loc}")))
        object.__setattr__(self, "d0", s("d0"))
        object.__setattr__(self, "d1", s("d1"))


#: The canonical per-cell correction-symbol instances. Shared across
#: every consumer so the cached per-cell solves are computed ONCE.
COR_I = CorrectionSymbols("celli")
COR_IP1 = CorrectionSymbols("cellip1")


def derive_correction_system(
    cell: tuple[sp.Expr, ...], mom: CellSymbols, cor: CorrectionSymbols
) -> dict[str, sp.Expr]:
    r"""Steps 3–4 for one cell: promote, subtract, prove the lag terms drop.

    Step 3 (Larsen (18), the "discrete P1 approximation" slot
    assignment): promote every explicit :math:`\ell \le 1` moment
    (:math:`\phi_n \to \phi_n + f_n` at edges and average) and the
    scattering source (:math:`\phi_{ni}^l \to \phi_{ni}^{l+1/2} +
    f_{ni}`); the :math:`\phi_2` terms and the opaque closure
    functionals :math:`L_0[\gamma\psi]`, :math:`L_0[\beta\psi]` stay at
    :math:`l+1/2` — they are lagged, structurally untouchable here
    because they are opaque symbols.

    Step 4: promoted minus original, with :math:`\phi_{ni}^l =
    \phi_{ni}^{l+1/2} - d_n`. The derivation-grade check: the result
    must be free of EVERY absolute symbol — moments, lagged
    functionals, sources — leaving a linear system in
    :math:`\{f, d\}` alone (Larsen (20)–(21)).
    """
    original = _explicit_moment_equations(*cell, mom)

    promo: dict[sp.Symbol, sp.Expr] = {}
    for loc in ("left", "right", "bar"):
        phi, f = getattr(mom, loc), getattr(cor, loc)
        for n in range(2):
            promo[phi[n]] = phi[n] + f[n]
    promo[mom.phi0_lag] = mom.bar[0] + cor.bar[0]
    promo[mom.phi1_lag] = mom.bar[1] + cor.bar[1]

    lag_realization = {
        mom.phi0_lag: mom.bar[0] - cor.d0,
        mom.phi1_lag: mom.bar[1] - cor.d1,
    }
    absolute_symbols = (
        [getattr(mom, loc)[n] for loc in ("left", "right", "bar") for n in (0, 1)]
        + [mom.phi2_left, mom.phi2_right, mom.g_left, mom.g_right,
           mom.b_left, mom.b_right, mom.s0, mom.s1]
    )

    corrections: dict[str, sp.Expr] = {}
    for name, residual in original.items():
        promoted = residual.subs(promo, simultaneous=True)
        delta = sp.expand(promoted - residual).subs(lag_realization)
        delta = sp.expand(delta)
        for sym in absolute_symbols:
            _require(
                _is_zero(sp.diff(delta, sym)),
                f"step 4 [{name}]: every absolute symbol must cancel in the "
                f"subtraction — {sym} survived",
            )
        corrections[name] = delta
    return corrections


@lru_cache(maxsize=None)
def eliminate_to_edge_pair(
    cell: tuple[sp.Expr, ...], cor: CorrectionSymbols
) -> tuple[tuple[sp.Expr, sp.Expr], dict[sp.Symbol, sp.Expr]]:
    r"""Eliminate the cell-average corrections: Larsen (22)/(24).

    Substitutes the closure-derived cell averages (21a/21b) into the
    balance corrections (20a/20b), producing the cell's two edge-only
    equations (residual form, == 0) — and returns the cell-average
    solutions alongside (the raw material of the updates (28)).
    """
    mom = CellSymbols(f"tmp_{cor.tag}")
    system = derive_correction_system(cell, mom, cor)
    bar_solution = sp.solve(
        [system["P0_closure"], system["P1_closure"]],
        [cor.bar[0], cor.bar[1]],
        dict=True,
    )
    _require(
        len(bar_solution) == 1,
        "(21): the closure corrections must determine the cell averages "
        "uniquely",
    )
    sub = bar_solution[0]
    return (
        (
            sp.expand(system["P0_balance"].subs(sub)),
            sp.expand(system["P1_balance"].subs(sub)),
        ),
        sub,
    )


@lru_cache(maxsize=None)
def solve_cell_for_edge_f1(
    cell: tuple[sp.Expr, ...], cor: CorrectionSymbols
) -> dict[sp.Symbol, sp.Expr]:
    r"""Solve one cell's edge pair for its edge :math:`f_1` values.

    Larsen (25)/(26): the cell's two edge-only equations, solved as a
    2×2 linear system for :math:`(f_{1,\mathrm{left}},
    f_{1,\mathrm{right}})` in terms of the cell's :math:`f_0` edge
    values and its displacement sources. Instantiated at a shared edge
    from the two adjacent cells, the two expressions' equality IS the
    consistency condition that becomes the tridiagonal (27).
    """
    (p0_eq, p1_eq), _bars = eliminate_to_edge_pair(cell, cor)
    solution = sp.solve(
        [p0_eq, p1_eq], [cor.left[1], cor.right[1]], dict=True
    )
    _require(
        len(solution) == 1,
        "(25)/(26): the edge pair must determine the edge f1 values "
        "uniquely",
    )
    return solution[0]


# ═══════════════════════════════════════════════════════════════════════
# The main theorem — the interior tridiagonal row ≡ Larsen (27)/(23a–f)
# ═══════════════════════════════════════════════════════════════════════

#: The three interior edge corrections around the shared edge i+1/2 and
#: the two cells' displacement sources — the unknowns of the row.
F0_LEFT_EDGE, F0_SHARED, F0_RIGHT_EDGE = sp.symbols(
    "f0_im12 f0_ip12 f0_ip32"
)
D0_I, D1_I, D0_IP1, D1_IP1 = sp.symbols("d0_i d1_i d0_ip1 d1_ip1")


def larsen_23_coefficients(
    h: sp.Expr,
    sig_t: sp.Expr,
    sig_s0: sp.Expr,
    sig_s1: sp.Expr,
    a_wd: sp.Expr,
) -> dict[str, sp.Expr]:
    r"""The TRANSCRIBED coefficient set (23a–f) — the comparison target.

    This is the module's only transcription from [Larsen1982a]_ (printed
    p. 53, cross-checked against the literature memo §2.2). The
    derivation must REPRODUCE the row built from these; nothing else
    consumes them.

    .. math::

        \hat\sigma_{R} &= \frac{\sigma_T - \sigma_{S0}}
            {1 + \tfrac{3}{2}\rho(\sigma_T-\sigma_{S0})h}, \quad
        \hat\sigma_{S} = \frac{\sigma_{S0}}
            {1 + \tfrac{3}{2}\rho(\sigma_T-\sigma_{S0})h}, \\
        D &= \frac{1}{3(\sigma_T-\sigma_{S1})} + \frac{\rho h}{2}, \quad
        a = \frac{\sigma_{S1}}{\sigma_T-\sigma_{S1}}, \\
        g_0 &= \hat\sigma_{S}\,h\,d_0, \qquad g_1 = a\,d_1 .
    """
    rho = rho_of(a_wd)
    sig_r = sig_t - sig_s0
    denom = 1 + sp.Rational(3, 2) * rho * sig_r * h
    return {
        "sigma_R_hat": sig_r / denom,
        "sigma_S_hat": sig_s0 / denom,
        "D": 1 / (3 * (sig_t - sig_s1)) + rho * h / 2,
        "a": sig_s1 / (sig_t - sig_s1),
    }


def larsen_27_row() -> sp.Expr:
    r"""The TRANSCRIBED interior row (27) (residual form, == 0).

    .. math::

        -\frac{D_{i+1}}{h_{i+1}}(f_{0,i+3/2} - f_{0,i+1/2})
        + \frac{D_i}{h_i}(f_{0,i+1/2} - f_{0,i-1/2})
        + \tfrac14\bigl[\hat\sigma_{R,i+1}h_{i+1}(f_{0,i+3/2}+f_{0,i+1/2})
                       + \hat\sigma_{R,i}h_i(f_{0,i+1/2}+f_{0,i-1/2})\bigr]
        - \tfrac12(g_{0,i+1} + g_{0,i}) + (g_{1,i+1} - g_{1,i})
    """
    ci = larsen_23_coefficients(H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I)
    cp = larsen_23_coefficients(
        H_IP1, SIG_T_IP1, SIG_S0_IP1, SIG_S1_IP1, A_WD_IP1
    )
    g0_i, g0_p = ci["sigma_S_hat"] * H_I * D0_I, cp["sigma_S_hat"] * H_IP1 * D0_IP1
    g1_i, g1_p = ci["a"] * D1_I, cp["a"] * D1_IP1
    return (
        -cp["D"] / H_IP1 * (F0_RIGHT_EDGE - F0_SHARED)
        + ci["D"] / H_I * (F0_SHARED - F0_LEFT_EDGE)
        + sp.Rational(1, 4)
        * (
            cp["sigma_R_hat"] * H_IP1 * (F0_RIGHT_EDGE + F0_SHARED)
            + ci["sigma_R_hat"] * H_I * (F0_SHARED + F0_LEFT_EDGE)
        )
        - sp.Rational(1, 2) * (g0_p + g0_i)
        + (g1_p - g1_i)
    )


def derive_interior_row() -> sp.Expr:
    r"""THE MAIN THEOREM: shared-edge :math:`f_1` continuity ≡ Larsen (27).

    Solves cell :math:`i` and cell :math:`i+1` for their edge
    :math:`f_1` values, instantiates both at the SHARED edge
    :math:`i+1/2` on the global edge unknowns, and forms the continuity
    residual :math:`f_{1,i+1/2}^{(i)} - f_{1,i+1/2}^{(i+1)}`. The
    theorem: under Larsen's printed convention (:math:`W_0 = 1`,
    :math:`W_2 = 1/3`) this residual is a nonzero multiple of the
    transcribed row (27) with the coefficient set (23a–f) — i.e. the
    four-step, executed symbolically, PRODUCES the printed consistent
    low-order operator.

    Returns the derived row (the general-quadrature form, before the
    printed-convention substitution) for downstream reuse.
    """
    cor_i, cor_p = COR_I, COR_IP1

    f1_i = solve_cell_for_edge_f1(
        (H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I), cor_i
    )
    f1_p = solve_cell_for_edge_f1(
        (H_IP1, SIG_T_IP1, SIG_S0_IP1, SIG_S1_IP1, A_WD_IP1), cor_p
    )

    into_global_i = {
        cor_i.left[0]: F0_LEFT_EDGE,
        cor_i.right[0]: F0_SHARED,
        cor_i.d0: D0_I,
        cor_i.d1: D1_I,
    }
    into_global_p = {
        cor_p.left[0]: F0_SHARED,
        cor_p.right[0]: F0_RIGHT_EDGE,
        cor_p.d0: D0_IP1,
        cor_p.d1: D1_IP1,
    }
    f1_shared_from_i = f1_i[cor_i.right[1]].subs(into_global_i)
    f1_shared_from_p = f1_p[cor_p.left[1]].subs(into_global_p)
    derived_row = f1_shared_from_i - f1_shared_from_p

    # ── the printed-convention comparison ────────────────────────────
    # The continuity residual determines the row up to a nonzero scalar
    # FIELD (it is one equation, not a normalized operator row), so the
    # theorem is exact proportionality: every 2×2 minor of the two rows'
    # coefficient vectors over the row unknowns must vanish. Checking
    # minors avoids ever forming (or simplifying) the scale quotient.
    unknowns = (
        F0_LEFT_EDGE, F0_SHARED, F0_RIGHT_EDGE, D0_I, D1_I, D0_IP1, D1_IP1
    )

    def coeffs(row: sp.Expr) -> list[sp.Expr]:
        collected = {u: sp.diff(row, u) for u in unknowns}
        constant = row - sum(c * u for u, c in collected.items())
        _require(
            _rational_zero(constant),
            "(27): the row must be homogeneous-linear in the unknowns",
        )
        return [collected[u] for u in unknowns]

    derived_c = [
        _to_printed_convention(c) for c in coeffs(sp.expand(derived_row))
    ]
    target_c = [
        _to_printed_convention(c) for c in coeffs(sp.expand(larsen_27_row()))
    ]
    pivot = unknowns.index(F0_SHARED)
    _require(
        not _rational_zero(derived_c[pivot]),
        "(27): the derived row must actually involve the shared edge",
    )
    for k, u in enumerate(unknowns):
        minor = (
            target_c[k] * derived_c[pivot] - target_c[pivot] * derived_c[k]
        )
        _require(
            _rational_zero(minor),
            f"(27): coefficient proportionality must hold at {u} — the "
            f"derived shared-edge continuity row must be a scalar multiple "
            f"of the transcribed Larsen row with coefficients (23a-f)",
        )
    return derived_row


def derive_update_relations() -> None:
    r"""The cell-average updates ≡ Larsen (28a/28b).

    Substituting the cell's solved edge :math:`f_1` values into its
    closure-derived averages (21a/21b) must reproduce the printed
    updates

    .. math::

        f_{0i} &= \bigl(\tfrac12 - \tfrac34\rho_i\hat\sigma_{Ri}h_i\bigr)
            (f_{0,i+1/2} + f_{0,i-1/2}) + \tfrac32\rho_i g_{0i}, \\
        f_{1i} &= \bigl(\tfrac12\rho_i - D_i/h_i\bigr)
            (f_{0,i+1/2} - f_{0,i-1/2}) + g_{1i},

    under the printed convention. (:math:`\phi_n^{l+1} = \phi_n^{l+1/2}
    + f_{ni}` is then the accelerated iterate.)
    """
    cor = COR_I
    cell = (H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I)
    _pair, bars = eliminate_to_edge_pair(cell, cor)
    f1 = solve_cell_for_edge_f1(cell, cor)

    rho = rho_of(A_WD_I)
    co = larsen_23_coefficients(H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I)
    g0 = co["sigma_S_hat"] * H_I * cor.d0
    g1 = co["a"] * cor.d1
    edge_sum = cor.right[0] + cor.left[0]
    edge_diff = cor.right[0] - cor.left[0]

    targets = {
        cor.bar[0]: (
            (sp.Rational(1, 2)
             - sp.Rational(3, 4) * rho * co["sigma_R_hat"] * H_I) * edge_sum
            + sp.Rational(3, 2) * rho * g0
        ),
        cor.bar[1]: (rho / 2 - co["D"] / H_I) * edge_diff + g1,
    }
    for which, (bar_sym, target) in zip(("28a", "28b"), targets.items()):
        derived = bars[bar_sym].subs(f1)
        _require(
            _rational_zero(
                _to_printed_convention(derived)
                - _to_printed_convention(target)
            ),
            f"({which}): the update relation must match the printed form",
        )


def derive_boundary_rows() -> dict[str, sp.Expr]:
    r"""The correction-equation boundary rows: Marshak (38a) + reflecting (39a).

    The correction edge angular flux is the :math:`\ell \le 1` synthesis
    :math:`\Psi_m = f_0 + (\mu_m/W_2) f_1` (Larsen (33) with the honest
    :math:`1/W_2`; his "3"). At a left vacuum boundary the prescribed
    zero-incident condition is imposed on :math:`\Psi` in the Marshak
    (half-range current) sense:

    .. math::

        0 = \sum_{\mu_m > 0} \mu_m \omega_m \Psi_m
          = \gamma_N f_{0,1/2} + \frac{W_2^+}{W_2} f_{1,1/2},

    whose coefficients are DERIVED here from the half-range reduction
    and checked against the printed (38a) pair :math:`(\gamma_N, 1/2)`
    under the printed convention — the :math:`1/2` is
    :math:`3\sum_{\mu>0}\mu_m^2\omega_m`, exact for symmetric rules with
    :math:`W_2 = 1/3`, NOT a continuum constant. Reflecting (39a) is
    :math:`f_{1,1/2} = 0`.

    Returns the two rows (in the boundary-edge :math:`f_0, f_1`
    symbols) plus the CLOSED vacuum row: :math:`f_{1,1/2}` eliminated
    via the boundary cell's one-sided relation (the (26) analog), i.e.
    the actual first row of the consistent low-order system in
    :math:`f_0` alone.
    """
    f0_b, f1_b = sp.symbols("f0_bnd f1_bnd")

    psi_of = lambda mu: f0_b + mu / W2 * f1_b  # noqa: E731 — (33)
    marshak = half_range_sum(lambda mu: mu * psi_of(mu), positive=True)
    marshak = sp.expand(marshak)

    _require(
        _rational_zero(marshak.coeff(f0_b) - GAMMA_HALF),
        "(38a): the f0 coefficient must be the half-range current γ_N",
    )
    _require(
        _rational_zero(
            _to_printed_convention(marshak.coeff(f1_b)) - sp.Rational(1, 2)
        ),
        "(38a): the f1 coefficient must be W2⁺/W2 — printed value 1/2",
    )

    # Close the row: f1 at the boundary edge from the boundary cell's
    # solved pair (the one-sided (26) analog; cell i plays the boundary
    # cell with its LEFT edge on the boundary).
    cor = COR_I
    cell = (H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I)
    f1 = solve_cell_for_edge_f1(cell, cor)
    one_sided = {f0_b: cor.left[0], f1_b: f1[cor.left[1]]}
    closed_vacuum_row = sp.together(marshak.subs(one_sided))

    _require(
        not _rational_zero(sp.diff(closed_vacuum_row, cor.right[0])),
        "the closed vacuum row must couple to the interior edge",
    )
    return {
        "marshak_open": marshak,
        "vacuum_closed": closed_vacuum_row,
        "reflecting": f1[cor.left[1]],  # == 0 is the reflecting row
    }


def derive_dd_instance() -> dict[str, sp.Expr]:
    r"""The diamond-difference member: coefficients (23a–f) at α = 0.

    The DD low-order operator — the object Phase 3b wires — with the
    quadrature moments still symbolic where they enter:

    .. math::

        \hat\sigma_R = \sigma_T - \sigma_{S0}, \quad
        \hat\sigma_S = \sigma_{S0}, \quad
        D = \frac{1}{3(\sigma_T - \sigma_{S1})}, \quad
        a = \frac{\sigma_{S1}}{\sigma_T - \sigma_{S1}},

    i.e. Alcouffe's edge-centered scheme with the ¼-weighted three-point
    removal — recovered here as the :math:`a_{\rm wd} \to 0` member of
    the proven WD family (nothing DD-specific is re-derived, and none of
    Alcouffe's errata-bearing printed forms are consumed).
    """
    dd = {A_WD_I: sp.Integer(0), A_WD_IP1: sp.Integer(0)}
    co = {
        k: sp.simplify(v.subs(dd))
        for k, v in larsen_23_coefficients(
            H_I, SIG_T_I, SIG_S0_I, SIG_S1_I, A_WD_I
        ).items()
    }
    _require(
        _is_zero(co["sigma_R_hat"] - (SIG_T_I - SIG_S0_I)),
        "DD: σ̂_R must reduce to the plain removal σ_T − σ_S0",
    )
    _require(
        _is_zero(co["D"] - 1 / (3 * (SIG_T_I - SIG_S1_I))),
        "DD: D must reduce to 1/[3(σ_T − σ_S1)]",
    )
    return co


def _self_check_steps34() -> None:
    derive_interior_row()
    derive_update_relations()
    derive_boundary_rows()
    derive_dd_instance()


def run_all() -> None:
    """Run the complete derivation battery (the test entry point)."""
    _self_check_steps12()
    _self_check_steps34()


if __name__ == "__main__":  # pragma: no cover — manual derivation run
    _self_check_steps12()
    print("steps 1-2 OK: annihilation, P-recursion, (16a-d) split exact")
    derive_interior_row()
    print(
        "steps 3-4 OK: THE MAIN THEOREM — the four-step's shared-edge f1 "
        "continuity reproduces Larsen (27) with (23a-f) at W0=1, W2=1/3"
    )
    derive_update_relations()
    print("updates OK: (28a/28b) reproduced")
    derive_boundary_rows()
    print("boundary OK: Marshak (38a) coefficients derived; row closed")
    derive_dd_instance()
    print("DD instance OK: (23a-f) at α = 0")
