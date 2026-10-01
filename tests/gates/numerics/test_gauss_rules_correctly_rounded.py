r"""Gate (a) of the platform-independence remedy: every Gauss rule is CORRECTLY ROUNDED.

**The claim.** ``GeneratingMeasure.gauss(n)`` returns, for every shipped
family and every ``n``, the nodes and weights of the exact ``n``-point Gauss
rule of that measure, each real number rounded ONCE to the nearest float64.
A correctly rounded value is unique, so this rule is the same bytes on every
platform, every LAPACK and every thread count; it is the property whose
absence let a macOS update (Accelerate, 2026-09-30) move the Golub–Welsch
rules by 1–445 ULP and redden 24 downstream pins (``[M]``
``scratch/platform_drift/memo.md``).

**Claim kind: REFERENCE.** The reference is computed here, at
``_DPS = 80`` decimal digits, by a route that shares nothing with production
above the trusted-library line (``instrument-doctrine`` X4):

* production evaluates each family's MONIC recurrence coefficients
  :math:`(\alpha_k, \beta_k)` (Gautschi's convention) and diagonalises the
  Jacobi matrix (Golub–Welsch);
* the reference finds each node by NEWTON's method on the classical
  orthogonal polynomial in its standard (non-monic) normalisation — Bonnet's
  recurrence for Legendre (DLMF 18.9.1), the three-term Jacobi recurrence of
  DLMF 18.9.2, mpmath's own ``laguerre`` / ``hermite`` (hypergeometric
  definitions) — and takes each weight from the CHRISTOFFEL closed form in
  the polynomial's derivative (Szegő, *Orthogonal Polynomials*, 1975,
  Eqs. 15.3.1, 15.3.5 and 15.3.6; Abramowitz & Stegun 25.4.29, .38, .45, .46);
  the two Chebyshev families use their trigonometric closed forms
  (A&S 25.4.38–39), written through ``sin`` so the centre node of an odd rule
  is an exact zero.

The Newton seeds come from ``scipy.special.roots_*`` (an independent float
implementation); a seed only selects WHICH root Newton converges to, and the
rows assert ``n`` distinct roots, so a bad seed reddens rather than passes.
The reference itself is verified twice, below production (the
``foundation`` rows): (i) the hand-written recurrences agree with mpmath's
hypergeometric definitions at a generic point; (ii) every reference rule
integrates the monomials :math:`x^0..x^{2n-1}` against its own weight to 70
digits, with the moments computed in closed form (Gamma/Beta functions), the
independent ground the Christoffel formulas are checked against.

**Rounding.** Each 80-digit value is rounded with mpmath's explicit
round-to-nearest (``to_float(..., rnd=round_nearest)``; ``float(mpf)`` follows
the context's rounding mode, which is the wrong thing to rely on). A value
within :math:`10^{-70}` relative of a rounding midpoint RAISES
(``_NearTie``), because there the 80-digit value cannot certify which side the
exact root is on. ``[M]`` 2026-10-01: no value in the population comes close.

**First red at HEAD** (``0a5a23fa``, Golub–Welsch ``eigh``): ``[M]``
``.venv/bin/python -O -m pytest tests/gates/numerics/test_gauss_rules_correctly_rounded.py``
— the equality rows red for every family at ``n >= 3`` (LEGENDRE n=4: nodes 1,
weights 7 ULP; n=64: 5 / 445 ULP; HERMITE n=64: weights off by 6.7e16 ULP,
i.e. tiny tail weights with no relative accuracy) and for ``jacobi(1.5, 2.25)``
already at ``n = 1`` (its mass ``exp(gammaln(...))`` is 2 ULP from
:math:`2^{a+b+1}B(a+1,b+1)`). The symmetry rows are green at HEAD (production
IMPOSES the mirror) and are the law that must SURVIVE the retirement of that
imposition; their tooth is the un-imposed Golub–Welsch rule (gates_spec.md §a).

Runtime ``[M]`` about 6 s serial; the cost is the reference, not production.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import mpmath as mp
import numpy as np
import pytest
from mpmath.libmp import round_nearest, to_float
from scipy import special

from orpheus.numerics.generating_measure import (
    CHEBYSHEV_T,
    CHEBYSHEV_U,
    HERMITE,
    LEGENDRE,
    GeneratingMeasure,
    jacobi,
    laguerre,
)

_HERE = "tests/gates/numerics/test_gauss_rules_correctly_rounded.py"

#: An mpmath real (``mpmath.mpf``). mpmath ships no type stubs and its ``mpf``
#: is a per-context class, so the annotation names the role, not the class.
_Real = Any

#: Working precision of the reference, in decimal digits. Production is
#: specified at ~40 and ~60 digits; the reference sits ABOVE both so that a
#: production rounding error cannot be shared with it.
_DPS = 80

#: Newton stops when the step is below this relative size; one more step is
#: then taken, so the converged root is good to ~2x this many digits minus
#: the polynomial's conditioning. 10**-(_DPS-5) keeps every root at >= 70
#: correct digits for n <= 128 ([M]: the moment rows close to 1e-70).
_NEWTON_TOL_DIGITS = _DPS - 5

#: A rounding decision is trusted only if the 80-digit value is farther than
#: this (relative) from a float64 midpoint.
_TIE_MARGIN_DIGITS = 70


class _NearTie(ArithmeticError):
    """The 80-digit value sits too close to a float64 rounding midpoint."""


# ── the reference: Newton on the classical polynomial + Christoffel weights ──


def _newton(
    value_and_slope: Callable[[_Real], tuple[_Real, _Real]], seeds: np.ndarray
) -> list[_Real]:
    tolerance = mp.mpf(10) ** (-_NEWTON_TOL_DIGITS)
    roots = []
    for seed in seeds:
        x = mp.mpf(float(seed))
        for _ in range(200):
            value, slope = value_and_slope(x)
            step = value / slope
            x -= step
            if abs(step) <= tolerance * max(1, abs(x)):
                value, slope = value_and_slope(x)
                x -= value / slope
                break
        else:  # pragma: no cover - a reference failure, not a production one
            raise ArithmeticError(f"Newton did not converge from seed {seed!r}")
        roots.append(x)
    return roots


def _legendre_bonnet(n: int, x: _Real) -> tuple[_Real, _Real]:
    """P_n(x) and P_n'(x) by Bonnet's recurrence (DLMF 18.9.1, a=b=0)."""
    p_prev, p = mp.mpf(1), x
    for k in range(1, n):
        p_prev, p = p, ((2 * k + 1) * x * p - k * p_prev) / (k + 1)
    return p, n * (x * p - p_prev) / (x * x - 1)


def _jacobi_dlmf(n: int, a: _Real, b: _Real, x: _Real) -> _Real:
    """P_n^{(a,b)}(x) by the three-term recurrence of DLMF 18.9.2 (non-monic)."""
    if n == 0:
        return mp.mpf(1)
    p_prev, p = mp.mpf(1), ((a + b + 2) * x + (a - b)) / 2
    for k in range(1, n):
        s = 2 * k + a + b
        lead = 2 * (k + 1) * (k + a + b + 1) * s
        p_prev, p = p, (
            (s + 1) * ((s + 2) * s * x + a * a - b * b) * p
            - 2 * (k + a) * (k + b) * (s + 2) * p_prev
        ) / lead
    return p


def _rule_legendre(n: int) -> tuple[list[_Real], list[_Real]]:
    if n == 1:
        return [mp.mpf(0)], [mp.mpf(2)]
    nodes = _newton(lambda x: _legendre_bonnet(n, x), special.roots_legendre(n)[0])
    # Christoffel numbers, Szego (15.3.1) / A&S 25.4.29: w = 2 / ((1-x^2) P_n'(x)^2)
    return nodes, [2 / ((1 - x * x) * _legendre_bonnet(n, x)[1] ** 2) for x in nodes]


def _rule_jacobi(a_: float, b_: float) -> Callable[[int], tuple[list, list]]:
    def rule(n: int) -> tuple[list[_Real], list[_Real]]:
        a, b = mp.mpf(a_), mp.mpf(b_)
        if n == 1:
            return [(b - a) / (a + b + 2)], [_mass_jacobi(a, b)]

        def value_and_slope(x: _Real) -> tuple[_Real, _Real]:
            # d/dx P_n^{(a,b)} = (n+a+b+1)/2 P_{n-1}^{(a+1,b+1)}  (DLMF 18.9.15)
            return _jacobi_dlmf(n, a, b, x), (n + a + b + 1) / 2 * _jacobi_dlmf(
                n - 1, a + 1, b + 1, x
            )

        nodes = _newton(value_and_slope, special.roots_jacobi(n, a_, b_)[0])
        # Szego (15.3.5): w = 2^{a+b+1} G(n+a+1) G(n+b+1) / (G(n+a+b+1) n! (1-x^2) P_n'(x)^2)
        c = (
            2 ** (a + b + 1)
            * mp.gamma(n + a + 1)
            * mp.gamma(n + b + 1)
            / (mp.gamma(n + a + b + 1) * mp.factorial(n))
        )
        return nodes, [c / ((1 - x * x) * value_and_slope(x)[1] ** 2) for x in nodes]

    return rule


def _rule_laguerre(a_: float) -> Callable[[int], tuple[list, list]]:
    def rule(n: int) -> tuple[list[_Real], list[_Real]]:
        a = mp.mpf(a_)
        if n == 1:
            return [a + 1], [mp.gamma(a + 1)]

        def value_and_slope(x: _Real) -> tuple[_Real, _Real]:
            # d/dx L_n^{(a)} = -L_{n-1}^{(a+1)}  (DLMF 18.9.23)
            return mp.laguerre(n, a, x), -mp.laguerre(n - 1, a + 1, x)

        nodes = _newton(value_and_slope, special.roots_genlaguerre(n, a_)[0])
        # Szego (15.3.6) / A&S 25.4.45: w = G(n+a+1) / (n! x L_n^{(a)}'(x)^2)
        c = mp.gamma(n + a + 1) / mp.factorial(n)
        return nodes, [c / (x * value_and_slope(x)[1] ** 2) for x in nodes]

    return rule


def _rule_hermite(n: int) -> tuple[list[_Real], list[_Real]]:
    if n == 1:
        return [mp.mpf(0)], [mp.sqrt(mp.pi)]
    # H_n' = 2n H_{n-1}  (DLMF 18.9.25)
    nodes = _newton(
        lambda x: (mp.hermite(n, x), 2 * n * mp.hermite(n - 1, x)),
        special.roots_hermite(n)[0],
    )
    # A&S 25.4.46: w = 2^{n-1} n! sqrt(pi) / (n^2 H_{n-1}(x)^2)
    c = 2 ** (n - 1) * mp.factorial(n) * mp.sqrt(mp.pi) / n**2
    return nodes, [c / mp.hermite(n - 1, x) ** 2 for x in nodes]


def _rule_chebyshev_t(n: int) -> tuple[list[_Real], list[_Real]]:
    # x_i = cos((2i-1)pi/(2n)) = sin(pi (n+1-2i) / (2n)); w_i = pi/n  (A&S 25.4.38)
    nodes = [mp.sin(mp.pi * (n + 1 - 2 * i) / (2 * n)) for i in range(1, n + 1)]
    return nodes, [mp.pi / n] * n


def _rule_chebyshev_u(n: int) -> tuple[list[_Real], list[_Real]]:
    # x_i = cos(i pi/(n+1)); w_i = pi/(n+1) sin^2(i pi/(n+1))  (A&S 25.4.40)
    angles = [mp.pi * (n + 1 - 2 * i) / (2 * (n + 1)) for i in range(1, n + 1)]
    return [mp.sin(t) for t in angles], [mp.pi / (n + 1) * mp.cos(t) ** 2 for t in angles]


def _mass_jacobi(a: _Real, b: _Real) -> _Real:
    return 2 ** (a + b + 1) * mp.beta(a + 1, b + 1)


def _round_once(value: _Real) -> float:
    """Round an 80-digit value to the nearest float64, refusing a near-tie."""
    nearest = to_float(value._mpf_, rnd=round_nearest)
    if nearest != 0.0:
        below = math.nextafter(nearest, -math.inf)
        above = math.nextafter(nearest, math.inf)
        margin = min(
            abs(value - (mp.mpf(nearest) + mp.mpf(below)) / 2),
            abs(value - (mp.mpf(nearest) + mp.mpf(above)) / 2),
        )
        if margin <= abs(value) * mp.mpf(10) ** (-_TIE_MARGIN_DIGITS):
            raise _NearTie(f"{mp.nstr(value, 30)} is within 1e-{_TIE_MARGIN_DIGITS} of a midpoint")
    return nearest


def reference_rule(
    rule: Callable[[int], tuple[list, list]], n: int
) -> tuple[np.ndarray, np.ndarray]:
    """The correctly rounded (nodes, weights) of ``rule(n)``, ascending nodes."""
    with mp.workdps(_DPS):
        nodes, weights = rule(n)
        order = sorted(range(n), key=lambda i: nodes[i])
        return (
            np.array([_round_once(nodes[i]) for i in order]),
            np.array([_round_once(weights[i]) for i in order]),
        )


# ── the population ──────────────────────────────────────────────────────────

#: (id, production measure, reference rule, exact zeroth moment, even weight?)
_FAMILIES: list[tuple[str, GeneratingMeasure, Callable, Callable[[], _Real], bool]] = [
    ("legendre", LEGENDRE, _rule_legendre, lambda: mp.mpf(2), True),
    ("chebyshev_t", CHEBYSHEV_T, _rule_chebyshev_t, lambda: +mp.pi, True),
    ("chebyshev_u", CHEBYSHEV_U, _rule_chebyshev_u, lambda: mp.pi / 2, True),
    ("hermite", HERMITE, _rule_hermite, lambda: mp.sqrt(mp.pi), True),
    ("jacobi(2,2)", jacobi(2.0, 2.0), _rule_jacobi(2.0, 2.0),
     lambda: _mass_jacobi(mp.mpf(2), mp.mpf(2)), True),
    ("jacobi(1.5,2.25)", jacobi(1.5, 2.25), _rule_jacobi(1.5, 2.25),
     lambda: _mass_jacobi(mp.mpf(1.5), mp.mpf(2.25)), False),
    ("jacobi(-0.5,0.75)", jacobi(-0.5, 0.75), _rule_jacobi(-0.5, 0.75),
     lambda: _mass_jacobi(mp.mpf(-0.5), mp.mpf(0.75)), False),
    ("laguerre(0)", laguerre(), _rule_laguerre(0.0), lambda: mp.mpf(1), False),
    ("laguerre(1.5)", laguerre(1.5), _rule_laguerre(1.5), lambda: mp.gamma(mp.mpf(2.5)), False),
]

#: n = 1 (the eigensolve-free branch), 2 (fixed by symmetry + mass at HEAD),
#: odd orders (an exact centre node for an even weight), the orders production
#: reaches (4, 8, 16, 64) and 33 (prime, odd). 128 is Legendre-only (the
#: plan's largest measured order; the other families' references cost 0.5 s
#: each there).
_ORDERS = (1, 2, 3, 4, 5, 8, 16, 33, 64)
_LEGENDRE_ORDERS = _ORDERS + (128,)


def _cases(even_only: bool = False) -> list:
    out = []
    for label, measure, rule, mass, even in _FAMILIES:
        if even_only and not even:
            continue
        orders = _LEGENDRE_ORDERS if label == "legendre" else _ORDERS
        out += [pytest.param(measure, rule, mass, n, id=f"{label}-n{n}") for n in orders]
    return out


def _ulps(a: np.ndarray, b: np.ndarray) -> int:
    return int(np.max(np.abs(a.view(np.int64) - b.view(np.int64))))


# ── foundation: the reference is right, below production ───────────────────

_REF_SELF_CHECK = f"{_HERE}::test_reference_recurrences_agree_with_mpmath_definitions"
_REF_MOMENTS = f"{_HERE}::test_reference_rules_integrate_their_moments_exactly"


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 5, 16, 64])
def test_reference_recurrences_agree_with_mpmath_definitions(n: int) -> None:
    """The two hand-typed recurrences of the reference (Bonnet, DLMF 18.9.1; Jacobi, 18.9.2)
    agree with mpmath's hypergeometric ``legendre`` / ``jacobi`` at a generic
    point to 70 digits — so a typo in a recurrence coefficient reddens HERE,
    not as a disagreement with production."""
    with mp.workdps(_DPS):
        x = mp.mpf("0.3141592653589793238462643383279502884197")
        p, _ = _legendre_bonnet(n, x)
        if abs(p - mp.legendre(n, x)) > mp.mpf(10) ** -70 * max(1, abs(p)):
            pytest.fail(f"Bonnet P_{n} disagrees with mpmath.legendre at x={x}")
        for a, b in ((0, 0), (2, 2), (1.5, 2.25), (-0.5, 0.75)):
            pj = _jacobi_dlmf(n, mp.mpf(a), mp.mpf(b), x)
            exact = mp.jacobi(n, a, b, x)
            if abs(pj - exact) > mp.mpf(10) ** -70 * max(1, abs(exact)):
                pytest.fail(f"DLMF 18.9.2 P_{n}^({a},{b}) disagrees with mpmath.jacobi")


def _moment(label: str, k: int) -> _Real:
    """Closed-form k-th moment of each family's weight."""
    if label == "legendre":
        return mp.mpf(2) / (k + 1) if k % 2 == 0 else mp.mpf(0)
    if label == "chebyshev_t":  # int x^k (1-x^2)^(-1/2) = B((k+1)/2, 1/2) for even k
        return mp.beta(mp.mpf(k + 1) / 2, mp.mpf(1) / 2) if k % 2 == 0 else mp.mpf(0)
    if label == "chebyshev_u":
        return mp.beta(mp.mpf(k + 1) / 2, mp.mpf(3) / 2) if k % 2 == 0 else mp.mpf(0)
    if label == "hermite":
        return mp.gamma(mp.mpf(k + 1) / 2) if k % 2 == 0 else mp.mpf(0)
    if label.startswith("laguerre"):
        a = mp.mpf(label[len("laguerre("):-1])
        return mp.gamma(a + k + 1)
    a, b = (mp.mpf(s) for s in label[len("jacobi("):-1].split(","))
    # x^k = sum_j C(k,j) (1+x)^j (-1)^(k-j); int (1-x)^a (1+x)^(b+j) = 2^(a+b+j+1) B(a+1, b+j+1)
    return mp.fsum(
        mp.binomial(k, j) * (-1) ** (k - j) * 2 ** (a + b + j + 1) * mp.beta(a + 1, b + j + 1)
        for j in range(k + 1)
    )


@pytest.mark.foundation
@pytest.mark.parametrize("family", [f[0] for f in _FAMILIES])
@pytest.mark.parametrize("n", [2, 3, 5, 8])
def test_reference_rules_integrate_their_moments_exactly(family: str, n: int) -> None:
    """Every reference rule (80 digits, before rounding) integrates
    x^0..x^{2n-1} against its OWN weight to 60 digits relative to the moment
    scale — the defining property of the Gauss rule, checked against
    closed-form moments that share nothing with the node or weight formulas."""
    rule = next(f[2] for f in _FAMILIES if f[0] == family)
    with mp.workdps(_DPS):
        nodes, weights = rule(n)
        if len({mp.nstr(x, 40) for x in nodes}) != n:
            pytest.fail(f"{family} n={n}: Newton returned repeated roots")
        for k in range(2 * n):
            quad = mp.fsum(w * x**k for x, w in zip(nodes, weights))
            exact = _moment(family, k)
            scale = mp.fsum(abs(w) * abs(x) ** k for x, w in zip(nodes, weights))
            if abs(quad - exact) > mp.mpf(10) ** -60 * scale:
                pytest.fail(f"{family} n={n}: moment x^{k} = {quad} != {exact}")


# ── the gate: production is the correctly rounded rule, bit for bit ─────────


@pytest.mark.l1
@pytest.mark.rests_on(_REF_SELF_CHECK, _REF_MOMENTS)
@pytest.mark.parametrize(("measure", "rule", "mass", "n"), _cases())
def test_rule_is_the_correctly_rounded_rule(measure, rule, mass, n: int) -> None:
    """``measure.gauss(n)`` == the 80-digit reference rounded once, bit for bit."""
    produced = measure.gauss(n)
    nodes, weights = reference_rule(rule, n)
    same_nodes = np.array_equal(produced.nodes, nodes)
    same_weights = np.array_equal(produced.weights, weights)
    if not (same_nodes and same_weights):
        pytest.fail(
            f"{measure.name} n={n} is not correctly rounded: nodes off by up to "
            f"{_ulps(produced.nodes, nodes)} ULP, weights by up to "
            f"{_ulps(produced.weights, weights)} ULP (max|dw|/max|w| = "
            f"{np.max(np.abs(produced.weights - weights)) / np.max(weights):.2e})"
        )


@pytest.mark.l1
@pytest.mark.rests_on(f"{_HERE}::test_rule_is_the_correctly_rounded_rule")
@pytest.mark.parametrize(("measure", "rule", "mass", "n"), _cases(even_only=True))
def test_even_weight_rule_is_exactly_symmetric(measure, rule, mass, n: int) -> None:
    """THEOREM: for an even weight the exact rule satisfies x_i = -x_{n-1-i},
    w_i = w_{n-1-i}, and rounding to nearest commutes with negation, so the
    correctly rounded rule is EXACTLY mirror-symmetric with an exact 0.0 at
    the centre of an odd rule — with no imposition step. This row is what
    licenses retiring the mirror average in ``gauss``: an angular quadrature's
    reflection closure (``Mirror("x")``) is then integer index arithmetic."""
    produced = measure.gauss(n)
    x, w = produced.nodes, produced.weights
    if not (np.array_equal(x, -x[::-1]) and np.array_equal(w, w[::-1])):
        pytest.fail(
            f"{measure.name} n={n}: rule not exactly symmetric "
            f"(node defect {np.max(np.abs(x + x[::-1])):.2e}, "
            f"weight defect {np.max(np.abs(w - w[::-1])):.2e})"
        )
    if n % 2 == 1 and x[n // 2] != 0.0:
        pytest.fail(f"{measure.name} n={n}: centre node {x[n // 2]!r} is not exactly 0")


@pytest.mark.l1
@pytest.mark.rests_on(f"{_HERE}::test_rule_is_the_correctly_rounded_rule")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 8, 16, 33, 64])
@pytest.mark.parametrize(
    ("a", "b", "constant"),
    [(0.0, 0.0, LEGENDRE), (-0.5, -0.5, CHEBYSHEV_T), (0.5, 0.5, CHEBYSHEV_U)],
    ids=["jacobi(0,0)=legendre", "jacobi(-.5,-.5)=chebyshev_t", "jacobi(.5,.5)=chebyshev_u"],
)
def test_jacobi_specialisation_is_its_constant_bit_for_bit(a, b, constant, n: int) -> None:
    """THEOREM of the construction: a correctly rounded real is unique, so the
    general Jacobi recurrence at (0,0) / (-1/2,-1/2) / (1/2,1/2) and the
    specialised recurrence of the constant produce the SAME bytes (the module
    docstring's "bit-identically" claim, which Golub–Welsch made true only by
    accident of identical float inputs)."""
    general, special_ = jacobi(a, b).gauss(n), constant.gauss(n)
    if not (
        np.array_equal(general.nodes, special_.nodes)
        and np.array_equal(general.weights, special_.weights)
    ):
        pytest.fail(
            f"jacobi({a},{b}).gauss({n}) != {constant.name}.gauss({n}): nodes "
            f"{_ulps(general.nodes, special_.nodes)} ULP, weights "
            f"{_ulps(general.weights, special_.weights)} ULP apart"
        )


@pytest.mark.l1
@pytest.mark.rests_on(f"{_HERE}::test_rule_is_the_correctly_rounded_rule")
@pytest.mark.parametrize(("measure", "rule", "mass", "n"), _cases())
def test_weight_sum_is_within_one_float_step_of_the_mass(measure, rule, mass, n: int) -> None:
    r"""THEOREM (why the zeroth-moment renormalisation retires). With
    :math:`\hat w_i = w_i(1+\delta_i)`, :math:`|\delta_i|\le u = 2^{-53}`, the
    EXACT sum of the rounded weights satisfies
    :math:`|\sum\hat w_i-\mu_0|\le u\mu_0<\operatorname{ulp}(\mu_0)`, and
    ``math.fsum`` returns that sum correctly rounded; so ``fsum(w)`` lies
    within ONE float step of the correctly rounded mass. A renormalisation
    ``w * (mu0 / sum(w))`` would re-round every weight and break the equality
    row above. RECORD on top of the theorem (``[M]`` 2026-10-01, the
    reference): for LEGENDRE the sum is EXACTLY 2.0 at every order in the
    population; for the others it is within the theorem's one step (e.g.
    chebyshev_t n=3: one step above pi)."""
    weights = measure.gauss(n).weights
    with mp.workdps(_DPS):
        mu0 = _round_once(mass())
    total = math.fsum(weights)
    steps = {mu0, math.nextafter(mu0, math.inf), math.nextafter(mu0, -math.inf)}
    if total not in steps:
        pytest.fail(f"{measure.name} n={n}: fsum(w) = {total!r}, mass = {mu0!r}")
    if measure == LEGENDRE and total != 2.0:
        pytest.fail(f"legendre n={n}: fsum(w) = {total!r}, recorded exactly 2.0")


@pytest.mark.foundation
@pytest.mark.parametrize("measure", [LEGENDRE, jacobi(1.5, 2.25)], ids=["legendre", "jacobi(1.5,2.25)"])
def test_a_returned_rule_cannot_corrupt_the_next_one(measure: GeneratingMeasure) -> None:
    """The (family, n) cache (ruled 2026-10-01) must hand out FRESH or
    read-only arrays: a caller writing into one returned rule must not change
    the rule the next caller gets. Green before the cache exists; its tooth is
    a cache that returns its stored arrays by reference."""
    first = measure.gauss(8)
    pristine = first.weights.copy(), first.nodes.copy()
    for array in (first.weights, first.nodes):
        if array.flags.writeable:
            array[0] = 99.0
    second = measure.gauss(8)
    if not (np.array_equal(second.weights, pristine[0]) and np.array_equal(second.nodes, pristine[1])):
        pytest.fail(f"{measure.name}.gauss(8): a write into one returned rule changed the next call's rule")
