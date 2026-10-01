r"""Gate (c): the homogeneous k∞, flux and condensed cross sections against the EXACT answer of their float inputs, within a DERIVED forward-error bound.

Why this replaces a byte pin
----------------------------
``test_byte_stability`` froze ``k_inf`` (as hex), the flux (as bytes) and
``sig_prod``/``sig_abs`` (as hex) for 8 mixtures. A frozen byte of a LAPACK
output pins the platform: macOS 27.0.1's Accelerate moved ``geev`` by 1 ULP
on 2 of the 8 cases with no ORPHEUS change (``[M]``
``scratch/platform_drift/memo.md``). The property worth pinning is that
production's answer is the EXACT answer to the problem it was handed, to
within the rounding its algorithm is entitled to, and that is what this gate
asserts. The exact answer is computed in rational arithmetic by
:mod:`orpheus.derivations.common.exact_homogeneous` (L0; ``fractions`` only;
self-certified against :math:`F\varphi = kA\varphi` with zero residual).

The premise: what production computes (the rank-one route, ruled 2026-10-01)
-----------------------------------------------------------------------------
1. :math:`\hat A = \mathrm{fl}(\operatorname{diag}\Sigma_t - (\Sigma_{s0}^T + 2\Sigma_{2,0}^T))`,
   the hub's ``loss.as_matrix()``;
2. :math:`\hat u` = one LU solve (``getrf`` with partial pivoting, then
   ``getrs``) of :math:`\hat A u = \chi` — ``MatrixInverseOperator(loss).apply(chi)``;
3. :math:`\hat k = \mathrm{fl}\langle\nu\Sigma_f, \hat u\rangle` (a :math:`G`-term reduction);
4. :math:`\hat\varphi = \hat u \cdot \mathrm{fl}(100 / \hat n)`, with :math:`\hat n` a
   :math:`G`-term reduction of :math:`\langle\nu\Sigma_f,\hat u\rangle` (the ``ScaleGauge``);
5. :math:`\hat\sigma_x = \mathrm{fl}(\hat N_x / \hat D)`,
   :math:`\hat N_x = \mathrm{fl}\langle\Sigma_x,\hat\varphi\rangle`,
   :math:`\hat D = \mathrm{fl}\langle 1,\hat\varphi\rangle`.

The bound (every constant is derived; :math:`u = 2^{-53}`, :math:`\gamma_m = mu/(1-mu)`)
------------------------------------------------------------------------------------
*Assembly.* :math:`E = \hat A - A` is computed EXACTLY (both are rationals).

*Solve* — Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.
(2002), Theorem 9.4: the computed solution satisfies
:math:`(\hat A + \Delta)\hat u = \chi` with
:math:`|\Delta| \le \gamma_{3G}\,P^T|\hat L||\hat U|`, where :math:`\hat L, \hat U`
are the COMPUTED factors. The gate reads them from ``scipy.linalg.lu_factor``
of the same :math:`\hat A` in the same process (the factors production's
``MatrixInverseOperator`` holds), so the bound is a posteriori in the actual
factors, valid on any LAPACK that implements the standard algorithm. Hence
:math:`(A + E + \Delta)\hat u = \chi` and, with
:math:`D = |E| + \gamma_{3G}P^T|\hat L||\hat U|` and :math:`W = |A^{-1}|D \ge 0`,

.. math::

   |\hat u - u| = |A^{-1}(E+\Delta)\hat u| \le W(|u| + |\hat u - u|)
   \;\Longrightarrow\;
   e := |\hat u - u| \le (I-W)^{-1}W|u| \le W|u| + \tfrac{r}{1-r}\,\|W|u|\|_\infty\mathbf 1,

with :math:`r = \|W\|_\infty < 1` (asserted; the Neumann series of a
non-negative matrix of norm :math:`r` majorises the inverse). This is the
componentwise form of Higham's Theorem 7.4 (forward error from a componentwise
backward error), with no first-order truncation.

*Dot product* — Higham Eq. (3.5): a :math:`G`-term inner product in any
summation order satisfies :math:`|\mathrm{fl}(a^Tb) - a^Tb| \le
\gamma_G|a|^T|b|`. The gate uses :math:`\gamma_{G+1}`, admitting one extra
rounding per term for the counting measure's weight multiply (exact when the
weight is 1.0, as it is here; the extra rounding is the stated price of not
depending on that). So

.. math::

   |\hat k - k| \le B_k := |\nu\Sigma_f|^T e + \gamma_{G+1}|\nu\Sigma_f|^T(|u|+e).

*Gauge.* With :math:`\eta = B_k/k` and two roundings in
:math:`\hat u_i\cdot\mathrm{fl}(100/\hat n)`:
:math:`|\hat\varphi_i - \varphi_i| \le B_{\varphi,i} := \tfrac{100}{k}\big[e_i\tfrac{1+\gamma_2}{1-\eta}
+ |u_i|\,\rho\big]`, :math:`\rho = \max\big(\tfrac{1+\gamma_2}{1-\eta}-1,\;
1-\tfrac{1-\gamma_2}{1+\eta}\big)`.

*Condensed cross sections.* Every term is non-negative (the cone), so the
relative errors add: with :math:`r_N = B_N/N`, :math:`r_D = B_D/D` from the
dot-product bound on :math:`\hat\varphi`, the ratio's relative error is at most
:math:`\max\big(\tfrac{(1+r_N)(1+u)}{1-r_D}-1,\;1-\tfrac{(1-r_N)(1-u)}{1+r_D}\big)`.

Evaluated per case (``[M]`` 2026-10-01, bound / ULP of the exact value):

==================  =====  ===============  ======  ======
case                 k      flux (per group)  σ_prod  σ_abs
==================  =====  ===============  ======  ======
homo_1eg, A_1g       3.75   5.2              18.8    12.5
homo_2eg (+eg, A_2g) 18.4   24.0 34.4        71.1    75.2
homo_2eg_n2n         17.7   23.7 30.5        45.1    76.9
homo_4eg, A_4g       43.2   47.3 79.5 52.4 69.8  173.4  115.1
==================  =====  ===============  ======  ======

Measured deviations on the same population: HEAD's ``geev`` route
:math:`\le 1.55` ULP in k; the rank-one route :math:`\le 0.86` ULP; every flux
component within 7 % of its bound. A worst-case bound is loose by
construction (it holds for EVERY rounding pattern); its value is that a
platform cannot move it.

Resolution, and what is out of it (honestly)
--------------------------------------------
A 1-ULP perturbation of one χ entry moves the exact k by at most ~1 ULP,
which is INSIDE every case's bound: this gate cannot see it, and no gate
built on a rigorous bound for a float LU solve can (the rounding of the
solve alone is entitled to more). ``[M]`` battery (``gates_spec.md`` §c): a
χ scale of :math:`(1 + m\,u)` reds at the smallest ``m`` listed there, per
case; the structural mutations (a dropped dot-product term, the (n,2n) term
dropped from :math:`A`, a transposed scattering matrix, a wrong gauge
target) red by orders of magnitude. The bit-level platform witness for the
homogeneous route is therefore NOT this gate; it is the absence of any
platform-dependent primitive on the path except one LU solve, whose bytes
remain platform-dependent and are deliberately not pinned.

Claim kind: REFERENCE (exact rational arithmetic, structurally independent:
the reference shares only the float inputs with production).
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import lu_factor

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.exact_homogeneous import (
    ExactInfiniteMedium,
    exact_infinite_medium_of,
)
from orpheus.homogeneous.solver import (
    HomogeneousProblem,
    HomogeneousResult,
    solve_homogeneous_infinite,
)

from ._homogeneous_population import mixture_cases

_HERE = "tests/gates/homogeneous/test_kinf_exact_reference.py"

#: Unit roundoff of IEEE binary64, exact.
_U = Fraction(1, 2**53)


def _gamma(m: int) -> Fraction:
    return m * _U / (1 - m * _U)


def _exact(x: float) -> Fraction:
    return Fraction(float(x))


def _ulp(x: Fraction) -> Fraction:
    return Fraction(math.ulp(float(x)))


def _lu_abs_product(a_hat: np.ndarray) -> list[list[Fraction]]:
    """:math:`P^T|\\hat L||\\hat U|` of the COMPUTED LU factors, exactly."""
    lu, piv = lu_factor(np.array(a_hat, dtype=float))
    g = lu.shape[0]
    lower = [[Fraction(1) if i == j else (_exact(lu[i, j]) if i > j else Fraction(0)) for j in range(g)] for i in range(g)]
    upper = [[_exact(lu[i, j]) if j >= i else Fraction(0) for j in range(g)] for i in range(g)]
    product = [[sum((abs(lower[i][m]) * abs(upper[m][j]) for m in range(g)), Fraction(0)) for j in range(g)] for i in range(g)]
    row_of = list(range(g))  # LAPACK's ipiv: row i was swapped with row piv[i], in order
    for i, p in enumerate(piv):
        row_of[i], row_of[p] = row_of[p], row_of[i]
    out = [[Fraction(0)] * g for _ in range(g)]
    for i in range(g):
        out[row_of[i]] = product[i]
    return out


class _Bounds:
    """The derived forward-error bounds of one case (module docstring)."""

    def __init__(self, ref: ExactInfiniteMedium, a_hat: np.ndarray) -> None:
        g = len(ref.chi)
        assembly = [[abs(_exact(a_hat[i, j]) - ref.loss[i][j]) for j in range(g)] for i in range(g)]
        lu_term = _lu_abs_product(a_hat)
        backward = [[assembly[i][j] + _gamma(3 * g) * lu_term[i][j] for j in range(g)] for i in range(g)]
        w = [
            [sum((abs(ref.loss_inverse[i][m]) * backward[m][j] for m in range(g)), Fraction(0)) for j in range(g)]
            for i in range(g)
        ]
        self.r = max(sum(row, Fraction(0)) for row in w)
        if self.r >= 1:
            raise ArithmeticError(f"||W||_inf = {float(self.r)} >= 1: the bound does not apply")
        u_abs = [abs(v) for v in ref.u]
        w_u = [sum((w[i][j] * u_abs[j] for j in range(g)), Fraction(0)) for i in range(g)]
        tail = self.r / (1 - self.r) * max(w_u)
        self.e = [w_u[i] + tail for i in range(g)]
        nsf = [abs(v) for v in ref.nu_sigma_f]
        dot = _gamma(g + 1)
        self.k = sum((nsf[i] * self.e[i] for i in range(g)), Fraction(0)) + dot * sum(
            (nsf[i] * (u_abs[i] + self.e[i]) for i in range(g)), Fraction(0)
        )
        eta = self.k / ref.k_inf
        g2 = _gamma(2)
        rho = max((1 + g2) / (1 - eta) - 1, 1 - (1 - g2) / (1 + eta))
        self.flux = [
            Fraction(100) / ref.k_inf * (self.e[i] * (1 + g2) / (1 - eta) + u_abs[i] * rho)
            for i in range(g)
        ]
        self.sig_prod = self._ratio(ref, ref.nu_sigma_f, dot)
        self.sig_abs = self._ratio(ref, ref.sigma_a, dot)

    def _ratio(self, ref: ExactInfiniteMedium, sigma: list[Fraction], dot: Fraction) -> Fraction:
        phi, b = ref.flux, self.flux
        s = [abs(v) for v in sigma]
        num = sum((si * p for si, p in zip(s, phi)), Fraction(0))
        den = sum(phi, Fraction(0))
        b_num = sum((si * bi for si, bi in zip(s, b)), Fraction(0)) + dot * sum(
            (si * (p + bi) for si, p, bi in zip(s, phi, b)), Fraction(0)
        )
        b_den = sum(b, Fraction(0)) + dot * sum((p + bi for p, bi in zip(phi, b)), Fraction(0))
        r_num, r_den = b_num / num, b_den / den
        rel = max((1 + r_num) * (1 + _U) / (1 - r_den) - 1, 1 - (1 - r_num) * (1 - _U) / (1 + r_den))
        return rel * num / den


def _case(name: str) -> tuple[ExactInfiniteMedium, _Bounds, HomogeneousResult]:
    mix = mixture_cases()[name]
    if not isinstance(mix, Mixture):
        raise TypeError(f"{name}: the population helper returned {type(mix).__name__}, not a Mixture")
    ref = exact_infinite_medium_of(mix)
    a_hat = np.asarray(HomogeneousProblem(mix).loss.as_matrix(), dtype=float)
    return ref, _Bounds(ref, a_hat), solve_homogeneous_infinite(mix)


def _within(value: float, exact: Fraction, bound: Fraction) -> bool:
    return abs(_exact(value) - exact) <= bound


def _report(label: str, value: float, exact: Fraction, bound: Fraction) -> str:
    unit = _ulp(exact)
    return (
        f"{label} = {value!r} is {float((_exact(value) - exact) / unit):+.2f} ULP from the "
        f"exact {float(exact)!r}; the derived bound is {float(bound / unit):.2f} ULP"
    )


_CASES = sorted(mixture_cases())
_REFERENCE_ROW = f"{_HERE}::test_the_exact_reference_certifies_itself"


@pytest.mark.foundation
@pytest.mark.parametrize("case", _CASES)
def test_the_exact_reference_certifies_itself(case: str) -> None:
    """``certify`` re-checks A A⁻¹ = I, F φ = k A φ, tr(A⁻¹F) = k, the gauge and the
    cone in exact arithmetic (it also runs at construction; this row makes it a
    reported gate, and asserts the bound's precondition r = ‖W‖∞ < 1/2)."""
    ref, bounds, _ = _case(case)
    ref.certify()
    if bounds.r >= Fraction(1, 2):
        pytest.fail(f"{case}: ||W||_inf = {float(bounds.r):.3e}: the solve is too ill-conditioned for the bound")


@pytest.mark.foundation
@pytest.mark.parametrize("case", _CASES)
def test_the_acceptance_predicate_accepts_the_rounded_exact_value_and_rejects_twice_the_bound(case: str) -> None:
    """Decoder control (X1), two legs per quantity: the correctly rounded exact
    value is ACCEPTED (its error, half an ULP, is below every bound), and the
    float nearest ``exact ± 2·bound`` is REJECTED (it is at least
    ``2·bound − ulp/2 > bound`` away, since every bound exceeds one ULP)."""
    ref, bounds, _ = _case(case)
    for label, exact, bound in (
        ("k_inf", ref.k_inf, bounds.k),
        ("sig_prod", ref.sig_prod, bounds.sig_prod),
        ("sig_abs", ref.sig_abs, bounds.sig_abs),
        *((f"flux[{g}]", x, b) for g, (x, b) in enumerate(zip(ref.flux, bounds.flux))),
    ):
        if bound <= _ulp(exact):
            pytest.fail(f"{case}.{label}: bound {float(bound / _ulp(exact)):.2f} ULP is below one ULP")
        if not _within(float(exact), exact, bound):
            pytest.fail(f"{case}.{label}: the rounded exact value is rejected")
        for side in (1, -1):
            if _within(float(exact + side * 2 * bound), exact, bound):
                pytest.fail(f"{case}.{label}: a value two bounds away is accepted")


@pytest.mark.l1
@pytest.mark.rests_on(_REFERENCE_ROW)
@pytest.mark.parametrize("case", _CASES)
def test_kinf_is_the_exact_eigenvalue_within_the_derived_bound(case: str) -> None:
    ref, bounds, result = _case(case)
    if not _within(result.k_inf, ref.k_inf, bounds.k):
        pytest.fail(f"{case}: " + _report("k_inf", result.k_inf, ref.k_inf, bounds.k))


@pytest.mark.l1
@pytest.mark.rests_on(_REFERENCE_ROW)
@pytest.mark.parametrize("case", _CASES)
def test_flux_is_the_exact_gauged_eigenvector_within_the_derived_bound(case: str) -> None:
    ref, bounds, result = _case(case)
    flux = np.asarray(result.flux, dtype=float).ravel()
    if flux.shape != (len(ref.flux),):
        pytest.fail(f"{case}: flux shape {flux.shape} is not ({len(ref.flux)},)")
    for g, (value, exact, bound) in enumerate(zip(flux, ref.flux, bounds.flux)):
        if not _within(float(value), exact, bound):
            pytest.fail(f"{case}: " + _report(f"flux[{g}]", float(value), exact, bound))


@pytest.mark.l1
@pytest.mark.rests_on(_REFERENCE_ROW)
@pytest.mark.parametrize("case", _CASES)
def test_condensed_cross_sections_are_exact_within_the_derived_bound(case: str) -> None:
    ref, bounds, result = _case(case)
    for label, value, exact, bound in (
        ("sig_prod", result.sig_prod, ref.sig_prod, bounds.sig_prod),
        ("sig_abs", result.sig_abs, ref.sig_abs, bounds.sig_abs),
    ):
        if not _within(float(value), exact, bound):
            pytest.fail(f"{case}: " + _report(label, float(value), exact, bound))
