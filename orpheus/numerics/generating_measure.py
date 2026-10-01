r"""The measure that GENERATES a Gauss rule — and that its exactness is ABOUT.

Every Gaussian quadrature in the classical zoo — Legendre, Chebyshev,
Jacobi, Laguerre, Hermite — is **one construction** applied to different
data. This module is that construction, plus the data.

The construction (Golub & Welsch 1969)
--------------------------------------

A measure :math:`d\lambda = w(x)\,dx` on an interval fixes a family of
orthogonal polynomials, which obey a three-term recurrence

.. math::
   :label: generating-measure-three-term

   p_{k+1}(x) = (x - \alpha_k)\, p_k(x) - \beta_k\, p_{k-1}(x),
   \qquad p_{-1} \equiv 0, \quad p_0 \equiv 1.

Assemble the symmetric tridiagonal **Jacobi matrix**

.. math::
   :label: generating-measure-jacobi-matrix

   J_n = \begin{pmatrix}
     \alpha_0 & \sqrt{\beta_1} & & \\
     \sqrt{\beta_1} & \alpha_1 & \ddots & \\
     & \ddots & \ddots & \sqrt{\beta_{n-1}} \\
     & & \sqrt{\beta_{n-1}} & \alpha_{n-1}
   \end{pmatrix}

then the :math:`n`-point Gauss rule for :math:`d\lambda` is

.. math::
   :label: generating-measure-golub-welsch

   x_i = \lambda_i(J_n), \qquad
   w_i = \mu_0 \, \bigl[v_i\bigr]_1^{\,2},

the eigenvalues of :math:`J_n` and the squared **first components** of
its unit eigenvectors, scaled by the zeroth moment
:math:`\mu_0 = \int d\lambda`.

The rule is correctly rounded
-----------------------------

The equation above DEFINES the rule; it is not how the floats are produced.
Every node and weight shipped is the float64 nearest its exact value, so the
rule is a property of the measure alone. Until 2026-10-01 the floats were the
platform ``eigh``'s output, and `[M]` a macOS update (Accelerate, 2026-09-30)
moved the GL-4 rule by 1-2 ULP and reddened 24 bit pins on an unchanged tree;
Golub-Welsch output also sits up to 445 ULP from the correctly rounded weights
at :math:`n = 64`. The construction now:

1. takes float Golub-Welsch eigenvalues as SEEDS only;
2. refines each to the root of the orthonormal polynomial :math:`\hat p_n` by
   Newton on the recurrence, in ``mpmath`` arithmetic, and reads its weight as
   the Christoffel number :math:`w_i = 1/\sum_{k<n}\hat p_k(x_i)^2`
   (Gautschi 2004, Thm 1.46: the same quantity as
   :math:`\mu_0[v_i]_1^2`);
3. does so at two working precisions (40 and 60 digits) and rounds each once;
   the two must agree bit for bit, and every value must round to ONE float64
   across its whole error interval (a radius :math:`10^{10-d}` relative at
   :math:`d` digits) — that interval check is what certifies "correctly
   rounded", since two precisions can agree on the same wrong side of a tie.
   Either failure raises: a value that cannot be certified is refused, never
   guessed;
4. for an even weight computes the non-negative half and mirrors it, with an
   odd rule's centre at exactly zero: the exact rule is symmetric, and
   rounding commutes with negation.

The recurrence coefficients are written ONCE, in arbitrary precision
(:attr:`GeneratingMeasure.exact_recurrence`); the float coefficients
(:meth:`GeneratingMeasure.recurrence`) are their correctly rounded image. The
zeroth-moment renormalisation the float construction ended with is gone: it
would un-round correctly rounded weights. Rules are cached per
``(measure, recurrence, n)``: keying on the recurrence too keeps two
constructions of one measure (``jacobi(0, 0)`` and :data:`LEGENDRE`) computed
separately, so their agreement stays a cross-check rather than a cache read.

So a family is nothing but :math:`(\alpha_k, \beta_k, \mu_0)`. It has no
behaviour to override, which is why the families here are **values, not
subclasses** — a subclass whose entire content is returning three arrays
is ceremony around data (see ``coding-elegance``, "the instance with
extra ceremony"). One generic body consumes them all, and it is
measured 2–6× *faster* than the specialised routines it replaces, so
there is not even a performance argument for an override hook.

.. list-table:: The classical families
   :header-rows: 1
   :widths: 20 18 30 32

   * - Family
     - Support
     - Weight :math:`w(x)`
     - :math:`\mu_0 = \int w\,dx`
   * - Legendre
     - :math:`[-1,1]`
     - :math:`1`
     - :math:`2`
   * - Chebyshev-1 (:math:`T`)
     - :math:`[-1,1]`
     - :math:`(1-x^2)^{-1/2}`
     - :math:`\pi`
   * - Chebyshev-2 (:math:`U`)
     - :math:`[-1,1]`
     - :math:`(1-x^2)^{+1/2}`
     - :math:`\pi/2`
   * - Jacobi :math:`(a,b)`
     - :math:`[-1,1]`
     - :math:`(1-x)^a (1+x)^b`
     - :math:`2^{a+b+1} B(a{+}1, b{+}1)`
   * - Laguerre :math:`(a)`
     - :math:`[0,\infty)`
     - :math:`x^a e^{-x}`
     - :math:`\Gamma(a+1)`
   * - Hermite
     - :math:`\mathbb{R}`
     - :math:`e^{-x^2}`
     - :math:`\sqrt{\pi}`

Jacobi **generalises** the first three: ``jacobi(0, 0)`` is Legendre,
``jacobi(-1/2, -1/2)`` is Chebyshev-1, ``jacobi(1/2, 1/2)`` is
Chebyshev-2. That is not a remark — it is a *verification instrument*.
Two independent constructions of the same measure must agree, and
``tests/gates/numerics/test_generating_measure.py`` asserts they do: both
round to the same correctly rounded rule, so nodes AND weights agree bit for
bit. The campaign's
acceptance rule is that ≥2 realizations *prove* an implementation rather
than merely pinning it; here the second realization is free.

Why the exactness claim lives here
----------------------------------

A Gauss rule's degree of exactness is a claim **with respect to a
measure**, and the measure is this one:

.. math::

   \sum_{i=1}^{n} w_i\, q(x_i) \;=\; \int q(x)\, w(x)\, dx
   \qquad \text{for all } \deg q \le 2n-1.

The weight function is *in the integral*, not in the quadrature. Before
2026-08-02 ``gauss_legendre`` and ``gauss_chebyshev`` both shipped
``degree_of_exactness = 2n-1`` in the same field meaning different
things — Chebyshev's docstring said so in prose while the type erased
it. `[M]` at :math:`n=4`, :math:`q(x)=x^6`: the Chebyshev rule
reproduces :math:`\int q (1-x^2)^{-1/2} dx` to machine precision and
misses the *unweighted* integral by **0.696**.

Hence :attr:`~orpheus.numerics.measure.DiscreteMeasure.generating_measure`:
a rule built from its measure carries that measure, so its exactness
claim states what it is about. This is ``coding-elegance`` Pattern 4
applied to quadrature — **a rule built from its measure cannot lie about
its exactness space**, because the exactness follows from the
construction rather than being an integer someone typed.

The converse is the diagnostic. A rule with no generating measure
(``lebedev_sphere``, ``level_symmetric_sn``) is not thereby wrong — its
claim simply rests on external authority (a published table plus a
citation) rather than on a construction in this codebase. That is a real
distinction and it is why the field is optional. It is also how
issue #327 was possible: ``level_symmetric_sn`` assigned one weight to
every ordinate by hand, so nothing constrained its advertised degree,
and `[M]` it WAS degree-3 at *every* order while claiming :math:`N-1`
(fixed 2026-08-06 — the weights are now per-orbit moment-solved and
the advertised degree is gated both directions).

What is deliberately NOT here
-----------------------------

**Node constraints** (Gauss / Radau / Lobatto) are *not* families. They
are the same construction with prescribed nodes, obtained by modifying
the trailing entries of :math:`J_n` (Golub 1973). They belong on a
second, orthogonal axis — a closed sum type applied as a morphism on the
rule — and are not built here because nothing consumes them yet.

References
----------

* Golub, G.H. and Welsch, J.H. (1969). "Calculation of Gauss quadrature
  rules." *Mathematics of Computation* **23**(106), 221-230.
* Golub, G.H. (1973). "Some modified matrix eigenvalue problems."
  *SIAM Review* **15**(2), 318-334. (Radau / Lobatto as a rank-one
  modification of :math:`J_n`.)
* Gautschi, W. (2004). *Orthogonal Polynomials: Computation and
  Approximation*. Oxford. §1.3 for the recurrence coefficients of every
  family above; the :math:`\beta_0 := \mu_0` storage convention used
  here is his.
* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. §3.6 (Gauss quadrature, exactness theorem 3.6.20).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import Callable

import mpmath
import numpy as np

from .exactness import ExactnessClaim, OrthogonalSystem
from .measure import DiscreteMeasure
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    HALF_LINE,
    Interval,
    Manifold,
    REAL_LINE,
)

# Type of a recurrence-coefficient generator: (n, digits) -> (alpha, beta),
# each a list of n ``mpmath.mpf`` evaluated at ``digits`` decimal digits, with
# the storage convention beta[0] = mu_0. The precision is an ARGUMENT, not
# mpmath's ambient setting, so a bare call cannot return 15-digit values under
# the name "exact". The closed forms are written once, in arbitrary
# precision; the float coefficients and every Gauss rule are derived from
# them, so the two cannot drift.
RecurrenceCoefficients = Callable[[int, int], "tuple[list, list]"]

# The two working precisions (decimal digits) a rule is computed at. Each
# result is rounded once to float64 and the two must agree bit for bit (a
# disagreement means some exact value is unresolved at the coarser
# precision), and the finer result must also round to one float64 across its
# whole error interval (see ``_CERTIFIED_DIGITS``). Either failure is refused,
# never resolved silently.
_WORKING_DIGITS = (40, 60)

# The precision the float views (``recurrence``, ``mass``) and the symmetry
# probe evaluate the closed forms at: far beyond float64, so each float is
# the correctly rounded image of its closed form.
_VIEW_DIGITS = 40

# How many of the fine run's digits are certified: Newton stops when its step
# is below ``10**(5 - digits)`` (relative), and converges quadratically, so
# the node is accurate far below that; the Christoffel sum at that node loses
# at most a few digits more. A radius of ``10**(_CERTIFIED_DIGITS - digits)``
# (relative) is therefore a generous error interval, still ~35 orders of
# magnitude below a float64 spacing at 60 digits. `[M]` 2026-10-01 (qa
# review): 40-digit values differ from 100-digit recomputations by at most
# 1.0e-37 relative, against the 1e-30 radius, over Laguerre n=128,
# Laguerre(1.5) n=200, Hermite n=128, Jacobi(-1/2, 3/4) n=128 and Legendre
# n=256. A value is shipped only when both ends of its interval round to the
# same float64.
_CERTIFIED_DIGITS = 10

# How many Newton steps a node may take from its float seed. A seed from the
# float Golub-Welsch rule carries ~16 correct digits and Newton doubles them
# per step, so 60 digits need three; the cap only bounds a pathological seed.
_NEWTON_STEP_LIMIT = 50

# How many recurrence coefficients to inspect when deciding whether a
# measure is symmetric. `alpha_k == 0` either holds for all k or fails at
# the first: for every classical family alpha is identically zero or
# grows monotonically in |k|, so a short prefix decides it. Eight is
# generous — jacobi(a, b) with a != b already separates at k = 0.
_SYMMETRY_PROBE_ORDER = 8


@dataclass(frozen=True)
class GeneratingMeasure:
    r"""A continuous measure :math:`d\lambda = w(x)\,dx` that generates a
    Gauss rule, and that the rule's exactness claim is about.

    The defining data is the three-term recurrence
    :eq:`generating-measure-three-term`. Everything else — nodes,
    weights, the degree of exactness, the total mass — follows from it,
    which is the point: none of them can drift from the measure they
    describe.

    Attributes
    ----------
    name : str
        Canonical mathematical identity, e.g. ``"legendre"``,
        ``"jacobi(a=1.5, b=2.25)"``. This is what equality compares —
        two *constructions* of the same measure are the same measure, so
        ``jacobi(0, 0) == LEGENDRE`` is ``True`` even though the two
        carry different :attr:`exact_recurrence` callables and reach the
        coefficients by different code paths. (That the two paths agree
        is a separate, testable claim, and it is tested: both round to the
        same correctly rounded rule, bit for bit.)
    support : Manifold
        The interval the measure lives on.
    exact_recurrence : callable
        ``(n, digits) -> (alpha, beta)``, each a list of ``n``
        ``mpmath.mpf`` evaluated at ``digits`` decimal digits, following Gautschi's
        storage convention :math:`\beta_0 := \mu_0`. The closed forms live
        here ONCE, in arbitrary precision; :meth:`recurrence` (the float
        view) and :meth:`gauss` (the rule) are derived from it. Carrying
        the zeroth moment *inside* the recurrence rather than beside it
        means the mass and the family that defines it cannot disagree.
        Excluded from equality and ``repr``: it is the implementation of
        the identity that :attr:`name` states, not part of it.
    """

    name: str
    support: Manifold
    exact_recurrence: RecurrenceCoefficients = field(compare=False, repr=False)

    # -- derived quantities -------------------------------------------

    def recurrence(self, n: int) -> "tuple[np.ndarray, np.ndarray]":
        r"""The first ``n`` recurrence coefficients as float64, each
        correctly rounded from :attr:`exact_recurrence`.

        The float view of the one definition, so a consumer reading
        coefficients and a consumer reading the rule read the same numbers.
        """
        alpha, beta = self.exact_recurrence(n, _VIEW_DIGITS)
        return (
            np.array([float(a) for a in alpha]),
            np.array([float(b) for b in beta]),
        )

    @property
    def is_symmetric(self) -> bool:
        r"""Is the weight even, :math:`w(-x) = w(x)`?

        **Derived, never declared** — it is a theorem that

        .. math::

           \alpha_k \equiv 0 \quad \Longleftrightarrow \quad
           w \text{ is even about the origin},

        because :math:`\alpha_k = \langle x p_k, p_k\rangle /
        \langle p_k, p_k\rangle` is the first moment of an even
        function against an odd integrand when :math:`w` is even.

        scipy carries this as a hand-set ``symmetrize`` boolean passed
        to its generic routine (``_gen_roots_and_weights``). Reading it
        off the recurrence instead costs nothing and cannot fall out of
        step with the family it describes — the same reason
        :attr:`mass` is read rather than stored. `[M]` It agrees with
        scipy's hand-set flag on every family shipped here, including
        the parameterised ones: ``jacobi(a, b)`` is symmetric exactly
        when ``a == b``, which the derivation gets right without being
        told.
        """
        alpha, _ = self.exact_recurrence(_SYMMETRY_PROBE_ORDER, _VIEW_DIGITS)
        return all(a == 0 for a in alpha)

    @property
    def mass(self) -> float:
        r""":math:`\mu_0 = \int w\,dx`, the zeroth moment, correctly rounded.

        Read from the recurrence rather than stored, so it is the mass
        of *this* family by construction.
        """
        _, beta = self.recurrence(1)
        return float(beta[0])

    @property
    def orthogonal_system(self) -> OrthogonalSystem:
        r"""Always :attr:`~orpheus.numerics.exactness.OrthogonalSystem.ALGEBRAIC`
        — and that is a **theorem of the construction, not a choice**.

        A three-term recurrence
        :eq:`generating-measure-three-term` generates a sequence of
        *polynomials*, degree :math:`k` at index :math:`k`, orthogonal
        with respect to :math:`w`. So a measure defined by such a
        recurrence has algebraic polynomials as its orthogonal system by
        definition; there is no measure of this class whose degree could
        index anything else.

        This property is what makes :class:`GeneratingMeasure` satisfy
        :class:`~orpheus.numerics.exactness.ReferenceMeasure` — the
        broader protocol an exactness claim is typed against. Systems
        with no such recurrence (the Fourier basis on the circle, the
        spherical harmonics) are reference measures that are **not**
        generating measures, which is exactly why the claim is typed
        against the protocol rather than against this class.
        """
        return OrthogonalSystem.ALGEBRAIC


    # -- the construction ---------------------------------------------

    def gauss(self, n: int) -> DiscreteMeasure:
        r"""The :math:`n`-point Gauss rule for this measure, correctly rounded.

        Each node and weight is the float64 nearest the exact value
        (:eq:`generating-measure-golub-welsch` names them), so the rule is
        a property of the measure and of nothing else: not of the platform's
        linear-algebra library, not of its thread count, not of its version.
        The construction is :func:`_correctly_rounded_rule`; this method
        wraps its cached numbers in a fresh measure, so no caller can mutate
        a shared array. The cache key includes the recurrence itself, not
        only the measure's identity: ``jacobi(0, 0) == LEGENDRE``, but the
        two constructions are computed separately, so their agreement stays
        a cross-check rather than a cache read.

        The returned measure carries ``generating_measure=self``, so its
        ``degree_of_exactness = 2n - 1`` states which integral it is
        exact for. For a weighted family that matters: the claim is
        about :math:`\int q\,w\,dx`, never about :math:`\int q\,dx`.

        Parameters
        ----------
        n : int
            Number of nodes; must be :math:`\ge 1`.

        Returns
        -------
        DiscreteMeasure
            On this measure's :attr:`support`, with
            ``degree_of_exactness = 2n - 1`` and
            ``generating_measure = self``.

        Raises
        ------
        ArithmeticError
            When the two working precisions round some node or weight
            differently, when a value's error interval straddles a float64
            rounding boundary (it cannot be certified), or when Newton fails
            to converge from the float seed.
        """
        if n < 1:
            raise ValueError(f"{self.name}.gauss requires n >= 1, got n={n}")
        nodes, weights = _correctly_rounded_rule(self, self.exact_recurrence, n)
        return DiscreteMeasure(
            nodes=np.array(nodes),
            weights=np.array(weights),
            support=self.support,
            # The claim names its own reference: this rule is exact to
            # algebraic degree 2n-1 against THIS measure, and (for a
            # weighted family) against no other. A rule built from its
            # measure cannot over-claim, because the degree follows from
            # the construction rather than being typed in beside it.
            exactness=ExactnessClaim(reference=self, degree=2 * n - 1),
        )

    # -- morphisms ----------------------------------------------------

    def on(self, a: float, b: float) -> "GeneratingMeasure":
        r"""Push this measure forward along the affine map
        :math:`[-1,1] \to [a,b]`.

        Under :math:`x = \tfrac{1}{2}\bigl[(b-a)t + (a+b)\bigr]` the
        recurrence transforms as

        .. math::

           \alpha'_k = \tfrac{1}{2}\bigl[(b-a)\alpha_k + (a+b)\bigr],
           \qquad
           \beta'_0 = \tfrac{b-a}{2}\,\beta_0,
           \qquad
           \beta'_k = \Bigl(\tfrac{b-a}{2}\Bigr)^2 \beta_k \;\; (k \ge 1),

        the shift and scale of the nodes, the Jacobian on the mass, and
        the square of the scale on the off-diagonal (which enters
        :math:`J` under a square root, so it scales linearly there).
        The transform is applied in the working precision, so the remapped
        rule is correctly rounded too (the float endpoints are exact
        binary numbers and enter exactly).

        Defined only for measures on :math:`[-1,1]`; the unbounded
        families have no finite interval to remap onto.
        """
        return _affine_remap(self, float(a), float(b))


@functools.lru_cache(maxsize=None)
def _affine_remap(measure: GeneratingMeasure, a: float, b: float) -> GeneratingMeasure:
    """:meth:`GeneratingMeasure.on`, cached per ``(measure, a, b)`` so a
    repeated remap returns the same measure and hits the rule cache."""
    if measure.support != COSINE_INTERVAL:
        raise ValueError(
            f"affine remap is defined for measures on "
            f"{COSINE_INTERVAL.name}, but {measure.name} lives on "
            f"{measure.support.name}"
        )
    if not a < b:
        raise ValueError(f"require a < b, got a={a}, b={b}")
    inner = measure.exact_recurrence

    def remapped(n: int, digits: int) -> "tuple[list, list]":
        alpha, beta = inner(n, digits)
        with mpmath.workdps(digits):
            scale = (mpmath.mpf(b) - mpmath.mpf(a)) / 2
            shift = (mpmath.mpf(a) + mpmath.mpf(b)) / 2
            alpha = [scale * x + shift for x in alpha]
            beta = [scale * beta[0]] + [scale**2 * x for x in beta[1:]]
        return alpha, beta

    return GeneratingMeasure(
        name=f"{measure.name}_on[{a},{b}]",
        # The remapped family lives on the interval itself, not on a string
        # that spells one — ``Interval`` is the same type the source support
        # already is, so ``affine_remap`` composes with itmeasure.
        support=Interval(a, b),
        exact_recurrence=remapped,
    )


# ---------------------------------------------------------------------------
# The construction: Newton on the orthonormal recurrence, rounded once
# ---------------------------------------------------------------------------


def _orthonormal_at(x, alpha, beta, n: int):
    r"""The orthonormal :math:`\hat p_n(x)`, its derivative, and the
    Christoffel sum :math:`\sum_{k<n} \hat p_k(x)^2`, by the recurrence.

    :math:`\sqrt{\beta_{k+1}}\,\hat p_{k+1} = (x - \alpha_k)\,\hat p_k -
    \sqrt{\beta_k}\,\hat p_{k-1}`, :math:`\hat p_0 = \beta_0^{-1/2}`
    (Gautschi 2004, Eq. 1.3.13); the derivative obeys the differentiated
    recurrence. Needs ``alpha``, ``beta`` of length ``n + 1``.
    """
    p_previous, p = mpmath.mpf(0), 1 / mpmath.sqrt(beta[0])
    dp_previous, dp = mpmath.mpf(0), mpmath.mpf(0)
    christoffel_sum = p * p
    for k in range(n):
        coupling_next = mpmath.sqrt(beta[k + 1])
        coupling = mpmath.sqrt(beta[k]) if k > 0 else mpmath.mpf(0)
        p_next = ((x - alpha[k]) * p - coupling * p_previous) / coupling_next
        dp_next = (p + (x - alpha[k]) * dp - coupling * dp_previous) / coupling_next
        p_previous, p, dp_previous, dp = p, p_next, dp, dp_next
        if k < n - 1:
            christoffel_sum += p * p
    return p, dp, christoffel_sum


def _float_seeds(measure: GeneratingMeasure, n: int) -> np.ndarray:
    r"""Golub-Welsch in float64: the eigenvalues of :math:`J_n`, ascending.

    Only a STARTING POINT for Newton — accurate to a few ULP on any
    platform, never returned. What the platform's ``eigh`` gives here cannot
    reach the rule, because Newton converges every seed to the same exact
    root and the result is rounded from that.
    """
    alpha, beta = measure.recurrence(n)
    if n == 1:
        return alpha[:1].copy()
    off_diagonal = np.sqrt(beta[1:n])
    jacobi_matrix = (
        np.diag(alpha) + np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    )
    return np.linalg.eigvalsh(jacobi_matrix)


def _rule_at(
    measure: GeneratingMeasure, n: int, seeds: np.ndarray, digits: int
) -> "tuple[tuple[float, ...], tuple[float, ...]]":
    r"""The rule computed at ``digits`` working digits, rounded once.

    Each node is the root of :math:`\hat p_n` reached by Newton from its
    seed; its weight is the Christoffel number
    :math:`w_i = 1 / \sum_{k<n} \hat p_k(x_i)^2` (Gautschi 2004, Thm 1.46),
    the same quantity Golub-Welsch reads off the eigenvector. For an even
    weight only the non-negative half is computed and the other half is its
    mirror, with the centre of an odd rule at exactly zero: the exact rule
    has that symmetry, and rounding commutes with negation, so the rounded
    rule has it bit for bit.
    """
    alpha, beta = measure.exact_recurrence(n + 1, digits)
    with mpmath.workdps(digits):
        tolerance = mpmath.mpf(10) ** (5 - digits)
        symmetric = measure.is_symmetric
        roots = []
        for seed in (seeds[n // 2:] if symmetric else seeds):
            if symmetric and n % 2 == 1 and len(roots) == 0:
                roots.append(mpmath.mpf(0))
                continue
            x = mpmath.mpf(float(seed))
            for _ in range(_NEWTON_STEP_LIMIT):
                p, dp, _ = _orthonormal_at(x, alpha, beta, n)
                step = p / dp
                x -= step
                if abs(step) <= tolerance * max(1, abs(x)):
                    break
            else:
                raise ArithmeticError(
                    f"{measure.name}.gauss({n}): Newton did not converge from "
                    f"the seed {float(seed)!r} in {_NEWTON_STEP_LIMIT} steps"
                )
            roots.append(x)
        weights = [1 / _orthonormal_at(x, alpha, beta, n)[2] for x in roots]
        if symmetric:
            mirrored = [-x for x in reversed(roots) if x != 0]
            roots = mirrored + roots
            weights = list(reversed(weights[len(weights) - len(mirrored):])) + weights
        nodes = tuple(_certified_float(x, digits, measure, n) for x in roots)
        if any(b <= a for a, b in zip(nodes, nodes[1:])):
            raise ArithmeticError(
                f"{measure.name}.gauss({n}): Newton converged two seeds to one "
                f"root, or out of order: {nodes}"
            )
        return nodes, tuple(_certified_float(w, digits, measure, n) for w in weights)


def _certified_float(value, digits: int, measure: "GeneratingMeasure", n: int) -> float:
    r"""``value`` rounded to float64, certified: both ends of its error
    interval :math:`[v - \rho|v|, v + \rho|v|]`,
    :math:`\rho = 10^{\texttt{_CERTIFIED_DIGITS} - \texttt{digits}}`, round to
    the same float64, so the exact value does too. A value within
    :math:`\rho` of a rounding boundary is refused rather than guessed."""
    with mpmath.workdps(digits):
        radius = abs(value) * mpmath.mpf(10) ** (_CERTIFIED_DIGITS - digits)
        low, high = float(value - radius), float(value + radius)
    if low != high:
        raise ArithmeticError(
            f"{measure.name}.gauss({n}): {float(value)!r} lies within its error "
            f"radius of a float64 rounding boundary, so it cannot be certified "
            f"correctly rounded"
        )
    return low


@functools.lru_cache(maxsize=None)
def _correctly_rounded_rule(
    measure: GeneratingMeasure, recurrence: RecurrenceCoefficients, n: int
) -> "tuple[tuple[float, ...], tuple[float, ...]]":
    r"""The correctly rounded :math:`n`-point rule, as tuples (immutable, so
    the cache can hand them out).

    Computed at each of the two working precisions and rounded once; the
    two must agree bit for bit, and each value is certified against its error
    interval. Cached per ``(measure, recurrence, n)``: ``recurrence`` is
    ``measure.exact_recurrence``, passed so the key names the CONSTRUCTION and
    not only the measure's identity (equal measures reached by different
    recurrences are computed separately, which is what keeps the
    ``jacobi(0, 0)`` / Legendre agreement a cross-check).
    """
    seeds = _float_seeds(measure, n)
    coarse, fine = (_rule_at(measure, n, seeds, d) for d in _WORKING_DIGITS)
    if coarse != fine:
        raise ArithmeticError(
            f"{measure.name}.gauss({n}): the rule rounds differently at "
            f"{_WORKING_DIGITS[0]} and {_WORKING_DIGITS[1]} digits, so some exact "
            f"node or weight sits on a float64 rounding tie"
        )
    return fine


# ---------------------------------------------------------------------------
# The classical families
# ---------------------------------------------------------------------------
#
# Coefficients: Gautschi 2004 §1.3, each written once in arbitrary precision
# (mpmath arithmetic at the precision ``_at_digits`` sets). Each family's rule
# is verified against an independent correctly rounded reference in
# tests/gates/numerics/test_gauss_rules_correctly_rounded.py, and the
# parameterised families are cross-checked against the constants they
# specialise to in tests/gates/numerics/test_generating_measure.py.


def _at_digits(closed_form: "Callable[[int], tuple[list, list]]") -> RecurrenceCoefficients:
    """Lift a closed form written in mpmath arithmetic to a
    :data:`RecurrenceCoefficients`: evaluate it at the precision the caller
    names, never at mpmath's ambient one."""

    @functools.wraps(closed_form)
    def recurrence(n: int, digits: int) -> "tuple[list, list]":
        with mpmath.workdps(digits):
            return closed_form(n)

    return recurrence


@_at_digits
def _legendre_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.mpf(2)] + [mpmath.mpf(k * k) / (4 * k * k - 1) for k in range(1, n)]
    return alpha, beta[:n]


@_at_digits
def _chebyshev_t_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    # beta_1 = 1/2 breaks the otherwise-constant 1/4 pattern; this is the
    # only irregular coefficient in the family.
    beta = [mpmath.pi, mpmath.mpf(1) / 2] + [mpmath.mpf(1) / 4] * max(n - 2, 0)
    return alpha, beta[:n]


@_at_digits
def _chebyshev_u_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.pi / 2] + [mpmath.mpf(1) / 4] * (n - 1)
    return alpha, beta


@_at_digits
def _hermite_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.sqrt(mpmath.pi)] + [mpmath.mpf(k) / 2 for k in range(1, n)]
    return alpha, beta


#: Weight :math:`w(x) = 1` on :math:`[-1,1]`. The unweighted family —
#: its Gauss rule is exact for plain polynomial integration.
LEGENDRE = GeneratingMeasure(
    name="legendre",
    support=COSINE_INTERVAL,
    exact_recurrence=_legendre_recurrence,
)

#: Weight :math:`w(x) = (1-x^2)^{-1/2}` on :math:`[-1,1]`.
#: **Its exactness is about the WEIGHTED integral** — the rule does not
#: integrate bare polynomials on :math:`[-1,1]`.
CHEBYSHEV_T = GeneratingMeasure(
    name="chebyshev_t",
    support=COSINE_INTERVAL,
    exact_recurrence=_chebyshev_t_recurrence,
)

#: Weight :math:`w(x) = (1-x^2)^{+1/2}` on :math:`[-1,1]`.
CHEBYSHEV_U = GeneratingMeasure(
    name="chebyshev_u",
    support=COSINE_INTERVAL,
    exact_recurrence=_chebyshev_u_recurrence,
)

#: Weight :math:`w(x) = e^{-x^2}` on :math:`\mathbb{R}`.
HERMITE = GeneratingMeasure(
    name="hermite",
    support=REAL_LINE,
    exact_recurrence=_hermite_recurrence,
)


def jacobi(a: float, b: float) -> GeneratingMeasure:
    """The Jacobi measure for exponents ``(a, b)``; see :func:`_jacobi`.

    Cached per ``(float(a), float(b))``, so a repeated call returns the same
    measure, with the same recurrence, and hits the rule cache.
    """
    return _jacobi(float(a), float(b))


@functools.lru_cache(maxsize=None)
def _jacobi(a: float, b: float) -> GeneratingMeasure:
    r"""Weight :math:`w(x) = (1-x)^a (1+x)^b` on :math:`[-1,1]`.

    The parent of the three :math:`[-1,1]` constants above:
    ``jacobi(0, 0)`` is :data:`LEGENDRE`, ``jacobi(-1/2, -1/2)`` is
    :data:`CHEBYSHEV_T`, ``jacobi(1/2, 1/2)`` is :data:`CHEBYSHEV_U`.
    Those specialisations carry the **canonical name** of the family
    they equal, so they compare equal to the constants — while still
    running the general recurrence, which is what makes the agreement a
    genuine cross-check rather than an alias.

    Parameters
    ----------
    a, b : float
        Exponents; both must exceed :math:`-1` for the weight to be
        integrable. They enter the working precision exactly (a float64 is
        a binary rational).
    """
    if a <= -1.0 or b <= -1.0:
        raise ValueError(
            f"jacobi requires a > -1 and b > -1 for an integrable weight, "
            f"got a={a}, b={b}"
        )

    @_at_digits
    def recurrence(n: int) -> "tuple[list, list]":
        A, B = mpmath.mpf(a), mpmath.mpf(b)
        AB = A + B
        # mu_0 = 2^(a+b+1) B(a+1, b+1).
        alpha = [(B - A) / (AB + 2)]
        beta = [mpmath.power(2, AB + 1) * mpmath.beta(A + 1, B + 1)]
        if n > 1:
            # k = 1 is a REMOVABLE singularity of the general beta
            # formula below: its numerator carries (k + ab) and its
            # denominator (2k + ab - 1), and at k = 1 both are (1 + ab).
            # They cancel. Evaluating the general form there would divide
            # by zero whenever a + b = -1 — which is exactly Chebyshev-1,
            # jacobi(-1/2, -1/2), the single most common member.
            alpha.append((B * B - A * A) / ((2 + AB) * (4 + AB)))
            beta.append(4 * (1 + A) * (1 + B) / ((2 + AB) ** 2 * (3 + AB)))
        for k in range(2, n):
            two_k_ab = 2 * k + AB
            alpha.append((B * B - A * A) / (two_k_ab * (two_k_ab + 2)))
            beta.append(
                4 * k * (k + A) * (k + B) * (k + AB)
                / (two_k_ab**2 * (two_k_ab + 1) * (two_k_ab - 1))
            )
        return alpha[:n], beta[:n]

    return GeneratingMeasure(
        name=_jacobi_name(a, b),
        support=COSINE_INTERVAL,
        exact_recurrence=recurrence,
    )


def _jacobi_name(a: float, b: float) -> str:
    """Canonical name: the classical families keep their own names.

    Equality on :class:`GeneratingMeasure` is equality of mathematical
    identity, and ``jacobi(0, 0)`` *is* Legendre — so it must not
    advertise itself as something else merely because it was reached
    through the general constructor.
    """
    for known_a, known_b, known in (
        (0.0, 0.0, LEGENDRE),
        (-0.5, -0.5, CHEBYSHEV_T),
        (0.5, 0.5, CHEBYSHEV_U),
    ):
        if a == known_a and b == known_b:
            return known.name
    return f"jacobi(a={a}, b={b})"


def laguerre(a: float = 0.0) -> GeneratingMeasure:
    """The generalised Laguerre measure for exponent ``a``; see :func:`_laguerre`.

    Cached per ``float(a)``, so a repeated call returns the same measure,
    with the same recurrence, and hits the rule cache.
    """
    return _laguerre(float(a))


@functools.lru_cache(maxsize=None)
def _laguerre(a: float = 0.0) -> GeneratingMeasure:
    r"""Weight :math:`w(x) = x^a e^{-x}` on :math:`[0,\infty)`.

    Parameters
    ----------
    a : float, optional
        Exponent; must exceed :math:`-1`. Default ``0.0`` is the plain
        Laguerre weight :math:`e^{-x}`.
    """
    if a <= -1.0:
        raise ValueError(
            f"laguerre requires a > -1 for an integrable weight, got a={a}"
        )

    @_at_digits
    def recurrence(n: int) -> "tuple[list, list]":
        A = mpmath.mpf(a)
        alpha = [2 * k + A + 1 for k in range(n)]
        beta = [mpmath.gamma(A + 1)] + [k * (k + A) for k in range(1, n)]
        return alpha, beta

    return GeneratingMeasure(
        name="laguerre" if a == 0.0 else f"laguerre(a={a})",
        support=HALF_LINE,
        exact_recurrence=recurrence,
    )


__all__ = [
    "CHEBYSHEV_T",
    "CHEBYSHEV_U",
    "HERMITE",
    "LEGENDRE",
    "GeneratingMeasure",
    "RecurrenceCoefficients",
    "jacobi",
    "laguerre",
]
