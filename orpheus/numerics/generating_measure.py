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
``tests/numerics/test_generating_measure.py`` asserts they do (measured:
nodes agree **bit-identically**, weights to 2.8e-17). The campaign's
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
issue #327 was possible: ``level_symmetric_sn`` assigns one weight to
every ordinate by hand, so nothing constrains its advertised degree, and
`[M]` it is degree-3 at *every* order while claiming :math:`N-1`.

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

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.special import gammaln

from .exactness import OrthogonalSystem
from .measure import (
    SPACE_HALF_LINE,
    SPACE_INTERVAL_M11,
    SPACE_R,
    DiscreteMeasure,
    Space,
)

# Type of a recurrence-coefficient generator: n -> (alpha, beta), each of
# shape (n,), with the storage convention beta[0] = mu_0.
RecurrenceCoefficients = Callable[[int], "tuple[np.ndarray, np.ndarray]"]

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
        carry different :attr:`recurrence` callables and reach the
        coefficients by different code paths. (That the two paths agree
        numerically is a separate, testable claim, and it is tested.)
    support : Space
        The interval the measure lives on.
    recurrence : callable
        ``n -> (alpha, beta)``, each of shape ``(n,)``, following
        Gautschi's storage convention :math:`\beta_0 := \mu_0`. Carrying
        the zeroth moment *inside* the recurrence rather than beside it
        means the mass and the family that defines it cannot disagree.
        Excluded from equality and ``repr``: it is the implementation of
        the identity that :attr:`name` states, not part of it.
    """

    name: str
    support: Space
    recurrence: RecurrenceCoefficients = field(compare=False, repr=False)

    # -- derived quantities -------------------------------------------

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
        alpha, _ = self.recurrence(_SYMMETRY_PROBE_ORDER)
        return bool(np.all(alpha == 0.0))

    @property
    def mass(self) -> float:
        r""":math:`\mu_0 = \int w\,dx`, the zeroth moment.

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
        r"""The :math:`n`-point Gauss rule for this measure.

        Golub-Welsch :eq:`generating-measure-golub-welsch`: diagonalise
        the Jacobi matrix, take eigenvalues as nodes and
        :math:`\mu_0 \, [v_i]_1^2` as weights.

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
        """
        if n < 1:
            raise ValueError(f"{self.name}.gauss requires n >= 1, got n={n}")
        alpha, beta = self.recurrence(n)
        if n == 1:
            # J_1 = (alpha_0): one node at the mean, carrying all the mass.
            nodes = np.array([alpha[0]], dtype=float)
            weights = np.array([beta[0]], dtype=float)
        else:
            off_diagonal = np.sqrt(beta[1:n])
            jacobi_matrix = (
                np.diag(alpha)
                + np.diag(off_diagonal, 1)
                + np.diag(off_diagonal, -1)
            )
            # eigh, not eig: J is symmetric by construction, and the
            # symmetric driver returns REAL eigenvalues in ascending
            # order with orthonormal eigenvectors — all three are
            # properties the Gauss rule needs and none is guaranteed by
            # the general driver.
            nodes, eigenvectors = np.linalg.eigh(jacobi_matrix)
            weights = beta[0] * eigenvectors[0, :] ** 2

            if self.is_symmetric:
                # IMPOSE the reflection symmetry the measure has, rather
                # than inheriting it to within round-off. For an even
                # weight the exact rule satisfies x_i = -x_{n-1-i} and
                # w_i = w_{n-1-i}; averaging the computed rule against
                # its own mirror enforces that BIT-EXACTLY, and for odd
                # n puts an exact 0.0 at the centre.
                #
                # `[M]` The defect goes from 3.3e-16 (n=8) / 8.6e-16
                # (n=32) to exactly 0.0, and the rule gets *more*
                # accurate as a bonus: worst relative error over degrees
                # 0..2n-1 improves 4.4e-16 -> 1.1e-16 at n=8 and
                # 7.9e-16 -> 3.3e-16 at n=32, because the mirror average
                # cancels the antisymmetric part of the eigensolver's
                # error.
                #
                # This matters here more than it would elsewhere: an
                # angular quadrature's invariance group is load-bearing
                # (a reflecting boundary is exactly representable only
                # if the node set is closed under the face reflection),
                # and an EXACTLY symmetric rule makes that a matter of
                # integer index arithmetic rather than of tolerance.
                nodes = (nodes - nodes[::-1]) / 2.0
                weights = (weights + weights[::-1]) / 2.0

        # Impose the zeroth moment too — degree-0 exactness is the one
        # coefficient we know in closed form, so there is no reason to
        # let the eigensolver's round-off decide it. numpy and scipy
        # both end their Gauss routines this way.
        weights = weights * (beta[0] / weights.sum())

        return DiscreteMeasure(
            nodes=nodes,
            weights=weights,
            support=self.support,
            degree_of_exactness=2 * n - 1,
            generating_measure=self,
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

        Defined only for measures on :math:`[-1,1]`; the unbounded
        families have no finite interval to remap onto.
        """
        if self.support != SPACE_INTERVAL_M11:
            raise ValueError(
                f"affine remap is defined for measures on "
                f"{SPACE_INTERVAL_M11}, but {self.name} lives on "
                f"{self.support}"
            )
        if not a < b:
            raise ValueError(f"require a < b, got a={a}, b={b}")
        scale = (b - a) / 2.0
        shift = (a + b) / 2.0
        inner = self.recurrence

        def remapped(n: int) -> "tuple[np.ndarray, np.ndarray]":
            alpha, beta = inner(n)
            alpha = scale * alpha + shift
            beta = beta.copy()
            beta[0] *= scale
            beta[1:] *= scale**2
            return alpha, beta

        return GeneratingMeasure(
            name=f"{self.name}_on[{a},{b}]",
            support=f"[{a},{b}]",
            recurrence=remapped,
        )


# ---------------------------------------------------------------------------
# The classical families
# ---------------------------------------------------------------------------
#
# Coefficients: Gautschi 2004 §1.3.  Each is verified against an
# independent oracle in tests/numerics/test_generating_measure.py
# (numpy leggauss / laggauss / hermgauss, scipy roots_jacobi /
# roots_chebyu, and the Chebyshev-1 closed form), and the parameterised
# families are additionally cross-checked against the constants they
# specialise to.


def _legendre_recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
    alpha = np.zeros(n)
    beta = np.empty(n)
    beta[0] = 2.0
    if n > 1:
        k = np.arange(1, n, dtype=float)
        beta[1:] = k**2 / (4.0 * k**2 - 1.0)
    return alpha, beta


def _chebyshev_t_recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
    alpha = np.zeros(n)
    beta = np.full(n, 0.25)
    beta[0] = np.pi
    if n > 1:
        # beta_1 = 1/2 breaks the otherwise-constant 1/4 pattern; this is
        # the only irregular coefficient in the family.
        beta[1] = 0.5
    return alpha, beta


def _chebyshev_u_recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
    alpha = np.zeros(n)
    beta = np.full(n, 0.25)
    beta[0] = np.pi / 2.0
    return alpha, beta


def _hermite_recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
    alpha = np.zeros(n)
    beta = np.empty(n)
    beta[0] = np.sqrt(np.pi)
    if n > 1:
        beta[1:] = np.arange(1, n, dtype=float) / 2.0
    return alpha, beta


#: Weight :math:`w(x) = 1` on :math:`[-1,1]`. The unweighted family —
#: its Gauss rule is exact for plain polynomial integration.
LEGENDRE = GeneratingMeasure(
    name="legendre",
    support=SPACE_INTERVAL_M11,
    recurrence=_legendre_recurrence,
)

#: Weight :math:`w(x) = (1-x^2)^{-1/2}` on :math:`[-1,1]`.
#: **Its exactness is about the WEIGHTED integral** — the rule does not
#: integrate bare polynomials on :math:`[-1,1]`.
CHEBYSHEV_T = GeneratingMeasure(
    name="chebyshev_t",
    support=SPACE_INTERVAL_M11,
    recurrence=_chebyshev_t_recurrence,
)

#: Weight :math:`w(x) = (1-x^2)^{+1/2}` on :math:`[-1,1]`.
CHEBYSHEV_U = GeneratingMeasure(
    name="chebyshev_u",
    support=SPACE_INTERVAL_M11,
    recurrence=_chebyshev_u_recurrence,
)

#: Weight :math:`w(x) = e^{-x^2}` on :math:`\mathbb{R}`.
HERMITE = GeneratingMeasure(
    name="hermite",
    support=SPACE_R,
    recurrence=_hermite_recurrence,
)


def jacobi(a: float, b: float) -> GeneratingMeasure:
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
        integrable.
    """
    if a <= -1.0 or b <= -1.0:
        raise ValueError(
            f"jacobi requires a > -1 and b > -1 for an integrable weight, "
            f"got a={a}, b={b}"
        )
    ab = a + b

    def recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
        alpha = np.empty(n)
        beta = np.empty(n)
        # mu_0 = 2^(a+b+1) B(a+1, b+1), via log-gamma so large exponents
        # do not overflow on the way to a finite answer.
        beta[0] = np.exp(
            (ab + 1.0) * np.log(2.0)
            + gammaln(a + 1.0)
            + gammaln(b + 1.0)
            - gammaln(ab + 2.0)
        )
        alpha[0] = (b - a) / (ab + 2.0)
        if n > 1:
            # k = 1 is a REMOVABLE singularity of the general beta
            # formula below: its numerator carries (k + ab) and its
            # denominator (2k + ab - 1), and at k = 1 both are (1 + ab).
            # They cancel. Evaluating the general form there would divide
            # by zero whenever a + b = -1 — which is exactly Chebyshev-1,
            # jacobi(-1/2, -1/2), the single most common member.
            alpha[1] = (b * b - a * a) / ((2.0 + ab) * (4.0 + ab))
            beta[1] = (
                4.0 * (1.0 + a) * (1.0 + b)
                / ((2.0 + ab) ** 2 * (3.0 + ab))
            )
        for k in range(2, n):
            two_k_ab = 2.0 * k + ab
            alpha[k] = (b * b - a * a) / (two_k_ab * (two_k_ab + 2.0))
            beta[k] = (
                4.0 * k * (k + a) * (k + b) * (k + ab)
                / (two_k_ab**2 * (two_k_ab + 1.0) * (two_k_ab - 1.0))
            )
        return alpha, beta

    return GeneratingMeasure(
        name=_jacobi_name(a, b),
        support=SPACE_INTERVAL_M11,
        recurrence=recurrence,
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

    def recurrence(n: int) -> "tuple[np.ndarray, np.ndarray]":
        k = np.arange(n, dtype=float)
        alpha = 2.0 * k + a + 1.0
        beta = np.empty(n)
        beta[0] = float(np.exp(gammaln(a + 1.0)))
        if n > 1:
            kk = k[1:]
            beta[1:] = kk * (kk + a)
        return alpha, beta

    return GeneratingMeasure(
        name="laguerre" if a == 0.0 else f"laguerre(a={a})",
        support=SPACE_HALF_LINE,
        recurrence=recurrence,
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
