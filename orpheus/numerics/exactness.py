r"""What a quadrature rule's exactness claim is ABOUT.

A degree of exactness on its own is half a claim. The whole claim is

.. math::
   :label: exactness-claim

   \sum_i w_i\, f(x_i) \;=\; \int_{\mathcal{X}} f \, d\lambda
   \qquad \text{for every } f \in V,

and it names **two** things the bare integer does not: the reference
measure :math:`d\lambda` the two sides are compared against, and the
finite-dimensional space :math:`V` the claim quantifies over. The
degree is an *index into* :math:`V`, meaningless without it.

Why this module exists
----------------------

The tree carried the degree and the reference as two independent fields
on :class:`~orpheus.numerics.measure.DiscreteMeasure`, so a degree could
exist without saying what it was about. Two measured consequences:

* **Different measures, same number.** `[M]` Gauss-Legendre and
  Gauss-Chebyshev both ship ``degree_of_exactness = 2n - 1``. At
  :math:`n = 4`, :math:`q = x^6`, GC reproduces
  :math:`\int q\,(1-x^2)^{-1/2}dx` exactly and misses the *unweighted*
  integral by **0.696**. One field, two incompatible claims.
* **Different SPACES, same number.** `[M]` The :math:`n`-node equispaced
  rule on an interval is the midpoint rule, of *algebraic* degree
  :math:`1`. The same nodes on the **circle** are the periodic
  trapezoid, exact for *trigonometric* polynomials of degree
  :math:`n - 1`. Same nodes, same weights, and the integer cannot tell
  the two claims apart — which is why composing a polar rule with an
  azimuthal one by ``min()`` on the raw integers yields **degree 1**
  when the truth is :math:`\min(2n_\mu - 1,\, n_\varphi - 1)`.

The second is the one that blocks a principled product rule, and it is
not fixable by correcting a number: **what is missing is the space
attached to the degree.**

The unification
---------------

:math:`V` is not free-floating. For every rule in this tree it is the
span of the first :math:`m` functions of the reference measure's own
**orthogonal system** — and *which* system that is follows from the
measure:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - reference measure
     - orthogonal system
     - degree indexes
   * - weight :math:`w` on an interval
     - orthogonal polynomials
     - algebraic degree
   * - uniform on the circle
     - the Fourier basis
     - trigonometric degree
   * - uniform on :math:`S^2`
     - spherical harmonics
     - :math:`\ell`

So "the generating measure IS the exactness space" holds in general
once the system is read off the measure rather than assumed to be
polynomials. The three rows are genuinely non-isomorphic function
systems, which is why :class:`OrthogonalSystem` is a value a claim
carries and not a detail a caller may assume.

Reference vs generator — a deliberate separation
------------------------------------------------

:class:`~orpheus.numerics.generating_measure.GeneratingMeasure` is
defined by a three-term recurrence: it is the object that **builds** a
Gauss rule (Golub-Welsch). That is strictly more than an exactness
claim needs, and strictly less than every rule has:

* the circle's Fourier system satisfies no three-term recurrence of
  that kind, yet the periodic trapezoid has a perfectly good exactness
  claim;
* Lebedev's claim is about the uniform measure on :math:`S^2`, which
  generates nothing by Golub-Welsch — its authority is a published
  table.

So this module types the claim against :class:`ReferenceMeasure` — *a
measure plus its orthogonal system* — and ``GeneratingMeasure``
satisfies that protocol as the sub-case that **also** knows how to
construct a rule from it. Widening ``GeneratingMeasure`` instead would
have welded "what the claim is about" to "how the rule was built",
which is the same welding
:func:`~orpheus.numerics.quadrature.rules_product.product_mu_phi`
already demonstrates the cost of.

See Also
--------
:class:`orpheus.numerics.generating_measure.GeneratingMeasure`
    The Golub-Welsch generator; the sub-case that builds as well as
    certifies.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol, runtime_checkable

from .measure import SPACE_CIRCLE, SPACE_SPHERE, Space


class OrthogonalSystem(Enum):
    r"""Which family of functions a degree of exactness indexes.

    A **value**, not a hierarchy: the systems differ in what they span,
    not in behaviour, and nothing dispatches on them beyond deciding
    whether two claims are comparable. (The project rule is to mint a
    type only when there are :math:`\ge 2` non-isomorphic realizations
    *and* a non-identity morphism is applied to it; here the morphism —
    the product theorem — acts on the *claim*, not on the system.)

    Attributes
    ----------
    ALGEBRAIC
        Orthogonal polynomials of the reference weight. Degree
        :math:`p` means every polynomial of degree :math:`\le p` is
        integrated exactly *against that weight*. The Gauss families.
    TRIGONOMETRIC
        The Fourier basis on the circle. Degree :math:`p` means every
        trigonometric polynomial :math:`\sum_{|k| \le p} c_k e^{ik\phi}`
        is integrated exactly. The periodic trapezoid.
    SPHERICAL_HARMONIC
        Real spherical harmonics on :math:`S^2`. Degree :math:`\ell`
        means every :math:`Y_{\ell'}^m` with :math:`\ell' \le \ell` is
        integrated exactly. Lebedev, and what a level-symmetric rule
        claims.
    """

    ALGEBRAIC = "algebraic"
    TRIGONOMETRIC = "trigonometric"
    SPHERICAL_HARMONIC = "spherical_harmonic"


@runtime_checkable
class ReferenceMeasure(Protocol):
    r"""The continuous measure an exactness claim is **about**.

    Deliberately narrower than
    :class:`~orpheus.numerics.generating_measure.GeneratingMeasure`:
    a claim needs to know *what* it is compared against and *which*
    functions it quantifies over, and needs nothing at all about how
    the rule was constructed. Structural (``Protocol``) rather than
    inherited, so a generator satisfies it by having the attributes
    rather than by being told to.
    """

    @property
    def name(self) -> str:
        """Canonical mathematical identity — what equality compares."""
        ...

    @property
    def support(self) -> Space:
        """The domain the measure lives on."""
        ...

    @property
    def orthogonal_system(self) -> OrthogonalSystem:
        """Which function family a degree against this measure indexes."""
        ...


@dataclass(frozen=True)
class UniformMeasure:
    r"""Lebesgue measure on a domain — a reference that generates nothing.

    The second realization of :class:`ReferenceMeasure`, and the one that
    shows why the protocol is not just
    :class:`~orpheus.numerics.generating_measure.GeneratingMeasure` under
    another name. A uniform measure is what the non-Gauss rules are exact
    against:

    * the **midpoint / equispaced** rule on an interval — exact for
      algebraic polynomials of degree 1, against :math:`dx`;
    * the **periodic trapezoid** on the circle — exact for trigonometric
      polynomials of degree :math:`n-1`, against :math:`d\varphi`;
    * **Lebedev** on :math:`S^2` — exact for spherical harmonics, against
      :math:`d\Omega`.

    None of the three is produced by a three-term recurrence, and the
    last two are not polynomial claims at all. The orthogonal system is
    passed explicitly rather than inferred from :attr:`support`, because
    inferring it would mean parsing a string tag to decide a
    mathematical fact.

    Attributes
    ----------
    support : Space
        The domain :math:`\mathcal{X}`.
    orthogonal_system : OrthogonalSystem
        Which family the harmonic analysis of :math:`\mathcal{X}`
        provides — polynomials on an interval, Fourier on the circle,
        spherical harmonics on :math:`S^2`.
    """

    support: Space
    orthogonal_system: OrthogonalSystem

    @property
    def name(self) -> str:
        """Canonical identity — two uniform measures on the same domain
        with the same system ARE the same measure, so the name is
        derived rather than stored and cannot disagree with the fields."""
        return f"uniform({self.support})"


#: Lebesgue measure on the circle :math:`[0, 2\pi)`. The reference the
#: periodic trapezoid is exact against, and the one whose orthogonal
#: system is the **Fourier basis** — so a degree against it is a
#: *trigonometric* degree, not the algebraic degree the same nodes would
#: carry when read as a midpoint rule on an interval.
UNIFORM_ON_CIRCLE = UniformMeasure(
    support=SPACE_CIRCLE, orthogonal_system=OrthogonalSystem.TRIGONOMETRIC,
)

#: Lebesgue measure on :math:`S^2`. The reference for cubatures whose
#: claim is about spherical harmonics (Lebedev, and what a
#: level-symmetric rule claims) rather than about a weight on an
#: interval.
UNIFORM_ON_SPHERE = UniformMeasure(
    support=SPACE_SPHERE,
    orthogonal_system=OrthogonalSystem.SPHERICAL_HARMONIC,
)


@dataclass(frozen=True)
class ProductMeasure:
    r"""The product reference :math:`\lambda_1 \otimes \lambda_2`.

    A tensor product of rules lands on the **product space**, so its
    claim is about the product measure — not about either factor's. That
    distinction is what separates a tensor product from a direct sum, and
    conflating them is a live bug this type exists to prevent: a shared
    "combined degree" helper served both operations, so a product could
    inherit a factor's reference and thereby claim exactness against a
    measure it is not exact against.

    Only defined when the factors share an :class:`OrthogonalSystem`. A
    mixed product — algebraic :math:`\times` trigonometric, which is
    exactly the polar :math:`\times` azimuthal case — spans a *mixed
    tensor* system whose degree is not a minimum of the factors' but a
    theorem about the target space. Such a product legitimately has no
    generic claim here, and the theorem belongs with the rule that
    applies it.
    """

    factors: tuple[ReferenceMeasure, ...]

    def __post_init__(self) -> None:
        if len(self.factors) < 2:
            raise ValueError(
                f"a product measure needs >= 2 factors, got "
                f"{len(self.factors)}"
            )
        systems = {f.orthogonal_system for f in self.factors}
        if len(systems) != 1:
            raise ValueError(
                f"a ProductMeasure's factors must share an orthogonal "
                f"system; got {sorted(s.value for s in systems)}. A mixed "
                f"product spans a mixed tensor system whose degree is a "
                f"theorem about the target space, not a minimum — build "
                f"that claim where the theorem lives."
            )

    @property
    def name(self) -> str:
        return " ⊗ ".join(f.name for f in self.factors)

    @property
    def support(self) -> Space:
        return " × ".join(f.support for f in self.factors)

    @property
    def orthogonal_system(self) -> OrthogonalSystem:
        """The shared system — guaranteed unique by construction."""
        return self.factors[0].orthogonal_system


@dataclass(frozen=True)
class ExactnessClaim:
    r"""A quadrature rule's exactness, stated whole: :eq:`exactness-claim`.

    Bundling the degree with its reference is the point. As two loose
    fields either could be read without the other, and both failure
    modes in this module's header are exactly that: a number compared
    against the wrong measure, or against the wrong function space.
    Here the degree is unreachable without the thing it indexes.

    Attributes
    ----------
    reference : ReferenceMeasure
        The measure :math:`d\lambda` the claim is against, carrying the
        orthogonal system the degree indexes.
    degree : int
        The largest :math:`p` such that every element of the reference's
        orthogonal system up to index :math:`p` is integrated exactly.
        Must be :math:`\ge 0` — a rule that integrates not even the
        constant exactly has no claim to make, and should carry ``None``
        rather than a negative degree.
    """

    reference: ReferenceMeasure
    degree: int

    def __post_init__(self) -> None:
        if self.degree < 0:
            raise ValueError(
                f"an exactness claim needs degree >= 0, got {self.degree}; "
                f"a rule with no exact subspace carries no claim (None), "
                f"rather than a negative one"
            )

    @property
    def system(self) -> OrthogonalSystem:
        """The function system this claim quantifies over — read from the
        reference, never stored separately, so the two cannot drift."""
        return self.reference.orthogonal_system

    def tensor_with(self, other: "ExactnessClaim") -> "ExactnessClaim | None":
        r"""The claim of a TENSOR PRODUCT of the two rules, or ``None``.

        The product lands on the **product space**, so its reference is
        :class:`ProductMeasure` — never either factor's. For separable
        integrands of degree :math:`p` per axis both factors must be
        exact to :math:`p`, so the degree is :math:`\min(p_1, p_2)`.

        ⭐ **This makes a correct claim representable that the tree
        previously had to suppress.** The old shared "combined degree"
        helper refused any pair with different references, on the
        strength of a measurement that `[M]`
        ``gauss_legendre(4) * gauss_chebyshev(4)`` advertised degree 7
        while integrating the constant :math:`1` to **6.2832** instead
        of :math:`4`. Re-measured with the reference named: the product
        **is** exact to degree 7 per axis against
        :math:`dx \otimes (1-y^2)^{-1/2}dy` — worst error **4.16e-13**
        over every :math:`(a,b) \le 7`, with the degree-8 control
        missing by 1.5e-2, so the bound is tight. And :math:`2\pi =
        6.2832` **is** the correct mass of ``legendre ⊗ chebyshev_t``;
        the old expectation of :math:`4` was the *Lebesgue* product,
        which is not this product's reference. The refusal was a
        conservative workaround for a missing type, not a law.

        Returns ``None`` when the factors' systems differ — a mixed
        product (algebraic :math:`\times` trigonometric, i.e. polar
        :math:`\times` azimuthal) spans a mixed tensor system whose
        degree comes from a theorem about the target space, not from a
        minimum. That theorem belongs with the rule that applies it.

        .. note:: Why there is no ``combined_with``

           An earlier draft carried a meet "for operations on one space",
           intended for the direct sum. It is **not** the direct sum's
           law: :meth:`~orpheus.numerics.measure.DiscreteMeasure.__add__`
           requires equal supports, so summing two rules for
           :math:`\lambda` gives a rule for :math:`2\lambda` — keeping
           the shared reference would assert exactness against a measure
           the sum is not exact against, the very error this type exists
           to prevent. The direct sum therefore carries no claim, and the
           meet had no other consumer.
        """
        if self.system is not other.system:
            return None
        return ExactnessClaim(
            reference=ProductMeasure(
                factors=(self.reference, other.reference),
            ),
            degree=min(self.degree, other.degree),
        )

    def __str__(self) -> str:
        return (
            f"exact to {self.system.value} degree {self.degree} "
            f"against {self.reference.name}"
        )


__all__ = [
    "UNIFORM_ON_CIRCLE",
    "UNIFORM_ON_SPHERE",
    "ExactnessClaim",
    "OrthogonalSystem",
    "ReferenceMeasure",
    "UniformMeasure",
]
