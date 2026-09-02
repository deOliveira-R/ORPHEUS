r"""1-D quadrature rules on the polar cosine interval :math:`[-1, 1]`,
and their adoption onto the orbit space :math:`S^2/O(2)_a`.

Two objects live here, and the distinction is the point (tracker 2.4 of
the angular-spaces campaign, #429, 2026-09-01):

* :func:`gauss_legendre_on_mu` — the CHART-level rule. A Gauss rule for
  Lebesgue measure on :math:`[-1, 1]`, tagged with the :math:`\sigma_x`
  invariance group (the reflection :math:`\mu \to -\mu`, which under the
  tree's canonical :math:`(\mu, 0, 0)` embedding is the mirror with
  normal :math:`\hat e_x`). It names NO orbit space, because it serves
  two: it is the polar factor of every product rule in
  :mod:`.rules_product`, whose :math:`\mu` is the cosine against
  :math:`z`, and the raw material of the slab's rule, whose :math:`\mu`
  is the cosine against :math:`x`. A factor that declared one axis inside
  a product about the other would be a false claim, so the chart-level
  rule declares neither.
* :func:`gauss_legendre_on_polar_orbit` — the same atoms READ on
  :math:`S^2/O(2)_a` for a named axis :math:`a` (the orbit space of the
  axis's STABILISER — rotations about it and mirrors through it; named
  by that group since #432), via
  :meth:`~orpheus.numerics.measure.DiscreteMeasure.on_orbit_space`. This
  is what the SN slab solver consumes, through
  :meth:`~orpheus.numerics.quadrature.Quadrature.gauss_legendre`, and
  what the quadrature registry's ``GaussLegendre1D`` spec builds — so
  its space is ``L2[S^2/O2_x]``, not ``L2[[-1,1]]``, and it can no
  longer be confused with an 8-node SPATIAL rule on the same interval
  (`[M]` the two were ``==`` and hash-equal before this).

The group the polar marginal was quotiented **by** — the stabiliser
:math:`O(2)_a` of the axis, whose rotation half's fiber action is the
curvilinear :math:`\alpha` term — is thus carried by the orbit space
itself (:attr:`Quotient.by`), and the
geometry table's
:attr:`~orpheus.numerics.quadrature.registry.AngularSymmetry.continuous_isotropy`
names the same group; stage 0 of quadrature selection compares the two.
Until 2026-08-03 the spent group was carried in the invariance TAG; see
the note at :func:`gauss_legendre_on_mu`'s return statement for why that
had to move, :class:`~orpheus.numerics.symmetry.SO2` for why the
group now names its axis, and :class:`~orpheus.numerics.symmetry.O2` for
why the entry is named by the stabiliser rather than the rotation half.

There is **no adapter class** in between — ``Quadrature`` is one type
with named classmethod factories.

References
----------

* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. Theorem 3.6.20 (Gauss-Legendre rule
  is exact for polynomials of degree :math:`\le 2N - 1`).
"""

from __future__ import annotations

from ..measure import (
    DiscreteMeasure,
)
from ..generating_measure import LEGENDRE
from ..manifold import SPHERE
from ..symmetry import SubgroupOfO3


def gauss_legendre_on_mu(n: int) -> DiscreteMeasure:
    r"""N-point Gauss-Legendre rule on :math:`\mu \in [-1, 1]`.

    The nodes :math:`\mu_i` are the roots of the Legendre polynomial
    :math:`P_n` and the weights are the Christoffel numbers. Built by
    :meth:`~orpheus.numerics.generating_measure.GeneratingMeasure.gauss`
    on :data:`~orpheus.numerics.generating_measure.LEGENDRE` — the same
    generic construction as every other Gauss family in the tree.
    Agrees with :func:`numpy.polynomial.legendre.leggauss` to 1-4 ULP,
    not bit-for-bit: numpy Newton-refines. Neither is "the" answer,
    since the exact nodes are irrational.

    Polynomial exactness (Stoer & Bulirsch 2002, Theorem 3.6.20):
    integrates every polynomial of degree :math:`\le 2n - 1` exactly,
    so the returned measure carries
    ``degree_of_exactness = 2*n - 1``.

    Weight sum: :math:`\sum_i w_i = 2`, matching the Lebesgue measure
    of :math:`[-1, 1]` (so the rule is unweighted in the classical
    sense).

    Symmetry: the tag is :math:`\sigma_x`, the reflection
    :math:`\mu \to -\mu` — the property SN consumers actually use, and
    the one this measure actually has. It holds **bit-exactly**: the
    generic construction imposes it rather than inheriting it (`[M]`
    defect 0.0 at every :math:`n`, where the previous route left
    ~1e-16). So the slab's two sweep senses pair by integer index —
    ``partner(i) = n - 1 - i`` is exact arithmetic, not a nearest-node
    search.

    Parameters
    ----------
    n : int
        Number of quadrature nodes; must be :math:`\ge 1`. SN slab
        sweeps require ``n`` even (one ordinate per direction sense),
        but this primitive accepts any :math:`n \ge 1`.

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(n,)``, weights shape ``(n,)``, on
        ``support=COSINE_INTERVAL`` — the chart, naming no axis — with
        ``invariance_group=SubgroupOfO3.Mirror("x")`` and an exactness
        claim of degree :math:`2n-1` against ``LEGENDRE``.

    See Also
    --------
    :func:`gauss_legendre_on_polar_orbit` — the same atoms declared on
    :math:`S^2/O(2)_a`; what the slab consumes.
    :meth:`orpheus.numerics.quadrature.Quadrature.gauss_legendre` — the
    named factory SN calls. It is a factory on the single
    ``Quadrature`` type, not an adapter class: the four per-family
    wrappers that once cached ``mu_x`` / ``mu_y`` / ``weights`` views
    were retired into classmethod factories.
    """
    # ONE construction. This is the Gauss rule for the Legendre measure,
    # so it is built by asking that measure for it — the same generic
    # Golub-Welsch body every other family goes through.
    #
    # It used to call `numpy.leggauss` instead. Both are correct; the two
    # differ by 1-4 ULP in the nodes because numpy Newton-refines, and
    # neither is "the" answer since the exact nodes are irrational. The
    # tie-breaker is that the generic body imposes the reflection
    # symmetry (see `GeneratingMeasure.gauss`), so `mu -> -mu` maps the
    # node set onto itself BIT-EXACTLY and the slab's two sweep senses
    # pair by index rather than by tolerance. Consolidating cost a
    # deliberate re-baselining of the SN slab snapshots, taken
    # 2026-08-02 (user ruling), and bought a single source of truth for
    # every Gauss rule in the tree.
    return LEGENDRE.gauss(n).with_metadata(
        # σ_x: μ -> -μ, and the generic Golub-Welsch body imposes it
        # BIT-EXACTLY (see the note above), which is what lets the slab's
        # two sweep senses pair by index rather than by tolerance.
        #
        # This field used to read `SubgroupOfO3.SO2`, defended as "the tag
        # names the group that was INTEGRATED OUT to produce the μ-marginal
        # — the spent half of the geometry's symmetry". That is a true and
        # useful fact, but it is a DIFFERENT fact from the one this field's
        # name promises, and a field whose values correlate with something
        # other than its name is doing two jobs. Its own defence gave the
        # tell: "SO(2) acts trivially on μ, so EVERY measure on [-1,1]
        # satisfies it" — a tag satisfied by every possible value carries no
        # information and cannot be wrong, so it could never have been
        # checked. `[M]` It also made the μ-marginal answer O(3) to the
        # invariance walk, and let an ASYMMETRIC μ-set certify as
        # SO(2)-invariant.
        #
        # The spent half is NOT lost: `gauss_legendre_on_polar_orbit`
        # declares the orbit space S²/O(2)_a this rule is read on, and
        # `AngularSymmetry.continuous_isotropy` (registry.py) names the
        # same group for the geometry — stage 0 compares the two.
        invariance_group=SubgroupOfO3.Mirror("x"),
    )


def gauss_legendre_on_polar_orbit(n: int, axis: str) -> DiscreteMeasure:
    r"""The :math:`n`-point Gauss-Legendre rule READ on :math:`S^2/O(2)_a`.

    The polar marginal about ``axis``: the same nodes and weights as
    :func:`gauss_legendre_on_mu`, whose :math:`\mu` is now known to be
    the cosine against a named axis — :math:`\mu = \hat\Omega \cdot
    \hat e_a`. The slab and sphere geometries take ``axis="x"``: the
    slab's streaming axis is the spatial :math:`x`, the 1-D embedding is
    :math:`(\mu, 0, 0)`, and `[M]` the real spherical-harmonic basis
    takes its pole there (``cos θ = μ_x``), so :math:`P_\ell(\mu)` are
    its :math:`m = 0` members.

    What the declaration buys, and it is a repair rather than wiring:
    the measure's space is ``L2[S^2/O2_x]``, so it no longer compares
    equal to a spatial rule on :math:`[-1, 1]` with the same node count
    (`[M]` 2026-09-01: an 8-node slab angular space and an 8-node spatial
    space were ``==`` AND hash-equal, 2.1's energy/spatial collision one
    level up), and the registry's stage-0 admission
    (:meth:`~orpheus.numerics.quadrature.registry.AngularSymmetry.admits_domain`)
    compares an orbit space that names its group against a geometry that
    names the group it spends — one fact on both sides, so a chart-level
    rule, or a rule about the WRONG axis, is refused there.

    The residual :math:`\{e, \sigma_a\}` — the reflection :math:`\mu \to
    -\mu` the slab owes its two sweep senses — descends to the orbit
    space (it normalises :math:`O(2)_a`) and the Gauss nodes are closed
    under it bit-exactly (see :func:`gauss_legendre_on_mu`), so it is
    re-tagged here for the named axis; :meth:`DiscreteMeasure.on_orbit_space`
    drops the chart's tag because a subgroup of :math:`O(3)` is a
    statement about an embedding.

    Parameters
    ----------
    n : int
        Number of nodes, as for :func:`gauss_legendre_on_mu`.
    axis : str
        The axis the marginal's cosine is measured against — the axis its
        stabiliser :math:`O(2)_a` fixes — one of ``"x"`` / ``"y"`` / ``"z"``.

    Returns
    -------
    DiscreteMeasure
        On ``support=SPHERE.quotient(SubgroupOfO3.O2(axis))``, with
        ``invariance_group=SubgroupOfO3.Mirror(axis)``; ``phase`` is
        ``"angular"`` from the manifold alone and ``quotient_group`` is
        the stabiliser :math:`O(2)_a` (#432: the orbit space is named by
        the largest group with its orbits, not by its rotation half).
    """
    orbit_space = SPHERE.quotient(SubgroupOfO3.O2(axis))
    return gauss_legendre_on_mu(n).on_orbit_space(orbit_space).with_metadata(
        invariance_group=SubgroupOfO3.Mirror(axis),
    )


__all__ = ["gauss_legendre_on_mu", "gauss_legendre_on_polar_orbit"]
