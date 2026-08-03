r"""1-D quadrature rules on the polar cosine interval :math:`[-1, 1]`.

These rules return :class:`~orpheus.numerics.measure.DiscreteMeasure`
instances tagged with the :math:`\sigma_x` invariance group — the
reflection :math:`\mu \to -\mu`, which under the tree's canonical
:math:`(\mu, 0, 0)` embedding is the mirror with normal :math:`\hat e_x`.
It is the symmetry the slab and sphere geometries actually owe, and the
one their two sweep senses consume.

The group the polar marginal was quotiented **by** — :math:`SO(2)`, whose
fiber action is the curvilinear :math:`\alpha` term — is a different fact
and lives in
:attr:`~orpheus.numerics.quadrature.registry.AngularSymmetry.continuous_isotropy`,
where it derives the support :math:`S^2/SO(2) = [-1, 1]`. It was carried
in the invariance tag until 2026-08-03; see the note at
:func:`gauss_legendre_on_mu`'s return statement for why that had to move.

The Gauss-Legendre rule is the canonical 1-D primitive, consumed by the
SN slab solver through
:meth:`~orpheus.numerics.quadrature.Quadrature.gauss_legendre` and by the
product rule in :mod:`.rules_product`. There is **no adapter class** in
between — ``Quadrature`` is one type with named classmethod factories.

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
        ``space="[-1,1]"``, with
        ``invariance_group=SubgroupOfO3.Mirror("x")`` and
        ``degree_of_exactness=2n-1``.

    See Also
    --------
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
        # The spent half is NOT lost: it lives in
        # `AngularSymmetry.continuous_isotropy` (registry.py), which is
        # where it belongs and where it does real work — deriving the
        # support S²/SO(2) = [-1,1] that this measure carries.
        invariance_group=SubgroupOfO3.Mirror("x"),
    )


__all__ = ["gauss_legendre_on_mu"]
