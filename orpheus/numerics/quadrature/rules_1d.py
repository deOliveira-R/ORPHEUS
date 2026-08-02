r"""1-D quadrature rules on the polar cosine interval :math:`[-1, 1]`.

These rules return :class:`~orpheus.numerics.measure.DiscreteMeasure`
instances tagged with the :math:`SO(2)` invariance group — the maximal
:math:`O(3)` subgroup that fixes a 1-D measure on the polar cosine,
because rotations about the polar axis act trivially on
:math:`\mu = \cos\theta`.

The Gauss-Legendre rule is the canonical 1-D primitive consumed by
the SN slab solver via the
:class:`~orpheus.sn.quadrature.GaussLegendre1D` adapter (and by the
product rule in :mod:`.rules_product`).

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

    Symmetry: the tag is :math:`SO(2)`, the group the polar marginal
    was quotiented BY (see the note at the return statement). The
    property SN consumers actually use is the reflection
    :math:`\mu \to -\mu`, and that one holds **bit-exactly**: the
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
        ``invariance_group=SubgroupOfO3.SO2`` and
        ``degree_of_exactness=2n-1``.

    See Also
    --------
    :class:`orpheus.sn.quadrature.GaussLegendre1D` — the SN-side
    adapter that caches ``mu_x`` / ``mu_y`` / ``weights`` views of
    this measure.
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
        # SO(2) here is CORRECT, and correct for a different reason
        # than a rule on S² would be: this measure lives on the
        # QUOTIENT S²/SO(2) = [-1,1]. The tag names the group that was
        # integrated out to produce the μ-marginal — the "spent" half
        # of the geometry's symmetry, whose fiber action is the
        # curvilinear α term. Do NOT "correct" this to a finite group
        # by analogy with :func:`product_mu_phi` (whose SO(2) tag WAS
        # false): SO(2) acts trivially on μ, so every measure on
        # [-1,1] satisfies it, and here that triviality is the point.
        invariance_group=SubgroupOfO3.SO2,
    )


__all__ = ["gauss_legendre_on_mu"]
