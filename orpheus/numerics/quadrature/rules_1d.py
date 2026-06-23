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

import numpy as np

from ..measure import (
    SPACE_INTERVAL_M11,
    DiscreteMeasure,
)
from ..symmetry import SubgroupOfO3


def gauss_legendre_on_mu(n: int) -> DiscreteMeasure:
    r"""N-point Gauss-Legendre rule on :math:`\mu \in [-1, 1]`.

    The nodes :math:`\mu_i` are the roots of the Legendre polynomial
    :math:`P_n` and the weights are the Christoffel numbers — exactly
    what :func:`numpy.polynomial.legendre.leggauss` returns.

    Polynomial exactness (Stoer & Bulirsch 2002, Theorem 3.6.20):
    integrates every polynomial of degree :math:`\le 2n - 1` exactly,
    so the returned measure carries
    ``degree_of_exactness = 2*n - 1``.

    Weight sum: :math:`\sum_i w_i = 2`, matching the Lebesgue measure
    of :math:`[-1, 1]` (so the rule is unweighted in the classical
    sense).

    Symmetry: :math:`SO(2)` — the rule is invariant under rotations
    about the polar axis (these act trivially on
    :math:`\mu`). Stronger symmetry would require nodes paired
    :math:`(\mu_i, -\mu_i)` symmetrically; this happens to hold for
    Gauss-Legendre but is not the maximal property we tag here, since
    SN consumers care about :math:`SO(2)` for slab geometry.

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
    if n < 1:
        raise ValueError(f"gauss_legendre_on_mu requires n >= 1, got n={n}")
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=SPACE_INTERVAL_M11,
        invariance_group=SubgroupOfO3.SO2,
        degree_of_exactness=2 * n - 1,
    )


__all__ = ["gauss_legendre_on_mu"]
