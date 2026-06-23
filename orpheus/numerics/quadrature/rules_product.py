r"""Product quadrature on the sphere: Gauss-Legendre :math:`(\mu)`
:math:`\times` equispaced :math:`(\phi)`.

The polar angle is discretised by Gauss-Legendre on
:math:`\mu = \cos\theta \in [-1, 1]` and the azimuthal angle is
sampled at the :math:`n_\phi` left-endpoints of an even partition of
:math:`[0, 2\pi)` (matching the long-standing convention in
:class:`orpheus.sn.quadrature.ProductQuadrature` — the bit-identical
contract enforced by the regression snapshots requires this exact
convention).

Direction cosines per ordinate:

.. math::
   :label: product-mu-phi-cosines

   \mu_x &= \sin\theta \cos\phi, \\
   \mu_y &= \sin\theta \sin\phi, \\
   \mu_z &= \cos\theta = \mu.

Weight per ordinate: :math:`w = w_{\text{GL}}(\mu) \cdot (2\pi /
n_\phi)`. Total weight sum: :math:`4\pi`.

Symmetry: :math:`SO(2)` about the :math:`\mu_z` axis. The polar
factor is :math:`SO(2)`-invariant trivially (it does not see
:math:`\phi`); the azimuthal factor's :math:`n_\phi`-fold cyclic
symmetry is a discrete subgroup of :math:`SO(2)`. The maximal
continuous-rotation symmetry the discrete rule respects is
:math:`SO(2)` only when :math:`n_\phi \to \infty`; for finite
:math:`n_\phi` the rule is :math:`C_{n_\phi}`-invariant strictly.
For the algebra-of-record, we tag :math:`SO(2)` as a conservative
upper bound (the discrete rule is invariant under the continuous
group only in a moment sense; a stricter consumer can use
:meth:`SubgroupOfO3.is_invariant` to verify which exact group the
discrete grid satisfies).

References
----------

* Bailey, T.S., Adams, M.L., Yang, B., and Zika, M.R. (2009). "A
  piecewise linear discontinuous finite element spatial discretization
  of the transport equation." *Annals of Nuclear Energy* **35**,
  1929-1936. Eq. 50 — α-recursion convention used by the cylindrical
  sweep that consumes this product rule.
"""

from __future__ import annotations

import numpy as np

from ..measure import SPACE_SPHERE, DiscreteMeasure
from ..symmetry import SubgroupOfO3
from .rules_sphere import LevelStructure


def product_mu_phi(
    n_mu: int, n_phi: int
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Build the GL :math:`\times` equispaced product rule on :math:`S^2`.

    Implements the layout in
    :eq:`product-mu-phi-cosines`. Ordinate ordering matches the
    long-standing :class:`orpheus.sn.quadrature.ProductQuadrature`
    convention bit-for-bit: outer loop over the :math:`n_\mu`
    Gauss-Legendre nodes (axial cosine), inner loop over the
    :math:`n_\phi` azimuthal samples. Per-level indexing lists are
    sorted by increasing :math:`\eta = \mu_x` to match the
    cylindrical-sweep convention from Bailey et al. (2009) Eq. 50.

    Weight sum: :math:`\sum_i w_i = 4\pi`.

    Polynomial exactness: the polar factor is exact at
    degree :math:`2 n_\mu - 1` (Stoer-Bulirsch); the azimuthal
    midpoint rule is exact for trigonometric polynomials of degree
    :math:`< n_\phi`. We tag the conservative
    ``degree_of_exactness = min(2*n_mu - 1, n_phi - 1)`` as the
    weakest of the two — for general polynomial integrands on
    :math:`S^2`, both factors must be exact simultaneously.

    Parameters
    ----------
    n_mu : int
        Number of Gauss-Legendre points in :math:`\mu`. Must be
        :math:`\ge 1`.
    n_phi : int
        Number of equispaced points in :math:`\phi`. Must be
        :math:`\ge 1`.

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(n_mu * n_phi, 3)``, weights shape
        ``(n_mu * n_phi,)``, on ``space="S^2"``.
        ``invariance_group=SubgroupOfO3.SO2``,
        ``degree_of_exactness=min(2*n_mu-1, n_phi-1)``.
    LevelStructure
        Per-:math:`\mu`-level indexing metadata used by the
        cylindrical SN sweep.

    See Also
    --------
    :class:`orpheus.sn.quadrature.ProductQuadrature` — the SN-side
    adapter caching this rule's outputs for hot-path access.
    """
    if n_mu < 1:
        raise ValueError(f"product_mu_phi requires n_mu >= 1, got n_mu={n_mu}")
    if n_phi < 1:
        raise ValueError(f"product_mu_phi requires n_phi >= 1, got n_phi={n_phi}")

    # GL points in μ = cos(θ).
    mu_gl, w_gl = np.polynomial.legendre.leggauss(n_mu)

    # Equispaced φ in [0, 2π) — left-endpoints, matching the legacy
    # ProductQuadrature contract pinned by the regression snapshots.
    phi_pts = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    w_phi = 2.0 * np.pi / n_phi

    n_total = n_mu * n_phi
    mu_x = np.empty(n_total)
    mu_y = np.empty(n_total)
    mu_z = np.empty(n_total)
    weights = np.empty(n_total)
    level_indices: list[np.ndarray] = []

    idx = 0
    for p in range(n_mu):
        mu_val = mu_gl[p]
        sin_theta = np.sqrt(1.0 - mu_val**2)
        level_idx: list[int] = []
        for m in range(n_phi):
            mu_x[idx] = sin_theta * np.cos(phi_pts[m])
            mu_y[idx] = sin_theta * np.sin(phi_pts[m])
            mu_z[idx] = mu_val
            weights[idx] = w_gl[p] * w_phi
            level_idx.append(idx)
            idx += 1
        # Sort by increasing η (mu_x) for cylindrical sweep convention.
        level_arr = np.array(level_idx)
        order = np.argsort(mu_x[level_arr])
        level_indices.append(level_arr[order])

    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.SO2,
        degree_of_exactness=min(2 * n_mu - 1, n_phi - 1),
    )
    structure = LevelStructure(
        n_levels=n_mu,
        level_indices=level_indices,
        level_mu=mu_gl,
    )
    return measure, structure


__all__ = ["product_mu_phi"]
