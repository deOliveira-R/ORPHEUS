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

Symmetry: :math:`D_{n_\phi h}`, of order :math:`4 n_\phi`. The
:math:`n_\phi` equispaced azimuths are preserved by the cyclic
rotations :math:`C_{n_\phi}` and by the vertical mirrors at
:math:`k\pi/n_\phi`; the Gauss-Legendre :math:`\mu` nodes are
symmetric about :math:`\mu = 0`, which supplies :math:`\sigma_h`.
Together these generate :math:`D_{n_\phi h}` — and
:func:`~orpheus.numerics.symmetry.maximal_invariance_groups`
**computes** exactly that from the nodes, so the tag above is pinned
by a permanent test rather than merely asserted. It is not verified
at construction: the computed check runs in ``tests/``, not here.

.. caution::

   Both sides of that check are ``cos``/``sin`` evaluations — the
   :math:`\phi` grid here, and the checker's own :math:`C_n` and
   :math:`\sigma_v` operators — so the agreement holds to ~1e-16, not
   bit-exactly, and the match window
   (``symmetry._NODE_WINDOW_FACTOR``) is what absorbs the difference.
   Lebedev's :math:`O_h` claim, by contrast, is exact on both sides
   because signed permutations are exact in IEEE. Making this one
   exact needs :mod:`orpheus.numerics.roots_of_unity` on the generator
   AND on the checker's cyclic/mirror operators (issue #325).

.. note::

   Until 2026-08-02 this rule tagged :math:`SO(2)` and this
   paragraph defended it as "a conservative upper bound". Both
   halves were wrong, in opposite directions. **No finite point set
   on** :math:`S^2` **is** :math:`SO(2)`-**closed**, so the tag was
   simply false — and for an *invariance* claim a larger group is a
   **stronger** claim, so an upper bound is an over-claim, never a
   conservative one. (Contrast ``degree_of_exactness``, where the
   ``min`` genuinely *is* conservative: a lower degree is a weaker
   claim.) The prose also undersold the truth in the other
   direction, calling the rule ":math:`C_{n_\phi}`-invariant
   strictly" — :math:`C_{n_\phi}` is an index-4 subgroup of the real
   answer. The false tag was load-bearing: it was the only reason
   the registry's geometry gate admitted this rule for a cylinder
   (see :func:`~orpheus.numerics.quadrature.registry.select_quadrature`).

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
from .rules_1d import gauss_legendre_on_mu
from .rules_sphere import LevelStructure, PolarInvariant


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
        ``invariance_group=SubgroupOfO3.Dnh(n_phi)``,
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

    # GL points in μ = cos(θ). The polar factor IS the slab rule — same
    # construction, one source — rather than a third independent call to
    # a Gauss-Legendre routine. Consolidated 2026-08-02; it had been
    # `np.polynomial.legendre.leggauss(n_mu)`, which the audit found as
    # a surviving twin of the two rules the campaign was retiring.
    polar = gauss_legendre_on_mu(n_mu)
    mu_gl, w_gl = polar.nodes, polar.weights

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
        # D_{n_φ h}, NOT SO(2). The φ grid is a FINITE set of n_φ
        # equispaced azimuths, so no continuous rotation preserves it —
        # only the n_φ-fold cyclic rotations, the vertical mirrors at
        # kπ/n_φ, and (because the GL μ nodes are symmetric) σ_h. The
        # rule carried ``SO2`` until 2026-08-02: a claim no finite point
        # set on S² can satisfy. See ``maximal_invariance_groups``, which
        # computes this group from the nodes and pins the declaration.
        invariance_group=SubgroupOfO3.Dnh(n_phi),
        degree_of_exactness=min(2 * n_mu - 1, n_phi - 1),
    )
    structure = LevelStructure(
        n_levels=n_mu,
        level_indices=level_indices,
        level_mu=mu_gl,
        # SIGNED: each GL polar node is its own level, so the levels run
        # over all of [-1, 1] and a level's fiber is ONE circle. The
        # level-symmetric rule fibers over |mu_z| instead, and its
        # levels carry two.
        polar_invariant=PolarInvariant.SIGNED_MU_Z,
        # The fiber coordinate, taken from the phi GRID rather than
        # recovered from the cosines: phi_pts[m] is exactly the angle
        # that generated ordinate m.
        azimuth=np.tile(phi_pts, n_mu),
        hemisphere=np.sign(mu_z).astype(np.int64),
    )
    return measure, structure


__all__ = ["product_mu_phi"]
