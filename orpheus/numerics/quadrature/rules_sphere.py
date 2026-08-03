r"""Sphere quadrature rules on :math:`S^2`.

Two families:

* **Lebedev rules** — :math:`O_h`-invariant grids on the unit sphere
  with high polynomial exactness for the order chosen (Lebedev 1976).
  Returned by :func:`lebedev_sphere`.
* **Level-symmetric :math:`S_N` rules** — Carlson-Lathrop's
  :math:`O_h`-invariant triangular discrete-ordinate grids organised
  into :math:`N/2` polar levels per hemisphere with permutations
  :math:`(\eta, \xi, \mu)` per octant (Carlson & Lathrop 1968, Lewis &
  Miller §4.2). Returned by :func:`level_symmetric_sn`.

Both families are :math:`O_h`-invariant, so the returned measures
carry ``invariance_group=SubgroupOfO3.OctahedralOh``.

The level-symmetric rule additionally carries per-level metadata
inside the ``levels`` field of the result (a small adjacent
namedtuple), needed by the cylindrical SN sweep for azimuthal
redistribution. The :math:`O(1)`-attribute SN-side adapter
(:class:`~orpheus.sn.quadrature.LevelSymmetricSN`) caches that
metadata directly.

References
----------

* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Comp. Math.*
  **16**(2), 10-24.
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan et al., eds., Gordon & Breach.
* Lewis, E.E. and Miller, W.F. (1993). *Computational Methods of
  Neutron Transport*. Wiley. §4.2 (level-symmetric construction
  and degree of exactness).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from ..exactness import UNIFORM_ON_SPHERE, ExactnessClaim
from ..measure import SPACE_SPHERE, DiscreteMeasure
from ..symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Lebedev sphere quadrature
# ---------------------------------------------------------------------------


def lebedev_sphere(order: int) -> DiscreteMeasure:
    r"""Lebedev :math:`O_h`-invariant rule on :math:`S^2` of given order.

    Wraps :func:`scipy.integrate.lebedev_rule`. The ``order`` parameter
    is the maximum degree of the polynomial integrated exactly, so the
    returned measure carries ``degree_of_exactness=order``. The number
    of nodes is determined by SciPy's tabulated Lebedev grid for that
    order — typical sizes are 6, 14, 26, 38, 50, …, 590.

    Weight sum: :math:`\sum_i w_i = 4\pi`, matching the area of
    :math:`S^2`.

    Symmetry: :math:`O_h` (full octahedral group, 48 elements). The
    construction is :math:`O_h`-invariant by design (Lebedev 1976) —
    every node sits in an :math:`O_h`-orbit and the corresponding
    weights are equal across the orbit.

    Nodes are shape ``(N, 3)`` carrying the Cartesian direction
    cosines :math:`(\mu_x, \mu_y, \mu_z)` per row.

    Parameters
    ----------
    order : int
        Lebedev order; must match a tabulated value (SciPy raises
        otherwise).

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(N, 3)``, weights shape ``(N,)``, on
        ``space="S^2"``, with ``invariance_group=SubgroupOfO3.OctahedralOh`` and
        ``degree_of_exactness=order``.

    See Also
    --------
    :class:`orpheus.sn.quadrature.LebedevSphere` — the SN-side adapter
    that caches ``mu_x`` / ``mu_y`` / ``mu_z`` views and reflection
    indices.
    """
    from scipy.integrate import lebedev_rule

    pts, w = lebedev_rule(order)
    # ``pts`` from scipy has shape ``(3, N)``; transpose to ``(N, 3)``
    # for the canonical ``DiscreteMeasure`` (N, d) layout.
    nodes = np.ascontiguousarray(pts.T)  # (N, 3)
    return DiscreteMeasure(
        nodes=nodes,
        weights=w,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        # SPHERICAL-HARMONIC degree, against Lebesgue measure on S^2 —
        # not an algebraic degree, and not against a weight on an
        # interval. Lebedev's reference generates nothing by
        # Golub-Welsch: the claim's authority is the published table
        # (Lebedev 1976), which is exactly why the reference here is a
        # ``UniformMeasure`` rather than a ``GeneratingMeasure``.
        exactness=ExactnessClaim(reference=UNIFORM_ON_SPHERE, degree=order),
    )


# ---------------------------------------------------------------------------
# Level-symmetric S_N quadrature
# ---------------------------------------------------------------------------


class PolarInvariant(Enum):
    r"""Which invariant a quadrature's levels are the fibers OF.

    A level is a **fiber of an invariant**, not an orbit, and the two
    shipped producers fiber over *different* invariants — which the
    single ``LevelStructure`` type used to erase.

    * :attr:`SIGNED_MU_Z` — :func:`product_mu_phi`. Each Gauss-Legendre
      polar node is its own level, so the levels run over the whole
      range :math:`\mu_z \in [-1, 1]` and the fiber over a level is a
      full circle.
    * :attr:`ABS_MU_Z` — :func:`level_symmetric_sn`. Levels are indexed
      by :math:`|\mu_z|`, so each carries BOTH hemispheres and the
      fiber is two circles.

    Consequence, and the reason this has to be on the type: "the
    ordinates of level :math:`p`" means a different set in the two
    cases, so any consumer that reads ``level_indices`` without knowing
    the invariant is reading an ambiguous object.
    """

    SIGNED_MU_Z = "signed_mu_z"
    ABS_MU_Z = "abs_mu_z"


@dataclass(frozen=True)
class LevelStructure:
    r"""Per-level metadata for a triangular SN-style sphere quadrature.

    Captured alongside the :class:`DiscreteMeasure` returned by
    :func:`level_symmetric_sn` and :func:`product_mu_phi`. The cylindrical
    SN sweep needs this structure to compute the azimuthal redistribution
    coefficients (Bailey et al. 2009, Eq. 50).

    Attributes
    ----------
    n_levels : int
        Number of polar levels per hemisphere (or polar axis, for the
        product-quadrature variant).
    level_indices : list[np.ndarray]
        For each level :math:`p`, the indices into the flattened node
        array of ordinates on that level, sorted by increasing
        :math:`\eta = \mu_x` (the radial cosine — the cylindrical
        sweep convention from Bailey et al. 2009 Eq. 50).

        .. warning::

           That sort key is **2-to-1 on the fiber**. A level is a
           circle, and :math:`\eta = \sin\theta\cos\varphi` is even in
           :math:`\varphi`, so ordering by it is an ordering of the
           circle *modulo the mirror* rather than of the circle.
           `[M]` On one level of ``product(2, 8)``, 8 ordinates give
           only 5 distinct :math:`\eta`, and the resulting order does
           not determine :math:`\operatorname{sign}(\xi)` — that
           information is not in the key. Use :attr:`azimuth` when the
           ordering must be a function of the fiber. This is the
           mechanism behind issue #326.
    level_mu : np.ndarray
        The invariant's value per level — read
        :attr:`polar_invariant` to know *which* invariant.
    polar_invariant : PolarInvariant
        Which invariant the levels fiber over. Not decoration: it
        decides whether level :math:`p` holds one hemisphere's circle
        or both.
    azimuth : np.ndarray
        Per-ordinate azimuthal angle :math:`\varphi \in [0, 2\pi)`.
    hemisphere : np.ndarray
        Per-ordinate :math:`\operatorname{sign}(\mu_z)`, as ``-1`` /
        ``0`` / ``+1``.

        Together with :attr:`azimuth` this is **the fiber's own
        coordinate**, and it takes both because the fiber is not always
        one circle: under :attr:`PolarInvariant.ABS_MU_Z` a level
        carries BOTH hemispheres, so :math:`\varphi` alone repeats. The
        pair is injective on a level under either invariant, which is
        exactly what an ordering of the fiber needs.
    """

    n_levels: int
    level_indices: list[np.ndarray]
    level_mu: np.ndarray
    polar_invariant: PolarInvariant
    azimuth: np.ndarray
    hemisphere: np.ndarray

    def fiber(self, level: int) -> np.ndarray:
        """Ordinate indices of ``level``, ordered by the FIBER coordinate.

        Lexicographic in ``(hemisphere, azimuth)`` — an ordering that is
        a function of the fiber, unlike :attr:`level_indices`, which is
        ordered by a projection that is 2-to-1 on it.

        Kept as a separate accessor rather than replacing the stored
        order: the stored order is what the cylindrical sweep consumes
        today, and changing it moves results (issue #326). This gives
        the correct ordering a name and a home first.
        """
        members = self.level_indices[level]
        keys = np.lexsort((self.azimuth[members], self.hemisphere[members]))
        return members[keys]


def _build_level_symmetric_arrays(
    sn_order: int,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, np.ndarray,
    list[np.ndarray], np.ndarray, np.ndarray,
]:
    """Construct level-symmetric :math:`S_N` quadrature arrays.

    This implementation mirrors the existing
    :func:`orpheus.sn.quadrature._build_level_symmetric` byte-for-byte
    (verified by the bit-identical foundation tests at
    ``tests/numerics/test_rules_sphere.py``).

    Returns
    -------
    mu_x, mu_y, mu_z, weights : np.ndarray
        Direction cosines and weights, shape ``(N,)`` each.
    n_levels : int
        Number of polar levels per hemisphere.
    level_mu : np.ndarray
        Polar cosine per level.
    level_indices : list[np.ndarray]
        Per-level ordinate indices, sorted by increasing ``mu_x``.
    azimuth, hemisphere : np.ndarray
        The fiber coordinate (Q4): :math:`\\varphi \\in [0, 2\\pi)` and
        :math:`\\operatorname{sign}(\\mu_z)` per ordinate. The annotation
        omitted these two from 2026-08-02 (`3afb52c2`) until the exactness
        carve, while the body returned all nine — a stale signature the
        Sphinx build could not see and the tests could not fail on.
    """
    if sn_order % 2 != 0 or sn_order < 2:
        raise ValueError(f"S_N order must be positive even, got {sn_order}")

    n_half = sn_order // 2

    if n_half == 1:
        mu2_levels = np.array([1.0 / 3.0])
    else:
        mu1_sq = 1.0 / (sn_order * (sn_order + 2) / 4)
        delta = 2.0 * (1.0 - 3.0 * mu1_sq) / (sn_order - 2)
        mu2_levels = mu1_sq + np.arange(n_half) * delta

    mu_levels = np.sqrt(mu2_levels)

    # The defining property of a level-symmetric set: every ordinate's
    # three direction cosines are drawn from the SAME level array, so
    # the set is closed under permuting the axes. The level index of the
    # third cosine is not free — it is fixed by the other two.
    #
    #   mu2[p] + mu2[k] + mu2[j] = 3*mu1_sq + (p + k + j)*delta = 1
    #   and delta = 2(1 - 3*mu1_sq)/(N - 2)
    #   =>  p + k + j = (1 - 3*mu1_sq)/delta = (N - 2)/2 = n_half - 1.
    #
    # So j is INDEX ARITHMETIC. Until 2026-08-02 this loop instead
    # recovered the third cosine numerically as
    # ``sqrt(1 - mu_z**2 - eta**2)``, which lands within ~1e-16 of
    # ``mu_levels[j]`` but not ON it: `[M]` at N=16 the y-axis then
    # carried 22 distinct magnitudes where the level array has 8, and
    # only 8 of the 48 O_h operators mapped the node set onto itself
    # bit-exactly — the 8 pure sign flips, since negation is exact in
    # IEEE while a coordinate swap is only exact if the values match.
    # The rule advertised O_h invariance and realized D_2h.
    #
    # Reading the level value instead makes all 48 exact, so the
    # octahedral symmetry becomes an integer permutation of ordinate
    # indices rather than a question about tolerances. It also retires
    # the `xi_sq < -1e-14` guard: `j` is provably in range for every
    # admissible (p, k), because p + k <= n_half - 1 by the loop bound.
    octant_dirs: list[tuple[float, float, float]] = []
    for p in range(n_half):
        mu_z = mu_levels[p]
        for k in range(n_half - p):
            j = n_half - 1 - p - k
            octant_dirs.append((mu_levels[k], mu_levels[j], mu_z))

    n_octant = len(octant_dirs)
    w_octant = 4.0 * np.pi / (8.0 * n_octant)

    all_eta: list[float] = []
    all_xi: list[float] = []
    all_mu: list[float] = []
    all_w: list[float] = []
    for eta, xi, mu_z in octant_dirs:
        for s_eta in (-1, 1):
            for s_xi in (-1, 1):
                for s_mu in (-1, 1):
                    all_eta.append(s_eta * eta)
                    all_xi.append(s_xi * xi)
                    all_mu.append(s_mu * mu_z)
                    all_w.append(w_octant)

    mu_x = np.array(all_eta)
    mu_y = np.array(all_xi)
    mu_z_arr = np.array(all_mu)
    weights = np.array(all_w)

    n_levels = n_half
    level_mu_vals = mu_levels
    level_indices: list[np.ndarray] = []
    for p in range(n_levels):
        # EXACT membership. The 8-fold sign replication above copies
        # mu_z straight out of `mu_levels`, so |mu_z| IS the level value
        # bit-for-bit and the comparison is an equality, not a
        # neighbourhood test. This carried `tol = 1e-12` until
        # 2026-08-02 — a symmetry question answered by float comparison
        # in a loop that already knew the answer exactly.
        idx = np.where(np.abs(mu_z_arr) == level_mu_vals[p])[0]
        order = np.argsort(mu_x[idx])
        level_indices.append(idx[order])

    # The fiber coordinate: azimuth about the polar axis. Distinct
    # ordinates on a level sit at distinct phi, which is what makes an
    # ordering by it a function OF the fiber (unlike mu_x, which is even
    # in phi and therefore 2-to-1 on it).
    azimuth = np.mod(np.arctan2(mu_y, mu_x), 2.0 * np.pi)
    hemisphere = np.sign(mu_z_arr).astype(np.int64)

    return (mu_x, mu_y, mu_z_arr, weights, n_levels, level_mu_vals,
            level_indices, azimuth, hemisphere)


def level_symmetric_sn(
    sn_order: int,
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Carlson-Lathrop level-symmetric :math:`S_N` rule on :math:`S^2`.

    Standard triangular discrete-ordinate quadrature with :math:`N/2`
    polar levels per hemisphere. The construction is
    :math:`O_h`-invariant: every octant carries the same set of
    direction-cosine triples up to sign permutations, with equal
    weights inside an octant.

    Weight sum: :math:`\sum_i w_i = 4\pi`.

    Symmetry: :math:`O_h` — invariant under all 48 rotation /
    reflection elements of the octahedral group.

    Polynomial exactness: depends on :math:`N`. For the simple
    equal-weight construction implemented here (matching the existing
    :func:`orpheus.sn.quadrature._build_level_symmetric`), the rule is
    exact at degree :math:`1` for :math:`N = 2` (zeroth and first
    moments) and reaches degree :math:`N - 1` for higher orders under
    the moment-conditions construction in Carlson & Lathrop 1968.
    Conservatively, ``degree_of_exactness=N-1`` is recorded; consumers
    that need a tighter guarantee should refer to Lewis & Miller
    Table 4-2.

    Returns the :class:`DiscreteMeasure` and a
    :class:`LevelStructure` capturing the per-level grouping needed by
    the cylindrical SN sweep.

    Parameters
    ----------
    sn_order : int
        Even :math:`N \ge 2`. Common values: 4, 8, 16.

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(N_total, 3)``, weights shape ``(N_total,)``.
        ``invariance_group=SubgroupOfO3.OctahedralOh``,
        ``degree_of_exactness=sn_order-1``.
    LevelStructure
        Per-level indexing metadata.

    See Also
    --------
    :class:`orpheus.sn.quadrature.LevelSymmetricSN` — the SN-side
    adapter that caches all of this for hot-path attribute access.
    """
    mu_x, mu_y, mu_z, w, n_levels, level_mu, level_indices, azimuth, hemisphere = (
        _build_level_symmetric_arrays(sn_order)
    )
    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=w,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        # ⛔ This claim is FALSE and known to be — issue #327. `[M]` the
        # measured degree is 3 at EVERY order, an over-claim of 12 at
        # S_16. It is recorded unchanged here so the carve stays
        # behaviour-neutral; #327 is where it gets corrected.
        #
        # Naming the reference makes the defect legible rather than
        # fixing it: this rule has NO generating measure (it hand-assigns
        # one weight to every ordinate), so nothing constrains the
        # integer, and the claim rests on an authority the construction
        # does not actually implement. That is the mechanism of #327.
        exactness=ExactnessClaim(
            reference=UNIFORM_ON_SPHERE, degree=sn_order - 1,
        ),
    )
    structure = LevelStructure(
        n_levels=n_levels,
        level_indices=level_indices,
        level_mu=level_mu,
        polar_invariant=PolarInvariant.ABS_MU_Z,
        azimuth=azimuth,
        hemisphere=hemisphere,
    )
    return measure, structure


__all__ = [
    "LevelStructure",
    "lebedev_sphere",
    "level_symmetric_sn",
]
