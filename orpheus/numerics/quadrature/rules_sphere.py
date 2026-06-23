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

import numpy as np

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
        space=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        degree_of_exactness=order,
    )


# ---------------------------------------------------------------------------
# Level-symmetric S_N quadrature
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LevelStructure:
    """Per-level metadata for a triangular SN-style sphere quadrature.

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
        :math:`\\eta = \\mu_x` (the radial cosine — the cylindrical
        sweep convention from Bailey et al. 2009 Eq. 50).
    level_mu : np.ndarray
        Polar cosine value :math:`\\mu_p \\ge 0` per level.
    """

    n_levels: int
    level_indices: list[np.ndarray]
    level_mu: np.ndarray


def _build_level_symmetric_arrays(
    sn_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, np.ndarray, list[np.ndarray]]:
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

    octant_dirs: list[tuple[float, float, float]] = []
    for p in range(n_half):
        mu_z = mu_levels[p]
        sin_theta_sq = 1.0 - mu_z**2
        n_azi = n_half - p
        for k in range(n_azi):
            eta = mu_levels[k]
            xi_sq = sin_theta_sq - eta**2
            if xi_sq < -1e-14:
                continue
            xi = np.sqrt(max(xi_sq, 0.0))
            octant_dirs.append((eta, xi, mu_z))

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
        tol = 1e-12
        idx = np.where(np.abs(np.abs(mu_z_arr) - level_mu_vals[p]) < tol)[0]
        order = np.argsort(mu_x[idx])
        level_indices.append(idx[order])

    return mu_x, mu_y, mu_z_arr, weights, n_levels, level_mu_vals, level_indices


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
    mu_x, mu_y, mu_z, w, n_levels, level_mu, level_indices = (
        _build_level_symmetric_arrays(sn_order)
    )
    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=w,
        space=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        degree_of_exactness=sn_order - 1,
    )
    structure = LevelStructure(
        n_levels=n_levels,
        level_indices=level_indices,
        level_mu=level_mu,
    )
    return measure, structure


__all__ = [
    "LevelStructure",
    "lebedev_sphere",
    "level_symmetric_sn",
]
