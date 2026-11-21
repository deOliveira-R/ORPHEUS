"""Test-side oracles of the geometry layer: the four methods that left production in wave W0, unit ``dead``.

Rows 6, 32, 33 and 36 of the dead census had no production reader; each is a closed-form or enumeration oracle the geometry gates
read, and each is a free function of its former receiver (decision D2: the object is the first argument, a property becomes a
function of the object). The helpers a moved body calls (``_sine``, ``_realize``, ``permutation_between``) stay in production and are
imported from it.

* :func:`tangencies` (was ``DirectionDomain.tangencies``): the closed-form break set of a direction domain at a level.
* :func:`_group_elements` (was ``symmetry._group_elements``): every element of a finite tag's group, ``None`` for a continuous one.
* :func:`stabilizer_order` (was the property ``OrbitCertificate.stabilizer_order``): the orbit-stabilizer oracle.
* :func:`permutes` (was ``RigidMotion.permutes``): the permutation a motion induces on a point set, or ``None``.
"""

from __future__ import annotations

from math import pi
from typing import TYPE_CHECKING, assert_never

import numpy as np

from orpheus.geometry.chart import DirectionDomain, DirectionShape, _sine
from orpheus.geometry.symmetry import OrbitCertificate, _realize
from orpheus.geometry.transformation import Permutation, RigidMotion, permutation_between

if TYPE_CHECKING:
    from numpy.typing import ArrayLike


def tangencies(direction_domain: DirectionDomain, level: float) -> np.ndarray:
    r"""The values of the first axis where the line through ``direction_domain.point`` is tangent to :math:`c` = ``level``.

    Defined where the chart's group acts on the kept space, for
    :math:`0 < \ell \le c`, where :math:`b = \ell` has solutions:
    :math:`\Omega_x = \pm\sqrt{(1 - \ell/c)(1 + \ell/c)}` on the sphere,
    :math:`\alpha = \arcsin(\ell/c)` and :math:`\pi - \arcsin(\ell/c)` on
    the cylinder, one value (the grazing direction) at :math:`\ell = c`.
    Sorted; empty elsewhere, on the slab and on a stratum.
    """
    r, l = direction_domain.orbit_coordinate, float(level)
    if not direction_domain.chart.acts_on_kept_space or not 0.0 < l <= r:
        return np.empty(0)
    ratio = l / r
    match direction_domain.shape:
        case DirectionShape.COSINE:
            m = float(_sine(np.asarray(ratio)))
            return np.unique([-m, m])
        case DirectionShape.ANGLE_AXIAL:
            a = float(np.arcsin(ratio))
            return np.unique([a, pi - a])
        case DirectionShape.WHOLE | DirectionShape.AXIAL_COSINE:
            return np.empty(0)
        case unreachable:
            assert_never(unreachable)


def _group_elements(tag) -> "list[RigidMotion] | None":
    """Every element of a FINITE tag's group (memoised through production's ``_realize``), or ``None`` for a continuous one."""
    realization = _realize(tag)
    return list(realization.elements) if realization.is_finite else None


def stabilizer_order(orbit_certificate: OrbitCertificate) -> np.ndarray:
    r""":math:`|\mathrm{Stab}(x_i)|` for every node ``i`` of the certificate's point set.

    ``1`` exactly on the regular (free) points; the singular set is ``> 1``. Meaningful only when the certificate was built from the
    FULL group rather than a generating set. By orbit-stabilizer the orbit length is :math:`|G|` divided by it.
    """
    order = np.zeros(orbit_certificate.n_points, dtype=np.int64)
    for perm in orbit_certificate.permutations:
        order[perm.fixed_points] += 1
    return order


def permutes(rigid_motion: RigidMotion, points: "ArrayLike", *, atol: float) -> "Permutation | None":
    r"""The permutation :math:`\pi` with :math:`g(x_i) = x_{\pi(i)}`, or ``None`` if ``g`` does not map the set onto itself.

    Returns the permutation rather than a ``bool``: every downstream reading (orbits, stabilisers, the singular set) is a reading of
    :math:`\pi`, and a returned permutation makes its own bijectivity assertable.

    Two guards, and the second is the one that is easy to omit: every image is within ``atol`` of some point of the set; and the match is
    a **bijection** (nearest-neighbour matching alone proves only that every image has *a* partner; a bit-identical duplicate of one
    point is the minimal witness). Both live in production's ``permutation_between`` (ERR-073), which this function calls.
    """
    x = np.asarray(points, dtype=float)
    return permutation_between(x, rigid_motion.on_points(x), atol=atol)
