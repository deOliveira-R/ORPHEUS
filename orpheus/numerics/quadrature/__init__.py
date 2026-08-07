r"""Composable quadrature rule library returning :class:`DiscreteMeasure`.

This subpackage hosts the *publication-grade* quadrature rules used by
the project. Each function returns a
:class:`~orpheus.numerics.measure.DiscreteMeasure` carrying its
``invariance_group`` (per Issue 3 of
``.claude/plans/sn_reshape.md``) and its
``degree_of_exactness`` so downstream code can reason about polynomial
exactness and symmetry without recomputing it.

The ``invariance_group`` tag is a **construction-time** fact — the
factory built the node set, so it knows its symmetry — and a permanent
gate (``test_every_registry_rule_declares_a_symmetry_it_actually_has``)
pins every shipped rule's tag against
:func:`~orpheus.numerics.symmetry.maximal_invariance_groups`, which
computes the group from the nodes. Anything that *selects* on symmetry
asks the nodes rather than the tag: a tag can be true without being
maximal, and until 2026-08-02 one of them was not even true.

The canonical directional :class:`Quadrature`
---------------------------------------------

R-1 Phase A detour-C retired the per-family SN-side adapter classes
(GaussLegendre1D / LebedevSphere / LevelSymmetricSN /
ProductQuadrature) and the AngularQuadrature Protocol. The four
families collapse into named classmethod factories on a single
:class:`Quadrature` class:

* :meth:`Quadrature.gauss_legendre(n)` — slab 1-D polar.
* :meth:`Quadrature.lebedev(order)` — :math:`O_h`-invariant sphere.
* :meth:`Quadrature.level_symmetric(sn_order)` — Carlson-Lathrop S_N.
* :meth:`Quadrature.product(n_mu, n_phi)` — GL × equispaced.

Each factory returns a :class:`Quadrature` wrapping the canonical
:class:`DiscreteMeasure` produced by the rule functions
(:func:`gauss_legendre_on_mu`, :func:`lebedev_sphere`,
:func:`level_symmetric_sn`, :func:`product_mu_phi`) — the mathematical
ground truth — plus the SN-side derived-data surface (ordinate
permutations via :meth:`Quadrature.ordinate_permutation`, the octant
partition, the level structure).

The legacy ``mu_x`` / ``mu_y`` / ``mu_z`` / ``weights`` / ``N``
attribute surface survives on :class:`Quadrature` as ``@property``
views over the underlying ``measure.nodes`` / ``measure.weights``,
preserving back-compat for the ~150 consumer call sites that still
read by these names. The Pattern 7 violation that made the legacy
surface dangerous (denormalised dataclass *fields* in four parallel
adapter classes) is closed by construction here: there is exactly
ONE producer of the ordinate data, the ``measure``; the named
accessors are derived views with no separate storage.

References
----------

* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. Theorem 3.6.20 (Gauss-Legendre
  exactness).
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan et al., eds., Gordon & Breach.
  Level-symmetric :math:`S_N` construction; :math:`O_h`-invariant.
* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Comp. Math.*
  **16**(2), 10-24. The :math:`O_h`-invariant sphere quadrature.

See Also
--------
:ref:`discrete-measures` (theory page) — composition algebra,
metadata-propagation table, the bridge pattern documented here.
"""

from .directional import Quadrature
from .registry import (
    GEOMETRY_ANGULAR_SYMMETRY,
    AngularSymmetry,
    QuadratureSelectionError,
    QuadratureSpec,
    SelectionLog,
    quadrature_registry,
    select_quadrature,
)
from .rules_1d import gauss_legendre_on_mu
from .rules_circle import NODE_ALIGNED, STAGGERED, periodic_trapezoid
from .rules_product import product_mu_phi, spherical_product
from .rules_sphere import (
    LevelStructure,
    lebedev_sphere,
    level_symmetric_sn,
)

__all__ = [
    "GEOMETRY_ANGULAR_SYMMETRY",
    "NODE_ALIGNED",
    "STAGGERED",
    "AngularSymmetry",
    "LevelStructure",
    "Quadrature",
    "QuadratureSelectionError",
    "QuadratureSpec",
    "SelectionLog",
    "gauss_legendre_on_mu",
    "lebedev_sphere",
    "level_symmetric_sn",
    "periodic_trapezoid",
    "product_mu_phi",
    "quadrature_registry",
    "select_quadrature",
    "spherical_product",
]
