r"""Composable quadrature rule library returning :class:`DiscreteMeasure`.

This subpackage hosts the *publication-grade* quadrature rules used by
the project. Each function returns a
:class:`~orpheus.numerics.measure.DiscreteMeasure` carrying its
``invariance_group`` (per Issue 3 of
``.claude/plans/sn_reshape.md``) and its
``degree_of_exactness`` so downstream code can reason about polynomial
exactness and symmetry without recomputing it.

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
ground truth — plus the SN-side derived data (reflection partners,
octant partition, level structure) cached at construction time.

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
    GEOMETRY_GROUPS,
    QuadratureSelectionError,
    QuadratureSpec,
    SelectionLog,
    quadrature_registry,
    select_quadrature,
)
from .rules_1d import gauss_legendre_on_mu
from .rules_product import product_mu_phi
from .rules_sphere import (
    LevelStructure,
    lebedev_sphere,
    level_symmetric_sn,
)

__all__ = [
    "GEOMETRY_GROUPS",
    "LevelStructure",
    "Quadrature",
    "QuadratureSelectionError",
    "QuadratureSpec",
    "SelectionLog",
    "gauss_legendre_on_mu",
    "lebedev_sphere",
    "level_symmetric_sn",
    "product_mu_phi",
    "quadrature_registry",
    "select_quadrature",
]
