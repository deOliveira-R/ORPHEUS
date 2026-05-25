r"""Composable quadrature rule library returning :class:`DiscreteMeasure`.

This subpackage hosts the *publication-grade* quadrature rules used by
the project. Each function returns a
:class:`~orpheus.numerics.measure.DiscreteMeasure` carrying its
``invariance_group`` (per Issue 3 of
``.claude/plans/sn_reshape.md``) and its
``degree_of_exactness`` so downstream code can reason about polynomial
exactness and symmetry without recomputing it.

The bridge pattern
------------------

Domain-specific consumers (e.g. SN) need fast attribute access to
``mu_x`` / ``mu_y`` / ``mu_z`` / ``weights`` arrays in hot inner
loops, plus precomputed reflection-index arrays. The SN-side adapters
at :mod:`orpheus.sn.quadrature` cache numpy views of the
:class:`DiscreteMeasure` returned here, preserving O(1) attribute
access for the ~50 consumer call sites (sweeps, BiCGSTAB operators,
solvers, geometry meshes) that depend on the legacy attribute API.

This separation gives:

* **One algebra-of-record** — node and weight construction lives
  here; :class:`DiscreteMeasure` is the canonical mathematical
  object.
* **Composability** — rules can be tensor-producted
  (:meth:`DiscreteMeasure.__mul__`), restricted, pushed forward, etc.,
  without inheriting any SN-specific surface.
* **Backward-compatible SN** — the :class:`AngularQuadrature`
  Protocol and its four implementations
  (:class:`~orpheus.sn.quadrature.GaussLegendre1D`,
  :class:`~orpheus.sn.quadrature.LebedevSphere`,
  :class:`~orpheus.sn.quadrature.LevelSymmetricSN`,
  :class:`~orpheus.sn.quadrature.ProductQuadrature`) keep their
  attribute-access surface bit-identical, with Wave B's regression
  snapshots at ``tests/sn/regression/snapshots/`` enforcing the
  contract.

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
