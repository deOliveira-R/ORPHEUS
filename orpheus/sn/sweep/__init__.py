"""The discrete-ordinates sweep machinery: closures, seeds, caches, the scan.

This package is the sweep's kernel substrate — the pieces the walk
executors (:mod:`orpheus.sn.loss_representation`) and the sweep-family
operators (:mod:`orpheus.sn.operators`) consume per ordinate.  The
method-generic spatial closure layer (the ``DiscretizationScheme``
contract, Diamond Difference, Linear Discontinuous, the UBLD kernels,
and the per-cell WDD balance) lives in :mod:`orpheus.transport.spatial`
since #272:

* :mod:`~orpheus.sn.sweep.pole_angular_closure` — the curvilinear
  angular-redistribution closure family (Morel–Montry, identity).
* :mod:`~orpheus.sn.sweep.psi_half_angle_seed` — the Hébert §3.9.4
  starting-direction DD march (the #282 route-(a) direct ψ½ solver).
* :mod:`~orpheus.sn.sweep.cache` — the two-stratum hot-path
  cache keyed to solver/mesh × quadrature lifecycles.
* :mod:`~orpheus.sn.sweep.scan` — the 1-D Blelloch ordinate scan and
  the Cartesian row-march.
* :mod:`~orpheus.sn.sweep.pairing` — pairing-validity predicates for
  the (spatial ⊗ angular) discretization product (#236).

History: this package was ``orpheus.sn.spatial`` until the R9 estate
rename (task #54, 2026-07-13) — the old name described the pre-#272
contents (the spatial trial-space layer), not this residual.
"""

from .pairing import pair_diffusion_limit_consistent
from .pole_angular_closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    PoleAngularClosureBase,
    default_angular_closure_class,
)
from .psi_half_angle_seed import carlson_inward_sweep_from_source
from .scan import ordinate_scan

# Issue #168 Phase C retired the BoundaryFaceFlux Protocol entirely;
# the sweep-frame apply matvec subsumed the boundary closure into
# the WDD propagation chain.
#
# Issue #168 Phase D shipped a PsiHalfAngleSeed strategy family;
# ERR-058 (#195) flipped its default to the angular-edge
# extrapolation.  #282 route (a) (#280 Phase 2.5d, 2026-07-04)
# RETIRED the whole strategy zoo: the starting-direction flux is
# first-class composite STATE (R12a presence predicate), the solve
# marches it directly via ``carlson_inward_sweep_from_source`` on the
# TRUE q½ source, the apply reads the given carrier block, and the
# non-carrying cylinder levels inline the edge extrapolation on the
# closure (``MorelMontryAngularSweep.edge_extrapolated_seed``).

__all__ = [
    "IdentityAngularClosure",
    "MorelMontryAngularSweep",
    "PoleAngularClosureBase",
    "carlson_inward_sweep_from_source",
    "default_angular_closure_class",
    "ordinate_scan",
    "pair_diffusion_limit_consistent",
]
