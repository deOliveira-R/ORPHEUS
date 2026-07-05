"""SN sweep-walk support: angular closures, seeds, caches, and the scan.

What remains here after #272 promoted the method-generic spatial
closure layer (the ``DiscretizationScheme`` contract, Diamond
Difference, Linear Discontinuous, the UBLD kernels, and the per-cell
WDD balance) to :mod:`orpheus.transport.spatial` is the
discrete-ordinates sweep machinery proper:

* :mod:`~orpheus.sn.spatial.pole_angular_closure` — the curvilinear
  angular-redistribution closure family (Morel–Montry, identity).
* :mod:`~orpheus.sn.spatial.psi_half_angle_seed` — the Hébert §3.9.4
  starting-direction DD march (the #282 route-(a) direct ψ½ solver).
* :mod:`~orpheus.sn.spatial.sweep_cache` — the two-stratum hot-path
  cache keyed to solver/mesh × quadrature lifecycles.
* :mod:`~orpheus.sn.spatial.scan` — the 1-D Blelloch ordinate scan and
  the Cartesian row-march.
* :mod:`~orpheus.sn.spatial.pairing` — pairing-validity predicates for
  the (spatial ⊗ angular) discretization product (#236).

NOTE (R9, 2026-07-04): the package NAME is stale for this residual —
what lives here is sweep-walk + angular machinery, not spatial trial
space.  The rename is deliberately deferred to the #280
orientation×kernel walk unification (Phase 2.5), which owns the
sweep-layer layout.
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
