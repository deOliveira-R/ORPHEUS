"""SN sweep-walk support: angular closures, seeds, caches, and the scan.

What remains here after #272 promoted the method-generic spatial
closure layer (the ``DiscretizationScheme`` contract, Diamond
Difference, Linear Discontinuous, the UBLD kernels, and the per-cell
WDD balance) to :mod:`orpheus.transport.spatial` is the
discrete-ordinates sweep machinery proper:

* :mod:`~orpheus.sn.spatial.pole_angular_closure` — the curvilinear
  angular-redistribution closure family (Morel–Montry, identity).
* :mod:`~orpheus.sn.spatial.psi_half_angle_seed` — starting-direction
  seed strategies for the curvilinear half-angle level.
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
from .psi_half_angle_seed import (
    AngularEdgeExtrapolation,
    CarlsonInwardSweep,
    CarlsonSweepContext,
    PsiHalfAngleSeed,
    PsiHalfAngleSeedBase,
    ZeroSeed,
)
from .scan import ordinate_scan

# Issue #168 Phase C retired the BoundaryFaceFlux Protocol entirely;
# the sweep-frame apply matvec subsumed the boundary closure into
# the WDD propagation chain.
#
# Issue #168 Phase D shipped the PsiHalfAngleSeed strategy family —
# composed into MorelMontryAngularSweep via Option α (composition,
# not sibling Protocol).  ERR-058 (#195, 2026-06-12) flipped the
# default from CarlsonInwardSweep (proxy-source seed, exact only at
# flat-flux equilibrium) to AngularEdgeExtrapolation (the input
# field extrapolated in μ to the level's angular edge — the
# operator-consistent seed).  CarlsonInwardSweep stays registered as
# the Hébert §3.9.4 recurrence host for future TRUE-source-driven
# sweep-side seeding.

__all__ = [
    "AngularEdgeExtrapolation",
    "CarlsonInwardSweep",
    "CarlsonSweepContext",
    "IdentityAngularClosure",
    "MorelMontryAngularSweep",
    "PoleAngularClosureBase",
    "PsiHalfAngleSeed",
    "PsiHalfAngleSeedBase",
    "ZeroSeed",
    "default_angular_closure_class",
    "ordinate_scan",
    "pair_diffusion_limit_consistent",
]
