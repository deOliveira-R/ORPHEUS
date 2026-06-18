"""Per-cell spatial-update strategies for the SN sweep.

This package hosts the strategy contract and the concrete
strategies — Diamond Difference (here) and (Wave C-extension)
Linear Discontinuous, Exponential Characteristic, Step — that a
1-D SN sweep uses to march through a cell and produce its
average flux + outgoing states.

The contract itself lives in :mod:`orpheus.sn.spatial.scheme`:

* :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` — the
  ``@runtime_checkable Protocol``.
* :class:`~orpheus.sn.spatial.scheme.UpstreamState` —
  per-cell input state (spatial face flux + optional angular
  half-flux).
* :class:`~orpheus.sn.spatial.scheme.CellResult` — per-cell
  output state (average flux, outgoing spatial flux, outgoing
  angular state).

Round 1 (Issue #157) shipped the contract.
Round 2 (Issue #158, this file's ``DiamondDifference`` re-export)
ships the first concrete strategy as a bit-identical extraction of
the existing inlined sweep math at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``).  Wave
C-extension will add Linear Discontinuous, Step, and Exponential
Characteristic.

See the SN reshape campaign plan at
``.claude/plans/sn_reshape.md`` and the Wave C plan at
``.claude/plans/mossy-mapping-pine.md`` for context.
"""

from .scheme import CellResult, DiscretizationScheme, UpstreamState
from .diamond import DiamondDifference
from .linear_discontinuous import LinearDiscontinuous
from .pairing import pair_diffusion_limit_consistent
from .pole_angular_closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    PoleAngularClosure,
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
    "CellResult",
    "DiscretizationScheme",
    "DiamondDifference",
    "IdentityAngularClosure",
    "LinearDiscontinuous",
    "MorelMontryAngularSweep",
    "PoleAngularClosure",
    "PoleAngularClosureBase",
    "PsiHalfAngleSeed",
    "PsiHalfAngleSeedBase",
    "UpstreamState",
    "ZeroSeed",
    "default_angular_closure_class",
    "ordinate_scan",
    "pair_diffusion_limit_consistent",
]
