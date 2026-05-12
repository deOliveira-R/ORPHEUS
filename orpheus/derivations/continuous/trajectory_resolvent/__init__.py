"""Peierls Variant α Green's-function references.

This sub-package hosts the angle-resolved Green's-function references
produced via the Peierls Variant α formulation (Plan-2). The Variant α
form interpolates between vacuum (:math:`\\alpha = 0`) and specular
(:math:`\\alpha = 1`) boundary conditions on a single rank-1 resolvent,
sharing the closure operator across geometry kernels.

The package is the continuous-reference twin of
:mod:`orpheus.derivations.continuous.peierls_nystrom`: both descend from
the Peierls integral form, but the Nyström family discretizes the
*scalar-flux* integral equation while the Variant α family discretizes
the *angle-resolved* Green's function and reduces it analytically via
the rank-1 trick.

Sub-modules:

- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core` —
  the shared rank-1 Variant α resolvent and surface-flux closure
  primitives consumed by both sphere and cylinder solvers.
- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function` —
  sphere production solver (fixed-source Green's function and
  multi-region eigenvalue/k-inf evaluation).
- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder` —
  cylinder production solver (mirror of the sphere driver, mounted on
  the same Variant α core).
- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.origins` —
  symbolic *origins* (SymPy V_α1..V_α3 operator-level identities) for
  sphere and cylinder.
- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.billiard` —
  the :class:`Billiard` class (R2 hindsight refactor) — a math-rich
  facade over every Variant α geometry that returns the SHARED
  cross-method :class:`CriticalSolution` /
  :class:`FluxSolution` result types from
  :mod:`orpheus.derivations.common.solution_types`.
- :mod:`~orpheus.derivations.continuous.trajectory_resolvent.chord_oracle` —
  the :class:`ChordOracle` Protocol + per-geometry implementations
  (R3 hindsight refactor) — extracts the chord-arithmetic primitives
  (the *base atlas* in the fiber-bundle frame) into one Protocol with
  six concrete frozen-dataclass instances, one per geometry.

Public API entry point: :class:`Billiard`. Construct directly via
``Billiard(geometry=structured_geometry, materials={0: mix}, alpha=...)``,
then call :meth:`Billiard.solve_critical` for *k*-eigenproblems or
:meth:`Billiard.solve_fixed_source` for fixed-source problems. The
legacy ``solve_greens_function_*`` entry points remain available for
back-compat.
"""

from orpheus.derivations.continuous.trajectory_resolvent.billiard import (
    Billiard,
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    AnnulusChordOracle,
    ChordOracle,
    CylinderChordOracle,
    HollowSphereChordOracle,
    MultiRegionCylinderChordOracle,
    MultiRegionSphereChordOracle,
    SlabAsymmetricChordOracle,
    SphereChordOracle,
)

__all__ = [
    "Billiard",
    "CriticalSolution",
    "FluxSolution",
    # R3 ChordOracle Protocol + per-geometry concrete implementations.
    "ChordOracle",
    "SphereChordOracle",
    "MultiRegionSphereChordOracle",
    "CylinderChordOracle",
    "MultiRegionCylinderChordOracle",
    "SlabAsymmetricChordOracle",
    "HollowSphereChordOracle",
    "AnnulusChordOracle",
]
