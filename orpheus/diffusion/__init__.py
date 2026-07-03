"""Multigroup diffusion on the operator algebra (#290).

The CANONICAL diffusion API — the scalar-composite operator family
(P4), the boundary realizer (P3), and the k-eigenvalue solver on the
shared engines (P5):

* ``solve_diffusion_1d(materials, mesh)`` — the modern driver;
* ``DiffusionSolver`` / ``DiffusionResult`` — the
  ``EigenvalueSolver``-protocol implementer and its typed result;
* ``LeakageOperator`` (``L``, FULL) / ``DiffusionBoundaryOperator``
  (``B``, BOUNDARY) — the two diffusion-owned operator leaves
  (``C``/``S``/``F`` are the shared transport operators);
* ``DiffusionBoundaryRealizer`` / ``DiffusionMethodSpace`` — the
  albedo-family law realization ``J⁻ = 𝒜·J⁺`` (closes #182).

Importing the package auto-registers the realizer.

The legacy MATLAB-port island (``CoreGeometry`` / ``TwoGroupXS`` / the
BiCGSTAB solver) was retired at #290 P6; ``orpheus.diffusion.solver``
now IS the modern module (family naming parity with sn/cp/homogeneous).
"""

# #290 P5 — the modern k-eigenvalue solver on the operator algebra
# (renamed k_eigenvalue.py → solver.py at P6 when the island retired).
from .solver import DiffusionResult, DiffusionSolver, solve_diffusion_1d

# #290 P3 (closes #182) -- the FUNCTIONAL DiffusionBoundaryRealizer
# (albedo family J⁻ = 𝒜·J⁺ on the scalar partial-current trace) +
# its method space. Importing the module auto-registers the realizer.
from .boundary_realizer import DiffusionBoundaryRealizer
from .method_space import DiffusionMethodSpace

# #290 P4 -- the diffusion operator family's two new leaves on the
# scalar composite: the elliptic leakage L (FULL) and the realized
# boundary block B (BOUNDARY). C/S/F are the shared transport operators
# (MultiplicationOperator / K_iso / FissionOperator), which gained
# their scalar-composite arms in the same phase.
from .operators import DiffusionBoundaryOperator, LeakageOperator

__all__ = [
    "DiffusionBoundaryOperator",
    "DiffusionBoundaryRealizer",
    "DiffusionMethodSpace",
    "DiffusionResult",
    "DiffusionSolver",
    "LeakageOperator",
    "solve_diffusion_1d",
]
