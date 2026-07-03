"""Multigroup diffusion on the operator algebra (#290).

The CANONICAL diffusion API — the scalar-composite operator family
(P4), the boundary realizer (P3), and the k-eigenvalue solver on the
shared engines (P5):

* ``solve_diffusion_1d(materials, mesh)`` — the modern driver;
* ``DiffusionMesh`` — the diffusion phase space (P7a): MaterialMesh +
  scalar ``(J⁺, J⁻)`` trace + boundary laws realized at construction;
* ``DiffusionSolver`` / ``DiffusionResult`` — the
  ``EigenvalueSolver``-protocol implementer and its typed result;
* ``LeakageOperator`` (``L``, FULL) / ``DiffusionBoundaryOperator``
  (``B``, BOUNDARY) — the two diffusion-owned operator leaves
  (``C``/``S``/``F`` are the shared transport operators);
* ``DiffusionBoundaryRealizer`` / ``DiffusionMethodSpace`` — the
  albedo-family law realization ``J⁻ = 𝒜·J⁺`` (closes #182).

``DiffusionMesh`` realizes its boundary laws at construction through the
shared :func:`orpheus.transport.method.resolve_boundary_conditions` body
(the string-keyed realizer registry was dissolved at #290 P7b).

The legacy MATLAB-port island (``CoreGeometry`` / ``TwoGroupXS`` / the
BiCGSTAB solver) was retired at #290 P6; ``orpheus.diffusion.solver``
now IS the modern module (family naming parity with sn/cp/homogeneous).
"""

# #290 P7a -- the diffusion method-mesh (MaterialMesh + scalar trace +
# realized boundary laws; the SNMesh sibling).
from .augmented_mesh import DiffusionMesh

# #290 P5 — the modern k-eigenvalue solver on the operator algebra
# (renamed k_eigenvalue.py → solver.py at P6 when the island retired).
from .solver import DiffusionResult, DiffusionSolver, solve_diffusion_1d

# #290 P3 (closes #182) -- the FUNCTIONAL DiffusionBoundaryRealizer
# (albedo family J⁻ = 𝒜·J⁺ on the scalar partial-current trace) +
# its method space. Owned directly by DiffusionMesh (P7b dissolved the
# realizer registry — no import-side-effect registration remains).
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
    "DiffusionMesh",
    "DiffusionMethodSpace",
    "DiffusionResult",
    "DiffusionSolver",
    "LeakageOperator",
    "solve_diffusion_1d",
]
