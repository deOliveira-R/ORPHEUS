from .solver import (
    CoreGeometry,
    DiffusionResult,
    DiffusionSolver,
    TwoGroupXS,
    solve_diffusion_1d,
)

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
