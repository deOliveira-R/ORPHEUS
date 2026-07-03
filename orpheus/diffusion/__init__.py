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
