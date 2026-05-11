from .solver import (
    CoreGeometry,
    DiffusionResult,
    DiffusionSolver,
    TwoGroupXS,
    solve_diffusion_1d,
)

# Wave 5 -- auto-register the (stub) DiffusionBoundaryRealizer.
from . import boundary_realizer  # noqa: F401
