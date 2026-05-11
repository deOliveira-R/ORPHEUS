from .solver import CPMesh, CPParams, CPResult, CPSolver, solve_cp

# Wave 5 -- auto-register the (stub) CPBoundaryRealizer.
from . import boundary_realizer  # noqa: F401
