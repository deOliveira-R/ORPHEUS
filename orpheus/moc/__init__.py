from .solver import MoCResult, solve_moc
from .geometry import MOCMesh, Segment, Track
from .quadrature import MOCQuadrature

# Wave 5 -- auto-register the (stub) MoCBoundaryRealizer so
# `BoundaryRealizerRegistry.get("MoC")` succeeds after a plain
# `import orpheus.moc`. The stub raises NotImplementedError on
# realize(); the registration wires the architecture end-to-end.
from . import boundary_realizer  # noqa: F401
