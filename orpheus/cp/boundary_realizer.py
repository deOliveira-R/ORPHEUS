r"""CP BoundaryRealizer — STUB.

Stub realizer for Collision Probabilities. Registered in
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` under
``method_name = "CP"`` so the architecture is wired end-to-end, but
:meth:`realize` raises :class:`NotImplementedError`. Functional
implementation lands when CP adopts the unified BC architecture —
see follow-up issue ``BC: CPBoundaryRealizer functional implementation``
filed at the end of the 12-wave refactor.

The CP realisation maps each legacy BC onto the appropriate
boundary-to-region or region-to-boundary block of the
collision-probability matrix: vacuum is the zero block; specular
folds rays back into the kernel via the reflected partner regions;
white redistributes uniformly across all entering boundary
segments. The CP "method space" will carry the region/segment
indexing the realizer needs.
"""

from __future__ import annotations

from typing import Any

from orpheus.geometry.boundary import BoundaryRealizerRegistry


__all__ = ["CPBoundaryRealizer"]


@BoundaryRealizerRegistry.register("CP")
class CPBoundaryRealizer:
    """Stub realizer — raises :class:`NotImplementedError`."""

    method_name: str = "CP"

    def realize(self, law: Any, method_space: Any):
        raise NotImplementedError(
            f"CPBoundaryRealizer.realize({type(law).__name__}, ...) "
            f"is not yet implemented. The stub is registered so the "
            f"architecture is wired end-to-end. File issue when CP "
            f"adopts the unified BC architecture; see "
            f".claude/plans/transient-giggling-cake.md Wave 5 / C5.4."
        )
