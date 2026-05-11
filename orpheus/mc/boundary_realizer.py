r"""MC BoundaryRealizer — STUB.

Stub realizer for Monte Carlo. Registered in
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` under
``method_name = "MC"`` so the architecture is wired end-to-end, but
:meth:`realize` raises :class:`NotImplementedError`. Functional
implementation lands when MC adopts the unified BC architecture —
see follow-up issue ``BC: MCBoundaryRealizer functional implementation``
filed at the end of the 12-wave refactor.

The MC realisation differs from SN/MoC at a categorical level: MC
doesn't have a linear operator acting on a flux vector; it has a
*particle event* (kill / reflect / re-emit) that fires when a
random walk reaches a boundary. The "realized operator" returned by
:meth:`realize` will therefore be a thin wrapper around an event
callback rather than a sparse-matrix-like
:class:`~orpheus.numerics.operator.LinearOperator`. The unified
interface is meaningful precisely because the operator algebra
type system carries the abstraction; MC ships against the same
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances as
the deterministic methods.
"""

from __future__ import annotations

from typing import Any

from orpheus.geometry.boundary import BoundaryRealizerRegistry


__all__ = ["MCBoundaryRealizer"]


@BoundaryRealizerRegistry.register("MC")
class MCBoundaryRealizer:
    """Stub realizer — raises :class:`NotImplementedError`."""

    method_name: str = "MC"

    def realize(self, law: Any, method_space: Any):
        raise NotImplementedError(
            f"MCBoundaryRealizer.realize({type(law).__name__}, ...) "
            f"is not yet implemented. The stub is registered so the "
            f"architecture is wired end-to-end. File issue when MC "
            f"adopts the unified BC architecture; see "
            f".claude/plans/transient-giggling-cake.md Wave 5 / C5.4."
        )
