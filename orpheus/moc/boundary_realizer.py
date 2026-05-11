r"""MoC BoundaryRealizer — STUB.

Stub realizer for the Method of Characteristics. Registered in
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` under
``method_name = "MoC"`` so the architecture is wired end-to-end, but
:meth:`realize` raises :class:`NotImplementedError`. Functional
implementation lands when MoC adopts the unified BC architecture —
see follow-up issue ``BC: MoCBoundaryRealizer functional implementation``
filed at the end of the 12-wave refactor.

The MoC realisation of each legacy BC will differ structurally from
the SN realizer: vacuum kills entering track boundary fluxes (no
ordinate concept); specular reflects the track direction at the
incident face; white redistributes outgoing track current
isotropically into the entering tracks of the same face. The
stub registration today is a forward-reference; the functional
implementation lands during the MoC modernisation phase.
"""

from __future__ import annotations

from typing import Any

from orpheus.geometry.boundary import BoundaryRealizerRegistry


__all__ = ["MoCBoundaryRealizer"]


@BoundaryRealizerRegistry.register("MoC")
class MoCBoundaryRealizer:
    """Stub realizer — raises :class:`NotImplementedError`."""

    method_name: str = "MoC"

    def realize(self, law: Any, method_space: Any):
        raise NotImplementedError(
            f"MoCBoundaryRealizer.realize({type(law).__name__}, ...) "
            f"is not yet implemented. The stub is registered so the "
            f"architecture is wired end-to-end. File issue when MoC "
            f"adopts the unified BC architecture; see "
            f".claude/plans/transient-giggling-cake.md Wave 5 / C5.4."
        )
