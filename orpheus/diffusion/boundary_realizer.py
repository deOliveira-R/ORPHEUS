r"""Diffusion BoundaryRealizer — STUB.

Stub realizer for the diffusion solver. Registered in
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` under
``method_name = "diffusion"`` so the architecture is wired end-to-end,
but :meth:`realize` raises :class:`NotImplementedError`. Functional
implementation lands when diffusion adopts the unified BC
architecture — see follow-up issue ``BC: DiffusionBoundaryRealizer
functional implementation`` filed at the end of the 12-wave refactor.

Diffusion BCs are Marshak-extrapolated forms of the transport
BCs (Marshak 1947; Bell & Glasstone 1970 §3.4): vacuum becomes
:math:`\phi + (2/3)\,\nabla\phi \cdot \hat n = 0` at the surface
(albedo ≈ 0); reflective becomes :math:`\nabla\phi \cdot \hat n = 0`;
albedo and white share a generalised mixed-extrapolation form
parametrised by the half-current ratio. The realizer here will map
each legacy BC onto the appropriate Robin-coefficient block of the
diffusion finite-volume matrix.
"""

from __future__ import annotations

from typing import Any

from orpheus.geometry.boundary import BoundaryRealizerRegistry


__all__ = ["DiffusionBoundaryRealizer"]


@BoundaryRealizerRegistry.register("diffusion")
class DiffusionBoundaryRealizer:
    """Stub realizer — raises :class:`NotImplementedError`."""

    method_name: str = "diffusion"

    def realize(self, law: Any, method_space: Any):
        raise NotImplementedError(
            f"DiffusionBoundaryRealizer.realize({type(law).__name__}, ...) "
            f"is not yet implemented. The stub is registered so the "
            f"architecture is wired end-to-end. File issue when "
            f"diffusion adopts the unified BC architecture; see "
            f".claude/plans/transient-giggling-cake.md Wave 5 / C5.4."
        )
