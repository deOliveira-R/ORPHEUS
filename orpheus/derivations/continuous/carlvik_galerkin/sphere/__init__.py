r"""Sphere Carlvik-Galerkin solver — odd-Legendre basis (for r·φ)."""
from __future__ import annotations

from .one_group_anisotropic import (
    CarlvikGalerkinSphereResult,
    solve_carlvik_galerkin_sphere,
)

__all__ = [
    "CarlvikGalerkinSphereResult",
    "solve_carlvik_galerkin_sphere",
]
