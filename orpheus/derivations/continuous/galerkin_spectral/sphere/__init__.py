r"""Sphere Carlvik-Galerkin solver — odd-Legendre basis (for r·φ)."""
from __future__ import annotations

from .one_group_anisotropic import (
    CarlvikGalerkinSphereResult,
    solve_galerkin_spectral_sphere,
)

__all__ = [
    "CarlvikGalerkinSphereResult",
    "solve_galerkin_spectral_sphere",
]
