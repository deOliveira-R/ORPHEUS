r"""Slab Carlvik-Galerkin solver — even-Legendre basis."""
from __future__ import annotations

from .one_group_anisotropic import (
    CarlvikGalerkinSlabResult,
    solve_galerkin_spectral_slab,
)

__all__ = [
    "CarlvikGalerkinSlabResult",
    "solve_galerkin_spectral_slab",
]
