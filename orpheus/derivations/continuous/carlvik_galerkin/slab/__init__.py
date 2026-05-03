r"""Slab Carlvik-Galerkin solver — even-Legendre basis."""
from __future__ import annotations

from .one_group_anisotropic import (
    CarlvikGalerkinSlabResult,
    solve_carlvik_galerkin_slab,
)

__all__ = [
    "CarlvikGalerkinSlabResult",
    "solve_carlvik_galerkin_slab",
]
