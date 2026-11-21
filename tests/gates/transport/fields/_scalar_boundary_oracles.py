"""Test-side oracles of :class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`.

``p1_boundary_scalar_flux`` is the P1 closure of the Cauchy pair, ``φ_Γ = 2(J⁺ + J⁻)``, written from the
public ``outflow_view`` and ``inflow_view``. It was a method of the field with no production caller (retired
2026-11-20); the two gate consumers, one transport and one diffusion, read it here.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux


def p1_boundary_scalar_flux(scalar_boundary_flux: ScalarBoundaryFlux, face: str) -> NDArray:
    r"""P1 boundary scalar flux ``φ_Γ = 2(J⁺ + J⁻)`` on ``face`` (a copy)."""
    return np.asarray(2.0 * (scalar_boundary_flux.outflow_view(face) + scalar_boundary_flux.inflow_view(face)))
