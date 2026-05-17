r"""Comparison test — unified matvec vs legacy spherical matvec (PR-TYPED-6c Step 2).

Issue #197 PR-TYPED-6c Step 2 — verifies that the new
:func:`~orpheus.sn.operator.transport_operator_matvec_unified` produces
output bit-exact (within ULP) to the legacy
:func:`~orpheus.sn.operator.transport_operator_matvec_spherical` for
spherical geometry.

The two functions consume DIFFERENT input layouts:
  - legacy: ``psi_packed`` shape ``(n_unknowns,)`` (F-order ravel of
    ``(ng, n_eq)``)
  - unified: ``psi_view`` shape ``(N, ng, nx, ny)`` canonical layout

Equivalence chain:
  1. Start with a canonical ``psi_view`` (N, ng, nx, ny).
  2. Decode to packed via the same path ``solution_to_angular_flux_spherical``
     uses (sparse scatter): ``psi_packed[g + ng*k] = psi_view[ord[k], g, ix[k], iy[k]]``.
  3. Run legacy matvec on ``psi_packed`` → ``m_packed (n_unknowns,)``.
  4. Decode ``m_packed`` to canonical ``m_legacy_view (N, ng, nx, ny)``
     using ``solution_to_angular_flux_spherical``.
  5. Run unified matvec on ``psi_view`` → ``m_unified (N, ng, nx, ny)``.
  6. Compare ``m_unified`` against ``m_legacy_view`` at the unknown
     slots — equivalence holds iff bit-exact within reduction-tree ULP.

The comparison is performed ONLY at unknown slots (eq_map.ordinate/ix/iy
positions). BC-resolved slots are filled by the legacy decoder via
BC.apply and may differ structurally — those are not part of the matvec's
domain.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map_spherical,
    solution_to_angular_flux_spherical,
    transport_operator_matvec_spherical,
    transport_operator_matvec_unified,
)
from orpheus.sn.quadrature import GaussLegendre1D
from tests.sn._test_helpers import placeholder_materials


def _build_sphere(n_cells: int = 5, n_ord: int = 4) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, n_cells + 1),
        mat_ids=np.zeros(n_cells, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials())


def _canonical_to_packed(
    psi_view: np.ndarray, eq_map, ng: int,
) -> np.ndarray:
    """Convert canonical (N, ng, nx, ny) to legacy F-order packed format."""
    # packed[g + ng*k] = psi_view[ord[k], g, ix[k], iy[k]]
    # Reshape via (ng, n_eq) and F-order ravel matches the legacy encoder.
    flux = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        flux[:, k] = psi_view[
            eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
        ]
    return flux.ravel(order='F')


def _bc_fill_outer(
    psi_view: np.ndarray, sn_mesh, eq_map,
) -> np.ndarray:
    """Make psi_view BC-consistent: BC-fill the inflow-at-outer-boundary slots.

    The legacy ``solution_to_angular_flux_spherical`` decoder calls
    ``bc_outer.apply(outgoing)`` to populate the BC-resolved incoming-
    ordinate slots at i=nx-1. The unified matvec reads psi_view as-is,
    so for bit-exact comparison the test input must already be BC-filled
    at those positions.
    """
    # Identify BC-resolved positions: outer cell (i=nx-1) for incoming ordinates.
    # For sphere with reflective BC, these are mu_x < 0 at i=nx-1.
    quad = sn_mesh.quad
    incoming_mask = quad.mu_x < -1e-15
    if not incoming_mask.any():
        return psi_view
    # Apply BC: outgoing → incoming at the outer face.
    outer_face = psi_view[:, :, -1, 0]  # (N, ng) — current full outer face
    inflow_full = sn_mesh.bc_right.apply(outer_face)  # (N, ng)
    psi_view = psi_view.copy()
    psi_view[incoming_mask, :, -1, 0] = inflow_full[incoming_mask, :]
    return psi_view


def _extract_at_unknown_slots(
    field_4d: np.ndarray, eq_map,
) -> np.ndarray:
    """Gather field_4d at the eq_map's (ordinate, ix, iy) slots → (ng, n_eq)."""
    ng = field_4d.shape[1]
    out = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        out[:, k] = field_4d[
            eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
        ]
    return out


pytestmark = [pytest.mark.l0]


class TestUnifiedMatvecSphere:
    """Unified matvec equivalent to legacy spherical matvec at ULP scale."""

    @pytest.mark.parametrize("n_cells", [5, 10, 20])
    @pytest.mark.parametrize("n_ord", [4, 8])
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_unified_matches_legacy_spherical(
        self, n_cells: int, n_ord: int, seed: int,
    ) -> None:
        sn_mesh = _build_sphere(n_cells=n_cells, n_ord=n_ord)
        ng = 1
        eq_map = build_equation_map_spherical(n_cells, sn_mesh.quad, ng)

        # Random psi in canonical (N, ng, nx, ny) layout.
        rng = np.random.default_rng(seed)
        psi_view = rng.standard_normal(
            (sn_mesh.quad.N, ng, n_cells, 1)
        ).astype(np.float64)
        # BC-fill the inflow-at-outer slots to match the legacy decoder's
        # behavior. Both unified and legacy will see the SAME 4-D input.
        psi_view = _bc_fill_outer(psi_view, sn_mesh, eq_map)

        sigma_t = np.full((ng, n_cells, 1), 2.0)

        # Legacy path: canonical → packed → legacy matvec → decode to (N, ng, nx, ny).
        psi_packed = _canonical_to_packed(psi_view, eq_map, ng)
        reduced = sn_mesh.reduced
        m_packed = transport_operator_matvec_spherical(
            psi_packed, eq_map, sn_mesh.quad, sigma_t, n_cells, ng,
            reduced.face_areas, sn_mesh.volumes,
            reduced.alpha_half, reduced.redist_dAw, reduced.tau_mm,
            sn_mesh=sn_mesh, bc_outer=sn_mesh.bc_right,
            pole_angular_closure=sn_mesh.pole_angular_closure,
        )
        # Decode m_packed to (N, ng, nx, ny) via the same path the
        # legacy uses on the read side.
        m_legacy_view = solution_to_angular_flux_spherical(
            m_packed, eq_map, sn_mesh.quad, n_cells, ng,
        )

        # Unified path: canonical → unified matvec.
        m_unified = transport_operator_matvec_unified(
            psi_view, sn_mesh, sigma_t,
        )

        # Compare at unknown slots only (BC-resolved slots may differ
        # because solution_to_angular_flux_spherical fills them via
        # BC.apply but the unified matvec doesn't write there).
        m_legacy_unknowns = _extract_at_unknown_slots(m_legacy_view, eq_map)
        m_unified_unknowns = _extract_at_unknown_slots(m_unified, eq_map)

        # Bit-exact at ULP scale.
        np.testing.assert_allclose(
            m_unified_unknowns, m_legacy_unknowns,
            rtol=1e-13, atol=1e-14,
            err_msg=(
                f"unified vs legacy spherical mismatch at n_cells={n_cells} "
                f"n_ord={n_ord} seed={seed}"
            ),
        )

    def test_unified_constant_psi_gives_sigma_t(self) -> None:
        """At ψ = constant on homogeneous reflective sphere, unified
        matvec returns σ_t · ψ (= 2.0 for σ_t = 2.0 here). Sanity check."""
        sn_mesh = _build_sphere(n_cells=5, n_ord=4)
        ng = 1
        sigma_t_val = 2.0
        sigma_t = np.full((ng, sn_mesh.nx, 1), sigma_t_val)
        psi_view = np.ones((sn_mesh.quad.N, ng, sn_mesh.nx, 1))

        m_unified = transport_operator_matvec_unified(
            psi_view, sn_mesh, sigma_t,
        )
        # At constant ψ = 1: (L+C)·1 ≈ σ_t · 1 = 2.0 everywhere.
        eq_map = build_equation_map_spherical(sn_mesh.nx, sn_mesh.quad, ng)
        m_at_unknowns = _extract_at_unknown_slots(m_unified, eq_map)
        np.testing.assert_allclose(
            m_at_unknowns, sigma_t_val, rtol=1e-13, atol=1e-14,
        )

    def test_unified_zero_psi_gives_zero(self) -> None:
        """Linear operator: zero input → zero output."""
        sn_mesh = _build_sphere(n_cells=5, n_ord=4)
        ng = 1
        sigma_t = np.full((ng, sn_mesh.nx, 1), 2.0)
        psi_view = np.zeros((sn_mesh.quad.N, ng, sn_mesh.nx, 1))

        m_unified = transport_operator_matvec_unified(
            psi_view, sn_mesh, sigma_t,
        )
        np.testing.assert_array_equal(
            m_unified, np.zeros_like(m_unified),
        )
