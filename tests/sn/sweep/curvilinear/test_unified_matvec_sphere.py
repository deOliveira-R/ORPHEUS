r"""Sanity gates for the spherical ``(L + C)`` matvec.

.. warning::

   **This file NO LONGER carries a unified-vs-legacy comparison.**  It
   was created for Issue #197 PR-TYPED-6c Step 2 as exactly that: a
   bit-identity (within-ULP) gate of the then-new
   ``_transport_operator_matvec_unified`` against the legacy
   ``transport_operator_matvec_spherical``, via a six-step encode /
   run-legacy / decode / compare chain at the equation-bearing slots.
   PR-TYPED-6c Step 7 (2026-05-18) RETIRED
   ``transport_operator_matvec_spherical``, and the bit-identity test
   retired with it (see :class:`TestUnifiedMatvecSphere`); the Wave T
   T.5 matvec retirement then deleted
   ``_transport_operator_matvec_unified`` itself.  Do not read the file
   name — nor the orphaned ``_canonical_to_packed`` / ``_bc_fill_outer``
   helpers, which are un-called residue of the retired encode/decode
   chain — as evidence of a cross-implementation comparison.  BOTH
   sides of that comparison are gone.

What survives here are two L0 sanity gates on the ONE production
matvec, reached through the operator algebra
(:meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.apply`,
whose 1-D kernel is
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`):

1. ``ψ ≡ const`` on a homogeneous reflective sphere → ``σ_t·ψ`` at
   every equation-bearing slot (the per-ordinate flat-flux invariant).
2. ``ψ ≡ 0`` → ``0`` (linearity sentinel).

``legacy_proxy_matvec`` is the pre-B1'' *boundary-fill convention*
bridge (cell-centre face proxy), NOT retired code: it lets these gates
feed the bare ``(N, ng, nx, ny)`` arrays the L0 hand references were
built around.  Both gates assert ONLY at the equation-bearing slots
(all ``(n, ix, 0)`` except inward ordinates at ``ix == nx - 1``, whose
values the reflective BC determines).

The matvec's structurally-independent L0 reference — and its L1
``k_∞`` / trajectory-resolvent anchors — live in
``test_unified_matvec_cylinder.py``; the slab L1 ``k_∞`` anchor is in
``test_unified_matvec_slab.py``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import _LC_matvec
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import legacy_proxy_matvec, placeholder_materials


def _build_sphere(n_cells: int = 5, n_ord: int = 4) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, n_cells + 1),
        mat_ids=np.zeros(n_cells, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
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
    inflow_full = sn_mesh.bc["xmax"].apply(outer_face)  # (N, ng)
    psi_view = psi_view.copy()
    psi_view[incoming_mask, :, -1, 0] = inflow_full[incoming_mask, :]
    return psi_view


def _extract_at_unknown_slots(
    field_4d: np.ndarray, sn_mesh: SNMesh,
) -> np.ndarray:
    """Gather field_4d at the curvilinear equation-bearing slots → (ng, n_eq).

    D-J (2026-05-30) — replaces the legacy :class:`EquationMap`-driven
    slot map.  Curvilinear 1-D equation set: all ``(n, ix, 0)`` slots
    EXCEPT inward ordinates at the outermost cell ``ix == nx - 1``
    (the reflective BC determines those values; they are NOT unknowns).
    The quad-derived mask is legitimised by the (now-retired)
    ``test_native_matvec.TestQuadDerivedMaskEqualsLegacySlotMap`` pin
    which proved ``quad.mu_x > 0`` ≡ ``eq_map.face_outer_ordinate``.
    """
    quad = sn_mesh.quad
    nx = sn_mesh.nx
    ng = field_4d.shape[1]
    inflow_outer = quad.mu_x < -1e-15  # (N,)
    cols = []
    for ix in range(nx):
        for n in range(quad.N):
            if ix == nx - 1 and inflow_outer[n]:
                continue
            cols.append(field_4d[n, :, ix])
    return np.stack(cols, axis=1)  # (ng, n_eq)


pytestmark = [pytest.mark.l0]


class TestUnifiedMatvecSphere:
    """Unified matvec sphere sanity gates.

    PR-TYPED-6c Step 7 (2026-05-18) retired
    ``transport_operator_matvec_spherical``; the pre-retirement
    ``test_unified_matches_legacy_spherical`` ULP-bit-identity gate
    served its Step 2 purpose (verify the unified matvec matches the
    legacy helper at ULP) and retires with the helper.  Matvec
    correctness now relies on the constant-flux and zero-input
    sanity gates here plus the broader L0/L1 anchors elsewhere
    (the unified matvec's L1 reference for cylinder is the
    trajectory_resolvent test in ``test_unified_matvec_cylinder.py``).
    """

    def test_unified_constant_psi_gives_sigma_t(self) -> None:
        """At ψ = constant on homogeneous reflective sphere, unified
        matvec returns σ_t · ψ (= 2.0 for σ_t = 2.0 here). Sanity check."""
        sn_mesh = _build_sphere(n_cells=5, n_ord=4)
        ng = 1
        sigma_t_val = 2.0
        sigma_t = np.full((ng, sn_mesh.nx), sigma_t_val)
        psi_view = np.ones((sn_mesh.quad.N, ng, sn_mesh.nx))

        m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
        # At constant ψ = 1: (L+C)·1 ≈ σ_t · 1 = 2.0 everywhere.
        m_at_unknowns = _extract_at_unknown_slots(m_unified, sn_mesh)
        np.testing.assert_allclose(
            m_at_unknowns, sigma_t_val, rtol=1e-13, atol=1e-14,
        )

    def test_unified_zero_psi_gives_zero(self) -> None:
        """Linear operator: zero input → zero output."""
        sn_mesh = _build_sphere(n_cells=5, n_ord=4)
        ng = 1
        sigma_t = np.full((ng, sn_mesh.nx), 2.0)
        psi_view = np.zeros((sn_mesh.quad.N, ng, sn_mesh.nx))

        m_unified = legacy_proxy_matvec(psi_view, sn_mesh, sigma_t)
        np.testing.assert_array_equal(
            m_unified, np.zeros_like(m_unified),
        )
