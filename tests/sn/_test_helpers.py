"""Shared test helpers for the SN test suite.

Issue #197 PR-TYPED-0 introduced ``materials`` as a REQUIRED parameter
on :class:`SNMesh`.  Many geometry-only tests don't consume cross-
section values — they exercise sweep DAGs, BC realisation,
quadrature, cache structure, etc.  ``placeholder_materials`` provides
a minimal-but-valid :class:`Mixture` dict that those tests can hand
to :class:`SNMesh` so the construction succeeds without inviting any
real cross-section semantics into the test.

Issue #197 PR-TYPED-2 introduced :class:`BoundaryFlux` as the typed
replacement for the stringly-typed ``psi_bc: dict``.  Test fixtures
that previously passed ``{}`` to :func:`transport_sweep` should now
build a zero-initialised :class:`BoundaryFlux` via
``sn_mesh.zeros_boundary_flux()`` (or :func:`make_boundary_flux_zero`
below for non-SNMesh callers).

Tests that DO need realistic cross sections continue to use
``orpheus.derivations.common.xs_library.get_mixture`` etc. — this
helper is for the geometry-only call sites.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    from orpheus.sn.angular_flux import AngularFlux
    from orpheus.sn.boundary_flux import BoundaryFlux
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.scalar_flux import ScalarFlux


def placeholder_materials(
    ng: int = 1, mat_ids: tuple[int, ...] = (0,),
) -> dict:
    """Build a placeholder ``{mat_id: Mixture}`` dict for SNMesh tests.

    Parameters
    ----------
    ng
        Number of energy groups.  All Mixtures will report this value
        via :attr:`Mixture.ng`.
    mat_ids
        Material ids to include in the dict.  Default ``(0,)`` covers
        the common case where the mesh's ``mat_ids`` / ``mat_map`` is
        zeros.

    Returns
    -------
    dict[int, Mixture]
        Each entry has ``SigT = ones(ng)`` and all other cross sections
        zero.  Suitable for SNMesh tests that don't compute physical
        quantities from the materials.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return {
        int(mid): Mixture(
            SigC=z.copy(),
            SigL=z.copy(),
            SigF=z.copy(),
            SigP=z.copy(),
            SigT=np.ones(ng),
            SigS=[z_mat],
            Sig2=z_mat,
            chi=z.copy(),
        )
        for mid in mat_ids
    }


def make_boundary_flux_zero(sn_mesh: "SNMesh") -> "BoundaryFlux":
    """Build a zero-initialised :class:`BoundaryFlux` for ``sn_mesh``.

    Issue #197 PR-TYPED-2 — typed replacement for ``psi_bc = {}``.
    Allocates only the buffers the mesh's geometry consumes (slab gets
    two 1-D faces; curvilinear gets one outer face; 2-D Cartesian gets
    the persistent ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)``
    buffers).  Per-geometry dispatch lives inside
    :meth:`SNMesh.zeros_boundary_flux`; this helper is a clean alias
    so test fixtures don't have to chain through ``sn_mesh``.
    """
    return sn_mesh.zeros_boundary_flux()


def make_angular_flux_zero(sn_mesh: "SNMesh") -> "AngularFlux":
    """Build a zero-initialised :class:`AngularFlux` for ``sn_mesh``."""
    return sn_mesh.zeros_angular_flux()


def make_scalar_flux_zero(sn_mesh: "SNMesh") -> "ScalarFlux":
    """Build a zero-initialised :class:`ScalarFlux` for ``sn_mesh``."""
    return sn_mesh.zeros_scalar_flux()
