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
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.sn.geometry import SNMesh
    from orpheus.transport.fields.scalar_flux import ScalarFlux


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


def legacy_proxy_matvec(
    psi_view: "np.ndarray", sn_mesh: "SNMesh", sigma_t: "np.ndarray",
    *, bc_outer=None, pole_angular_closure=None,
) -> "np.ndarray":
    """Call :func:`_transport_operator_matvec_unified` with the
    cell-centre-proxy boundary fill semantics (pre-B1'' convention).

    Tests that compare against L0 hand-derived references built BEFORE
    the B1'' face-aware architecture (i.e. references constructed
    around "no face state, fall back to cell-centre proxy" semantics)
    feed bare ``psi_view`` ``(N, ng, nx, ny)`` ndarrays and expect a
    bare ``(N, ng, nx, ny)`` cell-output ndarray.  This helper bridges:

    1. Build a :class:`BoundaryFlux` whose face buffers carry
       ``psi_view``'s cell-centre value at the outer (and slab-inner)
       face — the cell-centre-proxy fill.
    2. Wrap into a :class:`TimedFullField`.
    3. Call :func:`_transport_operator_matvec_unified`.
    4. Return ``result.bulk.values`` as a bare ndarray.

    The "legacy" prefix refers to the BOUNDARY-FILL CONVENTION (the
    pre-B1'' cell-centre proxy), not to retired code.  Production
    code uses the B1'' face-aware path via
    :class:`InvertibleOperator` (= ``L + C``); this helper exists only
    for L0 tests that pin the legacy convention's behaviour against
    closed-form hand references.
    """
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.timed_full_field import TimedFullField
    from orpheus.sn.operator import (
        CollisionOperator,
        StreamingOperator,
    )

    # Wave T T.5 close-out (matvec retirement): route through the
    # public operator-algebra path `(L + C).apply`.  The legacy
    # `_transport_operator_matvec_unified` helper was DELETED;
    # `_MSpatialOperatorSum._compute_decomposition` is the canonical
    # dual-emission body, and `(L + C).apply` = `L.apply + C.apply`
    # = `(M_spat + M_ang - sigma_t * psi) + sigma_t * psi` = (L+C)
    # bit-exact for the legacy semantic.  The ``bc_outer`` /
    # ``pole_angular_closure`` override parameters are not used by
    # any production caller of this helper today (all call sites pass
    # `bc_outer=None, pole_angular_closure=None`) — kept in the
    # function signature for legacy back-compat but ignored.
    del bc_outer, pole_angular_closure  # explicitly mark unused
    boundary = BoundaryFlux.zeros_for_sn_mesh(sn_mesh)
    boundary.face_view("xmax")[:] = psi_view[:, :, -1, 0]
    if "xmin" in boundary.layout.faces:
        boundary.face_view("xmin")[:] = psi_view[:, :, 0, 0]
    composite = TimedFullField(
        bulk=AngularFlux.from_mesh(psi_view, sn_mesh),
        boundary=boundary,
        _history=(),
        history_depth=2,
    )
    L_op = StreamingOperator(sn_mesh, sigma_t)
    C_op = CollisionOperator(sn_mesh, sigma_t)
    result = (L_op + C_op).apply(composite)
    return result.bulk.values


def _LC_matvec(
    psi: "TimedFullField", sigma_t: "np.ndarray",
) -> "TimedFullField":
    r"""Test-helper shim: returns ``(L + C).apply(psi)`` as a TimedFullField.

    Wave T T.5 close-out (matvec retirement, post-T.5.2): the module-
    level helper ``_transport_operator_matvec_unified`` was DELETED;
    its body lives as
    :meth:`_MSpatialOperatorSum._compute_decomposition`.  This shim
    constructs the canonical ``(L + C)`` operator-algebra composite
    and delegates to its public :meth:`apply` — the migration target
    for tests that previously called the deleted helper directly.

    Bit-identical to the legacy call ``_transport_operator_matvec_unified(psi, sigma_t)``
    for default ``bc_outer`` / ``pole_angular_closure`` (the only call
    convention any test actually exercised).
    """
    from orpheus.sn.operator import CollisionOperator, StreamingOperator
    sn_mesh = psi.bulk.mesh
    L = StreamingOperator(sn_mesh, sigma_t)
    C = CollisionOperator(sn_mesh, sigma_t)
    return (L + C).apply(psi)


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


def make_scalar_flux_zero(sn_mesh: "SNMesh") -> "ScalarFlux":
    """Build a zero-initialised :class:`ScalarFlux` for ``sn_mesh``."""
    return sn_mesh.zeros_scalar_flux()
