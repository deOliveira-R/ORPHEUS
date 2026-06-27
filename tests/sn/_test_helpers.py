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
``BoundaryFlux.zeros_on(sn_mesh)`` (or :func:`make_boundary_flux_zero`
below for non-SNMesh callers).

Tests that DO need realistic cross sections continue to use
``orpheus.derivations.common.xs_library.get_mixture`` etc. — this
helper is for the geometry-only call sites.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux

if TYPE_CHECKING:
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields.scalar_flux import ScalarFlux


# Anchor for shared, sn-root-relative test data (regression snapshots,
# the sweep reference vector, the Wave-T fixture .npz files).  The
# capability-taxonomy reorg nests tests several directories deep
# (e.g. ``sweep/curvilinear/``); a ``Path(__file__).parent``-relative
# data lookup would break on every move.  Every consumer resolves data
# through this single anchor so the data store stays at the sn-root
# regardless of where the test that reads it lives.
SN_TESTS_ROOT = Path(__file__).resolve().parent
"""Absolute path to ``tests/sn/`` — the anchor for shared test data."""


def stamp_capability_marker(items, conftest_file: str, capability: str) -> None:
    """Apply ``@pytest.mark.cap(<capability>)`` to every test under a dir.

    The capability-taxonomy reorg encodes the SN-capability tier of a
    test as the *directory* it lives in (single source of truth). Rather
    than decorate every test file with a ``cap(...)`` marker — which can
    drift from the directory it documents — each capability directory
    carries a one-line ``conftest.py`` that delegates here. Every item
    collected at or below the conftest's directory gets the marker; the
    existing ``l0/l1/l2/foundation/verifies/catches`` markers on each
    test are untouched (``cap`` is orthogonal and composable).

    Parameters
    ----------
    items
        The collected items passed to ``pytest_collection_modifyitems``.
    conftest_file
        ``__file__`` of the calling capability-directory conftest.
    capability
        The capability name (one of the DAG nodes; see the ``cap``
        marker description in ``pyproject.toml``).
    """
    import pytest

    here = Path(conftest_file).resolve().parent
    marker = pytest.mark.cap(capability)
    for item in items:
        try:
            item_path = Path(str(item.fspath)).resolve()
        except Exception:
            continue
        # ``here`` is the item's parent OR an ancestor — covers both a
        # flat capability dir and a nested one (e.g. sweep/core under
        # sweep). The sweep/ conftest does NOT stamp because each leaf
        # (core, slab, curvilinear, cartesian_2d) carries its own.
        if item_path.parent == here:
            item.add_marker(marker)



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


# ── Shared curvilinear / slab mesh builders ──────────────────────────
#
# Lifted verbatim from the legacy ``test_cylindrical.py`` /
# ``test_spherical.py`` modules during the SN taxonomy reorg so the
# cylinder + sphere split files (eigenvalue/, sweep/curvilinear/,
# verification/analytical/) share ONE definition rather than each
# carrying a copy.

_COORD_TO_TAG = {
    "CARTESIAN": "SLB",
    "CYLINDRICAL": "CYL",
    "SPHERICAL": "SPH",
}


def _bcs_for(tag: str, bc):
    """BC tuple matching the geometry tag's endpoint count."""
    if tag == "SLB":
        return (bc, bc)
    return (bc,)


def curvilinear_homogeneous_mesh(
    n_cells: int,
    total_width: float,
    mat_id: int = 0,
    coord=None,
    bc=None,
):
    """Single-region uniform mesh in any coordinate system.

    SN tests default to ``BC.reflective`` (the eigenvalue / lattice
    convention). CP tests must override to ``BC.white`` because CP
    only supports ``"vacuum"`` / ``"white"``.
    """
    from orpheus.geometry import (
        BC, CoordSystem, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    if coord is None:
        coord = CoordSystem.CARTESIAN
    if bc is None:
        bc = BC.reflective
    tag = _COORD_TO_TAG[coord.name]
    geom = StructuredGeometry(
        geometry=tag,
        regions=(Region(mat_id=mat_id, outer_thickness_cm=total_width),),
        bcs=_bcs_for(tag, bc),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def curvilinear_two_region_mesh(
    outers,
    mat_ids,
    n_cells,
    coord,
    bc=None,
):
    """Two-region mesh with absolute outer-edge convention."""
    from orpheus.geometry import (
        BC, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.geometry import Mesh1D
    if bc is None:
        bc = BC.reflective
    tag = _COORD_TO_TAG[coord.name]
    geom = StructuredGeometry(
        geometry=tag,
        regions=(
            Region(mat_id=mat_ids[0], outer_thickness_cm=outers[0]),
            Region(mat_id=mat_ids[1], outer_thickness_cm=outers[1] - outers[0]),
        ),
        bcs=_bcs_for(tag, bc),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=n_cells[0]),
        RegionMesh(n_cells=n_cells[1]),
    ))


def cart2d_2g_nonsquare(nx: int = 5, ny: int = 7) -> "SNMesh":
    """2-D Cartesian, reflective, 2G, NON-SQUARE (the x↔y-swap moat).

    The discriminating config for structural operator/representation
    tests: the octant frame, the pure-z branch, and the interior
    recurrence all degenerate on a 1G flat square box (vv §H1/§H2).
    Promoted from ``test_one_octant_walk.py`` when the S6.5
    one-instance tests became its second consumer.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh2D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.mesh.augmented_mesh import SNMesh

    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2))


def het_operands(sn: "SNMesh"):
    """Heterogeneous σ_t + a non-flat random state (≥2G, non-degenerate).

    Returns ``(sig_t, psi)`` where ``sig_t`` is a random per-group
    per-cell total cross section and ``psi`` is a
    :class:`~orpheus.transport.timed_full_field.TimedFullField` with
    random bulk AND boundary values — every term of the loss action is
    activated (nothing nulled by a flat or zero state).
    """
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.timed_full_field import TimedFullField

    rng = np.random.default_rng(20260611)
    sig_t = rng.uniform(0.3, 3.0, size=(sn.ng, *sn.spatial_shape))
    psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)
    psi.bulk.values[...] = rng.standard_normal(psi.bulk.values.shape)
    for face in psi.boundary.layout.faces:
        fv = psi.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)
    return sig_t, psi


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
    from orpheus.sn.operators.streaming import (
        StreamingOperator,
    )
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

    # Wave T T.5 close-out (matvec retirement): route through the
    # public operator-algebra path `(L + C).apply`.  The legacy
    # `_transport_operator_matvec_unified` helper was DELETED; the
    # canonical 1-D matvec body is now `_OneDimScanWalk._apply_walk`
    # (the fused `(L+C)ψ`), and `(L + C).apply` = `L.apply + C.apply`
    # = `((L+C)ψ - sigma_t * psi) + sigma_t * psi` = (L+C)
    # bit-exact for the legacy semantic.  The ``bc_outer`` /
    # ``pole_angular_closure`` override parameters are not used by
    # any production caller of this helper today (all call sites pass
    # `bc_outer=None, pole_angular_closure=None`) — kept in the
    # function signature for legacy back-compat but ignored.
    del bc_outer, pole_angular_closure  # explicitly mark unused
    boundary = BoundaryFlux.zeros_on(sn_mesh)
    boundary.face_view("xmax")[:] = psi_view[:, :, -1]
    if "xmin" in boundary.layout.faces:
        boundary.face_view("xmin")[:] = psi_view[:, :, 0]
    composite = TimedFullField(
        bulk=AngularFlux.from_mesh(psi_view, sn_mesh),
        boundary=boundary,
        _history=(),
        history_depth=2,
    )
    L_op = StreamingOperator(sn_mesh)
    C_op = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
    result = (L_op + C_op).apply(composite)
    return result.bulk.values


def _LC_matvec(
    psi: "TimedFullField", sigma_t: "np.ndarray",
) -> "TimedFullField":
    r"""Test-helper shim: returns ``(L + C).apply(psi)`` as a TimedFullField.

    Wave T T.5 close-out (matvec retirement, post-T.5.2): the module-
    level helper ``_transport_operator_matvec_unified`` was DELETED;
    the canonical 1-D matvec body is now
    :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
    (the fused ``(L+C)ψ``).  This shim constructs the canonical
    ``(L + C)`` operator-algebra composite and delegates to its public
    :meth:`apply` — the migration target for tests that previously
    called the deleted helper directly.

    Bit-identical to the legacy call ``_transport_operator_matvec_unified(psi, sigma_t)``
    for default ``bc_outer`` / ``pole_angular_closure`` (the only call
    convention any test actually exercised).
    """
    from orpheus.sn.operators.streaming import StreamingOperator
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
    sn_mesh = psi.bulk.mesh
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
    return (L + C).apply(psi)


def make_boundary_flux_zero(sn_mesh: "SNMesh") -> "BoundaryFlux":
    """Build a zero-initialised :class:`BoundaryFlux` for ``sn_mesh``.

    Issue #197 PR-TYPED-2 — typed replacement for ``psi_bc = {}``.
    Allocates only the buffers the mesh's geometry consumes (slab gets
    two 1-D faces; curvilinear gets one outer face; 2-D Cartesian gets
    the persistent ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)``
    buffers).  Per-geometry dispatch lives inside
    :meth:`~orpheus.transport.fields.boundary_flux.BoundaryFlux.zeros_on`; this helper is a clean alias
    so test fixtures don't have to chain through ``sn_mesh``.
    """
    return BoundaryFlux.zeros_on(sn_mesh)


def make_scalar_flux_zero(sn_mesh: "SNMesh") -> "ScalarFlux":
    """Build a zero-initialised :class:`ScalarFlux` for ``sn_mesh``."""
    return ScalarFlux.zeros_on(sn_mesh)
