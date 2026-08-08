"""Foundation tests — SNMesh consumes ReducedStreamingOperator.

Round 1.1 of Wave D of the SN reshape campaign (Issue #159).
Updated for Wave E Round 2 (Issue #164): the six curvature-specific
deprecated accessors (``alpha_half`` / ``redist_dAw`` / ``tau_mm``
sphere; ``alpha_per_level`` / ``redist_dAw_per_level`` /
``tau_mm_per_level`` cylinder) were retired along with the BiCGSTAB
FD-operator API surface they were the only consumer of.  Only
``face_areas`` and ``delta_A`` retain a transitional deprecation-
warning surface.

Pins three software invariants on the post-refactor :class:`SNMesh`:

1. ``self.reduced`` is a :class:`ReducedStreamingOperator` instance with
   the matching :class:`CoordSystem` for slab / sphere / cylinder.
2. The remaining transitional ``@property`` accessors (``face_areas``,
   ``delta_A``) emit :class:`DeprecationWarning` and route to the
   matching attribute on ``self.reduced``.
3. The values returned via the deprecated properties are the *same
   array objects* held by ``self.reduced`` — no copy, no recompute.

The bit-identical regression contract (11 frozen snapshots at
``tests/sn/regression/snapshots/``) is what gates the *math* of the
refactor; these tests gate the *contract* of the new
canonical accessor + transitional deprecation surface.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.geometry.reduced_operator import ReducedStreamingOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _slab_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = Quadrature.gauss_legendre(4)
    return SNMesh(mesh, quad, placeholder_materials())


def _sphere_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = Quadrature.gauss_legendre(4)
    return SNMesh(mesh, quad, placeholder_materials())


def _cylinder_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=4)
    return SNMesh(mesh, quad, placeholder_materials())


# ---------------------------------------------------------------------------
# 1. self.reduced is a ReducedStreamingOperator with matching coord
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_slab_reduced_is_reduced_streaming_operator() -> None:
    sn = _slab_mesh()
    assert isinstance(sn.reduced, ReducedStreamingOperator)
    assert sn.reduced.coord is CoordSystem.CARTESIAN
    # Slab carries no curvature math.
    assert sn.reduced.requires_upstream_angular_state is False
    assert sn.reduced.angular_marching_axis is None


@pytest.mark.foundation
def test_sphere_reduced_is_reduced_streaming_operator() -> None:
    sn = _sphere_mesh()
    assert isinstance(sn.reduced, ReducedStreamingOperator)
    assert sn.reduced.coord is CoordSystem.SPHERICAL
    assert sn.reduced.requires_upstream_angular_state is True
    assert sn.reduced.angular_marching_axis == "mu"
    # Sphere populates alpha_half / redist_dAw.  (Issue #236 Step C
    # retired the geometry-side tau_mm — the M-M angular weight is now
    # closure-owned.)
    assert sn.reduced.alpha_half is not None
    assert sn.reduced.redist_dAw is not None


@pytest.mark.foundation
def test_cylinder_reduced_is_reduced_streaming_operator() -> None:
    sn = _cylinder_mesh()
    assert isinstance(sn.reduced, ReducedStreamingOperator)
    assert sn.reduced.coord is CoordSystem.CYLINDRICAL
    assert sn.reduced.requires_upstream_angular_state is True
    assert sn.reduced.angular_marching_axis == "mu"
    # Cylinder populates per-level structures.  (Issue #236 Step C
    # retired the geometry-side tau_mm_per_level — the M-M angular
    # weight is now closure-owned.)
    assert sn.reduced.alpha_per_level is not None
    assert sn.reduced.redist_dAw_per_level is not None


# ---------------------------------------------------------------------------
# 2. Deprecated property access emits DeprecationWarning
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("attr", ["face_areas", "delta_A"])
def test_sphere_deprecated_properties_warn(attr: str) -> None:
    sn = _sphere_mesh()
    with pytest.warns(DeprecationWarning, match=f"SNMesh.{attr} is deprecated"):
        getattr(sn, attr)


@pytest.mark.foundation
@pytest.mark.parametrize("attr", ["face_areas", "delta_A"])
def test_cylinder_deprecated_properties_warn(attr: str) -> None:
    sn = _cylinder_mesh()
    with pytest.warns(DeprecationWarning, match=f"SNMesh.{attr} is deprecated"):
        getattr(sn, attr)


# ---------------------------------------------------------------------------
# 3. Deprecated properties return the same arrays as self.reduced
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_sphere_deprecated_properties_route_to_reduced() -> None:
    sn = _sphere_mesh()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        # Identity check, not just equality — the property returns the
        # exact array held by ``self.reduced``, not a copy.
        assert sn.face_areas is sn.reduced.face_areas
        assert sn.delta_A is sn.reduced.delta_A


@pytest.mark.foundation
def test_cylinder_deprecated_properties_route_to_reduced() -> None:
    sn = _cylinder_mesh()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert sn.face_areas is sn.reduced.face_areas
        assert sn.delta_A is sn.reduced.delta_A


# ---------------------------------------------------------------------------
# 4. Slab keeps the Cartesian streaming stencils + has self.reduced
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_slab_keeps_cartesian_streaming_arrays() -> None:
    """``_setup_cartesian`` stays — DD-denominator arrays are SN-specific.

    HISTORY (C3.6): the slab carries exactly ONE per-axis streaming array
    — the phantom ``ny=1`` ``streaming_y`` padding is gone (the stencil is
    built over ``range(ndim)``); ``streaming(1)`` on a 1-D mesh is an
    IndexError, not a silently-shaped phantom.
    """
    sn = _slab_mesh()
    assert sn.streaming(0) is not None
    with pytest.raises(IndexError, match="out of range for ndim=1"):
        sn.streaming(1)
    # No curvature on Cartesian.
    assert sn.curvature is None


@pytest.mark.foundation
def test_curvilinear_streaming_accessor_raises() -> None:
    """G-b4 (C3.6): ``streaming(axis)`` is Cartesian-only — curvilinear
    meshes carry streaming in ``reduced.streaming_terms`` (the chain-scan
    substrate), never the per-axis DD stencil."""
    for sn in (_sphere_mesh(), _cylinder_mesh()):
        with pytest.raises(AttributeError, match="Cartesian-only"):
            sn.streaming(0)


# ---------------------------------------------------------------------------
# 5. dag_walk — DAG-topological traversal
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_sphere_dag_walk_outward_order() -> None:
    """Outward sweep (μ ≥ 0): cell 0 → nx-1, downstream = outer face."""
    sn = _sphere_mesh()
    quad = sn.quad
    # Pick a positive-μ ordinate.
    n = quad.N - 1
    assert quad.mu_x[n] > 0
    visits = list(sn.dag_walk(ordinate_idx=n))
    assert [v.cell_idx for v in visits] == list(range(sn.nx))
    # Each visit's downstream face is the OUTER face of the cell.
    for v in visits:
        assert v.face_area_downstream == v.streaming_terms.face_area_outer
        # Geometric labels are direction-independent — anchor.
        assert v.streaming_terms.face_area_inner < v.streaming_terms.face_area_outer


@pytest.mark.foundation
def test_sphere_dag_walk_inward_order() -> None:
    """Inward sweep (μ < 0): cell nx-1 → 0, downstream = inner face."""
    sn = _sphere_mesh()
    quad = sn.quad
    # Pick a negative-μ ordinate.
    n = 0
    assert quad.mu_x[n] < 0
    visits = list(sn.dag_walk(ordinate_idx=n))
    assert [v.cell_idx for v in visits] == list(range(sn.nx - 1, -1, -1))
    # Each visit's downstream face is the INNER face.
    for v in visits:
        assert v.face_area_downstream == v.streaming_terms.face_area_inner


@pytest.mark.foundation
def test_slab_dag_walk_no_face_areas() -> None:
    """Slab: forward / backward iteration with face_area_downstream = 1.0.

    Issue #196 Step 2.5: slab carries ``face_area_downstream = 1.0``
    (neutral curvature) so the unified cell-balance helper consumes
    one geometry-blind number.  Pre-Step-2.5 slab carried ``None``.
    """
    sn = _slab_mesh()
    quad = sn.quad
    # Forward (μ > 0).
    n_pos = quad.N - 1
    visits_pos = list(sn.dag_walk(ordinate_idx=n_pos))
    assert [v.cell_idx for v in visits_pos] == list(range(sn.nx))
    for v in visits_pos:
        assert v.face_area_downstream == 1.0
    # Backward (μ < 0).
    n_neg = 0
    visits_neg = list(sn.dag_walk(ordinate_idx=n_neg))
    assert [v.cell_idx for v in visits_neg] == list(range(sn.nx - 1, -1, -1))


@pytest.mark.foundation
def test_cylinder_dag_walk_per_level() -> None:
    """Cylindrical: per-level traversal; downstream face follows η sign."""
    sn = _cylinder_mesh()
    quad = sn.quad
    # Pick a level with both signs of η represented.
    for p, level_idx in enumerate(quad.level_indices):
        eta_signs = [np.sign(quad.eta[g]) for g in level_idx]
        if 1 in eta_signs and -1 in eta_signs:
            chosen_level = p
            break
    else:  # pragma: no cover — safety
        pytest.skip("No level has both signs of eta")

    level_idx = quad.level_indices[chosen_level]
    for m_local in range(len(level_idx)):
        global_n = int(level_idx[m_local])
        eta_n = quad.eta[global_n]
        visits = list(
            sn.dag_walk(
                ordinate_idx=m_local, mu_level_idx=chosen_level,
            )
        )
        # Each visit carries the GLOBAL ordinate's |η| in abs_mu —
        # bug 2 fix anchor.
        for v in visits:
            assert v.streaming_terms.abs_mu == abs(eta_n)
        if abs(eta_n) < 1e-15:
            # Degenerate: no spatial flow, forward order.  Issue #196
            # Step 2.5: ``face_area_downstream == 0.0`` (geometric
            # truth) replaces the ``None`` sentinel.
            assert all(v.face_area_downstream == 0.0 for v in visits)
            assert [v.cell_idx for v in visits] == list(range(sn.nx))
        elif eta_n >= 0:
            # Outward: forward, downstream = outer.
            assert [v.cell_idx for v in visits] == list(range(sn.nx))
            for v in visits:
                assert v.face_area_downstream == (
                    v.streaming_terms.face_area_outer
                )
        else:
            # Inward: backward, downstream = inner.
            assert [v.cell_idx for v in visits] == list(
                range(sn.nx - 1, -1, -1)
            )
            for v in visits:
                assert v.face_area_downstream == (
                    v.streaming_terms.face_area_inner
                )


@pytest.mark.foundation
def test_cylinder_dag_walk_requires_level_idx() -> None:
    """Cylindrical dag_walk without mu_level_idx raises."""
    sn = _cylinder_mesh()
    with pytest.raises(ValueError, match="mu_level_idx"):
        list(sn.dag_walk(ordinate_idx=0))
