"""Foundation tests — SNMesh consumes ReducedStreamingOperator.

Round 1.1 of Wave D of the SN reshape campaign (Issue #159).

Pins three software invariants on the post-refactor :class:`SNMesh`:

1. ``self.reduced`` is a :class:`ReducedStreamingOperator` instance with
   the matching :class:`CoordSystem` for slab / sphere / cylinder.
2. The deprecated ``@property`` accessors (``alpha_half``,
   ``redist_dAw``, ``tau_mm``, ``alpha_per_level``,
   ``redist_dAw_per_level``, ``tau_mm_per_level``, ``face_areas``,
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
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _slab_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = GaussLegendre1D.create(4)
    return SNMesh(mesh, quad)


def _sphere_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = GaussLegendre1D.create(4)
    return SNMesh(mesh, quad)


def _cylinder_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = ProductQuadrature.create(n_mu=4, n_phi=4)
    return SNMesh(mesh, quad)


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
    # Sphere populates alpha_half / redist_dAw / tau_mm.
    assert sn.reduced.alpha_half is not None
    assert sn.reduced.redist_dAw is not None
    assert sn.reduced.tau_mm is not None


@pytest.mark.foundation
def test_cylinder_reduced_is_reduced_streaming_operator() -> None:
    sn = _cylinder_mesh()
    assert isinstance(sn.reduced, ReducedStreamingOperator)
    assert sn.reduced.coord is CoordSystem.CYLINDRICAL
    assert sn.reduced.requires_upstream_angular_state is True
    assert sn.reduced.angular_marching_axis == "mu"
    # Cylinder populates per-level structures.
    assert sn.reduced.alpha_per_level is not None
    assert sn.reduced.redist_dAw_per_level is not None
    assert sn.reduced.tau_mm_per_level is not None


# ---------------------------------------------------------------------------
# 2. Deprecated property access emits DeprecationWarning
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize(
    "attr",
    ["face_areas", "delta_A", "alpha_half", "redist_dAw", "tau_mm"],
)
def test_sphere_deprecated_properties_warn(attr: str) -> None:
    sn = _sphere_mesh()
    with pytest.warns(DeprecationWarning, match=f"SNMesh.{attr} is deprecated"):
        getattr(sn, attr)


@pytest.mark.foundation
@pytest.mark.parametrize(
    "attr",
    [
        "face_areas",
        "delta_A",
        "alpha_per_level",
        "redist_dAw_per_level",
        "tau_mm_per_level",
    ],
)
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
        assert sn.alpha_half is sn.reduced.alpha_half
        assert sn.redist_dAw is sn.reduced.redist_dAw
        assert sn.tau_mm is sn.reduced.tau_mm


@pytest.mark.foundation
def test_cylinder_deprecated_properties_route_to_reduced() -> None:
    sn = _cylinder_mesh()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert sn.face_areas is sn.reduced.face_areas
        assert sn.delta_A is sn.reduced.delta_A
        assert sn.alpha_per_level is sn.reduced.alpha_per_level
        assert sn.redist_dAw_per_level is sn.reduced.redist_dAw_per_level
        assert sn.tau_mm_per_level is sn.reduced.tau_mm_per_level


# ---------------------------------------------------------------------------
# 4. Slab keeps the Cartesian streaming stencils + has self.reduced
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_slab_keeps_cartesian_streaming_arrays() -> None:
    """``_setup_cartesian`` stays — DD-denominator arrays are SN-specific."""
    sn = _slab_mesh()
    # 2-D stencils populated for 1-D slab (ny=1 padding).
    assert sn.streaming_x is not None
    assert sn.streaming_y is not None
    # No curvature on Cartesian.
    assert sn.curvature is None
