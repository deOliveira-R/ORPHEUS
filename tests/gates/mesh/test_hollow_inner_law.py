"""A hollow curvilinear body's inner law is refused unless it is reflective (#511).

The radial axis (:class:`~orpheus.mesh.axis.RadialAxisMesh`) has one law
slot, the outer surface, so the declared inner law of a hollow cylinder or
sphere was dropped by :func:`~orpheus.mesh.axis.axes_from_legacy_mesh` and
every method built on the axes (SN, diffusion) computed a reflective cavity
instead, silently: before the guard, an inner vacuum and an inner reflective
law gave bit-identical SN fluxes on a cylinder and a sphere with
:math:`r_0 = 0.5` (``[M]`` 2026-09-25, the P1 specification's probe
``hollow_inner_law.py``). P1 step 2 makes a hollow body declarable on a
:class:`~orpheus.geometry.StructuredGeometry`, so the guard lands with it
(the specification's gate S3.9, moved into step 2 by the user's ruling of
2026-09-29).

The refusal is a declared scope boundary (``SCOPE-BOUNDARY[guard]``): the
machinery that would remove it is an inner-surface trace on the radial axis.
Reflective is admitted because it is what the methods compute, and that
is the void-cavity answer (the last gate). ``None`` is no law: since P1
step 3b a mesh refuses it on every face.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import solve_sn_fixed_source

pytestmark = pytest.mark.foundation

_CURVILINEAR = (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL)
_REFUSED_INNER_LAWS = {
    "vacuum": BC.vacuum,
    "white": BC("white"),
    "albedo": BC("albedo", {"albedo": 0.5}),
}


def _materials():
    return {0: get_mixture("B", "2g")}


def _quadrature(coord: CoordSystem) -> Quadrature:
    if coord is CoordSystem.CYLINDRICAL:
        return Quadrature.folded_product(n_mu=4, n_phi=8)
    return Quadrature.gauss_legendre(n_ordinates=8)


def _hollow_mesh(coord: CoordSystem, inner_law) -> Mesh1D:
    """A hollow body [0.5, 2.0] declared on a geometry, 8 cells."""
    geometry = StructuredGeometry(
        coord=coord,
        breakpoints=(0.5, 2.0),
        mat_ids=(0,),
        boundaries=(inner_law, BC.vacuum),
    )
    return Mesher(geometry).partition(CellsByCount.uniform_volume(8)).mesh


def _direct_hollow_mesh(coord: CoordSystem, inner_law) -> Mesh1D:
    """The same body built by the bare constructor, with no geometry."""
    edges = np.linspace(0.5, 2.0, 9)
    return Mesh1D(
        coord=coord, edges=edges, volumes=coord.measure(edges),
        mat_ids=np.zeros(8, int), face_laws={"xmin": inner_law, "xmax": BC.vacuum},
    )


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
@pytest.mark.parametrize("law", _REFUSED_INNER_LAWS, ids=str)
def test_sn_refuses_a_declared_inner_law_it_would_drop(coord, law):
    mesh = _hollow_mesh(coord, _REFUSED_INNER_LAWS[law])
    with pytest.raises(NotImplementedError, match="#511"):
        SNProblem(mesh, _quadrature(coord), _materials())


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
def test_a_direct_mesh_is_refused_too(coord):
    """The guard sits at the one adapter, so a mesh built without a
    geometry (the bare constructor) is covered."""
    mesh = _direct_hollow_mesh(coord, BC.vacuum)
    with pytest.raises(NotImplementedError, match="#511"):
        SNProblem(mesh, _quadrature(coord), _materials())


def test_diffusion_refuses_it_too():
    """Diffusion reads the same axes, so it dropped the same law."""
    mesh = _hollow_mesh(CoordSystem.SPHERICAL, BC.vacuum)
    with pytest.raises(NotImplementedError, match="#511"):
        DiffusionMesh(mesh, _materials())


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
def test_an_undeclared_inner_law_is_the_reflective_one(coord):
    """RE-POSED at P1 step 3b. Until then this was a RECORD that an
    undeclared (``None``) inner law and a declared reflective one gave the
    same flux bit for bit. ``None`` is no longer a law: the mesh refuses it,
    so the only way to spell that cavity is ``inner=BC.reflective``, whose
    answer the next gate verifies."""
    with pytest.raises(TypeError, match="None is not a boundary law"):
        _direct_hollow_mesh(coord, None)


def _core_mixture(sigma: float):
    """A pure absorber of total cross section ``sigma`` in both groups."""
    return make_mixture(
        sig_t=np.full(2, sigma), sig_c=np.full(2, sigma), sig_f=np.zeros(2),
        nu=np.zeros(2), chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def _shell_flux_with_core(coord: CoordSystem, core, n_core: int = 4) -> np.ndarray:
    """The shell's scalar flux of a SOLID body whose core [0, 0.5] holds ``core``."""
    quadrature = _quadrature(coord)
    geometry = StructuredGeometry(
        coord=coord, breakpoints=(0.0, 0.5, 2.0), mat_ids=(1, 0),
        boundaries=(BC.vacuum,),
    )
    mesh = Mesher(geometry).partition(
        (CellsByCount.uniform_volume(n_core), CellsByCount.uniform_volume(8)),
    ).mesh
    source = np.concatenate(
        [np.zeros((quadrature.N, 2, n_core)), np.ones((quadrature.N, 2, 8))], axis=2,
    )
    flux = np.asarray(solve_sn_fixed_source(
        {0: get_mixture("B", "2g"), 1: core}, mesh, quadrature, source,
    ).scalar_flux.values)
    return flux[:, n_core:]


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
def test_a_reflective_inner_law_is_a_void_cavity(coord):
    r"""The admitted law, verified as if derived: what SN computes for a
    hollow body with a reflective inner law is the transport answer for a
    VOID cavity.

    In 1-D spherical or cylindrical symmetry a ray entering a void core
    leaves it at the same impact parameter with the radial direction cosine
    reversed, so a void cavity reflects specularly. A solid body whose core
    is a pure absorber of cross section :math:`\varepsilon` therefore
    differs from the hollow reflective body by the core's absorption, a gap
    linear in :math:`\varepsilon`. The gate reads the gap at two values of
    :math:`\varepsilon` and extrapolates it to :math:`\varepsilon = 0`: the
    intercept is the discrete void-cavity answer's distance from the hollow
    reflective answer. Measured at this gate's own fixture (``[M]``
    2026-09-29: an 8-cell shell [0.5, 2.0], a 4-cell core, material B 2G,
    GL-8 on the sphere, the folded product 4 x 8 on the cylinder): gap / ε
    = 0.3263 (sphere) and 0.9983 (cylinder), the same at ε = 1e-6 and 1e-4;
    intercepts 2.0e-11 and 1.4e-10. On the sphere the gap per ε is also
    independent of the core's refinement (2 and 8 cells) and of the
    quadrature (GL-8, GL-16) in a 12-cell-shell probe, so the sweep across
    the void core adds nothing. Control: a black core (σ = 50) differs by
    0.360 (sphere) and 0.546 (cylinder).
    """
    hollow = np.asarray(solve_sn_fixed_source(
        _materials(), _hollow_mesh(coord, BC.reflective), _quadrature(coord),
        np.ones((_quadrature(coord).N, 2, 8)),
    ).scalar_flux.values)
    scale = np.max(np.abs(hollow))

    def gap(sigma: float) -> float:
        return float(np.max(np.abs(_shell_flux_with_core(coord, _core_mixture(sigma)) - hollow)) / scale)

    small, large = 1e-6, 1e-4
    gap_small, gap_large = gap(small), gap(large)
    slope = (gap_large - gap_small) / (large - small)
    intercept = gap_small - slope * small
    assert gap_large > 10.0 * gap_small, (
        f"the gap does not grow with the core's absorption "
        f"({gap_small:.3e} at {small}, {gap_large:.3e} at {large}): it is a floor, "
        f"and the reflective inner law is not the void cavity"
    )
    assert abs(intercept) < 1e-3 * gap_large, (
        f"the gap extrapolated to a void core is {intercept:.3e}, not zero "
        f"(gaps {gap_small:.3e}, {gap_large:.3e})"
    )
    black = float(np.max(np.abs(_shell_flux_with_core(coord, _core_mixture(50.0)) - hollow)) / scale)
    assert black > 0.1, f"control: a black core should differ at O(1), got {black:.3e}"
