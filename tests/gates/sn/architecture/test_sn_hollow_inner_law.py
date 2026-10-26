"""A hollow curvilinear body's inner condition is refused unless it is the unit specular wall (#511).

The radial axis (:class:`~orpheus.mesh.axis.RadialAxisMesh`) has one trace slot, the outer surface, so the declared inner
condition of a hollow cylinder or sphere has no operator to land in, and every method built on the axes (SN, diffusion) computes
a reflective cavity instead, silently (until the wave-2 carve the geometry adapter refused it with ``NotImplementedError``; the
inner face is now a TAG of the carrier, with a condition bound to it, and the refusal is the realizer table's: SN's table admits
the inner site only at a unit specular return, ``Refused(PLUMBING, #511)`` otherwise): before the guard, an inner vacuum and an inner reflective
law gave bit-identical SN fluxes on a cylinder and a sphere with
:math:`r_0 = 0.5` (``[M]`` 2026-09-25, the P1 specification's probe
``hollow_inner_law.py``). P1 step 2 makes a hollow body declarable on a
:class:`~orpheus.geometry.StructuredGeometry`, so the guard lands with it
(the specification's gate S3.9, moved into step 2 by the user's ruling of
2026-09-29).

The refusal is a declared scope boundary (``SCOPE-BOUNDARY[guard]``): the
machinery that would remove it is an inner-surface trace on the radial axis.
The perfect mirror, ``WallCondition.specular(1.0)`` (rule R2: a curved face has no mirror deck), is admitted because it is what the
methods compute, and that is the void-cavity answer (the last gate). ``None`` is no role: a mesh refuses it on every face.

Moved from ``tests/gates/mesh/test_hollow_inner_law.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26): the refusal is
SN's, so the rows live with SN's gates; the bodies are unchanged. The hollow body, its tags and the refused conditions shared
with diffusion are in ``tests/gates/transport/_hollow_inner_helpers.py``; the ``None``-role refusal of the mesh stayed in
``tests/gates/mesh/test_hollow_inner_law.py``.

Diffusion's refusal of the same condition is ``tests/gates/diffusion/test_diffusion_hollow_inner_law.py``."""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture, make_mixture
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, WallCondition
from tests.gates.transport._hollow_inner_helpers import (
    _CURVILINEAR,
    _INNER,
    _OUTER,
    _REFUSED_INNER_LAWS,
    _conditions,
    _hollow_mesh,
    _materials,
)

pytestmark = pytest.mark.foundation


def _quadrature(coord: CoordSystem) -> Quadrature:
    if coord is CoordSystem.CYLINDRICAL:
        return Quadrature.folded_product(n_mu=4, n_phi=8)
    return Quadrature.gauss_legendre(n_ordinates=8)


def _direct_hollow_mesh(coord: CoordSystem, inner_role) -> Mesh1D:
    """The same body built by the bare constructor, with no geometry."""
    edges = np.linspace(0.5, 2.0, 9)
    return Mesh1D(
        coord=coord, edges=edges, volumes=coord.measure(edges),
        region_ids=np.zeros(8, int), region_materials=(0,), face_roles={"xmin": inner_role, "xmax": _OUTER},
    )


@pytest.mark.catches("ERR-120")
@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
@pytest.mark.parametrize("law", _REFUSED_INNER_LAWS, ids=str)
def test_sn_refuses_a_declared_inner_law_it_would_drop(coord, law):
    mesh = _hollow_mesh(coord)
    with pytest.raises(BoundaryRefused) as info:
        SNProblem(mesh, _quadrature(coord), _materials(), boundary_conditions=_conditions(_REFUSED_INNER_LAWS[law]))
    assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 511 and info.value.face == "xmin"


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
def test_a_direct_mesh_is_refused_too(coord):
    """The refusal is the realizer table's, read at the binding of the carrier's faces, so a mesh built without a
    geometry (the bare constructor) is covered."""
    mesh = _direct_hollow_mesh(coord, _INNER)
    with pytest.raises(BoundaryRefused) as info:
        SNProblem(mesh, _quadrature(coord), _materials(), boundary_conditions=_conditions(WallCondition.vacuum()))
    assert info.value.admission.issue == 511


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
        face_roles=(_OUTER,),
    )
    mesh = Mesher(geometry).partition(
        (CellsByCount.uniform_volume(n_core), CellsByCount.uniform_volume(8)),
    ).mesh
    source = np.concatenate(
        [np.zeros((quadrature.n_ordinates, 2, n_core)), np.ones((quadrature.n_ordinates, 2, 8))], axis=2,
    )
    flux = np.asarray(solve_sn_fixed_source(
        {0: get_mixture("B", "2g"), 1: core}, mesh, quadrature, source, boundary_conditions={_OUTER: WallCondition.vacuum()},
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
        _materials(), _hollow_mesh(coord), _quadrature(coord),
        np.ones((_quadrature(coord).n_ordinates, 2, 8)), boundary_conditions=_conditions(WallCondition.specular(1.0)),
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
