r"""L1 cross-check, leg 3: one PRODUCTION SN sweep applies the algebra of record's face transmission.

Moved from ``tests/gates/transport/spatial/test_face_transmission_xverif.py`` by layer DAG Phase 4 (2026-11-01): legs 1,
1b and 2 there exercise the transport schemes' cell kernel and verdict; this leg drives one production sweep on a
one-cell mesh through :func:`orpheus.sn.loss_representation.default_for`, an SN path. The Branch-1 reference, the
tolerance and the helpers are that module's (imported, not copied); its docstring states the path table this leg
covers (diamond d = 1 ``CumprodScan``, d = 2 ``ScanMarch``, d = 3 the kernel wavefront; LD d = 1 ``CumprodScan``,
d = 2 the kernel) and why legs 1 and 3 are not redundant.

Markers: ``l1``; the diamond rows carry ``verifies("dd-face-transmission-spectrum")``; the d = 3 LD row is a strict
xfail (#503).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.derivations.discrete.sn.face_transmission import Closure
from orpheus.mesh import AxisMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import loss_representation
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from tests.gates._exterior import conditions_for
from tests.gates.transport.spatial.test_face_transmission_xverif import (
    _SCHEME,
    _SYM,
    LABEL,
    _assert_transmission,
    _fail_unless,
)

#: Exterior conditions of the tags this file's meshes carry.
_CONDITIONS = {0: WallCondition.vacuum()}


#: Two pure absorbers, so one sweep carries two different optical thicknesses
#: and a group-index slip (Mode 2) moves the answer.
_SIGMA_T = np.array([0.9, 3.7])


_WIDTHS = (0.7, 1.3, 0.4)


_AXES = "xyz"


def _one_cell_sweep_maps(scheme, d: int):
    """Per group: the octant map the production sweep applies, and its g."""
    mixture = make_mixture(
        sig_t=_SIGMA_T, sig_c=_SIGMA_T.copy(), sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )
    vacuum = 0
    axes = tuple(
        AxisMesh(edges=np.array([0.0, w]), role_low=vacuum, role_high=vacuum) for w in _WIDTHS[:d]
    )
    quadrature = Quadrature.level_symmetric(4) if d > 1 else Quadrature.gauss_legendre(n_ordinates=4)
    problem = SNProblem.from_axes(
        axes, quadrature, {0: mixture}, mat_map=np.zeros((1,) * d, dtype=int), scheme=scheme, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )
    representation = loss_representation.default_for(problem, problem.scheme, problem.angular_closure)
    q = problem.quad
    cosines = (q.mu_x, q.mu_y, q.mu_z)[:d]
    first_octant = np.logical_and.reduce([mu > 0 for mu in cosines])
    n = int(np.argmax(first_octant))
    _fail_unless(bool(first_octant[n]), f"d={d}: no first-octant ordinate")
    stratum = representation.bind_sigma(problem.sigma_t_cell)
    cell_moments = (2**d,) if isinstance(scheme, LinearDiscontinuous) else ()
    ng = _SIGMA_T.size
    per_face = AngularBoundaryFlux.zeros(problem.angular_trace).face_view("xmin")[n].size // ng
    maps = []
    for group in range(ng):
        columns = []
        for a in range(d):
            for moment in range(per_face):
                flux = AngularBoundaryFlux.zeros(problem.angular_trace)
                flux.face_view(f"{_AXES[a]}min")[n].reshape(ng, per_face)[group, moment] = 1.0
                representation.sweep(np.zeros((q.n_ordinates, ng) + (1,) * d + cell_moments), stratum, flux)
                columns.append(np.concatenate([
                    flux.face_view(f"{_AXES[b]}max")[n].reshape(ng, per_face)[group] for b in range(d)
                ]))
        g = [sp.Float(abs(float(cosines[a][n])) / (_WIDTHS[a] * _SIGMA_T[group]), 40) for a in range(d)]
        maps.append((np.array(columns).T, g))
    return maps, type(representation).__name__


_SWEEP_ROWS = [
    pytest.param(Closure.DIAMOND, 1, "CumprodScan", marks=pytest.mark.verifies(LABEL), id="diamond-1"),
    pytest.param(Closure.DIAMOND, 2, "ScanMarch", marks=pytest.mark.verifies(LABEL), id="diamond-2"),
    pytest.param(Closure.DIAMOND, 3, "FullFieldWavefront", marks=pytest.mark.verifies(LABEL), id="diamond-3"),
    pytest.param(Closure.LINEAR_DISCONTINUOUS, 1, "CumprodScan", id="ld-1"),
    pytest.param(Closure.LINEAR_DISCONTINUOUS, 2, "MovingFrontierWindow", id="ld-2"),
]


@pytest.mark.l1
@pytest.mark.rests_on(
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_diamond_kernel_applies_the_branch1_transmission",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_ld_kernel_applies_the_branch1_transmission",
    f"{_SYM}::test_diamond_transmission_is_the_page_closed_form",
)
@pytest.mark.parametrize(("closure", "d", "representation"), _SWEEP_ROWS)
def test_production_sweep_applies_the_branch1_transmission(closure, d, representation) -> None:
    r"""One sweep through a one-cell mesh applies Branch 1's octant map, per group.

    The representation is asserted because the row's coverage claim is about
    that path: if ``default_for`` changes its choice, this docstring and the
    module docstring's path table change with it.
    """
    maps, ran = _one_cell_sweep_maps(_SCHEME[closure](), d)
    _fail_unless(ran == representation, f"d={d}: default_for chose {ran}, the row covers {representation}")
    for group, (produced, g) in enumerate(maps):
        _assert_transmission(produced, closure, d, g, f"{closure.value} sweep d={d} group {group}")


@pytest.mark.l1
@pytest.mark.xfail(
    strict=True,
    raises=NotImplementedError,
    reason="#503: production LD sweeps at d=3 raise in _ubld.assemble_inflow_axis (interior axis)",
)
def test_production_ld_sweep_at_d3_applies_the_branch1_transmission() -> None:
    """The d = 3 LD sweep, red today for the named reason only (``raises=``)."""
    maps, _ = _one_cell_sweep_maps(LinearDiscontinuous(), 3)
    for group, (produced, g) in enumerate(maps):
        _assert_transmission(produced, Closure.LINEAR_DISCONTINUOUS, 3, g, f"LD sweep d=3 group {group}")
