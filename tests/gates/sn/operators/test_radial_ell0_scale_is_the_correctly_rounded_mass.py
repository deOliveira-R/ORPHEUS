r"""**#708 / G2** -- the radial characteristic's constant-reproducing fold weight ``1/M`` reads the angular section's mass.

``RadialCharacteristicReconstruction._ell0_scale`` is the arc family's fold weight for the CONSTANT function, ``1/M``. ``M`` is the
mass of the angular axis, defined once, by the section ``AxisSectionOperator.of(angular_bulk_space, "angular").mass``: the correctly
rounded sum (``math.fsum``) of the stored weights (``ontology_w0_spec.md`` 5.7.2 item 6, E5, decision A1). Before the unit the scale
summed the weights itself (``1.0 / float(quadrature.weights.sum())``), a fifth spelling of the mass that differs from the section's by
one or two ULP on 5 of the 6 rules below.

The constructor of the unit takes the problem's axis-built angular space: ``RadialCharacteristicReconstruction(field_space, *,
coord, angular_bulk_space, n_moments=1)``. This file calls it in that spelling.

Kinds (``vv-principles``). Leg 1 is a REFERENCE (the platform-free ``1 / math.fsum`` of the stored weights, a different route from
the section's ``math.fsum(covector * unit)``); leg 2 is a THEOREM of the routing (halving under a doubled section mass is exact);
the sphere leg is a RECORD of the non-arc arm (``_ell0_scale is None``), the control that the change touches nothing there.

First red at the base tip ``954c4a4f`` (``-O``): every node fails at the constructor, ``TypeError: ... unexpected keyword argument
'angular_bulk_space'``. That is the committed red; the value-level red (the same inputs with the call spelled ``quadrature=sn.quad``:
``op._ell0_scale == 1/np.sum`` True and ``== 1/fsum`` False on the five rules where the sums differ, both True on ``(4, 8)``) is
recorded in the commit body from an uncommitted probe, so the gate is red for the signature and the value separately.

Activation. The rules ``(2, 12)``, ``(4, 12)``, ``(4, 16)``, ``(6, 12)`` are four of the 29 rules of the 120-rule population where
``np.sum != fsum``; ``(2, 6)`` is outside that population and also differs; ``(4, 8)`` is the control on which both agree.
The fixture asserts the difference (or the agreement) before reading, so a leg cannot pass vacuously.

Mutation (in-process, the mechanism of 1.3, applied by the reviewers): wrap ``RadialCharacteristicReconstruction.__init__`` so that
after the original body it assigns ``self._ell0_scale = 1.0 / float(np.asarray(sn.quad.weights).sum())`` when
``self._ell0_scale is not None`` (never on the sphere arm). Leg 1 reddens on the five differing rules and not on ``(4, 8)`` (which
is why the control exists); leg 2 reddens on all six, the control included; the sphere leg stays green.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.axis_marginal import AxisSectionOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.radial_characteristic import RadialCharacteristicReconstruction
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

pytestmark = pytest.mark.l0

_CONDITIONS = {0: WallCondition.specular(1.0)}
#: (folded_product arguments, whether ``np.sum`` and ``math.fsum`` of the stored weights differ)
_RULES = {(2, 12): True, (4, 12): True, (4, 16): True, (6, 12): True, (2, 6): True, (4, 8): False}


def _mixture():
    return make_mixture(
        sig_t=np.array([1.0]), sig_c=np.array([1.0]), sig_f=np.zeros(1), nu=np.zeros(1), chi=np.zeros(1), sig_s=np.array([[0.0]]),
    )


def _cylinder(rule: tuple[int, int]) -> SNProblem:
    geom = StructuredGeometry.cylinder((0.0, 2.0), (0,), outer=0)
    mesh = Mesher(geom).partition(CellsByCount.uniform_width(10)).mesh
    return SNProblem(
        mesh, Quadrature.folded_product(*rule), {0: _mixture()}, boundary_conditions=conditions_for(mesh, _CONDITIONS),
    )


def _sphere() -> SNProblem:
    geom = StructuredGeometry.sphere((0.0, 2.0), (0,), outer=0)
    mesh = Mesher(geom).partition(CellsByCount.uniform_width(10)).mesh
    return SNProblem(
        mesh, Quadrature.gauss_legendre(n_ordinates=8), {0: _mixture()}, boundary_conditions=conditions_for(mesh, _CONDITIONS),
    )


def _operator(sn: SNProblem) -> RadialCharacteristicReconstruction:
    return RadialCharacteristicReconstruction(
        sn.radial_characteristic_field_space, coord=sn.coord, angular_bulk_space=sn.angular_bulk_space,
    )


@pytest.mark.rests_on("tests/gates/numerics/test_axis_marginal_lift.py::test_g5_the_nodal_mass_is_correctly_rounded")
@pytest.mark.parametrize("rule", list(_RULES), ids=lambda r: f"folded{r}")
def test_the_scale_is_one_over_the_correctly_rounded_mass(rule):
    """``op._ell0_scale == 1 / math.fsum(stored weights)``, the platform-free reference (a private attribute: the one stored
    constant the operator folds with; no test read it before this one)."""
    sn = _cylinder(rule)
    weights = np.asarray(sn.quad.weights, dtype=float)
    assert (float(weights.sum()) != math.fsum(weights)) is _RULES[rule], "fixture activation: the sums differ exactly where declared"
    assert _operator(sn)._ell0_scale == 1.0 / math.fsum(weights)


@pytest.mark.rests_on("tests/gates/sn/operators/test_radial_ell0_scale_is_the_correctly_rounded_mass.py::test_the_scale_is_one_over_the_correctly_rounded_mass")
@pytest.mark.parametrize("rule", list(_RULES), ids=lambda r: f"folded{r}")
def test_the_scale_reads_the_sections_mass_object(rule, monkeypatch):
    """Doubling the section's mass halves the scale bit for bit (division by ``2M`` is division by ``M`` then by 2): a scale that
    sums the weights itself does not move. Red on the control rule ``(4, 8)`` too, which leg 1 cannot see."""
    sn = _cylinder(rule)
    section = AxisSectionOperator.of(sn.angular_bulk_space, "angular")
    before = _operator(sn)._ell0_scale
    assert before is not None and before > 0.0
    with monkeypatch.context() as m:
        m.setattr(section, "_mass", 2.0 * section.mass)  # a plain assignment inside the context is not tracked
        after = _operator(sn)._ell0_scale
    assert after == before / 2.0


@pytest.mark.rests_on("tests/gates/sn/operators/test_radial_ell0_scale_is_the_correctly_rounded_mass.py::test_the_scale_is_one_over_the_correctly_rounded_mass")
def test_the_sphere_arm_carries_no_scale():
    """The non-arc arm is untouched: ``_ell0_scale is None`` on a sphere (green before and after the change, but for the
    constructor's signature)."""
    assert _operator(_sphere())._ell0_scale is None
