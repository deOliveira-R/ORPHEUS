r"""**#708 / G1** -- the angular mass has ONE definition, the angular section's, and every reader of the producer-side ``1/W`` reads it.

The law (``ontology_w0_spec.md`` 5.7.4). The angular axis carries the discrete measure of the STORED weights ``w_n``; the section
``E`` is ``phi -> (u / M) phi`` and ``R E = id`` iff ``M = sum_n w_n`` exactly. ``M`` is therefore the correctly rounded sum of the
stored weights (``math.fsum``), which does not depend on the order of the nodes; ``float(np.sum(w))`` is a pairwise sum in storage
order and is not that value on 29 of the 120 rules of ``test_axis_marginal_lift._RULES`` (GL7 and GL12 give ``1.9999999999999998``
where ``fsum`` gives ``2.0``). The law of the unit: ``lift.total_weight == AxisSectionOperator.of(lift._codomain_interior,
"angular").mass == math.fsum(stored weights)`` bit for bit, invariant under a permutation of the stored weights, and every reader of
``total_weight`` divides by the section's mass object, not by a sum of its own.

Kinds (``vv-principles``). Test 1 is a THEOREM (a law of the stored weights, against the platform-free ``math.fsum``); test 2 is a THEOREM
(order invariance of the mass); test 3 is a THEOREM of the routing (each reader's output scales exactly with the section's mass
object, a pure power-of-two scaling, so a reader that sums the weights itself does not move). No row is a RECORD.

First reds at the base tip ``954c4a4f`` (command and failing lines are in the commit body; ``-O``):

* ``test_every_lift_divides_by_the_correctly_rounded_mass[7|12]``: ``assert 1.9999999999999998 == 2.0``.
* ``test_the_divisor_does_not_depend_on_the_enumeration_order``: ``{'F': 1.9999999999999998, 'S_L1': 1.9999999999999998}`` against
  ``{'F': 2.0, 'S_L1': 2.0}``.
* ``test_every_reader_reads_the_sections_mass_object``: every one of the 6 (lift, reader) legs reads ``False`` (none divides by
  the section's mass) on GL12; each leg is its own node, so the red set names the reader.

Activation. GL7 and GL12 are rules on which ``np.sum != fsum``, so the legs can fail; GL4 and GL8 are the controls on which both agree
(``2.0``) and the first test passes today. The permuted enumeration is the stable argsort of the weights.

Mutations (the mechanism of 1.3: in-process ``monkeypatch``; the reviewers apply them, this file does not):

(a) ``monkeypatch.setattr(AngularLift, "total_weight", property(lambda self: float(np.asarray(self.frame.measure.weights).sum())))``
    with ``_combine`` already on the section: tests 1 and 2 redden on both lifts, and the reader legs redden on ``_interior_transpose``
    (the ``apply_transpose`` legs) and on the two transfer redistributions.
(b) ``monkeypatch.setattr(TransferOperator, "_redistribute_ordinates", f)`` where ``f`` is a copy of its body dividing by
    ``float(np.asarray(self.frame.measure.weights).sum())``: only the ``redistribute_ordinates`` leg of ``S_L1`` reddens (the red set
    is the nodes of that reader, which is why the readers are separate parametrised nodes and not a loop).

The retirement half of the unit (the mixed sum, ``CONTAINS``) is G3 and ``tests/gates/_retired/rows_w0_angular_mass.py``; DSA is not
a leg (it stores ``self._section.mass`` already; a gate on it is green before and after, X1).
"""
from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.axis_marginal import AxisSectionOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields import AngularSourceSink, ScalarSourceSink
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from tests.gates.transport._phase_space_cases import angular_full_field_space, posing_quadrature, posing_scheme

pytestmark = pytest.mark.foundation

_RULES_WHERE_THE_SUMS_DIFFER = (7, 12)  # float(np.sum(w)) != math.fsum(w): 1.9999999999999998 against 2.0
_RULES_WHERE_THE_SUMS_AGREE = (4, 8)  # both 2.0: the controls, green before and after


def _mesh() -> MaterialMesh:
    """One slab, 20 cells, the two-group mixture ``A``: the fixture of ``tests/gates/transport/test_angular_lift.py``."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(CellsByCount.uniform_width(20)).mesh
    return MaterialMesh(mesh, {0: get_mixture("A", "2g")})


def _permuted(quadrature: Quadrature) -> Quadrature:
    """The same rule enumerated in order of increasing weight (nodes and weights permuted together)."""
    order = np.argsort(np.asarray(quadrature.weights), kind="stable")
    measure = dataclasses.replace(
        quadrature.measure, nodes=np.asarray(quadrature.nodes)[order], weights=np.asarray(quadrature.weights)[order],
    )
    return Quadrature(measure=measure)


def _posed(n: int, *, permuted: bool = False):
    """``(mesh, quadrature, space, lifts)``: the lifts are ``F`` (fission) and ``S_L1`` (scattering at order 1)."""
    mesh = _mesh()
    quadrature = posing_quadrature(mesh, n)
    if permuted:
        quadrature = _permuted(quadrature)
    space = angular_full_field_space(mesh, quadrature, posing_scheme(mesh))
    lifts = {
        "F": FissionOperator.from_solver_data(mat_xs=mesh.mat_xs, space=space),
        "S_L1": ScatteringOperator.from_solver_data(mat_xs=mesh.mat_xs, scattering_order=1, space=space),
    }
    return mesh, quadrature, space, lifts


def _mass_of(lift) -> float:
    return AxisSectionOperator.of(lift._codomain_interior, "angular").mass


@pytest.mark.rests_on("tests/gates/numerics/test_axis_marginal_lift.py::test_g5_the_nodal_mass_is_correctly_rounded")
@pytest.mark.parametrize("n", _RULES_WHERE_THE_SUMS_DIFFER)
def test_every_lift_divides_by_the_correctly_rounded_mass(n):
    """``total_weight == section.mass == math.fsum(w)`` on the two rules where ``np.sum`` and ``fsum`` differ.

    Red today: ``F.total_weight`` is ``1.9999999999999998`` and the section's mass ``2.0``. The first assertion is the fixture's
    activation: on a rule where the two sums agree the leg would pass vacuously."""
    _, quadrature, _, lifts = _posed(n)
    weights = np.asarray(quadrature.weights, dtype=float)
    assert float(weights.sum()) != math.fsum(weights), "fixture must be a rule on which the two sums differ"
    for key, lift in lifts.items():
        assert lift.total_weight == _mass_of(lift) == math.fsum(weights), key


@pytest.mark.rests_on("tests/gates/numerics/test_axis_marginal_lift.py::test_g5_the_nodal_mass_is_correctly_rounded")
@pytest.mark.parametrize("n", _RULES_WHERE_THE_SUMS_AGREE)
def test_the_control_rules_agree_in_both_spellings(n):
    """Control: on GL4 and GL8 the two spellings agree (``2.0``), so this leg is green before and after the change."""
    _, quadrature, _, lifts = _posed(n)
    weights = np.asarray(quadrature.weights, dtype=float)
    assert float(weights.sum()) == math.fsum(weights)
    for key, lift in lifts.items():
        assert lift.total_weight == math.fsum(weights), key


@pytest.mark.rests_on("tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py::test_every_lift_divides_by_the_correctly_rounded_mass")
def test_the_divisor_does_not_depend_on_the_enumeration_order():
    """The divisor on GL7 equals the divisor on the same rule enumerated in sorted order: a mass that depends on the enumeration
    is not the measure's mass. Red today: ``{'F': 1.9999999999999998, 'S_L1': 1.9999999999999998}`` against ``2.0`` for both."""
    natural = {k: lift.total_weight for k, lift in _posed(7)[3].items()}
    sorted_order = {k: lift.total_weight for k, lift in _posed(7, permuted=True)[3].items()}
    assert natural == sorted_order


# ── the readers: one node per (rule, lift, reader) ───────────────────────────────────────────────────────────────────────────


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def _read_combine(lift, mesh, space):
    iso_space = lift._scalar_interior_space
    iso = ScalarSourceSink(values=_rng(1).standard_normal(iso_space.shape), space=iso_space)
    return np.asarray(lift._combine(iso, None).values)


def _read_transpose(lift, mesh, space):
    interior = AngularSourceSink(
        values=_rng(2).standard_normal(lift._codomain_interior.shape), space=lift._codomain_interior,
    )
    cotangent = FullField(interior=interior, boundary=AngularBoundaryFlux.zeros(space.trace_space))
    return np.asarray(lift.apply_transpose(cotangent).interior.values)


def _read_redistribute_ordinates(lift, mesh, space):
    domain = lift._domain_interior
    flux = AngularFlux(values=_rng(3).standard_normal(domain.shape), space=domain)
    return np.asarray(lift._redistribute_ordinates(flux).values)


def _read_redistribute_moments(lift, mesh, space):
    on_moments = lift.on_moment_domain()
    moment_space = on_moments._domain_interior
    moments = HarmonicMomentFlux(values=_rng(4).standard_normal(moment_space.shape), space=moment_space, L=1)
    assert moments.head_layout.L == 1  # the leg holds under a MomentField refusal that derives L from the space
    return np.asarray(on_moments._redistribute_moments(moments).values)


#: (lift key, reader name) -> reader; ``F`` has no scattering transfer, so it carries two readers and ``S_L1`` four.
_READERS = {
    ("F", "combine"): _read_combine,
    ("F", "transpose"): _read_transpose,
    ("S_L1", "combine"): _read_combine,
    ("S_L1", "transpose"): _read_transpose,
    ("S_L1", "redistribute_ordinates"): _read_redistribute_ordinates,
    ("S_L1", "redistribute_moments"): _read_redistribute_moments,
}


@pytest.mark.rests_on("tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py::test_every_lift_divides_by_the_correctly_rounded_mass")
@pytest.mark.parametrize("n", _RULES_WHERE_THE_SUMS_DIFFER)
@pytest.mark.parametrize(("key", "reader"), list(_READERS), ids=[f"{k}-{r}" for k, r in _READERS])
def test_every_reader_reads_the_sections_mass_object(n, key, reader, monkeypatch):
    """Doubling the section's mass halves the reader's output exactly (division by ``2M`` is division by ``M`` then by 2, exact in
    binary floating point); a reader that sums the weights itself does not move. Each leg first asserts activation (a non-zero
    output, else a zero reading passes). Red today: all 6 legs, on GL12, read ``before`` where ``before / 2`` is owed."""
    mesh, _, space, lifts = _posed(n)
    lift = lifts[key]
    read = _READERS[(key, reader)]
    section = AxisSectionOperator.of(lift._codomain_interior, "angular")
    before = read(lift, mesh, space)
    assert np.abs(before).max() > 0.0
    with monkeypatch.context() as m:
        m.setattr(section, "_mass", 2.0 * section.mass)  # a plain assignment inside the context is not tracked
        after = read(lift, mesh, space)
    np.testing.assert_array_equal(after, before / 2.0)
