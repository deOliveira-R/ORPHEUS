r"""**#708 / G3** -- a scalar source does not add to an angular one; the isotropic source enters through the angular section.

``ScalarSourceSink`` and ``AngularSourceSink`` carry different units (``SCALAR_RATE_UNITS`` against ``ANGULAR_RATE_UNITS``, differing by
``1/sr``), so their sum is a dimensional error until the scalar is carried across the angular axis by the section ``E``, which
divides by the mass: ``E`` is the only map from the first to the second (``ontology_w0_spec.md`` 5.7.4). Before the unit the cross-class
arm ``AngularSourceSink.__add__`` / ``__radd__`` (with ``CONTAINS``) summed them with no division. After it no leaf overrides ``__add__``
and ``Field._check_partner`` refuses the sum with ``TypeError: ... requires a same-class partner``, a ``raise`` that survives
``python -O``. The refusal is not a new guard (no ``ELEGANCE-DEBT`` tag): the same-class gate already exists.

Legs, with the first red at the base tip ``954c4a4f`` (``-O``):

1. ``test_the_mixed_sum_is_refused[iso+aniso | aniso+iso]``: a CHILD interpreter (``python -O -P -c``, ``PYTHONPATH=<root>:<repository>``)
   builds both operands and evaluates the sum; it must exit non-zero with ``same-class partner`` on stderr. RED today: both orders
   construct and return an ``AngularSourceSink`` (exit 0). The child prints ``orpheus.__file__`` first (asserted under ``<root>``) and
   the sentinel ``IMPORTED`` after its last import (asserted), so a mutant that does not import is never read as red evidence.
2. ``test_the_isotropic_source_enters_through_the_sections_mass``: ``F._combine(iso, aniso).values`` equals
   ``broadcast(iso / math.fsum(w)) + aniso`` bit for bit on GL7 (a rule on which ``np.sum != fsum``). RED today: the divisor is the
   ``np.sum`` one. The within-class sum ``AngularSourceSink + AngularSourceSink`` is the control (green before and after).
3. ``test_the_coherence_refusal_moved_into_the_section``: an ``iso`` on a stretched slab's marginal space (the same shape, unequal
   content) is refused by ``_combine`` with a ``ValueError`` matching ``is not a member``. RED today on the message: the arm refuses
   with its own ``ValueError`` (``requires the iso operand's space ...``). The ``.values`` form of the rewrite would accept it.

The retirement of the names ``CONTAINS``, ``__add__`` and ``__radd__`` is held by the rows of ``tests/gates/_retired/rows_w0_angular_mass.py``
(and by no leg here: ``hasattr`` is not read by the engine).

Kinds (``vv-principles``): leg 1 a THEOREM (units: a sum of rates of different units is refused), leg 2 a REFERENCE (platform-free
``math.fsum``), leg 3 a THEOREM (space coherence). Layer: ``orpheus.data``, ``geometry``, ``mesh``, ``numerics``, ``transport`` only.

Mutations, by the mechanism of 1.3 (the reviewers apply them; ``_run_child(order, root)`` takes ``root`` for that purpose):

(i) the ``-O`` leg: in a ``.py``-only copy of ``orpheus/`` under ``T/orpheus/``, replace the search text
    ``if not self._is_same_class(other):\n            raise TypeError(`` (once in ``numerics/field.py``) by
    ``if not self._is_same_class(other):\n            assert False, (`` and run leg 1 with ``root=T``. The child imports (``IMPORTED``);
    under ``-O`` the ``assert`` is stripped, the sum then raises ``ValueError: ... equal space`` from the NEXT check of
    ``_check_partner``, whose stderr has no ``same-class partner``: the fragment assertion reddens, which is why the leg asserts the
    fragment and not only the exit status.
(ii) the arm re-planted: after the unit, in the copy of ``transport/fields/angular_source_sink.py`` restore the three imports and the
    constant ``ANGULAR_SOURCE_SINK_ARM`` of ``tests/gates/_retired/_controls.py`` after ``UNITS``; the ``aniso + iso`` child then exits 0
    and prints ``AngularSourceSink`` (only that order constructs: ``Field.__add__`` raises for ``iso + aniso``). (ii') plus
    ``BULK_FIELD_CONTAINS`` in ``_bases.py``: both children exit 0, and the ``BulkField.CONTAINS`` and ``AngularSourceSink.__add__`` rows
    of the engine redden.
(iii) leg 3, in process: ``monkeypatch.setattr(AngularLift, "_combine", f)`` with ``f`` the three-line body using
    ``section.apply(iso.values)``; only leg 3 reddens (the stretched operand passes the shape check of the bare-array path).
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.axis_marginal import AxisRetractionOperator
from orpheus.transport.fields import AngularSourceSink, ScalarSourceSink
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from tests.gates.transport._phase_space_cases import angular_full_field_space, posing_quadrature, posing_scheme

pytestmark = pytest.mark.foundation

_REPO = Path(__file__).resolve().parents[4]

#: The child: builds the refused input (the imports of 5.7.2 item 1), prints the sentinel, then evaluates the sum.
_CHILD = """
import orpheus
print(orpheus.__file__)
import numpy as np
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import Mesher, CellsByCount
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.numerics.axis_marginal import AxisRetractionOperator
from orpheus.transport.fields import ScalarSourceSink, AngularSourceSink
from tests.gates.transport._phase_space_cases import angular_full_field_space, posing_quadrature, posing_scheme
print('IMPORTED')
sn = MaterialMesh(Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(CellsByCount.uniform_width(20)).mesh, {{0: get_mixture("A", "2g")}})
interior = angular_full_field_space(sn, posing_quadrature(sn, 4), posing_scheme(sn)).interior_space
scalar_space = AxisRetractionOperator.of(interior, "angular").codomain
iso = ScalarSourceSink(values=np.ones(scalar_space.shape), space=scalar_space)
aniso = AngularSourceSink.zeros(interior)
print(type({expr}).__name__)
"""


def _run_child(expr: str, root: Path = _REPO) -> subprocess.CompletedProcess:
    """Run the child under ``python -O -P`` with ``<root>`` (a directory CONTAINING ``orpheus/``) ahead of the repository root.

    ``-P`` drops the working-directory entry that ``-c`` puts ahead of ``PYTHONPATH``; the CALLER asserts the printed
    ``orpheus.__file__`` lies under ``root``, else a run that imported the main tree would read as the mutant's."""
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([str(root), str(_REPO)])}
    return subprocess.run(
        [sys.executable, "-O", "-P", "-c", _CHILD.format(expr=expr)], env=env, capture_output=True, text=True, cwd=str(_REPO),
    )


@pytest.mark.rests_on("tests/gates/numerics/test_field.py::test_cross_class_addition_rejected")
@pytest.mark.parametrize("expr", ["iso + aniso", "aniso + iso"])
def test_the_mixed_sum_is_refused(expr):
    """``ScalarSourceSink + AngularSourceSink`` raises ``TypeError ... same-class partner`` in both orders, under ``python -O``."""
    done = _run_child(expr)
    lines = done.stdout.splitlines()
    assert lines and Path(lines[0]).is_relative_to(_REPO), f"child imported another tree: {lines[:1]}"
    assert "IMPORTED" in lines, f"the child never reached the expression:\n{done.stderr}"
    assert done.returncode != 0, f"{expr!r} constructed: {lines[-1:]}"
    assert "same-class partner" in done.stderr, done.stderr


def _slab(width: float = 1.0) -> MaterialMesh:
    mesh = Mesher(StructuredGeometry.from_homogeneous(width, 0)).partition(CellsByCount.uniform_width(20)).mesh
    return MaterialMesh(mesh, {0: get_mixture("A", "2g")})


def _fission_lift(mesh: MaterialMesh, n: int):
    quadrature = posing_quadrature(mesh, n)
    space = angular_full_field_space(mesh, quadrature, posing_scheme(mesh))
    return FissionOperator.from_solver_data(mat_xs=mesh.mat_xs, space=space), quadrature


@pytest.mark.rests_on("tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py::test_every_lift_divides_by_the_correctly_rounded_mass")
def test_the_isotropic_source_enters_through_the_sections_mass():
    """``F._combine(iso, aniso).values == broadcast(iso / fsum(w)) + aniso`` bit for bit on GL7, where ``np.sum != fsum``; the
    leg CALLS ``_combine``, so it survives a later move of the combine to a named entry. The within-class sum is the control."""
    lift, quadrature = _fission_lift(_slab(), 7)
    weights = np.asarray(quadrature.weights, dtype=float)
    assert float(weights.sum()) != math.fsum(weights)
    rng = np.random.default_rng(7)
    iso_values = rng.standard_normal(lift._scalar_interior_space.shape)
    aniso_values = rng.standard_normal(lift._codomain_interior.shape)
    out = lift._combine(
        ScalarSourceSink(values=iso_values, space=lift._scalar_interior_space),
        AngularSourceSink(values=aniso_values, space=lift._codomain_interior),
    )
    expected = np.broadcast_to((iso_values / math.fsum(weights))[None], aniso_values.shape) + aniso_values
    np.testing.assert_array_equal(np.asarray(out.values), expected)


def test_a_within_class_angular_sum_is_unchanged():
    """Control, green before and after: ``AngularSourceSink + AngularSourceSink`` adds the values."""
    lift, _ = _fission_lift(_slab(), 4)
    space = lift._codomain_interior
    rng = np.random.default_rng(3)
    a, b = rng.standard_normal(space.shape), rng.standard_normal(space.shape)
    total = AngularSourceSink(values=a, space=space) + AngularSourceSink(values=b, space=space)
    np.testing.assert_array_equal(np.asarray(total.values), a + b)


@pytest.mark.rests_on("tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py::test_every_lift_divides_by_the_correctly_rounded_mass")
def test_the_coherence_refusal_moved_into_the_section():
    """An isotropic operand on a stretched slab's marginal space (the same shape as the lift's, unequal content) is refused by
    ``_combine`` as not a member of the section's domain. Red today on the message only: the retired arm refused with its own
    ``ValueError`` (``requires the iso operand's space ...``)."""
    lift, _ = _fission_lift(_slab(), 4)
    stretched = _slab(2.0)
    stretched_interior = angular_full_field_space(
        stretched, posing_quadrature(stretched, 4), posing_scheme(stretched),
    ).interior_space
    foreign = AxisRetractionOperator.of(stretched_interior, "angular").codomain
    assert foreign.shape == lift._scalar_interior_space.shape and foreign != lift._scalar_interior_space
    with pytest.raises(ValueError, match="is not a member"):
        lift._combine(ScalarSourceSink(values=np.ones(foreign.shape), space=foreign), None)
