r"""#710 (W1): on a measure arm an ``Axis``'s weights are the origin measure's.

``Axis.weights`` is an init field stored beside the origin that generated it, while the axis's ``nodes`` are a
derived copy read off the origin; ``__post_init__`` checks shape and finiteness only. The refusal is ONE
``raise ValueError`` at the end of ``Axis.__post_init__``, tagged ``ELEGANCE-DEBT[guard] #710 (W1)``, over the
fifth read ``determined_weights`` of the origin arms (``FromMeasure`` and ``FromQuadrature`` return their
measure's weights, the other three ``None``). Specification: ``.claude/plans/ontology_w0_spec.md`` 5.1 to 5.3.

Law: for ``FromMeasure(m)`` and ``FromQuadrature(Q)`` the canonical weights (``None`` read as ones) are
bit-equal to ``m.weights`` / ``Q.measure.weights``. A weight-less axis over a non-counting measure is refused
too (decision 10): ``None`` there dresses the measure as the counting measure, which is the defect.

Activation: the weights 7.0 differ from the Gauss-Legendre weights at every node, so a build comparing any
subset fails; the RANK-2 admitted axis is the leg a comparison that drops the reshape cannot pass (with rank 1
every admitted input is identical with or without it). The ``FromBasis`` and ``Declared`` arms are NOT covered
(95 constructions on them carry their own weights).

Mutations: the refusal replaced by ``assert False, (`` in a copy (leg 1 constructs under ``-O``); the
``np.reshape`` of the comparison dropped in a copy (the rank-2 admitted leg is refused, no rank-1 leg is);
``FromQuadrature.determined_weights`` made ``None`` in process (only the ``FromQuadrature`` leg reddens).

Every assertion is a function call or an ``if ... pytest.fail`` (``python -O`` safe).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, Declared, FromMeasure, FromQuadrature
from orpheus.numerics.basis import LegendreBasis
from orpheus.numerics.quadrature import Quadrature, gauss_legendre_on_polar_orbit
from tests.gates import _refusal_gates as rg

pytestmark = pytest.mark.foundation

FRAGMENT = "disagree with the measure its origin records"
_AXIS_FILE = "orpheus/numerics/axis.py"

_CHILD_IMPORTS = (
    "import numpy as np\n"
    "from orpheus.numerics.axis import Axis, FromMeasure, FromQuadrature\n"
    "from orpheus.numerics.quadrature import Quadrature, gauss_legendre_on_polar_orbit\n"
    "m = gauss_legendre_on_polar_orbit(4, 'x')\n"
    "Q = Quadrature.gauss_legendre(n_ordinates=4)\n"
)
_CHILD_INPUTS = {
    "FromMeasure, weights 7": (
        "Axis('angular', (m.n_points,), weights=np.full(m.n_points, 7.0), origin=FromMeasure(m))"
    ),
    "FromQuadrature, weights 7": (
        "Axis('angular', (Q.n_ordinates,), weights=np.full(Q.n_ordinates, 7.0), origin=FromQuadrature(Q))"
    ),
    "FromMeasure, weights None": "Axis('angular', (m.n_points,), origin=FromMeasure(m))",
}
_RANK2_ADMITTED = (
    "Axis('angular', (2, 2), weights=np.reshape(m.weights, (2, 2)), origin=FromMeasure(m))"
)


@pytest.fixture(scope="module")
def measure():
    return gauss_legendre_on_polar_orbit(4, "x")


@pytest.fixture(scope="module")
def quadrature():
    return Quadrature.gauss_legendre(n_ordinates=4)


# ----------------------------------------------------------------- admitted legs


def test_the_measure_the_axis_is_minted_from_is_not_the_counting_measure(measure):
    """CONTROL: the Gauss-Legendre weights are not all ones, so ``weights=None`` over them is a real disguise."""
    if np.array_equal(measure.weights, np.ones(measure.n_points)):
        pytest.fail("CONTROL INVALID: the fixture measure is the counting measure")


def test_the_two_minting_factories_are_admitted(measure, quadrature):
    Axis.of_measure(measure, "angular")
    Axis.of_quadrature(quadrature)


def test_the_origins_own_weights_passed_explicitly_are_admitted(measure, quadrature):
    """``Axis(..., weights=m.weights, origin=FromMeasure(m))`` equals what the factory stores, so it passes."""
    Axis("angular", (measure.n_points,), weights=measure.weights, origin=FromMeasure(measure))
    Axis(
        "angular", (quadrature.n_ordinates,), weights=quadrature.measure.weights,
        origin=FromQuadrature(quadrature),
    )


def test_a_rank_two_axis_over_the_reshaped_weights_is_admitted(measure):
    """It constructs today, and ``np.array_equal(np.reshape(w, (2, 2)), w)`` is False: the leg that sees a
    comparison dropping the reshape."""
    Axis("angular", (2, 2), weights=np.reshape(measure.weights, (2, 2)), origin=FromMeasure(measure))


def test_the_arms_that_carry_their_own_weights_are_untouched():
    """``FromBasis`` (``of_basis`` stores the basis's coefficient weights, ``[4 pi, 4 pi / 3, 4 pi / 5]`` at L = 2)
    and ``Declared`` keep ``weights`` as an init field: the refusal ranges over the two measure arms only."""
    axis = Axis.of_basis(LegendreBasis(L=2, axis="z"))
    if axis.weights is None:
        pytest.fail("CONTROL INVALID: of_basis stores no weights")
    Axis("x", (3,), weights=np.full(3, 7.0), origin=Declared(BasisKind.NODAL))


# ----------------------------------------------------------------- refused inputs, in process


def test_FromMeasure_with_foreign_weights_is_refused(measure):
    """First red: weights 7 against the origin's 0.3478548451... construct today."""
    with pytest.raises(ValueError, match=FRAGMENT):
        Axis("angular", (measure.n_points,), weights=np.full(measure.n_points, 7.0), origin=FromMeasure(measure))


def test_FromQuadrature_with_foreign_weights_is_refused(quadrature):
    with pytest.raises(ValueError, match=FRAGMENT):
        Axis(
            "angular", (quadrature.n_ordinates,), weights=np.full(quadrature.n_ordinates, 7.0),
            origin=FromQuadrature(quadrature),
        )


def test_FromMeasure_with_no_weights_over_a_non_counting_measure_is_refused(measure):
    """Decision 10: ``None`` reads as ones, which silently dresses a non-counting measure as the counting one."""
    with pytest.raises(ValueError, match=FRAGMENT):
        Axis("angular", (measure.n_points,), origin=FromMeasure(measure))


# ----------------------------------------------------------------- cold child, -O


def test_the_unmutated_child_constructs_the_rank_two_axis():
    """CONTROL of the child legs: an admitted input imports and constructs under ``python -O``."""
    child = rg.run_child(_CHILD_IMPORTS, _RANK2_ADMITTED)
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.constructed(child, FRAGMENT):
        pytest.fail(f"the admitted rank-2 axis was refused in the child:\n{child.stderr}")


@pytest.mark.parametrize("which", sorted(_CHILD_INPUTS))
def test_the_refusal_survives_python_O(which):
    child = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which])
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.refused(child, FRAGMENT):
        pytest.fail(f"{which!r} was not refused under -O (exit {child.returncode}):\n{child.stderr}")


# ----------------------------------------------------------------- mutations


@pytest.mark.parametrize("which", sorted(_CHILD_INPUTS))
def test_mutation_the_refusal_replaced_by_an_assert_constructs_under_python_O(which, tmp_path):
    """Leg (1) mutation: the ``raise`` after the #710 tag becomes ``assert False, (``; the refused input constructs."""
    honest = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which])
    if not rg.refused(honest, FRAGMENT):
        pytest.fail(f"the unmutated tree does not refuse {which!r}, so the mutation proves nothing:\n{honest.stderr}")
    root = rg.copy_orpheus_py(rg.import_root(), tmp_path)
    rg.remove_refusal(root, _AXIS_FILE, 710)
    mutant = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which], root=root)
    if not rg.imported_from(mutant, root):
        pytest.fail(f"the mutant did not import from the copy (never red evidence):\n{mutant.stdout}\n{mutant.stderr}")
    if not rg.constructed(mutant, FRAGMENT):
        pytest.fail(f"the mutant still refused {which!r} (exit {mutant.returncode}):\n{mutant.stderr}")


def test_mutation_dropping_the_reshape_refuses_the_rank_two_axis(tmp_path):
    """The comparison's ``np.reshape(np.asarray(wanted, dtype=float), shape)`` replaced, as a whole expression, by
    ``np.asarray(wanted, dtype=float)`` (deleting ``np.reshape(`` alone leaves an unmatched parenthesis, a
    ``SyntaxError`` that would redden every leg for the wrong reason): the admitted rank-2 axis is then refused
    with the refusal's own fragment, which no rank-1 leg gives."""
    honest = rg.run_child(_CHILD_IMPORTS, _RANK2_ADMITTED)
    if not rg.constructed(honest, FRAGMENT):
        pytest.fail(f"the unmutated tree refuses the rank-2 axis, so the mutation proves nothing:\n{honest.stderr}")
    root = rg.copy_orpheus_py(rg.import_root(), tmp_path)
    rg.replace_once(
        root, _AXIS_FILE,
        "np.reshape(np.asarray(wanted, dtype=float), shape)", "np.asarray(wanted, dtype=float)",
    )
    mutant = rg.run_child(_CHILD_IMPORTS, _RANK2_ADMITTED, root=root)
    if not rg.imported_from(mutant, root):
        pytest.fail(f"the mutant did not import from the copy (never red evidence):\n{mutant.stdout}\n{mutant.stderr}")
    if not rg.refused(mutant, FRAGMENT):
        pytest.fail(f"the mutant still admits the rank-2 axis (exit {mutant.returncode}):\n{mutant.stderr}")


def test_mutation_FromQuadrature_determining_nothing_reddens_only_its_leg(measure, quadrature, monkeypatch):
    """In process: ``FromQuadrature.determined_weights`` returning ``None`` makes the FromQuadrature input construct;
    the FromMeasure input is still refused (the legs are disjoint)."""
    if not hasattr(FromQuadrature, "determined_weights"):
        raise rg.Uninstallable("FromQuadrature has no determined_weights: the fifth read is not installed yet")
    with pytest.raises(ValueError, match=FRAGMENT):
        Axis(
            "angular", (quadrature.n_ordinates,), weights=np.full(quadrature.n_ordinates, 7.0),
            origin=FromQuadrature(quadrature),
        )
    monkeypatch.setattr(FromQuadrature, "determined_weights", property(lambda self: None))
    Axis(
        "angular", (quadrature.n_ordinates,), weights=np.full(quadrature.n_ordinates, 7.0),
        origin=FromQuadrature(quadrature),
    )
    with pytest.raises(ValueError, match=FRAGMENT):
        Axis("angular", (measure.n_points,), weights=np.full(measure.n_points, 7.0), origin=FromMeasure(measure))
