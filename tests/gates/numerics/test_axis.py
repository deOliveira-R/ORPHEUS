r"""Intrinsic laws of the space-factor axis (campaign 1, CS1 step 1).

Every DEFINING law of the ``Axis`` / ``EnergyAxis`` concept gets a direct
test here, not merely a usage test — the project's "test a math concept's
intrinsic properties" standard. The concept: an axis is
(index shape, factor measure, basis kind, generator identity), frozen,
structurally identified PER SUBCLASS, with ``weights=None`` meaning the
COUNTING measure deliberately (an axis has no "unbound" state), and
weights stored CANONICALLY (all-ones collapses to ``None``; ``-0.0``
normalized; non-finite refused) so one measure has one spelling and one
identity.

Gate ids A1–A12 refer to the CS1 battery of record
(``scratch/cs1_verification_plan.md`` §2); the canonicalization gates
realize the Q-T1/Q-T3 rulings recorded in
``.claude/plans/cs1_energy_space_design.md`` §T-R.

No theory ``:label:`` — these are software/math invariants of a type, so
``foundation``, never ``verifies(...)``.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.axis import Declared, Axis, BasisKind, EnergyAxis
from orpheus.numerics.axis import Modal
from orpheus.numerics.metric import DiagonalMetric, NoForm

pytestmark = pytest.mark.foundation

#: Two groups, strictly DESCENDING (the canonical fast-first convention,
#: ``EnergyGrid`` refuses anything else).
_EDGES_2G = np.array([1.0e7, 1.0e3, 1.0e-3])


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def test_axis_is_frozen() -> None:
    """A1 — the axis is a VALUE: rebinding a field raises.

    LOADED leg (there is no control): the mutation is the witness. An axis
    that is hashable and mutable is an illegal state — its identity could
    change after it has been used as a dict key or compared inside a
    space's ``__eq__``.
    """
    ax = Axis("x", (3,), origin=Declared(BasisKind.NODAL))
    with pytest.raises(dataclasses.FrozenInstanceError):
        ax.label = "y"  # type: ignore[misc]


def test_no_form_is_the_counting_form() -> None:
    """A2 — a declared origin with no form IS the counting form, deliberately.

    The counting measure has ONE spelling: ``Declared(kind)`` with no form is
    ``DiagonalMetric.counting(shape)``, so the no-form-spelled and the
    ones-spelled counting axes are not merely equivalent, they are EQUAL (one
    identity), and no consumer can ever observe two spellings of the
    identity metric.
    """
    none_spelled = Axis("x", (3,), origin=Declared(BasisKind.NODAL))
    ones_spelled = Axis("x", (3,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.ones(3))))
    _require(np.array_equal(none_spelled.form.weights, np.ones(3)), "no form must read the counting form over the shape")
    _require(
        np.array_equal(ones_spelled.form.weights, np.ones(3)),
        "an all-ones form is the counting form (one spelling per measure)",
    )
    _require(
        none_spelled == ones_spelled and hash(none_spelled) == hash(ones_spelled),
        "the two spellings of the counting measure must be ONE identity",
    )


def test_weights_are_stored_canonically_and_read_only() -> None:
    """Q-T1/Q-T3 — the stored measure is canonical bytes, immutable.

    ``-0.0`` normalizes to ``+0.0`` (one measure, one byte pattern — the
    name digest and the hash read these bytes), and the stored array
    refuses in-place writes.
    """
    neg_zero = Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([-0.0, 2.0]))))
    pos_zero = Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([0.0, 2.0]))))
    _require(
        neg_zero == pos_zero and hash(neg_zero) == hash(pos_zero),
        "-0.0 and +0.0 spell ONE measure and must be ONE identity",
    )
    assert isinstance(neg_zero.form, DiagonalMetric)  # narrowing for the type checker
    w = neg_zero.form.weights
    _require(
        w.tobytes() == np.array([0.0, 2.0]).tobytes(),
        "stored weights must be the canonical (+0.0) byte pattern",
    )
    _require(not w.flags.writeable, "stored weights must be read-only")
    with pytest.raises((ValueError, RuntimeError)):
        w[0] = 5.0


def test_non_finite_weights_are_refused() -> None:
    """Q-T1 — a factor measure has FINITE weights; nan/inf are refused.

    (No non-negativity guard, deliberately: CS2's quadrature axes legally
    carry signed weights.)
    """
    for bad in (np.array([np.nan, 1.0]), np.array([np.inf, 1.0])):
        with pytest.raises(ValueError, match="finite"):
            Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(bad)))
    # The signed-weight door stays OPEN (the anti-claim of the guard):
    signed = Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([-0.5, 1.5]))))
    assert isinstance(signed.form, DiagonalMetric)
    w = signed.form.weights
    _require(bool(w[0] == -0.5), "signed weights are legal (quadrature axes)")




def test_synthetic_equality_is_ng_only() -> None:
    """A4 — ``synthetic(ng)`` identity is exactly ``ng``."""
    _require(EnergyAxis.synthetic(2) == EnergyAxis.synthetic(2), "same ng must be equal")
    _require(
        hash(EnergyAxis.synthetic(2)) == hash(EnergyAxis.synthetic(2)),
        "equal synthetic axes must hash equal",
    )
    _require(EnergyAxis.synthetic(2) != EnergyAxis.synthetic(3), "ng differs => unequal")


def test_synthetic_and_of_edges_at_the_same_ng_are_UNEQUAL() -> None:
    r"""A5 ⭐ — the sharp pair: same ``ng``, same shape, different axes.

    ``synthetic(2)`` and ``of_edges(<2-group edges>)`` describe the same
    INDEX SET and different PARTITIONS. Q2 rules identity = ``ng`` + edges
    CONTENT, so they must differ; if identity collapsed to ``ng`` the
    derived space NAME would collide and the axes tuples would compare
    equal — since the identity flip (CS4c step 6) the axes tuple IS the
    space identity, and until it the derived name was — so two physically
    different spaces would compare equal and compose silently either way.
    """
    synthetic = EnergyAxis.synthetic(2)
    gridded = EnergyAxis.of_edges(_EDGES_2G)
    _require(synthetic.shape == gridded.shape, "precondition: same index set")
    _require(synthetic != gridded, "identity collapsed to ng — the partition was dropped")


def test_energy_axis_hash_agrees_with_eq() -> None:
    """A6 — the A3/A4/A5 population round-trips through a ``set`` and a
    ``dict`` with exactly the expected cardinality.

    The A3 pair is two axes over two DISTINCT, equal edge arrays: equality is the
    edges' CONTENT, not the array's identity (the numerics half of A3; the data half,
    two fresh ``EnergyGrid`` objects of one mixture, is
    ``tests/gates/data/test_energy_grid.py::test_energy_axis_equality_is_edges_CONTENT_not_grid_identity``,
    moved by the layer DAG Phase 2 refresh)."""
    edges_a, edges_b = _EDGES_2G.copy(), _EDGES_2G.copy()
    _require(edges_a is not edges_b, "the A3 pair is built over two distinct arrays")
    population = [
        EnergyAxis.of_edges(edges_a),
        EnergyAxis.of_edges(edges_b),  # equal to the previous (A3)
        EnergyAxis.synthetic(2),
        EnergyAxis.synthetic(2),  # equal to the previous (A4)
        EnergyAxis.synthetic(3),
        EnergyAxis.of_edges(np.array([2.0e7, 1.0, 1.0e-5])),
    ]
    distinct = {population[0], population[2], population[4], population[5]}
    _require(
        len(set(population)) == 4,
        f"expected 4 distinct axes in the population, got {len(set(population))}",
    )
    _require(len({ax: None for ax in population}) == 4, "dict keying must agree with set")
    _require(len(distinct) == 4, "the four representatives must stay distinct")


def test_axis_is_immune_to_mutation_of_the_source_array() -> None:
    r"""A7 — the axis defensively copies its edges/weights.

    ``Mixture`` is a plain (NON-frozen) ``@dataclass`` and ``eg`` is a bare
    ndarray, so the caller's array is live. A hashable value object that
    aliases a mutable buffer is an illegal state: its ``__hash__`` silently
    changes after it has been stored in a dict.

    The energy leg mints through ``EnergyAxis.of_edges`` directly, so no
    grid stands between the caller's buffer and the axis: ``of_edges``
    reads the buffer with ``np.asarray`` (an alias, not a copy, for a
    contiguous float array), and the axis's own defensive copy in
    ``__post_init__`` is the sole guard on the path. (Through the old
    ``EnergyGrid`` route the guard was the same: the grid's
    ``np.ascontiguousarray`` also aliases such a buffer.)
    """
    edges = _EDGES_2G.copy()
    axis_before = EnergyAxis.of_edges(edges)
    reference = EnergyAxis.of_edges(_EDGES_2G)
    hash_before = hash(axis_before)
    edges[0] = 9.9e9  # mutate the caller's buffer
    _require(hash(axis_before) == hash_before, "hash moved under source mutation")
    _require(axis_before == reference, "equality moved under source mutation")

    weights_src = np.array([2.0, 5.0])
    weighted = Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(weights_src)))
    weighted_ref = Axis("x", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([2.0, 5.0]))))
    weights_src[0] = 7.0
    _require(weighted == weighted_ref, "weights identity moved under source mutation")


def test_the_data_reading_constructors_are_retired_from_the_numerics_axis() -> None:
    """layer_dag.md P2 — the numerics axis names nothing in data.

    ``EnergyAxis.from_grid`` (read an ``EnergyGrid``) and
    ``EnergyAxis.from_materials`` (the energy-arm rule over ``Mixture``)
    were re-homed to data as ``EnergyGrid.as_axis`` and
    ``orpheus.data.materials.energy_axis_of``; their return here would put
    the ``numerics -> data`` edge back.
    """
    for retired in ("from_grid", "from_materials"):
        _require(not hasattr(EnergyAxis, retired), f"EnergyAxis.{retired} is back")


def test_of_edges_derives_the_group_count_and_the_label() -> None:
    """``of_edges`` takes CONTENT alone: ng = len(edges) - 1, label "energy"."""
    axis = EnergyAxis.of_edges(_EDGES_2G)
    _require(
        axis == EnergyAxis("energy", (2,), edges=_EDGES_2G),
        f"of_edges did not derive ng and the label: {axis!r}",
    )


@pytest.mark.parametrize("edges", [np.array([1.0]), np.array([[2.0, 1.0], [1.0, 0.5]])], ids=["one-edge", "rank-2"])
def test_of_edges_refuses_fewer_than_two_edges_or_a_non_vector(edges: np.ndarray) -> None:
    """A single boundary bounds no group; a rank-2 array is no edge vector."""
    with pytest.raises(ValueError, match="of_edges needs a 1-D array"):
        EnergyAxis.of_edges(edges)


def test_rank_zero_shape_is_admitted() -> None:
    """A8 — an axis may have index shape ``()``: one point, no array axis (the face-spatial
    factor of a 1-D face), with the counting form of a single weight."""
    point = Axis("x", (), origin=Declared(BasisKind.NODAL))
    _require(point.shape == () and np.array_equal(point.form.weights, np.ones(())), "a rank-0 axis carries the counting form of one weight")


def test_weights_shape_must_match_shape() -> None:
    """A9 — the factor measure lives over exactly ``shape``."""
    with pytest.raises(ValueError, match="shape"):
        Axis("x", (3,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.ones(4) * 2.0)))


def test_basis_kind_participates_in_identity() -> None:
    """A10 — NODAL and MODAL are different axes.

    The cone metadata reads ``kind``; if kind were metadata-only, a modal
    space and a nodal space would compare equal and the CS1 step-4 cone
    refusal would be reachable through the wrong space.
    """
    nodal = Axis("x", (3,), origin=Declared(BasisKind.NODAL))
    modal = Axis("x", (3,), origin=Declared(Modal(np.eye(3)[0])))
    _require(nodal != modal, "basis kind must participate in identity")


def test_identity_is_structural_PER_SUBCLASS() -> None:
    """A11 — an ``EnergyAxis`` never equals a generic ``Axis`` carrying the
    same tuple of fields.

    "Structural per subclass, from day one" (§A). Free if ``EnergyAxis`` is
    its own dataclass (dataclass ``__eq__`` requires ``other.__class__ is
    self.__class__``); NOT free if ``EnergyAxis`` is a factory returning a
    bare ``Axis`` — which is exactly the design decision this pins.
    """
    energy = EnergyAxis.synthetic(2)
    generic = Axis("energy", (2,), origin=Declared(BasisKind.NODAL))
    _require(energy != generic, "subclass identity collapsed (EnergyAxis == Axis)")
    _require(generic != energy, "subclass identity must be symmetric")


def test_weights_are_part_of_identity() -> None:
    r"""A12 — same label, same shape, different measure ⟹ different axis.

    The §B clause-2 retrodiction at the axis tier ("a genuine one-cell slab
    keeps its axis with weight V != 1, distinguished from the quotient point
    by MEASURE") and the Q8 migration doctrine ("metric differences imply
    space differences"). ⚠ Per F2 this distinction is INVISIBLE to ``.H``
    on a one-cell space (a scalar metric commutes with everything), so
    identity is the only instrument that can carry it.
    """
    counting = Axis("spatial", (1,), origin=Declared(BasisKind.NODAL))  # quotient point, weight 1
    one_cell = Axis("spatial", (1,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([2.0]))))
    other_cell = Axis("spatial", (1,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array([3.0]))))
    _require(counting != one_cell, "V != 1 must differ from the quotient point")
    _require(one_cell != other_cell, "different measures must be different axes")


def test_energy_axis_has_no_form_and_no_modal_kind() -> None:
    """The counting theorem and the faces reading, enforced at construction.

    An ``EnergyAxis`` carries no weights at all, so a weighted one is
    unspellable (the deliberately non-physical D4b control toy uses a generic
    ``Axis`` instead): its origin is fixed to ``Declared(NODAL, NoForm())``, and
    a ``weights=`` or ``form=`` keyword is a ``TypeError``; a MODAL ``EnergyAxis``
    would contradict the groups-are-cells reading.
    """
    for keyword in ("weights", "form"):
        with pytest.raises(TypeError):
            EnergyAxis("energy", (2,), **{keyword: np.array([2.0, 5.0])})  # pyright: ignore[reportCallIssue]
    # NODAL by construction: the kind is fixed by the class, so a modal
    # energy axis is not refused but unspellable (S7, ruling 82).
    with pytest.raises(TypeError, match="origin"):
        EnergyAxis("energy", (2,), origin=Declared(Modal(np.eye(2)[0])))  # pyright: ignore[reportCallIssue]
    _require(isinstance(EnergyAxis("energy", (2,)).form, NoForm), "an EnergyAxis carries no form")
    _require(EnergyAxis("energy", (2,)) == EnergyAxis.synthetic(2), "two spellings of the synthetic axis are one axis")
