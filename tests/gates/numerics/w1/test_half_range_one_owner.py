r"""G02a (value row and call arm), G02d, G02f and G02g: the half-range has one owner and the owner is exact (N02, spec 6.3 and 7.3).

``numerics/tangential.py::half_range(quadrature, axis_index, *, sign=+1)`` reads ``sign * quadrature.mean_axis_cosine(axis_index)`` and
returns a frozen ``HalfRange(outflow, inflow, tangential, n)`` of three ``Coord`` parts that partition the ordinate index set. The
comparisons are exact: ``c > 0``, ``c < 0``, ``c == 0``.

Every row here reads integer index sets or the ``ast`` of a source file: no float of a LAPACK or libm routine is compared, so each row is
exact on any host. ``half_range`` and ``HalfRange`` are imported inside the test bodies, so that the ``ast`` rows are red for their own
reason (a consumer that never calls the owner, a compare of a cosine array) and not as a collection error.

First red on the tip ``37a82d78`` (each recorded with its command in the unit's report):

* value row, refusals, partition, ``sign()``: import-only (``half_range`` does not exist); the behavioural witness is that today's three
  selectors give ``[2]``, ``[1, 2]`` and ``[-1, 0, 1]`` on ``[-0.5, 9.5e-16, 0.5]`` and accept ``[nan, 0.5]``;
* the call arm: 0 of 8 consumers call ``half_range``;
* G02g: 13 ``Compare`` nodes read a cosine array in the scopes.

Each instrument is run on its mutant (``_half_range_helpers``: a tolerance classifier, an unguarded classifier, a copy of a consumer
with its call deleted, a selector with a compare put back) in a row that is green on any tree: a gate that cannot read a wrong answer
as wrong is not a gate.
"""

from __future__ import annotations

import ast
import dataclasses

import numpy as np
import pytest

from tests.gates import _half_range_helpers as H

pytestmark = pytest.mark.foundation


def _owner():
    from orpheus.numerics.tangential import half_range

    return half_range


def _half_range_class():
    from orpheus.numerics.tangential import HalfRange

    return HalfRange


# ── G02a value row ───────────────────────────────────────────────────


def test_the_owner_is_exact_on_a_round_off_cosine() -> None:
    r"""``[-0.5, 9.5e-16, 0.5]`` has outflow ``(1, 2)`` and inflow ``(0,)``; the subnormal stub reads ``5e-324`` as genuine, ``-0.0`` as zero."""
    bad = H.value_row_violations(_owner())
    assert not bad, f"half_range is not the exact classifier: {bad}"


@pytest.mark.parametrize(
    "band, reddens",
    [
        (1e-15, {"[-0.5, 9.5e-16, 0.5]", "[-0.5, 9.5e-16, 0.5] sign=-1", "[0.0, -0.0, 5e-324, -5e-324]"}),
        (1e-12, {"[-0.5, 9.5e-16, 0.5]", "[-0.5, 9.5e-16, 0.5] sign=-1", "[0.0, -0.0, 5e-324, -5e-324]"}),
        # 4 eps = 8.9e-16 < 9.5e-16: the old band already read the round-off cosine as genuine; only the subnormal stub reddens it.
        (4.0 * H.EPS, {"[0.0, -0.0, 5e-324, -5e-324]"}),
    ],
    ids=["1e-15", "1e-12", "4eps"],
)
def test_the_value_row_reads_a_tolerance_classifier_as_wrong(band: float, reddens: set[str]) -> None:
    r"""Control of the value row: each banded copy of the owner reddens it, and the ``4 eps`` copy only on the subnormal stub."""
    bad = H.value_row_violations(H.banded_owner(band))
    assert {b.split(": ")[0] for b in bad} == reddens, bad


def test_the_value_row_reads_the_exact_classifier_as_right() -> None:
    r"""Control of the control: the exact classifier with refusals reads clean, so the rows above are not red for any other reason."""
    assert H.value_row_violations(H.exact_owner) == []


# ── G02a call arm ────────────────────────────────────────────────────


def test_every_consumer_calls_the_owner() -> None:
    r"""In each of the 8 consumer functions the body contains a call named ``half_range``."""
    rows = H.consumer_calls(H.orpheus_root(), "half_range")
    assert len(rows) == 8, f"the consumer table has {len(rows)} rows, wanted 8"
    not_found = [(f, q) for f, q, calls in rows if calls is None]
    assert not not_found, f"consumer functions not found on this tree (a named object that does not exist): {not_found}"
    missing = [(f, q) for f, q, calls in rows if not calls]
    assert not missing, f"{len(missing)} of 8 consumers do not call half_range: {missing}"


def test_the_call_predicate_reads_a_deleted_call_as_missing_on_synthetic_text() -> None:
    r"""Control of the call predicate on text that exists on any tree: the call is seen, and its deletion is seen."""
    src = "def f(q):\n    return half_range(q, 0).outflow\n"
    assert H.calls_name(ast.parse(src), "half_range")
    assert H.calls_name(ast.parse("def f(q):\n    return mod.half_range(q, 0)\n"), "half_range")
    assert not H.calls_name(ast.parse(src.replace("half_range", "_deleted")), "half_range")
    assert not H.calls_name(ast.parse("def f(q):\n    return q.outflow\n"), "half_range")


@pytest.mark.parametrize("rel, qual", H.CONSUMERS, ids=[q for _, q in H.CONSUMERS])
def test_deleting_the_call_from_a_consumer_reddens_the_arm(rel: str, qual: str) -> None:
    r"""Mutation, 8 copies: the consumer's own source with its ``half_range`` calls renamed away reads as not calling the owner.

    The first assertion is the precondition (the honest consumer calls the owner), so this row is red until the consumer is re-pointed.
    """
    honest, mutant = H.consumer_source_without_call(H.orpheus_root(), rel, qual, "half_range")
    assert H.calls_name(ast.parse(honest), "half_range"), f"{qual} does not call half_range: nothing to delete"
    assert not H.calls_name(ast.parse(mutant), "half_range"), "the mutant still calls the owner: the mutation did not bite"


# ── G02d ─────────────────────────────────────────────────────────────


def test_half_range_refuses_what_is_not_a_partition() -> None:
    r"""``HalfRange`` refuses a part off the ordinate set, an overlap and a non-cover; ``half_range`` refuses a non-finite cosine, a
    two-dimensional cosine array and a ``sign`` other than ``+1`` or ``-1``; the refusals are ``raise`` and survive ``python -O``."""
    from orpheus.numerics.declaration import Coord

    bad = H.partition_refusal_violations(_half_range_class(), Coord) + H.refusal_violations(_owner())
    assert not bad, f"refusals missing: {bad}"


def test_a_half_range_is_frozen_and_equal_partitions_are_one_value() -> None:
    half = _owner()(H.CosineStub([-0.5, 0.0, 0.5]), 0)
    with pytest.raises(dataclasses.FrozenInstanceError):
        half.n = 5  # type: ignore[misc]
    again = _owner()(H.CosineStub([-0.5, 0.0, 0.5]), 0)
    assert half == again and hash(half) == hash(again)
    other = _owner()(H.CosineStub([-0.5, 0.0, 0.5]), 0, sign=-1)
    assert half != other


def test_the_refusal_instrument_reads_an_unguarded_classifier_as_wrong() -> None:
    r"""Control of the refusal row: dropping the finite, rank and sign checks leaves all six bad inputs accepted."""
    bad = H.refusal_violations(H.unguarded_owner)
    assert len(bad) == 6, bad
    assert H.refusal_violations(H.exact_owner) == []


def test_the_refusal_instrument_reads_a_dropped_finite_check_as_wrong() -> None:
    r"""Mutation: ``drop the isfinite check`` leaves the three non-finite rows (and only those) accepted."""

    def no_finite(quadrature, axis_index, *, sign=+1):
        c = sign * np.asarray(quadrature.mean_axis_cosine(axis_index), dtype=np.float64)
        if sign not in (+1, -1):
            raise ValueError("sign")
        if c.ndim != 1:
            raise ValueError("cosines must be one-dimensional")
        return H.exact_owner(H.CosineStub(np.nan_to_num(c, nan=1.0, posinf=1.0, neginf=-1.0)), axis_index)

    bad = H.refusal_violations(no_finite)
    assert [b.split(":")[0] for b in bad] == ["nan", "+inf", "-inf"], bad


def test_the_partition_instrument_reads_a_dropped_cover_check_as_wrong() -> None:
    r"""Mutation: ``HalfRange`` without its cover check accepts the non-cover case (and only that case)."""
    from orpheus.numerics.declaration import Coord

    bad = H.partition_refusal_violations(H.no_cover_partition, Coord)
    assert [b.split(":")[0] for b in bad] == ["parts that do not cover"], bad


def test_the_extent_is_parsed_like_a_coord_extent() -> None:
    r"""``HalfRange(..., n=3.0)`` is refused as ``Coord`` refuses a float extent, and ``n=3`` is accepted (review finding: ``3.0 == 3`` hashed alike)."""
    from orpheus.numerics.declaration import Coord

    cls = _half_range_class()
    parts = (Coord(3, (0,)), Coord(3, (1,)), Coord(3, (2,)))
    assert cls(*parts, 3).n == 3
    with pytest.raises(TypeError):
        cls(*parts, 3.0)
    with pytest.raises(TypeError):
        cls(*parts, True)


# ── G02f ─────────────────────────────────────────────────────────────


def test_octant_labels_are_the_owners_signs() -> None:
    r"""``Quadrature.octants`` labels equal ``np.sign`` of the quadrature's own cosines, and the owner's ``sign()`` equals it too, on the 17 rules and the nudged rule.

    The expected value is read from ``mean_axis_cosine`` by ``np.sign`` against exact zero, so it shares no code with ``half_range``.

    The behavioural red comes first: on the nudged ``product(2, 4)`` the x label is ``+1`` on the 4 nudged ordinates; today it is ``0``.
    """
    nudged = H.nudged_product(9.5e-16)
    mu_x = np.asarray(nudged.mean_axis_cosine(0))
    nudge = np.flatnonzero((mu_x > 0) & (mu_x < 1e-10))
    assert nudge.size == 4, f"the fixture is not nudged: {nudge.tolist()}"
    assert H.octant_labels(nudged)[nudge, 0].tolist() == [1, 1, 1, 1], "the octant label reads a 9.5e-16 cosine as tangential"

    half_range = _owner()
    inputs = [("nudged product(2, 4)", nudged)] + [(name, make()) for name, make in H.census_rules()]
    for name, quad in inputs:
        dim = int(quad.dim)
        signs = [half_range(quad, i).sign() for i in range(dim)]
        for i, s in enumerate(signs):
            assert s.dtype == np.int8, f"{name} axis {i}: sign() dtype {s.dtype}, wanted int8"
            assert set(np.unique(s).tolist()) <= {-1, 0, 1}, f"{name} axis {i}: sign() takes values {np.unique(s).tolist()}"
            tangential = np.asarray(half_range(quad, i).tangential.indices, dtype=np.int64)
            assert np.array_equal(np.flatnonzero(s == 0), tangential), f"{name} axis {i}: sign() is 0 off the tangential set"
        independent = [np.sign(np.asarray(quad.mean_axis_cosine(i), dtype=np.float64)).astype(np.int8) for i in range(dim)]
        for i in range(dim):
            assert np.array_equal(signs[i], independent[i]), f"{name} axis {i}: the owner's sign() differs from np.sign of the cosines"
        want = independent[0] if dim == 1 else np.stack(independent, axis=1)
        got = H.octant_labels(quad)
        assert np.array_equal(got, want), f"{name}: octants labels differ from the owner's signs on {int(np.sum(got != want))} entries"


# ── G02g ─────────────────────────────────────────────────────────────


def test_no_selector_compares_a_cosine_array() -> None:
    r"""In the scopes S no ``Compare`` reads a cosine array, and the two trace selectors read ``half_ranges`` and compare nothing."""
    root = H.orpheus_root()
    hits = H.scope_cosine_compares(root)
    assert not hits, f"{len(hits)} Compare nodes read a cosine array in the scopes (file, line): {hits}"

    tree = ast.parse((root / "numerics/spaces/angular_trace_space.py").read_text())
    for name in ("inflow_indices_for_face", "outflow_indices_for_face"):
        fn = H.find_function(tree, f"AngularTraceSpace.{name}")
        assert fn is not None, f"AngularTraceSpace.{name} is not found"
        assert not [n for n in ast.walk(fn) if isinstance(n, ast.Compare)], f"{name} contains a Compare"
        assert any(isinstance(n, ast.Attribute) and n.attr == "half_ranges" for n in ast.walk(fn)), f"{name} does not read half_ranges"


def test_the_cosine_compare_predicate_reads_a_put_back_compare_on_synthetic_text() -> None:
    r"""Control of G02g: a selector body with ``np.flatnonzero(row < 0.0)`` reads one hit, the clean body none; a compare of ``mu_x``
    in a copy of ``_dag_legs`` reads; a rank question (``nodes.ndim == 1``) does not."""
    clean = "class A:\n    def inflow_indices_for_face(self, face):\n        row = self.half_ranges[0]\n        return row\n"
    assert H.cosine_compares(ast.parse(clean), ("inflow_indices_for_face",)) == []
    bad = clean.replace("return row", "return np.flatnonzero(row < 0.0)")
    assert H.cosine_compares(ast.parse(bad), ("inflow_indices_for_face",)) == [4]
    legs = "def _dag_legs(self):\n    return mu_x > 0.0\n"
    assert H.cosine_compares(ast.parse(legs), ("_dag_legs",)) == [2]
    assert H.cosine_compares(ast.parse("def g(nodes):\n    return nodes.ndim == 1\n")) == []
