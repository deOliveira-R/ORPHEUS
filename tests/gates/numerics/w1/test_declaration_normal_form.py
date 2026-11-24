r"""G01f: the declaration algebra: ``Coord``, ``Spanned``, ``Ann`` and the verbs ``contains`` and ``intersect`` (N01, spec 3.1, 2.3, 7.2).

Laws read here: ``Ann(Ann(d)) = d``; ``Ann(Coord(n, J)) = Coord(n, complement)``; ``d ∩ Ann(d) = Coord(n, ())`` for a ``Spanned`` ``d``;
``contains(outer, inner)`` (``inner ⊆ outer``, the outer operand first) is decided on ``Coord`` in ``Coord``, the empty ``Coord`` in
anything, anything in the whole ``Coord`` and ``d`` in ``d``, and RAISES elsewhere (it never guesses); ``intersect(a, b)`` is symmetric,
exact on ``Coord`` and on the two exact cases ``Coord(n, all) ∩ d = d`` and ``Coord(n, ()) ∩ d = Coord(n, ())``, and RAISES elsewhere
(it never returns ``None``). The refusals are pinned by their fragments (spec 3.1).

A second spelling of one subspace is bounded one-sidedly, not forbidden: two ``Spanned`` with different keys may name one subspace
(a key is a generator's content address, not a canonical form), so they compare UNEQUAL; ``test_one_subspace_two_keys_are_two_values``
records that boundary and carries no ``verifies``. The roster legs (digests) are ``test_content_identity_declaration.py``.
"""

from __future__ import annotations

import pytest

from orpheus.numerics.declaration import Ann, Coord, Declaration, Spanned, contains, intersect

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


N = 4
EMPTY = Coord(N, ())
WHOLE = Coord(N, (0, 1, 2, 3))
C13 = Coord(N, (1, 3))
C1 = Coord(N, (1,))
S = Spanned(N, ("k", 1), 2)
S2 = Spanned(N, ("k", 2), 2)
A = Ann(N, S)


def test_annihilator_normal_form() -> None:
    _require(S.ann() == A and isinstance(S.ann(), Ann), "Spanned.ann() is the Ann of it")
    _require(S.ann().ann() == S, "Ann(Ann(d)) = d for a Spanned")
    _require(A.ann() == S and A.dim == N - S.dim, "Ann.ann() is its Spanned and dim = n - dim(of)")
    _require(C13.ann() == Coord(N, (0, 2)) and isinstance(C13.ann(), Coord), "Coord.ann() is the Coord complement")
    _require(C13.ann().ann() == C13, "Coord.ann().ann() = Coord")
    _require(EMPTY.ann() == WHOLE and WHOLE.ann() == EMPTY, "the complement of the empty Coord is the whole Coord")
    _require(C13.dim == 2 and S.dim == 2, "dim reads the declaration")


def test_ann_is_built_only_from_a_spanned_of_the_same_n() -> None:
    with pytest.raises((ValueError, TypeError), match="annihilator of a Spanned"):
        Ann(N, C13)  # type: ignore[arg-type]
    with pytest.raises((ValueError, TypeError), match="annihilator of a Spanned"):
        Ann(N + 1, S)


def test_spanned_refusals() -> None:
    with pytest.raises(ValueError, match="a spanned subspace is not the zero subspace"):
        Spanned(N, ("k",), 0)
    with pytest.raises(ValueError, match="a spanned subspace is not the whole space"):
        Spanned(N, ("k",), N)
    with pytest.raises((ValueError, TypeError)):
        Spanned(N, (), 1)  # an empty key names nothing
    with pytest.raises((ValueError, TypeError)):
        Spanned(N, ([1],), 1)  # an unhashable key


def test_coord_refusals() -> None:
    for bad in ((2, 1), (1, 1), (0, 0, 1)):
        with pytest.raises(ValueError, match="indices must be sorted ascending and unique"):
            Coord(N, bad)
    with pytest.raises(ValueError):
        Coord(N, (0, N))
    with pytest.raises(ValueError):
        Coord(N, (-1, 0))
    with pytest.raises(ValueError):
        Coord(-1, ())


def test_equal_declarations_are_one_value() -> None:
    _require(Coord(N, (1, 3)) == C13 and hash(Coord(N, (1, 3))) == hash(C13), "two builds of one Coord are one value")
    _require(Spanned(N, ("k", 1), 2) == S and hash(Spanned(N, ("k", 1), 2)) == hash(S), "two builds of one Spanned are one value")
    _require(Ann(N, Spanned(N, ("k", 1), 2)) == A, "two builds of one Ann are one value")
    _require(C13 != Coord(N, (1, 2)) and S != Spanned(N, ("k", 1), 1), "a moved part is another value")


def test_one_subspace_two_keys_are_two_values() -> None:
    """The one-sided boundary of 3.1: a key is a generator's address, so two keys over one subspace are NOT identified."""
    _require(S != S2, "two Spanned of different keys must compare unequal: the declaration cannot decide they name one subspace")


def test_contains_decides_exactly_where_it_can() -> None:
    _require(contains(C13, C1) is True, "Coord in Coord: {1} in {1, 3}")
    _require(contains(C1, C13) is False, "Coord in Coord: {1, 3} not in {1}")
    _require(contains(S, EMPTY) is True and contains(A, EMPTY) is True, "the empty Coord is in anything")
    _require(contains(WHOLE, S) is True and contains(WHOLE, A) is True, "anything is in the whole Coord")
    _require(contains(S, S) is True and contains(A, A) is True and contains(C13, C13) is True, "d in d")
    _require(contains(EMPTY, EMPTY) is True and contains(WHOLE, WHOLE) is True, "the corners")
    _require(contains(EMPTY, C1) is False, "a non-empty Coord is not in the empty Coord")


@pytest.mark.parametrize(
    "outer,inner",
    [(S, S2), (S, A), (A, S), (S, C1), (C13, S), (A, C13), (C1, A)],
    ids=["spanned_spanned", "spanned_ann", "ann_spanned", "spanned_coord", "coord_spanned", "ann_coord", "coord_ann"],
)
def test_contains_raises_where_it_cannot_decide(outer: Declaration, inner: Declaration) -> None:
    with pytest.raises(ValueError, match="subspace inclusion is undecided"):
        contains(outer, inner)


def test_intersect_is_exact_and_symmetric() -> None:
    _require(intersect(C13, Coord(N, (0, 1))) == C1, "Coord ∩ Coord is the common indices")
    _require(intersect(Coord(N, (0, 1)), C13) == C1, "symmetric on Coord")
    _require(intersect(C1, Coord(N, (0, 2))) == EMPTY, "disjoint Coords meet in the empty Coord")
    for d in (S, A, C13):
        _require(intersect(WHOLE, d) == d and intersect(d, WHOLE) == d, f"Coord(n, all) ∩ {d} = {d}, both orders")
        _require(intersect(EMPTY, d) == EMPTY and intersect(d, EMPTY) == EMPTY, f"Coord(n, ()) ∩ {d} = the empty Coord, both orders")


def test_a_spanned_meets_its_annihilator_in_the_zero_subspace() -> None:
    _require(intersect(S, S.ann()) == EMPTY and intersect(S.ann(), S) == EMPTY, "d ∩ Ann(d) = Coord(n, ()) for a Spanned d, both orders")


@pytest.mark.parametrize(
    "a,b",
    [(S, S2), (S, C13), (C13, S), (A, C13), (S, Ann(N, S2)), (C1, A)],
    ids=["spanned_spanned", "spanned_coord", "coord_spanned", "ann_coord", "spanned_ann_other", "coord_ann"],
)
def test_intersect_raises_where_it_cannot_decide(a: Declaration, b: Declaration) -> None:
    with pytest.raises(ValueError, match="subspace intersection is undecided"):
        intersect(a, b)
