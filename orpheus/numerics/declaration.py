r"""Declarations of a subspace of :math:`\mathbb{R}^n`: the language in which a form states its kernel.

A *declaration* names a subspace of the index space :math:`\mathbb{R}^n` without carrying a numerical basis. There are three
spellings, closed under the annihilator:

* :class:`Coord` is an exact index subset: the coordinate subspace spanned by the unit vectors of the listed indices. The empty
  ``Coord`` is the zero subspace and the ``Coord`` of every index is the whole space, each with one spelling.
* :class:`Spanned` is a subspace of a given dimension that a generator knows by a content key (for example the digest of the Gram
  whose kernel it is). The numerical basis is never stored in the declaration, so two keys over one subspace are two values: a key
  is a generator's content address, not a canonical form.
* :class:`Ann` is the annihilator of a :class:`Spanned`; it is built through :meth:`Spanned.ann` only.

The annihilator is an involution: ``Ann(Ann(d)) = d``, and the annihilator of ``Coord(n, J)`` is the ``Coord`` of the complement of
``J``. :func:`contains` and :func:`intersect` decide inclusion and intersection exactly where the arms allow it and refuse
everywhere else, so a rule that needs an undecidable answer must be given the result by the generator that knows it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from orpheus.values.content import ContentIdentity

__all__ = ["Coord", "Spanned", "Ann", "Declaration", "contains", "intersect"]


def _parse_extent(name: str, value: object) -> int:
    r"""``value`` as a plain ``int`` (a ``bool`` or a non-integer is refused)."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an int, got {type(value).__name__}")
    return value


@dataclass(frozen=True, eq=False)
class Coord(ContentIdentity):  # concept: declaration | role: arm | kind: set
    r"""The coordinate subspace of :math:`\mathbb{R}^n` spanned by the unit vectors of ``indices``.

    ``indices`` is sorted ascending, unique and inside ``[0, n)``.
    """

    n: int
    indices: tuple[int, ...]

    def __post_init__(self) -> None:
        n = _parse_extent("Coord.n", self.n)
        if n < 0:
            raise ValueError(f"Coord needs n >= 0, got {n}")
        idx = tuple(_parse_extent("a Coord index", i) for i in self.indices)
        if any(b <= a for a, b in zip(idx, idx[1:])):
            raise ValueError(f"indices must be sorted ascending and unique, got {idx}")
        if idx and not (0 <= idx[0] and idx[-1] < n):
            raise ValueError(f"indices must lie in [0, {n}), got {idx}")
        object.__setattr__(self, "indices", idx)

    @property
    def dim(self) -> int:
        r"""The dimension of the subspace: the number of listed indices."""
        return len(self.indices)

    def ann(self) -> "Coord":
        r"""The annihilator: the ``Coord`` of the complementary indices."""
        listed = set(self.indices)
        return Coord(self.n, tuple(i for i in range(self.n) if i not in listed))


@dataclass(frozen=True, eq=False)
class Spanned(ContentIdentity):  # concept: declaration | role: arm | kind: set
    r"""A subspace of :math:`\mathbb{R}^n` of dimension ``dim`` that a generator knows by the content ``key``.

    The zero subspace is ``Coord(n, ())`` and the whole space is the ``Coord`` of every index, so ``1 <= dim <= n - 1``.
    """

    n: int
    key: tuple
    dim: int

    def __post_init__(self) -> None:
        n = _parse_extent("Spanned.n", self.n)
        dim = _parse_extent("Spanned.dim", self.dim)
        if not isinstance(self.key, tuple):
            raise TypeError(f"a Spanned key must be a tuple, got {type(self.key).__name__}")
        if not self.key:
            raise ValueError("a spanned subspace needs a non-empty content key")
        hash(self.key)
        if dim < 1:
            raise ValueError(f"a spanned subspace is not the zero subspace (dim = {dim}); the zero subspace is Coord(n, ())")
        if dim > n - 1:
            raise ValueError(f"a spanned subspace is not the whole space (dim = {dim}, n = {n}); the whole space is the Coord of every index")

    def ann(self) -> "Ann":
        r"""The annihilator, an :class:`Ann` of this subspace."""
        return Ann(self.n, self)


@dataclass(frozen=True, eq=False)
class Ann(ContentIdentity):  # concept: declaration | role: arm | kind: set
    r"""The annihilator of the :class:`Spanned` ``of``, in :math:`\mathbb{R}^n`.

    Build it through :meth:`Spanned.ann`: the annihilator of a ``Coord`` is a ``Coord`` and an ``Ann`` of one is unspellable.
    """

    n: int
    of: Spanned

    def __post_init__(self) -> None:
        n = _parse_extent("Ann.n", self.n)
        if not isinstance(self.of, Spanned) or self.of.n != n:
            raise ValueError(f"Ann is the annihilator of a Spanned of the same n ({n}), got {self.of!r}")

    @property
    def dim(self) -> int:
        r"""The dimension of the annihilator: ``n`` minus the dimension of ``of``."""
        return self.n - self.of.dim

    def ann(self) -> Spanned:
        r"""The annihilator of the annihilator: the :class:`Spanned` itself."""
        return self.of


Declaration: TypeAlias = Coord | Spanned | Ann
r"""A declared subspace: an exact index subset, a keyed subspace, or the annihilator of a keyed subspace."""


def _same_n(outer: Declaration, inner: Declaration, verb: str) -> None:
    if outer.n != inner.n:
        raise ValueError(f"{verb}: the declarations live in R^{outer.n} and R^{inner.n}")


def _is_whole(d: Declaration) -> bool:
    return isinstance(d, Coord) and d.dim == d.n


def contains(outer: Declaration, inner: Declaration) -> bool:
    r"""Whether ``inner`` is a subspace of ``outer``.

    Decided exactly when both are ``Coord``, or ``inner`` is the empty ``Coord``, or ``outer`` is the ``Coord`` of every index, or
    the two declarations are equal, or ``inner`` has a larger dimension than ``outer`` (then it is not inside: ``False``, exact
    because no basis is needed). Anywhere else the declarations cannot decide it and the call raises ``ValueError``
    ("subspace inclusion is undecided"): the generator that knows the answer must declare it.
    """
    _same_n(outer, inner, "subspace inclusion")
    if isinstance(inner, Coord) and inner.dim == 0:
        return True
    if _is_whole(outer):
        return True
    if inner.dim > outer.dim:
        return False
    if outer == inner:
        return True
    if isinstance(outer, Coord) and isinstance(inner, Coord):
        return set(inner.indices) <= set(outer.indices)
    raise ValueError(f"subspace inclusion is undecided: {inner!r} in {outer!r}")


def intersect(a: Declaration, b: Declaration) -> Declaration:
    r"""The intersection of two declarations, symmetric in its operands.

    Exact on two ``Coord``, on the whole space (``Coord`` of every index) with anything, on the empty ``Coord`` with anything, and
    on a ``Spanned`` with its own annihilator (the empty ``Coord``). Anywhere else it raises ``ValueError`` ("subspace
    intersection is undecided"); it never returns ``None``.
    """
    _same_n(a, b, "subspace intersection")
    if isinstance(a, Coord) and isinstance(b, Coord):
        return Coord(a.n, tuple(sorted(set(a.indices) & set(b.indices))))
    for x, y in ((a, b), (b, a)):
        if _is_whole(x):
            return y
        if isinstance(x, Coord) and x.dim == 0:
            return x
    for x, y in ((a, b), (b, a)):
        if isinstance(x, Spanned) and isinstance(y, Ann) and y.of == x:
            return Coord(a.n, ())
    # SCOPE-BOUNDARY[guard]
    # machinery: the meet of a coordinate subspace with a keyed subspace needs the key's basis, which the generator owns (W8a)
    # ruling: SD-25 and scratch/ontology/w1/open/N03/rulings.md (R5), the orchestrator, 2026-11-24
    # revisit: W8a, which types the generators
    raise ValueError(f"subspace intersection is undecided: {a!r} and {b!r}")
