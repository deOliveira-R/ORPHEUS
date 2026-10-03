r"""Parsers for scalar and array inputs: one definition of "a real number".

One definition each of "a real number", "a finite real number", "an
integer" and "a positive quantity", shared by every layer that admits
numbers: the content encoder (:mod:`orpheus.numerics.content`), the
mesh-free functions and question values of ``orpheus.numerics``, the
geometry and the mesh.

The conversion every one of them makes is :func:`exact_double` (and
:func:`canonical_reals`, the same rule entrywise): an integer beyond
2**53 is refused (a double would round it, so two unequal integers would
become one value), a magnitude beyond the range of a double is refused,
NaN is refused (it is not a number, and not equal to itself), and ``-0.0``
becomes ``+0.0`` (:func:`canonical_real`; the two compare equal, so they
are one value). The stored bits and the digest therefore agree (#405 P1
step 5; the earlier copies were unified by #559). The parsers also refuse
``bool`` wherever a number is expected; the encoder admits it, because
``True == 1`` and its content follows ``==``.

The module imports nothing from ``orpheus``: it is a leaf, below the
content encoder that uses it.
"""
from __future__ import annotations

import math
from collections.abc import Iterable
from numbers import Real

import numpy as np


EXACT_INTEGER_LIMIT = 2**53
"""The largest integer magnitude a double carries exactly (every integer up to it)."""


def canonical_real(value: float, where: str) -> float:
    """The canonical ``float`` of one real value: NaN refused, ``-0.0`` made ``+0.0``.

    ``where`` names the value in the refusal. Infinities pass; a caller
    needing a finite value uses :func:`parse_finite_real`.
    """
    if math.isnan(value):
        raise ValueError(f"{where} is NaN, which is not a number")
    return value + 0.0  # -0.0 + 0.0 is +0.0


def exact_double(value: object, where: str) -> float:
    """The double that carries a real scalar exactly as ``==`` sees it, or a keyed refusal.

    An integer beyond :data:`EXACT_INTEGER_LIMIT` is refused (a double would
    round it, so two unequal integers would become one value), a magnitude
    beyond the range of a double is refused, and the result passes through
    :func:`canonical_real`. ``bool`` is admitted here, as ``True == 1``; the
    parsers that refuse it check the type first.
    """
    if isinstance(value, (int, np.integer)) and abs(int(value)) > EXACT_INTEGER_LIMIT:
        raise ValueError(
            f"{where}: the integer {int(value)} lies beyond 2**53, where a double cannot carry it exactly"
        )
    try:
        double = float(value)  # type: ignore[arg-type]
    except OverflowError:
        raise ValueError(f"{where}: {value!r} lies beyond the range of a double") from None
    return canonical_real(double, where)


def canonical_reals(values: np.ndarray, where: str) -> np.ndarray:
    """:func:`exact_double` entrywise: a ``float64`` copy of a boolean, integer or real array.

    A refused entry is named by its index; the ``-0.0`` fold is
    :func:`canonical_real`'s, applied to the whole array at once.
    """
    if values.dtype.kind in "iu" and values.size:
        for position in (int(np.argmin(values)), int(np.argmax(values))):
            index = np.unravel_index(position, values.shape)
            exact_double(int(values[index]), f"{where}, the entry {tuple(int(i) for i in index)}")
    doubles = np.array(values, dtype=np.float64)
    for index in zip(*np.nonzero(np.isnan(doubles))):
        canonical_real(math.nan, f"{where}, the entry {tuple(int(i) for i in index)}")
    return doubles + 0.0  # -0.0 + 0.0 is +0.0


def parse_real(value: object, where: str) -> float:
    """A real scalar as a ``float``, or a keyed refusal.

    NaN is refused: it is not a number, and not equal to itself, so a value
    holding one has no content identity (#405 P1 step 5). Infinities are
    real and pass; a caller needing a finite value checks it.
    """
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(f"{where} must be a real number, got {type(value).__name__}")
    return exact_double(value, where)


def parse_finite_real(value: object, where: str) -> float:
    """A real, finite scalar as a ``float``, or a keyed refusal."""
    parsed = parse_real(value, where)
    if math.isinf(parsed):
        raise ValueError(f"{where} is infinite, not a finite real number")
    return parsed


def parse_finite_reals(value: object, where: str) -> np.ndarray:
    """A real array of any rank as a read-only ``float`` copy, every entry finite.

    Every entry is parsed by :func:`parse_finite_real`, so an array obeys the
    scalar rule exactly: a ``bool`` or a non-real entry is refused, a NaN or
    infinite entry is refused naming its index, and ``-0.0`` becomes ``+0.0``.
    An ``ndarray`` whose dtype is not real (``bool`` included) is refused as a
    whole; an object array and a nested sequence are parsed entry by entry.
    The copy is taken, so the caller's array can change afterwards
    without moving the parsed one.
    """
    if isinstance(value, np.ndarray) and value.dtype.kind != "O" and value.dtype.kind not in "iuf":
        raise TypeError(f"{where} must be real numbers, got an array of dtype {value.dtype}")
    entries = np.asarray(value, dtype=object)
    array = np.empty(entries.shape, dtype=float)
    for index, entry in np.ndenumerate(entries):
        array[index] = parse_finite_real(entry, f"{where}, the entry {index}")
    array.flags.writeable = False
    return array


def parse_positive_real(value: object, where: str, noun: str) -> float:
    """A positive, finite real scalar as a ``float``, or a keyed refusal."""
    parsed = parse_real(value, where)
    if not (math.isfinite(parsed) and parsed > 0.0):
        raise ValueError(f"{where}: {noun} is positive and finite, got {value!r}")
    return parsed


def parse_integer(value: object, where: str, noun: str) -> int:
    """An integer scalar as an ``int``, or a keyed refusal."""
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError(f"{where}: {noun} is an int, got {type(value).__name__}")
    return int(value)


def parse_index(value: object, where: str, noun: str) -> int:
    """A non-negative integer scalar (an index) as an ``int``, or a keyed refusal."""
    parsed = parse_integer(value, where, noun)
    if parsed < 0:
        raise ValueError(f"{where}: {noun} is a non-negative index, got {parsed}")
    return parsed


def parse_member(value: object, kinds: tuple[type, ...], where: str, noun: str, what: str) -> object:
    """A value that is an instance of one of ``kinds`` (the members of a closed sum), or a keyed refusal.

    The refusal names the owner, the role and the sum, and lists its members:
    ``"FixedSource: the source is a mesh-free function (RegionwiseConstant or
    Symbolic), got a numpy.ndarray"``.
    """
    if not isinstance(value, kinds):
        names = " or ".join(kind.__name__ for kind in kinds)
        raise TypeError(
            f"{where}: {noun} is {what} ({names}), got a {type(value).__module__}.{type(value).__qualname__}"
        )
    return value


def parse_positive_integer(value: object, where: str, noun: str) -> int:
    """An integer scalar of at least 1 as an ``int``, or a keyed refusal."""
    parsed = parse_integer(value, where, noun)
    if parsed < 1:
        raise ValueError(f"{where}: {noun} is at least 1, got {parsed}")
    return parsed


def parse_entries(value: object, where: str, expected: str) -> tuple[object, ...]:
    """The entries of a sequence as a tuple, or a keyed refusal."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise TypeError(f"{where} must be a sequence of {expected}, got {type(value).__name__}")
    return tuple(value)


def parse_positions(value: object, where: str) -> np.ndarray:
    """A 1-D sequence of real, finite positions as a read-only float array."""
    entries = parse_entries(value, where, "real numbers")
    array = np.array(
        [parse_finite_real(v, f"{where}[{k}]") for k, v in enumerate(entries)], dtype=float,
    )
    array.flags.writeable = False
    return array


__all__ = [
    "EXACT_INTEGER_LIMIT",
    "canonical_real",
    "canonical_reals",
    "exact_double",
    "parse_entries",
    "parse_finite_real",
    "parse_index",
    "parse_member",
    "parse_finite_reals",
    "parse_integer",
    "parse_positions",
    "parse_positive_integer",
    "parse_positive_real",
    "parse_real",
]
