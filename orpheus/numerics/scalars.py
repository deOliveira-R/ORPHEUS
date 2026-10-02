r"""Parsers for scalar and array inputs: one definition of "a real number".

One definition each of "a real number", "a finite real number", "an
integer" and "a positive quantity", shared by every layer that admits
numbers: the content encoder (:mod:`orpheus.numerics.content`), the
mesh-free functions and question values of :mod:`orpheus.numerics`, the
geometry and the mesh. ``bool`` is refused wherever a number is expected
(``True`` is an ``int``), NaN is refused (it is not a number, and not equal
to itself), and ``-0.0`` becomes ``+0.0``: the two compare equal, so they
are one value. :func:`canonical_real` is that rule on a ``float``; every
parser here and the encoder call it, so the stored bits and the digest
agree (#405 P1 step 5; the three earlier copies were unified by #559).

The module imports nothing from ``orpheus``: it is a leaf, below the
content encoder that uses it.
"""
from __future__ import annotations

import math
from collections.abc import Iterable
from numbers import Real

import numpy as np


def canonical_real(value: float, where: str) -> float:
    """The canonical ``float`` of one real value: NaN refused, ``-0.0`` made ``+0.0``.

    ``where`` names the value in the refusal. Infinities pass; a caller
    needing a finite value uses :func:`parse_finite_real`.
    """
    if math.isnan(value):
        raise ValueError(f"{where} is NaN, which is not a number")
    return value + 0.0  # -0.0 + 0.0 is +0.0


def parse_real(value: object, where: str) -> float:
    """A real scalar as a ``float``, or a keyed refusal.

    NaN is refused: it is not a number, and not equal to itself, so a value
    holding one has no content identity (#405 P1 step 5). Infinities are
    real and pass; a caller needing a finite value checks it.
    """
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(f"{where} must be a real number, got {type(value).__name__}")
    return canonical_real(float(value), where)


def parse_finite_real(value: object, where: str) -> float:
    """A real, finite scalar as a ``float``, or a keyed refusal."""
    parsed = parse_real(value, where)
    if math.isinf(parsed):
        raise ValueError(f"{where} is infinite, not a finite real number")
    return parsed


def parse_finite_reals(value: object, where: str) -> np.ndarray:
    """A real array of any rank as a read-only ``float`` copy, every entry finite.

    The rule of :func:`parse_finite_real`, applied entrywise: a non-real
    dtype (``bool`` included) is refused, and a NaN or infinite entry is
    refused naming its index. The copy is taken, so the caller's array can
    change afterwards without moving the parsed one.
    """
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf":
        raise TypeError(f"{where} must be real numbers, got an array of dtype {raw.dtype}")
    array = np.array(raw, dtype=float) + 0.0  # -0.0 + 0.0 is +0.0
    for index in zip(*np.nonzero(~np.isfinite(array))):
        parse_finite_real(array[index], f"{where}, the entry {tuple(int(i) for i in index)}")
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
        [parse_real(v, f"{where}[{k}]") for k, v in enumerate(entries)], dtype=float,
    )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{where} must be finite; got {array}")
    array.flags.writeable = False
    return array


__all__ = [
    "canonical_real",
    "parse_entries",
    "parse_finite_real",
    "parse_finite_reals",
    "parse_integer",
    "parse_positions",
    "parse_positive_integer",
    "parse_positive_real",
    "parse_real",
]
