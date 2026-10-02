r"""Parsers for the scalar and sequence inputs of the geometry and the mesh.

One definition each of "a real number", "an integer" and "a positive
quantity", shared by :mod:`orpheus.geometry` and :mod:`orpheus.mesh` so the
two layers canonicalise the same values the same way. ``bool`` is
refused wherever a number is expected (``True`` is an ``int``), NaN is
refused (it is not a number), and ``-0.0`` becomes ``+0.0``: the two
compare equal, so they are one position. The content encoder
(:mod:`orpheus.numerics.content`, #405 P1 step 5) applies the same
``-0.0`` rule to every value it digests, so the stored bits and the digest
agree.
"""
from __future__ import annotations

import math
from collections.abc import Iterable
from numbers import Real

import numpy as np


def parse_real(value: object, where: str) -> float:
    """A real scalar as a ``float``, or a keyed refusal.

    NaN is refused: it is not a number, and not equal to itself, so a value
    holding one has no content identity (#405 P1 step 5). Infinities are
    real and pass; a caller needing a finite value checks it.
    """
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(f"{where} must be a real number, got {type(value).__name__}")
    parsed = float(value)
    if math.isnan(parsed):
        raise ValueError(f"{where} is NaN, which is not a number")
    return parsed + 0.0


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
    "parse_entries",
    "parse_integer",
    "parse_positions",
    "parse_positive_integer",
    "parse_positive_real",
    "parse_real",
]
