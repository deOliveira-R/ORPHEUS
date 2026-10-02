r"""Mesh-free functions a specification states: a per-region table and a symbolic function.

A specification states its source and its detector before any mesh exists,
so neither can be an array on a discretised space. Two intensional forms
cover the cases (#405 P1 step 6):

* :class:`RegionwiseConstant`: one real value per (region, energy group), a
  function on the **angle-integrated** space. A region is a positional index
  of a geometry's intervals; a group is an index of the energy axis.
* :class:`Symbolic`: one SymPy expression :math:`q_g(r, \mu, \varphi)` per
  energy group, a function on **phase space**, stored as ``srepr`` text.

**The two arrows into phase space** (the user's ruling of 2026-10-02: the
:math:`4\pi` is the measure, not a convention). A per-region table lives on
the space :math:`R` maps to, where :math:`R\psi = \int \psi\,\mathrm d\Omega`
is the angular retraction (fibre integration). It enters phase space by one
of two arrows, and the ROLE the specification gives it picks the arrow:

* a source of rate :math:`Q` enters the right-hand side through the section
  :math:`E`, defined by :math:`R \circ E = \mathrm{id}`, so the rate is kept;
* a detector :math:`\Sigma_d` is the functional
  :math:`\psi \mapsto \langle \Sigma_d, R\psi\rangle`, whose Riesz
  representative is the adjoint :math:`R^\dagger \Sigma_d`, the pullback.

The two differ by :math:`R \circ R^\dagger`, the mass of the angular measure
(:math:`4\pi` on the sphere), which is therefore never written: neither type
carries a role or a density, and the field of the specification that holds
the value is the role. A :class:`Symbolic` is already a function on phase
space and takes no lift. The discrete arrows are
:class:`~orpheus.numerics.operator.AxisSectionOperator` and
:class:`~orpheus.numerics.operator.AxisPullbackOperator`; the continuous
ones are Branch 1's, :mod:`orpheus.derivations.common.angular_measure`.

**The chart.** :math:`\mu` and :math:`\varphi` are read in the chart the
geometry's coordinate system declares,
:attr:`orpheus.geometry.coord.CoordSystem.angular_chart`; this module is
free of geometry, so a :class:`Symbolic` is chart-relative until a
specification pairs it with a geometry.

**SymPy** is imported inside the functions that need it: importing this
module, or building a :class:`RegionwiseConstant`, loads no SymPy.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import TYPE_CHECKING, Any

import numpy as np

from orpheus.numerics.content import ContentIdentity

if TYPE_CHECKING:
    import sympy


@dataclass(frozen=True, eq=False)
class RegionwiseConstant(ContentIdentity):
    r"""One real value per (region, energy group): a function on the angle-integrated space.

    ``values`` has shape ``(regions, groups)``, at least one of each, and is
    stored as a read-only float copy, so the caller's array can change
    afterwards without moving the value or its digest. Every entry is finite:
    NaN is not a value, and an infinite rate is not a function value.

    The table carries no role and no density. As a source it enters phase
    space through the angular section :math:`E` (the rate is kept); as a
    detector, through the retraction's adjoint :math:`R^\dagger` (the module
    docstring).
    """

    values: np.ndarray

    def __post_init__(self) -> None:
        raw = np.asarray(self.values)
        if raw.dtype.kind not in "biuf" or raw.dtype.kind == "b":
            raise TypeError(
                f"RegionwiseConstant: values must be real numbers, got an array of dtype {raw.dtype}"
            )
        if raw.ndim != 2:
            raise ValueError(
                f"RegionwiseConstant: values is a (regions, groups) table of rank 2, got rank {raw.ndim}"
            )
        for axis, noun in enumerate(("regions", "groups")):
            if raw.shape[axis] == 0:
                raise ValueError(f"RegionwiseConstant: the table has no {noun} (shape {raw.shape})")
        table = np.array(raw, dtype=float) + 0.0
        for index in zip(*np.nonzero(np.isnan(table))):
            raise ValueError(f"RegionwiseConstant: the entry {tuple(int(i) for i in index)} is NaN, which is not a number")
        for index in zip(*np.nonzero(np.isinf(table))):
            raise ValueError(
                f"RegionwiseConstant: the entry {tuple(int(i) for i in index)} is infinite; "
                f"a rate or a response is a finite function value"
            )
        table.flags.writeable = False
        object.__setattr__(self, "values", table)

    @property
    def n_regions(self) -> int:
        """The number of regions, the table's row count."""
        return int(self.values.shape[0])

    @property
    def n_groups(self) -> int:
        """The number of energy groups, the table's column count."""
        return int(self.values.shape[1])


class _OwnedSymbol:
    """A coordinate symbol of :class:`Symbolic`, built on first access (no import-time SymPy)."""

    def __init__(self, name: str) -> None:
        self._name = name

    def __get__(self, instance: object, owner: type) -> "sympy.Symbol":
        return _owned_symbols()[self._name]


@cache
def _owned_symbols() -> dict[str, "sympy.Symbol"]:
    import sympy

    return {name: sympy.Symbol(name, real=True) for name in ("r", "mu", "phi")}


@dataclass(frozen=True, eq=False)
class Symbolic(ContentIdentity):
    r"""One SymPy expression :math:`q_g(r, \mu, \varphi)` per energy group, a function on phase space.

    **The coordinates** are the class's own symbols :attr:`r`, :attr:`mu`
    and :attr:`phi`, each ``Symbol(name, real=True)``: the position, the
    cosine to the chart's polar axis and the azimuth about it, in the chart
    the geometry's coordinate system declares
    (:attr:`orpheus.geometry.coord.CoordSystem.angular_chart`). An expression
    with any other free symbol is refused, naming the symbol and its
    assumptions; a symbol with an owned name but other assumptions
    (``Symbol("mu", positive=True)``) is a different symbol to SymPy and is
    refused under its own message.

    **Storage is text.** ``srepr`` holds one canonical ``sympy.srepr`` string
    per group, and the content is that text, so identity is by spelling:
    :math:`(r + 1)^2` and its expansion are two values (a cache miss, never a
    wrong hit). ``sympy_version`` is the version of the SymPy that wrote the
    text, and it is part of the content: a SymPy that parses the same text
    into a different expression would otherwise return a stored answer for
    a different function, so every SymPy upgrade invalidates every key
    instead. An expression holding ``nan``, ``zoo`` or an infinity is
    refused: it is not a function value.

    Build one from expressions with :meth:`of`, or from stored text with
    :meth:`from_srepr`.
    """

    srepr: tuple[str, ...]
    sympy_version: str

    r = _OwnedSymbol("r")
    mu = _OwnedSymbol("mu")
    phi = _OwnedSymbol("phi")

    def __post_init__(self) -> None:
        import sympy

        texts = tuple(self.srepr)
        if not texts:
            raise ValueError("Symbolic: a function has at least one energy group")
        owned = _owned_symbols()
        canonical = []
        for g, text in enumerate(texts):
            if not isinstance(text, str):
                raise TypeError(f"Symbolic: group {g} is stored as srepr text, got {type(text).__name__}")
            expression = sympy.sympify(text)
            _admit(expression, g, owned)
            canonical.append(sympy.srepr(expression))
        object.__setattr__(self, "srepr", tuple(canonical))
        if not isinstance(self.sympy_version, str):
            raise TypeError(f"Symbolic: sympy_version is a version string, got {type(self.sympy_version).__name__}")

    @classmethod
    def of(cls, *expressions: Any) -> "Symbolic":
        """One expression per group, written by this process's SymPy."""
        import sympy

        return cls(srepr=tuple(sympy.srepr(sympy.sympify(e)) for e in expressions), sympy_version=sympy.__version__)

    @classmethod
    def from_srepr(cls, srepr: tuple[str, ...], sympy_version: str | None = None) -> "Symbolic":
        """Stored text, with the version of the SymPy that wrote it (this process's when omitted)."""
        import sympy

        return cls(srepr=tuple(srepr), sympy_version=sympy.__version__ if sympy_version is None else sympy_version)

    @property
    def expressions(self) -> tuple["sympy.Basic", ...]:
        """The expressions, one per group, parsed from the stored text."""
        import sympy

        return tuple(sympy.sympify(text) for text in self.srepr)

    @property
    def n_groups(self) -> int:
        """The number of energy groups, one expression each."""
        return len(self.srepr)

    @property
    def is_isotropic(self) -> bool:
        r"""Whether no group depends on the direction.

        Isotropic iff, in every group, ``simplify`` reduces both
        :math:`\partial q_g/\partial\mu` and :math:`\partial q_g/\partial\varphi`
        to 0. A derivative ``simplify`` cannot reduce counts as a dependence:
        the undecided case falls on the anisotropic side, which is the side a
        consumer refusing anisotropy refuses (``sin(φ)**2 + cos(φ)**2`` is
        decided isotropic; a free-symbols test would call it anisotropic).
        """
        import sympy

        return all(
            sympy.simplify(sympy.diff(q, coordinate)) == 0
            for q in self.expressions
            for coordinate in (self.mu, self.phi)
        )


def _admit(expression: "sympy.Basic", group: int, owned: dict[str, "sympy.Symbol"]) -> None:
    """Refuse a stray coordinate or a non-value in one group's expression."""
    import sympy

    for symbol in sorted(expression.free_symbols, key=str):
        name = str(symbol)
        if owned.get(name) == symbol:
            continue
        assumptions = {k: v for k, v in sorted(symbol.assumptions0.items()) if v is not None}
        if name in owned:
            raise ValueError(
                f"Symbolic: group {group} uses a symbol {name!r} with assumptions {assumptions}: "
                f"same name, different assumptions: use Symbolic.{name} (real=True)"
            )
        raise ValueError(
            f"Symbolic: group {group} has the free symbol {name!r} with assumptions {assumptions}; "
            f"the coordinates are Symbolic.r, Symbolic.mu and Symbolic.phi only"
        )
    for non_value in (sympy.nan, sympy.zoo, sympy.oo, -sympy.oo):
        if expression.has(non_value):
            raise ValueError(f"Symbolic: group {group} contains {non_value}, which is not a function value")


__all__ = ["RegionwiseConstant", "Symbolic"]
