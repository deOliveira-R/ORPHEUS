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
the space carries, which is therefore never written: :math:`4\pi` for
:math:`\mathrm d\Omega` on the sphere, 2 on the orbit space of :math:`\mu`
alone, where a 1-D rule's ordinates live. Neither type carries a role or a
density, and the field of the specification that holds the value is the
role.

A :class:`Symbolic` is a density with respect to :math:`\mathrm d\Omega`,
so it needs no lift onto the sphere; onto a rule whose ordinates are points
of an orbit space (a 1-D rule) it is pushed forward, integrated over each
orbit (:math:`\int_0^{2\pi} q\,\mathrm d\varphi`, which is
:math:`2\pi q` when q does not depend on :math:`\varphi`). That
pushforward is part of projecting a specification onto a method's unknowns
(#405 phase P4), not of this module. The discrete arrows are
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
from typing import TYPE_CHECKING, Any, TypeAlias, cast, get_args

import numpy as np

from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.scalars import parse_member, parse_finite_reals

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
        table = parse_finite_reals(self.values, "RegionwiseConstant: values")
        if table.ndim != 2:
            raise ValueError(
                f"RegionwiseConstant: values is a (regions, groups) table of rank 2, got rank {table.ndim}"
            )
        for axis, noun in enumerate(("regions", "groups")):
            if table.shape[axis] == 0:
                raise ValueError(f"RegionwiseConstant: the table has no {noun} (shape {table.shape})")
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
    """A coordinate symbol of :class:`Symbolic`, named by its attribute and
    built on first access (no import-time SymPy)."""

    def __set_name__(self, owner: type, name: str) -> None:
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
    instead. Stored text is parsed through a whitelist (calls to SymPy
    classes, SymPy singletons, literals and keywords, nothing else), never
    by a bare ``eval``.

    **Only real scalar functions are admitted.** Refused, each keyed: an
    object that is not a scalar expression (a relation, a matrix); a
    ``nan``, ``zoo`` or infinity; the imaginary unit; an undefined function
    (``f(r)``); a ``Piecewise`` with no otherwise branch, which has no value
    outside its conditions.

    Build one from expressions with :meth:`of`, or from stored text with
    :meth:`from_srepr`.
    """

    srepr: tuple[str, ...]
    sympy_version: str

    r = _OwnedSymbol()
    mu = _OwnedSymbol()
    phi = _OwnedSymbol()

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
            expression = _parse(text, g)
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
    def expressions(self) -> tuple["sympy.Expr", ...]:
        """The expressions, one per group, parsed from the stored text (each a
        scalar ``Expr``: the constructor admitted nothing else)."""
        return tuple(cast("sympy.Expr", _parse(text, g)) for g, text in enumerate(self.srepr))

    @property
    def n_groups(self) -> int:
        """The number of energy groups, one expression each."""
        return len(self.srepr)

    def depends_on(self, *coordinates: "sympy.Symbol") -> bool:
        r"""Whether some group's value changes when the given coordinates do.

        ``coordinates`` are the owned symbols (:attr:`r`, :attr:`mu`,
        :attr:`phi`). The function does NOT depend on them iff, in every
        group, ``simplify`` reduces :math:`q_g - q_g|_{c \to c'}` to 0, with
        every named coordinate :math:`c` replaced by a fresh real symbol
        :math:`c'` at once. A difference ``simplify`` cannot reduce counts as
        a dependence: the undecided case falls on the dependent side, which is
        the side a consumer refusing the dependence refuses
        (``sin(φ)**2 + cos(φ)**2`` is decided independent of φ; a free-symbols
        test would call it dependent). A derivative test is wrong here: a
        step, ``Piecewise((1, μ > 0), (0, True))``, has zero derivative
        wherever it is defined.
        """
        import sympy

        owned = (self.r, self.mu, self.phi)
        if not coordinates or any(c not in owned for c in coordinates):
            raise ValueError(
                f"Symbolic.depends_on: name at least one owned coordinate (Symbolic.r, Symbolic.mu, "
                f"Symbolic.phi), got {coordinates!r}; a symbol the function does not own is not a coordinate of it"
            )
        fresh = [(c, sympy.Symbol(f"{c.name}_other", real=True)) for c in coordinates]
        return not all(sympy.simplify(q - q.subs(fresh, simultaneous=True)) == 0 for q in self.expressions)

    @property
    def is_isotropic(self) -> bool:
        """Whether no group depends on the direction (``μ`` and ``φ``, :meth:`depends_on`)."""
        return not self.depends_on(self.mu, self.phi)


@cache
def _sympy_names() -> dict[str, Any]:
    """The names ``srepr`` text may call: every subclass of ``sympy.Basic`` by
    its class name (``ExprCondPair`` is not a top-level export), and SymPy's
    top-level singletons (``pi``, ``true``, ``oo``); a top-level export wins
    a name two classes share."""
    import sympy

    names: dict[str, Any] = {}
    stack = [sympy.Basic]
    while stack:
        cls = stack.pop()
        names.setdefault(cls.__name__, cls)
        stack.extend(cls.__subclasses__())
    for name in dir(sympy):
        value = getattr(sympy, name)
        if (isinstance(value, type) and issubclass(value, sympy.Basic)) or isinstance(value, sympy.Basic):
            names[name] = value
    return names


def _parse(text: str, group: int) -> "sympy.Basic":
    """``srepr`` text as an expression, through a whitelist and never a bare ``eval``.

    Admitted: calls whose callee is a name of a SymPy class or a SymPy
    singleton (``pi``, ``true``, ``oo``), numeric and string literals, a
    unary sign on a literal, keywords, tuples. Anything else (an attribute,
    a subscript, a lambda, a name SymPy does not export as a class) is
    refused before evaluation, which runs with no builtins.
    """
    import ast

    import sympy

    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as err:
        raise ValueError(f"Symbolic: group {group} is not srepr text ({err.msg})") from None
    namespace: dict[str, Any] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Expression, ast.Call, ast.keyword, ast.Tuple, ast.Load, ast.USub, ast.UAdd)):
            continue
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, str, bool)):
            continue
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)) and isinstance(node.operand, ast.Constant):
            continue
        if isinstance(node, ast.Name):
            target = _sympy_names().get(node.id)
            if target is not None:
                namespace[node.id] = target
                continue
            raise ValueError(f"Symbolic: group {group} names {node.id!r}, which is not a SymPy class or constant")
        raise ValueError(f"Symbolic: group {group} holds a {type(node).__name__}, which srepr text never does")
    return eval(compile(tree, "<srepr>", "eval"), {"__builtins__": {}}, namespace)


def _admit(expression: "sympy.Basic", group: int, owned: dict[str, "sympy.Symbol"]) -> None:
    """Refuse anything in one group that is not a real scalar function of the owned coordinates."""
    import sympy
    from sympy.core.function import AppliedUndef

    if not isinstance(expression, sympy.Expr):
        raise ValueError(
            f"Symbolic: group {group} is a {type(expression).__name__}, not a scalar expression"
        )

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
    if expression.has(sympy.I):
        raise ValueError(f"Symbolic: group {group} contains the imaginary unit; a function value is real")
    undefined = sorted(str(f.func) for f in expression.atoms(AppliedUndef))
    if undefined:
        raise ValueError(f"Symbolic: group {group} applies the undefined function(s) {undefined}")
    for piecewise in expression.atoms(sympy.Piecewise):
        # The last (expression, condition) pair's condition is the otherwise.
        if piecewise.args[-1].args[1] != sympy.true:
            raise ValueError(
                f"Symbolic: group {group} has a Piecewise with no otherwise branch, "
                f"so it has no value outside its conditions"
            )


MeshFreeFunction: TypeAlias = RegionwiseConstant | Symbolic
"""A function on phase space stored without a mesh: the closed set of the two spellings."""



def parse_mesh_free_function(value: Any, where: str, noun: str) -> MeshFreeFunction:
    """A mesh-free function in the named role (a source, a detector, a weight), or a refusal naming its owner."""
    return cast(MeshFreeFunction, parse_member(value, get_args(MeshFreeFunction), where, noun, "a mesh-free function"))


__all__ = ["MeshFreeFunction", "RegionwiseConstant", "Symbolic", "parse_mesh_free_function"]
