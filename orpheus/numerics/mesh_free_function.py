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
from typing import TYPE_CHECKING, Any, TypeAlias, assert_never, cast, get_args

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
            expression = parse_srepr(text, f"Symbolic: group {g}")
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
        return tuple(cast("sympy.Expr", parse_srepr(text, f"Symbolic: group {g}")) for g, text in enumerate(self.srepr))

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

    def without(self, *coordinates: "sympy.Symbol") -> "Symbolic":
        """The same function with the named coordinates eliminated, refused if it depends on one.

        The one definition of "independent of a coordinate" is
        :meth:`depends_on`; this is its constructive face: any named coordinate
        in an expression (a dependence that cancels, ``sin(φ)**2 + cos(φ)**2``)
        is set to 0, which changes no value because the function does not
        depend on it. The expression is NOT simplified first: ``simplify``
        rewrote ``Piecewise((1, sin(3r) > 0), (0, True))`` as the first period
        alone, ``0 < r < π/3``, dropping every later one (qa of #405 P2 step
        7b.2.1, `[M]` the value at r = 2.5 went from 1 to 0). A substitution
        alone is the identity on the function's values. A reader that needs
        the function on a space without those coordinates (the infinite
        medium has none) reads this instead of re-deciding constancy.
        """
        if self.depends_on(*coordinates):
            names = ", ".join(c.name for c in coordinates)
            raise ValueError(f"Symbolic.without: the function depends on {names}, so it cannot be read without them")
        eliminated = (q.subs({c: 0 for c in coordinates}, simultaneous=True) for q in self.expressions)
        return Symbolic.of(*eliminated)

    def steps(self, r_range: tuple[float, float]) -> tuple[float, ...]:
        r"""Where some group's value is not smooth in :math:`r` inside ``r_range``: sorted, each once.

        The one definition of a weight's step locations, read by every
        integrator that must split there (the mesh's exact cell integrals,
        :meth:`orpheus.mesh.structured.Mesh1D.cell_integrals`; a reference's
        float quadrature, #405 P2 step 7b.2.2). A step of a non-smooth
        construct sits where its argument changes sign: a ``Piecewise``
        condition's ``lhs − rhs``, the argument of ``Heaviside``, ``Abs`` and
        ``sign``, and the pairwise differences of ``Max`` and ``Min``'s
        arguments. Its locations are the real roots of that argument in the
        closed interval ``r_range``.

        SCOPE-BOUNDARY[guard] machinery: a root finder for a transcendental step argument, and the jump set of any other non-smooth head.
        ruling: the orchestrator, #405 P2 steps 7b.2.1 and 7b.2.2 (qa: `[M]` SymPy 1.14 integrates a sin step wrongly).
        revisit: when a weight with a transcendental step, or another non-smooth head, is needed.
        An ALLOW-list (elegance review of 7b.2.2, S1): every node of the
        expression is a smooth head (:func:`_is_smooth_head`) or a step
        construct whose locations this finds; any other head (``floor``,
        ``arg``, ``atan2``, ``Contains``, ...) is refused by name, since a
        reader that missed its jumps would integrate across them silently. A
        step is located only where its argument is polynomial in :math:`r`
        (``Piecewise((1, sin(3r) > 0), (0, True))`` is refused), with real
        coefficients of any kind (``r < pi/4``: exact roots where SymPy finds
        them, else 30-digit numerical roots). A step whose location depends on
        the direction (an argument in :math:`\mu` or :math:`\varphi`) is
        refused as well.
        """
        import sympy

        r = self.r
        low, high = (float(bound) for bound in r_range)
        locations: set[float] = set()
        for g, expression in enumerate(self.expressions):
            for node in sympy.preorder_traversal(expression):
                if not (_is_smooth_head(node) or isinstance(node, _step_heads())):
                    raise ValueError(
                        f"Symbolic.steps: group {g} has {node} (head {type(node).__name__}), which is neither smooth nor a "
                        f"step located here, so it cannot be integrated over the cells or split for a quadrature "
                        f"(a scope boundary)"
                    )
            for argument in _step_arguments(expression):
                if not argument.is_polynomial(r):
                    raise ValueError(
                        f"Symbolic.steps: group {g} steps where {argument} changes sign; a step is integrated "
                        f"only where its argument is polynomial in r (a scope boundary of the exact integration)"
                    )
                if argument.free_symbols - {r}:
                    raise ValueError(
                        f"Symbolic.steps: group {g} steps where {argument} changes sign, a location that "
                        f"depends on the direction; a step is located in r alone"
                    )
                locations |= _real_roots(argument, r)
        return tuple(sorted(x for x in locations if low <= x <= high))

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


def parse_srepr(text: str, where: str) -> "sympy.Basic":
    """``srepr`` text as an expression, through a whitelist and never a bare ``eval``; refusals begin with ``where``.

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
        raise ValueError(f"{where} is not srepr text ({err.msg})") from None
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
            raise ValueError(f"{where} names {node.id!r}, which is not a SymPy class or constant")
        raise ValueError(f"{where} holds a {type(node).__name__}, which srepr text never does")
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


def _is_smooth_head(node: "sympy.Basic") -> bool:
    """Whether ``node`` is smooth wherever it is defined: a number, a symbol, a sum or product, a power with a
    rational exponent (or a positive constant base), ``exp``, ``log``, a trigonometric or hyperbolic function."""
    import sympy
    from sympy.functions.elementary.hyperbolic import HyperbolicFunction
    from sympy.functions.elementary.trigonometric import TrigonometricFunction

    if isinstance(node, sympy.Pow):
        base, exponent = node.args
        return bool(exponent.is_Rational) or bool(base.is_number and base.is_positive)
    return isinstance(node, (sympy.Number, sympy.NumberSymbol, sympy.Symbol, sympy.Add, sympy.Mul, sympy.exp,
                             sympy.log, TrigonometricFunction, HyperbolicFunction))


@cache
def _step_heads() -> tuple[type, ...]:
    """The step constructs :meth:`Symbolic.steps` locates, and the logic a ``Piecewise`` condition is built from."""
    import sympy
    from sympy.functions.elementary.piecewise import ExprCondPair
    from sympy.logic.boolalg import BooleanAtom

    return (sympy.Piecewise, ExprCondPair, sympy.Heaviside, sympy.Abs, sympy.sign, sympy.Max, sympy.Min,
            sympy.core.relational.Relational, sympy.And, sympy.Or, sympy.Not, BooleanAtom)


def _real_roots(argument: "sympy.Expr", r: "sympy.Symbol") -> set[float]:
    """The real roots of a polynomial in ``r`` with real coefficients: exact where SymPy solves it, else numerical."""
    import sympy

    exact = sympy.solveset(argument, r, sympy.Reals)
    if isinstance(exact, sympy.FiniteSet):
        return {float(cast("sympy.Expr", root)) for root in exact}
    if exact is sympy.S.EmptySet or exact == sympy.S.Reals:  # no root, or an identically zero argument (no step)
        return set()
    numerical = cast("list[sympy.Expr]", sympy.Poly(argument, r).nroots(n=30))
    return {float(root) for root in numerical if root.is_real}


def _step_arguments(expression: "sympy.Expr") -> list["sympy.Expr"]:
    """The expressions whose sign changes are the steps of ``expression``'s non-smooth constructs."""
    import sympy

    arguments: list[sympy.Expr] = []
    for piecewise in expression.atoms(sympy.Piecewise):
        for branch in piecewise.args:
            for relation in branch.args[1].atoms(sympy.core.relational.Relational):
                arguments.append(cast("sympy.Expr", relation.lhs) - cast("sympy.Expr", relation.rhs))
    for kind in (sympy.Heaviside, sympy.Abs, sympy.sign):
        for atom in expression.atoms(kind):
            arguments.append(cast("sympy.Expr", atom.args[0]))
    for kind in (sympy.Max, sympy.Min):
        for atom in expression.atoms(kind):
            operands = [cast("sympy.Expr", operand) for operand in atom.args]
            arguments += [a - b for i, a in enumerate(operands) for b in operands[i + 1:]]
    return arguments


MeshFreeFunction: TypeAlias = RegionwiseConstant | Symbolic
"""A function on phase space stored without a mesh: the closed set of the two spellings."""



def parse_mesh_free_function(value: Any, where: str, noun: str) -> MeshFreeFunction:
    """A mesh-free function in the named role (a source, a detector, a weight), or a refusal naming its owner."""
    return cast(MeshFreeFunction, parse_member(value, get_args(MeshFreeFunction), where, noun, "a mesh-free function"))


def values_without_position(function: MeshFreeFunction, n_groups: int) -> tuple["sympy.Expr", ...]:
    """The function's value per group on a medium with no position (the infinite medium, a 0-D answer).

    The one reading of a mesh-free function where there is no coordinate to
    read it at: a regionwise-constant table must have exactly one region and
    is read EXACTLY (every double is a dyadic rational); a symbolic function is
    read :meth:`Symbolic.without` every owned coordinate, the same
    independence its admission decided. The group count must be ``n_groups``.
    Production reads these values as floats, a reference keeps them exact.
    """
    import sympy

    if function.n_groups != n_groups:
        raise ValueError(f"the function has {function.n_groups} groups; the medium has {n_groups}")
    match function:
        case RegionwiseConstant():
            if function.n_regions != 1:
                raise ValueError(f"the function has {function.n_regions} regions; a medium with no position has one")
            return tuple(sympy.Rational(float(value)) for value in function.values[0])
        case Symbolic():
            return function.without(Symbolic.r, Symbolic.mu, Symbolic.phi).expressions
        case _:
            assert_never(function)


__all__ = ["MeshFreeFunction", "RegionwiseConstant", "Symbolic", "parse_mesh_free_function", "parse_srepr", "values_without_position"]
