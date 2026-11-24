r"""Inputs and instruments shared by the half-range gates of unit N02 (``tests/gates/numerics/w1`` and ``tests/gates/mesh/w1``).

Every function here RETURNS data (index tuples, lists of violation messages, parse results) and asserts nothing, because an
``assert`` in a helper imported by a test is stripped under ``python -O``; the collected test files assert.

Three kinds of object live here:

* the inputs: :class:`CosineStub` (a quadrature that exposes only ``mean_axis_cosine``, the ``AxisCosineSource`` protocol of
  ``numerics/tangential.py``), the 17 shipped rules of the census, and the nudged ``product(2, 4)``;
* the instruments: functions that read an owner (a callable with the signature of ``half_range``) and return what it gets wrong,
  and the ``ast`` predicates of the call, cosine-compare and cutoff arms;
* the mutants: owners that are deliberately wrong (a tolerance classifier, a classifier that reads a non-finite cosine). A mutant
  returns light namespaces, not ``HalfRange``, so that its instrument runs on any tree, before the owner exists: each mutant is the
  positive control of the instrument that reads it, and every instrument is run on its mutant in the collected tests.

Nothing here classifies a cosine for production: the mutants are controls and are built from a tolerance that the owner must not
have.
"""

from __future__ import annotations

import ast
import re
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Callable

import numpy as np

EPS = float(np.finfo(np.float64).eps)

# ─────────────────────────────────────────────────────────────────────
# Inputs
# ─────────────────────────────────────────────────────────────────────


class CosineStub:
    r"""A quadrature that exposes only ``mean_axis_cosine``: the ``AxisCosineSource`` protocol and nothing else."""

    def __init__(self, cosines) -> None:
        self._cosines = np.asarray(cosines, dtype=np.float64)

    def mean_axis_cosine(self, axis_index: int) -> np.ndarray:
        return self._cosines


def census_rules() -> list[tuple[str, Callable[[], object]]]:
    r"""The 17 shipped rules of the census: ``(name, factory)``. Each is read at axes 0, 1, 2, so 51 (rule, axis) pairs."""
    from orpheus.numerics.quadrature import Quadrature as Q

    rules: list[tuple[str, Callable[[], object]]] = []
    rules += [(f"gl{n}", lambda n=n: Q.gauss_legendre(n)) for n in (2, 4, 8, 16)]
    rules += [(f"leb{o}", lambda o=o: Q.lebedev(o)) for o in (5, 17)]
    rules += [(f"ls{s}", lambda s=s: Q.level_symmetric(s)) for s in (2, 4, 6, 8)]
    rules += [(f"prod{a}x{b}", lambda a=a, b=b: Q.product(a, b)) for a, b in ((4, 4), (8, 8), (4, 3), (4, 5))]
    rules += [(f"fold{a}x{b}", lambda a=a, b=b: Q.folded_product(a, b)) for a, b in ((2, 4), (4, 8), (8, 8))]
    return rules


#: The rules of the qa probe ``q1`` population (a subset of the census): 11 rules, 33 pairs, 29 not identically zero.
Q1_RULES = ("gl4", "gl16", "leb5", "leb17", "ls4", "ls8", "prod4x4", "prod8x8", "prod4x3", "prod4x5", "fold4x8")


def nudged_product(value: float):
    r"""The shipped ``product(2, 4)`` with its exactly-zero ``mu_x`` replaced by ``value``, through the measure's nodes."""
    from orpheus.numerics.quadrature import Quadrature

    shipped = Quadrature.product(n_mu=2, n_phi=4)
    nodes = np.array(shipped.measure.nodes, copy=True)
    nodes[nodes[:, 0] == 0.0, 0] = value
    return replace(shipped, measure=replace(shipped.measure, nodes=nodes))


def octant_labels(quad) -> np.ndarray:
    r"""The octant label of every ordinate, ``(N,)`` for a one-dimensional rule and ``(N, dim)`` otherwise, read off ``Quadrature.octants``."""
    n = int(quad.n_ordinates)
    width = len(quad.octants[0].label)
    out = np.zeros((n, width), dtype=np.int64)
    for entry in quad.octants:
        out[np.asarray(entry.indices, dtype=np.int64)] = np.asarray(entry.label, dtype=np.int64)
    return out[:, 0] if int(quad.dim) == 1 else out


# ─────────────────────────────────────────────────────────────────────
# Mutants: owners that are wrong on purpose
# ─────────────────────────────────────────────────────────────────────


def _parts(n: int, out, inn, tan) -> SimpleNamespace:
    def part(idx) -> SimpleNamespace:
        return SimpleNamespace(n=n, indices=tuple(int(i) for i in idx))

    def sign() -> np.ndarray:
        s = np.zeros(n, dtype=np.int8)
        s[list(part(out).indices)] = 1
        s[list(part(inn).indices)] = -1
        return s

    return SimpleNamespace(n=n, outflow=part(out), inflow=part(inn), tangential=part(tan), sign=sign)


def banded_owner(band: float) -> Callable[..., SimpleNamespace]:
    r"""A tolerance classifier: outflow ``c > band``, inflow ``c < -band``, tangential ``|c| <= band`` (the retired spelling)."""

    def owner(quadrature, axis_index: int, *, sign: int = +1) -> SimpleNamespace:
        c = sign * np.asarray(quadrature.mean_axis_cosine(axis_index), dtype=np.float64)
        return _parts(c.size, np.flatnonzero(c > band), np.flatnonzero(c < -band), np.flatnonzero(np.abs(c) <= band))

    return owner


def unguarded_owner(quadrature, axis_index: int, *, sign: int = +1) -> SimpleNamespace:
    r"""An exact classifier that never refuses: no finite check, no rank check, no sign check."""
    c = sign * np.asarray(quadrature.mean_axis_cosine(axis_index), dtype=np.float64)
    c = c.reshape(-1)
    return _parts(c.size, np.flatnonzero(c > 0.0), np.flatnonzero(c < 0.0), np.flatnonzero(c == 0.0))


def exact_owner(quadrature, axis_index: int, *, sign: int = +1) -> SimpleNamespace:
    r"""The exact classifier with the refusals (finite, one-dimensional, sign), as light namespaces: the control of the instruments."""
    if sign not in (+1, -1):
        raise ValueError(f"sign must be +1 or -1, got {sign!r}")
    c = sign * np.asarray(quadrature.mean_axis_cosine(axis_index), dtype=np.float64)
    if c.ndim != 1:
        raise ValueError("cosines must be one-dimensional")
    if not np.all(np.isfinite(c)):
        raise ValueError("non-finite direction cosine")
    return _parts(c.size, np.flatnonzero(c > 0.0), np.flatnonzero(c < 0.0), np.flatnonzero(c == 0.0))


def no_cover_partition(*, outflow, inflow, tangential, n):
    r"""A ``HalfRange`` stand-in that checks the part extent and the overlap and never the cover (the dropped-check mutant)."""
    sets = [set(p.indices) for p in (outflow, inflow, tangential)]
    if any(p.n != n for p in (outflow, inflow, tangential)):
        raise ValueError("a part is not on the ordinate set")
    if sets[0] & sets[1] or sets[0] & sets[2] or sets[1] & sets[2]:
        raise ValueError("the three parts are not disjoint")
    return SimpleNamespace(outflow=outflow, inflow=inflow, tangential=tangential, n=n)


# ─────────────────────────────────────────────────────────────────────
# Instruments over an owner
# ─────────────────────────────────────────────────────────────────────


def _idx(part) -> tuple[int, ...]:
    return tuple(int(i) for i in part.indices)


def _reading(owner, cosines, sign: int = +1) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    h = owner(CosineStub(cosines), 0, sign=sign)
    return _idx(h.outflow), _idx(h.inflow), _idx(h.tangential)


def value_row_violations(owner) -> list[str]:
    r"""G02a value row: what ``owner`` gets wrong on the round-off cosine and on the subnormal stub (an empty list is correct)."""
    bad: list[str] = []
    expect = [
        ("[-0.5, 9.5e-16, 0.5]", [-0.5, 9.5e-16, 0.5], +1, ((1, 2), (0,), ())),
        ("[-0.5, 9.5e-16, 0.5] sign=-1", [-0.5, 9.5e-16, 0.5], -1, ((0,), (1, 2), ())),
        ("[0.0, -0.0, 5e-324, -5e-324]", [0.0, -0.0, 5e-324, -5e-324], +1, ((2,), (3,), (0, 1))),
    ]
    for label, cosines, sign, want in expect:
        got = _reading(owner, cosines, sign)
        if got != want:
            bad.append(f"{label}: (outflow, inflow, tangential) = {got}, wanted {want}")
    return bad


def refusal_violations(owner) -> list[str]:
    r"""G02d, the refusals of ``half_range``: each bad input that ``owner`` accepts is a violation (an empty list is correct)."""
    bad: list[str] = []
    for label, cosines, sign, fragment in [
        ("nan", [np.nan, 0.5], +1, "non-finite direction cosine"),
        ("+inf", [np.inf, 0.5], +1, "non-finite direction cosine"),
        ("-inf", [-np.inf, 0.5], +1, "non-finite direction cosine"),
        ("a two-dimensional cosine array", [[0.1, 0.2], [0.3, 0.4]], +1, "one-dimensional"),
        ("sign 2", [-0.5, 0.5], 2, None),
        ("sign 0", [-0.5, 0.5], 0, None),
    ]:
        try:
            owner(CosineStub(cosines), 0, sign=sign)
        except ValueError as exc:
            if fragment is not None and fragment not in str(exc):
                bad.append(f"{label}: refused with {str(exc)!r}, wanted a message containing {fragment!r}")
        except Exception as exc:  # a wrong exception class is a violation, not an error of the instrument
            bad.append(f"{label}: raised {type(exc).__name__}, wanted ValueError")
        else:
            bad.append(f"{label}: accepted without a refusal")
    return bad


def partition_violations(half, n: int, label: str) -> list[str]:
    r"""``half`` is a partition of ``range(n)`` into three disjoint parts, each on the ordinate set of extent ``n``."""
    bad: list[str] = []
    parts = {"outflow": half.outflow, "inflow": half.inflow, "tangential": half.tangential}
    if half.n != n:
        bad.append(f"{label}: n = {half.n}, wanted {n}")
    seen: list[int] = []
    for name, part in parts.items():
        if part.n != n:
            bad.append(f"{label}: {name}.n = {part.n}, wanted {n}")
        seen += list(_idx(part))
    if sorted(seen) != list(range(n)):
        bad.append(f"{label}: the three parts do not partition range({n}) (union with multiplicity has {len(seen)} entries)")
    return bad


def partition_refusal_violations(half_range_cls, coord_cls) -> list[str]:
    r"""G02d, the refusals of ``HalfRange.__post_init__``: a part off the ordinate set, an overlap, a non-cover (each with its fragment)."""
    bad: list[str] = []
    n = 3
    cases = [
        (
            "a part whose Coord.n differs from n",
            dict(outflow=coord_cls(n + 1, (n,)), inflow=coord_cls(n, (0, 1)), tangential=coord_cls(n, (2,)), n=n),
            "not on the ordinate set",
        ),
        (
            "overlapping parts",
            dict(outflow=coord_cls(n, (0, 1)), inflow=coord_cls(n, (1, 2)), tangential=coord_cls(n, ()), n=n),
            "not disjoint",
        ),
        (
            "parts that do not cover",
            dict(outflow=coord_cls(n, (0,)), inflow=coord_cls(n, (1,)), tangential=coord_cls(n, ()), n=n),
            "do not cover",
        ),
    ]
    for label, kwargs, fragment in cases:
        try:
            half_range_cls(**kwargs)
        except ValueError as exc:
            if fragment not in str(exc):
                bad.append(f"{label}: refused with {str(exc)!r}, wanted a message containing {fragment!r}")
        except Exception as exc:
            bad.append(f"{label}: raised {type(exc).__name__}, wanted ValueError")
        else:
            bad.append(f"{label}: accepted without a refusal")
    return bad


# ─────────────────────────────────────────────────────────────────────
# ``ast`` arms
# ─────────────────────────────────────────────────────────────────────

#: The 10 consumer functions of G02a's call arm: ``(path under orpheus/, dotted qualified name)``.
CONSUMERS: tuple[tuple[str, str], ...] = (
    ("mesh/axis.py", "face_outflow_ordinates"),
    ("numerics/quadrature/ordinate_pairing.py", "ordinate_pairing"),
    ("numerics/spaces/angular_trace_space.py", "AngularTraceSpace.from_quadrature_and_layout"),
    ("numerics/quadrature/directional.py", "Quadrature.octants"),
    ("sn/operators/loss_kernel_gauge.py", "_reflection_orbits"),
    ("sn/loss_representation/__init__.py", "_OneDimScanWalk._dag_legs"),
    ("sn/loss_representation/__init__.py", "_OneDimScanWalk._degenerate_positions"),
    ("sn/sweep/cache.py", "StreamingCoefficientCache.from_mesh_and_quad"),
    ("sn/problem.py", "SNProblem._representative_ordinate"),
    ("sn/problem.py", "SNProblem._iter_cylindrical_visits"),
)


def orpheus_root() -> Path:
    r"""The directory of the ``orpheus`` package that this process imports (read as a path, never imported from a layer above)."""
    import orpheus

    return Path(orpheus.__file__).resolve().parent


def find_function(tree: ast.AST, qualname: str) -> ast.FunctionDef | None:
    r"""The ``FunctionDef`` named by the dotted ``qualname`` (``Class.method`` or ``function``) in ``tree``, or ``None``."""
    scope: ast.AST = tree
    node: ast.AST | None = None
    for part in qualname.split("."):
        node = next(
            (
                n
                for n in getattr(scope, "body", [])
                if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name == part
            ),
            None,
        )
        if node is None:
            return None
        scope = node
    return node if isinstance(node, ast.FunctionDef) else None


def calls_name(node: ast.AST, name: str) -> bool:
    r"""``node`` contains a call whose callee is the name ``name`` or an attribute named ``name``."""
    for n in ast.walk(node):
        if isinstance(n, ast.Call):
            f = n.func
            if (isinstance(f, ast.Name) and f.id == name) or (isinstance(f, ast.Attribute) and f.attr == name):
                return True
    return False


def consumer_calls(root: Path, name: str) -> list[tuple[str, str, bool | None]]:
    r"""``(file, qualname, calls)`` for each of the 10 consumers: ``calls`` is ``None`` when the function is not found in its file."""
    out: list[tuple[str, str, bool | None]] = []
    for rel, qual in CONSUMERS:
        fn = find_function(ast.parse((root / rel).read_text()), qual)
        out.append((rel, qual, None if fn is None else calls_name(fn, name)))
    return out


def consumer_source_without_call(root: Path, rel: str, qual: str, name: str) -> tuple[str, str]:
    r"""``(function_source, source_with_the_call_renamed)``: the consumer's source text and a copy whose calls to ``name`` are gone."""
    src = (root / rel).read_text()
    fn = find_function(ast.parse(src), qual)
    if fn is None:
        raise LookupError(f"{rel}::{qual} is not found")
    seg = ast.get_source_segment(src, fn)
    if seg is None:
        raise LookupError(f"{rel}::{qual} has no source segment")
    import textwrap

    seg = textwrap.dedent(seg)
    return seg, re.sub(rf"\b{name}\b", "_the_owner_call_was_deleted", seg)


#: Names that hold a direction cosine in the scopes of G02g.
COSINE_NAMES = frozenset({"omega_dot_n", "mu", "mu_x", "mu_level", "mu_axis", "row", "nodes", "partner_mu", "abs_mu"})
_RANK_QUESTIONS = frozenset({"ndim", "shape", "size"})

#: The scopes of G02g: whole files, then the functions that share a file with sites that are not this unit's.
COSINE_SCOPE_FILES = (
    "mesh/axis.py",
    "numerics/quadrature/ordinate_pairing.py",
    "numerics/quadrature/directional.py",
    "numerics/spaces/angular_trace_space.py",
    "sn/sweep/cache.py",
)
COSINE_SCOPE_FUNCTIONS = {
    "sn/loss_representation/__init__.py": ("_dag_legs", "_degenerate_positions"),
    "sn/operators/loss_kernel_gauge.py": ("_reflection_orbits",),
}


def _reads_cosine(operand: ast.AST) -> bool:
    skipped = {
        id(y)
        for x in ast.walk(operand)
        if isinstance(x, ast.Attribute) and x.attr in _RANK_QUESTIONS
        for y in ast.walk(x.value)
    }
    return any(
        id(x) not in skipped
        and ((isinstance(x, ast.Name) and x.id in COSINE_NAMES) or (isinstance(x, ast.Attribute) and x.attr in COSINE_NAMES))
        for x in ast.walk(operand)
    )


def cosine_compares(tree: ast.AST, functions: tuple[str, ...] | None = None) -> list[int]:
    r"""Line numbers of the ``Compare`` nodes with an operand that reads a cosine array, in ``tree`` or in the named functions of it."""
    roots: list[ast.AST] = (
        [tree]
        if functions is None
        else [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in functions]
    )
    lines: set[int] = set()
    for root in roots:
        for n in ast.walk(root):
            if isinstance(n, ast.Compare) and any(_reads_cosine(s) for s in (n.left, *n.comparators)):
                lines.add(n.lineno)
    return sorted(lines)


def scope_cosine_compares(root: Path) -> list[tuple[str, int]]:
    r"""``(file, line)`` of every cosine-array ``Compare`` in the scopes of G02g."""
    hits: list[tuple[str, int]] = []
    for rel in COSINE_SCOPE_FILES:
        hits += [(rel, ln) for ln in cosine_compares(ast.parse((root / rel).read_text()))]
    for rel, funcs in COSINE_SCOPE_FUNCTIONS.items():
        hits += [(rel, ln) for ln in cosine_compares(ast.parse((root / rel).read_text()), funcs)]
    return hits


#: The six cutoff names of G02c (i).
CUTOFF_NAMES = ("TANGENTIAL_EPS", "_OUTFLOW_COSINE_TOL", "_TANGENTIAL_MU", "_MU_DIRECTION_EPS", "_OCTANT_SIGN_EPS", "_DEGENERATE_ABS_ETA_THRESHOLD")

#: The 9 edited production files that held a cutoff (G02c (ii)).
CUTOFF_FILES = (
    "mesh/axis.py",
    "numerics/tangential.py",
    "numerics/quadrature/ordinate_pairing.py",
    "numerics/quadrature/directional.py",
    "numerics/spaces/angular_trace_space.py",
    "sn/loss_representation/__init__.py",
    "sn/operators/loss_kernel_gauge.py",
    "sn/problem.py",
    "sn/sweep/cache.py",
)


def cutoff_name_hits(source: str, names: tuple[str, ...] | frozenset[str] = CUTOFF_NAMES) -> list[tuple[str, int]]:
    r"""``(name, line)`` of every definition, import and read of one of ``names`` in a parse of ``source``."""
    hits: set[tuple[str, int]] = set()
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Name) and n.id in names:
            hits.add((n.id, n.lineno))
        elif isinstance(n, ast.Attribute) and n.attr in names:
            hits.add((n.attr, n.lineno))
        elif isinstance(n, ast.alias) and (n.name in names or n.asname in names):
            hits.add((n.name, getattr(n, "lineno", 0)))
        elif isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name in names:
            hits.add((n.name, n.lineno))
    return sorted(hits)


def _tiny(v: object) -> bool:
    return isinstance(v, float) and 0 < abs(v) <= 1e-6


def _has_finfo(n: ast.AST) -> bool:
    return any(
        (isinstance(x, ast.Attribute) and x.attr == "finfo") or (isinstance(x, ast.Name) and x.id == "finfo") for x in ast.walk(n)
    )


def _literal(n: ast.AST) -> object:
    try:
        return ast.literal_eval(n)
    except Exception:
        return None


def tolerance_compares(source: str) -> list[tuple[int, str]]:
    r"""``(line, source segment)`` of each ``Compare`` that carries a tolerance operand: a float literal ``0 < v <= 1e-6``, a
    name bound to one or to a ``finfo`` expression at module level, in a class body or in a function (read bare or as an attribute),
    an imported name ending ``EPS`` or ``TOL``, or ``.eps``."""
    tree = ast.parse(source)
    tolerance_names: set[str] = set()
    for st in ast.walk(tree):  # every binding: module level, a class body (a ``ClassVar``), a function local
        target: str | None = None
        value: ast.AST | None = None
        if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name):
            target, value = st.targets[0].id, st.value
        elif isinstance(st, ast.AnnAssign) and isinstance(st.target, ast.Name) and st.value is not None:
            target, value = st.target.id, st.value
        if target is not None and value is not None and (_tiny(_literal(value)) or _has_finfo(value)):
            tolerance_names.add(target)
        if isinstance(st, ast.ImportFrom):
            for a in st.names:
                nm = a.asname or a.name
                if nm.endswith("EPS") or nm.endswith("TOL"):
                    tolerance_names.add(nm)
    out: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        hit = False
        for operand in (node.left, *node.comparators):
            for x in ast.walk(operand):
                if isinstance(x, ast.Constant) and _tiny(x.value):
                    hit = True
                elif isinstance(x, ast.Name) and x.id in tolerance_names:
                    hit = True
                elif isinstance(x, ast.Attribute) and (x.attr == "eps" or x.attr in tolerance_names):
                    hit = True
        if hit:
            out.append((node.lineno, ast.get_source_segment(source, node) or ""))
    return sorted(out)


#: The one comparison with a tolerance operand that survives in the 9 files: ``(file, source segment)``.
TOLERANCE_ALLOW_LIST = frozenset({("sn/operators/loss_kernel_gauge.py", "reading <= _GAUGE_AUDIBLE_FLOOR")})
