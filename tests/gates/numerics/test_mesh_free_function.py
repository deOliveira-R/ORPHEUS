r"""The laws of the mesh-free functions (#405 P1 step 6, S6.8, S6.11-S6.17, S6.21).

:class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant` is a
function on the angle-integrated space, one finite real per (region, group);
:class:`~orpheus.numerics.mesh_free_function.Symbolic` is a function on
phase space, one SymPy expression per group in the class's own coordinates
``r, mu, phi``, stored as ``srepr`` text with the version of the SymPy that
wrote it. Their content-identity rows are
``test_content_identity_mesh_free.py`` (S6.9); the Branch-1 lifts are
``tests/gates/derivations/test_angular_measure_symbolic.py`` (S6.10).

The specification is ``.claude/plans/reference_p1_spec.md`` §1.6.
"""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]

r, mu, phi = Symbolic.r, Symbolic.mu, Symbolic.phi


# ── S6.8: RegionwiseConstant's construction ─────────────────────────────


def test_s6_8_the_table_is_a_read_only_copy() -> None:
    caller = np.array([[1.0, 2.0], [3.0, 4.0]])
    table = RegionwiseConstant(caller)
    digest = table.content_digest
    caller[0, 0] = 99.0
    npt.assert_array_equal(table.values, [[1.0, 2.0], [3.0, 4.0]])
    if table.values.flags.writeable:
        pytest.fail("the stored table is writeable")
    if RegionwiseConstant(np.array([[1.0, 2.0], [3.0, 4.0]])).content_digest != digest:
        pytest.fail("the digest moved with the caller's array")


@pytest.mark.parametrize(
    "values,fragment",
    [
        (np.array(1.0), "of rank 2, got rank 0"),  # qa 2026-10-02: the #559 parse first raised numpy's unkeyed 0-d error
        (np.array([1.0, 2.0]), "of rank 2, got rank 1"),
        (np.ones((2, 2, 2)), "of rank 2, got rank 3"),
        (np.ones((0, 2)), "has no regions"),
        (np.ones((2, 0)), "has no groups"),
        (np.array([[1.0, np.nan]]), r"entry \(0, 1\) is NaN"),
        (np.array([[1.0], [-np.inf]]), r"entry \(1, 0\) is infinite"),
    ],
    ids=["rank 0", "rank 1", "rank 3", "no regions", "no groups", "NaN", "inf"],
)
def test_s6_8_the_refusals(values: np.ndarray, fragment: str) -> None:
    with pytest.raises(ValueError, match=fragment):
        RegionwiseConstant(values)


def test_s6_8_a_non_real_table_is_refused() -> None:
    for values in (np.array([[1 + 2j]]), np.array([[True]]), np.array([["a"]])):
        with pytest.raises(TypeError, match="must be real numbers"):
            RegionwiseConstant(values)


# ── S6.11 / S6.12: no role, no density; the group axis ───────────────────


def test_s6_11_neither_type_carries_a_role_or_a_density() -> None:
    """The role (source or detector) is the specification's field, and it
    picks the arrow (E or R†); a type field naming it would let a value carry
    a lift the field contradicts."""
    if [f.name for f in dataclasses.fields(RegionwiseConstant)] != ["values"]:
        pytest.fail(f"RegionwiseConstant fields: {[f.name for f in dataclasses.fields(RegionwiseConstant)]}")
    if [f.name for f in dataclasses.fields(Symbolic)] != ["srepr", "sympy_version"]:
        pytest.fail(f"Symbolic fields: {[f.name for f in dataclasses.fields(Symbolic)]}")


def test_s6_12_the_group_count() -> None:
    table = RegionwiseConstant(np.ones((3, 4)))
    if (table.n_regions, table.n_groups) != (3, 4):
        pytest.fail(f"(regions, groups) = {(table.n_regions, table.n_groups)}, not (3, 4)")
    if Symbolic.of(r, mu, 1).n_groups != 3:
        pytest.fail("a 3-expression Symbolic does not report 3 groups")


# ── S6.13: the owned coordinates ─────────────────────────────────────────


def test_s6_13_the_owned_symbols() -> None:
    for name, symbol in (("r", r), ("mu", mu), ("phi", phi)):
        if symbol != sp.Symbol(name, real=True):
            pytest.fail(f"Symbolic.{name} is {symbol!r} with {symbol.assumptions0}, not Symbol({name!r}, real=True)")


@pytest.mark.parametrize(
    "stray,fragment",
    [
        (sp.Symbol("a0"), "free symbol 'a0' with assumptions"),
        (sp.Symbol("x"), "free symbol 'x' with assumptions"),
        (sp.Symbol("mu"), "same name, different assumptions: use Symbolic.mu"),
        # The tree's own MMS builder spells its coordinates this way
        # (derivations/continuous/mms/sn.py), so a check by name admits it.
        (sp.Symbol("mu", positive=True, real=True), "same name, different assumptions: use Symbolic.mu"),
        (sp.Symbol("r", positive=True), "same name, different assumptions: use Symbolic.r"),
    ],
    ids=["a0", "x", "mu bare", "mu positive", "r positive"],
)
def test_s6_13_a_stray_symbol_is_refused_by_name_and_assumptions(stray: sp.Symbol, fragment: str) -> None:
    with pytest.raises(ValueError, match=fragment):
        Symbolic.of(1 + stray * r)


# ── S6.14: srepr round trip and process stability ────────────────────────


def _population() -> dict[str, tuple[sp.Expr, ...]]:
    many = sum(sp.Rational(k + 1, 7) * r**k * mu ** (k % 3) for k in range(12))
    return {
        "polynomial": (1 + 2 * r - 3 * r**2 * mu + mu**3,),
        "trigonometric": (sp.sin(sp.pi * r) * sp.exp(-mu) + sp.cos(phi) * sp.sqrt(1 - mu**2) * r,),
        "piecewise": (sp.Piecewise((sp.Float("1.234567890123456789012345678901", 30), r < sp.Rational(1, 2)),
                                   (sp.Rational(3, 7) * mu, True)),),
        "twelve terms": (many, many / 2),
    }


@pytest.mark.parametrize("name", list(_population()))
def test_s6_14_the_text_round_trips(name: str) -> None:
    value = Symbolic.of(*_population()[name])
    back = Symbolic.from_srepr(value.srepr, value.sympy_version)
    if back != value or back.content_digest != value.content_digest:
        pytest.fail(f"{name}: from_srepr(s.srepr) is not s")
    for original, parsed in zip(_population()[name], back.expressions):
        if sp.simplify(original - parsed) != 0:
            pytest.fail(f"{name}: the stored text parses to {parsed}, not {original}")
        if parsed.free_symbols - {r, mu, phi}:
            pytest.fail(f"{name}: the assumptions did not survive: {parsed.free_symbols}")


def test_s6_14_the_digest_is_the_same_in_every_process() -> None:
    """The digests printed under PYTHONHASHSEED 0 and 1 agree. A digest over
    ``hash(expr)`` would differ (control: ``hash`` of a Symbol differs)."""
    script = textwrap.dedent('''
        import sympy as sp
        from orpheus.numerics.mesh_free_function import Symbolic
        r, mu, phi = Symbolic.r, Symbolic.mu, Symbolic.phi
        s = Symbolic.of(sp.sin(sp.pi * r) * sp.exp(-mu) + sp.cos(phi) * r, 1 + mu)
        print(s.content_digest.hex(), hash(sp.Symbol("q")))
    ''')
    outputs = []
    for seed in ("0", "1"):
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(_ROOT)}
        result = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, env=env, cwd=_ROOT, check=True)
        outputs.append(result.stdout.split())
    if outputs[0][1] == outputs[1][1]:
        pytest.fail("control: hash() of a Symbol did not change between seeds, so this row cannot see a hash() digest")
    if outputs[0][0] != outputs[1][0]:
        pytest.fail(f"the digest moved between seeds: {outputs[0][0]} vs {outputs[1][0]}")


# ── S6.15: the anisotropy predicate ─────────────────────────────────────


@pytest.mark.parametrize(
    "expression,isotropic",
    [
        (1 + mu, False),
        (sp.cos(phi), False),
        (sp.Integer(2), True),
        (r**2 + 1, True),
        (mu - mu, True),
        # The traps a free-symbols test calls anisotropic.
        (sp.sin(phi) ** 2 + sp.cos(phi) ** 2, True),
        (mu**2 + (1 - mu**2) * sp.cos(phi) ** 2 + (1 - mu**2) * sp.sin(phi) ** 2, True),
        # Steps in the direction: zero derivative wherever defined, so a
        # derivative test calls them isotropic (qa, 2026-10-02).
        (sp.Piecewise((1, mu > 0), (0, True)), False),
        (sp.Piecewise((1, phi < sp.pi), (0, True)), False),
        # A step in position only is isotropic.
        (sp.Piecewise((1, r < 1), (2, True)), True),
    ],
    ids=["1+mu", "cos phi", "2", "r^2+1", "mu-mu", "sin^2+cos^2", "Omega.Omega", "step in mu", "step in phi", "step in r"],
)
def test_s6_15_the_anisotropy_predicate(expression: sp.Expr, isotropic: bool) -> None:
    if Symbolic.of(expression).is_isotropic is not isotropic:
        pytest.fail(f"{expression}: is_isotropic should be {isotropic}")


def test_s6_15_one_anisotropic_group_makes_the_function_anisotropic() -> None:
    if Symbolic.of(r, 1 + mu).is_isotropic:
        pytest.fail("a function with one anisotropic group read isotropic")


# ── S6.16: storage is text; identity is by spelling; non-values refused ─


def test_s6_16_storage_is_text() -> None:
    value = Symbolic.of(1 + mu)
    if not all(isinstance(part, str) for part in value.srepr) or not isinstance(value.sympy_version, str):
        pytest.fail(f"the content parts are not text: {value.content_parts()}")


def test_s6_16_identity_is_by_spelling() -> None:
    """(r+1)² and its expansion are one function and two values: a cache
    miss, never a wrong hit."""
    folded, expanded = (r + 1) ** 2, r**2 + 2 * r + 1
    if sp.simplify(folded - expanded) != 0:
        pytest.fail("control: the two spellings are not one function")
    if Symbolic.of(folded) == Symbolic.of(expanded):
        pytest.fail("two spellings of one function compared equal; identity is by spelling")


@pytest.mark.parametrize("non_value", [sp.nan, sp.zoo, sp.oo, -sp.oo], ids=["nan", "zoo", "oo", "-oo"])
def test_s6_16_a_non_value_is_refused(non_value: sp.Expr) -> None:
    with pytest.raises(ValueError, match="is not a function value"):
        Symbolic.of(r + non_value * mu)


@pytest.mark.parametrize(
    "candidate,fragment",
    [
        (r > 1, "is a StrictGreaterThan, not a scalar expression"),
        (sp.I * mu, "contains the imaginary unit"),
        (sp.Function("f")(r), r"applies the undefined function\(s\) \['f'\]"),
        (sp.Piecewise((1, r < 1)), "Piecewise with no otherwise branch"),
    ],
    ids=["relation", "imaginary", "undefined function", "piecewise without otherwise"],
)
def test_s6_16_only_real_scalar_functions_are_admitted(candidate, fragment: str) -> None:
    """The cases qa found admitted (2026-10-02): none is a real function value
    everywhere on phase space."""
    with pytest.raises(ValueError, match=fragment):
        Symbolic.of(candidate)


@pytest.mark.parametrize(
    "text,fragment",
    [
        ("__import__('os').getcwd()", "holds a Attribute"),
        ("getattr(Integer(1), 'func')", "names 'getattr'"),
        ("Integer(1)[0]", "holds a Subscript"),
        ("(lambda: Integer(1))()", "holds a Lambda"),
        ("Matrix([[Integer(1)]])", "names 'Matrix'"),
        ("Integer(1", "is not srepr text"),
    ],
    ids=["import", "builtin", "subscript", "lambda", "a matrix", "syntax"],
)
def test_s6_16_stored_text_is_parsed_through_a_whitelist(text: str, fragment: str) -> None:
    """``from_srepr`` evaluated stored text with ``sympify`` (``eval``), so a
    stored string could run code (qa, 2026-10-02: one wrote a file). Now only
    calls to SymPy classes and constants, literals and keywords are parsed,
    with no builtins."""
    with pytest.raises(ValueError, match=fragment):
        Symbolic.from_srepr((text,))


# ── S6.17: the SymPy version is content ─────────────────────────────────


def test_s6_17_the_version_is_this_sympys() -> None:
    if Symbolic.of(1 + mu).sympy_version != sp.__version__:
        pytest.fail(f"the version part is not sympy.__version__ = {sp.__version__}")


def test_s6_17_a_new_sympy_writes_a_new_value(monkeypatch) -> None:
    before = Symbolic.of(1 + mu)
    monkeypatch.setattr(sp, "__version__", "99.0.0")
    after = Symbolic.of(1 + mu)
    if after == before or after.content_digest == before.content_digest:
        pytest.fail("a value written by another SymPy has the same content")


# Written by SymPy 1.14.0 on 2026-10-02 (the producer pin: printed by
# ``Symbolic.of`` on the population of S6.14, never typed).
_SREPR_RECORD = {
    "polynomial": (
        "Add(Pow(Symbol('mu', real=True), Integer(3)), Mul(Integer(-1), Integer(3), Symbol('mu', real=True), "
        "Pow(Symbol('r', real=True), Integer(2))), Mul(Integer(2), Symbol('r', real=True)), Integer(1))",
    ),
    "trigonometric": (
        "Add(Mul(Symbol('r', real=True), Pow(Add(Integer(1), Mul(Integer(-1), Pow(Symbol('mu', real=True), "
        "Integer(2)))), Rational(1, 2)), cos(Symbol('phi', real=True))), Mul(exp(Mul(Integer(-1), "
        "Symbol('mu', real=True))), sin(Mul(pi, Symbol('r', real=True)))))",
    ),
    "piecewise": (
        "Piecewise(ExprCondPair(Float('1.23456789012345678901234567890097', precision=103), "
        "StrictLessThan(Symbol('r', real=True), Rational(1, 2))), ExprCondPair(Mul(Rational(3, 7), "
        "Symbol('mu', real=True)), true))",
    ),
}


def test_s6_17_record_the_srepr_text() -> None:
    """RECORD. Red means SymPy changed its ``srepr`` format: every ``Symbolic``
    cache key is invalidated (the version part already forces the miss);
    re-pin with the version."""
    for name, expected in _SREPR_RECORD.items():
        written = Symbolic.of(*_population()[name]).srepr
        if written != expected:
            pytest.fail(
                f"SymPy {sp.__version__} changed its srepr format for {name!r}: every Symbolic cache key "
                f"is invalidated (the version part already forces the miss); re-pin with the version.\n"
                f"written:  {written}\nrecorded: {expected}"
            )


# ── S6.21: the layer ─────────────────────────────────────────────────────


def test_s6_21_importing_the_module_loads_no_sympy() -> None:
    """In a fresh interpreter: importing the module and building a table load
    no SymPy; building a Symbolic does (the control that the probe can see it)."""
    script = textwrap.dedent('''
        import sys
        import numpy as np
        from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
        after_import = "sympy" in sys.modules
        RegionwiseConstant(np.ones((2, 2))).content_digest
        after_table = "sympy" in sys.modules
        Symbolic.of(1)
        print(after_import, after_table, "sympy" in sys.modules)
    ''')
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, env=env, cwd=_ROOT, check=True).stdout.split()
    if out != ["False", "False", "True"]:
        pytest.fail(f"(sympy after import, after a table, after a Symbolic) = {out}")


def test_s6_21_the_module_imports_no_geometry() -> None:
    """The module is chart-relative: it names the chart's home in prose and
    imports nothing from ``orpheus.geometry`` (its own import statements, by
    AST; the ``orpheus.numerics`` package already loads
    ``geometry.transformation`` through other modules)."""
    import ast

    import orpheus.numerics.mesh_free_function as module

    tree = ast.parse(Path(module.__file__).read_text())
    imported = [
        node.module if isinstance(node, ast.ImportFrom) else alias.name
        for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    ]
    if "orpheus.numerics.content" not in imported:
        pytest.fail(f"activation: the AST pass did not see the module's known import; it saw {imported}")
    if any(str(name).startswith("orpheus.geometry") for name in imported):
        pytest.fail(f"the module imports geometry: {imported}")
