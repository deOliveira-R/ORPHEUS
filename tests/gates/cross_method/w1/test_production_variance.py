r"""G04h and G04v: production mints each factor with its variance stated, and the function role is derived (N04, spec 7.5; spec 3.4a).

G04h cell check: the flux space of the SN slab, the SN sphere and the diffusion slab have angular ``Variance(MEASURE, DENSITY)``, energy ``FormlessVariance.MEASURE`` and spatial
``Variance(MEASURE, DENSITY)`` (the diffusion space has no angular slot). ``ast`` leg: every CALL in ``orpheus/`` of a ``TensorSpace``/``FunctionSpace`` ``of_measure``, ``of_basis`` or ``declared``, of
``Factor(`` and of ``Factor.whole_on(``, outside ``numerics/factor.py`` and outside the bodies of the constructors ``of_axes``, ``of_measure``, ``of_basis``, ``declared`` in ``numerics/space.py`` (which mint
``Factor(axis, variance, ...)`` with their parameter), states its variance as ``Variance(Kind.<member>, VarianceChart.<member>)``, ``FormlessVariance.<member>``, a module-level constant bound to
one, ``<factor>.variance``, or ``<factor>.variance.dual()`` / ``.riesz()`` (the re-mints of ``dual``, ``riesz`` and ``completed``); ``Factor.whole_on`` (the one branch on a form outside the factor's own invariants, ``ELEGANCE-DEBT[guard] #757``; decision ``[ON04-4]``, amended) takes a literal ``Kind.<member>`` and a literal
``VarianceChart.<member>``; a call with no variance argument is a red. G04v: ``bulk_space.dual().riesz(on={"spatial"})`` has energy ``FormlessVariance.FUNCTION`` and spatial
``(FUNCTION, DENSITY)`` and equals the importance row built from the bulk space's own axes; ``riesz(on={"energy"})`` is refused.

First red at the base tip ``cf924aaa``: the cell check is import-/attribute-only (``FunctionSpace`` has no ``factor``), G04v likewise (no ``riesz``); the ``ast`` leg is behavioural: it reads 3 builder calls
in ``numerics/frame_faces.py`` that state no variance (``of_basis`` at line 442, ``of_measure`` at 488 and a third; the opening's "0 builder calls" counted calls with a non-literal variance argument). Mutations: the angular mint with ``Variance(FUNCTION, DENSITY)`` reddens the cell check;
the gate's own planted ``Factor(`` call passing a variable, a ``Factor.whole_on(`` call passing a variable kind, a ``Factor.whole_on(`` call with no kind, a ``Factor(`` call passing ``Variance(f.kind, f.chart)`` and ``Variance(kind_variable, chart)`` at a builder each redden the ``ast`` predicate
(positive controls inside the file); ``dual`` as the identity and ``riesz`` as the identity each redden G04v. Exact on any host (the cells are variances and identities, the ``ast`` leg a parse).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.gates.cross_method.w1 import _fixtures as fx
from tests.gates.cross_method.w1 import _snapshot as snap

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
BUILDERS = {"of_measure": 2, "of_basis": 1, "declared": 2}  # name -> position of the variance argument
CONSTRUCTOR_BODIES = {"of_axes", "of_measure", "of_basis", "declared"}


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _getter(rows: tuple, name: str):
    return [g for n, g, _ in rows if n == name][0]


def _cells(space) -> dict[str, object]:
    return {f.axis.label: f.variance for f in space.factors}


def test_the_flux_spaces_carry_the_declared_variance() -> None:
    from orpheus.numerics.factor import FormlessVariance, Kind, Variance, VarianceChart

    measure_density = Variance(Kind.MEASURE, VarianceChart.DENSITY)
    for kind in ("slab", "sphere"):
        problem, _ = fx.sn_problem(kind, 2)
        cells = _cells(_getter(snap.SN_ROWS, "angular_bulk_space")(problem))
        _require(cells == {"angular": measure_density, "energy": FormlessVariance.MEASURE, "spatial": measure_density}, f"SN {kind} flux space variances: {cells}")
    mesh, _solver, _ = fx.diffusion_problem(2)
    cells = _cells(_getter(snap.DIFFUSION_ROWS, "bulk_space")(mesh))
    _require(cells == {"energy": FormlessVariance.MEASURE, "spatial": measure_density}, f"diffusion slab bulk space variances: {cells}")


def test_sd35b_function_role_is_derived() -> None:
    from orpheus.numerics.factor import Factor, FormlessVariance, Kind, Subquotient, Variance, VarianceChart
    from orpheus.numerics.space import TensorSpace

    mesh, _solver, _ = fx.diffusion_problem(2)
    bulk = _getter(snap.DIFFUSION_ROWS, "bulk_space")(mesh)
    role = bulk.dual().riesz(on={"spatial"})
    cells = _cells(role)
    _require(cells["energy"] == FormlessVariance.FUNCTION, f"the energy slot of the function role is {cells['energy']}")
    _require(cells["spatial"] == Variance(Kind.FUNCTION, VarianceChart.DENSITY), f"the spatial slot of the function role is {cells['spatial']}")
    e, s = bulk.factor("energy"), bulk.factor("spatial")
    own = TensorSpace.of_axes(
        Factor(e.axis, FormlessVariance.FUNCTION, Subquotient.whole(e.axis.shape[0])),
        Factor(s.axis, Variance(Kind.FUNCTION, VarianceChart.DENSITY), Subquotient.whole(s.axis.shape[0])),
    )
    _require(role == own, "the derived function role differs from the importance row built from the bulk space's own axes")
    with pytest.raises(Exception, match="riesz is not defined: "):
        bulk.dual().riesz(on={"energy"})


# ---- the ``ast`` predicate ------------------------------------------------------------------------------------------------------------------


def _is_literal_variance(node: ast.expr, constants: set[str]) -> bool:
    if isinstance(node, ast.Name):
        return node.id in constants
    if isinstance(node, ast.Attribute):
        if node.attr == "variance":
            return True
        return isinstance(node.value, (ast.Name, ast.Attribute)) and _tail(node.value) == "FormlessVariance"
    if isinstance(node, ast.Call):
        f = node.func
        if _tail(f) == "Variance":
            return len(node.args) == 2 and _is_member(node.args[0], "Kind") and _is_member(node.args[1], "VarianceChart")
        if isinstance(f, ast.Attribute) and f.attr in ("dual", "riesz"):
            return isinstance(f.value, ast.Attribute) and f.value.attr == "variance"
    return False


def _tail(node: ast.expr) -> str:
    return node.id if isinstance(node, ast.Name) else node.attr if isinstance(node, ast.Attribute) else ""


def _is_member(node: ast.expr, enum: str) -> bool:
    return isinstance(node, ast.Attribute) and _tail(node.value) == enum


def _arg(call: ast.Call, position: int, keyword: str) -> ast.expr | None:
    if len(call.args) > position:
        return call.args[position]
    return next((k.value for k in call.keywords if k.arg == keyword), None)


def violations(source: str, *, skip_constructor_bodies: bool = False) -> list[str]:
    tree = ast.parse(source)
    constants = {t.id for n in tree.body if isinstance(n, ast.Assign) and _is_literal_variance(n.value, set()) for t in n.targets if isinstance(t, ast.Name)}
    skipped: set[int] = set()
    if skip_constructor_bodies:
        for n in ast.walk(tree):
            if isinstance(n, ast.FunctionDef) and n.name in CONSTRUCTOR_BODIES:
                skipped |= {id(c) for c in ast.walk(n)}
    out: list[str] = []
    for c in (n for n in ast.walk(tree) if isinstance(n, ast.Call) and id(n) not in skipped):
        name = _tail(c.func)
        owner = _tail(c.func.value) if isinstance(c.func, ast.Attribute) else ""
        if name in BUILDERS and owner in ("TensorSpace", "FunctionSpace"):
            a = _arg(c, BUILDERS[name], "variance")
        elif name == "Factor":
            a = _arg(c, 1, "variance")
        elif name == "whole_on" and owner == "Factor":
            kind, chart = _arg(c, 1, "kind"), _arg(c, 2, "chart")
            if kind is None or chart is None or not (_is_member(kind, "Kind") and _is_member(chart, "VarianceChart")):
                out.append(f"line {c.lineno}: whole_on with a non-literal or missing kind or chart")
            continue
        else:
            continue
        if a is None:
            out.append(f"line {c.lineno}: {name}( states no variance")
        elif not _is_literal_variance(a, constants):
            out.append(f"line {c.lineno}: {name}( states its variance as {ast.unparse(a)}")
    return out


def test_the_ast_predicate_finds_each_planted_violation() -> None:
    ok = "Factor(axis, Variance(Kind.MEASURE, VarianceChart.DENSITY), sk)\nFactor(a, f.variance, sk)\nFactor(a, FormlessVariance.MEASURE, sk)\nFactor.whole_on(axis, Kind.MEASURE, VarianceChart.DENSITY)\nFactor(a, g.variance.dual(), sk)"
    _require(violations(ok) == [], f"a clean source reads {violations(ok)}")
    _require(len(violations("Factor(axis, v, sk)")) == 1, "a Factor( call passing a variable is not found")
    _require(len(violations("Factor.whole_on(axis, k, VarianceChart.DENSITY)")) == 1, "a Factor.whole_on( call passing a variable kind is not found")
    _require(len(violations("Factor.whole_on(axis)")) == 1, "a Factor.whole_on( call with no kind is not found")
    _require(len(violations("Factor(axis, Variance(f.kind, f.chart), sk)")) == 1, "a mint passing Variance(f.kind, f.chart) is not found")
    _require(len(violations("Factor(axis)")) == 1, "a Factor( mint with no variance is not found")
    _require(len(violations("Factor(axis, Variance(kind_var, VarianceChart.DENSITY), sk)")) == 1, "Variance(kind_variable, chart) at a builder is not found")
    _require(len(violations("TensorSpace.of_measure(m, 'x')")) == 1, "an of_measure( call with no variance is not found")


def test_every_production_mint_states_its_variance_literally() -> None:
    found: list[str] = []
    population = 0
    for path in sorted((ROOT / "orpheus").rglob("*.py")):
        rel = path.relative_to(ROOT).as_posix()
        if rel == "orpheus/numerics/factor.py":
            continue
        population += 1
        for v in violations(path.read_text(), skip_constructor_bodies=(rel == "orpheus/numerics/space.py")):
            found.append(f"{rel}: {v}")
    _require(population > 200, f"the census read {population} files")
    _require(not found, f"{len(found)} production mints do not state a literal variance: {found[:6]}")
