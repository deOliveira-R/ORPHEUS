r"""**G5.1** — no operator of the transport family dispatches on its
operand's CARRIER (CS4c step 5: *each binding acts through the body its
ends select*).

The AST census the step's done-when names, with its predicate STATED
(F-2 of the verification plan: a lexical predicate reads 0 while a
family parse lives one frame out, so the census walks the whole
package, not the verbs alone):

* **no ``singledispatchmethod``** anywhere under
  ``orpheus/transport/operators/`` (`[M]` 3 dispatch tables / 13 arms at
  ``f90f7914``);
* **no ``isinstance(·, <carrier>)``** anywhere under the package, where a
  carrier is a transport field leaf, a field family base, ``Field``,
  ``FullField`` or ``np.ndarray`` — EXCEPT the three declared carve-outs,
  the family's ONLY carrier parses:
  :func:`~orpheus.transport.operators.lift.admit_composite` (the
  composite admission), :func:`~orpheus.transport.operators.lift.admit_array`
  (the plain admission), and the constructors' SPACE parse
  (``isinstance(space, FullFieldSpace)`` — a space, not a carrier; not
  counted). `[M]` 12 carrier ``isinstance`` inside the verbs + 3 in
  helpers they call at ``f90f7914``;
* **the composite is assembled in ONE place** — every ``FullField(...)``
  construction under the package is inside
  :func:`~orpheus.transport.operators.lift.lift_bulk_action` (`[M]` 9
  hand spellings of the zero-trace emission before the carve), and no
  boundary leaf's ``.zeros(`` is spelled by name anywhere under it.

⚠ Two filters, validated: the census runs on a SYNTHETIC module carrying
each shape it must catch (a dispatcher, a lexical carrier arm, a helper
carrier arm, a hand-spelled zero trace) before its production zero is
believed (vv #17 — the harness lies before the code does).
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

_PACKAGE = Path("orpheus/transport/operators")

#: The carrier vocabulary an ``isinstance`` may not parse (the 20 leaves,
#: the family bases, the ABC, the composite, the bare array).
_CARRIERS = {
    "AngularFlux", "ScalarFlux", "HarmonicMomentFlux", "AngularBoundaryFlux",
    "ScalarBoundaryFlux", "RadialCharacteristicInteriorFlux",
    "RadialCharacteristicBoundaryFlux",
    "AngularSourceSink", "ScalarSourceSink", "HarmonicMomentSourceSink",
    "AngularBoundarySourceSink", "ScalarBoundarySourceSink",
    "RadialCharacteristicInteriorSourceSink",
    "RadialCharacteristicBoundarySourceSink",
    "AngularResidual", "ScalarResidual", "AngularBoundaryResidual",
    "RadialCharacteristicInteriorResidual",
    "RadialCharacteristicBoundaryResidual",
    # (CrossSectionField is a COEFFICIENT leaf, never an operand — the
    # multiplier's from_mesh parses it as an INPUT, which is not dispatch.)
    "BulkField", "AngularField", "ScalarField", "MomentField", "FaceField",
    "BoundaryField", "Field", "FullField", "TimedFullField", "ndarray",
}

#: The declared carve-outs: (module, enclosing class or None, function) →
#: the family's only carrier parses — the two admission verbs, and the
#: moment factor's typed arms (R-5 of the step-5 design round: the
#: role-changing edge ``HarmonicMomentFlux → HarmonicMomentSourceSink`` is
#: kept legible in Λ's signature; its consumer is the moment end's typed
#: route, and it is a parse of the MOMENT carrier by the operator that
#: owns that carrier, not a dispatch on which body to run).
_CARVE_OUTS = {
    ("lift.py", None, "admit_composite"),
    ("lift.py", None, "admit_array"),
    ("transfer.py", "LegendreMomentTransfer", "apply"),
    ("transfer.py", "LegendreMomentTransfer", "apply_transpose"),
}
#: The one place the composite is assembled.
_LIFT_SITE = ("lift.py", None, "lift_bulk_action")


def _name_of(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return ""


def _isinstance_targets(call: ast.Call) -> set[str]:
    """The class names an ``isinstance(x, T)`` / ``isinstance(x, (T, U))`` tests."""
    if _name_of(call.func) != "isinstance" or len(call.args) != 2:
        return set()
    target = call.args[1]
    parts = target.elts if isinstance(target, ast.Tuple) else [target]
    return {_name_of(t) for t in parts}


def _census(source: str, module: str) -> dict[str, list[tuple[str, int]]]:
    """Walk one module; report every offence by kind with (function, line)."""
    tree = ast.parse(source)
    found: dict[str, list[tuple[str, int]]] = {
        "dispatcher": [], "carrier_isinstance": [], "full_field_ctor": [],
        "hand_zero_trace": [],
    }
    for func in [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
        for deco in func.decorator_list:
            if _name_of(deco) == "singledispatchmethod" or (
                isinstance(deco, ast.Attribute) and deco.attr == "register"
            ):
                found["dispatcher"].append((func.name, func.lineno))
    # module-scope statements too (a dispatcher can be assigned, not decorated)
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id == "singledispatchmethod":
            found["dispatcher"].append(("<module>", node.lineno))
    # (class or None, function) enclosing each line — the carve-out key.
    enclosing: dict[int, tuple[str | None, str]] = {}

    def _index(body: list[ast.stmt], cls: str | None) -> None:
        for node in body:
            if isinstance(node, ast.ClassDef):
                _index(node.body, node.name)
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for sub in ast.walk(node):
                    if hasattr(sub, "lineno"):
                        enclosing.setdefault(sub.lineno, (cls, node.name))

    _index(tree.body, None)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        cls, func_name = enclosing.get(node.lineno, (None, "<module>"))
        where = (module, cls, func_name)
        label = f"{cls}.{func_name}" if cls else func_name
        targets = _isinstance_targets(node) & _CARRIERS
        if targets and where not in _CARVE_OUTS:
            found["carrier_isinstance"].append((label, node.lineno))
        if _name_of(node.func) == "FullField" and where != _LIFT_SITE:
            found["full_field_ctor"].append((label, node.lineno))
        if (
            isinstance(node.func, ast.Attribute) and node.func.attr == "zeros"
            and (_name_of(node.func.value).endswith("BoundarySourceSink")
                 or _name_of(node.func.value).endswith("BoundaryFlux"))
        ):
            found["hand_zero_trace"].append((label, node.lineno))
    return found


_SYNTHETIC = '''
from functools import singledispatchmethod
class Op:
    @singledispatchmethod
    def _apply_impl(self, psi): ...
    @_apply_impl.register
    def _(self, psi: FullField):
        if isinstance(psi.interior, AngularFlux):
            return FullField(interior=psi.interior, boundary=AngularBoundarySourceSink.zeros(psi.boundary.space))
def helper(op, psi):
    if isinstance(psi, (ScalarField, np.ndarray)):
        return psi
'''


def test_the_census_catches_every_shape_it_exists_for():
    """vv #17 — the positive control: one synthetic module, four shapes."""
    found = _census(_SYNTHETIC, "synthetic.py")
    assert len(found["dispatcher"]) >= 2, found          # decorator + register
    assert sorted(f for f, _ in found["carrier_isinstance"]) == ["Op._", "helper"], found
    assert len(found["full_field_ctor"]) == 1, found
    assert len(found["hand_zero_trace"]) == 1, found


def test_the_carve_outs_are_exempt_only_where_declared():
    """The carve-out filter is keyed on (module, function): the same
    ``isinstance`` in another function of ``lift.py`` is an offence."""
    src = (
        "def admit_array(op, x):\n    return isinstance(x, FullField)\n"
        "def other(x):\n    return isinstance(x, FullField)\n"
        "class LegendreMomentTransfer:\n    def apply(self, m):\n        return isinstance(m, HarmonicMomentFlux)\n"
        "class Other:\n    def apply(self, m):\n        return isinstance(m, HarmonicMomentFlux)\n"
    )
    assert [f for f, _ in _census(src, "lift.py")["carrier_isinstance"]] == ["other", "LegendreMomentTransfer.apply", "Other.apply"]
    assert [f for f, _ in _census(src, "transfer.py")["carrier_isinstance"]] == ["admit_array", "other", "Other.apply"]


@pytest.mark.parametrize(
    "module", sorted(p.name for p in _PACKAGE.glob("*.py")),
)
def test_no_carrier_dispatch_in(module):
    found = _census((_PACKAGE / module).read_text(), module)
    assert found["dispatcher"] == [], f"{module}: a singledispatch table regrew: {found['dispatcher']}"
    assert found["carrier_isinstance"] == [], (
        f"{module}: an isinstance carrier parse outside the two admission "
        f"verbs — the ends select the body, the carrier is not parsed: "
        f"{found['carrier_isinstance']}"
    )
    assert found["full_field_ctor"] == [], (
        f"{module}: the composite is assembled outside lift_bulk_action "
        f"(a second spelling of the zero-trace emission): {found['full_field_ctor']}"
    )
    assert found["hand_zero_trace"] == [], (
        f"{module}: a boundary leaf's .zeros( spelled by name — the trace's "
        f"zero rides the operand's role partner: {found['hand_zero_trace']}"
    )


def test_the_package_is_the_denominator_the_census_walked():
    """State the denominator: the modules walked, and that the three named
    carve-out sites exist (a carve-out naming a dead function exempts nothing)."""
    modules = sorted(p.name for p in _PACKAGE.glob("*.py"))
    assert "lift.py" in modules and "angular_lift.py" in modules
    tree = ast.parse((_PACKAGE / "lift.py").read_text())
    defined = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert {f for m, _, f in _CARVE_OUTS if m == "lift.py"} <= defined and _LIFT_SITE[2] in defined
    lmt = ast.parse((_PACKAGE / "transfer.py").read_text())
    classes = {n.name for n in ast.walk(lmt) if isinstance(n, ast.ClassDef)}
    assert "LegendreMomentTransfer" in classes
