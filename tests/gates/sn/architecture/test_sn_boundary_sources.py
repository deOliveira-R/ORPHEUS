r"""Source-level gates of the S2b carve: one body per builder, and no string dispatch in SN (wave 2, slice S2b).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.6 (decision 45), section 6.2 (S2b gates 11 and 19) and
section 7.2 (row R20). These gates READ FILES (``pathlib`` and ``ast``) and import no production module, so at the base
tip they fail for their NAMED reason (the retired builders still stand, the registry is still there) and not for the
``BC`` import chain that reddens every gate that executes SN. Each scanner is validated by a planted positive control
before its zero is believed (``instrument-doctrine`` X1), and each reports its input count (X2).

Gate 19, ``test_one_body_per_builder``
======================================

An AST scan of ``orpheus/sn/boundary/`` finds ``_deck_kernel``, ``_narrowed_zero_operator``,
``_checked_angular_average`` and ``_specular_crossed_diagnosis`` defined 0 times, ``_attenuated_kernel_operator``
defined once (in ``realized_operators.py``), ``realized_operators.py`` importing nothing from ``realizer.py``, and
``realizer.py`` keeping at most the two module-level helpers ``_outflow_restriction`` and ``_partner_outflow_restriction``
beside its classes (S2b deletes either when the typed operators' spaces make it dead: a subset, not an equality). Why each: the typed operators OWN the one body of each builder (a window in which both a builder and
its typed operator stood would be two bodies of one quantity, and the Lambertian divisor change of R-11 would land in one
copy); and an import in the other direction would red the WF2.4 intra-package cycle gate.

Mutations for the reviewers: a copy of ``realizer.py`` that keeps ``_narrowed_zero_operator`` (reddens the zero-definition
row); a copy of ``realized_operators.py`` with ``from orpheus.sn.boundary.realizer import _outflow_restriction`` (reddens
the import row, whose scanner the planted control proves can see it); a second ``_attenuated_kernel_operator`` in
``realizer.py`` (reddens the defined-once row).

Gate 11, ``test_no_string_dispatch_in_sn``
==========================================

No ``BOUNDARY_OPERATOR_REGISTRY`` on ``SNProblem`` or ``TransportMethod``, no ``.kind`` read (the sweep schedule's own
``schedule.kind`` is a different concept and is the one allowed read), and none of the welded questions the typed
operators now answer by type: ``law_permutes_ordinates``, ``response_kernel``, ``geometry_map``, ``.amplitude``,
``.law``; and no consumer imports the retiring shim ``bound_compat``, the ``SNBoundaryRealizer`` or the retired ``BC``
(``orpheus/sn/boundary/bound_compat.py`` itself stays until WF2.4, section 6.2; it is excluded from the scan and from
nothing else).
"""

from __future__ import annotations

import ast
import os
from collections import Counter
from pathlib import Path

import pytest

import orpheus

pytestmark = pytest.mark.foundation

_THIS = "tests/gates/sn/architecture/test_sn_boundary_sources.py"
#: The tree under test: the imported ``orpheus``'s parent, NOT resolved (a symlinked tree must stay itself).
ROOT = Path(os.path.abspath(orpheus.__file__)).parent
BOUNDARY = ROOT / "sn" / "boundary"
RETIRED_BUILDERS = ("_deck_kernel", "_narrowed_zero_operator", "_checked_angular_average", "_specular_crossed_diagnosis")
KEPT_HELPERS = frozenset({"_outflow_restriction", "_partner_outflow_restriction"})
TYPED_OPERATORS = ("RealizedBoundary", "ZeroMap", "ScaledOrdinatePermutation", "RankOneReemission", "RealizedSum")


def definitions(source: str) -> Counter[str]:
    """How many times each function or class name is DEFINED in ``source`` (any depth)."""
    return Counter(n.name for n in ast.walk(ast.parse(source)) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)))


def imports_realizer(source: str) -> list[str]:
    """The import statements of ``source`` that reach ``orpheus.sn.boundary.realizer`` (absolute or relative)."""
    found = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.ImportFrom):
            module = n.module or ""
            names = {a.name for a in n.names}
            if module.endswith("realizer") or (n.level and module == "realizer") or ("realizer" in names and module.endswith("boundary")) or (n.level and module == "" and "realizer" in names):
                found.append(ast.unparse(n))
        elif isinstance(n, ast.Import):
            found += [ast.unparse(n) for a in n.names if a.name.endswith(".realizer")]
    return found


def module_level_functions(source: str) -> set[str]:
    return {n.name for n in ast.parse(source).body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}


def _source(name: str) -> str:
    path = BOUNDARY / name
    assert path.exists(), f"{path} does not exist"
    return path.read_text()


# ═════════════════════════════════════════════════════════════════════════════
# Gate 19
# ═════════════════════════════════════════════════════════════════════════════


def test_the_scanners_see_what_they_must() -> None:
    r"""The positive controls of gate 19 (a negative reading is evidence only after the instrument has read a positive):
    a planted source defining the five builders reads one definition each; a planted ``from
    orpheus.sn.boundary.realizer import _outflow_restriction`` (and its relative and module spellings) is found; a
    planted source with an extra module-level helper is found by the module-level scan. Passes today; reddens on a
    scanner that reads nothing.
    """
    planted = "def _deck_kernel(): pass\ndef _narrowed_zero_operator(): pass\ndef _checked_angular_average(): pass\ndef _specular_crossed_diagnosis(): pass\ndef _attenuated_kernel_operator(): pass\nclass K:\n    def _deck_kernel(self): pass\n"
    counts = definitions(planted)
    assert counts["_deck_kernel"] == 2 and all(counts[n] == 1 for n in RETIRED_BUILDERS[1:] + ("_attenuated_kernel_operator",))
    for spelling in (
        "from orpheus.sn.boundary.realizer import _outflow_restriction",
        "from .realizer import _outflow_restriction",
        "from orpheus.sn.boundary import realizer",
        "import orpheus.sn.boundary.realizer",
    ):
        assert imports_realizer(spelling), f"the scanner cannot see {spelling!r}"
    assert not imports_realizer("from orpheus.sn.boundary.angular import PartialCurrentOperator")
    assert module_level_functions("def a(): pass\nclass C:\n    def b(self): pass\n") == {"a"}


@pytest.mark.rests_on(f"{_THIS}::test_the_scanners_see_what_they_must")
def test_one_body_per_builder() -> None:
    r"""The retired builders are defined 0 times in ``orpheus/sn/boundary/``; ``_attenuated_kernel_operator`` is defined
    once, in ``realized_operators.py``; ``realized_operators.py`` imports nothing from ``realizer.py``; ``realizer.py``
    keeps at most the two helpers ``_outflow_restriction`` and ``_partner_outflow_restriction`` at module level
    (the Gamma-to-Gamma_+ restrictions the realizer composes around the typed operators; S2b deletes either if the typed
    operators' spaces make it dead, hence a subset); the typed operators are defined in ``realized_operators.py`` and exported by the
    package ``__init__``.

    Red input at the base tip (``[M]`` the census of the count, 2026-10-23): each of the four builders and
    ``_attenuated_kernel_operator`` is defined once in ``realizer.py``, ``realized_operators.py`` does not exist.
    """
    files = sorted(BOUNDARY.glob("*.py"))
    total = Counter()
    where: dict[str, list[str]] = {}
    for path in files:
        for name, count in definitions(path.read_text()).items():
            total[name] += count
            where.setdefault(name, []).append(path.name)
    assert len(files) >= 3, f"scanned {len(files)} files of orpheus/sn/boundary/"
    for name in RETIRED_BUILDERS:
        assert total[name] == 0, f"{name} is still defined {total[name]} time(s), in {where[name]}"
    assert total["_attenuated_kernel_operator"] == 1 and where["_attenuated_kernel_operator"] == ["realized_operators.py"], f"_attenuated_kernel_operator: {where.get('_attenuated_kernel_operator')}"

    operators = _source("realized_operators.py")
    assert not imports_realizer(operators), f"realized_operators.py imports the realizer: {imports_realizer(operators)}"
    defined = definitions(operators)
    for name in TYPED_OPERATORS:
        assert defined[name] == 1, f"{name} is defined {defined[name]} times in realized_operators.py"
    realizer = _source("realizer.py")
    extra = module_level_functions(realizer) - KEPT_HELPERS
    assert not extra, f"realizer.py keeps module-level functions beyond the two Gamma-to-Gamma_+ restrictions: {sorted(extra)}"
    exported = {
        str(elt.value)
        for n in ast.parse(_source("__init__.py")).body
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "__all__" for t in n.targets)
        for elt in getattr(n.value, "elts", [])
        if isinstance(elt, ast.Constant)
    }
    assert set(TYPED_OPERATORS) <= exported, f"orpheus.sn.boundary exports {sorted(exported)}, missing {sorted(set(TYPED_OPERATORS) - exported)}"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 11
# ═════════════════════════════════════════════════════════════════════════════

#: Names whose mere use is a string-dispatch or a welded question about the declared law (section 2.6: "five SN functions
#: ask the law's factors"). ``amplitude`` and ``law`` and ``kind`` are ATTRIBUTE reads; the others are any name or attribute.
FORBIDDEN_NAMES = frozenset({
    "BOUNDARY_OPERATOR_REGISTRY", "law_permutes_ordinates", "response_kernel", "geometry_map", "from_problem_laws",
    "realize_boundary_law", "_BoundBoundaryOperator", "SNBoundaryRealizer", "BoundaryTraceLaw", "face_laws", "bound_compat",
})
FORBIDDEN_ATTRIBUTES = frozenset({"amplitude", "law", "kind"})
#: Allowed reads of a forbidden attribute (the receiver's source, ending): the sweep schedule's own ``kind``.
ALLOWED_RECEIVERS = ("schedule",)


def forbidden_uses(source: str) -> list[str]:
    """The uses of a forbidden name or attribute read in ``source`` (a ``bound_compat`` import included)."""
    hits = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Name) and n.id in FORBIDDEN_NAMES:
            hits.append(f"{n.lineno}: name {n.id}")
        elif isinstance(n, ast.Attribute):
            if n.attr in FORBIDDEN_NAMES:
                hits.append(f"{n.lineno}: attribute {n.attr}")
            elif n.attr in FORBIDDEN_ATTRIBUTES and not ast.unparse(n.value).endswith(ALLOWED_RECEIVERS):
                hits.append(f"{n.lineno}: attribute read {ast.unparse(n)}")
        elif isinstance(n, ast.ImportFrom):
            imported = {a.name for a in n.names} | {(n.module or "").rsplit(".", 1)[-1]}
            if imported & FORBIDDEN_NAMES or (n.module == "orpheus.geometry" and "BC" in imported):
                hits.append(f"{n.lineno}: import {ast.unparse(n)}")
        elif isinstance(n, ast.FunctionDef) and n.name in FORBIDDEN_NAMES:
            hits.append(f"{n.lineno}: def {n.name}")
        elif isinstance(n, ast.Constant) and False:  # pragma: no cover
            pass
    return hits


def _sn_files() -> list[Path]:
    sn = [p for p in sorted((ROOT / "sn").rglob("*.py")) if p.name != "bound_compat.py"]
    transport = [ROOT / "transport" / "method.py", ROOT / "transport" / "fields" / "angular_boundary_source_sink.py", ROOT / "transport" / "mesh" / "material_mesh.py"]
    return sn + transport


def test_the_dispatch_scanner_sees_what_it_must() -> None:
    r"""The positive controls of gate 11: each planted use is found, the sweep schedule's ``kind`` is not. Passes today."""
    for planted in (
        "x = problem.bc['xmin'].kind",
        "ok = law_permutes_ordinates(law)",
        "a = law.response_kernel.amplitude",
        "class P:\n    BOUNDARY_OPERATOR_REGISTRY = {}",
        "from orpheus.geometry import BC",
        "from orpheus.sn.boundary.bound_compat import _BoundBoundaryOperator",
        "face = bc[f].geometry_map.domain_face(f)",
        "y = op.law",
    ):
        assert forbidden_uses(planted), f"the scanner cannot see {planted!r}"
    assert not forbidden_uses("k = self.schedule.kind")
    assert not forbidden_uses("from orpheus.geometry import CoordSystem")


@pytest.mark.rests_on(f"{_THIS}::test_the_dispatch_scanner_sees_what_it_must")
def test_no_string_dispatch_in_sn() -> None:
    r"""No SN file (and no transport consumer of the binding) holds a forbidden name, a forbidden attribute read or a
    ``bound_compat``/``BC`` import: the five welded questions are answered by the typed operators' properties, the
    registry is the table, and the declaration is read from the binding. ``bound_compat.py`` itself is the one file
    excluded (it retires at WF2.4 and nothing imports it).

    Red input at the base tip: ``BOUNDARY_OPERATOR_REGISTRY`` stands on ``SNProblem`` (``problem.py:192``) and on
    ``TransportMethod`` (``method.py``), ``from orpheus.geometry import BC`` stands in ``method.py`` and ``solver.py``,
    and the five welded questions are asked in ``dsa.py``, ``boundary.py``, ``sweep_schedule.py``, ``solver.py`` and
    ``coupled_system.py``. Mutations: a copy of ``dsa.py`` that reads ``.law.response_kernel.is_zero`` again; a
    ``BOUNDARY_OPERATOR_REGISTRY`` class attribute restored on ``SNProblem``.
    """
    files = _sn_files()
    assert len(files) >= 30 and all(p.exists() for p in files), f"scanned {len(files)} files"
    offences = {str(p.relative_to(ROOT)): forbidden_uses(p.read_text()) for p in files}
    offences = {k: v for k, v in offences.items() if v}
    assert not offences, f"{sum(map(len, offences.values()))} uses in {len(offences)} of {len(files)} files: " + "; ".join(f"{k}: {v[:3]}" for k, v in sorted(offences.items())[:6])
