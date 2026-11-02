"""The traced memo's per-statement import-time pins (the cone store's increment B2, issue #689): precision and recall.

Spec: ``scratch/verification_tree/cone_store/b2/gate_spec.md`` (§2 the memo-level rows, §3 the declared boundaries).
The design under test: ``.claude/plans/verification_tree_refinement.md``, "B2 — the design", items 1-4.

Before B2 (B1, `e9df7ed9`) a run that executed any code of a first-party module pins the module's skeleton (``ModulePin``), so an edit to
ANY import-time statement of that module makes the manifest stale. B2 pins only the statements the run could observe
(``StatementsPin``). The rows here are one fixture each from :mod:`._statement_pin_cases`, by kind:

* ``precision`` rows are Red at B1 (`e9df7ed9`) (the skeleton moves) and green after B2;
* ``recall`` rows are green at B1 and must stay green after B2;
* ``under`` rows (the design's declared boundaries) assert the documented under-selection: Red at B1 (`e9df7ed9`), green after B2;
* ``over`` rows assert the documented over-selection: green at B1 and after;
* ``outcome`` rows (ruling R3 amended) change the import outcome: stale, green at B1 and after B2;
* ``QA-*`` rows are qa's 38 attack rows on the first design (``_statement_pin_qa_rows``), recall except PC2;
  the six findings F1-F6 of this spec are recall rows too, each closed by one ruling (R1, R2D, R4, R5).

Two instruments per row, each in a fresh process (the memo's source reads are kept per process): the trace
(``trace_call`` on the case's entry, with the package imported inside the call, as a generation imports it) and the
validation (``validate`` on the manifest's JSON after the edit). The validation process then runs the entry on the
edited tree, so a ``recall`` row also shows that its edit changes the answer (the edit is observable; X1).

A ``recall`` row's edit must be visible to an import-time statement pin ONLY: the row asserts that no ``DefPin`` of
the edited file is stale, so a row that a def pin would catch anyway cannot pass as statement-pin evidence.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ._statement_pin_cases import ALL_CASES, BY_ID, Case

pytestmark = pytest.mark.foundation

_THIS = "tests/gates/numerics/test_traced_memo_statement_pins.py"
_TABLE = f"{_THIS}::test_the_case_table_spans_every_kind"

REPO = Path(__file__).resolve().parents[3]

_TRACE = r'''
import json, sys
from orpheus.numerics.traced_memo import trace_call

import fxrun

value, manifest = trace_call(fxrun.run)
print("STMT-TRACE " + json.dumps({"value": repr(value), "manifest": manifest.to_json()}))
'''

_VALIDATE = r'''
import json, sys
from pathlib import Path
from orpheus.numerics.traced_memo import Manifest, validate

rows = json.loads(Path(sys.argv[1]).read_text())
reasons = list(validate(Manifest.from_json(rows)))
value = None
if sys.argv[2] == "run":
    import fxrun
    try:
        value = repr(fxrun.run())
    except Exception as error:  # the edit may make the entry raise: that is a changed answer too
        value = f"raised {type(error).__name__}"
print("STMT-VALID " + json.dumps({"reasons": reasons, "value": value}))
'''


def _environment(root: Path, on_path: bool = True) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("PYTEST_", "COV_"))}
    # The fixture root and this repository first; then whatever path the gate process received (a battery arm's
    # directory reaches the drivers this way, as the cone store's ``worker_environment`` passes it on).
    env["PYTHONPATH"] = os.pathsep.join([*([str(root)] if on_path else []), str(REPO),
                                         *filter(None, [os.environ.get("PYTHONPATH")])])
    env["OMP_NUM_THREADS"] = "1"
    # An edit of the same size within the second of the trace leaves a timestamp ``.pyc`` that the validating
    # process would run instead of the edited source (``[M]`` 2026-11-02: six recall rows read the old answer after
    # their edit when run in sequence, and passed alone). The pins read the files themselves; only the entry's
    # re-run needs this.
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def _run(root: Path, script: str, *arguments: str, tag: str, on_path: bool = True) -> dict:
    # The driver lives in a directory of its own: the script's directory is ``sys.path[0]``, and a pin's path is
    # relative to the first entry holding the file, so a driver beside the tree would pin ``tree/fx/...``.
    drivers = root.parent / "drivers"
    drivers.mkdir(exist_ok=True)
    driver = drivers / f"{root.name}.{tag.lower()}.py"
    driver.write_text(script)
    run = subprocess.run([sys.executable, "-O", str(driver), *arguments], cwd=root, env=_environment(root, on_path),
                         capture_output=True, text=True, timeout=300)
    line = next((l for l in run.stdout.splitlines() if l.startswith(f"STMT-{tag} ")), None)
    if run.returncode != 0 or line is None:
        pytest.fail(f"the {tag.lower()} driver exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    return json.loads(line[len(f"STMT-{tag} "):])


def build(root: Path, case: Case) -> Path:
    for relative, source in case.tree().items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    return root


def trace(root: Path) -> tuple[str, list]:
    found = _run(root, _TRACE, tag="TRACE")
    return found["value"], found["manifest"]


def validate(root: Path, manifest: list, *, run: bool) -> tuple[list[str], str | None]:
    rows = root.parent / f"{root.name}.manifest.json"
    rows.write_text(json.dumps(manifest))
    found = _run(root, _VALIDATE, str(rows), "run" if run else "no-run", tag="VALID")
    return found["reasons"], found["value"]


def apply(root: Path, case: Case) -> None:
    path = root / case.edit.file
    text = path.read_text()
    if case.edit.old is None:  # qa's rows: the file's whole new text
        if text == case.edit.new:
            pytest.fail(f"{case.id}: the edit leaves {case.edit.file} unchanged")
        path.write_text(case.edit.new)
        return
    if text.count(case.edit.old) != 1:
        pytest.fail(f"{case.id}: expected one {case.edit.old!r} in {case.edit.file}, found {text.count(case.edit.old)}")
    path.write_text(text.replace(case.edit.old, case.edit.new))


def pinned_files(manifest: list) -> set[str]:
    """The files the manifest pins by def or by module (or statement) grain, as their path strings."""
    return {str(row[1]).replace("\\", "/") for row in manifest if row[0] not in {"DefPin", "DataPin", "ChildPin"}
            and len(row) > 1 and isinstance(row[1], str) and row[1].endswith(".py")}


def about(reasons: list[str], file: str) -> tuple[list[str], list[str]]:
    """(every reason naming ``file``, the ones a ``DefPin`` gave): the def pin's reasons read ``function <path>:...``."""
    named = [r for r in reasons if file in r.replace("\\", "/")]
    return named, [r for r in named if r.startswith("function ")]


def not_by_def(reasons: list[str]) -> list[str]:
    """The reasons no ``DefPin`` gave (a module or statement pin, or the import-outcome check, whose wording the build
    chooses)."""
    return [r for r in reasons if not r.startswith("function ")]


def judged(tmp_path: Path, case: Case) -> tuple[list[str], list[str], str, str | None, list[str]]:
    """Trace, check the negative leg (no edit: valid) and that the edited file is pinned at all, edit, validate.
    Returns (reasons naming the edited file, those from a DefPin, the value before, the value after, every reason)."""
    root = build(tmp_path / "tree", case)
    before, manifest = trace(root)
    clean, _ = validate(root, manifest, run=False)
    if clean:
        pytest.fail(f"{case.id}: the manifest is stale with no edit (the negative leg): {clean}")
    if not any(path.endswith(case.edit.file) for path in pinned_files(manifest)):
        pytest.fail(f"{case.id}: the edited file {case.edit.file} is pinned by no module or statement pin, so a valid "
                    f"reading would be vacuous; pinned: {sorted(pinned_files(manifest))}")
    apply(root, case)
    reasons, after = validate(root, manifest, run=True)
    named, by_def = about(reasons, case.edit.file)
    return named, by_def, before, after, reasons


def _ids(kind: str) -> list[str]:
    return [case.id for case in ALL_CASES if case.kind == kind]


@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("case_id", _ids("precision"))
def test_precision_an_edit_the_run_cannot_observe_leaves_the_manifest_valid(tmp_path: Path, case_id: str) -> None:
    """An edit to a statement no code the run executed reads (by name, by literal, through a selected statement's
    reads) leaves the manifest valid. Red at B1 (`e9df7ed9`): the module skeleton moves. Mutation that reddens it after B2: the
    projection keeps every statement (``StatementsPin.names`` replaced by every identifier of the source)."""
    case = BY_ID[case_id]
    named, _, before, after, reasons = judged(tmp_path, case)
    assert not reasons, f"{case_id}: an edit the run cannot observe made the manifest stale: {reasons}"
    if after != before:
        pytest.fail(f"{case_id}: the fixture is wrong: the 'unobservable' edit changed the answer {before} -> {after}")


@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("case_id", _ids("recall"))
def test_recall_an_edit_the_run_observes_makes_the_manifest_stale(tmp_path: Path, case_id: str) -> None:
    """An edit to an import-time statement whose effect the run observes makes the manifest stale, through a
    module or statement pin (never only a def pin). Green at B1; the mutation that reddens each row after B2 is
    the spec's §2 table (drop the clause of the design the row names)."""
    case = BY_ID[case_id]
    named, by_def, before, after, reasons = judged(tmp_path, case)
    if case.before is not None:
        assert (before, after) == (case.before, case.after), f"{case_id}: the edit's effect is not the one written"
    else:
        assert before != after, f"{case_id}: the edit did not change the answer ({before}); the row is not a recall row"
    assert not_by_def(reasons), f"{case_id}: the run observes the edit ({before} -> {after}) and the manifest is still valid"
    assert not by_def, f"{case_id}: a def pin saw the edit ({by_def}); the row does not witness a statement pin"


@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("case_id", _ids("outcome"))
def test_outcome_a_changed_import_outcome_makes_the_manifest_stale(tmp_path: Path, case_id: str) -> None:
    """Ruling R3 (amended): the run records its import outcome (its first-party import set imported, and the warnings
    it emitted); when a byte of that set changes, the set is re-imported once and the pins are trusted only if the
    outcome is equal. An unselected statement that starts raising (I2), or a new import-time warning that leaves the
    answer unchanged (I1), makes the manifest stale. Green at B1 (the skeleton moves); the mutation that reddens it
    after B2: skip the outcome comparison (trust the statement pins alone)."""
    case = BY_ID[case_id]
    _, by_def, before, after, reasons = judged(tmp_path, case)
    assert (before, after) == (case.before, case.after), f"{case_id}: the edit's effect is not the one written"
    assert not_by_def(reasons), f"{case_id}: the import outcome changed and the manifest is still valid"
    assert not by_def, f"{case_id}: a def pin saw the edit ({by_def}); the row does not witness the outcome check"


@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("case_id", _ids("under"))
def test_boundary_the_documented_under_selection(tmp_path: Path, case_id: str) -> None:
    """A declared boundary of the design (item 4): the run observes the edit, and the manifest stays valid. The
    row asserts the documented limit so it is a measured fact. Red at B1 (`e9df7ed9`) (the module skeleton moves). If a later
    change closes the boundary the row reddens, and it is then re-posed as a recall row."""
    case = BY_ID[case_id]
    named, _, before, after, _ = judged(tmp_path, case)
    assert (before, after) == (case.before, case.after), f"{case_id}: the edit's effect is not the one written"
    assert not named, f"{case_id}: the boundary is closed (the manifest saw the edit: {named}); re-pose the row as recall"


@pytest.mark.rests_on(_TABLE)
@pytest.mark.parametrize("case_id", _ids("over"))
def test_boundary_the_documented_over_selection(tmp_path: Path, case_id: str) -> None:
    """A declared over-selection (the statement is the unit; names are not module-qualified; an effect statement
    keeps module grain): the run cannot observe the edit, and the manifest goes stale. Green at B1 and after B2."""
    case = BY_ID[case_id]
    named, _, before, after, _ = judged(tmp_path, case)
    assert before == after, f"{case_id}: the fixture is wrong: the edit changed the answer {before} -> {after}"
    assert named, f"{case_id}: the documented over-selection did not happen; the spec's boundary list is stale"


def test_the_case_table_spans_every_kind() -> None:
    """The table's population: every kind holds rows, and every id is unique (X2: an empty parametrisation
    reads as a green gate)."""
    kinds = {case.kind for case in ALL_CASES}
    assert kinds == {"precision", "recall", "outcome", "under", "over"}, kinds
    assert len(BY_ID) == len(ALL_CASES)
    assert sum(case.id.startswith("QA-") for case in ALL_CASES) == 38, "qa's 38 non-pytest rows"



# ── control 1 on the production tree, at the memo level ─────────────────────────────────────────────────────

_AXIS = "orpheus/mesh/axis.py"

_AXIS_TRACE = r'''
import json
import numpy as np
import orpheus.mesh.axis as axis
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import AxisMesh, FaceLabel
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.traced_memo import trace_call

axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 5), role_low=MirrorDeck(), role_high=MirrorDeck()),)
quad = Quadrature.gauss_legendre(n_ordinates=8)
value, manifest = trace_call(axis.face_outflow_ordinates, axes, FaceLabel(0, "max"), quad)
print("STMT-TRACE " + json.dumps({"value": repr(list(value)), "manifest": manifest.to_json()}))
'''

_SHADOW_VALIDATE = r'''
import json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[2])  # the shadow root: its orpheus/mesh/axis.py is the one the pin's path names now
from orpheus.numerics.traced_memo import Manifest, validate

print("STMT-VALID " + json.dumps({"reasons": list(validate(Manifest.from_json(json.loads(Path(sys.argv[1]).read_text())))),
                                  "value": None}))
'''


def _axis_judged(tmp_path: Path, old: str, new: str) -> list[str]:
    """Trace ``face_outflow_ordinates`` on the production tree (``orpheus/mesh/axis.py`` imported BEFORE the trace,
    so only the def that ran pins the file), write an edited copy of ``axis.py`` under a shadow root that the
    validating process puts first on ``sys.path`` (the production file is never touched), and return the reasons
    (every reason, whatever pin gave it)."""
    root = tmp_path / "axis"
    root.mkdir()
    found = _run(root, _AXIS_TRACE, tag="TRACE")
    manifest = found["manifest"]
    if not any(row[0] == "DefPin" and row[1].endswith(_AXIS) and row[2] == "face_outflow_ordinates" for row in manifest):
        pytest.fail(f"the trace did not run face_outflow_ordinates: {[r for r in manifest if r[0] == 'DefPin'][:5]}")
    shadow = tmp_path / "shadow"
    target = shadow / _AXIS
    target.parent.mkdir(parents=True)
    text = (REPO / _AXIS).read_text()
    if text.count(old) != 1:
        pytest.fail(f"expected one {old!r} in {_AXIS}, found {text.count(old)}")
    rows = tmp_path / "axis.manifest.json"
    rows.write_text(json.dumps(manifest))
    target.write_text(text)
    clean = _run(root, _SHADOW_VALIDATE, str(rows), str(shadow), tag="VALID")["reasons"]
    if clean:
        pytest.fail(f"an unedited shadow copy of axis.py made the manifest stale (the negative leg): {clean}")
    target.write_text(text.replace(old, new))
    return _run(root, _SHADOW_VALIDATE, str(rows), str(shadow), tag="VALID")["reasons"]


def test_control_1_memo_the_read_constant_is_pinned(tmp_path: Path) -> None:
    """Control 1 (recall leg, memo level). ``_OUTFLOW_COSINE_TOL`` is read by ``face_outflow_ordinates``: editing it
    makes the manifest stale through a module or statement pin (the def's own digest names the constant, not its
    value). Green at B1 and after B2."""
    reasons = _axis_judged(tmp_path, "_OUTFLOW_COSINE_TOL = 1e-15", "_OUTFLOW_COSINE_TOL = 2e-15")
    named = [r for r in reasons if _AXIS in r.replace("\\", "/")]
    assert not_by_def(named) and not [r for r in named if r.startswith("function ")], reasons


def test_control_1_memo_an_unread_constant_of_the_same_module_is_not_pinned(tmp_path: Path) -> None:
    """Control 1 (precision leg, memo level). ``_ENDPOINT_SUFFIX`` (``axis.py`` line ~93) is bound at import and read
    by no code ``face_outflow_ordinates`` runs: editing its value leaves the manifest valid. Red at B1 (`e9df7ed9`) (the skeleton
    of ``axis.py`` moves). ``[M]`` ``ta/probe_axis_e.py`` (the first design's prototype): ``_ENDPOINT_SUFFIX`` is not in
    the run's closure names (|E'| = 73), while ``_OUTFLOW_COSINE_TOL`` and ``_OUTWARD_ENDPOINTS`` are. Every reason
    counts, so a misfiring import-outcome check (whose reason names no file) cannot pass the row."""
    reasons = _axis_judged(tmp_path, '_ENDPOINT_SUFFIX = {"min": "min", "max": "max", "outer": "max"}',
                         '_ENDPOINT_SUFFIX = {"min": "min", "max": "max", "outer": "outer"}')
    assert not reasons, f"an edit to a constant the run never reads made the manifest stale: {reasons}"


@pytest.mark.rests_on(f"{_THIS}::test_the_case_table_spans_every_kind")
def test_a_run_that_evaluates_type_alias_values_is_pinned_not_refused(tmp_path: Path) -> None:
    """A run that evaluates a PEP 695 alias value, at module level and in a class body, is pinned: the alias value's
    code object (named ``V``, or ``W`` inside ``C``) is module-level code whose ``type`` statement pins it, so
    ``trace_call`` returns a manifest instead of refusing with ``NoDefAtLine``, and an edit to the alias value makes
    the manifest stale. Reddening mutation: drop ``ast.TypeAlias`` from the source index's scope walk
    (``traced_memo._parsed_index``): the trace raises ``Unpinnable(NoDefAtLine)`` (``[M]`` before the fix, row QA-R5)."""
    root = tmp_path / "tree"
    (root / "fx").mkdir(parents=True)
    (root / "fx" / "__init__.py").write_text("")
    (root / "fx" / "aliases.py").write_text("type PsiV = int\n\n\nclass PsiC:\n    type PsiW = float\n")
    (root / "fxrun.py").write_text("def run():\n    from fx import aliases\n"
                                   "    return (aliases.PsiV.__value__.__name__, aliases.PsiC.PsiW.__value__.__name__)\n")
    found = _run(root, _TRACE, tag="TRACE")
    assert found["value"] == "('int', 'float')", found["value"]
    assert any(row[1].endswith("fx/aliases.py") for row in found["manifest"] if row[0] != "DefPin"), "aliases.py is unpinned"
    path = root / "fx" / "aliases.py"
    path.write_text(path.read_text().replace("type PsiW = float", "type PsiW = complex"))
    reasons, after = validate(root, found["manifest"], run=True)
    assert after == "('int', 'complex')", after
    assert [r for r in reasons if "fx/aliases.py" in r and not r.startswith("function ")], reasons


def test_the_import_outcome_is_not_served_across_a_change_of_the_memos_own_source(tmp_path: Path,
                                                                                   monkeypatch: pytest.MonkeyPatch) -> None:
    """The import outcome is kept by the digest of its inputs, the memo's own source among them (it holds the
    measuring program): an outcome kept under one memo digest is measured again under another, and served under the
    same one. Instrument: the files in the outcome store (one per measurement), with a counter on the measuring
    subprocess. Reddening mutation: drop ``_INDEX_MEMO_DIGEST`` from the key of ``traced_memo._import_outcome`` (the
    second memo digest is served the first one's outcome: 1 file, 1 measurement)."""
    import subprocess as subprocess_module

    from orpheus.numerics import traced_memo as memo

    (tmp_path / "pkgo").mkdir()
    (tmp_path / "pkgo" / "__init__.py").write_text("")
    (tmp_path / "pkgo" / "m.py").write_text("OMICRON_VALUE = 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(memo, "_IMPORT_OUTCOMES", tmp_path / "outcomes")
    measured: list[int] = []
    real_run = subprocess_module.run

    def counting(*args, **kwargs):
        measured.append(1)
        return real_run(*args, **kwargs)

    monkeypatch.setattr(memo.subprocess, "run", counting)
    modules = (("pkgo/m.py", memo._digest((tmp_path / "pkgo" / "m.py").read_bytes()), "pkgo.m"),)
    monkeypatch.setattr(memo, "_INDEX_MEMO_DIGEST", "a" * 64)
    first = memo._import_outcome(modules, [str(tmp_path)]).digest
    monkeypatch.setattr(memo, "_INDEX_MEMO_DIGEST", "b" * 64)
    second = memo._import_outcome(modules, [str(tmp_path)]).digest
    monkeypatch.setattr(memo, "_INDEX_MEMO_DIGEST", "a" * 64)
    again = memo._import_outcome(modules, [str(tmp_path)]).digest
    kept = sorted((tmp_path / "outcomes").glob("*.json"))
    assert len(measured) == 2 and len(kept) == 2, f"{len(measured)} measurements, {len(kept)} kept outcomes"
    assert first == second == again, "one import set measured twice gave two outcomes"


#: A package found by a meta-path finder and not through ``sys.path`` (an editable install's route, or a run started
#: outside the repository): its files are pinned by ABSOLUTE path (qa's build review, F1).
_FINDER = r'''
import sys, importlib.machinery as _machinery
_ROOT = sys.argv[1]


class _Editable:
    @classmethod
    def find_spec(cls, name, path=None, target=None):
        if name.split(".")[0] in ("pkgf", "fxrunf"):
            return _machinery.PathFinder.find_spec(name, [_ROOT] if path is None else path)
        return None


sys.meta_path.append(_Editable)
'''
_FINDER_TRACE = _FINDER + r'''
import json
from orpheus.numerics.traced_memo import Unpinnable, trace_call
import fxrunf
try:
    value, manifest = trace_call(fxrunf.run)
    found = {"value": repr(value), "manifest": manifest.to_json(), "refusal": None}
except Unpinnable as refused:
    found = {"value": None, "manifest": None, "refusal": [type(refused.refusal).__name__, str(refused.refusal)]}
print("STMT-TRACE " + json.dumps(found))
'''


@pytest.mark.rests_on(f"{_THIS}::test_outcome_a_changed_import_outcome_makes_the_manifest_stale")
def test_a_package_found_by_the_runs_own_finder_is_measured_from_its_real_root(tmp_path: Path) -> None:
    """qa's F1 (build review), re-posed after the roots repair: a package found by a meta-path finder installed in the
    run's OWN process, so its files are pinned by absolute path and no ``sys.path`` entry names its root. The import
    outcome is measured with each module's root derived from its REAL file (the directory that, joined with its dotted
    name's path, gives the file), so the fresh measuring interpreter imports it: the trace is pinned (no refusal, an
    import-outcome pin present), and the S2k′ edit (a def deleted that only an unread table names, so the import now
    raises) makes the manifest stale. Reddening mutations: roots from the pin path (the package's root is wrong: the
    trace is refused, ``UnmeasurableImport``); drop the measurement-errors refusal together with the roots (the edit
    reads valid). History: before the roots repair this row asserted the refusal, the only safe answer then."""
    root = tmp_path / "abs"
    (root / "pkgf").mkdir(parents=True)
    (root / "pkgf" / "__init__.py").write_text("")
    (root / "pkgf" / "c.py").write_text("C = 0\n")
    (root / "pkgf" / "b.py").write_text("from . import c\nVALUE = 1\n\n\ndef old():\n    return 0\n\n\n_LEGACY = {'old': old}\n")
    (root / "fxrunf.py").write_text("def run():\n    from pkgf import b\n    return b.VALUE\n")
    found = _run(root, _FINDER_TRACE, str(root), tag="TRACE", on_path=False)
    assert found["refusal"] is None, f"the trace was refused: {found['refusal']}"
    manifest = found["manifest"]
    assert found["value"] == "1", found["value"]
    assert any(row[0] == "ImportOutcomePin" for row in manifest), "no import-outcome pin"
    if not any(row[0] == "StatementsPin" and os.path.isabs(row[1]) and row[1].endswith("pkgf/b.py") for row in manifest):
        pytest.fail(f"the route does not pin by absolute path: {[r[:2] for r in manifest if r[0] == 'StatementsPin']}")
    clean, _ = validate(root, manifest, run=False)
    assert not clean, f"stale with no edit (the negative leg): {clean}"
    path = root / "pkgf" / "b.py"
    path.write_text(path.read_text().replace("def old():\n    return 0\n", ""))
    reasons, _ = validate(root, manifest, run=False)
    assert not_by_def(reasons), "the import now raises and the manifest is still valid"


def test_a_kept_set_missing_once_is_read_when_it_appears_and_a_corrupt_one_is_rewritten(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """qa's build review K1 and K2. K1: a kept set missing at a process's first read is found by a later read in
    the same process once a generation keeps it (a miss is not remembered). K2: a kept-set file whose bytes do not
    digest to its name is rewritten by the next ``_keep_set`` of that value, and then read. Reddening mutations: K1,
    remember a miss (cache ``None`` in ``_read_kept_set``); K2, write only when the file is absent (``_keep_set``)."""
    from orpheus.numerics import traced_memo as memo

    monkeypatch.setattr(memo, "_KEPT_SETS", tmp_path / "sets")
    monkeypatch.setattr(memo, "_KEPT_SETS_READ", {})
    value = ["kappa_a", "kappa_b"]
    name = memo._digest(json.dumps(value, separators=(",", ":")).encode())
    assert memo._read_kept_set(name) is None, "a set never kept was read"
    assert memo._keep_set(value) == name
    assert memo._read_kept_set(name) == tuple(value), "K1: a set missing once stays missing in this process"
    other = ["kappa_c"]
    other_name = memo._digest(json.dumps(other, separators=(",", ":")).encode())
    (tmp_path / "sets").mkdir(exist_ok=True)
    (tmp_path / "sets" / f"{other_name}.json").write_bytes(b'["truncated"')
    assert memo._read_kept_set(other_name) is None, "a corrupt set was served"
    memo._keep_set(other)
    assert memo._read_kept_set(other_name) == tuple(other), "K2: a corrupt set file was not rewritten"


def test_the_blas_tag_is_pinned_by_its_def_and_by_the_statement_it_reads() -> None:
    """qa's final F3. The BLAS name keys an entry (``platform_tag``); it is extracted by ``_blas``, a def recorded on
    every call, from ``_NUMPY_BUILD_CONFIG``, the module statement that caches numpy's build configuration and that
    ``_blas`` reads by name. So an edit to ``_blas``'s body moves its def's digest, and an edit to that statement moves
    the projection on ``_blas``'s code names; and a traced call of ``platform_tag`` made after the cache is warm still
    pins ``_blas`` (its code starts on every call). Reddening mutation: cache ``_blas`` itself (``functools.cache``):
    the warm call pins no ``_blas`` def."""
    from orpheus.numerics import _import_time
    from orpheus.numerics import traced_memo as memo

    source = Path(memo.__file__).read_bytes()
    index = memo._source_index(source)
    names = frozenset(_import_time.code_names(memo._blas.__code__))
    assert "_NUMPY_BUILD_CONFIG" in names, sorted(names)
    body = source.replace(b'.get("name", "unknown") if isinstance', b'.get("name", "openblas") if isinstance', 1)
    statement = source.replace(b'functools.partial(np.show_config, mode="dicts")',
                               b'functools.partial(np.show_config, mode="stdout")', 1)
    assert body != source and statement != source, "an edit target moved: re-aim the row"
    assert memo._source_index(body).digest("_blas", 0) != index.digest("_blas", 0)
    assert memo._source_index(statement).projection(names) != index.projection(names)
    memo.platform_tag()  # warm every cache
    _, manifest = memo.trace_call(memo.platform_tag)
    assert any(p.qualname == "_blas" and p.path.endswith("traced_memo.py") for p in manifest.functions), \
        sorted(p.qualname for p in manifest.functions)
