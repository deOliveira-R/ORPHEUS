"""Gates G2 and G3 of unit ``memo-digest`` (specification part B 4.1): the traced memo is in no old home, and no file spells the old home.

G2a: ``orpheus/numerics/`` holds no file named for a memo module and no class named ``Recording`` or ``TracedMemo``.
G2b: the engine's module rows are read red by a planted import and a planted leftover file (the positive control; the 21 rows of
``tests/gates/_retired/rows_w0_memo.py`` are asserted green on the real tree by ``test_retirements.py``, the one engine).
G3: the old path is spelled nowhere as TEXT (the three patterns of 2.2, the test-location pattern, the segment form and the
path-expression form), over every tracked file of ``orpheus tools tests examples`` except ``tests/gates/_retired/`` and this file,
and over ``docs/`` with an exemption rule for a git path at a past commit (``git show <hash>:<path>`` must keep working).
This file builds the old spellings from pieces, so it does not spell what it forbids. Every walk takes a root (the planted legs).

First red in the tree at base (measured, see the commit body): G2a 3 files and 2 classes in the old home; G3 324 lines in 54 files
outside ``.claude/`` (the docs half reads 153 lines in 7 files and is red until the builder and the orchestrator re-point or exempt
them: ``_PAST_ROWS`` is the orchestrator's tense table, empty at the gate's landing).
"""
from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[2]
_OLD = "orpheus" + "." + "numerics"
_OLD_DIR = "numerics"
_MEMO_NAMES = ("traced_memo", "_traced_memo_boot", "_import_time")
_M = r"(?:traced_memo|_traced_memo_boot|_import_time)"
_TESTS = (
    "test_traced_memo_apparatus", "test_traced_memo_compiled_bytes", "test_traced_memo_data", "test_traced_memo_findings",
    "test_traced_memo_listing_ignores_bytecode", "test_traced_memo_manifest", "test_traced_memo_nested_recording",
    "test_traced_memo_ownership", "test_traced_memo_process", "test_traced_memo_records_each_file", "test_traced_memo_source_index",
    "test_traced_memo_statement_pins", "test_traced_memo_template", "test_traced_memo_validation",
    "_traced_memo_api", "_traced_memo_synthetic", "_statement_pin_cases", "_statement_pin_qa_rows",
)
_PATTERNS = {
    "P1": re.compile(re.escape(_OLD) + r"\." + _M + r"\b"),
    "P2": re.compile(_OLD_DIR + "/" + _M + r"(\.py)?\b"),
    "P3": re.compile(r"from\s+" + re.escape(_OLD) + r"\s+import\s+[^\n]*\b" + _M + r"\b"),
    "TESTPATH": re.compile(r"gates[/.]" + _OLD_DIR + r"[/.](?:" + "|".join(_TESTS) + r")\b"),
}
_P4 = re.compile(r'"orpheus",\s*"' + _OLD_DIR + '"')
_GIT_PAST = re.compile(r"[0-9a-f]{7,40}:orpheus/" + _OLD_DIR + "/")
SCOPE = ("orpheus", "tools", "tests", "examples")

#: The orchestrator's tense table for ``docs/``: ``(file, exact line text)`` of a line that names the old path as a fact about the
#: past and stays. The two ``git show <hash>:path`` lines are exempt by rule (``_GIT_PAST``); the four rows are a dated case and three history or git-at-a-past-commit lines.
_PAST_ROWS: tuple[tuple[str, str], ...] = (
    ('docs/development/evidence/vv-anti-patterns.md',
     "**Surprise.** #405 P3 step 2c (C1) memoised one group's characteristic transport block on disk (`group_transport`, a traced memo, `orpheus/numerics/traced_memo.py`, the path at the time). A miss generates the block in a fresh interpreter from the source ON DISK, and a hit serves what that source computed, so a monkeypatch or a `-p` plugin rebinding a kernel symbol in the pytest process never reaches a block read through the memo. `[M]` 2026-10-14, the test-architect's C1 battery (`scratch/characteristic_architecture/p3/c1/gate_spec.md` §4.3, §4.4, §6): with the system's blocks built at construction (arm `eager-blocks`), the old spy on `LineRule.of` reddened 4 construction rows and stayed green on the other 4 (F5 on three charts, D11), because those blocks' `LineRule.of` ran in the generating interpreter; re-pointed onto the system module's binding of `group_transport`, the same arm reddened all 8. And B1, which compares a served block with one computed under `bypass()`, reddens under ANY in-process kernel mutation, because only the bypassed side sees it: the control arm `control-kernel-garbage` reddened B1 on 4 of 4 fixtures and B6b."),
    ('docs/theory/verification/cone_store.rst',
     '       the statement model is ``orpheus/numerics/_import_time.py``; the'),
    ('docs/theory/verification/error_catalog.rst',
     '   HEAD:orpheus/numerics/traced_memo.py`` on ``f0d47585``). Full account:'),
    ('docs/theory/verification/reference_cache.rst',
     '       ``orpheus/numerics/_import_time.py``, a pinned memo source; the'),
)


def _tracked(*roots: str) -> list[str]:
    out = subprocess.run(["git", "ls-files", "-z", *roots], capture_output=True, text=True, cwd=REPO, check=True).stdout
    return [f for f in out.split("\0") if f]


def _text(root: Path, rel: str) -> str | None:
    try:
        return (root / rel).read_text(encoding="utf-8")
    except (UnicodeDecodeError, FileNotFoundError):
        return None


# ---------------------------------------------------------------- G2a: the old home holds no memo module and no memo class
def _old_home_findings(numerics: Path) -> list[str]:
    if not any(numerics.glob("*.py")):
        raise AssertionError(f"input count 0: {numerics} holds no module (a vacuous walk)")
    found = [f"file {n}.py" for n in _MEMO_NAMES if (numerics / f"{n}.py").exists()]
    for path in sorted(numerics.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        found += [f"class {c.name} in {path.name}" for c in ast.walk(tree) if isinstance(c, ast.ClassDef) and c.name in ("Recording", "TracedMemo")]
    return found


def test_the_old_numerics_home_holds_no_memo_module_and_no_memo_class() -> None:
    numerics = REPO / "orpheus" / _OLD_DIR
    assert (numerics / "instance_cache.py").is_file(), "positive control: the walk's population holds the module that stays"
    assert _old_home_findings(numerics) == []


def test_the_old_home_walk_finds_a_planted_file_and_a_planted_class(tmp_path: Path) -> None:
    (tmp_path / "instance_cache.py").write_text("x = 1\n")
    assert _old_home_findings(tmp_path) == []
    (tmp_path / "traced_memo.py").write_text("x = 1\n")
    (tmp_path / "other.py").write_text("class Recording:\n    pass\n")
    assert sorted(_old_home_findings(tmp_path)) == ["class Recording in other.py", "file traced_memo.py"]
    with pytest.raises(AssertionError, match="input count 0"):
        _old_home_findings(tmp_path / "empty")


# ---------------------------------------------------------------- G2b: the module rows' positive control
def test_the_module_rows_are_red_for_a_planted_import_and_a_planted_leftover_file(tmp_path: Path) -> None:
    """Positive control of G2b (runs at base): a planted import statement, a planted import_module string and a planted leftover file."""
    from tests.gates._retired import RetiredModule, _walk

    (tmp_path / "tools").mkdir()
    (tmp_path / "tools" / "x.py").write_text(
        f"from {_OLD} import traced_memo\nimport tests.gates.{_OLD_DIR}.test_traced_memo_data\n"
        f"import importlib\nimportlib.import_module({_OLD + '._import_time'!r})\n")
    (tmp_path / "orpheus" / _OLD_DIR).mkdir(parents=True)
    (tmp_path / "orpheus" / _OLD_DIR / "_traced_memo_boot.py").write_text("x = 1\n")
    rows = [RetiredModule(_OLD + "." + n) for n in _MEMO_NAMES] + [RetiredModule(f"tests.gates.{_OLD_DIR}.test_traced_memo_data"),
                                                                    RetiredModule(f"tests.gates.{_OLD_DIR}.test_traced_memo_manifest")]
    red = {row.module for (_k, row), found in _walk.walk(tmp_path, rows=rows).items() if found}
    assert red == {_OLD + "." + n for n in _MEMO_NAMES} | {f"tests.gates.{_OLD_DIR}.test_traced_memo_data"}


# ---------------------------------------------------------------- G3: the old path as text
def _path_expression_hits(source: str) -> list[int]:
    """Lines where a ``/`` chain or ``os.path.join`` call has the constant segment ``numerics`` IMMEDIATELY followed by a memo file name."""
    hits: list[int] = []

    def leaves(node: ast.expr) -> list[ast.expr]:
        return leaves(node.left) + leaves(node.right) if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div) else [node]

    def check(parts: list[ast.expr], line: int) -> None:
        consts = [p.value if isinstance(p, ast.Constant) and isinstance(p.value, str) else None for p in parts]
        if any(a == _OLD_DIR and b in tuple(f"{n}.py" for n in _MEMO_NAMES) for a, b in zip(consts, consts[1:])):
            hits.append(line)

    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            check(leaves(node), node.lineno)
        elif isinstance(node, ast.Call) and ast.unparse(node.func).endswith("path.join"):
            check(list(node.args), node.lineno)
    return sorted(set(hits))


def _text_findings(root: Path, files: list[str]) -> list[tuple[str, int, str]]:
    """``(file, line, pattern)`` over the text of ``files`` under ``root``; ``files`` is the walk's population."""
    if not files:
        raise AssertionError("input count 0: the walk's file list is empty (a vacuous zero)")
    found: list[tuple[str, int, str]] = []
    for rel in files:
        text = _text(root, rel)
        if text is None:
            continue
        for n, line in enumerate(text.splitlines(), 1):
            found += [(rel, n, name) for name, rx in _PATTERNS.items() if rx.search(line)]
        if _P4.search(text) and re.search(_M, text):
            found += [(rel, n, "P4") for n, line in enumerate(text.splitlines(), 1) if _P4.search(line)]
        if rel.endswith(".py"):
            found += [(rel, n, "PATHEXPR") for n in _path_expression_hits(text)]
    return found


def _code_files() -> list[str]:
    return [f for f in _tracked(*SCOPE)
            if f.endswith(".py") and not f.startswith("tests/gates/_retired/") and f != "tests/gates/test_memo_old_path.py"]


def _docs_findings(root: Path, files: list[str]) -> list[tuple[str, int, str]]:
    out = []
    for rel, n, pattern in _text_findings(root, files):
        line = (_text(root, rel) or "").splitlines()[n - 1]
        if _GIT_PAST.search(line) or (rel, line) in _PAST_ROWS:
            continue
        out.append((rel, n, pattern))
    return out


def test_no_tracked_code_file_spells_the_old_memo_path() -> None:
    files = _code_files()
    assert len(files) > 1000 and (REPO / "tools/cone_store/_recorder.py").is_file(), "the population is the tracked code of four trees"
    found = _text_findings(REPO, files)
    by_file = sorted({f for f, _n, _p in found})
    assert not found, f"{len(found)} findings in {len(by_file)} files (P1-P3, the test-location pattern, P4, the path expression): {by_file[:12]}"


def test_no_docs_page_spells_the_old_memo_path_outside_the_past_tense_table() -> None:
    files = [f for f in _tracked("docs") if not f.startswith("docs/_build/")]
    assert len(files) > 100 and (REPO / "docs/api/numerics.rst").is_file()
    found = _docs_findings(REPO, files)
    assert not found, f"{len(found)} docs lines in {len({f for f, _n, _p in found})} files name an old path and are in no tense row: {sorted({f for f, _n, _p in found})}"


def test_every_past_tense_row_names_a_live_line() -> None:
    dead = [(f, line) for f, line in _PAST_ROWS if line not in (_text(REPO, f) or "").splitlines()]
    assert not dead, f"dead exemption rows: {dead}"


def test_the_text_walk_finds_each_planted_spelling_and_reads_the_clean_sources_as_zero(tmp_path: Path) -> None:
    """Positive controls (run at base): each of P1, P2, P3, the test-location pattern, P4 and the path expression finds its planted
    witness; the new path, ``instance_cache`` and ``ORPHEUS_ROOT / "numerics"`` read 0."""
    planted = {
        "p1.py": f"import {_OLD}.traced_memo\n",
        "p2.py": f'x = "{_OLD_DIR}/_import_time.py"\n',
        "p3.py": f"from {_OLD} import space, traced_memo\n",
        "tp.py": f"import tests.gates.{_OLD_DIR}._traced_memo_api\n",
        "p4.py": f'a = ("orpheus", "{_OLD_DIR}")\nb = "_traced_memo_boot"\n',
        "pe.py": f'from pathlib import Path\nx = Path(".") / "{_OLD_DIR}" / "_traced_memo_boot.py"\n',
        "pj.py": f'import os\nx = os.path.join(t, "orpheus", "{_OLD_DIR}", "traced_memo.py")\n',
    }
    clean = {
        "c1.py": "import orpheus.memo.traced_memo\nfrom orpheus.memo import traced_memo\n",
        "c2.py": f"import {_OLD}.instance_cache\nfrom {_OLD} import traced_memo_x\n",
        "c3.py": f'ROOT = R / "{_OLD_DIR}"\ny = R / "{_OLD_DIR}" / "space.py"\n',
    }
    for name, text in {**planted, **clean}.items():
        (tmp_path / name).write_text(text)
    found = _text_findings(tmp_path, sorted(planted) + sorted(clean))
    assert {f for f, _n, _p in found} == set(planted), found
    assert {p for f, _n, p in found if f == "tp.py"} == {"TESTPATH"}
    assert {p for _f, _n, p in found} >= {"P1", "P2", "P3", "TESTPATH", "P4", "PATHEXPR"}
    with pytest.raises(AssertionError, match="input count 0"):
        _text_findings(tmp_path, [])


def test_a_git_path_at_a_past_commit_is_exempt_and_a_dead_exemption_is_not_taken(tmp_path: Path) -> None:
    (tmp_path / "d.rst").write_text(f"see ``git show 99897d10:orpheus/{_OLD_DIR}/traced_memo.py``\nlive: orpheus/{_OLD_DIR}/traced_memo.py\n")
    assert _docs_findings(tmp_path, ["d.rst"]) == [("d.rst", 2, "P2")]
