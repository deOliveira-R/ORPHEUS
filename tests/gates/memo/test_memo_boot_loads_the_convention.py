"""Gates G5, G6, G16, G17, G18 of unit ``memo-digest`` (specification part B 3.3, 4.2): how the boot loads the hash convention.

G5  the boot imports nothing outside the standard library before ``adopt_convention`` and ``_boot`` (two non-stdlib imports, in a
    stated order); ``values/digest.py`` is standard-library only. The predicate is run against a MODEL boot (a positive control that
    runs at base) and against the real file, with seven textual mutations of the model.
G6  the by-path module is THE module of its dotted name: one module object in a cold generating process and in an ordinary one.
G16 a recorded call of the convention is pinned as a def of ``values/digest.py``; the sealed bytecode entry obeys its law whole.
G17 the loader reads each pinned source once.
G18 the two ``__init__.py`` that ``import orpheus.values`` loads import nothing.

Imports of ``orpheus.memo`` and ``orpheus.values.digest`` are inside test bodies (one red per gate at base, not a collection error).
Every child is built by ``_memo_tree`` (store root asserted under the temporary tree, ``orpheus.__file__`` asserted under it).
A mutation arm raises ``_memo_tree.Uninstallable`` (a distinct error) when its anchor is not found exactly once: an arm that cannot
install is not a green arm.
"""
from __future__ import annotations

import ast
import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

from tests.gates.memo import _memo_tree
from tests.gates.memo._memo_tree import REPO, Uninstallable, build_tree, generating_child, mutate

pytestmark = pytest.mark.foundation

BOOT = REPO / "orpheus" / "memo" / "_traced_memo_boot.py"
DIGEST = REPO / "orpheus" / "values" / "digest.py"

# ================================================================ G5: the imports of the boot
_STDLIB = sys.stdlib_module_names


def _imports(tree: ast.AST) -> list[tuple[str | None, str, int]]:
    """``(enclosing function or None, module, line)`` of every import statement at any depth (``from . import x`` reads ``.``)."""
    found: list[tuple[str | None, str, int]] = []

    def visit(node: ast.AST, scope: str | None) -> None:
        for child in ast.iter_child_nodes(node):
            inner = child.name if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and scope is None else scope
            if isinstance(child, ast.Import):
                found.extend((scope, alias.name, child.lineno) for alias in child.names)
            elif isinstance(child, ast.ImportFrom):
                found.append((scope, "." * child.level + (child.module or ""), child.lineno))
            visit(child, inner)

    visit(tree, None)
    return found


def boot_findings(source: str) -> list[str]:
    """Violations of the boot's import law (empty when the law holds)."""
    tree = ast.parse(source)
    foreign = [(s, m, n) for s, m, n in _imports(tree) if m.split(".")[0] not in _STDLIB]
    out: list[str] = []
    expected = {("adopt_convention", "orpheus.values"), ("_boot", "orpheus.memo")}
    got = {(s, m) for s, m, _n in foreign}
    out += [f"unexpected non-stdlib import {m} in {s or 'module level'}" for s, m in sorted(got - expected, key=str)]
    out += [f"missing the import {m} in {s}" for s, m in sorted(expected - got)]
    if len(foreign) != len(got):
        out.append(f"{len(foreign)} non-stdlib imports for {len(got)} distinct (a second import of one module)")
    functions = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    boot, adopt = functions.get("_boot"), functions.get("adopt_convention")
    if boot is not None:
        path_line = next((n.lineno for n in ast.walk(boot) if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]).startswith("sys.path[")), None)
        call_line = next((n.lineno for n in ast.walk(boot) if isinstance(n, ast.Call) and ast.unparse(n.func) == "adopt_convention"), None)
        memo_line = next((n for s, m, n in foreign if s == "_boot" and m == "orpheus.memo"), None)
        if None in (path_line, call_line, memo_line):
            out.append(f"_boot lacks a step: sys.path {path_line}, adopt_convention() {call_line}, memo import {memo_line}")
        else:
            if not path_line < call_line < memo_line:
                out.append(f"_boot order: sys.path at {path_line}, adopt_convention() at {call_line}, memo import at {memo_line}")
            if any(s == "_boot" and n > memo_line for s, _m, n in foreign):
                out.append("a non-stdlib import follows the memo import in _boot")
    if adopt is not None:
        import_line = next((n for s, m, n in foreign if s == "adopt_convention"), None)
        assigns = [n.lineno for n in ast.walk(adopt) if isinstance(n, ast.Assign)]
        if import_line is not None and any(line < import_line for line in assigns):
            out.append("adopt_convention registers the module before `import orpheus.values` (decided D1: load, import, register)")
    return out


def convention_findings(source: str) -> list[str]:
    return [f"values/digest.py imports {m}" for _s, m, _n in _imports(ast.parse(source)) if m.split(".")[0] not in _STDLIB]


_MODEL = '''\
import importlib
import os
import sys
import types

if "orpheus" in sys.modules:
    _CONVENTION = _convention()
else:
    _CONVENTION = types.ModuleType("orpheus.values.digest")


def adopt_convention() -> None:
    import orpheus.values

    sys.modules["orpheus.values.digest"] = _CONVENTION
    orpheus.values.digest = _CONVENTION


def _boot():
    sys.modules["orpheus.memo._traced_memo_boot"] = sys.modules[__name__]
    sys.path[:] = ["x"]
    adopt_convention()
    from orpheus.memo import traced_memo

    return traced_memo
'''
_MUTATIONS = {
    "module-level import numpy": lambda s: "import numpy\n" + s,
    "module-level from orpheus.values import digest": lambda s: "from orpheus.values import digest\n" + s,
    "a second import orpheus.values inside _boot": lambda s: s.replace("    adopt_convention()\n", "    adopt_convention()\n    import orpheus.values\n"),
    "adopt_convention() before the sys.path assignment": lambda s: s.replace(
        '    sys.path[:] = ["x"]\n    adopt_convention()\n', '    adopt_convention()\n    sys.path[:] = ["x"]\n'),
    "memo import before adopt_convention()": lambda s: s.replace(
        "    adopt_convention()\n    from orpheus.memo import traced_memo\n", "    from orpheus.memo import traced_memo\n    adopt_convention()\n"),
    "import numpy in a helper": lambda s: s + "\n\ndef _helper():\n    import numpy\n",
    "registration above import orpheus.values": lambda s: s.replace(
        '    import orpheus.values\n\n    sys.modules["orpheus.values.digest"] = _CONVENTION\n',
        '    sys.modules["orpheus.values.digest"] = _CONVENTION\n    import orpheus.values\n\n'),
}


def test_the_import_law_holds_for_the_model_boot_and_each_of_seven_mutations_breaks_it() -> None:
    """Positive control of G5 (runs at base): the model satisfies the predicate; each mutation (anchored, changed text) is red."""
    assert boot_findings(_MODEL) == []
    for label, mutation in _MUTATIONS.items():
        mutated = mutation(_MODEL)
        assert mutated != _MODEL, f"Uninstallable: {label} did not change the model"
        assert boot_findings(mutated), f"the predicate is blind to: {label}"


def test_the_real_boot_imports_nothing_outside_the_standard_library_before_the_adoption() -> None:
    """G5. Red at base: the boot has ONE non-stdlib import (the memo, at its old path) and no ``adopt_convention``; a scratch probe of
    this predicate on the base file reads 1 where it demands 2 (commit body)."""
    assert BOOT.is_file(), f"{BOOT} does not exist (input count 0)"
    assert boot_findings(BOOT.read_text(encoding="utf-8")) == []
    assert DIGEST.is_file() and convention_findings(DIGEST.read_text(encoding="utf-8")) == []
    assert convention_findings("import hashlib\nimport numpy\nfrom orpheus import x\n") == ["values/digest.py imports numpy", "values/digest.py imports orpheus"]


# ================================================================ G18: the two __init__ import nothing
def _import_count(path: Path) -> int:
    return len(_imports(ast.parse(path.read_text(encoding="utf-8"))))


def test_the_two_init_files_that_import_orpheus_values_import_nothing(tmp_path: Path) -> None:
    """G18. Green at base (test-only first red: the planted ``import numpy`` in a copy). Positive control: the counter reads 1 on a
    planted ``import x``."""
    planted = tmp_path / "init.py"
    planted.write_text("import x\n")
    assert _import_count(planted) == 1
    assert (_import_count(REPO / "orpheus" / "__init__.py"), _import_count(REPO / "orpheus" / "values" / "__init__.py")) == (0, 0)


# ================================================================ G6: one module object of the convention
_GENERATING_BODY = '''
import gc, hashlib, importlib, types
import orpheus.values.digest as x
from orpheus.values.digest import digest as d
from orpheus.values import digest as y
import orpheus.values.content as c
dpath = os.path.realpath(os.path.join(R, "orpheus", "values", "digest.py"))
objs = [o for o in gc.get_objects() if isinstance(o, types.ModuleType) and getattr(o, "__file__", None) and os.path.realpath(o.__file__) == dpath]
facts.update(
    is_boot=x is module._CONVENTION, func=d is x.digest, attr=y is x, gc=len(objs), file=os.path.realpath(x.__file__) == dpath,
    content=c._convention is x, memo=traced_memo.hexdigest is x.hexdigest,
    pinned=module.PINNED_SOURCES[dpath] == module._CONVENTION_SOURCE,
    code=any(isinstance(k, types.CodeType) and k.co_name == "digest" and k == x.digest.__code__
             for k in compile(module._CONVENTION_SOURCE, dpath, "exec").co_consts),
    hexok=x.hexdigest(b"x") == hashlib.blake2b(b"x", digest_size=32).hexdigest())
'''
_ONE_OBJECT = ("is_boot", "func", "attr", "file", "content", "memo", "pinned", "code", "hexok")


def _good(facts: dict[str, object]) -> bool:
    return facts["gc"] == 1 and all(facts[k] is True for k in _ONE_OBJECT) and facts["orpheus_before"] is False


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_the_two_init_files_that_import_orpheus_values_import_nothing')
def test_a_generating_process_holds_one_module_object_of_the_convention_read_once(tmp_path: Path) -> None:
    """G6 (i). The facts of 3.1 (i) on the real ``_boot`` in a copy of the tree. Red at base: ``digest.py`` is absent."""
    tree = build_tree(tmp_path / "t")
    facts = generating_child(tree, _GENERATING_BODY, tmp_path / "w")
    assert _good(facts), facts


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_the_two_init_files_that_import_orpheus_values_import_nothing')
def test_an_ordinary_process_holds_one_module_object_and_no_source_bytes(tmp_path: Path) -> None:
    """G6 (ii): ``import orpheus.memo.traced_memo`` from the repository root; the name ``_CONVENTION_SOURCE`` exists in the by-path
    branch only (decided C5)."""
    program = (
        "import gc, json, os, sys, types\nimport orpheus.memo.traced_memo\nimport orpheus.values.digest as x\n"
        "boot = sys.modules['orpheus.memo._traced_memo_boot']\n"
        "dpath = os.path.realpath(x.__file__)\n"
        "n = sum(1 for o in gc.get_objects() if isinstance(o, types.ModuleType) and getattr(o, '__file__', None) and os.path.realpath(o.__file__) == dpath)\n"
        "print(json.dumps({'gc': n, 'same': boot._CONVENTION is x, 'nosource': not hasattr(boot, '_CONVENTION_SOURCE')}))\n")
    import json

    assert json.loads(_memo_tree.cold_child(REPO, program).splitlines()[-1]) == {"gc": 1, "same": True, "nosource": True}


def test_the_by_path_load_with_no_boot_run_gives_the_same_file_and_equal_outputs() -> None:
    """G6 (iii): the cone-store capture's mode (decision 9.Q8 (B), #734): only what holds under either ruling."""
    program = (
        "import importlib.util, json, os, sys\n"
        f"R = {str(REPO)!r}\n"
        "p = os.path.join(R, 'orpheus', 'memo', '_traced_memo_boot.py')\n"
        "s = importlib.util.spec_from_file_location('orpheus.memo._traced_memo_boot', p)\n"
        "m = importlib.util.module_from_spec(s); sys.modules[s.name] = m; s.loader.exec_module(m)\n"
        "if 'orpheus' in sys.modules: raise SystemExit('orpheus imported by the boot')\n"
        "sys.path.insert(0, R)\n"
        "import orpheus.values.digest as x\n"
        "print(json.dumps([os.path.realpath(m._CONVENTION.__file__) == os.path.realpath(x.__file__), m._CONVENTION.hexdigest(b'q') == x.hexdigest(b'q')]))\n")
    import json

    done = subprocess.run([sys.executable, "-P", "-c", program], capture_output=True, text=True, cwd=REPO.parent)
    assert done.returncode == 0, done.stderr
    assert json.loads(done.stdout.splitlines()[-1]) == [True, True]


@pytest.mark.parametrize(("label", "anchor", "replacement"), [
    ("no adopt_convention() call", "    adopt_convention()\n", ""),
    ("no sys.modules registration (the attribute assignment stays)", '    sys.modules["orpheus.values.digest"] = _CONVENTION\n', ""),
])
@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_a_generating_process_holds_one_module_object_of_the_convention_read_once')
def test_g6_mutation_a_boot_that_does_not_adopt_leaves_two_module_objects(tmp_path: Path, label: str, anchor: str, replacement: str) -> None:
    """G6 mutations (1) and (2): the honest tree reads good (computed first); the mutant's facts are not good. [R] at base: the
    anchors do not exist until sub-unit D, so the arm raises ``Uninstallable`` (a distinct error), never a quiet green."""
    honest = generating_child(build_tree(tmp_path / "ok"), _GENERATING_BODY, tmp_path / "w1")
    assert _good(honest), honest
    mutant = generating_child(build_tree(tmp_path / "bad", memo_edit=("_traced_memo_boot.py", anchor, replacement)), _GENERATING_BODY, tmp_path / "w2")
    assert not _good(mutant), (label, mutant)


# ================================================================ G6 (iv), (v): the by-path module has a spec and is not served stale bytecode
_SPEC_BODY = '''
import importlib.util
import orpheus.values.digest as x
found = importlib.util.find_spec("orpheus.values.digest")
facts.update(spec_is_own=found is x.__spec__, spec_loader=repr(x.__spec__.loader),
             named=os.path.realpath(x.__file__) in traced_memo._module_names())
'''


def test_a_generating_process_answers_find_spec_for_the_convention_and_names_its_file(tmp_path: Path) -> None:
    """G6 (iv). The by-path module carries a spec, so ``find_spec`` returns the module's own ``__spec__`` (it raised ``ValueError``
    on a module with no ``__spec__``) and the memo's ``_module_names`` keeps ``digest.py``; the spec has NO loader, so nothing can
    re-read the file through it. Red before the fix: ``find_spec`` raised and ``named`` was False."""
    facts = generating_child(build_tree(tmp_path / "t"), _SPEC_BODY, tmp_path / "w")
    assert facts["spec_is_own"] is True and facts["named"] is True and facts["spec_loader"] == "None", facts


_STALE_BODY = '''
import re
x_text = open(os.path.join(R, "orpheus", "values", "digest.py"), encoding="utf-8").read()
import orpheus.values.digest as x
facts.update(running=x.NAME_SIZE, disk=int(re.search(r"NAME_SIZE: Final = (\\d+)", x_text).group(1)))
'''
_CONTROL = ("import orpheus.values.digest as x\nprint(x.NAME_SIZE)\n")


def _tree_with_a_stale_pyc(dest: Path, *, memo_edit: tuple[str, str, str] | None = None) -> Path:
    """A copy of the tree whose ``digest.py`` has a timestamp ``.pyc`` (both optimisation levels) of ``NAME_SIZE = 8`` while the file
    says ``9``: the same-size edit, the file's mtime restored, so the interpreter finds the ``.pyc`` valid."""
    import os
    import py_compile

    tree = build_tree(dest, memo_edit=memo_edit)
    path = tree / "orpheus" / "values" / "digest.py"
    before = path.stat()
    for level in (0, 1):
        py_compile.compile(str(path), doraise=True, optimize=level, invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP)
    mutate(path, "NAME_SIZE: Final = 8", "NAME_SIZE: Final = 9")
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert path.stat().st_size == before.st_size
    return tree


def _stale_readings(tree: Path, tmp_path: Path) -> tuple[str, dict[str, object]]:
    control = _memo_tree.cold_child(tree, _CONTROL).strip()
    return control, generating_child(tree, _STALE_BODY, tmp_path / "w")


def test_a_generating_process_computes_with_the_disk_value_when_a_stale_pyc_is_live(tmp_path: Path) -> None:
    """G6 (v). The positive control: an ordinary import in the same tree reads the stale ``.pyc`` (``8``) while the disk says ``9``;
    a generating process compiles the bytes it pins and so runs ``9``. Green before and after the spec fix, by design: its teeth
    are the next test's mutation (an ordinary import in place of the by-path load)."""
    control, facts = _stale_readings(_tree_with_a_stale_pyc(tmp_path / "t"), tmp_path)
    assert control == "8", control  # the stale .pyc is live
    assert (facts["running"], facts["disk"]) == (9, 9), facts


def test_g6_mutation_an_ordinary_import_in_place_of_the_by_path_load_runs_the_stale_pyc(tmp_path: Path) -> None:
    """G6 (v) mutation: replace the by-path load by ``importlib.import_module`` (with the tree on ``sys.path``): the generating
    process then runs the ``.pyc``'s ``8`` against a disk of ``9``, and the honest leg above reads red."""
    text = BOOT.read_text(encoding="utf-8")
    start = text.index("    with open(_CONVENTION_PATH")
    stop = text.index("def adopt_convention")
    anchor = text[start:stop].rstrip("\n") + "\n"
    replacement = (
        "    sys.path.insert(0, os.path.dirname(_PACKAGE))\n    _CONVENTION = importlib.import_module(\"orpheus.values.digest\")\n"
        "    sys.path.remove(os.path.dirname(_PACKAGE))\n    _CONVENTION_SOURCE = open(_CONVENTION_PATH, \"rb\").read()\n")
    tree = _tree_with_a_stale_pyc(tmp_path / "t", memo_edit=("_traced_memo_boot.py", anchor, replacement))
    _, facts = _stale_readings(tree, tmp_path)
    assert (facts["running"], facts["disk"]) == (8, 9), facts


# ================================================================ G16: a recorded call is pinned; the sealed entry's law
def _api():
    from tests.gates.memo import _traced_memo_api

    return _traced_memo_api


def _boot_module():
    from orpheus.memo import _traced_memo_boot

    return _traced_memo_boot


def _convention():
    from orpheus.values import digest

    return digest


def test_a_recorded_call_of_the_convention_is_pinned_as_a_def_of_values_digest() -> None:
    """G16 leg (1). The recorder is in-process (``trace_call``); the store is read, never written (asserted under the worktree)."""
    boot = _boot_module()
    if not Path(boot.CACHE).resolve().is_relative_to(REPO.resolve()):
        raise Uninstallable(f"the store root {boot.CACHE} is outside the worktree")
    convention = _convention()
    _value, manifest = _api().trace_call(convention.hexdigest, b"x")
    names = {(p.path, p.qualname) for p in manifest.functions}
    assert any(path.endswith("orpheus/values/digest.py") and qualname == "hexdigest" for path, qualname in names), sorted(names)[:5]


def seal_laws(boot, key: bytes, payload: bytes) -> dict[str, bool]:
    entry = boot.seal_entry(key, payload)
    return {
        "roundtrip": boot.unseal_entry(key, entry) == payload,
        "key": entry[:32] == key,
        "digest": entry[32:64] == hashlib.blake2b(payload, digest_size=32).digest(),
        "payload": entry[64:] == payload,
    }


def test_the_sealed_bytecode_entry_round_trips_and_nothing_else_unseals() -> None:
    """G16 legs (2) and (3), the law whole: round trip, its converse on five bad entries, and the entry's bytes against ``hashlib``."""
    boot = _boot_module()
    key = hashlib.blake2b(b"key", digest_size=32).digest()
    other = hashlib.blake2b(b"other", digest_size=32).digest()
    for payload in (b"", b"p", bytes(range(256)) * 40):
        assert all(seal_laws(boot, key, payload).values()), payload[:4]
        entry = boot.seal_entry(key, payload)
        flipped = entry[:40] + bytes([entry[40] ^ 1]) + entry[41:]
        for bad in (flipped, entry[:-1] if payload else entry[:63], b"", entry + b"x"):
            assert boot.unseal_entry(key, bad) is None
        assert boot.unseal_entry(other, entry) is None
        assert boot.unseal_entry(key, entry[:31]) is None


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_the_sealed_bytecode_entry_round_trips_and_nothing_else_unseals')
def test_g16_mutation_a_self_consistent_wrong_digest_passes_the_round_trip_and_fails_the_content_leg(monkeypatch: pytest.MonkeyPatch) -> None:
    """The pair stays self-consistent under a constant digest; only leg (3) sees it (mutate INSIDE the law's class)."""
    boot = _boot_module()
    key = hashlib.blake2b(b"key", digest_size=32).digest()
    assert all(seal_laws(boot, key, b"payload").values())
    monkeypatch.setattr(boot._CONVENTION, "digest", lambda data: b"\0" * 32)
    laws = seal_laws(boot, key, b"payload")
    assert laws["roundtrip"] and not laws["digest"], laws


# ================================================================ G17: the loader reads each pinned source once
def test_the_loader_compiles_the_pinned_bytes_and_keeps_the_pin(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """G17. Red at base: the loader overwrites the pin and compiles the disk's bytes (``f()`` is 2). The store of this test is the
    tmp directory (``_BYTECODE`` is monkeypatched; the live bytecode cache is never written)."""
    boot = _boot_module()
    source = tmp_path / "m.py"
    source.write_bytes(b"def f(): return 2\n")
    real = str(source.resolve())
    monkeypatch.setattr(boot, "PINNED_SOURCES", {real: b"def f(): return 1\n"})
    monkeypatch.setattr(boot, "_BYTECODE", str(tmp_path / "bytecode"))
    namespace: dict[str, object] = {}
    exec(boot._PinnedSourceLoader("m", str(source)).get_code("m"), namespace)
    assert namespace["f"]() == 1  # type: ignore[operator]
    assert boot.PINNED_SOURCES[real] == b"def f(): return 1\n"
    monkeypatch.setattr(boot, "PINNED_SOURCES", {})
    exec(boot._PinnedSourceLoader("m", str(source)).get_code("m"), namespace)
    assert namespace["f"]() == 2 and boot.PINNED_SOURCES[real] == b"def f(): return 2\n"  # type: ignore[operator]


@pytest.mark.rests_on('tests/gates/memo/test_memo_boot_loads_the_convention.py::test_the_loader_compiles_the_pinned_bytes_and_keeps_the_pin')
def test_g17_mutation_a_loader_that_reads_again_compiles_the_disk_and_overwrites_the_pin(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """The arm: a COPY of the boot with ``real in PINNED_SOURCES`` replaced by ``False``, loaded by path under a private name. [R] at
    base (the anchor appears with the build): ``Uninstallable``."""
    import importlib.util

    copy = tmp_path / "orpheus" / "memo" / "_traced_memo_boot.py"
    copy.parent.mkdir(parents=True)
    copy.write_text(BOOT.read_text(encoding="utf-8"), encoding="utf-8")
    mutate(copy, "real in PINNED_SOURCES", "False")
    spec = importlib.util.spec_from_file_location("_g17_mutant_boot", copy)
    assert spec is not None and spec.loader is not None
    mutant = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mutant)
    source = tmp_path / "m.py"
    source.write_bytes(b"def f(): return 2\n")
    real = str(source.resolve())
    monkeypatch.setattr(mutant, "PINNED_SOURCES", {real: b"def f(): return 1\n"})
    monkeypatch.setattr(mutant, "_BYTECODE", str(tmp_path / "bytecode"))
    namespace: dict[str, object] = {}
    exec(mutant._PinnedSourceLoader("m", str(source)).get_code("m"), namespace)
    assert namespace["f"]() == 2 and mutant.PINNED_SOURCES[real] == b"def f(): return 2\n"  # type: ignore[operator]
