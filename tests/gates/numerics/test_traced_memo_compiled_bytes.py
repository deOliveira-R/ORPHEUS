r"""#602, the two defects found by round two's specification: a generation runs the bytes it pins.

Spec ``scratch/characteristic_architecture/p3/m602/round2_gate_spec.md`` §0, rows P1–P9 and A1–A2.

**P (the stale timestamp ``.pyc``).** CPython accepts a timestamp ``.pyc`` when the source's whole-second mtime and
its size match. A same-size edit within that second therefore ran the OLD code in a generating process, while the
manifest pinned the NEW bytes. The regenerated entry served the old answer and validated (the probe
``ta/p8_pyc_through_memo.py``: 4.5 where 6.0 is right, then a ``Hit``). The fix: the boot script installs a
loader that compiles every first-party module from the bytes it reads, records those bytes for the manifest, and
keeps its bytecode in a content-keyed cache apart from ``__pycache__``.

**A (a docstring-only def).** The docstring strip put a ``pass`` with no location in a body that was only a
docstring, and ``_def_pin`` read that body's first line. An ``__annotate__`` of such a def (a ``Protocol`` stub
whose annotations are read) raised ``AttributeError`` in every generation that ran it.

Every generator is written under ``tmp_path``, and every edit is to that copy. The bytecode cache is the
production one (``.cache/bytecode``, gitignored): it is content-keyed, so a test's entries are never another
module's.
"""
from __future__ import annotations

import ast
import hashlib
import importlib
import marshal
import os
import subprocess
import sys
import textwrap
import threading
import time
from pathlib import Path

import pytest

from orpheus.numerics import _traced_memo_boot

from . import _traced_memo_api as api
from .test_traced_memo_findings import _Package

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]
_M31 = 'tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never'
_M32 = 'tests/gates/numerics/test_traced_memo_process.py::test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does'


@pytest.fixture
def make(tmp_path):
    made: list[_Package] = []

    def build(**modules: str) -> _Package:
        made.append(_Package(tmp_path, **modules))
        return made[-1]

    with api.cache_root(tmp_path / "cache"):
        yield build
    for package in made:
        package.close()


def _boot():
    return _traced_memo_boot


def _bytecode_dir(optimize: int = sys.flags.optimize) -> Path:
    """The bytecode cache of one ``-O`` level, spelled from its documented layout (not read off the boot)."""
    return REPO / ".cache" / "bytecode" / f"{sys.implementation.cache_tag}-O{optimize}"


def _bytecode_entry(source: Path, optimize: int = sys.flags.optimize) -> Path:
    """The cache entry of ``source`` as a generating process names it: digest(path, interpreter, ``-O``, bytes),
    spelled here independently. The loader's path is the import's (a ``sys.path`` entry joined with the module's
    relative path), which may or may not be the resolved one; the existing candidate is returned."""
    data = source.read_bytes()
    candidates = []
    for path in dict.fromkeys((str(source), os.path.realpath(source))):
        head = f"{path}\0{sys.implementation.cache_tag}\0{optimize}\0".encode()
        candidates.append(_bytecode_dir(optimize) / hashlib.blake2b(head + data, digest_size=32).hexdigest())
    found = [c for c in candidates if c.is_file()]
    assert found, f"no bytecode entry for {source} among {candidates}"
    return found[0]


def _entry_is_whole(data: bytes, key: bytes | None = None) -> bool:
    """An entry is its key (32 bytes), the digest of its payload (32 bytes) and the marshalled payload."""
    return (len(data) > 64 and (key is None or data[:32] == key)
            and data[32:64] == hashlib.blake2b(data[64:], digest_size=32).digest())


GENERATOR = '''
    from orpheus.numerics.traced_memo import traced_memo

    SCALE = 3.0


    @traced_memo
    def generate(x: float) -> float:
        return x * SCALE
'''


@pytest.mark.rests_on(_M31, _M32)
@pytest.mark.catches("ERR-115")
def test_p1_a_same_size_edit_in_the_same_second_regenerates_the_new_answer(make):
    """P1 (the probe ``ta/p8`` as a gate): a same-size edit whose mtime is put back in the same nanosecond as the
    previous source, the state a scripted edit within one second produces and a timestamp ``.pyc`` cannot tell
    apart, is ``Stale``. The regeneration returns the NEW answer, and the entry it writes then validates.
    First red: the loader not installed, where the child compiled the stale ``.pyc`` and returned 4.5."""
    package = make(m=GENERATOR)
    m = package.module()  # the parent's import writes the timestamp .pyc of the old source
    assert m.generate(1.5) == 4.5
    source = package.dir / "m.py"
    before = os.stat(source)
    package.edit("SCALE = 3.0", "SCALE = 4.0")
    os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert os.stat(source).st_size == before.st_size, "the edit must keep the size"
    assert api.verdict_kind(m.generate.lookup(1.5)) == "Stale"
    assert m.generate(1.5) == 6.0, "the generating process ran the old code (a stale .pyc)"
    assert api.verdict_kind(m.generate.lookup(1.5)) == "Hit"


@pytest.mark.timeout(300)
@pytest.mark.rests_on(_M31, _M32)
@pytest.mark.catches("ERR-115")
def test_p2_an_edit_landing_during_a_generation_leaves_the_entry_stale(make, tmp_path):
    """P2: the generator waits on a barrier after its module was compiled. The test then edits the module and
    releases it, so the generation runs the old code. The entry must pin the bytes that were COMPILED, so it is
    ``Stale`` against the edited file and never served. First red: bytes re-read after the run, which pin the
    new bytes beside the old answer and make the entry a ``Hit``."""
    ready, go = tmp_path / "ready", tmp_path / "go"
    package = make(m=f'''
        import os
        import time

        from orpheus.numerics.traced_memo import traced_memo

        SCALE = 3.0


        @traced_memo
        def slow(x: float) -> float:
            open({str(ready)!r}, "w").close()
            while not os.path.exists({str(go)!r}):
                time.sleep(0.01)
            return x * SCALE
    ''')
    m = package.module()
    outcome: list[object] = []
    thread = threading.Thread(target=lambda: outcome.append(_call(m.slow, 1.5)))
    thread.start()
    deadline = time.monotonic() + 120
    while not ready.exists():
        assert time.monotonic() < deadline, "the generation never reached its barrier"
        time.sleep(0.01)
    package.edit("SCALE = 3.0", "SCALE = 4.0")
    go.touch()
    thread.join(timeout=120)
    assert not thread.is_alive()
    (result,) = outcome
    # The generation ran the old code; what the caller receives is either that answer, or the memo's refusal to
    # serve an entry that no longer validates. A Hit of the old answer is the defect.
    assert result == 4.5 or (isinstance(result, RuntimeError) and "does not validate" in str(result)), result
    assert api.verdict_kind(m.slow.lookup(1.5)) == "Stale"
    assert m.slow(1.5) == 6.0


def _call(function, *args):
    try:
        return function(*args)
    except Exception as error:  # noqa: BLE001 - the outcome is asserted by the gate
        return error


@pytest.mark.rests_on(_M31)
def test_p3_the_loaders_own_reads_are_not_the_runs(make, tmp_path):
    """P3: the loader reads every source file and its bytecode-cache entry, and none of those reads is the
    run's dependency. No ``DataPin`` names a module the run imported, and no data or presence pin lies under
    the bytecode cache. First red: the loader's body run outside ``unrecorded()``, which pins every cache
    entry as data."""
    package = make(m=GENERATOR)
    package.module().generate(1.5)
    (entry,) = api.entry_files(tmp_path / "cache")
    record = api.read_entry(entry)
    cache = os.path.realpath(REPO / ".cache" / "bytecode")
    data = [path for path, _ in api.manifest_rows(record, "DataPin")]
    probed = [path for path, _ in api.manifest_rows(record, "PresencePin")]
    assert not [p for p in data if p.endswith(".py")], data
    assert not [p for p in data + probed if os.path.realpath(p).startswith(cache)], (data, probed)
    assert api.manifest_rows(record, "ModulePin"), "the instrument read no module pins: the entry is not a generation's"


LOOSE_INIT = '''
    import importlib.util
    import sys

    LOOSE = {loose!r}


    class _LooseFinder:
        """A meta-path finder, appended as an editable install's is: it alone finds ``{name}``."""

        @classmethod
        def find_spec(cls, fullname, path=None, target=None):
            if fullname != {name!r}:
                return None
            return importlib.util.spec_from_file_location(fullname, LOOSE)


    sys.meta_path.append(_LooseFinder)
'''


@pytest.mark.rests_on(_M31, _M32)
def test_p4_a_module_found_by_a_meta_path_finder_is_compiled_from_pinned_bytes(make, tmp_path):
    """P4: a module that only a meta-path finder appended after the path finder knows (as an editable install's
    finder is) is compiled by the pinned loader in the generating process. It is pinned by its absolute path,
    an edit to it is ``Stale``, and the regeneration returns the new answer. First red: a pinned-source finder
    that sees only path-finder imports, which leaves the module with no recorded bytes and refuses the run."""
    loose = tmp_path / "loose" / "loose_p4.py"
    loose.parent.mkdir()
    loose.write_text("def f(x):\n    return x + 1.0\n")
    name = f"loose_p4_{os.urandom(4).hex()}"
    package = make(__init__=LOOSE_INIT.format(loose=str(loose), name=name), m=f'''
        from orpheus.numerics.traced_memo import traced_memo


        @traced_memo
        def via_loose(x: float) -> float:
            import {name}
            return {name}.f(x)
    ''')
    m = package.module()
    assert m.via_loose(1.0) == 2.0
    (entry,) = api.entry_files(tmp_path / "cache")
    defs = {(path, qualname) for path, qualname, *_ in api.manifest_rows(api.read_entry(entry), "DefPin")}
    assert (str(loose.resolve()), "f") in defs, defs
    loose.write_text("def f(x):\n    return x + 5.0\n")
    assert api.verdict_kind(m.via_loose.lookup(1.0)) == "Stale"
    assert m.via_loose(1.0) == 6.0


_P4B_DRIVER = '''
import sys
sys.path.insert(0, {root!r})
import orpheus  # through the editable install's finder: the repository root is on no sys.path entry
assert not any(p and orpheus.__file__.startswith(p.rstrip("/") + "/") for p in sys.path), (orpheus.__file__, sys.path)
from orpheus.numerics import traced_memo as tm
import {name}.m as m
with tm.cache_root({cache!r}):
    print(m.generate(1.5))
'''


@pytest.mark.rests_on(_M31)
def test_p4b_a_process_that_imports_orpheus_through_the_editable_finder_generates(make, tmp_path):
    """P4b: a script whose ``sys.path`` does not hold the repository (``orpheus`` found by the editable install's
    finder, which ``site`` imports before the boot script runs) generates and returns the right answer. In its
    generating process the finder module and the boot script were compiled before the loader existed, and
    their bytes are recorded when it is installed. First red: the first cut of the fix, which refused the
    editable finder module (no bytes recorded for it).

    Scope: the script imports the MAIN tree's ``orpheus`` (the editable finder points there), so a battery arm
    applied to a copy of the tree cannot reach this row."""
    package = make(m=GENERATOR)
    script = tmp_path / "outside" / "driver.py"
    script.parent.mkdir()
    script.write_text(_P4B_DRIVER.format(root=str(package.root), name=package.name, cache=str(tmp_path / "cache_p4b")))
    run = subprocess.run([sys.executable, *(["-O"] if sys.flags.optimize else []), str(script)], cwd=script.parent,
                         capture_output=True, text=True, env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"})
    assert run.returncode == 0, run.stderr[-3000:]
    assert run.stdout.split()[-1] == "4.5", run.stdout


@pytest.mark.rests_on(_M31)
def test_p5_the_generation_writes_no_pycache_and_its_bytecode_entry_carries_its_key(make):
    """P5: a module the generation imports first (``n``, never imported by the parent) leaves no ``.pyc`` in
    ``__pycache__``. Its compiled code is in the bytecode cache under its content key, as the key (32 bytes),
    the digest of the payload (32 bytes) and the marshalled code. First red: the loader not installed, which
    writes a timestamp ``.pyc`` for ``n``."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo


        @traced_memo
        def lazily(x: float) -> float:
            from . import n
            return n.SHIFT + x
    ''', n="SHIFT = 0.25\n")
    package.module().lazily(1.0)
    pycache = package.dir / "__pycache__"
    assert not [p.name for p in pycache.glob("n.*")] if pycache.is_dir() else True
    data = _bytecode_entry(package.dir / "n.py").read_bytes()
    assert _entry_is_whole(data), "the entry does not carry the digest of its payload"
    assert data[:32].hex() == _bytecode_entry(package.dir / "n.py").name, "the entry does not carry its key"
    assert isinstance(marshal.loads(data[64:]), type(compile("", "", "exec")))


_P6_DRIVER = '''
import sys
sys.path[:0] = [{repo!r}, {root!r}]  # the tree under test first, so a copy's boot is the one generating
from orpheus.numerics import traced_memo as tm
import {name}.m as m
with tm.cache_root({cache!r}):
    print(m.debug_level(0.0))
'''


@pytest.mark.rests_on(_M31)
def test_p6_each_optimisation_level_has_its_own_bytecode(make, tmp_path):
    """P6: ``__debug__`` is a constant the compiler folds, so the same bytes compile to different code under
    ``-O`` and without it. A generation at the OTHER level (a subprocess) caches its code first. A generation at
    this level then still answers with this level's value, and the two levels' entries lie in two directories.
    First red: a cache key and directory without the ``-O`` level, which serves the other level's code."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo


        @traced_memo
        def debug_level(x: float) -> float:
            return x + (1.0 if __debug__ else 2.0)
    ''')
    other = 0 if sys.flags.optimize else 1
    script = tmp_path / "driver.py"
    script.write_text(_P6_DRIVER.format(repo=str(REPO), root=str(package.root), name=package.name, cache=str(tmp_path / "cache_other")))
    run = subprocess.run([sys.executable, *(["-O"] if other else []), str(script)], cwd=REPO, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-3000:]
    assert float(run.stdout.split()[-1]) == (1.0 if other == 0 else 2.0)
    assert package.module().debug_level(0.0) == (2.0 if sys.flags.optimize else 1.0)
    source = package.dir / "m.py"
    assert _bytecode_entry(source, other).parent != _bytecode_entry(source).parent


_P7_SOURCE = '''
    from orpheus.numerics.traced_memo import traced_memo

    SIBLING = {sibling!r}


    @traced_memo
    def run_sibling(x: float) -> float:
        namespace = {{}}
        with open(SIBLING) as handle:
            exec(compile(handle.read(), SIBLING, "exec"), namespace)
        return namespace["f"](x)
'''


@pytest.mark.rests_on(_M31)
def test_p7_code_of_a_file_the_loader_did_not_compile_is_refused(make, tmp_path):
    """P7: a generator that ``exec``\\ s a sibling source file runs code from bytes no loader recorded. A re-read
    after the run could pin bytes that never ran, so the run is refused, by a message naming exec and runpy,
    and nothing is written. First red: today's re-read, which pins the sibling's bytes after the run."""
    package = make(m=_P7_SOURCE.format(sibling="SIBLING_PLACEHOLDER"), sib="def f(x):\n    return x + 3.0\n")
    package.edit("SIBLING_PLACEHOLDER", str(package.dir / "sib.py"))
    m = package.module()
    with pytest.raises(api.name("Unpinnable"), match="exec or runpy"):
        m.run_sibling(1.0)
    assert api.entry_files(tmp_path / "cache") == []


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_compiled_bytes.py::test_p7_code_of_a_file_the_loader_did_not_compile_is_refused')
def test_p7b_trace_call_still_pins_the_exec_d_file_by_a_re_read(make):
    """P7b: :func:`trace_call`, the in-process instrument, has no loader; it keeps the re-read, so the same
    generator's sibling def is pinned. First red: the refusal applied to ``trace_call`` too."""
    package = make(m=_P7_SOURCE.format(sibling="SIBLING_PLACEHOLDER"), sib="def f(x):\n    return x + 3.0\n")
    package.edit("SIBLING_PLACEHOLDER", str(package.dir / "sib.py"))
    value, manifest = api.trace_call(package.module().run_sibling.__wrapped__, 1.0)
    assert value == 4.0
    assert any(path.endswith("sib.py") and qualname == "f" for path, qualname, *_ in manifest.functions), manifest.functions


@pytest.mark.timeout(300)
def test_p8_concurrent_writers_never_leave_a_partial_bytecode_entry(tmp_path):
    """P8: two processes each compile one module 200 times through the pinned loader, deleting the entry before
    each compile, while this process reads the entry as fast as it can. Every read that finds the file finds a
    whole entry: its key, the digest of its payload, and the payload. The loader recompiles a partial entry
    silently, so the answer is blind to a write in place; the file itself is the instrument. First red: a write
    straight to the entry's name, with no ``os.replace``."""
    source = tmp_path / f"p8_{os.urandom(4).hex()}.py"
    source.write_text("VALUE = " + " + ".join(str(i) for i in range(400)) + "\n")
    boot = _boot()
    key = boot._bytecode_key(str(source), source.read_bytes())
    entry = Path(boot._BYTECODE) / key.hex()
    writer = textwrap.dedent(f'''
        import importlib.util, os
        spec = importlib.util.spec_from_file_location("boot", {boot.__file__!r})
        boot = importlib.util.module_from_spec(spec); spec.loader.exec_module(boot)
        for _ in range(200):
            try:
                os.remove({str(entry)!r})
            except FileNotFoundError:
                pass
            boot._PinnedSourceLoader("p8", {str(source)!r}).get_code("p8")
        print("done")
    ''')
    procs = [subprocess.Popen([sys.executable, *(["-O"] if sys.flags.optimize else []), "-c", writer], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(2)]
    reads = partial = 0
    while any(p.poll() is None for p in procs):
        try:
            data = entry.read_bytes()
        except FileNotFoundError:
            continue
        reads += 1
        partial += not _entry_is_whole(data, key)
    for p in procs:
        out, err = p.communicate()
        assert p.returncode == 0 and "done" in out, err[-2000:]
    assert reads > 0, "the reader never saw the entry: the instrument read nothing"
    assert partial == 0, f"{partial} of {reads} reads found a partial bytecode entry"


def _corruptions(data: bytes, other: bytes) -> dict[str, bytes]:
    flip = lambda b, i: b[:i] + bytes([b[i] ^ 0x01]) + b[i + 1:]  # noqa: E731
    return {"empty": b"", "half": data[: len(data) // 2], "payload-bit": flip(data, len(data) - 5),
            "digest-bit": flip(data, 40), "other-file": other}


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_compiled_bytes.py::test_p5_the_generation_writes_no_pycache_and_its_bytecode_entry_carries_its_key')
@pytest.mark.parametrize("kind", ["empty", "half", "payload-bit", "digest-bit", "other-file"])
def test_p9_a_corrupt_bytecode_entry_is_recompiled_never_run(make, kind):
    """P9: the cache entry of the generator's module is corrupted (truncated to nothing or to half, one bit of
    its code flipped, one bit of its stored digest flipped, or replaced by another source's whole entry). The
    next generation answers from the module's own bytes, and the entry is whole again afterwards. First red:
    the key and digest checks skipped, so the other source's code runs (``other-file``)."""
    package = make(m=GENERATOR)
    m = package.module()
    assert m.generate(1.0) == 3.0
    path = _bytecode_entry(package.dir / "m.py")
    code = compile("SCALE = 99.0\nfrom orpheus.numerics.traced_memo import traced_memo\n"
                   "@traced_memo\ndef generate(x: float) -> float:\n    return x * SCALE\n", str(package.dir / "m.py"), "exec")
    payload = marshal.dumps(code)
    other = hashlib.blake2b(b"another source", digest_size=32).digest() + hashlib.blake2b(payload, digest_size=32).digest() + payload
    path.write_bytes(_corruptions(path.read_bytes(), other)[kind])
    assert m.generate(2.0) == 6.0, f"a corrupt bytecode entry ({kind}) ran"
    assert _entry_is_whole(path.read_bytes(), bytes.fromhex(path.name)), "the corrupt entry was not replaced"


@pytest.mark.rests_on(_M31)
def test_a1_a_protocol_stub_whose_annotations_are_read_generates(make, tmp_path):
    """A1: a generator reads the annotations of a ``Protocol`` method whose body is only a docstring, which runs
    its ``__annotate__``. The generation succeeds. Its manifest pins the module's skeleton (which holds the
    signature) and no body of ``Proto.method``, since that body never ran. First red: the strip's ``pass`` with
    no location (``AttributeError: 'Pass' object has no attribute 'lineno'`` crossed back)."""
    package = make(m='''
        import typing

        from orpheus.numerics.traced_memo import traced_memo


        class Proto(typing.Protocol):
            def method(self, x: int) -> int:
                """Only a docstring."""


        @traced_memo
        def hints(x: float) -> float:
            return x + len(typing.get_type_hints(Proto.method))
    ''')
    assert package.module().hints(1.0) == 3.0
    (entry,) = api.entry_files(tmp_path / "cache")
    record = api.read_entry(entry)
    assert any(path.endswith("m.py") for path, _ in api.manifest_rows(record, "ModulePin"))
    assert not [q for _, q, *_ in api.manifest_rows(record, "DefPin") if q == "Proto.method"]


def _tracked_python() -> list[Path]:
    listed = subprocess.run(["git", "ls-files", "*.py"], cwd=REPO, capture_output=True, text=True, check=True).stdout.split()
    return [REPO / p for p in listed if (REPO / p).is_file()]


def _docstring_only_bodies(tree: ast.AST) -> list[ast.Expr]:
    """The docstring of every def or class whose body is ONLY that docstring."""
    found = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and len(node.body) == 1:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                found.append(first)
    return found


def _with_pass_for_docstring_only_bodies(data: bytes, docstrings: list[ast.Expr]) -> bytes:
    """``data`` with each docstring-only body's docstring replaced, in the TEXT, by ``pass`` at its indentation."""
    lines = data.decode().splitlines(keepends=True)
    for doc in sorted(docstrings, key=lambda d: d.lineno, reverse=True):
        indent = lines[doc.lineno - 1][: doc.col_offset]
        if indent.strip():  # the docstring shares its line with the header (``def f(): """d"""``)
            line = lines[doc.lineno - 1]
            lines[doc.lineno - 1 : doc.end_lineno] = [line[: doc.col_offset] + "pass\n"]
        else:
            lines[doc.lineno - 1 : doc.end_lineno] = [indent + "pass\n"]
    return "".join(lines).encode()


def _all_digests(module, data: bytes) -> tuple:
    text = module._parsed_index(data)
    return text.skeleton, tuple((d.qualname, d.ordinal, d.digest) for d in text.defs)


def test_a2_a_docstring_only_body_digests_as_a_body_of_pass():
    """A2: over every tracked Python file with a docstring-only def or class body, the skeleton and every def's
    digest equal those of the same file with each such docstring replaced, IN THE TEXT, by ``pass``: the strip
    reads a docstring-only body as ``pass``, and the location the fix gives that ``pass`` moves no digest
    (a digest excludes locations). The reference side is a textual edit parsed by ``ast``, not the strip.
    First red: a fix that substitutes another node (``Expr(Constant(None))``)."""
    memo = api.module()
    files = _tracked_python()
    assert len(files) > 1000, len(files)
    moved, witnesses = [], 0
    for path in files:
        data = path.read_bytes()
        docstrings = _docstring_only_bodies(ast.parse(data))
        if not docstrings:
            continue
        witnesses += 1
        if _all_digests(memo, data) != _all_digests(memo, _with_pass_for_docstring_only_bodies(data, docstrings)):
            moved.append(str(path.relative_to(REPO)))
    assert witnesses, "no file has a docstring-only body: the row compares nothing"
    assert not moved, f"{len(moved)} of {witnesses} files digest a docstring-only body unlike pass: {moved[:5]}"
