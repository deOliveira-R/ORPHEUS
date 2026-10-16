r"""#602 round two, B: what a manifest and a validation read from a source file, kept on disk by its content.

Spec ``scratch/characteristic_architecture/p3/m602/round2_gate_spec.md`` §1, rows B1–B14.

A miss used to parse every source file the generation touched, and a process's first validation parsed every
pinned file. The index (:class:`~orpheus.numerics.traced_memo._SourceIndex`) is what those reads take from a
tree: the skeleton digest, each outermost def with its span and digest, and each node a code object can be
compiled from, with its span. It is kept under ``.cache/sources`` and keyed by the file's bytes,
the interpreter and the memo's own source.

**The oracle is the test's own parse (X4).** B1–B4 compare the index a process reads from disk with a parse
this file makes for itself: the docstring strip, the outermost-def walk, the skeleton and the digests are
re-spelled below from their definitions, and B3/B4's ``_def_pin`` tree walk is the retired
``_def_pin``/``_holds`` (``main`` before B). Production's builder is on neither side, so the rows pin the
builder, the serialisation and the queries together; only ``ast`` and ``hashlib`` are shared.

**Only a generating process writes an index** (the ruling after qa's round-two review: a process whose memo
came from a timestamp ``.pyc``, or was edited after the import, wrote an index under new code's name). A row that
needs an index written in THIS process puts it in the generating condition (:func:`_generating`: the running
memo's digest set, as a generation's pinned bytes set it); B16 and B17 gate the rule itself.

Rows that need an isolated index directory monkeypatch ``_SOURCE_INDEXES`` to ``tmp_path`` in THIS process. A
generating child computes its own path from its ``__file__`` (the tree's ``.cache/sources``), so B13 and B14,
which observe a child, read that fixed directory and keep their witnesses unique by content (a module named by
a fresh uuid).
"""
from __future__ import annotations

import ast
import copy
import hashlib
import os
import subprocess
import sys
import textwrap
import types
import uuid
from pathlib import Path
from typing import NamedTuple

import pytest

from . import _traced_memo_api as api
from .test_traced_memo_findings import _Package

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]
_M13B = 'tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3b_nested_code_maps_to_its_outermost_def'
_M21 = 'tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed'
_B1 = 'tests/gates/numerics/test_traced_memo_source_index.py::test_b1_an_index_survives_its_json_round_trip'
_B3 = 'tests/gates/numerics/test_traced_memo_source_index.py::test_b3_def_pin_from_the_index_is_def_pin_from_the_tree'

#: qa's adversarial fixture of 26 defs (``scratch/characteristic_architecture/p3/m602/qa/shapes.py``, ``FIX``):
#: a def under every block, decorators with lambdas, async, ``try*``, ``match`` with a guard, PEP 695, overloads.
QA_SHAPES = '''
"""module doc"""
import functools
@functools.cache
@(lambda f: f)
def deco(x: int = (lambda: 1)()) -> int:
    """d"""
    def inner():
        """i"""
        class InDef:
            """c"""
            def m(self):
                """m"""
                return 1
        return InDef
    return inner
async def coro():
    """a"""
    async with ctx() as c:
        async def nested_async():
            return 1
    async for i in it():
        pass
    return await c
if A:
    def cond(): return 1
elif B:
    def cond(): return 2
else:
    def cond(): return 3
for i in range(3):
    def in_for(): return i
else:
    def for_else(): return 0
while W:
    def in_while(): return 1
else:
    def while_else(): return 2
with open("x") as f, open("y") as g:
    def in_with(): return 1
try:
    def in_try(): return 1
except ValueError as e:
    def in_except(): return 2
else:
    def try_else(): return 3
finally:
    def try_finally(): return 4
try:
    def in_trystar(): return 1
except* ValueError:
    def in_exceptstar(): """s"""; return 2
match y:
    case [1, *rest] if rest:
        class InCase:
            """k"""
            def meth(self): """q"""; return 3
    case _:
        def in_default(): return 4
class Outer:
    """o"""
    X: int = 1
    class Inner:
        """in"""
        @property
        def p(self): """p"""; return 1
        if T:
            def cond_meth(self): return 1
        try:
            def try_meth(self): return 1
        except* E:
            pass
    def generic[T](self, v: T) -> T:
        """only a docstring"""
    def doc_only(self, x: int) -> int:
        """only a docstring"""
squares = [v * v for v in range(3)]
total = sum(v for v in range(3))
'''


# ── the oracle: the test's OWN parse and the retired tree walk (``main`` before B) ─────────────────────
#
# Independent of production's builder (X4): the docstring strip, the outermost-def walk, the skeleton and the
# digests are re-spelled here from their definitions (the module docstring of ``traced_memo``), and
# ``_def_pin``'s tree walk is the retired ``_def_pin``/``_holds``. Only ``ast`` and ``hashlib`` are shared.

_DEFS = (ast.FunctionDef, ast.AsyncFunctionDef)
_STATEMENT_HOLDERS = (ast.stmt, ast.excepthandler, ast.match_case)
_SIGNATURE_SCOPES = ("__annotate__", "<generic parameters of ", "<TypeVar bound of ", "<TypeVar constraint of ",
                     "<TypeAlias ", "<type parameters of ")
_ANONYMOUS = {"<lambda>": ast.Lambda, "<genexpr>": ast.GeneratorExp, "<listcomp>": ast.ListComp,
              "<setcomp>": ast.SetComp, "<dictcomp>": ast.DictComp}


def _hex(node: ast.AST) -> str:
    return hashlib.blake2b(ast.dump(node, include_attributes=False).encode(), digest_size=32).hexdigest()


def _span(node) -> range:
    first = min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])
    return range(first, node.end_lineno + 1)


class _OracleDef(NamedTuple):
    qualname: str
    ordinal: int
    span: range
    node: ast.FunctionDef | ast.AsyncFunctionDef
    digest: str


class _Oracle:
    """One file parsed by the test: the stripped tree, its outermost defs and its skeleton digest."""

    def __init__(self, data: bytes) -> None:
        tree = ast.parse(data)
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, *_DEFS)) and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                    node.body = node.body[1:] or [ast.copy_location(ast.Pass(), first)]
        self.tree = tree
        self.defs: list[_OracleDef] = []
        seen: dict[str, int] = {}

        def walk(parent: ast.AST, prefix: str) -> None:
            for node in ast.iter_child_nodes(parent):
                if isinstance(node, _DEFS):
                    qualname = prefix + node.name
                    seen[qualname] = seen.get(qualname, -1) + 1
                    self.defs.append(_OracleDef(qualname, seen[qualname], _span(node), node, _hex(node)))
                elif isinstance(node, ast.ClassDef):
                    walk(node, f"{prefix}{node.name}.")
                elif isinstance(node, _STATEMENT_HOLDERS):
                    walk(node, prefix)

        walk(tree, "")
        hollow = copy.deepcopy(tree)
        hollow_defs: list[ast.FunctionDef | ast.AsyncFunctionDef] = []

        def collect(parent: ast.AST) -> None:
            for node in ast.iter_child_nodes(parent):
                if isinstance(node, _DEFS):
                    hollow_defs.append(node)
                elif isinstance(node, (ast.ClassDef, *_STATEMENT_HOLDERS)):
                    collect(node)

        collect(hollow)
        for node in hollow_defs:
            node.body = [ast.Pass()]
        self.skeleton = _hex(hollow)

    def enclosing(self, line: int):
        return next((d for d in self.defs if line in d.span), None)


def _tree_holds(scope: ast.AST, line: int, name: str) -> bool:
    if name.startswith(_SIGNATURE_SCOPES):
        return True
    for node in ast.walk(scope):
        if isinstance(node, (*_DEFS, ast.ClassDef)) and node.name == name and line in _span(node):
            return True
        if name in _ANONYMOUS and isinstance(node, _ANONYMOUS[name]) and line in _span(node):
            return True
    return False


def _tree_def_pin(memo, oracle: _Oracle, path: str, code):
    line, name = code.co_firstlineno, code.co_name
    found = oracle.enclosing(line)
    if name == "<module>" or (found is not None and name.startswith(_SIGNATURE_SCOPES) and line < found.node.body[0].lineno):
        return None
    if found is None:
        if _tree_holds(oracle.tree, line, name):
            return None
    elif _tree_holds(found.node, line, name):
        return memo.DefPin(path, found.qualname, found.ordinal, found.digest)
    raise memo.Unpinnable(f"{path}:{line} {code.co_qualname}: no def of that name at its line")


def _outcome(function):
    try:
        return ("pin", function())
    except Exception as error:  # noqa: BLE001 - the outcome (a pin or a refusal) is the compared value
        return (type(error).__name__, str(error))


def _codes(code: types.CodeType):
    yield code
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType):
            yield from _codes(constant)


def _generating(memo, patch, digest: str | None = None) -> None:
    """Put this process in the generating condition: the memo it runs is named (``_RUNNING_MEMO_DIGEST``, which
    a generation sets from the bytes its loader pinned) and its indexes are named by the same digest."""
    digest = digest if digest is not None else memo._memo_digest()
    patch.setattr(memo, "_RUNNING_MEMO_DIGEST", digest)
    patch.setattr(memo, "_INDEX_MEMO_DIGEST", digest)


@pytest.fixture(scope="module")
def population(tmp_path_factory):
    """``(name, bytes, oracle, index)`` for every tracked Python file and qa's fixture. The index is the one a
    process reads from DISK: built and written into a directory of this module's own, the in-memory cache
    cleared, then read back."""
    memo = api.module()
    listed = subprocess.run(["git", "ls-files", "*.py"], cwd=REPO, capture_output=True, text=True, check=True).stdout.split()
    sources = [(p, (REPO / p).read_bytes()) for p in listed if (REPO / p).is_file()]
    sources.append(("<qa shapes>", QA_SHAPES.encode()))
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(memo, "_SOURCE_INDEXES", tmp_path_factory.mktemp("sources"))
        _generating(memo, patch)
        _fresh(memo)
        for _, data in sources:
            memo._source_index(data)
        _fresh(memo)
        rows = [(name, data, _Oracle(data), memo._source_index(data)) for name, data in sources]
        _fresh(memo)
    print(f"population: {len(rows)} files")
    assert len(rows) > 1000, len(rows)
    return rows


def test_b1_an_index_survives_its_json_round_trip(population):
    """B1: the index of every tracked file (and qa's fixture), read back from disk, is rebuilt equal from its
    own JSON. First red: a field written in one order and read in another."""
    memo = api.module()
    index_type = type(population[0][3])
    broken = [name for name, _, _, index in population if index_type.from_json(index.to_json()) != index]
    assert not broken, f"{len(broken)} of {len(population)} indexes changed through JSON: {broken[:5]}"
    assert memo  # the module the population was read through


@pytest.mark.rests_on(_B1)
def test_b2_the_index_digests_are_the_oracles(population):
    """B2: the skeleton and the digest of every outermost def ``(qualname, ordinal)``, read from the index on
    disk, equal the test's own parse (the oracle above). First red: a def's digest wrong or missing in the
    stored index (a key without the bytes serves another file's)."""
    wrong = []
    for name, _, oracle, index in population:
        if index.skeleton != oracle.skeleton:
            wrong.append((name, "<skeleton>"))
        wrong += [(name, d.qualname) for d in oracle.defs if index.digest(d.qualname, d.ordinal) != d.digest]
    assert not wrong, f"{len(wrong)} digests differ: {wrong[:5]}"


@pytest.mark.rests_on(_B1, _M13B)
def test_b3_def_pin_from_the_index_is_def_pin_from_the_tree(population):
    """B3: over every code object ``compile`` makes from every tracked file and qa's fixture (defs, lambdas,
    generator expressions, class bodies, ``__annotate__``, type-parameter scopes), ``_def_pin`` read from the
    index gives the pin, or the refusal, the retired tree walk over the test's own parse gives. First red: a
    query that ignores which def holds a scope."""
    memo = api.module()
    n_codes, wrong = 0, []
    for name, data, oracle, index in population:
        for code in _codes(compile(data, name, "exec")):
            n_codes += 1
            ours = _outcome(lambda: memo._def_pin(index, name, code))
            theirs = _outcome(lambda: _tree_def_pin(memo, oracle, name, code))
            if ours != theirs:
                wrong.append((name, code.co_firstlineno, code.co_name, ours, theirs))
    print(f"B3: {len(population)} files, {n_codes} code objects")
    assert n_codes > 10_000, n_codes
    assert not wrong, f"{len(wrong)} of {n_codes} code objects pinned differently: {wrong[:3]}"


@pytest.mark.rests_on(_B3)
def test_b4_the_span_boundaries_agree_with_the_tree(population):
    """B4: every def, class and anonymous scope of the test's own parse, its name probed at its first line minus
    one, its first, its last and its last plus one (lines no compiled code object need report), gives the
    tree's answer. A span is inclusive in the index and a ``range`` in the tree, so this is where an
    off-by-one lives. First red: spans written exclusive, or the def holding a scope ignored."""
    memo = api.module()
    anonymous = {kind: name for name, kind in _ANONYMOUS.items()}
    probes, wrong = 0, []
    for name, _, oracle, index in population:
        named = {(node.name if isinstance(node, (*_DEFS, ast.ClassDef)) else anonymous[type(node)], _span(node))
                 for node in ast.walk(oracle.tree) if isinstance(node, (*_DEFS, ast.ClassDef, *anonymous))}
        for scope_name, span in named:
            for line in {span.start - 1, span.start, span.stop - 1, span.stop}:
                code = types.SimpleNamespace(co_firstlineno=line, co_name=scope_name, co_qualname=scope_name)
                probes += 1
                ours = _outcome(lambda: memo._def_pin(index, name, code))
                theirs = _outcome(lambda: _tree_def_pin(memo, oracle, name, code))
                if ours != theirs:
                    wrong.append((name, line, scope_name, ours, theirs))
    print(f"B4: {probes} probes")
    assert not wrong, f"{len(wrong)} of {probes} boundary probes differ: {wrong[:3]}"


@pytest.mark.rests_on(_B3)
def test_b5_a_real_manifest_is_the_same_from_a_cold_and_a_warm_index(tmp_path, monkeypatch):
    """B5: one recording of five cheap characteristic blocks makes the same manifest with an empty index
    directory (every index parsed and written) and with the indexes it wrote (every index read, 0 parses).
    First red: an index that does not survive the disk."""
    memo = api.module()
    boot = __import__("orpheus.numerics._traced_memo_boot", fromlist=["Recording"])
    import numpy as np

    from orpheus.derivations.continuous.characteristic import PanelBasis, TransportResolution, Wall, Walls
    from orpheus.derivations.continuous.characteristic.assembly import group_transport
    from orpheus.geometry.chart import Chart
    from orpheus.geometry.chord import ConcentricPartition
    from orpheus.geometry.coord import CoordSystem

    sphere = CoordSystem.SPHERICAL
    walls = lambda n: Walls((Wall(n, 0.0, 1.0, n),), n, Chart(sphere))  # noqa: E731
    cases = [(PanelBasis.of(ConcentricPartition(Chart(sphere), (0.0, 1.0)), 3, 2, 0.4), walls(1), (0.7 + 0.01 * k,), TransportResolution(4)) for k in range(3)]
    cases += [(PanelBasis.of(ConcentricPartition(Chart(sphere), (0.0, 0.5, 1.5, 2.0)), 3, 2, 0.4), walls(3), (0.6 + 0.01 * k, 1.3, 0.45), TransportResolution(4)) for k in range(2)]
    recording = boot.Recording().start()
    try:
        with api.bypass():
            for case in cases:
                group_transport(*case, np.asarray(case[2]) > 0)
    finally:
        recording.stop()
    monkeypatch.setattr(memo, "_SOURCE_INDEXES", tmp_path / "sources")
    _generating(memo, monkeypatch)
    store = memo.Store(tmp_path / "store")
    _fresh(memo)
    cold = memo._manifest(recording, store)
    written = sorted((tmp_path / "sources").glob("*.json"))
    _fresh(memo)
    parses: list[int] = []
    real_parse = ast.parse
    monkeypatch.setattr(ast, "parse", lambda *a, **k: (parses.append(1), real_parse(*a, **k))[1])
    warm = memo._manifest(recording, store)
    monkeypatch.setattr(ast, "parse", real_parse)
    _fresh(memo)
    assert len(cold.functions) > 50 and written, (len(cold.functions), len(written))
    assert parses == [], f"a warm index directory still parsed {len(parses)} files"
    assert (set(cold.functions), set(cold.modules)) == (set(warm.functions), set(warm.modules))


# ── the key, the file and the route ─────────────────────────────────────────────────────────────────


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    """The index directory moved to ``tmp_path`` for this process, the in-memory caches cleared, and a counter
    of ``ast.parse`` calls (the route a claim about the path is gated on, ``vv-principles`` #26)."""
    memo = api.module()
    monkeypatch.setattr(memo, "_SOURCE_INDEXES", tmp_path / "sources")
    _generating(memo, monkeypatch)
    _fresh(memo)
    parses: list[int] = []
    real_parse = ast.parse
    monkeypatch.setattr(ast, "parse", lambda *a, **k: (parses.append(1), real_parse(*a, **k))[1])
    yield memo, parses
    _fresh(memo)


SOURCE = b'def f(x):\n    """Doc."""\n    return x + 1\n\n\nLIMIT = 3\n'


def _fresh(memo) -> None:
    """Every in-memory cache a source text is read through, so the next read is the disk's or a parse."""
    for name in ("_source_index", "_parsed_index"):
        if hasattr(cached := getattr(memo, name, None), "cache_clear"):
            cached.cache_clear()  # type: ignore[union-attr]


def _index_files(memo) -> list[Path]:
    return sorted(Path(memo._SOURCE_INDEXES).glob("*.json"))


def test_b6_an_edit_to_the_bytes_is_another_index(isolated):
    """B6: an index is found again for the same bytes (0 parses) and never for edited bytes (one parse, the new
    digest). First red: a key without the bytes, which serves the first file's index for every file."""
    memo, parses = isolated
    first = memo._source_index(SOURCE)
    _fresh(memo)
    parses.clear()
    assert memo._source_index(SOURCE) == first and parses == []
    edited = SOURCE.replace(b"x + 1", b"x + 2")
    (oracle,) = [d for d in _Oracle(edited).defs if d.qualname == "f"]
    _fresh(memo)
    parses.clear()
    assert memo._source_index(edited).digest("f") == oracle.digest != first.digest("f")
    assert len(parses) == 1, parses
    assert len(_index_files(memo)) == 2


def test_b7_another_interpreter_is_another_index(isolated, monkeypatch):
    """B7: an index written under one interpreter identity is not read under another (one parse, a second
    file). No realisable cross-interpreter state is known: 0 of 390 ``orpheus/`` files differ in ``ast.dump``
    between 3.13.12 and 3.14.3 (``[M]`` 2026-10-15, ``ta/p2_ast_across_pythons.py``); the component stays
    because the dump format is not a stable interface, and this row is mechanical. First red: a key and a
    header without the identity."""
    memo, parses = isolated
    memo._source_index(SOURCE)
    _fresh(memo)
    parses.clear()
    monkeypatch.setattr(memo, "_python_identity", lambda: "3.99.0 (another interpreter)|cpython-399")
    memo._source_index(SOURCE)
    assert len(parses) == 1, "an index written by another interpreter was read"
    assert len(_index_files(memo)) == 2


def test_b8_another_normaliser_is_another_index(isolated, monkeypatch):
    """B8: the memo's source keys the index, since the normaliser that built it lives there. A normaliser that
    keeps docstrings (with the memo digest it would carry) serves its own digest of ``f``, which moves with the
    docstring, and never the index the old normaliser wrote. First red: a key and a header without the memo."""
    memo, parses = isolated
    old = memo._source_index(SOURCE).digest("f")
    _fresh(memo)
    _generating(memo, monkeypatch, "a-normaliser-that-keeps-docstrings")
    monkeypatch.setattr(memo, "_without_docstrings", lambda tree: tree)
    new = memo._source_index(SOURCE).digest("f")
    assert new != old, "the old normaliser's index was served"


def test_b9_an_index_is_the_same_under_minus_o(tmp_path):
    """B9: the index of qa's fixture built in a process with ``-O`` and in one without is the same JSON, which is
    why ``-O`` is not in the key. First red: a build that depends on ``-O``."""
    data = tmp_path / "shapes.py"
    data.write_bytes(QA_SHAPES.encode())
    bodies = []
    for flag in ([], ["-O"]):
        directory = tmp_path / f"sources{len(bodies)}"  # each process builds its own: nothing is read back
        script = ("import sys; from pathlib import Path; from orpheus.numerics import traced_memo as tm; "
                  f"tm._SOURCE_INDEXES = Path({str(directory)!r}); "
                  f"sys.stdout.buffer.write(tm._source_index(Path({str(data)!r}).read_bytes()).to_json())")
        bodies.append(subprocess.run([sys.executable, *flag, "-c", script], cwd=REPO, capture_output=True, check=True).stdout)
    assert bodies[0] and bodies[0] == bodies[1]


def _header(raw: bytes) -> dict[str, str]:
    """The key fields an index file's first line names (its body's digest left out)."""
    import json

    return {k: v for k, v in json.loads(raw.split(b"\n", 1)[0]).items() if k != "body"}


def _verdict(memo, path: Path, header: dict[str, str]) -> str:
    """What ``_read_index`` makes of ``path``: ``"index"``, ``"none"`` or ``"corrupt"`` (anything else it returns)."""
    found = memo._read_index(path, header)
    return "none" if found is None else "index" if hasattr(found, "skeleton") else "corrupt"


def _corrupt(memo, data: bytes) -> dict[str, bytes]:
    """The index file of ``data`` as written, and each way it is corrupted."""
    memo._source_index(data)
    (path,) = _index_files(memo)
    raw = path.read_bytes()
    head, body = raw.split(b"\n", 1)
    skeleton = memo._source_index(data).skeleton.encode()
    flipped = skeleton[:5] + (b"0" if skeleton[5:6] != b"0" else b"1") + skeleton[6:]
    other = SOURCE.replace(b"LIMIT = 3", b"LIMIT = 4")
    _fresh(memo)
    memo._source_index(other)
    (other_path,) = [p for p in _index_files(memo) if p != path]
    copied = other_path.read_bytes()
    other_path.unlink()
    import json

    header = json.loads(head)
    return {"as-written": raw, "empty": b"", "one-byte": raw[:1], "half": raw[: len(raw) // 2], "last-byte-missing": raw[:-1],
            "digest-hex-flipped": head + b"\n" + body.replace(skeleton, flipped, 1), "other-file": copied,
            "header-identity": json.dumps({**header, "python": "another"}).encode() + b"\n" + body}


@pytest.mark.parametrize("kind", ["empty", "one-byte", "half", "last-byte-missing", "digest-hex-flipped", "other-file", "header-identity"])
def test_b10_a_corrupt_index_is_rebuilt_never_read(isolated, kind):
    """B10: a truncated index file, one hex digit of its skeleton digest flipped (still valid JSON), the index of
    other bytes copied under this key, and a header naming another interpreter are each a ``_CorruptIndex``; the
    next read gives the parse's digests, and the file holds a valid index again. First red: the body's digest
    not checked (the flipped digit is served), or the header not checked (the copied index is served)."""
    memo, _ = isolated
    variants = _corrupt(memo, SOURCE)
    _fresh(memo)
    (path,) = _index_files(memo)
    path.write_bytes(variants[kind])
    header = _header(variants["as-written"])
    assert _verdict(memo, path, header) == "corrupt", kind
    oracle = _Oracle(SOURCE)
    index = memo._source_index(SOURCE)
    assert (index.skeleton, index.digest("f")) == (oracle.skeleton, oracle.defs[0].digest), "a corrupt index was served"
    assert _verdict(memo, path, header) == "index", "the corrupt index was not rewritten"


@pytest.mark.timeout(300)
def test_b11_concurrent_writers_never_expose_a_partial_index(tmp_path, monkeypatch):
    """B11: four processes each write one key's index 200 times (removing it before each write) while this
    process reads it. Every read is the index or nothing; 0 reads are a ``_CorruptIndex``. A corrupt read is
    rebuilt silently, so the pins are blind to a write in place, and the typed verdict is the instrument.
    First red: an index written straight to its name, with no ``os.replace``."""
    memo = api.module()
    directory = tmp_path / "sources"
    data = (SOURCE * 40) + b"\n".join(b"def g%d(x):\n    return x * %d\n" % (i, i) for i in range(300))
    (tmp_path / "data.py").write_bytes(data)
    monkeypatch.setattr(memo, "_SOURCE_INDEXES", directory)
    _generating(memo, monkeypatch)
    _fresh(memo)
    memo._source_index(data)
    _fresh(memo)
    (path,) = sorted(directory.glob("*.json"))
    header = _header(path.read_bytes())
    writer = textwrap.dedent(f'''
        import os
        from pathlib import Path
        from orpheus.numerics import traced_memo as tm
        tm._SOURCE_INDEXES = Path({str(directory)!r})
        tm._RUNNING_MEMO_DIGEST = tm._INDEX_MEMO_DIGEST  # the generating condition: only a generation writes
        data = Path({str(tmp_path / "data.py")!r}).read_bytes()
        for _ in range(200):
            tm._source_index.cache_clear()
            try:
                os.remove({str(path)!r})
            except FileNotFoundError:
                pass
            tm._source_index(data)
        print("done")
    ''')
    procs = [subprocess.Popen([sys.executable, "-c", writer], cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(4)]
    reads = corrupt = 0
    while any(p.poll() is None for p in procs):
        verdict = _verdict(memo, path, header)
        if verdict != "none":
            reads += 1
            corrupt += verdict == "corrupt"
    for p in procs:
        out, err = p.communicate()
        assert p.returncode == 0 and "done" in out, err[-2000:]
    assert reads > 0, "the reader never found the index: the instrument read nothing"
    assert corrupt == 0, f"{corrupt} of {reads} reads found a partial index"


_B12_VALIDATE = '''
import ast, json, sys
from pathlib import Path
sys.path.insert(0, {root!r})
from orpheus.numerics import traced_memo as tm
tm._SOURCE_INDEXES = Path({sources!r})
parses = []
real = ast.parse
ast.parse = lambda *a, **k: (parses.append(1), real(*a, **k))[1]
reasons = tm.validate(tm.Manifest.from_json(json.loads(Path({manifest!r}).read_text())))
print(json.dumps({{"parses": len(parses), "reasons": list(reasons)}}))
'''


@pytest.mark.rests_on(_M21)
def test_b12_a_processs_first_validation_reads_indexes_instead_of_parsing(tmp_path, monkeypatch):
    """B12: a fresh process validating a manifest parses nothing when the indexes of the pinned files exist, and
    parses each distinct pinned file once when they do not; both validations are clean. Both sides of the
    branch are gated. First red: a validation that does not consult the index."""
    memo = api.module()
    package = _Package(tmp_path, m="def f(x):\n    return g(x) + 1\n\n\ndef g(x):\n    return x * 2\n")
    try:
        _, manifest = api.trace_call(package.module().f, 1.0)
        files = {path for path, *_ in manifest.functions} | {path for path, _ in manifest.modules}
        assert files, "the manifest pins no source"
        (tmp_path / "manifest.json").write_text(__import__("json").dumps(manifest.to_json()))
        warm, cold = tmp_path / "warm", tmp_path / "cold"
        monkeypatch.setattr(memo, "_SOURCE_INDEXES", warm)
        _generating(memo, monkeypatch)  # written under the disk source's name, which a fresh process reads by
        _fresh(memo)
        assert api.validate(manifest) == ()
        _fresh(memo)
        results = []
        for directory in (warm, cold):
            script = _B12_VALIDATE.format(root=str(package.root), sources=str(directory), manifest=str(tmp_path / "manifest.json"))
            run = subprocess.run([sys.executable, "-c", script], cwd=REPO, capture_output=True, text=True)
            assert run.returncode == 0, run.stderr[-3000:]
            results.append(__import__("json").loads(run.stdout.splitlines()[-1]))
    finally:
        package.close()
    assert results[0] == {"parses": 0, "reasons": []}, results[0]
    assert results[1] == {"parses": len(files), "reasons": []}, (results[1], sorted(files))


def _inodes(directory: Path) -> dict[str, int]:
    return {p.name: p.stat().st_ino for p in directory.glob("*.json")} if directory.is_dir() else {}


@pytest.mark.rests_on(_B3)
def test_b13_a_miss_indexes_only_bytes_it_has_not_indexed(tmp_path):
    """B13, through the generating process (which uses the tree's fixed ``.cache/sources``): the first miss of
    a module with new bytes writes that module's index; a second miss under another key reads every index it
    needs and writes or replaces none (the directory's names and inodes are unchanged; ``os.replace`` would
    give a replaced file a new inode). First red: an index the child writes and never reads."""
    memo = api.module()
    directory = REPO / ".cache" / "sources"  # where a generating process of this tree writes
    witness = f"# a content no other module has: {uuid.uuid4().hex}\n"
    package = _Package(tmp_path, m=witness + "from orpheus.numerics.traced_memo import traced_memo\n\n\n"
                       "@traced_memo\ndef f(x: float) -> float:\n    return x + 1.0\n")
    try:
        with api.cache_root(tmp_path / "cache"):
            m = package.module()
            assert m.f(1.0) == 2.0
            digest = hashlib.blake2b((package.dir / "m.py").read_bytes(), digest_size=32).hexdigest()
            assert any(digest in p.read_bytes().split(b"\n", 1)[0].decode() for p in directory.glob("*.json")), \
                "the first miss did not index the module's new bytes"
            before = _inodes(directory)
            assert m.f(2.0) == 3.0
            after = _inodes(directory)
    finally:
        package.close()
    assert after == before, f"the second miss wrote {len(set(after.items()) - set(before.items()))} index files"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_findings.py::test_a_parent_does_not_pin_its_childs_sources_as_data')
def test_b14_the_index_and_the_memos_own_source_are_no_dependency_of_a_run(tmp_path):
    """B14: a generation that validates a child entry (so reads indexes) has no pin under ``.cache/sources``,
    and no data pin on the memo's own source (which the running memo digests at its import, inside the
    generation). First red: that import-time digest recorded as the run's read."""
    memo = api.module()
    package = _Package(tmp_path, m='''
        from orpheus.numerics.traced_memo import traced_memo
        from .n import inner

        @traced_memo
        def outer(x: float) -> float:
            return inner(x) + 1.0
    ''', n='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def inner(x: float) -> float:
            return x * 5.0
    ''')
    try:
        with api.cache_root(tmp_path / "cache"):
            m = package.module()
            m.inner(2.0)
            assert m.outer(2.0) == 11.0
    finally:
        package.close()
    parent = next(f for f in api.entry_files(tmp_path / "cache") if api.entry_function(f).endswith(":outer"))
    record = api.read_entry(parent)
    sources = os.path.realpath(REPO / ".cache" / "sources")
    memo_sources = {os.path.realpath(REPO / "orpheus" / "numerics" / f) for f in ("traced_memo.py", "_traced_memo_boot.py", "content.py")}
    paths = [p for kind in ("DataPin", "PresencePin", "ListingPin") for p, *_ in api.manifest_rows(record, kind)]
    assert not [p for p in paths if os.path.realpath(p).startswith(sources)], paths
    data = [p for p, *_ in api.manifest_rows(record, "DataPin")]
    # A presence pin on the memo's file is the memo's import resolving its own path (``Path(__file__).resolve()``
    # probes each ancestor); it holds a kind, never bytes, so it does not move with an edit. A DATA pin would.
    assert not [p for p in data if os.path.realpath(p) in memo_sources], data
    assert api.manifest_rows(record, "ChildPin"), "the generation read no child entry: the row exercised nothing"



_MEMO_EDIT_DRIVER = '''
import os, sys, threading, time
sys.path[:0] = [{tree!r}]
import orpheus
assert orpheus.__file__.startswith({tree!r}), orpheus.__file__
from orpheus.numerics import traced_memo as tm
import memo_edit_fixture.m as m
outcome = []
def call():
    try:
        outcome.append(m.slow(1.5))
    except Exception as error:
        outcome.append(error)
with tm.cache_root({cache!r}):
    thread = threading.Thread(target=call)
    thread.start()
    deadline = time.monotonic() + 120
    while not os.path.exists({ready!r}):
        assert time.monotonic() < deadline
        time.sleep(0.01)
    with open(os.path.join({tree!r}, "orpheus", "numerics", "traced_memo.py"), "a") as handle:
        handle.write("# an edit to the memo while a generation runs\\n")
    open({go!r}, "w").close()
    thread.join(120)
    print("OUTCOME", type(outcome[0]).__name__, "VERDICT", type(m.slow.lookup(1.5)).__name__)
'''


@pytest.mark.timeout(600)
@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_compiled_bytes.py::test_p2_an_edit_landing_during_a_generation_leaves_the_entry_stale')
def test_b15_the_memos_own_source_edited_during_a_generation_leaves_the_entry_stale(tmp_path):
    """B15 (the elegance review's F4, the shape of P2): the memo's own source is edited while a generation runs.
    The generation ran the old memo, so its entry pins the memo as that process imported it, and the entry is
    ``Stale`` against the edited source, never served. Run on a copy of the package's Python files (the memo
    under test is the copy's), in a driver process. First red: the memo pinned by re-reading its source after
    the run, which pins the edited bytes and makes the entry a ``Hit``."""
    import shutil

    tree = tmp_path / "tree"
    shutil.copytree(REPO / "orpheus", tree / "orpheus",
                    ignore=lambda d, names: [n for n in names if not n.endswith(".py") and not (Path(d) / n).is_dir() or n == "__pycache__"])
    ready, go = tmp_path / "ready", tmp_path / "go"
    fixture = tree / "memo_edit_fixture"
    fixture.mkdir()
    (fixture / "__init__.py").write_text("")
    (fixture / "m.py").write_text(textwrap.dedent(f'''
        import os
        import time

        from orpheus.numerics.traced_memo import traced_memo


        @traced_memo
        def slow(x: float) -> float:
            open({str(ready)!r}, "w").close()
            while not os.path.exists({str(go)!r}):
                time.sleep(0.01)
            return x * 3.0
    '''))
    driver = tmp_path / "driver.py"
    driver.write_text(_MEMO_EDIT_DRIVER.format(tree=str(tree), cache=str(tmp_path / "cache"), ready=str(ready), go=str(go)))
    run = subprocess.run([sys.executable, *(["-O"] if sys.flags.optimize else []), str(driver)], cwd=tmp_path,
                         capture_output=True, text=True, timeout=500)
    assert run.returncode == 0, run.stderr[-3000:]
    line = run.stdout.split("OUTCOME", 1)[1].split()
    outcome, verdict = line[0], line[2]
    assert verdict == "Stale", f"the entry generated by the old memo is {verdict} against the edited memo (outcome {outcome})"


# ── who writes an index (qa's round-two review, F1) and what an unreadable one is ────────────────────


def test_b16_a_process_that_is_not_generating_never_writes_an_index(tmp_path, monkeypatch):
    """B16: outside a generation (no running memo digest: the memo may have come from a timestamp ``.pyc`` or
    have been edited after the import) a source text's index is parsed in memory and served, and no file is
    written; the same call in the generating condition writes one. Both sides of the rule are gated. First
    red: an index written whatever the process."""
    memo = api.module()
    monkeypatch.setattr(memo, "_SOURCE_INDEXES", tmp_path / "sources")
    monkeypatch.setattr(memo, "_RUNNING_MEMO_DIGEST", None)
    _fresh(memo)
    oracle = _Oracle(SOURCE)
    assert (memo._source_index(SOURCE).skeleton, memo._source_index(SOURCE).digest("f")) == (oracle.skeleton, oracle.defs[0].digest)
    assert _index_files(memo) == [], "a process that is not generating wrote an index"
    _generating(memo, monkeypatch)
    _fresh(memo)
    memo._source_index(SOURCE)
    assert len(_index_files(memo)) == 1, "the generating condition wrote no index: the row's other side is blind"
    _fresh(memo)


_B17_PROCESS = '''
import sys
sys.path.insert(0, sys.argv[1])
import orpheus.numerics.traced_memo as tm
assert tm.__file__.startswith(sys.argv[1]), tm.__file__
data = open(sys.argv[2], "rb").read()
served, parsed = tm._source_index(data).digest("f"), tm._parsed_index(data).digest("f")
print("SERVED", served, "PARSED", parsed)
'''


def test_b17_old_memo_code_cannot_write_an_index_under_new_codes_name(tmp_path):
    """B17, qa's p2 scenario as a gate. On a copy of the memo's four files: process A imports it (a timestamp
    ``.pyc``); the normaliser is edited at the same size to keep docstrings, with the mtime put back; process B
    imports the OLD code from the ``.pyc`` while the memo's source on disk is the NEW one, and reads a sample's
    index; process C compiles the NEW code fresh and reads the sample's index. C is served its own
    normalisation's digest. First red: B, running old code, writes its index under the new source's name, and
    C reads the old normalisation's digest as its own."""
    import shutil

    tree = tmp_path / "tree"
    package = tree / "orpheus" / "numerics"
    package.mkdir(parents=True)
    (tree / "orpheus" / "__init__.py").write_text("")
    (package / "__init__.py").write_text("")
    for name in ("traced_memo.py", "_traced_memo_boot.py", "content.py", "scalars.py"):
        shutil.copy2(REPO / "orpheus" / "numerics" / name, package)
    memo_file = package / "traced_memo.py"
    old, new = b"isinstance(first.value.value, str)", b"isinstance(first.value.value, int)"
    source = memo_file.read_bytes()
    assert source.count(old) == 1, "the edit target moved: re-aim the row"
    sample = tmp_path / "sample.py"
    sample.write_text('def f():\n    "a docstring"\n    return 1\n')
    script = tmp_path / "process.py"
    script.write_text(_B17_PROCESS)

    def run() -> tuple[str, str]:
        done = subprocess.run([sys.executable, *(["-O"] if sys.flags.optimize else []), str(script), str(tree), str(sample)],
                              cwd=tmp_path, capture_output=True, text=True)
        assert done.returncode == 0, done.stderr[-3000:]
        words = done.stdout.split()
        return words[words.index("SERVED") + 1], words[words.index("PARSED") + 1]

    run()  # A: writes the timestamp .pyc of the original memo
    stat = memo_file.stat()
    memo_file.write_bytes(source.replace(old, new))
    os.utime(memo_file, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert memo_file.stat().st_size == stat.st_size
    shutil.rmtree(tree / ".cache", ignore_errors=True)
    b_served, b_parsed = run()  # B: the old code, from the .pyc, under the new source's name
    shutil.rmtree(package / "__pycache__")
    c_served, c_parsed = run()  # C: the new code, compiled fresh
    assert b_parsed != c_parsed, "the two normalisations agree on the sample: the row cannot see a wrong index"
    assert c_served == c_parsed, "a process running old memo code wrote an index under the new code's name"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_source_index.py::test_b12_a_processs_first_validation_reads_indexes_instead_of_parsing')
def test_b18_a_validation_run_inside_a_recording_does_not_pin_the_index(isolated):
    """B18: the one caller that reaches the index while a recording is active outside the store's own lookup is a
    recorded run that validates a manifest itself (``trace_call`` of a function calling :func:`validate`). Its
    read of the index directory is the memo's, not the run's: no pin names a file under it. First red: the
    index read outside ``unrecorded()``."""
    memo, _ = isolated
    package = _Package(Path(memo._SOURCE_INDEXES).parent, m="def f(x):\n    return x + 1\n")
    try:
        _, manifest = api.trace_call(package.module().f, 1.0)
        _fresh(memo)
        assert api.validate(manifest) == ()  # writes the indexes (the generating condition)
        assert _index_files(memo), "no index was written: the row reads nothing"
        _fresh(memo)
        _, outer = api.trace_call(lambda: memo.validate(manifest))
    finally:
        package.close()
    directory = os.path.realpath(memo._SOURCE_INDEXES)
    pinned = [p for p, *_ in [*outer.data, *outer.of(memo.PresencePin), *outer.of(memo.ListingPin)]]
    assert not [p for p in pinned if os.path.realpath(p).startswith(directory)], pinned


def test_b19_an_unreadable_index_is_corrupt_and_the_parse_is_served(isolated):
    """B19: an index file that cannot be read (``PermissionError``) is a corrupt index, never an error the caller
    sees, and the digests served are the parse's. First red: an ``OSError`` other than a missing file escaping
    the read."""
    memo, _ = isolated
    memo._source_index(SOURCE)
    (path,) = _index_files(memo)
    header = _header(path.read_bytes())
    _fresh(memo)
    path.chmod(0)
    try:
        if os.access(path, os.R_OK):
            pytest.skip("this user can read a mode-0 file (root): no unreadable state exists here")
        assert _verdict(memo, path, header) == "corrupt"
        oracle = _Oracle(SOURCE)
        index = memo._source_index(SOURCE)
        assert (index.skeleton, index.digest("f")) == (oracle.skeleton, oracle.defs[0].digest)
    finally:
        path.chmod(0o644)
