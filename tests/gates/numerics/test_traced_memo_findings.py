r"""The qa review of the traced memo (#405 P3, 2026-10-04): one gate per finding, each red on ``99897d10``.

qa reproduced 13 findings against the first build of :mod:`orpheus.numerics.traced_memo` (the list is
``scratch/reference_architecture/p3/qa/findings.md``; the reproducers ``probes/probe_findings.py``). Seven
served a stale or a wrong value, which is the one failure a cache must never have. Each gate below is that
reproducer turned into an assertion of the right behaviour; the finding's number is in the gate's name.
Every gate writes its own generator package under ``tmp_path`` and edits only that copy.
"""
from __future__ import annotations

import importlib
import importlib.util
import os
import subprocess
import sys
import textwrap
import threading
import uuid
from pathlib import Path

import numpy as np
import pytest

from . import _traced_memo_api as api

pytestmark = pytest.mark.foundation


class _Package:
    """A generator package written from source strings under ``tmp_path``, on the front of ``sys.path``."""

    def __init__(self, tmp_path: Path, **modules: str) -> None:
        self.root = tmp_path / "src"
        self.name = f"memo_finding_{uuid.uuid4().hex[:12]}"
        self.dir = self.root / self.name
        self.dir.mkdir(parents=True)
        (self.dir / "__init__.py").write_text("")
        for module, source in modules.items():
            (self.dir / f"{module}.py").write_text(textwrap.dedent(source))
        sys.path.insert(0, str(self.root))
        importlib.invalidate_caches()

    def module(self, name: str = "m"):
        return importlib.import_module(f"{self.name}.{name}")

    def edit(self, old: str, new: str, module: str = "m") -> None:
        path = self.dir / f"{module}.py"
        text = path.read_text()
        if text.count(old) != 1:
            raise AssertionError(f"edit target {old!r} occurs {text.count(old)} times")
        path.write_text(text.replace(old, new))

    def close(self) -> None:
        while str(self.root) in sys.path:
            sys.path.remove(str(self.root))
        for key in [k for k in sys.modules if k == self.name or k.startswith(self.name + ".")]:
            del sys.modules[key]


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


def _kind(verdict) -> str:
    return api.verdict_kind(verdict)


def test_q1_the_def_that_ran_is_pinned_among_defs_of_one_name(make):
    """Finding 1: a property getter beside its setter, and alternative defs under ``if``/``else``, share a
    qualified name; the pin is the def that RAN (by its rank), so an edit to it is ``Stale`` (it was a ``Hit``
    serving 11.0 where a fresh call gives 12.0)."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        class Box:
            @property
            def x(self):
                return 1.0

            @x.setter
            def x(self, value):
                raise AttributeError

        if True:
            def branch():
                return 10.0
        else:
            def branch():
                return -10.0

        @traced_memo
        def gen(a: float) -> float:
            return Box().x + branch() + a
    ''')
    m = package.module()
    assert m.gen(0.0) == 11.0
    package.edit("return 1.0", "return 2.0")
    assert _kind(m.gen.lookup(0.0)) == "Stale" and m.gen(0.0) == 12.0
    package.edit("return 10.0", "return 20.0")
    assert _kind(m.gen.lookup(0.0)) == "Stale" and m.gen(0.0) == 22.0


def test_q2_a_memo_called_from_a_worker_thread_is_pinned_as_a_child(make, tmp_path):
    """Finding 2: the recording and the store are process state, so a memo called from a worker thread of a
    generation is the generation's child and is written under the caller's store (it was neither)."""
    package = make(m='''
        from concurrent.futures import ThreadPoolExecutor
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def inner(x: float) -> float:
            return x * 5.0

        @traced_memo
        def outer(x: float) -> float:
            with ThreadPoolExecutor(1) as pool:
                return pool.submit(inner, x).result() + 1.0
    ''')
    m = package.module()
    assert m.outer(2.0) == 11.0
    functions = {api.entry_function(f).rsplit(":", 1)[-1] for f in api.entry_files(tmp_path / "cache")}
    assert functions == {"inner", "outer"}, functions
    package.edit("return x * 5.0", "return x * 7.0")
    assert _kind(m.outer.lookup(2.0)) == "Stale" and m.outer(2.0) == 15.0


@pytest.mark.parametrize("function", ["pair", "mode", "masked"])
def test_q3_a_subclass_of_a_payload_type_is_refused_never_written_as_its_base(make, function):
    """Finding 3: a ``NamedTuple``, an ``IntEnum`` and a masked array are subclasses of types the payload
    writes; written as their bases they came back as a plain tuple, an ``int`` and an array without its mask
    (mean 2.0 directly, -3.3e29 through the memo). They are refused, and nothing is written."""
    package = make(m='''
        import enum
        from typing import NamedTuple
        import numpy as np
        from orpheus.numerics.traced_memo import traced_memo

        class Pair(NamedTuple):
            k: float
            n: int

        class Mode(enum.IntEnum):
            FAST = 1

        @traced_memo
        def pair(a: float) -> Pair:
            return Pair(a, 3)

        @traced_memo
        def mode(a: float) -> Mode:
            return Mode.FAST

        @traced_memo
        def masked(a: float) -> np.ndarray:
            return np.ma.masked_less(np.array([a, -1e30, 3.0]), -1.0)
    ''')
    memo = getattr(package.module(), function)
    with pytest.raises(api.name("Unencodable")):
        memo(1.0)
    assert _kind(memo.lookup(1.0)) == "Absent"


def test_q4_the_key_is_exact_below_the_top_level(make):
    """Finding 4: the key followed ``==``, so ``sign(-0.0)`` was served ``sign(0.0)``'s 1.0, and
    ``count(Spec(8.0))`` was served ``Spec(8)``'s entry where the direct call raises. The key is the exact
    encoding: the sign of zero and a field's type are part of it."""
    package = make(m='''
        import dataclasses, math
        from orpheus.numerics.traced_memo import traced_memo

        @dataclasses.dataclass(frozen=True)
        class Spec:
            n: int

        @traced_memo
        def sign(x: float) -> float:
            return math.copysign(1.0, x)

        @traced_memo
        def count(spec: Spec) -> int:
            return len(range(spec.n))
    ''')
    m = package.module()
    assert m.sign(0.0) == 1.0 and m.sign(-0.0) == -1.0
    assert m.count.key(m.Spec(8)) != m.count.key(m.Spec(8.0))
    assert m.count(m.Spec(8)) == 8
    with pytest.raises(TypeError):
        m.count(m.Spec(8.0))


def test_q5_absence_a_listing_and_a_bytes_path_are_dependencies(make):
    """Finding 5: an answer that depended on a file's ABSENCE (an ``open`` that failed), on a directory
    listing, or on a file opened by a bytes path was served stale after the change; each is now pinned. (An
    absence probed with ``exists()`` raises no audit event: the declared blind spot of the module docstring;
    this generator probes by opening.)"""
    package = make(m='''
        import os
        from pathlib import Path
        from orpheus.numerics.traced_memo import traced_memo

        HERE = Path(__file__).parent

        @traced_memo
        def override(a: float) -> float:
            try:
                return float((HERE / "override.txt").read_text())
            except FileNotFoundError:
                return a

        @traced_memo
        def listing(a: float) -> float:
            return a + len([f for f in os.listdir(HERE / "tables") if f.endswith(".dat")])

        @traced_memo
        def bytes_open(a: float) -> float:
            with open(os.fsencode(HERE / "t.dat"), "rb") as handle:
                return a + float(handle.read())
    ''')
    (package.dir / "tables").mkdir()
    (package.dir / "tables" / "a.dat").write_text("1")
    (package.dir / "t.dat").write_text("1")
    m = package.module()
    assert (m.override(1.0), m.listing(1.0), m.bytes_open(1.0)) == (1.0, 2.0, 2.0)
    (package.dir / "override.txt").write_text("99")
    (package.dir / "tables" / "b.dat").write_text("1")
    (package.dir / "t.dat").write_text("50")
    for memo in (m.override, m.listing, m.bytes_open):
        assert _kind(memo.lookup(1.0)) == "Stale", memo.__name__
    assert (m.override(1.0), m.listing(1.0), m.bytes_open(1.0)) == (99.0, 3.0, 51.0)


def test_q5_a_relative_read_pins_the_working_directory(make, tmp_path, monkeypatch):
    """Finding 5, the working directory: ``open("table.txt")`` read from one directory was served in another.
    A run that read a relative path pins its working directory."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def table(a: float) -> float:
            with open("table.txt") as handle:
                return a + float(handle.read())
    ''')
    one, two = tmp_path / "one", tmp_path / "two"
    for directory, value in ((one, "1"), (two, "2")):
        directory.mkdir()
        (directory / "table.txt").write_text(value)
    m = package.module()
    monkeypatch.chdir(one)
    assert m.table(0.0) == 1.0
    monkeypatch.chdir(two)
    assert _kind(m.table.lookup(0.0)) == "Stale" and m.table(0.0) == 2.0


def test_q6_a_def_under_any_module_level_block_is_pinned(make):
    """Finding 6: defs under a module-level ``for`` and ``match`` were pinned by nothing (the skeleton holds no
    body), so an edit was served stale (101.0 instead of 105.0)."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        for _i in range(1):
            def looped(x):
                return x + 1.0

        match 1:
            case 1:
                def matched(x):
                    return x + 100.0

        @traced_memo
        def gen(a: float) -> float:
            return looped(a) + matched(a)
    ''')
    m = package.module()
    assert m.gen(0.0) == 101.0
    package.edit("return x + 1.0", "return x + 5.0")
    assert _kind(m.gen.lookup(0.0)) == "Stale" and m.gen(0.0) == 105.0


def test_q7_a_hit_serves_the_bytes_it_verified(tmp_path):
    """Finding 7: the value served was read a second time after the digest was checked, so a write in between
    served a value no generation produced. A hit carries the verified payload."""
    store = api.name("Store")(tmp_path / "cache")
    manifest = api.name("Manifest")(frozenset())
    store.write("f:x", "key", manifest, (1.0, np.array([1.0, 1.0])), set())
    hit = store.lookup("f:x", "key")
    store.write("f:x", "key", manifest, (2.0, np.array([2.0, 2.0])), set())
    first, array = api.decode_payload(hit.tree, hit.arrays, set())
    assert first == 1.0 and array.tolist() == [1.0, 1.0]


def test_q8_the_memo_pins_its_own_source(make, tmp_path):
    """Finding 8: the memo writes the entry after the recording stops, so its own code was in no manifest; an
    entry written by a defective memo stayed a hit after the fix. Every entry pins the memo's source."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def third(a: float) -> float:
            return a / 3.0
    ''')
    m = package.module()
    m.third(1.0)
    (entry,) = api.entry_files(tmp_path / "cache")
    (memo_pin,) = api.manifest_rows(api.read_entry(entry), "MemoPin")
    manifest = api.name("Manifest").from_json(api.read_entry(entry)["manifest"])
    pin = api.name("MemoPin")
    reasons = api.validate(manifest.replacing(pin, [pin("0" * 64)]))
    assert memo_pin and any("own source" in r for r in reasons), reasons


def test_q9_a_print_in_a_generator_never_corrupts_the_answer(make):
    """Finding 9: standard output was the answer's channel, so a ``print(..., flush=True)`` in a generator
    (``data/macro_xs/mixture.py`` prints so) made the caller raise although the entry was written, and a
    raising generator that printed lost its own exception."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def gen(a: float) -> float:
            print("  Sigma-zero iterations...", end=" ", flush=True)
            return a * 3.0

        @traced_memo
        def bad(a: float) -> float:
            print("progress", flush=True)
            raise ValueError("the real refusal")
    ''')
    m = package.module()
    assert m.gen(1.0) == 3.0
    with pytest.raises(ValueError, match="the real refusal"):
        m.bad(1.0)


def test_q10_concurrent_writes_of_one_key_are_atomic(tmp_path):
    """Finding 10: two writers of one key failed 458 of 800 times (``rmtree`` then ``os.replace`` of a
    directory) and leaked their staging directories. An entry is one file replaced atomically."""
    store = api.name("Store")(tmp_path / "cache")
    manifest = api.name("Manifest")(frozenset())
    errors: list[BaseException] = []

    def writer() -> None:
        for _ in range(50):
            try:
                store.write("f:x", "key", manifest, (np.zeros(20000),), set())
            except BaseException as error:  # noqa: BLE001 - every failure is the finding
                errors.append(error)

    threads = [threading.Thread(target=writer) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors, errors[:2]
    assert [p.name for p in (tmp_path / "cache" / "f:x").iterdir()] == ["key.npz"]
    assert _kind(store.lookup("f:x", "key")) == "Hit"


def test_q11_an_exception_that_cannot_be_rebuilt_crosses_described(make):
    """Finding 11 (a): an exception that pickles in the generating process but cannot be rebuilt here (its
    ``__init__`` takes other arguments than its ``args``) arrived anonymous; it crosses as a ``RuntimeError``
    carrying the original type, message and traceback."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        class Refusal(Exception):
            def __init__(self, what, why):
                super().__init__(f"{what}: {why}")

        @traced_memo
        def gen(a: float) -> float:
            raise Refusal("gen", "outside the domain")
    ''')
    with pytest.raises(RuntimeError, match=r"Refusal: gen: outside the domain"):
        package.module().gen(1.0)


def test_q11_an_undeclared_dataclass_return_is_refused_before_it_is_written(make):
    """Finding 11 (b): a dataclass returned by a function whose return annotation does not name it was
    written, then refused on every load. It is refused at the write, and nothing is cached."""
    package = make(m='''
        import dataclasses
        from orpheus.numerics.traced_memo import traced_memo

        @dataclasses.dataclass(frozen=True)
        class R:
            k: float

        @traced_memo
        def gen(a: float):
            return R(a)
    ''')
    m = package.module()
    with pytest.raises(api.name("Unencodable")):
        m.gen(1.0)
    assert _kind(m.gen.lookup(1.0)) == "Absent"


def test_q12_a_module_on_no_sys_path_entry_is_pinned_by_its_absolute_path(tmp_path):
    """Finding 12: a process importing ``orpheus`` through the editable install's finder (a script, a
    notebook) has the package on no ``sys.path`` entry, and the manifest refused it. Such a file is pinned by
    its absolute path, validates, and goes stale when edited."""
    source = tmp_path / "loose" / "loose_module.py"
    source.parent.mkdir()
    source.write_text("def f(x):\n    return x + 1.0\n")
    spec = importlib.util.spec_from_file_location(f"loose_{uuid.uuid4().hex[:8]}", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _, manifest = api.trace_call(module.f, 1.0)
    assert str(source.resolve()) in {path for path, *_ in manifest.functions}
    assert api.validate(manifest) == ()
    source.write_text("def f(x):\n    return x + 2.0\n")
    assert any("changed" in r for r in api.validate(manifest))


@pytest.mark.timeout(600)
def test_q13_a_memo_calling_itself_with_its_own_arguments_is_refused(make):
    """Finding 13: ``f(x)`` calling ``f(x)`` would start interpreters without end; the generating process
    knows the calls it is answering and refuses the cycle."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def f(x: float) -> float:
            return f(x)
    ''')
    with pytest.raises(RecursionError, match="calls itself"):
        package.module().f(1.0)


def test_a_started_program_is_pinned_and_a_started_python_is_refused(make, tmp_path):
    """A process the run starts: a declared machine query is pinned by its executable's bytes
    (``platform``'s processor query runs ``uname``, so every real client starts one); a Python interpreter
    is refused, since the code it ran is code no recording saw."""
    package = make(m=f'''
        import subprocess, sys
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def external(a: float) -> float:
            subprocess.run(["uname"], capture_output=True, check=True)
            return a

        @traced_memo
        def python(a: float) -> float:
            subprocess.run([sys.executable, "-c", "pass"], check=True)
            return a
    ''')
    m = package.module()
    assert m.external(1.0) == 1.0
    (entry,) = api.entry_files(tmp_path / "cache")
    assert any(path.endswith("/uname") for path, _ in api.manifest_rows(api.read_entry(entry), "DataPin"))
    with pytest.raises(api.name("Unpinnable"), match="may start only"):
        m.python(1.0)


# ── qa's second review (2026-10-04, ``scratch/reference_architecture/p3/qa2/``) ─────────────────────────


def test_r1_the_key_is_exact_all_the_way_up(make):
    """Second review, finding 1: the key was exact at its leaves only, so a masked array and its plain twin,
    two dict orders, a sparse matrix's explicit zero and its format, a ``NamedTuple`` and its tuple, and a
    ``compare=False`` constructor field each shared one key and served the other call's answer. The key is
    the exact encoding of what the generating process RECEIVES (each object's constructor form)."""
    package = make(m='''
        import dataclasses
        from typing import NamedTuple
        import numpy as np
        from orpheus.numerics.traced_memo import traced_memo

        class P(NamedTuple):
            a: float
            b: float

        @dataclasses.dataclass(frozen=True)
        class Spec:
            n: int
            label: str = dataclasses.field(default="a", compare=False)

        @traced_memo
        def anything(x) -> str:
            return type(x).__name__
    ''')
    import scipy.sparse as sp

    from orpheus.numerics.content import ContentlessError, FrozenMapping

    m = package.module()
    key = m.anything.key
    with pytest.raises(ContentlessError):  # its state holds a type, which has no content
        key(np.ma.masked_array(np.array([1.0, 2.0]), mask=[0, 1]))
    from orpheus.data.emission_spectrum import EmissionSpectrum

    assert key(EmissionSpectrum(np.array([0.6, 0.4]))) != key(np.array([0.6, 0.4]))  # a stateless subclass: its type
    assert key({"fuel": 1, "mod": 2}) != key({"mod": 2, "fuel": 1})
    assert key(FrozenMapping([("fuel", 1), ("mod", 2)])) != key(FrozenMapping([("mod", 2), ("fuel", 1)]))
    plain = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 0.0]]))
    stored_zero = sp.csr_matrix((np.array([1.0, 0.0]), np.array([0, 1]), np.array([0, 1, 2])), shape=(2, 2))
    assert key(plain) != key(stored_zero)
    assert key(plain) != key(plain.tocsc())
    assert key((1.0, 2.0)) != key(m.P(1.0, 2.0))
    assert key(m.Spec(1, "a")) != key(m.Spec(1, "b"))
    assert m.anything(m.P(1.0, 2.0)) == "P" and m.anything((1.0, 2.0)) == "tuple"


@pytest.mark.parametrize("starter", ["pool", "shell", "cat"])
def test_r2_r3_only_a_declared_program_may_be_started(make, starter):
    """Second review, findings 2 and 3: a spawn-context process pool ran Python nothing recorded, and a
    started program was pinned by its launcher's bytes while the files it read went unpinned (``cat`` served
    1.0 against 50.0; ``shell=True`` and ``env python3`` slipped past the Python refusal). A generation may
    start only a declared machine query (``uname``); every other start is refused."""
    package = make(m='''
        import multiprocessing, subprocess
        from orpheus.numerics.traced_memo import traced_memo

        def _double(x):
            return 2.0 * x

        @traced_memo
        def pool(a: float) -> float:
            with multiprocessing.get_context("spawn").Pool(1) as workers:
                return workers.map(_double, [a])[0]

        @traced_memo
        def shell(a: float) -> float:
            return a + float(subprocess.run("echo 1", shell=True, capture_output=True, text=True).stdout)

        @traced_memo
        def cat(a: float) -> float:
            return a + float(subprocess.run(["cat", __file__], capture_output=True, text=True).stdout.count("def"))
    ''')
    memo = getattr(package.module(), starter)
    with pytest.raises(api.name("Unpinnable"), match="may start only"):
        memo(1.0)
    assert _kind(memo.lookup(1.0)) == "Absent"


def test_r4_a_listing_made_by_a_module_body_during_its_import_is_pinned(make):
    """Second review, finding 4: a directory listed by a module body while it was imported was exempt as the
    import system's own (an import was on the stack), and served 11.0 against 22.0. Only the import system's
    OWN calls are exempt."""
    package = make(m='''
        import os
        from pathlib import Path
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def gen(a: float) -> float:
            from . import lister
            return a + lister.COUNT
    ''', lister='''
        import os
        from pathlib import Path
        COUNT = 10.0 * len(os.listdir(Path(__file__).parent / "tables"))
    ''')
    (package.dir / "tables").mkdir()
    (package.dir / "tables" / "a.dat").write_text("1")
    m = package.module()
    assert m.gen(1.0) == 11.0
    (package.dir / "tables" / "b.dat").write_text("1")
    assert _kind(m.gen.lookup(1.0)) == "Stale" and m.gen(1.0) == 21.0


def test_r5_a_python_file_read_as_data_is_pinned(make):
    """Second review, finding 5: every ``.py`` read was skipped as code, so a generator reading a ``.py`` file
    as data served 3.0 against 30.0. Only the import system's own reads are exempt."""
    package = make(m='''
        from pathlib import Path
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def gen(a: float) -> float:
            return a * float((Path(__file__).parent / "table.py").read_text().split("=")[1])
    ''')
    (package.dir / "table.py").write_text("VALUE = 3")
    m = package.module()
    assert m.gen(1.0) == 3.0
    (package.dir / "table.py").write_text("VALUE = 30")
    assert _kind(m.gen.lookup(1.0)) == "Stale" and m.gen(1.0) == 30.0


def test_r7_a_run_that_changes_directory_then_reads_a_relative_path_is_a_hit(make, tmp_path):
    """Second review, finding 7: the working-directory pin held the run's FINAL directory, so a generator that
    changed directory and read a relative path was stale on every call. It holds the directory the run
    started in."""
    package = make(m='''
        import os
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def gen(a: float, where: str) -> float:
            os.chdir(where)
            with open("table.txt") as handle:
                return a + float(handle.read())
    ''')
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "table.txt").write_text("4")
    m = package.module()
    assert m.gen(1.0, str(tmp_path / "data")) == 5.0
    assert _kind(m.gen.lookup(1.0, str(tmp_path / "data"))) == "Hit"


def test_q5_absence_probed_by_exists_is_a_dependency(make):
    """Finding 5, the half the first fix declared blind: an answer that depended on ``exists()`` finding no
    file was served stale (1.0 against 99.0). The recorder wraps ``os.stat`` while a run is recorded, so a
    probe's finding is pinned."""
    package = make(m='''
        from pathlib import Path
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def override(a: float) -> float:
            path = Path(__file__).parent / "override.txt"
            return float(path.read_text()) if path.exists() else a
    ''')
    m = package.module()
    assert m.override(1.0) == 1.0
    (package.dir / "override.txt").write_text("99")
    assert _kind(m.override.lookup(1.0)) == "Stale" and m.override(1.0) == 99.0


def test_q8_the_memo_pin_covers_the_write_path():
    """Finding 8, its second half: the payload's constructor form (``content.constructor_arguments``) runs
    on the write path after the recording stops, so the memo's pin covers ``content.py`` too."""
    names = {path.name for path in api.name("_MEMO_SOURCES")}
    assert names == {"traced_memo.py", "_traced_memo_boot.py", "content.py"}, names


def test_q10_a_reader_never_sees_a_half_written_entry(tmp_path):
    """Finding 10's reader leg (the second review: the writer-only gate was blind to an in-place write): while
    four threads rewrite one key, every lookup is a ``Hit`` or ``Absent``, never ``Corrupt``."""
    store = api.name("Store")(tmp_path / "cache")
    manifest = api.name("Manifest")(frozenset())
    stop = threading.Event()
    seen: list[str] = []

    def writer() -> None:
        while not stop.is_set():
            store.write("f:x", "key", manifest, (np.arange(20000.0),), set())

    threads = [threading.Thread(target=writer) for _ in range(4)]
    for thread in threads:
        thread.start()
    try:
        for _ in range(300):
            seen.append(_kind(store.lookup("f:x", "key")))
    finally:
        stop.set()
        for thread in threads:
            thread.join()
    assert "Corrupt" not in seen and "Hit" in seen, sorted(set(seen))


@pytest.mark.timeout(600)
def test_q13_a_chain_of_generations_through_distinct_keys_is_bounded(make):
    """Finding 13's bound (the second review: removing the cycle guard hung instead of failing): a chain
    through DISTINCT keys (``f(n)`` calling ``f(n + 1)``) has no repeated key for the ancestry to catch, so
    the depth of nested generations is bounded."""
    package = make(m='''
        from orpheus.numerics.traced_memo import traced_memo

        @traced_memo
        def f(n: int) -> int:
            return f(n + 1)
    ''')
    with pytest.raises(RecursionError, match="deep"):
        package.module().f(0)
