"""The fixture tree of the cone store's gates, and the helpers that drive ``tools.cone_store`` on it.

Spec: ``scratch/verification_tree/cone_store/spec_b1.md`` (§1 the surface, §2 the tree). The expected sets
in the gates are read off the sources below BY HAND (the spec's §5): no gate derives its expectation from a
recorder. A gate that edits the tree edits its own copy (:func:`copy_captured`), never a shared one.

The command runs in a subprocess (``python -O -m tools.cone_store``), so ``PYTHONPATH`` carries the fixture
root (for ``pkg`` and ``tests``) and this repository (for ``tools`` and ``orpheus``). Every ``select`` that a
gate reads runs in a subprocess too (:func:`select`): the memo's ``ChildPin`` validation resolves its pins
through the validating process's ``sys.path``, and the memo's ``Kept`` reads are process-global, so an
in-process ``select`` in the gate process would read a stale or foreign state.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]

#: How long a written file must age before its identity is trusted (the memo's ``TRUSTED_AGE_NS``, 2 s),
#: plus a margin. Read from the memo, never restated: a change there moves this.
def _trusted_age_s() -> float:
    from orpheus.numerics._traced_memo_boot import TRUSTED_AGE_NS

    return TRUSTED_AGE_NS / 1e9 + 0.2


#: The files of the base tree. Each value is dedented source.
_BASE: dict[str, str] = {
    "pyproject.toml": """
        [tool.pytest.ini_options]
        addopts = "-p no:cacheprovider"
        testpaths = ["tests"]
    """,
    "pkg/__init__.py": "",
    "pkg/consts.py": "TOL = 1e-15\n",
    "pkg/conf_dep.py": """
        def conftest_helper():
            return 3 if FLAG else 0


        FLAG = True
    """,
    "pkg/user.py": """
        import functools
        from pkg.consts import TOL


        def _build(n):
            return [i * TOL for i in range(n)]


        @functools.cache
        def table(n):
            return _build(n)


        def above(x):
            return x > TOL


        def untouched():
            return 0
    """,
    "pkg/hidden.py": """
        import functools


        def _fill_body(n):
            return list(range(n))


        @functools.cache
        def cached_f(n):
            return _fill_body(n)
    """,
    "pkg/twins.py": """
        if False:
            def pick():
                return 2
        else:
            def pick():
                return 1


        def outer():
            return (lambda y: y + 1)(1)
    """,
    "tests/__init__.py": "",
    "tests/conftest.py": """
        from pkg.conf_dep import conftest_helper

        conftest_helper()
    """,
    "tests/_helpers.py": """
        def helper():
            return 7
    """,
    "tests/data/ref.txt": "0.5\n",
    "tests/test_a.py": """
        import pytest
        from pkg.user import table


        @pytest.mark.parametrize("t", [table(3)])
        def test_a(t):
            assert len(t) == 3
    """,
    "tests/test_b.py": """
        from pkg import user


        def test_b():
            assert len(user.table(3)) == 3
            assert user.above(1.0)
    """,
    "tests/test_c.py": """
        from pkg import hidden


        def test_c1():
            assert hidden.cached_f(4) == [0, 1, 2, 3]


        def test_c2():
            assert hidden.cached_f(4) == [0, 1, 2, 3]


        def test_c3():
            # Run after test_c1 in ONE process this is red; in a child forked from a parent that never ran
            # test_c1 it is green (spec T4: state leaking between tests of one module).
            assert hidden.cached_f.cache_info().currsize == 0
    """,
    "tests/test_d.py": """
        def test_d():
            from tests._helpers import helper

            assert helper() == 7


        def test_plain():
            assert 1 + 1 == 2
    """,
    "tests/test_h.py": """
        from pathlib import Path


        def test_reads_data():
            assert Path(__file__).with_name("data").joinpath("ref.txt").read_text() == "0.5\\n"


        def test_tmp(tmp_path):
            (tmp_path / "x").write_text("1")
            assert (tmp_path / "x").read_text() == "1"
    """,
    "tests/test_f.py": """
        from pkg import twins


        def test_twins():
            assert twins.pick() == 1
            assert twins.outer() == 2
    """,
    "tests/test_g.py": """
        from pkg import user


        def test_fail():
            assert user.above(0.0)
    """,
    "tests/test_m.py": """
        from pkg import user
        from orpheus.numerics.traced_memo import DefPin, trace_call


        def test_nested():
            value, manifest = trace_call(user.above, 1.0)
            assert value is True
            defs = {(p.path.replace("\\\\", "/"), p.qualname) for p in manifest.pins if isinstance(p, DefPin)}
            assert any(path.endswith("pkg/user.py") and name == "above" for path, name in defs)
            assert not any("tools/cone_store" in path for path, _ in defs)
    """,
}

#: The memo's part (spec §2, test_e): a stored reference, its validation test and a reader.
_MEMO: dict[str, str] = {
    "pkg/oracle.py": """
        from orpheus.numerics.traced_memo import traced_memo


        def _oracle_body(n):
            return float(n) * 0.5


        @traced_memo
        def oracle(n: int) -> float:
            return _oracle_body(n)
    """,
    "tests/test_e.py": """
        from orpheus.numerics.traced_memo import bypass
        from pkg.oracle import oracle


        def test_reader():
            assert oracle(2) == 1.0


        def test_oracle_validation():
            with bypass():
                assert oracle(2) == 1.0
    """,
}

#: A test that reads a memo entry from a store of its own under ``tmp_path`` (the traced memo's own gates do), whose
#: generation imports a module the test's process never imports (B3 acceptance X1.3, 2026-11-03: 70 such readers were
#: selected on every ``select``, their ``ChildPin`` judged against the default store). The second test also runs the
#: generated ``_deep_body`` in its own process, so its own cone holds that def although the read's entry holds it too.
_OWN_STORE: dict[str, str] = {
    "pkg/deep.py": """
        def _deep_body(n):
            return float(n) * 0.25
    """,
    "pkg/own_oracle.py": """
        from orpheus.numerics.traced_memo import traced_memo


        @traced_memo
        def own_oracle(n: int) -> float:
            from pkg.deep import _deep_body

            return _deep_body(n)
    """,
    "tests/test_own.py": """
        from orpheus.numerics.traced_memo import cache_root
        from pkg.own_oracle import own_oracle


        def test_own_store(tmp_path):
            with cache_root(tmp_path / "own"):
                assert own_oracle(4) == 1.0


        def test_own_store_and_its_generator(tmp_path):
            from pkg.deep import _deep_body

            with cache_root(tmp_path / "own"):
                assert own_oracle(4) == _deep_body(4)
    """,
}

_SPAWNERS: dict[str, str] = {
    "tests/test_spawn.py": """
        import multiprocessing
        import os
        import subprocess
        import sys


        def test_subprocess_run():
            subprocess.run([sys.executable, "-c", "pass"], check=True)


        def test_os_system():
            assert os.system("true") == 0


        def test_multiprocessing():
            process = multiprocessing.get_context("spawn").Process(target=os.getpid)
            process.start()
            process.join()


        def test_os_fork():
            pid = os.fork()
            if pid == 0:
                os._exit(0)
            os.waitpid(pid, 0)


        def test_posix_spawn():
            pid = os.posix_spawn("/usr/bin/true", ["true"], dict(os.environ))
            os.waitpid(pid, 0)
    """,
    "tests/test_spawn_c.py": """
        import subprocess
        import sys

        subprocess.run([sys.executable, "-c", "pass"], check=True)


        def test_c_one():
            pass


        def test_c_two():
            pass
    """,
}

#: A test that reads a first-party source file as TEXT while its module also ran (spec S11, the battery's
#: finding on the store's drop of a DataPin for a ``.py`` file that is also a pinned module).
_READS_SOURCE: dict[str, str] = {
    "tests/test_k.py": """
        from pathlib import Path

        import pkg.user


        def test_reads_source():
            assert "The table." not in Path(pkg.user.__file__).read_text()
    """,
}

#: Review round 1 (qa's under-selections U3, U7, U9): a test sub-directory, a test that reads a repository file
#: and rewrites it, and a test that holds still while the gate edits code it ran.
_SUBDIR: dict[str, str] = {
    "tests/sub/__init__.py": "",
    "tests/sub/test_x.py": """
        def test_x():
            assert True
    """,
}
_READ_REWRITE: dict[str, str] = {
    "tests/data/counter.txt": "0\n",
    "tests/test_rw.py": """
        from pathlib import Path


        def test_rw():
            counter = Path(__file__).with_name("data") / "counter.txt"
            value = int(counter.read_text())
            assert value < 100
            counter.write_text(f"{value + 1}\\n")
    """,
}
_RACE: dict[str, str] = {
    "tests/test_race.py": """
        import os
        import time
        from pathlib import Path

        from pkg import user


        def test_race():
            assert user.above(1.0)
            Path(os.environ["CONE_RACE_MARK"]).write_text("ran")
            deadline = time.monotonic() + 30
            while not Path(os.environ["CONE_RACE_MARK"] + ".edited").exists() and time.monotonic() < deadline:
                time.sleep(0.05)
    """,
}

#: qa's second pass (D2, D3, D4/D5): an off-path conftest with a session hook, a user ``pytest_ignore_collect``
#: that fills a cache a test reads again, and a local distribution that can grow a ``pytest11`` entry point.
_OFFPATH: dict[str, str] = {
    "tests/zother/__init__.py": "",
    "tests/zother/conftest.py": """
        def pytest_configure(config):
            pass
    """,
    "tests/zother/test_z.py": """
        def test_z():
            assert True
    """,
}
_USER_IGNORE: dict[str, str] = {
    "tests/data/flag.txt": "ok\n",
    "tests/flagread.py": """
        import functools
        from pathlib import Path


        @functools.cache
        def flag():
            return (Path(__file__).parent / "data" / "flag.txt").read_text()
    """,
    "tests/test_flag.py": """
        from tests.flagread import flag


        def test_flag():
            assert flag() == "ok\\n"
    """,
}
_USER_IGNORE_CONFTEST = """

def pytest_ignore_collect(collection_path, config):
    from tests.flagread import flag

    flag()
    return None
"""
_LOCALPLUG: dict[str, str] = {
    "localplug_mod.py": """
        def pytest_configure(config):
            pass
    """,
    "localplug-1.0.dist-info/METADATA": """
        Metadata-Version: 2.1
        Name: localplug
        Version: 1.0
    """,
    "localplug-1.0.dist-info/RECORD": "",
}

#: B3 (spec ``scratch/verification_tree/cone_store/b3/spec_b3.md`` §2.2): a test that calls, in its own run, the def
#: the root conftest already ran at import (``conftest_helper``). Its OWN cone must hold that def; ``test_b``'s must not.
_PHASES: dict[str, str] = {
    "tests/test_p.py": """
        from pkg.conf_dep import conftest_helper


        def test_p():
            assert conftest_helper() == 3
    """,
}

#: B3 placement (spec §2.3-§2.5, rows PL1-PL11): a layer ``pkg/lay`` ranked by :data:`PLACEMENT_ORDER`, with one
#: import of each kind, an upstream layer ``pkg/up``, a higher one ``pkg/top``, and a conftest that runs ``alpha.fd`` (so ``fd``, ``fb``, ``fa``
#: are in every test's inherited floor). Module names are chosen so that rank order, alphabetical order and def
#: count disagree: the top of the order is ``alpha``, the bottom ``zeta``, and ``zeta`` holds two defs.
PLACEMENT_ORDER: tuple[tuple[str, ...], ...] = (
    ("zeta", "s"), ("beta", "beta2"), ("gamma",), ("alpha",))
_PLACEMENT: dict[str, str] = {
    "pkg/lay/__init__.py": "",
    "pkg/lay/zeta.py": """
        def fa():
            return _fa_inner()


        def _fa_inner():
            return 1
    """,
    "pkg/lay/s.py": "".join(f"def f{i}():\n    return {i}\n\n\n" for i in range(20))
                    + "ALL = (" + ", ".join(f"f{i}" for i in range(20)) + ")\n",
    "pkg/lay/beta.py": """
        from pkg.lay import s, zeta


        def fb():
            return zeta.fa() + 1
    """,
    "pkg/lay/beta2.py": """
        def fb2():
            return 2
    """,
    "pkg/lay/gamma.py": """
        from typing import TYPE_CHECKING

        if TYPE_CHECKING:
            from pkg.lay.beta import fb


        def fc():
            from pkg.lay.zeta import fa

            return fa()
    """,
    "pkg/lay/alpha.py": """
        from pkg.lay import beta


        def fd():
            return beta.fb() + 1
    """,
    "pkg/up/__init__.py": "",
    "pkg/up/u.py": """
        def u():
            return 0
    """,
    "pkg/top/__init__.py": "",
    "pkg/top/t.py": """
        def t():
            return 0
    """,
    "tests/lay/__init__.py": "",
    "tests/lay/conftest.py": """
        from pkg.lay.alpha import fd

        fd()
    """,
    "tests/lay/test_pl.py": """
        import subprocess
        import sys

        import pytest

        from pkg.lay import alpha, beta, beta2, gamma, s, zeta
        from pkg.top import t
        from pkg.up import u


        def test_low():
            assert zeta.fa() == 1


        def test_d_only():
            assert alpha.fd() == 3


        def test_none():
            assert alpha.fd() + beta2.fb2() == 5


        def test_lazy_typ():
            assert gamma.fc() + beta.fb() == 3


        def test_tie():
            assert beta.fb() + beta2.fb2() == 4


        def test_layers():
            assert zeta.fa() + u.u() + t.t() == 1


        def test_nodef():
            assert True


        def test_fails():
            assert zeta.fa() == 2


        def test_spawn():
            assert zeta.fa() == 1
            subprocess.run([sys.executable, "-c", "pass"], check=True)


        @pytest.mark.parametrize("k", [1, 2])
        def test_s_lower(k):
            assert sum(f() for f in s.ALL) == 190


        def test_s_upper95():
            assert sum(f() for f in s.ALL[:19]) + beta.fb() == 173


        def test_s_upper90():
            assert sum(f() for f in s.ALL[:18]) == 153


        def test_s_super():
            assert sum(f() for f in s.ALL) + alpha.fd() == 193


        def test_tiny():
            assert beta2.fb2() == 2


        def test_dup1():
            assert sum(f() for f in s.ALL[:6]) + zeta.fa() == 16


        def test_dup2():
            assert sum(f() for f in s.ALL[:6]) + zeta.fa() > 0
    """,
}

#: B5 reads (spec ``scratch/verification_tree/cone_store/b5/spec_b5.md`` §3.0): a layer ``pkg/rd`` ranked by
#: :data:`READS_ORDER`, one test per shape of read the resolution must name (or must not), an apparatus helper
#: ``tests/rd/_helper.py``, a conftest whose pytest hook reads another package (``pkg/other``) and whose fixture builds
#: a ``low.Box``, and source audits of one file, of the layer and of the package. The fixture's tests raise
#: ``AssertionError`` by hand rather than ``assert``, except ``test_getattr``: the assertion rewriter adds loads of
#: its own (a bare module load per operand), which only R11 is about. No fixture test carries a docstring (its string
#: constants would be paired with the modules it reads, C3).
READS_ORDER: tuple[tuple[str, ...], ...] = (("low",), ("mid", "rex"), ("high",))
_READS: dict[str, str] = {
    "pkg/rd/__init__.py": """
        from pkg.rd.low import Rec
    """,
    "pkg/rd/low.py": """
        from dataclasses import dataclass
        from typing import Protocol, runtime_checkable

        TOL = 1e-12
        ALIAS = TOL


        def _wrap(function):
            def wrapped():
                return function()

            return wrapped


        def f():
            return 1


        f = _wrap(f)


        @dataclass(frozen=True)
        class Rec:
            value: float


        @runtime_checkable
        class Proto(Protocol):
            def apply(self): ...


        class Box:
            LIMIT = 3


        from pkg.other.alt import REBOUND

        REBOUND = 5

        globals()["LATE2"] = 2

        try:
            from pkg.rd._fast import SPEED
        except ImportError:
            SPEED = 0
    """,
    "pkg/rd/mid.py": """
        from pkg.rd import low

        BUILT = low.Rec(low.TOL)
        FROM_LATE = low.LATE2


        def g():
            return low.f()
    """,
    "pkg/rd/rex.py": """
        from pkg.rd.low import TOL
    """,
    "pkg/rd/high.py": """
        HIGH = 7


        def h():
            return 1
    """,
    "pkg/other/__init__.py": "",
    "pkg/other/o.py": "OTHER = 1\n",
    "pkg/other/alt.py": "REBOUND = 99\n",
    "tests/rd/__init__.py": "",
    "tests/rd/_helper.py": """
        from pkg.rd import high
        from pkg.rd.mid import BUILT

        TABLE = (BUILT,)


        def make():
            return high.HIGH
    """,
    "tests/rd/conftest.py": """
        import pytest

        import pkg.other.o as other
        from pkg.rd import low


        def pytest_runtest_setup(item):
            if other.OTHER != 1:
                raise RuntimeError("the hook's read")


        @pytest.fixture
        def box():
            return low.Box()
    """,
    "tests/rd/test_rd.py": """
        import pkg.rd.low
        import pytest

        import pkg.other.o as other_module
        from pkg.rd import Rec, high, low, mid
        from pkg.rd.rex import TOL
        from tests.rd import _helper
        from tests.rd._helper import BUILT

        _G = (1.0, 20.0)
        CASES = (low.Rec(1.0),)


        def test_global():
            if not low.TOL > 0:
                raise AssertionError


        def test_chain():
            if not pkg.rd.low.TOL > 0:
                raise AssertionError


        def test_in_def():
            from pkg.rd import low as m

            if not m.TOL > 0:
                raise AssertionError


        def test_cell():
            from pkg.rd import low as m

            def inner():
                return m.TOL

            if not inner() > 0:
                raise AssertionError


        def test_reexport():
            if not TOL > 0:
                raise AssertionError


        def test_package_reexport():
            if Rec(1.0).value != 1.0:
                raise AssertionError


        def test_fallback():
            if low.SPEED != 0:
                raise AssertionError


        def test_rebound():
            if low.REBOUND != 5:
                raise AssertionError


        def test_superinstruction():
            import pkg as pk; value = pk.rd.high.h()  # noqa: E702
            if value != 1:
                raise AssertionError


        def test_closure_unbound():
            if mid.FROM_LATE + high.h() != 3:
                raise AssertionError


        def test_two_cells():
            def first():
                from pkg.rd import low as m

                def inner():
                    return m.TOL

                return inner()

            def second():
                from pkg.rd import high as m

                def inner():
                    return m.HIGH

                return inner()

            if second() + high.h() != 8 or first() <= 0:
                raise AssertionError


        def test_fast_pair():
            from pkg.rd import low as m
            x = 1.0
            if max(x, m.TOL) != x or low.f() != 1:
                raise AssertionError


        def test_def_and_assign():
            if low.f() != 1:
                raise AssertionError


        def test_dataclass():
            if low.Rec(1.0) != low.Rec(1.0):
                raise AssertionError


        def test_protocol():
            if isinstance(object(), low.Proto):
                raise AssertionError


        def test_getattr():
            assert getattr(low, "TOL") > 0
            assert hasattr(low, "ALIAS")


        def test_absent():
            if hasattr(low, "GONE"):
                raise AssertionError


        def test_dunder():
            if not low.__file__.endswith(".py"):
                raise AssertionError


        def test_instance_attr(box):
            if box.LIMIT != 3:
                raise AssertionError


        def test_helper_reexport():
            if BUILT.value <= 0:
                raise AssertionError


        class TestParam:
            @pytest.mark.parametrize("tol", [low.TOL])
            def test_parametrized(self, tol):
                if tol <= 0:
                    raise AssertionError


        def test_premise():
            if not max(_G) / min(_G) > 10:
                raise AssertionError


        def test_case_table():
            if len(CASES) != 1:
                raise AssertionError


        def test_helper_table():
            if len(_helper.TABLE) != 1:
                raise AssertionError


        def test_legendre_like():
            if mid.BUILT.value <= 0:
                raise AssertionError


        def test_helper_def():
            if _helper.make() != 7:
                raise AssertionError


        def test_uses_premise():
            if low.f() * max(_G) <= 0:
                raise AssertionError


        def test_input_value():
            if low.f() + high.HIGH != 8:
                raise AssertionError


        def test_two_reads():
            if mid.BUILT.value != low.TOL:
                raise AssertionError


        def test_two_reads_apart():
            if mid.BUILT.value + high.HIGH <= 0:
                raise AssertionError


        def test_other():
            if other_module.OTHER != 1:
                raise AssertionError


        @pytest.fixture
        def late():
            setattr(low, "LATE", 1)


        def test_unresolved(late):
            if low.LATE != 1:
                raise AssertionError
    """,
    "tests/rd/test_star.py": """
        from pkg.rd.low import *  # noqa: F403


        def test_star():
            if not TOL > 0:  # noqa: F405
                raise AssertionError
    """,
    "tests/rd/test_audit.py": """
        from pathlib import Path

        from pkg.rd import mid


        def test_audit_one():
            if "HIGH" not in Path("pkg/rd/high.py").read_text():
                raise AssertionError


        def test_audit_layer():
            if len([p.read_text() for p in sorted(Path("pkg/rd").glob("*.py"))]) != 5:
                raise AssertionError


        def test_audit_package():
            if len([p.read_text() for p in sorted(Path("pkg").rglob("*.py"))]) < 10:
                raise AssertionError


        def test_audit_run():
            if mid.g() != 1 or len([p.read_text() for p in sorted(Path("pkg/rd").glob("*.py"))]) != 5:
                raise AssertionError
    """,
}


def _defs(prefix: str, n: int) -> str:
    """``n`` one-line defs ``<prefix>0`` ... returning their index (a fixture module's source)."""
    return "".join(f"def {prefix}{i}():\n    return {i}\n\n\n" for i in range(n))


def _base(prefix: str) -> str:
    """A base module of 40 defs and their tuple ``ALL``."""
    return _defs(prefix, 40) + "ALL = (" + ", ".join(f"{prefix}{i}" for i in range(40)) + ")\n"


#: Ruling 50 (a support's lower SUBJECT, the function's own defs at its home nodes, must be contained too): a layer
#: ``pkg/sj`` ranked :data:`SUBJECT_ORDER`, with one base of 40 defs per scenario so no two scenarios share a def.
#: (A, base ``a``) ``test_lower`` runs the 40 and ``mid.m0`` (41 defs, home ``mid``, subject ``{m0}``);
#: ``test_upper_misses`` runs the 40 and ``top.t0``-``t9`` (50; cone 40/41, subject 0); ``test_upper_runs`` runs the 40,
#: ``m0`` and ``top.u0``-``u4`` (46; 41/41, 1/1). (B, base ``b``, a TIE) ``test_tie_lower`` runs the 40, ``mid.p0`` and
#: ``mid2.n0`` (42, tied at ``mid`` and ``mid2``, subject ``{p0, n0}``); ``test_tie_partial`` the 40, ``p0`` and
#: ``top.w0``-``w4`` (46; 41/42, 1/2); ``test_tie_full`` the 40, ``p0``, ``n0`` and ``top.x0``-``x4`` (47; 42/42, 2/2).
#: (C, base ``c``, AGAINST the order) ``test_against_lower`` runs the 40 and ``mid.q0`` (home ``mid``);
#: ``test_against_upper`` the 40, ``q0`` and ``side.s0``-``s4`` (home ``side``, which imports nothing: 41/41, 1/1).
#: (D, base ``f``, an EDGE) ``test_edge_lower`` runs the 40 and ``mid.r0``-``r9`` (50, subject the ten ``r``);
#: ``test_edge_upper`` the 40, ``r0``-``r8`` and ``top.y0`` (cone 49/50 = 0.98, subject 9/10 = 0.9 exactly).
#: (E, base ``g``, TWO HOMES) ``test_two_homes[1]`` runs the 40 and ``mid.v0`` (home ``mid``); ``test_two_homes[2]``
#: the 40, ``mid.v1`` and ``top.h0`` (home ``top``): the function's homes are ``mid`` and ``top`` and its subject is
#: ``{v0, v1, h0}`` (case 2's ``v1`` is at the function's home ``mid``, not at the case's own); ``test_two_homes_upper``
#: runs the 40, ``v0``, ``h0`` and ``top.k0``-``k4`` (cone 42/43, subject 2/3; the per-case reading, ``{v0, h0}``,
#: would give 1). Every other pair holding at 0.8 has a subject of 0 (the largest such cone, 49/50, is
#: ``test_edge_lower`` over ``test_edge_upper``). ``test_reads_a_constant`` runs no def and reads ``mid.CONST``;
#: ``test_reads_and_runs_above`` reads it and runs ``pkg/sjx/o.o0``, a def ABOVE the layer: placed by statement, with
#: a cone of one def.
SUBJECT_ORDER: tuple[tuple[str, ...], ...] = (("a", "b", "c", "f", "g"), ("mid", "mid2"), ("side", "top"))
_SUBJECT: dict[str, str] = {
    "pkg/sj/__init__.py": "",
    "pkg/sj/a.py": _base("a"),
    "pkg/sj/b.py": _base("b"),
    "pkg/sj/c.py": _base("c"),
    "pkg/sj/f.py": _base("f"),
    "pkg/sj/g.py": _base("g"),
    "pkg/sj/side.py": _defs("s", 5),
    "pkg/sjx/__init__.py": "",
    "pkg/sjx/o.py": _defs("o", 1),
    "pkg/sj/mid.py": textwrap.dedent("""
        from pkg.sj import a, b

        CONST = 3


        def m0():
            return 0


        def p0():
            return 0
    """) + "\n\n" + "".join(_defs(prefix, n) for prefix, n in (("q", 1), ("r", 10), ("v", 2))),
    "pkg/sj/mid2.py": """
        from pkg.sj import b


        def n0():
            return 0
    """,
    "pkg/sj/top.py": "from pkg.sj import mid, mid2\n\n\n"
                      + "".join(_defs(prefix, n) for prefix, n in (("t", 10), ("u", 5), ("w", 5), ("x", 5), ("y", 1),
                                                                 ("h", 1), ("k", 5))),
    "tests/sj/__init__.py": "",
    "tests/sj/test_sj.py": """
        import pytest

        from pkg.sj import a, b, c, f, g, mid, mid2, side, top
        from pkg.sjx import o


        def test_lower():
            if sum(f() for f in a.ALL) + mid.m0() != 780:
                raise AssertionError


        def test_upper_misses():
            if sum(f() for f in a.ALL) + sum(getattr(top, f"t{i}")() for i in range(10)) != 825:
                raise AssertionError


        def test_upper_runs():
            if sum(f() for f in a.ALL) + mid.m0() + sum(getattr(top, f"u{i}")() for i in range(5)) != 790:
                raise AssertionError


        def test_tie_lower():
            if sum(f() for f in b.ALL) + mid.p0() + mid2.n0() != 780:
                raise AssertionError


        def test_tie_partial():
            if sum(f() for f in b.ALL) + mid.p0() + sum(getattr(top, f"w{i}")() for i in range(5)) != 790:
                raise AssertionError


        def test_tie_full():
            if sum(f() for f in b.ALL) + mid.p0() + mid2.n0() + sum(getattr(top, f"x{i}")() for i in range(5)) != 790:
                raise AssertionError


        def test_reads_a_constant():
            if mid.CONST != 3:
                raise AssertionError


        def test_reads_and_runs_above():
            if mid.CONST + o.o0() != 3:
                raise AssertionError


        def test_against_lower():
            if sum(fn() for fn in c.ALL) + mid.q0() != 780:
                raise AssertionError


        def test_against_upper():
            if sum(fn() for fn in c.ALL) + mid.q0() + sum(getattr(side, f"s{i}")() for i in range(5)) != 790:
                raise AssertionError


        def test_edge_lower():
            if sum(fn() for fn in f.ALL) + sum(getattr(mid, f"r{i}")() for i in range(10)) != 825:
                raise AssertionError


        def test_edge_upper():
            if sum(fn() for fn in f.ALL) + sum(getattr(mid, f"r{i}")() for i in range(9)) + top.y0() != 816:
                raise AssertionError


        @pytest.mark.parametrize("case", [1, 2])
        def test_two_homes(case):
            total = sum(fn() for fn in g.ALL) + (mid.v0() if case == 1 else mid.v1() + top.h0())
            if total != 780 + (case - 1):
                raise AssertionError


        def test_two_homes_upper():
            if sum(fn() for fn in g.ALL) + mid.v0() + top.h0() + sum(getattr(top, f"k{i}")() for i in range(5)) != 790:
                raise AssertionError
    """,
}

_THREAD_LINES = """
import threading as _threading
import time as _time

_threading.Thread(target=_time.sleep, args=(3600,), daemon=True).start()
"""

_UNMEASURED_PLATFORM_LINES = """
from orpheus.numerics import traced_memo as _traced_memo

_traced_memo.fork_platform = lambda: ("linux", "x86_64", "openblas")
"""

_THREAD_MODULE = """
from pkg import user
""" + _THREAD_LINES + """

def test_thread():
    assert user.above(1.0)
"""

#: Node ids of the base tree (memo part included), written by hand.
BASE_IDS = frozenset({
    "tests/test_a.py::test_a[t0]",
    "tests/test_b.py::test_b",
    "tests/test_c.py::test_c1",
    "tests/test_c.py::test_c2",
    "tests/test_c.py::test_c3",
    "tests/test_d.py::test_d",
    "tests/test_d.py::test_plain",
    "tests/test_h.py::test_reads_data",
    "tests/test_h.py::test_tmp",
    "tests/test_f.py::test_twins",
    "tests/test_g.py::test_fail",
    "tests/test_m.py::test_nested",
})
MEMO_IDS = frozenset({"tests/test_e.py::test_reader", "tests/test_e.py::test_oracle_validation"})
SPAWNER_IDS = {
    "subprocess.run": "tests/test_spawn.py::test_subprocess_run",
    "os.system": "tests/test_spawn.py::test_os_system",
    "multiprocessing": "tests/test_spawn.py::test_multiprocessing",
    "os.fork": "tests/test_spawn.py::test_os_fork",
    "os.posix_spawn": "tests/test_spawn.py::test_posix_spawn",
}
SPAWN_AT_COLLECTION_IDS = frozenset({"tests/test_spawn_c.py::test_c_one", "tests/test_spawn_c.py::test_c_two"})
#: What ``select`` lists with no edit after a capture of the base tree: the tests that did not pass
#: (spec S1). No base test spawns.
BASELINE = frozenset({"tests/test_g.py::test_fail"})


#: Ruling 53 (a declaration's verdict: ``holds`` when SOME lower under its target holds, a lower with a layer subject
#: by ``c >= threshold`` AND ``c_home >= threshold``, a lower with none by ``c_items >= threshold``): a layer ``pkg/vd``
#: ranked :data:`VERDICT_ORDER`, one base of 40 defs per scenario. Counted by hand off the sources below.
#: (1, base ``a``, the SUBJECT WITHOUT THE CONE) ``test_v_lower`` runs the 40 and ``mid.m0`` (41, home ``mid``, subject
#: ``{m0}``); ``test_v_subject_only`` runs ``a0``-``a29``, ``m0`` and ``top.t0`` (cone 31/41, subject 1/1).
#: (2, base ``b``, SOME LOWER) ``TestPair::test_held`` runs the 40 and ``mid.p0``; ``TestPair::test_missed`` the 40 and
#: ``mid.n0`` (each 41, home ``mid``); ``test_v_some`` runs the 40, ``p0`` and ``top.x0`` (41/41 and 1/1 over the
#: first, 40/41 and 0/1 over the second).
#: (3, a STRUCTURAL lower) ``test_v_reads`` runs no def and reads ``mid.CONST``; ``test_v_reads_too`` reads it and runs
#: ``top.u0`` (it holds every item the lower read); ``test_v_subject_only`` reads nothing of ``mid``.
#: (4, base ``f``, the THRESHOLD) ``test_v_edge_lower`` runs the 40 and ``mid.r0``-``r9`` (50, subject the ten ``r``);
#: ``test_v_edge_upper`` the 40, ``r0``-``r8`` and ``top.y0`` (cone 49/50 = 0.98, subject 9/10 = 0.9 exactly).
#: (5, a target holding the UPPER ALONE) ``TestAlone::test_alone`` runs ``top.t0``; declared on its own class, the
#: target's only function is the upper itself.
VERDICT_ORDER: tuple[tuple[str, ...], ...] = (("a", "b", "f"), ("mid",), ("top",))
_VERDICT: dict[str, str] = {
    "pkg/vd/__init__.py": "",
    "pkg/vd/a.py": _base("a"),
    "pkg/vd/b.py": _base("b"),
    "pkg/vd/f.py": _base("f"),
    "pkg/vd/mid.py": "from pkg.vd import a, b, f\n\nCONST = 3\n\n\n" + "".join(
        _defs(prefix, n) for prefix, n in (("m", 1), ("p", 1), ("n", 1), ("r", 10))),
    "pkg/vd/top.py": "from pkg.vd import mid\n\n\n" + "".join(_defs(prefix, 1) for prefix in ("t", "x", "u", "y")),
    "tests/vd/__init__.py": "",
    "tests/vd/test_vd.py": """
        from pkg.vd import a, b, f, mid, top


        def test_v_lower():
            if sum(fn() for fn in a.ALL) + mid.m0() != 780:
                raise AssertionError


        def test_v_subject_only():
            if sum(fn() for fn in a.ALL[:30]) + mid.m0() + top.t0() != 435:
                raise AssertionError


        class TestPair:
            def test_held(self):
                if sum(fn() for fn in b.ALL) + mid.p0() != 780:
                    raise AssertionError

            def test_missed(self):
                if sum(fn() for fn in b.ALL) + mid.n0() != 780:
                    raise AssertionError


        def test_v_some():
            if sum(fn() for fn in b.ALL) + mid.p0() + top.x0() != 780:
                raise AssertionError


        def test_v_reads():
            if mid.CONST != 3:
                raise AssertionError


        def test_v_reads_too():
            if mid.CONST + top.u0() != 3:
                raise AssertionError


        def test_v_edge_lower():
            if sum(fn() for fn in f.ALL) + sum(getattr(mid, f"r{i}")() for i in range(10)) != 825:
                raise AssertionError


        def test_v_edge_upper():
            if sum(fn() for fn in f.ALL) + sum(getattr(mid, f"r{i}")() for i in range(9)) + top.y0() != 816:
                raise AssertionError


        class TestAlone:
            def test_alone(self):
                if top.t0() != 0:
                    raise AssertionError
    """,
}


def build(root: Path, *, memo: bool = True, spawners: bool = False, thread_at_conftest: bool = False,
          thread_at_module: bool = False, unmeasured_platform: bool = False, reads_source: bool = False,
          subdir: bool = False, read_rewrite: bool = False, race: bool = False, offpath: bool = False,
          user_ignore: bool = False, localplug: bool = False, phases: bool = False, placement: bool = False,
          own_store: bool = False, reads: bool = False, subject: bool = False, verdict: bool = False,
          age: bool = True) -> Path:
    """Write the fixture tree under ``root`` and return it; ``age`` waits until every file's identity is
    settled, so a capture taken next records settled identities (spec F1)."""
    files = dict(_BASE)
    if memo:
        files |= _MEMO
    if spawners:
        files |= _SPAWNERS
    if reads_source:
        files |= _READS_SOURCE
    if subdir:
        files |= _SUBDIR
    if read_rewrite:
        files |= _READ_REWRITE
    if race:
        files |= _RACE
    if offpath:
        files |= _OFFPATH
    if user_ignore:
        files |= _USER_IGNORE
    if localplug:
        files |= _LOCALPLUG
    if phases:
        files |= _PHASES
    if placement:
        files |= _PLACEMENT
    if own_store:
        files |= _OWN_STORE
    if reads:
        files |= _READS
    if subject:
        files |= _SUBJECT
    if verdict:
        files |= _VERDICT
    conftest = textwrap.dedent(files["tests/conftest.py"])
    if thread_at_conftest:
        conftest += _THREAD_LINES
    if unmeasured_platform:
        conftest += _UNMEASURED_PLATFORM_LINES
    if user_ignore:
        conftest += _USER_IGNORE_CONFTEST
    for relative, source in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(conftest if relative == "tests/conftest.py" else textwrap.dedent(source))
    if thread_at_module:
        (root / "tests/test_thread.py").write_text(_THREAD_MODULE)
    if age:
        time.sleep(_trusted_age_s())
    return root


def ids_of(root: Path) -> frozenset[str]:
    """The node ids pytest collects on the tree (an instrument independent of the store)."""
    run = subprocess.run([sys.executable, "-O", "-m", "pytest", "--collect-only", "-q", "tests"], cwd=root,
                         env=environment(root), capture_output=True, text=True, timeout=120)
    return frozenset(line.strip() for line in run.stdout.splitlines() if "::" in line)


def environment(root: Path, **overrides: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("PYTEST_", "COV_"))}
    # The fixture root and this repository first, then the path the gate process received (passed on as the cone
    # store's ``worker_environment`` passes it, so an instrument installed on it reaches every child).
    env["PYTHONPATH"] = os.pathsep.join([str(root), str(REPO), *filter(None, [os.environ.get("PYTHONPATH")])])
    env["OMP_NUM_THREADS"] = "1"
    env.update(overrides)
    return env


def paths(root: Path) -> dict[str, Path]:
    return {"store": root.parent / f"{root.name}.cone.sqlite", "memo": root.parent / f"{root.name}.memo"}


def cli(root: Path, command: str, *arguments: str, env: dict[str, str] | None = None,
        timeout: float = 600) -> subprocess.CompletedProcess[str]:
    """``python -O -m tools.cone_store <command> --repo <root> --store ... --package pkg --memo-root ...``."""
    where = paths(root)
    argv = [sys.executable, "-O", "-m", "tools.cone_store", command, "--repo", str(root), "--store", str(where["store"]),
            "--package", "pkg", "--memo-root", str(where["memo"]), *arguments]
    run = subprocess.run(argv, cwd=root, env=environment(root, **(env or {})), capture_output=True, text=True,
                         timeout=timeout)
    if run.returncode != 0:
        raise AssertionError(f"{' '.join(argv)} exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    return run


def capture(root: Path, *arguments: str, **env: str) -> subprocess.CompletedProcess[str]:
    return cli(root, "capture", *arguments, env=env)


def store_of(root: Path, store: Path | None = None):
    """The store as ``tools.cone_store.open_store`` reads it (spec §1)."""
    import tools.cone_store as cone

    return cone.open_store(store if store is not None else paths(root)["store"])


_SELECT_DRIVER = r'''
import json, os, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
opened = []


def hook(event, args):
    if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        path = Path(os.fsdecode(args[0]))
        try:
            resolved = path.resolve()
        except OSError:
            return
        if resolved.suffix == ".py" and resolved.is_relative_to(root):
            opened.append(resolved.relative_to(root).as_posix())


import tools.cone_store as cone

sys.addaudithook(hook)
memo = Path(sys.argv[3]) if sys.argv[3] else None
selected = cone.select(root, Path(sys.argv[2]), package="pkg", memo_root=memo)
sys.addaudithook(lambda *a: None)
print("CONE-SELECT " + json.dumps({"selected": sorted(selected), "opened": sorted(set(opened))}))
'''


def select(root: Path, **env: str) -> tuple[frozenset[str], frozenset[str]]:
    """``tools.cone_store.select`` called in a fresh process: (the node ids selected, the ``.py`` files under the
    root that the selecting process opened). The audit hook sees the store's own reads, not those of the
    ``--collect-only`` subprocess the selection starts (spec §1)."""
    where = paths(root)
    driver = root.parent / f"{root.name}.select_driver.py"
    driver.write_text(_SELECT_DRIVER)
    run = subprocess.run([sys.executable, "-O", str(driver), str(root), str(where["store"]), str(where["memo"])],
                         cwd=root, env=environment(root, **env), capture_output=True, text=True, timeout=600)
    line = next((l for l in run.stdout.splitlines() if l.startswith("CONE-SELECT ")), None)
    if run.returncode != 0 or line is None:
        raise AssertionError(f"select driver exited {run.returncode}:\n{run.stdout[-3000:]}\n{run.stderr[-3000:]}")
    found = json.loads(line[len("CONE-SELECT "):])
    return frozenset(found["selected"]), frozenset(found["opened"])


def copy_captured(source: Path, target: Path) -> Path:
    """A copy of a captured tree with its store and memo store, for a gate that edits. The copy has new
    inodes and change times, so every file goes to stage two: right for selection rows, wrong for the
    freshness route rows (F1-F6 build their own tree)."""
    shutil.copytree(source, target)
    for kind, path in paths(source).items():
        if path.is_dir():
            shutil.copytree(path, paths(target)[kind])
        elif path.exists():
            shutil.copy2(path, paths(target)[kind])
    return target


def edit(root: Path, relative: str, old: str, new: str) -> None:
    """Replace ``old`` by ``new`` exactly once in a fixture file, refusing when ``old`` is not there once (an
    edit that did not apply would read as a correct empty selection)."""
    path = root / relative
    text = path.read_text()
    if text.count(old) != 1:
        raise AssertionError(f"{relative}: expected one {old!r}, found {text.count(old)}")
    path.write_text(text.replace(old, new))


def pin_key(pin) -> tuple:
    """A pin as (kind, path, qualname, ordinal): the identity of the node, without its digest."""
    kind = type(pin).__name__
    path = getattr(pin, "path", getattr(pin, "function_id", ""))
    return kind, str(path).replace("\\", "/"), getattr(pin, "qualname", ""), getattr(pin, "ordinal", -1)


def code_keys(pins) -> frozenset[tuple]:
    """The ``DefPin`` and ``StatementsPin`` keys of a row set (the statement pin replaced the skeleton's ``ModulePin``
    at B2: one per first-party module that ran, keyed by its path)."""
    return frozenset(pin_key(p) for p in pins if type(p).__name__ in {"DefPin", "StatementsPin"})


IN_SUITE_MONITOR = r'''
"""The in-suite monitor of round 1 (cone_plugin_rearm.py): one process, PY_START, DISABLE, re-armed per phase."""
import json, os, sys
from pathlib import Path

import pytest

ROOT = str(Path(os.environ["CONE_FIXTURE_ROOT"]).resolve()) + os.sep
TOOL = 5
_current = None
_cones = {}


def _on_start(code, offset):
    if not code.co_filename.startswith(ROOT):
        return sys.monitoring.DISABLE
    if _current is not None:
        _cones.setdefault(_current, set()).add((code.co_filename[len(ROOT):], code.co_qualname))
    return sys.monitoring.DISABLE


def pytest_configure(config):
    sys.monitoring.use_tool_id(TOOL, "in-suite")
    sys.monitoring.register_callback(TOOL, sys.monitoring.events.PY_START, _on_start)
    sys.monitoring.set_events(TOOL, sys.monitoring.events.PY_START)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_setup(item):
    global _current
    sys.monitoring.restart_events(); _current = item.nodeid; yield; _current = None


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    global _current
    sys.monitoring.restart_events(); _current = item.nodeid; yield; _current = None


def pytest_unconfigure(config):
    Path(os.environ["CONE_MONITOR_OUT"]).write_text(json.dumps({k: sorted(map(list, v)) for k, v in _cones.items()}))
'''


def in_suite_cones(root: Path, *files: str) -> dict[str, frozenset[tuple[str, str]]]:
    """The cones the in-suite monitor records for ``files`` run in one pytest process."""
    plugin_dir = root.parent / f"{root.name}.monitor"
    plugin_dir.mkdir(exist_ok=True)
    (plugin_dir / "in_suite_monitor.py").write_text(IN_SUITE_MONITOR)
    out = plugin_dir / "cones.json"
    env = environment(root, CONE_FIXTURE_ROOT=str(root), CONE_MONITOR_OUT=str(out))
    env["PYTHONPATH"] = os.pathsep.join([str(plugin_dir), env["PYTHONPATH"]])
    subprocess.run([sys.executable, "-O", "-m", "pytest", "-q", "-p", "in_suite_monitor", *files], cwd=root, env=env,
                   capture_output=True, text=True, timeout=300)
    return {k: frozenset(map(tuple, v)) for k, v in json.loads(out.read_text()).items()}


class Captured:
    """A base tree captured twice: once cold (the memo generates ``oracle(2)``) and once warm (the reader hits).
    ``cold_unknown`` keeps each test's unknown reasons from the cold capture (spec U8)."""

    def __init__(self, root: Path, cold_unknown: dict[str, tuple[str, ...]]) -> None:
        self.root = root
        self.cold_unknown = cold_unknown

    @property
    def store(self):
        return store_of(self.root)


def capture_base(root: Path, **variants: bool) -> Captured:
    build(root, **variants)
    capture(root)
    cold = {test: tuple(record.unknown) for test, record in store_of(root).tests().items()}
    capture(root)
    return Captured(root, cold)
