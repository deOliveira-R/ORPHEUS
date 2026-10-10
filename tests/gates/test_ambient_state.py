r"""Step 5 of #405 P3: no answer under ``orpheus/`` depends on ambient state a traced memo cannot key.

Spec ``.claude/plans/reference_p3_spec.md`` §1.5, gates M5.1–M5.4. A traced memo keys a call on its arguments
and validates it against the code that ran. Two kinds of state escape both: the process environment (read
at import or at call time) and a GLOBAL precision written by one function and read by the next (``mpmath.mp``).
The first enters a child unchanged and is in no key; the second makes an in-process (bypass or legacy) answer
depend on the order of calls. These gates read the tree by AST, with a positive control per shape (X1, X2).
"""
from __future__ import annotations

import ast
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[2]
ORPHEUS = REPO / "orpheus"

#: The one environment read allowed: the withdrawal switch decides whether a WITHDRAWN generator may run at
#: all (it refuses otherwise); it selects no value. A traced memo must therefore never wrap a withdrawn
#: generator (spec §1.5, rider R5.1).
ENVIRONMENT_READERS_ALLOWED = {
    "orpheus/derivations/common/withdrawal.py",
    "orpheus/numerics/traced_memo.py",
    # The recorder snapshots the environment a run RECEIVED, to pin its declared variables (#405 P3).
    "orpheus/numerics/_traced_memo_boot.py",
}
#: The census's positive control: a reader that exists on today's tree. (The traced memo reads the environment
#: only to hand a generating process its declared variables, and its recorder to pin them, M3.16; neither
#: selects a value.)
ENVIRONMENT_READER_KNOWN = "orpheus/derivations/common/withdrawal.py"

_ENV_ATTRIBUTES = {"environ", "environb", "getenv", "getenvb", "putenv", "unsetenv"}
_PRECISION_ATTRIBUTES = {"dps", "prec"}


def environment_reads(source: str) -> list[int]:
    """Lines reading or writing the process environment: ``os.environ``/``os.getenv`` (any receiver), and the
    same names imported from ``os`` (``from os import environ, getenv as g``) used bare."""
    tree = ast.parse(source)
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "os":
            imported |= {a.asname or a.name for a in node.names if a.name in _ENV_ATTRIBUTES}
    lines = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in _ENV_ATTRIBUTES:
            lines.append(node.lineno)
        elif isinstance(node, ast.Name) and node.id in imported:
            lines.append(node.lineno)
    return sorted(set(lines))


def global_precision_writes(source: str) -> list[int]:
    """Lines assigning a ``.dps``/``.prec`` attribute (``mp.dps = 30``, ``mpmath.mp.prec += 4``), or
    ``setattr(x, "dps"|"prec", ...)``: a global precision write. ``with mp.workdps(30):`` is the local form and
    is not one."""
    lines = []
    for node in ast.walk(ast.parse(source)):
        targets: list[ast.expr] = []
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AugAssign, ast.AnnAssign)):
            targets = [node.target]
        if isinstance(node, ast.stmt) and any(isinstance(t, ast.Attribute) and t.attr in _PRECISION_ATTRIBUTES for t in targets):
            lines.append(node.lineno)
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "setattr"
                and len(node.args) >= 2 and isinstance(node.args[1], ast.Constant) and node.args[1].value in _PRECISION_ATTRIBUTES):
            lines.append(node.lineno)
    return sorted(set(lines))


ENV_CONTROLS = {
    "os.environ.get": ("import os\nx = os.environ.get('A')\n", [2]),
    "os.getenv": ("import os\nx = os.getenv('A')\n", [2]),
    "environ-imported": ("from os import environ\nx = environ['A']\n", [2]),
    "getenv-aliased": ("from os import getenv as g\nx = g('A')\n", [2]),
    "aliased-module": ("import os as _os\nx = _os.environ.get('A', '0') != '1'\n", [2]),
    "none": ("import os\nx = os.path.join('a', 'b')\n", []),
}

PRECISION_CONTROLS = {
    "mp.mp.dps": ("import mpmath as mp\nmp.mp.dps = 30\n", [2]),
    "mpmath.mp.prec+=": ("import mpmath\nmpmath.mp.prec += 10\n", [2]),
    "setattr": ("import mpmath\nsetattr(mpmath.mp, 'dps', 30)\n", [2]),
    "workdps-is-local": ("import mpmath\nwith mpmath.workdps(30):\n    pass\n", []),
}


@pytest.mark.parametrize("shape", sorted(ENV_CONTROLS))
def test_m5_1a_the_environment_census_finds_every_shape(shape):
    """M5.1 (a), the census's positive controls: each spelling of an environment read is found on its line,
    and an ``os.path`` use is not."""
    source, lines = ENV_CONTROLS[shape]
    assert environment_reads(source) == lines


@pytest.mark.rests_on('tests/gates/test_ambient_state.py::test_m5_1a_the_environment_census_finds_every_shape')
def test_m5_1b_no_answer_reads_the_environment():
    """M5.1 (b): no module under ``orpheus/`` reads the process environment except the withdrawal switch. A
    generator reading it would answer by a variable no key holds (``ORPHEUS_SLAB_VIA_E1``, read at import by
    ``peierls_nystrom/cases.py``: ``[M]`` 2026-10-04, the first red). The input count is asserted, and the
    allowed reader must be found (the census's own control)."""
    files = sorted(ORPHEUS.rglob("*.py"))
    assert len(files) > 300, f"the census read {len(files)} files"
    found = {str(p.relative_to(REPO)): environment_reads(p.read_text()) for p in files}
    found = {k: v for k, v in found.items() if v}
    assert ENVIRONMENT_READER_KNOWN in found, f"the census missed the known reader: {sorted(found)}"
    offenders = {k: v for k, v in found.items() if k not in ENVIRONMENT_READERS_ALLOWED}
    assert not offenders, f"{len(offenders)} of {len(files)} modules read the environment: {offenders}"


@pytest.mark.parametrize("shape", sorted(PRECISION_CONTROLS))
def test_m5_2a_the_precision_census_finds_every_shape(shape):
    """M5.2 (a): each spelling of a global precision write is found; the local ``workdps`` context is not."""
    source, lines = PRECISION_CONTROLS[shape]
    assert global_precision_writes(source) == lines


@pytest.mark.rests_on('tests/gates/test_ambient_state.py::test_m5_2a_the_precision_census_finds_every_shape')
def test_m5_2b_no_function_writes_the_global_precision():
    """M5.2 (b): no module under ``orpheus/`` writes ``mpmath.mp``'s precision globally (``[M]`` 2026-10-04,
    the first red: ``cylinder_derivations.py:404`` and ``greens_function_slab.py:456``, both ``= 30``)."""
    files = sorted(ORPHEUS.rglob("*.py"))
    assert len(files) > 300
    offenders = {str(p.relative_to(REPO)): v for p in files if (v := global_precision_writes(p.read_text()))}
    assert not offenders, f"{len(offenders)} of {len(files)} modules write a global precision: {offenders}"


@pytest.mark.rests_on('tests/gates/test_ambient_state.py::test_m5_2b_no_function_writes_the_global_precision')
def test_m5_3_the_two_origins_derivations_leave_the_precision_as_they_found_it():
    """M5.3, the behaviour the census predicts: with ``mp.dps`` set to 17, each of the two derivations that
    wrote it leaves it 17 (``[M]`` 2026-10-04: both leave 30; ``probes/probe_mp_dps.log``). Run in a fresh
    interpreter so no earlier test's precision leaks in."""
    script = textwrap.dedent('''
        import mpmath
        from orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations import derive_bessel_wronskian_identity
        from orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab import derive_T00_equals_P_ss_slab
        for f in (derive_bessel_wronskian_identity, derive_T00_equals_P_ss_slab):
            mpmath.mp.dps = 17
            result = f()
            print(f.__name__, mpmath.mp.dps, bool(result["pass"]))
    ''')
    run = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, cwd=str(REPO))
    rows = [line.split() for line in run.stdout.splitlines()]
    assert len(rows) == 2, run.stderr[-2000:]
    for name, dps, passed in rows:
        assert dps == "17", f"{name} left mp.dps = {dps}"
        assert passed == "True", f"{name}: its own identity no longer passes"


@pytest.mark.rests_on('tests/gates/test_ambient_state.py::test_m5_1b_no_answer_reads_the_environment')
def test_m5_4_the_slab_routing_does_not_follow_the_environment():
    """M5.4: the Peierls slab routing is a setting of the call, not of the environment: two fresh interpreters,
    one with ``ORPHEUS_SLAB_VIA_E1=1`` and one without, import the cases module and agree on every module-level
    value it exposes (today they disagree on ``_SLAB_VIA_UNIFIED``: the first red). The two existing rows that
    pin the environment behaviour as the contract (``test_peierls_multigroup.py``,
    ``test_default_flag_is_unified`` and ``test_env_var_forces_native``) are INVERTED gates after the carve and
    are re-posed onto the explicit setting in the same commit (spec §1.5)."""
    script = textwrap.dedent('''
        import json
        from orpheus.derivations.continuous.peierls_nystrom import cases
        print(json.dumps({k: v for k, v in vars(cases).items() if isinstance(v, (bool, int, float, str)) and not k.startswith("__")}, sort_keys=True))
    ''')
    env = {k: v for k, v in os.environ.items() if k != "ORPHEUS_SLAB_VIA_E1"}
    plain = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, cwd=str(REPO), env=env)
    forced = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, cwd=str(REPO),
                            env=env | {"ORPHEUS_SLAB_VIA_E1": "1"})
    assert plain.stdout and forced.stdout, plain.stderr[-1500:] + forced.stderr[-1500:]
    assert plain.stdout == forced.stdout, f"the environment moved the module: {plain.stdout} != {forced.stdout}"
