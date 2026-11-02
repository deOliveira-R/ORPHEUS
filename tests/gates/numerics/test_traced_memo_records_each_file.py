r"""The recorder keeps every file's code, even when two files hold equal code objects.

Code objects compare equal across files: ``co_filename`` takes no part in ``==`` or the hash, so an empty
``__init__.py`` equals every other empty module, and a def equals an identical def at the same line of another
file. A recording that kept its code objects in a set kept only the first of each equal pair, and the manifest
lost the other file's pins: an edit to that file left the entry warm (found by the cone store's builder,
2026-11-01, on two empty package ``__init__`` modules). The recording now keys each code object by its file.
"""
from __future__ import annotations

import importlib
import sys
import textwrap

import pytest

from orpheus.numerics import traced_memo
from orpheus.numerics.traced_memo import DefPin, StatementsPin

pytestmark = pytest.mark.foundation

_TWIN = textwrap.dedent(
    """\
    def twin():
        return 1
    """
)


@pytest.fixture
def two_equal_files(tmp_path, monkeypatch):
    """Two packages with empty ``__init__`` modules, each holding a module whose ``twin`` def is the same text at
    the same line, so the two defs' code objects compare equal."""
    for name in ("twin_left", "twin_right"):
        package = tmp_path / name
        package.mkdir()
        (package / "__init__.py").write_text("")
        (package / "mod.py").write_text(_TWIN)
    monkeypatch.syspath_prepend(str(tmp_path))
    yield tmp_path
    for name in ("twin_left", "twin_left.mod", "twin_right", "twin_right.mod"):
        sys.modules.pop(name, None)


def _run_both() -> int:
    left = importlib.import_module("twin_left.mod")
    right = importlib.import_module("twin_right.mod")
    return left.twin() + right.twin()


def test_the_premise_equal_code_objects_from_two_files(two_equal_files):
    """The premise this file guards: the two defs' code objects are equal, as are two empty modules."""
    left = importlib.import_module("twin_left.mod")
    right = importlib.import_module("twin_right.mod")
    if left.twin.__code__ != right.twin.__code__:
        pytest.fail("premise lost: the two twin defs no longer compare equal, so this file guards nothing")
    if compile("", "a/__init__.py", "exec") != compile("", "b/__init__.py", "exec"):
        pytest.fail("premise lost: empty modules no longer compare equal")


@pytest.mark.catches("ERR-124")
@pytest.mark.rests_on("tests/gates/numerics/test_traced_memo_records_each_file.py::test_the_premise_equal_code_objects_from_two_files")
def test_both_files_of_an_equal_pair_are_pinned(two_equal_files):
    """A traced run through two equal defs in two files, importing two empty packages, pins all four files."""
    value, manifest = traced_memo.trace_call(_run_both)
    if value != 2:
        pytest.fail(f"the traced run returned {value}, not 2")
    module_paths = {pin.path for pin in manifest.pins if isinstance(pin, StatementsPin)}
    def_paths = {pin.path for pin in manifest.pins if isinstance(pin, DefPin) and pin.qualname == "twin"}
    wanted_modules = {"twin_left/__init__.py", "twin_right/__init__.py", "twin_left/mod.py", "twin_right/mod.py"}
    if not wanted_modules <= module_paths:
        pytest.fail(f"module pins lost to code equality: missing {sorted(wanted_modules - module_paths)}")
    if def_paths != {"twin_left/mod.py", "twin_right/mod.py"}:
        pytest.fail(f"the twin def is pinned in {sorted(def_paths)}, not in both files")
