"""Shared instruments of the memo-digest gates G5 to G7, G16 to G18 (specification part B 1.3, 4.2): a temporary copy of the memo and
values packages, a textual mutation that refuses to run when its anchor is absent, and a generating-process child.

A child is ``[sys.executable, "-P", "-c", program]`` run from a directory OUTSIDE the tree, so that the envelope's ``sys_path`` alone
decides what ``orpheus`` is. The program loads the boot by path under its dotted name (the way ``tools/cone_store/_recorder.py`` does),
exits at once unless the store root the boot computes lies under the temporary directory (the live ``.cache`` is never written, the
store root is three directories above the boot file), runs ``_boot()`` with the envelope on standard input, and prints the facts of
its ``BODY`` to standard error as one JSON line (``_boot`` points standard output at standard error).

This is a helper module, outside any collected test module: it carries no bare ``assert`` (the suite runs under ``python -O``).
"""
from __future__ import annotations

import json
import pickle
import shutil
import struct
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ORPHEUS = REPO / "orpheus"


class Uninstallable(Exception):
    """A mutation or a tree could not be built (an anchor absent or present more than once): the arm did not run."""


def build_tree(dest: Path, *, extra_values: dict[str, str] | None = None, memo_edit: tuple[str, str, str] | None = None,
               values_edit: tuple[str, str, str] | None = None, source: Path = ORPHEUS) -> Path:
    """Copy ``<source>/__init__.py``, ``<source>/memo/*.py`` and ``<source>/values/*.py`` (``source`` is the real ``orpheus`` package by default) to ``dest/orpheus`` and return ``dest``.

    ``extra_values`` adds ``orpheus/values/<name>``; an edit is ``(file name, anchor, replacement)`` and raises :class:`Uninstallable`
    unless the anchor occurs exactly once in the copied file."""
    package = dest / "orpheus"
    sources = {"memo": sorted((source / "memo").glob("*.py")), "values": sorted((source / "values").glob("*.py"))}
    if not sources["memo"] or not sources["values"]:
        raise Uninstallable(f"the tree to copy is empty: memo {len(sources['memo'])} files, values {len(sources['values'])} files")
    package.mkdir(parents=True)
    shutil.copy(ORPHEUS / "__init__.py", package / "__init__.py")
    for name, files in sources.items():
        (package / name).mkdir()
        for path in files:
            shutil.copy(path, package / name / path.name)
    for name, text in (extra_values or {}).items():
        (package / "values" / name).write_text(text)
    for subdir, edit in (("memo", memo_edit), ("values", values_edit)):
        if edit is not None:
            mutate(package / subdir / edit[0], edit[1], edit[2])
    return dest


def mutate(path: Path, anchor: str, replacement: str) -> None:
    """Replace ``anchor`` by ``replacement`` in the COPY ``path``; :class:`Uninstallable` unless the anchor occurs exactly once."""
    if not path.is_file():
        raise Uninstallable(f"{path} does not exist")
    text = path.read_text(encoding="utf-8")
    if text.count(anchor) != 1:
        raise Uninstallable(f"anchor occurs {text.count(anchor)} times in {path.name} (exactly one required): {anchor!r}")
    path.write_text(text.replace(anchor, replacement), encoding="utf-8")


_PREFIX = '''
import importlib.util, json, os, sys
R = {tree!r}
path = os.path.join(R, "orpheus", "memo", "_traced_memo_boot.py")
spec = importlib.util.spec_from_file_location("orpheus.memo._traced_memo_boot", path)
module = importlib.util.module_from_spec(spec)
sys.modules["orpheus.memo._traced_memo_boot"] = module
spec.loader.exec_module(module)
if not os.path.realpath(module.CACHE).startswith(os.path.realpath(R)):
    raise SystemExit("store root outside the tree: " + module.CACHE)
facts = {{"orpheus_before": "orpheus" in sys.modules, "cache": module.CACHE}}
recording, answers, jobs, traced_memo = module._boot()
import orpheus
if not os.path.realpath(orpheus.__file__).startswith(os.path.realpath(R)):
    raise SystemExit("orpheus resolves outside the tree: " + orpheus.__file__)
facts["orpheus_file_in_tree"] = True
'''
_SUFFIX = "\nsys.stderr.write('\\nFACTS ' + json.dumps(facts, default=str) + '\\n')\nsys.stderr.flush()\nos._exit(0)\n"


def generating_child(tree: Path, body: str, workdir: Path) -> dict[str, object]:
    """Run ``body`` after the real ``_boot()`` of the boot copy in ``tree``; return the ``facts`` dict it filled.

    ``body`` sees ``module`` (the boot), ``facts``, ``traced_memo``, ``R``. The envelope is ``{"sys_path": [tree, *this process's
    sys.path]}``. A non-zero exit raises ``RuntimeError`` with the child's standard error."""
    workdir.mkdir(parents=True, exist_ok=True)
    envelope = pickle.dumps({"sys_path": [str(tree), *sys.path]})
    program = _PREFIX.format(tree=str(tree)) + body + _SUFFIX
    done = subprocess.run([sys.executable, "-P", "-c", program], input=struct.pack(">Q", len(envelope)) + envelope,
                          capture_output=True, cwd=workdir)
    err = done.stderr.decode(errors="replace")
    if done.returncode != 0 or "FACTS " not in err:
        raise RuntimeError(f"generating child exited {done.returncode}:\n{err[-3000:]}")
    return json.loads(err.rsplit("FACTS ", 1)[1].splitlines()[0])


def cold_child(tree: Path, program: str) -> str:
    """Run ``program`` with ``-c`` from ``tree`` (so ``tree`` is first on ``sys.path``); its first statements assert that ``orpheus``
    resolves inside ``tree``; return its standard output."""
    guard = (f"import os, orpheus\nif not os.path.realpath(orpheus.__file__).startswith(os.path.realpath({str(tree)!r})):\n"
             "    raise SystemExit('orpheus resolves outside the tree: ' + orpheus.__file__)\n")
    done = subprocess.run([sys.executable, "-c", guard + program], capture_output=True, text=True, cwd=tree)
    if done.returncode != 0:
        raise RuntimeError(f"cold child exited {done.returncode}:\n{done.stderr[-3000:]}")
    return done.stdout
