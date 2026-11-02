"""Which tests pytest would run: one ``--collect-only`` subprocess (the floor of ``select`` and ``capture``)."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections.abc import Sequence
from typing import NamedTuple
from pathlib import Path

from tools.cone_store.plan import Config, worker_command, worker_environment


def normalise(repo: Path, argument: str) -> str:
    """A pytest argument as it is spelled relative to the repository (an absolute path under it is made relative)."""
    path, separator, rest = argument.partition("::")
    if os.path.isabs(path):
        path = os.path.relpath(path, repo)
    return f"{path}{separator}{rest}"


class Collected(NamedTuple):
    """What a collection yields: the node ids pytest would run, in its order, and the modules whose collection failed
    (their tests are in no id: a module that does not import has none to list)."""

    ids: list[str]
    errors: list[str]
    #: Every conftest a full run loads, relative to the repository (:mod:`tools.cone_store.conftests`).
    conftests: list[str]


def collect_ids(repo: Path, package: str, arguments: Sequence[str] = ()) -> Collected:
    """The node ids pytest collects for ``arguments`` (default: the project's ``testpaths``), and the modules that failed to."""
    with tempfile.TemporaryDirectory(prefix="cone_list_") as results:
        config = Path(results, "config.json")
        Config(str(repo), package, None, results).dump(config)
        run = subprocess.run(worker_command("list", config, *(normalise(repo, a) for a in arguments)), cwd=repo,
                             env=worker_environment(repo), capture_output=True, text=True)
        listed, errors, found = Path(results, "ids.json"), Path(results, "errors.json"), Path(results, "conftests.json")
        if not listed.exists():
            raise RuntimeError(f"cone_store: pytest collected nothing (exit {run.returncode}):\n{run.stdout[-2000:]}{run.stderr[-2000:]}")
        return Collected(json.loads(listed.read_text()), json.loads(errors.read_text()), json.loads(found.read_text()))
