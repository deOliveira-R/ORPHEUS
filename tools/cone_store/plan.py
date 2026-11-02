"""What a driver tells a worker process: one JSON file, read by ``tools.cone_store.worker``."""
from __future__ import annotations

import json
import os
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path


def pycache_dir(repo: Path | str) -> str:
    """The bytecode cache the workers run with (``sys.pycache_prefix``), under the repository's own ``.cache``: a
    cache copied from another checkout holds code objects named after the file they were compiled from there
    (pytest's assertion-rewriting cache does not repair the name on load), and a cone is attributed by file name.
    The one place that path is spelled: each worker sets it, the cone drops what lives in it."""
    return os.path.join(os.path.realpath(repo), ".cache", "cone_store", "pycache")


@dataclass(frozen=True)
class Config:
    """The repository, its package, the memo's store (or ``None`` for the memo's default), the directory the
    workers write their results to, and the plan: for each test module, in order, the pytest arguments that
    collect its tests (a module path, or the node ids the caller named in it). ``fresh_only`` runs every test in
    its own interpreter (``--route fresh``). ``since_ns`` is the capture's start and ``recent`` the identity,
    at that moment, of each first-party file written within the memo's ``TRUSTED_AGE_NS`` before it (a path
    relative to the repository, as a list): a file whose identity differs from it at the end of a test changed
    during the capture."""

    repo: str
    package: str
    memo_root: str | None
    results: str
    groups: tuple[tuple[str, tuple[str, ...]], ...] = ()
    fresh_only: bool = False
    since_ns: int = 0
    recent: dict[str, list[int]] = field(default_factory=dict)

    def dump(self, path: Path) -> None:
        path.write_text(json.dumps(asdict(self)))

    @classmethod
    def load(cls, path: str | Path) -> Config:
        raw = json.loads(Path(path).read_text())
        raw["groups"] = tuple((module, tuple(args)) for module, args in raw["groups"])
        return cls(**raw)


def interpreter() -> list[str]:
    """This interpreter with this process's ``-O`` level: the environment key names it, so a worker must run under it."""
    return [sys.executable, *(["-O"] * sys.flags.optimize)]


def worker_command(role: str, config: Path, *arguments: str) -> list[str]:
    return [*interpreter(), "-m", "tools.cone_store.worker", role, str(config), *arguments]


def worker_environment(repo: Path) -> dict[str, str]:
    """The environment a worker starts in: the repository and the tools' own root on ``PYTHONPATH``, so
    ``tools.cone_store`` imports whatever the repository is (a fixture tree has no ``tools``). The bytecode cache
    (:func:`pycache_dir`) is set inside the worker (``sys.pycache_prefix``), never exported: an exported
    ``PYTHONPYCACHEPREFIX`` reached every interpreter a test starts and changed its verdict (M1, 2026-11-01:
    ``test_b17`` of the memo's source-index suite flipped)."""
    tools_root = str(Path(__file__).resolve().parents[2])
    roots = [str(repo), tools_root, *filter(None, [os.environ.get("PYTHONPATH")])]
    return {**os.environ, "PYTHONPATH": os.pathsep.join(roots)}
