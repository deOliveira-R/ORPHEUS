"""Repository files named in code spans, checked for existence.

A code span that names a file of this repository claims the file exists, and a
rename or a deletion makes the claim false with no build warning: a code span
is not a link, and the ID check (``ids.py``) skips code by design. Three
spellings are read, outside fences:

* a **path**: a span with a ``/`` whose first segment is a top-level directory
  holding tracked files, or whose last segment carries a file extension
  (``orpheus/sn/solver.py``, ``sn/sweep/scan.py``, a skill's ``scripts/x.md``;
  a ``::test``, ``:line`` or ``#anchor`` suffix is dropped), which must be a
  tracked file or a directory holding one, or the tail of one (so a path
  relative to a page's own directory, or shortened from the root, resolves);
* a **dotted module** whose first segment is such a directory
  (``tools.harness.ids``), of which the longest prefix that is a module or a
  package must reach two segments or more (what follows it is attributes,
  unchecked);
* a **bare file name** with an extension (``error_catalog.rst``), which must name
  some tracked file, anywhere, so it cannot tell two files of one name apart.

A file the reader is to create, and an example, carry a placeholder
(``diag_<NN>_<step>.py``, ``tests/<method>/``): a span holding ``<``, ``>``,
``*``, ``{``, ``}``, an ellipsis or whitespace is not read. Not read either,
stated so the zero is read for what it is: fences; a bare identifier (``SNMesh``), which no
existence test can tell from prose. A path git ignores (the literature folder)
is local by design and is not read; a path that is neither tracked nor ignored
is dead for every other checkout, whatever this one holds.
"""
from __future__ import annotations

import re
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from .ids import _CODE, _FENCE
from .source import REPO_ROOT, Page

_PLACEHOLDER = re.compile(r"[<>*{}\s]|\.\.\.|…")
_SUFFIX = re.compile(r"::.*$|:\d+(?:-\d+)?$|#.*$")
_BARE = re.compile(r"[\w.-]+\.(?:md|rst|py|json|toml|ya?ml|sh|txt|cfg|ini)")
_DOTTED = re.compile(r"\w+(?:\.\w+)+")


@dataclass(frozen=True)
class Tree:
    files: frozenset[str]
    dirs: frozenset[str]
    names: frozenset[str]
    top: frozenset[str]
    tails: frozenset[str]   # every trailing run of segments of a tracked file or directory

    @classmethod
    def tracked(cls) -> Tree:
        out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout
        files = frozenset(f for f in out.split("\0") if f)
        dirs = frozenset("/".join(f.split("/")[:i]) for f in files for i in range(1, f.count("/") + 1))
        tails = frozenset("/".join(p.split("/")[i:]) for p in files | dirs for i in range(p.count("/") + 1))
        return cls(files, dirs, frozenset(PurePosixPath(f).name for f in files),
                   frozenset(d for d in dirs if "/" not in d), tails)

    def resolves(self, span: str) -> bool | None:
        """Whether ``span`` names something tracked; None when it is not read as a file."""
        span = span.strip().rstrip(".,;:)")
        if not span or _PLACEHOLDER.search(span):
            return None
        span = _SUFFIX.sub("", span)
        if "/" in span:
            path = span.rstrip("/")
            if path.split("/", 1)[0] not in self.top and not _BARE.fullmatch(path.rsplit("/", 1)[-1]):
                return None
            return path in self.tails
        if _DOTTED.fullmatch(span) and span.split(".", 1)[0] in self.top and not _BARE.fullmatch(span):
            segs = span.split(".")
            return any(f"{'/'.join(segs[:n])}.py" in self.files or f"{'/'.join(segs[:n])}/__init__.py" in self.files
                       or "/".join(segs[:n]) in self.dirs for n in range(len(segs), 1, -1))
        if _BARE.fullmatch(span):
            return span in self.names
        return None


def spans(text: str) -> list[tuple[str, int]]:
    """Every code span outside a fence, with its line."""
    masked = _FENCE.sub(lambda m: "\n" * m.group(0).count("\n"), text)
    return [(m.group(0)[1:-1], masked.count("\n", 0, m.start()) + 1) for m in _CODE.finditer(masked)]


def ignored(paths: Iterable[str], root: Path = REPO_ROOT) -> frozenset[str]:
    """The paths git ignores: local by design, so a citation of one is not dead. Each is asked as
    written and as a directory: a pattern with a trailing slash (``scratch/literature/``) matches a
    path git knows to be a directory, and on a checkout where the directory is absent (CI's) only the
    slashed spelling tells it so (`[M]` 2026-09-22: the first push of this check read the literature
    folder dead in CI and ignored here, where the folder exists)."""
    paths = [p.rstrip("/") for p in paths]
    query = "\0".join(q for p in paths for q in (p, p + "/"))
    if not query:
        return frozenset()
    out = subprocess.run(["git", "check-ignore", "--no-index", "--stdin", "-z"], cwd=root, input=query + "\0",
                         capture_output=True, text=True, check=False).stdout
    hits = {h for h in out.split("\0") if h}
    return frozenset(p for p in paths if p in hits or p + "/" in hits)


def check(pages: Iterable[Page], tree: Tree | None = None) -> tuple[int, list[str]]:
    """(the number of spans that resolved, the problems): a span naming a file the tree does not hold."""
    tree = tree or Tree.tracked()
    resolved = 0
    dead: list[tuple[str, int, str]] = []
    for page in pages:
        offset = page.path.read_text(encoding="utf-8").count("\n") - page.body.count("\n")  # the front matter's lines
        for span, line in spans(page.body):
            line += offset
            match tree.resolves(span):
                case True:
                    resolved += 1
                case False:
                    dead.append((page.rel, line, span))
                case None:
                    pass
    local = ignored(_SUFFIX.sub("", s.strip().rstrip(".,;:)")).rstrip("/") for _, _, s in dead if "/" in s)
    problems = [f"{where}:{line}: `{span}` names no tracked file (write a file to be created, or an example, with a <placeholder>)"
                for where, line, span in dead if _SUFFIX.sub("", span.strip().rstrip(".,;:)")).rstrip("/") not in local]
    return resolved, problems
