"""How a relative link moves from the source's directory to the target's, and
which anchors MyST mints for a page.

Every linked file must exist, and an anchor into a ``.md`` page must be one
MyST mints for it: the anchors are computed with MyST's own slug function and
de-duplication over headings up to ``HEADING_ANCHOR_LEVELS`` (which
``docs/conf.py`` reads from here), so the check agrees with the docs build by
construction. A link to a page that is itself generated is written to that
page's generated copy (``aliases``): the harness copy is the one already in an
agent's context.
"""
from __future__ import annotations

import functools
import os
import re
from pathlib import Path

from markdown_it import MarkdownIt
from markdown_it.tree import SyntaxTreeNode
from myst_parser.mdit_to_docutils.base import compute_unique_slug

HEADING_ANCHOR_LEVELS = 4  # docs/conf.py: myst_heading_anchors = this
_LINK = re.compile(r"\]\((?!https?://|mailto:)([^)\s#]*)(#[^)]*)?\)")
_FRONT = re.compile(r"\A---\n.*?\n---\n", re.S)


@functools.lru_cache(maxsize=None)
def anchors(md: Path) -> frozenset[str]:
    """The anchors MyST mints for ``md``: parsed, not regexed, so a ``#`` inside
    a fence is not a heading; slugged and de-duplicated by MyST's own function."""
    tree = SyntaxTreeNode(MarkdownIt("commonmark").parse(_FRONT.sub("", md.read_text(encoding="utf-8"))))
    slugs: list[str] = []
    for node in tree.children:
        if node.type == "heading" and int(node.tag[1]) <= HEADING_ANCHOR_LEVELS:
            slugs.append(compute_unique_slug(node, slugs))
    return frozenset(slugs)


def relink(body: str, src: Path, dst: Path, aliases: dict[Path, Path], problems: list[str], where: str) -> str:
    """Re-point every relative link in ``body`` from ``src``'s directory to
    ``dst``'s. A missing file or a missing anchor is appended to ``problems``
    under ``where``, and that link is left as written."""
    def sub(m: re.Match) -> str:
        path, anchor = m.group(1), m.group(2) or ""
        target = (src.parent / path).resolve() if path else src.resolve()
        if not target.is_file():
            problems.append(f"{where}: link to missing file {path}")
            return m.group(0)
        if anchor and target.suffix == ".md" and anchor[1:] not in anchors(target):
            problems.append(f"{where}: link to missing heading {path}{anchor}")
        written = aliases.get(target, target)
        new = Path(os.path.relpath(written, dst.parent.resolve())).as_posix() if path else ""
        return f"]({new}{anchor})"
    return _LINK.sub(sub, body)
