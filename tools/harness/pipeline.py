"""Generation for one harness, then the two checks against the tree.

Every page the harness realises becomes an output at the path the harness
names, shaped as the harness says (a whole file written as-is; a block
spliced into whatever the target holds today), each measured against its
budget. **Drift** is an output that differs from disk. An **orphan** is a file
the harness owns that carries this generator's stamp, or a block marker, and
that no page produces any more.
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from . import budget
from .links import anchors, relink
from typing import assert_never

from .render import BEGIN_PREFIX, GENERATOR, STAMP_PREFIX, Block, WholeFile, block_text, splice
from .source import Page, rel
from .targets.base import Harness

# A stamp of ours names this generator; a stamp naming another generator (the
# error-catalogue index under the skills directory) is that generator's to reclaim.
_OURS = re.compile(rf"(?m)^(?:{re.escape(STAMP_PREFIX)}{re.escape(GENERATOR)} |{re.escape(BEGIN_PREFIX)})")


@dataclass(frozen=True)
class Output:
    page: Page
    path: Path
    text: str      # the file as it will stand on disk; what a session loads, so what the always-on sum counts
    budgeted: str  # what the budget bounds: the whole file, or the block with its markers (the hand-maintained rest is not the page's)


def generate(harness: Harness, pages: Sequence[Page]) -> tuple[dict[Path, Output], list[str]]:
    anchors.cache_clear()
    problems: list[str] = []
    mine = [p for p in pages if p.kind in harness.kinds]
    aliases = {p.path.resolve(): harness.target(p).resolve() for p in mine}
    outputs: dict[Path, Output] = {}
    for page in mine:
        dst = harness.target(page)
        if not any(dst == r or dst.is_relative_to(r) for r in harness.roots):
            problems.append(f"{page.rel}: the harness put it at {rel(dst)}, outside the roots it owns; an orphan there could never be reclaimed")
            continue
        if dst in outputs:
            problems.append(f"{rel(dst)}: produced by both {outputs[dst].page.rel} and {page.rel}")
            continue
        body = relink(page.body, page.path, dst, aliases, problems, page.rel)
        rendered = harness.render(page, body)
        match rendered:
            case WholeFile(text):
                out = budgeted = text
            case Block() as block:
                current = dst.read_text(encoding="utf-8") if dst.exists() else ""
                out, problem = splice(block, current, page.rel, rel(dst))
                if problem:
                    problems.append(problem)
                budgeted = block_text(block, page.rel)
            case _:
                assert_never(rendered)
        if problem := budget.check(budgeted, page.budget_tokens, f"{rel(dst)} ({page.kind.value})"):
            problems.append(problem)
        outputs[dst] = Output(page, dst, out, budgeted)
    return outputs, problems


def drift(outputs: dict[Path, Output]) -> list[Path]:
    return [p for p, o in outputs.items() if not p.exists() or p.read_text(encoding="utf-8") != o.text]


def orphans(harness: Harness, outputs: dict[Path, Output]) -> list[str]:
    candidates: list[Path] = []
    for root in harness.roots:
        if root.is_dir():
            candidates += [p for p in root.rglob("*.md") if "worktrees" not in p.relative_to(root).parts]
        elif root.is_file():
            candidates.append(root)
    return [f"{rel(p)}: carries a GENERATED stamp but no source page produces it"
            for p in candidates if p not in outputs and _OURS.search(p.read_text(encoding="utf-8"))]
