"""The brief-rules block: the one place a copy of the rules is structurally forced.

The three Support agents that are launched without the project rules load no
rule, so a brief to one of them is the only channel a project rule has; the template on
``docs/development/workflows.md`` therefore carries a "Rules that apply to
you" list. Written by hand, that list is a copy of the rules that drifts
(review 2026-09-21, concept C7). Here it is assembled instead: every rule page
that declares ``harness.brief`` contributes one item, in page-name order, and
the block is spliced into the workflows page between GENERATED markers like a
role block, so ``--check`` reports a hand edit as drift. The block is a
source-to-source derivation, harness-independent: the page is docs, and what a
rule says to a rule-less agent is a fact of the rule.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from .render import Block, splice
from .source import DOCS_DEV, Kind, Page, rel

BRIEF_PAGE = DOCS_DEV / "workflows.md"
LABEL = "brief rules"
SOURCE = "the harness.brief of every page under docs/development/rules/"


def assemble(pages: Sequence[Page]) -> Block:
    """One list item per rule that declares a brief, in page-name order."""
    briefs = sorted((p.name, p.brief) for p in pages if p.kind is Kind.RULE and p.brief)
    items = [f"- `{name}`: {' '.join(text.split())}" for name, text in briefs]
    return Block(LABEL, "\n".join(items))


def render(pages: Sequence[Page], page: Path = BRIEF_PAGE) -> tuple[str, str | None]:
    """The workflows page as it should stand, and the splice problem if any."""
    current = page.read_text(encoding="utf-8") if page.exists() else ""
    return splice(assemble(pages), current, SOURCE, rel(page))
