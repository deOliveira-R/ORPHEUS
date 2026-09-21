"""What a page costs, and what a budget is.

Tokens are estimated as characters / ``CHARS_PER_TOKEN``. Every generated
target carries a ``budget_tokens`` in its source page's ``harness:`` block,
and exceeding it fails generation, because the budget is the point of the
whole arrangement: the always-on cost every dispatch pays for a rule or an
agent's role block, the per-load cost for a skill.

How a budget is set: a round figure above the measured size after a review
pass, never more than ``SLACK_MAX`` tokens above it — generation fails on more
slack, so an over-generous budget cannot let growth in unreviewed; an agent's
role block takes one flat value per role (Key 500, Support 300) inside the
same slack. A raise is made in the same commit as the text that needs it,
which the front-matter placement makes literal. The budget is a size gate,
never a quality gate: the review is.

What a budget bounds is the page's own output: a whole file, or a spliced block
with its markers (the hand-maintained rest of that file is not the page's to
budget). What the always-on sum in ``__main__`` counts is the file as a session
loads it, whole.
"""
from __future__ import annotations

CHARS_PER_TOKEN = 3.6  # [M] 2026-09-20: 3.51 and 3.53 chars/token on the API tokenizer, keep-minus-omit probe (plan, "T4 executed")
SLACK_MAX = 400  # a budget more than this above the measured size lets growth in unreviewed


def tokens(text: str) -> int:
    return int(len(text) / CHARS_PER_TOKEN)


def check(text: str, budget: int, what: str) -> str | None:
    """The problem with ``text`` against ``budget``, or None: over it, or more
    than ``SLACK_MAX`` under it."""
    measured = tokens(text)
    if measured > budget:
        return f"{what}: ≈{measured} tokens > budget {budget}"
    if budget - measured > SLACK_MAX:
        return (f"{what}: budget {budget} is {budget - measured} tokens above ≈{measured}, "
                f"more than SLACK_MAX {SLACK_MAX}; lower it")
    return None
