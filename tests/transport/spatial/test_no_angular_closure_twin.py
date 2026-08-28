r"""The Morel–Montry relation is UNSPELLABLE in ``orpheus/transport/`` — permanently.

P4.9a's done-when, made a gate (the M8 anti-twin arm of
``scratch/p4_9a_verification_plan.md`` §4; the same promotion pattern P4.2
used for the L0 ladder).  Before the un-weld, ``DiamondDifference.update``
carried an inline copy of the closure's march — a Pattern-2 twin FORCED by
the layer (L2 may not import L3), which meant it could silently reappear
the same way.  This gate makes the regression grep-loud: the transport
layer's source may not spell the update relation's ``1 − τ`` kernel, nor
name the closure's march entry points.

[M] 2026-08-28, the pre-carve reading this gate would have caught:
``diamond.py:229`` spelled ``(psi_avg - (1.0 - tau) * ...) / tau`` — one
of six production spellings of the relation, and the only one in L2.

Two surfaces, scanned differently: the RELATION's ``1 − τ`` spelling is
forbidden in ALL text (comments and docstrings included — a
commented-out copy is a re-introduction waiting for an uncomment; the
twin this catches is an inline expression no import audit can see,
`plan-authoring` §6d's intra-file lesson).  The owner's entry-point
NAMES are forbidden only in CODE (AST names/imports/attributes) —
docstrings legitimately POINT readers at the owner
(``:func:`…march_psi_half_step```), and pointing is not shopping.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

import orpheus.transport as _transport_pkg

pytestmark = [pytest.mark.foundation]

_TRANSPORT_ROOT = Path(_transport_pkg.__file__).parent

# The relation's kernel in its spelling variants — forbidden in all text.
_FORBIDDEN_TEXT = (
    re.compile(r"1\.0\s*-\s*tau\b"),
    re.compile(r"\b1\s*-\s*tau\b"),
)
# The owner's entry points — forbidden as CODE names (an L2 module naming
# them in code went shopping for the closure again).
_FORBIDDEN_NAMES = frozenset({"advance_psi_half", "march_psi_half_step"})


def _code_names(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            yield node.id, node.lineno
        elif isinstance(node, ast.Attribute):
            yield node.attr, node.lineno
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            yield node.name, node.lineno
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                yield alias.name, node.lineno


def test_transport_source_cannot_spell_the_morel_montry_march() -> None:
    """[foundation] No transport module spells the M-M update or names its owner.

    Positive controls first (vv #17 / L61: a filter must prove it can see
    a known member of every shape before its zero means anything) — each
    pattern against its own witness: the OWNER module for the spellings
    the owner carries, a synthetic pre-carve line for the ``1 - tau``
    variant the owner legitimately never uses.
    """
    owner_path = _TRANSPORT_ROOT.parent / "sn" / "angular" / "closure.py"
    owner = owner_path.read_text()
    text_witnesses = {
        r"1\.0\s*-\s*tau\b": owner,
        r"\b1\s*-\s*tau\b": "psi = (avg - (1 - tau) * up) / tau",
    }
    for pat in _FORBIDDEN_TEXT:
        assert pat.search(text_witnesses[pat.pattern]) is not None, (
            f"positive control failed: {pat.pattern!r} matched neither its "
            f"witness — the filter is broken, not the tree clean"
        )
    owner_names = {n for n, _ in _code_names(ast.parse(owner))}
    assert _FORBIDDEN_NAMES <= owner_names, (
        "positive control failed: the owner module does not define the "
        "forbidden names — the AST filter is broken, not the tree clean"
    )

    offenders: list[str] = []
    for mod in sorted(_TRANSPORT_ROOT.rglob("*.py")):
        text = mod.read_text()
        rel = mod.relative_to(_TRANSPORT_ROOT.parent.parent)
        for i, line in enumerate(text.splitlines(), 1):
            for pat in _FORBIDDEN_TEXT:
                if pat.search(line):
                    offenders.append(f"{rel}:{i}: {line.strip()[:80]}")
        for name, lineno in _code_names(ast.parse(text)):
            if name in _FORBIDDEN_NAMES:
                offenders.append(f"{rel}:{lineno}: code names {name!r}")
    if offenders:
        pytest.fail(
            "the Morel-Montry relation leaked back into orpheus/transport/ "
            "(P4.9a made the angular axis the closure's alone):\n  "
            + "\n  ".join(offenders)
        )
