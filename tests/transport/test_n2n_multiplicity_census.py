r"""XD-2 — the (n,2n) multiplicity's ONE home (CS4c step 3; design record
§12 F6, ruled option (i); verification plan §6.1).

The multiplicity :math:`\nu_{2n} = 2` is a physics constant of the
channel with exactly one production spelling —
:attr:`orpheus.transport.kernels.N2NKernel.multiplicity`. Every solver
family consumes THAT (datum consumption — recycling the vocabulary per
the sharpening-order law), so a thirteenth literal home is unspellable
without reddening this census.

Predicate (AST — ⛔ never grep for this: `.claude/rules/nexus-tools.md`
records ugrep's silent-zero anchor hazard): every ``ast.Constant`` with
value in ``{2, 2.0}`` (bool excluded) that is an operand of a
``BinOp(Mult)`` or the value of an ``AugAssign(Mult)``, inside a
function whose body text matches ``n2n|sig2|sig_2n|_2n\b``
(case-insensitive), outside ``kernels.py`` — validated against all four
historical spellings (``2.0 *``, integer ``2 *``, ``w *= 2.0``, and the
reversed ``* 2.0``) by the mandatory positive control (vv #17).
"""
from __future__ import annotations

import ast
import pathlib

import pytest

pytestmark = pytest.mark.foundation

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_NAME_NET = ("n2n", "sig2", "sig_2n", "_2n")

#: Named false-positive exclusions, each with its measured reason.
_EXCLUDED: set[tuple[str, int]] = set()
# orpheus/mc/solver.py — `phi = 2.0 * np.pi * rng.random()`: an
# azimuthal-angle sample that lands in the net only because the same
# `_random_walk` body also spells `sig_2n`. [M] the one false positive
# tree-wide (plan §6.1). Registered by (file, function) rather than a
# brittle line number:
_EXCLUDED_BY_FUNCTION = {("orpheus/mc/solver.py", "_random_walk", "pi")}


def _census(tree: ast.AST, rel: str):
    hits = []
    for func in ast.walk(tree):
        if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        body_text = ast.unparse(func).lower()
        if not any(tok in body_text for tok in _NAME_NET):
            continue
        for node in ast.walk(func):
            twos = []
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
                twos = [node.left, node.right]
            elif isinstance(node, ast.AugAssign) and isinstance(node.op, ast.Mult):
                twos = [node.value]
            for cand in twos:
                if (
                    isinstance(cand, ast.Constant)
                    and not isinstance(cand.value, bool)
                    and cand.value in (2, 2.0)
                ):
                    # π-sample exclusion: a 2·π literal is an angle, not a
                    # multiplicity — keyed by (file, function, "pi") so the
                    # exclusion cannot silently widen.
                    parent_src = ast.unparse(node)
                    if (
                        (rel, func.name, "pi") in _EXCLUDED_BY_FUNCTION
                        and "pi" in parent_src
                    ):
                        continue
                    hits.append((rel, cand.lineno, func.name, parent_src[:60]))
    return hits


def test_multiplicity_literal_census_is_empty_outside_kernels():
    hits = []
    for f in sorted((_ROOT / "orpheus").rglob("*.py")):
        rel = str(f.relative_to(_ROOT))
        if rel == "orpheus/transport/kernels.py":
            continue  # the ONE home (the ClassVar's own docstring row)
        tree = ast.parse(f.read_text())
        hits.extend(h for h in _census(tree, rel) if (h[0], h[1]) not in _EXCLUDED)
    if hits:
        listing = "\n".join(f"  {r}:{ln} in {fn}(): {src}" for r, ln, fn, src in hits)
        pytest.fail(
            f"{len(hits)} production multiplicity literal(s) outside "
            f"N2NKernel.multiplicity (XD-2, ruled option (i)):\n{listing}"
        )


def test_census_positive_control():
    """vv #17 — the filter finds all four historical spellings on a
    synthetic source; without this a broken filter and a clean tree
    print the same thing (L61)."""
    src = '''
def n2n_source_assembly(sigma, w, phi):
    a = 2.0 * sigma          # the classic spelling
    b = 2 * sigma            # the integer spelling (moc/core.py:316 class)
    w *= 2.0                 # the AugAssign spelling (mc/solver.py:447 class)
    c = sigma * 2.0          # the reversed operand
    return a + b + c + w
'''
    hits = _census(ast.parse(src), "synthetic.py")
    lines = sorted(h[1] for h in hits)
    if len(hits) != 4:
        pytest.fail(
            f"positive control: expected 4 synthetic hits "
            f"(got {len(hits)} at lines {lines}) — the census filter is broken"
        )
