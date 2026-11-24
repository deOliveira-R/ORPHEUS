r"""G01e: the frame's rank cutoff has one address, and the frame's count equals ``pinv``'s nullity (N01, spec 7.2, 3.3).

Structure (``ast`` over ``numerics/frame_faces.py``): the float literal ``1e-12`` appears exactly once, as the value assigned to
``FRAME_GRAM_RANK_RTOL``, within three lines of ``machinery:``, ``ruling:`` and ``revisit:`` (the window of the ledger gate
``test_elegance_debt_is_tagged.py``); ``frame_gram_kernel`` reads ``np.linalg.svd(..., hermitian=True)`` WITH the vectors, never
``compute_uv=False`` and never ``eigvalsh``. The predicates are run on source copies with a forbidden form (the positive control).

Behaviour: on 201 Grams with one eigenvalue at ``f * 1e-12`` (``f`` in ``[0.9, 1.1]``, the tie of ``test_quotient_metric.py::test_B1c``),
``frame_gram_kernel(g).dim`` equals the nullity ``5 - round(trace(pinv(g, hermitian=True, rcond=FRAME_GRAM_RANK_RTOL) @ g))`` of an
oracle computed here (G01b scans ``metric.py`` only), not "the nullity of the inverse the body builds", which holds for any declared
count. Both nullities 0 and 1 occur. First red today: the literal appears 0 times in ``frame_faces.py`` and the rank cutoff lives in
``metric.py``.
"""

from __future__ import annotations

import ast
import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics import frame_faces as FF

pytestmark = pytest.mark.foundation

# tests/gates/numerics/test_quotient_metric.py:93 and :96 (copied: a gate owns its inputs).
_Q5 = np.linalg.qr(np.random.default_rng(0).standard_normal((5, 5)))[0]


def _gram(eigs):
    G = _Q5 @ np.diag(eigs) @ _Q5.T
    return (G + G.T) / 2


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def frame_faces_path() -> Path:
    spec = importlib.util.find_spec("orpheus.numerics.frame_faces")
    _require(spec is not None and spec.origin is not None, "orpheus.numerics.frame_faces is not found")
    assert spec is not None and spec.origin is not None
    return Path(spec.origin)


def tolerance_literals(source: str) -> list[int]:
    """The line of each float constant equal to ``1e-12`` in a parse of ``source``."""
    return sorted(n.lineno for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Constant) and isinstance(n.value, float) and n.value == 1e-12)


def literal_is_the_named_constant(source: str) -> bool:
    """True iff the one literal is the value of an assignment to ``FRAME_GRAM_RANK_RTOL``."""
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, (ast.Assign, ast.AnnAssign)) and isinstance(n.value, ast.Constant) and n.value.value == 1e-12:
            targets = n.targets if isinstance(n, ast.Assign) else [n.target]
            if any(isinstance(t, ast.Name) and t.id == "FRAME_GRAM_RANK_RTOL" for t in targets):
                return True
    return False


def tag_window_has_the_ledger_fields(source: str, line: int) -> bool:
    lines = source.splitlines()
    window = "\n".join(lines[max(0, line - 4): line + 3])
    return all(token in window for token in ("machinery:", "ruling:", "revisit:"))


def count_reads(source: str) -> dict[str, bool]:
    """How ``frame_gram_kernel`` reads the spectrum: hermitian ``svd`` with the vectors, a ``compute_uv`` keyword, ``eigvalsh``."""
    fn = next((n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.FunctionDef) and n.name == "frame_gram_kernel"), None)
    out = {"found": fn is not None, "svd_with_vectors": False, "compute_uv_keyword": False, "eigvalsh": False}
    if fn is None:
        return out
    for c in ast.walk(fn):
        if isinstance(c, ast.Call):
            name = c.func.attr if isinstance(c.func, ast.Attribute) else c.func.id if isinstance(c.func, ast.Name) else None
            kws = {k.arg for k in c.keywords}
            if name == "svd":
                out["compute_uv_keyword"] |= "compute_uv" in kws
                out["svd_with_vectors"] |= "hermitian" in kws and "compute_uv" not in kws
            if name == "eigvalsh":
                out["eigvalsh"] = True
    return out


def test_frame_tolerance_has_one_address() -> None:
    source = frame_faces_path().read_text()
    lines = tolerance_literals(source)
    _require(len(lines) == 1, f"the literal 1e-12 appears {len(lines)} times in frame_faces.py (lines {lines}), expected exactly 1")
    _require(literal_is_the_named_constant(source), "the one 1e-12 literal is not the value of FRAME_GRAM_RANK_RTOL")
    _require(tag_window_has_the_ledger_fields(source, lines[0]), "FRAME_GRAM_RANK_RTOL is not within three lines of machinery:, ruling: and revisit:")
    _require(FF.FRAME_GRAM_RANK_RTOL == 1e-12, f"FRAME_GRAM_RANK_RTOL = {FF.FRAME_GRAM_RANK_RTOL!r}")


def test_the_count_reads_the_decomposition_the_body_reads() -> None:
    r = count_reads(frame_faces_path().read_text())
    _require(r["found"], "frame_gram_kernel is not defined in frame_faces.py")
    _require(r["svd_with_vectors"], "frame_gram_kernel must read np.linalg.svd(..., hermitian=True) with the vectors")
    _require(not r["compute_uv_keyword"], "frame_gram_kernel passes compute_uv: the count must read the body's own call (with the vectors)")
    _require(not r["eigvalsh"], "frame_gram_kernel calls eigvalsh: a different routine reads 1 of 201 tie Grams across the cutoff")


_EIGVALSH = (r"np\.linalg\.svd\(\s*(\w+)\s*,\s*hermitian=True\s*\)", r"np.linalg.eigvalsh(\1)")
_COMPUTE_UV = (r"hermitian=True\s*\)", "hermitian=True, compute_uv=False)")


@pytest.mark.parametrize(
    "mutate,predicate",
    [
        (lambda s: s + "\n_EXTRA = 1e-12\n", lambda s: len(tolerance_literals(s)) == 1),
        (lambda s: re.sub(*_COMPUTE_UV, s), lambda s: not count_reads(s)["compute_uv_keyword"]),
        (lambda s: re.sub(*_EIGVALSH, s), lambda s: not count_reads(s)["eigvalsh"]),
    ],
    ids=["second_literal", "compute_uv_false", "eigvalsh"],
)
def test_the_predicates_read_red_on_a_forbidden_form(mutate, predicate) -> None:
    """Positive control: each predicate reads True on the real source and False with its forbidden form introduced."""
    source = frame_faces_path().read_text()
    _require(predicate(source), "the predicate must read True on the real source (so the control is a change of reading)")
    mutated = mutate(source)
    _require(mutated != source, "the mutation did not change the source: the control is a null control")
    _require(not predicate(mutated), "the predicate still reads True with the forbidden form introduced")


def test_the_frame_count_equals_the_pinv_nullity_across_the_cutoff() -> None:
    bad, counts = [], set()
    for f in np.linspace(0.9, 1.1, 201):
        g = _gram([1.0, 0.5, 0.3, 0.2, f * FF.FRAME_GRAM_RANK_RTOL])
        declared = FF.frame_gram_kernel(g).dim
        oracle = 5 - int(round(float(np.trace(np.linalg.pinv(g, hermitian=True, rcond=FF.FRAME_GRAM_RANK_RTOL) @ g))))
        counts.add(declared)
        if declared != oracle:
            bad.append((round(float(f), 4), declared, oracle))
    _require(not bad, f"the frame's count differs from pinv's nullity at {len(bad)} of 201: {bad[:5]}")
    _require(counts == {0, 1}, f"activation: the sweep must straddle the cutoff with declared dimensions {{0, 1}}, got {counts}")
