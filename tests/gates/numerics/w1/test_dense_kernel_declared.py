r"""G01c, G01d and G01e: the rank of a dense form is declared, the declaration is exact, and the frame's cutoff has one address (N01, spec 7.2, 3.3).

**G01c.** ``DenseMetric(matrix, *, kernel)``: ``kernel`` is required and is ``Coord(n, ())`` (full rank) or a ``Spanned``. The rank on
elements of shape ``shape`` is ``(n - kernel.dim) * (size // n)``, ``size`` the product of ``shape``; the kernel basis has
``kernel.dim`` columns (the signature ``(n - kernel.dim, kernel.dim, 0)`` is declared, and read here through these two counts).
The eigenvalue ``5e-13`` of ``diag([1, 5e-13])`` is BELOW the retired cutoff ``1e-12`` and is a kept direction when the form
declares full rank: the "too tight" defect of #634 reading it as dropped is the red the inverse row catches. First red today:
``DenseMetric(np.diag([1.0, 5e-13])).rank((2,))`` reads 1 (the cutoff), and the ``kernel`` keyword does not exist.

**G01d.** The form constructor does not test ``G K = 0`` in floating point (3.3); this gate does, with its tolerance in the gate. For
each of the 14 deficient Grams of the harness group (c), with ``w, v = np.linalg.eigh(g_sym)`` (ascending, so the kernel is the FIRST
``kernel.dim`` columns): ``max |G K| <= 1e-14 * norm(G)`` on the 13 frame-built or exactly deficient ones and
``<= FRAME_GRAM_RANK_RTOL * norm(G)`` on the synthetic near-cutoff Gram (its dropped value is 5e-13 by construction, so it is the one
member of the second tolerance); ``w[0] >= -1e-14 * norm(G)`` on all 21 (no negative direction). Wrong-declaration arms: the declared
dimension must separate the singular spectrum, ``s[n - dim] < 1e-6 s_max`` (what a declaration inverts) and ``s[n - dim - 1] > 1e-6 s_max``
(what it zeroes); the arms ``dim - 1`` (14 of 14 deficient) and ``dim + 1`` (21 of 21) must read red. The indefinite control
``diag([1, -1])`` reads red on the smallest-eigenvalue predicate. The declared kernel is read LIVE through ``frame_gram_kernel`` (and the
harness' own declaration for the four hand-written Grams), so a production mutation of the frame's count reddens these rows.

**G01e.** Structure (``ast`` over ``numerics/frame_faces.py``): the float literal ``1e-12`` appears exactly once, as the value assigned
to ``FRAME_GRAM_RANK_RTOL``, within three lines of ``machinery:``, ``ruling:`` and ``revisit:`` (the window of the ledger gate
``test_elegance_debt_is_tagged.py``); ``frame_gram_kernel`` reads ``np.linalg.svd(..., hermitian=True)`` WITH the vectors, never
``compute_uv=False`` and never ``eigvalsh``. The predicates are run on source copies with a forbidden form (the positive control).
Behaviour: on 201 Grams with one eigenvalue at ``f * 1e-12`` (``f`` in ``[0.9, 1.1]``, the tie of ``test_quotient_metric.py::test_B1c``),
``frame_gram_kernel(g).dim`` equals the nullity ``5 - round(trace(pinv(g, hermitian=True, rcond=FRAME_GRAM_RANK_RTOL) @ g))`` of an oracle
computed here (G01b scans ``metric.py`` only), not "the nullity of the inverse the body builds", which holds for any declared count; both
nullities 0 and 1 occur. First red today: the literal appears 0 times in ``frame_faces.py`` and the rank cutoff lives in ``metric.py``.
"""

from __future__ import annotations

import ast
import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics import frame_faces as FF
from orpheus.numerics.declaration import Ann, Coord, Spanned
from orpheus.numerics.metric import DenseMetric
from tests.gates._dense import declared_kernel
from tests.gates.numerics.w1 import _stored_grams as sg

pytestmark = pytest.mark.foundation

_RANK2_3X3 = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


@pytest.mark.rests_on("tests/gates/numerics/w1/test_declaration_normal_form.py::test_annihilator_normal_form")
def test_rank_is_declared() -> None:
    near = DenseMetric(np.diag([1.0, 5e-13]), kernel=Coord(2, ()))
    _require(near.rank((2,)) == 2, f"a declared-full-rank form with a 5e-13 value has rank {near.rank((2,))}, expected 2 (no cutoff)")
    _require(near.rank((2, 3)) == 6, f"trailing multiplicity: rank((2, 3)) = {near.rank((2, 3))}, expected 6")
    _require(near.kernel_basis((2,)).shape == (2, 0), f"kernel basis shape {near.kernel_basis((2,)).shape}, expected (2, 0)")
    deficient = DenseMetric(_RANK2_3X3, kernel=declared_kernel(_RANK2_3X3, 1))
    _require(deficient.rank((3,)) == 2, f"rank((3,)) = {deficient.rank((3,))}, expected 2")
    _require(deficient.rank((3, 2)) == 4, f"rank((3, 2)) = {deficient.rank((3, 2))}, expected (3 - 1) * (6 // 3) = 4")
    _require(deficient.kernel_basis((3,)).shape == (3, 1), f"kernel basis shape {deficient.kernel_basis((3,)).shape}, expected (3, 1)")
    _require(deficient.kernel_basis((3, 2)).shape == (6, 2), f"kernel basis shape {deficient.kernel_basis((3, 2)).shape}, expected (6, 2)")
    k = deficient.kernel_basis((3,))[:, 0]
    _require(float(np.max(np.abs(_RANK2_3X3 @ k))) <= 1e-14, "the kernel basis column is not in the kernel of the matrix")


def test_kept_direction_below_the_old_cutoff_is_inverted() -> None:
    """A kept ``5e-13`` is inverted byte for byte: ``inverse_matrix`` equals ``np.diag([1, 1 / 5e-13])`` (the #634 defect reversed)."""
    m = DenseMetric(np.diag([1.0, 5e-13]), kernel=Coord(2, ()))
    expected = np.diag([1.0, 1.0 / 5e-13])
    _require(np.asarray(m.inverse_matrix).tobytes() == expected.tobytes(), f"inverse face {np.asarray(m.inverse_matrix).tolist()} differs from {expected.tolist()}")


@pytest.mark.rests_on("tests/gates/numerics/w1/test_declaration_normal_form.py::test_annihilator_normal_form")
def test_kernel_declaration_refusals() -> None:
    d2 = np.diag([1.0, 2.0])
    with pytest.raises(ValueError, match="kernel is not declared by a spanned subspace"):
        DenseMetric(d2, kernel=Coord(2, (1,)))
    with pytest.raises(ValueError, match="kernel is not declared by a spanned subspace"):
        DenseMetric(d2, kernel=Coord(2, (0, 1)))  # the whole space: the zero form stays unspellable
    with pytest.raises(ValueError, match="kernel is declared on"):
        DenseMetric(d2, kernel=Coord(3, ()))
    with pytest.raises(ValueError, match="kernel is declared on"):
        DenseMetric(d2, kernel=Spanned(3, ("k",), 1))
    with pytest.raises(ValueError, match="kernel is not declared by a spanned subspace"):
        DenseMetric(d2, kernel=Ann(2, Spanned(2, ("k",), 1)))
    with pytest.raises(TypeError):
        DenseMetric(d2)  # type: ignore[call-arg]


def test_a_declared_kernel_is_part_of_the_identity() -> None:
    """The kernel is a compare field (P2-6): equal matrices with different declared kernels are different forms."""
    from orpheus.values.content import content_digest

    a = DenseMetric(_RANK2_3X3, kernel=declared_kernel(_RANK2_3X3, 1))
    b = DenseMetric(_RANK2_3X3.copy(), kernel=declared_kernel(_RANK2_3X3, 1))
    c = DenseMetric(_RANK2_3X3, kernel=Spanned(3, ("another-key",), 1))
    _require(content_digest(a) == content_digest(b), "one form built twice must have one digest")
    _require(content_digest(a) != content_digest(c), "two kernels over one matrix must be two forms (two digests)")


# ============================================================
# G01d
# ============================================================

_TIGHT = 1e-14
_GAP = 1e-6
_NEAR_CUTOFF = "syn.near_cutoff_5x5"


def kernel_residual(g_sym: np.ndarray, dim: int) -> float:
    """``max |G K| / norm(G)`` over the first ``dim`` eigenvectors (ascending eigenvalues) of the symmetric ``g_sym``."""
    _, v = np.linalg.eigh(g_sym)
    return float(np.max(np.abs(g_sym @ v[:, :dim])) / np.linalg.norm(g_sym)) if dim else 0.0


def smallest_eigenvalue_over_norm(g_sym: np.ndarray) -> float:
    return float(np.linalg.eigvalsh(g_sym)[0] / np.linalg.norm(g_sym))


def spectrum_is_separated_by(g_sym: np.ndarray, dim: int) -> bool:
    """True iff declaring ``dim`` neither inverts a value below ``1e-6 s_max`` nor zeroes one above it (exact NaN-safe form)."""
    s = np.linalg.svd(g_sym, hermitian=True)[1]
    n, top = s.shape[0], s[0]
    inverts_small = n - dim > 0 and not (s[n - dim - 1] > _GAP * top)  # the smallest KEPT value is not above the gap
    zeroes_large = dim > 0 and not (s[n - dim] < _GAP * top)  # the largest DROPPED value is not below the gap
    return not inverts_small and not zeroes_large


@pytest.mark.parametrize("family", sg.families())
@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_the_frame_count_equals_the_pinv_nullity_across_the_cutoff")
def test_declared_kernel_is_exact_for_the_frames(family: str) -> None:
    e = sg.stored(family)
    g = e.g_sym
    dim = sg.declared_for(e).dim
    if dim:
        tol = FF.FRAME_GRAM_RANK_RTOL if family == _NEAR_CUTOFF else _TIGHT
        r = kernel_residual(g, dim)
        _require(r <= tol, f"[{family}] max|G K|/norm(G) = {r:.3e} over the declared {dim}-dimensional kernel, above {tol:g}")
    low = smallest_eigenvalue_over_norm(g)
    _require(low >= -_TIGHT, f"[{family}] smallest eigenvalue / norm(G) = {low:.3e}: a negative direction")
    _require(spectrum_is_separated_by(g, dim), f"[{family}] the declared dimension {dim} does not separate the singular spectrum")


def test_the_arms_read_red_on_a_wrong_dimension() -> None:
    """A wrong declaration reads red: ``dim - 1`` on 14 of 14 deficient Grams, ``dim + 1`` on 21 of 21."""
    fams = sg.families()
    deficient = [f for f in fams if sg.declared_for(sg.stored(f)).dim > 0]
    _require(len(deficient) == 14, f"{len(deficient)} deficient families, expected 14")
    low_red = [f for f in deficient if not spectrum_is_separated_by(sg.stored(f).g_sym, sg.declared_for(sg.stored(f)).dim - 1)]
    high_red = [f for f in fams if not spectrum_is_separated_by(sg.stored(f).g_sym, sg.declared_for(sg.stored(f)).dim + 1)]
    _require(len(low_red) == 14, f"the arm dim - 1 read red on {len(low_red)} of 14, not red on {sorted(set(deficient) - set(low_red))}")
    _require(len(high_red) == 21, f"the arm dim + 1 read red on {len(high_red)} of 21, not red on {sorted(set(fams) - set(high_red))}")


def test_the_indefinite_control_reads_red() -> None:
    g = np.diag([1.0, -1.0])
    low = smallest_eigenvalue_over_norm(g)
    _require(low < -_TIGHT, f"the control diag([1, -1]) reads {low:.3e}, not negative")
    _require(not (low >= -_TIGHT), "the smallest-eigenvalue predicate must read red on the indefinite control")


# ============================================================
# G01e
# ============================================================

# tests/gates/numerics/test_quotient_metric.py:93 and :96 (copied: a gate owns its inputs).
_Q5 = np.linalg.qr(np.random.default_rng(0).standard_normal((5, 5)))[0]


def _gram(eigs):
    G = _Q5 @ np.diag(eigs) @ _Q5.T
    return (G + G.T) / 2


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


def test_the_rank_rule_counts_an_exact_tie_as_dropped() -> None:
    """The rank rule is ``<=``: a value exactly at ``FRAME_GRAM_RANK_RTOL`` of the largest is dropped, as ``pinv`` at that ``rcond`` does.

    ``diag([1, 1e-12])`` is the exact tie (the product ``1e-12 * 1.0`` is the stored value). A mutation to ``<`` declares ``Coord(2, ())``.
    """
    g = np.diag([1.0, FF.FRAME_GRAM_RANK_RTOL])
    declared = FF.frame_gram_kernel(g)
    oracle = 2 - int(round(float(np.trace(np.linalg.pinv(g, hermitian=True, rcond=FF.FRAME_GRAM_RANK_RTOL) @ g))))
    _require(oracle == 1, f"the pinv oracle reads nullity {oracle} on the tie, expected 1 (the control that the tie is a tie)")
    _require(declared.dim == 1, f"the tie is declared with dimension {declared.dim}, expected 1 (a value at the cutoff is dropped)")


def test_the_frame_kernel_key_is_the_digest_of_the_symmetrised_gram() -> None:
    """The content key is ``('frame-gram-kernel', digest)``, and two different deficient Grams of equal ``n`` and ``dim`` are unequal."""
    from orpheus.values.content import content_digest

    g1 = _gram([1.0, 0.5, 0.3, 0.0, 0.0])
    g2 = _gram([1.0, 0.5, 0.7, 0.0, 0.0])
    k1, k2 = FF.frame_gram_kernel(g1), FF.frame_gram_kernel(g2)
    _require(isinstance(k1, Spanned) and isinstance(k2, Spanned) and k1.dim == k2.dim == 2 and k1.n == k2.n == 5, "the two Grams are not a same-(n, dim) pair")
    sym = DenseMetric.symmetrised(g1)
    sym.setflags(write=False)
    _require(k1.key == ("frame-gram-kernel", content_digest(sym)), f"the kernel key {k1.key} is not ('frame-gram-kernel', digest of the symmetrised Gram)")
    _require(k1 != k2, "two different deficient Grams of equal n and dim declared equal kernels")
    _require(k1 == FF.frame_gram_kernel(g1.copy()), "one Gram declared twice must declare one kernel")


def test_inverse_of_honours_a_declaration_that_differs_from_the_retired_count() -> None:
    """The production door reads the declared count and no cutoff: ``inverse_of(diag([1, 5e-13]), kernel=Coord(2, ()))`` inverts ``5e-13`` (#634)."""
    m = DenseMetric.inverse_of(np.diag([1.0, 5e-13]), kernel=Coord(2, ()))
    expected = np.diag([1.0, 1.0 / 5e-13])
    _require(np.asarray(m.matrix).tobytes() == expected.tobytes(), f"the metric face {np.asarray(m.matrix).tolist()} differs from {expected.tolist()}")
    _require(m.rank((2,)) == 2, f"rank {m.rank((2,))}, expected 2")


def test_both_doors_share_one_parse() -> None:
    """One bad input is one refusal on both doors, in every mode: the kernel and the shape are parsed before any decomposition."""
    cases = [
        (np.diag([1.0, 0.0, 0.0]), Coord(5, ()), ValueError, "kernel is declared on R"),
        (np.diag([1.0, 2.0]), 3, TypeError, "must be a Declaration"),
        (np.diag([1.0, 2.0]), None, TypeError, "must be a Declaration"),
        (np.diag([1.0, 2.0]), Ann(2, Spanned(2, ("k",), 1)), ValueError, "kernel is not declared by a spanned subspace"),
        (np.ones((2, 3)), Coord(2, ()), ValueError, "square 2-D"),
    ]
    for matrix, kernel, exc, msg in cases:
        with pytest.raises(exc, match=msg):
            DenseMetric(matrix, kernel=kernel)  # type: ignore[arg-type]
        with pytest.raises(exc, match=msg):
            DenseMetric.inverse_of(matrix, kernel=kernel)  # type: ignore[arg-type]


def test_the_frame_refuses_the_zero_gram_in_its_own_words() -> None:
    with pytest.raises(ValueError, match="frame Gram has no value above"):
        FF.frame_gram_kernel(np.zeros((3, 3)))


def test_inverse_of_runs_two_decompositions_per_mint(monkeypatch: pytest.MonkeyPatch) -> None:
    """3.3: ``inverse_of`` decomposes the Gram once, and the bare door's second face is not recomputed (a spy on ``np.linalg.svd``)."""
    calls: list[int] = []
    real = np.linalg.svd

    def spy(*a, **k):
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(np.linalg, "svd", spy)
    g = _gram([1.0, 0.5, 0.3, 0.2, 0.1])
    DenseMetric.inverse_of(g, kernel=Coord(5, ()))
    _require(len(calls) == 1, f"inverse_of made {len(calls)} decompositions of the Gram, expected 1")
    calls.clear()
    DenseMetric(g, kernel=Coord(5, ()))
    _require(len(calls) == 1, f"the bare door made {len(calls)} decompositions, expected 1")
