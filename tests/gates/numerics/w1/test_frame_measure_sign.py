r"""G03r: a frame refuses an under-resolved basis, the Gram of which is indefinite on a measure with a negative weight (N03 R1', re-ruled by N03b, spec 7.4b).

In this gate's own process, run under ``python -O`` (so that a bare ``assert`` could not carry it and a refusal under ``if __debug__:`` is not read):
``angular_frame(Quadrature.lebedev(13), 7)`` and ``lebedev(27), 14`` raise ``ValueError`` at the frame's mint (``basis_space``, which reads
``FromFrame.form``), the message naming the Krein arm, #706 and the count of negative kept directions; ``lebedev(13)`` at L = 6, ``lebedev(11)`` at L = 7 and
``lebedev(25)`` at L = 13 build and ``basis_space`` reads their forms without a refusal; ``Axis.of_quadrature(Quadrature.lebedev(13))`` builds (the same rule as
a diagonal axis, signature ``(66, 0, 8)``); ``numerics/frame.py`` carries the ``SCOPE-BOUNDARY[guard]`` tag with its three fields (the ledger's own predicate)
and names #706, and has ONE ``svd`` call (the private body that counts the rank and reads the signs).

The count in the message is compared with an independent count from ``np.linalg.eigvalsh`` of the symmetrised Gram, at the threshold ``1e-9`` times
the largest eigenvalue; ``test_the_refusal_margin_is_far_from_every_threshold`` asserts the spectral gap that makes the count exact on any LAPACK: no eigenvalue has
a magnitude between ``1e-13`` and ``1e-6`` of the largest, and the smallest negative and the smallest positive kept direction exceed ``1e-6``
of it (measured on this host, ``[M]`` this gate's helper: 1.5e-2 and 6.7e-2 for ``lebedev(13)`` at L = 7, 1.5e-3 and 1.8e-2 for ``lebedev(27)`` at L = 14). First red at ``ff305e0f``: ``lebedev(13)`` at L = 7 builds under
``-O`` and has 4 negative eigenvalues; the Parseval mint does not refuse. Mutations: drop the refusal (reads ``lebedev(13) L=7 reads accepted``); move it
under ``if __debug__:`` (the same line, read in the ``-O`` process).
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _mint(order: int, L: int) -> tuple[str, str]:
    """``(verdict, message)``: the frame's mint through ``basis_space``."""
    from orpheus.numerics.frame_faces import angular_frame, basis_space
    from orpheus.numerics.quadrature.directional import Quadrature

    try:
        basis_space(angular_frame(Quadrature.lebedev(order), L))
    except ValueError as e:
        return "refused", str(e)
    return "accepted", ""


def _independent_negative_count(order: int, L: int) -> tuple[int, float, float]:
    """``(negative kept directions, smallest |negative| / s_max, smallest positive-kept / s_max)`` by ``eigvalsh`` of the symmetrised Gram."""
    from orpheus.numerics.frame_faces import angular_frame
    from orpheus.numerics.quadrature.directional import Quadrature

    g = np.asarray(angular_frame(Quadrature.lebedev(order), L).discrete_gram, dtype=float)
    ev = np.linalg.eigvalsh((g + g.T) / 2)
    smax = float(np.max(np.abs(ev)))
    neg = ev[ev < -1e-9 * smax]
    pos = ev[ev > 1e-9 * smax]
    return int(neg.size), float(np.min(np.abs(neg))) / smax if neg.size else float("nan"), float(np.min(pos)) / smax


def test_a_frame_refuses_an_indefinite_gram() -> None:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.quadrature.directional import Quadrature

    _require(not __debug__, "this gate must run under python -O so that only a real raise can carry the refusal")
    for order, L in ((13, 7), (27, 14)):
        verdict, msg = _mint(order, L)
        _require(verdict == "refused", f"lebedev({order}) L={L} reads {verdict}")
        n_neg, _, _ = _independent_negative_count(order, L)
        _require(n_neg > 0, f"lebedev({order}) L={L}: the independent count finds no negative direction")
        _require("Krein" in msg and "#706" in msg, f"lebedev({order}) L={L}: the message does not name the Krein arm and #706: {msg}")
        _require(re.search(rf"\b{n_neg}\b", msg) is not None, f"lebedev({order}) L={L}: the message {msg!r} does not carry the count {n_neg}")
    for order, L in ((13, 6), (11, 7), (25, 13)):
        verdict, msg = _mint(order, L)
        _require(verdict == "accepted", f"lebedev({order}) L={L} reads {verdict}: {msg}")
    _require(Axis.of_quadrature(Quadrature.lebedev(13)).form.signature == (66, 0, 8), "Axis.of_quadrature(lebedev(13)) does not read signature (66, 0, 8)")


def test_the_refusal_margin_is_far_from_every_threshold() -> None:
    """The count is exact on any LAPACK only if the negative directions sit far from zero and the kept ones far from the rank cutoff."""
    from orpheus.numerics.frame_faces import angular_frame
    from orpheus.numerics.quadrature.directional import Quadrature

    for order, L in ((13, 7), (27, 14)):
        g = np.asarray(angular_frame(Quadrature.lebedev(order), L).discrete_gram, dtype=float)
        ev = np.abs(np.linalg.eigvalsh((g + g.T) / 2))
        rel = ev / float(np.max(ev))
        gap = rel[(rel > 1e-13) & (rel < 1e-6)]
        _require(gap.size == 0, f"lebedev({order}) L={L}: {gap.size} eigenvalues sit between 1e-13 and 1e-6 of the largest: the rank cutoff would decide them")
        n_neg, neg_margin, pos_margin = _independent_negative_count(order, L)
        _require(n_neg > 0 and neg_margin > 1e-6, f"lebedev({order}) L={L}: the smallest negative direction is {neg_margin:.3e} of the largest: the count is host-dependent")
        _require(pos_margin > 1e-6, f"lebedev({order}) L={L}: the smallest positive direction is {pos_margin:.3e} of the largest")


def test_the_site_is_tagged_and_reads_one_decomposition() -> None:
    from orpheus.numerics import frame as FR
    from tests.gates import test_elegance_debt_is_tagged as ledger

    path = Path(FR.__file__)
    src = path.read_text()
    _require("SCOPE-BOUNDARY[guard]" in src, "numerics/frame.py carries no SCOPE-BOUNDARY[guard]")
    _require(not ledger._boundary_violations_in(src, label="frame.py"), "a boundary tag of frame.py lacks machinery, ruling or revisit")
    _require("706" in src, "numerics/frame.py does not name #706")
    n_svd = sum(1 for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "svd")
    _require(n_svd == 1, f"numerics/frame.py has {n_svd} svd calls, expected the one decomposition that counts the rank and reads the signs")
    _require(hasattr(FR, "frame_gram_kernel") and hasattr(FR, "frame_kernel"), "frame.py lacks the two verbs frame_gram_kernel and frame_kernel")


def test_a_negative_gram_diagonal_takes_the_dense_verdict_and_is_refused_by_the_kernel_door() -> None:
    """The diagonal verdict needs no refusal of its own: ``FrameBase.discrete_gram_structure`` returns DENSE when any entry of the Gram's diagonal is
    negative (its first test, ``np.any(diag < 0.0)``), so no frame with a negative diagonal reaches the diagonal arm of ``FromFrame.form``. Real
    ``GalerkinFrame`` instances. Leg one: ``lebedev(13)`` at L = 7 to 15 and ``lebedev(27)`` at L = 14, 15 have a negative Gram diagonal entry (11 of the
    48 frames of ``lebedev`` 13, 25, 27 at L = 0 to 15 with one; their off-diagonals are large, so the rank rule alone would also route them DENSE). Leg two
    isolates the diagonal test: the negated ``lebedev(13)`` rule at L = 0 has the one-mode Gram ``[-4 pi]``, no off-diagonal at all, so only the negative
    diagonal entry routes it DENSE. Every such frame is refused by ``frame_kernel`` naming the Krein arm and #706. Mutation: drop the diagonal test of the
    verdict (leg two reads DIAGONAL and reddens)."""
    from dataclasses import replace

    from orpheus.numerics.basis.base import GramStructure
    from orpheus.numerics.frame import frame_kernel
    from orpheus.numerics.frame_faces import angular_frame
    from orpheus.numerics.quadrature.directional import Quadrature

    negated = Quadrature(measure=replace(Quadrature.lebedev(13).measure, weights=-np.asarray(Quadrature.lebedev(13).measure.weights)))
    cases = [(f"lebedev(13) L={L}", Quadrature.lebedev(13), L) for L in range(7, 16)]
    cases += [(f"lebedev(27) L={L}", Quadrature.lebedev(27), L) for L in (14, 15)]
    cases += [("negated lebedev(13) L=0", negated, 0)]
    for label, rule, L in cases:
        frame = angular_frame(rule, L)
        _require(bool(np.any(np.diagonal(frame.discrete_gram) < 0.0)), f"{label}: no negative Gram diagonal entry")
        _require(frame.discrete_gram_structure is GramStructure.DENSE, f"{label}: a negative Gram diagonal reads {frame.discrete_gram_structure}, not DENSE")
        try:
            frame_kernel(frame)
        except ValueError as e:
            _require("Krein" in str(e) and "#706" in str(e), f"{label}: the refusal does not name the Krein arm and #706: {e}")
        else:
            pytest.fail(f"{label}: frame_kernel accepted a Gram with a negative kept direction")
