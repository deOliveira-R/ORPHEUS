r"""G03m and G03n: the one reader of the axis form in ``numerics/space.py`` gives today's entry (N03, spec 7.4).

G03m. On five frames built by ``angular_frame`` (numerics only, verdicts read): DIAGONAL ``level_symmetric(4)`` L = 2 and ``folded_product(2, 4)`` L = 1;
DENSE ``folded_product(2, 4)`` L = 2, ``product(4, 4)`` L = 2 and ``level_symmetric(4)`` L = 3. ``FromFrame.form`` is ``Diag(diag G)`` on the diagonal
verdict and ``DenseMetric(symmetrised G, kernel, face=METRIC)`` on the dense one, byte-equal to the frame's ``discrete_gram`` symmetrised; the reader's
entry of a diagonal head equals a hand-written masked reciprocal ``np.where(d > 0, 1 / np.where(d > 0, d, 1), 0)`` of the diagonal byte for byte (an
independent expression, not the production body); a dense head gives no entry and its form with the face flipped equals the Gram stated with ``face=GRAM``
in both faces and in content digest; ``origin.frame is frame``.

G03n. ``_axes_metric`` over a space of a counting axis and a ``NoForm`` axis (an ``EnergyAxis``): every entry is ``None`` and the metric is ``None``; a
weighted axis gives its ``DiagonalMetric`` entry. Green at the base tip (a refactor guard: the tip's ``_axes_metric`` returns ``None`` for an
``EnergyAxis`` and a counting axis); the counting and ``NoForm`` rules move no byte and no public observable, so the private structure is their only
witness, and it retires with ``_axis_entries`` at N04.

Exact on any host: byte comparisons of a computation with an independent expression, and structure. First red at ``ff305e0f`` (G03m): ``FromFrame`` has no
``form`` (``AttributeError``) and ``numerics.frame.frame_gram_kernel`` does not exist. Mutations (G03m): drop the reciprocal in the reader; return the dense
form as the entry; take the dense head's integrated form as the axis form's stated face with no flip. (G03n): the counting rule dropped (counting is an
entry); ``NoForm`` read as a counting entry.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation

_DIAGONAL = (("level_symmetric", (4,), 2), ("folded_product", (2, 4), 1))
_DENSE = (("folded_product", (2, 4), 2), ("product", (4, 4), 2), ("level_symmetric", (4,), 3))


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _frame_head(factory: str, args: tuple[int, ...], L: int):
    from orpheus.numerics.frame_faces import angular_frame, basis_space
    from orpheus.numerics.quadrature.directional import Quadrature

    frame = angular_frame(getattr(Quadrature, factory)(*args), L)
    return frame, basis_space(frame).axes[0]


def _masked_reciprocal(d: np.ndarray) -> np.ndarray:
    return np.where(d > 0, 1 / np.where(d > 0, d, 1), 0)


def test_from_frame_form_and_its_entry() -> None:
    from orpheus.numerics import frame as FR
    from orpheus.numerics import space as S
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, Face
    from orpheus.values.content import content_digest

    for factory, args, L in _DIAGONAL:
        frame, ax = _frame_head(factory, args, L)
        g = DenseMetric.symmetrised(np.asarray(frame.discrete_gram, dtype=float))
        _require(isinstance(ax.form, DiagonalMetric), f"{factory}{args} L={L}: the form is {type(ax.form).__name__}, expected the diagonal verdict's Diag(diag G)")
        _require(ax.form.weights.tobytes() == np.ascontiguousarray(np.diag(g)).tobytes(), f"{factory}{args} L={L}: the form is not diag G")
        entry = S._entry_of_axis(ax)
        _require(isinstance(entry, DiagonalMetric), f"{factory}{args} L={L}: the reader's entry is {entry!r}")
        want = np.ascontiguousarray(_masked_reciprocal(np.diag(g)))
        _require(np.asarray(entry.weights).tobytes() == want.tobytes(), f"{factory}{args} L={L}: the entry is not the masked reciprocal of diag G")
        _require(ax.origin.frame is frame, f"{factory}{args} L={L}: origin.frame is not the frame")
    for factory, args, L in _DENSE:
        frame, ax = _frame_head(factory, args, L)
        g = DenseMetric.symmetrised(np.asarray(frame.discrete_gram, dtype=float))
        kernel = FR.frame_gram_kernel(np.asarray(frame.discrete_gram, dtype=float))
        _require(isinstance(ax.form, DenseMetric), f"{factory}{args} L={L}: the form is {type(ax.form).__name__}, expected the dense verdict's DenseMetric")
        _require(ax.form.face is Face.METRIC and ax.form.stated.tobytes() == g.tobytes(), f"{factory}{args} L={L}: the axis form is not the symmetrised Gram stated as the metric")
        _require(ax.form.kernel == kernel, f"{factory}{args} L={L}: the declared kernel differs")
        _require(S._entry_of_axis(ax) is None, f"{factory}{args} L={L}: a dense head gives an entry")
        head = DenseMetric(g, kernel=kernel, face=Face.GRAM)
        flipped = ax.form.flipped()
        _require(
            flipped.matrix.tobytes() == head.matrix.tobytes() and flipped.inverse_matrix.tobytes() == head.inverse_matrix.tobytes(),
            f"{factory}{args} L={L}: the axis form flipped is not the head's integrated metric in both faces",
        )
        _require(content_digest(flipped) == content_digest(head), f"{factory}{args} L={L}: the flipped axis form and the head have different content digests")
        _require(ax.origin.frame is frame, f"{factory}{args} L={L}: origin.frame is not the frame")


def test_counting_and_noform_give_no_entry() -> None:
    from orpheus.numerics import space as S
    from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
    from orpheus.numerics.metric import DiagonalMetric
    from orpheus.numerics.quadrature.directional import Quadrature

    counting = Axis("s", (3,), origin=Declared(BasisKind.NODAL))
    energy = EnergyAxis.synthetic(2)
    entries = S._axis_entries((counting, energy))
    _require(len(entries) == 2 and all(e is None for _, e in entries), f"a counting axis and an energy axis give entries {entries}")
    _require(S._axes_metric((counting, energy)) is None, "the metric of a counting and an energy axis is not None")
    q = Quadrature.gauss_legendre(4)
    weighted = Axis.of_measure(q.measure, "x")
    wentries = S._axis_entries((weighted, energy))
    _require(isinstance(wentries[0][1], DiagonalMetric) and wentries[1][1] is None, f"a weighted axis gives {wentries}")
    _require(np.array_equal(np.asarray(wentries[0][1].weights), np.asarray(q.measure.weights)), "the weighted axis's entry is not its measure")
