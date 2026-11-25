r"""G03m: the form of a frame head's factor is today's entry (N03, spec 7.4; re-pointed at N04 from the reader to ``Factor.form``).

G03m. On five frames built by ``angular_frame`` (numerics only, verdicts read): DIAGONAL ``level_symmetric(4)`` L = 2 and ``folded_product(2, 4)`` L = 1;
DENSE ``folded_product(2, 4)`` L = 2, ``product(4, 4)`` L = 2 and ``level_symmetric(4)`` L = 3. ``FromFrame.form`` is ``Diag(diag G)`` on the diagonal
verdict and ``DenseMetric(symmetrised G, kernel, face=METRIC)`` on the dense one, byte-equal to the frame's ``discrete_gram`` symmetrised; the reader's
entry of a diagonal head equals a hand-written masked reciprocal ``np.where(d > 0, 1 / np.where(d > 0, d, 1), 0)`` of the diagonal byte for byte (an
independent expression, not the production body); a dense head gives no entry and its form with the face flipped equals the Gram stated with ``face=GRAM``
in both faces and in content digest; ``origin.frame is frame``.

(G03n, the reader's counting and ``NoForm`` entries, retired with ``_axis_entries`` at N04: a counting form is a form that is applied, a ``NoForm`` slot is skipped by the walk, gated by G04e and G04f. The refusal of a ``Declared`` dense axis, ELEGANCE-DEBT #755, retired with the overlay: a ``Declared`` axis holding a ``DenseMetric`` is a factor on a dense form like any other.)

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
    from orpheus.numerics.frame_faces import basis_space
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, Face
    from orpheus.values.content import content_digest

    for factory, args, L in _DIAGONAL:
        frame, ax = _frame_head(factory, args, L)
        g = DenseMetric.symmetrised(np.asarray(frame.discrete_gram, dtype=float))
        _require(isinstance(ax.form, DiagonalMetric), f"{factory}{args} L={L}: the form is {type(ax.form).__name__}, expected the diagonal verdict's Diag(diag G)")
        _require(ax.form.weights.tobytes() == np.ascontiguousarray(np.diag(g)).tobytes(), f"{factory}{args} L={L}: the form is not diag G")
        entry = basis_space(frame).factors[0].form
        _require(isinstance(entry, DiagonalMetric), f"{factory}{args} L={L}: the head factor's form is {entry!r}")
        want = np.ascontiguousarray(_masked_reciprocal(np.diag(g)))
        _require(np.asarray(entry.weights).tobytes() == want.tobytes(), f"{factory}{args} L={L}: the form is not the masked reciprocal of diag G")
        _require(ax.origin.frame is frame, f"{factory}{args} L={L}: origin.frame is not the frame")
    for factory, args, L in _DENSE:
        frame, ax = _frame_head(factory, args, L)
        g = DenseMetric.symmetrised(np.asarray(frame.discrete_gram, dtype=float))
        kernel = FR.frame_gram_kernel(np.asarray(frame.discrete_gram, dtype=float))
        _require(isinstance(ax.form, DenseMetric), f"{factory}{args} L={L}: the form is {type(ax.form).__name__}, expected the dense verdict's DenseMetric")
        _require(ax.form.face is Face.METRIC and ax.form.stated.tobytes() == g.tobytes(), f"{factory}{args} L={L}: the axis form is not the symmetrised Gram stated as the metric")
        _require(ax.form.kernel == kernel, f"{factory}{args} L={L}: the declared kernel differs")
        _require(isinstance(basis_space(frame).factors[0].form, DenseMetric), f"{factory}{args} L={L}: a dense head's factor has no dense form")
        head = DenseMetric(g, kernel=kernel, face=Face.GRAM)
        flipped = ax.form.flipped()
        _require(
            flipped.matrix.tobytes() == head.matrix.tobytes() and flipped.inverse_matrix.tobytes() == head.inverse_matrix.tobytes(),
            f"{factory}{args} L={L}: the axis form flipped is not the head's integrated metric in both faces",
        )
        _require(content_digest(flipped) == content_digest(head), f"{factory}{args} L={L}: the flipped axis form and the head have different content digests")
        _require(ax.origin.frame is frame, f"{factory}{args} L={L}: origin.frame is not the frame")


def test_a_diagonal_frame_head_factor_reads_the_space_metric_bytes() -> None:
    """The factor's integrated form on a diagonal frame head (the masked reciprocal of the head's diagonal Gram, ELEGANCE-DEBT #756) is a
    ``DiagonalMetric`` whose apply is byte-equal to the metric the dressed space applies (``basis_space(frame).apply_metric``), on
    each of the two diagonal-verdict frames; the dual of the density form is NOT that metric on a dead slot (it divides by zero there). Exact on any host
    (``.tobytes()`` of two reads of one arithmetic, no decomposition). Mutation: the factor's branch replaced by ``density.dual()`` reddens
    the type leg on the first frame (a ``DualMetric``)."""
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.factor import Factor, Kind, Subquotient, Variance, VarianceChart
    from orpheus.numerics.frame_faces import basis_space
    from orpheus.numerics.metric import DiagonalMetric

    for factory, args, L in _DIAGONAL:
        frame, ax = _frame_head(factory, args, L)
        n = int(np.prod(ax.shape))
        dead = ax.form.kernel  # the head's dead slots: a density factor on the whole space must quotient by them
        density = Factor(ax, Variance(Kind.MEASURE, VarianceChart.DENSITY), Subquotient(Coord(n, tuple(range(n))), dead))
        form = density.riesz().form  # the integrated end: the form induced from its density partner
        _require(isinstance(form, DiagonalMetric), f"{factory}{args} L={L}: the integrated form is {type(form).__name__}")
        space = basis_space(frame)
        x = np.random.default_rng(7).standard_normal(space.shape)
        _require(
            np.asarray(form.apply(x)).tobytes() == np.asarray(space.apply_metric(x)).tobytes(),
            f"{factory}{args} L={L}: the factor's integrated form does not apply byte-equal to the space's metric",
        )


def test_a_declared_dense_axis_is_a_factor_on_a_dense_form() -> None:
    """A ``Declared`` axis holding a ``DenseMetric`` is a factor like any other: the space pairs by the Gram (the pairing of ``x = (1, 2, 3)`` with
    ``y = (1, -1, 2)`` reads ``x^T G y`` by hand, 22.0, where the Euclidean form reads 5.0). Control: a dense frame head's factor holds a dense form.
    Exact on any host up to the last bit of a 3-term sum, compared within the relative tolerance of one hand product."""
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.factor import Kind, Variance, VarianceChart
    from orpheus.numerics.metric import DenseMetric, Face
    from orpheus.numerics.space import TensorSpace

    gram = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]])
    axis = Axis("d", (3,), origin=Declared(BasisKind.NODAL, DenseMetric(gram, kernel=Coord(3, ()), face=Face.METRIC)))
    space = TensorSpace.declared("d", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=axis.form)
    x, y = np.array([1.0, 2.0, 3.0]), np.array([1.0, -1.0, 2.0])
    _require(abs(space.inner_product(x, y) - float(x @ gram @ y)) <= 1e-12 * 22.0, f"the pairing is {space.inner_product(x, y)}, by hand {float(x @ gram @ y)}")
    _, head = _frame_head("folded_product", (2, 4), 2)
    _require(isinstance(head.form, DenseMetric), "the control head is not a dense form")
