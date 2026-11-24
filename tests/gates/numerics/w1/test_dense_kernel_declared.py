r"""G01c: the rank of a dense form is declared, never read from a cutoff (N01, spec 7.2, 3.3).

``DenseMetric(matrix, *, kernel)``: ``kernel`` is required and is ``Coord(n, ())`` (full rank) or a ``Spanned``. The rank on
elements of shape ``shape`` is ``(n - kernel.dim) * (size // n)``, ``size`` the product of ``shape``; the kernel basis has
``kernel.dim`` columns (the signature ``(n - kernel.dim, kernel.dim, 0)`` is declared, and read here through these two counts).
The eigenvalue ``5e-13`` of ``diag([1, 5e-13])`` is BELOW the retired cutoff ``1e-12`` and is a kept direction when the form
declares full rank: the "too tight" defect of #634 reading it as dropped is the red the last row catches.

First red today: ``DenseMetric(np.diag([1.0, 5e-13])).rank((2,))`` reads 1 (the cutoff), and the ``kernel`` keyword does not exist.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.declaration import Ann, Coord, Spanned
from orpheus.numerics.metric import DenseMetric
from tests.gates._dense import declared_kernel

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
