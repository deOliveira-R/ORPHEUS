r"""G01g: the subject of #634: the inverse face is derived, never supplied, never admitted by a tolerance (N01, spec 7.2).

Defect 1 (too tight): ``DenseMetric.inverse_of`` REFUSED a legal Gram whose smallest eigenvalue is 2e-12 or 1e-10 (residuals 3.5e+05
and 9.9e+02 of the retired Penrose admission). Defect 2 (too loose): a supplied inverse ``X = G+ + v w^T`` was ADMITTED
(``|P - P^T| = 1.0``). After N01 the first Gram builds with rank 5 and the second spelling does not exist: there is no public
supplied-inverse keyword, and the field has ``init=False`` so ``dataclasses.replace`` cannot supply it either.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.declaration import Coord
from orpheus.numerics.metric import DenseMetric

pytestmark = pytest.mark.foundation

# tests/gates/numerics/test_quotient_metric.py:93 and :96 (copied: a gate owns its inputs).
_Q5 = np.linalg.qr(np.random.default_rng(0).standard_normal((5, 5)))[0]


def _gram(eigs):
    G = _Q5 @ np.diag(eigs) @ _Q5.T
    return (G + G.T) / 2


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


@pytest.mark.parametrize("smallest", [2e-12, 1e-10], ids=["2e-12", "1e-10"])
@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_rank_is_declared")
def test_a_legal_gram_with_a_small_eigenvalue_builds(smallest: float) -> None:
    g = _gram([1.0, 0.5, 0.3, 0.2, smallest])
    m = DenseMetric.inverse_of(g, kernel=Coord(5, ()))
    _require(m.rank((5,)) == 5, f"rank {m.rank((5,))}, expected 5 for a declared full-rank Gram of smallest eigenvalue {smallest:g}")
    g_sym = (g + g.T) / 2.0
    _require(np.asarray(m.inverse_matrix).tobytes() == g_sym.tobytes(), "the gram face of inverse_of is the exact symmetrisation of the Gram")
    # the faces are inverses to the conditioning of the Gram: a legal Gram is inverted, not refused and not zeroed.
    residual = float(np.linalg.norm(np.asarray(m.matrix) @ g_sym - np.eye(5)))
    _require(residual <= 1e-3, f"|G+ G - I| = {residual:.3e}: the face is not the inverse")


@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_rank_is_declared")
def test_no_public_supplied_inverse() -> None:
    g = _gram([1.0, 0.5, 0.3, 0.2, 0.1])
    gp = np.linalg.pinv(g, hermitian=True)
    rng = np.random.default_rng(3)
    x = gp + np.outer(rng.standard_normal(5), rng.standard_normal(5))  # not symmetric: P = X G is not an orthogonal projector
    with pytest.raises(TypeError):
        DenseMetric(g, kernel=Coord(5, ()), inverse_matrix=x)  # type: ignore[call-arg]
    m = DenseMetric(g, kernel=Coord(5, ()))
    with pytest.raises(TypeError):
        dataclasses.replace(m, inverse_matrix=x)
    # control: replacing the matrix re-derives the second face (the field is derived data).
    again = dataclasses.replace(m, matrix=2.0 * g)
    _require(np.asarray(again.inverse_matrix).tobytes() != np.asarray(m.inverse_matrix).tobytes(), "replace(matrix=) must rebuild the second face")
