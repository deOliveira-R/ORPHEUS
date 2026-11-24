r"""G01g and G01h. **G01g.** The subject of #634: the inverse face is derived, never supplied, never admitted by a tolerance (N01, spec 7.2).

Defect 1 (too tight): ``DenseMetric.inverse_of`` REFUSED a legal Gram whose smallest eigenvalue is 2e-12 or 1e-10 (residuals 3.5e+05
and 9.9e+02 of the retired Penrose admission). Defect 2 (too loose): a supplied inverse ``X = G+ + v w^T`` was ADMITTED
(``|P - P^T| = 1.0``). After N01 the first Gram builds with rank 5 and the second spelling does not exist: there is no public
supplied-inverse keyword, and the field has ``init=False`` so ``dataclasses.replace`` cannot supply it either.

**G01h.** The Moore-Penrose identities on the 21 harness Grams (N01, spec 7.2).

The replacement for the retired Penrose admission, which was a runtime check and left no gate (structure check L2). For each
stored Gram ``G`` (symmetrised) and ``G+ = DenseMetric(G, kernel=k).inverse_matrix`` built LIVE by the bare path:
``G G+ G = G``, ``G+ G G+ = G+`` and ``G+ G`` symmetric, each at relative Frobenius ``<= 1e-10`` (the maxima over the 21 are of
order 1e-12, ``syn.kappa1e4_8x8`` and ``syn.near_cutoff_5x5`` the largest). A wrong declaration reads red: ``dim + 1`` on 21 of 21
(a kept direction is zeroed), ``dim - 1`` on the 14 deficient Grams (a null direction is inverted: ``inf`` or a huge value).
A refusal of the wrong construction by the constructor under ``__debug__`` counts as read red.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.declaration import Coord, Declaration
from orpheus.numerics.metric import DenseMetric
from tests.gates._dense import declared_kernel
from tests.gates.numerics.w1 import _stored_grams as sg

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
    """A supplied inverse ``X = G+ + k w^T`` (``k`` in the kernel of ``G``) satisfied the retired admission ``G X G = G`` and is not the
    Moore-Penrose inverse (``X G`` is not symmetric): it is unspellable, as a keyword and through ``dataclasses.replace``."""
    g = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])  # rank 2, kernel spanned by (1, -1, 0)
    kernel = declared_kernel(g, 1)
    k = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    x = np.linalg.pinv(g, hermitian=True) + np.outer(k, np.array([0.0, 0.0, 1.0]))
    _require(float(np.max(np.abs(g @ x @ g - g))) <= 1e-14, "activation: X must satisfy the retired admission G X G = G")
    _require(float(np.max(np.abs(x @ g - (x @ g).T))) > 0.5, "activation: X G must be far from symmetric (X is not the Moore-Penrose inverse)")
    with pytest.raises(TypeError):
        DenseMetric(g, kernel=kernel, inverse_matrix=x)  # type: ignore[call-arg]
    m = DenseMetric(g, kernel=kernel)
    with pytest.raises(TypeError):
        dataclasses.replace(m, inverse_matrix=x)
    # control: replacing the matrix re-derives the second face (the field is derived data).
    again = dataclasses.replace(m, matrix=2.0 * g)
    _require(np.asarray(again.inverse_matrix).tobytes() != np.asarray(m.inverse_matrix).tobytes(), "replace(matrix=) must rebuild the second face")


# ============================================================
# G01h
# ============================================================

_TOL = 1e-10


def penrose_residuals(g: np.ndarray, gp: np.ndarray) -> tuple[float, float, float]:
    """Relative Frobenius residuals of ``G G+ G = G``, ``G+ G G+ = G+`` and the symmetry of ``G+ G``."""
    p = gp @ g
    return (
        float(np.linalg.norm(g @ gp @ g - g) / np.linalg.norm(g)),
        float(np.linalg.norm(gp @ g @ gp - gp) / np.linalg.norm(gp)),
        float(np.linalg.norm(p - p.T) / np.linalg.norm(p)),
    )


def reads_red(g: np.ndarray, kernel: Declaration) -> bool:
    """True iff the declaration is refused or its inverse breaks an identity at ``_TOL`` (NaN reads red)."""
    try:
        gp = np.asarray(DenseMetric(g, kernel=kernel).inverse_matrix)
    except ValueError:
        return True
    with np.errstate(all="ignore"):
        return not all(r <= _TOL for r in penrose_residuals(g, gp))


@pytest.mark.parametrize("family", sg.families())
@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_declared_kernel_is_exact_for_the_frames")
@pytest.mark.rests_on("tests/gates/cross_method/w1/test_dense_inverse_bytes.py::test_bare_path_equals_the_stored_metric")
def test_moore_penrose_identities_on_the_harness(family: str) -> None:
    e = sg.stored(family)
    g = e.g_sym
    gp = np.asarray(DenseMetric(g, kernel=sg.declared_for(e)).inverse_matrix)
    r1, r2, r3 = penrose_residuals(g, gp)
    _require(r1 <= _TOL, f"[{family}] G G+ G = G: relative residual {r1:.3e}")
    _require(r2 <= _TOL, f"[{family}] G+ G G+ = G+: relative residual {r2:.3e}")
    _require(r3 <= _TOL, f"[{family}] G+ G symmetric: relative asymmetry {r3:.3e}")


def test_a_wrong_declaration_reads_red() -> None:
    fams = sg.families()
    dims = {f: sg.declared_for(sg.stored(f)).dim for f in fams}
    up = [f for f in fams if reads_red(sg.stored(f).g_sym, declared_kernel(sg.stored(f).g_sym, dims[f] + 1))]
    down = [f for f in fams if dims[f] > 0 and reads_red(sg.stored(f).g_sym, declared_kernel(sg.stored(f).g_sym, dims[f] - 1))]
    _require(len(up) == 21, f"the arm dim + 1 read red on {len(up)} of 21")
    _require(len(down) == sum(1 for f in fams if dims[f] > 0) == 14, f"the arm dim - 1 read red on {len(down)} of 14 deficient Grams")
