r"""G01h: the Moore-Penrose identities on the 21 harness Grams (N01, spec 7.2).

The replacement for the retired Penrose admission, which was a runtime check and left no gate (structure check L2). For each
stored Gram ``G`` (symmetrised) and ``G+ = DenseMetric(G, kernel=k).inverse_matrix`` built LIVE by the bare path:
``G G+ G = G``, ``G+ G G+ = G+`` and ``G+ G`` symmetric, each at relative Frobenius ``<= 1e-10`` (the maxima over the 21 are of
order 1e-12, ``syn.kappa1e4_8x8`` and ``syn.near_cutoff_5x5`` the largest). A wrong declaration reads red: ``dim + 1`` on 21 of 21
(a kept direction is zeroed), ``dim - 1`` on the 14 deficient Grams (a null direction is inverted: ``inf`` or a huge value).
A refusal of the wrong construction by the constructor under ``__debug__`` counts as read red.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.declaration import Declaration
from orpheus.numerics.metric import DenseMetric
from tests.gates._dense import declared_kernel
from tests.gates.numerics.w1 import _stored_grams as sg

pytestmark = pytest.mark.foundation

_TOL = 1e-10


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


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
@pytest.mark.rests_on("tests/gates/numerics/w1/test_declared_kernel_is_exact.py::test_declared_kernel_is_exact_for_the_frames")
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
