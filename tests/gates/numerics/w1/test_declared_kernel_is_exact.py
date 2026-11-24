r"""G01d: the declared kernel of each harness Gram is the exact kernel of the matrix (N01, spec 7.2, 3.3).

The form constructor does not test ``G K = 0`` in floating point (3.3); this gate does, with its tolerance in the gate. For each
of the 14 deficient Grams of the harness group (c), with ``w, v = np.linalg.eigh(g_sym)`` (ascending, so the kernel is the FIRST
``kernel.dim`` columns): ``max |G K| <= 1e-14 * norm(G)`` on the 13 frame-built or exactly deficient ones and
``<= FRAME_GRAM_RANK_RTOL * norm(G)`` on the synthetic near-cutoff Gram (its dropped value is 5e-13 by construction, so it is the
one member of the second tolerance); ``w[0] >= -1e-14 * norm(G)`` on all 21 (no negative direction). Wrong-declaration arms:
the declared dimension must separate the singular spectrum, ``s[n - dim] < 1e-6 s_max`` (what a declaration inverts) and
``s[n - dim - 1] > 1e-6 s_max`` (what it zeroes); the arms ``dim - 1`` (14 of 14 deficient) and ``dim + 1`` (21 of 21) must read
red. The indefinite control ``diag([1, -1])`` reads red on the smallest-eigenvalue predicate.

The declared kernel is read LIVE through ``frame_gram_kernel``, so a production mutation of the frame's count reddens these rows.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.frame_faces import FRAME_GRAM_RANK_RTOL
from tests.gates.numerics.w1 import _stored_grams as sg

pytestmark = pytest.mark.foundation

_TIGHT = 1e-14
_GAP = 1e-6
_NEAR_CUTOFF = "syn.near_cutoff_5x5"


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


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
@pytest.mark.rests_on("tests/gates/numerics/w1/test_frame_tolerance_has_one_address.py::test_the_frame_count_equals_the_pinv_nullity_across_the_cutoff")
def test_declared_kernel_is_exact_for_the_frames(family: str) -> None:
    e = sg.stored(family)
    g = e.g_sym
    dim = sg.declared_for(e).dim
    if dim:
        tol = FRAME_GRAM_RANK_RTOL if family == _NEAR_CUTOFF else _TIGHT
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
