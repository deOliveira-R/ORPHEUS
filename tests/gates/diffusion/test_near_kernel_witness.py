r"""U3b rows D-W1/D-W2: the near-kernel witness against diffusion's predicted kernel and the dense SVD of its loss.

Spec ``scratch/membership_covector/u3b/gate_spec.md`` §2.2. The family covers the absorbing-chain predicate's shapes:
a one-group class, a two-group class, a class with a feeder (``g2_only``: the left co-vector reads the feeder's
absorbing-chain probability, not the flat unit), and two regular posings (up-scatter closes no class; a white wall at
0.999 leaks). The witness's dense route reads the flattened loss the solver factorises; the oracle is numpy's SVD of
the same loss assembled column by column through its matvec.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.diffusion.solver import DiffusionSolver
from orpheus.numerics.near_kernel import NearKernel
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.diffusion.test_fixed_source_verdict import _columns, _dense, _dm
from tests.gates.numerics._near_kernel_oracles import dense_kernels, residual
from tests.gates.transport._singular_loss_fixtures import BOTH_2G, G2_ONLY, G2_UPSCATTER, PURE_1G

pytestmark = [pytest.mark.l1, pytest.mark.rests_on("tests/gates/numerics/test_near_kernel.py")]

_DIFFUSION = [("pure1g", dict(mixtures={0: PURE_1G}), 1), ("both2g", dict(mixtures={0: BOTH_2G}), 1),
              ("g2_only", dict(mixtures={0: G2_ONLY}), 1), ("upscatter", dict(mixtures={0: G2_UPSCATTER}), 0),
              ("white0.999", dict(mixtures={0: PURE_1G}, left=0, right=0, bcs={0: WallCondition.white(0.999)}), 0)]


@pytest.mark.parametrize("label,kw,dim", _DIFFUSION, ids=[r[0] for r in _DIFFUSION])
def test_DW_the_witness_agrees_with_the_prediction_and_the_dense_svd(label, kw, dim):
    """Measured dim == dense dim == predicted dim; with a kernel, the predicted right span contains the measured one, the
    predicted co-vectors' Riesz representatives contain the measured left vectors, and the gap is >= 1e10 (`[M]` the
    test-architect's p2: 2.87e15, 6.54e14, 1.71e15); without one, the measured sigma_min/sigma_max is the dense one (`[M]`
    9.86e-3 and 2.37e-5)."""
    solver = DiffusionSolver(_dm(**kw))
    A = _dense(solver)
    g = np.asarray(solver.predicted_kernel.space.metric_diagonal(), dtype=float).ravel()
    sv, R, LH, k = dense_kernels(A, g)
    w = NearKernel.dense(solver.pencil.lhs)
    kernel = solver.predicted_kernel
    assert w.dimension == k == dim == kernel.dimension
    if dim:
        n = A.shape[0]
        assert residual(_columns(kernel.right, n), w.right) <= 1e-12
        assert residual(_columns(kernel.H.right, n), w.left) <= 1e-12
        assert w.singular_gap is not None and w.singular_gap >= 1e10
    else:
        assert w.sigma_min is not None
        np.testing.assert_allclose(w.sigma_min / w.sigma_max, sv[-1] / sv[0], rtol=1e-12)
