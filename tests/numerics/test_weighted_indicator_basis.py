r"""Intrinsic-property tests for :class:`WeightedIndicatorBasis` — the PG test side.

The weighted cell indicator is the **test** basis of the spatial homogenisation
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`: its analysis functional is the
:math:`w`-reweighted region integral :math:`\langle\chi_R, f\rangle_W = \int_R w f\,
\mathrm{d}V`, with the geometric support coming from the wrapped trial
:class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`. These tests pin its
defining laws (the project standard: every math-bearing type ships a test of its
defining property):

* :meth:`evaluate` is the **plain membership** table — the weight is an *analysis*
  weight, NOT a tabulation (so a per-group flux rides ``analyze``'s trailing axes
  instead of inflating the table);
* :meth:`analyze` folds the weight into the region integral;
* the **degenerate** ``w ≡ 1`` reduces ``analyze`` bit-identically to the plain
  indicator (the basis-level analogue of the frame-level PG→Galerkin degenerate);
* the synthesis side **raises** (a test-only basis — never a half-consistent one).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis import IndicatorBasis, WeightedIndicatorBasis

pytestmark = [pytest.mark.foundation]


# 4 fine nodes on [0, 4], 2 coarse cells [0,2],[2,4]; non-uniform measure weights.
_EDGES = (np.array([0.0, 2.0, 4.0]),)
_NODES = np.array([0.5, 1.5, 2.5, 3.5])
_WEIGHTS = np.array([0.5, 1.5, 1.0, 1.0])      # volume measure V_i


def _trial() -> IndicatorBasis:
    return IndicatorBasis(_EDGES)


def test_evaluate_is_the_plain_membership_table_weight_free():
    """``evaluate`` returns the geometric membership — the weight does NOT enter it.

    The weight is an analysis weight; the support (which node is in which cell) is
    weight-independent. A weight-in-evaluate would double-count it in every
    contraction.
    """
    trial = _trial()
    w = np.array([2.0, 3.0, 5.0, 7.0])
    wb = WeightedIndicatorBasis(trial, w)
    np.testing.assert_array_equal(wb.evaluate(_NODES), trial.evaluate(_NODES))


def test_analyze_folds_the_weight_into_the_region_integral():
    r"""``analyze(f) = \sum_{i\in R} w_i\,V_i\,f_i`` — pinned against the hand loop."""
    trial = _trial()
    w = np.array([2.0, 3.0, 5.0, 7.0])
    wb = WeightedIndicatorBasis(trial, w)
    table = wb.evaluate(_NODES)
    f = np.array([10.0, 20.0, 30.0, 40.0])

    out = wb.analyze(f, table, _WEIGHTS)
    regions = [[0, 1], [2, 3]]
    expected = np.array([sum(w[i] * _WEIGHTS[i] * f[i] for i in sel) for sel in regions])
    np.testing.assert_allclose(out, expected, rtol=1e-14)


def test_per_group_weight_aligns_to_the_source_group_of_a_matrix_channel():
    r"""A per-group weight ``w (n, ng)`` weights a ``[g_from, g_to]`` matrix by ``g_from``.

    The leading-aligned broadcast aligns the weight's group axis to the FIRST trailing
    (source) axis of the analysed field — the source-group weighting homogenisation's
    scattering channels require.
    """
    trial = _trial()
    ng = 2
    w = np.array([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0], [7.0, 8.0]])   # (n, ng)
    wb = WeightedIndicatorBasis(trial, w)
    table = wb.evaluate(_NODES)
    rng = np.random.default_rng(7)
    mat = rng.standard_normal((4, ng, ng))                          # (n, g_from, g_to)

    out = wb.analyze(mat, table, _WEIGHTS)                          # (n_R, g_from, g_to)
    regions = [[0, 1], [2, 3]]
    expected = np.zeros((2, ng, ng))
    for R, sel in enumerate(regions):
        for gf in range(ng):
            for gt in range(ng):
                expected[R, gf, gt] = sum(
                    w[i, gf] * _WEIGHTS[i] * mat[i, gf, gt] for i in sel
                )
    np.testing.assert_allclose(out, expected, rtol=1e-13)


def test_unit_weight_reduces_to_the_plain_indicator_analyze_bit_identically():
    """``w ≡ 1`` ⟹ ``analyze`` is BIT-IDENTICAL to the plain indicator's.

    The basis-level analogue of the frame-level PG→Galerkin degenerate: the weighted
    test reduces to its trial side exactly when the weight is unity (``array_equal``).
    """
    trial = _trial()
    wb = WeightedIndicatorBasis(trial, np.ones(_NODES.shape[0]))
    table = trial.evaluate(_NODES)
    rng = np.random.default_rng(9)
    f = rng.standard_normal((4, 3))
    np.testing.assert_array_equal(
        wb.analyze(f, table, _WEIGHTS), trial.analyze(f, table, _WEIGHTS),
    )


@pytest.mark.parametrize(
    "op",
    ["synthesize", "reconstruct", "reconstruct_transpose", "analyze_transpose"],
)
def test_synthesis_side_raises_test_only_basis(op):
    """The synthesis-side operations raise — a test-only basis, never half-consistent.

    A :class:`PetrovGalerkinFrame` takes its reconstruction from the TRIAL basis, so
    the weighted test's synthesis has no consumer. Rather than silently delegate to
    the *unweighted* indicator (``analyze`` folds ``w``, an unweighted synthesis would
    not — a landmine), the synthesis side raises until a real adjoint consumer builds
    the weighted transposes consistently.
    """
    wb = WeightedIndicatorBasis(_trial(), np.array([2.0, 3.0, 5.0, 7.0]))
    table = wb.evaluate(_NODES)
    with pytest.raises(NotImplementedError, match="no consumer"):
        if op == "synthesize":
            wb.synthesize(np.zeros((2,)), table)
        elif op == "reconstruct":
            wb.reconstruct(np.zeros((2,)), table)
        elif op == "reconstruct_transpose":
            wb.reconstruct_transpose(np.zeros((4,)), table)
        elif op == "analyze_transpose":
            wb.analyze_transpose(np.zeros((2,)), table, _WEIGHTS)


def test_mass_matrix_raises_test_only_basis():
    """``mass_matrix`` raises — the test-only synthesis side has no consumer."""
    from orpheus.numerics.measure import DiscreteMeasure

    wb = WeightedIndicatorBasis(_trial(), np.array([2.0, 3.0, 5.0, 7.0]))
    measure = DiscreteMeasure(nodes=_NODES, weights=_WEIGHTS, support="spatial_R1")
    with pytest.raises(NotImplementedError, match="no consumer"):
        wb.mass_matrix(measure)
