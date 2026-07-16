"""Intrinsic tests for the ORPHEUS canonical energy-group convention.

The canonical convention is **group 0 = FASTEST** (highest energy): ``eg`` is
strictly DESCENDING and ``SigS[g_from, g_to]`` downscatter sits in the UPPER
triangle (``g_to > g_from``). NJOY/GENDF stores the opposite (group 0 = thermal,
ascending energy); the single normalisation at the data-ingest boundary
(:func:`~orpheus.data.micro_xs.gendf._to_canonical_group_order`, via
:func:`~orpheus.data.micro_xs.gendf._reverse_groups_2d` for the matrix channels)
flips it once, so every downstream consumer is order-transparent.

See ``docs/theory/foundations/cross_section_data.rst`` (canonical group convention).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.micro_xs import load_isotope
from orpheus.data.micro_xs.gendf import _reverse_groups_2d

pytestmark = pytest.mark.foundation


class TestReverseGroups2D:
    """The sparse both-axes group reversal — the matrix-channel primitive."""

    def test_moves_downscatter_lower_to_upper_triangle(self) -> None:
        # A lower-triangular downscatter matrix [g_from, g_to] (thermal-first):
        # row g_from scatters down to columns g_to <= g_from.
        lower = csr_matrix(
            np.array(
                [
                    [0.5, 0.0, 0.0, 0.0],
                    [0.3, 0.5, 0.0, 0.0],
                    [0.1, 0.3, 0.5, 0.0],
                    [0.0, 0.1, 0.3, 0.5],
                ]
            )
        )
        rev = np.asarray(_reverse_groups_2d(lower).todense())
        # After reversal the off-diagonal mass sits in the UPPER triangle.
        assert np.triu(rev, 1).sum() > 0.0
        assert np.tril(rev, -1).sum() == 0.0
        # In-group self-scatter (the diagonal) is order-invariant.
        assert np.allclose(np.diag(rev), [0.5, 0.5, 0.5, 0.5])

    def test_is_an_involution(self) -> None:
        rng = np.random.default_rng(0)
        m = csr_matrix(rng.random((5, 5)))
        twice = np.asarray(_reverse_groups_2d(_reverse_groups_2d(m)).todense())
        assert np.array_equal(twice, np.asarray(m.todense()))

    def test_reverses_both_axes_as_a_full_permutation(self) -> None:
        m = np.arange(9, dtype=float).reshape(3, 3)
        rev = np.asarray(_reverse_groups_2d(csr_matrix(m)).todense())
        assert np.array_equal(rev, m[::-1, ::-1])


class TestLoadedIsotopeIsCanonical:
    """The production XS data obeys the canonical fast-first convention."""

    def test_eg_is_strictly_descending(self) -> None:
        iso = load_isotope("U_235", 294)
        assert np.all(np.diff(iso.eg) < 0), "eg must be DESCENDING (group 0 = highest E)"
        assert iso.eg[0] > iso.eg[-1]

    def test_group_zero_is_the_fastest_group(self) -> None:
        # U-235 fission is huge at thermal, tiny at fast ⟹ group 0 (fast) ≪ group -1 (thermal).
        iso = load_isotope("U_235", 294)
        sigF = np.asarray(iso.sigF)[-1]  # most-dilute base point
        assert sigF[0] < sigF[-1], "group 0 must be the FAST group (small fission XS)"

    def test_downscatter_is_upper_triangular(self) -> None:
        iso = load_isotope("U_235", 294)
        sig_s0 = np.asarray(iso.sigS[0][-1].todense())
        assert np.triu(sig_s0, 1).sum() > np.tril(sig_s0, -1).sum(), (
            "downscatter must dominate the UPPER triangle (g_to > g_from) in fast-first order"
        )
