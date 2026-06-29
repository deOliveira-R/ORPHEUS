r"""Intrinsic-property gate for the energy-condensation value objects.

P5.0 (GitHub #274). Energy condensation collapses fine-group cross
sections onto a coarse group structure, spectrum-weighted. The value
objects that carry its structure:

* ``EnergyGrid`` — the energy analogue of a coarse ``Mesh``: a strictly
  DESCENDING boundary array (the canonical fast-first convention — see
  ``tests/data/test_gendf_canonical_order.py`` and
  ``docs/theory/cross_section_data.rst``). N+1 boundaries → N groups. It
  yields the binary fine→coarse mismatch treatment ``fine.overlap_to(
  coarse)``.
* ``OverlapBasis`` (returned by ``overlap_to``) — the fine→coarse trial
  carrying the fractional membership table ``overlap_table`` (a partition
  of unity, rows summing to 1; nested → one-hot). Its ``dominant_column``
  is the per-fine-group dominant coarse index (the former
  ``GroupCondensation.coarse_of_fine``), and ``fractional_columns`` reports
  the coarse groups that leaned on the within-group model (the former
  ``GroupCondensation.locally_interpolated``).

This file pins the DEFINING LAWS of those types (the user's standing
rule: every math-bearing type ships a test of its intrinsic properties).
Foundation tests — software invariants, no theory ``:label:``.

The orientation trap (P5.0 memo FLAG 1)
---------------------------------------
``IndicatorBasis.evaluate`` buckets points with ``searchsorted(side=
"right")``, which assumes ASCENDING edges. The canonical ``eg`` is
DESCENDING. Feeding descending coarse edges to ``IndicatorBasis`` does
NOT break the partition (row-sums stay 1, because searchsorted+clip
always returns an in-range bin) — it silently REVERSES the coarse-group
order (fast fine groups land in the thermal coarse column). So the
``overlap_to`` law that MUST be pinned is the coarse-group ORDER
(fast-first), not merely partition completeness. ``test_orientation_trap``
*demonstrates* the reversal directly against the live ``IndicatorBasis``
(it needs no SUT), so the implementer cannot ship a condensation whose
coarse groups are silently thermal-first.

vv Mode-8: every assertion is ``np.testing.*`` / ``pytest.raises`` /
``pytest.fail`` (the suite runs ``-O``; bare ``assert`` is stripped).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.energy_grid import EnergyGrid, InverseEnergySpectrum
from orpheus.numerics.basis.indicator_basis import IndicatorBasis

pytestmark = pytest.mark.foundation


# A small DESCENDING energy grid: 4 groups, fast-first (group 0 fastest).
_EG_DESC = np.array([1.0e7, 1.0e6, 1.0e3, 1.0e0, 1.0e-5])  # 5 edges → 4 groups


# ══════════════════════════════════════════════════════════════════════
# The orientation trap — demonstrated against the LIVE IndicatorBasis
# (needs NO SUT; pins the FLAG-1 silent coarse-group-order reversal)
# ══════════════════════════════════════════════════════════════════════


class TestOrientationTrap:
    """The descending-edge / ``searchsorted`` coarse-group-order reversal.

    These run unconditionally — they exercise the existing
    ``IndicatorBasis`` to prove WHY ``GroupCondensation`` must pin the
    coarse-group ORDER, not just partition completeness.
    """

    @staticmethod
    def _fine_points(eg_desc: np.ndarray) -> np.ndarray:
        """Geometric-midpoint representative energy of each fine group."""
        return np.sqrt(eg_desc[:-1] * eg_desc[1:])

    def test_descending_edges_keep_partition_but_reverse_order(self) -> None:
        """Descending coarse edges → valid partition BUT reversed coarse order.

        Fast fine groups (low fine-index, high energy) land in the HIGH
        coarse-column when edges are descending — the silent reversal the
        condensation MUST NOT inherit.
        """
        fine_pts = self._fine_points(_EG_DESC)            # 4 pts, descending energy
        # Coarse: fast {g0,g1}, thermal {g2,g3}. DESCENDING coarse edges.
        coarse_desc = np.array([1.0e7, 1.0e3, 1.0e-5])
        table = IndicatorBasis(edges_per_axis=(coarse_desc,)).evaluate(fine_pts)
        # Partition is "valid" (every fine group → exactly one coarse cell):
        np.testing.assert_array_equal(
            table.sum(axis=1), np.ones(4),
            err_msg="descending edges should still give a one-hot partition",
        )
        # ...but the column assignment is REVERSED: fast fine groups (0,1)
        # land in coarse column 1, thermal (2,3) in column 0.
        assignment = table.argmax(axis=1)
        np.testing.assert_array_equal(
            assignment, np.array([1, 1, 0, 0]),
            err_msg="descending edges silently reverse the coarse-group order",
        )

    def test_ascending_edges_required_for_fast_first(self) -> None:
        """Reversing edges to ASCENDING is what the searchsorted bucketing needs.

        With ascending coarse edges the fast fine groups map to the HIGH
        ascending-cell index; the implementer must reverse that cell axis
        back to recover fast-first coarse groups (or build the partition map
        explicitly, orientation-free).
        """
        fine_pts = self._fine_points(_EG_DESC)
        coarse_asc = np.array([1.0e-5, 1.0e3, 1.0e7])     # ascending
        table = IndicatorBasis(edges_per_axis=(coarse_asc,)).evaluate(fine_pts)
        np.testing.assert_array_equal(table.sum(axis=1), np.ones(4))
        # Fast fine groups (high energy) → high ascending-cell index (1).
        np.testing.assert_array_equal(
            table.argmax(axis=1), np.array([1, 1, 0, 0]),
            err_msg="ascending-edge cell index puts fast groups in the high cell",
        )


# ══════════════════════════════════════════════════════════════════════
# EnergyGrid — intrinsic properties (strict-descending, positive, partition)
# ══════════════════════════════════════════════════════════════════════


class TestEnergyGridIntrinsic:
    """The DEFINING laws of an EnergyGrid (the #265 invariant slice)."""

    def test_valid_descending_grid_constructs(self) -> None:
        """A strictly-descending, all-positive boundary array constructs."""
        grid = EnergyGrid(_EG_DESC)
        # N+1 boundaries → N groups (the defining count law).
        np.testing.assert_array_equal(
            grid.n_groups, len(_EG_DESC) - 1,
            err_msg="N+1 boundaries must give N groups",
        )

    def test_boundaries_strictly_descending(self) -> None:
        """The stored boundaries are strictly decreasing (fast-first)."""
        grid = EnergyGrid(_EG_DESC)
        edges = np.asarray(grid.edges)
        np.testing.assert_array_less(
            edges[1:], edges[:-1],
            err_msg="EnergyGrid boundaries must be strictly DESCENDING",
        )

    def test_ascending_grid_raises(self) -> None:
        """An ASCENDING boundary array violates the fast-first law → raises."""
        with pytest.raises((ValueError, AssertionError)):
            EnergyGrid(_EG_DESC[::-1])

    def test_nonmonotone_grid_raises(self) -> None:
        """A non-monotone boundary array (a repeat / inversion) → raises."""
        bad = np.array([1.0e7, 1.0e6, 1.0e6, 1.0e0, 1.0e-5])  # repeat → not strict
        with pytest.raises((ValueError, AssertionError)):
            EnergyGrid(bad)

    def test_nonpositive_energy_raises(self) -> None:
        """A non-positive boundary (energy ≤ 0) → raises (energies are positive).

        NOTE: WIMS-69 stores a bottom edge of exactly 0.0 eV. If EnergyGrid
        admits a 0.0 floor, the implementer should relax this to ``< 0``
        ONLY — keep the strict-positive-interior law and pin the floor
        convention explicitly. This test pins the strict-positive default;
        adjust to the chosen floor convention at implementation time.
        """
        bad = np.array([1.0e7, 1.0e6, 1.0e3, 0.0, -1.0])  # negative floor
        with pytest.raises((ValueError, AssertionError)):
            EnergyGrid(bad)

    def test_partition_completeness(self) -> None:
        """The groups tile the energy axis with no gaps/overlaps.

        Each group g occupies [edges[g+1], edges[g]) (descending), and the
        union of the groups is exactly [floor, ceiling] — a complete
        partition. Pinned via the group-width sum equalling the total span.
        """
        grid = EnergyGrid(_EG_DESC)
        edges = np.asarray(grid.edges)
        widths = edges[:-1] - edges[1:]                   # per-group width > 0
        np.testing.assert_array_less(
            0.0, widths, err_msg="every group must have positive width"
        )
        np.testing.assert_allclose(
            widths.sum(), edges[0] - edges[-1], rtol=1e-12,
            err_msg="group widths must sum to the total energy span (complete)",
        )


class TestEnergyGridGeometry:
    """The per-group geometry properties (#276 P4-F): representative energy,
    energy widths, lethargy widths — the diagnostics the homogeneous / MoC
    spectra read off the grid instead of re-deriving the group structure."""

    # A small hand-computable descending grid: 2 groups over [1, 100] eV, each
    # spanning exactly one decade.
    _GRID = EnergyGrid(np.array([100.0, 10.0, 1.0]))

    def test_energy_widths_are_upper_minus_lower(self) -> None:
        r"""\Delta E_g = E_up - E_lo (hand-checked, positive)."""
        np.testing.assert_allclose(self._GRID.energy_widths, [90.0, 9.0], rtol=1e-12)

    def test_energy_widths_sum_to_total_span(self) -> None:
        r"""Completeness in energy: \sum \Delta E_g = ceiling - floor."""
        g = self._GRID
        np.testing.assert_allclose(
            g.energy_widths.sum(), g.edges[0] - g.edges[-1], rtol=1e-12,
        )

    def test_lethargy_widths_are_log_ratio(self) -> None:
        r"""\Delta u_g = ln(E_up / E_lo) (hand-checked; each group spans a
        decade so \Delta u = ln 10)."""
        np.testing.assert_allclose(
            self._GRID.lethargy_widths, [np.log(10.0), np.log(10.0)], rtol=1e-12,
        )

    def test_lethargy_widths_sum_to_total_lethargy(self) -> None:
        r"""Completeness in lethargy: \sum \Delta u_g = ln(ceiling / floor)."""
        g = self._GRID
        np.testing.assert_allclose(
            g.lethargy_widths.sum(), np.log(g.edges[0] / g.edges[-1]), rtol=1e-12,
        )

    def test_representative_energy_is_geometric_not_arithmetic(self) -> None:
        r"""The group centre is the GEOMETRIC mean :math:`\sqrt{E_{up} E_{lo}}` —
        the natural centre on the log energy axis — STRICTLY below the
        arithmetic midpoint (AM-GM). This is the P4-F correctness property: the
        homogeneous / MoC spectrum abscissa uses THIS, not the old
        ``0.5*(E_up+E_lo)``; a regression to the arithmetic midpoint reddens here.
        """
        g = self._GRID
        upper, lower = g.edges[:-1], g.edges[1:]
        np.testing.assert_allclose(
            g.representative_energy, np.sqrt(upper * lower), rtol=1e-12,
        )
        np.testing.assert_allclose(
            g.representative_energy, [np.sqrt(1000.0), np.sqrt(10.0)], rtol=1e-12,
        )
        np.testing.assert_array_less(g.representative_energy, 0.5 * (upper + lower))


# ══════════════════════════════════════════════════════════════════════
# overlap_to — intrinsic properties (containing-interval + fast-first)
# ══════════════════════════════════════════════════════════════════════


class TestOverlapToIntrinsic:
    """The DEFINING containing-interval law + the fast-first order pin.

    ``EnergyGrid.overlap_to(coarse)`` returns the ``OverlapBasis`` carrying
    the fine→coarse mismatch table; its ``dominant_column`` is the per-fine
    dominant coarse index (the former ``GroupCondensation.coarse_of_fine``).
    For a NESTED map every fine group maps wholly to one coarse group, so
    ``dominant_column`` IS the exact containing-interval partition.
    """

    @staticmethod
    def _fine_grid() -> np.ndarray:
        """8 fine groups, strictly descending (9 boundaries)."""
        return np.array([1e7, 5e6, 1e6, 1e3, 1e1, 1e0, 1e-1, 1e-2, 1e-3])

    @staticmethod
    def _coarse_grid() -> np.ndarray:
        """3 coarse groups: {g0,g1}, {g2,g3,g4}, {g5,g6,g7}. Descending."""
        return np.array([1e7, 1e6, 1e0, 1e-3])  # 4 edges → 3 coarse groups

    def _build(self):
        """Build the OverlapBasis from fine→coarse EnergyGrids."""
        fine = EnergyGrid(self._fine_grid())
        coarse = EnergyGrid(self._coarse_grid())
        return fine.overlap_to(coarse)

    @staticmethod
    def _coarse_of_fine(overlap) -> np.ndarray:
        """The per-fine-group dominant coarse index as a 1-D int array.

        ``OverlapBasis.dominant_column`` is the containing-interval map for a
        nested grid (and the argmax-fraction for a straddle).
        """
        return np.asarray(overlap.dominant_column, dtype=int).ravel()

    def test_every_fine_maps_to_exactly_one_coarse(self) -> None:
        """Containing-interval: each fine group → exactly one coarse group."""
        overlap = self._build()
        cof = self._coarse_of_fine(overlap)
        np.testing.assert_array_equal(
            cof.shape, (8,), err_msg="one coarse index per fine group"
        )
        # Every coarse index is a valid group id in [0, 3).
        np.testing.assert_array_less(
            cof.max(), 3, err_msg="coarse indices must be < n_coarse"
        )
        np.testing.assert_array_less(
            -1, cof.min(), err_msg="coarse indices must be ≥ 0"
        )

    def test_partition_is_contiguous_no_gaps_overlaps(self) -> None:
        """Each coarse group's fine-group set is a contiguous run (no gaps)."""
        overlap = self._build()
        cof = self._coarse_of_fine(overlap)
        # The expected fast-first contiguous partition of 8 fine into 3 coarse.
        np.testing.assert_array_equal(
            cof, np.array([0, 0, 1, 1, 1, 2, 2, 2]),
            err_msg="fine→coarse must be the contiguous fast-first partition",
        )
        # Contiguity: the coarse index is non-decreasing along fine index
        # (fast-first → coarse 0 first), and each value appears in one run.
        np.testing.assert_array_less(
            np.diff(cof), 2, err_msg="coarse index must step by ≤1 (contiguous)"
        )
        np.testing.assert_array_equal(
            np.all(np.diff(cof) >= 0), True,
            err_msg="coarse index must be non-decreasing (fast-first order)",
        )

    def test_coarse_group_order_is_fast_first(self) -> None:
        """The coarse-group ORDER is fast-first (FLAG-1 orientation pin).

        Coarse group 0 must contain the FASTEST fine groups (lowest fine
        index / highest energy). A condensation built on raw descending
        edges would silently reverse this — coarse 0 would be thermal.
        """
        overlap = self._build()
        cof = self._coarse_of_fine(overlap)
        # Fine group 0 (fastest) MUST land in coarse group 0.
        np.testing.assert_array_equal(
            cof[0], 0,
            err_msg="fastest fine group must map to coarse group 0 (fast-first)",
        )
        # Fine group 7 (slowest) MUST land in the LAST coarse group.
        np.testing.assert_array_equal(
            cof[-1], 2,
            err_msg="slowest fine group must map to the last coarse group",
        )

    def test_mismatched_grids_raise(self) -> None:
        """A coarse grid not nested in the fine grid (non-coincident span) raises.

        The containing-interval rule requires the coarse boundaries to lie
        within the fine span. A coarse ceiling above the fine ceiling is
        ill-posed (``overlap_to`` raises the span guard).
        """
        fine = EnergyGrid(self._fine_grid())
        bad_coarse = EnergyGrid(np.array([2e7, 1e6, 1e-3]))  # ceiling 2e7 > fine 1e7
        with pytest.raises((ValueError, AssertionError)):
            fine.overlap_to(bad_coarse)


# ══════════════════════════════════════════════════════════════════════
# FRACTIONAL-OVERLAP re-binning (#274 follow-up) — the NON-NESTED case
# ══════════════════════════════════════════════════════════════════════
#
# The production case (421-group → WIMS-69/172) is NON-NESTED: a coarse
# boundary can fall INSIDE a fine group, so that fine group must apportion
# a FRACTION of its rate to each coarse group it overlaps (conservative
# fractional re-binning; Hébert §3.5, the rank-0 case of Generalized Energy
# Condensation). ``GroupCondensation.table`` becomes FRACTIONAL T[g,G]∈[0,1]
# with rows summing to 1 (partition of unity); nested → one-hot (the
# regression-safe degenerate). The fraction is the flux-weighted overlap
#
#     f_{g,G} = (∫_{g∩G} w(E) dE) / (∫_g w(E) dE)
#
# with w(E) the SELECTABLE within-group flux model. For 1/E (flat-in-
# lethargy, the default ``InverseEnergySpectrum``) this is the LETHARGY
# ratio ln(hi_ov/lo_ov)/ln(hi_g/lo_g).
#
# THE STRADDLE FIXTURE (the discriminating geometry — a non-straddle makes
# fractional ≡ one-hot, so every fractional gate uses THIS):
#   Fine edges (descending): [1e6, 1e4, 1e2, 1e0, 1e-2] → 4 fine groups
#     g0 [1e4,1e6)  g1 [1e2,1e4)  g2 [1e0,1e2)  g3 [1e-2,1e0)
#   Coarse boundary at 1e1 eV falls INSIDE g2 [1e0,1e2). Coarse edges
#     [1e6, 1e1, 1e-2] → G0 [1e1,1e6) (fast)  G1 [1e-2,1e1) (thermal).
#   ⇒ g0,g1 → G0 fully; g3 → G1 fully; g2 STRADDLES (1/E split [0.5,0.5]).
#
# The oracle table is computed INDEPENDENTLY of the SUT (hand-coded energy
# overlap + the within-group model), NEVER read back from the SUT — the vv
# L11 structural-independence requirement.

_FINE_STRADDLE = np.array([1.0e6, 1.0e4, 1.0e2, 1.0e0, 1.0e-2])  # 4 fine groups
_COARSE_STRADDLE = np.array([1.0e6, 1.0e1, 1.0e-2])              # 2 coarse, cut @1e1
_NG_FINE_S, _NG_COARSE_S = 4, 2


def _inv_e_weight(lo: float, hi: float) -> float:
    """1/E within-group integrated weight ∫_lo^hi dE/E = ln(hi/lo) (lethargy)."""
    return float(np.log(hi / lo))


def _flat_e_weight(lo: float, hi: float) -> float:
    """A DIFFERENT within-group model: flat-in-ENERGY ∫_lo^hi dE = hi-lo.

    Constructed in-test (the F4 negative flavour) to prove the 1/E model is
    genuinely SELECTABLE / load-bearing: it gives a different straddle split.
    """
    return float(hi - lo)


def _fractional_table_oracle(fine_e: np.ndarray, coarse_e: np.ndarray, weight) -> np.ndarray:
    """Hand-build the FRACTIONAL membership T[g,G] from energy overlaps.

    Structurally independent of the SUT: pure energy arithmetic +the chosen
    within-group ``weight`` model. T[g,G] = (∫_{g∩G} w)/(∫_g w); rows sum to 1
    for a coarse grid spanning the fine grid (partition of unity).
    """
    ng_f = fine_e.size - 1
    ng_c = coarse_e.size - 1
    T = np.zeros((ng_f, ng_c))
    for g in range(ng_f):
        g_hi, g_lo = float(fine_e[g]), float(fine_e[g + 1])     # descending
        denom = weight(g_lo, g_hi)
        for G in range(ng_c):
            G_hi, G_lo = float(coarse_e[G]), float(coarse_e[G + 1])
            ov_lo, ov_hi = max(g_lo, G_lo), min(g_hi, G_hi)
            if ov_hi > ov_lo:
                T[g, G] = weight(ov_lo, ov_hi) / denom
    return T


def _build_fractional(fine_e: np.ndarray, coarse_e: np.ndarray, within_group=None):
    """Build the fractional OverlapBasis onto the given grids.

    ``fine.overlap_to(coarse, within_group)`` IS the fractional mismatch
    treatment — the within-group model defaults to 1/E when ``None``.
    """
    fine = EnergyGrid(fine_e)
    coarse = EnergyGrid(coarse_e)
    if within_group is None:
        return fine.overlap_to(coarse)
    return fine.overlap_to(coarse, within_group)


# ══════════════════════════════════════════════════════════════════════
# F4 (oracle half) — the 1/E within-group model IS the lethargy ratio
# (pure foundation: runs without ANY SUT — it pins the math the SUT must
#  reproduce, and proves 1/E ≠ flat-energy for the straddle)
# ══════════════════════════════════════════════════════════════════════


class TestF4WithinGroupModelOracle:
    """The 1/E (InverseEnergySpectrum) overlap fraction = lethargy ratio.

    Two unconditional legs: (a) the 1/E straddle split is exactly the
    hand-computed lethargy ratio; (b) a DIFFERENT within-group model
    (flat-in-energy) gives a DIFFERENT split — so the model is genuinely
    load-bearing (vv #11 positive + negative). These run with no SUT; the
    SUT-side discriminator (that ``InverseEnergySpectrum`` reproduces leg-a)
    lives in tests/data/test_mixture_condense.py.
    """

    def test_inv_e_straddle_fraction_is_lethargy_ratio(self) -> None:
        """g2 [1e0,1e2) split by the coarse cut at 1e1 → ln(10)/ln(100) = 0.5."""
        T = _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, _inv_e_weight)
        # g2's overlap with G0 is [1e1,1e2): lethargy ln(1e2/1e1)=ln(10);
        # full g2 lethargy ln(1e2/1e0)=ln(100). Ratio = 0.5.
        np.testing.assert_allclose(
            T[2], np.array([0.5, 0.5]), atol=1e-12,
            err_msg="1/E straddle split must be the lethargy ratio ln(10)/ln(100)",
        )

    def test_flat_energy_model_gives_different_split(self) -> None:
        """flat-in-ENERGY gives g2→G0 = 90/99 ≈ 0.909, NOT 0.5 — model matters.

        Negative flavour: a wrong within-group model changes the straddle
        apportionment, so the model is not droppable / not interchangeable.
        """
        T_inv = _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, _inv_e_weight)
        T_flat = _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, _flat_e_weight)
        # flat-energy: (1e2-1e1)/(1e2-1e0) = 90/99.
        np.testing.assert_allclose(
            T_flat[2, 0], 90.0 / 99.0, atol=1e-12,
            err_msg="flat-energy straddle fraction must be (1e2-1e1)/(1e2-1e0)",
        )
        # The two models DISAGREE on the straddle (the load-bearing point):
        if abs(float(T_flat[2, 0] - T_inv[2, 0])) < 1e-3:
            pytest.fail(
                "within-group model is NOT load-bearing: 1/E and flat-energy give "
                "the same straddle split — the F4 SUT discriminator would be blind"
            )

    def test_both_models_are_partitions_of_unity(self) -> None:
        """Either within-group model yields rows summing to 1 (conservation).

        The fraction normalisation Σ_G f_{g,G} = (Σ_G ∫_{g∩G} w)/(∫_g w) = 1
        holds for ANY positive w because the coarse grid tiles the fine span.
        """
        for weight in (_inv_e_weight, _flat_e_weight):
            T = _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, weight)
            np.testing.assert_allclose(
                T.sum(axis=1), np.ones(_NG_FINE_S), atol=1e-12,
                err_msg="any within-group model must give a partition of unity",
            )


# ══════════════════════════════════════════════════════════════════════
# F2 — Partition of unity (fractional rows sum to 1; nested → one-hot {0,1})
# ══════════════════════════════════════════════════════════════════════


class TestF2PartitionOfUnity:
    """table.sum(axis=1) == 1 (fractional straddle) AND {0,1} bit-exact (nested)."""

    def test_fractional_rows_sum_to_one(self) -> None:
        """The straddle table's rows sum to 1 (a genuinely fractional row exists)."""
        overlap = _build_fractional(
            _FINE_STRADDLE, _COARSE_STRADDLE,
            within_group=InverseEnergySpectrum(),
        )
        T = np.asarray(overlap.overlap_table, dtype=float)
        np.testing.assert_array_equal(
            T.shape, (_NG_FINE_S, _NG_COARSE_S),
            err_msg="fractional table must be (n_fine, n_coarse)",
        )
        np.testing.assert_allclose(
            T.sum(axis=1), np.ones(_NG_FINE_S), atol=1e-12,
            err_msg="fractional table rows must sum to 1 (partition of unity)",
        )
        # The straddle row (g2) MUST be genuinely fractional (not one-hot) —
        # otherwise the fixture is not exercising the fractional path.
        straddle = T[2]
        is_one_hot = bool(np.any(np.isclose(straddle, 1.0)) and np.sum(straddle > 1e-12) == 1)
        if is_one_hot:
            pytest.fail(
                "g2 row is one-hot — the straddle fixture is not fractional; "
                "the partition-of-unity test is not exercising the fractional path"
            )
        np.testing.assert_array_less(
            1e-12, straddle,
            err_msg="g2 straddles → BOTH coarse columns must receive a fraction",
        )

    def test_nested_table_is_bit_exact_one_hot(self) -> None:
        """A NESTED partition gives a {0,1} table, bit-identical to the one-hot.

        ``array_equal`` against the explicit one-hot (built from the
        containing-interval ``dominant_column``): a fractional impl that ALWAYS
        splits (never snaps to {0,1} on alignment) would drift → RED. The
        nested grid puts every coarse boundary ON a fine boundary so no group
        straddles.
        """
        # Nested: coarse edges are a SUBSET of the fine edges → no straddle.
        fine_e = _FINE_STRADDLE                                   # 4 fine groups
        coarse_e = np.array([1.0e6, 1.0e2, 1.0e-2])              # cut @1e2 (a fine edge)
        overlap = _build_fractional(fine_e, coarse_e, within_group=InverseEnergySpectrum())
        T = np.asarray(overlap.overlap_table, dtype=float)
        # Independent one-hot oracle from the containing-interval map.
        cof = np.asarray(overlap.dominant_column, dtype=int).ravel()
        one_hot = np.zeros_like(T)
        one_hot[np.arange(T.shape[0]), cof] = 1.0
        np.testing.assert_array_equal(
            T, one_hot,
            err_msg="nested fractional table must be the bit-exact {0,1} one-hot",
        )
        # And it really is in {0,1} (no fractional entries crept in):
        np.testing.assert_array_equal(
            np.isin(T, (0.0, 1.0)), np.ones_like(T, dtype=bool),
            err_msg="nested table entries must be exactly 0 or 1",
        )


# ══════════════════════════════════════════════════════════════════════
# F3 — Nested degeneracy / regression (the fractional table == the old
#      one-hot containing-interval result on the easy case)
# ══════════════════════════════════════════════════════════════════════


class TestF3NestedDegeneracy:
    """A nested fine→coarse fractional table reproduces the one-hot exactly.

    The 'no regression on the easy case' pin. The EXISTING nested intrinsic
    gates (TestOverlapToIntrinsic) still pass unchanged — this adds the
    specific claim that the fractional ``overlap_table`` equals the one-hot
    ``dominant_column`` map for a coarse grid whose boundaries coincide with
    fine boundaries.
    """

    def test_fractional_reduces_to_one_hot_when_nested(self) -> None:
        """8 fine → 3 coarse, all coarse edges ON fine edges → one-hot table."""
        fine_e = np.array([1e7, 5e6, 1e6, 1e3, 1e1, 1e0, 1e-1, 1e-2, 1e-3])  # 8 groups
        coarse_e = np.array([1e7, 1e6, 1e0, 1e-3])  # cuts at fine edges 1e6,1e0
        overlap = _build_fractional(fine_e, coarse_e, within_group=InverseEnergySpectrum())
        T = np.asarray(overlap.overlap_table, dtype=float)
        # Expected fast-first contiguous one-hot: {g0,g1}→0, {g2,g3,g4}→1, {g5,g6,g7}→2
        expected_cof = np.array([0, 0, 1, 1, 1, 2, 2, 2])
        one_hot = np.zeros((8, 3))
        one_hot[np.arange(8), expected_cof] = 1.0
        np.testing.assert_array_equal(
            T, one_hot,
            err_msg="nested fractional table must equal the contiguous fast-first one-hot",
        )

    def test_coarse_of_fine_unchanged_for_nested(self) -> None:
        """The dominant-coarse report (argmax of table) is the contiguous map."""
        fine_e = np.array([1e7, 5e6, 1e6, 1e3, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
        coarse_e = np.array([1e7, 1e6, 1e0, 1e-3])
        overlap = _build_fractional(fine_e, coarse_e, within_group=InverseEnergySpectrum())
        np.testing.assert_array_equal(
            np.asarray(overlap.dominant_column, dtype=int).ravel(),
            np.array([0, 0, 1, 1, 1, 2, 2, 2]),
            err_msg="nested dominant_column must be the contiguous fast-first partition",
        )


# ══════════════════════════════════════════════════════════════════════
# F5 — Upscaling guard (condensation is downsampling-ONLY)
# ══════════════════════════════════════════════════════════════════════


class TestF5UpscalingGuard:
    """coarse.n_groups > fine.n_groups MUST raise — condensation cannot
    fabricate finer structure. Positive control: a valid coarser target does
    NOT raise. The guard lives in ``EnergyGrid.overlap_to``.
    """

    def test_upscaling_target_raises(self) -> None:
        """A coarse grid FINER than the fine grid (more groups) → ValueError.

        64 fine groups → a 200-group 'coarse' target is upscaling, which
        condensation forbids (the few-group structure must be coarser). The
        guard belongs at the ``overlap_to`` (orientation-/span-check) site.
        """
        # 64 fine groups (65 descending edges, log-spaced 1e7..1e-3).
        fine = EnergyGrid(np.logspace(7, -3, 65))
        # 200 'coarse' groups within the SAME span — strictly finer ⇒ upscaling.
        coarse = EnergyGrid(np.logspace(7, -3, 201))
        with pytest.raises((ValueError, AssertionError)):
            fine.overlap_to(coarse)

    def test_valid_coarser_target_does_not_raise(self) -> None:
        """Positive control: a genuinely coarser target constructs cleanly.

        Pins that F5's guard is the upscaling DIRECTION, not a blanket reject
        of unequal group counts (vv #11 positive leg — the guard must not be
        a no-op-that-always-raises nor an always-passes).
        """
        fine = EnergyGrid(np.logspace(7, -3, 65))     # 64 fine
        coarse = EnergyGrid(np.logspace(7, -3, 9))    # 8 coarse, nested
        # MUST NOT raise — the downsampling direction is allowed.
        fine.overlap_to(coarse)
        # n_coarse < n_fine — the downsampling direction is the valid one.
        np.testing.assert_array_less(
            coarse.n_groups, fine.n_groups,
            err_msg="a coarser target must construct (downsampling is allowed)",
        )


# ══════════════════════════════════════════════════════════════════════
# F6 — Local-interpolation report (which coarse groups leaned on w(E))
# ══════════════════════════════════════════════════════════════════════


class TestF6LocalInterpolationReport:
    """``fractional_columns`` is EMPTY for a nested map, NON-EMPTY for a
    straddle (the coarse groups whose data leans on the within-group model).
    """

    def test_nested_reports_no_interpolation(self) -> None:
        """A nested condensation interpolated nothing → empty report."""
        fine_e = np.array([1e7, 5e6, 1e6, 1e3, 1e1, 1e0, 1e-1, 1e-2, 1e-3])
        coarse_e = np.array([1e7, 1e6, 1e0, 1e-3])  # all cuts on fine edges
        overlap = _build_fractional(fine_e, coarse_e, within_group=InverseEnergySpectrum())
        li = np.asarray(overlap.fractional_columns, dtype=int).ravel()
        np.testing.assert_array_equal(
            li.size, 0,
            err_msg="nested condensation must report NO locally-interpolated groups",
        )

    def test_straddle_reports_interpolated_coarse_groups(self) -> None:
        """The straddle leans on w(E) → the affected coarse groups are reported.

        g2 straddles the 1e1 cut, fanning fractionally into BOTH G0 and G1, so
        BOTH coarse groups received a fractional (straddle) contribution.
        """
        overlap = _build_fractional(
            _FINE_STRADDLE, _COARSE_STRADDLE,
            within_group=InverseEnergySpectrum(),
        )
        li = set(np.asarray(overlap.fractional_columns, dtype=int).ravel().tolist())
        np.testing.assert_array_less(
            0, len(li),
            err_msg="a straddle must report ≥1 locally-interpolated coarse group",
        )
        # Both coarse groups that g2 fans into are interpolation-dependent.
        T = np.asarray(overlap.overlap_table, dtype=float)
        fractional_cols = set(np.nonzero((T > 1e-12) & (T < 1.0 - 1e-12))[1].tolist())
        np.testing.assert_array_equal(
            fractional_cols.issubset(li), True,
            err_msg="every coarse column that received a FRACTION must be reported "
            f"as locally-interpolated; fractional cols={fractional_cols}, report={li}",
        )
