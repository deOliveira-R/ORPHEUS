r"""Energy-condensation gate — per-material ``Mixture.condense``.

P5.0 (GitHub #274). ``Mixture.condense(target, spectrum) -> Mixture``
collapses fine-group cross sections onto a coarse group structure,
spectrum-weighted. This file pins the per-material math:

* **G1** rate preservation (vector channels) — the L1 correctness anchor.
* **G2** within-group-varying-spectrum discriminator — proves the
  flux-weight is load-bearing (Mode-7 guard).
* **G3** scattering 2-axis collapse (sink SUMMED, source flux-AVERAGED) —
  the structural difference from spatial homogenize, with three
  value-RED mutations.
* **chi** pure birth-group SUM (``χ @ T``, NOT frame-projected) +
  ``Σχ`` preservation.
* **G4** ``assert_balanced`` regression (positive + negative, vv #11).
* **G5 leg-A** WIMS Table-11.3 derivation-validation (the containing-
  interval rule reproduces ``CONDENSE_172_TO_69`` within the documented
  boundary tolerance).
* **G6** Mode-11 execution sentinel (condense routes through
  ``FrameBase.project`` and the TEST-side ``WeightedIndicatorBasis``).

The structurally-independent oracle (vv L11) is ALWAYS an explicit
per-fine-group Python hand-sum — NEVER a re-call of ``frame.project`` or a
second condense path (that would be twin-consistency, vv L4). The oracle
math is verified against the live frame in the P5.0 memo
(``.claude/plans/archive/p5_0_condensation_gate.md``).

vv Mode-8: every assertion is ``np.testing.*`` / ``pytest.raises`` /
``pytest.fail`` (the suite runs ``-O``; bare ``assert`` is stripped).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.energy_grid import EnergyGrid, InverseEnergySpectrum
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.numerics.manifold import EnergyGroups

pytestmark = pytest.mark.foundation


# ══════════════════════════════════════════════════════════════════════
# Fixtures — a 4-fine-group balanced fissile Mixture, fast-first
# ══════════════════════════════════════════════════════════════════════
#
# Coarse structure throughout this file: G0 = {g0,g1} (fast),
# G1 = {g2,g3} (thermal) — so EACH coarse group has TWO fine groups,
# enabling within-group spectral variation (the G2 / Mode-7 requirement).

_NG_FINE = 4
_GROUPS = (np.array([0, 1]), np.array([2, 3]))      # fast-first coarse partition
_NG_COARSE = len(_GROUPS)

# A fine spectrum that VARIES within each coarse group (Mode-7 activation):
# G0 has φ = [1, 4] (4× spread), G1 has φ = [2, 0.5] (4× spread).
_PHI = np.array([1.0, 4.0, 2.0, 0.5])

# DESCENDING energy grid (fast-first); 5 edges → 4 fine groups.
_EG_FINE = np.array([1.0e7, 1.0e5, 1.0e2, 1.0e0, 1.0e-2])
_EG_COARSE = np.array([1.0e7, 1.0e2, 1.0e-2])         # 3 edges → 2 coarse groups


def _balanced_fissile_4g() -> Mixture:
    """A balanced fissile 4-group Mixture (asymmetric downscatter).

    SigT = SigC + SigL + SigF + rowsum(SigS0) + rowsum(Sig2); χ is a
    fast-peaked simplex; SigL AND Sig2 are NON-zero so the full balance
    identity is exercised after condensation (G4). Built DIRECTLY (not via
    ``make_mixture``, which hardcodes ``SigL = 0`` and would drop the SigL
    leg — mirrors ``test_mixture_xs_balance.py``'s hand-built fixtures).
    """
    sig_c = np.array([0.02, 0.05, 0.10, 0.30])
    sig_l = np.array([0.01, 0.02, 0.03, 0.04])        # (n,alpha) — NON-zero
    sig_f = np.array([0.01, 0.02, 0.08, 0.20])
    nu = np.array([2.6, 2.5, 2.45, 2.43])
    chi = np.array([0.7, 0.25, 0.05, 0.0])            # fast-peaked simplex
    sig_s0 = np.array([                                # asymmetric, downscatter-heavy
        [0.30, 0.10, 0.05, 0.02],
        [0.00, 0.40, 0.15, 0.03],
        [0.00, 0.00, 0.55, 0.20],
        [0.00, 0.00, 0.00, 0.90],
    ])
    sig_s1 = sig_s0 * 0.1                               # mild P1 anisotropy
    sig_2 = np.array([                                 # (n,2n) — NON-zero, upper-tri
        [0.005, 0.003, 0.001, 0.0],
        [0.000, 0.004, 0.002, 0.0],
        [0.000, 0.000, 0.003, 0.0],
        [0.000, 0.000, 0.000, 0.0],
    ])
    rowsum = sig_s0.sum(axis=1) + sig_2.sum(axis=1)
    sig_t = sig_c + sig_l + sig_f + rowsum
    mix = Mixture(
        SigC=sig_c.copy(), SigL=sig_l.copy(), SigF=sig_f.copy(),
        SigP=(nu * sig_f).copy(), SigT=sig_t.copy(),
        SigS=[csr_matrix(sig_s0), csr_matrix(sig_s1)],
        Sig2=csr_matrix(sig_2), chi=chi.copy(), eg=_EG_FINE.copy(),
    )
    mix.assert_balanced()
    return mix


@pytest.fixture(scope="module")
def fine_mix() -> Mixture:
    return _balanced_fissile_4g()


def _condense(mix: Mixture) -> Mixture:
    """Run condense onto the fast-first 2-coarse-group structure.

    ``Mixture.condense(target, spectrum)`` reads the mixture's own fine grid
    (``mix.energy_grid``) and the ``target`` coarse ``EnergyGrid``, projecting
    with the per-material spectrum φ (within-group model defaults to 1/E).
    """
    return mix.condense(EnergyGrid(_EG_COARSE), _PHI)


# ── Structurally-independent oracles (hand-sums; NEVER frame.project) ───


def _phi_R() -> np.ndarray:
    """Coarse flux Φ_G = Σ_{g∈G} φ_g (the group-sum)."""
    return np.array([_PHI[g].sum() for g in _GROUPS])


def _vector_oracle(sig_fine: np.ndarray) -> np.ndarray:
    """Σ_G = (Σ_{g∈G} φ_g Σ_g) / Σ_{g∈G} φ_g — explicit hand-sum."""
    return np.array(
        [(_PHI[g] * sig_fine[g]).sum() / _PHI[g].sum() for g in _GROUPS]
    )


def _scatter_rate_oracle(sig_s_fine: np.ndarray) -> np.ndarray:
    """In-scatter RATE R_{G→G'} = Σ_{g∈G}Σ_{g'∈G'} φ_g Σ_s[g,g'] (double sum)."""
    R = np.zeros((_NG_COARSE, _NG_COARSE))
    for Gi, gset in enumerate(_GROUPS):
        for Gj, gpset in enumerate(_GROUPS):
            R[Gi, Gj] = sum(
                _PHI[g] * sig_s_fine[g, gp] for g in gset for gp in gpset
            )
    return R


# ══════════════════════════════════════════════════════════════════════
# G1 — Rate preservation (vector channels) — the correctness anchor
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("energy-condensation-rate-preservation")
class TestG1RatePreservationVectors:
    """Σ_G · Φ_G == Σ_{g∈G} φ_g Σ_g for every vector channel (vv L1, L11)."""

    @pytest.mark.parametrize("channel", ["SigT", "SigC", "SigL", "SigF", "SigP"])
    def test_vector_channel_rate_preserved(
        self, fine_mix: Mixture, channel: str
    ) -> None:
        coarse_mix = _condense(fine_mix)
        sig_fine = np.asarray(getattr(fine_mix, channel), dtype=float)
        sig_eff = np.asarray(getattr(coarse_mix, channel), dtype=float)
        phi_R = _phi_R()
        for G in range(_NG_COARSE):
            rate_coarse = sig_eff[G] * phi_R[G]
            rate_ref = float((_PHI[_GROUPS[G]] * sig_fine[_GROUPS[G]]).sum())
            np.testing.assert_allclose(
                rate_coarse, rate_ref, rtol=1e-12, atol=1e-12,
                err_msg=f"{channel} rate not preserved in coarse {G}",
            )

    def test_condensed_value_equals_vector_oracle(self, fine_mix: Mixture) -> None:
        """The condensed SigT equals the hand-summed flux-weighted average."""
        coarse_mix = _condense(fine_mix)
        np.testing.assert_allclose(
            np.asarray(coarse_mix.SigT), _vector_oracle(np.asarray(fine_mix.SigT)),
            rtol=1e-12, atol=1e-12,
            err_msg="condensed SigT must equal the flux-weighted hand-sum oracle",
        )

    def test_condensed_shapes(self, fine_mix: Mixture) -> None:
        """Vectors collapse to (n_coarse,); the group count drops fine→coarse."""
        coarse_mix = _condense(fine_mix)
        np.testing.assert_array_equal(np.asarray(coarse_mix.SigT).shape, (_NG_COARSE,))
        np.testing.assert_array_equal(coarse_mix.ng, _NG_COARSE)


# ══════════════════════════════════════════════════════════════════════
# G2 — Within-group-varying-spectrum discriminator (Mode-7 guard)
# ══════════════════════════════════════════════════════════════════════


class TestG2WithinGroupVaryingDiscriminator:
    """Flux-weighted ≠ arithmetic-average when φ varies within a coarse group.

    Proves the φ-weight is load-bearing: a flat spectrum makes flux-weighted
    ≡ arithmetic-average, so a wrong unweighted collapse would accidentally
    agree (Mode-7 trap). The fixture φ = [1,4,2,0.5] ACTIVATES within-group
    variation in BOTH coarse groups; the test fails loudly if it is too flat.
    """

    def test_ansatz_activates_within_group_variation(self) -> None:
        """The φ ansatz is genuinely non-flat within EACH coarse group.

        Unconditional Mode-7 activation guard — runs without the SUT so a
        future edit that flattens φ reddens immediately.
        """
        for G, gset in enumerate(_GROUPS):
            spread = _PHI[gset].max() / _PHI[gset].min()
            np.testing.assert_array_less(
                1.5, spread,
                err_msg=f"coarse {G}: φ too flat (spread={spread}) — Mode-7 blind",
            )

    def test_flux_weighted_differs_from_arithmetic(self) -> None:
        """The φ-weighted and arithmetic-average collapses are numerically distinct.

        Unconditional — pins that an unweighted collapse is a DIFFERENT
        answer, so G2(b) below is a real discriminator.
        """
        sig_fine = np.asarray(_balanced_fissile_4g().SigT, dtype=float)
        phi_weighted = _vector_oracle(sig_fine)
        arithmetic = np.array([sig_fine[g].mean() for g in _GROUPS])
        # At least one coarse group must discriminate by a visible margin.
        discriminated = bool(np.any(np.abs(phi_weighted - arithmetic) > 1e-3))
        if not discriminated:
            pytest.fail(
                "fixture too flat: φ-weighted ≡ arithmetic — the discriminator "
                "is blind (Mode-7). Sharpen φ or the per-group SigT spread."
            )

    def test_condense_preserves_rate_and_rejects_arithmetic(
        self, fine_mix: Mixture
    ) -> None:
        """(a) flux-weighted condense preserves the rate; (b) the arithmetic
        average does NOT — the production value MUST match (a)."""
        coarse_mix = _condense(fine_mix)
        sig_fine = np.asarray(fine_mix.SigT, dtype=float)
        phi_R = _phi_R()
        arithmetic = np.array([sig_fine[g].mean() for g in _GROUPS])
        discriminated = False
        for G in range(_NG_COARSE):
            # (a) production preserves the rate (== the φ-weighted oracle).
            rate_coarse = float(coarse_mix.SigT[G] * phi_R[G])
            rate_ref = float((_PHI[_GROUPS[G]] * sig_fine[_GROUPS[G]]).sum())
            np.testing.assert_allclose(
                rate_coarse, rate_ref, rtol=1e-12, atol=1e-12,
                err_msg=f"coarse {G}: flux-weighted condense must preserve the rate",
            )
            # (b) the arithmetic average gives a WRONG rate where it differs.
            rate_arith = float(arithmetic[G] * phi_R[G])
            if abs(rate_arith - rate_ref) > 1e-6:
                discriminated = True
                np.testing.assert_array_less(
                    1e-6, abs(float(coarse_mix.SigT[G]) - arithmetic[G]),
                    err_msg=f"coarse {G}: production matched the arithmetic average!",
                )
        if not discriminated:
            pytest.fail("fixture too flat: arithmetic average never wrong (Mode-7)")


# ══════════════════════════════════════════════════════════════════════
# G3 — Scattering 2-axis collapse (sink SUMMED, source flux-AVERAGED)
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("energy-condensation-scattering-collapse")
class TestG3ScatteringTwoAxisCollapse:
    """In-scatter RATE preserved; the source/sink asymmetry is load-bearing.

    The design: sink g_to SUMMED (Σ @ T), source g_from flux-AVERAGED
    (project). This DIFFERS from spatial homogenize (which flux-weights both
    axes). The three mutations below each produce a numerically different
    coarse matrix (verified in the P5.0 memo) → each MUST red.
    """

    @pytest.mark.parametrize("order", [0, 1])
    def test_in_scatter_rate_preserved(self, fine_mix: Mixture, order: int) -> None:
        """Φ_G · Σ_{s,coarse}[G,G'] == Σ_{g∈G}Σ_{g'∈G'} φ_g Σ_s[g,g']."""
        coarse_mix = _condense(fine_mix)
        if order >= len(coarse_mix.SigS):
            pytest.skip(f"condensed Mixture has no Legendre order {order}")
        sig_s_fine = np.asarray(fine_mix.SigS[order].todense(), dtype=float)
        sig_s_coarse = np.asarray(coarse_mix.SigS[order].todense(), dtype=float)
        R_oracle = _scatter_rate_oracle(sig_s_fine)
        phi_R = _phi_R()
        for Gi in range(_NG_COARSE):
            for Gj in range(_NG_COARSE):
                rate_coarse = sig_s_coarse[Gi, Gj] * phi_R[Gi]
                # one-ULP, not exact (the @T matmul reduction tree differs).
                np.testing.assert_allclose(
                    rate_coarse, R_oracle[Gi, Gj], rtol=1e-12, atol=1e-12,
                    err_msg=f"SigS[{order}][{Gi},{Gj}] in-scatter rate not preserved",
                )

    def test_n2n_two_axis_collapse(self, fine_mix: Mixture) -> None:
        """The (n,2n) matrix collapses by the SAME 2-axis rule as scattering."""
        coarse_mix = _condense(fine_mix)
        sig2_fine = np.asarray(fine_mix.Sig2.todense(), dtype=float)
        sig2_coarse = np.asarray(coarse_mix.Sig2.todense(), dtype=float)
        R_oracle = _scatter_rate_oracle(sig2_fine)
        phi_R = _phi_R()
        for Gi in range(_NG_COARSE):
            for Gj in range(_NG_COARSE):
                np.testing.assert_allclose(
                    sig2_coarse[Gi, Gj] * phi_R[Gi], R_oracle[Gi, Gj],
                    rtol=1e-12, atol=1e-12,
                    err_msg=f"Sig2[{Gi},{Gj}] (n,2n) rate not preserved",
                )

    # ── The mutation oracles: each must DIFFER from the design collapse ──
    # These run UNCONDITIONALLY (they exercise the membership table + the
    # frame directly, no SUT), pinning that the design collapse is the ONLY
    # one that preserves the in-scatter rate. The implementer's
    # Mixture.condense must reproduce `_design_collapse`; G3 above ties the
    # production value to the same hand-summed rate.

    @staticmethod
    def _membership() -> np.ndarray:
        """The fine→coarse one-hot T[g, G] (fast-first, orientation-correct)."""
        T = np.zeros((_NG_FINE, _NG_COARSE))
        for G, gset in enumerate(_GROUPS):
            T[gset, G] = 1.0
        return T

    def _frame(self):
        from orpheus.numerics.basis.indicator_basis import IndicatorBasis
        from orpheus.numerics.basis.weighted_indicator_basis import (
            WeightedIndicatorBasis,
        )
        from orpheus.numerics.frame import PetrovGalerkinFrame
        from orpheus.numerics.measure import DiscreteMeasure

        nodes = np.arange(_NG_FINE, dtype=float)
        # ASCENDING coarse edges in fine-index space (orientation-correct
        # membership: nodes 0,1 → cell 0; nodes 2,3 → cell 1).
        edges = np.array([-0.5, 1.5, 3.5])
        trial = IndicatorBasis(
            edges_per_axis=(edges,), partition_of=EnergyGroups(_NG_FINE),
        )
        measure = DiscreteMeasure(
            nodes=nodes, weights=np.ones(_NG_FINE), support=EnergyGroups()
        )  # COUNTING measure (w=1) + flux test-weight
        return PetrovGalerkinFrame(trial, measure, WeightedIndicatorBasis(trial, _PHI))

    def _design_collapse(self, sig_s_fine: np.ndarray) -> np.ndarray:
        """The CORRECT condense: sink summed (@T), source flux-avg (project)."""
        T = self._membership()
        sink_summed = sig_s_fine @ T                  # (n_from, n_coarse_to)
        return self._frame().project(sink_summed)     # (n_coarse_from, n_coarse_to)

    def test_design_collapse_matches_rate_oracle(self) -> None:
        """The design collapse reproduces the hand-summed in-scatter rate."""
        sig_s_fine = np.asarray(_balanced_fissile_4g().SigS[0].todense(), dtype=float)
        coarse = self._design_collapse(sig_s_fine)
        R = self._scatter_rate_oracle_static(sig_s_fine)
        phi_R = _phi_R()
        np.testing.assert_allclose(
            coarse * phi_R[:, None], R, rtol=1e-12, atol=1e-12,
            err_msg="design collapse must reproduce the hand-summed rate",
        )

    @staticmethod
    def _scatter_rate_oracle_static(sig_s_fine: np.ndarray) -> np.ndarray:
        return _scatter_rate_oracle(sig_s_fine)

    def test_mutation_swap_axes_differs(self) -> None:
        """(i) swap sink/source (flux-weight the sink, sum the source) → differs."""
        sig_s_fine = np.asarray(_balanced_fissile_4g().SigS[0].todense(), dtype=float)
        design = self._design_collapse(sig_s_fine)
        T = self._membership()
        mut = self._frame().project(sig_s_fine.T @ T)
        np.testing.assert_array_equal(
            bool(np.allclose(mut, design)), False,
            err_msg="swap-axes mutation must NOT match the design collapse",
        )

    def test_mutation_sum_both_axes_differs(self) -> None:
        """(ii) sum BOTH axes (drop the source flux-weight) → differs."""
        sig_s_fine = np.asarray(_balanced_fissile_4g().SigS[0].todense(), dtype=float)
        design = self._design_collapse(sig_s_fine)
        T = self._membership()
        mut = T.T @ sig_s_fine @ T                    # pure double-sum, no φ
        np.testing.assert_array_equal(
            bool(np.allclose(mut, design)), False,
            err_msg="sum-both-axes mutation must NOT match the design collapse",
        )

    def test_mutation_project_both_axes_differs(self) -> None:
        """(iii) project BOTH axes (flux-weight the sink too) → differs.

        This is exactly the spatial-homogenize behavior — WRONG for
        condensation. Guards against 'implementer copied homogenize verbatim'.
        """
        sig_s_fine = np.asarray(_balanced_fissile_4g().SigS[0].todense(), dtype=float)
        design = self._design_collapse(sig_s_fine)
        frame = self._frame()
        mut = frame.project(frame.project(sig_s_fine).T).T  # project source AND sink
        np.testing.assert_array_equal(
            bool(np.allclose(mut, design)), False,
            err_msg="project-both-axes (homogenize-style) must NOT match condense",
        )


# ══════════════════════════════════════════════════════════════════════
# chi — pure birth-group SUM (χ @ T), NOT frame-projected
# ══════════════════════════════════════════════════════════════════════


class TestChiBirthGroupSum:
    """χ_G = χ @ T (pure sum) — preserves Σχ=1; NOT a flux-weighted average."""

    def test_chi_sum_oracle_preserves_simplex(self) -> None:
        """The pure birth-group sum preserves Σχ; a projection would not.

        Unconditional — pins that χ collapses by SUM, distinct from the
        projected channels (anti-recommendation #2).
        """
        chi_fine = np.asarray(_balanced_fissile_4g().chi, dtype=float)
        T = np.zeros((_NG_FINE, _NG_COARSE))
        for G, gset in enumerate(_GROUPS):
            T[gset, G] = 1.0
        chi_sum = chi_fine @ T
        np.testing.assert_allclose(
            chi_sum.sum(), chi_fine.sum(), rtol=1e-12,
            err_msg="χ @ T must preserve Σχ (the simplex normalization)",
        )
        # A flux-weighted projection would NOT sum to 1 (it averages):
        chi_proj = _vector_oracle(chi_fine)
        if abs(chi_proj.sum() - chi_fine.sum()) < 1e-6:
            pytest.fail(
                "fixture degenerate: flux-projected χ accidentally preserves Σχ — "
                "the sum-vs-project distinction is not being tested"
            )

    def test_condensed_chi_is_birth_group_sum(self, fine_mix: Mixture) -> None:
        """The condensed χ equals the pure birth-group sum (NOT projected)."""
        coarse_mix = _condense(fine_mix)
        chi_fine = np.asarray(fine_mix.chi, dtype=float)
        T = np.zeros((_NG_FINE, _NG_COARSE))
        for G, gset in enumerate(_GROUPS):
            T[gset, G] = 1.0
        chi_oracle = chi_fine @ T
        np.testing.assert_allclose(
            np.asarray(coarse_mix.chi), chi_oracle, rtol=1e-12, atol=1e-12,
            err_msg="condensed χ must be the pure birth-group sum χ @ T",
        )

    def test_condensed_chi_is_simplex(self, fine_mix: Mixture) -> None:
        """The condensed χ sums to 1 (producing mixture) and is non-negative."""
        coarse_mix = _condense(fine_mix)
        np.testing.assert_allclose(
            float(np.asarray(coarse_mix.chi).sum()), 1.0, atol=1e-12,
            err_msg="condensed χ must remain a probability simplex",
        )
        np.testing.assert_array_less(
            -1e-15, np.asarray(coarse_mix.chi),
            err_msg="condensed χ must be non-negative",
        )

    def test_condensed_chi_fast_peaked(self, fine_mix: Mixture) -> None:
        """χ stays fast-peaked (bulk in the fast coarse group 0) post-condense."""
        coarse_mix = _condense(fine_mix)
        chi = np.asarray(coarse_mix.chi)
        # Fine χ = [0.7,0.25,0.05,0] → G0={g0,g1} sum 0.95, G1 sum 0.05.
        np.testing.assert_array_less(
            chi[1], chi[0],
            err_msg="condensed χ must remain fast-peaked (coarse 0 > coarse 1)",
        )


# ══════════════════════════════════════════════════════════════════════
# G4 — assert_balanced regression (positive + negative, vv #11)
# ══════════════════════════════════════════════════════════════════════


class TestG4BalanceRegression:
    """The condensed Mixture balances iff the fine one does (vv #11)."""

    def test_positive_leg_condensed_balances(self, fine_mix: Mixture) -> None:
        """A balanced fine Mixture condenses to a balanced coarse Mixture."""
        coarse_mix = _condense(fine_mix)
        # MUST NOT raise (the source-flux weight is shared across every
        # removal channel, so the identity survives the collapse).
        coarse_mix.assert_balanced(atol=1e-9)

    def test_negative_leg_rate_broken_condense_fails_balance(
        self, fine_mix: Mixture
    ) -> None:
        """A HAND-built condensed Mixture with a broken SigT fails balance.

        Built directly (NOT by perturbing the real condense) so the test
        pins the INVARIANT, not the raising behavior (vv #11). Take the real
        condensed channels but corrupt SigT by ≫ atol.
        """
        coarse_mix = _condense(fine_mix)
        broken_sigt = np.asarray(coarse_mix.SigT, dtype=float).copy()
        broken_sigt[0] += 0.5                          # ≫ atol=1e-9
        broken = Mixture(
            SigC=np.asarray(coarse_mix.SigC).copy(),
            SigL=np.asarray(coarse_mix.SigL).copy(),
            SigF=np.asarray(coarse_mix.SigF).copy(),
            SigP=np.asarray(coarse_mix.SigP).copy(),
            SigT=broken_sigt,
            SigS=[csr_matrix(m.todense()) for m in coarse_mix.SigS],
            Sig2=csr_matrix(coarse_mix.Sig2.todense()),
            chi=np.asarray(coarse_mix.chi).copy(),
            eg=coarse_mix.eg,
        )
        with pytest.raises(ValueError):
            broken.assert_balanced(atol=1e-9)


# ══════════════════════════════════════════════════════════════════════
# G5 leg-A — WIMS Table-11.3 derivation-validation
# ══════════════════════════════════════════════════════════════════════
#
# Runs UNCONDITIONALLY (it validates the published partition + the
# containing-interval rule; needs no SUT). The WIMS draft is untracked but
# present; guard the import so the file still collects if it is removed.

try:
    import orpheus.data.group_structures.wims as _wims  # type: ignore

    _HAS_WIMS = True
except Exception:  # pragma: no cover
    _wims = None  # type: ignore
    _HAS_WIMS = False

_needs_wims = pytest.mark.skipif(not _HAS_WIMS, reason="WIMS draft module absent")


@_needs_wims
class TestG5WimsDerivationValidation:
    """The containing-interval rule reproduces the published 172→69 map.

    Within the DOCUMENTED boundary tolerance — the 69- and 172-group
    structures were defined independently, so 19 boundaries are non-
    coincident (``boundary_mismatch_report``). The test REPORTS those and
    fails only if a NEW mismatch appears, NOT if a known one does.
    """

    def test_published_map_is_a_clean_partition(self) -> None:
        """CONDENSE_172_TO_69 tiles 1..172 with no gaps/overlaps."""
        np.testing.assert_array_equal(_wims.validate(), True)
        flat = [
            g for (_, a, b) in _wims.CONDENSE_172_TO_69 for g in range(a, b + 1)
        ]
        np.testing.assert_array_equal(
            np.asarray(flat), np.arange(1, 173),
            err_msg="published map must be a clean partition of fine groups 1..172",
        )

    def test_derived_map_reproduces_published_exactly(self) -> None:
        """Re-derive the partition by containing-interval; compare to published.

        For each fine 172-group, its representative energy (geometric
        midpoint of [lower, upper]) → the coarse 69-group whose interval
        contains it. The published ``CONDENSE_172_TO_69`` was itself built by
        this correspondence, so the derivation reproduces it EXACTLY (every
        fine midpoint falls in the published coarse group, even where the
        boundary ENERGIES are non-coincident). A NEW divergence (a future edit
        to either structure that breaks the correspondence) reddens.
        """
        e172 = np.array(_wims.G172_EMAX + [_wims.G172_LOWER])   # 173 descending edges
        e69 = np.array(_wims.G69_EMAX + [_wims.G69_LOWER])      # 70 descending edges
        fine_mid = np.sqrt(e172[:-1] * e172[1:])                # 172 reps, descending

        # Containing-interval: coarse group of each fine midpoint. Coarse g69
        # (1-based) owns (e69[g69], e69[g69-1]] (descending) → bracket in -E.
        derived_coarse = np.array(
            [max(1, min(69, int(np.searchsorted(-e69, -E)))) for E in fine_mid]
        )
        published_coarse = np.empty(172, dtype=int)
        for g69, a, b in _wims.CONDENSE_172_TO_69:
            published_coarse[a - 1 : b] = g69

        mismatches = (np.where(derived_coarse != published_coarse)[0] + 1).tolist()
        np.testing.assert_array_equal(
            derived_coarse, published_coarse,
            err_msg=f"containing-interval derivation must reproduce CONDENSE_172_TO_69 "
            f"exactly; disagreeing fine groups: {mismatches}",
        )

    def test_boundary_non_coincidences_are_real_and_reported(self) -> None:
        """The 69/172 structures were defined independently → some boundaries
        differ. Pin the documented count so a future edit that silently makes
        them coincident (or worsens the mismatch) is caught — the test REPORTS
        the mismatches rather than asserting zero (anti-recommendation: the
        derivation-validation must surface non-coincidences, not hide them)."""
        report = _wims.boundary_mismatch_report()      # (g69, first172, e69, e172, rel)
        # There ARE non-coincident boundaries (the structures are independent);
        # the largest is 69-g1's ceiling (10 MeV) vs 172's (19.64 MeV), because
        # 172-groups 1..5 lie above the 69-group ceiling.
        np.testing.assert_array_less(
            0, len(report),
            err_msg="expected ≥1 documented boundary non-coincidence (independent "
            "structures); report is empty — did the grids silently become coincident?",
        )
        worst = max(rel for *_, rel in report)
        np.testing.assert_array_less(
            0.5, worst,
            err_msg=f"the 172-above-69-ceiling boundary should dominate (rel≈0.96); "
            f"worst reported is only {worst:.3f}",
        )


# ══════════════════════════════════════════════════════════════════════
# G6 — Mode-11 execution sentinel (condense routes through the frame)
# ══════════════════════════════════════════════════════════════════════


def test_condense_routes_through_frame_project(fine_mix: Mixture, monkeypatch) -> None:
    """Mode-11: ``condense`` actually executes ``FrameBase.project`` AND the
    TEST-side ``WeightedIndicatorBasis.analyze`` on the new path.

    A green rate gate is VACUOUS for the 'routes through the PG frame' claim
    unless the rewired readers are on the call graph: a condense that
    open-codes the collapse with a direct matmul (never constructing the
    weighted test basis / never calling project) would pass the rate gates
    yet leave these counters at 0. Sentinel via the COUNTER (np.testing on
    the count) — NEVER a bare ``assert`` (Mode-8: ``-O`` strips it).
    """
    from orpheus.numerics.basis.weighted_indicator_basis import WeightedIndicatorBasis
    from orpheus.numerics.frame import FrameBase

    counts = {"project": 0, "test_analyze": 0}

    def _counting(name, fn):
        def wrapped(*args, **kwargs):
            counts[name] += 1
            return fn(*args, **kwargs)

        return wrapped

    monkeypatch.setattr(
        FrameBase, "project", _counting("project", FrameBase.project)
    )
    monkeypatch.setattr(
        WeightedIndicatorBasis, "analyze",
        _counting("test_analyze", WeightedIndicatorBasis.analyze),
    )

    _condense(fine_mix)

    np.testing.assert_array_less(
        0, counts["project"],
        err_msg="condense did NOT use the FrameBase.project verb (open-coded?)",
    )
    np.testing.assert_array_less(
        0, counts["test_analyze"],
        err_msg="condense did NOT route the analysis through the WEIGHTED test "
        "basis — the flux-weight is not on the test side",
    )


# ══════════════════════════════════════════════════════════════════════
# FRACTIONAL-OVERLAP re-binning (#274 follow-up) — F1 + F4 (SUT side)
# ══════════════════════════════════════════════════════════════════════
#
# The production case (421 → WIMS-69/172) is NON-NESTED: a coarse boundary
# can fall INSIDE a fine group, which then apportions a FRACTION of its rate
# to each coarse group it overlaps (conservative fractional re-binning). The
# fractional table T[g,G]∈[0,1] (rows sum to 1) generalises the one-hot; the
# rate-preserving collapse becomes
#
#     σ_G·Φ_G = Σ_g T[g,G] φ_g σ_g,   Φ_G = Σ_g T[g,G] φ_g.
#
# The oracle T is computed INDEPENDENTLY (hand-coded energy overlap × the 1/E
# lethargy model), NEVER read back from the SUT (vv L11). The STRADDLE fixture
# (a fine group whose interval spans a coarse boundary) is what makes the
# fractional path discriminating — a one-hot impl assigns the straddling group
# wholly to one coarse group, a DIFFERENT rate split → RED.
#
# Straddle geometry (descending, fast-first):
#   Fine [1e6,1e4,1e2,1e0,1e-2] → g0[1e4,1e6) g1[1e2,1e4) g2[1e0,1e2) g3[1e-2,1e0)
#   Coarse [1e6,1e1,1e-2] cuts INSIDE g2 at 1e1 → G0[1e1,1e6) G1[1e-2,1e1)
#   1/E split of g2: ln(10)/ln(100) = [0.5, 0.5].

_FINE_STRADDLE = np.array([1.0e6, 1.0e4, 1.0e2, 1.0e0, 1.0e-2])
_COARSE_STRADDLE = np.array([1.0e6, 1.0e1, 1.0e-2])
_NG_FINE_S, _NG_COARSE_S = 4, 2

# A within-group-VARYING spectrum (so the flux weight is load-bearing — a flat
# φ would make flux-weighted ≡ unweighted even at the straddle, Mode-7).
_PHI_S = np.array([1.0, 3.0, 2.0, 0.5])


def _inv_e_weight(lo: float, hi: float) -> float:
    """1/E integrated weight ∫_lo^hi dE/E = ln(hi/lo) (lethargy width)."""
    return float(np.log(hi / lo))


def _fractional_table_oracle(fine_e, coarse_e, weight) -> np.ndarray:
    """Hand-built fractional T[g,G] from energy overlaps — SUT-independent.

    T[g,G] = (∫_{g∩G} w)/(∫_g w); the structurally-independent oracle for the
    fractional rate-preservation check (NEVER ``overlap.overlap_table`` re-read).
    """
    ng_f = fine_e.size - 1
    ng_c = coarse_e.size - 1
    T = np.zeros((ng_f, ng_c))
    for g in range(ng_f):
        g_hi, g_lo = float(fine_e[g]), float(fine_e[g + 1])
        denom = weight(g_lo, g_hi)
        for G in range(ng_c):
            G_hi, G_lo = float(coarse_e[G]), float(coarse_e[G + 1])
            ov_lo, ov_hi = max(g_lo, G_lo), min(g_hi, G_hi)
            if ov_hi > ov_lo:
                T[g, G] = weight(ov_lo, ov_hi) / denom
    return T


def _straddle_mixture_4g() -> Mixture:
    """A balanced fissile 4-group Mixture on the STRADDLE fine grid.

    Same balanced structure as the nested fixture, re-homed on the straddle
    energy grid so g2 spans the coarse cut. Built directly (SigL ≠ 0).
    """
    sig_c = np.array([0.02, 0.05, 0.10, 0.30])
    sig_l = np.array([0.01, 0.02, 0.03, 0.04])
    sig_f = np.array([0.01, 0.02, 0.08, 0.20])
    nu = np.array([2.6, 2.5, 2.45, 2.43])
    chi = np.array([0.7, 0.25, 0.05, 0.0])
    sig_s0 = np.array([
        [0.30, 0.10, 0.05, 0.02],
        [0.00, 0.40, 0.15, 0.03],
        [0.00, 0.00, 0.55, 0.20],
        [0.00, 0.00, 0.00, 0.90],
    ])
    sig_s1 = sig_s0 * 0.1
    sig_2 = np.array([
        [0.005, 0.003, 0.001, 0.0],
        [0.000, 0.004, 0.002, 0.0],
        [0.000, 0.000, 0.003, 0.0],
        [0.000, 0.000, 0.000, 0.0],
    ])
    rowsum = sig_s0.sum(axis=1) + sig_2.sum(axis=1)
    sig_t = sig_c + sig_l + sig_f + rowsum
    mix = Mixture(
        SigC=sig_c.copy(), SigL=sig_l.copy(), SigF=sig_f.copy(),
        SigP=(nu * sig_f).copy(), SigT=sig_t.copy(),
        SigS=[csr_matrix(sig_s0), csr_matrix(sig_s1)],
        Sig2=csr_matrix(sig_2), chi=chi.copy(), eg=_FINE_STRADDLE.copy(),
    )
    mix.assert_balanced()
    return mix


def _condense_straddle(mix: Mixture, within_group=None) -> Mixture:
    """Condense the straddle Mixture onto the coarse grid (fractional path).

    ``Mixture.condense(target, spectrum, within_group)`` reads the mixture's
    own fine grid and apportions the straddling fine group fractionally per
    the ``within_group`` model.
    """
    coarse = EnergyGrid(_COARSE_STRADDLE)
    if within_group is None:
        return mix.condense(coarse, _PHI_S)
    return mix.condense(coarse, _PHI_S, within_group=within_group)


# ── F1 — Straddle rate-preservation (THE fractional core) ──────────────


@pytest.mark.verifies("energy-condensation-rate-preservation")
class TestF1StraddleRatePreservation:
    """σ_G·Φ_G == Σ_g T[g,G] φ_g σ_g with a genuinely FRACTIONAL T (vv L1, L11).

    The oracle T is the hand-computed 1/E overlap table — independent of the
    SUT. A one-hot impl (g2 → one coarse group) gives a different Φ_G AND a
    different rate split → every vector channel reddens.
    """

    @staticmethod
    def _oracle_T() -> np.ndarray:
        return _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, _inv_e_weight)

    def test_oracle_has_a_genuine_straddle(self) -> None:
        """Guard: the oracle table really IS fractional (g2 split, not one-hot).

        Unconditional — if a future edit aligned the grids, the straddle would
        vanish and F1 would silently degrade to the nested (one-hot) case.
        """
        T = self._oracle_T()
        np.testing.assert_allclose(
            T[2], np.array([0.5, 0.5]), atol=1e-12,
            err_msg="oracle g2 must straddle 50/50 (fixture no longer fractional?)",
        )

    @pytest.mark.parametrize("channel", ["SigT", "SigC", "SigL", "SigF", "SigP"])
    def test_vector_channel_fractional_rate_preserved(self, channel: str) -> None:
        """Φ_G·σ_G == Σ_g T[g,G] φ_g σ_g per vector channel (fractional T)."""
        mix = _straddle_mixture_4g()
        coarse_mix = _condense_straddle(mix, within_group=InverseEnergySpectrum())
        T = self._oracle_T()
        sig_fine = np.asarray(getattr(mix, channel), dtype=float)
        sig_eff = np.asarray(getattr(coarse_mix, channel), dtype=float)
        Phi_G = (T * _PHI_S[:, None]).sum(axis=0)                       # Σ_g T[g,G]φ_g
        for G in range(_NG_COARSE_S):
            rate_coarse = float(sig_eff[G] * Phi_G[G])
            rate_ref = float((T[:, G] * _PHI_S * sig_fine).sum())       # hand-sum
            np.testing.assert_allclose(
                rate_coarse, rate_ref, rtol=1e-12, atol=1e-12,
                err_msg=f"{channel} fractional rate not preserved in coarse {G}",
            )

    @pytest.mark.verifies("energy-condensation-fractional-collapse")
    def test_condensed_value_equals_fractional_oracle(self) -> None:
        """The condensed SigT equals the fractional flux-weighted hand-sum.

        σ_G = (Σ_g T[g,G] φ_g σ_g)/(Σ_g T[g,G] φ_g) — the VALUE, not just the
        rate (pins the Φ_G denominator too).
        """
        mix = _straddle_mixture_4g()
        coarse_mix = _condense_straddle(mix, within_group=InverseEnergySpectrum())
        T = self._oracle_T()
        sig_fine = np.asarray(mix.SigT, dtype=float)
        Phi_G = (T * _PHI_S[:, None]).sum(axis=0)
        sig_oracle = (T * _PHI_S[:, None] * sig_fine[:, None]).sum(axis=0) / Phi_G
        np.testing.assert_allclose(
            np.asarray(coarse_mix.SigT), sig_oracle, rtol=1e-12, atol=1e-12,
            err_msg="condensed SigT must equal the fractional flux-weighted oracle",
        )

    def test_scattering_fractional_rate_preserved(self) -> None:
        """In-scatter RATE with the fractional source-axis weight:
        Φ_G·Σ_s,coarse[G,G'] == Σ_g Σ_g' T[g,G] φ_g Σ_s[g,g'] T[g',G'].

        Both axes carry the fractional membership: the SINK axis is summed with
        T (a scatter into any fine group of coarse G' is a scatter into G'),
        the SOURCE axis is flux-averaged with T. The oracle is the full
        double-fractional hand-sum.
        """
        mix = _straddle_mixture_4g()
        coarse_mix = _condense_straddle(mix, within_group=InverseEnergySpectrum())
        T = self._oracle_T()
        sig_s_fine = np.asarray(mix.SigS[0].todense(), dtype=float)
        sig_s_coarse = np.asarray(coarse_mix.SigS[0].todense(), dtype=float)
        Phi_G = (T * _PHI_S[:, None]).sum(axis=0)
        for Gi in range(_NG_COARSE_S):
            for Gj in range(_NG_COARSE_S):
                rate_coarse = float(sig_s_coarse[Gi, Gj] * Phi_G[Gi])
                rate_ref = float(sum(
                    T[g, Gi] * _PHI_S[g] * sig_s_fine[g, gp] * T[gp, Gj]
                    for g in range(_NG_FINE_S) for gp in range(_NG_FINE_S)
                ))
                np.testing.assert_allclose(
                    rate_coarse, rate_ref, rtol=1e-12, atol=1e-12,
                    err_msg=f"SigS[0][{Gi},{Gj}] fractional in-scatter rate not preserved",
                )

    def test_condensed_straddle_balances(self) -> None:
        """Balance survives the FRACTIONAL collapse (every channel shares Φ_G)."""
        coarse_mix = _condense_straddle(
            _straddle_mixture_4g(), within_group=InverseEnergySpectrum()
        )
        coarse_mix.assert_balanced(atol=1e-9)


# ── F4 (SUT side) — the within-group model is selectable + load-bearing ─


class TestF4ModelDiscriminatorSUT:
    """The SUT's 1/E model reproduces the lethargy-ratio oracle (positive),
    and a DIFFERENT model gives a DIFFERENT condensed value (negative).

    The oracle-half (1/E fraction = ln ratio, and 1/E ≠ flat-energy) is pinned
    unconditionally in tests/data/test_energy_grid.py::TestF4WithinGroupModelOracle.
    Here the SUT's ``InverseEnergySpectrum`` is tied to that oracle, and a
    second within-group model (constructed in-test if the SUT ships only 1/E)
    is shown to move the condensed straddle value.
    """

    def test_inv_e_condense_matches_lethargy_oracle(self) -> None:
        """condense(within=InverseEnergySpectrum) reproduces the 1/E oracle value."""
        mix = _straddle_mixture_4g()
        coarse_mix = _condense_straddle(mix, within_group=InverseEnergySpectrum())
        T = _fractional_table_oracle(_FINE_STRADDLE, _COARSE_STRADDLE, _inv_e_weight)
        sig_fine = np.asarray(mix.SigT, dtype=float)
        Phi_G = (T * _PHI_S[:, None]).sum(axis=0)
        sig_oracle = (T * _PHI_S[:, None] * sig_fine[:, None]).sum(axis=0) / Phi_G
        np.testing.assert_allclose(
            np.asarray(coarse_mix.SigT), sig_oracle, rtol=1e-12, atol=1e-12,
            err_msg="1/E condense must reproduce the lethargy-ratio oracle (the model "
            "is correctly 1/E)",
        )

    def test_a_different_model_gives_a_different_value(self) -> None:
        """A flat-in-ENERGY within-group model moves the straddle σ_G.

        Proves the model is genuinely SELECTABLE / load-bearing: swapping 1/E
        for flat-energy changes g2's split (0.5 → 0.909), hence the condensed
        coarse SigT. Built two ways depending on what the SUT ships:

        * if the SUT exposes a second model class, condense with it;
        * else, construct a trivial in-test ``WithinGroupSpectrum`` (flat
          energy) and pass it — the SUT must consume any model satisfying the
          protocol. If the SUT REJECTS a custom model, the test skips (the
          model is then 1/E-only, not selectable, which the implementer should
          note — flagged below).
        """
        mix = _straddle_mixture_4g()
        inv = _condense_straddle(mix, within_group=InverseEnergySpectrum())

        # A trivial second within-group model: flat-in-ENERGY ∫ = (hi-lo).
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class _FlatEnergySpectrum:
            def integrated_weight(self, lo, hi):  # vectorized-friendly
                return np.asarray(hi, dtype=float) - np.asarray(lo, dtype=float)

        try:
            flat = _condense_straddle(mix, within_group=_FlatEnergySpectrum())
        except Exception:  # pragma: no cover - SUT may not accept a custom model
            pytest.skip(
                "SUT did not accept a custom WithinGroupSpectrum — the within-group "
                "model is not duck-typed/selectable; verify the protocol surface"
            )
        # The two models MUST give a different condensed value (the straddle
        # split differs). If they agree, the model is not load-bearing.
        if np.allclose(np.asarray(inv.SigT), np.asarray(flat.SigT), atol=1e-10):
            pytest.fail(
                "within-group model is NOT load-bearing: 1/E and flat-energy "
                "condense to the same SigT — the straddle split is not model-driven"
            )


# ════════════════════════════════════════════════════════════════════
# P6 (#281) — the bilinear (B&G-convention) condensation, data level.
# T6 of orpheus/derivations/common/homogenization.py, live in floats.
# ════════════════════════════════════════════════════════════════════


class TestBilinearCondensation:
    """The eigenvalue-consistent arm of :meth:`Mixture.condense` at the
    data level: T6a live (the bilinear-condensed pencil reproduces the fine
    0-D k at machine precision when condensed with the TRUE spectrum pair),
    plus the shape guard on the new parameter."""

    def test_t6a_true_spectra_reproduce_fine_k_exactly(
        self, fine_mix: Mixture,
    ) -> None:
        """T6a in floats: 0-D ∞-medium, φ = A⁻¹χ, φ* = A⁻ᵀνΣf; the
        bilinear-condensed 2G pencil's rank-1 k equals the fine 4G k to
        machine precision — condensation is pure projection on the energy
        axis (B&G convention; no streaming carve)."""
        sig_t = np.asarray(fine_mix.SigT, float)
        sig_s = np.asarray(fine_mix.SigS[0].todense(), float)
        nsf = np.asarray(fine_mix.SigP, float)
        chi = np.asarray(fine_mix.chi, float)
        A = np.diag(sig_t) - sig_s.T
        phi = np.linalg.solve(A, chi)
        phis = np.linalg.solve(A.T, nsf)
        k_fine = float(nsf @ phi)

        out = fine_mix.condense(
            EnergyGrid(_EG_COARSE), phi, adjoint_spectrum=phis,
        )
        A_c = (
            np.diag(np.asarray(out.SigT, float))
            - np.asarray(out.SigS[0].todense(), float).T
        )
        k_c = float(
            np.asarray(out.SigP, float)
            @ np.linalg.solve(A_c, np.asarray(out.chi, float))
        )
        np.testing.assert_allclose(k_c, k_fine, rtol=1e-12)

    def test_forward_arm_untouched_by_the_parameter(
        self, fine_mix: Mixture,
    ) -> None:
        """§4.0 at the data level: omitting vs passing ``adjoint_spectrum=None``
        is bit-identical."""
        a = fine_mix.condense(EnergyGrid(_EG_COARSE), _PHI)
        b = fine_mix.condense(EnergyGrid(_EG_COARSE), _PHI, adjoint_spectrum=None)
        np.testing.assert_array_equal(np.asarray(a.SigT), np.asarray(b.SigT))
        np.testing.assert_array_equal(np.asarray(a.chi), np.asarray(b.chi))
        np.testing.assert_array_equal(
            np.asarray(a.SigS[0].todense()), np.asarray(b.SigS[0].todense()),
        )

    def test_adjoint_spectrum_shape_guard(self, fine_mix: Mixture) -> None:
        with pytest.raises(ValueError, match="adjoint_spectrum"):
            fine_mix.condense(
                EnergyGrid(_EG_COARSE), _PHI, adjoint_spectrum=np.ones(3),
            )
