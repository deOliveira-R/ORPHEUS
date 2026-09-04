"""S10c — Mixture cross-section balance invariant gate.

``Mixture`` gains a SCOPED balance invariant (NOT ``__post_init__``):

    Mixture.balance_residual -> (NG,)  # |SigT - derived| per group
    Mixture.assert_balanced(atol=1e-9) -> None  # raises ValueError if max > atol

where the balance is the definitional removal identity (per group)::

    SigT == SigC + SigL + SigF + rowsum(SigS[0]) + rowsum(Sig2)

This is VERBATIM the SigT-derivation line in ``compute_macro_xs`` (the
production code derives SigT this way, so a real mixture ALWAYS balances —
the guard catches a FUTURE derivation bug, it does not break real data).
NOTE: ``SigP`` (production = νΣf) is a fission-SOURCE multiplier, NOT a
removal channel — it is deliberately ABSENT from the identity. ``SigF``
(fission removal) IS present. Confusing the two is the Mode-2 role-swap
trap; gate 1 has an explicit leg pinning the SigP exclusion.

SCOPED design (the L20 ruling): the law lives on the TYPE but is INVOKED
selectively. ``compute_macro_xs`` calls it (free real-path regression
guard); the verification harness's PHYSICAL synthetic tables are swept
(the value gate / typo-catcher). There is NO unconditional ``__post_init__``
enforcement — the Atalay criticality-parameter Mixtures, the structural
test scaffolds, and the billiard SigP-carrier are INTENTIONALLY imbalanced
and build directly, bypassing ``assert_balanced``.

vv anti-pattern #11 (L11): gate 1 has BOTH a positive leg (balanced →
no raise) AND a negative leg (imbalanced → raises). The broken fixture is
built BY HAND with SigT off by ≫ atol — never by perturbing the production
builder — so the test cannot be satisfied by ``assert_balanced`` re-spelling
a wrong rowsum.

Mode-8: every assertion is ``np.testing.*`` / ``pytest.raises`` /
``pytest.fail`` (the suite runs ``-O``; bare ``assert`` is stripped).

Foundation tests (software invariant — no theory ``:label:``).

PRE-IMPLEMENTATION: the ``Mixture.assert_balanced`` / ``balance_residual``
API does not exist yet. Gates that exercise it are guarded by
``_HAS_BALANCE_API`` and ``pytest.mark.skipif`` so this file collects clean
and does NOT false-green. Remove the guard when the production code lands.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture

pytestmark = pytest.mark.foundation

# Flip-to-live guard: the SCOPED balance API does not exist pre-impl.
_HAS_BALANCE_API = hasattr(Mixture, "assert_balanced") and hasattr(
    Mixture, "balance_residual"
)
_needs_api = pytest.mark.skipif(
    not _HAS_BALANCE_API,
    reason="S10c: Mixture.assert_balanced / balance_residual not yet implemented",
)

_ATOL = 1e-9


# ══════════════════════════════════════════════════════════════════════
# Hand-built fixtures (L11 structural independence — NOT via any factory)
# ══════════════════════════════════════════════════════════════════════


def _mixture(
    *,
    sig_t: np.ndarray,
    sig_c: np.ndarray,
    sig_l: np.ndarray,
    sig_f: np.ndarray,
    sig_p: np.ndarray,
    sig_s0: np.ndarray,
    sig_2: np.ndarray,
    chi: np.ndarray,
) -> Mixture:
    """Build a Mixture DIRECTLY from N-group arrays.

    Built without ``make_mixture`` / ``compute_macro_xs`` so the test does
    not depend on the builders the balance guard protects (which would make
    the test self-referential). All removal channels are explicit so a leg
    can dial SigL / Sig2 to NON-zero and exercise those terms of the
    identity (every PHYSICAL fixture in the codebase has SigL = Sig2 = 0).
    """
    return Mixture(
        SigC=sig_c.copy(),
        SigL=sig_l.copy(),
        SigF=sig_f.copy(),
        SigP=sig_p.copy(),
        SigT=sig_t.copy(),
        SigS=[csr_matrix(sig_s0)],
        Sig2=[csr_matrix(sig_2)],
        chi=chi.copy(),
    )


def _balanced_2g_with_sig2_and_sigl() -> Mixture:
    """A hand-balanced 2-group Mixture with NON-zero SigL AND Sig2.

    This is the ONLY fixture that activates the ``+ SigL`` and
    ``+ rowsum(Sig2)`` terms of the identity. Per-group SigT is set to make
    the balance hold EXACTLY:

        SigT_g = SigC_g + SigL_g + SigF_g + rowsum(SigS0)_g + rowsum(Sig2)_g
    """
    sig_c = np.array([0.10, 0.20])
    sig_l = np.array([0.03, 0.05])  # (n,alpha) — NON-zero, exercised here
    sig_f = np.array([0.01, 0.08])
    sig_s0 = np.array([[0.40, 0.10], [0.00, 0.90]])
    sig_2 = np.array([[0.02, 0.01], [0.00, 0.04]])  # (n,2n) — NON-zero
    rowsum_s0 = sig_s0.sum(axis=1)
    rowsum_2 = sig_2.sum(axis=1)
    sig_t = sig_c + sig_l + sig_f + rowsum_s0 + rowsum_2
    # χ is a simplex on the producing group set; νΣf carried in SigP.
    return _mixture(
        sig_t=sig_t,
        sig_c=sig_c,
        sig_l=sig_l,
        sig_f=sig_f,
        sig_p=np.array([0.025, 0.20]),
        sig_s0=sig_s0,
        sig_2=sig_2,
        chi=np.array([1.0, 0.0]),
    )


def _imbalanced_2g(*, delta: float = 0.5) -> Mixture:
    """A hand-imbalanced 2-group Mixture: SigT off by ``delta`` (≫ atol).

    The imbalance is injected ONLY into SigT (not a removal component) so the
    test reads as 'the balance is broken' and pins the invariant CLAIM, not
    merely that the method can raise. ``delta`` defaults to ~0.5 (≈ 9 orders
    above atol).
    """
    sig_c = np.array([0.10, 0.20])
    sig_l = np.zeros(2)
    sig_f = np.array([0.01, 0.08])
    sig_s0 = np.array([[0.40, 0.10], [0.00, 0.90]])
    sig_2 = csr_matrix((2, 2)).toarray()
    rowsum_s0 = sig_s0.sum(axis=1)
    sig_t = sig_c + sig_l + sig_f + rowsum_s0
    sig_t[0] += delta  # break the identity in group 0 by a macroscopic amount
    return _mixture(
        sig_t=sig_t,
        sig_c=sig_c,
        sig_l=sig_l,
        sig_f=sig_f,
        sig_p=np.array([0.025, 0.20]),
        sig_s0=sig_s0,
        sig_2=sig_2,
        chi=np.array([1.0, 0.0]),
    )


# ══════════════════════════════════════════════════════════════════════
# Gate 1 — assert_balanced intrinsic-property (vv #11, L11)
# ══════════════════════════════════════════════════════════════════════


@_needs_api
class TestAssertBalancedIntrinsic:
    """The DEFINING balance law: positive + negative + a Sig2/SigL leg."""

    def test_balanced_does_not_raise_positive_leg(self) -> None:
        """A hand-balanced Mixture passes assert_balanced (MUST NOT raise)."""
        mix = _balanced_2g_with_sig2_and_sigl()
        # MUST NOT raise — bare call (if it raised, pytest reports the error).
        mix.assert_balanced(atol=_ATOL)

    def test_balanced_residual_near_zero_positive_leg(self) -> None:
        """balance_residual ~0 for a balanced Mixture (per-group, ≤ FP floor)."""
        mix = _balanced_2g_with_sig2_and_sigl()
        np.testing.assert_array_less(
            np.asarray(mix.balance_residual),
            1e-12,
            err_msg="balanced Mixture residual should be at FP floor",
        )

    def test_imbalanced_raises_negative_leg(self) -> None:
        """A hand-imbalanced Mixture (SigT off ≫ atol) raises ValueError."""
        mix = _imbalanced_2g(delta=0.5)
        with pytest.raises(ValueError):
            mix.assert_balanced(atol=_ATOL)

    def test_imbalanced_residual_reports_right_magnitude(self) -> None:
        """balance_residual reports the injected defect (~delta) in the right group."""
        delta = 0.5
        mix = _imbalanced_2g(delta=delta)
        res = np.asarray(mix.balance_residual)
        np.testing.assert_allclose(
            res[0], delta, rtol=1e-12,
            err_msg=f"group-0 residual={res[0]} should equal injected delta={delta}",
        )
        np.testing.assert_allclose(
            res[1], 0.0, atol=1e-12,
            err_msg="group-1 residual should be ~0 (defect injected only in group 0)",
        )

    def test_sig2_and_sigl_terms_are_actually_tested(self) -> None:
        """The Sig2/SigL leg activates BOTH the +SigL and +rowsum(Sig2) terms.

        Guard against a fixture whose SigL / Sig2 are silently zero (which
        would leave those two identity terms UNTESTED — every physical
        fixture in the codebase has them zero, so gate 2 is blind to them).
        Then confirm the balanced Mixture passes (the terms are correctly
        included in assert_balanced's recompute).
        """
        mix = _balanced_2g_with_sig2_and_sigl()
        np.testing.assert_array_less(
            0.0, np.asarray(mix.SigL).max(),
            err_msg="SigL leg must be NON-zero to exercise the +SigL term",
        )
        np.testing.assert_array_less(
            0.0, float(np.array(mix.Sig2[0].sum(axis=1)).ravel().max()),
            err_msg="Sig2 leg must be NON-zero to exercise the +rowsum(Sig2) term",
        )
        mix.assert_balanced(atol=_ATOL)  # MUST NOT raise

    def test_sigp_excluded_from_identity_anti_mode2(self) -> None:
        """SigP (production = νΣf) is NOT part of the balance identity.

        Anti Mode-2 role-swap (SigF↔SigP): build a balanced Mixture, then
        scale SigP arbitrarily large. The balance MUST be unaffected — if
        assert_balanced erroneously included SigP, this would red.
        """
        mix = _balanced_2g_with_sig2_and_sigl()
        mix.SigP = mix.SigP * 1000.0  # production has no removal meaning here
        mix.assert_balanced(atol=_ATOL)  # MUST still pass — SigP not in identity
        np.testing.assert_array_less(
            np.asarray(mix.balance_residual), 1e-12,
            err_msg="SigP must not enter balance_residual",
        )


# ══════════════════════════════════════════════════════════════════════
# Gate 2 — Physical-table sweep (THE VALUE GATE / typo-catcher)
# ══════════════════════════════════════════════════════════════════════
#
# Mode-7 honest scope: this pins the PHYSICAL synthetic tables only —
# xs_library A/B/C/D × 1g/2g/4g, the Sood LA13511 registry, and the
# homogeneous derive builders. The Atalay criticality-parameter Mixtures,
# the structural scaffolds, and the billiard carrier are deliberately
# imbalanced and are EXCLUDED (Atalay is structurally absent from
# LA13511_CASES — it lives in ATALAY_ALL_CASES — so no filter is needed).


def _xs_library_mixtures():
    from orpheus.derivations.common.xs_library import get_mixture

    return [
        pytest.param(get_mixture(region, ng), id=f"xslib-{region}-{ng}")
        for region in "ABCD"
        for ng in ("1g", "2g", "4g")
    ]


def _la13511_mixtures():
    from orpheus.derivations.continuous.sood_registry import LA13511_CASES

    out = []
    for case_id, case in LA13511_CASES.items():
        for mat_id, mix in case.materials.items():
            out.append(pytest.param(mix, id=f"la13511-{case_id}-m{mat_id}"))
    return out


def _homogeneous_mixtures():
    from orpheus.derivations.continuous.analytical.homogeneous import (
        derive_1g,
        derive_2g,
        derive_4g,
    )

    out = []
    for fn in (derive_1g, derive_2g, derive_4g):
        case = fn()
        for mat_id, mix in case.materials.items():
            out.append(pytest.param(mix, id=f"homog-{fn.__name__}-m{mat_id}"))
    return out


@_needs_api
class TestPhysicalTableSweep:
    """Every PHYSICAL synthetic Mixture must balance. The typo-catcher."""

    @pytest.mark.parametrize("mix", _xs_library_mixtures())
    def test_xs_library_balances(self, mix: Mixture) -> None:
        mix.assert_balanced(atol=_ATOL)

    @pytest.mark.parametrize("mix", _la13511_mixtures())
    def test_la13511_balances(self, mix: Mixture) -> None:
        mix.assert_balanced(atol=_ATOL)

    @pytest.mark.parametrize("mix", _homogeneous_mixtures())
    def test_homogeneous_balances(self, mix: Mixture) -> None:
        mix.assert_balanced(atol=_ATOL)


# ══════════════════════════════════════════════════════════════════════
# Gate 3 — Real-path guard (compute_macro_xs invokes assert_balanced)
# ══════════════════════════════════════════════════════════════════════
#
# The guard is a regression CATCHER, not a breaker: a real isotope-built
# mixture must construct WITHOUT raising (compute_macro_xs derives SigT via
# the exact identity, so it always balances). A minimal real-isotope mixture
# suffices; a full recipe is too heavy for a foundation test.


@_needs_api
@pytest.mark.slow
class TestRealPathGuard:
    """A real GENDF/load_isotope-built mixture constructs without raising."""

    def test_minimal_real_isotope_mixture_constructs(self) -> None:
        """A real isotope-built mixture constructs (the guard catches, not breaks).

        ``compute_macro_xs`` invokes ``assert_balanced`` after building the
        Mixture; because it DERIVES SigT via the balance identity, a real
        mixture always balances. This pins that the in-function guard does
        NOT raise on real production data.
        """
        from pathlib import Path

        micro_xs = (
            Path(__file__).resolve().parents[2] / "orpheus" / "data" / "micro_xs"
        )
        if not (micro_xs / "U_235.h5").exists():
            pytest.skip("real U_235.h5 micro-XS library absent")

        from orpheus.data.macro_xs.mixture import compute_macro_xs
        from orpheus.data.micro_xs import load_isotope

        # Minimal real-isotope mixture (1 fissile + 1 moderator); the
        # compute_macro_xs internal assert_balanced MUST NOT raise.
        u235 = load_isotope("U_235", 900)
        o16 = load_isotope("O_016", 900)
        mix = compute_macro_xs([u235, o16], np.array([0.03, 2.0]), fissile_indices=[0])
        # belt-and-braces: the derived-SigT mixture balances to FP floor.
        mix.assert_balanced()
        np.testing.assert_array_less(
            np.asarray(mix.balance_residual),
            _ATOL,
            err_msg="real derived-SigT mixture must balance well within atol",
        )


# ══════════════════════════════════════════════════════════════════════
# Gate 4 — Exemption integrity (the scoping genuinely leaves them alone)
# ══════════════════════════════════════════════════════════════════════


class TestExemptionIntegrity:
    """Intentionally-imbalanced Mixtures still CONSTRUCT (no __post_init__ raise)."""

    def test_direct_construction_of_imbalanced_does_not_raise(self) -> None:
        """An imbalanced Mixture builds DIRECTLY — pins 'no unconditional enforcement'.

        This must hold WHETHER OR NOT the balance API exists (the whole point
        of the SCOPED design is that __post_init__ stays silent on the
        balance). NOT guarded by _needs_api.
        """
        # MUST NOT raise on construction.
        mix = _imbalanced_2g(delta=0.5)
        np.testing.assert_array_equal(
            np.asarray(mix.SigT).shape, (2,),
            err_msg="imbalanced Mixture must construct fine (no __post_init__ balance raise)",
        )

    def test_atalay_criticality_mixture_constructs(self) -> None:
        """Atalay νΣf=(c-1)Σt c>1 mixture (residual≈SigF) still builds.

        Atalay is in ATALAY_ALL_CASES, NOT in LA13511_CASES, so gate 2 never
        touches it. Confirm it constructs and IS imbalanced (the exemption is
        real, not an accident of a balanced fixture).
        """
        from orpheus.derivations.continuous.sood_registry.atalay1997 import (
            _mix_iso_at_c,
        )

        mix = _mix_iso_at_c(1.30)  # MUST NOT raise on construction
        rowsum_s0 = np.array(mix.SigS[0].sum(axis=1)).ravel()
        rowsum_2 = np.array(mix.Sig2[0].sum(axis=1)).ravel()
        derived = mix.SigC + mix.SigL + mix.SigF + rowsum_s0 + rowsum_2
        residual = float(np.max(np.abs(mix.SigT - derived)))
        # It is INTENTIONALLY imbalanced (the c>1 criticality encoding).
        np.testing.assert_array_less(
            1e-3, residual,
            err_msg=f"Atalay c=1.30 must be (intentionally) imbalanced; residual={residual}",
        )

    def test_atalay_not_in_la13511_registry(self) -> None:
        """Structural exemption: no Atalay case_id appears in LA13511_CASES.

        This is WHY gate 2 needs no defensive filter — the swept registry
        simply does not contain the exempt set. If a future edit adds Atalay
        to LA13511_CASES, this reddens BEFORE gate 2 starts failing.
        """
        from orpheus.derivations.continuous.sood_registry import LA13511_CASES

        atalay_ids = [cid for cid in LA13511_CASES if "atalay" in cid.lower()]
        np.testing.assert_array_equal(
            len(atalay_ids), 0,
            err_msg=f"Atalay must stay out of LA13511_CASES; found {atalay_ids}",
        )


# ══════════════════════════════════════════════════════════════════════
# Gate 5 — atol band justification (residuals separate balancers/violators)
# ══════════════════════════════════════════════════════════════════════


class TestAtolBand:
    """atol=1e-9 sits in the dead-band between balancers (≤4.4e-16) and
    violators (≥1.6e-2). This test pins the band so a future tightening of
    the synthetic tables (or a regressed Atalay encoding) is caught."""

    @staticmethod
    def _residual(mix: Mixture) -> float:
        rowsum_s0 = np.array(mix.SigS[0].sum(axis=1)).ravel()
        rowsum_2 = np.array(mix.Sig2[0].sum(axis=1)).ravel()
        derived = mix.SigC + mix.SigL + mix.SigF + rowsum_s0 + rowsum_2
        return float(np.max(np.abs(mix.SigT - derived)))

    def test_physical_tables_below_1e_12(self) -> None:
        from orpheus.derivations.common.xs_library import get_mixture

        for region in "ABCD":
            for ng in ("1g", "2g", "4g"):
                r = self._residual(get_mixture(region, ng))
                np.testing.assert_array_less(
                    r, 1e-12,
                    err_msg=f"physical {region}-{ng} residual={r} should be at FP floor",
                )

    def test_atalay_violator_above_1e_3(self) -> None:
        from orpheus.derivations.continuous.sood_registry.atalay1997 import (
            _mix_iso_at_c,
        )

        r = self._residual(_mix_iso_at_c(1.30))
        np.testing.assert_array_less(
            1e-3, r,
            err_msg=f"Atalay c=1.30 residual={r} should be ≫ atol (intentional)",
        )

    def test_atol_sits_in_dead_band(self) -> None:
        """1e-9 is above every balancer floor and below every violator."""
        np.testing.assert_array_less(1e-12, _ATOL)  # above FP floor
        np.testing.assert_array_less(_ATOL, 1e-3)   # below the smallest violator
