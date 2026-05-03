r"""Foundation tests for the F_N method LA-13511 :math:`k_\infty` cases.

This file gates the first slice of the
:mod:`orpheus.derivations.continuous.fn_method` package — the two
LA-13511 cases that need only Sood's own equations
(``PUa-1-0-IN`` and ``PU-2-0-IN``):

* **Branch-1 SymPy gate** — every ``derive_*()`` function in
  :mod:`...fn_method.origins.k_inf_derivations` has exactly one
  ``@pytest.mark.foundation`` test pinning its symbolic identity.
  These are pure algebra; they do NOT depend on any numerical
  benchmark value.
* **Branch-2 reference-value gate** — :func:`compute_kinf_*` reproduces
  the published LA-13511 reference values to 5 significant figures
  (the precision Sood et al. reports the cases at).
* **Branch-1 ↔ Branch-2 agreement gate** — the symbolic closed forms
  evaluated in numpy must match the SymPy ``sp.lambdify`` result, and
  must match each other across the 1G/2G/MG entry points in their
  domain of overlap.
* **Cross-implementation gate** — the F_N module's ``compute_kinf_mg``
  must agree with the existing
  :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` to
  machine precision. This is a structurally-independent check: F_N
  evaluates the closed form Sood Eq 76 directly, while
  ``kinf_homogeneous`` solves the dominant eigenvalue problem of
  :math:`A^{-1}F`. They share only ``numpy`` (above the trusted-library
  line) — a disagreement would point at a real bug in one or the other.

Tag policy: all tests are marked ``@pytest.mark.foundation`` because
the underlying claims are algebra/software-invariants of the
infinite-medium k_inf machinery — there is no spatial/angular
discretisation, hence no L0/L1/L2 ladder applies. The test count
equals 1 per ``derive_*`` + 1 per Branch-2 entry point + the
agreement and reduction tests.

References
----------

* Sood, Forster, Parsons (1999), LA-13511, Appendix A.
* Literature memo:
  ``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``.
* Closeout memo:
  ``.claude/agent-memory/method-implementer/fn_method_kinf_first_slice.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_spectrum_homogeneous,
    kinf_homogeneous,
)
from orpheus.derivations.continuous.sood_registry import (
    PU_2_0_IN,
    PUA_1_0_IN,
)
from orpheus.derivations.continuous.fn_method.multi_group import (
    compute_flux_ratio_2g_no_upscatter,
    compute_kinf_1g,
    compute_kinf_2g_general,
    compute_kinf_2g_no_upscatter,
    compute_kinf_mg,
)
from orpheus.derivations.continuous.fn_method.origins import (
    derive_kinf_1g_eq_19,
    derive_kinf_1g_eq_20_simplifies_to_eq_19,
    derive_kinf_2g_general_from_matrix,
    derive_kinf_2g_no_upscatter,
    derive_kinf_mg_matrix_form,
    derive_kinf_mg_reduces_to_1g,
    derive_phi_ratio_2g_no_upscatter,
)


# ═══════════════════════════════════════════════════════════════════
# Branch-1 SymPy gate — one test per derive_*() function
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_fn1_1_kinf_1g_eq_19():
    r"""V_fn1.1 — 1G balance equation reduces to :math:`k_\infty = \nu\Sigma_f / \Sigma_a`."""
    result = derive_kinf_1g_eq_19()
    assert result["pass"], (
        f"V_fn1.1 failed: derived k = {result['k_derived']}, "
        f"expected = {result['k_eq19']}, diff = {result['diff']}"
    )


@pytest.mark.foundation
def test_v_fn1_2_kinf_eq_20_simplifies_to_eq_19():
    """V_fn1.2 — Sood Eq 20 algebraically equals Eq 19 (the c factor cancels)."""
    result = derive_kinf_1g_eq_20_simplifies_to_eq_19()
    assert result["pass"], (
        f"V_fn1.2 failed: Eq 20 - Eq 19 = {result['diff']} (expected 0)"
    )


@pytest.mark.foundation
def test_v_fn2_1_kinf_2g_general_from_matrix():
    """V_fn2.1 — 2G general k_inf from det(M)=0 matches Eq 29 in no-upscatter limit;
    Sood's Eq 28 typo is confirmed algebraically."""
    result = derive_kinf_2g_general_from_matrix()
    assert result["pass_eq29_match"], (
        "V_fn2.1: derived 2G k_inf does not reduce to Eq 29 at no-upscatter; "
        f"diff = {result['diff_eq29']}"
    )
    assert result["pass_eq28_corrected"], (
        "V_fn2.1: derived 2G k_inf does not match the *corrected* Eq 28; "
        f"diff = {result['diff_eq28_corrected']}"
    )
    assert result["pass_eq28_typo_confirmed"], (
        "V_fn2.1: Eq 28 as printed appears to algebraically match the "
        "derived form — but it should NOT (we expect a typo). "
        f"diff_printed = {result['diff_eq28_printed']}"
    )


@pytest.mark.foundation
def test_v_fn2_2_kinf_2g_no_upscatter_makes_det_zero():
    """V_fn2.2 — Substituting Eq 29 into det(M(no upscatter)) gives zero."""
    result = derive_kinf_2g_no_upscatter()
    assert result["pass"], (
        f"V_fn2.2 failed: det(M)|k=Eq 29 = {result['detM_at_eq29']} "
        f"(expected 0)"
    )


@pytest.mark.foundation
def test_v_fn2_3_phi_ratio_eq_32():
    """V_fn2.3 — phi_2/phi_1 derivation (chi-sum + balance) matches Sood Eq 32."""
    result = derive_phi_ratio_2g_no_upscatter()
    assert result["pass_chi_eliminates"], (
        "V_fn2.3: chi_1 did not eliminate from Eq 23 + Eq 24 sum"
    )
    assert result["pass"], (
        f"V_fn2.3: derived ratio = {result['ratio_derived']}, "
        f"Eq 32 = {result['ratio_eq32']}, diff = {result['diff']}"
    )


@pytest.mark.foundation
def test_v_fn_mg_1_eq_76_g2_form():
    """V_fnMG.1 — Sood Eq 76 for G=2 is the trace of the rank-1 fission matrix."""
    result = derive_kinf_mg_matrix_form()
    assert result["pass_M_is_rank1"], (
        "V_fnMG.1: A^{-1} χ νΣ_f is not rank-1 (eigenvalue check failed)"
    )
    assert result["pass_trace_equals_eq76"], (
        f"V_fnMG.1: trace(M) ≠ Eq 76; diff = {result['diff_trace']}"
    )
    assert result["pass_eigenvalue_equals_eq76"], (
        "V_fnMG.1: dominant eigenvalue ≠ Eq 76"
    )


@pytest.mark.foundation
def test_v_fn_mg_2_reduces_to_1g():
    """V_fnMG.2 — Eq 76 with G=1 reduces to Eq 19 bit-equal."""
    result = derive_kinf_mg_reduces_to_1g()
    assert result["pass"], (
        f"V_fnMG.2 failed: Eq 76 (G=1) - Eq 19 = {result['diff']}"
    )


# ═══════════════════════════════════════════════════════════════════
# Branch-2 reference-value gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_PUa_1_0_IN_matches_sood_table():
    r"""Case 1 — PUa-1-0-IN k_inf reproduces Sood Table reference to 6 digits."""
    case = PUA_1_0_IN
    k_inf = compute_kinf_1g(
        float(case.sigma_t[0]),
        float(case.sigma_s[0, 0]),
        float(case.nu_sigma_f[0]),
    )
    # Sood publishes 2.612903 (6 digits + trailing 3); a 1e-5 absolute
    # tolerance is appropriate.
    assert k_inf == pytest.approx(case.k_eff_or_kinf, abs=1e-5), (
        f"PUa-1-0-IN: computed k_inf = {k_inf}, Sood = {case.k_eff_or_kinf}, "
        f"diff = {abs(k_inf - case.k_eff_or_kinf):.3e}"
    )


@pytest.mark.foundation
def test_PU_2_0_IN_kinf_matches_sood_table():
    r"""Case 5 — PU-2-0-IN k_inf reproduces Sood Table reference to 6 digits."""
    case = PU_2_0_IN
    k_inf = compute_kinf_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_inf == pytest.approx(case.k_eff_or_kinf, abs=1e-5), (
        f"PU-2-0-IN: computed k_inf = {k_inf}, Sood = {case.k_eff_or_kinf}, "
        f"diff = {abs(k_inf - case.k_eff_or_kinf):.3e}"
    )


@pytest.mark.foundation
def test_PU_2_0_IN_flux_ratio_matches_sood_table():
    r"""Case 5 — PU-2-0-IN flux ratio reproduces Sood Eq 32 reference to 6 digits."""
    case = PU_2_0_IN
    # Sood publishes phi_2/phi_1 = phi_fast/phi_slow = 0.675229.
    # Cross-check with return_in_orpheus_order=False (Sood-style).
    sood_ratio_published = 0.675229
    fast_over_slow = compute_flux_ratio_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
        return_in_orpheus_order=False,
    )
    assert fast_over_slow == pytest.approx(sood_ratio_published, abs=1e-5), (
        f"PU-2-0-IN: computed phi_fast/phi_slow = {fast_over_slow}, "
        f"Sood = {sood_ratio_published}, "
        f"diff = {abs(fast_over_slow - sood_ratio_published):.3e}"
    )

    # And the catalogue-stored ORPHEUS-order ratio matches as well.
    slow_over_fast = compute_flux_ratio_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
        return_in_orpheus_order=True,
    )
    assert case.flux_ratio_groupwise is not None
    assert slow_over_fast == pytest.approx(
        case.flux_ratio_groupwise[1], abs=1e-5
    )


# ═══════════════════════════════════════════════════════════════════
# Branch-1 ↔ Branch-2 agreement gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_branch1_branch2_kinf_1g_agreement():
    """Branch-1 SymPy Eq 19 evaluated symbolically agrees with Branch-2 numpy."""
    import sympy as sp

    # Pull the SymPy Eq 19 result and lambdify with the PUa-1-0-IN values.
    result = derive_kinf_1g_eq_19()
    Sigma_t, Sigma_s, nu_Sigma_f = sp.symbols(
        "Sigma_t Sigma_s nu_Sigma_f", positive=True
    )
    fn = sp.lambdify(
        (Sigma_t, Sigma_s, nu_Sigma_f),
        result["k_eq19"],
        modules="numpy",
    )

    case = PUA_1_0_IN
    k_branch1 = float(fn(
        float(case.sigma_t[0]),
        float(case.sigma_s[0, 0]),
        float(case.nu_sigma_f[0]),
    ))
    k_branch2 = compute_kinf_1g(
        float(case.sigma_t[0]),
        float(case.sigma_s[0, 0]),
        float(case.nu_sigma_f[0]),
    )
    assert k_branch1 == pytest.approx(k_branch2, abs=1e-12), (
        f"Branch-1 ({k_branch1}) vs Branch-2 ({k_branch2}) disagree at 1e-12"
    )


@pytest.mark.foundation
def test_branch1_branch2_kinf_2g_agreement():
    """Branch-1 SymPy Eq 29 evaluated symbolically agrees with Branch-2 numpy
    on the PU-2-0-IN case."""
    import sympy as sp

    case = PU_2_0_IN
    # Sood-side names; map ORPHEUS arrays back.
    (
        Sigma_1, Sigma_2, Sigma_11s, Sigma_22s, Sigma_12s,
        nu_1_Sigma_1f, nu_2_Sigma_2f, chi_1, chi_2,
    ) = sp.symbols(
        "Sigma_1 Sigma_2 Sigma_11s Sigma_22s Sigma_12s "
        "nu_1_Sigma_1f nu_2_Sigma_2f chi_1 chi_2",
        positive=True,
    )
    Sigma_1rem = Sigma_1 - Sigma_11s
    Sigma_2rem = Sigma_2 - Sigma_22s
    k_eq29 = (
        chi_1 * nu_1_Sigma_1f / Sigma_1rem
        + chi_2 * (
            nu_1_Sigma_1f * Sigma_12s / (Sigma_1rem * Sigma_2rem)
            + nu_2_Sigma_2f / Sigma_2rem
        )
    )
    fn = sp.lambdify(
        (Sigma_1, Sigma_2, Sigma_11s, Sigma_22s, Sigma_12s,
         nu_1_Sigma_1f, nu_2_Sigma_2f, chi_1, chi_2),
        k_eq29,
        modules="numpy",
    )

    # Sood ↔ ORPHEUS mapping: 2 = ORPHEUS 0 (fast), 1 = ORPHEUS 1 (slow).
    k_branch1 = float(fn(
        float(case.sigma_t[1]),       # Sigma_1 (slow)
        float(case.sigma_t[0]),       # Sigma_2 (fast)
        float(case.sigma_s[1, 1]),    # Sigma_11s (slow self)
        float(case.sigma_s[0, 0]),    # Sigma_22s (fast self)
        float(case.sigma_s[0, 1]),    # Sigma_12s (from fast to slow)
        float(case.nu_sigma_f[1]),    # nu_1·Sigma_1f
        float(case.nu_sigma_f[0]),    # nu_2·Sigma_2f
        float(case.chi[1]),           # chi_1 (slow)
        float(case.chi[0]),           # chi_2 (fast)
    ))
    k_branch2 = compute_kinf_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_branch1 == pytest.approx(k_branch2, abs=1e-12), (
        f"Branch-1 ({k_branch1}) vs Branch-2 ({k_branch2}) disagree at 1e-12"
    )


@pytest.mark.foundation
def test_2g_general_equals_no_upscatter_when_sigma_21s_zero():
    """compute_kinf_2g_general (corrected Eq 28) agrees with
    compute_kinf_2g_no_upscatter (Eq 29) bit-equal when Σ_21s = 0."""
    case = PU_2_0_IN  # has sigma_s[1, 0] = 0
    assert case.sigma_s[1, 0] == 0.0, "PU-2-0-IN should have no upscatter"
    k_general = compute_kinf_2g_general(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    k_no_upscatter = compute_kinf_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_general == pytest.approx(k_no_upscatter, abs=1e-14)


# ═══════════════════════════════════════════════════════════════════
# Cross-implementation agreement gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_kinf_mg_reduces_to_kinf_1g_at_n_groups_1():
    """compute_kinf_mg with G=1 inputs reduces to compute_kinf_1g exactly."""
    case = PUA_1_0_IN
    k_mg = compute_kinf_mg(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    k_1g = compute_kinf_1g(
        float(case.sigma_t[0]),
        float(case.sigma_s[0, 0]),
        float(case.nu_sigma_f[0]),
    )
    # G=1 reduction is pure scalar algebra — should be exactly bit-equal.
    assert k_mg == k_1g, (
        f"compute_kinf_mg ({k_mg}) vs compute_kinf_1g ({k_1g}) "
        f"differ at G=1; expected bit-equal"
    )


@pytest.mark.foundation
@pytest.mark.parametrize(
    "case",
    [PUA_1_0_IN, PU_2_0_IN],
    ids=lambda c: c.case_id,
)
def test_kinf_mg_agrees_with_existing_orpheus_kinf_homogeneous(case):
    """compute_kinf_mg agrees with existing kinf_homogeneous to 1e-12.

    This is the structural-independence gate. The two solvers share
    only ``numpy`` (above the trusted-library line per
    ``algebra-of-record``). Disagreement would expose a bug in one
    or the other.
    """
    k_fn = compute_kinf_mg(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    k_orph = kinf_homogeneous(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_fn == pytest.approx(k_orph, abs=1e-12), (
        f"{case.case_id}: compute_kinf_mg ({k_fn}) vs "
        f"kinf_homogeneous ({k_orph}) disagree at 1e-12; "
        f"diff = {abs(k_fn - k_orph):.3e}"
    )


@pytest.mark.foundation
def test_PU_2_0_IN_flux_spectrum_matches_kinf_and_spectrum_homogeneous():
    """The Eq 32 flux ratio agrees with the dominant eigenvector of
    ``kinf_and_spectrum_homogeneous`` to 1e-12.

    A second structural-independence gate: Eq 32 is derived by adding
    the chi-summed balance equations, while ``kinf_and_spectrum_homogeneous``
    extracts the dominant eigenvector via :func:`numpy.linalg.eig`.
    Two completely different code paths to the same physical ratio.
    """
    case = PU_2_0_IN
    k_orph, phi_orph = kinf_and_spectrum_homogeneous(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    # phi_orph is L2-normalised non-negative; ratio is invariant under norm.
    fast_over_slow_eig = float(phi_orph[0] / phi_orph[1])
    fast_over_slow_eq32 = compute_flux_ratio_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
        return_in_orpheus_order=False,
    )
    assert fast_over_slow_eq32 == pytest.approx(fast_over_slow_eig, abs=1e-12)
