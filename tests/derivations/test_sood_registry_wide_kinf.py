r"""Phase B3 wide-enumeration k_inf gate — Sood/Forster/Parsons LA-13511.

This file covers the **20 k_inf cases** in
:data:`orpheus.derivations.continuous.sood_registry.WIDE_SLICE_KINF`.
Each case carries published reference values from LA-13511 (Tables 5,
12, 16, 20, 24, 28, 30-31, 33-34, 36-37, 39-40, 43-44, 46-47, 49-50,
59-67); the matching :func:`compute_kinf_*` Branch-2 function reproduces
them at machine precision, modulo Sood's published 5-7 digit truth.

Tag policy: all tests are :func:`pytest.mark.foundation` because the
underlying claims are pure linear-algebra invariants of the multi-group
infinite-medium k_inf machinery — there is no spatial/angular
discretisation, hence no L0/L1/L2 ladder applies.

Cross-implementation gate: each case is computed via TWO structurally
independent paths (Sood Eq 19/29/76 ↔ ``kinf_homogeneous``'s
:func:`numpy.linalg.eig` of :math:`A^{-1}F`). They share only ``numpy``
above the trusted-library line per ``algebra-of-record``.

References
----------

* Sood, Forster, Parsons (1999), LA-13511.
* Phase B3 closeout memo:
  ``.claude/agent-memory/method-implementer/sood_registry_wide_enumeration_phase_b3.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_spectrum_homogeneous,
    kinf_homogeneous,
)
from orpheus.derivations.continuous.fn_method.multi_group import (
    compute_kinf_1g,
    compute_kinf_2g_general,
    compute_kinf_2g_no_upscatter,
    compute_kinf_mg,
)
from orpheus.derivations.continuous.sood_registry import (
    LA13511_CASES,
    WIDE_SLICE_KINF,
)


# ═══════════════════════════════════════════════════════════════════
# Group cases by ngroups and upscatter for parametrize fan-out
# ═══════════════════════════════════════════════════════════════════


KINF_1G_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_KINF if case.n_groups == 1
]

KINF_2G_NO_UPSCATTER_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_KINF
    if case.n_groups == 2 and float(case.sigma_s[1, 0]) == 0.0
]

KINF_2G_WITH_UPSCATTER_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_KINF
    if case.n_groups == 2 and float(case.sigma_s[1, 0]) > 0.0
]

KINF_MG_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_KINF if case.n_groups >= 3
]

ALL_KINF_CASE_IDS = [case.case_id for case in WIDE_SLICE_KINF]


# ═══════════════════════════════════════════════════════════════════
# 1G k_inf gate (11 cases including 4 P_1-anisotropic which reduce to
# isotropic for k_inf since P_1 doesn't affect infinite-medium balance)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", KINF_1G_CASE_IDS)
def test_kinf_1g_matches_sood(case_id: str) -> None:
    r"""1G k_inf via Sood Eq 19 reproduces published table value to ≤ 1e-5."""
    case = LA13511_CASES[case_id]
    truth = case.truth.k_eff_or_kinf
    k = compute_kinf_1g(
        float(case.sigma_t[0]),
        float(case.sigma_s[0, 0]),
        float(case.nu_sigma_f[0]),
    )
    assert k == pytest.approx(truth, abs=1e-5), (
        f"{case_id}: Eq 19 k_inf = {k}, Sood truth = {truth}, "
        f"diff = {abs(k - truth):.3e}"
    )


# ═══════════════════════════════════════════════════════════════════
# 2G k_inf gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", KINF_2G_NO_UPSCATTER_CASE_IDS)
def test_kinf_2g_no_upscatter_matches_sood(case_id: str) -> None:
    r"""2G k_inf (no upscatter) via Sood Eq 29 reproduces published table value to ≤ 1e-5."""
    case = LA13511_CASES[case_id]
    truth = case.truth.k_eff_or_kinf
    k = compute_kinf_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k == pytest.approx(truth, abs=1e-5), (
        f"{case_id}: Eq 29 k_inf = {k}, Sood truth = {truth}, "
        f"diff = {abs(k - truth):.3e}"
    )


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", KINF_2G_WITH_UPSCATTER_CASE_IDS)
def test_kinf_2g_with_upscatter_matches_sood(case_id: str) -> None:
    r"""2G k_inf (with upscatter) via Sood Eq 28 (general) reproduces published to ≤ 1e-5.

    The no-upscatter Eq-29 specialisation does NOT apply; this gate
    pins :func:`compute_kinf_2g_general` for the URRb/URRc cases that
    have :math:`\Sigma_{21s} > 0`.
    """
    case = LA13511_CASES[case_id]
    truth = case.truth.k_eff_or_kinf
    assert float(case.sigma_s[1, 0]) > 0.0, (
        f"{case_id}: this fixture is for upscatter cases only"
    )
    k = compute_kinf_2g_general(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k == pytest.approx(truth, abs=1e-5), (
        f"{case_id}: Eq 28 (general) k_inf = {k}, Sood truth = {truth}, "
        f"diff = {abs(k - truth):.3e}"
    )


# ═══════════════════════════════════════════════════════════════════
# 3G/6G k_inf gate (compute_kinf_mg via Eq 76)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", KINF_MG_CASE_IDS)
def test_kinf_mg_matches_sood(case_id: str) -> None:
    r"""Multi-group k_inf via Sood Eq 76 reproduces published table value to ≤ 1e-5."""
    case = LA13511_CASES[case_id]
    truth = case.truth.k_eff_or_kinf
    k = compute_kinf_mg(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k == pytest.approx(truth, abs=1e-5), (
        f"{case_id}: Eq 76 k_inf = {k}, Sood truth = {truth}, "
        f"diff = {abs(k - truth):.3e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Cross-implementation structural-independence gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", ALL_KINF_CASE_IDS)
def test_kinf_mg_agrees_with_kinf_homogeneous_eig(case_id: str) -> None:
    r"""F_N ``compute_kinf_mg`` agrees with ``kinf_homogeneous`` to 1e-12.

    Two structurally independent paths to the same physical k_inf:

    * ``compute_kinf_mg`` evaluates Sood Eq 76 directly via
      ``numpy.linalg.solve`` (closed-form result).
    * ``kinf_homogeneous`` builds the transfer matrix
      :math:`A^{-1}\,\chi\,(\nu\Sigma_f)^T` and extracts the dominant
      eigenvalue via ``numpy.linalg.eig``.

    They share only ``numpy`` above the trusted-library line; a
    disagreement would expose a bug in one or the other (algebra-of-record
    cross-check pillar).
    """
    case = LA13511_CASES[case_id]
    k_fn = compute_kinf_mg(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    k_orph = kinf_homogeneous(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_fn == pytest.approx(k_orph, abs=1e-12), (
        f"{case_id}: compute_kinf_mg ({k_fn}) vs "
        f"kinf_homogeneous ({k_orph}) disagree at 1e-12; "
        f"diff = {abs(k_fn - k_orph):.3e}"
    )


# ═══════════════════════════════════════════════════════════════════
# 2G no-upscatter ↔ general agreement gate
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case_id", KINF_2G_NO_UPSCATTER_CASE_IDS)
def test_kinf_2g_general_equals_no_upscatter_when_sigma_21s_zero(case_id: str) -> None:
    r"""For :math:`\Sigma_{21s} = 0`, the general Eq-28 result equals the
    no-upscatter Eq-29 result bit-for-bit.

    This is a sanity check that the corrected Eq-28 algebra reduces
    cleanly to Eq-29 in the no-upscatter limit (asserted abstractly
    in :func:`derive_kinf_2g_general_from_matrix`; this test pins it
    on every concrete no-upscatter 2G Sood case).
    """
    case = LA13511_CASES[case_id]
    assert float(case.sigma_s[1, 0]) == 0.0
    k_general = compute_kinf_2g_general(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    k_no_upscatter = compute_kinf_2g_no_upscatter(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    assert k_general == pytest.approx(k_no_upscatter, abs=1e-14), (
        f"{case_id}: general ({k_general}) ≠ no_upscatter ({k_no_upscatter})"
    )


# ═══════════════════════════════════════════════════════════════════
# Multi-group flux-spectrum gate (URR-3-0-IN, URR-6-0-IN)
# ═══════════════════════════════════════════════════════════════════
#
# These two cases are the only Phase-B3 entries that publish multi-group
# flux ratios at the case spec (LA-13511 Eqs 64-65, Section IX). Verifies
# that the dominant-eigenvector spectrum from kinf_and_spectrum_homogeneous
# matches Sood's published ratios.


@pytest.mark.foundation
def test_URR_3_0_IN_flux_spectrum() -> None:
    r"""URR-3-0-IN: dominant-eigenvector spectrum reproduces Sood Eq 64-65 ratios.

    Sood publishes :math:`\phi_2/\phi_3 = 0.480` (mid/fast) and
    :math:`\phi_1/\phi_3 = 0.150` (slow/fast), both exact rational
    constructions. ORPHEUS convention: g=0 fast, g=1 mid, g=2 slow.
    """
    case = LA13511_CASES["URR-3-0-IN"]
    _, phi = kinf_and_spectrum_homogeneous(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    phi_normalized = phi / phi[0]   # normalise to fast = 1
    # phi_normalized[0] = 1.0 (fast, ORPHEUS 0 = Sood g3)
    # phi_normalized[1] = phi_mid/phi_fast = 0.480 (Sood g2/g3)
    # phi_normalized[2] = phi_slow/phi_fast = 0.150 (Sood g1/g3)
    assert phi_normalized[0] == pytest.approx(1.0, abs=1e-12)
    assert phi_normalized[1] == pytest.approx(0.480, abs=1e-6)
    assert phi_normalized[2] == pytest.approx(0.150, abs=1e-6)


@pytest.mark.foundation
def test_URR_6_0_IN_flux_spectrum_mirror() -> None:
    r"""URR-6-0-IN: 6-group spectrum has the published mirror structure.

    Sood designed the 6-group case as 2 coupled URR-3-0-IN blocks with
    χ distributing fission across both. The expected spectrum:

    * φ(g=6 fast) = φ(g=1 fast) — mirror pair
    * φ(g=5 mid)  = φ(g=2 mid)  — mirror pair
    * φ(g=4 slow) = φ(g=3 slow) — mirror pair
    * Within each block: φ_mid/φ_fast = 0.480, φ_slow/φ_fast = 0.150.

    ORPHEUS index convention: 0=g6, 1=g5, 2=g4, 3=g3, 4=g2, 5=g1.
    """
    case = LA13511_CASES["URR-6-0-IN"]
    _, phi = kinf_and_spectrum_homogeneous(
        case.sigma_t, case.sigma_s, case.nu_sigma_f, case.chi,
    )
    # Mirror property: phi[i] == phi[5-i]
    for i in range(3):
        assert phi[i] == pytest.approx(phi[5 - i], abs=1e-12), (
            f"URR-6-0-IN mirror broken: phi[{i}] = {phi[i]} vs phi[{5-i}] = {phi[5-i]}"
        )
    # Ratios within the upper block (ORPHEUS 0,1,2 = Sood 6,5,4):
    phi_norm_upper = phi[:3] / phi[0]
    assert phi_norm_upper[1] == pytest.approx(0.480, abs=1e-6)
    assert phi_norm_upper[2] == pytest.approx(0.150, abs=1e-6)


# ═══════════════════════════════════════════════════════════════════
# Sanity: no test-case duplication
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_wide_slice_kinf_count() -> None:
    """Phase B3 wide-slice k_inf has the expected case count (20)."""
    assert len(WIDE_SLICE_KINF) == 20, (
        f"Expected 20 wide-slice k_inf cases, got {len(WIDE_SLICE_KINF)}"
    )


@pytest.mark.foundation
def test_wide_slice_kinf_partition_complete() -> None:
    """Every WIDE_SLICE_KINF case is in exactly one of the parametrize buckets."""
    bucketed = (
        set(KINF_1G_CASE_IDS)
        | set(KINF_2G_NO_UPSCATTER_CASE_IDS)
        | set(KINF_2G_WITH_UPSCATTER_CASE_IDS)
        | set(KINF_MG_CASE_IDS)
    )
    expected = {case.case_id for case in WIDE_SLICE_KINF}
    assert bucketed == expected, (
        f"Bucket partition mismatch:\n"
        f"  in buckets but not WIDE_SLICE_KINF: {bucketed - expected}\n"
        f"  in WIDE_SLICE_KINF but no bucket: {expected - bucketed}"
    )
    # Disjointness:
    buckets = (
        set(KINF_1G_CASE_IDS),
        set(KINF_2G_NO_UPSCATTER_CASE_IDS),
        set(KINF_2G_WITH_UPSCATTER_CASE_IDS),
        set(KINF_MG_CASE_IDS),
    )
    for i, bi in enumerate(buckets):
        for j, bj in enumerate(buckets):
            if i < j:
                assert bi.isdisjoint(bj), (
                    f"Buckets {i} and {j} overlap on {bi & bj}"
                )
