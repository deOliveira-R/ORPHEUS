r"""L1 cross-check: NM 1980 reflected-slab F_N vs Sood problem 4 (2003 Table 7; 1999 Table 10).

Promoted from ``derivations/diagnostics/diag_01_sood_table_geometry.py``
on 2026-05-03 (numerics-investigator).

Background
----------

Wave 2-A (2026-05-03) reported a "34% disagreement between the NM 1980
F_N solver and Sood PUa-H2O(1)-1-0-SL Table 10" (problem 3,
PUa-H2O(1)-1-0-SL, is 1999 Table 9 and 2003 Table 6) — a hypothesis that
turned out to be a Sood-table geometry-mismatch artefact, NOT a
method bug.

Sood problem 3 (2003 Table 6; 1999 Table 9) is a **NON-symmetric**
two-region slab: Pu-239 core + 1 mfp H2O reflector on **one side
only**, Pu+H2O total radius 4.542175 cm. Pu r_c = 0.482566 mfp
(2003 edition, Table 6; the 1999 LA-13511 Table 9 printed 0.48255).

Sood problem 4 (2003 Table 7; 1999 Table 10) is a **symmetric three-region**
slab: Pu-239 core + 0.5 mfp H2O reflector **each side**, total Pu+H2O
radius 2.849725 cm. Pu r_c = 0.43015 mfp (2003; 1999 printed 0.43014).
**This is the natural NM
1980 comparator** (NM also models a symmetric reflector around a
critical Pu core).

The Wave 2-A diagnostic compared Sood problem 3 (2003 Table 6; one-sided, 1 mfp)
against NM Case 6 (two-sided, 1 mfp each = 2 mfp total) — different
geometries, hence the 34% gap.

Reference
---------

* Sood, A., Forster, R.A., Parsons, D.K. (2003). *Prog. Nucl. Energy* **42**,
  55, Tables 6 and 7 (Tables 9 and 10 of the 1999 report LA-13511).
* Neshat, K., Maiorino, J.R. (1980). *Ann. Nucl. Energy* **7**, 79-81.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_reflected_critical,
)


@pytest.mark.l1
@pytest.mark.verifies("nm1980-eq15-critical-condition")
def test_sood_table10_problem4_symmetric_pu_h2o_05():
    """Sood problem 4 (2003 Table 7; 1999 Table 10) = `PUa-H2O(0.5)-1-0-SL`:
    SYMMETRIC reflected slab, Pu c=1.50, H2O c=0.90, Δ=0.5 mfp each side.
    Sood literature value (2003 edition): Pu r_c = 0.43015 mfp.
    """
    result = solve_fn_slab_reflected_critical(
        c_core=1.50, c_reflector=0.90,
        reflector_half_thickness=0.5, n_modes=7,
    )
    assert result.converged

    sood_truth = 0.43015  # mfp, Sood 2003 Table 7, problem 4
    rel_diff = abs(result.tau_critical_mfp - sood_truth) / sood_truth
    print(
        f"\nSood 2003 Table 7 / problem 4 (PUa-H2O(0.5)-1-0-SL):"
        f"\n  Sood r_c = {sood_truth:.5f} mfp"
        f"\n  F_7 r_c  = {result.tau_critical_mfp:.5f} mfp"
        f"\n  rel diff = {rel_diff:.4e}"
    )
    # Sood publishes 5 sig figs; tolerance: 1e-3 absolute.
    assert abs(result.tau_critical_mfp - sood_truth) < 1e-3, (
        f"F_7 = {result.tau_critical_mfp:.5f}, "
        f"Sood = {sood_truth:.5f}, "
        f"diff = {abs(result.tau_critical_mfp - sood_truth):.4e}"
    )


@pytest.mark.foundation
def test_nm_case6_is_not_table10():
    """NM 1980 Case 6 has Δ=1.0 each side; Sood problem 4 (2003 Table 7) has Δ=0.5
    each side. They are NOT the same problem."""
    nm_case6 = solve_fn_slab_reflected_critical(
        c_core=1.50, c_reflector=0.90,
        reflector_half_thickness=1.0, n_modes=7,
    )
    sood_table10 = solve_fn_slab_reflected_critical(
        c_core=1.50, c_reflector=0.90,
        reflector_half_thickness=0.5, n_modes=7,
    )
    print(
        f"\n  NM Case 6 (Δ=1.0): r_c = {nm_case6.tau_critical_mfp:.5f} mfp"
        f"\n  Sood 2003 Table 7 (Δ=0.5): r_c = {sood_table10.tau_critical_mfp:.5f} mfp"
        f"\n  ratio = {sood_table10.tau_critical_mfp / nm_case6.tau_critical_mfp:.4f}"
    )
    # NM Δ=1.0 has more reflector; expect smaller core (more flux returned)
    assert nm_case6.tau_critical_mfp < sood_table10.tau_critical_mfp


@pytest.mark.foundation
def test_wave2a_memo_table9_nonsymmetric_geometry_mismatch():
    """Wave 2-A memo compared NM Case 6 (Δ=1 SYMMETRIC two-sided)
    to Sood 2003 Table 6 problem 3 (`PUa-H2O(1)-1-0-SL`,
    H2O thickness = 1 mfp NON-SYMMETRIC, single-sided).

    These are fundamentally different geometries:
        - NM Case 6: reflector on BOTH sides, 1 mfp each = 2 mfp total
        - 2003 Table 6 #3: reflector on ONE side only, 1 mfp total

    Even if the convention were identical, comparing these two would
    be an ill-posed cross-check. The disagreement (34%) is consistent
    with "single-sided reflector returns less flux to the core ⇒
    larger critical core needed", giving problem 3's r_c = 0.482566
    > NM Case 6's 0.3597.

    The NM-comparable case is 2003 Table 7 problem 4 (Δ=0.5 each side,
    symmetric).
    """
    # Sood 2003 Table 6 problem 3 (NON-symmetric):
    sood_table9_truth = 0.482566  # mfp, ONE-SIDED 1 mfp H2O reflector (2003 Table 6)
    # NM Case 6 (SYMMETRIC, Δ=1 each side):
    nm_case6_truth = 0.3597  # mfp

    rel_diff = abs(sood_table9_truth - nm_case6_truth) / nm_case6_truth
    # The Wave 2-A "34% gap" — but it's a geometry mismatch, not a method bug
    assert rel_diff > 0.30, (
        "Wave 2-A's reported 34% gap reproduced — but the geometries "
        "are incompatible (one-sided vs two-sided reflector)."
    )
    print(
        f"\n  Wave 2-A geometry-mismatch verification:"
        f"\n  Sood 2003 Table 6 #3 (1-sided, 1 mfp H2O): r_c = {sood_table9_truth}"
        f"\n  NM Case 6 (2-sided, 1 mfp each):    r_c = {nm_case6_truth}"
        f"\n  reported gap = {rel_diff*100:.1f}%   (geometry-induced, NOT method bug)"
    )


if __name__ == "__main__":
    test_sood_table10_problem4_symmetric_pu_h2o_05()
    test_nm_case6_is_not_table10()
    test_wave2a_memo_table9_nonsymmetric_geometry_mismatch()


@pytest.mark.l1
@pytest.mark.parametrize("n_modes", [11, 13, 15])
def test_sood_problem4_at_the_published_digit(n_modes):
    """Sood problem 4 at the precision it is published to, which decides the
    edition (P1 step 2b; the user ruled the 2003 edition the reference).

    The 2003 edition prints 0.43015 mfp and the 1999 edition 0.43014; a
    five-digit value stands for the interval of half a unit in its last
    place, 5e-6. F_N converged in N (measured 2026-09-29: 0.4301459,
    0.4301477 and 0.4301459 at N = 11, 13 and 15, an oscillation of
    1.9e-6) lies inside 2003's interval and outside 1999's at every one of
    these N. The 1e-3 gates above cannot tell the two editions apart.
    """
    tau = solve_fn_slab_reflected_critical(
        c_core=1.50, c_reflector=0.90, reflector_half_thickness=0.5, n_modes=n_modes,
    ).tau_critical_mfp
    half_unit = 0.5e-5
    assert abs(tau - 0.43015) <= half_unit, f"N={n_modes}: tau = {tau!r} outside 0.43015 +/- 5e-6"
    assert abs(tau - 0.43014) > half_unit, f"N={n_modes}: tau = {tau!r} cannot tell 1999 from 2003"
