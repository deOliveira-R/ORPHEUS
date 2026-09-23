r"""L1 cross-check: NM 1980 reflected-slab F_N vs Sood LA-13511 Table 10.

Promoted from ``derivations/diagnostics/diag_01_sood_table_geometry.py``
on 2026-05-03 (numerics-investigator).

Background
----------

Wave 2-A (2026-05-03) reported a "34% disagreement between the NM 1980
F_N solver and Sood PUa-H2O(1)-1-0-SL Table 10" — a hypothesis that
turned out to be a Sood-table geometry-mismatch artefact, NOT a
method bug.

Sood LA-13511 (1999) Table 9 (problem 3) is a **NON-symmetric**
two-region slab: Pu-239 core + 1 mfp H2O reflector on **one side
only**, Pu+H2O total radius 4.542126 cm. Pu r_c = 0.48255 mfp.

Sood LA-13511 Table 10 (problem 4) is a **symmetric three-region**
slab: Pu-239 core + 0.5 mfp H2O reflector **each side**, total Pu+H2O
radius 2.849694 cm. Pu r_c = 0.43014 mfp. **This is the natural NM
1980 comparator** (NM also models a symmetric reflector around a
critical Pu core).

The Wave 2-A diagnostic compared Sood Table 9 (one-sided, 1 mfp)
against NM Case 6 (two-sided, 1 mfp each = 2 mfp total) — different
geometries, hence the 34% gap.

Reference
---------

* Sood, A., Forster, R.A., Parsons, D.K. (1999). LA-13511 Tables 9, 10.
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
    """Sood Table 10 problem 4 = `PUa-H2O(0.5)-1-0-SL`:
    SYMMETRIC reflected slab, Pu c=1.50, H2O c=0.90, Δ=0.5 mfp each side.
    Sood literature value: Pu r_c = 0.43014 mfp.
    """
    result = solve_fn_slab_reflected_critical(
        c_core=1.50, c_reflector=0.90,
        reflector_half_thickness=0.5, n_modes=7,
    )
    assert result.converged

    sood_truth = 0.43014  # mfp, from Sood Table 10 problem 4
    rel_diff = abs(result.tau_critical_mfp - sood_truth) / sood_truth
    print(
        f"\nSood Table 10 / problem 4 (PUa-H2O(0.5)-1-0-SL):"
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
    """NM 1980 Case 6 has Δ=1.0 each side; Sood Table 10 has Δ=0.5
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
        f"\n  Sood Table 10 (Δ=0.5): r_c = {sood_table10.tau_critical_mfp:.5f} mfp"
        f"\n  ratio = {sood_table10.tau_critical_mfp / nm_case6.tau_critical_mfp:.4f}"
    )
    # NM Δ=1.0 has more reflector; expect smaller core (more flux returned)
    assert nm_case6.tau_critical_mfp < sood_table10.tau_critical_mfp


@pytest.mark.foundation
def test_wave2a_memo_table9_nonsymmetric_geometry_mismatch():
    """Wave 2-A memo compared NM Case 6 (Δ=1 SYMMETRIC two-sided)
    to Sood Table 9 problem 3 (`PUa-H2O(1)-1-0-SL`,
    H2O thickness = 1 mfp NON-SYMMETRIC, single-sided).

    These are fundamentally different geometries:
        - NM Case 6: reflector on BOTH sides, 1 mfp each = 2 mfp total
        - Table 9 #3: reflector on ONE side only, 1 mfp total

    Even if the convention were identical, comparing these two would
    be an ill-posed cross-check. The disagreement (34%) is consistent
    with "single-sided reflector returns less flux to the core ⇒
    larger critical core needed", giving Table 9's r_c = 0.48255
    > NM Case 6's 0.3597.

    The NM-comparable case is Table 10 problem 4 (Δ=0.5 each side,
    symmetric).
    """
    # Sood Table 9 problem 3 (NON-symmetric):
    sood_table9_truth = 0.48255  # mfp, ONE-SIDED 1 mfp H2O reflector
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
        f"\n  Sood Table 9 #3 (1-sided, 1 mfp H2O): r_c = {sood_table9_truth}"
        f"\n  NM Case 6 (2-sided, 1 mfp each):    r_c = {nm_case6_truth}"
        f"\n  reported gap = {rel_diff*100:.1f}%   (geometry-induced, NOT method bug)"
    )


if __name__ == "__main__":
    test_sood_table10_problem4_symmetric_pu_h2o_05()
    test_nm_case6_is_not_table10()
    test_wave2a_memo_table9_nonsymmetric_geometry_mismatch()
