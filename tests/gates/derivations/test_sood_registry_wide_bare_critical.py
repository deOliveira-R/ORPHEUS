r"""Phase B3 wide-enumeration 1G bare-critical gate — Sood/Forster/Parsons LA-13511.

This file covers the **5 1G bare-critical cases** in
:data:`orpheus.derivations.continuous.sood_registry.WIDE_SLICE_BARE_CRITICAL_1G`
that the existing slab/sphere F_N solvers can solve TODAY:

* **Slab (3 cases)**: PUa-1-0-SL (c=1.50), PUb-1-0-SL (c=1.40),
  UD2O-1-0-SL (c=1.02).
* **Sphere (2 cases)**: PUb-1-0-SP (c=1.40), UD2O-1-0-SP (c=1.02).

The Phase A first-slice case ``Ua-1-0-SL`` (c=1.30) and ``Ua-1-0-SP``
(c=1.30) are NOT re-tested here — they live in
``test_fn_la13511_slab.py`` / ``test_fn_la13511_sphere.py``.

STUB cases (cylinders, 2G bare-critical) are tracked in the registry
but explicitly skipped here with ``pytest.skip`` markers — see the
``test_*_stub_explicitly_skipped`` cases below.

Tag policy: ``@pytest.mark.foundation`` for the convergence/registry
bookkeeping invariants, ``@pytest.mark.l1`` for the F_N reference-value
gates (these implement an analytical-reference-pillar L1 verification
per ``vv-principles``).

Tolerance: ≤ 1e-5 absolute on the critical dimension in mfp. Sood
publishes 6-10 digits; F_N at the chosen N reaches at most 3e-6
(slab N=12) or 5e-8 (sphere N=10) — both well within budget.

References
----------

* Sood, Forster, Parsons (1999), LA-13511, Tables 6, 7, 17.
* Phase B3 closeout memo:
  ``.claude/agent-memory/method-implementer/sood_registry_wide_enumeration_phase_b3.md``.
"""
from __future__ import annotations

import pytest

# Suppress numpy divide-by-zero warnings from the F_N bracket scan.
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]

from orpheus.derivations.continuous.fn_method.slab.one_group import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.fn_method.sphere.one_group import (
    solve_fn_sphere_bare_critical,
)
from orpheus.derivations.continuous.sood_registry import (
    LA13511_CASES,
    WIDE_SLICE_BARE_CRITICAL_1G,
    WIDE_SLICE_STUBS,
)


# ═══════════════════════════════════════════════════════════════════
# Group cases by geometry for parametrize fan-out
# ═══════════════════════════════════════════════════════════════════


SLAB_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_BARE_CRITICAL_1G
    if case.geometry_kind == "slab"
]

SPHERE_CASE_IDS = [
    case.case_id for case in WIDE_SLICE_BARE_CRITICAL_1G
    if case.geometry_kind == "sphere"
]


def _compute_c(case) -> float:
    """Compute the secondaries ratio c = (Σ_s + νΣ_f)/Σ_t for a 1G case."""
    sigma_t = float(case.materials[0].SigT[0])
    sigma_s = float(case.materials[0].SigS[0][0, 0])
    nu_sigma_f = float(case.materials[0].SigP[0])
    return (sigma_s + nu_sigma_f) / sigma_t


# ═══════════════════════════════════════════════════════════════════
# Slab F_N L1 reference-value gate
# ═══════════════════════════════════════════════════════════════════
#
# N=12 chosen empirically: hits ≤ 3e-6 absolute for all three cases.
# Lower N (8/10) gives 1.6e-5/4.8e-6 — borderline; higher N (14)
# loses bracket for low c (UD2O at c=1.02).


SLAB_N_MODES = 12


@pytest.mark.l1
@pytest.mark.parametrize("case_id", SLAB_CASE_IDS)
def test_slab_fn_critical_dimension_matches_sood(case_id: str) -> None:
    r"""1G bare slab F_N reproduces Sood :math:`a_c` to ≤ 1e-5 absolute.

    Solver: :func:`solve_fn_slab_bare_critical` at N=12 modes.
    Reference: LA-13511 Tables 6, 7, 17 (slab cases 2, 6, 22).
    """
    case = LA13511_CASES[case_id]
    truth = case.truth.critical_dimension_mfp
    c = _compute_c(case)
    res = solve_fn_slab_bare_critical(c=c, n_modes=SLAB_N_MODES)
    err = abs(res.a_critical_mfp - truth)
    assert err < 1e-5, (
        f"{case_id}: F_{SLAB_N_MODES} a_c = {res.a_critical_mfp}, "
        f"Sood truth = {truth}, err = {err:.3e}"
    )
    # Also pin determinant residual sanity (matrix solved at the root).
    assert abs(res.determinant_residual) < 1e-10, (
        f"{case_id}: det M({res.a_critical_mfp}) residual = "
        f"{res.determinant_residual} not converged"
    )


# ═══════════════════════════════════════════════════════════════════
# Sphere F_N L1 reference-value gate
# ═══════════════════════════════════════════════════════════════════
#
# N=10 chosen empirically: hits ≤ 5e-8 absolute (way better than 1e-5).
# Spheres converge much faster than slabs in F_N.


SPHERE_N_MODES = 10


@pytest.mark.l1
@pytest.mark.parametrize("case_id", SPHERE_CASE_IDS)
def test_sphere_fn_critical_dimension_matches_sood(case_id: str) -> None:
    r"""1G bare sphere F_N reproduces Sood :math:`R_c` to ≤ 1e-5 absolute.

    Solver: :func:`solve_fn_sphere_bare_critical` at N=10 modes.
    Reference: LA-13511 Tables 7, 17 (sphere cases 8, 24).
    """
    case = LA13511_CASES[case_id]
    truth = case.truth.critical_dimension_mfp
    c = _compute_c(case)
    res = solve_fn_sphere_bare_critical(c=c, n_modes=SPHERE_N_MODES)
    err = abs(res.R_critical_mfp - truth)
    assert err < 1e-5, (
        f"{case_id}: F_{SPHERE_N_MODES} R_c = {res.R_critical_mfp}, "
        f"Sood truth = {truth}, err = {err:.3e}"
    )


# ═══════════════════════════════════════════════════════════════════
# Stub coverage tracker — cylinder + 2G bare-critical
# ═══════════════════════════════════════════════════════════════════
#
# These cases live in WIDE_SLICE_STUBS — registered but no solver
# tests until the dispatched solver lands. Marked as skip so pytest
# reports them rather than silently passing.


CYLINDER_STUB_IDS = [
    case.case_id for case in WIDE_SLICE_STUBS
    if case.geometry_kind == "cylinder"
]

TWO_G_BARE_CRITICAL_STUB_IDS = [
    case.case_id for case in WIDE_SLICE_STUBS
    if case.materials[0].ng == 2 and case.geometry_kind in ("slab", "sphere")
]


@pytest.mark.l1
@pytest.mark.parametrize("case_id", CYLINDER_STUB_IDS)
def test_cylinder_stub_solver_pending_b1_dispatch(case_id: str) -> None:
    """STUB: B1 dispatch will activate this case once the cylinder F_N solver lands."""
    case = LA13511_CASES[case_id]
    assert case.geometry_kind == "cylinder"
    assert case.truth.critical_dimension_mfp is not None  # truth registered
    pytest.skip(
        f"{case_id}: cylinder F_N solver pending — B1 dispatch will "
        f"activate (Westfall-Metcalf 1973 NSE 52, 1)"
    )


@pytest.mark.l1
@pytest.mark.parametrize("case_id", TWO_G_BARE_CRITICAL_STUB_IDS)
def test_2g_bare_critical_stub_solver_pending_2g_fn(case_id: str) -> None:
    """STUB: needs Siewert-Thomas 1986 2G F_N machinery (matrix dispersion law)."""
    case = LA13511_CASES[case_id]
    assert case.materials[0].ng == 2
    assert case.geometry_kind in ("slab", "sphere")
    assert case.truth.critical_dimension_mfp is not None  # truth registered
    pytest.skip(
        f"{case_id}: 2G bare-critical F_N solver pending — needs "
        f"Siewert-Thomas 1986 2G machinery (matrix Λ dispersion law, "
        f"complex 2x2 Case eigenvalues). See "
        f".claude/agent-memory/method-implementer/"
        f"fn_method_sphere_fn_proper.md for the geometry-sign abstraction."
    )


# ═══════════════════════════════════════════════════════════════════
# Bookkeeping
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_wide_slice_bare_critical_count() -> None:
    """Phase B3 wide-slice 1G bare-critical has the expected case count (5)."""
    assert len(WIDE_SLICE_BARE_CRITICAL_1G) == 5, (
        f"Expected 5 wide-slice 1G bare-critical cases, "
        f"got {len(WIDE_SLICE_BARE_CRITICAL_1G)}"
    )


@pytest.mark.foundation
def test_wide_slice_stubs_count() -> None:
    """Phase B3 wide-slice STUBS has the expected case count (12)."""
    # 2 cylinder stubs (PUb, UD2O) + 5 PU/U/UAL/URRa/UD2O × {SL, SP} = 12
    assert len(WIDE_SLICE_STUBS) == 12, (
        f"Expected 12 wide-slice stubs, got {len(WIDE_SLICE_STUBS)}"
    )
