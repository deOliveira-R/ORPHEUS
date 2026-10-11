r"""L1 cross-check: WM-72 cylinder vs the characteristic reference's cylinder at Sood
``Ua-1-0-CY`` configuration.

This test pins the **structurally-independent** agreement between the
two cylinder critical-radius solvers in ORPHEUS:

* :func:`orpheus.derivations.continuous.singular_eigenfunction.cylinder.one_group.solve_singular_eigenfunction_cylinder_bare_critical`
  — Westfall-Metcalf 1973 singular-eigenfunction expansion via direct
  Nyström discretisation of the Mitsis cylindrical integral transport
  equation (modified Bessel kernel
  :math:`K_0(\max/\mu)\,I_0(\min/\mu)/\mu^2`).

* :class:`orpheus.derivations.continuous.characteristic.CharacteristicDerivation`
  — transport integrated along lines through the cylinder, a Galerkin
  pencil on graded panels (no Bickley-Naylor / :math:`\mathrm{Ki}_n`
  integrals; structurally distinct from the modified-Bessel kernel of
  WM-72). Until step (e2) of the characteristic-reference campaign this
  side was the trajectory-resolvent cylinder (Variant α), deleted with
  its family.

These two methods share **only** the dispersion-root primitive
(:func:`orpheus.derivations.continuous.fn_method.core.dispersion.case_nu0`,
which is a medium property — the dispersion function
:math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu) = 0` is the same
across all 1G isotropic-scattering singular-eigenfunction expansions).
Above the trusted-library line, the methods are entirely disjoint:

* WM-72: integral transport equation in `(r, t)` space with modified
  Bessel kernel.
* The characteristic reference: transport along lines, the closure as
  the boundary resolvent, a dense Galerkin pencil.

Per ``algebra-of-record`` § "Structural independence applies above
the trusted-library line", agreement at Sood ``Ua-1-0-CY`` is a true
structurally-independent L1 cross-check.

Accuracy floor — post-hardening
--------------------------------

The hardened WM-72 implementation (full Mitsis-WM Fredholm method
with Mitsis-Zweifel singular subtraction + Lagrangian derivative)
reaches **≤ 3e-7 relative** at Sood ``Ua-1-0-CY`` at
:math:`n_{\rm grid} = 24` — comparable to the 6-7 digit precision of
the published WM-72 Table II values. The cross-check now uses a
**1e-5 relative tolerance** (the target set by the brief), with
~30× margin to platform variation.

V&V triangle for Sood ``Ua-1-0-CY``:

* The characteristic reference at Sood's printed radius:
  ``tests/gates/derivations/test_characteristic_independent_references.py``
  (``test_k_is_one_at_the_one_group_cylinders_published_critical_radius``).
* WM-72 via singular-eigenfunction Fredholm: ≤ 3e-7 (this module).
* Cross-check WM-72 ↔ the characteristic reference: ≤ 1e-5 (this test;
  ``[M]`` 2026-10-10, k − 1 = −1.7e-8 at rung 3,
  ``scratch/characteristic_architecture/p1_step_e/ta_e2/probe_xrows.py``).

Two structurally-independent paths, both anchored at the published
Sood truth value to ≤ 1e-5. A third leg via ``peierls_nystrom``
(Bickley-Naylor :math:`\mathrm{Ki}_3`) is available for future
expansion.
"""
from __future__ import annotations

import pytest

from orpheus.geometry import BC, StructuredGeometry
from orpheus.derivations.continuous.singular_eigenfunction import (
    solve_singular_eigenfunction_cylinder_bare_critical,
)
from orpheus.derivations.continuous.sood_registry import SOOD2003_CASES


pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:.*roundoff error.*:scipy.integrate.IntegrationWarning"
    ),
    pytest.mark.filterwarnings("ignore::RuntimeWarning"),
]


def _characteristic_k(geometry, mixture, degree: int) -> float:
    """The characteristic reference's k for one mixture on ``geometry``, at rung ``degree`` of the joint ladder, in
    this process (:func:`~orpheus.numerics.traced_memo.bypass`)."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.derivations.continuous.characteristic import CharacteristicDerivation
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.question import Eigen
    from orpheus.numerics.traced_memo import bypass
    from orpheus.specification.specification import GeometrySpecification
    from tests.gates.derivations._characteristic_ladders import rung

    spec = GeometrySpecification(Materials({0: mixture}), geometry, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))
    with bypass():
        return float(CharacteristicDerivation(spec, rung(degree)).evaluate(Eigenvalue()).value)


@pytest.mark.l1
@pytest.mark.slow
def test_wm72_vs_the_characteristic_cylinder_at_sood_ua_1_0_cy():
    r"""L1 cross-check — WM-72 r_c agrees with the characteristic reference's r_c at the
    Sood ``Ua-1-0-CY`` benchmark configuration to ≤ 1e-5 relative.

    Both solvers reproduce the published :math:`r_c = 1.72500292` mfp
    (= 5.284935 cm) at :math:`c = 1.30` to ≤ 1e-5 relative. Their
    agreement at this precision is the structurally-independent V&V
    cross-check anchor.

    Procedure:

    1. Run the WM-72 hardened Fredholm solver to compute
       :math:`r_c^{\rm WM}` mfp at :math:`n_{\rm grid} = 24`.
    2. Convert to cm via :math:`R = r_c^{\rm WM} / \Sigma_t`.
    3. Run the characteristic reference at the WM-72-converted radius
       with vacuum walls (rung 3); the eigenvalue must be ≈ 1 to within
       :math:`10^{-5}` (the old row's 1e-3 was the trajectory resolvent's
       margin; this reference reads 1.7e-8 there).
    4. Assert WM-72 R_c agrees with Sood truth to ≤ 1e-5 (not 2%).
    """
    case = SOOD2003_CASES["Ua-1-0-CY"]
    truth_mfp = case.truth.critical_dimension_mfp  # 1.72500292
    sigma_t = float(case.materials[0].SigT[0])  # 0.32640 cm⁻¹

    # WM-72 hardened path.
    res_wm = solve_singular_eigenfunction_cylinder_bare_critical(
        c=1.30, n_grid=24, sigma_t=sigma_t,
    )
    err_wm_rel = abs(res_wm.r_c_mfp - truth_mfp) / truth_mfp
    assert err_wm_rel < 1.0e-5, (
        f"WM-72 R_c agreement with Sood truth = {err_wm_rel:.3e} > 1e-5; "
        f"R_c = {res_wm.r_c_mfp:.9f} mfp, truth = {truth_mfp}."
    )

    # The characteristic reference at WM-72's R (in cm), rung 3 (about 30 s).
    assert res_wm.r_c_cm is not None
    k_eff = _characteristic_k(StructuredGeometry.cylinder((0.0, res_wm.r_c_cm), (0,), outer=BC.vacuum), case.materials[0], 3)
    err_va = abs(k_eff - 1.0)
    assert err_va < 1.0e-5, (
        f"characteristic k_eff = {k_eff:.8f} at WM-72's R = {res_wm.r_c_cm:.6f} cm; "
        f"expected k ≈ 1 to ≤ 1e-5. Got |k - 1| = {err_va:.3e}."
    )
