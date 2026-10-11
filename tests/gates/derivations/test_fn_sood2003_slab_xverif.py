r"""L1 cross-check: F_N slab vs the characteristic reference's slab at Sood ``Ua-1-0-SL``.

This is the load-bearing **structural-independence** verification
between two methods that share only ``numpy``/``scipy`` (above the
trusted-library line):

* **F_N method** (this slice's Branch-2 solver
  :func:`solve_fn_slab_bare_critical`) — Case singular eigenfunctions
  + Wiener-Hopf factorization (slab F_N collocation per Siewert-
  Benoist Part I + Grandjean-Siewert Part II).
* **The characteristic reference**
  (:class:`~orpheus.derivations.continuous.characteristic.CharacteristicDerivation`,
  vacuum walls) — transport integrated along lines, a Galerkin pencil
  on graded panels. Until step (e2) of the characteristic-reference
  campaign the second method was the trajectory-resolvent slab
  (Variant α), deleted with its family.

Both methods produce the same Sood/KLL truth :math:`a_c =
0.93772556` mfp via mathematically disjoint paths; their agreement is
the strongest available L1 evidence for either solver.

The agreement test uses the Sood Ua-1-0-SL XS at the published
critical half-thickness :math:`a = 2.872934` cm (so the slab full
width is :math:`L = 5.745868` cm), and verifies that BOTH methods
return :math:`k_{\rm eff} \approx 1.0` to within their respective
quadrature floors.

Quadrature notes
----------------

* F_N at :math:`N = 10` reaches :math:`a_c` truth to ~5e-6.
* The characteristic reference at rung 4 reads :math:`k = 1 - 3.5\times 10^{-6}`
  at F_N's thickness, converged (rungs 3 and 4 agree to 1.2e-8,
  ``[M]`` 2026-10-10, ``scratch/characteristic_architecture/p1_step_e/ta_e2/probe_xrows.py``):
  the residual is F_N's own ~5e-6 in :math:`a_c`.

Cross-check tolerance is set to 5e-5 to accommodate both floors.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.sood_registry import UA_1_0_SL_STUB
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)
from orpheus.geometry import BC, StructuredGeometry
from tests.gates.derivations.test_characteristic_system import _mixture


# Suppress the bracket-scan divide-by-zero warnings (see slab tests).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
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
def test_fn_slab_vs_the_characteristic_slab_at_sood_ua_1_0_sl():
    r"""Both methods return :math:`k_{\rm eff} = 1` at the Sood Ua-1-0-SL
    truth dimension :math:`L = 5.745868` cm.

    F_N runs at :math:`N = 10` (its inherent convergence floor for this
    case is ~5e-6 in :math:`a_c`); the characteristic reference runs at
    rung 4 (converged to 1.2e-8 here). The agreement is required to
    ≤ 5e-5, as before the re-point.
    """
    case = UA_1_0_SL_STUB
    sigma_t = float(case.materials[0].SigT[0])
    sigma_s = float(case.materials[0].SigS[0][0, 0])
    nu_sigma_f = float(case.materials[0].SigP[0])

    # F_N: returns a_c in mfp; convert to cm via /sigma_t.
    res_fn = solve_fn_slab_bare_critical(c=(sigma_s + nu_sigma_f) / sigma_t,
                                          n_modes=10)
    a_fn_mfp = res_fn.a_critical_mfp
    a_fn_cm = a_fn_mfp / sigma_t  # 2.872xx cm
    L_full_cm = 2.0 * a_fn_cm

    # F_N self-consistency: should match Sood truth.
    truth_mfp = case.truth.critical_dimension_mfp
    err_fn = abs(a_fn_mfp - truth_mfp)
    assert err_fn < 5e-5, (
        f"F_N internal: a_c={a_fn_mfp:.10f} vs truth {truth_mfp}, "
        f"err={err_fn:.2e}"
    )

    # The characteristic reference at the F_N-truth thickness should give k_eff = 1.
    geometry = StructuredGeometry.slab((0.0, L_full_cm), (0,), left=BC.vacuum, right=BC.vacuum)
    k_eff_at_truth = _characteristic_k(geometry, case.materials[0], 4)
    err_va = abs(k_eff_at_truth - 1.0)
    assert err_va < 5e-5, (
        f"characteristic reference at F_N truth thickness: k_eff={k_eff_at_truth:.8f}, "
        f"err={err_va:.2e}"
    )

    # The two-method cross-check: both report critical-thickness
    # agreement at the level of their respective floors. If one method
    # is biased by O(1e-3) the assertion above will fail and the bug
    # source can be pinned to whichever method moved.


@pytest.mark.l1
def test_fn_slab_vs_the_characteristic_slab_at_grandjean_siewert_c150():
    r"""Cross-check at a different :math:`c` value — Grandjean-Siewert
    Table XI :math:`c = 1.50` row. Both methods should agree at the
    same critical thickness to ≤ 1e-4 (looser at higher :math:`c` due
    to thinner slab + steeper boundary gradient).
    """
    c = 1.50
    # Solve via F_N (unit XS).
    res_fn = solve_fn_slab_bare_critical(c=c, n_modes=10)
    a_fn = res_fn.a_critical_mfp  # in mfp = cm at sigma_t=1.

    # The characteristic reference at L = 2*a_fn using unit XS (sigma_s = 0, nu sigma_f = c).
    geometry = StructuredGeometry.slab((0.0, 2.0 * a_fn), (0,), left=BC.vacuum, right=BC.vacuum)
    k_eff = _characteristic_k(geometry, _mixture([1.0], [[0.0]], None, [c], [1.0]), 4)
    err = abs(k_eff - 1.0)
    assert err < 1e-4, (
        f"GS c=1.50: F_N a={a_fn:.10f}, characteristic k_eff at that "
        f"thickness = {k_eff:.8f}, err={err:.2e}"
    )
