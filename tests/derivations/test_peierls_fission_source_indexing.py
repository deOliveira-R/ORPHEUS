r"""L1 — Peierls multi-group fission source-indexing gate (ERR-063).

INTRINSIC PROPERTY (the structural teeth that distinguish source- from
sink-indexed :math:`\chi`): *a non-fissile region's fission emission
spectrum* :math:`\chi` *must not affect* :math:`k_{\rm eff}`. A non-fissile
region (:math:`\nu\Sigma_f = 0`) emits zero fission neutrons, so its
:math:`\chi` — the energy distribution OF the fission neutrons it would
emit — is physically meaningless and must drop out of the eigenvalue
problem entirely. This holds for the FORWARD problem regardless of
geometry, group count, or boundary condition.

Physics ground (structurally independent of the ORPHEUS theory page, which
documented the wrong convention — reference contamination, vv §6): Hébert
(2009) *Applied Reactor Physics* Eq. 3.57/3.58 — the fission emission
density :math:`\chi_g(r)\,\sum_{g'}\nu\Sigma_{f,g'}(r)\,\varphi(r)` is a
single LOCAL quantity at the fission point :math:`r`; when
:math:`\nu\Sigma_f(r)=0` the whole emission is zero and :math:`\chi(r)`
multiplies nothing. The transport kernel :math:`K(r_i,r_j)` is the sole
carrier of spatial coupling and never touches :math:`\chi`.

This is the promotion target of
``derivations/diagnostics/diag_err063_probe_e_intrinsic_property.py``.
It is RED under the bug (sink-indexed :math:`\chi_i` leaks the non-fissile
region's spectrum into the fission birth at fissile-node→non-fissile-node
kernel couplings) and GREEN under the fix (source-indexed :math:`\chi_j`).
Mutating ``geometry.py`` / ``slab.py`` back to ``chi[i]`` reddens it — the
``catches("ERR-063")`` teeth, mutation-verified.

vv #11 both-legs structure:
  - POSITIVE leg (``test_fissile_region_chi_genuinely_affects_keff``):
    the FISSILE region's :math:`\chi` DOES move :math:`k_{\rm eff}` (the
    property is non-degenerate — guards against a vacuous "χ never matters"
    reading of the negative leg).
  - NEGATIVE leg / teeth
    (``test_nonfissile_region_chi_does_not_affect_keff``): mutating ONLY a
    non-fissile region's :math:`\chi` must NOT move :math:`k_{\rm eff}`.

Mode-8: the suite runs ``python -O`` (bare ``assert`` is stripped); every
assertion routes through ``np.testing.assert_array_less`` (a function call
that fires under ``-O``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_nystrom.geometry import (
    CYLINDER_1D,
    SPHERE_1D,
    solve_peierls_mg,
)

# k_eff must be invariant to a non-fissile region's χ to solver tolerance,
# with generous headroom. The bug produces O(1e-2) drift (4.7% measured).
_INVARIANCE_TOL = 1e-7

pytestmark = [
    pytest.mark.l1,
    pytest.mark.catches("ERR-063"),
    pytest.mark.verifies("peierls-mg-operator"),
]


def _solve_two_region(geometry, *, chi: np.ndarray) -> float:
    r"""Solve the fixed 2-region test problem for a given per-region χ.

    Inner region (0): fissile, distinct birth spectrum. Outer region (1):
    pure absorber/scatterer (:math:`\nu\Sigma_f = 0`). 2-group, downscatter
    g0->g1. Only ``chi`` (the ``(n_regions, ng)`` per-region emission
    spectrum) varies across legs — every other input is fixed, so a k_eff
    move is attributable to ``chi`` alone. Cheap, discriminating quadrature:
    we measure STABILITY under a χ mutation, not a converged absolute k_eff,
    so low order suffices. Kwargs are passed explicitly (no ``**dict`` splat)
    so each flows into its typed parameter.
    """
    sol = solve_peierls_mg(
        geometry,
        radii=np.array([0.5, 1.0]),
        sig_t=np.array([[1.0, 1.2], [0.9, 1.1]]),
        sig_s=np.array([
            [[0.30, 0.20], [0.00, 0.80]],   # inner
            [[0.25, 0.15], [0.00, 0.70]],   # outer
        ]),
        nu_sig_f=np.array([[0.40, 0.60], [0.00, 0.00]]),   # outer NON-fissile
        chi=chi,
        boundary="vacuum",
        n_panels_per_region=2, p_order=3, n_angular=12, n_rho=12,
        n_surf_quad=12, dps=15, max_iter=300, tol=1e-10,
    )
    if sol.k_eff is None:  # the eigenvalue solve always converges here
        pytest.fail("solve_peierls_mg returned k_eff=None (no eigenvalue)")
    return sol.k_eff


@pytest.mark.parametrize("geometry", [
    pytest.param(CYLINDER_1D, id="cylinder-1d"),
    pytest.param(SPHERE_1D, id="sphere-1d"),
])
def test_nonfissile_region_chi_does_not_affect_keff(geometry):
    """NEGATIVE leg / teeth: a non-fissile region's χ must not move k_eff.

    Three outer-region χ choices — physical [1,0], swapped [0,1], and
    zeroed [0,0] — must all give the same k_eff because the outer region
    emits no fission neutrons. Sink-indexed χ (ERR-063) fails this; the
    source-indexed fix passes it.
    """
    k = {}
    for tag, outer_chi in [
        ("phys", np.array([1.0, 0.0])),
        ("swap", np.array([0.0, 1.0])),
        ("zero", np.array([0.0, 0.0])),
    ]:
        k[tag] = _solve_two_region(
            geometry, chi=np.array([[1.00, 0.00], outer_chi]),
        )

    base = k["phys"]
    for tag in ("swap", "zero"):
        drift = abs(k[tag] - base) / max(abs(base), 1e-30)
        np.testing.assert_array_less(
            drift, _INVARIANCE_TOL,
            err_msg=(
                f"[{geometry.kind}] non-fissile outer χ ({tag}) moved "
                f"k_eff: k_phys={base:.10f}, k_{tag}={k[tag]:.10f}, "
                f"rel drift={drift:.3e} > {_INVARIANCE_TOL:.1e}. "
                f"A region with νΣ_f=0 emits no fission neutrons; its χ "
                f"must be inert (Hébert 2009 Eq. 3.57/3.58). Sink-indexed "
                f"χ (ERR-063) leaks it into the fissile region's birth."
            ),
        )


def test_fissile_region_chi_genuinely_affects_keff():
    """POSITIVE leg / non-degeneracy: the FISSILE region's χ DOES matter.

    Guards against a degenerate "χ never affects k_eff" reading of the
    negative leg — proving the negative leg has teeth only because the
    mutated region is non-fissile, not because χ is globally inert. The
    outer region is non-fissile (χ pinned to [0,0]); only the INNER
    (fissile) region's χ varies.
    """
    k_a = _solve_two_region(CYLINDER_1D, chi=np.array([[1.0, 0.0], [0.0, 0.0]]))
    k_b = _solve_two_region(CYLINDER_1D, chi=np.array([[0.3, 0.7], [0.0, 0.0]]))
    rel = abs(k_a - k_b) / max(abs(k_a), 1e-30)
    np.testing.assert_array_less(
        1e-3, rel,
        err_msg=(
            f"The FISSILE region's χ must affect k_eff (non-degeneracy "
            f"guard): χ_inner=[1,0]->{k_a:.6f}, χ_inner=[0.3,0.7]->"
            f"{k_b:.6f}, rel={rel:.3e}. If equal, the negative leg's "
            f"invariance is vacuous."
        ),
    )
