r"""C2 — the same-mesh XS-replacement contrast ladder (P6, #281). L2 tier.

The one L2 gate of the adjoint-weighted collapse battery, in its own file so
the V&V level marker does not conflict with the L0 file-level marker of
``test_homogenization.py`` (qa NIT-1). Everything else about the battery —
§4.0 pins, C1, C3, Cχ — lives there; the spec is
``.claude/plans/p6_adjoint_verification_spec.md`` §4 (B3 addendum g).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

from tests.sn.test_homogenization import _balanced_fissile

pytestmark = [pytest.mark.l2, pytest.mark.cap("solve")]


def _contrast_materials(eps: float):
    """The base material + a contrast-scaled partner: m1(ε) = m0 + ε·Δ.

    ε is the SMALLNESS parameter of the XS-collapse perturbation (δA = O(ε)),
    so first-order perturbation theory predicts gap_fwd = O(ε) and
    gap_adj = O(ε²) — the order signature the C2 ladder measures.
    """
    base_c = np.array([0.20, 0.30]); base_f = np.array([0.10, 0.20])
    base_s = np.array([[0.60, 0.10], [0.00, 0.90]])
    d_c = np.array([0.15, 0.15]); d_f = np.array([-0.05, 0.10])
    d_s = np.array([[-0.10, -0.08], [0.00, -0.05]])
    chi = [1.0, 0.0]; nu = [2.4, 2.4]
    m0 = _balanced_fissile(base_c, base_f, nu, chi, base_s)
    m1 = _balanced_fissile(base_c + eps * d_c, base_f + eps * d_f, nu, chi,
                           base_s + eps * d_s)
    return {0: m0, 1: m1}


@pytest.mark.verifies("sn-homogenization-adjoint-weighted")
class TestC2ComparativeKeffOrder:
    """C2 (redesigned TWICE — see the spec §4 correction notes): the
    SAME-MESH XS-replacement contrast ladder.

    The worth theorems (T0/T5) speak about replacing the fine per-cell XS
    by region-collapsed constants ON THE SAME DISCRETE SYSTEM — so the
    gate re-solves the FINE 16-cell mesh with the per-region collapsed
    materials (no coarse re-discretization: the coarse-mesh DD error is a
    confounder that at coarse partitions swamps and even inverts the worth
    delta — measured, spec §4 note). The smallness knob is the MATERIAL
    CONTRAST ε (m1 = m0 + ε·Δ): δA = O(ε) ⇒ the forward gap is O(ε) while
    the adjoint-weighted gap is O(ε²) (first-order worth identically
    zeroed). A partition ladder is NOT the knob — the alternating-material
    pattern keeps within-region heterogeneity constant at every P, and
    single-material regions null the weight entirely.

    Measured signature (2026-07-26): fwd ratios 2.05/2.01 (first order),
    adj ratios 6.08/9.24 (≥ second order), adjoint gap smaller on every
    rung. k_fine(ε) is each rung's own L1-anchored fine reference
    (anti-#5 pairing)."""

    _EDGES16 = np.linspace(0.0, 4.0, 17)
    _MAT16 = np.tile([0, 1], 8)          # every 2-cell window mixes materials

    def _gaps(self, eps: float) -> tuple[float, float]:
        from orpheus.sn.solver import solve_sn_adjoint

        mats = _contrast_materials(eps)
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        fine = Mesh1D(
            edges=self._EDGES16, mat_ids=self._MAT16,
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"), bc_right=BC("reflective"),
        )
        fwd = solve_sn(mats, fine, quad, scattering_order=0)
        adj = solve_sn_adjoint(mats, fine, quad, scattering_order=0)
        k_fine = fwd.keff
        assert k_fine is not None

        P = 2
        coarse = Mesh1D(
            edges=np.linspace(0.0, 4.0, P + 1), mat_ids=np.zeros(P, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"), bc_right=BC("reflective"),
        )
        centers = 0.5 * (self._EDGES16[:-1] + self._EDGES16[1:])
        region_of = np.clip(
            np.searchsorted(coarse.edges, centers, side="right") - 1, 0, P - 1,
        )
        gaps = {}
        for tag, mm in (
            ("adj", fwd.homogenize(coarse, adjoint=adj)),
            ("fwd", fwd.homogenize(coarse)),
        ):
            # SAME-MESH replacement: the fine geometry, region-constant XS.
            replaced = Mesh1D(
                edges=self._EDGES16, mat_ids=region_of.astype(int),
                coord=CoordSystem.CARTESIAN,
                bc_left=BC("vacuum"), bc_right=BC("reflective"),
            )
            k = solve_sn(dict(mm.materials), replaced, quad,
                         scattering_order=0).keff
            assert k is not None
            gaps[tag] = abs(k - k_fine)
        return gaps["adj"], gaps["fwd"]

    def test_contrast_ladder_orders_discriminate(self):
        gaps = {eps: self._gaps(eps) for eps in (1.0, 0.5, 0.25)}
        for eps, (ga, gf) in gaps.items():
            assert ga < gf, (
                f"eps={eps}: adjoint gap {ga:.3e} not smaller than forward {gf:.3e}"
            )
        for hi, lo in ((1.0, 0.5), (0.5, 0.25)):
            ratio_adj = gaps[hi][0] / gaps[lo][0]
            ratio_fwd = gaps[hi][1] / gaps[lo][1]
            assert ratio_fwd < 3.0, (
                f"forward gap not first-order: ratio {ratio_fwd:.2f} (expect ~2)"
            )
            assert ratio_adj > 3.0, (
                f"adjoint gap not higher-order: ratio {ratio_adj:.2f} "
                f"(expect ≳4; first-order contamination if ~2)"
            )
            assert ratio_adj > ratio_fwd, "adjoint must shrink strictly faster"
