r"""The curvilinear MMS suite CANNOT adjudicate the per-level ordinate
ordering — this file is that Mode-7 declaration made executable.

Promoted 2026-08-01 from
``derivations/diagnostics/diag_326_mms_ordering_blindness.py``
(numerics-investigator, GitHub issue **#326**).

WHY IT LIVES HERE
-----------------
``test_mms_curvilinear.py`` and ``test_curvilinear_aniso_convergence.py`` are
the gates a fresh session will reach for when asked "does the curvilinear MMS
see X?".  For X = *the per-level ordinate ordering* the answer is **no, and
exactly no** — but nothing in those files says so, and their Mode-7 honest-scope
note exists only as prose.  This module is the executable form of it, so the
claim cannot decay silently.

THE SITUATION
-------------
``rules_product.py`` orders each mu-level by ``np.argsort(mu_x)``.
``eta = mu_x = sin(theta) cos(phi)`` is 2-to-1 over ``phi`` in ``[0, 2*pi)``:
the azimuthal mirror pair ``(phi, 2*pi - phi)`` shares it.  The level was never
totally ordered, and the tie-break — decided today by rounding noise — moves a
heterogeneous cylindrical solve by 0.6 %–7.2 %.

WHAT THIS FILE PROVES
---------------------
1. ``alpha`` and ``tau`` are BIT-IDENTICAL across tie-breaks.  ``alpha`` is a
   partial sum of ``w_m * eta_m`` and ``w`` is constant within a product-rule
   level, so permuting equal ``eta`` cannot move any partial sum; ``tau`` reads
   only the sorted ``eta`` sequence.  The tie-break moves ONLY which ordinate
   sits at which position — i.e. which member of each mirror pair receives
   ``tau = 1`` and which the clamped ``tau = 1/2``.

2. **BOTH production cylindrical MMS ansatzes are functions of ``eta_n`` and
   ``xi_n**2`` ONLY** (``mms/sn.py``: ``psi_n = (A + B*eta_n)/W``;
   ``Q_n = eta A' + eta**2 B' + xi**2 B/r + (SigT-SigS) A + SigT eta B``).  The
   mirror pair shares ``eta`` — that IS the sort-key degeneracy — and shares
   ``xi**2``.  So it carries identical ``psi_ref``, identical ``Q``, identical
   ``w``: swapping it is a PURE RELABELING of two ordinates carrying identical
   data, and ``phi = sum_n w_n psi_n`` is invariant.

   **Therefore the curvilinear MMS is structurally blind to the entire
   ``xi``-ODD sector**, of which the ordering degeneracy is one instance.  That
   is ``vv-principles`` Mode 12 in exact form — the measured functional's
   invariance group CONTAINS the mutation class.  No tolerance, no mesh
   refinement, no quadrature order can expose it through these ansatzes.

3. A ``xi``-ODD companion ansatz ``psi_n = (A(r) + B(r) xi_n)/W`` leaves the
   symmetric sector and DOES see the tie-break (20.6 % apart at nx=10) — but it
   does NOT adjudicate: both orderings converge to the SAME angular floor from
   opposite sides (0.28 % apart at nx=80) at zero spatial order.  There is no
   "right" ordering to find; the closure itself is what is under-determined.

WHAT THIS FILE IS BLIND TO
--------------------------
It verifies no transport equation and carries no ``verifies(...)`` edge.  Every
row is a claim about the TEST DESIGN (an algebraic reduction invariant, a
property of an MMS case object, a blindness measurement), which is what
``foundation`` means.  The criterion that DOES bite the underlying defect is
the geometric ``xi``-mirror symmetry —
``tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py``.

MEASURED 2026-08-01 (branch ``refactor/operator-strategy-layers``)
------------------------------------------------------------------
Isotropic cylindrical MMS, nx [20,40,80,160], the production gate's ladder::

    lexsort  [2.160863081428e-03 5.389276744521e-04 1.346506256551e-04 3.365754646309e-05]
    stable   [2.160863081428e-03 5.389276744522e-04 1.346506256550e-04 3.365754646321e-05]
    max rel diff 3.4e-12, orders identical [2.0034 2.0009 2.0002]

Anisotropic cylindrical MMS, nx [10,20,40,80]::

    both orderings  [0.022096038325 0.019504300884 0.019116118385 0.019037806484]
    max rel diff 8.6e-15
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    _make_1g_mixture,
    build_cylindrical_mms_case,
    build_cylindrical_anisotropic_mms_case,
)
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source

from tests.sn._test_helpers import product_level_ordering, volume_weighted_l2

# Every claim here is about test design / quadrature algebra, not about a
# transport equation — so the whole module is `foundation` and NOTHING carries
# a `verifies(...)` edge.
pytestmark = pytest.mark.foundation


# ── 0. the swap is NON-VACUOUS (control legs) ───────────────────────────

def test_the_two_tie_breaks_really_order_the_level_differently():
    """Control: without this the invariance results below prove nothing."""
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break):
            got[tie_break] = [
                list(map(int, a))
                for a in Quadrature.product(n_mu=4, n_phi=8).level_indices
            ]
    print(f"level0 lexsort={got['lexsort'][0]}  stable={got['stable'][0]}")
    assert got["lexsort"] != got["stable"], (
        f"the two tie-breaks produced the SAME ordering {got} — the "
        f"blindness measurements below would be vacuous"
    )


def test_trig_nodes_hide_the_tie_break_behind_rounding_noise():
    """#325 <-> #326 coupling: today's nodes make the tie-break unreachable.

    With ``np.cos(np.linspace(0, 2pi, n_phi, endpoint=False))`` the mirror
    pair's ``eta`` differ by ~1 ULP, so EVERY sort rule agrees and the
    degeneracy is invisible.  Algebraically-exact nodes (#325) make the ties
    exact and hand the decision to the tie-break — which is why #326 blocks
    #325 rather than the reverse, and why every "ordering alone" measurement
    must be made on exact nodes.
    """
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break, exact_nodes=False):
            got[tie_break] = [
                list(map(int, a))
                for a in Quadrature.product(n_mu=4, n_phi=8).level_indices
            ]
    assert got["lexsort"] == got["stable"], (
        f"trig-evaluated nodes unexpectedly left a REACHABLE tie: {got}"
    )
    production = Quadrature.product(n_mu=4, n_phi=8)      # today's nodes
    level = np.asarray(production.level_indices[0])
    gaps = np.abs(np.diff(production.eta[level]))
    ulp_gaps = gaps[gaps < 1e-12]
    print(f"production eta gaps that should be exact ties: {ulp_gaps}")
    assert ulp_gaps.size > 0 and np.all(ulp_gaps > 0.0), (
        f"expected non-zero ULP-scale gaps standing in for exact ties, "
        f"got {ulp_gaps}"
    )


# ── 1. the tie-break cannot move alpha or tau ───────────────────────────

def test_alpha_and_tau_are_bit_identical_across_tie_breaks():
    """The geometry coefficients do not see the tie-break at all."""
    from orpheus.sn.sweep.pole_angular_closure import morel_montry_tau_per_level

    def coefficients(tie_break):
        with product_level_ordering(tie_break):
            quad = Quadrature.product(n_mu=4, n_phi=8)
            alphas = []
            for level in quad.level_indices:
                eta, w = quad.mu_x[level], quad.weights[level]
                alpha = np.zeros(len(level) + 1)
                for m in range(len(level)):
                    alpha[m + 1] = alpha[m] - w[m] * eta[m]
                alphas.append(alpha)
            taus = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)
        return np.array(alphas), np.array(taus)

    alpha_lex, tau_lex = coefficients("lexsort")
    alpha_stable, tau_stable = coefficients("stable")
    np.testing.assert_array_equal(alpha_lex, alpha_stable)   # bit-identical
    np.testing.assert_array_equal(tau_lex, tau_stable)


def test_the_degenerate_eta_pair_splits_tau_into_one_and_one_half():
    """The MECHANISM: a mirror pair's shared ``eta`` collapses its angular cell.

    ``eta_edge[m+1] = (eta[m] + eta[m+1])/2 = eta[m]`` when the two are equal,
    so ``tau_raw`` is 1 for the lower partner and 0 for the upper, and the
    structural ``[1/2, 1]`` clamp turns that into ``{1, 1/2}``.
    """
    from orpheus.sn.sweep.pole_angular_closure import (
        morel_montry_tau_per_level, morel_montry_tau_raw_per_level,
    )
    quad = Quadrature.product(n_mu=2, n_phi=8)
    raw = morel_montry_tau_raw_per_level(quad, CoordSystem.CYLINDRICAL)[0]
    clamped = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)[0]
    print(f"tau_raw={raw}\ntau_clamped={clamped}")
    np.testing.assert_allclose(raw, [0, 1, 0, 1, 0, 1, 0, 1], atol=1e-14)
    np.testing.assert_allclose(clamped, [0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1])


# ── 2. the ansatzes live inside the symmetric sector ────────────────────

@pytest.mark.parametrize(
    "builder",
    [build_cylindrical_mms_case, build_cylindrical_anisotropic_mms_case],
    ids=["isotropic", "anisotropic"],
)
def test_production_mms_source_is_identical_on_the_mirror_pair(builder):
    """``Q_n`` and ``psi_ref`` depend on ``eta`` and ``xi**2`` only.

    This is the Mode-7 declaration made executable: the ansatz does not merely
    make the ordering error small, it makes the two ordinates the ordering
    permutes CARRY THE SAME DATA.  Everything downstream follows from here.
    """
    case = builder()
    quad = case.quadrature
    partner = quad.reflection_index("y")
    Q = case.external_source(case.build_mesh(20))           # (N, ng, nx)
    np.testing.assert_allclose(Q, Q[partner], rtol=0, atol=1e-15)
    np.testing.assert_array_equal(quad.weights, quad.weights[partner])


# ── 3. the MMS ladders are invariant (the headline negative result) ─────

def _ladder(case, n_cells):
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature,
            case.external_source(mesh), max_inner=800, inner_tol=1e-13,
        )
        errors.append(volume_weighted_l2(
            np.asarray(result.scalar_flux.values)[0, :],
            case.phi_exact(mesh.centers), mesh.volumes,
        ))
    return np.asarray(errors)


@pytest.mark.slow
@pytest.mark.parametrize(
    "builder,n_cells",
    [
        (build_cylindrical_mms_case, [20, 40, 80]),
        (build_cylindrical_anisotropic_mms_case, [10, 20, 40]),
    ],
    ids=["isotropic", "anisotropic"],
)
def test_curvilinear_mms_ladder_is_blind_to_the_tie_break(builder, n_cells):
    """The production MMS gate cannot adjudicate issue #326's ordering."""
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break):
            got[tie_break] = _ladder(builder(), n_cells)
    rel = float(np.max(np.abs(got["lexsort"] - got["stable"])
                       / np.abs(got["lexsort"])))
    print(f"lexsort={got['lexsort']}\nstable ={got['stable']}\nrel={rel:.3e}")
    assert rel < 1e-9, (
        f"the curvilinear MMS unexpectedly SEES the tie-break (rel {rel:.3e}) "
        f"— the Mode-7 blindness declaration this module exists to pin is no "
        f"longer true; re-adjudicate #326 against it"
    )


# ── 4. the xi-ODD companion: sees it, but still cannot adjudicate ───────

@dataclass
class _CylXiOddMMSCase:
    r"""``psi_n(r) = (A(r) + B(r) xi_n)/W`` — the same ``A``, ``B`` shapes as
    the production anisotropic case, but ODD in the azimuthal cosine, so the
    mirror pair carries DIFFERENT data.

    Substituting into ``eta d_r psi - (xi/r) d_omega psi + SigT psi
    = Qext + SigS phi/W`` with ``d_omega xi = eta`` and
    ``phi = sum_n w_n psi_n = A`` (since ``sum_n w_n xi_n = 0``)::

        Qext_n = [eta A' + eta xi B' - eta xi B/r
                  + (SigT - SigS) A + SigT xi B] / W

    ``A(R) = B(R) = 0`` and ``A(0) = B(0) = 0``, so psi vanishes at both r = 0
    (reflective) and r = R (vacuum) exactly as the production case does.
    """

    sigma_t: float
    sigma_s: float
    radius: float
    materials: dict
    mat_id: int
    quadrature: Quadrature

    def A(self, r):
        return np.sin(np.pi * np.asarray(r) / self.radius)

    def Ap(self, r):
        R = self.radius
        return (np.pi / R) * np.cos(np.pi * np.asarray(r) / R)

    def B(self, r):
        rr = np.asarray(r) / self.radius
        return rr * (1.0 - rr) * np.cos(np.pi * rr)

    def Bp(self, r):
        R = self.radius
        rr = np.asarray(r) / R
        return ((1.0 - 2.0 * rr) * np.cos(np.pi * rr) / R
                - (np.pi / R) * rr * (1.0 - rr) * np.sin(np.pi * rr))

    def phi_exact(self, r):
        return self.A(r)

    def build_mesh(self, n_cells: int) -> Mesh1D:
        return Mesh1D(
            edges=np.linspace(0.0, self.radius, n_cells + 1),
            mat_ids=np.full(n_cells, self.mat_id, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"), bc_right=BC("vacuum"),
        )

    def external_source(self, mesh: Mesh1D) -> np.ndarray:
        r = mesh.centers
        A_, Ap_, B_, Bp_ = self.A(r), self.Ap(r), self.B(r), self.Bp(r)
        eta, xi = self.quadrature.eta, self.quadrature.xi
        W = float(self.quadrature.weights.sum())
        Q = (
            eta[:, None] * Ap_[None, :]
            + (eta * xi)[:, None] * Bp_[None, :]
            - (eta * xi)[:, None] * (B_ / r)[None, :]
            + (self.sigma_t - self.sigma_s) * A_[None, :]
            + self.sigma_t * xi[:, None] * B_[None, :]
        ) / W
        return Q[:, None, :]


def _build_xi_odd_case(n_mu: int = 4, n_phi: int = 8) -> _CylXiOddMMSCase:
    return _CylXiOddMMSCase(
        sigma_t=1.0, sigma_s=0.5, radius=5.0,
        materials={1: _make_1g_mixture(1.0, 0.5)}, mat_id=1,
        quadrature=Quadrature.product(n_mu=n_mu, n_phi=n_phi),
    )


def test_xi_odd_companion_source_DOES_differ_on_the_mirror_pair():
    """Control: the companion ansatz really leaves the symmetric sector."""
    case = _build_xi_odd_case()
    partner = case.quadrature.reflection_index("y")
    Q = case.external_source(case.build_mesh(20))
    rel = float(np.max(np.abs(Q - Q[partner])) / np.max(np.abs(Q)))
    print(f"xi-odd companion: mirror-pair source difference {rel:.3e}")
    assert rel > 1e-2, (
        f"the xi-odd companion source is (nearly) mirror-symmetric "
        f"({rel:.3e}) — it would be as blind as the production ansatz"
    )


@pytest.mark.slow
def test_xi_odd_companion_sees_the_tie_break_but_does_not_adjudicate_it():
    """The companion ansatz separates the orderings — then reunites them.

    Measured 2026-08-01 (``product(4, 8)``, nx [10,20,40,80]):
      lexsort [0.0727, 0.0790, 0.0806, 0.0810]  orders [-0.12, -0.03, -0.007]
      stable  [0.0877, 0.0827, 0.0815, 0.0812]  orders [+0.08, +0.02, +0.005]
    20.6 % apart at nx=10, 0.28 % apart at nx=80, spatial order ~0 for both:
    they converge to the SAME angular floor from opposite sides.  Neither is
    "correct"; the floor is the defect.
    """
    n_cells = [10, 20, 40, 80]
    got = {}
    for tie_break in ("lexsort", "stable"):
        with product_level_ordering(tie_break):
            got[tie_break] = _ladder(_build_xi_odd_case(), n_cells)
    coarse = abs(got["lexsort"][0] - got["stable"][0]) / got["lexsort"][0]
    fine = abs(got["lexsort"][-1] - got["stable"][-1]) / got["lexsort"][-1]
    orders = {k: np.log2(v[:-1] / v[1:]) for k, v in got.items()}
    print(f"lexsort={got['lexsort']} orders={orders['lexsort']}")
    print(f"stable ={got['stable']} orders={orders['stable']}")
    print(f"coarse gap={coarse:.4e}  fine gap={fine:.4e}")

    assert coarse > 1e-2, (
        f"the xi-odd companion did NOT separate the orderings ({coarse:.3e}); "
        f"it is then no better an adjudicator than the production ansatz"
    )
    assert fine < coarse / 10.0, (
        f"the orderings did NOT reunite under refinement "
        f"(coarse {coarse:.3e}, fine {fine:.3e}) — one ordering may genuinely "
        f"converge to a different limit; RE-ADJUDICATE #326"
    )
    for tie_break, order in orders.items():
        assert abs(order).max() < 0.5, (
            f"{tie_break}: the xi-odd ladder acquired a spatial convergence "
            f"order {order} — the angular floor lifted; re-adjudicate #326"
        )
