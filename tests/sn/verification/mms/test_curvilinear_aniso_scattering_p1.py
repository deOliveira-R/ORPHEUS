r"""#9 — path-(II) P1 Legendre SCATTERING in curvilinear SN (NEW coverage).

The operator-admits trick at L0: feed a KNOWN anisotropic angular flux
:math:`\psi_{\text{ref},n} = (A(r) + \zeta_n B(r))/W` to the within-group
scattering operator :math:`S` at ``scattering_order=1``, isolate the
:math:`\ell=1` contribution as the operator DIFFERENCE
:math:`S_1.\text{apply}(\psi) - S_0.\text{apply}(\psi)`, and demand it
match a STRUCTURALLY-INDEPENDENT hand-reference **per ordinate** (NOT
weight-summed — the α-dome telescopes under the weight sum, masking
half-angle defects; vv-principles anti-pattern #8).

Why this is the right gate for #9
=================================

The two curvilinear "anisotropic" paths are UNRELATED (see the W4
verification memo + ``discrete_ordinates.rst`` :ref:`pn-scattering`):

* **Path (I)** — geometric angular redistribution
  :math:`(1-\mu^2)/r\,\partial_\mu\psi`, threaded by the Morel–Montry
  Carlson α-dome.  P0-only; the "anisotropy" lives in the angular-flux
  ANSATZ.  The existing curvilinear aniso MMS cases
  (:func:`build_spherical_mms_case` etc.) exercise ONLY this.

* **Path (II)** — the :math:`P_\ell` Legendre SCATTERING moments
  :math:`R\,\Lambda\,M` (:eq:`pn-scatter` / :eq:`flux-moments`), wired
  via :func:`~orpheus.sn.coupled_system.build_within_group_system` (the ``S`` in the ``(L+C), S, B``
  triple carries :math:`P_1` when ``scattering_order=1``).  This is
  geometry-AGNOSTIC, but NO curvilinear test exercised it before #9 —
  the only prior :eq:`pn-scatter` tests are 2-D Cartesian
  (``test_keff_2d.py::TestBicgstabPnScattering``).  This file closes
  that coverage gap by verifying :eq:`pn-scatter`/:eq:`flux-moments`
  in spherical AND cylindrical geometry.

Why this verifies the OPERATOR and may use 1 energy group
=========================================================

This is a flux-shape / OPERATOR claim, NOT an eigenvalue claim: the
per-ordinate :math:`P_1` source reads the :math:`\ell=1` flux moment
:math:`\phi_1` (flux-shape-dependent by construction — a flat flux makes
:math:`\phi_1 = 0`, so the term is INACTIVE; the anisotropic ansatz
:math:`(A + \zeta B)/W` with :math:`B \ne 0` ACTIVATES it).  The
Cardinal Rule bars 1-GROUP **eigenvalue** tests (``k = \nu\Sigma_f /
\Sigma_a`` is flux-shape independent); it does NOT bar a 1-group
operator-source claim.  A single group keeps the hand-reference
transparent and the structural independence sharp.  The directional
EIGENVALUE claim (P1 lowers k_eff) lives at L1 in
``tests/sn/eigenvalue/test_keff_curvilinear.py`` with a 2-group fissile
mixture.

Structural independence of the hand-references
==============================================

* SPHERE — fully SH-table-INDEPENDENT (strongest): with the radial
  cosine :math:`\mu_n = \zeta_n`, the :math:`\ell=1` Legendre source is

  .. math::

     q_n^{P_1} = \frac{1}{W}\,3\,\mu_n\,\Sigma_{s,1}\,\phi_1,
     \qquad
     \phi_1 = B(r)\,\frac{\sum_n w_n \mu_n^2}{W},

  derived from :math:`P_1(\mu) = \mu` directly — it never touches the
  code's ``spherical_harmonics`` table.

* CYLINDER — explicit real-:math:`Y_1^m` moment sum, independent of the
  ``R \Lambda M`` einsum the production code composes:

  .. math::

     q_n^{P_1} = \frac{1}{W}\,3\,\Sigma_{s,1}
       \sum_m Y_1^m(\Omega_n)\,\phi_1^m,
     \qquad
     \phi_1^m = \sum_n w_n\,Y_1^m(\Omega_n)\,\psi_n.

  This reuses ``quad.spherical_harmonics(1)`` (a verified primitive,
  below the trusted-library line — see ``algebra-of-record`` §structural
  independence) but assembles the moment sum by an INDEPENDENT explicit
  Python loop, not the production frame analysis/reconstruction faces /
  ``LegendreMomentScattering`` einsums.

Negative control (vv-principles anti-pattern #11): a dropped-P1
regression makes :math:`S_1 \equiv S_0` so the operator difference is
identically zero — which would fail the non-zero hand-ref match anyway,
but we PIN ``max|S_1 - S_0| > 1e-6`` explicitly so the test is not a
vacuous tautology (measured: sphere ~8.5e-3, cyl ~1.3e-3).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.mms.sn import (
    build_cylindrical_mms_case,
    build_spherical_mms_case,
)
from orpheus.geometry import CoordSystem
from orpheus.sn.solver import SNSolver, _as_sn_mesh
from orpheus.transport.fields.angular_flux import AngularFlux
from tests.sn._test_helpers import curvilinear_homogeneous_mesh


# ── shared builders ─────────────────────────────────────────────────────

#: P1 self-scatter of the 1g "A" mixture (``SigS[1] = [[0.025]]``) — a
#: single scalar so the hand-reference stays transparent.  ``A`` is the
#: only mixture whose 1g entry carries a non-zero :math:`\ell=1` block.
_SIGMA_S1 = 0.025


def _scattering_op(coord: CoordSystem, quad, scattering_order: int):
    """Build the curvilinear ``ScatteringOperator`` at the requested order.

    Routes through the SAME ``_as_sn_mesh`` + ``SNSolver`` path
    ``solve_sn_fixed_source`` uses, so ``solver.scattering_op`` IS the
    ``S`` of the production ``(L+C), S, B`` within-group triple.
    """
    mat = get_mixture("A", "1g")
    mesh = curvilinear_homogeneous_mesh(8, 2.0, mat_id=0, coord=coord)
    sn_mesh = _as_sn_mesh(mesh, quad, {0: mat}, "vacuum", mat_map=None)
    solver = SNSolver(
        sn_mesh,
        inner_solver="source_iteration",
        scattering_order=scattering_order,
        max_inner=5,
        inner_tol=1e-12,
    )
    return solver.scattering_op, sn_mesh, mesh.centers


def _anisotropic_flux(sn_mesh, centers):
    r"""ψ_ref,n = (A(r) + ζ_n B(r))/W with A,B > 0 and ζ = radial cosine.

    ``A(r) > 0`` and ``B(r) != 0`` are both load-bearing: ``A`` keeps the
    isotropic part non-trivial; ``B`` ACTIVATES the :math:`\ell=1` moment
    (a flat ``B = 0`` would null the very term under test).
    """
    quad = sn_mesh.quad
    W = float(quad.weights.sum())
    r = centers
    A = 1.0 + 0.5 * r
    B = 0.3 + 0.2 * r
    zeta = quad.mu_x  # radial cosine in both spherical (dim=1) and cyl
    vals = np.zeros((quad.N, sn_mesh.ng, *sn_mesh.spatial_shape))
    for n in range(quad.N):
        vals[n, 0, :] = (A + zeta[n] * B) / W
    return AngularFlux(values=vals, space=sn_mesh.angular_bulk_space), A, B, W


def _isolated_p1(coord: CoordSystem, quad):
    """Return ``(P1, psi, sn_mesh, A, B, W)`` — the isolated ℓ=1 source.

    ``P1 = S_1.apply(psi) - S_0.apply(psi)`` is the operator difference
    that subtracts off the P0 in-scatter, leaving ONLY the :math:`\\ell=1`
    Galerkin reconstruction.
    """
    S1, sn_mesh, centers = _scattering_op(coord, quad, scattering_order=1)
    S0, _, _ = _scattering_op(coord, quad, scattering_order=0)
    psi, A, B, W = _anisotropic_flux(sn_mesh, centers)
    P1 = S1.apply(psi).values - S0.apply(psi).values
    return P1, psi, sn_mesh, A, B, W


# ── L0 tests ────────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("pn-scatter", "flux-moments", "sn-p1-sphere-hand-ref")
def test_spherical_p1_source_matches_hand_reference():
    r"""[L0] Spherical curvilinear ℓ=1 scattering source == hand-ref.

    #9 path-(II) coverage — verifies :eq:`pn-scatter` / :eq:`flux-moments`
    in SPHERICAL geometry, the first curvilinear exercise of the
    :math:`R\,\Lambda\,M` scattering composition.

    The hand-reference is fully SH-table-INDEPENDENT — it uses
    :math:`P_1(\mu) = \mu` directly with the radial cosine, never the
    production ``spherical_harmonics`` table.  A sign-flipped /
    transposed / mis-normalised :math:`\ell=1` block would land far from
    this reference (this is the curvilinear sibling of the 2-D
    ``TestBicgstabPnScattering`` checks).
    """
    quad = build_spherical_mms_case().quadrature
    P1, _psi, _sn, _A, B, W = _isolated_p1(CoordSystem.SPHERICAL, quad)

    # Negative control (anti-pattern #11): a dropped-P1 makes S1 == S0
    # so the difference is identically zero — the gate must not be a
    # tautology that a no-op P1 satisfies.
    peak = float(np.max(np.abs(P1)))
    assert peak > 1e-6, (
        f"spherical S1-S0 peak={peak:.3e} <= 1e-6 — the P1 source is "
        f"inert; the dropped-P1 regression class (S1 == S0) is NOT "
        f"caught by the hand-ref match alone."
    )

    # Hand-reference: q_n^P1 = (1/W) 3 mu_n Sigma_s1 phi_1,
    #   phi_1 = B(r) * sum_n(w_n mu_n^2) / W.
    mu = quad.mu_x
    wts = quad.weights
    phi1 = B * float(np.sum(wts * mu**2)) / W
    hand = np.zeros_like(P1)
    for n in range(quad.N):
        hand[n, 0, :] = (1.0 / W) * 3.0 * mu[n] * _SIGMA_S1 * phi1

    np.testing.assert_allclose(
        P1, hand, rtol=1e-12, atol=1e-13,
        err_msg=(
            "spherical curvilinear P1 scattering source disagrees with the "
            "SH-table-independent hand-reference per ordinate "
            f"(peak|P1|={peak:.3e}, max abs err="
            f"{float(np.max(np.abs(P1 - hand))):.3e})."
        ),
    )


@pytest.mark.l0
@pytest.mark.verifies("pn-scatter", "flux-moments", "sn-p1-cylinder-hand-ref")
def test_cylindrical_p1_source_matches_hand_reference():
    r"""[L0] Cylindrical curvilinear ℓ=1 scattering source == hand-ref.

    #9 path-(II) coverage — verifies :eq:`pn-scatter` / :eq:`flux-moments`
    in CYLINDRICAL geometry (a genuine 3-D quadrature, dim=3, so all
    three real :math:`Y_1^m` are active).  The hand-reference assembles
    the moment sum by an INDEPENDENT explicit loop over the
    ``spherical_harmonics`` table, NOT the production frame analysis/reconstruction
    faces / ``LegendreMomentScattering`` einsums — so a transposed or
    mis-shaped einsum in the production path is detectable.
    """
    quad = build_cylindrical_mms_case().quadrature
    P1, psi, sn_mesh, _A, _B, W = _isolated_p1(CoordSystem.CYLINDRICAL, quad)

    peak = float(np.max(np.abs(P1)))
    assert peak > 1e-6, (
        f"cylindrical S1-S0 peak={peak:.3e} <= 1e-6 — the P1 source is "
        f"inert; the dropped-P1 regression class (S1 == S0) is NOT "
        f"caught by the hand-ref match alone."
    )

    # Hand-reference: q_n^P1 = (1/W) 3 Sigma_s1 sum_m Y_1^m(Omega_n) phi_1^m,
    #   phi_1^m = sum_n w_n Y_1^m(Omega_n) psi_n.
    # Y shape (N, L+1, 2L+1); the ell=1 block lives at index [:, 1, :],
    # with the m-slot running 0,1,2 (the addition-theorem-shifted index).
    Y = quad.spherical_harmonics(1)
    wts = quad.weights
    psi_vals = psi.values[:, 0, :]  # (N, nx)
    nx = sn_mesh.spatial_shape[0]
    hand = np.zeros_like(P1)
    for n in range(quad.N):
        acc = np.zeros(nx)
        for m_slot in range(3):
            phi1m = np.einsum("n,nx->x", wts * Y[:, 1, m_slot], psi_vals)
            acc += Y[n, 1, m_slot] * phi1m
        hand[n, 0, :] = (1.0 / W) * 3.0 * _SIGMA_S1 * acc

    np.testing.assert_allclose(
        P1, hand, rtol=1e-12, atol=1e-13,
        err_msg=(
            "cylindrical curvilinear P1 scattering source disagrees with "
            "the explicit-SH-moment-sum hand-reference per ordinate "
            f"(peak|P1|={peak:.3e}, max abs err="
            f"{float(np.max(np.abs(P1 - hand))):.3e})."
        ),
    )
