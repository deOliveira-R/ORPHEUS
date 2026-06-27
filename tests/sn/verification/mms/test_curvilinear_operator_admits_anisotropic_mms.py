r"""ERR-026 / ERR-058 catcher (anisotropic complement) — the curvilinear
``(L+C).apply`` matvec admits a NON-FLAT-in-μ manufactured solution per ordinate.

Promoted 2026-06-26 from
``derivations/diagnostics/diag_sphere_matvec_nonflat_per_ordinate.py``
(numerics-investigator) during the curvilinear snapshot-cleanup gate that
re-baselined the stale pre-ERR-058 SPH apply snapshots (#195/#209).

This is the ANISOTROPIC complement to
:mod:`tests.sn.verification.mms.test_curvilinear_operator_admits_mms` (the
isotropic gate).  The isotropic ansatz :math:`\psi_{\text{ref},n}=A(r)/W` is
FLAT in μ, so it identically nulls the spherical angular-redistribution term
:math:`(1-\mu^2)B(r)/r` — the M-M closure's hardest math (Mode-7 blind spot,
lessons L6/L27).  This gate uses a NON-FLAT ansatz
:math:`\psi_{\text{ref},n}(r)=(A(r)+B(r)\mu_n)/W` that ACTIVATES that term,
applies the PRODUCTION ``(L+C).apply`` DIRECTLY (not via the solver/inverse),
and demands the PER-ORDINATE residual against a hand-derived CONTINUOUS
reference decay under mesh refinement.

Structural independence (vv-principles L11): the reference is the continuous
spherical SN streaming+removal operator evaluated analytically,

.. math::

   [(L+C)\psi]_n(r) = \mu_n A'(r) + \mu_n^2 B'(r) + (1-\mu_n^2)\,B(r)/r
                      + \Sigma_t\,\bigl(A(r) + \mu_n B(r)\bigr)
                      \quad (\text{all}\ /W),

derived from the PDE (MMS pillar) — NO shared discrete ``c_in``/``c_out``/τ/dAw
recurrence or reduction tree with the M-M closure under test.

Why per-ordinate + volume-weighted: the redistribution α-dome telescopes under
the angular weight sum, so a scalar residual is structurally blind to a
half-angle-thread defect (vv anti-pattern #8); the pole-adjacent cells carry a
bounded non-decaying POINTWISE residual (closure truncation × the
:math:`\Delta A/V \sim 1/h` geometry on :math:`r^2\,dr`-vanishing cells), so the
volume-weighted per-ordinate L2 is the stability-relevant norm.  The rate is the
documented inherent ~O(h) of the M-M half-angle-thread interpolation on a
per-ordinate-imposed non-flat field (the #229 angular floor); the SOLVER's
converged scalar flux still reaches O(h²) (V∼r² suppresses the floor — pinned by
``test_curvilinear_aniso_convergence.py``), so this gate asserts CLEAN MONOTONE
DECAY + small magnitude, NOT O(h²).

Teeth (mutation-verified 2026-06-26): a ``c_in`` sign-flip drives the residual to
~7.5 (order ~0); a ``c_out``×2 factor error to ~3.8 (order ~0) — 3-4 orders above
the ~1.7e-3 baseline and non-decaying, so a wrong redistribution sign/factor
(ERR-026 / ERR-058 class) reddens the decay + magnitude assertions below.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_anisotropic_mms_case,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.timed_full_field import TimedFullField


def _lc_apply_on_psi_ref(case, nc: int):
    """Build ψ_ref = (A + B·μ)/W, apply the production ``(L+C).apply``, and
    return ``(residual_per_ordinate, V)`` where
    ``residual = (L+C)ψ_ref − hand_continuous_q``.
    """
    mesh = case.build_mesh(nc)
    sn_mesh = SNMesh(mesh, case.quadrature, case.materials)

    r = mesh.centers                                   # (nx,)
    mu = case.quadrature.mu_x                           # (N,)
    W = float(case.quadrature.weights.sum())
    N = case.quadrature.N
    nx = mesh.N
    ng = 1

    # ── ψ_ref,n(r) = (A(r) + B(r)·μ_n) / W  (NON-FLAT in μ) ──
    A_ = case.A(r)
    B_ = case.B(r)
    Ap_ = case.Ap(r)
    Bp_ = case.Bp(r)
    psi_vals = (A_[None, :] + B_[None, :] * mu[:, None]) / W   # (N, nx)
    vals = np.zeros((N, ng, nx))
    vals[:, 0, :] = psi_vals

    # zero boundary trace — the matvec's bulk action is what we probe.
    zero = TimedFullField.zeros(
        bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
    )
    psi_ref = TimedFullField(
        bulk=AngularFlux.from_mesh(vals, sn_mesh), boundary=zero.boundary,
    )

    # ── Production (L+C).apply on ψ_ref ──
    sigma_t = np.full((ng, nx), case.sigma_t)
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
    lc_psi = (L + C).apply(psi_ref).bulk.values         # (N, ng, nx)

    # ── Hand CONTINUOUS reference: [(L+C)ψ]_n(r) per ordinate, /W ──
    #   streaming  : μ A' + μ² B'
    #   redistrib  : (1 − μ²) B / r        <-- the spherical angular term
    #   removal    : Σ_t (A + μ B)          <-- full collision (NO scatter)
    stream = mu[:, None] * Ap_[None, :] + (mu[:, None] ** 2) * Bp_[None, :]
    redist = (1.0 - mu[:, None] ** 2) * (B_ / r)[None, :]
    removal = case.sigma_t * (A_[None, :] + mu[:, None] * B_[None, :])
    hand_q = ((stream + redist + removal) / W)[:, None, :]   # (N, ng, nx)

    residual = lc_psi - hand_q                          # (N, ng, nx)
    return residual, sn_mesh.volumes


def _vol_weighted_per_ordinate_rms(residual: np.ndarray, V: np.ndarray) -> float:
    """RMS volume-weighted per-ordinate residual (matches the isotropic
    operator-admits gate's norm)."""
    return float(np.sqrt(np.einsum("x,ngx->", V, residual ** 2) / V.sum()))


@pytest.mark.l1
@pytest.mark.catches("ERR-026", "ERR-058")
def test_sphere_lc_matvec_nonflat_per_ordinate_decays_under_refinement():
    """PRODUCTION ``(L+C).apply`` on a NON-FLAT-in-μ ψ must converge to the
    continuous spherical operator PER ORDINATE, under mesh refinement.

    This activates the ``(1-μ²)B/r`` angular-redistribution term (zero for
    flat-flux).  The volume-weighted per-ordinate residual must DECAY
    MONOTONELY at ≥ O(h) and land small; a wrong angular-redistribution
    sign/factor/index lands as a NON-decaying O(1)-class residual (ERR-026 /
    ERR-058 class — mutation-verified: ``c_in`` flip → ~7.5 order~0,
    ``c_out``×2 → ~3.8 order~0).  See the module docstring for the rate
    rationale (inherent ~O(h) of the half-angle-thread interpolation).
    """
    case = build_spherical_anisotropic_mms_case(n_ordinates=16)

    r40, V40 = _lc_apply_on_psi_ref(case, 40)
    r80, V80 = _lc_apply_on_psi_ref(case, 80)
    r160, V160 = _lc_apply_on_psi_ref(case, 160)

    e40 = _vol_weighted_per_ordinate_rms(r40, V40)
    e80 = _vol_weighted_per_ordinate_rms(r80, V80)
    e160 = _vol_weighted_per_ordinate_rms(r160, V160)
    order_lo = float(np.log2(e40 / e80))
    order_hi = float(np.log2(e80 / e160))

    # Clean monotone decay at ≥ O(h) (a wrong closure gives order ~0).
    assert order_lo > 0.85 and order_hi > 0.85, (
        f"Sphere (L+C).apply per-ordinate residual on a NON-FLAT-in-μ field "
        f"does not decay at ≥ O(h) under refinement "
        f"(e40={e40:.3e}, e80={e80:.3e}, e160={e160:.3e}; "
        f"orders {order_lo:.2f}, {order_hi:.2f}).  A wrong angular-"
        f"redistribution sign/factor/index in the M-M closure (ERR-026 / "
        f"ERR-058 class) lands here as a NON-decaying O(1)-class residual "
        f"(mutation-verified: c_in flip → ~7.5 order~0; c_out×2 → ~3.8 order~0)."
    )
    # Magnitude band: baseline lands ~1.7e-3; a wrong closure is O(1)-class.
    assert e160 < 5.0e-3, (
        f"Sphere (L+C).apply finest per-ordinate residual e160={e160:.3e} "
        f"exceeds 5e-3 — a wrong redistribution closure lands O(1)-class "
        f"(mutation-verified ~3.8-7.5)."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-026")
def test_sphere_lc_matvec_redistribution_term_is_active():
    """Guard the Mode-7 declaration: the ``(1-μ²)B/r`` redistribution term must
    be ACTIVE (non-trivial) in this ansatz, else the per-ordinate check above is
    a vacuous tautology that a no-op closure would satisfy.
    """
    case = build_spherical_anisotropic_mms_case(n_ordinates=16)
    mesh = case.build_mesh(40)
    r = mesh.centers
    mu = case.quadrature.mu_x
    W = float(case.quadrature.weights.sum())
    B_ = case.B(r)

    # The continuous redistribution contribution that ONLY the B·μ term sources.
    redist = (1.0 - mu[:, None] ** 2) * (B_ / r)[None, :] / W   # (N, nx)
    peak = float(np.max(np.abs(redist)))
    assert peak > 1e-3, (
        f"spherical redistribution term (1-μ²)B/r peak={peak:.3e} ≤ 1e-3 — "
        f"the non-flat ansatz does NOT activate the angular-redistribution "
        f"path; the per-ordinate matvec check is a vacuous tautology."
    )
