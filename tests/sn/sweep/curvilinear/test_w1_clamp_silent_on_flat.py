r"""W1 — spherical Morel–Montry τ-clamp removal: invariance + positivity gates.

W1 replaces the spherical M-M weight clamp
``tau_mm[n] = max(0.5, min(1.0, tau_raw))`` with the unclamped
Bailey-Morel-Chang weight ``tau_mm[n] = tau_raw`` (the unique weight
exact-on-linear-in-:math:`\mu`).  The cylinder KEEPS its clamp
(structural ``τ_raw=0`` division-by-zero block), so W1 touches the
sphere closure ONLY.

The load-bearing invariant
==========================

On a **flat-in-:math:`\mu`** angular field (every ordinate equal at a
cell — the converged isotropic / infinite-medium solution), the M-M
half-angle thread is the IDENTITY regardless of :math:`\tau`.  The
closure is ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ``; at the flat fixed
point ``ψ^a_in == ψ_avg`` so ``ψ^a_out == ψ_avg`` for any ``τ ≠ 0``.
Equivalently, the net angular redistribution coefficient
``dA_w·(c_out − c_in) = dA_w·(α_out − α_in)`` and the flat-field cell
coefficient ``denom − dA_w·c_in = 2|μ|A_down + dA_w·(α_out − α_in) +
Σ_t·V`` are both :math:`\tau`-INDEPENDENT (the ``α_out/τ`` and
``(1−τ)/τ·α_out`` terms cancel exactly in real arithmetic).

CRITICAL FP NUANCE (measured 2026-06-13):  the cancellation is exact
in real arithmetic but NOT bit-identical at IEEE-754 — the closure
``(ψ − (1−τ)ψ)/τ`` returns ``ψ`` exactly only ~81 % of the time and
within 1 ULP otherwise (the reduction order ``α_out/τ − (1−τ)/τ·α_out``
differs from ``α_out``).  Over an iterative solve this accumulates to
~1e-12 on a flat-flux sphere.  So:

* Gate :func:`test_flat_field_coefficient_tau_independent` asserts the
  exact-arithmetic cancellation to a **few-ULP** FP bound (NOT zero) —
  this is the unit-level proof that the clamp is silent on flat-in-μ.
* Gate :func:`test_homogeneous_reflective_sphere_iso_unchanged` asserts
  the converged homogeneous-reflective (flat-flux) sphere flux is
  unchanged to the iterative FP tail — NOT bit-identity.  This is the
  converged-solve sibling of the unit invariant.

Together they are the iso-invariance gate: W1 does not disturb the
isotropic curvilinear-sphere physics beyond FP-non-associativity noise.
The aniso-improvement gate lives separately in
``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``
(W1 changes — improves at coarse mesh — the ANISOTROPIC solution).

Why these gates are intrinsic (hold post-change by construction)
===============================================================
The unit gate recomputes BOTH the clamped and the raw :math:`\tau`
from the SAME :func:`spherical_streaming` geometry and compares the two
closures — it does not depend on which :math:`\tau` production ships,
so it is valid both before and after W1 lands.  The converged-solve
gate compares the live solver against a frozen pre-W1 reference flux
captured in this module (the homogeneous-reflective sphere is flat-flux
so its converged value is W1-invariant to FP).
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.geometry.reduced_operator import spherical_streaming
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from tests.sn._test_helpers import (
    curvilinear_homogeneous_mesh as _homogeneous_mesh,
)


# ── tau reconstruction (clamped + raw) from the SAME geometry ──────────


def _clamped_and_raw_tau(quad: Quadrature) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(tau_clamped, tau_raw)`` for a GL quadrature.

    Mirrors :func:`spherical_streaming` exactly so the unit test is
    independent of which :math:`\\tau` production currently ships.
    """
    mu, w, N = quad.mu_x, quad.weights, quad.N
    mu_edge = np.zeros(N + 1)
    mu_edge[0] = -1.0
    for n in range(N):
        mu_edge[n + 1] = mu_edge[n] + w[n]
    tau_raw = np.empty(N)
    tau_clamped = np.empty(N)
    for n in range(N):
        dmu = mu_edge[n + 1] - mu_edge[n]
        raw = (mu[n] - mu_edge[n]) / dmu if abs(dmu) > 1e-15 else 0.5
        tau_raw[n] = raw
        tau_clamped[n] = max(0.5, min(1.0, raw))
    return tau_clamped, tau_raw


def _flat_field_coeff(
    st, A_down: float, total_xs: float, tau: float,
    *, alpha_in: float, alpha_out: float,
) -> float:
    """``denom − dA_w·c_in`` — the cell coefficient on ψ_avg at the flat
    fixed point (where ``ψ^a_in == ψ_avg``).

    Issue #236 Step C: the M-M α dome is no longer packed on
    ``StreamingTerms`` (the τ/α packing was retired); the caller passes
    ``alpha_in`` / ``alpha_out`` from the surviving
    ``ReducedStreamingOperator.alpha_half`` dome.
    """
    abs_mu = st.abs_mu
    dA_w = st.delta_A_over_w
    V = st.volume
    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in
    denom = 2.0 * abs_mu * A_down + dA_w * c_out + total_xs * V
    return denom - dA_w * c_in


def _redist_net(
    st, tau: float, *, alpha_in: float, alpha_out: float,
) -> float:
    """``dA_w·(c_out − c_in)`` — the net flat-field angular redistribution."""
    dA_w = st.delta_A_over_w
    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in
    return dA_w * (c_out - c_in)


# ═══════════════════════════════════════════════════════════════════════
# Gate 1a — closure-unit invariant: τ is irrelevant on a flat-in-μ field
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
def test_flat_field_coefficient_tau_independent():
    r"""W1 [closure invariant] The flat-field M-M cell coefficient + the net
    angular redistribution are :math:`\tau`-INDEPENDENT.

    This is the unit-level proof that "the clamp is silent on flat-in-μ":
    for the 8 ordinates where the clamp actually bites (``tau_raw < 0.5``),
    the flat-field coefficient ``denom − dA_w·c_in`` and the net
    redistribution ``dA_w·(c_out − c_in) = dA_w·(α_out − α_in)`` computed
    with the CLAMPED τ equal those computed with the RAW τ to a
    few-ULP FP-non-associativity bound.

    Catches: a future closure refactor that lets τ leak into the
    flat-field cell algebra (re-floors the isotropic solution — the
    ERR-026 redistribution class).  Measured worst FP residual
    2026-06-13: flat coeff 1.1e-13, redist net 1.8e-14.
    """
    R = 5.0
    nx = 20
    edges = np.linspace(0.0, R, nx + 1)
    mesh = Mesh1D(
        edges=edges, mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(16)
    op = spherical_streaming(mesh, quad)
    tau_clamped, tau_raw = _clamped_and_raw_tau(quad)

    changed = np.where(tau_raw != tau_clamped)[0]
    # Guard the guard: the clamp MUST bite somewhere, else the test is
    # vacuous (it would compare τ against itself).
    assert changed.size > 0, (
        "clamp never bites for S16 — the τ-independence test is vacuous; "
        f"tau_raw={tau_raw}"
    )

    total_xs = 0.8
    mu = quad.mu_x
    worst_coeff = 0.0
    worst_redist = 0.0
    for cell_idx in range(nx):
        for n in changed:
            st = op.streaming_terms(cell_idx, int(n))
            # α dome from the surviving operator array (Step C).
            alpha_in = float(op.alpha_half[int(n)])
            alpha_out = float(op.alpha_half[int(n) + 1])
            mu_n = float(mu[n])
            A_down = (
                float(op.face_areas[cell_idx + 1]) if mu_n > 0
                else float(op.face_areas[cell_idx])
            )
            cc = _flat_field_coeff(
                st, A_down, total_xs, float(tau_clamped[n]),
                alpha_in=alpha_in, alpha_out=alpha_out,
            )
            cr = _flat_field_coeff(
                st, A_down, total_xs, float(tau_raw[n]),
                alpha_in=alpha_in, alpha_out=alpha_out,
            )
            rc = _redist_net(
                st, float(tau_clamped[n]),
                alpha_in=alpha_in, alpha_out=alpha_out,
            )
            rr = _redist_net(
                st, float(tau_raw[n]),
                alpha_in=alpha_in, alpha_out=alpha_out,
            )
            worst_coeff = max(worst_coeff, abs(cc - cr))
            worst_redist = max(worst_redist, abs(rc - rr))

    # Few-ULP bound: the cancellation is exact in real arithmetic; the
    # residual is reduction-order FP-non-associativity in
    # ``α_out/τ − (1−τ)/τ·α_out``.  Bound by ``(reduction depth ≈ a few)
    # × typical magnitude (~1) × eps``; 1e-12 is ~4500 ULP headroom over
    # the measured 1.1e-13, but still 4 orders below any algorithmic
    # (clamp-leak) change which would be O(0.1) per the clamp gap.
    assert worst_coeff < 1e-12, (
        f"flat-field cell coefficient is τ-DEPENDENT (worst Δ={worst_coeff:.3e} "
        f"≥ 1e-12) — the clamp is NOT silent on flat-in-μ; the M-M closure "
        f"leaks τ into the isotropic solution (ERR-026 redistribution class)"
    )
    assert worst_redist < 1e-12, (
        f"net angular redistribution is τ-DEPENDENT (worst Δ={worst_redist:.3e} "
        f"≥ 1e-12) — c_out−c_in ≠ α_out−α_in; τ leaks into the flat-field "
        f"redistribution (ERR-026 class)"
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 1b — converged-solve iso invariant: homogeneous-reflective sphere
# (flat flux) is unchanged by W1 to the iterative FP tail
# ═══════════════════════════════════════════════════════════════════════
#
# The homogeneous-reflective sphere converges to the flat infinite-medium
# solution: φ(r) ≡ const, ψ flat-in-μ.  By Gate 1a the M-M closure is
# τ-irrelevant on that field, so W1 leaves the converged value invariant
# up to the FP tail (the ~1e-12 reduction-order jitter Gate 1a quantifies,
# accumulated over the power iteration).  The reference below is the
# CLAMPED converged value (pre-W1); a frozen literal so this gate holds
# regardless of which τ ships — and so the strong claim "iso is unchanged"
# is checked against an actual number, not against a re-run of the same
# code.  Anchored to the closed-form k_inf = νΣ_f/Σ_a = 1.875 (2G mixture
# A) — a structurally-independent ground (vv-principles L11), NOT a
# solver-vs-solver tautology.

# Frozen CLAMPED reference (pre-W1), GL-8, n=20, 2G mixture A reflective
# sphere.  φ is spatially flat; report the flat value + k_eff.
_ISO_SPHERE_KEFF_CLAMPED = 1.8750000000103633
_ISO_SPHERE_FLUX_FLAT_CLAMPED = 0.1326291192  # per-group flat value (both g)
_K_INF_2G_A = 1.875  # closed-form νΣ_f/Σ_a — structurally-independent anchor


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
def test_homogeneous_reflective_sphere_iso_unchanged():
    r"""W1 [iso invariance] The homogeneous-reflective (flat-flux) sphere
    eigenpair is unchanged to the iterative FP tail.

    This is the converged-solve sibling of Gate 1a.  The flux is
    spatially flat (max/min == 1) so the field is flat-in-μ at the fixed
    point; W1 may shift only at the reduction-order FP tail (~1e-12).
    The gate asserts:

    * ``k_eff`` matches the frozen clamped reference within ``1e-9``
      (the iterative tail; bug-era re-floor would be O(1e-3)+);
    * ``k_eff`` matches the closed-form ``k_inf = 1.875``
      (structurally-independent anchor — proves both the clamped
      reference AND the live value are CORRECT, not merely close);
    * the flux is spatially flat (the flat-in-μ precondition — else the
      invariance claim is vacuous);
    * every group's flat value matches the frozen clamped reference
      within ``1e-9``.

    NOT a bit-identity gate: W1's τ-cancellation reduction order differs
    at IEEE-754 (measured homogeneous-sphere drift 2026-06-13:
    |Δk|=2.3e-13, max|Δφ|=4.4e-13 — both well inside 1e-9).
    """
    mix = get_mixture("A", "2g")
    mesh = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = Quadrature.gauss_legendre(8)
    result = solve_sn(
        {0: mix}, mesh, quad,
        max_outer=500, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=300, inner_tol=1e-10,
    )

    flux = np.asarray(result.scalar_flux.values, dtype=np.float64)  # (ng, nx)

    # Flat-in-μ precondition: the converged scalar flux must be spatially
    # flat, else this homogeneous-reflective case does not exercise the
    # iso invariant.
    for g in range(flux.shape[0]):
        row = flux[g]
        np.testing.assert_allclose(
            row.max(), row.min(), rtol=1e-9,
            err_msg=(
                f"group {g} flux not spatially flat (max/min="
                f"{row.max() / row.min():.6f}); homogeneous-reflective sphere "
                f"should converge to the flat infinite-medium solution — the "
                f"iso-invariance precondition fails"
            ),
        )

    # k_eff: unchanged from the frozen clamped reference (FP tail only) AND
    # correct against the closed-form anchor.
    assert abs(result.keff - _ISO_SPHERE_KEFF_CLAMPED) < 1e-9, (
        f"k_eff={result.keff:.13f} drifted from the frozen clamped reference "
        f"{_ISO_SPHERE_KEFF_CLAMPED:.13f} by "
        f"{abs(result.keff - _ISO_SPHERE_KEFF_CLAMPED):.2e} ≥ 1e-9 — W1 "
        f"disturbed the ISOTROPIC sphere physics beyond the FP tail "
        f"(the clamp was supposed to be silent on flat-in-μ)"
    )
    assert abs(result.keff - _K_INF_2G_A) < 1e-6, (
        f"k_eff={result.keff:.10f} ≠ closed-form k_inf={_K_INF_2G_A} "
        f"(structurally-independent anchor) — the reference itself is wrong"
    )

    # Per-group flat flux: unchanged from the frozen clamped reference.
    for g in range(flux.shape[0]):
        flat = float(flux[g].mean())
        assert abs(flat - _ISO_SPHERE_FLUX_FLAT_CLAMPED) < 1e-9, (
            f"group {g} flat flux={flat:.12e} drifted from clamped reference "
            f"{_ISO_SPHERE_FLUX_FLAT_CLAMPED:.10e} by "
            f"{abs(flat - _ISO_SPHERE_FLUX_FLAT_CLAMPED):.2e} ≥ 1e-9"
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate 4 — positivity on stress configs post-W1
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "groups, radius, n_ord",
    [
        ("1g", 2.0, 16),   # thick absorber, fine angular
        ("2g", 2.0, 64),   # 2G, S64 high-order angular
    ],
    ids=["1g_R2_S16", "2g_R2_S64"],
)
def test_unclamped_sphere_flux_strictly_positive(groups, radius, n_ord):
    r"""W1 [positivity] The converged sphere scalar flux is strictly
    positive on stress configs after the τ-clamp removal.

    The clamp's stated purpose was to keep the M-M weighting positive
    (prevent zero/negative half-angle thread).  W1 removes it; this gate
    proves the converged scalar flux stays strictly positive (and finite)
    on a thick absorber + a high-order (S64) angular config — the
    configurations most likely to expose a negativity if the unclamp were
    unsafe.  Measured 2026-06-13 (unclamped): flux_min 3.98e-2 (1g),
    1.33e-1 (2g S64), both finite.
    """
    mix = get_mixture("A", groups)
    mesh = _homogeneous_mesh(40, radius, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = Quadrature.gauss_legendre(n_ord)
    result = solve_sn(
        {0: mix}, mesh, quad, max_inner=800, inner_tol=1e-11,
    )

    fv = np.asarray(result.scalar_flux.values, dtype=np.float64)
    assert np.all(np.isfinite(fv)), "non-finite flux post-unclamp"
    assert np.all(fv > 0.0), (
        f"non-positive flux post-unclamp: min={fv.min():.4e} — the M-M "
        f"weight unclamp produced a zero/negative half-angle thread "
        f"(the failure mode the clamp guarded against)"
    )
