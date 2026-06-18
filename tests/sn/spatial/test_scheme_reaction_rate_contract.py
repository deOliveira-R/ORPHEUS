r"""Diffusion-readiness contract gate — the scheme is a Σ-stateless generic
advection–reaction discretization (#240 Phase 2 Step D4).

Why this gate exists
====================

Step/DD/LD are generic advection–reaction spatial discretizations (Step ↔
first-order upwind, DD ↔ central / Keller box, LD ↔ DG-P1-upwind).  Their
coefficient triple ``(a, inverse_denom, w)`` is a pure function of the
wave-speed (``|μ| → a``) and the **reaction-rate** (``Σ_t → r``); the
``w(Σ)`` blend is the CFD Péclet / κ-scheme blend, NOT an SN artefact.  This
``DiscretizationScheme`` layer is the **model-agnostic advection–reaction
spatial discretization** the diffusion solver will consume (standalone AND
the consistent-DSA preconditioner, #2).  For that consumption to be clean,
the scheme MUST be a *reaction-rate-parameterized* interface that holds NO
cross-section state of its own — a diffusion/DSA consumer supplies an
arbitrary reaction-rate (an off-diagonal removal ``Σ_r = Σ_t − Σ_s0``, a
diffusion-removal, or a pure-advection ``Σ → 0``) and gets coefficients
matching the closed form evaluated AT THAT reaction-rate.

This is a **software-contract gate** (the diffusion-readiness contract), not
an equation verification, so it carries ``@pytest.mark.foundation`` and NO
``verifies(...)`` label.

What this gate pins
===================

1. **POSITIVE — the closed form holds at an ARBITRARY reaction-rate.**  For
   ``DiamondDifference`` AND ``LinearDiscontinuous``, calling
   :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.affine_scan_coefficients`
   with an arbitrary ``sig_t`` (a removal-like ``Σ_r < Σ_t`` value, a
   near-zero pure-advection ``Σ → 0`` value) returns ``(a, inverse_denom, w)``
   matching the closed form evaluated at that ``sig_t``:

   * **DD** — ``w = ½`` always (central / box scheme, Σ-independent);
     ``denom = Σ·V + 2|μ|·A_down`` (slab ``A_down = 1``), so
     ``inverse_denom = 1/denom`` and ``a = 2|μ|·A_total/denom − 1``.
   * **LD** — ``w = 1/(1+k)`` with ``k = (|μ|/θ)/(Σ·V + |μ|/θ)`` (``m = |μ|·
     A_down``, ``p = m/θ``, ``D₂ = Σ·V + p``, ``k = p/D₂``, ``S = (Σ·V + m)
     + m·p/D₂``).  The **Péclet limits**: ``w → ½`` as ``Σ → 0`` (the
     pure-advection central limit) and ``w → 1`` as ``Σ → ∞`` (the full
     upwind limit).

2. **NEGATIVE / statelessness teeth — the scheme holds NO Σ.**  A
   single reused scheme instance, called twice with DIFFERENT ``sig_t``,
   returns coefficients that differ ONLY per the formula (determinism +
   Σ-responsiveness; no hidden Σ memory between calls).  And the constructed
   scheme carries NO cross-section attribute (``vars(scheme)`` is empty) — the
   regression guard against a future hidden-Σ smuggle (the "frame-leak
   parameter" the cross-domain-attacker flagged).

Mode-8 safe: every assertion is a function call (``np.testing.assert_*`` /
``pytest.fail``), never a bare ``assert`` — so the gate fires under ``-O``
(the canonical ORPHEUS invocation strips bare ``assert`` to a no-op;
see ``vv-principles`` failure-mode 8).

Epistemic boundary — what this gate does NOT verify
===================================================

The closed forms above **deliberately re-encode** the production
``affine_scan_coefficients`` algebra; they are NOT a structurally-independent
derivation of the coefficient *values*.  So these POSITIVE tests are NOT a
value-verification of the discretization — that would be circular (vv-principles
L11).  They are a **parameterization** oracle: a fixed re-encoded formula is the
correct instrument for the thing this gate actually pins — that the coefficients
are a pure function of the *passed* reaction-rate, with NO instance Σ state
(the diffusion-readiness contract).  Coefficient-VALUE correctness is inherited
upstream from the structurally-independent oracles that already gate this code:
the LD linear-exactness MMS (``test_linear_flux_recovered_exactly`` /
``test_mms_ld_slab`` O(h²)) and the ``group3 ≡ group2`` two-paths cross-check.
The load-bearing assertions here are the **statelessness teeth** (the NEGATIVE
item) — those genuinely discriminate a stateful from a stateless scheme.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous


# ── Slab-neutral geometry fixtures (the wave-speed / geometry inputs held
#    fixed so the reaction-rate Σ is the only varying axis) ────────────────
#
# One ordinate (|μ| = 0.6), one group, one cell.  Slab neutral element:
# A_down = 1, A_total = 2, ΔA/w = 0, c_out = 0 (the `slab_streaming` neutral
# values — LD's curvilinear guard demands dA_w == 0 and c_out == 0 exactly).
# An arbitrary non-power-of-2 cell volume so any FP mismatch surfaces.

_ABS_MU = 0.6
_V = 0.7
_A_DOWN = 1.0
_A_TOTAL = 2.0
_N, _NG, _NX = 1, 1, 1


def _slab_geometry() -> dict[str, np.ndarray]:
    """The fixed wave-speed / geometry inputs (everything but ``sig_t``)."""
    return {
        "abs_mu": np.array([_ABS_MU]),
        "A_down": np.full((_N, _NX), _A_DOWN),
        "A_total": np.full((_N, _NX), _A_TOTAL),
        "dA_w": np.zeros((_N, _NX)),
    }


def _sig_t(value: float) -> np.ndarray:
    """A ``(N, ng, nx)`` reaction-rate field at the given Σ."""
    return np.full((_N, _NG, _NX), value)


def _dd_coeffs(scheme: DiamondDifference, sigma: float):
    geo = _slab_geometry()
    return scheme.affine_scan_coefficients(
        c_out=np.zeros((_N,)),          # DD: c_out is (N,)
        V=np.full((_N, _NX), _V),
        reaction_xs=_sig_t(sigma),
        **geo,
    )


def _ld_coeffs(scheme: LinearDiscontinuous, sigma: float):
    geo = _slab_geometry()
    return scheme.affine_scan_coefficients(
        c_out=np.zeros((_N, _NX)),      # LD: c_out is (N, nx)
        V=np.full((_N, _NX), _V),
        reaction_xs=_sig_t(sigma),
        **geo,
    )


# ── Closed-form references (the reaction-rate-parameterized formula, written
#    out independently of the production reduction tree) ─────────────────────

def _dd_closed_form(sigma: float) -> tuple[float, float, float]:
    r"""DD ``(a, inverse_denom, w)`` at reaction-rate ``sigma`` — slab.

    ``denom = Σ·V + 2|μ|·A_down``; ``w = ½``; ``a = 2|μ|·A_total/denom − 1``.
    """
    denom = sigma * _V + 2.0 * _ABS_MU * _A_DOWN
    a = 2.0 * _ABS_MU * _A_TOTAL / denom - 1.0
    return a, 1.0 / denom, 0.5


def _ld_closed_form(sigma: float, theta: float) -> tuple[float, float, float]:
    r"""LD ``(a, inverse_denom, w)`` at reaction-rate ``sigma`` — slab.

    ``m = |μ|·A_down``; ``p = m/θ``; ``D₂ = Σ·V + p``; ``k = p/D₂``;
    ``S = (Σ·V + m) + m·p/D₂``; ``w = 1/(1+k)``; ``a = m(1+k)²/S − k``.
    """
    m = _ABS_MU * _A_DOWN
    t = sigma * _V
    p = m / theta
    d2 = t + p
    k = p / d2
    s = (t + m) + m * p / d2
    a = m * (1.0 + k) ** 2 / s - k
    return a, 1.0 / s, 1.0 / (1.0 + k)


# Reaction-rate probe values: a removal-like Σ_r < Σ_t, and a near-zero
# pure-advection Σ → 0.  Σ_t itself is not a fixture — the scheme never holds
# it; these are arbitrary reaction-rates a diffusion/DSA consumer could pass.
_REMOVAL_LIKE = 0.35            # Σ_r = Σ_t − Σ_s0 < Σ_t
_PURE_ADVECTION = 1.0e-13       # Σ → 0


# ═══════════════════════════════════════════════════════════════════════
# POSITIVE — the closed form holds at an arbitrary reaction-rate
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.foundation
@pytest.mark.parametrize(
    "sigma", [_REMOVAL_LIKE, _PURE_ADVECTION], ids=["removal_like", "pure_advection"],
)
def test_dd_coefficients_match_closed_form_at_arbitrary_reaction_rate(sigma: float) -> None:
    r"""DD ``(a, inverse_denom, w)`` equals the reaction-rate closed form.

    The scheme is supplied an ARBITRARY reaction-rate ``sigma`` (a value the
    mesh never owns — a removal-like Σ_r, then a pure-advection Σ → 0) and the
    returned coefficients match ``denom = Σ·V + 2|μ|·A_down``, ``w = ½``,
    ``a = 2|μ|·A_total/denom − 1``.  This is the diffusion-readiness positive:
    coefficients are a pure function of the passed reaction-rate.
    """
    a, inv, w = (x.item() for x in _dd_coeffs(DiamondDifference(), sigma))
    a_ref, inv_ref, w_ref = _dd_closed_form(sigma)
    np.testing.assert_allclose(a, a_ref, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(inv, inv_ref, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(w, w_ref, rtol=0.0, atol=0.0)  # DD w = ½ exactly


@pytest.mark.foundation
@pytest.mark.parametrize(
    "sigma", [_REMOVAL_LIKE, _PURE_ADVECTION], ids=["removal_like", "pure_advection"],
)
def test_ld_coefficients_match_closed_form_at_arbitrary_reaction_rate(sigma: float) -> None:
    r"""LD ``(a, inverse_denom, w)`` equals the reaction-rate closed form.

    The scheme is supplied an ARBITRARY reaction-rate ``sigma`` and the
    returned coefficients match the Schur-reduced slab LD closed form:
    ``w = 1/(1+k)``, ``k = (|μ|/θ)/(Σ·V + |μ|/θ)``, ``inverse_denom = 1/S``,
    ``a = m(1+k)²/S − k``.  Pure function of the passed reaction-rate.
    """
    scheme = LinearDiscontinuous()
    a, inv, w = (x.item() for x in _ld_coeffs(scheme, sigma))
    a_ref, inv_ref, w_ref = _ld_closed_form(sigma, scheme.theta)
    np.testing.assert_allclose(a, a_ref, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(inv, inv_ref, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(w, w_ref, rtol=1e-14, atol=0.0)


@pytest.mark.foundation
def test_ld_blend_weight_peclet_limits() -> None:
    r"""LD ``w(Σ)`` realizes the CFD Péclet / κ-scheme blend.

    ``w → ½`` as ``Σ → 0`` (the pure-advection central limit) and ``w → 1`` as
    ``Σ → ∞`` (the full upwind limit).  This is the κ-scheme blend
    ``w = 1/(1+k)`` with ``k = (|μ|/θ)/(Σ·V + |μ|/θ)``: as ``Σ → 0`` the
    reaction term vanishes (``k → 1``, ``w → ½``); as ``Σ → ∞`` the reaction
    dominates (``k → 0``, ``w → 1``).  This is the generic advection–reaction
    structure, NOT an SN artefact — diffusion will consume the SAME blend.
    """
    scheme = LinearDiscontinuous()
    # Σ → 0: the pure-advection / central limit, w → ½.
    _, _, w_zero = (x.item() for x in _ld_coeffs(scheme, 1.0e-13))
    np.testing.assert_allclose(w_zero, 0.5, rtol=0.0, atol=1e-12)
    # Σ → ∞: the full-upwind limit, w → 1.
    _, _, w_inf = (x.item() for x in _ld_coeffs(scheme, 1.0e9))
    np.testing.assert_allclose(w_inf, 1.0, rtol=0.0, atol=1e-8)
    # Monotone in between: an intermediate Σ sits strictly between the limits.
    _, _, w_mid = (x.item() for x in _ld_coeffs(scheme, 1.0))
    if not (0.5 < w_mid < 1.0):
        pytest.fail(
            f"LD blend weight w={w_mid} at Σ=1 is not strictly between the "
            f"Péclet limits (½, 1) — the κ-scheme blend is broken."
        )


@pytest.mark.foundation
def test_dd_blend_weight_is_reaction_rate_independent() -> None:
    r"""DD ``w = ½`` for EVERY reaction-rate — the central / box scheme.

    DD is the central (Keller-box) advection scheme: its cell-average blend is
    the symmetric diamond mean ``ψ̄ = ½(ψ_in + ψ_out)`` regardless of the
    reaction-rate.  Unlike LD's Péclet-blended ``w(Σ)``, DD's ``w`` is the
    fixed ``½`` at Σ → 0, at a removal-like Σ_r, and at Σ → ∞ — the
    Σ-independence that makes its negative-flux pathology the central ``Pe > 2``
    wiggle.
    """
    scheme = DiamondDifference()
    for sigma in (1.0e-13, _REMOVAL_LIKE, 1.0e9):
        _, _, w = (x.item() for x in _dd_coeffs(scheme, sigma))
        np.testing.assert_allclose(w, 0.5, rtol=0.0, atol=0.0)


# ═══════════════════════════════════════════════════════════════════════
# NEGATIVE / statelessness teeth — the scheme holds NO Σ
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.foundation
@pytest.mark.parametrize(
    "make_scheme, coeffs_at",
    [
        (DiamondDifference, _dd_coeffs),
        (LinearDiscontinuous, _ld_coeffs),
    ],
    ids=["DiamondDifference", "LinearDiscontinuous"],
)
def test_scheme_holds_no_cross_section_state(make_scheme, coeffs_at) -> None:
    r"""The constructed scheme carries NO cross-section attribute.

    The lock-in regression guard against a future hidden-Σ smuggle: a scheme
    is a *stateless* discretization strategy.  ``vars(scheme)`` MUST be empty —
    if a future edit stores a ``self.sigma_t`` / ``self._materials`` on the
    instance, this gate fails and the diffusion-readiness contract (the scheme
    is parameterized by the passed reaction-rate, never holds one) is restored
    at review time.
    """
    scheme = make_scheme()
    instance_state = vars(scheme)
    if instance_state:
        pytest.fail(
            f"{make_scheme.__name__} carries instance state {instance_state!r}; "
            f"the diffusion-readiness contract requires a Σ-stateless scheme "
            f"(coefficients a pure function of the passed reaction-rate)."
        )


@pytest.mark.foundation
@pytest.mark.parametrize(
    "make_scheme, coeffs_at, closed_form",
    [
        (DiamondDifference, _dd_coeffs, lambda s: _dd_closed_form(s)),
        (
            LinearDiscontinuous,
            _ld_coeffs,
            lambda s: _ld_closed_form(s, LinearDiscontinuous().theta),
        ),
    ],
    ids=["DiamondDifference", "LinearDiscontinuous"],
)
def test_one_instance_is_pure_in_reaction_rate(make_scheme, coeffs_at, closed_form) -> None:
    r"""A reused instance is a pure function of the passed reaction-rate.

    Call the SAME scheme instance twice with DIFFERENT reaction-rates A then B;
    the two coefficient triples must (a) DIFFER (Σ-responsiveness — the scheme
    is not ignoring its reaction-rate argument), and (b) each match the closed
    form evaluated at its own Σ (determinism + no hidden Σ memory carried
    between calls).  This is the statelessness lock: a hidden ``self.sigma``
    would either leak A's value into B's call (breaking (b)) or freeze the
    result (breaking (a)).
    """
    scheme = make_scheme()
    sigma_a, sigma_b = 0.35, 1.25

    a_a, inv_a, w_a = (x.item() for x in coeffs_at(scheme, sigma_a))
    a_b, inv_b, w_b = (x.item() for x in coeffs_at(scheme, sigma_b))

    # (a) Σ-responsiveness — DIFFERENT reaction-rates ⇒ DIFFERENT inverse_denom.
    #     (inverse_denom = 1/denom is strictly monotone in Σ for both schemes.)
    if inv_a == inv_b:
        pytest.fail(
            f"{make_scheme.__name__}.affine_scan_coefficients returned the same "
            f"inverse_denom={inv_a} for Σ={sigma_a} and Σ={sigma_b}; the scheme "
            f"is not responding to its reaction-rate argument (a hidden frozen Σ?)."
        )

    # (b) Determinism / no hidden memory — each call matches the closed form at
    #     its OWN Σ, with no contamination from the prior call's Σ.
    a_ref_a, inv_ref_a, w_ref_a = closed_form(sigma_a)
    a_ref_b, inv_ref_b, w_ref_b = closed_form(sigma_b)
    np.testing.assert_allclose([a_a, inv_a, w_a], [a_ref_a, inv_ref_a, w_ref_a], rtol=1e-14, atol=0.0)
    np.testing.assert_allclose([a_b, inv_b, w_b], [a_ref_b, inv_ref_b, w_ref_b], rtol=1e-14, atol=0.0)

    # Re-call at Σ_a a SECOND time — must reproduce the first call bit-for-bit
    # (no state mutated by the intervening Σ_b call).
    a_a2, inv_a2, w_a2 = (x.item() for x in coeffs_at(scheme, sigma_a))
    np.testing.assert_array_equal([a_a2, inv_a2, w_a2], [a_a, inv_a, w_a])
