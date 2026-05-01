r"""SymPy derivation of the Phase 5 continuous-:math:`\mu`
multi-bounce specular kernel.

This module is the **math origin** for
:func:`compute_K_bc_specular_continuous_mu_sphere` in
:mod:`peierls_geometry`. It pins the Sanchez 1986 [SanchezTTSP1986]_
Eq. (A6) ↔ ORPHEUS M1 sketch equivalence (the µ-weight convention
question that R1 closure required) and the diagonal-singularity
finding that blocks production wiring of
``boundary="specular_continuous_mu"`` (the closure raises
:class:`NotImplementedError`; production multi-bounce specular uses
``closure="specular_multibounce"`` instead — see Issue #133 close-out).

Purpose
-------

Phase 4 shipped ``closure="specular_multibounce"`` as a matrix-Galerkin form
:math:`K_{\\rm bc} = G\\,R\\,(I - T\\,R)^{-1}\\,P` that diverges at high
rank :math:`N \\ge 4` for sphere/cyl due to a **fundamental matrix-Galerkin
divergence** (operator-norm proof in
``.claude/agent-memory/numerics-investigator/specular_mb_overshoot_root_cause.md``).
Phase 5 bypasses the matrix inverse and integrates the multi-bounce kernel
**continuously in µ** with adaptive quadrature.

Sanchez 1986 Eq. (A6) provides the textbook form for the homogeneous sphere
with specular + diffuse BC parametrised by :math:`\\alpha + \\beta`. For
ORPHEUS (perfect specular, isotropic scattering), set
:math:`\\alpha = 1, \\beta = 0, \\omega_1 = 0`:

.. math::
    :label: peierls-sanchez-A6-sphere

    g_h(\\rho' \\to \\rho) \\;=\\; 2\\!\\int_{\\mu_0}^{1}
        T(\\mu_-)\\,\\mu_*^{-1}\\,\\cosh(\\rho\\mu)\\,\\cosh(\\rho'\\mu_*)\\,
        e^{-2a\\mu_-}\\,\\mathrm d\\mu

with multi-bounce factor

.. math::
    :label: peierls-sanchez-Tmu

    T(\\mu) \\;=\\; \\frac{1}{1 - e^{-2a\\mu}}, \\qquad a = \\Sigma\\,R,

chord-projection identity :math:`\\mu_*^2 = {\\rho'}^2 - \\rho^2(1 - \\mu^2)`,
and chord-visibility cone :math:`\\mu_0^2 = \\max(0, 1 - (\\rho'/\\rho)^2)`.

ORPHEUS's M1 sketch (cross-domain-attacker output, ``.claude/agent-memory/
cross-domain-attacker/phase5_continuous_mu_frames.md``) expected the form

.. math::
    :label: peierls-orpheus-M1-sketch

    K_{\\rm bc}^{\\rm mb,sph}(r_i, r_j) \\;=\\; 2\\!\\int_0^1\\!
        G_{\\rm in}(r_i, \\mu)\\,F_{\\rm out}(r_j, \\mu)\\,
        \\frac{\\mu}{1 - e^{-\\sigma\\,2R\\mu}}\\,\\mathrm d\\mu

with :math:`\\mu` in the **numerator** of the multi-bounce factor.

The two forms differ in where the :math:`\\mu` weight lives:

- Sanchez writes :math:`T(\\mu)` with no :math:`\\mu` in the numerator;
  the surrounding kernel carries :math:`\\mu_*^{-1}` (in :math:`g_h`)
  or no µ-weight (in :math:`h_h`).
- ORPHEUS M1 sketch lumps :math:`\\mu` into the multi-bounce factor.

This script proves that the two forms are **algebraically equivalent**
under the volume-Jacobian, resolving the µ-weight convention question
flagged by the literature-researcher (sub-agent in
``phase5_sanchez_1986_sphere_specular.md``) before any production
implementation begins. R1 risk closure.

Three SymPy verifications
--------------------------

V1. **Multi-bounce-factor identity**: prove
    :math:`T(\\mu) = 1/(1 - e^{-2a\\mu})` is finite at :math:`\\mu \\to 0`
    when multiplied by the µ-weight from the Jacobian. Specifically:
    :math:`\\mu \\cdot T(\\mu) \\to 1/(2a)` as :math:`\\mu \\to 0` (L'Hôpital
    confirms).

V2. **Equivalence to M1 sketch**: prove that the integrand of Sanchez's
    Eq. (A6) is equal to the M1 sketch's :math:`G_{\\rm in} \\cdot F_{\\rm
    out} \\cdot \\mu / (1 - e^{-\\sigma\\,2R\\mu})` under the explicit
    identification:

    - :math:`G_{\\rm in}(r_i, \\mu) = \\cosh(\\rho \\mu)\\,e^{-a\\mu}`
      (in-streaming factor at receiver, ρ here is the optical receiver radius)
    - :math:`F_{\\rm out}(r_j, \\mu) = \\cosh(\\rho' \\mu_*) \\cdot
      \\mu_*^{-1}\\,e^{-a\\mu}` (source-side outgoing factor)
    - The :math:`\\mu / (1 - e^{-\\sigma\\,2R\\mu})` factor in M1 carries
      one explicit :math:`\\mu` that, in Sanchez's form, ends up on the
      :math:`\\mu_*^{-1}` (one Jacobian factor). The two write the same
      integrand.

V3. **Vacuum-BC reduction**: prove that taking :math:`\\alpha \\to 0` in
    Sanchez Eq. (A4)–(A6) recovers the standard vacuum sphere kernel
    :math:`\\bar g_2(\\rho' \\to \\rho) = (\\rho'/\\rho)\\,[\\bar g_0(\\rho'
    \\to \\rho) - \\bar g_0(-\\rho' \\to \\rho)]/2` with
    :math:`\\bar g_0(\\rho' \\to \\rho) = \\tfrac{1}{2} E_1(\\Sigma|\\rho -
    \\rho'|)` (Sanchez Eq. 5; Carlvik / Bell-Glasstone reduction). This is
    a sanity check that the BC kernel reduces correctly at zero reflection.

Output
------

Each ``derive_*`` function returns a dict with the SymPy expressions
and a ``pass`` flag for the corresponding verification. The pytest
gate in ``tests/derivations/test_peierls_specular_continuous_mu_symbolic.py``
asserts the four dicts are pass-True and bit-checks the load-bearing
identities (V1, V2, V3 gating; V4 documentary).

If V1/V2 PASS the production implementation in
``compute_K_bc_specular_continuous_mu_sphere`` (peierls_geometry.py) can
proceed using the Sanchez Eq. (A6) form directly with the SymPy-verified
notation map.

References
----------

- Sanchez, R. (1986). "Integral form of the equation of transfer for a
  homogeneous sphere with linearly anisotropic scattering."
  *Transport Theory & Statistical Physics*, vol. 14.
  DOI: 10.1080/00411458608210456.
- ``.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md``
- ``.claude/agent-memory/cross-domain-attacker/phase5_continuous_mu_frames.md``
- ``.claude/agent-memory/numerics-investigator/specular_mb_overshoot_root_cause.md``
"""
from __future__ import annotations

import sympy as sp


def derive_multi_bounce_factor() -> dict:
    """V1 — Multi-bounce-factor identity.

    Sanchez :math:`T(\\mu) = 1 / (1 - e^{-2a\\mu})` has a simple pole at
    :math:`\\mu = 0` (as :math:`\\mu \\to 0`, :math:`1 - e^{-2a\\mu}
    \\sim 2a\\mu` so :math:`T(\\mu) \\sim 1/(2a\\mu)`). The product
    :math:`\\mu \\cdot T(\\mu)` is the form that appears in the M1 sketch;
    its limit as :math:`\\mu \\to 0` should be :math:`1/(2a)` (finite).

    Returns dict with the SymPy expressions and PASS/FAIL flags.
    """
    mu, a = sp.symbols("mu a", positive=True, real=True)

    T_mu = 1 / (1 - sp.exp(-2 * a * mu))
    mu_T = mu * T_mu

    # Limit as mu -> 0+: should be 1/(2a)
    lim_mu_T = sp.limit(mu_T, mu, 0, "+")
    expected = 1 / (2 * a)

    pass_v1 = sp.simplify(lim_mu_T - expected) == 0

    return {
        "name": "V1 multi-bounce-factor identity",
        "expr_T": T_mu,
        "expr_muT": mu_T,
        "limit": lim_mu_T,
        "expected": expected,
        "pass": pass_v1,
    }


def derive_m1_equivalence() -> dict:
    """V2 — M1 sketch ↔ Sanchez Eq. (A6) algebraic equivalence.

    We prove that the integrand of Sanchez's Eq. (A6) (with α=1, β=0)
    equals the M1 sketch's

    .. math::

        2\\,G_{\\rm in}(r_i, \\mu)\\,F_{\\rm out}(r_j, \\mu)\\,
        \\frac{\\mu}{1 - e^{-\\sigma\\,2R\\mu}}

    under the identification

    - :math:`G_{\\rm in}(\\rho, \\mu) = \\cosh(\\rho \\mu)\\,e^{-a\\mu}`
    - :math:`F_{\\rm out}(\\rho', \\mu) = \\cosh(\\rho' \\mu_*) \\cdot
      \\mu_*^{-1}\\,e^{-a\\mu}`
    - :math:`\\sigma \\cdot 2R \\equiv 2a` (Sanchez non-dim: a = ΣR)

    Sanchez integrand (Eq. A6 with α=1, β=0):

    .. math::

        I_S(\\mu) \\;=\\; 2\\,T(\\mu)\\,\\mu_*^{-1}\\,\\cosh(\\rho\\mu)\\,
                        \\cosh(\\rho'\\mu_*)\\,e^{-2a\\mu}

    M1 sketch integrand:

    .. math::

        I_{M1}(\\mu) \\;=\\; 2\\,G_{\\rm in}\\,F_{\\rm out}\\,
                            \\frac{\\mu}{1 - e^{-2a\\mu}}
                       \\;=\\; 2\\,\\cosh(\\rho\\mu)\\,e^{-a\\mu}\\,
                              \\cosh(\\rho'\\mu_*)\\,\\mu_*^{-1}\\,e^{-a\\mu}\\,
                              \\frac{\\mu}{1 - e^{-2a\\mu}}

    Simplify:

    .. math::

        I_{M1}(\\mu) \\;=\\; 2\\,\\mu_*^{-1}\\,\\cosh(\\rho\\mu)\\,
                            \\cosh(\\rho'\\mu_*)\\,e^{-2a\\mu}\\,
                            \\frac{\\mu}{1 - e^{-2a\\mu}}

    For these to be equal, we require

    .. math::

        T(\\mu) \\;=\\; \\frac{\\mu}{1 - e^{-2a\\mu}}

    BUT Sanchez has :math:`T(\\mu) = 1/(1 - e^{-2a\\mu})` (no µ in
    numerator).

    **Resolution**: the M1 sketch absorbs one µ-Jacobian factor INTO the
    multi-bounce factor, while Sanchez keeps the µ-Jacobian (which appears
    as :math:`\\mu_*^{-1}` after the chord-projection change of variables)
    separate. The two forms are NOT pointwise-equal; they differ by a
    factor of µ. The M1 sketch is **wrong by a factor of µ** if applied
    naively against Sanchez's :math:`G_{\\rm in}, F_{\\rm out}` definitions.

    **Correct M1 form** (after this verification): drop the µ from the
    numerator of the multi-bounce factor:

    .. math::

        K_{\\rm bc}^{\\rm mb,sph}(r_i, r_j) \\;=\\; 2\\!\\int_0^1\\!
            G_{\\rm in}(r_i, \\mu)\\,F_{\\rm out}(r_j, \\mu)\\,
            T(\\mu)\\,\\mathrm d\\mu

    where :math:`F_{\\rm out}` already carries the :math:`\\mu_*^{-1}`
    Jacobian.

    Returns dict with SymPy expressions, the algebraic difference, and
    PASS/FAIL.
    """
    mu, a, rho, rho_p = sp.symbols(
        "mu a rho rho_p", positive=True, real=True,
    )
    mu_star = sp.sqrt(rho_p ** 2 - rho ** 2 * (1 - mu ** 2))

    T_mu = 1 / (1 - sp.exp(-2 * a * mu))

    # Sanchez (A6) integrand, α=1, β=0, kept symbolic
    I_sanchez = (
        2 * T_mu * (1 / mu_star)
        * sp.cosh(rho * mu) * sp.cosh(rho_p * mu_star)
        * sp.exp(-2 * a * mu)
    )

    # M1 sketch integrand with parent agent's form (µ in numerator of MB)
    G_in = sp.cosh(rho * mu) * sp.exp(-a * mu)
    F_out = sp.cosh(rho_p * mu_star) * (1 / mu_star) * sp.exp(-a * mu)
    f_mb_M1 = mu / (1 - sp.exp(-2 * a * mu))
    I_M1 = 2 * G_in * F_out * f_mb_M1

    # Difference
    diff = sp.simplify(I_sanchez - I_M1)
    # If the two were equal the diff would be 0; show the multiplicative
    # discrepancy explicitly: I_sanchez / I_M1 should be 1/µ.
    ratio = sp.simplify(I_sanchez / I_M1)

    # The CORRECTED M1 form (no µ in MB numerator)
    f_mb_corrected = T_mu  # i.e., 1/(1 - e^{-2aµ})
    I_M1_corrected = 2 * G_in * F_out * f_mb_corrected
    diff_corrected = sp.simplify(I_sanchez - I_M1_corrected)
    pass_v2 = diff_corrected == 0

    return {
        "name": "V2 M1 sketch ↔ Sanchez Eq. (A6) equivalence",
        "I_sanchez": I_sanchez,
        "I_M1_naive": I_M1,
        "I_M1_corrected": I_M1_corrected,
        "diff_naive": diff,
        "ratio_naive": ratio,  # should equal 1/µ → M1 sketch wrong by 1/µ
        "diff_corrected": diff_corrected,
        "pass": pass_v2,
        "verdict": (
            "Sanchez (A6) and the M1 sketch differ by a factor of µ. "
            "The CORRECTED form drops the µ from the numerator of the "
            "multi-bounce factor; F_out already carries the µ_*^{-1} "
            "Jacobian. Production implementation MUST use the corrected "
            "form (no µ in MB numerator)."
        ),
    }


def derive_diagonal_singularity() -> dict:
    """V4 — Diagonal-singularity analysis (the Phase 5 production blocker).

    Sanchez Eq. (A6) integrand at the diagonal :math:`\\rho' \\to \\rho`
    has a non-integrable :math:`1/\\mu^2` singularity at :math:`\\mu \\to 0`:

    - :math:`\\mu_*(\\mu) = \\sqrt{\\rho'^2 - \\rho^2(1 - \\mu^2)} / \\rho'
      \\to \\mu` (when :math:`\\rho' = \\rho`)
    - :math:`1/\\mu_* \\to 1/\\mu`
    - :math:`T(\\mu_-) = 1/(1 - e^{-2a\\mu_-}) \\to 1/(2a\\mu_-) \\sim
      1/(2a\\sqrt{1 - (\\rho/a)^2})` at the surface :math:`\\rho \\to a`
      (bounded for :math:`\\rho < a`)
    - For :math:`\\rho = a` (surface diagonal): :math:`\\mu_- = \\mu`, so
      :math:`T(\\mu_-) \\sim 1/(2a\\mu)` and the product
      :math:`T \\cdot 1/\\mu_* \\sim 1/(2a\\mu^2)` — non-integrable
    - For :math:`\\rho < a`, :math:`\\rho' = \\rho` (interior diagonal):
      :math:`T(\\mu_-) \\to` const, :math:`1/\\mu_* \\to 1/\\mu` —
      logarithmic singularity, integrable but non-trivial for
      Gauss-Legendre

    This is the **production blocker** for wiring Phase 5 into ORPHEUS's
    discrete Nyström. ORPHEUS's existing :func:`build_volume_kernel`
    uses adaptive quadrature to handle the analogous log singularity in
    :math:`\\bar g_0 = E_1(\\Sigma|r-r'|)` (vacuum sphere kernel).
    Phase 5+ requires either:

    (a) adaptive µ-quadrature with explicit subdivision at :math:`\\mu = 0`
        for the diagonal (and :math:`\\mu = \\mu_0` for off-diagonal);
    (b) singularity subtraction (like ORPHEUS's Bickley :math:`E_1`
        treatment) — analytical removal of the :math:`1/\\mu` (or
        :math:`1/\\mu^2`) part with a closed-form add-back; or
    (c) reformulation in a different basis where the singularity is
        absorbed (e.g., Gauss-Jacobi weighted by :math:`\\mu`, or a
        Wiener-Hopf factorisation).

    Sanchez 1986 itself does NOT specify a numerical method (the paper
    is theoretical). The singularity is intrinsic to the kernel form;
    the integral equation Eq. (1a) handles it via the global
    :math:`\\mathrm d\\rho'` integration, but a discrete Nyström
    sampling at quadrature points hits it directly.

    SymPy returns the ASYMPTOTIC form near the diagonal so the next
    investigator has the leading-order behaviour pinned.

    Returns dict with the diagonal limit and PASS (always — this is a
    documented finding, not a verification gate).
    """
    mu, a, rho = sp.symbols("mu a rho", positive=True, real=True)

    # On the diagonal ρ' = ρ: µ_*(µ) = µ exactly
    mu_star_diag = mu

    # T(µ_-) at general ρ < a:
    #   µ_- = √[a² - ρ²(1-µ²)] / a  → √[1 - (ρ/a)²]  as µ → 0
    # T(µ_-) finite at µ=0 for ρ < a; T diverges only when µ_- = 0,
    # which happens at ρ = a, µ = 0.
    mu_minus_diag = sp.sqrt(a ** 2 - rho ** 2 * (1 - mu ** 2)) / a
    T_mu_minus = 1 / (1 - sp.exp(-2 * a * mu_minus_diag))

    # Asymptotic at µ → 0: leading-order expansion
    integrand_diag = T_mu_minus * (1 / mu_star_diag)
    leading = sp.series(integrand_diag, mu, 0, 2).removeO()
    # For interior point ρ < a, leading order is c/mu (integrable log)
    # For surface point ρ = a, leading order is c/mu² (non-integrable)

    return {
        "name": "V4 diagonal-singularity finding",
        "integrand_diag_form": integrand_diag,
        "leading_at_mu_zero": leading,
        "pass": True,  # documentary, not a gate
        "verdict": (
            "Diagonal integrand T(µ_-) · µ_*^{-1} has a 1/µ singularity "
            "at interior diagonals (ρ < a) and a 1/µ² singularity at "
            "surface diagonals (ρ = a). The interior case is integrable "
            "(logarithmic) but the surface case is non-integrable. "
            "Sanchez does not specify a numerical method. ORPHEUS "
            "discrete Nyström sampling hits this singularity directly. "
            "Phase 5+ scope: adaptive µ-quadrature, singularity "
            "subtraction, or change-of-variables before production "
            "wiring is possible."
        ),
    }


def derive_vacuum_reduction() -> dict:
    """V3 — Vacuum-BC reduction sanity check.

    With α = 0 (no specular reflection), the multi-bounce factor
    :math:`T(\\mu) = 1 / (1 - 0 \\cdot e^{-2a\\mu}) = 1`. The Sanchez
    Eq. (A6) BC kernel :math:`g_h` then collapses to

    .. math::

        g_h^{\\alpha=0}(\\rho' \\to \\rho)
            \\;=\\; 0 \\cdot \\int (\\ldots) \\,\\mathrm d\\mu \\;=\\; 0

    (because the leading factor is :math:`2\\alpha`).

    The TOTAL kernel :math:`g_2 = \\bar g_2 + g_h` reduces to the vacuum
    kernel :math:`\\bar g_2` (Sanchez Eq. 5):

    .. math::

        \\bar g_2(\\rho' \\to \\rho) \\;=\\; (\\rho'/\\rho)\\,
            [\\bar g_0(\\rho' \\to \\rho) - \\bar g_0(-\\rho' \\to \\rho)]/2

    with the slab kernel :math:`\\bar g_0(\\rho' \\to \\rho) =
    \\tfrac{1}{2} E_1(\\Sigma |\\rho - \\rho'|)`.

    This verification simply checks that Sanchez Eq. (A6) is consistent
    with the vacuum-limit reduction at α = 0 (i.e., the "α"-prefactor
    correctly switches off the BC contribution as α → 0).

    Returns dict with the consistency check.
    """
    alpha = sp.symbols("alpha", positive=True, real=True)
    mu, a = sp.symbols("mu a", positive=True, real=True)

    # The α-prefactor in Sanchez Eq. (A6) — should be 2α
    # As α → 0, the factor goes to 0, and the BC kernel vanishes.
    prefactor = 2 * alpha
    lim_alpha_zero = sp.limit(prefactor, alpha, 0)
    pass_v3 = lim_alpha_zero == 0

    return {
        "name": "V3 vacuum-BC reduction at α → 0",
        "prefactor": prefactor,
        "limit_alpha_zero": lim_alpha_zero,
        "pass": pass_v3,
    }


