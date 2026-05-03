r"""SymPy derivations of the LA-13511 closed-form :math:`k_\infty` cases.

This module is the **algebra-of-record** for the infinite-medium
fission-eigenvalue formulae from Sood/Forster/Parsons LA-13511 (1999),
Appendix A (Eqs 18-32 + 72-76). Branch-1 SymPy proves that:

* The 1G closed form (Eq 20 with the leading factor :math:`c`)
  algebraically reduces to the simpler form (Eq 19) :math:`k_\infty
  = \nu\Sigma_f / (\Sigma_t - \Sigma_s)`. The :math:`c` factor cancels
  identically.
* The 2G general formula (Eq 28, derived here directly from
  :math:`\det(M)=0` of Eq 25) and the no-upscatter reduction (Eq 29)
  follow from the 2G matrix balance equation.
* The 2G flux-ratio formula (Eq 32) follows from adding Eq 23 and
  Eq 24 with :math:`\chi_1 + \chi_2 = 1`.
* The general multi-group formula (Eq 76)
  :math:`k = \nu\Sigma_f^T (\Sigma_t - \Sigma_s)^{-1} \chi \nu\Sigma_f \phi`
  follows from the matrix balance :math:`\Sigma_t \phi = \Sigma_s
  \phi + (1/k)\, \chi (\nu\Sigma_f^T \phi)`.
* The MG formula at :math:`G=1` reduces bit-for-bit to Eq 19.

A note on Sood's typo
=====================

The published Eq 28 contains a typo: when the no-upscatter limit
:math:`\Sigma_{21s} \to 0` is taken and the result simplified, it
does NOT reduce to the printed Eq 29. The typo is in which
:math:`\Sigma_g^{\rm rem}` factor multiplies which :math:`\chi_g`
in the numerator: as printed, the :math:`\chi_1` numerator has
:math:`\Sigma_1^{\rm rem}\,\nu_1\Sigma_{1f}` (and :math:`\chi_2`
has :math:`\Sigma_2^{\rm rem}\,\nu_2\Sigma_{2f}`); the correct form
that matches Eq 29 swaps these so :math:`\chi_1` carries
:math:`\Sigma_2^{\rm rem}` and :math:`\chi_2` carries
:math:`\Sigma_1^{\rm rem}`. The SymPy derivation here computes
:math:`\det(M)=0` from Eq 25 directly and proves the corrected form
algebraically, then verifies the printed Eq 29 against the
corrected general form. The numerical reference values in
LA-13511 (k_inf = 2.683767, phi_ratio = 0.675229) match Eq 29 and
the corrected Eq 28; they do NOT match Eq 28 as printed.

Sood's group convention
=======================

Sood numbers groups :math:`g=N` (fast) → :math:`g=1` (slow), the
reverse of typical nuclear-engineering convention. The SymPy
derivations here use **Sood's symbols verbatim** (so equations match
the report letter-for-letter), but Branch-2 numpy code in
:mod:`..multi_group.k_inf` uses ORPHEUS's :math:`g=0` (fast) →
:math:`g=N-1` (slow) convention. The conversion is purely a relabeling
— no algebra changes — so the SymPy results apply equally to either
side.

References
----------

* Sood, Forster, Parsons (1999), LA-13511, Appendix A.
* Sood, Forster, Parsons (2003), *Progress in Nuclear Energy* 42, 55
  — journal-published condensation.
"""
from __future__ import annotations

import sympy as sp


# ═══════════════════════════════════════════════════════════════════
# 1G Infinite Medium — Sood Eqs 18-20
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_1g_eq_19() -> dict:
    r"""V_fn1.1 — 1G infinite-medium :math:`k_\infty` from balance.

    Starting from Sood Eq 18 (the 1G integrated transport equation
    for an infinite isotropic homogeneous medium):

    .. math::

        \Sigma_t\,\phi = \Sigma_s\,\phi + \frac{\nu\Sigma_f}{k_\infty}\,\phi

    we factor :math:`\phi` and solve for :math:`k_\infty`:

    .. math::

        k_\infty = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s}

    which is Sood Eq 19. The flux :math:`\phi` cancels (constant
    everywhere in the infinite medium), confirming the eigenvalue is
    flux-shape independent — the 1G degeneracy that makes 1G
    eigenvalue claims insufficient as L1 verification of any operator
    other than absorption/production ratios.

    Returns
    -------
    dict
        Includes ``"pass"`` (bool), the symbolic balance equation,
        and the solved-for :math:`k_\infty`.
    """
    # Symbols: positive reals, generic (no algebraic prejudice).
    Sigma_t, Sigma_s, nu_Sigma_f, phi, k = sp.symbols(
        "Sigma_t Sigma_s nu_Sigma_f phi k", positive=True
    )

    # Sood Eq 18: balance equation.
    eq18 = sp.Eq(Sigma_t * phi, Sigma_s * phi + nu_Sigma_f * phi / k)

    # Solve for k.
    k_solutions = sp.solve(eq18, k)
    assert len(k_solutions) == 1, f"Expected single k root, got {k_solutions}"
    k_derived = sp.simplify(k_solutions[0])

    # Expected closed form (Sood Eq 19).
    k_eq19 = nu_Sigma_f / (Sigma_t - Sigma_s)

    diff = sp.simplify(k_derived - k_eq19)
    pass_eq19 = (diff == 0)

    return {
        "name": "V_fn1.1: 1G k_inf from balance equation reduces to Eq 19",
        "eq18": eq18,
        "k_derived": k_derived,
        "k_eq19": k_eq19,
        "diff": diff,
        "pass": pass_eq19,
    }


def derive_kinf_1g_eq_20_simplifies_to_eq_19() -> dict:
    r"""V_fn1.2 — Sood Eq 20 algebraically equals Eq 19.

    Sood states the same 1G result two ways:

    * Eq 19: :math:`k_\infty = \nu\Sigma_f / (\Sigma_t - \Sigma_s)`.
    * Eq 20: :math:`k_\infty = c \cdot \nu\Sigma_f \Sigma_t / [(\Sigma_t
      - \Sigma_s)(\Sigma_s + \nu\Sigma_f)]` where
      :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`.

    Substituting :math:`c` into Eq 20:

    .. math::

        k_\infty = \frac{(\Sigma_s + \nu\Sigma_f)}{\Sigma_t}
                   \cdot \frac{\nu\Sigma_f \Sigma_t}
                              {(\Sigma_t - \Sigma_s)(\Sigma_s + \nu\Sigma_f)}
                = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s}

    so the :math:`c` and :math:`\Sigma_t` factors cancel and Eq 20
    collapses to Eq 19 exactly.

    This identity is the trivial Branch-1 anchor for Case 1
    (PUa-1-0-IN). It also shows why SymPy is the right tool here: a
    by-hand reader needs to mentally cancel four factors; the
    symbolic engine does it mechanically and the test confirms the
    cancellation closes to zero.
    """
    Sigma_t, Sigma_s, nu_Sigma_f = sp.symbols(
        "Sigma_t Sigma_s nu_Sigma_f", positive=True
    )
    c = (Sigma_s + nu_Sigma_f) / Sigma_t

    k_eq19 = nu_Sigma_f / (Sigma_t - Sigma_s)
    k_eq20 = c * nu_Sigma_f * Sigma_t / (
        (Sigma_t - Sigma_s) * (Sigma_s + nu_Sigma_f)
    )

    diff = sp.simplify(k_eq20 - k_eq19)
    pass_id = (diff == 0)

    return {
        "name": "V_fn1.2: Eq 20 simplifies to Eq 19 (c factor cancels)",
        "c_definition": c,
        "k_eq19": k_eq19,
        "k_eq20": k_eq20,
        "diff": diff,
        "pass": pass_id,
    }


# ═══════════════════════════════════════════════════════════════════
# 2G Infinite Medium — Sood Eqs 21-32
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_2g_general_from_matrix() -> dict:
    r"""V_fn2.1 — 2G general :math:`k_\infty` from :math:`\det(M)=0`.

    Sood Eqs 23-24 rearrange the 2G balance equations Eq 21-22 into
    a 2x2 homogeneous linear system :math:`M(k_\infty)\,\vec\phi = 0`
    (Sood Eq 25). Critical fission balance requires :math:`\det(M) = 0`,
    which is a quadratic in :math:`k_\infty`. One root is :math:`k=0`
    (trivial); the other is the desired :math:`k_\infty`.

    Sood notation (g=2 fast, g=1 slow, :math:`\Sigma_{g}^{\rm rem}
    = \Sigma_g - \Sigma_{ggs}`):

    .. math::

        M = \begin{pmatrix}
              -(\Sigma_{21s} + \frac{\chi_2}{k}\nu_1\Sigma_{1f}) &
              \Sigma_{2}^{\rm rem} - \frac{\chi_2}{k}\nu_2\Sigma_{2f} \\[2pt]
              \Sigma_{1}^{\rm rem} - \frac{\chi_1}{k}\nu_1\Sigma_{1f} &
              -(\Sigma_{12s} + \frac{\chi_1}{k}\nu_2\Sigma_{2f})
            \end{pmatrix}

    This SymPy derivation:

    1. Sets up :math:`M` symbolically.
    2. Computes :math:`\det(M)`, factors out :math:`1/k^2`, and solves
       the resulting quadratic for :math:`k`.
    3. Discards the :math:`k=0` root.
    4. Returns the surviving root as the **derived general 2G
       formula**, which is the corrected Eq 28.

    The PASS flag verifies that the derived root, when restricted
    to :math:`\Sigma_{21s} = 0` (no upscatter), reduces to Eq 29
    (which IS printed correctly in LA-13511 and matches the
    numerical reference). This indirectly proves Eq 28 as printed
    has a typo.
    """
    # Sood symbols verbatim. positive=True for reals; we let k range
    # over the reals since at quadratic level we will pick the
    # non-zero root manually.
    (
        Sigma_1, Sigma_2, Sigma_11s, Sigma_22s, Sigma_12s, Sigma_21s,
        nu_1, Sigma_1f, nu_2, Sigma_2f, chi_1, chi_2,
    ) = sp.symbols(
        "Sigma_1 Sigma_2 Sigma_11s Sigma_22s Sigma_12s Sigma_21s "
        "nu_1 Sigma_1f nu_2 Sigma_2f chi_1 chi_2",
        positive=True,
    )
    k = sp.symbols("k", positive=True)

    Sigma_2rem = Sigma_2 - Sigma_22s
    Sigma_1rem = Sigma_1 - Sigma_11s

    # Sood Eq 25 — matrix M acting on (phi_1, phi_2)^T:
    #
    #   row 1 (from Eq 23, the phi_2 balance rearranged):
    #     [-(Σ_{21s} + (χ_2/k)·ν_1·Σ_{1f}),  Σ_2^rem - (χ_2/k)·ν_2·Σ_{2f}]
    #   row 2 (from Eq 24, the phi_1 balance rearranged):
    #     [Σ_1^rem - (χ_1/k)·ν_1·Σ_{1f},   -(Σ_{12s} + (χ_1/k)·ν_2·Σ_{2f})]
    #
    # Writing it from the LA-13511 PDF page 43 (Eq 25) verbatim:
    M = sp.Matrix([
        [-(Sigma_21s + chi_2 / k * nu_1 * Sigma_1f),
         Sigma_2rem - chi_2 / k * nu_2 * Sigma_2f],
        [Sigma_1rem - chi_1 / k * nu_1 * Sigma_1f,
         -(Sigma_12s + chi_1 / k * nu_2 * Sigma_2f)],
    ])

    detM = sp.expand(M.det())

    # detM is a rational function in k of the form (A + B/k + C/k^2)
    # where C is the k=0 root contribution. Multiply through by k^2
    # to land on a polynomial.
    detM_poly = sp.expand(detM * k**2)
    poly_in_k = sp.Poly(detM_poly, k)
    coeffs = poly_in_k.all_coeffs()  # leading coefficient first

    # Quadratic in k: ax^2 + bx + c = 0. We expect c = 0 (so k=0 is a
    # root) and the surviving root is k = -c'/a where c' is the linear
    # coefficient. Sympy will hand it back as one of two solve() roots.
    k_roots = sp.solve(detM_poly, k)
    # Filter out the trivial k=0 root.
    k_roots_nonzero = [r for r in k_roots if sp.simplify(r) != 0]
    assert len(k_roots_nonzero) == 1, (
        f"Expected exactly one non-trivial root, got {k_roots_nonzero}"
    )
    k_general_derived = sp.simplify(k_roots_nonzero[0])

    # Cross-check against the no-upscatter limit (Sigma_21s -> 0).
    # Expected: Sood Eq 29 (printed correctly).
    k_no_upscatter_derived = sp.simplify(
        k_general_derived.subs(Sigma_21s, 0)
    )

    # Sood Eq 29 (verbatim transcription):
    #   k = chi_1·nu_1·Sigma_1f / Sigma_1^rem
    #     + chi_2·[nu_1·Sigma_1f·Sigma_12s / (Sigma_1^rem·Sigma_2^rem)
    #               + nu_2·Sigma_2f / Sigma_2^rem]
    k_eq29 = (
        chi_1 * nu_1 * Sigma_1f / Sigma_1rem
        + chi_2 * (
            nu_1 * Sigma_1f * Sigma_12s / (Sigma_1rem * Sigma_2rem)
            + nu_2 * Sigma_2f / Sigma_2rem
        )
    )
    k_eq29_simplified = sp.simplify(k_eq29)

    diff_29 = sp.simplify(k_no_upscatter_derived - k_eq29_simplified)
    pass_eq29_match = (diff_29 == 0)

    # ALSO check Eq 28 as printed in the PDF — which we expect to FAIL
    # (i.e. NOT match k_general_derived). This documents the typo
    # algebraically rather than just by numerical example.
    # Eq 28 as printed:
    #   k = [chi_1·(nu_2·Sigma_2f·Sigma_21s + Sigma_1^rem·nu_1·Sigma_1f)
    #        + chi_2·(nu_1·Sigma_1f·Sigma_12s + Sigma_2^rem·nu_2·Sigma_2f)]
    #     / [Sigma_1^rem·Sigma_2^rem - Sigma_12s·Sigma_21s]
    k_eq28_as_printed = (
        chi_1 * (nu_2 * Sigma_2f * Sigma_21s + Sigma_1rem * nu_1 * Sigma_1f)
        + chi_2 * (nu_1 * Sigma_1f * Sigma_12s + Sigma_2rem * nu_2 * Sigma_2f)
    ) / (Sigma_1rem * Sigma_2rem - Sigma_12s * Sigma_21s)

    # The corrected Eq 28 (what the SymPy derivation should produce):
    k_eq28_corrected = (
        chi_1 * (nu_2 * Sigma_2f * Sigma_21s + Sigma_2rem * nu_1 * Sigma_1f)
        + chi_2 * (nu_1 * Sigma_1f * Sigma_12s + Sigma_1rem * nu_2 * Sigma_2f)
    ) / (Sigma_1rem * Sigma_2rem - Sigma_12s * Sigma_21s)

    diff_28_printed = sp.simplify(k_general_derived - k_eq28_as_printed)
    diff_28_corrected = sp.simplify(k_general_derived - k_eq28_corrected)

    pass_eq28_typo_confirmed = (diff_28_printed != 0)
    pass_eq28_corrected = (diff_28_corrected == 0)

    pass_overall = bool(
        pass_eq29_match and pass_eq28_typo_confirmed and pass_eq28_corrected
    )

    return {
        "name": "V_fn2.1: 2G general k_inf derived from det(M)=0; "
                "Eq 29 verified, Eq 28 typo confirmed",
        "M": M,
        "k_general_derived": k_general_derived,
        "k_no_upscatter_derived": k_no_upscatter_derived,
        "k_eq29": k_eq29_simplified,
        "diff_eq29": diff_29,
        "pass_eq29_match": pass_eq29_match,
        "k_eq28_as_printed": k_eq28_as_printed,
        "k_eq28_corrected": k_eq28_corrected,
        "diff_eq28_printed": diff_28_printed,
        "diff_eq28_corrected": diff_28_corrected,
        "pass_eq28_typo_confirmed": pass_eq28_typo_confirmed,
        "pass_eq28_corrected": pass_eq28_corrected,
        "pass": pass_overall,
    }


def derive_kinf_2g_no_upscatter() -> dict:
    r"""V_fn2.2 — Eq 29 closed form solves :math:`\det(M_{\Sigma_{21s}=0})=0`.

    Standalone Branch-1 verification of the no-upscatter 2G formula
    (Sood Eq 29):

    .. math::

        k_\infty = \frac{\chi_1 \nu_1 \Sigma_{1f}}{\Sigma_1^{\rm rem}}
                 + \chi_2 \left[
                     \frac{\nu_1\Sigma_{1f}\Sigma_{12s}}
                          {\Sigma_1^{\rm rem}\Sigma_2^{\rm rem}}
                     + \frac{\nu_2\Sigma_{2f}}{\Sigma_2^{\rm rem}}
                   \right]

    Substituting Eq 29 into the 2G balance system (Sood Eqs 21-22 with
    :math:`\Sigma_{21s} = 0`) should make the system consistent (LHS
    of Eq 25 with :math:`\Sigma_{21s}=0` becomes zero).

    Equivalently: substitute :math:`k = k_{\rm Eq 29}` into
    :math:`\det(M(\Sigma_{21s}=0))` and verify it simplifies to zero.

    This is independent of V_fn2.1 (which derives the formula from
    scratch) — V_fn2.2 verifies the published formula directly without
    going through the quadratic-root machinery.
    """
    (
        Sigma_1, Sigma_2, Sigma_11s, Sigma_22s, Sigma_12s,
        nu_1, Sigma_1f, nu_2, Sigma_2f, chi_1, chi_2,
    ) = sp.symbols(
        "Sigma_1 Sigma_2 Sigma_11s Sigma_22s Sigma_12s "
        "nu_1 Sigma_1f nu_2 Sigma_2f chi_1 chi_2",
        positive=True,
    )

    Sigma_2rem = Sigma_2 - Sigma_22s
    Sigma_1rem = Sigma_1 - Sigma_11s

    # Sood Eq 29 (no upscatter):
    k_eq29 = (
        chi_1 * nu_1 * Sigma_1f / Sigma_1rem
        + chi_2 * (
            nu_1 * Sigma_1f * Sigma_12s / (Sigma_1rem * Sigma_2rem)
            + nu_2 * Sigma_2f / Sigma_2rem
        )
    )

    # Build M with Sigma_21s = 0:
    #   row 1: [-(0 + (χ_2/k)·ν_1·Σ_{1f}),  Σ_2^rem - (χ_2/k)·ν_2·Σ_{2f}]
    #   row 2: [Σ_1^rem - (χ_1/k)·ν_1·Σ_{1f},  -(Σ_{12s} + (χ_1/k)·ν_2·Σ_{2f})]
    k_sym = sp.symbols("k_sym", positive=True)
    M = sp.Matrix([
        [-chi_2 / k_sym * nu_1 * Sigma_1f,
         Sigma_2rem - chi_2 / k_sym * nu_2 * Sigma_2f],
        [Sigma_1rem - chi_1 / k_sym * nu_1 * Sigma_1f,
         -(Sigma_12s + chi_1 / k_sym * nu_2 * Sigma_2f)],
    ])

    detM = M.det()
    detM_at_eq29 = sp.simplify(detM.subs(k_sym, k_eq29))

    pass_eq29_zero = (detM_at_eq29 == 0)

    return {
        "name": "V_fn2.2: Sood Eq 29 makes det(M)=0 at no-upscatter limit",
        "k_eq29": sp.simplify(k_eq29),
        "detM_at_eq29": detM_at_eq29,
        "pass": pass_eq29_zero,
    }


def derive_phi_ratio_2g_no_upscatter() -> dict:
    r"""V_fn2.3 — Sood Eq 32 :math:`\phi_2/\phi_1` from balance + chi-sum.

    Adding Sood Eq 23 + Eq 24 with the constraint :math:`\chi_1 +
    \chi_2 = 1` eliminates the :math:`\chi_g` from the resulting
    relation (Sood Eq 30):

    .. math::

        \left[\Sigma_1^{\rm rem} - \Sigma_{21s} - \frac{\nu_1\Sigma_{1f}}{k_\infty}\right] \phi_1
      + \left[\Sigma_2^{\rm rem} - \Sigma_{12s} - \frac{\nu_2\Sigma_{2f}}{k_\infty}\right] \phi_2
        = 0

    For the no-upscatter case (:math:`\Sigma_{21s} = 0`) this becomes
    Sood Eq 32:

    .. math::

        \frac{\phi_2}{\phi_1}
        = \frac{\Sigma_1^{\rm rem} - \nu_1\Sigma_{1f}/k_\infty}
               {\nu_2\Sigma_{2f}/k_\infty - \Sigma_2^{\rm rem} + \Sigma_{12s}}

    This SymPy derivation:

    1. Builds the chi-sum-equals-one relation symbolically by adding
       Sood Eqs 23 + 24.
    2. Solves the resulting linear equation for :math:`\phi_2/\phi_1`.
    3. Compares to the printed Eq 32.
    """
    (
        Sigma_1, Sigma_2, Sigma_11s, Sigma_22s, Sigma_12s, Sigma_21s,
        nu_1, Sigma_1f, nu_2, Sigma_2f, k,
    ) = sp.symbols(
        "Sigma_1 Sigma_2 Sigma_11s Sigma_22s Sigma_12s Sigma_21s "
        "nu_1 Sigma_1f nu_2 Sigma_2f k",
        positive=True,
    )
    chi_1, chi_2, phi_1, phi_2 = sp.symbols("chi_1 chi_2 phi_1 phi_2", positive=True)

    Sigma_2rem = Sigma_2 - Sigma_22s
    Sigma_1rem = Sigma_1 - Sigma_11s

    # Sood Eq 23: [Σ_2 - Σ_22s - (χ_2/k)·ν_2·Σ_2f]·φ_2
    #             - [Σ_21s + (χ_2/k)·ν_1·Σ_1f]·φ_1 = 0
    eq23_lhs = (
        (Sigma_2rem - chi_2 / k * nu_2 * Sigma_2f) * phi_2
        - (Sigma_21s + chi_2 / k * nu_1 * Sigma_1f) * phi_1
    )

    # Sood Eq 24: [Σ_1 - Σ_11s - (χ_1/k)·ν_1·Σ_1f]·φ_1
    #             - [Σ_12s + (χ_1/k)·ν_2·Σ_2f]·φ_2 = 0
    eq24_lhs = (
        (Sigma_1rem - chi_1 / k * nu_1 * Sigma_1f) * phi_1
        - (Sigma_12s + chi_1 / k * nu_2 * Sigma_2f) * phi_2
    )

    # Sum + apply χ_1 + χ_2 = 1.
    sum_eq = sp.expand(eq23_lhs + eq24_lhs)
    sum_eq_chi_eliminated = sp.simplify(sum_eq.subs(chi_2, 1 - chi_1))

    # The chi_1 dependence should drop out.
    chi_1_coeff = sp.simplify(sum_eq_chi_eliminated.coeff(chi_1))
    pass_chi_eliminates = (sp.simplify(chi_1_coeff) == 0)

    # The remaining equation is Sood Eq 30:
    eq30 = sp.simplify(sum_eq_chi_eliminated.subs(chi_1, 0))  # chi-free residue

    # Now restrict to no-upscatter (Σ_21s = 0).
    eq30_nou = eq30.subs(Sigma_21s, 0)

    # Solve for the ratio φ_2/φ_1.
    ratio = sp.symbols("ratio", positive=True)
    eq30_nou_in_ratio = sp.simplify(
        eq30_nou.subs(phi_2, ratio * phi_1) / phi_1
    )
    ratio_solutions = sp.solve(eq30_nou_in_ratio, ratio)
    assert len(ratio_solutions) == 1, (
        f"Expected single ratio root, got {ratio_solutions}"
    )
    ratio_derived = sp.simplify(ratio_solutions[0])

    # Sood Eq 32:
    ratio_eq32 = (
        (Sigma_1rem - nu_1 * Sigma_1f / k)
        / (nu_2 * Sigma_2f / k - Sigma_2rem + Sigma_12s)
    )
    ratio_eq32_simplified = sp.simplify(ratio_eq32)

    diff = sp.simplify(ratio_derived - ratio_eq32_simplified)
    pass_eq32 = (diff == 0)

    return {
        "name": "V_fn2.3: phi_2/phi_1 derivation matches Sood Eq 32",
        "ratio_derived": ratio_derived,
        "ratio_eq32": ratio_eq32_simplified,
        "diff": diff,
        "pass_chi_eliminates": pass_chi_eliminates,
        "pass": pass_eq32 and pass_chi_eliminates,
    }


# ═══════════════════════════════════════════════════════════════════
# General multi-group infinite medium — Sood Eqs 72-76
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_mg_matrix_form() -> dict:
    r"""V_fnMG.1 — General MG balance reduces to single matrix inversion.

    Sood Eq 72 is the matrix balance for an infinite medium:

    .. math::

        \overline{\overline{\Sigma_t}}\,\bar\phi
        = \overline{\overline{\Sigma_s}}\,\bar\phi
        + \frac{1}{k}\, \bar\chi\, \overline{\nu\Sigma_f}\,\bar\phi

    Rearranging (Eq 73):

    .. math::

        (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})\,\bar\phi
        = \frac{1}{k}\, \bar\chi\, (\overline{\nu\Sigma_f}\,\bar\phi)

    Inverting and projecting onto the production vector
    :math:`\overline{\nu\Sigma_f}` (Eqs 74-76):

    .. math::

        k = \overline{\nu\Sigma_f}\,
            (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})^{-1}\,
            \bar\chi

    is a scalar (single matrix inversion).

    This SymPy verification derives this for **G=2** symbolically as
    proof of concept, and confirms the result equals the dominant
    eigenvalue of :math:`A^{-1}\,\chi\,(\nu\Sigma_f)^T` (the form
    used by :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`).

    For G ≥ 5, symbolic eigenvalue closed forms break (Abel-Ruffini);
    Sood Eq 76 is the right form for the G=2 case where the formula
    *is* a clean closed form. We verify Eq 76 symbolically for G=2;
    the algebraic structure is identical for general G but the
    eigenvalue is no longer a closed-form expression in the cross
    sections for G ≥ 5.

    Note on convention
    ------------------

    Sood writes the scattering operator as :math:`\Sigma_s\phi`, which
    in matrix form is :math:`(\Sigma_s)_{gg'}` = scattering FROM g'
    TO g. ORPHEUS stores ``sigma_s[g, h]`` = scattering FROM g TO h,
    which is :math:`(\Sigma_s^T)_{gg'}`. The two forms are related by
    transpose and the dominant eigenvalue is invariant under this
    transposition since :math:`A^{-1}F` and :math:`F^T A^{-T}` have
    the same spectrum. This subtlety is consequential at the
    flux-eigenvector level; for k_inf alone the convention is invisible.
    """
    # G=2 case symbolic.
    # Use generic positive symbols (no upscatter restriction, since
    # the matrix form handles upscatter naturally).
    Sigma_t1, Sigma_t2 = sp.symbols("Sigma_t1 Sigma_t2", positive=True)
    Sigma_s11, Sigma_s12, Sigma_s21, Sigma_s22 = sp.symbols(
        "Sigma_s11 Sigma_s12 Sigma_s21 Sigma_s22", nonnegative=True
    )
    nuSf1, nuSf2, chi1, chi2 = sp.symbols(
        "nuSf1 nuSf2 chi1 chi2", positive=True
    )
    phi1, phi2, k = sp.symbols("phi1 phi2 k", positive=True)

    # Sood-style: Σ_s acts as (Σ_s_ij = scattering TO i FROM j).
    # i.e. the scattering source into group i is sum_j Σ_s_ij·φ_j.
    Sigma_t = sp.Matrix([[Sigma_t1, 0], [0, Sigma_t2]])
    Sigma_s = sp.Matrix([
        [Sigma_s11, Sigma_s12],
        [Sigma_s21, Sigma_s22],
    ])
    chi = sp.Matrix([[chi1], [chi2]])
    nuSf = sp.Matrix([[nuSf1, nuSf2]])  # row vector
    phi = sp.Matrix([[phi1], [phi2]])

    # Sood Eq 72: Σ_t·φ = Σ_s·φ + (1/k)·χ·(νΣ_f·φ)
    A = Sigma_t - Sigma_s
    fission_source = chi * (nuSf * phi)  # (2x1) · (1x1) = 2x1
    eq72 = sp.Eq(A * phi, fission_source / k)

    # Sood Eq 76: k = νΣ_f · A^{-1} · χ
    A_inv = A.inv()
    k_eq76_matrix = nuSf * A_inv * chi  # (1x2)·(2x2)·(2x1) = (1x1) scalar
    k_eq76 = sp.simplify(k_eq76_matrix[0, 0])

    # Cross-check: dominant eigenvalue of M = A^{-1} · χ · νΣ_f
    M = A_inv * chi * nuSf  # (2x2) outer-product matrix, rank 1
    eigvals = list(M.eigenvals().keys())
    # Rank-1 matrix has one zero eigenvalue and one nonzero; trace = nonzero ev.
    M_trace = sp.simplify(sp.trace(M))
    # Trace of a rank-1 matrix outer(u, v) is dot(v, u). Same as Eq 76.
    diff = sp.simplify(M_trace - k_eq76)
    pass_trace = (diff == 0)

    # Also verify rank-1 structure (one zero eigenvalue).
    nonzero_evs = [e for e in eigvals if sp.simplify(e) != 0]
    pass_rank1 = (len(nonzero_evs) == 1)

    # Verify the surviving eigenvalue equals Eq 76.
    if pass_rank1:
        diff_ev = sp.simplify(nonzero_evs[0] - k_eq76)
        pass_ev = (diff_ev == 0)
    else:
        diff_ev = None
        pass_ev = False

    return {
        "name": "V_fnMG.1: Sood Eq 76 for G=2 — k = nuSf · A^{-1} · chi",
        "A": A,
        "k_eq76": k_eq76,
        "M_trace": M_trace,
        "diff_trace": diff,
        "pass_trace_equals_eq76": pass_trace,
        "pass_M_is_rank1": pass_rank1,
        "pass_eigenvalue_equals_eq76": pass_ev,
        "pass": pass_trace and pass_rank1 and pass_ev,
    }


def derive_kinf_mg_reduces_to_1g() -> dict:
    r"""V_fnMG.2 — MG formula with G=1 reduces to 1G Eq 19 bit-equal.

    The general MG formula Sood Eq 76:

    .. math::

        k = \overline{\nu\Sigma_f}\,
            (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})^{-1}\,
            \bar\chi

    at G=1 collapses to scalars: :math:`A` is the :math:`1\times 1`
    matrix :math:`(\Sigma_t - \Sigma_s)`, :math:`A^{-1} = 1/(\Sigma_t -
    \Sigma_s)`, :math:`\chi = (1)`, :math:`\nu\Sigma_f` is a scalar.
    Therefore :math:`k_{\rm Eq 76} = \nu\Sigma_f / (\Sigma_t -
    \Sigma_s) = k_{\rm Eq 19}`.

    This is the trivial dimensional-reduction check: the MG
    infrastructure must reproduce the 1G result exactly when only
    one group is present.
    """
    Sigma_t, Sigma_s, nu_Sigma_f = sp.symbols(
        "Sigma_t Sigma_s nu_Sigma_f", positive=True
    )

    # G=1 case as 1x1 matrices.
    A = sp.Matrix([[Sigma_t - Sigma_s]])
    chi = sp.Matrix([[1]])
    nuSf = sp.Matrix([[nu_Sigma_f]])
    k_eq76_g1 = sp.simplify((nuSf * A.inv() * chi)[0, 0])

    k_eq19 = nu_Sigma_f / (Sigma_t - Sigma_s)

    diff = sp.simplify(k_eq76_g1 - k_eq19)
    pass_id = (diff == 0)

    return {
        "name": "V_fnMG.2: Eq 76 with G=1 reduces to Eq 19",
        "k_eq76_g1": k_eq76_g1,
        "k_eq19": k_eq19,
        "diff": diff,
        "pass": pass_id,
    }
