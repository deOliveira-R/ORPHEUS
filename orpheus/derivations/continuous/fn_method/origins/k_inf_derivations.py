r"""SymPy derivations of Sood's closed-form :math:`k_\infty` cases.

This module is the **algebra-of-record** for the infinite-medium
fission-eigenvalue formulae of Sood, Forster & Parsons (2003),
Appendix A (Eqs (A.1)-(A.15) + (A.55)-(A.59)). Branch-1 SymPy proves that:

* The 1G closed form (Eq (A.3) with the leading factor :math:`c`)
  algebraically reduces to the simpler form (Eq (A.2)) :math:`k_\infty
  = \nu\Sigma_f / (\Sigma_t - \Sigma_s)`. The :math:`c` factor cancels
  identically.
* The 2G general formula (Eq (A.11)) is the non-trivial root of
  :math:`\det(M)=0` of Eq (A.8), as printed, and its no-upscatter
  reduction is Eq (A.12).
* The 2G flux-ratio formula (Eq (A.15)) follows from adding Eq (A.6) and
  Eq (A.7) with :math:`\chi_1 + \chi_2 = 1`.
* The general multi-group formula (Eq (A.59))
  :math:`k = \nu\Sigma_f^T (\Sigma_t - \Sigma_s)^{-1} \chi \nu\Sigma_f \phi`
  follows from the matrix balance :math:`\Sigma_t \phi = \Sigma_s
  \phi + (1/k)\, \chi (\nu\Sigma_f^T \phi)`.
* The MG formula at :math:`G=1` reduces bit-for-bit to Eq (A.2).

Eq (A.11) is printed correctly (a refuted typo claim)
=====================================================

Until 2026-09-29 this module asserted that the published Eq 28 (the
1999 numbering of Eq (A.11)) has a typo, and its gate certified the
claim. The claim came from a
mis-transcription: the "as printed" form typed here paired
:math:`\chi_1` with :math:`\Sigma_1^{\rm rem}\,\nu_1\Sigma_{1f}` and
:math:`\chi_2` with :math:`\Sigma_2^{\rm rem}\,\nu_2\Sigma_{2f}`,
while the report prints :math:`\chi_1` with
:math:`\Sigma_2^{\rm rem}\,\nu_1\Sigma_{1f}` and :math:`\chi_2` with
:math:`\Sigma_1^{\rm rem}\,\nu_2\Sigma_{2f}` (the 1999 page image and
the text of both editions, which print the same body). The printed
Eq (A.11) is the non-trivial root of
:math:`\det(M)=0`, and on PU-2-0-IN it gives the published
k_inf = 2.683767. The mis-transcribed form gives 2.862. The gate that
certified the typo compared the derivation with the mis-transcription,
so it would have passed for any wrong transcription (ERR-092).

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

* Sood, A., Forster, R.A., Parsons, D.K. (2003), "Analytical benchmark
  test set for criticality code verification", *Progress in Nuclear
  Energy* 42(1), 55-106, Appendix A: the edition cited here.
* The same authors' 1999 report LA-13511 numbers these equations
  18-32 and 72-76: its appendix Eq N is Eq (A.N-17) here.
"""
from __future__ import annotations

import sympy as sp


# ═══════════════════════════════════════════════════════════════════
# 1G Infinite Medium — Sood Eqs (A.1)-(A.3)
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_1g_eq_a2() -> dict:
    r"""V_fn1.1 — 1G infinite-medium :math:`k_\infty` from balance.

    Starting from Sood Eq (A.1) (the 1G integrated transport equation
    for an infinite isotropic homogeneous medium):

    .. math::

        \Sigma_t\,\phi = \Sigma_s\,\phi + \frac{\nu\Sigma_f}{k_\infty}\,\phi

    we factor :math:`\phi` and solve for :math:`k_\infty`:

    .. math::

        k_\infty = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s}

    which is Sood Eq (A.2). The flux :math:`\phi` cancels (constant
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

    # Sood Eq (A.1): balance equation.
    eq_a1 = sp.Eq(Sigma_t * phi, Sigma_s * phi + nu_Sigma_f * phi / k)

    # Solve for k.
    k_solutions = sp.solve(eq_a1, k)
    assert len(k_solutions) == 1, f"Expected single k root, got {k_solutions}"
    k_derived = sp.simplify(k_solutions[0])

    # Expected closed form (Sood Eq (A.2)).
    k_eq_a2 = nu_Sigma_f / (Sigma_t - Sigma_s)

    diff = sp.simplify(k_derived - k_eq_a2)
    pass_eq_a2 = (diff == 0)

    return {
        "name": "V_fn1.1: 1G k_inf from balance equation reduces to Eq (A.2)",
        "eq_a1": eq_a1,
        "k_derived": k_derived,
        "k_eq_a2": k_eq_a2,
        "diff": diff,
        "pass": pass_eq_a2,
    }


def derive_kinf_1g_eq_a3_simplifies_to_eq_a2() -> dict:
    r"""V_fn1.2 — Sood Eq (A.3) algebraically equals Eq (A.2).

    Sood states the same 1G result two ways:

    * Eq (A.2): :math:`k_\infty = \nu\Sigma_f / (\Sigma_t - \Sigma_s)`.
    * Eq (A.3): :math:`k_\infty = c \cdot \nu\Sigma_f \Sigma_t / [(\Sigma_t
      - \Sigma_s)(\Sigma_s + \nu\Sigma_f)]` where
      :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`.

    Substituting :math:`c` into Eq (A.3):

    .. math::

        k_\infty = \frac{(\Sigma_s + \nu\Sigma_f)}{\Sigma_t}
                   \cdot \frac{\nu\Sigma_f \Sigma_t}
                              {(\Sigma_t - \Sigma_s)(\Sigma_s + \nu\Sigma_f)}
                = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s}

    so the :math:`c` and :math:`\Sigma_t` factors cancel and Eq (A.3)
    collapses to Eq (A.2) exactly.

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

    k_eq_a2 = nu_Sigma_f / (Sigma_t - Sigma_s)
    k_eq_a3 = c * nu_Sigma_f * Sigma_t / (
        (Sigma_t - Sigma_s) * (Sigma_s + nu_Sigma_f)
    )

    diff = sp.simplify(k_eq_a3 - k_eq_a2)
    pass_id = (diff == 0)

    return {
        "name": "V_fn1.2: Eq (A.3) simplifies to Eq (A.2) (c factor cancels)",
        "c_definition": c,
        "k_eq_a2": k_eq_a2,
        "k_eq_a3": k_eq_a3,
        "diff": diff,
        "pass": pass_id,
    }


# ═══════════════════════════════════════════════════════════════════
# 2G Infinite Medium — Sood Eqs (A.4)-(A.15)
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_2g_general_from_matrix() -> dict:
    r"""V_fn2.1 — 2G general :math:`k_\infty` from :math:`\det(M)=0`.

    Sood Eqs (A.6)-(A.7) rearrange the 2G balance equations Eq (A.4)-(A.5) into
    a 2x2 homogeneous linear system :math:`M(k_\infty)\,\vec\phi = 0`
    (Sood Eq (A.8)). Critical fission balance requires :math:`\det(M) = 0`,
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
       formula**.

    The PASS flag verifies that the derived root equals Eq (A.11) as
    printed, and that its restriction to :math:`\Sigma_{21s} = 0`
    (no upscatter) is Eq (A.12). An Eq (A.11) transcription that swaps the
    :math:`\Sigma_g^{\rm rem}` factors between the two
    :math:`\chi_g` terms fails the first check (ERR-092).
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

    # Sood Eq (A.8) — matrix M acting on (phi_1, phi_2)^T:
    #
    #   row 1 (from Eq (A.6), the phi_2 balance rearranged):
    #     [-(Σ_{21s} + (χ_2/k)·ν_1·Σ_{1f}),  Σ_2^rem - (χ_2/k)·ν_2·Σ_{2f}]
    #   row 2 (from Eq (A.7), the phi_1 balance rearranged):
    #     [Σ_1^rem - (χ_1/k)·ν_1·Σ_{1f},   -(Σ_{12s} + (χ_1/k)·ν_2·Σ_{2f})]
    #
    # Writing it from Sood 2003 p. 100 (Eq (A.8)) verbatim:
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
    # Expected: Sood Eq (A.12) (printed correctly).
    k_no_upscatter_derived = sp.simplify(
        k_general_derived.subs(Sigma_21s, 0)
    )

    # Sood Eq (A.12) (verbatim transcription):
    #   k = chi_1·nu_1·Sigma_1f / Sigma_1^rem
    #     + chi_2·[nu_1·Sigma_1f·Sigma_12s / (Sigma_1^rem·Sigma_2^rem)
    #               + nu_2·Sigma_2f / Sigma_2^rem]
    k_eq_a12 = (
        chi_1 * nu_1 * Sigma_1f / Sigma_1rem
        + chi_2 * (
            nu_1 * Sigma_1f * Sigma_12s / (Sigma_1rem * Sigma_2rem)
            + nu_2 * Sigma_2f / Sigma_2rem
        )
    )
    k_eq_a12_simplified = sp.simplify(k_eq_a12)

    diff_a12 = sp.simplify(k_no_upscatter_derived - k_eq_a12_simplified)
    pass_eq_a12_match = (diff_a12 == 0)

    # Sood Eq (A.11) (verbatim transcription; read on the 1999 page
    # image, where it is Eq 28):
    #   k = [chi_1·(nu_2·Sigma_2f·Sigma_21s + Sigma_2^rem·nu_1·Sigma_1f)
    #        + chi_2·(nu_1·Sigma_1f·Sigma_12s + Sigma_1^rem·nu_2·Sigma_2f)]
    #     / [Sigma_1^rem·Sigma_2^rem - Sigma_12s·Sigma_21s]
    k_eq_a11 = (
        chi_1 * (nu_2 * Sigma_2f * Sigma_21s + Sigma_2rem * nu_1 * Sigma_1f)
        + chi_2 * (nu_1 * Sigma_1f * Sigma_12s + Sigma_1rem * nu_2 * Sigma_2f)
    ) / (Sigma_1rem * Sigma_2rem - Sigma_12s * Sigma_21s)

    diff_a11 = sp.simplify(k_general_derived - k_eq_a11)
    pass_eq_a11_match = (diff_a11 == 0)

    pass_overall = bool(pass_eq_a12_match and pass_eq_a11_match)

    return {
        "name": "V_fn2.1: 2G general k_inf derived from det(M)=0; "
                "Eq (A.11) and Eq (A.12) verified",
        "M": M,
        "k_general_derived": k_general_derived,
        "k_no_upscatter_derived": k_no_upscatter_derived,
        "k_eq_a12": k_eq_a12_simplified,
        "diff_eq_a12": diff_a12,
        "pass_eq_a12_match": pass_eq_a12_match,
        "k_eq_a11": k_eq_a11,
        "diff_eq_a11": diff_a11,
        "pass_eq_a11_match": pass_eq_a11_match,
        "pass": pass_overall,
    }


def derive_kinf_2g_no_upscatter() -> dict:
    r"""V_fn2.2 — Eq (A.12) closed form solves :math:`\det(M_{\Sigma_{21s}=0})=0`.

    Standalone Branch-1 verification of the no-upscatter 2G formula
    (Sood Eq (A.12)):

    .. math::

        k_\infty = \frac{\chi_1 \nu_1 \Sigma_{1f}}{\Sigma_1^{\rm rem}}
                 + \chi_2 \left[
                     \frac{\nu_1\Sigma_{1f}\Sigma_{12s}}
                          {\Sigma_1^{\rm rem}\Sigma_2^{\rm rem}}
                     + \frac{\nu_2\Sigma_{2f}}{\Sigma_2^{\rm rem}}
                   \right]

    Substituting Eq (A.12) into the 2G balance system (Sood Eqs (A.4)-(A.5) with
    :math:`\Sigma_{21s} = 0`) should make the system consistent (LHS
    of Eq (A.8) with :math:`\Sigma_{21s}=0` becomes zero).

    Equivalently: substitute :math:`k = k_{\rm Eq (A.12)}` into
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

    # Sood Eq (A.12) (no upscatter):
    k_eq_a12 = (
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
    detM_at_eq_a12 = sp.simplify(detM.subs(k_sym, k_eq_a12))

    pass_eq_a12_zero = (detM_at_eq_a12 == 0)

    return {
        "name": "V_fn2.2: Sood Eq (A.12) makes det(M)=0 at no-upscatter limit",
        "k_eq_a12": sp.simplify(k_eq_a12),
        "detM_at_eq_a12": detM_at_eq_a12,
        "pass": pass_eq_a12_zero,
    }


def derive_phi_ratio_2g_no_upscatter() -> dict:
    r"""V_fn2.3 — Sood Eq (A.15) :math:`\phi_2/\phi_1` from balance + chi-sum.

    Adding Sood Eq (A.6) + Eq (A.7) with the constraint :math:`\chi_1 +
    \chi_2 = 1` eliminates the :math:`\chi_g` from the resulting
    relation (Sood Eq (A.13)):

    .. math::

        \left[\Sigma_1^{\rm rem} - \Sigma_{21s} - \frac{\nu_1\Sigma_{1f}}{k_\infty}\right] \phi_1
      + \left[\Sigma_2^{\rm rem} - \Sigma_{12s} - \frac{\nu_2\Sigma_{2f}}{k_\infty}\right] \phi_2
        = 0

    For the no-upscatter case (:math:`\Sigma_{21s} = 0`) this becomes
    Sood Eq (A.15):

    .. math::

        \frac{\phi_2}{\phi_1}
        = \frac{\Sigma_1^{\rm rem} - \nu_1\Sigma_{1f}/k_\infty}
               {\nu_2\Sigma_{2f}/k_\infty - \Sigma_2^{\rm rem} + \Sigma_{12s}}

    This SymPy derivation:

    1. Builds the chi-sum-equals-one relation symbolically by adding
       Sood Eqs (A.6) + (A.7).
    2. Solves the resulting linear equation for :math:`\phi_2/\phi_1`.
    3. Compares to the printed Eq (A.15).
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

    # Sood Eq (A.6): [Σ_2 - Σ_22s - (χ_2/k)·ν_2·Σ_2f]·φ_2
    #             - [Σ_21s + (χ_2/k)·ν_1·Σ_1f]·φ_1 = 0
    eq_a6_lhs = (
        (Sigma_2rem - chi_2 / k * nu_2 * Sigma_2f) * phi_2
        - (Sigma_21s + chi_2 / k * nu_1 * Sigma_1f) * phi_1
    )

    # Sood Eq (A.7): [Σ_1 - Σ_11s - (χ_1/k)·ν_1·Σ_1f]·φ_1
    #             - [Σ_12s + (χ_1/k)·ν_2·Σ_2f]·φ_2 = 0
    eq_a7_lhs = (
        (Sigma_1rem - chi_1 / k * nu_1 * Sigma_1f) * phi_1
        - (Sigma_12s + chi_1 / k * nu_2 * Sigma_2f) * phi_2
    )

    # Sum + apply χ_1 + χ_2 = 1.
    sum_eq = sp.expand(eq_a6_lhs + eq_a7_lhs)
    sum_eq_chi_eliminated = sp.simplify(sum_eq.subs(chi_2, 1 - chi_1))

    # The chi_1 dependence should drop out.
    chi_1_coeff = sp.simplify(sum_eq_chi_eliminated.coeff(chi_1))
    pass_chi_eliminates = (sp.simplify(chi_1_coeff) == 0)

    # The remaining equation is Sood Eq (A.13):
    eq_a13 = sp.simplify(sum_eq_chi_eliminated.subs(chi_1, 0))  # chi-free residue

    # Now restrict to no-upscatter (Σ_21s = 0).
    eq_a13_nou = eq_a13.subs(Sigma_21s, 0)

    # Solve for the ratio φ_2/φ_1.
    ratio = sp.symbols("ratio", positive=True)
    eq_a13_nou_in_ratio = sp.simplify(
        eq_a13_nou.subs(phi_2, ratio * phi_1) / phi_1
    )
    ratio_solutions = sp.solve(eq_a13_nou_in_ratio, ratio)
    assert len(ratio_solutions) == 1, (
        f"Expected single ratio root, got {ratio_solutions}"
    )
    ratio_derived = sp.simplify(ratio_solutions[0])

    # Sood Eq (A.15):
    ratio_eq_a15 = (
        (Sigma_1rem - nu_1 * Sigma_1f / k)
        / (nu_2 * Sigma_2f / k - Sigma_2rem + Sigma_12s)
    )
    ratio_eq_a15_simplified = sp.simplify(ratio_eq_a15)

    diff = sp.simplify(ratio_derived - ratio_eq_a15_simplified)
    pass_eq_a15 = (diff == 0)

    return {
        "name": "V_fn2.3: phi_2/phi_1 derivation matches Sood Eq (A.15)",
        "ratio_derived": ratio_derived,
        "ratio_eq_a15": ratio_eq_a15_simplified,
        "diff": diff,
        "pass_chi_eliminates": pass_chi_eliminates,
        "pass": pass_eq_a15 and pass_chi_eliminates,
    }


# ═══════════════════════════════════════════════════════════════════
# General multi-group infinite medium — Sood Eqs (A.55)-(A.59)
# ═══════════════════════════════════════════════════════════════════


def derive_kinf_mg_matrix_form() -> dict:
    r"""V_fnMG.1 — General MG balance reduces to single matrix inversion.

    Sood Eq (A.55) is the matrix balance for an infinite medium:

    .. math::

        \overline{\overline{\Sigma_t}}\,\bar\phi
        = \overline{\overline{\Sigma_s}}\,\bar\phi
        + \frac{1}{k}\, \bar\chi\, \overline{\nu\Sigma_f}\,\bar\phi

    Rearranging (Eq (A.56)):

    .. math::

        (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})\,\bar\phi
        = \frac{1}{k}\, \bar\chi\, (\overline{\nu\Sigma_f}\,\bar\phi)

    Inverting and projecting onto the production vector
    :math:`\overline{\nu\Sigma_f}` (Eqs (A.57)-(A.59)):

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
    Sood Eq (A.59) is the right form for the G=2 case where the formula
    *is* a clean closed form. We verify Eq (A.59) symbolically for G=2;
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

    # Sood Eq (A.55): Σ_t·φ = Σ_s·φ + (1/k)·χ·(νΣ_f·φ)
    A = Sigma_t - Sigma_s
    fission_source = chi * (nuSf * phi)  # (2x1) · (1x1) = 2x1
    eq_a55 = sp.Eq(A * phi, fission_source / k)

    # Sood Eq (A.59): k = νΣ_f · A^{-1} · χ
    A_inv = A.inv()
    k_eq_a59_matrix = nuSf * A_inv * chi  # (1x2)·(2x2)·(2x1) = (1x1) scalar
    k_eq_a59 = sp.simplify(k_eq_a59_matrix[0, 0])

    # Cross-check: dominant eigenvalue of M = A^{-1} · χ · νΣ_f
    M = A_inv * chi * nuSf  # (2x2) outer-product matrix, rank 1
    eigvals = list(M.eigenvals().keys())
    # Rank-1 matrix has one zero eigenvalue and one nonzero; trace = nonzero ev.
    M_trace = sp.simplify(sp.trace(M))
    # Trace of a rank-1 matrix outer(u, v) is dot(v, u). Same as Eq (A.59).
    diff = sp.simplify(M_trace - k_eq_a59)
    pass_trace = (diff == 0)

    # Also verify rank-1 structure (one zero eigenvalue).
    nonzero_evs = [e for e in eigvals if sp.simplify(e) != 0]
    pass_rank1 = (len(nonzero_evs) == 1)

    # Verify the surviving eigenvalue equals Eq (A.59).
    if pass_rank1:
        diff_ev = sp.simplify(nonzero_evs[0] - k_eq_a59)
        pass_ev = (diff_ev == 0)
    else:
        diff_ev = None
        pass_ev = False

    return {
        "name": "V_fnMG.1: Sood Eq (A.59) for G=2 — k = nuSf · A^{-1} · chi",
        "A": A,
        "k_eq_a59": k_eq_a59,
        "M_trace": M_trace,
        "diff_trace": diff,
        "pass_trace_equals_eq_a59": pass_trace,
        "pass_M_is_rank1": pass_rank1,
        "pass_eigenvalue_equals_eq_a59": pass_ev,
        "pass": pass_trace and pass_rank1 and pass_ev,
    }


def derive_kinf_mg_reduces_to_1g() -> dict:
    r"""V_fnMG.2 — MG formula with G=1 reduces to 1G Eq (A.2) bit-equal.

    The general MG formula Sood Eq (A.59):

    .. math::

        k = \overline{\nu\Sigma_f}\,
            (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})^{-1}\,
            \bar\chi

    at G=1 collapses to scalars: :math:`A` is the :math:`1\times 1`
    matrix :math:`(\Sigma_t - \Sigma_s)`, :math:`A^{-1} = 1/(\Sigma_t -
    \Sigma_s)`, :math:`\chi = (1)`, :math:`\nu\Sigma_f` is a scalar.
    Therefore :math:`k_{\rm Eq (A.59)} = \nu\Sigma_f / (\Sigma_t -
    \Sigma_s) = k_{\rm Eq (A.2)}`.

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
    k_eq_a59_g1 = sp.simplify((nuSf * A.inv() * chi)[0, 0])

    k_eq_a2 = nu_Sigma_f / (Sigma_t - Sigma_s)

    diff = sp.simplify(k_eq_a59_g1 - k_eq_a2)
    pass_id = (diff == 0)

    return {
        "name": "V_fnMG.2: Eq (A.59) with G=1 reduces to Eq (A.2)",
        "k_eq_a59_g1": k_eq_a59_g1,
        "k_eq_a2": k_eq_a2,
        "diff": diff,
        "pass": pass_id,
    }
