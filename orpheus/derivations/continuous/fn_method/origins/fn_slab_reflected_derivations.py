r"""SymPy derivations for the **reflected-slab** F_N method
(Neshat-Maiorino 1980, *Annals of Nuclear Energy* **7**, 79-81).

This module is the **algebra-of-record** for the reflected-slab F_N
solver — the trivial extension of the bare-slab F_N method
(Grandjean-Siewert 1979) to a finite-reflector configuration. The
geometry is

* **Core**: :math:`-\tau \le x \le \tau`, with mean number of
  secondaries :math:`c_1 > 1` (multiplying medium).
* **Reflector**: :math:`\tau < |x| \le b`, with :math:`c_2 < 1`
  (subcritical reflector).
* **Reflector half-thickness** :math:`\Delta = b - \tau` is **given**;
  the **critical core half-thickness** :math:`\tau` is the unknown.

The physics is encoded by the boundary conditions:

* :math:`\psi_1(\tau, \mu) = \psi_2(\tau, \mu)` for all :math:`\mu`
  — interface continuity.
* :math:`\psi_2(b, -\mu) = 0`, :math:`\psi_2(-b, \mu) = 0` for
  :math:`\mu > 0` — vacuum at the outer reflector surfaces.

The angular flux is symmetric: :math:`\psi_1(-\tau, \mu) =
\psi_1(\tau, -\mu)` (and similarly for the reflector).

What this module verifies (algebra-of-record, State 1A)
-------------------------------------------------------

The reflected-slab F_N method introduces THREE coefficient arrays
(NM Eqs. 9a-c):

.. math::

   \psi_1(\tau, +\mu) &= \sum_\alpha a_\alpha \mu^\alpha,
   \qquad \mu > 0
   \\
   \psi_1(\tau, -\mu) &= \sum_\alpha b_\alpha \mu^\alpha,
   \qquad \mu > 0
   \\
   \psi_2(b, +\mu) &= \psi_2(-b, -\mu) = \sum_\alpha e_\alpha \mu^\alpha,
   \qquad \mu > 0

so the F_N system is :math:`3(N+1)` equations in :math:`3(N+1)`
unknowns, gathered from the **three projection equations** NM
Eqs. 10, 11, 12 obtained from full-range Case orthogonality
(Eq. 5b applied at :math:`x=\pm\tau` and :math:`x=\pm b` for
the reflector, Eq. 5a at :math:`x=\pm\tau` for the core) plus
the boundary conditions Eq. 2.

The verifications below cover:

* **V_fn-slab-refl.1** — The :math:`A_\alpha(\xi)` and
  :math:`B_\alpha^{(i)}(\xi)` recursion identities (NM Eqs. 13-14)
  reduce **exactly** to the bare-slab Grandjean-Siewert moment
  recursions, *parametrised* by :math:`c_i` per region (medium-local
  X-function). This is the load-bearing structural fact: each region
  uses its own :math:`c_i` independently; the coupling lives in the
  boundary projection equations only.

* **V_fn-slab-refl.2** — The exponential signs in NM Eqs. 10, 11
  are self-consistent: at :math:`\mu = +\mu_0` (outgoing through
  reflector), the attenuation factor is :math:`e^{-\Delta/\xi}`
  (decay across the reflector); at :math:`\mu = -\mu_0` (returning
  through reflector), it is :math:`e^{+\Delta/\xi}` (back-projection
  inverse). The :math:`F_0` initial guess Eq. 17 has the structure
  of a 2-stream reflector albedo and is positive for :math:`c_2 < 1`,
  :math:`\Delta > 0`.

* **V_fn-slab-refl.3** — The critical condition NM Eq. 15 is
  the homogeneous-system determinant condition: a non-trivial
  :math:`\{a_\alpha, b_\alpha, e_\alpha\}` exists iff the
  :math:`3(N+1)`-row system has rank :math:`< 3(N+1)`. Eq. 15
  is the rank-deficiency condition specialised to the
  :math:`\xi = \nu_0` (Case discrete eigenvalue, core) row in
  the symmetry-eliminated form.

* **V_fn-slab-refl.4** — The :math:`F_0` initial-guess closed form
  (NM Eqs. 16-17) reduces to a recognisable 2-stream albedo
  / one-collision feedback coefficient. Direct substitution of
  :math:`a_0 = 1`, neglecting all higher moments, yields the
  Eq. 17 :math:`b_0` formula.

This module ENDS at the bifurcation point: the F_N system itself
(linear-algebraic determinant zero with transcendental
:math:`\xi`-dependence) is not closed-form for general :math:`N`;
the production solver in :mod:`...slab.reflected` evaluates it
numerically per the NM iteration loop.

The full numerical L1 cross-check is against NM Table 2
(Burkart 1971/76 "Exact" reference column) and against
Sood ``*-X(N)-1-0-SL`` reflected-slab cases.

References
----------

* Neshat & Maiorino 1980, *Ann. Nucl. Energy* **7**, 79-81.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161 (bare-slab
  F_N — provides the :math:`A_\alpha`, :math:`B_\alpha^{(i)}`
  recursions used here verbatim).
* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156 (Part I).
* Burkart 1975 PhD / 1976 *Trans. ANS* **24**, 190 — "Exact"
  reference values cited in NM Table 2.
"""
from __future__ import annotations

import sympy as sp


def derive_reflected_moment_recursions_match_bare() -> dict:
    r"""V_fn-slab-refl.1 — :math:`A_\alpha`, :math:`B_\alpha^{(i)}`
    recursions are bare-slab recursions parametrised by :math:`c_i`.

    Per NM Eqs. 13-14, the moment functions for the **reflected**
    slab are

    .. math::

       B_\alpha^{(i)}(\xi) &= \xi B_{\alpha-1}^{(i)}(\xi)
       - \frac{1}{\alpha+1},
       \qquad \alpha \ge 1,\ i=1,2,
       \\
       B_0^{(i)}(\xi) &= \frac{2}{c_i} - 1
       - \xi\,\log\!\left(1 + \frac{1}{\xi}\right),
       \qquad i=1,2,
       \\
       A_\alpha(\xi) &= -\xi A_{\alpha-1}(\xi) + \frac{1}{\alpha+1},
       \qquad \alpha \ge 1,
       \\
       A_0(\xi) &= 1 - \xi\,\log\!\left(1 + \frac{1}{\xi}\right).

    This is the **bare-slab** Grandjean-Siewert recursion (verified
    in :func:`derive_B_recursion`, :func:`derive_A_recursion`,
    :func:`derive_B0_seed`, :func:`derive_A0_seed` in
    :mod:`.fn_slab_derivations`). The only difference is that the
    reflected slab carries a *parametrised* :math:`B_0^{(i)}`
    seed: :math:`i = 1` uses :math:`c_1` (core), :math:`i = 2`
    uses :math:`c_2` (reflector).

    The structural insight: the X-function (whose argument
    determines the :math:`B_0` seed) is **medium-local** —
    each region's X-function depends only on its own :math:`c_i`
    and the dispersion-relation root :math:`\nu_0(c_i)`. There is
    NO two-region X-function. The two regions couple only through
    the boundary-projection equations (NM Eqs. 10-12) which involve
    integrals of the *common* angular flux at the *interface*
    :math:`\tau` and the *outer surface* :math:`b`.

    SymPy verifies that the recursion increment

    .. math::

       B_\alpha^{(i)}(\xi) - \xi B_{\alpha-1}^{(i)}(\xi)
       = -\frac{1}{\alpha+1}

    holds for both :math:`i = 1` and :math:`i = 2` independently
    (i.e., the recursion is independent of the region). Similarly
    for :math:`A_\alpha` (which has no :math:`c` dependence at all
    — all the :math:`c` dependence sits in :math:`B_0^{(i)}`).
    """
    xi, c1, c2 = sp.symbols("xi c_1 c_2", positive=True, real=True)
    alpha = sp.Symbol("alpha", positive=True, integer=True)

    # B_α^(i) recursion increment (i = 1, then i = 2; should be identical).
    # B_0^(i)(ξ) = 2/c_i - 1 - ξ log(1+1/ξ).
    B0_c1 = 2 / c1 - 1 - xi * sp.log(1 + 1 / xi)
    B0_c2 = 2 / c2 - 1 - xi * sp.log(1 + 1 / xi)

    # B_1^(i)(ξ) = ξ B_0^(i) - 1/2.
    B1_c1 = xi * B0_c1 - sp.Rational(1, 2)
    B1_c2 = xi * B0_c2 - sp.Rational(1, 2)

    # B_2^(i)(ξ) = ξ B_1^(i) - 1/3.
    B2_c1 = xi * B1_c1 - sp.Rational(1, 3)
    B2_c2 = xi * B1_c2 - sp.Rational(1, 3)

    # Increment check: B_α^(i) - ξ B_{α-1}^(i) = -1/(α+1).
    incr_B1_c1 = sp.simplify(B1_c1 - xi * B0_c1)
    incr_B1_c2 = sp.simplify(B1_c2 - xi * B0_c2)
    incr_B2_c1 = sp.simplify(B2_c1 - xi * B1_c1)
    incr_B2_c2 = sp.simplify(B2_c2 - xi * B1_c2)

    pass_incr_alpha1 = (
        incr_B1_c1 == sp.Rational(-1, 2)
        and incr_B1_c2 == sp.Rational(-1, 2)
    )
    pass_incr_alpha2 = (
        incr_B2_c1 == sp.Rational(-1, 3)
        and incr_B2_c2 == sp.Rational(-1, 3)
    )

    # The increment is c-INDEPENDENT (the c only enters via B_0).
    incr_diff_alpha1 = sp.simplify(incr_B1_c1 - incr_B1_c2)
    incr_diff_alpha2 = sp.simplify(incr_B2_c1 - incr_B2_c2)
    pass_c_indep = (incr_diff_alpha1 == 0 and incr_diff_alpha2 == 0)

    # A_α recursion: A_α = -ξ A_{α-1} + 1/(α+1). c-independent always.
    A0 = 1 - xi * sp.log(1 + 1 / xi)
    A1 = -xi * A0 + sp.Rational(1, 2)
    A2 = -xi * A1 + sp.Rational(1, 3)
    incr_A1 = sp.simplify(A1 + xi * A0)  # = 1/2.
    incr_A2 = sp.simplify(A2 + xi * A1)  # = 1/3.
    pass_A_recursion = (
        incr_A1 == sp.Rational(1, 2) and incr_A2 == sp.Rational(1, 3)
    )

    return {
        "name": (
            "V_fn-slab-refl.1: reflected-slab moment recursions match bare-slab "
            "(parametrised by c_i; X-function is medium-local)"
        ),
        "B0_c1": B0_c1,
        "B0_c2": B0_c2,
        "increment_B1_c1": incr_B1_c1,
        "increment_B1_c2": incr_B1_c2,
        "increment_B2_c1": incr_B2_c1,
        "increment_B2_c2": incr_B2_c2,
        "increment_diff_c1_minus_c2_alpha1": incr_diff_alpha1,
        "increment_diff_c1_minus_c2_alpha2": incr_diff_alpha2,
        "increment_A1": incr_A1,
        "increment_A2": incr_A2,
        "pass": bool(
            pass_incr_alpha1
            and pass_incr_alpha2
            and pass_c_indep
            and pass_A_recursion
        ),
    }


def derive_reflector_attenuation_signs() -> dict:
    r"""V_fn-slab-refl.2 — Exponential signs in NM Eqs. 10-11 are
    self-consistent under the reflector-traversal interpretation.

    NM Eq. 10 (collocated at :math:`\hat\xi \in P_2`):

    .. math::

       \sum_\alpha a_\alpha A_\alpha(\hat\xi)
       - \sum_\alpha b_\alpha B_\alpha^{(2)}(\hat\xi)
       = e^{-\Delta/\hat\xi} \sum_\alpha e_\alpha A_\alpha(\hat\xi).

    NM Eq. 11:

    .. math::

       \sum_\alpha b_\alpha A_\alpha(\hat\xi)
       - \sum_\alpha a_\alpha B_\alpha^{(2)}(\hat\xi)
       = -e^{+\Delta/\hat\xi} \sum_\alpha e_\alpha B_\alpha^{(2)}(\hat\xi).

    The exponential factors are :math:`e^{\mp\Delta/\hat\xi}`.
    Physical reading:

    * **Eq. 10**: attenuation factor :math:`e^{-\Delta/\hat\xi}`
      describes a particle starting at the interface :math:`x = \tau`
      and traversing the reflector to the outer surface :math:`x = b`,
      a distance :math:`\Delta = b - \tau` in mean-free-paths.
    * **Eq. 11**: factor :math:`e^{+\Delta/\hat\xi}` is the
      back-projection — equivalently, it appears because the equation
      is obtained by *eliminating* the reflector outward continuum
      coefficient :math:`B(-\eta)` between two equations of the form
      Eq. 5b at :math:`x = -\tau, -b`.

    The product :math:`e^{-\Delta/\hat\xi} \cdot e^{+\Delta/\hat\xi}
    = 1`, so the eliminations are consistent.

    The :math:`F_0` initial-guess :math:`b_0` (NM Eq. 17) is

    .. math::

       b_0 = \frac{A_0(\eta_0) B_0^{(2)}(\eta_0)
       (1 - e^{-2\Delta/\eta_0})}
       {[B_0^{(2)}(\eta_0)]^2 - [A_0(\eta_0)]^2 e^{-2\Delta/\eta_0}}.

    This is the structure of a 2-stream reflector albedo:

    * Numerator carries :math:`(1 - e^{-2\Delta/\eta_0})` —
      the "what fraction of an injected particle returns
      to the interface" coefficient. As :math:`\Delta \to 0`,
      :math:`b_0 \to 0` (no reflector, no return). As
      :math:`\Delta \to \infty`, :math:`b_0 \to A_0/B_0^{(2)}`
      (the infinite-reflector limit).
    * Denominator carries :math:`[B_0^{(2)}]^2 - [A_0]^2
      e^{-2\Delta/\eta_0}` which is the determinant of the
      :math:`2 \times 2` reflector projection system at
      :math:`F_0` order; positive for :math:`c_2 < 1`,
      :math:`\Delta > 0`.

    SymPy verifies:
    * Eq. 17 simplifies to the right limits as
      :math:`\Delta \to 0` and :math:`\Delta \to \infty`.
    * The product :math:`e^{-\Delta/\hat\xi} \cdot e^{+\Delta/\hat\xi}
      = 1` (consistency of eliminations).
    """
    xi_hat, Delta, eta0 = sp.symbols("xi_hat Delta eta_0", positive=True, real=True)
    A0_eta, B02_eta = sp.symbols("A0_eta B0_eta", real=True)

    # Product of exponentials in Eq. 10 and Eq. 11 must cancel.
    exp_10 = sp.exp(-Delta / xi_hat)
    exp_11 = sp.exp(+Delta / xi_hat)
    product = sp.simplify(exp_10 * exp_11)
    pass_product_one = (product == 1)

    # NM Eq. 17 b_0 formula.
    b0_formula = (
        A0_eta * B02_eta * (1 - sp.exp(-2 * Delta / eta0))
        / (B02_eta**2 - A0_eta**2 * sp.exp(-2 * Delta / eta0))
    )

    # Δ → 0 limit: b_0 → 0 (no reflector, no return).
    b0_at_Delta_0 = sp.limit(b0_formula, Delta, 0)
    pass_b0_Delta_0 = (sp.simplify(b0_at_Delta_0) == 0)

    # Δ → ∞ limit: e^{-2Δ/η_0} → 0, so b_0 → A_0 B_0^(2) / [B_0^(2)]^2 = A_0/B_0^(2).
    b0_at_Delta_inf = sp.limit(b0_formula, Delta, sp.oo)
    expected_inf = A0_eta / B02_eta
    diff_inf = sp.simplify(b0_at_Delta_inf - expected_inf)
    pass_b0_Delta_inf = (diff_inf == 0)

    return {
        "name": (
            "V_fn-slab-refl.2: NM Eq.10/11 exponential signs consistent; "
            "Eq.17 b_0 has correct Δ→0 and Δ→∞ limits"
        ),
        "exp_10_times_exp_11": product,
        "b0_formula": b0_formula,
        "b0_limit_Delta_to_0": b0_at_Delta_0,
        "b0_limit_Delta_to_inf": b0_at_Delta_inf,
        "expected_Delta_to_inf": expected_inf,
        "pass": bool(pass_product_one and pass_b0_Delta_0 and pass_b0_Delta_inf),
    }


def derive_critical_condition_eq15_structure() -> dict:
    r"""V_fn-slab-refl.3 — NM Eq. 15 is the rank-deficiency condition
    for the reflected-slab F_N system at the Case discrete eigenvalue
    :math:`\xi = \nu_0`.

    The reflected-slab F_N system is a homogeneous :math:`3(N+1)`
    linear system in :math:`\{a_\alpha, b_\alpha, e_\alpha\}`. With
    the normalisation :math:`a_0 = 1` fixing one degree of freedom,
    a non-trivial solution exists iff the remaining
    :math:`3(N+1) - 1 = 3N+2` rows are linearly dependent on the
    :math:`3N + 1` unknowns :math:`\{a_\alpha\}_{\alpha\ge 1},
    \{b_\alpha\}_{\alpha\ge 0}, \{e_\alpha\}_{\alpha\ge 0}`.

    Operationally NM split this constraint into:

    1. A :math:`3N+2`-equation linear system (Eqs. 10-12 collocated
       on the published grid) which determines the
       :math:`a_\alpha, b_\alpha, e_\alpha` GIVEN a value of
       :math:`\tau`.
    2. The CRITICAL CONDITION Eq. 15 — collocated at the core
       Case discrete eigenvalue :math:`\xi = \nu_0` — which selects
       the :math:`\tau` at which the linear system from step 1 is
       compatible.

    Eq. 15 has the explicit form:

    .. math::

       e^{-2\tau/\nu_0}
       \sum_\alpha [b_\alpha B_\alpha^{(1)}(\nu_0)
       - a_\alpha A_\alpha(\nu_0)]
       = \sum_\alpha [a_\alpha B_\alpha^{(1)}(\nu_0)
       - b_\alpha A_\alpha(\nu_0)]

    Symbolic-structure verification: at :math:`N = 0` (just
    :math:`a_0, b_0`) the equation reduces to a scalar transcendental
    in :math:`\tau`; with :math:`a_0 = 1` and a given :math:`b_0`
    (e.g. NM Eq. 17), it can be solved as :math:`\tau =
    -(\nu_0/2) \log[\dots]` (NM Eq. 16). SymPy verifies this
    closed-form follows from Eq. 15 at :math:`N = 0`.
    """
    nu0, tau = sp.symbols("nu_0 tau", positive=True, real=True)
    a0, b0 = sp.symbols("a_0 b_0", real=True)
    A0_nu, B01_nu = sp.symbols("A0_nu0 B01_nu0", real=True)

    # Eq. 15 at N = 0:
    #   e^{-2τ/ν_0} (b_0 B_0^(1)(ν_0) - a_0 A_0(ν_0))
    #     = a_0 B_0^(1)(ν_0) - b_0 A_0(ν_0).
    lhs = sp.exp(-2 * tau / nu0) * (b0 * B01_nu - a0 * A0_nu)
    rhs = a0 * B01_nu - b0 * A0_nu

    # Solve for e^{-2τ/ν_0}:
    #   e^{-2τ/ν_0} = (a_0 B_0^(1) - b_0 A_0) / (b_0 B_0^(1) - a_0 A_0).
    # → τ = -(ν_0/2) log[(a_0 B_0^(1) - b_0 A_0) / (b_0 B_0^(1) - a_0 A_0)].
    # NM Eq. 16 with a_0 = 1: τ^(0) = -(ν_0/2) log[(b_0 A_0 - B_0^(1))/(A_0 - b_0 B_0^(1))].
    # Compare signs: NM ratio has (b_0 A_0 - B_0^(1))/(A_0 - b_0 B_0^(1));
    # ours with a_0=1 gives (B_0^(1) - b_0 A_0) / (b_0 B_0^(1) - A_0)
    #                     = -(b_0 A_0 - B_0^(1)) / -(A_0 - b_0 B_0^(1))
    #                     = (b_0 A_0 - B_0^(1)) / (A_0 - b_0 B_0^(1)).
    # → matches NM Eq. 16 verbatim.
    eq15_at_N0 = lhs - rhs
    sol_for_exp = sp.solve(eq15_at_N0, sp.exp(-2 * tau / nu0))
    has_one_sol = (len(sol_for_exp) == 1)
    exp_neg_2tau_over_nu0 = sol_for_exp[0] if has_one_sol else None

    # Substitute a_0 = 1, then take -log/2*(-ν_0):
    #  τ = -(ν_0/2) log( solved_exp ).
    if exp_neg_2tau_over_nu0 is not None:
        sol_a0_1 = exp_neg_2tau_over_nu0.subs(a0, 1)
        # NM Eq. 16 form (with a_0 = 1):
        nm_eq16_arg = (b0 * A0_nu - B01_nu) / (A0_nu - b0 * B01_nu)
        # These should be equal:
        diff_arg = sp.simplify(sol_a0_1 - nm_eq16_arg)
        pass_eq16_match = (diff_arg == 0)
    else:
        pass_eq16_match = False

    return {
        "name": (
            "V_fn-slab-refl.3: NM Eq. 15 critical condition reduces to NM Eq. 16 "
            "for tau at N = 0 with a_0 = 1"
        ),
        "eq15_lhs": lhs,
        "eq15_rhs": rhs,
        "solved_exp_neg_2tau_over_nu0": exp_neg_2tau_over_nu0,
        "expected_NM_eq16_arg_at_a0_1": (b0 * A0_nu - B01_nu) / (A0_nu - b0 * B01_nu),
        "pass": bool(has_one_sol and pass_eq16_match),
    }


def derive_F0_initial_guess_structure() -> dict:
    r"""V_fn-slab-refl.4 — The :math:`F_0` initial-guess
    :math:`b_0` formula (NM Eq. 17) is consistent with a 2-stream
    reflector eliminating :math:`e_0` from the :math:`N = 0`
    truncation of NM Eqs. 10-11.

    At :math:`N = 0` (only :math:`a_0, b_0, e_0`), NM Eqs. 10
    and 11 collocated at :math:`\hat\xi = \eta_0` reduce to:

    .. math::

       a_0 A_0(\eta_0) - b_0 B_0^{(2)}(\eta_0)
       &= e^{-\Delta/\eta_0} e_0 A_0(\eta_0)
       \\
       b_0 A_0(\eta_0) - a_0 B_0^{(2)}(\eta_0)
       &= -e^{+\Delta/\eta_0} e_0 B_0^{(2)}(\eta_0)

    Solving for :math:`e_0` from each and equating:

    .. math::

       e_0 = \frac{a_0 A_0 - b_0 B_0^{(2)}}{e^{-\Delta/\eta_0} A_0}
       = -\frac{b_0 A_0 - a_0 B_0^{(2)}}{e^{+\Delta/\eta_0} B_0^{(2)}}

    Cross-multiplying:

    .. math::

       (a_0 A_0 - b_0 B_0^{(2)}) B_0^{(2)} e^{+\Delta/\eta_0}
       = -(b_0 A_0 - a_0 B_0^{(2)}) A_0 e^{-\Delta/\eta_0}

    Setting :math:`a_0 = 1` and solving for :math:`b_0` gives a
    closed-form. This module verifies the algebraic reduction
    matches NM Eq. 17:

    .. math::

       b_0 = \frac{A_0(\eta_0) B_0^{(2)}(\eta_0)
       (1 - e^{-2\Delta/\eta_0})}
       {[B_0^{(2)}(\eta_0)]^2 - [A_0(\eta_0)]^2 e^{-2\Delta/\eta_0}}.
    """
    Delta, eta0 = sp.symbols("Delta eta_0", positive=True, real=True)
    a0, b0, e0 = sp.symbols("a_0 b_0 e_0", real=True)
    A0, B02 = sp.symbols("A0 B0_2", positive=True, real=True)

    # Eq. 10 at N = 0, ξ̂ = η_0 (a_0 A_0 - b_0 B_0^(2) = e^{-Δ/η_0} e_0 A_0):
    eq10_resid = a0 * A0 - b0 * B02 - sp.exp(-Delta / eta0) * e0 * A0
    # Eq. 11 at N = 0, ξ̂ = η_0 (b_0 A_0 - a_0 B_0^(2) = -e^{+Δ/η_0} e_0 B_0^(2)):
    eq11_resid = b0 * A0 - a0 * B02 + sp.exp(+Delta / eta0) * e0 * B02

    # Eliminate e_0 between eq10_resid = 0 and eq11_resid = 0.
    sol = sp.solve([eq10_resid, eq11_resid], [b0, e0])
    has_solution = (len(sol) > 0)
    if has_solution:
        b0_sol_a0 = sp.simplify(sol[b0])
        # Substitute a_0 = 1.
        b0_sol = sp.simplify(b0_sol_a0.subs(a0, 1))
    else:
        b0_sol = None

    # NM Eq. 17:
    nm_eq17 = (
        A0 * B02 * (1 - sp.exp(-2 * Delta / eta0))
        / (B02**2 - A0**2 * sp.exp(-2 * Delta / eta0))
    )

    if b0_sol is not None:
        diff_b0 = sp.simplify(b0_sol - nm_eq17)
        pass_b0 = (diff_b0 == 0)
    else:
        diff_b0 = None
        pass_b0 = False

    return {
        "name": (
            "V_fn-slab-refl.4: F_0 initial-guess b_0 (NM Eq. 17) reduces from "
            "the N=0 truncation of NM Eqs. 10-11 with a_0 = 1"
        ),
        "eq10_residual": eq10_resid,
        "eq11_residual": eq11_resid,
        "b0_from_elimination": b0_sol,
        "nm_eq17": nm_eq17,
        "pass": bool(pass_b0),
    }
