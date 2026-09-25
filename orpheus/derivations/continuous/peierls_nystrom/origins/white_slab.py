r"""SymPy origin of the slab white-boundary closed form
:eq:`peierls-white-bc-slab` and of the two exponential-integral
antiderivatives it rests on.

A homogeneous pure-absorber slab :math:`[0, L]` with a uniform unit
source and a Mark (isotropic re-entry, albedo 1) white boundary on both
faces has the flat flux :math:`\varphi \equiv 1/\Sigma_t`. The closed
form is implemented, in arbitrary precision, by
:func:`~orpheus.derivations.continuous.peierls_nystrom.reference.slab_uniform_source_white_bc_analytical`;
its derivation is on the theory page, section
``peierls-slab-white-bc-analytical``. This module carries that
derivation symbolically, one step per ``derive_*`` function.

The load-bearing step is the antiderivative

.. math::

   \int_0^{\tau} E_2(u)\,\mathrm du \;=\; \tfrac12 - E_3(\tau).

ERR-032 (``docs/theory/verification/error_catalog.rst``) shipped the
closed form derived with :math:`1 - E_3(\tau)` instead. The two
candidates have the SAME derivative, :math:`E_2(\tau)`, because
:math:`E_3'(\tau) = -E_2(\tau)`; they differ only in the integration
constant, which is fixed by the value at :math:`\tau = 0`
(:math:`E_3(0) = \tfrac12`). The identity is therefore proven here by
the fundamental theorem of calculus: the derivative residual is zero
and the limit at zero is zero. ``sympy.integrate(expint(2, u))`` is
deliberately NOT used: it returns an ``Ei(tau*exp_polar(I*pi))`` form
that ``simplify`` does not reduce to the right-hand side (measured
2026-09-25, ``scratch/reference_architecture/p0probe/sympy_e2.py``), so
an integrate-based identity would be a proof that cannot close.

Scope. This module proves the ALGEBRA of the closed form; it evaluates
no production code, so no defect in a code path can make it read
differently (the ``instrument-doctrine`` X3 point: it carries no
``catches`` marker). The code-side catchers of ERR-032 are the
Wigner–Seitz identity row and the quadrature partial-current balance
row in ``tests/gates/derivations/test_peierls_reference.py``.
"""

from __future__ import annotations

import sympy as sp

__all__ = [
    "derive_en_antiderivatives",
    "derive_slab_white_bc_flux",
    "slab_white_bc_flux",
]


def _ftc_residuals(antiderivative: sp.Expr, integrand: sp.Expr,
                   tau: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    r"""The two fundamental-theorem-of-calculus conditions for
    :math:`F(\tau) = \int_0^\tau f(u)\,\mathrm du`: the derivative
    residual :math:`F'(\tau) - f(\tau)` (simplified) and the value
    :math:`\lim_{\tau \to 0^+} F(\tau)`. Both are zero iff ``F`` is the
    antiderivative anchored at 0."""
    derivative_residual = sp.simplify(sp.diff(antiderivative, tau) - integrand)
    value_at_zero = sp.limit(antiderivative, tau, 0, "+")
    return derivative_residual, value_at_zero


def derive_en_antiderivatives() -> dict:
    r"""V_ws1 — the antiderivatives of :math:`E_1` and :math:`E_2`
    anchored at zero, and the ERR-032 candidate refuted.

    Proves:

    - :math:`\int_0^\tau E_2(u)\,\mathrm du = \tfrac12 - E_3(\tau)`
      (derivative residual 0, value at 0 equal to 0);
    - :math:`\int_0^\tau E_1(u)\,\mathrm du = 1 - E_2(\tau)` (the same
      two conditions; this is the volume term of the slab Peierls
      equation);
    - the ERR-032 candidate :math:`1 - E_3(\tau)` has derivative
      residual 0 as well, and value :math:`\tfrac12 \neq 0` at zero: it
      is an antiderivative of :math:`E_2`, anchored at the wrong
      constant. This negative leg is what makes the positive leg
      discriminating: a derivative check alone passes both candidates.

    Returns a dict with the residuals, the values at zero and ``pass``
    (True iff every positive condition holds AND the negative leg is
    refuted by its value at zero while agreeing in derivative).
    """
    tau = sp.Symbol("tau", positive=True)
    e1, e2, e3 = (sp.expint(n, tau) for n in (1, 2, 3))

    int_e2 = sp.Rational(1, 2) - e3
    int_e1 = 1 - e2
    int_e2_err032 = 1 - e3

    e2_residual, e2_at_zero = _ftc_residuals(int_e2, e2, tau)
    e1_residual, e1_at_zero = _ftc_residuals(int_e1, e1, tau)
    err032_residual, err032_at_zero = _ftc_residuals(int_e2_err032, e2, tau)

    positive = (e2_residual == 0 and e2_at_zero == 0
                and e1_residual == 0 and e1_at_zero == 0)
    negative_refuted = (err032_residual == 0
                        and err032_at_zero == sp.Rational(1, 2))
    return {
        "name": "V_ws1: int_0^tau E_2 = 1/2 - E_3 and int_0^tau E_1 = 1 - E_2; "
                "the ERR-032 candidate 1 - E_3 is off by the constant 1/2",
        "tau": tau,
        "int_E2": int_e2,
        "int_E1": int_e1,
        "int_E2_err032": int_e2_err032,
        "E2_derivative_residual": e2_residual,
        "E2_value_at_zero": e2_at_zero,
        "E1_derivative_residual": e1_residual,
        "E1_value_at_zero": e1_at_zero,
        "err032_derivative_residual": err032_residual,
        "err032_value_at_zero": err032_at_zero,
        "pass": bool(positive and negative_refuted),
    }


def slab_white_bc_flux(int_e2_of: sp.Lambda, sig_t: sp.Symbol,
                       L: sp.Symbol, x: sp.Symbol) -> sp.Expr:
    r"""The slab white-boundary scalar flux assembled from the Peierls
    equation and the Mark partial-current balance, for a given
    antiderivative :math:`\tau \mapsto \int_0^\tau E_2`.

    .. math::

       \varphi(x) = \frac{1}{2\Sigma_t}\bigl[(1 - E_2(\Sigma_t x))
                    + (1 - E_2(\Sigma_t(L-x)))\bigr]
                    + 2 J^-\bigl[E_2(\Sigma_t x) + E_2(\Sigma_t(L-x))\bigr],
       \qquad
       J^- = \frac{J^+_{\rm vol}}{1 - 2E_3(\Sigma_t L)},
       \quad
       J^+_{\rm vol} = \frac{1}{2\Sigma_t}\int_0^{\Sigma_t L} E_2.

    The volume term is :math:`\tfrac12\int_0^L E_1(\Sigma_t|x-x'|)\,\mathrm dx'`
    split at :math:`x` and changed to :math:`u = \Sigma_t|x - x'|` (the
    one step done by hand), then closed with
    :math:`\int_0^\tau E_1 = 1 - E_2(\tau)` from
    :func:`derive_en_antiderivatives`. The antiderivative of
    :math:`E_2` is a parameter so that the ERR-032 candidate can be
    assembled through the same algebra.
    """
    tau_left = sig_t * x
    tau_right = sig_t * (L - x)
    tau_L = sig_t * L
    volume = (2 - sp.expint(2, tau_left) - sp.expint(2, tau_right)) / (2 * sig_t)
    j_plus_volume = int_e2_of(tau_L) / (2 * sig_t)
    j_minus = j_plus_volume / (1 - 2 * sp.expint(3, tau_L))
    return volume + 2 * j_minus * (sp.expint(2, tau_left) + sp.expint(2, tau_right))


def derive_slab_white_bc_flux() -> dict:
    r"""V_ws2 — :eq:`peierls-white-bc-slab`: the flux assembled by
    :func:`slab_white_bc_flux` with :math:`\int_0^\tau E_2 = \tfrac12 - E_3`
    is identically :math:`1/\Sigma_t`; with the ERR-032 antiderivative
    :math:`1 - E_3` it is not.

    The collapse is the cancellation
    :math:`(\tfrac12 - E_3)/(1 - 2E_3) = \tfrac12`, so
    :math:`2J^- = 1/(2\Sigma_t)` and the :math:`E_2` terms of the volume
    and the boundary contributions cancel. The negative leg's residual
    is evaluated at the ERR-032 arm's probe point
    :math:`(x, L, \Sigma_t) = (0.3, 1, 1)`, where the catalogue's wrong
    closed form reads 1.45097 against the true 1, so its non-vanishing
    is a number, not a failure of ``simplify``.

    Returns a dict with both residuals and ``pass``.
    """
    sig_t, L = sp.symbols("Sigma_t L", positive=True)
    x = sp.Symbol("x", positive=True)
    u = sp.Symbol("u", positive=True)
    derived = derive_en_antiderivatives()
    tau = derived["tau"]
    int_e2 = sp.Lambda(u, derived["int_E2"].subs(tau, u))
    int_e2_err032 = sp.Lambda(u, derived["int_E2_err032"].subs(tau, u))

    phi = slab_white_bc_flux(int_e2, sig_t, L, x)
    residual = sp.simplify(phi - 1 / sig_t)

    phi_err032 = slab_white_bc_flux(int_e2_err032, sig_t, L, x)
    probe = {x: sp.Rational(3, 10), L: 1, sig_t: 1}
    err032_residual_at_probe = sp.N(
        (phi_err032 - 1 / sig_t).subs(probe), 30)

    return {
        "name": "V_ws2: slab white-BC flux is identically 1/Sigma_t "
                "(eq. peierls-white-bc-slab); the ERR-032 form is not",
        "phi": phi,
        "residual": residual,
        "err032_residual_at_probe": err032_residual_at_probe,
        "pass": bool(derived["pass"] and residual == 0
                     and abs(err032_residual_at_probe) > sp.Rational(1, 10)),
    }
