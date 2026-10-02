r"""The continuous angular measure and its two lifts, in SymPy (Branch 1; #405 P1 step 6).

The direction sphere :math:`S^2` carries the measure
:math:`\mathrm d\Omega = \mathrm d\mu\,\mathrm d\varphi` over
:math:`\mu \in [-1, 1]`, :math:`\varphi \in [0, 2\pi)` (:data:`SPHERE`, the one
place the chart's ranges are written). On it:

* the **retraction** :math:`R q = \int_{S^2} q\,\mathrm d\Omega`, fibre
  integration over the direction (:func:`retraction`);
* its **section** :math:`E Q = Q / m`, with :math:`m = R 1` the mass of the
  measure, so that :math:`R \circ E = \mathrm{id}` (:func:`section`);
* the retraction's **adjoint** :math:`R^\dagger \Sigma = \Sigma`, the
  pullback, constant in direction (:func:`pullback`):
  :math:`\langle R^\dagger\Sigma, \psi\rangle_{S^2} = \langle \Sigma, R\psi\rangle`.

A source of rate :math:`Q` enters phase space as :math:`E Q`; a detector
:math:`\Sigma_d` as :math:`R^\dagger \Sigma_d`. The two differ by
:math:`R \circ R^\dagger = m`, and :math:`m = 4\pi` is DERIVED here by
integrating 1 over :data:`SPHERE`, never typed (the user's ruling of
2026-10-02: the :math:`4\pi` is the measure, not a convention).

This is the continuous counterpart of the discrete arrows
:class:`~orpheus.numerics.operator.AxisSectionOperator` and
:class:`~orpheus.numerics.operator.AxisPullbackOperator`, built independently
of them: it reads no quadrature and no production operator (X4). It is the
Branch-1 measure object #557 asks the manufactured sources and the
Green's-function sphere to lift through.

Verification functions, each pinned by a foundation test in
``tests/gates/derivations/test_angular_measure_symbolic.py``:

* :func:`derive_section_identity`: :math:`R(E Q) = Q`, exactly, per region and
  group of a ``Rational`` table;
* :func:`derive_adjoint_identity`: :math:`\langle R^\dagger\Sigma, \psi\rangle =
  \langle\Sigma, R\psi\rangle`, exactly, for a ψ anisotropic in both angles.
"""

from __future__ import annotations

import sympy as sp

R_POS = sp.Symbol("r", real=True)
MU = sp.Symbol("mu", real=True)
PHI = sp.Symbol("phi", real=True)

SPHERE = ((MU, -1, 1), (PHI, 0, 2 * sp.pi))
"""The integration domain of :math:`\\mathrm d\\Omega = \\mathrm d\\mu\\,\\mathrm d\\varphi`."""


def retraction(q: sp.Expr) -> sp.Expr:
    r""":math:`R q = \int_{S^2} q\,\mathrm d\Omega`."""
    return sp.integrate(q, *SPHERE)


def mass() -> sp.Expr:
    r"""The mass of the measure, :math:`m = R 1`."""
    return retraction(sp.Integer(1))


def section(rate: sp.Expr) -> sp.Expr:
    r""":math:`E Q = Q / m`: the isotropic function whose integral over the sphere is ``rate``."""
    return rate / mass()


def pullback(response: sp.Expr) -> sp.Expr:
    r""":math:`R^\dagger \Sigma = \Sigma`, constant in direction (the retraction's adjoint)."""
    return response


def derive_section_identity() -> dict:
    r"""V_E: :math:`R(E Q) = Q` exactly, for every entry of a 3-region, 2-group ``Rational`` table."""
    table = [[sp.Rational(7, 3), sp.Rational(1, 5)], [sp.Integer(0), sp.Rational(-2, 7)], [sp.Rational(11, 13), 1]]
    residuals = [sp.simplify(retraction(section(Q)) - Q) for row in table for Q in row]
    return {
        "name": "V_E: the section is a right inverse of the retraction",
        "mass": mass(),
        "residuals": residuals,
        "pass": all(res == 0 for res in residuals),
    }


def derive_adjoint_identity() -> dict:
    r"""V_R†: :math:`\int_0^1\!\!\int_{S^2} (R^\dagger\Sigma)\,\psi\,\mathrm d\Omega\,\mathrm dr
    = \int_0^1 \Sigma\,(R\psi)\,\mathrm dr`, exactly, with
    :math:`\Sigma(r) = 1 + r` and :math:`\psi = (1 + \mu + \mu^2\cos\varphi)\,e^{-r}`."""
    sigma = 1 + R_POS
    psi = (1 + MU + MU**2 * sp.cos(PHI)) * sp.exp(-R_POS)
    lhs = sp.integrate(sp.integrate(pullback(sigma) * psi, *SPHERE), (R_POS, 0, 1))
    rhs = sp.integrate(sigma * retraction(psi), (R_POS, 0, 1))
    wrong_arrow = sp.integrate(sp.integrate(section(sigma) * psi, *SPHERE), (R_POS, 0, 1))
    return {
        "name": "V_R†: the pullback is the adjoint of the retraction",
        "lhs": lhs,
        "rhs": rhs,
        "ratio_if_lifted_by_the_section": sp.simplify(wrong_arrow / rhs),
        "pass": sp.simplify(lhs - rhs) == 0,
    }


__all__ = [
    "MU",
    "PHI",
    "R_POS",
    "SPHERE",
    "derive_adjoint_identity",
    "derive_section_identity",
    "mass",
    "pullback",
    "retraction",
    "section",
]
