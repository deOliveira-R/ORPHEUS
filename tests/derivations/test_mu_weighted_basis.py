r"""µ-weighted orthonormal polynomial basis on [0, 1] (Sanchez-McCormick
§III.F.1 Eq. 165).

Surface modes :math:`f^\rho_+(\mu)` form an orthonormal basis on
:math:`[0, 1]` under the µ-weighted inner product

.. math::

   \int_0^1 f^\rho(\mu)\,f^\nu(\mu)\,\mu\,\mathrm d\mu = \delta_{\rho\nu}

Equivalently — after scaling — these are Jacobi polynomials
:math:`P^{(0,1)}_n` shifted to :math:`[0, 1]`. The first mode is
:math:`f^0 = \sqrt{2}` (constant: :math:`\int_0^1 \mu\,\mathrm d\mu = 1/2`,
so :math:`C^2 \cdot 1/2 = 1 \Rightarrow C = \sqrt 2`).

This file pins three identities derived symbolically:

1. **Orthonormality** (symbolic SymPy: :math:`\langle f^i, f^j\rangle_\mu
   = \delta_{ij}` exactly).
2. **Numerical orthonormality** under high-order Gauss-Legendre quadrature
   to :math:`10^{-12}` (cross-checks the symbolic derivation against
   independent numerical integration).
3. **Cross-check vs shifted Jacobi** (:math:`f^n / P^{(0,1)}_n(2\mu - 1)`
   is constant for :math:`n = 0, 1, 2, 3`).

This is a math-origin promote — a general mathematical identity
verification, not tied to any particular Sphinx theory equation
``:label:`` (the basis is *used* in Sanchez-McCormick mode
construction but the identity itself stands on its own as a Gram-
Schmidt / Jacobi cross-check). Foundation marker, no ``verifies(...)``.

Promoted from ``derivations/diagnostics/derive_mu_weighted_basis.py``
(2026-04-30 triage).
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp


def _gram_schmidt_mu_weighted(max_n):
    """Gram-Schmidt of {1, µ, µ², ..., µ^max_n} with weight µ on [0, 1].

    Returns ``(mu_symbol, list_of_orthonormal_polys)`` where
    ``list_of_orthonormal_polys[i]`` is the i-th basis function,
    normalised so that :math:`\\int_0^1 \\mu\\,(f^i)^2\\,\\mathrm d\\mu = 1`.
    """
    mu = sp.Symbol("mu", positive=True)
    basis = [sp.Integer(1)]
    for k in range(1, max_n + 1):
        basis.append(mu**k)

    def inner(f, g):
        return sp.integrate(f * g * mu, (mu, 0, 1))

    ortho = []
    for v in basis:
        u = v
        for q in ortho:
            u = u - inner(v, q) / inner(q, q) * q
        ortho.append(sp.expand(u))

    f = []
    for u in ortho:
        nrm_sq = sp.simplify(inner(u, u))
        f.append(sp.simplify(u / sp.sqrt(nrm_sq)))
    return mu, f


@pytest.mark.foundation
def test_mu_weighted_basis():
    r"""Three-pronged check on the µ-weighted orthonormal basis.

    Bug it would catch: any arithmetic error in the Gram-Schmidt
    derivation, an incorrect normalisation constant (e.g. ``f^0 = 1``
    instead of :math:`\sqrt{2}`), or a wrong mapping to shifted Jacobi
    (which would surface in mode-recombination identities elsewhere
    in the rank-N Marshak / DP_N machinery).
    """
    mu, f = _gram_schmidt_mu_weighted(3)

    # 1. Symbolic orthonormality.
    for i in range(4):
        for j in range(4):
            val = sp.integrate(f[i] * f[j] * mu, (mu, 0, 1))
            expected = 1 if i == j else 0
            assert sp.simplify(val - expected) == 0, (
                f"f^{i}·f^{j} symbolic orthonormality failed: got {val}"
            )

    # 2. f^0 = √2.
    assert sp.simplify(f[0] - sp.sqrt(2)) == 0, (
        f"Expected f^0 = √2, got {f[0]}"
    )

    # 3. Numerical orthonormality at 1e-12 (independent of symbolic
    # integration — uses 64-pt Gauss-Legendre on [0, 1]).
    f_num = [sp.lambdify(mu, fi, "numpy") for fi in f]
    n_pts = 64
    x_leg, w_leg = np.polynomial.legendre.leggauss(n_pts)
    mu_pts = 0.5 * (x_leg + 1)
    w_pts = 0.5 * w_leg
    for i in range(4):
        for j in range(4):
            f_i_vals = f_num[i](mu_pts)
            f_j_vals = f_num[j](mu_pts)
            integ = float(np.sum(w_pts * mu_pts * f_i_vals * f_j_vals))
            expected = 1.0 if i == j else 0.0
            assert abs(integ - expected) < 1e-12, (
                f"⟨f^{i}, f^{j}⟩_µ = {integ:.3e} vs {expected} "
                f"(diff = {integ - expected:.3e})"
            )

    # 4. Shifted Jacobi cross-check: f^n ∝ P^{(0,1)}_n(2µ - 1) for
    # n = 0, 1, 2, 3 (proportionality constant absorbs the
    # normalisation difference between f^n and the standard Jacobi
    # convention).
    x = sp.Symbol("x")
    jac = [sp.jacobi(n, 0, 1, x) for n in range(4)]
    jac_shifted = [sp.simplify(p.subs(x, 2 * mu - 1)) for p in jac]
    for n in range(4):
        ratio = sp.simplify(f[n] / jac_shifted[n])
        assert sp.simplify(sp.diff(ratio, mu)) == 0, (
            f"f^{n} not proportional to shifted Jacobi P^{{(0,1)}}_{{{n}}}: "
            f"ratio = {ratio} is not constant"
        )
