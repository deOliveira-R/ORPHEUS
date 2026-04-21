"""Derive the µ-weighted orthonormal polynomial basis on [0, 1].

Created by numerics-investigator on 2026-04-21 for Issue #119.

Sanchez-McCormick §III.F.1 Eq. 165 defines surface modes
f^ρ_+(µ) as an orthonormal basis on [0, 1] under the µ-weighted inner
product:

    ∫_0^1 f^ρ(µ) · f^ν(µ) · µ dµ = δ_ρν

Equivalently (after scaling), these are Jacobi polynomials P^{(0,1)}_n
shifted to [0, 1]. The first mode f^0 is constant:
f^0 = √2 (since ∫_0^1 µ dµ = 1/2, so C² · (1/2) = 1 ⇒ C = √2).

We derive f^0, f^1, f^2, f^3 symbolically via Gram-Schmidt, then
cross-check against shifted Jacobi P^{(0,1)}_n (which is the
standard-table identification).

Promotion: this is a general mathematical verification and should go
into tests/derivations/ if the recipe closes.
"""
import numpy as np
import pytest
import sympy as sp


def _gram_schmidt_mu_weighted(max_n):
    """Gram-Schmidt {1, µ, µ², ..., µ^max_n} with weight µ on [0,1]."""
    mu = sp.Symbol("mu", positive=True)
    basis = [sp.Integer(1)]
    for k in range(1, max_n + 1):
        basis.append(mu**k)

    # Orthogonalize under ⟨f, g⟩ = ∫_0^1 f·g·µ dµ.
    def inner(f, g):
        return sp.integrate(f * g * mu, (mu, 0, 1))

    ortho = []
    for v in basis:
        u = v
        for q in ortho:
            u = u - inner(v, q) / inner(q, q) * q
        ortho.append(sp.expand(u))

    # Normalize: f^n = u_n / √⟨u_n, u_n⟩.
    f = []
    for u in ortho:
        nrm_sq = sp.simplify(inner(u, u))
        f.append(sp.simplify(u / sp.sqrt(nrm_sq)))
    return mu, f


def test_mu_weighted_basis_orthonormality():
    """f^ρ must satisfy ∫ µ f^ρ f^ν dµ = δ_ρν to machine precision."""
    mu, f = _gram_schmidt_mu_weighted(3)
    for i in range(4):
        for j in range(4):
            val = sp.integrate(f[i] * f[j] * mu, (mu, 0, 1))
            expected = 1 if i == j else 0
            assert sp.simplify(val - expected) == 0, (
                f"f^{i}·f^{j} orthonormality failed: got {val}"
            )


def test_mu_weighted_basis_f0_is_sqrt2():
    """f^0 must be the constant √2 under µ-weighted normalisation."""
    mu, f = _gram_schmidt_mu_weighted(0)
    assert sp.simplify(f[0] - sp.sqrt(2)) == 0, (
        f"Expected f^0 = √2, got {f[0]}"
    )


def test_mu_weighted_basis_matches_shifted_jacobi_0_1():
    """Verify f^n matches the Jacobi P^{(0,1)}_n shifted to [0, 1]
    up to a normalization constant.

    Standard-table Jacobi P^{(alpha,beta)}_n(x) on [-1, 1] has weight
    (1-x)^alpha (1+x)^beta. Transform x → 2µ - 1 maps to [0, 1] with
    weight (2(1-µ))^0 · (2µ)^1 = 2µ → proportional to weight µ.
    So P^{(0,1)}_n(2µ - 1) is the shifted Jacobi; both sides
    orthogonal under ∫ µ · (·) dµ, up to a scalar constant.
    """
    mu, f_mu = _gram_schmidt_mu_weighted(3)

    x = sp.Symbol("x")
    # Jacobi polynomials P^{(0,1)}_n(x) on [-1, 1]
    jac = [sp.jacobi(n, 0, 1, x) for n in range(4)]
    # Shift to [0, 1]: x = 2µ - 1
    jac_shifted = [sp.simplify(p.subs(x, 2 * mu - 1)) for p in jac]

    # Check proportionality: ratio f^n / J^n must be constant.
    for n in range(4):
        ratio = sp.simplify(f_mu[n] / jac_shifted[n])
        # Constant check: derivative w.r.t. mu is zero.
        d = sp.diff(ratio, mu)
        assert sp.simplify(d) == 0, (
            f"f^{n} not proportional to shifted Jacobi P^{{(0,1)}}_{{{n}}}: "
            f"ratio = {ratio} is not constant"
        )


def test_mu_weighted_basis_numerical_orthonormality_1e12():
    """Numerical orthonormality check to 1e-12."""
    mu, f = _gram_schmidt_mu_weighted(3)

    # Convert to numerical callables.
    f_num = [sp.lambdify(mu, fi, "numpy") for fi in f]

    # High-order Gauss-Legendre on [0, 1].
    n_pts = 64
    x_leg, w_leg = np.polynomial.legendre.leggauss(n_pts)
    mu_pts = 0.5 * (x_leg + 1)  # map [-1,1] → [0,1]
    w_pts = 0.5 * w_leg  # map weight

    for i in range(4):
        for j in range(4):
            f_i_vals = f_num[i](mu_pts)
            f_j_vals = f_num[j](mu_pts)
            integ = np.sum(w_pts * mu_pts * f_i_vals * f_j_vals)
            expected = 1.0 if i == j else 0.0
            assert abs(integ - expected) < 1e-12, (
                f"⟨f^{i}, f^{j}⟩_µ = {integ:.3e} vs {expected} "
                f"(diff = {integ-expected:.3e})"
            )


def print_symbolic_basis():
    """Print the first 4 µ-weighted orthonormal polynomials."""
    mu, f = _gram_schmidt_mu_weighted(3)
    print("µ-weighted orthonormal basis on [0, 1]:")
    print("  ⟨f^ρ, f^ν⟩_µ = ∫_0^1 f^ρ(µ) · f^ν(µ) · µ dµ = δ_ρν")
    for n, fn in enumerate(f):
        print(f"  f^{n}(µ) = {sp.nsimplify(fn, rational=False)}")
        print(f"         = {sp.simplify(fn)}")
    print()
    # Check values at µ=0 and µ=1
    for n, fn in enumerate(f):
        v0 = float(fn.subs(mu, 0))
        v1 = float(fn.subs(mu, 1))
        print(f"  f^{n}(0) = {v0:.6f},  f^{n}(1) = {v1:.6f}")


if __name__ == "__main__":
    print_symbolic_basis()
