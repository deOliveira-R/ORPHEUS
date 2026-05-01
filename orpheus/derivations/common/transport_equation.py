"""Symbolic transport equation — the math root for both verification paths.

This module is the **single source of truth** for the symbolic form of
the linear Boltzmann transport equation that ORPHEUS verifies. Both
the discrete (Path 1) and continuous (Path 2) verification paths
branch from the equations defined here:

- Path 1 — *discrete production solvers* (``orpheus.derivations.discrete.*``):
  start from the canonical first-order form below and apply the
  specific spatial / angular discretisation each method commits to
  (S\\ :sub:`N` discrete-ordinate sweep, MOC characteristic integration,
  CP collision-probability projection, finite-volume diffusion).

- Path 2 — *continuous reference derivations*
  (``orpheus.derivations.continuous.*``):
  start from the **same** canonical form and apply
  *analytical* / *semi-analytical* techniques that bypass spatial
  discretisation entirely (homogeneous infinite-medium algebra,
  flat-source CP integral identities, Peierls integral form via
  Nyström / collocation, Method of Manufactured Solutions).

Equality of the two paths at the eigenvalue level is what verifies
the production solvers — the references in Path 2 are the truth
that Path 1 commits to converge to.

Canonical form
--------------

The mono-energetic, time-independent, isotropic-scattering linear
Boltzmann equation in 1-D Cartesian geometry, written in *first-order*
operator form:

.. math::
   :label: math-root-transport

   \\mu \\frac{\\partial \\psi}{\\partial r}(r, \\mu)
   + \\Sigma_t(r) \\, \\psi(r, \\mu)
   = \\frac{1}{2} \\Sigma_s(r) \\, \\phi(r)
   + \\frac{\\chi}{2 k} \\nu \\Sigma_f(r) \\, \\phi(r)
   + \\frac{Q(r)}{2}

with the scalar flux

.. math::

   \\phi(r) = \\int_{-1}^{1} \\psi(r, \\mu)\\, d\\mu .

In multi-group form, group-to-group scattering and fission spectrum
:math:`\\chi_g` enter as sums over the source group :math:`g'`; the
scalar-flux integral runs over the same :math:`(\\mu, E)` domain.
Curvilinear geometries (cylinder, sphere) replace the streaming
operator's :math:`\\mu \\partial_r` with the appropriate angular
redistribution + radial advection terms; this module's
:func:`streaming_operator` returns the Cartesian form by default and
documents the curvilinear extensions in its docstring.

Boundary conditions
-------------------

Three BC families are exposed symbolically:

- vacuum (no incoming current): :math:`\\psi(r_b, \\mu) = 0` for
  :math:`\\mu` pointing *into* the domain at boundary :math:`r_b`.
- reflective (specular mirror): :math:`\\psi(r_b, \\mu) =
  \\psi(r_b, -\\mu)` at the symmetry plane.
- white (isotropic re-emission): the incoming current at :math:`r_b`
  is uniformly redistributed in the inward hemisphere.

Self-test
---------

Run ``python -m orpheus.derivations.common.transport_equation``
to print every symbolic form this module defines.

See also
--------

- ``orpheus.derivations.discrete`` — production solver discretisations.
- ``orpheus.derivations.continuous`` — continuous reference derivations.
- ``docs/theory/verification.rst`` for the architectural treatment
  of the two paths and why this module is the contract between them.
"""

from __future__ import annotations

import sympy as sp


# ─────────────────────────────────────────────────────────────────────────
# Symbols
# ─────────────────────────────────────────────────────────────────────────
# Independent variables
r = sp.Symbol("r", real=True, nonnegative=True)
mu = sp.Symbol("mu", real=True)
E = sp.Symbol("E", real=True, positive=True)

# Scalar functions of (r) — the unknowns of the *integrated* form
phi = sp.Function("phi")(r)         # scalar flux φ(r)
Q = sp.Function("Q")(r)             # external (fixed) source Q(r)

# Angular flux ψ(r, µ) — the unknown of the *first-order* form
psi = sp.Function("psi")(r, mu)

# Macroscopic cross-sections — functions of r so the symbolic form
# admits heterogeneous (region-piecewise) materials
Sigma_t = sp.Function("Sigma_t")(r)   # total
Sigma_s = sp.Function("Sigma_s")(r)   # scattering (single-group form)
Sigma_f = sp.Function("Sigma_f")(r)   # fission
nu = sp.Symbol("nu", positive=True)   # neutrons per fission
chi = sp.Symbol("chi", real=True)     # fission spectrum (single-group: 1)

# Eigenvalue
k = sp.Symbol("k", positive=True)


# ─────────────────────────────────────────────────────────────────────────
# Operators
# ─────────────────────────────────────────────────────────────────────────
def streaming_operator(geometry: str = "cartesian"):
    """Return the streaming + collision operator :math:`T \\psi`.

    Parameters
    ----------
    geometry : {"cartesian", "cylindrical", "spherical"}
        Cartesian: :math:`\\mu \\partial_r \\psi + \\Sigma_t \\psi`.
        Curvilinear forms add angular-redistribution terms; only the
        Cartesian form is returned in closed symbolic form here. The
        curvilinear extensions are derived inside the per-geometry
        modules (``continuous/peierls/{cylinder,sphere}.py``,
        ``discrete/sn/balance.py``).
    """
    if geometry == "cartesian":
        return mu * sp.diff(psi, r) + Sigma_t * psi
    raise NotImplementedError(
        f"Curvilinear streaming operator for {geometry!r} is "
        "documented in the per-geometry modules; this root form is "
        "Cartesian-only."
    )


def scattering_source_isotropic():
    """Isotropic scattering source: :math:`\\frac{1}{2} \\Sigma_s \\phi`."""
    return sp.Rational(1, 2) * Sigma_s * phi


def fission_source_isotropic():
    """Fission source: :math:`\\frac{\\chi}{2 k} \\nu \\Sigma_f \\phi`."""
    return chi * nu * Sigma_f * phi / (2 * k)


def external_source_isotropic():
    """Isotropic external source: :math:`Q(r) / 2`."""
    return Q / 2


def transport_equation(geometry: str = "cartesian", with_fission: bool = True,
                       with_external: bool = True):
    """Assemble the canonical first-order transport equation as ``LHS - RHS``.

    Returns a SymPy ``Eq(LHS, RHS)``. Toggling ``with_fission`` /
    ``with_external`` produces the eigenvalue / fixed-source / pure
    streaming sub-cases that the verification cases in this package
    target individually.
    """
    lhs = streaming_operator(geometry)
    rhs = scattering_source_isotropic()
    if with_fission:
        rhs = rhs + fission_source_isotropic()
    if with_external:
        rhs = rhs + external_source_isotropic()
    return sp.Eq(lhs, rhs)


def scalar_flux_definition():
    """Return :math:`\\phi(r) = \\int_{-1}^{1} \\psi(r, \\mu) d\\mu` symbolically."""
    return sp.Eq(phi, sp.Integral(psi, (mu, -1, 1)))


# ─────────────────────────────────────────────────────────────────────────
# Multi-group operators
# ─────────────────────────────────────────────────────────────────────────
def multigroup_scattering_source(num_groups: int):
    """Return the symbolic group-to-group scattering source for group :math:`g`.

    Returns a tuple ``(rhs_per_group, phi_g_symbols, Sigma_s_gp_to_g_symbols)``
    where ``rhs_per_group[g] = (1/2) Σ_{g'} Σ_s^{g'→g} φ_{g'}``.

    The returned objects are pure ``sp.IndexedBase`` symbols suitable
    for substitution by either path's discretisation routine.
    """
    g = sp.Symbol("g", integer=True, nonnegative=True)
    gp = sp.Symbol("gprime", integer=True, nonnegative=True)
    phi_g = sp.IndexedBase("phi")
    Sigma_s_gp_to_g = sp.IndexedBase("Sigma_s")
    rhs = [
        sp.Rational(1, 2) * sp.Sum(
            Sigma_s_gp_to_g[gp, g_idx] * phi_g[gp], (gp, 0, num_groups - 1)
        )
        for g_idx in range(num_groups)
    ]
    return rhs, phi_g, Sigma_s_gp_to_g


def multigroup_fission_source(num_groups: int):
    """Return the symbolic multi-group fission source for group :math:`g`.

    Returns ``rhs_per_group[g] = (χ_g / 2k) Σ_{g'} νΣ_f^{g'} φ_{g'}``.
    """
    gp = sp.Symbol("gprime", integer=True, nonnegative=True)
    phi_g = sp.IndexedBase("phi")
    chi_g = sp.IndexedBase("chi")
    nuSigma_f_gp = sp.IndexedBase("nuSigma_f")
    rhs = [
        chi_g[g_idx] / (2 * k) * sp.Sum(
            nuSigma_f_gp[gp] * phi_g[gp], (gp, 0, num_groups - 1)
        )
        for g_idx in range(num_groups)
    ]
    return rhs, phi_g, chi_g, nuSigma_f_gp


# ─────────────────────────────────────────────────────────────────────────
# Boundary conditions (symbolic forms)
# ─────────────────────────────────────────────────────────────────────────
def vacuum_bc(r_boundary: sp.Expr, mu_in_sign: int):
    """Vacuum BC: :math:`\\psi(r_b, \\mu) = 0` for :math:`\\mu` in the inward
    hemisphere defined by ``mu_in_sign in {+1, -1}``."""
    if mu_in_sign not in (-1, +1):
        raise ValueError("mu_in_sign must be +1 or -1")
    return sp.Eq(psi.subs(r, r_boundary), 0)


def reflective_bc(r_boundary: sp.Expr):
    """Reflective (specular) BC: :math:`\\psi(r_b, \\mu) = \\psi(r_b, -\\mu)`."""
    return sp.Eq(psi.subs(r, r_boundary),
                 psi.subs([(r, r_boundary), (mu, -mu)]))


def white_bc(r_boundary: sp.Expr):
    """White (isotropic re-emission) BC.

    The incoming partial current at :math:`r_b` is uniformly
    redistributed across the inward hemisphere. Symbolically:

    .. math::
       \\psi^{\\mathrm{in}}(r_b, \\mu) = \\frac{2}{\\pi}
       \\int_0^1 |\\mu'| \\psi^{\\mathrm{out}}(r_b, \\mu')\\, d\\mu' .
    """
    mu_p = sp.Symbol("muprime", real=True)
    psi_out = sp.Function("psi_out")(r_boundary, mu_p)
    psi_in = sp.Function("psi_in")(r_boundary, mu)
    rhs = sp.Rational(2, 1) / sp.pi * sp.Integral(
        sp.Abs(mu_p) * psi_out, (mu_p, 0, 1)
    )
    return sp.Eq(psi_in, rhs)


# ─────────────────────────────────────────────────────────────────────────
# Self-test
# ─────────────────────────────────────────────────────────────────────────
def _print_summary():
    """Print every symbolic form this module defines (for self-test)."""
    print("=" * 72)
    print("Canonical first-order transport equation (1-D Cartesian)")
    print("=" * 72)
    print(transport_equation())
    print()
    print("Scalar flux definition:")
    print(scalar_flux_definition())
    print()
    print("Streaming operator T·ψ (cartesian):")
    print(streaming_operator("cartesian"))
    print()
    print("Sources:")
    print("  scattering :", scattering_source_isotropic())
    print("  fission    :", fission_source_isotropic())
    print("  external   :", external_source_isotropic())
    print()
    print("Multi-group scattering RHS, 2 groups:")
    rhs_s, _, _ = multigroup_scattering_source(2)
    for g_idx, expr in enumerate(rhs_s):
        print(f"  g={g_idx}: {expr}")
    print()
    print("Multi-group fission RHS, 2 groups:")
    rhs_f, _, _, _ = multigroup_fission_source(2)
    for g_idx, expr in enumerate(rhs_f):
        print(f"  g={g_idx}: {expr}")
    print()
    print("Boundary conditions:")
    r_b = sp.Symbol("r_b", positive=True)
    print("  vacuum (mu>0 inward):", vacuum_bc(r_b, +1))
    print("  reflective          :", reflective_bc(r_b))
    print("  white (isotropic)   :", white_bc(r_b))
    print("=" * 72)


if __name__ == "__main__":
    _print_summary()
