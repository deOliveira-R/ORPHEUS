r"""Paired symbolic-vs-closed-form contract test for the specular BC.

Math-origin pattern (mirroring :func:`orpheus.derivations.sn_balance.
derive_cumprod_recurrence`): the SymPy derivation in
:mod:`orpheus.derivations.peierls_specular` is the **source of truth**
for the rank-:math:`N` specular reflection matrix
:math:`R_{\rm spec} = \tfrac{1}{2}\,M^{-1}`. This test consumes the
SymPy origin as a contract and pins three identities:

1. **Closed-form M tridiagonal == direct symbolic integration of M_nm**
   (verifies the analytical tridiagonal closed form derived from the
   :math:`x\,P_n(x)` Bonnet recurrence).
2. **The :math:`2\,M\,R_{\rm spec} = I` partial-current contract**
   (verifies the construction enforces :math:`J^{-}_m = J^{+}_m` for
   :math:`m = 0,\ldots,N-1`).
3. **Production numpy** :func:`reflection_specular` **matches** the
   lambdified symbolic :math:`R_{\rm spec}` to numerical precision
   (the only discrepancy is float64 inversion roundoff vs exact
   rational arithmetic — the symbolic origin pins the production
   implementation bit-for-bit modulo this floor).

This is the rank-N specular companion to the L0 sn_balance derivation
chain. End-to-end :math:`k_{\rm eff}` ladders for the BC closure live
in :file:`test_peierls_specular_bc.py` (foundation level — they exercise
the full pipeline rather than the primitive contract).

The verified equation is :eq:`peierls-specular-bc-defn` (specular BC
definition, rank-N projection).
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.peierls_geometry import reflection_specular
from orpheus.derivations.peierls_specular import (
    build_M_closed_form,
    build_M_symbolic,
    build_R_specular_symbolic,
)


@pytest.mark.l1
@pytest.mark.verifies("peierls-specular-bc-defn")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_M_closed_form_matches_direct_symbolic_integration(n):
    r"""The closed-form tridiagonal :math:`M` from the Bonnet
    :math:`x\,P_n(x)` recurrence equals the direct symbolic
    integration :math:`M_{nm} = \int_0^1 \mu\,\tilde P_n\,\tilde P_m\,
    \mathrm d\mu` at every rank :math:`N = 1, \ldots, 6`.
    """
    M_sym = build_M_symbolic(n)
    M_closed = build_M_closed_form(n)
    diff = sp.simplify(M_sym - M_closed)
    assert diff == sp.zeros(n, n), (
        f"closed form != direct integration at N={n}: diff={diff}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-specular-bc-defn")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_R_specular_satisfies_2MR_eq_I(n):
    r"""The partial-current contract :math:`2\,M\,R_{\rm spec} = I`
    is the algebraic statement of the rank-:math:`N` specular condition
    :math:`J^{-}_m = J^{+}_m` for :math:`m = 0,\ldots,N-1`. This is
    the symbolic origin's load-bearing identity — every consumer of
    :func:`reflection_specular` relies on it.
    """
    M = build_M_symbolic(n)
    R = build_R_specular_symbolic(n)
    diff = sp.simplify(2 * M * R - sp.eye(n))
    assert diff == sp.zeros(n, n), (
        f"2 M R != I at N={n}: residual={diff}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-specular-bc-defn")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_reflection_specular_matches_symbolic(n):
    r"""Production numpy :func:`reflection_specular(n)` matches the
    lambdified symbolic :math:`R_{\rm spec}` to within float64
    inversion roundoff. The production path inverts a float64 :math:`M`;
    the symbolic path inverts an exact rational matrix. Their diff at
    :math:`N \le 6` is bounded by :math:`O(\kappa(M)\,\epsilon_{\rm
    mach})` — about :math:`1.5 \times 10^{-14}` at :math:`N = 6`.
    """
    R_sym = build_R_specular_symbolic(n)
    R_exact = np.array(R_sym.tolist(), dtype=float)
    R_production = reflection_specular(n)
    np.testing.assert_allclose(
        R_production, R_exact, rtol=1e-12, atol=1e-13,
        err_msg=(
            f"production reflection_specular(n={n}) drift from "
            f"symbolic R_spec exceeds float64 inversion floor"
        ),
    )
