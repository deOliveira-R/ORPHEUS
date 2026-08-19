r"""The cone predicate against PRODUCTION output — the DD negative-flux pair.

Campaign 1 CS3 step 4 (`.claude/plans/space_and_kernel_binding_campaign.md`
§4): cone membership is an element PREDICATE because the shipped
realization does not preserve K — diamond difference can march a positive
source and positive inflow into negative outgoing flux
(``DiscretizationScheme.is_positivity_preserving = False``, honestly).
This module exercises :meth:`~orpheus.numerics.field.Field.cone_violations`
on CONVERGED output of the public entry ``solve_sn_fixed_source``, both ways
(``vv`` #11):

* **positive leg** — optically moderate slab: the converged ψ and φ are in K;
* **negative leg** — the SAME materials / quadrature / BC / cell size
  (``Δx·Σ_t = 100`` in BOTH legs — ``nx`` and ``width`` move together), a
  DEEPER domain: past enough decay, DD's dome overshoot drives entries of ψ
  (and φ) negative. `[M]` frozen at the CS3 verification plan §5.3
  (2026-08-19): ``min ψ = −6.399383e-01`` with 2 of 8 entries negative,
  ``min φ = −8.438399e-01``; the positive sibling ``min ψ = +2.181405e-01``.

Because the per-cell optical thickness is IDENTICAL across the legs, the
negative leg cannot be explained by a materials, quadrature, or cell-size
difference — it is the discretization's marching behaviour in depth, which
is the ruling's own argument for predicate-not-invariant (a ψ≥0 type would
REFUSE this legitimate production output). The cell-size mechanism itself
is documented on the cone chapter's scan tables (``field_algebra.rst``):
``Δx·Σ_t = 1`` stays in K, ``= 2`` leaves it.

**What this module does NOT claim** (the predicate observes; nothing
enforces): production does not keep fields in K, a violation is not handled
or warned, DD is not repaired, and a violating solve is not refused.

Marks ``foundation`` — a software-predicate wiring gate through the public
entry; no physics ``:label:`` claim. Explicit raises / ``np.testing`` only
(``-O`` strips bare ``assert``).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source


pytestmark = pytest.mark.foundation

_SIG_T = 10.0
_C = 0.5


def _mix() -> Mixture:
    z = np.zeros(1)
    return Mixture(
        SigC=np.array([(1.0 - _C) * _SIG_T]),
        SigL=z.copy(), SigF=z.copy(), SigP=z.copy(),
        SigT=np.array([_SIG_T]),
        SigS=[csr_matrix(np.array([[_C * _SIG_T]]))],
        Sig2=csr_matrix((1, 1)), chi=z.copy(),
    )


def _solve(nx: int, width: float):
    mesh = Mesh1D(
        edges=np.linspace(0.0, width, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=2)
    q = np.zeros((quad.N, 1, nx))
    # Per-ordinate source of 100 in cell 0 — the asymmetry that drives the
    # DD dome; magnitude matches the CS3 plan §5.3 freeze so the [M] numbers
    # cross-reference verbatim (min ψ = −6.399383e-01 on the thick leg).
    q[:, 0, 0] = 100.0
    sol = solve_sn_fixed_source({0: _mix()}, mesh, quad, external_source=q)
    history = sol.history
    if history is None:
        raise AssertionError("the solve returned no iteration history")
    np.testing.assert_equal(bool(history.converged), True)
    return sol


def test_converged_moderate_slab_is_in_the_cone() -> None:
    """POSITIVE (vv #11 pairing) — the benign sibling: same materials,
    quadrature, BC, source and CELL SIZE; half the domain depth."""
    sol = _solve(nx=2, width=20.0)
    psi = sol.angular_flux.interior
    if float(np.min(psi.values)) <= 0.0:
        raise AssertionError(
            "activation lost: the benign sibling's ψ is no longer strictly "
            "positive — the fixture pair stopped discriminating."
        )
    if psi.cone_violations() != []:
        raise AssertionError("the predicate rejected a cone member (ψ)")
    if sol.scalar_flux.cone_violations() != []:
        raise AssertionError("the predicate rejected a cone member (φ)")


def test_dd_thick_cell_output_violates_and_the_report_says_where() -> None:
    """NEGATIVE — same cell size (Δx·Σ_t = 100), doubled domain depth: DD
    marches the same positive source into negative ψ entries, the converged
    solve reports success, and the predicate says exactly WHERE."""
    sol = _solve(nx=4, width=40.0)
    psi = sol.angular_flux.interior
    values = np.asarray(psi.values)
    if not float(np.min(values)) < -0.5:
        raise AssertionError(
            "activation lost: the thick-cell fixture no longer produces the "
            "deep negative entries it was frozen with (min ψ ≈ −0.64)."
        )
    got = psi.cone_violations()
    expected = {tuple(int(i) for i in idx) for idx in np.argwhere(values < 0.0)}
    if set(got) != expected:
        raise AssertionError(
            f"the violation report {got} is not exactly the negative-entry "
            f"set {sorted(expected)}"
        )
    if got[0] != tuple(int(i) for i in
                       np.unravel_index(int(values.argmin()), values.shape)):
        raise AssertionError("the worst violation is not reported first")
    if sol.scalar_flux.cone_violations() == []:
        raise AssertionError("φ inherited no violation — the reduction leg "
                             "of the witness went silent")
