"""The source-iteration count on an all-reflective pure absorber (#340).

With no scattering and no leakage, the only coupling left between sweeps is
the reflective boundary and the only damping is absorption.  The slow mode is
the diamond-difference face sawtooth (the ``-1`` eigenvalue of
:eq:`dd-face-transmission-spectrum`): it has zero cell average, so the
collision term barely sees it, and it decays only through the absorption
suffered around the reflective loop (the theory page's
``sn-boundary-gs-rate-regime`` section, whose control row of 1631 sweeps this
module reproduces).  Two properties of the METHOD follow, and one claim of the
budget law rests on them:

* the count grows by roughly an order per added dimension (d=1 < d=2 < d=3);
* the count is absorption-limited, ``Sigma_t x n_inner`` roughly invariant;
* the derived default budget (:func:`~orpheus.numerics.convergence.default_iteration_budget`)
  converges the d=3 case that opened #340.

These three are a RECORD of how source iteration behaves, not a THEOREM: a
better splitting or an accelerator would legitimately change them, and each
test's failure message says to re-derive the budget guidance when that
happens.

Fixture: two groups, ``Sigma_t = (s, 2 s)``, ``Sigma_s = 0``, a flat
isotropic source ``Q = (1.0, 0.5) / W``, level-symmetric S4, every face
reflective, ``inner_tol = 1e-13``, the budget left to the entry's default.
The exact answer is the flat field ``psi = Q_g / (W Sigma_t,g)`` on every
ordinate; each solve is checked against it before its count is read.  The
fixture nulls scattering, leakage and spatial heterogeneity on purpose: those
are the channels that damp the slow mode, and the claim is about the case
where none of them does.

Origin: ``derivations/diagnostics/diag_d3_absorber_02_si_rate_scaling.py``
(retired, R19), the probe of the 2026-08-08 #340 budget study, which no
collected test ran while a theory page cited it as the pin of the 1631
control.

Runtime: `[M]` 2026-09-22, 14.8 s for the module (five solves, the
d=3 ``s = 0.8`` solve shared by two tests through a module cache).
"""

from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.convergence import default_iteration_budget
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

from tests.gates.numerics.test_default_iteration_budget import _MEASURED_D3_ABSORBER

pytestmark = [
    pytest.mark.l1,
    # Every face reflective closes the DD face mode into an exact kernel;
    # the entry gauge-fixes the TRACE and says so.  The bulk read here is
    # unaffected (sn-loss-kernel-gauge), so the notice is expected.
    pytest.mark.filterwarnings(
        "ignore::orpheus.sn.operators.loss_kernel_gauge.GaugeFreedomWarning"
    ),
]

_INNER_TOL = 1e-13
_Q_G = np.array([1.0, 0.5])
#: SAFETY x the solver's own convergence tolerance: the bound on the
#: converged answer's distance from the flat field.
_FLAT_FIELD_TOL = 10.0 * _INNER_TOL


def _absorber(sig: tuple[float, float]):
    s = np.asarray(sig, float)
    return make_mixture(
        sig_t=s, sig_c=s, sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def _mesh(d: int):
    if d == 1:
        return Mesh1D(edges=np.linspace(0.0, 1.0, 4),
                      mat_ids=np.zeros(3, dtype=int),
                      coord=CoordSystem.CARTESIAN), (3,)
    if d == 2:
        return Mesh2D(edges_x=np.linspace(0.0, 1.0, 4),
                      edges_y=np.linspace(0.0, 2.0, 5),
                      mat_map=np.zeros((3, 4), dtype=int)), (3, 4)
    return tuple(
        AxisMesh(edges=np.linspace(0.0, ext, n + 1), bc_low=None, bc_high=None)
        for ext, n in zip((1.0, 2.0, 3.0), (3, 4, 5))
    ), (3, 4, 5)


@functools.lru_cache(maxsize=None)
def _solve(d: int, s: float):
    """``(record, max relative distance from the flat field)``, default budget."""
    sig = (s, 2.0 * s)
    quad = Quadrature.level_symmetric(sn_order=4)
    mesh, shape = _mesh(d)
    W = float(np.sum(quad.weights))
    q = np.broadcast_to(
        (_Q_G / W).reshape(1, 2, *([1] * len(shape))), (quad.N, 2, *shape),
    ).copy()
    sol = solve_sn_fixed_source(
        {0: _absorber(sig)}, mesh, quad, external_source=q,
        boundary_condition="reflective", inner_tol=_INNER_TOL,
    )
    psi = np.asarray(sol.angular_flux.interior.values)
    flat = _Q_G / (W * np.asarray(sig))
    err = max(float(np.max(np.abs(psi[:, g] - flat[g]) / flat[g])) for g in range(2))
    return sol.record, err


def _converged_count(d: int, s: float) -> int:
    record, err = _solve(d, s)
    assert record.converged, f"d={d}, Sigma_t={s}: {record.report()}"
    assert err < _FLAT_FIELD_TOL, (
        f"d={d}, Sigma_t={s}: the converged answer is {err:.3e} from the flat "
        f"field Q/(W Sigma_t) (bound {_FLAT_FIELD_TOL:.1e}); the count of a wrong "
        f"answer means nothing"
    )
    return record.n_iterations


@pytest.mark.rests_on(
    "tests/gates/sn/solve/test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact",
)
def test_the_count_grows_with_dimension() -> None:
    """d=1 < d=2 < d=3 sweeps to the same tolerance on the same material.

    What it catches: a change to source iteration's reflective coupling
    (a schedule, a splitting, an accelerator, a boundary fold) that changes
    how the slow mode's decay scales with the number of reflective axis
    pairs.  `[M]` 2026-09-22: 33, 259 and 1632 sweeps.  The claim is the
    ORDERING with at least a factor of 2 per added dimension (the measured
    factors are 7.8 and 6.3), which is what makes a dimension-blind fixed
    cap wrong.
    """
    n1, n2, n3 = (_converged_count(d, 0.8) for d in (1, 2, 3))
    assert 2 * n1 < n2 < n3 and 2 * n2 < n3, (
        f"reflective-absorber SI counts d=1,2,3 = {n1}, {n2}, {n3}: no longer "
        f"growing by a factor >= 2 per dimension.  Did the iteration change "
        f"(an accelerator, a new splitting)?  Re-derive this gate and the "
        f"default_iteration_budget guidance with it."
    )


@pytest.mark.rests_on(
    "tests/gates/sn/solve/test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact",
)
def test_the_count_is_absorption_limited() -> None:
    """``Sigma_t x n_inner`` is invariant: absorption is the only damping.

    What it catches: a change that damps the slow mode by a mechanism other
    than absorption (a leak in the reflective law, an accelerator, a
    splitting that folds the face sawtooth), after which every budget
    heuristic derived from the absorption-limited rate is void.

    The discriminator is derived from the two hypotheses.  Absorption-limited:
    ``n ~ 1/Sigma_t``, so the products are equal.  Not absorption-limited:
    ``n`` independent of ``Sigma_t``, so the products spread by the full
    ``Sigma_t`` range, 3.2 / 0.8 = 4.  The gate is the log-midpoint,
    ``spread < sqrt(4) = 2``.  `[M]` 2026-09-22: counts 1632, 851, 438 at
    ``Sigma_t = 0.8, 1.6, 3.2``; products 1306, 1362, 1402; spread 1.07.
    """
    sigmas = (0.8, 1.6, 3.2)
    counts = [_converged_count(3, s) for s in sigmas]
    assert counts[0] > counts[1] > counts[2], (
        f"counts {counts} at Sigma_t = {sigmas} do not fall as Sigma_t rises"
    )
    products = np.array([s * n for s, n in zip(sigmas, counts)])
    spread = float(products.max() / products.min())
    bound = math.sqrt(max(sigmas) / min(sigmas))
    assert spread < bound, (
        f"Sigma_t x n_inner = {products} spreads {spread:.3f} >= {bound:.3f}: "
        f"the reflective-SI slow mode is no longer absorption-limited"
    )


@pytest.mark.rests_on(
    "tests/gates/numerics/test_default_iteration_budget.py::TestTheBudgetTracksItsTolerance::test_it_covers_the_configuration_that_opened_the_issue",
    "tests/gates/sn/solve/test_convergence_contract.py::TestEveryEntryDerivesItsBudget::test_the_solver_resolves_none_to_the_derived_budget",
    "tests/gates/sn/solve/test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact",
)
def test_the_derived_budget_converges_the_case_that_opened_340() -> None:
    """The d=3 case converges inside the DEFAULT budget, and the budget
    table's count for it is still the one production takes.

    What it catches, two things:

    * a default budget that no longer covers the founding case (the
      retired constant 1000 against the 1631 sweeps this case needs); the
      solve then exits truncated and ``converged`` is False;
    * a stale ``_MEASURED_D3_ABSORBER`` table in
      ``tests/gates/numerics/test_default_iteration_budget.py``.  That table is a
      solve-free RECORD ("reproduce by solving at max_inner=6000"), so a
      sweep change (#337 moved 1369 to 1631) silently leaves it describing
      another tree.  This leg re-measures its ``1e-13`` row on every run and
      reds when production's count moves, so that all four rows are
      re-measured together.

    The table's ``needed`` is ``len(record.trajectory)``, the number of
    residual readings, which is ``n_iterations - 1``.
    """
    record, _ = _solve(3, 0.8)
    _converged_count(3, 0.8)
    assert record.budget.limit == default_iteration_budget(_INNER_TOL), record.report()
    assert not record.exhausted_budget, record.report()
    tabled = dict(_MEASURED_D3_ABSORBER)[_INNER_TOL]
    needed = len(record.trajectory)
    assert needed == tabled, (
        f"production now needs {needed} residual readings at inner_tol="
        f"{_INNER_TOL:g}; _MEASURED_D3_ABSORBER in "
        f"tests/gates/numerics/test_default_iteration_budget.py records {tabled}.  "
        f"Re-measure all four of its rows (the same case, max_inner=6000)."
    )
