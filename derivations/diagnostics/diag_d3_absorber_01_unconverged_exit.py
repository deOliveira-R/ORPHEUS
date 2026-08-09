"""Diagnostic: the d=3 all-reflective pure-absorber fixed source exits
``converged=False`` at the default ``max_inner``, and the gate that reads it
never checks the flag.

Created by numerics-investigator on 2026-08-08.

Root cause of the red in
``tests/sn/solve/test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact``.
The exact uniform field IS the discrete solution (``||A psi - q||/||q|| =
1.06e-15``); the solve simply stops 632 iterations short of it because
source iteration on an ALL-reflective pure absorber needs ~1631 sweeps at
``inner_tol=1e-13`` while ``solve_sn_fixed_source`` defaults to
``max_inner=1000``.

Three legs, each a separate regression gate:

* ``test_d3_absorber_default_max_inner_does_not_converge`` -- the STATE OF
  THE TREE today (strict xfail once ``max_inner`` is raised at the entry or
  the gate; deleting the marker is then the todo signal).
* ``test_d3_absorber_converges_when_given_enough_iterations`` -- the
  positive control: with headroom the answer IS exact.  This is the leg
  that must be promoted; it pins the physics claim the red gate was
  trying to make.
* ``test_d3_absorber_exact_uniform_field_is_the_discrete_solution`` -- the
  reference-independent structural pin: the closed-form field has a
  machine-zero equation residual through the PRODUCTION operator, so the
  discretization is exonerated no matter what the solver does.

If this catches a real bug, promote to ``tests/sn/solve/test_d3_admission.py``
(the same module the red gate lives in).
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    _bare_loss_arm,
    _build_fixed_source_rhs,
    evaluate_residual,
    solve_sn_fixed_source,
)
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.mesh.axis import AxisMesh

# The gate's own fixture, restated locally so the diagnostic is standalone.
_EXTENTS = (1.0, 2.0, 3.0)
_CELLS = (3, 4, 5)
_SIG_T = np.array([0.8, 1.6])
_Q_G = np.array([1.0, 0.5])


def _axes():
    return tuple(
        AxisMesh(edges=np.linspace(0.0, ext, n + 1), bc_low=None, bc_high=None)
        for ext, n in zip(_EXTENTS, _CELLS)
    )


def _absorber():
    return make_mixture(
        sig_t=_SIG_T, sig_c=_SIG_T,
        sig_f=np.zeros(2), nu=np.zeros(2), chi=np.zeros(2),
        sig_s=np.zeros((2, 2)),
    )


def _source(quad):
    W = float(np.sum(quad.weights))
    return np.broadcast_to(
        (_Q_G / W)[None, :, None, None, None], (quad.N, 2, *_CELLS),
    ).copy(), W


def _solve(max_inner, inner_schedule="gauss_seidel"):
    quad = Quadrature.level_symmetric(sn_order=4)
    mix = _absorber()
    q, W = _source(quad)
    sol = solve_sn_fixed_source(
        {0: mix}, _axes(), quad, external_source=q,
        boundary_condition="reflective", inner_tol=1e-13,
        max_inner=max_inner, inner_schedule=inner_schedule,
    )
    return sol, quad, mix, q, W


def _max_rel_err(sol, W):
    psi = np.asarray(sol.angular_flux.interior.values)
    return max(
        float(np.max(np.abs(psi[:, g] - _Q_G[g] / (W * _SIG_T[g]))
                     / (_Q_G[g] / (W * _SIG_T[g]))))
        for g in range(2)
    )


@pytest.mark.l1
def test_d3_absorber_default_max_inner_does_not_converge() -> None:
    """The default ``max_inner=1000`` truncates the solve mid-descent.

    Characterizes the tree as measured 2026-08-08.  The solver is HONEST
    about it (``history.converged is False``); the defect is that a
    best-effort answer is indistinguishable from a converged one at the
    call site.  Flip to ``xfail(strict=True)`` when the entry default is
    raised or the certificate learns to refuse a silent best-effort exit.
    """
    sol, _quad, _mix, _q, W = _solve(max_inner=1000)
    h = sol.history

    np.testing.assert_equal(bool(h.converged), False)
    np.testing.assert_equal(h.n_inner, 999)
    # The running residual is still four orders above the requested tol.
    last = float(list(h.flux_residuals)[-1])
    if not (last > 1e-10):
        pytest.fail(
            f"running residual {last:.3e} is no longer above 1e-10 — the "
            f"iteration count or the rate changed; re-derive this gate."
        )
    # ...and the answer is correspondingly off the closed form.
    err = _max_rel_err(sol, W)
    if not (err > 1e-11):
        pytest.fail(
            f"best-effort error {err:.3e} is now small — the truncation no "
            f"longer bites; this leg has served its purpose, retire it."
        )


@pytest.mark.l1
@pytest.mark.verifies("transport-cartesian")
def test_d3_absorber_converges_when_given_enough_iterations() -> None:
    """THE physics claim: with iteration headroom, psi_n = Q_g/(W Sigma_t,g).

    Pure absorber + all-reflective ⇒ the flat field satisfies every DD face
    difference and every reflective pairing, so the discrete solution IS the
    closed form.  ``max_inner=4000`` clears the measured 1631 sweeps with
    ~2.5x headroom (the count scales as ~1300/Sigma_t — see
    ``diag_d3_absorber_02_si_rate_scaling``).
    """
    sol, quad, mix, _q, W = _solve(max_inner=4000)
    np.testing.assert_equal(bool(sol.history.converged), True)

    psi = np.asarray(sol.angular_flux.interior.values)
    sig_t = np.asarray(mix.SigT)
    for g in range(2):
        np.testing.assert_allclose(
            psi[:, g], _Q_G[g] / (W * sig_t[g]), rtol=1e-10,
            err_msg=f"per-ordinate flat-flux identity, group {g}",
        )
    phi = np.asarray(sol.scalar_flux.values)
    for g in range(2):
        np.testing.assert_allclose(phi[g], _Q_G[g] / sig_t[g], rtol=1e-10)


@pytest.mark.l0
@pytest.mark.verifies("transport-cartesian")
def test_d3_absorber_exact_uniform_field_is_the_discrete_solution() -> None:
    """Solver-independent: the closed form has a MACHINE-ZERO residual.

    Builds the exact uniform state (bulk AND boundary trace) and pushes it
    through the production loss operator ``A = L + C - S - B``.  This
    exonerates the discretization without reference to any iterate, and is
    the measurement that discriminates 'the solve failed to arrive' from
    'the discrete system does not admit the closed form'.

    Control leg: the same evaluation on a DELIBERATELY perturbed uniform
    field must be O(perturbation) — otherwise the residual is blind and the
    machine-zero reading above carries no information (vv anti-pattern #19).
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    mix = _absorber()
    q, W = _source(quad)

    # A throwaway solve only to obtain a production SNMesh + a state template
    # on the SAME instance (the operators enforce mesh identity).
    sol = solve_sn_fixed_source(
        {0: mix}, _axes(), quad, external_source=q,
        boundary_condition="reflective", inner_tol=1e-13, max_inner=1,
    )
    sn = sol.mesh
    solver = SNSolver(sn, inner_solver="source_iteration",
                      scattering_order=0, max_inner=1, inner_tol=1e-13)
    system = build_within_group_system(
        sn, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    loss = _bare_loss_arm(system)
    qc = _build_fixed_source_rhs(q, sn)
    q_norm = float(np.linalg.norm(np.asarray(qc.to_flat())))

    want = np.array([_Q_G[g] / (W * _SIG_T[g]) for g in range(2)])
    tmpl = sol.angular_flux

    def uniform_state(scale=1.0, perturb_group=None, perturb=0.0):
        vals = np.empty_like(np.asarray(tmpl.interior.values))
        for g in range(2):
            vals[:, g] = want[g] * scale
        if perturb_group is not None:
            vals[:, perturb_group] *= 1.0 + perturb
        faces = {}
        for fn in tmpl.boundary.space.face_names:
            fv = np.asarray(tmpl.boundary.face_view(fn))
            arr = np.empty_like(fv)
            for g in range(2):
                arr[:, g] = want[g] * scale
            if perturb_group is not None:
                arr[:, perturb_group] *= 1.0 + perturb
            faces[fn] = arr
        return dataclasses.replace(
            tmpl,
            interior=dataclasses.replace(tmpl.interior, values=vals),
            boundary=AngularBoundaryFlux.from_face_arrays(sn, faces),
        )

    def rel_resid(state):
        r = evaluate_residual(loss, state, qc)
        return float(np.linalg.norm(np.asarray(r.to_flat()))) / q_norm

    exact = rel_resid(uniform_state())
    np.testing.assert_array_less(
        exact, 1e-13,
        err_msg=(f"the closed-form uniform field must satisfy the discrete "
                 f"system to machine precision; got ||A psi - q||/||q|| = "
                 f"{exact:.3e}"),
    )

    # CONTROL: a 1e-6 perturbation must move the residual O(1e-6) relative,
    # i.e. at least 1e3 x the exact reading.  Without this the machine-zero
    # above could be a blind functional.
    perturbed = rel_resid(uniform_state(perturb_group=0, perturb=1e-6))
    if not (perturbed > 1e3 * max(exact, 1e-16)):
        pytest.fail(
            f"residual control failed: perturbing group 0 by 1e-6 moved "
            f"||A psi - q||/||q|| only from {exact:.3e} to {perturbed:.3e} "
            f"— the residual is not loaded on the bulk field."
        )
