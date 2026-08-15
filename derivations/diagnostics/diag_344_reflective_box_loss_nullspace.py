"""Diagnostic: is ``A = L + C - S - B``'s singularity on an all-reflective
Cartesian box tangential-slot bookkeeping, or real trace underdetermination?

Created by numerics-investigator on 2026-08-14 for GitHub issue #344
(the structural by-product recorded in ``scratch/issue_341_gs_jacobi_mechanism.md``
§1.3 and ``scratch/d3_absorber_diagnosis.md`` §4).

**Answer, measured**: BOTH, and they are separable.  ``ker A`` is the direct
sum of

  T  the tangential trace slots (``Omega.n == 0``), whose dimension equals
     the tangential trace-DOF count exactly, present ONLY for quadratures
     that place a vanishing direction cosine (``product``, ``lebedev``); and
  R  a genuine trace underdetermination on CURRENT-CARRYING ordinates,
     present for ``level_symmetric`` — which has ZERO tangential slots.

On the fixture #341/#344 were measured on (``level_symmetric``) ``T = 0`` and
the whole kernel is ``R``, so the "benign tangential bookkeeping" reading is
refuted there.

If this test catches a real bug, promote to ``tests/sn/`` — the natural home
is ``tests/sn/operators/`` (the loss-operator rank contract) or
``tests/sn/solve/`` (the convergence-functional blindness legs).  These are
CHARACTERIZATION tests of a known structural property, not red gates: they
pin what the operator does today so a change to the boundary closure or the
spatial scheme cannot move it silently.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    _as_sn_mesh,
    _balance_projection,
    _unwindowed_cold_start,
)
from orpheus.transport.mesh.axis import AxisMesh
from orpheus.transport.spatial.linear_discontinuous import (
    LinearDiscontinuous,
)

V, R = BC("vacuum"), BC("reflective")

#: relative singular-value threshold.  Justified by the MEASURED gap: at
#: d=2 / S4 the twelve smallest singular values lie in 7.2e-17 .. 1.2e-14
#: against ``||A||_2 = 27.5`` and the thirteenth is 1.16e-01 — a ratio of
#: 9.5e+12, so any threshold in [1e-13, 1e-2] returns the same rank.
_RANK_RTOL = 1e-10


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def _axes(extents, cells, bcs):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=lo, bc_high=hi)
        for e, n, (lo, hi) in zip(extents, cells, bcs)
    )


def _absorber(ng=2, sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t[:ng], dtype=float)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t.copy(), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=np.zeros((ng, ng)),
    )


def _build(extents, cells, bcs, mixture, quad):
    sn_mesh = _as_sn_mesh(_axes(extents, cells, bcs), quad, {0: mixture})
    solver = SNSolver(sn_mesh, inner_solver="source_iteration")
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    template = _unwindowed_cold_start(sn_mesh, history_depth=0)
    return sn_mesh, system, template


def _apply_A(system):
    """``A = implicit_operator - sum(explicit_gains)`` — the composition the
    SI driver iterates (``M - N``), built from the PRODUCTION splitting."""
    implicit = system.implicit_operator
    gains = list(system.explicit_gains)

    def apply(x):
        out = implicit.apply(x)
        for g in gains:
            out = out - g.apply(x)
        return out

    return apply


def _dense(apply_fn, template):
    n = template.to_flat().size
    out = np.empty((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        out[:, j] = apply_fn(type(template).from_flat(e, template)).to_flat()
        e[j] = 0.0
    return out


def _tangential_dof_count(sn_mesh):
    odn = sn_mesh.angular_trace.omega_dot_n
    layout = sn_mesh.angular_trace.layout
    total = 0
    for f_idx, f in enumerate(layout.faces):
        slot = layout.faces[f]
        per_ord = int(np.prod(slot.shape[1:]))
        total += int(np.sum(np.abs(odn[f_idx]) < 1e-14)) * per_ord
    return total


def _null_basis(A):
    _u, s, vt = np.linalg.svd(A)
    mask = s < _RANK_RTOL * s[0]
    return vt[mask, :].T, s


# --------------------------------------------------------------------------
# 1. the counting law, and the singular-value GAP that justifies the rank
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "sn_order, ng, expected",
    [(2, 1, 2), (2, 2, 4), (4, 1, 6), (4, 2, 12), (6, 2, 24)],
)
def test_d2_all_reflective_loss_operator_kernel_is_ng_times_N_over_4(
    sn_order, ng, expected,
):
    """``dim ker A = n_groups * N / 4`` on a 2-D all-reflective box.

    Mesh- and scattering-independent (measured over cells (2,2)..(6,8) and
    c = 0, 0.5, 0.9); the rank is read off a 9.5e+12 singular-value gap.
    """
    quad = Quadrature.level_symmetric(sn_order=sn_order)
    sn_mesh, system, template = _build(
        (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(ng), quad,
    )
    A = _dense(_apply_A(system), template)
    Vn, s = _null_basis(A)
    dim_ker = Vn.shape[1]
    assert dim_ker == expected, (
        f"dim ker A = {dim_ker} (expected {expected} = ng*N/4 with "
        f"ng={ng}, N={sn_mesh.angular_trace.omega_dot_n.shape[1]}); "
        f"sigma_max={s[0]:.6e}, sigma[-{dim_ker}]={s[-dim_ker]:.6e}, "
        f"sigma[-{dim_ker + 1}]={s[-dim_ker - 1]:.6e}"
    )
    # the gap is what makes the rank threshold non-arbitrary
    gap = s[-dim_ker - 1] / s[-dim_ker]
    assert gap > 1e8, f"no clean rank gap: sigma ratio = {gap:.3e}"


def test_d2_kernel_vanishes_the_moment_one_face_stops_being_reflective():
    """The singularity is a property of the CLOSED reflective box.

    Measured: all-reflective 12; x-vacuum 0; y-vacuum 0; all-vacuum 0;
    even the mixed ``xmin`` reflective / ``xmax`` vacuum row is 0.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    for bcs, expected in (
        ([(R, R), (R, R)], 12),
        ([(V, V), (R, R)], 0),
        ([(R, R), (V, V)], 0),
        ([(V, V), (V, V)], 0),
        ([(R, V), (R, R)], 0),
    ):
        sn_mesh, system, template = _build(
            (1.0, 2.0), (3, 4), bcs, _absorber(2), quad,
        )
        A = _dense(_apply_A(system), template)
        Vn, s = _null_basis(A)
        assert Vn.shape[1] == expected, (
            f"bcs={bcs}: dim ker A = {Vn.shape[1]}, expected {expected} "
            f"(sigma_min/sigma_max = {s[-1] / s[0]:.3e})"
        )


# --------------------------------------------------------------------------
# 2. ⭐ THE #344 DISCRIMINATOR — the kernel is NOT the tangential slots
# --------------------------------------------------------------------------
def test_level_symmetric_has_zero_tangential_slots_yet_a_nontrivial_kernel():
    """The refutation of the "benign tangential bookkeeping" reading.

    ``level_symmetric`` places its cosines on shells ``|mu| >= mu_1 > 0``,
    so ``Omega.n`` never vanishes — yet ``ker A`` is 12-dimensional and
    every null direction is carried by ordinates with ``|Omega.n| > 0``.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh, system, template = _build(
        (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2), quad,
    )
    odn = sn_mesh.angular_trace.omega_dot_n
    n_tan_pairs = int(np.sum(np.abs(odn) < 1e-14))
    assert n_tan_pairs == 0, (
        f"fixture precondition broken: {n_tan_pairs} tangential "
        f"(face, ordinate) pairs; min |Omega.n| = {np.min(np.abs(odn)):.6e}"
    )

    A = _dense(_apply_A(system), template)
    Vn, s = _null_basis(A)
    assert Vn.shape[1] > 0, "fixture no longer singular — re-derive #344"

    nb = template.interior.values.size
    P_diag = np.einsum("ij,ij->i", Vn, Vn)   # basis-independent null mass
    bulk_mass = float(P_diag[:nb].sum())
    assert bulk_mass < 1e-20, (
        f"the kernel is not pure-trace: bulk mass = {bulk_mass:.6e}"
    )

    layout = sn_mesh.angular_trace.layout
    tan_mass = nontan_mass = 0.0
    for f_idx, f in enumerate(layout.faces):
        slot = layout.faces[f]
        per_ord = int(np.prod(slot.shape[1:]))
        for m in range(slot.shape[0]):
            lo = nb + slot.offset + m * per_ord
            w = float(P_diag[lo:lo + per_ord].sum())
            if abs(odn[f_idx, m]) < 1e-14:
                tan_mass += w
            else:
                nontan_mass += w
    assert tan_mass < 1e-20 and nontan_mass == pytest.approx(
        Vn.shape[1], rel=1e-9
    ), (
        f"null mass split: tangential {tan_mass:.6e}, "
        f"non-tangential {nontan_mass:.6f} of dim ker {Vn.shape[1]}"
    )


def test_product_quadrature_kernel_splits_into_tangential_plus_a_remainder():
    """``dim ker A = (#tangential trace DOFs) + R``, exactly, per quadrature.

    Measured at d=2, cells (3,4), ng=2, all reflective:

    ==================  ========  ===========  ===
    quadrature          tan DOFs  dim ker A      R
    ==================  ========  ===========  ===
    level_symmetric(4)         0           12   12
    product(4,4)             224          224    0
    product(8,8)             448          464   16
    lebedev(11)              224          242   18
    ==================  ========  ===========  ===
    """
    expectations = {
        "level_symmetric(4)": (Quadrature.level_symmetric(sn_order=4), 0, 12),
        "product(4,4)": (Quadrature.product(n_mu=4, n_phi=4), 224, 0),
        "lebedev(11)": (Quadrature.lebedev(order=11), 224, 18),
    }
    for name, (quad, exp_tan, exp_rem) in expectations.items():
        sn_mesh, system, template = _build(
            (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2), quad,
        )
        n_tan = _tangential_dof_count(sn_mesh)
        A = _dense(_apply_A(system), template)
        Vn, s = _null_basis(A)
        dim_ker = Vn.shape[1]
        assert n_tan == exp_tan, (
            f"{name}: tangential DOFs {n_tan}, expected {exp_tan}"
        )
        assert dim_ker - n_tan == exp_rem, (
            f"{name}: dim ker A = {dim_ker}, tangential DOFs = {n_tan}, "
            f"remainder = {dim_ker - n_tan} (expected {exp_rem})"
        )


# --------------------------------------------------------------------------
# 3. the blindness of both convergence functionals, WITH a positive control
# --------------------------------------------------------------------------
def test_both_convergence_functionals_are_blind_to_ker_A_but_not_to_a_control():
    """Neither the SI stop nor the balance defect can see a ``ker A`` shift.

    (a) ``||A psi - q|| / ||q||`` is the SI loop's OWN stopping quantity —
        it computes ``||sum_i g_i (psi_{n-1} - psi_n)|| / ||q||``, which
        equals the equation residual under an exact ``M``
        (``orpheus/numerics/iteration.py``, the rho-honest stop).
    (b) ``||R_g(A psi - q)|| / ||R_g(q)||`` is ``_exit_balance_defect``,
        through the production ``_balance_projection``.

    The POSITIVE CONTROL (vv anti-pattern #17) is a NON-null perturbation of
    the SAME flat 2-norm: both functionals must move by many orders.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh, system, template = _build(
        (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2), quad,
    )
    n = template.to_flat().size
    nb = template.interior.values.size
    apply_A = _apply_A(system)
    A = _dense(apply_A, template)
    Vn, _s = _null_basis(A)
    assert Vn.shape[1] > 0

    W = float(sn_mesh.quad.weights.sum())
    sig_t = np.asarray(sn_mesh.materials[0].SigT, dtype=float)
    Q_g = np.array([1.0, 2.0])
    psi_g = Q_g / (W * sig_t)

    layout = sn_mesh.angular_trace.layout
    flat_sol = np.zeros(n)
    interior = np.zeros(template.interior.values.shape)
    q_int = np.zeros(template.interior.values.shape)
    for g in range(2):
        interior[:, g, ...] = psi_g[g]
        q_int[:, g, ...] = Q_g[g] / W
    flat_sol[:nb] = interior.ravel()
    for f, slot in layout.faces.items():
        view = np.zeros(slot.shape)
        for g in range(2):
            view[:, g, ...] = psi_g[g]
        flat_sol[nb + slot.offset:
                 nb + slot.offset + slot.flat_size] = view.ravel()
    q_flat = np.zeros(n)
    q_flat[:nb] = q_int.ravel()

    def fld(flat):
        return type(template).from_flat(flat, template)

    den_bal = float(np.linalg.norm(
        np.asarray(_balance_projection(fld(q_flat), sn_mesh=sn_mesh))))

    def si_stop(flat):
        r = apply_A(fld(flat)).to_flat() - q_flat
        return float(np.linalg.norm(r) / np.linalg.norm(q_flat))

    def balance_defect(flat):
        r = apply_A(fld(flat)).to_flat() - q_flat
        return float(np.linalg.norm(
            np.asarray(_balance_projection(fld(r), sn_mesh=sn_mesh)))) / den_bal

    # the analytic uniform state IS the discrete solution
    assert si_stop(flat_sol) < 1e-13, si_stop(flat_sol)
    assert balance_defect(flat_sol) < 1e-13, balance_defect(flat_sol)

    rng = np.random.default_rng(0)
    coef = rng.standard_normal(Vn.shape[1])
    v = Vn @ (coef / np.linalg.norm(coef))
    u = rng.standard_normal(n)
    u -= Vn @ (Vn.T @ u)
    u /= np.linalg.norm(u)

    def trace_dev(flat):
        worst = 0.0
        for f, slot in layout.faces.items():
            view = flat[nb + slot.offset:
                        nb + slot.offset + slot.flat_size].reshape(slot.shape)
            for g in range(2):
                worst = max(worst, float(
                    np.max(np.abs(view[:, g, ...] / psi_g[g] - 1.0))))
        return worst

    v_pert = v * (0.1126 / trace_dev(flat_sol + v))
    u_pert = u * np.linalg.norm(v_pert)

    dev = trace_dev(flat_sol + v_pert)
    assert dev == pytest.approx(0.1126, rel=1e-6), dev

    # BLIND to the null shift
    assert si_stop(flat_sol + v_pert) < 1e-13, (
        f"the SI stop MOVED on a ker A shift: "
        f"{si_stop(flat_sol + v_pert):.6e} (trace deviation {dev:.4%})"
    )
    assert balance_defect(flat_sol + v_pert) < 1e-13, (
        f"the balance defect MOVED on a ker A shift: "
        f"{balance_defect(flat_sol + v_pert):.6e}"
    )
    # POSITIVE CONTROL — same flat norm, outside ker A: both must move
    assert si_stop(flat_sol + u_pert) > 1e-3, (
        f"positive control did not move the SI stop — the instrument is "
        f"broken: {si_stop(flat_sol + u_pert):.6e}"
    )
    assert balance_defect(flat_sol + u_pert) > 1e-5, (
        f"positive control did not move the balance defect: "
        f"{balance_defect(flat_sol + u_pert):.6e}"
    )


def test_the_trace_metric_SEES_ker_A_even_though_it_annihilates_tangential_rows():
    """The remedy-relevant asymmetry (``vv-principles`` anti-pattern #18).

    A tangential slot carries ``G = |Omega.n| w_n == 0`` EXACTLY, so no
    G-weighted functional can ever observe it.  The #344 kernel is on
    current-carrying ordinates, so ``<v,v>_G > 0`` — a partial-current-norm
    gate CAN see it.  What it does NOT move is the net current.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh, system, template = _build(
        (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2), quad,
    )
    n = template.to_flat().size
    nb = template.interior.values.size
    A = _dense(_apply_A(system), template)
    Vn, _s = _null_basis(A)
    G = sn_mesh.angular_trace.partial_current_metric

    rng = np.random.default_rng(0)
    coef = rng.standard_normal(Vn.shape[1])
    v = Vn @ (coef / np.linalg.norm(coef))
    g_norm = float(np.sqrt(np.sum(G * v[nb:] ** 2)))
    assert g_norm > 1e-2, (
        f"the kernel is metric-annihilated after all: <v,v>_G^(1/2) = "
        f"{g_norm:.6e} — that would make it the tangential class"
    )
    # but it carries no net current, which is why R_g is blind (test above)
    vfield = type(template).from_flat(v, template)
    for f in sn_mesh.angular_trace.layout.faces:
        J = vfield.boundary.net_current(f)
        assert np.max(np.abs(J)) < 1e-14, (
            f"face {f}: the null mode carries net current "
            f"max|J| = {np.max(np.abs(J)):.6e}"
        )

    # CONTRAST: on product(4,4) the tangential slots really are G == 0
    quad_p = Quadrature.product(n_mu=4, n_phi=4)
    m2, _s2, _t2 = _build((1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2),
                          quad_p)
    G2 = m2.angular_trace.partial_current_metric
    odn2 = m2.angular_trace.omega_dot_n
    lay2 = m2.angular_trace.layout
    tan_max = 0.0
    for f_idx, f in enumerate(lay2.faces):
        slot = lay2.faces[f]
        per_ord = int(np.prod(slot.shape[1:]))
        for m in range(slot.shape[0]):
            if abs(odn2[f_idx, m]) < 1e-14:
                lo = slot.offset + m * per_ord
                tan_max = max(tan_max, float(np.max(np.abs(G2[lo:lo + per_ord]))))
    assert tan_max == 0.0, (
        f"product(4,4) tangential slots should carry G == 0 exactly; "
        f"max G = {tan_max:.6e}"
    )


# --------------------------------------------------------------------------
# 4. PART II — the DISPOSITION gates: the deviation is the SPLITTING's, the
#    trace error is O(h), and the G-orthogonal gauge projection is EXACT
# --------------------------------------------------------------------------
def _uniform_fixture(cells, ng=2, sig_t=(0.8, 1.6), ndim=2):
    """All-reflective box + uniform iso source; returns everything needed to
    compare a solve against the CLOSED FORM ``psi = Q / (sigma_t * sum w)``."""
    from orpheus.sn.solver import _build_fixed_source_rhs

    extents = tuple([1.0, 2.0, 3.0][:ndim])
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh, system, template = _build(
        extents, cells, [(R, R)] * ndim, _absorber(ng, sig_t), quad,
    )
    n = template.to_flat().size
    nb = template.interior.values.size
    W = float(sn_mesh.quad.weights.sum())
    st = np.asarray(sn_mesh.materials[0].SigT, dtype=float)
    psi_g = 1.0 / (W * st)
    q_int = np.zeros(template.interior.values.shape)
    for g in range(ng):
        q_int[:, g, ...] = 1.0 / W
    q_field = _build_fixed_source_rhs(np.asarray(q_int, dtype=float), sn_mesh)
    exact = np.zeros(n)
    it = np.zeros(template.interior.values.shape)
    for g in range(ng):
        it[:, g, ...] = psi_g[g]
    exact[:nb] = it.ravel()
    for f, slot in sn_mesh.angular_trace.layout.faces.items():
        vw = np.zeros(slot.shape)
        for g in range(ng):
            vw[:, g, ...] = psi_g[g]
        exact[nb + slot.offset:
              nb + slot.offset + slot.flat_size] = vw.ravel()
    return sn_mesh, system, template, n, nb, q_field, exact


def _solve(system, sn_mesh, template, q_field, n, schedule, tol=1e-13):
    from orpheus.numerics.iteration import SourceIteration
    from orpheus.sn.solver import _select_si_splitting

    base, gains = _select_si_splitting(
        system.implicit_operator, *system.explicit_gains, sn_mesh, schedule,
    )
    si = SourceIteration(base.inverse(), *gains, max_iter=400000, tol=tol)
    zero = type(template).from_flat(np.zeros(n), template)
    psi, rec = si.solve(q_field, initial_guess=zero)
    assert rec.converged, f"{schedule} did not converge"
    return psi.to_flat()


@pytest.mark.parametrize("cells", [(3, 4), (3, 3), (5, 4)])
def test_the_returned_boundary_trace_depends_on_the_SPLITTING(cells):
    """A splitting cannot change the equation — but it selects which member
    of the SOLUTION MANIFOLD is returned.

    CHARACTERIZATION of the #344 disposition.  Same operator, same source,
    same ZERO cold start: ``jacobi`` returns the analytic uniform trace to
    solver tolerance; ``gauss_seidel`` (the production default for Cartesian
    d>=2) returns a point 6-10 % away, and the difference is in ``ker A``.

    ``vv-principles`` Mode 9 says a splitting must not move the fixed POINT.
    On this configuration there is no fixed point — there is a manifold — so
    that guarantee is void, and this gate records it.
    """
    sn_mesh, system, template, n, nb, q_field, exact = _uniform_fixture(cells)
    apply_A = _apply_A(system)
    devs = {}
    for sched in ("gauss_seidel", "jacobi"):
        flat = _solve(system, sn_mesh, template, q_field, n, sched)
        d = flat - exact
        devs[sched] = (
            float(np.max(np.abs(d[nb:] / exact[nb:]))),
            float(np.max(np.abs(d[:nb] / exact[:nb]))),
            float(np.linalg.norm(
                apply_A(type(template).from_flat(d, template)).to_flat())
                / max(np.linalg.norm(d), 1e-300)),
        )
    gs_trace, gs_bulk, gs_Ad = devs["gauss_seidel"]
    j_trace, j_bulk, _j_Ad = devs["jacobi"]
    assert j_trace < 1e-10, (
        f"jacobi no longer returns the exact trace: {j_trace:.6e}"
    )
    assert gs_trace > 1e-2, (
        f"gauss_seidel trace deviation {gs_trace:.6e} — if this has become "
        f"small the defect was fixed; delete this characterization"
    )
    assert gs_Ad < 1e-9, (
        f"the gauss_seidel deviation is NOT a ker A component "
        f"(||A d||/||d|| = {gs_Ad:.3e}) — the mechanism changed"
    )
    # the BULK is unaffected under BOTH schedules — the control
    assert gs_bulk < 1e-10 and j_bulk < 1e-10, (gs_bulk, j_bulk)


def test_the_exact_solution_is_the_minimum_G_norm_member_of_the_manifold():
    """The gauge is CANONICAL, and projecting onto it is an EXACT fix.

    ``psi_exact`` is G-orthogonal to ``ker A`` to machine precision, so the
    analytic answer IS the minimum-``||.||_G`` point of the solution manifold;
    removing the returned iterate's ``ker A`` content recovers it.
    """
    cells = (3, 4)
    sn_mesh, system, template, n, nb, q_field, exact = _uniform_fixture(cells)
    A = _dense(_apply_A(system), template)
    Vn, _s = _null_basis(A)
    assert Vn.shape[1] > 0

    w = np.asarray(sn_mesh.quad.weights, dtype=float)
    vols = np.asarray(sn_mesh.volumes, dtype=float)
    Gf = np.empty(n)
    Gb = np.zeros(template.interior.values.shape)
    Gb[...] = (w.reshape((w.size, 1) + (1,) * (Gb.ndim - 2))
               * vols.reshape((1, 1) + template.interior.values.shape[2:]))
    Gf[:nb] = Gb.ravel()
    Gf[nb:] = sn_mesh.angular_trace.partial_current_metric

    S = np.sqrt(Gf)
    Qb, _r = np.linalg.qr(S[:, None] * Vn)

    def P_G(x):
        return (Qb @ (Qb.T @ (S * x))) / S

    frac = float(np.sqrt(np.sum(Gf * P_G(exact) ** 2))
                 / np.sqrt(np.sum(Gf * exact ** 2)))
    assert frac < 1e-12, (
        f"psi_exact is NOT G-orthogonal to ker A ({frac:.6e}) — the "
        f"minimum-norm gauge would no longer be the right representative"
    )

    flat = _solve(system, sn_mesh, template, q_field, n, "gauss_seidel")
    before = float(np.max(np.abs((flat[nb:] - exact[nb:]) / exact[nb:])))
    fixed = flat - P_G(flat)
    after = float(np.max(np.abs((fixed[nb:] - exact[nb:]) / exact[nb:])))
    assert before > 1e-2, before
    assert after < 1e-10, (
        f"the G-orthogonal gauge projection did NOT recover the exact trace: "
        f"{before:.6e} -> {after:.6e}"
    )


def test_the_gauss_seidel_trace_error_is_first_order_in_h():
    """``err * n`` is CONSTANT on the ODD-cell family — the error is O(h).

    The error CONVERGES — so this is not a "converges to the wrong limit"
    defect.  The ladder must stay inside ONE PARITY CLASS: at even ``n_x`` the
    deviation is identically zero (``vv-principles`` #13 — a 4/8/16/32 ladder
    reports nothing).

    ⚠ The COEFFICIENT is fixture-dependent, the ORDER is not.  Measured
    ``err * n`` over ``n = 5, 7, 9, 11, 15, 21, 31``:

    * extents ``(1.0, 1.0)``  ->  ``0.311671`` (8 s.f., 6.2x mesh range)
    * extents ``(1.0, 2.0)``  ->  ``0.293872`` (12 s.f. over n = 5, 7, 9)

    So the assertion below pins the ORDER tightly (the spread) and the
    coefficient only loosely — the order is the claim, the coefficient is a
    property of this box.
    """
    prods = []
    for n_cells in (5, 7, 9):
        sn_mesh, system, template, n, nb, q_field, exact = _uniform_fixture(
            (n_cells, n_cells)
        )
        flat = _solve(system, sn_mesh, template, q_field, n, "gauss_seidel",
                      tol=1e-13)
        dev = float(np.max(np.abs((flat[nb:] - exact[nb:]) / exact[nb:])))
        prods.append(dev * n_cells)
    spread = (max(prods) - min(prods)) / np.mean(prods)
    assert spread < 1e-5, (
        f"err * n is not constant -> the trace error is not O(h): "
        f"{prods} (relative spread {spread:.3e})"
    )
    assert abs(np.mean(prods) - 0.293872) < 1e-4, (
        f"the O(h) coefficient for extents (1.0, 2.0) moved: "
        f"{np.mean(prods):.6f} (was 0.293872)"
    )


# --------------------------------------------------------------------------
# 5. PART III — is boundary Gauss-Seidel a SPLITTING of A, or INCOHERENT?
#    (the ERR-056 question).  The last test is the POSITIVE CONTROL and is
#    what makes the three before it carry information.
# --------------------------------------------------------------------------
_KERNEL_FREE = {
    # label: (extents, cells, bcs, scheme) — four INDEPENDENT ways to remove
    # ker A, so the coherence conclusion is not an artefact of how it was
    # removed (vv-principles L11 structural independence).
    "x-vacuum": ((1.0, 2.0), (3, 4), [(V, V), (R, R)], None),
    "xmin-refl/xmax-vac + y-refl": ((1.0, 2.0), (3, 4), [(R, V), (R, R)], None),
    "LD on the ALL-reflective box": ((1.0, 2.0), (3, 4), [(R, R)] * 2, "LD"),
    "d3 one reflective pair": ((1.0, 2.0, 3.0), (2, 2, 2),
                               [(V, V), (V, V), (R, R)], None),
}


def _scatterer(ng=2, c=0.5, sig_t=(0.8, 1.6)):
    st = np.asarray(sig_t[:ng], dtype=float)
    s = np.zeros((ng, ng))
    for g in range(ng):
        s[g, g] = c * st[g]
    return make_mixture(
        sig_t=st, sig_c=st - s.sum(axis=1), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=s,
    )


def _kernel_free_build(label):
    extents, cells, bcs, scheme_tag = _KERNEL_FREE[label]
    scheme = LinearDiscontinuous() if scheme_tag == "LD" else None
    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh = _as_sn_mesh(_axes(extents, cells, bcs), quad,
                          {0: _scatterer()}, scheme=scheme)
    solver = SNSolver(sn_mesh, inner_solver="source_iteration")
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    template = _unwindowed_cold_start(sn_mesh, history_depth=0)
    return sn_mesh, system, template


def _iso_source(sn_mesh, template):
    from orpheus.sn.solver import _build_fixed_source_rhs

    W = float(sn_mesh.quad.weights.sum())
    q_int = np.zeros(template.interior.values.shape)
    for g in range(template.interior.values.shape[1]):
        q_int[:, g, ...] = (1.0 + g) / W
    return _build_fixed_source_rhs(np.asarray(q_int, dtype=float), sn_mesh)


def _both_schedules(sn_mesh, system, template):
    """Return {schedule: converged flat iterate} from the SAME zero cold start."""
    from orpheus.numerics.iteration import SourceIteration
    from orpheus.sn.solver import _select_si_splitting

    n = template.to_flat().size
    q_field = _iso_source(sn_mesh, template)
    out = {}
    for sched in ("gauss_seidel", "jacobi"):
        M_op, gains = _select_si_splitting(
            system.implicit_operator, *system.explicit_gains, sn_mesh, sched,
        )
        si = SourceIteration(M_op.inverse(), *gains, max_iter=400000,
                             tol=1e-13)
        psi, rec = si.solve(
            q_field,
            initial_guess=type(template).from_flat(np.zeros(n), template),
        )
        assert rec.converged, f"{sched} did not converge"
        out[sched] = psi.to_flat()
    return out


@pytest.mark.parametrize(
    "label",
    ["d2 all-reflective absorber", "d2 all-reflective c=0.9", "x-vacuum"],
)
def test_both_schedules_are_splittings_of_the_SAME_A(label):
    """``M - N == A`` for BOTH schedules — the definition of a splitting.

    Measured BIT-EXACTLY zero on 20 of 20 fixtures (absorber, c=0.5, c=0.9,
    heterogeneous mat_map, LD, d=3, kernel-free variants).  If this ever
    fails, boundary Gauss-Seidel is not a splitting of ``A`` and #344's
    DETERMINISM verdict reverses to a correctness bug.
    """
    from orpheus.sn.solver import _select_si_splitting

    quad = Quadrature.level_symmetric(sn_order=4)
    if label == "x-vacuum":
        sn_mesh, system, template = _kernel_free_build("x-vacuum")
    else:
        mixture = (_absorber(2) if "absorber" in label
                   else _scatterer(2, c=0.9))
        sn_mesh, system, template = _build(
            (1.0, 2.0), (3, 4), [(R, R)] * 2, mixture, quad,
        )
    A = _dense(_apply_A(system), template)
    nA = np.linalg.norm(A)
    for sched in ("jacobi", "gauss_seidel"):
        M_op, gains = _select_si_splitting(
            system.implicit_operator, *system.explicit_gains, sn_mesh, sched,
        )
        M = _dense(lambda x: M_op.apply(x), template)

        def n_apply(x, gains=gains):
            out = None
            for g in gains:
                gx = g.apply(x)
                out = gx if out is None else out + gx
            return out

        N = _dense(n_apply, template)
        rel = float(np.linalg.norm(M - N - A) / nA)
        assert rel < 1e-13, (
            f"{label}/{sched}: ||M - N - A||/||A|| = {rel:.6e} — the "
            f"schedule is NOT a splitting of A"
        )


def test_the_gauss_seidel_inverse_is_exact_on_the_DRIVER_RHS_subspace():
    """``M_GS.inverse()`` is a SUBSPACE inverse — exact where it is used.

    ⛔ It is NOT a full-space inverse: ``||M M^-1 - I||`` is O(1), and the
    defect is EXACTLY the (inflow-row, outflow-column) block.  The SI driver's
    right-hand side ``q + S.psi + B_upper.psi`` has identically ZERO
    outflow-trace content, so on that subspace the inverse is exact to ~1e-15.

    ⚠ HAZARD this gate exists to pin: a **Krylov preconditioner** built from
    ``M_GS.inverse()`` WOULD feed it arbitrary vectors and is outside the
    contract.  ``M_J.inverse()`` is exact in full space.
    """
    from orpheus.sn.solver import _select_si_splitting

    quad = Quadrature.level_symmetric(sn_order=4)
    sn_mesh, system, template = _build(
        (1.0, 2.0), (3, 4), [(R, R)] * 2, _absorber(2), quad,
    )
    n = template.to_flat().size
    nb = template.interior.values.size
    odn = sn_mesh.angular_trace.omega_dot_n
    layout = sn_mesh.angular_trace.layout
    outflow = np.zeros(n, dtype=bool)
    for f_idx, f in enumerate(layout.faces):
        slot = layout.faces[f]
        per = int(np.prod(slot.shape[1:]))
        for k in range(slot.shape[0]):
            if odn[f_idx, k] > 0:
                lo = nb + slot.offset + k * per
                outflow[lo:lo + per] = True

    M_op, gains = _select_si_splitting(
        system.implicit_operator, *system.explicit_gains, sn_mesh,
        "gauss_seidel",
    )
    M = _dense(lambda x: M_op.apply(x), template)

    def n_apply(x, gains=gains):
        out = None
        for g in gains:
            gx = g.apply(x)
            out = gx if out is None else out + gx
        return out

    N = _dense(n_apply, template)
    inv_op = M_op.inverse()
    zero = type(template).from_flat(np.zeros(n), template)
    q = _iso_source(sn_mesh, template).to_flat()

    rng = np.random.default_rng(1)
    for _ in range(3):
        r = q + N @ rng.standard_normal(n)
        assert np.max(np.abs(r[outflow])) == 0.0, (
            "the driver rhs grew outflow-trace content — the subspace "
            "argument no longer holds and this contract must be re-derived"
        )
        y = inv_op.apply(
            type(template).from_flat(r, template), initial_guess=zero,
        ).to_flat()
        rel = float(np.max(np.abs(M @ y - r)) / np.max(np.abs(r)))
        assert rel < 1e-12, (
            f"M_GS^-1 is not an inverse on the driver rhs subspace: {rel:.6e}"
        )


@pytest.mark.parametrize("label", list(_KERNEL_FREE))
def test_kernel_free_configs_give_the_SAME_TRACE_under_both_schedules(label):
    """⭐ THE COHERENCE GATE.  With ``dim ker A == 0`` there is nowhere to hide.

    Four INDEPENDENT ways of removing the kernel.  If the boundary
    Gauss-Seidel schedule were incoherent (ERR-056's failure mode — reflecting
    a face whose trace is inconsistent with the bulk state), at least one of
    these would disagree.  Measured: trace ``<= 1.7e-12``, bulk ``<= 2.2e-13``.

    Teeth verified by ``test_the_err056_first_group_reflect_mutation_...``
    below, which drives this same comparison to 1.0 (trace) and 0.39-0.80
    (bulk).
    """
    sn_mesh, system, template = _kernel_free_build(label)
    A = _dense(_apply_A(system), template)
    Vn, _s = _null_basis(A)
    assert Vn.shape[1] == 0, (
        f"{label} is NOT kernel-free (dim ker A = {Vn.shape[1]}) — it cannot "
        f"serve as the coherence control"
    )
    nb = template.interior.values.size
    sol = _both_schedules(sn_mesh, system, template)
    d = sol["gauss_seidel"] - sol["jacobi"]
    ref = sol["jacobi"]
    trace_rel = float(np.max(np.abs(d[nb:])) / np.max(np.abs(ref[nb:])))
    bulk_rel = float(np.max(np.abs(d[:nb])) / np.max(np.abs(ref[:nb])))
    assert trace_rel < 1e-10 and bulk_rel < 1e-10, (
        f"{label}: the two schedules disagree with NO kernel to hide in — "
        f"trace {trace_rel:.6e}, bulk {bulk_rel:.6e}.  Boundary G-S is "
        f"INCOHERENT and #344 reverses to a correctness bug."
    )


def test_the_err056_first_group_reflect_mutation_reddens_the_coherence_gate():
    """⭐ THE POSITIVE CONTROL (``vv-principles`` #17) for the gate above.

    Reflect each face after the FIRST outflowing octant group instead of the
    LAST — exactly the failure ``SweepSchedule.gauss_seidel``'s deferred-reflect
    rule exists to prevent.  On a kernel-free configuration this MUST move both
    the trace and the BULK; if it does not, the coherence gate is a dud and its
    green carries no information.

    Measured: trace 1.0000e+00 / bulk 7.17e-01 (x-vacuum) against the baseline
    7.81e-13 / 1.67e-13 — twelve orders of dynamic range.
    """
    from orpheus.sn.loss_representation import sweep_schedule as ss

    label = "x-vacuum"
    sn_mesh, system, template = _kernel_free_build(label)
    nb = template.interior.values.size
    base = _both_schedules(sn_mesh, system, template)
    base_trace = float(
        np.max(np.abs(base["gauss_seidel"][nb:] - base["jacobi"][nb:]))
        / np.max(np.abs(base["jacobi"][nb:]))
    )
    assert base_trace < 1e-10, base_trace

    original = ss.SweepSchedule.gauss_seidel.__func__

    def first_group_gauss_seidel(cls, mesh):
        reflective = ss._reflective_faces(mesh)
        ordered, by_label = [], {}
        for entry in mesh.quad.octants:
            sweep = ss._octant_sweep(entry, mesh.ndim)
            if sweep.label not in by_label:
                by_label[sweep.label] = []
                ordered.append(sweep.label)
            by_label[sweep.label].append(sweep)
        first = {}
        for gi, lab in enumerate(ordered):
            for f in ss._outgoing_faces(lab):
                if f in reflective and f not in first:
                    first[f] = gi          # ⛔ the ERR-056 mutation
        by_group = {gi: [] for gi in range(len(ordered))}
        for face, gi in first.items():
            by_group[gi].append(face)
        return cls(
            groups=tuple(
                ss.OctantSweepGroup(
                    sweeps=tuple(by_label[lab]),
                    reflect_faces=tuple(sorted(by_group[gi])),
                )
                for gi, lab in enumerate(ordered)
            ),
            kind="gauss_seidel",
        )

    # the mutation must BITE: the schedule really changes
    shipped = tuple(g.reflect_faces
                    for g in ss.SweepSchedule.gauss_seidel(sn_mesh).groups)
    mutated = tuple(g.reflect_faces
                    for g in first_group_gauss_seidel(ss.SweepSchedule,
                                                      sn_mesh).groups)
    assert shipped != mutated, (
        f"the ERR-056 mutation did not change the schedule ({shipped}) — "
        f"the control is inert"
    )

    ss.SweepSchedule.gauss_seidel = classmethod(first_group_gauss_seidel)
    try:
        bad = _both_schedules(sn_mesh, system, template)
    finally:
        ss.SweepSchedule.gauss_seidel = classmethod(original)

    bad_trace = float(
        np.max(np.abs(bad["gauss_seidel"][nb:] - bad["jacobi"][nb:]))
        / np.max(np.abs(bad["jacobi"][nb:]))
    )
    bad_bulk = float(
        np.max(np.abs(bad["gauss_seidel"][:nb] - bad["jacobi"][:nb]))
        / np.max(np.abs(bad["jacobi"][:nb]))
    )
    assert bad_trace > 1e-2 and bad_bulk > 1e-2, (
        f"the ERR-056 mutation did NOT redden the coherence gate "
        f"(trace {bad_trace:.6e}, bulk {bad_bulk:.6e}) — the gate is a dud "
        f"and every coherence green in this module is uninformative"
    )
    # and the revert must restore the baseline
    again = _both_schedules(sn_mesh, system, template)
    assert np.array_equal(again["gauss_seidel"], base["gauss_seidel"]), (
        "the monkeypatch did not revert cleanly"
    )
