r"""The measured ladders the trajectory-resolvent cross-check's tolerances are derived from, and the command that re-measures each.

Every tolerance of the SN-against-trajectory-resolvent rows
(``tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py``,
``tests/gates/sn/sweep/curvilinear/test_unified_matvec_cylinder.py``,
``tests/gates/sn/verification/analytical/test_l1_standoff_slab_cylinder.py``)
and of the Garcia 2021 rows
(``tests/gates/derivations/test_peierls_greens_function_garcia2021.py``) is
COMPUTED here from a table of measured rungs by the functions at the bottom;
no tolerance is typed. A ladder ESTIMATES a reference's error: it can refute a
reference, never certify one (the user's step-5 ruling of 2026-10-03, #405
P2), so no number here is a bound, and the rows comparing against the
uncertified trajectory resolvent are ``compare_uncertified`` comparisons or
strict xfails (#566, #516). Each table names the command that re-measures it::

    python -O -m tests.gates.derivations._trajectory_resolvent_ladders <ladder>

with ``<ladder>`` one of ``sphere-reference``, ``sphere-sn``,
``cylinder-reference``, ``cylinder-sn``, ``records``, ``garcia``, ``gate4``. The
reference rungs take minutes to an hour each (the cylinder's 128-node
azimuthal rung about 40 minutes on the host). Run from the repository root.

The problem of the SN rows: fuel A | moderator B | fuel A at outer radii 0.5,
1.5, 2.0 cm, 2 groups, reflective at r = R (the ``*_3region`` snapshots of
``tests/gates/sn/regression/_generate_snapshots.py``). All values ``[M]``
2026-09-26 on the ERR-090 reference (one emission-density spline per region).
"""

from __future__ import annotations

import math
import sys

# ── the sphere reference, (n_r, n_mu) -> k ─────────────────────────────────
#: ``sphere-reference``. n_traj_quad = 64, tol = 1e-9.
SPHERE_3REG_REFERENCE_K: dict[tuple[int, int], float] = {
    (24, 24): 1.3837371812114572,
    (36, 24): 1.3819174916441415,
    (48, 24): 1.3812925738847415,
    (72, 24): 1.3813805536790806,
    (36, 48): 1.3806829855617522,
    (36, 96): 1.381169542129444,
    (36, 192): 1.3810415940162106,
    (72, 96): 1.3809323414728516,
}
#: The shape metric (fission-gauged cell averages on the SN mesh, fast group governing) between the fixture
#: rung (36, 96) and each finer rung. ``[M]`` 2026-09-26 on the reference's NODAL φ read through a per-region
#: spline, the reading retired at #405 P2 step 7b.2.3; ``sphere-reference`` now re-measures it through the
#: reference's natural extension (``_aba_reference.shape_observables``), not re-run since.
SPHERE_3REG_REFERENCE_SHAPE_STEP = {"n_r": 6.7e-4, "n_mu": 4.8e-4}
#: The sphere's radial order: the per-region ladder at n_mu = 24 has increment
#: ratio 0.34 on n_r = 24, 36, 48, the second-order value (0.35); the n_r
#: step 36 -> 72 at the fixture is extrapolated at that order.
SPHERE_3REG_RADIAL_ORDER = 2

# ── the SN side, relative residuals at the fixtures ───────────────────────
#: ``sphere-sn``: Gauss-Legendre 32, 40 cells (the fixture) against 160 cells
#: (mesh step) and 64 ordinates at 160 cells (angular step), k relative and
#: shape metric.
SPHERE_3REG_SN_STEPS = {
    "k": {"mesh": 5.5e-6, "angle": 9.7e-6},
    "shape": {"mesh": 3.74e-3, "angle": 3.52e-3},
}
#: ``cylinder-sn``: folded 16x32, 40 cells (the fixture) against 32x64 at 40
#: cells (angular step) and 16x32 at 160 cells (mesh step).
CYLINDER_3REG_SN_STEPS = {
    "k": {"mesh": 1.11e-5, "angle": 1.88e-5},
    "shape": {"mesh": 1.59e-3, "angle": 3.44e-3},
}
#: ``cylinder-sn``: the folded 4x8, 40-cell solve's k against its 32x64
#: limit (1.2310213 against 1.2317490): the unified and standoff rows' SN side.
CYLINDER_3REG_SN_4X8_K_STEP = 5.9e-4

# ── the cylinder reference: no bound ──────────────────────────────────────
#: ``cylinder-reference``: (n_r, n_mu_axial, n_phi_az) -> k. At 8 axial nodes
#: the azimuthal steps 16 -> 32 -> 64 -> 128 are +9.6e-3, +2.3e-3, -1.7e-3
#: (relative 7.8e-3, 1.9e-3, 1.4e-3): not monotone, so no bound follows from
#: them; the radial step 24 -> 36 at the fixture is +8.8e-4 and couples to
#: the azimuth (the kink angles move with the node radius).
CYLINDER_3REG_REFERENCE_K: dict[tuple[int, int, int], float] = {
    (24, 8, 16): 1.2213077635577716,
    (24, 8, 32): 1.2309315590544474,
    (24, 8, 64): 1.23326036,
    (24, 8, 128): 1.23157514,
    (24, 4, 32): 1.2306654201728575,
    (24, 16, 32): 1.231036749830859,
    (24, 32, 32): 1.23103956,
    (36, 16, 32): 1.23211676,
}

# ── Garcia 2021 Case 1 (fixed source, vacuum) ─────────────────────────────
#: ``garcia``: the trajectory resolvent's largest relative change at the 14
#: interior table radii between the fixture (n_r, n_mu) = (48, 24) and its
#: finer rungs, and at the outer surface r = R.
GARCIA_CASE1_RESOLVENT_STEP = {"interior": 1.43e-3, "surface": 1.65e-2}


# ── the rules ─────────────────────────────────────────────────────────────


def richardson_error(step: float, refinement_ratio: float, order: float) -> float:
    r"""The coarse rung's error from one step at a known order.

    With :math:`e(n) = C n^{-p}` and the step :math:`s = e(n) - e(\rho n)`,
    :math:`e(n) = s / (1 - \rho^{-p})`.
    """
    return abs(step) / (1.0 - refinement_ratio ** (-order))


def alternating_error(step: float) -> float:
    """The coarse rung's error when the sequence alternates with shrinking steps: the limit lies within one step."""
    return abs(step)


def ceil_one_significant_figure(x: float) -> float:
    """The smallest one-significant-figure value at least ``x`` (1.5e-3 -> 2e-3, 3.2e-3 -> 4e-3)."""
    magnitude = 10.0 ** math.floor(math.log10(x))
    return math.ceil(round(x / magnitude, 9)) * magnitude


def tolerance_for(sut_residual: float, reference_error: float | None) -> float:
    r"""The tolerance a comparison is held to.

    The smallest one-significant-figure value :math:`T` with
    :math:`T \ge 10\,b` (the verification floor) and
    :math:`T \ge 2\,(e + b)` (room for both errors), for the SUT's residual
    :math:`e` and the reference's error :math:`b`: a derived bound once the
    family is certified, a ladder ESTIMATE until then (the sphere today, which
    makes its rows uncertified comparisons, not verification). With no figure
    for the reference it is assumed at the floor, :math:`b = T/10`, which gives
    :math:`T \ge 2.5\,e`: the tolerance the row will hold once a reference
    certifies a tenth of it.
    """
    if reference_error is None:
        return ceil_one_significant_figure(2.5 * sut_residual)
    return ceil_one_significant_figure(
        max(10.0 * reference_error, 2.0 * (sut_residual + reference_error))
    )


def sphere_3reg_reference_ladder_estimate() -> dict[str, float]:
    """The sphere reference's relative error ESTIMATE at (36, 96): radial step extrapolated at second order, plus the
    alternating mu step. Not a bound: no ladder certifies a reference (the user's step-5 ruling of 2026-10-03), so
    this is the documented provenance of the sphere rows' tolerances, never a certificate (#405 P2 step 7b.2.3;
    named ``sphere_3reg_reference_bound`` until then)."""
    k = SPHERE_3REG_REFERENCE_K
    fixture = k[(36, 96)]
    radial = richardson_error(k[(72, 96)] - fixture, 2.0, SPHERE_3REG_RADIAL_ORDER)
    angular = alternating_error(k[(36, 192)] - fixture)
    shape = SPHERE_3REG_REFERENCE_SHAPE_STEP
    return {
        "k": (radial + angular) / fixture,
        "shape": richardson_error(shape["n_r"], 2.0, SPHERE_3REG_RADIAL_ORDER)
        + alternating_error(shape["n_mu"]),
    }


def sn_residual(steps: dict[str, float]) -> float:
    """The SUT's residual at its fixture: the sum of its per-axis steps."""
    return float(sum(steps.values()))


# ── re-measurement ────────────────────────────────────────────────────────


def _sphere_reference() -> None:
    """Each rung's k, and the shape steps from the fixture, read through the trajectory-resolvent reference."""
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue
    from tests.gates.sn.verification.analytical import _aba_reference as aba
    from tests.gates.sn.verification.analytical import test_phase_c_crosscheck as pc
    observables = aba.shape_observables(pc._sphere_3reg_sn_gl32()[1])
    shapes = {}
    for n_r, n_mu in SPHERE_3REG_REFERENCE_K:
        ref = aba.aba_reference_at(CoordSystem.SPHERICAL, {"n_r": n_r, "n_mu": n_mu, "n_traj_quad": 64})
        print(f"({n_r}, {n_mu}): k = {ref.read(Eigenvalue()).value!r}", flush=True)
        if (n_r, n_mu) in ((36, 96), (72, 96), (36, 192)):
            shapes[(n_r, n_mu)] = {(i, g): ref.read(ratio).value for i, g, ratio in observables}
    fixture = shapes[(36, 96)]
    for other in ((72, 96), (36, 192)):
        scale = max(abs(v) for v in shapes[other].values())
        step = max(abs(fixture[key] - shapes[other][key]) for key in fixture) / scale
        print(f"shape step (36, 96) -> {other}: {step:.3e}")


def _sn(geometry: str) -> None:
    """The SN rungs' k, and the shape steps from the fixture in the gauged cell-average metric."""
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    from tests.gates.sn.verification.analytical import _aba_reference as aba
    if geometry == "sphere":
        rungs = {"fixture": (40, 32), "mesh": (160, 32), "angle": (160, 64)}
        build = lambda n, q: {**generator._sphere_3region("2g", n),
                              "quadrature": Quadrature.gauss_legendre(n_ordinates=q), "max_inner": 2000}
    else:
        rungs = {"fixture": (40, (16, 32)), "angle": (40, (32, 64)), "mesh": (160, (16, 32)), "4x8": (40, (4, 8))}
        build = lambda n, q: {**generator._cylinder_3region("2g", n, "folded_4x8"),
                              "quadrature": Quadrature.folded_product(n_mu=q[0], n_phi=q[1]), "max_inner": 2000}
    observables = aba.shape_observables(build(40, rungs["fixture"][1])["mesh"])
    solved = {}
    for name, (n, q) in rungs.items():
        result = generator.run_case(build(n, q))
        k = result.read(Eigenvalue()).value
        solved[name] = (k, {(i, g): result.read(ratio).value for i, g, ratio in observables})
        print(f"{name} {n} {q}: k = {k!r}", flush=True)
    k0, shape0 = solved["fixture"]
    for name in ("mesh", "angle", "4x8"):
        if name in solved:
            k, shape = solved[name]
            scale = max(abs(v) for v in shape.values())
            print(f"{name} step: k {abs(k - k0) / k:.3e}, shape {max(abs(shape0[key] - shape[key]) for key in shape0) / scale:.3e}")


def _cylinder_reference() -> None:
    """Each rung's k through the trajectory-resolvent reference (the rows' factory, its one spelling of the solve)."""
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue
    from tests.gates.sn.verification.analytical import _aba_reference as aba
    for n_r, n_mu, n_phi in CYLINDER_3REG_REFERENCE_K:
        quadrature = {"n_r": n_r, "n_mu_axial": n_mu, "n_phi_az": n_phi, "n_traj_quad": 64}
        k = aba.aba_reference_at(CoordSystem.CYLINDRICAL, quadrature).read(Eigenvalue()).value
        print(f"({n_r}, {n_mu}, {n_phi}): k = {k!r}", flush=True)


def _records() -> None:
    from orpheus.numerics.observable import Eigenvalue
    from tests.gates.sn.sweep.curvilinear import test_unified_matvec_cylinder as unified
    from tests.gates.sn.verification.analytical import test_l1_standoff_slab_cylinder as standoff
    from tests.gates.sn.verification.analytical import test_phase_c_crosscheck as pc
    print("phase C:", pc._cylinder_readings(), flush=True)
    print("unified k:", repr(unified._unified_cylinder_solution().read(Eigenvalue()).value), flush=True)
    print("standoff sweep k (nx=40):", repr(standoff._solve_cyl_via_sweep(nx=40).read(Eigenvalue()).value), flush=True)


def _garcia() -> None:
    import numpy as np
    from tests.gates.derivations import test_peierls_greens_function_garcia2021 as g
    rungs = ((48, 24), (48, 48), (96, 48), (96, 96))
    flux = {}
    for n_r, n_mu in rungs:
        res = g.solve_case1(n_r=n_r, n_mu=n_mu)
        flux[(n_r, n_mu)] = np.array([g._phi_at(res, r) for r in g.GARCIA_2021_CASE1_R])
        print(f"({n_r}, {n_mu}) relative error against the table:",
              np.abs(flux[(n_r, n_mu)] - g.GARCIA_TO_VARIANT_ALPHA_FACTOR * g.GARCIA_2021_CASE1_PHI)
              / (g.GARCIA_TO_VARIANT_ALPHA_FACTOR * g.GARCIA_2021_CASE1_PHI), flush=True)
    fixture = flux[rungs[0]]
    steps = np.max([np.abs(fixture - flux[r]) / np.abs(flux[r]) for r in rungs[1:]], axis=0)
    print("largest step from (48, 24), interior:", steps[:-1].max(), "surface:", steps[-1])


def _gate4() -> None:
    """The interface jump of ``test_mr_interface_continuity_3region`` along n_r = 24, 36, 48, 72."""
    import numpy as np
    from scipy.interpolate import CubicSpline
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
        solve_greens_function_cylinder_mr,
    )
    radii = np.array([2.0, 3.0, 4.0])
    for n_r in (24, 36, 48, 72):
        res = solve_greens_function_cylinder_mr(
            radii=radii, sigma_t=np.array([2.0, 0.5, 1.5])[:, None],
            sigma_s=np.array([1.8, 0.3, 1.0]).reshape(3, 1, 1),
            nu_sigma_f=np.array([0.1, 0.0, 0.05])[:, None], alpha=1.0,
            n_r=n_r, n_mu_axial=16, n_phi_az=32, n_traj_quad=48, max_iter=300, tol=1e-9,
        )
        region = np.asarray(res.region_at_node)
        phi = res.phi_g[0]
        jumps = []
        for k, r_k in enumerate(radii[:-1]):
            left = CubicSpline(res.r_nodes[region == k], phi[region == k])(r_k)
            right = CubicSpline(res.r_nodes[region == k + 1], phi[region == k + 1])(r_k)
            jumps.append(float(abs(left - right) / max(abs(left), abs(right))))
        print(f"n_r = {n_r}: relative jumps {jumps}", flush=True)


_LADDERS = {
    "sphere-reference": _sphere_reference,
    "sphere-sn": lambda: _sn("sphere"),
    "cylinder-reference": _cylinder_reference,
    "cylinder-sn": lambda: _sn("cylinder"),
    "records": _records,
    "garcia": _garcia,
    "gate4": _gate4,
}

if __name__ == "__main__":
    _LADDERS[sys.argv[1]]()
