r"""The measured ladders of the trajectory-resolvent reference (the old family), and the command that re-measures each.

Since P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md``
the SN rows compare against the characteristic reference, and their tolerances
come from :mod:`tests.gates.derivations._characteristic_ladders` (the SN side's
residuals moved there with them). What remains here describes the OLD
reference only, for its last consumers until step (e) deletes the family:
the corroboration rows
(``tests/gates/derivations/test_characteristic_reference_corroboration.py``,
G1-G3 and G5) and the Garcia 2021 rows
(``tests/gates/derivations/test_peierls_greens_function_garcia2021.py``). A
ladder ESTIMATES a reference's error: it can refute a reference, never certify
one (the user's step-5 ruling of 2026-10-03, #405 P2). The rules that turn a
ladder into a tolerance are :mod:`tests.gates.derivations._ladder_rules`.

Each table names the command that re-measures it::

    python -O -m tests.gates.derivations._trajectory_resolvent_ladders <ladder>

with ``<ladder>`` one of ``garcia``, ``gate4``. The A|B|A rungs
(``SPHERE_3REG_REFERENCE_K``, ``CYLINDER_3REG_REFERENCE_K``) are re-measured by
``python -O -m tests.gates.derivations._characteristic_ladders
old-sphere-reference`` (``old-cylinder-reference``): they build the A|B|A
specification from ``_aba_reference``, which reaches the new reference since
step (d), and this module is on the old side of the corroboration rows, whose
independence leg refuses an old-side module whose imports reach it. The
reference rungs take minutes to an hour each (the cylinder's 128-node
azimuthal rung about 40 minutes on the host). Run from the repository root.

The problem: fuel A | moderator B | fuel A at outer radii 0.5, 1.5, 2.0 cm, 2
groups, reflective at r = R. All values ``[M]`` 2026-09-26 on the ERR-090
reference (one emission-density spline per region).
"""

from __future__ import annotations

import sys

from tests.gates.derivations._ladder_rules import alternating_error, richardson_error

# ── the sphere reference, (n_r, n_mu) -> k ─────────────────────────────────
#: ``old-sphere-reference`` (module docstring). n_traj_quad = 64, tol = 1e-9.
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
#: spline, the reading retired at #405 P2 step 7b.2.3; ``old-sphere-reference`` now re-measures it through the
#: reference's natural extension (``_aba_reference.shape_observables``), not re-run since.
SPHERE_3REG_REFERENCE_SHAPE_STEP = {"n_r": 6.7e-4, "n_mu": 4.8e-4}
#: The sphere's radial order: the per-region ladder at n_mu = 24 has increment
#: ratio 0.34 on n_r = 24, 36, 48, the second-order value (0.35); the n_r
#: step 36 -> 72 at the fixture is extrapolated at that order.
SPHERE_3REG_RADIAL_ORDER = 2

# ── the cylinder reference: no bound ──────────────────────────────────────
#: ``old-cylinder-reference`` (module docstring): (n_r, n_mu_axial, n_phi_az) -> k. At 8 axial nodes
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


# ── the estimate ──────────────────────────────────────────────────────────


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


# ── re-measurement ────────────────────────────────────────────────────────


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
    "garcia": _garcia,
    "gate4": _gate4,
}

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="python -O -m tests.gates.derivations._trajectory_resolvent_ladders",
                                     description="Re-measure one of the old reference's ladders.")
    parser.add_argument("ladder", choices=sorted(_LADDERS))
    _LADDERS[parser.parse_args().ladder]()
