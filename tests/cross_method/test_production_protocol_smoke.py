r"""Foundation smoke tests for the production CP / SN Protocol entry points.

This file is the foundation-tier gate for the
:meth:`~orpheus.cp.solver.CPSolver.from_problem` /
:meth:`~orpheus.sn.solver.SNSolver.from_problem` factories that
mint the production discrete-mesh solvers onto the
:class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol. It pins three properties of each factory:

1. **Constructibility** — both classes can be built from a
   ``(materials, geometry_spec, discretization)`` triple without
   crashing.
2. **Non-crashing solve** — both classes' :meth:`solve_critical`
   returns a :class:`CriticalSolution` with a positive eigenvalue
   on a basic case (the "smoke" check, NOT a correctness check).
3. **Protocol conformance** — both classes satisfy
   ``isinstance(x, TransportSolver)``.

History
-------

This file replaced an earlier scaffold-based variant
(``tests.cross_method.discrete_adapters``) that wrapped the
production CP / SN solvers behind ad-hoc ``CPMeshAdapter`` /
``SNMeshAdapter`` facades. Step 4 of the input-cleanup track
(2026-05-04) deleted the scaffold and promoted the production
classes onto the Protocol natively. The smoke tests stayed in
place — they're the same shape contracts, now wired against the
production factory directly.

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants on the
production-class Protocol surface, not equation labels.
"""
from __future__ import annotations

import contextlib
import io

import numpy as np
import pytest

from orpheus.cp.solver import CPSolver
from orpheus.derivations.common.discretization_spec import DiscretizationSpec
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.common.solver_protocol import TransportSolver
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.mesh import BC
from orpheus.sn.solver import SNSolver


pytestmark = [pytest.mark.foundation]


def _smoke_mixture():
    """Pu-like 1G fissile mixture (νσ_f / σ_a > 1)."""
    return make_mixture(
        sig_t=np.array([0.5]),
        sig_c=np.array([0.05]),
        sig_f=np.array([0.05]),
        nu=np.array([2.5]),
        chi=np.array([1.0]),
        sig_s=np.array([[0.40]]),
    )


def _smoke_sphere_spec_cp(R_cm: float = 4.0) -> GeometrySpec:
    """Sphere spec compatible with the CP solver's outer-BC registry.

    CPMesh supports outer BC ∈ {white, vacuum}; reflective is not
    in its registry. The smoke test uses (reflective centre, white
    outer) — the standard CP infinite-lattice convention.
    """
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=0.5 * R_cm,
        critical_dimension_cm=R_cm,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC("white"),
    )


def _smoke_sphere_spec_sn(R_cm: float = 4.0) -> GeometrySpec:
    """Sphere spec compatible with the SN solver's BC registry.

    SN supports {vacuum, reflective}; we use (reflective centre,
    vacuum outer) — the standard SN bare-critical convention.
    """
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=0.5 * R_cm,
        critical_dimension_cm=R_cm,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )


def test_cp_solver_constructible():
    """CPSolver.from_problem builds without crash."""
    cp = CPSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec_cp(),
        discretization=DiscretizationSpec(n_cells=20, n_chord_quad=32),
    )
    assert cp.method_name == "cp"
    assert cp._discretization.n_cells == 20
    assert cp._discretization.n_chord_quad == 32


def test_sn_solver_constructible():
    """SNSolver.from_problem builds without crash."""
    sn = SNSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec_sn(),
        discretization=DiscretizationSpec(n_cells=20, n_angular=8),
    )
    assert sn.method_name == "sn"
    assert sn._discretization.n_cells == 20
    assert sn._discretization.n_angular == 8


def test_cp_solver_solve_critical_runs():
    """CPSolver.solve_critical returns a positive eigenvalue (smoke)."""
    cp = CPSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec_cp(R_cm=8.0),
        discretization=DiscretizationSpec(n_cells=20, n_chord_quad=16),
    )
    with contextlib.redirect_stdout(io.StringIO()):
        sol = cp.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue > 0.0
    assert sol.parameter_value == 8.0
    assert sol.parameter_kind == "domain_extent_cm"
    assert "raw_result" in sol.metadata


def test_sn_solver_solve_critical_runs():
    """SNSolver.solve_critical returns a positive eigenvalue (smoke).

    Closed slab (reflective both sides) → k_eff = k_inf =
    νσ_f / σ_a.
    """
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=2.0,
        critical_dimension_cm=4.0,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    sn = SNSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=spec,
        discretization=DiscretizationSpec(n_cells=20, n_angular=8),
    )
    sol = sn.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue > 0.0


def test_cp_solver_satisfies_protocol():
    """CPSolver satisfies the TransportSolver Protocol."""
    cp = CPSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec_cp(),
    )
    assert isinstance(cp, TransportSolver)


def test_sn_solver_satisfies_protocol():
    """SNSolver satisfies the TransportSolver Protocol."""
    sn = SNSolver.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec_sn(),
    )
    assert isinstance(sn, TransportSolver)


@pytest.mark.skip(
    reason="L4 cross-check pending follow-up dispatch — see "
           "spec § 'Deliberately NOT verified'."
)
def test_cross_method_discrete_vs_reference_TODO():
    r"""Placeholder for the L4 cross-check between discrete and reference.

    The Protocol enables a uniform comparator: feed the same
    ``(materials, geometry_spec)`` to ``Billiard.from_problem`` (or
    ``MomentSpace.from_problem``) for the L1-truth-backed continuous
    reference, and to ``CPSolver.from_problem`` (or
    ``SNSolver.from_problem``) for the production discrete solver.
    Compare ``solve_critical().eigenvalue``.
    """
    pass
