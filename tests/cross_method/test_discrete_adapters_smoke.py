r"""Foundation smoke tests for the discrete-mesh Protocol adapters.

This file is the spec § "File 6" gate for the
:mod:`tests.cross_method.discrete_adapters` scaffold. It pins three
properties of :class:`CPMeshAdapter` and :class:`SNMeshAdapter`:

1. **Constructibility** — both adapters can be built from a
   ``(materials, geometry_spec, discretization)`` triple without
   crashing.
2. **Non-crashing solve** — both adapters' :meth:`solve_critical`
   returns a :class:`CriticalSolution` with a positive eigenvalue
   on a basic case (the "smoke" check, NOT a correctness check).
3. **Protocol conformance** — both adapters satisfy
   ``isinstance(x, TransportSolver)``.

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants on the adapter
scaffold, not equation labels.
"""
from __future__ import annotations

import contextlib
import io

import numpy as np
import pytest

from orpheus.derivations.common.discretization_spec import DiscretizationSpec
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.common.solver_protocol import TransportSolver
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.mesh import BC

from .discrete_adapters import CPMeshAdapter, SNMeshAdapter


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


def _smoke_sphere_spec(R_cm: float = 4.0) -> GeometrySpec:
    return GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=0.5 * R_cm,
        critical_dimension_cm=R_cm,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )


def test_cp_mesh_adapter_constructible():
    """CPMeshAdapter.from_problem builds without crash."""
    adapter = CPMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec(),
        discretization=DiscretizationSpec(n_cells=20, n_chord_quad=32),
    )
    assert adapter.method_name == "cp"
    assert adapter.discretization.n_cells == 20
    assert adapter.discretization.n_chord_quad == 32


def test_sn_mesh_adapter_constructible():
    """SNMeshAdapter.from_problem builds without crash."""
    adapter = SNMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec(),
        discretization=DiscretizationSpec(n_cells=20, n_angular=8),
    )
    assert adapter.method_name == "sn"
    assert adapter.discretization.n_cells == 20
    assert adapter.discretization.n_angular == 8


def test_cp_mesh_adapter_solve_critical_runs():
    """CPMeshAdapter.solve_critical returns a positive eigenvalue (smoke)."""
    adapter = CPMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec(R_cm=8.0),
        discretization=DiscretizationSpec(n_cells=20, n_chord_quad=16),
    )
    with contextlib.redirect_stdout(io.StringIO()):
        sol = adapter.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue > 0.0
    assert sol.parameter_value == 8.0
    assert sol.parameter_kind == "domain_extent_cm"
    assert "raw_result" in sol.metadata


def test_sn_mesh_adapter_solve_critical_runs():
    """SNMeshAdapter.solve_critical returns a positive eigenvalue (smoke).

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
    adapter = SNMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=spec,
        discretization=DiscretizationSpec(n_cells=20, n_angular=8),
    )
    sol = adapter.solve_critical()
    assert isinstance(sol, CriticalSolution)
    assert sol.eigenvalue_kind == "k_eff"
    assert sol.eigenvalue > 0.0


def test_cp_mesh_adapter_satisfies_protocol():
    """CPMeshAdapter satisfies the TransportSolver Protocol."""
    adapter = CPMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec(),
    )
    assert isinstance(adapter, TransportSolver)


def test_sn_mesh_adapter_satisfies_protocol():
    """SNMeshAdapter satisfies the TransportSolver Protocol."""
    adapter = SNMeshAdapter.from_problem(
        materials={0: _smoke_mixture()},
        geometry_spec=_smoke_sphere_spec(),
    )
    assert isinstance(adapter, TransportSolver)


@pytest.mark.skip(
    reason="L4 cross-check pending follow-up dispatch — see "
           "spec § 'Deliberately NOT verified'."
)
def test_cross_method_discrete_vs_reference_TODO():
    r"""Placeholder for the L4 cross-check between discrete and reference.

    The Protocol enables a uniform comparator: feed the same
    ``(materials, geometry_spec)`` to ``Billiard.from_problem`` (or
    ``MomentSpace.from_problem``) for the L1-truth-backed continuous
    reference, and to ``CPMeshAdapter.from_problem`` (or
    ``SNMeshAdapter.from_problem``) for the production discrete
    solver. Compare ``solve_critical().eigenvalue``.
    """
    pass
