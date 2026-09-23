"""A σ-variant hub reaches the diffusion solver (C3a, 2026-09-13).

The total cross section is a DATUM of the Problem (``MaterialMesh.sigma_t_cell``)
and the diffusion solver reads the hub's ONE field (``mesh.mat_xs``), so a
``with_cross_sections`` Problem changes the removal term of the diffusion
eigenvalue.  ``[M]`` (test-architect C3 delta §A.6, probe ``p11``): a ×3 override
moves this gate's own fixture ``keff 0.930946184 → 0.030008461`` (``D``
following the datum; ``0.029084534`` with ``D`` frozen — the C3a state).

RULED 2026-09-13 (fork 4 (a)): ``D`` FOLLOWS the overridden σ_t — ``[M]``
``Mixture.diffusion_coefficient = 1/(3·Σ_tr)`` with ``Σ_tr = Σ_t − Σ_s1``-outflow
is always DERIVED, so the field's per-cell view derives ``D`` from the hub's
datum and the material's P1 outflow (C3b-2): bit-identical on a non-overridden
hub, consistent on a σ-variant — the removal AND the leakage terms read one σ_t.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion import DiffusionMesh, solve_diffusion_1d
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _keff(hub: DiffusionMesh) -> float:
    """The entry's own drive (``solve_diffusion_1d``), over a HUB instead of raw data."""
    from orpheus.diffusion import DiffusionSolver
    from orpheus.numerics.eigenvalue import power_iteration

    solver = DiffusionSolver(hub, keff_tol=1e-10, flux_tol=1e-9)
    return float(power_iteration(solver, max_iter=1000, budget_name="max_outer").keff)


def test_a_sigma_variant_hub_reaches_the_diffusion_removal_term() -> None:
    mats = {0: get_mixture("A", "2g")}
    mesh = Mesh1D(edges=np.linspace(0.0, 4.0, 9), mat_ids=np.zeros(8, dtype=int),
                  bc_left=BC("reflective"), bc_right=BC("vacuum"))
    hub = DiffusionMesh(mesh, mats)
    hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
    _require(hub_b != hub and hub_b.same_phase_space(hub), "a σ-variant: another Problem, the same phase space")
    _require(np.array_equal(hub_b.mat_xs.total_cross_section, 3.0 * hub.mat_xs.total_cross_section),
             "the diffusion hub's field reads its own datum")
    D, D_b = hub.mat_xs.diffusion_coefficient, hub_b.mat_xs.diffusion_coefficient
    p1 = np.asarray(mats[0].p1_outflow, dtype=float).reshape(-1, 1)
    _require(np.array_equal(D, 1.0 / (3.0 * (hub.sigma_t_cell - p1))),
             "D derives from the hub's datum: 1/(3(σ_t − Σ_s1-outflow)) per cell")
    _require(np.array_equal(D_b, 1.0 / (3.0 * (3.0 * hub.sigma_t_cell - p1))),
             "…and FOLLOWS a σ-variant's datum (fork 4 (a))")
    _require(np.array_equal(D, np.asarray(mats[0].diffusion_coefficient, dtype=float).reshape(-1, 1) * np.ones_like(D)),
             "bit-identical to the material's own D on a non-overridden hub")
    k, k_b = _keff(hub), _keff(hub_b)
    _require(abs(k - float(solve_diffusion_1d(mats, mesh).keff)) < 1e-9,
             "the hub-driven solve reproduces the raw-data entry (the fixture is honest)")
    _require(abs(k_b - k) / abs(k) > 1e-2, f"the σ-variant must move the diffusion keff (k={k:.6f}, k_b={k_b:.6f})")
