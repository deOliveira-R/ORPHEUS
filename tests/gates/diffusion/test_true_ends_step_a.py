r"""#707 Step A, the diffusion rows: the loss is V → V*, and k integrates its output ONCE.

Spec: ``scratch/true_ends/stepA/gate_spec.md`` rows DF-1 (RED today) and DF-2 (GREEN today, must stay green: the census
defect ``diffusion/solver.py:431`` is a SPELLING that becomes a double integration the moment the loss returns
integrated values; this row is the regression guard that catches it).
"""

from __future__ import annotations

from functools import cache

import numpy as np
import pytest

from tests.gates.cross_method import _source_chart_cases as cc

pytestmark = pytest.mark.foundation


@cache
def _solver():
    return cc.diffusion_problem()


@cache
def _anchor() -> dict[str, np.ndarray]:
    with np.load(cc.ANCHOR) as data:
        return {k: data[k] for k in data.files}


@pytest.mark.parametrize("label", ["loss", "fission"])
def test_df1_flux_to_source_ends(label: str) -> None:
    """DF-1 RED. The diffusion loss and fission run V → V* (census verdict c: ``LeakageOperator``,
    ``DiffusionBoundaryOperator``, ``MultiplicationOperator`` declare one space at both ends today)."""
    op = getattr(_solver(), label)
    assert op.codomain != op.domain, f"diffusion {label}: codomain == domain"


@pytest.mark.rests_on("tests/gates/cross_method/test_source_chart_law.py::test_diffusion_output_is_integrated")
def test_df2_k_integrates_the_loss_once() -> None:
    r"""DF-2 GREEN, must stay green. At a NON-eigen state, ``compute_keff(φ) = P(φ) / ⟨1, (Aφ)_bulk⟩``, the oracle
    built from the anchored density output of the loss, ``Σ_i V_i Σ_g (Aφ)_{g,i}``, and the flux-side production
    ``Σ_i V_i Σ_g νΣf_g φ_{g,i}`` from raw cross sections. A carve that stores the loss integrated and keeps
    ``_integration.evaluate`` at ``solver.py:431`` multiplies by V twice and reds this row by the ratio
    ``Σ V² (Aφ) / Σ V (Aφ)``.
    """
    solver = _solver()
    a = _anchor()
    phi_flat = a["diffusion/state"]
    member = type(solver.template).from_flat(phi_flat, solver.template)
    V = cc.diffusion_volumes(solver)
    mats = cc._mixtures()
    mat_of_cell = np.asarray(cc._slab_mesh().mat_ids, dtype=int)  # per cell, read from the mesh
    nusf = np.stack([np.asarray(mats[m].SigP, float) for m in mat_of_cell], axis=1)  # (ng, nx)
    phi = np.asarray(member.interior.values, float)
    production = float(np.sum(V[None, :] * nusf * phi))
    loss = float(np.sum(V[None, :] * a["diffusion/loss/bulk"]))
    got = solver.compute_keff(phi_flat)
    assert abs(got - production / loss) <= 1e-13 * abs(production / loss), (got, production / loss)
