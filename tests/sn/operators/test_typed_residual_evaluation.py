r"""Piece 3 (#208 box 7) — typed equation-residual evaluation via from_balance.

The residual mint :meth:`AngularResidual.from_balance` /
:meth:`BoundaryResidual.from_balance` now has a production-reachable consumer:
:func:`orpheus.sn.solver.evaluate_residual` types the within-group balance
defect :math:`r = (L+C-S-B)\psi - q` as a composite
``TimedFullField(bulk=AngularResidual, boundary=BoundaryResidual)`` — a
diagnostic (``balance_map`` / ``boundary_vs_interior_split`` / ``relative_to``)
and the consistent-DSA (`#2`) low-order-correction substrate.

Config: a 1-D slab, 2 groups, fuel|moderator heterogeneous, P1-anisotropic
(Cardinal-6 ≥2G + heterogeneous + anisotropic — the downscatter / redistribution
the per-ordinate balance probe must capture). 1-D never windows, so the
converged iterate ``bulk`` is a full-angular ``AngularFlux`` the operators
consume directly.

Marks: ``foundation`` for the mint-guard + split + relative_to invariants;
``l0`` for the balance-map physics row (it verifies the discrete balance term
``(L+C-S-B)ψ-q``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import (
    SNSolver,
    _within_group_si,
    _within_group_triple,
    boundary_vs_interior_split,
    evaluate_residual,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.residuals import AngularResidual
from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField


def _converged_slab_2g(nx: int = 24, n_ord: int = 8):
    """Converge a 2G het (fuel|moderator) P1 slab via SI; return
    ``(solver, loss_op, q_ext, psi)`` with ``psi`` the full-angular flux."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    mat_ids = np.zeros(nx, dtype=int)
    mat_ids[: nx // 2] = 2  # fuel | moderator split
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1), mat_ids=mat_ids,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    sn_mesh = SNMesh(mesh, quad, {2: fuel, 0: mod})
    solver = SNSolver(sn_mesh, inner_solver="source_iteration", scattering_order=1)
    LC, S, B = _within_group_triple(solver)
    si, _base, _gains, windowed = _within_group_si(
        LC, S, B, sn_mesh, inner_schedule=solver.inner_schedule,
        max_iter=600, tol=1e-12,
    )
    if windowed:
        raise AssertionError("1-D slab must not window")
    q_ext = TimedFullField(
        bulk=AngularSourceSink.from_isotropic(
            np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0), sn_mesh,
        ),
        boundary=BoundarySourceSink.zeros_on(sn_mesh),
        _history=(), history_depth=2,
    )
    ig = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    psi, _ = si.solve(q_ext, initial_guess=ig)
    return solver, (LC - S - B), q_ext, psi


# ── T4.1 — from_balance mints (positive) / mis-typed operand raises (negative)


@pytest.mark.foundation
def test_from_balance_mints_residual_with_correct_type_units_space():
    r"""[L11 paired] well-formed same-class source operands → AngularResidual
    on the ``"angular_residual"`` space; a flux operand (wrong units) RAISES."""
    fuel = get_mixture("A", "2g")
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 9), mat_ids=np.zeros(8, dtype=int),
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    sn_mesh = SNMesh(mesh, quad, {0: fuel})
    rng = np.random.default_rng(208)
    shape = (quad.N, sn_mesh.ng, *sn_mesh.spatial_shape)
    a_psi = AngularSourceSink.from_mesh(rng.standard_normal(shape), sn_mesh)
    q = AngularSourceSink.from_mesh(rng.standard_normal(shape), sn_mesh)
    r = AngularResidual.from_balance(lhs=a_psi, rhs=q)  # POSITIVE
    if type(r) is not AngularResidual:
        raise AssertionError(f"from_balance returned {type(r).__name__}")
    np.testing.assert_array_equal(r.values, a_psi.values - q.values)
    if r.space.name != "angular_residual":
        raise AssertionError(f"residual on wrong space {r.space.name!r}")
    flux = AngularFlux.from_mesh(rng.standard_normal(shape), sn_mesh)
    with pytest.raises(TypeError):  # NEGATIVE — flux operand (wrong units/class)
        _ = AngularResidual.from_balance(lhs=flux, rhs=q)  # type: ignore[arg-type]


# ── T4.2 — balance_map ≈ 0 at convergence / detectably ≠ 0 on perturbed ψ


@pytest.mark.l0
def test_balance_map_zero_at_convergence_nonzero_on_perturbation():
    r"""[term: (L+C-S-B)ψ-q] balance_map() ≈ 0 per-cell at the SI fixed point,
    detectably ≠ 0 when one interior cell is perturbed.

    The typed per-ordinate flat-flux residual probe (Signature 1 / §H3): the
    per-cell defect global conservation HIDES. L11-paired: converged → ≈0
    (POSITIVE); perturbed → ≠0 localised (NEGATIVE — proves the map reads the
    defect, not a latent dud)."""
    _solver, loss_op, q_ext, psi = _converged_slab_2g()
    r = evaluate_residual(loss_op, psi, q_ext)            # POSITIVE
    bmap = r.bulk.balance_map()                            # (ng, *spatial)
    q_scale = float(np.abs(q_ext.bulk.values).max())
    rel = float(np.abs(bmap).max()) / max(q_scale, 1e-30)
    if not (rel < 1e-7):
        raise AssertionError(f"balance_map not ≈0 at convergence: rel={rel:.2e}")

    # NEGATIVE — perturb one interior cell (group 0, all ordinates) by 10 %.
    bad_vals = psi.bulk.values.copy()
    ix = psi.bulk.values.shape[2] // 2
    bad_vals[:, 0, ix] *= 1.1
    psi_bad = TimedFullField(
        bulk=AngularFlux.from_mesh(bad_vals, psi.bulk.mesh),
        boundary=psi.boundary, _history=(), history_depth=psi.history_depth,
    )
    r_bad = evaluate_residual(loss_op, psi_bad, q_ext)
    bmap_bad = r_bad.bulk.balance_map()
    if not (np.abs(bmap_bad).max() > 100.0 * max(np.abs(bmap).max(), 1e-30)):
        raise AssertionError(
            "balance_map did not detect a deliberate ψ perturbation — latent dud."
        )
    # Locality: the perturbed cell's balance is broken (its residual is large).
    # The defect ALSO spreads to streaming neighbours (the L stencil couples
    # cells), so the GLOBAL peak may sit at ix±1 — assert the perturbed cell
    # itself, not that it is the global argmax.
    if not (np.abs(bmap_bad[:, ix]).max()
            > 100.0 * max(np.abs(bmap).max(), 1e-30)):
        raise AssertionError(
            "balance_map shows no defect AT the perturbed cell — not localised"
        )


# ── T4.3 — boundary_vs_interior_split quadrature-sums to ‖r‖


@pytest.mark.foundation
def test_boundary_vs_interior_split_quadrature_sum():
    r"""sqrt(b² + i²) == ‖r‖ (the composite flat L2 — the typed residual is
    bulk ⊕ boundary; the split is the norms of the two members)."""
    _solver, loss_op, q_ext, psi = _converged_slab_2g()
    r = evaluate_residual(loss_op, psi, q_ext)
    b, i = boundary_vs_interior_split(r)
    total = float(np.linalg.norm(r.to_flat()))
    np.testing.assert_allclose(np.hypot(b, i), total, rtol=1e-12)


# ── T4.4 — relative_to(source) is the norm ratio


@pytest.mark.foundation
def test_relative_to_source_is_norm_ratio():
    r"""relative_to(q) == ‖r‖/‖q‖ — the tolerance-portable residual criterion."""
    _solver, loss_op, q_ext, psi = _converged_slab_2g()
    r = evaluate_residual(loss_op, psi, q_ext)
    np.testing.assert_allclose(
        r.bulk.relative_to(q_ext.bulk), r.bulk.l2 / q_ext.bulk.l2, rtol=1e-12,
    )
