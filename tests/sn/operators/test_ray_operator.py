r"""L1 convergence of ``A_BB`` — the radial-characteristic resolvent vs the closed-form ODE.

The operator-level L1 anchor for :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
(``A_BB``, campaign step 1c of the coupled-block-operator campaign,
``.claude/plans/archive/coupled_block_operator_campaign.md``). ``A_BB.solve`` is the
two-leg Carlson march (the resolvent action :math:`A_{BB}^{-1}`); its
``μ = -1`` inward leg solves the starting-direction ODE
:math:`-\,\mathrm d\phi/\mathrm dr + \sigma\phi = q` on ``[0, R]`` with the
prescribed ``r = R`` inflow datum :math:`\phi(R) = \phi_R`.

**Claim layer + pillar (vv-principles §1.5).** This is a *convergence-order*
claim (the lowest layer — pure math), verified against a **closed-form**
reference: the exact attenuation integral

.. math::

   \phi(r) \;=\; \frac{q}{\sigma} \;+\;
                 \Bigl(\phi_R - \frac{q}{\sigma}\Bigr)\,
                 e^{-\sigma\,(R - r)} .

**Structural independence (vv-principles §1).** The reference is an analytic
exponential — it does NO marching, so it is structurally independent of the
diamond-difference recurrence the operator wraps (Hébert 3.434/3.435). It is
NOT a hand-re-execution of that recurrence (the ERR-032 procedural-independence
hazard). The engine ``carlson_inward_sweep_from_source`` carries its OWN
closed-form convergence gate in
``tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py::TestB1DirectSolverClosedForm``;
THIS module is the operator-level anchor — it exercises the WRAP layer (the
``RadialCharacteristicField`` interior/boundary member views, the
source-block reading, the composite flux output) that the bare-engine gate
never touches.

**Config discipline (§0.6).** The fixture is genuinely ≥2G — DISTINCT
per-group :math:`(\sigma, q, \phi_R)`, both groups anchored — so a group-axis
transpose in the view layer is caught (the max-over-groups error is the gate
measurand). Both the MANDATORY graded-mesh leg (vv Mode 5 — a uniform mesh is
blind to a per-cell ``dr[k]`` index drift) AND the uniform leg run. The
``r = R`` inflow is a genuinely-exponential :math:`\phi_R \neq q/\sigma`
(the flat leg ``φ_R = q/σ`` — pinned exactly by the equilibrium gate in
``test_psi_half_coupling.py`` — nulls the dynamics).

**Runtime discipline (vv Mode 8).** Every gate reddens via
:func:`pytest.fail` (a function call — fires under the canonical
``python -O``), never a bare ``assert``.
"""
from __future__ import annotations

import numpy as np
import pytest

import orpheus.sn.operators.radial_characteristic as _rc_mod
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.radial_characteristic import RadialCharacteristicOperator
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)

pytestmark = pytest.mark.l1

# ── The manufactured 2G convergence problem (constant σ, constant q per group) ──
#
# DISTINCT per-group data so the max-over-groups error catches a group-axis
# transpose in the WRAP.  R = 1.0 keeps σ·R modest (≤ 1.3) so the asymptotic
# O(Δr²) regime is reached by nx = 64; φ_R ≠ q/σ makes the exponential genuinely
# non-flat (the flat leg nulls the dynamics — §0.6).
_R = 1.0
_SIGMA = np.array([1.3, 0.7])       # (ng,) constant-in-r total XS per group
_Q = np.array([0.7, 1.1])           # (ng,) constant q̄ per group
_PHI_R = np.array([2.0, 0.5])       # (ng,) r = R inflow datum per group
_NG = 2


def _mixture(ng: int):
    """A group-graded diagonal-scatter mixture (the operator reads its OWN σ_t,
    so the mesh mixture only has to build a seed-carrying sphere)."""
    st = np.array([1.0 + 0.4 * g for g in range(ng)])
    ss = np.diag([0.4 * (1.0 + 0.4 * g) for g in range(ng)])
    return make_mixture(sig_t=st, sig_c=st - ss.sum(axis=0), sig_f=np.zeros(ng),
                        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=ss)


def _sphere(edges: np.ndarray, ng: int = _NG):
    """A seed-carrying sphere (GL S4) on the given radial edges."""
    nx = edges.size - 1
    mesh = Mesh1D(edges=edges, mat_ids=np.zeros(nx, dtype=int),
                  coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"))
    return SNMesh(mesh, Quadrature.gauss_legendre(4), {0: _mixture(ng)})


def _uniform_edges(nx: int) -> np.ndarray:
    return np.linspace(0.0, _R, nx + 1)


def _graded_edges(nx: int, p: float = 1.5) -> np.ndarray:
    r"""Self-similar power grading ``r_j = R·(j/nx)^p`` — genuinely non-uniform
    ``dr`` (the per-cell ``dr[k]``-indexing leg; a uniform mesh is BLIND to a
    ``dr[k] → dr[k±1]`` index drift, vv Mode 5)."""
    j = np.arange(nx + 1, dtype=float)
    return _R * (j / nx) ** p


def _radial_characteristic_exact(r: np.ndarray, q: float, sigma: float,
                                 phi_R: float) -> np.ndarray:
    r"""φ(r) = q/σ + (φ_R − q/σ)·exp(−σ·(R − r)) — the closed attenuation
    integral of −dφ/dr + σφ = q with φ(R) = φ_R (constant σ, q).

    Structurally independent of the DD recurrence (no marching — pure exp)."""
    return q / sigma + (phi_R - q / sigma) * np.exp(-sigma * (_R - r))


def _solve_inward_cells(edges: np.ndarray) -> np.ndarray:
    """Build the operator + the uniform-q½ source, return ``solve``'s inward-leg
    cells ``(ng, nx)`` = the ψ½ marched on the ``μ = -1`` ray."""
    sn = _sphere(edges)
    nx = sn.nx
    # σ_t constant in r, distinct per group — the C_ray collision coefficient,
    # a typed mesh-bound CrossSectionField (the operator's mesh-identity guard).
    sigma_t = np.stack([np.full(nx, _SIGMA[g]) for g in range(_NG)], axis=0)
    op = RadialCharacteristicOperator(sn, CrossSectionField.from_mesh(sigma_t, sn))
    # System B's SOURCE composite (4e native — the block boundary speaks the
    # split member directly; no unified bridge).
    src = RadialCharacteristicField.source_zeros(sn.radial_characteristic_field_space)
    for g in range(_NG):
        # q½ built DIRECTLY (uniform q̄ per group) — NOT via the R14 fold (step 2).
        src.interior.cells(0, -1)[g, :] = _Q[g]
        src.boundary.corner(0, -1)[g] = _PHI_R[g]   # r = R Dirichlet datum
    return op.solve(src).interior.cells(0, -1)


def _max_group_inf_error(edges: np.ndarray) -> float:
    """‖solve inward cells − φ_exact(r_center)‖_∞, maxed over both groups."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    cells = _solve_inward_cells(edges)
    return max(
        float(np.max(np.abs(
            cells[g] - _radial_characteristic_exact(centers, _Q[g], _SIGMA[g],
                                                     _PHI_R[g]))))
        for g in range(_NG)
    )


# Measured (fresh operator, distinct-group): uniform err(64) ≈ 7.4e-5 with
# ratios 3.75→3.93 (order 1.91→1.98); graded err(64) ≈ 1.6e-4 with ratios
# 3.53→3.87 (order 1.82→1.95).  The band brackets O(Δr²) (ratio → 4) with
# cross-machine headroom; a wrong-limit mutant plateaus at O(1e-1) (≫ 5e-4).
_RATIO_LO, _RATIO_HI = 3.3, 4.7
_FINEST_TOL = 5.0e-4
_REFINEMENTS = (8, 16, 32, 64)


def _order_ok(errs: list[float]) -> bool:
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    return all(_RATIO_LO <= rat <= _RATIO_HI for rat in ratios)


def _check_convergence(edges_of, *, expect_pass: bool, label: str) -> None:
    r"""The O(Δr²) criterion: ``err_∞`` shrinks by ≈4× per halving AND the
    finest error is genuinely small.  ``expect_pass=False`` inverts it (the
    mutation-teeth leg)."""
    errs = [_max_group_inf_error(edges_of(nx)) for nx in _REFINEMENTS]
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    passed = _order_ok(errs) and errs[-1] < _FINEST_TOL
    if expect_pass and not passed:
        pytest.fail(
            f"A_BB.solve[{label}] lost O(Δr²) convergence to the closed-form "
            f"ODE φ=q/σ+(φ_R−q/σ)e^{{−σ(R−r)}} — errs={errs!r}, ratios={ratios!r} "
            f"(band [{_RATIO_LO}, {_RATIO_HI}], finest < {_FINEST_TOL})."
        )
    if not expect_pass and passed:
        pytest.fail(
            f"A_BB.solve[{label}] STILL converges after the mutation — the "
            f"convergence gate has no teeth (errs={errs!r}, ratios={ratios!r})."
        )


def _dd_mutant(kind: str):
    r"""A mutated twin of ``carlson_inward_sweep_from_source`` (the DD engine
    A_BB wraps) — same ``(phi_cells, phi_face_final)`` signature so the operator
    consumes it unchanged.  Monkeypatched into the operator's module namespace
    (never a ``git checkout``).

    ``kind``:
    - ``closure_sign``   — Hébert (3.435) ``2φ − f`` → ``2φ + f`` (Mode 1).
    - ``denom_sign``     — DD denominator ``Δr·σ + 2`` → ``Δr·σ − 2`` (Mode 3).
    - ``diamond_factor`` — the ``2·f_in`` numerator weight → ``1·f_in`` (Mode 3).
    - ``index_drift``    — ``dr[k]`` → ``dr[k−1]`` (Mode 5 — uniform-blind).
    """

    def solver(Q_bar, sigma_t, dr, bc_outer_value):
        ng, nx = Q_bar.shape
        phi_aux = np.zeros((ng, nx), dtype=Q_bar.dtype)
        phi_face = bc_outer_value.copy()
        for k in range(nx - 1, -1, -1):
            dk = dr[k - 1] if kind == "index_drift" else dr[k]
            denom = (dk * sigma_t[:, k] - 2.0) if kind == "denom_sign" \
                else (dk * sigma_t[:, k] + 2.0)
            two = 1.0 if kind == "diamond_factor" else 2.0
            phi_cell = (dk * Q_bar[:, k] + two * phi_face) / denom
            phi_aux[:, k] = phi_cell
            phi_face = (2.0 * phi_cell + phi_face) if kind == "closure_sign" \
                else (2.0 * phi_cell - phi_face)
        return phi_aux, phi_face

    return solver


class TestA_BB_RadialBVP:
    r"""L1 — ``A_BB.solve``'s inward leg converges O(Δr²) to the exponential ODE.

    The convergence-order claim on the resolvent's ``μ = -1`` march, anchored on
    the closed-form attenuation integral (structurally independent of the DD
    recurrence). Necessary-but-not-sufficient on its own — paired with the
    foundation invariants (adjoint consistency, WRAP bit-identity, pole
    continuation, Dirichlet propagation, fixed-source equilibrium) in
    ``test_psi_half_coupling.py::TestA_BB_RadialBVP`` that pin the CONVERGED
    value and the two-leg structure.
    """

    @pytest.mark.verifies("hebert-3-434")
    @pytest.mark.verifies("hebert-3-435")
    def test_solve_converges_to_closed_form_uniform(self) -> None:
        r"""Uniform mesh: err_∞ shrinks at O(Δr²) toward the exact exponential."""
        _check_convergence(_uniform_edges, expect_pass=True, label="uniform")

    @pytest.mark.verifies("hebert-3-434")
    @pytest.mark.verifies("hebert-3-435")
    def test_solve_converges_to_closed_form_graded(self) -> None:
        r"""The MANDATORY per-cell ``dr[k]``-indexing leg (vv Mode 5) — a graded
        mesh drives the per-cell width index out of the uniform blind spot."""
        _check_convergence(_graded_edges, expect_pass=True, label="graded")

    # ── Mutation teeth (monkeypatch the engine — never git-revert) ─────────

    @pytest.mark.parametrize("kind", ["closure_sign", "denom_sign",
                                      "diamond_factor"])
    def test_teeth_dd_coefficient_mutations_red_both_meshes(
        self, kind, monkeypatch,
    ) -> None:
        r"""A wrong DD coefficient (Hébert 3.434/3.435) → convergence to the
        WRONG limit — the operator gate REDs on BOTH uniform and graded (the
        error is O(1) in the coefficient, so no mesh hides it)."""
        monkeypatch.setattr(
            _rc_mod, "carlson_inward_sweep_from_source", _dd_mutant(kind))
        _check_convergence(_uniform_edges, expect_pass=False,
                           label=f"{kind}/uniform")
        _check_convergence(_graded_edges, expect_pass=False,
                           label=f"{kind}/graded")

    def test_teeth_index_drift_uniform_blind_graded_red(self, monkeypatch) -> None:
        r"""``dr[k] → dr[k−1]`` (Mode 5): the config-blindness keystone — a
        uniform mesh CANNOT see the drift (all widths equal, arithmetically
        identical), the graded mesh REDs it. This asymmetry IS the evidence that
        the graded leg is load-bearing (§0.6)."""
        monkeypatch.setattr(
            _rc_mod, "carlson_inward_sweep_from_source", _dd_mutant("index_drift"))
        _check_convergence(_uniform_edges, expect_pass=True,
                           label="index_drift/uniform (blind — expected)")
        _check_convergence(_graded_edges, expect_pass=False,
                           label="index_drift/graded")
