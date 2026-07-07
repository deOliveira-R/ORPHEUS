r"""L0/L1 tests for the Hébert §3.9.4 starting-direction DD march.

Issue #168 Phase D shipped the curvilinear starting-direction seed as a
swappable strategy family composed onto
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`.
The full retirement narrative lives in the production module docstring of
:mod:`orpheus.sn.spatial.psi_half_angle_seed` (its "Development history").

#282 route (a) (#280 Phase 2.5d, 2026-07-04) **retired the whole strategy
zoo**: the starting-direction flux is now first-class STATE (the
``RadialCharacteristicField`` block of the augmented composite), the SOLVE
marches it directly from the TRUE q½ source, and the APPLY reads the given
carrier.  What SURVIVES as the direct engine is the free function
:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
— the Hébert (3.434)-(3.435) march, now returning the
``(phi_cells, phi_face_final)`` tuple.

Test layers
-----------

* **L0** (term verification via flat-ψ algebraic identity + multi-region
  + multi-group probes): the canonical Hébert §3.9.4 (3.432)-(3.435)
  algebra is verified directly against the source-driven march.
* **L1** (§16.B direct-solver gates, 2.5d d2): ``carlson_inward_sweep_from_source``
  vs the structurally-independent closed form of the μ = −1
  starting-direction ODE, plus the R14 Legendre fold helper.

The former foundation block (Protocol conformance / registry / immutability
/ zero-seed bit-identity / structural-independence probe / default-seed
field) pinned the retired strategy machinery and was dropped with the zoo
(API-smoke + retired-strategy-only tests — no successor).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.spaces.radial_characteristic_space import (
    fold_moments_to_radial_characteristic,
)
from orpheus.sn.spatial.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
)


# ═══════════════════════════════════════════════════════════════════════
# L0 term verification — flat-ψ algebraic identity
# ═══════════════════════════════════════════════════════════════════════


class TestCarlsonFlatPsiAlgebraicIdentity:
    r"""Hébert §3.9.4 flat-ψ trace (literature memo §3 algebraic proof).

    With flat ψ = C, consistent source ``Q̄ = Σ_t · C``, and
    consistent BC ``bc_outer = C``, the Hébert (3.434)-(3.435) inward
    recurrence is self-similar: every cell and face stays at C.

    The L0 verification:

    .. math::

       \\bar\\phi_i = \\frac{\\Delta r \\cdot \\Sigma \\cdot C + 2 C}
                              {\\Delta r \\cdot \\Sigma + 2}
                    = C \\cdot \\frac{\\Delta r \\Sigma + 2}{\\Delta r \\Sigma + 2}
                    = C

    #282 route (a): the ``Q̄`` construction moved OUT of the retired
    strategy (which folded ``Q̄ = ½·Σ_t·φ_0`` from a per-ordinate ψ field
    at weights summing to 2, i.e. ``Q̄ = Σ_t·C`` for flat ψ = C) — the
    tests now build the flat/isotropic source explicitly and call the
    surviving :func:`carlson_inward_sweep_from_source` directly.
    """

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_flat_psi_identity_reflective(self) -> None:
        r"""Flat Q̄ = Σ_t·C, reflective bc_outer = C → all cells = C."""
        ng, nx = 1, 10
        C = 1.0
        sigma_t_value = 0.5
        sigma_t = np.full((ng, nx), sigma_t_value)
        dr = np.full(nx, 2.0 / nx)
        # Flat consistent source Q̄ = Σ_t · C (the L=0 fold of flat ψ = C
        # the retired strategy produced internally).  bc_outer = C
        # (reflective: outer face at μ = −1 mirrors μ = +1 = C).
        Q_bar = np.full((ng, nx), sigma_t_value * C)
        bc_outer = np.full((ng,), C)
        cells, _face = carlson_inward_sweep_from_source(
            Q_bar, sigma_t, dr, bc_outer,
        )
        # Every cell should be C = 1 to machine precision.
        np.testing.assert_allclose(cells, C * np.ones((ng, nx)), atol=1e-15)

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_flat_psi_identity_varying_C(self) -> None:
        r"""Identity for C ≠ 1 — algebraic invariance under scaling."""
        for C in (0.7, 2.3, 100.0):
            ng, nx = 1, 8
            sigma_t_value = 0.5
            sigma_t = np.full((ng, nx), sigma_t_value)
            dr = np.full(nx, 2.0 / nx)
            Q_bar = np.full((ng, nx), sigma_t_value * C)
            bc_outer = np.full((ng,), C)
            cells, _face = carlson_inward_sweep_from_source(
                Q_bar, sigma_t, dr, bc_outer,
            )
            np.testing.assert_allclose(
                cells, C * np.ones((ng, nx)), atol=1e-13 * abs(C),
                err_msg=f"flat-source identity failed at C = {C}",
            )

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_vacuum_BC_flat_source_nx_3(self) -> None:
        r"""Vacuum BC + flat source Q̄ = 0.5, nx = 3: hand calculation.

        Per literature memo §3.b (vacuum trace) the recurrence marches
        inward from the vacuum outer face (φ̄_{nx+1/2} = 0).  The Q̄
        construction (½·Σ_t·φ_0 = ½·0.5·2 = 0.5 for flat ψ = 1 at
        weights summing to 2) is now explicit; the hand calc below is
        unchanged (the recurrence math is identical — only the Q̄
        assembly moved external).
        """
        nx = 3
        ng = 1
        dr_val = 2.0 / nx  # 0.667
        sigma_t_value = 0.5
        sigma_t = np.full((ng, nx), sigma_t_value)
        dr = np.full(nx, dr_val)
        # Q̄ = ½·Σ_t·φ_0 = ½·0.5·2.0 = 0.5 (flat isotropic source).
        Q_bar = np.full((ng, nx), 0.5)
        bc_outer = np.array([0.0])  # vacuum
        # Hand calc:
        # phi_face = 0 (outer)
        # k=2: denom = 0.667·0.5 + 2 = 2.333
        #      phi_cell = (0.667 · 0.5 + 0) / 2.333 = 0.333/2.333 = 0.1429
        #      phi_face = 2·0.1429 - 0 = 0.2857
        # k=1: phi_cell = (0.667·0.5 + 2·0.2857) / 2.333 = 0.388
        #      phi_face = 2·0.388 - 0.2857 = 0.490
        # k=0: phi_cell = (0.333 + 2·0.490) / 2.333 = 0.563
        expected_2 = (dr_val * 0.5) / (dr_val * 0.5 + 2)
        face_at_2 = 2 * expected_2
        expected_1 = (dr_val * 0.5 + 2 * face_at_2) / (dr_val * 0.5 + 2)
        face_at_1 = 2 * expected_1 - face_at_2
        expected_0 = (dr_val * 0.5 + 2 * face_at_1) / (dr_val * 0.5 + 2)
        expected = np.array([[expected_0, expected_1, expected_2]])

        cells, _face = carlson_inward_sweep_from_source(
            Q_bar, sigma_t, dr, bc_outer,
        )
        np.testing.assert_allclose(cells, expected, rtol=1e-13)

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_multi_region_sigma_t_step(self) -> None:
        r"""Σ_t jumps at midpoint (heterogeneous medium).

        Per Hébert §3.9.4 the recurrence handles per-cell Σ_t
        trivially (no special treatment for material discontinuities
        as long as the radial mesh respects boundaries).
        """
        ng, nx = 1, 6
        # σ_t step: first 3 cells at 0.3, last 3 at 1.0.
        sigma_t = np.zeros((ng, nx))
        sigma_t[:, :3] = 0.3
        sigma_t[:, 3:] = 1.0
        dr = np.full(nx, 1.0 / nx)
        # Flat isotropic source per cell: Q̄[k] = ½·Σ_t[k]·φ_0 = ½·Σ_t[k]·2
        #   = Σ_t[k]  (flat ψ = 1, weights summing to 2).
        Q_bar = sigma_t.copy()
        bc_outer = np.array([0.0])  # vacuum
        cells, _face = carlson_inward_sweep_from_source(
            Q_bar, sigma_t, dr, bc_outer,
        )
        # Hand-trace: phi_face_outer = 0; Q̄[k] = σ_t[k].
        phi_face = 0.0
        expected_cells = np.zeros(nx)
        for k in range(nx - 1, -1, -1):
            sig = sigma_t[0, k]
            Q_bar_k = sig  # = 0.5 · sig · 2.0
            phi_cell = (dr[k] * Q_bar_k + 2.0 * phi_face) / (dr[k] * sig + 2.0)
            expected_cells[k] = phi_cell
            phi_face = 2.0 * phi_cell - phi_face
        np.testing.assert_allclose(cells[0], expected_cells, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# L0 term verification — linearity (foundation of operator algebra)
# ═══════════════════════════════════════════════════════════════════════


class TestSeedLinearity:
    """The source-driven march is linear in ``Q_bar``.

    Required by the operator algebra: the M-M closure is linear, and
    a non-linear seed would propagate non-linearity into the apply
    matvec.  Linearity is a pre-condition for
    :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply_transpose`
    correctness (dense-matrix probing assumes linearity).
    """

    @pytest.mark.foundation
    def test_carlson_inward_sweep_is_linear(self) -> None:
        ng, nx = 1, 5
        sigma_t = np.full((ng, nx), 0.5)
        dr = np.full(nx, 2.0 / nx)
        # bc_outer = 0 so the recurrence is purely linear in Q_bar (a
        # non-zero bc_outer adds an affine offset — its linearity is
        # pinned separately in test_sweep_vs_apply_consistency.py).
        bc_outer = np.zeros(ng)
        rng = np.random.default_rng(seed=42)
        Q1 = rng.standard_normal((ng, nx))
        Q2 = rng.standard_normal((ng, nx))
        alpha = 1.7
        beta = -0.3
        cells_combined, _ = carlson_inward_sweep_from_source(
            alpha * Q1 + beta * Q2, sigma_t, dr, bc_outer,
        )
        cells_a, _ = carlson_inward_sweep_from_source(Q1, sigma_t, dr, bc_outer)
        cells_b, _ = carlson_inward_sweep_from_source(Q2, sigma_t, dr, bc_outer)
        np.testing.assert_allclose(
            cells_combined, alpha * cells_a + beta * cells_b, atol=1e-13,
        )


# ═══════════════════════════════════════════════════════════════════════
# §16.B — the DIRECT-SEED value gates (2.5d d2, #282 route (a))
# ═══════════════════════════════════════════════════════════════════════
#
# B1 pins carlson_inward_sweep_from_source — the SOURCE-DRIVEN direct
# entry (route (a)'s solve engine) — against the structurally-
# independent closed form of the μ = −1 starting-direction ODE.  Every
# pre-existing pin above either runs the FLAT identity (exact but
# nulls the dynamics) or re-executes the same recurrence by hand
# (procedurally, not structurally, independent — the ERR-032 hazard);
# B1b is the missing structural pillar.  B2 pins the R14 fold helper
# (the q½ source construction) on the Legendre identity P_ℓ(−1) =
# (−1)^ℓ — the anisotropic case is MANUFACTURED here because an
# all-isotropic suite is silently blind to a mis-signed φ_{ℓ≥1} fold
# (§0.6 isotropic-snapshot blindness).


def _radial_characteristic_exact(r, q, sigma, phi_R, R):
    r"""φ(r) = q/σ + (φ_R − q/σ)·exp(−σ·(R − r)) — the closed attenuation
    integral of −dφ/dr + σφ = q with φ(R) = φ_R (constant σ, q).
    Structurally independent of the DD recurrence (no marching)."""
    return q / sigma + (phi_R - q / sigma) * np.exp(-sigma * (R - r))


def _uniform_edges(nx: int, R: float = 1.0) -> np.ndarray:
    return np.linspace(0.0, R, nx + 1)


def _graded_edges(nx: int, R: float = 1.0, p: float = 1.5) -> np.ndarray:
    r"""Self-similar power grading r_j = R·(j/nx)^p — genuinely nonuniform
    dr (the per-cell dr[k]-indexing leg; a uniform mesh is BLIND to a
    dr[k] → dr[k±1] index drift, vv Mode 5)."""
    j = np.arange(nx + 1, dtype=float)
    return R * (j / nx) ** p


def _direct_solver_inf_error(solver, edges, *, q=0.7, sigma=1.3, phi_R=2.0):
    r"""‖solver output − φ_exact(r_center)‖_∞ on one mesh (ng = 1).

    ``solver`` returns the ``(phi_cells, phi_face_final)`` tuple (both the
    real :func:`carlson_inward_sweep_from_source` and the mutation twins
    honour it), so the cells are unpacked uniformly."""
    R = float(edges[-1])
    dr = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    nx = dr.size
    Q_bar = np.full((1, nx), q)
    sig = np.full((1, nx), sigma)
    bc = np.array([phi_R])
    cells, _face = solver(Q_bar, sig, dr, bc)
    exact = _radial_characteristic_exact(centers, q, sigma, phi_R, R)
    return float(np.max(np.abs(cells[0] - exact)))


def _b1_order_check(solver, edges_of, *, expect_pass: bool, label: str):
    r"""The B1b criterion: err_∞ shrinks at O(Δr²) across nx ∈ {8..64}
    (ratio ∈ [3.4, 4.6]) AND the finest error is genuinely small.
    ``expect_pass=False`` inverts it (the mutation-teeth leg)."""
    errs = [
        _direct_solver_inf_error(solver, edges_of(nx))
        for nx in (8, 16, 32, 64)
    ]
    ratios = [errs[i] / errs[i + 1] for i in range(3)]
    order_ok = all(3.4 <= rat <= 4.6 for rat in ratios)
    # Measured true-solver scale: err(64) ≈ 7.4e-5 (uniform); a wrong-limit
    # mutant plateaus at O(1e-1) — 2e-4 sits 3 orders below the plateau.
    small_ok = errs[-1] < 2.0e-4
    if expect_pass and not (order_ok and small_ok):
        pytest.fail(
            f"B1b[{label}]: direct solver lost O(Δr²) convergence to the "
            f"closed form — errs={errs!r}, ratios={ratios!r}"
        )
    if not expect_pass and (order_ok and small_ok):
        pytest.fail(
            f"B1b-teeth[{label}]: the mutated solver still passes the "
            f"convergence criterion — the gate has no teeth "
            f"(errs={errs!r}, ratios={ratios!r})"
        )


class TestB1DirectSolverClosedForm:
    r"""§16.B B1 — carlson_inward_sweep_from_source vs the exponential ODE.

    B1a (the exact flat leg ``φ_R = q/σ ⟹ output ≡ q/σ``) is already
    pinned above by ``TestCarlsonFlatPsiAlgebraicIdentity`` /
    ``test_carlson_vacuum_BC_flat_source_nx_3`` — cited, not
    duplicated.  B1b is the genuinely-exponential closed-form leg
    (φ_R ≠ q/σ) the flat pins null.
    """

    @pytest.mark.l1
    @pytest.mark.verifies("hebert-3-434")
    @pytest.mark.verifies("hebert-3-435")
    def test_b1b_uniform_mesh_second_order(self) -> None:
        _b1_order_check(
            carlson_inward_sweep_from_source, _uniform_edges,
            expect_pass=True, label="uniform",
        )

    @pytest.mark.l1
    @pytest.mark.verifies("hebert-3-434")
    @pytest.mark.verifies("hebert-3-435")
    def test_b1b_graded_mesh_second_order(self) -> None:
        r"""The MANDATORY per-cell dr[k]-indexing leg (vv Mode 5)."""
        _b1_order_check(
            carlson_inward_sweep_from_source, _graded_edges,
            expect_pass=True, label="graded",
        )

    # ── Mutation teeth (in-process twins — never git-revert) ─────────

    @staticmethod
    def _mutant(kind: str):
        r"""A mutated twin of the direct solver (the §16.B B1 matrix).

        Returns the SAME ``(phi_cells, phi_face_final)`` tuple shape as
        :func:`carlson_inward_sweep_from_source` so
        :func:`_direct_solver_inf_error` unpacks it uniformly."""

        def solver(Q_bar, sigma_t, dr, bc_outer_value):
            ng, nx = Q_bar.shape
            phi_aux = np.zeros((ng, nx), dtype=Q_bar.dtype)
            phi_face = bc_outer_value.copy()
            for k in range(nx - 1, -1, -1):
                dk = dr[k - 1] if kind == "index_drift" else dr[k]
                if kind == "denom_sign":
                    denom = dk * sigma_t[:, k] - 2.0
                else:
                    denom = dk * sigma_t[:, k] + 2.0
                two = 1.0 if kind == "diamond_factor" else 2.0
                phi_cell = (dk * Q_bar[:, k] + two * phi_face) / denom
                phi_aux[:, k] = phi_cell
                if kind == "closure_sign":
                    phi_face = 2.0 * phi_cell + phi_face  # Hébert (3.435) sign
                else:
                    phi_face = 2.0 * phi_cell - phi_face
            return phi_aux, phi_face

        return solver

    @pytest.mark.l1
    @pytest.mark.parametrize(
        "kind", ["closure_sign", "denom_sign", "diamond_factor"],
    )
    def test_teeth_wrong_limit_mutations_red_both_rows(self, kind) -> None:
        r"""Mutations 1–3: convergence to the WRONG limit — both rows RED."""
        _b1_order_check(
            self._mutant(kind), _uniform_edges,
            expect_pass=False, label=f"{kind}/uniform",
        )
        _b1_order_check(
            self._mutant(kind), _graded_edges,
            expect_pass=False, label=f"{kind}/graded",
        )

    @pytest.mark.l1
    def test_teeth_index_drift_uniform_green_graded_red(self) -> None:
        r"""Mutation 4 (dr[k] → dr[k−1]): the config-blindness keystone —
        a uniform mesh CANNOT see the drift (all widths equal, the mutant
        is arithmetically identical), the graded row REDS it."""
        _b1_order_check(
            self._mutant("index_drift"), _uniform_edges,
            expect_pass=True, label="index_drift/uniform (blind — expected)",
        )
        _b1_order_check(
            self._mutant("index_drift"), _graded_edges,
            expect_pass=False, label="index_drift/graded",
        )


class TestB2RadialCharacteristicSourceFold:
    r"""§16.B B2 — the R14 fold helper ``Q̄(±1) = Σ_ℓ (2ℓ+1)/2·Q_ℓ·(±1)^ℓ``.

    The single source of the q½ construction (the S/F seed arms and the
    d3 q_ext factories all fold through it).  B2b manufactures the
    anisotropic case the isotropic suite is blind to (§0.6); production
    reach today is ℓ = 0 (the arms feed the iso moment only — the
    honest-scope note on the helper), so B2b is the live foundation pin
    on the helper itself.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("hebert-3-432-source")
    def test_b2a_isotropic_fold_is_half_q0(self) -> None:
        r"""ℓ=0: Q̄(±1) = ½·Q₀ — matching the direct solver's documented
        ``Q_bar`` convention ((2·0+1)/2 with P₀(±1) = 1)."""
        rng = np.random.default_rng(20260704)
        Q0 = rng.normal(size=(2, 5))
        np.testing.assert_array_equal(
            fold_moments_to_radial_characteristic(Q0[None], -1), 0.5 * Q0,
        )
        np.testing.assert_array_equal(
            fold_moments_to_radial_characteristic(Q0[None], +1), 0.5 * Q0,
        )

    @pytest.mark.foundation
    @pytest.mark.verifies("hebert-3-432-source")
    def test_b2b_two_term_legendre_hand_check(self) -> None:
        r"""ℓ ∈ {0, 1}: Q̄(−1) = ½Q₀ − (3/2)Q₁ and Q̄(+1) = ½Q₀ + (3/2)Q₁
        (P₁(−1) = −1 — the Mode-1/6 sign this gate exists to pin)."""
        rng = np.random.default_rng(7)
        Q = rng.normal(size=(2, 3, 4))  # (ℓ=0..1, ng, nx)
        np.testing.assert_allclose(
            fold_moments_to_radial_characteristic(Q, -1),
            0.5 * Q[0] - 1.5 * Q[1],
            rtol=1e-14,
        )
        np.testing.assert_allclose(
            fold_moments_to_radial_characteristic(Q, +1),
            0.5 * Q[0] + 1.5 * Q[1],
            rtol=1e-14,
        )

    @pytest.mark.foundation
    def test_b2b_teeth_dropped_sign_detected(self) -> None:
        r"""MUTATION: a fold with P_ℓ(−1) = +1 (the dropped (−1)^ℓ) gives
        ½Q₀ + (3/2)Q₁ at μ = −1 — assert the hand identity REDS on it."""
        rng = np.random.default_rng(8)
        Q = rng.normal(size=(2, 3, 4))

        def mutant_fold(moments, sign):
            ell = np.arange(moments.shape[0])
            coeff = (2.0 * ell + 1.0) / 2.0  # P_ℓ(sign) → +1: the Mode-1 drop
            return np.tensordot(coeff, moments, axes=(0, 0))

        wrong = mutant_fold(Q, -1)
        expected = 0.5 * Q[0] - 1.5 * Q[1]
        if np.allclose(wrong, expected, rtol=1e-14):
            pytest.fail("the sign-drop mutant matched the hand identity — "
                        "the B2b pin has no teeth")

    @pytest.mark.foundation
    def test_fold_input_validation(self) -> None:
        with pytest.raises(ValueError, match="sign must be"):
            fold_moments_to_radial_characteristic(np.ones((1, 2, 2)), 0)
        with pytest.raises(ValueError, match="leading ℓ axis"):
            fold_moments_to_radial_characteristic(np.float64(1.0), -1)
