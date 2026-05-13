r"""Phase G Step 1 — verification gates for :class:`AngularRedistribution`.

Step 1 promotes ``MorelMontryAngularSweep`` →
``AngularRedistribution(LinearOperator)``. The promotion wraps:

- The M-M angular recurrence ``MorelMontryAngularSweep.__call__``
  (currently at ``orpheus/sn/spatial/pole_angular_closure.py``
  lines 466+), and
- The Carlson seed (``carlson_inward_sweep_from_source`` /
  ``CarlsonInwardSweep``) composition.

Capability set declared: ``frozenset({CAP_APPLY})``.  ``CAP_SOLVE``
is **dropped at Step 1** per Pattern 4 (advertise only what works);
back-substitution is deferred to a later step.

The ``apply(psi_cells, *, ...) → R`` method delegates to
:meth:`MorelMontryAngularSweep.__call__`.  Input ``psi_cells`` is
``(ng, M, nx)`` — the natural M-M recurrence shape.

These tests extend (not replace) the existing 57 function-level
Carlson + twin-path tests at
``tests/sn/spatial/test_sweep_vs_apply_consistency.py``. The
existing tests pin the FUNCTION contract; this module pins the
OPERATOR contract on the same algebraic identities.

Gate design lives in
``.claude/agent-memory/test-architect/issue_196_phase_g_step1_verification_gates.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    cylindrical_streaming,
    spherical_streaming,
)
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    MissingCapability,
)
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.spatial import AngularRedistribution, MorelMontryAngularSweep
from orpheus.sn.spatial.psi_half_angle_seed import (
    CarlsonInwardSweep,
    CarlsonSweepContext,
    carlson_inward_sweep_from_source,
)


# ═══════════════════════════════════════════════════════════════════════
# Shared fixtures
# ═══════════════════════════════════════════════════════════════════════


def _spherical_redist_inputs(nx=8, radius=1.0, quad_order=4):
    """Build the M-M recurrence inputs for a 1-D sphere.

    Returns a dict with ``psi_cells_shape``, ``alpha_half``,
    ``redist_dAw``, ``tau_mm``, ``volume``, ``level_indices=None``,
    plus the underlying ``sn_op`` and ``quad`` for context building.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(quad_order)
    op = spherical_streaming(mesh, quad)
    # Spherical M-M inputs: single level, flat arrays (no per-level wrapping).
    # Volume per cell from streaming_terms (direction-independent geometric scalar).
    volume = np.array([op.streaming_terms(i, 0).volume for i in range(nx)])
    return {
        "quad": quad,
        "op": op,
        "mesh": mesh,
        "alpha_half": op.alpha_half,        # (M+1,)
        "redist_dAw": op.redist_dAw,        # (nx, M)
        "tau_mm": op.tau_mm,                # (M,)
        "volume": volume,                   # (nx,)
        "dr": mesh.widths,
        "level_indices": None,
    }


def _cylindrical_redist_inputs(nx=8, radius=1.0):
    """Build the M-M recurrence inputs for a 1-D cylinder."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = ProductQuadrature.create(n_mu=4, n_phi=4)
    op = cylindrical_streaming(mesh, quad)
    # Volume per cell from streaming_terms (direction-independent geometric scalar).
    volume = np.array([
        op.streaming_terms(i, 0, mu_level_idx=0).volume for i in range(nx)
    ])
    return {
        "quad": quad,
        "op": op,
        "mesh": mesh,
        "alpha_half": op.alpha_per_level,  # list of per-level
        "redist_dAw": op.redist_dAw_per_level,  # list
        "tau_mm": op.tau_mm_per_level,  # list
        "volume": volume,
        "dr": mesh.widths,
        "level_indices": quad.level_indices,
    }


def _build_inputs(geometry, nx=8):
    if geometry == "sphere":
        return _spherical_redist_inputs(nx=nx)
    elif geometry == "cylinder":
        return _cylindrical_redist_inputs(nx=nx)
    raise ValueError(geometry)


# ═══════════════════════════════════════════════════════════════════════
# Gate 3 — Capability surface for AngularRedistribution
# ═══════════════════════════════════════════════════════════════════════


class TestAngularRedistributionCapabilities:
    """Gate 3 — capability surface for ``AngularRedistribution``.

    Step 1 declares ``frozenset({CAP_APPLY})``. ``CAP_SOLVE`` is
    **dropped at Step 1** per Pattern 4 (advertise only what works):
    back-substitution is not implemented in this Step.

    See open-questions §4 in the gate-design memo; the
    method-implementer's Step 1 decision is to scope ``CAP_SOLVE``
    out of advertising.
    """

    @pytest.mark.foundation
    def test_cap_apply_declared(self):
        """``CAP_APPLY`` MUST advertise the M-M recurrence transform."""
        op = AngularRedistribution()
        assert CAP_APPLY in op.capabilities

    @pytest.mark.foundation
    def test_cap_solve_NOT_declared_at_step_1(self):
        """``CAP_SOLVE`` MUST NOT be declared at Step 1.

        Per Pattern 4 — advertise only what works.  The M-M recurrence
        is upper-triangular in ordinate index so its inverse exists
        algebraically, but no working back-substitution code ships at
        Step 1.  Advertising the capability without the implementation
        would violate illegal-states-unrepresentable.
        """
        op = AngularRedistribution()
        assert CAP_SOLVE not in op.capabilities

    @pytest.mark.foundation
    def test_cap_apply_transpose_NOT_declared_at_step_1(self):
        """``CAP_APPLY_TRANSPOSE`` MUST NOT be declared at Step 1.

        The M-M adjoint is the reverse-order recurrence — ``.H``
        propagation is Step 5 scope.
        """
        op = AngularRedistribution()
        assert CAP_APPLY_TRANSPOSE not in op.capabilities


# ═══════════════════════════════════════════════════════════════════════
# Gate 6 #5 — Linearity of AngularRedistribution.apply
# ═══════════════════════════════════════════════════════════════════════


class TestAngularRedistributionApplyLinearity:
    """Gate 6 #5 — ``AngularRedistribution.apply`` is linear in input.

    The M-M recurrence (Hébert §3.9.4 Eqs. 3.437 / 3.439) is an
    affine transformation in the cell-averaged scalar flux.  When
    the Carlson seed (linear in ``psi_cells``) is composed as the
    inner half-angle seed, the M-M recurrence is purely linear.

    Test: ``op.apply(α·ψ_a + β·ψ_b) ≡ α·op.apply(ψ_a) + β·op.apply(ψ_b)``
    at ``rtol=1e-12``.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("pole-mm-recurrence")
    @pytest.mark.parametrize(
        "geometry", ["sphere", "cylinder"],
        ids=["sphere", "cylinder"],
    )
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize(
        "alpha,beta",
        [(1.0, 1.0), (2.5, -1.0), (0.7, 0.3)],
        ids=["sum", "diff", "convex"],
    )
    def test_apply_linear_in_psi_input(
        self, geometry, n_groups, alpha, beta,
    ):
        """Linearity: ``op.apply(α·ψ_a + β·ψ_b) = α·apply(ψ_a) + β·apply(ψ_b)``.

        Input shape is ``(ng, M, nx)`` — the natural M-M recurrence
        input.  We feed ``carlson_context=None`` here so the M-M
        recurrence runs with the seed under the operator's bound
        :class:`CarlsonInwardSweep` strategy reading the same context
        regardless of input — but since context is None, the
        recurrence falls back to ZeroSeed-equivalent behaviour, which
        is unambiguously linear in ψ.
        """
        inputs = _build_inputs(geometry, nx=8)
        quad = inputs["quad"]
        nx = inputs["volume"].shape[0]
        M = quad.N

        op = AngularRedistribution()
        rng = np.random.default_rng(
            seed=20260512 + n_groups + hash((geometry, alpha, beta)) % 100000,
        )
        psi_a = rng.standard_normal((n_groups, M, nx))
        psi_b = rng.standard_normal((n_groups, M, nx))

        R_a = op.apply(
            psi_cells=psi_a,
            alpha_half=inputs["alpha_half"],
            redist_dAw=inputs["redist_dAw"],
            tau_mm=inputs["tau_mm"],
            volume=inputs["volume"],
            level_indices=inputs["level_indices"],
            carlson_context=None,
        )
        R_b = op.apply(
            psi_cells=psi_b,
            alpha_half=inputs["alpha_half"],
            redist_dAw=inputs["redist_dAw"],
            tau_mm=inputs["tau_mm"],
            volume=inputs["volume"],
            level_indices=inputs["level_indices"],
            carlson_context=None,
        )
        R_lin = op.apply(
            psi_cells=alpha * psi_a + beta * psi_b,
            alpha_half=inputs["alpha_half"],
            redist_dAw=inputs["redist_dAw"],
            tau_mm=inputs["tau_mm"],
            volume=inputs["volume"],
            level_indices=inputs["level_indices"],
            carlson_context=None,
        )
        np.testing.assert_allclose(
            R_lin, alpha * R_a + beta * R_b,
            rtol=1e-12, atol=1e-13,
            err_msg=(
                f"AngularRedistribution.apply not linear on {geometry} "
                f"ng={n_groups} (α={alpha}, β={beta})"
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate 6 #4 — Flat-flux closure (Σw = 2 Hébert convention)
# ═══════════════════════════════════════════════════════════════════════


class TestAngularRedistributionFlatFluxClosure:
    r"""Gate 6 #4 — flat-ψ closure of the Carlson seed via the operator API.

    Operator-level dual of the Phase F Gate 1.6 identity at
    :file:`tests/sn/test_phase_c_gates.py:577`.  With Σw = 2 (Hébert
    convention) and ``Q_bar = Σ_t · ψ_const``, the inward Carlson
    sweep reproduces ``ψ_const`` at every cell.

    Mode 7 caveat — this gate uses a flat-ψ ansatz (isotropic-in-μ
    by construction) which nulls angular variation.  The
    angularly-non-trivial companion is the Phase F twin-path defense
    gate at
    :file:`tests/sn/spatial/test_sweep_vs_apply_consistency.py`,
    which probes ``sphere_2g_3reg`` (heterogeneous MR multi-group);
    both must ship together per the gate-design memo.

    We pin the seed identity by accessing the operator's bound
    :class:`CarlsonInwardSweep` strategy directly: the operator's
    ``angular_sweep.psi_half_seed`` is the inner LinearOperator that
    Step 1 promotes alongside the M-M recurrence; the algebraic
    identity holds on the seed itself.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("phase-f-carlson-seed-source-driven")
    @pytest.mark.verifies("phase-f-q-bar-twin-forms")
    @pytest.mark.verifies("mm-weights")
    @pytest.mark.catches("ERR-026")
    @pytest.mark.parametrize(
        "geometry", ["sphere", "cylinder"],
        ids=["sphere", "cylinder"],
    )
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize("sigma_t_value", [0.5, 1.5])
    @pytest.mark.parametrize("psi_const", [1.0, 2.5])
    def test_flat_psi_carlson_seed_identity_via_operator(
        self, geometry, n_groups, sigma_t_value, psi_const,
    ):
        """Operator-level flat-ψ Carlson seed identity.

        Setup mirrors the Phase F Gate 1.6 test at
        :file:`tests/sn/test_phase_c_gates.py:577` — Σw = 2 Hébert
        convention, ``Q_bar = Σ_t · ψ_const``, reflective ``bc_outer
        = ψ_const`` → seed value ψ_const at every cell.

        We compose ``AngularRedistribution()`` with its default
        ``MorelMontryAngularSweep`` (which carries the default
        ``CarlsonInwardSweep`` seed) and access the seed via the
        operator's ``angular_sweep.psi_half_seed`` attribute — same
        composition path the operator uses at ``apply`` time.
        """
        inputs = _build_inputs(geometry, nx=8)
        nx = inputs["volume"].shape[0]

        # Σw = 2 convention: M=2 ordinates, weights [1, 1], μ = [-1, +1].
        M_apply = 2
        weights_apply = np.array([1.0, 1.0])
        mu_apply = np.array([-1.0, 1.0])
        psi_level_flat = np.full((n_groups, M_apply, nx), psi_const)
        sigma_t_gx = np.full((n_groups, nx), sigma_t_value)
        bc_outer = np.full((n_groups,), psi_const)
        ctx = CarlsonSweepContext(
            sigma_t=sigma_t_gx,
            dr=inputs["dr"],
            mu_quad=mu_apply,
            weights=weights_apply,
            bc_outer_value=bc_outer,
        )

        op = AngularRedistribution()
        # The bound Carlson seed strategy on the M-M sweep.
        seed_strategy = op.angular_sweep.psi_half_seed
        seed = seed_strategy(psi_level_flat, ctx)

        expected_const = np.full((n_groups, nx), psi_const)
        np.testing.assert_allclose(
            seed, expected_const, rtol=1e-13, atol=1e-13,
            err_msg=(
                f"Carlson seed flat-ψ identity broken on {geometry} "
                f"(ng={n_groups}, Σ_t={sigma_t_value}, ψ_const={psi_const}). "
                "With Σw=2 Hébert convention + reflective bc_outer=ψ_const, "
                "the inward Carlson sweep MUST reproduce ψ_const at every "
                "cell.  Regression would indicate the operator's bound "
                "CarlsonInwardSweep strategy drifted from the canonical "
                "Hébert §3.9.4 (3.434)-(3.435) algebra."
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate 6 #3 — Carlson seed equivalence (operator-level)
# ═══════════════════════════════════════════════════════════════════════


class TestAngularRedistributionCarlsonSeedEquivalence:
    r"""Gate 6 #3 — operator-level equivalence between the apply-path
    Carlson seed (``CarlsonInwardSweep.__call__``) and the
    source-driven helper (``carlson_inward_sweep_from_source``) when
    BOTH are composed inside ``AngularRedistribution``.

    The existing 57 function-level tests at
    :file:`tests/sn/spatial/test_sweep_vs_apply_consistency.py` pin
    equivalence at the BARE-FUNCTION level.  This class pins that the
    OPERATOR-level composition preserves the equivalence (the
    promotion didn't introduce a path that ignores the canonical
    helper).
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("phase-f-carlson-seed-source-driven")
    @pytest.mark.catches("ERR-026")
    @pytest.mark.parametrize("nx", [4, 8, 16])
    @pytest.mark.parametrize("sigma_t_value", [0.1, 0.5, 1.5])
    @pytest.mark.parametrize("bc_outer_value", [0.0, 0.3, 1.0])
    @pytest.mark.parametrize("psi_const", [1.0, 2.5])
    def test_operator_carlson_seed_matches_source_driven_helper(
        self, nx, sigma_t_value, bc_outer_value, psi_const,
    ):
        """``AngularRedistribution``'s bound Carlson seed strategy
        matches ``carlson_inward_sweep_from_source`` on matched
        flat-ψ inputs.

        Mirrors the parametrize matrix from the existing
        ``test_carlson_seed_apply_sweep_equivalence_flat_psi`` (57
        tests at ``test_sweep_vs_apply_consistency.py:46+``) — lifts
        each function-level assertion to the operator-level API.
        """
        sigma_t_gx = np.full((1, nx), sigma_t_value)
        dr = np.full(nx, 1.0 / nx)
        bc_outer = np.full((1,), bc_outer_value)

        # Apply path (operator-level): bound CarlsonInwardSweep on
        # AngularRedistribution's M-M sweep.
        M = 2
        weights = np.array([1.0, 1.0])
        mu = np.array([-1.0, 1.0])
        psi_level = np.full((1, M, nx), psi_const)
        ctx = CarlsonSweepContext(
            sigma_t=sigma_t_gx, dr=dr,
            mu_quad=mu, weights=weights,
            bc_outer_value=bc_outer,
        )
        op = AngularRedistribution()
        seed_via_op = op.angular_sweep.psi_half_seed(psi_level, ctx)

        # Source-driven path: the helper directly.
        Q_bar = np.full((1, nx), sigma_t_value * psi_const)
        seed_via_helper = carlson_inward_sweep_from_source(
            Q_bar=Q_bar, sigma_t=sigma_t_gx, dr=dr,
            bc_outer_value=bc_outer,
        )

        np.testing.assert_allclose(
            seed_via_op, seed_via_helper,
            rtol=1e-14, atol=1e-14,
            err_msg=(
                f"Operator-level Carlson seed ≠ source-driven helper "
                f"(nx={nx}, Σ_t={sigma_t_value}, bc_outer={bc_outer_value}, "
                f"ψ_const={psi_const}). The AngularRedistribution promotion "
                "must keep the seed strategy bound and consume it via "
                "the same algebraic path."
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate 2 — apply-solve round-trip — SKIPPED at Step 1 per Q4 decision
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.skip(
    reason=(
        "AngularRedistribution drops CAP_SOLVE at Step 1 per Pattern 4 "
        "(advertise only what works).  Back-substitution is deferred to "
        "a later step; round-trip test re-enables when CAP_SOLVE is "
        "advertised.  See open-questions §4 in the gate-design memo."
    ),
)
class TestAngularRedistributionApplySolveRoundTrip:
    """Gate 2 — DEFERRED at Step 1 (no ``CAP_SOLVE``)."""

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "geometry", ["sphere", "cylinder"],
        ids=["sphere", "cylinder"],
    )
    @pytest.mark.parametrize("n_groups", [1, 2])
    def test_round_trip_apply_solve_recovers_input(
        self, geometry, n_groups,
    ):
        """``op.apply(op.solve(ψ)) ≡ ψ`` — DEFERRED, no CAP_SOLVE."""
        pass


__all__: list[str] = []
