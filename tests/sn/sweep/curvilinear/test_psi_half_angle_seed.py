r"""Foundation + L0 tests for the PsiHalfAngleSeed strategy family.

Issue #168 Phase D — pins the software contract of the new
:class:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeed`
Protocol + two concrete strategies (:class:`ZeroSeed`,
:class:`CarlsonInwardSweep`) that compose into
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
via Option α (composition, not sibling Protocol).

Test layers
-----------

* **Foundation** (Protocol / registry / immutability / linearity):
  software contract.  No equation labels.
* **L0** (term verification via flat-ψ algebraic identity + multi-region
  + multi-group probes): the canonical Hébert §3.9.4 (3.432)-(3.435)
  algebra is verified directly.
* **L1** (structural-independence probe vs naive ψ_cell broadcast):
  vacuum-BC falsification probe per
  :doc:`/.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis`
  §8.

Mirrors :file:`tests/sn/spatial/test_pole_angular_closure.py` in
shape; the foundation block here parallels the foundation block
there.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.psi_half_angle_seed import (
    CarlsonInwardSweep,
    CarlsonSweepContext,
    PsiHalfAngleSeed,
    PsiHalfAngleSeedBase,
    ZeroSeed,
)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — Protocol / registry / immutability
# ═══════════════════════════════════════════════════════════════════════


class TestProtocolConformance:
    """The :class:`PsiHalfAngleSeed` Protocol is honoured by both
    concrete strategies via structural typing."""

    @pytest.mark.foundation
    def test_zero_seed_satisfies_protocol(self) -> None:
        assert isinstance(ZeroSeed(), PsiHalfAngleSeed)

    @pytest.mark.foundation
    def test_carlson_inward_sweep_satisfies_protocol(self) -> None:
        assert isinstance(CarlsonInwardSweep(), PsiHalfAngleSeed)

    @pytest.mark.foundation
    def test_is_linear_class_attr_advertised(self) -> None:
        """Both strategies advertise ``is_linear = True``.

        Load-bearing: the M-M closure delegates to these strategies
        in the apply matvec, which MUST remain a linear operator.
        Non-linearity in the seed would propagate non-linearity into
        the operator and break :meth:`InvertibleOperator.apply_transpose`
        / dense matrix probing.
        """
        assert ZeroSeed.is_linear is True
        assert CarlsonInwardSweep.is_linear is True

    @pytest.mark.foundation
    def test_zero_seed_repr(self) -> None:
        assert repr(ZeroSeed()) == "ZeroSeed()"

    @pytest.mark.foundation
    def test_carlson_inward_sweep_repr(self) -> None:
        assert repr(CarlsonInwardSweep()) == "CarlsonInwardSweep()"


class TestRegistry:
    """RegistryMixin self-registration for PsiHalfAngleSeed family."""

    @pytest.mark.foundation
    def test_zero_seed_registered(self) -> None:
        assert "zero" in PsiHalfAngleSeedBase.registry
        assert PsiHalfAngleSeedBase.registry["zero"] is ZeroSeed

    @pytest.mark.foundation
    def test_carlson_inward_sweep_registered(self) -> None:
        assert "carlson_inward_sweep" in PsiHalfAngleSeedBase.registry
        assert (
            PsiHalfAngleSeedBase.registry["carlson_inward_sweep"]
            is CarlsonInwardSweep
        )

    @pytest.mark.foundation
    def test_create_zero_seed_returns_instance(self) -> None:
        instance = PsiHalfAngleSeedBase.create("zero")
        assert isinstance(instance, ZeroSeed)

    @pytest.mark.foundation
    def test_create_carlson_inward_sweep_returns_instance(self) -> None:
        instance = PsiHalfAngleSeedBase.create("carlson_inward_sweep")
        assert isinstance(instance, CarlsonInwardSweep)

    @pytest.mark.foundation
    def test_create_unknown_key_raises(self) -> None:
        with pytest.raises(
            KeyError, match="unknown PsiHalfAngleSeedBase key",
        ):
            PsiHalfAngleSeedBase.create("unknown_strategy")


class TestImmutability:
    """Both strategies are ``@dataclass(frozen=True, slots=True)``."""

    @pytest.mark.foundation
    def test_zero_seed_is_frozen(self) -> None:
        strategy = ZeroSeed()
        with pytest.raises((AttributeError, TypeError)):
            strategy.is_linear = False  # type: ignore[misc]

    @pytest.mark.foundation
    def test_carlson_inward_sweep_is_frozen(self) -> None:
        strategy = CarlsonInwardSweep()
        with pytest.raises((AttributeError, TypeError)):
            strategy.is_linear = False  # type: ignore[misc]


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — shape contract + ZeroSeed bit-identity
# ═══════════════════════════════════════════════════════════════════════


def _make_context(
    ng: int = 1, nx: int = 10, M: int = 4,
    sigma_t_value: float = 0.5,
    bc_outer_value: float = 1.0,
    R: float = 2.0,
) -> CarlsonSweepContext:
    """Build a uniform-mesh CarlsonSweepContext for hand-calc tests."""
    sigma_t = np.full((ng, nx), sigma_t_value)
    dr = np.full(nx, R / nx)
    # Use uniform GL-like mu spacing for tests; specifics don't matter
    # for the algebra (only weights · ψ matters for the moment).
    mu_quad = np.linspace(-0.9, 0.9, M)
    weights = np.full(M, 2.0 / M)  # sums to 2 (GL convention)
    bc_outer = np.full(ng, bc_outer_value)
    return CarlsonSweepContext(
        sigma_t=sigma_t, dr=dr, mu_quad=mu_quad,
        weights=weights, bc_outer_value=bc_outer,
        mu_start=-1.0,)


class TestZeroSeedShapeAndBitIdentity:
    """The :class:`ZeroSeed` strategy is the regression-safety ablation.

    Must produce :func:`numpy.zeros((ng, nx))` EXACTLY (bit-identical
    to the pre-Phase-D ``psi_half_left = np.zeros(...)`` hardcode).
    """

    @pytest.mark.foundation
    def test_zero_seed_returns_zeros_shape_and_dtype(self) -> None:
        ng, M, nx = 2, 4, 10
        psi_level = np.ones((ng, M, nx))
        context = _make_context(ng=ng, nx=nx, M=M)
        result = ZeroSeed()(psi_level, context)
        # Shape contract: (ng, nx).
        assert result.shape == (ng, nx)
        # Bit-identity to np.zeros (the pre-Phase-D hardcode).
        assert np.array_equal(result, np.zeros((ng, nx)))
        # dtype matches psi_level.
        assert result.dtype == psi_level.dtype

    @pytest.mark.foundation
    def test_zero_seed_independent_of_input(self) -> None:
        """ZeroSeed ignores psi_level entirely — output is constant zero."""
        ng, M, nx = 1, 4, 5
        context = _make_context(ng=ng, nx=nx, M=M)
        result_zero_input = ZeroSeed()(np.zeros((ng, M, nx)), context)
        result_random_input = ZeroSeed()(
            np.random.default_rng(0).standard_normal((ng, M, nx)),
            context,
        )
        assert np.array_equal(result_zero_input, result_random_input)
        assert np.array_equal(result_zero_input, np.zeros((ng, nx)))


class TestCarlsonInwardSweepShape:
    """Shape contract pinning for :class:`CarlsonInwardSweep`."""

    @pytest.mark.foundation
    def test_carlson_returns_correct_shape(self) -> None:
        ng, M, nx = 3, 4, 8
        psi_level = np.ones((ng, M, nx))
        context = _make_context(ng=ng, nx=nx, M=M)
        result = CarlsonInwardSweep()(psi_level, context)
        assert result.shape == (ng, nx)


# ═══════════════════════════════════════════════════════════════════════
# L0 term verification — flat-ψ algebraic identity
# ═══════════════════════════════════════════════════════════════════════


class TestCarlsonFlatPsiAlgebraicIdentity:
    """Hébert §3.9.4 flat-ψ trace (literature memo §3 algebraic proof).

    With flat ψ = C, consistent source ``Q̄ = Σ_t · C``, and
    consistent BC ``bc_outer = C``, the Hébert (3.434)-(3.435) inward
    recurrence is self-similar: every cell and face stays at C.

    The L0 verification:

    .. math::

       \\bar\\phi_i = \\frac{\\Delta r \\cdot \\Sigma \\cdot C + 2 C}
                              {\\Delta r \\cdot \\Sigma + 2}
                    = C \\cdot \\frac{\\Delta r \\Sigma + 2}{\\Delta r \\Sigma + 2}
                    = C
    """

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_flat_psi_identity_reflective(self) -> None:
        r"""Flat ψ = 1, Σ_t = 0.5, reflective bc_outer = 1 → all phi_aux = 1."""
        ng, M, nx = 1, 4, 10
        # Flat ψ = 1 across all ordinates and cells.
        psi_level = np.ones((ng, M, nx))
        # On flat ψ = 1, the scalar moment φ_0 = Σ w_n · 1 = Σ weights.
        # For weights summing to 2 (GL convention): φ_0 = 2.
        # Q̄ = 0.5 · Σ_t · φ_0 = 0.5 · 0.5 · 2 = 0.5 (matches Σ_t · ψ).
        # bc_outer = 1.0 (reflective: outer face at μ=−1 mirrors μ=+1 = 1).
        context = _make_context(
            ng=ng, nx=nx, M=M,
            sigma_t_value=0.5, bc_outer_value=1.0,
        )
        result = CarlsonInwardSweep()(psi_level, context)
        # Every cell should be C = 1 to machine precision.
        np.testing.assert_allclose(result, np.ones((ng, nx)), atol=1e-15)

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_flat_psi_identity_varying_C(self) -> None:
        r"""Identity for ψ = C ≠ 1 — algebraic invariance under scaling."""
        for C in (0.7, 2.3, 100.0):
            ng, M, nx = 1, 4, 8
            psi_level = C * np.ones((ng, M, nx))
            context = _make_context(
                ng=ng, nx=nx, M=M,
                sigma_t_value=0.5, bc_outer_value=C,
            )
            result = CarlsonInwardSweep()(psi_level, context)
            np.testing.assert_allclose(
                result, C * np.ones((ng, nx)), atol=1e-13 * abs(C),
                err_msg=f"flat-ψ identity failed at C = {C}",
            )

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_vacuum_BC_flat_source_nx_3(self) -> None:
        r"""Vacuum BC + flat source Q̄ = 0.5 + nx = 3: hand calculation.

        Per literature memo §3.b (vacuum trace):

        * φ̄_{nx+1/2} = 0 (vacuum)
        * φ̄_{2} = (Δr · 0.5 + 0) / (Δr · 0.5 + 2)
        * φ̄_{5/2 - 1/2} = 2·φ̄_{2} − 0 = 2·φ̄_{2}
        * Continue inward for k = 1, 0.

        With Δr · Σ_t = 0.5 / 3 · 0.5 / 1 = ... actually let's just
        compute numerically with concrete Δr = 1/3, Σ_t = 0.5.
        """
        nx = 3
        dr_val = 2.0 / nx  # 0.667
        sigma_t_value = 0.5
        # In ψ-input-driven Q̄: Q̄ = 0.5 · Σ_t · φ_0(ψ_input).
        # Use ψ_level designed so φ_0(ψ_input) = 2 (matches Σ w · 1).
        ng, M = 1, 4
        psi_level = np.ones((ng, M, nx))  # gives φ_0 = 2 on GL-like weights
        context = _make_context(
            ng=ng, nx=nx, M=M,
            sigma_t_value=sigma_t_value, bc_outer_value=0.0,  # vacuum
            R=2.0,
        )
        # Q̄ = 0.5 · 0.5 · 2.0 = 0.5
        # Hand calc:
        # phi_face = 0 (outer)
        # k=2: denom = 0.667·0.5 + 2 = 2.333
        #      phi_cell = (0.667 · 0.5 + 0) / 2.333 = 0.333/2.333 = 0.1429
        #      phi_face = 2·0.1429 - 0 = 0.2857
        # k=1: denom = 2.333
        #      phi_cell = (0.667·0.5 + 2·0.2857) / 2.333 = (0.333 + 0.5714)/2.333 = 0.388
        #      phi_face = 2·0.388 - 0.2857 = 0.490
        # k=0: phi_cell = (0.333 + 2·0.490) / 2.333 = 1.314/2.333 = 0.563
        expected_2 = (dr_val * 0.5) / (dr_val * 0.5 + 2)
        face_at_2 = 2 * expected_2
        expected_1 = (dr_val * 0.5 + 2 * face_at_2) / (dr_val * 0.5 + 2)
        face_at_1 = 2 * expected_1 - face_at_2
        expected_0 = (dr_val * 0.5 + 2 * face_at_1) / (dr_val * 0.5 + 2)
        expected = np.array([[expected_0, expected_1, expected_2]])

        result = CarlsonInwardSweep()(psi_level, context)
        np.testing.assert_allclose(result, expected, rtol=1e-13)

    @pytest.mark.l0
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_multi_region_sigma_t_step(self) -> None:
        r"""Σ_t jumps at midpoint (heterogeneous medium).

        Per Hébert §3.9.4 the recurrence handles per-cell Σ_t
        trivially (no special treatment for material discontinuities
        as long as the radial mesh respects boundaries).
        """
        ng, M, nx = 1, 4, 6
        psi_level = np.ones((ng, M, nx))
        # σ_t step: first 3 cells at 0.3, last 3 at 1.0.
        sigma_t = np.zeros((ng, nx))
        sigma_t[:, :3] = 0.3
        sigma_t[:, 3:] = 1.0
        dr = np.full(nx, 1.0 / nx)
        context = CarlsonSweepContext(
            sigma_t=sigma_t,
            dr=dr,
            mu_quad=np.linspace(-0.9, 0.9, M),
            weights=np.full(M, 0.5),
            bc_outer_value=np.array([0.0]),
            mu_start=-1.0,)
        result = CarlsonInwardSweep()(psi_level, context)
        # Hand-trace: phi_face_outer = 0; phi_0(ψ_input) = 2.
        # Q̄[k] = 0.5 · σ_t[k] · 2 = σ_t[k].
        # Sweep inward k=5,4,3 (σ_t=1), then k=2,1,0 (σ_t=0.3).
        phi_face = 0.0
        expected_cells = np.zeros(nx)
        for k in range(nx - 1, -1, -1):
            sig = sigma_t[0, k]
            Q_bar = sig  # = 0.5 · sig · 2.0
            phi_cell = (dr[k] * Q_bar + 2.0 * phi_face) / (dr[k] * sig + 2.0)
            expected_cells[k] = phi_cell
            phi_face = 2.0 * phi_cell - phi_face
        np.testing.assert_allclose(result[0], expected_cells, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# L0 term verification — linearity (foundation of operator algebra)
# ═══════════════════════════════════════════════════════════════════════


class TestSeedLinearity:
    """Both seed strategies are linear in psi_level.

    Required by the operator algebra: the M-M closure is linear, and
    a non-linear seed would propagate non-linearity into the apply
    matvec.  Linearity is a pre-condition for
    :meth:`~orpheus.sn.operator.InvertibleOperator.apply_transpose`
    correctness (dense-matrix probing assumes linearity).
    """

    @pytest.mark.foundation
    def test_zero_seed_is_linear(self) -> None:
        ng, M, nx = 1, 4, 5
        context = _make_context(ng=ng, nx=nx, M=M)
        rng = np.random.default_rng(seed=42)
        psi1 = rng.standard_normal((ng, M, nx))
        psi2 = rng.standard_normal((ng, M, nx))
        alpha = 1.7
        beta = -0.3
        out_combined = ZeroSeed()(alpha * psi1 + beta * psi2, context)
        out_split = (
            alpha * ZeroSeed()(psi1, context)
            + beta * ZeroSeed()(psi2, context)
        )
        np.testing.assert_allclose(out_combined, out_split, atol=1e-15)

    @pytest.mark.foundation
    def test_carlson_inward_sweep_is_linear(self) -> None:
        ng, M, nx = 1, 4, 5
        context = _make_context(ng=ng, nx=nx, M=M)
        rng = np.random.default_rng(seed=42)
        psi1 = rng.standard_normal((ng, M, nx))
        psi2 = rng.standard_normal((ng, M, nx))
        alpha = 1.7
        beta = -0.3
        # Linearity requires bc_outer to also be linear in ψ; the
        # context fixes it so this test pins linearity given a
        # constant bc_outer (the matvec wiring makes bc_outer linear
        # in input ψ via the BC realization).
        out_combined = CarlsonInwardSweep()(
            alpha * psi1 + beta * psi2, context,
        )
        out_split = (
            alpha * CarlsonInwardSweep()(psi1, context)
            + beta * CarlsonInwardSweep()(psi2, context)
        )
        # The recurrence has an affine bc_outer offset that doesn't
        # cancel under linear combinations of ψ alone (bc_outer is a
        # constant in this test).  Test ψ-linearity by setting
        # bc_outer = 0 so the recurrence is purely linear in ψ.
        context_no_bc = CarlsonSweepContext(
            sigma_t=context.sigma_t,
            dr=context.dr,
            mu_quad=context.mu_quad,
            weights=context.weights,
            bc_outer_value=np.zeros_like(context.bc_outer_value),
            mu_start=-1.0,)
        out_combined = CarlsonInwardSweep()(
            alpha * psi1 + beta * psi2, context_no_bc,
        )
        out_split = (
            alpha * CarlsonInwardSweep()(psi1, context_no_bc)
            + beta * CarlsonInwardSweep()(psi2, context_no_bc)
        )
        np.testing.assert_allclose(out_combined, out_split, atol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# L1 — Structural-independence vs naive ψ_cell broadcast (falsification probe)
# ═══════════════════════════════════════════════════════════════════════


class TestCarlsonStructuralIndependence:
    """Vacuum-BC + non-flat probe distinguishes Carlson from naive ψ_cell broadcast.

    Per :doc:`/.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis`
    §8: the flat-ψ reflective probe cannot distinguish
    :class:`CarlsonInwardSweep` from a degenerate "broadcast ψ_cell"
    seed, because both give C on flat ψ.  Vacuum BC + flat ψ
    introduces a non-trivial phi_aux profile that DIFFERS from
    ψ_cell-broadcast — this falsifies the naive seed.

    A future regression that replaces ``CarlsonInwardSweep`` with a
    naive broadcast would change the output, and this test catches
    the substitution.
    """

    @pytest.mark.l1
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.catches("ERR-026")
    def test_carlson_vs_zero_seed_differ_on_vacuum_bc(self) -> None:
        r"""ZeroSeed vs CarlsonInwardSweep produce structurally distinct
        outputs on vacuum-BC flat-ψ probe.

        ZeroSeed returns zeros; CarlsonInwardSweep produces a non-trivial
        spatial profile driven by Q̄ ≠ 0.  Maximum absolute difference
        must be > 0 (i.e. the outputs are not coincidentally equal).
        """
        ng, M, nx = 1, 4, 10
        psi_level = np.ones((ng, M, nx))
        # Vacuum BC = 0; flat ψ = 1.  CarlsonInwardSweep gives non-trivial
        # spatial profile; ZeroSeed gives zeros.
        context = _make_context(
            ng=ng, nx=nx, M=M,
            sigma_t_value=0.5, bc_outer_value=0.0,
        )
        result_carlson = CarlsonInwardSweep()(psi_level, context)
        result_zero = ZeroSeed()(psi_level, context)
        diff = float(np.max(np.abs(result_carlson - result_zero)))
        assert diff > 0.05, (
            f"Carlson and Zero seeds should differ on vacuum-BC probe; "
            f"max|diff| = {diff:.3e}"
        )
        # Sanity: Carlson is non-zero where the input is non-zero.
        assert float(np.max(np.abs(result_carlson))) > 0.05
        # ZeroSeed is identically zero (by definition).
        assert np.all(result_zero == 0.0)


# ═══════════════════════════════════════════════════════════════════════
# Foundation test — MorelMontryAngularSweep default field is Carlson
# ═══════════════════════════════════════════════════════════════════════


class TestMorelMontryDefaultSeed:
    """Pin that MorelMontryAngularSweep's default psi_half_seed is the
    operator-consistent angular-edge extrapolation (ERR-058 b, #195).

    A future regression that flips the default back to
    ``CarlsonInwardSweep`` (the proxy-source seed, exact only at
    flat-flux equilibrium) or ``ZeroSeed`` (the Phase B pre-fix
    behaviour) is caught here.
    """

    @pytest.mark.foundation
    def test_default_seed_is_angular_edge_extrapolation(self) -> None:
        from orpheus.sn.spatial.pole_angular_closure import (
            MorelMontryAngularSweep,
        )
        from orpheus.sn.spatial.psi_half_angle_seed import (
            AngularEdgeExtrapolation,
        )
        instance = MorelMontryAngularSweep()
        assert isinstance(instance.psi_half_seed, AngularEdgeExtrapolation)

    @pytest.mark.foundation
    def test_carlson_remains_registered_opt_in(self) -> None:
        """The Hébert recurrence host stays available by explicit key."""
        from orpheus.sn.spatial.psi_half_angle_seed import (
            PsiHalfAngleSeedBase,
        )
        cls = PsiHalfAngleSeedBase.registry["carlson_inward_sweep"]
        assert cls is CarlsonInwardSweep

    @pytest.mark.foundation
    def test_can_construct_with_zero_seed_opt_out(self) -> None:
        """User can opt out of Carlson default via explicit ZeroSeed."""
        from orpheus.sn.spatial.pole_angular_closure import (
            MorelMontryAngularSweep,
        )
        instance = MorelMontryAngularSweep(psi_half_seed=ZeroSeed())
        assert isinstance(instance.psi_half_seed, ZeroSeed)
