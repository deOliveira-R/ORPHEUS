r"""Foundation / L0 tests for :class:`MultiplicationOperator` (#257 S3b).

The §5.7 promotion ``M[f]`` of a
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
to a diagonal :class:`~orpheus.numerics.operator.LinearOperator`. Two
test families:

1. **Broadcast oracle** — the generalised N-D
   :class:`~orpheus.numerics.operator.DiagonalOperator` engine on the
   σ_t leading-ordinate form reduces EXACTLY to the legacy
   ``sigma[None] * psi`` at the VALUES level. The discriminating regime
   is a 2-D mesh carrier ``(N_ord, ng, nx, ny)`` with **nx ≠ ny** and
   **ng = 2** — the only shape where a transposed broadcast axis would
   mis-scale rather than silently agree (variable-swap failure mode #2).
   Primary assertion :func:`numpy.testing.assert_array_equal` (it IS the
   same operation); a ``nulp=1`` fallback is documented only against a
   future re-association (per the test-architect spec, gate item 1).

2. **Multiplier-algebra law-suite** (the intrinsic-properties directive)
   — the laws that make "a coefficient IS a diagonal operator" literally
   true: :math:`M[1]=I`, :math:`M[0]=0` (codomain-aware: collision's zero
   is a SOURCE, not a flux), linearity :math:`M[af+bg]=aM[f]+bM[g]`
   (tested on ≥2G ASYMMETRIC HETEROGENEOUS coefficients — anti-pattern
   #3/#4: a 1G / homogeneous coefficient is blind here), self-adjointness
   :math:`M[f]^*=M[f]`, the spectrum law
   :math:`\mathrm{spec}(M[f])=\mathrm{ess\,range}(f)` (``CAP_SOLVE`` iff
   ``min|f|>0``), and the homomorphism :math:`M[f]M[g]=M[fg]` tested at
   the VALUES level via the numerics engine on the raw product array
   (``σ·σ'`` has units ``cm⁻²``, so the homomorphism is a numerics-engine
   property — built directly, NOT as a ``CrossSectionField``).

The resolvent-assembly leg (``L + M[σ] → InvertibleOperator``, the
sweep/matvec matching the legacy resolvent and the analytical
``k_∞ = νΣf/Σa`` / streaming-equilibrium references) is verified by the
EXISTING structurally-independent gates that stay green:
``tests/sn/verification/analytical/test_kinf_homogeneous.py`` and
``tests/sn/operators/test_invertible_operator.py`` (the
``test_si_carve_recovers_analytical_kinf`` + the
``..._recovers_q_over_sigma_composite`` streaming-equilibrium leg). This
file adds an L0 bare-engine streaming-equilibrium values check so the
broadcast itself is pinned without a full solve.

``foundation`` — software invariants (the multiplier-algebra laws), not
a theory-page equation claim.

References
----------

* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S3b.
* ``.claude/agent-memory/test-architect/issue_257_s3_multiplication_operator_verification.md``
  — the 5 gate verdicts + failure-mode exposure.
* Grand Report v3 §5.5–5.7 (the operator promotion).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    DiagonalOperator,
    IdentityOperator,
    MissingCapability,
    ZeroOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.multiplication_operator import MultiplicationOperator
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
# strips bare ``assert`` to a NO-OP. These structural gates (isinstance /
# capability-set / membership) MUST fire under -O, so they route through
# ``pytest.fail`` (a function call) rather than a bare ``assert`` (a
# false-green tripwire under -O).
def _require(condition: bool, message: str) -> None:
    """A -O-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — a 1-D slab and the discriminating 2-D nx≠ny carrier.
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cartesian_2d_mesh(nx: int = 5, ny: int = 3, ng: int = 2) -> SNMesh:
    """The discriminating regime: nx ≠ ny, ng = 2.

    A transposed broadcast axis (the N-D generalisation's NEW risk —
    variable-swap mode #2) mis-scales or raises ONLY when the carrier
    has ≥2 complement axes of DIFFERENT sizes. nx ≠ ny + ng = 2 makes a
    wrong ``broadcast_axes`` observable; a square mesh would silently
    agree.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _random_state(sn_mesh: SNMesh, ng: int = 2, seed: int = 42) -> TimedFullField:
    """Random :class:`TimedFullField` whose bulk has shape ``(N, ng, *spatial)``."""
    rng = np.random.default_rng(seed)
    N = sn_mesh.quad.N
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    return replace(
        state,
        bulk=replace(
            state.bulk,
            values=rng.standard_normal((N, ng, *sn_mesh.spatial_shape)),
        ),
    )


def _positive_sigma(sn_mesh: SNMesh, ng: int = 2, seed: int = 11) -> np.ndarray:
    """Heterogeneous positive σ ``(ng, *spatial)``, bounded away from 0."""
    rng = np.random.default_rng(seed)
    return 0.3 + 0.5 * rng.random((ng, *sn_mesh.spatial_shape))


def _multiplier(sn_mesh: SNMesh, sigma: np.ndarray) -> MultiplicationOperator:
    """``M[σ]`` from a raw ndarray (wrapped into a CrossSectionField)."""
    return MultiplicationOperator(
        coefficient=CrossSectionField.from_mesh(sigma, sn_mesh),
    )


# ═══════════════════════════════════════════════════════════════════════
# 1. Broadcast oracle — the generalised engine ≡ legacy sigma[None]*psi.
# ═══════════════════════════════════════════════════════════════════════


class TestBroadcastOracle:
    """The σ_t leading-ordinate broadcast reduces EXACTLY to ``σ[None]·ψ``."""

    def test_2d_nxneqny_engine_equals_legacy_broadcast(self):
        """nx≠ny, ng=2: the engine result IS ``sigma[None] * psi`` at VALUES level.

        The discriminating regime for axis-ordering bugs (variable-swap
        mode #2). ``assert_array_equal`` is the PRIMARY gate — both forms
        are ``expand_dims`` on the leading axis, so they are the SAME
        operation (0 ULP expected, reduction_depth = 1 pure
        broadcast-multiply). A ``nulp=1`` fallback would be documented
        only against a future re-association; none exists today.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = _positive_sigma(sn, ng=2, seed=101)
        psi = _random_state(sn, ng=2, seed=202)

        M = _multiplier(sn, sigma)
        engine_result = M.apply(psi).bulk.values
        legacy_result = sigma[None] * psi.bulk.values

        np.testing.assert_array_equal(engine_result, legacy_result)

    def test_2d_codomain_is_source_not_flux(self):
        """The codomain of ``M[σ]`` is a SOURCE (collision rate), not a flux."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        M = _multiplier(sn, _positive_sigma(sn, ng=2, seed=303))
        out = M.apply(_random_state(sn, ng=2, seed=404))
        _require(
            isinstance(out.bulk, AngularSourceSink),
            f"M[σ].apply codomain bulk must be AngularSourceSink (a SOURCE), "
            f"got {type(out.bulk).__name__}",
        )
        # Multiplier has no face-trace action → boundary is implicit zero.
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    def test_solve_codomain_returns_to_flux(self):
        """``M[σ]^{-1} q = q/σ`` returns a flux (inverse of flux→source)."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = _positive_sigma(sn, ng=2, seed=505)
        q = _random_state(sn, ng=2, seed=606)
        M = _multiplier(sn, sigma)
        out = M.solve(q)
        _require(
            isinstance(out.bulk, AngularFlux),
            f"M[σ].solve codomain bulk must be AngularFlux (a flux), "
            f"got {type(out.bulk).__name__}",
        )
        np.testing.assert_array_equal(out.bulk.values, q.bulk.values / sigma[None])

    def test_round_trip_solve_apply_identity(self):
        """``solve(apply(ψ)) == ψ`` for a coefficient bounded away from 0."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = _positive_sigma(sn, ng=2, seed=707)
        psi = _random_state(sn, ng=2, seed=808)
        M = _multiplier(sn, sigma)
        round_trip = M.solve(M.apply(psi))
        np.testing.assert_allclose(
            round_trip.bulk.values, psi.bulk.values, rtol=1e-14, atol=1e-15,
        )


# ═══════════════════════════════════════════════════════════════════════
# 2. The multiplier-algebra law-suite (intrinsic-properties directive).
# ═══════════════════════════════════════════════════════════════════════


class TestMultiplierAlgebraLaws:
    r"""The §5.7 embedding ``M : L^∞ → B(L²)`` is a faithful unital
    *-homomorphism onto the diagonal subalgebra.
    """

    def test_M_one_is_identity(self):
        r""":math:`M[1] = I` — the constant-one coefficient acts as the identity."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        ones = np.ones((2, *sn.spatial_shape))
        M1 = _multiplier(sn, ones)
        psi = _random_state(sn, ng=2, seed=909)
        I = IdentityOperator()
        # On the VALUES level (the engine domain): M[1]·x == I·x == x.
        np.testing.assert_array_equal(
            M1.apply(psi).bulk.values, I.apply(psi.bulk.values),
        )

    def test_M_zero_is_codomain_aware_zero_operator(self):
        r""":math:`M[0] = 0` — the zero coefficient promotes to a ZeroOperator.

        Codomain-aware: collision's zero is a *zero SOURCE*, not a zero
        flux. The structurally-independent reference is the numerics
        :class:`ZeroOperator` with the codomain-zero hook returning the
        SAME zero AngularSourceSink composite ``M[0].apply`` produces.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        zeros = np.zeros((2, *sn.spatial_shape))
        M0 = _multiplier(sn, zeros)
        psi = _random_state(sn, ng=2, seed=111)

        out = M0.apply(psi)
        # The zero coefficient maps every flux to the zero SOURCE.
        _require(
            isinstance(out.bulk, AngularSourceSink),
            f"M[0].apply codomain bulk must be AngularSourceSink (a zero "
            f"SOURCE, not a zero flux), got {type(out.bulk).__name__}",
        )
        np.testing.assert_array_equal(out.bulk.values, 0.0)

        # Structurally-independent: a codomain-aware ZeroOperator whose
        # codomain_zero builds the zero source-composite from psi.
        def _zero_source(p: TimedFullField) -> TimedFullField:
            return replace(
                p,
                bulk=AngularSourceSink.from_mesh(
                    np.zeros_like(p.bulk.values), p.bulk.mesh,
                ),
            )

        Z = ZeroOperator(codomain_zero=_zero_source)
        np.testing.assert_array_equal(
            out.bulk.values, Z.apply(psi).bulk.values,
        )
        # M[0] is NOT invertible — its spectrum is {0}.
        _require(
            CAP_SOLVE not in M0.capabilities,
            "M[0] must NOT advertise CAP_SOLVE (spectrum is {0}); "
            f"got {sorted(M0.capabilities)}",
        )

    def test_linearity_M_af_plus_bg(self):
        r""":math:`M[a f + b g] = a\,M[f] + b\,M[g]` on ≥2G asymmetric het.

        The coefficient fields ``f``, ``g`` are 2-group HETEROGENEOUS and
        per-group ASYMMETRIC (NOT 1G / homogeneous — anti-pattern #3/#4,
        which would null the group/spatial structure). Verified at the
        VALUES level: the promotion is a linear map on the coefficient
        vector space.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        f = _positive_sigma(sn, ng=2, seed=131)
        g = _positive_sigma(sn, ng=2, seed=141)
        # Asymmetry between the two groups (so a g0↔g1 swap is observable).
        f[1] *= 2.5
        g[0] *= 0.4
        a, b = 1.7, -0.6

        psi = _random_state(sn, ng=2, seed=151)
        M_combo = _multiplier(sn, a * f + b * g)
        M_f = _multiplier(sn, f)
        M_g = _multiplier(sn, g)

        lhs = M_combo.apply(psi).bulk.values
        # a·M[f]·ψ + b·M[g]·ψ at the values level.
        rhs = (
            a * M_f.apply(psi).bulk.values
            + b * M_g.apply(psi).bulk.values
        )
        np.testing.assert_allclose(lhs, rhs, rtol=1e-14, atol=1e-15)

    def test_self_adjoint_M_H_equals_M(self):
        r""":math:`M[f]^* = M[f]` — a real coefficient is self-adjoint.

        ``apply_transpose == apply``, and the metric-blind Euclidean
        ``.H`` (domain ``None``) reduces to the representation transpose,
        which equals ``apply``.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        M = _multiplier(sn, _positive_sigma(sn, ng=2, seed=161))
        psi = _random_state(sn, ng=2, seed=171)

        forward = M.apply(psi).bulk.values
        np.testing.assert_array_equal(
            M.apply_transpose(psi).bulk.values, forward,
        )
        # .H action (metric-blind, Euclidean) equals the forward action.
        np.testing.assert_array_equal(M.H.apply(psi).bulk.values, forward)

    def test_spectrum_cap_solve_iff_min_abs_positive(self):
        r""":math:`\mathrm{spec}(M[f]) = \mathrm{ess\,range}(f)`.

        ``CAP_SOLVE`` advertised iff ``min|f| > 0``. A coefficient with a
        zero entry has 0 in its spectrum → solve REVOKED (raises
        MissingCapability), while apply still works.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        psi = _random_state(sn, ng=2, seed=181)

        # Bounded away from 0 → invertible.
        pos = _positive_sigma(sn, ng=2, seed=191)
        M_pos = _multiplier(sn, pos)
        _require(
            M_pos.capabilities
            == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE}),
            f"min|f|>0 must advertise all three caps; got "
            f"{sorted(M_pos.capabilities)}",
        )

        # A single zero entry → 0 in the spectrum → solve revoked.
        with_zero = pos.copy()
        with_zero[1, 2, 1] = 0.0
        M_zero = _multiplier(sn, with_zero)
        _require(
            M_zero.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE}),
            f"a coefficient with a zero entry must drop CAP_SOLVE; got "
            f"{sorted(M_zero.capabilities)}",
        )
        _require(
            CAP_SOLVE not in M_zero.capabilities,
            "0 in the spectrum must revoke CAP_SOLVE",
        )
        with pytest.raises(MissingCapability):
            M_zero.solve(psi)
        # apply still works (the multiplier is defined; only the inverse
        # is undefined at the zero).
        np.testing.assert_array_equal(
            M_zero.apply(psi).bulk.values, with_zero[None] * psi.bulk.values,
        )

    def test_homomorphism_M_f_compose_M_g_equals_M_fg(self):
        r""":math:`M[f]\,M[g] = M[f\,g]` at the VALUES level.

        The composition of two multipliers is the multiplier of the
        product. Tested via the numerics engine on the raw arrays — the
        product ``σ·σ'`` has units ``cm⁻²`` (NOT a cross section), so the
        homomorphism is a property of the engine, built directly on the
        product array (NOT a ``CrossSectionField``, which would mis-state
        the units). The composite carrier round-trips through
        AngularSourceSink between the two applies, so the engine product
        is the structurally-independent reference.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        f = _positive_sigma(sn, ng=2, seed=211)
        g = _positive_sigma(sn, ng=2, seed=221)
        f[1] *= 1.9  # 2G asymmetry.
        psi = _random_state(sn, ng=2, seed=231)

        # M[f] ∘ M[g] applied via the two engines (values level).
        eng_g = DiagonalOperator(g, broadcast_axes=(0,))
        eng_f = DiagonalOperator(f, broadcast_axes=(0,))
        composed = eng_f.apply(eng_g.apply(psi.bulk.values))

        # M[f·g] — the engine on the raw product array (units cm⁻², so
        # NOT a CrossSectionField; the homomorphism is a numerics fact).
        eng_fg = DiagonalOperator(f * g, broadcast_axes=(0,))
        product = eng_fg.apply(psi.bulk.values)

        np.testing.assert_allclose(composed, product, rtol=1e-14, atol=1e-15)


# ═══════════════════════════════════════════════════════════════════════
# 3. L0 streaming-equilibrium — the bare-engine values check.
#
# The full resolvent-assembly leg (L + M[σ] → InvertibleOperator, sweep
# matching kinf / streaming-equilibrium) is the EXISTING gate set that
# stays green:
#   tests/sn/verification/analytical/test_kinf_homogeneous.py
#   tests/sn/operators/test_invertible_operator.py::
#     TestInvertibleSolveBridgeRegression::{test_si_carve_recovers_analytical_kinf,
#       test_fixed_source_homogeneous_reflective_recovers_q_over_sigma_composite}
# This pins the BARE multiplier algebra without a solve.
# ═══════════════════════════════════════════════════════════════════════


class TestStreamingEquilibriumValuesLeg:
    r"""On a flat-σ collision-only balance ``C ψ = Q``, ``ψ = Q/σ``.

    A degenerate streaming-equilibrium check at the multiplier level: with
    no streaming term the within-group balance is purely ``M[σ_t]ψ = Q``,
    whose solution is ``ψ = M[σ_t]^{-1}Q = Q/σ_t``. This is the L0 leg
    that the bare engine is correct; the full ``φ = Q/(Σ_t(1−c))``
    streaming-equilibrium with scattering is the InvertibleOperator gate.
    """

    def test_collision_only_balance_recovers_q_over_sigma(self):
        sn = _slab_mesh(nx=4, ng=2)
        sigma_t = 0.8 * np.ones((2, *sn.spatial_shape))  # flat Σ_t.
        Q = _random_state(sn, ng=2, seed=241)
        M = _multiplier(sn, sigma_t)
        psi = M.solve(Q)
        # M[σ_t]·ψ == Q (the balance is satisfied).
        np.testing.assert_allclose(
            M.apply(psi).bulk.values, Q.bulk.values, rtol=1e-14, atol=1e-15,
        )
        # And ψ == Q/σ_t explicitly.
        np.testing.assert_array_equal(
            psi.bulk.values, Q.bulk.values / sigma_t[None],
        )
