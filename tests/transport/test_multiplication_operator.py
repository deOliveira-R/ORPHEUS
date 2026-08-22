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
   :math:`\mathrm{spec}(M[f])=\mathrm{ess\,range}(f)` (``is_invertible``
   iff ``min|f|>0``), and the homomorphism :math:`M[f]M[g]=M[fg]` tested at
   the VALUES level via the numerics engine on the raw product array
   (``σ·σ'`` has units ``cm⁻²``, so the homomorphism is a numerics-engine
   property — built directly, NOT as a ``CrossSectionField``).

The resolvent-assembly leg (``L + M[σ] → StreamingCollisionOperator``, the
sweep/matvec matching the legacy resolvent and the analytical
``k_∞ = νΣf/Σa`` / streaming-equilibrium references) is verified by the
EXISTING structurally-independent gates that stay green:
``tests/sn/verification/analytical/test_kinf_homogeneous.py`` and
``tests/sn/operators/test_streaming_collision_operator.py`` (the
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
    DiagonalOperator,
    IdentityOperator,
    IncompatibleOperatorComposition,
    NotInvertible,
    OperatorSum,
    ZeroOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
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
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
    return replace(
        state,
        interior=replace(
            state.interior,
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
        engine_result = M.apply(psi).interior.values
        legacy_result = sigma[None] * psi.interior.values

        np.testing.assert_array_equal(engine_result, legacy_result)

    def test_2d_codomain_is_source_not_flux(self):
        """The codomain of ``M[σ]`` is a SOURCE (collision rate), not a flux."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        M = _multiplier(sn, _positive_sigma(sn, ng=2, seed=303))
        out = M.apply(_random_state(sn, ng=2, seed=404))
        _require(
            isinstance(out.interior, AngularSourceSink),
            f"M[σ].apply codomain bulk must be AngularSourceSink (a SOURCE), "
            f"got {type(out.interior).__name__}",
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
            isinstance(out.interior, AngularFlux),
            f"M[σ].solve codomain bulk must be AngularFlux (a flux), "
            f"got {type(out.interior).__name__}",
        )
        np.testing.assert_array_equal(out.interior.values, q.interior.values / sigma[None])

    def test_round_trip_solve_apply_identity(self):
        """``solve(apply(ψ)) == ψ`` for a coefficient bounded away from 0."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = _positive_sigma(sn, ng=2, seed=707)
        psi = _random_state(sn, ng=2, seed=808)
        M = _multiplier(sn, sigma)
        round_trip = M.solve(M.apply(psi))
        np.testing.assert_allclose(
            round_trip.interior.values, psi.interior.values, rtol=1e-14, atol=1e-15,
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
            M1.apply(psi).interior.values, I.apply(psi.interior.values),
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
            isinstance(out.interior, AngularSourceSink),
            f"M[0].apply codomain bulk must be AngularSourceSink (a zero "
            f"SOURCE, not a zero flux), got {type(out.interior).__name__}",
        )
        np.testing.assert_array_equal(out.interior.values, 0.0)

        # Structurally-independent: a codomain-aware ZeroOperator whose
        # codomain_zero builds the zero source-composite from psi.
        def _zero_source(p: TimedFullField) -> TimedFullField:
            return replace(
                p,
                interior=AngularSourceSink(
                    values=np.zeros_like(p.interior.values),
                    space=p.interior.space,
                ),
            )

        Z = ZeroOperator(codomain_zero=_zero_source)
        np.testing.assert_array_equal(
            out.interior.values, Z.apply(psi).interior.values,
        )
        # M[0] is NOT invertible — its spectrum is {0}.
        _require(
            not M0.is_invertible,
            "M[0] must NOT be invertible (spectrum is {0})",
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

        lhs = M_combo.apply(psi).interior.values
        # a·M[f]·ψ + b·M[g]·ψ at the values level.
        rhs = (
            a * M_f.apply(psi).interior.values
            + b * M_g.apply(psi).interior.values
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

        forward = M.apply(psi).interior.values
        np.testing.assert_array_equal(
            M.apply_transpose(psi).interior.values, forward,
        )
        # .H action (metric-blind, Euclidean) equals the forward action.
        np.testing.assert_array_equal(M.H.apply(psi).interior.values, forward)

    def test_spectrum_invertible_iff_min_abs_positive(self):
        r""":math:`\mathrm{spec}(M[f]) = \mathrm{ess\,range}(f)`.

        ``is_invertible`` iff ``min|f| > 0``. A coefficient with a zero
        entry has 0 in its spectrum → the inverse axis REVOKED (solve
        raises NotInvertible), while apply still works.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        psi = _random_state(sn, ng=2, seed=181)

        # Bounded away from 0 → invertible (and always adjointable).
        pos = _positive_sigma(sn, ng=2, seed=191)
        M_pos = _multiplier(sn, pos)
        _require(
            M_pos.is_invertible and M_pos.is_adjointable,
            "min|f|>0 must advertise both structural axes",
        )

        # A single zero entry → 0 in the spectrum → invertibility revoked.
        with_zero = pos.copy()
        with_zero[1, 2, 1] = 0.0
        M_zero = _multiplier(sn, with_zero)
        _require(
            M_zero.is_adjointable and not M_zero.is_invertible,
            "a coefficient with a zero entry must revoke invertibility "
            "(and keep the self-adjoint transpose)",
        )
        with pytest.raises(NotInvertible):
            M_zero.solve(psi)
        # apply still works (the multiplier is defined; only the inverse
        # is undefined at the zero).
        np.testing.assert_array_equal(
            M_zero.apply(psi).interior.values, with_zero[None] * psi.interior.values,
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
        composed = eng_f.apply(eng_g.apply(psi.interior.values))

        # M[f·g] — the engine on the raw product array (units cm⁻², so
        # NOT a CrossSectionField; the homomorphism is a numerics fact).
        eng_fg = DiagonalOperator(f * g, broadcast_axes=(0,))
        product = eng_fg.apply(psi.interior.values)

        np.testing.assert_allclose(composed, product, rtol=1e-14, atol=1e-15)


# ═══════════════════════════════════════════════════════════════════════
# 3. L0 streaming-equilibrium — the bare-engine values check.
#
# The full resolvent-assembly leg (L + M[σ] → StreamingCollisionOperator, sweep
# matching kinf / streaming-equilibrium) is the EXISTING gate set that
# stays green:
#   tests/sn/verification/analytical/test_kinf_homogeneous.py
#   tests/sn/operators/test_streaming_collision_operator.py::
#     TestStreamingCollisionSolveBridgeRegression::{test_si_carve_recovers_analytical_kinf,
#       test_fixed_source_homogeneous_reflective_recovers_q_over_sigma_composite}
# This pins the BARE multiplier algebra without a solve.
# ═══════════════════════════════════════════════════════════════════════


class TestStreamingEquilibriumValuesLeg:
    r"""On a flat-σ collision-only balance ``C ψ = Q``, ``ψ = Q/σ``.

    A degenerate streaming-equilibrium check at the multiplier level: with
    no streaming term the within-group balance is purely ``M[σ_t]ψ = Q``,
    whose solution is ``ψ = M[σ_t]^{-1}Q = Q/σ_t``. This is the L0 leg
    that the bare engine is correct; the full ``φ = Q/(Σ_t(1−c))``
    streaming-equilibrium with scattering is the StreamingCollisionOperator gate.
    """

    def test_collision_only_balance_recovers_q_over_sigma(self):
        sn = _slab_mesh(nx=4, ng=2)
        sigma_t = 0.8 * np.ones((2, *sn.spatial_shape))  # flat Σ_t.
        Q = _random_state(sn, ng=2, seed=241)
        M = _multiplier(sn, sigma_t)
        psi = M.solve(Q)
        # M[σ_t]·ψ == Q (the balance is satisfied).
        np.testing.assert_allclose(
            M.apply(psi).interior.values, Q.interior.values, rtol=1e-14, atol=1e-15,
        )
        # And ψ == Q/σ_t explicitly.
        np.testing.assert_array_equal(
            psi.interior.values, Q.interior.values / sigma_t[None],
        )


# ═══════════════════════════════════════════════════════════════════════
# 3. #261 — the composition-guard space metadata.
#
# The gap the collapse of ``CollisionOperator`` folded UP into the base:
# a multiplier may carry the optional composite ``FunctionSpace`` and so
# JOIN the W-D ``OperatorSum`` composition guard (``L + C``). The retired
# subclass added NOTHING the base lacked once the base gained ``space`` —
# these gates pin that the base now does the job (the collapse's premise).
# ═══════════════════════════════════════════════════════════════════════


class TestSpaceMetadataAndGuardJoin:
    """``space`` → ``domain``/``codomain`` → the multiplier joins the guard."""

    def test_space_anonymous_by_default(self):
        """No ``space`` → ``domain``/``codomain`` are ``None`` (the mesh-free
        default — the guard SKIPS a ``None``-spaced operand, the pre-W-D
        behaviour preserved for non-composite multipliers)."""
        sn = _slab_mesh()
        M = _multiplier(sn, _positive_sigma(sn))
        _require(M.domain is None, "default domain must be None (space-anonymous)")
        _require(M.codomain is None, "default codomain must be None")

    def test_space_threads_to_domain_and_codomain(self):
        """A supplied ``space`` is reported as BOTH ``domain`` AND ``codomain``
        — the multiplier is space-endomorphic (flux block → source block, same
        shape), so ``codomain == domain``."""
        sn = _slab_mesh()
        space = sn.full_field_space
        M = MultiplicationOperator(
            coefficient=CrossSectionField.from_mesh(_positive_sigma(sn), sn),
            space=space,
        )
        _require(M.domain is space, "domain must be the supplied space")
        _require(M.codomain is space, "codomain must equal domain (endomorphic)")

    def test_from_mesh_defaults_space_from_mesh(self):
        """``from_mesh(σ, sn_mesh)`` defaults ``space`` to the mesh's
        ``full_field_space`` — the faithful drop-in for the retired
        ``CollisionOperator(sn_mesh, σ)``, which reached the same composite
        space through ``sn_mesh.full_field_space`` (#261)."""
        sn = _slab_mesh()
        M = MultiplicationOperator.from_mesh(_positive_sigma(sn), sn)
        _require(
            M.space is sn.full_field_space,
            "from_mesh must default space to mesh.full_field_space",
        )

    def test_from_mesh_chain_short_circuits_on_full_field_space(self):
        """D7 (CS1, path gate; re-derived at CS4b S2) — on a method mesh the
        OPERATOR's space resolves to ``full_field_space`` by identity, never
        to the scalar bulk.

        The chain order is load-bearing: ``SNMesh.bulk_space`` is the
        scalar ``(ng, *spatial)`` bulk, NOT the angular composite, so a
        flipped chain would silently re-space every SN multiplier. Until
        CS4b this was pinned with a poison on ``bulk_space`` ("the chain
        NEVER consults it") — that instrument died when the space-source
        flip made ``bulk_space`` the shared arm of every carrier mint (the
        COEFFICIENT field's space is now correctly sourced from it, and
        the energy-arm single source routes ``angular_bulk_space`` through
        it too). The claim that survives is the OUTCOME, pinned by
        identity on both legs: the operator's space is the composite, and
        the coefficient's is the carrier's scalar bulk.
        """
        sn = _slab_mesh()
        sigma = _positive_sigma(sn)

        M = MultiplicationOperator.from_mesh(sigma, sn)
        _require(
            M.space is sn.full_field_space,
            "the chain must resolve the composite by identity, first arm",
        )
        _require(
            M.space is not sn.bulk_space,
            "the operator's space must never be the scalar bulk",
        )
        _require(
            M.coefficient.space is sn.bulk_space,
            "the coefficient's space is the carrier's cached scalar bulk "
            "(CS4b: one mint per carrier)",
        )

    def test_matching_spaces_compose(self):
        """Two multipliers on the SAME space compose into an ``OperatorSum``
        (the guard VALIDATES — equal domains AND codomains)."""
        sn = _slab_mesh()
        space = sn.full_field_space
        a = MultiplicationOperator(
            coefficient=CrossSectionField.from_mesh(_positive_sigma(sn, seed=1), sn),
            space=space,
        )
        b = MultiplicationOperator(
            coefficient=CrossSectionField.from_mesh(_positive_sigma(sn, seed=2), sn),
            space=space,
        )
        _require(
            isinstance(a + b, OperatorSum),
            "matching-space multipliers must compose into an OperatorSum",
        )

    def test_mismatched_spaces_red_the_guard(self):
        """Two multipliers on DIFFERENT spaces RED the ``OperatorSum``
        composition guard — the gap-fill's whole purpose. A formerly
        ``None``-spaced ``CollisionOperator`` skipped this check; the
        space-carrying ``MultiplicationOperator`` no longer can mis-compose
        silently."""
        sn_a = _slab_mesh(nx=4)
        sn_b = _slab_mesh(nx=6)  # distinct spatial shape → distinct full_field_space
        a = MultiplicationOperator(
            coefficient=CrossSectionField.from_mesh(_positive_sigma(sn_a), sn_a),
            space=sn_a.full_field_space,
        )
        b = MultiplicationOperator(
            coefficient=CrossSectionField.from_mesh(_positive_sigma(sn_b), sn_b),
            space=sn_b.full_field_space,
        )
        with pytest.raises(IncompatibleOperatorComposition):
            _ = a + b


# ═══════════════════════════════════════════════════════════════════════
# 4. The meshless bare-ndarray apply arm (#276 — homogeneous C − K_iso).
#
# The singledispatch arm that lets ``M[σ]`` act on a bare ``(ng, *spatial)``
# group block (NO ordinate axis, NO mesh), so the collision operator
# composes on a mesh-free MaterialMesh and realizes as a dense matrix via
# apply-to-basis-vectors. The load-bearing gate is CROSS-ARM CONSISTENCY:
# the bare arm must agree with the validated FullField-engine arm per
# ordinate (structurally independent — raw ``coeff*x`` vs the engine's
# leading-axis broadcast). The single-cell homogeneous config is
# axis-/group-blind, so these gates run on the discriminating 2-D nx≠ny
# carrier (a wrong axis or a group-drop is observable here, NOT meshless).
# ═══════════════════════════════════════════════════════════════════════


class TestMeshlessBareArm:
    """``M[σ].apply(bare (ng,*spatial))`` — the meshless diagonal multiply."""

    @staticmethod
    def _asym_sigma(ng: int, nx: int, ny: int) -> np.ndarray:
        """Deterministic σ asymmetric across groups AND the two unequal
        spatial axes — a group-drop or an x↔y swap is observable."""
        return (
            np.arange(1, 1 + ng * nx * ny, dtype=float).reshape(ng, nx, ny)
        ) * 0.1

    def test_bare_arm_is_sigma_times_x(self):
        r""":math:`M[\sigma]\,\phi = \sigma \odot \phi` on a bare ``(ng,*spatial)``.

        The explicit ``out.shape == sigma.shape`` is load-bearing:
        ``assert_array_equal`` BROADCASTS, so a spurious leading axis (the
        engine-style ``(1,ng,nx,ny)`` mutation) would pass the value compare
        — only the shape gate reds it.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = self._asym_sigma(2, 5, 3)
        x = np.random.default_rng(0).uniform(0.1, 1.0, size=(2, 5, 3))
        C = _multiplier(sn, sigma)

        out = C.apply(x)
        _require(
            out.shape == sigma.shape,
            f"bare apply must return (ng,*spatial)={sigma.shape}, got {out.shape}",
        )
        np.testing.assert_array_equal(out, sigma * x)

    def test_cross_arm_consistency_bare_equals_fullfield_per_ordinate(self):
        r"""The bare arm equals the FullField-engine arm on EVERY ordinate.

        The load-bearing structurally-independent pin: feed the SAME block
        ``x`` to the bare arm AND broadcast it across all N ordinates into a
        :class:`FullField`; the bare result must equal the FullField bulk at
        every ordinate. The two arms reach the value by DIFFERENT code paths
        (raw ``coeff*x`` vs the ``DiagonalOperator`` leading-axis broadcast),
        so a wrong axis / broadcast in the new arm reds against the validated
        206-caller engine path.
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        sigma = self._asym_sigma(2, 5, 3)
        x = np.random.default_rng(1).uniform(0.1, 1.0, size=(2, 5, 3))
        C = _multiplier(sn, sigma)

        N = sn.quad.N
        bulk_vals = np.broadcast_to(x[None], (N, 2, 5, 3)).copy()
        state = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        )
        psi_bcast = replace(state, interior=replace(state.interior, values=bulk_vals))

        bare_out = C.apply(x)
        ff_out = C.apply(psi_bcast).interior.values
        _require(
            bare_out.shape == (2, 5, 3),
            f"bare apply shape must be (ng,nx,ny)=(2,5,3), got {bare_out.shape}",
        )
        for n in range(N):
            np.testing.assert_array_equal(bare_out, ff_out[n])

    def test_bare_arm_is_self_adjoint(self):
        r""":math:`M[\sigma]^* = M[\sigma]` on a bare block (real coefficient)."""
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        C = _multiplier(sn, self._asym_sigma(2, 5, 3))
        x = np.random.default_rng(2).uniform(0.1, 1.0, size=(2, 5, 3))
        np.testing.assert_array_equal(C.apply_transpose(x), C.apply(x))

    def test_dispatch_raises_on_unsupported_type(self):
        """ndarray-only scope: an unregistered carrier raises ``TypeError``.

        Pins that the bare arm is ndarray-ONLY — a typed :class:`ScalarFlux`
        is NOT silently accepted (false symmetry with ``FissionOperator``'s
        ScalarFlux arm, deferred until a real scalar-flux collision consumer
        exists; ``coding-elegance`` Pattern 6).
        """
        sn = _cartesian_2d_mesh(nx=5, ny=3, ng=2)
        C = _multiplier(sn, self._asym_sigma(2, 5, 3))
        with pytest.raises(TypeError):
            C.apply(object())
        with pytest.raises(TypeError):
            C.apply(ScalarFlux.from_mesh(self._asym_sigma(2, 5, 3), sn))


class TestInverseOperatorFace:
    """#226 taxonomy step 1 (§12.1/§12.4): the typed ``.inverse()`` face."""

    @pytest.mark.foundation
    def test_inverse_apply_is_solve_bit_identical_and_involutive(self):
        """``M.inverse().apply ≡ M.solve`` bit-id (delegation, not a
        reciprocal-field twin) + I1 round-trip on the bulk values +
        ``(M⁻¹)⁻¹ is M`` by object identity."""
        mesh = _slab_mesh()
        sigma = _positive_sigma(mesh)
        M = _multiplier(mesh, sigma)
        psi = _random_state(mesh)
        q = M.apply(psi)  # a genuine source-typed rhs
        inv = M.inverse()
        np.testing.assert_array_equal(
            inv.apply(q).interior.values, M.solve(q).interior.values,
            err_msg="M.inverse().apply is not M.solve bit-identically",
        )
        np.testing.assert_array_almost_equal_nulp(  # I1: M⁻¹(M ψ) = ψ
            inv.apply(q).interior.values, np.asarray(psi.interior.values), nulp=2
        )
        assert inv.inverse() is M  # (M⁻¹)⁻¹ — object identity

    @pytest.mark.foundation
    def test_singular_coefficient_inverse_raises(self):
        """NEGATIVE (§12.4): a TRUE-zero entry (spectrum law) refuses."""
        mesh = _slab_mesh()
        sigma = _positive_sigma(mesh).copy()
        sigma[0, 0] = 0.0
        with pytest.raises(NotInvertible, match="zero"):
            _multiplier(mesh, sigma).inverse()
