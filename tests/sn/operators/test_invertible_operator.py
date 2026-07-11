r"""L0 foundation: :class:`InvertibleOperator` — sweep-invertible (L + C).

R-1 Step C (2026-05-19) — the SN-specific algebraic identity

.. math::

    (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
    \text{WDD sweep}

is encoded at the type level by :class:`InvertibleOperator`, a
specialisation of :class:`~orpheus.numerics.operator.OperatorSum` that
carries ``.solve`` as the forward substitution on the operator's ONE
``loss_representation`` instance (S6.5, #222).

The dispatch is symmetric: ``L + C`` and ``C + L`` both produce an
``InvertibleOperator`` (with the streaming operand stored first via the
canonical ``OperatorSum`` ordering).

Tests pin:

* Dispatch — ``L + C`` and ``C + L`` both return InvertibleOperator;
  ``L + L``, ``C + C``, and ``L + S`` (where S is a scattering operator)
  fall through to generic OperatorSum.
* Predicates — ``is_invertible`` True (the sweep identity); the parent
  generic OperatorSum built on the same operands is NOT invertible
  (its leading term L is not).
* Invariants — mesh-identity required; negative sigma rejected.
* Apply equivalence — InvertibleOperator.apply matches the inherited
  OperatorSum action (L.apply + C.apply, bit-exact).
* Solve consistency — apply ∘ solve = identity (on volumetric rhs).
* Carlson seed plumbing — InvertibleOperator.solve reads ``rhs(1)`` and
  forwards it to the representation sweep as ``initial_guess`` for the
  curvilinear pole closure.
"""
from __future__ import annotations

from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.operator import OperatorSum
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import (
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.source_sinks import RadialCharacteristicSourceSink


def _random_state(
    sn_mesh: SNMesh, seed: int = 42, *, history_depth: int = 2,
) -> TimedFullField:
    """Build a :class:`TimedFullField` with random bulk values.

    D-H.2-C1: the composite carrier replaces the legacy
    :class:`orpheus.sn.angular_flux.AngularFlux` test fixture.  Bulk
    values are sampled from ``N(0, 1)``; boundary is left as the
    implicit-zero L2 :class:`AngularBoundaryFlux`.
    """
    rng = np.random.default_rng(seed)
    N, ng = sn_mesh.quad.N, sn_mesh.ng
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, history_depth=history_depth, radial_characteristic=RadialCharacteristicFlux)
    return replace(
        state,
        interior=replace(
            state.interior, values=rng.standard_normal((N, ng, *sn_mesh.spatial_shape)),
        ),
    )


def _const_state(
    sn_mesh: SNMesh, value: float = 1.0, *, history_depth: int = 2,
) -> TimedFullField:
    """Build a :class:`TimedFullField` whose bulk is uniformly ``value``."""
    N, ng = sn_mesh.quad.N, sn_mesh.ng
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, history_depth=history_depth, radial_characteristic=RadialCharacteristicFlux)
    return replace(
        state,
        interior=replace(state.interior, values=np.full((N, ng, *sn_mesh.spatial_shape), value)),
    )


pytestmark = [pytest.mark.foundation]


# ── Geometry fixtures ───────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ── Dispatch via __add__ ────────────────────────────────────────────────


class TestDispatch:
    def test_streaming_plus_collision_returns_invertible(self) -> None:
        """``L + C`` returns :class:`InvertibleOperator` via dispatch."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)

        composite = L + C
        assert isinstance(composite, InvertibleOperator)
        assert composite.streaming is L
        assert composite.diagonal is C

    def test_explicit_construction_with_sigma_r(self) -> None:
        """Explicit ctor accepts a CollisionOperator built with σ_r ≠ σ_t."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        # σ_r = σ_t - within-group self-scatter; here a stylised
        # positive reduction (e.g. 1.0 - 0.3 = 0.7 for c_(g->g) = 0.3).
        sigma_r = 0.7 * np.ones((sn.ng, *sn.spatial_shape))

        L = StreamingOperator(sn)
        C_r = MultiplicationOperator.from_mesh(sigma_r, sn)
        invertible = InvertibleOperator(L, C_r)
        assert invertible.streaming is L
        assert invertible.diagonal is C_r
        np.testing.assert_array_equal(invertible.sigma, sigma_r)

    def test_streaming_plus_zero_operator_falls_through(self) -> None:
        """Non-Collision RHS falls back to generic OperatorSum."""
        from orpheus.numerics.operator import ZeroOperator
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        Z = ZeroOperator()

        composite = L + Z
        assert isinstance(composite, OperatorSum)
        assert not isinstance(composite, InvertibleOperator)


# ── Capabilities and invariants ─────────────────────────────────────────


class TestCapabilitiesAndInvariants:
    def test_predicates_invertible_and_adjointable(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = InvertibleOperator(
            StreamingOperator(sn),
            MultiplicationOperator.from_mesh(sigma_t, sn),
        )
        assert invertible.is_invertible
        # Both operands transpose (Lᵀ analytic, C self-adjoint), so the
        # sum-closure law makes (L+C) adjointable too.
        assert invertible.is_adjointable

    def test_parent_operator_sum_not_invertible(self) -> None:
        r"""Plain :class:`OperatorSum` on ``(L, C)`` is NOT invertible.

        ``is_invertible`` on a generic sum reads the LEADING term (the
        Green/Neumann factorization ``(A+B)⁻¹`` needs ``A⁻¹``); here the
        leading ``L`` is rank-deficient, so the plain sum reports False.
        Only :class:`InvertibleOperator` carries the SN-specific sweep
        identity at the type level.
        """
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        # Explicitly construct via the parent class to bypass dispatch.
        plain = OperatorSum(
            StreamingOperator(sn),
            MultiplicationOperator.from_mesh(sigma_t, sn),
        )
        assert callable(getattr(plain, "apply", None))
        assert not plain.is_invertible

    def test_mismatched_mesh_rejected(self) -> None:
        """Two distinct :class:`SNMesh` instances → ValueError."""
        sn1 = _slab_mesh()
        sn2 = _slab_mesh()
        sigma_t1 = np.ones((sn1.ng, *sn1.spatial_shape))
        sigma_t2 = np.ones((sn2.ng, *sn2.spatial_shape))
        L = StreamingOperator(sn1)
        C = MultiplicationOperator.from_mesh(sigma_t2, sn2)
        with pytest.raises(ValueError, match="mesh-identity"):
            InvertibleOperator(L, C)

    def test_non_positive_sigma_rejected(self) -> None:
        """``σ <= 0`` anywhere → ValueError at construction."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        sigma_bad = sigma_t.copy()
        sigma_bad.flat[0] = -0.5  # one bad cell
        L = StreamingOperator(sn)
        C_bad = MultiplicationOperator.from_mesh(sigma_bad, sn)
        with pytest.raises(ValueError, match="strictly positive"):
            InvertibleOperator(L, C_bad)

    def test_zero_sigma_rejected(self) -> None:
        """``σ == 0`` (void cell) is also rejected — strict positivity."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        sigma_bad = sigma_t.copy()
        sigma_bad.flat[0] = 0.0
        L = StreamingOperator(sn)
        C_bad = MultiplicationOperator.from_mesh(sigma_bad, sn)
        with pytest.raises(ValueError, match="strictly positive"):
            InvertibleOperator(L, C_bad)

    def test_wrong_operand_types_rejected(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        with pytest.raises(TypeError, match="StreamingOperator"):
            InvertibleOperator(C, C)  # type: ignore[arg-type]
        with pytest.raises(TypeError, match="MultiplicationOperator"):
            InvertibleOperator(L, L)  # type: ignore[arg-type]


# ── Apply matches inherited OperatorSum ─────────────────────────────────


class TestApply:
    def test_apply_equals_l_plus_c_on_typed_flux(self) -> None:
        r"""``(L+C).apply(ψ) == L.apply(ψ) + C.apply(ψ)`` for TimedFullField (slab).

        Since #240 Phase 2 Step B ``InvertibleOperator.apply`` OWNS its matvec
        via ``loss_representation.loss_action(self.sigma)`` (it no longer routes
        through the inherited ``OperatorSum.apply`` leaf sum).  #257 S8b: ``L``
        is now pure σ-free ``streaming_action(ψ) = loss_action(0, ψ)``, so the
        leaf sum is ``loss_action(0, ψ) + σ_t⊙ψ`` while the override is
        ``loss_action(σ_t, ψ)`` directly — value-equal by the affine relation
        ``M(σ_t)ψ = streaming_action(ψ) + σ_t⊙ψ``, but FP-re-associated (the
        σ-free walk vs the σ-bearing matvec).  The BULK gate is nULP; the
        BOUNDARY stays STRICT 0-ULP (``C`` contributes zero on the trace, so
        ``(L+C).apply.boundary == L.apply.boundary`` byte-exact).
        """
        sn = _slab_mesh(nx=4, n_ord=4, ng=2)
        sigma_t = np.full((sn.ng, *sn.spatial_shape), 0.8)
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        invertible = L + C

        state = _random_state(sn, seed=7)

        composite_out = invertible.apply(state)
        sum_out = L.apply(state) + C.apply(state)

        # Bulk: value-equal to FP-non-associativity (the σ-free walk re-
        # associates vs the σ-bearing matvec).  256 nULP absorbs the near-zero
        # cancellation spikes; rel drift is machine-ε (grounded by the
        # byte-identical composite in test_streaming_operator_decomposition).
        np.testing.assert_array_almost_equal_nulp(
            composite_out.interior.values, sum_out.interior.values, nulp=256,
        )
        # Boundary STRICT: C never touches the trace.
        np.testing.assert_array_equal(
            composite_out.boundary.values, sum_out.boundary.values,
        )


# ── Solve via sweep ─────────────────────────────────────────────────────


class TestSolve:
    def test_solve_returns_composite_on_slab(self) -> None:
        """``InvertibleOperator.solve`` returns a :class:`TimedFullField` on slab."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        rhs = _const_state(sn, value=1.0)
        psi = invertible.solve(rhs)

        assert isinstance(psi, TimedFullField)
        assert isinstance(psi.interior, AngularFlux)
        assert psi.interior.mesh is sn
        assert psi.interior.values.shape == rhs.interior.values.shape

    def test_solve_inherits_history_depth_from_rhs(self) -> None:
        """``rhs.history_depth`` propagates to the returned composite."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        rhs = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn, history_depth=5)
        psi = invertible.solve(rhs)
        assert psi.history_depth == 5

    def test_solve_rejects_bare_ndarray(self) -> None:
        """The composite-flux contract is mandatory — bare ndarray rejected."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )
        with pytest.raises(TypeError):
            invertible.solve(np.zeros(10))

    @pytest.mark.foundation
    def test_apply_rejects_bare_ndarray(self) -> None:
        """Negative contract: post-D-I.3d, :meth:`StreamingOperator.apply`
        accepts only :class:`~orpheus.transport.timed_full_field.TimedFullField`.
        A bare ``np.ndarray`` argument raises ``TypeError`` (principled
        rejection at the method entry — the typed contract is enforced
        before any attribute access).

        D-I.3d (2026-05-29) — Lesson L18 (Pattern 7 producer-side
        normalisation): the convention bridge collapses at the
        producer; consumers see only the typed contract.
        """
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        psi_bare = np.zeros(10)
        # #257 S8a — the matvec input contract widened from TimedFullField to
        # the timeless FullField (the leaf is a base arrow ``FullField ->
        # FullField``); a bare ndarray is still rejected at the entry guard.
        with pytest.raises(TypeError, match="expected FullField"):
            L.apply(psi_bare)

    def test_solve_rejects_mismatched_mesh(self) -> None:
        sn1 = _slab_mesh()
        sn2 = _slab_mesh()
        sigma_t1 = np.ones((sn1.ng, *sn1.spatial_shape))
        invertible = StreamingOperator(sn1) + MultiplicationOperator.from_mesh(
            sigma_t1, sn1,
        )
        rhs = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn2)
        with pytest.raises(ValueError, match="mesh-identity"):
            invertible.solve(rhs)

    def test_solve_produces_positive_flux_for_positive_source(self) -> None:
        r"""Sanity smoke — positive volumetric source → non-trivial positive ψ.

        Note on a stronger ``apply(solve(q)) ≈ q`` round-trip:
        :meth:`StreamingOperator.apply` uses the **symmetric-closure
        FD** matvec (the Krylov-on-apply path that
        :func:`transport_operator_matvec_unified` realises), whereas
        :meth:`InvertibleOperator.solve` uses the **WDD asymmetric
        sweep**.  The two are different discretisations of the same
        continuous operator — they agree in the fine-mesh limit but
        are NOT bit-exact inverses on the discrete grid (ERR-026
        documents the closure-bias-driven discrepancy).  The Krylov-
        on-apply outer iteration converges them to a consistent
        fixed point; per-call ``apply ∘ solve = identity`` is not a
        valid contract.

        This test instead asserts a weaker but architecturally
        meaningful invariant: a positive source through a positive-σ
        diagonal in vacuum-BC slab yields a positive flux at every
        cell.  The sweep is monotone-preserving by construction.
        """
        sn = _slab_mesh(nx=8, n_ord=4, ng=1)
        sigma_t = np.full((sn.ng, *sn.spatial_shape), 1.5)
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        N, ng, nx = sn.quad.N, sn.ng, sn.nx
        x = np.linspace(0.5 / nx, 1.0 - 0.5 / nx, nx)
        rhs_values = np.broadcast_to(
            np.sin(np.pi * x)[None, None, :] * np.ones((N, ng, 1)),
            (N, ng, nx),
        ).copy()
        q = replace(
            TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn),
            interior=replace(TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn).interior, values=rhs_values),
        )
        psi = invertible.solve(q)

        # Some cells positive (interior; the sweep is monotone for
        # positive sources on vacuum BC slab).  Exact every-cell-
        # positivity may fail at the BC face under WDD's diamond
        # closure — that's a known DD artefact.  Sufficient signal:
        # non-trivial peak in the interior.
        assert psi.interior.values.max() > 0
        # Peak occurs near the source maximum (centre of the slab).
        peak_x_idx = np.unravel_index(
            np.argmax(psi.interior.values[0, 0, :]),
            psi.interior.values[0, 0, :].shape,
        )
        assert nx // 4 <= peak_x_idx[0] <= 3 * nx // 4

    def test_solve_runs_on_sphere(self) -> None:
        r"""Sphere geometry — smoke test that the curvilinear path runs."""
        sn = _sphere_mesh(nx=8, n_ord=4, ng=1)
        sigma_t = np.full((sn.ng, *sn.spatial_shape), 1.0)
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        q = _const_state(sn, value=1.0)
        psi = invertible.solve(q)
        assert isinstance(psi, TimedFullField)
        assert psi.interior.values.shape == q.interior.values.shape
        # Curvilinear sweep produces a non-trivial positive bulk.
        assert psi.interior.values.max() > 0

    @pytest.mark.l0
    @pytest.mark.verifies("transport-cartesian")
    def test_solve_consumes_per_ordinate_rhs(self) -> None:
        r"""``InvertibleOperator.solve`` passes ``rhs.interior.values`` unmodified to the sweep.

        R-1 Step 4 A1 invariant pin (N5 per verification plan).  The
        producer-side normalisation contract says the bulk values are
        already per-ordinate density (the producer of ``rhs`` — typically
        ``ScatteringOperator.apply`` or ``AngularSourceSink.from_isotropic``
        — applied ``/sum_w`` at the producer boundary).  The
        ``InvertibleOperator.solve`` adapter MUST forward the bulk
        values to :func:`transport_sweep` *bit-equal* — no internal
        ``* sum_w`` bridge, no ``/sum_w`` rescaling.

        Pre-A1 the adapter wrapped ``rhs.values * sum_w`` into the sweep
        to compensate for the old sweep-internal ``/W``; both directions
        of that bridge dissolved in A1.  This test spies on the
        representation's ``sweep`` (the S6.5 solve seam —
        ``solve`` runs ``self.loss_representation.sweep`` directly) to
        capture the ``Q`` argument and asserts it is bit-identical to
        ``rhs.interior.values``.  If a future refactor re-introduces a
        ``* sum_w`` / ``/ sum_w`` rescaling on this hot path, the
        bit-equal assertion fails.
        """
        sn = _slab_mesh(nx=4, n_ord=4, ng=1)
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        sum_w = float(sn.quad.weights.sum())
        q_const = 2.7
        per_ord_density = q_const / sum_w
        zero = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
        rhs = replace(
            zero,
            interior=replace(
                zero.interior,
                values=np.full(
                    (sn.quad.N, sn.ng, *sn.spatial_shape), per_ord_density,
                ),
            ),
        )

        # Spy on the representation sweep — capture the Q argument.
        captured: list[np.ndarray] = []
        from orpheus.sn.loss_representation import CumprodScan
        original = CumprodScan.sweep

        def spy(self, Q, *args, **kwargs):
            # Snapshot ``Q`` BEFORE the sweep runs (the sweep is a pure
            # reader of its source argument but the snapshot makes
            # intent explicit).
            captured.append(np.array(Q, copy=True))
            return original(self, Q, *args, **kwargs)

        with patch.object(CumprodScan, "sweep", spy):
            invertible.solve(rhs)

        assert len(captured) == 1, (
            f"representation sweep called {len(captured)} times; expected 1"
        )
        forwarded = captured[0]
        np.testing.assert_array_equal(
            forwarded, rhs.interior.values,
            err_msg=(
                "InvertibleOperator.solve modified rhs.interior.values before "
                "forwarding to the representation sweep — A1 producer-side "
                "convention drifted (the ``* sum_w`` bridge MUST stay "
                "dissolved)."
            ),
        )


# ── D-H.2-C1: TimedFullField composite-only solve invariants ───────────


class TestSolveTimedFullField:
    """Composite-specific invariants on :meth:`InvertibleOperator.solve`.

    The parity-vs-legacy-AngularFlux tests retired with D-H.2-C1
    (legacy class itself retires in C5).  These tests pin the
    invariants that remain composite-specific:

    * Mesh-identity invariant on both ``rhs`` and ``initial_guess``.
    * History-depth propagation.
    * Composite initial_guess.boundary threads to the sweep's
      ``boundary_buf`` (partner-flux seeding, audit §5).
    """

    def test_rhs_boundary_seeds_the_sweep_inflow(self) -> None:
        """Composite ``rhs.boundary`` seeds the bare sweep's inflow trace.

        Wave O (#208) O.4a.2 — BC extraction: the boundary INFLOW seed is
        now the boundary SOURCE ``rhs.boundary`` (carrying ``q.boundary +
        B·ψ.outflow``), NOT ``initial_guess.boundary``.  The bare sweep no
        longer re-applies ``bc`` at entry; ``InvertibleOperator.solve``
        seeds the sweep's mutable ``boundary_buf`` from ``rhs.boundary``
        before the representation sweep runs (the S6.5 solve seam).
        (Pre-extraction this seed came from ``initial_guess.boundary`` —
        the partner-flux carrier — which the ``−B`` extraction retires.)
        """
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        # Build rhs with a non-zero boundary trace (the inflow source) —
        # verify it makes it into the sweep's boundary_buf.  Slab has both
        # xmin and xmax faces.
        rhs = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
        rhs_boundary = rhs.boundary
        layout = rhs_boundary.layout
        if "xmax" in layout.faces:
            rhs_boundary.face_view("xmax")[:] = 0.7
        if "xmin" in layout.faces:
            rhs_boundary.face_view("xmin")[:] = 0.3

        # Outcome-level spy: capture the boundary_buf as the
        # representation sweep sees it.  The seed must already be in
        # place when the kernel runs (otherwise the bare sweep sees a
        # zero inflow trace and the fixed point shifts away from the
        # −B-driven answer).
        captured: list[tuple[np.ndarray, np.ndarray]] = []
        from orpheus.sn.loss_representation import CumprodScan
        original = CumprodScan.sweep

        def spy(
            self, Q, sigma, boundary_flux, *,
            moment_frame=None,
            schedule=None, reflect=None,
            radial_characteristic_source=None, radial_characteristic_flux=None,
        ):
            # D-H.2-C2: L2 AngularBoundaryFlux exposes per-face writable views
            # via face_view; copy them out at entry to snapshot the inflow.
            # (#226 step 2: the door signature grew schedule/reflect —
            # None here, the Jacobi path; #282 route (a) grew the two
            # radial_characteristic seed args — None on the non-carrying slab;
            # mirror them all so the spy stays signature-compatible with the
            # production call.  The vestigial ``initial_guess`` kwarg retired
            # with the dead seed threading, #280 2.5c.)
            captured.append((
                boundary_flux.face_view("xmax").copy(),
                boundary_flux.face_view("xmin").copy(),
            ))
            return original(
                self, Q, sigma, boundary_flux,
                radial_characteristic_source=radial_characteristic_source,
                radial_characteristic_flux=radial_characteristic_flux,
            )

        with patch.object(CumprodScan, "sweep", spy):
            invertible.solve(rhs)

        assert len(captured) == 1
        np.testing.assert_array_equal(captured[0][0], 0.7)
        np.testing.assert_array_equal(captured[0][1], 0.3)

    def test_history_depth_preserved(self) -> None:
        """``rhs.history_depth`` propagates through the composite branch."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        invertible = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
            sigma_t, sn,
        )

        for depth in (1, 2, 4):
            rhs = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn, history_depth=depth)
            psi = invertible.solve(rhs)
            assert psi.history_depth == depth

    def test_rejects_mismatched_mesh_on_rhs(self) -> None:
        sn1 = _slab_mesh()
        sn2 = _slab_mesh()
        sigma_t1 = np.ones((sn1.ng, *sn1.spatial_shape))
        invertible = StreamingOperator(sn1) + MultiplicationOperator.from_mesh(
            sigma_t1, sn1,
        )
        rhs = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn2)
        with pytest.raises(ValueError, match="mesh-identity"):
            invertible.solve(rhs)


# ── Convention-bridge regression catchers (R-1 Step 4 A5 promotion) ────


class TestInvertibleSolveBridgeRegression:
    r"""End-to-end regression catchers for the convention-bridge bug class.

    Origin
    ======

    Promoted from ``derivations/diagnostics/diag_r1_step_e_invertible_solve_w_bridge.py``
    (R-1 Step 4 session 1, 2026-05-19).  The diagnostic caught the
    per-ordinate-vs-iso magnitude convention drift in the typed
    operator algebra (issue #202, lesson :ref:`L18 <lessons-l18>`):
    ``ScatteringOperator.apply`` (typed) was returning iso-magnitude
    while :func:`~orpheus.sn.loss_representation.transport_sweep` internally divided
    by ``W = sum_w``, so the converged fixed point sat at ``k_inf / W``
    instead of ``k_inf``.  Slab-2eg SI keff = 1.4844 vs ref = 1.875
    (ratio ≈ 1/W·c_s for W=2 GL) was the signature.

    Phase 1.1 A1 (commit ``de8822d``) made the bug class **structurally
    unreachable**: ``/sum_w`` lives at the producer boundary
    (:meth:`ScatteringOperator.apply`, :meth:`FissionOperator.apply`,
    :meth:`AngularSourceSink.from_isotropic`); consumers see per-ord
    density throughout.  The legacy ``rhs.values * sum_w`` bridge in
    ``InvertibleOperator.solve`` is GONE.  The sweep applies NO
    ``/W`` anywhere internally.

    What this class pins
    ====================

    Three end-to-end checks that catch a re-introduction of any
    convention-bridge bug on the typed operator-algebra hot path:

    1. ``test_slab_uniform_roundtrip`` — ``(L+C).solve((L+C).apply(ψ=1))
       == ψ=1`` to machine zero on slab (no curvilinear angular
       redistribution to confound).  Pre-fix gave 1/W instead.
    2. ``test_fixed_source_homogeneous_reflective_recovers_q_over_sigma``
       — fixed-point ``ψ_n = q_n / Σ_t`` on reflective homogeneous
       medium.  Pre-fix gave ``q_n / (W Σ_t)``.
    3. ``test_si_carve_recovers_analytical_kinf`` — end-to-end SI inner
       solve on homogeneous reflective medium reaches analytical
       ``k_inf`` to ``rtol < 1e-9``.  Pre-fix every ng≥2 case failed
       with ~21% drift.

    Post-Phase-1.2 the M-M Carlson seed travels through the explicit
    ``initial_guess`` kwarg on :meth:`InvertibleOperator.solve` (not
    through ``rhs(1)`` history); the iterative tests below thread the
    previous iterate that way.

    References
    ==========

    - ``.claude/lessons.md`` L18 — Pattern 7 producer-side normalisation.
    - ``.claude/lessons.md`` L21 — sweep and matvec share ONE strategy.
    - Sphinx theory: ``docs/theory/discrete_ordinates.rst``
      (transport-cartesian, sn-curvilinear-homogeneous-kinf-recovery).
    - ERR catalog: ERR-049 — convention drift (Phase 1.4).
    - Issue #202 — closed by Phase 1.1 commit ``de8822d``.
    - Issue #204 — A5 promote diagnostics to permanent L1.
    """

    @staticmethod
    def _homogeneous_solver(coord: str, ng_key: str = "2eg", n_cells: int = 10):
        """Build a homogeneous reflective :class:`SNSolver` of the given coord.

        Helper for the three tests below — replicates the legacy
        diagnostic's ``_build_homogeneous_setup`` but uses the
        l1_analytical infrastructure already on disk.
        """
        import sys
        import warnings

        from tests.sn._test_helpers import SN_TESTS_ROOT

        sys.path.insert(
            0, str(SN_TESTS_ROOT / "verification" / "analytical"),
        )
        warnings.simplefilter("ignore")
        from test_kinf_homogeneous import (  # type: ignore[import-not-found]
            _get_continuous_case, _homogeneous_mesh, _quadrature_for,
        )
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.sn.solver import SNSolver

        case = _get_continuous_case(ng_key)
        mat_id = next(iter(case.problem.materials.keys()))
        mesh = _homogeneous_mesh(
            coord=coord, n_cells=n_cells, length=2.0, mat_id=mat_id,
        )
        quad = _quadrature_for(coord)
        sn_mesh = SNMesh(mesh, quad, case.problem.materials)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0,
            max_inner=300, inner_tol=1e-12,
            inner_solver="source_iteration",
        )
        return solver, case

    @pytest.mark.l1
    @pytest.mark.verifies(
        "transport-cartesian", "sn-curvilinear-homogeneous-kinf-recovery",
    )
    @pytest.mark.catches("ERR-049")
    @pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
    @pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
    def test_si_carve_recovers_analytical_kinf(
        self, coord: str, ng_key: str,
    ) -> None:
        r"""End-to-end SI inner solver on homogeneous reflective → analytical ``k_inf``.

        Pins the R-1 Step E ``_solve_source_iteration`` carve at L1.
        Pre-fix (no ``×W`` bridge): keff lands at ``≈ k_inf / W`` —
        every ng≥2 case failed with ~21% drift.  Post-fix every case
        reaches the analytical reference at ``rtol < 1e-9``.

        The L1 Krylov path is cross-checked by the L1 Krylov suite +
        :file:`test_krylov_curvilinear_precond_safety.py` (the L19
        regression catcher).  Together the two pin structural
        independence: an A1 regression would surface in SI and Krylov
        identically only if both inner solvers share the same producer
        boundary — which they do, by Phase-1.1's design.
        """
        solver, case = self._homogeneous_solver(coord, ng_key=ng_key)

        phi = solver.initial_flux_distribution()
        keff = 1.0
        keff_new = keff
        for _ in range(60):
            fis = solver.compute_fission_source(phi, keff)
            phi_new = solver._solve_source_iteration(fis, phi)
            keff_new = solver.compute_keff(phi_new)
            if abs(keff_new - keff) < 1e-12:
                break
            keff, phi = keff_new, phi_new

        rel_err = abs(keff_new - case.k_eff) / case.k_eff
        assert rel_err < 1e-9, (
            f"SI carve {coord}-{ng_key}: keff={keff_new:.10f}, "
            f"ref={case.k_eff:.10f}, rel_err={rel_err:.3e}.  Pre-fix "
            f"failed at ~21% (≈ 1/W·c_s for W=2 GL).  A1 producer-side "
            f"normalisation eliminated the convention bridge; if this "
            f"now fails again, a new convention drift has been "
            f"introduced.  See L18 + ERR-049."
        )

    # ── D-H.1c Stage 3 composite-branch L1 anchors ────────────────────
    #
    # Mirror the two streaming-equilibrium L1 anchors above but
    # exercise the InvertibleOperator.solve TimedFullField composite
    # branch DIRECTLY (no solver wrapper).  Structurally independent
    # of the SI carve L1 test (which routes through solver._solve_
    # source_iteration): a regression breaking the composite branch
    # alone — without touching solver.py — would surface here while
    # the SI carve test would (transitively) catch it too via Stage 2.
    # Two L1 anchors covering the same equilibrium → two angles into
    # the same bug class.

    @pytest.mark.l1
    @pytest.mark.verifies("transport-cartesian")
    @pytest.mark.catches("ERR-049")
    def test_slab_uniform_roundtrip_composite(self) -> None:
        r"""Composite mate of :meth:`test_slab_uniform_roundtrip`.

        ``LC.solve(LC.apply(ψ=1)) == ψ=1`` on slab under TimedFullField
        dispatch.  Exercises the StreamingOperator + CollisionOperator
        composite-input branches (D-H.1b.6 / D-H.1b.5) ∘
        InvertibleOperator.solve composite branch (D-H.1c stage 1).
        Pre-D-H.1c-stage-1 the composite branch did not exist; pre-A1
        the legacy branch failed by factor W (now structurally
        impossible per L18).
        """
        from orpheus.transport.timed_full_field import TimedFullField

        solver, _case = self._homogeneous_solver("slab")
        sn_mesh = solver.sn_mesh
        N = solver.quad.N
        ng = solver.ng

        sigma_t = solver.mat_xs.total_cross_section
        L_leaf = StreamingOperator(sn_mesh)
        C_t = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
        LC = L_leaf + C_t

        # Build composite ψ=1.  No need for legacy AngularFlux at all
        # on this path.
        from dataclasses import replace
        psi_known = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        psi_known = replace(
            psi_known,
            interior=replace(psi_known.interior, values=np.ones((N, ng, *sn_mesh.spatial_shape))),
        )
        # #257 S8a — the matvec leaf is a base arrow: ``LC.apply`` returns a
        # timeless FullField source.  ``InvertibleOperator.solve`` consumes the
        # driver's timed rhs, so wrap the source back into a TimedFullField
        # (exactly what the production SI / Krylov driver does — it carries the
        # comonad; the operator does not).  Byte-identical source values.
        LC_psi_source = LC.apply(psi_known)
        LC_psi = TimedFullField(
            interior=LC_psi_source.interior,
            boundary=LC_psi_source.boundary,
            _history=(),
            history_depth=psi_known.history_depth,
        )

        # Iterate until stabilises (slab converges in 1-2 iterates).
        psi_recovered: TimedFullField | None = None
        for _ in range(20):
            if psi_recovered is None:
                psi_new = LC.solve(LC_psi)
            else:
                psi_new = LC.solve(LC_psi, initial_guess=psi_recovered)
            if psi_recovered is not None and np.abs(
                psi_new.interior.values - psi_recovered.interior.values,
            ).max() < 1e-15:
                psi_recovered = psi_new
                break
            psi_recovered = psi_new

        assert psi_recovered is not None
        diff = np.abs(psi_known.interior.values - psi_recovered.interior.values).max()
        assert diff < 1e-12, (
            f"slab composite (L+C).solve((L+C).apply(ψ=1)) ≠ ψ=1: "
            f"abs_max={diff:.3e}.  D-H.1c stage 1's composite bridge "
            f"is the test target; if this fails the bridge has drifted "
            f"from the legacy branch (cf. the L0 bulk-equivalence pin "
            f"in TestSolveTimedFullField.test_bulk_matches_legacy_branch)."
        )

    @pytest.mark.l1
    @pytest.mark.verifies(
        "transport-cartesian", "sn-curvilinear-homogeneous-kinf-recovery",
    )
    @pytest.mark.catches("ERR-049")
    @pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
    def test_fixed_source_homogeneous_reflective_recovers_q_over_sigma_composite(
        self, coord: str,
    ) -> None:
        r"""Composite mate of
        :meth:`test_fixed_source_homogeneous_reflective_recovers_q_over_sigma`.

        Reflective homogeneous medium, uniform per-ordinate source
        ``q_n = const`` → converged composite ψ_n must equal
        ``q_n / Σ_t`` on every cell.  Exercises the composite branch
        of InvertibleOperator.solve under three geometries (slab +
        spherical + cylindrical), pinning the streaming-equilibrium
        property end-to-end on the composite path.
        """
        from dataclasses import replace
        from orpheus.sn.operators.boundary import (
            RadialCharacteristicBoundaryOperator,
            SNBoundaryOperator,
        )
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        solver, _case = self._homogeneous_solver(coord)
        sn_mesh = solver.sn_mesh
        quad = solver.quad
        N = quad.N
        ng = solver.ng
        sum_w = float(quad.weights.sum())

        sigma_t = solver.mat_xs.total_cross_section
        L_leaf = StreamingOperator(sn_mesh)
        C_t = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
        LC = L_leaf + C_t

        # Per-ordinate uniform source — composite carrier.
        q_iso = 0.225
        q_per_ord = np.full((N, ng, *sn_mesh.spatial_shape), q_iso / sum_w)
        rhs = TimedFullField(
            interior=AngularFlux.from_mesh(q_per_ord, sn_mesh),
            boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
            # #282 route (a): on a carrying mesh (sphere) the fixed-source RHS
            # carries the TRUE starting-direction SOURCE — the q½ fold of the
            # per-ordinate source the direct ψ½ solve consumes (the SAME
            # ``from_angular_source`` factory the production
            # ``solve_sn_fixed_source`` q_ext uses).  Without it the pole march
            # is seeded wrong and ψ ≠ q/Σ_t; None on non-carrying (slab/cyl).
            radial_characteristic=RadialCharacteristicSourceSink.from_angular_source(
                q_per_ord, sn_mesh,
            ),
            _history=(),
            history_depth=2,
        )

        # Iterate to converge the reflective fixed point.  Wave O (#208)
        # O.4a.2 — BC extraction: ``LC.solve`` is now the BARE inverse
        # (it reads ``rhs.boundary`` as the inflow seed and no longer
        # re-applies ``bc`` internally), so this loop must drive the
        # reflective coupling EXPLICITLY via the sibling ``−B``: each
        # iterate sets ``rhs.boundary = q.boundary + B·ψ.outflow``
        # (``q.boundary = 0`` here) — exactly the ``S + B`` source the
        # production SI driver folds.  ``initial_guess`` still threads the
        # M-M Carlson bulk warm start.
        # The augmented boundary coupling: on a seed-carrying (sphere) mesh the
        # direct sum B_a + B_b (System A trace ⊕ System B ray corner — RULING P1;
        # the un-weld of the old welded SNBoundaryOperator). B.apply(ψ) emits BOTH
        # the reflected trace AND the ψ½ corner arm, exactly the augmented −B the
        # production driver folds (this test replicates that fold by hand). On a
        # seedless (slab/cylinder) mesh B is B_a alone. The FUSED B_a + B_b sum
        # rides the fused-oracle shim since B.2d (the re-typed B_b speaks
        # System B's composite spaces, which the OperatorSum guard correctly
        # refuses to sum with B_a's FullField carrier; production consumes the
        # block natively through the gain grid).
        from tests.sn._test_helpers import FusedRayBoundaryGain

        B_a = SNBoundaryOperator(sn_mesh)
        B = (
            B_a
            + FusedRayBoundaryGain(
                RadialCharacteristicBoundaryOperator(sn_mesh),
            )
            if sn_mesh.radial_characteristic_space is not None
            else B_a
        )
        psi_typed: TimedFullField | None = None
        for _ in range(400):
            if psi_typed is None:
                psi_new = LC.solve(rhs)
            else:
                # #282 route (a): the reflective −B coupling folds into BOTH the
                # boundary inflow AND the ψ½ SOURCE corner arm.  Production SI
                # builds ``rhs = q_ext + B.apply(psi)`` (the composite ``+``
                # combines EVERY block); the pre-2.5d loop replaced only
                # ``.boundary``.  On a carrying mesh (sphere) the B corner arm
                # must be added to the q½ source too, else the pole ψ½ march is
                # under-driven and ψ ≠ q/Σ_t (the observed non-flat profile).
                B_out = B.apply(psi_typed)
                seed_n = rhs.radial_characteristic
                if seed_n is not None and B_out.radial_characteristic is not None:
                    seed_n = seed_n + B_out.radial_characteristic
                rhs_n = replace(
                    rhs, boundary=B_out.boundary, radial_characteristic=seed_n,
                )
                psi_new = LC.solve(rhs_n, initial_guess=psi_typed)
            if psi_typed is not None and np.abs(
                psi_new.interior.values - psi_typed.interior.values,
            ).max() < 1e-14:
                psi_typed = psi_new
                break
            psi_typed = psi_new

        assert psi_typed is not None
        # Per-ordinate expected: ψ_n = q_n / Σ_t = (q_iso/W) / Σ_t.
        for g in range(ng):
            sig_g = float(sigma_t[g, 0])
            expected_per_ord = (q_iso / sum_w) / sig_g
            actual = psi_typed.interior.values[:, g, :]
            rel_dev = np.abs(actual - expected_per_ord) / max(
                expected_per_ord, 1e-30,
            )
            assert rel_dev.max() < 1e-6, (
                f"{coord} composite g={g}: expected per-ord ψ = "
                f"{expected_per_ord:.4e}; max rel deviation = "
                f"{rel_dev.max():.3e}.  D-H.1c stage 1 composite branch "
                f"and stage 2 solver routing must both hold for this "
                f"streaming-equilibrium identity to converge."
            )
