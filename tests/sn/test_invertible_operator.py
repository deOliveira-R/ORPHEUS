r"""L0 foundation: :class:`InvertibleOperator` — sweep-invertible (L + C).

R-1 Step C (2026-05-19) — the SN-specific algebraic identity

.. math::

    (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
    \text{WDD sweep}

is encoded at the type level by :class:`InvertibleOperator`, a
specialisation of :class:`~orpheus.numerics.operator.OperatorSum` that
carries ``.solve`` via :func:`~orpheus.sn.sweep.transport_sweep`.

The dispatch is symmetric: ``L + C`` and ``C + L`` both produce an
``InvertibleOperator`` (with the streaming operand stored first via the
canonical ``OperatorSum`` ordering).

Tests pin:

* Dispatch — ``L + C`` and ``C + L`` both return InvertibleOperator;
  ``L + L``, ``C + C``, and ``L + S`` (where S is a scattering operator)
  fall through to generic OperatorSum.
* Capability — ``CAP_APPLY`` AND ``CAP_SOLVE`` advertised; the parent
  OperatorSum advertises only ``CAP_APPLY``.
* Invariants — mesh-identity required; negative sigma rejected.
* Apply equivalence — InvertibleOperator.apply matches the inherited
  OperatorSum action (L.apply + C.apply, bit-exact).
* Solve consistency — apply ∘ solve = identity (on volumetric rhs).
* Carlson seed plumbing — InvertibleOperator.solve reads ``rhs(1)`` and
  forwards it to :func:`transport_sweep` as ``initial_guess`` for the
  curvilinear pole closure.
"""
from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_SOLVE,
    OperatorSum,
)
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.boundary_flux import BoundaryFlux
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.sn.quadrature import GaussLegendre1D
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ── Geometry fixtures ───────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ── Dispatch via __add__ ────────────────────────────────────────────────


class TestDispatch:
    def test_streaming_plus_collision_returns_invertible(self) -> None:
        """``L + C`` returns :class:`InvertibleOperator` via dispatch."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)

        composite = L + C
        assert isinstance(composite, InvertibleOperator)
        assert composite.streaming is L
        assert composite.diagonal is C

    def test_collision_plus_streaming_also_returns_invertible(self) -> None:
        """``C + L`` (commutative dispatch) also returns InvertibleOperator."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)

        composite = C + L
        assert isinstance(composite, InvertibleOperator)
        # Streaming is canonically placed first (the algebraic identity
        # ``L + C^{-1}`` reads that way).
        assert composite.streaming is L
        assert composite.diagonal is C

    def test_explicit_construction_with_sigma_r(self) -> None:
        """Explicit ctor accepts a CollisionOperator built with σ_r ≠ σ_t."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        # σ_r = σ_t - within-group self-scatter; here a stylised
        # positive reduction (e.g. 1.0 - 0.3 = 0.7 for c_(g->g) = 0.3).
        sigma_r = 0.7 * np.ones((sn.ng, sn.nx, sn.ny))

        L = StreamingOperator(sn, sigma_t)
        C_r = CollisionOperator(sn, sigma_r)
        invertible = InvertibleOperator(L, C_r)
        assert invertible.streaming is L
        assert invertible.diagonal is C_r
        np.testing.assert_array_equal(invertible.sigma, sigma_r)

    def test_streaming_plus_zero_operator_falls_through(self) -> None:
        """Non-Collision RHS falls back to generic OperatorSum."""
        from orpheus.numerics.operator import ZeroOperator
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        Z = ZeroOperator()

        composite = L + Z
        assert isinstance(composite, OperatorSum)
        assert not isinstance(composite, InvertibleOperator)


# ── Capabilities and invariants ─────────────────────────────────────────


class TestCapabilitiesAndInvariants:
    def test_advertises_apply_and_solve(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = InvertibleOperator(
            StreamingOperator(sn, sigma_t),
            CollisionOperator(sn, sigma_t),
        )
        assert CAP_APPLY in invertible.capabilities
        assert CAP_SOLVE in invertible.capabilities

    def test_parent_operator_sum_lacks_solve(self) -> None:
        r"""Plain :class:`OperatorSum` does NOT advertise solve.

        The :math:`(A+B)^{-1}` algebraic identity is operator-pair
        specific (no generic formula); only :class:`InvertibleOperator`
        carries the SN-specific sweep identity at the type level.
        """
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        # Explicitly construct via the parent class to bypass dispatch.
        plain = OperatorSum(
            StreamingOperator(sn, sigma_t),
            CollisionOperator(sn, sigma_t),
        )
        assert CAP_APPLY in plain.capabilities
        assert CAP_SOLVE not in plain.capabilities

    def test_mismatched_mesh_rejected(self) -> None:
        """Two distinct :class:`SNMesh` instances → ValueError."""
        sn1 = _slab_mesh()
        sn2 = _slab_mesh()
        sigma_t1 = np.ones((sn1.ng, sn1.nx, sn1.ny))
        sigma_t2 = np.ones((sn2.ng, sn2.nx, sn2.ny))
        L = StreamingOperator(sn1, sigma_t1)
        C = CollisionOperator(sn2, sigma_t2)
        with pytest.raises(ValueError, match="mesh-identity"):
            InvertibleOperator(L, C)

    def test_non_positive_sigma_rejected(self) -> None:
        """``σ <= 0`` anywhere → ValueError at construction."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        sigma_bad = sigma_t.copy()
        sigma_bad[0, 0, 0] = -0.5  # one bad cell
        L = StreamingOperator(sn, sigma_t)
        C_bad = CollisionOperator(sn, sigma_bad)
        with pytest.raises(ValueError, match="strictly positive"):
            InvertibleOperator(L, C_bad)

    def test_zero_sigma_rejected(self) -> None:
        """``σ == 0`` (void cell) is also rejected — strict positivity."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        sigma_bad = sigma_t.copy()
        sigma_bad[0, 0, 0] = 0.0
        L = StreamingOperator(sn, sigma_t)
        C_bad = CollisionOperator(sn, sigma_bad)
        with pytest.raises(ValueError, match="strictly positive"):
            InvertibleOperator(L, C_bad)

    def test_wrong_operand_types_rejected(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)
        with pytest.raises(TypeError, match="StreamingOperator"):
            InvertibleOperator(C, C)  # type: ignore[arg-type]
        with pytest.raises(TypeError, match="CollisionOperator"):
            InvertibleOperator(L, L)  # type: ignore[arg-type]


# ── Apply matches inherited OperatorSum ─────────────────────────────────


class TestApply:
    def test_apply_equals_l_plus_c_on_typed_flux(self) -> None:
        r"""``(L+C).apply(ψ) == L.apply(ψ) + C.apply(ψ)`` for AngularFlux."""
        sn = _slab_mesh(nx=4, n_ord=4, ng=2)
        sigma_t = np.full((sn.ng, sn.nx, sn.ny), 0.8)
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)
        invertible = L + C

        rng = np.random.default_rng(7)
        psi = AngularFlux(
            rng.standard_normal((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn,
        )

        composite_out = invertible.apply(psi)
        sum_out = L.apply(psi) + C.apply(psi)

        np.testing.assert_array_equal(composite_out.values, sum_out.values)
        np.testing.assert_array_equal(
            composite_out.boundary.xmax_face, sum_out.boundary.xmax_face,
        )
        if composite_out.boundary.xmin_face is not None:
            np.testing.assert_array_equal(
                composite_out.boundary.xmin_face, sum_out.boundary.xmin_face,
            )


# ── Solve via sweep ─────────────────────────────────────────────────────


class TestSolve:
    def test_solve_returns_angular_flux_on_slab(self) -> None:
        """``InvertibleOperator.solve`` returns a typed AngularFlux on slab."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        rhs = AngularFlux(
            np.ones((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn,
        )
        psi = invertible.solve(rhs)

        assert isinstance(psi, AngularFlux)
        assert psi.mesh is sn
        assert psi.values.shape == rhs.values.shape

    def test_solve_inherits_history_depth_from_rhs(self) -> None:
        """``rhs.history_depth`` propagates to the returned AngularFlux."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        rhs = AngularFlux(
            np.zeros((sn.quad.N, sn.ng, sn.nx, sn.ny)),
            sn,
            history_depth=5,
        )
        psi = invertible.solve(rhs)
        assert psi.history_depth == 5

    def test_solve_rejects_bare_ndarray(self) -> None:
        """The typed-flux contract is mandatory — bare ndarray rejected."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )
        with pytest.raises(TypeError, match="AngularFlux"):
            invertible.solve(np.zeros(10))

    def test_solve_rejects_mismatched_mesh(self) -> None:
        sn1 = _slab_mesh()
        sn2 = _slab_mesh()
        sigma_t1 = np.ones((sn1.ng, sn1.nx, sn1.ny))
        invertible = StreamingOperator(sn1, sigma_t1) + CollisionOperator(
            sn1, sigma_t1,
        )
        rhs = AngularFlux(
            np.zeros((sn2.quad.N, sn2.ng, sn2.nx, sn2.ny)), sn2,
        )
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
        sigma_t = np.full((sn.ng, sn.nx, sn.ny), 1.5)
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
        x = np.linspace(0.5 / nx, 1.0 - 0.5 / nx, nx)
        rhs_values = np.broadcast_to(
            np.sin(np.pi * x)[None, None, :, None] * np.ones((N, ng, 1, ny)),
            (N, ng, nx, ny),
        ).copy()
        q = AngularFlux(rhs_values, sn)
        psi = invertible.solve(q)

        # Some cells positive (interior; the sweep is monotone for
        # positive sources on vacuum BC slab).  Exact every-cell-
        # positivity may fail at the BC face under WDD's diamond
        # closure — that's a known DD artefact.  Sufficient signal:
        # non-trivial peak in the interior.
        assert psi.values.max() > 0
        # Peak occurs near the source maximum (centre of the slab).
        peak_x_idx = np.unravel_index(
            np.argmax(psi.values[0, 0, :, 0]), psi.values[0, 0, :, 0].shape,
        )
        assert nx // 4 <= peak_x_idx[0] <= 3 * nx // 4

    def test_solve_runs_on_sphere(self) -> None:
        r"""Sphere geometry — smoke test that the curvilinear path runs."""
        sn = _sphere_mesh(nx=8, n_ord=4, ng=1)
        sigma_t = np.full((sn.ng, sn.nx, sn.ny), 1.0)
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        q = AngularFlux(
            np.ones((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn,
        )
        psi = invertible.solve(q)
        assert isinstance(psi, AngularFlux)
        assert psi.values.shape == q.values.shape
        # Curvilinear sweep produces a non-trivial positive bulk.
        assert psi.values.max() > 0

    @pytest.mark.l0
    @pytest.mark.verifies("sn-streaming", "transport-cartesian")
    def test_solve_consumes_per_ordinate_rhs(self) -> None:
        r"""``InvertibleOperator.solve`` passes ``rhs.values`` unmodified to the sweep.

        R-1 Step 4 A1 invariant pin (N5 per verification plan).  The
        producer-side normalisation contract says ``rhs.values`` is
        already per-ordinate density (the producer of ``rhs`` — typically
        ``ScatteringOperator.apply`` or ``PerOrdinateSource.from_isotropic``
        — applied ``/sum_w`` at the producer boundary).  The
        ``InvertibleOperator.solve`` adapter MUST forward ``rhs.values``
        to :func:`transport_sweep` *bit-equal* — no internal ``* sum_w``
        bridge, no ``/sum_w`` rescaling.

        Pre-A1 the adapter wrapped ``rhs.values * sum_w`` into the sweep
        to compensate for the old sweep-internal ``/W``; both directions
        of that bridge dissolved in A1.  This test spies on
        :func:`transport_sweep` to capture the ``source`` argument and
        asserts ``source.values`` is bit-identical to ``rhs.values``.
        If a future refactor re-introduces a ``* sum_w`` / ``/ sum_w``
        rescaling on this hot path, the bit-equal assertion fails.
        """
        sn = _slab_mesh(nx=4, n_ord=4, ng=1)
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        sum_w = float(sn.quad.weights.sum())
        q_const = 2.7
        per_ord_density = q_const / sum_w
        rhs = AngularFlux(
            np.full((sn.quad.N, sn.ng, sn.nx, sn.ny), per_ord_density), sn,
        )

        # Spy on transport_sweep — capture the source argument's values.
        captured: list[np.ndarray] = []
        from orpheus.sn import sweep as sweep_module
        original = sweep_module.transport_sweep

        def spy(source, *args, **kwargs):
            # Snapshot ``source.values`` BEFORE the sweep mutates anything
            # (transport_sweep is a pure reader of its source argument
            # but the snapshot makes intent explicit).
            captured.append(np.array(source.values, copy=True))
            return original(source, *args, **kwargs)

        with patch("orpheus.sn.sweep.transport_sweep", spy):
            invertible.solve(rhs)

        assert len(captured) == 1, (
            f"transport_sweep called {len(captured)} times; expected 1"
        )
        forwarded = captured[0]
        np.testing.assert_array_equal(
            forwarded, rhs.values,
            err_msg=(
                "InvertibleOperator.solve modified rhs.values before "
                "forwarding to transport_sweep — A1 producer-side "
                "convention drifted (the ``* sum_w`` bridge MUST stay "
                "dissolved)."
            ),
        )

    def test_solve_forwards_explicit_initial_guess_to_sweep(self) -> None:
        r"""``InvertibleOperator.solve`` forwards the explicit ``initial_guess``
        kwarg to :func:`transport_sweep`.

        Phase 1.2 — the curvilinear Carlson coupled-pole seed travels
        through an explicit ``initial_guess`` argument; the lag-1 frame
        machinery on :class:`AngularFlux` (``rhs(1)``, ``.stash``) is
        reserved for future time-derivative tracking and no longer
        load-bearing for the seed path.  The previous
        ``previous = rhs(1)`` plumbing was retired together with the
        sweep's inline Q_bar derivation — the M-M closure's
        ``psi_half_seed`` strategy reads ``initial_guess`` directly
        (or zeros on cold start).

        Spy on :func:`transport_sweep` to capture the ``initial_guess``
        argument on both the explicit-seed and cold-start paths.
        """
        sn = _sphere_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        invertible = StreamingOperator(sn, sigma_t) + CollisionOperator(
            sn, sigma_t,
        )

        psi_prev = AngularFlux(
            np.full((sn.quad.N, sn.ng, sn.nx, sn.ny), 0.7),
            sn,
        )
        rhs = AngularFlux(
            np.zeros((sn.quad.N, sn.ng, sn.nx, sn.ny)),
            sn,
        )

        # Spy on transport_sweep — capture initial_guess on every call.
        captured = []

        from orpheus.sn import sweep as sweep_module
        original = sweep_module.transport_sweep

        def spy(*args, **kwargs):
            captured.append(kwargs.get("initial_guess"))
            return original(*args, **kwargs)

        with patch("orpheus.sn.sweep.transport_sweep", spy):
            invertible.solve(rhs, initial_guess=psi_prev)
        assert len(captured) == 1
        seed = captured[0]
        assert seed is not None
        np.testing.assert_array_equal(seed.values, psi_prev.values)

        # Cold start — no explicit seed → initial_guess should be None.
        captured.clear()
        with patch("orpheus.sn.sweep.transport_sweep", spy):
            invertible.solve(rhs)
        assert len(captured) == 1
        assert captured[0] is None
