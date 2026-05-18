r"""Foundation + L0 tests for ``MorelMontryAngularSweep.compute_psi_half_per_level``.

Issue #197 PR-TYPED-6b — the public method that exposes the
:math:`\phi_{m\pm 1/2, i, g}` half-angle grid the unified matvec
needs to consume.  The method is the load-bearing exposure of the
M-M recurrence's intermediate state; this test gate pins:

* **Foundation**: the method exists on
  :class:`MorelMontryAngularSweep` and is invocable.
* **L0 — shape contract**: returns ``(ng, M+1, nx)``.
* **L0 — recurrence formula**: the half-angle grid satisfies the
  Hébert §3.9.4 Eqs. 3.437 / 3.439 recurrence
  :math:`\phi_{m+1/2,i,g} = (\phi_{m,i,g} - (1-\tau_m)\,
  \phi_{m-1/2,i,g})/\tau_m` exactly.
* **L0 — seed contract**: ``carlson_context=None`` reproduces the
  Phase B zero seed; ``carlson_context=ctx`` invokes
  :attr:`psi_half_seed` to build the recurrence's :math:`\phi_{1/2}`.
* **L0 — redistribution consistency**: composing the public method
  with the geometry-redistribution coefficient ``(ΔA/w)/V``
  reproduces the redistribution output of the existing ``__call__``
  bit-for-bit (Pattern 2 round-trip).
* **L0 — linearity in input ψ**: the method is linear in
  ``psi_level`` (when ``carlson_context`` is held fixed) — preserves
  operator linearity for the matvec's consumption.

These tests live in ``tests/sn/spatial/`` next to the existing
``test_psi_half_angle_seed.py`` foundation gates.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.pole_angular_closure import (
    MorelMontryAngularSweep,
    _MMHalfGrid,
    _mm_psi_half_grid_single_level,
    _mm_weighted_angular_recurrence_single_level,
)
from orpheus.sn.spatial.psi_half_angle_seed import (
    CarlsonInwardSweep,
    CarlsonSweepContext,
    ZeroSeed,
)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — method exists, signature contract
# ═══════════════════════════════════════════════════════════════════════


class TestMethodExists:
    """The public method is a stable public surface."""

    @pytest.mark.foundation
    def test_method_attached_to_morel_montry(self) -> None:
        sweep = MorelMontryAngularSweep()
        assert hasattr(sweep, "compute_psi_half_per_level")
        assert callable(sweep.compute_psi_half_per_level)

    @pytest.mark.foundation
    def test_method_consumes_keyword_only_carlson_context(self) -> None:
        """``carlson_context`` MUST be keyword-only so future positional
        args (e.g. spatial-axis stride) do not collide with it."""
        import inspect

        sweep = MorelMontryAngularSweep()
        sig = inspect.signature(sweep.compute_psi_half_per_level)
        carlson_param = sig.parameters.get("carlson_context")
        assert carlson_param is not None
        assert carlson_param.kind == inspect.Parameter.KEYWORD_ONLY


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — shape contract
# ═══════════════════════════════════════════════════════════════════════


class TestShapeContract:
    """The half-angle grid is shape ``(ng, M+1, nx)`` regardless of
    the carlson_context branch."""

    @pytest.mark.parametrize("ng,M,nx", [(1, 4, 8), (2, 6, 16), (3, 8, 32)])
    @pytest.mark.l0
    def test_shape_no_carlson(self, ng: int, M: int, nx: int) -> None:
        sweep = MorelMontryAngularSweep()
        rng = np.random.default_rng(seed=42)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        # PR-TYPED-6c Step 1.5: return is _MMHalfGrid; underlying faces
        # ndarray has the same (ng, M+1, nx) shape.
        assert psi_half.faces.shape == (ng, M + 1, nx)
        assert psi_half.n_groups == ng
        assert psi_half.n_ordinates == M
        assert psi_half.n_cells == nx

    @pytest.mark.parametrize("ng,M,nx", [(1, 4, 8), (2, 6, 16)])
    @pytest.mark.l0
    def test_shape_with_carlson(self, ng: int, M: int, nx: int) -> None:
        sweep = MorelMontryAngularSweep()
        rng = np.random.default_rng(seed=43)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        # Build a valid carlson context.
        ctx = CarlsonSweepContext(
            sigma_t=np.full((ng, nx), 1.0),
            dr=np.full(nx, 0.1),
            mu_quad=np.linspace(-1.0, 1.0, M),
            weights=np.full(M, 2.0 / M),
            bc_outer_value=np.zeros(ng),
        )
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=ctx,
        )
        assert psi_half.faces.shape == (ng, M + 1, nx)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — _MMHalfGrid typed accessor contracts (PR-TYPED-6c Step 1.5)
# ═══════════════════════════════════════════════════════════════════════


class TestMMHalfGridAccessors:
    """The :class:`_MMHalfGrid` wrapper exposes named accessors so the
    matvec consumer cannot off-by-one between upstream and downstream
    half-faces.

    Pattern 4 — illegal states unrepresentable. The matvec's
    ``psi_angular_upstream`` argument is populated via
    :meth:`_MMHalfGrid.upstream` (semantic, not raw indexing).
    """

    def _build_grid(self) -> tuple[_MMHalfGrid, int, int, int]:
        sweep = MorelMontryAngularSweep()
        ng, M, nx = 2, 5, 10
        rng = np.random.default_rng(seed=100)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.7)
        grid = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        return grid, ng, M, nx

    @pytest.mark.l0
    def test_return_type_is_mm_half_grid(self) -> None:
        """``compute_psi_half_per_level`` returns an :class:`_MMHalfGrid`,
        not a raw ndarray."""
        grid, _, _, _ = self._build_grid()
        assert isinstance(grid, _MMHalfGrid)

    @pytest.mark.l0
    def test_shape_properties(self) -> None:
        """``n_groups``, ``n_ordinates``, ``n_cells`` decode the
        underlying ``faces.shape``."""
        grid, ng, M, nx = self._build_grid()
        assert grid.n_groups == ng
        assert grid.n_ordinates == M
        assert grid.n_cells == nx
        assert grid.faces.shape == (ng, M + 1, nx)

    @pytest.mark.l0
    def test_upstream_per_ordinate_drops_trailing_face(self) -> None:
        """``upstream_per_ordinate`` is shape ``(ng, M, nx)`` — the
        first ``M`` faces, which are the upstream-per-ordinate set."""
        grid, ng, M, nx = self._build_grid()
        upstream = grid.upstream_per_ordinate
        assert upstream.shape == (ng, M, nx)
        np.testing.assert_array_equal(upstream, grid.faces[:, :-1, :])

    @pytest.mark.l0
    def test_downstream_per_ordinate_drops_leading_face(self) -> None:
        """``downstream_per_ordinate`` is shape ``(ng, M, nx)`` — the
        last ``M`` faces, which are the downstream-per-ordinate set."""
        grid, ng, M, nx = self._build_grid()
        downstream = grid.downstream_per_ordinate
        assert downstream.shape == (ng, M, nx)
        np.testing.assert_array_equal(downstream, grid.faces[:, 1:, :])

    @pytest.mark.l0
    def test_upstream_indexing(self) -> None:
        """``upstream(m)`` returns ``faces[:, m, :]`` — the upstream
        face of ordinate ``m`` (= downstream of ``m-1``)."""
        grid, _, M, _ = self._build_grid()
        for m in range(M):
            np.testing.assert_array_equal(
                grid.upstream(m), grid.faces[:, m, :],
            )

    @pytest.mark.l0
    def test_downstream_indexing(self) -> None:
        """``downstream(m)`` returns ``faces[:, m+1, :]`` — the
        downstream face of ordinate ``m`` (= upstream of ``m+1``)."""
        grid, _, M, _ = self._build_grid()
        for m in range(M):
            np.testing.assert_array_equal(
                grid.downstream(m), grid.faces[:, m + 1, :],
            )

    @pytest.mark.l0
    def test_upstream_downstream_chain(self) -> None:
        """``downstream(m) == upstream(m+1)`` for adjacent ordinates —
        the half-face identity that makes the M-M recurrence consistent."""
        grid, _, M, _ = self._build_grid()
        for m in range(M - 1):
            np.testing.assert_array_equal(
                grid.downstream(m), grid.upstream(m + 1),
            )

    @pytest.mark.l0
    def test_mm_half_grid_is_frozen(self) -> None:
        """``_MMHalfGrid`` is a frozen dataclass — its ``faces`` attribute
        cannot be reassigned (Pattern 4: prevent accidental mutation)."""
        grid, _, _, _ = self._build_grid()
        with pytest.raises((AttributeError, TypeError)):
            grid.faces = np.zeros_like(grid.faces)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — recurrence formula (Hébert Eqs. 3.437 / 3.439)
# ═══════════════════════════════════════════════════════════════════════


class TestRecurrenceFormula:
    """The half-angle grid satisfies the M-M recurrence exactly.

    For each :math:`m \\in \\{0, \\ldots, M-1\\}`:
        :math:`\\phi_{m+1/2,i,g} = (\\phi_{m,i,g} - (1-\\tau_m)
        \\phi_{m-1/2,i,g}) / \\tau_m`.
    """

    @pytest.mark.l0
    def test_recurrence_at_tau_half_zero_seed(self) -> None:
        """At τ=1/2 (Hébert canonical) the recurrence is pure DD:
        ``ψ_{m+1/2} = 2·ψ_m − ψ_{m-1/2}``."""
        sweep = MorelMontryAngularSweep()
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=1)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        # Zero seed: φ_{1/2} = 0.
        assert np.array_equal(psi_half.faces[:, 0, :], np.zeros((ng, nx)))
        # Verify pure-DD recurrence φ_{m+1/2} = 2φ_m − φ_{m-1/2}.
        for m in range(M):
            expected = 2.0 * psi_level[:, m, :] - psi_half.faces[:, m, :]
            np.testing.assert_allclose(psi_half.faces[:, m + 1, :], expected, rtol=1e-14)

    @pytest.mark.l0
    @pytest.mark.parametrize("tau_value", [0.5, 0.6, 0.75, 0.9, 1.0])
    def test_recurrence_at_arbitrary_tau(self, tau_value: float) -> None:
        """For any τ ∈ [1/2, 1] the M-M weighted DD recurrence
        holds: ``ψ_{m+1/2} = (ψ_m − (1−τ)·ψ_{m-1/2})/τ``."""
        sweep = MorelMontryAngularSweep()
        ng, M, nx = 2, 6, 12
        rng = np.random.default_rng(seed=2)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, tau_value)
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        for m in range(M):
            expected = (
                psi_level[:, m, :] - (1.0 - tau_value) * psi_half.faces[:, m, :]
            ) / tau_value
            np.testing.assert_allclose(
                psi_half.faces[:, m + 1, :], expected, rtol=1e-13,
            )


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — seed contract (carlson_context vs None)
# ═══════════════════════════════════════════════════════════════════════


class TestSeedContract:
    """The ``carlson_context`` parameter switches between Phase B
    zero seed (when ``None``) and the Carlson coupled-pole seed
    (when supplied)."""

    @pytest.mark.l0
    def test_none_context_gives_zero_seed(self) -> None:
        sweep = MorelMontryAngularSweep()
        ng, M, nx = 2, 4, 8
        rng = np.random.default_rng(seed=3)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        # Phase B zero seed: ψ_{1/2} = 0.
        np.testing.assert_array_equal(
            psi_half.faces[:, 0, :], np.zeros((ng, nx)),
        )

    @pytest.mark.l0
    def test_carlson_context_uses_psi_half_seed_strategy(self) -> None:
        """``carlson_context=ctx`` triggers
        ``self.psi_half_seed(psi_level, ctx)`` for the seed."""
        sweep = MorelMontryAngularSweep()
        # Use ZeroSeed in this sweep so the carlson-context path also
        # produces zero (sanity-check the wiring without depending on
        # CarlsonInwardSweep numerical specifics).
        sweep_zero = MorelMontryAngularSweep(psi_half_seed=ZeroSeed())
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=4)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        ctx = CarlsonSweepContext(
            sigma_t=np.full((ng, nx), 1.0),
            dr=np.full(nx, 0.1),
            mu_quad=np.linspace(-1.0, 1.0, M),
            weights=np.full(M, 2.0 / M),
            bc_outer_value=np.zeros(ng),
        )
        psi_half_with_ctx = sweep_zero.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=ctx,
        )
        psi_half_without_ctx = sweep_zero.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        # ZeroSeed returns zero regardless of context, so both branches
        # produce the same seed → bit-identical recurrence output.
        np.testing.assert_array_equal(
            psi_half_with_ctx.faces, psi_half_without_ctx.faces,
        )


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — Pattern 2 round-trip: __call__ output ≡ method composition
# ═══════════════════════════════════════════════════════════════════════


class TestPattern2Roundtrip:
    """Composing ``compute_psi_half_per_level`` with the geometry-
    redistribution coefficient ``(ΔA/w)/V`` reproduces the
    ``__call__`` redistribution output bit-for-bit.

    This is the load-bearing Pattern 2 assertion: BOTH the public
    method AND the redistribution body inside ``__call__`` route
    through the same kernel ``_mm_psi_half_grid_single_level``.
    Pre-refactor both produced the same numbers because they shared
    the recurrence body inline; post-refactor they share the
    function-call boundary; the test pins the equivalence.
    """

    @pytest.mark.l0
    def test_redistribution_from_psi_half_matches_call(self) -> None:
        sweep = MorelMontryAngularSweep(psi_half_seed=ZeroSeed())
        ng, M, nx = 2, 6, 12
        rng = np.random.default_rng(seed=5)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        # Build α-dome with the dome-positivity invariant.
        alpha_half = np.zeros(M + 1)
        for m in range(M):
            alpha_half[m + 1] = alpha_half[m] - 0.1
        # Re-base so α_{M+1/2} = 0 (dome zero at the closing boundary).
        alpha_half -= alpha_half[M]
        dAw_level = rng.random((nx, M))
        volume = rng.uniform(0.5, 1.5, nx)

        # Method path: get half-angle grid, then compose with coef.
        psi_half = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        redist_from_method = np.empty((ng, M, nx))
        for m in range(M):
            redist_from_method[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_half[m + 1] * psi_half.faces[:, m + 1, :]
                   - alpha_half[m] * psi_half.faces[:, m, :])
                / volume.reshape(1, nx)
            )

        # Free-function path: bypass strategy and call the recurrence
        # helper directly (this is what __call__ uses internally).
        redist_from_helper = _mm_weighted_angular_recurrence_single_level(
            psi_level, alpha_half, dAw_level, tau_level, volume,
            psi_half_seed=None,
        )
        np.testing.assert_allclose(
            redist_from_method, redist_from_helper, rtol=1e-14,
        )

    @pytest.mark.l0
    def test_method_delegates_to_psi_half_grid_helper(self) -> None:
        """``compute_psi_half_per_level`` returns the same grid as
        ``_mm_psi_half_grid_single_level``.  Pattern 2 — same kernel.
        """
        sweep = MorelMontryAngularSweep(psi_half_seed=ZeroSeed())
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=6)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        from_method = sweep.compute_psi_half_per_level(
            psi_level, tau_level, carlson_context=None,
        )
        from_helper = _mm_psi_half_grid_single_level(
            psi_level, tau_level, psi_half_seed=None,
        )
        # _MMHalfGrid wraps an ndarray; .faces returns the raw grid which
        # MUST be bit-identical to the free-function helper's output.
        np.testing.assert_array_equal(from_method.faces, from_helper)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — linearity in psi_level (operator linearity preservation)
# ═══════════════════════════════════════════════════════════════════════


class TestLinearity:
    """The method is linear in ``psi_level`` when ``carlson_context``
    is held fixed.  Required because the matvec consumes this output
    in a linear operator chain (the streaming + collision algebra)."""

    @pytest.mark.l0
    def test_linearity_in_psi_level_no_context(self) -> None:
        sweep = MorelMontryAngularSweep()
        ng, M, nx = 2, 4, 8
        rng = np.random.default_rng(seed=7)
        psi_a = rng.random((ng, M, nx))
        psi_b = rng.random((ng, M, nx))
        alpha, beta = 1.3, -0.7
        tau_level = np.full(M, 0.5)
        # ψ_combined = α·ψ_a + β·ψ_b
        psi_combined = alpha * psi_a + beta * psi_b
        result_combined = sweep.compute_psi_half_per_level(
            psi_combined, tau_level, carlson_context=None,
        )
        result_a = sweep.compute_psi_half_per_level(
            psi_a, tau_level, carlson_context=None,
        )
        result_b = sweep.compute_psi_half_per_level(
            psi_b, tau_level, carlson_context=None,
        )
        expected = alpha * result_a.faces + beta * result_b.faces
        np.testing.assert_allclose(result_combined.faces, expected, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — call-output bit-identity preservation (refactor regression)
# ═══════════════════════════════════════════════════════════════════════


class TestCallOutputUnchanged:
    """The refactor that factored ``_mm_psi_half_grid_single_level``
    out of ``_mm_weighted_angular_recurrence_single_level`` must not
    perturb the latter's output.  This pin catches an accidental
    reordering of the formula or an off-by-one in the half-angle
    storage."""

    @pytest.mark.l0
    @pytest.mark.parametrize("seed", [10, 11, 12])
    def test_recurrence_output_random_seed(self, seed: int) -> None:
        """Random-seed sanity probe — the recurrence body is unchanged."""
        rng = np.random.default_rng(seed=seed)
        ng, M, nx = 2, 6, 16
        psi_level = rng.random((ng, M, nx))
        tau_level = rng.uniform(0.5, 1.0, M)
        # Build a valid α-dome.
        alpha_half = np.zeros(M + 1)
        for m in range(M):
            alpha_half[m + 1] = alpha_half[m] - rng.uniform(0.05, 0.2)
        alpha_half -= alpha_half[M]
        dAw_level = rng.random((nx, M))
        volume = rng.uniform(0.5, 1.5, nx)

        # Direct call to the recurrence helper.
        redist = _mm_weighted_angular_recurrence_single_level(
            psi_level, alpha_half, dAw_level, tau_level, volume,
            psi_half_seed=None,
        )
        # Reconstruct via grid + redistribution coefficient.
        psi_half = _mm_psi_half_grid_single_level(
            psi_level, tau_level, psi_half_seed=None,
        )
        redist_check = np.empty((ng, M, nx))
        for m in range(M):
            redist_check[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_half[m + 1] * psi_half[:, m + 1, :]
                   - alpha_half[m] * psi_half[:, m, :])
                / volume.reshape(1, nx)
            )
        np.testing.assert_array_equal(redist, redist_check)
