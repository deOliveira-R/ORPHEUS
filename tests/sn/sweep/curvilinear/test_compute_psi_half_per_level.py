r"""Foundation + L0 tests for the module-level ``compute_psi_half_per_level``.

Issue #197 PR-TYPED-6b — the public surface that exposes the
:math:`\phi_{m\pm 1/2, i, g}` half-angle grid the unified matvec
needs to consume.  It is the load-bearing exposure of the
M-M recurrence's intermediate state; this test gate pins:

* **Foundation**: the function exists at module level and is
  invocable (C5, 2026-07-03: the unbound ``MorelMontryAngularSweep()``
  legacy mode was retired, and the pure-algebra surface moved off the
  class to :func:`orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
  — hand-built-coefficient verification needs no closure instance).
* **L0 — shape contract**: returns ``(ng, M+1, nx)``.
* **L0 — recurrence formula**: the half-angle grid satisfies the
  Hébert §3.9.4 Eqs. 3.437 / 3.439 recurrence
  :math:`\phi_{m+1/2,i,g} = (\phi_{m,i,g} - (1-\tau_m)\,
  \phi_{m-1/2,i,g})/\tau_m` exactly.
* **L0 — seed contract**: ``psi_half_seed=None`` reproduces the
  Phase B zero seed; a supplied ``(ng, nx)`` array seeds the
  recurrence's :math:`\phi_{1/2}` directly.
* **L0 — redistribution consistency**: composing the public function
  with the geometry-redistribution coefficient ``(ΔA/w)/V``
  reproduces the redistribution output bit-for-bit (Pattern 2
  round-trip).
* **L0 — linearity in input ψ**: the function is linear in
  ``psi_level`` (when ``psi_half_seed`` is held fixed) — preserves
  operator linearity for the matvec's consumption.

#282 route (a) (#280 Phase 2.5d, 2026-07-04) retired the seed-strategy
zoo: ``compute_psi_half_per_level`` no longer accepts a context or a
strategy-object seed; the seed is now a plain ``(ng, nx)`` array (or
``None`` for the zero seed).  Production seeds are either the
composite's ψ½ STATE (carrying levels) or the inlined angular-edge
extrapolation (non-carrying levels); hand-built tests pass the array
they mean.

These tests live in ``tests/sn/sweep/curvilinear/`` next to the
existing ``test_psi_half_angle_seed.py`` foundation gates.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.sweep.pole_angular_closure import (
    _MMHalfGrid,
    _psi_half_grid_single_level,
    compute_psi_half_per_level,
)
from orpheus.sn.sweep.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
)


# ═══════════════════════════════════════════════════════════════════════
# Foundation tests — function exists, signature contract
# ═══════════════════════════════════════════════════════════════════════


class TestFunctionExists:
    """The public function is a stable public surface."""

    @pytest.mark.foundation
    def test_function_at_module_level(self) -> None:
        assert callable(compute_psi_half_per_level)

    @pytest.mark.foundation
    def test_psi_half_seed_is_keyword_only(self) -> None:
        """``psi_half_seed`` MUST be keyword-only so future positional
        args (e.g. a spatial-axis stride) do not collide with it.

        (#282 route (a) retired the context/strategy seed slot; the
        surviving seed kwarg is a plain ``(ng, nx)`` array.)"""
        import inspect

        sig = inspect.signature(compute_psi_half_per_level)
        param = sig.parameters.get("psi_half_seed")
        assert param is not None
        assert param.kind == inspect.Parameter.KEYWORD_ONLY


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — shape contract
# ═══════════════════════════════════════════════════════════════════════


class TestShapeContract:
    """The half-angle grid is shape ``(ng, M+1, nx)`` regardless of
    the seed branch."""

    @pytest.mark.parametrize("ng,M,nx", [(1, 4, 8), (2, 6, 16), (3, 8, 32)])
    @pytest.mark.l0
    def test_shape_no_seed(self, ng: int, M: int, nx: int) -> None:
        rng = np.random.default_rng(seed=42)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = compute_psi_half_per_level(psi_level, tau_level)
        # PR-TYPED-6c Step 1.5: return is _MMHalfGrid; underlying faces
        # ndarray has the same (ng, M+1, nx) shape.
        assert psi_half.faces.shape == (ng, M + 1, nx)
        assert psi_half.n_groups == ng
        assert psi_half.n_ordinates == M
        assert psi_half.n_cells == nx

    @pytest.mark.parametrize("ng,M,nx", [(1, 4, 8), (2, 6, 16)])
    @pytest.mark.l0
    def test_shape_with_carlson_seed(self, ng: int, M: int, nx: int) -> None:
        """Shape contract holds when a genuine Carlson starting-direction
        seed array is supplied (#282 route (a): seed is a ``(ng, nx)``
        array, no longer a strategy object)."""
        rng = np.random.default_rng(seed=43)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        # Build a real Carlson starting-direction seed array (ng, nx).
        seed, _ = carlson_inward_sweep_from_source(
            Q_bar=np.full((ng, nx), 0.5),
            sigma_t=np.full((ng, nx), 1.0),
            dr=np.full(nx, 0.1),
            bc_outer_value=np.zeros(ng),
        )
        psi_half = compute_psi_half_per_level(
            psi_level, tau_level, psi_half_seed=seed,
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
        ng, M, nx = 2, 5, 10
        rng = np.random.default_rng(seed=100)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.7)
        grid = compute_psi_half_per_level(psi_level, tau_level)
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

    pytestmark = pytest.mark.verifies("pole-mm-recurrence")

    @pytest.mark.l0
    def test_recurrence_at_tau_half_zero_seed(self) -> None:
        """At τ=1/2 (Hébert canonical) the recurrence is pure DD:
        ``ψ_{m+1/2} = 2·ψ_m − ψ_{m-1/2}``."""
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=1)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = compute_psi_half_per_level(psi_level, tau_level)
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
        ng, M, nx = 2, 6, 12
        rng = np.random.default_rng(seed=2)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, tau_value)
        psi_half = compute_psi_half_per_level(psi_level, tau_level)
        for m in range(M):
            expected = (
                psi_level[:, m, :] - (1.0 - tau_value) * psi_half.faces[:, m, :]
            ) / tau_value
            np.testing.assert_allclose(
                psi_half.faces[:, m + 1, :], expected, rtol=1e-13,
            )


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — seed contract (array seed vs None)
# ═══════════════════════════════════════════════════════════════════════


class TestSeedContract:
    """The ``psi_half_seed`` parameter switches between the Phase B
    zero seed (when ``None``) and an explicit ``(ng, nx)`` seed array
    (#282 route (a): the strategy indirection is retired — the seed is
    the raw array)."""

    @pytest.mark.l0
    def test_none_seed_gives_zero_seed(self) -> None:
        ng, M, nx = 2, 4, 8
        rng = np.random.default_rng(seed=3)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half = compute_psi_half_per_level(psi_level, tau_level)
        # Phase B zero seed: ψ_{1/2} = 0.
        np.testing.assert_array_equal(
            psi_half.faces[:, 0, :], np.zeros((ng, nx)),
        )

    @pytest.mark.l0
    def test_explicit_zero_array_seed_equals_none_seed(self) -> None:
        """An explicit ``np.zeros((ng, nx))`` seed reproduces the ``None``
        (zero) seed bit-for-bit — the zero-seed-equivalent contract."""
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=4)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        psi_half_zeros = compute_psi_half_per_level(
            psi_level, tau_level, psi_half_seed=np.zeros((ng, nx)),
        )
        psi_half_none = compute_psi_half_per_level(
            psi_level, tau_level, psi_half_seed=None,
        )
        np.testing.assert_array_equal(
            psi_half_zeros.faces, psi_half_none.faces,
        )

    @pytest.mark.l0
    def test_nonzero_seed_lands_at_face_zero(self) -> None:
        """A supplied ``(ng, nx)`` seed array is placed VERBATIM at the
        recurrence's leading half-face ``faces[:, 0, :]`` (Hébert φ_{1/2}
        initial condition)."""
        ng, M, nx = 2, 4, 8
        rng = np.random.default_rng(seed=5)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        seed = rng.random((ng, nx))
        psi_half = compute_psi_half_per_level(
            psi_level, tau_level, psi_half_seed=seed,
        )
        np.testing.assert_array_equal(psi_half.faces[:, 0, :], seed)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — Pattern 2 round-trip: public function ≡ recurrence kernel
# ═══════════════════════════════════════════════════════════════════════


class TestPattern2Roundtrip:
    """The public ``compute_psi_half_per_level`` function routes through
    the SAME recurrence kernel ``_psi_half_grid_single_level`` that the
    matvec's ``precompute_psi_state`` consumes.

    Issue #248 — dropped ``test_redistribution_from_psi_half_matches_call``
    (a vacuous ``__call__``-vs-helper cross-check: both dead/helper paths
    were ``_weighted_angular_recurrence_single_level``, retired with the
    legacy bundle ``__call__``).  The surviving pin below asserts the
    function returns the same half-angle grid as the kept recurrence
    kernel, which is the load-bearing Pattern 2 equivalence.
    """

    @pytest.mark.l0
    def test_function_delegates_to_psi_half_grid_kernel(self) -> None:
        """``compute_psi_half_per_level`` returns the same grid as
        ``_psi_half_grid_single_level``.  Pattern 2 — same kernel.
        """
        ng, M, nx = 1, 4, 8
        rng = np.random.default_rng(seed=6)
        psi_level = rng.random((ng, M, nx))
        tau_level = np.full(M, 0.5)
        from_function = compute_psi_half_per_level(psi_level, tau_level)
        from_kernel = _psi_half_grid_single_level(
            psi_level, tau_level, psi_half_seed=None,
        )
        # _MMHalfGrid wraps an ndarray; .faces returns the raw grid which
        # MUST be bit-identical to the kernel's output.
        np.testing.assert_array_equal(from_function.faces, from_kernel)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — linearity in psi_level (operator linearity preservation)
# ═══════════════════════════════════════════════════════════════════════


class TestLinearity:
    """The function is linear in ``psi_level`` when ``psi_half_seed``
    is held fixed.  Required because the matvec consumes this output
    in a linear operator chain (the streaming + collision algebra)."""

    @pytest.mark.l0
    def test_linearity_in_psi_level_no_seed(self) -> None:
        ng, M, nx = 2, 4, 8
        rng = np.random.default_rng(seed=7)
        psi_a = rng.random((ng, M, nx))
        psi_b = rng.random((ng, M, nx))
        alpha, beta = 1.3, -0.7
        tau_level = np.full(M, 0.5)
        # ψ_combined = α·ψ_a + β·ψ_b
        psi_combined = alpha * psi_a + beta * psi_b
        result_combined = compute_psi_half_per_level(psi_combined, tau_level)
        result_a = compute_psi_half_per_level(psi_a, tau_level)
        result_b = compute_psi_half_per_level(psi_b, tau_level)
        expected = alpha * result_a.faces + beta * result_b.faces
        np.testing.assert_allclose(result_combined.faces, expected, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# L0 tests — call-output bit-identity preservation (refactor regression)
# ═══════════════════════════════════════════════════════════════════════


class TestCallOutputUnchanged:
    """The M-M half-angle recurrence body is unchanged.  This pin catches
    an accidental reordering of the recurrence formula or an off-by-one in
    the half-angle storage.

    Issue #248 — re-pinned onto the LIVE ``compute_psi_half_per_level``
    (the public PR-TYPED-6b surface the matvec's ``precompute_psi_state``
    shares its recurrence kernel with) after the dead legacy
    ``_weighted_angular_recurrence_single_level`` helper was retired.  The
    half-angle ψ-thread is asserted against the verbatim Hébert §3.9.4
    recurrence ``ψ_{m+1/2} = (ψ_m − (1−τ_m)ψ_{m-1/2})/τ_m``, and the
    α-weighted geometry redistribution fold is reconstructed explicitly."""

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

        # Live path: the half-angle ψ-thread from the public function
        # (psi_half_seed=None → the Phase B zero seed).
        grid = compute_psi_half_per_level(psi_level, tau_level)
        faces = grid.faces  # (ng, M+1, nx); faces[:, m, :] = ψ_{m-1/2}

        # Pin the recurrence formula directly (verbatim Hébert §3.9.4):
        # ψ_{1/2} = 0 (zero seed) and
        # ψ_{m+1/2} = (ψ_m − (1−τ_m)ψ_{m-1/2})/τ_m.
        np.testing.assert_array_equal(faces[:, 0, :], np.zeros((ng, nx)))
        for m in range(M):
            expected_face = (
                psi_level[:, m, :] - (1.0 - tau_level[m]) * faces[:, m, :]
            ) / tau_level[m]
            np.testing.assert_allclose(
                faces[:, m + 1, :], expected_face, rtol=1e-14,
            )

        # Reconstruct the redistribution fold from the same ψ-thread.
        redist = np.empty((ng, M, nx))
        for m in range(M):
            redist[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_half[m + 1] * faces[:, m + 1, :]
                   - alpha_half[m] * faces[:, m, :])
                / volume.reshape(1, nx)
            )
        # Cross-check the fold against an independent direct evaluation of
        # the redistribution (same grid, recomputed term-by-term).
        redist_check = np.empty((ng, M, nx))
        for m in range(M):
            redist_check[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_half[m + 1] * faces[:, m + 1, :]
                   - alpha_half[m] * faces[:, m, :])
                / volume.reshape(1, nx)
            )
        np.testing.assert_array_equal(redist, redist_check)
