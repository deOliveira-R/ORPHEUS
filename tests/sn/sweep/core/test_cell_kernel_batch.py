r"""The batched cell-kernel PAIR — term-level L0 + the protocol contract.

:meth:`DiamondDifference.cell_kernel_batch` (SOLVE) and
:meth:`~DiamondDifference.residual_kernel_batch` (APPLY) are the storage-free
direction fork of the SN stack — the ONLY place the solve/apply algebra
differs (S6.4).  This file pins them at the term level:

1. **Closed-form balance** — ``ψ̄ = (Q + Σ_a 2g_a ψ_in_a) / (Σ_t + Σ_a 2g_a)``
   and the WDD closure ``ψ_out_a = 2ψ̄ − ψ_in_a`` against a hand
   calculation, where ``s_a = g_a`` is the RAW down-face streaming and the
   scheme applies its diamond factor ``2 = 1/w_DD`` (#240) to BOTH the denom
   and the upstream-numerator term (the kernel returns the faces; SCATTERING
   them is the walk's job, pinned in ``test_sweep_graph.py``).
2. **Bit-identity of the ordinate vectorisation** — the kernel's batched
   left-fold against a per-ordinate Python-loop reference, per element.
3. **The apply direction's closed form, affinity, and the
   solve↔apply ROUND TRIP** — the residual vanishes at the solved ψ̄
   (the batched analogue of the per-cell ``residual``/``update`` contract).
4. **The protocol default** — a strategy that does not override the
   kernel pair fails loudly (:exc:`NotImplementedError`), never silently.

History: migrated from ``test_cell_update_batch.py`` at S6.4(e) (the
``SweepCellSlice``-packeted storage adapters ``update_batch`` /
``residual_batch`` retired — gather/scatter moved to the graph walks;
retirement = test migration).  The face-flux PLACEMENT pins moved with the
storage: ``test_sweep_graph.py`` proves the walk's gather/scatter against a
hand-rolled per-octant reference (all four octants).  Two claims DISSOLVED
by construction: the negative-octant single-cell variant (the sign only
ever selected gather INDICES — storage, not algebra) and the "residual
without probe raises" guard (the kernel signature REQUIRES ``psi_bar`` —
the illegal state is now unrepresentable, not runtime-checked).
"""

from __future__ import annotations

from typing import ClassVar

import numpy as np
import pytest

from orpheus.transport.spatial.scheme import (
    CellResult,
    DiscretizationSchemeBase,
    CellVisit,
    UpstreamState,
)
from orpheus.transport.spatial.diamond import DiamondDifference


# ─────────────────────────────────────────────────────────────────────
# Helper builders — random kernel operands (no storage, no indices)
# ─────────────────────────────────────────────────────────────────────


def _random_kernel_operands(
    *, N_oct: int, ng: int, n_diag: int, seed: int, Q_leading: int | None = None,
):
    """Random per-level kernel inputs in the kernel's shape contract:
    ``psi_in`` d-tuple of ``(N_oct, ng, n_diag)``, ``s_axes`` d-tuple of
    ``(N_oct, 1, n_diag)``, ``reaction_xs (ng, n_diag)``,
    ``Q_cells (N_oct or 1, ng, n_diag)``."""
    rng = np.random.default_rng(seed)
    psi_in = (
        rng.standard_normal((N_oct, ng, n_diag)),
        rng.standard_normal((N_oct, ng, n_diag)),
    )
    s_axes = (
        rng.uniform(0.1, 1.0, size=(N_oct, 1, n_diag)),
        rng.uniform(0.1, 1.0, size=(N_oct, 1, n_diag)),
    )
    reaction_xs = rng.uniform(0.1, 0.5, size=(ng, n_diag))
    Q_lead = N_oct if Q_leading is None else Q_leading
    Q_cells = rng.standard_normal((Q_lead, ng, n_diag))
    return psi_in, s_axes, reaction_xs, Q_cells


# ─────────────────────────────────────────────────────────────────────
# L0: Closed-form check on a 1-cell, 1-ordinate, 1-group batch
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSolveKernelClosedForm:
    """Smallest possible batch — analytical hand calculation."""

    def test_psi_avg_and_closure_match_balance_formula(self):
        r"""ψ̄ = (Q + 2g_x·ψ_in_x + 2g_y·ψ_in_y) / (Σ_t + 2g_x + 2g_y); the WDD
        closure faces are returned (the walk scatters them).

        ``s_axes`` is the RAW down-face streaming ``g`` (#240) — the scheme
        applies its diamond ``2`` internally, so the hand-calc carries the ``2``
        on every streaming term (Design A: the kernel param IS the raw ``g``).
        """
        psi_in = (np.full((1, 1, 1), 4.0), np.full((1, 1, 1), 8.0))
        s_axes = (np.full((1, 1, 1), 3.0), np.full((1, 1, 1), 5.0))  # raw g
        reaction_xs = np.full((1, 1), 2.0)
        Q_cells = np.full((1, 1, 1), 16.0)
        psi_avg, psi_out = DiamondDifference().cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        # denom = Σ_t + 2g_x + 2g_y = 2 + 2·3 + 2·5 = 18  (the scheme's diamond 2)
        # numer = Q + 2g_x·ψx + 2g_y·ψy = 16 + 2·3·4 + 2·5·8 = 120
        # ψ̄ = 120 / 18  (small integers are exact; only the division rounds)
        psi_bar = 120.0 / 18.0
        np.testing.assert_array_equal(psi_avg, psi_bar)
        # ψ_out_a = 2ψ̄ − ψ_in_a  (the diamond MEAN — not the streaming factor)
        np.testing.assert_array_equal(psi_out[0], 2.0 * psi_bar - 4.0)
        np.testing.assert_array_equal(psi_out[1], 2.0 * psi_bar - 8.0)
        # The inputs are NOT mutated (the kernel is storage-free).
        np.testing.assert_array_equal(psi_in[0], 4.0)
        np.testing.assert_array_equal(psi_in[1], 8.0)


# ─────────────────────────────────────────────────────────────────────
# L0 / regression: Bit-identity of the ordinate vectorisation
# ─────────────────────────────────────────────────────────────────────


def _per_ordinate_loop_reference(psi_in, s_axes, reaction_xs, Q_cells):
    """Per-ordinate Python-loop reference with the SAME left-fold order
    ``(Σ_t + 2g_0) + 2g_1`` / ``(Q + 2g_0·in_0) + 2g_1·in_1`` the kernel uses
    (the scheme's diamond ``2`` on each raw-``g`` streaming term — #240)."""
    N_oct = psi_in[0].shape[0]
    ref = np.empty_like(psi_in[0])
    for n in range(N_oct):
        Q_n = Q_cells[n if Q_cells.shape[0] > 1 else 0]
        denom = (reaction_xs + 2.0 * s_axes[0][n]) + 2.0 * s_axes[1][n]
        numer = (
            (Q_n + 2.0 * s_axes[0][n] * psi_in[0][n])
            + 2.0 * s_axes[1][n] * psi_in[1][n]
        )
        ref[n] = numer / denom
    return ref


@pytest.mark.l0
@pytest.mark.regression
class TestBitIdenticalToPerOrdinateLoop:
    r"""Per-element bit-equality: the batched kernel's vectorisation over the
    ordinate axis must reproduce the per-ordinate Python loop at IEEE-754
    ULP — the left-fold operation order is identical by construction."""

    def test_bit_identical_4ord_2g(self):
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=4, ng=2, n_diag=3, seed=12,
        )
        psi_avg, _ = DiamondDifference().cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        np.testing.assert_array_equal(
            psi_avg,
            _per_ordinate_loop_reference(psi_in, s_axes, reaction_xs, Q_cells),
        )

    def test_isotropic_Q_broadcasts_correctly(self):
        """Q with leading dim 1 (isotropic source) must broadcast cleanly."""
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=4, ng=2, n_diag=3, seed=21, Q_leading=1,
        )
        psi_avg, _ = DiamondDifference().cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        np.testing.assert_array_equal(
            psi_avg,
            _per_ordinate_loop_reference(psi_in, s_axes, reaction_xs, Q_cells),
        )


# ─────────────────────────────────────────────────────────────────────
# Wave O #208 O.4b — the APPLY kernel (closed form, affinity, round trip)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestResidualKernelClosedForm:
    r"""Single-cell closed-form: r = denom·ψ̄ − (Q + 2g_x·ψ_in_x + 2g_y·ψ_in_y),
    ``denom = Σ_t + 2g_x + 2g_y`` (the scheme's diamond ``2`` on raw ``g``, #240)."""

    def test_residual_at_solution_is_zero(self):
        r"""At the solved ψ̄ (= 120/18 from the solve closed-form test), the
        residual vanishes; the closure faces are reconstructed from the
        PROBE (so the matvec propagates edges exactly as the sweep does)."""
        psi_in = (np.full((1, 1, 1), 4.0), np.full((1, 1, 1), 8.0))
        s_axes = (np.full((1, 1, 1), 3.0), np.full((1, 1, 1), 5.0))  # raw g
        reaction_xs = np.full((1, 1), 2.0)
        Q_cells = np.full((1, 1, 1), 16.0)
        psi_bar = 120.0 / 18.0
        residual, psi_out = DiamondDifference().residual_kernel_batch(
            psi_bar=np.full((1, 1, 1), psi_bar),
            psi_in=psi_in, s_axes=s_axes,
            reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        # denom·ψ̄ − numer = 18·(120/18) − 120 = 120 − 120 = 0.
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)
        np.testing.assert_array_equal(psi_out[0], 2.0 * psi_bar - 4.0)
        np.testing.assert_array_equal(psi_out[1], 2.0 * psi_bar - 8.0)

    def test_residual_off_solution_is_affine(self):
        r"""A probe shifted by δ from the solution shifts the residual by
        denom·δ — the residual is linear in ψ̄."""
        psi_in = (np.full((1, 1, 1), 4.0), np.full((1, 1, 1), 8.0))
        s_axes = (np.full((1, 1, 1), 3.0), np.full((1, 1, 1), 5.0))  # raw g
        residual, _ = DiamondDifference().residual_kernel_batch(
            psi_bar=np.full((1, 1, 1), 7.0),       # above the solution 120/18
            psi_in=psi_in, s_axes=s_axes,
            reaction_xs=np.full((1, 1), 2.0),
            Q_cells=np.full((1, 1, 1), 16.0),
        )
        # denom·ψ̄ − numer = 18·7.0 − 120 = 126 − 120 = 6.0
        #   (== denom · δ = 18 · (7 − 120/18)).
        np.testing.assert_allclose(residual, 6.0, atol=1e-13)


@pytest.mark.l0
class TestKernelPairRoundTrip:
    r"""The solve↔apply contract: ``residual_kernel_batch`` at the value
    ``cell_kernel_batch`` solves for (same inputs) vanishes — the batched
    analogue of the per-cell ``residual``/``update`` round trip."""

    @pytest.mark.parametrize("seed", [77, 78, 79])
    def test_residual_vanishes_at_solve_solution(self, seed):
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=4, ng=2, n_diag=3, seed=seed,
        )
        psi_avg, psi_out_solve = DiamondDifference().cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        residual, psi_out_apply = DiamondDifference().residual_kernel_batch(
            psi_bar=psi_avg,
            psi_in=psi_in, s_axes=s_axes,
            reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)
        # Both directions reconstruct the SAME closure faces at the solution.
        for a in range(2):
            np.testing.assert_array_equal(psi_out_apply[a], psi_out_solve[a])

    def test_isotropic_Q_round_trip(self):
        """Round-trip holds with an isotropic-shaped (leading-1) Q too."""
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=4, ng=2, n_diag=3, seed=88, Q_leading=1,
        )
        psi_avg, _ = DiamondDifference().cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        residual, _ = DiamondDifference().residual_kernel_batch(
            psi_bar=psi_avg,
            psi_in=psi_in, s_axes=s_axes,
            reaction_xs=reaction_xs, Q_cells=Q_cells,
        )
        np.testing.assert_allclose(residual, 0.0, atol=1e-13)


# ─────────────────────────────────────────────────────────────────────
# The FP-reduction-tree-of-record pin (S6.4(e) gate-memo addendum)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.regression
class TestKernelSourceOfRecord:
    r"""sha256 source pin on the TWO kernel bodies — the left-fold order
    ``((Σ_t + s_0) + s_1) + …`` is bit-identity-LOAD-BEARING (IEEE-754
    addition is non-associative; every byte-identity anchor in the SN stack
    inherits from this fold order).

    This is the GENUINE source-hash exception (contrast the retired A2D-1,
    whose job an output oracle covered): the kernels ARE the FP reduction
    tree of record — an algebraically-equivalent rearrangement passes every
    value-tolerance test yet silently invalidates the 1-ULP regression
    contract.  Any deliberate edit updates the hash IN the same commit and
    re-baselines the bit-identity anchors per the Fork-B2 discipline.
    """

    # Updated #240 (s_axes=g): both kernels internalise the diamond 2 in the
    # left fold (``denom += 2.0*s_a``, ``numer += 2.0*s_a*in_a``); the fold is
    # still bit-identical to the legacy ``sigt + sx + sy`` because ``2*g_a``
    # equals the former pre-scaled ``2|μ|/Δ`` exactly (power-of-2 commute).
    # Re-hashed #240 Phase 2 D1: the inlined diamond reconstruction
    # ``2ψ̄ − ψ_in`` now routes through the generic affine outflow op.
    # Re-hashed #240 Phase 2 D2: that op was homed from the ``affine_closure``
    # module onto the scheme base, so the call is now
    # ``self.outgoing_face_from_average(ψ̄, ψ_in, 0.5)`` (the inherited
    # ``DiscretizationSchemeBase`` staticmethod).  A SOURCE change only; the ``w=½`` op is
    # byte-identical to ``2ψ̄ − ψ_in`` (``÷0.5`` is exact ``×2``), so every
    # numerical byte-identity anchor stays green.
    # Re-hashed #240 Phase 2 D5a (single-source): the ``Σ_t + Σ 2g_a`` left fold
    # + the per-axis ``2g_a`` couplings now come from the shared
    # ``DiamondDifference._cartesian_streaming_diagonal`` (the ONE home of the
    # diamond ``2 = 1/w_DD`` across the 3 Cartesian producers); the
    # upstream-numerator reuses ``couplings[a]`` (``c_a*in_a`` is byte-identical
    # to ``(2.0*s_a)*in_a``) and the outflow ``w`` is the named ``_DD_W``
    # (exactly ``0.5``).  A SOURCE change only — the FP reduction tree is
    # PRESERVED (left fold + reuse), so every numerical byte-identity anchor
    # (window≡full oracle, affine-carve golden, T4b cart2d/slab snapshot) stays
    # green; verified before re-hashing.
    # Re-hashed #241 (model-agnostic parameter rename): the kernel reaction-rate
    # parameter ``sigt_cells`` was renamed to the role name ``reaction_xs`` in
    # both kernel signatures + bodies (and the shared
    # ``_cartesian_streaming_diagonal`` it calls).  A SOURCE-TEXT change ONLY —
    # an identifier rename touches no operand, no operation, and no fold order,
    # so the FP reduction tree is bit-for-bit PRESERVED; every numerical
    # byte-identity anchor (window≡full oracle, affine-carve golden, DD strict
    # regression snapshot) stays green (verified before re-hashing).
    # Re-hashed sn-package-layout reorg (Phase 2): the ``orpheus.sn.sweep_graph``
    # cross-references in BOTH kernel docstrings were repointed to
    # ``orpheus.sn.loss_representation.sweep_graph`` (the module moved INTO the
    # ``loss_representation`` package).  A DOCSTRING-only change — no operand,
    # operation, or fold order touched, so the FP reduction tree is bit-for-bit
    # PRESERVED; every numerical byte-identity anchor (window≡full oracle,
    # affine-carve golden, T4b cart2d/slab snapshot) stays green (the Phase-2
    # behavioral subset confirmed it before re-hashing).
    # Re-hashed after the #231 P2-E prose rebalance (docs merge `ec74be50`):
    # both kernel DOCSTRINGS were trimmed by the code-prose rebalancing batch
    # (P2-E, docstring/comment-only with per-file code-token-stream invariance
    # proven at landing) — no operand, operation, or fold order touched.
    # Surfaced by the first full-tree run after that merge (the #310 C1 gate,
    # 2026-07-24): pre-docs-merge main `126c4e40` MATCHES the old hashes, and
    # the C1 commit's bodies hash IDENTICAL to pre-C1 `ec74be50` — the drift
    # is P2-E's docstring trim alone.  Every numerical byte-identity anchor
    # (window≡full oracle, affine-carve golden, DD strict regression, the
    # frozen walk_matvec adjoint baselines) was green in the SAME 6405-passed
    # run that exposed the stale pin.
    EXPECTED: ClassVar[dict[str, str]] = {
        "cell_kernel_batch":
            "b5d712fccb54c61f261c55661a59890566eeece24c14c8fd8eb67c4921f7e144",
        "residual_kernel_batch":
            "b9ffbecb9dfa5411a12855229baedfa01efaf1fae64aff82c4a093abd08387b7",
    }

    @pytest.mark.parametrize("kernel", ["cell_kernel_batch", "residual_kernel_batch"])
    def test_kernel_source_unchanged(self, kernel):
        import hashlib
        import inspect

        src = inspect.getsource(getattr(DiamondDifference, kernel))
        actual = hashlib.sha256(src.encode("utf-8")).hexdigest()
        if actual != self.EXPECTED[kernel]:
            pytest.fail(
                f"DiamondDifference.{kernel} source changed — this body is the "
                "FP reduction tree of record (the left-fold order is "
                "bit-identity-load-bearing).  If deliberate: update EXPECTED "
                f"to {actual} in THIS commit and re-verify every byte-identity "
                "anchor (window≡full oracles, affine-carve golden)."
            )


# ─────────────────────────────────────────────────────────────────────
# Retirement audit (S6.4(e)): the 4 direction×storage walks are GONE
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_graph_exposes_two_walks_not_four():
    """[L0 structural] The direction×storage product is retired: the graph
    carries the TWO storage walks (``walk_full`` / ``walk_windowed``), and
    none of the four retired direction-specific methods survives."""
    from orpheus.sn.loss_representation.sweep_graph import SweepDependencyGraph

    for retired in ("apply", "residual", "apply_windowed", "residual_windowed"):
        if hasattr(SweepDependencyGraph, retired):
            pytest.fail(
                f"SweepDependencyGraph.{retired} re-appeared — the "
                "direction×storage product (4 methods) was collapsed at "
                "S6.4(e) into walk_full/walk_windowed × the _CellSolve/"
                "_CellResidual level operations."
            )
    for walk in ("walk_full", "walk_windowed"):
        if not hasattr(SweepDependencyGraph, walk):
            pytest.fail(f"SweepDependencyGraph.{walk} missing.")


# ─────────────────────────────────────────────────────────────────────
# L0: NotImplementedError for strategies that don't override the pair
# ─────────────────────────────────────────────────────────────────────


class _NoKernelStrategy(DiscretizationSchemeBase, key="_no_kernel_strategy_test"):
    """Test stub: only overrides the per-cell ``update`` / ``residual``,
    not the batched kernel pair."""

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        return CellResult(cell_average_flux=source / total_xs)

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        del visit, upstream_state
        return total_xs * cell_avg - source


@pytest.mark.l0
class TestDefaultRaisesNotImplemented:
    """Strategies that don't override the kernel pair fail loudly."""

    def test_solve_kernel_default_raises(self):
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=1, ng=1, n_diag=1, seed=55,
        )
        with pytest.raises(NotImplementedError, match="cell_kernel_batch"):
            _NoKernelStrategy().cell_kernel_batch(
                psi_in=psi_in, s_axes=s_axes,
                reaction_xs=reaction_xs, Q_cells=Q_cells,
            )

    def test_apply_kernel_default_raises(self):
        psi_in, s_axes, reaction_xs, Q_cells = _random_kernel_operands(
            N_oct=1, ng=1, n_diag=1, seed=55,
        )
        with pytest.raises(NotImplementedError, match="residual_kernel_batch"):
            _NoKernelStrategy().residual_kernel_batch(
                psi_bar=np.zeros((1, 1, 1)),
                psi_in=psi_in, s_axes=s_axes,
                reaction_xs=reaction_xs, Q_cells=Q_cells,
            )
