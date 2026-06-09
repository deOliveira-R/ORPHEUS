r"""ERR-054 regression catcher — ``ordinate_scan`` chain-reset (``a = 0``).

GitHub issue: `#209 <https://github.com/deOliveira-R/ORPHEUS/issues/209>`_.

Promoted from the numerics-investigator diagnostic cascade
(``derivations/diagnostics/diag_si_cyl_20cell_nan_step1_characterize.py``
+ ``..._step5_root_cause.py``) per the ERR-054 catalog entry.  This
file pins the load-bearing contract for the bug class

    "Blelloch closed-form first-order linear-recurrence scan divides by
    the cumulative product ``cumprod_a``; an exact chain reset
    ``a[i] = 0`` collapses the cumprod to zero, the division produces
    ``Inf``, ``cumsum`` propagates NaN, and ``0 · NaN = NaN``."

The math is finite at the same chain.  At ``a[i] = 0`` the recurrence
``ψ[i+1] = a[i]·ψ[i] + b[i]`` degenerates to ``ψ[i+1] = b[i]`` — the
chain *forgets its history*, a fully-defined "fully attenuated" exit
flux.  The pole cell of a 1-D cylindrical mesh hits this exactly: the
inner radial face has zero area, so ``A_down = 0`` and the cache's
``a = 2|μ|·A_total / (dA_w·c_out + Σ_t·V) − 1`` lands bit-exactly on 0
at the ``(thick=2, n=20, μ_x=-1/√20, Σ_t=1)`` resonance.

Two structurally-independent grounds (``vv-principles`` §1/§6):

* **Scan-form contract** (``TestOrdinateScanReset``).  The reference is
  the explicit Python recurrence loop — a structurally-independent
  ground for the closed-form scan: a different algorithm (sequential
  fold, no division) for the same recurrence.  This catches the bug
  WITHOUT any solver in the loop.
* **Solver-level eigenvalue pin** (``TestSICylinderResonance``).  The
  reference is the analytical infinite-medium eigenvalue
  ``k_inf = νΣ_f / Σ_a = 1.875`` for mixture A 2G — a closed-form
  reference (pillar: closed-form).  NaN is the one failure mode the
  ``k = νΣ_f/Σ_a`` flux-shape invariance cannot smuggle past the
  eigenvalue identity, so it propagates to ``keff = NaN`` pre-fix.

The fast-path bit-identity gate (``test_fast_path_bit_identical``)
guards the other half of the fix contract: when the chain has NO
reset, ``ordinate_scan`` must run the UNCHANGED Blelloch closed form,
bit-for-bit, so the hot SN sweep path neither regresses numerically
nor in performance.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from orpheus.sn.spatial.scan import ordinate_scan
from tests.sn._test_helpers import curvilinear_homogeneous_mesh

pytestmark = [pytest.mark.verifies("blelloch-1990-eq-1-5")]


def _explicit_scan_loop(
    a: np.ndarray, b: np.ndarray, psi_0: np.ndarray,
) -> np.ndarray:
    r"""Sequential reference for ``ψ[i+1] = a[i]·ψ[i] + b[i]``.

    Structurally independent of :func:`ordinate_scan`: a serial fold
    with no division, exact at ``a = 0`` (the in-cell update reduces to
    ``ψ = b``).  ``out[i]`` is the state AFTER cell ``i``.
    """
    nx = a.shape[0]
    out = np.empty_like(a)
    psi = (
        np.broadcast_to(psi_0, a.shape[1:]).astype(float).copy()
        if a.ndim > 1
        else np.asarray(psi_0, dtype=float)
    )
    for i in range(nx):
        psi = a[i] * psi + b[i]
        out[i] = psi
    return out


def _blelloch_closed_form(
    a: np.ndarray, b: np.ndarray, psi_0: np.ndarray,
) -> np.ndarray:
    r"""The legacy division-based Blelloch form (no-reset fast path).

    Reference for the bit-identity gate: the exact expression that
    :func:`ordinate_scan` MUST still execute when the chain has no
    reset.
    """
    cumprod_a = np.cumprod(a, axis=0)
    return cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))


# ═══════════════════════════════════════════════════════════════════════
# 1.  Scan-form contract (no solver — structurally independent ground)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
class TestOrdinateScanReset:
    r"""``ordinate_scan`` must be finite and exact at every ``a = 0``.

    The reference is the explicit recurrence loop — a different
    algorithm (serial fold, no division) for the same math.  Pre-fix
    these all FAIL (the Blelloch division NaNs from the reset onward).
    """

    @pytest.mark.catches("ERR-054")
    def test_ordinate_scan_at_a_zero_returns_finite_via_loop(self) -> None:
        r"""A single exact reset ``a[i] = 0`` gives the loop's value, finite.

        At the reset the chain restarts: ``ψ[i+1] = b[i]``, independent
        of all prior history.  Entries AFTER the reset must NOT inherit
        NaN.  Bit-exact against the explicit loop.
        """
        a = np.array([0.5, 0.0, 0.5, 2.0])
        b = np.array([1.0, 3.0, 0.5, 0.1])
        psi_0 = 0.7

        result = ordinate_scan(a, b, psi_0)
        reference = _explicit_scan_loop(a, b, psi_0)

        assert np.all(np.isfinite(result)), (
            f"ordinate_scan returned non-finite at a=0: {result}"
        )
        # The reset cell forgets history: ψ[1] = b[1] = 3.0 exactly.
        assert result[1] == b[1], (
            f"reset cell ψ[1] = {result[1]!r}, expected b[1] = {b[1]!r}"
        )
        np.testing.assert_array_equal(result, reference)

    @pytest.mark.catches("ERR-054")
    def test_ordinate_scan_multiple_and_consecutive_resets(self) -> None:
        r"""Multiple, consecutive, chain-leading, per-lane resets all finite.

        Stresses the reset handling beyond a single zero: lane 0 has two
        *consecutive* resets (``(0, b)`` annihilating ``(0, b')`` to its
        left under ``⊕``); lane 1 has one interior reset; lane 2 has a
        reset at the chain *start* and one at the *end*.  Each lane must
        match the serial loop to machine precision.
        """
        rng = np.random.default_rng(seed=209)
        a = rng.uniform(0.3, 1.4, size=(8, 3))
        b = rng.uniform(-0.4, 0.7, size=(8, 3))
        a[1, 0] = 0.0
        a[2, 0] = 0.0  # consecutive reset, lane 0
        a[5, 1] = 0.0  # interior reset, lane 1
        a[0, 2] = 0.0  # chain-leading reset, lane 2
        a[7, 2] = 0.0  # chain-trailing reset, lane 2
        psi_0 = np.array([0.7, 0.3, 1.1])

        result = ordinate_scan(a, b, psi_0)
        reference = _explicit_scan_loop(a, b, psi_0)

        assert np.all(np.isfinite(result)), (
            "ordinate_scan non-finite with multiple/consecutive resets"
        )
        np.testing.assert_allclose(result, reference, rtol=1e-13, atol=1e-15)

    def test_fast_path_bit_identical(self) -> None:
        r"""NO reset ⇒ the UNCHANGED Blelloch closed form, bit-for-bit.

        The other half of the fix contract (``vv-principles``
        §"Bit-identity vs principled-equivalence"): the hot SN sweep
        path takes a reset only at the rare pole-cell resonance, so the
        common no-reset path MUST be the legacy expression byte-for-byte
        — no FP-reduction-tree change, no performance regression.
        Parametrised over scalar / multi-group / multi-ordinate shapes.
        """
        rng = np.random.default_rng(seed=11)
        for shape in [(30,), (30, 4), (40, 8, 3)]:
            a = rng.uniform(0.3, 1.4, size=shape)
            b = rng.uniform(-0.4, 0.7, size=shape)
            psi_0 = (
                rng.uniform(0.1, 0.5, size=shape[1:])
                if len(shape) > 1
                else 0.3
            )
            new = ordinate_scan(a, b, psi_0)
            legacy = _blelloch_closed_form(a, b, psi_0)
            np.testing.assert_array_equal(
                new, legacy,
                err_msg=f"fast path not bit-identical at shape {shape}",
            )


# ═══════════════════════════════════════════════════════════════════════
# 2.  Solver-level eigenvalue pin (closed-form k_inf ground)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
class TestSICylinderResonance:
    r"""The Issue #209 minimal reproducer at the exact resonance.

    Cylinder (thick=2, n=20, reflective) × LevelSymmetric S8 × mixture
    A 2G.  The smallest-|μ| ordinate ``μ_x = -1/√20`` hits the pole-cell
    ``a = 0`` resonance bit-exactly.  Pre-fix ``inner_solver=
    "source_iteration"`` returns ``keff = NaN``; ``"krylov"`` (which
    never touches ``ordinate_scan``) returns the analytical ``k_inf``.
    """

    _K_INF = 1.875  # νΣ_f / Σ_a for mixture A, 2 groups (closed form).

    @staticmethod
    def _build():
        mix = get_mixture("A", "2g")
        mesh = curvilinear_homogeneous_mesh(
            20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL,
        )
        quad = Quadrature.level_symmetric(8)
        return {0: mix}, mesh, quad

    def test_resonant_ordinate_present(self) -> None:
        r"""``μ_x = -1/√20`` is in S8 — the resonance is reachable."""
        _, _, quad = self._build()
        assert np.any(np.isclose(quad.mu_x, -1.0 / np.sqrt(20))), (
            "LevelSymmetric S8 must carry the resonant ordinate "
            "μ_x = -1/√20; without it the reproducer is inert"
        )

    @pytest.mark.catches("ERR-054")
    def test_si_agrees_with_kinf_at_resonance(self) -> None:
        r"""SI returns the analytical ``k_inf``, finite — not NaN.

        THE L1 correctness pin.  Pre-fix this FAILS: SI returns
        ``keff = NaN`` because ``ordinate_scan`` NaNs at the pole-cell
        reset and the NaN propagates through every source iteration.
        """
        materials, mesh, quad = self._build()
        result = solve_sn(
            materials, mesh, quad,
            inner_solver="source_iteration",
            max_inner=500, inner_tol=1e-10,
        )
        assert np.isfinite(result.keff), (
            f"SI keff = {result.keff!r} (NaN ⇒ ERR-054 regressed: "
            "ordinate_scan divides by a zero cumprod at the pole cell)"
        )
        assert abs(result.keff - self._K_INF) < 1e-6, (
            f"SI keff = {result.keff:.10f} vs k_inf = {self._K_INF}"
        )

    @pytest.mark.catches("ERR-054")
    def test_si_agrees_with_krylov_at_resonance(self) -> None:
        r"""SI and Krylov agree at the resonance to tight tolerance.

        Krylov is the structurally-independent solver ground (the
        matvec ``apply`` path never consumes ``ordinate_scan``).  Both
        paths must converge to the same eigenvalue.
        """
        materials, mesh, quad = self._build()
        si = solve_sn(
            materials, mesh, quad, inner_solver="source_iteration",
            max_inner=500, inner_tol=1e-10,
        )
        kr = solve_sn(
            materials, mesh, quad, inner_solver="krylov",
            max_inner=500, inner_tol=1e-10,
        )
        assert np.isfinite(si.keff) and np.isfinite(kr.keff)
        assert abs(si.keff - kr.keff) < 1e-7, (
            f"SI keff = {si.keff:.10f} vs Krylov keff = {kr.keff:.10f}"
        )
