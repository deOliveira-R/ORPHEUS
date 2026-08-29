r"""``ordinate_scan`` conditioning-fallback regression catcher (ERR-054 + ERR-057).

GitHub issues: `#209 <https://github.com/deOliveira-R/ORPHEUS/issues/209>`_
(ERR-054 — the exact-reset ``a = 0`` cause) and
`#222 <https://github.com/deOliveira-R/ORPHEUS/issues/222>`_
(ERR-057 — the denormal-cumprod underflow cause).

Promoted from the numerics-investigator diagnostic cascade
(``derivations/diagnostics/diag_si_cyl_20cell_nan_step1_characterize.py``
+ ``..._step5_root_cause.py``) per the ERR-054 catalog entry.  The
``_step5_root_cause`` half was **retired 2026-08-09** (#347): the
division-free backend below inverted its central "the scan NaNs"
assertion — exactly the retirement trigger its own docstring named —
and its ``μ_x = 1/√20`` fixture no longer exists in any shipped rule.
This file and ``test_si_cyl_20cell_nan_regression.py`` are its
successors.  This one pins the load-bearing contract for the bug class

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

ERR-057 (`#222 <https://github.com/deOliveira-R/ORPHEUS/issues/222>`_) is
the SIBLING conditioning gap of the same closed-form scan, with a
DIFFERENT cause: a long optically-thick chain decays ``cumprod_a`` into
the IEEE-754 denormals (``~1e-312``) WITHOUT reaching an exact zero, so
``b / cumprod_a`` overflows to ``inf`` and ``cumprod_a · inf = NaN``.
The original ``cumprod_a[-1] == 0`` guard tested a PROXY for the
exact-reset cause and silently MISSED this denormal band — a
guard-predicate-incompleteness bug.  The fix dispatches on the TRUE
failure condition (a non-finite closed form, which catches BOTH causes);
``TestOrdinateScanDenormalUnderflow`` pins it against the explicit serial
loop.

Two structurally-independent grounds (``vv-principles`` §1/§6):

* **Scan-form contract** (``TestOrdinateScanReset``).  The reference is
  the explicit Python recurrence loop — a structurally-independent
  ground for the closed-form scan: a different algorithm (sequential
  fold, no division) for the same recurrence.  This catches the bug
  WITHOUT any solver in the loop.
* **Solver-level tier** (``TestSICylinderResonance``).  Since the 6.3
  flip the admitted cylinder family is the σ_y fold, where the
  pole-cell resonance is UNREACHABLE at physical Σ_t (``[M]``
  2026-08-08 — the per-ordinate resonant ``Σ_t*`` is ≤ 0 everywhere);
  the class now pins that unreachability as the 6.4 tripwire (the
  absorber retirement changes ``c_out``) plus a plain SI-vs-Krylov
  agreement row.  The closed-form ``k_inf`` pin lives in the ``@slow``
  solver gate
  ``tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py``
  (its own fixture migrates with the curvilinear leg).

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
from orpheus.sn.sweep.scan import ordinate_scan
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
    def test_ordinate_scan_multiple_and_consecutive_resets(self) -> None:
        r"""Multiple, consecutive, chain-leading, per-lane resets all finite.

        Subsumes the single-reset case (the solver-level gate
        ``test_si_cyl_20cell_nan_regression`` also pins a single trailing
        zero against the explicit loop). Stresses the reset handling
        beyond a single zero: lane 0 has two *consecutive* resets
        (``(0, b)`` annihilating ``(0, b')`` to its left under ``⊕``);
        lane 1 has one interior reset; lane 2 has a reset at the chain
        *start* and one at the *end*.  Each lane must match the serial
        loop to machine precision.
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


@pytest.mark.l0
class TestOrdinateScanDenormalUnderflow:
    r"""``ordinate_scan`` must stay finite in the denormal-cumprod band.

    ERR-057 / Issue #222.  This is the SIBLING of the ``a = 0`` reset
    (ERR-054): the chain has NO exact reset (``|a| < 1`` throughout), but
    a long optically-thick run decays ``cumprod_a`` into the IEEE-754
    denormals (``~1e-312``) so the Blelloch division ``b / cumprod_a``
    overflows to ``inf`` and ``cumprod_a · inf = NaN``.  The old
    ``cumprod_a[-1] == 0`` proxy guard caught the exact-reset cause but
    MISSED this denormal band — a guard-predicate-incompleteness bug.
    The finite-closed-form dispatch routes it to the division-free
    pair-monoid; the reference is the explicit serial loop (no division),
    which stays finite throughout.

    Pre-fix this FAILS (the closed form leaks a NaN); the assertions use
    ``pytest.fail`` (not bare ``assert``) so they fire under ``python -O``
    (``vv-principles`` Mode 8).
    """

    @pytest.mark.catches("ERR-057")
    def test_denormal_cumprod_underflow_stays_finite(self) -> None:
        r"""A denormal-band cumprod (no exact zero) stays finite + exact."""
        nx = 312
        a = np.full(nx, 0.1)
        b = np.full(nx, 1.0)
        psi_0 = 1.0

        # Precondition: the chain must land IN the denormal band — a
        # nonzero cumprod below `tiny` — or the old exact-zero guard would
        # already have caught it and this regime would be inert.  Kept as
        # a function call so it fires under -O (Mode 8).
        cumprod_tail = float(np.cumprod(a)[-1])
        if not 0.0 < cumprod_tail < np.finfo(float).tiny:
            pytest.fail(
                f"chain mis-tuned: cumprod[-1]={cumprod_tail:.3e} is not in "
                "the denormal band (0, tiny); the ERR-057 regime is inert"
            )

        result = ordinate_scan(a, b, psi_0)
        reference = _explicit_scan_loop(a, b, psi_0)

        if not np.all(np.isfinite(result)):
            n_bad = int(np.sum(~np.isfinite(result)))
            pytest.fail(
                "ordinate_scan leaked a non-finite value in the denormal-"
                f"cumprod band (ERR-057, {n_bad} non-finite): the closed-form "
                "guard missed the underflow"
            )
        np.testing.assert_allclose(result, reference, rtol=1e-12, atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# 2.  Solver-level eigenvalue pin (closed-form k_inf ground)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
class TestSICylinderResonance:
    r"""Issue #209's pole-cell resonance — UNREACHABLE on the admitted family.

    History in three states. (1) The pre-#337 fixture FOUND the
    ``a = 0`` resonance at a numerical coincidence: the OLD
    project-convention LS8 seed carried ``μ_x = −1/√20``, which solved
    the pole-cell identity ``2|μ|·A_total = dA_w·c_out + Σ_t·V``
    (ERR-054) bit-exactly at ``(thick=2, n=20, Σ_t=1)``.  (2) #337's
    moment-matched seed retired that node (its μ₁² is 1/21) — the old
    ordinate-present tripwire fired exactly as designed.  (3) At the
    6.3 flip the admitted cylinder family is the σ_y fold, and there
    the resonance is UNREACHABLE at physical cross sections: ``[M]``
    (2026-08-08) the per-ordinate resonant ``Σ_t*`` — from the affine
    ``1/(a+1)`` in ``Σ_t`` — is ≤ 0 for EVERY inward ordinate at
    every probed ``(n, R)``, i.e. ``2|μ|·A_total ≤ dA_w·c_out`` on
    every clamped folded arc, so ``a < 0`` for all ``Σ_t > 0``.

    The unreachability test IS the tripwire: the [½,1] absorber's
    retirement (6.4) changes ``c_out`` — if that makes the resonance
    reachable again, the gate reds and 6.4 owes a live reproducer.
    The ERR-054 class itself stays pinned by the quadrature-free
    scan-form contract test above (``a[·] = 0`` chains fed directly);
    the solver-level ``catches`` marker moved there with the
    reachability.
    """

    @staticmethod
    def _build():
        mix = get_mixture("A", "2g")
        mesh = curvilinear_homogeneous_mesh(
            20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL,
        )
        quad = Quadrature.folded_product(n_mu=8, n_phi=16)
        return {0: mix}, mesh, quad

    def test_pole_resonance_unreachable_on_admitted_family(self) -> None:
        r"""``Σ_t* ≤ 0`` for every inward ordinate — the ``a = 0``
        pole-cell resonance cannot be assembled at physical Σ_t on the
        folded rule.  REDS if a closure change (the 6.4 absorber
        retirement) makes it reachable — then a live reproducer is
        owed here again."""
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.sn.sweep.cache import (
            CollisionCache,
            StreamingCoefficientCache,
        )

        materials, mesh, quad = self._build()
        probe = SNMesh(mesh, quad, materials)
        geom = StreamingCoefficientCache.from_mesh_and_quad(probe, probe.angular_closure)
        mu = np.asarray(quad.mu_x)
        inward = np.flatnonzero(mu < 0)

        def _pole_a(ord_i: int, sig1: float) -> float:
            sig_t = np.ones((probe.ng, probe.nx))
            sig_t[0, :] = sig1
            cache = CollisionCache.from_geometry(
                geom, sig_t, probe.scheme,
            )
            return float(cache.a_attenuation[ord_i, 0, -1])

        worst = -np.inf
        for ord_i in map(int, inward):
            f1 = 1.0 / (_pole_a(ord_i, 1.0) + 1.0)
            f2 = 1.0 / (_pole_a(ord_i, 2.0) + 1.0)
            worst = max(worst, 1.0 + (1.0 - f1) / (f2 - f1))
        assert worst <= 1e-9, (
            f"an inward ordinate admits a PHYSICAL resonant sigma_t* = "
            f"{worst:.6g} on the admitted folded rule — the ERR-054 "
            f"pole-cell resonance is reachable again (absorber "
            f"retirement?); restore a live solver-level reproducer here"
        )

    def test_si_agrees_with_krylov(self) -> None:
        r"""SI and Krylov agree on the folded cylinder to tight
        tolerance — Krylov is the structurally-independent solver
        ground (the matvec ``apply`` path never consumes
        ``ordinate_scan``)."""
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
