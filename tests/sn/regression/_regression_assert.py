"""Principled regression assertion + drift tripwire (single source of truth).

This module replaces the per-test magic tolerance floors (``rtol=1e-12``,
``atol=1e-13``, ``rtol=5e-6``) that ``test_dd_regression.py`` historically
carried with ONE principled correctness gate and ONE informational drift
tripwire. Coding-elegance Pattern 5 (build the right primitive, not the
right product): every regression case calls :func:`assert_regression`; no
case re-derives a tolerance.

Two layers, in priority order:

1. **Correctness gate (hard-fail).** The tolerance is tied to what the
   computation *promises*, not to a magic number:

   * ``kind="iterative"`` — k_eff / scalar flux produced by an outer
     power iteration or fixed-source source-iteration. The solver only
     converged to its own stopping criterion ``conv_tol`` (``keff_tol``
     for k_eff, ``flux_tol`` / ``inner_tol`` for flux). Asserting tighter
     than the solver promised is unphysical, so the gate is
     ``SAFETY × conv_tol``. ``conv_tol`` is read off the actual solver
     config used to produce the snapshot (see ``_generate_snapshots.py``
     ``RUN_CONFIG``) — it is NEVER hardcoded here.
   * ``kind="direct"`` — a single sweep with no outer iteration. The only
     drift possible is FP-non-associativity in the reduction tree, bounded
     by ``reduction_depth × ULP``. The gate is
     :func:`numpy.testing.assert_array_almost_equal_nulp` at
     ``nulp=reduction_depth``.

2. **Drift tripwire (informational).** If the value passed the correctness
   gate BUT moved beyond bit-identity (ULP distance > 0), emit a
   :class:`DriftWarning`. pytest summarises warnings at run end (deduped,
   visible) where a ``logging.debug`` line would be buried; the warning is
   *escalatable* — ``-W error::DriftWarning`` turns ANY sub-tolerance bit
   drift into a HARD fail, which is the strict bit-identity gate a
   "pure-refactor, zero-numerical-change" PR should run. ``logging.debug``
   carries the forensic per-element ULP breakdown for ``--log-cli-level``.

**``-O`` safety.** ``python -O`` strips bare ``assert`` statements in any
module that pytest's assertion-rewriter did not touch — and helper modules
imported by tests are NOT rewritten. Every correctness check in this module
therefore uses :mod:`numpy.testing` (which raises ``AssertionError``
unconditionally, independent of ``__debug__``) or an explicit ``raise``.
There is no bare ``assert`` anywhere below; the canonical ORPHEUS
invocation ``python -O -m pytest`` keeps the gate live. See
``vv-principles`` Mode 8.
"""
from __future__ import annotations

import logging
import warnings

import numpy as np


logger = logging.getLogger(__name__)


# A 10× headroom over the solver's own stopping criterion. The solver
# converged ``|x_k − x_{k−1}| < conv_tol`` (k_eff) or
# ``‖φ_k − φ_{k−1}‖/‖φ_k‖ < conv_tol`` (flux); the *gap to the true fixed
# point* is the increment amplified by ``ρ/(1−ρ)`` where ρ is the spectral
# radius of the iteration map. For the regression mixtures ρ ≲ 0.9 so the
# amplification ≲ 10; SAFETY=10 covers that amplification AND the
# cross-run FP-non-associativity jitter at the converged tail. It is a
# *principled* bound (one iteration-map amplification factor), not a fudge.
SAFETY = 10.0


class DriftWarning(UserWarning):
    """A regression value passed the correctness gate but is no longer bit-identical.

    Informational by default (pytest deduplicates + summarises it at run
    end). Escalate to a hard failure on a pure-refactor PR with
    ``-W error::DriftWarning`` — then any sub-tolerance bit drift fails
    the build, giving a strict bit-identity contract for changes that
    claim zero numerical effect.
    """


def _ulp_distance(actual: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """Per-element ULP distance between two float64 arrays.

    The ULP distance is the count of representable float64 values between
    ``actual`` and ``expected`` — the natural unit for FP-non-associativity
    drift. Implemented via the monotone integer ordering of the IEEE-754
    bit pattern (the standard "almost-equal-ulp" trick that
    :func:`numpy.testing.assert_array_almost_equal_nulp` uses internally).
    """
    a = np.asarray(actual, dtype=np.float64)
    e = np.asarray(expected, dtype=np.float64)
    ai = a.view(np.int64).copy()
    ei = e.view(np.int64).copy()
    # Map the two's-complement int view to a monotone ordering: negative
    # floats get reflected so that adjacent floats are adjacent ints.
    neg_a = ai < 0
    neg_e = ei < 0
    ai[neg_a] = np.int64(np.iinfo(np.int64).min) - ai[neg_a]
    ei[neg_e] = np.int64(np.iinfo(np.int64).min) - ei[neg_e]
    return np.abs(ai - ei).astype(np.int64)


def assert_regression(
    actual,
    expected,
    *,
    conv_tol: float,
    case_name: str,
    kind: str,
    reduction_depth: int | None = None,
    quantity: str = "scalar_flux",
) -> None:
    """Assert ``actual`` reproduces the frozen ``expected`` snapshot.

    Parameters
    ----------
    actual, expected
        Scalars or arrays. Coerced to float64.
    conv_tol
        The solver's own convergence stopping criterion for *this quantity*
        (``keff_tol`` for k_eff; ``flux_tol`` / ``inner_tol`` for flux).
        Required for ``kind="iterative"``; ignored for ``kind="direct"``.
        Read off the run config — NEVER hardcoded by the caller.
    case_name
        Snapshot case name, for the failure / drift message.
    kind
        ``"iterative"`` (outer-iterated result → ``SAFETY × conv_tol``) or
        ``"direct"`` (single sweep → ``nulp=reduction_depth``).
    reduction_depth
        FP-reduction depth for ``kind="direct"`` (≈ #cells or #ordinates
        summed). Required for ``kind="direct"``; ignored otherwise.
    quantity
        Human label ("k_eff" / "scalar_flux") for the messages.
    """
    actual_arr = np.asarray(actual, dtype=np.float64)
    expected_arr = np.asarray(expected, dtype=np.float64)

    if actual_arr.shape != expected_arr.shape:
        raise AssertionError(
            f"{case_name}: {quantity} shape regression — "
            f"actual {actual_arr.shape} != expected {expected_arr.shape}"
        )

    # ── Layer 1: principled correctness gate (hard-fail) ──────────────
    if kind == "iterative":
        if conv_tol is None:
            raise ValueError(
                f"{case_name}: kind='iterative' requires conv_tol "
                "(the solver's own stopping criterion for this quantity)."
            )
        tol = SAFETY * conv_tol
        # rtol bound on the iterative converged value. ``assert_allclose``
        # raises unconditionally (not gated on __debug__), so it survives -O.
        np.testing.assert_allclose(
            actual_arr, expected_arr, rtol=tol, atol=tol,
            err_msg=(
                f"{case_name}: {quantity} regression EXCEEDS the principled "
                f"correctness tolerance SAFETY×conv_tol = {SAFETY:g}×{conv_tol:.1e} "
                f"= {tol:.1e}. This is past the solver's own convergence floor — "
                f"investigate as a real drift, not FP noise."
            ),
        )
    elif kind == "direct":
        if reduction_depth is None or reduction_depth < 1:
            raise ValueError(
                f"{case_name}: kind='direct' requires reduction_depth >= 1 "
                "(the FP-reduction-tree depth bounding the ULP drift)."
            )
        # nulp = reduction_depth: a single sweep's only drift source is
        # FP-non-associativity over the reduction tree. assert_*_nulp
        # raises unconditionally, surviving -O.
        np.testing.assert_array_almost_equal_nulp(
            actual_arr, expected_arr, nulp=int(reduction_depth),
        )
    else:
        raise ValueError(
            f"{case_name}: unknown kind={kind!r}; expected "
            "'iterative' or 'direct'."
        )

    # ── Layer 2: drift tripwire (informational; escalatable) ──────────
    ulp = _ulp_distance(actual_arr, expected_arr)
    n_ulp = int(ulp.max()) if ulp.size else 0
    if n_ulp > 0:
        denom = np.maximum(np.abs(expected_arr), np.finfo(np.float64).tiny)
        rel = float(np.max(np.abs(actual_arr - expected_arr) / denom))
        tol_repr = (
            f"{SAFETY * conv_tol:.1e}" if kind == "iterative"
            else f"nulp={reduction_depth}"
        )
        warnings.warn(
            DriftWarning(
                f"{case_name}: {quantity} drifted {n_ulp} ULP / "
                f"{rel:.2e} rel (within tol {tol_repr})"
            ),
            stacklevel=2,
        )
        # Forensic layer — surfaced via --log-cli-level=DEBUG. Per-element
        # breakdown of WHERE the drift lives + the reduction-depth bound.
        if ulp.ndim == 0:
            logger.debug(
                "%s drift detail: scalar %r vs %r → %d ULP",
                case_name, float(actual_arr), float(expected_arr), n_ulp,
            )
        else:
            worst = np.unravel_index(int(np.argmax(ulp)), ulp.shape)
            logger.debug(
                "%s drift detail: max %d ULP at index %s "
                "(actual=%.17g expected=%.17g); n_drifted=%d/%d elements; "
                "kind=%s tol=%s",
                case_name, n_ulp, worst,
                float(actual_arr[worst]), float(expected_arr[worst]),
                int(np.count_nonzero(ulp)), ulp.size, kind, tol_repr,
            )
