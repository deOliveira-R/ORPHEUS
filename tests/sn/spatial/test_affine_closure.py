r"""Unit gate for the generic affine-scheme reconstruction ops (#240 Phase 2 D1).

:func:`~orpheus.sn.spatial.affine_closure.outgoing_face_from_average` is the
INVERSE of :func:`~orpheus.sn.spatial.affine_closure.cell_average`: given the
cell-average ``ψ̄`` and the upstream face ``ψ_in``, recover the downstream face
``ψ_out = (ψ̄ − (1−w)·ψ_in)/w``.  The two together form the affine cell-average
blend ``ψ̄ = (1−w)·ψ_in + w·ψ_out``.  This is the single-source reconstruction
that DD (``w=½``) and LD (``w=1/(1+k)``) both route through (D1 collapses the
inlined ``2ψ̄ − ψ_in`` / ``ψ̄ + (g/θ)(ψ̄ − ψ_in)/D₂`` duplicates).

The gate pins three facts:

1. **Exact inverse of** :func:`cell_average` — round-trip to FP tol for any
   ``w ∈ (0, 1]``.
2. **DD ``w=½`` byte-identity** — bit-for-bit equal to the inlined ``2ψ̄ − ψ_in``
   (``÷0.5`` is the exact power-of-2 ``×2``; round-to-nearest commutes with
   exact doubling).
3. **LD ``w=1/(1+k)`` algebraic equality** — equal (to FP tol, a DIFFERENT
   reduction tree) to the inlined ``ψ̄ + (g/θ)(ψ̄ − ψ_in)/D₂``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.affine_closure import (
    cell_average,
    outgoing_face_from_average,
)


@pytest.mark.foundation
def test_outgoing_face_is_exact_inverse_of_cell_average() -> None:
    r"""``outgoing_face_from_average(cell_average(in, out, w), in, w) == out``.

    The convex face blend and its inverse round-trip to FP noise for any
    ``w ∈ (0, 1]`` and any face pair.
    """
    rng = np.random.default_rng(240)
    face_in = rng.standard_normal((4, 3, 7))
    face_out = rng.standard_normal((4, 3, 7))
    # w spans the affine range: ½ (DD), 1 (Step/donor), and an LD-like blend.
    for w in (0.5, 1.0, 0.3, np.full((4, 3, 7), 0.42), rng.uniform(0.05, 1.0, (4, 3, 7))):
        psi_bar = cell_average(face_in, face_out, w)
        recon = outgoing_face_from_average(psi_bar, face_in, w)
        np.testing.assert_allclose(recon, face_out, rtol=1e-12, atol=1e-13)


@pytest.mark.foundation
def test_outgoing_face_dd_half_is_byte_identical() -> None:
    r"""DD ``w=½``: ``(ψ̄ − 0.5·ψ_in)/0.5`` is BIT-IDENTICAL to ``2ψ̄ − ψ_in``.

    ``÷0.5`` is the exact power-of-2 ``×2`` and round-to-nearest commutes with
    exact doubling, so the generic op matches the inlined DD reconstruction
    bit-for-bit — DD scan/matvec snapshots stay byte-identical under D1.
    """
    rng = np.random.default_rng(31415)
    psi_bar = rng.standard_normal((5, 4, 9))
    face_in = rng.standard_normal((5, 4, 9))
    inlined = 2.0 * psi_bar - face_in
    generic = outgoing_face_from_average(psi_bar, face_in, 0.5)
    if not np.array_equal(inlined, generic):  # Mode-8: fires under -O
        pytest.fail(
            "DD w=½ generic outflow is NOT byte-identical to 2ψ̄ − ψ_in; "
            f"max |Δ| = {np.max(np.abs(inlined - generic)):.3e}"
        )


@pytest.mark.foundation
def test_outgoing_face_ld_blend_matches_inlined_schur() -> None:
    r"""LD ``w=1/(1+k)``: generic op equals the inlined ``ψ̄ + (g/θ)(ψ̄−ψ_in)/D₂``.

    Algebraically identical (``k = (g/θ)/D₂``, ``w = 1/(1+k)``) but a DIFFERENT
    reduction tree → a principled ~1-ULP re-baseline, not byte-identity.
    """
    rng = np.random.default_rng(2718)
    psi_bar = rng.standard_normal((6, 2, 5))
    face_in = rng.standard_normal((6, 2, 5))
    g_over_theta = rng.uniform(0.01, 5.0, (6, 2, 5))
    sigt = rng.uniform(0.01, 5.0, (6, 2, 5))
    d2 = g_over_theta + sigt
    inlined = psi_bar + g_over_theta * (psi_bar - face_in) / d2
    w = 1.0 / (1.0 + g_over_theta / d2)
    generic = outgoing_face_from_average(psi_bar, face_in, w)
    np.testing.assert_allclose(generic, inlined, rtol=1e-12, atol=1e-13)
