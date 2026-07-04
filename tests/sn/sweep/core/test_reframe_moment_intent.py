r"""#246 — ``_reframe`` keys the moment-frame involution on INTENT, not on a
coincidental trailing-axis length (the S4-hazard negative control).

SKELETON (test-architect, PRE-IMPL for #246). The production fix is NOT written.
These tests pin the CONTRACT the fix must satisfy and will FAIL/ERROR until the
``is_moment_valued`` thread lands on ``_reframe`` (see the spec
``.claude/agent-memory/test-architect/issue_246_moment_axis_predicate_spec.md``).

The hazard
==========
``_reframe`` (``orpheus/sn/sweep_graph.py``) today guards with::

    if frame_signs is None or arr.shape[-1] != frame_signs.shape[0]:
        return arr
    return arr * frame_signs

The second clause keys the "does ``arr`` carry the ``2^d`` moment axis?" question
on a COINCIDENTAL trailing-axis length. On a d=2 mesh ``2^d == 4``; a genuine
NON-moment array whose trailing axis is coincidentally 4 collides with
``frame_signs.shape[0] == 4`` → the moment-frame sign-flip involution
(``octant_moment_frame_signs``) mis-fires on a non-moment array. Latent today (no
flat-source d≥2 LD path reaches ``_reframe`` with such an array), so the fix is
prophylactic: thread an ``is_moment_valued`` bool the caller sets from the
array's TYPED ORIGIN, keep the already-typed ``frame_signs is None ⟺ DD/Step``
half.

vv-principles: this is the negative control that proves the predicate keys on
intent, not coincidence (anti-pattern #11 positive/negative pair; Mode-6
convention drift at the test-design level).

``-O``-safe (Mode 8): ``np.testing.assert_*`` / ``pytest.fail`` only — NO bare
``assert``.

IMPLEMENTATION-SHAPE NOTE
=========================
This skeleton assumes the recommended signature ``_reframe(arr, frame_signs, *,
is_moment_valued)``. If the implementer instead expresses "non-moment" by
passing ``frame_signs=None`` at the caller (the redundant-with-``frame_signs is
not None`` route), retarget rows 1/2/4 to that semantics — the ASSERTIONS
(pass-through on non-moment, involution on moment, None short-circuit on DD) are
unchanged; only the call signature differs. Coordinate at impl time.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.spatial._ubld import octant_moment_frame_signs

# NOTE: import the production helper once the signature is settled. Until then
# this import pins the symbol under test.
from orpheus.sn.loss_representation.sweep_graph import _reframe


def _backward_d2_frame_signs() -> np.ndarray:
    r"""A ``(4,)`` involution for a per_axis=2, d=2 scheme with at least one -1.

    A BACKWARD octant (``(1, -1)``) is MANDATORY: an all-forward octant returns
    all-ones, ``arr * signs == arr``, and the collision is INVISIBLE (Mode-10
    activated-but-unconstrained). The premise guard below fails loudly if the
    involution does not actually flip a slot.
    """
    signs = octant_moment_frame_signs((1, -1), 2)
    if signs is None or -1.0 not in signs:
        pytest.fail(
            "test premise broken: frame_signs must carry a -1 so the involution "
            f"actually flips a moment slot (got {signs}) — otherwise the S4 "
            "collision is invisible (Mode-10)"
        )
    return signs


@pytest.mark.l0
def test_reframe_non_moment_4trailing_passthrough_when_not_moment_valued() -> None:
    r"""Row 1 (the centerpiece): a NON-moment array whose trailing axis is
    coincidentally length ``2^d == 4`` passes through UNCHANGED when
    ``is_moment_valued=False``.

    The OLD probe ``arr.shape[-1] != frame_signs.shape[0]`` is FALSE here (both
    are 4) → OLD would multiply by ``frame_signs`` and corrupt the non-moment
    array. The NEW typed-origin thread keeps it untouched.
    """
    frame_signs = _backward_d2_frame_signs()                  # (4,)
    ng = 2
    # A non-moment buffer (a flat/DD-shaped array) that HAPPENS to have a
    # trailing axis of length 4 — the d=2 2^d collision.
    arr = np.arange(ng * 4, dtype=float).reshape(ng, 4)
    out = _reframe(arr, frame_signs, is_moment_valued=False)  # NEW signature
    np.testing.assert_array_equal(
        out, arr,
        err_msg=(
            "non-moment array with trailing axis 4 was reframed — the predicate "
            "keyed on the coincidental axis length, not on is_moment_valued (S4)"
        ),
    )


@pytest.mark.l0
def test_reframe_real_moment_applies_involution_when_moment_valued() -> None:
    r"""Row 2 (positive half): a genuine moment array DOES get the involution
    when ``is_moment_valued=True`` — the fix is not a blanket 'never reframe'."""
    frame_signs = _backward_d2_frame_signs()
    ng = 2
    arr = np.arange(ng * 4, dtype=float).reshape(ng, 4)       # real 2^d moment
    out = _reframe(arr, frame_signs, is_moment_valued=True)
    np.testing.assert_array_equal(
        out, arr * frame_signs,
        err_msg="moment array was NOT reframed at a moment closure",
    )


@pytest.mark.l0
@pytest.mark.parametrize("is_moment_valued", [True, False])
def test_reframe_dd_step_passthrough_frame_signs_none(is_moment_valued) -> None:
    r"""Row 3 (backward-compat invariant): ``frame_signs is None`` (DD/Step)
    short-circuits REGARDLESS of ``is_moment_valued`` → byte-identical DD/Step."""
    arr = np.arange(2 * 4, dtype=float).reshape(2, 4)
    out = _reframe(arr, None, is_moment_valued=is_moment_valued)
    np.testing.assert_array_equal(
        out, arr,
        err_msg="DD/Step (frame_signs is None) must always pass through",
    )


@pytest.mark.l0
def test_reframe_old_probe_would_misfire_on_collision() -> None:
    r"""Row 4 (the documenting negative-control PAIR, anti-pattern #11): the OLD
    array-length probe WOULD have flipped signs on the collision, while the NEW
    thread does not. Demonstrates both the broken instance and the fix."""
    frame_signs = _backward_d2_frame_signs()
    arr = np.arange(2 * 4, dtype=float).reshape(2, 4)

    # Inline replica of the OLD guard logic.
    old_passes = (frame_signs is None) or (arr.shape[-1] != frame_signs.shape[0])
    out_old = arr if old_passes else arr * frame_signs
    if np.array_equal(out_old, arr):
        pytest.fail(
            "the OLD probe did NOT misfire on this input — the collision regime "
            "is not reproduced (pick a frame_signs with a -1 and matching length)"
        )

    out_new = _reframe(arr, frame_signs, is_moment_valued=False)
    np.testing.assert_array_equal(
        out_new, arr,
        err_msg="the NEW thread should pass the non-moment array through",
    )


@pytest.mark.l0
@pytest.mark.parametrize("d", [2, 3])
def test_reframe_keys_on_intent_for_d2_and_d3(d) -> None:
    r"""Blind-spot coverage (spec Gate-4 note 2): the predicate-keys-on-intent
    claim holds for the d=2 ``2^d==4`` AND the d=3 ``2^d==8`` collision. The d=3
    row is the only coverage 3-D LD gets for this until a 3-D consumer arrives.

    Builds a backward octant (all axes -1 → guaranteed a -1 in the involution),
    a non-moment array with the matching trailing length, and asserts
    pass-through under ``is_moment_valued=False``."""
    signs = octant_moment_frame_signs((-1,) * d, 2)           # (2^d,)
    if signs is None or -1.0 not in signs:
        pytest.fail(f"premise: frame_signs for d={d} must carry a -1 (got {signs})")
    width = 2 ** d
    arr = np.arange(2 * width, dtype=float).reshape(2, width)
    out = _reframe(arr, signs, is_moment_valued=False)
    np.testing.assert_array_equal(
        out, arr,
        err_msg=f"d={d}: non-moment array of trailing length {width} was reframed",
    )
