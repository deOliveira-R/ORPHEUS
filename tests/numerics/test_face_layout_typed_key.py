r"""The typed-key ``FaceLayout[K]`` — the ψ½ space's ``(level, sign, part)`` layout.

C2a of the codim-1 ``FaceField`` carve migrated
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
off its hand-rolled ``_leg_offset``/``cells_slice``/``corner_slice`` arithmetic
onto a real :class:`~orpheus.numerics.face_layout.FaceLayout`, keyed by the
structured ``(level, sign, part)`` tuple instead of the spatial trace's ``str``
face names. The offsets are *identical* — the point of the carve is a single
flat-buffer discipline shared across key regimes, not a behavior change.

**Why this gate exists (test-architect deliverable-1 gap).** The production
sphere-GL mesh carries exactly **one** seed-carrying μ-level, so every existing
pin exercises the layout at ``n_levels == 1``. At a single level the flat offset
of leg ``(pos, sign)`` is ``sign_pos · per_leg`` — the level-stride term
``2 · pos · per_leg`` is identically zero (``pos == 0`` always) and therefore
**unconstrained**. A typed-key regression in the level-index term (a wrong
emission order, a dropped ``2·`` stride, a re-hand-rolled offset) is invisible to
every single-level test. This module pins the layout at a **multi-level**
synthetic (``levels = (0, 2, 5)``), where ``pos`` runs ``0, 1, 2`` and the stride
term is live — and mutation-proves the golden has teeth.

Pure numerics — no mesh: :meth:`RadialCharacteristicSpace.for_levels` needs only
``levels``, ``ng``, ``nx``, ``cell_volumes``. All assertions are
``np.testing.assert_*`` / :func:`pytest.fail` (function calls — they fire under
the canonical ``python -O``, unlike bare ``assert``; vv Mode 8).
"""

from __future__ import annotations

import numpy as np
import pytest

import orpheus.numerics.spaces.radial_characteristic_space as sds_mod
from orpheus.numerics.spaces.radial_characteristic_space import RadialCharacteristicSpace

# The canonical layout conventions, hard-coded HERE (independent of the
# production module) so a production emission-order mutation diverges from the
# golden rather than silently tracking it.
_SIGN_ORDER = (-1, +1)  # sign_pos: -1 → 0 (inward, first), +1 → 1 (outward)
_PARTS = ("cells", "corner")  # cells first, corner second, within each leg


def _closed_form_slot(
    levels: tuple[int, ...], ng: int, nx: int, level: int, sign: int, part: str,
) -> tuple[int, tuple[int, ...], int]:
    r"""The ``(offset, shape, flat_size)`` a correct flat layout MUST assign.

    Derived from first principles, NOT from the production code:

    * one leg = cells ``(ng, nx)`` then corner ``(ng,)`` ⇒ ``per_leg = ng·nx + ng``;
    * legs walk ``(level ascending, sign ∈ (-1, +1))`` ⇒ leg index
      ``2·pos + sign_pos`` (``pos`` the level's rank, ``sign_pos`` from
      :data:`_SIGN_ORDER`);
    * cells at the leg base, corner at ``base + ng·nx``.
    """
    per_leg = ng * nx + ng
    pos = levels.index(level)
    sign_pos = _SIGN_ORDER.index(sign)
    base = (2 * pos + sign_pos) * per_leg
    if part == "cells":
        return base, (ng, nx), ng * nx
    return base + ng * nx, (ng,), ng


def _assert_layout_matches_closed_form(
    space: RadialCharacteristicSpace, levels: tuple[int, ...], ng: int, nx: int,
) -> None:
    r"""Assert every ``(level, sign, part)`` slot matches the closed form.

    Raises :class:`AssertionError` (via ``np.testing``) on the first mismatch —
    so callers can both USE it as a gate and CATCH it to prove teeth.
    """
    for level in levels:
        for sign in _SIGN_ORDER:
            for part in _PARTS:
                key = (level, sign, part)
                slot = space.layout.faces[key]
                exp_offset, exp_shape, exp_flat = _closed_form_slot(
                    levels, ng, nx, level, sign, part,
                )
                np.testing.assert_equal(
                    slot.offset, exp_offset,
                    err_msg=f"offset mismatch at {key}",
                )
                np.testing.assert_equal(
                    slot.shape, exp_shape, err_msg=f"shape mismatch at {key}",
                )
                np.testing.assert_equal(
                    slot.flat_size, exp_flat,
                    err_msg=f"flat_size mismatch at {key}",
                )


# ── The multi-level fixture ──────────────────────────────────────────────
_LEVELS = (0, 2, 5)
_NG, _NX = 2, 4
_CELL_VOLUMES = np.array([1.5, 2.5, 3.5, 4.5])  # shape (nx,), SPD


def _multilevel_space() -> RadialCharacteristicSpace:
    return RadialCharacteristicSpace.for_levels(
        _LEVELS, ng=_NG, nx=_NX, cell_volumes=_CELL_VOLUMES,
    )


# ── The pin ──────────────────────────────────────────────────────────────


def test_multilevel_offsets_match_closed_form():
    r"""Every ``(level, sign, part)`` slot lands at the closed-form offset.

    THE gate the single-level production pins cannot provide: at ``levels =
    (0, 2, 5)`` the level-stride term ``2·pos·per_leg`` is exercised for
    ``pos ∈ {0, 1, 2}`` (offsets 0, 20, 40 for the inward-cells slots — stride
    ``2·per_leg = 20``).
    """
    space = _multilevel_space()
    _assert_layout_matches_closed_form(space, _LEVELS, _NG, _NX)
    # total_size = n_legs · per_leg = (3 levels · 2 signs) · (ng·nx + ng).
    np.testing.assert_equal(space.layout.total_size, 3 * 2 * (_NG * _NX + _NG))
    np.testing.assert_equal(space.shape, (3 * 2 * (_NG * _NX + _NG),))


def test_multilevel_metric_weights_correspond_to_the_layout_slots():
    r"""Each slot's G_sd = V_cell weights sit at the slot's flat region.

    The layout offsets AND the state-metric weights are emitted in ONE walk
    over the legs (C2c single-source), so a slot's flat region and its V_cell
    weights must correspond: cells → ``V_cell`` per group; corner → ``V(r=R)``.
    A regression that re-derived the metric order independently of the layout
    walk (the retired parallel ``np.tile``) would break this correspondence
    without moving any offset — so the offset golden alone cannot catch it.
    """
    space = _multilevel_space()
    weights = np.asarray(space.inner_product_weights)
    for level in _LEVELS:
        for sign in _SIGN_ORDER:
            cells = space.layout.faces[(level, sign, "cells")]
            corner = space.layout.faces[(level, sign, "corner")]
            cells_w = weights[cells.offset : cells.offset + cells.flat_size]
            corner_w = weights[corner.offset : corner.offset + corner.flat_size]
            # cells: V_cell per group (ng tiled copies of the nx-vector).
            np.testing.assert_array_equal(cells_w, np.tile(_CELL_VOLUMES, _NG))
            # corner: the r = R gauge weight V[-1] on every group.
            np.testing.assert_array_equal(
                corner_w, np.full(_NG, float(_CELL_VOLUMES[-1])),
            )


def test_multilevel_views_roundtrip_without_aliasing():
    r"""Per-leg writes through ``cells_view``/``corner_view`` do not alias.

    The behavioral consequence of correct, non-overlapping offsets: a distinct
    sentinel written into each of the 3 · 2 legs (cells and corner) reads back
    exactly, with no cross-leg bleed. A single-level layout cannot exercise the
    inter-level stride this catches.
    """
    space = _multilevel_space()
    buf = np.zeros(space.shape)
    # Write a unique sentinel per (level, sign) leg.
    sentinels: dict[tuple[int, int], tuple[float, float]] = {}
    tag = 1.0
    for level in _LEVELS:
        for sign in _SIGN_ORDER:
            cell_val, corner_val = tag, tag + 0.5
            space.cells_view(buf, level, sign)[:] = cell_val
            space.corner_view(buf, level, sign)[:] = corner_val
            sentinels[(level, sign)] = (cell_val, corner_val)
            tag += 10.0
    # Read back — each leg must recover EXACTLY its own sentinel.
    for level in _LEVELS:
        for sign in _SIGN_ORDER:
            cell_val, corner_val = sentinels[(level, sign)]
            np.testing.assert_array_equal(
                space.cells_view(buf, level, sign), np.full((_NG, _NX), cell_val),
            )
            np.testing.assert_array_equal(
                space.corner_view(buf, level, sign), np.full((_NG,), corner_val),
            )


# ── Why multi-level: the stride term is single-level-blind ───────────────


def test_level_stride_term_is_single_level_blind():
    r"""``2·pos`` vs ``pos`` coincide at one level, diverge from ``pos = 1``.

    A pure-arithmetic demonstration (no production call) that the level-index
    coefficient is unobservable at ``n_levels == 1`` — the reason a single-level
    fixture is an insufficient pin and this module builds a 3-level one.
    """
    per_leg = _NG * _NX + _NG

    def leg_offset(pos: int, sign_pos: int, *, correct: bool) -> int:
        stride = (2 * pos) if correct else pos  # the mutant drops the factor 2
        return (stride + sign_pos) * per_leg

    # Single level (pos = 0): correct and mutant are IDENTICAL for both signs.
    for sign_pos in (0, 1):
        np.testing.assert_equal(
            leg_offset(0, sign_pos, correct=True),
            leg_offset(0, sign_pos, correct=False),
        )
    # From pos = 1 the two DIVERGE — only a multi-level fixture sees it.
    np.testing.assert_(leg_offset(1, 0, correct=True) != leg_offset(1, 0, correct=False))
    np.testing.assert_(leg_offset(2, 0, correct=True) != leg_offset(2, 0, correct=False))


# ── Teeth: a production emission-order mutation reds the golden ───────────


def test_offset_golden_reds_under_production_emission_mutation(monkeypatch):
    r"""Corrupting the production leg-emission order REDs the multi-level golden.

    Mutation-proof that :func:`test_multilevel_offsets_match_closed_form` is not
    vacuous. Reversing the module's ``_SIGNS`` constant to ``(+1, -1)`` flips the
    within-leg sign order that :meth:`for_levels` emits — so ``(level, -1)`` and
    ``(level, +1)`` legs swap offsets, diverging from the hard-coded
    :data:`_SIGN_ORDER` golden. ``monkeypatch`` auto-restores ``_SIGNS`` (never a
    ``git checkout`` on an uncommitted file).
    """
    monkeypatch.setattr(sds_mod, "_SIGNS", (+1, -1))
    space = _multilevel_space()
    try:
        _assert_layout_matches_closed_form(space, _LEVELS, _NG, _NX)
    except AssertionError:
        return  # expected — the golden caught the emission-order corruption.
    pytest.fail(
        "the offset golden did NOT red under the _SIGNS emission-order "
        "mutation — the gate is vacuous (a wrong layout would ship green)."
    )
