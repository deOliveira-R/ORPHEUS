r"""#337 — the level-symmetric node seed: pre-carve captures and node gates.

Step 0 of the #337 verification plan
(``scratch/issue_337_verification_plan.md`` §8): the frozen pre-carve
data below is captured from the tree BEFORE the moment-matched seed
lands, in its own commit — a post-change capture would be a
self-referential snapshot certifying nothing.

At step 2 (the atomic landing) this module gains the G1–G8 gate set
(the μ₁ table against the 50-digit re-solve, the degree table both
directions, the smallest-root selection rule, root-find totality, the
``∫μ_z^N`` naming condition, the frontier flip). Today it asserts the
captures against the current tree — trivially green, which is exactly
what a recording commit should be.
"""

from __future__ import annotations

import numpy as np

from orpheus.numerics.quadrature import Quadrature

# ── Frozen pre-carve data ───────────────────────────────────────────
# Captured 2026-08-08 from the PRE-#337 tree (seed μ₁² = 4/(N(N+2)),
# commit bad5223e) — step 0 of the #337 verification plan, BEFORE the
# seed change lands. repr() round-trips float64 exactly.

_S2_NODES_PRE337 = np.array([
    [-0.5773502691896257, -0.5773502691896257, -0.5773502691896257],
    [-0.5773502691896257, -0.5773502691896257, 0.5773502691896257],
    [-0.5773502691896257, 0.5773502691896257, -0.5773502691896257],
    [-0.5773502691896257, 0.5773502691896257, 0.5773502691896257],
    [0.5773502691896257, -0.5773502691896257, -0.5773502691896257],
    [0.5773502691896257, -0.5773502691896257, 0.5773502691896257],
    [0.5773502691896257, 0.5773502691896257, -0.5773502691896257],
    [0.5773502691896257, 0.5773502691896257, 0.5773502691896257],
])

_S2_WEIGHTS_PRE337 = np.array([
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
    1.5707963267948966,
])

# level_indices per order — `[M]` the #337 plan measured these
# BIT-IDENTICAL old-seed vs new-seed at S4..S12, so this literal is the
# permutation-stability gate for the carve (G7-family): the seed change
# moves node VALUES, never the level membership or intra-level order.
_PRE337_LEVEL_INDICES = {
    4: ((8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15), (16, 17, 18, 19, 20, 21, 22, 23)),
    6: ((19, 18, 17, 16, 11, 10, 9, 8, 0, 3, 2, 1, 6, 5, 4, 7, 12, 13, 14, 15, 22, 20, 21, 23), (32, 33, 34, 35, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37, 38, 39), (40, 41, 42, 43, 44, 45, 46, 47)),
    8: ((24, 25, 26, 27, 19, 17, 16, 18, 11, 10, 9, 8, 3, 2, 1, 0, 7, 6, 5, 4, 15, 14, 13, 12, 21, 22, 23, 20, 30, 28, 29, 31), (51, 50, 49, 48, 43, 42, 41, 40, 32, 35, 34, 33, 38, 37, 36, 39, 44, 45, 46, 47, 54, 52, 53, 55), (64, 65, 66, 67, 56, 57, 58, 59, 60, 61, 62, 63, 68, 69, 70, 71), (72, 73, 74, 75, 76, 77, 78, 79)),
    10: ((35, 34, 33, 32, 27, 26, 25, 24, 19, 18, 17, 16, 8, 9, 10, 11, 0, 1, 2, 3, 7, 6, 5, 4, 14, 13, 12, 15, 20, 21, 22, 23, 28, 30, 29, 31, 38, 36, 37, 39), (64, 65, 66, 67, 59, 57, 56, 58, 51, 50, 49, 48, 43, 42, 41, 40, 47, 46, 45, 44, 55, 54, 53, 52, 61, 62, 63, 60, 70, 68, 69, 71), (91, 90, 89, 88, 83, 82, 81, 80, 72, 75, 74, 73, 78, 77, 76, 79, 84, 85, 86, 87, 94, 92, 93, 95), (104, 105, 106, 107, 96, 97, 98, 99, 100, 101, 102, 103, 108, 109, 110, 111), (112, 113, 114, 115, 116, 117, 118, 119)),
    12: ((40, 43, 41, 42, 32, 34, 35, 33, 27, 26, 25, 24, 17, 16, 19, 18, 10, 8, 11, 9, 0, 3, 2, 1, 7, 6, 5, 4, 13, 12, 14, 15, 23, 22, 21, 20, 30, 29, 28, 31, 36, 37, 38, 39, 46, 44, 45, 47), (83, 82, 81, 80, 75, 74, 73, 72, 67, 66, 65, 64, 56, 57, 58, 59, 48, 49, 50, 51, 55, 54, 53, 52, 62, 61, 60, 63, 68, 69, 70, 71, 76, 78, 77, 79, 86, 84, 85, 87), (112, 113, 114, 115, 107, 105, 104, 106, 99, 98, 97, 96, 91, 90, 89, 88, 95, 94, 93, 92, 103, 102, 101, 100, 109, 110, 111, 108, 118, 116, 117, 119), (139, 138, 137, 136, 131, 130, 129, 128, 120, 123, 122, 121, 126, 125, 124, 127, 132, 133, 134, 135, 142, 140, 141, 143), (152, 153, 154, 155, 144, 145, 146, 147, 148, 149, 150, 151, 156, 157, 158, 159), (160, 161, 162, 163, 164, 165, 166, 167)),
}


# ── G6 — the S2 bit-identity control ────────────────────────────────

def test_S2_is_BIT_IDENTICAL_to_the_pre_337_rule():
    """The regression control: S2 has no freedom (μ² = 1/3 forced by
    the diffusion condition), so a node change there would mean the
    RECURSION itself moved, not the seed.

    Frozen pre-carve arrays, captured BEFORE the change lands (a
    self-recorded post-change snapshot would certify nothing).
    """
    q = Quadrature.level_symmetric(2)
    np.testing.assert_array_equal(
        np.asarray(q.measure.nodes), _S2_NODES_PRE337
    )
    np.testing.assert_array_equal(
        np.asarray(q.measure.weights), _S2_WEIGHTS_PRE337
    )


def test_S2_nodes_are_the_forced_diffusion_value():
    """WHY the row above is exact: every S2 coordinate is ±√(1/3) —
    the forced value — so the gate survives any node re-ordering."""
    q = Quadrature.level_symmetric(2)
    np.testing.assert_array_equal(
        np.abs(np.asarray(q.measure.nodes)),
        np.full((8, 3), np.sqrt(1.0 / 3.0)),
    )


# ── G7-family — the level structure is seed-independent ─────────────

def test_level_indices_are_the_pre_337_permutation():
    """`[M]` (#337 plan §0): ``level_indices`` is bit-identically the
    same permutation under the old and the moment-matched seed at
    every order S4..S12 — level membership is |μ_z|-keyed index
    arithmetic and the intra-level order is the stored sort, neither
    of which reads the seed VALUE. This literal pins that: the seed
    carve moves node values, never the level structure.
    """
    for n, expected_levels in _PRE337_LEVEL_INDICES.items():
        q = Quadrature.level_symmetric(n)
        assert len(q.level_indices) == len(expected_levels), f"S{n}"
        for p, (got, exp) in enumerate(
            zip(q.level_indices, expected_levels, strict=True)
        ):
            np.testing.assert_array_equal(
                np.asarray(got), np.asarray(exp, dtype=np.asarray(got).dtype),
                err_msg=f"S{n} level {p}",
            )
