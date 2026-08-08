r"""#337 — the moment-matched level-symmetric seed: the node gate set.

The seed :math:`\mu_1` is the smallest root of "the rule integrates
:math:`\mu_z^N` exactly" (G5's law) — the construction behind the
published LA-3186 tables. This module carries the G1–G8 gates of the
#337 verification plan (``scratch/issue_337_verification_plan.md``)
plus the pre-carve captures its step 0 froze BEFORE the change landed.

⭐ **Honest scope (the plan's R5).** The degree table (G2) is
structurally BLIND to the seed at N = 12, 16 and 18: the achieved
degree is :math:`N-1` under BOTH the retired convention seed and the
moment-matched one (`[M]` 11/11, 15/15, 17/17), so no reading of the
degree can distinguish them there. Those orders are covered by G1 (the
μ₁ value against the 50-digit re-solve) and G5 (the axial residual)
alone — and G1's own tolerance ladder is loosest exactly there (the
residual slope collapses 1.03 → 1.7e-4 across the family), so a
sub-1e-10 μ₁ error at S16/S18 is ungated by every row in this module.
That is §3.3's arithmetic floor, not a gate weakness to fix.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.quadrature.rules_sphere import (
    _greedy_orbit_solve,
    _octant_directions,
    _sphere_monomial_integral,
)

pytestmark = pytest.mark.foundation

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


# ── G1 — the μ₁ table (the value keystone) ──────────────────────────

#: `[M]` 2026-08-08. Derived, never transcribed: an independent
#: re-implementation of the level recursion + greedy solve, RE-SOLVED
#: at 50 digits in mpmath and rounded to float64 (the plan's
#: scratch/_probe_337_prec2.py). These are the MATHEMATICALLY CORRECT
#: values, not our float64 root-finder's output. Corroboration
#: (render-verified against the LA-3186 scan, 2026-08-08 — NOT the
#: reference): Table I tabulates S4/S6/S8/S12/S16/S20 and agrees at
#: every one, with THREE one-ulp last-digit print slips (S6/S12/S16 —
#: these values are the correctly-rounded ones, the print is off);
#: S10/S14/S18 are untabulated and corroborated by an independent
#: level-system reproduction of the C-L construction
#: (scratch/issue_337_la3186_la4058_extraction.md) — the two
#: derivations share the level system, so that leg is procedural
#: corroboration, not structural.
_MU1 = {
    4: 0.3500211745815407,    # == sqrt((5 - sqrt(10)) / 15), 1.0 ULP
    6: 0.26663540151670473,
    8: 0.2182178902359924,    # mu1**2 == 1/21 exactly
    10: 0.18932132647801048,
    12: 0.16721265282271328,
    14: 0.1519858614610319,
    16: 0.13895687506778034,
    18: 0.12934450454592483,
}

#: ⭐ PER-ORDER, and DERIVED (plan §3.3). The float64 root-find's
#: distance from the 50-digit truth grows 1.0 → 40 653 ULP across the
#: family because the residual slope collapses (1.03 → 1.7e-4). A
#: single rtol is either a false red at S14+ or a thrown-away gate at
#: S4. `[M]` measured error ×10, decade-rounded.
_MU1_RTOL = {4: 1e-14, 6: 1e-13, 8: 1e-12, 10: 1e-11,
             12: 1e-11, 14: 1e-10, 16: 1e-10, 18: 1e-10}


def test_S4_mu1_is_the_EXACT_RADICAL():
    r"""⭐ The only order with a closed form — the fully-structural anchor.

    One :math:`O_h` orbit ⟹ :math:`\Sigma w = 4\pi` forces the weight,
    so μ₁ is the SOLE freedom and the μ⁴ condition closes in radicals:
    :math:`15x^2 - 10x + 1 = 0`, :math:`x = \mu_1^2 = (5-\sqrt{10})/15`.
    Structurally independent of the SUT — a quadratic solved by
    radicals, not by the builder's linear algebra. Every other order is
    a numerical re-solve (procedurally independent only), so this row
    carries the family's one genuine structural leg.
    """
    want = math.sqrt((5.0 - math.sqrt(10.0)) / 15.0)
    got = float(
        np.min(np.abs(np.asarray(Quadrature.level_symmetric(4).measure.nodes)))
    )
    # `[M]` 1.0 ULP against the 50-digit solve; 5e-15 is ~9x that.
    np.testing.assert_allclose(got, want, rtol=0.0, atol=5e-15)


@pytest.mark.parametrize("n", sorted(_MU1))
def test_mu1_matches_the_50_DIGIT_solve(n):
    got = float(
        np.min(np.abs(np.asarray(Quadrature.level_symmetric(n).measure.nodes)))
    )
    np.testing.assert_allclose(got, _MU1[n], rtol=_MU1_RTOL[n], atol=0.0)


# ── G2 — the degree table's frozen third corner ─────────────────────

#: `[M]` 2026-08-08, closed-form monomial sweep at atol 1e-10,
#: cross-checked at 50 digits. NOT a formula: N+1 at 4/6/8/10/14 and
#: N-1 at 12/16/18 (see the module docstring's honest-scope note).
_LS_DEGREE = {2: 3, 4: 5, 6: 7, 8: 9, 10: 11, 12: 11, 14: 15, 16: 15, 18: 17}


def _own_degree_sweep(nodes: np.ndarray, weights: np.ndarray) -> int:
    """This module's OWN monomial sweep — the triangle's second edge.

    With a build-measured stamp, "measured == advertised" alone is two
    re-implementations of one gamma identity agreeing (the plan's R3);
    the frozen literal plus this independent sweep make it a three-way
    consistency triangle: a wrong production sweep reddens the stamp
    row only, a wrong test sweep reddens this row only, and a genuine
    family change reddens both — the signal the literals need
    re-deriving.
    """
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    degree = 0
    while degree < 64:
        d = degree + 1
        for a in range(d + 1):
            for b in range(d + 1 - a):
                c = d - a - b
                s = float(np.sum(weights * x**a * y**b * z**c))
                if abs(s - _sphere_monomial_integral(a, b, c)) > 1e-10:
                    return degree
        degree = d
    return degree


@pytest.mark.parametrize("n", sorted(_LS_DEGREE))
def test_the_stamp_matches_the_FROZEN_table(n):
    """Production's stamp == the literal (catches a wrong build-time
    measurement without dragging this module's own sweep into it)."""
    got = Quadrature.level_symmetric(n).measure.exactness.degree
    assert got == _LS_DEGREE[n], (
        f"S{n}: stamped {got}, frozen literal {_LS_DEGREE[n]}. If the "
        f"OWN-sweep row reddened too, the family changed and the "
        f"literals need re-deriving; if only this row, the build-time "
        f"measurement is wrong."
    )


@pytest.mark.parametrize("n", sorted(_LS_DEGREE))
def test_the_MEASUREMENT_matches_the_FROZEN_table(n):
    """This module's own sweep == the literal (catches a wrong test
    measure)."""
    q = Quadrature.level_symmetric(n)
    got = _own_degree_sweep(
        np.asarray(q.measure.nodes), np.asarray(q.measure.weights)
    )
    assert got == _LS_DEGREE[n], f"S{n}: swept {got} != {_LS_DEGREE[n]}"


# ── G3 — the smallest-root selection rule ───────────────────────────

#: `[M]` 2026-08-08 — the OTHER root of the (0,0,N) residual on
#: (0, 1/3) in μ₁. Exists only at N = 6, 10, 14, 18. Every one yields
#: a strongly NEGATIVE weight, which is why "the smallest root" is the
#: rule (plan §10 — "smallest positive-weight root" was REJECTED:
#: positivity belongs to the refusal contract, not the selection).
_MU1_SECOND_ROOT = {6: 0.4225186537611119, 10: 0.3471231056,
                    14: 0.3008785706, 18: 0.2689383416}


def _weights_at_seed(n: int, mu1_sq: float) -> np.ndarray:
    """The production per-orbit solve at a HAND-FED seed (no root-find,
    no positivity contract) — what the builder would have shipped had
    it taken this seed."""
    _levels, dirs, orbit = _octant_directions(n, mu1_sq)
    per_orbit, _labels, _conds, _system, _targets = _greedy_orbit_solve(
        n, dirs, orbit
    )
    return np.asarray(per_orbit)


@pytest.mark.parametrize("n", sorted(_MU1_SECOND_ROOT))
def test_the_SMALLEST_root_is_the_one_taken(n):
    r"""⭐ The μ^N residual has TWO roots on (0, 1/3); only one is a rule.

    `[M]` the second root at each of these orders yields a solve with a
    strongly NEGATIVE weight (−0.597 / −1.087 / −3.428 / −7.669), so
    taking it would either ship a negative-weight rule or make the
    family refuse an order it can serve. A construction that says only
    "root-find on (0, 1/3)" does not distinguish them — and a global
    bracketing solve cannot even be CALLED there (same sign at both
    ends).

    Two legs, because either alone is satisfiable by accident:
    (a) the shipped μ₁ is the SMALLER of the two roots; (b) the LARGER
    root, constructed by hand, has a negative weight — the reason (a)
    is the right rule and not a coincidence.
    """
    shipped = float(
        np.min(np.abs(np.asarray(Quadrature.level_symmetric(n).measure.nodes)))
    )
    assert shipped == pytest.approx(_MU1[n], rel=1e-9)
    assert shipped < _MU1_SECOND_ROOT[n]
    w2 = _weights_at_seed(n, _MU1_SECOND_ROOT[n] ** 2)
    assert float(np.min(w2)) < 0.0, "the second root is not the discarded one"


# ── G4 — the root-find is TOTAL over the orders it serves ───────────

@pytest.mark.parametrize("n", [4, 6, 8, 10, 12, 14, 16, 18])
def test_the_construction_does_not_raise_a_solver_error(n):
    """The plan's R2: 31 % of the S20 bracket is rank-deficient for the
    greedy solve and the residual is discontinuous where its chosen set
    switches. The family must fail ONLY through the positivity channel
    with the documented message — never through LinAlgError,
    RuntimeError, or a silent NaN out of the inner solve. This row is
    why the outer root-find brackets LOCALLY (stepping over NaN trial
    points) rather than solving over (0, 1/3) blind.
    """
    q = Quadrature.level_symmetric(n)
    assert np.all(np.isfinite(np.asarray(q.measure.nodes)))
    assert np.all(np.isfinite(np.asarray(q.measure.weights)))


def test_the_refusal_channel_is_the_ONLY_failure_channel():
    """`match=` is load-bearing: a bare raises(ValueError) is satisfied
    by a rank-deficiency raise from the inner solve — exactly the
    failure R2 predicts — and the gate would be green on the wrong bug.
    """
    for n in (20, 22, 24):
        with pytest.raises(ValueError, match="no POSITIVE solution"):
            Quadrature.level_symmetric(n)


# ── G5 — the seed's LAW: ∫μ_z^N is exact ────────────────────────────

@pytest.mark.parametrize("n", [4, 6, 8, 10, 12, 14, 16, 18])
def test_the_axial_monomial_of_degree_N_is_exact(n):
    r"""⭐ The law the seed is chosen to satisfy, stated directly.

    `[M]` at every order 4..20 the first ladder condition the weight
    solve cannot satisfy on its own is exactly (0, 0, N) — the pure
    axial monomial. So "the moment-matched seed" has a NAME: the rule
    integrates :math:`\mu_z^N` exactly. This row makes the construction
    checkable without re-running its own root-find; it is the equation,
    not the algorithm. (Both roots satisfy it by definition — the
    selection claim belongs to G3, not here.)

    ⭐ This is the gate that carries S12/S16/S18, where the degree table
    is blind (module docstring). Its detection floor per order tracks
    the residual slope: ~1e-13 relative μ₁ error detectable at S4,
    degrading to ~5.8e-10 at S18.

    Tolerance 1e-13 RELATIVE, derived: the sum is N(N+2) POSITIVE
    terms (κ = 1), so the FP floor is n_terms·eps; `[M]` the realised
    residual at the exact seed is 3.5e-16 … 7.7e-15 relative across
    the family — 13× headroom at the worst order.
    """
    q = Quadrature.level_symmetric(n)
    nodes = np.asarray(q.measure.nodes)
    w = np.asarray(q.measure.weights)
    got = float(np.sum(w * nodes[:, 2] ** n))
    want = _sphere_monomial_integral(0, 0, n)
    assert abs(got - want) <= 1e-13 * want, (
        f"S{n}: sum w*mu_z^{n} = {got!r} vs closed form {want!r} — the "
        f"seed does not satisfy its own defining condition"
    )


# ── G8 — the frontier flip, two-sided ───────────────────────────────

def test_S18_BUILDS_and_S20_REFUSES():
    """The frontier moved from S12 to S18 — the point of #337.

    `[M]` min weight, VERIFIED AT 50 DIGITS (plan §3.2):
    S18 = +0.0001750014795257 POSITIVE, S20 = −0.0002191020803263
    NEGATIVE; float64 agrees to every printed digit and the float64
    root's margin is 2.3e-11 from the exact one — the sign has ~7
    orders of margin, so the frontier is a property of the FAMILY, not
    of our arithmetic. (The builder still COMPUTES it from its own
    solution's sign — this row pins the behaviour at 18 and 20, not a
    hardcoded bound.)
    """
    q18 = Quadrature.level_symmetric(18)
    assert float(np.min(np.asarray(q18.measure.weights))) > 0.0
    with pytest.raises(ValueError, match="no POSITIVE solution") as exc:
        Quadrature.level_symmetric(20)
    assert "lebedev" in str(exc.value) or "product" in str(exc.value)


@pytest.mark.parametrize("n,want", [(14, 1.298982584513e-2),
                                    (16, 1.6300087256084e-2),
                                    (18, 1.750014795257e-4)])
def test_the_frontier_MARGIN_is_pinned_not_just_its_sign(n, want):
    """⭐ A sign decided by 2e-4 goes stale silently; a VALUE does not.

    S18's +1.75e-4 is 93× smaller than S16's +1.63e-2. If a future
    change to pivoting / BLAS / the greedy rank tolerance perturbs the
    solve at the 1e-4 level, the SIGN row above flips with no warning
    and the family silently sheds an order. This row reddens FIRST,
    with the number in the message.

    Tolerance rel=1e-5: `[M]` the min weight at the float64 root
    differs from the 50-digit one by 1.33e-7 relative at S18 (1.6e-11
    at S16, 5.7e-11 at S14) — 75× headroom at the worst order. The
    ``want`` values are the 50-DIGIT ones.
    """
    got = float(np.min(np.asarray(Quadrature.level_symmetric(n).measure.weights)))
    assert got == pytest.approx(want, rel=1e-5)
