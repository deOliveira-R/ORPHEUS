r"""Intrinsic-property gates for :func:`orpheus.numerics.roots_of_unity.roots_of_unity`.

This module tests the primitive's **defining laws**, not its usage — the
project standard for any math-bearing type (see the ``test-architect``
memory ``feedback_test_intrinsic_properties``). ``roots_of_unity``
generates the :math:`q`-th roots of unity by the action of the dihedral
group :math:`D_q`, so the laws under test are that group's laws: the
mirror symmetries, the fixed points of the octant swap, closure of the
node set under the action, and the single-orbit structure of the
rotation.

Why an index-permutation mirror and not ``cos(-x) == cos(x)``
============================================================

The mirror MUST be applied as an **index permutation on the generated
set** (``k -> (q - k) % q``). Testing ``cos(-theta) == cos(theta)``
proves nothing: libm's cosine is even *by construction*, so that identity
holds for any implementation whatsoever, including the parameterized one
this module replaces. What production consumes is the node ARRAY, and
what must be bit-exact is the relation between its entries.

The Mode-12 stabiliser note (read before crediting any gate)
============================================================

Two of the laws below are satisfied by implementations that are grossly
wrong, and are therefore kept as *properties*, never cited as catchers:

* :func:`test_unit_norm` — the functional :math:`c^2 + s^2 - 1` has
  invariance group :math:`O(2)` acting on the pair, so it is blind to a
  ``cos``/``sin`` swap, to every sign flip, to conjugation, and even to
  the constant map :math:`(1, 0)`. Measured: it passes 11 of 11 mutants
  of the shipped body. Its ONE tooth is a gross magnitude error
  (measured: reds at ``_SQRT_HALF = 0.707``, green at :math:`\pm1` ULP).
* :func:`test_x_mirror_is_an_exact_index_permutation` and its ``y``
  sibling — a mirror law is satisfied by the constant map and by global
  conjugation (:math:`p \to -p` is an automorphism of :math:`D_q`).
  They are closed by :func:`test_nodes_are_distinct` (kills the constant
  map) and by :func:`test_agrees_with_arbitrary_precision_to_one_ulp`
  (the only pointwise anchor that sees conjugation).

Every ``==`` comparison in this file is blind to signed zero
(``-0.0 == 0.0``), so the ``+ 0.0`` canonicalisation has exactly ONE
catcher: :func:`test_no_negative_zero`, written on :func:`numpy.signbit`.

See Also
--------
:mod:`orpheus.numerics.roots_of_unity`
    The construction and its two integer-decided fixed points.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.roots_of_unity import _SQRT_HALF, roots_of_unity

pytestmark = pytest.mark.foundation


# ``q`` classes, each activating a branch the others null:
#   1, 2       -- r == 0 for EVERY p; the only members that enter the
#                 ``quad == 2`` select lane where ``-sin_quad`` produces
#                 the -0.0 that ``test_no_negative_zero`` pins.
#   4 | q      -- the on-axis fixed point.
#   8 | q      -- the 45-degree diagonal fixed point.
#   odd q      -- the generic octant fold; NO axis beyond p = 0.
#   q = 2 mod 4-- a y-mirror with no x-axis-crossing node.
Q_SWEEP = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 32, 40, 48, 64, 96, 128]
Q_DIVISIBLE_BY_8 = [q for q in Q_SWEEP if q % 8 == 0]
Q_EVEN = [q for q in Q_SWEEP if q % 2 == 0]
Q_ODD = [q for q in Q_SWEEP if q % 2 == 1]

ULP1 = float(np.spacing(1.0))


def _quadrant_and_residual(q: int) -> tuple[np.ndarray, np.ndarray]:
    """``divmod(4p, q)`` for ``p = 0 .. q-1`` — the construction's own
    integer split, recomputed here so the gate does not import it."""
    return np.divmod(4 * np.arange(q), q)


# ---------------------------------------------------------------------------
# L1 -- the axis is exact, and it falls out with NO branch of its own
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_SWEEP)
def test_on_axis_nodes_are_exactly_zero_and_unit(q: int) -> None:
    r"""``4p == 0 (mod q)`` ⟹ the node is exactly ``(±1, 0)`` or ``(0, ±1)``.

    The axis needs no special case: ``r == 0`` gives ``theta == 0.0`` and
    IEEE-754 makes ``cos(0.0) == 1.0``, ``sin(0.0) == 0.0`` exactly. A
    future "simplification" that perturbs ``theta`` by even one ULP
    silently loses this, which is why the law is pinned rather than
    assumed.
    """
    cos, sin = roots_of_unity(np.arange(q), q)
    quad, r = _quadrant_and_residual(q)
    on_axis = r == 0
    lanes = [quad == 0, quad == 1, quad == 2, quad == 3]
    want_cos = np.select(lanes, [1.0, 0.0, -1.0, 0.0])
    want_sin = np.select(lanes, [0.0, 1.0, 0.0, -1.0])
    np.testing.assert_array_equal(
        cos[on_axis], want_cos[on_axis],
        err_msg=f"q={q}: on-axis cosines are not exactly 0/±1",
    )
    np.testing.assert_array_equal(
        sin[on_axis], want_sin[on_axis],
        err_msg=f"q={q}: on-axis sines are not exactly 0/±1",
    )


# ---------------------------------------------------------------------------
# L2 / L3 -- the mirrors, as INDEX PERMUTATIONS on the generated set
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_SWEEP)
def test_x_mirror_is_an_exact_index_permutation(q: int) -> None:
    r"""``k -> (q - k) % q`` maps ``cos`` to itself and ``sin`` to its
    negation, **bit-exactly**.

    This is the reflection across the :math:`x` axis. It is bit-exact
    because the mirror partners fold to the SAME first-octant residual
    and therefore share the two transcendental evaluations — they are not
    computed independently.

    NOT sufficient on its own: the constant map and a global conjugation
    both satisfy it (see the module docstring).
    """
    p = np.arange(q)
    cos, sin = roots_of_unity(p, q)
    mirror = (q - p) % q
    np.testing.assert_array_equal(cos, cos[mirror], err_msg=f"q={q}: x-mirror moved cos")
    np.testing.assert_array_equal(sin, -sin[mirror], err_msg=f"q={q}: x-mirror sin sign")


@pytest.mark.parametrize("q", Q_EVEN)
def test_y_mirror_is_an_exact_index_permutation(q: int) -> None:
    r"""``k -> (q/2 - k) % q`` negates ``cos`` and preserves ``sin``,
    bit-exactly. Defined only for even ``q`` (``pi`` must be on the grid)."""
    p = np.arange(q)
    cos, sin = roots_of_unity(p, q)
    mirror = (q // 2 - p) % q
    np.testing.assert_array_equal(cos, -cos[mirror], err_msg=f"q={q}: y-mirror cos sign")
    np.testing.assert_array_equal(sin, sin[mirror], err_msg=f"q={q}: y-mirror moved sin")


# NOTE — a tautological guard was REMOVED here (2026-08-01).
#
# An "activation guard" asserting that odd ``q`` carries no node at ``pi``
# (``not any(2 * p == q)``) looks like a companion check but is a THEOREM
# ABOUT PARITY: ``2p`` is even and ``q`` is odd, so no input can make it
# fail. It executed under every runtime mode and could never red —
# vv-principles Mode 8, class 2 (tautological companion guard). The
# genuine content, that ``Q_EVEN`` and ``Q_ODD`` partition ``Q_SWEEP``, is
# structural and needs no test.


# ---------------------------------------------------------------------------
# L4 -- the 45-degree diagonal, VALUE-PINNED (not merely |cos| == |sin|)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_DIVISIBLE_BY_8)
def test_diagonal_components_are_the_correctly_rounded_sqrt_half(q: int) -> None:
    r"""``2r == q`` ⟹ both components have magnitude exactly ``sqrt(0.5)``.

    The value is pinned, NOT merely the equality ``|cos| == |sin|``. That
    weaker form has a one-parameter invariance group (any common value
    ``v`` satisfies it) and is measurably blind to a 1-ULP-wrong constant:
    ``1/np.sqrt(2)`` and ``np.sin(np.pi/4)`` are both 1 ULP below
    ``np.sqrt(0.5)`` and pass the equality form.
    """
    cos, sin = roots_of_unity(np.arange(q), q)
    _, r = _quadrant_and_residual(q)
    diagonal = 2 * r == q
    n_diag = int(diagonal.sum())
    if n_diag == 0:
        pytest.fail(f"q={q} is divisible by 8 but has no 45-degree node")
    want = np.full(n_diag, np.sqrt(0.5))
    np.testing.assert_array_equal(
        np.abs(cos[diagonal]), want,
        err_msg=f"q={q}: diagonal |cos| is not the correctly-rounded sqrt(2)/2",
    )
    np.testing.assert_array_equal(
        np.abs(sin[diagonal]), want,
        err_msg=f"q={q}: diagonal |sin| is not the correctly-rounded sqrt(2)/2",
    )


def test_sqrt_half_constant_is_correctly_rounded() -> None:
    r"""``_SQRT_HALF`` is ``np.sqrt(0.5)`` — pinned against a respelling.

    ``sqrt`` is correctly rounded per IEEE-754 and ``0.5`` is exact, so
    ``np.sqrt(0.5)`` is the correctly-rounded :math:`\sqrt2/2`. Two
    natural-looking respellings are 1 ULP LOW and would break the octant
    swap at its own fixed point; two others are bit-equal and safe. The
    enumeration is the point of this test.
    """
    np.testing.assert_array_equal(_SQRT_HALF, np.sqrt(0.5))
    # Bit-equal, therefore safe respellings.
    np.testing.assert_array_equal(np.cos(np.pi / 4), np.sqrt(0.5))
    np.testing.assert_array_equal(np.sqrt(2.0) / 2.0, np.sqrt(0.5))
    # 1 ULP LOW — using either would break the swap's fixed point.
    if float(np.sin(np.pi / 4)) == float(np.sqrt(0.5)):
        pytest.fail("np.sin(pi/4) is no longer 1 ULP below sqrt(0.5) — "
                    "the diagonal patch's rationale needs re-measuring")
    if float(1.0 / np.sqrt(2.0)) == float(np.sqrt(0.5)):
        pytest.fail("1/sqrt(2) is no longer 1 ULP below sqrt(0.5) — "
                    "the diagonal patch's rationale needs re-measuring")


# ---------------------------------------------------------------------------
# L5 -- unit norm.  KEPT AS A PROPERTY; NEVER CREDITED AS A CATCHER.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_SWEEP)
def test_unit_norm(q: int) -> None:
    r"""Every node lies on the unit circle to within 2 ULP.

    HONEST SCOPE (vv-principles Mode 12): the functional
    :math:`c^2 + s^2 - 1` is invariant under the whole of :math:`O(2)`
    acting on the pair. Measured, it passes 11 of 11 mutants of the
    shipped body — including the CONSTANT map :math:`(1, 0)`. Its one
    tooth is a gross magnitude error (reds at ``_SQRT_HALF = 0.707``,
    green at :math:`\pm 1` ULP). Do NOT cite it as evidence for any
    symmetry, orientation, or permutation claim.
    """
    cos, sin = roots_of_unity(np.arange(q), q)
    np.testing.assert_allclose(
        cos * cos + sin * sin, 1.0, rtol=0.0, atol=2 * ULP1,
        err_msg=f"q={q}: nodes left the unit circle",
    )


# ---------------------------------------------------------------------------
# L6 -- the pointwise anchor: agreement with the TRUE value
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", [1, 2, 3, 4, 5, 7, 8, 12, 16, 24, 32, 48, 64, 128])
def test_agrees_with_arbitrary_precision_to_one_ulp(q: int) -> None:
    r"""Every component is within 1 ULP of 1.0 of the arbitrary-precision value.

    The reference is :mod:`mpmath` at 50 digits — structurally
    independent of both this construction and the ``np.linspace`` +
    ``np.cos`` parameterization it replaces.

    **The reference is NOT** ``np.cos(2 * np.pi * p / q)``. That
    expression is not the true value: the argument ``fl(2*pi*p/q)``
    already carries rounding, so the parameterized form is up to
    ``3.72 ulp(1.0)`` from the truth while this construction is within
    ``0.57 ulp(1.0)``. Comparing against numpy would gate the new code
    against the error it exists to remove.

    This is one of only two catchers for a global conjugation
    (:math:`p \to -p`), which every mirror and closure law admits.
    """
    mpmath = pytest.importorskip("mpmath")
    with mpmath.workdps(50):
        cos, sin = roots_of_unity(np.arange(q), q)
        for k in range(q):
            theta = 2 * mpmath.pi * mpmath.mpf(k) / mpmath.mpf(q)
            for got, ref, name in (
                (cos[k], mpmath.cos(theta), "cos"),
                (sin[k], mpmath.sin(theta), "sin"),
            ):
                err = abs(mpmath.mpf(float(got)) - ref)
                if err > ULP1:
                    pytest.fail(
                        f"q={q} k={k} {name}: |{float(got)!r} - exact| = "
                        f"{float(err):.6e} > 1 ulp(1.0) = {ULP1:.6e}"
                    )


# ---------------------------------------------------------------------------
# L7 / L8 / L9 -- the group structure of the node SET
# ---------------------------------------------------------------------------


def _node_set(q: int) -> set[tuple[float, float]]:
    cos, sin = roots_of_unity(np.arange(q), q)
    return {(float(a), float(b)) for a, b in zip(cos, sin)}


@pytest.mark.parametrize("q", Q_SWEEP)
def test_node_set_is_closed_under_the_dihedral_action(q: int) -> None:
    r"""The node set is invariant under the generators of :math:`D_q`:
    the rotation :math:`p \to p+1` and the reflection :math:`p \to -p`.

    This is the group law the whole construction exists to respect. Like
    the mirrors it admits the constant map and global conjugation; it is
    closed by the two tests below plus the arbitrary-precision anchor.
    """
    p = np.arange(q)
    base = _node_set(q)
    cos_r, sin_r = roots_of_unity((p + 1) % q, q)
    rotated = {(float(a), float(b)) for a, b in zip(cos_r, sin_r)}
    cos_m, sin_m = roots_of_unity((-p) % q, q)
    reflected = {(float(a), float(-b)) for a, b in zip(cos_m, sin_m)}
    if rotated != base:
        pytest.fail(f"q={q}: the node set is not closed under rotation")
    if reflected != base:
        pytest.fail(f"q={q}: the node set is not closed under reflection")


@pytest.mark.parametrize("q", Q_SWEEP)
def test_nodes_are_distinct(q: int) -> None:
    r"""There are exactly ``q`` distinct roots.

    The non-degeneracy companion to every symmetry law in this file: the
    constant map :math:`(1, 0)` satisfies both mirrors AND dihedral
    closure, and this is what rules it out. Without it those laws are
    satisfiable by a map that discards all angular information.
    """
    n_distinct = len(_node_set(q))
    if n_distinct != q:
        pytest.fail(f"q={q}: only {n_distinct} distinct nodes — the map is degenerate")


@pytest.mark.parametrize("q", Q_SWEEP)
def test_rotation_acts_as_a_single_q_cycle(q: int) -> None:
    r"""``p -> p+1`` permutes the node set as ONE orbit of length ``q``.

    Closure (previous test) says the rotation maps the set into itself;
    this says it does so as a *generator* — a strictly stronger statement
    that a sign-flipped quadrant lane violates.
    """
    p = np.arange(q)
    cos, sin = roots_of_unity(p, q)
    index = {(float(a), float(b)): i for i, (a, b) in enumerate(zip(cos, sin))}
    if len(index) != q:
        pytest.fail(f"q={q}: node set is degenerate; the orbit test is vacuous")
    cos_r, sin_r = roots_of_unity((p + 1) % q, q)
    try:
        perm = [index[(float(a), float(b))] for a, b in zip(cos_r, sin_r)]
    except KeyError:
        pytest.fail(f"q={q}: rotation left the node set")
    seen: set[int] = set()
    cur = 0
    for _ in range(q):
        if cur in seen:
            break
        seen.add(cur)
        cur = perm[cur]
    if len(seen) != q:
        pytest.fail(
            f"q={q}: rotation orbit has length {len(seen)}, expected {q} "
            f"(the permutation is not a single cycle)"
        )


# ---------------------------------------------------------------------------
# L10 -- signed zero.  The ONLY catcher for the `+ 0.0` canonicalisation.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_SWEEP)
def test_no_negative_zero(q: int) -> None:
    r"""No component is a negative zero.

    Written on :func:`numpy.signbit` on purpose. ``-0.0 == 0.0`` is
    ``True``, so EVERY other comparison in this file — and every
    ``array_equal`` / ``allclose`` in the wider suite — is blind to the
    ``+ 0.0`` canonicalisation. Measured: dropping it leaves ten of the
    thirteen laws green. This is its sole catcher.

    It matters beyond tidiness: a byte-level snapshot hash DOES see the
    sign bit while ``np.array_equal`` does not, so an uncanonicalised
    node array compares equal by value and unequal by ``sha256``.
    """
    cos, sin = roots_of_unity(np.arange(q), q)
    for name, arr in (("cos", cos), ("sin", sin)):
        bad = np.signbit(arr) & (arr == 0.0)
        if np.any(bad):
            pytest.fail(
                f"q={q}: {name} carries negative zero at indices "
                f"{np.flatnonzero(bad).tolist()}"
            )


# ---------------------------------------------------------------------------
# L11 / L12 / L13 -- periodicity, identity, and the group's own product law
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("q", Q_SWEEP)
def test_periodic_in_the_numerator_mod_q(q: int) -> None:
    """``f(p + q) == f(p) == f(p - 3q)`` bit-exactly.

    Negative and out-of-range numerators are accepted and reduced, which
    is what lets ``-p`` name the mirror partner of ``p`` directly.
    """
    p = np.arange(q)
    cos, sin = roots_of_unity(p, q)
    for shift in (q, -3 * q, 7 * q):
        cos_s, sin_s = roots_of_unity(p + shift, q)
        np.testing.assert_array_equal(cos, cos_s, err_msg=f"q={q} shift={shift}: cos")
        np.testing.assert_array_equal(sin, sin_s, err_msg=f"q={q} shift={shift}: sin")


@pytest.mark.parametrize("q", Q_SWEEP)
def test_zeroth_root_is_exactly_one(q: int) -> None:
    r""":math:`\zeta_q^0 = (1, 0)` exactly, for every ``q``."""
    cos, sin = roots_of_unity(np.array([0]), q)
    np.testing.assert_array_equal(cos, np.array([1.0]), err_msg=f"q={q}: zeta^0 cos")
    np.testing.assert_array_equal(sin, np.array([0.0]), err_msg=f"q={q}: zeta^0 sin")


@pytest.mark.parametrize("q", Q_SWEEP)
def test_addition_theorem(q: int) -> None:
    r""":math:`\zeta^p \cdot \zeta^1 = \zeta^{p+1}` as a complex product.

    The genuine group law, and the reason it earns a row of its own: an
    arbitrary relabelling of the circle can satisfy every symmetry law
    above while violating this, because this one constrains the node
    VALUES against each other multiplicatively rather than by symmetry.
    """
    p = np.arange(q)
    cos_p, sin_p = roots_of_unity(p, q)
    cos_1, sin_1 = roots_of_unity(np.array([1]), q)
    cos_next, sin_next = roots_of_unity((p + 1) % q, q)
    prod_cos = cos_p * cos_1[0] - sin_p * sin_1[0]
    prod_sin = sin_p * cos_1[0] + cos_p * sin_1[0]
    np.testing.assert_allclose(
        prod_cos, cos_next, rtol=0.0, atol=8 * ULP1,
        err_msg=f"q={q}: addition theorem failed on the real part",
    )
    np.testing.assert_allclose(
        prod_sin, sin_next, rtol=0.0, atol=8 * ULP1,
        err_msg=f"q={q}: addition theorem failed on the imaginary part",
    )


# ---------------------------------------------------------------------------
# Guards -- the negative tests.  Each calls the PRODUCTION entry point.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "bad_numerators",
    [
        pytest.param(np.array([0.0, 1.0]), id="float64-array"),
        pytest.param(np.array([0.5]), id="fractional"),
        pytest.param(1.5, id="python-float-scalar"),
        pytest.param(np.array([True, False]), id="bool-array"),
    ],
)
def test_non_integer_numerators_are_refused(bad_numerators) -> None:
    """A non-integer numerator RAISES rather than being truncated.

    The exactness guarantee rests on integer arithmetic deciding the two
    fixed points; a silently truncated float would land on the wrong root
    with no diagnostic.

    ``bool`` is included deliberately: ``np.issubdtype(np.bool_,
    np.integer)`` is ``False``, so booleans are refused today. Pinning it
    makes that a decision rather than an accident.
    """
    with pytest.raises(ValueError, match="integer numerators"):
        roots_of_unity(bad_numerators, 8)


@pytest.mark.parametrize("bad_q", [0, -1, -4])
def test_non_positive_denominator_is_refused(bad_q: int) -> None:
    """``denominator < 1`` RAISES — there is no zeroth root of unity."""
    with pytest.raises(ValueError, match="denominator >= 1"):
        roots_of_unity(np.arange(3), bad_q)


@pytest.mark.parametrize(
    "numerators",
    [
        pytest.param(3, id="python-int-scalar"),
        pytest.param(np.int64(3), id="numpy-scalar"),
        pytest.param(np.array(3), id="zero-d-array"),
        pytest.param(np.arange(4), id="one-d"),
        pytest.param(np.arange(6).reshape(2, 3), id="two-d"),
    ],
)
def test_output_shape_matches_input_shape(numerators) -> None:
    """The docstring promises "shape matches ``numerators``" — pin it,
    including the 0-d case where ``np.select`` could plausibly promote."""
    cos, sin = roots_of_unity(numerators, 8)
    want = np.asarray(numerators).shape
    if cos.shape != want or sin.shape != want:
        pytest.fail(
            f"input shape {want} produced cos{cos.shape} / sin{sin.shape}"
        )
