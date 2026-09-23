r"""Foundation gates for :mod:`orpheus.geometry.transformation` — campaign step G2.

The 42 gates of the G2 verification plan (``scratch/g2_verification_plan.md``),
written against **pure mathematics**, never against
:mod:`orpheus.numerics.symmetry`. That prohibition is the whole reason G2 runs
before G7 (test migration): once the primitive is verified against
domain-independent ground truth, a system-under-test and its oracle both using
it is not contamination — independence lives in the INPUT ASSEMBLY, not in the
shared routine. Verified against ``symmetry.py`` instead, the tree would lose
its independent check silently.

**Forbidden references in this file** (grep-checkable): ``numerics.symmetry``,
``_rotation_z``, ``_reflections``, ``_orbit_closure``, ``_close_group``,
``Quadrature``, ``DiscreteMeasure``, ``roots_of_unity``. The G3 bit-identity
gate ("the migration is a no-op") is a *different* claim and belongs in a
different file.

The local fixture helper is deliberately named ``_reflection_family`` rather
than ``_reflections``: the audit above is meant to be a MECHANICAL grep, and a
local symbol sharing a forbidden name would defeat it.

The reference pillars actually available for a group-theory primitive — there
is no MMS row and no semi-analytical row, because there is no differential
operator to manufacture a solution for and nothing reduces to a quadrature:

===  =========================================================================
P1   SymPy symbolic identity under an EXPLICIT unit parameterisation
     (``n̂ = (sinθcosφ, sinθsinφ, cosθ)``). Exact algebra, no floats.
P2   ``scipy.spatial.transform.Rotation.from_rotvec`` — an independently
     maintained, quaternion-based implementation: a different ALGORITHM,
     not a re-spelling.
P3   ``scipy.linalg.expm(θ(vuᵀ − uvᵀ))`` — the Lie definition ``R = exp(θJ)``
     by Padé scaling-and-squaring. The only DIMENSION-GENERIC reference.
P4   Closed forms written independently in the test (the projector form of a
     Householder; ``[[c,−s],[s,c]]``; ``arange(n)[::-1]``).
P5   Published tables: group orders, and the cube's Wyckoff site symmetries
     (Hamermesh 1962 §9.4; International Tables Vol. A).
P6   Exact integer arithmetic — a permutation's bijectivity,
     ``π(gh) = π(g)∘π(h)``, ``ord(g^k) = n/gcd(k,n)``. No tolerance possible.
===  =========================================================================

Gate map (plan §4). ``O`` is the orbit/action group, lettered so no gate ID
collides with a campaign STEP (``G1``…``G8``):

* ``A1``–``A6`` group axioms; ``B1``–``B4`` orthogonality/determinant
* ``C1``–``C5`` element order and the fixed set
* ``D1``–``D5`` independent constructions; ``E1``–``E3`` conjugation
* ``F1``–``F4`` seating; ``O1``–``O8`` ``permutes``/``preserves``/orbits
* ``X1``–``X3`` gates the landed design earns; ``H1``–``H3`` dimension genericity

Every tolerance below is MEASURED (see ``_TOL_*`` provenance comments), never
guessed, and the bit-exact-versus-tolerance choice is per-law: identity and
signed-permutation associativity are bit-exact 500/500, while general-element
associativity is bit-exact 0/500. A uniform choice would be either a false red
or a thrown-away gate.

Tagged ``@pytest.mark.foundation`` per :mod:`tests._harness` — these are
mathematical invariants of a primitive with no theory-page ``:label:`` today,
so they carry no ``verifies(...)``.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest
import sympy as sp
from numpy.typing import ArrayLike
from scipy.linalg import expm
from scipy.spatial.transform import Rotation as ScipyRotation

from orpheus.geometry import transformation as _txn
from orpheus.geometry.transformation import (
    NotAFinitePointGroupError,
    Permutation,
    RigidMotion,
    close_group,
)

pytestmark = pytest.mark.foundation


# ---------------------------------------------------------------------------
# Measured tolerances — each carries its provenance
# ---------------------------------------------------------------------------

#: General-element residuals, **relative to the quantity's own scale**.
#:
#: An absolute bound is the wrong contract for a translation: the residual of
#: a product is ``O(ops × ‖t‖ × eps)``, so a fixture drawing ``t ~ N(0, 1)``
#: produces residuals that track ``‖t‖`` and a fixed ``atol`` reds on correct
#: code for the large draws. MEASURED over 3000-4000 draws per law at
#: d = 1, 2, 3, normalising by ``max(1, ‖desired‖_∞)``: associativity ≤7.5e-16,
#: inverse ≤4.6e-15, anti-homomorphism ≤5.2e-15, affine conjugation ≤6.9e-15,
#: seating ≤2.4e-16. Worst overall 6.9e-15, so this leaves ~14× margin over
#: the FP floor while every mutation the gates target is O(1).
_TOL_GENERAL = 1e-13

#: Cross-check against ``scipy…Rotation.from_rotvec``. MEASURED max abs
#: 1.1e-15 over 300 draws against the SHIPPED ``rotation_about_axis``.
#: ABSOLUTE, never ULP: the same data reads "5120 ULP" because the matrix
#: entries cross zero, and ULP is the wrong metric for arrays crossing zero.
_TOL_SCIPY_ROTVEC = 1e-14

#: Cross-check against ``scipy.linalg.expm`` of the planar generator.
#: MEASURED 3.9e-14 over 300 draws at d = 2, 3, 4 against the SHIPPED
#: ``rotation``. Padé scaling-and-squaring is not exact; a 1e-15 gate here
#: would be a false red.
_TOL_EXPM = 1e-12

#: The ``permutes`` match window for exactly-constructed point sets (cube
#: vertices, square, centred lattices). MEASURED: a seated reflection of a
#: uniform cell-centre lattice lands its images 5.6e-17 – 2.2e-16 from their
#: partners — never bit-exact — while the induced PERMUTATION is exactly
#: ``arange(n)[::-1]``. Integers exact, coordinates not.
_TOL_MATCH = 1e-12

#: Angles bounded away from zero for the fixed-subspace-dimension gate.
#: ``fixed_subspace_dimension`` uses ``matrix_rank(Q − I, tol=1e-10)``, an
#: ABSOLUTE singular-value threshold, so a rotation by an angle small compared
#: with it is indistinguishable from the identity — a documented property,
#: pinned as an honest-scope control in ``test_c4_...degenerate_control``.
_ROTATION_ANGLES = (0.1, 0.7, 1.3, 2.0, 3.0)


def _rng(seed: int = 20260803) -> np.random.Generator:
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# Fixtures — built from pure math, never from orpheus.numerics
# ---------------------------------------------------------------------------


def _assert_close(
    actual: "ArrayLike", desired: "ArrayLike", *, err_msg: str = ""
) -> None:
    """``assert_allclose`` with the tolerance scaled by the quantity's own size.

    See :data:`_TOL_GENERAL`: the residual of a rigid-motion product scales
    with ``‖t‖``, so a fixed absolute bound is a false red on large draws and
    an over-loose gate on small ones.
    """
    target = np.asarray(desired, dtype=float)
    scale = max(1.0, float(np.max(np.abs(target))) if target.size else 1.0)
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float),
        target,
        atol=_TOL_GENERAL * scale,
        rtol=0.0,
        err_msg=err_msg,
    )


def _orthonormal_plane(rng: np.random.Generator, d: int) -> np.ndarray:
    """A random orthonormal pair ``(u, v)`` spanning a 2-plane of R^d.

    Built by QR, not by one Gram-Schmidt pass. MEASURED: a single
    orthogonalisation of a random ``v`` against ``u`` loses enough
    orthogonality on near-parallel draws to reach ``max|GGᵀ − I| = 2.1e-12``
    — ABOVE the shipped ``_ORTHOGONALITY_ATOL`` — so the production
    constructor (correctly) refuses the fixture. The refusal is right; the
    fixture was wrong.
    """
    q, _ = np.linalg.qr(rng.normal(size=(d, 2)))
    return q.T


def _random_rotation(rng: np.random.Generator, d: int) -> RigidMotion:
    return RigidMotion.rotation(
        plane=_orthonormal_plane(rng, d), angle=rng.uniform(-np.pi, np.pi)
    )


def _random_linear(rng: np.random.Generator, d: int) -> RigidMotion:
    """A random origin-fixing element, valid at every ``d`` including 1.

    ``SO(1) = {e}``, so at d = 1 the whole of ``O(1)`` is ``{±1}`` and the
    only non-trivial element is the reflection. Reaching for a rotation there
    is the 1-D/3-D arm split creeping back into the *test* fixture.
    """
    if d == 1:
        return (
            RigidMotion.reflection(normal=[1.0])
            if rng.integers(2)
            else RigidMotion.identity(1)
        )
    return _random_rotation(rng, d)


def _random_motion(rng: np.random.Generator, d: int) -> RigidMotion:
    """A general element of E(d): a random linear part with a random ``t``."""
    return RigidMotion(_random_linear(rng, d).linear, rng.normal(size=d))


#: Reflection normals per dimension. **Every reflection gate must carry at
#: least one NON-COORDINATE normal**: with ``n̂ = e_x`` the Householder
#: degenerates to a diagonal sign flip, which is exactly the degenerate case
#: this campaign is dissolving — a coordinate-only fixture is structurally
#: blind to the general ``I − 2n̂n̂ᵀ`` path. (Plan §7-V10; MEASURED: ``σ² == I``
#: is bit-exact for a coordinate normal and 0/300 for general normals.)
_NORMALS: dict[int, tuple[tuple[str, list[float]], ...]] = {
    1: (("e0", [1.0]), ("neg", [-2.5])),
    2: (("e0", [1.0, 0.0]), ("diagonal", [1.0, 1.0]), ("generic", [1.0, 2.0])),
    3: (
        ("e2", [0.0, 0.0, 1.0]),
        ("diagonal", [1.0, 1.0, 0.0]),
        ("generic", [1.0, 2.0, 3.0]),
        ("body", [1.0, 1.0, 1.0]),
    ),
}


def _reflection_family(d: int) -> list[tuple[str, RigidMotion]]:
    return [
        (name, RigidMotion.reflection(normal=n)) for name, n in _NORMALS[d]
    ]


def _hyperoctahedral(d: int) -> list[RigidMotion]:
    """All ``2^d · d!`` signed permutations — ``O_h`` at d = 3, ``D_4`` at d = 2.

    Exactly orthogonal (integer entries), and NON-ABELIAN for ``d >= 2``,
    which several gates require: on an abelian group the composition-order
    mutation moves nothing and the homomorphism law is vacuous (plan §7-V3).
    """
    return [
        RigidMotion.signed_permutation(permutation=perm, signs=signs)
        for signs in itertools.product([-1.0, 1.0], repeat=d)
        for perm in itertools.permutations(range(d))
    ]


def _cube_vertices() -> np.ndarray:
    return np.array(list(itertools.product([-1.0, 1.0], repeat=3))) / np.sqrt(3.0)


def _cube_face_centres() -> np.ndarray:
    return np.vstack([np.eye(3), -np.eye(3)])


def _cube_edge_midpoints() -> np.ndarray:
    raw = (
        [[a, b, 0.0] for a in (-1.0, 1.0) for b in (-1.0, 1.0)]
        + [[a, 0.0, b] for a in (-1.0, 1.0) for b in (-1.0, 1.0)]
        + [[0.0, a, b] for a in (-1.0, 1.0) for b in (-1.0, 1.0)]
    )
    return np.array(raw) / np.sqrt(2.0)


def _square_vertices() -> np.ndarray:
    return np.array([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]) / np.sqrt(2.0)


def _square_edge_midpoints() -> np.ndarray:
    return np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])


def _cell_centres(n: int, *, origin: float = 0.0, h: float = 0.1) -> np.ndarray:
    """A 1-D uniform cell-centre lattice, as ``Mesh1D`` lays one out."""
    return (origin + (np.arange(n) + 0.5) * h).reshape(-1, 1)


#: Wyckoff decomposition of the cube's special positions under the full
#: hyperoctahedral group (P5 — International Tables Vol. A; Hamermesh 1962
#: §9.4). ``(|orbit|, |Stab|)``, whose product is |G|.
_WYCKOFF_CUBE = {
    "vertices (site C_3v)": (_cube_vertices, 8, 6),
    "face centres (site C_4v)": (_cube_face_centres, 6, 8),
    "edge midpoints (site C_2v)": (_cube_edge_midpoints, 12, 4),
}
_WYCKOFF_SQUARE = {
    "vertices (site C_s)": (_square_vertices, 4, 2),
    "edge midpoints (site C_s)": (_square_edge_midpoints, 4, 2),
    "centre (site D_4)": (lambda: np.zeros((1, 2)), 1, 8),
}


# ===========================================================================
# A — group axioms of the rigid-motion algebra
# ===========================================================================


@pytest.mark.parametrize("d", [1, 2, 3])
def test_a1_associativity_is_bit_exact_on_signed_permutations(d: int) -> None:
    r"""A1, tier 1 — ``(g∘h)∘k = g∘(h∘k)`` BIT-IDENTICALLY on the integer subgroup.

    The two tiers are MEASURED, not assumed: signed permutations seated at
    integer centres carry only integer-valued entries, so every product is
    exact and ``array_equal`` is the honest gate. General elements are
    bit-exact 0/500 (tier 2 below).

    Falsifier: composing as ``(Q₁Q₂, Q₂t₁ + t₂)`` — the operand swap in the
    translation rule — breaks associativity by O(1). Nothing else in this file
    pins the ORDER inside the translation rule except A4 and O2.
    """
    seat = np.arange(1, d + 1, dtype=float)  # integer centre: (I−Q)c stays exact
    elements = [g.seated_at(seat) for g in _hyperoctahedral(d)]
    rng = _rng()
    for _ in range(60):
        a, b, c = (elements[i] for i in rng.integers(len(elements), size=3))
        left, right = (a @ b) @ c, a @ (b @ c)
        np.testing.assert_array_equal(
            left.linear, right.linear, err_msg="associativity moved the linear part"
        )
        np.testing.assert_array_equal(
            left.translation,
            right.translation,
            err_msg="associativity moved the translation part",
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_a1_associativity_holds_for_general_elements(d: int) -> None:
    """A1, tier 2 — general elements, to a MEASURED tolerance.

    MEASURED over 500 random triples: ``‖ΔQ‖ ≤ 2.2e-16``, ``‖Δt‖ ≤ 1.8e-15``,
    and bit-exact in **0/500**. Gating this tier at ``array_equal`` would be a
    false red on correct code.
    """
    rng = _rng(11)
    for _ in range(80):
        a, b, c = (_random_motion(rng, d) for _ in range(3))
        left, right = (a @ b) @ c, a @ (b @ c)
        _assert_close(left.linear, right.linear)
        _assert_close(left.translation, right.translation)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_a2_identity_acts_bit_exactly_on_both_sides(d: int) -> None:
    """A2 — ``e∘g = g∘e = g``, BIT-EXACT.

    MEASURED 500/500 bit-exact, so ``array_equal`` is earned rather than
    aspirational. Falsifier: an identity built at the wrong dimension, or one
    carrying any non-zero translation.
    """
    e = RigidMotion.identity(d)
    rng = _rng(2)
    for _ in range(40):
        g = _random_motion(rng, d)
        for label, product in (("e@g", e @ g), ("g@e", g @ e)):
            np.testing.assert_array_equal(
                product.linear, g.linear, err_msg=f"{label} moved Q"
            )
            np.testing.assert_array_equal(
                product.translation, g.translation, err_msg=f"{label} moved t"
            )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_a3_inverse_on_both_sides_and_against_the_homogeneous_matrix(d: int) -> None:
    r"""A3 — ``g∘g⁻¹ = g⁻¹∘g = e``, AND ``(Qᵀ, −Qᵀt)`` equals the closed-form
    inverse of the ``(d+1)×(d+1)`` homogeneous matrix (P4).

    The fixture MUST carry general rotations with a general ``t``: every
    reflection has ``Q = Qᵀ``, and for a *seated* reflection ``t ∈ span(n̂)``
    so ``−Qᵀt = t`` — a ``t ↦ −t`` bug is invisible on the whole reflection
    family (plan §7-V12).

    MEASURED: both sides residual 1.1e-15 / 1.8e-15, bit-exact **0/500** —
    the two sides are genuinely different computations.
    """
    e = RigidMotion.identity(d)
    rng = _rng(3)
    for _ in range(60):
        g = _random_motion(rng, d)
        for label, product in (
            ("g∘g⁻¹", g @ g.inverse()),
            ("g⁻¹∘g", g.inverse() @ g),
        ):
            _assert_close(
                product.linear, e.linear, err_msg=f"{label} linear part"
            )
            # The residual of g∘g⁻¹'s translation scales with ‖t‖, so the
            # bound is taken against the INPUT's scale, not the (zero) output.
            np.testing.assert_allclose(
                product.translation,
                e.translation,
                atol=_TOL_GENERAL * max(1.0, float(np.max(np.abs(g.translation)))),
                rtol=0.0,
                err_msg=f"{label} translation part",
            )
        # Independent construction: the homogeneous embedding, inverted by
        # a general-purpose solver that knows nothing about orthogonality.
        homogeneous = np.eye(d + 1)
        homogeneous[:d, :d] = g.linear
        homogeneous[:d, d] = g.translation
        reference = np.linalg.inv(homogeneous)
        _assert_close(g.inverse().linear, reference[:d, :d])
        _assert_close(g.inverse().translation, reference[:d, d])


@pytest.mark.parametrize("d", [2, 3])
def test_a4_inverse_is_an_anti_homomorphism(d: int) -> None:
    r"""A4 — ``(g∘h)⁻¹ = h⁻¹∘g⁻¹``, and the SWAPPED order is measurably wrong.

    The companion leg is what makes this non-vacuous: on an abelian group both
    orders agree and the law says nothing (plan §7-V3). MEASURED: the swapped
    order differs on 500/500 random draws at d ≥ 2.
    """
    rng = _rng(4)
    swapped_differs = 0
    for _ in range(60):
        g, h = _random_motion(rng, d), _random_motion(rng, d)
        correct = h.inverse() @ g.inverse()
        _assert_close((g @ h).inverse().linear, correct.linear)
        _assert_close((g @ h).inverse().translation, correct.translation)
        swapped = g.inverse() @ h.inverse()
        if not swapped.approximately_equals(correct, atol=1e-8):
            swapped_differs += 1
    assert swapped_differs > 0, (
        "the SWAPPED inverse order agreed on every draw — the fixture is "
        "effectively abelian and this gate is vacuous"
    )


def test_a5_generated_groups_have_the_textbook_order() -> None:
    r"""A5 — closure orders against PUBLISHED tables (P5).

    ``|⟨σ⟩| = 2``; ``|⟨σ_x,σ_y,σ_z⟩| = |D_2h| = 8``; ``|⟨C_4, σ_x⟩| = |C_4v| = 8``;
    ``|⟨C_n⟩| = n``; the hyperoctahedral generators close to ``2^d·d!`` — 48 at
    d = 3, which is ``|O_h|``.

    Falsifier: seeding the closure with the identity but never adding the
    generators themselves (orders collapse to 1), or dropping the ``+0.0``
    ``-0.0`` canonicalisation from the element key (orders inflate).
    """
    sx = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    sy = RigidMotion.reflection(normal=[0.0, 1.0, 0.0])
    sz = RigidMotion.reflection(normal=[0.0, 0.0, 1.0])
    quarter_turn = RigidMotion.rotation_from_circle_point(
        plane=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], point=(0.0, 1.0)
    )
    cases: list[tuple[str, list[RigidMotion], int]] = [
        ("<sigma> (d=1)", [RigidMotion.reflection(normal=[1.0])], 2),
        ("<sigma_z> (d=3)", [sz], 2),
        ("<sigma_x, sigma_y, sigma_z> = D_2h", [sx, sy, sz], 8),
        ("<C_4, sigma_x> = C_4v", [quarter_turn, sx], 8),
        ("hyperoctahedral d=2 = D_4", _hyperoctahedral(2), 8),
        ("hyperoctahedral d=3 = O_h", _hyperoctahedral(3), 48),
    ]
    for n in (2, 3, 4, 6):
        generator = RigidMotion.rotation(
            plane=[[1.0, 0.0], [0.0, 1.0]], angle=2.0 * np.pi / n
        )
        cases.append((f"<C_{n}> (d=2)", [generator], n))
    for label, generators, expected in cases:
        order = len(close_group(generators))
        assert order == expected, f"|{label}| = {order}, expected {expected}"
    assert len(close_group([], dimension=3)) == 1
    with pytest.raises(ValueError, match="explicit dimension"):
        close_group([])


def test_a6_a_generating_set_of_infinite_order_is_refused() -> None:
    r"""A6 — an infinite generated group RAISES, from the production entry point.

    Two configurations, both ordinary rather than exotic: a pure translation,
    and **two PARALLEL mirrors** — whose product is the pure translation by
    twice their separation, and which are the two reflective faces of a slab.

    The raising call is ``close_group`` itself, never a locally-constructed
    exception (plan §7-V8: a ``pytest.raises`` around a ``raise`` you wrote is
    green forever and pins nothing about production). Falsifier: remove the
    ``_MAX_GROUP_ORDER`` cap.
    """
    translation = RigidMotion.translation_by([1.0, 0.0, 0.0])
    with pytest.raises(NotAFinitePointGroupError, match="finite group"):
        close_group([translation])

    near = RigidMotion.reflection(normal=[1.0, 0.0, 0.0], offset=0.0)
    far = RigidMotion.reflection(normal=[1.0, 0.0, 0.0], offset=0.5)
    # The mechanism, asserted rather than assumed: the product of two parallel
    # mirrors IS the pure translation by twice their separation.
    product = far @ near
    _assert_close(product.linear, np.eye(3))
    _assert_close(product.translation, [1.0, 0.0, 0.0])
    assert product.element_order() is None, "a pure translation has infinite order"
    with pytest.raises(NotAFinitePointGroupError, match="no common seat"):
        close_group([near, far])


# ===========================================================================
# B — orthogonality and determinant
# ===========================================================================


def test_b1_orthogonality_and_determinant_are_symbolic_identities() -> None:
    r"""B1 — SymPy (P1), exact algebra with no floats in the chain.

    Five identities: the Householder's ``QᵀQ − I = 0`` and ``det = −1`` at
    d = 2 and d = 3, and Rodrigues' ``QᵀQ − I = 0``, ``det = +1``,
    ``trace = 1 + 2cos t``.

    The unit constraint MUST be imposed by an explicit parameterisation.
    MEASURED: substituting ``Σnᵢ² → 1`` after expansion does NOT fire and
    leaves ``4n₀²(Σnᵢ² − 1)`` residuals — the naive spelling silently proves
    nothing. Timings 0.1 / 0.3 / 0.5 / 0.8 / 5.0 s; kept in-file and NOT
    marked slow, because a ``slow``-marked foundation gate is deselected by
    the project's canonical ``-m "not slow"`` sweep and would stop guarding.
    """
    theta, phi, t = sp.symbols("theta phi t", real=True)

    normal_2d = sp.Matrix([sp.cos(phi), sp.sin(phi)])
    householder_2d = sp.eye(2) - 2 * normal_2d * normal_2d.T
    assert sp.simplify(householder_2d.T * householder_2d - sp.eye(2)) == sp.zeros(2, 2)
    assert sp.simplify(householder_2d.det()) == -1

    normal_3d = sp.Matrix(
        [sp.sin(theta) * sp.cos(phi), sp.sin(theta) * sp.sin(phi), sp.cos(theta)]
    )
    householder_3d = sp.eye(3) - 2 * normal_3d * normal_3d.T
    assert sp.simplify(householder_3d.T * householder_3d - sp.eye(3)) == sp.zeros(3, 3)
    assert sp.simplify(householder_3d.det()) == -1

    axis = normal_3d
    cross = sp.Matrix(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ]
    )
    rodrigues = (
        sp.eye(3) * sp.cos(t) + sp.sin(t) * cross + (1 - sp.cos(t)) * (axis * axis.T)
    )
    assert sp.simplify(rodrigues.det()) == 1
    assert sp.simplify(sp.trace(rodrigues) - (1 + 2 * sp.cos(t))) == 0
    assert sp.simplify(rodrigues.T * rodrigues - sp.eye(3)) == sp.zeros(3, 3)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_b2a_the_types_own_orthogonality_invariant(d: int) -> None:
    r"""B2a — the invariant the TYPE promises, and its rejection threshold.

    ``__post_init__`` accepts any ``Q`` with ``max|QᵀQ − I| ≤
    _ORTHOGONALITY_ATOL``, so an element carrying a sub-threshold shear is a
    **legal value of this type**. Asserting ``1e-14`` on an arbitrary element
    would assert a property ``RigidMotion`` does not promise — the
    constructors' far better quality is a separate, stronger claim (B2b).

    Both legs are mandatory (anti-pattern #11): the positive leg (a
    sub-threshold shear is ACCEPTED and satisfies the invariant) and the
    negative leg (a supra-threshold shear is REFUSED by the production
    constructor). Falsifier: delete the ``__post_init__`` guard.
    """
    atol = _txn._ORTHOGONALITY_ATOL
    if d >= 2:
        shear = np.eye(d)
        shear[0, 1] = 5.0e-13  # below the threshold -> a LEGAL element
        legal = RigidMotion(shear)
        departure = float(np.max(np.abs(legal.linear.T @ legal.linear - np.eye(d))))
        assert departure <= atol, departure
        assert departure > 1e-14, (
            "the 'legal shear' witness is inside 1e-14, so it does not "
            "demonstrate that the type is looser than the constructors"
        )
        too_sheared = np.eye(d)
        too_sheared[0, 1] = 2.0e-12  # above the threshold
        with pytest.raises(ValueError, match="not orthogonal"):
            RigidMotion(too_sheared)
    with pytest.raises(ValueError, match="not orthogonal"):
        RigidMotion(2.0 * np.eye(d))  # a pure scaling is not an isometry
    with pytest.raises(ValueError, match="square matrix"):
        RigidMotion(np.zeros((d, d + 1)))

    # Every product of legal elements is itself legal.
    rng = _rng(5)
    for _ in range(40):
        product = _random_motion(rng, d) @ _random_motion(rng, d)
        assert (
            float(np.max(np.abs(product.linear.T @ product.linear - np.eye(d)))) <= atol
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_b2b_the_named_constructors_are_far_better_than_the_invariant(d: int) -> None:
    """B2b — the constructors' ACTUAL quality: ``1e-14``, and exactly ``0`` for
    the integer-entried signed permutations.

    This is the row that catches a constructor regression; B2a cannot, because
    a sub-threshold shear is a legal element. Falsifier: build the Householder
    as ``I − n̂n̂ᵀ`` (missing the 2), or skip the normalisation of ``n̂``.
    """
    built: list[tuple[str, RigidMotion]] = [
        ("identity", RigidMotion.identity(d)),
        ("inversion", RigidMotion.inversion(d)),
        ("translation", RigidMotion.translation_by(np.arange(1.0, d + 1))),
    ]
    built += _reflection_family(d)
    if d >= 2:
        rng = _rng(6)
        for k in range(3):
            built.append((f"rotation[{k}]", _random_rotation(rng, d)))
    if d == 3:
        built.append(
            ("rotation_about_axis", RigidMotion.rotation_about_axis(axis=[1, 2, 3], angle=0.9))
        )
    for label, g in built:
        departure = float(np.max(np.abs(g.linear.T @ g.linear - np.eye(d))))
        assert departure <= 1e-14, f"{label}: max|QᵀQ − I| = {departure:.3e}"

    for g in _hyperoctahedral(d):
        np.testing.assert_array_equal(
            g.linear.T @ g.linear,
            np.eye(d),
            err_msg="signed permutations have integer entries — QᵀQ must be exactly I",
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_b3_determinant_values(d: int) -> None:
    r"""B3 — ``det`` per construction, dimension-generically.

    Reflection ``−1``; rotation ``+1``; translation ``+1``; and inversion
    ``(−1)^d``. The inversion row is why d = 1 and d = 3 must both be present:
    at d = 2 ``det(−I) = +1``, so that dimension is structurally BLIND to an
    inversion-versus-half-turn confusion.
    """
    for label, g in _reflection_family(d):
        assert g.determinant == pytest.approx(-1.0, abs=1e-14), label
        assert not g.is_proper, label
    assert RigidMotion.translation_by(np.ones(d)).determinant == pytest.approx(
        1.0, abs=1e-14
    )
    assert RigidMotion.inversion(d).determinant == pytest.approx(
        float((-1) ** d), abs=1e-14
    )
    if d >= 2:
        rng = _rng(7)
        for _ in range(10):
            rot = _random_rotation(rng, d)
            assert rot.determinant == pytest.approx(1.0, abs=1e-14)
            assert rot.is_proper


@pytest.mark.parametrize("d", [1, 2, 3])
def test_b4_det_is_a_homomorphism_and_seating_preserves_the_linear_part(d: int) -> None:
    r"""B4 — ``det(g∘h) = det(g)·det(h)``, and ``seated_at`` leaves ``Q``
    BIT-IDENTICAL.

    The homomorphism leg is VACUOUS on proper-only draws (``1 = 1·1``), so the
    fixture MUST mix parities: reflection×reflection (→ +1),
    reflection×rotation (→ −1), and ``−I`` (plan §7-V4). The assertion below
    checks that both signs actually occur.
    """
    pool: list[RigidMotion] = [g for _, g in _reflection_family(d)]
    pool.append(RigidMotion.inversion(d))
    pool.append(RigidMotion.identity(d))
    if d >= 2:
        rng = _rng(8)
        pool += [_random_rotation(rng, d) for _ in range(3)]

    seen_signs: set[int] = set()
    for g, h in itertools.product(pool, repeat=2):
        product = g @ h
        assert product.determinant == pytest.approx(
            g.determinant * h.determinant, abs=1e-12
        )
        seen_signs.add(int(round(product.determinant)))
    assert seen_signs == {-1, 1}, (
        f"only {seen_signs} occurred — a proper-only fixture makes det(gh) = "
        "det(g)det(h) the vacuous identity 1 = 1·1"
    )

    centre = np.arange(1.0, d + 1) * 0.37
    for g in pool:
        np.testing.assert_array_equal(
            g.seated_at(centre).linear,
            g.linear,
            err_msg="seating must not touch the linear part",
        )
        assert g.seated_at(centre).determinant == pytest.approx(
            g.determinant, abs=1e-14
        )


# ===========================================================================
# C — element order, and the fixed set
# ===========================================================================


@pytest.mark.parametrize("d", [1, 2, 3])
def test_c1_a_reflection_is_an_involution(d: int) -> None:
    r"""C1 — ``σ² = e``, two tiers, WITH ITS BLINDNESS DECLARED.

    ⚠ **This gate is Mode-12 designed-green on the affine offset and MUST
    NEVER be credited for it.** ``(σ, t)² = (I, σt + t)``, and ``σt = −t`` for
    *every* ``t ∈ span(n̂)``, so every offset whatsoever is an involution.
    MEASURED: ``t ∈ {d n̂, 2d n̂, 4d n̂, −2d n̂}`` are all involutions while the
    true mirror plane moves by 0.37 / 0.74 / 1.48. The catcher for the offset
    is **C3**, not this row.

    Two tiers, MEASURED: ``σ² == I`` is bit-exact for a coordinate normal and
    **0/300** for a general one.
    """
    e = RigidMotion.identity(d)
    for label, sigma in _reflection_family(d):
        assert sigma.element_order() == 2, label
        squared = sigma @ sigma
        np.testing.assert_allclose(
            squared.linear, e.linear, atol=1e-15, rtol=0.0, err_msg=label
        )
        seated = sigma.seated_at(np.arange(1.0, d + 1) * 0.21)
        assert seated.element_order() == 2, f"{label} seated"

    coordinate = RigidMotion.reflection(normal=np.eye(d)[0])
    np.testing.assert_array_equal(
        (coordinate @ coordinate).linear,
        np.eye(d),
        err_msg="a coordinate-normal Householder is a signed diagonal — σ² is exact",
    )


@pytest.mark.parametrize("n", [2, 3, 4, 6, 8, 12])
def test_c2_the_cyclic_generator_has_order_exactly_n(n: int) -> None:
    r"""C2 — ``ord(R(2π/n)) = n`` EXACTLY, and ``ord(R(2πk/n)) = n/gcd(k,n)``.

    ``g^n = e`` alone is satisfied by every element whose order DIVIDES ``n``,
    so only the smallest such power pins the element.

    The parametrisation MUST include even ``n``: MEASURED, the ``k = 2``
    generator mutation gives ``n/gcd(2,n)`` — caught at n = 4, 6, 8, 12 (the
    order halves) and NOT at n = 3 (the order stays 3). An odd-only fixture is
    blind (plan §7-V2).

    Not tolerance-fragile: MEASURED at n = 24, ``‖Qⁿ − I‖ ≈ 1e-15`` against
    ``min_{k<n}‖Q^k − I‖ ≥ 1.6e-1`` — six orders of margin.
    """
    plane = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    generator = RigidMotion.rotation(plane=plane, angle=2.0 * np.pi / n)
    assert generator.element_order() == n

    for k in range(1, n + 1):
        power = RigidMotion.rotation(plane=plane, angle=2.0 * np.pi * k / n)
        assert power.element_order() == n // math.gcd(k, n), f"k={k}"


@pytest.mark.parametrize("d", [1, 2, 3])
def test_c3_an_element_fixes_its_named_fixed_set_pointwise(d: int) -> None:
    r"""C3 — ⭐ the ONLY catcher for the affine offset.

    ``σ_{n̂,offset}`` must fix every point of ``{x : n̂·x = offset}``; a seated
    rotation must fix its centre; a pure translation must fix nothing.

    Falsifier — the offset factor: ``t = d n̂``, ``t = 4d n̂``, ``t = −2d n̂``.
    MEASURED, all three remain involutions (C1 stays green) and fail this gate
    by 0.37 / 0.74 / 1.48 respectively. The correct ``t = 2d n̂`` gives residual
    exactly 0.0.
    """
    rng = _rng(9)
    for label, raw_normal in _NORMALS[d]:
        normal = np.asarray(raw_normal, dtype=float)
        normal /= np.linalg.norm(normal)
        offset = 0.37
        sigma = RigidMotion.reflection(normal=raw_normal, offset=offset)
        # Points ON the mirror: project a random cloud onto {n̂·x = offset}.
        cloud = rng.normal(size=(24, d))
        on_plane = cloud - np.outer(cloud @ normal - offset, normal)
        np.testing.assert_allclose(
            on_plane @ normal, offset, atol=1e-14, rtol=0.0,
            err_msg="fixture error: the cloud is not on the plane",
        )
        assert sigma.fixes(on_plane, atol=1e-13).all(), (
            f"{label}: the seated reflection does not fix its own mirror plane"
        )
        # And it moves points OFF the plane by twice their signed distance.
        off_plane = on_plane + np.outer(np.full(len(on_plane), 0.5), normal)
        np.testing.assert_allclose(
            sigma.on_points(off_plane),
            on_plane - np.outer(np.full(len(on_plane), 0.5), normal),
            atol=1e-13,
            rtol=0.0,
        )

    centre = np.arange(1.0, d + 1) * 0.6
    if d >= 2:
        seated_rotation = _random_rotation(_rng(10), d).seated_at(centre)
        assert bool(seated_rotation.fixes(centre, atol=_TOL_GENERAL))
    assert not bool(
        RigidMotion.translation_by(np.ones(d)).fixes(centre, atol=1e-6)
    ), "a pure translation fixes nothing"


@pytest.mark.parametrize("d", [1, 2, 3])
def test_c4_fixed_subspace_dimension_is_the_structural_claim(d: int) -> None:
    r"""C4 — ``dim Fix(σ) = d−1`` and ``dim Fix(R) = d−2`` (ruling R2).

    This is what makes reflections and rotations dimension-generic, and it is
    read off ``ker(Q − I)`` — the LINEAR part, so seating cannot change it.

    Angles are bounded away from zero because the shipped implementation uses
    an ABSOLUTE singular-value threshold; the degenerate end is pinned
    separately as an honest-scope control below.
    """
    centre = np.arange(1.0, d + 1) * 0.4
    for label, sigma in _reflection_family(d):
        assert sigma.fixed_subspace_dimension == d - 1, label
        assert sigma.seated_at(centre).fixed_subspace_dimension == d - 1, label
    if d >= 2:
        plane = _orthonormal_plane(_rng(12), d)
        for angle in _ROTATION_ANGLES:
            rot = RigidMotion.rotation(plane=plane, angle=angle)
            assert rot.fixed_subspace_dimension == d - 2, angle
            assert rot.seated_at(centre).fixed_subspace_dimension == d - 2, angle
    assert RigidMotion.identity(d).fixed_subspace_dimension == d
    assert RigidMotion.translation_by(np.ones(d)).fixed_subspace_dimension == d


def test_c4_degenerate_control_a_tiny_rotation_reads_as_the_identity() -> None:
    r"""C4 control — the documented limit of ``fixed_subspace_dimension``.

    ``matrix_rank(Q − I, tol=1e-10)`` cannot distinguish a rotation by an
    angle small compared with that threshold from the identity. That is a
    property of the numerics, not a bug — pinned here so a future reader
    cannot over-claim the gate above, and so the threshold cannot silently
    move without something reddening.
    """
    plane = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    assert RigidMotion.rotation(plane=plane, angle=1.0e-6).fixed_subspace_dimension == 1
    assert RigidMotion.rotation(plane=plane, angle=1.0e-13).fixed_subspace_dimension == 3


def test_c5_cartan_dieudonne_sharp_form() -> None:
    r"""C5 — ``k = d − dim Fix(Q)`` reflections suffice, and ``det = (−1)^k``.

    The sharp form of Cartan–Dieudonné, and the row that ties B3 to C4:
    neither ``det`` nor ``dim Fix`` can be wrong alone without this reddening.
    MEASURED on all eight cases below, including ``−I`` at d = 1, 2, 3 (where
    ``k`` = 1, 2, 3 and ``det`` = −1, +1, −1), a rotoreflection, and the
    identity (``k = 0``).
    """
    cases: list[tuple[str, RigidMotion]] = [
        ("d=1 reflection", RigidMotion.reflection(normal=[1.0])),
        ("d=1 -I", RigidMotion.inversion(1)),
        ("d=2 rotation", RigidMotion.rotation(plane=[[1.0, 0.0], [0.0, 1.0]], angle=0.6)),
        ("d=2 -I", RigidMotion.inversion(2)),
        ("d=3 reflection", RigidMotion.reflection(normal=[1.0, 2.0, 3.0])),
        ("d=3 rotation", RigidMotion.rotation_about_axis(axis=[1, 1, 1], angle=0.9)),
        ("d=3 -I", RigidMotion.inversion(3)),
        (
            "d=3 rotoreflection",
            RigidMotion.reflection(normal=[0.0, 0.0, 1.0])
            @ RigidMotion.rotation_about_axis(axis=[0.0, 0.0, 1.0], angle=0.7),
        ),
        ("d=3 identity", RigidMotion.identity(3)),
    ]
    for label, g in cases:
        k = g.dimension - g.fixed_subspace_dimension
        assert g.determinant == pytest.approx(float((-1) ** k), abs=1e-12), (
            f"{label}: dim Fix = {g.fixed_subspace_dimension}, k = {k}, "
            f"det = {g.determinant:+.6f}, expected {(-1) ** k:+d}"
        )


# ===========================================================================
# D — independent constructions cross-checked
# ===========================================================================


@pytest.mark.parametrize("d", [2, 3])
def test_d1_householder_equals_the_projector_construction(d: int) -> None:
    r"""D1 — ``I − 2n̂n̂ᵀ`` equals the SPECTRAL form built in the test (P4).

    The reference is assembled from the definition — ``+1`` on an orthonormal
    basis of the hyperplane, ``−1`` on the normal — via a QR the SUT never
    calls. One is a rank-1 update, the other a spectral sum: structurally
    independent constructions of the same operator.
    """
    for label, raw in _NORMALS[d]:
        normal = np.asarray(raw, dtype=float)
        normal /= np.linalg.norm(normal)
        basis, _ = np.linalg.qr(np.column_stack([normal, np.eye(d)[:, : d - 1]]))
        # Re-align the first column with n̂ (QR may flip its sign).
        basis = basis * np.sign(basis[:, 0] @ normal if d else 1.0)
        tangent = basis[:, 1:]
        reference = tangent @ tangent.T - np.outer(normal, normal)
        np.testing.assert_allclose(
            RigidMotion.reflection(normal=raw).linear,
            reference,
            atol=1e-14,
            rtol=0.0,
            err_msg=label,
        )


def test_d2_rotation_about_axis_matches_scipy_from_rotvec() -> None:
    r"""D2 — cross-check against ``scipy…Rotation.from_rotvec`` (P2).

    An independently maintained, quaternion-based implementation: a different
    ALGORITHM, not a re-spelling of the same closed form. MEASURED max abs
    ``1.11e-15`` over 300 random (axis, angle) draws.

    Gated on ``atol``, NEVER on ULP: the same data reads "5120 ULP" because
    the matrix entries cross zero (plan §7-V7).
    """
    rng = _rng(13)
    worst = 0.0
    for _ in range(300):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = rng.uniform(-np.pi, np.pi)
        ours = RigidMotion.rotation_about_axis(axis=axis, angle=angle).linear
        theirs = ScipyRotation.from_rotvec(axis * angle).as_matrix()
        worst = max(worst, float(np.max(np.abs(ours - theirs))))
    assert worst <= _TOL_SCIPY_ROTVEC, f"max |ours − scipy| = {worst:.3e}"


@pytest.mark.parametrize("d", [2, 3, 4])
def test_d3_planar_rotation_matches_the_lie_exponential(d: int) -> None:
    r"""D3 — cross-check against ``expm(θ(vuᵀ − uvᵀ))`` (P3).

    The Lie-theoretic definition ``R = exp(θJ)``, evaluated by Padé
    scaling-and-squaring — a THIRD structurally independent route, and the
    only DIMENSION-GENERIC one, which is why d = 4 appears here even though
    no consumer needs it: it is the cheapest evidence that the ``(d−2)``
    formula is not a d ≤ 3 coincidence.

    MEASURED 3.9e-14 against the shipped constructor; a 1e-15 gate would red
    on Padé noise.
    """
    rng = _rng(14 + d)
    worst = 0.0
    for _ in range(100):
        plane = _orthonormal_plane(rng, d)
        angle = rng.uniform(-np.pi, np.pi)
        ours = RigidMotion.rotation(plane=plane, angle=angle).linear
        u, v = plane[0], plane[1]
        theirs = expm(angle * (np.outer(v, u) - np.outer(u, v)))
        worst = max(worst, float(np.max(np.abs(ours - theirs))))
    assert worst <= _TOL_EXPM, f"d={d}: max |ours − expm| = {worst:.3e}"


def test_d4_cartan_dieudonne_composition_symbolic_and_numeric() -> None:
    r"""D4 — ``σ_α ∘ σ_β = R(2(α−β))``: two reflections whose planes meet at
    ``θ/2`` compose to a rotation by ``θ``.

    Symbolic at d = 2 (P1 — MEASURED residual exactly 0 in 0.3 s) and numeric
    at d = 3 with NON-COORDINATE normals, so the row is not secretly a d = 2
    problem embedded in d = 3.

    This is the ONLY gate that pins the rotation's SENSE against an
    independent construction: swapping ``σ_α ∘ σ_β`` for ``σ_β ∘ σ_α`` gives
    the rotation by ``−θ``, and the companion leg below asserts the two
    genuinely differ.
    """
    alpha, beta = sp.symbols("alpha beta", real=True)

    def _householder_2d(angle: sp.Symbol) -> sp.Matrix:
        n = sp.Matrix([sp.cos(angle), sp.sin(angle)])
        return sp.eye(2) - 2 * n * n.T

    product = sp.simplify(_householder_2d(alpha) * _householder_2d(beta))
    turn = 2 * (alpha - beta)
    reference = sp.Matrix(
        [[sp.cos(turn), -sp.sin(turn)], [sp.sin(turn), sp.cos(turn)]]
    )
    assert sp.simplify(sp.expand_trig(product - reference)) == sp.zeros(2, 2)

    # Numeric, d = 3, normals off the coordinate planes.
    for theta in (0.3, 1.1, 2.0):
        axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
        # Two normals in the plane perpendicular to `axis`, separated by θ/2.
        u = np.cross(axis, [1.0, 0.0, 0.0])
        u /= np.linalg.norm(u)
        v = np.cross(axis, u)
        n_beta = u
        n_alpha = np.cos(theta / 2.0) * u + np.sin(theta / 2.0) * v
        composed = RigidMotion.reflection(normal=n_alpha) @ RigidMotion.reflection(
            normal=n_beta
        )
        reference_rotation = RigidMotion.rotation(plane=(u, v), angle=theta)
        np.testing.assert_allclose(
            composed.linear, reference_rotation.linear, atol=1e-14, rtol=0.0
        )
        assert composed.determinant == pytest.approx(1.0, abs=1e-14)
        swapped = RigidMotion.reflection(normal=n_beta) @ RigidMotion.reflection(
            normal=n_alpha
        )
        assert not swapped.approximately_equals(composed, atol=1e-8), (
            "the reversed reflection product agreed — this fixture cannot see "
            "the rotation's sense"
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_d5_inversion_is_a_product_of_d_orthogonal_reflections(d: int) -> None:
    r"""D5 — ``−I = σ_{e_1} ∘ … ∘ σ_{e_d}``, BIT-EXACT; and at d = 3 it is
    NEITHER a reflection NOR a rotation.

    MEASURED bit-exact (the coordinate Householders are signed diagonals).
    The second leg is the API gate: ``−I`` at d = 3 has ``dim Fix = 0``, which
    is neither ``d−1 = 2`` (reflection) nor ``d−2 = 1`` (rotation), so it
    cannot be typed as either — which is exactly why the shipped design has
    ONE ``RigidMotion`` with a COMPUTED determinant rather than a
    ``Reflection``/``Rotation`` type pair.
    """
    product = RigidMotion.identity(d)
    for axis in range(d):
        product = product @ RigidMotion.reflection(normal=np.eye(d)[axis])
    np.testing.assert_array_equal(product.linear, -np.eye(d))
    np.testing.assert_array_equal(
        product.linear, RigidMotion.inversion(d).linear
    )
    inversion = RigidMotion.inversion(d)
    assert inversion.fixed_subspace_dimension == 0
    if d == 3:
        assert inversion.fixed_subspace_dimension not in (d - 1, d - 2)
        assert inversion.determinant == pytest.approx(-1.0, abs=1e-14)


# ===========================================================================
# E — conjugation ("seat it elsewhere" is "re-orient it")
# ===========================================================================


def test_e1_conjugating_a_reflection_rotates_its_normal() -> None:
    r"""E1 — ``R σ_n̂ R⁻¹ = σ_{R n̂}``.

    Symbolic for ``R = R_z(t)`` against an in-plane normal ``n(φ)`` (P1 —
    MEASURED residual exactly 0 in 0.6 s, and it yields the bonus identity
    ``R_z(t) n(φ) = n(φ + t)``), plus numeric for a general ``R`` (MEASURED
    1.2e-15 over 500 draws).

    ⚠ The FULLY symbolic version — a general Rodrigues ``R`` against a general
    Householder — CHOKES SymPy (MEASURED: killed at 120 s). Specialising one
    factor closes in 0.6 s and proves the same law shape; the general case is
    covered numerically.

    Falsifier: conjugate the wrong way (``R⁻¹σR``) — the normal rotates by
    ``−t``.
    """
    t, phi = sp.symbols("t phi", real=True)
    rot_z = sp.Matrix(
        [[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]]
    )
    normal = sp.Matrix([sp.cos(phi), sp.sin(phi), 0])
    sigma = sp.eye(3) - 2 * normal * normal.T
    conjugated = (rot_z * sigma * rot_z.T).applyfunc(
        lambda entry: sp.simplify(sp.expand_trig(sp.expand(entry)))
    )
    moved_normal = (rot_z * normal).applyfunc(sp.simplify)
    reference = (sp.eye(3) - 2 * moved_normal * moved_normal.T).applyfunc(sp.simplify)
    assert (conjugated - reference).applyfunc(
        lambda entry: sp.simplify(sp.expand_trig(sp.expand(entry)))
    ) == sp.zeros(3, 3)
    # The bonus identity the specialisation exposes: the mirror plane rotates.
    assert sp.simplify(sp.expand_trig(moved_normal[0] - sp.cos(phi + t))) == 0
    assert sp.simplify(sp.expand_trig(moved_normal[1] - sp.sin(phi + t))) == 0

    rng = _rng(15)
    for _ in range(200):
        raw = rng.normal(size=3)
        rotation = _random_rotation(rng, 3)
        conjugate = RigidMotion.reflection(normal=raw).conjugated_by(rotation)
        expected = RigidMotion.reflection(normal=rotation.on_directions(raw))
        np.testing.assert_allclose(
            conjugate.linear, expected.linear, atol=_TOL_GENERAL, rtol=0.0
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_e2_affine_conjugation_moves_the_mirror_plane(d: int) -> None:
    r"""E2 — ⭐ ``g ∘ σ_{n̂,δ} ∘ g⁻¹ = σ_{Qn̂,\; δ + (Qn̂)·t}`` for ANY rigid
    motion ``g = (Q, t)``.

    This is the law that makes "seat this mirror somewhere else" and
    "re-orient it" ONE operation. MEASURED residual 2.7e-15 over 300 draws.

    Falsifier: drop the ``(Qn̂)·t`` term — i.e. conjugate only the linear part
    — and the mirror lands at the wrong offset. Without this row the affine
    part is unverified under the group action, and the companion leg below
    asserts the dropped-term spelling genuinely differs.
    """
    rng = _rng(16)
    dropped_term_differs = 0
    for _ in range(120):
        raw = rng.normal(size=d)
        normal = raw / np.linalg.norm(raw)
        offset = float(rng.normal())
        g = _random_motion(rng, d)
        conjugate = RigidMotion.reflection(normal=raw, offset=offset).conjugated_by(g)

        moved_normal = g.on_directions(normal)
        moved_offset = offset + float(moved_normal @ g.translation)
        expected = RigidMotion.reflection(normal=moved_normal, offset=moved_offset)
        _assert_close(conjugate.linear, expected.linear)
        np.testing.assert_allclose(
            conjugate.translation, expected.translation, atol=_TOL_GENERAL, rtol=0.0
        )
        naive = RigidMotion.reflection(normal=moved_normal, offset=offset)
        if not naive.approximately_equals(expected, atol=1e-8):
            dropped_term_differs += 1
    assert dropped_term_differs > 0, (
        "dropping the (Qn̂)·t term changed nothing on any draw — this fixture "
        "cannot see the affine half of conjugation"
    )


@pytest.mark.parametrize("d", [2, 3])
def test_e3_conjugation_is_a_class_function(d: int) -> None:
    """E3 — conjugation preserves ``element_order``, ``determinant`` and
    ``fixed_subspace_dimension``.

    Any implementation of ``conjugated_by`` that is not a group automorphism
    breaks at least one of the three.
    """
    rng = _rng(17)
    conjugators = [_random_motion(rng, d) for _ in range(4)]
    subjects: list[RigidMotion] = [g for _, g in _reflection_family(d)]
    subjects.append(RigidMotion.inversion(d))
    subjects.append(
        RigidMotion.rotation(plane=_orthonormal_plane(rng, d), angle=2 * np.pi / 6)
    )
    for h in conjugators:
        for g in subjects:
            moved = g.conjugated_by(h)
            assert moved.element_order() == g.element_order()
            assert moved.determinant == pytest.approx(g.determinant, abs=1e-12)
            assert moved.fixed_subspace_dimension == g.fixed_subspace_dimension


# ===========================================================================
# F — the seating relation (ruling R1)
# ===========================================================================


@pytest.mark.parametrize("d", [1, 2, 3])
def test_f1_seating_puts_the_centre_in_the_fixed_set(d: int) -> None:
    """F1 — ``seat(Q, c) = (Q, (I−Q)c)`` fixes ``c``.

    MEASURED residual 8.9e-16, bit-exact only 58/500 — a tolerance is
    required. Falsifier: ``t = (Q−I)c`` (sign flip) gives ``g(c) = 2Qc − c``.
    """
    rng = _rng(18)
    for _ in range(60):
        linear = _random_linear(rng, d)
        centre = rng.normal(size=d)
        seated = linear.seated_at(centre)
        tol = _TOL_GENERAL * max(1.0, float(np.max(np.abs(centre))))
        assert bool(seated.fixes(centre, atol=tol)), (
            f"max|g(c) − c| = {np.max(np.abs(seated.on_points(centre) - centre)):.3e}"
        )
        np.testing.assert_array_equal(seated.linear, linear.linear)

    # HONEST-SCOPE CONTROL: `seated_at` is the general conjugation, so on an
    # element that ALREADY carries a translation it does NOT fix its argument
    # and is not idempotent. Pinned so nobody reads the gate above as
    # "seated_at always fixes its centre" — MEASURED, a glide seated at a
    # random c misses it by O(1) (4.1 on a d = 1 draw).
    glide = RigidMotion(
        _random_linear(rng, d).linear, np.full(d, 1.5)
    )
    stray = np.full(d, 0.75)
    if not np.allclose(glide.linear, np.eye(d)):
        assert not bool(glide.seated_at(stray).fixes(stray, atol=1e-8)), (
            "seated_at is documented as the general conjugation — an element "
            "already carrying a translation must NOT be re-seated into fixing c"
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_f2_every_row_of_the_seating_table(d: int) -> None:
    r"""F2 — ⭐ each row of ruling R1's table, recomputed from ``(I−Q)c``.

    ==============  =======================  =========  ===============
    element         geometric data           ``Q``      ``t``
    ==============  =======================  =========  ===============
    reflection      normal ``n̂``, offset δ   I − 2n̂n̂ᵀ   ``2δ n̂``
    rotation        plane, angle, centre c   R          ``(I − R)c``
    inversion       centre ``c``             −I         ``2c``
    translation     vector ``v``             I          ``v``
    ==============  =======================  =========  ===============

    The reflection and inversion rows are the two carrying a FACTOR OF 2 —
    exactly the factor C1 is structurally blind to (see C3).
    """
    rng = _rng(19)
    raw = rng.normal(size=d)
    normal = raw / np.linalg.norm(raw)
    offset = 0.37
    centre = rng.normal(size=d)

    reflection = RigidMotion.reflection(normal=raw, offset=offset)
    np.testing.assert_allclose(
        reflection.translation, 2.0 * offset * normal, atol=1e-15, rtol=0.0
    )
    np.testing.assert_allclose(
        reflection.linear, np.eye(d) - 2.0 * np.outer(normal, normal),
        atol=1e-15, rtol=0.0,
    )

    inversion = RigidMotion.inversion(d).seated_at(centre)
    _assert_close(inversion.translation, 2.0 * centre)
    np.testing.assert_array_equal(inversion.linear, -np.eye(d))

    translation = RigidMotion.translation_by(centre)
    np.testing.assert_array_equal(translation.linear, np.eye(d))
    np.testing.assert_array_equal(translation.translation, centre)

    if d >= 2:
        rotation = _random_rotation(rng, d)
        seated = rotation.seated_at(centre)
        _assert_close(seated.translation, (np.eye(d) - rotation.linear) @ centre)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_f3_the_offset_has_one_source_of_truth(d: int) -> None:
    r"""F3 — the FORMULA and the ROUTING, gated separately and deliberately.

    Leg 1 pins the formula: ``reflection(n̂, δ).translation`` is BIT-IDENTICAL
    to ``2δn̂`` recomputed here. Leg 2 pins the routing: the same element
    reached via ``reflection(n̂).seated_at(δn̂)`` agrees to 1.11e-16
    (MEASURED — the two are genuinely different code paths, ``2δn̂`` directly
    versus ``(I−Q)c``).

    They MUST NOT be merged. If ``reflection`` is ever refactored to delegate
    to ``seated_at``, leg 2 silently becomes ``X == X`` while leg 1 keeps
    pinning the arithmetic.
    """
    rng = _rng(20)
    raw = rng.normal(size=d)
    normal = raw / np.linalg.norm(raw)
    for offset in (0.0, 0.37, -1.5):
        direct = RigidMotion.reflection(normal=normal, offset=offset)
        np.testing.assert_array_equal(
            direct.translation,
            2.0 * offset * normal,
            err_msg="the offset formula t = 2 δ n̂ is not exact",
        )
        via_seating = RigidMotion.reflection(normal=normal).seated_at(offset * normal)
        np.testing.assert_allclose(
            via_seating.translation, direct.translation, atol=1e-15, rtol=0.0
        )


@pytest.mark.parametrize("d", [1, 2, 3])
def test_f4_seating_at_the_origin_is_the_linear_element(d: int) -> None:
    """F4 — ``seat(Q, 0) = (Q, 0)`` bit-exactly, with no ``-0.0`` leaking in.

    ``-0.0`` compares equal to ``0.0`` but hashes differently unless
    canonicalised, so the sign bit is asserted explicitly: this is what keeps
    ``is_linear`` and the closure key from disagreeing about the same element.
    """
    rng = _rng(21)
    pool: list[RigidMotion] = [g for _, g in _reflection_family(d)]
    pool.append(RigidMotion.inversion(d))
    if d >= 2:
        pool.append(_random_rotation(rng, d))
    for g in pool:
        seated = g.seated_at(np.zeros(d))
        np.testing.assert_array_equal(seated.linear, g.linear)
        np.testing.assert_array_equal(seated.translation, np.zeros(d))
        assert seated.is_linear
        assert not np.signbit(seated.translation).any(), "a -0.0 leaked into t"


# ===========================================================================
# O — permutes / preserves, and the orbit action (ruling R4)
# ===========================================================================


def _point_set_cases() -> list[tuple[str, np.ndarray, list[RigidMotion]]]:
    """(label, points, group elements) triples that are genuinely invariant."""
    return [
        ("cube vertices / O_h", _cube_vertices(), _hyperoctahedral(3)),
        ("cube face centres / O_h", _cube_face_centres(), _hyperoctahedral(3)),
        ("cube edge midpoints / O_h", _cube_edge_midpoints(), _hyperoctahedral(3)),
        ("square vertices / D_4", _square_vertices(), _hyperoctahedral(2)),
        (
            "centred 1-D lattice n=6",
            _cell_centres(6, origin=-0.3),
            [RigidMotion.identity(1), RigidMotion.reflection(normal=[1.0])],
        ),
        (
            "centred 1-D lattice n=7",
            _cell_centres(7, origin=-0.35),
            [RigidMotion.identity(1), RigidMotion.reflection(normal=[1.0])],
        ),
    ]


def test_o1_permutes_returns_a_genuine_permutation() -> None:
    r"""O1 — POSITIVE: an invariant set yields a real ``Permutation``.

    The positive leg is what makes O3/O4 mean anything (anti-pattern #11): a
    negative-only suite validates the *refusing*, not the claim. Here the
    returned object is re-checked independently — ``sort(π) == arange(n)`` and
    ``‖g(x) − x[π]‖ ≤ atol`` — and the 1-D rows pin the closed form
    ``π(σ) == arange(n)[::-1]`` EXACTLY (P4/P6).
    """
    for label, points, elements in _point_set_cases():
        for g in elements:
            pi = g.permutes(points, atol=_TOL_MATCH)
            assert pi is not None, f"{label}: an invariant set was refused"
            assert isinstance(pi, Permutation)
            np.testing.assert_array_equal(
                np.sort(pi.indices), np.arange(len(points)), err_msg=label
            )
            np.testing.assert_allclose(
                g.on_points(points), points[pi.indices], atol=_TOL_MATCH, rtol=0.0,
                err_msg=label,
            )

    for n in (6, 7):
        lattice = _cell_centres(n, origin=-0.05 * n)
        mirror = RigidMotion.reflection(normal=[1.0])
        pi = mirror.permutes(lattice, atol=_TOL_MATCH)
        assert pi is not None
        np.testing.assert_array_equal(pi.indices, np.arange(n)[::-1])


@pytest.mark.parametrize(
    "label,points,elements",
    [
        ("cube / O_h (non-abelian)", _cube_vertices(), _hyperoctahedral(3)),
        ("square / D_4 (non-abelian)", _square_vertices(), _hyperoctahedral(2)),
    ],
)
def test_o2_the_induced_permutation_is_a_group_homomorphism(
    label: str, points: np.ndarray, elements: list[RigidMotion]
) -> None:
    r"""O2 — ⭐⭐ ``π_{g∘h} = π_g ∘ π_h``, plus ``π_e = id`` and
    ``π_{g⁻¹} = π_g⁻¹``.

    The deepest cheap law available on this object: it cross-checks TWO
    independent layers — the element algebra (matrix composition) and the
    action (nearest-neighbour matching on points) — using only exact integer
    arithmetic, no tolerance possible.

    It simultaneously pins the composition order, the row-versus-column action
    convention, and the ``π``-versus-``π⁻¹`` choice. **None of the three is
    checkable alone.**

    VACUOUS ON AN ABELIAN GROUP, so both fixtures are non-abelian and the
    companion leg asserts the reversed order genuinely disagrees. MEASURED:
    the reversed order violates 102/144 pairs on the cube and 42/144 on a
    12-element sample; on the abelian ``C_4`` it agrees 16/16.
    """
    subset = elements[:12] if len(elements) > 12 else elements
    reversed_disagreements = 0
    for a, b in itertools.product(subset, repeat=2):
        pi_a = a.permutes(points, atol=_TOL_MATCH)
        pi_b = b.permutes(points, atol=_TOL_MATCH)
        pi_ab = (a @ b).permutes(points, atol=_TOL_MATCH)
        assert pi_a is not None and pi_b is not None and pi_ab is not None
        assert pi_ab == pi_a @ pi_b, (
            f"{label}: π(a∘b) != π(a)∘π(b) — "
            f"{pi_ab.indices.tolist()} vs {(pi_a @ pi_b).indices.tolist()}"
        )
        if pi_ab != pi_b @ pi_a:
            reversed_disagreements += 1
    assert reversed_disagreements > 0, (
        f"{label}: the REVERSED composition order agreed on every pair — the "
        "fixture is abelian and this gate is vacuous"
    )

    n = len(points)
    identity = RigidMotion.identity(points.shape[1])
    assert identity.permutes(points, atol=_TOL_MATCH) == Permutation.identity(n)
    for g in subset:
        pi = g.permutes(points, atol=_TOL_MATCH)
        pi_inv = g.inverse().permutes(points, atol=_TOL_MATCH)
        assert pi is not None and pi_inv is not None
        assert pi_inv == pi.inverse()
        assert pi @ pi.inverse() == Permutation.identity(n)


def test_o3_a_non_invariant_set_is_refused() -> None:
    """O3 — NEGATIVE-1: a genuinely asymmetric set gets ``None``.

    The image simply is not a point of the set. Falsifier: remove the
    distance-window guard.
    """
    asymmetric = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.3, 0.4, 0.9]])
    sigma = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    assert sigma.permutes(asymmetric, atol=_TOL_MATCH) is None
    assert (
        sigma.preserves(
            asymmetric, np.ones(3), atol=_TOL_MATCH, weight_atol=_TOL_MATCH
        )
        is None
    )
    # …and a set that IS invariant is not refused, so the guard is not simply
    # returning None for everything.
    assert sigma.permutes(_cube_vertices(), atol=_TOL_MATCH) is not None


@pytest.mark.catches("ERR-073")
def test_o4_a_duplicated_point_is_refused_and_guard_two_is_the_one_that_fires() -> None:
    r"""O4 — NEGATIVE-2 (ERR-073), WITH GUARD ISOLATION.

    Appending a bit-identical duplicate of one cube vertex makes every match
    map many-to-one: the duplicated position carries twice the mass of its
    image, so the set is NOT invariant — yet an injectivity-free check
    certifies it, which is exactly how ERR-073 shipped.

    ``permutes`` returns ``None`` for three different reasons, so ``is None``
    alone proves only that *one of* the guards fired. The isolation is
    therefore established from the INPUT side, in plain numpy: MEASURED over
    all 48 ``O_h`` operators, guard 1 (image-is-a-point) passes 48/48, guard 3
    (weights) passes 48/48, and guard 2 (bijection) passes **0/48**.

    Note for the mutation record: deleting the ``np.unique(π).size == n``
    check makes the ``Permutation`` CONSTRUCTOR raise rather than the function
    return ``None``, so this test reds with a ``ValueError`` instead of an
    ``AssertionError``. Still a red — but for a substituted mechanism.
    """
    vertices = _cube_vertices()
    duplicated = np.vstack([vertices, vertices[0]])
    weights = np.full(len(duplicated), 0.125)
    n = len(duplicated)

    guard_one = guard_two = guard_three = 0
    for g in _hyperoctahedral(3):
        assert g.permutes(duplicated, atol=_TOL_MATCH) is None
        assert (
            g.preserves(
                duplicated, weights, atol=_TOL_MATCH, weight_atol=_TOL_MATCH
            )
            is None
        )
        # Input-side isolation, computed here rather than trusted.
        moved = g.on_points(duplicated)
        dist = np.linalg.norm(moved[:, None, :] - duplicated[None, :, :], axis=2)
        match = np.argmin(dist, axis=1)
        guard_one += bool(np.all(dist[np.arange(n), match] <= _TOL_MATCH))
        guard_two += bool(np.unique(match).size == n)
        guard_three += bool(np.all(np.abs(weights[match] - weights) <= _TOL_MATCH))
    assert (guard_one, guard_two, guard_three) == (48, 0, 48), (
        "the duplicate witness no longer ISOLATES the bijectivity guard: "
        f"guard1={guard_one}/48 guard2={guard_two}/48 guard3={guard_three}/48"
    )

    # The type makes the illegal state unrepresentable, not merely guarded.
    with pytest.raises(ValueError, match="not a permutation"):
        Permutation(np.array([0, 0, 2]))


def test_o5_the_match_window_is_an_explicit_parameter_and_it_bites() -> None:
    r"""O5 — ⭐ the window is a first-class correctness parameter.

    The ERR-073 bijectivity guard does NOT make a nearest-neighbour certifier
    sound. MEASURED: a two-point set off-symmetry by ``1e-9`` certifies under
    a ``1e-7`` window with a **perfectly injective** ``π``. Window and
    bijection are INDEPENDENT failure modes.

    ``atol`` is keyword-only and REQUIRED on the production signature, so this
    gate is not signature-tautological: the varying input is genuinely
    reachable from production. The first assertion pins that.
    """
    import inspect

    signature = inspect.signature(RigidMotion.permutes)
    atol_param = signature.parameters["atol"]
    assert atol_param.kind is inspect.Parameter.KEYWORD_ONLY
    assert atol_param.default is inspect.Parameter.empty, (
        "atol must have no default — a module constant would make the "
        "'the window bites' claim signature-tautological"
    )

    sigma = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    for delta in (1e-15, 1e-12, 1e-9, 1e-6):
        perturbed = np.array([[1.0, 0.0, 0.0], [-1.0 + delta, 0.0, 0.0]])
        for window in (1e-13, 1e-7):
            pi = sigma.permutes(perturbed, atol=window)
            if delta <= window:
                assert pi is not None, (delta, window)
                # …and when it certifies, the map really is a bijection: the
                # window, not injectivity, is what admitted it.
                np.testing.assert_array_equal(pi.indices, [1, 0])
            else:
                assert pi is None, (delta, window)


def test_o6_preserves_is_permutes_plus_exactly_one_weight_guard() -> None:
    r"""O6 — the weight guard, positive AND negative.

    On equal weights ``preserves`` and ``permutes`` answer identically. The
    NEGATIVE leg needs UNEQUAL weights arranged non-invariantly, which is
    impossible on any shipped quadrature (weights are equal within an orbit) —
    so the witness is hand-built: two antipodal points with weights ``(1, 2)``
    under ``σ_x``. MEASURED: ``π = [1, 0]`` is a *valid* permutation while
    ``w[π] = (2, 1) ≠ w``, so ``permutes`` SUCCEEDS and ``preserves`` FAILS.
    That is guard 3 in isolation.
    """
    points = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    sigma = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])

    pi = sigma.permutes(points, atol=_TOL_MATCH)
    assert pi is not None
    np.testing.assert_array_equal(pi.indices, [1, 0])
    assert (
        sigma.preserves(
            points, np.array([1.0, 2.0]), atol=_TOL_MATCH, weight_atol=_TOL_MATCH
        )
        is None
    ), "unequal weights on a swapped pair must NOT be preserved"
    assert (
        sigma.preserves(
            points, np.array([1.0, 1.0]), atol=_TOL_MATCH, weight_atol=_TOL_MATCH
        )
        == pi
    )
    # Equal weights: preserves ≡ permutes on every shipped-style fixture.
    for label, pts, elements in _point_set_cases():
        weights = np.full(len(pts), 0.31)
        for g in elements:
            assert g.preserves(
                pts, weights, atol=_TOL_MATCH, weight_atol=_TOL_MATCH
            ) == g.permutes(pts, atol=_TOL_MATCH), label
    with pytest.raises(ValueError, match="weights must have shape"):
        sigma.preserves(
            points, np.ones(3), atol=_TOL_MATCH, weight_atol=_TOL_MATCH
        )


def test_o7_orbit_stabiliser_against_the_published_wyckoff_table() -> None:
    r"""O7 — ⭐ ``|orbit|·|Stab| = |G|`` AND the ABSOLUTE site symmetries (P5).

    The relation alone is only a self-consistency check: it is invariant under
    conjugating the whole realization by an arbitrary ``R``, so it survives a
    uniformly-rotated group. Only the ABSOLUTE orbit types catch that class —
    hence the published Wyckoff decomposition of the cube (site symmetries
    ``C_3v`` / ``C_4v`` / ``C_2v``) and of the square.

    The 1-D rows are the cheapest non-trivial stabiliser in the file: an
    ODD-``n`` centred lattice has exactly one fixed point (the centre cell),
    an even-``n`` one has none.
    """
    for d, table in ((3, _WYCKOFF_CUBE), (2, _WYCKOFF_SQUARE)):
        group = close_group(_hyperoctahedral(d))
        order = len(group)
        assert order == 2**d * math.factorial(d)
        for label, (build, expected_orbit, expected_stabiliser) in table.items():
            points = np.asarray(build(), dtype=float)
            stabiliser = np.zeros(len(points), dtype=int)
            for g in group:
                pi = g.permutes(points, atol=_TOL_MATCH)
                assert pi is not None, f"{label}: not invariant under the full group"
                stabiliser[pi.fixed_points] += 1
            assert len(points) == expected_orbit, label
            assert set(stabiliser.tolist()) == {expected_stabiliser}, (
                f"{label}: |Stab| = {sorted(set(stabiliser.tolist()))}, "
                f"expected {{{expected_stabiliser}}}"
            )
            assert expected_orbit * expected_stabiliser == order, label

    mirror_group = [RigidMotion.identity(1), RigidMotion.reflection(normal=[1.0])]
    for n, expected_fixed in ((6, 0), (7, 1)):
        lattice = _cell_centres(n, origin=-0.05 * n)
        stabiliser = np.zeros(n, dtype=int)
        for g in mirror_group:
            pi = g.permutes(lattice, atol=_TOL_MATCH)
            assert pi is not None
            stabiliser[pi.fixed_points] += 1
        assert int((stabiliser > 1).sum()) == expected_fixed, (n, stabiliser.tolist())
        for i in range(n):
            orbit_length = len(mirror_group) // int(stabiliser[i])
            assert orbit_length * int(stabiliser[i]) == len(mirror_group)


def test_o8_the_seat_is_forced_by_the_centroid_theorem() -> None:
    r"""O8 — ⭐⭐ the centroid of a ``G``-preserved weighted point set is
    ``G``-FIXED, so the seat is a THEOREM rather than a convention.

    Proof in three lines: ``g(c) = Qc + t = Σwᵢ g(xᵢ)/Σw = Σwᵢ x_{π(i)}/Σw``,
    and re-indexing with ``w_{π(i)} = wᵢ`` returns ``c``. That is ruling R10's
    pure-math backing — and the operational consequence is the campaign's
    motivating defect: a mesh starting at ``origin = 0.0`` fails a mirror test
    only because the mirror is seated at the origin instead of the centroid.

    MEASURED on a cube shifted to ``(2.5, −1.25, 0.75)``: all **48/48**
    *seated* elements preserve it (``max‖g(c) − c‖ = 4.4e-16``) while only
    **1/48** *unseated* ones do — the identity.
    """
    shift = np.array([2.5, -1.25, 0.75])
    points = _cube_vertices() + shift
    weights = np.full(len(points), 0.37)
    centroid = (weights[:, None] * points).sum(axis=0) / weights.sum()
    np.testing.assert_allclose(centroid, shift, atol=_TOL_GENERAL, rtol=0.0)

    seated_hits = unseated_hits = 0
    worst_centroid_drift = 0.0
    for g in _hyperoctahedral(3):
        seated = g.seated_at(centroid)
        if (
            seated.preserves(
                points, weights, atol=_TOL_MATCH, weight_atol=_TOL_MATCH
            )
            is not None
        ):
            seated_hits += 1
        if (
            g.preserves(points, weights, atol=_TOL_MATCH, weight_atol=_TOL_MATCH)
            is not None
        ):
            unseated_hits += 1
        worst_centroid_drift = max(
            worst_centroid_drift,
            float(np.max(np.abs(seated.on_points(centroid) - centroid))),
        )
    assert seated_hits == 48, f"only {seated_hits}/48 seated elements preserve the set"
    assert unseated_hits == 1, (
        f"{unseated_hits}/48 UNSEATED elements preserved the set — expected only "
        "the identity; the affine part is not being exercised"
    )
    assert worst_centroid_drift <= _TOL_GENERAL, worst_centroid_drift

    # The 1-D analogue: the production-mesh case.
    for n in (6, 7):
        origin, h = 0.0, 0.1
        lattice = _cell_centres(n, origin=origin, h=h)
        cell_volumes = np.full(n, h)
        centre = float(
            (cell_volumes[:, None] * lattice).sum() / cell_volumes.sum()
        )
        assert centre == pytest.approx(origin + n * h / 2.0, abs=1e-15)
        mirror = RigidMotion.reflection(normal=[1.0])
        assert mirror.permutes(lattice, atol=_TOL_MATCH) is None
        seated_pi = mirror.seated_at([centre]).permutes(lattice, atol=_TOL_MATCH)
        assert seated_pi is not None
        np.testing.assert_array_equal(seated_pi.indices, np.arange(n)[::-1])


# ===========================================================================
# X — gates the landed design earns
# ===========================================================================


@pytest.mark.parametrize("d", [1, 2, 3])
def test_x1_the_two_actions_are_distinct_and_consistent(d: int) -> None:
    r"""X1 — ⭐ ``on_points`` and ``on_directions`` are separate for a reason.

    A direction has no position, so ``t`` must not act on it; applying the
    affine map to a direction silently DENORMALISES unit vectors. Three legs:
    the difference is exactly ``t``; the direction action is norm-preserving;
    and — the sharp one — ``on_directions`` is INVARIANT under ``seated_at``
    for every centre. Where a mirror sits does not change which way it points,
    and that is what makes the boundary law
    ``ψ_in(x,Ω) = ψ_out(g⁻¹x, Q_g⁻¹Ω)`` spellable.
    """
    rng = _rng(22)
    vectors = rng.normal(size=(16, d))
    for _ in range(20):
        g = _random_motion(rng, d)
        # Stated in the direction that IS a float theorem: `on_points` is
        # `x @ Qᵀ + t` and `on_directions` is `x @ Qᵀ`, so adding `t` back
        # recomputes the same expression. MEASURED bit-exact 6000/6000.
        # (The subtraction `on_points − on_directions == t` is NOT exact —
        # `fl(a + t) − a ≠ t` in general — and asserting it would be a false
        # red at 2.2e-16.)
        np.testing.assert_array_equal(
            g.on_points(vectors), g.on_directions(vectors) + g.translation
        )
        _assert_close(
            np.linalg.norm(g.on_directions(vectors), axis=1),
            np.linalg.norm(vectors, axis=1),
        )
        for _ in range(3):
            centre = rng.normal(size=d)
            np.testing.assert_array_equal(
                g.seated_at(centre).on_directions(vectors),
                g.on_directions(vectors),
                err_msg="seating changed the direction action",
            )


def test_x2_the_circle_point_constructor_preserves_exactness() -> None:
    r"""X2 — the exact-turn path EXISTS, and the angle path cannot reproduce it.

    ``np.cos(np.pi/2) = 6.1e-17``, not ``0``. A caller holding an exact circle
    point — from a root-of-unity construction — keeps that exactness through
    ``rotation_from_circle_point``; ``rotation(angle=π/2)`` does not.

    This is the API-side gate for campaign step G4, deliberately written
    WITHOUT importing ``roots_of_unity`` (a forbidden reference here): G2 pins
    that the exact path exists and agrees with the angle path to ``1e-14``;
    G4 supplies the exact circle points.
    """
    plane = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    exact = RigidMotion.rotation_from_circle_point(plane=plane, point=(0.0, 1.0))
    via_angle = RigidMotion.rotation(plane=plane, angle=np.pi / 2.0)

    np.testing.assert_array_equal(
        exact.linear,
        np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    assert np.count_nonzero(exact.linear) == 3, "the exact quarter turn has exact zeros"
    assert np.count_nonzero(via_angle.linear) > 3, (
        "np.cos(pi/2) is no longer inexact — this gate's premise has changed"
    )
    np.testing.assert_allclose(via_angle.linear, exact.linear, atol=1e-14, rtol=0.0)

    with pytest.raises(ValueError, match="unit circle"):
        RigidMotion.rotation_from_circle_point(plane=plane, point=(0.5, 0.5))
    with pytest.raises(ValueError, match="ORTHONORMAL"):
        RigidMotion.rotation_from_circle_point(
            plane=[[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], point=(0.0, 1.0)
        )


def test_x3_a_glide_has_infinite_order_though_its_linear_part_does_not() -> None:
    r"""X3 — ⭐ an element translating off its own fixed subspace generates an
    INFINITE cyclic group.

    ``σ_x`` has order 2; compose it with a translation ALONG the mirror and
    the result is a glide reflection, whose square is a pure translation — so
    the element has infinite order and ``element_order()`` must report
    ``None``. Reporting the linear part's order would be a lie about the
    element, and it is precisely the error a "just take the order of ``Q``"
    implementation makes.

    The element-level analogue of A6, and the control leg pins the asymmetry:
    translating ALONG the normal instead keeps the order finite (it is a
    seated reflection).
    """
    mirror = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    assert mirror.element_order() == 2

    glide = RigidMotion.translation_by([0.0, 1.0, 0.0]) @ mirror
    assert glide.element_order() is None, "a glide reflection has infinite order"
    squared = glide @ glide
    np.testing.assert_allclose(squared.linear, np.eye(3), atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(
        squared.translation, [0.0, 2.0, 0.0], atol=1e-15, rtol=0.0
    )

    # Control: a translation ALONG the normal is a seated reflection — finite.
    seated = RigidMotion.translation_by([1.0, 0.0, 0.0]) @ mirror
    assert seated.element_order() == 2

    # And a screw: a rotation carrying a translation along its own axis.
    screw = RigidMotion.translation_by([0.0, 0.0, 1.0]) @ RigidMotion.rotation_about_axis(
        axis=[0.0, 0.0, 1.0], angle=np.pi / 2.0
    )
    assert screw.element_order() is None


# ===========================================================================
# H — dimension genericity (why R^1 is first-class)
# ===========================================================================


def test_h1_r1_is_complete_and_correct() -> None:
    r"""H1 — ⭐ ``O(1) = {+1, −1}`` and ``SO(1) = {e}``.

    ``R^1`` is the case that dissolves the 1-D/3-D arm split, so it is
    first-class rather than an afterthought: ``reflection(normal=[1])`` is
    ``[[−1]]`` with ``det = −1``, ``dim Fix = 0`` and order 2; seated at ``c``
    it is ``x ↦ 2c − x``, fixing ``c`` **exactly**; and a 1-D rotation raises
    naming ``SO(1) = {e}`` rather than blaming the caller's plane.
    """
    sigma = RigidMotion.reflection(normal=[1.0])
    np.testing.assert_array_equal(sigma.linear, [[-1.0]])
    assert sigma.determinant == pytest.approx(-1.0, abs=0.0)
    assert sigma.fixed_subspace_dimension == 0
    assert sigma.element_order() == 2
    assert not sigma.is_proper

    centre = 0.3
    seated = sigma.seated_at([centre])
    np.testing.assert_array_equal(seated.translation, [2.0 * centre])
    assert bool(seated.fixes(np.array([[centre]]), atol=0.0)[0]), (
        "x -> 2c - x must fix c EXACTLY in R^1"
    )
    np.testing.assert_array_equal(
        seated.on_points(np.array([[1.0], [0.0]])), [[2.0 * centre - 1.0], [2.0 * centre]]
    )

    # O(1) has exactly two elements; SO(1) exactly one.
    assert len(close_group([sigma])) == 2
    assert len(close_group([], dimension=1)) == 1
    np.testing.assert_array_equal(RigidMotion.inversion(1).linear, sigma.linear)

    with pytest.raises(ValueError, match=r"SO\(1\) = \{e\}"):
        RigidMotion.rotation(plane=[[1.0], [1.0]], angle=0.5)


@pytest.mark.parametrize("source,target", [(1, 2), (1, 3), (2, 3), (2, 5), (3, 4)])
def test_h2_the_embedding_is_a_homomorphism(source: int, target: int) -> None:
    r"""H2 — ⭐ ``ι: E(d) → E(D)``, ``(Q,t) ↦ (diag(Q,I), (t,0))``.

    The operation that makes "dimension generic" a property OF THE TYPE rather
    than a habit of its callers. Without it, relating a 1-D element to its 3-D
    counterpart is ``np.eye(3)`` padding at each call site — and then the
    1-D/3-D split this module exists to delete has merely MOVED.

    Four laws: it is a group homomorphism; it preserves the determinant; it
    raises ``dim Fix`` by exactly ``D − d`` (so a reflection stays a reflection
    and keeps fixing a hyperplane); and it acts identically on the first ``d``
    coordinates while fixing the complement pointwise. Lowering dimension
    raises — there is no canonical restriction.
    """
    rng = _rng(23 + target)
    pool: list[RigidMotion] = [g for _, g in _reflection_family(source)]
    pool.append(RigidMotion.inversion(source))
    pool.append(RigidMotion.translation_by(rng.normal(size=source)))
    if source >= 2:
        pool.append(_random_motion(rng, source))

    for g in pool:
        lifted = g.embedded_in(target)
        assert lifted.dimension == target
        assert lifted.determinant == pytest.approx(g.determinant, abs=1e-12)
        assert (
            lifted.fixed_subspace_dimension
            == g.fixed_subspace_dimension + (target - source)
        )
        # Acts identically on the first `source` coordinates; fixes the rest.
        points = rng.normal(size=(12, target))
        moved = lifted.on_points(points)
        _assert_close(moved[:, :source], g.on_points(points[:, :source]))
        np.testing.assert_array_equal(moved[:, source:], points[:, source:])

    for a, b in itertools.product(pool, repeat=2):
        composed = (a @ b).embedded_in(target)
        lifted_product = a.embedded_in(target) @ b.embedded_in(target)
        np.testing.assert_allclose(
            composed.linear, lifted_product.linear, atol=1e-15, rtol=0.0
        )
        np.testing.assert_allclose(
            composed.translation, lifted_product.translation, atol=1e-15, rtol=0.0
        )

    with pytest.raises(ValueError, match="cannot lower dimension"):
        RigidMotion.identity(target).embedded_in(source)


def test_h2_the_1d_mirror_embeds_as_the_coordinate_mirror() -> None:
    r"""H2, the load-bearing instance — ``ι(σ_{[1]}) = σ_x`` at d = 3.

    The same element, read in a bigger room: in ``R^1`` it fixes a
    0-dimensional set (the origin); embedded in ``R^3`` it fixes the plane
    ``x = 0``. The complementary axes are fixed pointwise — which is exactly
    the content of the ``(μ, 0, 0)`` embedding the angular layer performs on
    its data, and the reason the 1-D/3-D arm split is deletable rather than
    merely reconcilable.

    Bit-exact, so nothing can drift here.
    """
    lifted = RigidMotion.reflection(normal=[1.0]).embedded_in(3)
    np.testing.assert_array_equal(lifted.linear, np.diag([-1.0, 1.0, 1.0]))
    np.testing.assert_array_equal(
        lifted.linear, RigidMotion.reflection(normal=[1.0, 0.0, 0.0]).linear
    )
    assert lifted.fixed_subspace_dimension == 2
    assert lifted.determinant == pytest.approx(-1.0, abs=0.0)
    # Embedding twice equals embedding once.
    twice = RigidMotion.reflection(normal=[1.0]).embedded_in(2).embedded_in(3)
    np.testing.assert_array_equal(twice.linear, lifted.linear)


def test_h3_mixing_dimensions_raises() -> None:
    """H3 — a d = 1 element composed with a d = 3 element is unspellable.

    Numpy would happily broadcast some of these. The raising calls are the
    PRODUCTION ``__matmul__`` and ``close_group``, never a locally-constructed
    exception.
    """
    one = RigidMotion.reflection(normal=[1.0])
    three = RigidMotion.reflection(normal=[1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="different dimensions"):
        _ = one @ three
    with pytest.raises(ValueError, match="different dimensions"):
        _ = three @ one
    with pytest.raises(ValueError, match="same dimension"):
        close_group([one, three])
    with pytest.raises(ValueError, match="contradicts"):
        close_group([three], dimension=2)
    with pytest.raises(ValueError, match="trailing dimension"):
        three.on_points(np.zeros((4, 2)))
    with pytest.raises(ValueError, match="trailing dimension"):
        three.on_directions(np.zeros((4, 2)))
    # `permutes` owns no width check of its own since 2026-09-02 (#429
    # tracker 2.2b): the match is `permutation_between` on `on_points`'s
    # image, so the refusal is `on_points`'s — one home for the contract.
    with pytest.raises(ValueError, match="trailing dimension"):
        three.permutes(np.zeros((4, 2)), atol=_TOL_MATCH)
    with pytest.raises(ValueError, match="compose permutations"):
        _ = Permutation.identity(3) @ Permutation.identity(4)
