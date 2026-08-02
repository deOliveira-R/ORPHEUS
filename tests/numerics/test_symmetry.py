r"""Foundation tests for :mod:`orpheus.numerics.symmetry`.

Every test here is tagged ``foundation`` because the assertions verify
software invariants of the static :math:`O(3)`-subgroup lattice and the
node-permutation orbit-closure check, not a physics-equation claim with
an L0..L3 verification ladder. The mathematical content (the lattice
itself) is fixed by group theory (Hamermesh 1962 §2.5, §9.4) and the
quadrature-construction theorems (Lebedev 1976; Carlson & Lathrop 1968)
that *guarantee* the invariance the tests confirm.

The tests fall into two groups:

1. **Containment lattice** — exhaustive coverage of every named
   relation plus the parameterised-family rules
   (:math:`C_n \subseteq C_m` iff :math:`n \mid m`, etc.). Reflexivity
   is asserted for every named entry.

2. **Invariance against existing quadratures** — a
   :class:`~orpheus.numerics.measure.DiscreteMeasure` is built from
   each quadrature's ``mu_x``/``mu_y``/``mu_z``/``weights`` arrays
   (the "bridge" pattern that the Issue 4 adapter will formalise) and
   the appropriate group's :meth:`is_invariant` is asserted true.
   Includes a negative case: Lebedev grids are :math:`O_h`-invariant
   *not* :math:`I_h`-invariant.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.numerics.quadrature import Quadrature


# ============================================================================
# Helpers
# ============================================================================


def _measure_from_sphere_quad(q) -> DiscreteMeasure:
    """Build an :math:`S^2` DiscreteMeasure from any quadrature exposing
    ``mu_x``/``mu_y``/``mu_z``/``weights``.

    Mirrors the bridge pattern that the Issue 4 quadrature adapter will
    formalise; tests can use it today without depending on Issue 4.
    """
    nodes = np.column_stack([q.mu_x, q.mu_y, q.mu_z])
    return DiscreteMeasure(nodes=nodes, weights=q.weights, support="S^2")


def _measure_from_1d_quad(q) -> DiscreteMeasure:
    """Build a 1-D DiscreteMeasure on :math:`[-1, 1]` from a 1-D quadrature."""
    return DiscreteMeasure(
        nodes=q.mu_x, weights=q.weights, support="[-1,1]",
    )


# ============================================================================
# Containment lattice — named entries
# ============================================================================


_NAMED = [
    SubgroupOfO3.Trivial,
    SubgroupOfO3.Z2,
    SubgroupOfO3.SO2,
    SubgroupOfO3.O2,
    SubgroupOfO3.OctahedralOh,
    SubgroupOfO3.IcosahedralIh,
    SubgroupOfO3.SO3,
    SubgroupOfO3.O3,
]


@pytest.mark.foundation
@pytest.mark.parametrize("G", _NAMED, ids=lambda g: g.name)
def test_named_reflexivity(G: SubgroupOfO3) -> None:
    """Every named subgroup contains itself."""
    assert G.contains(G)
    assert G.is_subgroup_of(G)


@pytest.mark.foundation
@pytest.mark.parametrize("G", _NAMED, ids=lambda g: g.name)
def test_trivial_inside_every_named_group(G: SubgroupOfO3) -> None:
    """The trivial group sits at the bottom of the lattice."""
    assert G.contains(SubgroupOfO3.Trivial)


@pytest.mark.foundation
@pytest.mark.parametrize(
    "G",
    [
        SubgroupOfO3.Z2,
        SubgroupOfO3.SO2,
        SubgroupOfO3.O2,
        SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh,
        SubgroupOfO3.SO3,
    ],
    ids=lambda g: g.name,
)
def test_o3_contains_every_proper_subgroup(G: SubgroupOfO3) -> None:
    """The orthogonal group sits at the top of the lattice."""
    assert SubgroupOfO3.O3.contains(G)


@pytest.mark.foundation
def test_so2_chain() -> None:
    """``SO(2) ⊂ O(2) ⊂ O(3)`` is the canonical axisymmetric tower."""
    assert SubgroupOfO3.O2.contains(SubgroupOfO3.SO2)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.O2)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.SO2)  # transitive


@pytest.mark.foundation
def test_z2_chain_via_octahedral() -> None:
    """``Z_2 ⊂ O_h ⊂ O(3)`` — the Cartesian-symmetry tower."""
    assert SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Z2)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.OctahedralOh)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Z2)


@pytest.mark.foundation
def test_octahedral_not_in_so3() -> None:
    """``O_h`` contains improper rotations (inversion, reflections), so
    it is **not** a subgroup of the proper-rotation group ``SO(3)``.

    This matches the chemistry/Lebedev convention used for the
    project's quadratures (Lebedev grids carry the full :math:`O_h`,
    not just the chiral subgroup :math:`O`).
    """
    assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.OctahedralOh)


@pytest.mark.foundation
def test_icosahedral_not_in_so3() -> None:
    """Same reasoning as above — :math:`I_h` includes inversion."""
    assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.IcosahedralIh)


@pytest.mark.foundation
def test_so3_in_o3_only() -> None:
    """``SO(3)`` is contained only in ``O(3)`` (and itself)."""
    for G in _NAMED:
        if G in (SubgroupOfO3.SO3, SubgroupOfO3.O3):
            assert G.contains(SubgroupOfO3.SO3), G
        else:
            assert not G.contains(SubgroupOfO3.SO3), G


# ============================================================================
# Containment — parameterised families (Cn, Dnh)
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 3, 4, 6, 8])
def test_cn_reflexivity(n: int) -> None:
    """``C_n ⊃ C_n`` for every ``n``."""
    G = SubgroupOfO3.Cn(n)
    assert G.contains(G)


@pytest.mark.foundation
@pytest.mark.parametrize(
    "outer_n, inner_n, expected",
    [
        (6, 2, True),     # C_2 ⊂ C_6
        (6, 3, True),     # C_3 ⊂ C_6
        (6, 6, True),     # reflexivity
        (4, 2, True),     # C_2 ⊂ C_4
        (6, 4, False),    # 4 ∤ 6
        (3, 2, False),    # 2 ∤ 3
        (5, 2, False),    # 2 ∤ 5
    ],
)
def test_cyclic_containment_by_divisibility(
    outer_n: int, inner_n: int, expected: bool,
) -> None:
    """``C_n ⊆ C_m`` iff ``n ∣ m``."""
    assert (
        SubgroupOfO3.Cn(outer_n).contains(SubgroupOfO3.Cn(inner_n))
    ) is expected


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 3, 4, 6, 8])
def test_cn_in_so2(n: int) -> None:
    r""":math:`SO(2) = \bigcup_n C_n`, so every :math:`C_n` is inside."""
    assert SubgroupOfO3.SO2.contains(SubgroupOfO3.Cn(n))


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 3, 4, 6])
def test_cn_in_dnh_principal_axis(n: int) -> None:
    """``C_n ⊆ D_{nh}`` — principal-axis rotations sit inside the
    dihedral group."""
    assert SubgroupOfO3.Dnh(n).contains(SubgroupOfO3.Cn(n))


@pytest.mark.foundation
def test_dnh_reflection_in_dnh() -> None:
    """``Z_2`` (a single reflection) sits inside every ``D_{nh}``."""
    for n in (1, 2, 3, 4, 6):
        assert SubgroupOfO3.Dnh(n).contains(SubgroupOfO3.Z2)


@pytest.mark.foundation
def test_dnh_in_o3() -> None:
    """``D_{nh} ⊂ O(3)`` for every ``n`` (and ``D_{nh} ⊂ O(2)``)."""
    for n in (1, 2, 3, 4, 6):
        assert SubgroupOfO3.O3.contains(SubgroupOfO3.Dnh(n))
        assert SubgroupOfO3.O2.contains(SubgroupOfO3.Dnh(n))
        # Not in SO(2)/SO(3) — D_nh contains improper rotations.
        assert not SubgroupOfO3.SO2.contains(SubgroupOfO3.Dnh(n))
        assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.Dnh(n))


# ============================================================================
# is_subgroup_of — readability synonym
# ============================================================================


@pytest.mark.foundation
def test_is_subgroup_of_reverses_contains() -> None:
    """``A.is_subgroup_of(B)`` is equivalent to ``B.contains(A)``."""
    assert SubgroupOfO3.Z2.is_subgroup_of(SubgroupOfO3.OctahedralOh)
    assert SubgroupOfO3.SO2.is_subgroup_of(SubgroupOfO3.O3)
    assert SubgroupOfO3.OctahedralOh.is_subgroup_of(SubgroupOfO3.O3)
    assert not SubgroupOfO3.O3.is_subgroup_of(SubgroupOfO3.OctahedralOh)


# ============================================================================
# Invariance — Lebedev sphere quadrature (Oh by construction)
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("order", [11, 17, 23])
def test_lebedev_is_octahedral_invariant(order: int) -> None:
    r"""Lebedev (1976) constructs sphere quadratures by enforcing
    :math:`O_h` symmetry through the choice of generating points
    (face-centres, edge-centres, vertices, and free-parameter orbits).
    The construction *guarantees* :math:`O_h`-invariance — this test
    confirms the project's wrapper preserves it.
    """
    q = Quadrature.lebedev(order=order)
    mu = _measure_from_sphere_quad(q)
    assert SubgroupOfO3.OctahedralOh.is_invariant(mu)


@pytest.mark.foundation
def test_lebedev_is_NOT_icosahedral_invariant() -> None:
    r"""Lebedev grids are :math:`O_h`-symmetric, **not**
    :math:`I_h`-symmetric (the icosahedron and the cube have
    incompatible vertex symmetries — :math:`O_h \cap I_h` is just the
    common subgroup, not either one).

    A generic Lebedev grid fails the icosahedral fingerprint. The test
    is therefore a **negative** case: it would erroneously pass if the
    invariance check were too permissive (e.g. if it only checked
    weight sums rather than orbit closure).
    """
    q = Quadrature.lebedev(order=17)
    mu = _measure_from_sphere_quad(q)
    assert not SubgroupOfO3.IcosahedralIh.is_invariant(mu)


# ============================================================================
# Invariance — level-symmetric SN (Oh by construction)
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("sn_order", [2, 4, 6, 8])
def test_level_symmetric_is_octahedral_invariant(sn_order: int) -> None:
    r"""Level-symmetric :math:`S_N` quadrature (Carlson & Lathrop 1968)
    is built by picking direction cosines on each octant and reflecting
    through every coordinate-plane combination — the construction
    enforces :math:`O_h`-invariance directly.
    """
    q = Quadrature.level_symmetric(sn_order=sn_order)
    mu = _measure_from_sphere_quad(q)
    assert SubgroupOfO3.OctahedralOh.is_invariant(mu)


# ============================================================================
# Invariance — 1-D Gauss-Legendre (SO2-invariant by degeneracy)
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("n", [4, 8, 16])
def test_gauss_legendre_1d_so2_invariant(n: int) -> None:
    """A 1-D measure on :math:`[-1, 1]` is trivially
    :math:`SO(2)`-invariant: there is no azimuthal coordinate to rotate.

    This is the "axisymmetric slab" tag — slab geometry uses GL1D
    because the slab's azimuthal symmetry is automatically inherited.
    """
    q = Quadrature.gauss_legendre(n_ordinates=n)
    mu = _measure_from_1d_quad(q)
    assert SubgroupOfO3.SO2.is_invariant(mu)


@pytest.mark.foundation
@pytest.mark.parametrize("n", [4, 8, 16])
def test_gauss_legendre_1d_z2_reflective(n: int) -> None:
    r"""Gauss-Legendre nodes are symmetric about :math:`x = 0` (roots
    of an even polynomial come in :math:`\pm` pairs with equal
    Christoffel weights), so the measure is :math:`Z_2`-invariant
    under :math:`x \to -x`.
    """
    q = Quadrature.gauss_legendre(n_ordinates=n)
    mu = _measure_from_1d_quad(q)
    assert SubgroupOfO3.Z2.is_invariant(mu)


# ============================================================================
# Invariance — Trivial group is universal
# ============================================================================


@pytest.mark.foundation
def test_trivial_invariance_is_universal() -> None:
    """Every measure is invariant under the identity group.

    Sanity check on the dispatch logic: the trivial group's
    ``is_invariant`` short-circuit must return ``True`` regardless
    of the measure's structure.
    """
    # An asymmetric, non-symmetric, non-axisymmetric measure.
    nodes = np.array([[0.5, 0.3, 0.1], [0.7, -0.2, 0.4]])
    weights = np.array([1.0, 2.0])
    mu = DiscreteMeasure(nodes=nodes, weights=weights, support="S^2")
    assert SubgroupOfO3.Trivial.is_invariant(mu)


@pytest.mark.foundation
def test_asymmetric_measure_not_octahedral_invariant() -> None:
    """An asymmetric measure on :math:`S^2` fails the :math:`O_h`
    invariance check.

    Constructed by hand: two arbitrary points on the sphere with
    different weights. A correctly-implemented orbit-closure check
    rejects this.
    """
    p1 = np.array([1.0, 0.0, 0.0]) / np.sqrt(1.0)
    p2 = np.array([0.7, 0.5, 0.5]) / np.linalg.norm([0.7, 0.5, 0.5])
    nodes = np.vstack([p1, p2])
    weights = np.array([1.0, 2.0])
    mu = DiscreteMeasure(nodes=nodes, weights=weights, support="S^2")
    assert not SubgroupOfO3.OctahedralOh.is_invariant(mu)


# ============================================================================
# Repr / equality smoke
# ============================================================================


@pytest.mark.foundation
def test_repr_and_equality() -> None:
    """Singleton named entries compare equal to themselves; the
    parameterised families compare by their ``n``."""
    assert SubgroupOfO3.Z2 == SubgroupOfO3.Z2
    assert SubgroupOfO3.Cn(6) == SubgroupOfO3.Cn(6)
    assert SubgroupOfO3.Cn(4) != SubgroupOfO3.Cn(6)
    assert SubgroupOfO3.Dnh(4) != SubgroupOfO3.Dnh(6)
    # Repr should not raise and should mention the group name.
    assert "Z2" in repr(SubgroupOfO3.Z2)
    assert "Cn(6)" in repr(SubgroupOfO3.Cn(6))
    assert "Dnh(4)" in repr(SubgroupOfO3.Dnh(4))


# ============================================================================
# ERR-072 — a CONTINUOUS group cannot be decided by a FINITE sample
# ============================================================================
#
# `is_invariant` asks node-set closure: "does every g in G permute the nodes
# among themselves, weights matched?"  For a continuous G that question has an
# exact, decidable answer that no sampling can approximate -- a finite sample
# generates a finite SUBgroup, and closure under a subgroup is strictly weaker
# than closure under G.  The old code sampled {0,90,180,270} degrees about z,
# i.e. C_4, and called it SO(2); the certification it produced was therefore a
# function of `n_phi mod 4`.  These gates pin the exact criterion instead.


def _production_sphere_rules() -> list[tuple[str, object]]:
    """Every S^2 rule the tree actually ships, across its parameter range."""
    rules: list[tuple[str, object]] = []
    for n_mu, n_phi in ((2, 4), (4, 4), (4, 8), (4, 12), (4, 16), (6, 12)):
        rules.append((f"product({n_mu},{n_phi})",
                      Quadrature.product(n_mu=n_mu, n_phi=n_phi)))
    for sn_order in (4, 8, 16):
        rules.append((f"level_symmetric({sn_order})",
                      Quadrature.level_symmetric(sn_order=sn_order)))
    for order in (5, 11, 17):
        rules.append((f"lebedev({order})", Quadrature.lebedev(order=order)))
    return rules


@pytest.mark.foundation
@pytest.mark.catches("ERR-072")
def test_no_discrete_cubature_is_so2_invariant() -> None:
    r"""No finite angular cubature is :math:`SO(2)`-invariant.

    The :math:`SO(2)` orbit of any node off the :math:`z`-axis is a whole
    circle, so a finite set containing that node cannot be closed. Every
    shipped rule places nodes off-axis, hence every one must answer ``False``.
    """
    for name, q in _production_sphere_rules():
        mu = _measure_from_sphere_quad(q)
        assert not SubgroupOfO3.SO2.is_invariant(mu), (
            f"{name} certified SO(2)-invariant, but a finite point set with "
            f"off-axis nodes has infinite SO(2) orbits and cannot be closed"
        )


@pytest.mark.foundation
@pytest.mark.catches("ERR-072")
def test_so2_verdict_is_not_a_function_of_n_phi_mod_4() -> None:
    r"""The ERR-072 signature: the old answer tracked ``n_phi mod 4``.

    Sampling :math:`\{0, 90, 180, 270\}^\circ` generates :math:`C_4`, so a
    product rule whose azimuthal count is a multiple of 4 closed under the
    sample and was certified. ``n_phi`` in 4/8/12/16 read ``True`` while
    2/3/5/6/7 read ``False`` -- a verdict that depends on the sample, not on
    the group. The correct verdict is constant across the family.
    """
    verdicts = {
        n_phi: SubgroupOfO3.SO2.is_invariant(
            _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        )
        for n_phi in (2, 3, 4, 5, 6, 7, 8, 12, 16)
    }
    assert set(verdicts.values()) == {False}, (
        f"SO(2) verdict varies across the product family: {verdicts}. A "
        f"correct criterion cannot depend on n_phi mod 4 -- that dependence "
        f"IS the sampled-subgroup bug"
    )


def _icosahedron_vertices() -> np.ndarray:
    r"""The 12 icosahedron vertices, built from the golden ratio.

    Constructed here rather than from :mod:`~orpheus.numerics.symmetry`'s own
    operator tables, so the fixture is structurally independent of the module
    under test (``vv-principles`` L11).
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = []
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            raw.extend([
                [0.0, s1 * 1.0, s2 * phi],
                [s1 * 1.0, s2 * phi, 0.0],
                [s2 * phi, 0.0, s1 * 1.0],
            ])
    v = np.array(raw)
    return v / np.linalg.norm(v, axis=1, keepdims=True)


@pytest.mark.foundation
@pytest.mark.catches("ERR-072")
def test_icosahedral_vertex_set_is_not_so3_invariant() -> None:
    r"""The discriminating fixture for the :math:`SO(3)` half of ERR-072.

    The old code tested :math:`I_h` closure and *called it* :math:`SO(3)`. On
    every rule the tree ships that happened to give the right answer for the
    wrong reason -- none of them is :math:`I_h`-invariant, so the check said
    ``False`` by accident. The 12 icosahedron vertices ARE :math:`I_h`-closed,
    which is exactly where the two criteria part company: the old check
    certifies them :math:`SO(3)`-invariant, and they are not. Every
    :math:`SO(3)` orbit of a non-origin point is a whole 2-sphere, so only a
    measure supported at the origin can be :math:`SO(3)`-closed.

    Also pins the second edge: since :math:`-I \in I_h`, the old ``SO3`` and
    ``O3`` branches ran the SAME operator set and were identically equal for
    every input -- and 60 of those 120 matrices are improper, hence not in
    :math:`SO(3)` at all.
    """
    verts = _icosahedron_vertices()
    ico = DiscreteMeasure(
        nodes=verts, weights=np.full(len(verts), 1.0), support="S^2",
    )
    # Sanity: the fixture really is the I_h orbit it claims to be.
    assert len(verts) == 12
    assert SubgroupOfO3.IcosahedralIh.is_invariant(ico), (
        "fixture is not I_h-invariant -- it cannot discriminate the criteria"
    )
    assert not SubgroupOfO3.SO3.is_invariant(ico), (
        "the 12 icosahedron vertices certified SO(3)-invariant; a finite set "
        "cannot be, since every SO(3) orbit of a non-origin point is S^2"
    )
    assert not SubgroupOfO3.O3.is_invariant(ico)


@pytest.mark.foundation
def test_no_sphere_cubature_is_so3_or_o3_invariant() -> None:
    r"""Regression floor for :math:`SO(3)` / :math:`O(3)` on the shipped rules.

    Deliberately NOT marked ``catches("ERR-072")``: measured, the pre-fix code
    also answered ``False`` for every one of these rules, so this gate would
    not have caught the defect. It pins the answer against future drift; the
    catcher is
    :func:`test_icosahedral_vertex_set_is_not_so3_invariant`.
    """
    for name, q in _production_sphere_rules():
        mu = _measure_from_sphere_quad(q)
        assert not SubgroupOfO3.SO3.is_invariant(mu), f"{name} certified SO(3)"
        assert not SubgroupOfO3.O3.is_invariant(mu), f"{name} certified O(3)"


@pytest.mark.foundation
def test_axis_supported_measure_is_so2_invariant_but_not_so3() -> None:
    r"""Positive control + discriminator, so the SO(2) gate is not vacuous.

    A two-pole measure IS :math:`SO(2)`-invariant (both nodes are fixed by
    every rotation about :math:`z`) and IS :math:`O(2)`-invariant
    (:math:`\sigma_h` swaps them at equal weight) -- but is NOT
    :math:`SO(3)`-invariant, since the orbit of a pole is the whole sphere.
    A criterion that simply returned ``False`` would fail this test.
    """
    poles = DiscreteMeasure(
        nodes=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
        weights=np.array([1.0, 1.0]),
        support="S^2",
    )
    assert SubgroupOfO3.SO2.is_invariant(poles)
    assert SubgroupOfO3.O2.is_invariant(poles)
    assert not SubgroupOfO3.SO3.is_invariant(poles)

    # Unequal weights break sigma_h, hence O(2), while SO(2) survives.
    lopsided = DiscreteMeasure(
        nodes=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
        weights=np.array([1.0, 2.0]),
        support="S^2",
    )
    assert SubgroupOfO3.SO2.is_invariant(lopsided)
    assert not SubgroupOfO3.O2.is_invariant(lopsided)


@pytest.mark.foundation
def test_product_rule_honest_group_is_dnh_of_n_phi() -> None:
    r"""The product rule's true azimuthal group is the FINITE
    :math:`D_{n_\varphi h}` -- parameter-dependent, never :math:`SO(2)`.

    Deliberately NOT marked ``catches("ERR-072")``: measured, this gate is
    green under the pre-fix code too (the :math:`C_n` / :math:`D_{nh}`
    branches were always sound), so it never caught the defect. It states
    positively what the rule's symmetry IS -- the fact the registry's
    constant ``invariance_group=SO2`` tag contradicts.

    Verified in both directions: the rule IS :math:`D_{n_\varphi h}`- and
    :math:`C_{n_\varphi}`-invariant, and is NOT invariant under the
    twice-as-fine :math:`C_{2 n_\varphi}`. The second half is what forbids
    over-claiming: an azimuthal group of order :math:`2 n_\varphi` would
    require a node at every half-spacing, and there is none.
    """
    for n_phi in (4, 8, 12, 16):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        assert SubgroupOfO3.Cn(n_phi).is_invariant(mu), (
            f"product(4,{n_phi}) is not C_{n_phi}-invariant"
        )
        assert SubgroupOfO3.Dnh(n_phi).is_invariant(mu), (
            f"product(4,{n_phi}) is not D_{n_phi}h-invariant"
        )
        assert not SubgroupOfO3.Cn(2 * n_phi).is_invariant(mu), (
            f"product(4,{n_phi}) certified C_{2 * n_phi}-invariant -- that "
            f"would need a node at every half-spacing"
        )


# ============================================================================
# ERR-073 — the orbit match must be a BIJECTION, not merely a partner map
# ============================================================================


@pytest.mark.foundation
@pytest.mark.catches("ERR-073")
def test_duplicated_node_breaks_oh_invariance() -> None:
    r"""A bit-identical duplicate node destroys invariance, and must be seen.

    Appending a copy of one node to an :math:`O_h`-invariant rule leaves the
    node SET unchanged but not the MEASURE: the duplicated position now
    carries twice the mass of its mirror image, so :math:`M_\# \mu \neq \mu`.
    Nearest-neighbour matching still finds a same-weight partner for every
    image, so an injectivity-free check certifies it. The match must be
    required to be a bijection.
    """
    q = Quadrature.level_symmetric(sn_order=4)
    nodes = np.column_stack([q.mu_x, q.mu_y, q.mu_z])
    weights = np.asarray(q.weights)

    duplicated = DiscreteMeasure(
        nodes=np.vstack([nodes, nodes[0]]),
        weights=np.concatenate([weights, [weights[0]]]),
        support="S^2",
    )

    # Independent ground truth: O_h contains sigma_x, so the mass carried at
    # p0 and at sigma_x(p0) must agree. It does not.
    p0 = nodes[0]
    image = np.array([-p0[0], p0[1], p0[2]])
    dup_nodes = np.vstack([nodes, nodes[0]])
    dup_w = np.concatenate([weights, [weights[0]]])
    mass_p0 = dup_w[np.linalg.norm(dup_nodes - p0, axis=1) < 1e-12].sum()
    mass_image = dup_w[np.linalg.norm(dup_nodes - image, axis=1) < 1e-12].sum()
    assert not np.isclose(mass_p0, mass_image), (
        "fixture is not actually non-invariant -- pick a node off the mirror"
    )

    assert not SubgroupOfO3.OctahedralOh.is_invariant(duplicated), (
        f"duplicated measure certified O_h-invariant, but mass at p0 is "
        f"{mass_p0!r} against {mass_image!r} at its mirror image"
    )


@pytest.mark.foundation
def test_bijection_requirement_keeps_every_shipped_rule_certified() -> None:
    """Positive control: the bijection requirement adds no false negative.

    Each shipped rule must still certify under the group it is built to
    carry. A stricter check that rejected a genuine rule would be a
    regression, not a fix.
    """
    for sn_order in (4, 8, 16):
        mu = _measure_from_sphere_quad(Quadrature.level_symmetric(sn_order=sn_order))
        assert SubgroupOfO3.OctahedralOh.is_invariant(mu)
    for order in (5, 11, 17):
        mu = _measure_from_sphere_quad(Quadrature.lebedev(order=order))
        assert SubgroupOfO3.OctahedralOh.is_invariant(mu)
    for n_phi in (4, 8, 16):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        assert SubgroupOfO3.Dnh(n_phi).is_invariant(mu)
