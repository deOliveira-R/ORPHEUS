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
    return DiscreteMeasure(nodes=nodes, weights=q.weights, space="S^2")


def _measure_from_1d_quad(q) -> DiscreteMeasure:
    """Build a 1-D DiscreteMeasure on :math:`[-1, 1]` from a 1-D quadrature."""
    return DiscreteMeasure(
        nodes=q.mu_x, weights=q.weights, space="[-1,1]",
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
    mu = DiscreteMeasure(nodes=nodes, weights=weights, space="S^2")
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
    mu = DiscreteMeasure(nodes=nodes, weights=weights, space="S^2")
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
