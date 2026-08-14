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
    SubgroupOfO3.Mirror("z"),
    SubgroupOfO3.SO2,
    SubgroupOfO3.Dinfh,
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
        SubgroupOfO3.Mirror("z"),
        SubgroupOfO3.SO2,
        SubgroupOfO3.Dinfh,
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
    r""":math:`SO(2) \subset D_{\infty h} \subset O(3)` — the axisymmetric tower."""
    assert SubgroupOfO3.Dinfh.contains(SubgroupOfO3.SO2)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Dinfh)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.SO2)  # transitive


@pytest.mark.foundation
def test_mirror_chain_via_octahedral() -> None:
    r""":math:`\sigma_z \subset O_h \subset O(3)` — the Cartesian tower."""
    assert SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Mirror("z"))
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.OctahedralOh)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Mirror("z"))


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
def test_dnh_contains_sigma_h_at_every_n_but_sigma_x_only_at_EVEN_n() -> None:
    r"""WHICH reflection sits inside :math:`D_{nh}` is ``n``-dependent.

    ⚠ Re-posed 2026-08-02. This asserted ``Dnh(n).contains(Z2)`` for
    ``n in (1, 2, 3, 4, 6)`` under the docstring *"Z_2 (a single
    reflection) sits inside every D_nh"*. The prose was a FALSE
    GENERALISATION that the parameter-free tag made unfalsifiable: ``Z2``
    was realized as :math:`\sigma_z`, so the test only ever asked about
    the horizontal mirror — which every :math:`D_{nh}` does carry.

    `[M]` Named, the claim splits: :math:`\sigma_h` is in every
    :math:`D_{nh}`, but :math:`\sigma_x` is a VERTICAL mirror and
    :math:`D_{nh}`'s vertical mirrors sit at :math:`k\pi/n` — so the
    x-normal plane is one of them only for even ``n``. Two of the five
    orders this test was parametrized over answer False.
    """
    for n in (1, 2, 3, 4, 6):
        assert SubgroupOfO3.Dnh(n).contains(SubgroupOfO3.Mirror("z")), n
        assert SubgroupOfO3.Dnh(n).contains(SubgroupOfO3.Mirror("x")) == (
            n % 2 == 0
        ), n


@pytest.mark.foundation
def test_dnh_in_o3() -> None:
    r""":math:`D_{nh} \subset O(3)` and :math:`D_{nh} \subset D_{\infty h}`.

    The second relation is the one this test was always asserting, and it is
    TRUE — but it was written against an entry named ``O2`` whose realization
    is :math:`C_{\infty h}` (axial rotations + :math:`\sigma_h`). Neither
    :math:`C_{\infty h}` nor the true :math:`O(2) = C_{\infty v}` contains
    :math:`D_{nh}`, because :math:`D_{nh}` carries :math:`C_2` axes lying IN
    the plane. Renaming the entry to :math:`D_{\infty h}` — the full
    cylindrical group, and what a cylinder actually carries — makes the
    assertion correct rather than deleting it.
    """
    for n in (1, 2, 3, 4, 6):
        assert SubgroupOfO3.O3.contains(SubgroupOfO3.Dnh(n))
        assert SubgroupOfO3.Dinfh.contains(SubgroupOfO3.Dnh(n))
        # Not in SO(2)/SO(3) — D_nh contains improper rotations.
        assert not SubgroupOfO3.SO2.contains(SubgroupOfO3.Dnh(n))
        assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.Dnh(n))


# ============================================================================
# is_subgroup_of — readability synonym
# ============================================================================


@pytest.mark.foundation
def test_is_subgroup_of_reverses_contains() -> None:
    """``A.is_subgroup_of(B)`` is equivalent to ``B.contains(A)``."""
    assert SubgroupOfO3.Mirror("z").is_subgroup_of(SubgroupOfO3.OctahedralOh)
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
def test_gauss_legendre_1d_is_NOT_so2_invariant(n: int) -> None:
    r"""A polar marginal is **not** :math:`SO(2)`-invariant, and the
    inverted claim is what made a 1-D measure answer :math:`O(3)`.

    This gate is the INVERSION of a test that asserted the opposite
    ("trivially :math:`SO(2)`-invariant: there is no azimuthal coordinate
    to rotate"). The premise was a category error. A measure on
    :math:`[-1,1]` embeds as :math:`(\mu, 0, 0)` — the tree's canonical
    embedding, written down in :meth:`Quadrature.axis_cosines` and used by
    ``spherical_harmonics`` internally — and a rotation about :math:`z`
    **moves** :math:`(\mu, 0, 0)` off the x-axis. The old claim was only
    true for a :math:`(0, 0, \mu)` embedding, i.e. it silently used a
    different convention from the mirror branch sitting beside it.

    What "there is no azimuthal coordinate to rotate" really described is
    the group the marginal was quotiented BY. That is a real and useful
    fact, but it is a property of the REDUCTION, not of this point set,
    and it lives in ``AngularSymmetry.continuous_isotropy`` where it
    derives the support.

    The old claim also could not fail: `[M]` it held for EVERY measure on
    :math:`[-1,1]`, including a deliberately asymmetric one — certifying a
    continuous group that was never tested, which is ERR-072's shape and
    contradicted :func:`_is_axis_supported`, the exact criterion the same
    module applies to three-dimensional nodes.
    """
    q = Quadrature.gauss_legendre(n_ordinates=n)
    mu = _measure_from_1d_quad(q)
    assert not SubgroupOfO3.SO2.is_invariant(mu)
    # And the property the slab actually owes DOES hold — otherwise this
    # gate would be satisfied by a measure with no symmetry at all.
    assert SubgroupOfO3.Mirror("x").is_invariant(mu)


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
    # sigma_x, not "a reflection": the polar marginal embeds as
    # (mu, 0, 0), so mu -> -mu is the mirror normal to x. sigma_y and
    # sigma_z fix the embedded nodes POINTWISE, so they hold trivially —
    # asserted alongside, because a lone True for sigma_x would not show
    # that the plane is the thing being tested.
    assert SubgroupOfO3.Mirror("x").is_invariant(mu)
    assert SubgroupOfO3.Mirror("y").is_invariant(mu)
    assert SubgroupOfO3.Mirror("z").is_invariant(mu)


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
    assert SubgroupOfO3.Mirror("z") == SubgroupOfO3.Mirror("z")
    assert SubgroupOfO3.Mirror("x") != SubgroupOfO3.Mirror("z")
    assert SubgroupOfO3.Cn(6) == SubgroupOfO3.Cn(6)
    assert SubgroupOfO3.Cn(4) != SubgroupOfO3.Cn(6)
    assert SubgroupOfO3.Dnh(4) != SubgroupOfO3.Dnh(6)
    # Repr should not raise and should mention the group name.
    # The repr MUST carry the axis: the lattice walk keys `visited`
    # and `_GROUP_CACHE` on it, so a plane-dropping repr would
    # silently collapse the three mirrors into one entry.
    assert "Mirror('z')" in repr(SubgroupOfO3.Mirror("z"))
    assert repr(SubgroupOfO3.Mirror("x")) != repr(SubgroupOfO3.Mirror("z"))
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
    for sn_order in (4, 8, 12):
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
    assert SubgroupOfO3.Dinfh.is_invariant(poles)
    assert not SubgroupOfO3.SO3.is_invariant(poles)

    # Unequal weights break sigma_h, hence O(2), while SO(2) survives.
    lopsided = DiscreteMeasure(
        nodes=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
        weights=np.array([1.0, 2.0]),
        support="S^2",
    )
    assert SubgroupOfO3.SO2.is_invariant(lopsided)
    assert not SubgroupOfO3.Dinfh.is_invariant(lopsided)


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
    for sn_order in (4, 8, 12):
        mu = _measure_from_sphere_quad(Quadrature.level_symmetric(sn_order=sn_order))
        assert SubgroupOfO3.OctahedralOh.is_invariant(mu)
    for order in (5, 11, 17):
        mu = _measure_from_sphere_quad(Quadrature.lebedev(order=order))
        assert SubgroupOfO3.OctahedralOh.is_invariant(mu)
    for n_phi in (4, 8, 16):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        assert SubgroupOfO3.Dnh(n_phi).is_invariant(mu)


# ============================================================================
# Containment computed from the realized operator sets
# ============================================================================
#
# A hand-maintained lattice is a claim with no construction behind it, and
# this module shipped two such claims that were false. Finite-vs-finite
# containment is now COMPUTED from the same operator sets `is_invariant`
# applies, so the two verbs answer questions about the same group.
#
# The sense is LITERAL containment in the standard setting (principal axis
# along z, vertex on x). That is the question the selection gate asks, since
# a rule and its geometry are used in one frame. The setting-INDEPENDENT
# relation -- subconjugacy, exists g with g H g^-1 <= K -- is a different
# question and is not what `contains` answers.


@pytest.mark.foundation
def test_realized_group_orders_are_textbook() -> None:
    """The realizations generate groups of the right ORDER.

    If a realization were wrong, every containment computed from it would be
    wrong in the same direction and silently. Orders are the cheapest
    independent check: |D_nh| = 4n, |O_h| = 48, |I_h| = 120.
    """
    from orpheus.numerics.symmetry import _close_group, _realized_ops

    expected = {
        "Trivial": (SubgroupOfO3.Trivial, 1),
        "sigma_z": (SubgroupOfO3.Mirror("z"), 2),
        "C_3": (SubgroupOfO3.Cn(3), 3),
        "C_4": (SubgroupOfO3.Cn(4), 4),
        "D_2h": (SubgroupOfO3.Dnh(2), 8),
        "D_3h": (SubgroupOfO3.Dnh(3), 12),
        "D_4h": (SubgroupOfO3.Dnh(4), 16),
        "O_h": (SubgroupOfO3.OctahedralOh, 48),
        "I_h": (SubgroupOfO3.IcosahedralIh, 120),
    }
    for label, (tag, order) in expected.items():
        ops = _realized_ops(tag._tag)
        assert ops is not None, f"{label} has no finite realization"
        assert len(_close_group(ops)) == order, (
            f"|{label}| computed as {len(_close_group(ops))}, expected {order}"
        )


@pytest.mark.foundation
def test_octahedral_contains_exactly_its_z_axis_subgroups() -> None:
    r""":math:`O_h` contains :math:`C_1, C_2, C_4` and :math:`D_{1h},
    D_{2h}, D_{4h}` -- and no others among these families.

    ``O`` is isomorphic to :math:`S_4`, whose element orders are 1, 2, 3, 4:
    so there is no 6-fold proper rotation and :math:`C_6` is excluded under
    ANY setting. :math:`C_3` is excluded under THIS setting because the
    cube's 3-fold axis is the body diagonal, not :math:`z`.

    These relations were all ``False`` before containment was computed --
    a deliberate gap (the source said "we do not encode this until a
    consumer needs it"), whose accompanying note also mis-stated the answer
    as ``n in {1,2,3,4,6}``.
    """
    for n in (1, 2, 4):
        assert SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Cn(n)), f"C_{n}"
        assert SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Dnh(n)), f"D_{n}h"
    for n in (3, 5, 6):
        assert not SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Cn(n)), f"C_{n}"
    for n in (3, 6):
        assert not SubgroupOfO3.OctahedralOh.contains(SubgroupOfO3.Dnh(n)), f"D_{n}h"


@pytest.mark.foundation
def test_computed_containment_obeys_the_order_relation_laws() -> None:
    """Reflexive, antisymmetric, transitive -- and consistent with Lagrange.

    A computed lattice can be checked against the laws an order relation
    must satisfy; a hand-written one cannot be checked at all.
    """
    from orpheus.numerics.symmetry import _close_group, _realized_ops

    finite = [
        SubgroupOfO3.Trivial, SubgroupOfO3.Mirror("z"),
        SubgroupOfO3.Cn(1), SubgroupOfO3.Cn(2), SubgroupOfO3.Cn(3),
        SubgroupOfO3.Cn(4), SubgroupOfO3.Dnh(1), SubgroupOfO3.Dnh(2),
        SubgroupOfO3.Dnh(4), SubgroupOfO3.OctahedralOh,
    ]
    order = {
        g: len(_close_group(_realized_ops(g._tag)))  # type: ignore[arg-type]
        for g in finite
    }

    for a in finite:
        assert a.contains(a), f"not reflexive at {a.name}"

    for a in finite:
        for b in finite:
            if a.contains(b) and b.contains(a):
                assert order[a] == order[b], (
                    f"{a.name} and {b.name} contain each other but have "
                    f"orders {order[a]} != {order[b]}"
                )
            if a.contains(b):
                # Lagrange: a subgroup's order divides the group's.
                assert order[a] % order[b] == 0, (
                    f"{a.name} contains {b.name} but {order[b]} does not "
                    f"divide {order[a]} -- impossible for a subgroup"
                )

    for a in finite:
        for b in finite:
            for c in finite:
                if a.contains(b) and b.contains(c):
                    assert a.contains(c), (
                        f"transitivity broken: {a.name} > {b.name} > {c.name}"
                    )


@pytest.mark.foundation
def test_vertical_mirror_planes_follow_the_dnh_setting() -> None:
    r"""The :math:`\sigma_v` PLANES sit at :math:`k\pi/n`, not their normals.

    The standard :math:`D_{nh}` setting puts a vertex on the x-axis, so the
    :math:`\varphi = 0` plane is a mirror. Placing the k-th *normal* at
    :math:`k\pi/n` instead rotates every plane by :math:`\pi/2` -- which
    maps the set onto itself for EVEN ``n`` and is therefore invisible
    there, while giving genuinely different planes for ODD ``n``.

    Orthogonality, determinant, closure and group order all survive a
    rotated mirror set, so only comparing angles against the setting can
    see this. Pinned on odd ``n``, where the product rule is demonstrably
    closed under its :math:`\varphi = 0` mirror.
    """
    for n_phi in (3, 5, 7):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        assert SubgroupOfO3.Dnh(n_phi).is_invariant(mu), (
            f"product(4,{n_phi}) is not D_{n_phi}h-invariant -- the sigma_v "
            f"planes are rotated off the setting"
        )


# ============================================================================
# The subgroup graph, and walking it to compute a measure's symmetry
# ============================================================================
#
# Crystallography walks the Hasse diagram of the subgroup lattice downward
# from high symmetry to find the symmetry a structure actually has. Doing the
# same here replaces a DECLARED invariance group with a COMPUTED one -- a
# declaration is a claim with no construction behind it, and two such claims
# already shipped false in this module.


def _walk_rules() -> list[tuple[str, object]]:
    """A representative rule set spanning all three families.

    Deliberately not the full parameter sweep: ``_orbit_closure`` is an
    O(N^2) Python loop per operator, so a brute-force scan over every
    candidate of a 110-node Lebedev rule costs minutes. This set keeps the
    discriminating cases -- odd and even ``n_phi``, the two ``O_h`` families
    -- at a cost the foundation tier can carry on every run.
    """
    rules: list[tuple[str, object]] = []
    for n_phi in (3, 4, 8):
        rules.append((f"product(4,{n_phi})",
                      Quadrature.product(n_mu=4, n_phi=n_phi)))
    for sn_order in (4, 8):
        rules.append((f"level_symmetric({sn_order})",
                      Quadrature.level_symmetric(sn_order=sn_order)))
    for order in (5, 11):
        rules.append((f"lebedev({order})", Quadrature.lebedev(order=order)))
    return rules


@pytest.mark.foundation
def test_walk_and_bruteforce_agree() -> None:
    """The two realizations of ``Sym`` must agree on every rule.

    This is the verification instrument, not a fast-path check: the pruned
    Hasse walk and the exhaustive scan compute the same definition by
    structurally different routes, so agreement across the whole production
    family is what PROVES the walk rather than merely pinning it.
    """
    from orpheus.numerics.symmetry import maximal_invariance_groups

    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        walk = sorted(g.name for g in maximal_invariance_groups(mu, method="walk"))
        brute = sorted(
            g.name for g in maximal_invariance_groups(mu, method="bruteforce")
        )
        assert walk == brute, f"{name}: walk={walk} bruteforce={brute}"


@pytest.mark.foundation
def test_computed_symmetry_matches_the_construction() -> None:
    r"""The walk recovers each rule's symmetry from its NODES alone.

    ``product(n_mu, n_phi)`` is built with equispaced azimuth from
    :math:`\varphi = 0` over Gauss-Legendre polar levels, so it carries
    :math:`D_{n_\varphi h}` -- a group that depends on a PARAMETER.
    Level-symmetric and Lebedev rules are built :math:`O_h`-invariant.
    Nothing is declared here; the measure is asked.
    """
    from orpheus.numerics.symmetry import maximal_invariance_groups

    for n_phi in (3, 4, 5, 8, 12):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        got = maximal_invariance_groups(mu)
        assert [g.name for g in got] == [SubgroupOfO3.Dnh(n_phi).name], (
            f"product(4,{n_phi}): expected D_{n_phi}h, got "
            f"{[g.name for g in got]}"
        )

    for sn_order in (4, 8):
        mu = _measure_from_sphere_quad(
            Quadrature.level_symmetric(sn_order=sn_order)
        )
        assert [g.name for g in maximal_invariance_groups(mu)] == [
            SubgroupOfO3.OctahedralOh.name
        ]
    for order in (5, 11):
        mu = _measure_from_sphere_quad(Quadrature.lebedev(order=order))
        assert [g.name for g in maximal_invariance_groups(mu)] == [
            SubgroupOfO3.OctahedralOh.name
        ]


@pytest.mark.foundation
def test_invariance_is_downward_closed() -> None:
    r"""The law the walk's pruning rests on.

    If :math:`\mu` is :math:`G`-invariant and :math:`H \le G`, then
    :math:`\mu` is :math:`H`-invariant -- so the invariant set is an order
    ideal, a failing node kills every supergroup, and a passing node implies
    every subgroup. The walk is unsound without it.

    Simultaneously a test OF the lattice: a wrong containment edge shows up
    here as a violation. This law measured 68 violations before the lattice
    was computed from the operator sets rather than tabulated.
    """
    from orpheus.numerics.symmetry import candidate_groups

    violations: list[str] = []
    checked = 0
    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        cands = candidate_groups(mu)
        invariant = {repr(g): g.is_invariant(mu) for g in cands}
        for outer in cands:
            for inner in cands:
                if outer == inner or not outer.contains(inner):
                    continue
                checked += 1
                if invariant[repr(outer)] and not invariant[repr(inner)]:
                    violations.append(
                        f"{name}: {outer.name}-invariant but not "
                        f"{inner.name}-invariant, though {inner.name} <= "
                        f"{outer.name}"
                    )
    assert checked > 500, f"only {checked} pairs checked -- gate is too thin"
    assert not violations, (
        f"{len(violations)} downward-closure violations:\n  "
        + "\n  ".join(violations[:10])
    )


@pytest.mark.foundation
def test_every_registry_rule_declares_a_symmetry_it_actually_has() -> None:
    r"""Every shipped rule's ``invariance_group`` is a TRUE claim.

    This gate replaces (2026-08-02) a characterization test that pinned the
    opposite: the product rule declared :math:`SO(2)` while the walk said
    :math:`D_{n_\varphi h}`, and that test existed to hold the contradiction
    visible until the registry was re-posed. Its docstring predicted it would
    "start passing for the opposite reason once the registry is re-posed" --
    this is that re-posing, so the gate is inverted rather than deleted.

    The claim tested is deliberately WEAKER than equality with the computed
    maximal group, and the justification is the CONTRACT, not an example.
    ``DiscreteMeasure.invariance_group`` is a claim the measure makes about
    its own nodes; this gate exists to forbid a claim the nodes do NOT
    satisfy. It must not also forbid a claim that is true but modest, because
    the type permits that and a gate may not assert more than the type
    promises (``vv-principles`` #16). Strengthening to equality would convict
    a legitimate value.

    ⭐ ``gauss_legendre_on_mu`` shows why equality is not even well-POSED
    here: `[M]` its maximal invariance groups are
    :math:`\{\sigma_x, \sigma_y, \sigma_z\}` -- **three** of them, because the
    polar marginal embeds as :math:`(\mu, 0, 0)` so :math:`\sigma_y` and
    :math:`\sigma_z` fix every node pointwise. There is no single "the
    maximal group" to compare against, so an equality gate would have to pick
    one arbitrarily.

    ⛔ This paragraph justified the weakness differently until 2026-08-14:
    *"one shipped rule genuinely is [true-but-not-maximal]:
    ``gauss_legendre_on_mu`` declares* :math:`SO(2)`\ *, the group its domain
    was quotiented BY."* Both halves are now false. `[M]`
    ``gauss_legendre_on_mu(8).invariance_group`` is
    ``SubgroupOfO3.Mirror('x')`` -- the 2026-08-02 correction that replaced
    the slab/sphere residual :math:`Z_2` with :math:`\sigma_x` reached the
    rule's declared tag too -- and `[M]` walking
    :func:`~orpheus.numerics.symmetry.maximal_invariance_groups` over all four
    registry rules, **0 of 4** are true-but-not-maximal (GL declares
    :math:`\sigma_x`, one of its three; Lebedev and LS_N declare
    :math:`O_h`; the product rule declares :math:`D_{6h}`).
    ⟹ the gate is unchanged and still correct; what was repaired is its
    EVIDENCE. A weakness argued from a single example silently becomes
    unjustified the moment that example moves, and the shape of this failure
    is that nothing goes red -- the gate kept passing for the whole time its
    stated reason was false.
    """
    from orpheus.numerics.quadrature.registry import quadrature_registry

    checked: list[str] = []
    for spec in quadrature_registry:
        params = spec.degree_of_exactness_for(5)
        assert params is not None, f"{spec.name} cannot reach degree 5"
        measure = spec.build(params)
        declared = measure.invariance_group
        assert declared is not None, (
            f"{spec.name} carries no invariance_group; DiscreteMeasure.phase "
            f"reads it to classify the measure as angular"
        )
        assert declared.is_invariant(measure), (
            f"{spec.name} at {params} declares {declared.name}, but its own "
            f"nodes are not invariant under it -- a declared group with no "
            f"construction behind it"
        )
        checked.append(spec.name)

    assert set(checked) == {
        "GaussLegendre1D",
        "LebedevSphere",
        "LevelSymmetricSN",
        "ProductQuadrature",
    }


@pytest.mark.foundation
def test_product_declares_the_group_the_walk_computes() -> None:
    r"""The product rule's declaration is now the computed answer itself.

    The one rule whose group is parameter-dependent is the one that lied, and
    it lied because the spec field's type could not spell the dependence. With
    the declaration moved onto the measure -- where ``n_phi`` is in scope --
    declared and computed agree at every order, including the odd ones the
    selection gate now rejects for a cylinder.
    """
    from orpheus.numerics.symmetry import maximal_invariance_groups

    for n_phi in (2, 3, 4, 5, 6, 8):
        mu = _measure_from_sphere_quad(
            Quadrature.product(n_mu=4, n_phi=n_phi)
        )
        computed = maximal_invariance_groups(mu)
        assert [g.name for g in computed] == [SubgroupOfO3.Dnh(n_phi).name], (
            f"n_phi={n_phi}: walk says {[g.name for g in computed]}"
        )
        assert not SubgroupOfO3.SO2.is_invariant(mu), (
            f"n_phi={n_phi}: no finite point set on S^2 is SO(2)-closed"
        )


# ============================================================================
# The orbit certificate and the singular set Sigma
# ============================================================================
#
# `_orbit_closure` always computed pi_M(i) with M x_i = x_{pi(i)} and threw it
# away to return a bool. A `-> bool` predicate that internally builds the
# permutation IS the missing primitive (L-013). Widening the return type gives
# Sigma for free, because pi_M(i) == i MEANS x_i is in Fix(M).


def _cert_rules() -> list[tuple[str, object]]:
    return [
        ("product(4,4)", Quadrature.product(n_mu=4, n_phi=4)),
        ("product(4,8)", Quadrature.product(n_mu=4, n_phi=8)),
        ("level_symmetric(8)", Quadrature.level_symmetric(sn_order=8)),
        ("lebedev(11)", Quadrature.lebedev(order=11)),
    ]


@pytest.mark.foundation
def test_singular_set_membership_is_exact_not_thresholded() -> None:
    r""":math:`\Sigma` does not move when the tolerance moves.

    Membership is :math:`\pi_M(i) = i`, an integer identity. The ad-hoc
    float comparisons this replaces would each shift with their epsilon;
    this must not, across four orders of magnitude.
    """
    from orpheus.numerics.symmetry import maximal_invariance_groups, singular_set

    for name, q in _cert_rules():
        mu = _measure_from_sphere_quad(q)
        group = maximal_invariance_groups(mu)[0]
        answers = {
            atol: tuple(singular_set(mu, group, atol=atol).tolist())
            for atol in (1e-15, 1e-13, 1e-11)
        }
        assert len(set(answers.values())) == 1, (
            f"{name}: Sigma varies with atol -- {[(k, len(v)) for k, v in answers.items()]}"
        )


@pytest.mark.foundation
def test_singular_set_under_d2h_reproduces_the_epsilon_detectors() -> None:
    r"""The retirement warrant for the three ad-hoc epsilon detectors.

    ``_OCTANT_SIGN_EPS``, ``_MU_DIRECTION_EPS`` and ``TANGENTIAL_EPS`` all ask
    "is a direction cosine zero", i.e. does this node lie on a COORDINATE
    mirror. That is :math:`\Sigma` under :math:`D_{2h}` — the group generated
    by the three coordinate reflections, isomorphic to
    :math:`(\mathbb{Z}_2)^3`, whose chambers are exactly the octants. Naming
    the group is what makes the question precise; the epsilons could not say
    WHICH mirrors they meant.

    Exact agreement on every production rule is what licenses retiring them.
    """
    from orpheus.numerics.symmetry import singular_set

    for name, q in _cert_rules() + [
        ("lebedev(17)", Quadrature.lebedev(order=17)),
        # S12, the top of the family's defined range since #327 (the
        # moment-matched solve has no positive solution above it). The row
        # wants a HIGH order with many levels, not the number 16.
        ("level_symmetric(12)", Quadrature.level_symmetric(sn_order=12)),
    ]:
        mu = _measure_from_sphere_quad(q)
        nodes = np.column_stack([q.mu_x, q.mu_y, q.mu_z])
        by_epsilon = np.flatnonzero((np.abs(nodes) < 1e-15).any(axis=1))
        by_group = singular_set(mu, SubgroupOfO3.Dnh(2))
        assert np.array_equal(np.sort(by_group), np.sort(by_epsilon)), (
            f"{name}: Sigma|D_2h has {len(by_group)} points, the epsilon "
            f"detector {len(by_epsilon)}"
        )


@pytest.mark.foundation
def test_full_group_sigma_sees_diagonal_mirrors_the_epsilons_cannot() -> None:
    r"""And why the retirement must NAME :math:`D_{2h}`, not just swap spelling.

    :math:`\Sigma` under the rule's FULL symmetry group is strictly larger
    whenever that group carries mirrors off the coordinate planes. A node on
    the diagonal plane :math:`x = y` is singular with no zero coordinate, so
    a cosine-magnitude test is structurally blind to it. Most starkly:
    ``level_symmetric`` has NO node with a zero cosine, so the epsilon answer
    is empty while 32 of its ordinates lie on :math:`O_h` mirrors.
    """
    from orpheus.numerics.symmetry import maximal_invariance_groups, singular_set

    mu = _measure_from_sphere_quad(Quadrature.level_symmetric(sn_order=8))
    coordinate_only = singular_set(mu, SubgroupOfO3.Dnh(2))
    full = singular_set(mu, maximal_invariance_groups(mu)[0])
    assert len(coordinate_only) == 0
    assert len(full) > 0, "expected O_h diagonal mirrors to fix some ordinate"
    assert set(coordinate_only.tolist()) <= set(full.tolist())

    # product(4, n_phi) has n_phi vertical mirrors, only 2 of them coordinate
    # planes, so Sigma grows with n_phi while the epsilon answer does not.
    sizes = {}
    for n_phi in (4, 8, 16):
        m = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        sizes[n_phi] = (
            len(singular_set(m, SubgroupOfO3.Dnh(2))),
            len(singular_set(m, maximal_invariance_groups(m)[0])),
        )
    assert [s[0] for s in sizes.values()] == [16, 16, 16], sizes
    assert [s[1] for s in sizes.values()] == [16, 32, 64], sizes


@pytest.mark.foundation
def test_singular_set_requires_an_invariant_measure() -> None:
    r""":math:`\Sigma` is unrepresentable without the closure proof.

    A quotient needs something to quotient: the singular set is defined only
    on a :math:`G`-invariant measure. Since the certificate exists only when
    closure holds, the precondition is enforced by construction rather than
    documented in a comment.
    """
    from orpheus.numerics.symmetry import singular_set

    lopsided = DiscreteMeasure(
        nodes=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        weights=np.array([1.0, 2.0]),
        support="S^2",
    )
    with pytest.raises(ValueError, match="invariant"):
        singular_set(lopsided, SubgroupOfO3.OctahedralOh)


@pytest.mark.foundation
def test_orbit_stabilizer_theorem_holds_on_every_certificate() -> None:
    r"""The certificate's defining law: :math:`|Gx|\cdot|\mathrm{Stab}(x)| = |G|`.

    An intrinsic property of the object, not of its use — if the
    permutations were not a genuine group action this would fail, so it
    independently re-checks the bijection requirement (ERR-073) on every
    element rather than only on the duplicate-node witness.
    """
    from orpheus.numerics.symmetry import (
        _close_group, _realized_ops, maximal_invariance_groups, orbit_certificate,
    )

    for name, q in _cert_rules():
        mu = _measure_from_sphere_quad(q)
        group = maximal_invariance_groups(mu)[0]
        cert = orbit_certificate(mu, group)
        assert cert is not None, name
        order = len(_close_group(_realized_ops(group._tag)))  # type: ignore[arg-type]
        stab = cert.stabilizer_order
        for orbit in cert.orbits():
            for i in orbit:
                assert len(orbit) * int(stab[i]) == order, (
                    f"{name}: orbit-stabilizer fails at node {i} -- "
                    f"|orbit|={len(orbit)} |Stab|={stab[i]} |G|={order}"
                )

        # Stabilizer order is 1 exactly off Sigma.
        sigma = set(cert.singular_set.tolist())
        for i in range(cert.n_points):
            assert (int(stab[i]) > 1) is (i in sigma), (
                f"{name}: node {i} has |Stab|={stab[i]} but "
                f"{'is' if i in sigma else 'is not'} in Sigma"
            )


# ============================================================================
# Proper vs improper, and the 1-D polar-marginal action
# ============================================================================


@pytest.mark.foundation
def test_a_mirror_is_improper_so_it_is_not_inside_the_rotation_groups() -> None:
    r""":math:`Z_2` is realized as :math:`\sigma_z`, ``det = -1``.

    An improper element cannot lie in a proper-rotation group, so
    :math:`Z_2 \not\le SO(2)` and :math:`Z_2 \not\le SO(3)`. The lattice
    asserted the second until 2026-08-02, and it broke monotonicity on any
    measure that is :math:`SO(3)`-invariant without being
    reflection-symmetric.

    The proper order-2 sibling is ``Cn(2)``, which IS inside both — that
    distinction is the reason the two spellings exist, and it is exactly
    what the false edge erased.
    """
    from orpheus.numerics.symmetry import _close_group, _realized_ops

    dets = [
        M.determinant
        for M in _close_group(_realized_ops(SubgroupOfO3.Mirror("z")._tag))  # type: ignore[arg-type]
    ]
    assert sorted(round(d) for d in dets) == [-1, 1], dets

    assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.Mirror("z"))
    assert not SubgroupOfO3.SO2.contains(SubgroupOfO3.Mirror("z"))
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Mirror("z"))
    # The proper sibling behaves oppositely.
    assert SubgroupOfO3.SO3.contains(SubgroupOfO3.Cn(2))
    assert SubgroupOfO3.SO2.contains(SubgroupOfO3.Cn(2))


@pytest.mark.foundation
def test_so3_on_a_polar_marginal_requires_reflection_symmetry() -> None:
    r""":math:`SO(3)` does not act trivially on a 1-D :math:`\mu`-measure.

    :math:`SO(2)` and :math:`C_n` rotate about :math:`z` and leave the polar
    cosine alone, so a 1-D measure is trivially invariant under them. But
    :math:`R_x(\pi) \in SO(3)` induces :math:`\mu \to -\mu`, so an
    :math:`SO(3)`-invariant polar marginal must be reflection-symmetric.

    Returning ``True`` unconditionally made an asymmetric node set read
    :math:`SO(3)`-invariant and :math:`Z_2`-non-invariant at once.

    **Re-posed 2026-08-03, and the premise in the paragraph above is now
    part of the history rather than the claim.** ":math:`SO(2)` and
    :math:`C_n` rotate about :math:`z` and leave the polar cosine alone"
    is true only under a :math:`(0,0,\mu)` embedding; the tree embeds a
    polar marginal as :math:`(\mu, 0, 0)`, where a rotation about
    :math:`z` **moves** the node. The two readings sat in adjacent
    branches of one function, and the rotational one is what let a 1-D
    measure answer :math:`O(3)`.

    So all three rotational groups now read ``False`` on both fixtures,
    for two different exact reasons — :math:`SO(2)` and :math:`SO(3)`
    because a finite set closed under a CONTINUOUS group must be
    axis-/origin-supported, and :math:`C_n` because it genuinely moves
    :math:`(\mu,0,0)`. The discriminating content moves to
    :math:`\sigma_x`, which is where it always belonged: it is what the
    slab and sphere rows of ``GEOMETRY_ANGULAR_SYMMETRY`` require.
    """
    asymmetric = DiscreteMeasure(
        nodes=np.array([-0.9, -0.1, 0.3]),
        weights=np.ones(3),
        support="[-1,1]",
    )
    symmetric = DiscreteMeasure(
        nodes=np.array([-0.6, -0.2, 0.2, 0.6]),
        weights=np.ones(4),
        support="[-1,1]",
    )

    # No rotational group is invariant on a polar marginal, symmetric or
    # not — and crucially these are NOT vacuous "always False" rows: the
    # sigma_x pair below shows the same fixtures ARE distinguishable.
    for rotational in (SubgroupOfO3.SO2, SubgroupOfO3.Cn(3), SubgroupOfO3.SO3):
        assert not rotational.is_invariant(asymmetric), rotational
        assert not rotational.is_invariant(symmetric), rotational

    # The reflection is the discriminator, and it separates the fixtures.
    assert not SubgroupOfO3.Mirror("x").is_invariant(asymmetric)
    assert SubgroupOfO3.Mirror("x").is_invariant(symmetric)


@pytest.mark.foundation
def test_invariance_is_downward_closed_on_polar_marginals() -> None:
    """The monotonicity law again, on the 1-D path.

    The sphere-side gate cannot see 1-D bugs: that path never calls
    ``_orbit_closure`` and has its own dispatch. Both defects fixed above
    surface here as violations.
    """
    tags = [
        # All THREE mirrors. On a polar marginal (embedded (mu, 0, 0))
        # sigma_y and sigma_z hold trivially — they fix every node
        # pointwise — so sigma_x is the only one that can constrain
        # anything here, and listing only sigma_z would leave the
        # reflection family untested on this path.
        SubgroupOfO3.Mirror("x"), SubgroupOfO3.Mirror("y"),
        SubgroupOfO3.Trivial, SubgroupOfO3.Mirror("z"), SubgroupOfO3.SO2,
        SubgroupOfO3.Dinfh, SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
        SubgroupOfO3.Cn(1), SubgroupOfO3.Cn(2), SubgroupOfO3.Cn(3),
        SubgroupOfO3.Dnh(1), SubgroupOfO3.Dnh(2), SubgroupOfO3.Dnh(4),
    ]
    measures = {
        "asymmetric": DiscreteMeasure(
            nodes=np.array([-0.9, -0.1, 0.3]), weights=np.ones(3),
            support="[-1,1]",
        ),
        "symmetric": DiscreteMeasure(
            nodes=np.array([-0.6, -0.2, 0.2, 0.6]), weights=np.ones(4),
            support="[-1,1]",
        ),
        "symmetric-nodes-unequal-weights": DiscreteMeasure(
            nodes=np.array([-0.6, 0.6]), weights=np.array([1.0, 2.0]),
            support="[-1,1]",
        ),
    }

    violations: list[str] = []
    checked = 0
    for label, mu in measures.items():
        invariant = {repr(t): t.is_invariant(mu) for t in tags}
        for outer in tags:
            for inner in tags:
                if outer == inner or not outer.contains(inner):
                    continue
                checked += 1
                if invariant[repr(outer)] and not invariant[repr(inner)]:
                    violations.append(
                        f"{label}: {outer.name}-invariant but not "
                        f"{inner.name}-invariant"
                    )
    assert checked > 150, f"only {checked} pairs checked"
    assert not violations, (
        f"{len(violations)} violations:\n  " + "\n  ".join(violations[:10])
    )


# ============================================================================
# The mirror family — Q5.0.2 (naming the plane)
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_a_mirror_reaches_a_REAL_arm_for_every_continuous_outer_group(
    axis: str,
) -> None:
    r"""``_contains`` ends in a bare ``return False``, so a tag it has no
    arm for gets a wrong-but-SILENT answer rather than an error.

    `[M]` That is not hypothetical: a prototype mirror with only the
    realization and invariance arms extended answered
    ``O3.contains(Mirror('x')) = False`` and
    ``Dinfh.contains(Mirror('x')) = False`` — both flatly wrong, and
    ``O(3) ⊉`` anything breaks the soundness precondition
    :func:`~orpheus.numerics.symmetry.maximal_invariance_groups` states
    for the lattice walk.

    The finite outer groups are decided by computed matrix containment
    and cannot go silently wrong this way; it is exactly the CONTINUOUS
    ones that fall through to the table, and the table is typed
    enum-to-enum so a parameterised tag can never be in it.
    """
    sigma = SubgroupOfO3.Mirror(axis)
    # Improper (det = -1) ⇒ inside the full groups, outside the proper ones.
    assert SubgroupOfO3.O3.contains(sigma)
    assert SubgroupOfO3.Dinfh.contains(sigma)
    assert not SubgroupOfO3.SO3.contains(sigma)
    assert not SubgroupOfO3.SO2.contains(sigma)


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_a_mirror_has_exactly_two_subgroups(axis: str) -> None:
    """:math:`\\{e, \\sigma\\}` has order 2, so its only proper subgroup is
    the trivial one — and no other mirror is inside it."""
    sigma = SubgroupOfO3.Mirror(axis)
    assert sigma.contains(SubgroupOfO3.Trivial)
    assert sigma.contains(sigma)
    for other in ("x", "y", "z"):
        if other != axis:
            assert not sigma.contains(SubgroupOfO3.Mirror(other))
    assert not sigma.contains(SubgroupOfO3.SO2)


@pytest.mark.foundation
def test_the_axis_is_validated_at_CONSTRUCTION() -> None:
    """There is no unnamed reflection — which is what the retired tag
    pretended there was."""
    with pytest.raises(ValueError, match="axis in x/y/z"):
        SubgroupOfO3.Mirror("w")


@pytest.mark.foundation
def test_the_two_invariance_ARMS_agree_on_the_canonical_embedding() -> None:
    r"""⭐ The reason the plane had to be named.

    A polar marginal is dispatched to the 1-D arm; the SAME data
    embedded as :math:`(\mu, 0, 0)` — the tree's canonical embedding,
    written down in ``Quadrature.axis_cosines`` — is dispatched to the
    3-D arm. Before the plane was named the two arms asked DIFFERENT
    questions (plane-free :math:`x \to -x` vs :math:`\sigma_z`), and
    `[M]` on an asymmetric :math:`\mu`-set the 3-D arm CERTIFIED a set
    that violates :math:`\mu \to -\mu` — a false certification in the
    dangerous direction, the ERR-072 family.

    The asymmetric leg is the discriminating one; the symmetric leg is
    the control that keeps the agreement from being vacuous (on a
    symmetric set every plane answers True and the two arms would agree
    even while asking different questions — which is precisely how this
    went unnoticed).
    """
    weights = np.full(4, 0.5)

    def both_arms(mu: np.ndarray, axis: str) -> tuple[bool, bool]:
        one_d = DiscreteMeasure(nodes=mu, weights=weights, support="[-1,1]")
        embedded = DiscreteMeasure(
            nodes=np.column_stack([mu, np.zeros(4), np.zeros(4)]),
            weights=weights,
            support="S^2",
        )
        g = SubgroupOfO3.Mirror(axis)
        return g.is_invariant(one_d), g.is_invariant(embedded)

    asymmetric = np.array([-0.9, -0.3, 0.3, 0.7])   # violates mu -> -mu
    symmetric = np.array([-0.7, -0.3, 0.3, 0.7])

    # DISCRIMINATING: sigma_x must be False on both arms, y/z True on both.
    assert both_arms(asymmetric, "x") == (False, False)
    assert both_arms(asymmetric, "y") == (True, True)
    assert both_arms(asymmetric, "z") == (True, True)

    # CONTROL: on a symmetric set every plane is True on both arms, so
    # the agreement above is not an artefact of everything being True.
    for axis in ("x", "y", "z"):
        assert both_arms(symmetric, axis) == (True, True)


@pytest.mark.foundation
def test_a_shipped_rule_already_discriminates_the_planes() -> None:
    r"""`[M]` ``product(4, 3)`` is :math:`\sigma_z`-closed and NOT
    :math:`\sigma_x`-closed.

    This is what falsifies the retired tag's docstring claim that "any
    single reflection works; the choice is convention". No synthetic
    fixture is needed — an odd-:math:`n_\varphi` product rule the tree
    ships already tells the planes apart, and the parameter-free tag
    answered ``True``.
    """
    measure = Quadrature.product(n_mu=4, n_phi=3).measure
    assert not SubgroupOfO3.Mirror("x").is_invariant(measure)
    assert SubgroupOfO3.Mirror("y").is_invariant(measure)
    assert SubgroupOfO3.Mirror("z").is_invariant(measure)
