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

from orpheus.numerics.manifold import COSINE_INTERVAL, SPHERE
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.numerics.quadrature import Quadrature

# 2.2b gates (test-architect, 2026-09-02) — the draft's imports, verbatim
import dataclasses
import itertools

import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    SPHERE,
    Ball,
    Quotient,
    barycentre,
    quotient_onto,
)
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import GEOMETRY_ANGULAR_SYMMETRY
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
)
from orpheus.numerics.symmetry import (IdentityComponent, Realization, SubgroupOfO3, _ELEMENT_ATOL,)
from orpheus.numerics.invariance import candidate_groups


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
    return DiscreteMeasure(nodes=nodes, weights=q.weights, support=SPHERE)


def _measure_from_1d_quad(q) -> DiscreteMeasure:
    """Build a 1-D DiscreteMeasure on :math:`[-1, 1]` from a 1-D quadrature."""
    return DiscreteMeasure(
        nodes=q.mu_x, weights=q.weights, support=COSINE_INTERVAL,
    )


# ============================================================================
# Containment lattice — named entries
# ============================================================================


_NAMED = [
    SubgroupOfO3.Trivial,
    SubgroupOfO3.Mirror("z"),
    SubgroupOfO3.SO2("z"),
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
        SubgroupOfO3.SO2("x"),
        SubgroupOfO3.SO2("z"),
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
    r""":math:`SO(2)_z \subset D_{\infty h} \subset O(3)` — the axisymmetric
    tower, about the axis :math:`D_{\infty h}` is realized on. About any
    other axis the middle link breaks and the outer one holds: a rotation
    about :math:`x` is a proper rotation (so it is in :math:`O(3)`) that does
    not preserve the :math:`z`-axis (so it is not in :math:`D_{\infty h}`)."""
    assert SubgroupOfO3.Dinfh.contains(SubgroupOfO3.SO2("z"))
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Dinfh)
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.SO2("z"))  # transitive
    for other in ("x", "y"):
        assert not SubgroupOfO3.Dinfh.contains(SubgroupOfO3.SO2(other))
        assert SubgroupOfO3.O3.contains(SubgroupOfO3.SO2(other))
        assert SubgroupOfO3.SO3.contains(SubgroupOfO3.SO2(other))


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
    r""":math:`SO(2)_z = \bigcup_n C_n`, so every :math:`C_n` is inside —
    and inside no other axis's rotation group, since :math:`C_n` is
    realized about :math:`z`. Except :math:`C_1 = \{e\}`, which is inside
    everything: ⛔ until 2026-09-02 the tabulated axial arm answered
    ``SO2('x') ⊉ C_1`` while ``SO2('x') ⊇ {e}`` — one group under two
    spellings, two answers — and this test pinned the wrong one; the
    relation is now COMPUTED from the realization (since R1 of #434,
    ``Realization.contains`` — until then ``_fixes_axis``) and cannot
    disagree with itself (#432).  The ``n = 1`` row is therefore a
    ``Trivial``-containment row since R1 (``Cn(1)`` IS ``Trivial``); that is
    the merge, not a gap."""
    assert SubgroupOfO3.SO2("z").contains(SubgroupOfO3.Cn(n))
    assert SubgroupOfO3.O2("z").contains(SubgroupOfO3.Cn(n))
    for other in ("x", "y"):
        assert SubgroupOfO3.SO2(other).contains(SubgroupOfO3.Cn(n)) is (n == 1)
        assert SubgroupOfO3.O2(other).contains(SubgroupOfO3.Cn(n)) is (n == 1)


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
        assert not SubgroupOfO3.SO2("z").contains(SubgroupOfO3.Dnh(n))
        assert not SubgroupOfO3.SO3.contains(SubgroupOfO3.Dnh(n))


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
    assert mu.is_invariant_under(SubgroupOfO3.OctahedralOh)


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
    assert not mu.is_invariant_under(SubgroupOfO3.IcosahedralIh)


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
    assert mu.is_invariant_under(SubgroupOfO3.OctahedralOh)


# ============================================================================
# Invariance — 1-D Gauss-Legendre: invariant under its OWN axial group only
# ============================================================================


@pytest.mark.foundation
@pytest.mark.parametrize("n", [4, 8, 16])
def test_a_polar_marginal_is_invariant_under_its_OWN_axial_group_only(
    n: int,
) -> None:
    r"""A polar marginal about :math:`x` is :math:`SO(2)_x`-invariant and
    :math:`SO(2)_z`-NON-invariant — the two answers the bare ``SO2`` tag
    could not tell apart, and the whole reason the group names its axis.

    ⚠ This gate has been INVERTED twice, and both earlier forms were
    right about what they measured. The first asserted the marginal is
    *"trivially SO(2)-invariant: there is no azimuthal coordinate to
    rotate"* — true of the group the marginal was quotiented BY, false as
    a node-set statement under the retired z-realized ``SO2``. The second
    asserted it is **not** :math:`SO(2)`-invariant — true, because that
    ``SO2`` was about :math:`z` and a rotation about :math:`z` moves
    :math:`(\mu, 0, 0)` off the x-axis, and it recorded the reason as an
    axis-convention collision (the angular-spaces plan's Part IV obstacle
    1). Tracker 2.4 (2026-09-01) resolved the collision by naming the axis
    rather than by fiat, and the derivation is one line: a finite point
    set is :math:`SO(2)_a`-closed iff every node lies ON axis :math:`a`
    (:meth:`IdentityComponent.fixes` — until R1 of #434 ``_is_axis_supported``);
    the marginal embeds along :math:`x`; so
    the verdict is ``True`` for :math:`a = x` and ``False`` otherwise.

    Both the chart-level measure (a bare interval, embedded along column
    0 by convention) and the slab's DECLARED measure (support
    ``S^2/O2_x``, embedded along the axis its orbit space names) answer
    the same way — and a marginal declared about :math:`z` answers the
    opposite way, which is what makes the axis load-bearing rather than
    decorative.

    The old "cannot fail" defect is still guarded one section down:
    ``test_invariance_is_exact_on_polar_marginals_not_vacuous`` shows the
    same fixtures ARE distinguishable, by :math:`\sigma_x`.
    """
    from orpheus.numerics.quadrature import gauss_legendre_on_polar_orbit

    q = Quadrature.gauss_legendre(n_ordinates=n)
    for mu in (_measure_from_1d_quad(q), q.measure):
        assert mu.is_invariant_under(SubgroupOfO3.SO2("x"))
        assert not mu.is_invariant_under(SubgroupOfO3.SO2("y"))
        assert not mu.is_invariant_under(SubgroupOfO3.SO2("z"))
        # And the property the slab actually owes DOES hold — otherwise
        # this gate would be satisfied by a measure with no symmetry at all.
        assert mu.is_invariant_under(SubgroupOfO3.Mirror("x"))

    about_z = gauss_legendre_on_polar_orbit(n, axis="z")
    assert about_z.is_invariant_under(SubgroupOfO3.SO2("z"))
    assert not about_z.is_invariant_under(SubgroupOfO3.SO2("x"))
    assert about_z.is_invariant_under(SubgroupOfO3.Mirror("z"))


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
    assert mu.is_invariant_under(SubgroupOfO3.Mirror("x"))
    assert mu.is_invariant_under(SubgroupOfO3.Mirror("y"))
    assert mu.is_invariant_under(SubgroupOfO3.Mirror("z"))


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
    mu = DiscreteMeasure(nodes=nodes, weights=weights, support=SPHERE)
    assert mu.is_invariant_under(SubgroupOfO3.Trivial)


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
    mu = DiscreteMeasure(nodes=nodes, weights=weights, support=SPHERE)
    assert not mu.is_invariant_under(SubgroupOfO3.OctahedralOh)


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
    # The repr MUST carry the axis: it is what a reader sees in a failure
    # message, and the three mirrors are three groups.  (Until R1 of #434
    # the lattice walk keyed `visited` and a `_GROUP_CACHE` on it; both now
    # key on the VALUE, whose dataclass `__eq__`/`__hash__` read the tag.)
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
        for axis in ("x", "y", "z"):
            assert not mu.is_invariant_under(SubgroupOfO3.SO2(axis)), (
                f"{name} certified SO(2)_{axis}-invariant, but a finite point "
                f"set with off-axis nodes has infinite SO(2) orbits and cannot "
                f"be closed"
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
        n_phi: _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi)).is_invariant_under(SubgroupOfO3.SO2("z"))
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
        nodes=verts, weights=np.full(len(verts), 1.0), support=SPHERE,
    )
    # Sanity: the fixture really is the I_h orbit it claims to be.
    assert len(verts) == 12
    assert ico.is_invariant_under(SubgroupOfO3.IcosahedralIh), (
        "fixture is not I_h-invariant -- it cannot discriminate the criteria"
    )
    assert not ico.is_invariant_under(SubgroupOfO3.SO3), (
        "the 12 icosahedron vertices certified SO(3)-invariant; a finite set "
        "cannot be, since every SO(3) orbit of a non-origin point is S^2"
    )
    assert not ico.is_invariant_under(SubgroupOfO3.O3)


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
        assert not mu.is_invariant_under(SubgroupOfO3.SO3), f"{name} certified SO(3)"
        assert not mu.is_invariant_under(SubgroupOfO3.O3), f"{name} certified O(3)"


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
        support=SPHERE,
    )
    assert poles.is_invariant_under(SubgroupOfO3.SO2("z"))
    assert poles.is_invariant_under(SubgroupOfO3.Dinfh)
    assert not poles.is_invariant_under(SubgroupOfO3.SO3)
    # The axis discriminates: the z-poles are moved by a rotation about x.
    assert not poles.is_invariant_under(SubgroupOfO3.SO2("x"))
    assert not poles.is_invariant_under(SubgroupOfO3.SO2("y"))

    # Unequal weights break sigma_h, hence O(2), while SO(2) survives.
    lopsided = DiscreteMeasure(
        nodes=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
        weights=np.array([1.0, 2.0]),
        support=SPHERE,
    )
    assert lopsided.is_invariant_under(SubgroupOfO3.SO2("z"))
    assert not lopsided.is_invariant_under(SubgroupOfO3.Dinfh)


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
        assert mu.is_invariant_under(SubgroupOfO3.Cn(n_phi)), (
            f"product(4,{n_phi}) is not C_{n_phi}-invariant"
        )
        assert mu.is_invariant_under(SubgroupOfO3.Dnh(n_phi)), (
            f"product(4,{n_phi}) is not D_{n_phi}h-invariant"
        )
        assert not mu.is_invariant_under(SubgroupOfO3.Cn(2 * n_phi)), (
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
        support=SPHERE,
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

    assert not duplicated.is_invariant_under(SubgroupOfO3.OctahedralOh), (
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
        assert mu.is_invariant_under(SubgroupOfO3.OctahedralOh)
    for order in (5, 11, 17):
        mu = _measure_from_sphere_quad(Quadrature.lebedev(order=order))
        assert mu.is_invariant_under(SubgroupOfO3.OctahedralOh)
    for n_phi in (4, 8, 16):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        assert mu.is_invariant_under(SubgroupOfO3.Dnh(n_phi))


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
        SubgroupOfO3.Cn(2), SubgroupOfO3.Cn(3),
        SubgroupOfO3.Cn(4), SubgroupOfO3.Dnh(1), SubgroupOfO3.Dnh(2),
        SubgroupOfO3.Dnh(4), SubgroupOfO3.OctahedralOh,
    ]
    # Nine DISTINCT groups: `Cn(1)` IS `Trivial` since R1 of #434, and a
    # roster listing one group twice would run the laws over a dict that
    # silently holds one key fewer than the list (vv-principles #20).
    assert len(set(finite)) == len(finite) == 9
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
        assert mu.is_invariant_under(SubgroupOfO3.Dnh(n_phi)), (
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

    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        walk = sorted(g.name for g in mu.symmetry_groups(method="walk"))
        brute = sorted(
            g.name for g in mu.symmetry_groups(method="bruteforce")
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

    for n_phi in (3, 4, 5, 8, 12):
        mu = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        got = mu.symmetry_groups()
        assert [g.name for g in got] == [SubgroupOfO3.Dnh(n_phi).name], (
            f"product(4,{n_phi}): expected D_{n_phi}h, got "
            f"{[g.name for g in got]}"
        )

    for sn_order in (4, 8):
        mu = _measure_from_sphere_quad(
            Quadrature.level_symmetric(sn_order=sn_order)
        )
        assert [g.name for g in mu.symmetry_groups()] == [
            SubgroupOfO3.OctahedralOh.name
        ]
    for order in (5, 11):
        mu = _measure_from_sphere_quad(Quadrature.lebedev(order=order))
        assert [g.name for g in mu.symmetry_groups()] == [
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
    from orpheus.numerics.invariance import candidate_groups

    violations: list[str] = []
    checked = 0
    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        cands = candidate_groups(mu)
        invariant = {repr(g): mu.is_invariant_under(g) for g in cands}
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
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.symmetry_groups` over all four
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
        assert measure.is_invariant_under(declared), (
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

    for n_phi in (2, 3, 4, 5, 6, 8):
        mu = _measure_from_sphere_quad(
            Quadrature.product(n_mu=4, n_phi=n_phi)
        )
        computed = mu.symmetry_groups()
        assert [g.name for g in computed] == [SubgroupOfO3.Dnh(n_phi).name], (
            f"n_phi={n_phi}: walk says {[g.name for g in computed]}"
        )
        for axis in ("x", "y", "z"):
            assert not mu.is_invariant_under(SubgroupOfO3.SO2(axis)), (
                f"n_phi={n_phi}: no finite point set on S^2 is SO(2)_{axis}-closed"
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
    from orpheus.numerics.invariance import singular_set_under

    for name, q in _cert_rules():
        mu = _measure_from_sphere_quad(q)
        group = mu.symmetry_groups()[0]
        answers = {
            atol: tuple(mu.singular_set_under(group, atol=atol).tolist())
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
    from orpheus.numerics.invariance import singular_set_under

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
        by_group = mu.singular_set_under(SubgroupOfO3.Dnh(2))
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
    from orpheus.numerics.invariance import singular_set_under

    mu = _measure_from_sphere_quad(Quadrature.level_symmetric(sn_order=8))
    coordinate_only = mu.singular_set_under(SubgroupOfO3.Dnh(2))
    full = mu.singular_set_under(mu.symmetry_groups()[0])
    assert len(coordinate_only) == 0
    assert len(full) > 0, "expected O_h diagonal mirrors to fix some ordinate"
    assert set(coordinate_only.tolist()) <= set(full.tolist())

    # product(4, n_phi) has n_phi vertical mirrors, only 2 of them coordinate
    # planes, so Sigma grows with n_phi while the epsilon answer does not.
    sizes = {}
    for n_phi in (4, 8, 16):
        m = _measure_from_sphere_quad(Quadrature.product(n_mu=4, n_phi=n_phi))
        sizes[n_phi] = (
            len(m.singular_set_under(SubgroupOfO3.Dnh(2))),
            len(m.singular_set_under(m.symmetry_groups()[0])),
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
    from orpheus.numerics.invariance import singular_set_under

    lopsided = DiscreteMeasure(
        nodes=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        weights=np.array([1.0, 2.0]),
        support=SPHERE,
    )
    with pytest.raises(ValueError, match="invariant"):
        lopsided.singular_set_under(SubgroupOfO3.OctahedralOh)


@pytest.mark.foundation
def test_orbit_stabilizer_theorem_holds_on_every_certificate() -> None:
    r"""The certificate's defining law: :math:`|Gx|\cdot|\mathrm{Stab}(x)| = |G|`.

    An intrinsic property of the object, not of its use — if the
    permutations were not a genuine group action this would fail, so it
    independently re-checks the bijection requirement (ERR-073) on every
    element rather than only on the duplicate-node witness.
    """
    from orpheus.numerics.symmetry import (_close_group, _realized_ops,)
    from orpheus.numerics.invariance import certificate_under

    for name, q in _cert_rules():
        mu = _measure_from_sphere_quad(q)
        group = mu.symmetry_groups()[0]
        cert = mu.certificate_under(group)
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
    assert not SubgroupOfO3.SO2("z").contains(SubgroupOfO3.Mirror("z"))
    assert SubgroupOfO3.O3.contains(SubgroupOfO3.Mirror("z"))
    # The proper sibling behaves oppositely.
    assert SubgroupOfO3.SO3.contains(SubgroupOfO3.Cn(2))
    assert SubgroupOfO3.SO2("z").contains(SubgroupOfO3.Cn(2))


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
        support=COSINE_INTERVAL,
    )
    symmetric = DiscreteMeasure(
        nodes=np.array([-0.6, -0.2, 0.2, 0.6]),
        weights=np.ones(4),
        support=COSINE_INTERVAL,
    )

    # No rotational group is invariant on a polar marginal, symmetric or
    # not — and crucially these are NOT vacuous "always False" rows: the
    # sigma_x pair below shows the same fixtures ARE distinguishable.
    for rotational in (
        SubgroupOfO3.SO2("z"), SubgroupOfO3.SO2("y"),
        SubgroupOfO3.Cn(3), SubgroupOfO3.SO3,
    ):
        assert not asymmetric.is_invariant_under(rotational), rotational
        assert not symmetric.is_invariant_under(rotational), rotational
    # ⭐ The one rotational group that IS invariant is the marginal's OWN —
    # SO(2)_x acts trivially on the x-axis the marginal embeds along — and
    # it is invariant on the asymmetric fixture too, which is exactly why
    # it cannot be the discriminator: it says nothing about the node set.
    assert asymmetric.is_invariant_under(SubgroupOfO3.SO2("x"))
    assert symmetric.is_invariant_under(SubgroupOfO3.SO2("x"))

    # The reflection is the discriminator, and it separates the fixtures.
    assert not asymmetric.is_invariant_under(SubgroupOfO3.Mirror("x"))
    assert symmetric.is_invariant_under(SubgroupOfO3.Mirror("x"))


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
        SubgroupOfO3.Trivial, SubgroupOfO3.Mirror("z"),
        SubgroupOfO3.SO2("x"), SubgroupOfO3.SO2("y"), SubgroupOfO3.SO2("z"),
        SubgroupOfO3.Dinfh, SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
        SubgroupOfO3.Cn(1), SubgroupOfO3.Cn(2), SubgroupOfO3.Cn(3),
        SubgroupOfO3.Dnh(1), SubgroupOfO3.Dnh(2), SubgroupOfO3.Dnh(4),
    ]
    measures = {
        "asymmetric": DiscreteMeasure(
            nodes=np.array([-0.9, -0.1, 0.3]), weights=np.ones(3),
            support=COSINE_INTERVAL,
        ),
        "symmetric": DiscreteMeasure(
            nodes=np.array([-0.6, -0.2, 0.2, 0.6]), weights=np.ones(4),
            support=COSINE_INTERVAL,
        ),
        "symmetric-nodes-unequal-weights": DiscreteMeasure(
            nodes=np.array([-0.6, 0.6]), weights=np.array([1.0, 2.0]),
            support=COSINE_INTERVAL,
        ),
    }

    violations: list[str] = []
    checked = 0
    for label, mu in measures.items():
        invariant = {repr(t): mu.is_invariant_under(t) for t in tags}
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
    r"""The retired ``_contains`` (R1 of #434) ended in a bare ``return
    False``, so a tag it had no arm for got a wrong-but-SILENT answer rather
    than an error; the realization has no fall-through — this row pins that
    a mirror against every CONTINUOUS outer group gets a real answer.

    `[M]` That is not hypothetical: a prototype mirror with only the
    realization and invariance arms extended answered
    ``O3.contains(Mirror('x')) = False`` and
    ``Dinfh.contains(Mirror('x')) = False`` — both flatly wrong, and
    ``O(3) ⊉`` anything breaks the soundness precondition
    :func:`~orpheus.numerics.invariance.symmetry_groups` states
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
    for rot_axis in ("x", "y", "z"):
        assert not SubgroupOfO3.SO2(rot_axis).contains(sigma)


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
    for rot_axis in ("x", "y", "z"):
        assert not sigma.contains(SubgroupOfO3.SO2(rot_axis))


@pytest.mark.foundation
def test_the_rotation_axis_is_validated_at_CONSTRUCTION() -> None:
    """There is no unnamed axial rotation group — which is what the retired
    bare ``SO2`` entry pretended there was (tracker 2.4, 2026-09-01)."""
    for axis in ("x", "y", "z"):
        g = SubgroupOfO3.SO2(axis)
        assert g.name == f"SO2_{axis}"
        assert repr(g) == f"SubgroupOfO3.SO2({axis!r})"
    with pytest.raises(ValueError, match="axis in x/y/z"):
        SubgroupOfO3.SO2("w")
    with pytest.raises(ValueError, match="unnamed axial rotation"):
        SubgroupOfO3.SO2("")


@pytest.mark.foundation
def test_so2_axes_are_three_incomparable_groups() -> None:
    r"""Three axes, three groups: :math:`SO(2)_a \cap SO(2)_b = \{e\}` for
    :math:`a \ne b`, so neither contains the other, and they are not equal
    — the property that keeps ``S^2/O2_x`` and ``S^2/O2_z`` two different
    orbit spaces in the catalogue."""
    axes = ("x", "y", "z")
    for a in axes:
        for b in axes:
            g, h = SubgroupOfO3.SO2(a), SubgroupOfO3.SO2(b)
            assert g.contains(h) is (a == b)
            assert (g == h) is (a == b)
            assert (hash(g) == hash(h)) is (a == b)
        # The only finite subgroup of any of them is the trivial one.
        assert SubgroupOfO3.SO2(a).contains(SubgroupOfO3.Trivial)
        for n in (1, 2, 3, 4):
            assert not SubgroupOfO3.SO2(a).contains(SubgroupOfO3.Dnh(n))


@pytest.mark.foundation
def test_rotation_axis_is_the_continuous_dual_of_mirror_axis() -> None:
    """``rotation_axis`` answers the axial family alone — exactly as
    ``mirror_axis`` answers the reflection family alone — so an orbit-space
    derivation can read the invariant coordinate off the group and nothing
    that merely CONTAINS axial rotations answers by accident."""
    for axis, index in (("x", 0), ("y", 1), ("z", 2)):
        assert SubgroupOfO3.SO2(axis).rotation_axis == index
        assert SubgroupOfO3.SO2(axis).mirror_axis is None
        assert SubgroupOfO3.Mirror(axis).rotation_axis is None
    for contains_rotations in (
        SubgroupOfO3.Dinfh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
        SubgroupOfO3.Cn(4), SubgroupOfO3.Dnh(2), SubgroupOfO3.Trivial,
    ):
        assert contains_rotations.rotation_axis is None


@pytest.mark.foundation
def test_candidate_groups_offer_all_three_axial_rotation_groups() -> None:
    r"""The walk is offered every axis of BOTH axial families, so a polar
    marginal along x reports :math:`O(2)_x` — its full stabiliser — rather
    than reading as carrying no continuous symmetry (what the retired
    z-only ``SO2`` amounted to) or stopping at the rotations (what the walk
    reported before #432). Walk and brute force must agree on it, as on
    everything.

    `[M]` 2026-09-03 (R2 of #434) the slab's maximal elements are
    ``{O2_x, D_2h}``: :math:`SO(2)_x`, :math:`\sigma_y` and :math:`\sigma_z`
    all sit inside :math:`O(2)_x`; :math:`\sigma_x` — which FLIPS the axis —
    sits in neither, and since the candidate set is derived from the orbit
    BARYCENTRES ``(mu, 0, 0)`` (two z-azimuths) :math:`D_{2h}` is offered,
    is invariant (each of its elements normalises :math:`O(2)_x` and acts on
    :math:`\mu` by :math:`\pm`), and ABSORBS :math:`\sigma_x`.  (Until R2 the
    candidate set was read off the STORED width — a 1-D rule offered no
    :math:`C_n`/:math:`D_{nh}` at all — and the answer was
    ``{O2_x, sigma_x}``.)  The group the two generate is :math:`D_{\infty h}`
    about :math:`x`, which the lattice cannot spell (``Dinfh`` is about z),
    so the two incomparable maxima are the honest answer."""
    from orpheus.numerics.invariance import candidate_groups

    q = Quadrature.gauss_legendre(n_ordinates=8)
    offered = {g.name for g in candidate_groups(q.measure)}
    assert {"SO2_x", "SO2_y", "SO2_z", "O2_x", "O2_y", "O2_z"} <= offered
    walk = {g.name for g in q.measure.symmetry_groups()}
    brute = {
        g.name for g in q.measure.symmetry_groups(method="bruteforce")
    }
    assert walk == brute
    assert walk == {"O2_x", "D_2h"}
    # The groups the maximum absorbs: the rotations and the two mirrors
    # THROUGH the axis are invariant AND inside O2_x, hence not maximal.
    for absorbed in ("SO2_x", "sigma_y", "sigma_z"):
        g = next(c for c in candidate_groups(q.measure) if c.name == absorbed)
        assert q.measure.is_invariant_under(g) and SubgroupOfO3.O2("x").contains(g)
    assert not SubgroupOfO3.O2("x").contains(SubgroupOfO3.Mirror("x"))


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
        one_d = DiscreteMeasure(nodes=mu, weights=weights, support=COSINE_INTERVAL)
        embedded = DiscreteMeasure(
            nodes=np.column_stack([mu, np.zeros(4), np.zeros(4)]),
            weights=weights,
            support=SPHERE,
        )
        g = SubgroupOfO3.Mirror(axis)
        return one_d.is_invariant_under(g), embedded.is_invariant_under(g)

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
    assert not measure.is_invariant_under(SubgroupOfO3.Mirror("x"))
    assert measure.is_invariant_under(SubgroupOfO3.Mirror("y"))
    assert measure.is_invariant_under(SubgroupOfO3.Mirror("z"))


# ============================================================================
# The axial STABILISER O(2)_a = C_∞v — #432, 2026-09-02
# ============================================================================
#
# An orbit space is named by the LARGEST group with its orbits. The
# constant-μ circles on S² are the orbits of SO(2)_a AND of its stabiliser
# O(2)_a = {g ∈ O(3) : g ê_a = ê_a} — a vertical mirror carries each circle
# onto itself — so the two have one invariant ring, one derivation and one
# orbit space, and the entry records the bigger group. This section gates the
# new member: its lattice edges, its exact invariance criterion, and the
# two-component generic-image set the isotypic probe reads.
#
# Claim layer: **term-level (L0)**, closed form. Every reference below is a
# group-theoretic fact written from the definition (`vv-principles` §pillars,
# and §8 of the test-architect lessons: a pure-math primitive has no MMS row
# and no semi-analytical row — every row is closed form or exact integer
# arithmetic).


#: The groups the lattice can spell today, in one list so the order-relation
#: laws below run over an ENUMERATED denominator rather than a sample
#: (`plan-authoring` §2). 23 members: the named six, three mirrors, three
#: SO(2), three O(2), C_1..C_4, D_1h..D_4h.
_SPELLABLE = (
    [
        SubgroupOfO3.Trivial, SubgroupOfO3.Dinfh, SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
    ]
    + [SubgroupOfO3.Mirror(a) for a in ("x", "y", "z")]
    + [SubgroupOfO3.SO2(a) for a in ("x", "y", "z")]
    + [SubgroupOfO3.O2(a) for a in ("x", "y", "z")]
    + [SubgroupOfO3.Cn(n) for n in (2, 3, 4)]  # C_1 IS Trivial (R1 of #434)
    + [SubgroupOfO3.Dnh(n) for n in (1, 2, 3, 4)]
)


@pytest.mark.foundation
def test_the_O2_axis_is_validated_at_CONSTRUCTION() -> None:
    r"""There is no unnamed axial GROUP, for the reason there is no unnamed
    axial ROTATION group (tracker 2.4) and no unnamed reflection (2026-08-02).

    The name and the ``repr`` both carry the axis, and both are load-bearing
    rather than cosmetic: :attr:`SubgroupOfO3.name` keys the orbit-space
    catalogue (``S^2/O2_x``), and the VALUE — whose dataclass ``__eq__`` /
    ``__hash__`` read the tag — keys the realization cache and the lattice
    walk's ``visited`` set (R1 of #434; until then ``repr`` keyed both), so
    an axis-blind spelling would silently merge three groups into one entry.
    """
    for axis in ("x", "y", "z"):
        g = SubgroupOfO3.O2(axis)
        assert g.name == f"O2_{axis}"
        assert repr(g) == f"SubgroupOfO3.O2({axis!r})"
    # three names, three reprs — the collapse the cache would suffer.
    assert len({SubgroupOfO3.O2(a).name for a in ("x", "y", "z")}) == 3
    assert len({repr(SubgroupOfO3.O2(a)) for a in ("x", "y", "z")}) == 3

    with pytest.raises(ValueError, match="axis in x/y/z"):
        SubgroupOfO3.O2("w")
    with pytest.raises(ValueError, match="unnamed axial group"):
        SubgroupOfO3.O2("")


@pytest.mark.foundation
def test_rotation_axis_answers_for_the_STABILISER_too_and_mirror_axis_does_not() -> None:
    r"""``rotation_axis`` is the *axial-family* accessor, not the
    *rotation-family* one — it answers for :math:`O(2)_a` as well as
    :math:`SO(2)_a`, because both fix :math:`\hat e_a` pointwise and both
    name the axis whose polar interval their orbit space IS.

    That is exactly what ``manifold._sphere_mod_o2`` reads to pick the
    surviving invariant AND, since R4 of #434 (2026-09-03), to call
    ``_coordinate_chart`` for the entry's chart/lift pair — the lift is a
    FIELD the builder populates, no longer a branch that read this property
    at lift time (until 2026-09-02 a separate ``_polar_axis_of`` read it to
    pick an embedding axis; #429 tracker 2.2b retired that into the entry's
    lift; R4 retired the lift's own branch). A property that answered
    ``None`` for :math:`O(2)_a` would send the builder through the ``None``
    branch.

    The companion of the existing SO(2)-only row
    :func:`test_rotation_axis_is_the_continuous_dual_of_mirror_axis`; the
    NEGATIVE half is shared with it — nothing that merely *contains* axial
    rotations may answer.
    """
    for axis, index in (("x", 0), ("y", 1), ("z", 2)):
        assert SubgroupOfO3.O2(axis).rotation_axis == index
        assert SubgroupOfO3.O2(axis).mirror_axis is None
    # NEGATIVE: D_∞h contains every rotation about z and answers None,
    # because sigma_h FLIPS the axis — the discriminating case, since a
    # naive "does it contain the rotations?" reading would answer 2.
    assert SubgroupOfO3.Dinfh.contains(SubgroupOfO3.SO2("z"))
    assert SubgroupOfO3.Dinfh.rotation_axis is None
    for other in (
        SubgroupOfO3.SO3, SubgroupOfO3.O3, SubgroupOfO3.Cn(4),
        SubgroupOfO3.Dnh(2), SubgroupOfO3.Trivial,
        SubgroupOfO3.Mirror("x"), SubgroupOfO3.OctahedralOh,
    ):
        assert other.rotation_axis is None, other.name


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_the_stabiliser_sits_strictly_above_its_rotation_half(axis: str) -> None:
    r"""The identity component: :math:`SO(2)_a \subsetneq O(2)_a`, and the
    containment does NOT run the other way — the vertical mirrors are
    improper, so they lie in no rotation group.

    Against the continuous named entries, the two families answer
    DIFFERENTLY exactly where properness bites: both sit inside
    :math:`O(3)`, only the proper half sits inside :math:`SO(3)`. That
    asymmetry is the whole content of "stabiliser, not rotation group", and
    it is the row a mutation dropping the determinant test from
    ``IdentityComponent.contains_element`` (the proper/improper distinction,
    one body since R1 of #434) must redden.
    """
    rot, stab = SubgroupOfO3.SO2(axis), SubgroupOfO3.O2(axis)

    assert stab.contains(rot)
    assert not rot.contains(stab)
    assert stab != rot
    # ⚠ ``hash`` does NOT discriminate them and must not be asserted to:
    # a frozen dataclass's generated ``__hash__`` hashes the FIELD TUPLE and
    # not the class, so `[M]` 2026-09-02 ``hash(SO2(a)) == hash(O2(a)) ==
    # hash(Mirror(a))`` on every axis (and ``hash(Cn(n)) == hash(Dnh(n))``).
    # Equality still discriminates — the generated ``__eq__`` requires
    # ``other.__class__ is self.__class__`` — so dicts and sets are correct
    # and this is a collision, not a bug. The property that matters is the
    # one asserted here.
    assert len({stab, rot, SubgroupOfO3.Mirror(axis)}) == 3

    assert SubgroupOfO3.O3.contains(stab)
    assert SubgroupOfO3.O3.contains(rot)
    # the discriminating pair: O(2)_a carries reflections, SO(2)_a does not.
    assert not SubgroupOfO3.SO3.contains(stab)
    assert SubgroupOfO3.SO3.contains(rot)

    # D_∞h = O(2)_z × {e, σ_h} is realized about z, so it contains BOTH
    # axial families about z and NEITHER about any other axis.
    assert SubgroupOfO3.Dinfh.contains(stab) is (axis == "z")
    assert SubgroupOfO3.Dinfh.contains(rot) is (axis == "z")
    # ... and no axial group contains D_∞h, SO(3) or O(3): none fixes an axis.
    for continuous in (SubgroupOfO3.Dinfh, SubgroupOfO3.SO3, SubgroupOfO3.O3):
        assert not stab.contains(continuous)
        assert not rot.contains(continuous)


@pytest.mark.foundation
def test_O2_axes_are_three_incomparable_groups() -> None:
    r"""Three axes, three groups — the companion of
    :func:`test_so2_axes_are_three_incomparable_groups`.

    :math:`O(2)_a \cap O(2)_b = \{e, \sigma_c\}` for :math:`a \ne b`, which
    is neither of them, so the two are incomparable; and equality and hash
    agree with containment. That is the property keeping ``S^2/O2_x`` and
    ``S^2/O2_z`` two different orbit spaces in the catalogue, and the one a
    repr or a ``__eq__`` that dropped the axis would destroy.
    """
    axes = ("x", "y", "z")
    for a in axes:
        for b in axes:
            g, h = SubgroupOfO3.O2(a), SubgroupOfO3.O2(b)
            assert g.contains(h) is (a == b)
            assert (g == h) is (a == b)
            assert (hash(g) == hash(h)) is (a == b)
            # cross-family: O(2)_a ⊇ SO(2)_b iff a == b, never conversely.
            assert g.contains(SubgroupOfO3.SO2(b)) is (a == b)
            assert SubgroupOfO3.SO2(b).contains(g) is False


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_exactly_the_mirrors_THROUGH_the_axis_lie_in_its_stabiliser(
    axis: str,
) -> None:
    r"""``O2(a) ⊇ Mirror(b)`` **iff** :math:`b \ne a` — computed from
    :math:`\sigma_b`'s own realization, not tabulated.

    :math:`\sigma_b` has normal :math:`\hat e_b`, so its plane CONTAINS
    :math:`\hat e_a` for :math:`b \ne a` (a vertical mirror, which fixes the
    axis pointwise) and is PERPENDICULAR to it for :math:`b = a` (the
    horizontal mirror :math:`\sigma_h`, which sends :math:`\hat e_a \to
    -\hat e_a`). Two of three planes in, one out — the fact that makes the
    Legendre basis admissible on a :math:`\sigma_b`-fold and refused on a
    :math:`\sigma_a`-fold (``tests/numerics/test_frame.py``'s nine-pairing
    table is the same statement one tier up).

    The rotation half admits NO mirror on any axis (``det = -1``), which is
    the leg that keeps this from reading as "any reflection is fine".
    """
    stab, rot = SubgroupOfO3.O2(axis), SubgroupOfO3.SO2(axis)
    for plane in ("x", "y", "z"):
        assert stab.contains(SubgroupOfO3.Mirror(plane)) is (plane != axis), (
            f"O2({axis}) vs sigma_{plane}"
        )
        assert not rot.contains(SubgroupOfO3.Mirror(plane))
    # exactly two of the three coordinate mirrors, never three, never one.
    assert sum(
        stab.contains(SubgroupOfO3.Mirror(p)) for p in ("x", "y", "z")
    ) == 2


@pytest.mark.foundation
def test_the_finite_relations_of_the_stabiliser_are_COMPUTED_from_realizations() -> None:
    r"""Every :math:`(\text{axial}, \text{finite})` relation follows from one
    question — *does every element of the finite group fix* :math:`\hat e_a`? —
    so the answers are n-dependent and axis-dependent in ways no enum-to-enum
    table can spell.

    `[M]` 2026-09-02, over the four orders and three axes below:

    * :math:`C_n \subseteq O(2)_z` for every :math:`n` (realized about z),
      and inside another axis's stabiliser only for :math:`n = 1`;
    * :math:`D_{1h} = \{e, \sigma_z, \sigma_y, C_2(x)\}` — order 4 — is
      inside :math:`O(2)_x` **and no other axis's**, because every one of
      its four elements fixes :math:`\hat e_x`;
    * :math:`D_{nh} \not\subseteq O(2)_a` for :math:`n \ge 2` on every axis
      (its in-plane :math:`C_2` axes flip whichever axis is not theirs);
    * :math:`O_h` and :math:`I_h` are inside no axial group (they move every
      axis), and no axial group is inside them (a finite group cannot
      contain a continuous one).

    The :math:`D_{1h} \subseteq O(2)_x`-only row is the sharp one: read as
    "a dihedral group is never inside an axial group" it would be a false
    universal, and read as "D_1h is inside every O(2)_a" it would be false
    on two of three axes. Only the computed relation gets both right.
    """
    for n in (1, 2, 3, 4, 6):
        assert SubgroupOfO3.O2("z").contains(SubgroupOfO3.Cn(n)), n
        for other in ("x", "y"):
            assert SubgroupOfO3.O2(other).contains(SubgroupOfO3.Cn(n)) is (n == 1)

    for axis in ("x", "y", "z"):
        assert SubgroupOfO3.O2(axis).contains(SubgroupOfO3.Dnh(1)) is (axis == "x")
        for n in (2, 3, 4, 6):
            assert not SubgroupOfO3.O2(axis).contains(SubgroupOfO3.Dnh(n)), (axis, n)
        for polyhedral in (
            SubgroupOfO3.OctahedralOh, SubgroupOfO3.IcosahedralIh,
        ):
            assert not SubgroupOfO3.O2(axis).contains(polyhedral)
            assert not polyhedral.contains(SubgroupOfO3.O2(axis))
        # ... and the trivial group is inside every one of them — under ONE
        # spelling: since R1 of #434 `Cn(1)` IS `Trivial` (the normalisation
        # is the fact; asserting containment of both would be one assertion
        # written twice — see `test_cn_in_so2`'s history note).
        assert SubgroupOfO3.Cn(1) == SubgroupOfO3.Trivial
        assert SubgroupOfO3.O2(axis).contains(SubgroupOfO3.Trivial)


@pytest.mark.foundation
def test_the_whole_spellable_lattice_obeys_the_order_relation_laws() -> None:
    r"""``vv-principles`` #15's compatibility discipline applied to the
    lattice itself, over EVERY group the module can spell — continuous
    members included.

    The existing :func:`test_computed_containment_obeys_the_order_relation_laws`
    runs over nine FINITE groups, so the axial families (which have no finite
    realization to compare) are outside it entirely. This row widens the
    denominator to all 22 spellable groups — `[M]` 2026-09-03: **484** ordered
    pairs, **108** strict edges, of which **20 touch an** :math:`O(2)_a`
    (23 members / 529 / 131 / 23 until R1 of #434 merged ``C_1`` into
    ``Trivial``: every edge the duplicate entered was an edge the merged value
    already entered) — and checks the three laws an order relation must
    satisfy:

    * reflexivity on every member;
    * antisymmetry, stated as *mutual containment implies the same group* —
      by ORDER for a finite pair and by ``==`` when either side is
      continuous. Until R1 of #434 ``C_1`` was a second value for
      :math:`\{e\}` that compared UNEQUAL to ``Trivial`` while containing
      it, which is why the finite arm compares orders rather than identity;
      the arm survives because two distinct finite groups can still share an
      order, and the merge is asserted by ``TestR1TheTrivialGroupHasExactlyOneSpelling``;
    * transitivity, over all :math:`22^3 = 10648` triples.

    Lagrange (a subgroup's order divides the group's) is asserted only on the
    finite × finite pairs, where an order exists.

    The containment MATRIX is computed once and the laws are then read off
    it — 484 ``contains`` calls rather than 10648, which is what makes the
    triple loop affordable (`[M]` 0.43 s for the matrix before R1's memo).
    """
    from orpheus.numerics.symmetry import _group_elements

    n = len(_SPELLABLE)
    assert n == 22, f"the spellable set moved: {n} members"
    assert len(set(_SPELLABLE)) == n, "the spellable set lists one group twice"
    below = [[a.contains(b) for b in _SPELLABLE] for a in _SPELLABLE]
    elements = [_group_elements(g._tag) for g in _SPELLABLE]
    order = [None if e is None else len(e) for e in elements]

    for i, g in enumerate(_SPELLABLE):
        assert below[i][i], f"not reflexive at {g.name}"

    for i in range(n):
        for j in range(n):
            if i == j or not (below[i][j] and below[j][i]):
                continue
            if order[i] is None or order[j] is None:
                assert _SPELLABLE[i] == _SPELLABLE[j], (
                    f"{_SPELLABLE[i].name} and {_SPELLABLE[j].name} contain "
                    f"each other but are not equal, and one is continuous"
                )
            else:
                assert order[i] == order[j], (
                    f"{_SPELLABLE[i].name} and {_SPELLABLE[j].name} contain "
                    f"each other but have orders {order[i]} != {order[j]}"
                )

    for i in range(n):
        for j in range(n):
            outer_order, inner_order = order[i], order[j]
            if not below[i][j] or outer_order is None or inner_order is None:
                continue
            assert outer_order % inner_order == 0, (
                f"{_SPELLABLE[i].name} contains {_SPELLABLE[j].name} but "
                f"{inner_order} does not divide {outer_order} — impossible"
            )

    broken = [
        (_SPELLABLE[i].name, _SPELLABLE[j].name, _SPELLABLE[k].name)
        for i in range(n) for j in range(n) for k in range(n)
        if below[i][j] and below[j][k] and not below[i][k]
    ]
    assert not broken, f"{len(broken)} transitivity violations: {broken[:5]}"

    # NON-VACUITY (`vv-principles` #20): the widening must actually reach the
    # new member, or the row is the old gate wearing a bigger list.
    axial = {i for i, g in enumerate(_SPELLABLE) if g.name.startswith("O2_")}
    touching = sum(
        1 for i in range(n) for j in range(n)
        if i != j and below[i][j] and (i in axial or j in axial)
    )
    assert touching == 20, f"O(2)-touching strict edges: {touching}"
    strict = sum(1 for i in range(n) for j in range(n) if i != j and below[i][j])
    assert strict == 108, f"strict edges: {strict}"


# ---------------------------------------------------------------------------
# Exact invariance for O(2)_a — decided, never sampled (ERR-072)
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_O2_invariance_is_axis_support_and_needs_NO_mu_reflection(
    axis: str,
) -> None:
    r"""A finite point set is :math:`O(2)_a`-closed **iff** every node lies
    ON axis :math:`a` — the same exact criterion as :math:`SO(2)_a`, because
    the rotation half already forces axis support and a point on the axis is
    fixed by every vertical mirror. Decided, never sampled (ERR-072).

    ⭐ The ASYMMETRIC leg is the discriminating one and it is the reason this
    row is not a restatement of the :math:`\sigma`-mirror gates: a
    :math:`\mu`-set with no :math:`\mu \to -\mu` symmetry is
    :math:`O(2)_a`-invariant (``True``) while :math:`\sigma_a` — the mirror
    that FLIPS the axis, the slab's own residual — reads ``False`` on the
    same data. `[M]` 2026-09-02 on ``mu = [-0.9, -0.3, 0.3, 0.7]`` declared
    on ``S^2/O2_a``: :math:`O(2)_a` ``True``, :math:`\sigma_a` ``False``,
    both other axes' stabilisers ``False``.

    A criterion that had (wrongly) required the reflection would read
    ``False`` on the asymmetric leg; one that had dropped the axis would read
    ``True`` on all three axes. The symmetric leg is the CONTROL that keeps
    the agreement from being the trivial "everything is True".
    """
    from orpheus.numerics.manifold import SPHERE

    entry = SPHERE.quotient(SubgroupOfO3.O2(axis))
    weights = np.full(4, 0.5)
    asymmetric = np.array([-0.9, -0.3, 0.3, 0.7])   # violates mu -> -mu
    symmetric = np.array([-0.7, -0.3, 0.3, 0.7])

    for mu, mirror_holds in ((asymmetric, False), (symmetric, True)):
        measure = DiscreteMeasure(nodes=mu, weights=weights, support=entry)
        for other in ("x", "y", "z"):
            assert measure.is_invariant_under(SubgroupOfO3.O2(other)) is (other == axis)
            assert measure.is_invariant_under(SubgroupOfO3.SO2(other)) is (other == axis)
        # the axis-flipping mirror is a DIFFERENT question, and it separates
        # the two legs — which is what makes the O(2) answer non-vacuous.
        assert measure.is_invariant_under(SubgroupOfO3.Mirror(axis)) is mirror_holds


@pytest.mark.foundation
def test_no_shipped_sphere_cubature_is_O2_invariant_on_ANY_axis() -> None:
    r"""The intended consequence of deciding the continuous groups exactly:
    a rule with genuinely off-axis nodes has infinite orbits and cannot be
    closed under any :math:`O(2)_a`.

    `[M]` 2026-09-02, **4 rules × 3 axes = 12 verdicts, all ``False``**:
    ``product(4, 8)``, ``level_symmetric(4)``, ``lebedev(11)``,
    ``folded_product(4, 8)``. The 1-D polar marginal — the one measure that
    DOES pass, and only for its own axis — is the positive control, so this
    is a two-sided statement rather than a blanket refusal.

    The sibling :func:`test_no_discrete_cubature_is_so2_invariant` makes the
    same claim for the rotation half; both are needed, because a criterion
    keyed on the wrong half would satisfy one and not the other.
    """
    rules = {
        "product(4,8)": Quadrature.product(n_mu=4, n_phi=8).measure,
        "level_symmetric(4)": Quadrature.level_symmetric(sn_order=4).measure,
        "lebedev(11)": Quadrature.lebedev(order=11).measure,
        "folded_product(4,8)": Quadrature.folded_product(n_mu=4, n_phi=8).measure,
    }
    verdicts = 0
    for name, measure in rules.items():
        for axis in ("x", "y", "z"):
            verdicts += 1
            assert not measure.is_invariant_under(SubgroupOfO3.O2(axis)), (
                f"{name} certified O(2)_{axis}-invariant, but a finite point "
                f"set with off-axis nodes has infinite O(2) orbits"
            )
    assert verdicts == 12

    # POSITIVE CONTROL: the polar marginal passes, for its own axis alone.
    from orpheus.numerics.quadrature import gauss_legendre_on_polar_orbit

    for axis in ("x", "y", "z"):
        marginal = gauss_legendre_on_polar_orbit(8, axis=axis)
        assert marginal.is_invariant_under(SubgroupOfO3.O2(axis))


@pytest.mark.foundation
def test_the_downward_closure_law_actually_EXERCISES_the_new_member() -> None:
    r"""``vv-principles`` #15 / #20 — the compatibility law's DENOMINATOR.

    :func:`test_invariance_is_downward_closed` runs over
    ``candidate_groups``, which since #432 offers the three
    :math:`O(2)_a`. A widened candidate list makes that gate look stronger
    without proving it reaches the new member: a pair
    :math:`(\text{outer}, \text{inner})` only enters the law when
    ``outer.contains(inner)`` holds, so an O(2) row that is inside nothing
    and contains nothing would be checked ZERO times and the law would still
    report "no violations".

    `[M]` 2026-09-02, over the SEVEN rules :func:`_walk_rules` supplies
    (``product(4, n_phi)`` for ``n_phi`` in 3/4/8, ``level_symmetric`` 4/8,
    ``lebedev`` 5/11): **919** ordered pairs enter the law, of which **145**
    have an :math:`O(2)_a` on one side — and BOTH directions are populated
    (**117** with the stabiliser as the outer group, **28** as the inner).
    Zero violations on the exercised subset alone.  (1090 / 166 / 138 / 28
    until R1 of #434 the same day: ``C_1`` was a second candidate for
    :math:`\{e\}`, and every pair it entered — it lies inside everything —
    was a pair the merged ``Trivial`` already entered.)
    """
    from orpheus.numerics.invariance import candidate_groups

    checked = 0
    touching_o2 = 0
    o2_outer = 0
    o2_inner = 0
    violations: list[str] = []
    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        cands = candidate_groups(mu)
        invariant = {repr(g): mu.is_invariant_under(g) for g in cands}
        for outer in cands:
            for inner in cands:
                if outer == inner or not outer.contains(inner):
                    continue
                checked += 1
                out_is, in_is = (
                    outer.name.startswith("O2_"), inner.name.startswith("O2_"),
                )
                touching_o2 += int(out_is or in_is)
                o2_outer += int(out_is)
                o2_inner += int(in_is)
                if invariant[repr(outer)] and not invariant[repr(inner)]:
                    violations.append(f"{name}: {outer.name} > {inner.name}")
    assert not violations, violations[:10]
    assert checked == 919, f"pairs entering the law: {checked}"
    assert touching_o2 == 145, f"pairs touching an O(2): {touching_o2}"
    assert o2_outer == 117 and o2_inner == 28, (o2_outer, o2_inner)


# ---------------------------------------------------------------------------
# generic_images — the two components the isotypic probe samples
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_generic_images_of_the_stabiliser_carry_BOTH_components(
    axis: str,
) -> None:
    r"""``generic_images(O2(a))`` is :math:`SO(2)_a`'s six incommensurate
    rotations **plus** each of them composed with one vertical mirror — a
    generic element of each of the group's two connected components.

    The determinant is the component label and it is exact: `[M]` 2026-09-02
    the twelve reconstructed linear parts read :math:`+1` six times and
    :math:`-1` six times, against :math:`SO(2)_a`'s six :math:`+1`. A
    mutation that drops the mirrored half moves the count 12 → 6 and the
    determinant multiset; one that composes with the HORIZONTAL mirror
    :math:`\sigma_a` instead keeps both, and is caught by the next row.

    ⚠ The count is not arbitrary — six is
    ``len(symmetry._INCOMMENSURATE_ANGLES)``, read here rather than written,
    so a future angle added to the sample does not red this row for the
    wrong reason.
    """
    from orpheus.numerics.symmetry import _INCOMMENSURATE_ANGLES

    points = _generic_directions()
    m = len(_INCOMMENSURATE_ANGLES)
    rotations = SubgroupOfO3.SO2(axis).generic_images(points)
    both = SubgroupOfO3.O2(axis).generic_images(points)
    assert len(rotations) == m
    assert len(both) == 2 * m

    dets = [round(float(np.linalg.det(_linear_part(points, im))), 9) for im in both]
    assert dets.count(1.0) == m and dets.count(-1.0) == m, dets
    assert all(
        round(float(np.linalg.det(_linear_part(points, im))), 9) == 1.0
        for im in rotations
    )
    # the proper half of O(2)_a IS SO(2)_a's sample, element for element.
    for proper, rotation in zip(both[:m], rotations):
        assert np.array_equal(proper, rotation)


@pytest.mark.foundation
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_every_generic_image_of_the_stabiliser_FIXES_THE_AXIS(axis: str) -> None:
    r"""The DEFINING property, as a gate: :math:`O(2)_a = \{g : g\hat e_a =
    \hat e_a\}`, so every sampled element — rotation and mirrored alike —
    must leave :math:`\hat e_a` exactly where it is, and must be an isometry.

    `[M]` 2026-09-02: :math:`\max|g\hat e_a - \hat e_a| = 0.000e+00` over
    all twelve images on all three axes (the coordinate mirror is a signed
    diagonal, hence bit-exact, and the rotation fixes the axis by
    construction), and the Gram matrix of nine generic directions moves by
    at most ``5.6e-16``.

    ⭐ The NEGATIVE control is the mirror that is NOT in the group: the
    horizontal :math:`\sigma_a` sends :math:`\hat e_a \to -\hat e_a`, so a
    ``generic_images`` composing with it instead of a vertical mirror would
    read :math:`|g\hat e_a - \hat e_a| = 2`. Asserted here, because without
    it the row would pass for a sample that fixed the axis by accident —
    and ``2.0`` is what the mutation actually produces.
    """
    from orpheus.geometry.transformation import RigidMotion
    from orpheus.numerics.symmetry import _axis_vector

    e_a = _axis_vector(axis)
    images = SubgroupOfO3.O2(axis).generic_images(e_a.reshape(1, 3))
    assert len(images) == 12
    assert max(float(np.max(np.abs(im[0] - e_a))) for im in images) == 0.0

    points = _generic_directions()
    gram = points @ points.T
    for image in SubgroupOfO3.O2(axis).generic_images(points):
        np.testing.assert_allclose(image @ image.T, gram, atol=1e-14)

    # NEGATIVE control — the horizontal mirror is outside O(2)_a and moves
    # the axis by the full diameter.
    sigma_h = RigidMotion.reflection(normal=e_a)
    assert float(np.max(np.abs(sigma_h.linear @ e_a - e_a))) == 2.0
    assert not SubgroupOfO3.O2(axis).contains(SubgroupOfO3.Mirror(axis))


@pytest.mark.foundation
def test_the_MIRRORED_half_is_inert_for_the_incommensurate_probe_and_that_is_a_THEOREM() -> None:
    r"""⛔ **A declared blindness, measured** — and the refutation of the
    obvious control for the previous two rows.

    The natural companion to "``generic_images(O2)`` carries both
    components" is *"a* :math:`\sigma_v`-odd function is constant across the
    rotation half and not across the mirrored half"*. **No such function
    exists.** If :math:`f` is constant on every :math:`SO(2)_a` orbit and
    :math:`\sigma_v` carries each constant-:math:`\mu` circle onto itself,
    then :math:`f(\sigma_v x) = f(x)` for every :math:`x` — which is exactly
    the theorem ``manifold._sphere_mod_o2`` records,
    :math:`\mathbb{R}[x]^{SO(2)_a} = \mathbb{R}[x]^{O(2)_a}`, and the reason
    ONE derivation serves both groups.

    `[M]` 2026-09-02, over **3 axes × L = 1..6 = 18 rows**: the descending
    mask computed from ``generic_images(O2(a))`` is ``array_equal`` to the
    one computed from ``generic_images(SO2(a))``. So dropping the mirrored
    half is production-INERT at the shipped incommensurate angles, and no
    ``descending_slots`` row can catch it.

    ⟹ Where the mirrored half IS observable is the DEGENERATE-sample
    regime, and it is gated elsewhere:
    ``test_manifold.py::TestDescendingSlots::test_a_right_angle_sample_falsely_admits_slots_from_L_four``
    measures right angles generating :math:`C_{4v}` rather than
    :math:`C_4` (6 rather than 7 live slots at :math:`L = 4`; 8 rather than
    10 at :math:`L = 5`). This row exists so the inertness is written down
    rather than discovered by a battery arm that reddens nothing.
    """
    from orpheus.numerics.basis.spherical_harmonic_basis import (
        SphericalHarmonicBasis,
    )

    points = _generic_directions()

    def mask(basis, group) -> np.ndarray:
        reference = np.asarray(basis.evaluate(points))
        out = np.ones(reference.shape[1:], dtype=bool)
        for image in group.generic_images(points):
            out &= np.all(
                np.isclose(basis.evaluate(image), reference, rtol=0.0, atol=1e-12),
                axis=0,
            )
        return out

    rows = 0
    for axis in ("x", "y", "z"):
        for L in (1, 2, 3, 4, 5, 6):
            basis = SphericalHarmonicBasis(L=L)
            both = mask(basis, SubgroupOfO3.O2(axis))
            proper = mask(basis, SubgroupOfO3.SO2(axis))
            assert np.array_equal(both, proper), (
                f"axis {axis}, L={L}: the mirrored half changed the mask — "
                f"which the invariant-ring theorem forbids"
            )
            rows += 1
    assert rows == 18

    # NON-VACUITY: the masks are not "everything" — about x they select one
    # real slot per degree, so the equality above is a real coincidence of
    # two selections rather than of two universal answers.
    basis = SphericalHarmonicBasis(L=4)
    live = mask(basis, SubgroupOfO3.O2("x")) & basis.live_slot_mask
    assert int(live.sum()) == 5 and int(basis.live_slot_mask.sum()) == 25


def _generic_directions(seed: int = 20260902, n: int = 9) -> np.ndarray:
    """Deterministic generic unit directions — no component zero, no pair
    related by a coordinate symmetry."""
    raw = np.random.default_rng(seed).normal(size=(n, 3))
    return raw / np.linalg.norm(raw, axis=1)[:, None]


def _linear_part(points: np.ndarray, image: np.ndarray) -> np.ndarray:
    """Recover ``g`` from ``points @ g.T = image`` — the test reads the
    MATRIX back out of the action rather than importing the constructor the
    production code used."""
    solved, *_ = np.linalg.lstsq(points, image, rcond=None)
    return solved.T


# ---------------------------------------------------------------------------
# The orbit STABILISER — the group an orbit space is named by
# ---------------------------------------------------------------------------


class TestOrbitStabiliser:
    r"""The defining laws of :attr:`SubgroupOfO3.orbit_stabiliser` — the
    largest subgroup of :math:`O(3)` with this group's orbits (a math concept
    ships with a test of its intrinsic properties, not only of its use).

    Two laws over the whole spellable lattice, and the census of WHICH
    members grow: exactly the two connected groups, whose improper
    extension fixes every invariant.  (That the two axial groups then share
    one orbit space at the tier the tree reads — identical descending-slot
    masks — is ``test_the_MIRRORED_half_is_inert_…_THEOREM`` above.)
    """

    @staticmethod
    def _lattice() -> list[SubgroupOfO3]:
        """Every family the module realizes, one or more members each — the
        DENOMINATOR of the two universals below, written out."""
        G = SubgroupOfO3
        return (
            [G.Trivial, G.Dinfh, G.OctahedralOh, G.IcosahedralIh, G.SO3, G.O3]
            + [G.Cn(n) for n in (1, 2, 3, 4, 6)]
            + [G.Dnh(n) for n in (1, 2, 3, 4)]
            + [G.Mirror(a) for a in "xyz"]
            + [G.SO2(a) for a in "xyz"]
            + [G.O2(a) for a in "xyz"]
        )

    def test_it_CONTAINS_the_group_and_is_IDEMPOTENT(self) -> None:
        """Closure and idempotence: ``G ⊆ Stab(G)`` and ``Stab(Stab(G)) =
        Stab(G)`` — a stabiliser is its own stabiliser — on 24 of 24."""
        groups = self._lattice()
        assert len(groups) == 24
        for g in groups:
            s = g.orbit_stabiliser
            assert s.contains(g), (g, s)
            assert s.orbit_stabiliser == s, (g, s)

    def test_exactly_the_two_CONNECTED_groups_grow(self) -> None:
        """`[M]` the non-identity cases are ``SO2(a) → O2(a)`` on each axis and
        ``SO3 → O3``; every finite group, ``D_∞h``, ``O2(a)`` and ``O3`` is
        its own stabiliser (24 of 24 enumerated above)."""
        grown = {
            repr(g): repr(g.orbit_stabiliser)
            for g in self._lattice()
            if g.orbit_stabiliser != g
        }
        assert grown == {
            "SubgroupOfO3.SO2('x')": "SubgroupOfO3.O2('x')",
            "SubgroupOfO3.SO2('y')": "SubgroupOfO3.O2('y')",
            "SubgroupOfO3.SO2('z')": "SubgroupOfO3.O2('z')",
            "SubgroupOfO3.SO3": "SubgroupOfO3.O3",
        }


# =============================================================================
# 2.2b — the Γ-slot: gates drafted by the test-architect (2026-09-02)
# =============================================================================

_g22b_AXES = ("x", "y", "z")


def _g22b_trivially_quotiented_product() -> DiscreteMeasure:
    from dataclasses import replace

    return replace(
        Quadrature.product(4, 8).measure,
        support=SPHERE.quotient(SubgroupOfO3.Trivial),
    )




def _g22b_rot(axis: str, theta: float) -> RigidMotion:
    return RigidMotion.rotation_about_axis(axis=np.eye(3)["xyz".index(axis)], angle=theta)


def _g22b_mirror(axis: str) -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(3)["xyz".index(axis)])


#: Pairwise-incommensurate angles.  A right-angle sample of a continuous
#: family generates C_4 and certifies what it never tested (ERR-072,
#: ``vv-principles`` #13); the CONTROL below is allowed to sample because it
#: can only REFUTE, and it is compared against a criterion that never does.
_g22b_INCOMMENSURATE = (1.0, float(np.sqrt(2.0)), float(np.e), 2.5, float(np.sqrt(7.0)))


# =============================================================================


# G5 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheNormaliserCriterionIsEXACTAndAgreesWithBruteConjugation:
    r"""(iii) ``is_normalised_by`` / ``normalises``, positive AND negative
    legs per family, with a brute-conjugation CONTROL.

    The production criterion is **exact and never sampled** (ERR-072): a
    motion normalises the axial families iff it maps ``ê_a`` to ``±ê_a``, and
    the finite families are decided by conjugating the ELEMENT SET.

    ⛔ A 45-row "brute-conjugation control" stood here until 2026-09-03 and
    was RETIRED by R1 of #434: the review AST-unparsed it and found it
    α-identical to production's finite arm, with its element list read from
    production's own realization — a control that could not disagree
    (``vv-principles`` #22, the shared-EXPRESSION case).  The independent
    reference — every element set rebuilt from the group's DEFINITION in
    plain numpy, 31 members × 46 motions — is
    ``TestR1TheNormaliserIsExactAndComputed``.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_an_axial_group_is_normalised_exactly_by_the_axis_preserving_motions(
        self, axis
    ):
        for H in (SubgroupOfO3.O2(axis), SubgroupOfO3.SO2(axis)):
            assert H.is_normalised_by(_g22b_rot(axis, float(np.e)))       # rotation about a
            assert H.is_normalised_by(_g22b_mirror(axis))                 # flips ê_a -> -ê_a
            for other in _g22b_AXES:
                if other == axis:
                    continue
                assert not H.is_normalised_by(_g22b_rot(other, 1.0))
                assert H.is_normalised_by(_g22b_mirror(other))            # fixes ê_a

    @pytest.mark.foundation
    def test_the_GROUP_level_normaliser_positive_and_negative_legs(self):
        r"""`[M]` every row measured 2026-09-02 on the shadow of the design.

        ⚠ ``O2_x.normalises(σ_y)`` is **False** although ``O2_x CONTAINS
        σ_y`` — containment does not imply normalisation for a subgroup that
        is not normal, and a rotation about x conjugates σ_y to a mirror in a
        tilted plane.  That pair is the reason step 1 (normalise) must run
        BEFORE step 2 (contains), and it is the sharpest row here.
        """
        O2x, SO2x = SubgroupOfO3.O2("x"), SubgroupOfO3.SO2("x")
        sx, sy, sz = (SubgroupOfO3.Mirror(a) for a in _g22b_AXES)
        assert O2x.contains(sy) and not O2x.normalises(sy)      # <- the sharp row
        assert O2x.normalises(sx) and O2x.normalises(SO2x) and O2x.normalises(O2x)
        assert O2x.normalises(SubgroupOfO3.Trivial)
        assert not O2x.normalises(SubgroupOfO3.Cn(4))
        assert SubgroupOfO3.Dinfh.normalises(SubgroupOfO3.O2("z"))
        assert SubgroupOfO3.Dinfh.normalises(sz)
        assert not SubgroupOfO3.Dinfh.normalises(O2x)
        assert SubgroupOfO3.SO3.normalises(SubgroupOfO3.Trivial)
        assert not SubgroupOfO3.SO3.normalises(sz)
        assert not SubgroupOfO3.SO3.normalises(SubgroupOfO3.O2("z"))
        for a in _g22b_AXES:                       # D_2h normalises each of its mirrors
            assert SubgroupOfO3.Dnh(2).normalises(SubgroupOfO3.Mirror(a))
        assert SubgroupOfO3.Dnh(2).normalises(O2x)
        assert not SubgroupOfO3.Cn(4).normalises(sy)
        assert SubgroupOfO3.Cn(4).normalises(sz)
        assert SubgroupOfO3.SO2("y").normalises(sy)
        assert not SubgroupOfO3.SO2("z").normalises(sy)

    @pytest.mark.foundation
    def test_normalisation_is_REFLEXIVE_on_every_spellable_group(self):
        r"""⛔ **The design's spec, read literally, gets this WRONG for the
        continuous families.**  It prescribes ``H ⊆ {e, −I}`` for ``SO(3)``,
        which makes ``SO3.normalises(SO3)`` and ``SO3.normalises(O3)``
        answer ``False`` — and every group normalises itself, by theorem.

        `[M]` the shadow written from the design text answers ``False`` on
        both.  It is LATENT for the shipped inputs (a ``Quotient.by`` is only
        ever ``O2_a``, ``σ_a`` or ``Trivial``), which is exactly why it needs
        a gate rather than a fix note: nothing else in the tree will ever ask.
        """
        for g in _every_spellable_group():
            assert g.normalises(g), g.name

    @pytest.mark.foundation
    def test_the_two_spellings_of_the_normaliser_agree_where_both_apply(self):
        """``G.normalises(H)`` must be ``all(H.is_normalised_by(g))`` over a
        finite ``G`` — one relation, two entry points, no drift."""
        for G in (
            SubgroupOfO3.Dnh(2),
            SubgroupOfO3.Cn(4),
            SubgroupOfO3.Mirror("x"),
            SubgroupOfO3.OctahedralOh,
        ):
            for H in (
                SubgroupOfO3.Mirror("y"),
                SubgroupOfO3.O2("x"),
                SubgroupOfO3.Trivial,
            ):
                assert G.normalises(H) is all(
                    H.is_normalised_by(g) for g in _closed(G)
                ), (G.name, H.name)

    @pytest.mark.foundation
    def test_a_CONTINUOUS_normaliser_is_not_contradicted_by_any_sampled_element(self):
        """A necessary condition on the exact criterion: if ``G.normalises(H)``
        then no element of ``G`` — sampled at incommensurate angles, which a
        certificate may not do but a REFUTATION may — can fail
        ``H.is_normalised_by``."""
        for axis in _g22b_AXES:
            for G in (SubgroupOfO3.SO2(axis), SubgroupOfO3.O2(axis)):
                for H in (SubgroupOfO3.O2(axis), SubgroupOfO3.Mirror(axis)):
                    if not G.normalises(H):
                        continue
                    for theta in _g22b_INCOMMENSURATE:
                        assert H.is_normalised_by(_g22b_rot(axis, theta)), (G.name, H.name)
                    for other in _g22b_AXES:
                        assert H.is_normalised_by(_g22b_mirror(other)) or other == axis


def _every_spellable_group() -> list[SubgroupOfO3]:
    named = [
        SubgroupOfO3.Trivial,
        SubgroupOfO3.Dinfh,
        SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh,
        SubgroupOfO3.SO3,
        SubgroupOfO3.O3,
    ]
    named += [SubgroupOfO3.Mirror(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.SO2(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.O2(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.Cn(n) for n in (1, 2, 3, 4)]
    named += [SubgroupOfO3.Dnh(n) for n in (1, 2, 3, 4)]
    return named


# =============================================================================


# G6 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheTrivialQuotientAnswersExactlyLikeTheBareSphere:
    r"""(v) The bit-identity witness that the orbit-space route CONTAINS the
    ambient kernel.

    ``S²/{e} = S²`` is a theorem (``_mod_trivial``: ``P = I``), so routing a
    trivially-quotiented measure through the kernel (``invariance.is_invariant_under``
    since R2 of #434; ``_invariance_on_orbit_space`` until then) must
    reproduce ``_invariance_on_points`` exactly — for the finite families
    (step 5 becomes plain orbit closure on the base's own coordinates) AND
    for the continuous ones (step 4 becomes the axis-support / origin-support
    rule the ambient kernel already applies).

    `[M]` 2026-09-02: **85 rows** (5 shipped sphere rules × 17 groups), **0
    differ**.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "name,factory",
        [
            ("product(4,8)", lambda: Quadrature.product(4, 8).measure),
            ("product(2,4)", lambda: Quadrature.product(2, 4).measure),
            ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4).measure),
            ("lebedev(5)", lambda: Quadrature.lebedev(5).measure),
            ("lebedev(11)", lambda: Quadrature.lebedev(11).measure),
        ],
    )
    def test_every_candidate_group_answers_identically(self, name, factory):
        from dataclasses import replace

        bare = factory()
        assert bare.support == SPHERE
        trivial = replace(bare, support=SPHERE.quotient(SubgroupOfO3.Trivial))
        assert isinstance(trivial.support, Quotient)
        checked = 0
        for g in _every_spellable_group():
            checked += 1
            assert bare.is_invariant_under(g) is trivial.is_invariant_under(g), (name, g.name)
        assert checked >= 17, checked


# =============================================================================


# G7 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheSlabAndSphereFamiliesAnswerExactlyAsBeforeTheCarve:
    r"""(vi) The PRE-CHANGE baseline, as literals.

    Measured on the pristine tree at ``4b7d24c3`` on 2026-09-02, over each
    rule's OWN ``candidate_groups`` set — not a hand list, which is how the
    first pass of this measurement under-counted the fold's flips by two
    (``D_1h`` and ``C_2`` were absent from the hand list).

    `[M]` the whole carve moves **20 of 415** rows over 21 shipped rules,
    every one of them on a σ_y FOLD.  Nothing on a slab rule, a polar
    marginal on any axis, a product rule, a level-symmetric rule or a
    Lebedev rule moves at all.

    ⭐ R2 of #434 (2026-09-03): the candidate set is derived from the orbit
    BARYCENTRES, so a 1-D rule is offered the polyhedral groups its
    barycentres admit — the slab's ``(mu, 0, 0)`` have two z-azimuths and
    gain ``C_2``/``D_1h``/``D_2h``, the z-marginal's ``(0, 0, mu)`` lie ON
    the axis and gain ``D_1h``.  `[M]` every one of the new rows reads
    ``True`` and **no shared key moved** (15 of 15 identical on each rule);
    the literals below carry the new keys, marked.
    """

    _SLAB_GL8 = {
        "O2_x": True, "SO2_x": True, "Trivial": True,
        "sigma_x": True, "sigma_y": True, "sigma_z": True,
        "C_2": True, "D_1h": True, "D_2h": True,          # offered since R2
        "Dinfh": False, "Ih": False, "O2_y": False, "O2_z": False,
        "O3": False, "Oh": False, "SO2_y": False, "SO2_z": False, "SO3": False,
    }
    _POLAR_Z = {
        "Dinfh": True, "O2_z": True, "SO2_z": True, "Trivial": True,
        "sigma_x": True, "sigma_y": True, "sigma_z": True,
        "D_1h": True,                                     # offered since R2
        "Ih": False, "O2_x": False, "O2_y": False, "O3": False, "Oh": False,
        "SO2_x": False, "SO2_y": False, "SO3": False,
    }

    @pytest.mark.foundation
    def test_the_slab_rule_is_bit_identical_to_the_pre_change_answers(self):
        m = Quadrature.gauss_legendre(8).measure
        got = {g.name: m.is_invariant_under(g) for g in candidate_groups(m)}
        assert got == self._SLAB_GL8

    @pytest.mark.foundation
    def test_a_polar_marginal_about_z_is_bit_identical_too(self):
        r"""⭐ The ``S^2/O2_z`` marginal is the ONLY shipped input that
        reaches step 4 with a CONTINUOUS quotienting group: ``D_∞h``
        normalises ``O(2)_z`` and is not contained in it, so the verdict is
        decided by the identity component fixing every chart node and then by
        the coset representatives — not by the step-2 short circuit.  `[M]`
        ``Dinfh`` reads ``True`` here and ``False`` on an ASYMMETRIC marginal
        on the same entry, so the row discriminates."""
        m = gauss_legendre_on_polar_orbit(8, "z")
        got = {g.name: m.is_invariant_under(g) for g in candidate_groups(m)}
        assert got == self._POLAR_Z
        asym = DiscreteMeasure(
            nodes=np.array([-0.9, -0.3, 0.3, 0.7]),
            weights=np.full(4, 0.5),
            support=SPHERE.quotient(SubgroupOfO3.O2("z")),
        )
        assert asym.is_invariant_under(SubgroupOfO3.O2("z"))
        assert not asym.is_invariant_under(SubgroupOfO3.Dinfh)

    @pytest.mark.foundation
    @pytest.mark.parametrize("n", [2, 3, 4, 8, 16])
    @pytest.mark.parametrize("axis", _g22b_AXES)
    def test_no_axial_marginal_on_any_axis_moves(self, axis, n):
        """`[M]` 3 axes × 21 groups = 63 rows, 0 changed."""
        m = gauss_legendre_on_polar_orbit(n, axis)
        assert m.is_invariant_under(SubgroupOfO3.O2(axis))
        assert m.is_invariant_under(SubgroupOfO3.SO2(axis))
        assert m.is_invariant_under(SubgroupOfO3.Mirror(axis))   # GL nodes are ±-paired
        for other in _g22b_AXES:
            if other == axis:
                continue
            assert not m.is_invariant_under(SubgroupOfO3.O2(other))
            assert not m.is_invariant_under(SubgroupOfO3.SO2(other))
            assert m.is_invariant_under(SubgroupOfO3.Mirror(other))  # fixes the axis


# =============================================================================


# G8 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheFoldIsInvariantUnderExactlyTheGroupsThatDescendAndClose:
    r"""The carve's own consequence — and the DENOMINATOR the design memo did
    not carry.

    ⛔ **The design's §8 names two flips (``Mirror('y')`` and ``Dnh(2)``).
    `[M]` there are FOUR**, on every shipped fold: ``σ_y``, ``C_2``,
    ``D_1h`` and ``D_2h``, all ``False → True``.  Measured over each rule's
    own candidate set: 5 folds × 4 = **20 of 415** rows.

    ⚠ **And the headline flip is the least falsifiable one.**  ``σ_y`` is
    answered by the step-2 short circuit (``H.contains(G)`` — G acts
    trivially on ``M/H``), which never looks at a node.  `[M]` it stays
    ``True`` on a fold with a node DELETED and on a fold with a weight
    perturbed by 1.5×, while ``D_2h`` and ``σ_x`` both go ``False`` on both.
    So ``σ_y`` gates the WIRING and ``D_2h`` gates the node-level closure;
    citing the σ_y row as evidence of the latter is ``vv-principles`` #19.
    """

    # C_1 IS Trivial since 2026-09-02 (R1 of #434): one row for {e}, not two.
    # C_4 / D_4h are NOT offered since 2026-09-03 (R2 of #434): the candidate
    # set is read off the barycentres, whose y column the fold's projection
    # zeroes — two z-azimuths, not four (18 candidates, 20 until R2).
    _FOLD_AFTER = {
        "C_2": True, "D_1h": True, "D_2h": True, "Trivial": True,
        "sigma_x": True, "sigma_y": True, "sigma_z": True,
        "Dinfh": False, "Ih": False,
        "O2_x": False, "O2_y": False, "O2_z": False, "O3": False, "Oh": False,
        "SO2_x": False, "SO2_y": False, "SO2_z": False, "SO3": False,
    }

    @pytest.mark.foundation
    def test_the_whole_verdict_set_on_the_shipped_fold(self):
        m = Quadrature.folded_product(4, 8).measure
        assert isinstance(m.support, Quotient)
        got = {g.name: m.is_invariant_under(g) for g in candidate_groups(m)}
        assert got == self._FOLD_AFTER

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8), (6, 8), (8, 8)])
    def test_all_four_flips_hold_on_every_shipped_fold(self, n_mu, n_phi):
        m = Quadrature.folded_product(n_mu, n_phi).measure
        for name in ("sigma_y", "C_2", "D_1h", "D_2h"):
            g = next(c for c in candidate_groups(m) if c.name == name)
            assert m.is_invariant_under(g), name

    @pytest.mark.foundation
    def test_the_node_reading_flips_are_falsifiable_and_sigma_y_is_NOT(self):
        r"""The discriminating half of the flip, and the honest admission
        about the other half.  `[M]` both perturbations below move ``D_2h``,
        ``C_2``, ``D_1h`` and ``σ_x`` to ``False`` and leave ``σ_y`` ``True``
        — because σ_y acts trivially on ``S²/σ_y`` for EVERY measure there,
        which is a theorem about the space, not a fact about the rule."""
        m = Quadrature.folded_product(4, 8).measure
        nodes = np.asarray(m.nodes, dtype=float)
        weights = np.asarray(m.weights, dtype=float)

        bumped = weights.copy()
        bumped[0] *= 1.5
        perturbed = DiscreteMeasure(nodes=nodes, weights=bumped, support=m.support)
        dropped = DiscreteMeasure(
            nodes=nodes[1:], weights=weights[1:], support=m.support
        )
        for broken in (perturbed, dropped):
            assert not broken.is_invariant_under(SubgroupOfO3.Dnh(2))
            assert not broken.is_invariant_under(SubgroupOfO3.Cn(2))
            assert not broken.is_invariant_under(SubgroupOfO3.Mirror("x"))
            # unfalsifiable BY CONSTRUCTION — stated, not hidden
            assert broken.is_invariant_under(SubgroupOfO3.Mirror("y"))

    @pytest.mark.foundation
    def test_a_group_that_does_not_NORMALISE_the_fold_is_refused_whatever_the_nodes(
        self,
    ):
        r"""Step 1 is a statement about the GROUPS, not the measure: ``C_4``
        about z conjugates ``σ_y`` to ``σ_x``, so it does not act on
        ``S²/σ_y`` at all and the verdict is ``False`` for every rule on that
        space — including the one whose PARENT is ``C_4``-invariant.  `[M]`
        ``product(4, 8)`` is ``C_4``-invariant; its σ_y fold is not."""
        parent = Quadrature.product(4, 8).measure
        fold = Quadrature.folded_product(4, 8).measure
        assert parent.is_invariant_under(SubgroupOfO3.Cn(4))
        assert not fold.is_invariant_under(SubgroupOfO3.Cn(4))
        assert SubgroupOfO3.Dnh(2).normalises(SubgroupOfO3.Mirror("y"))
        assert not SubgroupOfO3.Cn(4).normalises(SubgroupOfO3.Mirror("y"))


# =============================================================================


# G9 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheCompatibilityLawHoldsOnTheOrbitSpaceRoute:
    r"""(iv) ``vv-principles`` #15 — ``A ⊆ B ∧ P(B, μ) ⟹ P(A, μ)``, re-run
    over the lattice × the ORBIT-SPACE route.

    ⛔ **The module's existing downward-closure gates cannot reach this
    carve.**  `[M]` ``_walk_rules()`` supplies 3 products + 2 level-symmetric
    + 2 Lebedev and `_measure_from_sphere_quad` re-declares
    ``support=SPHERE`` on every one — so the whole existing instrument runs
    on the ambient kernel and is structurally blind to the orbit-space one.
    This class is the missing denominator.

    `[M]` 2026-09-02, 0 violations on the fold, the slab rule, the three
    axial marginals and a trivial-quotient sphere rule.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "name,factory",
        [
            ("folded_product(4,8)", lambda: Quadrature.folded_product(4, 8).measure),
            ("folded_product(2,4)", lambda: Quadrature.folded_product(2, 4).measure),
            ("folded_product(6,8)", lambda: Quadrature.folded_product(6, 8).measure),
            ("slab gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8).measure),
            ("polar_orbit(8,z)", lambda: gauss_legendre_on_polar_orbit(8, "z")),
            ("trivial-quotient product(4,8)", _g22b_trivially_quotiented_product),
        ],
    )
    def test_invariance_is_downward_closed_on_the_orbit_space_route(self, name, factory):
        measure = factory()
        assert isinstance(measure.support, Quotient), name
        cands = candidate_groups(measure)
        invariant = {repr(g): measure.is_invariant_under(g) for g in cands}
        checked = 0
        violations = []
        for outer in cands:
            for inner in cands:
                if outer == inner or not outer.contains(inner):
                    continue
                checked += 1
                if invariant[repr(outer)] and not invariant[repr(inner)]:
                    violations.append(f"{name}: {outer.name} > {inner.name}")
        assert checked > 40, f"{name}: only {checked} pairs — the gate is too thin"
        assert not violations, violations[:8]


# =============================================================================


# G10 — tests/numerics/test_symmetry.py
# =============================================================================
class TestTheWalkOnAFoldReportsTheGroupTheGeometryOwes:
    r"""(x) The Hasse walk on a fold — measured BEFORE, and predicted AFTER
    from the Γ-closure on the orbit space.

    `[M]` before the carve: ``['sigma_z', 'sigma_x']`` — two incomparable
    maxima, because the ambient kernel could not see that σ_y is free and
    therefore could not close the pair up to ``D_2h``.
    `[M]` after: ``['D_2h']``, on all five shipped folds, with the walk and
    the bruteforce scan AGREEING (the module's own two-realization
    instrument, which is why the answer is derivable and not merely pinned).

    ⭐ ``candidate_groups`` reads ``_distinct_azimuths`` of the orbit
    BARYCENTRES (R2 of #434 — until then of the stored nodes, so one fold had
    two candidate sets by spelling) and never asks the invariance predicate —
    so the move lives in the walk's verdicts and, on a fold, in a search
    space that no longer offers ``C_4`` / ``D_4h`` (the barycentre projection
    zeroes the fold's ``y`` column: two azimuths, not four).

    ⚠ The walk (``symmetry_groups``; ``maximal_invariance_groups`` until R2 of
    #434) has **0 production callers** (`[M]` 2026-09-03 the only ``orpheus/``
    hits are its ``def`` in ``invariance`` and the measure verb that delegates
    to it), so this change cannot move a solve.  Nothing in the tree pinned the
    walk on a fold before this class.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8), (6, 8), (8, 8)])
    def test_the_walk_and_the_bruteforce_scan_agree_and_report_D2h(self, n_mu, n_phi):
        m = Quadrature.folded_product(n_mu, n_phi).measure
        walk = sorted(g.name for g in m.symmetry_groups(method="walk"))
        brute = sorted(
            g.name for g in m.symmetry_groups(method="bruteforce")
        )
        assert walk == brute, (walk, brute)
        assert walk == ["D_2h"]

    @pytest.mark.foundation
    def test_the_candidate_SET_is_derived_from_the_barycentres(self):
        """`[M]` 20 names until R2 of #434; the barycentre projection zeroes
        the fold's y column, so its azimuth count falls 4 -> 2 and ``C_4`` /
        ``D_4h`` are no longer offered.  (Named "untouched by the carve" until
        2026-09-03, which R2 made false.)"""
        m = Quadrature.folded_product(4, 8).measure
        names = sorted(g.name for g in candidate_groups(m))
        assert names == sorted(
            [
                # C_1 IS Trivial since 2026-09-02 (R1 of #434): one group, one spelling.
                "C_2", "D_1h", "D_2h", "Dinfh", "Ih",
                "O2_x", "O2_y", "O2_z", "O3", "Oh", "SO2_x", "SO2_y", "SO2_z",
                "SO3", "Trivial", "sigma_x", "sigma_y", "sigma_z",
            ]
        )


# =============================================================================


# G11 — tests/numerics/test_symmetry.py
# =============================================================================
class TestOrbitCertificateAndInvarianceNowAnswerDifferentQuestions:
    r"""⛔ **A divergence the carve CREATES, and the message that goes
    present-tense-false with it.**

    ``certificate_under`` reads ``measure.nodes`` in the AMBIENT space, so
    after the carve ``is_invariant`` and ``certificate_under`` answer
    different questions on a fold: `[M]` ``σ_y`` / ``C_2`` / ``D_2h`` are
    ``is_invariant = True`` with ``certificate_under = None``, while ``σ_x``
    has both.  That is CORRECT — a group acting trivially on the orbit space
    permutes no ambient node — but two things in the tree assert otherwise:

    * ``certificate_under``'s docstring gives ``None`` exactly two causes
      ("not invariant" / "continuous"); there is now a THIRD;
    * ``DiscreteMeasure.quotient``'s refusal says *"this measure is not
      σ_y-invariant"*, which the carve makes FALSE.

    Both are ``coding-standards`` MUST-FIXes, in the carve's own commit.  The
    honest message names the mechanism: *no finite permutation of the NODES
    realizes σ_y; the measure IS σ_y-invariant on its orbit space, where σ_y
    acts trivially, so a second fold by the same group is the identity.*
    """

    @pytest.mark.foundation
    def test_the_certificate_and_the_verdict_cannot_disagree_on_a_fold(self):
        r"""✅ **Re-posed 2026-09-02.**  The draft predicted a DIVERGENCE
        (``is_invariant`` True with ``certificate_under`` None) because the
        certificate read the AMBIENT nodes.  `[M]` the coordinator closed it
        the better way instead: ``certificate_under`` now builds on the
        CHART, through ``induced_action``, via ``invariance._orbit_closure(...,
        images_of=)`` (the action REQUIRED since R2 of #434 — its ambient
        default was dead) — so there is ONE invariance notion, which is R1.

        What is owed is therefore the AGREEMENT, gated: on a fold, for every
        candidate group, ``is_invariant`` and ``certificate_under is not
        None`` must answer the same — except for the continuous groups, which
        have no finite realization to permute anything with.
        """
        from orpheus.numerics.symmetry import _realized_ops
        from orpheus.numerics.invariance import certificate_under

        m = Quadrature.folded_product(4, 8).measure
        checked = 0
        for g in candidate_groups(m):
            if _realized_ops(g._tag) is None:   # continuous: no realization
                assert m.certificate_under(g) is None
                continue
            checked += 1
            assert (m.certificate_under(g) is not None) is m.is_invariant_under(g), g.name
        # 9 finite candidates since R2 of #434 (C_4 / D_4h left the fold's set).
        assert checked >= 9, checked
        # and the permutations are the induced ones: sigma_y acts TRIVIALLY on
        # S^2/sigma_y, so its chart permutation is the identity.
        cert = m.certificate_under(SubgroupOfO3.Mirror("y"))
        assert cert is not None
        for perm in cert.permutations:
            assert perm.indices.tolist() == list(range(m.n_points))
        # the control: sigma_x MOVES the chart, so its permutation is not id
        cert_x = m.certificate_under(SubgroupOfO3.Mirror("x"))
        assert cert_x is not None
        assert any(
            p.indices.tolist() != list(range(m.n_points)) for p in cert_x.permutations
        )

    @pytest.mark.foundation
    def test_the_refold_refusal_no_longer_claims_the_measure_is_not_invariant(self):
        r"""✅ Already repaired by the coordinator in flight: `[M]`
        ``tests/numerics/test_measure.py::test_the_fold_consumes_the_symmetry_idempotent_only_on_a_trivial_action``
        now expects ``"lies in the spent group sigma_y"``.  This row is the
        SIBLING claim on the ``Quadrature``-tier fold, which that test does
        not reach — keep both, or the tier with no witness is the one that
        drifts back."""
        m = Quadrature.folded_product(4, 8).measure
        with pytest.raises(ValueError) as excinfo:
            m.quotient(SubgroupOfO3.Mirror("y"))
        message = str(excinfo.value)
        assert "not sigma_y-invariant" not in message, (
            "the refusal still asserts non-invariance, which tracker 2.2b "
            f"made false: {message!r}"
        )
        assert "spent" in message

    @pytest.mark.foundation
    def test_orbit_certificates_docstring_gives_None_its_THIRD_cause(self):
        r"""⛔ A docstring MUST-FIX the carve creates.  ``certificate_under``
        documents ``None`` as meaning *"not group-invariant, or the group is
        CONTINUOUS"* — after 2.2b there is a third cause, and it is the one a
        fold hits.  A gate on prose is unusual; it is here because the claim
        is a COVERAGE claim (an audit reads that docstring as the function's
        contract) and nothing else can see it."""
        from orpheus.numerics.invariance import certificate_under

        doc = certificate_under.__doc__ or ""
        assert "orbit space" in doc or "chart" in doc, (
            "certificate_under's docstring still describes the AMBIENT "
            "reading; tracker 2.2b moved its body onto the CHART (through "
            "Quotient.induced_action), and added a third cause of None (the "
            "group does not normalise the quotienting group). Its two "
            "downstream readers -- `singular_set_under` and "
            "`DiscreteMeasure.quotient` -- inherit the new semantics."
        )


# =============================================================================


def _closed(group: SubgroupOfO3) -> list[RigidMotion]:
    """The group's realized ELEMENTS (the closure), for the brute-force
    conjugation controls above — the draft defined it beside the manifold
    gates; the symmetry gates need it too."""
    from orpheus.numerics.symmetry import _group_elements

    elems = _group_elements(group._tag)
    assert elems is not None
    return list(elems)


# =============================================================================
# ===== R1 gates ==============================================================
# =============================================================================
#
# R1 — every question about a group is computed from its realization (#434).
# The reference constructions below are INDEPENDENT of the module: see this
# file's header for what "independent" means here and where it stops.


_r1_ATOL = 1e-9
_r1_AXES = ("x", "y", "z")

#: Pairwise-incommensurate angles.  A right-angle sample of a continuous family
#: generates ``C_4`` and certifies what it never tested (ERR-072,
#: ``vv-principles`` #13); every continuous-group SAMPLE below is drawn from
#: these, and it is used only where a sample can REFUTE.
_r1_INCOMMENSURATE = (1.0, float(np.sqrt(2.0)), float(np.e), 2.5,
                      float(np.sqrt(7.0)), 0.37)


# --------------------------------------------------------------------------
# Reference constructions — plain numpy, from each group's definition
# --------------------------------------------------------------------------


def _r1_axis(letter: str) -> np.ndarray:
    return np.eye(3)["xyz".index(letter)]


def _r1_skew(u) -> np.ndarray:
    u = np.asarray(u, float)
    return np.array([[0.0, -u[2], u[1]], [u[2], 0.0, -u[0]], [-u[1], u[0], 0.0]])


def _r1_rot(u, theta: float) -> np.ndarray:
    """Rodrigues: ``exp(theta [u]_x)`` — not the module's circle-point form."""
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    K = _r1_skew(u)
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def _r1_rot_z(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _r1_refl(normal) -> np.ndarray:
    n = np.asarray(normal, float)
    n = n / np.linalg.norm(n)
    return np.eye(3) - 2.0 * np.outer(n, n)


def _r1_half_turn(u) -> np.ndarray:
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    return 2.0 * np.outer(u, u) - np.eye(3)


def _r1_cn(n: int) -> np.ndarray:
    return np.array([_r1_rot_z(2.0 * np.pi * k / n) for k in range(n)])


def _r1_dnh(n: int) -> np.ndarray:
    r""":math:`D_{nh} = D_n \times \{e, \sigma_h\}`, order :math:`4n`.

    Written from the DEFINITION: the ``n`` proper rotations about z, the ``n``
    half-turns about the in-plane axes at azimuth :math:`k\pi/n` (the standard
    setting puts a vertex on x, so :math:`k = 0` is the x-axis), and each of
    those times :math:`\sigma_h` — which turns the half-turns into the vertical
    mirrors.  Production builds the same group as ``C_n + [σ_z] + C_n σ_0``
    closed under composition; neither construction can be derived from the
    other by inspection, and `[M]` they agree to 1.166e-15.
    """
    sigma_h = np.diag([1.0, 1.0, -1.0])
    rotations = [_r1_rot_z(2.0 * np.pi * k / n) for k in range(n)]
    out = list(rotations) + [r @ sigma_h for r in rotations]
    for k in range(n):
        alpha = np.pi * k / n
        out.append(_r1_half_turn([np.cos(alpha), np.sin(alpha), 0.0]))
        out.append(_r1_refl([-np.sin(alpha), np.cos(alpha), 0.0]))
    return np.array(out)


def _r1_oh() -> np.ndarray:
    """The 48 signed permutation matrices, written out by index assignment."""
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            M = np.zeros((3, 3))
            for i, j in enumerate(perm):
                M[i, j] = signs[i]
            out.append(M)
    return np.array(out)


def _r1_icosahedron() -> np.ndarray:
    """The 12 vertices in the golden-ratio convention, unit-normalised."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = np.array([
        [0, 1, phi], [0, 1, -phi], [0, -1, phi], [0, -1, -phi],
        [1, phi, 0], [1, -phi, 0], [-1, phi, 0], [-1, -phi, 0],
        [phi, 0, 1], [phi, 0, -1], [-phi, 0, 1], [-phi, 0, -1],
    ], dtype=float)
    return raw / np.linalg.norm(raw[0])


def _r1_ih() -> np.ndarray:
    r""":math:`I_h` by the FLAG construction — a different algorithm from
    production's generator closure.

    A rotation of the icosahedron is determined by the image of one flag
    (a vertex together with one of its five neighbours): the frame
    :math:`[v, w, v\times w]` maps to :math:`[v', w', v'\times w']`, so the
    rotation is :math:`F' F_0^{-1}`.  Twelve vertices times five neighbours is
    the rotation group :math:`I` (60); :math:`I_h = I \sqcup (-I)\,I` (120),
    because :math:`-I` is central and not a rotation.
    """
    V = _r1_icosahedron()
    d2 = ((V[:, None, :] - V[None, :, :]) ** 2).sum(-1)
    nearest = np.sort(d2[0])[1]
    adjacency = [np.flatnonzero(np.abs(d2[i] - nearest) < 1e-9) for i in range(12)]
    v0, w0 = V[0], V[adjacency[0][0]]
    F0_inverse = np.linalg.inv(np.column_stack([v0, w0, np.cross(v0, w0)]))
    rotations = np.array([
        np.column_stack([V[i], V[j], np.cross(V[i], V[j])]) @ F0_inverse
        for i in range(12) for j in adjacency[i]
    ])
    return np.concatenate([rotations, -rotations])


def _r1_finite_reference() -> "dict[str, np.ndarray]":
    """Every finite member's element set, by construction.  Memoised at module
    scope (`[M]` 2.9 ms to build; the ``I_h`` half is the only non-trivial cost)."""
    if _r1_FINITE_CACHE:
        return _r1_FINITE_CACHE
    _r1_FINITE_CACHE["Trivial"] = np.array([np.eye(3)])
    _r1_FINITE_CACHE["Oh"] = _r1_oh()
    _r1_FINITE_CACHE["Ih"] = _r1_ih()
    for a in _r1_AXES:
        _r1_FINITE_CACHE[f"sigma_{a}"] = np.array([np.eye(3), _r1_refl(_r1_axis(a))])
    for n in range(1, 9):
        _r1_FINITE_CACHE[f"C_{n}"] = _r1_cn(n)
        _r1_FINITE_CACHE[f"D_{n}h"] = _r1_dnh(n)
    return _r1_FINITE_CACHE


_r1_FINITE_CACHE: "dict[str, np.ndarray]" = {}

#: The ROSTER — every member the module can spell today, one row per distinct
#: value the constructors accept.  ``C_1`` is kept in the list deliberately: R1
#: makes it EQUAL to ``Trivial``, so its row must reproduce ``Trivial``'s
#: answers exactly, and a table that dropped it could not see that.
_r1_MEMBERS = (
    ["Trivial", "Dinfh", "Oh", "Ih", "SO3", "O3"]
    + [f"sigma_{a}" for a in _r1_AXES]
    + [f"SO2_{a}" for a in _r1_AXES]
    + [f"O2_{a}" for a in _r1_AXES]
    + [f"C_{n}" for n in range(1, 9)]
    + [f"D_{n}h" for n in range(1, 9)]
)


def _r1_group(name: str) -> SubgroupOfO3:
    """The production group a roster label names."""
    G = SubgroupOfO3
    fixed = {"Trivial": G.Trivial, "Dinfh": G.Dinfh, "Oh": G.OctahedralOh,
             "Ih": G.IcosahedralIh, "SO3": G.SO3, "O3": G.O3}
    if name in fixed:
        return fixed[name]
    if name.startswith("sigma_"):
        return G.Mirror(name[-1])
    if name.startswith("SO2_"):
        return G.SO2(name[-1])
    if name.startswith("O2_"):
        return G.O2(name[-1])
    if name.startswith("C_"):
        return G.Cn(int(name[2:]))
    if name.startswith("D_") and name.endswith("h"):
        return G.Dnh(int(name[2:-1]))
    raise ValueError(name)  # pragma: no cover - roster guard


def _r1_is_member(name: str, matrices: np.ndarray, atol: float = _r1_ATOL) -> np.ndarray:
    r"""``(m,)`` boolean: which of the ``(m, 3, 3)`` matrices lie in the group.

    Analytic for every continuous member — :math:`O(3)` is everything,
    :math:`SO(3)` is :math:`\det > 0`, :math:`SO(2)_a` / :math:`O(2)_a` are the
    stabiliser of :math:`\hat e_a` (with / without the properness cut) and
    :math:`D_{\infty h}` is the stabiliser of the z-LINE — exact for a set with
    no finite element list.  Exact matrix comparison for the finite ones.
    """
    matrices = np.asarray(matrices, float)
    if name == "O3":
        return np.ones(len(matrices), bool)
    determinant = np.linalg.det(matrices)
    if name == "SO3":
        return determinant > 0
    if name.startswith("SO2_") or name.startswith("O2_"):
        e = _r1_axis(name[-1])
        fixes_axis = np.abs(matrices @ e - e).max(axis=1) <= atol
        return fixes_axis & (determinant > 0) if name.startswith("SO2_") else fixes_axis
    if name == "Dinfh":
        e = _r1_axis("z")
        image = matrices @ e
        return ((np.abs(image - e).max(axis=1) <= atol)
                | (np.abs(image + e).max(axis=1) <= atol))
    reference = _r1_finite_reference()[name]
    distance = np.abs(matrices[:, None] - reference[None]).max(axis=(2, 3))
    return distance.min(axis=1) <= atol


def _r1_sample(name: str) -> np.ndarray:
    """Every element (finite) or a GENERIC incommensurate sample (continuous).

    A sample of a continuous group can only REFUTE a universal, never certify
    one — which is exactly how it is used below (``contains`` and ``normalises``
    with a continuous inner/self argument).  Drawing it at INCOMMENSURATE
    angles is what stops it generating a finite subgroup and certifying by
    accident (ERR-072 / ``vv-principles`` #13).
    """
    finite = _r1_finite_reference()
    if name in finite:
        return finite[name]
    if name.startswith("SO2_") or name.startswith("O2_"):
        a = name[-1]
        out = [_r1_rot(_r1_axis(a), t) for t in _r1_INCOMMENSURATE]
        if name.startswith("O2_"):
            b, c = [i for i in range(3) if i != "xyz".index(a)]
            for t in _r1_INCOMMENSURATE:
                v = np.zeros(3)
                v[b], v[c] = np.cos(t), np.sin(t)
                out.append(_r1_refl(v))
        return np.array(out)
    if name == "Dinfh":
        sigma_h = np.diag([1.0, 1.0, -1.0])
        rotations = [_r1_rot_z(t) for t in _r1_INCOMMENSURATE]
        out = list(rotations) + [sigma_h @ r for r in rotations]
        for t in _r1_INCOMMENSURATE:
            out.append(_r1_refl([np.cos(t), np.sin(t), 0.0]))
            out.append(_r1_half_turn([np.cos(t), np.sin(t), 0.0]))
        return np.array(out)
    rng = np.random.default_rng(20260902)
    out = []
    while len(out) < 24:
        Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
        Q = Q * np.sign(np.diag(R))
        if np.linalg.det(Q) < 0:
            Q = Q @ np.diag([-1.0, 1.0, 1.0])
        assert np.linalg.det(Q) > 0        # the sampler's own positive control
        out.append(Q)
    if name == "O3":
        out = out + [-m for m in out]
    return np.array(out)


def _r1_contains(outer: str, inner: str) -> bool:
    """``inner ⊆ outer`` — the independent reference."""
    finite = _r1_finite_reference()
    if inner in finite and outer in finite:
        distance = np.abs(finite[outer][None] - finite[inner][:, None]).max(axis=(2, 3))
        return bool((distance.min(axis=1) <= _r1_ATOL).all())
    if inner not in finite and outer in finite:
        return False        # a finite group holds no infinite one
    return bool(_r1_is_member(outer, _r1_sample(inner)).all())


def _r1_normalised_by(name: str, motions: np.ndarray) -> np.ndarray:
    r"""``(m,)`` boolean: does each motion normalise the group?

    By CONJUGATION and membership — :math:`gHg^{-1} \subseteq H` element by
    element — never by an axis rule.  That is the structural independence
    that matters here: production decides the continuous families through
    :math:`\mathrm{Ad}_g\mathfrak g = \mathfrak g`, a statement about the Lie
    ALGEBRA, and this reference never touches an algebra.
    """
    S = _r1_sample(name)
    motions = np.asarray(motions, float)
    inverses = np.linalg.inv(motions)
    conjugated = np.einsum("mij,kjl,mln->mkin", motions, S, inverses)
    m, k = conjugated.shape[:2]
    inside = _r1_is_member(name, conjugated.reshape(m * k, 3, 3)).reshape(m, k)
    return np.asarray(inside.all(axis=1))


def _r1_stratified_motions() -> np.ndarray:
    r"""A motion set chosen to SEPARATE the criteria, not to enumerate a group.

    46 matrices: the identity and inversion; per coordinate axis a mirror, a
    half-turn, a quarter-turn, two incommensurate rotations and an
    :math:`S_4`; the four diagonal mirrors and diagonal half-turns; a
    three-fold about :math:`(1,1,1)` and a five-fold about an icosahedral
    vertex; and eight seeded generic rotations with their improper partners.
    `[M]` over the 31-member roster this gives **578 True / 848 False** — both
    truth values populated on 27 of 31 members.  The four CONSTANT columns are
    ``Trivial``, ``C_1``, ``SO(3)`` and ``O(3)``, all-True because those groups
    are NORMAL in :math:`O(3)`; that is a theorem, not a coverage gap
    (``vv-principles`` #20), and §G4's gate names them.
    """
    out = [np.eye(3), -np.eye(3)]
    for a in _r1_AXES:
        e = _r1_axis(a)
        out.append(_r1_refl(e))
        out.append(_r1_half_turn(e))
        out.append(_r1_rot(e, np.pi / 2))
        out.append(_r1_rot(e, 1.0))
        out.append(_r1_rot(e, float(np.sqrt(2.0))))
        out.append(-_r1_rot(e, np.pi / 2))
    for v in ([1, 1, 0], [1, -1, 0], [1, 0, 1], [0, 1, 1]):
        out.append(_r1_refl(v))
        out.append(_r1_half_turn(v))
    out.append(_r1_rot([1, 1, 1], 2 * np.pi / 3))
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    out.append(_r1_rot([0, 1, phi], 2 * np.pi / 5))
    rng = np.random.default_rng(20260902)
    for _ in range(8):
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        R = _r1_rot(u, float(rng.uniform(0.3, 3.0)))
        out.append(R)
        out.append(-R)
    return np.array(out)


def _r1_motion(matrix: np.ndarray) -> RigidMotion:
    return RigidMotion(linear=np.asarray(matrix, float), translation=np.zeros(3))


#: Every finite member's textbook order — the denominator of §G2, written out
#: so a construction that silently lost an element cannot pass by agreeing with
#: itself.  |C_n| = n, |D_nh| = 4n, |sigma| = 2, |O_h| = 48, |I_h| = 120.
_r1_ORDERS = {"Trivial": 1, "Oh": 48, "Ih": 120,
              **{f"sigma_{a}": 2 for a in _r1_AXES},
              **{f"C_{n}": n for n in range(1, 9)},
              **{f"D_{n}h": 4 * n for n in range(1, 9)}}

#: `[M]` 2026-09-02 — the corrected table.  Every FINITE member's identity
#: component is the TRIVIAL group (a discrete group's identity component is a
#: point); before R1 the property returned the group ITSELF, which is wrong on
#: 21 of the roster's 22 finite members (all but ``Trivial``).
_r1_IDENTITY_COMPONENT = {
    "Trivial": "Trivial", "Dinfh": "SO2_z", "Oh": "Trivial", "Ih": "Trivial",
    "SO3": "SO3", "O3": "SO3",
    **{f"sigma_{a}": "Trivial" for a in _r1_AXES},
    **{f"SO2_{a}": f"SO2_{a}" for a in _r1_AXES},
    **{f"O2_{a}": f"SO2_{a}" for a in _r1_AXES},
    **{f"C_{n}": "Trivial" for n in range(1, 9)},
    **{f"D_{n}h": "Trivial" for n in range(1, 9)},
}

#: `[M]` the dimension of each member as a Lie group.  Note what is ABSENT:
#: no member has dimension 2, and none can — ``so(3)`` has no two-dimensional
#: subalgebra, which is the theorem the whole realization design rests on.
_r1_DIM = {"Trivial": 0, "Dinfh": 1, "Oh": 0, "Ih": 0, "SO3": 3, "O3": 3,
           **{f"sigma_{a}": 0 for a in _r1_AXES},
           **{f"SO2_{a}": 1 for a in _r1_AXES},
           **{f"O2_{a}": 1 for a in _r1_AXES},
           **{f"C_{n}": 0 for n in range(1, 9)},
           **{f"D_{n}h": 0 for n in range(1, 9)}}


# --------------------------------------------------------------------------
# G1 — the realization is the group
# --------------------------------------------------------------------------


class TestR1TheRealizationIsWhatTheGroupIS:
    r"""The structural laws of :class:`Realization` / :class:`IdentityComponent`.

    A math-bearing type ships a test of its DEFINING laws, not only of its use.
    :math:`G = \bigsqcup_r r\,G^0` is a partition claim, ``generators`` is a
    basis claim, and both are assertable without any reference solver.
    """

    @pytest.mark.foundation
    def test_every_member_realizes_and_no_member_has_DIMENSION_TWO(self) -> None:
        r"""``dim`` over the whole roster, and the theorem behind the design.

        :math:`\mathfrak{so}(3)` is simple and three-dimensional, so its
        subalgebras are :math:`\{0\}`, the lines :math:`\mathbb R[\hat a]_\times`
        and the whole algebra — **there is no two-dimensional subalgebra**.
        That is what makes ``(identity component, coset representatives)`` a
        COMPLETE representation of a closed subgroup of :math:`O(3)`, and it is
        why ``generators`` may only ever have length 0, 1 or 3.

        `[M]` 22 members at dim 0, 7 at dim 1 (:math:`D_{\infty h}` and the six
        axial groups), 2 at dim 3 — 31 of 31, and 0 at dim 2.
        """
        seen = {}
        for name in _r1_MEMBERS:
            group = _r1_group(name)
            realization = group.realization
            assert isinstance(realization, Realization), name
            assert isinstance(realization.component, IdentityComponent), name
            n_generators = len(realization.component.generators)
            assert n_generators in (0, 1, 3), (
                f"{name}: so(3) has no {n_generators}-dimensional subalgebra"
            )
            assert realization.component.dim == n_generators == group.dim, name
            assert group.dim == _r1_DIM[name], (name, group.dim, _r1_DIM[name])
            seen[group.dim] = seen.get(group.dim, 0) + 1
        assert seen == {0: 22, 1: 7, 3: 2}, seen
        assert 2 not in seen, "a two-dimensional subalgebra of so(3) is impossible"

    @pytest.mark.foundation
    def test_the_generators_are_SKEW_and_span_what_dim_claims(self) -> None:
        r"""Each generator is in :math:`\mathfrak{so}(3)` (:math:`X^{\mathsf T} =
        -X`) and the tuple is linearly INDEPENDENT — ``dim`` is a rank, not a
        length.  A duplicated generator would pass a length check and fail this.
        """
        for name in _r1_MEMBERS:
            generators = _r1_group(name).realization.component.generators
            for X in generators:
                assert X.shape == (3, 3), name
                np.testing.assert_allclose(X, -X.T, atol=1e-14,
                                           err_msg=f"{name}: generator not skew")
            if generators:
                flat = np.array([X.ravel() for X in generators])
                assert np.linalg.matrix_rank(flat, tol=1e-9) == len(generators), (
                    f"{name}: {len(generators)} generators of rank "
                    f"{np.linalg.matrix_rank(flat, tol=1e-9)}"
                )

    @pytest.mark.foundation
    def test_the_identity_comes_first_and_the_cosets_PARTITION(self) -> None:
        r""":math:`G = \bigsqcup_r r\,G^0` is a DISJOINT union, and the module's
        contract says the identity is representative zero.

        Both halves are load-bearing.  The order matters because consumers read
        ``representatives[0]`` as the identity coset; the disjointness matters
        because ``contains_element`` is an ``any`` over cosets, and a repeated
        coset would make the walk quadratic while hiding a wrong representative.
        """
        identity = RigidMotion.identity(3)
        for name in _r1_MEMBERS:
            realization = _r1_group(name).realization
            representatives = realization.representatives
            assert representatives, name
            assert representatives[0].approximately_equals(identity, atol=_ELEMENT_ATOL), (
                f"{name}: representative 0 is not the identity"
            )
            component = realization.component
            for i, r in enumerate(representatives):
                for j, s in enumerate(representatives):
                    if i == j:
                        continue
                    assert not component.contains_element(r.inverse() @ s), (
                        f"{name}: representatives {i} and {j} share a coset"
                    )

    @pytest.mark.foundation
    def test_a_FINITE_realization_is_a_group_and_a_continuous_one_refuses(self) -> None:
        r"""For a finite member the representatives ARE the elements: closed
        under composition and inverse, containing each of its own elements, and
        of the textbook order.  For a continuous member ``elements`` refuses
        rather than returning the representatives — which would be a coset list
        wearing an element list's name.
        """
        for name in _r1_MEMBERS:
            realization = _r1_group(name).realization
            if not realization.is_finite:
                with pytest.raises(ValueError, match="continuous"):
                    _ = realization.elements
                continue
            elements = realization.elements
            assert len(elements) == _r1_ORDERS[name], (
                f"{name}: |G| = {len(elements)}, textbook {_r1_ORDERS[name]}"
            )
            matrices = np.array([e.linear for e in elements])
            products = np.einsum("aij,bjk->abik", matrices, matrices).reshape(-1, 3, 3)
            assert _r1_is_member(name, products).all(), f"{name} is not closed"
            inverses = matrices.transpose(0, 2, 1)      # orthogonal
            assert _r1_is_member(name, inverses).all(), f"{name} lacks an inverse"
            for e in elements:
                assert realization.contains_element(e), f"{name}: an element is not its own"

    @pytest.mark.foundation
    def test_the_component_axis_exists_EXACTLY_for_the_torus(self) -> None:
        r"""``axis`` is defined for ``dim == 1`` and refused otherwise — a
        ``{e}`` or an :math:`SO(3)` has no distinguished direction, and a
        property that answered anyway would hand a fiction to every consumer
        that keys on it.  For the six axial members the axis is the coordinate
        axis the tag NAMES, and it is the kernel of the generator.
        """
        for name in _r1_MEMBERS:
            component = _r1_group(name).realization.component
            if component.dim != 1:
                with pytest.raises(ValueError, match="axis"):
                    _ = component.axis
                continue
            expected = _r1_axis("z" if name == "Dinfh" else name[-1])
            axis = component.axis
            np.testing.assert_allclose(np.abs(axis), np.abs(expected), atol=1e-12,
                                       err_msg=f"{name}: axis {axis}")
            (X,) = component.generators
            np.testing.assert_allclose(X @ axis, np.zeros(3), atol=1e-12,
                                       err_msg=f"{name}: the axis is not ker X")


# --------------------------------------------------------------------------
# G2 — the finite realizations ARE the textbook groups
# --------------------------------------------------------------------------


class TestR1TheFiniteRealizationsAreTheTextbookGroups:
    r"""The keystone reference: an element set built from each group's
    DEFINITION, in plain numpy, agrees with production's realization.

    Production closes a GENERATING set; this closes nothing — every element is
    written down.  ``I_h`` is the sharpest row: production runs a BFS closure of
    :math:`\{R_5, R_3, -I\}`, the reference solves 60 little linear systems for
    the flag maps.  Nothing but the shared standard SETTING is common to the two.
    """

    @pytest.mark.foundation
    def test_the_reference_construction_is_ITSELF_a_group(self) -> None:
        """The reference's validity control (``vv-principles`` #17): the sets
        below are orthogonal, closed, of the textbook order, and contain the
        identity — measured WITHOUT touching production, so a later row's
        agreement is evidence about the SUT rather than about my arithmetic."""
        reference = _r1_finite_reference()
        assert set(reference) == set(_r1_ORDERS), sorted(set(reference) ^ set(_r1_ORDERS))
        for name, M in reference.items():
            assert len(M) == _r1_ORDERS[name], (name, len(M), _r1_ORDERS[name])
            np.testing.assert_allclose(
                M @ M.transpose(0, 2, 1), np.broadcast_to(np.eye(3), M.shape),
                atol=1e-14, err_msg=f"{name}: not orthogonal",
            )
            assert (np.abs(M - np.eye(3)).max(axis=(1, 2)) <= 1e-14).sum() == 1, name
            products = np.einsum("aij,bjk->abik", M, M).reshape(-1, 3, 3)
            distance = np.abs(products[:, None] - M[None]).max(axis=(2, 3))
            assert (distance.min(axis=1) <= 1e-12).all(), f"{name}: not closed"

    @pytest.mark.foundation
    def test_production_realizes_exactly_the_independently_built_element_set(self) -> None:
        r"""22 of 22 finite members, set equality both ways.

        `[M]` 2026-09-02: the worst per-element residual over all 22 members is
        **1.166e-15**, i.e. ``_ELEMENT_ATOL`` (1e-9) is slack by six orders —
        the band absorbs no disagreement here, it only tolerates round-off.  The
        assertion is at 1e-9 rather than at the measured 1.2e-15 because the
        module's own contract is ``_ELEMENT_ATOL`` and a gate may not assert
        more than the type promises (``vv-principles`` #16); the measured
        headroom is the sentence above.
        """
        reference = _r1_finite_reference()
        worst = 0.0
        for name, expected in reference.items():
            realization = _r1_group(name).realization
            assert realization.is_finite, name
            realized = np.array([e.linear for e in realization.elements])
            assert len(realized) == len(expected), (name, len(realized), len(expected))
            distance = np.abs(realized[:, None] - expected[None]).max(axis=(2, 3))
            forward, backward = distance.min(axis=1), distance.min(axis=0)
            assert forward.max() <= _ELEMENT_ATOL, (
                f"{name}: a realized element is {forward.max():.3e} from every "
                f"reference element"
            )
            assert backward.max() <= _ELEMENT_ATOL, (
                f"{name}: a reference element is {backward.max():.3e} from every "
                f"realized element"
            )
            worst = max(worst, float(forward.max()))
        assert worst <= 1e-12, f"worst residual {worst:.3e} — investigate before relaxing"


# --------------------------------------------------------------------------
# G3 — containment
# --------------------------------------------------------------------------


class TestR1ContainmentIsComputedFromTheRealizations:
    r"""``contains`` over the whole roster against an independent membership
    predicate — the row that replaces the retired ``_NAMED_LATTICE``.

    The table it replaced was a hand-written claim with no construction behind
    it, and it had shipped two false edges.  This gate is what makes "computed"
    checkable: neither side can be wrong alone.
    """

    @pytest.mark.foundation
    def test_the_whole_table_agrees_with_an_independent_membership_predicate(self) -> None:
        r"""`[M]` 2026-09-02 — **961 of 961** ordered pairs over the 31-member
        roster.  Both truth values are populated and the census is asserted, so
        a predicate that collapsed to a constant cannot pass.

        Independence, per argument kind:

        * finite ⊆ finite — exact matrix-set containment against the §G2 sets;
        * continuous ⊆ finite — ``False`` by the theorem that a finite group
          holds no infinite one, never by a table;
        * anything ⊆ continuous — the analytic membership predicate
          (:func:`_r1_is_member`), applied to every element / to an
          incommensurate sample.  A sample can only refute, and incommensurate
          angles are what stop it certifying by generating a finite subgroup.
        """
        disagreements = []
        table = {}
        for outer in _r1_MEMBERS:
            for inner in _r1_MEMBERS:
                got = _r1_group(outer).contains(_r1_group(inner))
                expected = _r1_contains(outer, inner)
                table[(outer, inner)] = got
                if got != expected:
                    disagreements.append(f"{outer} ⊇ {inner}: code={got} reference={expected}")
        assert not disagreements, (
            f"{len(disagreements)} of {len(table)} pairs disagree:\n  "
            + "\n  ".join(disagreements[:12])
        )
        trues = sum(table.values())
        assert len(table) == 961, len(table)
        assert 100 < trues < 861, f"only {trues} True of 961 — the table has collapsed"

    @pytest.mark.foundation
    def test_C_1_answers_EXACTLY_as_Trivial_in_every_direction(self) -> None:
        """R1 merges the two spellings of :math:`\\{e\\}`, so the ``C_1`` row and
        column of the containment table must be bit-identical to ``Trivial``'s.

        The roster keeps ``C_1`` precisely so this is assertable; a roster that
        dropped it would make the merge invisible.  Before R1 the two were
        distinct values that contained each other — one group, two answers to
        ``==`` — and a committed test pinned the wrong half of that.
        """
        for other in _r1_MEMBERS:
            g = _r1_group(other)
            assert _r1_group("C_1").contains(g) is _r1_group("Trivial").contains(g), other
            assert g.contains(_r1_group("C_1")) is g.contains(_r1_group("Trivial")), other


# --------------------------------------------------------------------------
# G4 — the normaliser, exact, one body
# --------------------------------------------------------------------------


class TestR1TheNormaliserIsExactAndComputed:
    r"""``is_normalised_by`` and ``normalises`` against a CONJUGATION reference.

    Production decides a continuous member through :math:`\mathrm{Ad}_g\,
    \mathfrak g = \mathfrak g` and the bracket criterion — statements about the
    Lie ALGEBRA.  The reference never touches an algebra: it conjugates elements
    and asks membership.  Two structurally different routes to one answer.
    """

    @pytest.mark.foundation
    def test_is_normalised_by_agrees_with_a_conjugation_reference(self) -> None:
        r"""`[M]` 2026-09-02 — **1426 of 1426** over 31 members × 46 stratified
        motions, **578 True / 848 False**.

        ⚠ Four columns are CONSTANT (all-True) and that is a theorem, not a
        coverage gap (``vv-principles`` #20): ``Trivial``, ``C_1``, :math:`SO(3)`
        and :math:`O(3)` are NORMAL in :math:`O(3)`, so no motion can separate
        them.  The gate asserts that census so a future criterion that made one
        of them separable — or made a separable member constant — is caught.

        The denominator is 46 motions rather than the behaviour contract's 336
        because production costs `[M]` ~78 ms per motion across the roster
        (17.8 s for the full set); the stratified set is chosen to SEPARATE
        (mirrors, half-turns, quarter-turns, incommensurate rotations, improper
        partners, diagonal and polyhedral axes), not to enumerate.  The full
        336-motion table is the main agent's ``scratch/_r1_behaviour.py``
        contract, run once per carve.
        """
        motions = _r1_stratified_motions()
        assert len(motions) == 46, len(motions)
        wrapped = [_r1_motion(m) for m in motions]
        constant, disagreements, trues = [], [], 0
        for name in _r1_MEMBERS:
            expected = _r1_normalised_by(name, motions)
            group = _r1_group(name)
            got = np.array([group.is_normalised_by(m) for m in wrapped])
            trues += int(got.sum())
            if (got != expected).any():
                where = np.flatnonzero(got != expected)[:4]
                disagreements.append(f"{name}: motions {list(where)} disagree")
            if got.all() or not got.any():
                constant.append(name)
        assert not disagreements, "\n  ".join(disagreements[:12])
        assert sorted(constant) == ["C_1", "O3", "SO3", "Trivial"], constant
        assert trues == 578, f"True count moved: {trues} (was 578 on 2026-09-02)"

    @pytest.mark.foundation
    def test_normalises_agrees_with_the_elementwise_reference(self) -> None:
        r"""`[M]` **931 of 961** ordered pairs; the 30 excluded are named.

        EXCLUDED, with the reason: ``normalises(self, other)`` asks
        ``other.is_normalised_by`` once per element of ``self``, and
        ``is_normalised_by`` on :math:`O_h` / :math:`I_h` is
        :math:`O(|G|^2)` coset arithmetic — `[M]` the 44 pairs with those two as
        ``other`` cost 19.9 s of a 23.0 s full table.  Sixteen of the 44 are
        kept (every continuous member, the mirrors, ``C_2``, ``C_4``, ``D_2h``
        and ``Trivial``), which covers both the exact-component leg and the
        representative leg; the 30 dropped rows are finite selves whose only
        distinctive content is a longer element list.
        """
        cheap = [m for m in _r1_MEMBERS if m not in ("Oh", "Ih")]
        polyhedral_selves = (["Trivial", "Dinfh", "SO3", "O3"]
                             + [f"sigma_{a}" for a in _r1_AXES]
                             + [f"SO2_{a}" for a in _r1_AXES]
                             + [f"O2_{a}" for a in _r1_AXES]
                             + ["C_2", "C_4", "D_2h"])
        pairs = ([(a, b) for a in _r1_MEMBERS for b in cheap]
                 + [(a, b) for a in polyhedral_selves for b in ("Oh", "Ih")])
        assert len(pairs) == 931, len(pairs)
        disagreements, trues = [], 0
        for a, b in pairs:
            expected = bool(_r1_normalised_by(b, _r1_sample(a)).all())
            got = _r1_group(a).normalises(_r1_group(b))
            trues += got
            if got != expected:
                disagreements.append(f"{a}.normalises({b}): code={got} reference={expected}")
        assert not disagreements, (
            f"{len(disagreements)} of {len(pairs)}:\n  " + "\n  ".join(disagreements[:12])
        )
        assert 100 < trues < len(pairs) - 100, f"{trues} True of {len(pairs)} — collapsed"

    @pytest.mark.foundation
    def test_the_component_criterion_answers_every_CASE_of_the_theorem(self) -> None:
        r"""The five cases of :math:`\mathrm{Lie}\,N(H) = \{X : [X,\mathfrak h]
        \subseteq \mathfrak h,\; X - \mathrm{Ad}_s X \in \mathfrak h\}`, each
        with a positive AND a negative leg (``vv-principles`` #11).

        The retired ``_identity_component_normalises`` had five hand-written
        arms; one of them was reached ONCE in 670 tests.  These are the same
        five answers demanded of the one body, so a future arm-collapse cannot
        quietly lose a case.

        * **finite other** — :math:`SO(2)_a` normalises a finite :math:`H` iff
          every element commutes with :math:`[\hat a]_\times` (:math:`\mathfrak
          h = 0`, so only the second clause bites).  `[M]` that separates
          :math:`C_n` from :math:`D_{nh}` about the SAME axis: the cyclic
          rotations commute, the in-plane :math:`C_2'` axes and the vertical
          mirrors are rotated to axes the group does not contain.
        * **axial other, SAME axis** — the bracket vanishes and
          :math:`X - \mathrm{Ad}_{\sigma_v}X = 2X \in \mathfrak h`.
        * **axial other, DIFFERENT axis** — :math:`[[\hat a]_\times,[\hat
          b]_\times] = [\hat a\times\hat b]_\times \notin \mathbb R[\hat
          b]_\times`.
        * **D∞h other** — the same, with the axis pinned to :math:`z`.
        * **component SO(3)** — :math:`\mathrm{Ad}` is transitive on the
          non-central elements, so only :math:`H \subseteq \{e, -I\}` survives.
        """
        G = SubgroupOfO3
        # (1) finite other: C_n and D_nh are realized about z, and the criterion
        #     SEPARATES them — this pair is the case's discriminating content.
        for n in (2, 3, 4, 6):
            assert G.SO2("z").normalises(G.Cn(n)), n
            assert not G.SO2("z").normalises(G.Dnh(n)), n
            assert not G.SO2("x").normalises(G.Cn(n)), n
        assert not G.SO2("z").normalises(G.Dnh(1))
        #     the mirror table is the diagonal: sigma_b commutes with [a]_x iff b == a
        for a in _r1_AXES:
            for b in _r1_AXES:
                assert G.SO2(a).normalises(G.Mirror(b)) is (a == b), (a, b)
        # (2)/(3) axial other, same and different axis.
        for a in _r1_AXES:
            assert G.SO2(a).normalises(G.O2(a)) and G.SO2(a).normalises(G.SO2(a))
            assert G.O2(a).normalises(G.SO2(a)) and G.O2(a).normalises(G.O2(a))
            for b in _r1_AXES:
                if b == a:
                    continue
                assert not G.SO2(a).normalises(G.O2(b)), (a, b)
                assert not G.SO2(a).normalises(G.SO2(b)), (a, b)
        # (4) D_inf_h other — pinned to z.
        assert G.SO2("z").normalises(G.Dinfh) and G.Dinfh.normalises(G.Dinfh)
        assert not G.SO2("x").normalises(G.Dinfh)
        assert not G.SO2("y").normalises(G.Dinfh)
        # (5) component SO(3): only the centre {e, -I} survives.
        assert G.SO3.normalises(G.Trivial) and G.SO3.normalises(G.SO3)
        assert G.SO3.normalises(G.O3) and G.O3.normalises(G.SO3)
        for finite in (G.Mirror("z"), G.Cn(2), G.Dnh(2), G.OctahedralOh, G.IcosahedralIh):
            assert not G.SO3.normalises(finite), finite.name

    @pytest.mark.foundation
    def test_CONTAINMENT_and_NORMALISATION_are_independent_in_BOTH_directions(self) -> None:
        r"""The row that forces step 1 (normalise) to run BEFORE step 2
        (contains) in the invariance kernel — and it separates the two
        predicates in **both** directions on one group and one family.

        `[M]` 2026-09-02, :math:`O(2)_x` against the three coordinate mirrors:

        =========  ===========  ============
        mirror     contained?   normalised?
        =========  ===========  ============
        σ_x        **no**       **yes**
        σ_y        **yes**      **no**
        σ_z        **yes**      **no**
        =========  ===========  ============

        Reading it: :math:`\sigma_x` has NORMAL :math:`\hat e_x`, so it flips
        the axis and is not in the stabiliser — yet it commutes with every
        rotation about :math:`x`, so it normalises.  :math:`\sigma_y` is IN
        :math:`O(2)_x` and a rotation about :math:`x` carries it to a mirror in
        a tilted plane, which is not one of the two elements of
        :math:`\{e, \sigma_y\}`.  Neither predicate implies the other, so a
        kernel that asked ``contains`` first would answer a question about a
        group that does not act.
        """
        O2x = SubgroupOfO3.O2("x")
        expected = {"x": (False, True), "y": (True, False), "z": (True, False)}
        for axis, (contained, normalised) in expected.items():
            mirror = SubgroupOfO3.Mirror(axis)
            assert O2x.contains(mirror) is contained, axis
            assert O2x.normalises(mirror) is normalised, axis


# --------------------------------------------------------------------------
# G5 — IdentityComponent.fixes: the ERR-072 exactness witnesses, re-keyed
# --------------------------------------------------------------------------


def _r1_sphere_measure(q) -> DiscreteMeasure:
    return DiscreteMeasure(
        nodes=np.column_stack([q.mu_x, q.mu_y, q.mu_z]), weights=q.weights, support=SPHERE,
    )


class TestR1AConnectedGroupIsDecidedByItsLieAlgebra:
    r"""``IdentityComponent.fixes`` — the exact criterion that replaced the
    sampled one, re-keyed off the retired ``_is_axis_supported`` /
    ``_is_origin_supported`` / ``_fixes_every_point`` triple.

    A connected group's orbits are connected, and a connected orbit inside a
    finite set is a point; so :math:`G^0` fixes a finite set iff :math:`Xp = 0`
    for every generator and every point.  The retired sample tested closure
    under :math:`C_4` and called it :math:`SO(2)` — unsound in the DANGEROUS
    direction, certifying non-invariant rules as a function of ``n_phi mod 4``
    (ERR-072).
    """

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-072")
    @pytest.mark.parametrize("axis", _r1_AXES)
    def test_a_torus_fixes_exactly_the_points_ON_its_axis(self, axis: str) -> None:
        r""":math:`[\hat a]_\times p = \hat a \times p`, whose norm IS the point's
        distance from the axis — so ``fixes`` is the axis-support criterion,
        stated once, with no branch on the axis letter.

        Positive leg: points on the axis, at any signed coordinate, including
        an ASYMMETRIC set (the criterion must NOT secretly require a mirror).
        Negative leg: one point off the axis by :math:`\rho`, at
        :math:`\rho = 10^{-6}` down to :math:`\rho = 1` — the refusal is on the
        DISTANCE, so the smallest one that clears the window must already fail.
        """
        component = SubgroupOfO3.SO2(axis).realization.component
        e = _r1_axis(axis)
        on_axis = np.array([0.9 * e, -0.3 * e, np.zeros(3), 1.0 * e])
        assert component.fixes(on_axis, atol=1e-13)
        other = _r1_axis(_r1_AXES[( "xyz".index(axis) + 1) % 3])
        for rho in (1e-6, 1e-3, 0.5, 1.0):
            off_axis = np.vstack([on_axis, (0.2 * e + rho * other)[None, :]])
            assert not component.fixes(off_axis, atol=1e-13), rho
        # the O(2)_a component is the SAME torus — one criterion, two members
        assert SubgroupOfO3.O2(axis).realization.component.fixes(on_axis, atol=1e-13)

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-072")
    def test_no_shipped_azimuthal_rule_is_axis_supported_at_ANY_n_phi(self) -> None:
        r"""The ERR-072 signature, in the regime that produced it: the retired
        check answered as a function of ``n_phi mod 4`` — ``True`` for 4/8/12/16
        and ``False`` for 2/3/5/6/7 — because four right angles generate
        :math:`C_4`.  ``fixes`` answers ``False`` at every ``n_phi``, because a
        product rule's nodes are off every coordinate axis.

        The ladder deliberately BREAKS its own arithmetic pattern
        (``vv-principles`` #13): two residues mod 4, an odd, a prime.
        """
        for n_phi in (3, 4, 5, 6, 8, 12, 16):
            nodes = np.column_stack([
                (q := Quadrature.product(4, n_phi)).mu_x, q.mu_y, q.mu_z,
            ])
            for axis in _r1_AXES:
                component = SubgroupOfO3.SO2(axis).realization.component
                assert not component.fixes(nodes, atol=1e-13), (n_phi, axis)
                assert not _r1_sphere_measure(q).is_invariant_under(SubgroupOfO3.SO2(axis)), (
                    n_phi, axis,
                )

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-072")
    def test_the_SO3_component_fixes_only_the_ORIGIN(self) -> None:
        r"""Every non-origin :math:`SO(3)`-orbit is a whole 2-sphere, so a finite
        invariant set can hold only the origin.  The icosahedral vertex set is
        the witness the retired check got WRONG: it sampled :math:`I_h` and
        called the answer :math:`SO(3)`, and 60 of those 120 matrices are
        improper, hence not in :math:`SO(3)` at all.
        """
        component = SubgroupOfO3.SO3.realization.component
        assert component.dim == 3
        assert component.fixes(np.zeros((1, 3)), atol=1e-13)
        assert component.fixes(np.zeros((4, 3)), atol=1e-13)
        assert not component.fixes(_r1_icosahedron(), atol=1e-13)
        assert not component.fixes(np.array([[1e-6, 0.0, 0.0]]), atol=1e-13)
        assert SubgroupOfO3.O3.realization.component.dim == 3

    @pytest.mark.foundation
    @pytest.mark.catches("ERR-072")
    def test_the_SAMPLED_criterion_would_certify_what_fixes_REFUSES(self) -> None:
        r"""The bug, re-introduced HERE as a reference rather than in production
        — the honest way to show the exact criterion is not merely
        self-consistent.

        ``product(4, 8)`` is closed under the four right-angle rotations about
        :math:`z`, so the retired sample certified it :math:`SO(2)_z`-invariant.
        `[M]` the sampled predicate says **True**, ``fixes`` says **False**, and
        the true azimuthal group is the finite :math:`D_{8h}` — which the
        module does report.
        """
        q = Quadrature.product(4, 8)
        nodes = np.column_stack([q.mu_x, q.mu_y, q.mu_z])

        def sampled_so2_z_invariant(points: np.ndarray) -> bool:
            """The retired {0, 90, 180, 270} degree check — C_4, not SO(2)."""
            for k in range(4):
                image = points @ _r1_rot_z(np.pi * k / 2.0).T
                distance = np.abs(image[:, None] - points[None]).max(axis=2)
                if not (distance.min(axis=1) <= 1e-9).all():
                    return False
            return True

        assert sampled_so2_z_invariant(nodes) is True
        assert SubgroupOfO3.SO2("z").realization.component.fixes(nodes, atol=1e-13) is False
        assert _r1_sphere_measure(q).is_invariant_under(SubgroupOfO3.SO2("z")) is False
        assert _r1_sphere_measure(q).is_invariant_under(SubgroupOfO3.Dnh(8)) is True


# --------------------------------------------------------------------------
# G6/G7 — one group, one value
# --------------------------------------------------------------------------


class TestR1TheTrivialGroupHasExactlyOneSpelling:
    r"""``Cn(1)`` IS :math:`\{e\}`, so it must be the same VALUE as ``Trivial``.

    Before R1 they were two values that contained each other, and the door they
    hit was the orbit-space catalogue: ``SPHERE.quotient(Cn(1))`` reported "no
    catalogue entry" — the one message the door promises never to give, for the
    group every catalogue has.
    """

    @pytest.mark.foundation
    def test_both_spellings_are_ONE_value_through_every_constructor(self) -> None:
        r"""``==``, ``hash``, the container, ``name`` and ``repr`` — through the
        classmethod AND through the bare constructor, because
        ``__post_init__`` is what normalises and a factory-only fix would leave
        ``SubgroupOfO3(Cn(1))`` spelling the old value.

        ⚠ **``is`` is deliberately NOT asserted.** `[M]` 2026-09-02
        ``SubgroupOfO3.Cn(1) is SubgroupOfO3.Trivial`` is **False** — the
        constructor returns a fresh instance carrying the normalised tag, not
        the pre-instantiated singleton.  The plan's done-when spells the
        predicate with ``is``; that leg would be a FALSE RED against a correct
        implementation, and the claim that matters is value identity.
        """
        from orpheus.numerics.symmetry import Cn as CnTag

        trivial = SubgroupOfO3.Trivial
        for spelling in (SubgroupOfO3.Cn(1), SubgroupOfO3(CnTag(1))):
            assert spelling == trivial
            assert hash(spelling) == hash(trivial)
            assert spelling.name == trivial.name == "Trivial"
            assert repr(spelling) == repr(trivial)
            assert spelling.dim == 0 and spelling.orbit_stabiliser == trivial
        assert len({SubgroupOfO3.Cn(1), SubgroupOfO3(CnTag(1)), trivial}) == 1
        # and the merge must not swallow a real member
        assert len({SubgroupOfO3.Cn(2), SubgroupOfO3.Cn(1), trivial}) == 2

    @pytest.mark.foundation
    def test_the_orbit_space_door_answers_for_the_merged_value(self) -> None:
        """``S^2/C_1`` is ``S^2/{e}`` — the catalogue lookup must find the entry
        it always had, under either spelling, rather than reporting a missing
        derivation for the group every catalogue contains."""
        assert SPHERE.quotient(SubgroupOfO3.Cn(1)) == SPHERE.quotient(SubgroupOfO3.Trivial)


class TestR1TheGroupIsAFrozenValue:
    """A group is a value: frozen, hashable, compared by its tag.

    Mutability was not theoretical — ``g._tag = …`` succeeded and MOVED
    ``hash(Quotient)``, so a group already used as a dictionary key could be
    edited out from under the dictionary.
    """

    @pytest.mark.foundation
    def test_assignment_to_the_only_field_is_refused(self) -> None:
        """Modelled on ``TestManifold::test_variants_are_frozen_and_compare_by_value``.

        The field is discovered through ``dataclasses.fields`` rather than
        named, so the gate also asserts the class IS a dataclass with exactly
        one field — which is the structural half of "the tag is the identity".
        """
        assert dataclasses.is_dataclass(SubgroupOfO3)
        fields = dataclasses.fields(SubgroupOfO3)
        assert len(fields) == 1, [f.name for f in fields]
        group = SubgroupOfO3.Cn(4)
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(group, fields[0].name, None)

    @pytest.mark.foundation
    def test_distinct_members_are_distinct_KEYS(self) -> None:
        r"""Separation is asserted through the CONTAINER, never through
        ``hash(a) != hash(b)``.

        ⚠ A frozen dataclass hashes its FIELD TUPLE, so members of different
        tag families can share a hash while ``__eq__`` still separates them
        (``__eq__`` opens with a class check).  `[M]` 2026-09-02
        ``hash(SO2('x')) == hash(O2('x')) == hash(Mirror('x'))``; a
        ``hash !=`` leg reds on CORRECT code and reads as extra rigour.
        """
        members = [_r1_group(name) for name in _r1_MEMBERS]
        # C_1 and Trivial are ONE value by design; every other pair is distinct.
        assert len(set(members)) == len(_r1_MEMBERS) - 1
        assert len({SubgroupOfO3.Mirror(a) for a in _r1_AXES}) == 3
        assert len({SubgroupOfO3.SO2("x"), SubgroupOfO3.O2("x"),
                    SubgroupOfO3.Mirror("x")}) == 3
        lookup = {g: g.name for g in members}
        for g in members:
            assert lookup[g] == g.name


# --------------------------------------------------------------------------
# G8 — identity_component, corrected
# --------------------------------------------------------------------------


class TestR1TheIdentityComponentIsTheConnectedHalf:
    r"""``identity_component`` is :math:`G^0`, the CONNECTED component of the
    identity — which for a discrete group is a point, i.e. :math:`\{e\}`.

    Before R1 the property returned the group ITSELF for every finite member.
    `[M]` that is wrong on **21 of the roster's 22** finite members (all but
    ``Trivial``), and it had **zero** readers tree-wide — a property with no
    consumer and no test, which is exactly how a wrong answer survives.
    """

    @pytest.mark.foundation
    def test_the_named_table_over_the_whole_roster(self) -> None:
        """`[M]` 2026-09-02: every finite member answers ``Trivial``; the axial
        pair answers its own ``SO2_a``; ``D_∞h`` answers ``SO2_z``; ``O(3)``
        answers ``SO(3)``.  31 of 31 rows written out, because a table is the
        only way to state which rows the carve MOVED."""
        for name in _r1_MEMBERS:
            got = _r1_group(name).identity_component
            assert got.name == _r1_IDENTITY_COMPONENT[name], (
                f"{name}: identity_component is {got.name}, "
                f"expected {_r1_IDENTITY_COMPONENT[name]}"
            )

    @pytest.mark.foundation
    def test_it_obeys_the_three_laws_G0_must_obey(self) -> None:
        r"""Independent of the table: :math:`G^0 \subseteq G`,
        :math:`\dim G^0 = \dim G` (the component carries the whole dimension),
        and :math:`(G^0)^0 = G^0` (idempotence — a connected group is its own
        identity component).  A table can be transcribed wrongly; these cannot.
        """
        for name in _r1_MEMBERS:
            group = _r1_group(name)
            component = group.identity_component
            assert group.contains(component), name
            assert component.dim == group.dim, (name, component.dim, group.dim)
            assert component.identity_component == component, name
            assert component.realization.component.dim == group.realization.component.dim


# --------------------------------------------------------------------------
# G9 — the orbit stabiliser is a genuine MAXIMUM
# --------------------------------------------------------------------------


def _r1_orbit_probe_points() -> np.ndarray:
    """A generic direction, the six coordinate poles' half, the face and vertex
    diagonals, and five seeded generic directions — the loci a wrong stabiliser
    would have to agree on."""
    p0 = np.array([0.3708, -0.5921, 0.7154])
    p0 /= np.linalg.norm(p0)
    rng = np.random.default_rng(7)
    generic = [(lambda v: v / np.linalg.norm(v))(rng.normal(size=3)) for _ in range(5)]
    return np.array([p0, _r1_axis("x"), _r1_axis("y"), _r1_axis("z"),
                     np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
                     np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
                     np.array([0.6, 0.8, 0.0]), np.array([0.0, 0.6, 0.8]), *generic])


def _r1_in_orbit(name: str, images: np.ndarray, p: np.ndarray) -> np.ndarray:
    r"""``(m,)`` boolean: is each row of ``images`` in the :math:`H`-orbit of
    ``p``?  Analytic per family — the orbits of the continuous members are the
    latitude circles / circle-pairs / spheres, so membership is two scalar
    comparisons and needs no element list."""
    finite = _r1_finite_reference()
    if name in finite:
        orbit = finite[name] @ p
        return np.abs(images[:, None, :] - orbit[None]).max(axis=2).min(axis=1) <= _r1_ATOL
    same_radius = np.abs(np.linalg.norm(images, axis=1) - np.linalg.norm(p)) <= _r1_ATOL
    if name in ("SO3", "O3"):
        return same_radius
    if name.startswith("SO2_") or name.startswith("O2_"):
        e = _r1_axis(name[-1])
        return same_radius & (np.abs(images @ e - p @ e) <= _r1_ATOL)
    if name == "Dinfh":
        e = _r1_axis("z")
        return same_radius & (np.abs(np.abs(images @ e) - abs(p @ e)) <= _r1_ATOL)
    raise ValueError(name)  # pragma: no cover - roster guard


class TestR1TheOrbitStabiliserIsAGenuineMaximum:
    r"""``orbit_stabiliser`` is the LARGEST subgroup of :math:`O(3)` with this
    group's orbits — a maximum, not a lookup.

    The search is made finite by an argument that costs nothing: if
    :math:`g\,p_0 \in H p_0` then :math:`h^{-1}g \in \mathrm{Stab}(p_0)` for
    some :math:`h \in H`, so the whole candidate set lies in
    :math:`H\cdot O(2)_{p_0}` — enumerable by sampling the one-parameter
    stabiliser of a generic point.
    """

    @staticmethod
    def _stabiliser_of(p0: np.ndarray, n: int = 90) -> np.ndarray:
        r""":math:`O(2)_{p_0}` sampled: ``n`` rotations about :math:`p_0` and
        ``n`` reflections in planes through it."""
        out = [_r1_rot(p0, t) for t in np.linspace(0.0, 2 * np.pi, n, endpoint=False)]
        e1 = np.cross(p0, [0.0, 0.0, 1.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(p0, e1)
        out += [_r1_refl(np.cos(t) * e1 + np.sin(t) * e2)
                for t in np.linspace(0.0, np.pi, n, endpoint=False)]
        return np.array(out)

    @pytest.mark.foundation
    def test_nothing_OUTSIDE_the_stabiliser_shares_the_orbits(self) -> None:
        r"""The MAXIMALITY half — `[M]` 2026-09-02, **31 of 31** members, 0.6 s.

        Every candidate :math:`g \in H\cdot O(2)_{p_0}` (180-point stabiliser
        sample, 13 probe points) that carries every probe point inside its own
        :math:`H`-orbit is shown to lie INSIDE ``orbit_stabiliser(H)``.  So no
        group strictly larger than the reported one can have :math:`H`'s
        orbits, which is the half a "same orbits" check alone cannot give.
        """
        probes = _r1_orbit_probe_points()
        stabiliser = self._stabiliser_of(probes[0])
        for name in _r1_MEMBERS:
            H = _r1_sample(name)
            candidates = np.einsum("hij,sjk->hsik", H, stabiliser).reshape(-1, 3, 3)
            keep = np.ones(len(candidates), bool)
            for p in probes:
                live = np.flatnonzero(keep)
                if live.size == 0:
                    break
                inside = _r1_in_orbit(name, candidates[live] @ p, p)
                keep[live[~inside]] = False
            survivors = candidates[keep]
            reported = _r1_group(name).orbit_stabiliser
            label = next(m for m in _r1_MEMBERS if _r1_group(m) == reported)
            assert _r1_is_member(label, survivors).all(), (
                f"{name}: a motion with {name}'s orbits lies OUTSIDE the "
                f"reported stabiliser {reported.name}"
            )

    @pytest.mark.foundation
    def test_the_stabiliser_really_HAS_the_group_s_orbits(self) -> None:
        r"""The correctness half, and the one that witnesses the two members
        that GROW: every element of ``orbit_stabiliser(H)`` must carry each
        probe point inside its :math:`H`-orbit.

        This is where :math:`SO(2)_a \to O(2)_a` and :math:`SO(3) \to O(3)` are
        earned rather than asserted — the vertical mirrors and the inversion
        are elements of the reported stabiliser that are NOT in :math:`H`, and
        they must still preserve every orbit.
        """
        probes = _r1_orbit_probe_points()
        grew = {}
        for name in _r1_MEMBERS:
            group = _r1_group(name)
            reported = group.orbit_stabiliser
            label = next(m for m in _r1_MEMBERS if _r1_group(m) == reported)
            K = _r1_sample(label)
            for p in probes:
                assert _r1_in_orbit(name, K @ p, p).all(), (
                    f"{name}: {reported.name} moves {p} out of the orbit"
                )
            assert reported.contains(group), name
            assert reported.orbit_stabiliser == reported, name
            if reported != group:
                grew[name] = reported.name
                outside = ~_r1_is_member(name, K)
                assert outside.any(), (
                    f"{name}: {reported.name} is reported LARGER but every "
                    f"sampled element already lies inside {name} — the growth "
                    f"is not witnessed"
                )
        assert grew == {"SO2_x": "O2_x", "SO2_y": "O2_y", "SO2_z": "O2_z",
                        "SO3": "O3"}, grew

    @pytest.mark.foundation
    def test_a_WRONG_stabiliser_claim_fails_the_orbit_check(self) -> None:
        """The control (``vv-principles`` #17): the two halves above are only
        evidence if a wrong claim would fail them.  ``sigma_x`` does NOT have
        ``O(2)_x``'s orbits — a rotation about x moves a generic point off its
        two-element mirror orbit — so claiming it would red the correctness
        half; and ``O(2)_x`` does not have ``O(3)``'s orbits, so claiming a
        bigger stabiliser for it would red too.

        ⚠ **A DECLARED non-catcher.** This row calls no production predicate:
        it asserts a property of the REFERENCE, so no mutation of the module
        can redden it, and `[M]` none of the 20-arm battery does.  Its green is
        not coverage of ``orbit_stabiliser``; it is the licence to read the two
        rows above as coverage.  The same holds for
        ``test_the_reference_construction_is_ITSELF_a_group`` in §G2."""
        probes = _r1_orbit_probe_points()
        for subject, false_claim in (("sigma_x", "O2_x"), ("O2_x", "O3"),
                                     ("C_4", "D_4h"), ("Trivial", "sigma_z")):
            K = _r1_sample(false_claim)
            broken = any(
                not _r1_in_orbit(subject, K @ p, p).all() for p in probes
            )
            assert broken, (
                f"the control is inert: {false_claim} appears to share "
                f"{subject}'s orbits on every probe point"
            )


class TestR1TheRealizationRefusesIllegalStates:
    r"""The guards travel with the theorem (elegance review, 2026-09-03).

    Every predicate on :class:`IdentityComponent` reads ``dim`` as
    ``len(generators)``, where the mathematics means RANK, and
    :math:`\mathfrak{so}(3)` has subalgebras of dimension 0, 1 and 3 only. So
    a dependent, non-skew or two-element basis must be refused where it is
    written: `[M]` before the guard, ``IdentityComponent((X, X, X))`` read as
    :math:`\mathfrak{so}(3)` and admitted every proper motion, and
    ``Realization(IdentityComponent(()), ())`` was a "group" whose
    ``contains_element(identity)`` was ``False``.  The mutation battery's
    two-generator arm is now UNINSTALLABLE — refused at construction — which
    is the strongest verdict a guard can earn, and this class is its direct
    witness (``vv-principles`` #11: one positive leg, the negative legs).
    """

    @pytest.mark.foundation
    def test_a_two_dimensional_or_dependent_or_non_skew_basis_is_refused(self) -> None:
        skew_z = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        skew_x = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
        # positive leg: the three legal dimensions construct
        for generators in ((), (skew_z,), (skew_x, skew_z, skew_x @ skew_z - skew_z @ skew_x)):
            assert IdentityComponent(generators).dim == len(generators)
        with pytest.raises(ValueError, match="dimension 0, 1 or 3, not 2"):
            IdentityComponent((skew_x, skew_z))
        with pytest.raises(ValueError, match="linearly dependent"):
            IdentityComponent((skew_z, skew_z, skew_z))
        with pytest.raises(ValueError, match="skew-symmetric"):
            IdentityComponent((np.eye(3),))

    @pytest.mark.foundation
    def test_a_zero_generator_is_refused_as_dependent_not_as_non_skew(self) -> None:
        with pytest.raises(ValueError, match="linearly dependent"):
            IdentityComponent((np.zeros((3, 3)),))

    @pytest.mark.foundation
    def test_an_empty_or_identity_less_coset_list_is_refused(self) -> None:
        trivial = IdentityComponent(())
        identity = RigidMotion.identity(3)
        sigma_y = RigidMotion.reflection(normal=np.array([0.0, 1.0, 0.0]))
        # positive leg
        assert Realization(trivial, (identity, sigma_y)).contains_element(sigma_y)
        with pytest.raises(ValueError, match="must not be empty"):
            Realization(trivial, ())
        with pytest.raises(ValueError, match="first coset representative must be the identity"):
            Realization(trivial, (sigma_y, identity))
