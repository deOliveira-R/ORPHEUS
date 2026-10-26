r"""Invariance verdicts of a discrete measure under the subgroup lattice of :math:`O(3)`.

Re-filed 2026-10-25 from ``tests/gates/geometry/test_symmetry.py`` (its group 2)
by the layer-DAG campaign (ruling 23): the asserted verbs —
``DiscreteMeasure.is_invariant_under``, ``symmetry_groups``,
``singular_set_under``, ``certificate_under``, ``permutation_under`` — belong
to a numerics object, so the claims are homed in numerics; the lattice they
are asked of stays in geometry (``tests/gates/geometry/test_symmetry.py``),
and the helpers both files read live in
``tests/gates/geometry/_lattice_fixtures.py``.  One row joined from
``tests/gates/geometry/test_manifold.py``
(``TestR4OrbitBarycentresIsOneConceptOnBothWidths``, the harness fold answer).
Every test here is tagged ``foundation``, as it was at its old home.

**Invariance against existing quadratures** — a
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

from orpheus.geometry.transformation import RigidMotion
from orpheus.geometry.manifold import COSINE_INTERVAL, SPHERE, Ball, Quotient
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_polar_orbit
from orpheus.geometry.symmetry import SubgroupOfO3
from orpheus.numerics.measure import WEIGHT_ATOL, _node_window
from orpheus.geometry.symmetry import candidate_groups

from tests.gates.geometry._lattice_fixtures import _g22b_AXES, _every_spellable_group, _R4_I3


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


def _polar_x(nodes, weights) -> DiscreteMeasure:
    r"""A :math:`\mu`-rule on the polar orbit space :math:`S^2/O(2)_x`.

    Re-posed at S5 (ruling 44, 2026-10-19): until then these fixtures lived on the
    bare chart ``[-1, 1]``, which the old kernel zero-padded to :math:`(\mu, 0, 0)` and
    asked O(3) questions of. The chart carries only the 1-D motions, so an O(3)
    question there is now a ``TypeError``; the orbit space is where the slab's rule
    lives in production (tracker 2.4), and where the same questions have answers.
    """
    return DiscreteMeasure(
        nodes=np.asarray(nodes, dtype=float), weights=np.asarray(weights, dtype=float),
        support=COSINE_INTERVAL,
    ).on_orbit_space(SPHERE.quotient(SubgroupOfO3.O2("x")))


def _measure_from_1d_quad(q) -> DiscreteMeasure:
    """A 1-D quadrature's rule on :math:`S^2/O(2)_x` (see :func:`_polar_x`)."""
    return _polar_x(q.mu_x, q.weights)


def _candidates_of(measure: DiscreteMeasure):
    """The candidate groups of a measure: ``symmetry.candidate_groups`` on its ambient
    points, at the old default position window (S5: until then the free
    ``invariance.candidate_groups(measure)``)."""
    return candidate_groups(
        measure.support.ambient_points(measure.nodes), atol=_node_window(WEIGHT_ATOL)
    )


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

    The slab's measure (support ``S^2/O2_x``, embedded along the axis its
    orbit space names) answers this way (until S5 a chart-level twin on the
    bare interval, zero-padded along column 0, answered beside it; ruling 44
    retired the padding, and the twin re-posed onto ``S^2/O2_x`` IS this
    measure, so the loop keeps one arm) — and a marginal declared about :math:`z` answers the
    opposite way, which is what makes the axis load-bearing rather than
    decorative.

    The old "cannot fail" defect is still guarded one section down:
    ``test_invariance_is_exact_on_polar_marginals_not_vacuous`` shows the
    same fixtures ARE distinguishable, by :math:`\sigma_x`.
    """
    from orpheus.numerics.quadrature import gauss_legendre_on_polar_orbit

    q = Quadrature.gauss_legendre(n_ordinates=n)
    for mu in (q.measure,):  # S5: the chart arm re-posed onto S^2/O2_x IS q.measure
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

    Constructed here rather than from :mod:`~orpheus.geometry.symmetry`'s own
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
# Invariance — the D_nh product rule's vertical mirror planes (its section at the
# old home, "Containment computed from the realized operator sets", stays in geometry)
# ============================================================================


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

    violations: list[str] = []
    checked = 0
    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        cands = _candidates_of(mu)
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
        measure = spec.build(params).measure
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
    from orpheus.geometry.symmetry import _close_group, _realized_ops

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
    asymmetric = _polar_x(np.array([-0.9, -0.1, 0.3]), np.ones(3))
    symmetric = _polar_x(np.array([-0.6, -0.2, 0.2, 0.6]), np.ones(4))

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
        "asymmetric": _polar_x(np.array([-0.9, -0.1, 0.3]), np.ones(3)),
        "symmetric": _polar_x(np.array([-0.6, -0.2, 0.2, 0.6]), np.ones(4)),
        "symmetric-nodes-unequal-weights": _polar_x(np.array([-0.6, 0.6]), np.array([1.0, 2.0])),
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

    q = Quadrature.gauss_legendre(n_ordinates=8)
    offered = {g.name for g in _candidates_of(q.measure)}
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
        g = next(c for c in _candidates_of(q.measure) if c.name == absorbed)
        assert q.measure.is_invariant_under(g) and SubgroupOfO3.O2("x").contains(g)
    assert not SubgroupOfO3.O2("x").contains(SubgroupOfO3.Mirror("x"))


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

    Re-posed at S5 (ruling 44): there is one arm now, the composition. The
    "1-D" fixture lives on ``S^2/O2_x`` (the chart refuses an O(3) question) and
    the embedded one on ``Ball(3)``, where its barycentres lie; the claim that
    survives is that the orbit-space verdict equals the barycentres' verdict.

    The asymmetric leg is the discriminating one; the symmetric leg is
    the control that keeps the agreement from being vacuous (on a
    symmetric set every plane answers True and the two arms would agree
    even while asking different questions — which is precisely how this
    went unnoticed).
    """
    weights = np.full(4, 0.5)

    def both_arms(mu: np.ndarray, axis: str) -> tuple[bool, bool]:
        one_d = _polar_x(mu, weights)
        # S5 (ruling 44): (mu, 0, 0) is the orbit's BARYCENTRE, a point of the ball,
        # not of the sphere; until S5 this fixture declared SPHERE for interior
        # points (ERR-080's shape). The surviving claim: the orbit-space verdict
        # equals the verdict on the barycentres in the ball.
        embedded = DiscreteMeasure(
            nodes=np.column_stack([mu, np.zeros(4), np.zeros(4)]),
            weights=weights,
            support=Ball(3),
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
    from orpheus.geometry.manifold import SPHERE

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

    checked = 0
    touching_o2 = 0
    o2_outer = 0
    o2_inner = 0
    violations: list[str] = []
    for name, q in _walk_rules():
        mu = _measure_from_sphere_quad(q)
        cands = _candidates_of(mu)
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


# =============================================================================
# 2.2b — the Γ-slot: gates drafted by the test-architect (2026-09-02)
# =============================================================================


def _g22b_trivially_quotiented_product() -> DiscreteMeasure:
    from dataclasses import replace

    return replace(
        Quadrature.product(4, 8).measure,
        support=SPHERE.quotient(SubgroupOfO3.Trivial),
    )


# =============================================================================


# G6 — tests/gates/numerics/test_invariance_verdicts.py
# =============================================================================


class TestTheTrivialQuotientAnswersExactlyLikeTheBareSphere:
    r"""(v) The bit-identity witness that the orbit-space route CONTAINS the
    ambient kernel.

    ``S²/{e} = S²`` is a theorem (``_mod_trivial``: ``P = I``), so routing a
    trivially-quotiented measure through the kernel (``DiscreteMeasure.is_invariant_under``,
    a composition since S5; ``invariance.is_invariant_under`` from R2 of #434;
    ``_invariance_on_orbit_space`` until then) must
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


# G7 — tests/gates/numerics/test_invariance_verdicts.py
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
        got = {g.name: m.is_invariant_under(g) for g in _candidates_of(m)}
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
        got = {g.name: m.is_invariant_under(g) for g in _candidates_of(m)}
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


# G8 — tests/gates/numerics/test_invariance_verdicts.py
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
        got = {g.name: m.is_invariant_under(g) for g in _candidates_of(m)}
        assert got == self._FOLD_AFTER

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 4), (4, 8), (6, 8), (8, 8)])
    def test_all_four_flips_hold_on_every_shipped_fold(self, n_mu, n_phi):
        m = Quadrature.folded_product(n_mu, n_phi).measure
        for name in ("sigma_y", "C_2", "D_1h", "D_2h"):
            g = next(c for c in _candidates_of(m) if c.name == name)
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


# G9 — tests/gates/numerics/test_invariance_verdicts.py
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
        cands = _candidates_of(measure)
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


# G10 — tests/gates/numerics/test_invariance_verdicts.py
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
        names = sorted(g.name for g in _candidates_of(m))
        assert names == sorted(
            [
                # C_1 IS Trivial since 2026-09-02 (R1 of #434): one group, one spelling.
                "C_2", "D_1h", "D_2h", "Dinfh", "Ih",
                "O2_x", "O2_y", "O2_z", "O3", "Oh", "SO2_x", "SO2_y", "SO2_z",
                "SO3", "Trivial", "sigma_x", "sigma_y", "sigma_z",
            ]
        )


# =============================================================================


# G11 — tests/gates/numerics/test_invariance_verdicts.py
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
        CHART, through the orbit space's action (``induced_action`` until S5,
        ``Quotient.action`` since), composed with ``relabelling_to`` since S5
        (``invariance._orbit_closure`` before it) — so there is ONE invariance
        notion, which is R1.

        What is owed is therefore the AGREEMENT, gated: on a fold, for every
        candidate group, ``is_invariant`` and ``certificate_under is not
        None`` must answer the same — except for the continuous groups, which
        have no finite realization to permute anything with.
        """
        from orpheus.geometry.symmetry import _realized_ops

        m = Quadrature.folded_product(4, 8).measure
        checked = 0
        for g in _candidates_of(m):
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
        ``tests/gates/numerics/test_measure.py::test_the_fold_consumes_the_symmetry_idempotent_only_on_a_trivial_action``
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

        # S5: the free function is retired; the claim lives on the method now.
        doc = DiscreteMeasure.certificate_under.__doc__ or ""
        assert "orbit space" in doc or "chart" in doc, (
            "certificate_under's docstring still describes the AMBIENT "
            "reading; tracker 2.2b moved its body onto the CHART (through "
            "Quotient.action, induced_action until S5), and added a third cause of None (the "
            "group does not normalise the quotienting group). Its two "
            "downstream readers -- `singular_set_under` and "
            "`DiscreteMeasure.quotient` -- inherit the new semantics."
        )


# =============================================================================
# Re-filed from tests/gates/geometry/test_manifold.py (layer DAG ruling 23, 2026-10-25)
# =============================================================================


class TestR4OrbitBarycentresIsOneConceptOnBothWidths:
    """The numerics half of geometry's ``TestR4OrbitBarycentresIsOneConceptOnBothWidths``.

    Re-filed 2026-10-25 from ``tests/gates/geometry/test_manifold.py`` (layer-DAG
    ruling 23): the row below asserts a numerics object's verb, so it lives here; the class keeps its geometry name so the node id
    changes only in its path.  The manifold rows of the class stay in geometry.
    """

    @pytest.mark.foundation
    def test_the_HARNESS_fold_answer_now_agrees_with_the_kernel(self) -> None:
        r"""⭐ ``tests/_harness/references.py``'s partner map is derived from
        ``_embedded_nodes`` (``GSpace.ambient_points`` since S5), so R4 moves ONE of its 33 (rule × axis) answers per
        fold — and it moves it INTO agreement with the kernel.

        ``[M]`` 2026-09-03: 31 of 33 answers unchanged; the 2 that move are
        ``axis="y"`` on the two folds, where the harness RAISED (residual
        ``1.19e+00``) and now returns the identity permutation (residual
        ``0.00e+00``) — which is exactly what ``permutation_under`` has
        answered since #429 tracker 2.2b.  ``[M]`` **no call site passes
        ``axis="y"`` on a fold** (both pass ``"x"``: a literal at
        ``test_snmesh_realizer_wiring.py:439,450``, and the CYL fold's only
        realized face carries ``ReflectiveBoundary(axis='x')`` at
        ``test_b3_domain_narrowing.py:164,175``), so no committed assertion
        moved — which is why this capability needs a gate of its own or it
        lands unrecorded.
        """
        from tests._harness.references import mirror_partner_indices

        quad = Quadrature.folded_product(n_mu=4, n_phi=8)
        reference = mirror_partner_indices(quad, "y")
        np.testing.assert_array_equal(reference, np.arange(quad.n_ordinates))

        kernel = quad.measure.permutation_under(RigidMotion.reflection(normal=_R4_I3[1]), atol=1e-9)
        assert kernel is not None
        np.testing.assert_array_equal(np.asarray(kernel.indices), reference)

        # the axis the two folded CALL SITES actually pass is untouched
        np.testing.assert_array_equal(
            mirror_partner_indices(quad, "x"),
            mirror_partner_indices(quad, "x"),
        )
