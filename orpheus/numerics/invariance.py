r"""Invariance of a discrete measure under a subgroup of :math:`O(3)` — asked
ON the measure's orbit space.

The question is the MEASURE's: a measure carries its nodes, its weights and the
:class:`~orpheus.numerics.manifold.Quotient` it lives on; a group carries its
realization.  So the verbs live on :class:`~orpheus.numerics.measure.DiscreteMeasure`
— :meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.certificate_under`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.permutation_under`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.singular_set_under`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.symmetry_groups` — and delegate
here, where ONE closure (:func:`_orbit_space_closure`) serves all of them: the
measure's nodes are carried to the base's ambient space as orbit barycentres,
read as CHART points of the orbit space, and each isometry's image is its
induced action (:meth:`~orpheus.numerics.manifold.Quotient.induced_action`),
which exists only for an isometry that normalises the quotienting group.

The import graph reads like the mathematics (R2 of #434, 2026-09-03):
``geometry.transformation`` ← ``symmetry`` (groups) ← ``manifold`` (an orbit
space is a manifold and a group) ← ``measure`` (a measure lives on a manifold)
← this module (a measure × a group).  Until R2 the kernel lived in
:mod:`~orpheus.numerics.symmetry`, which therefore imported ``measure`` at
module scope; ``manifold`` could reference a group only under ``TYPE_CHECKING``
and duck-typed ten of its members (24 sites); ``measure`` deferred its import of the
certificate to function scope; and the kernel inlined a copy of the closure
(`[M]` an identical lambda in two functions) while three docstrings claimed
"one closure".  None of those survives: this module imports ``measure`` under
``TYPE_CHECKING`` only (it READS measures, never builds one), so no cycle
exists to defer around.

**Continuous groups are decided EXACTLY, never sampled (ERR-072):** the
identity component must fix every node's orbit barycentre — a connected orbit
inside a finite set is a point — and the finitely many coset representatives
are permuted like finite elements.  **The match is a bijection (ERR-073)**, the
nodes are matched at :data:`_NODE_WINDOW_FACTOR` × the weight window, and the
certificate the closure returns (:class:`OrbitCertificate`) carries every
permutation, so the singular set is an INTEGER identity.

The lattice walk (:func:`symmetry_groups`) descends the Hasse diagram of the
expressible groups (:func:`candidate_groups`) with both prunings — invariance is
downward-closed — and returns the MAXIMAL invariant candidates; the true
:math:`\mathrm{Sym}(\mu)` is the group they generate (#437 computes it from the
node set instead).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Iterable, Literal

import numpy as np

from orpheus.geometry.transformation import (
    Permutation,
    RigidMotion,
    permutation_preserving,
)

from .manifold import Quotient, RealSpace
from .symmetry import SubgroupOfO3

if TYPE_CHECKING:  # this module READS measures; it never builds one
    from .measure import DiscreteMeasure

#: The WEIGHT window every verb defaults to — the absolute tolerance at which
#: two weights are one weight.  Spelled once: the measure verbs and the
#: quadrature's fold default to THIS name, so the day the window becomes
#: relative (see :data:`_NODE_WINDOW_FACTOR`) one edit moves every caller.
WEIGHT_ATOL = 1e-13

#: The node-match window is this multiple of ``atol``, which is the WEIGHT
#: window. The asymmetry is deliberate and worth naming rather than leaving
#: as a bare ``* 100``: a node coordinate is the accumulated result of a
#: matrix product against a constructed direction cosine, so it carries more
#: round-off than a weight, which is usually read straight from a table. Both
#: windows are ABSOLUTE — for rules whose weights are O(1e-3) the weight test
#: is correspondingly stricter in relative terms, which is a known
#: characteristic of this check and not an accident of it.
_NODE_WINDOW_FACTOR = 100.0


def _node_window(weight_atol: float) -> float:
    """The node-match window that goes with a weight window."""
    return weight_atol * _NODE_WINDOW_FACTOR


# ---------------------------------------------------------------------------
# Invariance check
# ---------------------------------------------------------------------------


def _embedded_nodes(measure: DiscreteMeasure) -> np.ndarray:
    r"""The measure's nodes as orbit BARYCENTRES in :math:`\mathbb{R}^3`.

    A measure on an orbit space answers through its entry: the nodes are
    carried to the base's ambient space by :meth:`Quotient.orbit_barycentres
    <orpheus.numerics.manifold.Quotient.orbit_barycentres>` — the Reynolds
    projector of a fold's representatives (since #434 R4, 2026-09-03; until
    then a fold's nodes passed through as representatives), the entry's
    :attr:`~orpheus.numerics.manifold.Quotient.lift` of a polar marginal's
    chart coordinates (:math:`\mu \mapsto \mu\,\hat e_a`, inside the ball,
    on the sphere only at the poles) — the honest point an invariance
    question wants, since every isometry that descends to the orbit space
    carries barycentres to barycentres. A measure that names no orbit space
    keeps the tree's zero-padding convention: a chart-level interval's
    :math:`\mu` becomes :math:`(\mu, 0, 0)` (column 0 — a bare interval
    names no axis), a planar rule :math:`(x, y)` becomes :math:`(x, y, 0)`,
    a sphere rule is itself.

    ⭐ Since 2026-09-02 (#429 tracker 2.2b) the axis is READ off the entry's
    lift rather than by this function: until then it embedded a polar
    marginal through :func:`~orpheus.numerics.manifold.barycentre` after
    reading the axis itself, and the invariance question was then asked
    in the AMBIENT space — right for a bare sphere rule, wrong for a fold,
    whose representatives are not closed under the group that folded them
    although that group acts trivially on the orbit space. The question
    is now asked ON the orbit space (:func:`is_invariant_under`);
    this function only supplies the representatives.
    """
    support = measure.support
    if isinstance(support, Quotient):
        return support.orbit_barycentres(measure.nodes)
    nodes = measure.nodes
    if nodes.ndim == 1:
        nodes = nodes[:, None]
    n, d = nodes.shape
    if d == 3:
        return nodes
    embedded = np.zeros((n, 3))
    embedded[:, :d] = nodes
    return embedded


def _as_columns(points: np.ndarray) -> np.ndarray:
    """``(n,)`` chart coordinates as the ``(n, 1)`` column a point-set match wants."""
    arr = np.asarray(points, dtype=float)
    return arr if arr.ndim == 2 else arr[:, None]


@functools.cache
def _ambient_orbit_space() -> Quotient:
    r"""The trivial orbit space :math:`\mathbb{R}^3/\{e\}` — where a measure
    that names no orbit space is asked. Its chart is the ambient space
    itself, its lift the identity, and every isometry descends to it, so
    the orbit-space kernel reduces on it to the ambient question. The base
    is :math:`\mathbb{R}^3` and not the sphere on purpose: a zero-padded
    interval or planar rule lands OFF the sphere, and the container must
    honestly contain what is put in it."""
    return RealSpace(3).quotient(SubgroupOfO3.Trivial)


def _orbit_space_of(measure: DiscreteMeasure) -> Quotient:
    r"""The orbit space a measure is asked ON — its support when that is a
    :class:`~orpheus.numerics.manifold.Quotient`, else the trivial orbit
    space :math:`\mathbb{R}^3/\{e\}` (:func:`_ambient_orbit_space`), where
    the two readings coincide.  ONE spelling, so the closure and the
    invariance kernel cannot disagree about which space they are on."""
    support = measure.support
    return support if isinstance(support, Quotient) else _ambient_orbit_space()


def _orbit_space_closure(
    measure: DiscreteMeasure,
    motions: Iterable[RigidMotion],
    atol: float,
) -> OrbitCertificate | None:
    r"""The certificate of ``motions`` acting ON the measure's orbit space —
    the one closure every "does this isometry permute these ordinates"
    question in the tree is answered by (#429 tracker 2.2b).

    The measure's nodes are carried to the base (:func:`_embedded_nodes`),
    read as CHART points of the orbit space (a bare support is the trivial
    orbit space :math:`\mathbb{R}^3/\{e\}`, where chart and base coincide),
    and each motion's image of the chart set is its
    :meth:`~orpheus.numerics.manifold.Quotient.induced_action` — which
    exists only for a motion that normalises the quotienting group, so a
    motion that does not act on the orbit space is ``None`` here, exactly
    as one that acts without permuting is. The match itself is
    :func:`_orbit_closure`: ERR-073's bijection guard, ERR-074's no-bare-
    ``argmin`` guard, the weight leg, and the two windows.
    """
    orbit_space = _orbit_space_of(measure)
    barycentres = _embedded_nodes(measure)
    chart = _as_columns(orbit_space.orbit_coordinates(barycentres))
    motions = list(motions)
    if not all(orbit_space.by.is_normalised_by(g) for g in motions):
        return None
    return _orbit_closure(
        chart,
        measure.weights,
        motions,
        atol,
        images_of=lambda g: _as_columns(orbit_space.induced_action(g)(barycentres)),
    )


def permutation_under(
    measure: DiscreteMeasure, motion: RigidMotion, *, atol: float,
) -> Permutation | None:
    r"""The permutation ``motion`` induces on the measure's weighted nodes,
    ON the measure's orbit space — or ``None`` if it is not a symmetry there.

    The single-motion face of :func:`certificate_under` and the kernel
    :func:`is_invariant_under` reads, so that a quadrature's
    :meth:`~orpheus.numerics.quadrature.directional.Quadrature.ordinate_permutation`
    and its invariance can never disagree about a fold: `[M]` until
    2026-09-02 they did — :math:`\sigma_y` on ``folded_product(4, 8)`` was
    *invariant* to one reading and *no permutation* to the other, because
    the second matched the fold's representatives in the ambient space where
    their :math:`\sigma_y`-mates are absent. On the orbit space it is the
    identity permutation, which is what a reflecting face on a folded rule
    realizes. ``atol`` is the WEIGHT window; the node window is
    :data:`_NODE_WINDOW_FACTOR` times it.
    """
    certificate = _orbit_space_closure(measure, [motion], atol)
    return None if certificate is None else certificate.permutations[0]


def is_invariant_under(
    measure: DiscreteMeasure, group: SubgroupOfO3, *, atol: float = WEIGHT_ATOL,
) -> bool:
    r"""``True`` iff ``measure`` is invariant under every element of ``group``,
    asked ON the measure's orbit space — the one kernel.

    Three facts, in the order they are cheapest and most decisive:

    1. :math:`G` must NORMALISE :math:`H` (the group the measure's orbit space
       :math:`M/H` was quotiented by), else it does not act on :math:`M/H` at
       all (:meth:`SubgroupOfO3.normalises
       <orpheus.numerics.symmetry.SubgroupOfO3.normalises>`, exact through the
       Lie algebra) — ``False``.
    2. A CONTINUOUS :math:`G` is decided exactly through its realization: the
       identity component :math:`G^0` has connected orbits, and a connected
       orbit inside a finite set is a point, so :math:`G^0` must FIX every
       node's orbit BARYCENTRE (:meth:`IdentityComponent.fixes
       <orpheus.numerics.symmetry.IdentityComponent.fixes>` — the axis-support
       / origin rule of ERR-072, stated once), at the NODE window
       :data:`_NODE_WINDOW_FACTOR` × ``atol`` (a position question; until R2
       of #434 it ran at the weight window while the node match a step later
       used 100× it).  Exact on the barycentre because the chart is injective
       on orbits, so :math:`G^0` fixes :math:`P_H p` iff it fixes :math:`[p]`;
       when :math:`G^0 \subseteq H` the barycentre lies on the axis by
       construction and the test passes as the theorem says it must.
    3. Every element (finite :math:`G`) or coset representative (continuous
       :math:`G`) must permute the weighted node set in CHART coordinates
       through its induced action (:func:`_orbit_space_closure` — ERR-073's
       bijection guard, ERR-074's no-bare-``argmin`` guard, the weight leg).

    ⚠ There is no "``G ⊆ H`` ⇒ ``True``" short circuit.  It was an
    OPTIMISATION, and the closure re-proves it — a group inside :math:`H`
    acts trivially on :math:`M/H`, so every element's induced action is the
    identity permutation (`[M]` deleting it reddened 0 of 2290 kernel calls;
    32 shipped (rule × group) rows stay green through the closure).  Do not
    restore it as a missing guard.  ⚠ And do not apply the same argument to
    step 1: for a FINITE group its representatives are its elements and the
    closure's per-motion guard re-proves normalisation, but for a CONTINUOUS
    group the representatives do not cover :math:`G^0`, and the
    ``component`` half of :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.normalises`
    is the only check that :math:`G^0` normalises :math:`H` — deleting it
    is a silent wrong ``True`` in ERR-072's direction.

    ``atol`` is the WEIGHT window; the node window is
    :data:`_NODE_WINDOW_FACTOR` times it.
    """
    orbit_space = _orbit_space_of(measure)
    if not group.normalises(orbit_space.by):
        return False
    realization = group.realization
    if not realization.is_finite:
        barycentres = _embedded_nodes(measure)
        if not realization.component.fixes(barycentres, atol=_node_window(atol)):
            return False
    return _orbit_space_closure(measure, realization.representatives, atol) is not None


# ---------------------------------------------------------------------------
# Orbit-closure check, and the certificate it produces
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OrbitCertificate:
    r"""The realized action of a group on a point set.

    Returned by the closure check instead of a bare ``True``. The check
    ALREADY computes, for each operator :math:`M` and node :math:`i`, the
    index :math:`\pi_M(i)` with :math:`M x_i = x_{\pi_M(i)}` — it simply
    discarded it. A ``-> bool`` predicate that internally builds the
    permutation IS the missing primitive; widening the return type is
    cheaper than minting a class to recompute it.

    Everything below is a *reading* of the permutations, not new work:

    * :attr:`singular_set` — the orbifold **singular set**
      :math:`\Sigma = \{x : \mathrm{Stab}(x) \neq \{e\}\}`. Because
      :math:`\pi_M(i) = i` means exactly :math:`M x_i = x_i`, i.e.
      :math:`x_i \in \mathrm{Fix}(M)`, membership is an **integer
      identity** — exact, with no tolerance. The only place a tolerance
      enters is matching nodes while BUILDING :math:`\pi`, which is the
      one place the question is honestly numerical.
    * :attr:`stabilizer_order` — :math:`|\mathrm{Stab}(x_i)|`, the orbit
      type; by orbit-stabilizer the orbit length is :math:`|G|` divided
      by it.
    * :meth:`orbits` — the :math:`G`-orbits themselves.

    A certificate exists only for a :math:`G`-INVARIANT set, so
    :math:`\Sigma` is unrepresentable without the closure proof — which is
    the precondition ("the quotient is defined only on a G-invariant
    measure") enforced by construction rather than by a comment.
    """

    operators: tuple[RigidMotion, ...]
    permutations: tuple[Permutation, ...]

    @property
    def n_points(self) -> int:
        return self.permutations[0].n if self.permutations else 0

    @property
    def _non_identity(self) -> list[Permutation]:
        """The permutations of operators other than the identity.

        The identity fixes every point, so including it would report the
        whole set as singular. The test asks the ELEMENT whether it is the
        identity of its own dimension, rather than comparing against a
        hard-coded ``np.eye(3)`` — the certificate is a statement about a
        group acting on a point set, and nothing in it is three-dimensional.
        """
        return [
            perm
            for M, perm in zip(self.operators, self.permutations)
            if not M.approximately_equals(RigidMotion.identity(M.dimension))
        ]

    @property
    def singular_set(self) -> np.ndarray:
        r"""Indices of the **singular points** — :math:`\Sigma`.

        A point is singular iff some non-identity group element fixes it.
        Under a reflection the fixed locus is a **mirror**; a point fixed
        by two mirrors is a corner reflector; a point on a rotation axis
        with no mirror is a cone point.
        """
        fixed = np.zeros(self.n_points, dtype=bool)
        for perm in self._non_identity:
            fixed[perm.fixed_points] = True
        return np.flatnonzero(fixed)

    @property
    def stabilizer_order(self) -> np.ndarray:
        r""":math:`|\mathrm{Stab}(x_i)|` for every node ``i``.

        ``1`` exactly on the regular (free) points; the singular set is
        ``> 1``. Meaningful only when the certificate was built from the
        FULL group rather than a generating set.
        """
        order = np.zeros(self.n_points, dtype=np.int64)
        for perm in self.permutations:
            order[perm.fixed_points] += 1
        return order

    def orbits(self) -> tuple[np.ndarray, ...]:
        """The :math:`G`-orbits, as arrays of node indices."""
        parent = list(range(self.n_points))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for perm in self.permutations:
            for i, j in enumerate(perm.indices):
                ra, rb = find(i), find(int(j))
                if ra != rb:
                    parent[ra] = rb

        buckets: dict[int, list[int]] = {}
        for i in range(self.n_points):
            buckets.setdefault(find(i), []).append(i)
        return tuple(np.array(v, dtype=np.int64) for v in buckets.values())


def _orbit_closure(
    nodes: np.ndarray,
    weights: np.ndarray,
    ops: Iterable[RigidMotion],
    atol: float,
    *,
    images_of: Callable[[RigidMotion], np.ndarray],
) -> OrbitCertificate | None:
    """Check that applying every operator ``M ∈ ops`` to ``nodes``
    yields the same multiset of (node, weight) pairs.

    The per-element work is
    :func:`~orpheus.geometry.transformation.permutation_preserving`; what
    remains here is the *measure*-level question — that EVERY element
    preserves it, and the certificate assembled from the results. The
    matching, the bijectivity requirement (ERR-073) and the weight guard all
    live in the core, so this module no longer carries a second copy of them.
    ``images_of`` supplies each operator's image of the node set — its
    INDUCED action on the orbit space, applied to chart coordinates (#429
    tracker 2.2b).  REQUIRED since R2 of #434: the ambient default was dead
    (every caller passed the action) and the one caller is
    :func:`_orbit_space_closure`, which is what makes "one closure" a
    structural fact rather than a docstring's claim.

    The two windows: ``atol`` is the WEIGHT tolerance, and the node-match
    window is :data:`_NODE_WINDOW_FACTOR` times larger. That relationship is
    a claim about *measures* — a coordinate accumulates matrix-product
    round-off while a weight is usually read from a table — so it is stated
    here rather than inside the geometry primitive, which takes the two
    windows separately and refuses to guess.

    For each :math:`M`, compute :math:`M(\\text{nodes})` and find a
    permutation :math:`\\pi` such that
    :math:`M(\\text{nodes})_i = \\text{nodes}_{\\pi(i)}` (within
    ``atol``). The match is verified element-wise, then the weight
    equality :math:`w_i = w_{\\pi(i)}` is checked.

    The match is required to be a **bijection**. Nearest-neighbour
    matching alone only proves every image has *a* same-weight partner,
    which is strictly weaker than "the action permutes the nodes": two
    distinct sources may land on one target, leaving some node
    unmatched entirely. Such a set is NOT :math:`G`-invariant — its
    mass is distributed differently from its image — yet it passes an
    injectivity-free check. Appending a bit-identical duplicate of any
    node to an :math:`O_h`-invariant rule is the minimal witness: the
    duplicated position then carries twice the mass of its mirror
    image, and every one of the 48 match maps is non-injective.
    Since :math:`\\pi` maps an :math:`n`-set to an :math:`n`-set,
    injectivity is equivalent to bijectivity, so counting distinct
    targets suffices. (ERR-073.)

    Returns ``None`` at the first failure, else the
    :class:`OrbitCertificate` carrying every :math:`\\pi_M`.
    """
    kept_ops: list[RigidMotion] = []
    perms: list[Permutation] = []

    for M in ops:
        images = images_of(M)
        pi = permutation_preserving(
            nodes,
            images,
            weights,
            atol=_node_window(atol),
            weight_atol=atol,
        )
        if pi is None:
            return None
        kept_ops.append(M)
        perms.append(pi)
    return OrbitCertificate(operators=tuple(kept_ops), permutations=tuple(perms))


# ---------------------------------------------------------------------------
# The subgroup graph, and walking it to find a measure's symmetry
# ---------------------------------------------------------------------------
#
# Crystallography does not ASK a structure what its symmetry group is — it
# walks the subgroup lattice downward from high symmetry until it reaches
# the symmetry the structure actually has. The graph is the Hasse diagram of
# the lattice: nodes are groups, edges are MAXIMAL-subgroup relations
# (International Tables Vol. A1 renders it as a Bärnighausen tree).
#
# That is the machinery this module needs. A DECLARED invariance group is a
# claim with no construction behind it, and such claims have already shipped
# false here twice. A COMPUTED one cannot lie about the object it was
# computed from.


def _distinct_azimuths(nodes: np.ndarray, atol: float) -> int:
    r"""Number of distinct azimuthal angles among the off-axis nodes.

    Bounds the cyclic families: a :math:`C_n` rotation (``n > 1``) fixes no
    azimuth, so the azimuths split into FREE orbits of size ``n`` and
    therefore ``n`` divides this count. That turns two infinite families
    into a handful of divisors.  ``atol`` is BOTH the off-axis cut and the
    angle-merge window (until R2 of #434 a literal ``1e-9`` shadowed the
    argument for the merge, so the function's own parameter could not move
    its answer — `[M]` inert on all 15 shipped rules, witnessed by a
    manufactured pair of near-coincident azimuths).
    """
    rho = np.hypot(nodes[:, 0], nodes[:, 1])
    off_axis = nodes[rho > atol]
    if off_axis.shape[0] == 0:
        return 0
    phi = np.sort(np.mod(np.arctan2(off_axis[:, 1], off_axis[:, 0]), 2.0 * np.pi))
    distinct = [phi[0]]
    for p in phi[1:]:
        if p - distinct[-1] > atol:
            distinct.append(float(p))
    # The first and last can be the same angle seen either side of the branch.
    if len(distinct) > 1 and (2.0 * np.pi - distinct[-1] + distinct[0]) < atol:
        distinct.pop()
    return len(distinct)


def candidate_groups(
    measure: DiscreteMeasure, *, atol: float = WEIGHT_ATOL,
) -> tuple[SubgroupOfO3, ...]:
    """The expressible groups worth testing against ``measure``.

    The named entries always, plus :math:`C_n` / :math:`D_{nh}` for each
    divisor of the measure's distinct-azimuth count (see
    :func:`_distinct_azimuths` for why divisors suffice).
    """
    named = [
        SubgroupOfO3.Trivial,
        SubgroupOfO3.Dinfh, SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh, SubgroupOfO3.SO3, SubgroupOfO3.O3,
    ]
    # All three mirrors, always. The parameter set of a reflection family
    # is {x, y, z} — FINITE BY CONSTRUCTION, so unlike Cn/Dnh (which need
    # the distinct-azimuth divisor bound to stay finite) it needs no bound
    # at all, and unlike the retired parameter-free Z2 it offers the walk
    # all three planes instead of silently only sigma_z.
    named += [SubgroupOfO3.Mirror(a) for a in ("x", "y", "z")]
    # And all three axial rotation groups, for the same reason: offering
    # only z (what the retired bare SO2 amounted to) made every x-pole
    # rule — the slab's own polar marginal — read as carrying no
    # continuous symmetry at all.
    named += [SubgroupOfO3.SO2(a) for a in ("x", "y", "z")]
    # And the three stabilisers O(2)_a above them (#432, 2026-09-02).
    # Finite-vs-finite the walk never sees them as equal, and on a polar
    # marginal they are what is MAXIMAL: the slab's rule reports O(2)_x
    # (with its μ → −μ mirror σ_x beside it, since σ_x flips the axis and
    # so lies in neither) where it used to stop at the rotations.
    named += [SubgroupOfO3.O2(a) for a in ("x", "y", "z")]
    # The azimuth count is read off the orbit BARYCENTRES in R^3 — the same
    # points the invariance kernel reads — never off the stored width: until
    # R2 of #434 this branched on `measure.nodes.shape[1]`, so ONE fold had
    # two candidate sets by spelling (`[M]` 20 names at ambient width, 15 at
    # chart width, reporting {D_2h} and {sigma_x, sigma_y, sigma_z}).  On the
    # barycentres a polar marginal's set is (mu, 0, 0), two z-azimuths, so
    # D_2h is offered and absorbs sigma_x in the walk's answer.
    n_az = _distinct_azimuths(_embedded_nodes(measure), _node_window(atol))
    families: list[SubgroupOfO3] = []
    for d in (d for d in range(1, n_az + 1) if n_az % d == 0) if n_az else (1,):
        if d > 1:  # C_1 IS Trivial, already in `named`; D_1h is a real group of order 4
            families.append(SubgroupOfO3.Cn(d))
        families.append(SubgroupOfO3.Dnh(d))
    return tuple(named + families)


def _maximal(groups: Iterable[SubgroupOfO3]) -> tuple[SubgroupOfO3, ...]:
    """Those members not STRICTLY contained in another member.

    Strict containment — ``h ⊇ g`` and not ``g ⊇ h`` — as the docstring has
    always said; until R2 of #434 the body tested ``h != g and h.contains(g)``,
    which drops BOTH members of an equal-but-distinct pair (`[M]` with
    ``Cn(1)`` a second value for ``{e}`` that emptied the walk's answer on a
    symmetry-free measure; R1 merged the pair, so no expressible pair can
    tell the two predicates apart today — the gate's witness is synthetic).
    """
    items = list(groups)
    return tuple(
        g for g in items
        if not any(h.contains(g) and not g.contains(h) for h in items)
    )


def maximal_subgroups(
    group: SubgroupOfO3, candidates: Iterable[SubgroupOfO3],
) -> tuple[SubgroupOfO3, ...]:
    """The Hasse edges below ``group`` — its maximal proper subgroups.

    Derived from :meth:`SubgroupOfO3.contains`, never tabulated: a
    hand-drawn Bärnighausen tree would be exactly the unverifiable claim
    the computed lattice exists to eliminate.
    """
    below = [h for h in candidates if h != group and group.contains(h)]
    return _maximal(below)


def symmetry_groups(
    measure: DiscreteMeasure,
    *,
    atol: float = WEIGHT_ATOL,
    candidates: Iterable[SubgroupOfO3] | None = None,
    method: Literal["walk", "bruteforce"] = "walk",
) -> tuple[SubgroupOfO3, ...]:
    r"""The maximal groups leaving ``measure`` invariant — its symmetry.

    Returns the maximal elements of
    :math:`\{G \in C : \mu \text{ is } G\text{-invariant}\}`. Several
    incomparable maxima can survive (a rule may carry both a rotation group
    and a reflection with neither inside the other); the true
    :math:`\mathrm{Sym}(\mu)` is the group they GENERATE, which need not be
    expressible in the candidate set — so the maximal elements, not a single
    tag, are the honest answer.

    **Soundness.** Invariance is DOWNWARD-CLOSED: if :math:`\mu` is
    :math:`G`-invariant and :math:`H \le G`, then :math:`\mu` is
    :math:`H`-invariant. So the invariant set is an order ideal, giving the
    walk both prunings — a failing node kills every supergroup, a passing
    node implies every subgroup. The walk therefore *requires* a correct
    lattice and is simultaneously a test of it.

    Parameters
    ----------
    method : {"walk", "bruteforce"}
        Two realizations of the same definition. ``"walk"`` descends the
        Hasse diagram with both prunings; ``"bruteforce"`` tests every
        candidate. They must agree on every input — that agreement is the
        verification instrument, not a mere fast-path check.

    Notes
    -----
    The answer is about the group's realization in the **standard setting**
    (principal axis along z for the finite families :math:`C_n` /
    :math:`D_{nh}` and for :math:`D_{\infty h}`), not up to conjugation. A
    rule whose symmetry axis is not z reports a smaller group from those
    families, which is correct for a gate comparing against a geometry in
    the same frame. The three families whose parameter IS the axis —
    :class:`Mirror`, :class:`SO2` and :class:`O2` — are offered on all three
    axes, so a polar marginal along :math:`x` reports :math:`O(2)_x`.
    """
    cands = tuple(candidates) if candidates is not None else candidate_groups(
        measure, atol=atol
    )

    if method == "bruteforce":
        return _maximal(
            g for g in cands if is_invariant_under(measure, g, atol=atol)
        )
    if method != "walk":
        raise ValueError(
            f"method must be 'walk' or 'bruteforce', got {method!r}"
        )

    accepted: list[SubgroupOfO3] = []
    visited: set[SubgroupOfO3] = set()
    stack = list(_maximal(cands))
    while stack:
        group = stack.pop()
        if group in visited:
            continue
        visited.add(group)
        if any(a.contains(group) for a in accepted):
            continue  # already implied by an accepted supergroup
        if is_invariant_under(measure, group, atol=atol):
            accepted.append(group)  # everything below is implied — do not descend
        else:
            stack.extend(maximal_subgroups(group, cands))
    return _maximal(accepted)


def certificate_under(
    measure: DiscreteMeasure,
    group: SubgroupOfO3,
    *,
    atol: float = WEIGHT_ATOL,
) -> OrbitCertificate | None:
    r"""The realized action of ``group`` on ``measure``, ON the measure's
    orbit space — or ``None``.

    ``None`` when the group is CONTINUOUS (no finite element set to
    permute nodes by), when it does not act on the measure's orbit space
    (it does not normalise the quotienting group), or when some element
    does not permute the weighted node set there. The certificate is built
    from the group's ELEMENTS, not a generating set: orbit closure only
    needs generators, but a point's stabilizer may be generated by a
    non-generator, so :math:`\Sigma` needs them all. Its ``operators`` are
    the ambient motions; its ``permutations`` are of the measure's nodes
    read as CHART points of the orbit space under each element's
    :meth:`~orpheus.numerics.manifold.Quotient.induced_action` — for a bare
    support (the trivial orbit space) exactly the ambient permutations,
    for a fold the permutations of the REPRESENTATIVES, for a polar
    marginal the permutations of :math:`\mu` (#429 tracker 2.2b,
    2026-09-02; until then a 1-D node set was refused here by SHAPE, the
    II.11 defect).  One closure (:func:`_orbit_space_closure`) serves
    this, :func:`is_invariant_under` and
    :func:`permutation_under`, so they cannot disagree.
    """
    realization = group.realization
    if not realization.is_finite:
        return None
    return _orbit_space_closure(measure, realization.elements, atol)


def singular_set_under(
    measure: DiscreteMeasure,
    group: SubgroupOfO3,
    *,
    atol: float = WEIGHT_ATOL,
) -> np.ndarray:
    r"""The **singular set** :math:`\Sigma` of ``measure`` under ``group``.

    Indices of the **singular points** — those whose stabilizer is
    non-trivial, i.e. fixed by some non-identity group element. Under a
    reflection the fixed locus is a **mirror**; a point on two mirrors is a
    corner reflector; a point on a rotation axis with no mirror through it
    is a cone point.

    Membership is decided by :math:`\pi_M(i) = i` — an **integer identity**,
    exact, with no tolerance. The three ad-hoc float comparisons the tree
    grew for this question (``_OCTANT_SIGN_EPS``, ``_MU_DIRECTION_EPS``,
    ``TANGENTIAL_EPS``) were all asking it numerically; measured across 29
    production rules, the separation between "zero" and "nonzero" cosines is
    a factor of :math:`2.7\times10^{13}`, so the tolerance was never doing
    real work. The only honestly-numerical step is matching nodes while
    BUILDING :math:`\pi`, which is ``atol``.

    Raises
    ------
    ValueError
        If ``measure`` is not ``group``-invariant. :math:`\Sigma` is defined
        only on an invariant set — a quotient needs something to quotient —
        so the illegal state is unrepresentable rather than silently
        returning a wrong answer.
    """
    cert = certificate_under(measure, group, atol=atol)
    if cert is None:
        raise ValueError(
            f"the singular set is defined only for a {group.name}-invariant "
            f"measure with a finite realization; no certificate exists "
            f"because {group.name} is continuous (no finite node "
            f"permutation), or does not ACT on this measure's orbit space "
            f"{measure.support.name} (it does not normalise the group "
            f"already spent there), or acts there without permuting the "
            f"weighted nodes — this measure is not {group.name}-invariant "
            f"on its orbit space."
        )
    return cert.singular_set


#: The free functions carry the measure verbs' names — ``certificate_under``,
#: ``permutation_under``, ``singular_set_under`` were ``orbit_certificate``,
#: ``induced_permutation``, ``singular_set`` until R2 of #434 landed
#: (2026-09-03) — so the family greps as ONE word order (a rename census
#: keyed on either spelling had missed the other's docstrings).
__all__ = [
    "OrbitCertificate",
    "candidate_groups",
    "permutation_under",
    "is_invariant_under",
    "maximal_subgroups",
    "certificate_under",
    "singular_set_under",
    "symmetry_groups",
]
