r"""Sphere quadrature rules on :math:`S^2`.

Two families:

* **Lebedev rules** — :math:`O_h`-invariant grids on the unit sphere
  with high polynomial exactness for the order chosen (Lebedev 1976).
  Returned by :func:`lebedev_sphere`.
* **Level-symmetric :math:`S_N` rules** — Carlson-Lathrop's
  :math:`O_h`-invariant triangular discrete-ordinate grids organised
  into :math:`N/2` polar levels per hemisphere with permutations
  :math:`(\eta, \xi, \mu)` per octant (Carlson & Lathrop 1968, Lewis &
  Miller §4.2). Returned by :func:`level_symmetric_sn`.

Both families are :math:`O_h`-invariant, so the returned measures
carry ``invariance_group=SubgroupOfO3.OctahedralOh``.

The level-symmetric rule additionally returns a
:class:`LevelStructure` *alongside* the measure — a frozen dataclass
carrying the per-level metadata the cylindrical SN sweep needs for
azimuthal redistribution.
:meth:`Quadrature.level_symmetric
<orpheus.numerics.quadrature.Quadrature.level_symmetric>` stores it on
the ``level_structure`` field, so the sweep reads it in :math:`O(1)`.
There is no per-family adapter class: the four SN-side wrappers this
docstring used to point at (``LevelSymmetricSN`` and its siblings)
were retired into classmethod factories on the one ``Quadrature``
type (R-1 Phase A detour-C).

References
----------

* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Comp. Math.*
  **16**(2), 10-24.
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan et al., eds., Gordon & Breach.
* Lewis, E.E. and Miller, W.F. (1993). *Computational Methods of
  Neutron Transport*. Wiley. §4.2 (level-symmetric construction
  and degree of exactness).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from enum import Enum
from functools import lru_cache

import numpy as np

from ..exactness import UNIFORM_ON_SPHERE, ExactnessClaim
from ..measure import SPACE_SPHERE, DiscreteMeasure
from ..symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Lebedev sphere quadrature
# ---------------------------------------------------------------------------


def lebedev_sphere(order: int) -> DiscreteMeasure:
    r"""Lebedev :math:`O_h`-invariant rule on :math:`S^2` of given order.

    Wraps :func:`scipy.integrate.lebedev_rule`. The ``order`` parameter
    is the maximum degree of the polynomial integrated exactly, so the
    returned measure carries ``degree_of_exactness=order``. The number
    of nodes is determined by SciPy's tabulated Lebedev grid for that
    order — typical sizes are 6, 14, 26, 38, 50, …, 590.

    Weight sum: :math:`\sum_i w_i = 4\pi`, matching the area of
    :math:`S^2`.

    Symmetry: :math:`O_h` (full octahedral group, 48 elements). The
    construction is :math:`O_h`-invariant by design (Lebedev 1976) —
    every node sits in an :math:`O_h`-orbit and the corresponding
    weights are equal across the orbit.

    Nodes are shape ``(N, 3)`` carrying the Cartesian direction
    cosines :math:`(\mu_x, \mu_y, \mu_z)` per row.

    Parameters
    ----------
    order : int
        Lebedev order; must match a tabulated value (SciPy raises
        otherwise).

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(N, 3)``, weights shape ``(N,)``, on
        ``space="S^2"``, with ``invariance_group=SubgroupOfO3.OctahedralOh`` and
        ``degree_of_exactness=order``.

    See Also
    --------
    :meth:`orpheus.numerics.quadrature.Quadrature.lebedev` — the named
    factory SN consumers call. It wraps this measure; mirror pairings
    are derived on demand via ``ordinate_permutation``. There is no per-family
    adapter class: the SN-side wrapper this docstring used to point at
    (``orpheus.sn.quadrature.LebedevSphere``) was retired into a
    classmethod factory on the one ``Quadrature`` type, whose
    ``mu_x`` / ``mu_y`` / ``mu_z`` are ``@property`` views over these
    nodes rather than separately-cached fields.
    """
    from scipy.integrate import lebedev_rule

    pts, w = lebedev_rule(order)
    # ``pts`` from scipy has shape ``(3, N)``; transpose to ``(N, 3)``
    # for the canonical ``DiscreteMeasure`` (N, d) layout.
    nodes = np.ascontiguousarray(pts.T)  # (N, 3)
    return DiscreteMeasure(
        nodes=nodes,
        weights=w,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        # SPHERICAL-HARMONIC degree, against Lebesgue measure on S^2 —
        # not an algebraic degree, and not against a weight on an
        # interval. Lebedev's reference generates nothing by
        # Golub-Welsch: the claim's authority is the published table
        # (Lebedev 1976), which is exactly why the reference here is a
        # ``UniformMeasure`` rather than a ``GeneratingMeasure``.
        exactness=ExactnessClaim(reference=UNIFORM_ON_SPHERE, degree=order),
    )


# ---------------------------------------------------------------------------
# Level-symmetric S_N quadrature
# ---------------------------------------------------------------------------


class PolarInvariant(Enum):
    r"""Which invariant a quadrature's levels are the fibers OF.

    A level is a **fiber of an invariant**, not an orbit, and the two
    shipped producers fiber over *different* invariants — which the
    single ``LevelStructure`` type used to erase.

    * :attr:`SIGNED_MU_Z` — :func:`product_mu_phi`. Each Gauss-Legendre
      polar node is its own level, so the levels run over the whole
      range :math:`\mu_z \in [-1, 1]` and the fiber over a level is a
      full circle.
    * :attr:`ABS_MU_Z` — :func:`level_symmetric_sn`. Levels are indexed
      by :math:`|\mu_z|`, so each carries BOTH hemispheres and the
      fiber is two circles.

    Consequence, and the reason this has to be on the type: "the
    ordinates of level :math:`p`" means a different set in the two
    cases, so any consumer that reads ``level_indices`` without knowing
    the invariant is reading an ambiguous object.
    """

    SIGNED_MU_Z = "signed_mu_z"
    ABS_MU_Z = "abs_mu_z"


@dataclass(frozen=True)
class LevelStructure:
    r"""Per-level metadata for a triangular SN-style sphere quadrature.

    Captured alongside the :class:`DiscreteMeasure` returned by
    :func:`level_symmetric_sn` and :func:`product_mu_phi`. The cylindrical
    SN sweep needs this structure to compute the azimuthal redistribution
    coefficients (Bailey et al. 2009, Eq. 50).

    Attributes
    ----------
    n_levels : int
        Number of polar levels per hemisphere (or polar axis, for the
        product-quadrature variant).
    level_indices : list[np.ndarray]
        For each level :math:`p`, the indices into the flattened node
        array of ordinates on that level, ordered by the **fiber's own
        coordinate**: primarily by increasing :math:`\eta = \mu_x` (the
        radial cosine — the cylindrical sweep convention from Bailey
        et al. 2009 Eq. 50), ties broken by increasing :math:`\varphi`
        then increasing :math:`\operatorname{sign}(\mu_z)`.

        :meth:`from_level_membership` is the constructor that
        establishes that order, and it REFUSES a level on which the key
        repeats — so the order is a property of the rule, with no sort
        algorithm contributing. Until 2026-08-13 the key was
        :math:`\eta` ALONE, which the warning below explains cannot
        order a full rule's fiber; the two producers were then obliged
        to break the resulting ties, and did so differently
        (``product_mu_phi`` named its tie-break ``kind="stable"``;
        ``level_symmetric_sn`` left it to :func:`numpy.argsort`'s
        introsort partition).

        .. warning::

           :math:`\eta` ALONE is **2-to-1 on the fiber** of a **full**
           rule — the reason the key needs its other two components.
           A level is a circle, and
           :math:`\eta = \sin\theta\cos\varphi` is even in
           :math:`\varphi`, so ordering by it is an ordering of the
           circle *modulo the mirror* rather than of the circle.
           `[M]` On one level of ``product(2, 8)``, 8 ordinates give
           only 5 distinct :math:`\eta`, and the resulting order does
           not determine :math:`\operatorname{sign}(\xi)` — that
           information is not in the key. This is the mechanism
           behind issue #326.

           On a **folded** rule (:meth:`quotient` of a
           :math:`\sigma_y` mirror fold) the disease is gone by
           construction: a level is a single ARC, :math:`\eta` is
           injective on it, and the stored order IS an ordering of
           the fiber — the arc traversed in strictly *decreasing*
           :math:`\omega`, which is the azimuthal march order (T22b:
           one order seen through two charts; :math:`\eta`-ascending
           and :math:`\omega`-descending are the same traversal
           because :math:`\eta = \sin\theta\cos\omega` is strictly
           decreasing in :math:`\omega` on the surviving arc). The
           campaign's second accessor ``fiber()`` — an
           :math:`\omega`-ascending ordering kept while #326 was
           being adjudicated — retired on this theorem: on the object
           the sweep will consume, the two orders differ only by
           orientation, and the stored one is the march.
    level_mu : np.ndarray
        The invariant's value per level — read
        :attr:`polar_invariant` to know *which* invariant.
    polar_invariant : PolarInvariant
        Which invariant the levels fiber over. Not decoration: it
        decides whether level :math:`p` holds one hemisphere's circle
        or both.
    azimuth : np.ndarray
        Per-ordinate azimuthal angle :math:`\varphi \in [0, 2\pi)`.
    hemisphere : np.ndarray
        Per-ordinate :math:`\operatorname{sign}(\mu_z)`, as ``-1`` /
        ``0`` / ``+1``.

        Together with :attr:`azimuth` this is **the fiber's own
        coordinate**, and it takes both because the fiber is not always
        one circle: under :attr:`PolarInvariant.ABS_MU_Z` a level
        carries BOTH hemispheres, so :math:`\varphi` alone repeats. The
        pair is injective on a level under either invariant, which is
        exactly what an ordering of the fiber needs.
    """

    n_levels: int
    level_indices: list[np.ndarray]
    level_mu: np.ndarray
    polar_invariant: PolarInvariant
    azimuth: np.ndarray
    hemisphere: np.ndarray

    @classmethod
    def from_level_membership(
        cls,
        membership: Sequence[np.ndarray],
        *,
        nodes: np.ndarray,
        level_mu: np.ndarray,
        polar_invariant: PolarInvariant,
        azimuth: np.ndarray,
        hemisphere: np.ndarray,
    ) -> LevelStructure:
        r"""Build from unordered level MEMBERSHIP; the type orders it.

        The producers of a levelled rule know which ordinates share a
        level — that is index arithmetic on their own construction. They
        do NOT get to decide what "the order of a level" means, because
        that is a property of the fiber, and the fiber's chart lives
        here (:attr:`azimuth`, :attr:`hemisphere`). So ``membership`` is
        a **partition** — its intra-level order is read and discarded —
        and the canonical order is established once, in this body.

        The order is the fiber's own coordinate, lexicographically:

        1. :math:`\eta = \mu_x` **ascending** — the cylindrical sweep
           convention (Bailey et al. 2009 Eq. 50), and the only
           component any consumer's arithmetic reads;
        2. :math:`\varphi` **ascending** — the fiber angle;
        3. :math:`\operatorname{sign}(\mu_z)` **ascending** — which
           circle of the fiber, needed only under
           :attr:`PolarInvariant.ABS_MU_Z`, where a level carries two.

        Components 2-3 are the pair the class docstring names as
        injective on a level, so **the key admits no tie** and no sort
        algorithm can contribute to the answer. That is the whole point:
        this replaces two producers that sorted on :math:`\eta` alone
        and therefore *had* to break ties somehow — the product rule
        named its tie-break (``kind="stable"``), the level-symmetric
        rule left it to :func:`numpy.argsort`'s introsort partition,
        which is neither documented nor stable across array sizes
        (`[M]` numpy falls to insertion sort at :math:`\le 16` elements,
        so the two agreed on small levels and diverged from 24 up).

        Injectivity is **certified, not assumed** — see
        :meth:`quotient` for the same discipline on the fiberwise
        precondition. A degenerate key would silently reinstate exactly
        the accident this constructor exists to remove, so it is refused
        rather than tolerated. Spelled as a ``raise`` because production
        ``assert`` vanishes under ``-O``, which is this project's
        canonical test invocation.

        Parameters
        ----------
        membership
            One index array per level, in level order. Treated as a
            SET per level; any incoming order is discarded.
        nodes
            The companion measure's ``(N, 3)`` direction cosines. Only
            column 0 (:math:`\eta`) is read — the primary sort key.
        level_mu, polar_invariant, azimuth, hemisphere
            Passed through to the corresponding attributes.

        Raises
        ------
        ValueError
            If the fiber key is not injective on some level, i.e. two
            ordinates of one level share :math:`(\eta, \varphi,
            \operatorname{sign}\mu_z)` — then "the order of that level"
            is not defined by the rule.
        """
        eta = nodes[:, 0]
        ordered: list[np.ndarray] = []
        for p, members in enumerate(membership):
            level = np.asarray(members)
            key = np.stack([eta[level], azimuth[level], hemisphere[level]])
            distinct_fiber_points = int(np.unique(key, axis=1).shape[1])
            if distinct_fiber_points != level.size:
                raise ValueError(
                    f"the fiber key is not injective on level {p}: "
                    f"{level.size} ordinates share only "
                    f"{distinct_fiber_points} distinct "
                    f"(eta, phi, sign mu_z) triples, so the level has no "
                    f"order the rule determines. A levelled rule must "
                    f"place distinct ordinates at distinct points of its "
                    f"fiber; duplicated nodes are the usual cause "
                    f"(cf. ERR-073)."
                )
            ordered.append(level[np.lexsort((key[2], key[1], key[0]))])

        return cls(
            # DERIVED, never passed: a structure whose n_levels
            # disagrees with its own level count is unspellable through
            # this constructor.
            n_levels=len(ordered),
            level_indices=ordered,
            level_mu=level_mu,
            polar_invariant=polar_invariant,
            azimuth=azimuth,
            hemisphere=hemisphere,
        )

    def quotient(
        self,
        *,
        parent: DiscreteMeasure,
        onto: DiscreteMeasure,
        mass_rtol: float = 1e-12,
    ) -> LevelStructure:
        r"""Descend this structure along a quotient ``parent`` has taken.

        The folded rule's structure, derived by **selection** — the
        companion of :meth:`DiscreteMeasure.quotient
        <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
        folds the measure but cannot know about levels. A quotient
        never moves a node: every atom of the folded measure is an
        orbit representative, a bit-copy of one of ``parent``'s nodes.
        So every per-ordinate field descends by pure selection, and the
        level order descends by **restriction**:

        * ``level_indices[p]`` — the parent's :math:`\eta`-sorted level,
          restricted to the surviving representatives. A subsequence of
          a sorted sequence is sorted, so the sort convention is
          spelled ONCE (in the producer) and inherited here, never
          re-derived.
        * ``azimuth`` / ``hemisphere`` — selections of the parent's
          charts, bit-identical. Nothing is recomputed, so nothing can
          drift.
        * ``n_levels`` / ``level_mu`` / ``polar_invariant`` — unchanged,
          because the descent is defined exactly when the group action
          is **fiberwise**: it may permute a level's circle but must
          not move weight between levels.

        The fiberwise precondition is *certified, not assumed*: a
        quotient conserves mass level by level iff no orbit crosses a
        level boundary, so each level's folded mass is checked against
        its parent mass. A level-merging fold (e.g. the
        :math:`\mu`-mirror :math:`\sigma_z`, which pairs the
        :math:`\pm\mu_z` levels) is refused loudly.

        Both measures arrive keyword-only: they share a type, and a
        positional swap would be silent.

        Parameters
        ----------
        parent : DiscreteMeasure
            The rule this structure describes.
        onto : DiscreteMeasure
            The folded rule — ``parent.quotient(G)`` for a fiberwise
            ``G``.
        mass_rtol : float, optional
            Relative tolerance of the per-level mass certificate. The
            folded weights are the parent's addends re-associated
            (orbit-stabilizer sums), so the honest gap is a few ULP.

        Returns
        -------
        LevelStructure
            The folded rule's structure. On a mirror fold its levels
            are single ARCS: the :math:`\eta`-key is injective there,
            so the stored order — 2-to-1 on a full level — becomes a
            genuine ordering of the fiber, traversed in the march
            orientation (see the :attr:`level_indices` warning).
        """
        parent_nodes = np.ascontiguousarray(parent.nodes)
        parent_index: dict[bytes, int] = {}
        for i in range(parent.n_points):
            key = parent_nodes[i].tobytes()
            if key in parent_index:
                raise ValueError(
                    f"structure descent needs distinct parent nodes; "
                    f"parent atoms {parent_index[key]} and {i} coincide"
                )
            parent_index[key] = i

        onto_nodes = np.ascontiguousarray(onto.nodes)
        old_of_new = np.empty(onto.n_points, dtype=np.intp)
        for j in range(onto.n_points):
            match = parent_index.get(onto_nodes[j].tobytes())
            if match is None:
                raise ValueError(
                    f"not a selection of this structure's rule: folded "
                    f"atom {j} matches no parent node bit-for-bit — a "
                    f"quotient never moves a node, so `onto` is not a "
                    f"quotient of `parent`"
                )
            old_of_new[j] = match
        new_of_old = {int(i): j for j, i in enumerate(old_of_new)}

        folded_levels: list[np.ndarray] = []
        for p, members in enumerate(self.level_indices):
            survivors = np.array(
                [new_of_old[int(i)] for i in members if int(i) in new_of_old],
                dtype=np.intp,
            )
            parent_mass = float(parent.weights[members].sum())
            folded_mass = float(onto.weights[survivors].sum())
            if not np.isclose(folded_mass, parent_mass, rtol=mass_rtol, atol=0.0):
                raise ValueError(
                    f"the quotient does not act on the fiber: level {p} "
                    f"(invariant value {self.level_mu[p]}) has folded mass "
                    f"{folded_mass} against parent mass {parent_mass} — an "
                    f"orbit crossed a level boundary (e.g. a μ-mirror "
                    f"fold), so the level decomposition does not descend"
                )
            folded_levels.append(survivors)

        return replace(
            self,
            level_indices=folded_levels,
            azimuth=self.azimuth[old_of_new],
            hemisphere=self.hemisphere[old_of_new],
        )


def _even_monomial_conditions(max_degree: int) -> list[tuple[int, int, int]]:
    r"""The independent moment conditions up to ``max_degree``, lowest first.

    Only EVEN triples constrain anything: an :math:`O_h`-invariant node set is
    closed under sign flips, so every odd monomial sums to zero on both sides
    identically. And :math:`O_h` contains the coordinate permutations, so
    :math:`(4,2,0)` and :math:`(0,4,2)` are the SAME equation — deduping by the
    sorted triple keeps the system square rather than padded with copies.
    """
    seen: set[tuple[int, int, int]] = set()
    out: list[tuple[int, int, int]] = []
    for degree in range(0, max_degree + 1, 2):
        for a in range(0, degree + 1, 2):
            for b in range(0, degree - a + 1, 2):
                c = degree - a - b
                if c % 2:
                    continue
                key = (a, b, c)
                canonical = tuple(sorted(key))
                if canonical in seen:
                    continue
                seen.add(canonical)  # type: ignore[arg-type]
                out.append(key)
    return out


def _sphere_monomial_integral(a: int, b: int, c: int) -> float:
    r""":math:`\int_{S^2} x^a y^b z^c \,\mathrm{d}\Omega` in closed form."""
    if a % 2 or b % 2 or c % 2:
        return 0.0
    return 2.0 * (
        math.gamma((a + 1) / 2) * math.gamma((b + 1) / 2)
        * math.gamma((c + 1) / 2) / math.gamma((a + b + c + 3) / 2)
    )


def _measured_exactness_degree(
    nodes: np.ndarray, weights: np.ndarray, *, atol: float = 1e-10,
) -> int:
    r"""The rule's polynomial exactness, MEASURED against the closed form.

    The largest ``d`` with every monomial of total degree ``<= d``
    integrated to ``atol`` — the same sweep, with the same tolerance,
    that ``tests/numerics/test_advertised_degree_is_measured.py``
    re-runs independently. Used to STAMP the level-symmetric claim
    (#337): under the moment-matched seed the achieved degree is not a
    clean formula of :math:`N`, so the honest stamp is a measurement
    of the built rule, not an integer the construction hopes for.
    """
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    degree = 0
    while degree < 64:
        d = degree + 1
        for a in range(d + 1):
            for b in range(d + 1 - a):
                c = d - a - b
                quadrature_sum = float(
                    np.sum(weights * x**a * y**b * z**c)
                )
                if abs(
                    quadrature_sum - _sphere_monomial_integral(a, b, c)
                ) > atol:
                    return degree
        degree = d
    return degree


def _moment_matched_octant_weights(
    sn_order: int,
    octant_dirs: "list[tuple[float, float, float]]",
    octant_orbit: "list[tuple[int, int, int]]",
) -> "list[float]":
    r"""Carlson–Lathrop weights: one free weight per :math:`O_h` orbit, solved.

    Issue **#327**. Until 2026-08-06 this was ``4π/(8·n_octant)`` — ONE weight
    for every ordinate — and the rule was therefore **degree 3 at every order**
    while advertising :math:`N-1`, an over-claim of 12 at :math:`S_{16}`. The
    degree-3 it did reach was free: any :math:`O_h`-symmetric node set with
    :math:`\sum w = 4\pi` reaches 3, so nothing about the weights was earning it.

    :math:`O_h`-invariance FORCES the weight to be constant on each orbit, so
    the free parameters are exactly one per orbit — no more, no fewer. Solving
    the lowest independent moment conditions for them is the classical
    construction the docstring has always cited.

    ⭐ **At** :math:`S_2` **this is provably a no-op.** The node set is a
    SINGLE orbit, so invariance plus :math:`\sum w = 4\pi` determine the
    weight uniquely — bit-identical across #327 AND #337 (the pre-carve
    capture in ``tests/numerics/test_level_symmetric_nodes.py`` pins it).
    (:math:`S_4` was also a forced no-op at #327; #337 then moved its
    NODES, so its weight — still forced — changed with the orbit sums.)

    `[M]` 2026-08-08, on the moment-matched nodes (#337); the #327
    convention-seed values are struck through in the corpus table:

    ====  ======  ======  ============
    N     orbits  degree  min weight
    ====  ======  ======  ============
    2     1       3       1.570796
    4     1       5       0.523599
    6     2       7       0.246940
    8     3       9       0.142535
    10    4       11      0.070755
    12    5       11      0.040607
    14    7       15      0.012990
    16    8       15      0.016300
    18    10      17      0.000175
    ====  ======  ======  ============

    Raises
    ------
    ValueError
        When the solve yields a **non-positive** weight — measured, from
        :math:`S_{20}` upward on the moment-matched nodes (`[M]`
        ``-2.191e-4`` at :math:`S_{20}`, 50-digit-confirmed — the sign
        has ~7 orders of margin over float64, so the flip is the
        family's, not the arithmetic's).

        ⛔ **Positivity is not a preference.** :math:`\phi = \sum_n w_n \psi_n`
        must be non-negative for a non-negative angular flux, and the boundary
        response kernels asserts it outright
        (``assert_response_positive_if_declared``; the Lambertian's sub-Markov
        bound). A rule that integrates a positive flux to a negative scalar
        flux is not a transport quadrature, so the family is **defined exactly
        where it is valid and refuses elsewhere** rather than shipping a
        silent second construction under one name.

        The frontier is **computed, never hardcoded**: it is read off the
        solution's own sign, so it tracks the node set instead of going stale
        beside it. Today that puts it between :math:`S_{18}` and
        :math:`S_{20}` — and it is the END of the road, not a solver
        limitation: `[M]` :math:`S_{20}` is servable only by the
        axis-weight decomposition at full degree 11 (dominated by our
        :math:`S_{14}`), and at :math:`S_{22}` an LP over the whole
        decomposition kernel certifies NO nonnegative point weights
        exist for the even-moment family at all.
    """
    per_orbit, labels, conditions, system, targets = _greedy_orbit_solve(
        sn_order, octant_dirs, octant_orbit,
    )

    if float(np.min(per_orbit)) <= 0.0:
        raise ValueError(
            f"level-symmetric S_{sn_order}: the moment-matched construction "
            f"has no POSITIVE solution on this node set (min weight "
            f"{float(np.min(per_orbit)):.6f}). phi = sum(w*psi) must stay "
            f"non-negative for a non-negative angular flux, and the boundary "
            f"response kernels assert it, so a negative weight is not a "
            f"tolerable trade for the extra degree. The level-symmetric family "
            f"is positive up to S_18 on the moment-matched levels; above it "
            f"use Quadrature.lebedev(order) or Quadrature.product(n_mu, "
            f"n_phi), both of which reach their advertised degree with "
            f"positive weights (issues #327, #337)."
        )

    weight_of = dict(zip(labels, (float(w) for w in per_orbit)))
    return [weight_of[label] for label in octant_orbit]


def _greedy_orbit_solve(
    sn_order: int,
    octant_dirs: "list[tuple[float, float, float]]",
    octant_orbit: "list[tuple[int, int, int]]",
) -> "tuple[np.ndarray, list[tuple[int, int, int]], list[tuple[int, int, int]], np.ndarray, np.ndarray]":
    r"""The per-orbit weight solve, WITHOUT the positivity contract.

    Returns ``(per_orbit, labels, conditions, system, targets)`` — the
    solved orbit weights plus the full condition ladder they were
    solved inside, so a caller can evaluate the residual of any ladder
    row NOT in the solved set. Two callers:

    * :func:`_moment_matched_octant_weights` — adds the positivity
      raise (the production contract);
    * :func:`_axial_condition_residual` — the μ₁ root-find's objective
      (#337), which must evaluate trial seeds whose weights may be
      transiently negative or whose system may be rank-short. A trial
      failure there is a NaN to step over, not a contract violation —
      which is exactly why the raise cannot live in this core.

    Raises
    ------
    ValueError
        When fewer independent conditions than orbits exist in the
        ladder at these nodes (a rank-short greedy — degenerate level
        coincidences at pathological seeds).
    """
    orbits: "dict[tuple[int, int, int], list[int]]" = {}
    for index, label in enumerate(octant_orbit):
        orbits.setdefault(label, []).append(index)
    labels = list(orbits)

    # One octant's contribution, replicated over the 8 sign octants: the 8-fold
    # replication below multiplies every orbit sum by the same factor, so it is
    # carried here rather than left for the caller to remember.
    directions = np.asarray(octant_dirs, dtype=float)
    conditions = _even_monomial_conditions(max(sn_order + 2, 4))
    system = np.array(
        [
            [
                8.0 * float(np.sum(
                    directions[members, 0] ** a
                    * directions[members, 1] ** b
                    * directions[members, 2] ** c
                ))
                for members in (orbits[label] for label in labels)
            ]
            for (a, b, c) in conditions
        ],
        dtype=float,
    )
    targets = np.array(
        [_sphere_monomial_integral(a, b, c) for (a, b, c) in conditions],
        dtype=float,
    )

    # Take the lowest INDEPENDENT conditions, one per free weight, and solve
    # exactly. ⛔ Not least squares over all of them: with more conditions than
    # orbits that is an overdetermined compromise satisfying NONE of them —
    # `[M]` it misses even ``Σw = 4π``, the degree-0 condition, and the rule
    # then integrates a constant wrongly.
    chosen: list[int] = []
    for row in range(system.shape[0]):
        trial = chosen + [row]
        if np.linalg.matrix_rank(system[trial], tol=1e-10) == len(trial):
            chosen = trial
        if len(chosen) == len(labels):
            break
    if len(chosen) != len(labels):
        raise ValueError(
            f"level-symmetric S_{sn_order}: only {len(chosen)} independent "
            f"moment conditions found for {len(labels)} orbit weights at "
            f"this seed (a rank-short greedy — degenerate level "
            f"coincidence)."
        )
    per_orbit = np.linalg.solve(system[chosen], targets[chosen])
    return per_orbit, labels, conditions, system, targets


def _axial_condition_residual(sn_order: int, mu1_sq: float) -> float:
    r"""The :math:`\int_{S^2}\mu_z^N` defect of the trial seed ``mu1_sq``.

    The moment-matched family's DEFINING equation (#337): with weights
    solved per orbit, one condition remains for the node seed, and
    `[M]` it is always the pure axial monomial ``(0, 0, N)`` — present
    in the ladder at every order, outside the greedy's chosen set, and
    generically unsatisfied. The seed therefore has a *name*: the rule
    integrates :math:`\mu_z^N` exactly.

    Returns NaN when the trial seed is unsolvable (rank-short greedy) —
    the root-find steps over such points rather than crashing, because
    `[M]` up to 31 % of the bracket is rank-deficient at S20.
    """
    try:
        _mu_levels, octant_dirs, octant_orbit = _octant_directions(
            sn_order, mu1_sq
        )
        per_orbit, _labels, conditions, system, targets = (
            _greedy_orbit_solve(sn_order, octant_dirs, octant_orbit)
        )
    except (ValueError, np.linalg.LinAlgError):
        return float("nan")
    row = conditions.index((0, 0, sn_order))
    return float(system[row] @ per_orbit - targets[row])


@lru_cache(maxsize=None)
def _moment_matched_mu1_sq(sn_order: int) -> float:
    r"""The moment-matched seed: the SMALLEST root of the axial defect.

    Root-finds :func:`_axial_condition_residual` over
    :math:`\mu_1^2 \in (0, 1/3)` (issue #337 — the seed is no longer a
    project convention but the choice that extends the exactly-held
    moment set by :math:`\int\mu_z^N`, the construction behind the
    published LA-3186 tables: `[M]` reproduces their :math:`\mu_1` to
    every printed digit at S4/S6/S8/S12/S16).

    Two measured hazards shape the algorithm:

    * **The root is not unique** — at N = 6, 10, 14, 18 a SECOND root
      exists in (0, 1/3), and its weight solve is strongly negative
      (−0.6 … −7.7). The rule is **the smallest root**: scan LEFT to
      right and bisect the first sign change.
    * **The objective is not total** — the inner greedy is
      rank-deficient on parts of the bracket (`[M]` 31 % at S20) and
      its chosen-set switches make the residual discontinuous there.
      NaN trial points are stepped over, never bracketed across.
    """
    lo, hi, n_scan = 1e-4, 1.0 / 3.0 - 1e-4, 512
    xs = np.linspace(lo, hi, n_scan)
    prev_x = prev_f = math.nan
    for x in xs:
        f = _axial_condition_residual(sn_order, float(x))
        if not math.isfinite(f):
            prev_x = prev_f = math.nan
            continue
        if math.isfinite(prev_f) and prev_f * f < 0.0:
            a, b, fa = prev_x, float(x), prev_f
            for _ in range(200):
                mid = 0.5 * (a + b)
                fm = _axial_condition_residual(sn_order, mid)
                if not math.isfinite(fm):
                    mid = a + 0.499 * (b - a)
                    fm = _axial_condition_residual(sn_order, mid)
                    if not math.isfinite(fm):
                        break
                if fa * fm <= 0.0:
                    b = mid
                else:
                    a, fa = mid, fm
                if b - a < 1e-17:
                    break
            return 0.5 * (a + b)
        prev_x, prev_f = float(x), f
    raise ValueError(
        f"level-symmetric S_{sn_order}: the moment-matched seed equation "
        f"(the integral of mu_z^{sn_order} exact) has no root in "
        f"(0, 1/3) — the family is not realizable at this order."
    )


def _octant_directions(
    sn_order: int, mu1_sq: float,
) -> "tuple[np.ndarray, list[tuple[float, float, float]], list[tuple[int, int, int]]]":
    r"""The triangular octant point set at a given seed.

    The ONE producer of level-symmetric octant geometry, shared by the
    builder and by the μ₁ root-find's trial evaluations (#337) — so the
    trial nodes ARE the production nodes at that seed, by construction
    rather than by a second spelling.

    Levels follow Carlson–Lathrop Eq. (3-52) (LA-3251-MS printed
    p. 32): :math:`\mu_i^2 = \mu_1^2 + (i-1)\Delta`,
    :math:`\Delta = 2(1-3\mu_1^2)/(N-2)` (:math:`S_2`'s single level
    is the seed itself — no recursion, no division by :math:`N-2`).
    """
    n_half = sn_order // 2
    if n_half == 1:
        mu_levels = np.sqrt(np.array([mu1_sq]))
    else:
        delta = 2.0 * (1.0 - 3.0 * mu1_sq) / (sn_order - 2)
        mu_levels = np.sqrt(mu1_sq + np.arange(n_half) * delta)

    # The defining property of a level-symmetric set: every ordinate's
    # three direction cosines are drawn from the SAME level array, so
    # the set is closed under permuting the axes. The level index of the
    # third cosine is not free — it is fixed by the other two.
    #
    #   mu2[p] + mu2[k] + mu2[j] = 3*mu1_sq + (p + k + j)*delta = 1
    #   and delta = 2(1 - 3*mu1_sq)/(N - 2)
    #   =>  p + k + j = (1 - 3*mu1_sq)/delta = (N - 2)/2 = n_half - 1.
    #
    # So j is INDEX ARITHMETIC. Until 2026-08-02 this loop instead
    # recovered the third cosine numerically as
    # ``sqrt(1 - mu_z**2 - eta**2)``, which lands within ~1e-16 of
    # ``mu_levels[j]`` but not ON it: `[M]` at N=16 the y-axis then
    # carried 22 distinct magnitudes where the level array has 8, and
    # only 8 of the 48 O_h operators mapped the node set onto itself
    # bit-exactly — the 8 pure sign flips, since negation is exact in
    # IEEE while a coordinate swap is only exact if the values match.
    # The rule advertised O_h invariance and realized D_2h.
    #
    # Reading the level value instead makes all 48 exact, so the
    # octahedral symmetry becomes an integer permutation of ordinate
    # indices rather than a question about tolerances. It also retires
    # the `xi_sq < -1e-14` guard: `j` is provably in range for every
    # admissible (p, k), because p + k <= n_half - 1 by the loop bound.
    octant_dirs: list[tuple[float, float, float]] = []
    # ⭐ The O_h ORBIT of each octant direction, as EXACT index arithmetic.
    #
    # O_h contains the coordinate permutations and all sign flips, so two
    # directions lie in the same orbit iff their |cosine| MULTISETS agree —
    # and every cosine here is ``mu_levels[·]``, so the multiset is the sorted
    # triple of LEVEL INDICES. That makes the orbit label an integer fact of
    # the construction rather than a float comparison discovered afterwards
    # (the same reasoning that made ``j`` index arithmetic above: a symmetry
    # question answered by a loop that already knows the answer exactly).
    octant_orbit: list[tuple[int, int, int]] = []
    for p in range(n_half):
        mu_z = mu_levels[p]
        for k in range(n_half - p):
            j = n_half - 1 - p - k
            octant_dirs.append((mu_levels[k], mu_levels[j], mu_z))
            octant_orbit.append(tuple(sorted((p, k, j))))  # type: ignore[arg-type]
    return mu_levels, octant_dirs, octant_orbit


def _build_level_symmetric_arrays(
    sn_order: int,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
    list[np.ndarray], np.ndarray, np.ndarray,
]:
    r"""Construct level-symmetric :math:`S_N` quadrature arrays.

    This body was lifted byte-for-byte from the private
    ``_build_level_symmetric`` of the now-retired
    ``orpheus/sn/quadrature.py`` (R-1 Phase A detour-C), and is today
    the ONE producer of level-symmetric ordinates in the tree. The
    old side-by-side bit-identity comparison is therefore no longer
    expressible, and **nothing witnesses that the carve preserved the
    node order** — the honest statement, which the companion gate
    ``tests/numerics/test_rules_sphere.py::
    test_level_symmetric_bit_identical_to_legacy_adapter`` already
    carries in its own docstring. What exists going forward is a
    POST-carve floor: the frozen literal in
    ``tests/numerics/test_level_symmetric_nodes.py`` (nodes at S2, the
    level partition at S4..S12).

    ⛔ Until 2026-08-13 this paragraph named
    ``tests/sn/regression/snapshots/cyl_*_LS4_*.npz`` as the pin. That
    was false twice over. Those three files were **deleted** at
    ``c39b7d44`` when the cylinder snapshot families were re-captured on
    the σ_y fold (Q5.6.3) — the survivors are ``cyl_*_folded_*``. And
    the mechanism it presumed is refused: a level-symmetric rule cannot
    reach the cylindrical sweep at all, because an ``ABS_MU_Z`` level
    carries both hemispheres, so :math:`\eta` is 4-fold degenerate and
    ``assert_carrying_quadrature`` refuses on ``eta_level[0] ==
    eta_level[1]`` (`[M]` at every order S2..S12, on the spherical arm
    too via :math:`\tau_{\rm raw} \notin [0, 1]`; #336 tracks the
    refuse-or-reduce design). The claim's correction had already landed
    on the TEST-side twin and never reached this producer.

    Returns
    -------
    mu_x, mu_y, mu_z, weights : np.ndarray
        Direction cosines and weights, shape ``(N,)`` each.
    level_mu : np.ndarray
        Polar cosine per level.
    level_membership : list[np.ndarray]
        Per-level ordinate indices as a PARTITION — unordered. The
        canonical intra-level order is established by
        :meth:`LevelStructure.from_level_membership`, which is where
        the fiber's chart lives; this function must not pre-empt it.
        The level COUNT is not returned either: it is ``len()`` of this,
        and a second spelling of it is a second thing to keep in sync.
    azimuth, hemisphere : np.ndarray
        The fiber coordinate (Q4): :math:`\varphi \in [0, 2\pi)` and
        :math:`\operatorname{sign}(\mu_z)` per ordinate. The annotation
        omitted these two from 2026-08-02 (`3afb52c2`) until the exactness
        carve, while the body returned all nine — a stale signature the
        Sphinx build could not see and the tests could not fail on.
    """
    if sn_order % 2 != 0 or sn_order < 2:
        raise ValueError(f"S_N order must be positive even, got {sn_order}")

    n_half = sn_order // 2

    if n_half == 1:
        # S_2 has no freedom: mu^2 = 1/3 is forced by the diffusion
        # condition (LA-3251-MS printed p. 45).
        mu1_sq = 1.0 / 3.0
    else:
        # The recursion is Carlson-Lathrop Eq. (3-52) (LA-3251-MS
        # printed p. 32), and since #337 the SEED is MOMENT-MATCHED:
        # the smallest root of "the rule integrates mu_z^N exactly" —
        # the construction behind the published LA-3186 tables (`[M]`
        # reproduces their mu_1 to every printed digit at
        # S4/S6/S8/S12/S16). Until 2026-08-08 this line read the
        # project convention 4/(N(N+2)), which limited the family by
        # the seed alone (degree ceiling N-1; positivity frontier
        # S_12); the corpus provenance note in angular_quadrature.rst
        # carries the verification against the primary sources
        # (#327, #337).
        mu1_sq = _moment_matched_mu1_sq(sn_order)

    mu_levels, octant_dirs, octant_orbit = _octant_directions(
        sn_order, mu1_sq
    )

    w_octant_per_dir = _moment_matched_octant_weights(
        sn_order, octant_dirs, octant_orbit,
    )

    all_eta: list[float] = []
    all_xi: list[float] = []
    all_mu: list[float] = []
    all_w: list[float] = []
    for (eta, xi, mu_z), w_dir in zip(octant_dirs, w_octant_per_dir):
        for s_eta in (-1, 1):
            for s_xi in (-1, 1):
                for s_mu in (-1, 1):
                    all_eta.append(s_eta * eta)
                    all_xi.append(s_xi * xi)
                    all_mu.append(s_mu * mu_z)
                    all_w.append(w_dir)

    mu_x = np.array(all_eta)
    mu_y = np.array(all_xi)
    mu_z_arr = np.array(all_mu)
    weights = np.array(all_w)

    n_levels = n_half
    level_mu_vals = mu_levels
    level_membership: list[np.ndarray] = []
    for p in range(n_levels):
        # EXACT membership. The 8-fold sign replication above copies
        # mu_z straight out of `mu_levels`, so |mu_z| IS the level value
        # bit-for-bit and the comparison is an equality, not a
        # neighbourhood test. This carried `tol = 1e-12` until
        # 2026-08-02 — a symmetry question answered by float comparison
        # in a loop that already knew the answer exactly.
        #
        # MEMBERSHIP ONLY: this loop answers "which ordinates share a
        # level", which is index arithmetic on the replication above.
        # It does NOT order them — that is a property of the fiber, and
        # `LevelStructure.from_level_membership` owns it. Until
        # 2026-08-13 the next line was a bare `np.argsort(mu_x[idx])`,
        # and eta is 4-fold degenerate on an ABS_MU_Z level (the ±xi,
        # ±mu_z replications), so the intra-level order was introsort's
        # partition detail — `[M]` reproduced verbatim in the frozen
        # literal of `tests/numerics/test_level_symmetric_nodes.py`,
        # where S6 level 0 read `..., 0, 3, 2, 1, 6, 5, 4, 7, ...`.
        level_membership.append(
            np.where(np.abs(mu_z_arr) == level_mu_vals[p])[0]
        )

    # The fiber coordinate: azimuth about the polar axis. Distinct
    # ordinates on a level sit at distinct phi, which is what makes an
    # ordering by it a function OF the fiber (unlike mu_x, which is even
    # in phi and therefore 2-to-1 on it).
    azimuth = np.mod(np.arctan2(mu_y, mu_x), 2.0 * np.pi)
    hemisphere = np.sign(mu_z_arr).astype(np.int64)

    return (mu_x, mu_y, mu_z_arr, weights, level_mu_vals,
            level_membership, azimuth, hemisphere)


def level_symmetric_sn(
    sn_order: int,
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Carlson-Lathrop level-symmetric :math:`S_N` rule on :math:`S^2`.

    Standard triangular discrete-ordinate quadrature with :math:`N/2`
    polar levels per hemisphere. The construction is
    :math:`O_h`-invariant: every octant carries the same set of
    direction-cosine triples up to sign permutations, with the weight
    constant on each :math:`O_h` orbit (one free weight per orbit,
    solved — see :func:`_moment_matched_octant_weights`; #327).

    Weight sum: :math:`\sum_i w_i = 4\pi`.

    Symmetry: :math:`O_h` — invariant under all 48 rotation /
    reflection elements of the octahedral group.

    Polynomial exactness: BUILD-MEASURED (`[M]` :math:`3` at
    :math:`S_2`; :math:`N+1` at :math:`S_4`–:math:`S_{10}` and
    :math:`S_{14}`; :math:`N-1` at :math:`S_{12}`/:math:`S_{16}`/
    :math:`S_{18}` — no clean formula of :math:`N`), gated three ways:
    the stamp, an independent sweep, and frozen 50-digit-checked
    literals (``tests/numerics/test_level_symmetric_nodes.py`` +
    ``test_advertised_degree_is_measured.py``).

    The node seed is **moment-matched** (#337, 2026-08-08): the
    smallest root of "the rule integrates :math:`\mu_z^N` exactly",
    which reproduces the published LA-3186 :math:`\mu_1` to every
    printed digit at each tabulated order (and beats the print at
    three — one-ulp last-digit slips documented in the corpus). Our
    POINT weights, however, are the per-orbit cross-moment solve, NOT
    LA-3186's axis-weight decomposition (:math:`p = a_i + a_j + a_k`):
    `[M]` on the same nodes the published decomposition reaches full
    3-D degree 11 at every order :math:`\geq 14` while this solve
    reaches 15/15/17 at 14/16/18 — and the two frontiers INTERLEAVE
    (theirs dies at 18, serves 20 at degree 11; ours serves 18 and
    dies at 20). :math:`S_{22}` is intrinsically dead for the whole
    even-moment family (LP certificate: no nonnegative decomposition
    exists). The corpus note in
    ``docs/theory/methods/sn/angular_quadrature.rst`` carries the
    full provenance and the comparison table.

    Returns the :class:`DiscreteMeasure` and a
    :class:`LevelStructure` capturing the per-level grouping needed by
    the cylindrical SN sweep.

    Parameters
    ----------
    sn_order : int
        Even :math:`N \ge 2`, supported through :math:`S_{18}`.
        Common values: 4, 8, 12. (Read ``S_{12}`` here until 2026-08-13
        — the pre-#337 convention seed's frontier, superseded by the
        moment-matched one; see ``Raises``.)

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(N_total, 3)``, weights shape ``(N_total,)``,
        ``invariance_group=SubgroupOfO3.OctahedralOh``. The
        ``degree_of_exactness`` is **build-measured**, not a formula of
        :math:`N` — see the stamp's own comment in the body.

        ⛔ Read ``max(3, sn_order-1)`` here until 2026-08-13. Same
        staleness as the frontier above and from the same cause: the
        formula was true for the pre-#337 convention seed, and under the
        moment-matched one it UNDER-claims by 2 wherever it is wrong.
        `[M]` measured / formula: S4 **5**/3, S6 **7**/5, S8 **9**/7,
        S14 **15**/13 — and it happens to agree at S2, S12, S16, S18, so
        a spot-check on the wrong order confirms it. The body has said
        the degree is build-measured since #337; this block did not.
    LevelStructure
        Per-level indexing metadata.

    Raises
    ------
    ValueError
        Above :math:`S_{18}` — the per-orbit weight solve has no
        POSITIVE solution on this node seed from :math:`S_{20}`, and
        positivity is not tradeable. `[M]` 2026-08-13:
        :math:`S_{14}/S_{16}/S_{18}` build with smallest weight
        ``0.012990`` / ``0.016300`` / ``1.75e-4``;
        :math:`S_{20}`/:math:`S_{22}` refuse. The frontier is computed
        from the solution's own sign, never hardcoded; the full doctrine
        lives on :func:`_moment_matched_octant_weights`.

        ⛔ This block read "Above :math:`S_{12}` … a negative weight …
        from :math:`S_{14}` (`[M]` ``-0.027``)" until 2026-08-13. That
        measurement was correct **for the pre-#337 convention seed**
        :math:`\mu_1^2 = 4/(N(N+2))`, whose positivity frontier really
        was :math:`S_{12}`; #337 replaced the seed with the
        moment-matched root and moved the frontier to :math:`S_{18}`,
        without updating this docstring. The number was never wrong —
        its CONFIGURATION stopped being the shipped one.

    See Also
    --------
    :meth:`orpheus.numerics.quadrature.Quadrature.level_symmetric` —
    the named factory SN consumers call. It wraps this measure and
    holds the :class:`LevelStructure` on its ``level_structure`` field
    for hot-path access. There is no per-family adapter class: the
    SN-side wrapper this docstring used to point at
    (``orpheus.sn.quadrature.LevelSymmetricSN``) was retired into a
    classmethod factory on the one ``Quadrature`` type.
    """
    (mu_x, mu_y, mu_z, w, level_mu, level_membership,
     azimuth, hemisphere) = _build_level_symmetric_arrays(sn_order)
    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=w,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.OctahedralOh,
        # ⭐ BUILD-MEASURED since #337: the achieved degree is not a
        # clean formula of N under the moment-matched seed (`[M]` N+1
        # at 4/6/8/10/14, N-1 at 12/16/18, 3 at S_2), so the stamp is
        # what THIS rule measures against the closed-form monomial
        # integral — the same sweep the both-directions gate
        # (``tests/numerics/test_advertised_degree_is_measured.py``)
        # re-runs independently, with the gate's third corner a table
        # of FROZEN literals cross-checked at 50 digits (so the stamp,
        # the gate's sweep, and the literals must all three agree).
        #
        # History: ``degree=sn_order - 1`` unconditionally until #327
        # (FALSE both ways — equal weights delivered 3 at every order
        # while S_2 over-delivered); ``max(3, sn_order-1)`` until #337
        # (true for the convention seed, an under-claim of 2 at
        # S_4..S_10 under the moment-matched one).
        exactness=ExactnessClaim(
            reference=UNIFORM_ON_SPHERE,
            degree=_measured_exactness_degree(nodes, w),
        ),
    )
    structure = LevelStructure.from_level_membership(
        level_membership,
        nodes=nodes,
        level_mu=level_mu,
        # ABS_MU_Z: levels are indexed by |mu_z|, so each carries BOTH
        # hemispheres and its fiber is two circles. That is exactly why
        # `sign(mu_z)` is load-bearing in the ordering key here and inert
        # for `product_mu_phi` — `[M]` every eta-tie on a level of
        # S2..S12 holds two phi-collisions that only the hemisphere
        # separates, against zero on the product side.
        polar_invariant=PolarInvariant.ABS_MU_Z,
        azimuth=azimuth,
        hemisphere=hemisphere,
    )
    return measure, structure


__all__ = [
    "LevelStructure",
    "lebedev_sphere",
    "level_symmetric_sn",
]
