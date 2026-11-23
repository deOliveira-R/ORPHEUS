r"""Tagged quadrature registry: the catalogue of rules and the geometry-to-symmetry table.

The module holds two catalogues. :data:`quadrature_registry` lists every
rule the tree can build, one :class:`QuadratureSpec` per family: the
:class:`~orpheus.numerics.quadrature.directional.Quadrature` factory, the
inversion of the rule's degree-of-exactness formula
(:attr:`QuadratureSpec.degree_of_exactness_for`), and the two facts a
parameter-free field can state truthfully (``axis_aligned``,
``half_range_clean``). :data:`GEOMETRY_ANGULAR_SYMMETRY` assigns each
geometry its :class:`AngularSymmetry`: the angular domain the dimensional
reduction left behind, the symmetry its solution keeps, and the closure a
reflecting face still owes. A rule is admitted for a geometry by the
domain leg (:meth:`AngularSymmetry.admits_domain`,
:meth:`AngularSymmetry.domain_refusal`) and the symmetry leg
(:meth:`AngularSymmetry.admits_symmetry`); the symmetry of a rule is
**computed from its nodes**
(:meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`),
never read from a declared tag, and its exactness claim is read off the
built measure
(:class:`~orpheus.numerics.exactness.ExactnessClaim`), never declared on
the spec. A caller names its rule through the
:class:`~orpheus.numerics.quadrature.directional.Quadrature` factories
with an explicit parameter, which is a modelling choice the user owns.

Geometry → angular-symmetry assignment
--------------------------------------

Each geometry declares three facts about its angular symmetry (see
:class:`AngularSymmetry`): the stabiliser :math:`K` the dimensional
reduction **spends**, the finite symmetry :math:`\Gamma` the solution
keeps **unspent** in the local frame (what a fold may spend), and the
closure :math:`R` still **owed** to a reflecting face.

.. list-table::
   :header-rows: 1
   :widths: 14 10 10 10 10 46

   * - Geometry
     - :math:`K` (spent)
     - Domain :math:`S^2/K`
     - :math:`\Gamma` (unspent)
     - :math:`R` (owed)
     - Rationale
   * - ``"slab"``
     - :math:`O(2)_x`
     - :math:`S^2/O(2)_x` (chart :math:`[-1,1]`)
     - :math:`\{e\}`
     - :math:`\sigma_x`
     - 1-D in :math:`x`, the slab normal and the tree's polar axis
       (tracker 2.4). Azimuthal rotation about the slab normal — and every
       reflection in a plane containing it — is integrated out
       analytically, so the angular variable is :math:`\mu = \cos\theta`
       alone. Nothing finite is left unspent: a slab solution is not
       even in :math:`\mu` for a general source. What remains owed is
       :math:`\mu \to -\mu`, the reflection that pairs the two sweep
       senses and that a reflecting end face consumes. Gauss-Legendre
       nodes are symmetric (Stoer-Bulirsch §3.6), so it holds.
   * - ``"sphere"``
     - :math:`O(2)_x`
     - :math:`S^2/O(2)_x` (chart :math:`[-1,1]`)
     - :math:`\{e\}`
     - :math:`\sigma_x`
     - The 1-D **radial** spherical SN reduces to GL on
       :math:`\mu_r = \cos\theta_r`, the cosine of the angle between
       the ordinate and the radial direction (Lewis & Miller 1993
       §4.4). The continuous problem has :math:`O(3)` symmetry; the
       radial reduction spends the stabiliser of :math:`\hat r` and
       leaves the same owed :math:`\sigma_r` as the slab. Here the spent
       half is not free — its fiber action reappears in the sweep as the
       angular-redistribution :math:`\alpha` term.
   * - ``"cylinder"``
     - trivial
     - :math:`S^2`
     - :math:`D_{1h} = \{e, \sigma_y, \sigma_z, C_2^x\}`
     - :math:`D_{2h}`
     - An axisymmetric cylinder is :math:`\phi`-independent in
       *space*, but that does not reduce the *angular* domain: both
       angular degrees of freedom survive, so the rule lives on
       :math:`S^2` or on a fold of it. In the local frame
       (:math:`\eta = \Omega\cdot\hat r`, :math:`\xi = \Omega\cdot\hat\phi`,
       :math:`\mu = \Omega\cdot\hat z` — the tree's columns 0, 1, 2) the
       solution IS even under :math:`\sigma_y` (:math:`\xi \to -\xi`, the
       azimuthal cosine: the problem is invariant under :math:`\phi \to
       -\phi`, which pairs :math:`+\xi` with :math:`-\xi`) and under
       :math:`\sigma_z` (:math:`\mu \to -\mu`, the axial cosine: the
       geometry is z-uniform, :math:`\partial_z` is absent, and :math:`\mu`
       enters the streaming operator only through :math:`|\eta| =
       \sqrt{1-\mu^2}`; a scattering kernel in :math:`\Omega\cdot\Omega'` is
       invariant under any joint orthogonal map) — the group they generate is
       :math:`D_{1h}` about the radial axis (`[M]` 2026-09-03 the
       realization's four elements). The shipped ``folded_product`` spends
       :math:`\sigma_y`. Owed is :math:`D_{2h}`, the coordinate-plane
       mirrors. The cylindrical SN sweep additionally requires
       per-:math:`\mu` polar-level structure for the azimuthal
       redistribution coefficients; the rule carries it as its
       ``level_structure`` (:class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`).
   * - ``"cartesian2d"``
     - trivial
     - :math:`S^2`
     - :math:`\sigma_z`
     - :math:`D_{2h}`
     - 2-D Cartesian (x-y), z-uniform: :math:`\psi` is even in
       :math:`\mu_z`, and in NOTHING else — the :math:`\sigma_y` fold that
       is legitimate on the cylinder is not licensed here (`[M]` until R3
       of #434 it was admitted at both stages, and 2 of the 4
       :math:`(\operatorname{sign}\mu_x, \operatorname{sign}\mu_y)` sweep
       quadrants were then empty). :math:`D_{2h} \cong (\mathbb{Z}_2)^3`
       is generated by the three coordinate-plane mirrors, and its
       chambers are exactly the octants the sweep decomposes into —
       which is precisely the symmetry a reflecting :math:`x` or
       :math:`y` face needs to be representable exactly.

.. note::

   Before 2026-08-02 this table held a single group per geometry and
   recorded the **spent** half in it: ``slab``/``sphere``/``cylinder``
   all read :math:`SO(2)`. That is a true statement *about the
   geometry* and a useless one *for selecting a quadrature*, because
   no finite point set on :math:`S^2` is :math:`SO(2)`-closed — the
   gate :math:`G_{\text{geom}} \subseteq G_Q` was unsatisfiable by any
   discrete azimuthal rule and could only ever pass on a false tag.
   ``cartesian2d`` read :math:`O_h`, a 6× over-claim demanding the
   :math:`x \leftrightarrow z` exchange, never a symmetry of a
   z-uniform problem. Splitting spent from owed is what makes the gate
   both satisfiable and discriminating — and splitting UNSPENT from owed
   (R3 of #434, 2026-09-03) is what stops the owed closure being read as
   a fold licence.

When 2-D / 3-D spherical SN lands, this table gains ``"sphere2d"`` /
``"sphere3d"`` entries; both spend nothing continuously, keep nothing
unspent, and owe at least :math:`D_{2h}`.

References
----------

* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the method
  of discrete ordinates." In *Computing Methods in Reactor Physics*,
  Greenspan, Kelber, Okrent, eds., Gordon & Breach. Level-symmetric
  :math:`S_N` quadratures.
* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Computational
  Mathematics and Mathematical Physics* **16**(2), 10-24. The
  octahedral-invariant sphere quadrature.
* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. Theorem 3.6.20 (Gauss-Legendre
  exactness :math:`\le 2n - 1`).
* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. Symmetry-group framework underlying the
  G-compatibility filter.
* Lewis, E.E. and Miller, W.F. (1993). *Computational Methods of
  Neutron Transport*. Wiley. §4.2 (level-symmetric construction);
  §4.4 (1-D spherical SN).
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010). "The Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes."
  *Nucl. Sci. Eng.* **165**(2), 149-169, doi:10.13182/NSE08-66.
  **Eq. (50)** is the cylindrical (R-Z) :math:`\alpha`-recursion whose
  per-level march is what needs ``level_structured``; **Eq. (52)** is
  the per-level edge-cosine accumulation, i.e. the statement that a
  level's ordinates run ascending in the radial cosine — which is the
  ORDER the ``LevelStructure`` index lists carry.

  ⛔ *Retracted citation (2026-08-27).* This entry read "Bailey,
  T.S., Adams, M.L., Yang, B., Zika, M.R. (2009). 'A piecewise linear
  discontinuous finite element spatial discretization of the transport
  equation.' Annals of Nuclear Energy 35, 1929-1936." **No such
  publication exists** - its title, author list, journal reference and
  year each trace to a different real source. The equation numbers were
  right and belong to BMC 2010 above. Full account, the BMC
  equation-number map, the ⚠ published-typo warning on Eq. (50) and the
  scoping note on Eq. (52): ``docs/theory/methods/sn/curvilinear_one_group.rst`` §bmc-equation-map.

"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Any, Callable

from ..exactness import UNIFORM_ON_SPHERE, ReferenceMeasure
from ..generating_measure import pushforward_reference
from ..measure import DiscreteMeasure
from orpheus.geometry.manifold import SPHERE, Manifold, Quotient, quotient_onto
from orpheus.geometry.symmetry import SubgroupOfO3
from .directional import Quadrature


# ---------------------------------------------------------------------------
# Available Lebedev orders — ASKED of SciPy, never tabulated
# ---------------------------------------------------------------------------

#: Upper bound of the discovery sweep in :func:`_lebedev_orders`. A **search
#: window**, not a claim about what SciPy serves — the difference matters,
#: because a claim can be wrong and a window can only be too small, which
#: :func:`_lebedev_orders` refuses rather than silently truncating. Generous
#: because it is free: `[M]` 2026-08-14 sweeping to 1001 costs **16.6 ms**,
#: identical to sweeping to 201, since an unavailable order raises
#: immediately without constructing anything.
_LEBEDEV_SEARCH_CEILING = 1001


@cache
def _lebedev_orders() -> tuple[int, ...]:
    r"""The Lebedev orders the **installed** SciPy actually serves.

    ``scipy.integrate.lebedev_rule`` serves a fixed, gappy set of odd orders
    and raises ``NotImplementedError`` outside it, with no public accessor
    for the set (`[M]` ``scipy.integrate._lebedev`` exposes only functions;
    the orders live inside their bodies). So the set is discovered the only
    way SciPy offers: by asking.

    ⭐ This is the same ruling :func:`_ls_sn_invert` already states for the
    level-symmetric frontier — *ask the family, never tabulate* — applied to
    the family next door, and it is applied here because the tabulated twin
    had already gone stale exactly as that ruling predicts.

    ⛔ Until 2026-08-14 this was a literal 18-tuple topping at **47**, whose
    own comment read *"probed at registry-construction time on the SciPy
    installed in this dev container"*. That probe is ``60f9fb29``,
    **2026-05-06**, and was never re-run. `[M]` the installed SciPy 1.17.1
    serves **32** orders topping at **131**, so the selector (retired 2026-11-20) refused degrees
    the tree could deliver: ``_lebedev_invert(53)`` returned ``None`` while
    ``lebedev_sphere(53)`` builds 974 nodes summing to :math:`4\pi`. A
    frozen probe is a measurement whose configuration (the installed SciPy)
    can change underneath it without anything going red.

    Cost is paid at most once per process, and only if an inverter is called
    (`[M]` the retired selector had zero production consumers). `[M]` 16.6 ms.

    Raises
    ------
    RuntimeError
        If the largest order found equals :data:`_LEBEDEV_SEARCH_CEILING` —
        the window is then binding and the answer would be a silent
        truncation, which is the failure this function exists to remove.
    """
    orders = tuple(
        n for n in range(3, _LEBEDEV_SEARCH_CEILING + 1, 2)
        if _lebedev_order_is_available(n)
    )
    if orders and orders[-1] == _LEBEDEV_SEARCH_CEILING:
        raise RuntimeError(
            f"the Lebedev discovery sweep hit its ceiling "
            f"({_LEBEDEV_SEARCH_CEILING}): SciPy serves that order, so the "
            f"window is binding and larger orders may exist but go unseen. "
            f"Raise _LEBEDEV_SEARCH_CEILING — it is a search bound, not a "
            f"claim, and widening it is free."
        )
    return orders


def _lebedev_order_is_available(order: int) -> bool:
    """Does the installed SciPy serve this Lebedev order?

    Asked by attempting the construction, which is how SciPy answers.
    `[M]` a miss costs ~0 s (it raises before building) and the largest
    hit — order 131, 5810 nodes — costs 1.85 ms.
    """
    try:
        Quadrature.lebedev(order)  # the registry's factory: one rule list (ruling 32)
    except NotImplementedError:
        return False
    return True


# ---------------------------------------------------------------------------
# Spec dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class QuadratureSpec:  # concept: quadrature | role: part | kind: element
    r"""Tagged registry entry for a single quadrature rule.

    Each spec captures what a caller needs to build a rule at a degree: the
    parameter-inversion logic, the two catalogue flags, and the factory.

    ⛔ It read "…needs to *rank* the rule **without instantiating it**:
    … the node-count formula …" until 2026-08-14, and carried an
    ``expected_node_count`` callable to serve that. The premise had not
    held since the 2026-08-02 stage-0/1 rework: the selector (retired on
    2026-11-20) **had to** instantiate every surviving candidate, because
    its stages 0 and 1 asked the nodes. So the cost was computed twice — once
    as a formula here, once implicitly by the measure the selector was already
    holding — and `[M]` the two agreed on all 25 shipped configurations, which
    is exactly what makes a twin source dangerous rather than safe. The rank read
    :attr:`~orpheus.numerics.measure.DiscreteMeasure.n_points` off the built
    measure, so a formula cannot drift from the rule it claims to describe.
    (The first family that would have broken it already exists:
    ``folded_product`` quotients by a mirror, so ``n_mu * n_phi``
    over-counts it by 2×.)

    It deliberately carries **no invariance group**. A rule's symmetry
    is parameter-dependent (``product_mu_phi``'s is
    :math:`D_{n_\phi h}`), so a parameter-free field cannot state it
    truthfully — and the field that used to sit here declared
    :math:`SO(2)`, which is false for every finite point set on
    :math:`S^2`. The group is now **computed** from the instantiated
    measure by
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`.
    The measure's own ``invariance_group`` metadata survives as the
    single source of that tag (``DiscreteMeasure.phase`` reads it);
    duplicating it here was a twin source of truth.

    Attributes
    ----------
    name : str
        Human-readable identifier (used in messages and Sphinx
        documentation).
    factory : callable
        The :class:`~orpheus.numerics.quadrature.directional.Quadrature`
        classmethod that makes the rule (``Quadrature.gauss_legendre``,
        ``.lebedev``, ``.level_symmetric``, ``.product``): the registry
        holds the same rule list the class does, not a second one over
        the module rule functions (ruling 32, 2026-10-19). It accepts
        the keyword parameters its signature names (the keys
        :attr:`degree_of_exactness_for` returns) and returns a
        :class:`~orpheus.numerics.quadrature.directional.Quadrature`,
        whose ``level_structure`` rides with the measure, so the
        :func:`level_symmetric_sn` and :func:`product_mu_phi` pairing
        of a measure with a
        :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
        is no longer unpacked and dropped here.
    degree_of_exactness_for : callable
        Inverse of the rule's degree-of-exactness formula:
        ``target_degree -> dict[str, Any]`` returning the smallest
        parameter dict :math:`(p_1, p_2, \ldots)` such that the rule
        evaluated at those parameters satisfies
        :math:`\deg(Q) \ge \text{target\_degree}`. Returns ``None`` if
        no parameter combination supported by this rule reaches the
        target degree (e.g., a Lebedev order higher than scipy's
        tabulated 47).
    axis_aligned : bool
        Are nodes aligned with coordinate axes? Lebedev's :math:`O_h`
        construction includes axis-aligned points; level-symmetric
        does *not* place ordinates on the axes. Cartesian-friendly
        consumers may prefer axis-aligned rules.
    half_range_clean : bool
        Does ``measure.restrict(lambda x: x > 0)`` give a valid
        quadrature without re-normalisation? Gauss-Legendre on
        :math:`[-1, 1]` is half-range-clean: the
        :math:`n`-point rule has :math:`n/2` positive nodes (when
        :math:`n` is even) carrying half the total weight. Lebedev
        and level-symmetric rules are not (they are designed to span
        the full sphere).

    Two flags are NOT attributes, because they are facts of the BUILT rule
    and a parameter-free field cannot state them truthfully:
    ``positive_weights`` (every weight :math:`> 0`; negative-weight rules
    amplify quadrature noise) and ``level_structured`` (the rule carries a
    :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`, which
    the cylindrical redistribution of Bailey-Morel-Chang 2010 Eq. (50)
    marches over). `[M]` 2026-10-19: the typed ``positive_weights=True`` on
    Lebedev was false for 3 of its 32 orders (13, 25, 27); the selector
    (retired with ERR-118 on 2026-11-20) returned order 13 for cartesian2d at
    degree 12, whose weights include -0.372.
    """

    name: str
    factory: Callable[..., Quadrature]
    degree_of_exactness_for: Callable[[int], dict[str, Any] | None]
    axis_aligned: bool
    half_range_clean: bool

    def build(self, parameters: dict[str, Any]) -> Quadrature:
        """Instantiate the rule at ``parameters``.

        The factory is a :class:`~orpheus.numerics.quadrature.directional.Quadrature`
        classmethod, so the result is the quadrature itself, level structure
        included. A consumer reads ``.measure`` where it asks the domain,
        symmetry and exactness questions.
        """
        return self.factory(**parameters)


# ---------------------------------------------------------------------------
# Per-rule degree-of-exactness inversion
# ---------------------------------------------------------------------------
#
# Each rule's degree formula is parameter-dependent. The inversion finds
# the smallest parameter dict satisfying deg(Q) >= target_degree, or
# returns None if no available parameter reaches that degree.
#
# We keep these as module-level helpers (not closures inside the spec
# constructor) to make stack traces clean and to make the implementations
# directly testable.


def _gl1d_invert(target_degree: int) -> dict[str, Any] | None:
    """Gauss-Legendre 1-D: ``deg = 2n - 1``, so ``n = ceil((d + 1) / 2)``.

    SN slab consumers traditionally require ``n`` even (one ordinate
    per direction sense). We round up to the next even value so the
    resulting rule is usable directly in slab sweeps without further
    adjustment.
    """
    if target_degree < 0:
        return None
    n_min = (target_degree + 1 + 1) // 2  # ceil((d+1)/2)
    if n_min < 1:
        n_min = 1
    if n_min % 2 == 1:
        n_min += 1
    return {"n_ordinates": n_min}


def _lebedev_invert(target_degree: int) -> dict[str, Any] | None:
    r"""Lebedev: ``deg = order``; pick the smallest order SciPy serves.

    The admissible set is gappy and is discovered by
    :func:`_lebedev_orders`, never tabulated here. If no available order
    reaches ``target_degree``, the rule cannot satisfy the constraint and
    the inversion returns ``None`` (a caller reads it as "this spec
    cannot reach the degree").

    ⛔ Read ``_LEBEDEV_ORDERS``, a frozen 18-tuple topping at 47, until
    2026-08-14. `[M]` that made this return ``None`` for every target in
    :math:`(47, 131]` — degrees the installed SciPy serves and
    :func:`~orpheus.numerics.quadrature.rules_sphere.lebedev_sphere`
    builds. See :func:`_lebedev_orders` for the measurement.
    """
    if target_degree < 0:
        return None
    for order in _lebedev_orders():
        if order >= target_degree:
            return {"order": order}
    return None



def _ls_sn_invert(target_degree: int) -> dict[str, Any] | None:
    r"""Level-symmetric :math:`S_N`: the CHEAPEST order meeting the target.

    Carlson-Lathrop's level-symmetric rule has parameter-dependent
    exactness with **no formula in** :math:`N`. `[M]` 2026-08-14 the realized
    degree by order is :math:`S_2 \to 3`, :math:`S_4 \to 5`,
    :math:`S_6 \to 7`, :math:`S_8 \to 9`, :math:`S_{10} \to 11`,
    :math:`S_{12} \to 11`, :math:`S_{14} \to 15`, :math:`S_{16} \to 15`,
    :math:`S_{18} \to 17` (the rule-side docstring on
    :func:`~orpheus.numerics.quadrature.rules_sphere.level_symmetric_sn`
    carries the same table and its provenance). So the degree is **asked of
    each candidate**, exactly as the frontier is.

    The search window is two orders wide, and both ends are theorems rather
    than guesses: :math:`\deg(N) \ge N - 1` always, so any order from
    :math:`d + 1` up that builds meets the target and there is no reason to
    look higher; and `[M]` :math:`\deg(N) \le N + 1` on every buildable
    order, so nothing below :math:`d - 1` can reach it. At most two builds.

    ⛔ This inverted :math:`\deg = N - 1` directly until 2026-08-14 —
    smallest even :math:`N \ge d + 1`, one build, no question asked. That is
    **safe**, since :math:`N - 1` is a true lower bound, and it was
    **expensive**: `[M]` it returned an order two higher than the one that
    achieves the target at **6 of the 9 buildable orders**. Target 3 was met
    by :math:`S_2` (8 nodes) and it returned :math:`S_4` (24) — 3×; target
    15 by :math:`S_{14}` (224) and it returned :math:`S_{16}` (288), which
    `[M]` achieves only degree 15 as well, so those 64 nodes bought nothing
    whatsoever. Stage 4 ranks by node count, so the over-shoot priced the
    family above its true cost and could lose it a minimum-points tie-break
    it should win.

    ⛔ The summary line read *"conservative ``deg = N - 1``"* and the body
    claimed the wrapper *"records ``degree_of_exactness = N - 1`` as a
    conservative lower bound"*. The wrapper has recorded a build-measured
    degree since #337; "conservative" named the over-shoot without
    recognising it as one.

    ⛔ **The family is BOUNDED, and the bound is asked of the family rather
    than tabulated here.** Since #327 the weights are solved per :math:`O_h`
    orbit, and the solve has no POSITIVE solution from :math:`S_{20}` on this
    node seed. `[M]` 2026-08-13: :math:`S_{14}/S_{16}/S_{18}` build with
    smallest weight ``0.012990`` / ``0.016300`` / ``1.75e-4``;
    :math:`S_{20}`/:math:`S_{22}` refuse. Above the frontier the rule does not
    exist, so this inverter returns ``None`` — the "cannot reach
    target_degree with any supported parameters" channel of this field — (the retired selector fell through to Lebedev or the product rule on it).

    ⚠ The bound is discovered by ATTEMPTING the construction, not by comparing
    against a constant. A literal here would be a second copy of a frontier
    that lives in the solve, and the two would drift the first time the node
    set changes — an inverter would then either refuse a rule that works or
    hand back one that raises. Inversion is not a hot path (`[M]` the retired
    selector had no production consumers at all), so the honest check is
    affordable.

    ⭐ **That ruling was VINDICATED by a real event, and its own prose was the
    casualty.** #337 (``59bb38a0``) moved the frontier :math:`S_{12} \to S_{18}`
    **without touching this file**, and the inverter tracked it correctly
    because it holds no literal. Meanwhile the paragraph above it read *"no
    POSITIVE solution above* :math:`S_{12}` *(*``[M]`` *min weight* ``-0.027``
    *at* :math:`S_{14}`*)"* until 2026-08-14 — a stale ``12`` sitting eight
    lines above the warning that a literal ``12`` would drift. The number was
    never wrong; it was measured on the **pre-#337 convention seed**
    :math:`\mu_1^2 = 4/(N(N+2))`, whose positivity frontier really was
    :math:`S_{12}`, and it stopped describing the shipped configuration.
    ⟹ when citing this ruling as precedent, cite the ⚠ paragraph **and** this
    one: the CODE half was right, the PROSE half became the failure it
    predicted, and the displayed rigour is exactly what makes a stale
    neighbouring measurement read as trustworthy.

    ⛔ This docstring previously read *"The construction supports any even
    N ≥ 2. There is no upper bound here, but in practice N ≥ 24 runs into
    moment-condition algebra that the simple equal-weight construction in this
    module does not capture."* Three things were wrong: there IS an upper
    bound, it is not 24, and the construction is no longer equal-weight.
    """
    if target_degree < 0:
        return None

    def _even_at_least(value: int) -> int:
        floor = max(2, value)
        return floor + 1 if floor % 2 == 1 else floor

    # The search window is two orders wide, and both ends are theorems
    # about the family rather than guesses. deg(N) >= N - 1 always (the
    # moment conditions the construction solves), so every order from
    # target + 1 upward that BUILDS meets the target — no need to look
    # past it. And `[M]` deg(N) <= N + 1 on every buildable order, so no
    # order below target - 1 can reach it.
    candidate = _even_at_least(target_degree - 1)
    ceiling = _even_at_least(target_degree + 1)

    while candidate <= ceiling:
        try:
            measure = Quadrature.level_symmetric(candidate).measure  # the registry's factory
        except ValueError:
            # Above the positivity frontier. Every larger order is too, so
            # the family cannot serve this target — exactly what the
            # ``None`` contract is for.
            return None
        degree = measure.degree_of_exactness
        if degree is not None and degree >= target_degree:
            return {"sn_order": candidate}
        candidate += 2
    return None


def _product_invert(target_degree: int) -> dict[str, Any] | None:
    """Product GL :math:`\\times` equispaced :math:`\\phi`.

    Two parameters: :math:`n_\\mu` (polar GL) gives
    :math:`2 n_\\mu - 1`; :math:`n_\\phi` (azimuthal midpoint /
    left-endpoint) gives :math:`n_\\phi - 1` for trigonometric
    polynomials. The product rule's polynomial-exactness on
    :math:`S^2` is the **minimum** of the two factors (see
    :func:`~orpheus.numerics.quadrature.rules_product.product_mu_phi`
    docstring). Match that minimum to ``target_degree`` by sizing
    each factor independently.

    The inverse picks the smallest pair satisfying the constraint:

    - :math:`n_\\mu = \\lceil (d + 1)/2 \\rceil`
    - :math:`n_\\phi = d + 1`

    For SN cylindrical-sweep compatibility the :math:`\\phi`-count
    must be at least :math:`4 n_\\mu` to avoid azimuthal-resolution
    pathology — but this is a *consumer* requirement, not a
    polynomial-exactness one. The catalogue has no path to
    request a larger ``n_phi``; a consumer that needs one builds the rule
    directly with :meth:`Quadrature.product
    <orpheus.numerics.quadrature.directional.Quadrature.product>`. The base
    inversion only meets the V criterion.
    """
    if target_degree < 0:
        return None
    n_mu_min = max(1, (target_degree + 1 + 1) // 2)
    n_phi_min = max(1, target_degree + 1)
    return {"n_mu": n_mu_min, "n_phi": n_phi_min}


# ---------------------------------------------------------------------------
# Registry population
# ---------------------------------------------------------------------------


def _registry(*specs: QuadratureSpec) -> tuple[QuadratureSpec, ...]:
    r"""Freeze the registry and refuse duplicate names, at import.

    Two guarantees the bare literal did not give:

    **Immutability.** `[M]` 2026-08-14 ``quadrature_registry`` was the only
    module-level *list*-shaped registry in ``orpheus/``; every sibling
    (``LOSS_REPRESENTATIONS``, ``SHIPPED_CLASS_A``, …) is a ``tuple``. A
    mutable global that a reader walks is an invitation to append at
    runtime and get a different answer per import order.

    **Unique names.** :attr:`QuadratureSpec.name` is not merely a label —
    it is the key a reader addresses a spec by, and a `dict` keyed on it keeps
    the LAST value for a repeated key, so two rules sharing a name would make
    one of them **silently vanish** from that view. Nothing would warn.

    ⟹ raise here, at import, rather than gate it in a test: the failure is
    a property of the registry's construction, and the earliest possible
    refusal is the one that cannot be skipped.
    """
    names = [spec.name for spec in specs]
    duplicated = sorted({n for n in names if names.count(n) > 1})
    if duplicated:
        raise ValueError(
            f"quadrature registry has duplicate spec name(s) {duplicated}; "
            f"names must be unique because a spec is addressed by "
            f"its name and a dict keyed on it would silently drop all but "
            f"the last"
        )
    return specs


quadrature_registry: tuple[QuadratureSpec, ...] = _registry(
    QuadratureSpec(
        name="GaussLegendre1D",
        # The rule the SLAB/SPHERE geometries admit is the polar marginal
        # about x — declared as such, so that the domain leg compares an orbit
        # space naming its group against a geometry naming the group it
        # spends. The chart-level `gauss_legendre_on_mu` (the product
        # rules' polar factor) is deliberately NOT registered: it names no
        # axis, and the domain leg refuses it for exactly that reason.
        factory=Quadrature.gauss_legendre,
        degree_of_exactness_for=_gl1d_invert,
        axis_aligned=True,        # 1-D nodes ARE the axis (μ-axis)
        half_range_clean=True,     # restrict(mu>0) gives n/2-point GL on [0,1]
    ),
    QuadratureSpec(
        name="LebedevSphere",
        factory=Quadrature.lebedev,
        degree_of_exactness_for=_lebedev_invert,
        axis_aligned=True,         # Lebedev grids include axis-aligned orbits
        half_range_clean=False,    # half-sphere restriction is not a valid Lebedev
    ),
    QuadratureSpec(
        name="LevelSymmetricSN",
        factory=Quadrature.level_symmetric,
        degree_of_exactness_for=_ls_sn_invert,
        axis_aligned=False,        # LS_N nodes never lie on coordinate axes
        half_range_clean=True,     # one hemisphere = exactly half the rule (Oh symmetry)
    ),
    QuadratureSpec(
        name="ProductQuadrature",
        factory=Quadrature.product,
        degree_of_exactness_for=_product_invert,
        axis_aligned=False,        # depends on phi grid; conservatively False
        half_range_clean=True,     # GL polar half restricted cleanly to mu>0
    ),
)


# ---------------------------------------------------------------------------
# Geometry → angular symmetry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AngularSymmetry:  # concept: group | role: weld | kind: element | welded: group, quadrature
    r"""What a geometry demands of its angular discretisation — a ledger.

    A geometry's symmetry does not act on the angular variable as one
    undifferentiated thing.  It splits by **what each part is used for**,
    and the three parts are three different facts that place three
    different demands on a quadrature:

    * :attr:`spent` :math:`K` — the stabiliser the dimensional reduction
      integrates away.  A slab is invariant under the full stabiliser
      :math:`O(2)_x` of its normal — every rotation about it and every
      reflection in a plane containing it — so :math:`\psi` depends on
      :math:`\Omega` only through :math:`\mu = \Omega\cdot\hat e_x`; the
      azimuth is integrated out analytically and never discretised.  What
      this part determines is the *domain*: the angular variable lives on
      the orbit space :math:`S^2/K` (:attr:`support`).  In curvilinear
      geometry its non-trivial fiber action is what appears in the sweep
      as the angular-redistribution (:math:`\alpha`) term — the reduction
      is "paid for" there.
    * :attr:`unspent` :math:`\Gamma` — the FINITE symmetry the solution
      still has in the geometry's local frame, which the reduction did not
      spend and which a FOLD may.  A z-uniform problem is even in
      :math:`\mu_z`; an axisymmetric cylinder is even under the local
      azimuthal reflection as well.  What this part determines is which folds the domain leg admits: a rule on :math:`S^2/H` is admissible iff
      :math:`H \subseteq \Gamma K` — the fold spent nothing the solution
      does not have.
    * :attr:`owed` :math:`R` — the closure a reflecting face needs,
      still owed by the quadrature: the face reflection :math:`\sigma_{\hat
      n}` maps ordinate :math:`m` to ordinate :math:`m'` *exactly* only if
      the node set is closed under :math:`\sigma_{\hat n}`.  It cannot be
      integrated away; it must be realized as a permutation of the
      ordinates (the symmetry leg).

    ⚠ Until R3 of #434 (2026-09-03) the ledger had two entries and the
    owed closure did the unspent symmetry's job: the domain leg admitted any fold
    by a subgroup of what the geometry OWED, so a :math:`\sigma_y` fold —
    legitimate on the cylinder — was admitted for ``cartesian2d``, whose
    solution is even in :math:`\mu_z` only, and 2 of the 4 sweep quadrants
    were empty.  A symmetry recorded for one job (closing a face) was
    spent on another (licensing a fold).  And the two entries were
    documented as :math:`G^0` and :math:`G/G^0`, a factorisation the
    slab's :math:`O(2)_x` — disconnected — never was.

    The first split (spent from owed, 2026-08-02) is what made the
    selection gate expressible at all: the table it replaced recorded only
    the spent group and compared it to a rule's declared group, which is
    unsatisfiable by construction — no finite point set on :math:`S^2` is
    :math:`SO(2)`-closed — so the gate could only ever pass on a false
    declaration, and did.

    Attributes
    ----------
    spent : SubgroupOfO3
        :math:`K`, the stabiliser spent by the reduction.  Determines
        :attr:`support` and :attr:`reference`.  Use
        :attr:`SubgroupOfO3.Trivial` when the reduction spends nothing (a
        2-D or 3-D geometry retains both angular degrees of freedom).
    unspent : SubgroupOfO3
        :math:`\Gamma`, the finite symmetry the solution keeps in the local
        frame — the fold licence.  Must be finite (its elements are
        enumerated by the coverage test).
    owed : SubgroupOfO3
        :math:`R`, the finite reflection closure a quadrature must realize as an
        ordinate permutation for a reflecting face.  Must be finite.
    """

    spent: SubgroupOfO3
    unspent: SubgroupOfO3
    owed: SubgroupOfO3

    def __post_init__(self) -> None:
        for role, group, job in (
            ("unspent", self.unspent, "its elements are enumerated by the fold licence"),
            ("owed", self.owed, "it must be realized as an ordinate permutation"),
        ):
            if not group.realization.is_finite:
                raise ValueError(
                    f"AngularSymmetry.{role} must be a finite group ({job}); "
                    f"{group.name} is continuous"
                )

    @property
    def support(self) -> Manifold:
        r"""The angular domain :math:`S^2/K` — the orbit space of the SPENT
        stabiliser (disconnected on the 1-D rows: :math:`O(2)_x`), derived,
        not declared.

        Returns the point set a rule's measure must live on to be admissible
        for this geometry. Deriving it (rather than storing a second,
        independent column) keeps the domain and the spent group from drifting
        apart: they are one fact.

        ⚠ Until 2026-09-01 this was annotated ``-> str`` while its comparison
        partner :attr:`DiscreteMeasure.support` was annotated ``-> Space``, and
        ``Space`` **was** ``str`` — so the two were the same type to a checker
        and different strings to every census. That is why an audit of the
        ``Space`` alias could not return this method even though
        :meth:`admits_domain` compares it against the alias on every rule.
        """
        spent = self.spent
        if spent == SubgroupOfO3.Trivial:
            # S²/{e} = S² is a theorem, and the sphere is the name every
            # 2-D / 3-D rule declares. The catalogue's derived identity
            # quotient (`_mod_trivial`) spells the same point set
            # `S^2/Trivial` as the derivation's OUTPUT and is pinned
            # against this row's realization; the geometry's domain is
            # the base itself.
            return SPHERE
        # S²/O(2)_a — the polar marginal about the spent axis. The orbits
        # of the stabiliser are the constant-μ circles, so the quotient is
        # parameterised by μ alone; the catalogue derives it, and a
        # group with no entry refuses there, naming the missing WORK.
        # Until 2026-09-01 this row was a hand-written `COSINE_INTERVAL`,
        # which is the CHART of that orbit space and not the orbit space —
        # exactly the axis-blind spelling a slab rule could share with a
        # spatial interval (tracker 2.4).
        return SPHERE.quotient(spent)

    @property
    def reference(self) -> ReferenceMeasure:
        r"""The measure a degree must be **against** — READ off the orbit space.

        Transport integrates :math:`\phi = \int \psi \, d\Omega`
        **unweighted**, so the reference is Lebesgue measure on whatever
        angular domain the dimensional reduction left behind — the
        pushforward of :math:`d\Omega` along the quotient map of the spent
        group, which is a derivation OUTPUT of the catalogue entry
        (:func:`~orpheus.numerics.generating_measure.pushforward_reference`,
        the measure layer's answer for the orbit space). Like
        :attr:`support`, it is a function of the spent group alone: the
        domain and the measure on it are one fact, both derived from
        the spent group's orbit space (``pushforward_reference`` is keyed on
        that orbit space's group).  Until #429 tracker 3.1 (2026-09-02) this
        property was the hat-box theorem's second home — ``LEGENDRE``
        tabulated here for any rotation axis — the same twin :attr:`support`
        was until tracker 2.4 (`[M]` one production reader, the retired selector's V
        stage).

        The one arm that is not a catalogue read is the geometry that spends
        NOTHING. Its domain is the bare sphere, not the trivial quotient
        (:attr:`support` returns ``SPHERE`` so that every 2-D / 3-D rule's
        declared support matches at the domain leg), and Lebesgue on :math:`S^2`
        is ``UNIFORM_ON_SPHERE`` — a property of the base, which is why
        ``pushforward_reference`` of the trivial entry is ``None`` (user-ruled
        2026-09-02).

        ⭐ Why this is not redundant with :attr:`support`. The support says
        *which space*; the reference says *which measure on that space*, and
        a rule can get the first right and the second wrong. That is not
        hypothetical — it is the defect this property exists to make
        unspellable: `[M]` Gauss-Legendre and Gauss-Chebyshev at :math:`n=4`
        agree on support (:math:`[-1,1]`), on orthogonal system
        (``ALGEBRAIC``) **and** on degree (7). They differ only in the
        reference, and integrating :math:`x^6` against the *unweighted*
        measure the transport equation actually asks for, Legendre gives
        :math:`0.285714` (exact :math:`2/7`) while Chebyshev gives
        :math:`0.981748` — off by **0.696**, at full advertised degree.

        ``LEGENDRE`` is Lebesgue measure on :math:`[-1,1]`: its weight is
        :math:`w(x) = 1`, and `[M]` its mass is exactly 2. It is named for
        the polynomial family it generates, not for a weighting; the
        hat-box constant :math:`2\pi` is no claim's business.

        Raises
        ------
        NotImplementedError
            If the spent group's orbit space has no reference — no
            shipped :class:`ReferenceMeasure` realization spells the
            pushforward (the mirrors' weighted disk measure). The message
            names the missing WORK. A spent group with no catalogue entry
            at all is refused one call earlier, by :attr:`support`.
        """
        support = self.support
        if isinstance(support, Quotient):
            reference = pushforward_reference(support)
            if reference is None:
                raise NotImplementedError(
                    f"the catalogue entry for {support.name} carries no "
                    f"exactness reference: no shipped ReferenceMeasure "
                    f"realization spells the pushforward of dOmega along its "
                    f"quotient map. Add one to orpheus/numerics/exactness.py "
                    f"and answer it in `pushforward_reference` in "
                    f"orpheus/numerics/generating_measure.py."
                )
            return reference
        # The bare sphere: a geometry that spends nothing integrates
        # against Lebesgue on S^2 itself, and the sphere is the base, not
        # a catalogue entry.
        return UNIFORM_ON_SPHERE

    def admits_domain(self, measure: DiscreteMeasure) -> bool:
        r"""The domain leg — the rule lives on this geometry's angular domain, or
        on a FOLD of it by a symmetry the solution still has.

        Two conditions, both read off the lattice and neither declared: the
        descent arrow :math:`S^2/K \to X` onto the rule's support must EXIST
        (:func:`~orpheus.geometry.manifold.quotient_onto` — the same arrow
        a frame's G0 reads, of which equality is the identity case), AND
        the group :math:`H` the rule's support was quotiented by
        (:attr:`DiscreteMeasure.quotient_group
        <orpheus.numerics.measure.DiscreteMeasure.quotient_group>`,
        :math:`\{e\}` for a bare rule) must lie in :math:`\Gamma K`
        (:meth:`SubgroupOfO3.is_subset_of_product
        <orpheus.geometry.symmetry.SubgroupOfO3.is_subset_of_product>`):
        a fold may spend only what the solution has.  Total — every
        (rule, geometry) pair answers, no arm raises.

        `[M]` the slab's own rule, :math:`O(2)_x \subseteq \{e\}\,O(2)_x`,
        admitted; the shipped :math:`\sigma_y` fold on the cylinder,
        :math:`\sigma_y \subseteq D_{1h}`, admitted; the same fold on
        ``cartesian2d``, :math:`\sigma_y \not\subseteq \sigma_z`, REFUSED.
        Until R3 of #434 (2026-09-03) the second condition read "what the
        arrow spends lies in the OWED closure" — the closure a reflecting
        face needs, read as a fold licence — and the ``cartesian2d`` fold
        was admitted; until 2026-09-02 this was ``measure.support ==
        self.support`` — equality, no lattice — and it refused the shipped
        cylinder configuration (#429 tracker 2.2b).
        """
        return self.domain_refusal(measure) is None

    def domain_refusal(self, measure: DiscreteMeasure) -> str | None:
        r"""The domain leg's verdict WITH its reason — the one failing clause, or
        ``None`` when the rule is admitted.

        :meth:`admits_domain` is this predicate's ``is None``; a caller
        appends the returned reason.  Spelled once, so the report cannot
        drift from the decision (`[M]` 2026-09-03, over 4 geometries × 7
        rules — the five ``Quadrature`` factories, ``gauss_legendre_on_mu``
        and a constructed σ_z fold — the 17 stage-0 refusals split 14 arrow / 3 coverage /
        0 both — a disjunctive message named a satisfied fact on all 17).
        """
        if quotient_onto(self.support, measure.support) is None:
            return (
                f"the rule's nodes live on {measure.support.name}, and "
                f"{self.support.name} has no descent arrow onto it"
            )
        # The rule's support is a quotient of this geometry's domain, so it
        # is angular and the group it was folded by is a subgroup of O(3);
        # ``None`` — no fold — is the trivial group HERE, honestly, because
        # the arrow check ran first (a non-angular support never reaches
        # this line, and on the SPENT side ``None`` may mean "nothing acts").
        folded_by = (
            SubgroupOfO3.Trivial if measure.quotient_group is None
            else measure.quotient_group
        )
        if not folded_by.is_subset_of_product(gamma=self.unspent, kappa=self.spent):
            return (
                f"the rule's nodes live on {measure.support.name}, a fold by "
                f"{folded_by.name}, which is not a symmetry the solution has "
                f"in this geometry's frame (unspent {self.unspent.name}, "
                f"spent {self.spent.name})"
            )
        return None

    def admits_symmetry(self, measure: DiscreteMeasure) -> bool:
        """The symmetry leg — is the rule closed under the owed discrete symmetry?

        Computed from the nodes, never read from a declared tag — the
        measure's own question (:meth:`DiscreteMeasure.is_invariant_under
        <orpheus.numerics.measure.DiscreteMeasure.is_invariant_under>`,
        R2 of #434).
        """
        return measure.is_invariant_under(self.owed)


# Static table — one entry per supported geometry. New geometries
# (hexagonal lattice, 2-D / 3-D spherical, …) are added here, not in
# a consumer.
#
# `slab` / `sphere`: the 1-D reductions. The stabiliser O(2) of the
#   local axis — the rotations about it and the mirrors through it — is
#   spent, leaving a μ-marginal; what remains owed is the
#   forward/backward reflection μ → −μ, which is what pairs the two
#   sweep senses and what a reflecting end face consumes.
#
#   That reflection is σ_x, and until 2026-08-02 these two rows said
#   `Z2` — a parameter-free tag REALIZED as σ_z. The polar marginal
#   embeds as (μ, 0, 0) (`Quadrature.axis_cosines`), so μ → −μ is the
#   mirror through the plane normal to **x**; σ_z is a different matrix.
#   The rows were harmless only because two independent guards kept the
#   σ_z realization unreachable (the 1-D dispatch never builds it, and
#   the domain leg rejected every S² rule for these geometries before the
#   symmetry leg asked, in the retired selector's pipeline) — an accident of
#   dispatch, not a designed invariant.
# `cylinder` / `cartesian2d`: nothing continuous is spent — both
#   angular degrees of freedom survive. Owed is D_2h ≅ (ℤ₂)³, the group
#   generated by the three coordinate-plane mirrors, whose chambers are
#   exactly the octants the sweep decomposes into.
#
# `cartesian2d` carried O_h before 2026-08-02, a 6× over-claim: O_h
# demands the x↔z exchange and the diagonal mirrors, which are
# symmetries of a *cube*, never of a z-uniform problem.

GEOMETRY_ANGULAR_SYMMETRY: dict[str, AngularSymmetry] = {
    # The spent group names its AXIS (2026-09-01): the polar marginal
    # embeds along x, the owed mirror's normal is x, and the real
    # spherical-harmonic pole is x — one axis, stated on both halves.
    # Until then `SO2` was a bare member realized about z, so the slab's
    # own rule read as NOT invariant under the group it was quotiented
    # by (the angular-spaces plan's Part IV obstacle 1).
    # And it is the FULL stabiliser O(2)_x (2026-09-02, #432), not its
    # rotation half: a slab is symmetric under y -> -y as well as under
    # every rotation about x, and with O(2)_x spent nothing finite is left
    # unspent (a slab solution is not even in mu); the owed closure is the
    # end face's sigma_x.  Under SO(2)_x the unspent group would be a
    # Klein four-group; the orbit space is the same either way, and its
    # catalogue entry is named by the stabiliser.
    "slab": AngularSymmetry(
        spent=SubgroupOfO3.O2("x"),
        unspent=SubgroupOfO3.Trivial,
        owed=SubgroupOfO3.Mirror("x"),
    ),
    "sphere": AngularSymmetry(
        spent=SubgroupOfO3.O2("x"),
        unspent=SubgroupOfO3.Trivial,
        owed=SubgroupOfO3.Mirror("x"),
    ),
    # z-uniform (even in the AXIAL cosine mu: sigma_z) AND azimuthally
    # symmetric (even in the azimuthal cosine xi: sigma_y) — D_1h about the
    # radial axis, {e, sigma_y, sigma_z, C_2^x}; the derivation is the
    # module docstring's ledger row.  The plan's [R] row said
    # Mirror('y') alone; the opener refuted it on the plan's own
    # z-uniformity argument for cartesian2d (2026-09-03).
    "cylinder": AngularSymmetry(
        spent=SubgroupOfO3.Trivial,
        unspent=SubgroupOfO3.Dnh(1),
        owed=SubgroupOfO3.Dnh(2),
    ),
    # z-uniform only: even in mu_z, and in nothing else.
    "cartesian2d": AngularSymmetry(
        spent=SubgroupOfO3.Trivial,
        unspent=SubgroupOfO3.Mirror("z"),
        owed=SubgroupOfO3.Dnh(2),
    ),
}


__all__ = [
    "AngularSymmetry",
    "GEOMETRY_ANGULAR_SYMMETRY",
    "QuadratureSpec",
    "quadrature_registry",
]
