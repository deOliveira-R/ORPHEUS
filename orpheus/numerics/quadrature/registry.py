r"""Tagged quadrature registry with G + V + structural selector.

This module installs the **automated quadrature selector** that closes
the loop on Issues 3 (symmetry tagging) and 4 (rule-side wrappers): given
a geometry and a target polynomial degree, walk the registry of all known
rules and return the cheapest measure whose tags satisfy the geometry's
demands.

Motivation
----------

Quadrature selection is *fundamentally* a five-stage filter, in priority
order:

0. **Domain compatibility.**
   :math:`\mathcal{D}_Q = S^2 / G^0_{\text{geom}}` — the rule's nodes
   must live on the angular domain the geometry's *dimensional
   reduction* left behind. A slab integrates :math:`\phi` out
   analytically and wants a :math:`\mu`-marginal on :math:`[-1,1]`; a
   cylinder retains both angular degrees of freedom and wants
   :math:`S^2`. This stage is **not** a refinement of the symmetry
   stage below — it is the other half of the same decomposition (see
   :class:`AngularSymmetry`), and without it the symmetry stage alone
   admits a Lebedev rule for a slab.
1. **G compatibility** (symmetry).
   :math:`\Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q)` — the
   rule must be invariant under the geometry's *discrete residual*
   symmetry, the half of its symmetry group that survives the
   reduction and must therefore be realized as an ordinate
   permutation. A quadrature with **less** symmetry imprints spurious
   low-order asymmetry on a symmetric problem (Lebedev 1976, §1;
   Stiefel & Fässler 1979 Ch. 5) and cannot represent a reflecting
   boundary exactly. A quadrature with **more** symmetry is fine —
   the extra symmetry is unused, not violated.

   :math:`\operatorname{Sym}(Q)` is **computed from the rule's nodes**
   (:meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant`),
   never read from a declared tag. A declared group is a claim with no
   construction behind it, and the claim that used to sit here was
   false: ``product_mu_phi`` advertised :math:`SO(2)`, which no finite
   point set on :math:`S^2` can satisfy, and that falsehood was the
   only reason this gate admitted the product rule for a cylinder.
2. **V compatibility** (exactness).
   :math:`\mathcal{E}(Q) \succeq (\lambda_{\text{geom}}, d)` — the rule's
   exactness **claim** must dominate the query's, where one claim dominates
   another iff they are against the **same reference measure** and the
   degree is at least as large.

   The reference is not decoration. A degree is an *index into the
   orthogonal system of a measure*, so the same integer means different
   things against different measures, and a rule can match on space, on
   system and on degree while being exact against the wrong thing: `[M]`
   Gauss-Legendre and Gauss-Chebyshev at :math:`n = 4` agree on all three
   and differ by **0.696** on :math:`\int_{-1}^{1} x^6 dx`. Transport
   integrates :math:`\phi = \int \psi \, d\Omega` unweighted, so
   :math:`\lambda_{\text{geom}}` is Lebesgue measure on the geometry's
   angular domain — derived from the spent group by
   :attr:`AngularSymmetry.reference`, exactly as the domain is.

   ⛔ This conjunct read :math:`\deg(Q) \ge d` until 2026-08-14 and was
   never *checked* at all — the stage only **inverted**, asking each spec
   for parameters and trusting the answer. Both halves are now verified
   against the claim the built rule carries: see the two-part evaluation
   note below.

   Most rules have a parameter-dependent degree (:math:`2n - 1` for
   Gauss-Legendre, ``order`` for Lebedev), and the selector inverts this to
   choose the smallest parameter set satisfying the constraint.
   ⚠ Level-symmetric :math:`S_N` has **no formula**: its
   realized degree is build-measured, and :math:`N - 1` is only a lower
   bound the inversion uses — see :func:`_ls_sn_invert` for the measured
   table and the over-shoot it causes. (This list read ":math:`N - 1` for
   level-symmetric :math:`S_N`" as though it were the degree, until
   2026-08-14.)
3. **Structural compatibility**.
   Boolean flags the consumer can request (``positive_weights``,
   ``axis_aligned``, ``level_structured``, ``half_range_clean``) — the
   rule's flag set must be a superset of the requested set.
4. **Minimum points**.
   Among rules passing 0-3, pick the smallest :math:`N`. This is the
   cost minimisation: each ordinate is a sweep direction in the SN
   solver, and sweeps dominate the runtime. The count is read off the
   **built** measure, which stages 0 and 1 have already forced into
   existence — never from a parallel formula that could drift from it.
   (Read "passing 1-3" here until 2026-08-14, from before stage 0.)

Formally,

.. math::
   :label: quadrature-selection-criterion

   Q^{\star} \;=\; \arg\min\Bigl\{\, n(Q) \;:\;\;
   \mathcal{D}_Q = S^2 / G^0_{\text{geom}}
   \;\wedge\; \Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q)
   \;\wedge\; \mathcal{E}(Q) \succeq (\lambda_{\text{geom}},\, d)
   \;\wedge\; F_{\text{req}} \subseteq F_Q
   \,\Bigr\},

where :math:`n(Q)` is the number of nodes,
:math:`\mathcal{D}_Q` is the domain the rule's nodes live on,
:math:`G^0_{\text{geom}}` and :math:`\Gamma_{\text{geom}}` are the
continuous and discrete halves of the geometry's angular symmetry
(:class:`AngularSymmetry`), :math:`\operatorname{Sym}(Q)` is the group
the rule's nodes are actually invariant under,
:math:`\mathcal{E}(Q)` is the rule's exactness claim
(:class:`~orpheus.numerics.exactness.ExactnessClaim` — a reference
measure together with a degree against it),
:math:`\lambda_{\text{geom}}` is the measure the geometry integrates
against, and
:math:`F_Q \subseteq \{\text{positive\_weights}, \text{axis\_aligned},
\text{level\_structured}, \text{half\_range\_clean}\}` is the rule's
structural-flag set. **Domination** is

.. math::

   \mathcal{E}(Q) \succeq (\lambda, d)
   \quad\Longleftrightarrow\quad
   \operatorname{ref}\mathcal{E}(Q) = \lambda
   \;\wedge\;
   \deg\mathcal{E}(Q) \ge d ,

a partial order, not a total one: two claims against different
references are simply **incomparable**, which is the whole point —
:math:`\deg = 7` against ``chebyshev_t`` is neither better nor worse
than :math:`\deg = 7` against ``legendre``, it is an answer to a
different question.

The conjunction is order-free, but the *evaluation* is not, and it
splits **one conjunct across two moments**:

* :math:`\mathcal{E}(Q)` is a property of an instantiated rule, and the
  nodes need parameters that only the V constraint determines. So the V
  stage runs **first as an inversion** (ask the spec for the smallest
  parameters reaching :math:`d`), the rule is built, and then the claim
  it actually carries is **verified** against
  :math:`(\lambda_{\text{geom}}, d)`.
* The verification runs *after* stages 0 and 1 rather than immediately,
  because :math:`\lambda_{\text{geom}}` encodes :math:`\mathcal{D}_Q`: a
  rule on the wrong space fails both, and the domain stage gives the
  better diagnosis. What the verification uniquely catches is a rule on
  the **right** space carrying the **wrong measure**.

This split is a real consequence of the theorem, not an implementation
detail. Likewise a rule's invariance group is **parameter-dependent**
(the product rule's is :math:`D_{n_\phi h}`), so no parameter-free field
on :class:`QuadratureSpec` can express it, and the older design that
tried to had to lie. The same is true of the claim: it is read off the
built measure, never declared on the spec.

The :func:`select_quadrature` function returns both the chosen measure
**and** a :class:`SelectionLog` listing every rejected candidate with
the stage-and-reason that disqualified it. The log makes the choice
explainable for the SN solver's diagnostic output and for teaching:
"why was rule X picked over rule Y?" answers itself.

Why a registry, not a hardcoded ladder
--------------------------------------

``solve_sn`` keeps explicit quadrature-passing as the canonical API
(per the design note in ``.claude/plans/sn_reshape.md`` Issue 5):
explicit is better than implicit, and a quadrature is a load-bearing
modelling choice the user must own. The registry adds
:func:`select_quadrature` as **opt-in convenience**: a default for
"run this geometry with target degree :math:`d`" prototyping.

⛔ This paragraph also claimed the registry was "a documentation artifact —
every entry's docstring becomes a Sphinx table row narrating the trade-off"
until 2026-08-14. That described a mechanism that **cannot exist**:
:class:`QuadratureSpec` is a dataclass, so its four *instances* share the class
docstring (`[M]` ``spec.__doc__ is type(spec).__doc__`` for all four) and an
instance cannot carry prose of its own — the entries at
:data:`quadrature_registry` carry inline ``#`` comments, which no builder reads.
And nothing here renders: this module has **no** ``automodule`` directive
(deliberately — see ``docs/api/discrete_ordinates.rst``, which cross-references
it rather than rendering it to avoid duplicate-label collisions with the theory
pages), so no docstring in this file becomes a table row at any severity.
⟹ the per-rule narration lives in the theory page
:doc:`/theory/foundations/discrete_measures`, and a rule's own docstring is on
its factory in ``rules_*.py``. If per-entry narration is ever wanted *here*, it
needs a field to live in — it cannot be inferred from a shared class docstring.

Geometry → angular-symmetry assignment
--------------------------------------

Each geometry declares the two halves of its angular symmetry (see
:class:`AngularSymmetry`): the continuous part :math:`G^0` that the
dimensional reduction **spends**, and the discrete residual
:math:`\Gamma` still **owed** to the quadrature.

.. list-table::
   :header-rows: 1
   :widths: 18 12 12 12 46

   * - Geometry
     - :math:`G^0` (spent)
     - Domain :math:`S^2/G^0`
     - :math:`\Gamma` (owed)
     - Rationale
   * - ``"slab"``
     - :math:`SO(2)`
     - :math:`[-1,1]`
     - :math:`Z_2`
     - 1-D in :math:`z`. Azimuthal rotation about the slab normal is
       integrated out analytically, so the angular variable is
       :math:`\mu = \cos\theta` alone. What remains owed is
       :math:`\mu \to -\mu`, the reflection that pairs the two sweep
       senses and that a reflecting end face consumes. Gauss-Legendre
       nodes are symmetric (Stoer-Bulirsch §3.6), so it holds.
   * - ``"sphere"``
     - :math:`SO(2)`
     - :math:`[-1,1]`
     - :math:`Z_2`
     - The 1-D **radial** spherical SN reduces to GL on
       :math:`\mu_r = \cos\theta_r`, the cosine of the angle between
       the ordinate and the radial direction (Lewis & Miller 1993
       §4.4). The continuous problem has :math:`O(3)` symmetry; the
       radial reduction spends the azimuth about :math:`\hat r` and
       leaves the same :math:`Z_2` as the slab. Here the spent half is
       not free — its fiber action reappears in the sweep as the
       angular-redistribution :math:`\alpha` term.
   * - ``"cylinder"``
     - trivial
     - :math:`S^2`
     - :math:`D_{2h}`
     - An axisymmetric cylinder is :math:`\phi`-independent in
       *space*, but that does not reduce the *angular* domain: both
       angular degrees of freedom survive, so the rule must live on
       :math:`S^2`. Owed is :math:`D_{2h}`, the coordinate-plane
       mirrors. The cylindrical SN sweep additionally requires
       per-:math:`\mu` polar-level structure for the azimuthal
       redistribution coefficients; request it via the
       ``level_structured=True`` structural flag.
   * - ``"cartesian2d"``
     - trivial
     - :math:`S^2`
     - :math:`D_{2h}`
     - 2-D Cartesian (x-y). :math:`D_{2h} \cong (\mathbb{Z}_2)^3` is
       generated by the three coordinate-plane mirrors, and its
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
   both satisfiable and discriminating.

When 2-D / 3-D spherical SN lands, this table gains ``"sphere2d"`` /
``"sphere3d"`` entries; both spend nothing continuously and owe at
least :math:`D_{2h}`.

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

See Also
--------
:ref:`quadrature-selection-algorithm` (theory page) — narrates the
five-stage precedence chain in prose, with the worked examples and the
independence witnesses. It is the RENDERED statement of this design; this
module is deliberately not ``automodule``-rendered, so the page is what a
reader sees. (Read "the four-stage precedence chain" here until
2026-08-14 — the count predates stage 0, and the anchor was the whole
page rather than the section.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from collections.abc import Sequence
from typing import Any, Callable

from ..exactness import UNIFORM_ON_SPHERE, ReferenceMeasure
from ..generating_measure import LEGENDRE
from ..measure import DiscreteMeasure
from orpheus.numerics.manifold import COSINE_INTERVAL, SPHERE, Manifold
from ..symmetry import SubgroupOfO3
from .rules_1d import gauss_legendre_on_mu
from .rules_product import product_mu_phi
from .rules_sphere import lebedev_sphere, level_symmetric_sn


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
    serves **32** orders topping at **131**, so the selector refused degrees
    the tree could deliver: ``_lebedev_invert(53)`` returned ``None`` while
    ``lebedev_sphere(53)`` builds 974 nodes summing to :math:`4\pi`. A
    frozen probe is a measurement whose configuration (the installed SciPy)
    can change underneath it without anything going red.

    Cost is paid at most once per process, and only if selection is used
    (`[M]` the selector has zero production consumers). `[M]` 16.6 ms.

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
        lebedev_sphere(order)
    except NotImplementedError:
        return False
    return True


# ---------------------------------------------------------------------------
# Spec dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class QuadratureSpec:
    r"""Tagged registry entry for a single quadrature rule.

    Each spec captures what :func:`select_quadrature` needs to reach a
    candidate rule: the parameter-inversion logic, the structural flags,
    and the factory.

    ⛔ It read "…needs to *rank* the rule **without instantiating it**:
    … the node-count formula …" until 2026-08-14, and carried an
    ``expected_node_count`` callable to serve that. The premise had not
    held since the 2026-08-02 stage-0/1 rework: :func:`select_quadrature`
    **must** instantiate every surviving candidate, because stages 0 and 1
    ask the nodes. So the cost was computed twice — once as a formula here,
    once implicitly by the measure the selector was already holding — and
    `[M]` the two agreed on all 25 shipped configurations, which is exactly
    what makes a twin source dangerous rather than safe. The rank now reads
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
    :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant`.
    The measure's own ``invariance_group`` metadata survives as the
    single source of that tag (``DiscreteMeasure.phase`` reads it);
    duplicating it here was a twin source of truth.

    Attributes
    ----------
    name : str
        Human-readable identifier (used in :class:`SelectionLog`
        messages and Sphinx documentation).
    factory : callable
        Constructor returning the measure. Must accept the parameter
        keys listed in :attr:`parameters` and return either a
        :class:`DiscreteMeasure` directly or a tuple
        ``(DiscreteMeasure, ...)`` — the selector unpacks the first
        element when a tuple is returned, matching the
        :func:`level_symmetric_sn` and :func:`product_mu_phi`
        signatures that pair the measure with a
        :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`.
    parameters : dict[str, type]
        Parameter names and their types. The selector calls
        :meth:`degree_of_exactness_for` to invert the degree formula
        and obtain the minimum parameters that meet the target degree.
    degree_of_exactness_for : callable
        Inverse of the rule's degree-of-exactness formula:
        ``target_degree -> dict[str, Any]`` returning the smallest
        parameter dict :math:`(p_1, p_2, \ldots)` such that the rule
        evaluated at those parameters satisfies
        :math:`\deg(Q) \ge \text{target\_degree}`. Returns ``None`` if
        no parameter combination supported by this rule reaches the
        target degree (e.g., a Lebedev order higher than scipy's
        tabulated 47).
    positive_weights : bool
        Are all weights non-negative? Standard for Gauss-Legendre,
        Lebedev, level-symmetric :math:`S_N`, product. Negative-weight
        rules (Newton-Cotes :math:`n \ge 9`) are rejected by physics
        consumers because they amplify quadrature noise.
    axis_aligned : bool
        Are nodes aligned with coordinate axes? Lebedev's :math:`O_h`
        construction includes axis-aligned points; level-symmetric
        does *not* place ordinates on the axes. Cartesian-friendly
        consumers may prefer axis-aligned rules.
    level_structured : bool
        Does the rule expose per-:math:`\mu` polar levels?
        Cylindrical SN sweeps need this for the azimuthal
        redistribution coefficients (Bailey-Morel-Chang 2010 Eq. (50)),
        which march per level in the order Eq. (52) fixes.
        :func:`level_symmetric_sn` and :func:`product_mu_phi` both
        return a :class:`LevelStructure`; :func:`lebedev_sphere`
        does not (Lebedev points sit on :math:`O_h` orbits, not on
        polar levels).
    half_range_clean : bool
        Does ``measure.restrict(lambda x: x > 0)`` give a valid
        quadrature without re-normalisation? Gauss-Legendre on
        :math:`[-1, 1]` is half-range-clean: the
        :math:`n`-point rule has :math:`n/2` positive nodes (when
        :math:`n` is even) carrying half the total weight. Lebedev
        and level-symmetric rules are not (they are designed to span
        the full sphere).
    """

    name: str
    factory: Callable[..., Any]
    parameters: dict[str, type]
    degree_of_exactness_for: Callable[[int], dict[str, Any] | None]
    positive_weights: bool
    axis_aligned: bool
    level_structured: bool
    half_range_clean: bool

    def structural_flags(self) -> dict[str, bool]:
        """Return the spec's structural flags as a dict.

        The selector uses this to verify
        :math:`F_{\\text{req}} \\subseteq F_Q`: every requested flag
        must hold ``True`` in the spec.
        """
        return {
            "positive_weights": self.positive_weights,
            "axis_aligned": self.axis_aligned,
            "level_structured": self.level_structured,
            "half_range_clean": self.half_range_clean,
        }

    def build(self, parameters: dict[str, Any]) -> DiscreteMeasure:
        """Instantiate the rule at ``parameters`` and return its measure.

        The registered factories return either a
        :class:`~orpheus.numerics.measure.DiscreteMeasure` or a
        ``(measure, LevelStructure)`` pair. This is the **one** place
        that unpacking happens — the selector needs a measure to ask
        the domain and symmetry questions, and needs one again to
        return the winner.
        """
        result = self.factory(**parameters)
        return result[0] if isinstance(result, tuple) else result


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
    return {"n": n_min}


def _lebedev_invert(target_degree: int) -> dict[str, Any] | None:
    r"""Lebedev: ``deg = order``; pick the smallest order SciPy serves.

    The admissible set is gappy and is discovered by
    :func:`_lebedev_orders`, never tabulated here. If no available order
    reaches ``target_degree``, the rule cannot satisfy the constraint and
    the inversion returns ``None`` (the selector then rejects this
    spec at the V-compatibility stage with a clear message).

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
    exist, so this inverter returns ``None`` — the selector's own "cannot reach
    target_degree with any supported parameters" channel — and the caller falls
    through to Lebedev or the product rule.

    ⚠ The bound is discovered by ATTEMPTING the construction, not by comparing
    against a constant. A literal here would be a second copy of a frontier
    that lives in the solve, and the two would drift the first time the node
    set changes — the selector would then either refuse a rule that works or
    hand back one that raises. Selection is not a hot path (`[M]` the selector
    has no production consumers at all today), so the honest check is
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
            measure, _ = level_symmetric_sn(candidate)
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

    The selector picks the smallest pair satisfying the constraint:

    - :math:`n_\\mu = \\lceil (d + 1)/2 \\rceil`
    - :math:`n_\\phi = d + 1`

    For SN cylindrical-sweep compatibility the :math:`\\phi`-count
    must be at least :math:`4 n_\\mu` to avoid azimuthal-resolution
    pathology — but this is a *consumer* requirement, not a
    polynomial-exactness one, so we leave it to the consumer to
    request a larger ``n_phi`` explicitly via the ``parameters``
    override path. The base inversion only meets the V criterion.
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
    mutable global that a selector reads is an invitation to append at
    runtime and get a different answer per import order.

    **Unique names.** :attr:`QuadratureSpec.name` is not merely a label —
    it is what :class:`SelectionLog` reports rejections under, and the
    theory page teaches ``dict(log.rejected)["ProductQuadrature"]`` as the
    way to read one. `dict` keeps the LAST value for a repeated key, so two
    rules sharing a name would make one rejection **silently vanish** from
    that view while the log itself still listed both. Nothing would warn.

    ⟹ raise here, at import, rather than gate it in a test: the failure is
    a property of the registry's construction, and the earliest possible
    refusal is the one that cannot be skipped.
    """
    names = [spec.name for spec in specs]
    duplicated = sorted({n for n in names if names.count(n) > 1})
    if duplicated:
        raise ValueError(
            f"quadrature registry has duplicate spec name(s) {duplicated}; "
            f"names must be unique because SelectionLog reports rejections "
            f"under them and dict(log.rejected) would silently drop all but "
            f"the last"
        )
    return specs


quadrature_registry: tuple[QuadratureSpec, ...] = _registry(
    QuadratureSpec(
        name="GaussLegendre1D",
        factory=gauss_legendre_on_mu,
        parameters={"n": int},
        degree_of_exactness_for=_gl1d_invert,
        positive_weights=True,
        axis_aligned=True,        # 1-D nodes ARE the axis (μ-axis)
        level_structured=False,    # no per-mu sub-levels (it IS the mu axis)
        half_range_clean=True,     # restrict(mu>0) gives n/2-point GL on [0,1]
    ),
    QuadratureSpec(
        name="LebedevSphere",
        factory=lebedev_sphere,
        parameters={"order": int},
        degree_of_exactness_for=_lebedev_invert,
        positive_weights=True,
        axis_aligned=True,         # Lebedev grids include axis-aligned orbits
        level_structured=False,    # no polar-level structure (Oh orbits)
        half_range_clean=False,    # half-sphere restriction is not a valid Lebedev
    ),
    QuadratureSpec(
        name="LevelSymmetricSN",
        factory=level_symmetric_sn,
        parameters={"sn_order": int},
        degree_of_exactness_for=_ls_sn_invert,
        positive_weights=True,
        axis_aligned=False,        # LS_N nodes never lie on coordinate axes
        level_structured=True,     # canonical N/2 polar levels per hemisphere
        half_range_clean=True,     # one hemisphere = exactly half the rule (Oh symmetry)
    ),
    QuadratureSpec(
        name="ProductQuadrature",
        factory=product_mu_phi,
        parameters={"n_mu": int, "n_phi": int},
        degree_of_exactness_for=_product_invert,
        positive_weights=True,
        axis_aligned=False,        # depends on phi grid; conservatively False
        level_structured=True,     # n_mu polar levels by construction
        half_range_clean=True,     # GL polar half restricted cleanly to mu>0
    ),
)


# ---------------------------------------------------------------------------
# Geometry → angular symmetry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AngularSymmetry:
    r"""What a geometry demands of its angular discretisation.

    A geometry's symmetry group :math:`G` does not act on the angular
    variable as one undifferentiated thing. It splits by **how the
    action is used**, and the two halves place two different demands
    on a quadrature:

    * :attr:`continuous_isotropy` :math:`G^0` — the continuous
      subgroup that the dimensional reduction **spends**. A slab is
      invariant under every rotation about :math:`z`, so
      :math:`\psi` depends on :math:`\Omega` only through
      :math:`\mu = \Omega\cdot\hat z`; the azimuth is integrated out
      analytically and never discretised. What this half determines is
      the *domain*: the angular variable lives on the quotient
      :math:`S^2/G^0` (:attr:`support`). In curvilinear geometry its
      non-trivial fiber action is what appears in the sweep as the
      angular-redistribution (:math:`\alpha`) term — the reduction is
      "paid for" there.
    * :attr:`discrete_residual` :math:`\Gamma = G/G^0` — the finite
      residual, still **owed**. It cannot be integrated away; it must
      be realized as a permutation of the ordinates. This is the half
      a reflecting boundary condition consumes: the face reflection
      :math:`\sigma_{\hat n}` maps ordinate :math:`m` to ordinate
      :math:`m'` *exactly* only if the node set is closed under
      :math:`\sigma_{\hat n}`.

    Splitting them is what makes the selection gate expressible. The
    table this replaced recorded only :math:`G^0` and compared it to a
    rule's declared group, which is unsatisfiable by construction — no
    finite point set on :math:`S^2` is :math:`SO(2)`-closed — so the
    gate could only ever pass on a false declaration, and did.

    Attributes
    ----------
    continuous_isotropy : SubgroupOfO3
        :math:`G^0`, the half spent by the reduction. Determines
        :attr:`support`. Use
        :attr:`SubgroupOfO3.Trivial` when the reduction spends nothing
        (a 2-D or 3-D geometry retains both angular degrees of
        freedom).
    discrete_residual : SubgroupOfO3
        :math:`\Gamma`, the finite half a quadrature must realize as
        an ordinate permutation.
    """

    continuous_isotropy: SubgroupOfO3
    discrete_residual: SubgroupOfO3

    @property
    def support(self) -> Manifold:
        r"""The angular domain :math:`S^2/G^0` — derived, not declared.

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
        spent = self.continuous_isotropy
        if spent == SubgroupOfO3.SO2:
            # S²/SO(2) — the polar marginal. The orbits of the axial
            # rotation are the constant-μ circles, so the quotient is
            # parameterised by μ alone.
            return COSINE_INTERVAL
        if spent == SubgroupOfO3.Trivial:
            return SPHERE
        raise NotImplementedError(
            f"no angular domain is defined for the quotient S^2/{spent.name}; "
            f"extend AngularSymmetry.support when a geometry first spends it"
        )

    @property
    def reference(self) -> ReferenceMeasure:
        r"""The measure a degree must be **against** — derived, not declared.

        Transport integrates :math:`\phi = \int \psi \, d\Omega`
        **unweighted**, so the reference is Lebesgue measure on whatever
        angular domain the dimensional reduction left behind. Like
        :attr:`support`, it is a function of the spent group alone: the
        domain and the measure on it are one fact, and storing a second
        column would let them drift.

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
        the polynomial family it generates, not for a weighting.
        """
        spent = self.continuous_isotropy
        if spent == SubgroupOfO3.SO2:
            return LEGENDRE
        if spent == SubgroupOfO3.Trivial:
            return UNIFORM_ON_SPHERE
        raise NotImplementedError(
            f"no exactness reference is defined for the quotient "
            f"S^2/{spent.name}; extend AngularSymmetry.reference when a "
            f"geometry first spends it"
        )

    def admits_domain(self, measure: DiscreteMeasure) -> bool:
        """Stage 0 — does the rule live on this geometry's angular domain?"""
        return measure.support == self.support

    def admits_symmetry(self, measure: DiscreteMeasure) -> bool:
        """Stage 1 — is the rule closed under the owed discrete symmetry?

        Computed from the nodes, never read from a declared tag.
        """
        return self.discrete_residual.is_invariant(measure)


# Static table — one entry per supported geometry. New geometries
# (hexagonal lattice, 2-D / 3-D spherical, …) are added here, not in
# the selector itself.
#
# `slab` / `sphere`: the 1-D reductions. Azimuthal symmetry about the
#   local axis is spent, leaving a μ-marginal; what remains owed is the
#   forward/backward reflection μ → −μ, which is what pairs the two
#   sweep senses and what a reflecting end face consumes.
#
#   That reflection is σ_x, and until 2026-08-02 these two rows said
#   `Z2` — a parameter-free tag REALIZED as σ_z. The polar marginal
#   embeds as (μ, 0, 0) (`Quadrature.axis_cosines`), so μ → −μ is the
#   mirror through the plane normal to **x**; σ_z is a different matrix.
#   The rows were harmless only because two independent guards kept the
#   σ_z realization unreachable (the 1-D dispatch never builds it, and
#   stage 0 rejects every S² rule for these geometries before stage 1
#   asks) — an accident of dispatch, not a designed invariant.
# `cylinder` / `cartesian2d`: nothing continuous is spent — both
#   angular degrees of freedom survive. Owed is D_2h ≅ (ℤ₂)³, the group
#   generated by the three coordinate-plane mirrors, whose chambers are
#   exactly the octants the sweep decomposes into.
#
# `cartesian2d` carried O_h before 2026-08-02, a 6× over-claim: O_h
# demands the x↔z exchange and the diagonal mirrors, which are
# symmetries of a *cube*, never of a z-uniform problem.

GEOMETRY_ANGULAR_SYMMETRY: dict[str, AngularSymmetry] = {
    "slab": AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.SO2,
        discrete_residual=SubgroupOfO3.Mirror("x"),
    ),
    "sphere": AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.SO2,
        discrete_residual=SubgroupOfO3.Mirror("x"),
    ),
    "cylinder": AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.Trivial,
        discrete_residual=SubgroupOfO3.Dnh(2),
    ),
    "cartesian2d": AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.Trivial,
        discrete_residual=SubgroupOfO3.Dnh(2),
    ),
}


# ---------------------------------------------------------------------------
# Selection log
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SelectionLog:
    """Explainability record of a single :func:`select_quadrature` call.

    Returned alongside the chosen measure so consumers can render the
    decision in diagnostic output ("why did the selector pick rule X?").
    The log is also the right place to inspect when a selection
    **fails**: every spec that was tried is listed with the precise
    stage-and-reason that disqualified it.

    Attributes
    ----------
    geometry : str
        The geometry tag from the original call.
    angular_symmetry : AngularSymmetry
        The geometry's angular symmetry — both halves — resolved from
        :data:`GEOMETRY_ANGULAR_SYMMETRY`. Carries the required
        ``support`` and the owed ``discrete_residual``, so a reader of
        the log can reconstruct stages 0 and 1 without re-deriving
        them.
    target_degree : int
        The polynomial-exactness target from the original call.
    requested_flags : dict[str, bool]
        The structural flags requested via ``**structural``. Only the
        keys with value ``True`` are stored — these are the *required*
        flags. (False values are interpreted as "don't care".)
    chosen_spec : QuadratureSpec | None
        The spec selected, or ``None`` if no spec passed every stage.
    chosen_parameters : dict[str, Any] | None
        The parameter dict used to instantiate the chosen rule.
    rejected : list[tuple[str, str]]
        For each rejected spec, ``(spec.name, reason)``. The reason
        string identifies the failing stage (domain / symmetry / V /
        structural) and the missing predicate.
    """

    geometry: str
    angular_symmetry: AngularSymmetry
    target_degree: int
    requested_flags: dict[str, bool]
    chosen_spec: QuadratureSpec | None
    chosen_parameters: dict[str, Any] | None
    rejected: list[tuple[str, str]] = field(default_factory=list)

    def summary(self) -> str:
        """One-line summary for log output / Sphinx tables."""
        if self.chosen_spec is None:
            return (
                f"select_quadrature(geometry={self.geometry!r}, "
                f"target_degree={self.target_degree}, "
                f"flags={self.requested_flags}) -> NO MATCH "
                f"(rejected {len(self.rejected)} candidates)"
            )
        return (
            f"select_quadrature(geometry={self.geometry!r}, "
            f"target_degree={self.target_degree}) -> "
            f"{self.chosen_spec.name}"
            f"({self.chosen_parameters}) "
            f"[{len(self.rejected)} rejected]"
        )


# ---------------------------------------------------------------------------
# Selector
# ---------------------------------------------------------------------------


class QuadratureSelectionError(ValueError):
    """Raised when no registered spec satisfies the constraints.

    The exception's :attr:`log` attribute is a fully-populated
    :class:`SelectionLog` so callers can introspect *why* selection
    failed without re-running the search.
    """

    def __init__(self, message: str, log: SelectionLog) -> None:
        super().__init__(message)
        self.log = log


def select_quadrature(
    geometry: str,
    target_degree: int,
    *,
    registry: Sequence[QuadratureSpec] | None = None,
    **structural: bool,
) -> tuple[DiscreteMeasure, SelectionLog]:
    r"""Pick the cheapest registered quadrature satisfying every constraint.

    Implements :eq:`quadrature-selection-criterion` over the global
    :data:`quadrature_registry` (or a caller-supplied ``registry`` for
    testing). The **five**-stage filter — domain, symmetry, V, structural,
    minimum-points — runs in priority order and produces an explainability
    log alongside the chosen measure.

    ⛔ This read *"The four-stage filter — G, V, structural,
    minimum-points"* until 2026-08-14. Both halves were wrong, and in the
    same way: **stage 0 (domain) is missing**. It was added 2026-08-02 as
    the other half of the geometry's angular-symmetry decomposition —
    :math:`G^0` spent by the dimensional reduction, :math:`\Gamma` still
    owed to the quadrature — and without it the symmetry stage alone admits
    a Lebedev rule for a slab. See the module docstring's stages 0-4, and
    :ref:`quadrature-selection-algorithm` for the measured witness.

    Parameters
    ----------
    geometry : str
        Geometry tag — one of the keys of
        :data:`GEOMETRY_ANGULAR_SYMMETRY`
        (``"slab"``, ``"sphere"``, ``"cylinder"``, ``"cartesian2d"``).
        Case-sensitive. Unknown tags raise :class:`KeyError` immediately
        (this is a programming error, not a selection failure).
    target_degree : int
        Minimum polynomial-exactness degree :math:`d` required.
        :math:`d = 0` accepts any rule that integrates constants
        exactly.
    registry : Sequence[QuadratureSpec], optional
        Override the global registry — primarily for unit tests.
        Defaults to :data:`quadrature_registry`.
    **structural : bool
        Required structural flags. Only flags passed with value
        ``True`` are required; ``False`` and missing keys are
        interpreted as "don't care." Recognised keys are the four flag
        names of :class:`QuadratureSpec`:
        ``positive_weights``, ``axis_aligned``, ``level_structured``,
        ``half_range_clean``.

    Returns
    -------
    DiscreteMeasure
        The chosen quadrature measure (or, for rules whose factory
        returns ``(measure, structure)``, just the measure — the
        structure is dropped because the registry is for "give me a
        measure" use cases; consumers needing per-level structure
        should call the factory directly).
    SelectionLog
        Decision record: chosen spec + parameters + rejection list.

    Raises
    ------
    KeyError
        If ``geometry`` is not in :data:`GEOMETRY_ANGULAR_SYMMETRY`.
    QuadratureSelectionError
        If no registered spec satisfies every stage. The exception's
        ``log`` attribute carries the full :class:`SelectionLog` for
        introspection.
    """
    if geometry not in GEOMETRY_ANGULAR_SYMMETRY:
        raise KeyError(
            f"unknown geometry {geometry!r}; expected one of "
            f"{sorted(GEOMETRY_ANGULAR_SYMMETRY.keys())}"
        )
    if registry is None:
        registry = quadrature_registry

    angular_symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]

    # Normalise the structural request: only True values are
    # constraints. (Passing flag=False is "don't care" — we never
    # filter rules out for *not* having a flag.)
    requested_flags = {k: True for k, v in structural.items() if v}

    # Reject unknown flag keys early — better to fail loudly than to
    # silently ignore a typo.
    valid_flags = {
        "positive_weights",
        "axis_aligned",
        "level_structured",
        "half_range_clean",
    }
    unknown = set(structural.keys()) - valid_flags
    if unknown:
        raise KeyError(
            f"unknown structural flag(s) {sorted(unknown)}; "
            f"valid flags are {sorted(valid_flags)}"
        )

    rejected: list[tuple[str, str]] = []
    # The measure rides along: stages 0/1 had to build it, so rebuilding
    # the winner afterwards would be a second, divergence-capable
    # construction of the same object.
    candidates: list[
        tuple[QuadratureSpec, dict[str, Any], int, DiscreteMeasure]
    ] = []

    for spec in registry:
        # ---- Stage 2: V compatibility (FIRST — it fixes the
        #      parameters, and stages 0 and 1 need the nodes) ----------
        params = spec.degree_of_exactness_for(target_degree)
        if params is None:
            rejected.append((
                spec.name,
                f"V mismatch: rule cannot reach target_degree="
                f"{target_degree} with any supported parameters",
            ))
            continue

        measure = spec.build(params)

        # ---- Stage 0: angular domain ----------------------------------
        if not angular_symmetry.admits_domain(measure):
            rejected.append((
                spec.name,
                f"domain mismatch: geometry {geometry!r} discretises "
                f"{angular_symmetry.support.name} (= S^2/"
                f"{angular_symmetry.continuous_isotropy.name}), but the rule's "
                f"nodes live on {measure.support.name}",
            ))
            continue

        # ---- Stage 1: the owed discrete symmetry ----------------------
        #
        # Computed from the instantiated nodes. Note this is NOT
        # `residual.is_subgroup_of(measure.invariance_group)`: a declared
        # tag is only required to be TRUE of the nodes, never maximal, so
        # the lattice route can reject a rule that satisfies the owed
        # symmetry perfectly well. Asking the nodes directly cannot go
        # wrong that way.
        #
        # ⛔ This cited "the 1-D rule declares SO(2), the group its domain
        # was quotiented BY" as the worked example until 2026-08-14. `[M]`
        # gauss_legendre_on_mu(8).invariance_group is Mirror('x') — the
        # 2026-08-02 slab/sphere correction reached the rule's declared tag
        # too — and `[M]` walking maximal_invariance_groups over all four
        # registry rules, 0 of 4 are true-but-not-maximal today. The
        # argument stands on the CONTRACT (a declaration may be modest);
        # the example it used to lean on no longer exists.
        if not angular_symmetry.admits_symmetry(measure):
            rejected.append((
                spec.name,
                f"symmetry mismatch: geometry {geometry!r} owes "
                f"{angular_symmetry.discrete_residual.name}, which the rule's "
                f"nodes at {params} are not invariant under",
            ))
            continue

        # ---- Stage 2 (second half): the V conjunct, checked WHOLE -----
        #
        # Stage 2 above only INVERTED — it asked the spec for parameters
        # and trusted the answer. A degree is not a claim on its own: it
        # is an index into the orthogonal system of a REFERENCE MEASURE,
        # and the same integer means different things against different
        # measures (see orpheus.numerics.exactness). So the claim the
        # built rule actually carries is verified here, against the claim
        # the query asked for.
        #
        # It runs AFTER stages 0 and 1, not before, because the reference
        # encodes the domain: a rule on the wrong space fails both, and
        # the domain stage gives the better message. What can only be
        # caught here is a rule on the RIGHT space carrying the WRONG
        # measure — Gauss-Chebyshev for a slab, which matches on support,
        # on orthogonal system and on degree.
        claim = measure.exactness
        wanted = angular_symmetry.reference
        if claim is None:
            rejected.append((
                spec.name,
                f"V mismatch: the rule carries no exactness claim, so "
                f"nothing certifies it integrates anything exactly against "
                f"{wanted.name}",
            ))
            continue
        if claim.reference != wanted:
            rejected.append((
                spec.name,
                f"V mismatch: the rule is exact against "
                f"{claim.reference.name}, but geometry {geometry!r} "
                f"integrates against {wanted.name}. A degree is an index "
                f"into its reference's orthogonal system, so "
                f"{claim.degree} against {claim.reference.name} is not "
                f"{claim.degree} against {wanted.name}",
            ))
            continue
        if claim.degree < target_degree:
            rejected.append((
                spec.name,
                f"V mismatch: the rule built at {params} is exact only to "
                f"degree {claim.degree} against {wanted.name}, short of "
                f"target_degree={target_degree} — the spec's inversion "
                f"promised a parameter set that its own rule does not "
                f"deliver",
            ))
            continue

        # ---- Stage 3: structural compatibility -----------------------
        spec_flags = spec.structural_flags()
        missing = [
            k for k in requested_flags
            if not spec_flags.get(k, False)
        ]
        if missing:
            rejected.append((
                spec.name,
                f"structural mismatch: rule does not satisfy "
                f"requested flag(s) {missing}",
            ))
            continue

        # ---- Stage 4: cost (minimum points) tracked here -------------
        # The cost is the measure's OWN node count. Stages 0 and 1 have
        # already built it, so a second, formula-shaped source for the same
        # number would be a twin that can only ever drift.
        candidates.append((spec, params, measure.n_points, measure))

    def _selection_log(
        chosen_spec: QuadratureSpec | None,
        chosen_parameters: dict[str, Any] | None,
    ) -> SelectionLog:
        """Build the SelectionLog with the shared selection context.

        One typed construction site for the two outcomes (no-match raise vs
        success): the four context fields are passed with their precise
        types instead of splatting a union-valued ``dict`` — which pyright
        cannot match against the dataclass's distinct field types.
        """
        return SelectionLog(
            chosen_spec=chosen_spec,
            chosen_parameters=chosen_parameters,
            rejected=rejected,
            geometry=geometry,
            angular_symmetry=angular_symmetry,
            target_degree=target_degree,
            requested_flags=requested_flags,
        )

    if not candidates:
        log = _selection_log(chosen_spec=None, chosen_parameters=None)
        raise QuadratureSelectionError(
            f"no registered quadrature satisfies geometry={geometry!r}, "
            f"target_degree={target_degree}, flags={requested_flags}. "
            f"Rejected {len(rejected)} candidates; "
            f"see exc.log.rejected for per-rule reasons.",
            log,
        )

    # Stage 4 tie-break: smallest node count wins. Equal costs are broken
    # by registry order, because `list.sort` is stable.
    #
    # ⛔ This comment claimed registry order "puts the most specialised rule
    # first" until 2026-08-14. Nothing establishes that: the order is
    # GaussLegendre1D, LebedevSphere, LevelSymmetricSN, ProductQuadrature,
    # and no docstring, test, or commit states a specialisation ranking it
    # is meant to encode. `[M]` no test exercises the tie-break at all — no
    # shipped pair of rules produces an equal node count on any geometry —
    # so the policy is both unstated and unpinned. Stability is a real
    # guarantee and is worth keeping; the rationale offered for it was not
    # a claim about this code. If registry order is to MEAN something, it
    # has to be said and gated (cf. LOSS_REPRESENTATIONS, whose own comment
    # states "the registry ORDER is the default-selection policy").
    candidates.sort(key=lambda c: c[2])
    chosen_spec, chosen_params, _, measure = candidates[0]

    log = _selection_log(
        chosen_spec=chosen_spec, chosen_parameters=chosen_params,
    )
    return measure, log


__all__ = [
    "AngularSymmetry",
    "GEOMETRY_ANGULAR_SYMMETRY",
    "QuadratureSelectionError",
    "QuadratureSpec",
    "SelectionLog",
    "quadrature_registry",
    "select_quadrature",
]
