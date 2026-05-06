r"""Tagged quadrature registry with G + V + structural selector.

This module installs the **automated quadrature selector** that closes
the loop on Issues 3 (symmetry tagging) and 4 (rule-side wrappers): given
a geometry and a target polynomial degree, walk the registry of all known
rules and return the cheapest measure whose tags satisfy the geometry's
demands.

Motivation
----------

Quadrature selection is *fundamentally* a four-stage filter, in priority
order:

1. **G compatibility** (symmetry).
   :math:`G_{\text{geom}} \subseteq G_{\text{quad}}` — the rule's
   invariance group must contain every symmetry of the geometry. A
   quadrature with **less** symmetry imprints spurious low-order
   asymmetry on a symmetric problem (Lebedev 1976, §1; Stiefel &
   Fässler 1979 Ch. 5). A quadrature with **more** symmetry is fine —
   the extra symmetry is unused, not violated.
2. **V compatibility** (polynomial exactness, "vague" / Galerkin sense).
   :math:`\deg(Q) \ge d` — the rule's degree of exactness must reach at
   least the target. Most rules have a parameter-dependent degree
   (:math:`2n - 1` for Gauss-Legendre, :math:`N - 1` for level-symmetric
   :math:`S_N`, ``order`` for Lebedev), and the selector inverts this
   to choose the smallest parameter set satisfying the constraint.
3. **Structural compatibility**.
   Boolean flags the consumer can request (``positive_weights``,
   ``axis_aligned``, ``level_structured``, ``half_range_clean``) — the
   rule's flag set must be a superset of the requested set.
4. **Minimum points**.
   Among rules passing 1-3, pick the smallest :math:`N`. This is the
   cost minimisation: each ordinate is a sweep direction in the SN
   solver, and sweeps dominate the runtime.

Formally,

.. math::
   :label: quadrature-selection-criterion

   Q^{\star} \;=\; \arg\min\Bigl\{\, n(Q) \;:\;\;
   G_{\text{geom}} \subseteq G_Q
   \;\wedge\; \deg(Q) \ge d
   \;\wedge\; F_{\text{req}} \subseteq F_Q
   \,\Bigr\},

where :math:`n(Q)` is the number of nodes,
:math:`G_Q \subseteq O(3)` is the invariance group,
:math:`\deg(Q)` is the polynomial-exactness degree, and
:math:`F_Q \subseteq \{\text{positive\_weights}, \text{axis\_aligned},
\text{level\_structured}, \text{half\_range\_clean}\}` is the rule's
structural-flag set.

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
"run this geometry with target degree :math:`d`" prototyping, and a
documentation artifact — every entry's docstring becomes a Sphinx
table row narrating the trade-off.

Geometry → symmetry-group assignment
------------------------------------

The geometries the registry handles today are tagged with the
following native symmetry groups:

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Geometry
     - :math:`G_{\text{geom}}`
     - Rationale
   * - ``"slab"``
     - :math:`SO(2)`
     - 1-D problem in :math:`z`; the angular dependence reduces to
       :math:`\mu = \cos\theta` only, with full azimuthal rotation
       symmetry. The slab's transverse :math:`Z_2` reflection is
       carried by the polar-cosine endpoints
       :math:`\mu \to -\mu` and is automatic for any 1-D rule on
       :math:`[-1, 1]` with symmetric nodes (Gauss-Legendre is
       symmetric by Stoer-Bulirsch §3.6).
   * - ``"sphere"``
     - :math:`SO(2)`
     - The 1-D **radial** spherical SN reduces to GL on
       :math:`\mu_r = \cos\theta_r`, the cosine of the angle between
       the ordinate and the radial direction (Lewis & Miller 1993
       §4.4). The continuous problem has :math:`O(3)` symmetry, but
       the discrete radial reduction collapses azimuthal dependence,
       leaving :math:`SO(2)` as the relevant group.
   * - ``"cylinder"``
     - :math:`SO(2)`
     - Axisymmetric cylinder: :math:`\phi`-independence is
       :math:`SO(2)` about the axial direction. The cylindrical SN
       sweep additionally requires the rule to have **per-:math:`\mu`
       polar-level structure** for the azimuthal redistribution
       coefficients (Bailey et al. 2009 Eq. 50); request this via the
       ``level_structured=True`` structural flag.
   * - ``"cartesian2d"``
     - :math:`O_h`
     - 2-D Cartesian (x-y) carries the full octahedral group when
       both reflection axes are interchangeable. Conservatively
       tagging :math:`O_h` (rather than just :math:`D_{4h}`) means
       the selector accepts any :math:`O_h`-invariant rule —
       Lebedev and level-symmetric :math:`S_N` both qualify.

The judgment call is the spherical case: continuous spherical
geometry is :math:`O(3)` but the **discrete-ordinate radial SN** is a
1-D problem in :math:`\mu_r`, so :math:`SO(2)` is the right group at
the discretisation level. When 2-D / 3-D spherical SN lands, this
table will gain a ``"sphere2d"`` / ``"sphere3d"`` entry tagged
:math:`O_h` or :math:`O(3)`.

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
* Bailey, T.S., Adams, M.L., Yang, B., Zika, M.R. (2009). "A piecewise
  linear discontinuous finite element spatial discretization of the
  transport equation." *Annals of Nuclear Energy* **35**, 1929-1936.
  Eq. 50 (cylindrical α-recursion that needs ``level_structured``).

See Also
--------
:ref:`discrete-measures` (theory page) — the "Quadrature selection
algorithm" section narrates the four-stage precedence chain in
prose.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from ..measure import DiscreteMeasure
from ..symmetry import SubgroupOfO3
from .rules_1d import gauss_legendre_on_mu
from .rules_product import product_mu_phi
from .rules_sphere import lebedev_sphere, level_symmetric_sn


# ---------------------------------------------------------------------------
# Available Lebedev orders (scipy's tabulated subset)
# ---------------------------------------------------------------------------
#
# scipy.integrate.lebedev_rule supports a fixed set of orders. Using an
# unsupported order raises NotImplementedError; the selector must avoid
# those gaps when inverting the degree-of-exactness constraint.
#
# Source: probed at registry-construction time on the SciPy installed
# in this dev container; matches the published Lebedev 1976 / Lebedev &
# Laikov 1999 tabulation.
_LEBEDEV_ORDERS: tuple[int, ...] = (
    3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47,
)


# ---------------------------------------------------------------------------
# Spec dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class QuadratureSpec:
    r"""Tagged registry entry for a single quadrature rule.

    Each spec captures everything :func:`select_quadrature` needs to
    reason about the rule **without** instantiating it: the invariance
    group, the polynomial-exactness degree formula, the structural
    flags, the parameter-inversion logic. The actual measure is built
    lazily by ``factory(**parameters)`` once the selector has decided
    which rule to use.

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
    invariance_group : SubgroupOfO3
        Group :math:`G_Q` under which the rule's measure is invariant.
        Drives the **G compatibility** filter:
        ``geometry_group.is_subgroup_of(invariance_group)``.
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
    expected_node_count : callable
        Maps a parameter dict to the resulting number of nodes. Used
        by the **minimum-points** tie-break stage. Always pure (does
        not call :attr:`factory`) so the selector stays cheap.
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
        redistribution coefficients (Bailey et al. 2009 Eq. 50).
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
    invariance_group: SubgroupOfO3
    parameters: dict[str, type]
    degree_of_exactness_for: Callable[[int], dict[str, Any] | None]
    expected_node_count: Callable[[dict[str, Any]], int]
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


def _gl1d_node_count(params: dict[str, Any]) -> int:
    return int(params["n"])


def _lebedev_invert(target_degree: int) -> dict[str, Any] | None:
    """Lebedev: ``deg = order``; pick the smallest tabulated order.

    SciPy supports a fixed set of Lebedev orders (see
    :data:`_LEBEDEV_ORDERS`). If no tabulated order reaches
    ``target_degree``, the rule cannot satisfy the constraint and
    the inversion returns ``None`` (the selector then rejects this
    spec at the V-compatibility stage with a clear message).
    """
    if target_degree < 0:
        return None
    for order in _LEBEDEV_ORDERS:
        if order >= target_degree:
            return {"order": order}
    return None


# Lebedev node count per order — taken from the SciPy tabulation so
# the selector knows how cheap each candidate is without instantiating.
_LEBEDEV_NODE_COUNTS: dict[int, int] = {
    3: 6, 5: 14, 7: 26, 9: 38, 11: 50, 13: 74, 15: 86, 17: 110,
    19: 146, 21: 170, 23: 194, 25: 230, 27: 266, 29: 302, 31: 350,
    35: 434, 41: 590, 47: 770,
}


def _lebedev_node_count(params: dict[str, Any]) -> int:
    return _LEBEDEV_NODE_COUNTS[int(params["order"])]


def _ls_sn_invert(target_degree: int) -> dict[str, Any] | None:
    r"""Level-symmetric :math:`S_N`: conservative ``deg = N - 1``.

    Carlson-Lathrop's level-symmetric rule has parameter-dependent
    exactness; the rule-side wrapper records ``degree_of_exactness =
    N - 1`` as a conservative lower bound (see
    :func:`~orpheus.numerics.quadrature.rules_sphere.level_symmetric_sn`
    docstring). To meet ``target_degree``, choose the smallest even
    :math:`N \ge \mathrm{target\_degree} + 1`.

    The construction supports any even :math:`N \ge 2`. There is no
    upper bound here, but in practice :math:`N \ge 24` runs into
    moment-condition algebra that the simple equal-weight construction
    in this module does not capture; consumers asking for very high
    target degrees should prefer Lebedev (which scales cleanly to
    order 47).
    """
    if target_degree < 0:
        return None
    n_min = target_degree + 1
    if n_min < 2:
        n_min = 2
    if n_min % 2 == 1:
        n_min += 1
    return {"sn_order": n_min}


def _ls_sn_node_count(params: dict[str, Any]) -> int:
    """Level-symmetric :math:`S_N` total node count.

    Carlson-Lathrop construction: per octant there are
    :math:`(N/2)(N/2 + 1)/2 = N(N + 2)/8` ordinates; total over the
    8 octants is :math:`N(N + 2)`. Verified against the rule-side
    constructor for :math:`N = 4, 8, 16` (4·6=24, 8·10=80, 16·18=288).
    """
    n = int(params["sn_order"])
    return n * (n + 2)


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


def _product_node_count(params: dict[str, Any]) -> int:
    return int(params["n_mu"]) * int(params["n_phi"])


# ---------------------------------------------------------------------------
# Registry population
# ---------------------------------------------------------------------------


quadrature_registry: list[QuadratureSpec] = [
    QuadratureSpec(
        name="GaussLegendre1D",
        factory=gauss_legendre_on_mu,
        invariance_group=SubgroupOfO3.SO2,
        parameters={"n": int},
        degree_of_exactness_for=_gl1d_invert,
        expected_node_count=_gl1d_node_count,
        positive_weights=True,
        axis_aligned=True,        # 1-D nodes ARE the axis (μ-axis)
        level_structured=False,    # no per-mu sub-levels (it IS the mu axis)
        half_range_clean=True,     # restrict(mu>0) gives n/2-point GL on [0,1]
    ),
    QuadratureSpec(
        name="LebedevSphere",
        factory=lebedev_sphere,
        invariance_group=SubgroupOfO3.OctahedralOh,
        parameters={"order": int},
        degree_of_exactness_for=_lebedev_invert,
        expected_node_count=_lebedev_node_count,
        positive_weights=True,
        axis_aligned=True,         # Lebedev grids include axis-aligned orbits
        level_structured=False,    # no polar-level structure (Oh orbits)
        half_range_clean=False,    # half-sphere restriction is not a valid Lebedev
    ),
    QuadratureSpec(
        name="LevelSymmetricSN",
        factory=level_symmetric_sn,
        invariance_group=SubgroupOfO3.OctahedralOh,
        parameters={"sn_order": int},
        degree_of_exactness_for=_ls_sn_invert,
        expected_node_count=_ls_sn_node_count,
        positive_weights=True,
        axis_aligned=False,        # LS_N nodes never lie on coordinate axes
        level_structured=True,     # canonical N/2 polar levels per hemisphere
        half_range_clean=True,     # one hemisphere = exactly half the rule (Oh symmetry)
    ),
    QuadratureSpec(
        name="ProductQuadrature",
        factory=product_mu_phi,
        invariance_group=SubgroupOfO3.SO2,
        parameters={"n_mu": int, "n_phi": int},
        degree_of_exactness_for=_product_invert,
        expected_node_count=_product_node_count,
        positive_weights=True,
        axis_aligned=False,        # depends on phi grid; conservatively False
        level_structured=True,     # n_mu polar levels by construction
        half_range_clean=True,     # GL polar half restricted cleanly to mu>0
    ),
]


# ---------------------------------------------------------------------------
# Geometry → group lookup
# ---------------------------------------------------------------------------
#
# Static table — see the module docstring "Geometry → symmetry-group
# assignment" subsection for the rationale of each entry. New geometries
# (hexagonal lattice, 2-D / 3-D spherical, …) are added here, not in the
# selector itself.

GEOMETRY_GROUPS: dict[str, SubgroupOfO3] = {
    "slab": SubgroupOfO3.SO2,
    "sphere": SubgroupOfO3.SO2,
    "cylinder": SubgroupOfO3.SO2,
    "cartesian2d": SubgroupOfO3.OctahedralOh,
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
    geometry_group : SubgroupOfO3
        The native symmetry group resolved from
        :data:`GEOMETRY_GROUPS`.
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
        string identifies the failing stage (G / V / structural) and
        the missing predicate.
    """

    geometry: str
    geometry_group: SubgroupOfO3
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
    registry: list[QuadratureSpec] | None = None,
    **structural: bool,
) -> tuple[DiscreteMeasure, SelectionLog]:
    r"""Pick the cheapest registered quadrature satisfying every constraint.

    Implements :eq:`quadrature-selection-criterion` over the global
    :data:`quadrature_registry` (or a caller-supplied ``registry`` for
    testing). The four-stage filter — G, V, structural, minimum-points —
    runs in priority order and produces an explainability log alongside
    the chosen measure.

    Parameters
    ----------
    geometry : str
        Geometry tag — one of the keys of :data:`GEOMETRY_GROUPS`
        (``"slab"``, ``"sphere"``, ``"cylinder"``, ``"cartesian2d"``).
        Case-sensitive. Unknown tags raise :class:`KeyError` immediately
        (this is a programming error, not a selection failure).
    target_degree : int
        Minimum polynomial-exactness degree :math:`d` required.
        :math:`d = 0` accepts any rule that integrates constants
        exactly.
    registry : list[QuadratureSpec], optional
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
        If ``geometry`` is not in :data:`GEOMETRY_GROUPS`.
    QuadratureSelectionError
        If no registered spec satisfies every stage. The exception's
        ``log`` attribute carries the full :class:`SelectionLog` for
        introspection.
    """
    if geometry not in GEOMETRY_GROUPS:
        raise KeyError(
            f"unknown geometry {geometry!r}; expected one of "
            f"{sorted(GEOMETRY_GROUPS.keys())}"
        )
    if registry is None:
        registry = quadrature_registry

    geom_group = GEOMETRY_GROUPS[geometry]

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
    candidates: list[tuple[QuadratureSpec, dict[str, Any], int]] = []

    for spec in registry:
        # ---- Stage 1: G compatibility ---------------------------------
        if not geom_group.is_subgroup_of(spec.invariance_group):
            rejected.append((
                spec.name,
                f"G mismatch: geometry {geom_group.name} is not a subgroup "
                f"of rule's invariance group {spec.invariance_group.name}",
            ))
            continue

        # ---- Stage 2: V compatibility ---------------------------------
        params = spec.degree_of_exactness_for(target_degree)
        if params is None:
            rejected.append((
                spec.name,
                f"V mismatch: rule cannot reach target_degree="
                f"{target_degree} with any supported parameters",
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
        n_nodes = spec.expected_node_count(params)
        candidates.append((spec, params, n_nodes))

    log_template = dict(
        geometry=geometry,
        geometry_group=geom_group,
        target_degree=target_degree,
        requested_flags=requested_flags,
    )

    if not candidates:
        log = SelectionLog(
            chosen_spec=None,
            chosen_parameters=None,
            rejected=rejected,
            **log_template,
        )
        raise QuadratureSelectionError(
            f"no registered quadrature satisfies geometry={geometry!r}, "
            f"target_degree={target_degree}, flags={requested_flags}. "
            f"Rejected {len(rejected)} candidates; "
            f"see exc.log.rejected for per-rule reasons.",
            log,
        )

    # Stage 4 tie-break: smallest node count wins. Ties broken by
    # registry order (stable sort), which puts the most specialised
    # rule first when costs match.
    candidates.sort(key=lambda c: c[2])
    chosen_spec, chosen_params, _ = candidates[0]

    # The registered factories return either a DiscreteMeasure or a
    # (DiscreteMeasure, LevelStructure) tuple. Unpack uniformly.
    result = chosen_spec.factory(**chosen_params)
    if isinstance(result, tuple):
        measure = result[0]
    else:
        measure = result

    log = SelectionLog(
        chosen_spec=chosen_spec,
        chosen_parameters=chosen_params,
        rejected=rejected,
        **log_template,
    )
    return measure, log


__all__ = [
    "GEOMETRY_GROUPS",
    "QuadratureSelectionError",
    "QuadratureSpec",
    "SelectionLog",
    "quadrature_registry",
    "select_quadrature",
]
