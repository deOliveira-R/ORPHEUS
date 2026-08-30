r"""The angular factor — :math:`\alpha`, the dome, and the redistribution structure.

The curvilinear redistribution operator **factors**,

.. math::

   \mathcal{R} \;=\; R_{\rm spatial} \;\otimes\; A_{\rm angular}(\tau,\ \alpha,\ w),

and this module is the **member-independent half of the angular factor**: the
:math:`\alpha` dome and the starting direction, which every angular-closure
member reads and none of them owns.  :mod:`~orpheus.sn.angular.closure` is the
member-DEPENDENT half — the Morel--Montry weight :math:`\tau` and the
half-angle grid it interpolates on.

:func:`angular_redistribution` is **THE** single producer.  `[M]` both of its
outputs are functions of ``(quadrature, coord)`` alone — no cell, no mesh, no
material enters — which is why the object can be built on the d≥2 Cartesian
path with no mesh at all
(:mod:`~orpheus.sn.mesh.augmented_mesh`, the ``else:`` arm).

**Why this is S**\ :sub:`N`\ **-specific, and therefore why it lives in**
``sn/``.  An angular redistribution term exists only where an angular unknown
is indexed in a *local rotating frame* **and** its angular derivative is
*collocated*.  MoC fails the middle condition — :math:`\Omega` is fixed in the
global frame, so :math:`\Omega\cdot\nabla = d/ds` is chart-free and the
curvature moves into track segmentation — and CP fails the first, since angle
is integrated out into the kernel before anything is discretised.  P\ :sub:`n`
is the sharp case: it needs the *spatial* factor :math:`R` but **no**
:math:`\alpha` at all, because its angular derivative folds into the Legendre
recurrence and no differencing coefficient is ever formed.  ⟹ :math:`\alpha` is
S\ :sub:`N`'s; :math:`R` is shared with any method that discretizes an angular
variable in a rotating frame.  Full argument and its falsification conditions:
:doc:`/theory/foundations/structured_geometry`, "Who needs a connection
coefficient — and who does not".

**History.**  This module was carved out of ``orpheus/geometry/reduced_operator.py``
on 2026-08-28 (the un-weld arc, P4.2), which is the *second* half of a move the
P2 phase could not make: while the three ``*_streaming`` factories still lived
in ``geometry/`` they CALLED :func:`angular_redistribution` at runtime, so
moving it created a ``geometry -> sn`` import — a declared ``FORBIDDEN_EDGES``
violation *and* a measured circular-import failure of ``import orpheus.geometry``.
P4.4 moved those factories to :mod:`orpheus.sn.mesh.reduced_operator` first,
which turned the call into ``sn -> sn`` and made this landing trivial.
⚠ :class:`AngularMeasure` travelled with :math:`\alpha` rather than with the
factories, because :class:`AngularRedistribution` names it as a field type —
leaving it behind would have been the same forbidden edge in miniature.

⚠ **Whether :class:`AngularMeasure` should EXIST is a separate, open question**
(CS5: should ``DiscreteMeasure`` branch into Spatial / Angular / Energy?).  This
module decides only where it *lives*; moving it prejudges nothing.

Mathematical content
====================

The :math:`\alpha`-dome recursion (sphere) — Hébert (2009)
§3.9.4 Eqs. 3.423-3.424, in the ORPHEUS factor-of-2-absorbed
normalization:

.. math::
   :label: alpha-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

(Hébert's Eq. 3.424 reads
:math:`\alpha^{H}_{n+1/2} = \alpha^{H}_{n-1/2} - 2\mathcal{W}_n\,\mu_n`
with the corresponding redistribution divisor :math:`\Delta S_i /
(2\,\mathcal{W}_n)` in Eq. 3.428; the ORPHEUS arrays carry
:math:`\alpha^{O} = \alpha^{H}/2` so the factor of 2 is absorbed
into the recurrence and the redistribution divisor reads
:math:`\Delta A_i / w_n`.  Both forms are mathematically equivalent.)

For Gauss--Legendre quadrature with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome that closes
back to zero at the upper boundary by GL antisymmetry.

The cylindrical analog is per-:math:`\mu`-level: each level :math:`p`
carries its own :math:`(M+1)`-tuple of :math:`\alpha` with
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine.

The Morel--Montry angular closure weight (Bailey-Morel-Chang 2010
Eq. 43 = Lathrop 2000 Eq. 23) is the *barycentric coordinate* of the
ordinate :math:`\mu_m` between the two edges of its own angular cell
:math:`[\mu_{m-1/2}, \mu_{m+1/2}]`:

.. math::

   \tau_m = \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}}.

It is the UNIQUE choice exact for a flux affine in the radial cosine,
with admissible range :math:`\tau \in [0, 1]`.  The formula carries **no
geometry**: the geometry lives entirely in the *angular cell partition*
it reads.  The canonical definition-of-record and the partition are on
the theory page — see :eq:`morel-montry-closure` (the closure) and
:eq:`angular-cell-partition` (the partition) in
:doc:`/theory/foundations/structured_geometry`, the retired-absorber
account at :ref:`sn-tau-absorber-retirement`, and the comprehensive W1
vindication at :ref:`sn-curvilinear-aniso-norm-reconciliation` in
:doc:`/theory/verification/sn`.

**No clamp, in either arm** (sphere unclamped at W1, 2026-06-13;
cylinder at Q5.6.4, 2026-08-11).  On Gauss--Legendre the sphere's
:math:`\tau_m \in [0.39, 0.61]`; on the shipped folded cylinder arc
:math:`\tau_m \in [\tfrac14, \tfrac34]`.  The former
:math:`[\tfrac12, 1]` clamp was an over-conservative positivity floor
(mis-cited to Lewis & Miller §4.5) on the sphere, and on the cylinder a
compensation for a *wrong cell partition* — the retired chord-midpoint
edges, which disagreed with :math:`\alpha`'s half-angle edges by
~17.5 % in :math:`\omega`-width.  ``morel_montry_tau_raw_per_level``
retired with it: there is one :math:`\tau`.

⚠ This module no longer produces :math:`\tau` at all — #236 Step C
excised the geometry-side producer.  :math:`\tau` is single-sourced in
:func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`,
reading
:func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`.
What stays here is genuinely geometric: the :math:`\alpha`-dome, the
redistribution factor, the face areas, and the starting-direction edges.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 — **§3.9.3
  (cylinder, printed pp. 137-141)** and **§3.9.4 (sphere, printed pp.
  141-144)**; the whole Eq. 3.418-3.439 range is *spherical*.  The
  authority for the curvilinear S\ :sub:`N` **cell balance, DD difference
  relations, sweep ordering and Carlson starting direction** — i.e. for
  everything this module still produces.  He defines **no**
  :math:`\tau` anywhere in chapter 3, in either geometry, so he is
  **not** an authority for the weighted angular closure above; that
  mis-citation produced a wrong cylinder :math:`\tau` and is recorded at
  ``docs/theory/methods/sn/curvilinear_one_group.rst
  §sn-citation-corrections``.  Local copy:
  ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Morel, J. E., & Montry, G. R. (1984).  *Analysis and Elimination of
  the Discrete-Ordinates Flux Dip*.  Transport Theory and Statistical
  Physics 13(5):615-633, doi:10.1080/00411458408211661.  **PRIMARY**
  source for the weighted angular closure :math:`\tau` above.  Local
  copy: ``scratch/literature/Morel-Montry(1984)Analysis and elimination
  of the discrete-ordinates flux dip.pdf``.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *The Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356; OA at
  https://www.osti.gov/servlets/purl/1020346).  **Eqs. (42)/(43)** are
  the form implemented: the Morel--Montry weight :math:`\tau_m` above as
  the unique weight exact for a flux affine in the radial cosine,
  admissible range :math:`\tau \in [0, 1]` — the W1 source for dropping
  the spherical clamp.  Their **Eq. (41)** is the first-order
  diffusion-limit condition :math:`\beta = 0`, and forcing it to zero is
  what DETERMINES the weights; the paper frames that
  :math:`\beta`-contamination analysis via formal-:math:`\varepsilon`
  expansion, and its Eq. 53 names :math:`\tau = \tfrac12` "the diamond
  scheme".
* Reed, W. H., & Lathrop, K. D. (1970).  *Truncation Error Analysis of
  Finite Difference Approximations to the Transport Equation*.  NSE
  41(2):237-248, doi:10.13182/NSE70-A20710.  Their **Eq. (13c)** IS BMC
  Eq. (43), forty years earlier; their **Eqs. (15)/(16)** give the
  sharpest accuracy criterion on :math:`\tau` — second order iff the
  ordinate is the :math:`\mu`-midpoint of its own cell to
  :math:`O(w^2)`, i.e. iff :math:`\tau = \tfrac12 + O(w)`.  Unlike
  :math:`\beta` that criterion is POINTWISE.  Local copy:
  ``scratch/literature/07-Truncation Error Analysis of Finite Difference
  Approximations to the Transport Equation.pdf``.

  *Citation correction (Issue #168 Phase B)*: this module's
  pre-Phase-B docstring cited "Bailey, T. S., Adams, M. L., Yang, B.,
  & Zika, M. R. (2009).  *A piecewise linear finite element
  discretization of the diffusion equation for arbitrary polyhedral
  grids*. JCP 227, 3738-3757", which is the **wrong** Bailey paper —
  Bailey-Adams-Yang-Zika is a piecewise-linear FE diffusion paper
  unrelated to curvilinear S\ :sub:`N` :math:`\alpha`-recursion.  The
  intended reference is Bailey-Morel-Chang 2010 (above).
* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §4.5 — angular redistribution closure.  (NB:
  the historical :math:`[1/2, 1]` clamp was MIS-cited to this §4.5 —
  L&M §4.5 does not prescribe it; W1 traced the exact-on-linear weight
  to Bailey-Morel-Chang 2010 Eq. 43.  ⭐ The interval's real origin is
  Grant, I. P. (1968), *J. Comp. Phys.* 2(4):381-402,
  doi:10.1016/0021-9991(68)90044-2, where :math:`[\tfrac12, 1]` bounds
  the **SPATIAL** weighted-diamond parameter; Reed & Lathrop's footnote
  8 says so and adds that "Grant does not determine angular weights".
  ⚠ Grant 1968 is not in the local library and has not been read.)


See also
========

* :mod:`orpheus.sn.angular.closure` — the member-dependent half of the angular
  factor: the Morel--Montry weight :math:`\tau`, the angular cell partition,
  and the closure family that consumes :class:`AngularRedistribution`.
* :mod:`orpheus.sn.mesh.reduced_operator` — the chart connection; its three
  factories are :func:`angular_redistribution`'s only production callers.
* :mod:`orpheus.derivations.discrete.sn.angular_differencing` — the L0 P0/P4
  predicate ladder, which grades measures whose dome does **not** close.  It is
  ``derivations`` (L0) and may not import ``sn`` (L3), so it ACCEPTS
  :math:`\alpha` as a keyword argument rather than reaching up for it — the
  same treatment ``tau`` and ``edges`` received on 2026-08-12.  That is also
  why :func:`alpha_dome` and :func:`_assert_alpha_dome_closes` are two
  functions: a guard welded into the recursion would make that analysis
  unspellable.
* :doc:`/theory/foundations/structured_geometry` — the architecture page.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import numpy as np

from orpheus.geometry.coord import CoordSystem

if TYPE_CHECKING:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.quadrature.rules_sphere import LevelStructure



# ═══════════════════════════════════════════════════════════════════════
# Angular-measure adapter (Issue 4 contract — local, structural)
# ═══════════════════════════════════════════════════════════════════════

class AngularMeasure(Protocol):
    """Minimal contract this module needs from an angular quadrature.

    The concrete implementer is
    :class:`~orpheus.numerics.quadrature.Quadrature`, which wraps a
    :class:`~orpheus.numerics.measure.DiscreteMeasure`; here we only
    depend on the attributes the connection-coefficient math reads, so
    the geometry layer needs no import from the quadrature package at
    all.  This structural Protocol outlived the type it was written
    against: the ``AngularQuadrature`` Protocol and its four
    per-family adapter classes at ``orpheus.sn.quadrature`` were
    retired into classmethod factories on the single ``Quadrature``
    type (R-1 Phase A detour-C), and the read contract below did not
    have to move.
    """

    # Read-only properties (not plain attributes): this module only READS
    # the angular measure, and concrete implementers (``Quadrature``) expose
    # these as ``@property``. A mutable attribute in the Protocol is invariant
    # and rejects a read-only ``@property`` ("property is not assignable to
    # ndarray"); a read-only Protocol member accepts both a property and a
    # plain attribute, matching the actual usage (covariant read contract).

    @property
    def mu_x(self) -> np.ndarray:
        """(N,) primary direction cosine — :math:`\\mu` for sphere,
        :math:`\\eta` (radial) for cylindrical product quadrature."""
        ...

    @property
    def weights(self) -> np.ndarray:
        """(N,) quadrature weights."""
        ...

    @property
    def N(self) -> int:
        """Number of ordinates."""
        ...

    @property
    def eta(self) -> np.ndarray:
        r"""(N,) cylindrical-frame radial direction cosine
        :math:`\eta = \Omega \cdot \hat{r}` (read by the curvilinear
        connection-coefficient recursion only)."""
        ...

    @property
    def mu_z(self) -> np.ndarray:
        """(N,) axial direction cosine (curvilinear / level-structure side)."""
        ...

    @property
    def level_indices(self) -> list[np.ndarray]:
        """Per-:math:`\\mu`-level ordinate-index partition: ``[arange(N)]`` for
        a slab quadrature (single level), one index array per level for
        cylindrical-compatible cubatures."""
        ...

    @property
    def level_structure(self) -> "LevelStructure | None":
        """The per-level side-channel (``None`` for slab / 2-D rules).

        Declared here since CS5 (2026-08-29): the cylinder factory was
        already reading it past the contract via
        ``getattr(angular_measure, "level_structure", None)`` — the
        Protocol under-declared what its consumer needs. The tolerant
        ``getattr`` at the probe survives (structural conformers may
        predate this declaration); the CONTRACT now tells the truth.
        """
        ...

    def axis(self, label: str = "angular") -> "Axis":
        """Mint the generator-stamped angular space-factor axis.

        Declared since the P4-remainder (2026-08-29), the same repair
        shape as ``level_structure`` above: the three streaming factories
        MINT the axis from the measure they are handed, so the mint is
        part of what this contract's consumers need — a Protocol that
        omitted it would under-declare again. The concrete implementer is
        :meth:`Quadrature.axis <orpheus.numerics.quadrature.directional.Quadrature.axis>`.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# The α-dome — ONE recursion, and its admission contract
# ═══════════════════════════════════════════════════════════════════════
#
# Both curvilinear factories below need the Lathrop--Carlson angular
# redistribution coefficients, and until 2026-08-12 each spelled the
# recursion out in its own body (sphere and cylinder, identical arithmetic
# — the sphere IS the cylinder's single-level case).  That twin path is
# why the closure contract was enforced on ONE arm only, and why the one
# check that did exist was a bare ``assert`` (see below).  Cardinal Rule 2:
# one concept, one body.
#
# The RECURSION and the CONTRACT are deliberately two functions, not one.
# ``orpheus.derivations.discrete.sn.angular_differencing`` explores
# hypothetical and deliberately-non-closing measures (its P0/P4 predicate
# ladder) and must be able to run the recursion on them; a guard welded
# into the recursion would make that analysis unspellable.  The contract
# belongs at the production ADMISSION point — the factories — where a
# non-closing measure is genuinely illegal.


def alpha_dome(
    cosines: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""The Lathrop--Carlson angular redistribution recursion.

    .. math::
       \alpha_{m+1/2} \;=\; \alpha_{m-1/2} \;-\; w_m\,\mu_m ,
       \qquad \alpha_{1/2} = 0

    in the ORPHEUS factor-of-2-absorbed normalization (Hébert 2009
    §3.9.4 Eqs. 3.423-3.424 for the sphere; his cylinder is §3.9.3, where
    the :math:`\eta_{p,q\pm1/2}` construction is credited to Alcouffe &
    O'Dell 1986 — not local, not read).  The recursion itself traces to
    Lathrop & Carlson, *J. Comp. Phys.* **1**:173 (1966) — likewise not
    local, not read; cited by Reed & Lathrop 1970 (their ref. 7) as "a
    requirement commonly invoked to define the α coefficients".

    ⭐ **The recursion is strictly ONE-SIDED**: :math:`\alpha_{1/2} = 0` is
    an axiom of the construction, but the closure at the far end,
    :math:`\alpha_{M+1/2} = 0`, is **not** — it is a *consequence* of the
    measure's antisymmetry, i.e. a genuine contract on every admitted
    quadrature, and nothing here enforces it.  Callers that require a
    physically admissible dome must run
    :func:`_assert_alpha_dome_closes` on the result.

    Parameters
    ----------
    cosines :
        Shape ``(M,)`` — the marching direction cosine per ordinate,
        **ordered by increasing cosine**: :math:`\mu` on the sphere, the
        radial :math:`\eta` within one :math:`\mu`-level on the cylinder.
    weights :
        Shape ``(M,)`` — the matching quadrature weights (the level's
        azimuthal share on the cylinder).

    Returns
    -------
    np.ndarray
        Shape ``(M+1,)``, the dome on the angular EDGES, with
        ``alpha[0] == 0.0`` exactly.
    """
    cosines = np.asarray(cosines, dtype=float)
    weights = np.asarray(weights, dtype=float)
    alpha = np.zeros(cosines.size + 1)
    for m in range(cosines.size):
        alpha[m + 1] = alpha[m] - weights[m] * cosines[m]
    return alpha


#: Closure tolerance for :func:`_assert_alpha_dome_closes`.  `[M]` 2026-08-12
#: the worst residual over the shipped rules is ``2.78e-16``
#: (``folded_product(4, 32)``); Gauss--Legendre reads ``5.6e-17`` (N=4)
#: … ``8.2e-17`` (N=64), i.e. ≈ 1 ULP of the dome peak and not drifting
#: with N.  ``1e-12`` therefore clears every shipped rule by ~4 orders
#: while still refusing any real antisymmetry violation.
_ALPHA_CLOSURE_ATOL = 1e-12


def _assert_alpha_dome_closes(
    alpha: np.ndarray,
    *,
    coord: CoordSystem,
    level: int | None = None,
) -> None:
    r"""The admission contract :math:`\alpha_{M+1/2} = 0`.

    A dome that does not close means the measure's marching cosines are
    not antisymmetric about the level's centre — a mis-ordered, truncated,
    duplicated or simply inadmissible rule.  The consequence is not a
    small error, and it is not local: the closure's coefficient
    precompute reads :math:`\alpha_{M+1/2}` as the LAST entry of its
    ``alpha_out`` slice, so it lands in **both** cell-balance
    coefficients of the final ordinate —
    ``c_out[M-1] = α_{M+1/2}/τ`` (a denominator term) and
    ``c_in[M-1] = (1-τ)/τ·α_{M+1/2} + α_{M-1/2}`` (an
    upstream-numerator term); see
    :class:`~orpheus.sn.angular.closure.MorelMontryAngularSweep`'s
    ``_c_in_per_level`` / ``_c_out_per_level`` precompute.  A closing dome
    makes that denominator term vanish — angular redistribution stops at
    the level's top edge, which is the whole point of the edge.  A
    non-zero value instead redistributes flux past it, into nothing.
    Refuse it here.

    ⛔ **This guard replaced a bare ``assert`` on 2026-08-12, and the
    ``assert`` was INERT where it mattered.** The sphere arm carried
    ``assert abs(alpha[N]) < 1e-12``; the cylinder arm carried nothing at
    all.  The project's canonical runner is ``python -O -m pytest``
    (``.claude/rules/vv-testing.md``), and ``-O`` sets ``__debug__ =
    False`` and strips every ``assert`` statement — so on the one arm that
    had a check, that check did not run in the canonical test suite.
    `[M]` measured on the verbatim recursion: a measure closing at
    ``alpha[N] = +0.2000`` is REFUSED under plain ``python`` and
    **ACCEPTED** under ``python -O``.  A contract that a compiler flag can
    remove is not a contract; this raises.

    Sibling of
    :func:`~orpheus.sn.angular.closure._assert_tau_within_unit_interval`
    — the same shape (an admission predicate on the angular scheme that
    refuses upstream rather than absorbing downstream), for the same
    reason.
    """
    residual = float(abs(alpha[-1]))
    if not residual < _ALPHA_CLOSURE_ATOL:
        where = "" if level is None else f" on level {level}"
        raise ValueError(
            f"the alpha dome does not close{where}: "
            f"alpha[M+1/2] = {float(alpha[-1])!r}, which exceeds the "
            f"{_ALPHA_CLOSURE_ATOL:g} admission tolerance for a "
            f"{coord.name} rule. The recursion "
            f"alpha_{{m+1/2}} = alpha_{{m-1/2}} - w_m*mu_m starts at "
            f"alpha_1/2 = 0 and can only return there if the marching "
            f"cosines are antisymmetric about the level's centre, so a "
            f"non-zero endpoint certifies an inadmissible measure — "
            f"mis-ordered, truncated, or duplicated ordinates upstream. "
            f"alpha[M+1/2] multiplies the last ordinate's outflow "
            f"half-face, so absorbing it here would leak flux out of the "
            f"level's top edge. Fix the quadrature upstream."
        )


# ═══════════════════════════════════════════════════════════════════════
# The angular factor — the redistribution structure of a measure in a chart
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class AngularRedistribution:
    r"""The angular measure's redistribution structure on a coordinate chart.

    **The object exists because the redistribution operator FACTORS.**  The
    curvilinear angular-redistribution term of the transport operator is a
    tensor product

    .. math::

       \mathcal{R} \;=\; R_{\rm spatial} \;\otimes\;
                          A_{\rm angular}(\tau,\ \alpha,\ w)

    — a spatial pairing against an angular operator — and this type is the
    **member-independent half of the angular factor**: the :math:`\alpha`
    dome and the starting direction, which every angular-closure member
    shares.  What is NOT here is the member's own choice: :math:`\tau` and
    the derived :math:`c_{\rm in}` / :math:`c_{\rm out}` belong to the
    closure (Morel--Montry's barycentric weight, plain diamond's
    :math:`\tfrac12`, the neutral :math:`\tau \equiv 1`), and a shared
    object holding them would forbid a second member by construction.

    (Derived 2026-08-26, ``scratch/tau_under_ld_dip_analysis.md`` §3 and
    §7: the factorization is why the angular closure's weight cannot
    acquire spatial content, and the member-independent /
    member-dependent split is that memo's §7.1/§7.2 table.)

    **One owner, two consumers.**  Until this type existed, :math:`\alpha`
    and :math:`\mu_{\rm start}` were stored as six ``Optional`` fields on
    :class:`ReducedStreamingOperator` — a per-coordinate union discharged
    by ``assert`` — and were read both by the angular closure and by
    :meth:`ReducedStreamingOperator.streaming_terms`.  They are produced
    HERE, once, by :func:`angular_redistribution`.

    **Cartesian is the NEUTRAL element, not a special case.**  A slab has
    no curvature, so its dome is identically zero and its starting
    direction is the diameter ray; spelling those values rather than
    ``None`` is what lets the per-coordinate union die (Pattern 4 — the
    "no redistribution" state is not separately representable).

    **Sphere is the one-level case.**  Every field is
    per-:math:`\mu`-level; the sphere carries exactly one level.  The
    closure's own constructor already normalized to this shape internally,
    so the normal form is moved here rather than invented.

    Attributes
    ----------
    coord :
        The chart the redistribution was built on.
    alpha_per_level :
        The Lathrop--Carlson dome :math:`\alpha_{m\pm 1/2}` per level,
        each of shape ``(M_p + 1,)`` and closing at
        :math:`\alpha_{M+1/2} = 0` (:func:`_assert_alpha_dome_closes`).
    mu_start_per_level :
        The starting-direction angular edge per level — the direction the
        half-angle thread's seed flux lives at.  Sphere: ``-1.0`` (the
        Hébert §3.9.4 diameter ray).  Cylinder:
        :math:`-\sin\theta_p = -\sqrt{1-\xi_p^2}`.
    """

    coord: CoordSystem
    alpha_per_level: tuple[np.ndarray, ...]
    mu_start_per_level: tuple[float, ...]

    @property
    def n_levels(self) -> int:
        r"""Number of :math:`\mu`-levels (sphere and slab: ``1``)."""
        return len(self.alpha_per_level)


def angular_redistribution(
    quadrature: AngularMeasure,
    coord: CoordSystem,
) -> AngularRedistribution:
    r"""Build the :class:`AngularRedistribution` for a measure on a chart.

    THE single producer of :math:`\alpha` and :math:`\mu_{\rm start}`.
    Both are functions of ``(quadrature, coord)`` alone — no cell, no
    mesh, no material enters — which is measured, not assumed: every
    value this returns is **bit-identical** to what the three streaming
    factories stored before the carve (`[M]` 2026-08-26, both geometries,
    all levels, ``max|delta| = 0.000e+00``).

    Per level, :math:`\alpha` is :func:`alpha_dome` over that level's
    cosines and weights, and the admission contract
    :math:`\alpha_{M+1/2} = 0` is checked here — at the one site that
    produces it, rather than on whichever arm remembered to.
    """
    if coord is CoordSystem.CARTESIAN:
        # The neutral element: no curvature ⇒ no redistribution.  The dome
        # is zero (so the closure's c_in = c_out = 0 fall out rather than
        # being special-cased) and the starting direction is the diameter.
        n = int(np.asarray(quadrature.mu_x).size)
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=(np.zeros(n + 1),),
            mu_start_per_level=(-1.0,),
        )

    mu_x = np.asarray(quadrature.mu_x)
    weights = np.asarray(quadrature.weights)

    if coord is CoordSystem.SPHERICAL:
        # One level: every ordinate marches on the single polar dome, and
        # the seed rides the diameter ray (Hébert §3.9.4).
        alpha = alpha_dome(mu_x, weights)
        _assert_alpha_dome_closes(alpha, coord=coord)
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=(alpha,),
            mu_start_per_level=(-1.0,),
        )

    if coord is CoordSystem.CYLINDRICAL:
        # One dome per μ-level; the level's seed rides its most-inward
        # azimuthal edge, η_{1/2} = -sin θ_p = -sqrt(1 - ξ_p²).
        alphas: list[np.ndarray] = []
        starts: list[float] = []
        xi = np.asarray(quadrature.mu_z)
        for p, level_idx in enumerate(quadrature.level_indices):
            level_idx = np.asarray(level_idx)
            alpha = alpha_dome(mu_x[level_idx], weights[level_idx])
            _assert_alpha_dome_closes(alpha, coord=coord, level=p)
            alphas.append(alpha)
            starts.append(-float(np.sqrt(1.0 - float(xi[level_idx[0]]) ** 2)))
        return AngularRedistribution(
            coord=coord,
            alpha_per_level=tuple(alphas),
            mu_start_per_level=tuple(starts),
        )

    raise ValueError(  # pragma: no cover — exhaustive match above
        f"Unknown coord system: {coord!r}"
    )


__all__ = [
    "AngularMeasure",
    "AngularRedistribution",
    "alpha_dome",
    "angular_redistribution",
]
