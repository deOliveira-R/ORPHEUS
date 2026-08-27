r"""Angular-redistribution closure strategies for the curvilinear FD operator.

The curvilinear S\ :sub:`N` cell-balance (Hébert 2009 §3.9.4 Eq. 3.428)
carries an angular-redistribution term in **half-angle face fluxes**
:math:`\phi_{n\pm 1/2,i}` — the flux between consecutive ordinate
sub-domains.  They are NOT cell-centre values and cannot be computed
from cell-centres without a **closure**; this module supplies that
closure as a strategy family.

The production closure is the per-cell Morel--Montry **weighted**-DD
angular recurrence
:math:`\phi_{n+1/2,i} = (\phi_{n,i} - (1-\tau_n)\phi_{n-1/2,i})/\tau_n`,
seeded at the Carlson starting direction :math:`\mu = -1` (where the
redistribution weight vanishes, :math:`\alpha_{1/2} = 0`).  It runs the
SAME algebra
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` runs inside
the sweep, lifted to operator level so the apply matvec and the sweep
solve the **same** discrete fixed point.

.. warning::

   ⛔ **The weighted recurrence above is NOT Hébert's.** Until 2026-08-11
   this docstring attributed it to "Hébert Eqs. 3.437 / 3.439", which was
   false three ways and is recorded here because the mis-citation is what
   produced a wrong cylinder :math:`\tau` (Q5.6.4):

   * §3.9.4 is Hébert's **SPHERE**; his cylinder is §3.9.3.  The whole
     3.418--3.439 range is spherical.
   * Eqs. 3.437 / 3.439 are **not weighted** — 3.439 reads
     :math:`\phi_{n+1/2,i} = 2\phi_{n,i} - \phi_{n-1/2,i}`, i.e. Eq. 3.431
     rearranged for the sweep, which is :math:`\tau \equiv \tfrac12`
     exactly.  The cylinder's azimuthal counterparts, 3.412 / 3.414, have
     the identical shape.
   * **Hébert defines no** :math:`\tau` **anywhere in chapter 3**, in
     either geometry.

   Hébert therefore ships the *plain* angular diamond, which this module
   deliberately does **not** use: Bailey--Morel--Chang 2010 prove
   (their Eq. 53 + §I) that the plain diamond is diffusion-limit
   consistent only to LEADING order, while the weighted diamond is the
   only member of the family correct through FIRST order — and that
   first-order consistency is what removes the flux dip in general.
   Never cite Hébert against BMC here.

Since Issue #282 route (a) the
seed :math:`\phi_{1/2,i}` is first-class STATE marched directly from the
source (carrying levels — the GL sphere and, since Q5.6, every level
of a σ_y-folded cylinder), or the inlined 2-point angular-edge
extrapolation on non-carrying cylinder levels (NODE_ALIGNED /
level-symmetric full-circle rules); the retired
pre-route-(a) :math:`\phi_{1/2,i} = 0` back edge survives only as the
:math:`\psi`-independent coefficient state (a ``None`` seed).

The full theory lives in the book (nothing derived here):

* ``docs/theory/methods/sn/curvilinear_one_group.rst §balance-curvilinear``
  — the redistribution term, the :math:`\alpha`-recursion, the
  :math:`\Delta A/w` flat-flux consistency proof, and the M-M weights
  (``mm-weights`` / ``wdd-face`` / ``pole-mm-recurrence``);
* ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-apply-sweep-equivalence``
  — why the apply matvec and the sweep solve the same discrete operator;
* ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-solve``
  and ``§sn-direct-seed-r12a`` — the route-(a) direct seed
  (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`,
  :meth:`MorelMontryAngularSweep.edge_extrapolated_seed`) and which
  levels carry a ψ½ block.

The single strategy contract — :class:`PoleAngularClosureBase`
==============================================================

The angular-closure family is a self-registering ABC,
:class:`PoleAngularClosureBase`, layered on
:class:`~orpheus.numerics.registry.RegistryMixin`.  Each strategy
inherits it (passing ``key="..."``), carries a class-level
``is_linear: bool`` trait, and implements the matvec/sweep contract
(:meth:`~PoleAngularClosureBase.precompute_psi_state` /
:meth:`~PoleAngularClosureBase.cell_contribution` /
:meth:`~PoleAngularClosureBase.angular_adjoint`, abstract on the ABC)
plus the closure-constant accessors
(:attr:`~PoleAngularClosureBase.c_in_per_ordinate` /
:attr:`~PoleAngularClosureBase.c_out_per_ordinate` /
:attr:`~PoleAngularClosureBase.tau_per_ordinate`).  The pole singularity
is **intrinsic geometry** (a coordinate-system singularity), not an
external boundary, so the closure stays a separate concern from the
boundary trace.  One strategy covers both curvilinear geometries via an
optional per-:math:`\mu`-level loop: sphere is the single-level case;
cylinder loops one azimuthal sub-problem per level (each with its own
:math:`\alpha` dome and :math:`\Delta A/w` factor), the per-cell DD
recurrence structurally identical, only the ordinate index list changing.

The contract's evolution — the retired Phase-A ``BoundaryFaceFlux``
Protocol it mirrored, the Protocol→ABC retype (#236 Phase 2 B2), and the
Issue #248 retirement of the divergent ``PoleAngularClosure`` Protocol +
the legacy ``__call__`` bundle + the ``LegacyTauSymmetricInterpolation``
/ ``BaileyFlatFluxRedist`` ablation strategies — is the record at
``docs/theory/methods/sn/curvilinear_one_group.rst §sn-pole-angular-closure-protocol``.
The citation record is a section of its own on the same page,
``§sn-citation-corrections`` — the wrong 2009 Bailey FE-diffusion paper the
pre-Phase-B geometry docstrings cited (Issue #168 Phase B), and the
Hébert-vs-BMC attribution of the weighted :math:`\tau` (Q5.6.4, the warning
above).  The per-object authority table is ``§sn-tau-source-of-record``.

References
==========

Sources are listed by WHAT they are the authority for; the weighted
:math:`\tau`, the :math:`\alpha` recursion and the sweep mechanics come
from three different places, and conflating them is the ERR-class this
module has already paid for once (see the warning above).

* **The weighted** :math:`\tau` **— PRIMARY:** Morel, J. E., & Montry,
  G. R. (1984). *Analysis and Elimination of the Discrete-Ordinates Flux
  Dip*. Transport Theory and Statistical Physics 13(5):615-633.
  doi:10.1080/00411458408211661.  Local copy:
  ``scratch/literature/Morel-Montry(1984)Analysis and elimination of the discrete-ordinates flux dip.pdf``.
* **The weighted** :math:`\tau` **— the form we implement:** Bailey,
  T. S., Morel, J. E., & Chang, J. H. (2010). *Asymptotic Diffusion-Limit
  Accuracy of Sn Angular Differencing Schemes*. NSE 165(2):149-169.
  **Eqs. (42)/(43)** define :math:`\tau_m` as the barycentric coordinate
  of the ordinate between its own cell's two edges **in the radial
  direction cosine**; their Eq. (41) is the first-order diffusion-limit
  condition :math:`\beta = 0`, and forcing it to zero is what DETERMINES
  these weights (so :math:`\tau` is derived, not chosen).  Eq. (12) is the
  sphere's cumulative-weight cell partition, implemented verbatim in
  :func:`angular_cell_edges_per_level`.
* **The same barycentric condition, 40 years earlier:** Reed, W. H., &
  Lathrop, K. D. (1970). *Truncation Error Analysis of Finite Difference
  Approximations to the Transport Equation*. NSE 41:237.  Their Eq. (13c)
  IS BMC Eq. (43).  They additionally impose Eq. (13b), which turns the
  system into a quadratic for the ORDINATE (edges in, ordinates out) — a
  DIFFERENT branch that fixes the quadrature, and not ours.  ⭐ Their
  **Eqs. (15)/(16)** are the sharpest available accuracy criterion on
  :math:`\tau`: the angular truncation error is second order iff the
  ordinate is the :math:`\mu`-MIDPOINT of its own cell to
  :math:`O(w^2)`, i.e. iff :math:`\tau = \tfrac12 + O(w)`.  Unlike
  :math:`\beta` this is POINTWISE, so a :math:`\sigma_y`-folded level does
  not annihilate it.
* **The** :math:`\alpha` **recursion** :math:`\alpha_{m+1/2} =
  \alpha_{m-1/2} - \mu_m w_m`: Lathrop, K., & Carlson, B. (1966). *J.
  Comp. Phys.* 1:173 — cited by Reed & Lathrop (their ref. 7) as "a
  requirement commonly invoked to define the :math:`\alpha`
  coefficients".  Hébert credits the cylindrical
  :math:`\eta_{p,q\pm1/2}` construction to Alcouffe, R. E., & O'Dell,
  R. D. (1986), *Transport Calculations for Nuclear Reactors*, in the CRC
  Handbook of Nuclear Reactors Calculations Vol. I (Y. Ronen, ed.).
  ⚠ **Neither is in the local library and neither has been read.**
* **The sweep mechanics and the Carlson starting direction:** Hébert, A.
  (2009). *Applied Reactor Physics*, Chapter 3 — **§3.9.3 (cylinder,
  printed pp. 137-141)** and **§3.9.4 (sphere, printed pp. 141-144)**.
  Local copy: ``scratch/literature/Hebert(2009)Chapter3.pdf``.  Authority
  for the cell-balance layout, the sweep ordering and the
  :math:`\alpha_{1/2} = 0` initialisation — **NOT** for the weighted
  :math:`\tau` (see the warning above).
* Grant, I. P. (1968). *J. Comp. Phys.* 2(4):381-402,
  doi:10.1016/0021-9991(68)90044-2 — the origin of the weighted-diamond
  ansatz that BOTH branches credit, and of the
  :math:`[\tfrac12, 1]` interval.  ⚠ **That interval is on the SPATIAL
  weight**, not on the angular :math:`\tau` (Reed & Lathrop footnote 8,
  read on the rendered page); transplanting it onto :math:`\tau` is
  exactly what the retired cylinder "absorber" did.  Not in the local
  library.
* ``docs/theory/methods/sn/curvilinear_one_group.rst`` — the pole
  angular closure (Issue #168 Phase B) and the α-recursion crosswalk
  (``docs/theory/conventions/normalization.rst §normalization-alpha-crosswalk``).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from orpheus.geometry import CoordSystem
from orpheus.geometry.reduced_operator import AngularRedistribution
from orpheus.numerics.registry import RegistryMixin

if TYPE_CHECKING:  # pragma: no cover
    from orpheus.transport.fields._bases import RadialCharacteristicInteriorField

    from ..mesh.augmented_mesh import SNMesh


# ═══════════════════════════════════════════════════════════════════════
# PoleAngularClosureBase — concrete ABC with self-registration
# ═══════════════════════════════════════════════════════════════════════
#
# This ABC is the SOLE angular-closure contract (Cardinal Rule 2 — one
# contract per concept). The retired ``@runtime_checkable
# PoleAngularClosure`` Protocol and the legacy ``__call__`` bundle (with its
# ``tau_mm`` argument) are the record's story:
# curvilinear_one_group.rst §sn-pole-angular-closure-protocol.



def _require_single_moment_pairing(pairing: np.ndarray, who: str) -> None:
    r"""Refuse a multi-moment spatial pairing — the cell SOLVE does not exist yet.

    The pairing's two axes are real and load-bearing (see
    :attr:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.redistribution_pairing`):
    ``n_mom`` is the scheme's spatial-moment count, ``n_thread`` is what the
    angular device propagates.  What does NOT exist yet is the other side of
    the contract: at ``n_mom > 1`` the closure's contribution to the cell
    balance stops being a scalar and becomes an ``(n_mom, n_thread)`` block,
    so ``psi_out = numer / denom`` becomes a linear SOLVE — which is the
    linear-discontinuous cell update itself (Issue #158, the curvilinear
    arm; the derivation is published, see that issue's record).

    ⚠ **Read its coverage precisely** (``plan-authoring`` §6c).  No SHIPPED
    scheme is multi-moment, so nothing in production reaches this guard —
    but a multi-moment pairing is a plain array, so a HAND-BUILT one does, and
    that is a real witness rather than a mutation of the code under test.
    The gate is
    ``tests/sn/sweep/curvilinear/test_pole_angular_closure.py``; it pins
    both shapes the literature actually proposes — a square ``(nx, 2, 2)``
    (Adams--Martin, closing per spatial moment) and a rank-1 ``(nx, 2, 1)``
    (ONETRAN, closing on the cell average only).  So this is a fail-loud
    marker for the step that adds a multi-moment scheme, AND it is gated.
    """
    pairing = np.asarray(pairing)
    if pairing.ndim != 3:
        raise ValueError(
            f"{who}: the redistribution pairing must be (nx, n_mom, n_thread); "
            f"got shape {pairing.shape}. Build it from "
            f"ReducedStreamingOperator.redistribution_pairing."
        )
    n_mom, n_thread = pairing.shape[1], pairing.shape[2]
    if n_mom != 1 or n_thread != 1:
        raise NotImplementedError(
            f"{who}: a multi-moment redistribution pairing "
            f"(n_mom={n_mom}, n_thread={n_thread}) needs the multi-moment "
            f"cell SOLVE, which is not implemented (Issue #158, the "
            f"curvilinear linear-discontinuous arm). The pairing's axes are "
            f"accepted here so the pairing is expressible; the cell balance "
            f"below still divides by a scalar."
        )


class PoleAngularClosureBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering pole-angular-closure strategies.

    Subclasses inherit this ABC and pass ``key="..."`` in the class
    statement to self-register; the registry is consulted via
    :meth:`PoleAngularClosureBase.create("morel_montry_angular_sweep")`
    (or any other registered key).

    Subclasses MUST declare:

    * ``is_linear: ClassVar[bool]`` — whether the closure is linear
      in the cell-centre angular flux.
    * :meth:`precompute_psi_state` / :meth:`cell_contribution` /
      :meth:`angular_adjoint` — the matvec/sweep strategy contract (the
      forward angular path, its per-cell cell-balance contribution, and
      its reverse-mode adjoint).
    * ``__init__`` constructible as ``cls(angular, pairing)`` — the family
      construction contract (abstract below): every closure binds to
      its :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` at
      construction (PR-TYPED-6.5 Phase 2.3), and the SNMesh
      default-closure dispatch instantiates through this signature
      (``default_angular_closure_class(coord)(angular, pairing)``). Concretes may
      widen (an optional mesh, extra keyword strategy slots).

    Notes
    -----
    Adding a new strategy is a one-line registration::

        class MyClosure(PoleAngularClosureBase, key="my_closure"):
            is_linear: ClassVar[bool] = True

            def precompute_psi_state(self, psi_view, *,
                                     radial_characteristic=None):
                ...

            def cell_contribution(self, psi_state, cell_idx, level_idx,
                                  within_positions):
                ...

            def angular_adjoint(self, numer_bar):
                ...

    No registry insert; ``PoleAngularClosureBase.create("my_closure")``
    is immediately callable.
    """

    registry: ClassVar[dict[str, type["PoleAngularClosureBase"]]] = {}

    is_linear: ClassVar[bool]

    # ── Instance state set by every concrete subclass' ``__init__`` ──
    # The per-:math:`\mu`-level ordinate partition and the precomputed
    # closure constants.  Annotated here (not assigned) so the polymorphic
    # ``c_in_per_ordinate`` / ``c_out_per_ordinate`` accessors below typecheck
    # against the shared contract; the concrete value is bound in each
    # subclass (M-M from its α-dome / τ, Identity to neutral zeros).  Every
    # closure binds its two tensor factors at construction (``cls(angular, pairing)``
    # family contract), so the state is always populated — the former
    # ``| None`` widenings served only the retired M-M unbound legacy mode.
    level_indices: "tuple[np.ndarray, ...]"
    _c_in_per_level: "tuple[np.ndarray, ...]"
    _c_out_per_level: "tuple[np.ndarray, ...]"
    # The Morel--Montry angular weight ``τ`` per μ-level (Issue #236 Phase 2):
    # the FUNDAMENTAL angular weight the closure owns, from which the two c
    # constants above are derived.  Annotated here so the polymorphic
    # ``tau_per_ordinate`` accessor below typechecks against the shared
    # contract; bound by every concrete ``__init__`` (M-M from
    # ``morel_montry_tau_per_level``, Identity to the neutral ``τ ≡ 1``).
    _tau_per_level: "tuple[np.ndarray, ...]"
    # Cached ``(N,)`` global-ordinate gathers of the per-level constants
    # (Issue #236 Phase 2 B2 Fix 1).  The per-level→global gather is a pure
    # permutation of the immutable ``_c_*_per_level`` / ``_tau_per_level``, so
    # it is computed ONCE in each concrete ``__init__`` (via
    # ``_build_per_ordinate_cache``) and
    # the public ``c_*_per_ordinate`` accessors return the cache — O(1) per
    # access instead of re-running the ``(N,)`` gather on every read (the
    # per-visit ``_make_cell_visit`` stamp would otherwise be O(N²·nx)).
    _c_in_per_ordinate_cache: np.ndarray
    _c_out_per_ordinate_cache: np.ndarray
    # Cached ``(N,)`` global-ordinate gather of the per-level Morel--Montry
    # angular weight ``τ`` (Issue #236 Phase 2 B3).  τ is the FUNDAMENTAL
    # angular weight from which ``c_out = α_out/τ`` and
    # ``c_in = (1−τ)/τ·α_out + α_in`` are derived; the live sweep
    # (``DiamondDifference.update``) and the CumprodScan fast path consume τ
    # (the closure's owned weight) rather than the former geometry-factory
    # ``StreamingTerms.tau_mm`` (retired in Step C).  Same caching rationale
    # as the c-caches: the
    # per-level→global gather is a pure permutation, computed ONCE in each
    # concrete ``__init__`` (via ``_build_per_ordinate_cache``).
    _tau_per_ordinate_cache: np.ndarray

    beta_first_order_consistent: ClassVar[bool] = False
    r"""Whether this angular redistribution closure satisfies the
    Bailey–Morel–Chang (2010) first-order diffusion-limit condition
    :math:`\beta = 0` (BMC Eq. (41), the angular functional
    :math:`\sum_m \mu_m[\alpha_{m+1/2}\mu_{m+1/2} - \alpha_{m-1/2}\mu_{m-1/2}]`).
    This is the ANGULAR half of the (spatial ⊗ angular) diffusion-limit
    condition; the SPATIAL half is the discretization scheme's
    ``diffusion_limit_consistent`` (Larsen–Morel–Miller 1987), and the PAIR's
    validity is their conjunction
    (:func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`).
    Opt-in (``False`` default): a redistribution closure is NOT assumed
    :math:`\beta`-consistent until it declares so with a citation.  ``True`` for
    Morel–Montry (BMC Eq. (42) sets :math:`\beta = 0`) and — vacuously — for the
    Cartesian identity closure (no redistribution term ⇒ all :math:`\alpha
    \equiv 0` ⇒ :math:`\beta \equiv 0`).  Read-only class attribute."""

    @classmethod
    def _registry_base(cls) -> type:
        return PoleAngularClosureBase

    @abstractmethod
    def __init__(
        self,
        angular: "AngularRedistribution",
        pairing: np.ndarray,
    ) -> None:
        r"""Family construction contract — a closure binds its two TENSOR FACTORS.

        The curvilinear redistribution operator factors as
        :math:`\mathcal{R} = R_{\rm spatial} \otimes A_{\rm angular}
        (\tau, \alpha, w)`, and a member of this family IS the angular
        factor.  So it takes exactly the two things that product is made
        of, and nothing else:

        * ``angular`` — the member-INDEPENDENT angular data (the
          :math:`\alpha`-dome, the starting direction, the measure), shared
          by every member;
        * ``pairing`` — the SPATIAL factor, ``(nx, n_mom, n_thread)``.

        What each member adds on top is its own :math:`\tau` and the
        derived :math:`c_{\rm in}` / :math:`c_{\rm out}` — which is
        precisely why those are NOT on the shared angular object: a
        second member (plain diamond, an angular-LD device) has to be
        able to choose differently.

        Every concrete strategy is constructible as ``cls(angular, pairing)``;
        the SNMesh default-closure dispatch
        (``default_angular_closure_class(coord)(angular, pairing)``)
        instantiates through this signature.  Concretes may WIDEN it with
        extra keyword slots — the two-positional call must stay valid.
        Abstract: declares the signature only; concrete ``__init__``
        bodies do not chain here.

        (Before the 2026-08-26 un-weld this read ``cls(sn_mesh)`` and every
        member reached into the mesh for its operands.  Nothing it reached
        for was a mesh fact: measured, the whole set was
        ``(quad, coord, ΔA)`` plus values derivable from them.)
        """
    # The M-M weighted-diamond closure derives two algebraic constants per
    # μ-level from its α-dome and τ weight (the index convention consumers
    # depend on — keep at point of use):
    #
    #   c_out[m] = α_{m+1/2} / τ_m
    #   c_in[m]  = (1−τ_m)/τ_m · α_{m+1/2} + α_{m−1/2}
    #
    # The closure is their CANONICAL owner (single source of truth): it
    # computes them once at construction and consumers read the gathered
    # ``(N,)`` views below; Cartesian IdentityAngularClosure returns neutral
    # zeros (α ≡ 0, τ ≡ 1 ⇒ c = 0) so slab consumers stay geometry-blind.
    # Ownership rationale: curvilinear_one_group.rst §sn-closure-c-constants-owned.

    def _gather_per_ordinate(
        self, per_level: "tuple[np.ndarray, ...]"
    ) -> np.ndarray:
        """Scatter per-:math:`\\mu`-level ``(M_p,)`` values to ``(N,)`` global order.

        Pure permutation (no arithmetic) keyed on :attr:`level_indices` — the
        per-level→global ordinate map the closure already owns.  Sphere is the
        single-level (``M_p = N``) case; cylinder has one block per level.
        """
        level_indices = self.level_indices
        N = sum(int(lvl.size) for lvl in level_indices)
        # ``np.zeros`` (not ``np.empty``): the level partition fully covers
        # [0, N) today so every slot is overwritten, but a future non-covering
        # partition then yields a deterministic 0 rather than uninitialised
        # garbage (elegance review, B1).
        out = np.zeros(N, dtype=np.float64)
        for level, values in zip(level_indices, per_level, strict=True):
            out[level] = values
        return out

    def _build_per_ordinate_cache(self) -> None:
        """Gather the three per-level constants (c_in / c_out / τ) to ``(N,)``
        ONCE at construction.

        The per-level→global gather is a pure permutation of immutable data,
        so caching it makes the public accessors O(1).  The cached arrays are
        marked READ-ONLY (``setflags(write=False)``) so a consumer holding a
        reference to the shared ``(N,)`` view (e.g. the ``StreamingCoefficientCache``
        populator) cannot corrupt the cache.

        Precondition: called as the LAST ``__init__`` step, after the
        per-level constants (``_c_*_per_level`` / ``_tau_per_level``) and
        ``level_indices`` are bound.
        """
        c_in_cache = self._gather_per_ordinate(self._c_in_per_level)
        c_out_cache = self._gather_per_ordinate(self._c_out_per_level)
        tau_cache = self._gather_per_ordinate(self._tau_per_level)
        c_in_cache.setflags(write=False)
        c_out_cache.setflags(write=False)
        tau_cache.setflags(write=False)
        self._c_in_per_ordinate_cache = c_in_cache
        self._c_out_per_ordinate_cache = c_out_cache
        self._tau_per_ordinate_cache = tau_cache

    #: Set by every concrete ``__init__`` (the family's two-tensor-factor
    #: contract); declared here so the base's :attr:`angular` accessor is
    #: typed against the contract rather than against one member.
    _angular: "AngularRedistribution"

    @property
    def angular(self) -> "AngularRedistribution":
        """The angular factor these coefficients were derived from.

        Provenance, and the diagnostic affordance R19 asks of a bound
        object: given a closure, recover the dome, the starting direction
        and the measure that produced its τ / c_in / c_out — without
        having to find the mesh that built it.
        """
        return self._angular

    @property
    def c_in_per_ordinate(self) -> np.ndarray:
        r"""Upstream-numerator closure constant :math:`c_{\rm in}` per global ordinate.

        ``(N,)`` array; :math:`(1-\tau_m)/\tau_m\,\alpha_{m+1/2} +
        \alpha_{m-1/2}` for M-M, all-zero for the Cartesian identity closure.
        Returns the read-only cache built once at construction (Fix 1, L16).
        """
        return self._c_in_per_ordinate_cache

    @property
    def c_out_per_ordinate(self) -> np.ndarray:
        r"""Denominator closure constant :math:`c_{\rm out}` per global ordinate.

        ``(N,)`` array; :math:`\alpha_{m+1/2}/\tau_m` for M-M, all-zero for the
        Cartesian identity closure.  Returns the read-only cache built once at
        construction (Fix 1, L16).
        """
        return self._c_out_per_ordinate_cache

    @property
    def tau_per_ordinate(self) -> np.ndarray:
        r"""Morel--Montry angular weight :math:`\tau` per global ordinate.

        ``(N,)`` array; the FUNDAMENTAL angular weight (Bailey--Morel--Chang
        2010 Eq. 43) the closure owns, from which
        :math:`c_{\rm out} = \alpha_{m+1/2}/\tau_m` and
        :math:`c_{\rm in} = (1-\tau_m)/\tau_m\,\alpha_{m+1/2} + \alpha_{m-1/2}`
        derive (:attr:`c_out_per_ordinate` / :attr:`c_in_per_ordinate`).
        ``1.0`` for every ordinate of the Cartesian identity closure (the
        neutral weight: the recurrence is then the identity).  Returns the
        read-only cache built once at construction.  The live sweep, matvec,
        and scan read THIS τ (via
        :attr:`~orpheus.transport.spatial.scheme.CellVisit.tau`), not the
        retired geometry-factory ``StreamingTerms.tau_mm`` — the τ ownership
        and Step-C retirement (with the Leg-1 structural-independence gate) are
        the record at
        ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-tau-step-c-closeout``.
        """
        return self._tau_per_ordinate_cache

    # ── Matvec strategy contract (the ABC's three abstract methods) ──
    #
    # The unified SN matvec (``loss_representation.py``) reads
    # ``sn_mesh.pole_angular_closure`` typed against THIS ABC and drives the
    # angular path through these three methods; declaring them abstract makes
    # the ABC the COMPLETE strategy contract.  Return types are deliberately
    # loose (``object`` for the per-level half-grid state) so the two
    # realizations vary: M-M returns a per-level ``_MMHalfGrid`` tuple from
    # :meth:`precompute_psi_state` (Identity returns ``None`` — no curvature
    # grid), while both honour the ``(denom, upstream_numer)`` / adjoint
    # array-pair shapes.

    @abstractmethod
    def precompute_psi_state(
        self,
        psi_view: np.ndarray,
        *,
        radial_characteristic: "RadialCharacteristicInteriorField | None" = None,
    ) -> object:
        r"""Precompute per-level closure state for one matvec pass.

        Returns opaque per-level state (M-M: a ``_MMHalfGrid`` tuple of
        seeded half-angle grids; Identity: ``None``) that
        :meth:`cell_contribution` consumes.  See the concrete overrides.

        ``radial_characteristic`` (#282 route (a), 2.5d): the composite's
        typed ψ½ block.  On a CARRYING level (R12a) the recurrence seed
        is read from its ``cells(level, -1)`` leg; ``None`` seeds those
        levels at zero (legitimate only for the ψ-independent
        COEFFICIENT use — the transpose walk's ``denom``-only state).
        Non-carrying levels inline the 2-point angular-edge
        extrapolation of ``psi_view`` regardless
        (:meth:`MorelMontryAngularSweep.edge_extrapolated_seed`).
        """
        ...

    @abstractmethod
    def cell_contribution(
        self,
        psi_state: Any,
        cell_idx: int,
        level_idx: int,
        within_positions: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Per-cell angular contribution to the cell-balance terms.

        Returns ``(denom_term, upstream_numer_term)`` — the closure's
        additions to the cell-balance denominator (``(ΔA/w)·c_out``) and
        upstream numerator (``(ΔA/w)·c_in·ψ_{m-1/2}``).  Zero arrays for
        the Cartesian identity closure.
        """
        ...

    @abstractmethod
    def angular_adjoint(
        self,
        numer_bar: "tuple[np.ndarray, ...]",
    ) -> "tuple[np.ndarray, dict[int, np.ndarray]]":
        r"""Reverse-mode adjoint of the matvec angular coupling (Wave O O.2b).

        Returns ``(psi_view_bar, seed_cells_bar)``: the ``(ng, N, nx)``
        bulk cotangent plus — #282 route (a) — one ``(ng, nx)``
        seed-cells cotangent per CARRYING level (keyed by level index),
        where the reverse recurrence STOPS: the seed is first-class
        state, so its cotangent lands on the output composite's
        ``radial_characteristic`` block instead of being scattered back
        onto the bulk (the retired strategy ``seed_adjoint``
        delegation).  Non-carrying levels scatter their (edge-
        extrapolation) seed cotangent onto the bulk internally.  Empty
        dict + zero array for the Cartesian identity closure.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# _MMHalfGrid — module-private typed accessor for the M-M half-angle grid
# ═══════════════════════════════════════════════════════════════════════
#
# The underscore prefix declares "module-private": consumers (matvec, sweep,
# tests) see only the public API of :class:`MorelMontryAngularSweep` and treat
# the half-grid as opaque strategy state; external code consumes via the
# :meth:`upstream` / :attr:`upstream_per_ordinate` / :attr:`trailing_face`
# accessors.
#
# ⛔ This comment claimed until 2026-08-13 that "the redistribution body inside
# the M-M class accesses the raw :attr:`faces` array directly".  That consumer
# does NOT exist and has not for some time: the paired ``(m, m+1)`` access it
# described was folded into the ``c_in`` / ``c_out`` constants (see the
# derivation at the head of this module), so :meth:`cell_contribution` reads
# :attr:`upstream_per_ordinate` ALONE.  Recorded rather than deleted because
# the claim was load-bearing — it was the stated justification for exposing the
# raw array at all, and a reader who believed it would have preserved that
# exposure on the next carve.


@dataclass(frozen=True, slots=True)
class _MMHalfGrid:
    r"""Typed accessor for the Morel-Montry half-angle face grid.

    Pattern 4 (illegal states unrepresentable) for the off-by-one trap the
    half-angle grid exposes.  The M-M recurrence produces :math:`M+1` face
    fluxes per level for :math:`M` ordinates:
    ``faces[g, 0, i] = ψ_{1/2, i, g}`` (Carlson seed, upstream of ordinate 0),
    ``faces[g, m, i] = ψ_{m-1/2, i, g}`` (upstream of ordinate m ≡ downstream
    of ordinate m-1), ``faces[g, M, i] = ψ_{M+1/2, i, g}`` (downstream of the
    last ordinate).

    Distinct consumers need DIFFERENT slices of this grid:

    * The **unified matvec** consumes the upstream-per-ordinate slice — one
      ``(ng, nx)`` block per ordinate for ``cell_balance_for_streaming``'s
      ``psi_angular_upstream`` argument. Use :meth:`upstream` (single
      ordinate) or :attr:`upstream_per_ordinate` (all ordinates).  This is
      the ONLY production consumer (:meth:`cell_contribution`).
    * The **endpoint defect** reads the far end alone. Use
      :attr:`trailing_face`.

    ⛔ A third bullet stood here until 2026-08-13: *"the redistribution fold
    — :math:`R_m = (\Delta A/w)/V (\alpha_{m+1/2}\phi_{m+1/2} -
    \alpha_{m-1/2}\phi_{m-1/2})` — uses the paired ``(m, m+1)`` access; use
    :attr:`faces` directly."*  **No such consumer exists.**  That pairing was
    folded into the ``c_in`` / ``c_out`` constants, so the redistribution
    never touches two faces at once and :meth:`cell_contribution` reads
    :attr:`upstream_per_ordinate` alone.  Kept as a record because it was the
    stated reason the raw array is exposed.

    The off-by-one trap (``faces[g, m, i]`` vs ``faces[g, m+1, i]``) is
    impossible by API design when consumers use :meth:`upstream` — the
    method's name AND signature enforce upstream-per-ordinate semantics. The
    raw :attr:`faces` array stays exposed so the redist body keeps its
    bit-identical paired access.

    Storage convention: ``faces`` shape ``(ng, M+1, nx)`` — group-leading,
    matching the M-M recurrence kernel.

    Attributes
    ----------
    faces :
        Shape ``(ng, M+1, nx)``. The raw half-angle grid produced by
        :func:`_psi_half_grid_single_level`. ``faces[g, 0, i]`` is the Carlson
        seed at ordinate 0 (= upstream of m=0); for ``m = 1, …, M``,
        ``faces[g, m, i]`` is :math:`\phi_{m-1/2, i, g}` (upstream of ordinate
        m, downstream of ordinate m-1); ``faces[g, M, i]`` is the downstream
        face of the last ordinate.
    """

    faces: np.ndarray

    @property
    def n_groups(self) -> int:
        """Number of energy groups (``faces.shape[0]``)."""
        return self.faces.shape[0]

    @property
    def n_ordinates(self) -> int:
        """Number of ordinates :math:`M` in this level (``faces.shape[1] − 1``)."""
        return self.faces.shape[1] - 1

    @property
    def n_cells(self) -> int:
        """Number of spatial cells (``faces.shape[2]``)."""
        return self.faces.shape[2]

    @property
    def upstream_per_ordinate(self) -> np.ndarray:
        r"""Shape ``(ng, M, nx)`` — upstream half-face flux for each ordinate.

        ``[g, m, i]`` is :math:`\phi_{m-1/2, i, g}` — the upstream face
        of ordinate ``m`` (equivalently, the downstream face of ordinate
        ``m-1``) in group ``g`` at cell ``i``. The matvec consumes this
        slice as ``psi_angular_upstream`` (one ``(ng, nx)`` block per
        ordinate). Excludes the trailing face ``faces[:, M, :]`` which
        is the downstream of the last ordinate (not consumed as anyone's
        upstream) — see :attr:`trailing_face`.
        """
        return self.faces[:, :-1, :]

    @property
    def trailing_face(self) -> np.ndarray:
        r"""Shape ``(ng, nx)`` — :math:`\psi_{M+1/2}`, the level's FAR
        angular endpoint.

        ``faces[:, M, :]``: the downstream face of the last ordinate, and
        the one slice of this grid that enters no cell balance.

        **Why nothing consumes it, and why that is correct.** The M-M
        closure is substituted INTO the balance — that substitution is
        precisely where the :math:`c_{\rm in}` / :math:`c_{\rm out}`
        constants come from — so a half-angle face appears only as some
        ordinate's UPSTREAM datum.  For the last ordinate the outgoing
        coefficient is :math:`c_{\rm out}[M-1] = \alpha_{M+1/2}/\tau = 0`,
        because the :math:`\alpha`-dome closes: a contract enforced at
        admission by
        :func:`~orpheus.geometry.reduced_operator._assert_alpha_dome_closes`
        and itself a CONSEQUENCE of the measure's antisymmetry, not an axiom
        of the (strictly one-sided) Lathrop--Carlson recursion.  So this face
        is computed and then annihilated.  :meth:`angular_adjoint` agrees
        independently: it seeds ``psi_half_bar[:, :M, :]`` only, leaving
        index ``M`` at zero, so the last ordinate has no path at all through
        the angular channel.

        **What it is FOR.** It is one of the level's TWO angular endpoints.
        The redistribution coefficient vanishes at both
        (:math:`\alpha_{1/2} = \alpha_{M+1/2} = 0`), so the balance decouples
        into a plain radial DD ODE at each, and production solves both — the
        near end becomes this grid's seed, the far end is marched directly
        into ``cells(p, +1)``.  Only one is imposed; their disagreement is
        the over-determination residual
        :meth:`MorelMontryAngularSweep.angular_endpoint_defect_per_level`.

        Reading this costs nothing: the recurrence already fills it
        (:func:`_psi_half_grid_single_level`'s loop writes ``m + 1`` for
        ``m = 0 … M-1``).
        """
        return self.faces[:, -1, :]

    @property
    def downstream_per_ordinate(self) -> np.ndarray:
        r"""Shape ``(ng, M, nx)`` — downstream half-face flux for each ordinate.

        ``[g, m, i]`` is :math:`\phi_{m+1/2, i, g}` — the downstream face
        of ordinate ``m`` (equivalently, the upstream face of ordinate
        ``m+1``). Excludes the leading face ``faces[:, 0, :]`` which is
        the Carlson seed (upstream of m=0, not any ordinate's downstream).
        """
        return self.faces[:, 1:, :]

    def upstream(self, m: int) -> np.ndarray:
        r"""Shape ``(ng, nx)`` — upstream half-face of ordinate ``m``.

        Returns :math:`\phi_{m-1/2, i, g}` per ``(g, i)``. The unified
        matvec consumes one of these per ordinate to populate
        ``cell_balance_for_streaming``'s ``psi_angular_upstream`` argument.
        """
        return self.faces[:, m, :]

    def downstream(self, m: int) -> np.ndarray:
        r"""Shape ``(ng, nx)`` — downstream half-face of ordinate ``m``.

        Returns :math:`\phi_{m+1/2, i, g}` per ``(g, i)``. Note that
        ``downstream(m) == upstream(m+1)`` for ``m < M-1``.
        """
        return self.faces[:, m + 1, :]


# ═══════════════════════════════════════════════════════════════════════
# Morel–Montry τ producer — the ANGULAR-scheme weight, owned by the closure
# ═══════════════════════════════════════════════════════════════════════
#
# τ — the Morel–Montry / weighted-diamond angular weight (BMC 2010 Eq. 43) —
# is a function of the quadrature ``(μ, w, levels)`` ALONE, an ANGULAR-scheme
# property the closure produces HERE (not read back from the streaming-geometry
# factory; Step C retired that twin).  Ownership + Step-C retirement:
# curvilinear_one_group.rst §sn-tau-step-c-closeout.
#
# STRUCTURAL INDEPENDENCE (vv-principles L11) — a constraint, not history: this
# is the closure's OWN replica of the factory arithmetic, NOT a call into
# ``contamination.morel_montry_weights``.  That function is the VERIFICATION
# reference for the Leg-1 cross-check; calling it in production would collapse
# the cross-check into a tautology (reference contamination).  Do NOT "tidy"
# the arithmetic below into a call to the reference.


@dataclass(frozen=True)
class MarchStart:
    r"""Structural facts about a μ-level's angular-march START edge (R12a).

    The two DISTINCT degeneracies that the raw first-ordinate
    Morel–Montry weight used to conflate into one float (campaign
    ruling T26 — ``τ_raw,0 = 0`` and ``τ_raw,0 = 1`` both read "not
    carrying", for unrelated structural reasons). Each fact is posed
    on the level's own realization, as a bit-exact identity — no
    tolerance, no derived arithmetic:

    Attributes
    ----------
    on_edge_node : bool
        The march-start edge direction IS an ordinate: the level's
        :math:`\eta`-minimum node lies on the singular set
        :math:`\Sigma = \{\xi = 0\}` (the :math:`\omega = \pi`
        most-inward direction; on the sphere, a node at
        :math:`\mu = -1`). The seed slot is then a rank-duplicate of
        :math:`\psi_0` — a bulk ordinate for free. [Cylinder
        NODE_ALIGNED even :math:`n_\varphi`; a folded NODE_ALIGNED
        rule; a future Gauss–Lobatto sphere.]
    degenerate : bool
        The :math:`\eta`-minimum is SHARED (:math:`\eta_0 = \eta_1`
        bit-exactly) — a double-cover tie: mirror partners on a full
        circle (products; odd :math:`n_\varphi` and STAGGERED
        included), hemisphere partners on level-symmetric rules. The
        midpoint edge collapses onto :math:`\eta_0`, so the seed's
        only consumption path — the recurrence's :math:`(1-\tau_0)`
        thread weight — vanishes. Dead state.

    The bit-equalities are honest *because the shipped rules are
    exact by construction*: roots-of-unity azimuths make mirror
    :math:`\eta`-ties and the :math:`\Sigma` node's :math:`\xi = 0`
    bit-exact (Q5.E/E3), and level-symmetric sign replication
    bit-copies :math:`\eta` across hemispheres. A rule built from
    sloppy trigonometry is classified as the realization it actually
    is — the cure for the pre-E3 5.6e-16 flip was making the rules
    exact, not tolerating the noise.
    """

    on_edge_node: bool
    degenerate: bool

    @property
    def consumes_independent_seed(self) -> bool:
        """Carrying ⟺ NEITHER degeneracy — a genuine off-node start.

        The R12a predicate: the M-M half-angle recurrence consumes an
        independent seed value iff the start edge is not itself an
        ordinate AND the thread weight that would carry the seed does
        not vanish. One name for the conjunction, two named conjuncts
        — never one float.
        """
        return not (self.on_edge_node or self.degenerate)


def march_start_structure_per_level(
    quad: Any,
    coord: CoordSystem,
) -> tuple[MarchStart, ...]:
    r"""The R12a facts per μ-level, read off the level's realization.

    The integer re-pose of the seed-presence predicate (T26): the
    consumer
    (:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_levels`)
    reads these two structural facts directly instead of deciding a
    structural question on the raw M-M float. The former encoding —
    ``τ_raw,0 ∈ (0,1)`` exclusive — is now a *theorem about the edge
    arithmetic*, gated bit-exactly per family
    (``on_edge_node ⟹ τ_raw,0 = 0``; ``degenerate ⟹ τ_raw,0 = 1``;
    neither ⟹ strict interior) rather than being the predicate itself.

    Level indexing matches :func:`morel_montry_tau_per_level`:
    the sphere is one level (the whole μ-ascending rule); cylinder
    levels index ``quad.level_indices``, each in stored η-ascending
    order.
    """
    if coord is CoordSystem.SPHERICAL:
        mu = quad.mu_x
        return (
            MarchStart(
                on_edge_node=bool(mu[0] == -1.0),
                degenerate=bool(len(mu) > 1 and mu[0] == mu[1]),
            ),
        )

    if coord is CoordSystem.CYLINDRICAL:
        eta = quad.mu_x
        xi = quad.mu_y
        out: list[MarchStart] = []
        for level_idx in quad.level_indices:
            first = level_idx[0]
            eta_level = eta[level_idx]
            out.append(
                MarchStart(
                    on_edge_node=bool(xi[first] == 0.0),
                    degenerate=bool(
                        len(level_idx) > 1 and eta_level[0] == eta_level[1]
                    ),
                )
            )
        return tuple(out)

    raise ValueError(
        f"march_start_structure_per_level supports SPHERICAL or "
        f"CYLINDRICAL coordinate systems; got {coord!r}. Cartesian has "
        f"no angular march and no seed."
    )


def non_carrying_levels(starts: tuple[MarchStart, ...]) -> tuple[int, ...]:
    r"""The level POSITIONS whose march start does not consume an
    independent seed.

    Returns the offenders — positions, not a bool (vv anti-#14: return
    the structure the docstring names) — so the caller reports WHICH
    level and WHY (R12a).  An empty tuple means every level is carrying.
    The MIXED case (some levels carrying, some not) is unreachable
    through every shipped factory but is first-class here: each level
    is judged on its own :class:`MarchStart` facts.
    """
    return tuple(
        p for p, s in enumerate(starts) if not s.consumes_independent_seed
    )


def assert_carrying_quadrature(quad: Any, coord: CoordSystem) -> None:
    r"""Raise unless EVERY mu-level of ``quad`` is carrying on ``coord``.

    The cylindrical ``SNMesh`` admission (Q5.6 step 6.3): the
    Morel--Montry azimuthal march is posed with an independent
    :math:`\psi_{1/2}` seed per level (route (a) — the
    :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
    march from the true :math:`q_{1/2}` source), so a rule whose march
    start is itself an ordinate, or whose seed thread weight vanishes,
    has no honest seed state to march and is refused at construction
    rather than silently mis-seeded.

    Admission is decided by STRUCTURE, not provenance: the facts come
    from :func:`march_start_structure_per_level` — the same producer
    :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_levels`
    reads — so a hand-built rule with carrying march starts is admitted
    with no ``folded_by`` tag, and a tagged quotient with a node on the
    mirror is refused.

    Raises
    ------
    ValueError
        Naming the FIRST offending level position and exactly the
        facts that are true on it.  The message deliberately does NOT
        contain the substring ``level structure`` — that fragment
        belongs to the earlier structure-less guard in
        :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`,
        which keeps ownership of slab/sphere cubatures.
    """
    starts = march_start_structure_per_level(quad, coord)
    offenders = non_carrying_levels(starts)
    if not offenders:
        return
    p = offenders[0]
    reasons = []
    if starts[p].on_edge_node:
        reasons.append(
            "starts on an ordinate (xi = 0 at the most-inward node, so "
            "the seed slot is a rank-duplicate of psi_0)"
        )
    if starts[p].degenerate:
        reasons.append(
            "has a shared eta-minimum (eta_0 == eta_1 bit-exactly, so "
            "the (1 - tau_0) thread weight vanishes and the seed is dead)"
        )
    raise ValueError(
        f"A cylindrical SNMesh admits only a quadrature whose every "
        f"mu-level is CARRYING (the R12a march-start predicate): level "
        f"{p} is non-carrying - it {' AND '.join(reasons)}. Use "
        f"Quadrature.folded_product(n_mu, n_phi) (the sigma_y quotient "
        f"of the staggered product) or any rule with the same "
        f"march-start structure. This rule HAS the required per-level "
        f"partitioning; the levels themselves are the problem."
    )


def _assert_tau_within_unit_interval(
    raw_levels: tuple[np.ndarray, ...],
) -> None:
    r"""The :math:`[0, 1]` membership guard on :math:`\tau_{\rm raw}` (Q5.5/T27).

    On a well-posed monotone march every ordinate lies inside its own
    angular cell, so :math:`\tau_{\rm raw} \in [0, 1]` by construction; a
    value outside (or a NaN) certifies an ill-posed march upstream —
    mis-ordered members (T22's ω-ordering), duplicated nodes, or a
    quadrature incompatible with the arm's edge convention (a raw 3-D
    rule on the 1-D spherical arm, #336).  Promoted from the
    :math:`[\tfrac12, 1]` absorber's silent absorption, which is how
    T22's mis-ordering laundered into a finite wrong answer instead of
    stopping here.
    """
    for p, tau_level in enumerate(raw_levels):
        inside = (tau_level >= 0.0) & (tau_level <= 1.0)
        if not np.all(inside):
            m = int(np.argmin(inside))
            raise ValueError(
                f"tau_raw[{m}] = {float(tau_level[m])!r} on level {p} lies "
                f"outside [0, 1]: the ordinate sits outside its own angular "
                f"cell, so the level is not a well-posed monotone march — "
                f"its members are mis-ordered (T22's ω-ordering), "
                f"duplicated, or the quadrature is not a μ-line rule "
                f"compatible with this closure's edge convention (a raw "
                f"3-D rule on the 1-D spherical arm reaches this; #336). "
                f"Fix the quadrature or the level ordering upstream; do "
                f"not absorb the value here."
            )


def angular_cell_edges_per_level(
    quad: Any,
    coord: CoordSystem,
) -> tuple[np.ndarray, ...]:
    r"""The ANGULAR CELL PARTITION of each μ-level, in radial cosine.

    One object, one producer.  A level's ordinates each own an angular
    *cell*; this returns that cell partition as ``(M+1,)`` edge arrays of
    the **radial** direction cosine (:math:`\mu` sphere, :math:`\eta`
    cylinder), so edge ``m`` is the boundary below ordinate ``m`` and
    edge ``M`` closes the level.

    Every coefficient that references "the boundary between cell
    :math:`m` and :math:`m+1`" must read it HERE — the partition is a
    *shared* structure on which the scheme then imposes **two different
    conditions**, and conflating those conditions is a live hazard:

    * :math:`\tau` imposes the **zeroth**-moment condition (P2, the
      barycentric coordinate — :func:`morel_montry_tau_per_level`);
    * :math:`\alpha` imposes the **first**-moment conservation recursion
      (P4, Hébert 3.397–3.399, after Alcouffe & O'Dell).

    ⚠ **Do NOT "unify" α and τ onto a single expression.** They share the
    partition, not the condition; forcing α to equal the geometric
    tangential cosine at these edges silently drives
    :math:`\delta \to 0`, i.e. :math:`\tau \to \tfrac12` — the angular
    *diamond* scheme, a different method (Hébert 3.406/3.431).

    Per-geometry conventions, and both are *derived*, not conventional:

    **Sphere** — cumulative-WEIGHT edges from :math:`\mu_{1/2} = -1`
    (:math:`\mu_{m+1/2} = \mu_{m-1/2} + w_m`, so a cell's width IS its
    quadrature weight).  This is Bailey–Morel–Chang 2010 Eq. (12)
    verbatim, corroborated independently by Lathrop 2000 p. 249
    (:math:`\sum \Delta\mu_m = 2`).  It works because a Gauss–Legendre
    weight *is* the cell's :math:`\mu`-measure.

    **Cylinder** — the **midpoint in** :math:`\omega`, with the level's
    endpoints at :math:`\omega = \pi` and :math:`\omega = 0` (hence
    :math:`\eta = \mp\sin\theta`).  The azimuthal march is a march in
    :math:`\omega`, arc by arc (T22b), so the cell boundary is the
    midpoint in :math:`\omega` and its radial cosine is
    :math:`\sin\theta\cos\omega_{m\pm1/2}`.  On an equispaced-:math:`\omega`
    rule this is exactly the half-angle boundary
    :math:`\omega_m \pm \Delta\omega/2`; taking the midpoint of the
    *stored* :math:`\omega` values keeps it correct for any monotone arc.

    ⛔ **Two conventions were REFUTED here, both measured (Q5.6.4,
    2026-08-11).**

    *The sphere's convention does not transplant.* Accumulating weights
    in :math:`\eta` (even correctly renormalised to
    :math:`\sum \bar{w} = 2\sin\theta`, BMC Eq. 52) violates P3 on this
    rule and **worsens with refinement**: ordinates outside their own
    cell go 0/4 → 4/8 → 12/16 → 28/32 at
    :math:`n_\varphi = 8/16/32/64`, and the solve diverges (NaN) from
    :math:`n_\varphi \ge 16`.  The reason is structural and is stated in
    ``derivations/discrete/sn/angular_differencing.py`` — *"weights are
    uniform in* :math:`\varphi` *-space, not* :math:`\eta` *-space"*: an
    arc cell's :math:`\eta`-measure is
    :math:`2\sin\theta\sin\omega_m\sin(\Delta\omega/2) \propto
    \sin\omega_m`, **not** constant, while the trapezoid weight is.  `[M]`
    the mismatch (cell :math:`\eta`-measure ÷ its mean, so a constant
    weight would need ``1.0`` everywhere) **WIDENS with refinement**:

    ==========  ==========================
    n_φ         ratio range over the level
    ==========  ==========================
    8           ``[0.5858, 1.4142]``
    16          ``[0.3045, 1.5307]``
    32          ``[0.1537, 1.5607]``
    64          ``[0.0770, 1.5683]``
    ==========  ==========================

    ⟹ BMC Eq. (52) is not a law; it is the statement that *in their*
    quadrature the weight equals the cell's :math:`\eta`-measure.  Ours
    does not — and increasingly does not — so we satisfy the same
    *predicate* by a different partition.  This is also the mechanism
    behind the worsening P3 violation above: the two disagree more at
    every order.

    *The η-midpoint (chord) partition, which this code used until
    2026-08-11, is this partition with its END CELLS STRETCHED.* Its
    interior edges are exactly :math:`\cos(\Delta\omega/2) \times` these
    (measured agreement 1e-16 — the identity is
    :math:`\tfrac12[\cos\omega_a + \cos\omega_b] =
    \cos(\tfrac{\omega_a+\omega_b}{2})\cos(\tfrac{\Delta\omega}{2})`,
    the :math:`\kappa` prefactor's sibling), while the endpoints stay
    pinned at :math:`\mp\sin\theta` **unscaled** — so the outermost cells
    stretch to absorb the shrink.  The :math:`\eta` error vanishes as
    :math:`\Delta\omega \to 0`, but the implied :math:`\omega`-width
    spread does **NOT**: it converges to :math:`\approx 17.45\,\%`
    (18.71 / 17.59 / 17.48 / 17.46 at
    :math:`n_\varphi = 8/16/32/64`) against a quadrature whose own cells
    are bit-exactly equal.  That O(1) inconsistency is what the retired
    :math:`[\tfrac12, 1]` absorber was compensating for.

    Returns
    -------
    tuple of ``(M_p + 1,)`` arrays, one per μ-level, ascending in the
    radial cosine.  Sphere is a single level; cylinder is one array per
    azimuthal level, indexed as ``quad.level_indices``.
    """
    if coord is CoordSystem.SPHERICAL:
        w = quad.weights
        N = quad.N
        mu_edge = np.zeros(N + 1)
        mu_edge[0] = -1.0
        for n in range(N):
            mu_edge[n + 1] = mu_edge[n] + w[n]
        return (mu_edge,)

    if coord is CoordSystem.CYLINDRICAL:
        mu_z = quad.mu_z
        edges: list[np.ndarray] = []
        for level_idx in quad.level_indices:
            eta = quad.mu_x[level_idx]
            xi = quad.mu_y[level_idx]
            M = len(level_idx)
            sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)
            # Stored order is η-ascending, i.e. ω-DESCENDING from near π
            # (most inward) to near 0 (most outward).
            omega = np.arctan2(xi, eta)
            # An ω-midpoint partition is only MEANINGFUL on a monotone
            # arc inside (0, π) — one traversal of the half-circle.  A
            # full-circle level carries ω of BOTH signs (the σ_y double
            # cover), for which "the midpoint in ω" is not defined and a
            # silent wrong answer is the alternative.  Refuse it HERE,
            # naming the cause, instead of letting the P3 guard report an
            # arbitrary out-of-range τ several frames later.
            if M > 1 and not np.all(np.diff(omega) < 0.0):
                raise ValueError(
                    f"angular_cell_edges_per_level: level {len(edges)} is "
                    f"not a monotone arc in omega (omega = {omega}), so its "
                    f"angular cells are not defined. The azimuthal march "
                    f"runs arc by arc, and a FULL-CIRCLE level covers each "
                    f"arc twice (omega of both signs — the sigma_y double "
                    f"cover). Use Quadrature.folded_product(n_mu, n_phi), "
                    f"or any rule whose levels are monotone half-circle "
                    f"arcs. A cylindrical SNMesh already refuses this rule "
                    f"at admission (assert_carrying_quadrature)."
                )
            edge_omega = np.empty(M + 1)
            edge_omega[0] = np.pi
            edge_omega[1:M] = 0.5 * (omega[:-1] + omega[1:])
            edge_omega[M] = 0.0
            edges.append(sin_theta * np.cos(edge_omega))
        return tuple(edges)

    raise ValueError(
        f"angular_cell_edges_per_level supports SPHERICAL or CYLINDRICAL "
        f"coordinate systems; got {coord!r}. Cartesian has no angular "
        f"march and hence no angular cell partition."
    )


def morel_montry_tau_per_level(
    quad: Any,
    coord: CoordSystem,
) -> tuple[np.ndarray, ...]:
    r"""The Morel–Montry angular closure weight :math:`\tau` per μ-level.

    ⚠ **``tau`` is overloaded three ways in this codebase.** THIS one is
    the **angular closure weight**: a dimensionless barycentric
    coordinate in :math:`[0, 1]`, one per ordinate. It is NOT the optical
    depth :math:`\Sigma_t s` (``peierls_nystrom``, MoC — also spelled
    ``tau``), and NOT the critical half-thickness in mean free paths
    (``fn_method``). See :mod:`orpheus.derivations.discrete.sn.angular_differencing`
    for the full nomenclature table.

    :math:`\tau_m = (\mu_m - \mu_{m-1/2})/(\mu_{m+1/2} - \mu_{m-1/2})`
    (Bailey–Morel–Chang 2010 Eq. 43 = Lathrop 2000 Eq. 23) — the
    predicate **P2**: :math:`\tau_m` is the *barycentric coordinate of
    the ordinate between the two edges of its own angular cell*,
    equivalently the UNIQUE closure weight exact for a flux **affine in
    the radial cosine**.  The cell edges come from the single partition
    producer :func:`angular_cell_edges_per_level`; P2 itself carries no
    geometry, so ONE body serves both arms.

    ⭐ **This is the whole τ — there is no "raw" and no "clamped".** The
    cylinder's :math:`[\tfrac12, 1]` absorber RETIRED at Q5.6.4
    (2026-08-11) and ``morel_montry_tau_raw_per_level`` retired with it
    (the distinction it named no longer exists); the ``[0, 1]``
    membership half became the P3 guard below.  Until Q5.4 this value
    was also the R12a seed-presence predicate's input (``τ_raw,0 ∈
    (0, 1)`` exclusive); the predicate is now posed on the structural
    facts directly (:func:`march_start_structure_per_level`), and the
    first-ordinate trichotomy — ``0`` where the start edge is a node,
    ``1`` where the start is η-degenerate, strict interior where the
    level genuinely consumes a seed — is a bit-exact gated CONSEQUENCE
    of those facts, not their definition. The full table and its
    rationale:
    ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-r12a``.

    Parameters
    ----------
    quad :
        Quadrature exposing ``mu_x`` (radial cosine — :math:`\mu` sphere,
        :math:`\eta` cylinder), ``weights``, ``mu_y``/``mu_z`` and
        ``level_indices`` (cylinder only) — the angular factor's measure.
    coord :
        :class:`~orpheus.geometry.CoordSystem` — ``SPHERICAL`` or
        ``CYLINDRICAL``.  Selects only the cell PARTITION (delegated);
        the τ formula is geometry-free.

    Returns
    -------
    tuple of ``(M_p,)`` arrays, one per μ-level.  Sphere is a single
    level ``(N,)``; cylinder is one ``(M_p,)`` array per azimuthal level,
    indexed as ``quad.level_indices``.

    Raises
    ------
    ValueError
        If any :math:`\tau_{\rm raw} \notin [0, 1]` — the membership
        guard (Q5.5/T27).  On a well-posed monotone march every
        ordinate lies inside its own angular cell, so an out-of-range
        value certifies an ILL-POSED march: mis-ordered members (T22's
        ω-ordering produced :math:`\tau_{\rm raw} = 1.079` here, then a
        NaN 400 lines downstream), or a quadrature incompatible with
        the arm's edge convention (a raw 3-D ``level_symmetric`` rule
        on the 1-D spherical arm measured
        :math:`\tau_{\rm raw} \in [-20.3,\, 1.13]`, 23 of 24 outside —
        consumed SILENTLY by the unclamped sphere closure until this
        guard landed; #336 tracks the refuse-or-reduce design).  The
        closed endpoints are legal: ``0`` is an edge-node start and
        ``1`` an η-degenerate tie
        (:func:`march_start_structure_per_level`).  ⚠ The guard does
        NOT catch the double cover — a full-circle level's
        ``[0, 1, 0, 1, …]`` fingerprint is entirely INSIDE
        :math:`[0, 1]`; that detector is the singular set
        :math:`\Sigma` / the fold criterion (T24), not this.
    """
    if coord not in (CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL):
        raise ValueError(
            f"morel_montry_tau_per_level supports SPHERICAL or CYLINDRICAL "
            f"coordinate systems; got {coord!r}. Cartesian uses the neutral "
            f"τ = 1.0 supplied by IdentityAngularClosure."
        )

    # ⭐ ONE generic body, both geometries (Q5.6.4, 2026-08-11).  P2 — the
    # barycentric coordinate of the ordinate between its own cell's two
    # edges — carries NO geometry; the geometry lives entirely in the cell
    # partition, which :func:`angular_cell_edges_per_level` owns.  Until
    # 2026-08-11 the sphere and cylinder each re-derived their own edges
    # inline here, which is how they came to disagree about what an
    # angular cell IS.
    cosines: tuple[np.ndarray, ...] = (
        (quad.mu_x,)
        if coord is CoordSystem.SPHERICAL
        else tuple(quad.mu_x[idx] for idx in quad.level_indices)
    )
    tau_per_level: list[np.ndarray] = []
    for mu, mu_edge in zip(
        cosines, angular_cell_edges_per_level(quad, coord)
    ):
        M = len(mu)
        tau = np.empty(M)
        for m in range(M):
            dmu = mu_edge[m + 1] - mu_edge[m]
            # The ½ fallback is a 0/0 regularization for a zero-width
            # cell (an η-degenerate tie), NOT a limiter on τ.
            tau[m] = (
                (mu[m] - mu_edge[m]) / dmu if abs(dmu) > 1e-15 else 0.5
            )
        tau_per_level.append(tau)
    levels = tuple(tau_per_level)
    _assert_tau_within_unit_interval(levels)
    return levels


# ⛔ RETIRED 2026-08-11 (Q5.6.4): ``morel_montry_tau_per_level`` used to be
# a WRAPPER that applied a cylinder-only ``max(0.5, min(1.0, τ))``
# absorber on top of a "raw" producer.  Both halves are gone and the two
# functions are now ONE (the name above is the survivor):
#
#   * the ``[0, 1]`` MEMBERSHIP half was already promoted at Q5.5 to the
#     raising guard ``_assert_tau_within_unit_interval`` — that guard IS
#     the predicate **P3** (an ordinate lies inside its own angular cell);
#   * the ``[½, 1]`` ABSORPTION half is RETIRED.  It never had a source:
#     `[M]` no reference prescribes any limiter on τ, the admissible range
#     is ``[0, 1]``, and BMC's own S₂ example gives
#     ``τ₁ = μ₁ + 1 = 1 − 1/√3 ≈ 0.4226 < ½`` (Eq. 47) — as does our own
#     sphere arm, where 4 of 8 τ sit below ½ at S₈ Gauss–Legendre.
#
#   ⭐ What it was actually compensating for: the η-midpoint cell
#   partition it sat on (see ``angular_cell_edges_per_level``), whose end
#   cells are stretched ~17.5 % in ω.  With the partition taken in ω the
#   compensation has nothing left to correct.  `[M]` and the absorber is
#   condemned on its own terms, with no MMS involved: the march implied
#   BY a clamped τ (ν-closure, BMC Eq. 43) OVERSHOOTS the level's own
#   endpoint by 1.6 % / 0.19 % / 0.024 % at n_φ = 8/16/32 where a derived
#   τ closes exactly — i.e. the clamped values correspond to no partition
#   of the level at all.
#
#   ⚠ Honest cost, ratified rather than hidden: on the anisotropic
#   cylindrical MMS at nx=320 the principled τ is BETTER at n_φ=8
#   (3.128e-3 vs 3.511e-3) and ~1.8–2× WORSE at n_φ = 16/32/64.
#   Principled ≠ more accurate; the scheme that satisfies P2/P3 wins over
#   one with a smaller number on a single manufactured fixture.


# ═══════════════════════════════════════════════════════════════════════
# The M-M recurrence kernel — pure algebra, module level
# ═══════════════════════════════════════════════════════════════════════
#
# The Morel--Montry weighted half-angle recurrence is pure algebra — all
# data (``ψ_level``, ``τ_level``, an optional seed) via arguments, no mesh
# state.  The mesh-bound strategy composes it (``_psi_half_grid_for_level``
# reads τ from ``self`` and delegates); algebraic-identity tests call
# :func:`compute_psi_half_per_level` with hand-built coefficient arrays — no
# closure instance (and hence no mesh) required.


def _psi_half_grid_single_level(
    psi_level: np.ndarray,
    tau_level: np.ndarray,
    psi_half_seed: np.ndarray | None = None,
) -> np.ndarray:
    r"""Pure-algebra M-M angular recurrence on a single level.

    Returns the half-angle grid :math:`\phi_{m\pm 1/2, i, g}` of
    shape ``(ng, M+1, nx)`` from cell-centres ``psi_level`` shape
    ``(ng, M, nx)`` and the angular closure weight ``tau_level`` shape
    ``(M,)``.
    ``psi_half[:, 0, :]`` is the recurrence seed (Phase D Carlson
    if supplied, else Phase B zero); subsequent slices are the
    downstream half-faces produced by the Morel--Montry weighted
    recurrence (BMC 2010 Eqs. (42)/(43); at :math:`\tau \equiv \tfrac12`
    this degenerates to Hébert's plain angular diamond, Eqs. 3.437/3.439
    sphere / 3.412/3.414 cylinder):
    ``ψ_{m+1/2,i,g} = (ψ_{m,i,g} - (1-τ_m)·ψ_{m-1/2,i,g}) / τ_m``.

    Pure kernel — accepts all data via arguments.  Used by the public
    :func:`compute_psi_half_per_level` exposure and by the mesh-bound
    :meth:`MorelMontryAngularSweep._psi_half_grid_for_level` (which
    reads ``tau_level`` from the strategy) via
    :meth:`MorelMontryAngularSweep.precompute_psi_state`.
    """
    ng, M, nx = psi_level.shape
    psi_half = np.empty((ng, M + 1, nx), dtype=psi_level.dtype)
    if psi_half_seed is None:
        psi_half[:, 0, :] = 0.0
    else:
        psi_half[:, 0, :] = psi_half_seed
    for m in range(M):
        tau_m = tau_level[m]
        psi_half[:, m + 1, :] = (
            psi_level[:, m, :] - (1.0 - tau_m) * psi_half[:, m, :]
        ) / tau_m
    return psi_half


def compute_psi_half_per_level(
    psi_level: np.ndarray,
    tau_level: np.ndarray,
    *,
    psi_half_seed: np.ndarray | None = None,
) -> _MMHalfGrid:
    r"""Return the half-angle grid :math:`\phi_{m\pm 1/2, i, g}`
    for one level under the M-M recurrence, wrapped in the private
    :class:`_MMHalfGrid` typed accessor.

    The hand-built-coefficient verification surface for the M-M
    recurrence — ``tau_level`` arrives via argument so
    algebraic-identity tests can verify the recurrence under
    hand-built :math:`\tau` (no closure instance, no mesh).
    Production code uses
    :meth:`MorelMontryAngularSweep.precompute_psi_state`, which reads
    :math:`\tau` from the mesh-bound strategy state and runs the SAME
    :func:`_psi_half_grid_single_level` kernel (Pattern 2 — single
    source of truth for the algebra).

    Parameters
    ----------
    psi_level :
        Shape ``(ng, M_p, nx)`` — the cell-centre angular flux at
        the :math:`M_p` ordinates of one level (sphere: every
        ordinate; cylinder: a per-:math:`\mu`-level azimuthal
        subset).
    tau_level :
        Shape ``(M_p,)``: the Morel-Montry angular closure weights
        :math:`\tau` for the level (there is no clamp — the
        :math:`[\tfrac12, 1]` absorber retired at Q5.6.4).
    psi_half_seed :
        The half-angle face flux seed VALUES :math:`\phi_{1/2,i,g}`,
        shape ``(ng, nx)``.  ``None`` seeds the recurrence at zero.
        Hand-built tests pass the array they mean; production seeds come
        from :meth:`MorelMontryAngularSweep.precompute_psi_state` (route (a)).

    Returns
    -------
    _MMHalfGrid
        Typed accessor wrapping the half-angle grid
        ``faces`` of shape ``(ng, M_p+1, nx)``.
    """
    faces = _psi_half_grid_single_level(
        psi_level, tau_level, psi_half_seed=psi_half_seed,
    )
    return _MMHalfGrid(faces=faces)


# ═══════════════════════════════════════════════════════════════════════
# MorelMontryAngularSweep — Phase B canonical (default)
# ═══════════════════════════════════════════════════════════════════════
#
# External code never reaches into the M-M algebra directly — every M-M
# consumer goes through the strategy's public surface:
# :meth:`precompute_psi_state` + :meth:`cell_contribution` (the live
# matvec/sweep path); hand-built-coefficient verification goes through
# the module-level :func:`compute_psi_half_per_level` (same kernel).


class MorelMontryAngularSweep(
    PoleAngularClosureBase, key="morel_montry_angular_sweep",
):  # noqa: E501  (#282 route (a) notes:)
    # #282 route (a): the M-M recurrence's half-angle seed ``ψ_{1/2,i,g}`` is
    # no longer a swappable strategy.  Seed dispatch (R12a): on a CARRYING
    # level (sphere) the seed is the composite's ψ½ STATE, read from the
    # ``radial_characteristic`` block; on the non-carrying cylinder levels the
    # 2-point angular-edge extrapolation of the input field is inlined
    # (:meth:`edge_extrapolated_seed`).  The trichotomy and the retired
    # ``PsiHalfAngleSeed`` zoo: curvilinear_one_group.rst §sn-direct-seed-r12a.
    r"""Canonical per-cell Morel--Montry weighted-DD angular recurrence.

    The weighted :math:`\tau` is Morel--Montry's (1984), in the
    Bailey--Morel--Chang 2010 Eqs. (42)/(43) form; Hébert supplies the cell
    balance and the sweep mechanics but **no** :math:`\tau` — see the
    module docstring's References and its warning.

    The Phase-B default for the curvilinear FD operator's angular
    redistribution.  Bound to an SNMesh at construction: all M-M coefficients
    (α-dome, ΔA/w, τ, c_in, c_out, level partition) are precomputed
    eagerly from the mesh's
    :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`, and
    the mesh-bound methods (:meth:`precompute_psi_state`,
    :meth:`cell_contribution`, :meth:`angular_adjoint`) read that state from
    ``self`` — callers never ship M-M data through arguments (Pattern 4).

    The half-angle recurrence
    :math:`\phi_{n+1/2,i,g} = (\phi_{n,i,g} - (1-\tau_n)\phi_{n-1/2,i,g})/\tau_n`
    (seeded at the Carlson direction — the #282 route-(a) direct seed, see the
    class-level note) and its redistribution fold
    :math:`R_{n,i,g} = (\Delta A/w)_{i,n}/V_i\,[\alpha_{n+1/2}\phi_{n+1/2,i,g}
    - \alpha_{n-1/2}\phi_{n-1/2,i,g}]` reduce to pure DD at :math:`\tau = 1/2`
    (the admissible range is :math:`[0, 1]` and BOTH arms run the derived
    weight unclamped; the cylinder's :math:`[\tfrac12, 1]` absorber retired
    at Q5.6.4).  The SAME recurrence runs
    inside :class:`~orpheus.transport.spatial.diamond.DiamondDifference`, so
    the apply matvec and the sweep solve the same discrete fixed point.

    ⚠ **That equivalence is PREVENTED, not detected** — say it that way,
    because the difference decides what a green suite is evidence of.  There
    is one recurrence body (:func:`_psi_half_grid_single_level`) and one
    coefficient pair (``c_in`` / ``c_out``, precomputed here and read through
    :func:`~orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
    by both consumers), so no input exists on which the two paths could
    disagree about the ANGULAR closure; there is no second definition to
    drift.  What still has teeth, and what to cite instead of an equivalence
    gate:

    * ``tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py`` — the
      scalar recurrence against its precomputed-split vectorized scan,
      bit-identical across the three consumer frames (two genuinely distinct
      float programs over the one algebra).
    * ``tests/sn/sweep/core/test_sweep_vs_apply_consistency.py`` — the
      end-to-end SI-vs-Krylov consistency leg: two different *solvers* over
      the one operator, which is where an apply-vs-sweep divergence would
      actually surface.
    * ``tests/sn/sweep/core/test_phase_c_gates.py`` — Gate 1.3 reciprocity
      and Gate 1.4 linearity on the composite.

    ⛔ This paragraph cited ``tests/sn/l1_analytical/test_pole_closure_sweep_equivalence.py``
    until 2026-08-13.  That file **has never existed** in this repository's
    history — a coverage claim with no gate behind it, which is worse than no
    claim because an audit trusts it.

    Derivation + apply↔sweep equivalence:
    ``docs/theory/methods/sn/curvilinear_one_group.rst §pole-mm-recurrence``
    and ``§sn-apply-sweep-equivalence``.  Cylinder loops the recurrence per
    :math:`\mu`-level (each with its own α-dome, ΔA/w, τ); sphere is the
    single-level (``M_p = N``) case.

    Parameters
    ----------
    sn_mesh : SNMesh
        The mesh + quadrature + materials bundle this strategy binds to
        (REQUIRED — the family's ``cls(sn_mesh)`` construction contract).
        M-M precomputes α-dome, ΔA/w, τ, c_in, c_out, level partition,
        μ_x, weights, Δr at construction; the strategy methods read these
        from ``self`` (no M-M data through arguments).  Tests that need the
        recurrence under hand-built coefficient arrays use the module-level
        :func:`compute_psi_half_per_level` (same kernel, all data via
        arguments).
    """

    is_linear: ClassVar[bool] = True
    """The M-M weighted DD angular recurrence is an affine combination
    of cell-centre values (constant α, ΔA/w, τ coefficients); the
    output is linear in ``psi_cells``."""

    beta_first_order_consistent: ClassVar[bool] = True
    r"""Morel–Montry is the UNIQUE weighted-diamond-in-angle closure that sets
    the Bailey–Morel–Chang first-order functional :math:`\beta = 0` (BMC 2010
    Eq. (42): :math:`\tau_m = (\mu_m - \mu_{m-1/2})/(\mu_{m+1/2} - \mu_{m-1/2})`,
    the weight that makes Eq. (41) vanish; BMC Table I shows the M-M sum is zero
    to round-off while step/diamond-in-angle are nonzero).  So an M-M curvilinear
    closure preserves the FIRST-order diffusion limit, not merely the
    leading-order one — the angular half of the pair-validity condition."""

    # ── M-M-only mesh-bound state (beyond the base's shared contract) ──
    # The α-dome / ΔA/w redistribution geometry per μ-level, the per-level
    # starting directions, and the R12a carrying-level set (which levels
    # own a first-class ψ½ STATE block — the seed-read dispatch of
    # ``precompute_psi_state`` / ``angular_adjoint``).  All bound eagerly
    # in ``__init__`` from the mesh's ``ReducedStreamingOperator``.
    _alpha_per_level: "tuple[np.ndarray, ...]"
    _dAw_per_level: "tuple[np.ndarray, ...]"
    _mu_start_per_level: "tuple[float, ...]"
    _mu_x: np.ndarray
    _carrying_levels: frozenset[int]

    def __init__(
        self,
        angular: "AngularRedistribution",
        pairing: np.ndarray,
    ) -> None:
        # The two TENSOR FACTORS, and nothing else (the un-weld arc's
        # Phase B): the angular factor carries the dome, the starting
        # direction and the measure; ``pairing`` is the spatial factor,
        # ``(nx, n_mom, n_thread)``.  All M-M coefficients are precomputed
        # here and the strategy methods read them from ``self``.
        _require_single_moment_pairing(pairing, type(self).__name__)
        self._gram = np.asarray(pairing, dtype=float)
        self._angular = angular
        # (Retained as PROVENANCE, and read through :attr:`angular` — the
        # coefficients below are derived from it, so a diagnostic that has
        # the closure can recover what they were built from.  It was
        # write-only when first written; an audit caught that, and a
        # stored-but-unread field is dead weight, not provenance.)

        coord = angular.coord
        quad = angular.quadrature
        N = int(np.asarray(quad.mu_x).size)
        # The single-moment contraction: today's per-cell spatial factor
        # is the one entry of a (1, 1) pairing — ΔA_i.
        delta_A = self._gram[:, 0, 0]

        # R12a (#282 route (a)): the carrying-level set — the levels whose
        # recurrence consumes independent starting-direction STATE.
        # Derived HERE from (quad, coord), which is where the facts live;
        # the mesh property that used to supply it reads the same producer
        # (``march_start_structure_per_level``), so this is the same value
        # from the same source, one hop shorter.
        # (The predicate: curvilinear_one_group.rst §sn-direct-seed-r12a.)
        self._carrying_levels = frozenset(
            p for p, s in enumerate(
                march_start_structure_per_level(quad, coord)
            )
            if s.consumes_independent_seed
        )

        # ── Per-level partition (M-M's concept, NOT the quadrature's)
        # Sphere: every ordinate is one level (M_p = N, n_levels = 1).
        # Cylinder: μ-levels from ProductQuadrature / LevelSymmetricSN.
        # The closure OWNS τ, produced HERE from the quadrature (see the
        # τ-producer note above for the structural-independence constraint).
        tau_per_level = morel_montry_tau_per_level(quad, coord)
        if coord is CoordSystem.SPHERICAL:
            # The angular factor is the shared AngularRedistribution
            # (one producer); the spatial factor is ΔA, and the closure
            # forms ΔA ⊗ 1/w itself rather than reading a stored product.
            self.level_indices = (np.arange(N),)
            self._alpha_per_level = angular.alpha_per_level
            self._dAw_per_level = (
                delta_A[:, None] / np.asarray(quad.weights)[None, :],
            )
            self._tau_per_level = tau_per_level
            self._mu_start_per_level = angular.mu_start_per_level
        elif coord is CoordSystem.CYLINDRICAL:
            # Same two factors, per μ-level (see the sphere arm).
            w = np.asarray(quad.weights)
            self.level_indices = tuple(
                np.asarray(lvl) for lvl in quad.level_indices
            )
            self._alpha_per_level = angular.alpha_per_level
            self._dAw_per_level = tuple(
                delta_A[:, None] / w[np.asarray(lvl)][None, :]
                for lvl in quad.level_indices
            )
            self._tau_per_level = tau_per_level
            self._mu_start_per_level = angular.mu_start_per_level
        else:
            raise ValueError(
                f"MorelMontryAngularSweep supports SPHERICAL or CYLINDRICAL "
                f"coordinate systems; got {coord!r}. Use "
                f"IdentityAngularClosure for Cartesian."
            )

        # ── Precompute c_in, c_out per level (algebraic constants) ──
        # c_out[m] = α_{m+1/2} / τ_m       — denominator coefficient
        # c_in[m]  = (1-τ_m)/τ_m · α_{m+1/2} + α_{m-1/2}
        #            — upstream-numerator coefficient
        c_in_levels: list[np.ndarray] = []
        c_out_levels: list[np.ndarray] = []
        for p in range(len(self.level_indices)):
            alpha = self._alpha_per_level[p]
            tau = self._tau_per_level[p]
            M_p = tau.size
            alpha_in = alpha[:M_p]
            alpha_out = alpha[1:M_p + 1]
            c_out_levels.append(alpha_out / tau)
            c_in_levels.append((1.0 - tau) / tau * alpha_out + alpha_in)
        self._c_in_per_level = tuple(c_in_levels)
        self._c_out_per_level = tuple(c_out_levels)
        # ── Cache the (N,) global-ordinate gather ONCE (Fix 1, L16) ──
        # Gathers c_in / c_out AND the owned τ (Issue #236 Phase 2 B3).
        self._build_per_ordinate_cache()

        # ── Angular-edge geometry (ψ-independent quadrature data) ──
        # μ_x feeds the per-level edge-extrapolation stencil (non-carrying
        # levels) and the level-slice views of the adjoint.
        self._mu_x = quad.mu_x

    # ── Recurrence kernels (PR-TYPED-6.5 Phase 2.3) ──────────────────
    # Mesh-bound instance methods. They read ``tau`` / ``alpha`` / ``dAw``
    # / ``V`` from ``self`` — no shipping of mesh data through
    # arguments. The only call-time inputs that vary across matvec
    # calls are ``psi_level`` (the angular-flux slice for one level)
    # and the seed VALUES ``psi_half_seed`` (carrier state on carrying
    # levels, the inlined edge extrapolation otherwise — #282 route (a)).

    def _psi_half_grid_for_level(
        self,
        psi_level: np.ndarray,           # (ng, M_p, nx) — ordinates within level p
        level_idx_p: int,
        *,
        psi_half_seed: np.ndarray | None = None,  # (ng, nx)
    ) -> np.ndarray:
        r"""Run the M-M angular recurrence on level ``level_idx_p``.

        Mesh-bound wrapper around :func:`_psi_half_grid_single_level`.
        Reads ``τ`` from ``self._tau_per_level[level_idx_p]`` — no
        coefficient data is shipped through arguments.

        Returns
        -------
        np.ndarray
            Half-angle grid, shape ``(ng, M_p+1, nx)``.
        """
        return _psi_half_grid_single_level(
            psi_level,
            self._tau_per_level[level_idx_p],
            psi_half_seed=psi_half_seed,
        )

    # ── Strategy Protocol surface (PR-TYPED-6.5 Phase 2.5+) ──────────
    # The matvec body consumes ``precompute_psi_state`` + ``cell_contribution``
    # to read M-M contributions per-cell without naming any M-M algebra
    # or coefficient (Pattern 1 — operator algebra reads as composition).

    def edge_extrapolated_seed(
        self,
        psi_level: np.ndarray,                      # (ng, M_p, nx)
        level_idx_p: int,
    ) -> np.ndarray:
        r"""The 2-point angular-edge extrapolation of a NON-carrying level.

        The recurrence seed :math:`\psi_{1/2,i}` is the field's value at the
        level's starting-direction edge :math:`\mu_{\rm start}`; on a level
        that carries NO independent ψ½ state (R12a: raw τ₀ ∈ {0, 1} — every
        production cylinder level) the operator-consistent seed is the input
        field extrapolated linearly in :math:`\mu` through the level's two
        most-inward distinct-μ ordinates:

        .. math::

           \psi_{1/2,i} \;=\; (1-t)\,\psi_{m_0,i} + t\,\psi_{m_1,i},
           \qquad
           t = \frac{\mu_{\rm start} - \mu_{m_0}}{\mu_{m_1} - \mu_{m_0}} .

        Exact on angle-flat and linear-in-μ fields, O(Δμ²)-consistent, linear
        in the input.  Bit-identical to the retired ``AngularEdgeExtrapolation``
        default on every NON-CARRYING cylinder (NODE_ALIGNED full-circle
        product rules hit its t = 0 degenerate; level-symmetric rules
        have a dead seed weight; a σ_y-folded cylinder never reaches
        this — its carrying levels march first-class ψ½ state, Q5.6);
        degenerate
        single-direction levels fall back to constant extrapolation (t = 0).
        The R12a trichotomy: curvilinear_one_group.rst §sn-direct-seed-r12a.
        """
        m0, m1, t = self._edge_seed_stencil(level_idx_p)
        return (1.0 - t) * psi_level[:, m0, :] + t * psi_level[:, m1, :]

    def _edge_seed_stencil(self, level_idx_p: int) -> tuple[int, int, float]:
        r"""The ``(m0, m1, t)`` stencil of :meth:`edge_extrapolated_seed`.

        Shared by the forward read and its adjoint scatter (Pattern 2 —
        the linear map and its transpose read ONE coefficient source).
        Degenerate single-direction levels return ``(m0, m0, 0.0)``.
        """
        mu = self._mu_x[np.asarray(self.level_indices[level_idx_p])]
        # ⚠ #361 — the LAST bare `argsort` on a level in `orpheus/`, and
        # it re-derives an ordering the producer already made: since the
        # 2026-08-13 fiber-key carve `level_indices` is eta-primary, so
        # `mu` arrives non-decreasing here. `[M]` `order[0] == 0` on
        # every shipped rule (so `m0` is right), but argsort is NOT the
        # identity on tied levels — it permutes within tied blocks, so
        # `m1` below is an arbitrary member of the next eta-group, and
        # m0/m1 index FLUX SLOTS. Not fixed with the carve because the
        # SPHERE arm passes `level_indices = (arange(N),)` (genuinely
        # unsorted), where the sort IS load-bearing; the fix needs the
        # trivial level's semantics decided with #336.
        #
        # ⛔ REACHABILITY IS NOW MEASURED (2026-08-26), and it refuted the
        # hypothesis this comment used to end on ("may be dead on the
        # cylinder arm, in which case the answer is retirement").  Dead on
        # the cylinder, yes — `[M]` 0 of 88 rules `assert_carrying_quadrature`
        # admits.  But that gate has ONE call site, inside
        # `case CoordSystem.CYLINDRICAL`; the SPHERE arm calls no admission
        # gate at all, and a Gauss-Lobatto polar rule builds a production
        # SNMesh(SPHERICAL) and reaches this line at 6 of 11 orders.
        # Retirement is off the table: it is the only seed path a
        # mu = -1-noded sphere rule has.
        #
        # ⭐ And where it IS reached the ambiguity is ANNIHILATED, by
        # theorem: `t` below is (mu_start - mu[m0]) / (mu[m1] - mu[m0]),
        # and a level is non-carrying precisely when `on_edge_node`, i.e.
        # mu_start == mu[m0] — so the numerator is exactly 0, `t` is
        # exactly 0.0, and the arbitrary `m1` is multiplied by zero.
        # `[M]` over 75 reachable non-carrying levels: tie-break live in
        # 0 of 75; perturbing slot m1 by 1e3 moves the seed by 0.000e+00
        # (control on m0: 1.0e+03).  So this is a latent-not-live defect.
        #
        # ⚠ IF THIS argsort IS EVER CHANGED, land a Gauss-Lobatto witness
        # row in the same step.  Without one the change ships with no gate
        # that can redden it (plan-authoring §6c).  Record + probes:
        # scratch/mu_start_reachability_census.md; #361 (CLOSED — a record
        # of a deliberate non-fix, not a work item) carries the analysis.
        order = np.argsort(mu)
        m0 = int(order[0])
        for cand in order[1:]:
            if abs(mu[int(cand)] - mu[m0]) > 1e-14:
                m1 = int(cand)
                t = float(
                    (self._mu_start_per_level[level_idx_p] - mu[m0])
                    / (mu[m1] - mu[m0])
                )
                return m0, m1, t
        return m0, m0, 0.0

    def precompute_psi_state(
        self,
        psi_view: np.ndarray,                       # (N, ng, nx)
        *,
        radial_characteristic: "RadialCharacteristicInteriorField | None" = None,
    ) -> "tuple[_MMHalfGrid, ...]":
        r"""Build per-level half-angle grids from the ψ½ state / edge seed.

        One ``_MMHalfGrid`` per level (sphere: 1 element; cylinder:
        ``n_levels`` elements).  Matvec consumes via :meth:`cell_contribution`.

        Seed dispatch (#282 route (a), R12a):

        * **carrying level** — the recurrence seed is the composite's
          FIRST-CLASS ψ½ state: ``radial_characteristic.cells(p, -1)``
          (the inward starting-direction leg).  The retired-strategy
          extrapolation of the ITERATE — the #282 back edge — is gone:
          the seed is upstream STATE in the augmented walk order.
          ``radial_characteristic=None`` seeds carrying levels at ZERO —
          legitimate ONLY for the ψ-independent coefficient use (the
          transpose walk's ``denom``-only state); value paths on a
          carrying mesh must hand the block in (the walk guards this).
        * **non-carrying level** — the inlined 2-point angular-edge
          extrapolation of ``psi_view`` (:meth:`edge_extrapolated_seed`)
          — bit-identical to the retired default (t = 0 exact on
          product rules; dead seed weight on level-symmetric rules).

        Parameters
        ----------
        psi_view :
            Current angular-flux iterate, canonical layout.
        radial_characteristic :
            The composite's typed ψ½ block (``None`` on non-carrying
            meshes and for coefficient-only state).
        """
        # (N, ng, nx) → (ng, N, nx) — reorder for level access.
        psi_g_first = psi_view.swapaxes(0, 1)
        per_level: list[_MMHalfGrid] = []
        for p, level_idx in enumerate(self.level_indices):
            psi_level = psi_g_first[:, level_idx, :]  # (ng, M_p, nx)
            if p in self._carrying_levels:
                if radial_characteristic is not None:
                    psi_half_seed_arr = radial_characteristic.cells(p, -1)
                else:
                    ng, _M, nx = psi_level.shape
                    psi_half_seed_arr = np.zeros((ng, nx))
            else:
                psi_half_seed_arr = self.edge_extrapolated_seed(psi_level, p)
            faces = self._psi_half_grid_for_level(
                psi_level, p, psi_half_seed=psi_half_seed_arr,
            )
            per_level.append(_MMHalfGrid(faces=faces))
        return tuple(per_level)

    def cell_contribution(
        self,
        psi_state: "tuple[_MMHalfGrid, ...]",
        cell_idx: int,
        level_idx: int,
        within_positions: np.ndarray,               # (n_mask,)
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Per-cell M-M contribution to the cell-balance ``(denom, upstream_numer)``.

        Reads ``ΔA/w``, ``c_in``, ``c_out`` and the level's
        precomputed half-angle grid from ``self`` / ``psi_state``.
        No M-M-specific arguments are passed in.

        Returns
        -------
        denom_term : np.ndarray
            Shape ``(n_mask,)``.  ``(ΔA/w) · c_out`` per ordinate-in-mask.
            Adds to the cell-balance denominator alongside streaming
            and collision terms.
        upstream_numer_term : np.ndarray
            Shape ``(ng, n_mask)``.  ``(ΔA/w) · c_in · ψ_{m-1/2,i,g}``
            per ``(group, ordinate-in-mask)``.  Adds to the cell-balance
            upstream-numerator alongside the spatial-streaming term.
        """
        p = level_idx
        dAw = self._dAw_per_level[p][cell_idx, within_positions]     # (n_mask,)
        c_in = self._c_in_per_level[p][within_positions]             # (n_mask,)
        c_out = self._c_out_per_level[p][within_positions]           # (n_mask,)
        denom_term = dAw * c_out                                     # (n_mask,)
        # psi_state[p].upstream_per_ordinate: (ng, M_p, nx) → (ng, n_mask).
        psi_ang = psi_state[p].upstream_per_ordinate[
            :, within_positions, cell_idx,
        ]
        upstream_numer_term = dAw[None, :] * c_in[None, :] * psi_ang  # (ng, n_mask)
        return denom_term, upstream_numer_term

    def angular_endpoint_defect_per_level(
        self,
        psi_state: "tuple[_MMHalfGrid, ...]",
        radial_characteristic: "RadialCharacteristicInteriorField",
    ) -> dict[int, np.ndarray]:
        r"""``D_p`` — the over-determination residual of level ``p``'s angular march.

        .. math::

           D_p \;:=\; \psi_{M+1/2}\big|_p \;-\; \psi^{\text{marched}}_p(+1)

        per carrying level, shape ``(ng, nx)`` each, keyed by level position.

        The level's angular march has **two** endpoints, not one.  The
        redistribution coefficient vanishes at each
        (:math:`\alpha_{1/2} = \alpha_{M+1/2} = 0`), so the cell balance
        decouples there into a plain radial DD ODE, and production solves
        **both** with one engine
        (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`,
        called twice from
        :meth:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve`):
        the inward leg becomes this grid's seed, the outward leg is stored as
        ``cells(p, +1)``.  The M-M recurrence, marched from the seed across
        all ``M`` ordinates, **also** predicts the far end — as
        :attr:`_MMHalfGrid.trailing_face`.  Only one of the two is imposed;
        ``D`` is what the unimposed one is violated by.

        Two structurally different discretizations of ONE point of phase
        space, previously computed and never compared.

        ⛔ **``D`` MUST NOT be used to correct the seed.**  The map's linear
        part is exactly :math:`(-1)^M I` (it follows from
        :math:`\prod_m (1-\tau_m)/\tau_m = 1`, gated on both arms in
        ``tests/sn/sweep/curvilinear/test_psi_half_positivity.py``), and BOTH
        endpoint values come from physics — so imposing both is an
        **over-determination**, a constraint on the interior solution, not an
        equation for a free parameter.  Zeroing ``D`` would merely force the
        marched endpoint onto the directly-marched one with no evidence the
        latter is the better of the two.

        ⛔ **``D`` IS NOT AN ERROR ESTIMATOR, AND MAY NOT VOTE ON**
        :math:`\tau`.  This is a measurement, not a caution.  `[M]` 2026-08-12
        against the **analytic** anisotropic cylindrical MMS
        (``build_cylindrical_anisotropic_mms_case``; nx = 80,
        ``inner_tol = 1e-13``, all 12 solves converged) the Pearson
        correlation of :math:`\log D` with :math:`\log` (true MMS error)
        across four :math:`\tau` variants runs ``+0.7515 / +0.2608 /
        +0.0630`` at :math:`n_\varphi` = 8 / 16 / 32, with 2/4 → 0/4 → 0/4
        rank agreement — i.e. it **degrades monotonically to zero** as angle
        refines.  Structurally it must: :math:`D = e_1 - e_2`, a DIFFERENCE
        of two truncation errors, hence small whenever both are large and
        equal.  ``D`` ranks the shipped partition first by 2.6–45× while
        being uncorrelated with accuracy; that ranking is therefore **not**
        evidence for the partition.  (Record:
        ``scratch/q65_endpoint_defect_findings.md`` §F7/§F10.)

        ⟹ What it legitimately is: a cheap, reference-free, pointwise
        **consistency** residual.  Its convergence under angular refinement
        is a property of the scheme worth pinning, and it is pinned — see
        ``tests/sn/sweep/curvilinear/test_angular_endpoint_defect.py``, which
        owns the measured ladder (do not copy those numbers here; that gate
        re-measures them).

        Parameters
        ----------
        psi_state :
            The per-level half-angle grids from :meth:`precompute_psi_state`.
        radial_characteristic :
            The composite's typed ψ½ block — the SAME object that seeded
            ``psi_state``.  Non-optional, unlike every sibling signature in
            this module: those accept ``None`` because they have a real
            fallback (a zero or edge-extrapolated seed), and this has none.
            Without the marched state there is no second endpoint, so the
            defect is **undefined**, not zero.  That is a typing contract,
            deliberately not a runtime ``raise`` — a "not ``None``" check
            here would be type-narrowing wearing a contract's clothes
            (``.claude/rules/coding-standards.md``), and the downstream
            ``AttributeError`` is immediate.

        Returns
        -------
        dict[int, np.ndarray]
            ``{level_position: D_p}`` over the **carrying** levels only.
            Empty on a mesh with no carrying level (no such level has two
            endpoints).  A non-carrying level seeds from
            :meth:`edge_extrapolated_seed`, which is an interpolation and not
            a solved endpoint, so no comparison is defined for it.
        """
        return {
            p: psi_state[p].trailing_face - radial_characteristic.cells(p, +1)
            for p in sorted(self._carrying_levels)
        }

    def angular_adjoint(
        self,
        numer_bar: "tuple[np.ndarray, ...]",
    ) -> "tuple[np.ndarray, dict[int, np.ndarray]]":
        r"""Adjoint of the matvec angular coupling (Wave O O.2b, #208).

        The Hilbert transpose ``Lᵀ`` of :meth:`StreamingOperator.apply`
        needs the reverse of the matvec's angular path, which the forward
        builds via :meth:`precompute_psi_state` + the per-cell
        :meth:`cell_contribution` ``angular_numer_upstream`` injection:

        .. math::

            \text{seed } \psi_{1/2} \;\to\;
            \text{M-M recurrence } \psi_{m-1/2} \;\to\;
            \text{angular\_numer} = (\Delta A/w)\,c_{\rm in}\,\psi_{m-1/2}.

        Given the cotangent of every ``angular_numer_upstream`` contribution
        (collected by the spatial reverse sweep, one ``(ng, M_p, nx)`` array
        per μ-level), reverse the level scatter and the recurrence down to
        the SEED cotangent, then route it by the R12a dispatch (#282
        route (a)):

        * **carrying level** — the seed is first-class ψ½ STATE, so the
          reverse recurrence STOPS here: the seed cotangent is returned
          in ``seed_cells_bar[p]`` and the caller lands it on the output
          composite's ``radial_characteristic`` block (the retired-strategy
          ``seed_adjoint`` delegation is gone with the zoo).
        * **non-carrying level** — the forward seed was the inlined
          edge extrapolation of the INPUT, so its cotangent scatters
          back onto the two stencil ordinates of the bulk
          (:meth:`_edge_seed_stencil` — the same coefficients as the
          forward, Pattern 2).  Bit-identical to the retired
          ``AngularEdgeExtrapolation.seed_adjoint``.

        The ``c_in``/``c_out`` roles do NOT change (the diagonal
        ``c_out`` lives in ``denom``, handled by the spatial diagonal).

        Parameters
        ----------
        numer_bar :
            One array per μ-level (matching :attr:`level_indices`), each shape
            ``(ng, M_p, nx)`` — the cotangent of that level's
            ``angular_numer_upstream``.

        Returns
        -------
        psi_view_bar : np.ndarray
            ``(ng, N, nx)`` g-first bulk cotangent.
        seed_cells_bar : dict[int, np.ndarray]
            Per CARRYING level (keyed by level index) the ``(ng, nx)``
            cotangent of the level's inward starting-direction cells —
            the walk adds it into the output seed block's ``cells(p, -1)``
            rows.  Empty on non-carrying meshes.
        """
        # Shapes from the cotangent arrays themselves (level partition
        # covers the quadrature; ``numer_bar`` matches ``level_indices``).
        ng = int(numer_bar[0].shape[0])
        nx = int(numer_bar[0].shape[2])
        N = int(self._mu_x.shape[0])
        psi_bar = np.zeros((ng, N, nx))
        seed_cells_bar: dict[int, np.ndarray] = {}
        for p, level_idx in enumerate(self.level_indices):
            level_idx = np.asarray(level_idx)
            M = level_idx.size
            # ── reverse the dAw·c_in level-scatter: upstream_bar = (ΔA/w)·c_in·numer_bar
            dAw_p = self._dAw_per_level[p]               # (nx, M)
            c_in_p = np.asarray(self._c_in_per_level[p])  # (M,)
            coeff = (dAw_p * c_in_p[None, :]).T           # (M, nx)
            upstream_bar = coeff[None, :, :] * numer_bar[p]   # (ng, M, nx)
            # ── reverse the M-M recurrence ψ_half[m+1] = (ψ[m] − (1−τ)ψ_half[m])/τ
            tau_p = np.asarray(self._tau_per_level[p])
            psi_half_bar = np.zeros((ng, M + 1, nx))
            psi_half_bar[:, :M, :] += upstream_bar        # upstream[:,m,:] = ψ_half[:,m,:]
            for m in range(M - 1, -1, -1):
                tau_m = tau_p[m]
                phb = psi_half_bar[:, m + 1, :]
                psi_bar[:, level_idx[m], :] += phb / tau_m
                psi_half_bar[:, m, :] += -((1.0 - tau_m) / tau_m) * phb
            seed_bar = psi_half_bar[:, 0, :]              # (ng, nx) seed cotangent
            # ── route the seed cotangent (R12a dispatch, #282 route (a)) ──
            if p in self._carrying_levels:
                # First-class ψ½ state: STOP — the caller lands this on
                # the output composite's radial_characteristic block.
                seed_cells_bar[p] = seed_bar
            else:
                # Inlined edge-extrapolation adjoint: scatter onto the
                # two stencil ordinates (the forward's coefficients).
                m0, m1, t = self._edge_seed_stencil(p)
                psi_bar[:, level_idx[m0], :] += (1.0 - t) * seed_bar
                psi_bar[:, level_idx[m1], :] += t * seed_bar
        return psi_bar, seed_cells_bar

    def __repr__(self) -> str:
        # The "()" repr is contractual: tests assert
        # ``repr(MorelMontryAngularSweep()) == "MorelMontryAngularSweep()"``.
        # For mesh-bound instances we keep the same shape; the mesh is
        # an implementation detail not surfaced in the repr.
        return "MorelMontryAngularSweep()"


# ═══════════════════════════════════════════════════════════════════════
# IdentityAngularClosure — slab default (no angular redistribution)
# ═══════════════════════════════════════════════════════════════════════
#
# Cartesian slab carries no angular redistribution term (no curvature → no
# Hébert §3.9.4 closure).  The Identity strategy makes the slab algebra a
# typed default (replacing a former ``pole_angular_closure is None`` matvec
# branch — a Pattern 4 leak): the strategy exists, has the same surface, and
# contributes zero to both the cell-balance denominator and upstream numerator.


class IdentityAngularClosure(PoleAngularClosureBase, key="identity_angular_closure"):
    r"""No-op pole-angular closure for Cartesian (slab + 2-D rectilinear).

    The Cartesian SN balance equation has no angular-redistribution term —
    Hébert §3.9.4's :math:`(\Delta A/w)` factor vanishes on flat geometry
    (cell faces are parallel, no curvature-driven coupling between consecutive
    ordinate sub-domains; see
    ``docs/theory/methods/sn/curvilinear_one_group.rst §balance-curvilinear``).
    This strategy returns ``(0, 0)`` from :meth:`cell_contribution`,
    contributing nothing to the matvec's per-cell denominator and upstream
    numerator, so the matvec body consumes the SAME
    ``cell_balance_for_streaming`` algebra for Cartesian as for sphere +
    cylinder — geometry-blind by data (Cardinal Rule 2), replacing a former
    ``pole_angular_closure is None`` matvec branch (a Pattern 4 leak).

    Parameters
    ----------
    sn_mesh : SNMesh
        Bound to the mesh so consumers have one uniform construction
        pattern.  Identity reads only the ordinate count off the angular
        factor's measure to size its zero-contribution returns; the
        ``level_indices`` attribute is the trivial single-level
        partition ``(arange(N),)``.
    """

    is_linear: ClassVar[bool] = True
    """Returning zero is the canonical linear operation."""

    beta_first_order_consistent: ClassVar[bool] = True
    r"""Vacuously :math:`\beta`-consistent: the Cartesian streaming operator
    :math:`\mu\partial_x` carries NO angular-redistribution term, so all
    Bailey–Morel–Chang :math:`\alpha` coefficients are identically zero (BMC 2010
    Eq. (41) is built entirely from the :math:`\alpha`'s; they arise only from
    the curvilinear :math:`(1-\mu^2)/r\,\partial_\mu` term — cf. BMC R–Z
    Eqs. (49)-(50)).  With no angular edge flux to close there is nothing to get
    wrong: :math:`\beta \equiv 0` term-by-term.  This is what makes the
    pair-validity predicate COLLAPSE to the spatial condition alone in Cartesian
    (:func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`)."""

    def __init__(
        self,
        angular: "AngularRedistribution",
        pairing: np.ndarray,
    ) -> None:
        # The same two tensor factors every member takes.  Cartesian's are
        # both NEUTRAL — a zero dome and a zero pairing — so this member needs
        # nothing from either beyond the ordinate count.
        _require_single_moment_pairing(pairing, type(self).__name__)
        self._angular = angular   # provenance; read through ``angular``
        self._N: int = int(np.asarray(angular.quadrature.mu_x).size)
        # Trivial single-level partition: every ordinate in one level.
        self.level_indices = (np.arange(self._N),)
        # Neutral M-M closure constants per level — Cartesian carries no
        # angular redistribution (all α ≡ 0 ⇒ c_in = c_out = 0), and the
        # neutral angular weight is ``τ ≡ 1`` (Issue #236 Phase 2 B3): the
        # recurrence ``(ψ̄ − (1−τ)ψ_in)/τ`` is then the identity.  Consumers
        # read these through the shared base contract
        # (``cell_contribution`` + the ``c_*``/``tau_per_ordinate``
        # accessors) — geometry-blind by data (Cardinal Rule 2).
        self._tau_per_level = (np.ones(self._N),)
        self._c_in_per_level = (np.zeros(self._N),)
        self._c_out_per_level = (np.zeros(self._N),)
        # ── Cache the (N,) global-ordinate gather ONCE (Fix 1, L16) ──
        # Cartesian always binds the neutral zeros; the cache is the read-only
        # ``(N,)`` zeros the slab visit-stamp reads (c_in == c_out == 0.0) plus
        # the neutral ``τ ≡ 1`` (Issue #236 Phase 2 B3) the slab visit reads as
        # ``CellVisit.tau`` — the identity M-M weight.
        self._build_per_ordinate_cache()

    # ── Strategy Protocol surface ────────────────────────────────────

    def precompute_psi_state(
        self,
        psi_view: np.ndarray,
        *,
        radial_characteristic: "RadialCharacteristicInteriorField | None" = None,
    ) -> None:
        """No state — Cartesian has no curvature half-grid to precompute.

        ``radial_characteristic`` is structurally ``None`` on Cartesian
        (R12a: no curvature ⇒ no starting-direction levels ⇒ the field
        cannot even be constructed on this mesh).
        """
        del psi_view, radial_characteristic
        return None

    def cell_contribution(
        self,
        psi_state: None,
        cell_idx: int,
        level_idx: int,
        within_positions: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Zero contribution to ``(denom, upstream_numer)``."""
        del psi_state, cell_idx, level_idx
        n = within_positions.size
        # ``(1, n)`` broadcasts over the group axis: "zero for every
        # group" is the honest object, and it is why this member does not
        # need to be told ``ng`` (an ENERGY fact) at construction.
        return np.zeros(n), np.zeros((1, n))

    def angular_adjoint(
        self,
        numer_bar: "tuple[np.ndarray, ...]",
    ) -> "tuple[np.ndarray, dict[int, np.ndarray]]":
        """Zero angular adjoint — Cartesian has no curvature coupling (O.2b).

        The seed-cotangent dict is empty by construction (no carrying
        levels on Cartesian, R12a).
        """
        ng, _, nx = numer_bar[0].shape
        return np.zeros((ng, self._N, nx)), {}

    def __repr__(self) -> str:
        return "IdentityAngularClosure()"


# ═══════════════════════════════════════════════════════════════════════
# default_angular_closure_class — coord-system → strategy factory dispatch
# ═══════════════════════════════════════════════════════════════════════


def default_angular_closure_class(coord: CoordSystem) -> "type[PoleAngularClosureBase]":
    """Return the default pole-angular-closure CLASS for a coordinate system.

    PR-TYPED-6.5 Phase 2.9.  The factory dispatch (instantiation with
    ``sn_mesh``) is the caller's job — typically ``SNMesh.__init__``
    after the ``match mesh.coord:`` block resolves geometry data.

    * ``CARTESIAN``    → :class:`IdentityAngularClosure`
    * ``SPHERICAL``    → :class:`MorelMontryAngularSweep`
    * ``CYLINDRICAL``  → :class:`MorelMontryAngularSweep`
    """
    if coord is CoordSystem.CARTESIAN:
        return IdentityAngularClosure
    if coord in (CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL):
        return MorelMontryAngularSweep
    raise ValueError(
        f"No default pole-angular-closure for coordinate system {coord!r}"
    )



__all__ = [
    "IdentityAngularClosure",
    "MorelMontryAngularSweep",
    "PoleAngularClosureBase",
    "angular_cell_edges_per_level",
    "assert_carrying_quadrature",
    "compute_psi_half_per_level",
    "default_angular_closure_class",
    "march_start_structure_per_level",
    "morel_montry_tau_per_level",
]
