r"""Angular-redistribution closure strategies for the curvilinear FD operator.

The curvilinear S\ :sub:`N` cell-balance (Hébert 2009 §3.9.4 Eq. 3.428)
carries an angular-redistribution term in **half-angle face fluxes**
:math:`\phi_{n\pm 1/2,i}` — the flux between consecutive ordinate
sub-domains.  They are NOT cell-centre values and cannot be computed
from cell-centres without a **closure**; this module supplies that
closure as a strategy family.

The production closure is the per-cell Morel--Montry weighted-DD angular
recurrence of Hébert Eqs. 3.437 / 3.439,
:math:`\phi_{n+1/2,i} = (\phi_{n,i} - (1-\tau_n)\phi_{n-1/2,i})/\tau_n`,
seeded at the Carlson starting direction :math:`\mu = -1` (where the
redistribution weight vanishes, :math:`\alpha_{1/2} = 0`).  It runs the
SAME algebra
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` runs inside
the sweep, lifted to operator level so the apply matvec and the sweep
solve the **same** discrete fixed point.  Since Issue #282 route (a) the
seed :math:`\phi_{1/2,i}` is first-class STATE marched directly from the
source (carrying levels — the sphere), or the inlined 2-point
angular-edge extrapolation on non-carrying cylinder levels; the retired
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
``docs/theory/methods/sn/curvilinear_one_group.rst §sn-pole-angular-closure-protocol``
(which also carries the Hébert citation correction: the primary source is
Hébert §3.9.4; Bailey--Morel--Chang 2010 is the auxiliary M-M-clamp
justification, not the wrong 2009 Bailey FE-diffusion paper the
pre-Phase-B geometry docstrings cited).

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Chapter 3 §3.9.4
  (pp. 141-144), Eqs. 3.418-3.439.  **Primary source** for the
  per-cell DD angular recurrence and the Carlson starting-direction
  initialisation.  Local copy: ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010). *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*. NSE
  165(2):149-169.  Auxiliary justification for the M-M clamp.
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
    * ``__init__`` constructible as ``cls(sn_mesh)`` — the family
      construction contract (abstract below): every closure binds to
      its :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` at
      construction (PR-TYPED-6.5 Phase 2.3), and the SNMesh
      default-closure dispatch instantiates through this signature
      (``default_angular_closure_class(coord)(mesh)``). Concretes may
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
    # closure binds to its SNMesh at construction (the ``cls(sn_mesh)``
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
    def __init__(self, sn_mesh: "SNMesh") -> None:
        """Family construction contract — a closure binds to its SNMesh.

        Every concrete strategy is constructible as ``cls(sn_mesh)``
        (PR-TYPED-6.5 Phase 2.3 mesh binding); the SNMesh
        default-closure dispatch
        (``default_angular_closure_class(coord)(mesh)``) instantiates
        through this signature. Concretes may WIDEN it with extra
        keyword slots — but the one-positional-mesh call must stay
        valid. Abstract: declares the signature only; concrete
        ``__init__`` bodies do not chain here.
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
        reference to the shared ``(N,)`` view (e.g. the ``GeometryCoefficients``
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
# the half-grid as opaque strategy state. The redistribution body inside the
# M-M class accesses the raw :attr:`faces` array directly; external code
# consumes via :meth:`upstream` / :attr:`upstream_per_ordinate` accessors.


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

    Two distinct consumers need DIFFERENT slices of this grid:

    * The **redistribution fold** — :math:`R_m = (\Delta A/w)/V \cdot
      (\alpha_{m+1/2} \phi_{m+1/2} - \alpha_{m-1/2} \phi_{m-1/2})` — uses the
      paired ``(m, m+1)`` access. Use :attr:`faces` directly.
    * The **unified matvec** consumes the upstream-per-ordinate slice — one
      ``(ng, nx)`` block per ordinate for ``cell_balance_for_streaming``'s
      ``psi_angular_upstream`` argument. Use :meth:`upstream` (single
      ordinate) or :attr:`upstream_per_ordinate` (all ordinates).

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
        upstream).
        """
        return self.faces[:, :-1, :]

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


def morel_montry_tau_raw_per_level(
    quad: Any,
    coord: CoordSystem,
) -> tuple[np.ndarray, ...]:
    r"""Produce the **UNCLAMPED** Morel–Montry raw weight :math:`\tau_{\rm raw}` per μ-level.

    :math:`\tau_{\rm raw,m} = (\mu_m - \mu_{m-1/2})/(\mu_{m+1/2} -
    \mu_{m-1/2})` — the raw Bailey–Morel–Chang Eq. 43 value BEFORE the
    cylinder's structural :math:`[\tfrac12, 1]` clamp.  Split out of
    :func:`morel_montry_tau_per_level` because the raw value carries
    structure the clamp destroys (it maps ``0 → ½``): it is the single
    source for BOTH the production τ (this, then the clamp) AND the **R12a
    seed-presence predicate**
    (:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_levels`):
    a level carries independent starting-direction state iff its
    first-ordinate ``τ_raw ∈ (0, 1)`` exclusive.  Bit-exact trichotomy:
    ``0`` on cylinder *product* rules, ``1`` on cylinder *level-symmetric*
    rules, ``∈ (0,1)`` on the sphere-GL dome — the full table and its
    rationale at
    ``docs/theory/methods/sn/curvilinear_one_group.rst §sn-direct-seed-r12a``.

    Parameters and per-geometry edge conventions are those of
    :func:`morel_montry_tau_per_level` (weight-sum edges from −1.0 for
    the sphere; η-midpoint edges with ±sinθ endpoints per level for the
    cylinder; the ½ degenerate fallback where an angular cell has zero
    width belongs to the RAW value — a 0/0 regularization, not the clamp).
    """
    if coord is CoordSystem.SPHERICAL:
        # Sphere τ (the former spherical_streaming producer, retired in
        # Step C): weight-sum edges from −1.0, unclamped τ_raw, ½ fallback.
        mu = quad.mu_x
        w = quad.weights
        N = quad.N
        mu_edge = np.zeros(N + 1)
        mu_edge[0] = -1.0
        for n in range(N):
            mu_edge[n + 1] = mu_edge[n] + w[n]
        tau_mm = np.empty(N)
        for n in range(N):
            dmu = mu_edge[n + 1] - mu_edge[n]
            tau_mm[n] = (
                (mu[n] - mu_edge[n]) / dmu if abs(dmu) > 1e-15 else 0.5
            )
        return (tau_mm,)

    if coord is CoordSystem.CYLINDRICAL:
        # Cylinder raw τ (the former cylindrical_streaming producer,
        # retired in Step C): η-midpoint edges with ±sinθ endpoints,
        # per μ-level.  The structural [½, 1] clamp is applied by
        # morel_montry_tau_per_level, NOT here.
        mu_z = quad.mu_z
        tau_per_level: list[np.ndarray] = []
        for level_idx in quad.level_indices:
            eta = quad.mu_x[level_idx]
            M = len(level_idx)
            sin_theta = np.sqrt(1.0 - mu_z[level_idx[0]] ** 2)
            eta_edge = np.zeros(M + 1)
            eta_edge[0] = -sin_theta
            for m in range(M - 1):
                eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])
            eta_edge[M] = sin_theta
            tau = np.empty(M)
            for m in range(M):
                deta = eta_edge[m + 1] - eta_edge[m]
                tau[m] = (
                    (eta[m] - eta_edge[m]) / deta if abs(deta) > 1e-15 else 0.5
                )
            tau_per_level.append(tau)
        return tuple(tau_per_level)

    raise ValueError(
        f"morel_montry_tau_raw_per_level supports SPHERICAL or CYLINDRICAL "
        f"coordinate systems; got {coord!r}. Cartesian uses the neutral "
        f"τ = 1.0 supplied by IdentityAngularClosure."
    )


def morel_montry_tau_per_level(
    quad: Any,
    coord: CoordSystem,
) -> tuple[np.ndarray, ...]:
    r"""Produce the Morel–Montry angular weight :math:`\tau` per μ-level.

    :math:`\tau_m = (\mu_m - \mu_{m-1/2})/(\mu_{m+1/2} - \mu_{m-1/2})`
    (Bailey–Morel–Chang 2010 Eq. 43) — the UNIQUE weight exact for a
    flux linear in :math:`\mu`.  This is a property of the angular
    quadrature only; it is produced HERE on the angular closure rather
    than on the streaming-geometry factory (Issue #236 Phase 2).

    The raw (unclamped) value comes from
    :func:`morel_montry_tau_raw_per_level` — the single source shared
    with the R12a seed-presence predicate (2.5d); this function applies
    the cylinder's structural clamp on top.

    Parameters
    ----------
    quad :
        Quadrature exposing ``mu_x`` (radial cosine / η), ``weights``,
        ``mu_z`` (axial cosine, cylinder only), and ``level_indices``
        (cylinder only).  This is ``sn_mesh.quad``.
    coord :
        :class:`~orpheus.geometry.CoordSystem` — ``SPHERICAL`` or
        ``CYLINDRICAL``.  The closure TYPE already dispatches on
        geometry; this branch only selects the edge convention.

    Returns
    -------
    tuple of ``(M_p,)`` arrays, one per μ-level.  Sphere is a single
    level ``(N,)``; cylinder is one ``(M_p,)`` array per azimuthal level.

    Notes
    -----
    The sphere weight is **UNCLAMPED** (W1; the sphere dome is
    non-singular on GL, :math:`\tau_\text{raw} \in \sim[0.39, 0.61]`).
    The cylinder weight is **CLAMPED** to :math:`[\tfrac12, 1]` — the
    most-inward azimuthal ordinate of a *product* rule sits exactly on
    the level boundary (:math:`\eta_0 = \eta_{1/2} = -\sin\theta`) so
    its raw weight is :math:`\tau_\text{raw} = 0` bit-exactly; the
    recurrence :math:`(\psi - (1-\tau)\psi)/\tau` would divide by zero
    unclamped (structural; #229).
    """
    raw = morel_montry_tau_raw_per_level(quad, coord)
    if coord is CoordSystem.SPHERICAL:
        return raw
    # Cylinder: the structural [½, 1] clamp, element-wise (bit-identical
    # to the pre-split inline ``max(0.5, min(1.0, tau_raw))``).
    clamped: list[np.ndarray] = []
    for tau_raw_level in raw:
        tau = np.empty_like(tau_raw_level)
        for m in range(tau_raw_level.size):
            tau[m] = max(0.5, min(1.0, float(tau_raw_level[m])))
        clamped.append(tau)
    return tuple(clamped)


# ═══════════════════════════════════════════════════════════════════════
# The M-M recurrence kernel — pure algebra, module level
# ═══════════════════════════════════════════════════════════════════════
#
# The Hébert Eqs. 3.437 / 3.439 half-angle recurrence is pure algebra — all
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
    ``(ng, M, nx)`` and τ clamp ``tau_level`` shape ``(M,)``.
    ``psi_half[:, 0, :]`` is the recurrence seed (Phase D Carlson
    if supplied, else Phase B zero); subsequent slices are the
    downstream half-faces produced by Hébert Eqs. 3.437 / 3.439:
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
        Shape ``(M_p,)``: Morel-Montry :math:`\tau` clamp values
        for the level.
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
    r"""Canonical Hébert §3.9.4 per-cell M-M weighted DD angular recurrence.

    The Phase-B default for the curvilinear FD operator's angular
    redistribution.  Bound to an SNMesh at construction: all M-M coefficients
    (α-dome, ΔA/w, τ clamp, c_in, c_out, level partition) are precomputed
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
    (the M-M clamp keeps :math:`\tau \in [1/2, 1]`).  The SAME recurrence runs
    inside :class:`~orpheus.transport.spatial.diamond.DiamondDifference`, so
    the apply matvec and the sweep solve the same discrete fixed point (pinned
    by :file:`tests/sn/l1_analytical/test_pole_closure_sweep_equivalence.py`;
    derivation + apply↔sweep equivalence:
    ``docs/theory/methods/sn/curvilinear_one_group.rst §pole-mm-recurrence``
    and ``§sn-apply-sweep-equivalence``).  Cylinder loops the recurrence per
    :math:`\mu`-level (each with its own α-dome, ΔA/w, τ clamp); sphere is the
    single-level (``M_p = N``) case.

    Parameters
    ----------
    sn_mesh : SNMesh
        The mesh + quadrature + materials bundle this strategy binds to
        (REQUIRED — the family's ``cls(sn_mesh)`` construction contract).
        M-M precomputes α-dome, ΔA/w, τ clamp, c_in, c_out, level partition,
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
        sn_mesh: "SNMesh",
    ) -> None:
        # The mesh binding is REQUIRED (the ``cls(sn_mesh)`` construction
        # contract): all M-M coefficients are precomputed here and the
        # strategy methods read them from ``self``.
        #
        # R12a (#282 route (a)): the carrying-level set — the levels whose
        # recurrence consumes independent starting-direction STATE.
        # Single-sourced from the mesh predicate (which reads the raw producer
        # ``morel_montry_tau_raw_per_level``); safe here because the predicate
        # needs only ``(quad, coord)``, both bound before the closure is built.
        # (The τ_raw ∈ (0,1) predicate: curvilinear_one_group.rst §sn-direct-seed-r12a.)
        self._carrying_levels = frozenset(sn_mesh.radial_characteristic_levels)

        coord = sn_mesh.coord
        quad = sn_mesh.quad
        reduced = sn_mesh.reduced
        N = quad.N

        # ── Per-level partition (M-M's concept, NOT the quadrature's)
        # Sphere: every ordinate is one level (M_p = N, n_levels = 1).
        # Cylinder: μ-levels from ProductQuadrature / LevelSymmetricSN.
        # The closure OWNS τ, produced HERE from the quadrature (see the
        # τ-producer note above for the structural-independence constraint).
        tau_per_level = morel_montry_tau_per_level(quad, coord)
        if coord is CoordSystem.SPHERICAL:
            # Factory contract: spherical_streaming populates the sphere
            # fields of ReducedStreamingOperator.
            assert reduced.mu_start is not None
            assert reduced.alpha_half is not None
            assert reduced.redist_dAw is not None
            self.level_indices = (np.arange(N),)
            self._alpha_per_level = (reduced.alpha_half,)
            self._dAw_per_level = (reduced.redist_dAw,)
            self._tau_per_level = tau_per_level
            self._mu_start_per_level = (float(reduced.mu_start),)
        elif coord is CoordSystem.CYLINDRICAL:
            # Factory contract: cylindrical_streaming populates the
            # per-level fields of ReducedStreamingOperator.
            assert reduced.mu_start_per_level is not None
            assert reduced.alpha_per_level is not None
            assert reduced.redist_dAw_per_level is not None
            self.level_indices = tuple(
                np.asarray(lvl) for lvl in quad.level_indices
            )
            self._alpha_per_level = tuple(reduced.alpha_per_level)
            self._dAw_per_level = tuple(reduced.redist_dAw_per_level)
            self._tau_per_level = tau_per_level
            self._mu_start_per_level = tuple(
                float(v) for v in reduced.mu_start_per_level
            )
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
        default on every production cylinder (product rules hit its t = 0
        degenerate; level-symmetric rules have a dead seed weight); degenerate
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
        psi_view: np.ndarray,                       # (N, ng, nx, 1) canonical
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
        pattern.  Identity reads only ``sn_mesh.quad.N`` and
        ``sn_mesh.ng`` to size its zero-contribution returns; the
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

    def __init__(self, sn_mesh: "SNMesh") -> None:
        self._sn_mesh = sn_mesh
        self._N: int = sn_mesh.quad.N
        self._ng: int = sn_mesh.ng
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
        return np.zeros(n), np.zeros((self._ng, n))

    def angular_adjoint(
        self,
        numer_bar: "tuple[np.ndarray, ...]",
    ) -> "tuple[np.ndarray, dict[int, np.ndarray]]":
        """Zero angular adjoint — Cartesian has no curvature coupling (O.2b).

        The seed-cotangent dict is empty by construction (no carrying
        levels on Cartesian, R12a).
        """
        nx = int(numer_bar[0].shape[2])
        return np.zeros((self._ng, self._N, nx)), {}

    def __repr__(self) -> str:
        return f"IdentityAngularClosure(sn_mesh=<{self._sn_mesh!r}>)"


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
    "compute_psi_half_per_level",
    "default_angular_closure_class",
    "morel_montry_tau_per_level",
    "morel_montry_tau_raw_per_level",
]
