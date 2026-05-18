r"""Angular-redistribution closure strategies for the curvilinear FD operator.

Why this abstraction exists
===========================

The curvilinear S\ :sub:`N` cell-balance equation (Hébert 2009 §3.9.4
Eq. 3.428) carries an angular redistribution term

.. math::
   :label: pole-redistribution-balance

   \frac{\Delta S_i}{2\,\mathcal{W}_n}
   \bigl[\alpha_{n+1/2}\,\phi_{n+1/2,i}
       - \alpha_{n-1/2}\,\phi_{n-1/2,i}\bigr]

evaluating :math:`\int \partial_\mu[(1-\mu^2)\psi]\,d\mu` per ordinate
sub-domain (Hébert Eq. 3.420).  The :math:`\phi_{n\pm 1/2,i}` are
**half-angle face fluxes** between consecutive ordinate sub-domains;
they are NOT cell-centre values and they are NOT computable from
cell-centres without a closure.

The pre-Phase-B :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
matvec evaluated the half-angle face fluxes via Morel--Montry
:math:`\tau`-weighted **symmetric interpolation**

.. math::
   :label: legacy-mm-symmetric-interpolation

   \phi_{n+1/2,i} \;\approx\; \tau_n \, \phi_{n+1,i}
                              + (1 - \tau_n) \, \phi_{n,i}

— exact when :math:`\psi` is constant in :math:`\mu` (the **flat-flux
collapse** the Bailey :math:`\Delta A/w` factor enforces by
construction), but only :math:`\mathcal{O}(1)` accurate on smooth
angularly-varying :math:`\psi`.  The factor-of-two truncation gap
identified in Issue #168 Defect 3 (memo §3, lines 198-244) is exactly
the gap between this collapsed evaluation and the canonical Hébert
form: at the sphere pole on an MMS that is **not** angle-flat, the
collapsed redistribution overcorrects by a factor of two.

Phase B (this module) ships the canonical fix: the **per-cell M-M
weighted DD angular recurrence** of Hébert Eqs. 3.437 / 3.439 with
the Morel--Montry :math:`\tau` clamp,

.. math::
   :label: dd-angular-recursion

   \phi_{n+1/2,i} \;=\;
   \frac{\phi_{n,i} \;-\; (1 - \tau_n)\,\phi_{n-1/2,i}}{\tau_n},
   \qquad \phi_{1/2,i} = 0,

initialised at the Carlson zero-weight starting direction
:math:`\mu = -1` (Hébert Eqs. 3.432-3.435 give the source-driven
sweep that fixes :math:`\phi_{1/2,i}` for the **inverse** problem;
for the **forward apply** matvec we adopt :math:`\phi_{1/2,i} = 0`,
the unique choice that makes the recursion's seed consistent with
:math:`\alpha_{1/2} = 0` and that the sweep converges to under
fixed-point iteration).  At :math:`\tau_n = 1/2` the recurrence
reduces to the pure DD form :math:`\phi_{n+1/2,i} = 2\,\phi_{n,i}
- \phi_{n-1/2,i}` (Hébert Eqs. 3.437 / 3.439).  For
:math:`\tau_n \in (1/2, 1]` the M-M clamp gives weighted-DD with
guaranteed positive M-M weighting (Bailey-Morel-Chang 2010).

The recursion runs once per (cell, group) across the GL-sorted
ordinates :math:`n = 1, \ldots, N` — the same algebra the
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` cell-update
runs inside the sweep, lifted up to operator level so the apply
matvec and the sweep solve the **same** discrete fixed point.

Architectural mirror of BoundaryFaceFlux (Phase A, retired in Phase C)
======================================================================

This module's Protocol shape was modelled after the Phase A
``BoundaryFaceFlux`` Protocol (which has since been **retired** by
Issue #168 Phase C — see commit ``3fd1302``).  ``PoleAngularClosure``
stays because the sphere centre is **intrinsic geometry** (a
coordinate-system singularity), not an external boundary; the
Phase C sweep-frame matvec subsumed ``BoundaryFaceFlux`` but the
pole-angular closure is a separate concern that remains
strategy-parametrised.  The Phase A architectural mirror at the
time of writing was:

* ``@runtime_checkable Protocol`` for the strategy contract;
* concrete ABC :class:`PoleAngularClosureBase` layered on
  :class:`~orpheus.numerics.registry.RegistryMixin` for self-registering
  strategies;
* class-level ``is_linear: bool`` trait;
* ``frozen=True, slots=True`` dataclass implementations with
  explicit ``__repr__``;
* legacy reproducer (:class:`BaileyFlatFluxRedist`) preserved for
  ablation studies / back-bisection.

The unification choice for cylindrical was: **one Protocol, optional
per-level loop**.  The strategy accepts an optional ``level_indices``
argument; when ``None``, the strategy treats ordinates as a single
sphere-style level (sphere geometry); when populated, the strategy
loops over the per-:math:`\mu`-level azimuthal sub-problems with each
level's own :math:`\alpha_{n\pm 1/2}` dome and :math:`\Delta A/w`
geometry factor.  This keeps the Protocol surface area small and lets
the same concrete strategy serve both geometries — the per-cell DD
angular recurrence is structurally identical in both cases, only the
ordinate index list changes.

Hébert citation correction (Issue #168 Phase B)
================================================

The pre-Phase-B :mod:`orpheus.geometry.reduced_operator` docstrings
cited "Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
*A piecewise linear finite element discretization of the diffusion
equation for arbitrary polyhedral grids*. JCP 227, 3738-3757."  This
is the **wrong Bailey paper** — it is a piecewise-linear FE diffusion
paper unrelated to curvilinear S\ :sub:`N` α-recursion.  The intended
reference is **Bailey, Morel & Chang (2010)**, NSE 165(2):149-169,
"Asymptotic Diffusion-Limit Accuracy of Sn Angular Differencing
Schemes" (LLNL preprint LLNL-JRNL-420356; OA at
https://www.osti.gov/servlets/purl/1020346).

Phase B corrects the citations across the operator and geometry
modules.  The math itself was always Hébert §3.9.4 — Bailey-Morel-Chang
2010 re-derives the M-M weighted-diamond clamp with formal-:math:`\varepsilon`
asymptotic-diffusion-limit analysis but does not discuss the
angular-redistribution closure.  Hébert is the canonical primary
source.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Chapter 3 §3.9.4
  (pp. 141-144), Eqs. 3.418-3.439.  **Primary source** for the
  per-cell DD angular recurrence and the Carlson starting-direction
  initialisation.  Local copy: ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010). *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*. NSE
  165(2):149-169.  Auxiliary justification for the M-M clamp.
* :doc:`/theory/discrete_ordinates` — "PoleAngularClosure (Issue
  #168 Phase B)".
* Issue #168 design memo — ``.claude/plans/issue_168_design.md``.
* Phase A closeout —
  ``.claude/agent-memory/method-implementer/issue_168_phase_a_closeout.md``.

See also
========

* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`` —
  **RETIRED Phase C** (Issue #168, commit ``3fd1302``).  The Phase A
  boundary-flux Protocol whose architecture this module mirrors; the
  sweep-frame matvec rewrite subsumed it.  Retained as a
  cross-reference for the architectural-mirror pattern only.
* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` —
  per-cell-update strategy contract; the curvilinear sweep also runs
  the DD angular recurrence inside its
  :class:`~orpheus.sn.spatial.diamond.DiamondDifference` strategy
  (with :math:`\tau = 1/2`).  Phase B aligns the apply matvec with
  the sweep's math at the operator level.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, ClassVar, Protocol, runtime_checkable

import numpy as np

from orpheus.geometry import CoordSystem
from orpheus.numerics.registry import RegistryMixin

from .psi_half_angle_seed import (
    CarlsonInwardSweep,
    CarlsonSweepContext,
    PsiHalfAngleSeed,
)

if TYPE_CHECKING:  # pragma: no cover
    from ..geometry import SNMesh


# ═══════════════════════════════════════════════════════════════════════
# PoleAngularClosure Protocol
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class PoleAngularClosure(Protocol):
    r"""Strategy contract for evaluating the curvilinear angular redistribution.

    The strategy consumes the cell-centre angular flux array
    :math:`\psi_{n,i,g}` plus the precomputed connection-coefficient
    arrays and returns the **per-(group, ordinate, cell)
    redistribution term** that the matvec adds to streaming +
    collision.

    The redistribution at cell :math:`i`, ordinate :math:`n`, group
    :math:`g` is

    .. math::

       R_{n,i,g} \;=\; \frac{(\Delta A/w)_{i,n}}{V_i}
                       \bigl[\alpha_{n+1/2}\,\phi_{n+1/2,i,g}
                           - \alpha_{n-1/2}\,\phi_{n-1/2,i,g}\bigr]

    where :math:`\phi_{n\pm 1/2,i,g}` are the **half-angle face
    fluxes** between consecutive ordinate sub-domains, evaluated by
    the strategy's closure choice.

    The unified Protocol covers spherical (``level_indices is None``)
    and cylindrical (``level_indices`` provided — one entry per
    :math:`\mu`-level, the strategy loops over levels).

    Parameters
    ----------
    psi_cells : np.ndarray
        Cell-centre angular flux, shape ``(ng, N, nx)``.
    alpha_half : np.ndarray or list[np.ndarray]
        :math:`\alpha`-dome values:

        * For ``level_indices is None`` (sphere): shape ``(N+1,)`` —
          :math:`\alpha_{n-1/2}` for :math:`n = 0, \ldots, N`.
        * For cylindrical (``level_indices`` provided): list of
          per-level arrays, each shape ``(M_p + 1,)`` for level
          :math:`p` of size :math:`M_p`.
    redist_dAw : np.ndarray or list[np.ndarray]
        :math:`\Delta A_i / w_n` geometry factor:

        * Sphere: shape ``(nx, N)``.
        * Cylindrical: list of per-level arrays, each shape
          ``(nx, M_p)``.
    tau_mm : np.ndarray or list[np.ndarray]
        Morel--Montry :math:`\tau` clamp values:

        * Sphere: shape ``(N,)``.
        * Cylindrical: list of per-level arrays, each shape
          ``(M_p,)``.

        Read by :class:`MorelMontryAngularSweep` to set the recurrence
        weighting; ignored by :class:`BaileyFlatFluxRedist` (which
        uses the flat-flux collapse, independent of :math:`\tau`).
    volume : np.ndarray
        Cell volumes :math:`V_i`, shape ``(nx,)``.
    level_indices : list[np.ndarray] or None
        Per-:math:`\mu`-level global ordinate index arrays for
        cylindrical geometry.  ``None`` for sphere — the strategy
        then treats all ordinates as a single level.

    Returns
    -------
    np.ndarray
        Redistribution term, shape ``(ng, N, nx)``.  The matvec adds
        this directly to the per-(group, ordinate, cell) streaming +
        collision contribution.

    Attributes
    ----------
    is_linear : bool
        Whether the closure is linear in ``psi_cells``.  Both
        :class:`MorelMontryAngularSweep` and :class:`BaileyFlatFluxRedist`
        are linear (affine combinations of cell-centre values).
        Read-only class attribute.

    Notes
    -----
    The strategy is **stateless** — the same instance can be reused
    across every matvec call without per-call allocation.  Per-call
    work is :math:`\mathcal{O}(N \cdot n_x \cdot n_g)` — the same
    asymptotic cost as the legacy inlined redistribution evaluation.
    """

    is_linear: bool

    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,
    ) -> np.ndarray:
        ...


# ═══════════════════════════════════════════════════════════════════════
# PoleAngularClosureBase — concrete ABC with self-registration
# ═══════════════════════════════════════════════════════════════════════


class PoleAngularClosureBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering pole-angular-closure strategies.

    Subclasses inherit this ABC and pass ``key="..."`` in the class
    statement to self-register; the registry is consulted via
    :meth:`PoleAngularClosureBase.create("morel_montry_angular_sweep")`
    (or any other registered key).

    Subclasses MUST declare:

    * ``is_linear: ClassVar[bool]`` — whether the closure is linear
      in ``psi_cells``.
    * :meth:`__call__` — the closure algorithm.

    Notes
    -----
    Adding a new strategy is a one-line edit::

        class MyClosure(PoleAngularClosureBase, key="my_closure"):
            is_linear: ClassVar[bool] = True

            def __call__(self, psi_cells, alpha_half, redist_dAw,
                          tau_mm, volume, level_indices=None):
                ...

    No registry insert; ``PoleAngularClosureBase.create("my_closure")``
    is immediately callable.
    """

    registry: ClassVar[dict[str, type["PoleAngularClosureBase"]]] = {}

    is_linear: ClassVar[bool]

    @classmethod
    def _registry_base(cls) -> type:
        return PoleAngularClosureBase

    @abstractmethod
    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,
    ) -> np.ndarray:
        ...


# ═══════════════════════════════════════════════════════════════════════
# Helper — single-level DD angular recurrence
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# _MMHalfGrid — module-private typed accessor for the M-M half-angle grid
# ═══════════════════════════════════════════════════════════════════════
#
# PR-TYPED-6.5 Phase 2: the underscore prefix declares "module-private".
# Consumers (matvec, sweep, tests) see only the public API of
# :class:`MorelMontryAngularSweep` and treat the half-grid as opaque
# strategy state. The redistribution body inside the M-M class accesses
# the raw :attr:`faces` array directly; external code consumes via
# :meth:`upstream` / :attr:`upstream_per_ordinate` accessors.


@dataclass(frozen=True, slots=True)
class _MMHalfGrid:
    r"""Typed accessor for the Morel-Montry half-angle face grid.

    Issue #197 PR-TYPED-6c Step 1.5 — Pattern 4 (illegal states
    unrepresentable) for the off-by-one trap the half-angle grid
    historically exposed.

    The M-M recurrence produces :math:`M+1` face fluxes per level for
    :math:`M` ordinates: ``faces[g, 0, i] = ψ_{1/2, i, g}`` (Carlson
    seed, upstream of ordinate 0), ``faces[g, m, i] = ψ_{m-1/2, i, g}``
    (upstream of ordinate m, equivalently downstream of ordinate m-1),
    ``faces[g, M, i] = ψ_{M+1/2, i, g}`` (downstream of last ordinate).

    Two distinct consumers need DIFFERENT slices of this grid:

    * The **redistribution body** (``MorelMontryAngularSweep._weighted_angular_recurrence_single_level``)
      consumes the paired ``(m, m+1)`` access for each ordinate, computing
      :math:`R_m = (\Delta A/w)/V \cdot (\alpha_{m+1/2} \phi_{m+1/2}
      - \alpha_{m-1/2} \phi_{m-1/2})`. Use :attr:`faces` for direct
      paired access.

    * The **unified matvec** (Issue #197 PR-TYPED-6c Step 2+) consumes the
      upstream-per-ordinate slice — one ``(ng, nx)`` block per ordinate
      to populate ``cell_balance_for_streaming``'s ``psi_angular_upstream``
      argument. Use :meth:`upstream` (single ordinate) or
      :attr:`upstream_per_ordinate` (all ordinates).

    The off-by-one trap (``faces[g, m, i]`` vs ``faces[g, m+1, i]``)
    is impossible by API design when consumers use :meth:`upstream` —
    the method's name AND signature enforce upstream-per-ordinate
    semantics. The raw :attr:`faces` array stays exposed so the redist
    body keeps its bit-identical paired access pattern.

    Storage convention (Step 1.5): ``faces`` shape ``(ng, M+1, nx)`` —
    group-leading, matching the existing M-M recurrence kernel. Step
    1.7 (packed-vector canonical alignment) will flip the storage to
    ``(M+1, ng, nx)`` ordinate-leading; this change is deferred so
    Step 1.5 stays purely additive.

    Attributes
    ----------
    faces :
        Shape ``(ng, M+1, nx)``. The raw half-angle grid as produced by
        :meth:`MorelMontryAngularSweep._psi_half_grid_single_level`. ``faces[g, 0, i]`` is
        the Carlson seed at ordinate 0 (= upstream of m=0); for
        ``m = 1, …, M``, ``faces[g, m, i]`` is the half-angle face flux
        :math:`\phi_{m-1/2, i, g}` (upstream of ordinate m, downstream
        of ordinate m-1). ``faces[g, M, i]`` is the downstream face of
        the last ordinate.
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
# MorelMontryAngularSweep — Phase B canonical (default)
# ═══════════════════════════════════════════════════════════════════════
#
# PR-TYPED-6.5 Phase 2.2: the M-M recurrence kernels (formerly the
# free module-level ``_mm_psi_half_grid_single_level`` and
# ``_mm_weighted_angular_recurrence_single_level``) live as private
# staticmethods on the class.  External code never reaches into the
# M-M algebra directly — every M-M consumer goes through the strategy's
# public surface (``__call__`` for the legacy bundle path that retires
# in PR-TYPED-6c Step 7, ``compute_psi_half_per_level`` for the matvec).


class MorelMontryAngularSweep(
    PoleAngularClosureBase, key="morel_montry_angular_sweep",
):  # noqa: E501  (Phase D notes:)
    # Phase D (Issue #168, ERR-026 closure): the M-M recurrence's
    # half-angle seed ``ψ_{1/2,i,g}`` is now sourced from a
    # :class:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeed`
    # strategy field.  The Phase D default :class:`CarlsonInwardSweep`
    # runs Hébert §3.9.4 Eqs. (3.432)-(3.435) inward μ = −1 sweep,
    # giving the canonical Carlson coupled-pole seed that makes the
    # M-M recurrence consistent with the apply matvec on sphere Gate
    # 1.1 MMS.  See the module docstring of
    # :mod:`orpheus.sn.spatial.psi_half_angle_seed` for the rationale.
    r"""Canonical Hébert §3.9.4 per-cell M-M weighted DD angular recurrence.

    PR-TYPED-6.5 Phase 2.3: the strategy is now **bound to an SNMesh at
    construction**.  All M-M coefficients (α-dome, ΔA/w, τ clamp,
    c_in, c_out, level partition) are precomputed from the mesh's
    :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
    eagerly.  The mesh-bound instance methods
    (:meth:`precompute_psi_state`, :meth:`cell_contribution`,
    :meth:`compute_psi_half_per_level`) read this state from ``self``
    — callers never ship M-M data through arguments
    (Pattern 4 — illegal states unrepresentable on "M-M strategy with
    inconsistent coefficients").

    Phase B default for the curvilinear FD operator's angular
    redistribution.  Implements Hébert Eqs. 3.437 / 3.439 with the
    Morel--Montry :math:`\tau` clamp: initialise :math:`\phi_{1/2,i,g}
    = 0` (Carlson zero-weight starting direction), then for :math:`n =
    1, \ldots, N`:

    .. math::

       \phi_{n+1/2,i,g} \;=\;
       \frac{\phi_{n,i,g} \;-\; (1 - \tau_n)\,\phi_{n-1/2,i,g}}{\tau_n}

    and evaluate

    .. math::

       R_{n,i,g} \;=\; \frac{(\Delta A/w)_{i,n}}{V_i}
                       \bigl[\alpha_{n+1/2}\,\phi_{n+1/2,i,g}
                           - \alpha_{n-1/2}\,\phi_{n-1/2,i,g}\bigr].

    At :math:`\tau_n = 1/2` the recurrence reduces to pure DD angular
    (Hébert Eqs. 3.437 / 3.439); the M-M clamp :math:`\tau \in [1/2,
    1]` keeps the M-M weighting positive (Bailey-Morel-Chang 2010).
    The same recurrence runs inside
    :class:`~orpheus.sn.spatial.diamond.DiamondDifference` (the sweep's
    cell update); applying this strategy in the matvec and running
    the sweep solve the **same** discrete fixed point — pinned by
    :file:`tests/sn/l1_analytical/test_pole_closure_sweep_equivalence.py`.

    For cylindrical geometry, the strategy loops over
    :math:`\mu`-levels (each level has its own :math:`\alpha`-dome,
    :math:`\Delta A/w` geometry factor, and :math:`\tau` clamp) and
    runs the recurrence independently per level.  Sphere is the
    single-level (``M_p = N``) case of the same algebra.

    Parameters
    ----------
    sn_mesh : SNMesh, optional
        The mesh + quadrature + materials bundle this strategy binds
        to.

        * When ``sn_mesh`` is supplied (**production path**), M-M
          precomputes α-dome, ΔA/w, τ clamp, c_in, c_out, level
          partition, μ_x, weights, Δr, V at construction.  The
          :meth:`precompute_psi_state` and :meth:`cell_contribution`
          strategy methods (PR-TYPED-6.5+) read these from ``self``
          — no M-M data is shipped through arguments.  This is the
          canonical production mode used by the unified matvec.
        * When ``sn_mesh`` is ``None`` (**legacy/test path**),
          precomputation is skipped.  The mesh-bound methods raise
          ``RuntimeError`` if called; only the legacy :meth:`__call__`
          interface works, taking α / ΔA/w / τ / V via arguments as
          before.  This mode supports tests that verify M-M algebra
          with hand-built coefficient arrays.  Phase 6 retires the
          legacy interface and ``sn_mesh`` becomes required.
    psi_half_seed : PsiHalfAngleSeed, optional
        Strategy producing the half-angle face flux seed
        :math:`\phi_{1/2,i,g}` for the M-M recurrence.  Default is
        :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
        — the canonical Hébert §3.9.4 (3.432)-(3.435) inward μ = −1
        sweep that closes ERR-026 on sphere Gate 1.1 MMS.  Set to
        :class:`~orpheus.sn.spatial.psi_half_angle_seed.ZeroSeed` for
        the regression-safety ablation reproducing Phase B's
        hardcoded zero behaviour.
    """

    is_linear: ClassVar[bool] = True
    """The M-M weighted DD angular recurrence is an affine combination
    of cell-centre values (constant α, ΔA/w, τ coefficients); the
    output is linear in ``psi_cells``."""

    def __init__(
        self,
        sn_mesh: "SNMesh | None" = None,
        *,
        psi_half_seed: PsiHalfAngleSeed | None = None,
    ) -> None:
        # Strategy state.  PR-TYPED-6.5 Phase 2: ``sn_mesh`` is OPTIONAL.
        # When supplied (production), precompute coefficients and bind
        # all mesh-aware methods.  When None (legacy / test), the
        # mesh-bound methods raise; only the legacy ``__call__`` works.
        self._sn_mesh = sn_mesh
        self.psi_half_seed: PsiHalfAngleSeed = (
            psi_half_seed if psi_half_seed is not None else CarlsonInwardSweep()
        )

        if sn_mesh is None:
            # Unbound — only the legacy ``__call__`` (takes args)
            # is functional.  Mesh-bound state is None; mesh-bound
            # methods (precompute_psi_state, cell_contribution,
            # _psi_half_grid_for_level, _redistribution_for_level)
            # raise on call.
            self.level_indices: "tuple[np.ndarray, ...] | None" = None
            self._alpha_per_level: "tuple[np.ndarray, ...] | None" = None
            self._dAw_per_level: "tuple[np.ndarray, ...] | None" = None
            self._tau_per_level: "tuple[np.ndarray, ...] | None" = None
            self._c_in_per_level: "tuple[np.ndarray, ...] | None" = None
            self._c_out_per_level: "tuple[np.ndarray, ...] | None" = None
            self._mu_x: "np.ndarray | None" = None
            self._weights: "np.ndarray | None" = None
            self._dr: "np.ndarray | None" = None
            self._V: "np.ndarray | None" = None
            self._N: "int | None" = None
            return

        coord = sn_mesh.mesh.coord
        quad = sn_mesh.quad
        reduced = sn_mesh.reduced
        N = quad.N

        # ── Per-level partition (M-M's concept, NOT the quadrature's)
        # Sphere: every ordinate is one level (M_p = N, n_levels = 1).
        # Cylinder: μ-levels from ProductQuadrature / LevelSymmetricSN.
        if coord is CoordSystem.SPHERICAL:
            self.level_indices = (np.arange(N),)
            self._alpha_per_level = (reduced.alpha_half,)
            self._dAw_per_level = (reduced.redist_dAw,)
            self._tau_per_level = (reduced.tau_mm,)
        elif coord is CoordSystem.CYLINDRICAL:
            self.level_indices = tuple(
                np.asarray(lvl) for lvl in quad.level_indices
            )
            self._alpha_per_level = tuple(reduced.alpha_per_level)
            self._dAw_per_level = tuple(reduced.redist_dAw_per_level)
            self._tau_per_level = tuple(reduced.tau_mm_per_level)
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

        # ── Carlson-sweep machinery (psi-independent mesh + quad data)
        self._mu_x = quad.mu_x
        self._weights = quad.weights
        self._dr = sn_mesh.dx
        self._V = sn_mesh.volumes[:, 0]
        self._N = N

    # ── Pure-algebra kernel staticmethods (test + legacy surface) ────
    # These were the M-M algebraic core pre-Phase-2.3 (as free
    # module-level functions; PR-TYPED-6.5 Phase 2.2 moved them onto
    # the class as ``@staticmethod``).  Phase 2.3 added the mesh-bound
    # instance methods (``_psi_half_grid_for_level``,
    # ``_redistribution_for_level``) that read from ``self`` —
    # production callers use those.  The staticmethods stay for
    # algebraic-identity tests that verify the kernel under
    # hand-constructed coefficient arrays (no mesh dependence) and
    # for the legacy :meth:`__call__` path.  Same algebra in either
    # surface (Pattern 2 — instance methods delegate to staticmethods
    # with self's data).

    @staticmethod
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

        Pure kernel — accepts all data via arguments.  Used by the
        legacy :meth:`__call__` interface and by algebraic-identity
        tests that verify the recurrence under hand-built coefficients.
        Production code uses :meth:`_psi_half_grid_for_level` which
        reads ``tau_level`` from ``self``.
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

    @staticmethod
    def _weighted_angular_recurrence_single_level(
        psi_level: np.ndarray,
        alpha_level: np.ndarray,
        dAw_level: np.ndarray,
        tau_level: np.ndarray,
        volume: np.ndarray,
        psi_half_seed: np.ndarray | None = None,
    ) -> np.ndarray:
        r"""Pure-algebra M-M redistribution on a single level.

        Composes :meth:`_psi_half_grid_single_level` with the
        ``(ΔA/w)/V`` geometry redistribution coefficient.  Returns the
        per-ordinate redistribution term :math:`R_{m,i,g} =
        (\Delta A/w)_{i,m}/V_i \,(α_{m+1/2}·ψ_{m+1/2} - α_{m-1/2}·
        ψ_{m-1/2})` of shape ``(ng, M, nx)``.

        Pure kernel — accepts all data via arguments.  Counterpart to
        :meth:`_redistribution_for_level` (which reads from self).
        """
        ng, M, nx = psi_level.shape
        psi_half = MorelMontryAngularSweep._psi_half_grid_single_level(
            psi_level, tau_level, psi_half_seed=psi_half_seed,
        )
        redist = np.empty((ng, M, nx), dtype=psi_level.dtype)
        for m in range(M):
            redist[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_level[m + 1] * psi_half[:, m + 1, :]
                   - alpha_level[m] * psi_half[:, m, :])
                / volume.reshape(1, nx)
            )
        return redist

    def _require_mesh_bound(self) -> None:
        """Raise if the strategy was constructed without an sn_mesh."""
        if self._sn_mesh is None:
            raise RuntimeError(
                "MorelMontryAngularSweep was constructed without an SNMesh "
                "(unbound mode for legacy/test compatibility). The "
                "mesh-bound methods (precompute_psi_state, cell_contribution, "
                "_psi_half_grid_for_level, _redistribution_for_level) require "
                "an SNMesh — instantiate via MorelMontryAngularSweep(sn_mesh)."
            )

    # ── Recurrence kernels (PR-TYPED-6.5 Phase 2.3) ──────────────────
    # Mesh-bound instance methods. They read ``tau`` / ``alpha`` / ``dAw``
    # / ``V`` from ``self`` — no shipping of mesh data through
    # arguments. The only call-time input that varies across matvec
    # calls is ``psi_level`` (the angular-flux slice for one level)
    # and the optional ``psi_half_seed`` produced by the Carlson sweep.

    def _psi_half_grid_for_level(
        self,
        psi_level: np.ndarray,           # (ng, M_p, nx) — ordinates within level p
        level_idx_p: int,
        *,
        psi_half_seed: np.ndarray | None = None,  # (ng, nx)
    ) -> np.ndarray:
        r"""Run the M-M angular recurrence on level ``level_idx_p``.

        Mesh-bound wrapper around :meth:`_psi_half_grid_single_level`.
        Reads ``τ`` from ``self._tau_per_level[level_idx_p]`` — no
        coefficient data is shipped through arguments.

        Returns
        -------
        np.ndarray
            Half-angle grid, shape ``(ng, M_p+1, nx)``.
        """
        self._require_mesh_bound()
        return self._psi_half_grid_single_level(
            psi_level,
            self._tau_per_level[level_idx_p],
            psi_half_seed=psi_half_seed,
        )

    def _redistribution_for_level(
        self,
        psi_level: np.ndarray,           # (ng, M_p, nx)
        level_idx_p: int,
        *,
        psi_half_seed: np.ndarray | None = None,  # (ng, nx)
    ) -> np.ndarray:
        r"""Hébert's per-cell M-M redistribution :math:`R_{m,i,g}` for level ``p``.

        Mesh-bound wrapper around
        :meth:`_weighted_angular_recurrence_single_level`.  Reads
        ``α``, ``ΔA/w``, ``τ``, ``V`` from ``self`` — no coefficient
        data is shipped.

        Returns
        -------
        np.ndarray
            Redistribution per level, shape ``(ng, M_p, nx)``.
        """
        self._require_mesh_bound()
        return self._weighted_angular_recurrence_single_level(
            psi_level,
            self._alpha_per_level[level_idx_p],
            self._dAw_per_level[level_idx_p],
            self._tau_per_level[level_idx_p],
            self._V,
            psi_half_seed=psi_half_seed,
        )

    # ── Strategy Protocol surface (PR-TYPED-6.5 Phase 2.5+) ──────────
    # The matvec body consumes ``precompute_psi_state`` + ``cell_contribution``
    # to read M-M contributions per-cell without naming any M-M algebra
    # or coefficient (Pattern 1 — operator algebra reads as composition).

    def precompute_psi_state(
        self,
        psi_view: np.ndarray,                       # (N, ng, nx, 1) canonical
        *,
        sigma_t: np.ndarray,                        # (ng, nx)
        bc_outer_inflow_estimate: np.ndarray,       # (N, ng)
    ) -> "tuple[_MMHalfGrid, ...]":
        r"""Build per-level half-angle grids with Carlson coupled-pole seeds.

        One ``_MMHalfGrid`` per level (sphere: 1 element; cylinder:
        ``n_levels`` elements).  Matvec consumes via :meth:`cell_contribution`.

        Parameters
        ----------
        psi_view :
            Current angular-flux iterate, canonical layout.
        sigma_t :
            Per-group per-cell total cross section ``(ng, nx)``.
            Needed by the Carlson coupled-pole sweep.
        bc_outer_inflow_estimate :
            Boundary-trace estimate at the outer face, ``(N, ng)``.
            One row per ordinate; ``CarlsonSweepContext.bc_outer_value``
            consumes the row matching the level's most-inward ordinate.
        """
        self._require_mesh_bound()
        # (N, ng, nx, 1) → (ng, N, nx) — drop y; reorder for level access.
        psi_g_first = psi_view.transpose(1, 0, 2, 3)[..., 0]
        per_level: list[_MMHalfGrid] = []
        for p, level_idx in enumerate(self.level_indices):
            mu_level = self._mu_x[level_idx]
            within_inward = int(np.argmin(mu_level))
            global_inward = int(level_idx[within_inward])
            bc_value = bc_outer_inflow_estimate[global_inward, :]
            psi_level = psi_g_first[:, level_idx, :]  # (ng, M_p, nx)
            carlson_ctx = CarlsonSweepContext(
                sigma_t=sigma_t,
                dr=self._dr,
                mu_quad=mu_level.copy(),
                weights=self._weights[level_idx].copy(),
                bc_outer_value=bc_value,
            )
            psi_half_seed_arr = self.psi_half_seed(psi_level, carlson_ctx)
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
        self._require_mesh_bound()
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

    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,
    ) -> np.ndarray:
        r"""Compute the Hébert §3.9.4 angular redistribution term.

        Legacy bundle interface — accepts all coefficient arrays via
        arguments.  Works regardless of whether the strategy was
        constructed with an SNMesh (production) or without (test).
        Dispatches by ``level_indices``: ``None`` → sphere (single
        level over all ordinates); populated → cylindrical (per-level
        loop).

        The production matvec consumes the new mesh-bound interface
        :meth:`precompute_psi_state` + :meth:`cell_contribution` and
        does NOT call this method; the sweep route and the legacy
        matvec helpers (``transport_operator_matvec_spherical`` /
        ``_cylindrical``) do until PR-TYPED-6c Step 7 / PR-TYPED-6.5
        Phase 6 retires them.

        Parameters
        ----------
        carlson_context :
            Phase D Carlson coupled-pole sweep inputs (Σ_t, Δr,
            μ_quad, weights, bc_outer_value).  See
            :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonSweepContext`.
            When supplied, the M-M recurrence's seed comes from
            ``self.psi_half_seed(psi_level, carlson_context)``.  When
            ``None``, the recurrence falls back to the Phase B
            hardcoded zero seed.  For cylindrical (``level_indices``
            populated), a list of per-level contexts is expected; for
            spherical (``level_indices is None``), a single context.
        """
        if level_indices is None:
            # Spherical: ordinates ARE a single level.
            assert isinstance(alpha_half, np.ndarray), (
                "Spherical pole closure expects alpha_half to be an "
                "ndarray; got list (cylindrical inputs)."
            )
            assert isinstance(redist_dAw, np.ndarray), (
                "Spherical pole closure expects redist_dAw to be an "
                "ndarray; got list (cylindrical inputs)."
            )
            assert isinstance(tau_mm, np.ndarray), (
                "Spherical pole closure expects tau_mm to be an "
                "ndarray; got list (cylindrical inputs)."
            )
            psi_half_seed_arr: np.ndarray | None = None
            if carlson_context is not None:
                assert isinstance(carlson_context, CarlsonSweepContext), (
                    "Spherical pole closure expects carlson_context "
                    "to be a single CarlsonSweepContext (not a list)."
                )
                psi_half_seed_arr = self.psi_half_seed(
                    psi_cells, carlson_context,
                )
            return self._weighted_angular_recurrence_single_level(
                psi_cells, alpha_half, redist_dAw, tau_mm, volume,
                psi_half_seed=psi_half_seed_arr,
            )

        # Cylindrical: loop over μ-levels.
        assert isinstance(alpha_half, list)
        assert isinstance(redist_dAw, list)
        assert isinstance(tau_mm, list)
        if carlson_context is not None:
            assert isinstance(carlson_context, list), (
                "Cylindrical pole closure expects carlson_context to "
                "be a list of per-level CarlsonSweepContext."
            )
            assert len(carlson_context) == len(level_indices), (
                "carlson_context length must match level_indices length"
            )
        ng, _N, nx = psi_cells.shape
        redist = np.zeros((ng, _N, nx), dtype=psi_cells.dtype)
        for p, level_idx in enumerate(level_indices):
            psi_level = psi_cells[:, level_idx, :]            # (ng, M_p, nx)
            psi_half_seed_arr = None
            if carlson_context is not None:
                psi_half_seed_arr = self.psi_half_seed(
                    psi_level, carlson_context[p],
                )
            redist_level = self._weighted_angular_recurrence_single_level(
                psi_level, alpha_half[p], redist_dAw[p], tau_mm[p],
                volume, psi_half_seed=psi_half_seed_arr,
            )
            redist[:, level_idx, :] = redist_level
        return redist

    # ── Issue #197 PR-TYPED-6b: ψ_half grid exposure ──────────────────

    def compute_psi_half_per_level(
        self,
        psi_level: np.ndarray,
        tau_level: np.ndarray,
        *,
        carlson_context: "CarlsonSweepContext | None" = None,
    ) -> _MMHalfGrid:
        r"""Return the half-angle grid :math:`\phi_{m\pm 1/2, i, g}`
        for one level under the M-M recurrence, wrapped in the private
        :class:`_MMHalfGrid` typed accessor.

        Legacy interface (PR-TYPED-6b) — accepts ``tau_level`` via
        argument so algebraic-identity tests can verify the recurrence
        under hand-built ``τ``.  Production code in the unified matvec
        uses :meth:`precompute_psi_state` (which reads ``τ`` from
        ``self._tau_per_level``).

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
        carlson_context :
            Optional Phase D Carlson coupled-pole seed context.  When
            supplied, the recurrence seeds at
            :math:`\phi_{1/2, i, g} = \mathrm{Carlson}(\psi_{\rm level},
            \mathrm{ctx})` via :attr:`psi_half_seed` (default
            :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`).
            When ``None`` the recurrence falls back to the Phase B
            hardcoded zero seed.

        Returns
        -------
        _MMHalfGrid
            Typed accessor wrapping the half-angle grid
            ``faces`` of shape ``(ng, M_p+1, nx)``.
        """
        psi_half_seed_arr: np.ndarray | None = None
        if carlson_context is not None:
            psi_half_seed_arr = self.psi_half_seed(psi_level, carlson_context)
        faces = self._psi_half_grid_single_level(
            psi_level, tau_level, psi_half_seed=psi_half_seed_arr,
        )
        return _MMHalfGrid(faces=faces)

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
# PR-TYPED-6.5 Phase 2.8.  Cartesian slab carries no angular
# redistribution term (no curvature → no Hébert §3.9.4 closure).
# Earlier code modelled this as ``sn_mesh.pole_angular_closure is None``
# inside the matvec body — a Pattern 4 violation (the matvec's
# ``if closure is None`` branch was an illegal-state check rather than
# a typed dispatch).  The Identity strategy makes the slab algebra a
# typed default: the strategy exists, has the same Protocol surface,
# and contributes zero to both the cell-balance denominator and the
# upstream numerator.


class IdentityAngularClosure(PoleAngularClosureBase, key="identity_angular_closure"):
    r"""No-op pole-angular closure for Cartesian (slab + 2-D rectilinear).

    PR-TYPED-6.5 Phase 2.8.  The Cartesian SN balance equation has no
    angular-redistribution term — Hébert §3.9.4's :math:`(\Delta A/w)`
    factor vanishes on flat geometry (cell faces are parallel, no
    curvature-driven coupling between consecutive ordinate sub-domains).
    This strategy returns ``(0, 0)`` from :meth:`cell_contribution`,
    contributing nothing to the matvec's per-cell denominator and
    upstream numerator.  The matvec body therefore consumes the SAME
    ``cell_balance_for_streaming`` algebra for Cartesian as for sphere
    + cylinder — geometry-blind by data (Cardinal Rule 2).

    Why it earns its keep
    ---------------------
    Before Phase 2.8 the matvec body's seed branch read

    .. code-block:: python

       if curvature == "cartesian":
           pole_face_seed = bc_inner.apply(...)
       else:
           pole_face_seed = psi_view[..., 0, 0].copy()

    The ``if curvature == "cartesian"`` was a Pattern-4 leak — geometry
    dispatch baked into the matvec body.  With Identity in place, the
    branch becomes ``if sn_mesh.bc_left is None`` (PR-TYPED-6.5 Phase
    3a) and ``sn_mesh.pole_angular_closure`` is ALWAYS a valid object
    (no ``is None`` test).

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

    def __init__(self, sn_mesh: "SNMesh") -> None:
        self._sn_mesh = sn_mesh
        self._N: int = sn_mesh.quad.N
        self._ng: int = sn_mesh.ng
        nx = sn_mesh.nx
        # Trivial single-level partition: every ordinate in one level.
        self.level_indices: tuple[np.ndarray, ...] = (np.arange(self._N),)
        # Neutral M-M coefficients per level — Cartesian carries no
        # angular redistribution, so α_{±1/2} = 0, τ = 1, ΔA/w = 0, and
        # the algebraic constants c_in = c_out = 0.  Exposed under the
        # same private names M-M uses so the matvec body can read from
        # the closure without geometry dispatch (PR-TYPED-6.5 Phase 3a.3).
        self._alpha_per_level: tuple[np.ndarray, ...] = (
            np.zeros(self._N + 1),
        )
        self._dAw_per_level: tuple[np.ndarray, ...] = (
            np.zeros((nx, self._N)),
        )
        self._tau_per_level: tuple[np.ndarray, ...] = (
            np.ones(self._N),
        )
        self._c_in_per_level: tuple[np.ndarray, ...] = (
            np.zeros(self._N),
        )
        self._c_out_per_level: tuple[np.ndarray, ...] = (
            np.zeros(self._N),
        )

    # ── Strategy Protocol surface ────────────────────────────────────

    def precompute_psi_state(
        self,
        psi_view: np.ndarray,
        *,
        sigma_t: np.ndarray | None = None,
        bc_outer_inflow_estimate: np.ndarray | None = None,
    ) -> None:
        """No state — Cartesian has no curvature half-grid to precompute."""
        del psi_view, sigma_t, bc_outer_inflow_estimate
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

    # ── Legacy bundle surface (Protocol-shape compatibility) ─────────

    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray] | None" = None,
        redist_dAw: "np.ndarray | list[np.ndarray] | None" = None,
        tau_mm: "np.ndarray | list[np.ndarray] | None" = None,
        volume: "np.ndarray | None" = None,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: Any = None,
    ) -> np.ndarray:
        """Slab has no redistribution — returns zeros shaped like ``psi_cells``."""
        del alpha_half, redist_dAw, tau_mm, volume, level_indices, carlson_context
        return np.zeros_like(psi_cells)

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


# ═══════════════════════════════════════════════════════════════════════
# BaileyFlatFluxRedist — legacy reproducer for ablation / back-bisection
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class BaileyFlatFluxRedist(
    PoleAngularClosureBase, key="bailey_flat_flux_redist",
):
    r"""Legacy flat-flux-collapsed redistribution — Issue #168 Defect 3.

    Reproduces the **pre-Phase-B** angular redistribution of the
    curvilinear FD operator.  The half-angle face fluxes were
    collapsed to cell-centre values by the legacy τ-weighted
    interpolation,

    .. math::

       \phi_{n\pm 1/2,i,g} \;\approx\; \phi_{n,i,g}

    (more precisely: the symmetric Morel--Montry interpolation
    :math:`\phi_{n+1/2} = \tau\,\phi_{n+1} + (1-\tau)\,\phi_n`
    collapses to :math:`\phi_n` on angularly-flat :math:`\psi`).
    The redistribution then becomes

    .. math::

       R_{n,i,g} \;=\; \frac{(\Delta A/w)_{i,n}}{V_i}
                       \bigl(\alpha_{n+1/2} - \alpha_{n-1/2}\bigr)
                       \phi_{n,i,g}
                       \;=\; -\frac{\mu_n\,\Delta A_i}{V_i}\,\phi_{n,i,g}

    using :math:`\alpha_{n+1/2} - \alpha_{n-1/2} = -w_n\,\mu_n`
    (Bailey 2009 dome recursion).  This is **exact** when
    :math:`\psi` is constant in :math:`\mu` (the per-ordinate
    flat-flux consistency invariant Bailey's :math:`\Delta A/w`
    factor enforces) but **only :math:`\mathcal{O}(1)` accurate** on
    angularly-varying :math:`\psi` — the factor-of-two truncation
    gap identified as Issue #168 Defect 3.

    Why this strategy exists
    ------------------------

    * **Ablation studies**: comparing
      :class:`MorelMontryAngularSweep` vs :class:`BaileyFlatFluxRedist`
      quantifies the Phase-B improvement on angularly-varying
      :math:`\psi`.
    * **Back-bisection**: if a future change to the matvec causes a
      regression, swapping in :class:`BaileyFlatFluxRedist`
      reproduces the historical pre-Phase-B operator output and
      pins which step in a multi-commit refactor introduced the
      drift.
    * **Documentation**: the strategy carries the bug it reproduces
      in its docstring.  Issue #168 Defect 3 is now a strategy
      pluggable into the matvec — anyone reading
      :mod:`orpheus.sn.spatial.pole_angular_closure` learns what the
      bug was, not just what the fix is.

    This is **not** a default; production callers MUST use
    :class:`MorelMontryAngularSweep`.

    Notes
    -----
    The legacy form built into :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
    pre-Phase-B used τ-weighted symmetric interpolation rather than
    the bare cell-centre substitution; on angularly-flat :math:`\psi`
    these give bit-identical results (the τ-weighted form collapses
    to :math:`\phi_n`).  This strategy reproduces the **collapse**
    semantics directly so that the Phase-B flat-flux identity test
    can pin the consistency between the two strategies on flat ψ
    without going through the matvec.
    """

    is_linear: ClassVar[bool] = True
    """Identity-like collapse on cell-centre values is trivially linear."""

    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",  # accepted but unused
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,  # noqa: E501
    ) -> np.ndarray:
        r"""Reproduce the flat-flux-collapsed redistribution.

        :math:`R_{n,i,g} = (\Delta A_i / w_n) (\alpha_{n+1/2} -
        \alpha_{n-1/2}) \phi_{n,i,g} / V_i`.

        ``tau_mm`` and ``carlson_context`` are accepted for
        Protocol-shape compatibility with
        :class:`MorelMontryAngularSweep` but are **not read** —
        :class:`BaileyFlatFluxRedist` collapses the half-angle face
        fluxes to cell centres unconditionally, independent of the
        M-M clamp and the Phase D Carlson seed.
        """
        del tau_mm  # explicitly unused — Protocol-shape compatibility only
        del carlson_context  # explicitly unused — Phase D compat only
        if level_indices is None:
            assert isinstance(alpha_half, np.ndarray)
            assert isinstance(redist_dAw, np.ndarray)
            return self._collapse_single_level(
                psi_cells, alpha_half, redist_dAw, volume,
            )

        assert isinstance(alpha_half, list)
        assert isinstance(redist_dAw, list)
        ng, N, nx = psi_cells.shape
        redist = np.zeros((ng, N, nx), dtype=psi_cells.dtype)
        for p, level_idx in enumerate(level_indices):
            psi_level = psi_cells[:, level_idx, :]
            redist_level = self._collapse_single_level(
                psi_level, alpha_half[p], redist_dAw[p], volume,
            )
            redist[:, level_idx, :] = redist_level
        return redist

    @staticmethod
    def _collapse_single_level(
        psi_level: np.ndarray,
        alpha_level: np.ndarray,
        dAw_level: np.ndarray,
        volume: np.ndarray,
    ) -> np.ndarray:
        ng, M, nx = psi_level.shape
        redist = np.empty((ng, M, nx), dtype=psi_level.dtype)
        for m in range(M):
            # R = (ΔA/w) × (α_out - α_in) × ψ_n / V
            d_alpha = alpha_level[m + 1] - alpha_level[m]
            redist[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * d_alpha
                * psi_level[:, m, :]
                / volume.reshape(1, nx)
            )
        return redist

    def __repr__(self) -> str:
        return "BaileyFlatFluxRedist()"


# ═══════════════════════════════════════════════════════════════════════
# LegacyTauSymmetricInterpolation — pre-Phase-B inlined form
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class LegacyTauSymmetricInterpolation(
    PoleAngularClosureBase, key="legacy_tau_symmetric",
):
    r"""Pre-Phase-B inlined :math:`\tau`-symmetric interpolation form.

    Reproduces the **exact** angular redistribution algebra that lived
    inline inside
    :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
    (and ``_cylindrical``) before Issue #168 Phase B.  The half-angle
    face fluxes were evaluated by τ-symmetric interpolation between
    consecutive cell-centres,

    .. math::

       \phi_{n+1/2,i,g} \;=\; \tau_n\,\phi_{n+1,i,g}
                              + (1 - \tau_n)\,\phi_{n,i,g}
                              \quad (\text{for } n < N - 1),

    .. math::

       \phi_{n-1/2,i,g} \;=\; \tau_{n-1}\,\phi_{n,i,g}
                              + (1 - \tau_{n-1})\,\phi_{n-1,i,g}
                              \quad (\text{for } n \ge 1),

    with the boundary conventions
    :math:`\phi_{N+1/2,i,g} = \phi_{N,i,g}` and
    :math:`\phi_{1/2,i,g} = \phi_{0,i,g}` (cell-centre extrapolation
    at the ordinate endpoints).  On angularly-flat :math:`\psi` this
    collapses to :math:`\phi_{n\pm 1/2,i,g} = \phi_{n,i,g}` —
    bit-identical to :class:`BaileyFlatFluxRedist`'s collapse — but
    on angularly-varying :math:`\psi` the τ-symmetric interpolation
    gives a quantitatively different result from either
    :class:`BaileyFlatFluxRedist` or :class:`MorelMontryAngularSweep`.

    This strategy IS the **Phase B default** — it preserves bit-for-bit
    agreement with the pre-Issue-#168 operator output, which is what
    the curvilinear regression-snapshot bit-identity contract demands.
    The Defect-3 truncation gap on angularly-varying :math:`\psi`
    survives under this strategy by design (the bug is *reproducible*
    so it can be cross-checked against Phase A's empirical behaviour);
    fixing it requires opting in to
    :class:`MorelMontryAngularSweep` and is gated on the spatial-closure
    follow-up beyond Phase B's scope.

    Why three strategies (not two)
    ------------------------------

    * :class:`LegacyTauSymmetricInterpolation` — reproduces the
      pre-Phase-B inlined math exactly.  **Default**.  Bit-identical
      regression preservation.  Carries Defect 3 by design.
    * :class:`BaileyFlatFluxRedist` — algebraic collapse equivalent
      (only on flat ψ); used for the ablation L1 test that pins the
      flat-flux equivalence between the legacy form and the Bailey
      collapse.
    * :class:`MorelMontryAngularSweep` — canonical Hébert §3.9.4
      form.  Opt-in.  Closes Defect 3 on angularly-varying ψ but
      requires spatial-closure alignment for full ERR-026 closure.

    Notes
    -----
    Frozen + slotted: instances are immutable and lightweight; a
    single :class:`LegacyTauSymmetricInterpolation` can be reused
    across every matvec call without per-call allocation.
    """

    is_linear: ClassVar[bool] = True
    """The τ-symmetric interpolation is an affine combination of
    cell-centre values; the redistribution output is linear in
    ``psi_cells``."""

    def __call__(
        self,
        psi_cells: np.ndarray,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,  # noqa: E501
    ) -> np.ndarray:
        r"""Reproduce the pre-Phase-B inlined τ-symmetric form.

        ``carlson_context`` is accepted for Protocol-shape
        compatibility with :class:`MorelMontryAngularSweep` but is
        **not read** — :class:`LegacyTauSymmetricInterpolation`
        evaluates half-angle face fluxes via τ-symmetric interpolation
        between cell centres, independent of any Carlson seed.
        """
        del carlson_context  # explicitly unused — Phase D compat only
        if level_indices is None:
            assert isinstance(alpha_half, np.ndarray)
            assert isinstance(redist_dAw, np.ndarray)
            assert isinstance(tau_mm, np.ndarray)
            return self._tau_symmetric_single_level(
                psi_cells, alpha_half, redist_dAw, tau_mm, volume,
            )

        assert isinstance(alpha_half, list)
        assert isinstance(redist_dAw, list)
        assert isinstance(tau_mm, list)
        ng, N, nx = psi_cells.shape
        redist = np.zeros((ng, N, nx), dtype=psi_cells.dtype)
        for p, level_idx in enumerate(level_indices):
            psi_level = psi_cells[:, level_idx, :]
            redist_level = self._tau_symmetric_single_level(
                psi_level, alpha_half[p], redist_dAw[p], tau_mm[p],
                volume,
            )
            redist[:, level_idx, :] = redist_level
        return redist

    @staticmethod
    def _tau_symmetric_single_level(
        psi_level: np.ndarray,           # (ng, M, nx)
        alpha_level: np.ndarray,         # (M+1,)
        dAw_level: np.ndarray,           # (nx, M)
        tau_level: np.ndarray,           # (M,)
        volume: np.ndarray,              # (nx,)
    ) -> np.ndarray:
        ng, M, nx = psi_level.shape
        redist = np.empty((ng, M, nx), dtype=psi_level.dtype)
        for m in range(M):
            tau_m = tau_level[m]
            # Pre-Phase-B inlined ψ_face_right (psi_angle_right):
            #   ψ_{m+1/2,i,g} = τ_m·ψ_{m+1,i,g} + (1-τ_m)·ψ_{m,i,g}
            # for m < M-1; at m = M-1 the legacy code substituted the
            # cell-centre (no upper neighbour available).
            if m < M - 1:
                psi_face_right = (
                    tau_m * psi_level[:, m + 1, :]
                    + (1.0 - tau_m) * psi_level[:, m, :]
                )
            else:
                psi_face_right = psi_level[:, m, :]
            # Pre-Phase-B inlined ψ_face_left (psi_angle_left):
            #   ψ_{m-1/2,i,g} = τ_{m-1}·ψ_{m,i,g}
            #                   + (1-τ_{m-1})·ψ_{m-1,i,g}
            # for m > 0; at m = 0 the legacy code substituted the
            # cell-centre (no lower neighbour available).
            if m > 0:
                tau_prev = tau_level[m - 1]
                psi_face_left = (
                    tau_prev * psi_level[:, m, :]
                    + (1.0 - tau_prev) * psi_level[:, m - 1, :]
                )
            else:
                psi_face_left = psi_level[:, m, :]
            # Redistribution:
            # R = (ΔA/w_m) × (α_{m+1/2}·ψ_{m+1/2} - α_{m-1/2}·ψ_{m-1/2}) / V
            redist[:, m, :] = (
                dAw_level[:, m].reshape(1, nx)
                * (alpha_level[m + 1] * psi_face_right
                   - alpha_level[m] * psi_face_left)
                / volume.reshape(1, nx)
            )
        return redist

    def __repr__(self) -> str:
        return "LegacyTauSymmetricInterpolation()"


__all__ = [
    "BaileyFlatFluxRedist",
    "CarlsonInwardSweep",
    "CarlsonSweepContext",
    "IdentityAngularClosure",
    "LegacyTauSymmetricInterpolation",
    "MorelMontryAngularSweep",
    "PoleAngularClosure",
    "PoleAngularClosureBase",
    "PsiHalfAngleSeed",
    "default_angular_closure_class",
]
