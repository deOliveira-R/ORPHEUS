r"""Typed return type for the SN transport solvers.

Issue #197 PR-TYPED-5 — :class:`Solution` and :class:`IterationHistory`
replaced the legacy bare-dataclass pair ``SNFixedSourceResult`` /
``SNResult`` (both formerly in ``orpheus.sn.solver``, deleted with the
migration).

The legacy types were data bags: they carried bare ``np.ndarray`` flux
fields, an opaque ``geometry`` / ``quadrature`` pair, and an ad-hoc set
of diagnostic scalars (``n_inner``, ``residual``, ``keff_history``,
``elapsed_seconds``).  Two distinct dataclasses for the two problem
kinds (fixed-source vs eigenvalue) duplicated the shape contract in two
places — a twin path waiting for drift (``coding-elegance`` anti-pattern
1: "two implementations of the same mathematical quantity").

Under PR-TYPED-5, **one carrier** (:class:`SolutionBase`, through both
role leaves) covers both problem kinds.  The kind discrimination lives
in two methods — :meth:`SolutionBase.is_eigenvalue` and
:meth:`SolutionBase.is_fixed_source` — that read the optional
:attr:`keff`.  The convergence-trajectory diagnostics live on a
separate :class:`IterationHistory`, exposed through method-style
accessors :meth:`SolutionBase.dominance_ratio` and
:meth:`SolutionBase.converged`.

Reads as the math (``coding-elegance`` Pattern 1 — match the algebra of
the domain)::

    sol = solve_sn(...)
    sol.is_eigenvalue()                    # True / False
    sol.dominance_ratio()                  # |k_n / k_{n-1}|
    sol.reaction_rate_density(sig_a)       # σ · φ
    sol.compare(other, rtol=1e-12)         # SolutionDiff

Two discrimination axes
=======================

The solution family discriminates along TWO independent axes, and the
axes deliberately use DIFFERENT mechanisms:

* **Problem kind** (fixed-source vs eigenvalue) is a **property** — one
  carrier covers both kinds via the optional :attr:`SolutionBase.keff`,
  read through :meth:`SolutionBase.is_eigenvalue` /
  :meth:`SolutionBase.is_fixed_source`.  The two kinds share every
  realization AND every operation (homogenizing a fixed-source flux is
  as meaningful as homogenizing an eigenmode), so a type here would be
  ceremony (the type-minting criterion fails on both prongs).

* **Solution role** (forward vs adjoint) is a **type** —
  :class:`SolutionBase` → {:class:`Solution`, :class:`AdjointSolution`}
  (campaign #276 A5 ruling, 2026-07-25).  The roles share the carrier
  (same fields, same packaging convention) but NOT the operation set:
  the forward-physics methods (:meth:`Solution.homogenize`,
  :meth:`Solution.condense`, :meth:`Solution.reaction_rate_density`)
  interpret :attr:`~SolutionBase.scalar_flux` as the flux :math:`\phi`
  and are physically meaningless on the importance :math:`\varphi^*` —
  homogenization/condensation collapse cross sections *preserving
  reaction rates*, an operation ON the forward flux; the adjoint enters
  only as the Petrov-Galerkin test weight that refines the collapse
  (the #281 P6-B2 parameter ``adjoint: AdjointSolution | None`` —
  LANDED, with the worth-zeroing taxonomy of
  :mod:`orpheus.derivations.common.homogenization`), never as its
  subject.  The type split makes the wrong physics UNSPELLABLE — an
  :class:`AdjointSolution` has no ``homogenize`` attribute at all
  (structural absence, not a runtime refusal) — and gives the adjoint
  machinery family (the landed adjoint-weighted collapse; perturbation
  theory :math:`\langle\varphi^*, \delta A\, \varphi\rangle` and
  generalized perturbation / response estimation to come) its
  signature-level carrier.

Mesh-binding consistency
========================

:class:`SolutionBase.__post_init__` validates that every typed field
(:attr:`angular_flux`, :attr:`scalar_flux`, :attr:`boundary_flux`)
carries the SAME :class:`SNMesh` instance as the
:attr:`SolutionBase.mesh`.  Mesh identity is checked via ``is`` — sharing
a mesh by value rather than by reference is forbidden by construction
(``coding-elegance`` Pattern 4 — illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Self

import numpy as np

if TYPE_CHECKING:
    from .mesh.augmented_mesh import SNMesh
    from orpheus.data.energy_grid import EnergyGrid, WithinGroupSpectrum
    from orpheus.numerics.convergence import IterationRecord
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry import Mesh1D, Mesh2D
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )
    from orpheus.transport.timed_full_field import TimedFullField


__all__ = [
    "AdjointSolution",
    "IterationHistory",
    "Solution",
    "SolutionBase",
    "SolutionDiff",
]


@dataclass(frozen=True)
class IterationHistory:
    r"""Convergence diagnostics for an SN solve — a **view** over the record.

    The solve's truth is an
    :class:`~orpheus.numerics.convergence.IterationRecord`: a tree of
    iteration levels, each carrying the quantities it stopped on with the
    tolerances it judged them against.  This type is the SN-facing reading
    of that tree — every scalar below is DERIVED, so there is exactly one
    source of truth and the flat surface cannot drift from the tree it
    summarises.

    ⛔ Until 2026-08-09 (#340 N2b-ii) these were six independent FIELDS and
    each producer filled them in by hand.  That is how ``n_inner`` came to
    mean ``len(residuals)`` on one path and ``len(residuals) + 1`` on the
    other, undocumented and exactly backwards from the truth; and how
    ``converged`` came to be written five times, once as a literal ``True``
    (#342).  A projection maintained by hand at N sites is N chances to
    disagree; a projection computed from one object is none.

    Parameters
    ----------
    record : IterationRecord
        What the solve measured — **required**, and the only input.  For a
        fixed-source solve this is the inner driver's own record (a leaf);
        for an eigenvalue solve it is the power iteration's outer record,
        whose ``children`` are the per-outer inner solves.
    keff_history : tuple of float
        Per-outer-iteration eigenvalue trajectory.  Empty for fixed-source
        problems.  This is a **field, not a derived reading**, because it is
        a physics output rather than a stopping criterion: what the outer
        stop test reads is its per-iteration INCREMENT, which lives in the
        record under the name the solver gave it (``dk``).

    Notes
    -----
    The trajectories are **tuple** rather than list — frozen dataclasses
    with mutable fields are an anti-pattern (``coding-elegance``
    Pattern 4: illegal states unrepresentable).  Callers that want a
    list pass through :func:`list` at the call site.

    ⭐ **Read :attr:`record` directly for anything this view flattens
    away** — which is most of it.  The flat readings below exist for
    established consumers; the record answers *which* criterion bound, at
    what rate it was closing, what budget would reach it, and which nested
    level actually failed.  Do not grow this surface: add the question to
    :class:`~orpheus.numerics.convergence.IterationRecord`, where the tree
    can answer it.
    """

    record: "IterationRecord"
    keff_history: tuple[float, ...] = ()

    # ── the verdict, derived ─────────────────────────────────────────────

    @property
    def converged(self) -> bool:
        """Did the solve's TOP level meet its own criteria?

        ⛔ A field until 2026-08-09, and before 2026-08-08 a field
        *defaulting to* ``True`` — the same defect as #342 with the
        assertion moved into the type: a history built by a producer that
        had not thought about convergence claimed it anyway.  Deriving it
        removes the question of who writes it.

        ⚠ Scoped to the top level, which for an eigenvalue solve is the
        OUTER.  A converged outer standing on a starved inner reads ``True``
        here and ``False`` at :attr:`fully_converged`; that gap is the #340
        headline, not a wart.  **A gate asserting physics wants
        :attr:`fully_converged`.**
        """
        return self.record.converged

    @property
    def fully_converged(self) -> bool:
        """Did EVERY level converge — the honest "can I trust this number?".

        New with the record (#340): the flat history had no way to express
        it, because it had already discarded the nested levels.
        """
        return self.record.fully_converged

    # ── the flat readings established consumers still take ───────────────

    @property
    def _is_outer(self) -> bool:
        """Does the top level drive nested solves, or IS it the solve?

        The discriminator the readings below need, taken from the tree's own
        STRUCTURE rather than from a label — a level with children is an
        outer over inners; a leaf is the within-group solve itself.  Reading
        it off ``label`` would be stringly-typed dispatch on a string chosen
        for humans.
        """
        return bool(self.record.children)

    @property
    def flux_residuals(self) -> tuple[float, ...]:
        r"""Per-iteration relative-flux residual trajectory
        :math:`\lVert\phi_n - \phi_{n-1}\rVert / \lVert\phi_n\rVert`.

        Empty on the eigenvalue path, whose top level stops on ``dk`` and
        ``dphi`` rather than on a within-group residual — read
        ``record.criteria`` (or ``record.children``) there.  That emptiness
        is the pre-existing behaviour, kept deliberately: widening this name
        to mean "``dphi`` when there is no residual" would silently re-point
        every consumer that branches on it.
        """
        if self._is_outer:
            return ()
        criterion = self.record.binding_criterion
        return () if criterion is None else criterion.trajectory

    @property
    def n_inner(self) -> int | None:
        """Inner (within-group) iterations consumed by THE inner solve.

        ``None`` on the eigenvalue path — there is no single inner solve to
        point at, which is why :attr:`total_inner_iterations` exists.
        """
        return None if self._is_outer else self.record.n_iterations

    @property
    def total_inner_iterations(self) -> int | None:
        """Inner iterations summed across ALL outer iterations.

        Populated by both paths: the eigenvalue path sums its children, the
        fixed-source path reports its single solve.  It is the measurand for
        the SI spectral-rate / Gauss-Seidel-recovery diagnostics.
        """
        if self._is_outer:
            return sum(child.n_iterations for child in self.record.children)
        return self.record.n_iterations

    @property
    def n_outer(self) -> int | None:
        """Outer (power) iterations consumed.  ``None`` for fixed-source."""
        return self.record.n_iterations if self._is_outer else None

    def dominance_ratio(self) -> float | None:
        r"""Return :math:`|k_n - k_{n-1}| / |k_{n-1}|` in the late-iteration limit.

        The dominance ratio approximates the spectral gap
        :math:`|k_1 / k_0|` that controls power-iteration convergence.
        Returns ``None`` when fewer than 2 entries are recorded — the
        ratio is undefined for a single-point history.
        """
        if len(self.keff_history) < 2:
            return None
        prev = self.keff_history[-2]
        if prev == 0.0:
            return None
        return abs(self.keff_history[-1] - prev) / abs(prev)

    def latest_keff(self) -> float | None:
        """Return the last eigenvalue, or ``None`` for empty history."""
        return self.keff_history[-1] if self.keff_history else None

    def latest_residual(self) -> float | None:
        """Return the last recorded flux residual, or ``None`` if absent."""
        return self.flux_residuals[-1] if self.flux_residuals else None


@dataclass(frozen=True)
class SolutionBase:
    r"""Shared carrier for SN transport solutions — the role-agnostic base.

    Bundles the typed flux fields + boundary state + eigenvalue (when
    eigenvalue problem) + iteration trajectory into one frozen
    dataclass.  Both solution ROLES share this carrier; the role is the
    concrete type (see the module docstring's "Two discrimination
    axes"):

    * :class:`Solution` — the FORWARD solve (:func:`solve_sn` /
      :func:`solve_sn_fixed_source`); :attr:`scalar_flux` is the scalar
      flux :math:`\phi`.  Carries the forward-physics operations
      (:meth:`~Solution.homogenize`, :meth:`~Solution.condense`,
      :meth:`~Solution.reaction_rate_density`).
    * :class:`AdjointSolution` — the DAGGERED solve
      (:func:`solve_sn_adjoint` / :func:`solve_sn_adjoint_fixed_source`);
      :attr:`scalar_flux` is the importance :math:`\varphi^*`.

    ``SolutionBase`` itself is NOT instantiable: the role set is closed
    ({forward, adjoint}) and a role-less solution is not a value that
    exists (:meth:`__post_init__` raises ``TypeError`` on the base).
    Role-agnostic consumers (convergence diagnostics, plotting, the
    :meth:`compare` regression check) type against the base; consumers
    that read the physics type against the leaf they mean.

    Parameters
    ----------
    angular_flux : TimedFullField
        Composite iteration state — pure-Field
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        bulk paired with pure-Field
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
        boundary trace plus time-derivative history.

        Under D-H.1b (2026-05-28), this field carries a
        :class:`~orpheus.transport.timed_full_field.TimedFullField`
        composite. The bulk (``angular_flux.interior``) is the
        per-ordinate angular flux :math:`\psi(\vec r, \hat\Omega_n, g)`
        on shape ``(N, ng, nx, ny)``; the boundary trace is exposed
        via the :attr:`SolutionBase.boundary_flux` delegate property.
    scalar_flux : ScalarFlux
        Scalar flux field
        :math:`\phi(\vec r, g) = \sum_n w_n \psi_n`,
        shape ``(ng, nx, ny)``.
    mesh : SNMesh
        The phase-space mesh.  Both typed flux fields above MUST carry
        this same :class:`SNMesh` instance (validated at construction).
    keff : float or None
        Multiplication eigenvalue.  ``None`` for fixed-source problems.
    history : IterationHistory or None
        Convergence-trajectory diagnostics.  ``None`` for paths that
        do not surface them.
    radial_characteristic : RadialCharacteristicField or None
        System B's converged ψ½ state (B.2d DP-Solution — its OWN typed
        member, never a block on :attr:`angular_flux`): the marched
        starting-direction flux composite on a carrying mesh (R12a — the
        sphere), ``None`` exactly when the mesh carries no seed level
        (presence is structural — validated as a biconditional at
        construction). Downstream System-A readers (:attr:`scalar_flux`,
        :attr:`boundary_flux`, :class:`SolutionDiff`) are untouched; a
        ray reader reads THIS member.

    Notes
    -----
    Pre-D-H, :attr:`angular_flux` was the legacy
    :class:`orpheus.transport.fields.angular_flux.AngularFlux` (a bulk Field that
    ALSO owned the boundary face state and iteration history). Under
    D-H.1b that conflation dissolves: the composite-state container
    :class:`TimedFullField` holds the bulk + boundary + history trio
    as a structured composite. The :attr:`SolutionBase.boundary_flux`
    delegate (line below) becomes a thin read-through to
    ``self.angular_flux.boundary`` (now a typed L2 AngularBoundaryFlux).
    """

    angular_flux: "TimedFullField"
    scalar_flux: "ScalarFlux"
    mesh: "SNMesh"
    keff: float | None = None
    history: IterationHistory | None = None
    radial_characteristic: "RadialCharacteristicField | None" = None

    def __post_init__(self) -> None:
        # Role closure: the role set is closed ({forward, adjoint} —
        # Solution / AdjointSolution) and a role-less carrier is not a
        # value that exists (Pattern 4).  The leaves inherit this
        # __post_init__; the guard fires only on the base itself.
        if type(self) is SolutionBase:
            raise TypeError(
                "SolutionBase is the role-agnostic carrier base and is "
                "not instantiable — construct Solution (forward) or "
                "AdjointSolution (adjoint)."
            )
        # Mesh-binding consistency: every typed field MUST share the
        # exact SNMesh instance (Pattern 4 — illegal states
        # unrepresentable; the typed field's invariants are leveraged
        # at Solution construction time so consumers downstream cannot
        # see a Solution with mismatched meshes).
        # D-H.1b: angular_flux is now a TimedFullField composite —
        # the mesh is on the bulk (and validated against boundary at
        # TimedFullField construction).
        if self.angular_flux.interior.mesh is not self.mesh:
            raise ValueError(
                f"{type(self).__name__}: angular_flux.interior.mesh is not "
                f"{type(self).__name__}.mesh (typed-field mesh-identity "
                "contract broken — every field must reference the same "
                "SNMesh instance)."
            )
        if self.scalar_flux.mesh is not self.mesh:
            raise ValueError(
                f"{type(self).__name__}: scalar_flux.mesh is not "
                f"{type(self).__name__}.mesh (typed-field mesh-identity "
                "contract broken — every field must reference the same "
                "SNMesh instance)."
            )
        # B.2d DP-Solution: System B's presence is STRUCTURAL — the member
        # exists exactly when the mesh carries seed levels (R12a). A
        # carrying-mesh Solution without its converged ψ½ state (or a
        # seedless one carrying a ray) is a wiring error, not a variant.
        carries = (
            getattr(self.mesh, "radial_characteristic_field_space", None) is not None
        )
        if carries != (self.radial_characteristic is not None):
            raise ValueError(
                f"{type(self).__name__}: radial_characteristic presence "
                f"must match the mesh's R12a predicate (mesh carries: "
                f"{carries}, member present: "
                f"{self.radial_characteristic is not None}) — System B's "
                "converged state is its own typed member on a carrying "
                "mesh, absent otherwise (B.2d)."
            )
        if (
            self.radial_characteristic is not None
            and self.radial_characteristic.mesh is not self.mesh
        ):
            raise ValueError(
                f"{type(self).__name__}: radial_characteristic.mesh is not "
                f"{type(self).__name__}.mesh (typed-field mesh-identity "
                "contract broken)."
            )

    # ── boundary_flux as a delegate property ─────────────────────────
    #
    # The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # owns the boundary face state via :attr:`TimedFullField.boundary`.
    # ``Solution.boundary_flux`` is a thin read-through that delegates
    # to the composite's owned boundary trace.

    @property
    def boundary_flux(self) -> "AngularBoundaryFlux":
        r"""Boundary face state — delegate to :attr:`TimedFullField.boundary`.

        The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
        owns the boundary face state via :attr:`TimedFullField.boundary`.
        This property provides the canonical ``sol.boundary_flux``
        access path.

        See :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
        for the per-geometry flat-layout contract.
        """
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

        # Role parse at the composite boundary: a Solution's composite IS a
        # flux composite, but the ``FullField.boundary`` slot erases the
        # role (the F2-sibling erasure — #289). A source-role trace here
        # means the flux-composite contract broke upstream — raise loudly.
        boundary = self.angular_flux.boundary
        if not isinstance(boundary, AngularBoundaryFlux):
            raise TypeError(
                f"Solution.boundary_flux: the solution composite carries "
                f"{type(boundary).__name__}, not AngularBoundaryFlux — the "
                f"flux-composite contract is broken."
            )
        return boundary

    # ── Discrimination ───────────────────────────────────────────────

    def is_eigenvalue(self) -> bool:
        """Return ``True`` when this is an eigenvalue-problem solution."""
        return self.keff is not None

    def is_fixed_source(self) -> bool:
        """Return ``True`` when this is a fixed-source-problem solution."""
        return self.keff is None

    # ── Convergence diagnostics (delegate to history) ────────────────

    def dominance_ratio(self) -> float | None:
        r"""Return the eigenvalue dominance ratio, or ``None`` if unavailable."""
        return self.history.dominance_ratio() if self.history else None

    def converged(self) -> bool:
        """Return whether the solver iteration met its convergence tolerance.

        ⚠ A solution with **no history at all** answers ``False``.  This
        used to answer ``True`` — the only production branch that asserted
        convergence rather than reading it, and the same optimistic-default
        defect as #342.  "Nobody recorded whether this converged" is not
        evidence that it did; a caller gating on this method should treat an
        unrecorded solve exactly as it treats a truncated one.
        """
        return self.history.converged if self.history else False

    def keff_history_list(self) -> list[float]:
        """Return the eigenvalue trajectory as a list (for plotting)."""
        return list(self.history.keff_history) if self.history else []

    @property
    def keff_history(self) -> list[float]:
        """Legacy alias for :meth:`keff_history_list`.

        The canonical access path is :attr:`history` + ``.keff_history``;
        this property keeps existing test fixtures (``len(result.keff_history)``,
        plotting code) working without a one-by-one migration of ~10 call
        sites.  The trajectory is a list, not a tuple, so ``len`` and
        slice operations behave as they did under the legacy
        ``SNResult.keff_history: list[float]``.
        """
        return self.keff_history_list()

    # ── Comparison (role-closed: Self vs Self) ──────────────────────

    def compare(self, other: Self, *, rtol: float = 1e-12) -> "SolutionDiff":
        r"""Return a field-by-field difference summary against ``other``.

        Compares :attr:`keff` (when both have it) and the
        :math:`L^\infty`-norms of the flux deltas.  Useful for
        regression / refactor consistency checks.

        Parameters
        ----------
        other : Self
            The reference solution to compare against — of the SAME
            role (``Self``-typed: a forward flux and an importance map
            are different physical quantities, so cross-role comparison
            is a type error statically and a ``TypeError`` at runtime).
        rtol : float
            Relative tolerance for the ``within_tolerance`` flag.

        Returns
        -------
        SolutionDiff
            Field-by-field summary.
        """
        if type(other) is not type(self):
            raise TypeError(
                f"{type(self).__name__}.compare: role mismatch — comparing "
                f"against {type(other).__name__}.  A forward flux and an "
                "importance map are different physical quantities; "
                "same-role comparison only (the Self-typed contract, "
                "enforced at runtime for untyped callers)."
            )
        if not self.mesh.is_same_phase_space(other.mesh):
            raise ValueError(
                f"{type(self).__name__}.compare: the solutions realize "
                "different discrete phase spaces — comparison is defined "
                "only across solves sharing the same constituents (see "
                "SNMesh.is_same_phase_space)."
            )

        if self.keff is not None and other.keff is not None:
            keff_abs = abs(self.keff - other.keff)
        else:
            keff_abs = None

        ang_diff = self.angular_flux.interior.values - other.angular_flux.interior.values
        sca_diff = self.scalar_flux.values - other.scalar_flux.values
        ang_linf = float(np.abs(ang_diff).max()) if ang_diff.size else 0.0
        sca_linf = float(np.abs(sca_diff).max()) if sca_diff.size else 0.0

        sca_norm = float(np.abs(other.scalar_flux.values).max())
        flux_ok = (sca_norm == 0.0) or (sca_linf <= rtol * sca_norm)
        keff_ok = (keff_abs is None) or (
            abs(other.keff or 1.0) > 0.0
            and keff_abs <= rtol * abs(other.keff or 1.0)
        )

        return SolutionDiff(
            keff_abs=keff_abs,
            angular_flux_linf=ang_linf,
            scalar_flux_linf=sca_linf,
            within_tolerance=bool(flux_ok and keff_ok),
        )


@dataclass(frozen=True)
class Solution(SolutionBase):
    r"""Canonical return type for the FORWARD solvers.

    The forward role of the :class:`SolutionBase` carrier — what
    :func:`solve_sn` / :func:`solve_sn_fixed_source` return:
    :attr:`~SolutionBase.scalar_flux` is the scalar flux
    :math:`\phi(\vec r, g) = \sum_n w_n \psi_n`, and the forward-physics
    operations live HERE and only here:

    * :meth:`reaction_rate_density` — :math:`\sigma \cdot \phi`;
    * :meth:`homogenize` — reaction-rate-preserving spatial collapse;
    * :meth:`condense` — reaction-rate-preserving energy collapse.

    All three interpret ``scalar_flux`` as the flux; none exists on
    :class:`AdjointSolution` (structural asymmetry — see the module
    docstring).  One type covers both PROBLEM KINDS (fixed-source and
    eigenvalue) via the optional :attr:`~SolutionBase.keff` /
    :attr:`~SolutionBase.history` (#197 PR-TYPED-5; replaces the legacy
    ``SNFixedSourceResult`` / ``SNResult`` pair).

    The adjoint-weighted refinement of :meth:`homogenize` /
    :meth:`condense` is the ratified #281 (P6-B2) API, **landed**: an
    optional keyword ``adjoint: AdjointSolution | None = None`` — ``None``
    keeps today's flux-weighted (Galerkin, :math:`\varphi^* = \phi`
    degenerate) collapse bit-identically; a real importance makes the
    collapse eigenvalue-consistent per the worth-zeroing taxonomy of
    :mod:`orpheus.derivations.common.homogenization` (spatial T1/T1b/T2/T3;
    energy = the B&G-convention bilinear, T6).

    Examples
    --------
    Reads as the math:

    >>> sol = solve_sn(materials, mesh, quadrature)             # doctest: +SKIP
    >>> sol.is_eigenvalue()                                     # doctest: +SKIP
    True
    >>> sol.dominance_ratio()                                   # doctest: +SKIP
    1.2e-08
    >>> sol.reaction_rate_density(materials[0].sig_a)           # doctest: +SKIP
    array(...)
    """

    # ── Reaction-rate accessor (Pattern 1 — read as math) ─────────────

    def reaction_rate_density(self, xs: np.ndarray) -> np.ndarray:
        r"""Compute the per-cell reaction-rate density :math:`\sigma \cdot \phi`.

        Parameters
        ----------
        xs : np.ndarray
            Cross-section array, shape ``(ng, *spatial)`` (a
            :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`-shaped
            ndarray).  Per Issue #197 Wave 1 design, ``ReactionRate`` is
            a ``NewType`` over ``np.ndarray`` and is NOT promoted to a
            dataclass — staying close to the natural ``ng × *spatial``
            layout the consumers (depletion, k-eff diagnostics, response
            functionals) want directly.

        Returns
        -------
        np.ndarray
            Per-cell rate density, shape ``(ng, *spatial)``.  Units:
            :math:`\rm 1/(cm^3 \cdot s)` — the reaction-rate density
            per group.
        """
        return xs * self.scalar_flux.values

    # ── Spatial homogenization (a domain operation on the solution) ──

    def homogenize(
        self,
        coarse: "Mesh1D | Mesh2D",
        *,
        adjoint: "AdjointSolution | None" = None,
    ) -> "MaterialMesh":
        r"""Flux·volume-weighted spatial homogenization onto a coarse mesh.

        Collapse this (fine) solution's per-cell cross sections onto the
        cells of a coarser mesh, preserving every reaction rate.  Each
        coarse cell :math:`R` becomes its own effective material whose
        cross sections are the flux·volume-weighted averages of the fine
        cells it contains:

        .. math::

            \Sigma_{R,g} \;=\;
            \frac{\sum_{i \in R} V_i\,\phi_{i,g}\,\Sigma_{i,g}}
                 {\sum_{i \in R} V_i\,\phi_{i,g}}
            \qquad
            \Phi_{R,g} \;=\; \sum_{i \in R} V_i\,\phi_{i,g}

        so that the volume-integrated reaction rate is preserved exactly:
        :math:`\Sigma_{R,g}\,\Phi_{R,g} = \sum_{i \in R} V_i\,\Sigma_{i,g}\,
        \phi_{i,g}`.  This is the **space-only**, **mesh-COUPLED** half of
        the condense/homogenize asymmetry law: the coarse cells carry the
        homogenized materials, so geometry + materials are born together
        (returned as a :class:`MaterialMesh`).  Energy is **not** condensed
        — the group structure (``eg``) is carried through unchanged.

        The matrix channels (per-Legendre :math:`\Sigma_{s,\ell}[g',g]` and
        :math:`\Sigma_{2n}[g',g]`, indexed ``[g_from, g_to]``) weight by the
        **source** group :math:`g'` flux (the group whose population drives
        the out-scatter), and :math:`\chi` is the **production-weighted**
        convex average (weight :math:`p_i = \sum_g \nu\Sigma_{f,i,g}\,
        \phi_{i,g}\,V_i`) — a convex combination of simplices, hence a
        simplex (validated by :class:`Mixture.__post_init__`).  Because every
        removal channel collapses with the *same* per-(R, g) weight, the
        definitional total-XS balance :math:`\Sigma_t = \Sigma_c + \Sigma_L
        + \Sigma_f + \mathrm{rowsum}(\Sigma_{s0}) + \mathrm{rowsum}(\Sigma_{2n})`
        is preserved cell-by-cell when the fine materials balance.

        With ``adjoint=`` (P6, #281) the collapse becomes **eigenvalue-
        consistent**: every channel takes the worth-zeroing rule of the
        algebra of record (:mod:`orpheus.derivations.common.homogenization`),
        so the coarse :math:`k` is first-order stationary in the flux shapes —
        the vector channels take the bilinear

        .. math::

            \Sigma_{R,g} \;=\;
            \frac{\sum_{i \in R} V_i\,\varphi^*_{i,g}\,\Sigma_{i,g}\,\varphi_{i,g}}
                 {\sum_{i \in R} V_i\,\varphi^*_{i,g}\,\varphi_{i,g}}

        (:eq:`sn-homogenization-adjoint-weighted` — the test weight is the
        PRODUCT :math:`\varphi^*\!\odot\varphi`, the Petrov-Galerkin lift the
        forward call is the :math:`\varphi^*{=}\,\text{flat}` degenerate of);
        :math:`\Sigma_t` takes the exact ANGULAR pairing, the matrices the
        per-pair sink×source rule, and the fission dyad the mixed-fold
        factored rule — the full taxonomy and its theorems live on
        :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through_bilinear`.

        .. warning:: An adjoint-weighted (worth-exact) collapse **breaks the
           total-XS balance identity** (the classical reactivity-vs-rates
           property; theorem T4).  Do NOT ``assert_balanced`` on the returned
           materials when ``adjoint`` was given.

        Parameters
        ----------
        coarse : Mesh1D or Mesh2D
            The coarse target mesh (any dimension matching the fine mesh).
            Must share this solution's outer boundary on every axis; its
            internal cell edges must align with fine-cell edges (each coarse
            cell is a contiguous union of fine cells).  Its own ``mat_ids`` /
            ``mat_map`` are ignored — homogenization assigns one fresh
            effective material per coarse cell.
        adjoint : AdjointSolution, optional
            The importance solution :math:`\psi^*` from
            :func:`~orpheus.sn.solver.solve_sn_adjoint` on the SAME mesh
            object (identity-checked).  ``None`` (default) keeps the forward
            flux-weighted (Galerkin-degenerate) collapse, bit-identical to
            the pre-P6 behaviour.

        Returns
        -------
        MaterialMesh
            The coarse mesh carrying the homogenized materials (one
            :class:`Mixture` per coarse cell, keyed by coarse-cell index).
            Promote to a solvable SN phase space with
            :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.from_material_mesh`.

        Notes
        -----
        Dimension-agnostic: the coarse cell-indicator basis and the fine
        volume measure are n-D (the membership table and the group axis ride
        the frame's contractions), so 1-D and 2-D meshes flow through the
        same body.  A group with zero region flux (:math:`\Phi_{R,g} = 0`)
        yields a zero effective XS for that (R, g) — there is no reaction
        rate to preserve there (the frame's Moore–Penrose Gram pseudo-inverse).
        """
        from orpheus.numerics.basis import WeightedIndicatorBasis
        from orpheus.numerics.frame import PetrovGalerkinFrame
        from orpheus.transport.mesh.material_mesh import MaterialMesh
        from orpheus.transport.mesh.material_xs_field import MaterialXSField

        fine = self.mesh
        trial = coarse.indicator_basis()       # coarse cell-indicator trial basis (n-D)
        if trial.ndim != fine.ndim:
            raise ValueError(
                f"Solution.homogenize: coarse mesh dimension {trial.ndim} must "
                f"match the fine mesh dimension {fine.ndim}."
            )

        # Each coarse cell is a contiguous union of fine cells, so the coarse mesh
        # must share the fine mesh's outer boundary on every axis.
        for axis in range(fine.ndim):
            fine_edges = np.asarray(fine.axes[axis].edges, dtype=float)
            coarse_edges = np.asarray(trial.edges_per_axis[axis], dtype=float)
            if not (
                np.isclose(coarse_edges[0], fine_edges[0])
                and np.isclose(coarse_edges[-1], fine_edges[-1])
            ):
                raise ValueError(
                    "Solution.homogenize: coarse mesh must share the fine mesh's "
                    f"outer boundary on axis {axis} "
                    f"[{fine_edges[0]}, {fine_edges[-1]}]; got "
                    f"[{coarse_edges[0]}, {coarse_edges[-1]}]."
                )

        measure = fine.volume_measure
        ng = fine.ng
        # (ng, *spatial) → (n_fine, ng) in the "ij"/C flat-cell order the measure
        # nodes and ``mat_map.ravel()`` share (1-D: a plain transpose; n-D: the
        # spatial axes collapse to one fine-cell axis).
        phi = np.asarray(self.scalar_flux.values, dtype=float).reshape(ng, -1).T

        if adjoint is None:
            # The two forward Petrov-Galerkin homogenisation frames. Trial = the
            # coarse cell indicators (the mesh YIELDS them); measure = the fine
            # geometric volume measure dV. The solution-weighting rides the TEST
            # side (the frame TYPE), NEVER folded into the measure: the measure
            # carries the axis + the fixed L² metric, the flux is a test-weighting
            # the solution emits. The two collapses preserve two different
            # conserved rates — Σ preserves the reaction rate (flux-weighted test
            # φ·1_R), χ preserves the emission rate (production-weighted test
            # p·1_R) — so each is its own frame. Both are the flat-φ* degenerate
            # of the eigenvalue-consistent adjoint-weighted (φ*≠φ) collapse below;
            # ``project`` is G⁻¹M with a diagonal (disjoint-indicator) Gram, whose
            # Moore–Penrose pseudo-inverse zeroes empty / zero-flux regions.
            sigma_frame = PetrovGalerkinFrame(
                trial, measure, WeightedIndicatorBasis(trial, phi),
            )

            mat_of_fine = np.asarray(fine.mat_map, dtype=int).ravel()  # (n_fine,)
            nu_sigma_f = np.array(
                [fine.materials[m].SigP for m in mat_of_fine]          # (n_fine, ng)
            )
            production = (nu_sigma_f * phi).sum(axis=1)                # p_i, (n_fine,)
            emission_frame = PetrovGalerkinFrame(
                trial, measure, WeightedIndicatorBasis(trial, production),
            )

            # Project the WHOLE cross-section field as one object: the field owns
            # the channel → weighting taxonomy and routes Σ → sigma_frame,
            # χ → emission_frame.
            homogenized = MaterialXSField.from_mesh(fine).project_through(
                sigma_frame, emission_frame,
            )
        else:
            # The eigenvalue-consistent arm (P6 #281). The role is a TYPE: only
            # an AdjointSolution can weight the test side (a forward Solution
            # here would silently compute the wrong physics), and the mesh must
            # be the SAME OBJECT (the σ↔geometry pairing tier — identity
            # guarantees shape AND the shared "ij" flat order).
            if not isinstance(adjoint, AdjointSolution):
                raise TypeError(
                    f"Solution.homogenize: adjoint must be an AdjointSolution "
                    f"(the importance is the test weight, never a forward flux); "
                    f"got {type(adjoint).__name__}."
                )
            if not fine.is_same_phase_space(adjoint.mesh):
                raise ValueError(
                    "Solution.homogenize: adjoint solves a different discrete "
                    "phase space — the importance must come from an adjoint "
                    "solve sharing this solution's constituents (same geometry "
                    "mesh, quadrature, and materials OBJECTS, same scheme; see "
                    "SNMesh.is_same_phase_space)."
                )
            phi_star = np.asarray(
                adjoint.scalar_flux.values, dtype=float,
            ).reshape(ng, -1).T
            # The exact collision pairing ρ_{i,g} = Σ_n w_n ψ*ψ (T1b — the
            # user-ruled angular rule; both solutions carry ψ).
            w = np.asarray(fine.quad.weights, dtype=float)             # (N,)
            psi = np.asarray(self.angular_flux.interior.values, dtype=float)
            psi_star = np.asarray(adjoint.angular_flux.interior.values, dtype=float)
            rho = np.einsum("n,n...->...", w, psi_star * psi)          # (ng, *spatial)
            rho = rho.reshape(ng, -1).T                                # (n_fine, ng)

            # The field owns the five-morphism bilinear taxonomy (T1/T1b/T2/T3).
            homogenized = MaterialXSField.from_mesh(fine).project_through_bilinear(
                trial, measure, phi=phi, phi_star=phi_star, rho=rho,
            )

        # The coarse geometry with each cell relabelled to its own material id
        # (dimension-agnostic — no Mesh1D/Mesh2D reconstruction branch), carrying the
        # one fresh effective material per coarse cell.
        return MaterialMesh(coarse.with_distinct_cell_ids(), homogenized)

    # ── Energy condensation (the energy-axis transpose of homogenize) ──

    def condense(
        self,
        coarse: "EnergyGrid",
        *,
        adjoint: "AdjointSolution | None" = None,
        within_group: "WithinGroupSpectrum | None" = None,
    ) -> dict[int, "Mixture"]:
        r"""Spectrum-weighted energy condensation onto a coarse group structure.

        Collapse this solution's per-material cross sections from the fine
        (solved) group structure onto a coarser :class:`~orpheus.data.energy_grid.EnergyGrid`,
        preserving every reaction rate.  Each material is condensed with its own
        **representative spectrum** — the flux·volume-weighted flux over the cells
        where the material appears:

        .. math::

            \varphi^{(m)}_g \;=\; \sum_{i:\,\mathrm{mat}(i)=m} V_i\,\phi_{i,g}

        used as the test weight in
        :meth:`orpheus.data.macro_xs.mixture.Mixture.condense` (which preserves the
        per-coarse-group reaction rate).  This is the **energy-only**,
        **mesh-DECOUPLED** half of the condense/homogenize asymmetry law: the
        result is **portable** few-group cross sections
        (``dict[material_id, Mixture]``), NOT bound to any mesh — geometry is
        untouched.  (Contrast :meth:`homogenize`, which collapses space and returns
        a mesh-COUPLED :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`.)

        Parameters
        ----------
        coarse : EnergyGrid
            The coarse target group structure (descending boundaries; no more
            groups than the fine structure — condensation only downsamples, see
            the upscaling guard in
            :meth:`~orpheus.data.energy_grid.EnergyGrid.overlap_to`).
        adjoint : AdjointSolution, optional
            The importance solution from
            :func:`~orpheus.sn.solver.solve_sn_adjoint` on the same discrete
            phase space (guarded via
            :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.is_same_phase_space`).
            ``None`` (default) keeps the flux-weighted collapse, bit-identical.
            Given, each material condenses with its representative SPECTRUM
            PAIR :math:`(\varphi^{(m)}, \varphi^{*(m)})` (the same flux·volume
            reduction applied to both solutions) through the **bilinear
            (eigenvalue-consistent) B&G-convention collapse** — see
            :meth:`Mixture.condense <orpheus.data.macro_xs.mixture.Mixture.condense>`
            (``adjoint_spectrum=``) and theorem T6 of
            :mod:`orpheus.derivations.common.homogenization`.  The bilinear
            constants do NOT satisfy the total-XS balance identity (the
            classical reactivity-vs-rates trade-off, T4) — do not
            ``assert_balanced`` on them.
        within_group : WithinGroupSpectrum, optional
            The sub-fine-group flux model for straddle apportionment, threaded
            through to :meth:`Mixture.condense` (default 1/E).

        Returns
        -------
        dict[int, Mixture]
            One condensed :class:`~orpheus.data.macro_xs.mixture.Mixture` per
            material id, carrying the coarse ``eg`` — portable few-group XS (e.g.
            to compare against a WIMS few-group library, or seed a coarse solve).

        Notes
        -----
        The fine→coarse :class:`~orpheus.numerics.basis.OverlapBasis` (via
        :meth:`~orpheus.data.energy_grid.EnergyGrid.overlap_to`) is built per material
        from its ``eg`` (the partition is identical across a uniform-grid mesh; only the
        spectrum differs).  A material with no flux in
        a fine group contributes zero weight there; the condense frame's
        Moore–Penrose Gram handles any empty coarse group.

        Raises
        ------
        ValueError
            If a material carries no energy grid (``eg is None`` — a synthetic
            mixture); condensation needs the fine library grid.
        """
        fine = self.mesh
        ng = fine.ng
        # (ng, *spatial) → (n_fine_cells, ng) in the "ij"/C flat-cell order the
        # volume measure and ``mat_map.ravel()`` share (same convention as homogenize).
        phi = np.asarray(self.scalar_flux.values, dtype=float).reshape(ng, -1).T
        volume = np.asarray(fine.volume_measure.weights, dtype=float)   # (n_cells,)
        mat_of_cell = np.asarray(fine.mat_map, dtype=int).ravel()       # (n_cells,)

        phi_star = None
        if adjoint is not None:
            # The role is a TYPE and the problems must match — the same guard
            # pair as the homogenize adjoint arm. (Collapse into a shared
            # helper when a THIRD adjoint-consuming verb lands — perturbation
            # worth / GPT are the anticipated instances.)
            if not isinstance(adjoint, AdjointSolution):
                raise TypeError(
                    f"Solution.condense: adjoint must be an AdjointSolution "
                    f"(the importance is the test weight, never a forward "
                    f"flux); got {type(adjoint).__name__}."
                )
            if not fine.is_same_phase_space(adjoint.mesh):
                raise ValueError(
                    "Solution.condense: adjoint solves a different discrete "
                    "phase space — the importance must come from an adjoint "
                    "solve sharing this solution's constituents (see "
                    "SNMesh.is_same_phase_space)."
                )
            phi_star = np.asarray(
                adjoint.scalar_flux.values, dtype=float,
            ).reshape(ng, -1).T

        condensed: dict[int, "Mixture"] = {}
        for mat_id, material in fine.materials.items():
            if material.eg is None:
                raise ValueError(
                    f"Solution.condense: material {mat_id} has no energy grid "
                    f"(eg is None); condensation needs the fine library grid."
                )
            cells = mat_of_cell == mat_id

            def _representative(field: np.ndarray) -> np.ndarray:
                # ONE reduction for the pair — the T6 carrier consistency
                # (B&G convention) requires φ and φ* be reduced by the SAME
                # operator; naming it once makes that structural.
                return (volume[cells, None] * field[cells]).sum(axis=0)

            spectrum = _representative(phi)                             # (ng,)
            if phi_star is None:
                condensed[mat_id] = material.condense(coarse, spectrum, within_group)
            else:
                condensed[mat_id] = material.condense(
                    coarse, spectrum, within_group,
                    adjoint_spectrum=_representative(phi_star),
                )
        return condensed


@dataclass(frozen=True)
class AdjointSolution(SolutionBase):
    r"""Canonical return type for the ADJOINT solvers.

    The adjoint role of the :class:`SolutionBase` carrier — what
    :func:`solve_sn_adjoint` / :func:`solve_sn_adjoint_fixed_source`
    return (#276 A4, the daggered posing):

    * :attr:`~SolutionBase.angular_flux` carries :math:`\psi^*` — the
      converged state of the daggered system: the exact DISCRETE
      transpose in the solution G-metric (``A.H`` — #280's swap-law
      adjoint), NOT a :math:`\mu`-reversed forward flux.
    * :attr:`~SolutionBase.scalar_flux` carries the **importance**
      :math:`\varphi^*(\vec r, g) = \sum_n w_n \psi^*_n` — the same
      :math:`w`-reduction as the forward scalar flux (the adjoint of
      the ISO source injection, not a new functional).
      :attr:`importance` is the domain-named alias.
    * :attr:`~SolutionBase.keff` is the eigenvalue of the daggered
      pencil :math:`(A^\dagger, F^\dagger)` — EXACTLY the forward
      :math:`k` in exact arithmetic (:math:`\operatorname{eig}(M^T) =
      \operatorname{eig}(M)`); the two power iterations agree to
      iteration tolerance (the P1.3 certification rows).

    Deliberately ABSENT: :meth:`Solution.homogenize`,
    :meth:`Solution.condense`, :meth:`Solution.reaction_rate_density`.
    Those interpret ``scalar_flux`` as the flux :math:`\phi` and
    preserve reaction rates — an importance map has no reaction rate to
    preserve.  The adjoint enters homogenization/condensation as the
    OPTIONAL TEST WEIGHT of the FORWARD collapse
    (``Solution.homogenize(..., adjoint=...)`` — the #281 P6-B2
    parameter, landed), never as its subject.  The absence is structural
    (no attribute exists), so the wrong physics cannot be spelled
    (``coding-elegance`` Pattern 4).
    """

    @property
    def importance(self) -> "ScalarFlux":
        r"""The importance map — the domain name for the adjoint scalar flux.

        :math:`\varphi^*(\vec r, g)`: the expected detector response
        per unit source particle introduced at :math:`(\vec r, g)` (the
        classical importance interpretation of the adjoint flux).  An
        alias of :attr:`~SolutionBase.scalar_flux` — one storage, two
        vocabularies.
        """
        return self.scalar_flux


@dataclass(frozen=True)
class SolutionDiff:
    r"""Result of :meth:`SolutionBase.compare` (both roles).

    Parameters
    ----------
    keff_abs : float or None
        Absolute eigenvalue difference :math:`|k_a - k_b|`, or
        ``None`` when at least one solution is fixed-source.
    angular_flux_linf : float
        :math:`L^\infty` norm of the angular-flux delta.
    scalar_flux_linf : float
        :math:`L^\infty` norm of the scalar-flux delta.
    within_tolerance : bool
        Aggregate verdict: True iff every available channel met the
        comparison rtol.
    """

    keff_abs: float | None
    angular_flux_linf: float
    scalar_flux_linf: float
    within_tolerance: bool
