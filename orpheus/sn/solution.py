r"""Typed return type for the SN transport solvers.

Issue #197 PR-TYPED-5 — :class:`Solution` and :class:`IterationHistory`
replace the legacy bare-dataclass pair
:class:`~orpheus.sn.solver.SNFixedSourceResult` /
:class:`~orpheus.sn.solver.SNResult`.

The legacy types were data bags: they carried bare ``np.ndarray`` flux
fields, an opaque ``geometry`` / ``quadrature`` pair, and an ad-hoc set
of diagnostic scalars (``n_inner``, ``residual``, ``keff_history``,
``elapsed_seconds``).  Two distinct dataclasses for the two problem
kinds (fixed-source vs eigenvalue) duplicated the shape contract in two
places — a twin path waiting for drift (``coding-elegance`` anti-pattern
1: "two implementations of the same mathematical quantity").

Under PR-TYPED-5, **one** typed :class:`Solution` covers both problem
kinds.  The discrimination lives in two methods —
:meth:`Solution.is_eigenvalue` and :meth:`Solution.is_fixed_source` —
that read the optional :attr:`keff`.  The convergence-trajectory
diagnostics live on a separate :class:`IterationHistory`, exposed
through method-style accessors :meth:`Solution.dominance_ratio` and
:meth:`Solution.converged`.

Reads as the math (``coding-elegance`` Pattern 1 — match the algebra of
the domain)::

    sol = solve_sn(...)
    sol.is_eigenvalue()                    # True / False
    sol.dominance_ratio()                  # |k_n / k_{n-1}|
    sol.reaction_rate_density(sig_a)       # σ · φ
    sol.compare(other, rtol=1e-12)         # SolutionDiff

Mesh-binding consistency
========================

:class:`Solution.__post_init__` validates that every typed field
(:attr:`angular_flux`, :attr:`scalar_flux`, :attr:`boundary_flux`)
carries the SAME :class:`SNMesh` instance as the
:attr:`Solution.mesh`.  Mesh identity is checked via ``is`` — sharing a
mesh by value rather than by reference is forbidden by construction
(``coding-elegance`` Pattern 4 — illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.timed_full_field import TimedFullField


__all__ = ["IterationHistory", "Solution", "SolutionDiff"]


@dataclass(frozen=True)
class IterationHistory:
    r"""Convergence trajectory diagnostics for an SN solve.

    Carries the per-outer / per-inner trajectory exposed by the legacy
    :class:`SNResult.keff_history` (list) plus the integer
    :attr:`SNFixedSourceResult.n_inner` + :attr:`residual` scalars — now
    bundled into one frozen container with method-style diagnostics.

    Parameters
    ----------
    keff_history : tuple of float
        Per-outer-iteration eigenvalue trajectory.  Empty tuple for
        fixed-source problems (no eigenvalue to track).
    flux_residuals : tuple of float
        Per-iteration relative-flux residual trajectory
        :math:`\lVert\phi_n - \phi_{n-1}\rVert / \lVert\phi_n\rVert`.
        Empty when the solver does not track this signal.
    n_inner : int or None
        Number of inner (within-group) iterations consumed.  ``None``
        for paths that do not surface this count (e.g. pure eigenvalue
        outer iteration).
    total_inner_iterations : int or None
        Total inner (within-group) iterations summed across ALL outer
        iterations.  Unlike :attr:`n_inner` (which is ``None`` on the
        eigenvalue outer path — there is no single inner solve to point
        at), this is populated by BOTH paths: the eigenvalue path
        accumulates it across the power-iteration outer loop, the
        fixed-source path sets it equal to :attr:`n_inner` (one inner
        solve).  It is the measurand for the SI spectral-rate /
        Gauss-Seidel-recovery diagnostics (Phase 3).  ``None`` only when
        no path surfaced a count.
    n_outer : int or None
        Number of outer (power) iterations consumed.  ``None`` for
        fixed-source problems.
    converged : bool
        ``True`` when the iteration met its convergence tolerance,
        ``False`` when ``max_iter`` was exhausted.

    Notes
    -----
    The trajectories are **tuple** rather than list — frozen dataclasses
    with mutable fields are an anti-pattern (``coding-elegance``
    Pattern 4: illegal states unrepresentable).  Callers that want a
    list pass through :func:`list` at the call site.
    """

    keff_history: tuple[float, ...] = ()
    flux_residuals: tuple[float, ...] = ()
    n_inner: int | None = None
    total_inner_iterations: int | None = None
    n_outer: int | None = None
    converged: bool = True

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
class Solution:
    r"""Canonical return type for :func:`solve_sn` / :func:`solve_sn_fixed_source`.

    Bundles the typed flux fields + boundary state + eigenvalue (when
    eigenvalue problem) + iteration trajectory into one frozen
    dataclass.  Replaces :class:`SNFixedSourceResult` + :class:`SNResult`
    (which were data bags without methods).  One type covers both
    fixed-source and eigenvalue problems via optional :attr:`keff` and
    :attr:`history`.

    Parameters
    ----------
    angular_flux : TimedFullField
        Composite iteration state — pure-Field
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        bulk paired with pure-Field
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        boundary trace plus time-derivative history.

        Under D-H.1b (2026-05-28), this field carries a
        :class:`~orpheus.transport.timed_full_field.TimedFullField`
        composite. The bulk (``angular_flux.bulk``) is the
        per-ordinate angular flux :math:`\psi(\vec r, \hat\Omega_n, g)`
        on shape ``(N, ng, nx, ny)``; the boundary trace is exposed
        via the :attr:`Solution.boundary_flux` delegate property.
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

    Notes
    -----
    Pre-D-H, :attr:`angular_flux` was the legacy
    :class:`orpheus.sn.angular_flux.AngularFlux` (a bulk Field that
    ALSO owned the boundary face state and iteration history). Under
    D-H.1b that conflation dissolves: the composite-state container
    :class:`TimedFullField` holds the bulk + boundary + history trio
    as a structured composite. The :attr:`Solution.boundary_flux`
    delegate (line below) becomes a thin read-through to
    ``self.angular_flux.boundary`` (now a typed L2 BoundaryFlux).

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

    angular_flux: "TimedFullField"
    scalar_flux: "ScalarFlux"
    mesh: "SNMesh"
    keff: float | None = None
    history: IterationHistory | None = None

    def __post_init__(self) -> None:
        # Mesh-binding consistency: every typed field MUST share the
        # exact SNMesh instance (Pattern 4 — illegal states
        # unrepresentable; the typed field's invariants are leveraged
        # at Solution construction time so consumers downstream cannot
        # see a Solution with mismatched meshes).
        # D-H.1b: angular_flux is now a TimedFullField composite —
        # the mesh is on the bulk (and validated against boundary at
        # TimedFullField construction).
        if self.angular_flux.bulk.mesh is not self.mesh:
            raise ValueError(
                "Solution: angular_flux.bulk.mesh is not Solution.mesh "
                "(typed-field mesh-identity contract broken — every "
                "field must reference the same SNMesh instance)."
            )
        if self.scalar_flux.mesh is not self.mesh:
            raise ValueError(
                "Solution: scalar_flux.mesh is not Solution.mesh "
                "(typed-field mesh-identity contract broken — every "
                "field must reference the same SNMesh instance)."
            )

    # ── boundary_flux as a delegate property ─────────────────────────
    #
    # The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # owns the boundary face state via :attr:`TimedFullField.boundary`.
    # ``Solution.boundary_flux`` is a thin read-through that delegates
    # to the composite's owned boundary trace.

    @property
    def boundary_flux(self) -> "BoundaryFlux":
        r"""Boundary face state — delegate to :attr:`TimedFullField.boundary`.

        The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
        owns the boundary face state via :attr:`TimedFullField.boundary`.
        This property provides the canonical ``sol.boundary_flux``
        access path.

        See :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        for the per-geometry flat-layout contract.
        """
        return self.angular_flux.boundary

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
        """Return whether the solver iteration met its convergence tolerance."""
        return self.history.converged if self.history else True

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
        :class:`SNResult.keff_history: list[float]`.
        """
        return self.keff_history_list()

    # ── Reaction-rate accessor (Pattern 1 — read as math) ─────────────

    def reaction_rate_density(self, xs: np.ndarray) -> np.ndarray:
        r"""Compute the per-cell reaction-rate density :math:`\sigma \cdot \phi`.

        Parameters
        ----------
        xs : np.ndarray
            Cross-section array, shape ``(ng, nx, ny)`` (a
            :class:`~orpheus.sn.material_xs_field.MaterialXSField`-shaped
            ndarray).  Per Issue #197 Wave 1 design, ``ReactionRate`` is
            a ``NewType`` over ``np.ndarray`` and is NOT promoted to a
            dataclass — staying close to the natural ``ng × nx × ny``
            layout the consumers (depletion, k-eff diagnostics, response
            functionals) want directly.

        Returns
        -------
        np.ndarray
            Per-cell rate density, shape ``(ng, nx, ny)``.  Units:
            :math:`\rm 1/(cm^3 \cdot s)` — the reaction-rate density
            per group.
        """
        return np.einsum("gxy,gxy->gxy", xs, self.scalar_flux.values)

    # ── Comparison ──────────────────────────────────────────────────

    def compare(self, other: "Solution", *, rtol: float = 1e-12) -> "SolutionDiff":
        r"""Return a field-by-field difference summary against ``other``.

        Compares :attr:`keff` (when both have it) and the
        :math:`L^\infty`-norms of the flux deltas.  Useful for
        regression / refactor consistency checks.

        Parameters
        ----------
        other : Solution
            The reference solution to compare against.
        rtol : float
            Relative tolerance for the ``within_tolerance`` flag.

        Returns
        -------
        SolutionDiff
            Field-by-field summary.
        """
        if self.mesh is not other.mesh:
            raise ValueError(
                "Solution.compare: meshes differ — cross-mesh comparison "
                "is not defined under the typed-field contract."
            )

        if self.keff is not None and other.keff is not None:
            keff_abs = abs(self.keff - other.keff)
        else:
            keff_abs = None

        ang_diff = self.angular_flux.bulk.values - other.angular_flux.bulk.values
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
class SolutionDiff:
    r"""Result of :meth:`Solution.compare`.

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
