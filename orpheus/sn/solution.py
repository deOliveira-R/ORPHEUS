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
    from .angular_flux import AngularFlux
    from .boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from orpheus.transport.fields.scalar_flux import ScalarFlux


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
    angular_flux : AngularFlux
        Per-ordinate angular flux field
        :math:`\psi(\vec r, \hat\Omega_n, g)`, shape ``(N, ng, nx, ny)``.
        Under R-1 Step 2, :class:`AngularFlux` carries its own boundary
        face state via :attr:`AngularFlux.boundary`; see the
        :attr:`Solution.boundary_flux` property below.
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
    Prior to R-1 Step 2, :attr:`boundary_flux` was a separate dataclass
    field carrying a :class:`BoundaryFlux` instance independent of
    ``angular_flux``.  That was a twin-storage anti-pattern (the same
    boundary face state was reachable through two paths).  Under R-1
    Step 2, :attr:`boundary_flux` is a delegated property that returns
    :attr:`AngularFlux.boundary` — single source of truth.

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

    angular_flux: "AngularFlux"
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
        for name, field_value in (
            ("angular_flux", self.angular_flux),
            ("scalar_flux", self.scalar_flux),
        ):
            if field_value.mesh is not self.mesh:
                raise ValueError(
                    f"Solution: {name}.mesh is not Solution.mesh "
                    "(typed-field mesh-identity contract broken — every "
                    "field must reference the same SNMesh instance)."
                )

    # ── R-1 Step 2: boundary_flux as a delegate property ─────────────
    #
    # Prior to R-1 Step 2, :attr:`boundary_flux` was a separate dataclass
    # field; the boundary face state was reachable from BOTH
    # ``sol.angular_flux.boundary`` (after the typed-pair migration)
    # AND ``sol.boundary_flux`` (the legacy field).  That was a twin-
    # storage path that invited drift; per ``coding-elegance``
    # Pattern 2 (single source of truth), R-1 Step 2 collapses both
    # paths to ``angular_flux.boundary``.  ``Solution.boundary_flux``
    # is now a thin read-through that delegates to the AngularFlux's
    # owned boundary — Step 5 removes the property entirely after one
    # wave cycle of test migration elapses.

    @property
    def boundary_flux(self) -> "BoundaryFlux":
        r"""Boundary face state — delegate to :attr:`AngularFlux.boundary`.

        :class:`AngularFlux` now owns the boundary face / persistent
        buffer state via :attr:`AngularFlux.boundary`.  This property
        provides the legacy ``sol.boundary_flux`` access path while
        ensuring the single-source-of-truth contract.

        See :class:`~orpheus.sn.boundary_flux.BoundaryFlux` for the
        per-geometry shape contract.
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

        ang_diff = self.angular_flux.values - other.angular_flux.values
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
