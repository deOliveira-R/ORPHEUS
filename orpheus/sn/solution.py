r"""Typed return type for the SN transport solvers.

Issue #197 PR-TYPED-5 — :class:`Solution` (and, until step 3 U6, ``IterationHistory``)
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
role leaves) covers both problem kinds.  Since step 3 of the consumers
campaign (2026-09-17) the carrier is the pair (Problem, posing) plus the
Strategy that produced it and the records — ``mesh``, a kind-typed
``outcome`` (the question, the returned STATE, the answer and the gauge that
picked the representative), ``strategy``, ``certificate`` and ``record`` —
and the KIND is the outcome's TYPE.  ⛔ Until step 3 the kind was read off
an optional ``keff`` through ``is_eigenvalue()`` / ``is_fixed_source()``, and
the flux members were stored fields; they are DERIVED readers of the state
now, and the convergence diagnostics read the :class:`IterationRecord`
directly (``IterationHistory``, the pre-step-3 view over it, survived as a
one-cycle reading through U2 and retired at U6 — its tree readings are the
record's own ``leaf_iterations`` / ``trajectory``, its magnitudes the
certificate's typed evidence).

Reads as the math (``coding-elegance`` Pattern 1 — match the algebra of
the domain)::

    sol = solve_sn(...)                    # Solution[EigenOutcome]
    sol.outcome.keff                       # λ under the k map
    sol.outcome.dominance_ratio()          # |k_n / k_{n-1}|
    sol.reaction_rate_density(sig_a)       # σ · φ
    sol.compare(other, rtol=1e-12)         # SolutionDiff

Two discrimination axes
=======================

The solution family discriminates along TWO independent axes, and the
axes deliberately use DIFFERENT mechanisms:

* **Problem kind** (source-driven vs eigen) is a **type PARAMETER** —
  ``SolutionBase[O]`` with ``O`` the kind-typed OUTCOME
  (:class:`~orpheus.numerics.outcome.EigenOutcome` /
  :class:`~orpheus.numerics.outcome.SourceOutcome`; consumers campaign
  step 3, RULED 2026-09-14).  ⛔ Until step 3 (2026-09-17) the kind was a
  *property*: one carrier covered both kinds via an optional ``keff``, read
  through ``is_eigenvalue()`` / ``is_fixed_source()`` — the kind derived
  from the ANSWER, one tier too late.  That ruling's reasoning was half
  right and is kept here as history: the two kinds DO share every
  operation (homogenizing a fixed-source flux is as meaningful as
  homogenizing an eigenmode — which is why the kind is a parameter and not
  a second pair of classes), but they do NOT share the realization: an eigen
  answer is a RAY plus λ plus the scale section that picked the
  representative, a source answer a COSET plus the kernel section, and their
  adjoints differ in arity (nullary vs a detector).  A parameter carries the
  first fact; the outcome's type carries the second; nothing is ``None``.

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

The state-on-domain law
=======================

:class:`SolutionBase.__post_init__` validates ONE invariant: the outcome's
state is an element of the Problem's coupled space (content identity —
campaign 1 step 6), which is the space the recorded question is posed on.
A cross-Problem pairing is refused; the flux members need no guard of their
own because they are READ off that one state (until step 3 each stored
field carried its own space-content check, and System B's presence was
re-derived by a hand-written biconditional against the mesh — the state's
ARITY answers it now, by construction).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Generic, Self, TypeVar, cast

import numpy as np

from orpheus.numerics.outcome import (
    EigenOutcome,
    Evidence,
    ExitCertificate,
    Measured,
    NotApplicable,
    SourceOutcome,
)

if TYPE_CHECKING:
    from .mesh.augmented_mesh import SNMesh
    from .splitting import Splitting
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.transport.fields.angular_flux import AngularFlux
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
    "Solution",
    "SolutionBase",
    "SolutionDiff",
]


O = TypeVar("O", EigenOutcome, SourceOutcome)
"""The KIND axis — a CONSTRAINED type variable with exactly the two kinds
(consumers campaign step 3, RULED 2026-09-14): a ``Solution[EigenOutcome]``
answers an eigen question (a ray, λ, a scale gauge), a ``Solution[SourceOutcome]``
an affine one (a coset, a kernel gauge).  Constrained rather than bounded so a
``Solution[Any]`` cannot smuggle a third kind in, and so pyright resolves
``sol.outcome.keff`` to a float on the first and to an ERROR on the second."""


@dataclass(frozen=True)
class SolutionBase(Generic[O]):
    r"""Shared carrier for SN transport solutions — the role-agnostic base.

    **A Solution is the pair (Problem, posing) plus the Strategy that produced
    it and the records** (R-cc2; step 3 of the consumers campaign, 2026-09-17):

    * :attr:`mesh` — the Problem (the hub; ``SNProblem`` at #412), the base
      point every other member is relative to;
    * :attr:`outcome` — the kind-typed answer, FUSED with the question it
      answered, the returned STATE and the gauge that picked the
      representative (:class:`~orpheus.numerics.outcome.EigenOutcome` /
      :class:`~orpheus.numerics.outcome.SourceOutcome`).  **The kind IS this
      member's type** — never a ``keff is not None`` read off the answer, which
      is what this class did until step 3 (see the module docstring's history);
    * :attr:`strategy` — the Strategy VALUE the solve drove: the
      :class:`~orpheus.sn.splitting.Splitting` (the labelled piece set, the
      schedule); the budgets and tolerances ride the record, per level;
    * :attr:`certificate` — what the exit MEASURED about the returned state,
      member by member, as typed evidence
      (:class:`~orpheus.numerics.outcome.ExitCertificate`);
    * :attr:`record` — the Strategy's path: the
      :class:`~orpheus.numerics.convergence.IterationRecord` tree.

    Both solution ROLES share this carrier; the role is the concrete type:
    :class:`Solution` (the FORWARD solve) carries the forward-physics
    operations (:meth:`~Solution.homogenize`, :meth:`~Solution.condense`,
    :meth:`~Solution.reaction_rate_density`); :class:`AdjointSolution` (the
    DAGGERED solve) does not — the verb set differs by ROLE, which is why the
    role is a class, and does NOT differ by KIND, which is why the kind is a
    type PARAMETER (``coding-standards`` "Type vs property": an axis that
    changes neither the arithmetic nor the shape may be a parameter — two
    leaves, not four).  ``SolutionBase`` itself is NOT instantiable.

    **The flux members are DERIVED, not stored** (RULED F3/F3b): the outcome's
    ``state`` is the returned iterate WHOLE — the one-system coupled state on a
    seedless mesh, the two-system one on a carrying (ray-bearing) mesh — and
    :attr:`angular_flux`, :attr:`boundary_flux`, :attr:`radial_characteristic`
    and :attr:`scalar_flux` are readers of it, under their historical names.
    So the presence of the ray member is the state's ARITY (no biconditional
    against the mesh to enforce), the scalar flux is :math:`\int\psi\,d\Omega`
    of the cell-average moment (no marginal-axes guard to keep it honest), and
    the one invariant left to check at construction is the STATE-ON-DOMAIN law:
    the state lives on the Problem's coupled space, which is the space the
    recorded question is posed on (a cross-Problem pairing is refused).

    ⚠ On the eigen path the derived scalar flux is :math:`\int\psi\,d\Omega` of
    the RETURNED ψ (one source-iteration step polished against the converged
    fission source, #448) — ``[M]`` 7.4e-11 relative from the power iteration's
    own converged φ that this member stored until step 3 (worst of 16 finalize
    cases; the artefacts were re-baselined with that ratio recorded, U2e).
    """

    mesh: "SNMesh"
    outcome: O
    strategy: "Splitting"
    certificate: ExitCertificate
    record: "IterationRecord"

    def __post_init__(self) -> None:
        # Role closure (unchanged): the role set is closed ({forward, adjoint})
        # and a role-less carrier is not a value that exists (Pattern 4).
        if type(self) is SolutionBase:
            raise TypeError(
                "SolutionBase is the role-agnostic carrier base and is "
                "not instantiable — construct Solution (forward) or "
                "AdjointSolution (adjoint)."
            )
        # The STATE-ON-DOMAIN law: the returned state is an element of the
        # Problem's coupled space (content identity — campaign 1 step 6), and
        # the recorded question is posed on that space.  A same-hub cross-kind
        # pairing is a LEGAL different solve (both kinds share the space); a
        # cross-hub one is refused here.  ``CoupledField.space`` is derived
        # from the members (step 3 U1).
        state_space = self.outcome.state.space
        if state_space != self.mesh.system.space:
            raise ValueError(
                f"{type(self).__name__}: the returned state lives on "
                f"{state_space!r}, not on this Problem's coupled space "
                f"{self.mesh.system.space!r} — a Solution's state is an element "
                "of its own Problem's space (the state-on-domain law)."
            )
        posed_on = _posed_domain(self.outcome)
        if posed_on is not None and posed_on != state_space:
            raise ValueError(
                f"{type(self).__name__}: the recorded question is posed on "
                f"{posed_on!r} but the state lives on {state_space!r} — the "
                "outcome's posing and its state must share one space."
            )

    # ── the state, and the flux members READ off it ──────────────────

    @property
    def state(self) -> "CoupledField":
        r"""The returned iterate WHOLE — the outcome's state (a one- or two-system coupled field)."""
        return self.outcome.state

    @property
    def angular_flux(self) -> "TimedFullField":
        r"""System A's composite: the per-ordinate angular flux ψ (bulk ⊕ trace).

        The arm's own convention, whole: a multi-moment (LD) closure's φ̂
        slopes ride the trailing moment axis on EVERY entry now (until step 3
        the eigen/adjoint tail stripped them to the cell average while the
        fixed-source arms kept them — two conventions for one member).
        """
        return cast("TimedFullField", self.state.systems[0])

    @property
    def radial_characteristic(self) -> "RadialCharacteristicField | None":
        r"""System B's converged ψ½ state on a carrying (R12a) mesh — the state's
        second member; ``None`` exactly when the state has one system.  Presence
        is the state's ARITY: no wiring guard is needed to keep it honest."""
        if self.state.n_systems == 1:
            return None
        return cast("RadialCharacteristicField", self.state.systems[1])

    @property
    def boundary_flux(self) -> "AngularBoundaryFlux":
        r"""Boundary face state — the composite's owned trace (``sol.boundary_flux``)."""
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

        boundary = self.angular_flux.boundary
        if not isinstance(boundary, AngularBoundaryFlux):
            raise TypeError(
                f"{type(self).__name__}.boundary_flux: the solution composite "
                f"carries {type(boundary).__name__}, not AngularBoundaryFlux — "
                "the flux-composite contract is broken."
            )
        return boundary

    @functools.cached_property
    def scalar_flux(self) -> "ScalarFlux":
        r"""The scalar flux :math:`\phi = \int\psi\,d\Omega` of the returned ψ's
        cell-average moment — DERIVED (RULED F3b): one quantity, one
        representation, computed once per Solution (a cached reading of a
        frozen state, so identity reads such as ``adj.importance is
        adj.scalar_flux`` hold).  On the eigen path this is the gauged,
        polished ψ's marginal — the section is applied to the returned state
        at the mint, so the derived φ sits on it too."""
        from orpheus.transport.fields.angular_flux import AngularFlux

        cell_average = AngularFlux(
            values=self.mesh.cell_average_moment(np.asarray(self.angular_flux.interior.values)),
            space=self.mesh.angular_bulk_space,
        )
        return cell_average.integrate_angular()

    # ── convergence (delegates to the record) ────────────────────────

    def converged(self) -> bool:
        """Can this answer be trusted — did EVERY level of the solve converge?
        (``IterationRecord.fully_converged``; the tree-wide question, #340 F1.)"""
        return self.record.fully_converged


    # ── the comparison (regression / refactor checks) ────────────────

    def compare(self, other: Self, *, rtol: float = 1e-12) -> "SolutionDiff":
        r"""Return a field-by-field difference summary against ``other``.

        Same ROLE (``Self``-typed — a forward flux and an importance map are
        different physical quantities), same KIND (an eigen answer and a
        source answer are different QUESTIONS — the kind is the outcome's
        type, so a cross-kind comparison is refused rather than silently
        skipping the eigenvalue channel, which is what the pre-step-3
        ``keff is not None`` branch did), same discrete phase space (the
        layout gate on the collapses — content identity on the
        constituents).  Compares λ (eigen kind) and the :math:`L^\infty`
        norms of the flux deltas.
        """
        if type(other) is not type(self):
            raise TypeError(
                f"{type(self).__name__}.compare: role mismatch — comparing "
                f"against {type(other).__name__}.  A forward flux and an "
                "importance map are different physical quantities; "
                "same-role comparison only (the Self-typed contract, "
                "enforced at runtime for untyped callers)."
            )
        if type(other.outcome) is not type(self.outcome):
            raise TypeError(
                f"{type(self).__name__}.compare: kind mismatch — a "
                f"{type(self.outcome).__name__} against a "
                f"{type(other.outcome).__name__}.  An eigen answer (a ray and "
                "λ) and a source answer (a coset) are different questions; "
                "same-kind comparison only."
            )
        if not self.mesh.same_phase_space(other.mesh):
            raise ValueError(
                f"{type(self).__name__}.compare: the solutions realize "
                "different discrete phase spaces — comparison is defined "
                "only across solves sharing the same constituents (see "
                "MaterialMesh.same_phase_space)."
            )
        keff_abs: Evidence
        if isinstance(self.outcome, EigenOutcome) and isinstance(other.outcome, EigenOutcome):
            keff_abs = Measured(abs(self.outcome.lam - other.outcome.lam))
            reference = abs(other.outcome.lam)
            keff_ok = reference > 0.0 and keff_abs.value <= rtol * reference
        else:
            keff_abs = NotApplicable("the source kind carries no eigenvalue")
            keff_ok = True
        ang_diff = self.angular_flux.interior.values - other.angular_flux.interior.values
        sca_diff = self.scalar_flux.values - other.scalar_flux.values
        ang_linf = float(np.abs(ang_diff).max()) if ang_diff.size else 0.0
        sca_linf = float(np.abs(sca_diff).max()) if sca_diff.size else 0.0
        sca_norm = float(np.abs(other.scalar_flux.values).max())
        flux_ok = (sca_norm == 0.0) or (sca_linf <= rtol * sca_norm)
        return SolutionDiff(
            keff_abs=keff_abs,
            angular_flux_linf=ang_linf,
            scalar_flux_linf=sca_linf,
            within_tolerance=bool(flux_ok and keff_ok),
        )


def _posed_domain(outcome: "EigenOutcome | SourceOutcome"):
    """The space the recorded question is posed on (``None`` when the operator declares no ends)."""
    if isinstance(outcome, EigenOutcome):
        return outcome.posing.pencil.lhs.domain
    return outcome.posing.operator.domain


@dataclass(frozen=True)
class Solution(SolutionBase[O]):
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
    docstring).  One generic type covers both PROBLEM KINDS — the kind is
    the type parameter, ``Solution[EigenOutcome]`` /
    ``Solution[SourceOutcome]`` (step 3 of the consumers campaign; until then
    an optional ``keff`` / ``history`` carried the kind by value — #197
    PR-TYPED-5, which replaced the legacy ``SNFixedSourceResult`` /
    ``SNResult`` pair).

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
    >>> sol.outcome.keff                                        # doctest: +SKIP
    True
    >>> sol.outcome.dominance_ratio()                           # doctest: +SKIP
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
            if not fine.same_phase_space(adjoint.mesh):
                raise ValueError(
                    "Solution.homogenize: adjoint solves a different discrete "
                    "phase space — the importance must come from an adjoint "
                    "solve sharing this solution's constituents (same geometry "
                    "mesh, quadrature, and materials OBJECTS, same scheme; see "
                    "MaterialMesh.same_phase_space)."
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
            :meth:`~orpheus.transport.mesh.material_mesh.MaterialMesh.same_phase_space`).
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
            if not fine.same_phase_space(adjoint.mesh):
                raise ValueError(
                    "Solution.condense: adjoint solves a different discrete "
                    "phase space — the importance must come from an adjoint "
                    "solve sharing this solution's constituents (see "
                    "MaterialMesh.same_phase_space)."
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
class AdjointSolution(SolutionBase[O]):
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
    * ``outcome.keff`` (an :class:`~orpheus.numerics.outcome.EigenOutcome`
      over the hub's ``eigen_posing.H()``) is the eigenvalue of the daggered
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
    keff_abs : Evidence
        :class:`~orpheus.numerics.outcome.Measured` — the absolute eigenvalue
        difference :math:`|k_a - k_b|` — for the eigen kind;
        :class:`~orpheus.numerics.outcome.NotApplicable` for the source kind
        (``compare`` refuses a cross-kind pair, so the two never mix).
    angular_flux_linf : float
        :math:`L^\infty` norm of the angular-flux delta.
    scalar_flux_linf : float
        :math:`L^\infty` norm of the scalar-flux delta.
    within_tolerance : bool
        Aggregate verdict: True iff every available channel met the
        comparison rtol.
    """

    keff_abs: Evidence
    angular_flux_linf: float
    scalar_flux_linf: float
    within_tolerance: bool
