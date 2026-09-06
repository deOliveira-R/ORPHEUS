r"""The diffusion k-eigenvalue solver on the operator algebra (#290 P5).

Poses the multigroup diffusion criticality problem

.. math::

    \underbrace{(L + C - S - B)}_{A}\,\psi \;=\; \frac{1}{k}\,F\,\psi

on the **scalar composite** ``FullField(interior=ScalarFlux,
boundary=ScalarBoundaryFlux)`` (campaign ruling 1) and drives it through
the ONE cross-method engine,
:func:`~orpheus.numerics.eigenvalue.power_iteration`, via the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` protocol — the
same boundary the SN / CP / homogeneous methods plug into. This
:class:`DiffusionSolver` replaced the MATLAB-port island that
previously lived at this module path (retired at #290 P6): where the
island hand-inlined :math:`(L + C - S)\varphi` into a scipy matvec with
a hardcoded 2-group ``[::-1]`` scatter flip and a BiCGSTAB inner
iteration, this solver composes the P4 operator leaves and is
ng-generic **by construction**. Before retirement the two were gated
equivalent on the CORE1D PWR problem (``|Δk| < 1e-8``, identical
discretizations in exact arithmetic — the P5 legacy-bridge gate, which
retired with the island it exercised).

Design decisions (the P5 record)
================================

**The protocol carrier is the flat composite vector.** The
``EigenvalueSolver`` boundary passes an opaque ``np.ndarray`` flux
distribution; here it is the composite's own flat serialization
(``FullField.to_flat`` — bulk C-ravel ⊕ trace buffer). The trace DOFs
``(J⁺, J⁻)`` are honestly PART of the eigenvector: the converged mode
carries its boundary partial currents, typed and inspectable, instead
of burying them in a post-processing gradient reconstruction (the
island's ``_boundary_gradient``). Conversion happens at exactly two
sites (:meth:`DiffusionSolver.unflatten` and the template-frozen
:class:`~orpheus.numerics.flat_operator.FlattenedOperator`).

**The inner solve is EXACT** — the campaign-ruled resolvent

.. code-block:: python

    resolvent = MatrixInverseOperator(FlattenedOperator(A, template))

one eager LU at construction, one back-substitution per outer
iteration. No scattering inner iteration exists at all: the loss ``A``
carries the full multigroup coupling, so ``solve_fixed_source`` is a
single application of :math:`A^{-1}`. (NEVER route through the
structure-keyed ``A.inverse()`` — the Green splitting diverges for
fine-mesh elliptic operators; campaign ruling.)

**The k-update is the integrated eigenvalue relation itself.** At any
iterate, :math:`k = P(\psi) / \langle 1, (A\psi)_{\rm bulk}\rangle_V`:
the production rate over the volume-integrated bulk loss rate. By the
P4-gated column-sum theorem :math:`\mathbf{1}^T(C - S) = \Sigma_a` and
the telescoping of ``L``'s conservative divergence (interior flows
cancel; boundary flows survive as the outward leakage), the
denominator decomposes as

.. math::

    \langle 1, (A\psi)_{\rm bulk}\rangle_V
        \;=\; \text{absorption} + \text{leakage}
              \;(-\; \text{net (n,2n) gain}),

i.e. the island's ``p_rate / (a_rate + leakage)`` — but derived through
the SAME operator that defines the fixed point, so no term can be
forgotten (the leakage term is structural, not hand-added). ``B``
contributes nothing to the bulk block, and the trace rows of
:math:`A\psi` (the constraint rows) are excluded by construction.

**The (n,2n) convention is loss-side** — mirroring the homogeneous
solver (:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`):
``S`` is the full K_iso pair
:math:`\Sigma_{s0}^T + 2\Sigma_2^T`
(:class:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering`
+ :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`),
and the production ``F`` (and :meth:`compute_production_rate`, the
ERR-052 renormalization anchor) is :math:`\nu\Sigma_f` ONLY. (SN poses
(n,2n) production-side instead — both are consistent posings of the
same balance; the posing table in :mod:`orpheus.numerics.eigenvalue`
documents the family.)

**Boundary conditions are the PHASE SPACE's contract, not the
solver's** (#290 P7a): the solver consumes a
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`, which
resolved and REALIZED its per-face laws at construction — each face's
:class:`~orpheus.geometry.mesh.BC` tag became a typed
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` and then the
albedo operator :math:`J^- = \mathcal{A} J^+` in ``mesh.bc``
(``SNMesh.bc`` parity; supported tags, ruling-3 semantics, and the
deliberate ``"white"`` absence are documented on
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`). ``B`` reads
``mesh.bc``; a solver on a phase space with unresolved BCs is
unrepresentable. The resolution itself is the ONE shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body
under the :class:`~orpheus.transport.method.TransportMethod` Protocol
(#290 P7b — it replaced the twin per-method ``_resolve_bcs`` loops;
each mesh keeps only its ``realize_boundary_law`` arm).

**Normalization**: :func:`power_iteration` renormalizes each iterate to
unit production rate (:meth:`compute_production_rate` — the
:class:`~orpheus.numerics.eigenvalue.ProductionRateSolver` contract,
ERR-052), so the returned mode carries the canonical convention
:math:`\int_V \nu\Sigma_f\,\phi\,dV = 1`. Rescaling to an absolute flux
at a target reactor power is a single downstream multiplication — the
island's hardcoded ``e_per_fission`` power window and its
``fi /= max(|fi|)`` conditioning hack are both retired by this
contract (#270 diffusion arm).

Verification (the P5 gates, ``tests/diffusion/test_solver.py``)
===============================================================

* cross-engine: ``power_iteration`` ≡
  :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue` on the
  materialized ``(A, F)`` pair at :math:`10^{-10}` in k (and the full
  composite eigenvector), across the albedo family;
* L2 infinite-medium: reflective diffusion k ≡ homogeneous
  :math:`k_\infty` (the shared-K_iso, independently-posed loss) —
  including the 3-group asymmetric case the island's flip trick
  structurally cannot represent;
* the legacy bridge: modern k ≡ island k on the CORE1D PWR problem
  under the ruling-4 bit-identical encoding (both discretizations are
  the same math in exact arithmetic);
* solution trace semantics per law (vacuum :math:`J^- = 0`, albedo
  :math:`J^- = \alpha J^+`, reflective :math:`J_{\rm net} = 0`,
  zero-flux :math:`\phi_\Gamma = 0`), the integrated balance identity,
  and k-monotonicity in the albedo.

.. seealso::

   :mod:`orpheus.diffusion.operators`
       The P4 operator family (``L``, ``B``) + the block structure and
       discretization derivation.
   ``.claude/plans/diffusion_crosswalk.md``
       The ``(J⁺, J⁻)`` convention contract.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

import numpy as np

from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
from orpheus.numerics.convergence import (
    IterationRecord,
    StoppingCriterion,
    warn_if_unconverged,
)
from orpheus.numerics.eigenvalue import power_iteration
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.lift import BulkLift
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicFission,
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry.mesh import Mesh1D
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["DiffusionResult", "DiffusionSolver", "solve_diffusion_1d"]


class DiffusionSolver:
    r"""Multigroup diffusion k-eigenvalue solver on the scalar composite.

    Implements the :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver`
    protocol (including the
    :class:`~orpheus.numerics.eigenvalue.ProductionRateSolver`
    production-rate extension) over the flat composite vector, with the
    EXACT resolvent inner solve — see the module docstring for the full
    design record.

    Parameters
    ----------
    mesh : DiffusionMesh
        The diffusion phase space (#290 P7a) — 1-D, bounded, with its
        per-face boundary laws already realized at construction
        (``mesh.bc``). Promote a plain
        :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
        via :meth:`DiffusionMesh.from_material_mesh
        <orpheus.diffusion.augmented_mesh.DiffusionMesh.from_material_mesh>`.
    keff_tol, flux_tol : float
        Outer-iteration convergence tolerances on ``|Δk|`` and the
        relative flux-distribution change.

    Attributes
    ----------
    loss : LinearOperator
        The composed :math:`A = L + C - S - B` (typed, un-materialized).
    fission : FissionOperator
        The shared rank-1 production :math:`F = \chi \otimes \nu\Sigma_f`.
    leakage : LeakageOperator
        The ``L`` leaf — also serves the net-current-profile
        reconstruction (:meth:`LeakageOperator.face_currents`).
    resolvent : MatrixInverseOperator
        :math:`A^{-1}` over the flat composite (eager LU).
    template : FullField
        The zero composite fixing the flat layout (shapes, spaces, mesh).
    """

    def __init__(
        self,
        mesh: "DiffusionMesh",
        *,
        keff_tol: float = 1e-10,
        flux_tol: float = 1e-9,
    ) -> None:
        self.mesh = mesh
        self.mat_xs: "MaterialXSField" = mesh.material_xs_field()
        self.keff_tol = float(keff_tol)
        self.flux_tol = float(flux_tol)

        # The operator family — P4 leaves + the shared transport
        # algebra, assembled on the phase space's own composite carrier.
        # Admission (1-D, bounded, realized BCs) is the MESH's
        # construction contract: a DiffusionMesh cannot exist otherwise.
        space = mesh.full_field_space
        self.leakage = LeakageOperator(mesh)
        collision = MultiplicationOperator(
            self.mat_xs.total_cross_section_field, domain=space, codomain=space,
        )
        # The full K_iso pair (loss-side (n,2n) — module docstring);
        # IsotropicN2N contributes exactly zero on a Σ₂-free mixture. The
        # energy bindings are PLAIN-bound on the mesh's scalar bulk (CS4c
        # step 5, R-4: the array carriers are theirs) and LIFTED onto the
        # scalar composite here, once, so the loss composes under the
        # OperatorSum ends guard and assembles for the exact resolvent
        # (the lift embeds the bulk block-diagonal in the composite flat
        # layout — index-identity, bit-for-bit).
        bulk = mesh.bulk_space
        scattering = BulkLift(
            IsotropicScattering.from_material_xs(self.mat_xs, space=bulk)
            + IsotropicN2N.from_material_xs(self.mat_xs, space=bulk),
            domain=space, codomain=space,
        )
        self.boundary = DiffusionBoundaryOperator(mesh)
        self.loss = self.leakage + collision - scattering - self.boundary
        # The fission ENERGY binding (CS4c step 4), lifted the same way:
        # diffusion consumes the scalar dyad on its scalar-bulk composite
        # — the angular composite binding (FissionOperator) is the SN
        # eigen-posing's.
        self.fission = BulkLift(
            IsotropicFission.from_material_xs(self.mat_xs, space=bulk),
            domain=space, codomain=space,
        )

        #: The zero composite freezing the flat layout.
        self.template = FullField.zeros(
            interior=ScalarFlux, boundary=ScalarBoundaryFlux, space=mesh.full_field_space,
        )
        #: The campaign-ruled exact resolvent (eager LU at construction).
        self.resolvent = MatrixInverseOperator(
            FlattenedOperator(self.loss, self.template)
        )
        self._n_flat = int(np.asarray(self.template.to_flat()).size)
        # #340 N4: one DIRECT record per exact inner solve, in order.
        # APPEND-ONLY across solves (``power_iteration`` slices from a
        # high-water mark it takes before its own loop).
        self._inner_records: list[IterationRecord] = []

    # ── Flat ↔ typed conversion ────────────────────────────────────────

    def unflatten(self, flux_distribution: np.ndarray) -> "FullField":
        r"""Rebuild the typed composite from a flat protocol vector."""
        return FullField.from_flat(
            np.asarray(flux_distribution, dtype=float).ravel(), self.template,
        )

    def _volume_integral(self, rate_density: np.ndarray) -> float:
        r"""``Σ_cells V·rate`` over all groups — the mesh's canonical
        ``volume_measure`` (the same reduction
        :class:`IntegratedReactionRate` rides)."""
        ng = rate_density.shape[0]
        flat = np.moveaxis(rate_density, 0, -1).reshape(-1, ng)
        return float(self.mesh.volume_measure(flat).sum())

    # ── EigenvalueSolver protocol ──────────────────────────────────────

    def initial_flux_distribution(self) -> np.ndarray:
        r"""Flat all-ones composite. (Only the bulk seeds the first
        fission source; the trace slots are overwritten by the first
        exact solve.)"""
        return np.ones(self._n_flat)

    def compute_fission_source(
        self, flux_distribution: np.ndarray, keff: float,
    ) -> np.ndarray:
        r"""Fission RHS :math:`F\psi / k` — the shared rank-1 dyad's
        scalar-composite arm (zero trace rows by construction)."""
        psi = self.unflatten(flux_distribution)
        fission = cast("FullField", self.fission.apply(psi))
        return np.asarray(fission.to_flat()) / keff

    @property
    def inner_records(self) -> "Sequence[IterationRecord]":
        r"""Every inner solve this instance has run, in order (#340 N4).

        Satisfies the optional
        :class:`~orpheus.numerics.eigenvalue.RecordingSolver` member. Each
        entry reads ``DIRECT`` — *did not iterate, DID converge* — because
        the inner is one LU back-substitution.

        Recording an exact level looks redundant and is not: without it the
        tree is one deep, and *"diffusion has no inner level"* is then
        indistinguishable from *"diffusion's inner was never recorded"*.
        The ``DIRECT`` child says which. It also means diffusion's
        ``fully_converged`` can only ever fail at the outer, which is the
        true statement about this method.
        """
        return self._inner_records

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""EXACT inner solve :math:`\psi = A^{-1} q` — one LU
        back-substitution; the initial guess is irrelevant."""
        # The default ``IterationBudget()`` has ``limit=0``: an exact solve
        # has no iteration budget, so ``exhausted_budget`` is False by
        # construction and the knob name is inert — nothing will ever advise
        # setting it.
        self._inner_records.append(
            IterationRecord(label="inner(exact resolvent, LU)")
        )
        return np.asarray(
            self.resolvent.apply(np.asarray(fission_source, dtype=float))
        )

    def compute_production_rate(self, flux_distribution: np.ndarray) -> float:
        r"""Total :math:`\int_V \langle\nu\Sigma_f, \phi\rangle\,dV` —
        the typed reaction-rate functional (#270 diffusion arm; the
        ERR-052 renormalization anchor). :math:`\nu\Sigma_f` ONLY:
        (n,2n) is loss-side here (module docstring)."""
        psi = self.unflatten(flux_distribution)
        return IntegratedReactionRate(
            self.mat_xs.fission_production_field
        ).evaluate(psi.interior.values)

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        r"""The integrated eigenvalue relation
        :math:`k = P(\psi) / \langle 1, (A\psi)_{\rm bulk}\rangle_V`
        (= production / (absorption + leakage), by the column-sum
        theorem + divergence telescoping — module docstring)."""
        psi = self.unflatten(flux_distribution)
        production = self.compute_production_rate(flux_distribution)
        loss_rate = self._volume_integral(self.loss.apply(psi).interior.values)
        return production / loss_rate

    def measure_stopping_criteria(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: np.ndarray,
        flux_old: np.ndarray,
    ) -> tuple[StoppingCriterion, ...]:
        r"""``|Δk|`` against ``keff_tol`` and the relative flux change against
        ``flux_tol`` (the SN convergence contract; no legacy ``print``).

        ⛔ Was ``converged(...) -> bool`` until 2026-08-09 (#340 N2b); the
        ``if iteration <= 2`` head is now
        :data:`~orpheus.numerics.eigenvalue.MINIMUM_OUTER_ITERATIONS`.
        """
        dk = abs(keff - keff_old)
        dphi = np.linalg.norm(flux_distribution - flux_old) / max(
            np.linalg.norm(flux_distribution), 1e-30,
        )
        return (
            StoppingCriterion.reading("dk", float(dk), self.keff_tol),
            StoppingCriterion.reading("dphi", float(dphi), self.flux_tol),
        )


@dataclass(frozen=True, eq=False)
class DiffusionResult:
    r"""Converged diffusion eigenmode (the modern, typed result).

    ``flux`` is the full composite — the bulk :class:`ScalarFlux` AND
    the boundary ``(J⁺, J⁻)`` trace — normalized to unit production
    rate (:math:`\int_V \nu\Sigma_f\,\phi\,dV = 1`; rescaling to a
    target power is one multiplication). ``current`` is the axis-signed
    net-current profile at every face, reconstructed by the leakage
    operator (interior faces from the condensed two-point currents,
    boundary faces from the trace).
    """

    keff: float
    flux: "FullField"
    current: np.ndarray            # (ng, nx+1) axis-signed net currents
    keff_history: "list[float]"
    mesh: "DiffusionMesh"

    record: IterationRecord
    """The iteration tree this solve actually ran (#340 N4).

    Ask :attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`
    to learn whether the returned eigenvalue can be trusted.  Diffusion's
    inner is one LU back-substitution, so each child reads ``DIRECT`` — it
    did not iterate and DID converge — and the only level that can fail is
    the outer power iteration.
    """


def solve_diffusion_1d(
    materials: "dict[int, Mixture]",
    mesh: "Mesh1D",
    *,
    keff_tol: float = 1e-10,
    flux_tol: float = 1e-9,
    max_outer: int = 1000,
) -> DiffusionResult:
    r"""Solve the 1-D multigroup diffusion k-eigenvalue problem.

    The modern driver (#290 P5, mesh-layered at P7a): ``mesh +
    materials →``
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (the
    data carrier) ``→`` :meth:`DiffusionMesh.from_material_mesh
    <orpheus.diffusion.augmented_mesh.DiffusionMesh.from_material_mesh>`
    (the promotion — trace built, boundary laws realized) ``→`` the
    operator family :math:`A = L + C - S - B`, :math:`F` ``→`` the
    shared :func:`~orpheus.numerics.eigenvalue.power_iteration` with
    the exact resolvent inner solve. Zoned cores are expressed by the
    mesh's own ``mat_ids`` (a 3-zone reflector–fuel–reflector slab is a
    ``Mesh1D`` with 3-valued ``mat_ids`` — the island's retired
    ``CoreGeometry`` container had no independent content); boundary
    conditions by the mesh's ``BC`` declarations (supported tags and
    ruling-3 semantics on
    :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`).

    Parameters
    ----------
    materials : dict[int, Mixture]
        Material id → :class:`~orpheus.data.macro_xs.mixture.Mixture`,
        keyed by the mesh's ``mat_ids`` values. The diffusion
        coefficient reads through the #290 P1 seam
        (``Mixture.diffusion_coefficient`` — transport-corrected when a
        P1 moment is present, :math:`1/(3\Sigma_t)` exactly otherwise).
    mesh : Mesh1D
        Geometry + zoning + BC declarations (1-D; slab / cylinder /
        sphere through its own areas + volumes).
    keff_tol, flux_tol : float
        Outer convergence tolerances (``|Δk|``, relative flux change).
    max_outer : int
        Outer power-iteration cap.

    Returns
    -------
    DiffusionResult
    """
    diffusion_mesh = DiffusionMesh.from_material_mesh(
        MaterialMesh(mesh, materials)
    )
    solver = DiffusionSolver(
        diffusion_mesh, keff_tol=keff_tol, flux_tol=flux_tol,
    )
    outcome = power_iteration(solver, max_iter=max_outer, budget_name="max_outer")
    keff, keff_history, flux_flat = (
        outcome.keff, outcome.keff_history, outcome.flux_distribution,
    )
    flux = solver.unflatten(flux_flat)
    current = solver.leakage.face_currents(flux)

    # ⭐ #340 N4.7 — see the note at the identical call in
    # :func:`~orpheus.cp.solver.solve_cp` for why this sits directly in the
    # entry and passes no ``balance_defect``.
    #
    # Diffusion's tree can only fail at the OUTER: its inner is an exact LU
    # resolvent, recorded with ``budget = 0`` (a DIRECT level, never
    # truncated).  So this warning is strictly about the power iteration —
    # which is also why ``max_outer`` is always the knob it names here.
    # ``balance_defect=None`` explicitly — see the note at ``solve_cp``.
    warn_if_unconverged(
        outcome.record, where="solve_diffusion_1d", balance_defect=None,
    )

    return DiffusionResult(
        keff=keff,
        flux=flux,
        current=current,
        keff_history=keff_history,
        mesh=diffusion_mesh,
        record=outcome.record,
    )
