r"""The diffusion k-eigenvalue solver on the operator algebra (#290 P5).

Poses the multigroup diffusion criticality problem

.. math::

    \underbrace{(L + C - S - B)}_{A}\,\psi \;=\; \frac{1}{k}\,F\,\psi

on the **scalar composite** ``FullField(bulk=ScalarFlux,
boundary=ScalarBoundaryFlux)`` (campaign ruling 1) and drives it through
the ONE cross-method engine,
:func:`~orpheus.numerics.eigenvalue.power_iteration`, via the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` protocol — the
same boundary the SN / CP / homogeneous methods plug into. The
:class:`DiffusionSolver` here is the modern successor of the MATLAB-port
island in :mod:`orpheus.diffusion.solver` (which P6 retires): where the
island hand-inlined :math:`(L + C - S)\varphi` into a scipy matvec with
a hardcoded 2-group ``[::-1]`` scatter flip, this solver composes the
P4 operator leaves and is ng-generic **by construction**.

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
(:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
+ :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`),
and the production ``F`` (and :meth:`compute_production_rate`, the
ERR-052 renormalization anchor) is :math:`\nu\Sigma_f` ONLY. (SN poses
(n,2n) production-side instead — both are consistent posings of the
same balance; the posing table in :mod:`orpheus.numerics.eigenvalue`
documents the family.)

**Boundary conditions resolve from the MESH's own declarations** — the
single source. Each face's :class:`~orpheus.geometry.mesh.BC` tag (from
``Mesh1D(..., bc_left=..., bc_right=...)``) is resolved to a typed
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw`, realized by the
P3 :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
into its albedo operator :math:`J^- = \mathcal{A} J^+`, and assembled
into ``B`` — mirroring ``SNMesh._resolve_bcs`` exactly (an undeclared
face defaults to ``BC("reflective")``, the infinite-lattice
convention). Supported tags: ``"vacuum"`` (Marshak :math:`J^- = 0`,
:math:`\mathcal{A} = 0` — ruling 3: vacuum MEANS zero incoming),
``"reflective"`` (:math:`\mathcal{A} = 1`), ``"albedo"``
(:math:`\mathcal{A} = \alpha` from ``params["albedo"]``), and
``"zero_flux"`` (the honestly-named Dirichlet idealization,
:math:`\mathcal{A} = -1` — what the legacy island mis-called
"vacuum"). A ``"white"`` tag is deliberately absent: at the P1 level
white coincides with reflective (the P3 realizer docstring's
coincidence note) — declare ``reflective`` or ``albedo``.
(P7 note: this per-method tag→law→realize resolution is the THIRD
spelling of the pattern — ``SNMesh._resolve_bcs`` +
homogenization's method layer + this — and is exactly what the
recorded ``TransportMethod`` trigger unifies.)

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
from typing import TYPE_CHECKING, ClassVar, Optional

import numpy as np

from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    ReflectiveBoundary,
    VacuumInflow,
    ZeroFluxBoundary,
)
from orpheus.geometry.mesh import BC
from orpheus.numerics.eigenvalue import power_iteration
from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.axis import face_labels
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry.boundary import BoundaryTraceLaw
    from orpheus.geometry.mesh import Mesh1D
    from orpheus.numerics.operator import LinearOperator
    from orpheus.transport.mesh.axis import FaceLabel
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
    mesh : MaterialMesh
        The mesh + materials carrier. Must be 1-D (slab / cylinder /
        sphere — refused otherwise by :class:`LeakageOperator`); its
        axes' :class:`~orpheus.geometry.mesh.BC` declarations are the
        boundary-condition source (``None`` defaults to reflective).
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

    #: BC-tag → boundary-law class, the diffusion method's registry
    #: (the SN precedent: ``SNMesh.BOUNDARY_OPERATOR_REGISTRY``).
    BOUNDARY_OPERATOR_REGISTRY: "ClassVar[dict[str, type]]" = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
        "albedo": AlbedoBoundary,
        "zero_flux": ZeroFluxBoundary,
    }

    def __init__(
        self,
        mesh: "MaterialMesh",
        *,
        keff_tol: float = 1e-10,
        flux_tol: float = 1e-9,
    ) -> None:
        self.mesh = mesh
        self.mat_xs: "MaterialXSField" = mesh.material_xs_field()
        self.keff_tol = float(keff_tol)
        self.flux_tol = float(flux_tol)

        # The operator family — P4 leaves + the shared transport algebra.
        # L constructs FIRST: it owns the honest 1-D refusal, before any
        # trace/space construction can trip on a multi-D mesh with a
        # duller error.
        self.leakage = LeakageOperator(mesh)          # refuses ndim != 1
        space = mesh.scalar_full_field_space
        collision = MultiplicationOperator(
            self.mat_xs.total_cross_section_field, space=space,
        )
        # The full K_iso pair (loss-side (n,2n) — module docstring);
        # IsotropicN2N contributes exactly zero on a Σ₂-free mixture.
        scattering = (
            IsotropicScattering(self.mat_xs, space=space)
            + IsotropicN2N(self.mat_xs, space=space)
        )
        self.boundary = DiffusionBoundaryOperator(mesh, self._resolve_bcs())
        self.loss = self.leakage + collision - scattering - self.boundary
        self.fission = FissionOperator.from_solver_data(
            mat_xs=self.mat_xs, full_field_space=space,
        )

        #: The zero composite freezing the flat layout.
        self.template = FullField.zeros(
            bulk=ScalarFlux, boundary=ScalarBoundaryFlux, mesh=mesh,
        )
        #: The campaign-ruled exact resolvent (eager LU at construction).
        self.resolvent = MatrixInverseOperator(
            FlattenedOperator(self.loss, self.template)
        )
        self._n_flat = int(np.asarray(self.template.to_flat()).size)

    # ── Boundary-condition resolution (mesh tags → realized laws) ─────

    def _resolve_bcs(self) -> "dict[str, LinearOperator]":
        r"""Resolve the mesh's per-face BC declarations into realized
        albedo operators — the diffusion mirror of ``SNMesh._resolve_bcs``
        (ONE loop over the same ``face_labels`` inventory the trace is
        built from; ``None`` defaults to reflective)."""
        default = BC("reflective")
        realizer = DiffusionBoundaryRealizer()
        realized: "dict[str, LinearOperator]" = {}
        for label in face_labels(self.mesh.axes):
            bc = self.mesh.axes[label.axis_index].bc[label.endpoint] or default
            law = self._law_from_tag(bc, label)
            method_space = DiffusionMethodSpace.for_face(
                mesh=self.mesh,
                face=label.face_name,
                trace=self.mesh.scalar_trace,
            )
            realized[label.face_name] = realizer.realize(law, method_space)
        return realized

    def _law_from_tag(self, bc: "BC", label: "FaceLabel") -> "BoundaryTraceLaw":
        r"""Construct the typed boundary law a :class:`BC` tag declares."""
        law_cls = self.BOUNDARY_OPERATOR_REGISTRY.get(bc.kind)
        if law_cls is None:
            supported = ", ".join(
                f"'{k}'" for k in sorted(self.BOUNDARY_OPERATOR_REGISTRY)
            )
            raise ValueError(
                f"Diffusion solver does not support boundary condition "
                f"'{bc.kind}' on face '{label.face_name}'. "
                f"Supported: {supported}."
            )
        if law_cls is ReflectiveBoundary:
            # Reflect across the face's own axis (the SN convention);
            # at the P1 level only the albedo=1 scalar survives.
            return ReflectiveBoundary(
                axis=AXIS_NAMES[label.axis_index], albedo=1.0,
            )
        if law_cls is AlbedoBoundary:
            try:
                return AlbedoBoundary(albedo=float(bc.params["albedo"]))
            except KeyError as exc:
                raise ValueError(
                    f"BC('albedo') on face '{label.face_name}' requires an "
                    f"'albedo' parameter; got params={bc.params!r}."
                ) from exc
        return law_cls()

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
        return np.asarray(self.fission.apply(psi).to_flat()) / keff

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""EXACT inner solve :math:`\psi = A^{-1} q` — one LU
        back-substitution; the initial guess is irrelevant."""
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
        ).evaluate(psi.bulk.values)

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        r"""The integrated eigenvalue relation
        :math:`k = P(\psi) / \langle 1, (A\psi)_{\rm bulk}\rangle_V`
        (= production / (absorption + leakage), by the column-sum
        theorem + divergence telescoping — module docstring)."""
        psi = self.unflatten(flux_distribution)
        production = self.compute_production_rate(flux_distribution)
        loss_rate = self._volume_integral(self.loss.apply(psi).bulk.values)
        return production / loss_rate

    def converged(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: np.ndarray,
        flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        r"""``|Δk| < keff_tol`` AND relative flux change ``< flux_tol``
        (the SN convergence contract; no legacy ``print``)."""
        if iteration <= 2:
            return False
        dk = abs(keff - keff_old)
        dphi = np.linalg.norm(flux_distribution - flux_old) / max(
            np.linalg.norm(flux_distribution), 1e-30,
        )
        return bool(dk < self.keff_tol and dphi < self.flux_tol)


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
    mesh: "MaterialMesh"


def solve_diffusion_1d(
    materials: "dict[int, Mixture]",
    mesh: "Mesh1D",
    *,
    keff_tol: float = 1e-10,
    flux_tol: float = 1e-9,
    max_outer: int = 1000,
) -> DiffusionResult:
    r"""Solve the 1-D multigroup diffusion k-eigenvalue problem.

    The modern driver (#290 P5): ``mesh + materials →``
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` ``→``
    resolved boundary laws ``→`` the operator family
    :math:`A = L + C - S - B`, :math:`F` ``→`` the shared
    :func:`~orpheus.numerics.eigenvalue.power_iteration` with the exact
    resolvent inner solve. Zoned cores are expressed by the mesh's own
    ``mat_ids`` (a 3-zone reflector–fuel–reflector slab is a ``Mesh1D``
    with 3-valued ``mat_ids`` — the legacy ``CoreGeometry`` is not
    needed); boundary conditions by the mesh's ``BC`` declarations
    (module docstring for the supported tags and ruling-3 semantics).

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
    material_mesh = MaterialMesh(mesh, materials)
    solver = DiffusionSolver(
        material_mesh, keff_tol=keff_tol, flux_tol=flux_tol,
    )
    keff, keff_history, flux_flat = power_iteration(
        solver, max_iter=max_outer,
    )
    flux = solver.unflatten(flux_flat)
    current = solver.leakage.face_currents(flux)
    return DiffusionResult(
        keff=keff,
        flux=flux,
        current=current,
        keff_history=keff_history,
        mesh=material_mesh,
    )
