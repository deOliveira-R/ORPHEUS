r"""The diffusion operator family on the scalar composite — ``L`` and ``B``.

#290 P4: the two genuinely-new leaves of the diffusion loss

.. math::

    A \;=\; L + C - S - B ,
    \qquad
    (A)\,x \;=\; \tfrac{1}{k}\,F\,x ,

acting on the scalar composite ``FullField(bulk=ScalarFlux,
boundary=ScalarBoundaryFlux)`` (user ruling 1). Everything else is the
SHARED algebra: ``C`` is
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
over :math:`\sigma_t`, ``S`` the K_iso kernels
(:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
[+ :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`]),
``F`` the shared rank-1
:class:`~orpheus.transport.operators.fission.FissionOperator` — all of
which gained their scalar-composite arms at P4. The in-group
:math:`\sigma_{gg}` cancellation between ``C`` and ``S`` REPRODUCES the
removal :math:`\Sigma_r = \sigma_t - \sigma_{s,gg}` as a THEOREM, not an
input (attacker Q2 convention ruling).

The block structure (the #208 partition, transferred verbatim)
==============================================================

On ``V = V_bulk ⊕ V_trace`` with trace unknowns the per-face partial
currents ``(J⁺, J⁻)`` (ruling 2; crosswalk)::

    C, S, F  →  [ A_bb  0 ]   (BULK)
                [ 0     0 ]

    L        →  [ A_bb  A_bs ] (FULL — the elliptic sibling of streaming)
                [ A_sb  A_ss^L]

    B        →  [ 0     0    ] (BOUNDARY — the realized albedo block)
                [ 0     A_ss^B]

``L``'s trace block mirrors the SN streaming convention exactly: the
**outflow-definition defect** on the
:attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.OUTFLOW_ROW`
plus the **inflow identity** on the ``INFLOW_ROW``, so that the loss's
trace rows read

.. math::

    \text{outflow row:}\quad
        J^+ - c_\phi\,\phi_e - c_{J^-}\,J^- \;=\; 0,
    \qquad
    \text{inflow row:}\quad
        J^- - \mathcal{A}\,J^+ \;=\; 0,

the second falling out of ``(L - B)`` because ``B`` emits
:math:`\mathcal{A} J^+` on the inflow row (the P3 realized law) and
``L`` supplies the :math:`J^-` identity.

The discretization (FD ≡ RT0 with mass lumping; attacker Q3)
============================================================

**Interior faces** carry the current-continuous two-point flux with the
series half-cell resistance (the harmonic-mean-of-``D`` form — equal to
the legacy island's arithmetic-mean-of-:math:`\Sigma_{\rm tr}` face
coefficient exactly, since :math:`D \propto 1/\Sigma_{\rm tr}` makes
``HM(D) = 1/(3·AM(Σ_tr))``):

.. math::

    J_f \;=\; -\,g_f\,(\phi_R - \phi_L),
    \qquad
    g_f \;=\; \Bigl(\frac{h_L}{2 D_L} + \frac{h_R}{2 D_R}\Bigr)^{-1} ,

and the bulk row is the conservative divergence
:math:`[(A J)_{i+1/2} - (A J)_{i-1/2}] / V_i` — written with the mesh's
face areas and cell volumes, so the SAME body is correct on slab
(``A ≡ 1``, ``V = h``: bit-identical to ``diff(J)/h``), cylinder, and
sphere (whose pole is not a face: the ``r = 0`` slot simply carries no
trace flow). Interior face currents stay CONDENSED (never trace DOFs);
only the boundary ``(J⁺, J⁻)`` pair survives as trace unknowns.

**Boundary faces** close through the P1 face algebra (Fick from the
edge-cell centre to the face + the P1 dictionary
:math:`\phi_\Gamma = 2(J^+ + J^-)`, :math:`J_{\rm net} = J^+ - J^-`):
with the half-cell resistance :math:`\rho = h_e / (2 D_e)`,

.. math::

    J^+ (\rho + 2) \;=\; \phi_e + (\rho - 2)\,J^- ,
    \qquad\Longrightarrow\qquad
    c_\phi = \frac{1}{\rho + 2}, \quad
    c_{J^-} = \frac{\rho - 2}{\rho + 2} ,

and the bulk edge row reads the net outward trace current
:math:`A_f\,(J^+ - J^-)/V_e` (face-local outward normal at EVERY face —
the crosswalk orientation row; the axis-signed conversion happens here,
the one consumer that needs the vector form).

Eliminating ``(J⁺, J⁻)`` from these rows (the Schur complement onto the
bulk) recovers the classic condensed closures exactly: zero-flux
(``𝒜 = −1``) → :math:`J_{\rm net} = \phi_e/\rho` — the legacy island's
``φ_0/(0.5·dz)`` vacuum-named arm; reflective (``𝒜 = 1``) →
:math:`J_{\rm net} = 0`; Marshak vacuum (``𝒜 = 0``) → the
:math:`d = 2D` extrapolation mixed condition. The composite is the
UN-condensed spelling of the same algebra, with the boundary law
factored into its own operator ``B`` (ruling 1).

The resolvent (P4 design ruling)
================================

:math:`A^{-1}` is the explicit dense direct inverse

.. code-block:: python

    template = FullField.zeros(bulk=ScalarFlux, boundary=ScalarBoundaryFlux, mesh=mesh)
    A_inv = MatrixInverseOperator(FlattenedOperator(A, template))

— NEVER the structure-keyed ``A.inverse()`` (the Green splitting
diverges for fine-mesh elliptic operators; campaign ruling). The
composite space's flat ``shape`` identity feeds ``as_matrix`` its basis
dimension automatically.

References
----------

* ``.claude/plans/diffusion_crosswalk.md`` — the ``(J⁺, J⁻)``
  convention contract (rows, signs, BC family, metric).
* Bell & Glasstone (1970) §3.4; Duderstadt & Hamilton (1976) §5 —
  the P1/Marshak boundary algebra.
* attacker memo ``diffusion_integration_frames.md`` Q2 (operator
  decomposition = shared loss + shared gains) / Q3 (FD-as-is IS RT0
  with mass lumping; Baliga–Patankar equivalence).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping, Optional

import numpy as np

from orpheus.numerics.operator import BlockRole, LinearOperator
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.axis import face_labels
from orpheus.transport.source_sinks import (
    ScalarBoundarySourceSink,
    ScalarSourceSink,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.mesh.material_mesh import MaterialMesh


__all__ = ["DiffusionBoundaryOperator", "LeakageOperator"]


# ─────────────────────────────────────────────────────────────────────
# The two discretization kernels — module-level pure functions so the
# formulas are stated once, readable, and mutation-testable (the P4
# stencil gate monkeypatches THESE to prove its teeth).
# ─────────────────────────────────────────────────────────────────────


def _interior_conductance(D: np.ndarray, h: np.ndarray) -> np.ndarray:
    r"""Current-continuous interior-face conductance ``g_f`` per group.

    .. math::

        g_f \;=\; \Bigl(\frac{h_L}{2 D_L} + \frac{h_R}{2 D_R}\Bigr)^{-1}

    — the series sum of the two half-cell resistances (the RT0 /
    harmonic-mean form; equal to the legacy arithmetic-mean-of-
    :math:`\Sigma_{\rm tr}` coefficient exactly on a uniform mesh).

    Parameters
    ----------
    D : (ng, nx) per-cell diffusion coefficients.
    h : (nx,) cell widths.

    Returns
    -------
    (ng, nx-1) conductances at the interior faces.
    """
    r_left = h[:-1] / (2.0 * D[:, :-1])
    r_right = h[1:] / (2.0 * D[:, 1:])
    return 1.0 / (r_left + r_right)


def _boundary_closure(
    D_edge: np.ndarray, h_edge: float,
) -> tuple[np.ndarray, np.ndarray]:
    r"""The P1 face-closure coefficients ``(c_φ, c_{J⁻})`` at one boundary face.

    From Fick between the edge-cell centre and the face plus the P1
    dictionary (module docstring): with :math:`\rho = h_e/(2 D_e)`,

    .. math::

        c_\phi = \frac{1}{\rho + 2}, \qquad
        c_{J^-} = \frac{\rho - 2}{\rho + 2} .

    Parameters
    ----------
    D_edge : (ng,) edge-cell diffusion coefficients.
    h_edge : the edge-cell width.

    Returns
    -------
    ``(c_phi, c_inflow)`` — each ``(ng,)``.
    """
    rho = h_edge / (2.0 * D_edge)
    return 1.0 / (rho + 2.0), (rho - 2.0) / (rho + 2.0)


@dataclass(frozen=True)
class _FaceClosure:
    r"""Precomputed geometry + closure data for ONE boundary face
    (named intermediates — Pattern 3; every entry has physical identity).
    """

    #: Edge-cell index along the axis (``0`` for the min face, ``-1`` for max).
    edge: int
    #: Index of this face's slot in the ``(nx+1,)`` face-flow buffer.
    flow_slot: int
    #: Axis sign of the OUTWARD normal (``-1`` at min, ``+1`` at max) —
    #: the one site converting face-local net current to the axis vector.
    axis_sign: float
    #: P1 closure coefficient on the edge-cell flux, ``(ng,)``.
    c_phi: np.ndarray
    #: P1 closure coefficient on the inflow partial current, ``(ng,)``.
    c_inflow: np.ndarray


class LeakageOperator(LinearOperator["FullField", "FullField"]):
    r"""The elliptic leakage leaf :math:`L = -\nabla\!\cdot D\nabla` — the
    FULL primitive of the diffusion algebra.

    The diffusion sibling of SN's
    :class:`~orpheus.sn.operators.streaming.StreamingOperator`: the ONE
    operator that couples bulk ↔ boundary (attacker Q2 — the #208
    partition transfers verbatim; the diffusion exception is the
    REALIZATION, elliptic-self-adjoint vs characteristic-triangular,
    never the algebra). Reads ``D`` per cell through the #290 P1 data
    seam (:attr:`MaterialXSField.diffusion_coefficient
    <orpheus.transport.mesh.material_xs_field.MaterialXSField.diffusion_coefficient>`
    = the per-cell gather of ``Mixture.diffusion_coefficient``).

    Action on ``FullField(bulk=ScalarFlux, boundary=ScalarBoundaryFlux)``
    (see the module docstring for the derivation):

    * **bulk** — the conservative FD divergence with condensed interior
      currents; the edge rows read the trace's net outward current
      :math:`A_f (J^+ - J^-)/V_e` (the ``A_bs`` coupling);
    * **outflow trace rows** — the P1 outflow-definition defect
      :math:`J^+ - c_\phi \phi_e - c_{J^-} J^-` (``A_sb`` + the
      outflow part of ``A_ss``);
    * **inflow trace rows** — the identity :math:`J^-`, so that
      ``(L − B)`` reads :math:`J^- - \mathcal{A} J^+` (the SN
      streaming-operator convention, mirrored).

    1-D meshes only today (slab / cylinder / sphere through the mesh's
    own areas + volumes); a multi-D mesh is refused at construction —
    the N-D stencil is a genuine extension, not a broadcast.

    Parameters
    ----------
    mesh : MaterialMesh
        The mesh + materials carrier; supplies ``D`` (through its XS
        field), the widths / areas / volumes, and the scalar trace the
        composite is bound to.
    """

    block_role = BlockRole.FULL

    def __init__(self, mesh: "MaterialMesh") -> None:
        if mesh.ndim != 1:
            raise ValueError(
                f"LeakageOperator supports 1-D meshes today (slab / "
                f"cylinder / sphere); got a {mesh.ndim}-D mesh. The N-D "
                f"diffusion stencil is a deliberate extension seam, not "
                f"a broadcast."
            )
        self.mesh = mesh
        D = np.asarray(mesh.material_xs_field().diffusion_coefficient, float)
        h = np.asarray(mesh.axis_widths[0], float)
        areas = np.asarray(mesh.areas, float)
        #: (ng, nx-1) interior-face conductances (the condensed currents).
        self._conductance = _interior_conductance(D, h)
        #: (nx+1,) face areas — interior slots weight the condensed
        #: currents in the conservative divergence.
        self._areas = areas
        #: (nx,) cell volumes (the divergence denominator).
        self._volumes = np.asarray(mesh.volumes, float)
        # Per boundary face: the closure bundle. Derived from the SAME
        # face_labels inventory the mesh's scalar trace was built from
        # (single source — a curvilinear pole is not a face and never
        # appears here, leaving its flow slot at zero: A(r=0)·J = 0).
        nx = h.size
        closures: dict[str, _FaceClosure] = {}
        for label in face_labels(mesh.axes):
            is_min = label.endpoint == "min"
            edge = 0 if is_min else nx - 1
            c_phi, c_inflow = _boundary_closure(D[:, edge], float(h[edge]))
            closures[label.face_name] = _FaceClosure(
                edge=edge,
                flow_slot=0 if is_min else nx,
                axis_sign=-1.0 if is_min else +1.0,
                c_phi=c_phi,
                c_inflow=c_inflow,
            )
        self._face_closures = closures

    @property
    def domain(self) -> "Optional[FunctionSpace]":
        return self.mesh.scalar_full_field_space

    @property
    def codomain(self) -> "Optional[FunctionSpace]":
        return self.mesh.scalar_full_field_space

    def _parse(self, psi: "FullField") -> "tuple[np.ndarray, ScalarBoundaryFlux]":
        r"""Guard + unpack a scalar composite → ``(φ values, trace)``.

        The shared entry of :meth:`apply` and :meth:`face_currents`
        (single source of the composite parses — Pattern 2)."""
        bulk = psi.bulk
        if not isinstance(bulk, ScalarFlux):
            raise TypeError(
                f"LeakageOperator: the input composite's bulk must be a "
                f"ScalarFlux; got {type(bulk).__name__}."
            )
        trace = psi.boundary
        if not isinstance(trace, ScalarBoundaryFlux):
            raise TypeError(
                f"LeakageOperator: the input composite's boundary must be "
                f"a ScalarBoundaryFlux (the (J⁺, J⁻) trace); got "
                f"{type(trace).__name__}."
            )
        if bulk.mesh is not self.mesh:
            raise ValueError(
                "LeakageOperator: input field and operator must "
                "share the same MaterialMesh instance (mesh-identity "
                f"invariant); got field mesh {bulk.mesh!r} vs operator "
                f"mesh {self.mesh!r}."
            )
        return bulk.values, trace

    def face_currents(self, psi: "FullField") -> np.ndarray:
        r"""Reconstruct the axis-signed net-current profile, ``(ng, nx+1)``.

        The current-density counterpart of :meth:`apply`'s flow assembly
        — and its single source (``apply`` consumes THIS reconstruction
        times the face areas): interior faces carry the condensed
        two-point current :math:`J_f = -g_f\,(\phi_R - \phi_L)`,
        boundary faces the trace's net outward current converted to the
        axis sign (the one vector-form consumer site — the crosswalk
        orientation row). Interior currents are never trace DOFs; this
        reconstruction is how the production current-profile output is
        served (#290 P5 — the modern successor of the legacy island's
        ``DiffusionResult.current``).
        """
        phi, trace = self._parse(psi)
        return self._face_currents_from(phi, trace)

    def _face_currents_from(
        self, phi: np.ndarray, trace: "ScalarBoundaryFlux",
    ) -> np.ndarray:
        ng, nx = phi.shape
        current = np.zeros((ng, nx + 1))
        current[:, 1:-1] = -self._conductance * (phi[:, 1:] - phi[:, :-1])
        for face, c in self._face_closures.items():
            current[:, c.flow_slot] = c.axis_sign * trace.net_current(face)
        return current

    def apply(self, psi: "FullField", /) -> "FullField":
        r"""Return :math:`L\,\psi` — leakage rate density ⊕ trace defect rows."""
        phi, trace = self._parse(psi)

        # ── Bulk block: the conservative divergence ────────────────────
        # Area-weighted, axis-signed flow at every face — the face-current
        # reconstruction (single source: :meth:`face_currents`) times the
        # face areas; interior slots carry the condensed currents,
        # boundary slots the trace's net outward current.
        flow = self._areas * self._face_currents_from(phi, trace)
        out_bulk = np.diff(flow, axis=1) / self._volumes  # (ng, nx)

        # ── Trace block: outflow-definition defect + inflow identity ──
        out_boundary = ScalarBoundarySourceSink.zeros_on(self.mesh)
        for face, c in self._face_closures.items():
            j_plus = trace.outflow_view(face)             # (ng,)
            j_minus = trace.inflow_view(face)             # (ng,)
            slot = out_boundary.face_view(face)
            slot[ScalarTraceSpace.OUTFLOW_ROW] = (
                j_plus - c.c_phi * phi[:, c.edge] - c.c_inflow * j_minus
            )
            slot[ScalarTraceSpace.INFLOW_ROW] = j_minus

        return FullField(
            bulk=ScalarSourceSink.from_mesh(out_bulk, self.mesh),
            boundary=out_boundary,
        )


class DiffusionBoundaryOperator(LinearOperator["FullField", "FullField"]):
    r"""The whole-trace boundary law ``B`` — the ``A_ss`` block of the
    diffusion algebra.

    The scalar mirror of
    :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`:
    block-diagonal over the mesh's boundary faces, ``B.apply(ψ)``
    returns a composite with **zero bulk** and, on each face's
    :attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.INFLOW_ROW`,
    the face's realized law applied to the outflow partial current —
    :math:`\mathcal{A}\,J^+` for the albedo family the
    :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
    produces (#290 P3). It composes as ``−B`` in ``(L + C − S − B)``;
    the loss's inflow rows then read :math:`J^- - \mathcal{A} J^+`
    (``L`` supplies the identity — the block-matrix table in
    :mod:`orpheus.diffusion.operators`).

    Parameters
    ----------
    mesh : MaterialMesh
        The mesh the trace (and the composite carrier) is bound to.
    face_laws : Mapping[str, LinearOperator]
        Per-face REALIZED boundary laws — the
        :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
        outputs, keyed by face name. Must cover exactly the mesh's
        boundary faces (a face without a law, or a law without a face,
        is an unrepresentable-by-construction error — the SN precedent
        where ``sn_mesh.bc`` and the trace layout share keys).
    """

    block_role = BlockRole.BOUNDARY

    def __init__(
        self,
        mesh: "MaterialMesh",
        face_laws: "Mapping[str, LinearOperator]",
    ) -> None:
        trace_faces = set(mesh.scalar_trace.face_names)
        given = set(face_laws)
        if given != trace_faces:
            raise ValueError(
                f"DiffusionBoundaryOperator: face_laws must cover exactly "
                f"the mesh's boundary faces {sorted(trace_faces)}; got "
                f"{sorted(given)} (missing: {sorted(trace_faces - given)}, "
                f"extra: {sorted(given - trace_faces)})."
            )
        self.mesh = mesh
        self.face_laws = dict(face_laws)

    @property
    def domain(self) -> "Optional[FunctionSpace]":
        # The composite carrier (NOT the bare trace) — matching the
        # L/C/S siblings for the OperatorSum composition guard (the
        # SNBoundaryOperator O.2b R5 precedent).
        return self.mesh.scalar_full_field_space

    @property
    def codomain(self) -> "Optional[FunctionSpace]":
        return self.mesh.scalar_full_field_space

    def apply(self, psi: "FullField", /) -> "FullField":
        r"""Return :math:`B\,\psi` — per-face :math:`\mathcal{A} J^+` on the
        inflow rows, zero everywhere else (zero bulk)."""
        trace = psi.boundary
        if not isinstance(trace, ScalarBoundaryFlux):
            raise TypeError(
                f"DiffusionBoundaryOperator: the input composite's "
                f"boundary must be a ScalarBoundaryFlux trace; got "
                f"{type(trace).__name__}."
            )
        if psi.bulk.mesh is not self.mesh:
            raise ValueError(
                "DiffusionBoundaryOperator.apply: input field and "
                "operator must share the same MaterialMesh instance "
                f"(mesh-identity invariant); got field mesh "
                f"{psi.bulk.mesh!r} vs operator mesh {self.mesh!r}."
            )
        out_boundary = ScalarBoundarySourceSink.zeros_on(self.mesh)
        for face, law in self.face_laws.items():
            out_boundary.face_view(face)[ScalarTraceSpace.INFLOW_ROW] = (
                law.apply(trace.outflow_view(face))
            )
        return FullField(
            bulk=ScalarSourceSink.zeros_on(self.mesh),
            boundary=out_boundary,
        )
