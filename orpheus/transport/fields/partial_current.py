r"""Partial-current boundary field :math:`(J^+, J^-)` on a material-mesh boundary.

The FLUX-role leaf of the SCALAR trace family (#290 P2) — the boundary
sibling of :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
exactly as the ℓ=0 moment structure dictates: where ``ScalarFlux`` is the
full-range moment :math:`\phi = \int_{4\pi}\psi\,d\Omega` of the bulk
angular flux, the partial currents are the two HALF-range first moments of
the boundary trace,

.. math::

   J^\pm_g \;=\; \int_{\Omega\cdot\hat n \gtrless 0}
                 |\Omega\cdot\hat n|\,\psi_g\,d\Omega ,

face-local against the OUTWARD normal: :math:`J^+` leaves the domain,
:math:`J^-` enters — at every face, ``xmin`` included. Under the P1
closure they are exactly the Cauchy data of the diffusion operator:

.. math::

   \phi_\Gamma = 2\,(J^+ + J^-), \qquad J = J^+ - J^-,
   \qquad J^\pm = \tfrac{\phi_\Gamma}{4} \pm \tfrac{J}{2},

(Bell & Glasstone 1970, the Marshak boundary treatment; Duderstadt &
Hamilton 1976). One input DOF (:math:`J^-`) and one output DOF
(:math:`J^+`) per face per group make an inconsistent Cauchy pair
unrepresentable, and every boundary law a member of the albedo family
:math:`J^- = \mathcal{A}\,J^+` — vacuum :math:`\mathcal{A}=0` (zero
incoming current, the Marshak realization), reflective
:math:`\mathcal{A}=1`, albedo :math:`\mathcal{A}=\alpha`, and the
zero-flux Dirichlet idealization :math:`\mathcal{A}=-1`
(:math:`\phi_\Gamma = 0`; positivity of :math:`J^\pm` is a property of
the PHYSICAL laws :math:`\mathcal{A}\in[0,1]`, not a type invariant).
The full convention contract is ``.claude/plans/diffusion_crosswalk.md``.

Storage follows the :class:`~orpheus.transport.fields._bases.TraceField`
locus discipline: ONE flat backing buffer on the mesh's cached
:class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
(:attr:`MaterialMesh.scalar_trace
<orpheus.transport.mesh.material_mesh.MaterialMesh.scalar_trace>`), face
slots of shape ``(2, ng, *face_spatial)`` with the component-row
convention owned by
:attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.OUTFLOW_ROW`
/ :attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.INFLOW_ROW`.

Units: :data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS`
(:math:`1/(\mathrm{cm^2\,s})`) — a partial current is an areal rate, the
same dimensional family as the scalar flux (eV-free, group-integrated by
construction).

DSA seam (#2): the SN→diffusion boundary restriction is the half-range
ℓ=0 moment of the angular ``BoundaryFlux`` under the SN trace metric
:math:`G_s = |\Omega\cdot\hat n|\odot w_n`; that reduction is owned by
the future moment-restriction operator — this leaf just carries its
codomain.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import ScalarTraceField
from orpheus.transport.fields._flux_role import FluxRole


__all__ = ["PartialCurrent"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class PartialCurrent(FluxRole, ScalarTraceField):
    r"""Partial-current pair :math:`(J^+, J^-)` per boundary face per group.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer of shape ``(layout.total_size,)``.
    space : ScalarTraceSpace
        The mesh's cached scalar trace
        (:attr:`MaterialMesh.scalar_trace
        <orpheus.transport.mesh.material_mesh.MaterialMesh.scalar_trace>`);
        construction via :meth:`zeros_on` / :meth:`from_mesh` /
        :meth:`from_face_arrays` is the canonical path.
    mesh : MaterialMesh
        The method-agnostic mesh+materials carrier the trace belongs to
        (inherited annotation — this family does NOT narrow it; the
        scalar trace is meaningful on any material mesh, including an
        :class:`SNMesh` when DSA restricts onto it).

    Notes
    -----
    Algebra (same-class ``+``/``-``, scalar ``*``/``/``, mesh-identity
    and layout guards) is inherited from :class:`TraceField`. The
    named-quantity accessors below own the P1 dictionary — no consumer
    re-derives it (Pattern 3 / the crosswalk contract).
    """

    #: Dimensional identity: areal rate ``1/(cm²·s)`` — the same family
    #: as ScalarFlux (a current IS an angle-weighted flux moment).
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS

    # ── Named component views (the P1 dictionary lives HERE) ──────────

    def outflow_view(self, face: str) -> NDArray:
        r"""``J⁺`` on ``face`` — shape ``(ng, *face_spatial)``, a VIEW.

        The outgoing partial current through the face's outward normal.
        Writes propagate to the backing buffer (the operator-side write
        target for the leakage operator's boundary emission).
        """
        return self.face_view(face)[ScalarTraceSpace.OUTFLOW_ROW]

    def inflow_view(self, face: str) -> NDArray:
        r"""``J⁻`` on ``face`` — shape ``(ng, *face_spatial)``, a VIEW.

        The incoming partial current (the BC datum: the albedo family
        sets ``J⁻ = 𝒜·J⁺``).
        """
        return self.face_view(face)[ScalarTraceSpace.INFLOW_ROW]

    def net_current(self, face: str) -> NDArray:
        r"""Net OUTWARD current ``J = J⁺ − J⁻`` on ``face`` (a copy).

        Positive = leakage out of the domain, at every face (face-local
        outward-normal convention; conversion to an axis-signed vector
        current happens at the consumer that needs it — the crosswalk's
        orientation row).
        """
        return np.asarray(self.outflow_view(face) - self.inflow_view(face))

    def boundary_scalar_flux(self, face: str) -> NDArray:
        r"""P1 boundary scalar flux ``φ_Γ = 2(J⁺ + J⁻)`` on ``face`` (a copy).

        The other half of the Cauchy pair, reconstructed through the P1
        closure — the single site of that dictionary (crosswalk).
        """
        return np.asarray(2.0 * (self.outflow_view(face) + self.inflow_view(face)))
