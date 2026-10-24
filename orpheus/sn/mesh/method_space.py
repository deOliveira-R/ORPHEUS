r"""SN method space: the quadrature and the trace ladder a realizer reads.

The :class:`SNMethodSpace` is the SN-specific argument to
:meth:`~orpheus.sn.boundary.realizer.SN_REALIZERS.realize`. It carries the two things an SN realizer needs to turn a declared
unit (a chain of the response network) into a
:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`:

* :attr:`quadrature` -- the angular quadrature (ordinate direction cosines and weights).
* :attr:`trace` -- the unified :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` precomputed for the
  quadrature and the face layout. The typed realized operators read every half-space of a face off it
  (``trace.inflow_space(face)``, ``trace.outflow_space(face)``).

The space names no face and holds no mesh. The face a realization acts on is the :class:`~orpheus.transport.boundary_conditions.FaceSite`
the realizer receives beside the space, so there is exactly one face value per realization. A method space that also carried a face
could disagree with the site without any check, which is ERR-041's own defect class (two independently supplied encodings of one
orientation); the field, the mesh the realizer chain never read, the index sets of ERR-041 and the optional trace all retired in
WF2.4. The ERR-041 invariant itself, that ``trace.inflow_space(face)`` holds exactly the ordinates with
:math:`\Omega\cdot n_f < 0`, is witnessed on the ladder
(``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``) and on the realizer
(``tests/gates/sn/operators/test_sn_boundary_realizer.py::TestVacuumRealizesOnTheTrace``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["SNMethodSpace"]


@dataclass(frozen=True)
class SNMethodSpace:
    r"""The method space of an SN face realization: the quadrature and the trace ladder over the problem's faces.

    See the module docstring for the field semantics. The canonical construction site is ``SNProblem.realize_boundary_role``,
    which builds the trace once per quadrature and layout and passes it through.
    """

    quadrature: "Quadrature"
    trace: "AngularTraceSpace"
