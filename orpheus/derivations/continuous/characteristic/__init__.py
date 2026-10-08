r"""The characteristic reference: transport integrated along the lines of a concentric body.

One reference posed from the specification (P1 of
``.claude/plans/characteristic_reference_architecture.md``), built on the
geometric kernel (:mod:`orpheus.geometry.chord`) and the reference
numerics kernel (:mod:`orpheus.derivations.common.dense_pencil`):

- :mod:`~orpheus.derivations.continuous.characteristic.walls` — each
  boundary point with what its law returns, read from the law's factors;
- :mod:`~orpheus.derivations.continuous.characteristic.closure` — the
  period of each line's unfolded path and the least solution of its cycle,
  the line part of the boundary resolvent, and the white walls' coupling,
  its diffuse part;
- :mod:`~orpheus.derivations.continuous.characteristic.basis` — the panel
  basis the emission density and the flux are represented in;
- :mod:`~orpheus.derivations.continuous.characteristic.transport` — the
  transport along each line on that basis: the traversals' source
  integrals, the Volterra block and the angular flux;
- :mod:`~orpheus.derivations.continuous.characteristic.lines` — the
  weighted sets of lines both test measures are built from: the impact,
  polar and cosine rules, their grading toward every tangency, each line's
  exact level, and the sources a line carries;
- :mod:`~orpheus.derivations.continuous.characteristic.assembly` — the
  Galerkin assembly over lines: one group's transport block, its line part
  and its diffuse walls' coupling, on a line rule graded from the group's
  optical scale;
- :mod:`~orpheus.derivations.continuous.characteristic.cross_sections` —
  each region's total cross section and emission matrices, and the
  regions each group emits in;
- :mod:`~orpheus.derivations.continuous.characteristic.system` — the
  multigroup Galerkin system on every group's block, its k pencil and its
  source pencil on the emission space: the fundamental and higher modes, the adjoint, the fixed
  source and the detector's adjoint flux;
- :mod:`~orpheus.derivations.continuous.characteristic.reading` — the
  reading at a point: the transported emission's scalar flux there, over
  the directions at the point as the line rule's own lines, and its angular
  flux;
- :mod:`~orpheus.derivations.continuous.characteristic.reference` — the
  door: the system posed from a specification and a resolution, the
  question it answers and the observables it reads;
- :mod:`~orpheus.derivations.continuous.characteristic.grading` — the
  geometric, exponential and hp gradings every rule of the package places
  its piece ends by.
"""

from .assembly import GroupTransport, LineRule, TransportResolution
from .basis import PanelBasis
from .closure import LinePeriod, TrappedSource, WallCoupling
from .reference import CharacteristicDerivation, Resolution, characteristic_reference
from .cross_sections import RegionCrossSections
from .system import EmissionSpace, GalerkinSystem
from .transport import TraversalRule
from .walls import Wall, Walls

__all__ = [
    "CharacteristicDerivation",
    "EmissionSpace",
    "GalerkinSystem",
    "GroupTransport",
    "LinePeriod",
    "LineRule",
    "PanelBasis",
    "RegionCrossSections",
    "Resolution",
    "TransportResolution",
    "TrappedSource",
    "TraversalRule",
    "Wall",
    "WallCoupling",
    "Walls",
    "characteristic_reference",
]
