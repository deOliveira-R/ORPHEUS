"""The exterior conditions of a transport problem: declaration, binding, admission and realization tables.

The package replaces the string-keyed boundary laws of ``orpheus.geometry.boundary`` with ONE shared exterior type
(:class:`WallCondition`, an affine pair of a response network and an inflow), a binding of declared conditions to a
carrier's tagged faces (:func:`bind`), and a per-method realizer table whose single door answers both "is it admitted"
and "realize it". Each module's docstring holds its mathematics and its gotchas; this file only re-exports.

* :mod:`.network`: the factors, the closed chain sum, :class:`ResponseNetwork`;
* :mod:`.inflow`: the affine part :math:`q`;
* :mod:`.condition`, :mod:`.displacement`: the conditions and their cone arithmetic;
* :mod:`.errors`, :mod:`.declared`, :mod:`.binding`: the refusals, the twin of ``Materials``, the binding;
* :mod:`.admission`, :mod:`.realization`: the admission lattice and the realizer tables.
"""

from orpheus.transport.boundary_conditions.admission import Admission, Approximated, FaceSite, Realized, RefusalAxis, Refused, SiteClass, hollow_inner_refusal
from orpheus.transport.boundary_conditions.binding import Binding, BoundaryCarrier, bind
from orpheus.transport.boundary_conditions.condition import BoundaryCondition, WallCondition
from orpheus.transport.boundary_conditions.declared import BoundaryConditions
from orpheus.transport.boundary_conditions.displacement import WallDisplacement
from orpheus.transport.boundary_conditions.errors import (
    BoundaryError,
    BoundaryRefused,
    BoundaryResponseNotPositiveError,
    NonInvariantConditionError,
    SubmarkovViolationError,
    UnboundTagError,
    UnusedConditionError,
)
from orpheus.transport.boundary_conditions.inflow import FaceFrameInflow, InflowProfile, InflowSampler, InflowTerms, LabFrameInflow
from orpheus.transport.boundary_conditions.network import (
    CHAIN_ORDER,
    Chain,
    CurrentRetraction,
    DiffuseChain,
    Factor,
    IsotropicSection,
    MirrorPermutation,
    ResponseNetwork,
    SpecularChain,
    ZeroNetwork,
    chain_composes,
)
from orpheus.transport.boundary_conditions.realization import (
    DeckRealizer,
    RealizerTable,
    RefusingRealizer,
    ResponseRealizer,
    SealedTableError,
    SourceRealizer,
    StateCone,
    TableNotSealedError,
    sample_terms,
    stamp_boundary_role,
)

__all__ = [
    "CHAIN_ORDER",
    "Admission",
    "Approximated",
    "Binding",
    "BoundaryCarrier",
    "BoundaryCondition",
    "BoundaryConditions",
    "BoundaryError",
    "BoundaryRefused",
    "BoundaryResponseNotPositiveError",
    "Chain",
    "CurrentRetraction",
    "DeckRealizer",
    "DiffuseChain",
    "Factor",
    "FaceFrameInflow",
    "FaceSite",
    "InflowProfile",
    "InflowSampler",
    "InflowTerms",
    "IsotropicSection",
    "LabFrameInflow",
    "MirrorPermutation",
    "NonInvariantConditionError",
    "RealizerTable",
    "Realized",
    "RefusalAxis",
    "Refused",
    "ResponseNetwork",
    "RefusingRealizer",
    "ResponseRealizer",
    "SealedTableError",
    "SiteClass",
    "SourceRealizer",
    "SpecularChain",
    "StateCone",
    "SubmarkovViolationError",
    "TableNotSealedError",
    "UnboundTagError",
    "UnusedConditionError",
    "WallCondition",
    "WallDisplacement",
    "ZeroNetwork",
    "bind",
    "chain_composes",
    "hollow_inner_refusal",
    "sample_terms",
    "stamp_boundary_role",
]
