"""Model-independent numerical methods for reactor physics."""

from .eigenvalue import EigenvalueSolver, ProductionRateSolver, power_iteration
from .field import Field
from .functional import Functional, InnerProductFunctional, R_co, V_contra
from .iteration import KEigenvalue, SourceIteration
from .measure import (
    BundleMeasure,
    DiscreteMeasure,
    DiscreteMeasurePartition,
    equispaced,
    gauss_chebyshev,
    gauss_legendre,
)
from .operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    DiagonalOperator,
    IdentityOperator,
    LinearOperator,
    MissingCapability,
    OperatorProduct,
    OperatorSum,
    RankOneOperator,
    ScaledOperator,
    SumOfTensorProductsOperator,
    TensorProductOperator,
    ZeroOperator,
    outer,
)
from orpheus.numerics.quadrature import gauss_legendre_on_mu, lebedev_sphere, level_symmetric_sn, product_mu_phi
from .basis import Basis, SphericalHarmonicBasis
from .frame import FrameBase, GalerkinFrame, PetrovGalerkinFrame
from .projection import AnalysisOperator, ReconstructionOperator
from .space import DualSpace, FunctionSpace, TensorProductSpace
from .symmetry import SubgroupOfO3
from .vector import V, Vector

__all__ = [
    "AnalysisOperator",
    "Basis",
    "BundleMeasure",
    "CAP_APPLY",
    "CAP_APPLY_TRANSPOSE",
    "CAP_SOLVE",
    "DiagonalOperator",
    "DiscreteMeasure",
    "DiscreteMeasurePartition",
    "DualSpace",
    "EigenvalueSolver",
    "Field",
    "FrameBase",
    "Functional",
    "FunctionSpace",
    "InnerProductFunctional",
    "GalerkinFrame",
    "ReconstructionOperator",
    "IdentityOperator",
    "KEigenvalue",
    "LinearOperator",
    "MissingCapability",
    "OperatorProduct",
    "OperatorSum",
    "PetrovGalerkinFrame",
    "ProductionRateSolver",
    "R_co",
    "RankOneOperator",
    "ScaledOperator",
    "SourceIteration",
    "SphericalHarmonicBasis",
    "SubgroupOfO3",
    "SumOfTensorProductsOperator",
    "TensorProductOperator",
    "TensorProductSpace",
    "V",
    "V_contra",
    "Vector",
    "ZeroOperator",
    "equispaced",
    "gauss_chebyshev",
    "gauss_legendre",
    "gauss_legendre_on_mu",
    "lebedev_sphere",
    "level_symmetric_sn",
    "outer",
    "power_iteration",
    "product_mu_phi",
]
