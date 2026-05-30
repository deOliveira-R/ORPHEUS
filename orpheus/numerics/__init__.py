"""Model-independent numerical methods for reactor physics."""

from .eigenvalue import EigenvalueSolver, power_iteration
from .field import Field
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
    LinearOperatorMixin,
    MissingCapability,
    OperatorProduct,
    OperatorSum,
    RankOneOperator,
    ScaledOperator,
    SumOfTensorProductsOperator,
    TensorProductOperator,
    ZeroOperator,
    as_scipy_linop,
)
from orpheus.numerics.quadrature import gauss_legendre_on_mu, lebedev_sphere, level_symmetric_sn, product_mu_phi
from .basis.spherical_harmonic_basis import SphericalHarmonicBasis
from .projection import (
    GalerkinProjection,
    HarmonicMomentReconstruction,
    MomentProjection,
    PetrovGalerkinProjection,
    ProjectionOperator,
    ReconstructionOperator,
)
from .space import DualSpace, FunctionSpace, TensorProductSpace
from .symmetry import SubgroupOfO3

__all__ = [
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
    "FunctionSpace",
    "GalerkinProjection",
    "HarmonicMomentReconstruction",
    "MomentProjection",
    "ReconstructionOperator",
    "IdentityOperator",
    "KEigenvalue",
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "OperatorProduct",
    "OperatorSum",
    "PetrovGalerkinProjection",
    "ProjectionOperator",
    "RankOneOperator",
    "ScaledOperator",
    "SourceIteration",
    "SphericalHarmonicBasis",
    "SubgroupOfO3",
    "SumOfTensorProductsOperator",
    "TensorProductOperator",
    "TensorProductSpace",
    "ZeroOperator",
    "as_scipy_linop",
    "equispaced",
    "gauss_chebyshev",
    "gauss_legendre",
    "gauss_legendre_on_mu",
    "lebedev_sphere",
    "level_symmetric_sn",
    "power_iteration",
    "product_mu_phi",
]
