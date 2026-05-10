"""Model-independent numerical methods for reactor physics."""

from .eigenvalue import EigenvalueSolver, power_iteration
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
    ScaledOperator,
    SumOfTensorProductsOperator,
    TensorProductOperator,
    ZeroOperator,
    as_scipy_linop,
)
from .quadrature import (
    gauss_legendre_on_mu,
    lebedev_sphere,
    level_symmetric_sn,
    product_mu_phi,
)
from .spherical_harmonics import evaluate_real_sh
from .symmetry import SubgroupOfO3

__all__ = [
    "BundleMeasure",
    "CAP_APPLY",
    "CAP_APPLY_TRANSPOSE",
    "CAP_SOLVE",
    "DiagonalOperator",
    "DiscreteMeasure",
    "DiscreteMeasurePartition",
    "EigenvalueSolver",
    "IdentityOperator",
    "KEigenvalue",
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "OperatorProduct",
    "OperatorSum",
    "ScaledOperator",
    "SourceIteration",
    "SubgroupOfO3",
    "SumOfTensorProductsOperator",
    "TensorProductOperator",
    "ZeroOperator",
    "as_scipy_linop",
    "equispaced",
    "evaluate_real_sh",
    "gauss_chebyshev",
    "gauss_legendre",
    "gauss_legendre_on_mu",
    "lebedev_sphere",
    "level_symmetric_sn",
    "power_iteration",
    "product_mu_phi",
]
