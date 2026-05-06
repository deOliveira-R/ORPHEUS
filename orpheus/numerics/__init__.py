"""Model-independent numerical methods for reactor physics."""

from .eigenvalue import EigenvalueSolver, power_iteration
from .measure import (
    BundleMeasure,
    DiscreteMeasure,
    equispaced,
    gauss_chebyshev,
    gauss_legendre,
)
from .operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    IdentityOperator,
    LinearOperator,
    LinearOperatorMixin,
    MissingCapability,
    OperatorProduct,
    OperatorSum,
    ScaledOperator,
    ZeroOperator,
    as_scipy_linop,
)

__all__ = [
    "BundleMeasure",
    "CAP_APPLY",
    "CAP_APPLY_TRANSPOSE",
    "CAP_SOLVE",
    "DiscreteMeasure",
    "EigenvalueSolver",
    "IdentityOperator",
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "OperatorProduct",
    "OperatorSum",
    "ScaledOperator",
    "ZeroOperator",
    "as_scipy_linop",
    "equispaced",
    "gauss_chebyshev",
    "gauss_legendre",
    "power_iteration",
]
