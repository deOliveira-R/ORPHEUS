"""Model-independent numerical methods for reactor physics."""

from .eigenvalue import EigenvalueSolver, power_iteration
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
    "EigenvalueSolver",
    "power_iteration",
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "OperatorSum",
    "OperatorProduct",
    "ScaledOperator",
    "IdentityOperator",
    "ZeroOperator",
    "as_scipy_linop",
    "CAP_APPLY",
    "CAP_SOLVE",
    "CAP_APPLY_TRANSPOSE",
]
