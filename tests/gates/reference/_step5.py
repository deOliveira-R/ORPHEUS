"""The step-5 gates' adapter and fixtures (#405 P2, ``.claude/plans/reference_p2_spec.md`` §1.5).

Every production name the step-5 gates reach is spelled ONCE here, resolved
on its module at call time (so a rebinding battery arm reaches every row,
lessons ``L100``), and a rename is one edit in this file.
"""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.citation import Citation
from orpheus.data.materials import Materials
from orpheus.numerics.mesh_free_function import RegionwiseConstant
from orpheus.numerics.question import Eigen, FixedSource
from tests.gates.specification import _fixtures as sf

PUBLISHED = "orpheus.reference.published"
CERTIFICATE = "orpheus.reference.certificate"
READING = "orpheus.reference.reading"
WITHDRAWAL = "orpheus.reference.withdrawal"
OBSERVABLE = "orpheus.numerics.observable"
ENCLOSURE = "orpheus.numerics.enclosure"
#: The one admission verb a publication and (step 6) the reader share (the orchestrator's ruling of
#: 2026-10-03: a module-level function, one exhaustive match over the observable).
ADMISSION = ("orpheus.specification.specification", "admit_observable")


def mod(name: str) -> Any:
    return importlib.import_module(name)


def c(module: str, name: str) -> Any:
    return getattr(mod(module), name)


# ── the observables and readings ────────────────────────────────────────────


def flux(values: Any = ((1.0, 0.5), (0.0, 2.0))) -> Any:
    return c(OBSERVABLE, "FluxIntegral")(RegionwiseConstant(np.array(values, float)))


def eigenvalue() -> Any:
    return c(OBSERVABLE, "Eigenvalue")()


def point(position: float = 0.5, group: int = 0) -> Any:
    return c(OBSERVABLE, "PointValue")(position, group)


def enclosure(value: float, bound: float) -> Any:
    return c(ENCLOSURE, "Enclosure")(value, bound)


def printed(text: str, locator: str = "Table 10", bibkey: str = "SoodForsterParsons2003") -> Any:
    return c(READING, "Printed")(text, Citation(bibkey, locator))


def withdrawal(reason: str = "an erratum moved the printed value", issue: int = 999) -> Any:
    return c(WITHDRAWAL, "Withdrawal")(reason, issue)


def current() -> Any:
    return c(PUBLISHED, "Current")()


# ── the specifications ──────────────────────────────────────────────────────


def eigen_medium() -> Any:
    """The infinite medium of the 2-group fuel under the k question (an eigen answer, no position)."""
    from orpheus.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, sf.fuel(), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def source_slab() -> Any:
    """A 2-region, 2-group slab under a fixed source (a field answer, no eigenvalue)."""
    from orpheus.specification import GeometrySpecification

    q = RegionwiseConstant(np.array([[1.0, 0.0], [0.0, 0.0]]))
    return GeometrySpecification(Materials({0: sf.fuel(), 1: sf.moderator()}), sf.slab2(), FixedSource(q))


def eigen_slab(region1: Any = None) -> Any:
    """The same slab under the k question (an eigen answer WITH a position); ``region1`` swaps the moderator."""
    from orpheus.specification import GeometrySpecification

    return GeometrySpecification(Materials({0: sf.fuel(), 1: sf.moderator() if region1 is None else region1}), sf.slab2(),
                                 Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


# ── the step-5 values ───────────────────────────────────────────────────────


def published(spec: Any = None, printed_map: Any = None, standing: Any = None) -> Any:
    spec = eigen_slab() if spec is None else spec
    printed_map = {eigenvalue(): printed("1.00000"), flux(): printed("0.6123", "Table 11")} if printed_map is None else printed_map
    return c(PUBLISHED, "PublishedSolution")(spec, printed_map, current() if standing is None else standing)


def exact(expression: Any = None, by: str = "derive_k_inf_rank_one") -> Any:
    """``Exact`` of a SymPy expression (default one third), stored as its ``srepr`` text."""
    import sympy

    expr = sympy.Rational(1, 3) if expression is None else sympy.sympify(expression)
    return c(CERTIFICATE, "Exact")(sympy.srepr(expr), by)


def derived(value: float = 1.0, bound: float = 1e-9, method: str = "Nystrom residual bound") -> Any:
    return c(CERTIFICATE, "DerivedBound")(enclosure(value, bound), method)


def claim(target: float, established: Any) -> Any:
    return c(CERTIFICATE, "Claim")(target, established)


def corroboration(anchor: Any, independence: str = "printed by an independent method (Sood 2003)") -> Any:
    return c(CERTIFICATE, "Corroboration")(anchor, independence)


def refinement(observable: Any, members: Any) -> Any:
    return c(CERTIFICATE, "Refinement")(observable, tuple(members))


def certificate(claims: Any, corroborations: Any = (), refinements: Any = (), standing: Any = None) -> Any:
    return c(CERTIFICATE, "ReferenceCertificate")(claims, tuple(corroborations), tuple(refinements),
                                                   current() if standing is None else standing)


def state_kind(cert: Any) -> str:
    """``Valid``, ``Invalid`` or ``Withdrawal``: the class name of the derived state."""
    return type(cert.state).__name__
