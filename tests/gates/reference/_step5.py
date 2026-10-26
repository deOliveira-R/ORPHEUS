"""The step-5 gates' adapter and fixtures (#405 P2, ``.claude/plans/reference_p2_spec.md`` §1.5).

Every production module the step-5 gates reach is imported ONCE here, as a
module, and each name is read off it at call time (``OBSERVABLE.Ratio``), so a
battery arm rebinding the module's attribute reaches every row (lessons
``L100``) and a module rename is one edit in this file. The imports are
statements, so the test-layer gate sees each edge (layer-DAG ruling 27,
2026-10-26: they were ``importlib.import_module`` of a string until then). The
module constants are upper-case: each is the module object a row reads a name
from, and the lower-case spellings are this file's builders (``enclosure``,
``certificate``, ``published``, ``withdrawal``).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.citation import Citation
from orpheus.data.materials import Materials
from orpheus.numerics.mesh_free_function import RegionwiseConstant
from orpheus.numerics.question import Eigen, FixedSource
from orpheus.numerics import enclosure as ENCLOSURE
from orpheus.numerics import observable as OBSERVABLE
from orpheus.numerics import outcome as OUTCOME
from orpheus.reference import certificate as CERTIFICATE
from orpheus.reference import published as PUBLISHED
from orpheus.reference import reading as READING
from orpheus.reference import withdrawal as WITHDRAWAL
from orpheus.specification import specification as SPECIFICATION
from tests.gates.specification import _fixtures as sf

#: The one admission verb a publication and (step 6) the reader share (the orchestrator's ruling of
#: 2026-10-03: a module-level function, one exhaustive match over the observable): its module and the
#: name every binding of it carries (a route gate walks ``sys.modules`` for that name).
ADMISSION = (SPECIFICATION, "admit_observable")


def untyped(value: Any) -> Any:
    """A value the type checker does not narrow (a row testing a RUNTIME refusal of a wrong type)."""
    return value


# ── the observables and readings ────────────────────────────────────────────


def flux(values: Any = ((1.0, 0.5), (0.0, 2.0))) -> Any:
    return OBSERVABLE.FluxIntegral(RegionwiseConstant(np.array(values, float)))


def eigenvalue() -> Any:
    return OBSERVABLE.Eigenvalue()


def point(position: float = 0.5, group: int = 0) -> Any:
    return OBSERVABLE.PointValue(position, group)


def enclosure(value: float, bound: float) -> Any:
    return ENCLOSURE.Enclosure(value, bound)


def printed(text: str, locator: str = "Table 10", bibkey: str = "SoodForsterParsons2003") -> Any:
    return READING.Printed(text, Citation(bibkey, locator))


def withdrawal(reason: str = "the publication was found unfit as a whole", issue: int = 999) -> Any:
    return WITHDRAWAL.Withdrawal(reason, issue)


def current() -> Any:
    return PUBLISHED.Current()


# ── the specifications ──────────────────────────────────────────────────────


def eigen_medium() -> Any:
    """The infinite medium of the 2-group fuel under the k question (an eigen answer, no position)."""
    from orpheus.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, sf.fuel(), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def source_slab() -> Any:
    """A 2-region, 2-group slab under a fixed source (a field answer, no eigenvalue)."""
    from orpheus.specification import GeometrySpecification

    q = RegionwiseConstant(np.array([[1.0, 0.0], [0.0, 0.0]]))
    return GeometrySpecification(Materials({0: sf.fuel(), 1: sf.moderator()}), sf.slab2(), FixedSource(q),
                                 exterior=sf.vacuum_exterior(sf.slab2()))


def eigen_slab(region1: Any = None) -> Any:
    """The same slab under the k question (an eigen answer WITH a position); ``region1`` swaps the moderator."""
    from orpheus.specification import GeometrySpecification

    return GeometrySpecification(Materials({0: sf.fuel(), 1: sf.moderator() if region1 is None else region1}), sf.slab2(),
                                 Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)), exterior=sf.vacuum_exterior(sf.slab2()))


# ── the step-5 values ───────────────────────────────────────────────────────


def published(spec: Any = None, printed_map: Any = None, standing: Any = None) -> Any:
    spec = eigen_slab() if spec is None else spec
    printed_map = {eigenvalue(): printed("1.00000"), flux(): printed("0.6123", "Table 11")} if printed_map is None else printed_map
    return PUBLISHED.PublishedSolution(spec, printed_map, current() if standing is None else standing)


def exact(expression: Any = None, by: str = "derive_k_inf_rank_one") -> Any:
    """``Exact`` of a SymPy expression (default one third), stored as its ``srepr`` text."""
    import sympy

    expr = sympy.Rational(1, 3) if expression is None else sympy.sympify(expression)
    return CERTIFICATE.Exact(sympy.srepr(expr), by)


def derived(value: float = 1.0, bound: float = 1e-9, method: str = "Nystrom residual bound") -> Any:
    return CERTIFICATE.DerivedBound(enclosure(value, bound), method)


def claim(target: float, established: Any) -> Any:
    return CERTIFICATE.Claim(target, established)


def corroboration(anchor: Any, independence: str = "printed by an independent method (Sood 2003)") -> Any:
    return CERTIFICATE.Corroboration(anchor, independence)


def refinement(observable: Any, members: Any) -> Any:
    return CERTIFICATE.Refinement(observable, tuple(members))


def certificate(claims: Any, corroborations: Any = (), refinements: Any = (), standing: Any = None) -> Any:
    return CERTIFICATE.ReferenceCertificate(claims, tuple(corroborations), tuple(refinements),
                                                   current() if standing is None else standing)


def state_kind(cert: Any) -> str:
    """``Valid``, ``Invalid`` or ``Withdrawal``: the class name of the derived state."""
    return type(cert.state).__name__
