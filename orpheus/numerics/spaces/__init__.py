r"""Function-space subclasses, organised by geometric / algebraic role.

The base :class:`~orpheus.numerics.space.FunctionSpace` lives in
``numerics/space.py``. This sub-package houses the specialised
subclasses that carry domain-specific metadata beyond
``(name, shape, inner_product_weights, units)``:

* :class:`SphericalHarmonicSpace` (P1.2) — moment-space carrier for SH
  coefficients with the ``MomentMassMatrix`` diagonal already broadcast
  to the storage layout.
* :class:`TraceSpace` (Depth B D-C) — ABC for boundary-trace function
  spaces. Concrete subclasses :class:`InflowTraceSpace` and
  :class:`OutflowTraceSpace` carry per-face inflow / outflow masks.
* Future: ``MeshFunctionSpace``, ``EnergyGroupSpace``,
  ``DiscreteAngularSpace`` per Grand Report v3 §5.3.

References
----------

* Grand Report v3 §5.3 — Space hierarchy.
* :mod:`orpheus.numerics.space` — the :class:`FunctionSpace` base.
"""

from __future__ import annotations

from orpheus.numerics.spaces.spherical_harmonic_space import (
    SphericalHarmonicSpace,
)
from orpheus.numerics.spaces.trace_space import (
    InflowTraceSpace,
    OutflowTraceSpace,
    TraceSpace,
)

__all__ = [
    "InflowTraceSpace",
    "OutflowTraceSpace",
    "SphericalHarmonicSpace",
    "TraceSpace",
]
