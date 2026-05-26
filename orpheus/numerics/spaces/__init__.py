r"""Function-space subclasses, organised by geometric / algebraic role.

The base :class:`~orpheus.numerics.space.FunctionSpace` lives in
``numerics/space.py`` (will move to ``numerics/spaces/function_space.py``
in P3.2 of the moment-space + layering plan). This sub-package houses
the specialised subclasses that carry domain-specific metadata beyond
``(name, shape, inner_product_weights)``:

* :class:`SphericalHarmonicSpace` (P1.2) — moment-space carrier for SH
  coefficients with the ``MomentMassMatrix`` diagonal already broadcast
  to the storage layout.
* Future: ``MeshFunctionSpace``, ``EnergyGroupSpace``,
  ``DiscreteAngularSpace`` per Grand Report v3 §5.3.

The existing :class:`~orpheus.numerics.trace_space.TraceSpace` is the
precedent for this sub-package pattern (frozen-dataclass inheritance
from ``FunctionSpace`` with ABC tag + per-face metadata). P3.2 absorbs
``trace_space.py`` into ``numerics/spaces/trace_space.py`` as part of
the broader reorganisation.

References
----------

* Grand Report v3 §5.3 — Space hierarchy.
* :mod:`orpheus.numerics.space` — the :class:`FunctionSpace` base.
* :mod:`orpheus.numerics.trace_space` — the inflow / outflow trace
  subclasses (precedent for this pattern).
"""

from __future__ import annotations

from orpheus.numerics.spaces.spherical_harmonic_space import (
    SphericalHarmonicSpace,
)

__all__ = ["SphericalHarmonicSpace"]
