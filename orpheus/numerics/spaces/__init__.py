r"""Function-space subclasses, organised by geometric / algebraic role.

The base :class:`~orpheus.numerics.space.FunctionSpace` lives in
``numerics/space.py``. This sub-package houses the specialised
subclasses that carry domain-specific metadata beyond
``(name, shape, inner_product_weights)``:

* :class:`SphericalHarmonicSpace` (P1.2) — moment-space carrier for SH
  coefficients with the ``MomentMassMatrix`` diagonal already broadcast
  to the storage layout.
* :class:`TraceSpace` (#205 / #201 unification) — the single
  whole-boundary trace function space. Carries the :class:`FaceLayout`
  + the signed :math:`\Omega\cdot\hat n` per face; inflow / outflow are
  *selectors* over it (no longer separate ``Inflow`` / ``Outflow``
  subclasses).
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
from orpheus.numerics.spaces.trace_space import TraceSpace

__all__ = [
    "SphericalHarmonicSpace",
    "TraceSpace",
]
