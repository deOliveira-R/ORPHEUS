r"""Discrete spectral bases on measure spaces (Grand Report v3 §5.4).

A *Basis* is a (possibly truncated) collection of functions on a measure
space that admits tabulation at sampling points and a discrete Gram
(mass) matrix against a quadrature rule. Derived operators — projection
:math:`M = Y^{*} W`, reconstruction :math:`R = (2\ell+1) Y` — are
constructed from the basis + a measure and live in
:mod:`orpheus.numerics.projection`.

The Basis family will grow — §5.4 of the architecture report lists nine
basis types (real spherical harmonics, Chebyshev, Lagrange polynomials,
finite-element shape functions, ...). This package houses them.
:class:`SphericalHarmonicBasis` (angular, on :math:`S^2`) is the first
member; :class:`IndicatorBasis` (piecewise-constant / P0, on a
tensor-product cell partition) is the second — the trial side of the
spatial / energy homogenisation Frame.

The :class:`~orpheus.numerics.basis.base.Basis` ABC (formerly deferred
"until a second concrete basis arrives") is now promoted — the forcing
consumer, the generic :class:`~orpheus.numerics.frame.FrameBase`, binds an
*abstract* basis to a measure, so it needs the interface rather than a
second instance. The contract is the three fundamental operations
:meth:`~orpheus.numerics.basis.base.Basis.evaluate` (tabulate),
:meth:`~orpheus.numerics.basis.base.Basis.synthesize` (naked :math:`S_0`),
and :meth:`~orpheus.numerics.basis.base.Basis.mass_matrix` (discrete Gram).
See :mod:`orpheus.numerics.basis.base` for the rationale.

References
----------

* Grand Report v3 §5.4 — Basis hierarchy.
* Grand Report v3 §19 — Harmonic projection.
"""

from __future__ import annotations

from orpheus.numerics.basis.base import Basis, GramStructure, TruncatedBasis
from orpheus.numerics.basis.descent import Descent
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.basis.overlap_basis import OverlapBasis
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.basis.weighted_indicator_basis import WeightedIndicatorBasis

__all__ = [
    "Basis",
    "Descent",
    "GramStructure",
    "IndicatorBasis",
    "LegendreBasis",
    "OverlapBasis",
    "SphericalHarmonicBasis",
    "TruncatedBasis",
    "WeightedIndicatorBasis",
]
