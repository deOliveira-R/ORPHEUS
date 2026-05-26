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
:class:`SphericalHarmonicBasis` is the first member.

A formal ``Basis`` ABC is deferred until a second concrete basis arrives
— see the project's ``feedback_unify_after_two_instances`` rule. For
now, every basis is expected to expose the methods below, but the
contract is enforced by tests / consumers, not by an abstract base.

Common method shape (informal)
==============================

.. code-block:: python

    class SomeBasis:
        # Tabulate basis functions at sampling points
        def evaluate(self, points) -> NDArray: ...

        # Discrete Gram matrix against a quadrature; approximates the
        # theoretical metric when the quadrature integrates degree
        # 2L exactly
        def discrete_mass_matrix(self, measure) -> NDArray: ...

References
----------

* Grand Report v3 §5.4 — Basis hierarchy.
* Grand Report v3 §19 — Harmonic projection.
"""

from __future__ import annotations

from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis

__all__ = ["SphericalHarmonicBasis"]
