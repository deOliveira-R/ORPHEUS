r"""The :class:`Basis` ABC — a discrete spectral basis on a measure space.

A *basis* is the **synthesis (trial) side** of a discrete frame
(:class:`orpheus.numerics.frame.Frame`): a collection of functions that can be
tabulated at sample points, reconstructed from coefficients, and weighed against
a quadrature to form a discrete Gram. It is the choice-free, measure-free half of
the analysis/synthesis pair — the :class:`~orpheus.numerics.measure.DiscreteMeasure`
supplies the analysis (test) side, and the :class:`Frame` binds the two.

Why the ABC is un-deferred on a single concrete basis
=====================================================

The package previously deferred a formal ``Basis`` ABC "until a second concrete
basis arrives" (``feedback_unify_after_two_instances``). It is promoted now — with
:class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
still the only concrete member — because the **forcing consumer has arrived**: the
generic :class:`~orpheus.numerics.frame.Frame` binds an *abstract* basis to a
measure and applies a non-identity morphism (analysis / synthesis) across the
interface. The Frame needs the interface, not a second basis; the interface itself
is math-rigid (Grand Report v3 §5.4 lists the nine eventual bases — real spherical
harmonics, Chebyshev, Lagrange, FE shape functions, ...) and every one of them
tabulates, reconstructs, and has a discrete Gram. Deferring the ABC would force the
Frame to name a concrete basis, which is exactly the coupling the abstraction exists
to break.

The contract (three fundamental operations)
===========================================

* :meth:`evaluate` — tabulate the basis functions at a set of sample points
  (the :math:`\Phi(\text{point}, \text{mode})` table).
* :meth:`synthesize` — the naked synthesis :math:`S_0(c) = \sum_k \phi_k\, c_k`
  (NO measure weights, NO dual-frame factor) — the choice-free reconstruction.
* :meth:`mass_matrix` — the discrete Gram :math:`\sum_n w_n\, \phi_j(x_n)\,
  \phi_k(x_n)` against a quadrature; equals the continuous Gram when the rule is
  exact to the basis's degree.

The *weighted* operations — the analysis :math:`M = \sum_n w_n\,\phi_k(x_n)\,(\cdot)`
and the canonical-dual reconstruction :math:`R` — bind the basis to a specific
measure (the choice), so they live on the :class:`Frame`, not here.

References
----------

* Grand Report v3 §5.4 — Basis hierarchy.
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed.
  Birkhäuser — the analysis/synthesis operator pair this ABC is the trial side of.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from numpy.typing import NDArray

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure


__all__ = ["Basis"]


class Basis(ABC):
    r"""Abstract discrete spectral basis — the synthesis side of a :class:`Frame`.

    Concrete bases (real spherical harmonics today; Legendre, Chebyshev,
    Lagrange/FE shape functions to come) implement the three operations below.
    A basis is *choice-free*: it knows its functions and their convention, but
    not which quadrature samples them — that choice is the
    :class:`~orpheus.numerics.measure.DiscreteMeasure`, bound to the basis by a
    :class:`~orpheus.numerics.frame.Frame`.
    """

    @abstractmethod
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""Tabulate the basis functions at ``points``.

        Returns the :math:`\Phi(\text{point}, \text{mode})` table — for the
        spherical-harmonic basis, :math:`Y[n, \ell, \ell+m] = Y_\ell^m(\hat\Omega_n)`
        at the ``(N, 3)`` direction cosines; a flat-mode basis returns ``(N, K)``.

        (Positional-only: the interface is by-position, so concrete bases may
        name the argument in their own domain vocabulary — e.g. ``directions``.)
        """
        ...

    @abstractmethod
    def synthesize(self, coefficients: NDArray, points: NDArray, /) -> NDArray:
        r"""Naked synthesis :math:`S_0(c) = \sum_k \phi_k(x)\, c_k` at ``points``.

        The pure (frame-theory *synthesis operator* :math:`T^*`) reconstruction
        with NO measure weights and NO dual-frame factor. The canonical-dual
        reconstruction :math:`R` (which DOES carry the dual factor) is a
        :class:`Frame` face built on top of this.
        """
        ...

    @abstractmethod
    def mass_matrix(self, measure: "DiscreteMeasure", /) -> NDArray:
        r"""Discrete Gram :math:`\sum_n w_n\, \phi_j(x_n)\, \phi_k(x_n)` over ``measure``.

        The frame operator :math:`S = T^* T` in discrete form. Equals the
        continuous Gram (the basis's intrinsic metric) when the quadrature is
        exact to the basis's degree; the residual is a quadrature-exactness
        diagnostic.
        """
        ...
