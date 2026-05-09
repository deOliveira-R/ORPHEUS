r"""Tensor-decomposed boundary conditions for transport solvers.

A boundary condition on a transport interface is, in full generality,
a linear operator :math:`R` mapping the *outgoing* angular flux at the
boundary to the *incoming* angular flux:

.. math::

    \psi_{\text{in}}(\Omega)
    = (R \psi_{\text{out}})(\Omega)
    = \sum_\alpha \bigl(G_\alpha \, \psi_{\text{out}}\bigr)(\Omega)
      \cdot A_\alpha,
    :label: bc-tensor-decomposition

where each term :math:`G_\alpha \otimes A_\alpha` is the tensor product
of a **geometric operator** :math:`G_\alpha` (a permutation, a
pushforward, an angular average, a spatial wrap-around, …) with a
**scalar amplitude** :math:`A_\alpha` (typically an albedo
:math:`\in [0, 1]`).

Most boundary conditions of practical interest are **rank-1**:

* :class:`VacuumBoundaryOperator` — :math:`R = 0` (the empty sum, rank 0; algebraically
  the trivial case of the decomposition).
* :class:`SpecularBoundaryOperator` — :math:`R = G_{\text{refl}} \cdot \alpha` where
  :math:`G_{\text{refl}}` is the angular-permutation operator that maps
  ordinate :math:`\Omega_n` to its reflected partner
  :math:`(\Omega_n \cdot \hat{n})` and :math:`\alpha \in [0, 1]` is the
  specular albedo. Equivalent to a
  :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` under the
  reflection map.
* :class:`WhiteBoundaryOperator` — :math:`R = G_{\text{diff}} \cdot \alpha` where
  :math:`G_{\text{diff}}` is the cosine-weighted angular average over the
  outgoing hemisphere, broadcast isotropically over the incoming
  hemisphere (Lambertian reflection). Rank-1 in *angle* even though the
  geometric operator is an integral, not a permutation.
* :class:`PeriodicBoundaryOperator` — :math:`R` is a *spatial* pushforward (wrap-around
  to the opposite face) with :math:`\alpha = 1`. Rank-1 in space; the
  angular structure is identity.
* :class:`AlbedoBoundaryOperator` — :math:`R = I \cdot \alpha` where :math:`I` is the
  angular identity. Rank-1; the geometric operator is trivial.

Mixed and partial-current boundaries are **rank-N** sums of the above
primitives:

* :class:`MixedBoundaryOperator` — a list of ``(weight, BoundaryOperator)`` pairs whose
  ``apply_to_incoming`` is the linear combination of the components.
  ``MixedBoundaryOperator([(0.3, SpecularBoundaryOperator), (0.7, WhiteBoundaryOperator)])`` realises the standard
  Marshak mixed boundary (Bell & Glasstone 1970, §1.5).

The protocol :class:`BoundaryOperator` is what production solvers consume.
Each solver's ``BC_REGISTRY`` factory returns a concrete
:class:`BoundaryOperator` instance from a solver-agnostic
:class:`~orpheus.geometry.mesh.BC` declaration; the sweep then calls
``resolved.apply_to_incoming(angular_flux_outgoing, quadrature)`` with
no string-kind branching at the call site.

The decomposition framing pays off architecturally because
multi-region interfaces, partial-current couplings, and asymmetric
left/right boundaries are all instances of the *same* algebra: pick
your geometric operator, pick your amplitude, sum. White, periodic,
and albedo BCs are added in this module as primitives even though
they are not yet wired into :func:`~orpheus.sn.solver.solve_sn` —
they ship as building blocks that downstream waves of the SN reshape
campaign can plumb into the sweep without further protocol work.

References
----------

* Lewis, E. E. & Miller, W. F. *Computational Methods of Neutron
  Transport*, §3.4 (boundary conditions in transport).
* Bell, G. I. & Glasstone, S. *Nuclear Reactor Theory*, §1.5
  (albedo, white, and Marshak boundary conditions).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol, Sequence, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


# ═══════════════════════════════════════════════════════════════════════
# Protocol
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class BoundaryOperator(Protocol):
    r"""Protocol for a tensor-decomposed boundary condition.

    A :class:`BoundaryOperator` is the runtime representation of the operator

    .. math::

        R = \sum_\alpha G_\alpha \otimes A_\alpha

    that consumes the outgoing angular flux at a face and produces the
    incoming angular flux. Concrete implementations may be **rank-1**
    (a single :math:`G \otimes A` term: vacuum, specular, white,
    periodic, albedo) or **rank-N** (:class:`MixedBoundaryOperator` — a sum of the
    rank-1 primitives).

    The :class:`~orpheus.sn.quadrature.AngularQuadrature` argument lets
    each concrete BC reach into the quadrature's structural metadata
    (``reflection_index(axis)``, weights, level structure) without the
    Protocol exposing those as required attributes — the responsibility
    for being able to query the quadrature lives in the BC, not in
    every consumer.
    """

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        r"""Compute the incoming angular flux from the outgoing.

        Parameters
        ----------
        angular_flux_outgoing : np.ndarray
            Angular flux at the boundary face, indexed over all
            ordinates of ``quadrature``. Shape ``(N_ord, ...)`` where
            the trailing axes are typically energy groups.
        quadrature : AngularQuadrature
            The angular quadrature; lets the BC query reflection
            partners, weights, and level structure.

        Returns
        -------
        np.ndarray
            Incoming angular flux at the boundary face, same shape as
            ``angular_flux_outgoing``. For ordinates whose direction
            cosine is *outgoing* at this face (and thus do not have an
            incoming value), the returned entries are zero — sweeps
            consume only the incoming entries.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# Concrete BCs
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class VacuumBoundaryOperator:
    r"""Vacuum boundary: :math:`R = 0`.

    The empty sum in the tensor decomposition: no incoming flux,
    irrespective of what leaks out. Algebraically a rank-0 case of
    :eq:`bc-tensor-decomposition`.
    """

    #: String tag for legacy string-kind comparisons. The Wave B
    #: refactor preserves the existing BC-resolution test contract
    #: (``sn_mesh.bc_right == "vacuum"`` continues to evaluate True)
    #: while consumers transition to direct
    #: :meth:`apply_to_incoming` calls. Wave C/D may drop this.
    kind: str = "vacuum"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, VacuumBoundaryOperator):
            return True
        return NotImplemented

    def __hash__(self) -> int:
        return hash(("VacuumBoundaryOperator",))

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        return np.zeros_like(angular_flux_outgoing)


@dataclass(frozen=True)
class SpecularBoundaryOperator:
    r"""Specular reflection with optional albedo.

    Tensor decomposition :math:`(G_{\text{refl}}, \alpha)` where
    :math:`G_{\text{refl}}` is the index permutation realising
    :math:`\Omega \mapsto \Omega - 2(\Omega \cdot \hat{n}) \hat{n}` on
    the quadrature ordinates and :math:`\alpha \in [0, 1]` is the
    specular albedo (1 = perfect reflection; the standard ``BC.reflective``
    case). The same operator is the
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` of
    the angular measure under the reflection map, with the Jacobian
    convention ``|R| = 1`` since reflections are isometries.

    Parameters
    ----------
    axis : str
        Axis of reflection: ``"x"``, ``"y"``, or ``"z"``. The
        :meth:`~orpheus.sn.quadrature.AngularQuadrature.reflection_index`
        method maps each ordinate to its reflected partner under this
        axis.
    albedo : float
        Specular albedo. Defaults to 1 (perfect reflection).
    """

    axis: str = "x"
    albedo: float = 1.0

    #: String tag for legacy string-kind comparisons. The default
    #: ``albedo == 1.0`` SpecularBoundaryOperator (the standard ``BC.reflective``
    #: case) compares equal to the string ``"reflective"``; tagged
    #: SpecularBoundaryOperator instances with ``albedo != 1`` compare equal to
    #: ``"partial"`` instead.
    @property
    def kind(self) -> str:
        return "reflective" if self.albedo == 1.0 else "partial"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, SpecularBoundaryOperator):
            return (
                self.axis == other.axis and self.albedo == other.albedo
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash(("SpecularBoundaryOperator", self.axis, self.albedo))

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        ref = quadrature.reflection_index(self.axis)
        return self.albedo * angular_flux_outgoing[ref]


@dataclass(frozen=True)
class WhiteBoundaryOperator:
    r"""White (Lambertian) boundary with optional albedo.

    Tensor decomposition :math:`(G_{\text{diff}}, \alpha)` where
    :math:`G_{\text{diff}}` is the cosine-weighted angular average
    over the outgoing hemisphere, broadcast isotropically over the
    incoming hemisphere:

    .. math::

        \psi_{\text{in}}(\Omega) = \frac{\alpha}{\pi}
            \sum_{\Omega' :\, \Omega' \cdot \hat{n} > 0}
            w(\Omega')\,|\Omega' \cdot \hat{n}|\,
            \psi_{\text{out}}(\Omega'),

    where the result is independent of the incoming :math:`\Omega`
    direction (Lambertian / cosine emission). The factor :math:`1/\pi`
    is the canonical Lambertian normalisation used in radiative
    transfer; the implementation here normalises by the outgoing
    cosine-weighted weight sum so the BC is conservative for any
    quadrature: the total returned current equals the incoming
    current (times :math:`\alpha`), which is the property the consumer
    actually needs — see Bell & Glasstone 1970 §1.5.

    Parameters
    ----------
    axis : str
        Boundary normal axis: ``"x"``, ``"y"``, or ``"z"``.
    outward_sign : int
        Sign of the outward normal along ``axis``: ``+1`` for the
        upper face (``xmax`` / ``ymax``) and ``-1`` for the lower face
        (``xmin`` / ``ymin``). Selects which ordinates are *outgoing*
        at this face.
    albedo : float
        Diffuse albedo. Defaults to 1 (perfectly reflecting).
    """

    axis: str = "x"
    outward_sign: int = +1
    albedo: float = 1.0

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Direction cosine along the boundary normal axis.
        if self.axis == "x":
            mu_n = quadrature.mu_x
        elif self.axis == "y":
            mu_n = quadrature.mu_y
        elif self.axis == "z":
            mu_n = getattr(quadrature, "mu_z", None)
            if mu_n is None:
                raise ValueError(
                    "WhiteBoundaryOperator(axis='z') requires a quadrature with mu_z "
                    "(2-D / 3-D adapters: Lebedev, level-symmetric, "
                    "product). The 1-D Gauss-Legendre adapter has no "
                    "mu_z attribute."
                )
        else:
            raise ValueError(f"Unknown axis: {self.axis!r}")

        weights = quadrature.weights
        # Outgoing ordinates at this face: those whose direction cosine
        # along the outward normal is positive.
        outgoing_mask = (self.outward_sign * mu_n) > 0.0
        cos_w = weights * (self.outward_sign * mu_n)
        # Cosine-weighted outgoing-current normalisation. ``np.where``
        # zeroes contributions from non-outgoing ordinates.
        cos_w = np.where(outgoing_mask, cos_w, 0.0)

        norm = cos_w.sum()
        if norm <= 0.0:
            # Degenerate quadrature — no outgoing ordinates. Return
            # zero rather than producing a NaN.
            return np.zeros_like(angular_flux_outgoing)

        # Cosine-weighted average of the outgoing flux.
        # Shape (N_ord, ...) → broadcast (..., ) average.
        psi_avg = (
            cos_w.reshape((-1,) + (1,) * (angular_flux_outgoing.ndim - 1))
            * angular_flux_outgoing
        ).sum(axis=0) / norm

        # Broadcast over all ordinates; sweeps consume only entries
        # whose direction is *incoming* at the face, but it is cheap
        # and conventional to fill the whole array uniformly.
        result = np.broadcast_to(
            psi_avg[None, ...] * self.albedo,
            angular_flux_outgoing.shape,
        ).copy()
        return result


@dataclass(frozen=True)
class PeriodicBoundaryOperator:
    r"""Periodic boundary: spatial pushforward to the partner face.

    Tensor decomposition :math:`(G_{\text{wrap}}, 1)` where
    :math:`G_{\text{wrap}}` is the spatial pushforward to the
    partner face of the domain (e.g. left ↔ right): the incoming flux
    at this face equals the outgoing flux at the partner face for
    every ordinate, with no angular permutation:

    .. math::

        \psi_{\text{in}}(\Omega, x_{\text{this}})
        = \psi_{\text{out}}(\Omega, x_{\text{partner}}).

    Realising periodicity at the *sweep* level requires coupling the
    two faces' boundary-flux buffers — which is a sweep-orchestration
    concern not modelled by ``apply_to_incoming`` alone (the
    :class:`BoundaryOperator` Protocol consumes one face's outgoing flux at
    a time). The primitive here therefore returns ``angular_flux_outgoing``
    unchanged: the contract is "the incoming side equals the outgoing
    flux you pass in", and the *spatial* coupling is handled by whoever
    instantiates :class:`PeriodicBoundaryOperator` and orchestrates the two-face
    plumbing. This is why periodic-BC support in :func:`solve_sn` is a
    downstream wave (it requires sweep changes); this primitive ships
    so that downstream code has the algebraic object to depend on.
    """

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Per the docstring above: the *partner-face* outgoing flux is
        # what this BC needs; the caller passes that array in via the
        # ``angular_flux_outgoing`` argument and we return it unchanged.
        # Angular structure is identity; spatial pushforward is the
        # caller's responsibility.
        return angular_flux_outgoing.copy()


@dataclass(frozen=True)
class AlbedoBoundaryOperator:
    r"""Pure albedo boundary: scalar multiple of the outgoing flux.

    Tensor decomposition :math:`(I, \alpha)` where :math:`I` is the
    angular identity and :math:`\alpha \in [0, 1]` is the albedo:

    .. math::

        \psi_{\text{in}}(\Omega) = \alpha \, \psi_{\text{out}}(\Omega).

    No angular redistribution. Useful as a *building block* for
    :class:`MixedBoundaryOperator` (where albedo and specular shares are independent
    parameters), and as a stand-alone primitive when the boundary is
    a pure attenuator with no angular structure.

    Parameters
    ----------
    albedo : float
        Albedo coefficient. ``0`` is vacuum, ``1`` is perfect
        same-direction return.
    """

    albedo: float = 0.0

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        return self.albedo * angular_flux_outgoing


@dataclass(frozen=True)
class MixedBoundaryOperator:
    r"""Linear combination of rank-1 BC primitives.

    Realises a rank-N tensor decomposition

    .. math::

        R = \sum_\alpha c_\alpha \, G_\alpha \otimes A_\alpha

    by storing a list of ``(coefficient, primitive)`` pairs and
    summing ``coefficient * primitive.apply_to_incoming(...)``. The
    coefficients :math:`c_\alpha` are typically convex (sum to 1) for
    a partial-current Marshak boundary
    (Bell & Glasstone 1970 §1.5: ``MixedBoundaryOperator([(0.3, SpecularBoundaryOperator()),
    (0.7, WhiteBoundaryOperator())])`` is "30% specular, 70% diffuse"); the linear-
    combination interface does not enforce this so other use cases
    (asymmetric weights, gain media) can also be expressed.

    Parameters
    ----------
    components : sequence of (coefficient, BoundaryOperator)
        The rank-N decomposition. Each component contributes
        ``coefficient * primitive.apply_to_incoming(...)`` to the
        incoming flux.
    """

    components: tuple[tuple[float, BoundaryOperator], ...] = field(default_factory=tuple)

    def __init__(
        self,
        components: Sequence[tuple[float, BoundaryOperator]],
    ) -> None:
        # Frozen-dataclass-with-Sequence-arg pattern: take a Sequence,
        # store as a tuple. ``object.__setattr__`` to bypass the frozen
        # guard during construction.
        object.__setattr__(self, "components", tuple(components))

    def apply_to_incoming(
        self,
        angular_flux_outgoing: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        result = np.zeros_like(angular_flux_outgoing)
        for coeff, primitive in self.components:
            result = result + coeff * primitive.apply_to_incoming(
                angular_flux_outgoing, quadrature,
            )
        return result


__all__ = [
    "BoundaryOperator",
    "VacuumBoundaryOperator",
    "SpecularBoundaryOperator",
    "WhiteBoundaryOperator",
    "PeriodicBoundaryOperator",
    "AlbedoBoundaryOperator",
    "MixedBoundaryOperator",
]
