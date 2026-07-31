r"""The affine boundary form's two factors, as typed **specifications**.

A boundary law is the affine map on the trace

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q,

and Grand Report v3 §16A.2 specifies its three members as typed fields rather
than as free-floating parameters:

.. code-block:: python

    @dataclass(frozen=True, slots=True)
    class BoundaryTraceLaw:
        geometry_map: BoundaryGeometryMap
        response: BoundaryResponseKernel
        source: BoundarySource

This module mints the first two. They were declared on
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` as ``-> Any`` properties
defaulting to ``None`` and **nothing ever populated them** — while five
production sites answered the questions they exist to answer by comparing
``law.kind`` strings.

Specification, not operator
===========================

The factors carry **what the geometry/response IS**, never a realized matrix.
That split is the design's, not an invention: §16A.2's realization step takes
the discretization as an argument —

.. code-block:: python

    def as_operator(self, trace_space): ...

— so the *law* owns the spec and the *trace space* produces the matrix. It is
also what makes the factors populatable at all: a specular mirror's realized
:math:`G` is a permutation **of ordinates**, which needs a quadrature the
method-agnostic law does not have; ``SpecularMirror(axis="x")`` needs nothing.

.. note::

   **No ``as_operator`` yet — deliberately.** These types are pure data in B1,
   which is a *pure addition*: nothing reads them, so nothing can regress.
   The realization method arrives in campaign phase **B4**, together with its
   first consumer (the realizer switching from ``isinstance``-dispatch to
   reading these factors), gated bit-identical against today's inline
   construction. Minting a method with no caller and no test is precisely the
   dead-capability pattern this campaign exists to remove — see the review's
   §4, which catalogues five instances of it in this subsystem alone.

Type-minting discipline — one deviation from §16A.2, stated
===========================================================

§16A.2 names a ``ZeroBoundaryResponse`` type and §16A.5 dispatches on
``isinstance(law.response, ZeroBoundaryResponse)``. This module ships **one**
:class:`ScalarResponse` instead, with :attr:`~ScalarResponse.is_zero` as a
property.

Reason: ``coding-standards.md`` mints a type **iff** there are ≥2
non-isomorphic realizations AND a non-identity morphism is applied. Measured,
**every** response in this law family is a scalar amplitude — both realizers
already reach it as a bare float (``float(law.albedo) * base`` in the SN arm,
``return float(law.albedo)`` in the diffusion arm). ``ZeroResponse`` and
``UnitResponse`` would be isomorphic singletons of ``ScalarResponse``, i.e.
ceremony around a value. The geometry maps below DO earn separate types —
a permutation, a hemispheric average and a spatial wrap are genuinely
non-isomorphic realizations.

If a future weak-form BC needs a non-scalar response (a full angular kernel,
which the theory page already anticipates as *"deferred"*), that is a second
non-isomorphic realization and the type split becomes earned at that point.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable


__all__ = [
    "BoundaryGeometryMap",
    "BoundaryResponseKernel",
    "HemisphericalAverage",
    "IdentityMap",
    "NullMap",
    "ScalarResponse",
    "SpatialWrap",
    "SpecularMirror",
]


# ═══════════════════════════════════════════════════════════════════════
# The two Protocols
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class BoundaryGeometryMap(Protocol):
    r""":math:`G : \Gamma_+ \to \Gamma_+` — the geometric relabeling.

    An **endomorphism of the outflow trace** (theory page
    :doc:`/theory/foundations/boundary_conditions`): it changes nothing about
    the physical interaction at the boundary, it only relabels which outgoing
    flux meets which incoming slot. The physics lives entirely in
    :class:`BoundaryResponseKernel`.

    The two predicates below are the structural questions production currently
    asks with string comparisons — ``bc[face] == "reflective"`` for the first,
    ``kind in _RULED_CORNER_KINDS`` for the second. Phase **B2** repoints those
    sites here.
    """

    @property
    def permutes_ordinates(self) -> bool:
        """Whether realizing this map permutes the ANGULAR index.

        ``True`` only for a specular mirror today. A spatial wrap permutes the
        *spatial* index and leaves angle alone, so it answers ``False`` — the
        distinction matters to the sweep schedule, which cares about angular
        coupling within a face.
        """
        ...

    @property
    def is_adjointable(self) -> bool:
        """Whether the realized map exposes an honest transpose.

        ``False`` for the hemispheric average today: its realized form is
        self-adjoint under the cosine-weighted inner product but NOT under the
        Euclidean one, and the codebase declines to advertise the ambiguous
        transpose. Phase **B5** types that map as the rank-one it is, at which
        point the adjoint becomes structurally available and this flips.
        """
        ...


@runtime_checkable
class BoundaryResponseKernel(Protocol):
    r""":math:`R : \Gamma_+ \to \Gamma_-` — the crossing.

    The physical amplitude with which outgoing flux returns as incoming. A
    scalar in :math:`[0, 1]` for every sub-Markov BC in this family (albedo,
    white, partial current).

    :class:`~orpheus.geometry.boundary.ZeroFluxBoundary` is expressible here but
    sits outside that range, at :math:`\mathcal{A} = -1` in the partial-current
    basis — which is exactly what the diffusion realizer builds
    (``ScaledOperator(-1.0, IdentityOperator)``). Note the distinction: its
    *realization* is affine, while its *posing* — :math:`\phi_\Gamma = 0`, i.e.
    :math:`A_-\gamma_- + A_+\gamma_+ = 0` — is a **relation**, a tier above the
    affine trace law (Grand Report v3; issue **#177**). Populating the factor is
    honest about the former without claiming the latter.
    """

    @property
    def scalar(self) -> float:
        """The amplitude, as the bare float both realizers already multiply by."""
        ...

    @property
    def is_zero(self) -> bool:
        """Whether this response returns nothing — the vacuum/prescribed case."""
        ...


# ═══════════════════════════════════════════════════════════════════════
# Concrete geometry maps
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class IdentityMap:
    r""":math:`G = I` — outgoing slot :math:`n` feeds incoming slot :math:`n`.

    The albedo family's geometry: no relabeling at all, the response scalar
    carries the entire law.
    """

    @property
    def permutes_ordinates(self) -> bool:
        return False

    @property
    def is_adjointable(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class NullMap:
    r""":math:`G = 0` — no outgoing flux reaches the inflow trace.

    The rank-0 geometry, shared by :class:`~orpheus.geometry.boundary.VacuumInflow`
    (whose response is also zero) and
    :class:`~orpheus.geometry.boundary.PrescribedInflow` (whose inflow comes
    entirely from :math:`q`). "Rank-0" is the existing in-code vocabulary —
    ``IncomingSourceOperator`` already describes itself as *"the rank-0 case
    where R = G = 0"*.
    """

    @property
    def permutes_ordinates(self) -> bool:
        return False

    @property
    def is_adjointable(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class SpecularMirror:
    r""":math:`G = G_{\text{refl}}` — mirror reflection about ``axis``.

    Realizes (B4) to the ordinate permutation
    ``quadrature.reflection_index(axis)``. The law carries only the axis; the
    permutation itself needs the quadrature, which is exactly why this is a
    spec and not an operator.
    """

    axis: str = "x"

    @property
    def permutes_ordinates(self) -> bool:
        return True

    @property
    def is_adjointable(self) -> bool:
        # A permutation's transpose is its inverse permutation — realized today
        # as ``argsort(perm)``, a genuine transpose rather than a
        # re-application. Measured exact: ‖T − Fᵀ‖∞ = 0.
        return True


@dataclass(frozen=True, slots=True)
class HemisphericalAverage:
    r""":math:`G = G_{\text{diff}}` — cosine-weighted Lambertian average.

    Contracts the outgoing hemisphere against :math:`|\Omega\cdot\hat n|\,w` and
    rebroadcasts the result isotropically — i.e. a **rank-one** map, currently
    realized as a hand-written contract-then-broadcast rather than typed as
    one (phase **B5**).
    """

    axis: str = "x"
    outward_sign: int = +1

    @property
    def permutes_ordinates(self) -> bool:
        # Every outgoing ordinate feeds every incoming one (rank-1 in angle) —
        # an all-to-all coupling, not a relabeling. This is why the production
        # Gauss-Seidel schedule excludes white from its octant split: a single
        # white face already couples all ordinates, so an octant-ordered G-S
        # degenerates to Jacobi.
        return False

    @property
    def is_adjointable(self) -> bool:
        # FALSE TODAY, and honestly so: the realized operator is self-adjoint
        # under the cosine-weighted inner product but not the Euclidean one,
        # and advertising the unweighted transpose would invite two different
        # ``.T`` semantics. B5 types it as ``u ⊗ v`` (transpose ``v ⊗ u``),
        # which makes the metric explicit rather than avoided — this flips
        # there, WITH its gate.
        return False


@dataclass(frozen=True, slots=True)
class SpatialWrap:
    r""":math:`G` = wrap-around along ``axis`` (periodic).

    Pushes the outflow of one face onto the inflow of the opposite face at the
    SAME angle — which is why periodic closes a sweep cycle from a *single* law
    while a lone reflecting face does not.

    **``axis``, not a partner face.** The first draft of this type carried
    ``partner_face``, which is wrong by the rule this campaign's B0 phase
    established: *a law carries only what is intrinsic to it, never what depends
    on the configuration or the discretization.* Which face is the partner
    depends on **where the law is installed** — configuration — whereas "wrap
    along x" is intrinsic, and it is the same shape as
    :class:`SpecularMirror`'s parameter. The realizer derives the partner from
    the installation face plus this axis.

    :class:`~orpheus.geometry.boundary.PeriodicBoundary` has never carried
    ANY field (issue #183), so the map it names was not expressible; ``axis`` is
    the parameter it was missing.

    Non-opposite gluing — a hex partner, a rotational quotient — is genuinely a
    different object and needs an explicit partner map. That is issue **#178**
    (``SymmetryBoundary``, "octant/quotient gluing distinct from physical
    mirror"), deliberately NOT this type.
    """

    axis: str = "x"

    @property
    def permutes_ordinates(self) -> bool:
        # Spatial, not angular: ordinate n at face A feeds ordinate n at face B.
        return False

    @property
    def is_adjointable(self) -> bool:
        # The face-pair swap is its own inverse, but the realized operator is
        # currently an angular identity with the spatial pushforward unbuilt
        # (#183), so there is no honest transpose to advertise yet.
        return False


# ═══════════════════════════════════════════════════════════════════════
# The response
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class ScalarResponse:
    r""":math:`R = \alpha` — a scalar crossing amplitude.

    The response of every law in this family. Both realizers already reach this
    number as a bare float; this type gives it a name and a home rather than a
    new computation.

    Construction does NOT clamp to :math:`[0, 1]`: the sub-Markov bound is a
    per-law invariant (``assert_albedo_in_unit_interval``), and the zero-flux
    idealization deliberately sits at :math:`-1`. Enforcing it here would
    duplicate the law-level invariant and make a legitimate value
    unrepresentable — issue #265 tracks the value-object treatment.
    """

    alpha: float = 1.0

    @property
    def scalar(self) -> float:
        return float(self.alpha)

    @property
    def is_zero(self) -> bool:
        # Exact compare, deliberately: this answers "is this law's response
        # structurally absent" (vacuum, prescribed inflow), not "is it small".
        # A near-zero albedo is a weak reflector, not a vacuum.
        return self.alpha == 0.0
