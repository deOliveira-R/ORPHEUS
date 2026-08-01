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

What separates :math:`G` from :math:`R` — the decidable criterion
=================================================================

The two factors are not "the first thing that happens" and "the second thing
that happens". They are **different kinds of mathematical object**, and the
distinction is decidable:

    :math:`G` is the **composition (Koopman) operator of a measure-preserving
    bijection of the boundary phase space** — :math:`(G\psi)(x) = \psi(g^{-1}x)`
    for some :math:`g` acting on :math:`\partial\Omega \times S^d`.

Such operators are invertible, preserve the trace measure
:math:`|\Omega\cdot\hat n|\,d\Omega\,dA`, form a **group**, and — the test that
decides membership — are **multiplicative**:

.. math::

    G(\psi\,\varphi) \;=\; (G\psi)\,(G\varphi).

A relabeling satisfies that identity. **An averaging operator never does.**
Anything that fails multiplicativity is a kernel, and kernels are :math:`R`.

Physically: *a change of direction caused by the* **geometry** *is* :math:`G`;
*a change of direction caused by the* **constitutive assumption of the BC** *is*
:math:`R`. Absorption, diffusivity and accommodation are :math:`R`. Mirrors,
translations and rotations are :math:`G`.

The crossing :math:`\Gamma_+ \to \Gamma_-` is itself geometric
--------------------------------------------------------------

The specular mirror :math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n`
is the unique ambient isometry fixing the face. It **exchanges the hemispheres**
and preserves :math:`|\Omega\cdot\hat n|`. So the crossing is not something the
*physics* does — it is something the *geometry* provides, and the honest typing
is

.. math::

    G : \Gamma_+ \to \Gamma_-, \qquad R : \Gamma_- \to \Gamma_-.

An earlier draft of this module typed ``G`` as an **endomorphism of the outflow
trace** (:math:`\Gamma_+ \to \Gamma_+`), which forces the mirror *into*
:math:`R` for a specular law — the same conflation running the other way.

Why the misassignment survived — a theorem
------------------------------------------

If :math:`R` is rank-one, :math:`R = u \otimes v`, then
:math:`R \circ G = u \otimes (G^{\mathsf T} v)`. The Lambertian's
:math:`v = |\Omega\cdot\hat n|` is **invariant** under both the mirror and the
periodic translation, so :math:`R \circ G = R`:

    :math:`G` **is unobservable exactly when** :math:`R` **is rank-one.**

White is precisely that case. Its :math:`G` slot was therefore free of any
observable consequence, and the physics drifted into it: the cosine-weighted
Lambertian average shipped as a ``BoundaryGeometryMap`` named
``HemisphericalAverage``, whose own docstring conceded it was *"an all-to-all
coupling, not a relabeling"* while implementing a Protocol documented as *"it
only relabels"*. Phase **B3** moved it to :class:`LambertianReemission`, where
it belongs.

This theorem is also *why the correction is safe*: on the one law whose factors
moved, the composite :math:`R\,G` is unchanged by construction.

The quotient picture — which BC is the orbifold
------------------------------------------------

:math:`G` is the **deck transformation of the quotient by which the physical
domain is represented**:

.. list-table::
   :header-rows: 1
   :widths: 18 30 20 32

   * - BC
     - quotient
     - fixed points
     - what it is
   * - periodic
     - :math:`\mathbb{R}^d/\Lambda` by a translation
     - none — **free**
     - a torus; a genuine **covering space**, a manifold
   * - reflective
     - by a reflection
     - the mirror plane
     - an **orbifold** (Thurston *reflector* boundary)
   * - rotational (⅛-core)
     - by a finite rotation
     - the axis
     - an **orbifold** (cone points)

So the orbifold label attaches to **reflective**, not to periodic — periodic is
the free / covering-space case. And :math:`R = I` **exactly when** the BC is a
pure symmetry statement adding no physics. Vacuum, white and albedo are not
symmetry statements at all, which is why their :math:`G` is the identity deck
element and all their content sits in :math:`R`.

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

   **No ``as_operator`` yet — deliberately.** These types are pure data, and
   nothing but the predicates below is read. The realization method arrives in
   campaign phase **B4**, together with its first consumer (the realizer
   switching from ``isinstance``-dispatch to reading these factors), gated
   bit-identical against today's inline construction. Minting a method with no
   caller and no test is precisely the dead-capability pattern this campaign
   exists to remove — see the review's §4, which catalogues five instances of
   it in this subsystem alone.

Type-minting discipline
=======================

``coding-standards.md`` mints a type **iff** there are ≥2 non-isomorphic
realizations AND a non-identity morphism is applied.

The **response tier** now has two members and earns the split: a scalar
amplitude and a rank-one angular kernel are not isomorphic, and B3 applies a
genuine non-identity morphism to each (a multiplication versus a
contract-then-broadcast). §16A.2's ``ZeroBoundaryResponse`` is still NOT minted
— it would be an isomorphic singleton of :class:`ScalarResponse`, so
:attr:`~ScalarResponse.is_zero` stays a property.

The **geometry tier** likewise: a mirror and a spatial wrap are genuinely
non-isomorphic deck transformations, and :class:`IdentityMap` is the group's
identity element rather than a fourth realization.

``NullMap`` was retired in B3 for the same reason the Lambertian moved. It
declared :math:`G = 0`, but **the zero map is not a bijection**, so it was never
a geometry map at all — vacuum's and prescribed-inflow's zero-ness is a property
of their *response*, which already spells it as ``ScalarResponse(0.0)``. Two
spellings of one fact, in the wrong two tiers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable


__all__ = [
    "BoundaryGeometryMap",
    "BoundaryResponseKernel",
    "IdentityMap",
    "LambertianReemission",
    "ScalarResponse",
    "SpatialWrap",
    "SpecularMirror",
]


# ═══════════════════════════════════════════════════════════════════════
# The two Protocols
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class BoundaryGeometryMap(Protocol):
    r""":math:`G : \Gamma_+ \to \Gamma_-` — the **deck transformation**.

    The composition operator of a measure-preserving bijection of the boundary
    phase space: it changes nothing about the physical interaction at the
    boundary, it only relabels which outgoing flux meets which incoming slot.
    The physics lives entirely in :class:`BoundaryResponseKernel`.

    Membership is decidable by **multiplicativity** — see
    :ref:`bc-factor-roles` on the theory page. The crossing
    from :math:`\Gamma_+` to :math:`\Gamma_-` is part of this factor, not of
    the response (:ref:`bc-factor-roles`).

    The two predicates below are the structural questions production used to
    ask with string comparisons — ``bc[face] == "reflective"`` for the first,
    ``kind in _RULED_CORNER_KINDS`` for the second. Phase **B2** repointed
    those sites here and retired both tag sets; the law-by-law equivalence is
    pinned in ``tests/geometry/test_boundary_factor_consumers.py``.
    """

    @property
    def permutes_ordinates(self) -> bool:
        """Whether realizing this map permutes the ANGULAR index.

        ``True`` only for a specular mirror today. A spatial wrap permutes the
        *spatial* index and leaves angle alone, so it answers ``False`` — the
        distinction matters to the sweep schedule, which cares about angular
        coupling within a face. The identity deck element answers ``False``:
        it relabels nothing.
        """
        ...

    @property
    def is_adjointable(self) -> bool:
        r"""Whether the realized map exposes an honest transpose.

        For a genuine deck transformation this is a **theorem, not a choice**:
        the composition operator of a bijection :math:`g` is invertible with
        :math:`G^{-1} = G_{g^{-1}}`, and measure-preservation makes that
        inverse the transpose. It is declared rather than derived only because
        :class:`SpatialWrap` answers ``False`` while its realized
        pushforward is still unbuilt (#183) — an implementation gap, not a
        property of the map.
        """
        ...


@runtime_checkable
class BoundaryResponseKernel(Protocol):
    r""":math:`R : \Gamma_- \to \Gamma_-` — the **constitutive response**.

    How much of the outgoing flux returns, and with what angular distribution.
    This is where every piece of *physics* in a boundary law lives: absorption,
    accommodation, diffusivity. It is a positive kernel operator, and — the
    discriminator — it is **not multiplicative**.

    Two realizations, genuinely non-isomorphic: :class:`ScalarResponse` (a bare
    amplitude, the whole story on a scalar trace where the angular distribution
    has no degrees of freedom) and :class:`LambertianReemission` (a rank-one
    angular kernel). The theory page anticipated exactly this tier as
    *"a full angular kernel in general weak-form BCs (deferred)"*; B3 mints it.

    :class:`~orpheus.geometry.boundary.ZeroFluxBoundary` is expressible here but
    sits outside the sub-Markov range, at :math:`\mathcal{A} = -1` in the
    partial-current basis — which is exactly what the diffusion realizer builds
    (``ScaledOperator(-1.0, IdentityOperator)``). Note the distinction: its
    *realization* is affine, while its *posing* — :math:`\phi_\Gamma = 0`, i.e.
    :math:`A_-\gamma_- + A_+\gamma_+ = 0` — is a **relation**, a tier above the
    affine trace law (Grand Report v3; issue **#177**). Populating the factor is
    honest about the former without claiming the latter.
    """

    @property
    def amplitude(self) -> float:
        r"""The scalar magnitude :math:`\alpha` with which flux returns.

        For :class:`ScalarResponse` this **is** the whole kernel. For a kernel
        with angular structure it is only the magnitude — but that is precisely
        the dimension-reduced view a **scalar trace** sees, which is why the
        diffusion realizer reads this one number and needs no geometric factor
        at all. Both realizers already reached it as a bare float before B1
        gave it a name.
        """
        ...

    @property
    def is_zero(self) -> bool:
        """Whether this response returns nothing — the vacuum/prescribed case."""
        ...

    @property
    def is_adjointable(self) -> bool:
        """Whether the realized kernel exposes an honest transpose.

        ``False`` for the Lambertian today: its realized form is self-adjoint
        under the cosine-weighted inner product but NOT under the Euclidean
        one, and the codebase declines to advertise the ambiguous transpose.
        Phase **B5** types that kernel as the rank-one it is, at which point
        the adjoint becomes structurally available and this flips.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# Concrete deck transformations
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class IdentityMap:
    r""":math:`G = \mathrm{id}` — the **identity element** of the deck group.

    The law imposes no geometric identification: the domain is represented by
    itself, not by a quotient. Every law carrying this map is one whose content
    is entirely constitutive — vacuum, prescribed inflow, white, albedo,
    zero-flux — and for all of them the crossing to :math:`\Gamma_-` is the
    face's own canonical one, which the realizer supplies.

    That this is *sound* rather than merely conventional is the rank-one
    theorem (:ref:`bc-factor-roles`): where the response destroys
    directional information, :math:`R \circ G = R` for **any** measure-
    preserving :math:`G`, so declaring the identity is not a lossy choice — it
    is the honest statement that the law fixes no geometry.
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

    The reflection :math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n` —
    the unique ambient isometry fixing the face. It exchanges the hemispheres
    (which is why it, and not the response, carries the
    :math:`\Gamma_+ \to \Gamma_-` crossing) and preserves
    :math:`|\Omega\cdot\hat n|`, so it is measure-preserving in the trace
    measure, as ERR-042 requires.

    Quotient reading: this is the deck transformation of a quotient **with
    fixed points** — the mirror plane itself — so a reflecting face makes the
    computational domain an **orbifold** with what Thurston calls a reflector
    boundary. (Contrast :class:`SpatialWrap`, whose action is free.)

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
class SpatialWrap:
    r""":math:`G` = wrap-around along ``axis`` (periodic).

    Pushes the outflow of one face onto the inflow of the opposite face at the
    SAME angle — which is why periodic closes a sweep cycle from a *single* law
    while a lone reflecting face does not.

    It is a genuine deck transformation: the translation :math:`x \mapsto x + L`
    carries face :math:`f'` onto face :math:`f`, and because the outward
    normals are opposite (:math:`\hat n_f = -\hat n_{f'}`) a direction that is
    *outgoing* at :math:`f'` is *incoming* at :math:`f` — so the crossing comes
    for free, with :math:`|\Omega\cdot\hat n|` preserved. Quotient reading: the
    action is **free** (no fixed points), so the quotient is a torus — a
    covering space and a manifold, NOT an orbifold.

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
        # The translation IS adjointable as a map (see the Protocol's note —
        # for a deck transformation this is a theorem). ``False`` reports an
        # IMPLEMENTATION gap: the realized operator is currently an angular
        # identity with the spatial pushforward unbuilt (#183). B3.4 builds it
        # and this flips, WITH its gate.
        return False


# ═══════════════════════════════════════════════════════════════════════
# Concrete response kernels
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class ScalarResponse:
    r""":math:`R = \alpha\,I` — a scalar amplitude, no angular redistribution.

    The flux returns attenuated but otherwise untouched, so whatever angular
    structure the deck transformation produced is what arrives. Paired with
    :class:`SpecularMirror` this is the specular albedo; paired with
    :class:`IdentityMap` on a **scalar** trace it is the partial-current albedo,
    where it is the complete story because the angular distribution has no
    degrees of freedom to fix.

    Construction does NOT clamp to :math:`[0, 1]`: the sub-Markov bound is a
    per-law invariant (``assert_albedo_in_unit_interval``), and the zero-flux
    idealization deliberately sits at :math:`-1`. Enforcing it here would
    duplicate the law-level invariant and make a legitimate value
    unrepresentable — issue #265 tracks the value-object treatment.
    """

    alpha: float = 1.0

    @property
    def amplitude(self) -> float:
        return float(self.alpha)

    @property
    def is_zero(self) -> bool:
        # Exact compare, deliberately: this answers "is this law's response
        # structurally absent" (vacuum, prescribed inflow), not "is it small".
        # A near-zero albedo is a weak reflector, not a vacuum.
        return self.alpha == 0.0

    @property
    def is_adjointable(self) -> bool:
        # A multiple of the identity is symmetric in every inner product.
        return True


@dataclass(frozen=True, slots=True)
class LambertianReemission:
    r""":math:`R = \alpha\,\frac{1}{\pi}\,|\Omega\cdot\hat n|
    \;\langle |\cdot\,\hat n|,\; \cdot\rangle` — diffuse re-emission.

    Contract the returning hemisphere against :math:`|\Omega\cdot\hat n|\,w` and
    rebroadcast the result isotropically: a **rank-one** kernel, currently
    realized as a hand-written contract-then-broadcast rather than typed as one
    (phase **B5**).

    **This is a response, not a geometry** — the distinction B3 corrected. An
    average is not multiplicative and is not a bijection, so it fails the
    membership test for :math:`G` (:ref:`bc-factor-roles`); and
    physically it is a *constitutive* statement about how the surface re-emits,
    not a symmetry of the domain. It shipped as a ``BoundaryGeometryMap`` named
    ``HemisphericalAverage`` because a rank-one response makes :math:`G`
    unobservable (:ref:`bc-factor-roles`), so the empty slot had no
    observable consequence and the physics settled into it.

    Watch the two thresholds. The realized operator's outward test must be the
    trace space's ``TANGENTIAL_EPS`` classification, NOT a strict ``> 0.0``:
    they disagree wherever a quadrature carries tangential ordinates, which is
    every production quadrature except ``gauss_legendre(4)``. Measured on a
    cylinder under ``product(2, 4)``, a strict compare reads two rows the trace
    calls tangential and the composite diverges by 6.1e-05.
    """

    alpha: float = 1.0
    axis: str = "x"
    outward_sign: int = +1

    @property
    def amplitude(self) -> float:
        return float(self.alpha)

    @property
    def is_zero(self) -> bool:
        # Same exact-compare rationale as ScalarResponse: a zero-amplitude
        # Lambertian returns nothing at all, which IS the vacuum answer.
        return self.alpha == 0.0

    @property
    def is_adjointable(self) -> bool:
        # FALSE TODAY, and honestly so: the realized operator is self-adjoint
        # under the cosine-weighted inner product but not the Euclidean one,
        # and advertising the unweighted transpose would invite two different
        # ``.T`` semantics. B5 types it as ``u ⊗ v`` (transpose ``v ⊗ u``),
        # which makes the metric explicit rather than avoided — this flips
        # there, WITH its gate.
        return False
