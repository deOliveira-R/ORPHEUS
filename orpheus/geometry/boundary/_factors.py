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
distinction is decidable — but it takes **two** tests, not one, and an earlier
draft of this module shipped only the first.

**The necessary test — multiplicativity.**

    :math:`G` is the **composition (Koopman) operator of a measure-preserving
    bijection of the boundary phase space** — :math:`(G\psi)(x) = \psi(g^{-1}x)`
    for some :math:`g` acting on :math:`\partial\Omega \times S^d`.

Such operators are invertible, preserve the trace measure
:math:`|\Omega\cdot\hat n|\,d\Omega\,dA`, form a **group**, and are
**multiplicative**:

.. math::

    G(\psi\,\varphi) \;=\; (G\psi)\,(G\varphi).

A relabeling satisfies that identity. **An averaging operator never does.**
Anything that *fails* multiplicativity is a kernel, and kernels are :math:`R`
— which is what disqualifies the Lambertian average, and is exactly the
argument phase B3 used to move it out of the geometry slot.

**Why that test is not sufficient.** A *specular kernel* is a permutation,
hence multiplicative too. So multiplicativity alone cannot separate

* a polished wall that returns :math:`\alpha` of the flux specularly, from
* a symmetry plane, across which the domain genuinely continues.

Both relabel; only one is geometry. Reading multiplicativity as *the* criterion
would put a surface's re-emission law in :math:`G`, and the two objects have
nothing in common but their matrix.

**The sufficient test — is it a quotient?**

    :math:`G` is the deck transformation **of an actual quotient of the
    physical domain**.

A physical surface is not a quotient: the domain does *not* continue on the
other side of it, so nothing is identified with anything and there is no deck
group to be an element of. A surface's specular pairing is therefore
**constitutive** — it is :math:`R`. This is the quotient table below, promoted
from an observation to a test.

**⇒ EXACTLY ONE of** :math:`G`, :math:`R` **is non-trivial**, which is this
module's own sentence *"*:math:`R = I` *exactly when the BC is a pure symmetry
statement adding no physics"* read as a law rather than as a remark. Its
contrapositive is the useful direction: a law that asserts any physics at all
has :math:`G = \mathrm{id}`.

Physically: *a change of direction caused by the* **geometry** *is* :math:`G`;
*a change of direction caused by the* **constitutive assumption of the BC** *is*
:math:`R`. Absorption, diffusivity and accommodation are :math:`R`. Mirrors,
translations and rotations are :math:`G` — but only when they are mirrors
*of the domain*, not of a wall standing in it.

.. list-table:: The law table this yields
   :header-rows: 1
   :widths: 34 22 22 22

   * - law
     - :math:`G`
     - :math:`R`
     - what it asserts
   * - ``ReflectiveBoundary(axis)``
     - :class:`SpecularMirror`
     - :math:`I`
     - a symmetry plane — a quotient, **zero physics**
   * - ``PeriodicBoundary(axis)``
     - :class:`SpatialWrap`
     - :math:`I`
     - a torus — a quotient
   * - ``AlbedoBoundary(α, SpecularReturn)``
     - :class:`IdentityMap`
     - :class:`SpecularReemission`
     - a **surface** returning :math:`\alpha` specularly
   * - ``AlbedoBoundary(α, IsotropicReturn)``
     - :class:`IdentityMap`
     - :class:`LambertianReemission`
     - a surface returning :math:`\alpha` diffusely
   * - ``VacuumInflow``
     - :class:`IdentityMap`
     - :math:`0`
     - a surface returning nothing

:class:`SpecularReturn` and :class:`SpecularMirror` realize to the *same*
permutation and are nonetheless different types. That is the **point**, not a
smell: two types make "put a surface's response in the geometry slot"
unspellable, which is the exact error this section corrects.

.. note::

   One row of the shipped code violates the law, deliberately and visibly:
   :class:`~orpheus.geometry.boundary.ReflectiveBoundary` still accepts an
   ``albedo`` parameter, so ``ReflectiveBoundary(axis, 0.7)`` has BOTH factors
   non-trivial. A symmetry plane cannot absorb — that object is
   ``AlbedoBoundary(0.7, SpecularReturn(axis))`` wearing the geometry costume.
   It is unreachable from a ``BC(...)`` tag (the tag parser hard-codes
   :math:`\alpha = 1`), so nothing production-facing rides on it; retiring the
   parameter is campaign phase **B5**.

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

The **response tier** now has three members and earns the split: a scalar
amplitude, a rank-one angular kernel and an angular permutation are pairwise
non-isomorphic, and each has a genuine non-identity morphism applied to it (a
multiplication, a contract-then-broadcast, a gather). §16A.2's
``ZeroBoundaryResponse`` is still NOT minted — it would be an isomorphic
singleton of :class:`ScalarResponse`, so :attr:`~ScalarResponse.is_zero` stays
a property.

The **geometry tier** likewise: a mirror and a spatial wrap are genuinely
non-isomorphic deck transformations, and :class:`IdentityMap` is the group's
identity element rather than a fourth realization.

The **closure tier** (:class:`ReemissionClosure`) is the newest and needs its
own defence, because it is one indirection above the kernel it produces.

A surface's re-emission law has **two independent degrees of freedom**: *how
much* comes back (:math:`\alpha`) and *in what angular shape*. They vary
independently — every shape admits every amplitude — so the elegant spelling
takes them as two parameters and multiplies, rather than offering a menu of
their product. :class:`SpecularReturn` and :class:`IsotropicReturn` are the
amplitude-**free** shapes; :meth:`~ReemissionClosure.kernel` is the
:math:`\alpha`-instantiation morphism that combines them. Two non-isomorphic
members, a non-identity morphism actually applied: the mint is earned.

The alternative — putting a fully-formed :class:`BoundaryResponseKernel` on the
law — was rejected because :attr:`~BoundaryResponseKernel.amplitude` is a
Protocol member the diffusion realizer already reads, so the kernel *must*
carry :math:`\alpha`; a law carrying both would hold **two sources of one
number**. Keeping the closure amplitude-free leaves :math:`\alpha` exactly
where the tag parser already puts it — on the law.

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
    "IsotropicReturn",
    "LambertianReemission",
    "ReemissionClosure",
    "ScalarResponse",
    "SpatialWrap",
    "SpecularMirror",
    "SpecularReemission",
    "SpecularReturn",
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


@dataclass(frozen=True, slots=True)
class SpecularReemission:
    r""":math:`R = \alpha\,P_{\text{mirror}}` — **specular** re-emission.

    A surface that returns :math:`\alpha` of the arriving flux into the
    mirror-partner direction. Realized as the ordinate permutation
    ``quadrature.reflection_index(axis)``, scaled by :math:`\alpha` — the same
    matrix :class:`SpecularMirror` realizes to.

    Why this is a **response** and not a geometry
    ---------------------------------------------

    It is the one kernel that passes the multiplicativity test — a permutation
    is multiplicative — and it is still :math:`R`, which is why that test alone
    was never sufficient (see this module's criterion section). The
    discriminator is the quotient one:

        A **symmetry plane** is a quotient of the domain. Nothing happens at
        it; the flux simply continues, relabeled, and :math:`R = I` records
        that no physics was asserted.

        A **polished wall** is not a quotient. The domain stops there. That the
        wall happens to return flux in the mirror direction is a *constitutive*
        statement about the wall — as constitutive as a Lambertian surface's
        diffuse return, and made in the same slot.

    So ``AlbedoBoundary(α, SpecularReturn(axis))`` and
    ``ReflectiveBoundary(axis, α)`` realize to the *same matrix* and assert
    *different things*: the first says "a wall here reflects :math:`\alpha`",
    the second says "the domain is symmetric about this plane, and (at
    :math:`\alpha < 1`, incoherently) also absorbs". Only the first survives
    the criterion at :math:`\alpha < 1`; the second's ``albedo`` parameter is
    retired in phase **B5**.

    The distinction is not academic. :attr:`permutes_ordinates` on the geometry
    tier drives the SN sweep schedule — a *quotient* couples ordinates across
    the sweep and can close a cycle; a *wall's* re-emission is an ordinary
    boundary source term. Reading one as the other is how a sweep gets
    needlessly lagged, or worse, wrongly not lagged.

    Parameters
    ----------
    alpha : float
        Fraction returned. Not clamped here — same rationale as
        :class:`ScalarResponse`: the sub-Markov bound is a per-law invariant.
    axis : str
        Mirror axis, ``"x"`` / ``"y"`` / ``"z"``. The permutation itself needs
        a quadrature, which is why this is a spec and not an operator.
    """

    alpha: float = 1.0
    axis: str = "x"

    @property
    def amplitude(self) -> float:
        return float(self.alpha)

    @property
    def is_zero(self) -> bool:
        # Same exact-compare rationale as ScalarResponse.
        return self.alpha == 0.0

    @property
    def is_adjointable(self) -> bool:
        # TRUE, unlike the Lambertian's, and for a structural reason rather
        # than an implementation one: a scaled permutation's transpose is the
        # inverse permutation scaled by the same α, and an axis reflection is
        # an involution, so the transpose is the forward action itself. No
        # inner-product ambiguity to decline — it is symmetric under both the
        # Euclidean and the cosine-weighted metric, because the permutation
        # preserves |Ω·n| pointwise.
        return True


# ═══════════════════════════════════════════════════════════════════════
# The re-emission closure tier — amplitude-free angular SHAPES
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class ReemissionClosure(Protocol):
    r"""The angular **shape** in which a surface returns flux — no amplitude.

    A re-emission law is a product of two independent choices: *how much*
    returns (:math:`\alpha`) and *in what direction*. This tier carries the
    second alone, and :meth:`kernel` combines it with the first.

    Why a surface needs one at all
    ------------------------------

    :math:`R = \alpha\,I` — a bare amplitude with no shape — is **complete on a
    scalar trace** and **incomplete on an angular one**. On a scalar trace
    there is one degree of freedom, :math:`J^- = \alpha J^+`, and the angular
    distribution has nothing to fix. On an angular trace :math:`\alpha\,I` is a
    map :math:`\Gamma_+ \to \Gamma_+` — an endomorphism of the *outgoing*
    hemisphere — while the law it must produce is
    :math:`\Gamma_+ \to \Gamma_-`. The identity supplies no crossing, and
    :class:`IdentityMap` (the geometry of every surface law) supplies none
    either. Something must say which outgoing direction feeds which incoming
    one; that something is this closure.

    Composing it anyway, without a closure, is not an approximation — it is a
    pairing by **array position**: incoming ordinate :math:`j` would receive
    :math:`\alpha` times the :math:`j`-th outgoing ordinate, an artefact of
    index order carrying no geometry. That is why the SN realizer refuses the
    closure-free spelling instead of choosing a default for you.

    This is the **angular-resolution** axis of the method-realizability
    taxonomy (``bc-method-realizability`` in this package's ``__init__``)
    biting for the first time: the same law is complete for one method and
    under-determined for another, because the methods resolve different
    coordinates.
    """

    def kernel(self, alpha: float) -> BoundaryResponseKernel:
        r"""Instantiate this shape at amplitude :math:`\alpha`.

        The law supplies :math:`\alpha` from its own field, so the number has
        exactly one home.
        """
        ...


@dataclass(frozen=True, slots=True)
class SpecularReturn:
    r"""Return into the **mirror-partner** direction — a polished surface.

    Deliberately a different type from :class:`SpecularMirror`, which realizes
    to the same permutation. The two are structurally identical and
    semantically disjoint, and keeping them apart is what makes "a wall's
    response in the geometry slot" unspellable — the exact conflation phase B3
    found in the white law and the user's 2026-08-01 ruling generalized.
    """

    axis: str = "x"

    def kernel(self, alpha: float) -> "SpecularReemission":
        return SpecularReemission(alpha=alpha, axis=self.axis)


@dataclass(frozen=True, slots=True)
class IsotropicReturn:
    r"""Return **diffusely** (Lambertian) — a matte surface.

    Produces the same :class:`LambertianReemission`
    :class:`~orpheus.geometry.boundary.WhiteBoundary` carries, which is what
    makes ``AlbedoBoundary(α, IsotropicReturn(axis, sign))`` and
    ``WhiteBoundary(axis, sign, α)`` the same law under two names — a genuine
    duplication in the *declaration* vocabulary, and the reason both realize
    through one shared body rather than two transcriptions.

    Needs ``outward_sign`` as well as ``axis`` because a cosine-weighted
    average is over a *hemisphere*, and which hemisphere is outgoing depends on
    which end of the axis the face sits at. A mirror needs no sign: it swaps
    the two hemispheres either way.
    """

    axis: str = "x"
    outward_sign: int = +1

    def kernel(self, alpha: float) -> "LambertianReemission":
        return LambertianReemission(
            alpha=alpha, axis=self.axis, outward_sign=self.outward_sign,
        )
