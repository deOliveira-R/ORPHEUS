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
     - :meth:`SelfPairedDeck.mirror`
     - :math:`I`
     - a symmetry plane — a quotient, **zero physics**
   * - ``PeriodicBoundary(axis)``
     - :class:`SpatialWrap`
     - :math:`I`
     - a torus — a quotient
   * - ``AlbedoBoundary(α, SpecularReturn)``
     - :meth:`SelfPairedDeck.identity`
     - :class:`SpecularReemission`
     - a **surface** returning :math:`\alpha` specularly
   * - ``AlbedoBoundary(α, IsotropicReturn)``
     - :meth:`SelfPairedDeck.identity`
     - :class:`LambertianReemission`
     - a surface returning :math:`\alpha` diffusely
   * - ``VacuumInflow``
     - :meth:`SelfPairedDeck.identity`
     - :math:`0`
     - a surface returning nothing

:class:`SpecularReturn` and the mirror deck element
(:meth:`SelfPairedDeck.mirror`) realize to the *same* permutation and are
nonetheless different types. That is the **point**, not a smell: two types make
"put a surface's response in the geometry slot" unspellable, which is the exact
error this section corrects.

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
and preserves :math:`|\Omega\cdot\hat n|`. So for a **quotient** law the
crossing is not something the *physics* does — it is something the *geometry*
provides, and the **classifying** typing is

.. math::

    G : \Gamma_+ \to \Gamma_-, \qquad R : \Gamma_- \to \Gamma_-.

An earlier draft of this module typed ``G`` as an **endomorphism of the outflow
trace** (:math:`\Gamma_+ \to \Gamma_+`), which forces the mirror *into*
:math:`R` for a specular law — the same conflation running the other way.

⛔ **And the typing above is a CLASSIFICATION, not the realized one — this
module presented it as both until 2026-08-04.** The sentence above generalises
an argument about the mirror to every law, but a wall is not a quotient, so at a
white or albedo face there is no isometry to provide the crossing and the
**response** performs it: `[M]` every realized response is
:math:`\Gamma_+ \to \Gamma_-`, never an endomorphism of :math:`\Gamma_-`. Read
the display as *"exactly one of the two is non-trivial, and it is the one that
crosses"* — which for :math:`G` present-and-non-trivial reduces to exactly the
sentence above. How each factor is EVALUATED is a separate question with a
per-kind answer (:class:`BoundaryGeometryMap` atomic,
:class:`BoundaryResponseKernel` composed); conflating classification with
factorization is what let this stand.

Why the misassignment survived — a theorem
------------------------------------------

If :math:`R` is rank-one, :math:`R = u \otimes v`, then
:math:`R \circ G = u \otimes (G^{\mathsf T} v)`. The Lambertian's
:math:`v = |\Omega\cdot\hat n|` is **invariant** under both the mirror and the
periodic translation, so :math:`G^{\mathsf T}v = v` **as a function** and the
composite does not depend on WHICH admissible :math:`G` was used:

    :math:`G` **is unobservable exactly when** :math:`R` **is rank-one.**

⛔ **This concluded** ":math:`R \circ G = R`" **until 2026-08-04, and that does
not type-check.** With :math:`G : \Gamma_+ \to \Gamma_-` the left side is
:math:`\Gamma_+ \to \Gamma_-` while :math:`R` is (classified)
:math:`\Gamma_- \to \Gamma_-`; the step silently identified :math:`G^{\mathsf T}v`
with :math:`v` by treating :math:`v` as a *function* without tracking which
half-trace it is restricted to. Harmless before **B3.2** made the halves distinct
spaces; a type error after. The theorem's content is unaffected — it is
:math:`G`-INDEPENDENCE of the composite, not an equality of operators.

⚠ **And the slogan's hypothesis is stronger than "rank-one".** The work is done
by :math:`v` being :math:`G`-**invariant**; a rank-one :math:`R` with
:math:`v = \delta_{\Omega_0}` makes :math:`G` fully observable. Rank-one is what
makes the *Lambertian's* :math:`v` a measure, and the measure is what :math:`G`
preserves. The slogan is kept verbatim because five code and test sites
cross-reference it; :ref:`bc-factor-roles` carries the sharpened form.

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
method-agnostic law does not have; ``SelfPairedDeck.mirror(axis="x")`` needs
nothing.

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

The **geometry tier** splits on the **pairing**, not on the individual motion.
:class:`SelfPairedDeck` is the face paired with *itself*
(:meth:`~SelfPairedDeck.domain_face` answers ``face``); :class:`SpatialWrap` is
the face paired with a *distinct* one (it answers the opposite face). Those two
are genuinely non-isomorphic — the self-paired guard rejects a translation
outright — so the tier earns exactly two types.

The trivial pairing and the mirror are **not** two more. They differ by a
*value*, not by structure: same guard, same ``domain_face``, and
:attr:`~SelfPairedDeck.permutes_ordinates` *derived* from the motion instead of
declared. So they are two constructors of one type
(:meth:`~SelfPairedDeck.identity`, :meth:`~SelfPairedDeck.mirror`). They were
two classes — ``IdentityMap`` and ``SpecularMirror`` — until phase **G5**, each
hand-declaring ``permutes_ordinates`` and ``is_adjointable``: a second source of
truth for two questions the deck element already answers about itself, the
second of them a theorem its construction guard proves.

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

import numpy as np

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.face_layout import AXIS_NAMES, face_normal, face_opposite


__all__ = [
    "BoundaryGeometryMap",
    "BoundaryResponseKernel",
    "IsotropicReturn",
    "LambertianReemission",
    "ReemissionClosure",
    "ScalarResponse",
    "SelfPairedDeck",
    "SpatialWrap",
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
    :ref:`bc-factor-roles` on the theory page.

    ⛔ This paragraph continued *"the crossing from* :math:`\Gamma_+` *to*
    :math:`\Gamma_-` *is part of this factor, not of the response"* until
    2026-08-04. That is proven for the case it argues — a specular mirror IS
    the unique ambient isometry fixing the face, so for a **quotient** law the
    geometry provides the crossing — but it was stated for *every* law, and it
    is false for a law with no isometry. A wall is not a quotient, so nothing
    provides the crossing geometrically and the **response** does it (see
    :class:`BoundaryResponseKernel`). The honest general form: **whichever
    factor is non-trivial carries the crossing**, which is well defined because
    exactly one of them ever is.

    This map, being a genuine deck transformation, is always the non-trivial
    factor when it is present at all — so for :math:`G` specifically the
    crossing IS geometric, and :meth:`domain_face` below is that statement at
    the level of which face's :math:`\Gamma_+` is consumed.

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
        inverse the transpose. Every map here answers ``True``; the declaration
        survives because a *future* map could be realized without one, and
        because :class:`SpatialWrap` answered ``False`` until campaign phase
        **B3.4c** built its pushforward — it was reporting an implementation
        gap (#183), never a property of the map.
        """
        ...

    def domain_face(self, face: str) -> str:
        r"""Which face's :math:`\Gamma_+` this map's domain is, installed on ``face``.

        Since **B3.2** an SN boundary law is typed :math:`\Gamma_+ \to
        \Gamma_-`; this names *whose* :math:`\Gamma_+`. The answer is ``face``
        itself for every map that acts within one face, and a DIFFERENT face
        exactly when the deck transformation carries the crossing between two
        faces of the fundamental domain — which today means
        :class:`SpatialWrap` alone.

        **Why the geometry tier owns this and the response tier cannot.** A
        response is *constitutive* — a property of the material surface at this
        face — so it answers to what arrives at this face and structurally
        cannot reach another one. Crossing between faces is an act of the deck
        group, which is B3.0's ruling (":math:`G` carries the crossing") read
        one level up: :math:`G` carries the crossing in ANGLE for a mirror and
        in SPACE for a wrap, and both are the same statement about which
        :math:`\Gamma_+` the law consumes.

        **Why it takes the installation face as an argument** rather than being
        a stored field: which face is the partner depends on where the law is
        installed — configuration — while the axis is intrinsic. That is the
        same B0 rule that kept ``partner_face`` off :class:`SpatialWrap`.

        A map whose declared axis cannot be reconciled with ``face`` MUST raise
        rather than answer: a wrap along ``y`` installed on an ``x`` face
        describes no identification at all, and silently returning ``face``
        would realize it as a bare identity on the wrong half-trace.
        """
        ...


@runtime_checkable
class BoundaryResponseKernel(Protocol):
    r"""The **constitutive response** — classified :math:`R`, realized
    :math:`\Gamma_+ \to \Gamma_-`.

    How much of the outgoing flux returns, and with what angular distribution.
    This is where every piece of *physics* in a boundary law lives: absorption,
    accommodation, diffusivity. It is a positive kernel operator, and — the
    discriminator — it is **not multiplicative**.

    .. important::

       ⭐ **The** :math:`R \circ G` **split is a TAXONOMY, not a computational
       factorization** — and this docstring typed :math:`R` as
       :math:`\Gamma_- \to \Gamma_-` until 2026-08-04, which conflated the two.

       As a *classification* the pair :math:`G : \Gamma_+ \to \Gamma_-`,
       :math:`R : \Gamma_- \to \Gamma_-` is coherent, and it answers the only
       question the taxonomy exists to answer: *is this law's content geometry
       or physics?* (multiplicativity + the quotient test — see
       :ref:`bc-factor-roles`).

       But **no realized response is an endomorphism of** :math:`\Gamma_-`.
       :class:`LambertianReemission` realizes as
       :class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator` composed on
       :class:`~orpheus.sn.boundary.angular.PartialCurrentOperator`, typed
       :math:`\Gamma_+ \to \Gamma_-`; :class:`SpecularReemission` realizes as a
       narrowed permutation, likewise :math:`\Gamma_+ \to \Gamma_-`;
       :class:`ScalarResponse` realizes as a commuting scale. **Whichever
       factor is non-trivial carries the crossing** — and since exactly one of
       them ever is, that is well defined. For a *quotient* law the crossing is
       geometric; for a *constitutive* law the physics does the crossing, by
       integrating the outgoing flux and re-emitting an incoming one. There is
       no ambient isometry at a wall to provide it.

       **How a response is actually evaluated** follows its kind: it is
       :math:`N` composed operations :math:`\Gamma_+ \to \dots \to \Gamma_-`.
       The Lambertian is an outflow **angle contraction** followed by an
       **isotropic broadcast**, through an intermediate state in the
       angle-integrated per-face scalar-current space — the cosine-weighted
       mean outgoing intensity, a real physical quantity. Contrast a
       :class:`BoundaryGeometryMap`, which is **atomic**: a measure-preserving
       bijection does not factor into two meaningful pieces.

       Why it matters rather than being pedantry: the adjoint. A deck
       transformation's is a *theorem* (:math:`G^{-1} = G_{g^{-1}}`, and
       measure-preservation makes that inverse the transpose). A response's
       exists iff the intermediate space carries a **non-degenerate** metric —
       and then :math:`R^* = C^*B^* = G_+^{-1}R^{\mathsf T}G_-`, with the
       intermediate metric CANCELLING (`[M]` verified over eleven orders of
       magnitude of that metric; broken only when it is singular).

    **Three** realizations, genuinely non-isomorphic: :class:`ScalarResponse`
    (a bare amplitude, the whole story on a scalar trace where the angular
    distribution has no degrees of freedom), :class:`LambertianReemission` (a
    rank-one angular kernel), and :class:`SpecularReemission` (a permutation —
    a polished wall's specular return, which is constitutive rather than
    geometric because a wall is not a quotient). The theory page anticipated
    exactly this tier as *"a full angular kernel in general weak-form BCs
    (deferred)"*; B3 mints it.

    (This paragraph read "Two realizations" and named only the first two until
    2026-08-04. **B3.4b** added the third and updated the module docstring's
    count at line 237 without updating this one — so the module asserted three
    and this Protocol asserted two, five hundred lines apart. Found while
    retiring a different symbol; a dead reference is a tripwire for a false
    claim, and this one was found by looking one line further.)

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
        r"""Whether the realized kernel exposes an honest transpose.

        ``True`` for every shipped kernel since **G6.3 step 3** (2026-08-04).

        It read *"``False`` for the Lambertian today … phase B5 types that
        kernel as the rank-one it is, at which point the adjoint becomes
        structurally available and this flips."* G6.3 absorbed that step by
        FACTORING the realization — :class:`PartialCurrentOperator` then
        :class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator` — which
        removes the two-``.T``-semantics ambiguity instead of resolving it:
        each link has one honest transpose, so the composite's Hilbert adjoint
        follows from its bound spaces.

        ⛔ This read *"its realized form is self-adjoint under the
        cosine-weighted inner product"* until 2026-08-04. **Type-incoherent
        since B3.4a**: self-adjointness requires domain :math:`=` codomain, and
        the realized Lambertian maps :math:`\Gamma_+ \to \Gamma_-` over disjoint
        index sets. It was true of the pre-B3.4a full-face endomorphism. What
        survives is a **structural symmetry** — :math:`R` and :math:`R^*` share
        one form with the half-traces exchanged — see
        :class:`~orpheus.sn.boundary.angular.PartialCurrentOperator` and its
        partner for the measured statement.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# Concrete deck transformations
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class SelfPairedDeck:
    r""":math:`G` for a face paired with **itself** — Poincaré's self-pairing case.

    A boundary law's :math:`G` is a **face pairing** of the fundamental domain
    (Poincaré). Two cases exist and they are structurally different:

    - **Self-paired** — the face is paired with itself, so
      :math:`\mathrm{dom} = \mathrm{cod}` and ``domain_face(face) is face``.
      This class. Its two inhabitants are the trivial pairing
      (:meth:`identity`, every constitutive law) and the mirror
      (:meth:`mirror`, a reflecting face).
    - **Paired with a DISTINCT face** — the codomain of one surface is the
      domain of its partner and vice versa. That is :class:`SpatialWrap`
      (periodic, a translation) and the unbuilt sector BC (#178, a rotation).
      **Not this class**: a pair cannot be named by one face, so it needs a
      surface-pair type that does not exist yet.

    The construction guard is the self-pairing condition itself, and it is what
    makes the distinction structural rather than documented:

    .. math::

        Q \text{ linear} \;\wedge\; \dim \mathrm{Fix}(Q) \ge d - 1

    **Involution is a THEOREM here, not a guard.** A linear orthogonal
    :math:`Q` fixing a hyperplane pointwise is :math:`I` or a reflection, so
    :math:`Q^2 = I` follows. Guarding on ``element_order() in (1, 2)``
    *instead* would be the converse of the insight and is **wrong**:
    :math:`E(3)` has FOUR involution families, and the other two — the
    half-turn (``order=2, det=+1, fix=1``) and the inversion
    (``order=2, det=-1, fix=0``) — map a face to its **opposite**, which is
    precisely :class:`SpatialWrap`'s deferred job. An order-only guard admits
    exactly the elements this type exists to exclude. It would also carry an
    ``atol`` into a type invariant, which the fixed-set spelling does not need.

    **Why the guard forbids an affine part** (``is_linear``), and this is the
    load-bearing half: :meth:`~orpheus.geometry.transformation.RigidMotion.on_directions`
    **drops the translation**, so a mirror's *offset* — where the plane sits —
    is bit-identically invisible at every level the tree can measure: the
    ordinate permutation, the realized image, the frozen snapshots, any scalar.
    A wrong mirror position is `vv-principles` Mode-12 **designed-green**, and
    no gate at any tolerance can ever catch it. Refusing the affine part makes
    the error *unspellable* instead of ungatable — the only closure available.
    (`[M]` ``on_directions`` identical at ``offset ∈ {0, 2.5, −17.0}``.)

    The angular/spatial split is the same statement one level up: :math:`G`
    carries the crossing **in ANGLE** for a self-paired face (the linear part
    acts, the position cannot matter because the face maps to itself) and
    **in SPACE** for a genuine pair (the translation is the whole content).
    """

    motion: RigidMotion

    def __post_init__(self) -> None:
        motion = self.motion
        if not isinstance(motion, RigidMotion):
            raise TypeError(
                f"SelfPairedDeck takes a RigidMotion, got {type(motion).__name__}. "
                f"The deck element is the transformation itself, not a name for "
                f"one — an axis LETTER cannot say which plane, and the letter "
                f"table had three live spellings when this type was minted."
            )
        if not motion.is_linear:
            raise ValueError(
                "a self-paired deck element must be LINEAR (no translation). "
                "A face paired with itself is fixed by its own pairing, so "
                "there is nothing for a translation to do — and `on_directions` "
                "drops it, which makes a wrong mirror POSITION invisible to "
                "every gate at every tolerance (vv-principles Mode 12). "
                "A translation belongs to a genuine face PAIR: see SpatialWrap."
            )
        fixed = motion.fixed_subspace_dimension
        if fixed < motion.dimension - 1:
            raise ValueError(
                f"a self-paired deck element must fix the face pointwise, i.e. "
                f"dim Fix >= {motion.dimension - 1}; this one fixes {fixed}. "
                f"An element fixing less carries the face to a DIFFERENT face "
                f"(a half-turn fixes a line, an inversion only the origin) — "
                f"that is a genuine face pairing and needs a surface-pair type, "
                f"not this one. Being an involution is NOT sufficient: E(3) has "
                f"four involution families and two of them are face-swapping."
            )

    @classmethod
    def identity(cls, dimension: int = 3) -> "SelfPairedDeck":
        r"""The trivial pairing — :math:`G = \mathrm{id}`, the deck group's unit.

        The law identifies no geometry: the domain is represented by itself,
        not by a quotient. Every law whose content is entirely constitutive
        lands here — vacuum, prescribed inflow, white, albedo, zero-flux — and
        for all of them the crossing to :math:`\Gamma_-` is the face's own
        canonical one, which the realizer supplies.

        That this is *sound* rather than merely conventional is the rank-one
        theorem (:ref:`bc-factor-roles`): where the response destroys
        directional information, the composite :math:`R \circ G` is the **same
        for every** measure-preserving :math:`G`, so declaring the identity is
        not a lossy choice — it is the honest statement that the law fixes no
        geometry.

        (This read ":math:`R \circ G = R` for **any** measure-preserving
        :math:`G`" until 2026-08-04. Same non-type-checking shorthand as the
        module docstring's: the composite is :math:`\Gamma_+ \to \Gamma_-` and
        cannot equal a :math:`\Gamma_- \to \Gamma_-` classification. The
        soundness argument is untouched — :math:`G`-independence is exactly what
        it needs.)
        """
        return cls(RigidMotion.identity(dimension))

    @classmethod
    def mirror(cls, *, axis: str, dimension: int = 3) -> "SelfPairedDeck":
        r"""The reflection whose plane is normal to ``axis``.

        The unique ambient isometry fixing the face. It exchanges the
        hemispheres — which is why :math:`G`, not the response, carries the
        :math:`\Gamma_+ \to \Gamma_-` crossing — and preserves
        :math:`|\Omega\cdot\hat n|`, so it is measure-preserving in the trace
        measure (ERR-042).

        Quotient reading: a deck transformation with **fixed points** (the
        mirror plane), so a reflecting face makes the computational domain an
        orbifold with what Thurston calls a *reflector boundary*. Contrast
        :class:`SpatialWrap`, whose action is **free** — that contrast IS the
        self-paired/genuinely-paired split this type is built on.

        The axis letter is resolved against
        :data:`~orpheus.numerics.face_layout.AXIS_NAMES` — the single home of
        that convention — and refused HERE rather than at realization, so a
        mis-named axis cannot be constructed and then fail late somewhere
        holding a quadrature.
        """
        try:
            index = AXIS_NAMES.index(axis)
        except ValueError:
            raise ValueError(
                f"axis must be one of {AXIS_NAMES}, got {axis!r}. The mirror "
                f"names the plane's NORMAL; validating here is what keeps a "
                f"bad axis from surviving construction and failing later "
                f"against a quadrature it never should have reached."
            ) from None
        if not 0 <= index < dimension:
            raise ValueError(
                f"axis {axis!r} is index {index}, out of range for a "
                f"{dimension}-dimensional deck element"
            )
        normal = np.zeros(dimension)
        normal[index] = 1.0
        return cls(RigidMotion.reflection(normal=normal))

    @property
    def permutes_ordinates(self) -> bool:
        """Whether realizing this map permutes the ANGULAR index.

        **Derived, not declared.** The linear part IS the action on
        directions, so "does it move angle?" is "is the linear part the
        identity?" — a question the motion answers about itself. Before this
        type the answer was a hand-written ``True``/``False`` per class, which
        is a second source of truth for something the element already knows.
        """
        return not np.array_equal(
            self.motion.linear, np.eye(self.motion.dimension)
        )

    @property
    def is_adjointable(self) -> bool:
        r"""Whether the realized map exposes an honest transpose.

        A **theorem** for a genuine deck transformation: the composition
        operator of a bijection :math:`g` is invertible with
        :math:`G^{-1} = G_{g^{-1}}`, and measure-preservation makes that
        inverse the transpose. Always ``True`` here — and unlike the
        per-class declaration it replaces, it cannot drift, because the
        construction guard already proved the element is an involution.
        """
        return True

    def domain_face(self, face: str) -> str:
        r"""Which face's :math:`\Gamma_+` this map's domain is.

        ``face`` itself — **by construction, not by choice**. That is the
        definition of self-paired, and it is why this is one line on the type
        rather than a per-class decision that could disagree with the guard.

        No axis-vs-face reconciliation happens here. A mirror's axis IS
        checked against the installation face, but at realization, where a
        mismatch is diagnosed against the actual reflection table (a mirror
        about ``y`` on an x-face relabels WITHIN each half-trace instead of
        exchanging them) and where the quadrature this method does not hold is
        available. A second, weaker name-only test here would be a twin that
        can disagree with it.
        """
        return face




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
    :meth:`SelfPairedDeck.mirror`'s ``axis``. The realizer derives the partner
    from the installation face plus this axis.

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
        # The translation IS adjointable as a map — for a deck transformation
        # that is a theorem (see the Protocol's note), and the realized form is
        # now a genuine one: an identity between the PARTNER's Γ₊ and this
        # face's Γ₋, whose transpose scatters the image back over the partner's
        # Γ₊. Until campaign phase **B3.4c** this answered ``False`` to report
        # the implementation gap (#183): the operator existed but was fed its
        # own face's outflow, so there was no partner leg for a transpose to
        # return along. B3.4c built the channel and this flipped, WITH its gate.
        return True

    def domain_face(self, face: str) -> str:
        r"""The PARTNER face — the one datum that makes periodic periodic.

        A wrap's law reads :math:`\gamma_-\psi|_f = \gamma_+\psi|_{f'}`: the
        inflow here is the outflow *across the domain*. It is the only
        :class:`BoundaryGeometryMap` whose answer is not ``face``, which is
        exactly the sense in which a torus is a quotient and a wall is not.

        The identification is well-posed because the two faces' outward normals
        are opposite (:math:`\hat n_f = -\hat n_{f'}`): a direction OUTGOING at
        :math:`f'` is INCOMING at :math:`f`, so :math:`\Gamma_+(f')` and
        :math:`\Gamma_-(f)` are the same index set and the crossing costs
        nothing. `[M]` measured 2026-08-01 as an exact set equality on
        ``gauss_legendre(8)``, ``product(2,4)``, ``level_symmetric(6)`` and
        ``lebedev(17)``, on both axis pairs. The realizer ASSERTS it rather
        than assuming it — the user's B3.4c ruling was that the quotient
        reading becomes a guard, not a mesh restructure.

        Raises
        ------
        BoundaryError
            If ``face`` does not lie on :attr:`axis`. A wrap along ``y``
            installed on an ``x`` face identifies nothing — the translation
            :math:`x \mapsto x + L_y` does not carry that face anywhere — and
            answering ``face`` would silently realize it as a bare identity on
            the wrong half-trace, which is the ERR-041 mis-declaration class
            (the diffuse arm's orientation cross-check is the same shape).
        """
        from ._errors import BoundaryError

        try:
            axis_index, _sign = face_normal(face)
        except ValueError as exc:
            raise BoundaryError(
                f"SpatialWrap(axis={self.axis!r}) cannot name a partner for "
                f"face {face!r}: {exc}",
                law="periodic",
            ) from exc
        if AXIS_NAMES[axis_index] != self.axis:
            raise BoundaryError(
                f"SpatialWrap declares axis={self.axis!r} but is installed on "
                f"face {face!r}, which lies on axis "
                f"{AXIS_NAMES[axis_index]!r}. A wrap identifies the two faces "
                f"NORMAL to its own axis; it carries an {self.axis}-face to "
                f"the opposite {self.axis}-face and says nothing about any "
                f"other face, so there is no partner to name here. Declare "
                f"PeriodicBoundary(axis={AXIS_NAMES[axis_index]!r}) if this "
                f"face's axis is the periodic one.",
                law="periodic",
            )
        return face_opposite(face)


# ═══════════════════════════════════════════════════════════════════════
# Concrete response kernels
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class ScalarResponse:
    r""":math:`R = \alpha\,I` — a scalar amplitude, no angular redistribution.

    The flux returns attenuated but otherwise untouched, so whatever angular
    structure the deck transformation produced is what arrives. Paired with
    :meth:`SelfPairedDeck.mirror` this is the specular albedo; paired with
    :meth:`SelfPairedDeck.identity` on a **scalar** trace it is the
    partial-current albedo, where it is the complete story because the angular
    distribution has no degrees of freedom to fix.

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
        # ⭐ FLIPPED 2026-08-04 at G6.3 step 3 (#330), which absorbed the B5
        # step this comment was waiting on. It read FALSE, "advertising the
        # unweighted transpose would invite two different ``.T`` semantics …
        # B5 types it as ``u ⊗ v`` (transpose ``v ⊗ u``) … this flips there,
        # WITH its gate."
        #
        # Factoring the realization into PartialCurrentOperator @
        # IsotropicEmissionOperator REMOVED the ambiguity rather than resolving
        # it: each link has ONE honest transpose (an outer product, a sum), so
        # there is no longer a choice of ``.T`` semantics to avoid. The
        # composite's Hilbert adjoint then follows from the bound spaces.
        #
        # The flip is REQUIRED, not cosmetic: TestFactorAdjointabilityMatches-
        # TheRealizedOperator holds this declaration against what the realizer
        # actually produces, and it caught the disagreement the moment the
        # realization gained its transpose — a landed capability staling its
        # own deferral contract.
        #
        # (Also said "the realized operator is self-adjoint under the
        # cosine-weighted inner product" until 2026-08-04. Type-incoherent
        # since B3.4a narrowed it to Γ₊ → Γ₋: self-adjointness needs
        # domain == codomain. What holds is that R and R* share ONE FORM with
        # the half-traces exchanged — see the property docstring above.)
        return True


@dataclass(frozen=True, slots=True)
class SpecularReemission:
    r""":math:`R = \alpha\,P_{\text{mirror}}` — **specular** re-emission.

    A surface that returns :math:`\alpha` of the arriving flux into the
    mirror-partner direction. Realized as the ordinate permutation
    ``quadrature.reflection_index(axis)``, scaled by :math:`\alpha` — the same
    matrix :meth:`SelfPairedDeck.mirror` realizes to.

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
        # TRUE for a STRUCTURAL reason rather than an implementation one:
        # (this read "TRUE, unlike the Lambertian's" until 2026-08-04. The
        # contrast expired in the SAME commit ~70 lines up, where the
        # Lambertian flipped to True once its realization was factored. Both
        # are adjointable now; only the REASONS still differ, which is what
        # the rest of this comment is about.) a scaled permutation's transpose is the
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
    :meth:`SelfPairedDeck.identity` (the geometry of every surface law) supplies
    none either. Something must say which outgoing direction feeds which
    incoming one; that something is this closure.

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

    Deliberately a different type from the mirror deck element
    :meth:`SelfPairedDeck.mirror` builds, which realizes to the same
    permutation. The two realize identically and are semantically disjoint, and
    keeping them apart is what makes "a wall's response in the geometry slot"
    unspellable — the exact conflation phase B3 found in the white law and the
    user's 2026-08-01 ruling generalized.
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
