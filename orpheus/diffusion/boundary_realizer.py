r"""The diffusion realizer table: ``DIFFUSION_REALIZERS`` (``boundary_realizer.py`` of ``orpheus/diffusion``).

**What this module is.** The ONE table through which diffusion admits and realizes a declared exterior. It replaces
``DiffusionBoundaryRealizer`` (a ``response_kernel.amplitude`` read behind two ``isinstance`` refusals, written for the
string-tagged boundary laws) with a :class:`~orpheus.transport.boundary_conditions.RealizerTable` keyed by UNIT and filled by six
realizer classes that register by class keyword; the table is sealed on the last line of this module, so anyone holding
``DIFFUSION_REALIZERS`` has executed every registration. It is owned by its method-mesh: :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`
holds it as ``boundary_realizers`` and the shared :func:`~orpheus.transport.method.realize_boundary` reads it.

Every diffusion boundary condition is an albedo
===============================================

Diffusion's boundary state is the scalar partial-current pair :math:`(J^+, J^-)` on the
:class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`: ONE outflow degree of freedom and ONE inflow degree of
freedom per face per group. A linear homogeneous boundary law therefore collapses to a single scalar response per face per group,

.. math::

    J^- \;=\; \mathcal{A}\, J^+ ,

and the realization of a unit is the composition of two maps: the unit to its scalar :math:`\mathcal{A}`, then the scalar to the
operator (:func:`_albedo_operator`, the structure-keyed collapse: ``ZeroOperator`` at :math:`\mathcal{A}=0`, ``IdentityOperator``
at :math:`\mathcal{A}=1`, :math:`\mathcal{A}\cdot I` otherwise; the operator is face-slot-shape-agnostic, a scalar multiple of the
identity on the ``(ng, *face_spatial)`` outflow data). The operator maps the outflow half of a face slot to the inflow half and is stamped
:attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` (the ``A_ss`` leaf of the composite state, consumed by
:class:`~orpheus.diffusion.operators.DiffusionBoundaryOperator`).

**Diffusion annihilates the shape.** It sees :math:`\Pi_{P_1}\colon (s, d) \mapsto s + d` on the response network of a wall: the
angular redistribution that tells a specular return from a Lambertian one (a permutation against a cosine-weighted average in SN)
is integrated out of the half-range :math:`\ell=0` moments. Both chains therefore realize as the scalar WEIGHT times the identity,
and a sum of chains is the scalar sum, which is law L-L (realization is a functor on the response network) read on this table:
the image of a mixed wall is :math:`\mathcal{A} = s + d`. The realization of one chain is EXACT and not faithful: this is the
tier-2 case of :ref:`bc-method-realizability` (the projection cannot tell the laws apart, which is a different statement from
approximating any of them), and it is stated as ``Realized(faithful=False)``.

The table (unit, admission, image)
==================================

=========================  =====================  ===========================  ========================================
Realizer                   Key (unit)             Admission                    Image
=========================  =====================  ===========================  ========================================
``DiffusionZero``          ``ZeroNetwork``        ``Realized(True)``, tier 1   ``ZeroOperator`` (:math:`\mathcal{A}=0`)
``DiffusionSpecularChain`` ``SpecularChain``      ``Realized(False)``, tier 2  :math:`w\cdot I` (the chain's weight)
``DiffusionDiffuseChain``  ``DiffuseChain``       ``Realized(False)``, tier 2  :math:`w\cdot I` (the chain's weight)
``DiffusionZeroFlux``      ``ZeroFluxCondition``  ``Realized(True)``, tier 1   :math:`-I` (:math:`\mathcal{A}=-1`)
``DiffusionSource``        ``InflowTerms``        ``Refused(PLUMBING, #656)``  none (a source is an array, no operator)
``DiffusionDecks``         the decks              mirror ``Realized(False)``,  mirror: :math:`I`; translation: refused
                                                  translation ``Refused(TOPOLOGY)``
=========================  =====================  ===========================  ========================================

* **Vacuum is tier 1, not an approximation (measured).** The empty network (vacuum, ``specular(0)``, ``white(0)``) has
  :math:`\mathcal{A} = 0` exactly, and a scalar commutes with every linear projection, so the projection loses nothing: the
  admission is ``Realized(True)``. ``[M]`` ``measure_diffusion_vacuum.md``: of nine law spellings the retired realizer was called
  on, seven realized, and the deciding line was ``return response.amplitude``. The Marshak closure (:math:`J^\pm = \phi/4 \pm J/2`
  with Fick's law, the classic Robin form :math:`\phi + 2D\,\nabla\phi\cdot\hat n = 0`; Marshak 1947, Bell & Glasstone 1970 section 3.4)
  is the METHOD's :math:`P_1` closure, applied identically to every albedo, so it belongs to the method and not to the admission of the
  vacuum wall. DRAFT 1 of the specification listed ``Approximated("Marshak")`` here; that was refuted by this measurement.
* **The two chains are tier 2.** ``specular(0.4)`` and ``white(0.4)`` realize the same operator :math:`0.4\,I`; the
  pair is the witness that diffusion cannot tell them apart (the specification's tier rule would call ``specular(0.4)`` tier 1, since
  the angular map is recovered, and the table states the lower tier because diffusion cannot see the shape).
* **Zero flux is tier 1 and diffusion-only** (:class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition`, :math:`\mathcal{A} = -1`,
  the Dirichlet idealization :math:`\phi_\Gamma = 0` under the :math:`P_1` dictionary :math:`\phi_\Gamma = 2(J^+ + J^-)`). It is outside the
  sub-Markov range, which is why the table's :class:`~orpheus.transport.boundary_conditions.StateCone` is ``SIGNED_CURRENT``.
* **A mirror deck is tier 2 with** :math:`\mathcal{A} = 1`: a mirror plane is a symmetry, whose image on a per-face scalar is the
  identity, the same operator as ``white(1.0)`` and ``specular(1.0)``.

Refusals
========

* **A source is plumbing** (``Refused(PLUMBING, #656)``). :math:`J^- = q` is the AFFINE part of the trace map: it is a vector of
  the inflow half-trace, which here is one number per face, and the trace carries it perfectly well. Nothing about the :math:`P_1`
  projection makes it unrepresentable, so the refusal disappears the day the fixed-source arm is wired, with no theory changing
  (#656 tracks diffusion's affine channel; the closed #290 is not the issue). ``DiffusionSource`` has no ``realize`` and no
  ``sample``: the door of :meth:`~orpheus.transport.boundary_conditions.RealizerTable.sample_inflow` raises before any body would run.

  Do not rebuild the SN design this used to cite: SN realizes the LINEAR factor of an inflow wall (the zero map) while :math:`q`
  travels the boundary-source channel (:ref:`bc-affine-source-channel`). If a diffusion source arm is built here, build it as a
  source on the trace block, never as an operator that ignores its input.
* **A translation deck is topological** (``Refused(TOPOLOGY)``): periodic's deck transformation is a translation, which acts on
  the spatial coordinate and leaves angle untouched, so P1 loses nothing angularly. What blocks it is that a PAIRED deck couples a
  face to a DIFFERENT face (the opposite one, for the wrap), and this table's codomain is a per-face scalar, a block-diagonal
  object with no slot for cross-face coupling. A method could resolve angle perfectly and still refuse this.
* **Any other condition** (another method's own condition) meets the table's unregistered-unit door,
  ``Refused(PLUMBING, #659)``: a signed-current table never refuses on the state cone (a negative return fraction is inside its cone).

**An albedo above one is unspellable.** The retired string registry realized ``albedo = 5.0`` (``k = 2.972`` on the two-region two-group
slab at ``5d04c296``). The shared wall's cone parse (``WallCondition`` refuses a return fraction above one at construction) removes that
spelling before it reaches this table: the albedo item of #265 (the other items of that issue are untouched, and it stays open).

**The retired door.** ``DiffusionBoundaryRealizer`` and ``_partial_current_albedo`` are gone (the amplitude read and the two
``isinstance`` refusals are the table's admissions now), together with ``DiffusionMesh.BOUNDARY_OPERATOR_REGISTRY`` (a string table
of admitted law kinds), in the same commit that mints this table.

References
----------

* Marshak, R. E. (1947). "Note on the spherical harmonic method as applied to the Milne problem for a sphere". *Phys. Rev.* 71, 443.
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*. Van Nostrand Reinhold. Section 3.4 (diffusion-theory boundary conditions).
* ``.claude/plans/diffusion_crosswalk.md``: the :math:`(J^+, J^-)` convention contract.
* ``.claude/plans/wave2_boundary_carve.md`` sections 2.7 and 4.2: the table and its pin.

Layer: ``diffusion`` (L3). Imports ``transport`` (the table machinery), ``geometry`` (the deck vocabulary), ``numerics`` (the operator
leaves) and its siblings.
"""

from __future__ import annotations

from typing import assert_never

from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.diffusion.zero_flux import ZeroFluxCondition
from orpheus.geometry.deck import DeckTransformation, FacePairing, MirrorDeck, TranslationDeck
from orpheus.numerics.operator import IdentityOperator, LinearOperator, ZeroOperator
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.boundary_conditions import (
    DeckRealizer,
    DiffuseChain,
    FaceSite,
    InflowTerms,
    Realized,
    RealizerTable,
    RefusalAxis,
    Refused,
    ResponseRealizer,
    SourceRealizer,
    SpecularChain,
    StateCone,
    ZeroNetwork,
    stamp_boundary_role,
)

__all__ = [
    "DIFFUSION_REALIZERS",
    "DiffusionDecks",
    "DiffusionDiffuseChain",
    "DiffusionSource",
    "DiffusionSpecularChain",
    "DiffusionZero",
    "DiffusionZeroFlux",
]


def _albedo_operator(albedo: float) -> LinearOperator:
    r"""The structure-keyed collapse of :math:`J^- = \mathcal{A} J^+`.

    Mirrors the SN realizer's scalar branch exactly: :math:`\mathcal{A} = 0` collapses to the structural :class:`ZeroOperator`,
    :math:`\mathcal{A} = 1` to the structural :class:`IdentityOperator`, anything else to ``𝒜 · IdentityOperator()`` (the one-term
    combination). The operator is face-slot-shape-agnostic: a scalar multiple of the identity on the ``(ng, *face_spatial)``
    outflow data, whatever the face.
    """
    if albedo == 0.0:
        return ZeroOperator()
    if albedo == 1.0:
        return IdentityOperator()
    return float(albedo) * IdentityOperator()


#: Diffusion's realizer table. ``SIGNED_CURRENT``: the trace state is a partial current, so a negative return fraction (zero flux) is
#: inside the cone. ``unregistered_issue=659``: the table registers every shared unit and meets its unregistered-unit door only with
#: another method's own condition (#659 tracks the SN and diffusion tables' unregistered units). The source space is the trace's
#: inflow half-space, a plain ``FunctionSpace`` here (the source verb is refused before any sampling).
DIFFUSION_REALIZERS: RealizerTable[DiffusionMethodSpace, LinearOperator, FunctionSpace] = RealizerTable(
    "diffusion", state_cone=StateCone.SIGNED_CURRENT, unregistered_issue=659,
)


class DiffusionZero(ResponseRealizer, table=DIFFUSION_REALIZERS, realizes=ZeroNetwork):
    r"""The empty network: :math:`\mathcal{A} = 0` exactly, ``ZeroOperator``. Tier 1: a scalar commutes with every projection."""

    status = Realized(True)

    def realize(self, unit: ZeroNetwork, weight: float, site: FaceSite, space: DiffusionMethodSpace) -> LinearOperator:
        return stamp_boundary_role(_albedo_operator(0.0))


class DiffusionSpecularChain(ResponseRealizer, table=DIFFUSION_REALIZERS, realizes=SpecularChain):
    r"""A specular chain of weight :math:`w`: :math:`\mathcal{A} = w`. Tier 2: the shape (a permutation) is annihilated by :math:`\Pi_{P_1}`."""

    status = Realized(False)

    def realize(self, unit: SpecularChain, weight: float, site: FaceSite, space: DiffusionMethodSpace) -> LinearOperator:
        return stamp_boundary_role(_albedo_operator(weight))


class DiffusionDiffuseChain(ResponseRealizer, table=DIFFUSION_REALIZERS, realizes=DiffuseChain):
    r"""A diffuse chain of weight :math:`w`: :math:`\mathcal{A} = w`. Tier 2: the shape (a rank-one re-emission) is annihilated by :math:`\Pi_{P_1}`."""

    status = Realized(False)

    def realize(self, unit: DiffuseChain, weight: float, site: FaceSite, space: DiffusionMethodSpace) -> LinearOperator:
        return stamp_boundary_role(_albedo_operator(weight))


class DiffusionZeroFlux(ResponseRealizer, table=DIFFUSION_REALIZERS, realizes=ZeroFluxCondition):
    r"""Zero flux, :math:`\phi_\Gamma = 0`: :math:`\mathcal{A} = -1`, the operator :math:`-I`. Tier 1."""

    status = Realized(True)

    def realize(self, unit: ZeroFluxCondition, weight: float, site: FaceSite, space: DiffusionMethodSpace) -> LinearOperator:
        return stamp_boundary_role(_albedo_operator(unit.return_fraction))


class DiffusionSource(SourceRealizer[FunctionSpace], table=DIFFUSION_REALIZERS):
    r"""The source channel: refused. The affine part :math:`q` of the trace map has no channel in diffusion yet (#656, plumbing).

    No ``realize`` (a source is an array, never an operator) and no ``sample`` override (the table's door raises before the
    default body, which is inherited, would run). Lifts with no theory change: see the module docstring.
    """

    status = Refused(RefusalAxis.PLUMBING, "diffusion's affine boundary source has no channel: J^- = q is a vector of the inflow half-trace the scalar trace carries, and the fixed-source arm is not wired", 656)


class DiffusionDecks(DeckRealizer, table=DIFFUSION_REALIZERS):
    r"""The decks: a mirror plane is :math:`\mathcal{A} = 1` (tier 2); a translation wrap is refused on the topological axis.

    The match is exhaustive over :data:`~orpheus.geometry.deck.DeckTransformation`.
    """

    def admit(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: DiffusionMethodSpace | None = None) -> Realized | Refused:
        match deck:
            case MirrorDeck():
                return Realized(False)
            case TranslationDeck():
                return Refused(
                    RefusalAxis.TOPOLOGY,
                    f"the wrap couples {pairing.face!r} to its opposite face {pairing.partner!r}: a per-face scalar has no slot for cross-face coupling",
                )
            case _:
                assert_never(deck)

    def realize(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: DiffusionMethodSpace) -> LinearOperator:
        match deck:
            case MirrorDeck():
                return stamp_boundary_role(_albedo_operator(1.0))
            case TranslationDeck():
                raise AssertionError("a translation deck is refused at admission, so it is never realized")  # noqa: S101
            case _:
                assert_never(deck)


DIFFUSION_REALIZERS.seal()
