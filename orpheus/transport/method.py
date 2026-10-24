r"""TransportMethod — the Protocol over the method-mesh layer.

A **method-mesh** is a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
(the method-agnostic mesh + materials DATA) augmented with one
transport method's BEHAVIOR. Two exist:
:class:`~orpheus.sn.problem.SNProblem` (quadrature + sweep
machinery + angular trace) and
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` (scalar trace
+ realized boundary operators). :class:`TransportMethod` is the structural
Protocol naming what every method-mesh exposes *as a method*, and
:func:`realize_boundary` is the ONE generic body that turns the method's
:class:`~orpheus.transport.boundary_conditions.Binding` (the declared exterior
conditions joined to the carrier's tagged faces, and its decks) into realized boundary
operators through the method's :class:`~orpheus.transport.boundary_conditions.RealizerTable`.

The mathematics the body walks
==============================

A face of the carrier has exactly one role: glued to a partner (a *deck*, a symmetry of the domain) or facing the exterior
(a *tag*, an integer colour of the true boundary, the twin of a material id). The tag-to-condition map is supplied where the
method is known: ``bind(carrier, conditions)`` joins the conditions to the carrier's faces once, refusing an unbound tag, an
unused condition and a condition the carrier's colouring symmetry does not fix
(:mod:`orpheus.transport.boundary_conditions.binding`). Realization is then a partial map from a declaration on a face to an
operator on the method's trace space, and its domain is the admission table of the method: the SAME table answers "is this
declaration admitted on this face" and "realize it", so the query and the realization cannot disagree
(:mod:`orpheus.transport.boundary_conditions.realization`).

:func:`realize_boundary` visits EVERY face of the binding, including a face with no trace slot (the inner face of a hollow
radial body is a tag, the radial axis has no inner-surface trace). A face with a trace is realized through the method's own
:meth:`~TransportMethod.realize_boundary_role` hook, which builds the method space of the face and calls the table's
``realize`` (``admit`` first, the raise of ``BoundaryRefused`` iff refused). A face without a trace has no operator to build
and meets the table through ``require`` alone, the type- and site-level admission, so a declaration the method cannot honour
there is refused at construction like any other.

The base of the method-meshes
=============================

:class:`MethodMesh` is the one place the two binding duties of the method-meshes that ARE a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` live (the key extended by the declared conditions, and the keywords a
re-spelling passes so the rebuilt problem keeps them); ``SNProblem`` and ``DiffusionMesh`` derive from it, so neither spells either
duty itself. The methods whose mesh wraps a ``Mesh1D`` and is not a ``MaterialMesh`` (collision probability, characteristics, Monte
Carlo) have no key and no re-spelling, hence no such duty; what the three share instead is the bind-and-realize step and the two hooks
of the shared body, and :class:`Mesh1DMethodMesh` is the one place those live.

Genesis (#290 P7b), and what the Protocol deliberately does NOT contain
======================================================================

The Protocol was deliberately **deferred until a second method-mesh existed** (the defer-until-two rule). Its two recorded
witnesses were the **homogenization method-layer** of :mod:`orpheus.transport.mesh.material_mesh` and the twin
``SNProblem._resolve_bcs`` / ``DiffusionMesh._resolve_bcs`` loops, each doing the same tag-to-law-to-realize walk with only the
leaf realization differing.

* **Promotion classmethods.** ``SNProblem.from_material_mesh(mm, quadrature, scheme=...)`` vs
  ``DiffusionMesh.from_material_mesh(mm)``: the signatures differ because the methods genuinely need different injections.
  The Protocol is **instance surface only**.
* **Method-space construction.** ``SNMethodSpace.for_face(mesh, quadrature, face, trace)`` vs
  ``DiffusionMethodSpace.for_face(trace, face)``: non-unifiable for the same reason. The per-method
  :meth:`~TransportMethod.realize_boundary_role` hook subsumes it.
* **The trace spaces.** ``SNProblem.angular_trace`` and ``DiffusionMesh.scalar_trace`` are locus-family SIBLINGS, not one concept.
* **A dispatch name.** There is no ``method_name`` member and no string registry: a method is its realizer table, held by the
  method-mesh (``boundary_realizers``). The Wave-5 ``BoundaryRealizerRegistry`` (string-keyed ``method_name`` to realizer) was
  retired at #290 P7b: you hold the method-mesh, therefore you hold its realizer.

What the wave-2 carve changed here
==================================

``BOUNDARY_OPERATOR_REGISTRY`` (a ``dict[str, type[BoundaryTraceLaw]]``, a string table of admitted law kinds) is gone: its
two jobs are the binding (which tags a carrier has, which condition each carries) and the realizer table (which declarations
the method admits). ``resolve_boundary_conditions`` and its tag-parse helper ``_law_from_tag`` (which patched an axis and an
outward sign into a law per face) are gone with the orientation arguments they patched: a shared condition is face-local
(:class:`~orpheus.transport.boundary_conditions.WallCondition` holds no axis), and the realizer reads the face's normal from
the face name once. This module imports no ``BC`` and no boundary law type, which is what lets the ``orpheus.transport``
package import without the retired ``orpheus.geometry.boundary``.

Conformance is structural: nobody imports this to conform. ``orpheus.sn`` and ``orpheus.diffusion`` import
:func:`realize_boundary` (the shared body), **never the Protocol class**. Static conformance is checked where it matters: at
each mesh's ``realize_boundary(self)`` call, pyright verifies the mesh structurally satisfies ``TransportMethod[OpT]`` and infers
the method's own realized-value type ``OpT``. The realized value is whatever the method's face needs, and it is not always an
operator: SN realizes a :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`, diffusion a bare
:class:`~orpheus.numerics.operator.LinearOperator`, collision probability a closure over its cell volumes, characteristics a
track-linking rule, Monte Carlo a periodic-wrap marker. ``OpT`` is therefore unbounded. The ratchet covers the ``realize_boundary``
call of SN and diffusion; ``pyproject``'s ignore list holds ``orpheus/cp``, ``orpheus/moc`` and ``orpheus/mc``, so the conformance of
those three heads is not checked by CI.

Layer (``tests/gates/test_layer_imports.py``): L2 ``transport``. Imports ``numerics`` (the operator base), ``mesh`` (the face
inventory) and its sibling :mod:`orpheus.transport.boundary_conditions`. The method packages (L3) import THIS module; never the
reverse.

References
----------

* ``.claude/plans/archive/diffusion_integration_290.md`` section P7b (the carve that minted the Protocol).
* ``.claude/plans/wave2_boundary_carve.md`` section 2.10 (the shared body and the method hooks).
* :doc:`/theory/foundations/boundary_conditions` — the three-layer decomposition (trace structure / physical law / method
  realisation).
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, ClassVar, Generic, Protocol, TypeVar, runtime_checkable

from orpheus.mesh import Mesh1D, axes_from_legacy_mesh, face_labels
from orpheus.transport.boundary_conditions import FaceSite, bind
from orpheus.transport.mesh.material_mesh import MaterialMesh, conditions_key

if TYPE_CHECKING:
    from collections.abc import Mapping

    from orpheus.geometry.deck import DeckTransformation
    from orpheus.mesh import Axis1D
    from orpheus.transport.boundary_conditions import BoundaryCondition, BoundaryConditions, Binding, RealizerTable


__all__ = ["Mesh1DMethodMesh", "MethodMesh", "TransportMethod", "realize_boundary"]


#: The method's realized-value type: what its ``realize_boundary_role`` arm returns and its ``bc`` values carry. Unbounded: SN's is a
#: :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`, diffusion's a bare ``LinearOperator``, CP's a closure, MoC's a
#: track-linking rule and MC's a wrap marker (the shared body only stores the values). Covariant, so a consumer typed against
#: ``TransportMethod[object]`` accepts every method-mesh.
OpT_co = TypeVar("OpT_co", covariant=True)

#: Invariant twin of :data:`OpT_co` for the generic *function* below (a function signature has no variance; the Protocol keeps
#: the covariant declaration).
OpT = TypeVar("OpT")


@runtime_checkable
class TransportMethod(Protocol[OpT_co]):
    r"""Structural Protocol of the method-mesh layer.

    A method-mesh, a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` subclass carrying one transport method's
    behavior, satisfies this Protocol **without importing it** (see the module docstring for the conformance mechanism and the
    deliberate exclusions).

    The members are exactly what the two witnesses expose generically:

    * :attr:`axes`, :attr:`binding`, :attr:`boundary_realizers`, :meth:`realize_boundary_role` and :meth:`face_site`, consumed by
      the shared :func:`realize_boundary` body;
    * :attr:`bc`, the realized-operator table every method-mesh ends up with.

    One universal member is deliberately NOT declared: ``full_field_space`` (the bulk + trace composite carrier). Both witnesses
    implement it as a ``functools.cached_property``, which the current pyright does not accept against a Protocol property
    member, and it has no method-generic consumer yet (the shared body never touches it). The trigger stands for the NEXT
    genuinely method-generic consumer; do NOT convert the witnesses' caching to appease the checker today.
    """

    #: Canonical spatial representation (the :class:`MaterialMesh` data block): one :class:`~orpheus.mesh.axis.Axis1D` per
    #: dimension, each carrying its endpoints' ROLES (a deck, or an integer tag).
    axes: tuple["Axis1D", ...]

    #: The declared exterior conditions joined to this carrier's faces: ``decks`` (faces glued to a partner or a mirror
    #: plane), ``tags`` (faces facing the exterior), the ``declared`` conditions and the carrier's colouring symmetry.
    binding: "Binding"

    #: The method's realizer table: the ONE door that answers "is this declaration admitted on this face" and "realize it".
    #: It is sealed by the method's defining module, so anyone holding the table has executed every registration.
    boundary_realizers: ClassVar["RealizerTable"]

    @property
    def bc(self) -> "Mapping[str, OpT_co]":
        """Per-face REALIZED boundary operators, keyed by face name.

        Built at construction by :func:`realize_boundary` over the faces of the binding that carry a trace: the
        method's trace derives from the same ``face_labels`` inventory, so coverage of faces by operators ≡ coverage of faces
        by trace slots by construction.
        """
        ...

    def realize_boundary_role(self, role: "DeckTransformation | BoundaryCondition", face: str) -> OpT_co:
        """Realize one declaration on one face with a trace: the per-method arm.

        The shared :func:`realize_boundary` body owns the face loop and the choice of the declaration (a deck, or the
        condition the face's tag is bound to); this hook owns everything genuinely method-specific: building the method space for
        ``face`` and calling the method's own :attr:`boundary_realizers` ``realize``.
        """
        ...

    def face_site(self, face: str) -> "FaceSite":
        """Where ``face`` sits (its name, the carrier's coordinate system and dimension): the admission table's site key.

        Computable from the carrier alone, before any method space exists, so the type- and site-level admission of a face
        with no trace slot (a hollow body's inner face) needs nothing the method has not built yet.
        """
        ...


class MethodMesh(MaterialMesh):
    r"""The base of the method-meshes: a :class:`MaterialMesh` that owns a ``binding``, and the two duties that follow from owning one.

    A method-mesh joins the declared exterior conditions to its faces (:attr:`binding`, set by the subclass's constructor through
    :func:`~orpheus.transport.boundary_conditions.bind`), and owes the identity and the re-spelling machinery of the
    :class:`MaterialMesh` two things that were spelled twice, identically, in ``SNProblem`` and ``DiffusionMesh`` until the S4
    carve (CP's, MoC's and MC's meshes wrap a ``Mesh1D``, are not material meshes, and owe neither duty):

    * **the key.** :attr:`_contractibility_key` extends the material mesh's by the declared conditions in tag order
      (:func:`~orpheus.transport.mesh.material_mesh.conditions_key`), so two problems with the same mesh and different
      conditions are different phase spaces, as the laws in the axes' key were before the carve;
    * **the re-spelling.** A method that rebuilds the problem over the same phase space (``with_cross_sections``, the scattering
      order) re-passes :meth:`_respelling_keywords`, because a rebuilt carrier whose tags are unbound is refused.

    The subclass's own key adds its data after this one (``SNProblem``: the quadrature's content and the scheme type), so the
    order of the key's members is the subclass chain's order.
    """

    #: The declared conditions joined to this carrier's faces; set by the subclass's constructor.
    binding: "Binding"

    @cached_property
    def _contractibility_key(self) -> tuple:
        return (*MaterialMesh._contractibility_key.func(self), conditions_key(self.binding.declared))

    def _respelling_keywords(self) -> "dict[str, BoundaryConditions]":
        """The keywords a re-spelling constructor call carries so the rebuilt problem keeps this one's conditions."""
        return {"boundary_conditions": self.binding.declared}


class Mesh1DMethodMesh(Generic[OpT]):
    r"""The base of the method-meshes that wrap a :class:`~orpheus.mesh.Mesh1D`: the mesh as the method's space, and the bind-and-realize step.

    Collision probability, characteristics and Monte Carlo each hold a ``Mesh1D`` (``self.mesh``) and no trace space: the face they
    realize is a closure over the cell volumes (CP), a track-linking rule (MoC) or a wrap marker (MC), built from the mesh itself, so
    the mesh is the *method space* the realizer table receives. Until this base existed the three heads spelled, identically, the
    same three things: the step ``axes_from_legacy_mesh``, then :func:`~orpheus.transport.boundary_conditions.bind`, then
    :func:`realize_boundary`; the site of a face (:class:`~orpheus.transport.boundary_conditions.FaceSite`, whose dimension is the
    length of the axis tuple, 1 here); and the per-face hook that hands ``self`` to the table. A change to any of the three (a field
    on ``FaceSite``, a second carrier) edited three copies; it edits this class.

    A subclass declares its realizer table as the class variable ``boundary_realizers`` and calls :meth:`_bind_and_realize` once its
    own scope guards have passed (a guard that refuses a mesh the method misreads runs BEFORE the binding, so the refusal names the
    mesh and not a wall). The base derives from nothing: these meshes are not :class:`MaterialMesh` subclasses, so they owe neither
    the key nor the re-spelling of :class:`MethodMesh`.

    Parameters of the type: ``OpT``, the method's realized value (an arm of its realizer table), as for :func:`realize_boundary`.
    """

    mesh: Mesh1D
    axes: tuple["Axis1D", ...]
    binding: "Binding"
    bc: dict[str, OpT]
    boundary_realizers: ClassVar["RealizerTable"]

    def _bind_and_realize(self, mesh: Mesh1D, conditions: "BoundaryConditions | Mapping[int, BoundaryCondition]") -> None:
        """Join the declared conditions to the mesh's faces and realize the binding through the shared body (a refusal raises here)."""
        self.axes = axes_from_legacy_mesh(mesh)
        self.binding = bind(mesh, conditions)
        self.bc = realize_boundary(self)

    def face_site(self, face: str) -> FaceSite:
        """Where ``face`` sits: its name, the mesh's coordinate system and the axis count (the table's site key)."""
        return FaceSite(face, self.mesh.coord, len(self.axes))

    def realize_boundary_role(self, role: "DeckTransformation | BoundaryCondition", face: str) -> OpT:
        """Realize one declaration on ``face``: the arm of the shared :func:`realize_boundary` body, the mesh being the method space."""
        return self.boundary_realizers.realize(role, self.face_site(face), self)


def realize_boundary(method: "TransportMethod[OpT]") -> dict[str, OpT]:
    r"""Realize a method-mesh's binding into per-face operators: the ONE generic body behind every method-mesh's ``bc`` table.

    Replaces ``resolve_boundary_conditions`` and its tag-to-law parse. For each face of ``method.binding`` (inventory order):

    1. the declaration is the face's deck if it is glued, else the condition its tag is bound to
       (``binding.condition_of(face)``);
    2. a face with a trace slot (a member of ``face_labels(method.axes)``) is realized through the method's
       :meth:`~TransportMethod.realize_boundary_role`, whose table call admits first and raises
       :class:`~orpheus.transport.boundary_conditions.BoundaryRefused` iff refused;
    3. a face with no trace slot (the inner face of a hollow radial body: a tag, with no inner-surface trace) has nothing to
       build and meets the table through ``require`` at the site level only (``space=None``), which raises the same refusal.

    The face inventory of ``bc`` IS the trace inventory by construction: a face that carries a trace gets exactly one entry,
    and the faces without one (the curvilinear pole is not a face at all; a hollow body's inner face has a tag and no trace)
    have no entry rather than a placeholder.

    Returns
    -------
    dict[str, OpT]
        Face name to the method's realized operator, in the method's own operator type.

    Raises
    ------
    BoundaryRefused
        A declaration the method's table does not admit on its face (the refusal carries the admission and the face).
    """
    table, binding = method.boundary_realizers, method.binding
    traced = {label.face_name for label in face_labels(method.axes)}
    realized: dict[str, OpT] = {}
    for face in binding.faces:
        declared = binding.decks[face] if face in binding.decks else binding.condition_of(face)
        if face in traced:
            realized[face] = method.realize_boundary_role(declared, face)
        else:
            table.require(declared, method.face_site(face))
    return realized
