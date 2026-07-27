r"""TransportMethod — the Protocol over the method-mesh layer.

A **method-mesh** is a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
(the method-agnostic mesh + materials DATA) augmented with one
transport method's BEHAVIOR. Two exist:
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (quadrature + sweep
machinery + angular trace) and
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` (scalar trace
+ realized boundary laws). :class:`TransportMethod` is the structural
Protocol naming what every method-mesh exposes *as a method*, and
:func:`resolve_boundary_conditions` is the ONE generic body that turns
the mesh axes' per-face :class:`~orpheus.geometry.mesh.BC` declarations
into realized boundary operators through it.

Genesis — the two recorded witnesses (#290 P7b)
===============================================

The Protocol was deliberately **deferred until a second method-mesh
existed** (the defer-until-≥2 rule). Its two recorded witnesses:

1. the **homogenization method-layer** —
   :mod:`orpheus.transport.mesh.material_mesh`'s docstring promised
   "a method-specific mesh IS a MaterialMesh that adds behavior,
   conforming structurally to a future ``TransportMethod`` Protocol";
2. the **boundary-realizer seam** — ``SNMesh._resolve_bcs`` and
   ``DiffusionMesh._resolve_bcs`` were twin loops (shape-identical
   after #290 P7a moved diffusion's resolution off the solver onto the
   phase space), each doing the same tag → law → realize walk with only
   the leaf realization differing.

#290 P7a minted ``DiffusionMesh`` (the missing second witness); this
module is P7b — the mint, plus the twin-loop collapse.

What the Protocol deliberately does NOT contain
===============================================

* **Promotion classmethods.** ``SNMesh.from_material_mesh(mm,
  quadrature, scheme=...)`` vs ``DiffusionMesh.from_material_mesh(mm)``
  — the signatures differ because the methods genuinely need different
  injections (SN has a quadrature; diffusion reads everything from the
  carrier). The asymmetry is honest; forcing a common signature would
  fabricate parameters. The Protocol is **instance surface only**.
* **Method-space construction.** ``SNMethodSpace.for_face(mesh,
  quadrature, face, trace)`` vs ``DiffusionMethodSpace.for_face(mesh,
  face)`` — non-unifiable for the same reason. The per-method
  :meth:`~TransportMethod.realize_boundary_law` hook subsumes it: the
  method space is an implementation detail *inside* each arm.
* **The trace spaces.** ``SNMesh.angular_trace`` (an
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`)
  and ``DiffusionMesh.scalar_trace`` (a
  :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`)
  are locus-family SIBLINGS, not one concept — a consumer that needs a
  trace needs the method's own trace type.
* **A dispatch name.** There is no ``method_name`` member: the one
  generic consumer that needs to name the method (the unsupported-tag
  refusal) uses ``type(method).__name__`` — self-truthing, and never
  usable for string dispatch (the anti-pattern the #290 P7b registry
  dissolution removed).

Conformance is structural — nobody imports this to conform
==========================================================

``orpheus.sn`` and ``orpheus.diffusion`` import
:func:`resolve_boundary_conditions` (the shared body), **never the
Protocol class**. Static conformance is checked exactly where it
matters: at each mesh's ``resolve_boundary_conditions(self)`` call,
pyright verifies the mesh structurally satisfies
``TransportMethod[OpT]`` and infers the method's own realized-operator
type ``OpT`` (SN: the ``_BoundBoundaryOperator`` kind-tagged shim;
diffusion: a bare :class:`~orpheus.numerics.operator.LinearOperator`).
Future generic consumers (the DSA driver #2 is the named one — it holds
an ``SNMesh`` and a promoted ``DiffusionMesh`` over the same carrier)
type against the Protocol under ``TYPE_CHECKING``.

The registry that this mint dissolved
=====================================

The Wave-5 ``BoundaryRealizerRegistry`` (string-keyed
``method_name → BoundaryRealizer``, with ``NotImplementedError`` stub
realizers auto-registered by ``import orpheus.{moc,mc,cp}``) was
retired at this carve. The registry existed to answer "given a method
NAME, find its realizer" — but no production consumer ever asked that
question: **you hold the method-mesh, therefore you hold its realizer**
(each :meth:`realize_boundary_law` arm calls its own realizer
directly). The string indirection carried real hazards for zero payoff
— registration-timing misses invisible to in-suite tests
(process-global state masked the empty-registry condition) and
import-side-effect coupling — and dissolving it deletes the hazard
class rather than gating it. Adding a transport method now means:
mint the method-mesh (conforms here structurally), its realizer, and
its ``BOUNDARY_OPERATOR_REGISTRY`` of admitted laws — no central
registration step.

Layer (``tests/test_layer_imports.py``): L2 ``transport``. Imports
``geometry`` (the BC tag + the law types) and ``numerics`` (the
operator base), like its sibling :mod:`~orpheus.transport.mesh.axis`.
The method packages (L3) import THIS module; never the reverse.

References
----------

* ``.claude/plans/diffusion_integration_290.md`` §P7b (the binding
  carve spec; supersedes the original §P7 and the banked
  ``realize_recursively_move_spec.md`` registry-routing design).
* :doc:`/theory/foundations/boundary_conditions` — the three-layer decomposition
  (trace structure / physical law / method realisation) this module's
  resolve body walks.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar, Protocol, TypeVar, runtime_checkable

from orpheus.geometry import BC
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryTraceLaw,
    ReflectiveBoundary,
)
from orpheus.numerics.operator import LinearOperator
from orpheus.transport.mesh.axis import AXIS_NAMES, face_labels

if TYPE_CHECKING:
    from collections.abc import Mapping

    from orpheus.transport.mesh.axis import Axis1D, FaceLabel


__all__ = ["TransportMethod", "resolve_boundary_conditions"]


#: The method's realized-operator type: what its ``realize_boundary_law``
#: arm returns and its ``bc`` values carry. SN narrows it to the
#: kind-tagged ``_BoundBoundaryOperator`` shim; diffusion uses the bare
#: :class:`LinearOperator`. Covariant — a consumer typed against
#: ``TransportMethod[LinearOperator]`` accepts every method-mesh.
OpT_co = TypeVar("OpT_co", bound=LinearOperator, covariant=True)

#: Invariant twin of :data:`OpT_co` for the generic *function* below
#: (a function signature has no variance; the Protocol keeps the
#: covariant declaration).
OpT = TypeVar("OpT", bound=LinearOperator)


@runtime_checkable
class TransportMethod(Protocol[OpT_co]):
    r"""Structural Protocol of the method-mesh layer.

    A method-mesh — a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
    subclass carrying one transport method's behavior — satisfies this
    Protocol **without importing it** (see the module docstring for the
    conformance mechanism and the deliberate exclusions).

    The members are exactly what the two witnesses expose generically:

    * :attr:`axes` + :attr:`BOUNDARY_OPERATOR_REGISTRY` +
      :meth:`realize_boundary_law` — consumed by the shared
      :func:`resolve_boundary_conditions` body;
    * :attr:`bc` — the realized-law table every method-mesh ends up
      with.

    One universal member is deliberately NOT declared:
    ``full_field_space`` (the bulk ⊕ trace composite carrier). Both
    witnesses implement it as a ``functools.cached_property``, which
    the current pyright does not accept against a Protocol property
    member (``cached_property[T]`` is compared as the descriptor, not
    its instance-access type) — and it has no method-generic consumer
    yet (the shared resolve body never touches it). The anticipated
    first consumer ("the DSA driver, #2") did NOT materialize: the R4
    ruling (2026-07-26) wired consistent DSA through an SN-side
    edge-centered low-order system consuming ``SNMesh`` directly —
    no ``DiffusionMesh`` enters the accelerated loop, so no generic
    ``full_field_space`` read exists. The trigger stands for the NEXT
    genuinely method-generic consumer; do NOT convert the witnesses'
    caching to appease the checker today.
    """

    #: Canonical spatial representation (the :class:`MaterialMesh`
    #: data block): one :class:`~orpheus.transport.mesh.axis.Axis1D`
    #: per dimension, each carrying its endpoints' ``BC`` declarations.
    axes: tuple["Axis1D", ...]

    #: The method's law admission table: ``BC.kind`` tag → the
    #: :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` subclass it
    #: declares. A tag absent here is REFUSED at construction with the
    #: supported list (each method admits exactly the laws its realizer
    #: can honestly realize — e.g. ``zero_flux`` is diffusion-only,
    #: and ``white`` is deliberately absent from diffusion where it
    #: coincides with reflective at P1).
    BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]]

    @property
    def bc(self) -> "Mapping[str, OpT_co]":
        """Per-face REALIZED boundary laws, keyed by face name.

        Built at construction by :func:`resolve_boundary_conditions`
        over the same ``face_labels`` inventory the method's trace
        derives from — law coverage ≡ face coverage by construction.
        """
        ...

    def realize_boundary_law(
        self,
        law: BoundaryTraceLaw,
        face: str,
    ) -> OpT_co:
        """Realize one typed law on one face — the per-method arm.

        The shared :func:`resolve_boundary_conditions` body owns the
        face loop, the reflective default, and the tag → law parse;
        this hook owns everything genuinely method-specific: building
        the method space for ``face`` and dispatching the method's own
        realizer.
        """
        ...


def resolve_boundary_conditions(
    method: "TransportMethod[OpT]",
) -> dict[str, OpT]:
    r"""Resolve a method-mesh's per-axis BC declarations into realized operators.

    The ONE generic body behind every method-mesh's ``bc`` table
    (#290 P7b — it replaced the twin ``SNMesh._resolve_bcs`` /
    ``DiffusionMesh._resolve_bcs`` loops; a third spelling on the
    diffusion solver died at P7a). Per face label of ``method.axes``:

    1. read the axis's declaration for that endpoint —
       ``axes[label.axis_index].bc[label.endpoint]``; ``None`` defaults
       to ``BC("reflective")`` (the infinite-lattice / eigenvalue
       convention, uniform across methods);
    2. parse the tag into its typed law via the method's
       :attr:`~TransportMethod.BOUNDARY_OPERATOR_REGISTRY` (see
       :func:`_law_from_tag` — unsupported tags and a parameter-less
       albedo refuse loudly HERE, at phase-space construction);
    3. hand the law to the method's
       :meth:`~TransportMethod.realize_boundary_law` arm.

    The face inventory IS the BC inventory by construction: a face
    that exists gets exactly one entry; a face that doesn't (the
    curvilinear pole r=0 — a regularity condition, not a BC face) is
    structurally absent. Both method traces derive from the SAME
    ``face_labels`` inventory, so law coverage ≡ trace-face coverage
    with no separate validation to drift.

    Returns
    -------
    dict[str, OpT]
        Face name → the method's realized operator, in the method's
        own operator type (``OpT`` infers from the mesh's
        ``realize_boundary_law`` return — the kind-tagged shim for SN,
        a bare :class:`LinearOperator` for diffusion).
    """
    default = BC("reflective")
    resolved: dict[str, OpT] = {}
    for label in face_labels(method.axes):
        tag = method.axes[label.axis_index].bc[label.endpoint] or default
        law = _law_from_tag(method, tag, label)
        resolved[label.face_name] = method.realize_boundary_law(
            law, label.face_name,
        )
    return resolved


def _law_from_tag(
    method: "TransportMethod[Any]",
    tag: BC,
    label: "FaceLabel",
) -> BoundaryTraceLaw:
    r"""Parse one ``BC`` tag into the typed law it declares.

    Method-generic given the method's admission table: the tag names
    the law class; the label supplies the face geometry. Two law
    families need construction context, resolved here in ONE place:

    * a **reflective** law reflects across the face's own axis —
      ``AXIS_NAMES[label.axis_index]`` — so the reflection partner is
      correct at any dimension by construction (a hand-listed
      face → axis map would silently build the wrong permutation for
      a z-face);
    * an **albedo** law carries its response as the tag parameter
      ``params["albedo"]``; a parameter-less ``BC("albedo")`` refuses
      with the face named.

    Every other admitted law is parameter-free (``law_cls()``). A law
    class whose tag needs new parameters extends this parse — one
    body, every method.
    """
    law_cls = method.BOUNDARY_OPERATOR_REGISTRY.get(tag.kind)
    if law_cls is None:
        supported = ", ".join(
            f"'{k}'" for k in sorted(method.BOUNDARY_OPERATOR_REGISTRY)
        )
        raise ValueError(
            f"{type(method).__name__} does not support boundary condition "
            f"'{tag.kind}' on face '{label.face_name}'. "
            f"Supported: {supported}."
        )
    if law_cls is ReflectiveBoundary:
        return ReflectiveBoundary(
            axis=AXIS_NAMES[label.axis_index], albedo=1.0,
        )
    if law_cls is AlbedoBoundary:
        try:
            return AlbedoBoundary(albedo=float(tag.params["albedo"]))
        except KeyError as exc:
            raise ValueError(
                f"BC('albedo') on face '{label.face_name}' requires an "
                f"'albedo' parameter; got params={tag.params!r}."
            ) from exc
    return law_cls()
