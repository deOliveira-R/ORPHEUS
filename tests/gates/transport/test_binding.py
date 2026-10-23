r"""The binding: the declared conditions joined to the faces of a carrier (wave 2, slice S2a, gate 14).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("``declared.py`` and ``binding.py``", decision 51 and the
carrier ruling of 2026-10-23), law L-C, 6.2 S2a gate 14.

``bind(carrier, conditions)`` is the ONE shared body every method head calls. The carrier owns the faces and their
colouring symmetry (a ``BoundaryCarrier``: ``face_roles`` as a ``Mapping[str, FaceRole]`` and ``symmetry``), so a
binding checked against another carrier's group is unspellable. It refuses, in this order: a tag a face carries with
no condition (``UnboundTagError``), a condition whose tag no face carries (``UnusedConditionError``), a condition
that is not invariant under the carrier's group (``NonInvariantConditionError``).

The claims (``vv-principles``):

* ``test_bind_takes_the_carrier`` (decision 51): STRUCTURE. The signature is exactly ``(carrier, conditions)``;
  the Protocol is satisfied by a test-local carrier, by a real ``Mesh1D`` and by a real ``Mesh2D``; a non-carrier is
  refused; the binding's symmetry IS the carrier's.
* ``test_bind_exact_cover`` and ``test_binding_invariants``: THEOREM on the cover (the tags faces carry equal the
  tags declared), its refusal order, and the type's own invariant (a ``Binding`` that joins a face to a condition of
  another tag cannot be built: the join is DERIVED from ``(tags, declared)``, never stored).
* ``test_bind_invariance`` (law L-C): THEOREM. For every non-identity ``h`` of the carrier's group and every tag whose
  faces ``h`` moves, the declared condition is invariant under ``h``. A lab-frame source on a tag a mirror swaps is
  refused; the controls bind. ``face_local`` implies invariant ONLY in one dimension (a 2-D mirror fixing a face
  reverses its tangent frame), which is why a ``UncomputedSymmetries`` carrier refuses a non-face-local condition on a
  multi-face tag, naming the missing machinery (``SCOPE-BOUNDARY[guard]``, the 2-D colouring group).

Activation: the symmetric slab ``(0, 1, 2)`` with tag 0 on both faces activates the mirror that the asymmetric slab
and the two-tag slab null; the lab-frame source activates the only condition kind that owes an invariance check (a
face-local source nulls it).

The authored first red (law L-C): the symmetric slab, tag 0 on both faces, ``framed(LabFrameInflow(f))`` raises
``NonInvariantConditionError``. At the base tip the equivalent typed ``PrescribedInflow`` realizes on both faces
silently (``[M]`` the plan's law table). Base-tip red here: ``ModuleNotFoundError``.

Mutations (each on a COPY of ``binding.py`` / ``condition.py``): replace ``invariant_under`` by ``True``; replace the
stabiliser by ``{e}`` (iterate the identity only); make ``bind`` take ``symmetry`` as a second argument; drop the
exact-cover check; check the cover after the invariance; read ``sources`` through ``isinstance(c, WallCondition)``.
"""

from __future__ import annotations

import inspect

import pytest

from orpheus.geometry.colouring import ColouringSymmetry, UncomputedSymmetries
from orpheus.geometry.deck import MirrorDeck
from orpheus.values.content import ContentlessError, FrozenMapping, content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import (
    Carrier,
    OpaqueSampler,
    bc,
    local,
    mesh2d,
    node,
    slab_mesh,
)

pytestmark = pytest.mark.foundation

_HERE = "test_binding.py"
_MAP = node("test_boundary_conditions_map.py", "test_conditions_map")
_ALGEBRA = node("test_wall_condition_algebra.py", "test_the_algebra_and_the_cone")
_FRAMES = node("test_wall_condition_algebra.py", "test_inflow_frames_and_terms")
_COLOURING = "tests/gates/geometry/test_colouring_symmetry.py::test_colouring_symmetry_reads_runs_and_roles"
_CARRIER = node(_HERE, "test_bind_takes_the_carrier")
_COVER = node(_HERE, "test_bind_exact_cover")


def _lab():
    return bc.WallCondition.framed(bc.LabFrameInflow(OpaqueSampler()))


def _symmetric():
    """The slab ``(0, 1, 2)``, material 0 throughout, tag 0 on both faces: the group is ``{e, mirror}``."""
    return slab_mesh({"xmin": 0, "xmax": 0})


# ── gate 14, part 1: the carrier ─────────────────────────────────────────────


@pytest.mark.rests_on(_MAP, _COLOURING)
def test_bind_takes_the_carrier() -> None:
    """``bind`` has exactly the parameters ``(carrier, conditions)`` (an ``inspect.signature`` read); the
    three-argument spelling ``bind(roles, conditions, symmetry)`` of DRAFT 2 is a ``TypeError``; a test-local
    carrier satisfying ``BoundaryCarrier``, a real ``Mesh1D`` and a real ``Mesh2D`` are carriers
    (``runtime_checkable``), a plain mapping of roles is not and is refused; the ``Binding`` holds the carrier's
    symmetry (``is``) and its roles; the same roles under another carrier's group give another answer (the mismatch
    that the two-value spelling admitted is unspellable).

    The ``MaterialMesh`` arm of the Protocol is slice S2b's gate 21.

    Reddens on: a copy of ``bind`` taking ``symmetry`` as a second argument (the signature row and the ``TypeError``
    row); one that re-reads the group from the roles instead of the carrier.
    """
    params = list(inspect.signature(bc.bind).parameters)
    require(params == ["carrier", "conditions"], f"bind takes {params}")
    sym_mesh = _symmetric()
    with pytest.raises(TypeError):
        bc.bind(dict(sym_mesh.face_roles), {0: bc.WallCondition.vacuum()}, sym_mesh.symmetry)  # type: ignore[call-arg]  # DRAFT 2's spelling
    with pytest.raises((TypeError, AttributeError)):
        bc.bind(dict(sym_mesh.face_roles), {0: bc.WallCondition.vacuum()})  # type: ignore[arg-type]  # roles are not a carrier
    require(not isinstance(dict(sym_mesh.face_roles), bc.BoundaryCarrier), "a plain mapping of roles is not a carrier")

    own = Carrier({"xmin": 0, "xmax": 0}, sym_mesh.symmetry)
    mesh_2d = mesh2d({"xmin": 0, "xmax": 0, "ymin": 1, "ymax": 2})
    for carrier in (own, sym_mesh, mesh_2d):
        require(isinstance(carrier, bc.BoundaryCarrier), f"{type(carrier).__name__} does not satisfy the BoundaryCarrier protocol")
    require(isinstance(mesh_2d.symmetry, UncomputedSymmetries) and isinstance(sym_mesh.symmetry, ColouringSymmetry), "control: the two meshes carry the two kinds of group")

    vac = bc.WallCondition.vacuum()
    binding = bc.bind(own, {0: vac})
    require(binding.symmetry is own.symmetry, "the binding holds another group than the carrier's")
    require(binding.faces == ("xmin", "xmax") and set(binding.tags) == {"xmin", "xmax"}, f"the binding's faces are {binding.faces}")

    trivial = bc.bind(Carrier({"xmin": 0, "xmax": 0}, slab_mesh({"xmin": 0, "xmax": 0}, materials=(0, 1)).symmetry), {0: _lab()})
    require(isinstance(trivial, bc.Binding) and len(trivial.symmetry.elements) == 1, "the asymmetric group does not bind the lab-frame source silently")
    with pytest.raises(bc.NonInvariantConditionError):
        bc.bind(own, {0: _lab()})  # the SAME roles, the symmetric carrier's group: refused


# ── gate 14, part 2: the exact cover ─────────────────────────────────────────


@pytest.mark.rests_on(_CARRIER)
def test_bind_exact_cover() -> None:
    """``UnboundTagError`` names the face and the tag; ``UnusedConditionError`` names the tag; a ``bool`` tag is
    refused; the refusals come in the order unbound, unused, non-invariant; ``bind`` coerces a bare mapping through
    ``BoundaryConditions.of`` once and accepts a ``BoundaryConditions``; a deck-only carrier binds the empty map.

    Reddens on: dropping the cover check (a typo'd tag then binds); a cover check after the invariance check (the
    last two rows); a ``bind`` that does not coerce.
    """
    vac = bc.WallCondition.vacuum()
    two_tags = slab_mesh({"xmin": 0, "xmax": 1})

    with pytest.raises(bc.UnboundTagError) as unbound:
        bc.bind(two_tags, {0: vac})
    require(unbound.value.tag == 1 and unbound.value.face == "xmax", f"UnboundTagError names tag {unbound.value.tag}, face {unbound.value.face!r}")
    require("xmax" in str(unbound.value) and "1" in str(unbound.value), f"the message does not name the face and the tag: {unbound.value}")
    require(isinstance(unbound.value, bc.BoundaryError), "UnboundTagError is a BoundaryError")

    with pytest.raises(bc.UnusedConditionError) as unused:
        bc.bind(_symmetric(), {0: vac, 5: vac})
    require(unused.value.tag == 5 and "5" in str(unused.value), f"UnusedConditionError names tag {unused.value.tag}")

    with pytest.raises(bc.UnboundTagError):
        bc.bind(two_tags, {0: vac, 7: vac})  # unbound (tag 1) is read before unused (tag 7)
    with pytest.raises(bc.UnusedConditionError):
        bc.bind(_symmetric(), {0: _lab(), 9: vac})  # unused (tag 9) is read before the invariance of tag 0
    with pytest.raises(TypeError):
        bc.bind(_symmetric(), {True: vac})

    as_mapping = bc.bind(_symmetric(), {0: vac})
    as_object = bc.bind(_symmetric(), bc.BoundaryConditions({0: vac}))
    require(isinstance(as_mapping.declared, bc.BoundaryConditions), "the binding holds the declaration object")
    require(as_mapping == as_object and content_digest(as_mapping) == content_digest(as_object), "a mapping and a declaration bind to one value")

    half_cell = slab_mesh({"xmin": MirrorDeck(), "xmax": MirrorDeck()})
    empty = bc.bind(half_cell, bc.BoundaryConditions.empty())
    require(empty.tags == FrozenMapping({}) and set(empty.decks) == {"xmin", "xmax"}, "a deck-only carrier binds the empty declaration")
    with pytest.raises(bc.UnusedConditionError):
        bc.bind(half_cell, {0: vac})


# ── gate 14, part 3: the invariance (law L-C) ────────────────────────────────


@pytest.mark.rests_on(_COVER, _FRAMES, _COLOURING)
def test_bind_invariance() -> None:
    """The authored first red and its controls, with the quantifier of the law.

    Refused: the symmetric slab ``(0, 1, 2)`` with materials ``(0, 0)`` and tag 0 on both faces, bound to
    ``framed(LabFrameInflow(f))``: ``NonInvariantConditionError`` naming the tag (0), the two faces and the
    non-identity isometry. Bind: the same with a uniform source, with a face-frame source, with a face-local
    response; an asymmetric slab (materials ``(0, 1)``) with the lab-frame source; the slab with two tags; the slab
    whose tag 0 has one face (a deck on the other); and ``UncomputedSymmetries`` (a real ``Mesh2D``) with a face-local
    condition on a two-face tag, and with a lab-frame source on a ONE-face tag.

    Refused with the missing-machinery message: ``UncomputedSymmetries`` with a two-face tag and a lab-frame source
    (``SCOPE-BOUNDARY[guard]``: the 2-D colouring group is not computed, so a non-face-local condition on a
    multi-face tag cannot be certified invariant).

    Reddens on: ``invariant_under`` replaced by ``True`` (the refused rows bind); the stabiliser replaced by the
    identity alone (the refused row binds); the ``UncomputedSymmetries`` arm deleted (the last refused row binds).
    """
    lab = _lab()
    with pytest.raises(bc.NonInvariantConditionError) as refused:
        bc.bind(_symmetric(), {0: lab})
    err = refused.value
    require(err.tag == 0 and tuple(err.faces) == ("xmin", "xmax"), f"the refusal names tag {err.tag}, faces {err.faces}")
    require(err.isometry is not None and not err.isometry.is_identity, "the refusal names the non-identity isometry that moves the faces")
    require(isinstance(err, bc.BoundaryError) and (err.declaration is lab or err.declaration is type(lab)), "the refusal carries the declaration (the object or its type)")

    sampler = OpaqueSampler()
    controls = {
        "uniform source, symmetric slab": (_symmetric(), {0: bc.WallCondition.uniform(1.0)}),
        "face-frame source, symmetric slab": (_symmetric(), {0: bc.WallCondition.framed(bc.FaceFrameInflow(sampler))}),
        "face-local response, symmetric slab": (_symmetric(), {0: bc.WallCondition.white(0.5)}),
        "lab-frame source, asymmetric slab": (slab_mesh({"xmin": 0, "xmax": 0}, materials=(0, 1)), {0: bc.WallCondition.framed(bc.LabFrameInflow(sampler))}),
        "lab-frame source, two tags": (slab_mesh({"xmin": 0, "xmax": 1}), {0: bc.WallCondition.framed(bc.LabFrameInflow(sampler)), 1: bc.WallCondition.vacuum()}),
        "lab-frame source, one-face tag": (slab_mesh({"xmin": MirrorDeck(), "xmax": 0}), {0: bc.WallCondition.framed(bc.LabFrameInflow(sampler))}),
        "face-local response, uncomputed group": (mesh2d({"xmin": 0, "xmax": 0, "ymin": 1, "ymax": 2}),
                                                  {0: bc.WallCondition.white(0.5), 1: bc.WallCondition.vacuum(), 2: bc.WallCondition.vacuum()}),
        "lab-frame source, uncomputed group, one-face tag": (mesh2d({"xmin": 0, "xmax": 0, "ymin": 1, "ymax": 2}),
                                                            {0: bc.WallCondition.vacuum(), 1: bc.WallCondition.framed(bc.LabFrameInflow(sampler)), 2: bc.WallCondition.vacuum()}),
    }
    for label, (carrier, conditions) in controls.items():
        bound = bc.bind(carrier, conditions)
        require(isinstance(bound, bc.Binding), f"{label}: did not bind")
    require(len(controls) >= 8, "the control population shrank")

    with pytest.raises(bc.NonInvariantConditionError) as machinery:
        bc.bind(mesh2d({"xmin": 0, "xmax": 0, "ymin": 1, "ymax": 2}),
                {0: lab, 1: bc.WallCondition.vacuum(), 2: bc.WallCondition.vacuum()})
    require(machinery.value.tag == 0 and "colouring" in str(machinery.value).lower(),
            f"the refusal does not name the colouring group as the missing machinery: {machinery.value}")


# ── the type's own invariant ─────────────────────────────────────────────────


@pytest.mark.rests_on(_COVER)
def test_binding_invariants() -> None:
    """``decks`` and ``tags`` partition the faces (disjoint, together the inventory in order); ``set(tags.values())``
    equals ``declared.tags`` EXACTLY, so a ``Binding`` joining a face to a condition of another tag cannot be built;
    ``condition_of``, ``faces_of`` and ``sources`` are derived from ``(tags, declared)``; ``sources`` reads the
    ``source`` accessor of EVERY condition type (a method-only condition with a nonzero source appears); a
    binding holding an opaque sampler has no digest (``ContentlessError``).

    Reddens on: a ``Binding`` that stores a face-to-condition map beside ``tags``; a ``sources`` that tests
    ``isinstance(c, WallCondition)``; a missing set equality (the direct-construction rows build).
    """
    vac, white = bc.WallCondition.vacuum(), bc.WallCondition.white(0.5)
    mesh = slab_mesh({"xmin": MirrorDeck(), "xmax": 0})
    binding = bc.bind(mesh, {0: white})
    require(set(binding.decks) == {"xmin"} and set(binding.tags) == {"xmax"} and set(binding.decks).isdisjoint(binding.tags), "decks and tags do not partition the faces")
    require(binding.decks["xmin"] == MirrorDeck() and binding.tags["xmax"] == 0, "the partition holds the wrong values")
    require(binding.faces == ("xmin", "xmax"), f"faces are {binding.faces}")
    require(binding.condition_of("xmax") == white, "condition_of(xmax) is the declared condition of its tag")
    require(binding.faces_of(0) == ("xmax",) and binding.faces_of(9) == (), "faces_of")
    require(len(binding.sources) == 0, "a source-free binding has sources")

    two = bc.bind(slab_mesh({"xmin": 0, "xmax": 0}), {0: bc.WallCondition.uniform(2.0)})
    require(set(two.sources) == {"xmin", "xmax"} and all(t.uniform == 2.0 for t in two.sources.values()), "sources lost a face or its weight")
    require(two.faces_of(0) == ("xmin", "xmax"), "faces_of keeps the inventory order")
    one = bc.bind(slab_mesh({"xmin": 0, "xmax": 1}), {0: vac, 1: bc.WallCondition.uniform(2.0)})
    require(set(one.sources) == {"xmax"}, f"sources holds {set(one.sources)}: only the faces whose condition has a nonzero source")
    sourced = local().SourcedProbe()
    through_accessor = bc.bind(slab_mesh({"xmin": 0, "xmax": 1}), {0: vac, 1: sourced})
    require(set(through_accessor.sources) == {"xmax"} and through_accessor.sources["xmax"].uniform == 3.0,
            "a method-only condition's source is not read: sources must go through BoundaryCondition.source")

    sym = _symmetric().symmetry
    empty_decks = FrozenMapping({})
    for tags, declared in (
        (FrozenMapping({"xmax": 1}), bc.BoundaryConditions({0: vac})),  # a face joined to a condition of another tag
        (FrozenMapping({"xmax": 0}), bc.BoundaryConditions({0: vac, 1: vac})),  # a declared tag no face carries
    ):
        with pytest.raises((ValueError, TypeError)):
            bc.Binding(decks=empty_decks, tags=tags, declared=declared, symmetry=sym)
    with pytest.raises((ValueError, TypeError)):
        bc.Binding(decks=FrozenMapping({"xmin": MirrorDeck()}), tags=FrozenMapping({"xmin": 0}), declared=bc.BoundaryConditions({0: vac}), symmetry=sym)

    opaque = bc.bind(_symmetric(), {0: bc.WallCondition.framed(bc.FaceFrameInflow(OpaqueSampler()))})
    with pytest.raises(ContentlessError):
        content_digest(opaque)
    plain = bc.bind(_symmetric(), {0: white})
    require(content_digest(plain) == content_digest(bc.bind(_symmetric(), {0: bc.WallCondition.white(0.5)})), "equal bindings have two digests")


# ── the type owns the laws the door states (S2a review) ──────────────────────


@pytest.mark.rests_on(_COVER)
def test_a_direct_binding_holds_the_laws_of_the_door() -> None:
    """``Binding(...)`` built directly refuses what ``bind`` refuses: a condition that is not invariant under the
    group (``NonInvariantConditionError``), a group acting on other faces than the binding's (``ValueError``: another
    carrier's group) and a face outside the face inventory (``ValueError``, no silent drop from ``faces`` and no bare
    ``KeyError`` from ``bind``). Controls: the same faces with a face-local condition build; the same faces with the
    group of the same faces build.

    Reddens on: the invariance check living only in ``bind`` (the first row builds); the isometry-domain check
    deleted (the second row builds); the face-name parse deleted (the third row builds and ``bind`` raises ``KeyError``).
    """
    sym = _symmetric().symmetry
    both = FrozenMapping({"xmin": 0, "xmax": 0})
    nodeck = FrozenMapping({})
    white = bc.BoundaryConditions({0: bc.WallCondition.white(0.5)})
    require(isinstance(bc.Binding(nodeck, both, white, sym), bc.Binding), "control: a face-local condition on the group's faces does not build")
    with pytest.raises(bc.NonInvariantConditionError):
        bc.Binding(nodeck, both, bc.BoundaryConditions({0: _lab()}), sym)
    with pytest.raises(ValueError, match="another carrier"):
        bc.Binding(nodeck, FrozenMapping({"xmax": 0}), white, sym)
    with pytest.raises(ValueError, match="not faces"):
        bc.Binding(nodeck, FrozenMapping({"outer": 0}), white, UncomputedSymmetries("Mesh2D"))
    with pytest.raises(ValueError, match="not faces"):
        bc.bind(Carrier({"outer": 0}, UncomputedSymmetries("Mesh2D")), {0: bc.WallCondition.white(0.5)})
