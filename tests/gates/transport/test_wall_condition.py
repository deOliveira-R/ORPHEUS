r"""The one shared exterior condition: its parse, its facades, its total chart (wave 2, slice S2a, gates 1, 2, 3, 11).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("``condition.py``"), section 3 (laws L-B, L-G), 6.2 S2a.

The claims, each with its layer and its kind (``vv-principles``):

* gate 1, ``test_wall_condition_parse``: THEOREM (law L-B). A wall is an affine pair ``(R, q)`` whose linear part
  lies in the sub-Markov cone ``K = {s, d >= 0, s + d <= 1}`` and whose source lies in the positive cone, so a
  value outside ``K`` is unconstructible, however it is spelled. The albedo twin is gone: ``specular(0.0)``,
  ``white(0.0)`` and ``vacuum()`` are ONE value (one digest).
* gate 2, ``test_orientation_is_unspellable`` and ``test_no_welded_surface_survives``: STRUCTURE (an introspection
  over the closed population, with a positive control). A wall is face-local: no orientation field exists on any
  condition, factor or chain, and none of the welds of the retired ``BoundaryTraceLaw`` (a string ``kind``, a
  ``response_kernel``, a deck slot ``geometry_map``, ``source.evaluate(space)``) survives on the new base.
* gate 3, ``test_white_is_a_facade``: THEOREM (law L-G, one definition per quantity). ``white(a)`` and
  ``from_amplitudes(0.0, a)`` are one value of one type, and ``a`` is the first positional argument.
* gate 11, ``test_from_amplitudes_is_total_on_the_cone``: THEOREM. ``from_amplitudes`` is the section of
  ``cone_point``: ``cone_point . from_amplitudes = id`` on ``K``, and nothing outside ``K`` is accepted. The edge
  takes two primitives, never a ``WallAmplitudes`` (the reference's value): references are insulated from production.

The module reads the package through ``bc`` at run time (see ``_boundary_conditions_helpers``), so at the slice's
base tip every row's red is its own ``ModuleNotFoundError`` line, and once the package exists each row's own
assertion.

First reds, at the base tip ``4ec95bc6``: every row fails with ``ModuleNotFoundError: No module named
'orpheus.transport.boundary_conditions'`` (rows 2 and 12 of the package file read their own assertions, after their
controls). The inputs that redden the rows once the module exists, ``[M]`` read at ``4ec95bc6`` from the retired
vocabulary (``orpheus.geometry.boundary``, which still imports there): ``WhiteBoundary(albedo=5.0)`` CONSTRUCTS and
keeps ``albedo == 5.0`` (and a diffusion realization of it is ``A = 5``, ``k = 2.972`` on a slab,
``scratch/layer_dag/wave2_spec/_probe_admission.out``), so the old tree admits what the new parse refuses;
``AlbedoBoundary(0.0, SpecularReturn(axis="x")) == VacuumInflow()`` is ``False`` (the twin), where
``specular(0.0) == vacuum()`` here; ``SpecularReturn(axis="y")`` constructs and keeps ``axis == "y"`` (gate 2). The
gates carry no import of that vocabulary: WF2.4 deletes it, and a row importing it would die with it.

Mutations (each applied to a COPY of the module, never by ``git checkout``; the reviewers' battery):

* gate 1: drop the cone parse from ``WallCondition.__post_init__`` (``specular(1.2)`` then constructs); drop the
  ``-0.0`` canonicalisation or the zero-weight drop of ``ResponseNetwork`` (``specular(0.0) != vacuum()``); accept a
  ``bool`` albedo; admit a negative inflow weight.
* gate 2: add a field ``axis: str = "x"`` to ``MirrorPermutation`` (introspection and ``MirrorPermutation(axis=)``
  both stop failing); restore ``kind`` on ``BoundaryCondition``.
* gate 3: make ``white`` return a ``WhiteCondition`` subclass; swap the argument order of ``white``.
* gate 11: make ``from_amplitudes`` drop the diffuse weight when ``specular > 0``; take a ``WallAmplitudes``.
"""

from __future__ import annotations

import inspect
import math
from typing import get_args

import pytest

from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import (
    ProfileSampler,
    bc,
    cone_grid,
    local,
    node,
    orientation_fields,
    walk_subclasses,
    wall,
)

pytestmark = pytest.mark.foundation

_HERE = "test_wall_condition.py"
_NETWORK = "test_response_network.py"
_PARSE = node(_HERE, "test_wall_condition_parse")


# ── gate 1: the parse (law L-B) ──────────────────────────────────────────────


@pytest.mark.rests_on(node(_NETWORK, "test_chains_compose_and_the_family_is_closed"))
def test_wall_condition_parse() -> None:
    """The cone is a parse: ``specular(0.0|0.5|1.0)`` construct; ``specular(1.2)`` raises
    ``SubmarkovViolationError``; ``white(-0.1)`` raises ``BoundaryResponseNotPositiveError``; NaN and a ``bool``
    albedo are refused; a negative inflow is refused; the facades return a ``WallCondition`` (``type(...) is``, not
    a subclass); ``specular(0.0) == white(0.0) == vacuum()`` with ONE digest.

    Reddens on: a copy without the cone parse (``specular(1.2)`` constructs); without the sign check
    (``white(-0.1)`` constructs); one that keeps a zero weight in the network or ``-0.0`` uncanonicalised (the three
    zero walls then differ in digest); one that accepts ``True`` as an albedo.
    """
    for albedo in (0.0, 0.5, 1.0):
        for factory in (bc.WallCondition.specular, bc.WallCondition.white):
            made = factory(albedo)
            require(type(made) is bc.WallCondition, f"{factory.__name__}({albedo}) is a {type(made).__name__}")
            require(made.return_fraction == albedo, f"{factory.__name__}({albedo}).return_fraction = {made.return_fraction}")
    require(bc.WallCondition.specular().cone_point == (1.0, 0.0), "the default albedo of specular is 1.0")
    require(bc.WallCondition.white().cone_point == (0.0, 1.0), "the default albedo of white is 1.0")

    for albedo in (1.2, 5.0, 1.0 + 2.0**-52):
        for factory in (bc.WallCondition.specular, bc.WallCondition.white):
            with pytest.raises(bc.SubmarkovViolationError):
                factory(albedo)
    for albedo in (-0.1, -5.0):
        for factory in (bc.WallCondition.specular, bc.WallCondition.white):
            with pytest.raises(bc.BoundaryResponseNotPositiveError):
                factory(albedo)
    require(issubclass(bc.SubmarkovViolationError, bc.BoundaryError) and issubclass(bc.BoundaryResponseNotPositiveError, bc.BoundaryError), "the cone errors are BoundaryErrors")
    require(issubclass(bc.BoundaryError, ValueError), "BoundaryError is a ValueError")
    for factory in (bc.WallCondition.specular, bc.WallCondition.white):
        with pytest.raises((TypeError, ValueError)):
            factory(float("nan"))
        with pytest.raises(TypeError):
            factory(True)
    with pytest.raises(bc.BoundaryError):
        bc.WallCondition.uniform(-1.0)

    zero = (bc.WallCondition.specular(0.0), bc.WallCondition.white(0.0), bc.WallCondition.vacuum(), bc.WallCondition.uniform(0.0))
    require(all(z == zero[0] for z in zero), f"the zero walls differ: {[z == zero[0] for z in zero]}")
    require(len({content_digest(z) for z in zero}) == 1, "the zero walls have more than one digest (the twin is not gone)")
    require(zero[0].cone_point == (0.0, 0.0) and zero[0].source.is_zero, "vacuum is the empty network over no source")
    require(bc.WallCondition.specular(-0.0) == bc.WallCondition.vacuum(), "-0.0 is canonical 0.0")

    profile = ProfileSampler()
    facades = (
        bc.WallCondition.vacuum(), bc.WallCondition.specular(0.5), bc.WallCondition.white(0.5), bc.WallCondition.uniform(2.0),
        bc.WallCondition.framed(bc.FaceFrameInflow(profile)), bc.WallCondition.framed(bc.LabFrameInflow(profile), 0.5),
        bc.WallCondition.from_amplitudes(0.25, 0.5),
    )
    for made in facades:
        require(type(made) is bc.WallCondition, f"a facade returned a {type(made).__name__}")
    inflow = bc.WallCondition.uniform(2.0)
    require(inflow.cone_point == (0.0, 0.0) and inflow.source.uniform == 2.0 and not inflow.source.sampled, "inflow(v): the zero network over uniform v")
    require(inflow.return_fraction == 0.0 and inflow.face_local is True, "a uniform inflow returns nothing and is face-local")
    framed = bc.WallCondition.framed(bc.FaceFrameInflow(profile), 0.5)
    require(framed.cone_point == (0.0, 0.0) and len(framed.source.sampled) == 1 and framed.source.sampled[0][0] == 0.5, "framed(p, w): the zero network over one sampled term")
    require(bc.WallCondition.framed(bc.FaceFrameInflow(profile), 0.0) == bc.WallCondition.vacuum(), "a zero-weight profile is not stored")


# ── gate 2: no orientation anywhere ──────────────────────────────────────────


def _population() -> list[type]:
    """The closed population gate 2 ranges over: the condition base and every subclass of it, the factors, the
    chains, the network."""
    classes = walk_subclasses(bc.BoundaryCondition) + list(get_args(bc.Factor)) + list(get_args(bc.Chain))
    classes += [bc.ResponseNetwork, bc.InflowTerms, bc.WallDisplacement, bc.FaceFrameInflow, bc.LabFrameInflow, bc.ZeroNetwork]
    return list(dict.fromkeys(classes))


@pytest.mark.rests_on(node(_NETWORK, "test_chains_compose_and_the_family_is_closed"))
def test_orientation_is_unspellable() -> None:
    """``MirrorPermutation(axis="y")`` and ``WallCondition.white(0.5).axis`` fail; no ``BoundaryCondition``
    subclass, factor, chain, network, displacement or profile has a field or attribute named ``axis`` or
    ``outward_sign`` (the realizer reads the normal from the face name, one parse, so one condition is correct on
    every face of its tag).

    Controls (the scanner can fire, X1): a planted dataclass with an ``axis`` field is reported; a planted one with
    ``outward_sign`` is reported; a clean planted dataclass is not. The population is listed and has at least the
    six classes the tables name (``WallCondition``, three factors, two chains) plus the network.

    Reddens on: a copy of ``network.py`` giving ``MirrorPermutation`` a field ``axis: str = "x"`` (the scan names
    it, and ``MirrorPermutation(axis="y")`` constructs); a copy giving ``WallCondition`` an ``outward_sign``
    property. At the base tip (old tree, ``[M]``): ``SpecularReturn(axis="y")`` constructs, and 97 test
    constructions of it exist.
    """
    import dataclasses

    @dataclasses.dataclass(frozen=True)
    class PlantedAxis:
        axis: str = "x"

    @dataclasses.dataclass(frozen=True)
    class PlantedSign:
        outward_sign: int = 1

    @dataclasses.dataclass(frozen=True)
    class Clean:
        weight: float = 1.0

    require(orientation_fields(PlantedAxis) == ["axis"], "control: a planted axis field is not reported")
    require(orientation_fields(PlantedSign) == ["outward_sign"], "control: a planted outward_sign field is not reported")
    require(orientation_fields(Clean) == [], "control: a clean dataclass is reported")

    population = _population()
    names = {cls.__name__ for cls in population}
    wanted = {"WallCondition", "MirrorPermutation", "CurrentRetraction", "IsotropicSection", "SpecularChain", "DiffuseChain", "ResponseNetwork"}
    require(wanted <= names, f"the population lacks {sorted(wanted - names)}: {sorted(names)}")
    # The local probes subclass BoundaryCondition too; the population is what the package and the session define.
    local()
    offenders = {cls.__name__: orientation_fields(cls) for cls in population if orientation_fields(cls)}
    require(not offenders, f"orientation is spellable on {offenders}")

    with pytest.raises(TypeError):
        bc.MirrorPermutation(axis="y")  # type: ignore[call-arg]  # the refusal input
    with pytest.raises(AttributeError):
        bc.WallCondition.white(0.5).axis  # type: ignore[attr-defined]  # noqa: B018  # the refusal input


@pytest.mark.rests_on(_PARSE)
def test_no_welded_surface_survives() -> None:
    """None of the welds of ``BoundaryTraceLaw`` survives on the new base (section 2.4's closing paragraph): no
    ``RegistryMixin`` in the MRO of a condition, no ``kind``, ``key``, ``registry``, ``geometry_map``,
    ``response_kernel`` or ``evaluate`` attribute on a condition, a displacement or the terms, and no ``assert_*``
    method taking a quadrature. Source is an ACCESSOR returning ``InflowTerms``, never ``source.evaluate(space)``.

    Controls: the retired surface is named by a planted class carrying the attributes, which the scan reports.

    Reddens on: a copy of ``condition.py`` restoring ``kind`` (a string tag) on ``BoundaryCondition``, or a base
    ``RegistryMixin``.
    """
    banned = {"kind", "key", "registry", "geometry_map", "response_kernel", "evaluate"}

    class Planted:
        kind = "marshak"
        geometry_map = None

    def welds(cls: type) -> list[str]:
        found = sorted(banned & {name for k in cls.__mro__ for name in vars(k)})
        found += [name for name in dir(cls) if name.startswith("assert_")]
        found += [k.__name__ for k in cls.__mro__ if "RegistryMixin" in k.__name__]
        return found

    require(welds(Planted) == ["geometry_map", "kind"], "control: the scan does not see the planted welds")
    population = [bc.BoundaryCondition, bc.WallCondition, bc.WallDisplacement, bc.InflowTerms]
    offenders = {cls.__name__: welds(cls) for cls in population if welds(cls)}
    require(not offenders, f"a weld of the retired base survives: {offenders}")
    terms = bc.WallCondition.uniform(1.0).source
    require(isinstance(terms, bc.InflowTerms), "the source accessor returns InflowTerms")
    require(not hasattr(terms, "evaluate"), "InflowTerms carries no evaluate(space): sampling is the method's verb")


# ── gate 3: white is a facade (law L-G) ──────────────────────────────────────


@pytest.mark.rests_on(_PARSE)
def test_white_is_a_facade() -> None:
    """``white(a) == from_amplitudes(0.0, a)``, equal digests, ``type(...) is WallCondition``, and the first
    positional argument is the albedo (the old ``WhiteBoundary`` took ``axis`` first).

    Reddens on: a copy returning a ``WhiteCondition`` subclass (``type is``); a copy swapping the arguments of
    ``white``; a copy whose ``white`` rebuilds the network through another route that differs in digest.
    """
    for first in (bc.WallCondition.white, bc.WallCondition.specular):
        params = list(inspect.signature(first).parameters)
        require(params == ["albedo"], f"{first.__name__} takes {params}, the albedo alone")
    for albedo in (0.0, 0.125, 0.3, 0.5, 1.0):
        made, edge = bc.WallCondition.white(albedo), bc.WallCondition.from_amplitudes(0.0, albedo)
        require(type(made) is bc.WallCondition and type(edge) is bc.WallCondition, "a facade returned a subclass")
        require(made == edge, f"white({albedo}) != from_amplitudes(0.0, {albedo})")
        require(content_digest(made) == content_digest(edge), f"white({albedo}) and the edge have two digests")
        require(made.cone_point == (0.0, albedo), f"white({albedo}).cone_point = {made.cone_point}")
        specular, edge = bc.WallCondition.specular(albedo), bc.WallCondition.from_amplitudes(albedo, 0.0)
        require(specular == edge and content_digest(specular) == content_digest(edge), f"specular({albedo}) is not the edge")
    require(bc.WallCondition.white(0.4).cone_point == (0.0, 0.4), "white(0.4): the first positional argument is the albedo")
    require(bc.WallCondition.white(0.4) != bc.WallCondition.specular(0.4), "white and specular are two walls")


# ── gate 11: from_amplitudes is total on the cone ────────────────────────────


@pytest.mark.rests_on(_PARSE)
def test_from_amplitudes_is_total_on_the_cone() -> None:
    """``from_amplitudes(0, 0)`` is ``vacuum()``; ``(s, 0)`` is ``specular(s)``; ``(0, d)`` is ``white(d)``;
    ``(0.3, 0.4)`` has ``cone_point == (0.3, 0.4)`` EXACTLY; ``(0.7, 0.7)`` raises the cone refusal; the round trip
    ``cone_point . from_amplitudes = id`` holds over a grid of the cone (no ``isinstance`` in the test: the chart is
    total because the chain family is closed); the cone sum is ``math.fsum`` (``0.9 + 0.1`` is accepted, one double
    step past 1 is refused); a negative amplitude is refused; the edge takes the two primitives and not an object
    (insulation of the reference: ``WallAmplitudes`` is the reference's value and this package imports no
    specification).

    Reddens on: a copy dropping the diffuse weight when ``specular > 0``; a copy that does not canonicalise ``-0.0``;
    one taking a ``WallAmplitudes``; a cone sum that rounds down or compares with ``<`` (the ``0.9, 0.1`` and
    ``1.0, 2**-52`` rows bound the correctly rounded sum from both sides).
    """
    require(bc.WallCondition.from_amplitudes(0, 0) == bc.WallCondition.vacuum(), "(0, 0) is vacuum")
    for s in (0.25, 0.5, 1.0):
        require(bc.WallCondition.from_amplitudes(s, 0.0) == bc.WallCondition.specular(s), f"({s}, 0) is specular({s})")
        require(bc.WallCondition.from_amplitudes(0.0, s) == bc.WallCondition.white(s), f"(0, {s}) is white({s})")
    mixed = bc.WallCondition.from_amplitudes(0.3, 0.4)
    require(mixed.cone_point == (0.3, 0.4), f"(0.3, 0.4) reads back as {mixed.cone_point}")
    require(type(mixed) is bc.WallCondition, "from_amplitudes returns a WallCondition")

    population = cone_grid()
    require(len(population) == 48, f"the cone grid has {len(population)} points, not 48")
    for s, d in population:
        point = wall(s, d).cone_point
        require(point == (s, d), f"cone_point . from_amplitudes is not the identity at {(s, d)}: {point}")

    require(math.fsum([0.9, 0.1]) <= 1.0, "control: the cone sum of 0.9 and 0.1 is correctly rounded to 1")
    require(bc.WallCondition.from_amplitudes(0.9, 0.1).cone_point == (0.9, 0.1), "(0.9, 0.1) is accepted: an intended unit wall")
    require(bc.WallCondition.from_amplitudes(0.25, 0.75).cone_point == (0.25, 0.75), "(0.25, 0.75) is a dyadic unit wall")
    for s, d in ((0.7, 0.7), (1.0, 2.0**-52), (1.0, 1.0), (2.0, 0.0), (0.0, 1.5)):
        with pytest.raises(bc.SubmarkovViolationError):
            bc.WallCondition.from_amplitudes(s, d)
    for s, d in ((-0.1, 0.5), (0.5, -0.1), (-0.5, -0.5)):
        with pytest.raises(bc.BoundaryResponseNotPositiveError):
            bc.WallCondition.from_amplitudes(s, d)

    params = list(inspect.signature(bc.WallCondition.from_amplitudes).parameters)
    require(params == ["specular", "diffuse"], f"from_amplitudes takes {params}")
    carrier = type("PlainAmplitudes", (), {"specular": 0.3, "diffuse": 0.4})()
    with pytest.raises(TypeError):
        bc.WallCondition.from_amplitudes(carrier, 0.0)  # type: ignore[arg-type]  # an object is not a primitive
