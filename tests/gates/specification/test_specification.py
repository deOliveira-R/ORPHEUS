r"""The specification's composition laws (#405 P1 step 8: S8.1, S8.10).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8, under the step-8 rulings and the user's ruling "The infinite medium is
the point in phase space"). The module is ``orpheus/specification``.

The layer a specification is posed at is its TYPE, ``Specification =
InfiniteMediumSpecification | GeometrySpecification``.
``InfiniteMediumSpecification(material_id, mixture, question)`` is posed on
energy alone: it holds one material and has no field a geometry could be given
in. ``GeometrySpecification(materials, geometry, question)`` keeps the
materials its geometry assigns (a spectator is dropped). At construction each
checks the question's datum and keys and stores the CANONICAL question: every
key (the ``Eigen`` parameter and each point key) is replaced by its resolved
explicit form, so ``spec.question`` may differ from, and be unequal to, the
question the caller passed (S8.10). A key is a ``CellCoefficient`` or a
``GeometryExtent``; a resolved ``CellCoefficient`` names no channel "in every
material".

Every refusal is keyed (``match=`` on a fragment the triggering argument
determines, the class name included), and the fragments are disjoint, asserted
once (``test_s8_1_the_refusal_fragments_are_disjoint``). Where two refusals can
fire on one input, a discrimination row pins the wiring order.
"""

from __future__ import annotations

import dataclasses
import pickle
import re
import typing
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
import sympy as sp

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import InconsistentMaterialsError, Materials
from orpheus.geometry import BC, CoordSystem, GeometryExtent, StructuredGeometry
from orpheus.geometry.boundary import PrescribedInflow
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.content import ContentlessError
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.question import Eigen, FixedSource, Nearest, Response
from orpheus.specification import (
    Coordinate,
    GeometrySpecification,
    InfiniteMediumSpecification,
    Specification,
)
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import (
    body,
    fission_only,
    fuel,
    moderator,
    n2n_only,
    slab2,
    slab3_repeated,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/specification/test_specification.py"
_Q = "tests/gates/numerics/test_question.py"
_MF = "tests/gates/numerics/test_mesh_free_function.py"
_CELLS = "tests/gates/data/test_cells.py"
_EXT = "tests/gates/geometry/test_geometry_extent.py"

F, S, N2 = Channel.FISSION_EMISSION, Channel.SCATTERING_EMISSION, Channel.N2N_EMISSION
r, mu, phi = Symbolic.r, Symbolic.mu, Symbolic.phi


def _k() -> CellCoefficient:
    return CellCoefficient.every(F)


def _table(regions: int = 2, groups: int = 2) -> RegionwiseConstant:
    return RegionwiseConstant(np.arange(1.0, 1.0 + regions * groups).reshape(regions, groups))


def _two_materials() -> Materials:
    return Materials({0: fuel(), 1: moderator()})


def _im(question: Any, mixture: Any = None, material_id: Any = 0) -> InfiniteMediumSpecification:
    """The infinite medium of one material (fuel unless given)."""
    return InfiniteMediumSpecification(material_id, fuel() if mixture is None else mixture, question)


def _gs(question: Any, materials: Any = None, geometry: Any = None) -> GeometrySpecification:
    """A geometry problem (the two-material slab unless given)."""
    return GeometrySpecification(_two_materials() if materials is None else materials,
                                 slab2() if geometry is None else geometry, question)


# ═════════════════════════════════════════════════════════════════════════════
# S8.1 — the refusals, keyed
# ═════════════════════════════════════════════════════════════════════════════
#
# One row per rule: (id, builder, the exception type, the fragment). The
# builder is CALLED in the test body, never at collection (vv #17(c): a
# mutation making production raise at collection kills the module).

_REFUSALS: list[tuple[str, Callable[[], Any], type[BaseException], str]] = [
    # (a) "the infinite medium has one material" is unspellable: its type holds one
    # mixture (the fields row below). (b) a geometry id missing from the materials
    # (the pinned fragment of Materials.restrict)
    ("b_geometry_id_not_declared",
     lambda: _gs(Eigen(_k()), materials=Materials({0: fuel()})),
     ValueError, r"references material ids \[1\]"),
    # (c) group counts differ across the materials the geometry assigns
    ("c_group_count",
     lambda: _gs(Eigen(_k()), materials=Materials({0: fuel(), 1: fuel(ng=3)})),
     InconsistentMaterialsError, r"uniform ng"),
    # (g) the datum against the problem's regions, groups and coordinates
    ("g_source_regions",
     lambda: _gs(FixedSource(_table(regions=3))),
     ValueError, r"GeometrySpecification: the source has 3 regions; the problem has 2$"),
    ("g_detector_regions",
     lambda: _gs(Response(_table(regions=1))),
     ValueError, r"GeometrySpecification: the detector has 1 regions; the problem has 2$"),
    # qa F1: three intervals over two distinct materials; the count is the intervals'
    ("g_regions_are_intervals_not_materials",
     lambda: _gs(FixedSource(_table(regions=2)), geometry=slab3_repeated()),
     ValueError, r"GeometrySpecification: the source has 2 regions; the problem has 3$"),
    ("g_source_groups_table",
     lambda: _gs(FixedSource(_table(groups=3))),
     ValueError, r"GeometrySpecification: the source has 3 groups; the materials have 2"),
    ("g_detector_groups_symbolic",
     lambda: _gs(Response(Symbolic.of(r))),
     ValueError, r"GeometrySpecification: the detector has 1 groups; the materials have 2"),
    ("g_infinite_medium_regions",
     lambda: _im(FixedSource(_table(regions=2))),
     ValueError, r"InfiniteMediumSpecification: the source has 2 regions; the problem has 1$"),
    ("g_infinite_medium_position",
     lambda: _im(FixedSource(Symbolic.of(r, 1))),
     ValueError, r"InfiniteMediumSpecification: the source depends on r, which this problem has no coordinate to read"),
    ("g_infinite_medium_direction",
     lambda: _im(Response(Symbolic.of(1, sp.cos(phi) * mu))),
     ValueError, r"InfiniteMediumSpecification: the detector depends on mu, phi, which this problem has no coordinate to read"),
    ("g_sphere_azimuth",
     lambda: _gs(FixedSource(Symbolic.of(1 + sp.cos(phi), 1)), geometry=body(CoordSystem.SPHERICAL)),
     ValueError, r"GeometrySpecification: the source depends on phi, which this problem has no coordinate to read"),
    # (h1) a parameter that is not a coordinate
    ("h1_parameter_str",
     lambda: _im(Eigen("fission-emission")),
     TypeError, r"InfiniteMediumSpecification: the parameter: 'fission-emission' \(a str\) is not a coordinate"),
    ("h1_parameter_none",
     lambda: _gs(Eigen(None)),
     TypeError, r"GeometrySpecification: the parameter: None \(a NoneType\) is not a coordinate"),
    # (h2) a point key that is not a coordinate, for every kind of question
    ("h2_point_key_eigen",
     lambda: _im(Eigen(_k(), {"boron": 0.1})),
     TypeError, r"the point key: 'boron' \(a str\) is not a coordinate"),
    ("h2_point_key_fixed_source",
     lambda: _im(FixedSource(_table(1), {"boron": 0.1})),
     TypeError, r"the point key: 'boron' \(a str\) is not a coordinate"),
    ("h2_point_key_response",
     lambda: _gs(Response(_table(), {("cell", 0): 0.1})),
     TypeError, r"the point key: \('cell', 0\) \(a tuple\) is not a coordinate"),
    # (h3) an extent on the infinite medium (as parameter and as point key), or out of range
    ("h3_extent_infinite_medium_parameter",
     lambda: _im(Eigen(GeometryExtent(0))),
     ValueError, r"InfiniteMediumSpecification: GeometryExtent\(interval=0\) names an extent, and the infinite medium has no geometry"),
    ("h3_extent_infinite_medium_point",
     lambda: _im(FixedSource(_table(1), {GeometryExtent(1): 0.5})),
     ValueError, r"InfiniteMediumSpecification: GeometryExtent\(interval=1\) names an extent, and the infinite medium has no geometry"),
    ("h3_extent_out_of_range",
     lambda: _gs(Eigen(GeometryExtent(2))),
     ValueError, r"interval 2 does not exist; the geometry has 2 interval"),
    # qa F1: the range is the intervals', not the distinct materials'
    ("h3_extent_out_of_range_repeated_ids",
     lambda: _gs(Eigen(GeometryExtent(3)), geometry=slab3_repeated()),
     ValueError, r"interval 3 does not exist; the geometry has 3 interval"),
    # (h4) a cell outside the problem's materials (a spectator included, qa F3), or a zero direction
    ("h4_material_not_in_the_problem",
     lambda: _gs(Eigen(CellCoefficient({(7, F)}))),
     ValueError, r"material 7 is not among the problem's materials \(ids: \[0, 1\]\)"),
    ("h4_spectator_cell",
     lambda: _gs(Eigen(CellCoefficient({(9, F)})), materials=Materials({0: fuel(), 1: moderator(), 9: fuel()})),
     ValueError, r"material 9 is not among the problem's materials \(ids: \[0, 1\]\)"),
    ("h4_k_on_a_non_producing_specification",
     lambda: _im(Eigen(_k()), mixture=moderator()),
     ValueError, r"zero direction"),
    # (h5) two point keys naming one coordinate
    ("h5_one_coordinate_twice",
     lambda: _im(FixedSource(_table(1), {_k(): 0.1, CellCoefficient({(0, F)}): 0.2})),
     ValueError, r"name the same coordinate"),
    # a specification exists to be keyed: a geometry with no content is refused eagerly
    ("contentless_geometry",
     lambda: _gs(Eigen(_k()), materials=Materials({0: fuel()}),
                 geometry=StructuredGeometry.slab((0.0, 1.0), (0,), left=PrescribedInflow(lambda space: np.zeros(space.shape)),  # pyright: ignore[reportArgumentType]  # a contentless source is the subject
                                                  right=BC.vacuum)),
     ContentlessError, r"GeometrySpecification\.geometry"),
    # the fields are typed
    ("question_not_a_question",
     lambda: _im(_table(1)),
     TypeError, r"InfiniteMediumSpecification: the question is an Eigen or FixedSource or Response, got a RegionwiseConstant"),
    ("geometry_none",
     lambda: GeometrySpecification(_two_materials(), None, Eigen(_k())),  # pyright: ignore[reportArgumentType]  # the refusal is the subject
     TypeError, r"GeometrySpecification: the geometry is a StructuredGeometry, got a NoneType"),
    ("geometry_tuple",
     lambda: GeometrySpecification(_two_materials(), (0.0, 1.0), Eigen(_k())),  # pyright: ignore[reportArgumentType]  # the refusal is the subject
     TypeError, r"GeometrySpecification: the geometry is a StructuredGeometry, got a tuple"),
    ("materials_dict",
     lambda: GeometrySpecification({0: fuel(), 1: moderator()}, slab2(), Eigen(_k())),  # pyright: ignore[reportArgumentType]  # the refusal is the subject
     TypeError, r"GeometrySpecification: the materials are a Materials, got a dict"),
    ("mixture_not_a_mixture",
     lambda: InfiniteMediumSpecification(0, Materials({0: fuel()}), Eigen(_k())),  # pyright: ignore[reportArgumentType]  # the refusal is the subject
     TypeError, r"InfiniteMediumSpecification: the mixture is a Mixture, got a Materials"),
    ("material_id_bool",
     lambda: InfiniteMediumSpecification(True, fuel(), Eigen(_k())),
     TypeError, r"InfiniteMediumSpecification: the material id is an int, got bool"),
]


@pytest.mark.rests_on(f"{_Q}::test_s7_3_the_signatures", f"{_CELLS}::test_s8_7_resolution_by_channel")
@pytest.mark.parametrize("build,error,fragment", [pytest.param(b, e, f, id=i) for i, b, e, f in _REFUSALS])
def test_s8_1_each_refusal_is_keyed(build, error: type[BaseException], fragment: str) -> None:
    with pytest.raises(error, match=fragment):
        build()


def test_s8_1_the_refusal_fragments_are_disjoint() -> None:
    """Each refusal's message matches its own fragment and NO other row's, so a
    row cannot be satisfied by a neighbouring rule (test-architect lessons §1)."""
    messages: dict[str, str] = {}
    for rid, build, error, _ in _REFUSALS:
        try:
            build()
        except error as err:
            messages[rid] = str(err)
        else:
            pytest.fail(f"{rid}: constructed")
    require(len(messages) == len(_REFUSALS), f"activation: {len(messages)} of {len(_REFUSALS)} refusals read")
    # Rows that share a rule share its fragment by design (one rule, several triggers).
    family = lambda rid: rid.split("_", 1)[0]  # noqa: E731
    overlaps = [
        f"{rid}'s message matches {other}'s fragment {frag!r}"
        for rid, msg in messages.items()
        for other, _, _, frag in _REFUSALS
        if other != rid and family(other) != family(rid) and re.search(frag, msg)
    ]
    require(not overlaps, "\n".join(overlaps))


# ── the wiring order: inputs that violate two rules ──────────────────────────


def test_s8_1_the_group_count_is_read_before_the_keys() -> None:
    """Two kept materials of different group counts AND a key naming a material
    outside the problem: (c) fires, (h4) does not."""
    with pytest.raises(InconsistentMaterialsError, match="uniform ng") as err:
        _gs(Eigen(CellCoefficient({(7, F)})), materials=Materials({0: fuel(), 1: fuel(ng=3)}))
    require("material 7" not in str(err.value), f"(h4) fired with (c): {err.value}")


def test_s8_1_a_spectator_of_another_group_count_is_dropped() -> None:
    """A declared material the geometry does not assign is dropped before the
    group count is read (the leak principle: a spectator changes no answer)."""
    spec = _gs(Eigen(_k()), materials=Materials({0: fuel(), 5: fuel(ng=3)}), geometry=slab2((0, 0)))
    require(sorted(spec.materials.ids) == [0], f"the spectator was kept: {sorted(spec.materials.ids)}")


def test_s8_1_the_composition_is_read_before_the_keys() -> None:
    """The geometry names an undeclared id AND the key names an undeclared
    material: (b) fires, (h4) does not (keys resolve against a valid declaration)."""
    with pytest.raises(ValueError, match=r"references material ids \[1\]") as err:
        _gs(Eigen(CellCoefficient({(7, F)})), materials=Materials({0: fuel()}))
    require("material 7" not in str(err.value), f"(h4) fired with (b): {err.value}")


# ── (d)-(f): the refused states of the old table are unspellable ─────────────


@pytest.mark.rests_on(f"{_Q}::test_s7_3_the_fields_are_exactly_the_roles")
def test_s8_1_the_layer_is_the_type() -> None:
    """The structural legs. The old (d)-(f): the datum lives in the question, so
    neither type has a ``source`` field. The old (a), "the infinite medium has
    one material", and "the infinite medium has no geometry": the infinite-medium
    type holds one mixture and no geometry field (the user's ruling, "The infinite
    medium is the point in phase space"). ``Specification`` is exactly the two
    types, and both expose the problem's ``materials``. Witness: a field
    re-added to either type reds this row."""
    infinite = tuple(f.name for f in dataclasses.fields(InfiniteMediumSpecification))
    geometric = tuple(f.name for f in dataclasses.fields(GeometrySpecification))
    require(infinite == ("material_id", "mixture", "question"), f"InfiniteMediumSpecification fields {infinite}")
    require(geometric == ("materials", "geometry", "question"), f"GeometrySpecification fields {geometric}")
    require(set(typing.get_args(Specification)) == {InfiniteMediumSpecification, GeometrySpecification},
            f"Specification is {typing.get_args(Specification)}")
    require(set(typing.get_args(Coordinate)) == {CellCoefficient, GeometryExtent}, f"Coordinate is {typing.get_args(Coordinate)}")
    medium = _im(Eigen(_k()), material_id=4)
    require(medium.materials == Materials({4: fuel()}), f"the infinite medium's materials are {medium.materials!r}")
    require(medium.n_regions == 1 and InfiniteMediumSpecification.n_regions == 1, "the infinite medium is one region")
    require(medium.unreadable == (r, mu, phi), f"unreadable {medium.unreadable}")


@pytest.mark.parametrize(
    "coord,unreadable",
    [(CoordSystem.CARTESIAN, ()), (CoordSystem.CYLINDRICAL, ()), (CoordSystem.SPHERICAL, (phi,))],
    ids=["slab", "cylinder", "sphere"],
)
def test_s8_1_a_geometry_problem_reads_what_its_chart_can(coord, unreadable) -> None:
    """``unreadable`` is read off the chart (``azimuth_reference is None`` on the sphere
    only), and ``n_regions`` is the interval count, not the distinct-material count (qa F1)."""
    spec = _gs(Eigen(_k()), geometry=body(coord))
    require(spec.unreadable == unreadable, f"{coord}: unreadable {spec.unreadable}")
    require(_gs(Eigen(_k()), geometry=slab3_repeated()).n_regions == 3, "n_regions counts the distinct materials")


# ── the fixtures are consistent (vv-testing: a hand-built mixture is gated) ──


@pytest.mark.parametrize("build", [fuel, moderator, n2n_only, fission_only], ids=lambda b: b.__name__)
def test_s8_1_each_fixture_mixture_balances(build) -> None:
    """sigma_t = sigma_c + sigma_l + sigma_f + sum_to SigS[0][g, :] + sum_to Sig2[0][g, :], bitwise as built."""
    m = build()
    removal = m.SigC + m.SigL + m.SigF + np.asarray(m.SigS[0].sum(axis=1)).ravel() + np.asarray(m.Sig2[0].sum(axis=1)).ravel()
    require(np.allclose(m.SigT, removal, rtol=4 * np.finfo(float).eps, atol=0), f"{build.__name__}: {m.SigT} vs {removal}")
    require(bool(np.all(m.SigC >= 0)), "a negative capture")


# ── the admissions: one positive leg per rule ────────────────────────────────

_ADMISSIONS: list[tuple[str, Callable[[], Specification]]] = [
    ("k on a fissile slab", lambda: _gs(Eigen(_k()))),
    ("k in the infinite medium", lambda: _im(Eigen(_k()))),
    ("c, every emission cell", lambda: _gs(Eigen(CellCoefficient.every(F, S, N2)))),
    ("a critical extent", lambda: _gs(Eigen(GeometryExtent(1), mode=Nearest(1.0)))),
    ("the last interval of a repeated-id geometry (qa F1)",
     lambda: _gs(Eigen(GeometryExtent(2)), geometry=slab3_repeated())),
    ("three regions on a repeated-id geometry (qa F1)",
     lambda: _gs(FixedSource(_table(regions=3)), geometry=slab3_repeated())),
    ("a spectator material", lambda: _gs(Eigen(_k()), materials=Materials({**_two_materials(), 9: moderator()}))),
    ("a fixed source on a non-producing slab, offset along scattering",
     lambda: _gs(FixedSource(_table(), {CellCoefficient.every(S): 0.1}), materials=Materials({0: moderator(), 1: moderator()}))),
    ("a response offset along an extent", lambda: _gs(Response(_table(), {GeometryExtent(0): -0.2}))),
    ("the parameter's own offset (step-7 ruling 3)", lambda: _im(Eigen(_k(), {_k(): 0.3}))),
    ("an explicit zero cell beside a non-zero one (ruling 5)", lambda: _gs(Eigen(CellCoefficient({(0, F), (1, F)})))),
    ("a constant Symbolic in the infinite medium", lambda: _im(FixedSource(Symbolic.of(2, 1)))),
    ("one region in the infinite medium", lambda: _im(Response(_table(1)))),
    ("any material id in the infinite medium", lambda: _im(Eigen(_k()), material_id=np.int64(7))),
]


@pytest.mark.parametrize("build", [pytest.param(b, id=i) for i, b in _ADMISSIONS])
def test_s8_1_each_rule_admits_its_legal_inputs(build) -> None:
    spec = build()
    require(isinstance(spec, Specification), "did not construct")


# ── (g): the phi-dependence on a sphere, by coordinate system ────────────────

_AZIMUTH = [
    # (coordinate system, expression, admitted)
    (CoordSystem.SPHERICAL, sp.cos(phi), False),
    (CoordSystem.SPHERICAL, sp.Piecewise((1, phi < sp.pi), (0, True)), False),
    (CoordSystem.SPHERICAL, mu, True),                              # mu is read against the radial direction
    (CoordSystem.SPHERICAL, sp.sin(phi) ** 2 + sp.cos(phi) ** 2, True),  # simplify decides it constant in phi
    (CoordSystem.SPHERICAL, r**2 * (1 + mu), True),
    (CoordSystem.CARTESIAN, sp.cos(phi), True),
    (CoordSystem.CYLINDRICAL, sp.cos(phi), True),
]


@pytest.mark.rests_on(f"{_MF}::test_s8_11_the_azimuth_predicate")
@pytest.mark.parametrize(
    "coord,expression,admitted", _AZIMUTH,
    ids=["sph-cos-phi", "sph-step-phi", "sph-mu", "sph-sin2cos2", "sph-r2mu", "slab-cos-phi", "cyl-cos-phi"],
)
@pytest.mark.parametrize("role", ["source", "detector"])
def test_s8_1_a_phi_dependent_function_is_refused_beside_a_sphere(coord, expression, admitted: bool, role: str) -> None:
    question = (FixedSource if role == "source" else Response)(Symbolic.of(expression, expression))
    build = lambda: _gs(question, geometry=body(coord))  # noqa: E731
    if admitted:
        build()
    else:
        with pytest.raises(ValueError, match=rf"the {role} depends on phi, which this problem has no coordinate to read"):
            build()


# ── (h4): the channel carriage matrix, through the specification ─────────────

_CARRIERS = {F: fission_only, S: moderator, N2: n2n_only}


@pytest.mark.rests_on(f"{_CELLS}::test_s8_6_the_carriage_predicate")
@pytest.mark.parametrize("asked", [F, S, N2], ids=lambda c: c.name)
@pytest.mark.parametrize("carried", [F, S, N2], ids=lambda c: c.name)
def test_s8_1_a_direction_resolves_iff_some_material_carries_it(asked: Channel, carried: Channel) -> None:
    """Each mixture carries exactly one emission channel, so the 3 x 3 matrix
    separates the three carriage predicates: a predicate read from the wrong
    stack reds an off-diagonal or a diagonal cell."""
    build = lambda: _im(Eigen(CellCoefficient.every(asked)), mixture=_CARRIERS[carried]())  # noqa: E731
    if asked is carried:
        spec = build()
        require(_parameter(spec) == CellCoefficient({(0, asked)}), f"resolved to {_parameter(spec)}")
    else:
        with pytest.raises(ValueError, match="zero direction"):
            build()


# ═════════════════════════════════════════════════════════════════════════════
# S8.10 — the stored question is the canonical question
# ═════════════════════════════════════════════════════════════════════════════


def _parameter(spec: Specification) -> Any:
    """The stored parameter of an eigen question (narrowed for the type checker)."""
    require(isinstance(spec.question, Eigen), f"the question is a {type(spec.question).__name__}")
    return getattr(spec.question, "parameter")


@pytest.mark.rests_on(f"{_CELLS}::test_s8_7_resolution_by_channel")
def test_s8_10_every_resolves_to_the_explicit_non_zero_cells() -> None:
    """``every(F)`` on {0: fuel, 1: moderator} is stored as {(0, F)}: the
    quantifier is gone and the moderator's zero cell is not a cell of the key."""
    q = Eigen(_k(), {CellCoefficient.every(S): 0.25})
    spec = _gs(q)
    require(_parameter(spec) == CellCoefficient({(0, F)}), f"parameter {_parameter(spec)}")
    require(dict(spec.question.point) == {CellCoefficient({(0, S), (1, S)}): 0.25}, f"point {spec.question.point}")
    keys = [_parameter(spec), *spec.question.point]
    require(all(k.channels_in_every_material == frozenset() for k in keys), f"a channel in every material survived in {keys}")
    require(spec.question != q, "activation: the stored question equals the written one")
    require(getattr(spec.question, "mode") == q.mode, "the mode moved")


@pytest.mark.parametrize(
    "written",
    [
        pytest.param(lambda: Eigen(CellCoefficient({(0, F)})), id="explicit"),
        pytest.param(lambda: Eigen(CellCoefficient({(0, F), (1, F)})), id="explicit with a zero cell"),
        pytest.param(lambda: Eigen(CellCoefficient([(0, F), (0, F)])), id="a duplicate cell"),  # pyright: ignore[reportArgumentType]  # the coercion is the subject
        pytest.param(lambda: Eigen(CellCoefficient({(np.int64(0), F)})), id="a numpy id"),  # pyright: ignore[reportArgumentType]  # the coercion is the subject
    ],
)
def test_s8_10_two_spellings_of_one_direction_are_one_specification(written) -> None:
    """Equal, one digest, one set member, and the stored key is the explicit cell set."""
    a = _gs(Eigen(_k()))
    b = _gs(written())
    require(a == b and a.content_digest == b.content_digest and len({a, b}) == 1, "two specifications")
    require(_parameter(b) == CellCoefficient({(0, F)}), f"stored {_parameter(b)}")


def test_s8_10_a_zero_cell_in_a_point_key_is_dropped() -> None:
    """qa F4: the canonical form drops a zero cell in a POINT key as in the parameter."""
    spec = _gs(FixedSource(_table(), {CellCoefficient({(0, F), (1, F)}): 0.1}))
    require(dict(spec.question.point) == {CellCoefficient({(0, F)}): 0.1}, f"point {dict(spec.question.point)}")


@pytest.mark.parametrize("kind", ["infinite medium", "geometry"])
def test_s8_10_canonicalisation_is_idempotent_and_survives_pickle(kind: str) -> None:
    """Re-posing from the STORED question returns an equal specification (resolution
    is idempotent). Re-posing with other materials must start from the caller's
    question (qa F2, documented in the module, not guarded)."""
    if kind == "geometry":
        spec = _gs(FixedSource(_table(), {CellCoefficient.every(S, N2): -0.5, GeometryExtent(1): 0.125}))
        again = dataclasses.replace(spec, question=spec.question)
        offsets = [-0.5, 0.125]
    else:
        spec = _im(Response(_table(1), {CellCoefficient.every(S, N2): -0.5}))
        again = dataclasses.replace(spec, question=spec.question)
        offsets = [-0.5]
    require(again == spec and again.question == spec.question, "re-posing the stored question moved it")
    back = pickle.loads(pickle.dumps(spec))
    require(back == spec and back.content_digest == spec.content_digest, "pickle moved the specification")
    require(list(spec.question.point.values()) == offsets, "the offsets moved")


def test_s8_10_the_two_layers_are_two_keys() -> None:
    """The infinite medium of a material and a reflective slab of it ask the same k
    and are two specifications (two types, two digests): the layer is content."""
    medium = _im(Eigen(_k()))
    slab = _gs(Eigen(_k()), materials=Materials({0: fuel()}),
               geometry=StructuredGeometry.from_homogeneous(1.0, BC.reflective))
    require(medium != slab and medium.content_digest != slab.content_digest and len({medium, slab}) == 2, "one key")
    require(medium.question == slab.question, "activation: the two do not ask one question")


def test_s8_10_the_materials_are_the_declaration_value() -> None:
    """One spelling: the materials are a ``Materials`` value, and a bare mapping
    is refused (the orchestrator's ruling 2 on the step-8 NEEDS, 2026-10-02)."""
    with pytest.raises(TypeError, match="the materials are a Materials, got a dict"):
        GeometrySpecification({0: fuel(), 1: moderator()}, slab2(), Eigen(_k()))  # pyright: ignore[reportArgumentType]  # the refusal is the subject


