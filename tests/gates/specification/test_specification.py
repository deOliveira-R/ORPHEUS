r"""The specification's composition laws (#405 P1 step 8: S8.1, S8.10).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8, under the step-8 rulings). The module is ``orpheus/specification``.

``Specification(materials, geometry | None, question)`` composes a materials
declaration, a geometry (``None`` is the infinite medium) and a question.
At construction it checks the composition, the question's functions and the
question's keys, and it stores the CANONICAL question: every key (the
``Eigen`` parameter and each point key) is replaced by its resolved explicit
form, so ``spec.question`` may differ from, and be unequal to, the question
the caller passed (S8.10). A key is a ``CellCoefficient`` (a set of
``(material id, Channel)`` cells) or a ``GeometryExtent(interval)``; the
quantifier ``CellCoefficient.every(...)`` never reaches a stored key.

Every refusal is keyed (``match=`` on a fragment the triggering argument
determines), and the fragments are disjoint, asserted once
(``test_s8_1_the_refusal_fragments_are_disjoint``). Where two refusals can
fire on one input, a discrimination row pins the wiring order.
"""

from __future__ import annotations

import dataclasses
import pickle
import re
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
from orpheus.specification import Specification
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import (
    body,
    fission_only,
    fuel,
    moderator,
    n2n_only,
    slab2,
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


# ═════════════════════════════════════════════════════════════════════════════
# S8.1 — the refusals, keyed
# ═════════════════════════════════════════════════════════════════════════════
#
# One row per rule: (id, builder, the exception type, the fragment). The
# builder is CALLED in the test body, never at collection (vv #17(c): a
# mutation making production raise at collection kills the module).

_REFUSALS: list[tuple[str, Callable[[], Any], type[BaseException], str]] = [
    # (a) the infinite medium has one material
    ("a_no_geometry_two_materials",
     lambda: Specification(materials=_two_materials(), geometry=None, question=Eigen(_k())),
     ValueError, r"the infinite medium has one material; 2 are declared"),
    # (b) a geometry id missing from the materials (the pinned fragment of Materials.restrict)
    ("b_geometry_id_not_declared",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=slab2(), question=Eigen(_k())),
     ValueError, r"references material ids \[1\]"),
    # (c) group counts differ across the materials the geometry assigns (the
    # orchestrator's ruling 3 on the step-8 NEEDS: a spectator is dropped first)
    ("c_group_count",
     lambda: Specification(materials=Materials({0: fuel(), 1: fuel(ng=3)}), geometry=slab2(), question=Eigen(_k())),
     InconsistentMaterialsError, r"uniform ng"),
    # (g) the functions against the geometry and the materials
    ("g_source_regions",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=FixedSource(_table(regions=3))),
     ValueError, r"the source has 3 regions; the geometry has 2 intervals"),
    ("g_detector_regions",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Response(_table(regions=1))),
     ValueError, r"the detector has 1 regions; the geometry has 2 intervals"),
    ("g_source_groups_table",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=FixedSource(_table(groups=3))),
     ValueError, r"the source has 3 groups; the materials have 2"),
    ("g_detector_groups_symbolic",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Response(Symbolic.of(r))),
     ValueError, r"the detector has 1 groups; the materials have 2"),
    ("g_infinite_medium_regions",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=FixedSource(_table(regions=2))),
     ValueError, r"the source has 2 regions; the infinite medium has 1 region"),
    ("g_infinite_medium_symbolic_position",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=FixedSource(Symbolic.of(r, 1))),
     ValueError, r"no coordinate system"),
    ("g_sphere_azimuth",
     lambda: Specification(materials=_two_materials(), geometry=body(CoordSystem.SPHERICAL),
                           question=FixedSource(Symbolic.of(1 + sp.cos(phi), 1))),
     ValueError, r"azimuth"),
    # (h1) a parameter that is not a coordinate
    ("h1_parameter_str",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen("fission-emission")),
     TypeError, r"parameter: 'fission-emission' \(a str\) is not a coordinate"),
    ("h1_parameter_none",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(None)),
     TypeError, r"parameter: None \(a NoneType\) is not a coordinate"),
    # (h2) a point key that is not a coordinate, for every kind of question
    ("h2_point_key_eigen",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(_k(), {"boron": 0.1})),
     TypeError, r"point key: 'boron' \(a str\) is not a coordinate"),
    ("h2_point_key_fixed_source",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=FixedSource(_table(1), {"boron": 0.1})),
     TypeError, r"point key: 'boron' \(a str\) is not a coordinate"),
    ("h2_point_key_response",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Response(_table(1), {("cell", 0): 0.1})),
     TypeError, r"point key: \('cell', 0\) \(a tuple\) is not a coordinate"),
    # (h3) an extent with no geometry, or out of range, as parameter and as point key
    ("h3_extent_no_geometry_parameter",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(GeometryExtent(0))),
     ValueError, r"parameter: a GeometryExtent names an extent, and the infinite medium \(no geometry\) has none"),
    ("h3_extent_no_geometry_point",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=FixedSource(_table(1), {GeometryExtent(0): 0.5})),
     ValueError, r"point key: a GeometryExtent names an extent, and the infinite medium \(no geometry\) has none"),
    ("h3_extent_out_of_range",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Eigen(GeometryExtent(2))),
     ValueError, r"interval 2 does not exist; the geometry has 2 interval"),
    # (h4) a cell coefficient naming an undeclared material, or a zero direction
    ("h4_undeclared_material",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Eigen(CellCoefficient({(7, F)}))),
     ValueError, r"material 7 is not declared"),
    ("h4_k_on_a_non_producing_specification",
     lambda: Specification(materials=Materials({0: moderator()}), geometry=None, question=Eigen(_k())),
     ValueError, r"zero direction"),
    # (h5) two point keys naming one coordinate
    ("h5_one_coordinate_twice",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None,
                           question=FixedSource(_table(1), {_k(): 0.1, CellCoefficient({(0, F)}): 0.2})),
     ValueError, r"name the same coordinate"),
    # a specification exists to be keyed: a geometry with no content is refused eagerly
    ("contentless_geometry",
     lambda: Specification(
         materials=Materials({0: fuel()}),
         geometry=StructuredGeometry.slab((0.0, 1.0), (0,), left=PrescribedInflow(lambda space: np.zeros(space.shape)),  # pyright: ignore[reportArgumentType]  # a contentless source is the subject
                                          right=BC.vacuum),
         question=Eigen(_k())),
     ContentlessError, r"Specification\.geometry"),
    # the three fields are typed
    ("question_not_a_question",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=_table(1)),  # pyright: ignore[reportArgumentType]
     TypeError, r"the question is an Eigen, a FixedSource or a Response"),
    ("geometry_not_a_geometry",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=(0.0, 1.0), question=Eigen(_k())),  # pyright: ignore[reportArgumentType]
     TypeError, r"the geometry is a StructuredGeometry or None"),
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


def test_s8_1_the_infinite_medium_rule_is_read_before_the_group_count() -> None:
    """No geometry, two materials of different group counts: (a) fires, (c) does
    not. The group count is read over the materials the specification KEEPS,
    and with no geometry that set is defined only once (a) has held."""
    with pytest.raises(ValueError, match="one material") as err:
        Specification(materials=Materials({0: fuel(), 1: fuel(ng=3)}), geometry=None, question=Eigen(_k()))
    require(not isinstance(err.value, InconsistentMaterialsError), f"(c) fired with (a): {err.value}")


def test_s8_1_a_spectator_of_another_group_count_is_dropped() -> None:
    """A declared material the geometry does not assign is dropped before the
    group count is read (the leak principle: a spectator changes no answer)."""
    spec = Specification(materials=Materials({0: fuel(), 5: fuel(ng=3)}), geometry=slab2((0, 0)), question=Eigen(_k()))
    require(sorted(spec.materials.ids) == [0], f"the spectator was kept: {sorted(spec.materials.ids)}")


def test_s8_1_the_composition_is_read_before_the_keys() -> None:
    """The geometry names an undeclared id AND the key names an undeclared
    material: (b) fires, (h4) does not (keys resolve against a valid declaration)."""
    with pytest.raises(ValueError, match=r"references material ids \[1\]") as err:
        Specification(materials=Materials({0: fuel()}), geometry=slab2(), question=Eigen(CellCoefficient({(7, F)})))
    require("material 7" not in str(err.value), f"(h4) fired with (b): {err.value}")


# ── (d)-(f): the refused states of the old table are unspellable ─────────────


@pytest.mark.rests_on(f"{_Q}::test_s7_3_the_fields_are_exactly_the_roles")
def test_s8_1_the_fields_are_exactly_materials_geometry_question() -> None:
    """(d)-(f) are structural: the datum lives in the question (``FixedSource.source``,
    ``Response.detector``), so the specification has no ``source`` field to hold a
    second copy (ruling of 2026-10-02 on the step-7 NEEDS). Witness: a ``source``
    field re-added reds this row."""
    names = tuple(f.name for f in dataclasses.fields(Specification))
    require(names == ("materials", "geometry", "question"), f"the fields are {names}")


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
    ("k on a fissile slab", lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Eigen(_k()))),
    ("k in the infinite medium", lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(_k()))),
    ("c, every emission cell", lambda: Specification(materials=_two_materials(), geometry=slab2(),
                                                      question=Eigen(CellCoefficient.every(F, S, N2)))),
    ("a critical extent", lambda: Specification(materials=_two_materials(), geometry=slab2(),
                                                question=Eigen(GeometryExtent(1), mode=Nearest(1.0)))),
    ("a spectator material", lambda: Specification(materials=Materials({**_two_materials(), 9: moderator()}), geometry=slab2(),
                                                   question=Eigen(_k()))),
    ("a fixed source on a non-producing slab, offset along scattering",
     lambda: Specification(materials=Materials({0: moderator(), 1: moderator()}), geometry=slab2(),
                           question=FixedSource(_table(), {CellCoefficient.every(S): 0.1}))),
    ("a response offset along an extent",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Response(_table(), {GeometryExtent(0): -0.2}))),
    ("the parameter's own offset (step-7 ruling 3)",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(_k(), {_k(): 0.3}))),
    ("an explicit zero cell beside a non-zero one (ruling 5)",
     lambda: Specification(materials=_two_materials(), geometry=slab2(), question=Eigen(CellCoefficient({(0, F), (1, F)})))),
    ("a constant Symbolic in the infinite medium",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=FixedSource(Symbolic.of(2, 1)))),
    ("one region in the infinite medium",
     lambda: Specification(materials=Materials({0: fuel()}), geometry=None, question=Response(_table(1)))),
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
    build = lambda: Specification(materials=_two_materials(), geometry=body(coord), question=question)  # noqa: E731
    if admitted:
        build()
    else:
        with pytest.raises(ValueError, match=rf"the {role} depends on the azimuth phi"):
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
    build = lambda: Specification(materials=Materials({0: _CARRIERS[carried]()}), geometry=None,  # noqa: E731
                                  question=Eigen(CellCoefficient.every(asked)))
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


def _cells_of(key: Any) -> frozenset:
    return frozenset(key.cells) if isinstance(key, CellCoefficient) else frozenset()


@pytest.mark.rests_on(f"{_CELLS}::test_s8_7_resolution_by_channel")
def test_s8_10_every_resolves_to_the_explicit_non_zero_cells() -> None:
    """``every(F)`` on {0: fuel, 1: moderator} is stored as {(0, F)}: the
    quantifier is gone and the moderator's zero cell is not a cell of the key."""
    q = Eigen(_k(), {CellCoefficient.every(S): 0.25})
    spec = Specification(materials=_two_materials(), geometry=slab2(), question=q)
    require(_parameter(spec) == CellCoefficient({(0, F)}), f"parameter {_parameter(spec)}")
    require(dict(spec.question.point) == {CellCoefficient({(0, S), (1, S)}): 0.25}, f"point {spec.question.point}")
    stored = _cells_of(_parameter(spec)) | frozenset().union(*(_cells_of(k) for k in spec.question.point))
    require(all(isinstance(m, int) for m, _ in stored), f"a quantifier survived in {stored}")
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
    a = Specification(materials=_two_materials(), geometry=slab2(), question=Eigen(_k()))
    b = Specification(materials=_two_materials(), geometry=slab2(), question=written())
    require(a == b and a.content_digest == b.content_digest and len({a, b}) == 1, "two specifications")
    require(_parameter(b) == CellCoefficient({(0, F)}), f"stored {_parameter(b)}")


def test_s8_10_canonicalisation_is_idempotent_and_survives_pickle() -> None:
    spec = Specification(materials=_two_materials(), geometry=slab2(),
                         question=FixedSource(_table(), {CellCoefficient.every(S, N2): -0.5, GeometryExtent(1): 0.125}))
    again = Specification(materials=spec.materials, geometry=spec.geometry, question=spec.question)
    require(again == spec and again.question == spec.question, "re-posing the stored question moved it")
    back = pickle.loads(pickle.dumps(spec))
    require(back == spec and back.content_digest == spec.content_digest, "pickle moved the specification")
    require([v for v in spec.question.point.values()] == [-0.5, 0.125], "the offsets moved")


def test_s8_10_the_materials_are_the_declaration_value() -> None:
    """One spelling: the materials are a ``Materials`` value, and a bare mapping
    is refused (the orchestrator's ruling 2 on the step-8 NEEDS, 2026-10-02)."""
    with pytest.raises(TypeError, match="the materials are a Materials, got a dict"):
        Specification(materials={0: fuel(), 1: moderator()}, geometry=slab2(), question=Eigen(_k()))  # pyright: ignore[reportArgumentType]  # the refusal is the subject


