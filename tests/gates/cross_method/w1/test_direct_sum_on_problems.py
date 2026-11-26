r"""G05a (runtime half), G05c, G05d, G05f (carrier half), G05g, G05j, G05l, G05o (trace-leaf arms), G05q and G05r: the direct sum on the harness problems (N05, spec 6.6 and 7.6).

The rows that build a harness problem live here (the N04 rule, 7.5); the rows that build only ``numerics`` objects live in ``tests/gates/numerics/w1/test_direct_sum.py``. The population is the 33
composites of the 15 harness fixtures (21 full fields: 12 SN, 6 System-B ray instances, 3 diffusion; 12 coupled spaces: 6 with one summand, 6 with two, the second a ray instance), reached through the
problems' public getters by the population builder of the composite record, which names neither class the unit retires.

Each name that exists only after the unit is read at call time, so on the base tip a row fails in its body. Where a row is a regression gate that today's classes already satisfy by construction (G05g's
``x.space == V`` checks, G05l's two charts) its first red is the mutation, and the row says so; its calls to the new surfaces (``carrier``, ``blocks``, ``with_blocks``, ``label=``) are red at the base tip.

Host exactness (CI: ubuntu with numpy 2.5.5; this host: macOS arm64 with Accelerate and numpy 2.4.4): every row is exact on any host, or compares a computation with a closed form within a stated
tolerance. The two charts are compared byte for byte within one process (two computations through one code path on one host; the comparison with stored bytes is G05m's and is host-keyed). The payload chart
is compared with the column-by-column flat chart within a relative 1e-13 of the largest entry (a dense block applies a matrix-matrix product where the flat chart applies a matrix-vector one, which differ in
the last ulps and not in the structure). Mutations are written for the reviewers and run on a copy; none is run by this file. Every row of this file is listed in the CI step "The space algebra gates (W1)".
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _composite_probe_population as pop
from tests.gates.cross_method.w1 import _fixtures as fx
from tests.gates.cross_method.w1 import _snapshot as snap
from tests.gates.numerics.w1 import _direct_sum_probe as dp
from tests.gates.numerics.w1 import _tensor_probe as tp

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
PAYLOAD_RTOL = 1e-13


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def sn_fixtures() -> list[tuple[str, object]]:
    return [(snap.sn_fixture_name(kind, ng), fx.sn_problem(kind, ng)[0]) for kind in fx.SN_KINDS for ng in fx.GROUPS]


def diffusion_fixtures() -> list[tuple[str, object]]:
    return [(snap.diffusion_fixture_name(ng), fx.diffusion_problem(ng)[0]) for ng in fx.GROUPS]


def composites() -> list[tuple[str, object]]:
    """The 33 composite spaces of the harness, each named ``<fixture>.<getter>``."""
    out: list[tuple[str, object]] = []
    for kind in fx.SN_KINDS:
        for ng in fx.GROUPS:
            out += list(pop._composites_of_sn(kind, ng).items())
    for name, mesh in diffusion_fixtures():
        out.append((f"{name}.full_field_space", mesh.full_field_space))
    return out


def composites_with_ng() -> list[tuple[str, object, int]]:
    """The 33 composites with the group count of their fixture."""
    out: list[tuple[str, object, int]] = []
    for kind in fx.SN_KINDS:
        for ng in fx.GROUPS:
            out += [(name, space, ng) for name, space in pop._composites_of_sn(kind, ng).items()]
    for ng in fx.GROUPS:
        out.append((f"{snap.diffusion_fixture_name(ng)}.full_field_space", fx.diffusion_problem(ng)[0].full_field_space, ng))
    return out


def test_the_population_is_33_composites() -> None:
    """The denominator of every row below (X2): 33 composites, 21 full fields (12 SN, 6 ray instances, 3 diffusion) and 12 coupled spaces. Exact on any host."""
    found = composites()
    kinds = [name.split(".", 1)[1] for name, _ in found]
    _require(len(found) == 33, f"the population holds {len(found)} composites, expected 33")
    _require(kinds.count("full_field_space") == 15 and kinds.count("radial_field_space") == 6 and kinds.count("system.space") == 12, f"the split of the population is {[(k, kinds.count(k)) for k in sorted(set(kinds))]}")


# ──────────────────────────────── G05a, the runtime half ────────────────────────────────


def test_the_production_sums_have_the_declared_shape() -> None:
    """G05a, the runtime half. The production full field is ``DirectSumSpace.of_pairs(((BULK, bulk), (TRACE, trace)))``, System B's ray instance is the same call with ``interior`` and ``boundary``, and the
    coupled system is a sum of the full field, nested. 15 SN and diffusion full fields, 6 ray instances, 12 coupled spaces. First red: import-only. Mutation: build the production full field as a flat
    sum of its two summands (the nested arm reddens)."""
    from orpheus.numerics import direct_sum as m

    for name, problem in sn_fixtures():
        field = problem.full_field_space
        _require(type(field) is m.DirectSumSpace and field.labels == (m.BULK, m.TRACE), f"{name}: the full field is {type(field).__name__} with labels {getattr(field, 'labels', None)}")
        coupled = problem.system.space
        _require(type(coupled) is m.DirectSumSpace and coupled.summand(coupled.labels[0]) == field, f"{name}: the coupled space does not nest the full field as its first summand")
        ray = problem.radial_characteristic_field_space
        if ray is not None:
            _require(type(ray) is m.DirectSumSpace and ray.labels == ("interior", "boundary"), f"{name}: the ray instance has labels {getattr(ray, 'labels', None)}")
            _require(len(coupled.labels) == 2 and coupled.summand(coupled.labels[1]) == ray, f"{name}: the coupled space does not nest the ray instance as its second summand")
        else:
            _require(len(coupled.labels) == 1, f"{name}: the coupled space has {len(coupled.labels)} summands without a ray instance")
    for name, mesh in diffusion_fixtures():
        _require(mesh.full_field_space.labels == (m.BULK, m.TRACE), f"{name}: the diffusion full field has labels {mesh.full_field_space.labels}")


# ──────────────────────────────── G05c ────────────────────────────────


def test_flip_addressed_to_a_label() -> None:
    """G05c, L6, two populations. (c1) tensor-only sums: a label held by one summand flips only that summand, and ``dual(on=L)`` likewise; a label held by no summand is a ``KeyError`` naming the labels; a
    flip whose addressed slot is formless in a summand that holds it (``riesz(on={"energy"})`` on a sum whose summand has a formless energy slot) is refused whole (``ValueError``) and nothing is returned;
    a summand label passed as ``on=`` is a ``KeyError`` (two namespaces). (c2) a mixed sum, a tensor summand and a legacy leaf: every non-empty ``on=`` flip and a whole ``riesz()`` is refused
    ("a legacy leaf holds no factor labels", SD-22 widened); whole ``dual()`` is the sum of the summands' duals. The production full field is (tensor + leaf) in 15 of 21 and (leaf + leaf) in the other
    6, and 12 of 12 coupled spaces are mixed, so (c2) runs on the 33. First red: import-only (a fixture with a whole ``sk`` on a singular diagonal form is unconstructible after N03, so the refusal
    the law needs is built from the formless summand). Mutations: (c1) apply the flip to the summands that accept it and skip the refusing one (the whole-refusal arm reddens); (c2) the pre-widening rule,
    flip only the tensor summand (the refusal arm reddens)."""
    S = dp.direct_sum()
    a, c, p = dp.space_of(dp.A_ARGS), dp.space_of(dp.C_ARGS), dp.positive_space()
    d = S.of(t=a, p=p)
    flipped = d.dual(on={"p"})
    _require(flipped.summand("t") == a and flipped.summand("p") == p.dual(on={"p"}) and flipped.summand("p") != p, "a label held by one summand did not flip only that summand")
    both = S.of(a=a, c=c).dual(on={"angular"})
    _require(both.summand("a") == a.dual(on={"angular"}) and both.summand("c") == c.dual(on={"angular"}), "a label held by two summands did not flip both")
    with pytest.raises(KeyError, match="nope"):
        d.dual(on={"nope"})
    with pytest.raises(KeyError):
        d.dual(on={"t"})  # a summand label is not a factor label (two namespaces)
    with pytest.raises(ValueError, match="riesz is not defined"):
        d.riesz(on={"energy"})
    _require(d.dual().summand("t") == a.dual() and d.dual().summand("p") == p.dual(), "whole dual is not the sum of the summands' duals")
    seen = 0
    for name, space in composites():
        seen += 1
        with pytest.raises((ValueError, TypeError), match="legacy leaf holds no factor labels"):
            space.dual(on={"angular"})
        with pytest.raises((ValueError, TypeError), match="legacy leaf holds no factor labels"):
            space.riesz()
        whole = space.dual()
        for label in space.labels:
            _require(whole.summand(label) == space.summand(label).dual(), f"{name}: the dual's summand {label!r} is not the dual of the summand")
    _require(seen == 33, f"c2 ran on {seen} composites, expected 33")


# ──────────────────────────────── G05d ────────────────────────────────


def _parts(space, rng: np.random.Generator, payload: tuple[int, ...]) -> dict[str, np.ndarray]:
    return {label: rng.standard_normal(space.summand(label).shape + payload) for label in space.labels}


def test_pack_unpack_and_slices() -> None:
    """G05d. (i) ``slice_of(label)`` slices cover ``range(n)`` exactly once in label order (parsed: each flat index in exactly one slice); (ii) ``unpack(pack(parts)) == parts`` for random parts with and
    without a trailing payload (byte-equal); (iii) ``pack(unpack(c)) == c``; (iv) on the 33 composites ``pack`` of the summand-wise ``apply_metric`` equals ``apply_metric`` of the flat chart byte for
    byte (two computations on one host); (v) ``ast`` over ``orpheus/``: no class named ``SpacePartition`` and no ``partition_by`` defined in a class derived from ``FunctionSpace`` (ON05-23; green today by
    construction, its first red is the mutation; ``DiscreteMeasure.partition_by`` is another object and is the positive control of the predicate). First red: import-only for (i) to (iv). Mutations: shift
    one slice's start by one (the cover check); swap two entries in ``unpack`` (the round trip); add ``partition_by`` to the base ``FunctionSpace`` (arm (v))."""
    S = dp.direct_sum()
    a, b, c = (dp.space_of(x) for x in (dp.A_ARGS, dp.B_ARGS, dp.C_ARGS))
    d = S.of(a=a, b=b, c=c)
    cover = np.zeros(d.n_coordinates, dtype=int)
    previous_stop = 0
    for label in d.labels:
        sl = d.slice_of(label)
        _require(sl.start == previous_stop, f"slice_of({label!r}) starts at {sl.start}, the previous slice stopped at {previous_stop}")
        cover[sl] += 1
        previous_stop = sl.stop
    _require(bool(np.all(cover == 1)), f"(i) the slices do not cover range({d.n_coordinates}) exactly once: {np.bincount(cover).tolist()}")
    rng = np.random.default_rng(31)
    for payload in ((), (2,)):
        parts = _parts(d, rng, payload)
        flat = d.pack(parts)
        _require(flat.shape == (d.n_coordinates,) + payload, f"(ii) pack returned shape {flat.shape} for payload {payload}")
        back = d.unpack(flat)
        _require(sorted(back) == sorted(parts) and all(np.asarray(back[k]).tobytes() == parts[k].tobytes() and np.asarray(back[k]).shape == parts[k].shape for k in parts), f"(ii) unpack(pack(parts)) differs from parts for payload {payload}")
        c0 = rng.standard_normal((d.n_coordinates,) + payload)
        _require(d.pack(d.unpack(c0)).tobytes() == c0.tobytes(), f"(iii) pack(unpack(c)) differs from c for payload {payload}")
    seen = 0
    for name, space in composites():
        v = snap._probe_vector(int(space.n_coordinates))
        parts = space.unpack(v)
        out = {}
        for label in space.labels:
            summand = space.summand(label)
            out[label] = snap._apply(summand, "apply_metric", np.asarray(parts[label]).reshape(-1)).reshape(summand.shape)
        got = space.pack(out)
        want = np.ascontiguousarray(np.asarray(space.apply_metric(v))).reshape(-1)
        _require(got.tobytes() == want.tobytes(), f"(iv) {name}: pack of the summand-wise apply_metric differs from apply_metric of the flat chart")
        seen += 1
    _require(seen == 33, f"(iv) ran on {seen} composites, expected 33")
    planted = "class FunctionSpace:\n    def partition_by(self, p):\n        return p\n"
    _require(partition_by_in_derived([("p.py", planted)]) == {"FunctionSpace": ["p.py"]}, "(v) the predicate does not read a planted partition_by on FunctionSpace")
    sources = dp.orpheus_sources()
    _require(partition_by_in_classes(sources).get("DiscreteMeasure"), "(v) positive control: DiscreteMeasure.partition_by is not read")
    found = partition_by_in_derived(sources)
    _require(not found, f"(v) partition_by is defined in classes derived from FunctionSpace: {found}")
    named = [p for p, s in sources if "SpacePartition" in dp.class_defs(s)]
    _require(not named, f"(v) a class named SpacePartition is defined in {named}")


def partition_by_in_classes(sources: list[tuple[str, str]]) -> dict[str, list[str]]:
    """``{class name: [files]}`` of every class that defines a method ``partition_by``."""
    out: dict[str, list[str]] = {}
    for path, source in sources:
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.ClassDef) and any(isinstance(m, ast.FunctionDef) and m.name == "partition_by" for m in node.body):
                out.setdefault(node.name, []).append(path)
    return out


def partition_by_in_derived(sources: list[tuple[str, str]]) -> dict[str, list[str]]:
    """The classes of ``partition_by_in_classes`` that are ``FunctionSpace`` or derive from it by base names within ``sources``."""
    bases: dict[str, set[str]] = {}
    for _path, source in sources:
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.ClassDef):
                bases.setdefault(node.name, set()).update(b.id if isinstance(b, ast.Name) else b.attr if isinstance(b, ast.Attribute) else "" for b in node.bases)
    derived = {"FunctionSpace"}
    changed = True
    while changed:
        changed = False
        for name, parents in bases.items():
            if name not in derived and parents & derived:
                derived.add(name)
                changed = True
    return {k: v for k, v in partition_by_in_classes(sources).items() if k in derived}


# ──────────────────────────────── G05f, the carrier half ────────────────────────────────


def test_carrier_pair_is_swapped_by_the_dual() -> None:
    """G05f, the carrier half, on the production coupled spaces (the pair class is the build's, so the row reads the pair through ``carrier``). The dual's carrier is the swapped pair:
    ``D.dual().carrier.zeros is D.carrier.dual_zeros`` and ``D.dual().carrier.dual_zeros is D.carrier.zeros``; the dual of a full field has no zero factory (21 of 21 duals: ``dual().carrier.zeros()`` raises
    ``RuntimeError``, the MB4 mechanism, which stays a strict xfail in ``test_member_seam_builders.py``, W1-P6). First red: import-only (the sum has no ``carrier`` today). Mutation: do not swap the pair
    in ``dual`` (the identity arms redden)."""
    seen = 0
    for name, problem in sn_fixtures():
        space = problem.system.space
        dual = space.dual()
        _require(dual.carrier.zeros is space.carrier.dual_zeros and dual.carrier.dual_zeros is space.carrier.zeros, f"{name}: the dual's carrier is not the swapped pair")
        _require(space.carrier.zeros is not None and space.carrier.dual_zeros is not None, f"{name}: the coupled space's pair is not installed (zeros and dual_zeros are wired at the 5 coupled sites)")
        seen += 1
    for name, space in composites():
        if name.endswith("full_field_space") or name.endswith("radial_field_space"):
            with pytest.raises((RuntimeError, TypeError, AttributeError)):  # the pair was never installed: its use raises (the exception type is the build's)
                space.dual().carrier.zeros()
    _require(seen == 12, f"the carrier half ran on {seen} coupled spaces, expected 12")


# ──────────────────────────────── G05g ────────────────────────────────


def test_carrier_space_equals_builder_space() -> None:
    """G05g. On the SN slab, box, sphere and cylinder and the diffusion slab (15 fixtures): the coupled states returned by ``_returned_state`` and ``_as_coupled`` (``sn/solver.py:710``, ``:717``) have
    ``.space == problem.system.space``; ``V.member(c)``, ``V.carrier.zeros()``, ``with_blocks`` and the restriction build ``CoupledField`` whose labels equal the space's (the formal paths cannot mismatch,
    ON05-21); ``problem.full_field_space.carrier.zeros().space == problem.full_field_space``; the coupled state ``problem.system.space.carrier.zeros().space == problem.system.space``; the restricted coupled
    field ``SystemRestrictionOperator(V, label=...).apply(z).space == r.codomain`` (addressed by label, ON05-24); the System-B field ``radial_characteristic_field_space.carrier.zeros().space ==`` the ray
    instance. Today's classes satisfy the equalities by construction (57 of 57 checks equal on the base tip: identity carries no label), so a regression row: its first red is the mutation, and its calls
    to ``carrier`` and ``label=`` are red on the base tip. Mutations: rename the label of one class constant (``FullField``'s trace label) or of the builder's coupled labels (the equality reddens on the
    affected fixtures); swap the two builder constants at ``sn/solver.py``'s two sites (the ``_returned_state`` arm reddens)."""
    from orpheus.numerics.coupled_system import SystemRestrictionOperator
    from orpheus.sn.solver import _as_coupled, _returned_state

    checks = 0
    for name, problem in sn_fixtures():
        field_space, V, ray_space = problem.full_field_space, problem.system.space, problem.radial_characteristic_field_space
        _require(field_space.carrier.zeros().space == field_space, f"{name}: the full field's zero does not live on the full field space")
        _require(V.carrier.zeros().space == V, f"{name}: the coupled zero does not live on the coupled space")
        member = V.member(snap._probe_vector(int(V.n_coordinates)))
        _require(member.space == V and member.labels == V.labels, f"{name}: V.member(c) has labels {member.labels}, expected {V.labels}")
        z = V.carrier.zeros()
        _require(z.labels == V.labels and z.with_blocks(z.blocks).labels == V.labels and z.with_blocks(z.blocks).space == V, f"{name}: zeros or with_blocks lost the labels")
        ray = None if ray_space is None else ray_space.carrier.zeros()
        _require(ray is None or ray.space == ray_space, f"{name}: the ray field's zero does not live on the ray instance")
        state = _returned_state(field_space.carrier.zeros(), ray)
        _require(state.space == V and _as_coupled(state).space == V, f"{name}: _returned_state or _as_coupled returns a state on another space than the hub's")
        if len(V.labels) == 1:
            _require(_as_coupled(field_space.carrier.zeros()).space == V, f"{name}: _as_coupled of a bare composite is not on the one-summand coupled space")
        for label in V.labels:
            restriction = SystemRestrictionOperator(V, label=label)
            restricted = restriction.apply(z)
            _require(restricted.space == restriction.codomain and restricted.labels == (label,), f"{name}: the restriction onto {label!r} returns a field on another space or with labels {restricted.labels}")
        checks += 1
    for name, mesh in diffusion_fixtures():
        _require(mesh.full_field_space.carrier.zeros().space == mesh.full_field_space, f"{name}: the diffusion full field's zero does not live on its space")
        checks += 1
    _require(checks == 15, f"G05g ran on {checks} fixtures, expected 15")


# ──────────────────────────────── G05j ────────────────────────────────


def _energy_summand(groups: int):
    from orpheus.numerics.axis import EnergyAxis

    axis = EnergyAxis.of_edges(np.linspace(1.0, 2.0, groups + 1))
    return tp.tensor_space().of_axes(tp.formless_factor(axis))


def test_energy_extent_of_a_sum() -> None:
    """G05j. A sum of tensor summands, one with an ``EnergyAxis`` of 2 groups, bound with data ``ng = 3``: refused with the fragment "energy extent"; with ``ng = 2`` it passes; a summand with no energy
    axis is skipped; a nested sum recurses; the 6 System-B ray instances of the harness pass at every ``ng``, because both their summands are legacy leaves with no energy axis until N06 (the row asserts that
    they pass, so that N06's removal of the edge is a red it must repair); the 27 tensor summands of the harness that hold an energy axis are refused at ``ng + 1`` and pass at the true ``ng``. The walk's
    leaf skip is tagged ``ELEGANCE-DEBT[guard] #759 legacy-leaf`` (pinned in G04o's tagged set by the unit's commit). First red: 33 of 33 harness composites pass unchecked today, at a wrong ``ng`` too
    (the guard returns on a sum), so the refusal arms are behavioural reds at the base tip (positive control: the bulk summand alone with ``ng = 3`` is refused today). Mutation: make the guard return on a
    sum again (today's behaviour)."""
    from orpheus.transport.operators._energy_conformity import assert_energy_extent_conforms as guard

    S = dp.direct_sum()
    energetic, plain = _energy_summand(2), dp.space_of(dp.A_ARGS)
    flat = S.of(e=energetic, p=plain)
    nested = S.of(n=S.of(e=energetic, p=plain), q=dp.positive_space())
    for name, space in (("flat", flat), ("nested", nested)):
        with pytest.raises(ValueError, match="energy extent"):
            guard(space, 3, operator="G05j")
        guard(space, 2, operator="G05j")
        del name
    guard(S.of(p=plain, q=dp.positive_space()), 7, operator="G05j")  # no energy axis anywhere: skipped
    refused = passed = 0
    for name, space, ng in composites_with_ng():
        guard(space, ng, operator="G05j")  # the true ng always passes
        if name.endswith("radial_field_space"):
            guard(space, ng + 1, operator="G05j")  # the edge: both summands are legacy leaves
            guard(space, ng + 7, operator="G05j")
            passed += 1
            continue
        with pytest.raises(ValueError, match="energy extent"):
            guard(space, ng + 1, operator="G05j")
        refused += 1
    _require((refused, passed) == (27, 6), f"the guard refused {refused} and passed {passed} of the 33 harness composites at ng + 1, expected 27 and 6")


# ──────────────────────────────── G05l ────────────────────────────────


def _class_names(blocks) -> list[str]:
    return [type(b).__name__ for b in blocks]


def test_the_two_charts_agree() -> None:
    """G05l. ``apply_metric`` and ``apply_inverse_metric`` of the sum on the typed carrier and on the flat coordinates (with and without a trailing payload) agree: ``to_flat(apply_metric(x)) ==
    apply_metric(to_flat(x))`` byte for byte on the 33 harness composites (the coupled spaces are the production nested shape); the payload chart equals the flat chart column by column within a relative
    1e-13 of the largest entry; the typed result is the carrier's own class, keeps the blocks' classes and ``.space``; ``inner_product`` and ``norm`` take the carrier only (flat arrays raise
    ``AttributeError``: 33 of 33 today); the sum reads and rebuilds its typed element through the carrier's ``blocks`` and ``with_blocks`` and, for a block whose summand is not a sum, through ONE helper
    that calls ``dataclasses.replace(block, values=...)`` (ON05-16): the ``ast`` arms are ``test_the_sum_acts_through_its_summands``. Today's classes accept both charts (33 of 33), so a regression row: its
    first red is the mutation; the calls to ``blocks`` are red at the base tip. Mutations: return the typed result as the base class; apply the metric to the flat chart in summand order reversed."""
    for name, space in composites():
        n = int(space.n_coordinates)
        v = snap._probe_vector(n)
        x = space.member(v)
        payload = dp_payload(n)
        for verb in ("apply_metric", "apply_inverse_metric"):
            typed = getattr(space, verb)(x)
            _require(type(typed) is type(x) and typed.space == x.space and _class_names(typed.blocks) == _class_names(x.blocks), f"{name}.{verb}: the typed result is {type(typed).__name__} with blocks {_class_names(typed.blocks)}")
            flat = np.ascontiguousarray(np.asarray(getattr(space, verb)(v))).reshape(-1)
            _require(np.ascontiguousarray(typed.to_flat()).reshape(-1).tobytes() == flat.tobytes(), f"{name}.{verb}: the typed chart differs from the flat chart in the bytes")
            both = np.asarray(getattr(space, verb)(payload))
            _require(both.shape == payload.shape, f"{name}.{verb}: the payload chart has shape {both.shape}, expected {payload.shape}")
            for k in range(payload.shape[1]):
                column = np.asarray(getattr(space, verb)(np.ascontiguousarray(payload[:, k]))).reshape(-1)
                scale = max(float(np.max(np.abs(column))), 1e-300)
                _require(float(np.max(np.abs(both[:, k] - column))) <= PAYLOAD_RTOL * scale, f"{name}.{verb}: payload column {k} differs from the flat chart beyond {PAYLOAD_RTOL:g}")
        with pytest.raises(AttributeError):
            space.inner_product(v, v)
        with pytest.raises(AttributeError):
            space.norm(v)


def dp_payload(n: int) -> np.ndarray:
    return pop._payload(n)


def _methods_of(source: str, cls_name: str) -> dict[str, ast.FunctionDef]:
    tree = ast.parse(source)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == cls_name)
    return {m.name: m for m in cls.body if isinstance(m, ast.FunctionDef)}


def forbidden_calls(function: ast.AST, names: frozenset[str]) -> list[int]:
    return [n.lineno for n in ast.walk(function) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr in names]


def replace_helpers(source: str) -> list[tuple[str, int]]:
    """``(function, line)`` of every call of ``replace`` (a name or an attribute) in a module, with the enclosing function."""
    tree = ast.parse(source)
    out: list[tuple[str, int]] = []
    for fn in [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
        for n in ast.walk(fn):
            if isinstance(n, ast.Call) and ((isinstance(n.func, ast.Name) and n.func.id == "replace") or (isinstance(n.func, ast.Attribute) and n.func.attr == "replace")):
                out.append((fn.name, n.lineno))
    return sorted(set(out))


def test_the_sum_acts_through_its_summands() -> None:
    """G05l, the ``ast`` arms over ``numerics/direct_sum.py``. In class ``DirectSumSpace``: no call of ``.apply``, ``.apply_inverse``, ``._walk`` or ``.apply_block`` (G04l reads class ``TensorSpace`` only, so the
    sum's delegation is read here), and each metric verb (``apply_metric``, ``apply_inverse_metric``, ``inner_product``, ``metric_diagonal``, ``rank``, ``kernel_basis``, ``negative``) reads the same
    member of a summand. The sum reads and rebuilds its typed element through ``blocks`` and ``with_blocks`` and, for a leaf block, through ONE helper that calls ``dataclasses.replace(block, values=...)``
    (ON05-16): no read of ``interior``, ``boundary`` or ``systems`` of a carrier, no call of ``_recombine``, no ``replace(`` with a ``systems`` keyword and no ``replace(`` outside that helper (scoped to the
    space's own access). The predicates are controlled by planted sources (the old accesses ``_rebuild`` of ``FullFieldSpace`` and ``_member_wise`` of ``CoupledSpace`` are the base tip's reads, retired
    with them). First red: import-only (``direct_sum.py`` does not exist). Mutation: rebuild through ``replace(x, systems=...)`` inside ``direct_sum.py``; call ``summand.apply`` in ``apply_metric``."""
    planted = "from dataclasses import replace\nclass DirectSumSpace:\n    def apply_metric(self, x):\n        y = x.interior\n        return replace(x, systems=(y,))\n    def rank(self):\n        return self._recombine(1).apply(2)\n"
    methods = _methods_of(planted, "DirectSumSpace")
    _require(forbidden_calls(methods["rank"], frozenset({"apply"})) != [], "the call predicate does not read a planted .apply")
    _require(replace_helpers(planted) == [("apply_metric", 5)], "the replace predicate does not read a planted replace")
    source = (ROOT / "orpheus/numerics/direct_sum.py").read_text()
    methods = _methods_of(source, "DirectSumSpace")
    for name, fn in methods.items():
        calls = forbidden_calls(fn, frozenset({"apply", "apply_inverse", "_walk", "apply_block", "_recombine"}))
        _require(not calls, f"DirectSumSpace.{name} calls a retired or forbidden verb at lines {calls}")
    for verb in ("apply_metric", "apply_inverse_metric", "inner_product", "metric_diagonal", "rank", "kernel_basis", "negative"):
        _require(verb in methods, f"DirectSumSpace defines no {verb}")
        reads = {n.attr for n in ast.walk(methods[verb]) if isinstance(n, ast.Attribute)}
        _require(verb in reads, f"DirectSumSpace.{verb} does not read {verb} of a summand (reads: {sorted(reads)})")
    tree = ast.parse(source)
    reads = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Attribute) and n.attr in {"interior", "boundary", "systems"}]
    _require(not reads, f"direct_sum.py reads interior, boundary or systems of a carrier at lines {reads}")
    helpers = replace_helpers(source)
    _require(len({name for name, _ in helpers}) == 1, f"replace( is spelled in {len({n for n, _ in helpers})} functions of direct_sum.py: {helpers}")
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and ((isinstance(n.func, ast.Name) and n.func.id == "replace") or (isinstance(n.func, ast.Attribute) and n.func.attr == "replace")):
            keywords = {k.arg for k in n.keywords}
            _require("systems" not in keywords and keywords == {"values"}, f"a replace( call of direct_sum.py has keywords {sorted(k for k in keywords if k)}, expected exactly values (line {n.lineno})")


# ──────────────────────────────── G05o, the trace-leaf arms ────────────────────────────────


def test_cone_conjunction_with_a_trace_leaf() -> None:
    """G05o, the arms with a legacy trace leaf (the numerics-only arms are in ``test_direct_sum.py``). A trace leaf answers None (precondition, control); on a sum with a trace leaf and nodal tensor summands
    ``has_coordinate_cone`` is None; on (modal head + trace leaf) it is False; on a sum nesting that sum beside a nodal tensor summand it is False (the order False before None, ON05-10); on the 33 harness
    composites it stays None (every production sum has a trace leaf that answers None until N06); on the 6 lifted full fields (moment bulk) it is False. First red at the base tip: the composite classes answer
    None on the modal-head sum. Mutation: evaluate None before False (the (modal head + trace leaf) and the nested arms redden; the older arms stay green, ``[M]`` check F4)."""
    S = dp.direct_sum()
    leaf = diffusion_fixtures()[0][1].scalar_trace
    head, _ = dp.modal_head_and_tail()
    nodal = dp.space_of(dp.A_ARGS)
    _require(leaf.has_coordinate_cone is None and head.has_coordinate_cone is False and nodal.has_coordinate_cone is True, "the controls of the cone arms do not read None, False and True")
    _require(S.of(t=leaf, a=nodal).has_coordinate_cone is None, "a trace leaf beside a nodal tensor summand is not None")
    mixed = S.of(h=head, t=leaf)
    _require(mixed.has_coordinate_cone is False, f"(modal head + trace leaf) reads {mixed.has_coordinate_cone!r}, expected False")
    _require(S.of(n=mixed, a=nodal).has_coordinate_cone is False, "the nested sum with a modal head reads other than False")
    for name, space in composites():
        _require(space.has_coordinate_cone is None, f"{name}: has_coordinate_cone is {space.has_coordinate_cone!r}, expected None")
    lifted = [pop._lifted_full_field(fx.sn_problem(kind, ng)[0]) for kind in ("slab", "box") for ng in fx.GROUPS]
    _require(len(lifted) == 6 and all(s.has_coordinate_cone is False for s in lifted), f"the lifted full fields read {[s.has_coordinate_cone for s in lifted]}, expected False on 6 of 6")


# ──────────────────────────────── G05q ────────────────────────────────


def test_typed_element_protocol() -> None:
    """G05q (ON05-1, ON05-16 to ON05-18). For every carrier class with a witness built for each (``FullField``, ``RadialCharacteristicField`` and ``CoupledField`` on the 33 harness composites, 0
    ``TimedFullField`` among them; a ``TimedFullField`` with a non-empty history is built in the gate): (q1) ``D.coordinates_of(D.member(c))`` is byte-equal to ``c`` for seeded ``c``; (q2) ``y =
    D.member(D.coordinates_of(x))`` has ``type(y) is type(x)``, ``y.space == x.space``, the classes of ``y.blocks`` equal those of ``x.blocks``, and ``D.coordinates_of(y)`` is byte-equal to
    ``D.coordinates_of(x)``; (q3) put-get: for blocks ``b'`` of the right classes ``x.with_blocks(b').blocks`` are ``b'`` element by element (``is``), and get-put: ``x.with_blocks(x.blocks).blocks`` are
    ``x.blocks`` (``is``), with the class kept; (q4) on a ``TimedFullField`` with a history, ``with_blocks`` keeps ``history_depth`` and empties the history, as ``_recombine`` does today. The carriers
    compare by identity (two mints of the same coordinates are unequal on 33 of 33), so no law is stated by ``==`` on a carrier. A timed field is not a member the space mints (``V.member`` returns a
    ``FullField``), so (q2) is not claimed of it; it satisfies (q3), (q4) and ``V.coordinates_of(t)`` equal to ``t.to_flat()``. First red: import-only (no carrier defines ``with_blocks`` today: ``blocks``
    there is a property of ``CoupledOperator``). Mutations: ``with_blocks`` returns ``self`` (q3); permute the blocks in ``coordinates_of`` (q1); keep the history (q4)."""
    classes: set[str] = set()
    for name, space in composites():
        n = int(space.n_coordinates)
        c = np.random.default_rng(101).standard_normal(n)
        x = space.member(c)
        classes.add(type(x).__name__)
        _require(np.ascontiguousarray(space.coordinates_of(space.member(c))).reshape(-1).tobytes() == c.tobytes(), f"(q1) {name}: coordinates_of(member(c)) differs from c")
        y = space.member(space.coordinates_of(x))
        _require(type(y) is type(x) and y.space == x.space and _class_names(y.blocks) == _class_names(x.blocks), f"(q2) {name}: the round trip returns {type(y).__name__} with blocks {_class_names(y.blocks)}")
        _require(np.ascontiguousarray(space.coordinates_of(y)).reshape(-1).tobytes() == np.ascontiguousarray(space.coordinates_of(x)).reshape(-1).tobytes(), f"(q2) {name}: the coordinates of the round trip differ")
        other = space.member(np.random.default_rng(102).standard_normal(n)).blocks
        put = x.with_blocks(other)
        _require(type(put) is type(x) and len(put.blocks) == len(other) and all(p is o for p, o in zip(put.blocks, other)), f"(q3) {name}: with_blocks does not put the blocks given")
        get_put = x.with_blocks(x.blocks)
        _require(type(get_put) is type(x) and all(p is o for p, o in zip(get_put.blocks, x.blocks)), f"(q3) {name}: with_blocks(x.blocks) does not return x's blocks")
    _require(classes == {"FullField", "RadialCharacteristicField", "CoupledField"}, f"the carrier classes of the 33 composites are {sorted(classes)}")
    problem = fx.sn_problem("slab", 2)[0]
    space, timed = pop._timed_pair(problem)
    _require(timed.history_length == 1 and timed.history_depth == 2, "the timed witness does not carry one frame in its history")
    rebuilt = timed.with_blocks(timed.blocks)
    _require(type(rebuilt) is type(timed) and rebuilt.history_depth == 2 and rebuilt.history_length == 0, f"(q4) with_blocks on a timed field gives depth {rebuilt.history_depth} and history {rebuilt.history_length}, expected 2 and 0")
    _require(all(p is o for p, o in zip(rebuilt.blocks, timed.blocks)), "(q3) the timed field's with_blocks does not keep the blocks")
    _require(np.ascontiguousarray(space.coordinates_of(timed)).reshape(-1).tobytes() == np.ascontiguousarray(timed.to_flat()).reshape(-1).tobytes(), "V.coordinates_of(t) differs from t.to_flat()")


# ──────────────────────────────── G05r ────────────────────────────────

COMPOSITE_CLASSES = (("orpheus/transport/full_field.py", "FullField"), ("orpheus/transport/timed_full_field.py", "TimedFullField"), ("orpheus/transport/radial_characteristic_field.py", "RadialCharacteristicField"))


def class_labels_constants(source: str, cls_name: str) -> list[ast.expr]:
    """The values assigned in the class body to a name containing ``labels`` (case-insensitive): the class constant of the label pair."""
    cls = next(n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.ClassDef) and n.name == cls_name)
    out = []
    for m in cls.body:
        if isinstance(m, ast.Assign) and any(isinstance(t, ast.Name) and "labels" in t.id.lower() for t in m.targets):
            out.append(m.value)
        if isinstance(m, ast.AnnAssign) and isinstance(m.target, ast.Name) and "labels" in m.target.id.lower() and m.value is not None:
            out.append(m.value)
    return out


def test_labels_are_class_constants_and_required() -> None:
    """G05r (ON05-2, ON05-20, ON05-21). On the 21 full fields of the harness the sum's labels are ``("bulk", "trace")`` for the 15 tensor-bulk instances and ``("interior", "boundary")`` for the 6 ray
    instances; the label pair is a class-level constant of each concrete ``Composite`` class (``FullField``, ``TimedFullField``: ``bulk`` and ``trace``; ``RadialCharacteristicField``: ``interior`` and
    ``boundary``), spelled from the builder's constants imported (names, never string literals), and no ``Composite`` class has a ``labels`` field; ``CoupledField(systems=...)`` without labels raises
    ``TypeError`` and no ``CoupledField`` field has a default label; ``bulk_and_trace`` of a ray instance raises. First red at the base tip: the 17 production ``CoupledField(`` constructors pass no
    label and the class accepts them (behavioural). Mutations: give ``CoupledField`` a default label tuple; set ``RadialCharacteristicField``'s constant to ``bulk`` and ``trace`` (the instance and
    the refusal arms redden)."""
    from orpheus.numerics import direct_sum as m
    from orpheus.numerics.coupled_system import CoupledField

    tensor_bulk = ray = 0
    for name, space in composites():
        if name.endswith("full_field_space"):
            _require(space.labels == (m.BULK, m.TRACE) == ("bulk", "trace"), f"{name}: labels {space.labels}")
            tensor_bulk += 1
        elif name.endswith("radial_field_space"):
            _require(space.labels == ("interior", "boundary"), f"{name}: labels {space.labels}")
            with pytest.raises(TypeError, match="composite DirectSumSpace"):
                m.bulk_and_trace(space, owner="G05r")
            ray += 1
    _require(tensor_bulk == 15 and ray == 6, f"the full-field population is {tensor_bulk} and {ray}, expected 15 and 6")
    for rel, cls_name in COMPOSITE_CLASSES:
        values = class_labels_constants((ROOT / rel).read_text(), cls_name)
        _require(len(values) == (0 if cls_name == "TimedFullField" else 1), f"{cls_name} has {len(values)} class-level label constants in its own body (TimedFullField inherits FullField's)")
        literals = [n for v in values for n in ast.walk(v) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
        _require(not literals, f"{cls_name}'s label constant spells a string literal {[l.value for l in literals]}: the builder's constants are imported (ON05-20)")
    planted = "class FullField:\n    LABELS = ('bulk', 'trace')\n"
    _require(len([n for n in ast.walk(class_labels_constants(planted, "FullField")[0]) if isinstance(n, ast.Constant)]) == 2, "the literal predicate does not read planted literals")
    from orpheus.transport.full_field import FullField
    from orpheus.transport.radial_characteristic_field import RadialCharacteristicField
    from orpheus.transport.timed_full_field import TimedFullField

    for cls in (FullField, TimedFullField, RadialCharacteristicField):
        _require("labels" not in {f.name for f in dataclasses.fields(cls)}, f"{cls.__name__} has a labels field")
    constants = {cls.__name__: [getattr(cls, k) for k in dir(cls) if "labels" in k.lower() and isinstance(getattr(cls, k), tuple)] for cls in (FullField, TimedFullField, RadialCharacteristicField)}
    _require(constants == {"FullField": [("bulk", "trace")], "TimedFullField": [("bulk", "trace")], "RadialCharacteristicField": [("interior", "boundary")]}, f"the class constants are {constants}")
    labels_field = {f.name: f for f in dataclasses.fields(CoupledField)}.get("labels")
    _require(labels_field is not None and labels_field.default is dataclasses.MISSING and labels_field.default_factory is dataclasses.MISSING, "CoupledField.labels is missing or has a default")
    z = fx.sn_problem("slab", 1)[0].full_field_space.carrier.zeros()
    with pytest.raises(TypeError):
        CoupledField(systems=(z,))  # type: ignore[call-arg]
