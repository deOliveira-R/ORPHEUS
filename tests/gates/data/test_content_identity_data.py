r"""Content identity of the data layer: ``Mixture`` and ``Materials`` (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5, gates S5.2 (equal content
gives equal ``==``, ``hash`` and digest), S5.3 (a change to any part changes
the digest, quantified over the parts the production instance reports) and
the data legs of S5.4 (signed zero canonicalised, NaN refused). The encoder
and the mixin are gated in ``tests/gates/numerics/test_content_identity.py``.

First red, measured on ``main`` ``1dc31163`` (spec §1.5, "First reds"): the
module fails to import (``orpheus.numerics.content`` does not exist). With the
encoder alone added to that tree, the rows below red on the legacy behaviour
each one rejects: ``Materials`` twins unequal (``eq=False``); ``Materials``
unpicklable (it holds a ``mappingproxy``); a ``Mixture`` whose ``SigL`` is
``-0.0`` unequal to the ``+0.0`` one; a ``Mixture`` with an explicitly stored
CSR zero unequal to the one without; a NaN cross section admitted.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.numerics.content import content_digest
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
    require,
)
from tests.gates.data.test_mixture_identity_anchors import _perturb

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/data/test_content_identity_data.py"
_ENCODER = "tests/gates/numerics/test_content_identity.py"
_S51 = f"{_ENCODER}::test_s5_1_digests_and_hashes_are_seed_stable"
_S54 = f"{_ENCODER}::TestS54EncoderCanonicalForms"


def _mix(region: str = "A", groups: str = "2g") -> Mixture:
    """A fresh library mixture: two calls are two objects of equal content."""
    return get_mixture(region, groups)


def _with_stored_zero(mix: Mixture) -> Mixture:
    """``mix`` with an explicit 0.0 stored at an unstored position of ``SigS[0]``.

    The same matrix: ``[M]`` mixture A 2g stores 3 of 4 entries of its P0 block
    and position (1, 0) is unstored; this builds a 4-entry CSR for it. Today
    ``Mixture`` canonicalises duplicates and index order but not stored zeros.
    """
    p0 = mix.SigS[0].tocoo()
    dense = p0.toarray()
    unstored = np.argwhere(dense == 0.0)
    require(len(unstored) >= 1, "activation: the P0 block needs an unstored position")
    r, c = (int(i) for i in unstored[0])
    rows = np.append(p0.row, r)
    cols = np.append(p0.col, c)
    vals = np.append(p0.data, 0.0)
    block = csr_matrix((vals, (rows, cols)), shape=p0.shape)
    require(block.nnz == p0.nnz + 1, "activation: the stored zero must be stored")
    return replace(mix, SigS=(block, *mix.SigS[1:]))


def _with_int64_indices(mix: Mixture) -> Mixture:
    p0 = mix.SigS[0]
    block = csr_matrix(
        (p0.data.copy(), p0.indices.astype(np.int64), p0.indptr.astype(np.int64)),
        shape=p0.shape,
    )
    return replace(mix, SigS=(block, *mix.SigS[1:]))


def _p1_moved(mix: Mixture) -> Mixture:
    return _perturb(mix, "SigS_high_order")


def _numpy_keyed() -> "dict[Any, Mixture]":
    """The declaration with numpy integer keys: the input the coercion parses."""
    return {np.int64(0): _mix("A"), np.int32(3): _mix("B")}


def _materials() -> Materials:
    return Materials({0: _mix("A"), 3: _mix("B")})


# ── The roster: the population S5.3 quantifies over ─────────────────────────

_MIXTURE = Entry(
    cls=Mixture,
    base=_mix,
    parts=("SigC", "SigL", "SigF", "SigP", "SigT", "SigS", "Sig2", "chi", "eg"),
    perturb={
        **{
            name: (leg("+0.125 on entry 0", lambda name=name: _perturb(_mix(), name)),)
            for name in ("SigC", "SigL", "SigF", "SigP", "SigT", "chi", "eg", "Sig2")
        },
        "SigS": (
            leg("P0 block", lambda: _perturb(_mix(), "SigS")),
            leg("P1 block alone", lambda: _p1_moved(_mix())),
        ),
    },
    pairs=(
        ("two library builds", _mix, _mix),
        ("CSR int64 indices", _mix, lambda: _with_int64_indices(_mix())),
        ("explicit stored CSR zero", _mix, lambda: _with_stored_zero(_mix())),
    ),
)

_MATERIALS = Entry(
    cls=Materials,
    base=_materials,
    parts=("mixtures",),
    fields_are_parts=False,
    perturb={
        "mixtures": (
            leg("a key relabelled", lambda: Materials({0: _mix("A"), 4: _mix("B")})),
            leg("a key added", lambda: Materials({0: _mix("A"), 3: _mix("B"), 5: _mix("C")})),
            leg("a value moved", lambda: Materials({0: _perturb(_mix("A"), "SigT"), 3: _mix("B")})),
            leg("the values swapped between keys", lambda: Materials({0: _mix("B"), 3: _mix("A")})),
        ),
    },
    pairs=(
        ("two builds", _materials, _materials),
        ("key insertion order", _materials, lambda: Materials({3: _mix("B"), 0: _mix("A")})),
        ("numpy integer keys", _materials, lambda: Materials(_numpy_keyed())),
        ("Materials.of a dict", _materials, lambda: Materials.of({0: _mix("A"), 3: _mix("B")})),
    ),
)

ROSTER: tuple[Entry, ...] = (_MIXTURE, _MATERIALS)


# ── S5.3: the population, then one row per perturbation ──────────────────────


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_3_population(entry: Entry) -> None:
    """S5.3's denominator: the parts are the type's own (``content_parts()``
    or the ``compare=True`` fields), and the table perturbs every one."""
    check_population(entry)


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s5_3_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
def test_s5_3_float32_rounding_is_a_value_change() -> None:
    """A float32-rounded ``SigT`` is a different value, so a different digest.

    The dtype is storage (float64 at construction), the VALUE is content.
    ``[M]`` mixture A 4g's ``SigT`` is not float32-representable, so the
    rounding moves it (asserted: the row is never vacuous).
    """
    base = _mix("A", "4g")
    rounded = replace(base, SigT=base.SigT.astype(np.float32))
    require(not np.array_equal(base.SigT, rounded.SigT), "activation: rounding must move SigT")
    require(rounded.SigT.dtype == np.float64, "the stored dtype is float64 either way")
    require(content_digest(rounded) != content_digest(base), "a rounded SigT must move the digest")
    require(rounded != base, "a rounded SigT is a different mixture")


# ── S5.2: equal content is one value ─────────────────────────────────────────


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s5_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    """First red with the encoder present: ``Materials`` twins unequal
    (``eq=False``); the stored CSR zero unequal (``[M]`` today nnz 4 vs 3)."""
    check_equal_pair(entry, pair)


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_2_pickle_round_trip(entry: Entry) -> None:
    """First red: ``pickle.dumps(Materials(...))`` raises ``TypeError: cannot
    pickle 'mappingproxy' object`` (``[M]`` 1dc31163)."""
    check_pickle(entry)


@pytest.mark.rests_on(_S51)
def test_s5_2_materials_keys_are_ints_and_order_is_not_content() -> None:
    """The keys are coerced to ``int``; the declared order is kept (iteration
    order is behaviour, which ``Materials.restrict`` documents: "preserves the
    given id order") and is not content (a reordered twin is equal, hashes
    alike). ``[M]`` before step 5 an ``np.int64`` key stayed ``np.int64``.

    Corrected 2026-10-02 by the main agent: this row first asserted sorted
    iteration (the plan's "key-sorted tuple"), which would break
    ``restrict``'s contract (``tests/gates/data/test_materials.py::TestRestrict``)."""
    declared: dict[Any, Mixture] = {np.int64(3): _mix("B"), 0: _mix("A")}
    mats = Materials(declared)
    keys = list(mats)
    require(keys == [3, 0], f"iteration order {keys} is not the declared order")
    require(all(type(k) is int for k in keys), f"key types {[type(k) for k in keys]}")
    require(3 in mats and np.int64(3) in mats, "membership by value")
    reordered = Materials({0: _mix("A"), 3: _mix("B")})
    require(mats == reordered and hash(mats) == hash(reordered), "id order is not content")


@pytest.mark.rests_on(_S51)
def test_s5_2_a_material_mesh_pickles() -> None:
    """The persisted value the reference cache keys on: ``MaterialMesh`` holds a
    ``Materials``, so it could not be pickled either (``[M]`` 1dc31163)."""
    import pickle

    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    geometry = StructuredGeometry(
        coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 1.0, 2.0), mat_ids=(0, 3),
        boundaries=(BC.reflective, BC.vacuum),
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh
    mm = MaterialMesh(mesh, _materials())
    twin = pickle.loads(pickle.dumps(mm))
    require(type(twin) is type(mm), "the reloaded type")
    require(twin.materials == mm.materials, "the reloaded declaration is equal")


# ── S5.4, the data legs: signed zero is one value; NaN is not a value ──────


@pytest.mark.rests_on(_S54)
def test_s5_4_signed_zero_is_one_cross_section() -> None:
    """``[M]`` today ``SigL = -0.0`` and ``SigL = +0.0`` give UNEQUAL mixtures
    (equality by bytes) while ``BC`` treats them as equal: the ruling (2026-09-25,
    NEEDS 1) is one convention, ``-0.0`` canonicalised."""
    plus = replace(_mix(), SigL=np.zeros(2))
    minus = replace(_mix(), SigL=-np.zeros(2))
    require(bool(np.signbit(-np.zeros(2)).all()), "activation: the input carries -0.0")
    require(content_digest(plus) == content_digest(minus), "signed zero must digest alike")
    require(plus == minus, "signed zero must compare equal")
    require(hash(plus) == hash(minus), "signed zero must hash equal")


def _nan_variant(name: str) -> "dict[str, Any]":
    """The ``replace`` kwargs putting one NaN into field ``name``."""
    mix = _mix()
    if name in ("SigS", "Sig2"):
        stack = [b.copy() for b in getattr(mix, name)]
        block = stack[0].tolil()
        block[0, 0] = np.nan
        stack[0] = csr_matrix(block)
        return {name: stack}
    if name == "eg":
        return {"eg": np.array([2.0e7, np.nan, 1.0e-5])}
    values = np.array(getattr(mix, name), dtype=float)
    values[0] = np.nan
    return {name: values}


@pytest.mark.rests_on(_S54)
@pytest.mark.parametrize("name", ["SigC", "SigL", "SigF", "SigP", "SigT", "chi", "eg", "SigS", "Sig2"])
def test_s5_4_a_nan_mixture_cannot_be_constructed(name: str) -> None:
    """Parse at the boundary (main's ruling, 2026-10-02): a NaN in any
    ``Mixture`` array (the five dense fields, ``chi``, ``eg``, a sparse block's
    data) is refused at CONSTRUCTION, a ``ValueError`` naming the field, so
    ``==`` and ``hash`` never raise on a constructible mixture. ``[M]`` at
    1dc31163 a NaN ``SigT`` constructs. No digest is taken: a refusal only at
    the encoder (its backstop) leaves this row red."""
    kwargs = _nan_variant(name)
    with pytest.raises(ValueError, match=name):
        replace(_mix(), **kwargs)
