r"""Consumers campaign step 1 — the PRE-CARVE anchors for ``Mixture``'s identity.

Landed on the UNMODIFIED tree (`31b632c3`), before the first production edit of the
identity step (plan ``.claude/plans/cs4c_binding_design.md`` §27, rulings
**R-cc3** / **R-cc8** / **R-cc9**; GitHub **#459**), so every row here is a
measurement of the tree the carve starts from.

``Mixture`` is the step's third fix and the LEAF of the whole identity
chain: the SN hub's materials entry, ``HomogeneousProblem``'s one field and
``Materials``' entries are all mixtures, so ONE content key serves all
three (§27's "one definition, three fixes").

**Two claim kinds live here and they have OPPOSITE fates — read the class
docstring before touching a row.**

* ``TestTodaysMixtureEqualityIsNotContentEquality`` — the **RECORD** of the
  state the carve deleted (the raising default equality, the ng = 1 false
  green, unhashability, mutability). It was green on the pre-carve tree
  and DELETED in the S1a commit with its subject, as its docstring asked.
* :class:`TestContentIdentity` and :class:`TestMutabilityIsTheHashHazard`
  — the RULED post-carve gates, shipped PRE-carve as ``xfail(strict=True)``
  (a self-retiring todo list) and turned GREEN in S1a, which deleted the
  markers; O-4 was ruled FROZEN, so the hazard row asserts the freeze.

⚠ **The ng = 1 trap, measured.** ``[M]`` ``Mixture(ng=1) == Mixture(ng=1)``
returns ``True`` TODAY — every field array holds one element, so
``bool(array)`` succeeds and the dataclass default is accidentally a
content comparison. At ng ≥ 2 the identical expression raises
``ValueError``. So a gate written on a 1-group fixture proves NOTHING about
the defect, and every row below that exercises equality asserts its own
group count (``vv`` #11 activation, ``lessons`` L67h).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture

pytestmark = [pytest.mark.foundation, pytest.mark.catches("ERR-084")]

#: The ruling this module's xfail rows retire with.
_RULING = "#459 / R-cc3 — ONE content-identity definition; Mixture is its leaf"


def _require(cond: object, msg: str) -> None:
    """``-O``-safe assertion (``coding-standards``: the canonical runner strips ``assert``)."""
    if not cond:
        raise AssertionError(msg)


def _fresh(region: str = "A", groups: str = "2g") -> Mixture:
    """A mixture built from the library — a NEW object on every call.

    ``[M]`` ``get_mixture`` re-derives; two calls are two objects carrying
    equal data, which is exactly the pair the content tier must equate and
    the constituent tier must not.
    """
    return get_mixture(region, groups)


def _perturb(mix: Mixture, field: str) -> Mixture:
    """Return a copy of ``mix`` with exactly ONE generating datum moved.

    The per-datum negative legs (``vv`` #11): every field of the dataclass
    gets its own flip, so a key that silently drops one is a red with a
    name rather than a green with a gap. The flip is the smallest one that
    survives ``__post_init__`` — ``chi`` is re-normalised by the emission
    law, so it is moved by PERMUTATION (mass-preserving) rather than by a
    scale, and the sparse stacks are moved in their ``.data``.
    """
    if field == "chi":
        chi = np.asarray(mix.chi, dtype=float).copy()
        if chi.size < 2:
            pytest.skip("the chi permutation leg needs ng >= 2")
        chi[[0, 1]] = chi[[1, 0]]
        return replace(mix, chi=chi)
    if field in ("SigS", "Sig2"):
        stack = [b.copy() for b in getattr(mix, field)]
        blk = stack[0].tolil()
        blk[0, 0] = blk[0, 0] + 0.125
        stack[0] = csr_matrix(blk)
        return replace(mix, **{field: stack})
    if field == "SigS_high_order":
        stack = [b.copy() for b in mix.SigS]
        if len(stack) < 2:
            pytest.skip("the high-order leg needs a P1 block in the library")
        blk = stack[1].tolil()
        blk[0, 0] = blk[0, 0] + 0.125
        stack[1] = csr_matrix(blk)
        return replace(mix, SigS=stack)
    if field == "eg":
        base = np.asarray(mix.eg) if mix.eg is not None else None
        if base is None:
            base = np.geomspace(1.0e-3, 2.0e7, mix.ng + 1)[::-1].copy()
            return replace(mix, eg=base)
        moved = base.copy()
        moved[0] = moved[0] * 1.5
        return replace(mix, eg=moved)
    arr = np.asarray(getattr(mix, field), dtype=float).copy()
    arr[0] = arr[0] + 0.125
    return replace(mix, **{field: arr})


#: Every generating datum of a ``Mixture``, one flip each. ``SigS_high_order``
#: is the eleventh row and is NOT a field: it moves the P1 block alone, which a
#: key that hashes only ``stack[0]`` would miss (``[M]`` the whole abstract
#: library ships ``len(SigS) == 2``, so the block EXISTS on every row here).
_PERTURBABLE = (
    "SigC", "SigL", "SigF", "SigP", "SigT",
    "chi", "eg", "SigS", "SigS_high_order", "Sig2",
)


# ════════════════════════════════════════════════════════════════════
# The RULED post-carve gates — strict xfail, a self-retiring todo list
# ════════════════════════════════════════════════════════════════════


class TestContentIdentity:
    """R-cc3: equal generating data ⟹ the same mixture, hashable.

    Every row is ``xfail(strict=True)``: it FAILS today (the comparison
    raises) and its XPASS is a failure, so the carve's own commit must
    delete the marker. Each row's body is shaped so exactly ONE statement
    can fail and it is the documented one (``vv`` Mode 8, fourth class).
    """

    @pytest.mark.parametrize("groups,ng", [("2g", 2), ("4g", 4)])
    def test_equal_data_compares_equal_and_hashes_equal(self, groups: str, ng: int) -> None:
        """The POSITIVE control — two independent builds of one mixture."""
        a, b = _fresh(groups=groups), _fresh(groups=groups)
        _require(a.ng == ng and a is not b, "activation: two distinct ng >= 2 objects")
        _require(a == b, "equal generating data must compare equal")
        _require(hash(a) == hash(b), "equal mixtures must hash equal")
        # The reportUnhashable diagnostic below IS the pre-carve state this
        # row exists to flip: a Mixture becomes hashable at #459.
        pair = {a, b}  # pyright: ignore[reportUnhashable]
        _require(len(pair) == 1, "a set of two equal mixtures holds one member")

    @pytest.mark.parametrize("field", _PERTURBABLE)
    def test_one_moved_datum_reads_unequal(self, field: str) -> None:
        """The per-datum NEGATIVE legs — one flip per generating datum.

        A key that drops any single field would pass the positive control
        and fail exactly one of these rows, which is the point of writing
        them one-per-field instead of as a loop inside one row (``vv``
        anti-#17's granularity clause).
        """
        base = _fresh(groups="2g")
        moved = _perturb(base, field)
        _require(base.ng == 2, "activation: the raise-arm group count")
        _require(not (base == moved), f"moving {field} must make the mixtures unequal")

    def test_the_p1_block_is_in_the_key(self) -> None:
        """``SigS[1]`` alone moved — the row a ``stack[0]``-only key fails.

        ``[M]`` every shipped ``xs_library`` mixture has ``len(SigS) == 2``,
        so the P1 block exists on every fixture in this file and the row is
        never vacuous. Asserted in-test.
        """
        base = _fresh(groups="2g")
        _require(len(base.SigS) >= 2, "activation: the library must ship a P1 block")
        _require(
            not (base == _perturb(base, "SigS_high_order")),
            "a key that hashes only the P0 block cannot see a P1 change",
        )

    def test_the_energy_grid_is_in_the_key(self) -> None:
        """``eg`` — the field ``[M]`` ``None`` on all 12 abstract fixtures.

        ``lessons`` L59b: a corpus uniform in a discriminating field leaves
        the arm witness-less. The witness here is MANUFACTURED (a grid is
        attached by ``replace``), exactly as the repo's only ``eg``-bearing
        homogeneous mixture is built.
        """
        base = _fresh(groups="2g")
        _require(base.eg is None, "activation: the abstract library ships eg = None")
        with_grid = _perturb(base, "eg")
        _require(with_grid.eg is not None, "the manufactured witness must carry a grid")
        _require(not (base == with_grid), "attaching an energy grid changes the mixture")


class TestCrossClassComparisonIsAlreadySafe:
    """MUST STAY GREEN — a property the carve must PRESERVE, not create.

    ⛔ I shipped this as an ``xfail`` row and the harness refuted it:
    ``[M]`` 2026-09-12 it ``XPASS(strict)``-ed on the unmodified tree. The
    generated ``__eq__`` opens with ``if other.__class__ is self.__class__``
    and returns ``NotImplemented`` otherwise, so a foreign comparison never
    reaches an ndarray field and never raises.

    ⟹ this is not a gap the carve closes; it is a REGRESSION PIN on the
    hand-written ``__eq__`` that replaces the generated one. The
    ``Axis`` precedent spells the same discipline explicitly
    (``numerics/axis.py:266`` — ``other.__class__ is not self.__class__``
    ⟹ ``NotImplemented``), and the reason it is load-bearing here is
    ``Solution``'s generated ``__eq__``: it reaches the mesh, hence the
    materials, hence every mixture, so a raising heterogeneous comparison
    would surface at a call site that never mentions ``Mixture``.
    """

    def test_comparison_across_classes_never_raises(self) -> None:
        a = _fresh(groups="2g")
        _require(a.ng == 2, "activation: the group count where equal-data == RAISES")
        _require(not (a == object()), "a mixture is not an arbitrary object")
        _require(not (a == None), "a mixture is not None")  # noqa: E711
        _require(a != object(), "the negated form must agree")


class TestMutabilityIsTheHashHazard:
    """The CACHED-key hazard, stated as a gate (open ruling O-4).

    ``[M]`` building the content key of a 421-group mixture costs **1022 µs**
    (3 214 496 bytes over 7 dense arrays and 14 CSR triples), while hashing a
    key already built costs **0.57 µs**. ``_GEOM_CACHE_INTERN`` reads the
    hash on every operator construction (``[M]`` 6 slab / 10 sphere per
    eigenvalue solve), so the key MUST be computed once — and a cached key
    over a MUTABLE datum is a silent lie the first time a test mutates one.

    The row asserts the outcome of whichever way O-4 is ruled: either
    ``Mixture`` is frozen (mutation raises) or the hash tracks mutation
    (no cache). It cannot be satisfied by a cached key on a mutable object,
    which is the state that must not ship.
    """

    def test_the_hash_cannot_go_stale(self) -> None:
        """O-4 RULED (2026-09-12): FROZEN. A field write and an in-place array
        write are both refused, so the cached key can never be stale."""
        import dataclasses

        mix = _fresh(groups="2g")
        before = hash(mix)
        moved = np.asarray(mix.SigT, dtype=float).copy()
        moved[0] = moved[0] + 0.125
        with pytest.raises(dataclasses.FrozenInstanceError):
            mix.SigT = moved  # type: ignore[misc]
        with pytest.raises(ValueError, match="read-only"):
            mix.SigT[0] = moved[0]
        with pytest.raises(ValueError, match="read-only"):
            mix.SigS[0].data[0] = 999.0
        _require(hash(mix) == before, "the key moved without the data moving")
        # the honest way to a variant: replace re-runs the laws and mints a NEW value
        twin = replace(mix, SigT=moved)
        _require(twin != mix and hash(twin) != before, "replace must mint a different value")
