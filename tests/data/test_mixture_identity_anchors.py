r"""Consumers campaign step 1 — the PRE-CARVE anchors for ``Mixture``'s identity.

Landed on the UNMODIFIED tree, before the first production edit of the
identity step (plan ``.claude/plans/cs4c_binding_design.md`` §27, rulings
**R-cc3** / **R-cc8** / **R-cc9**; GitHub **#459**), so every row here is a
measurement of the tree the carve starts from.

``Mixture`` is the step's third fix and the LEAF of the whole identity
chain: the SN hub's materials entry, ``HomogeneousProblem``'s one field and
``Materials``' entries are all mixtures, so ONE content key serves all
three (§27's "one definition, three fixes").

**Two claim kinds live here and they have OPPOSITE fates — read the class
docstring before touching a row.**

* :class:`TestTodaysMixtureEqualityIsNotContentEquality` — **RECORD** of a
  state the carve DELETES. Green today, designed to RED at the carve; the
  carve's commit deletes the class. Its value is that it makes the API
  change LOUD: without it, the ``xfail`` rows below could silently stay
  ``xfail`` if the new ``__eq__`` lands wrong (``vv`` Mode 8, the
  misattributed-strict-xfail class).
* :class:`TestContentIdentity` and :class:`TestMutabilityIsTheHashHazard`
  — the RULED post-carve gates, shipped as ``xfail(strict=True)`` so the
  marker set is a self-retiring todo list. Their XPASS is a FAILURE, which
  is what forces the marker's deletion in the carve's own commit.

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

pytestmark = pytest.mark.foundation

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
# RECORD — today's state, which the carve DELETES
# ════════════════════════════════════════════════════════════════════


class TestTodaysMixtureEqualityIsNotContentEquality:
    """RECORD of the defect #459 removes. Designed to RED at the carve.

    Delete this class in the commit that gives ``Mixture`` a content
    ``__eq__``/``__hash__``; do not "repair" it. Its whole job is to make
    the API change impossible to land quietly — the ``xfail`` rows below
    are silent if the new equality lands wrong, and these are not.
    """

    def test_equal_data_RAISES_at_two_groups_and_above(self) -> None:
        r"""``[M]`` 2026-09-12: ``ValueError: truth value of an array …`` at ng = 2 and 4.

        The generated ``__eq__`` compares the field TUPLE; the first ndarray
        field reaches ``bool(array)`` and raises. This is the defect: a
        comparison that cannot even RETURN, hidden behind a 1-group arm that
        accidentally works.
        """
        for groups, ng in (("2g", 2), ("4g", 4)):
            a, b = _fresh(groups=groups), _fresh(groups=groups)
            _require(a.ng == ng, f"activation: the {groups} fixture must have ng == {ng}")
            with pytest.raises(ValueError, match="truth value of an array"):
                a == b  # noqa: B015  # pyright: ignore[reportUnusedExpression]

    def test_the_one_group_arm_is_the_trap_not_the_contract(self) -> None:
        r"""``[M]`` ng = 1 returns ``True`` — accidentally, not by design.

        Every field array holds ONE element, so ``bool(array)`` succeeds and
        the dataclass default behaves like a content comparison. A gate
        written here would be green before AND after #459 and would prove
        nothing (``vv`` #19 — the reading that cannot change).
        """
        a, b = _fresh(groups="1g"), _fresh(groups="1g")
        _require(a.ng == 1, "activation: this row is ABOUT the degenerate group count")
        _require(a is not b, "the two mixtures must be distinct objects")
        _require(a == b, "today's ng = 1 arm returns True for equal data")
        _require(not (a == _fresh("C", "1g")), "and False for different data")

    def test_mixture_is_unhashable(self) -> None:
        r"""``[M]`` ``TypeError: unhashable type: 'Mixture'`` at every group count.

        ``@dataclass`` with the default ``eq=True`` sets ``__hash__ = None``.
        So a ``Mixture`` cannot be a dict key, a set member, or part of any
        other object's hash — which is why ``HomogeneousProblem.__hash__``
        currently hashes ``id(self.mixture)``.
        """
        for groups in ("1g", "2g", "4g"):
            with pytest.raises(TypeError, match="unhashable type"):
                hash(_fresh(groups=groups))

    def test_mixture_is_mutable_today(self) -> None:
        r"""``[M]`` a ``Mixture`` accepts post-construction assignment.

        Not a defect on its own — but it is the reason a CACHED content key
        is unsound without a freeze, and the reason
        :class:`TestMutabilityIsTheHashHazard` exists. ``[M]`` 33 sites in
        14 test files assign ``mix.SigS`` / ``.Sig2`` / ``.SigT`` after
        construction; 0 in ``orpheus/``.
        """
        mix = _fresh()
        before = float(mix.SigT[0])
        mix.SigT = np.asarray(mix.SigT, dtype=float).copy()
        mix.SigT[0] = before + 1.0
        _require(float(mix.SigT[0]) == before + 1.0, "Mixture is not frozen today")


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

    @pytest.mark.xfail(strict=True, reason=_RULING)
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

    @pytest.mark.xfail(strict=True, reason=_RULING)
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

    @pytest.mark.xfail(strict=True, reason=_RULING)
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

    @pytest.mark.xfail(strict=True, reason=_RULING)
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

    @pytest.mark.xfail(strict=True, reason=f"{_RULING}; open ruling O-4 (freeze vs live key)")
    def test_the_hash_cannot_go_stale(self) -> None:
        mix = _fresh(groups="2g")
        before = hash(mix)
        try:
            moved = np.asarray(mix.SigT, dtype=float).copy()
            moved[0] = moved[0] + 0.125
            mix.SigT = moved
        except Exception:
            return  # frozen — the other honest outcome
        _require(
            hash(mix) != before,
            "a cached content hash over a mutable Mixture is stale after a field write; "
            "freeze the datum (O-4) or derive the hash live",
        )
