r"""#434 R3 — the registry's ledger, and the coverage test that reads it: the gates.

Designed by the test-architect against the LANDED carve (2026-09-03; every ``[M]``
is a measurement on the landed code — the carve landed mid-dispatch, so the
pre-carve readings come from a ``git archive HEAD`` shadow tree, never a checkout).
Landed as ONE module rather than split across ``test_symmetry.py`` /
``test_registry.py`` / ``test_manifold.py`` because the classes share the
element-set construction helpers and the frozen stage-0 grid; each class states
the production home of the claim it pins.

* Groups A (the coverage predicate's laws — its two conjuncts separately, the
  finiteness refusal, a GENUINE product neither factor contains, the inverse
  direction, agreement with an independent element-set construction) are
  closed-form pure math; the structurally independent ground is the element set.
* Group B is the frozen pre-carve stage-0 grid
  (``tests/numerics/data/r3_stage_baseline.json``) with the RULED exceptions.
* Group C measures the ledger's ``unspent`` rows on the SOLVER — a fold licensed
  by the ledger must not change the flux; one the ledger refuses does.
* The fold consequence (the empty sweep quadrants ERR-081 catalogues), the
  ledger's finiteness guard, the retired names, and the stage-0 reason pinned
  where it is decided.

**Markers.**  Every row is a SOFTWARE INVARIANT or a closed-form law with no theory
``:label:``; ``foundation``, no ``verifies``.  The battery that proves these rows
bite is ``scratch/_r3_mut.py`` (17 arms; the two γ⁻¹ / product-order arms can bite
only on the genuine-product rows here, every shipped Γ being an involution group).
"""
from __future__ import annotations

import dataclasses
import json
import pathlib

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.manifold import SPHERE, Quotient, quotient_onto
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import (
    GEOMETRY_ANGULAR_SYMMETRY,
    AngularSymmetry,
    select_quadrature,
)
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.sn.solver import solve_sn_fixed_source
from tests.numerics import NUMERICS_DATA

ORPHEUS_ROOT = pathlib.Path(__file__).resolve().parents[2]

S = SubgroupOfO3
_I3 = np.eye(3)
_E = {"x": np.array([1.0, 0.0, 0.0]),
      "y": np.array([0.0, 1.0, 0.0]),
      "z": np.array([0.0, 0.0, 1.0])}
_AXIS = {"x": 0, "y": 1, "z": 2}


# ===========================================================================
# The INDEPENDENT construction — plain numpy, from the DEFINITIONS.
#
# `lessons` L72c, restated for this predicate: the reference SHARES the standard
# setting (a 3x3 orthogonal matrix in the same basis — "is this matrix in that
# group" is literal membership IN THAT SETTING and cannot be made independent of
# it).  What IS independent is the ALGORITHM: production decides `H ⊆ ΓK` by
# (Lie-algebra containment of the identity components) ∧ (a per-representative
# coset search); the reference decides it by ENUMERATING Γ·K's element set from
# each group's definition — a BFS closure for the finite factors, and for the one
# continuous factor that ships, `O(2)_a`, its CLOSED-FORM membership criterion
# `Q ê_a = ê_a`, which needs no sampling (ERR-072's trap: a finite sample of a
# continuous group generates a finite SUBgroup and certifies falsely).
# ===========================================================================


def _reflection(axis: str) -> np.ndarray:
    m = np.eye(3)
    m[_AXIS[axis], _AXIS[axis]] = -1.0
    return m


def _rotation(axis: str, theta: float) -> np.ndarray:
    k = _E[axis]
    cross = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return _I3 + np.sin(theta) * cross + (1.0 - np.cos(theta)) * (cross @ cross)


def _closure(generators, atol: float = 1e-9) -> list[np.ndarray]:
    """BFS closure — an ALGORITHM independent of production's `_close_group`
    only in the sense that it is written here, from the generators the
    definition names.  Its own guard: a runaway closure raises rather than
    silently returning a proper subgroup."""
    out, frontier = [_I3], [_I3]
    while frontier:
        nxt = []
        for a in frontier:
            for g in generators:
                b = g @ a
                if not any(np.abs(b - c).max() <= atol for c in out):
                    out.append(b)
                    nxt.append(b)
        frontier = nxt
        if len(out) > 200:
            raise RuntimeError("closure blew up — the generators are wrong")
    return out


@dataclasses.dataclass(frozen=True)
class _Tag:
    """The reference's own naming of a group — kept OFF the production tag type on
    purpose, so this construction cannot inherit a spelling from the thing it
    checks.  ``axis`` and ``order`` are the two parameters the shipped families
    carry; exactly one is set."""

    kind: str
    axis: str = ""
    order: int = 0

    def __str__(self) -> str:            # stable, readable parametrize ids
        return f"{self.kind}{self.axis or (self.order or '')}"


def _reference_elements(tag: _Tag) -> list[np.ndarray]:
    if tag.kind == "Trivial":
        return [_I3]
    if tag.kind == "Mirror":
        return [_I3, _reflection(tag.axis)]
    if tag.kind == "Cn":
        return [_rotation("z", 2.0 * np.pi * k / tag.order) for k in range(tag.order)]
    if tag.kind == "Dnh":
        # D_nh about x.  n == 1: generated by the two mirrors CONTAINING x
        # (`[M]` 2026-09-03 production realizes D_1h as {e, sigma_y, sigma_z,
        # C_2^x} — the group they generate).  n >= 2: the n-fold rotation about
        # x, the horizontal mirror sigma_x, and one vertical mirror.
        gens = ([_reflection("y"), _reflection("z")] if tag.order == 1
                else [_rotation("x", 2.0 * np.pi / tag.order),
                      _reflection("x"), _reflection("y")])
        return _closure(gens)
    raise KeyError(tag.kind)


def _reference_membership(tag: _Tag):
    if tag.kind == "O2":
        axis = _E[tag.axis]
        # O(2)_a is the stabiliser of the axis VECTOR: every rotation about a
        # fixes it, and every reflection in a plane CONTAINING a fixes it; the
        # elements that send a -> -a (sigma_a, C_2 about a perpendicular axis)
        # are in D_inf_h and NOT in O(2)_a.  Closed form, no sample.
        return lambda q: bool(np.abs(q @ axis - axis).max() <= 1e-9)
    elements = _reference_elements(tag)
    return lambda q: any(np.abs(q - c).max() <= 1e-9 for c in elements)


def _reference_covered(h_elements, gamma_elements, in_kappa) -> bool:
    r"""``H ⊆ Γ·K``: every ``h`` is ``γk``, i.e. ``∃γ: γ⁻¹h ∈ K``."""
    return all(any(in_kappa(g.T @ h) for g in gamma_elements) for h in h_elements)


_FINITE_TAGS: dict[_Tag, SubgroupOfO3] = {
    _Tag("Trivial"): S.Trivial,
    _Tag("Mirror", axis="x"): S.Mirror("x"),
    _Tag("Mirror", axis="y"): S.Mirror("y"),
    _Tag("Mirror", axis="z"): S.Mirror("z"),
    _Tag("Cn", order=2): S.Cn(2),
    _Tag("Cn", order=3): S.Cn(3),
    _Tag("Cn", order=4): S.Cn(4),
    _Tag("Dnh", order=1): S.Dnh(1),
    _Tag("Dnh", order=2): S.Dnh(2),
}
_KAPPA_TAGS: dict[_Tag, SubgroupOfO3] = dict(_FINITE_TAGS) | {
    _Tag("O2", axis="x"): S.O2("x"), _Tag("O2", axis="z"): S.O2("z"),
}


# ===========================================================================
# Group A — the coverage predicate's own laws.
#   claim layer: none (pure math).  pillar: closed-form + independent construction.
# ===========================================================================


class TestR3CoverageAgreesWithAnIndependentConstruction:
    r"""**A1, the keystone.** ``H.is_subset_of_product(Γ, K)`` against an element-set
    enumeration written from the groups' definitions.

    `[M]` 2026-09-03, **1188 triples, 0 disagreements** — 495 True / 693 False, so the
    reference is not a constant (a predicate that answered ``True`` everywhere would
    also agree with a reference that did).  Cost `[M]` **0.7 s** for all four sections
    of ``scratch/_r3_predicate_laws.py`` including interpreter start.

    ⛔ **Where the independence STOPS** (`lessons` L72c): both sides speak the same
    standard setting — a 3x3 orthogonal matrix in the same basis — because
    "is this matrix an element of that group" has no basis-free spelling.  The
    ALGORITHM is what is independent: production never enumerates ``Γ·K``, and the
    reference never touches a Lie algebra.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("h_tag", sorted(_FINITE_TAGS, key=str))
    def test_a1_every_gamma_kappa_pair_agrees(self, h_tag) -> None:
        h_group = _FINITE_TAGS[h_tag]
        h_elements = _reference_elements(h_tag)
        rows = trues = 0
        for gamma_tag, gamma in _FINITE_TAGS.items():
            for kappa_tag, kappa in _KAPPA_TAGS.items():
                want = _reference_covered(
                    h_elements, _reference_elements(gamma_tag),
                    _reference_membership(kappa_tag),
                )
                got = h_group.is_subset_of_product(gamma=gamma, kappa=kappa)
                assert got == want, (
                    f"{h_group.name} ⊆ {gamma.name}·{kappa.name}: "
                    f"production {got}, independent construction {want}"
                )
                rows += 1
                trues += int(want)
        # Non-vacuity (`vv-principles` #19): a reference that answers the same thing
        # everywhere agrees with ANY predicate, so the row would carry no
        # information.  `[M]` 2026-09-03 exactly ONE of the nine H's is constant,
        # and it is constant BY THEOREM — the trivial group lies in every product,
        # because e = e*e.  Stated as its own branch so the day another H goes
        # constant the suite says WHICH claim it lost (`lessons` L73a).
        assert rows == 99
        if h_group.is_trivial:
            assert trues == rows, "{e} ⊆ ΓK for every Γ, K — the theorem"
        else:
            assert 0 < trues < rows, (
                f"{h_group.name}: the reference answered {trues}/{rows} True — a "
                f"constant reference agrees with any predicate (vv-principles #19)"
            )

    @pytest.mark.foundation
    def test_a1b_a_continuous_H_is_covered_only_by_a_matching_continuous_K(self) -> None:
        r"""The first conjunct on a CONTINUOUS ``H``: no finite ``Γ`` can cover
        ``O(2)_x`` over a finite ``K``, and ``O(2)_z`` is not covered over ``O(2)_x``
        however large ``Γ`` is — ``ΓK`` is a finite union of cosets of ``K⁰``, so a
        1-parameter ``H⁰`` outside ``K⁰`` cannot fit in it."""
        for gamma in (S.Trivial, S.Mirror("x"), S.Dnh(2), S.Cn(4)):
            assert not S.O2("x").is_subset_of_product(gamma=gamma, kappa=S.Trivial)
            assert not S.O2("x").is_subset_of_product(gamma=gamma, kappa=S.Mirror("x"))
            assert not S.O2("z").is_subset_of_product(gamma=gamma, kappa=S.O2("x"))
            # the matching axis IS covered, so the row above is not "always False"
            assert S.O2("x").is_subset_of_product(gamma=gamma, kappa=S.O2("x"))


class TestR3CoverageIsASubsetOfAProduct:
    r"""**A2-A6** — the theorem's two halves, the finiteness refusal, the genuine
    product, the lattice compatibility law, and the two DECLARED-BLIND directions.

    The theorem (production's own docstring, re-derived here): ``ΓK = ⋃_γ γK`` is a
    finite union of cosets of ``K⁰``; those cosets are disjoint-or-equal, hence
    clopen in ``ΓK``, so the connected ``H⁰ ∋ e`` lands in the one containing ``e``,
    which is ``K⁰`` itself.  Given that, ``H = ⨆_r rH⁰ ⊆ ΓK ⟺ every r ∈ ΓK``.
    """

    # --- A2: the two conjuncts, isolated -----------------------------------
    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "label,h,gamma,kappa,component,reps,verdict",
        [
            ("component fails", S.O2("x"), S.Mirror("x"), S.Trivial, False, False, False),
            ("a rep fails", S.Mirror("y"), S.Mirror("z"), S.Trivial, True, False, False),
            ("both pass", S.Mirror("z"), S.Mirror("z"), S.Trivial, True, True, True),
            ("component alone", S.Cn(2), S.Trivial, S.O2("x"), True, False, False),
        ],
    )
    def test_a2_each_conjunct_has_a_discriminating_row(
        self, label, h, gamma, kappa, component, reps, verdict
    ) -> None:
        r"""`vv-principles` #11 + #17's granularity clause: a two-conjunct predicate
        is TWO claims, and mutating it as a unit certifies only the first reachable
        one.  ``component alone`` is the row that matters most — ``C_2`` (about z)
        has the trivial identity component, so the FIRST conjunct passes, and only
        the representative search refuses it."""
        got_component = kappa.realization.component.contains(h.realization.component)
        got_reps = all(
            any(kappa.realization.contains_element(g.inverse() @ r)
                for g in gamma.realization.elements)
            for r in h.realization.representatives
        )
        assert (got_component, got_reps) == (component, reps), label
        assert h.is_subset_of_product(gamma=gamma, kappa=kappa) is verdict, label

    # --- A3: Gamma must be FINITE ------------------------------------------
    @pytest.mark.foundation
    def test_a3_a_continuous_covering_factor_is_refused(self) -> None:
        r"""The elements of ``Γ`` are enumerated, so a continuous ``Γ`` has no
        answer.  ⛔ No SHIPPED table row reaches this arm — every geometry's
        ``unspent`` is finite by ``AngularSymmetry.__post_init__`` — so the witness
        is MANUFACTURED (`plan-authoring` §6c: the gate lands with the case it
        catches, and the case is built here)."""
        with pytest.raises(ValueError, match="must be finite"):
            S.Mirror("x").is_subset_of_product(gamma=S.O2("x"), kappa=S.Trivial)
        with pytest.raises(ValueError, match=r"O2_x is continuous"):
            S.Mirror("x").is_subset_of_product(gamma=S.O2("x"), kappa=S.Trivial)
        # positive control: the same call with a FINITE Gamma answers
        assert S.Mirror("x").is_subset_of_product(gamma=S.Mirror("x"), kappa=S.Trivial) is True
        # ... and a continuous K is fine — only the enumerated factor must be finite
        assert S.Mirror("y").is_subset_of_product(gamma=S.Trivial, kappa=S.O2("x")) is True

    # --- A4: the genuine product -------------------------------------------
    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "h,gamma,kappa",
        [
            (S.Mirror("x"), S.Mirror("x"), S.Mirror("y")),
            (S.Mirror("x"), S.Mirror("x"), S.O2("x")),
            (S.Mirror("x"), S.Cn(2), S.Mirror("y")),
            (S.Cn(2), S.Mirror("x"), S.O2("x")),
        ],
    )
    def test_a4_gamma_is_load_bearing_on_a_genuine_product(self, h, gamma, kappa) -> None:
        r"""⛔ **The shipped table never exercises the product.** `[M]` 2026-09-03 all
        four geometry rows have a TRIVIAL factor — slab/sphere ``(O(2)_x, {e})``,
        cylinder ``({e}, D_1h)``, cartesian2d ``({e}, σ_z)`` — so ``Γ·K`` is always
        ``Γ`` or ``K`` and the product structure has **no shipped witness**.  It IS
        reachable in the spellable vocabulary: `[M]` **137** ``(H, Γ, K)`` triples
        are covered WITH ``Γ`` and refused WITHOUT it.  Four of them ship here."""
        assert h.is_subset_of_product(gamma=gamma, kappa=kappa) is True
        assert h.is_subset_of_product(gamma=S.Trivial, kappa=kappa) is False, (
            "Gamma is not load-bearing on this row — it is not a genuine product"
        )

    # --- A5: the compatibility law with the existing lattice ----------------
    @pytest.mark.foundation
    def test_a5_a_trivial_gamma_reduces_to_containment(self) -> None:
        r"""``H ⊆ {e}·K ⟺ K ⊇ H`` — the new predicate and the shipped lattice
        (:meth:`SubgroupOfO3.contains`) must agree wherever both are defined
        (`vv-principles` #15: ship an order relation AND a predicate that must
        respect it, gate the compatibility law).  `[M]` 99 pairs; neither half can
        be wrong alone without this reddening."""
        rows = agreements = 0
        for h in _KAPPA_TAGS.values():
            for kappa in _KAPPA_TAGS.values():
                rows += 1
                agreements += int(
                    h.is_subset_of_product(gamma=S.Trivial, kappa=kappa) == kappa.contains(h)
                )
        assert agreements == rows, f"{rows - agreements} of {rows} pairs disagree"
        assert rows == 121

    @pytest.mark.foundation
    def test_a5b_the_law_is_monotone_in_gamma(self) -> None:
        r"""``Γ ⊆ Γ'`` ⟹ ``H ⊆ ΓK`` implies ``H ⊆ Γ'K`` — the fold licence can only
        WIDEN when the unspent symmetry does.  `[M]` the shipped instance:
        ``D_1h ⊇ σ_z``, and the cylinder admits every fold cartesian2d does."""
        assert S.Dnh(1).contains(S.Mirror("z"))
        for h in (S.Mirror("z"), S.Mirror("y"), S.Trivial, S.Cn(2)):
            if h.is_subset_of_product(gamma=S.Mirror("z"), kappa=S.Trivial):
                assert h.is_subset_of_product(gamma=S.Dnh(1), kappa=S.Trivial), h.name
        # and the widening is STRICT on sigma_y — the cylinder/cartesian2d split
        assert S.Mirror("y").is_subset_of_product(gamma=S.Dnh(1), kappa=S.Trivial)
        assert not S.Mirror("y").is_subset_of_product(gamma=S.Mirror("z"), kappa=S.Trivial)

    # --- A6: the DECLARED-BLIND directions ---------------------------------
    @pytest.mark.foundation
    def test_a6_the_inverse_direction_is_a_THEOREM_and_cannot_be_witnessed(self) -> None:
        r"""⛔ **A declared-blind row — its green is a LICENCE to read A1-A5 as
        coverage, not itself coverage** (`lessons` L72f).

        Spelling the search ``∃γ: γ r ∈ K`` instead of ``∃γ: γ⁻¹ r ∈ K`` cannot
        change any answer, because the existential ranges over a GROUP and a group
        is closed under inversion (substitute ``γ' = γ⁻¹``).  `[M]` 2026-09-03,
        **0 disagreements over 891** ``(H, Γ, K)`` triples — including ``Γ = C_4``,
        whose elements are NOT involutions, so the nullity is the group axiom and
        not an artefact of the shipped table's involution-only ``Γ``s.

        ⟹ the battery arm "γ instead of γ⁻¹" MUST stay green; a red there means the
        predicate stopped ranging over a group."""
        disagreements = 0
        for h in _FINITE_TAGS.values():
            for gamma in _FINITE_TAGS.values():
                for kappa in _KAPPA_TAGS.values():
                    left = all(
                        any(kappa.realization.contains_element(g.inverse() @ r)
                            for g in gamma.realization.elements)
                        for r in h.realization.representatives
                    )
                    right = all(
                        any(kappa.realization.contains_element(g @ r)
                            for g in gamma.realization.elements)
                        for r in h.realization.representatives
                    )
                    disagreements += int(left != right)
        assert disagreements == 0
        # the theorem's premise, gated so the day it stops holding the suite says WHICH
        # claim it lost (`lessons` L73a's premise-gating shape)
        for gamma in _FINITE_TAGS.values():
            inverses = [g.inverse() for g in gamma.realization.elements]
            for inv in inverses:
                assert any(np.abs(inv.linear - g.linear).max() <= 1e-12
                           for g in gamma.realization.elements), gamma.name

    @pytest.mark.foundation
    def test_a6b_the_SIDE_of_the_product_is_unwitnessable_in_this_vocabulary(self) -> None:
        r"""⛔ **The second declared-blind row.**  ``ΓK`` and ``KΓ`` differ in
        general (they coincide iff ``Γ`` normalises ``K``, or the product is a
        group).  `[M]` 2026-09-03 over the **81** finite ``(Γ, K)`` pairs spellable
        with the shipped vocabulary, ``ΓK == KΓ`` as SETS on **all 81** — every
        member is built from coordinate-axis rotations and coordinate-plane
        mirrors, and those products commute setwise.  A 45-degree mirror would
        separate them and :meth:`SubgroupOfO3.Mirror` admits only ``x``/``y``/``z``.

        ⟹ a "left coset instead of right coset" battery arm is a MEASURED null, not
        a blind gate.  Recorded here so nobody credits a green there as coverage."""
        differ = 0
        for gamma in _FINITE_TAGS.values():
            for kappa in _FINITE_TAGS.values():
                gk = {np.round(g.linear @ k.linear, 12).tobytes()
                      for g in gamma.realization.elements
                      for k in kappa.realization.elements}
                kg = {np.round(k.linear @ g.linear, 12).tobytes()
                      for g in gamma.realization.elements
                      for k in kappa.realization.elements}
                differ += int(gk != kg)
        assert differ == 0, (
            f"{differ} of 81 pairs now separate the two sides — the SIDE of the "
            f"product has become witnessable; give it its own gate"
        )


# ===========================================================================
# Group B — stage 0's answers against a FROZEN pre-carve baseline.
#   claim layer: none (a regression pin).  pillar: a table frozen before the edit.
# ===========================================================================

_BASELINE = NUMERICS_DATA / "r3_stage_baseline.json"


def _baseline_measures() -> list[tuple[str, DiscreteMeasure | str]]:
    """The 15 shipped rules of the #429 exit instrument + 7 MANUFACTURED folds.

    The manufactured rows exist because the shipped registry has no fold at all
    (`[M]` its four specs are GaussLegendre1D / LebedevSphere / LevelSymmetricSN /
    ProductQuadrature) and because ``folded_product`` is the ONLY fold the tree
    ships — so without them the carve's whole subject is one measure."""
    rules = [
        ("lebedev(5)", lambda: Quadrature.lebedev(5)),
        ("lebedev(9)", lambda: Quadrature.lebedev(9)),
        ("lebedev(13)", lambda: Quadrature.lebedev(13)),
        ("lebedev(17)", lambda: Quadrature.lebedev(17)),
        ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4)),
        ("level_symmetric(6)", lambda: Quadrature.level_symmetric(6)),
        ("level_symmetric(8)", lambda: Quadrature.level_symmetric(8)),
        ("product(4,8)", lambda: Quadrature.product(4, 8)),
        ("product(8,8)", lambda: Quadrature.product(8, 8)),
        ("gauss_legendre(2)", lambda: Quadrature.gauss_legendre(2)),
        ("gauss_legendre(8)", lambda: Quadrature.gauss_legendre(8)),
        ("gauss_legendre(16)", lambda: Quadrature.gauss_legendre(16)),
        ("product(4,4)", lambda: Quadrature.product(4, 4)),
        ("folded_product(2,4)", lambda: Quadrature.folded_product(2, 4)),
        ("folded_product(4,8)", lambda: Quadrature.folded_product(4, 8)),
    ]
    out: list[tuple[str, DiscreteMeasure | str]] = [(n, f().measure) for n, f in rules]
    base = Quadrature.product(4, 8).measure
    folds = [(f"product(4,8)/sigma_{a}", S.Mirror(a)) for a in "xyz"]
    folds += [("product(4,8)/Cn2", S.Cn(2)), ("product(4,8)/D1h", S.Dnh(1)),
              ("product(4,8)/D2h", S.Dnh(2)), ("product(4,8)/O2x", S.O2("x"))]
    for name, group in folds:
        try:
            out.append((name, base.quotient(group)))
        except Exception as exc:
            out.append((name, f"UNCONSTRUCTIBLE {type(exc).__name__}: {exc}"))
    gl8 = Quadrature.gauss_legendre(8).measure
    for tag, group in (("sigma_x", S.Mirror("x")), ("Dinfh", S.Dinfh), ("O3", S.O3)):
        try:
            out.append((f"gauss_legendre(8)/{tag}", gl8.quotient(group)))
        except Exception as exc:
            out.append((f"gauss_legendre(8)/{tag}",
                        f"UNCONSTRUCTIBLE {type(exc).__name__}: {exc}"))
    return out


#: The five cells R3 is RULED to move — everything else must be bit-equal to the
#: pre-carve table.  `[M]` 2026-09-03 over 21 constructible measures x 4 geometries
#: x 2 stages = **168 cells**; 5 moved, ALL stage 0, ALL True -> False.
_RULED_MOVES = {
    # D1 closed: the owed closure was read as a fold licence
    ("folded_product(2,4)", "cartesian2d.domain"): (True, False),
    ("folded_product(4,8)", "cartesian2d.domain"): (True, False),
    ("product(4,8)/sigma_y", "cartesian2d.domain"): (True, False),
    # NOT named in the plan: sigma_x is in NEITHER geometry's unspent group
    ("product(4,8)/sigma_x", "cylinder.domain"): (True, False),
    ("product(4,8)/sigma_x", "cartesian2d.domain"): (True, False),
}


class TestR3StageZeroBaselineIsUnchangedExceptWhereRuled:
    r"""**B1** — the carve's whole behaviour delta, against a table frozen BEFORE it.

    "the slab and sphere admissions unchanged" is a COMPARISON, not a claim
    (`AGENT.md` §0.5).  The baseline was captured from a ``git archive HEAD`` shadow
    tree, so it is reproducible after the carve landed on the working tree —
    ``scratch/_r3_capture.py <tree_dir|-> <out.json>``.

    `[M]` 2026-09-03 the delta is **5 of 168 cells**, and two of them are NOT in the
    plan's list: ``product(4,8)/σ_x`` is refused for the CYLINDER as well as for
    cartesian2d, because ``σ_x ∉ D_1h`` — a 1-D cylinder's radial cosine is not a
    symmetry direction (folding it would collapse the inward and outward sweeps).
    That is correct and is the second §6c witness for the ``Mirror('y') → D_1h``
    ruling; it is recorded here so it is a RULED move rather than an accident.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("geometry", ["slab", "sphere", "cylinder", "cartesian2d"])
    def test_b1_every_cell_matches_the_frozen_pre_carve_table(self, geometry) -> None:
        baseline = json.loads(_BASELINE.read_text())
        symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]
        checked = moved = 0
        for name, measure in _baseline_measures():
            want = baseline[name]
            if isinstance(measure, str):
                assert "note" in want, f"{name}: was constructible before the carve"
                assert want["note"].split(":")[0] == measure.split(":")[0], name
                continue
            assert "note" not in want, f"{name}: was UNCONSTRUCTIBLE before the carve"
            assert measure.support.name == want["support"], name
            assert int(measure.n_points) == want["n_points"], name
            for stage, verb in (("domain", symmetry.admits_domain),
                                ("symmetry", symmetry.admits_symmetry)):
                key = f"{geometry}.{stage}"
                got = bool(verb(measure))
                ruled = _RULED_MOVES.get((name, key))
                checked += 1
                if ruled is None:
                    assert got == want[key], (
                        f"{name} / {key}: UNRULED move {want[key]} -> {got}"
                    )
                else:
                    moved += 1
                    assert (want[key], got) == ruled, (
                        f"{name} / {key}: the ruled move is {ruled}, measured "
                        f"({want[key]}, {got})"
                    )
        # `[M]` 25 rows in the table, 7 of them UNCONSTRUCTIBLE (a fold of a
        # REDUCED domain, and the three folds with no catalogue entry) -> 18
        # constructible measures x 2 stages.  The count is asserted so a silently
        # SHRINKING denominator cannot make this row pass by measuring less.
        assert checked == 36, checked
        assert moved == sum(1 for _, k in _RULED_MOVES if k.startswith(geometry))

    @pytest.mark.foundation
    def test_b1b_every_ruled_move_is_a_TIGHTENING_on_a_FOLD(self) -> None:
        """No admission is WIDENED by R3, and every move is on a folded rule —
        so no shipped unfolded configuration can have changed."""
        for (name, key), (before, after) in _RULED_MOVES.items():
            assert (before, after) == (True, False), (name, key)
            assert key.endswith(".domain"), (name, key)
            assert "sigma" in name or "folded" in name, name

    @pytest.mark.foundation
    def test_b1c_admits_domain_is_TOTAL_over_every_pair(self) -> None:
        r"""⚠ **An honest-scope row.** ``admits_domain`` no longer has an arm that
        raises — but `[M]` 2026-09-03 the RETIRED arm had **no constructible
        witness** either: reaching ``spent_group``'s ``NotImplementedError`` needed a
        fold of a REDUCED domain, and both spellings are refused one layer down
        (``gauss_legendre(8).measure.quotient(σ_x)`` -> *"no catalogue entry for
        S^2/O2_x/sigma_x"*; ``.quotient(D_inf_h)`` -> *"a quotient is defined only
        for a Dinfh-invariant measure"*).  So totality is a claim about the CODE,
        witnessed by :class:`TestR3RetiredNames`, and this row is a REGRESSION pin —
        it would red if a future entry made such a fold constructible and the new
        predicate mishandled it."""
        for name, measure in _baseline_measures():
            if isinstance(measure, str):
                continue
            for symmetry in GEOMETRY_ANGULAR_SYMMETRY.values():
                assert isinstance(symmetry.admits_domain(measure), bool), name


# ===========================================================================
# Group C — what a gate CAN say about the PHYSICS claim behind `unspent`.
#   claim layer: FLUX SHAPE.  pillar: none needed — the reference is the SAME
#   solve read at two ordinates, i.e. an EXACT symmetry of the discrete operator.
# ===========================================================================


def _one_group_mixture(sigma_t: float = 1.0, sigma_s: float = 0.5) -> Mixture:
    zero = np.zeros(1)
    return Mixture(
        SigC=np.array([sigma_t - sigma_s]), SigL=zero.copy(), SigF=zero.copy(),
        SigP=zero.copy(), SigT=np.array([sigma_t]),
        SigS=[csr_matrix(np.array([[sigma_s]]))],
        Sig2=[csr_matrix(np.zeros((1, 1)))], chi=zero.copy(),
    )


def _ordinate_permutation(quad: Quadrature, matrix: np.ndarray, atol: float = 1e-10):
    """``n -> n'`` with ``Omega_{n'} = M Omega_n``, or ``None`` if the rule is not
    closed under ``M``."""
    omega = np.column_stack([quad.mu_x, quad.mu_y, quad.mu_z])
    mapped = omega @ np.asarray(matrix).T
    out = []
    for row in mapped:
        distance = np.abs(omega - row).max(axis=1)
        j = int(np.argmin(distance))
        if distance[j] > atol:
            return None
        out.append(j)
    perm = np.asarray(out)
    return perm if len(set(out)) == len(out) else None


def _interior(solution) -> np.ndarray:
    field = solution.angular_flux
    return np.asarray(field.interior.values if hasattr(field, "interior") else field.values)


def _evenness(psi: np.ndarray, perm: np.ndarray) -> float:
    axis = 0 if psi.shape[0] == len(perm) else 1
    return float(np.abs(np.take(psi, perm, axis=axis) - psi).max()
                 / max(np.abs(psi).max(), 1e-300))


@pytest.fixture(scope="module")
def _cylinder_solution():
    quad = Quadrature.folded_product(4, 8)
    mesh = Mesh1D(edges=np.linspace(0.0, 2.0, 9), mat_ids=np.zeros(8, dtype=int),
                  coord=CoordSystem.CYLINDRICAL,
                  bc_left=BC("reflective"), bc_right=BC("vacuum"))
    source = np.zeros((quad.N, 1, 8))
    source[:, :, :4] = 1.0                     # r-ASYMMETRIC: no flat-flux nulling
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000)
    return quad, _interior(solution)


@pytest.fixture(scope="module")
def _plane_solution():
    quad = Quadrature.product(4, 8)
    mesh = Mesh2D(edges_x=np.linspace(0.0, 2.0, 6), edges_y=np.linspace(0.0, 1.0, 4),
                  mat_map=np.zeros((5, 3), dtype=int), coord=CoordSystem.CARTESIAN,
                  bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
                  bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"))
    source = np.zeros((quad.N, 1, 5, 3))
    source[:, :, :2, :1] = 1.0                 # x- AND y-asymmetric corner source
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000)
    return quad, _interior(solution)


@pytest.fixture(scope="module")
def _slab_solution():
    quad = Quadrature.gauss_legendre(8)
    mesh = Mesh1D(edges=np.linspace(0.0, 2.0, 9), mat_ids=np.zeros(8, dtype=int),
                  coord=CoordSystem.CARTESIAN,
                  bc_left=BC("vacuum"), bc_right=BC("vacuum"))
    source = np.zeros((quad.N, 1, 8))
    source[:, :, :4] = 1.0
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000)
    return quad, _interior(solution)


_MIRRORS = {"sigma_x": np.diag([-1.0, 1.0, 1.0]),
            "sigma_y": np.diag([1.0, -1.0, 1.0]),
            "sigma_z": np.diag([1.0, 1.0, -1.0])}


class TestR3UnspentIsMeasuredOnTheSolver:
    r"""**C1-C3** — ``unspent`` is a claim about the SOLUTION, and this is the most a
    gate can say about it.

    ⛔ **What a gate CANNOT do.**  ``unspent`` is a statement about the CONTINUUM
    transport equation of a geometry — "psi is even under this finite group acting on
    directions alone".  No test verifies that; it is a derivation (the archivist's
    page), and a wrong derivation would make the geometry's whole reduction wrong,
    not just its registry row.  A test cannot even enumerate the claim's domain: it
    ranges over every source, every material, every boundary law.

    ⭐ **What a gate CAN do, and it is the load-bearing half.**  Every downstream fold
    rests on the DISCRETE solver's answer, and there the claim is exactly checkable:
    solve a deliberately asymmetric fixed-source problem and compare ``psi`` at
    ordinate ``n`` with ``psi`` at the ordinate the group element maps it to.  A
    positive leg per element of ``unspent``, a NEGATIVE leg per element outside it —
    without the second the row is `vv-principles` #19-blind, because a solver whose
    psi is even under EVERYTHING would pass the positive leg too.

    `[M]` 2026-09-03, all three fixtures in **1.0 s** total:

    ==============  ================  =====================  ==================
    geometry        rule              ``sigma_z``            ``sigma_y``
    ==============  ================  =====================  ==================
    cartesian2d     product(4,8)      **0.0** (EVEN)         6.043e-01 (NOT)
    cylinder        folded_product    **0.0** (EVEN)         no permutation
    slab            gauss_legendre    0.0 but VACUOUS        0.0 but VACUOUS
    ==============  ================  =====================  ==================

    and ``sigma_x`` moves every one of them: 8.175e-01 / 7.224e-01 / 6.493e-01.

    ⚠ **The trap the vacuity guard caught.**  On the SLAB the ``sigma_y`` and
    ``sigma_z`` permutations are the IDENTITY (a 1-D polar rule carries
    ``mu_y = mu_z = 0`` for every ordinate — ERR-080's territory), so their
    ``0.0`` readings are `vv-principles` Mode-8 tautologies, not evenness.  Reported
    without the guard they would have contradicted ``unspent = Trivial``.  Every
    positive leg here therefore asserts the permutation MOVES ordinates.
    """

    @pytest.mark.foundation
    def test_c1_cartesian2d_is_even_under_sigma_z_and_under_NOTHING_else(
        self, _plane_solution
    ) -> None:
        quad, psi = _plane_solution
        plane = GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"]
        assert plane.unspent == S.Mirror("z")

        perm = _ordinate_permutation(quad, _MIRRORS["sigma_z"])
        assert perm is not None and int(np.count_nonzero(perm != np.arange(len(perm)))) == 32
        assert _evenness(psi, perm) == 0.0, "z-uniformity is EXACT in the 2-D operator"

        for name in ("sigma_x", "sigma_y"):
            other = _ordinate_permutation(quad, _MIRRORS[name])
            assert other is not None
            assert int(np.count_nonzero(other != np.arange(len(other)))) == 24
            assert _evenness(psi, other) > 0.5, (
                f"{name} must MOVE the solution — it is outside cartesian2d's unspent "
                f"group, and a fold by it empties two sweep quadrants (D1)"
            )

    @pytest.mark.foundation
    def test_c2_the_cylinder_is_even_under_sigma_z_and_sigma_x_moves_it(
        self, _cylinder_solution
    ) -> None:
        r"""⛔ **The ``sigma_y`` half of ``D_1h`` is UN-WITNESSABLE on the solver
        path, by construction.**  `[M]` 2026-09-03 a cylindrical ``SNMesh`` admits
        only CARRYING quadratures (``assert_carrying_quadrature``): **15 of 15**
        ``folded_product`` rules pass, **0 of 20** ``product`` / ``lebedev`` /
        ``level_symmetric`` rules do.  Every admissible rule IS the ``sigma_y``
        quotient, so no cylindrical solve ever stores both signs of ``mu_y`` and no
        test can ask whether psi is even under the mirror the rule already spent.

        ⟹ ``sigma_y in unspent`` is, for the cylinder, a STRUCTURAL commitment the
        tree already makes one layer down — the shipped solver refuses every
        unfolded rule — and the assertion below is the strongest available form of
        it.  Read with `vv-principles` Mode 12: the functional (a folded solve)
        annihilates the error class, so no fixture can be found; saying so is the
        honest deliverable, not a green row."""
        quad, psi = _cylinder_solution
        cylinder = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        assert cylinder.unspent == S.Dnh(1)
        assert cylinder.unspent.contains(S.Mirror("y"))
        assert cylinder.unspent.contains(S.Mirror("z"))
        assert not cylinder.unspent.contains(S.Mirror("x"))

        # the sigma_y half: no permutation exists — the rule already spent it
        assert _ordinate_permutation(quad, _MIRRORS["sigma_y"]) is None
        assert quad.measure.quotient_group == S.Mirror("y")

        # the sigma_z half: measurable, and EXACT
        perm = _ordinate_permutation(quad, _MIRRORS["sigma_z"])
        assert perm is not None and int(np.count_nonzero(perm != np.arange(len(perm)))) == 16
        assert _evenness(psi, perm) == 0.0

        # the negative leg: sigma_x is outside D_1h and moves the solution 72 %
        other = _ordinate_permutation(quad, _MIRRORS["sigma_x"])
        assert other is not None
        assert _evenness(psi, other) > 0.5

    @pytest.mark.foundation
    def test_c3_the_slab_keeps_NOTHING_unspent(self, _slab_solution) -> None:
        r"""``unspent = Trivial`` says a slab solution is even under no non-trivial
        mirror.  ``sigma_x`` (the OWED closure) moves it by 65 %, and the other two
        mirrors act as the IDENTITY on a polar rule — so their zero readings carry
        no information and are asserted VACUOUS rather than read as evenness."""
        quad, psi = _slab_solution
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        assert slab.unspent == S.Trivial
        assert slab.owed == S.Mirror("x")

        perm = _ordinate_permutation(quad, _MIRRORS["sigma_x"])
        assert perm is not None
        assert int(np.count_nonzero(perm != np.arange(len(perm)))) == 8
        assert _evenness(psi, perm) > 0.5, (
            "a slab solution is not even in mu for a general source — which is why "
            "nothing finite is left UNSPENT and sigma_x is only OWED"
        )
        for name in ("sigma_y", "sigma_z"):
            vacuous = _ordinate_permutation(quad, _MIRRORS[name])
            assert vacuous is not None
            assert int(np.count_nonzero(vacuous != np.arange(len(vacuous)))) == 0, (
                f"{name} now MOVES ordinates on a polar rule — its 0.0 evenness "
                f"reading has stopped being a tautology and needs a real gate"
            )


class TestR3TheFoldConsequenceIsTheSweepQuadrant:
    r"""**C4** — the consequence tier of D1, at the tier the sweep consumes.

    The plan's D1 row says two of the four ``(sign mu_x, sign mu_y)`` sweep quadrants
    are EMPTY under the wrongly-admitted fold.  `[M]` 2026-09-03 on the ORDINATE
    cosines (``Quadrature.mu_x``/``mu_y`` — what a 2-D sweep partitions):

    ==========================  ====  ==========================================
    rule                        N     strictly-signed (mu_x, mu_y) quadrants
    ==========================  ====  ==========================================
    ``product(4,8)``            32    4 of 4 (4 each)
    ``product(4,8)/sigma_z``    16    **4 of 4** (2 each) — the LICENSED fold
    ``folded_product(4,8)``     16    **2 of 4** — ``mu_y > 0`` on all 16
    ``product(4,8)/sigma_y``    20    **2 of 4**
    ==========================  ====  ==========================================

    So the discriminator is exact and needs no solve: the fold cartesian2d admits
    keeps every sweep quadrant; the fold it now refuses empties half of them.
    """

    @staticmethod
    def _quadrants(nodes: np.ndarray) -> set[tuple[float, float]]:
        signs = np.sign(np.asarray(nodes)[:, :2])
        return {tuple(row) for row in signs if 0.0 not in row}

    @pytest.mark.foundation
    def test_c4_the_licensed_fold_keeps_all_four_sweep_quadrants(self) -> None:
        base = Quadrature.product(4, 8).measure
        z_fold = base.quotient(S.Mirror("z"))
        assert len(self._quadrants(base.nodes)) == 4
        assert len(self._quadrants(z_fold.nodes)) == 4
        assert GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"].admits_domain(z_fold)

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "label,build",
        [("folded_product(4,8)", lambda: Quadrature.folded_product(4, 8).measure),
         ("product(4,8)/sigma_y",
          lambda: Quadrature.product(4, 8).measure.quotient(S.Mirror("y")))],
    )
    def test_c4b_the_refused_fold_empties_two_of_them(self, label, build) -> None:
        fold = build()
        assert len(self._quadrants(fold.nodes)) == 2, label
        assert not GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"].admits_domain(fold), label
        # ... and the CYLINDER still admits it: the quadrant argument is about a
        # 2-D sweep, not about the fold being illegitimate everywhere
        assert GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain(fold), label


# ===========================================================================
# Group D — the ledger's construction invariants.
# ===========================================================================


class TestR3TheLedgerRefusesAContinuousFiniteRole:
    r"""**D1** — ``unspent`` and ``owed`` are enumerated, so both must be finite.

    ⛔ **No shipped row reaches either arm** — the four table rows are all finite by
    construction — so both witnesses are MANUFACTURED (`plan-authoring` §6c).  And
    the guard is a LOOP over two roles (`vv-principles` #17's granularity clause):
    mutating it as a unit certifies only ``unspent``, the first role checked, so the
    two arms are asserted SEPARATELY with the role name as the discriminating
    message fragment.
    """

    @pytest.mark.foundation
    def test_d1_a_continuous_unspent_is_refused_naming_that_role(self) -> None:
        with pytest.raises(ValueError, match=r"AngularSymmetry\.unspent must be a finite"):
            AngularSymmetry(spent=S.Trivial, unspent=S.O2("x"), owed=S.Mirror("x"))

    @pytest.mark.foundation
    def test_d1b_a_continuous_owed_is_refused_naming_THAT_role(self) -> None:
        """The second arm — reachable only with a finite ``unspent``, so a
        whole-guard mutation would have certified it through the first."""
        with pytest.raises(ValueError, match=r"AngularSymmetry\.owed must be a finite"):
            AngularSymmetry(spent=S.Trivial, unspent=S.Trivial, owed=S.O2("x"))

    @pytest.mark.foundation
    def test_d1c_control_a_continuous_SPENT_is_accepted(self) -> None:
        """The positive control, and the row that keeps the guard from widening:
        ``spent`` is NOT enumerated (it is asked ``contains_element``), and the
        slab's own row spends the continuous ``O(2)_x``."""
        row = AngularSymmetry(spent=S.O2("x"), unspent=S.Trivial, owed=S.Mirror("x"))
        assert row == GEOMETRY_ANGULAR_SYMMETRY["slab"]
        assert not row.spent.realization.is_finite


# ===========================================================================
# Group E — the retired names.
# ===========================================================================


class TestR3RetiredNames:
    r"""**E1** — ``manifold.spent_group`` and the two old field names are gone.

    This is the row that carries ``TestR3StageZeroBaselineIsUnchangedExceptWhereRuled``'s
    ``test_b1c`` totality claim: the ``NotImplementedError`` arm cannot fire because the
    function that raised it does not exist.
    """

    @pytest.mark.foundation
    def test_e1_spent_group_is_absent_from_manifold(self) -> None:
        import orpheus.numerics.manifold as manifold

        assert not hasattr(manifold, "spent_group")
        assert "spent_group" not in manifold.__all__

    @pytest.mark.foundation
    def test_e1b_the_ledger_has_exactly_the_three_roles(self) -> None:
        names = [f.name for f in dataclasses.fields(AngularSymmetry)]
        assert names == ["spent", "unspent", "owed"]
        row = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        for retired in ("continuous_isotropy", "discrete_residual"):
            assert not hasattr(row, retired), retired


# ===========================================================================
# Group F — the selection re-key.
# ===========================================================================


class TestR3TheStageZeroREASONIsPinnedWhereItIsDECIDED:
    r"""**F2/F3** — the wording half of stage 0, which the frozen table deliberately
    does NOT pin.

    ⛔ **The R2 gate did not cover R3 — it RED.**  ``TestR2SelectionIsUnchanged``
    pinned the whole frozen record including every ``rejected`` string, and `[M]`
    2026-09-03 **96 of 96 domain-mismatch strings moved** while **0 of 48** chosen
    specs / parameters / point counts did.  The repair that landed is the right one
    and is NOT duplicated here: that gate now pins stage 0 at its STAGE
    (``stage_of``) and keeps the frozen table as the value anchor.

    ⚠ So this class deliberately carries **no** whole-table row.  Re-pinning the 48
    verdicts here would be a Pattern-2 twin of
    ``tests/numerics/test_invariance.py::TestR2SelectionIsUnchanged``, and a twin
    that drifts is worse than no gate.  What is left is the half that gate hands
    over: the WORDING, pinned once, where the decision is made.
    """

    @pytest.mark.foundation
    def test_f2_the_stage_zero_refusal_names_the_ROLE_that_decided(self) -> None:
        """The message half.  Fragments, not the whole string — and the fragment is
        the ROLE that decided, so a message that reverts to naming the owed closure
        reds here (the D1 regression's own tell)."""
        # ⛔ The row I first wrote asserted the ARROW here and was WRONG, and the
        # landed sibling `test_the_rejection_message_names_the_missing_arrow` is
        # wrong the same way (`[M]` 2026-09-03 it is RED on the tree).  A 1-D rule
        # for a 2-D geometry passes the arrow — `quotient_onto(S^2, S^2/O2_x)` is
        # the entry's OWN quotient map — and is refused by COVERAGE.  Which clause
        # refuses which input is a measurement, not an intuition.
        _, log = select_quadrature("cartesian2d", 5)
        reasons = [r for name, r in log.rejected if name == "GaussLegendre1D"]
        assert reasons, "GaussLegendre1D was not rejected for cartesian2d"
        message = reasons[0]
        assert "domain mismatch" in message
        assert "unspent sigma_z" in message, message
        assert "no descent arrow" not in message, (
            "`domain_refusal` reports the ONE failing clause; naming the other is "
            f"the disjunctive message the refusal verb retired: {message}"
        )
        # and the ARROW clause's own selector witness, on the geometry that has it
        _, slab_log = select_quadrature("slab", 5)
        arrow = [r for name, r in slab_log.rejected if name == "LebedevSphere"]
        assert arrow and "no descent arrow" in arrow[0] and "unspent" not in arrow[0]

    @pytest.mark.foundation
    def test_f3_the_two_refusal_REASONS_are_disjoint_and_each_has_a_witness(self) -> None:
        r"""``domain_refusal`` is a ``-> str | None``, so it collapses TWO guards into
        ONE value and ``is None`` proves only that SOME guard fired (`lessons` L35f).
        Isolate from the INPUT side: one input that fails the ARROW and passes
        coverage, one that passes the arrow and fails COVERAGE, and the two
        fragments must not overlap (`lessons` L43c — a shared fragment makes either
        row a false attribution).

        `[M]` 2026-09-03 the shipped split over 4 geometries x 7 rules is
        **14 arrow / 3 coverage / 0 both** — so a disjunctive message named a
        SATISFIED fact on all 17 refusals before the verb was introduced."""
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        plane = GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"]
        y_fold = Quadrature.folded_product(4, 8).measure

        arrow = slab.domain_refusal(y_fold)          # S^2/O2_x has no arrow onto S^2/sigma_y
        coverage = plane.domain_refusal(y_fold)      # the arrow exists; sigma_y is not unspent
        assert arrow is not None and coverage is not None
        assert quotient_onto(slab.support, y_fold.support) is None
        assert quotient_onto(plane.support, y_fold.support) is not None

        assert "no descent arrow" in arrow and "unspent" not in arrow
        assert "unspent" in coverage and "no descent arrow" not in coverage
        # and the admitted case returns None, so `is None` is not vacuously true
        assert GEOMETRY_ANGULAR_SYMMETRY["cylinder"].domain_refusal(y_fold) is None
