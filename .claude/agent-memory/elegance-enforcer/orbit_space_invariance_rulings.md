---
name: orbit-space-invariance-rulings
description: "#429 tracker 2.2b (the Γ-slot) review — the ambient-vs-orbit-space TWIN CAMPS pattern, the lift/section codomain tell, and the replay-the-deleted-kernel regression instrument"
metadata:
  type: project
---

Reviewed 2026-09-02, uncommitted on `fix/angular-phantom-support` (base `4b7d24c3`):
`geometry/transformation.py`, `numerics/manifold.py`, `numerics/symmetry.py`,
`numerics/quadrature/registry.py`. Findings delivered to the coordinator; I did not edit.

## The durable, cross-transferable rulings

**1. ⭐⭐ A "re-ask the question in the right PLACE" carve produces TWO COHERENT CAMPS, and
the camp boundary is a shared low-level EMBEDDING helper — census that helper, not the
concept.** Here `is_invariant` + `orbit_certificate` moved to the orbit-space reading while
`Quadrature.ordinate_permutation` (3 production consumers) and `tests/_harness/references.py`
stayed on the ambient one; both camps are internally consistent, and `[M]` they contradict on
the shipped `folded_product(4,8)`: `is_invariant(σ_y)=True` vs `ordinate_permutation(σ_y)=None`
(pre-diff both said "no", so the carve CREATED the divergence). The tell that finds it in one
command: grep the shared helper (`_embedded_nodes`) — every consumer of the retired reading is
a consumer of that helper, and a `git grep` over it returns the camp roster while a grep over
the migrated method does not. **Generalise: when a carve changes what a question MEANS rather
than what it computes, the blast radius is the producers of its INPUT, not the callers of its
output.**

**2. The `lift` / `section` / `representative` family has a mechanical honesty tell: compare
`ManifoldMap.domain`/`.codomain` PER ARM.** `[M]` `Quotient.lift` returns three arrows whose
declared types disagree — axial `domain=S^2/O2_x, codomain=D^3` (a Ball, honest: the orbit
BARYCENTRE is off the sphere), mirror `domain=D^2, codomain=S^2`, trivial `domain=S^2,
codomain=S^2` — against its own docstring's "a right inverse of the quotient map on the
*realization*". `section_coordinates` then flattens all three into "points of the BASE", and
`[M]` `base.contains(...)` is **False** and `support.contains(...)` **raises** on the axial
family. That is ERR-080's shape (a chart coordinate dressed as a point of `S²`) re-minted
under a method whose NAME asserts the opposite. ⟹ on any typed-arrow property with per-family
arms, print `domain.name`/`codomain.name` for every arm before reading the docstring.

**3. ⭐ The cheapest decisive regression instrument for a KERNEL SWAP: re-instate the deleted
function in-process, verbatim from the diff, and replay both over the lattice × the shipped
rules.** `[M]` here: **3 of 230** answers changed, all three on the fold, all three the
intended flips (σ_y, C_2, D_2h) — every non-fold family, both axial entries and a bare
interval bit-identical to the retired kernel. Costs ~40 lines, needs no tree edit, and it is
the only thing that separates "the carve is surgical" from "the carve is green".

**4. Brute-force a GROUP-THEORY predicate over the whole lattice; and expect YOUR reference to
be the thing that is wrong first.** `[M]` 529 (G, H) pairs + 138 (H, out-of-lattice motion)
checks against a conjugation reference: 0 disagreements — the criterion (commutation with
`[ê_a]_×` for a connected `SO(2)_a`, `H ⊆ {e,−I}` for `SO(3)`, axis equality for a continuous
other) is correct AND complete. My FIRST run reported **46 disagreements**, all on `H ∈
{SO3, O3}` — my reference, not theirs (both are normal in `O(3)`; I had routed them into a
`Dinfh` membership test). `plan-authoring` §4's verify clause, working: a failed reproduction
is inconclusive until you have diagnosed whose failure it is.

**5. A guard applied to a CANONICAL point (a barycentre, a mean, a projection) can be
structurally unable to fail.** `[M]` `_is_axis_supported(section, a)` on the axial entry has
max off-axis component **0.000e+00 by construction** — an axial orbit's barycentre IS `μ ê_a`.
In-process mutation: forcing it `False` flips exactly one answer to a WRONG value; forcing it
`True` changes nothing. So it decides one answer with a predicate that cannot say no. **The
general tell: a position predicate applied downstream of a canonicalising map is a tautology;
the honest test is the CONTAINMENT that made the map canonical** (here `H ⊇ G⁰`, which also
makes the docstring's own "was answered at 2" clause true — `[M]` it is false today:
`Dinfh.normalises(O2_z)=True` but `O2_z.contains(Dinfh)=False`).

**6. An extraction can RE-duplicate what it factors out.** `[M]` `"points must have shape (n,"`
went **1 → 3** occurrences (`git show HEAD` vs the working tree): `permutes` used to call
`preserves`… (other way round) — `preserves` used to call `permutes` and inherit its validation;
after the carve both carry a byte-identical 6-line block. Check the message-string count across
`HEAD` and the tree whenever a shared helper is lifted out of two methods.

## Per-review inventory (this carve) — for the follow-up pass

VIOLATION ×3: the two-camps twin (§1); `section_coordinates`/`lift` type dishonesty (§2);
the tautological identity-component guard + its false docstring clause (§5).
CONCERN ×8: the re-duplicated shape validation (§6); `_as_columns` at 4 sites papering an
un-normalised `orbit_coordinates` (rank-1 for axial, rank-2 elsewhere — Pattern 7);
the `"SO3"` string sentinel discriminated at 2 sites (a missing type — the identity component
IS a `SubgroupOfO3`); `Quotient.lift`'s 3 arms with `[M]` 44/0/0 reach; `is_normalised_by`
reading `.linear` only on the continuous arms (`[M]` `Dinfh`/`SO3` answer True for a
z-translation, truth False; the finite arm gets it right); `orbit_certificate` spelling
finiteness twice; the spent-group door refusing `(S²/σ_y).quotient(Trivial)` while
`SPHERE.quotient(Trivial)` is admitted AND is what `_ambient_orbit_space()` builds;
`admits_domain` reading `X.by` instead of the group the descent ARROW spends.
Dead-name census: `_invariance_on_points` 0, `_section_nodes` 0, `M.preserves(` 0 —
`_polar_axis_of` **2 live present-tense claims** (`docs/theory/foundations/manifolds.rst:6285`
production docs, and `tests/numerics/test_symmetry.py:1698`).

Probes: `scratch/_22b_elegance_probe{1..10}.py`. Design memo: `scratch/_22b_design.md`.
