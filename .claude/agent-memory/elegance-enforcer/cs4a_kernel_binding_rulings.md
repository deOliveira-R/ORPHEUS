---
name: cs4a-kernel-binding-rulings
description: CS4a-R clear-context review of the "operators born bound" kernel core — the ruled-vs-accidental split, the raw-pairing regression, and the measured refutation of the ClassVar-on-frozen-dataclass trap.
metadata:
  type: project
---

Campaign 1 CS4a ("operators are born bound"), branch `feature/cs1-energy-space`,
`069e2caa`→`49b29391`. Reviewed clear-context 2026-08-21 (CS4a-R Phase 1).
Charter + ⏹ ledger: `.claude/plans/space_and_kernel_binding_campaign.md:688-868`.

**Why:** the campaign mints representation-free interaction kernels
(`orpheus/transport/kernels.py`) and makes the energy-family operators carry a
validated `space` at construction. The kernel/binding doctrine is the ground for
CS2's frame-at-binding and Campaign 2's resolvent.

**How to apply:** read this before reviewing CS4b / CS4c / CS2 — it records which
oddities are RULED (do not re-flag) and which are live.

## The recurring shape this campaign produced (transferable)

⭐ **A "bind it to the space" re-pose can DEMOTE a named typed functional into a
raw pairing, and the tests stay green because the two are coextensive on the
degenerate fixture.** CS4a's K2a replaced
`IntegratedReactionRate(field).evaluate(φ)` with
`space.inner_product(np.asarray(field.values), φ.reshape(ng,1))`. The RULING
(rates read the posing's measure, not the carrier's) is right; the SPELLING
un-named the co-vector and stripped two typed fields to bare ndarrays — in a
campaign whose whole doctrine is "born bound". `[M]` those three lines are the
**only** production sites in the tree that call `space.inner_product` with two
bare ndarrays; every other caller goes through `Field.inner_product`, which runs
`_check_partner`. ⟹ when a campaign re-poses a quantity onto a space, ask
*"which named object now carries the pairing?"* — if the answer is "none, the
solver spells it", the re-pose bound the operator and un-bound the functional.

⭐ **A doctrine can name a concept that has no TYPE.** CS4a's charter defines a
*binding* = "kernel × space [× assignment]" — but the tree has no binding type:
four operators each hand-call `assert_energy_extent_conforms(...)` from their own
`__post_init__`, opt-in and forgettable, with four different `ng` spellings
(three `self.mat_xs.ng`, one `self.coefficient.values.shape[0]`). The 5th binding
author simply omits the call and nothing tells anyone. ⟹ grep for *"is the
doctrine's noun a class?"* on any born-bound / born-typed campaign.

⭐ **A "the mint closes this hazard" claim is about CONSUMERS, and a module with
none closes nothing.** `[M]` 2026-08-21 the F4 writable-cache-alias hazard the
kernels exist to close is still fully open: `mx.sig_s_legendre(0)[0][0,0] = 999`
reaches the loss matrix (`apply_p0_in_scatter` → `[999.0, 1.0]`). All three
kernels have **zero** production and zero doc-source consumers. The producer-side
close (`setflags(write=False)` in `_build_dense_caches`; `[M]` every one of the
~20 consumers only READS) is one line and available today — Pattern 7, at the
definition site, for the 4 live consumers instead of for 0 future ones.

## Ruled — do NOT re-flag on CS4b/CS4c/CS2

- ng-guard live on only 4 of 13 bindings / 192 of 1022 constructions; the
  axis-keyed strengthening is CS2's. Shape-keying was measured UNRUNNABLE
  (`[M]` 275 rows destroyed as a mutation).
- `FissionKernel` has one gate + no production consumer until CS4c (Q2).
- C + the iso pair keep `Optional[space]` until CS4c (`[M]` M2.13: 84 reds today).
- `IntegratedReactionRate` survives for meshed diffusion/SN consumers (CS4b).
- The arity arm has no CS4a witness (first meshed binding = CS4b/CS4c).
- `_pose_space` vs `MaterialMesh.bulk_space` are two spellings of one space —
  and they are **structurally** coextensive, not coincidentally: `[M]`
  `Axis.__post_init__` canonicalizes all-ones weights to `None`, so
  `weights=volumes=[1.0]` and `weights=None` mint the same digest/name. G2.4+G2.5
  red on a swap. Not a twin worth flagging; the collapse is CS2's "single scalar
  bulk mint".

## ⛔ REFUTED — a standing AGENT.md institutional claim, measured false here

AGENT.md §"Institutional knowledge" #5 says: *"under `from __future__ import
annotations`, a stringized `ClassVar[...]` on a frozen dataclass leaf slips past
field-detection and becomes a dataclass field — so dataclass leaves tag with a
plain unannotated attr … do NOT 'fix' a missing ClassVar on a dataclass leaf."*

`[M]` 2026-08-21, Python 3.14, this tree's `.venv`: it does **not** slip past. All
four spellings — bare `ClassVar[int]`, `typing.ClassVar[int]`, a module-aliased
`t2.ClassVar[int]`, and a name-aliased `CV = ClassVar; CV[int]` — are correctly
excluded from `dataclasses.fields()` on a frozen dataclass under
`from __future__ import annotations`. `N2NKernel.multiplicity: ClassVar[int] = 2`
is a real ClassVar (`fields()` → `['matrix']`; `replace(k, multiplicity=7)`
raises `TypeError`).

⟹ the *mechanism* half of #5 is void on this interpreter; the "plain unannotated
attr" house style (`block_role = BlockRole.BULK`) is a **stylistic** choice, not a
requirement. What IS still real: an annotation **without** `ClassVar`
(`block_role: BlockRole = BlockRole.BULK`) genuinely becomes a field — which is
what the house comment is protecting against, imprecisely worded. Re-measure
before repeating either claim. See [[lessons-L-verify-before-you-flag]].
