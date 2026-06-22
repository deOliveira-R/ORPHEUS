---
name: issue-236-phase2-b1-c-fold
description: #236 Ph2 STEP B1 — fold P3 (sweep_cache GeometryCoefficients populator) onto the angular closure for the weighted-diamond c_in/c_out; mint-on-both ABC-vs-Protocol binding ruling + is_linear ClassVar
metadata:
  type: project
---

# #236 Ph2 STEP B1 — c_in/c_out fold P3→closure — PASS-WITH-NITS (none blocking, uncommitted)

Follows [[issue-236-phase2-step-a-tau-carve]] (Step A relocated τ PRODUCTION to the closure).
B1 folds the FIRST of FOUR c_in/c_out rebuild sites — the `GeometryCoefficients.from_mesh_and_quad`
populator (P3, `sweep_cache.py`) — onto the closure. B2/B3/C remain (crosswalk §7, lines 266-279,
VERIFIED tracked → anti-#11 satisfied).

## What B1 does (the 3 changes)
1. NEW public accessor `c_in_per_ordinate`/`c_out_per_ordinate` on BOTH the `@runtime_checkable`
   `PoleAngularClosure` Protocol (`:277-298` pure `...` stubs) AND the ABC `PoleAngularClosureBase`
   (`:422-450` real impl + shared `_gather_per_ordinate` helper). Returns the closure's per-μ-level
   `(M_p,)` c gathered to `(N,)` global-ordinate order (pure permutation). Identity→neutral zeros.
2. P3 rewire: `sweep_cache.py:319-321` reads `closure.c_out_per_ordinate`/`.c_in_per_ordinate`;
   inline `c_out=alpha_out/tau`/`c_in=(1-tau)/tau*alpha_out+alpha_in` DELETED; orphaned
   `alpha_in`/`alpha_out` locals dropped (`:236-237`); `tau` STAYS (feeds tau_inv + mm_a_in_coeff).
3. `is_linear: bool`→`ClassVar[bool]` on the Protocol (`:275`) — cleared a B1-triggered pyright
   assignability error at `geometry.py:465`.

## Producer SSOT verified
The ONE canonical c computation is M-M `__init__` `pole_angular_closure.py:905-920`. P3 now reads it.
Bit-id is STRUCTURAL not faith: consumer (`sweep_cache.py:248`) and producer (`:889-890`) BOTH derive
the cyl partition from the SAME `quad.level_indices`; sphere = `arange(N)` both sides; gather is
`out[level]=values` (pure perm in consumer `global_n` order); shared α slices; 0-ULP τ (Step-A Leg-1).

## ⭐ KEY RULING — mint-on-both verdict: CORRECT locally, but the two cases bind DIFFERENTLY (Cardinal-Rule-2 trigger)
The implementer cited `transverse_coupling_is_facewise` as the mint-on-both precedent. It IS one, but
the binding DIFFERS, and the difference is the whole reason the `is_linear` fix was needed:

| | Protocol member type | consumer `mesh.X` typed against |
|---|---|---|
| `scheme` (`scheme.py:453-459` Protocol / `:770,790` ABC) | **plain `bool`** | the **ABC** `DiscretizationSchemeBase` (`geometry.py:271`) |
| `pole_angular_closure` (B1) | now **`ClassVar[bool]`** | the **Protocol** `PoleAngularClosure` (`geometry.py:204/235/460`) |

- Precedent binds consumer to the **ABC** → ABC's `ClassVar` is the resolution target, assigns clean;
  Protocol's plain-`bool` is NEVER a `ClassVar`-assignment target → no error EVER. Protocol = pure
  `@runtime_checkable` (the test `test_pole_angular_closure.py:49` isinstance + `scheme` `:596`).
- B1 binds consumer to the **Protocol** → assigning an ABC instance (`is_linear: ClassVar`) into a
  Protocol slot declaring `is_linear: bool` (plain) is the pyright error. THAT is why mint-on-both
  here forces every new accessor onto BOTH, AND why `is_linear` had to flip.
- Verified: NO consumer is typed against the bare `DiscretizationScheme` Protocol anywhere; ALL use
  the ABC. The `pole_angular_closure` Protocol-binding is the anomaly.

CLEANER ARCH (the Architectural Opportunity, OUT OF SCOPE for B1): type `MeshSN.pole_angular_closure`
against the ABC `PoleAngularClosureBase` like `scheme` does → future accessors ABC-only (one decl),
`is_linear` fix unnecessary, two mint-on-both cases unified. Touches `geometry.py:204/235/460` sigs +
widens bit-id blast radius → FILE an issue (`module:sn`/`type:improvement`), don't do in B1.

## is_linear ClassVar fix verdict: CORRECT — KEEP, do NOT revert
Genuine class-level trait (ABC `:349` ClassVar, docstring `:265` "Read-only class attribute", every
subclass assigns ClassVar). The Protocol's old `is_linear: bool` was the anomaly (instance member
where impls give a class member). Reverting reintroduces the real `geometry.py:465` error. Caveat for
the record: this makes the `PoleAngularClosure` Protocol use `ClassVar[bool]` where `scheme.py:453`
Protocol uses plain `bool` — a symptom of the binding mismatch above, not a defect in the fix.

## NITS (none blocking)
- CONCERN (do-now, cheap): `_gather_per_ordinate` `pole_angular_closure.py:417` uses `np.empty(N)` +
  full-cover fancy-index overwrite → garbage on ANY ordinate not covered by a level. Holds today
  (sphere arange / cyl level_indices both cover) but the completeness invariant is SILENT + `np.empty`
  makes a future non-covering partition (degenerate-ordinate exclusion / mixed geometry) a Heisenbug.
  `strict=True` guards count, NOT coverage. FIX = `np.zeros(N)` (neutral gap) or assert
  `sum(l.size)==N`. anti-pattern #16 (unstated structural contract).
- NIT (record): `N` recomputed `sum(lvl.size...)` `:416` not `self._N` — defensible (`_N` is `|None`
  unbound; partition is the output-extent authority). Leave or read `_N` if the coverage assert lands.
- NIT (do-now cosmetic): consumer comments `sweep_cache.py:311-314` + `:393-397` narrate "the matvec"/
  "slab consumer" by prose not symbol; plan §7 HAS the anchors. Recurring "line# cites → symbol cites".

## Approval conditions
None block. Recommended: (1) `np.empty`→`np.zeros` or coverage assert; (2) keep `is_linear` ClassVar;
(3) file the ABC-binding-unification issue. Runtime already green + 0-ULP bit-id per brief.
