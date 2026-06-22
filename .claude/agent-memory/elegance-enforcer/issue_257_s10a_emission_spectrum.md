---
name: issue-257-s10a-emission-spectrum
description: S10a EmissionSpectrum(ndarray) value-object + Mixture/Isotope __post_init__ χ-simplex/null guard — PASS-WITH-NITS, but the is_fissile guard predicate keys on the WRONG quantity (SigF, must be SigP/production). Branch feature/field-typed-operator-algebra, UNCOMMITTED at review.
metadata:
  type: project
---

# #257 S10a — EmissionSpectrum value-object + container χ-guard (ELEGANCE REVIEW)

Branch `feature/field-typed-operator-algebra`, UNCOMMITTED. Additive + behavior-neutral
precursor (no χ value change; that is S10b multi-fissile production-weighting). Verdict
**PASS-WITH-NITS**; one item (the guard predicate) ruled a correctness-of-intent VIOLATION
to fix BEFORE commit.

## ⭐ RULING-1 (the load-bearing one) — guard predicate must key on PRODUCTION (SigP), not SigF.
The `__post_init__` χ-guard branches on `is_fissile = np.any(SigF>0)`. WRONG quantity.
χ is consumed wherever a FISSION SOURCE is built, and that source is **χ·νΣf·φ** EVERYWHERE
— I grepped every production consumer of `.chi`: ALL pair it with `SigP` (production rate),
NONE with `SigF`. Enumerated (verify before re-trusting): `sn/fission.py:263`
`RankOneOperator(chi, sig_p)`; `billiard.py:1031-1032` (SigP+chi, never SigF);
`fn_method/moment_space.py:399-400`; `homogeneous/solver.py:104`; `cp/solver.py:500,626`;
`moc/core.py:88`; `diffusion/solver.py:271`. The ONLY `SigF` reads in prod are
`mixture.absorption_xs:87` (SigF-as-XS, correct) + `mc/solver.py:423` (`sig_a`, absorption,
correct) — neither builds a fission source, neither gates χ. So the guard depends on
"is-producing" (SigP) but READS a differently-named coincidentally-agreeing quantity (SigF)
= textbook Pattern-3 quantity-mismatch. On REAL ENDF data `SigP=νΣf` ⇒ `SigP>0 ⟺ SigF>0`
(predicate identical); the divergence ONLY surfaces on synthetic SigP-direct fixtures —
which is the entire F_N / billiard / singular-eigenfunction idiom. **Required:** `Mixture.is_fissile
= bool(np.any(SigP>0))`; `Isotope.is_fissile = bool(np.any(nubar*sigF>0))` (the bool() wrap
ALSO clears the pre-existing `isotope.py:71` np.bool_ pyright red — mirror the Mixture wrap
per the implementer's own LESSON). STANDING TELL: when a guard/dispatch keys on a field,
ask "what quantity does the guarded behaviour actually DEPEND on?" — grep the consumers of
the guarded value; if they read a DIFFERENT field than the guard, the guard is the wrong
spelling (coincidence on real data masks it).

## ⭐ RULING-2 — the billiard SigF hack is a symptom, REVERT it.
`tests/derivations/test_trajectory_resolvent_billiard.py:105-114` `_mixture_from_xs` set
`SigF=nu_sf_arr.copy()` (dimensionally WRONG — conflates Σf with νΣf) + a 4-line apology
comment, purely to flip the SigF-keyed `is_fissile`→True so χ=[1,0] passes. anti-pattern
item-9 (code you must write a paragraph to excuse = the smell). With RULING-1's SigP-keyed
predicate, revert to `SigF=np.zeros(ng)` + delete the comment: `SigP>0` classifies the
multiplying medium correctly, χ stays the honest simplex, fixture is dimensionally honest.
Bug habitat if kept: `_mixture_from_xs` is the obvious copy-template for "build a multiplying
medium" → future contributor inherits `SigF==SigP` → wrong `absorption_xs`. Re-run DD byte-id
regression + 12 billiard tests after the flip+revert to confirm neutrality.

## NDARRAY-SUBCLASS vehicle — PASS (right call, gotcha-clean).
`EmissionSpectrum(np.ndarray)`, `__new__`=`np.asarray(values).view(cls)`,
`__array_finalize__` deliberate NO-OP. RIGHT over a frozen-dataclass-wrapper: a wrapper
forces `.values[None,:,:]` ripple at ~10 broadcast/einsum consumers for a guarantee χ
DOESN'T need (χ is NOT in the flux algebra — explicitly fenced from the #263 SpectrumField
in the module docstring). Subclass carries the law as methods (= operator-carries-.H analogy,
user's goal). GOTCHA AUDIT CLEAN: the no-op finalize is SAFE precisely because the law lives
in EXPLICITLY-INVOKED validators, not carried state — verified `assert_normalized`/`assert_null`
fire at EXACTLY 4 prod sites (`isotope.py:57,59`+`mixture.py:70,72`), all on the fresh
full-vector χ in `__post_init__`, NEVER on a slice/ufunc-result. A ufunc producing a
non-simplex EmissionSpectrum is harmless (nobody validates it). STANDING TELL for ndarray-
subclass value-objects: the no-op `__array_finalize__` is correct IFF validators are invoked-
not-automatic AND never called on a slice/arithmetic result — grep the validator call sites
to confirm they only hit the source full-vector.

## NIT (CONCERN, do-now) — `__post_init__` is a 2-instance twin → extract NOW.
The coerce-then-`if is_fissile: assert_normalized else assert_null` block is near-identical
in `isotope.py:50-59` + `mixture.py:63-72` = Cardinal-2 shared-concept-in-2-places. The
PREDICATE genuinely differs (SigF/SigP vs nubar·sigF) but the ENFORCEMENT (coerce+branch)
is identical and is the drift-prone part. Pattern-6 says defer-until-two — two instances are
IN THIS DIFF, so the trigger is MET (extract, don't defer). Remedy: free func in
emission_spectrum.py `enforce_emission_spectrum(chi, *, is_producing:bool)->EmissionSpectrum`
taking the ALREADY-DECIDED predicate as a value → each carrier keeps its own predicate
(polymorphic part) while the coerce+branch law single-sources next to the type.

## NIT (CONCERN) — `is_emitting` is dead weight + a THIRD spelling.
`emission_spectrum.py:107-110` `is_emitting` has ZERO consumers tree-wide (grep: only its
own unit test). Mild Pattern-6/item-11 (API ahead of 2nd instance). Real hazard = it's a
THIRD spelling of "is-producing" (`np.any(chi>0)`) alongside `is_fissile`(SigF) + the true
predicate(SigP) — three spellings of one concept, the Pattern-2 divergence in miniature
(χ=[1,0] is_emitting=True but SigF=0 ⇒ current is_fissile=False — the inconsistency RULING-1
fixes). Remedy: drop until a real consumer (lean), or doc it spectrum-level-≠-material-level.

## PASS — naming + 2-clause independence.
`assert_normalized`/`assert_null` read like the domain (normalized prob simplex / null
spectrum). 2-clause independence in assert_normalized (sum-clause ⊥ negativity-clause,
`[1.2,-0.2]` reds on negativity despite Σ=1) reads the intent literally. `assert_` prefix
matches the `_require` idiom.

## Architectural opportunity (filed-issue worthy, module:data type:improvement)
THREE circulating spellings of "neutron-producing material": `SigF>0`, `SigP>0`,
`chi.is_emitting`. Coincide on real ENDF, diverge on synthetic. Cardinal-2 concept-in-3-
places. Recommend an issue to single-source a canonical `is_producing`/`is_multiplying`
predicate on SigP across all fission-source sites + retire the SigF mis-spelling. The
implementer's closeout already documents the seam (the "new seam" LESSON) — the issue
captures it per Cardinal-4.

STANDING TELLS distilled: (1) guard-keys-on-field-X but consumers-read-field-Y ⇒ wrong
spelling masked by real-data coincidence — grep consumers to find the quantity the behaviour
DEPENDS on; (2) a hack that fakes a field to satisfy a guard + carries an apology comment ⇒
the GUARD is wrong, fix the predicate and delete the hack; (3) no-op `__array_finalize__`
on an ndarray-subclass value-object is safe IFF validators are invoked-not-automatic and
never hit a slice — confirm via grep of the validator call sites; (4) 2 instances in one
diff = the extract trigger, not the defer trigger.

Extends [[issue-257-s5-functional-category]] (same campaign, typed-field-algebra family).
