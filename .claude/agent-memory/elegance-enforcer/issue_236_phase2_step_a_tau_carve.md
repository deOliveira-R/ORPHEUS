---
name: issue-236-phase2-step-a-tau-carve
description: #236 Ph2 Step A — relocate Morel-Montry tau production onto the angular closure; PASS-WITH-NITS, the canonical staged-carve bit-identical bridge verdict
metadata:
  type: project
---

# #236 Phase 2 Step A — τ producer relocation onto the angular closure

**Verdict: PASS-WITH-NITS** (uncommitted; dispatching agent applies). Branch `feature/sn-spatial-angular-product`.

τ = `(μ_n−μ_{n−1/2})/(μ_{n+1/2}−μ_{n−1/2})` (BMC 2010 Eq.43) is an ANGULAR-scheme property historically
PRODUCED in the streaming-GEOMETRY factory (`reduced_operator.py`) and READ by the closure. Step A mints
producer `morel_montry_tau_per_level(quad, coord)` (`pole_angular_closure.py:514-600`) + flips
`MorelMontryAngularSweep.__init__` (`:775`/`:781`/`:790`) to source τ locally instead of `reduced.tau_mm`.
BIT-IDENTICAL + ADDITIVE (geometry factory still bakes identical τ for the sweep path; Steps B/C retire it).

## The KEY ruling — temporary Pattern-2 dup is the RIGHT staging
The producer REPLICATES `reduced_operator.py:681-688` (sphere) / `:798-815` (cyl) LINE-FOR-LINE → τ now in
TWO places. ACCEPTABLE as the canonical bit-identical parallel-run-and-compare de-risk bridge (same shape as
[[w1-mm-tau-clamp-unclamp]] and the windowing `window≡full` oracle). Clears the bar because:
1. Single-sourced at the CONCEPT (BMC Eq.43); twin PINNED 0-ULP by Leg-1 gate `test_tau_producer_equivalence.py`
   to BOTH the factory AND a structurally-independent ref (`contamination.morel_montry_weights`, L11). 9/9 -O.
2. Retirement COMMITTED (anti-#11): plan `.claude/plans/issue_236_phase2_tau_carve_crosswalk.md` §1 names all
   FOUR `c_in`/`c_out` rebuild sites (P0 closure `:805-816`, P1 `cell_balance.py:313-314`, P2 `diamond.py:306-307`,
   P3 `sweep_cache.py:309-310`) + assigns Step B to retire geometry producer + consolidate; theory stub
   `discrete_ordinates.rst:990-992` re-states Step B.
3. De-risk load-bearing not ceremony (4-site blast radius; can't 0-ULP-compare against a producer you just deleted).
Marked EXEMPLARY: comment `:500-511` states replication deliberate + explains the L11 reason it does NOT call
`contamination` in prod (would collapse cross-check to tautology). Gold-standard anti-#11 temporary marker.
No more elegant staging exists (one-commit flip+retire loses the oracle).

## PASS axes
- Cartesian neutrality by TYPE not branch: producer RAISES on CARTESIAN `:596-600`; `IdentityAngularClosure`
  supplies τ=1.0 structurally `:1410-1412` (Pattern 4). `coord`-branch inside producer selects only the EDGE
  CONVENTION (weight-sum-from-−1 vs η-midpoint+clamp) BELOW the type dispatch — right structure.
- Clean Pattern-5 primitive, domain-named intermediates (`mu_edge`/`eta_edge`/`dmu`/`deta`/`tau_raw`/`sin_theta`),
  reads like BMC Eq.43.
- Carve COMPLETE for P0: grep-clean — only `reduced.tau_mm` survivors are the explanatory comment `:770-771`.
- Replication verified EXACT: same accumulation, same `1e-15` guard, same `0.5` fallback (sphere); same `±sin_theta`
  + `max(0.5,min(1.0,…))` clamp (cyl); same source data.
- test_diamond core 53✓; combined diamond+curvilinear green through 85%+ (remaining `x` = #195 sphere xfails).

## NITS (none blocking)
- (#226 follow-up) pyright `:87`/`:148`: `Quadrature` doesn't statically satisfy `AngularMeasure` Protocol
  (`reduced_operator.py:157`). PRE-EXISTING SUITE-WIDE — `test_diamond.py:368/443/539/629/841/883` call the same
  factories with `quad`, NO `# type: ignore`. +13 cascade errors (`|None` narrowing on legacy-optional `_tau_per_level`).
  DO NOT fix in isolation (would make the test LESS consistent than its 6 siblings) → route to #226 suite-wide.
- (Step B, planned) `c_in`/`c_out` algebra `:813-814` is the 5th twin of a 4-site inline; owned+tracked by Step B
  (plan §5.2 drafts `test_mm_coefficient_single_source.py`).
- (optional do-now) line-number replication cites `:502/:554/:572` → use SYMBOL cites (recurring tell; helps Step-B
  handoff verify the replica against source after line drift).
