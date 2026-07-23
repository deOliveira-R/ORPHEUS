# QA Agent Memory — index

## 1. Lessons (the behavioral spine)

[lessons.md](lessons.md) — 56 behavioral lessons (`## L-NNN -- title`,
ascending, contiguous L-001..L-056). "What mistake did I make, what did I
learn that improved my review behaviour." Consult before every review;
sharpen in place after every task. The recurring spine:

- **Coverage ≠ count** (L-001, L-005, L-030) — demand a heterogeneous /
  multi-group / mesh-refinement gate; a per-axis invariant needs a fixture
  that VARIES along the non-reduced axis.
- **Mode-8 (-O strip)** (L-006, L-010, L-039, L-052) — bare asserts fire in
  collected `tests/` (rewriter), die only in `orpheus/` AND in your own probe
  SCRIPTS (L-052: a bare assert in a `python -O` throwaway prints PASSED while
  unequal — teeth-check through pytest or `np.testing.assert_*`, never a bare
  script assert); prove a gate fired by mutating it, not by inspection.
- **catches/verifies markers** (L-007) — a marker is a coverage CLAIM;
  mutation-verify the EXACT documented bug reddens THIS test or drop it.
- **Re-baseline integrity** (L-022, L-023, L-024, L-028, L-044, L-049) —
  grep the WHOLE tree for the old literal; pin a regen'd `.npy` to a
  STRUCTURALLY-INDEPENDENT reference, never old-vs-new ULP; prove the
  masking-check (old snapshot still hard-fails the new code).
- **Twin-path / Mode-11** (L-018, L-021, L-031, L-033, L-036, L-043) — a
  green "twin" must EXECUTE the rewired line; sentinel-instrument it,
  mutate the SHARED source not the dead-for-this-path method.
- **Circularity & structural independence** (L-026, L-029, L-046) — a
  re-encoded production formula is a circular VALUE check unless the other
  side is independently assembled; a weighted pin needs every factor in the
  hand-ref.
- **Bit-id micro-facts** (L-013, L-014, L-020) — resolve byte-id disputes
  with the IEEE/Python fact + git `.npy` status + a +1-ULP perturb, not a
  docstring.
- **API-migration rewire bit-id + phantom-symbol discipline** (L-051) — a
  "deleted class → new face" rewire's bit-id is VERIFIABLE (recompute the OLD
  einsum on a structurally-independent table, `np.array_equal`), never ASSUMED
  from the brief; confirm every brief-named symbol/file with find/grep (two
  phantoms caught) and byte-compile no-test generator scripts.
- **Adjoint-via-metric-composition** (L-052) — `A.H = G_dom⁻¹·Aᵀ·G_cod` is
  VERIFIABLE: re-derive the inner-product identity by hand, then prove it with a
  DENSE matrix built by loops + transposed directly + composed with metrics
  (zero shared code with the fused einsums); the "weight-free transpose" choice
  is provable (each transpose mirrors its OWN forward — asymmetric w/ a
  weight-baking sibling is correct, not a missing factor); teeth = drop-factor /
  bake-weight / reverse-factor / wrong-Gram-power, all RED via monkeypatch;
  per-term fold beats post-scale for a bit-faithful ref (0 vs 112 ULP).
- **Stale-snapshot triage** (L-034) — a HUGE-ULP red on ONE geometry while
  siblings pass = live-correct / frozen-stale; find the apply-changing
  commit that didn't re-capture.
- **Mode-10 honest-scope** (L-026, L-037, L-038) — an exercised term whose
  sign error is sub-floor is UNVERIFIED, not verified; structural teeth +
  no-op control when no dominant regime exists.
- **Protocol/category gates** (L-039, L-047) — runtime_checkable only checks
  member PRESENCE; the direct `not hasattr` negatives are the defense.
- **Behavior-neutral retype** (L-041, L-045, L-048, L-050) — re-prove
  inertness for EVERY consumer with a direct value comparison; a "neutral"
  claim holds only for the ONE contract it was proven against (ERR-063).
- **Doc-staleness adjudication** (L-055) — campaign-narration (`Phase X`/
  `Wave Y`) FIX bar = "provably lies about CURRENT code", verified via
  grep(symbol/wiring/workaround)+`gh issue view` BEFORE ruling; guards:
  a stale RUNTIME STRING is behavioral→KEEP, a "HALTs Phase X" load-bearing
  banner is a characterization record→KEEP.
- **Distillation-fidelity review** (L-056) — a skill→Sphinx doctrine page:
  the DOCTRINE is faithful (read vs preloaded skill), yield is in code-anchored
  SPECIFICS the build can't gate — `:mod:`/`:class:` roles are NOT `-W`-gated
  (dead target renders as plain text; grep corpus, the OUTLIER-count spelling
  is the bug), and the skill SOURCE carries stale specifics that propagate, so
  verify module-path/`mpmath`-vs-`scipy`/test-count against CODE + the consuming
  test's docstring, never the skill twin.

## 2. Active / in-flight state

None. All SN review campaigns whose lessons are recorded above are merged to
`main` (verified via `git merge-base --is-ancestor`). Only #236
(`feature/sn-spatial-angular-product`) remains open at the repo level, but it
carries no unresolved QA finding here. Git is authoritative for merge status —
re-verify before acting on any "in-flight" claim.

## 3. Durable reference (topic files)

- [field_role_typing_apply_sourcesink_contract.md](field_role_typing_apply_sourcesink_contract.md)
  — re-checkable SN role contract (`.apply`=AngularSourceSink,
  `.solve`=AngularFlux) + the A2D-1 source-hash-pin update procedure +
  affine-gate test-migration playbook. Cited by `qa/AGENT.md` enforcement
  rule #10 — **durable**.
- [phase1_moment_space_review.md](phase1_moment_space_review.md) — the
  ERR-039 moment-space P1.1–P1.7 verification-of-record; cited as the
  verification artifact by 3 plan/agent files — **durable**.
- [issue_247_legA_review.md](issue_247_legA_review.md) — full #247 Leg A
  (slope-source) review; distilled into L-037. **Retire candidate** (merged
  campaign; the reusable behavior is in the lesson).
- [issue_251_legB_review.md](issue_251_legB_review.md) — full #251 Leg B
  (boundary face-slope) review; distilled into L-038. **Retire candidate**
  (merged campaign; the reusable behavior is in the lesson).
