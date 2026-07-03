# Elegance Enforcer — Lessons

Review-PROCESS corrections only: "what did I mis-judge / catch, and what
sharpened my verdict?" The pattern/anti-pattern CATALOG lives in the preloaded
`coding-elegance` skill — never restate it here; cite it (`Pattern N`, `anti-#N`).
The cross-cutting institutional smells (twin-delivery plumbing, role-grid,
fuller-view-oracle exception, the recurring tells-to-grep) live in AGENT.md
§"Institutional knowledge" — do NOT duplicate them; these lessons are HOW I
verify a verdict, not WHAT the smells are.

The spine behind most of these — the three-leg VIOLATION standard (bug-habitat
edit named / coextensive-today→NIT / verify-the-live-tree-not-the-docstring) — is
**now a standing directive in AGENT.md** (§"The VIOLATION standard — three legs").
It is no longer recalled-when-relevant; it governs every verdict by definition.
Each lesson below is one face of that standard, kept for its forensic specifics.

---

## L-001 -- The single most recurring catch: "docstring names a primitive the code does not call"

Across nearly every typed-field / moment-layout carve I reviewed, the dominant
do-now NIT was the same shape: a docstring (or rationale comment) asserts the code
uses a shared single-source primitive (`is_moment_valued_by_rank`, `_n_face_moments`,
`face_moment_tail`), but the body OPEN-CODES the predicate as a fresh spelling
(`bulk.ndim == flat_ndim`, `per_axis ** (ndim-1)`). The two are provably
coextensive today, so it is a NIT not a VIOLATION — but the docstring is already
LYING, and a future maintainer greps for the named call, finds none, and is misled
about where the single source lives.

How to apply: when a docstring cites a primitive, GREP for the call. If absent,
the code open-codes a twin (anti-#1) AND the doc misstates an invariant (anti-#11)
— flag BOTH, and note they auto-resolve together when the code routes through the
primitive. The fix is almost always "add a shape-only overload of the primitive
(`is_moment_valued_by_flat_rank(arr, flat_ndim)`) so the consumer-with-no-reference-
array can join the single source." (#247 LegA NIT-1/2, #251 LegB NIT-2, #246, Step C.)

---

## L-002 -- A recurring inline predicate is a single-source DEBT, but the collapse trigger is the THIRD DIVERGENT consumer, not the second

The "X.ndim > ref.ndim + 1" / "per_axis^d" / "moment_tail != ()" family kept
re-appearing inline across `loss_representation` + `sweep_graph` + `geometry`. The
correct verdict is NOT "collapse on sight" — it is ACCEPTABLE-FOR-NOW when every
spelling is non-divergent (same convention, provably coextensive) AND I have
byte-verified the overlapping arithmetic. It STOPS being acceptable the moment a
THIRD consumer carrying a DIFFERENT rank lands (a `2^{d-1}` face variant, a moment
source with 2 extra axes). The do-now remedy at the second spelling is reciprocal
cross-ref comments + a tracked removal trigger naming the third-divergent-edit —
NOT premature unification (which is its own anti-#10).

How to apply: when you see the N-th inline copy of a layout predicate, count
DIVERGENT consumers, not raw copies. Co-locate the eventual primitive in the
convention-owning leaf (`numerics/moment_layout.py`, beside `face_moment_tail`).
State the collapse trigger explicitly in the verdict so a fresh session knows when
to pull it. (#246 ruling 2, #251 Arch-Opp, #240 D5b-S3 CONCERN.)

---

## L-003 -- "Two parallel predicates" is NOT a unify-trigger when they speak different return types to different consumers

Twice I was asked "these two trait/capability gates look parallel — co-locate them?"
and the correct answer was NO. Discriminator: do they return the SAME type to the
SAME kind of consumer? A `bool` physics-validity query (diffusion-limit pairing)
and a `Compatibility(ok, reason)` sweep-strategy-DISPATCHER query are parallel
IN SPIRIT (both per-axis-trait conjunctions) but live on DIFFERENT axes
(scheme⊗closure vs scheme×geometry) and speak DIFFERENT return types to DIFFERENT
layers. Co-locating CONFLATES two concepts — a real concept-merge, the opposite of
a single-source win. Different return types argue AGAINST premature unification, not
for it. Same logic retired the wrong instinct to extract a `_blend(require_convex:
bool)` over affine-combination vs convex-mix (anti-#3 flag dispatch + cross-role-
family dep) — the CONSTRAINT defines the operation's identity, so they are not one op.

How to apply: before recommending "unify these two," check return type + consumer
layer + whether the distinguishing axis is a genuine identity-change. If unifying
forces a boolean flag, a cross-family dependency, or a `Compatibility`-vs-`bool`
coercion, it is a concept-merge — KEEP them separate and SAY why. (#236 Ph1b ruling,
#257 S1 J4.)

---

## L-004 -- The stale-doc blast radius lands OUTSIDE the diff

→ **Now a standing directive in AGENT.md** (§"Scope of Review", the "diff boundary
is not the review boundary for a deletion/migration" paragraph). The general rule
(grep the whole tree for the deleted symbol AND its pre-deletion line numbers;
discriminate hits by tense) is now identity-level. Forensic discriminator retained
here for the line-cite case the directive compresses:
- DELETED LINE-NUMBER cross-reference ("Replicates X (file.py:681-688) EXACTLY") =
  MUST-FIX (asserts a live twin that no longer exists; prefer symbol-refs to line-cites).

Same blast-radius logic the method-implementer learned for ORACLE-READERS hiding
outside the audit's named files — it applies to DOC refs too. (#236 Step C standing
lesson; #222 S6.4(d) stale WavefrontFlux docstrings.)

---

## L-005 -- Prove the verification TEETH bite before crediting a re-baseline as non-laundering

The most decisive thing I do on a carve that re-baselines a snapshot or pins a
value: confirm the GATE that's supposed to catch a regression actually reddens under
the regression. The verdict "the re-baseline did not launder a bug" is earned only by:
1. The pin targets an INDEPENDENTLY-constructed structural reference, not the SUT's
   own output. `(L+C).apply == loss_action(σ)` byte-exact + `L.apply ≈ loss_action(σ)
   − C·ψ` against a SEPARATELY-built `CollisionOperator` is the gold standard; a
   snapshot regenerated by `np.save(out.values)` is a laundering habitat UNLESS an
   independent structural pin runs every session (verify it exists — it usually does
   elsewhere in the suite, which DOWNGRADES the snapshot smell from CONCERN to NIT).
2. A MUTATION/DISABLE proof: monkeypatch the term to re-introduce the bug (a σ-re-
   reading stub, a rename `apply`→`_DISABLED_apply` to let the leaf-sum take over) and
   confirm the teeth FAIL at the predicted ULP. For an affine-in-σ "value-correct-by-
   coincidence" twin, the teeth gate MUST be `array_equal` (0 ULP) — only bit-identity
   discriminates leak-from-override (both are value-equal).

How to apply: never accept "the gate is green so the re-baseline is safe." Ask what
the gate is pinned AGAINST (independent? or self?) and whether a mutation reddens it.
This is the same discipline as the method-implementer/qa "teeth bite" lessons — I
apply it to ADJUDICATE the carve's own verification, not just to write tests.
(#257 S8b rulings 1/3, #240 Step B.)

---

## L-006 -- A sub-floor (Mode-10) term is closed by STRUCTURAL teeth, NOT a tightened value band — do not flag the absence of a value-improvement leg

A term whose code path runs but whose error is O(h²)-small (a slope-source, a
boundary-trace transverse slope) is sub-floor: tightening the converged-flux value
band does NOT pin it. The complete closure is a producer-threads-at-machine-precision
proof + a consumed-source-row-sign-flip-moves-the-answer-≫-tol proof + a flat no-op
control pinning the asymmetry. A SHARPER case has NO value-improvement leg AT ALL
(a correctly-consumed slope can even make the converged value slightly WORSE) — and
there I must NOT demand a value-band gate, because manufacturing one would falsely RED
a correct term.

How to apply: when reviewing teeth for a localized/higher-order-small term, ACCEPT
structural teeth + no-op control as the whole proof; do not flag "but the value didn't
improve" as a coverage gap — for a sub-floor term that is the CORRECT signature.
(#247 LegA, #251 LegB, #257 S9; mirror of method-implementer L-007 / vv Mode-10.)

---

## L-007 -- Pattern-6 TRIM a named predicate whose OWN docstring says production will never use it

A clean recurring catch: a new named typed-query predicate (`has_spatial_moment_axis`)
minted ALONGSIDE a test pin, with ZERO production consumers (grep the whole tree), and
a docstring that explicitly says the inner walk MUST use a DIFFERENT predicate. That is
textbook over-abstraction (anti-#10 / Pattern-6 ≥2-instances): an abstraction whose own
documentation explains why prod will never call it. The invariant it was minted to pin
almost always has a PRE-EXISTING lower-level property to assert against
(`spatial_moments_per_axis == 1`) — find it and retarget the tests onto it.

How to apply: for any new predicate/Protocol/ABC in a diff, grep for production
consumers FIRST. Zero prod consumers + a docstring conceding prod won't use it =
TRIM-on-sight, retarget the tests to the pre-existing property. KEEP only if a real
production method-boundary consumer lands in the SAME commit (rule-of-two). Contrast:
a dunder-EMPTY role-marker mixin with its 2nd consumer IN THE SAME DIFF is NOT a TRIM —
the rule-of-two IS satisfied at landing, and a long docstring on an empty body is
ELEGANT there (the content is "absence of a gate," which cannot self-document via code).
(#246 ruling 1; contrast #257 S1 J1.)

---

## L-008 -- Prefer the emergent invariant-gate over a hand-written dunder when a construction invariant already exists

When a type has a construction invariant (`SpectrumField` simplex Σχ=2 enforced in
`__post_init__`), the elegant way to reject `χ+χ` is to let the inherited `__add__`
→ `replace()` on the frozen dataclass RE-RUN `__post_init__` (which raises) — NOT to
hand-write an explicit `__add__` override. The explicit override is a SECOND statement
of the same law (anti-#1: relax the tol in one, the other goes stale). The emergent
gate single-sources the invariant. (Contrast FluxRole, which MUST hand-write
`flux+flux→TypeError` because flux-addition is type-coherent-but-void — no construction
invariant catches it, so there is no emergent gate to lean on.)

How to apply: when reviewing a "should this type forbid operation X?" question, check
whether a construction invariant already exists that `replace()` will re-run. If yes,
the emergent gate is more elegant and I should flag an explicit dunder override as a
Pattern-2 duplication. If no invariant exists, the explicit gate is correct.
Caveat to flag: `@runtime_checkable` Protocol conformance (`isinstance(χ, Vector)` True)
is method-PRESENCE only — a gated-arithmetic leaf passes isinstance yet raises in a
generic accumulator loop; add a clarifying comment xref-ing the gate test. (#257 S1
J3/NIT-2.)

---

## L-009 -- The "honest typing" verdict (overload / TYPE_CHECKING) is adjudicated by a stash-baseline pyright DELTA + a runtime alias proof, not by error count

For carves that make a heteromorphic dispatch honestly typed (the first `@overload`-
on-`singledispatch` "Pattern M" precedent), I ratified the design by:
1. STASH-TO-BASELINE pyright delta: `git stash push -- <the N files>`, count errors by
   (file,line), `git stash pop`. The verdict needs NET errors AND zero-new — and
   because line numbers SHIFT under a diff, categorize the survivors SEMANTICALLY (read
   each surviving line) rather than by line number. "All errors pre-existing" is a
   claim to VERIFY this way, never to take on faith.
2. RUNTIME ALIAS proof: `Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']` (True =
   same descriptor, single source). Do NOT use `Cls.apply is Cls._apply_impl` — the
   descriptor `__get__` makes fresh bound wrappers per access (False is a test artifact,
   not a defect).
3. The design ranking is Beck's: "reveals intention" (rule 2) OUTRANKS "fewest elements"
   (rule 4). Pattern M (typing epilogue under `if TYPE_CHECKING:`, math at natural
   indentation) beats Pattern C (math buried in a double-indented `else:`) BECAUSE
   burying ~150 lines of real algebra to hide typing scaffold inverts the hierarchy.

How to apply: for any typing-honesty carve, run the stash-baseline pyright delta and
the `__dict__`-identity alias check; rank competing spellings by intention-revelation
before element-count. `cast()` is HONEST (not laundering) when the target is the EXACT
runtime invariant and a wrong carrier hits a fail-loud fallback. (#257 S8c.)

---

## L-010 -- When reviewing a sibling-repo (Nexus) tooling change, the latent crash is the MISSING guard a sibling method already has

On the sphinxcontrib-nexus overlay reviews, the highest-value findings were symmetry
breaks: one overlay method omits the `if node_id in self._g` stale-node guard that BOTH
its sibling methods carry — and that exact path (a node renamed by a rebuild between
ingest and query) is the sidecar design's whole reason to exist. The crash is silent
(`degree()` on a missing MultiDiGraph node returns a view, not an int, does not raise →
`TypeError: not JSON serializable` at the MCP boundary, violating "tools never raise out").

How to apply: on a tooling/library review, when N methods do "the same kind of thing,"
diff their guard clauses — a guard present in N-1 and absent in 1 is the latent bug
(symmetry-in-code). Also flag the self-disproving docstring ("exactly one metric family
populated" refuted by the same commit's heterogeneous-union `merge_runs`) — a false-
invariant doc IS the bug habitat (anti-#11), and "it's just a docstring" makes it more
insidious, not less, because it is the contract the next agent reads. (nexus #26 overlay.)

---

## L-011 -- The briefing's SCOPE is a claim to verify against live `git status`, not authority — a stale snapshot of WHICH files changed manufactures false "missing-gate" findings

Leg-3 of the VIOLATION standard ("verify the live tree, not the docstring") extends to
the DIFF SCOPE itself. On #226 inverse-as-operator the dispatch brief (and the harness
session-start `git status` snapshot it echoed) said "only production + ONE test file
changed." Trusting that, I built a near-certain SHOULD-FIX: "8 inverse types delivered,
0 round-trip/involution gates — violates §12 step-1's own 'pinned by the round-trip
gate'." Then I re-ran `git status` LIVE: three test files modified + an UNTRACKED
`test_inverse_universal.py` — the full §12 keystone suite (I1 both-ways + closed-form
anchor + `inverse().apply≡solve` bit-id + I2 per-type + negatives + the shim), with
deliberate mutation-teeth fixtures. The "ungated" finding was pure snapshot-staleness;
I retracted it and it became a PASS callout.

How to apply: at scope-identification, run `git status --short` + `git diff --stat HEAD`
FRESH — never trust the harness's session-start snapshot (it is frozen at launch) nor the
briefing's file list. Untracked test files (`??`) do NOT appear in `git diff` and are the
easiest coverage to miss — grep the working tree for the tested SYMBOLS (`grep -rn
InverseOperator tests/`), not just the diff. A "promise not backed by a test" verdict is
only earned after enumerating the LIVE test tree, including untracked files. Same family
as trust-git-not-a-frozen-claim (process-discipline) but for change-scope, not
merge-status. (#226 step-1 review, 2026-07-01.)

---

## L-012 -- To PROVE a `# type: ignore` is dead, enable `reportUnnecessaryTypeIgnoreComment` against the LIVE tree — the error-count ratchet is BLIND to it

A pyright retirement carve that tightens an Optional-away (`inflow_full: NDArray|None →
NDArray` by moving its binding into a narrowed scope) makes the site's `# type: ignore[index]`
UNNECESSARY — but the campaign's error-count ratchet (`pyright_ratchet`) counts errors, and
`reportUnnecessaryTypeIgnoreComment` defaults to OFF, so a dead ignore is INVISIBLE to the
ratchet and survives silently. On C5 the brief AND the dev-history doc both claimed "two
stale `type: ignore` retired"; the LIVE tree had only ONE removed (`quad.level_indices`) and
TWO `inflow_full[...]  # type: ignore[index]` survivors (loss_representation:3359,3389) that
the restructure had just made dead. Leg-3 catch (verify live tree, not the campaign's
self-description of what it retired).

How to apply: when a carve narrows a type such that an ignore SHOULD have become removable,
do NOT eyeball it — write a throwaway config at the repo ROOT (so `venvPath`/`ignore` resolve)
mirroring `pyrightconfig.json` + `"reportUnnecessaryTypeIgnoreComment": "error"` +
`"include": ["<the one file>"]`, run `npx pyright -p <tmp>.json`, grep the lines, then `rm`
the config. Read-only — never edit the file under review. Bug habitat = anti-#19: a dead
`[index]` ignore silences `reportIndexIssue` at that exact site, so the day a future edit
re-introduces an Optional there the real error is swallowed. Coextensive-today ⟹ NIT, but a
do-now NIT: the campaign's whole purpose is removing dead suppressions, and the doc that
CLAIMS they were retired is itself now false. (#226 C5 pyright-burndown review, 2026-07-03.)
