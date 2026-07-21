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
`"include": ["<the one file>"]`, run `/opt/homebrew/bin/pyright -p <tmp>.json` (repo-root path,
NOT `/tmp` — venvPath `.` resolves relative to the config's dir), grep the lines, then `rm`
the config. Read-only — never edit the file under review. Bug habitat = anti-#19: a dead
`[index]` ignore silences `reportIndexIssue` at that exact site, so the day a future edit
re-introduces an Optional there the real error is swallowed. Coextensive-today ⟹ NIT, but a
do-now NIT: the campaign's whole purpose is removing dead suppressions, and the doc that
CLAIMS they were retired is itself now false. (#226 C5 pyright-burndown review, 2026-07-03.)

**Sharpening (SN un-weld 1b, 2026-07-07):** the check applies to MOVED ignores, not just
NEW-on-narrowed ones. When a carve relocates a method to a new sibling class (here
`_reflect_radial_characteristic` → `RadialCharacteristicBoundaryOperator._reflect_corner`),
a `mesh=seed.mesh  # type: ignore[attr-defined]` it carries along can be dead — and was
ALREADY dead pre-carve (the field declared `mesh: SNMesh` all along; only 2 of 7
`RadialCharacteristicSourceSink(...)` sites carried the ignore, the tell). The carve's
"pyright transport:1 (+0)" gate was TRUE yet HOLLOW at both the moved copy AND the fresh
duplicate the un-weld minted. Precedent-setting carves make it worse: the dead ignore
becomes the template every future block operator copies. Run the check on any file that
gains/moves a `type: ignore`, not only on the Optional→required sites.

## L-013 -- A doc gather/split (monolith → chapter) is a RETIREMENT review; the recipe is mechanical and the spelling axis is project-specific

Reviewed the first SN mini-book chapter carve (`slab_one_group.rst` gathered from 8
regions of the 16k-line `index.rst`; branch `docs/sn-doc-architecture`, more chapters
coming — slab_multigroup / cartesian_multid / curvilinear_*). A doc extraction is a
`coding-elegance` Pattern-2 (single-source) + aggressive-retirement review wearing prose,
and the verification is fully mechanical — do NOT eyeball 700 lines:

1. **Twin content = label uniqueness + distinctive-prose grep.** Sphinx `-W` already
   guarantees no duplicate `:label:`/`.. _anchor:` across files (a genuine twin would red
   the gate the brief says is green) — so a moved label showing in BOTH files is either a
   real duplicate OR a `\b`-word-boundary false positive (`transport-cartesian` matched
   `transport-cartesian-2d`). Resolve by exact `grep -n`. For UNLABELED prose/eq (invisible
   to `-W`), grep distinctive fragments and assert count-in-chapter vs count-in-index is
   never "both".
2. **Math fidelity = extract-from-git-HEAD + char-diff, never by eye.** `git show
   HEAD:…/index.rst`, parse each `:label:` body out of both files with a tiny python
   extractor, assert byte-equality. Caught nothing wrong here (9/9 identical) but that IS
   the point — the char-diff is what lets me CERTIFY "verbatim travel", which eyeballing
   cannot.
3. **Account for EVERY deleted hunk** against {traveled | retired-as-named | stays}. The
   brief's retirement list is a CLAIM: grep the retired section's present-tense argument
   is ABSENT from the chapter (e.g. "Why not extract apply_transpose analytically" pro-dense
   argument — gone; the past-tense "that path retired in D-K" note is the KEEP, don't confuse
   them). A rewritten section (Krylov H1 → "The Krylov alternative") is faithful iff the
   stale story ("Consistency with the Sweep", "Round 2 deviation") is gone from BOTH files.
4. **The spelling axis is the project-specific discriminator (reusable every chapter):**
   the honest algebra is `A = L + C − S − B` (A = the FULL within-group operator); the sweep
   is `(L+C)^{-1}`, "the inner kernel of the full inverse, **never** the full inverse itself"
   — so **`A = L+C` and `A^{-1}`-for-the-sweep are the stale spelling to flag**. Scope: only
   NEW text (the chapter + inserted pointer paragraphs + rewritten blocks) must harmonize;
   RETAINED monolith sections keep the old local convention until their own chapter (brief
   said so — don't flag them). `A` as a generic bound variable in the reciprocity identity
   `⟨Aψ,φ⟩=⟨ψ,A*φ⟩` is fine (and char-identity-required). `A_wg^{-1}` inside a *refuted*
   approach is a defensible NIT, not a MUST-FIX.
5. **Diff boundary ≠ review boundary (L-004 again):** the retired blend-weight/Péclet
   spectrum's *correct home* (`discretization.rst`, the cross-link TARGET, outside the diff)
   had to be opened to confirm (a) the chapter only cross-links, doesn't restate, and (b) the
   monolith's retired version was genuinely WRONG (`k→∞,w→0` vs the code's `k=(g/θ)/(g/θ+Σt)
   ∈[0,1)`, `w∈[½,1]`). Also verify a CHANGED `:doc:` hint points where the `:ref:` actually
   lives — here the chapter *improved* on the source (`bc-extraction` hint corrected
   operator_algebra→boundary_conditions). Clean carves earn the credit; say so.

**Sharpening (slab_multigroup chapter 2, 2026-07-20):** when a gather/split ALSO declares a
convention-**harmonization** pass (here: retire local `A = L+C` → honest `A = L+C−S−B`, move
`F` to the eigenvalue RHS), the highest-value residue is a **stale-convention survivor in an
inline `:math:` inside a REWORDED-prose sentence.** The `.. math::`-block set-diff (extract all
`.. math::` bodies from HEAD-windows vs chapter, compare as sets) is necessary but BLIND to it:
the survivor is an *inline* `:math:` (not a block), so it lands in neither the "removed" nor the
"traveled-verbatim" set — it just rides along in prose. Prose-rewording and math-harmonization
are SEPARATE edits; the author reworded the KrylovAcceleration bullet (`inverter`→`preconditioner`,
"A⁻¹ step"→"step inverse") but left the embedded `M ≈ (A − S − F)^{-1}` in the retired convention
(`A=L+C` ⟹ `L+C−S−F`), which then CONTRADICTS the chapter's own harmonized posing `L+C−S−B`
(−F where the honest operator has −B; fission is cross-group, never a within-group loss term). The
L-013 point-4 targeted honest-algebra GREP (`A ?- ?S`, `\sum(A`, `A\psi`, `A = L \+ C` as identity)
is what catches it — run it over the WHOLE chapter, inline math included, NOT just the block
set-diff. Verdict = MUST-FIX (a math/physics contradiction in NEW scientific doc, Cardinal Rule 1),
but a one-line surgical drop of the stale term. Everything else (5/5 labeled eqs byte-identical,
13 unlabeled displays verbatim, all hunks accounted incl. two boundary-context headers that STAY
in the monolith, catalog compression + section retirement clean, both forward-promises landed) was
exemplary — a single inline-math harmonization miss in an otherwise clean 933-line carve.

**Sharpening 2 — the traveled-verbatim-vs-fresh discriminator decides MUST-FIX vs CONCERN on a
dual-`A` overload (cartesian_multid ch.3 Stage-1, 2026-07-21).** Same dual-meaning-of-`A` smell
recurred (fresh intro `A = L+C−S−B` juxtaposed with a §15A.2 display `A^{-1} = (L+C)^{-1}` that
only parses if `A = L+C`). Last chapter I called the analogue a MUST-FIX; here CONCERN. The
discriminator that flips the verdict is **did the author TOUCH the stale text, or splice it
byte-exact?** — settled by extracting the paragraph around the label from `git show HEAD` and
char-diffing. slab_multigroup:864 was FRESH-EDITED prose (author reworded the sentence, left the
math) + a real physics error (fission in a within-group operator) + in the declared harmonization
scope ⟹ MUST-FIX. cartesian_multid §15A.2 traveled BYTE-EXACT (verbatim splice, author never
touched it), is LOCALLY correct (`A=L+C` bound in-section by "the loss inverse (the sweep)"), and
the brief's declared staging method explicitly scopes traveled content as verbatim-verified /
harmonization-DEFERRED ⟹ CONCERN, deferrable to the Stage-2 harmonization of that section. Rule:
when a staged carve declares "traveled = byte-exact, fresh = harmonized," hold ONLY fresh content
to the honest-algebra standard; a stale spelling in a byte-exact span is out-of-scope-this-stage
(flag as forward-looking CONCERN, not MUST-FIX) — UNLESS the author edited it, which pulls it back
into scope. The fresh intro that NEWLY introduces the colliding binding is what makes it worth a
CONCERN at all. Also new this review: the scripted byte-exact splice made the WHOLE math layer
trivially verifiable — 4/4 labeled + 8/8 blocks verbatim, 0 fresh math blocks — so the review
collapses to (1) confirm the 3 declared transforms are the ONLY diffs [underline promotion, one
"above"→"" fix, W1/W2 title re-authoring], (2) fresh-prose honesty + dispatch-facts-vs-owning-page,
(3) label/hunk accounting. Verdict COMMIT-READY, zero MUST-FIX.

**Sharpening 3 — when the byte-exact core certifies clean, the ONLY live finding is the L-004
EXTERNAL-repoint blast radius; a doc-split's transform-5 that sweeps CODE-file qualifiers but
misses the TEST-file twins (cartesian_multid ch.3 Stage-2, 2026-07-21).** The two-span move
(W3a UBLD → cartesian_multid, W3b 2-D LD stress MMS → verification) certified PERFECTLY:
promotion counts 30/24 exact, underline census exactly 3 levels (no un-mapped 4th → no nesting
hazard), 31/31 labeled eqs char-identical, W3a = 2 declared annotations only + W3b = 0 drift,
section hierarchy nests with ZERO level jumps (the chapter's FLAT multi-`=` top level, set in
Stage-1, is why `-`→`=` correctly makes the span a top-level SIBLING — verify the promotion is
STRUCTURALLY right, not just count-right, by extracting the whole header tree and checking each
child sits exactly one level under its parent). Annotations truthful (every claimed code symbol
grep-confirmed: `spatial_moments=` threading, `D1ClosedForm.*`, `_moment_broadcast_sigma`,
`moment_scan_closure`). Dual-`A` recurred (`ld-ubld-cell-system` per-cell dense `A=G+F_out+ΣtM`
vs `ld-ubld-unified-moment-residual` operator `A=(L+C-S)`) — both traveled BYTE-EXACT + inline-
defined at point of use + the operator one is a legit L+C-S / F_in-on-RHS split (−B carried by
F_in, NOT forgotten) ⟹ CONCERN, exactly per Sharpening-2. So the whole in-diff review was PASS.
The decisive finding lived OUTSIDE the diff: `git grep methods/sn/index -- '*.py' '*.rst'` |
filter-to-moved-topics surfaced 4 un-repointed TEST files. Discriminate by CAUSATION with a
HEAD-vs-worktree count: `git show HEAD:index.rst | grep -c <label>` NONZERO + worktree `grep -c`
ZERO ⟹ THIS stage moved it ⟹ MUST-FIX (test_mms_2d.py:13, the 2-D Cartesian MMS `:doc:` still
naming index — the direct test twin of the mms/sn.py:646 repoint the author DID make). The other
3 (heterogeneous/nonvacuum, index=0 but moved in Tier-1) are pre-existing ⟹ SHOULD-FIX-for-
completeness (the author's own transform-5 already reached the mms/sn.py twins of that family, so
leaving the test twins is an internally-inconsistent half-cleanup). NOT-flagged control: the
curvilinear/aniso refs where the target is SPLIT (verification=34 / index=51) — content genuinely
still partly in index, so NOT cleanly stale; a HEAD-vs-worktree count keeps me from over-calling
those. Rule: a doc-split's transform-5 (external stale-qualifier repoints) is a RETIREMENT whose
blast radius is code + tests + docs (three searches, standing L-004) — the author almost always
greps the SOURCE tree and forgets `tests/`; the twin the MOVE ITSELF freshly staled is the hard
MUST-FIX, and the HEAD-nonzero/worktree-zero count is the causation oracle. Verdict: MUST-FIX
(1 in-scope + 3 completeness), everything else exemplary.

**Sharpening 9 — the CONCERN-DISCHARGE edit pass (harmonize a deferred notation): meaning-preservation + PROSE-twin hunt + the code-literal seam (dual-`A` one-pass, 2026-07-21, SHIP-WITH-FIXES).** After carrying the deferred dual-`A` CONCERN across all 4 splice stages (each traveled span kept a stale `A=L+C`/`A⁻¹`-for-sweep), the user ordered the one-pass harmonization — a SANCTIONED edit on traveled text (the deferral discharged). This is a DIFFERENT review shape from the byte-exact splices: not "did anything move?" but "is each edit meaning-preserving AND does it complete the law?" Method:
- **(a) The edits themselves: verify meaning-preservation against the surrounding math, then accept the honest-form judgment calls.** C6 de-bound `A=(L+C−S) ⟺ (L+C)ψ−Sψ=q` → `(L+C−S)ψ=q ⟺ (L+C)ψ−Sψ=q` — I certified it MEANING-PRESERVING *and* a type-improvement (equation⟺equation beats definition⟺equation) with −B correctly riding as F_in per the Stage-2 cartesian ruling. C8/U4 keep a LOCAL `A` (not respelled to `(L+C)⁻¹`) — ACCEPT: the swept object genuinely ranges over {L+C, M} (the reified G-S splitting), so `(L+C)⁻¹` would be WRONG; a local-`A`-bound-to-"the swept composite" + "inner kernel of the honest (L+C−S−B)⁻¹, never the full solve" is the strongest honest form (the §15A.2 precedent). U4's double-`A` ("A here = augmented loss L+C … the swept sub-composite of the honest A=L+C−S−B") is an INTENTIONAL local-vs-page-wide contrast, not a confusion — don't "fix" it.
- **(b) The miss is PROSE, and the author's machine scan is blind to it.** The pass's own scan "says the chapters are clean" — but it matched FORMULAS (`A = …`, `A⁻¹`), not PROSE. It renamed the Phase-C apply operator `A`→`L+C` at the formula site (U1) but left the SAME operator as bare "the operator :math:`A`" in the Key-Facts summary bullet AND the section-intro sentence of the SAME `phase-c-apply-sweep-equivalence` section — a pass-introduced internal inconsistency (A at the intro, L+C 26 lines down). RULE: on a notation-harmonization pass, grep the PROSE forms ("the operator A", "to the operator A", "discrete operator A") separately from the formula forms; the twins the formula-scan misses are the summary bullet + the section intro that paraphrase the renamed formula. Ground the rename first (here: `phase-c-cell-update`'s `(Tψ)=S+R+Σ_tψ` with `S`=SPATIAL-STREAMING face-flux — NOT scattering despite the letter — proved apply=L+C, so the operator IS L+C and the bare-A prose collides).
- **(c) The math-symbol/code-variable seam must not drift CODE-LITERALS.** The pass renamed the per-cell matrix `A`→`A_{\rm cell}` (math) but the CODE keeps dict-key/variable `A` (verified `return {"A": A, …}`). Correct for EQUATIONS and math-role prose. But it also rewrote a code-literal `` `A == A` `` pin-name → `` `A_cell == A_cell` `` — and the actual test spells the pin `` `A == A` `` (dict-key "A"). A double-backtick names a CODE construct; renaming it to a symbol the code doesn't have drifts the doc from the code (Cardinal Rule 3 / anti-#20). The brief's OWN seam rule ("code-literal ``A`` untouched") is the discriminator: math→A_cell, code-literal→A. When a rename pass has a math/code seam, grep every `` `A` `` (double-backtick) touched and confirm it still matches the live code symbol. Verdict: SHIP-WITH-FIXES (2 SHOULD-FIX — the prose twins + the code-literal drift); once fixed, the standing dual-`A` CONCERN from the 4 stage reviews CLOSES.

**Sharpening 10 — the FRESH-AUTHORED chapter (no splice): the mechanical char-identity recipe does NOT apply; the whole text is the honest-algebra + claim-truth surface, and the recurring catch is an OVERCLAIM the chapter contradicts ELSEWHERE IN ITSELF (curvilinear_multigroup ch. — the Part-B energy chapter, 2026-07-21, SHIP-WITH-FIXES, 1 MUST-FIX).** A ~400-line thin standalone chapter authored fresh (mints ZERO eq labels — pure composition-by-reference of `:eq:`multigroup`+`balance-general`+`alpha-recursion`+`mm-weights`+`wdd-closure`). Because there is no traveled-verbatim span, the L-013 Sharpening-2 "hold only fresh content to the standard" carve-out is VACUOUS — EVERY line is fresh, so the full honest-algebra + claim-truth standard applies to all of it (and it PASSED honest-algebra perfectly: `A=L+C−S−B` page-wide, sweep always `(L+C)⁻¹`, fission always external `q=Fψ/k` — the standing deferred dual-`A` CONCERN never even arose, 2nd fully-clean chapter after curvilinear_one_group S1). Method that replaces the char-diff:
- **(a) Claim-truth is the spine — spot-check the load-bearing claims against the LIVE code, not the diff's prose (leg-3).** The chapter's thesis was a three-tier decomposition (group-blind structure / group-diagonal data / group coupling in S,F); I verified the tier boundary EXACTLY: the unlabeled `denom[g,n]=2|μ|A_down+angular+Σt,g·V` display is byte-identical to `cell_balance_for_streaming`'s docstring math with `total_xs:(ng,)` the sole group input; `GeometryCoefficients` "No ng axis" + `CollisionCache (N,ng,nx)` confirmed in cache.py; `build_within_group_system` literally builds `A_AA=L+C−S−B_a` with the exact "within-group fission is zero (enters as q_ext)" comment; the for-g assembly loop confirmed. Verification pointers are the highest-value/most-falsifiable — I confirmed EVERY test name, tolerance, xfail, and guard: `test_kinf_homogeneous` (coord×ng×solver grid, xfail sphere-4g-krylov #200, rtol=1e-10, ρ/(1−ρ) amplification), the #196 gates (`test_si_krylov_eigenvalue_equivalence_sphere/_cylinder`, `_SI_KRYLOV_KEFF_TOL=1e-7`, floor ~1.9e-11, `assert maxmin>1.2` non-flatness), both multi-region classes + their 5 properties. A fresh chapter's verification section is a table of grep-checkable assertions — check them ALL.
- **(b) The recurring catch in a fresh composition-by-reference chapter is the OVERCLAIM, and the tell is that the chapter says the HONEST thing elsewhere in itself.** The sole MUST-FIX: line 130 "S, B, F **byte-for-byte the operators of** slab_multigroup" — but the SAME chapter says the correct weaker claim at line 257 ("the operators are consumed by composition, not re-derived") and line 118-119 (B "preserves energy just as in the slab" — energy-diagonality, NOT identity). "byte-for-byte" is a bit-identity claim: false for B (curvilinear reflection ≠ slab reflection at the angular level; code B_a carries the cosine-weighted |Ω·n|·w metric, geometry-specific), imprecise for S/F (same construction, not same array across meshes). Three-leg VIOLATION test passed: named future edit (a maintainer writes `array_equal(B_slab,B_sphere)` or caches B across geometries) + NOT-coextensive-today (diverges NOW for B, not latently) + verified-against-live-tree (the chapter's own two honest phrasings + the code). Fix = the brief's own "the same operators, consumed unchanged." RULE: on a fresh chapter, grep for strong identity/equality words ("byte-for-byte", "identical", "the same … verbatim") and check each against (i) the code and (ii) whether the chapter states a weaker/true version of the same claim elsewhere — the internal contradiction IS the proof it's an overclaim, not a judgment call.
- **(c) The §-number cross-ref gotcha: a bare "`:doc:`X` §7`" attribution is verified against SIBLING usage, because the numbering convention may count the title as §1.** "The inner-loop machinery of `:doc:`slab_one_group` §7`" LOOKS wrong (slab_one_group's 7th `=`-underline is "Verification hooks"; "Source iteration and its alternatives" is the 6th) — but slab_multigroup:624 cites the IDENTICAL "slab_one_group §7" for the same resolvent content, so the author consistently counts the TITLE as §1 (→ Source-iteration = §7). Consistent-with-sibling + the prose content-description self-disambiguates ⟹ CORRECT, not a finding. Don't over-call a §-number until you've checked how the sibling chapters cite the same section AND whether the surrounding prose names the content unambiguously. Everything else PASS: all 9 `:ref:` + 5 `:eq:` attributions resolve to correct live homes (the `choosing-inverse-realisation` parenthetical rightly targets slab_multigroup:844); the 3 sibling edits (toctree insert, router clause, two "What broadens next" `:doc:` links) correct/minimal; house style clean (Key Facts `:class: tip`, flat-top `=`/`=`/`-` matching siblings, no padding, WHY-monolithic rationale present); watch-items (b)/(c) from the brief were elegance WINS — the three-tier decomposition NEVER claims sweep coefficients are group-independent (it explicitly splits group-blind structure from per-group data), and the no-MMS admonition explicitly says the 1G MMS NULLS group coupling (covered by the 2G gates), not covers it. Verdict SHIP-WITH-FIXES on the single trivial overclaim swap.

**Sharpening 8 — DEMOTION stage (first of the campaign) + two-script double-application guard + the grep-null self-catch; the 4-stage Part-B campaign CLOSED (curvilinear floor→verification.rst ch. Stage-4, 2026-07-21, COMMIT-READY-after-1-MUST-FIX).** The finale moved a top-level H1 floor DOWN into verification.rst as an H2 — the REVERSE of every prior stage's promotion. Method points:
- **(a) Demotion direction reverses the map + introduces a DEEPER level.** `{=→-, -→~, ~→^}` census {1,8,6} (H1→H2, 8×H2→H3, 6×H3→H4). Char-identity certified perfectly (15 demotions exact, lengths preserved, + only F1/F2 = the two "above"→`(:doc:`curvilinear_numerics`)` fixes for the S3-departed ERR-058/#196 family). The NEW `^`/H4 level is the load-bearing check: extract the DESTINATION page's whole header tree, assert (i) the demoted subtree nests under the page's existing "Verification" H1 with ZERO global jumps, (ii) the `^` matches the page's H4 convention, (iii) the destination's next sibling (`ld-cartesian-2d` H2) stays intact at the right level directly after. All held.
- **(b) Two-script-split → double-application guard.** The splice aborted mid-run (a stale staying-landmark assert) and resumed via a second script; a resume CAN double-write. Single-homing (anchor once) rules out duplicate ANCHORS but NOT duplicate CONTENT (a paragraph can double with no anchor) — so the TRIPLE guard is: char-diff shows the span content exactly once (SequenceMatcher, no extra insert hunks) + single-homing (9/9 anchors verif≥1/idx=0) + `grep -c '^.. _anchor:'`==1 per anchor. All three confirmed no double-write.
- **(c) The grep-null SELF-catch (leg-3 turned on my own tooling).** My first `grep -cE 'Case singular-eigenfunction\|...'` returned 0 for verification.rst and nearly manufactured a false MUST-FIX ("cases/sn.py repoints a section that isn't there") — but the `\|` was a LITERAL under `-E`, a pattern bug. A corrected `grep -niE 'singular|eigenfunction'` found the section at verif:4173 (the brief's "~3493" was the PRE-splice line, +680 from the floor insert). RULE: a SURPRISING null/absence that would justify a MUST-FIX must be re-verified with a second, differently-spelled grep before flagging — verify the live tree INCLUDING my own tooling; a null can be a pattern bug, not a real absence.
- **(d) The foundation-doc straggler recurs EVEN with the lesson pre-applied.** Despite the author proactively applying my S2/S3 foundations-sweep lesson (reduced_operator.py floor→verification, cases/sn.py, 4 MMS tests all correctly repointed), `structured_geometry.rst:288` (`:ref:`sn-curvilinear-aniso-norm-reconciliation` in :doc:`index`` — the SAME floor anchor reduced_operator.py fixed) was MISSED → the lone MUST-FIX. The foundation docs (structured_geometry / boundary_conditions / coupled_block / operator_algebra) are the PERSISTENT blind spot of every doc-move; the per-anchor oracle over ALL of docs/theory/foundations/ is the catch. The bare-`:ref:`-vs-`:ref:…in :doc:index` discriminator held a 3rd time: index:128's bare `:ref:` to the moved floor resolves globally = NOT stale; structured_geometry:288's page-qualified one = stale = MUST-FIX.
- **Campaign arc (reusable playbook):** promotion (ch.1 Stage-1) → multi-span promotion w/ differential depth (ch.4 Stage-2) → chapter-MINT/gather (ch.5 Stage-3) → DEMOTION-to-verification (Stage-4). Across all: byte-exact splice is scriptable-verifiable (never eyeball); the live finding is ALWAYS the L-004 external blast radius (esp. foundation-doc inbound `:doc:index` twins the source-tree grep misses); the fresh skeleton/intro is the only honest-algebra surface (traveled dual-`A` = deferred CONCERN). Review loop demonstrably closes — all my S1/S2 findings were actioned by later stages.

**Sharpening 7 — chapter-MINT (gather) stage: the bycatch-twin, the umbrella supersession banner, and the review-loop-closed positive (curvilinear_numerics ch.5 Stage-3, 2026-07-21, COMMIT-READY, 0 MUST-FIX).** Final stage MINTED the 2nd Part-B chapter by gathering the Phase D/F/ERR-058/#196 chronicle out of the monolith (span A `~`→`=` {2}; span B `-`→`=`/`~`→`-`/`^`→`~` {3,40,9}). Char-identity certified PERFECTLY (both censuses exact, all lengths preserved, the ONE declared in-span fix FB1 = a retired-symbol `:func:`→literal isolated to one hunk, + the EOF-blank NIT). Header tree {6=,40-,9~} zero jumps; rump = Sweep-Algo H1 router → Unified-sweep-dispatch H2 with EXACTLY its 2 surviving H3s → floor H1, no dangling refs (grep the retained W6 head for `:ref:`/below to the CUT content — clean). Four method points:
- **(a) The bycatch-twin.** When a stage opportunistically fixes a PRIOR-stage straggler (the author split a test_krylov compound ref: `transport-cartesian`STAYS-index / `sn-curvilinear-homogeneous-kinf-recovery`MOVED-at-S2→curvilinear_one_group), that fix almost always has an IDENTICAL TWIN the author didn't trip over — here test_invertible_operator.py:687 carries the SAME compound docstring ref, unfixed. Grep the moved-anchor NAME tree-wide (not just this-stage anchors) when a bycatch fix appears; the twin is the SHOULD-FIX (S2-class completeness, not this-stage-caused). DISCRIMINATOR that avoids over-calling: a `@pytest.mark.verifies("anchor-name", ...)` registry marker is PAGE-AGNOSTIC (name→home resolved by the tagging registry, no page attribution) → NOT a straggler; only PROSE/docstring page-attributions (`index.rst (anchor)`) are. So 1 line flagged, 2 marker lines correctly not.
- **(b) The umbrella supersession banner — don't over-call a per-section gap.** A chronicle chapter's fresh intro claimed "each section carries its supersession banner"; my per-section heuristic flagged Phase F + ERR-058 as bannerless. FALSE ALARM: Phase D opens an `.. attention:: Superseded by #282` UMBRELLA banner that explicitly scopes "the ... claim in **the three sections below** is historical" (covers Phase F/ERR-058/#196), ERR-058 carries a `Status banner: #195 CLOSED` (it's the terminal CORRECT diagnosis, not superseded), and Phase F has an inline "*(since superseded.)*" marker. RULE: before flagging an intro's "every section is bannered" as an overclaim, grep upward for an umbrella banner that scopes "the N sections below" AND check whether an un-bannered section is a terminal-correct-conclusion (gets a Status/CLOSED banner, not a supersession one). The framing was honest.
- **(c) The "A/D/F above" cross-chapter honesty check.** A pre-applied prior-stage fix reworded a claim to "Phases A, D, and F above; B and C with the production machinery in :doc:`curvilinear_one_group`" — verify the RELATIVE-POSITION word ("above") lands true in the NEW file: A(ch.153)/D(191)/F(1153) ARE above the ERR-058 line (2052) that makes the claim; B/C genuinely in the sibling chapter. A relative-position claim that traveled must be re-checked against its NEW neighbourhood.
- **(d) The review loop closed — a positive worth logging.** All 5 of my Stage-2 MUST-FIX findings (bc.rst:655/738/743/3804 + coupled_block:65, the inbound foundation-doc `:doc:index`→moved-anchor stragglers) were RESOLVED this cycle (all now → curvilinear_one_group); the coordinator explicitly applied "your Stage-2 foundations-sweep lesson proactively." When a multi-stage campaign actions my prior findings, VERIFY they landed (grep the exact sites) and CREDIT it — the loop working is the signal the reviews are load-bearing. Fresh skeleton honest-algebra clean; span B byte-exact chronicle carries the standing deferred dual-`A` CONCERN (traveled-verbatim, harmonization-deferred), no new action.

**Sharpening 6 — multi-span stage: (a) differential promotion depth makes the header-tree nesting check load-bearing; (b) the inbound blast radius is FOUNDATION docs cross-referencing the moved anchors; (c) tombstoned-gate-chronicle travels with its campaign narrative, not to the living verification page (curvilinear_one_group ch.4 Stage-2, 2026-07-21, COMMIT-READY-after-5-MUST-FIX).** A three-span stage (W2a `-`→`=`; W2bc TWO-level `~`→`=`/`^`→`-`/`'`→`~`; W2d `-`→`=`/`~`→`-`) — DIFFERENT promotion maps per span, so classify PER SPAN (a union map is ambiguous: `~`→`=` in W2bc vs `~`→`-` in W2d). Char-identity certified perfectly: 3 spans, censuses {1}/{2,20,2}/{1,13} exact, all lengths preserved, F1–F5 each isolated to exactly one hunk (F1 W2a-end, F3/F4/F5 in W2bc, F2 W2d), only residual = the standing EOF trailing-blank trim (NIT). Three NEW method points:
- **(a) Differential-depth structural check.** When spans enter at different original H-depths but must land in the chapter's FLAT multi-`=` top level, the two-level promotion (W2bc H3→H1) is the load-bearing case — extract the WHOLE header tree and assert ZERO level jumps (here {13=,39-,14~} = Stage-1 {9,6,12} + Stage-2 {4,33,2}, no jumps; the 4 new machinery sections are flat-top siblings). Also verify the MONOLITH RUMP: the cut left the W6 H2 ("Unified sweep dispatch") with its 4 remaining H3s and the D/F/ERR-058 saga H2s intact, no orphan (going UP a level after a cut H3 is fine; a jump DOWN would be the smell).
- **(b) The inbound blast radius is the FOUNDATION-doc cross-refs, and it's the systematic miss.** The author repointed the obvious E1 (the moved symbol's own module docstring) + E2 (the direct gate test) — but MISSED five `:ref:`moved-anchor` in :doc:`/theory/methods/sn/index`` in FOUNDATION pages that cross-link the moved physics (boundary_conditions.rst ×4 naming `bare-sweep-extraction`/`sn-282-...`/`sn-pole-angular-closure-protocol`; coupled_block_operator.rst ×1 naming `sn-282-...`). These are the E1-CLASS silent wrong-target: the `:ref:` resolves globally to the chapter (primary link still works) but the `:doc:index` PAGE QUALIFIER is stale (index no longer hosts the anchor); no `-W` warning. Causation oracle nails this-stage (HEAD-idx-def=1 → WORK-idx-def=0, chapter-def=1). Verdict MUST-FIX (doc-correctness in heavily-xref'd foundation pages + internally-inconsistent-half-cleanup vs E1/E2 + fully-moved) though each is a one-line qualifier fix and the `:ref:` still works. The discriminator that keeps me from over-calling: bc:3271/3406 name Phase D/F anchors that STAY this stage (→ S3) and reduced_operator:74 names the aniso-floor that STAYS (→ S4) — all correctly NOT flagged, per the oracle. Rule: on a multi-anchor move, grep tree-wide `:doc:.../index` and run the per-anchor HEAD-vs-worktree oracle on EVERY named `:ref:`; the foundation docs (boundary_conditions, coupled_block_operator, operator_algebra) are where the inbound cross-refs to sweep/BC/route-a physics live and are the author's blind spot.
- **(c) The tombstoned-chronicle authoring adjudication.** A steer said the three-layer V&V ruling "governs the four nested MMS-gate subsections". The author kept them WITH Phase C (not verification.rst). CONCUR — and the concurrence is grep-decidable: the subsections carry explicit tombstone markers ("was **falsified** by W1–W4", "Gate 3.1 is therefore marked", "claims are **tombstoned** inline"), so they are RETIRED campaign chronicle, not LIVING case specs. Discriminator: living case → verification.rst (the W3b stress-MMS precedent); tombstoned/falsified/retracted gate → travels with the campaign narrative that explains its retirement (moving dead chronicle to a living-spec page pollutes it AND severs it from its explanation). Guard: concur only when a living-extraction commitment covers the still-living part (here the floor H1 → verification.rst at S4). The honest-algebra grep also found 2 traveled-verbatim dual-`A` spellings (`A^{-1}`≈sweep L2360; `A=L+C` augmented-loss L3545 — the latter even locally correct, the starting-direction march is scattering-free) → standing deferred CONCERN, not MUST-FIX (fresh Key Facts bullets clean).

**Sharpening 5 — the test-straggler MUST-vs-SHOULD discriminator is PER-REFERENCE / PER-LABEL, settled by
"does the moved content leave a correct anchor in the old page?" (curvilinear_one_group ch.4 Stage-1,
2026-07-21, COMMIT-READY-after-1-MUST-FIX).** The scripted byte-exact core certified PERFECTLY again:
31 non-equal hunks = exactly {26 promotions + F1–F5}; census by type {8→=,6→-,12→~} exact, all underline
LENGTHS preserved, F4/F5 SHAPES re-verified live on the migrated home (`ReducedStreamingOperator`
reduced_operator.py:421-428 — the fix kept the shape annotations, so keeping them is correct not a fresh
falsehood; VERIFY this, a transposed migration would inject a lie), F2 repoint's live consumers confirmed
(`cell_balance_terms`→`DiamondDifference.update` diamond.py:213; `_sweep_1d_*` retired = 0 prod matches),
header tree {9=,6-,12~} zero level-jumps, rump H1 router-only. FIRST fully-clean honest-algebra chapter —
ZERO dual-A anywhere (skeleton `A=L+C−S−B`, sweep never spelled `A^{-1}`), so the standing deferred CONCERN
didn't even arise. As Sharpening-3 predicted, the ONLY live findings were L-004 test stragglers the author's
source-tree sweep missed (E1 doc / E2 code-docstring / E3 test all repointed, but E3 was a CARTESIAN
straggler → the author never swept for CURVILINEAR test stragglers THIS stage's move creates; `grep -rn
curvilinear_one_group tests/` = 0). The NEW method beyond Sharpening-3: when a straggler is a COMPOUND
reference (one parenthetical citing ≥2 labels), run the HEAD-vs-worktree causation oracle PER LABEL, because
one can move while the other stays — and that partition IS the MUST-vs-SHOULD discriminator:
- **MUST-FIX** = the reference's target FULLY moved with NO correct anchor left in the old page — an
  ownership TRANSFER (a "canonical X gate" section whose title + `:label:` both left). Confirming query:
  `grep -c '<Section Title>' index` = 0 AND `grep -c ':label: <lab>' index` = 0 AND `=1` in the chapter.
  (test_streaming_equilibrium_curvilinear.py:48 — sole `:doc:index` pointer for "the canonical
  streaming-equilibrium gate"; the section+eq fully transferred → a maintainer following it lands on index,
  greps, finds the gate gone. Silent: `:doc:index` resolves, no `-W`.)
- **SHOULD-FIX** = a COMPOUND ref where one cited label moved and a sibling label STAYED and still anchors
  the prose (test_space_angle_separability.py:23 — ``index.rst`` :eq:`dd-curvilinear-scalar`(STAYS)
  / :eq:`mm-weights`(MOVED); both `:eq:` still resolve globally, only the page-attribution is half-stale).
  Don't over-call it MUST — the retained label makes the sentence partly-true and nothing is broken.
NOT-flagged controls (per-label oracle keeps me honest): test_phase_c_gates (#168 Phase C content is
Stage-2, retained now: phase-c-sweep-frame-matvec WORK=2), test_invertible_operator (transport-cartesian +
kinf-recovery, neither moved), test_curvilinear_pole_cell (dd-curvilinear-scalar STAYS), history.rst (bare
cross-chapter changelog refs, no page qualifier), index:275 "curvilinear geometries below extend" (points
at the RETAINED continuous-form Spherical/Cylindrical subsections, not the moved discrete balance — the
DISCRETE balance moved, the CONTINUOUS transport-equation forms stayed). Verdict: 1 MUST-FIX + 1 SHOULD-FIX,
both one-line test-docstring repoints, bundle into the same commit; the doc carve itself is exemplary.

**Sharpening 4 — two new verification wrinkles: DECLARED-IN-SPAN-FIXES as the expected residual,
and a WRAPPER DISSOLUTION (cartesian_multid ch.3 Stage-3, 2026-07-21, COMMIT-READY).** Final stage
moved W7 (angular windowing, IDENTITY promotion) + W8 (boundary G-S, `-`→`=`/`~`→`-`/`^`→`~`) and
DISSOLVED the `sn-iteration-primitives` H1 wrapper between them. Two method wrinkles beyond the
Stage-1/2 recipe: (a) **Fixes-as-residual**: when the brief declares N in-span fixes (here 7:
5×"drop stale same-page `above`/`:doc:.../index` qualifier", 2× repoint-to-slab_one_group), the
byte-exact `SequenceMatcher` diff must show EXACTLY those N hunks and NOTHING else — census (1/8/8,
1/4/4) + promotion-count (9) still exact, 5/5 labeled + 10/10 all math blocks verbatim (fixes were
100% prose, never touched an eq). Then EACH fix must repair a GENUINE falsehood, verified by
grepping the referenced content's ACTUAL home: `above`→cross-file was true (wavefront-flux-cochain
lives in foundations/wavefront_cochain.rst; the leakage-note + σ_r-fold live in slab_one_group.rst:
642/481); `:doc:.../index`→same-file-drop was true (si-gauss-seidel-reification + windowing-retyped
now in cartesian_multid via W8; inverse-application-driver stays foundations so F4 rightly kept
"below" only on windowing-retyped). A fix that DIDN'T map to a real home would be editorializing —
flag it. (b) **Wrapper dissolution** = a mini-retirement: verify (i) the info-bearing content is
preserved faithfully at its new home (the Wave-E-R1/#163 lift paragraph → slab_multigroup "Dev
history", tense-shifted "lifts"→"lifted" + sentence-joined, all 3 facts intact, placed before R2 to
complete R1→R2→R3), (ii) the rest was PURE POINTER (Para-2 "documented in slab_multigroup; this
section retains…" carried nothing unique), (iii) the anchor is gone from SOURCE with zero inbound —
GOTCHA: `grep -rl` returned 90 hits, ALL stale `docs/_build` HTML; grep source-only (`--include=*.rst
--include=*.py | grep -v _build`) → empty, and the `-W` build (catches broken `:ref:` as ref.ref)
confirms zero live inbound, (iv) no DUPLICATION at the destination (`grep -c "Round 1"` slab_multigroup
= 1, the new para only). The delicate seam is the dissolution point itself (W7→W8): confirm no
orphaned wrapper title / dangling "this section retains…" survives and the two now-adjacent sections
hand off sensibly (windowing-retyped closes "the reified M … is `si-gauss-seidel-reification`" → W8
IS the boundary-G-S schedule that DEFINES M — a better seam than the wrapper gave). L-004 sweep CLEAN
this time (no test stragglers — windowing/boundary-G-S content is internal, few external qualifiers;
contrast Stage-2's MMS). Fresh Key Facts bullet honest (`M=L+C−B_lower` correctly labeled the G-S
*regular splitting*, NOT conflated with `A=L+C−S−B`; "boundary-transient accelerator NOT scattering"
scopes S's absence). Dual-`A` (forward=L+C-or-M) purely traveled-verbatim + no fresh colliding
binding ⟹ the standing deferred CONCERN, no new action. Verdict COMMIT-READY, zero MUST-FIX — the
cleanest of the three cartesian_multid stages.
