# QA Lessons — hot digest

Read every dispatch. **Behavioral rules only** — imperative, standalone.

- **War stories, evidence, `file:line`, measured tables** live in
  `lessons_archive.md` (`## L-0NN`, ascending, L-001..L-067). Open it only for
  the exact `L-0NN` a rule points at. Never read it whole.
- **Doctrine is NOT restated here.** `vv-principles` (preloaded) owns
  anti-patterns #1–#17, Modes 7–12, the bit-identity criteria, 1-group
  degeneracy, the `catches` directive; `numerical-bug-signatures` owns Sig 1–10
  + H1–H5; `qa/AGENT.md` owns make-it-RED (#11) and the field-role contract
  (#10). A `[skill]` tag means the principle is already there — what is kept is
  the ORPHEUS mechanic or the procedural trap. §I is the map.
- New lesson ⟹ append `L-0NN` to the archive, land a 2–5 line rule here.
  Sharpen in place; never let this file grow narrative.

---

## A. Making a gate RED — mutation mechanics

**A1. Disable the OVERRIDE, not the value, when two paths are value-equal.** A
specialised `apply` overriding an inherited leaf-sum can agree to ≤2 ULP, so only
`array_equal` discriminates; rename the override away and every tooth must red. → L-024

**A2. Revert PRODUCTION ONLY (keep the new tests) to prove a fix's negatives
could have failed.** `git stash push -- <production files only>`; the red message
must NAME the original bug, not just `AttributeError`. → L-027

**A3. Mutate in-process (throwaway pytest plugin / monkeypatch); revert by
RE-EDITING** — never `git stash`/`checkout` a path with uncommitted state.
Untracked files make `git diff` empty, so the revert proof is gate-green-again +
zero mutation markers. A `-p <module>` plugin needs `PYTHONPATH`. → L-039, L-043, L-052

**A4. Your OWN mutation needs a bite check.** [skill: Mode-8 METHOD WARNING]
Residue: a capability REFUSAL is TWO-part — adding `apply_transpose` does not
lift `is_adjointable`/`is_invertible` (predicates defaulting False). And a
**0-call counter is a FINDING** (path unreachable), not an inert mutation. → L-061, L-062

**A5. When `simplify` is pathologically slow on the MUTATED expression, don't
wait on pytest** — call the `derive_*` builder with concrete Rationals and read
the residual directly. Seconds, and decisive. → L-029

**A6. Cripple a GENERATOR, not a value, for the sharpest coverage verdict** —
replacing `O_h`'s 48 ops with its 8 diagonal sign-flips (= `D_2h`) left a
182-test suite green. → L-062

**A8. Mutate INSIDE the object's algebraic class, or the reds are catching the
LAW you broke.** [skill: #18] A CONSTANT written into a linear operator's output
makes it affine → 60 Krylov/SI gates red; the realistic LINEAR bug (same 94k
rows) red exactly 1 of 5076. Ask of every red: *by what mechanism does THIS gate
see THIS property?* Over-power lies "richly caught" — the flattering direction. → L-063

**A7. ONE mutation direction is almost never enough — enumerate the leaks and
red each.** (a) *Capability default-OFF*: the factory AUTO-SELECTS the wider
shape, OR appends a PHANTOM length-1 axis (control:
`not hasattr(space,"factors")`). (b) *`xfail`→live flip*: red against the
re-introduced bug AND the EMULATED PRE-change behaviour — only the latter rules
out a gate already green at HEAD. (c) *Polymorphic hook*: override returns the
base type, AND override DROPPED (base `replace()` keeps state, so only the
empty-state tooth reds — and only if the test ADVANCES first). → L-032, L-038, L-041

**A10. Run the mutation over the WHOLE module tree in BOTH arms — the
symmetric difference turns "does an external pin exist?" from an argument into
a LIST.** `[M]` old-τ vs HEAD over `tests/sn`: 7 red only at HEAD, **32 red
only under old-τ**, 9 red in both (another agent's scope). The 32 named the
analytic-closed-form pins I had just concluded did not exist — my own draft,
refuted by my own measurement. Bite check first (the target gates must FLIP,
with a non-zero call count). → L-069

**A11. Two DUALS of A9, both flattering.** (a) **A normalized fingerprint reused as a
CHANGE detector inherits its deliberate BLINDNESSES** — `[M]` nexus `body_shingles` is
bit-identical under `rtol 1e-12→1e-6`, `max_inner 1000→50`, a re-baselined expected
value and a fixture-arg swap (it normalizes `Constant→"C"` for clone robustness), i.e.
blind to every Mode-8-class-7 decay cause; a ledger on it reports every decayed marker
FRESH. Intersect the fingerprint's invariance group with the change class — Mode 12,
asked of an INSTRUMENT. (b) **A recall counter placed DOWNSTREAM of a filter cannot
count what the filter dropped** — `[M]` `nexus runtime-ingest` printed
`nodes: 0 / unresolved: 0`, exit 0, on a real report, because all 339 file keys were
dropped by a path filter before reaching the resolver (absolute-vs-relative); the same
artifact joined **2892** nodes once normalized. Demand a per-REASON drop breakdown. → L-070

**A9. An AUDIT instrument needs one control PER STATE its predicate accepts.**
[skill: Mode-8 METHOD WARNING] A production predicate reused as a detector
inherits its OTHER meanings — `_claims_convergence` is also False for an EMPTY
history (= GMRES exited in 0 iterations = CONVERGED), so 44 of 90 census rows
were invented, in the flattering direction, while the positive control passed
(it only exercised the genuine branch). → L-067

---

## B. Where a gate is structurally blind (ORPHEUS shapes)

**B1. Mutate the SHARED source, not the dead-for-this-path method.** SWEEP and
MATVEC share only precomputed coefficients: three apply-path mutations gave
call-count 0 and identical error ladders (GREEN-BLIND on dead code). Instrument a
call-count FIRST to find which twin the gate runs. → L-036 [skill: Mode 11]

**B2. A round-trip / self-consistency test cannot pin a CONVENTION** — both arms
carry the stale input. Recurs at every scale (1-D `s_axes`, d≥2 `|μ_axis|`). → L-018, L-023, L-031

**B3. "The matvec is tested" needs an instrumented call-count, not a
round-trip.** SI sweeps never touch `loss_action`; the matvec runs only under
`inner_solver="krylov"` (measured 1600 / 0 kernel calls on an MMS solve). → L-018, L-021, L-033

**B4. A `@pytest.mark.slow` catcher is DESELECTED under the canonical
`-m "not slow"`** — a sibling of Mode-8. Never trust a plan's "test_X must red";
simulate the regression under the ACTUAL invocation. → L-053

**B5. A fixture SYMMETRIC in the axis under test cannot see that axis — make it
asymmetric.** Instances: a SQUARE `nx==ny` mesh hides axis-ORDERING, and the
algebra-law suite is swap-invariant anyway (linearity/homomorphism survive a
CONSISTENT transpose of all operands), so the catcher is a broadcast oracle on
`nx≠ny` → L-040; a UNIFORM fixture makes a per-cell and a global-mean check
indistinguishable, so one fixture must VARY along the non-reduced axis → L-030;
two SAME-AXIS faces make `|Ω·n|` bit-identical, annihilating the packing gate's
only knob-reader (`slot.slice_view(metric_flat)`) — `[M]` 0/10 red,
`changed=False` on every call; a y-face in the layout moves it 0.963 → L-065.

**B6. An A-vs-B INVARIANCE gate's coverage is the set of production lines that
READ the knob — grep them; it is usually ONE.** [skill: #23] Two runs of the same
code ⟹ blind to every non-knob-dependent mutation, so the CATASTROPHIC positive
control is INVALID (an identity kernel leaves it correctly 10/10 green); the
control must be knob-dependent — neuter the knob and the ACTIVATION leg must red.
Name the rows that structurally cannot see it (`vacuum`/`lambertian` never reach
the deck kernel) so their green is not counted. → L-065

**B7. A transpose/adjoint RECIPROCITY gate pins the transpose RELATIONSHIP, not
correctness** — green for ANY genuine `(S,Sᵀ)` pair, so Mode-12 blind to a
SYMMETRIC drop in both halves. Mutate BOTH ways and require the one-sided
`A∘A⁻¹≡I` companion; never let it be deleted on "reciprocity covers it". → L-060

**B8. The branch you are crediting may be dead under the shipped config.** A
"fires under quadrature Q" claim is a 3-line probe (`count_nonzero(|mu_x|<1e-15)`
→ zero at every LS order); likewise an instrument whose reflective arm never
runs. → L-016, L-059

**B9. A "the matrix says the operator is healthy" argument must cite a
certificate for the EXACT gated BC** — a sibling-BC certificate plus "same
mechanism" is inference. BUILD the fully-coupled matrix: `ρ_prod > ρ_matrix` ⟹
splitting/wall lag (honest); `ρ_matrix ≈ 1` ⟹ real consistency failure. → L-059

**B10. The headline category gate is usually the WEAK one.** `runtime_checkable`
checks member PRESENCE only — `isinstance` flips True on a monkeypatched attr and
stays False under a realistic PARTIAL leak. The direct `not hasattr(...)`
negatives are the defense; credit them, not the headline. → L-039, L-042, L-047

**B11. The MMS refinement ladder is BLIND to the diffusion limit — probe
`σ_t·h ≫ 1` on a COARSE mesh** (refinement drives `σ_t·h` thin, where flat-source
is fine; the failure lives where users run). Probe vs DD with an ε-scaled
diffusive material. A reflective `c≈1` probe is a TRAP — both schemes read ~82 %
wrong from non-convergence. → L-017

**B12. A CONVERGENCE flag at an eigenvalue entry is the OUTER fact only —
inner truncation under a power iteration is invisible BY CONSTRUCTION.**
`solve_sn`/`solve_sn_adjoint` warn on `max_outer`/`keff_tol`; a within-group
solve hitting `max_inner` never reaches the warning. Before crediting a
convergence sweep as complete, ask WHICH LOOP the flag belongs to, and wrap
`_certify_within_group_exit` for the other one. Companion holes: a suppressed
warning, an `xfail`-absorbed one, and `-m "not slow"` deselection. → L-067, L-053

**B14. The STATIC call graph cannot answer "did this test exercise it" — measure the
EXERCISED set with `coverage`, never with `calls`.** `[M]` 0 of 21 claiming tests reach
`Quadrature.ordinate_permutation` statically; **7** do at runtime — 100 % false-DEAD,
because every production call site is annotation-mediated (`quadrature.x()`,
`self.mesh.quad.x()`), which nexus #16 does not resolve. `nexus callers` → 0 and
`dead-functions` FLAGS the method. Recipe: `dynamic_context = test_function` +
`[json] show_contexts = True`, then join `contexts` to
`sphinxcontrib.nexus.runtime.build_node_index` spans (~15 lines; `[M]` 23/23 contexts,
1353 pairs, 1.45× runtime overhead). ⚠ It measures CO-EXECUTION, not co-constraint — a
candidate list to mutation-verify, never a licence. The rung ladder is
CLAIMED 21 → EXERCISED 7 → ASSERTED ≤2 → MUTATION-VERIFIED 0, and **no edge quality
separates rungs 2 and 3.** → L-070

**B13. A published COMMAND is a separate claim from the API it wraps — gate the
STRING.** [skill: Mode-8 EIGHTH class] `-W error::ConvergenceWarning` (4 doc
sites incl. the runtime message) does NOT parse — an undotted `-W` category
resolves against `builtins`; pytest exits ERROR, 0 collected. The file's own
"it is escalatable" test passes because it installs the filter
programmatically. Gate: `_pytest.config.parse_warning_filter(s, escape=False)`.
Also grep `pytest\.xfail(` — the CALL form can never XPASS, so its deferral is
immortal and its reason string rots unfalsifiably. → L-067

---

## C. Reference contamination & structural independence

**C1. The circularity test is: does the bug live on BOTH sides?** A re-encoded
production formula vs an INDEPENDENTLY-ASSEMBLED primitive is legitimate;
re-encoded-vs-re-encoded is circular. Two ORPHEUS tells that a *corroboration*
is only PROCEDURAL: (a) the two sides ride one **antiderivative identity** — a
discrete recursion summing `f` "confirms" a claim about the exact face value
`F=∫f` only because `F'=f`, so it is true whatever the claim is; (b) the
corroborating GATE's own docstring cites, as ITS reference, the very claim being
corroborated. `[M]` Q5.6.4: `α = κ·w_gl·ξ(e_arc)` to 1e-15 was offered as
evidence for "the exact face coefficient is `ξ`", and `test_alpha_closed_form`
names "Hebert 3.399 — α IS the tangential cosine at a half-angle boundary" as
its ground. → L-029, L-068

**C2. For a WEIGHTED value-pin, L11-independence is not enough.** The hand-ref
must carry EVERY weight factor AND the fixture must make a factor-BLIND formula
give a different answer. Prove both: hand-compute the blind number, then mutate
production blind and confirm only that gate reds. → L-046

**C3. Pin a re-baselined `.npy` to a STRUCTURALLY-INDEPENDENT reference** —
never to "whatever the changed leaf emits" (circular), never to old-vs-new ULP.
When the composite is byte-identical, the anchor is `composite − collision`. → L-049, L-034

**C4. Two independent implementations IS independence;
producer-vs-its-own-projector is not.** Compare production's emitted slot against
a TEST-SIDE `leggauss` ref, then separately pin the two projectors' agreement. → L-044

**C5. Recompute the OLD contraction on a structurally-independent table to
verify an API-migration's bit-identity** — a brief's "0 ULP" is a CLAIM. → L-051, L-052

**C6. For an adjoint, the independent reference is a DENSE matrix built by
LOOPS, transposed directly, composed with metrics by hand.** Re-derive the
inner-product identity first, then prove `(A⁻¹)ᵀ=(Aᵀ)⁻¹`. An ASYMMETRIC transpose
pair is CORRECT when each transpose mirrors its OWN forward. → L-052, L-060

**C7. A two-paths oracle's analytical anchor is often TRANSITIVE and in another
file** — confirm that file is green before crediting analytical grounding. → L-028

**C8. Independence has TWO axes — DERIVATION and INPUT RESOLUTION; a
single-sourcing retirement closes the second silently.** [skill: #22] Tell: the
test builds ONE domain object and hands it to the SUT *and* the reference, where
the SUT used to resolve it from its own tag. Prove per-axis by mutating the
shared RESOLUTION (`[M]` an axis-letter x↔y swap left the "genuinely independent
routes" file 15/15 GREEN and reddened 78 siblings). Prove a helper's independence
mechanically, not by reading: `dis.Bytecode(f).codeobj.co_names` — the forbidden
names usually appear only in the DOCSTRING. → L-064

---

## D. Re-baseline & bit-identity integrity

**D12. A retirement's CONCEPT grep must cover the retired FIELD's hyphenation,
not only the retired SYMBOL's — and the gutted package's OWN sibling docstrings
are where the survivors live.** `[M]` a grep for "reflection-index table" caught
3 sites and missed 4 "precomputes the reflection-partner map at construction"
claims in the very `numerics/quadrature/` files the commit rewrote. No Sphinx
severity can see them (no `automodule` for that package) and the xref checker
tests TARGETS, not prose truth — grep is the ONLY gate, so its vocabulary is the
whole audit. → L-064

**D13. MOVING a method to a sibling object can convert a SELF-consistency into a
cross-object coupling guarded only by EXTENT.** Ask what array the old owner's
callers relied on; if the new owner carries a COPY cross-checked by
shape/length only, a same-size-different-values pair is now ACCEPTED where it
used to RAISE. `[M]` `to_local` moved operator→space: the gather reads
`op.indices`, the remap now reads `space.ordinate_indices`. And a round-trip
gate that was harmless while ONE array existed (`searchsorted(op.indices,
op.indices[perm]) == perm`) becomes the gap the moment there are two. → L-065

**D14. Before judging a re-baseline's LEGITIMACY, `git log` the snapshot's own
directory for a commit that ALREADY made the decision** — the reds may be its
REMAINDER, and then the question is completeness. `[M]` `39b46a31`
re-baselined 2 cylinder artifacts with a sha256 sweep + a per-artefact
`τ:=0.7` screen, and its universal *"all 23 snapshots … the only two that
changed"* was scoped to ONE directory while 7 more references had moved.
⭐ And check its per-artefact NULL reasons against EVERY mechanism the commit
bundled: one conjunction in the subject = two mechanisms = two checks owed.
`[M]` *"at M=2 the partition is BIT-IDENTICAL, so this case's τ did not
change at all"* — partition true (`5e-17`), τ moved **2.071e-01**, because the
same commit also retired the `[½,1]` absorber. Right conclusion, void
argument, durable certificate of blindness. [skill: #25] → L-069

**D1. Grep the WHOLE tree for the OLD literal, not the diff's touched files.**
A cross-check against a derived value WILL break (genuine miss); a
self-consistency round-trip survives while feeding wrong physics (latent stale). → L-023, L-025

**D2. Run the MASKING-CHECK on any loosened gate or regenerated baseline.**
Loosened → re-run the untouched arms, confirm they STILL hard-fail ≫ the bound.
Regen → OLD-snapshot-vs-NEW-code must hard-fail (load-bearing) AND
NEW-snapshot-vs-OLD-code must hard-fail (live gate). → L-022, L-028

**D3. Characterize drift from the BINARY:** `git show <c>~1:x.npy` vs
`git show <c>:x.npy`, then ULP-diff. Live-code vs regenerated-snapshot is
necessarily 0 ULP and characterizes nothing. → L-022

**D4. Large ULP at small magnitude is a metric artifact** — inspect the worst
element's magnitude before calling 256 ULP a violation. Conversely ~1e15 ULP on
ONE geometry with siblings green is a stale snapshot. → L-022, L-024, L-028 [Sig-10]

**D5. A HARD nULP floor and a STRICT bit-identity floor are different gates —
verify WHICH invocation ran.** Strict = `-W error::DriftWarning` layered on top;
`tests/sn/regression/conftest.py` downgrades it for its own dir but does NOT leak
to siblings (measured — assume neither way). Prove a strict floor live by
perturbing the baseline 1 ULP (`np.nextafter`). → L-014, L-015

**D6. Settle a byte-identity dispute with the IEEE micro-fact + `git status
--short '**/*.npy'` — NEVER a docstring.** `0.5*(a+b)==0.5*a+0.5*b` bit-for-bit
for all doubles (a `w=½` affine closure IS byte-identical to DD, contra its own
docstring); `2*X/D ≠ 2*X*(1/D)` at 1 ULP (the real re-baseline trigger); an
einsum spectator lift `fc->gc` ⇒ `fc...->gc...` is `array_equal` at rank-2. → L-020, L-028, L-032

**D7. When a carve preserves the COMPOSITE and not the leaf, prove byte-identity
on the composite DIRECTLY** (both emitted against a read-only baseline worktree).
A brief's "≤16 ULP" can understate leaf drift ~7× — say which object is pinned. → L-049

**D8. Prove a "verbatim relocation" by NORMALIZED AST-diff, not by re-running
gates** — substitute into the old body, strip docstrings/imports/blanks,
`difflib`. A true move reduces to the signature line plus the declared fork. → L-013

**D9. For an ADDITIVE-only change, grep the tree for ANY importer of the new
module (excluding its own tests). Empty ⟹ it cannot perturb a pre-existing
outcome** — stronger than re-running the baseline reds. → L-042

**D10. Prove a `singledispatch` alias rename via
`Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']`** — `Cls.apply is
Cls._apply_impl` is False (fresh descriptor per access), a red herring. → L-050

**D11. Do NOT trust "byte-identical EXCEPT one LATENT collision" — PROBE it.**
Compute BOTH branches every call and `array_equal` them across the FULL gate
suite (plugin reassigns the symbol in EVERY importing module; attribute by
`item.nodeid`; read under `-s`). Measured 48 divergences at 70 % ⟹ REACHED. The
two-paths gate that found them is Mode-11 blind (shared callee). "Latent via the
public entry" can be TRUE while "latent everywhere" is FALSE. → L-035

---

## E. Markers, levels, and the ORPHEUS audit surface

**E7. NEVER read a Nexus V&V number as a coverage claim — the surface is a SEARCH
relation wearing a PROOF relation's name.** `[M]` all 2748 `tests` edges are
`test→equation` (there is **no** test→code edge); all 16624 `implements` edges are
`source="inferred"`, **81 % on ONE shared token** (`operator`/`method`/`case`/`cell`);
`verified` is set iff `len(tests)>0` with **no** confidence floor, so **351 of 692**
"verified" equations have no declared test and a **CP** test "verifies" an SN
cell-flatten identity via the token `"cell"`. `nexus provenance` compounds it —
`implemented_by` is really *"documented on the same page"* (10 printed / 1 real).
Ask of any status: **what predicate sets it, and what is its weakest admissible
evidence?** → L-070

**E8. The graph's marker surface is PARTIAL — census it before designing a
marker-driven query.** `[M]` `foundation` (1515 usages / 308 files) and `regression`
(10 files) have **no node attribute at all**; only `verifies`/`vv_level`/`catches`/`slow`
are lifted. `catches` is an ATTRIBUTE, not an edge, and no `ERR-NNN` node exists, so
its claims cannot be joined to the catalog by traversal. No `.npy`/`.npz` snapshot is a
node either — a frozen reference cannot even be NAMED. → L-070

**E1. `@foundation` stacked with `@verifies("<physics-eq>")` is silent level
conflation** — the harness records both, so Nexus credits a physics equation with
a foundation test's parametrizations. Tell: a `documented` equation whose ONLY
coverage is a foundation test. → L-007

**E2. A `catches("ERR-NNN")` marker DECAYS** [skill: Mode-8 class 7] — re-verify
on every review of the file, not once at authoring. Same-area misattribution is
constant (an `A==A` matrix pin credited with an INFLOW-factor bug). → L-031, L-054, L-061

**E3. An audit-MISSING `catches` has FOUR outcomes — grep the production RAISE
SITE first.** (1) genuine catcher → tag, mutation-verified; (2) the catalog's L0
test was RETIRED and the marker did not migrate → re-tag the successor; (3) the
typed error is exported but NEVER raised → dead scaffolding, NO CATCHER, do not
invent a marker; (4) `assert_X` delegates to a WEAKER sibling → NO CATCHER. → L-054

**E4. `xfail(strict=True)` is satisfied by ANY failure — verify the REASON**
(`--runxfail`, match the documented `reason=`, then re-run without it). A stale
array index made a gate a FALSE xfail. → L-008

**E5. Two more level-marker tells, distinct from E1.** `foundation` under a
module `pytestmark=l1` emits `conflicting V&V level markers` and the intended
level is SILENTLY DROPPED. A self-generated regression baseline wearing `l1` is
conflation — fix the marker, not the file; its `_load_or_skip` should HARD-FAIL. → L-058, L-061

**E6. Audit mechanics.** Orphan triage order D→B→A→C (class D = existing test
needs only the label; ~25 % of orphans). `matrix.rst` LAGS a label rename — re-run
`python -O -m tests._harness.audit --gaps` for live spelling. `vv-status`
rationale comments use (parens), never [brackets] (docutils reads citations). → L-002, L-003, L-004

---

## F. Claim-scope — the claim is broader than the evidence

**F1. A "behavior-neutral" claim holds only for the ONE contract it was proven
against.** [skill: #12 / ERR-063] Residue: the proxies that fooled the closeout
were "no guard errors" and "DD snapshots didn't move" — DD is the SN-only
consumer where the field IS inert. **Run the slow accuracy-band suites.** → L-045

**F2. Exercised ≠ constrained.** [skill: Mode 10] Three states — nulled /
exercised-but-unconstrained / verified; never collapse the middle. With NO
isolating regime the STRUCTURAL pair (machine-precision threading + sign-flip ≫
tol + a no-op control) is COMPLETE; never manufacture a value-improvement leg.
Calibrate the tol live — a deterministic SI re-solve floors at 0.0. → L-026, L-037, L-038

**F3. NEVER accept a "designed-green / blind to this mutation class" narrative by
inspection — RUN the mutation.** [skill: Mode 12] Both directions are real:
over-claiming a catcher, and UNDER-claiming blindness (a leaf-transpose DROP is a
NON-transpose operator; its k SHIFTS). When the pushback lands on the skill's own
example, flag the skill edit as a finding. → L-058

**F4. A cited mutation MAGNITUDE for a metric-adjoint SOLVE must be the
full-solve value — RUN it, never the angular-collapsed 0-D proxy.** Metric
conjugation of a MUTATED operator is not spectrum-preserving. A never-asserted
cited number is still a plausible-substitution error. → L-058

**F5. For a "no missed site" dedup claim, the PLAN is the scope authority, not
the closeout.** A residual hit is a defect only if (a) a direct reconstruction,
(b) not transitively routed one level deeper, AND (c) in declared scope. → L-025

**F6. A stress-ansatz mandated by the test-architect memo is a binding
contract** — shipping the canonical `sin(πx/L)` 1G homogeneous case instead is a
gate DOWNGRADE. Flag it even when all tests pass. → L-019

**F7. "BC X is load-bearing because k = k_∞" is TRUE only for HOMOGENEOUS.** On a
heterogeneous reflective sphere the flux is non-flat, so the term DOES move k
(measured larger than vacuum). Check the config. → L-012

**F8. Check what the test HELPER tolerates before crediting an enforcement
claim** — a `squeeze_density` helper made the suite agnostic to `keepdims`, so
the bit-identity claim held only up to a squeeze. → L-042

**F9. "Matvec twin verified" is KERNEL-level; end-to-end Krylov≡SI is a separate
claim.** A loud `NotImplementedError` on the deferred half is the CORRECT interim
— but say so, and don't let a spec's wording credit the un-shipped half. → L-031, L-033

**F10. Every brief-named symbol/file is a CLAIM — confirm with `find`/grep before
editing** (two phantoms in one brief). Byte-compile no-test generator SCRIPTS
after a rewire: a broken import there is a breakage no test run surfaces. → L-051

**F13. A tombstone naming N carriers may name ONE carrier PER LEG — split the
retired claim into legs and mutation-check each.** `[M]` a deleted binding row's
successors were "the split gate AND the periodic gate"; dropping the `codomain`
binding reddened only the periodic gate, because the split gate's own codomain
assertion is `a.codomain is b.codomain` — `None is None`-satisfiable, unlike the
`x is <concrete object>` domain leg beside it. → L-064

**F12. A retired claim quantified over a COMPLEMENT has one clause per partition
class — enumerate before ruling obsolete.** The SN face partition is THREE-way,
so "residual zero at non-outflow" = inflow (INVERTED by #208) ⊔ tangential (still
exactly true). Rewriting the sentence around the inverted clause silently retires
the survivor. Check whether the file's own fixtures can even EXPRESS the survivor
(all 3 there carry 0 tangential ordinates; GL-at-even-order is the only
production quadrature with none). → L-063

**F14. FILL a plan's ⏳PENDING decisive row yourself — and use the plan's own
anchors as the probe's positive control.** A "decide nothing until X is
measured" row is the highest-value thing you can produce; the plan usually
states what X must reproduce, and that IS the control (`[M]` Q5.6.4: my live
probe hit both anchors `6.5934e-02` / `1.2676e-01` to the printed digits, which
is what licensed reading its NEW row `τ≡½ = 1.0181e-01`). Cache the expensive
reference to disk — the arm sweep is then cheap. Grep `Solution.__dataclass_fields__`
before assuming an attribute name (`keff`, not `k_eff`). → L-068

**F15. An "honest cost" is a COMPARISON — measure it against the candidates on
the table, never in the abstract.** A caveat "X has the usual exposure to Y" is
inverted if every candidate has Y and X has the LEAST. [skill: #24(c)] → L-068

**F11. ACCEPT a floor-CHARACTER gate; do not demand a floor-REMOVAL gate.** When
a fix cleans a RATE but leaves a floor, the honest claim is a falsifiable scaling
pin (`err(S32) < err(S16)/2`): a closure-BUG floor is quadrature-independent
(ratio ≈ 1 → fails). Here the pushback is against over-demanding. → L-009

**F16. An inherited blast-radius number counts a NAME, not a TYPE — re-measure it
with an in-process wrap before it sizes (or DEFERS) the work.** `[M]` "~87 reads"
was `grep '\.converged' tests/derivations/`, of which **72** belonged to an
unrelated result family sharing the attribute name; wrapping the actual class's
`__init__` + `__getattribute__` gave **33 constructions / 0 without the field / 2
reads** — 43× over, in the direction that defers a zero-churn fix. Route:
**Nexus for producers** (an attribute read is not an edge — `degree: 1` is the
graph being right), **dynamic wrap for readers**, grep only to enumerate
candidates. Pair the dynamic `0` with a static no-other-path proof (no `**`
splat / `asdict` / `replace` on that class) or it is "not observed", not "none".
→ L-066

**F17. A hardcoded status constant is a defect only if the producer ITERATES —
triage one hop UP before a grep-driven sweep.** `[M]` 7 hardcoded
`converged=True`; 3 sat on direct `scipy.linalg.eig` / `np.linalg.solve`
producers where `True` is honest. A "fix every hardcode" pass mints **false
honesty** at those. The lies and the facts are grep-identical — which is how the
lies hid. (Same shape for any `success`/`valid`/`exact` flag.) → L-066

---

## G. Doc / prose correctness (Cardinal Rule 3 findings, not V&V)

**G1. Campaign-narration staleness: the FIX bar is "provably lies about CURRENT
code", VERIFIED before ruling** — grep the named symbol/wiring tree-wide,
`gh issue view N`. Default = KEEP. Guards: a stale line inside a RUNTIME STRING
is behavioral → KEEP; a "failure here HALTs Phase X" banner is a record → KEEP. → L-055

**G2. Reviewing a skill→Sphinx distillation: verify code-anchored specifics
against CODE, never the skill twin** (the source's stale specifics propagate
verbatim). Python-domain roles are NOT `-W`-gated — a dead `:mod:` renders as
plain text. Grep the corpus: the OUTLIER spelling count is the bug. → L-056

**G3. Reviewing a results-compilation page:** a count DE-FREEZE is CERTIFIABLE
(live `--collect-only` proves the old literal lied); a doc RETITLE can beat the
test's own stale name/docstring — verify against the live `assert` body; a
run-book delegating detail to a config file may point at a contradicting note. → L-057

**G5. A CAMPAIGN-STEP NAME in a forward-looking claim is a self-expiring token —
when the step lands, grep the step's own name.** `grep -rn 'G6\.5'` minus the
retrospective forms (`since|at|\(|—`) found the 2 survivors out of 33 hits in one
command; both said "G6.5 retires the lengths" and G6.5 deliberately did not — one
of them 146 lines from the SAME FILE's corrected twin (vv #21 aggravator: the
file can now be cited for either). Also: a brief declaring "the known baseline
reds" declares the reds of the batteries IT ran — widen scope, reconcile the new
ones against the PARENT commit in a worktree before attributing. → L-065

**G4. A test's own prose is the least reliable thing in the file.** A "frozen /
bit-identical to the pre-carve path" docstring stales SILENTLY when its `.npy` is
regenerated and the test file is untouched — on any regen, grep consumers for
"frozen". A cited issue number can be wrong (trust git-archaeology). A prose "the
ERR-NNN class" citation is a nit; the same string in `catches()` is a defect. → L-020, L-028, L-030, L-034

---

## H. Mechanics, environment, probe hygiene

**H1. NEVER use a bare `assert` in your OWN `python -O` probe script** — it is
stripped; a throwaway printed "PASSED" while the values were visibly unequal.
Run teeth-checks through pytest, `np.testing.*`, or an explicit `raise`. → L-052

**H2. Settle a brief's Mode-8 hypothesis about a `tests/` subtree in 2 min, then
PIVOT.** Synthetic control + a falsified COPY of a real file, both modes — the
premise is usually REFUTED (0/676 inert); the real surface is `orpheus/` plus
NON-COLLECTED helpers. Pivot to an AST census of *what the asserts assert*: only
~29 % of bare asserts pinned a VALUE. → L-006, L-010

**H3. Replicate the test's OWN solve helper before judging a value claim** — a
naive `solve_sn_fixed_source(...)` defaults to vacuum; a divergent
hand-replication usually means YOU dropped the BC. → L-011

**H4. Baseline via a READ-ONLY worktree + `PYTHONPATH`, and verify it took**
(`git worktree add -d <ref> /tmp/x`; the editable `.venv` otherwise resolves to
the MAIN tree — confirm via `orpheus.__file__` / `inspect.getsource`). A worktree
pyright count needs the main `.venv` symlinked into its root. → L-041, L-045, L-049

**H5. pyright deltas are apples-to-apples only after line-stripping and per-file
reconciliation** — a `(file, rule, msg)` diff gives FALSE net-new when a
type-RENDERING string shifts. Rule out a masked offset with: full-tree total
exactly the baseline, SUT isolation 0/0, EMPTY diff on any reverted seam file. → L-027, L-039, L-050

**H6. Removing a `NoReturn`-poisoned return UNMASKS every latent error downstream
of the first poisoned call** (pyright suppresses after a `Never`-returning call).
Expect net-new ≠ per-file delta; classify each as latent vs regression. → L-050

**H7. In a SHARED working tree, diff ONLY your own touched files** — other
agents' edits make every dirty file look like yours. → L-054, L-055

**H8. Locating slow/timeout tests: batch into runs that COMPLETE** (a SIGTERM'd
run writes no junit-xml and loses the `-rfE` reasons). Mark slow PARAMS with
`pytest.param(..., marks=...)`, not the function; verify with `--collect-only`. → L-005

**H9. zsh does NOT word-split an unquoted `$VAR`** — a
`for M in ...; do pytest $SUITE -p $M; done` loop passed the whole list as ONE
argument and collected zero tests. Use a function with `"$@"` and print the
baseline `N passed` inside the loop. → L-062

**H12. The SUBJECT of your review can move while you review it — re-`wc -l` and
`git log -1` the reviewed document before writing the verdict.** `[M]` the Q5.6.4
memo grew 721→879 lines mid-dispatch (`8db88596` added a §9bis.9 whose literature
argument was the strongest defence of the link I was refuting). Also re-check the
BRANCH: the harness's session-start git snapshot said `main`; git said
`refactor/operator-strategy-layers`. → L-068

**H13. `grep "^FAILED"` on COLOURED pytest output matches NOTHING — a false
all-green, in the flattering direction.** ANSI escapes precede the `F`. `[M]`
my extraction reported zero failures beside a `41 failed` summary; only the
warnings lines leaked, because `-W error::…` contains `error`. And **never
pipe a BACKGROUND command through `grep`** — the task file then holds only the
filtered output and the evidence cannot be re-extracted (17 min lost). Always
`--color=no`, redirect FULL output to a file, filter afterwards. → L-069

**H10. Run `-rs` and READ skip reasons on any suite with a non-zero skip count**
— a skip reason containing an exception message is a permanently-dead gate. → L-061

**H11. `full_output=True` does NOT make a scipy status readable — `disp=False` is
the load-bearing half.** `[M]` with `disp` defaulted `True`, a non-converged
`brentq`/`root_scalar` **raises** instead of returning `converged=False`, so the
`False` leg is an unreachable branch wearing an honest name. On any "we read
scipy's flag" claim, check `disp` first. (Omitting `full_output` entirely is
*honest-by-raising* — fine, but not *readable*, so it cannot serve a
warn-don't-raise contract.) → L-066

---

## I. Already in the skills — point, don't restate

Homes are `vv-principles` SKILL.md unless named otherwise.

| Doctrine | Home | Archive |
|---|---|---|
| Test count ≠ coverage; het + multi-group + refinement | #3/#4; `bug-signatures` H1–H5; AGENT.md #5 | L-001 |
| Mode-8 `-O` strip: rewriter boundary, `testpaths`, the 6 fires-but-cannot-fail classes, the bite-check METHOD WARNING | Mode 8 | L-006, L-010, L-061 |
| `catches` = coverage CLAIM; mutation-verify the EXACT bug; markers decay | §Log every caught bug | L-007, L-031, L-054 |
| Mode-10 activated-but-unconstrained, incl. no-isolating-regime | Mode 10 | L-026, L-037, L-038 |
| Mode-11 gate-never-executes-the-rewired-path + plugin sentinel | Mode 11 | L-018, L-031, L-033, L-043 |
| Mode-12 invariant-functional; metric repair; commutator criterion; k-tooth 0.171 | Mode 12 | L-058, L-060 |
| Behavior-neutral zeroing/retype needs per-consumer VALUE proof (ERR-063) | #12 | L-045, L-048 |
| Sample-generates-the-group (ERR-072); partner ≠ bijection (ERR-073); monotonicity law; positive control | #13/#14/#15/#17 | L-061, L-062 |
| ⭐ A **refinement ladder** `8/16/32/64` is ONE congruence class, not "every order" | #13 (the ladder sharpening) | L-068 |
| ⭐⭐ Validating an **ADJUDICATING instrument** (ranks designs, nothing to mutate): the BASIS check (probe modes vs the problem's symmetry + what the rule can represent — the WEIGHTING can be provably robust while the basis is wrong), the RANK-CORRELATION check (≥3 candidates: which mechanism is the metric ordered by?), the cost-against-alternatives check | #24 | L-068 |
| Bit-identity vs principled-equivalence; the 3 criteria; AMPLIFY | §Bit-identity | L-022, L-049 |
| Stale-snapshot huge-ULP triage; splitting verified in a degenerate regime | `bug-signatures` Sig-10; Mode 9 | L-034, L-036, L-041, L-053 |
| ⭐ A BUNDLED change's per-artefact NULL reason must be checked against EVERY mechanism it retired (else a durable false blindness certificate); a re-baseline's radius is the frozen REFERENCES, not one directory's `.npz` | #25 (added by L-069) | L-069 |
| ⚠ Sig-10's sibling-pass discriminator is VOID when the changed code is single-geometry — SLB/SPH green carries NO information about a cylinder-only carve; bisect instead | `bug-signatures` Sig-10 | L-069 |
| Green gate = nothing until RED; SN `.apply`/`.solve` role contract | `qa/AGENT.md` #11/#10 + the role memo | §A |
