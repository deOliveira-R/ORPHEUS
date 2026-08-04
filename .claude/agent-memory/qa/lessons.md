# QA Lessons — hot digest

Read every dispatch. **Behavioral rules only** — imperative, standalone.

- **War stories, evidence, `file:line`, measured tables** live in
  `lessons_archive.md` (`## L-0NN`, ascending, L-001..L-062). Open it only for
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

**B5. An algebra-law suite is swap-INVARIANT — blind to axis-ordering (Mode-2).**
Linearity and homomorphism survive a CONSISTENT transpose of all operands. The
catcher is a broadcast oracle on an `nx≠ny` mesh (square meshes agree in shape).
Never credit the law suite with swap coverage. → L-040

**B6. A per-AXIS invariant tested only with a uniform fixture is untested for
its per-axis-ness** — a per-cell check and a global-mean check are
indistinguishable there. One fixture must VARY along the non-reduced axis. → L-030

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

---

## C. Reference contamination & structural independence

**C1. The circularity test is: does the bug live on BOTH sides?** A re-encoded
production formula vs an INDEPENDENTLY-ASSEMBLED primitive is legitimate;
re-encoded-vs-re-encoded is circular. → L-029

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

---

## D. Re-baseline & bit-identity integrity

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

**F12. A retired claim quantified over a COMPLEMENT has one clause per partition
class — enumerate before ruling obsolete.** The SN face partition is THREE-way,
so "residual zero at non-outflow" = inflow (INVERTED by #208) ⊔ tangential (still
exactly true). Rewriting the sentence around the inverted clause silently retires
the survivor. Check whether the file's own fixtures can even EXPRESS the survivor
(all 3 there carry 0 tangential ordinates; GL-at-even-order is the only
production quadrature with none). → L-063

**F11. ACCEPT a floor-CHARACTER gate; do not demand a floor-REMOVAL gate.** When
a fix cleans a RATE but leaves a floor, the honest claim is a falsifiable scaling
pin (`err(S32) < err(S16)/2`): a closure-BUG floor is quadrature-independent
(ratio ≈ 1 → fails). Here the pushback is against over-demanding. → L-009

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

**H10. Run `-rs` and READ skip reasons on any suite with a non-zero skip count**
— a skip reason containing an exception message is a permanently-dead gate. → L-061

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
| Bit-identity vs principled-equivalence; the 3 criteria; AMPLIFY | §Bit-identity | L-022, L-049 |
| Stale-snapshot huge-ULP triage; splitting verified in a degenerate regime | `bug-signatures` Sig-10; Mode 9 | L-034, L-036, L-041, L-053 |
| Green gate = nothing until RED; SN `.apply`/`.solve` role contract | `qa/AGENT.md` #11/#10 + the role memo | §A |
