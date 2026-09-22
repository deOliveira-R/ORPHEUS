# Coding standards — evidence

Founding cases of [the coding-standards rule](../rules/coding-standards.md), moved verbatim (2026-09-20) from the `> [M]` blockquotes and dated parentheticals of `.claude/rules/coding-standards.md`. Each entry's first line names the clause it belongs to; the text below it is the original measurement, un-quoted and unchanged (its glyphs stay — it is history). The rule links here by heading anchor.

## Cases

### 2026-06-11 from axes roundtrip

Clause: Clean before extending

C5,
2026-06-11: `from_axes` round-tripped axes→legacy-mesh→axes, so 3-D admission would have
needed a new arm in the converter, constructor, AND trace gate — until the cleanup
inversion let 3-D flow through the one generic body.

### 2026-07-13 mass delete shadow

Clause: Retire as you go — mass-deletes are retirements too

2026-07-13: `15486f66` mass-deleted
`diag_cin_aware_split_basis_keff` while the CP rank-n protocol test's worker consumed it;
an untracked scratch copy masked the loss for ten weeks, then vanished — recovery had to
route through a surviving `.pyc`'s `co_filename` back into git history.

### 2026-07-15 sphinx severity

Clause: Retire as you go — the three searches, text-grep across code, tests AND docs (which Sphinx severities warn)

Measured 2026-07-15, Sphinx 9.1.0:
`:doc:` and `:ref:` **do** warn — `ref.doc` / `ref.ref` — so *page* moves and *label*
retirements ARE gated by `-W`; the silent class is the Python-domain roles, plus **raw path
strings** in prose/docstrings, which no build ever checks. A path assembled from segments —
`REPO_ROOT / "docs" / "theory"` — is invisible to a path-grep too; grep the **last
segment**.

This clause read
"unless the build runs `-n`" until 2026-08-03; that was false, and it told every
retirement audit it was covered when it was not.

### 2026-08-03 recaptured snapshots

Clause: Retire as you go — a rewire can silently DEMOTE a gate's claim class (check that the named replacement is real)

measured 2026-08-03: three `cyl_*` snapshots cited as the pre-carve anchor had all been
re-captured by the consolidation commit itself

### 2026-08-06 message fragment

Clause: Retire as you go — retiring a MESSAGE STRING, grep the shortest distinctive fragment

2026-08-06, G6.3 step
8.0: retiring `OperatorSum`'s inline check onto a shared helper reworded
`"OperatorSum requires equal domains"`. The audit grepped `requires equal domains` — which
matches the production line — and reported only the definition site. Two gates matched on
`"equal domains"` alone and went red in the wide run; a third reference was prose. The
correct pattern was the two-word fragment.

### 2026-08-06 matmul operand

Clause: Retire as you go — routing a call site through the ALGEBRA raises its operand requirement

2026-08-06,
G6.3 step 8: a `_NoTransposeLaw` stub with only `apply` hit
`TypeError: unsupported operand type(s) for @`.

### 2026-08-09 capability rows

Clause: Retire as you go — single-sourcing a duplicate demotes every gate that compared its copies

`[M]` 2026-08-09, #345. `capability_rows()` and the reference registry were two
hand-written enumerations, and their `r_0` name tags had *already* diverged
(`round(r0*100)` vs `round(r0/R_out*100)`, agreeing only because every shipped
`R_out` is `1.0`). Writing the promised row→registry join would have detected that;
hoisting the grid to one constant + one `reference_name()` made it **unspellable** —
and made the join tautological in the same commit. Step 2 saved it:
`test_builder_keyset_is_the_shipped_class_a_inventory` already pinned all 13 names
against a literal written independently of both, so the name set stayed anchored. The
join was kept, re-described in its docstring as testing the *discovery-and-registration
path* (which has no other catcher), and explicitly disclaimed as no longer able to
catch a spelling divergence.

### 2026-08-11 tau raw spellings

Clause: Retire as you go — a MATH symbol has three spellings, and a number is a fourth (and the concept grep's dual hazard)

2026-08-11, Q5.6.4: a retirement sweep briefed with `tau_raw` + `τ_raw`
reported clean. `docs/theory/methods/sn/angular_quadrature.rst:369` still
asserted *"`\tau_{\rm raw} \in [1/5, 4/5]` with the **bit-exact** reversal
identity"* — both halves present-tense-false after the ω-partition carve.
It was found only by grepping `tfrac15|tfrac45`, i.e. the NUMBER in its
LaTeX spelling. The page was also absent from the audit's own file list,
because that list was built by the same two-spelling grep.

⚠ Dual hazard, same measurement: the audit's file list was simultaneously
**over**-counted ~3× because `absorber` is also a *material* (`pure
absorber`, `cavity-absorber`) and `clamp` is also a GMRES `restart` clamp —
11 of 17 flagged pages were false positives. A concept grep needs its hits
triaged by MEANING before any of them is called a site.

### 2026-08-12 alpha dome assert

Clause: A bare `assert` in `orpheus/` is not a contract — the canonical runner strips it

`[M]` 2026-08-12. `α_{M+1/2} = 0` is a genuine admission contract on every curvilinear
quadrature (it is a *consequence* of the measure's antisymmetry, not an axiom of the
one-sided Lathrop–Carlson recursion). It was enforced on the **sphere** by
`assert abs(alpha[N]) < 1e-12` and on the **cylinder** by nothing at all. Demonstrated on
the verbatim recursion: a measure closing at `alpha[N] = +0.2000` is REFUSED under plain
`python` and **ACCEPTED** under `python -O`. Fixed at `bea6a367` — but note the fix was
*not* "add a check to the cylinder arm": the recursion had **three** copies, which is
precisely why the contract could live on one arm only. Cardinal Rule 2 first, then the
guard.

### 2026-08-12 blind pin twice

Clause: Retire as you go — the pin you reach for FIRST is the one most likely to be blind

`[M]` 2026-08-12, task #51 — **twice in one session, two unrelated blindness
mechanisms, both in the failing rows' own file.** (a) Cartesian octant: the
obvious licence was `test_2d_octant_sweep_closed_form_anchor` (`φ = Q/Σ_t`), green
at HEAD — but it is an all-reflective FLAT infinite medium, so it reads the
quadrature *only* through the total weight. `sum(w) = 4π` to **0.000e+00**
(bit-exact) at LS4/LS6/LS8 before AND after, while `μ₁` moved
`0.408248290463863 → 0.350021174581541` (14 %). It is a Σw-normalisation gate
(ERR-004/025), structurally blind to node PLACEMENT. The real licence was the LS
moment-exactness / advertised-degree suite. (b) Cylinder τ: every flat-flux L0
anchor is blind because the M-M recurrence gives `(ψ−(1−τ)ψ)/τ = ψ` for **every**
τ — including the `@verifies("streaming-equilibrium")` gate sitting in the same
file as three of the failing rows. The real licence was
`test_cyl_tau_equals_the_ANALYTIC_closed_form_not_the_chord_convention`, 1 of the
32 gates a whole-suite old-τ mutation actually reddens.

### 2026-08-13 pole mm recurrence

Clause: Retire as you go — a LABELLED EQUATION is an API

`[M]` 2026-08-13, task #67. `pole-mm-recurrence`'s first line read
`\phi_{1/2,i,g} = 0` while production marches the seed as an ODE. The page had ALREADY
condemned it twice — "replacing the hardcoded zero that Phase B had baked in" and a
`ZeroSeed` row reading "the pre-ERR-026 term-initialisation bug … wrong off flat flux" —
and a sibling page carried a subsection titled "The bug Phase B baked in". All ~2500
lines from the equation, none of it reaching it. Four sites inherited the zero,
including the page's own **Key Facts** card. Sphinx built clean throughout, and the
sole `verifies()` marker on the label was a suite that passes NO seed, so the matrix
reported the equation covered while the covering rows asserted the kernel default.

### 2026-08-14 q6c cross reference

Clause: Retire as you go — a CROSS-REFERENCE is a load-bearing dependency

`[M]` 2026-08-14, quadrature Q6-C. The theory page's stage-2 argument was rewritten to
say a degree means nothing without its reference measure, and cross-referenced the
*1-D primitive constructors* section for the Gauss rules. That section still stated
**both** rules as "`degree_of_exactness = 2n - 1`", with the Chebyshev one qualified
only by the prose "in the weighted sense" — i.e. exactly the bare-integer half-claim
the new argument exists to refute, sitting at the end of its own citation. Found by
the agent doing the repair, not by any build: `-W` was clean throughout, and the
section contains neither "reference" nor "claim" for a grep to catch.

### 2026-08-17 loss action promotion

Clause: Retire as you go — the MIRROR, a retirement can silently PROMOTE a gate's claim class

`[M]` 2026-08-17. `tests/sn/operators/test_loss_action_convention.py` asserted
`apply(ψ) == loss_action(σ_t, ψ) − C.apply(ψ)` and its own header called that check
*"tautological (`apply` is DEFINED as `loss_action − σ_t·ψ`)"* — true when written.
#257 S8b made `apply` σ-free (`loss_action(0, ψ)`), so the same line now reads
`loss_action(0,ψ) == loss_action(σ_t,ψ) − σ_t⊙ψ`: the **affinity of the walk in σ**,
a falsifiable property of two independently-evaluated walks. The gate went from
restating a definition to being the only check of an algebraic property, with **no
line of the test body changing**, and its docstring went on disclaiming it for months.

### 2026-08-26 curvature getattr

Clause: Retire as you go — a symbol grep cannot see a name that lives inside a STRING

`[M]` 2026-08-26, P1 item 8. Retiring `SNMesh.curvature` (whose `None`
**was** the Cartesian case), my residual grep returned only prose and I
called the set closed. `tests/sn/operators/test_native_matvec.py:392` read
it as `curv = getattr(sn_mesh, "curvature", None)` and branched on
`curv is None`, so after the retirement **every curvilinear mesh took the
slab branch** — 2 reds, sphere and cylinder. ⚠ The aggravator: I had run
exactly this string-form check for `mu_start` **one item earlier** and
confirmed it clean. The habit did not transfer across two commits by the
same author on the same afternoon, which is why this is a rule and not a
reminder. It failed loudly only because the assertion on the wrong branch
happened to be falsifiable; a `getattr` default that matches the common
case fails silently and green.

### 2026-08-27 p41a guard witness

Clause: Retire as you go — the ENFORCEMENT side, what KEPT THE COPIES EQUAL becomes load-bearing

`[M]` 2026-08-27, un-weld P4.1a. Retiring `ReducedStreamingOperator.coord` (a copy of
`mesh.coord`): what made it redundant was that each of the three factories *validates*
`mesh.coord` against the literal it then stored, so the identity held **by
construction**, not merely on 3/3 shipped fixtures. `grep "requires .* mesh"` returned
**3 hits, all three the production `raise` lines — zero witnesses tree-wide.** After
the retirement those guards are the only reason `op.mesh.coord` is the operator's
chart. The three `TestProperties` chart tests had been asserting the stored literal;
rewritten as the guards' witnesses (`vv-principles` #11 — one positive leg, two
negative legs each, matching the production message) they cost one edit and closed the
exposure in the commit that opened it.

### 2026-08-28 p44 dead references

Clause: Retire as you go — a COMPLETE import audit is still a PARTIAL audit; `dead_references` reads the docstring surface

`[M]` 2026-08-28, un-weld P4.4 (4 symbols, `geometry/` → `sn/mesh/`). The
import audit was done by **AST** (not grep), the residual check ran in
Python **with a positive control**, and returned **0**. The affected suite
was green, `pyright` 0, `sphinx -W` **clean**. `[M]`
`mcp__nexus__dead_references` then found **5 dead targets / 9 sites** —
`:class:`/`:func:`/`:attr:` cross-references to the old path sitting in
**docstrings** in `orpheus/sn/angular/closure.py`,
`transport/spatial/{scheme,cell_balance}.py` and two test modules. Nothing
else could see them: the module is not `automodule`'d, so `-W` is silent at
every severity. ⭐ The transferable half: the residual filter was
**validated and correct** — it was run over the wrong *surface* (import
statements), not with the wrong *pattern*. A positive control proves your
regex finds what you point it at; it says nothing about whether you pointed
it at the whole corpus. ⟹ **`dead_references` is the only instrument that
reads the docstring surface** — run it before calling any retirement or
re-home done, and again after the fix (this one went 9 → **0 dead / 52
checked**).

### 2026-09-07 guard is debt

Clause: A guard is elegance debt — tag it, and name what retires it

`[R]` user, 2026-09-07, ruling the R6 carrier guard: *"we're creating a protection, but
the ultimate state is to 'not need a guard'"*.

### 2026-09-18 underline scan

Clause: Retire as you go — a rename that changes a name's LENGTH owes the `.rst` corpus a section-underline scan

`[M]` 2026-09-18,
#412: two underlines short after pass 1, caught by the archivist's scan before the
build.

### 2026-06-22 naming rulings

Three rulings from June and July 2026, recorded in the main agent's memory until the `coding-standards` "Naming" section landed on 2026-09-22 (the memory notes `feedback_high_signal_names` and `feedback_naming_consistency_greppable` retired with it).

- **2026-06-22, `HarmonicMomentFlux` (the greppability law).** `HarmonicMomentField` was the one off-pattern member of a Domain-first `<Domain><Role>` family (`AngularFlux`, `ScalarFlux`, `BoundaryFlux`), so `grep Flux` missed it. Renamed `HarmonicMomentFlux`; `HarmonicMomentSourceSink` and `HarmonicMomentProjection` completed the matched set. A role-first proposal was refused because it introduced a THIRD word-order into the codebase. The user's ruling, verbatim: "if we're perfectly and predictably consistent, you will know exactly what to search for with a grep, and so will anyone else."
- **2026-06-23, the frame faces (native vocabulary; collisions; inputs name outputs).** Once the object was a discrete frame, its faces became `analysis` (`M`) and `reconstruction` (`R`, the dual-frame synthesis): not "projection", which already named the idempotent `R∘M` and so was the colliding member, and not `synthesis`, which in frame theory is the naked synthesis operator over the same family (`T*` in the Casazza–Lynch letters the corpus uses, `T` in Christensen's), a different object. `Frame(basis, measure)` has the basis space as codomain and the measure space as domain, so the inputs named the outputs. The FEM synonym went to the docstring. The user delegated the choice ("as long as it carries a strong signal"), and the name that emerged as understanding deepened was kept over the habitual one.
- **2026-07-10, `Composite[Interior, Boundary]` (structure over role).** The generic composite space was named by its structure, not by a domain role (`System[…]`): the domain reading is the consumer's and is carried by the specialization (`[R]` the user).

### 2026-09-22 memory file census

The blast-radius audit of 39 agent-memory topic files that the 2026-09-21 distillation had judged archaeology by name (`.claude/plans/archive/blast_radius_audit_2026-09-22/`). Three measurements changed the method mid-audit.

- **Agent-memory files are not graph nodes.** The distillation standard said they were, and every owner's brief repeated it. `[M]` `mcp__nexus__file_brief` on three memory files answered "not in the graph" (the orchestrator, the explorer and the cross-domain-attacker, independently); `.nexus/config.toml` sets `extra_source_dirs = ["tests"]` beside `orpheus/` and the docs. So `dead_references` is blind to a memory file's referrers, and the text census is the only instrument.
- **The stem census missed the wikilink spelling.** Memory files link each other as `[[hyphenated-slug]]`, the front-matter `name:`, while the census keyed on the underscore stem. `[M]` re-run in both spellings with a positive control per spelling (`[[issue-247-legA-review]]` found in its sibling): two referrers the stem census had missed, one for qa and one for the explorer.
- **Two censuses agreed because they shared an exclusion.** The orchestrator's census and the archivist's independent one (2 242 and 2 430 files) reproduced each other exactly, and both excluded every `scratch/` directory, which also excludes the TRACKED `.claude/scratch/`. `[M]` `git grep` after the deletions found `.claude/scratch/open_fronts_audit.md` citing two candidates, one of them an open tracker box. Agreement between two instruments with a shared blind spot is the X4 tautology in census form.
