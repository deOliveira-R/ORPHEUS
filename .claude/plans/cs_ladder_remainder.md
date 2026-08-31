# Every operator is born bound — the CS-ladder remainder (CS4c pulled forward + CS2 residue), the opening act of Campaign 2

> ## ▶ RESUME SURFACE — RULED 2026-08-30 (user), work NOT YET STARTED
>
> **This file is the entry point for the session that opens this work. It was
> written at a deliberate CLEAR-CONTEXT boundary — no compaction brief exists;
> everything a fresh session needs is here or pointed at from here.**
> Written at HEAD `067c4824` (main == origin/main, tree clean, Campaign 1
> closed out the same day).
>
> **The ruling (user, 2026-08-30, verbatim intent):** *"We will open the next
> session with the CS ladder to tighten the consumed operators before
> reshaping the consumers."* Concretely, per the readiness assessment the
> ruling accepted:
> - **CS4c pulls forward per its own charter clause** (the clause: *"may pull
>   earlier than post-CS2 at its chartering if frame-independence holds (C8
>   reserves only the frame MINT and the L/B + R6 rows)"* —
>   `space_and_kernel_binding_campaign.md:1838` §CS4c). The condition is now
>   `[M]`-checkable and looked satisfied at ruling time; the opener VERIFIES
>   it (see §Opener below).
> - **CS2 shrinks to its residue** (the S3 identity flip, the `*`-densifier
>   retirement, the L/B + R6 xfail rows) — landing before or interleaved with
>   CS4c's later steps; the sequencing detail is a design-round fork.
> - **CS1.5′ re-scopes toward the homogeneous coda** (per the F5 SN-first
>   ruling; its carrier-side design was already demoted to a contestant —
>   `space_and_kernel_binding_campaign.md:250`).
> - **This is Campaign 2's opening act** (clean-before-extend for the pencil:
>   `BoundOperator` + the binding recipes are exactly what
>   `GeneralizedEigenPencil` will consume). The pencil / resolvent /
>   partitioning phases follow AFTER the ladder; the **re-posed §5b build
>   fork and the O-3 split stay PENDING USER FORKS** at that later boundary —
>   this ladder work does not consume or decide them.

## 1. Goal, separately from means (plan-authoring §1)

**Goal.** Apply-time overloading retires from the operator layer: kernels are
representation-free DATA, a constructor binds Kernel × Space → the operator
(one domain, one codomain, ONE public `apply`), and the space supplies
everything the binding needs — so the consumers (the solver methods, reshaped
in later campaigns) receive operators that are already whole.

**Proposed means** (chartered 2026-08-20/21, NOT re-verified since — the
authoritative charter text is `space_and_kernel_binding_campaign.md`):
- **CS4c** (`:1838`): the per-instance feeding-normalization census first
  (vv-principles #29 method), the #205 ndarray keep-ruling re-litigated at
  the k-eigenvalue call sites, arm deletion per binding, C space-mandatory,
  S → `ScatteringKernel` shell with the iso pair + `LegendreMomentScattering`
  as ℓ=0/moment bindings of ONE datum (the O7 twin heal), the
  `BoundOperator(datum, space)` dataclass-ABC binding base with the
  mint-and-FORGET contract, the RESTRICTION verb joining retract/embed, and
  **the tightness gate landing WITH the binding**
  (`bind(K)† = bind(K†) ⟺ the frame is tight`, with a deliberately non-tight
  frame as the ERR-039-class negative control).
- **CS2 residue** (`:394`): SN's bulk space as first-class
  Energy ⊗ Spatial ⊗ Angular with STRUCTURAL identity (the "S3 flip" —
  `space.py` still says "until the S3 flip" at 4 sites, `[M]` 2026-08-30),
  the legacy `*`-product densifier retirement
  (`_dense_axes_weights`, `space.py:757` "retired with the densifier in
  CS2"), and the L/B mandatory-space flips + the R6 uniform-guard row.
- **The CS4c recorded-debts definition** (the fuller scope list):
  `.claude/plans/frame_square_recarve.md:321-329` — Riesz legs first-class
  (`riesz_lower`/`riesz_raise` space-minted; `InverseMetricOperator`
  retyped/absorbed); `DualSpace` (F1-B); `_AdjointOperator` retired into the
  leg composition; Λ's factor collection; ~~matrix-metric for the slab's
  non-diagonal Gram~~ ✅ **DISCHARGED by P7 2026-08-30**; `dual()` functor +
  declared-predicates-become-theorems gates; the mint-and-forget of S's frame
  accessor.

**Done when** (inherited from CS4c's charter + the ledger): the remaining
strict-xfail rows delete — **C and S at CS4c, L and B + R6 at CS2** — for the
campaign total of all 12 P1 strict-xfails deleted
(`tests/sn/architecture/test_monomorphic_leaves.py`, whose own header says
*"the strict-xfail set IS this file's todo list"*; `[M]` 2026-08-30 the F
rows are already deleted — CS4a K2b landed that flip); no
`singledispatchmethod` on S/F `apply`; the non-tight negative control REDs;
#359's M-M spellings re-checked (⚠ partially DISSOLVED — see §3 staleness).

## 2. Why the position is good — the `[M]` readiness facts (2026-08-30)

These are what made the pull-forward credible; the opener re-verifies each:

1. **The (kernel, space) collapse mechanism ALREADY LANDED — via CS5, not
   CS2.** CS2's design input 2 (`:400` block) chartered identity-content +
   accessor machinery *so that* "the scattering binding collapses to
   (kernel, space) — the separate `quadrature=` argument disappears". CS5
   landed a superseding-in-mechanism answer: `Axis.generator`
   (provenance-never-identity; `[M]` identity-inclusion RAISES) +
   `Axis.generator_as(Quadrature, consumer=…)` — production-proven since the
   P4-remainder (the closure mints + `streaming_terms` use it). ⚠ CS2's
   input-2 LETTER ("the axis never holds the grid"; "accessor rebuilds from
   stored edges") is superseded; its INTENT (no identity welding) is
   satisfied and gated. A §3 note now sits at the input in the binding plan.
2. **P7 pre-built the leg substrate and discharged the matrix-metric
   bullet.** The `HilbertMetric` family (`orpheus/numerics/metric.py`)
   single-sources the two faces the Riesz legs wrap; `_AdjointOperator`'s
   sandwich (`operator.py:1372/:1375`) is confirmed
   metric-representation-agnostic — its retirement into
   `A* = domain.riesz_raise ∘ A.dual() ∘ codomain.riesz_lower` (the target
   kept verbatim in `docs/theory/foundations/frame.rst`'s CS4c note) is a
   re-plumb of two lines plus composition; `DualSpace.of` threads metrics
   (P7 S2); `InverseMetricOperator` (`operator.py:3561`) sits ready to
   retype/absorb; the DENSE frames DRESS, so the tightness gate's
   "deliberately non-tight frame" negative control is cheap (model it on the
   D2/D4 loaded-not-blind idiom in `tests/numerics/test_frame.py`).
3. **The acceptance instrument is small and live**: `[M]` the ledger is down
   to C/S (CS4c) + L/B + R6 (CS2); per-ROW marks, so partial flips delete
   their own markers (the `_R1_ROWS`/`_R2_ROWS` pattern,
   `test_monomorphic_leaves.py:~253-290,:705,:1046`).
4. **The frame-independence condition** (the pull-forward's license): the
   frame mint needs the SPACE's angular measure — reachable via the CS5
   generator channel — and the SN angular bulk space is axis-built in
   production (`[M]` the harmonic mint reads `angular_space.axes` and raises
   otherwise; gated).

## 3. ⚠ STALENESS — every count in the CS charters predates ~10 landed phases

The charters date 2026-08-20/21. Landed since: CS4b, CS4a-R, CS5-NODAL, the
whole streaming un-weld arc (P4.x, P4.9a/b, P4b, P4-remainder), P6's
dissolution, P7. Numbers that MUST be re-censused before designing to them
(plan-authoring §2 — a count without a fresh predicate is a hypothesis):
- the **131/43-site C-mandatory migration** and the **909-site sugar-factory
  sweep** (both pre-un-weld);
- the **6-of-12 per-instance feeding census** (R-A);
- the **21 surviving consumer-side field-receiver reads** (the S4-deletion
  scoping at `:1823`'s preamble);
- **#359's "three M-M spellings"** — ⚠ partially DISSOLVED: P4.9a killed the
  Morel–Montry twin outright and the closure family was renamed
  (`angular_closure` / `AngularClosureBase`); re-derive what remains before
  citing the item;
- the L/B rows' surrounding claims (the boundary operator moved through the
  B-campaigns; R6's cited `boundary.py:343` line number is surely drifted).
The campaign pattern's record: **9 of 10 phase openers have corrected their
own charter section.** Assume this one will too.

✅ **RE-CENSUSED 2026-08-30 (the 10th correction, session 1; ground memos
`scratch/cs4c_opener_structural_ground.md` + `scratch/cs4c_opener_count_census.md`):**
C-migration HOLDS drifted (150/133/45 test-tree; production 3/3 space-bound;
only **7** sites mint truly space-anonymous objects — the 126 `from_mesh`
spellings resolve a space at construction); the 909-site sugar sweep is
**DISSOLVED** (S5 landed `2690a434`, 2026-08-24); the 21 field-receiver reads
are **DISSOLVED** (S4 landed; 0/0/0); the R-A roster is **INTACT** (0
added/removed; 2 HOMO sites re-spelled post-K2 — the census re-run must
expect its HOMO row to move); #359 is substantially dissolved (comment
posted); R6's substance holds at `boundary.py:714` (line rotted). The
frame-independence condition **HOLDS via PROVENANCE** (the `generator_as`
channel; the charter's "space's angular measure" wording corrected in the
ground memo §B.4).

## 4. What the OPENING SESSION does, in order

1. **Reground**: this file top-to-bottom; then the two charter sections
   (`space_and_kernel_binding_campaign.md:1838` CS4c — including its
   binding-base design inputs and the SN-FIRST F5 ruling — and `:394` CS2 —
   including the three 2026-08-20 adjudications and the §3 note at input 2);
   then the recorded-debts block (`frame_square_recarve.md:321-329`); then
   `kernel_and_medium_objectives.md` and `cs4b_fields_design.md` §Round 2
   (the binding-base inputs) as needed. Reconcile hashes/claims against git
   per plan-authoring §7 — never trust this file's tense over the tree.
2. **Opener ground re-measure** (dispatch the explorer; the mandated
   shelf-life pass): verify the frame-independence condition `[M]`; re-census
   every §3 number with stated predicates + positive controls; enumerate the
   current xfail ledger rows; map the current S/C/F/L/B leaf surfaces
   (`singledispatchmethod` sites, the #205/#276 keep-ruling call sites);
   locate what #359 still names.
3. **Design round (user)**: the step decomposition (CS4c is campaign-scale —
   "phases subdivide further as needed"); the CS2-residue interleaving; the
   formal order re-ruling (the pull-forward clause licenses proposing, the
   user rules); the binding-base shape against the round-2 inputs.
4. **test-architect BEFORE the carve** (the MUST proactive trigger — this is
   an operator-algebra carve crossing subsystem boundaries), then main-agent
   surgical implementation per `delegation.md` (NOT `method-implementer`).

> ✅ **STATUS 2026-08-30 (session 1):** steps 1–2 EXECUTED (reground + the
> two-explorer ground re-measure); step 3's design round SUBSTANTIVELY CLOSED
> the same day — **the design record is `.claude/plans/cs4c_binding_design.md`**
> (rulings: two-space `BoundOperator(domain, codomain)` base; frame handed in,
> faces minted, HUB interns; the three-tier construction discipline; kernels
> collapse to one datum; NO kernel dagger — operator-side adjoint by theorem;
> solution path deferred to its own consumers campaign; all six opener
> suggestions accepted incl. Riesz legs in-ladder and witness-1's re-rule).
> Step 4 (test-architect) + step 0 of the ladder await the user's go.

## 5. Standing constraints carried across the context boundary

- **Acceptance baseline: 9950 / 0 / 19 sk / 227 des / 70 xf, 13 trees rc=0**
  (measured 2026-08-30 at the P7 exit gate; the per-tree driver pattern is
  `scratch/_p7_gate_driver.sh` — canonical runner
  `.venv/bin/python -O -m pytest -p no:randomly -m "not slow" -q`, SERIAL,
  detached `Popen` + Monitor for the ~90-min full set).
- Predicted-then-measured full-gate deltas are the campaign norm (7
  consecutive reconciled-to-the-unit); ⚠ include CORPUS-COUPLED harness
  gates in the arithmetic (a registry-parametrized gate moves when docs add
  labels — the P7 +1; plan-authoring surprise log 2026-08-30).
- Mutation batteries: in-process monkeypatch plugins (`-p` + `PYTHONPATH`),
  per-arm verdict tables, positive control, `--continue-on-collection-errors`
  with `^ERROR` counted separately; NEVER disk mutation, NEVER
  `git checkout` on uncommitted paths. ⚠ The metric resolves PER CALL
  precisely so `object.__setattr__` batteries stay honest — do not re-cache
  it (`space.py:_resolved_metric`'s docstring carries the measured red).
- L37: never edit production (or gate-read tests) under a running gate.
  Never `git add -A` (stage by path; read porcelain). Branch + ff-merge only.
- ugrep: no `^`/`$` inside alternation groups; completeness claims re-run in
  Python with positive controls; enumerations never through `| head`.
- Pending at the LATER Campaign-2 boundary (untouched by this ladder): the
  re-posed §5b build fork (`scratch/p0_product_form_measurement.md` — the
  measured algebra is the SUM `L = D − E_sp·S − E_ang·A_ang`) · the O-3
  three-way split (deferred `ChartConnection` re-name + R→`transport/` +
  #407 + P5's content) · `_coll_cache`/`_pole_mirror_cache` untouched ·
  do NOT convert `CollisionCache`'s eager build to lazy (`_build_count`
  instruments the P4.9b COUNT gate).
- Open issues this work may touch: #158 (the curvilinear M/V VALUE — stays
  refused; E1 pins that #409 is not re-cited) · #423 (polynomial Basis —
  blocks nothing) · #359 (M-M spellings — re-derive first) · #407 (rides
  O-3, not this ladder).

---

**⏸ STATUS 2026-08-31:** CS4c steps 1+2 MERGED @ `9fc8bf04`; step 3
MERGED @ `600c5c80`; **step 4 MERGED @ `b25f9006`** (13-tree gate [M]
**10106/0/19sk/66xf**, reconciled exact against pre-registration; the F
rebind — FissionMaterialField + IsotropicFission + the frame-ℓ=0
FissionOperator, every consumer honest, N2N harmonized, G-F1/G-F2 +
battery B-4, the archivist's 19-page corpus pass). The LIVE resume
surface is `.claude/plans/cs4c_binding_design.md` — its **§18
COMPACTION POINT #4** is the LIVE one. §17 still holds the step-5 ▶
block (goal, standing constraints, the §6b members of any arm
deletion) and remains the resume surface FOR STEP 5 — but read §18
first: it discharges §17's census obligation, strikes §17's ruling 1
(the F ≡ N2N collapse is REFUTED by data), and records the user's
2026-08-31 ruling that **#426 and #428 are fixed BEFORE step 5
resumes**. §16.8 is the step-4 close-out; §15/§13 are superseded
history. Re-anchor at §18 → §17, never from a summary.

⚠ **The 13-tree baseline is UNMEASURED at HEAD `01ed1d79`.** Step 4's
**10106/0/19sk/66xf** is the last full gate; the 2026-08-31 session
added 10 rows and predicts **10116**, unrun. Run the gate before any
merge that claims a baseline.
