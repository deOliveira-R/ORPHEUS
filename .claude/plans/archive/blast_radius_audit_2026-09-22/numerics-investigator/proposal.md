# Blast-radius audit — numerics-investigator's three "archaeology" candidates

**Dispatch:** the harness campaign's memory-distillation close-out, owed item (1).
**Posture:** READ-ONLY on tracked files; this file is the only deliverable.
**Date:** 2026-09-21. Every status line below was verified against git / `gh`, never
against a memory's claim; the commands are quoted beside each.

## Verdicts

| file | verdict | reason (one line) | referrer edits |
|---|---|---|---|
| `phase5a_moment_consuming_scatter_derisk.md` | **RETIRE** | Both halves are in the corpus in richer form: the `Y_0^0 = 1` load-bearing fact as a named admonition on the angular-windowing theory page (and at its definition site), the methodology half as `vv-principles` anti-pattern #7 + the bit-identity criterion 2 + `instrument-doctrine` X4. | 1 line, in the distillation's own candidate list — see "Referrer edits" below. |
| `r1_step_d_sphere_preconditioner_oscillation.md` | **RETIRE** | ERR-050 catalogues the same defect with the same evidence table and MORE (how it hid, the catching test, the lesson) — and the file's stated root cause was structurally superseded three days after it was written, so the file is now stale-wrong about production. | 1 line, same. |
| `sn_sig_t_layout_drift_indexerror.md` | **RETIRE** | ERR-055 carries all three of its methodology lessons verbatim in substance; the fourth ("a deterministic crash skips the cascade") is the "use direct traceback" clause of both preloaded skills. | 1 line, same. |

**Counts:** RETIRE 3 · SALVAGE-then-RETIRE 0 · KEEP 0.
**Salvaged lines:** none (see "Salvage considered and declined").

## The blast radius, re-measured

`[M]` 2026-09-21, an independent re-run of the orchestrator's census with a positive
control, script kept at the session scratchpad (`census.py`): 2 704 files walked under
`/Users/rodrigo/git/nuclear/ORPHEUS` and the main agent's memory directory, excluding
`.git`, `docs/_build`, `.venv`, `scratch`, `__pycache__`; predicate = a line containing
the stem, the candidate itself excluded.

- `phase5a_moment_consuming_scatter_derisk` — **1** line.
- `r1_step_d_sphere_preconditioner_oscillation` — **1** line.
- `sn_sig_t_layout_drift_indexerror` — **1** line.
- Positive control `phase_f_step2_mesh_refinement` — **10** lines across **6** distinct
  files, which reconciles with `referrers.md`'s "6 referrers" (that census counts files,
  mine counts lines). The instrument can read positive.

Each candidate's single line is
`.claude/plans/archive/memory_distillation_2026-09-21/numerics-investigator/table.md`
(lines 73, 74, 75) — the distillation's own candidate list, which `referrers.md` already
classes as not a consumer. **No `repo`, `own-index`, `own-memory` or `other-memory`
consumer exists for any of the three.** In particular none of the three is listed in my
own `MEMORY.md` (grepped: zero hits inside the memory directory), so there is **no index
line to delete** for any of them.

**Nexus could not corroborate.** The brief says "agent-memory files are graph nodes"; in
the graph this session answers from, they are not. `query("agent-memory")` returns `[]`,
`query` on two candidate stems returns `[]`, and
`file_brief(".claude/agent-memory/numerics-investigator/r1_step_d_sphere_preconditioner_oscillation.md")`
returns `"... is not in the graph"`. The grep census above is therefore the whole
instrument, which is why it was re-run with a control rather than relayed.

## Referrer edits

One edit, shared by all three, and it is a **no-op by construction**:

- `.claude/plans/archive/memory_distillation_2026-09-21/numerics-investigator/table.md`
  lines 73–75 — **leave as history, no edit.** The file lives under
  `.claude/plans/archive/`, it is a dated record of what the 2026-09-21 distillation
  judged, and its rows are past-tense claims about candidates ("**0** referrers",
  "Phase-codename title"). Deleting the candidate makes those rows MORE true, not less.
  Per `vv-principles` anti-pattern #21's tense reconciliation, a past-tense history entry
  stays; only a present-tense claim would be a must-fix.

No other edit is owed anywhere: no index line, no cross-memory pointer, no docs
reference, no test or code reference.

## Why each file brings nothing forward — the trace, with its verification

### 1. `phase5a_moment_consuming_scatter_derisk.md`

Two sections, both superseded.

**(a) The LOAD-BEARING FACT (`Y_0^0 = 1`, unnormalized real harmonics, no `/√(4π)`
rescale).** This is a *convention*, and by Pattern 7 its home is its definition site, not
a copy in agent memory. It is at its definition site and on the consuming theory page:

- `orpheus/numerics/basis/spherical_harmonic_basis.py:58` (`Y_0^0 = 1` in the class
  docstring's math), `:127`, `:143`, `:191`, `:213` (the "no-prefactor normalisation"
  convention, carried by the basis object), `:591` (`Y[:, 0, 0] = 1.0`, the code).
- `docs/theory/foundations/spherical_harmonics.rst:66`, `:125`, `:260–265` — the
  definition of the no-prefactor convention.
- `docs/theory/methods/sn/cartesian_multid.rst:3111–3131` — a titled admonition, "The
  :math:`Y_0^0 = 1` convention — the scalar flux is read off", stating the full chain the
  memory file states (the `ℓ=0` moment IS the scalar flux, read "with no rescale", and
  "what makes the eigenvalue outer's scalar flux bit-identical to the full-angular
  `integrate_angular`").
- ERR-051's title carries the same convention ("no-prefactor SH convention") as a
  catalogued defect class, so a future violation has a catalogue entry too.

**(b) The methodology lesson and the de-risk evidence.** The whole Q-probe table the
memory file summarises is on the theory page, `docs/theory/methods/sn/cartesian_multid.rst:3218–3256`:
Q1/Q2/Q3 at 0 ULP, **Q2b** "vs INDEPENDENT Bell & Glasstone hand reconstruction (L11
structural-independence ground)" at `max rel = 3.4e-16`, and Q4's non-degeneracy row with
the same numbers (aniso max 0.49, `|Σ_{s,0} − Σ_{s,0}^T| = 0.1–0.18`). The prose beside
it even states the lesson in the memory file's own words: the independent reconstruction
"agreed at ~1.5 ULP (the expected floating-point distance for an independent reduction
order). This is the L11 structural-independence guard the bit-exact comparison alone
lacks."

The general form is rule-level, in two preloaded places: `vv-principles` anti-pattern #7
("NEVER treat 'two derivations agree' as proof — check *structural* independence", whose
**tell** is agreement that is too good) and #34 (α-normalised ASTs decide whether two
"independent" bodies are one), plus the "Bit-identity vs principled-equivalence" criterion
2 ("old-vs-new ULP distance is necessary, NEVER sufficient") and criterion 3 (the expected
drift is `reduction depth × ULP`); `instrument-doctrine` X4 is the always-on floor.

**Status verification.** The file claims "Phase 5a LANDED on origin/main (`93807aa` →
`63719a2`)". `git merge-base --is-ancestor` — **both ANCESTOR of HEAD**: `93807aa`
2026-06-07 "refactor(sn): factor scattering aniso onto a shared R·Λ moment→source
primitive", `63719a2` 2026-06-07 "docs(theory,sn): Phase 5c — document in-sweep moment
accumulation". The claim is true, and it is a claim about work that closed fifteen weeks
ago, which is the definition of archaeology.

The file's evidence pointer `derivations/diagnostics/diag_p5a_moment_consuming_scatter.py`
does not exist in the tree (the probe was consumed by the landing); the durable form is the
theory page's table.

### 2. `r1_step_d_sphere_preconditioner_oscillation.md`

This is the raw investigation note behind **ERR-050**, and the catalogue entry is a strict
superset — with one section the memory file cannot have, because it was written three days
too early.

- **Symptom.** The memory's oscillating-`keff` trace and "470× slowdown" are in ERR-050's
  "Empirical pre-fix evidence (R-1 Step D Probe B, 2026-05-19)" code-block, including the
  slab/sphere × default/identity comparison that is the memory's "Confirming evidence"
  table.
- **Root cause.** ERR-050's "Mechanism" states the same chain (the M-M Carlson coupled-pole
  closure, Hébert §3.9.4 Eqs. 3.432–3.435, `rhs(1)` returning the cold-frame default because
  GMRES feeds residual vectors with `history_depth == 0`, slab immune for want of a
  curvilinear pole) and adds a four-item "How it hid".
- **⚠ The memory file is now stale-wrong about production.** ERR-050's **Status** is
  "**CLOSED via structural supersession**": on 2026-05-22 Phase 1.2 (`c93355c`,
  **ANCESTOR** of HEAD, "refactor(sn): unify sweep+matvec through M-M's psi_half_seed
  strategy") made `StreamingCollisionOperator.solve` a pure function taking an explicit
  `initial_guess=` — *"the silent-fallback path no longer exists"*. Every file:line in the
  memory's "Files referenced" block describes the pre-`c93355c` tree
  (`solver.py:~608`, `iteration.py:596–601`, `operator.py:2570` "`carlson_seed = rhs(1)`").
  A future session reading the memory file would be diagnosing a mechanism that has not
  existed since 2026-05-22.
- **The fix it proposes did land, elsewhere.** `orpheus/sn/solver.py:959–961` now reads
  `preconditioner = lambda q: q  # noqa: E731` under the comment "explicit identity —
  issue #200 tracks the face-preconditioner re-enablement".
- **The regression catch it recommends did land.**
  `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:107` carries
  `pytest.mark.catches("ERR-050")`.
- **The live successor exists and is indexed.** `curvilinear_inverse_seed_taxonomy.md`
  (listed in my `MEMORY.md` §3) carries the current state of exactly this question —
  including "curvilinear GMRES ships the **IDENTITY** precond (`_within_group_krylov`
  solver.py:332 …), with #200 tracking the real one" and the sphere's seed-lag exit via
  route (a).
- **Issue status, verified.** `gh issue view 200 --json state` → **OPEN**
  ("SN: block-inverse preconditioner for Krylov on the typed AngularFlux algebra"); already
  carried by my `MEMORY.md` §2 open list. `gh issue view 203 --json state` → **CLOSED**
  ("numerics: KrylovAcceleration default precond requires CAP_STATELESS_INVERSE"), matching
  ERR-050's "Issue #203 closed by supersession".
- The one item with no home elsewhere is the file's own aside, "Long-term:
  `KrylovAcceleration` may want a sentinel like `preconditioner="identity"`". It is not
  worth carrying: the general principle is ERR-050's **Lesson** ("default values for
  behavioural parameters MUST either advertise their preconditions in the type system OR
  require explicit caller choice"), the production spelling is already explicit with a
  reason comment, and a stringly-typed `"identity"` sentinel is `coding-elegance`
  anti-pattern #4 in any case.
- The probe the file names, `derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py`,
  is gone from the tree — promoted, per ERR-050's "Which test catches it".

### 3. `sn_sig_t_layout_drift_indexerror.md`

**ERR-055** (`docs/theory/verification/error_catalog.rst:4901–4931`) is the durable home
and is strictly richer: it names all six affected tests by nodeid, gives the root cause
with the same file:line evidence, and its three-item **Lesson** is the memory's
methodology list nearly word for word —

1. "`coding-elegance` Pattern 7 applies to test fixtures, not just production … A test that
   calls an internal `_helper` directly with bare arrays is a Pattern-7 landmine."
2. "A layout/convention migration is incomplete until the direct-call tests are migrated
   too."
3. "ng=1 degeneracy hides axis-swap convention drift … another instance of cross-cutting
   hygiene rule H1 (1-group degeneracy) operating at the data-layout level."

ERR-055 also carries the masking twin the memory names (`test_unified_matvec_cylinder.py`
used the correct `(ng, n_cells, 1)` layout) under "How it hid".

The memory's fourth bullet — "Deterministic IndexError ⇒ skip the 8-step isolation cascade;
the traceback names file:line directly" — is the "When NOT to use" clause of both skills
that preload with me: `probe-cascade` SKILL.md:25 and `numerical-bug-signatures`
SKILL.md:38, "use direct traceback / `nexus-debugging` instead". Carrying it again in
memory would restate a preloaded skill, which the distillation standard forbids.

**Status verification.** `6cfdfd4`, the commit the file blames for the layout flip, is an
**ANCESTOR** of HEAD (2026-05-15, "perf(sn): CollisionCache (N, ng, nx) layout flip …"). The
file's own closing note says its two diagnostics were deleted post-fix; confirmed, neither
is in `derivations/diagnostics/`.

## Salvage considered and declined

One line was a genuine candidate and I am recording the decision rather than a bare
"considered and rejected", so the next reader does not re-open it:

> "For ANY bit-exact 'two paths agree' de-risk: add one probe whose math is hand-written
> from the textbook (different reduction order, different intermediate naming) and assert
> ~few-ULP, **NOT** ULP=0. Bit-exact-where-shared + few-ULP-where-independent is the honest
> signature; bit-exact everywhere is a tautology waiting to hide a shared upstream bug."
> — `phase5a_moment_consuming_scatter_derisk.md`, lines 42–49.

**Declined as a salvage**, because the SALVAGE bar is "nowhere else" and each half is
somewhere else, all of it preloaded: the "too-good agreement is the tell" half is
`vv-principles` anti-pattern #7's **tell** ("agreement at 1e-39 (ERR-032)") and #34; the
"an honest independent reduction reads a few ULP, not zero" half is the bit-identity
section's criterion 3 (drift is `reduction depth × ULP`, dimensionally explainable); the
refusal posture is criterion 2 and `instrument-doctrine` X4. The concrete instance is on
the theory page (§ above). Writing it into `lessons.md` would restate a preloaded skill,
which is the thing the distillation removed.

**Offered instead as an optional, separate uplift** (NOT part of this audit's deliverable,
and I do not think it is required): if the orchestrator ever wants the *positive*
discriminator stated once in the ceiling skill, the natural site is `vv-principles`
§"Bit-identity vs principled-equivalence", one clause on criterion 2 — *the independent
reference is expected to read a FEW ULP; a 0-ULP reading against a supposedly independent
probe is evidence the probe is not independent.* Its founding case would be the Q2b row
above. It is a skill edit, not a memory edit, so it is out of scope here.

## Findings outside the ask

Two, both read-only observations; neither blocks the retirements.

1. **A stale test path in a tracked doc.** `docs/theory/verification/error_catalog.rst`
   cites `tests/sn/test_krylov_curvilinear_precond_safety.py` twice — line 4500 (ERR-050,
   "Which test catches it") and line 4752 (ERR-053) — but the file lives at
   `tests/sn/solve/test_krylov_curvilinear_precond_safety.py`; `ls
   tests/sn/test_krylov_curvilinear_precond_safety.py` → "No such file or directory" (the
   only artefact left at the old path is a stale `__pycache__` entry). The marker itself is
   intact (`catches("ERR-050")` at `tests/sn/solve/…:107`), so this is a prose-side dead
   reference, not a coverage gap — exactly the class `dead_references` exists for, and the
   class that produces no build warning at any severity. Worth a one-line fix in whichever
   commit next touches the catalogue.

2. **Agent-memory files are not in the Nexus graph** (see "The blast radius, re-measured").
   If the campaign's later dispatches rely on `context()` over memory nodes, that premise
   needs re-checking before it is briefed again.

## NEEDS:

- **Nexus corroboration of the inbound edges was unobtainable**, for the reason in
  finding 2: the candidates are not graph nodes in the current build, so `context()` /
  `file_brief()` cannot speak to their inbound edges. The verdicts rest on the
  positive-controlled grep census above, which agrees with `referrers.md` on all three
  candidates and on the control.
- **No other blockers.** Every status line cited here was verified this session:
  `git merge-base --is-ancestor` for `93807aa`, `63719a2`, `6cfdfd4`, `c93355c`,
  `2c634ab`, `a29ab2d` (all ANCESTOR); `gh issue view` for #200 (OPEN) and #203 (CLOSED).
