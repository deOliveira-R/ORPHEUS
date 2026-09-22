# Proposal — explorer topic-file blast-radius audit (2026-09-21, HEAD d3794a2f)

Owner: explorer. Read-only on tracked files; this is the one deliverable. Every
status line below was re-verified this dispatch (`gh issue view N --json state`,
`git merge-base --is-ancestor <hash> HEAD`), never taken from a memory.

## 0. Verification ledger

Issues (`gh issue view N --json state,title`, 2026-09-21): #459 CLOSED · #331 CLOSED ·
#208 CLOSED · #168 CLOSED · #364 OPEN · #340 CLOSED · #325 CLOSED · #326 CLOSED ·
#280 CLOSED · #337 CLOSED · #281 CLOSED · #133 CLOSED ("Phase 5+: continuous-µ
specular sphere production wiring") · #137 CLOSED · #100/#103/#132 OPEN.

Hashes (`git merge-base --is-ancestor H HEAD` → all YES): 7aae9bf1 (2026-08-19),
4207c1d6 (2026-09-11), 7e9b6210 (2026-09-08), 5c39226b (2026-08-14), 82c0d441
(2026-09-08), 9ff1ee3f (2026-07-25). No `track-b`/`adjoint`/`consumers` branch
exists locally or remotely (`git branch -a | grep -i ...` → 0).

Referrer census, re-run: predicate = a line containing the stem, `grep -rn -I`
over the repository minus `.git`, `_build`, `.venv`, `scratch/`, `node_modules`,
excluding the file itself; result per stem is EXACTLY the set `referrers.md`
lists (see §2 per file). The `~/.claude/projects` tree was NOT counted: it holds
session transcripts (JSONL), which quote system prompts and are not a memory.

Nexus: `query("<stem>")` → `[]` for every stem; `file_brief` result recorded in
§2. (Whether agent-memory files are graph nodes in THIS build is stated there.)

Referrer census, two forms, both with a positive control:
- stem form (`grep -rn -I <stem>` over the repository minus `.git`, `_build`,
  `.venv`, `scratch/`, `__pycache__`, plus the main memory dir
  `~/.claude/projects/-Users-rodrigo-git-nuclear-ORPHEUS/memory/`; the file itself
  excluded). `[M]` positive control `phase_f_step2_mesh_refinement` → 6 files. Result
  per stem = exactly the set `referrers.md` lists; main memory → 0 for all five.
- **wikilink-slug form** (`grep -rn -F '[[<name-slug>]]'`, the hyphenated `name:`
  frontmatter field — the stem predicate CANNOT see it). `[M]` positive control
  `[[affine-operator-split-convention]]` → 2 files. Result: ONE referrer the
  orchestrator's census missed —
  `.claude/agent-memory/explorer/census_predicates_bound_reference_and_activation_traceback.md:66`
  → `[[sn-solve-exit-and-reflective-default]]`. The other three slugs → 0; the
  phase5 file has no slug (`name:` is a free-text title) and its title → 0.
- title form (`grep -F` on each index title) → only the index lines themselves, plus a
  HOMONYM: `.claude/agent-memory/cross-domain-attacker/problem_solution_split_frames.md`
  is the attacker's own file (`problem-solution-split-frames`), not a referrer.

Nexus (build `0e1d6079`, 2026-09-21, ancestor of HEAD, `changed_since_build: false`):
`file_brief("orpheus/sn/solver.py")` resolves (Python paths bridge); `query("field_algebra")`
resolves the doc node `std:file:theory/foundations/field_algebra`; but `query(<stem>)`
→ `[]` for all five, AND for the two positive controls `phase_f_step2_mesh_refinement`
and `affine_operator_split_convention` (both known-cited memory files), AND
`graph_query("file -contains-> * WHERE name=agent-memory*")` → `[]`. **`.claude/agent-memory`
is not indexed in this build**; the brief's "agent-memory files are graph nodes" and
the distillation standard's "Nexus indexes `.claude/agent-memory` ALONGSIDE `docs/`"
(2026-06-22) are claims with a vintage that this build refutes. The grep census above
is therefore the inbound-edge evidence (see NEEDS).

## 1. Verdict table

| file | verdict | reason (one line) | referrer edits |
|---|---|---|---|
| `flux_torsor_vs_cone_inventory.md` (vintage 2026-08-19, `5f6f6522`) | **RETIRE** | Maps retired code — `orpheus/transport/fields/_flux_role.py` absent; `[M]` `grep -rln 'FluxDisplacement\|class FluxRole\|affine_combination\|sibling_of' orpheus/` → 0 files (5 hits in `tests/`, all docstring prose about the retired design); its forward half (the cone-side footholds: `is_positivity_preserving` 0 readers, `realizer.py` ZeroFlux refusal, `TestCrossSectionConeAlgebra`, `power_iteration` ray-normalisation, the incumbent's concessions) is VERBATIM in `lessons_archive.md` §L-034 (lines 1281–1300); the overturned design's record is `docs/theory/foundations/field_algebra.rst` §`cone-the-overturned-affine-design` (line 1421); #208/#331/#168 all CLOSED. | `MEMORY.md:44` — delete. `table.md:111` — no edit (a dated judgment, true as written). |
| `identity_step_census_durable.md` (vintage 2026-09-11, `22ec2c67`) | **RETIRE** | The design it was to shape IS the tree: `MaterialMesh._contractibility_key`/`_identity_key`/`__eq__`/`__hash__`/`same_phase_space` (`orpheus/transport/mesh/material_mesh.py:363–438`) with a code comment block (`:305–333`) that restates the two-key design, the #459 history, AND resolves the file's "contested closure row" (R-cc8: identity INCLUDES the closure class + clamped order, contractibility EXCLUDES them — both sides were right about different keys); `SNProblem.scattering_order` clamped at `orpheus/sn/problem.py:325–330`; `[M]` `_GEOM_CACHE_INTERN` → 0 hits in `orpheus/`; the method lessons (clamp-vs-pad per entry, `La13511Case` homonym, `is`-on-Optional, dataclass hash lies) are digest L-041/L-042/L-043 + archive §L-043 (line 1579). #459 CLOSED. | `MEMORY.md:52` — delete. `table.md:112` — no edit. Its own line 58 (→ `problem_solution_split_census.md`) dies with the file (both retire). |
| `problem_solution_split_census.md` (vintage 2026-09-08, `9c93c97d`) | **RETIRE** | A `reference`-type pointer to an UNTRACKED scratch memo (`scratch/_consumers/explorer_problem_solution_census.md`, exists today, `??`), and 4 of 4 "durable" claims are false at HEAD: (1) the hub now OWNS bound operators — `SNProblem.loss_kernel_gauge`/`system`/`pencil`/`eigen_posing`/`fission` are `cached_property`s (`problem.py:1101–1210`); (2) `[M]` `_coll_cache\|_pole_mirror_cache\|_geom_cache` → 0 hits in `problem.py`/`solver.py`/`material_mesh.py`; (3) `scattering_order` lives on `SNProblem` (`problem.py:328`), not the solver; (4) `is_same_phase_space` renamed/repaired (#459 CLOSED). Its "How to apply" (`git log --since` the audit date) is M-2. | `MEMORY.md:51` — delete. `table.md:113` — no edit. `identity_step_census_durable.md:58` — no edit (referrer retires). |
| `sn_solve_exit_and_reflective_default.md` (vintage 2026-08-17, `144cdf51`) | **RETIRE** | §1 remedied — `[M]` `grep -n 'return Solution(\|return _package_solution(' orpheus/sn/solver.py` → 5 `_package_solution`, 0 direct (archive §L-030, line 1097, holds the story); §2 (unset BC ⟹ reflective at hub level; one-face-declared ⟹ the REST get reflective, not the entry's vacuum) is stated in THREE code docstrings — `resolve_boundary_conditions` (`orpheus/transport/method.py:229–234`, "``None`` defaults to ``BC("reflective")`` (the infinite-lattice / eigenvalue convention, uniform across methods)"), `_apply_default_bcs` (`orpheus/sn/solver.py:126–128`, "unchanged when it already carries ANY explicit BC … all-or-nothing"), `_as_problem` (`solver.py:163–171`, "fills faces only when the declaration carries no explicit BC … unset faces then resolve to the SNProblem-level reflective default") — a memory copy is a second copy that drifts; §3 (tangential bucket ⟹ half-range functionals miss it) is the module docstring (`orpheus/numerics/spaces/angular_trace_space.py:169`) + `docs/theory/foundations/spaces.rst:2080` + archive §L-010 (line 210); §4 is a perishable negative (`[M]` still 0 consumers of `.boundary_flux` in `orpheus/`/`examples/` today, but M-2: re-derive, never carry). | `MEMORY.md:36` — delete. `table.md:47`, `:114` — no edit. **`census_predicates_bound_reference_and_activation_traceback.md:66`** — re-point: replace `[[sn-solve-exit-and-reflective-default]]` with `L-030 (count the seam's return sites first; archive §L-030)` — the wikilink cited the "three exit sites" finding as a sibling of L-013's swap-and-run, and L-030 is that finding's durable home. |
| `phase5_mu_resolved_primitive_inventory.md` (vintage 2026-04-28, `f149c03f`) | **SALVAGE, then RETIRE** | Every referent moved: `orpheus/derivations/peierls_geometry.py` no longer exists (`bda76faf` "reorganize into common/discrete/continuous"; now `orpheus/derivations/continuous/peierls_nystrom/geometry.py`), so all ~40 `file:line` cells are dead; the campaign both LANDED (`_chord_tau_mu_sphere` exists at `geometry.py:2604` — the file's "extract" recommendation was executed; #133 CLOSED "Phase 5+: continuous-µ specular sphere production wiring") and was then RETREATED from (`docs/theory/references/trajectory_resolvent.rst:3781` "singularity that killed Phase 5", §`peierls-continuous-mu-retreat`); a design sketch for a dead campaign. ONE method correction (Risk #2) is nowhere else — see §2. | `MEMORY.md:49` — delete. `table.md:115` — no edit. |

Counts: RETIRE 4 · SALVAGE-then-RETIRE 1 · KEEP 0. Index lines deleted: `MEMORY.md`
36, 44, 49, 51, 52 (five). One own-memory referrer re-pointed (the wikilink the
stem census missed). No `repo` referrer outside the distillation's own archived table;
no other-agent or main-memory referrer for any of the five.

Optional (orchestrator's call): append one line to `table.md`'s "Orchestrator's notes
at apply time" recording the retirement hash, so the archived table's "not retired
here … a separate pass" sentence points at its outcome.

## 2. Salvaged

### From `phase5_mu_resolved_primitive_inventory.md` → digest `lessons.md`, new entry L-046 (or a rider on L-029 — orchestrator's placement)

Source lines (Risk #2, verbatim):

> The lesson: the explorer rule "the most-connected primitive is the canonical one" is
> FALSE here — `compute_P_esc_mode` (degree 131) is the rank-N Marshak primitive, NOT
> the canonical no-Jacobian form Phase 5 wants.

Still grounded at HEAD: `orpheus/derivations/continuous/peierls_nystrom/geometry.py:3984`
documents `compute_P_esc_mode`'s `(ρ_max/R)²` surface-to-observer Jacobian and `:5639`
names the "canonical no-Jacobian Mark-Lambert primitive" as a DIFFERENT helper. `[M]`
`grep -rn -i 'most.connected\|degree.*canonical\|hub.*convention'` over `lessons.md`,
`AGENT.md`, the two nexus skills and `.claude/rules/` → only the `god_nodes` tool
descriptions; the caveat exists nowhere.

Distilled entry (≤ 5 lines):

- **L-046** Graph DEGREE ranks connectivity, not canonicity: the most-connected
  primitive is the one most CALLERS chose, and it may encode a convention (a
  Jacobian, a normalisation, an index domain) the question at hand must NOT inherit
  (`compute_P_esc_mode`, degree 131, carried the rank-N Marshak `(ρ_max/R)²` Jacobian;
  the low-degree no-Jacobian sibling was the template). Before using a hub as a
  template, read the convention it encodes against the derivation of record and
  grep the siblings for "canonical". → M-6.

Archive: the three source lines above go to `lessons_archive.md` under a new
`## L-046` heading (the war story), so the digest entry has its pointer.

Nothing else in the five files is salvage: every other line is a `file:line` map
(dead or re-derivable), a count with a vintage, a negative claim (perishable by M-2),
or a design fact the tree now states in its own docstrings.

## 3. Discrepancies with the brief / referrers.md (findings)

1. `referrers.md`'s stem predicate misses the hyphenated `[[name-slug]]` wikilink form
   (filename `sn_solve_exit_…` vs `name: sn-solve-exit-…`); one own-memory referrer was
   missed. Recommend the orchestrator's census add the slug form for every owner.
2. "agent-memory files are graph nodes" is FALSE for build `0e1d6079`: two known-cited
   memory files are absent from the graph and `graph_query` finds no agent-memory file
   node. Either the build excludes the tree or the claim's vintage (2026-06-22) predates
   a config change; the blast radius here rests on grep, not the graph.
3. `problem_solution_split_census.md`'s headline ("the hub owns … no bound operators")
   has FLIPPED at HEAD (five operator-valued `cached_property`s on `SNProblem`); the
   `MEMORY.md:51` hook restates the flipped claim present-tense.
4. `phase5_…`'s referent module does not exist at its cited path.
5. **Tree moved during the audit (M-2 / L-007).** At close, `git status --short | grep -v '^??'`
   → ` M orpheus/sn/solver.py`, ` M docs/theory/verification/error_catalog.rst`, mtimes
   23:20:14 / 23:20:46 vs the brief's 23:12:54 — a sibling session's edits, NOT this
   dispatch's (every command here was a read or a write into `scratch/`). The solver hunk
   is docstring-only (`:3194–3200`, "One production caller" → "Two production callers")
   and touches no line this proposal cites (`:120–171` sit above it; the return-site
   counts are hunk-independent). Line numbers are stamped "at final read; tree moving".
   Never `git checkout`/`restore`/`stash` those two paths — they carry live work.

## NEEDS:

- Nexus inbound-edge check on memory files: unavailable in build `0e1d6079`
  (`.claude/agent-memory` not indexed; positive controls `[]`). Used the two-form grep
  census with positive controls instead. If the orchestrator wants the graph check, the
  build must index that tree first.
- Placement ruling for L-046 (new digest entry vs rider on L-029) and whether the
  archived `table.md` gets the one-line retirement note.
- UNVERIFIED: whether any test pins the §2 one-face-declared hazard (`Mesh2D(...,
  bc_xmin=BC("vacuum"))` + `boundary_condition="vacuum"` ⟹ other faces reflective). Not
  searched; if none exists it is a gate owed by the docstrings, out of this audit's scope.
