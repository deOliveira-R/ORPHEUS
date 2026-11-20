# W0 specification, part A: units insulation, dead, roles, defects, posing

Status: DRAFT for `spec-check` (`units`, then `check_only`) and the user's checkpoint. Written 2026-11-20 on `refactor/ontology` at `b68f8b8c`; AMENDED the same day by the orchestrator's decisions and the ten check files (the record is section 9; the unit sections carry no changelog). Author: the Draft stage (executor). Nothing in this file is a ruling; section 8 lists every choice the inputs did not settle.

Marker use follows `articulation` §5: `[M]` measured (the command beside it, or the measurement memo it was taken from), `[R]` reasoned and not measured, `[HYPOTHESIS]`. Every measurement cited here was taken by an earlier stage on tree `d2008d0e` or `refactor/ontology` before this file was written; those memos are untracked under `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/ontology/w0/` and a builder's worktree has no `scratch/`, so every row, path list and line number a builder needs is INLINED in this file. A line number holds only at the base hash named above and only before an earlier unit integrates (section 7); every line number is therefore given with a quoted anchor or a symbol, and the symbol decides.

## 1. Scope and bindings

### 1.1 What this specification covers

Five units of wave W0 of `.claude/plans/ontology_chart.md`: `insulation` (section 2), `dead` (section 3, four sub-units dead-A to dead-D), `roles` (section 4), `defects` (section 5, with its sub-unit `angular-mass`, 5.7, the root fix of #708 added by the third pass), `posing` (section 6). The units `registry`, `digest` and `memo` are another specification's. A precondition unit, `quadrature-ties` (2.0: the root fix of the rounding certificate that refuses an exact value on a float64 tie, so that `LEGENDRE.on(a, b).gauss(n)` is total), is a sixth unit; it integrates before every other unit and in particular before `insulation`, which calls the fixed code. Section 7 gives the order in which the units integrate, section 8 the decisions for the orchestrator.

### 1.1b What a builder's brief for one unit is made of

The build stage loads no CLAUDE.md, no rule and no skill. The brief for unit X is: section 1 whole (scope, bindings, conventions 1.3, the audit table 1.4, the shared gate 1.5, the unmeasured list 1.6); the section of unit X (2 for `insulation`, 3 for each dead sub-unit with its Appendix D-A to D-D, 4 for `roles` with Appendix R, 5 for `defects`, 6 for `posing`); section 7 (its place in the order and its conflicts); and the decisions of section 8 it cites. A dead sub-unit's brief additionally contains section 3 whole (the premise checks, the common procedure, the error-catalogue gate, the kept rows, the cycles). The brief for sub-unit `angular-mass` is section 1 whole, section 5.7 (self-contained: it needs no other part of section 5), section 7 and the items 35 and 36 of section 8.0 and the items 8.28 and 8.29 of section 8.4. The brief for the precondition `quadrature-ties` is section 1 whole, 2.0 whole (self-contained: it needs no other part of section 2), section 7, the items 37 to 41 of section 8.0, the rulings 8.42 to 8.46 of 8.4 and the open question 8.47 of 8.5.

### 1.2 The binding decisions (the orchestrator's, under ruling O2; not re-litigated here)

1. `orpheus/derivations/` is closed (rulings O8, O18/A2). It is never redesigned. A unit re-points mechanically the `derivations` imports its change breaks. No compatibility shim stays in production.
2. Insulation. The chord is `ConcentricPartition.chord` (bit-identical today). The composite Gauss-Legendre rule is spelled with the production quadrature package, `LEGENDRE.on(a, b).gauss(n)` per panel, concatenated inside `cp`; the three bit-exact pinned rows that move by 1 to 30 ulp are re-baselined with the reason in the commit. The Bickley `Ki_3` evaluator is copied into `orpheus/numerics/` as a special function with the SAME interpolant, so that no cylinder number moves; a more accurate evaluator is a filed issue, not W0. Production's imports of `orpheus.derivations` go to 0 and a gate forbids the edge from every production package in its OWN table (plain `FORBIDDEN_EDGES` rows redden 95 test files through `test_test_layers`).
3. Dead. A row is retired, moved to the test side or kept as `dead/reviewed.json` says (58 retire, 47 move, 41 keep; only these rows, never the unread census rows), with two departures: four of the 47 moves are deferred in W0 because their body reads private state (3.1 rule 2), and four of the 58 retires, rows 134 to 137, stay in production as catalogue data (3.0 item 2, kept for W8a). So the table W0 acts on is 54 retire, 43 move, and 49 kept in production (41 + 4 + 4). The 6 tests that reach `rayleigh_quotient_iteration` through a `getattr` string and skip when it is missing are deleted in the same sub-unit. Retiring `select_quadrature` migrates the ERR-081 marker to the surviving catalogue test and records ERR-118 as retired with its code in `docs/theory/verification/error_catalog.rst`. `thermal_hydraulics`, `kinetics` and `fuel` are NOT touched, so #713 is corrected as an issue and not fixed.
4. Roles. `BlockRole`, `SystemRole`, `_BlockRoleMeta`, the three marker classes, `stamp_boundary_role` and the four string fixtures retire with no replacement; the derived block pattern is a later wave's. Posing item u1.6 is inside this unit.
5. Defects. #709 is a refusal tagged W3, #710 a refusal tagged W1 on the `FromMeasure` and `FromQuadrature` arms, #711 a refusal tagged W8a. Each refusal is a `raise` that survives `python -O`, carries `ELEGANCE-DEBT[guard]` with the issue and the wave that retires it, and lands with the constructible input that it rejects. #708 is in W0 as a ROOT FIX, not a refusal (the orchestrator's ruling of 2026-11-20; sub-unit `angular-mass`, 5.7): one definition of the angular mass (the correctly rounded sum of the stored weights, owned by the angular section), `AngularLift.total_weight` and the fifth spelling in `radial_characteristic.py` read it, `_combine` injects through the section, and the mixed `ScalarSourceSink + AngularSourceSink` arm is deleted with `CONTAINS` (so `Field`'s existing same-class gate refuses the sum: no new guard, no tag, no wave). #392 is closed as already fixed.
6. Posing. Items u1.2, u1.3 and u0.3 only.
7. Branch mechanics (B9, O12). A unit branch is cut from `refactor/ontology`, rebased onto it and fast-forwarded; integrations are serial; a unit's file set is its call-site set (production, the tests that consume the changed interface, `derivations/`, `tools/`); an issue is closed by hand with the integrating hash.
8. The amendment of 2026-11-20 (the orchestrator's decisions, each inlined where it acts; a builder does not read the file that holds them): ONE retirement engine, `tests/gates/test_retirements.py`, with rows modules under `tests/gates/_retired/` and three row types (1.5); the engine binds by definition site and import resolution; every walk takes a root path and a mutation is a copy of the tree or an in-process monkeypatch (1.3); `tests/gates/test_layer_imports.py` is edited in the order roles, insulation, dead-B; a closure member that `dead/reviewed.json` does not rule is not touched (3.1 rule 3); a known defect is fixed before work is built on it (2.0).
9. The exact-tie unit (the orchestrator's rulings of 2026-11-20, 8.42 to 8.46): the root fix of `LEGENDRE.on(a, b).gauss(n)` is a precondition unit, `quadrature-ties` (2.0), that integrates before every other unit; an exact value on a float64 rounding tie is decided in rational arithmetic and rounded half to even; the closed-form recurrences are generic over the number type with no new stored field; the unit's gate is the investigation's 12 rows promoted into `tests/gates/numerics/test_generating_measure.py` plus the witness of the unequal single-boundary arm; an ERR entry is owed.

### 1.3 Conventions that hold in every unit

- **Tests run serially**: `.venv/bin/python -O -m pytest -p no:cacheprovider <file>`. Never a full-suite run inside a unit (about 75 minutes); the unit's gate files, the files it edits and the files named in its done-when. Python is `.venv/bin/python`.
- **`grep` in this environment is ugrep**, which can return a silent zero for an anchor (`^`, `$`) inside an alternation group. Every completeness claim in a done-when is therefore a Python `re` + `pathlib.rglob` (or `ast`) probe: write it to a file under a quoted heredoc (`<<'PY'`), and run a positive control (a symbol you know exists) beside every zero. Membership questions are answered by `ast`, never by a grep through a line window. Use `git grep` for tracked-file text searches (it reads only tracked files; it has no `\b`, so use `-P` with a lookbehind, or Python).
- **Stage files by path.** Never `git add -A`. Write a commit message to a file and use `git commit -F`; backticks in `-m` are command substitution. End the message with the attribution lines of the session. Close issues by hand with the integrating hash.
- **Never revert a mutation with `git checkout` / `git restore` / `git stash`** on a file with uncommitted edits. Copy the pristine file aside BEFORE the first mutation, mutate the copy or monkeypatch in-process, and `diff -q` against the copy after.
- **The layer gates** (`tests/gates/test_layer_imports.py`, `tests/gates/test_test_layers.py`) read `FORBIDDEN_EDGES`. A test file under `tests/gates/<dir>/` may import only what the package `<dir>` may import, and a `tests.gates.<other_dir>` helper is judged as the layer of `<other_dir>`. A new test helper lives in the directory of the layer it needs; a MOVED symbol's helper lives in the gate directory named for the package the symbol left, so that no consumer's package edge changes (3.1). Run `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_layer_imports.py tests/gates/test_test_layers.py` after every move.
- **Raises survive `python -O`.** A refusal is a `raise`, never an `assert`.
- **A mutation names its mechanism.** Every walk of a gate takes a root path parameter. A mutation of a source file is spelled: copy the tree (or the one file, at the same relative path) to a temporary directory, mutate the copy, call the walk on that root, require red; for a cold-interpreter leg, load the mutant copy with `importlib.util.spec_from_file_location`; for a behavioural gate the mutation is an in-process `monkeypatch`. "A temp copy" with neither mechanism is not a mutation instruction here. Every gate in a unit section states its input, its first red in today's tree (or says plainly that its first red is a test-only input) and its mutation with the mechanism.
- **A done-when grep excludes the rows modules.** `tests/gates/_retired/` names the retired symbols by design, so every done-when `git grep` over a retired name carries `':(exclude)tests/gates/_retired'`, and every count is restated over the full name list it claims (roles: 12 names, 117 hits in 17 production files, 275 in 19 test files).
- **Harness texts.** A retirement that makes a sentence in `CLAUDE.md`, `.claude/rules`, `.claude/skills` or an agent definition false edits its SOURCE under `docs/development/` and regenerates with `.venv/bin/python -m tools.harness` (B8). Generated files (`.claude/skills/vv-principles/error_index.md`, `docs/theory/verification/matrix.rst`) are regenerated ONCE at integration by the orchestrator, never by a unit builder (section 7).
- **The docs build.** `sphinx-build -E -W docs docs/_build/html` must stay clean. After any retirement, rebuild first (the Nexus graph is rebuilt by the docs build), then run `.venv/bin/nexus dead-references` and read the rows for the retired names. A false DEAD announces itself; a false ALIVE is silence, so a grep over the text sources (item B.6 below) accompanies it.

### 1.4 The retirement audit, item by item, with the command that discharges it

`.claude/skills/retirement-audit/SKILL.md` binds every retirement and move. A builder does not load it; the 25 items follow, each with its command. A unit's section says which items apply and what each finds in that unit.

| item | what it asks | command |
|---|---|---|
| A.1 | graph callers; necessary, not sufficient | `.venv/bin/nexus callers <symbol>` and `.venv/bin/nexus impact <symbol>` |
| A.2 | text surfaces: code, tests, `docs/` | `git grep -n -P '(?<![A-Za-z0-9_])NAME(?![A-Za-z0-9_])' -- orpheus tests tools examples docs .claude`; plus `sphinx-build -E -W -n` as a pre-edit-versus-post-edit set diff over the pages edited (nitpicky mode reports every dead `:class:`/`:func:`/`:meth:`/`:attr:` role on an `.rst` page) |
| A.3 | direct constructors of a guarded type | `git grep -n -P 'T\(' ...` per type; every `T(...)` caller, not only the factory |
| B.4 | docstring surface and `.claude/` | `.venv/bin/nexus dead-references` after a docs rebuild; `git grep -n -P '(?<![A-Za-z0-9_])NAME(?![A-Za-z0-9_])' -- .claude CLAUDE.md`, sorted by tense (E.19). `.claude/agents`, `.claude/skills`, `.claude/rules`, `CLAUDE.md` are edited through `docs/development/` and `python -m tools.harness`; `.claude/agent-memory/` and `.claude/plans/` are NOT in a unit's edit set (the orchestrator commits memory at close-out without curating it) |
| B.5 | a math symbol's spellings (ASCII, Unicode, LaTeX) and its number | grep the three spellings of any symbol the retired object's docs carry; read the numbers the claim carries |
| B.6 | a name inside a string (`getattr`, `hasattr`, `setattr`, `monkeypatch.setattr`, `__all__`, a dict key, an import fixture) | a Python `ast` probe over `Constant` strings equal to the name (non-docstring), over production, tests, tools, examples; `git grep -n -E "['\"]NAME['\"]"` as the cross-check. The shared gate of 1.5 holds it |
| B.7 | the concept, not only the symbol (paraphrase) | read the docs table rows and column headings of the pages in the tense table; grep hyphen and space variants of each name |
| B.8 | the paragraph that justified the retired thing | read the paragraph the retired name sits in; re-verify against the replacement |
| B.9 | a labelled equation is an API | `git grep -n -F ':eq:`LABEL`'` for each `:by:` or `verifies(...)` label a retired row's markers name; the tense tables list `:by:` targets (`eq/ref`) |
| B.10 | a cross-reference you add or lean on | read the target of every `:ref:`/`:doc:`/`:eq:` you add |
| B.11 | a message string is an API once a test pins a substring | `git grep -n -F '<shortest distinctive fragment>'` over `tests/` |
| C.12 | test migration: behavioural tests rewire, API-smoke tests delete | per row, the review's test fate (the appendix tables) |
| C.13 | marker migration: `catches(...)` / `verifies(...)` move with the test | `git grep -n -E 'catches\(|verifies\(' -- tests` for every marker the appendix table names; `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_error_catalogue_reconciles.py tests/gates/test_vv_harness_audit.py` |
| D.14 | demotion: a re-pointed gate may compare a value with itself | for every rewired comparison, one in-process mutation of the successor's body: still green means the gate never pinned anything |
| D.15 | promotion: a surviving gate may have become stronger; its docstring advertises the weaker claim | re-read the docstring of every gate whose subject was retired or single-sourced |
| D.16 | single-sourcing a duplicate demotes every gate that compared its copies; the demotion is correct | for each such gate: what input could still make the two sides differ; the replacement pin; the gate's own docstring says so |
| D.17 | a re-baselined pin: the first pin you reach for may be blind | one in-process mutation at the OLD value, which must redden the pin |
| D.18 | retiring a duplicate makes whatever kept the copies equal load-bearing | name the mechanism and its witness; write the witness in the same commit when none exists |
| E.19 | sort surviving references by tense; past stays, present is a must-fix | the tense tables of the appendices (`PRESENT` rewrite, `past` stays) |
| E.20 | routing through the algebra raises "has the verb" to "is the type" | read each rewired call site that now goes through `@` or `+` |
| E.21 | a retirement onto a shared helper moves a raise's provenance one frame out | any frame-asserting gate re-pointed at a helper also asserts the caller frame or the message owner |
| E.22 | mass deletes are retirements too | grep the module or script NAME, including string imports and `textwrap.dedent` bodies |
| F.23 | a length-changing rename owes an `.rst` underline scan | applies to no unit here (there is no rename); the moved names keep their spelling |
| G.24 | retire leaves first, in one commit | the closure column of the appendix tables; retire in leaf-first order in ONE commit per sub-unit; write the dependency table in the commit body |
| G.25 | a move's residues: read every tool that writes into the moved tree; one home per path constant | `git grep -n -E 'cone_store|generate_' -- tools` for any moved module; each moved symbol has one home that every consumer imports |

### 1.5 The shared retirement gate (created by unit `roles`, extended by every dead sub-unit)

**One engine, no wave in its name.** `tests/gates/test_retirements.py` is the only retirement census of the repository. It replaces the name walks of `tests/gates/test_wave2_retirement.py`: unit `roles` ports that file's `RETIRING_NAMES` into `RetiredName` rows and its `_HOMONYMS` survivors into `KEPT` entries of `tests/gates/_retired/rows_wave2.py`, then deletes that file's own walks (`retired_uses`, `_findings`, `test_nothing_under_orpheus_names_the_retired_vocabulary`, `test_no_test_example_or_tool_names_the_retired_vocabulary`, `test_every_excluded_homonym_survives_where_it_is_said_to`, the planted-census tests and the constants only they read). The tests of that file that are not name walks (the three retiring paths do not exist, the capture scripts, the cold-import list, the SN method space) stay, because they are not retirements by name; the builder reads the file whole before editing it. `tests/gates/test_w0_retirements.py` and `tests/gates/_w0_retired/` are NOT created. Rows live in modules under `tests/gates/_retired/` named per campaign and unit: `rows_wave2.py` (ported), `rows_w0_roles.py`, `rows_w0_dead_a.py`, `rows_w0_dead_b.py`, `rows_w0_dead_c.py`, `rows_w0_dead_d.py`; a later wave adds `rows_w<k>_<unit>.py`. The engine globs `rows_*.py`, so no two units edit one data file. `tests/gates/_retired/__init__.py` holds the row types and `tests/gates/_retired/_finder.py` the definition finder.

- **Row types: three NamedTuples; the gate dispatches on the type.** `Retired(file, qualname)`: the definition is gone from `file`. `Moved(file, qualname, home)`: gone from `file` and defined at `home`, a repository path (a new module named in 3.1, or the single consuming test file for an inline helper). `RetiredName(name)`: a name-level entry (such as `block_role`, written in about 25 files) checked by bare name. No verdict string and no empty home can be spelled.
- **The binding is by definition site and import resolution, never by bare spelling.** A `Retired` or `Moved` row is red when (a) `(file, qualname)` is still defined (the one finder of `_finder.py`: a `def`, a `class`, a module- or class-level `Assign` or `AnnAssign`, nested qualnames, and `self.x = ...` in a method body counts as defining `x`), or when (b) a module resolves an import or an attribute to that definition: an `ImportFrom` whose module (relative imports resolved against the importing file's package) is the row's module and which imports the row's last name; an `Attribute` whose value is a name bound to the row's module by `import m` or `from p import m`; an `__all__` member of the row's module or of a package `__init__` that re-exports it; a string constant that is the first or second argument of `getattr`, `hasattr`, `setattr`, `delattr` or `monkeypatch.setattr` whose first argument resolves to the row's module or class. A homonym, another object with the same name in `derivations/`, `tools/` or `tests/`, cannot redden the row and needs no exemption list. The gate scans every `.py` under `orpheus/` (including `derivations/`), `tests/` (except `tests/gates/_retired/`), `tools/`, `examples/` and `student_resources/`. For a method or attribute row (a dotted qualname) an attribute use through an instance cannot be resolved by name: those rows are red by (a) and by the string and import legs only, and the builder's red loop plus `.venv/bin/nexus dead-references` carry the instance uses (stated, not hidden).
- **A `RetiredName` row is the only bare-name check** and is used only where the name has no surviving homonym. The row's comment carries the measured count: `[M]` `scratch/ontology/w0/spec/amend/roles_count.py` (AST: `Name`, `Attribute`, `alias`, `keyword`, `arg`, `def`/`class` names and string constants over `orpheus`, `tests`, `tools`, `examples`, `student_resources`), over the 12 role names 392 hits in 36 distinct files (production 17 files, tests 19), 11 definitions, all of them the retired family; 0 homonyms. A bare-name walk checks every node kind listed in the previous sentence (`def`/`class` names, `arg` and `Name` in `Store` context included, so a leftover `def block_role` or a `block_role` parameter is seen). A `RetiredName` whose name is the last component of any `KEPT` qualname is a gate error.
- **String leg.** Every string constant that is valid Python source (`ast.parse` succeeds) is parsed and walked with the same finder, so an import written inside a fixture string is read like an import (a relative fixture `from ...numerics.eigenvalue import X` is matched by its module suffix). First red `[M]`: `tests/gates/test_layer_imports.py` at `b68f8b8c` holds six such strings, three that import `BlockRole` and `stamp_boundary_role` (lines 644, 667, 668; unit `roles`) and three that import `dominant_eigenpair` (lines 845, 847, 861; dead-B); no other gate reddens when the imported name is replaced by a name that does not exist, because `_check_source` only parses the string (`m_roles2.py`).
- **Positive controls, inside the gate body.** The finder runs on a synthetic source with one definition of each of five shapes (def, class, class attribute, module `Assign`, `AnnAssign`) and must find all five; the reference walk runs on a synthetic source holding `getattr(x, "NAME")`, `__all__ = ["NAME"]`, `monkeypatch.setattr(m, "NAME", f)`, `from m import NAME` and a fixture string `"from m import NAME\n"`, and must report all five; a negative control (`NAME = 1; f(NAME=True)` in a file that imports nothing from the row's module, and `np.polynomial.laguerre`) must report none.
- **First red** today, per row type: every `Retired` and `Moved` row is defined today `[M]` (`m_dead1.py`, over the file's verdicts before rows 134 to 137 were kept: 105 of 105 retire/move rows defined, 41 of 41 keep rows defined, 0 of 146 rows in `thermal_hydraulics`, `kinetics`, `fuel`); after the deferrals of 3.1 (D2) and the four rows kept by 3.0 item 2 the dead unit's rows are 92, and `roles` adds its own; the ONE further `Retired` entry that is not a row, `QuadratureSpec.structural_flags` (3.0 item 3), is defined today too. **Mutation** (G6, the root-path mechanism): the walk takes a root path; copy the tree, or the one file, to a temporary directory, add `def direct_eigenvalue(...)` to the copy of `numerics/eigenvalue.py` there, call the walk on that root and require red for that row; plant `monkeypatch.setattr(m, "direct_eigenvalue", f)` with `m` bound by `from orpheus.numerics import eigenvalue as m` in a copy of a test file and require red.
- **Rows modules are executable data (G8).** Every `Moved` row carries its concrete `home` path; every kept symbol appears in the module's `KEPT: tuple[tuple[str, str, str], ...]` of `(file, qualname, kept_for)`, assigned to the sub-unit that edits that file (3.3 gives the assignment of the 45 keep rows, among them rows 134 to 137 in dead-A; the four rows deferred by 3.1 D2 are in `KEPT` of their sub-unit with the reason `deferred, private state`); the engine asserts each `KEPT` entry is STILL defined (D1': a builder cannot retire a kept symbol by accident). `KEPT` of `rows_wave2.py` holds the wave-2 survivors that its `_HOMONYMS` listed. No placeholder string appears in a rows module. First red of D1': none (green today by construction); **mutation**: delete the definition of one KEPT row (`DiscreteMeasure.symmetry_groups`, row 58) in a temporary copy of `numerics/measure.py`, call the walk on that root, require red.
- **Exclusions stated in the claim.** The gate does not scan `.claude/`, `docs/` or comments (a text gate cannot sort tense; the tense tables and `dead_references` cover them); docstrings are covered by `tests/gates/test_docstring_xrefs.py`, which resolves roles by importing the target and reddens by design on a role to a retired name (about 34 of 105 rows are named by a docstring role, about 107 occurrences in 51 files `[M]` `m_dead9.py`, approximate: matched by name, not import resolution).
- **One definition finder for every leg**, and for the error-catalogue `Retired:` leg of 3.1: `tests/gates/_retired/_finder.py` (the catalogue gate imports it; `**Retired:**` names `file:qualname` pairs).

### 1.6 Unmeasured, with denominators (carried, not hidden)

- The full-suite red set of any unit was not run (the unit bound is a gate file, not 941). The only suites run: `tests/gates/test_layer_imports.py` 533 passed, `tests/gates/cp/test_cp_boundary_carry.py` 8 passed, plus in-process simulations.
- The defect spy covered 28 of 941 test files; "0 reddened" means 0 of those 28.
- The 73 suspect (file, method-name) pairs in 66 test files (method rows only) are unresolved homonyms until a builder's first red; with them counted 11 test files are shared between sub-units.
- Not read by this specification, listed in 8.3: the 73 suspect pairs and the tense tables of the appendices (the builder reads them); the closure members that `reviewed.json` does not rule; the unread census rows. The full-suite red set of any unit was not run.
- The 377 role lines of the tense table were sorted by a mechanical rule (a role marker within two lines of the name), not read one by one; the other 444 were read. Comments and non-docstring strings were not searched in the tense table.
- A hollow `CPMesh` cannot be built at the base hash (it is refused with `BoundaryRefused`/`UnboundTagError`, #513; one geometry, one tag measured, 2.2); the chord on a hollow body is measured at the function level only (breakpoints `[0.1, 0.3, 0.5]`, 7 values of y, cylinder and sphere), and a hollow body with an untagged inner face was not tried.
- The Ki3 copy was checked deterministic in-process on one platform only.
- Sub-unit `angular-mass` (5.7): the memo's count of 59 tests that read a divisor that differs, its 356-file swap run and its GL12 solve table were NOT re-run (the in-process fix of 5.7.3 was run over 21 test files, 343 tests, and 5 of those files again with an equality witness; 5.7.7 gives both populations); the `Fraction` exact-sum leg was not re-run; the full-suite red set of the change was not run; the half-trace mass at `orpheus/sn/boundary/realized_operators.py` line 522 (`np.sum(in_weights)`) is a sixth `np.sum` spelling of a mass, of another measure, found by the census and left (8.1).
- Unit `quadrature-ties` (2.0.12): the memo's 4643-pair before/after, its random census, its independent oracle, its gcd cost table and the per-file counts of the numerics gate files were NOT re-run (this pass re-ran the 12-row gate, its four mutations and the population probe, and ran the generic builders, the decision code of 2.0.5 and the rows of G4 and G5 in a prototype); the four call-site edits to `HighPrecision(...)` were not run (the prototype accepted an integer digit count at the unchanged sites); the full-suite red set of the unit was not run; the pyright count of the three edited production files was not measured.

## 2. Unit `insulation`: production CP stops importing the reference (#641, item 1)

**Outcome.** `orpheus/cp/solver.py` imports nothing from `orpheus.derivations`; the three symbols it used have production homes; a gate in its own table forbids the edge from every production module and from a cold import; the three carry rows that move by ulps are re-baselined with the reason. `derivations/` is not edited (decision 1): it keeps its copies of the chord, the rule and the Bickley interpolant, and the production copy of the interpolant is a TWIN of the reference's by construction (X4), so #641 item 2 (an independent verification chain for the kernel) stays OPEN: it is closed by hand only when its item 2 is a filed issue and production imports of `orpheus.derivations` are 0 (2.6.6). The unit integrates AFTER the precondition unit `quadrature-ties` (2.0).

### 2.0 Precondition unit `quadrature-ties`: an exact value on a float64 rounding tie is rounded half to even, never refused (decision I1; the orchestrator's rulings of 2026-11-20)

**Outcome.** `LEGENDRE.on(a, b).gauss(n)` (and the pushforward of every family on `[-1, 1]`) builds for every panel and every `n` where an exact node or weight is a rational number sitting on a float64 rounding tie, and each such value is the float64 nearest to it with the tie rounded to the EVEN neighbour (IEEE 754 `roundTiesToEven`, which is what "correctly rounded to nearest" means). The six closed-form recurrence builders become generic over the number type, so one definition is evaluated in high-precision floats (as today) or in rationals (new), and the affine remap maps the rational recurrence exactly. The unit integrates FIRST (7.1 step 0), before `insulation`, because the per-panel spelling of `insulation` calls this code on panels that raise today. Key `quadrature-ties` (the placeholder's earlier name, `quadrature-odd-n`, is retired: parity was never the law, 2.0.1). A builder's brief for this unit is: section 1 whole, this section 2.0 whole (self-contained: it needs no other part of section 2), section 7, and the items 37 to 41 of 8.0, the rulings 8.42 to 8.46 of 8.4 and 8.5.

#### 2.0.1 The premise, measured (every fact a builder needs, each with its command)

Instrument for every row of this list unless another is named: a probe copied to `scratch/ontology/w0/spec/amend/ties/` from `scratch/ontology/w0/oddn/` (the memo of the investigation) and re-run by this pass at `b68f8b8c` under `.venv/bin/python -O`; a builder has no `scratch/` and re-runs the gate of 2.0.7 instead.

- **Where it raises.** `orpheus/numerics/generating_measure.py`, `_certified_float` (lines 594 to 609 at the base hash), called by `_rule_at` at line 585 (nodes) and line 591 (weights), reached by `GeneratingMeasure.gauss` through `_correctly_rounded_rule`. It certifies a value only when both ends of its error interval `[v - rho|v|, v + rho|v|]`, `rho = 10^(10 - digits)`, round to one float64; otherwise it raises `ArithmeticError` (message tail `lies within its error radius of a float64 rounding boundary, so it cannot be certified correctly rounded`).
- **Why an exact value cannot pass.** The affine pushforward `on(a, b)` maps the reference rule to `x = c + s t`, `w = s omega`, `c = (a + b)/2`, `s = (b - a)/2` (exact dyadic rationals). Some values of the remapped rule are then EXACT rationals: the centre node `c` (odd `n`) and `s omega` for every rational reference weight `omega` (all weights at `n = 1, 2, 3`, the centre weight at odd `n`). An exact value with 54 significant bits is the midpoint of two adjacent float64: a TIE. Its error interval contains the tie, so the interval ends round to different floats, whatever the precision. The reference rule on `[-1, 1]` never meets the case: its only rational node is 0 (radius `|0| rho = 0`) and its rational weights (2, 1, 5/9, 8/9) are either representable or not dyadic.
- **Parity is not the law.** `[M]` the population probe `p1_population.py` (34 panels in 6 families, `n` in 1 to 64, 2176 (panel, n) pairs): 384 of 1088 odd-`n` pairs raise and 0 of 1088 even-`n` pairs; the predicate "raises iff `n` is odd and `(a + b)/2` is a tie (an exact `Fraction` test)" has 0 mismatches over the 2176 pairs; the tied panels are 6 of 12 on the carried cylinder and 4 of 12 on the carried sphere; the smallest reproducer is `n = 1` on `(0.7937005259840998, 0.8735804647362988)`. The even-`n` class exists: `[M]` (the memo's `p3_classes.log`, not re-run) the weight `(b - a)/2` at `n = 2` raises on the panel `(0.2615032494292951, 0.7646335475579331)`, whose width is a tie while its midpoint is not.
- **It is not the families' either.** `[M]` (`probes_ab.py`, run on today's tree) `CHEBYSHEV_T.on(0.2, 0.4).gauss(n)` raises for `n` = 1 and 3, `CHEBYSHEV_U.on(0.2, 0.4).gauss(3)`, `jacobi(0.5, 0.5).on(0.2, 0.4).gauss(3)`, `jacobi(-0.5, -0.5).on(0.2, 0.4).gauss(5)` raise, and so does `jacobi(0.0, 0.0).on(0.2615032494292951, 0.7646335475579331).gauss(1)`: 6 of 6 non-Legendre rows raise (and the 2 Legendre rows of the same probe, `(0.2, 0.4)` at `n` = 1 and 3); a seventh, `CHEBYSHEV_T.on(0.2615032494292951, 0.7646335475579331).gauss(2)`, builds (its weights carry the irrational mass `pi`, so no weight is rational and its nodes are not tied).
- **The contract's silence.** `docs/theory/foundations/discrete_measures.rst`, "The Gauss rules are correctly rounded" (the contract paragraph starting `**The contract.**`, line 836 at the base hash) says each value is "the float64 **nearest to its exact value**" and that the contract covers "the affine remap `on`"; step 5 (line 936) refuses a value so close to a midpoint that the approximation "could land on either side". The page states no tie rule (for an exact tie "nearest" names two floats). So the refusal of an exact tie is not a declared scope edge; it is a certificate correct for APPROXIMATE values (a Newton root, a Christoffel number known to about 1e-35) applied to an EXACT value, where there is no guess to make. The decision is ruled: round half to even.

#### 2.0.2 File set (by tree)

Hard (edited):
- production: `orpheus/numerics/generating_measure.py`, `orpheus/numerics/recurrence.py`.
- tests: `tests/gates/numerics/test_generating_measure.py` (the gate promoted into it, and the one existing test edited, 2.0.8), `tests/gates/numerics/test_concept_order.py` (a tier entry for the new module).
- docs: `docs/theory/foundations/discrete_measures.rst` (the contract paragraph and the sentences of 2.0.9), `docs/theory/verification/error_catalog.rst` (the entry stub of 2.0.10).

New: `orpheus/numerics/rational_polynomial.py` (2.0.5 E1), `tests/gates/numerics/test_rational_polynomial.py` (2.0.7 G7).

Run, not edited: `tests/gates/test_error_catalogue_reconciles.py` and `tests/gates/test_vv_harness_audit.py` (the marker's target must exist, 2.0.10), `tests/gates/test_layer_imports.py`, `tests/gates/numerics/test_gauss_rule_fingerprint.py`, `test_gauss_rules_correctly_rounded.py`, `test_exactness.py`, `test_registry.py`, `test_measure.py`, `test_legendre_basis.py`.

`[M]` the call-site set of the changed interface (an AST pass, `scratch/ontology/w0/spec/amend/ties/census.py`: `Name`, `Attribute`, `alias` and `keyword` nodes named `legendre_recurrence`, `chebyshev_t_recurrence`, `chebyshev_u_recurrence`, `hermite_recurrence`, `jacobi_recurrence`, `laguerre_recurrence`, `exact_recurrence`, `RecurrenceCoefficients` or `_at_digits` over every `*.py` under `orpheus/` (including `orpheus/derivations/`), `tests/`, `tools/`, `examples/`; positive control: the definition of `legendre_recurrence` in `recurrence.py` is found): 3 files of the 1391 `*.py` files under those four trees hold any of the names, `orpheus/numerics/recurrence.py` (the definitions), `orpheus/numerics/generating_measure.py` (the six builders imported at lines 203 to 211 and used at 655, 664, 671, 678, 719, 768; `exact_recurrence` at lines 300, 311, 342, 422, 477, 494, 558, 655, 664, 671, 678, 719, 768) and `tests/gates/numerics/test_generating_measure.py` (line 293, an identity assertion `general.exact_recurrence is not known.exact_recurrence`, which stays valid: `jacobi(a, b)` and the named constants stay different callables). 0 occurrences in `orpheus/derivations/`, `tools/`, `examples/` and any other test. `git grep -n -E "exact_recurrence|RecurrenceCoefficients|_recurrence|_at_digits" -- docs` finds two lines, `docs/theory/foundations/discrete_measures.rst` lines 896 and 898 (the construction's step 1; 2.0.9). Prose plan files (`.claude/plans/numerics_untangling.md`) name the symbols as history and are not edited.

Shared with other units (computed from the file sets of this specification by `scratch/ontology/w0/spec/amend/ties/overlap.py`, exact path matches read by section; 7.2 has the rows): `generating_measure.py` with dead-A (D-A.2, rows 54 and 55) and `defects` (5.1, #711); `test_generating_measure.py` with dead-A (D-A.3, rows 54 and 55 consumers); `discrete_measures.rst` with dead-A (D-A.4, D-A.6) and by prose with `defects` (5.1); `error_catalog.rst` with dead-A (ERR-081, ERR-118), `roles` (ERR-075) and by read-for-tense with `insulation` (2.1); `test_concept_order.py` with `insulation` (2.3, the tier entry of `bickley`). `recurrence.py` and the two new files are in no other unit.

#### 2.0.3 The law, whole, with its converse and its domain

Notation. For a measure `m` and `n >= 1`, `x_1 < ... < x_n` are the exact nodes (the roots of the monic orthogonal polynomial `pi_n`) and `w_i` the exact weights. For the exact value `v` of a node or weight and the working precision `d`, `rho = 10^(10 - d)` (`_CERTIFIED_DIGITS = 10`); `[low, high] = [float(v - rho|v|), float(v + rho|v|)]` (the existing computation, `float` of an `mpf` rounds to nearest-even). The boundaries of float64 rounding are the midpoints of adjacent floats.

**The law.** `m.gauss(n)` returns, for every node and every weight, the float64 nearest to the exact value, a tie going to the even neighbour (`roundTiesToEven`); or it raises `ArithmeticError`. Precisely, at the fine working precision of `_correctly_rounded_rule` (60 digits; the coarse run is certified the same way, and the two must agree, as today):
1. `low == high`: the value is `low` (unchanged: the new code is not consulted).
2. `low != high` and `low`, `high` are ADJACENT floats (`nextafter(low, +inf) == high`): the interval holds exactly ONE boundary `T = (low + high)/2`, a rational number, computed exactly as `(Fraction(low) + Fraction(high)) / 2`. Then, if the decision `v == T` is available and says EQUAL, the value is `float(T)` (`Fraction.__float__` is the int true division, which is correctly rounded, half to even); if it says UNEQUAL, or is unavailable, the value is refused exactly as today (same exception type, same message).
3. `low != high` and the ends are not adjacent: the interval holds two or more boundaries; refused exactly as today.

**The two decisions** (Q is the rational field; both read the measure's rational recurrence `(alpha_k, beta_k)`, the monic form `pi_{k+1}(x) = (x - alpha_k) pi_k(x) - beta_k pi_{k-1}(x)`, `pi_{-1} = 0`, `pi_0 = 1`, with the storage convention `beta_0 = mu_0`; Gautschi 2004, Eq. 1.3.13):
- **Node** `x_i`: `x_i = T` iff `pi_n(T) = 0`, evaluated exactly in Q. It reads `alpha_0 .. alpha_{n-1}` and `beta_1 .. beta_{n-1}`; it does NOT read `beta_0`. (Why it is an iff: `T` is within `rho|x_i|` of `x_i`, and the nodes of every rule built here are separated by far more than `2 rho |x|` `[R]` (the smallest gap of the 64-point Legendre rule is about 1e-3 relative; `rho` at 60 digits is 1e-50); so `T` is a root of `pi_n` iff it is the root `x_i`.)
- **Weight** `w_i`: `w_i = T` iff `x_i` is a root of `g = gcd(pi_n, K_n - 1/T)` in `Q[x]`, with `K_n(x) = sum_{k<n} pi_k(x)^2 / h_k`, `h_k = beta_0 beta_1 ... beta_k` (Christoffel numbers: `w_i = 1/K_n(x_i)`, Gautschi 2004 Thm 1.46 with Eq. 1.3.13). Proof: `w_i = T` iff `K_n(x_i) - 1/T = 0` iff `x_i` is a common root of `pi_n` and `K_n - 1/T` iff it is a root of their gcd (the roots of a gcd over Q are the common roots over C). `g` is constant: UNEQUAL (no node has weight `T`). Otherwise `x_i` is isolated by a sign change of `pi_n` on the rational bracket `x_i +- 10^-20 max(1, |x_i|)` (`x_i` is the working-precision `mpf` the rule holds, converted EXACTLY to a `Fraction`); since `g` divides `pi_n` and `pi_n` has simple roots, `g` changes sign (or vanishes) at a bracket end iff `g(x_i) = 0`. It reads `beta_0`.

**Domain, as part of the law (not a guard).** The decision exists for a measure iff its recurrence has a rational image: `alpha_k` and `beta_k` (`k >= 1`) are rationals for the six shipped families at every admitted parameter (the table of 2.0.4); the node decision exists for all of them; the weight decision exists iff `beta_0` is rational. Where the decision does not exist (`is_exactly is None`) the law's arm 2 is the refusal of today: the domain of the exact rounding is the set of measures with a rational recurrence, and outside it the construction is as before. A weight whose mass is `pi`, `sqrt(pi)` or a multiple of them (Chebyshev T and U, Hermite) is `beta_0` times an algebraic number, hence transcendental `[R]` (a nonzero algebraic multiple of `pi` is transcendental, Lindemann), so no such weight is ever a float: for those measures the unavailable weight decision refuses nothing that could be a tie.

**The converse, stated.** The law does NOT say that an UNEQUAL single-boundary value is refused forever: a higher working precision would round it (an escalation up to a cap that still raises is the specified improvement, out of scope: the unit refuses "exactly as today"). It does NOT claim the decision is cheap for every `n` (off a straddle it costs one rational recurrence of length `n + 1` per working precision, `[M]` unmeasurable in the cold rule times of 2.0.5 E3 item 3): `[M]` (the memo's `p9_gcd_cost.py`, not re-run) the worst case, a straddle with a trivial gcd at `n = 63`, takes 5.8 s; the decision runs once per working precision and only for a value whose interval holds exactly one boundary, which `[R]` never happens for an irrational node at `n >= 4`. It does NOT make `LEGENDRE.on(a, b)` symmetric (`is_symmetric` stays `False` on a panel; the mirror path is the reference rule's). The law has no parity clause and no panel clause.

**Raw inputs.** `v`, the `mpf` the rule computed; `digits`, `measure`, `n` (as today); `is_exactly`, `None` or a function `Fraction -> bool` built from `(measure, n)` and, for a weight, the node `x_i` (an `mpf`). The rational recurrence is NOT an input held beside the measure: it is `measure.exact_recurrence(n + 1, RATIONAL)`, the one definition evaluated in rationals (2.0.4).

#### 2.0.4 The generic form of the closed forms (the ruling: one definition, two number types, no new stored field)

`[M]` as spelled today (`orpheus/numerics/recurrence.py`, 160 lines): `RecurrenceCoefficients = Callable[[int, int], "tuple[list, list]"]` (line 60), `(n, digits) -> (alpha, beta)` with `n` `mpmath.mpf`; the private lift `_at_digits(closed_form)` (the decorator at line 63: a wrapper that enters `mpmath.workdps(digits)` and calls the closed form of `n`); four plain builders `legendre_recurrence`, `chebyshev_t_recurrence`, `chebyshev_u_recurrence`, `hermite_recurrence` (each `n -> (alpha, beta)`, decorated), and two builder factories `jacobi_recurrence(a, b)` and `laguerre_recurrence(a)` returning a `RecurrenceCoefficients` (`a`, `b` are float64, "enter the working precision exactly"). `GeneratingMeasure.exact_recurrence` (field, line 300, `compare=False`) carries one of them; the callers pass `digits` as an `int` (4 sites in `generating_measure.py`: line 311 `_VIEW_DIGITS`, line 342 `_VIEW_DIGITS`, line 480 inside the remap, line 558 `digits`). The census of 2.0.2 shows no other caller.

**The specified generic form.** The second argument of a builder changes from the digit count `digits: int` to an ARITHMETIC: an object that says which number type the closed form is evaluated in.

```python
# orpheus/numerics/recurrence.py  (imports: contextlib, dataclasses, fractions, math, typing, mpmath; nothing from orpheus)
Number = mpmath.mpf | Fraction            # a coefficient of the recurrence

class Irrational(ArithmeticError):
    """The rational arithmetic was asked for a quantity with no rational value (a transcendental mass)."""

class Arithmetic(Protocol):
    def working(self) -> ContextManager: ...                  # mpmath.workdps(d) | contextlib.nullcontext()
    def lift(self, x: int | float | Fraction) -> Number: ...  # the EXACT image of an int, a float64 or a Fraction
    def pi(self) -> Number: ...
    def sqrt(self, x: Number) -> Number: ...
    def gamma(self, x: Number) -> Number: ...
    def beta(self, a: Number, b: Number) -> Number: ...
    def power(self, base: int, exponent: Number) -> Number: ...

@dataclass(frozen=True)
class HighPrecision:                                          # today's arithmetic, `digits` decimal digits
    digits: int
    def working(self): return mpmath.workdps(self.digits)
    def lift(self, x): return mpmath.mpf(x) if not isinstance(x, Fraction) else mpmath.mpf(x.numerator) / x.denominator
    def pi(self): return +mpmath.pi                           # an mpf at the working precision, not the lazy constant
    def sqrt(self, x): return mpmath.sqrt(x)
    def gamma(self, x): return mpmath.gamma(x)
    def beta(self, a, b): return mpmath.beta(a, b)
    def power(self, base, e): return mpmath.power(base, e)

class _Rational:                                              # exact: Fraction; a transcendental quantity raises Irrational
    def working(self): return contextlib.nullcontext()
    def lift(self, x): return Fraction(x)                     # a float64 is a binary rational: exact
    def pi(self): raise Irrational("pi")
    def sqrt(self, x):                                        # rational only for a perfect-square ratio
        r = math.isqrt(x.numerator), math.isqrt(x.denominator)
        if r[0] ** 2 == x.numerator and r[1] ** 2 == x.denominator: return Fraction(*r)
        raise Irrational(f"sqrt({x})")
    def _positive_integer(self, x):
        x = Fraction(x)
        if x.denominator != 1 or x < 1: raise Irrational(f"gamma or beta of {x}")
        return int(x)
    def gamma(self, x): return Fraction(math.factorial(self._positive_integer(x) - 1))
    def beta(self, a, b):
        a, b = self._positive_integer(a), self._positive_integer(b)
        return Fraction(math.factorial(a - 1) * math.factorial(b - 1), math.factorial(a + b - 1))
    def power(self, base, e):
        e = Fraction(e)
        if e.denominator != 1: raise Irrational(f"{base}**{e}")
        return Fraction(base) ** int(e)

RATIONAL = _Rational()

RecurrenceCoefficients = Callable[[int, Arithmetic], "tuple[list[Number], list[Number | None]]"]
```

The wrapper `_at_digits` is replaced by `over_arithmetic`, public because `generating_measure.py` uses it for the remap (the same shape: enter `arithmetic.working()`, call `closed_form(n, arithmetic)`); a closed form is a function `(n, A) -> (alpha, beta)` of the arithmetic `A`. The mass is the ONLY place a transcendental function appears, and it is written through one helper:

```python
def _mass(thunk) -> Number | None:
    """beta_0: the value, or None when the arithmetic has no value for it (an irrational mass)."""
    try: return thunk()
    except Irrational: return None
```

The six builders, each written ONCE (this is the reference body; the law of 2.0.3 and the table below are the contract; a number leaves `lift` only through an exact image):

```python
@over_arithmetic
def legendre_recurrence(n, A):
    return [A.lift(0)] * n, ([A.lift(2)] + [A.lift(k * k) / (4 * k * k - 1) for k in range(1, n)])[:n]

@over_arithmetic
def chebyshev_t_recurrence(n, A):
    beta = [_mass(A.pi), A.lift(1) / 2] + [A.lift(1) / 4] * max(n - 2, 0)
    return [A.lift(0)] * n, beta[:n]

@over_arithmetic
def chebyshev_u_recurrence(n, A):
    return [A.lift(0)] * n, [_mass(lambda: A.pi() / 2)] + [A.lift(1) / 4] * (n - 1)

@over_arithmetic
def hermite_recurrence(n, A):
    return [A.lift(0)] * n, [_mass(lambda: A.sqrt(A.pi()))] + [A.lift(k) / 2 for k in range(1, n)]

def jacobi_recurrence(a: float, b: float) -> RecurrenceCoefficients:
    @over_arithmetic
    def recurrence(n, A):
        Aa, Bb = A.lift(a), A.lift(b)
        AB = Aa + Bb
        alpha = [(Bb - Aa) / (AB + 2)]
        beta = [_mass(lambda: A.power(2, AB + 1) * A.beta(Aa + 1, Bb + 1))]
        if n > 1:   # k = 1 is the removable singularity of the general formula (the comment of today's body is kept)
            alpha.append((Bb * Bb - Aa * Aa) / ((2 + AB) * (4 + AB)))
            beta.append(4 * (1 + Aa) * (1 + Bb) / ((2 + AB) ** 2 * (3 + AB)))
        for k in range(2, n):
            two_k_ab = 2 * k + AB
            alpha.append((Bb * Bb - Aa * Aa) / (two_k_ab * (two_k_ab + 2)))
            beta.append(4 * k * (k + Aa) * (k + Bb) * (k + AB) / (two_k_ab**2 * (two_k_ab + 1) * (two_k_ab - 1)))
        return alpha[:n], beta[:n]
    return recurrence

def laguerre_recurrence(a: float) -> RecurrenceCoefficients:
    @over_arithmetic
    def recurrence(n, A):
        Aa = A.lift(a)
        return [2 * k + Aa + 1 for k in range(n)], [_mass(lambda: A.gamma(Aa + 1))] + [k * (k + Aa) for k in range(1, n)]
    return recurrence
```

`[M]` (`compare_generic.py`, this body in `generic_recurrence.py`, over `legendre`, `chebyshev_t`, `chebyshev_u`, `hermite`, six `jacobi(a, b)` pairs `(0,0) (-0.5,-0.5) (1.5,2.25) (0.5,0.5) (2,3) (-0.5,0.5)` and four `laguerre(a)`, `a` in `0, 1.5, 2, -0.5`, `n` in 1, 2, 3, 8, 64, digits 40 and 60): 4368 of 4368 `mpf` entries equal today's builders by value (`mpmath.mpf(old) == new` under `workdps(digits)`); the one structural difference is `chebyshev_t_recurrence`'s `beta[0]`, the lazy `mpmath.pi` constant today and an `mpf` at the working precision after (its value at 40 and 60 digits is equal), and the same for the Chebyshev U and Hermite masses. Timing `[M]` (`time_builders.py`, min of 20, `n = 129`, 60 digits): Legendre 0.20 ms before and 0.22 ms after, `jacobi(1.5, 2.25)` 1.38 ms before and 1.38 ms after.

**What each family admits in rationals** (`[M]` `compare_generic.py`, the arithmetic `RATIONAL`, `n = 6`; "yes" is an all-entries test `isinstance(x, Fraction)`):

| family | `alpha_k` | `beta_k`, `k >= 1` | `beta_0` in `RATIONAL` | the node decision | the weight decision |
|---|---|---|---|---|---|
| Legendre | yes | yes | `2` | yes | yes |
| Chebyshev T, U | yes | yes | `None` (`pi`, `pi/2`) | yes | none (mass irrational) |
| Hermite | yes | yes | `None` (`sqrt(pi)`) | yes | none |
| Jacobi `(a, b)`, `a`, `b` non-negative integers | yes | yes | rational (`(0,0)`: 2; `(2,3)`: 16/15) | yes | yes |
| Jacobi, any other `(a, b)` (`(-0.5,-0.5)`, `(0.5,0.5)`, `(1.5,2.25)`, `(-0.5,0.5)`) | yes | yes | `None` | yes | none |
| Laguerre `a` a non-negative integer | yes | yes | rational (`a = 0`: 1; `a = 2`: 2) | yes | yes |
| Laguerre, any other `a` (`1.5`, `-0.5`) | yes | yes | `None` | yes | none |

So the decision code reads `beta[0] is None` as "the weight decision is unavailable"; the node decision never reads `beta[0]`. A family with non-rational parameters does not exist among the six (a float64 parameter is a binary rational, so `alpha_k`, `beta_k` are rational for every admitted `a`, `b`): what is non-rational is the MASS. A family whose `RATIONAL` evaluation raises `Irrational` outside the mass (a builder written later with a transcendental coefficient elsewhere) has no decision: `_rational_recurrence` catches `Irrational` and returns `None`, so both decisions are `None` (the law's domain, 2.0.3), so the value is refused as today.

**The remap, generic** (`_affine_remap`, `generating_measure.py` lines 465 to 495): `remapped(n, digits)` becomes `remapped(n, arithmetic)`, decorated by `over_arithmetic` (imported from `.recurrence`):

```python
def remapped(n, A):
    alpha, beta = inner(n, A)
    scale, shift = (A.lift(b) - A.lift(a)) / 2, (A.lift(a) + A.lift(b)) / 2        # a, b are the float64 endpoints, exact
    return ([scale * x + shift for x in alpha],
            [None if beta[0] is None else scale * beta[0]] + [scale**2 * x for x in beta[1:]])
```

so `alpha'_k = s alpha_k + c`, `beta'_0 = s beta_0`, `beta'_k = s^2 beta_k` (k >= 1) hold in `Fraction` (with `s`, `c` exact), and in `HighPrecision` the values are today's (the same operations: `scale` and `shift` are formed at the working precision from the same exact endpoints; `[M]` 4368 of 4368 above include no remap, the remap's equality is the gate rows of 2.0.7 and the fingerprint gate).

**Call sites that change** (the four of `generating_measure.py` named above, the whole set by the census of 2.0.2): line 311 `self.exact_recurrence(n, HighPrecision(_VIEW_DIGITS))`, line 342 likewise, line 558 `measure.exact_recurrence(n + 1, HighPrecision(digits))`, and the remap's closure. `HighPrecision` is imported from `.recurrence`. The identity assertion of `test_generating_measure.py` line 293 needs no edit. The cache key of `_correctly_rounded_rule(measure, recurrence, n)` is unchanged (`recurrence` is the callable `measure.exact_recurrence`, and a callable is a valid key); the key argument of `lru_cache` for the rational evaluation is not needed (the decision is evaluated only on a single-boundary straddle).

**Why this form (the structural reasons; the rejected candidates).** One definition per quantity (X4): the closed forms stay six, each evaluated in two number types. Rejected: (1) a second set of `Fraction` builders (the memo's prototype did this for Legendre: a second closed form per family, so the exact value and the float value would agree tautologically and could drift apart, which is the cardinal-2 stop signal); (2) a stored field `rational_recurrence` beside `exact_recurrence` (ruled out: a derived value stored beside its definition, and the identity of `compare=False` fields is the #711 problem again); (3) a type parameter that is the number CLASS (`mpf` or `Fraction`): `pi`, `gamma` and `beta` cannot be spelled from a class, they need the arithmetic object; (4) deciding in `mpf` at higher precision: an exact tie is not decidable by any finite precision (the whole defect).

#### 2.0.5 The exact edits

**E1. New module `orpheus/numerics/rational_polynomial.py`** (layer `numerics`, L1; phase: a pure function library over a type, no space, field or operator; imports only `fractions`; the module the decision needs because a polynomial over Q and its gcd are a concept in their own right, kept out of `generating_measure.py`; no re-export from `orpheus/numerics/__init__.py`). A polynomial is a tuple of `Fraction` coefficients, lowest degree first. Reference body (`[M]` this body, `rational_polynomial.py`, ran in the 12-row gate below):

```python
Polynomial = tuple[Fraction, ...]
ZERO: Polynomial = (Fraction(0),)
ONE: Polynomial = (Fraction(1),)

def trim(p) -> Polynomial:
    q = list(p)
    while len(q) > 1 and q[-1] == 0: q.pop()
    return tuple(q) if q else ZERO

def add(p, q) -> Polynomial:
    m = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m)])

def scale(p, c) -> Polynomial: return trim([c * x for x in p])

def multiply(p, q) -> Polynomial:
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q): out[i + j] += x * y
    return trim(out)

def remainder(p, q) -> Polynomial:
    """The remainder of p by the nonzero q, by long division over Q."""
    p, q = list(trim(p)), trim(q)
    if q == ZERO: raise ZeroDivisionError("polynomial division by zero")
    while len(p) >= len(q) and p != [0]:
        f = p[-1] / q[-1]
        shift = len(p) - len(q)
        for i, y in enumerate(q): p[shift + i] -= f * y
        p = list(trim(p[:-1])) if len(p) > 1 else [Fraction(0)]
    return trim(p)

def gcd(p, q) -> Polynomial:
    """The monic greatest common divisor in Q[x] (Euclid)."""
    p, q = trim(p), trim(q)
    while q != ZERO: p, q = q, remainder(p, q)
    return scale(p, 1 / p[-1]) if p != ZERO else ZERO

def evaluate(p, x) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p): acc = acc * x + c
    return acc
```

Its tier in `tests/gates/numerics/test_concept_order.py::TIERS`: the primitives tier (`("evidence", "face_layout", "gauge", "index_set", "instance_cache", "metric", "recurrence", "registry", "tangential", "units", "vector")`); add `"rational_polynomial"` immediately BEFORE `"recurrence"` (alphabetical; the line that `insulation` edits is the head of the tuple, `"bickley"`, 7.2). Without the entry `test_every_module_has_a_rank` is red (a by-design red of the unit, closed in the same commit). `generating_measure` (tier `("frame", "generating_measure")`) imports it: downstream of the primitives tier, so the edge follows the order.

**E2. `orpheus/numerics/recurrence.py`**: the Arithmetic, `Irrational`, `HighPrecision`, `RATIONAL`, `_mass`, `over_arithmetic` and the six builders of 2.0.4; `_at_digits` is deleted (0 other users, 2.0.2) and the module docstring's paragraph "Precision. Every coefficient is evaluated at an explicit number of decimal digits (the precision is an argument, never mpmath's ambient setting)..." and the sentence "each function returns `(alpha, beta)`, two lists of `n` `mpmath.mpf`" are rewritten (2.0.9). `__all__` gains `Arithmetic`, `HighPrecision`, `Irrational`, `RATIONAL`, `Number`, `over_arithmetic`.

**E3. `orpheus/numerics/generating_measure.py`**, six hunks, all anchored by quoted text (the symbol decides, line numbers are at `b68f8b8c`):

1. The call sites and the remap of 2.0.4 (lines 311, 342, 480 to 486, 558).
2. Add the decision code beside `_certified_float`, under a new header comment `# The exact tie` (after `_rule_at`, before `_certified_float`). Imports added: `from fractions import Fraction`, `from . import rational_polynomial as _rp` (or `from .rational_polynomial import ...`), and `Irrational`, `RATIONAL`, `HighPrecision` from `.recurrence`. Reference body (`[M]` this body, in the hoisted form of item 3, ran in the gate below, `final_decision2.py`):

```python
_ISOLATION = Fraction(1, 10**20)     # relative half-width of the rational bracket that isolates a node

def _exact(x) -> Fraction:
    """An ``mpf`` as the exact rational it is. NOT ``Fraction(str(x))`` and NOT ``mpf(x)``: re-rounding to the ambient precision
    would move a 136-bit node by 1e-16 and the bracket of 1e-20 would no longer isolate it."""
    mantissa, exponent = x.man_exp
    return Fraction(int(mantissa)) * Fraction(2) ** int(exponent)

def _monic_family(alpha, beta, n):
    """pi_0 .. pi_n as polynomials over Q (pi_{k+1} = (x - alpha_k) pi_k - beta_k pi_{k-1}; beta_0 is not read)."""
    pis, previous = [_rp.ONE], _rp.ZERO
    for k in range(n):
        step = _rp.multiply((-alpha[k], Fraction(1)), pis[-1])
        nxt = _rp.add(step, _rp.scale(previous, -beta[k])) if k > 0 else step
        previous = pis[-1]
        pis.append(nxt)
    return pis

def _pi_n_at(alpha, beta, n, x: Fraction) -> Fraction:
    previous, current = Fraction(0), Fraction(1)
    for k in range(n):
        previous, current = current, (x - alpha[k]) * current - (beta[k] * previous if k > 0 else 0)
    return current

def _rational_recurrence(measure, n):
    """``(alpha, beta)`` of ``measure`` as ``Fraction`` (length ``n + 1``; ``beta[0]`` is ``None`` for an irrational mass), or None
    when the rational arithmetic has no value for a coefficient: the domain of the exact decision (2.0.3)."""
    try:
        return measure.exact_recurrence(n + 1, RATIONAL)
    except Irrational:
        return None

def _node_is_exactly(rational, n):
    """``tie -> bool`` deciding ``x_i == tie`` (``pi_n(tie) == 0``), or None when there is no rational recurrence."""
    if rational is None:
        return None
    alpha, beta = rational
    return lambda tie: _pi_n_at(alpha, beta, n, tie) == 0

def _weight_is_exactly(rational, n, x, name):
    """``tie -> bool`` deciding ``w_i == tie`` for the node ``x`` (an ``mpf``), or None when the recurrence or the mass is not rational."""
    if rational is None or rational[1][0] is None:
        return None
    alpha, beta = rational

    def decide(tie: Fraction) -> bool:
        pis = _monic_family(alpha, beta, n)
        norms, h = [], Fraction(1)
        for k in range(n):
            h *= beta[k]
            norms.append(h)                                      # h_k = beta_0 ... beta_k
        reciprocal_christoffel = _rp.ZERO                        # K_n
        for k in range(n):
            reciprocal_christoffel = _rp.add(reciprocal_christoffel, _rp.scale(_rp.multiply(pis[k], pis[k]), 1 / norms[k]))
        common = _rp.gcd(pis[n], _rp.add(reciprocal_christoffel, (-1 / tie,)))
        if len(common) == 1:
            return False
        centre = _exact(x)
        half = _ISOLATION * max(1, abs(centre))
        low, high = centre - half, centre + half
        if (_rp.evaluate(pis[n], low) > 0) == (_rp.evaluate(pis[n], high) > 0):
            raise ArithmeticError(f"{name}.gauss({n}): the rational bracket does not isolate the node")
        return (_rp.evaluate(common, low) > 0) != (_rp.evaluate(common, high) > 0) or _rp.evaluate(common, low) == 0 or _rp.evaluate(common, high) == 0

    return decide
```

3. `_rule_at`: the rational recurrence is evaluated ONCE per call (not once per value), right after the line `alpha, beta = measure.exact_recurrence(n + 1, HighPrecision(digits))`: `rational = _rational_recurrence(measure, n)`. Line 585 `nodes = tuple(_certified_float(x, digits, measure, n) for x in roots)` becomes `nodes = tuple(_certified_float(x, digits, measure, n, is_exactly=_node_is_exactly(rational, n)) for x in roots)`; line 591 `return nodes, tuple(_certified_float(w, digits, measure, n) for w in weights)` becomes `return nodes, tuple(_certified_float(w, digits, measure, n, is_exactly=_weight_is_exactly(rational, n, x, measure.name)) for x, w in zip(roots, weights))` (the argument is passed by KEYWORD, as in the run that measured it; for a symmetric measure `roots` and `weights` are the mirrored lists of equal length, so `zip` pairs each weight with its node). The decision functions are built per value and the decision itself runs lazily, only on a single-boundary straddle. `[M]` this hoisted form (`final_decision2.py`) passes the 12 rows and the four mutations as before (12 passed; 10, 7, 3, 1 failed), and costs nothing measurable (`time_rule.py`, min of 3 cold runs, today against fixed: `LEGENDRE.gauss(128)` 0.643 s and 0.648 s, `jacobi(-0.5, 0.75).gauss(128)` 1.315 s and 1.316 s, `LEGENDRE.on(0.0, 1.0).gauss(63)` 0.319 s and 0.320 s).
4. `_certified_float(value, digits, measure, n, is_exactly=None)`: after the `low == high` return, the new arm; the message is unchanged. Reference body (`[M]`):

```python
def _certified_float(value, digits: int, measure: "GeneratingMeasure", n: int, is_exactly=None) -> float:
    with mpmath.workdps(digits):
        radius = abs(value) * mpmath.mpf(10) ** (_CERTIFIED_DIGITS - digits)
        low, high = float(value - radius), float(value + radius)
    if low == high:
        return low
    if is_exactly is not None and float(np.nextafter(low, np.inf)) == high:      # one boundary only: the tie
        tie = (Fraction(low) + Fraction(high)) / 2
        if is_exactly(tie):
            return float(tie)        # the exact value IS the tie: roundTiesToEven (Fraction -> float is correctly rounded)
    raise ArithmeticError(
        f"{measure.name}.gauss({n}): {float(value)!r} lies within its error "
        f"radius of a float64 rounding boundary, so it cannot be certified "
        f"correctly rounded"
    )
```
5. `_correctly_rounded_rule` (line 633): the tail of the two-precision refusal `so some exact node or weight sits on a float64 rounding tie` is made false by this unit (an exact tie is certified identically at both precisions); reword to `so some value is unresolved at the coarser precision`. The fragment pinned by the existing test is `rounds differently` (`test_generating_measure.py` line 757): unchanged; `git grep -n -F "rounding tie" -- tests docs orpheus` finds only line 633 itself, so nothing pins the old tail.
6. The docstrings: `GeneratingMeasure.gauss` Raises block (lines 412 to 418), the class docstring attribute `exact_recurrence` (lines 286 to 295), `GeneratingMeasure.on` (lines 437 to 462: add that the remap maps the rational recurrence exactly), the module docstring items 3 and the paragraph "The recurrence coefficients are written ONCE" (lines 61 to 66 and 71 to 77), `_certified_float`'s docstring (lines 595 to 599): 2.0.9. The new refusal "does not isolate the node" is added to the Raises block: `ArithmeticError` when the rational bracket of a node holds no sign change of `pi_n` (two nodes closer than `2 x 10^-20` relative, or a node that is not a root).

**`orpheus/geometry`, `transport`, `sn` and every other consumer are NOT edited**: `gauss` keeps its signature and its exception type.

#### 2.0.6 Layers, phases and what is created

`rational_polynomial.py`: layer `numerics` (L1), no `docs/api/numerics.rst` entry (`[M]` `git grep -n -E "recurrence|generating_measure" -- docs/api/numerics.rst` prints 0 `automodule` lines for the two sibling modules, so the module follows their precedent; adding an entry would make `docs/api/numerics.rst` a file shared with `insulation`, `dead` and `roles`), imports `fractions` only (the gate `FORBIDDEN_EDGES["numerics"]` is not touched); phase: no space, field or operator (a type's arithmetic: polynomials over Q as `tuple[Fraction, ...]`, the Euclidean ring of the decision); it is a primitive and has a tier entry (E1). `recurrence.py` gains `Arithmetic`, `HighPrecision`, `RATIONAL`, `Irrational`: layer unchanged (it imports `mpmath`, `fractions`, `math`, `contextlib`, nothing from `orpheus`); phase: the arithmetic is a value passed to a builder (the precision was an argument and is now an object; no ambient state). The concept registry entry for the number-type parameter is another unit's (`registry`). `generating_measure.py` stays at tier `("frame", "generating_measure")`.

#### 2.0.7 The gates, each with its input, its first red in today's tree, and the mutation that reddens it (with the mechanism)

All gates but G7 live in `tests/gates/numerics/test_generating_measure.py` (foundation block "the correctly rounded construction refuses what it cannot certify"). Every mutation below is an in-process `monkeypatch` of a module attribute (never an edit of a file; 1.3), except M7, which loads a mutated COPY of one file under another module name (1.3's `spec_from_file_location` mechanism); the numbers marked `[M]` were produced by the in-process plugins of `scratch/ontology/w0/spec/amend/ties/` (`plugin_final.py` installs the fix of this section, `plugin_mut.py` and `plugin_lenmut.py` the mutants) over the draft `test_exact_tie_gate.py`, run serially under `.venv/bin/python -O -m pytest -p no:cacheprovider`.

**G1. The 10 tie rows** (the memo's gate, `test_an_exact_value_on_a_rounding_tie_is_rounded_half_even`, one parametrized test, ids by `str`). Each row `(panel, n, field, index, exact)`: assert first that `exact` is a tie (a `_is_tie(Fraction)` helper: `f = float(v)`; `Fraction(f) != v` and `(Fraction(f) + Fraction(nextafter(f, toward v))) / 2 == v`), then `getattr(LEGENDRE.on(*panel).gauss(n), field)[index] == float(exact)`. The oracle is the CLOSED FORM of the tied value as a `Fraction`, independent of the fix's recurrence: `c = (a + b)/2`, `s = (b - a)/2`, `2s`, `s`, `5s/9`, `8s/9`. Rows: the centre at `(0.2, 0.4)` `n = 1` and `n = 3`, at `(0.55, 0.63)` `n = 63`; one centre where half-even picks the LOWER neighbour, `(0.5118216247002567, 0.9504636963259353)` `n = 3` (every other row picks the upper: the round-up mutant's witness); the width-tie panel `(0.2615032494292951, 0.7646335475579331)`: the weight `2s` at `n = 1`, the weights `s` at `n = 2` (indices 0 and 1: the even-`n` red); the non-dyadic reference weights at `n = 3`: `5s/9` on `(0.29134508037323065, 1.218729331071271)` (indices 0 and 2) and `8s/9` on `(0.08349338193695595, 2.0129640274100478)` (index 1), the class that refutes a catalogue of dyadic values. **First red** `[M]` (today's tree): 10 failed, 2 passed (G2 and G3); every red is `ArithmeticError ... cannot be certified correctly rounded`. **With the fix** `[M]`: 12 passed.

**G2. The untied control** (`test_the_untied_sibling_is_unchanged`): `LEGENDRE.on(0.0, 1.0).gauss(3).nodes[1] == 0.5` with `_is_tie((0+1)/2)` false. Green today and with the fix (it reads no new branch; it is the instrument's positive control that an untied odd panel is unchanged).

**G3. The anti-guesser, nodes** (`test_a_single_boundary_straddle_that_is_not_the_exact_value_is_still_refused`): `monkeypatch.setattr(gm, "_CERTIFIED_DIGITS", 23)` and `gm._correctly_rounded_rule.cache_clear()` (before and after), then `pytest.raises(ArithmeticError, match="cannot be certified")` on `LEGENDRE.on(0.0, 1.0).gauss(64)`. The radius `10^(23 - d)` makes 22 of the 128 values straddle exactly one boundary at 40 digits and 0 straddle several `[M]` (the memo's `p10_single_straddle.py`, not re-run). `[M]` (this pass, `probes_ab.py`): the refusal is raised by a NODE (the first node, `0.00449331426162784`; the node decision was invoked once and the weight decision 0 times), so G3 exercises the unequal arm of the NODE decision only. It is green today and with the fix; its mutant is "the decision always says equal" (M4).

**G4. The unequal arm of the WEIGHT decision, and the decision's domain** (new; the memo's list owes this witness: G3 stops on a node). A direct call of the private decision, `_weight_is_exactly(_rational_recurrence(m, 3), 3, x, m.name)` with `m = LEGENDRE.on(0.0, 1.0)` and `x` an `mpf` at 40 digits: the weights of that rule are `5/18, 4/9, 5/18`. `[M]` (`probes_final2.py`): for `x` the centre `1/2`, `decide(Fraction(4, 9))` is `True`; for `x` an outer node `1/2 + sqrt(3/5)/2`, `decide(Fraction(5, 18))` is `True`, `decide(Fraction(4, 9))` is `False` (the gcd is non-trivial and the node is not one of its roots: the arm "gcd non-trivial, this node not a root"), `decide(Fraction(1, 3))` is `False` (trivial gcd); for `x = mpf("0.3")` (not a node) `decide(Fraction(4, 9))` raises `ArithmeticError` matching `does not isolate`; and the same calls made OUTSIDE any `mpmath.workdps` context (the default 15 digits) give the same answers (the witness of `_exact`: a conversion that re-rounds the node to the ambient precision fails the bracket). The domain row: with `t = CHEBYSHEV_T.on(0.2, 0.4)` and `r = _rational_recurrence(t, 3)`, `_weight_is_exactly(r, 3, x, t.name)` is `None` and `_node_is_exactly(r, 3)` is not `None`. **First red**: `ImportError`/`AttributeError` (the private functions do not exist in today's tree: a test-only first red by construction; the consumer-level first red of the same decision is G1 row 5). **Mutations** `[M]`: (M5) `decide` returns `len(common) > 1` (the isolation test dropped): the 12 rows of G1 to G3 stay GREEN under it (the arm is witness-less there, measured: 12 passed), and the row `decide(Fraction(4, 9))` at an outer node reddens (the mutant returns `True`); (M8) `_exact` re-rounds (`mpmath.mpf(x).man_exp`): the outside-context rows redden (`does not isolate the node`, measured on the prototype with `mpf(x)`). Mechanism: `monkeypatch.setattr(gm, "_weight_is_exactly", mutated)` for the consumer rows, and `monkeypatch.setattr(gm, "_exact", mutated)` for M8 (`[M]` M5 over the 12 rows, `plugin_lenmut.py`: 12 passed; the mutant's value on the row `4/9` at an outer node, `probes_m5.py`: `True` against the required `False`; M8 `[M]` as measured on the prototype before it was corrected: `does not isolate the node`).

**G5. The families other than Legendre** (new; the first red of the generic claim): panels and `n` of `probes_ab.py`, each asserting the centre node equals `float(c)` (a `Fraction` closed form) and, for the Jacobi `(0, 0)` row, that the weight equals the Legendre rule's on the same panel (`jacobi(0.0, 0.0) == LEGENDRE` by name, so the equality is of the two constructions, `[M]` `on(W).gauss(1)` equal): `CHEBYSHEV_T.on(0.2, 0.4)` at `n` = 1 and 3, `CHEBYSHEV_U.on(0.2, 0.4)` at `n = 3`, `jacobi(0.5, 0.5).on(0.2, 0.4)` at `n = 3`, `jacobi(-0.5, -0.5).on(0.2, 0.4)` at `n = 5`, and `jacobi(0.0, 0.0).on(*WIDTH_TIE)` at `n = 1` (weight `2s`). **First red** `[M]`: the 6 of the 6 rows raise `ArithmeticError` on today's tree (`probes_ab.py` run without the fix) and the 6 of 6 build with the fix (with `probes_ab_final.py`, the decision code of 2.0.5). Mutation M6 (`[R]`, not run): `monkeypatch.setattr(_Rational, "pi", lambda self: Fraction(22, 7))` makes the Chebyshev mass rational, so G6's domain row reddens (and the weight decision would be wrongly available).

**G6. The generic builders' law** (new; in the same file): (a) the DOMAIN table of 2.0.4 as a parametrized test: for each family and parameter row, `builder(6, RATIONAL)` returns `Fraction` for every `alpha_k` and every `beta_k` (`k >= 1`) and `beta_0` is exactly the tabulated value (`2`, `None`, `16/15`, `1`, `2`); (b) `[R]` (not run) the agreement of the two arithmetics: for each row, `mpmath.mpf(q.numerator) / q.denominator` of each `Fraction` entry equals the `HighPrecision(60)` entry to `10^-55` relative (X4 says plainly: both evaluate one body, so this agreement is tautological in the body; what it pins is the TABLE, that the rational evaluation takes the same branch as the high-precision one); (c) `[R]` (not run) the remap's law against an INDEPENDENT oracle, the moments of the uniform measure on `[a, b]`: `LEGENDRE.on(a, b).exact_recurrence(3, RATIONAL)` has `alpha'_0 == (a + b)/2`, `beta'_0 == b - a`, `beta'_1 == (b - a)^2 / 12` as `Fraction`s (the mean, the mass and the variance of the uniform measure, not the recurrence), on `(0.2, 0.4)` and `(0.2615032494292951, 0.7646335475579331)`. **First red**: `ImportError` for `RATIONAL` (it does not exist today). **Mutations**: (M6) above for (a); (M7, `[R]`, not run) the remap with `scale**2` replaced by `scale` in `beta'_k`: a mutated copy of `generating_measure.py` loaded with `importlib.util.spec_from_file_location("orpheus.numerics._gm_mutant", tmp_path)` and `exec_module` (the relative imports resolve against the `orpheus.numerics` package), then `_mutant._affine_remap(LEGENDRE, 0.2, 0.4).exact_recurrence(3, RATIONAL)` must differ from `(b - a)^2 / 12` at `beta'_1`: red on (c).

**G7. `rational_polynomial`** (`tests/gates/numerics/test_rational_polynomial.py`): the Euclid on `p = 2(x-1)(x-2)` and `q = (x-1)(x-3)` (coefficient tuples `(4, -6, 2)` and `(3, -4, 1)`): `remainder(p, q) == (-2, 2)` and `gcd(p, q) == (-1, 1)` (`x - 1`, monic); `remainder((1, 0, 0, 1), (1, 0, 1)) == (1, -1)` (`x^3 + 1` by `x^2 + 1`); a degree-0 divisor leaves the zero polynomial (`remainder((1, 0, 1), (1,)) == ZERO`: the prototype crashed with `IndexError` here before it was corrected); `remainder(p, ZERO)` raises `ZeroDivisionError`; `evaluate` against `numpy.polynomial.Polynomial` at a dyadic point (an independent oracle). All entries are `Fraction`. **First red**: `ModuleNotFoundError`. **Mutation** (M9), `monkeypatch.setattr(rational_polynomial, "scale", lambda p, c: tuple(p))` (the monic normalisation removed): `[M]` (this pass, in-process) `gcd(p, q)` becomes `(-2, 2)` and the monic leg reddens.

**The mutations of the memo, all in-process (`monkeypatch`), measured over the 12 rows of G1 to G3** `[M]` (`plugin_final.py` + `plugin_mut.py`, `ODDN_MUT` = the mode): M1 the decision removed (the new arm never taken): 10 failed; M2 the tie rounded to the LOWER neighbour: 7 failed; M3 to the UPPER neighbour: 3 failed; M4 the decision always "equal": 1 failed (G3: the unequal single-boundary refusal vanishes); M5, M6, M7, M8 as above. The mechanism of M1 to M4 is `monkeypatch.setattr(gm, "_certified_float", mutated)` with a mutant that has the same signature as the shipped function (2.0.8).

#### 2.0.8 What stays bit-identical, and the one existing test edited

- **No rule that builds today moves.** The new arm is reached only after `low != high`, where today's code always raises; a rule that builds today therefore cannot change. `[M]` (the memo's `p7_fixed.py`, 4643 pairs: the 34 panels of `p1` and the width-tie panel for `n` in 1 to 64, 300 random panels (seed 7) for `n` in 1 to 8, `[-1, 1]` at `n` = 65, 96, 128; NOT re-run by this pass, and run on the memo's hand-written Legendre form, not on the generic builders): 817 raised before and 0 after; 0 of 3826 previously building rules moved under `np.array_equal`; the 825 tied closed-form values are bit-exact `float(Fraction)`. The generic-builder form of this section was run by this pass over the 12 rows, the 6 rows of G5 and the decision rows of G4 only; the builder re-runs the numerics gate files below, in which `test_gauss_rule_fingerprint.py` pins the bytes of the shipped rules, as the evidence that the generic builders move no bit.
- **The numerics files** `[M]` (the memo, `gates_base.log` and `gates_fixed.log`, each file run singly and serially before and after, with the fix installed by a pytest plugin; NOT re-run by this pass): `test_generating_measure` 260, `test_gauss_rules_correctly_rounded` 277, `test_gauss_rule_fingerprint` 21, `test_exactness` 21, `test_registry` 112, `test_measure` 68, `test_legendre_basis` 34 passed both ways, with ONE exception: `test_generating_measure.py::test_a_rule_that_rounds_differently_at_the_two_precisions_is_refused` fails with `TypeError` under the fix. `derivations/test_characteristic_basis.py` and `geometry/test_manifold.py`: 516 passed both ways.
- **The one existing test edited.** `tests/gates/numerics/test_generating_measure.py` line 754, in `test_a_rule_that_rounds_differently_at_the_two_precisions_is_refused`: `monkeypatch.setattr(gm, "_certified_float", lambda value, digits, measure, n: float(value))` becomes `lambda value, digits, measure, n, **_: float(value)` (the stand-in must accept the keyword `is_exactly`; the memo's form, run by the memo's plugin). `[M]` `git grep -n -E "_certified_float|_rule_at|is_exactly" -- tests` finds 1 test-side use, this line (and it is the grep's positive control). The stand-in keeps the test's purpose: it disarms the certification so that the two-precision comparison fires alone; the edit does not make it vacuous (the test still requires `rounds differently`).
- `test_a_value_its_error_interval_cannot_certify_is_refused` (line 725) stays green: its 1e-6 interval holds many boundaries and is refused whatever the decision says; it is NOT the anti-guesser control (G3 is).
- Gates that read source (concept order, layer imports, pyright) see the new `fractions` import and the new module: `fractions` is the standard library and the code stays in `numerics/`; the concept-order gate needs the tier entry of E1. `[R]` unaffected otherwise; the builder runs `tests/gates/test_layer_imports.py` and `tests/gates/numerics/test_concept_order.py` (2.0.11) and reports the pyright count of the three edited production files before and after (the number must not grow: the `Number | None` mass is the typed edge of the domain).

#### 2.0.9 Docs and docstrings made false by the unit, with their tense

Each is edited in the unit's commit (a present-tense falsehood is a bug; the derivation prose is the archivist's at the wave's documentation pass, 2.0.10):
- `docs/theory/foundations/discrete_measures.rst` line 126 (the "Every Gauss rule is correctly rounded" bullet: "the float64 nearest its exact value"): add "(an exact tie to the even neighbour)". The contract paragraph (line 836): add the tie rule as a sentence of the contract: "Where an exact value lies on the midpoint of two adjacent floats, the nearest is not unique and the rule is the IEEE 754 default: the even neighbour." Step 1 (lines 892 to 898: "written once, as closed forms evaluated at the caller's working precision ... one builder per family"): add that a builder is evaluated in the arithmetic it is given, `HighPrecision(digits)` or exact rationals, so the exact tie decision of step 5 reads the same definition. Step 5 (line 936): the refusal sentence ("the construction refuses to choose rather than choose silently") is true only for an approximate value; add the exact-tie decision (node: `pi_n(T) = 0`; weight: the gcd of `pi_n` and `K_n - 1/T`; equal rounds half to even; otherwise refused) and the pushforward argument; the full text with Gautschi's equation numbers, the population and the refuted candidates is the archivist's.
- `orpheus/numerics/generating_measure.py`: the module docstring, item 3 of the construction (lines 61 to 66: "a value that cannot be certified is refused, never guessed" gains "unless it is exactly a tie: then it is rounded half to even") and the paragraph at lines 71 to 77 ("The recurrence coefficients are written ONCE, in arbitrary precision (`exact_recurrence`)": "in an arithmetic: high precision or the rationals"); the class docstring attribute `exact_recurrence` (lines 286 to 295: "`(n, digits) -> (alpha, beta)`, each a list of `n` `mpmath.mpf`" becomes `(n, arithmetic) -> (alpha, beta)` with the mass possibly `None` in `RATIONAL`); `gauss`'s Raises block (lines 412 to 418: add the exact-tie rule and the non-isolation refusal); `on` (lines 437 to 462: "The transform is applied in the working precision" stays true; add "and exactly in the rationals"); `_certified_float`'s docstring (lines 595 to 599: "A value within rho of a rounding boundary is refused rather than guessed" gains the exact-value exception); the `_correctly_rounded_rule` message tail (E3 item 5).
- `orpheus/numerics/recurrence.py`: the module docstring ("each function returns `(alpha, beta)`, two lists of `n` `mpmath.mpf`"; "Precision. Every coefficient is evaluated at an explicit number of decimal digits (the precision is an argument, never mpmath's ambient setting)", lines 11 and 32 to 35) and the comment above `RecurrenceCoefficients` (lines 53 to 59): the precision argument is now an arithmetic.
- `tests/gates/numerics/test_generating_measure.py`: the docstring of the edited test names `_certified_float` by its old signature only in its body, which is the line-754 edit.
- Past tense that stays: the paragraph in `discrete_measures.rst` that says the platform eigensolver once produced the floats (the 2026-10-01 history).

#### 2.0.10 The ERR entry (a certificate correct for approximate values, applied to an exact value on a tie)

The defect class is "a certificate correct for approximate values applied to an exact value that sits on a tie" (failure mode 4/6 hybrid in the taxonomy of the catalogue); it hid because the reference rule's only exact node is 0, where `|0| rho = 0`, and its rational weights are representable or non-dyadic (2.0.1). The catalogue file is `docs/theory/verification/error_catalog.rst`; an entry is `.. error-entry:: ERR-NNN` with `:title:` and a body; the id is the NEXT FREE one at integration (`[M]` `ERR-127` at the base hash: 126 declared entries, the last `ERR-126`; no other unit of this specification mints one: dead-A records ERR-118 as retired and `roles` edits ERR-075; arm 4 below decides if another landing took it). The catching test carries `@pytest.mark.catches("ERR-NNN")` on G1 (the parametrized tie rows: the marker decorates the test function, so it covers every row).

`[M]` which gate needs the entry (read `tests/gates/test_error_catalogue_reconciles.py` at the base hash): arm 1 `test_every_marker_names_a_declared_entry` is RED for a `catches` marker whose id is not declared (and the build is `sphinx-build -E -W`, where nexus warns per unresolved marker), arm 2 `test_every_declared_entry_has_a_catching_test` needs the marker for a declared entry, arm 3 `test_the_generated_index_matches_the_corpus` is red until `.claude/skills/vv-principles/error_index.md` is regenerated (the orchestrator, 7.3: a by-design red of the unit that the generated file closes), arm 4 `test_the_id_sequence_has_no_gaps_or_duplicates` fixes the id, arm 5 `test_every_test_the_catalogue_cites_exists` requires every `tests/...py::name` cited in the entry to exist now. So the unit DOES add the entry's stub (the marker without it is red): `.. error-entry:: ERR-127`, `:title:` one sentence naming the defect, `**Status:** FIXED, merged @ <hash>` (filled at integration), `**Module:** numerics`, `**Failure mode:**` one line, and one `**Caught by:**` line citing `tests/gates/numerics/test_generating_measure.py::test_an_exact_value_on_a_rounding_tie_is_rounded_half_even` (the function exists after the unit). The full text (the defect, how it hid, the fix, the lesson) is the archivist's at the wave's documentation pass; the stub states that. The builder runs `tests/gates/test_error_catalogue_reconciles.py` and `tests/gates/test_vv_harness_audit.py` (2.0.11).

#### 2.0.11 Retirement-audit items and done-when

Items: A.1/A.2 (`_at_digits` and the six builders: the census of 2.0.2); B.4 (`nexus dead-references` after the docs rebuild: `_at_digits` is retired, so the docs of `recurrence` are read); B.6 (no string names these symbols); B.9 (a labelled equation: none is added); B.11 (the message tail of E3 item 5 has no pinning test: `git grep -n -F "rounding tie" -- tests` is empty); C.13 (the `catches` marker of 2.0.10, with its entry in the same commit); D.14/D.17 (every mutation of 2.0.7 reddens a gate; M5 is the one that the 12 rows cannot, and G4 carries it); D.18 (the mechanism that kept the exact value from ever needing a decision, the reference rule having only the node 0 as rational, is gone: G4 and G1 are the witness); E.19 (2.0.9). The two catalogue arms and the concept-order tier are the by-design reds, each closed in the same commit.

Done when (each is one command):
1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics/test_generating_measure.py tests/gates/numerics/test_rational_polynomial.py` exits 0, and the output names G1 (10 tie rows), G3, G4, G5, G6.
2. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics/test_gauss_rule_fingerprint.py tests/gates/numerics/test_gauss_rules_correctly_rounded.py tests/gates/numerics/test_exactness.py tests/gates/numerics/test_registry.py tests/gates/numerics/test_measure.py tests/gates/numerics/test_legendre_basis.py tests/gates/numerics/test_concept_order.py tests/gates/test_layer_imports.py tests/gates/test_test_layers.py` exits 0.
3. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_error_catalogue_reconciles.py tests/gates/test_vv_harness_audit.py` exits 0 except arm 3, which is red until the orchestrator regenerates the index (7.3); the builder reports exactly that one red.
4. `.venv/bin/python -O -c "from orpheus.numerics.generating_measure import LEGENDRE; r = LEGENDRE.on(0.2, 0.4).gauss(1); print(r.nodes[0])"` prints `0.30000000000000004` (today: `ArithmeticError`); and the 34-panel population probe of 2.0.1 (the builder writes the probe from the description: panels `(-1, 1)`, `(-0.3, 0.3)`, `(-2.5, 2.5)`, `(0, 1)`, `(1, 2)`, `(0.2, 0.4)`, `(0.55, 0.63)`, `(3, 7)` and the carried cylinder and sphere panels of `tests.gates.cp._cp_s4_helpers.cp_mesh(chart).edges`, `n` in 1 to 64) reports 0 raises in all of its 2176 pairs.
5. The AST census of 2.0.2, re-run, reports the six builder names and `exact_recurrence` in exactly the 3 files named there, and `_at_digits` in 0 files (positive control: `exact_recurrence` found in `generating_measure.py`).
6. `sphinx-build -E -W docs docs/_build/html` exits 0; `.venv/bin/nexus dead-references` has been read for `_at_digits`.
7. Mutation report: M1 to M9 of 2.0.7, each reddens the named gate (the builder lists the counts; M6 and M7 were not run by the specification).

#### 2.0.12 Unmeasured (carried, not hidden)

- The memo's 4643-pair before/after, its 12000-pair random census (`p6`), its independent mpmath oracle (`p8`: 0 of 4740 values differ for the untied values; the tied values are witnessed by the closed forms and G1), its gcd cost table (`p9`) and the per-file counts of 2.0.8 were NOT re-run by this pass; what this pass re-ran: the 12-row gate (today: 10 failed, 2 passed; fixed: 12 passed; mutations 10, 7, 3, 1 failed) and the population probe `p1` (384 of 1088 odd pairs raise, 0 of 1088 even, 0 mismatches of 2176), on a Python 3.14 `.venv`; plus the generic builders (4368 of 4368 entries), the table of 2.0.4, G4 and G5's rows and M5 (the prototype, `scratch/ontology/w0/spec/amend/ties/`).
- The unequal single-boundary arm for a WEIGHT has no consumer-level witness (no shipped rule straddles one boundary at the shipped radius with an unequal weight); G4 calls the decision directly. A non-equal value within 1e-30 relative of a tie is `[R]` vanishingly rare; the escalation of precision is not built.
- The decision's isolation bracket `10^-20` relative assumes node gaps above `2 x 10^-20` relative: `[R]`, true for every `n` up to a few thousand of the families here; beyond, a straddle raises `does not isolate the node`, a refusal and not a wrong answer.
- `[R]` no reference weight at an irrational node is rational for `n >= 4` (screened numerically for dyadic values by the memo's `p4`, not proved); the fix does not depend on it (the gcd decides every case).
- Timing of `tests/gates/numerics` and of a cold `gauss(63)` on a tied panel (0.43 s in the memo's `p7`): the builder reports the numerics suite's time before and after.

### 2.1 File set (by tree)

Existing, hard (edited):
- production: `orpheus/cp/solver.py`.
- tests: `tests/gates/cp/test_cp_boundary_carry.py` (3 rows re-baselined), `tests/gates/numerics/test_concept_order.py` (a tier entry), `tests/gates/test_layer_imports.py` (the production-insulation function and legs; the file is edited in the order roles, insulation, dead-B, and the later `memo` and `registry` units after all of part A; section 7).
- docs: `docs/architecture/layering.rst` (the sentence documenting the gate contract and the `derivations` rows; find with `git grep -n -E "derivations|FORBIDDEN_EDGES" docs/architecture/layering.rst`), `docs/api/numerics.rst` (the API entry of the new module), `docs/theory/methods/collision_probability.rst` (the cylinder kernel's provenance paragraph, 2.5).

New:
- `orpheus/numerics/bickley.py` (2.3).
- `tests/gates/cp/test_cp_composite_rule.py`, `tests/gates/cp/test_cp_chord_spelling.py`, `tests/gates/cp/test_cp_kernel_is_the_numerics_ki3.py` (2.4).
- `tests/gates/derivations/test_cp_ki3_copy.py` (2.4; NOT under `tests/gates/numerics/`: `FORBIDDEN_EDGES["numerics"]` forbids `derivations`, so a numerics-directory test cannot import the reference; `FORBIDDEN_EDGES["derivations"]` forbids L2, L3 and L4 only, so a derivations-directory test may import `orpheus.numerics` `[M]` measure_gates.md gate 3).

Read-only (listed so that nobody edits them): `orpheus/derivations/**` (0 edits); `orpheus/geometry/chord.py` (`ConcentricPartition.chord`, consumed as is); `tests/gates/test_test_layers.py` (reads `FORBIDDEN_EDGES`, must stay untouched); `tests/gates/transport/carry_head.npz` (see decision 2.6.1: the three rows are overridden in the test, the shared fixture is not edited); `orpheus/numerics/quadrature/__init__.py` (the composite rule is NOT exported from the quadrature package: that file belongs to dead-A, whose row 73 edits its import and `__all__`).

Prose that must be re-read for tense (no gate decides it; `retirement-audit` B.8 and E.19): the prose inside `orpheus/cp/solver.py` at four places, found by the quoted anchors: the comment under `class CPParams` beginning `# n_ki_table, ki_max: retained for backwards compatibility` (names `orpheus.derivations.continuous.flat_source_cp.geometry` and #94); the comment block beginning `# Q4 of the quadrature-architecture rollout` (names `orpheus.derivations.common.quadrature.composite_gauss_legendre`); the comment inside `_setup_radial_quadrature` beginning `# ``composite_gauss_legendre`` declares`; the docstring of `_setup_cylindrical` beginning `Cylindrical: canonical Ki_3 kernel` (says the solver "shares the Chebyshev interpolant ... bit-identical kernel evaluations"). The pages that may carry a present-tense claim that production CP shares the reference kernel (a candidate list by word match; read each, rewrite the present-tense ones, leave history): `docs/api/derivations.rst`, `docs/architecture/conceptual_view.rst`, `docs/theory/foundations/chart_and_chord.rst`, `docs/theory/methods/method_of_characteristics.rst`, `docs/theory/references/characteristic.rst`, `docs/theory/references/fn_method.rst`, `docs/theory/references/peierls.rst`, `docs/theory/references/peierls_nystrom.rst`, `docs/theory/verification/collision_probability.rst`, `docs/theory/verification/error_catalog.rst`, `docs/theory/verification/reference_solutions.rst`, `docs/theory/verification/sn.rst`. Command: `git grep -n -i -E "bit-identical|shares? the (chebyshev|interpolant)|production cp|cp/solver|composite_gauss_legendre|chord_half_lengths|_ki3" -- docs`. `tests/gates/cp/test_sphere_pss.py` (a docstring "agree by construction (X4)" naming `chord_half_lengths`) stays true after the swap (one shared upstream) and is not edited.

### 2.2 The production imports today, and the three edits in `orpheus/cp/solver.py`

`[M]` Instrument: AST over the 260 production files (every tracked `*.py` under `orpheus/` outside `orpheus/derivations/`, relative imports resolved): 3 imported names, 3 statements, 1 file, all at module level, none under `TYPE_CHECKING`:

```
orpheus/cp/solver.py:49   from orpheus.derivations.common.kernels import chord_half_lengths
orpheus/cp/solver.py:50   from orpheus.derivations.common.quadrature import composite_gauss_legendre
orpheus/cp/solver.py:51   from orpheus.derivations.continuous.flat_source_cp.geometry import _ki3_mp as _ki3_kernel
```
Cold import: one fresh interpreter per top-level package importing every module under it, counting loaded `orpheus.derivations*` modules: only `orpheus.cp` loads any (11 modules); 17 other packages load 0 with 0 import errors. No dynamic import with a `derivations` string literal exists in production (5 calls shaped like `import_module`/`__import__`, all in the two traced-memo files).

Edits (delete the three imports; add `from orpheus.geometry.chord import ConcentricPartition`, `from orpheus.geometry.line import Line`, `from orpheus.numerics.generating_measure import LEGENDRE`, `from orpheus.numerics.bickley import ki3`; `Chart` is already imported). `[M]` The method today is `CPMesh._setup_radial_quadrature` (`orpheus/cp/solver.py`, the lines `radii = self.mesh.edges[1:]`, `q = composite_gauss_legendre(self.mesh.edges.tolist(), p.n_quad_y)`, `in_plane = np.array([0.0, 1.0, 0.0])`, `self._y_pts = q.pts`, `self._y_wts = q.wts * Chart(self.mesh.coord).beam_density(q.pts, in_plane)`, `self._chords = chord_half_lengths(radii, self._y_pts)`); the in-plane direction is already spelled once there and is reused.

1. **The partition, the rule and the chord, in one method, with one parse of the edges.** `ConcentricPartition.__post_init__` is the only edge check: it refuses fewer than two breakpoints, a non-finite or non-increasing sequence (message fragments `strictly increasing` and `at least two`, from `orpheus/geometry/chord.py`) and, on a cylinder or sphere, `r_0 < 0`. The chart is built once. The rule reads `partition.breakpoints` and so does the chord. Replace the body after `p = self.params` by:

```python
chart = Chart(self.mesh.coord)
partition = ConcentricPartition(chart, tuple(self.mesh.edges.tolist()))     # the one parse of the edges
pts, wts = _composite_gauss_legendre(partition.breakpoints, p.n_quad_y)
# the measure on lines over the impact parameter y ... (the comment above `in_plane` stays)
in_plane = np.array([0.0, 1.0, 0.0])
self._y_pts = pts
self._y_wts = wts * chart.beam_density(pts, in_plane)
points = np.stack([pts, np.zeros_like(pts), np.zeros_like(pts)], axis=-1)
line = Line.through(points, np.broadcast_to(in_plane, points.shape))
n = partition.n_regions
self._chords = partition.chord(line).slot_length[..., -n:].T       # outbound half: regions 0..n-1, shape (n, n_y)
```
`radii` and the `.tolist()` comment are deleted. `[M]` this array equals `chord_half_lengths(edges[1:], y)` bit for bit, 1024 of 1024 entries, on the cylinder and the sphere (4 regions, 256 impact parameters; the inbound half reversed, `slot_length[..., :n][..., ::-1].T`, is also equal) (`scratch/ontology/w0/e2_chord_compare.py`, log `e2_chord_compare.log`). The reason is structural: `_radial_chord`, the body `ConcentricPartition.chord` runs, forms its slots from the same `region_half_lengths` of the same `half_chord` that `chord_half_lengths` calls, so the agreement is tautological (X4) and a chord swap has no behavioural first red. The decoding of the slot layout by position (`[..., -n:]`: n inbound slots, the cavity slot, then n outbound slots, `orpheus/geometry/chord.py`) is the producer's convention read at the consumer; it is pinned by G4 and recorded as a follow-up issue (8.1): select the outbound slots through `slot_region` and the sense, or an outbound accessor on `Chord`, in a later wave.

   **Domain, measured.** `[M]` (`scratch/ontology/w0/spec/amend/hollow_cp.py`, run under `-O`): a hollow cylinder `StructuredGeometry(coord=CYLINDRICAL, breakpoints=(0.5, 1.0, 2.0), mat_ids=(1, 0), face_roles=(0, 0))` meshed by `Mesher(g).partition(CellsByCount.uniform_width(2)).mesh` is refused at `CPMesh` construction: without a declared condition `UnboundTagError` (the inner face carries tag 0), with `{0: WallCondition.vacuum()}` `BoundaryRefused` ("CP reads only the outer face of a hollow body: the declared inner condition would be dropped", `_dropped_face`, a `SCOPE-BOUNDARY[guard]` of #513). So a hollow `CPMesh` is not constructible through the tag-bearing inner face at the base hash (one geometry, one tag; a hollow body with an untagged inner face was not tried), and no k exists to capture for it (#513 refuses the declaration of the inner face, so a hollow `CPMesh` cannot be built today). The discharge of decision I4 is therefore the function-level hollow leg of G4 below, and no `edges[0] != 0` refusal is added in `CPMesh`. The hollow case that exists is the function level: `[M]` (structure check, `hollow.py`) for breakpoints `[0.1, 0.3, 0.5]`, cylinder and sphere, 7 values of y in `[r_0, r_n]`, `ConcentricPartition.chord` equals `chord_half_lengths` bit for bit (`array_equal` True; the equality needs `y >= r_0`, and the y-rule starts at `edges[0]`, so no line below `r_0` is drawn). G4 carries that function-level hollow leg.

2. **The composite rule.** Delete the call `composite_gauss_legendre(...)` and its `.tolist()` comment (done in item 1). Add the module-level helper below, above `class CPMesh`. It is the realisation, until the primitive exists, of the direct sum of the panels' measures, `⊕_k LEGENDRE|[b_k, b_{k+1}]`: `DiscreteMeasure.__add__` (`orpheus/numerics/measure.py`, "Direct sum ... panelled Gauss-Legendre") is that sum but requires equal supports, and each panel's support is its own `Interval(a, b)` (`generating_measure.py`), so `[M]` (structure check) summing `LEGENDRE.on(a, b).gauss(8)` over abutting panels raises `ValueError: direct sum requires equal spaces, got '[0,0.2]' and '[0.2,0.4]'`. No per-panel cache is added: `_affine_remap` and `_correctly_rounded_rule` in `generating_measure.py` are already `lru_cache`d, and `gauss` returns fresh arrays built from cached tuples (the 84 s uncached run happened with both caches active), so a third cache would hold only an array allocation.

```python
# ELEGANCE-DEBT[guard] #<W8a-issue> (W8a): a composite rule is the direct sum of its panels' measures, and DiscreteMeasure.__add__ refuses abutting panels.
# Retires when DiscreteMeasure.__add__ admits a direct sum over abutting intervals (W8a): the helper is then `sum(LEGENDRE.on(a, b).gauss(n) for a, b in panels)`.
def _composite_gauss_legendre(breakpoints: Sequence[float], n_per_panel: int) -> tuple[NDArray, NDArray]:
    r"""Gauss-Legendre panels abutting on the breakpoints, as one rule: (nodes, weights). Laws: see 2.3."""
    panels = [LEGENDRE.on(lo, hi).gauss(n_per_panel) for lo, hi in zip(breakpoints[:-1], breakpoints[1:])]
    return np.concatenate([m.nodes for m in panels]), np.concatenate([m.weights for m in panels])
```
   The tag spelling: `ELEGANCE-DEBT[guard]` is the only token spelling the ledger gate `tests/gates/test_elegance_debt_is_tagged.py` reads (`[M]`: 38 occurrences of `ELEGANCE-DEBT[...]` in `orpheus/`, `.claude/` and `docs/development/`, every one `[guard]`; the gate requires `#\d+` on the token line and a `retires|retired|until|when` sentence within three lines). The helper is machinery debt rather than a guard and the grammar has no other token; the token is used so that the ledger lists it, and the issue number is the W8a issue filed in 8.1 BEFORE the unit's commit (so the line carries a real number). `import functools` is not needed; `from collections.abc import Sequence` and `from numpy.typing import NDArray` are added if absent. The final body is the builder's, the LAW (2.3) is the contract.
3. **The Bickley kernel.** `_setup_cylindrical` uses `ki3` in the two places that read `_ki3_kernel` today (`self._kernel = ...` and `self._kernel_zero = ...(np.zeros(n_y))`); delete the `_ki3_kernel` alias. Rewrite the four prose places of 2.1.

### 2.3 Every law whole, and the raw inputs, of the new functions

**`_composite_gauss_legendre(breakpoints, n_per_panel) -> (nodes, weights)`** (production, `orpheus/cp/solver.py`; layer L3, method-private; no phase: a quadrature helper). Raw inputs: `breakpoints`, the `partition.breakpoints` tuple of floats (already validated by `ConcentricPartition`: at least two, finite, strictly increasing, `r_0 >= 0`; the helper makes no check of its own, one parse per input); `n_per_panel`, an `int` (the user's `CPParams.n_quad_y`, default 64). The laws:
- (a) Shape: the refusals belong to the partition, not to the helper: `ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), (0.0,))` and `(0.0, 0.5, 0.5)` raise `ValueError` (fragments `at least two` and `strictly increasing`, from `orpheus/geometry/chord.py`). `[M]` no test pins that refusal today (`git grep -n "finite and strictly increasing\|at least two" -- tests/gates/geometry` prints nothing), and the unit now relies on it as the one edge check, so G3 carries the two refused legs (`retirement-audit` D.18: the check that kept `mesh.edges` valid is now load-bearing).
- (b) Size: `nodes.shape == weights.shape == ((len(breakpoints) - 1) * n_per_panel,)`.
- (c) Order: the nodes are strictly increasing and panel `k` lies in `[b_k, b_{k+1}]` (the concatenation preserves the panel order).
- (d) Exactness PER PANEL, the law the rule exists for (the global law "exact for degree <= 2n-1 on the whole interval" is weaker and a single global rule meets it): for a function that is a polynomial of degree at most `2 n_per_panel - 1` on each panel, `sum(weights * f(nodes))` equals the exact integral over `[b_0, b_last]` to rounding. The witness is the kinked integrand `f(x) = |x - 0.4|^3` on edges `[0, 0.2, 0.4, 0.55, 0.63]` at `n_per_panel = 8`, with the exact value `0.4^4/4 + 0.23^4/4` computed in the test from the closed form; `[M]` (structure check) the panel rule's error on it at n = 8 is 0.0 while one 32-point Gauss-Legendre rule has error 7.3e-7, so the leg separates a panelled rule from a global one (the former `x^7` leg is exact for both and is dropped).
- (e) Agreement with the reference rule (an oracle on an independent construction, `derivations.common.quadrature.composite_gauss_legendre`, numpy's `leggauss` mapped by float arithmetic): `np.allclose(..., rtol=1e-10, atol=0)` on nodes and weights. `[M]` on edges `[0, .2, .4, .55, .63]`: n = 8 max relative node difference 1.3e-15, weight 2.6e-15; n = 64 nodes 8.5e-14, weights 1.2e-12; not bit-equal (the reference's small end weights carry about 1e-12 relative error, the production rule is correctly rounded). The tolerance has about 100 times the measured margin.
- (f) Total in `n_per_panel`: this is the law of the precondition unit `quadrature-ties` (2.0.3), not a law of the helper; `LEGENDRE.on(a, b).gauss(n)` builds for every `n` on every panel, and there is no parity clause.

**`orpheus.numerics.bickley.ki3(tau) -> NDArray`**: `Ki_3(tau) = integral_0^{pi/2} cos^2(theta) exp(-tau / cos(theta)) dtheta` (the Bickley-Naylor function of order 3), evaluated as `exp(-tau)` times the degree-63 Chebyshev interpolant of the scaled kernel `f(tau) = exp(tau) Ki_3(tau)` on `[0, 50]`, and clamped to 0 for `tau > 50`. Raw input: `tau`, array-like of floats, read with `np.asarray(dtype=float)`; expected `tau >= 0` (it is an optical thickness). Output: a float64 array of the input's shape. Laws: (a) pointwise equal, bit for bit, to `orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` on every input (the SAME interpolant: same 64 Chebyshev-Lobatto nodes, same coefficients; a degree-62 mutant differs at 2003 of 2004 grid points `[M]` `m_ins3.py`); (b) the defining law, independent of the reference copy: for every `tau` in `[0, 49.9]`, `ki3(tau)` equals the mpmath integral `Ki_3(t) = integral_0^inf (1+u^2)^{-2} exp(-t sqrt(1+u^2)) du` at 30 decimal digits to relative error `3e-6`; `[M]` maximum relative error 2.6e-6, at `tau = 0`, over the 81 points `np.concatenate([[0.0], np.logspace(-8, np.log10(49.9), 80)])` (`scratch/ontology/w0/e2_ki3_accuracy.py`, `_ki3_mp interpolant` row; the exact values there are the reference's `ki_n_mp(3, t, 30)`, the gate writes the integral itself so that it does not import the reference); (c) `ki3(tau) = 0.0` exactly for `tau > 50` (the clamp; `Ki_3(50)` is about 3e-23); (d) the output of an empty input is empty; (e) a negative `tau` is not refused here (the original does not refuse it; `ki3` is the kernel evaluator of a solver whose inputs are optical thicknesses; refusing would add a guard on the hot path, and the derivation module's own `ki_n_mp(3, x)` refuses `x < 0` at build time only for the 64 build nodes, which are non-negative). The raw build input: the 64 Chebyshev-Lobatto nodes of `numpy.polynomial.Chebyshev.interpolate(func, deg=63, domain=[0.0, 50.0])`, evaluated with the mpmath integral at 30 decimal digits, times `float(np.exp(t))`.

Module `orpheus/numerics/bickley.py`. Contents: the constants `_KI3_TAU_MAX = 50.0`, `_KI3_DEG = 63`, `_KI3_DPS = 30`; the private mpmath integrand of order 3 only (what `n = 3` needs: the original's `_ki_integrand(n, x)` and `ki_n_mp` are general in `n` and only `n = 3` is consumed, so the copy is written for `n = 3`); the cached builder `ki3_scaled_chebyshev()` under `functools.lru_cache(maxsize=1)` (the original's idiom, about 0.3 s on first call: lazy, so importing the module costs nothing); and the public `ki3`. Stored coefficients are NOT used (decision 2.6.4): the interpolant is rebuilt from the integral exactly as the reference rebuilds it, and the witness compares the two coefficient arrays bit for bit. Layer: `numerics` (L1), imports only `functools`, `numpy`, `mpmath`; nothing from `orpheus` (so it imports no `derivations`, and `numerics` already imports `mpmath` in `generating_measure.py` and `recurrence.py`). Phase: a pure function with a cached constant; it has no space, field or operator: it is a special function, entered in the concept registry (another unit) as an extension concept `kernel: special function`; W9 re-spells collision probability as an integral operator and will own it. Tier in `tests/gates/numerics/test_concept_order.py::TIERS`: the primitives tier, the one listing `("evidence", "face_layout", "gauge", "index_set", "instance_cache", "metric", "recurrence", "registry", "tangential", "units", "vector")`; add `"bickley"` at the head of that tuple (alphabetical). The module is imported by path (`orpheus.numerics.bickley`) and is NOT re-exported from `orpheus/numerics/__init__.py`. No `verifies` marker for the defining law is attached: `Ki_3` has no labelled equation (`:eq:ki3-def` in `docs/theory` defines `Ki_2`, #348 open), and the new paragraph of 2.5 cites #348.

### 2.4 The gates, each with its input, its first red in today's tree, and the mutation that reddens it (with the mechanism)

Every walk below takes a root path parameter; a mutation of a source file is spelled by copying the file to a temporary directory, mutating the copy and calling the function on the copy's source (or `importlib.util.spec_from_file_location` for a cold-interpreter leg); a mutation of behaviour is an in-process `monkeypatch`.

**G1. The production-insulation table** (`tests/gates/test_layer_imports.py`, beside the reference-insulation block that starts at the constant `_REFERENCE_INSULATION_EXEMPT`; its model is `_reference_insulation_violations`, whose scan is `_imported_targets(rel, src) + _dynamic_imports(src)`, and `_INSULATION_REFUSED_LEGS` / `_INSULATION_ADMITTED_LEGS`). The refused set is the module constant `L0_PACKAGES` (`tests/gates/test_layer_imports.py`, a `frozenset` equal to `{"derivations"}` today; the new block spells no literal `"derivations"`); the function is `_production_insulation_violations(rel, src)`, which returns a violation for every import (`_imported_targets`) or string-literal dynamic import (`_dynamic_imports`, the existing scanner, called, not respelled) whose target's second dotted component is in `L0_PACKAGES`. Why a function beside the table and not new `FORBIDDEN_EDGES` rows: with plain rows (every row but `derivations` gains `L0_PACKAGES`), `tests/gates/test_test_layers.py` reports 100 new (file, package) rows in 95 test files (by directory `sn` 57, `cp` 13, `homogeneous` 7, `transport` 7, `mc` 5, `moc` 4, `diffusion` 2), the method gates reading their references `[M]` `m_ins9.py`; `FORBIDDEN_EDGES` has 19 rows keyed by top-level package, `derivations` is forbidden by 7 (`values numerics geometry data mesh specification reference`) and by none of `transport`, the L3 packages or `plotting`; and `_top_level_package` returns `"plotting.py"` for the module file `plotting.py`, so a per-row scheme cannot reach it (`ORPHEUS_TOP_LEVEL - FORBIDDEN_EDGES.keys() == {"plotting"}`). The population: every `*.py` under `orpheus/` whose top-level package is not in `L0_PACKAGES` (260 files `[M]` `m_ins1.py`, including `plotting.py` and `orpheus/__init__.py`), found by `_iter_python_modules`. The old coverage leg (a) is dropped: `_iter_python_modules` is an `rglob` and `ORPHEUS_TOP_LEVEL` is computed from the same directory, so both of its checks hold on every tree and no input can redden it (X1); coverage holds by construction, which the docstring states. The f-string-head import (`import_module(f"orpheus.derivations....")`) is a declared blind spot with the computed names; adding it to `_dynamic_imports` once would serve both insulation gates and is a follow-up, not this unit. Refused legs (each a synthetic source, run through the function): `cp/x.py` top-level `from orpheus.derivations.common.kernels import chord_half_lengths` (the input that is red at the base tip), `sn/x.py` late import inside a function, `transport/x.py` `from ..derivations import get` (relative), `sn/x.py` under `if TYPE_CHECKING:`, `plotting.py` top-level, `sn/x.py` `importlib.import_module('orpheus.derivations.common.kernels')`. Admitted legs: `sn/x.py` importing `orpheus.numerics`; `derivations/x.py` importing `orpheus.geometry` (the direction that stays legal: the target's package is not in `L0_PACKAGES`). Run-time edge not covered: `orpheus.memo`'s `_resolve` (`traced_memo.py`, census 1a) imports stored client names, and three `derivations` files are clients: a run-time production-to-`derivations` edge the AST cannot see; G2's cold import covers import-time executions only. **First red** `[M]` (`m_ins1.py`, AST independent of the gate's scanner): 3 statements in 1 file of 260: `cp/solver.py:49`, `:50`, `:51`. **Mutation**: call `_production_insulation_violations` on the source of `git show b68f8b8c:orpheus/cp/solver.py` (a copy of the base file, or a copy of the landed file with the first import line re-added, written to a temp directory), require red on that line.

**G2. The cold-import leg** (same file; idiom `_cold_loaded_modules(code)`, a fresh interpreter returning `sys.modules` after the code). For each name in `ORPHEUS_TOP_LEVEL - L0_PACKAGES`, one fresh interpreter imports every module under it and the gate requires that no loaded module name starts with `orpheus.<package in L0_PACKAGES>`. **First red** `[M]` `m_ins2.py`: `import orpheus.cp.solver` loads 11 `orpheus.derivations*` modules among 164 `orpheus` modules; `import orpheus.sn` 0 of 190; `import orpheus.numerics.quadrature` 0 of 59; control: `import orpheus.derivations.common.kernels` loads 6, so the filter sees members. **Mutation** `[M]` (`g2.py` of the structure check): load a mutant copy of `cp/solver.py` (the three imports re-added, in a temporary directory) as `orpheus.cp.solver` with `importlib.util.spec_from_file_location("orpheus.cp.solver", tmp_path)` and `exec_module` inside the cold interpreter, and require the loaded-module test red (it printed 6 `orpheus.derivations` modules); a shadow directory first on `sys.path` does NOT work (`orpheus.cp.solver.__file__` stays the tree's). What it does not catch: an import executed lazily inside a call.

**G3. The composite-rule witness** (`tests/gates/cp/test_cp_composite_rule.py`; a `tests/gates/cp/` test may import `orpheus.cp` and `orpheus.derivations`, a `tests/gates/derivations/` test may not import `orpheus.cp`, an L3 package). Calls `orpheus.cp.solver._composite_gauss_legendre` with the derivations rule as the oracle: laws (b), (c), (e) of 2.3 at n = 8 and n = 64 on `partition.breakpoints` of the 4-panel edges `[0, 0.2, 0.4, 0.55, 0.63]`; the per-panel exactness leg (d) with `|x - 0.4|^3` at n = 8; the odd-n leg of (f) at n = 3 and n = 63 on the same edges (agreement with the reference to `rtol=1e-10`); and the two refused partition legs of law (a). **First red**: `ImportError` (the helper does not exist before the unit lands); the odd-n legs additionally raise `ArithmeticError` on the base tree and are green once `quadrature-ties` has integrated (it integrates first, so a builder meets them green; the base-tree red is the evidence that the order matters). **Mutations**: drop the last panel in a monkeypatched copy of the helper (the size differs; measured); `n_per_panel = 63` passed as 64 (the size differs); replace the per-panel loop by one `LEGENDRE.on(b[0], b[-1]).gauss(n_per_panel * (len(b) - 1))` (the per-panel leg (d) reddens, `[M]` error 7.3e-7 against 0.0); for the refused legs, a monkeypatch of `ConcentricPartition.__post_init__` that skips the check turns them red.

**G4. The chord-spelling witness** (`tests/gates/cp/test_cp_chord_spelling.py`). Builds `CPMesh` for a cylinder and a sphere of 4 regions at the default `n_quad_y`, and requires `mesh._chords` equal bit for bit to `orpheus.derivations.common.kernels.chord_half_lengths(mesh.mesh.edges[1:], mesh._y_pts)` (`np.array_equal`), and the shape `(4, n_y)`; plus ONE function-level hollow leg: for breakpoints `[0.1, 0.3, 0.5]` on the cylinder and the sphere and 7 values of y in `[0.1, 0.5]`, `ConcentricPartition(...).chord(line).slot_length[..., -n:].T` equals `chord_half_lengths(breakpoints[1:], y)` bit for bit (`[M]` structure check, `hollow.py`). A hollow `CPMesh` cannot be built today (2.2 item 1, #513), so no hollow k exists to capture, and this function-level leg is the discharge of decision I4. State in its docstring that the agreement is tautological in its upstream (X4: both read `region_half_lengths` of `half_chord`) and has no behavioural first red; what it pins is the slot choice and the transposition (the outbound half, regions 0..n-1, `.T`). **First red**: none behavioural, by X4 (the gate lands green and says so: its first red is a test-only input, the mutation). **Mutation**, in-process: monkeypatch `orpheus.geometry.chord.Chord.slot_length` to return its array with the inbound half selected (`[..., :n]`, not reversed) for the `CPMesh` build, or, equivalently, run the leg's comparison against the `[..., :n]` slice; it differs from the reference. The behavioural pins of the chord are the carry rows that do NOT move, `sphere-white_1.0` among them, which passes bit-exact under the rule swap `[M]` (a wrong chord should redden it; `[R]`, no in-process chord mutation was run against the carry file; the builder runs one by `monkeypatch.setattr` of `CPMesh._setup_radial_quadrature` with `_chords` transposed and requires a red on `sphere-white_1.0`).

**G5. The Ki3 witness** (`tests/gates/derivations/test_cp_ki3_copy.py`). (a) `np.array_equal(orpheus.numerics.bickley.ki3(grid), orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp(grid))` on a 2004-point grid: `np.linspace(0, 50, 2000)` plus `[50.0, 50.0000001, 1e-300, 60.0]`; (b) `np.array_equal(ki3_scaled_chebyshev().coef, orpheus.derivations.continuous.flat_source_cp.geometry._ki3_scaled_cheb().coef)` and equal `.domain`; (c) the defining law of 2.3 (b): `ki3(t)` against the mpmath integral at 30 digits on the 81 points of 2.3, relative error `<= 3e-6` (this leg is independent of the reference copy; the two copies of one approximation agree by construction, X4, so legs (a) and (b) alone say nothing about accuracy); (d) `ki3(60.0) == 0.0`. Its docstring states what it proves (the two copies agree and the copy is accurate to 3e-6 against the integral) and does NOT prove (independence of the approximation: two copies of one 2.6e-6 approximation drift together, so #641 item 2 stays open). **First red**: `ModuleNotFoundError` for `orpheus.numerics.bickley` before the module exists. **Mutation** `[M]` `m_ins3.py`: a degree-62 interpolant (monkeypatch `_KI3_DEG` to 62 and clear the `lru_cache` of `ki3_scaled_chebyshev`) differs at 2003 of 2004 grid points (legs (a) and (b) redden; leg (c) reddens when the mutant is accurate to less than 3e-6, which the builder reports). Per-panel exactness, the kinked integrand, belongs to G3.

**G5b. The binding of the cylinder kernel** (`tests/gates/cp/test_cp_kernel_is_the_numerics_ki3.py`; a `tests/gates/derivations/` file may not import `orpheus.cp`, an L3 package, so this leg cannot live in G5's file): the cylinder `CPMesh` binds `orpheus.numerics.bickley.ki3` (`mesh._kernel is ki3`) and `mesh._kernel_zero` equals `ki3(np.zeros(n_y))` bit for bit; the sphere binds the exponential kernel. **First red**: `AttributeError`/`ImportError` for `orpheus.numerics.bickley` before the module exists. **Mutation**, in-process: `monkeypatch.setattr(orpheus.cp.solver, "ki3", lambda tau: np.exp(-tau))` before building the cylinder `CPMesh` (the production module reads the name at call time): red.

**G6. The re-baselined carry rows** (`tests/gates/cp/test_cp_boundary_carry.py`, `test_the_carry_of_the_cp_solves`). The three rows `[M]` `m_ins7.py` (in-process patch of `orpheus.cp.solver.composite_gauss_legendre` onto the production rule; the file goes 3 failed, 5 passed): `cylinder-vacuum` (fixture `keff_cp/cylinder/vacuum` 0.970601242717418, new 0.9706012427174181), `cylinder-white_1.0` (fixture `keff_cp/cylinder/white_1.0` 0.7349098770721597, new 0.7349098770721564), `sphere-vacuum` (fixture `keff_cp/sphere/vacuum` 0.6319877457593569, new 0.6319877457593566); the fixture values were read from `tests/gates/transport/carry_head.npz` (`np.load`, `repr(float(...))`), the new values are the ones the patched run produced (`m_ins7.py`); the builder re-captures all three with the landed code under `python -O` and requires them equal to these (the builder stops and returns the difference in NEEDS if they differ); `sphere-white_1.0` and `test_the_default_pin_cell_outer_face_is_bound_to_white` do not move. The moves are 1, 30 and 3 ulp. Edit: a module constant `_REBASELINED_BY_THE_RULE: dict[tuple, float]` keyed by `(chart, declaration)` holding the three new values with a comment giving the reason (the production Gauss-Legendre rule is correctly rounded, the reference's `numpy.polynomial.legendre.leggauss` is the platform eigensolver output; the k moves by 1 to 30 ulp; #641) and naming the trigger to delete it: a recapture of `carry_head.npz`; `expected = _REBASELINED_BY_THE_RULE.get((chart, declaration), H.fixture_k("cp", chart, declaration))`. A new leg asserts that each override DIFFERS from `H.fixture_k(...)` for its row, so a recapture of the fixture (which would make the override a stale second source of the expected k) reddens it. The carry test's docstring ("equals the HEAD `k` bit for bit", `test_cp_boundary_carry.py` lines 88 to 95 at base) and its failure message ("the HEAD k") are reworded for those three rows: they say "the HEAD k, or, for the three rows of the override, the k of the correctly rounded production rule". The shared fixture `carry_head.npz` is NOT edited (decision 2.6.1: 14 test files read it, and `test_the_fixture_holds_the_cp_rows` pins `meta/commit == HEAD_COMMIT`). The commit message gives the reason and the old and new values. **Check the pin is not blind** (D.17): with the three values restored to the old ones in a copy of the test (a temporary file copied from the real one and edited; never edit the real file for the mutation), the three rows must redden (they do: that is how the values were found); and the pins that do not move (`sphere-white_1.0`, the pin-cell test) remain bit-exact, which is the carry's evidence that the chord and the exponential kernel did not change. **First red** of the override leg: a test-only input (an override set equal to the fixture value reads red).

**By-design reds of the unit** (existing gates that move with it): `tests/gates/cp/test_cp_boundary_carry.py` 3 rows (G6); `tests/gates/numerics/test_concept_order.py::test_every_module_has_a_rank` for `bickley` until the tier entry lands. Nothing in `tests/` patches `orpheus.cp.solver.composite_gauss_legendre`, `_ki3_kernel` or `chord_half_lengths` by string or attribute `[M]` (`git grep -n -E "monkeypatch\.setattr\([^)]*cp\.solver|orpheus\.cp\.solver\.(composite|_ki3|chord)" -- tests` prints 0 lines; control: `tests/gates/cp/test_cp_first_collision.py` names `chord_half_lengths` in a `rests_on` test-id string, which stays valid because the derivations test stays). Outside `tests/gates/cp`, no gate moves under either swap `[M]` (99 passed in all three runs over the five files that call `solve_cp`: `cross_method/test_family_convergence_contract.py`, `mc/test_cross_verification.py`, `moc/test_verification.py`, `sn/verification/analytical/test_cp_standoff_curvilinear.py`, `derivations/test_sood_registry_compatibility.py::test_solve_cp_accepts_registry_sphere_case`). The suite `tests/gates/cp` runs 319 passed, 21 skipped, 6 xfailed at base (33 s). The 84 s run with the production rule uncached was taken with `_affine_remap` and `_correctly_rounded_rule` (both `lru_cache`d in `generating_measure.py`) active, so it was the first build of each distinct panel; the builder reports the time of `tests/gates/cp` after the unit (`[R]`, unmeasured here) and, if it exceeds the base by more than the first-build cost of the panels, files the cost as a follow-up (8.1) and does not add a cache.

### 2.5 Layers, bit-identity and re-baselines

- `orpheus/numerics/bickley.py`: layer `numerics` (L1), imports `functools`, `numpy`, `mpmath` only (the gate `FORBIDDEN_EDGES["numerics"]` forbids `data`, `mesh`, `specification`, `reference`, L2, L3, L0). `cp/solver.py` (L3) imports `orpheus.numerics.bickley`, `orpheus.numerics.generating_measure`, `orpheus.geometry.chord`, `orpheus.geometry.line`: allowed edges (L3 above numerics and geometry).
- Must stay bit-identical: the chord array (1024 of 1024), the exponential kernel and `sphere-white_1.0`, the whole cylinder kernel (the SAME interpolant), every other CP number except through the rule swap. The cylinder k values do not move through the kernel (a copy). The only moves: the three carry rows by 1, 30 and 3 ulp through the rule swap; the sphere-vacuum row moves through the rule although the sphere uses no Bickley kernel.
- `CPParams.n_ki_table` and `CPParams.ki_max` have 0 production readers (`[M]` structure check, `git grep`) and `tests/gates/cp/test_verification.py:1100` pins that they do nothing; after the unit `ki_max = 50.0` and `_KI3_TAU_MAX` are two homes of one value. Binding decision 3 retires only `reviewed.json` rows, so these two fields are NOT retired in W0 (unread rows, section 8 "unread" list); an issue is filed (8.1). The comment under `class CPParams` beginning `# n_ki_table, ki_max: retained for backwards compatibility` is rewritten for tense (2.1) and now says the fields are inert and the kernel's bound is `_KI3_TAU_MAX` of `orpheus.numerics.bickley`.
- Docs: `docs/theory/methods/collision_probability.rst` gets a paragraph (Cardinal Rule 3: derivation, rationale, evidence): the Bickley-Naylor function `Ki_3`, its tan-substitution integral, why the solver uses a degree-63 interpolant of the scaled kernel (the exponential decay becomes a slowly varying function), the measured accuracy (2.6e-6 relative at tau = 0 over the 81 points of 2.3), that no `verifies` marker is attached and #348 (the open issue for the `Ki` family's labelled equation) is the place to attach one, that production holds a COPY of the interpolant the reference `flat_source_cp` holds (so a defect in the shared approximation is invisible to the cross-check, #641 item 2 open), the correctly rounded composite rule, and the re-baseline of three carry rows. `docs/architecture/layering.rst`: the own-table edge and the cold-import leg.

### 2.6 Decisions inside this unit (each also listed in section 8)

1. **Where the carry re-baseline lives.** The test override (3 rows, with the reason, a leg asserting each override differs from the fixture value), not the shared binary fixture. Structural reason: the fixture is a capture of a named commit (`meta/commit == carry.HEAD_COMMIT`, asserted by `test_the_fixture_holds_the_cp_rows`) read by 14 test files; editing it would make the capture lie about its commit; the override touches 0 binary files and keeps the capture a faithful record.
2. **The root fix of the rounding certificate is a unit of its own, `quadrature-ties` (2.0; ruling 8.44).** This unit makes no parity refusal, carries no expected failure and defers nothing to an issue; the swap to `LEGENDRE.on(a, b).gauss(n)` is the consumer of that fix and integrates after it (7.1). The one open point is 8.47.
3. **The Ki3 copy rebuilds the interpolant from mpmath, it does not store 64 coefficients.** Structural reason: the coefficients are a derived value of one integral, and a stored table is an untraceable constant; the witness G5(b) compares the rebuilt coefficients with the reference's bit for bit. Cost: about 0.3 s once, lazily. The alternative (stored `repr` constants round-trip bit-exactly and remove the build) is refused because it adds a second source of truth to a unit whose purpose is a twin held equal by a witness.
4. **The composite rule stays a private helper in `cp`, tagged as debt.** Binding decision 2 keeps it inside `cp`; the missing primitive is a direct sum over abutting intervals in `DiscreteMeasure.__add__` (W8a), which is the trigger written in the tag; the issue is filed in 8.1 (orchestrator's decision I2).
5. **The hollow body.** A hollow `CPMesh` is refused upstream (`_dropped_face`, #513) for the tag-bearing inner face, so no `edges[0] != 0` refusal is added in `CPMesh`; the function-level equality is gated in G4 (2.2 item 1), and that leg discharges decision I4: a hollow `CPMesh` cannot be built today (#513), so no k exists to capture.
6. **`#641`.** Closed by hand when production's imports of `orpheus.derivations` are 0 (done-when 2) and the edge gate is green (done-when 1); item 2 of #641 (an independent verification chain) is a filed issue (an accurate Bickley evaluator, 8.1), and the issue's closing comment says what the unit did and did not do (orchestrator's decision I5).

### 2.7 Retirement-audit items for this unit

A.1/A.2: the three symbols' callers (`nexus callers chord_half_lengths` etc. agree with the AST census: 5, 8 and 17 callers; production 1 file each, derivations keeps its own, tests keep theirs); B.4: `.venv/bin/nexus dead-references` after the docs rebuild; B.8: the four prose places and the paragraph that justified sharing the interpolant; D.14: no gate compares production's chord with itself through a wrapper (the chord witness says why it is tautological); D.15: `tests/gates/cp/test_cp_first_collision.py` gives the cylinder an absolute band of `4e-6 / tau` "derived from the interpolant" and `test_verification.py` compares `keff` at `< 1e-5` and `rtol=1e-6`: they are tolerant of the interpolant's error by design and an accurate production kernel (a later issue) would let the cylinder band tighten; nothing to edit now; D.16: the reference copies stay, so no comparison is single-sourced; D.17: G6's old-value mutation; D.18: the shape refusals the helper no longer makes are the partition's, pinned by G3's two refused legs; E.19: the prose rewrites; G.25: the unit moves nothing out of `derivations/`.

### 2.8 Done when (each predicate is answered by one command)

1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_layer_imports.py tests/gates/test_test_layers.py tests/gates/numerics/test_concept_order.py` exits 0, and the output names the new G1 and G2 tests.
2. The AST probe (a file under a quoted heredoc, run with `.venv/bin/python -O`) over the 260 production files reports 0 imports of `orpheus.derivations`, with a positive control reporting 3 on a copy of `cp/solver.py` at the base hash (`git show b68f8b8c:orpheus/cp/solver.py`).
3. `.venv/bin/python -O -c "import sys, orpheus.cp.solver; print(sum(m.startswith('orpheus.derivations') for m in sys.modules))"` prints 0; the control `import orpheus.derivations.common.kernels` prints a number above 0.
4. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/cp tests/gates/derivations/test_cp_ki3_copy.py` exits 0 (no xfail is involved: the odd-n leg is a plain test that passes once `quadrature-ties` has integrated).
5. `sphinx-build -E -W docs docs/_build/html` exits 0, and the prose grep of 2.1 (`git grep -n -i -E "bit-identical|shares? the (chebyshev|interpolant)|production cp|cp/solver|composite_gauss_legendre|chord_half_lengths|_ki3" -- docs orpheus/cp`) has been run on `git show b68f8b8c:orpheus/cp/solver.py` written to a temporary directory as a positive control (the four anchors of 2.1 are found there) and every present-tense hit that survives in `docs/` or `orpheus/cp` after the unit has been read and is history or a reference-side claim (0 present-tense claims that production CP shares the reference kernel). `nexus dead-references` is also run, but the three names stay alive in `derivations/`, so it cannot read positive for them and is not the instrument of this item.
6. `git grep -n "orpheus.derivations" -- orpheus/cp` prints only lines the unit left on purpose (none expected; each remaining line is read).
7. `[M]`-style check of the ledger tag: `git grep -n "ELEGANCE-DEBT\[guard\] #" -- orpheus/cp` prints the helper's tag line with a real issue number (control: `orpheus/transport/fields/_bases.py` carries the `#662` tag).

## 3. Unit `dead`: retire, move to the test side, keep

**Outcome.** 92 rows leave production (49 retired, 43 moved to the test side) in four sub-units of disjoint code; 45 rows stay, each with the wave or consumer that keeps it (the 41 `keep` verdicts of `dead/reviewed.json`, and rows 134 to 137, which stay as catalogue data for W8a: 3.0 item 2); 4 further rows that `dead/reviewed.json` rules `move` (52, 99, 106, 116) are taken out of W0's moves because their body reads private state that no public accessor exposes (3.1, decision D2; section 8.3). The unit's rows are the READ verdicts of `scratch/ontology/w0/dead/reviewed.json` (146 rows: 58 retire, 47 move, 41 keep; 34 of the census's 146 actions were overturned by that read). Only a symbol that file RULES is retired or moved: a closure member it does not rule is not touched in W0 and is listed in section 8.3 as unread, for a later reviewed read (decision G9). The five role rows (64, 65, 66, 144, 145: `BulkOperator`, `FullOperator`, `BoundaryOperator`, `stamp_boundary_role`, `_BlockRoleMeta`) belong to unit `roles` (section 4); the other 100 are cut into four sub-units, of which 92 act in W0 (the 4 deferred moves and the 4 kept retires, rows 134 to 137, do not). The unread census rows (113 classed by rule, 96 rejected) are NOT in this specification, by binding decision 3.

### 3.0 Premise checks this specification ran on the review (each changes a builder's action)

1. `[M]` For every retire or move row with a production reference (38 of 105 rows have at least one), each reference is an import, an `__all__` string, a co-retired or co-moved caller in the same sub-unit, or a rewrite listed in the appendix table. The references that are neither an import nor a string: rows 44/45 (`evaluate_from_components`, each calls the other: both retire together, dead-B); rows 130 and 131 (`dominant_eigenpair` calls `_sign_normalised`, `direct_eigenvalue` calls `dominant_eigenpair`, `rayleigh_quotient_iteration` calls `_sign_normalised`: one cluster in `eigenvalue.py`, dead-B); rows 132 and 133 and the closure of row 73 (all inside `select_quadrature`, dead-A); rows 123, 124 and 139 (`foldable_part` and `residual_part` call `_sibling`, dead-D); rows 111 and 140 (`HarmonicMomentFlux.truncate` calls `truncated_head`, dead-D: one commit); rows 85, 86, 141, 142, 143 (the `assembly.py` module, dead-C); rows 144 and 145 (roles). The one reference that reads like a live production use and is not one: row 138, `TransferOperator.foldable_sigma` (retired) against `orpheus/sn/acceleration/dsa.py` `DSALowOrderSystem.from_problem`, which calls `mat_xs.foldable_sigma()` on a `MaterialXSField` (`orpheus/transport/mesh/material_xs_field.py:693`, KEPT): a homonym `[M]` (read at the quoted line `fold = mat_xs.foldable_sigma()`). 0 of 105 rows have a live production consumer outside the sets above.
2. **Rows 134 to 137 (`_gl1d_invert`, `_lebedev_invert`, `_ls_sn_invert`, `_product_invert`) leave W0: the four helpers and the field `QuadratureSpec.degree_of_exactness_for` stay as catalogue data, kept for W8a.** `[M]` `git grep -n degree_of_exactness_for -- orpheus`: the field is declared at `registry.py:545` as a REQUIRED field (`Callable[[int], dict[str, Any] | None]`, no default), assigned the four helpers at `registry.py:821, 828, 835, 842`, and read in production only at `registry.py:1384` inside `select_quadrature`; `git grep -c degree_of_exactness_for`: `registry.py` 8, `tests/gates/numerics/test_registry.py` 10, `test_registry_rule_list.py` 3, `test_invariance_verdicts.py` 1, `tests/gates/numerics/data/generate_s5_invariance_baseline.py` 1; docs 0. The four helpers are the values of a required field of a kept catalogue (G1-a), so they are catalogue data and not dead code, although the selector, the field's only production reader, retires. They go to the KEPT table (3.3) with the reason "a value of the required field `degree_of_exactness_for` of the kept catalogue; W8a", and the dead-A builder does not touch them, the field, `_lebedev_orders`, or the nested helpers `n_min`, `_even_at_least`, `ceiling`, `n_mu_min` and `n_phi_min`. No test line that calls the field is rewired and no literal parameter dict is written: a test that calls the field outside a test whose subject is the selector keeps calling it (`test_registry_rule_list.py::test_the_degree_inversion_is_the_cheapest_parameter_set_meeting_the_target`, `test_invariance_verdicts.py`, `generate_s5_invariance_baseline.py` and the `QuadratureSpec(..., degree_of_exactness_for=...)` builds of `test_registry.py` stay as they are), and the lines inside a test whose subject is the selector go with that test (3.6). The degree-to-parameter map of the four rules keeps its tests.
3. **`QuadratureSpec.structural_flags` retires with `select_quadrature` and ERR-118.** It is the subject of ERR-118 and has no consumer after the selector goes (`git grep -n structural_flags`: `registry.py` lines 536, 549 and 1481 (the class docstring, the method, the selector's call), `test_registry.py:127` and `:1825`, and one `:meth:` role in `docs/theory/verification/error_catalog.rst` at about line 11612, in ERR-118's `**Status:**` paragraph). The review's note on row 73 records the ERR as "retired with its code" and decision D4 adopts that record; the review rules the flag with row 73, so it is a ruled retirement and not an unruled closure member (the one exception to 3.1 rule 3). It is an ADDITIONAL `Retired` entry of the dead-A rows module (D-A.7, not counted among the rows). The ERR-118 witness is rewired as 3.1 gate D3 says, and the lines that remain to edit are: the class docstring sentence "`structural_flags` reads them off the rule", the test block "Structural-flags accessor must surface the four flags" of `test_registry.py` (about lines 126 to 133, the `flags = spec.structural_flags(...)` assertion and its set; the line `assert callable(spec.degree_of_exactness_for)` stays), and the `:meth:` role in ERR-118 (past tense). The two spec fields only it read (`axis_aligned` and `half_range_clean`) stay as catalogue data (they are rows of the kept catalogue G1-a); they become data read by tests only; that is recorded for the registry unit, not acted on.
4. **Row 104 `RealizerTable.catalogue`** is KEPT, for the method tables' admission matrix, asserted through it by 13 test files; nothing in W0 opens the realizer-table cycle (3.4, decision D3).
5. **Rows 54 and 55** (`HERMITE`, `laguerre`, move to the test side) and unit `defects` #711: `laguerre` and `HERMITE` are `GeneratingMeasure` factories whose construction runs the #711 refusal; the refusal is specified over the factories that REMAIN in production (`LEGENDRE`, `CHEBYSHEV_T`, `CHEBYSHEV_U`, `jacobi`, `_affine_remap`); a moved factory builds its measures through the same `GeneratingMeasure(...)` constructor and so passes through the refusal. Dead-A integrates before defects (section 7).

### 3.1 The sub-units and what each retires or moves

The cut is by package and closure so that a retired symbol and the helpers that exist only for it stay together (G.24). `[M]` `fs_final.py`: the pairwise intersections of the four sub-units' hard production-and-test sets are empty in 6 of 6 pairs. This is conditional on the suspect tier: with the 73 unresolved suspect (file, method-name) pairs counted, 11 distinct test files are shared (A and B 2, A and C 1, A and D 4, B and D 1, C and D 3), each reached only by a method-name hit (`.summary`, `.issues`, `.upstream`, `.discriminator`, `.dyad`, `.n2n_matrix`); the disjointness holds or fails with those files (section 7 serialises the sub-units, so a shared file is a textual rebase, never a lost edit). Docs pages are NOT cut file-disjointly (4 to 7 pages per pair, among them `error_catalog.rst`, `docs/api/numerics.rst`, `discrete_measures.rst`); each sub-unit edits the docs LINES of its own rows, and the orchestrator's docs pass reconciles.

| sub-unit | rows | retire | move | hard prod | hard test | hard docs (the tense table is the authority; the docs list is a subset) | suspect tests |
|---|---|---|---|---|---|---|---|
| dead-A: measure, quadrature selector, geometry, descent | 24 | 13 | 11 | 16 | 33 | 24 | 25 |
| dead-B: numerics rest and eigenvalue | 19 | 12 | 7 | 11 | 17 | 15 | 1 |
| dead-C: sn | 22 | 5 | 17 | 13 | 23 | 17 | 8 |
| dead-D: transport | 27 | 19 | 8 | 17 | 21 | 11 | 27 |

The rows, retire and move columns are the counts AFTER decision D2 (rows 52, 99, 106, 116 leave the moves: 3.1 below), after row 126 is given one home (dead-D; section 8.3) and after rows 134 to 137 are kept (3.0 item 2): 24 + 19 + 22 + 27 = 92 rows, 13 + 12 + 5 + 19 = 49 retire and 11 + 7 + 17 + 8 = 43 move (the dead-A column lost rows 134 to 137 to 3.0 item 2; command `.venv/bin/python scratch/ontology/w0/spec/amend/count_rows.py .claude/plans/ontology_w0_spec.md`, which counts the verdict column of the D-x.1 tables, minus the four deferred rows) The three "hard" and "suspect" columns were measured before those three changes and are carried as measured (`[M]` `fs_final.py`); the deferred rows' files are simply not edited for those rows, so the columns can only shrink. **Not read, stated plainly (decision D6):** the 73 suspect (file, method-name) pairs in 66 test files and the tense tables of the appendices are read by the builder, not by this specification; the pairwise file-disjointness of the four sub-units is conditional on those unresolved pairs, and the build order dead-A to dead-D is serial (section 7), so a test file shared by two sub-units is a textual rebase and never a lost edit (listed in 8.3).

Each sub-unit's appendix gives: every row with its file, symbol, verdict, destination or rewire, the review's test fate, markers, production references with their kind, closure and the test files naming it (Appendix D-A to D-D); its hard production, test and docs files; its suspect test files; and the tense table of every prose line naming one of its retired names (`PRESENT` is rewritten, `role` is a docstring role that is dead once the name is gone and is dropped or re-pointed, `eq/ref` is a `:by:` target or label to re-point, `past` stays, `collision` is not the symbol). The tense population `[M]` (`measure_tense.md`): 105 distinct bare names; corpus 1545 files (tracked `docs/` without `_build`, `.claude/agents`, `.claude/skills`, `.claude/rules`, `CLAUDE.md`, and the AST docstrings of tracked `orpheus/` and `tests/` `.py`); word-bounded Python `re`; 2343 raw hits, 826 after dropping ordinary-word uses of 28 generic names; of the 826: 377 dead roles (docs 178, orpheus 116, tests 83), 342 present-tense claims to fix (tests 157, docs 123, orpheus 60, `.claude` 2), 69 past (stay), 14 `:by:`/label targets, 24 lexical collisions; positive control 105 of 105 names match their own definition line; 18 names have 0 hits in docs and docstrings (zeros with that control). NOT searched: comments, non-docstring strings, `tools/`.

**Move rows** (43 act in W0). A moved symbol keeps its name on the test side. Three rules (decisions D1, D2, G9):

1. **One home per moved symbol.** A move with TWO OR MORE consuming test files listed in the appendix row's last column gets one module home, named below, under the gate directory of the package the symbol left; a move with ONE consuming file is inlined in that file. (The count is the appendix column; the builder first deletes the own-subject tests of the row's test fate, and if one consumer then survives the row keeps its module home: the row's `home` is fixed, so a rows module never changes shape. A function called by a sibling that has a module home goes with the sibling: `_MMHalfGrid.downstream_per_ordinate`, row 80.) Each row of 3.1's table below states which. Modules (all new, one per cluster):
   - dead-A: `tests/gates/geometry/_manifold_fixtures.py` (row 27 `UNIT_INTERVAL`, 6 consumers); `tests/gates/geometry/_geometry_oracles.py` (rows 6 `tangencies`, 32 `_group_elements`, 33 `stabilizer_order`, 36 `permutes`: 2, 2, 2 and 4 consumers); `tests/gates/numerics/_measure_fixtures.py` (rows 54 `HERMITE` 3, 55 `laguerre` 2, 56 `ENERGY` 4, 60 `equispaced` 7, 61 `from_srepr` 2); inline: row 34 `element_order` in `tests/gates/geometry/test_transformation.py` (1 consumer).
   - dead-B: `tests/gates/numerics/_dense_eigen_oracle.py` (rows 49 `direct_eigenvalue` 3, 130 `dominant_eigenpair` 2, 131 `_sign_normalised` 0, the dependency of the other two); `tests/gates/numerics/_overlap_basis_oracles.py` (row 43 `fractional_columns`, 2); `tests/gates/numerics/_space_oracles.py` (row 74 `project_live`, 2); inline: row 42 `dominant_column` in `tests/gates/data/test_energy_grid.py`, row 46 `multi_index` in `tests/gates/numerics/test_interval_legendre_basis.py` (1 consumer each).
   - dead-C: `tests/gates/sn/sweep/curvilinear/_psi_half.py` (rows 80, 81, 82); `tests/gates/sn/sweep/_ordinate_assembly.py` (the whole of `orpheus/sn/loss_representation/assembly.py`: rows 85, 86, 141, 142, 143; the module is deleted from production; its consumers `sn/sweep/test_assembly_mode.py`, `sn/sweep/core/test_multi_d_reverse_walk.py` and `tests/gates/derivations/test_dsa_production_tie.py`); `tests/gates/sn/sweep/_sweep_graph_faces.py` (rows 87 to 90, the four `face_*` accessors, 3 consumers each); `tests/gates/sn/operators/_kernel_dimension.py` (row 91 `predicted_kernel_dimension`, 3 consumers, all in `sn/operators`); `tests/gates/sn/architecture/_splitting_oracles.py` (row 98 `law_residual`, 2); inline: row 37 `n_unknowns_flat` and row 92 `SNProblem.n_unknowns_flat` in `tests/gates/sn/primitives/test_axis_primitive.py`, row 95 `boundary_vs_interior_split` in `tests/gates/sn/operators/test_typed_residual_evaluation.py`.
   - dead-D: `tests/gates/transport/_transfer_oracles.py` (rows 119 `dense_per_material` 3, 121 `on_basis` 7, 122 `on_frame` 2); `tests/gates/transport/fields/_scalar_boundary_oracles.py` (row 113 `p1_boundary_scalar_flux`, 2: one transport and one diffusion consumer); inline: row 105 `face_views` in `tests/gates/transport/fields/test_angular_boundary_flux.py`, rows 117 and 118 in `tests/gates/transport/test_kernels.py`, row 129 `kernel_rhs` in `tests/gates/transport/spatial/test_ld_ubld_primitive.py`.
2. **A moved method becomes a free function of its former receiver**, in its home, taking the object as the FIRST argument (a classmethod takes the class as first argument; a property becomes a function of the object), with every call site rewritten `x.m(a)` to `m(x, a)` and `x.p` to `p(x)`. The signature is in the rows module beside the row. Measured, the 43 moved rows read only public attributes of their receiver, with four exceptions that `[M]` (`scratch/ontology/w0/spec/amend/moves.py` and `priv.py`: an AST walk of each moved body for `self._x` and `cls._x`, then a read of each hit's owner) read private state that no public accessor exposes, and which are therefore TAKEN OUT of W0's moves (they stay in production, untouched, are listed in `KEPT` of their sub-unit with the reason `deferred, private state`, and are in 8.3 with the attribute names): row 52 `Field.where_largest` (calls `self._index_tuples`, a private method also read by the kept `cone_violations`), row 99 `CollisionCache.reset_build_count` (writes the private class variable `_build_count`, also incremented by `from_geometry`), row 106 `BoundaryField.from_face_arrays` (calls the private hook `cls._face_space_of`, which every leaf overrides), row 116 `MaterialXSField.n2n_matrix` (reads `self._n2n_dense` and calls `self._build_dense_caches`). Limit of the probe, stated: it sees `self` and `cls` receivers only; a private read through another receiver is the builder's first red.
3. **A closure member is not touched unless the review rules it (decision G9).** The closure column lists production definitions referenced in the retired body whose every production reference lies in that body. A member that is itself a row of `reviewed.json` is retired or moved with its row; a member that is not a row stays in production, even orphaned, and is listed in 8.3, with one exception: `QuadratureSpec.structural_flags`, which the review rules with row 73 through its ERR-118 note (3.0 item 3). Where a moved body needs such a member it imports it from production: the moved `laguerre` imports `_laguerre`, the moved `element_order` imports `_MAX_ELEMENT_ORDER`, the moved `predicted_kernel_dimension` imports `_anova_dimension`. `DiscreteMeasure.degree_of_exactness` and `AngularSymmetry.admits_symmetry` are KEPT (a native invariant of the measure and a leg of the ERR-081 rewire).

**The layer rule for a moved symbol: its test-side home is in the gate directory named for the package the symbol LEFT** (`numerics` symbols in `tests/gates/numerics/`, `geometry` symbols in `tests/gates/geometry/`, `sn` symbols in `tests/gates/sn/...`), and the home imports exactly what the symbol's body imported. Reason: the test-layer scan judges a consumer by the layers it reaches; a consumer that imported `orpheus.<pkg>` for the moved name now imports `tests.gates.<pkg>...`, the same layer, so its (file, package) edge set is unchanged and `OPEN_TEST_LAYER_VIOLATIONS` neither gains nor loses a row from the move. `[M]` two cases where a lowest-consumer-layer rule would have failed and this one does not: `geometry/test_manifold.py` already imports `orpheus.numerics` (`ENERGY`; a listed open violation) and keeps the same edge through `tests/gates/numerics/_measure_fixtures.py`; `derivations/test_dsa_production_tie.py` already imports `orpheus.sn` and `orpheus.transport` (both listed in `OPEN_TEST_LAYER_VIOLATIONS`) and keeps the same edge through `tests/gates/sn/sweep/_ordinate_assembly.py`; that test file is NOT moved. Simulate: `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_test_layers.py` must report 0 new and 0 stale rows after each sub-unit (a stale row means a consumer stopped violating: retire its row from `OPEN_TEST_LAYER_VIOLATIONS`, which makes `tests/gates/test_test_layers.py` part of that sub-unit's set; the first red decides). Disjointness: each sub-unit routes its helpers to its OWN new module or inline; two sub-units choosing one existing helper module (`tests/gates/sn/_test_helpers.py`, `tests/gates/diffusion/_wave2_helpers.py`) would collide. Do not use `tests/conftest.py` as a destination (every test reads it).

**Common per-sub-unit procedure** (all four):
1. Build the table of the sub-unit's rows from its appendix into `tests/gates/_retired/rows_w0_dead_<x>.py` (the data of the shared gate, 1.5; the rows module of Appendix D-<X>.7 is the starting text, in the three row types, with each `Moved` home concrete and the sub-unit's `KEPT` entries), then run `tests/gates/test_retirements.py`. The gate is red first (every row is defined today): that is the first red of the sub-unit.
2. Retire leaves first, in ONE commit per sub-unit, whose body carries the dependency table (the closure column) (G.24): the row's callers inside the set, then the helpers only the set called, then the top-level symbol; the `__init__`/`__all__` re-exports named in the production-reference column go in the same commit. A closure member that no row rules is not in the commit (3.1 rule 3).
3. For each test fate: rewire = keep the behavioural law and point it at the successor named in the "destination or rewire" column; delete = an API-smoke test of the retired object; keep markers (`catches`, `verifies`) by moving them with the behavioural test (C.13); never leave a `verifies` label whose only test is deleted (its equation turns orphan in `tests/_harness/audit.py`).
4. Edit the docs lines of the tense table: rewrite `PRESENT`, drop or re-point `role`, re-point `eq/ref` to the moved home, leave `past`. Re-run `.venv/bin/nexus dead-references` after a docs rebuild.
5. Run the sub-unit's test files (appendix list) one at a time, serially, plus `tests/gates/test_retirements.py`, `tests/gates/test_layer_imports.py`, `tests/gates/test_test_layers.py`, `tests/gates/test_docstring_xrefs.py`, `tests/gates/test_error_catalogue_reconciles.py`, `tests/gates/test_vv_harness_audit.py`. The suspect files (appendix) are opened only when a retired method name appears in them AND the class is the retired owner; a hit on a homonym is left alone (the builder's first red decides, never a name grep).
6. The string-keyed consumers that stay green when the name is gone (silent residue; B.6) and that the sub-unit must edit, by file (the quoted symbol decides, the line number is as of base): `tests/gates/mesh/test_partition.py` (dict key `compute_volumes_1d`, two places), `tests/gates/geometry/test_coord.py` (a set, `compute_volumes_1d`), `tests/gates/sn/acceleration/test_dsa_rate.py` (the allow-list set holding `foldable_sigma`), `tests/gates/sn/solve/test_si_diagnostic_trajectory.py` (`where_largest`, two places), `tests/gates/homogeneous/test_homogeneous.py` (a tuple holding `dominant_eigenpair`), `tests/gates/sn/primitives/test_method_space.py` (`minimal`), `tests/gates/transport/frames/test_harmonic_frame.py` (`from_measure`), `tests/gates/transport/frames/test_binding_tightness.py` (`equispaced`), `tests/gates/numerics/test_generating_measure.py` (`laguerre`), `tests/gates/cross_method/test_one_parse_on_production.py` (the dotted strings `"SNSolver.compute_group_absorption_rate"`, two places), `tests/gates/test_layer_imports.py` (three string fixtures that import `dominant_eigenpair`, in the refused legs of the reference-insulation block, the quoted lines naming `from orpheus.numerics.eigenvalue import dominant_eigenpair`), and the production `__all__`/re-export strings: `orpheus/geometry/__init__.py` (`compute_volumes_1d`), `orpheus/mesh/__init__.py` (`n_unknowns_flat`), `orpheus/numerics/quadrature/__init__.py` (`select_quadrature`, `SelectionLog`, `QuadratureSelectionError`), `orpheus/sn/__init__.py` (`build_coupled_system`), `orpheus/numerics/spaces/__init__.py` (`truncated_head`, two places). `[M]` `m_dead6.py` (the AST parent of each string) and `m_dead8.py` (word-bounded match inside non-docstring string constants shorter than 200 characters, production, tests, tools: 59 of 103 bare names have such a hit, the real residue being the sites above and many homonyms).
7. Every line number the appendices give is as of the base hash; an earlier unit may have moved it. The symbol and the quoted anchor decide.

**Gate D3, the error catalogue** (dead-A; the decision for ERR-081 and ERR-118; `tests/gates/test_error_catalogue_reconciles.py`, a CI-run file `[M]` listed in `.github/workflows/gates.yml`). `[M]` `m_dead3.py`: 126 declared entries, 125 marker ids, `uncaught = catalogue - markers - open` empty today; ERR-081 has exactly 1 `catches` marker (the class-level marker on `TestStageZeroIsTheDescentArrowPlusTheUnspentSymmetry`, `tests/gates/numerics/test_registry.py`) and ERR-118 exactly 1 (`test_a_rule_selected_for_positive_weights_has_positive_weights`, which calls `select_quadrature`). Deleting ERR-118's only marker makes arm 2 (`test_every_declared_entry_has_a_catching_test`) red with `uncaught == ['ERR-118']` (simulated in-process; control: removing ERR-075's six markers gives `['ERR-075']`). The gate has NO state for "retired with its code": `_open_ids` accepts only `**Status:** ... OPEN` with an issue and no marker. So the unit specifies the missing state (section 8, decision 4, lists the alternative):
- Catalogue: an entry whose catcher and subject are gone carries a paragraph `**Retired:** with its code: <name>[, <name>...]` (one paragraph; the entries are the `file:qualname` pairs of the retired definitions), added to ERR-118: `**Retired:** with its code: orpheus/numerics/quadrature/registry.py:select_quadrature, orpheus/numerics/quadrature/registry.py:QuadratureSpec.structural_flags` (the pairs are `file:qualname`, never bare names), followed by one sentence of the date (2026-11-20) and the reason (the selector had no production caller; the typed-flag defect it caught is unspellable without the selector and its flags accessor).
- Gate: `_RETIRED_RE = re.compile(r"^\s*\*\*Retired:\*\* with its code: (.+?)\s*$", re.M)`; `_retired_ids()` returns the ids whose entry body matches AND for which none of the named `file:qualname` pairs is DEFINED in its file (the ONE definition finder of 1.5, `tests/gates/_retired/_finder.py`, imported here: a second finder would be a second spelling of one predicate; and a bare-name search under `orpheus/` would let a later unrelated homonym of `summary` or `select_quadrature` silently un-retire ERR-118) AND which carry no `catches` marker; arm 2 becomes `uncaught = _catalogue_ids() - _marker_ids() - _open_ids() - _retired_ids()`; a new test `test_a_retired_entry_names_code_that_is_gone_and_has_no_catcher` with a positive control in its own body (a synthetic entry naming a defined function is NOT retired; one naming an undefined function is retired; one carrying both a retired line and a marker is stale) and a live assertion that `"ERR-118" in _retired_ids()`.
- **First red** of the new test: today's tree, where `select_quadrature` is defined (the gate reads ERR-118 as not retired, so the test that asserts `"ERR-118" in _retired_ids()` is red until the unit lands). **Mutation**: the retired-ids function takes the repository root as a parameter; copy `orpheus/numerics/quadrature/registry.py` to a temporary root at the same relative path, re-add `def select_quadrature` to the copy, call the function on that root, and require ERR-118 to leave the retired set.
- Arm 5/6 (`test_every_test_the_catalogue_cites_exists`): ERR-118's entry cites the deleted test by path; rewrite the citation to the module-and-directory form the gate's docstring prescribes (read lines 160 to 290 of the gate: a retired test is named by module and directory, never by path). Of 272 test citations in the catalogue, 1 names a test slated for deletion by the `dead` and `roles` units (this one).
- Other readers of the catalogue, read (`[M]` `tools/verification/generate_error_index.py`, `tests/_harness/audit.py`): the index generator lists every entry with no catcher under `## ⛔ Uncaught` with no notion of status (ERR-103 is listed there today), so ERR-118 will appear there after the unit; that is accurate (no test claims it) and needs no tool change; `audit.py` computes `dormant_errors` from catchers only and does not read the new paragraph. So no `tools/` file joins dead-A. Arm 1 stays green (the ERR-081 marker moves to another test, below). `test_the_generated_index_matches_the_corpus` compares id sets only and stays green; `.claude/skills/vv-principles/error_index.md` changes (ERR-081 and ERR-118 catcher counts) and is regenerated at integration (section 7) with `.venv/bin/python -m tools.verification.generate_error_index`.
- **ERR-081 migration.** The class-level marker stays on `TestStageZeroIsTheDescentArrowPlusTheUnspentSymmetry`. Its legs, bar one, already test `AngularSymmetry.admits_domain`, `domain_refusal` and `quotient_onto` (kept, G1-b) and do not call the selector. The one leg that calls it, `test_the_rejection_message_names_the_one_failing_clause`, is rewired onto the kept object: `GEOMETRY_ANGULAR_SYMMETRY["cylinder"].domain_refusal(Quadrature.gauss_legendre(8).measure)` must contain `S^2/O2_x`, `fold by O2_x` and `unspent D_1h` and must not contain `arrow`; `GEOMETRY_ANGULAR_SYMMETRY["slab"].domain_refusal(Quadrature.lebedev(5).measure)` must contain `no descent arrow` and must not contain `unspent` (the same fragments the leg asserts today against `log.rejected`, minus the selector's `domain mismatch` prefix, which belongs to the selector). Check D.14: mutate `domain_refusal`'s reason in-process and require the leg red.
- **ERR-118's witness `test_lebedev_order_13_is_not_positive_weighted`** (`tests/gates/numerics/test_registry.py`) asserts a fact about `structural_flags`, which retires with the selector (3.0 item 3); rewire it to the fact itself, on the built rule: `assert not np.all(Quadrature.lebedev(order=13).weights > 0.0)` (3 of 32 Lebedev orders, 13, 25 and 27, have a non-positive weight `[M]` ERR-118's entry), carrying NO `catches` marker (ERR-118 is retired).
- `docs/theory/verification/error_catalog.rst`: in ERR-118 the line `No production path calls ``select_quadrature``, so no solve met the rule.` (the entry's `**How it hid.**` list) is a PRESENT-tense claim that becomes false (the function does not exist); rewrite to the past ("had no production caller; retired with its code on 2026-11-20"); the `:meth:` role to `QuadratureSpec.structural_flags` in the `**Status:**` paragraph (about line 11612: "reads them off the rule") is a dead role once the method is gone: rewrite to the past ("read them off the rule; retired with the selector on 2026-11-20"); and the `Catching tests` paragraph per arm 5/6. ERR-081's entry text that mentions the marker `catches("ERR-081")` as a "REQUIRED companion edit" stays (history of the entry, past).

### 3.2 Gates of the unit (summary; each sub-unit's appendix lists its rows)

- **The definition leg, the reference leg and the keep leg** (1.5), per sub-unit, over the sub-unit's rows module. First red: every `Retired` and `Moved` row of the sub-unit is defined today (`[M]` `m_dead1.py`: 105 of 105 retire/move rows of the file, rows 134 to 137 among them; the sub-unit's share after D2 and 3.0 item 2). Mutations (the root-path mechanism of G6) in 1.5. The keep leg's first red is a test-only input (it is green today).
- **D3** (above), dead-A.
- **D4, no silent skip.** For dead-B: after the unit, `git grep -n "rayleigh_quotient_iteration" -- orpheus tests tools examples docs ':(exclude)tests/gates/_retired' ':(exclude)docs/theory/verification/error_catalog.rst'` returns 0 lines (the two exclusions are the rows module, which names the symbol by design, and the catalogue history; every other line is a failure; control: the same command without the exclusions at the base hash prints more than 0), and `TestRayleighQuotientIteration` (6 tests, `tests/gates/numerics/test_eigenvalue.py`, bound by `getattr(_eig, "rayleigh_quotient_iteration", None)`) is gone. `[M]` `m_dead2.py`: with the symbol present, 6 passed; with the attribute deleted in-process, 6 SKIPPED, 0 failed. So no gate reddens when the symbol alone is deleted. What forces the deletion of the six tests is the reference leg's string rule: the `getattr` string whose first argument resolves to `orpheus.numerics.eigenvalue` is a reference to the retired definition (the definition leg is red only while the definition exists, and forces no test deletion). The two `monkeypatch.setattr(_eig, "dominant_eigenpair", ...)` strings in `test_eigenvalue.py` raise `AttributeError` under pytest's default `raising=True` when `dominant_eigenpair` leaves the module, so they redden by design and are edited in dead-B.
- **D5, docstring cross-references.** `tests/gates/test_docstring_xrefs.py::test_no_docstring_names_a_missing_object[orpheus|tests|docs]` resolves every `:func:`/`:class:`/`:meth:`/`:attr:`/`:mod:` docstring role by importing its target (`tools/check_docstring_xrefs.py`, one subprocess per tree, an empty allowlist) and reddens by design on every surviving role to a retired name: `[M]` approximately 34 of 105 rows are named by a docstring role, about 107 occurrences in 51 files (`m_dead9.py`; some sit inside the retired object's own docstring and go with it). It is in the suite, not in the CI list of `gates.yml`; each sub-unit runs it.
- **D6, no demoted gate** (D.14, D.16): the rewired rows whose successor is a delegate chain (rows 24 `compute_volumes_1d` to `coord.measure(edges)`, 93 `loss_kernel_gauge` to `pencil.kernel.projector`, 94 `boundary_face_layout` to `angular_trace.layout`, 84 `build_coupled_system` to `build_within_group_system`, 2 `minimal` to the constructor, 28 `barycentre` to `Quotient.lift`) were delegates or aliases: a test that compared the delegate with its target compared a value with itself (X4); delete such a comparison with the row and say so in the test's module docstring; keep every behavioural law. Input and mutation per rewired successor (D.14; each mutation is an in-process `monkeypatch` of the successor's body, and the rewired test must go red): row 24 (`coord.measure`), row 28 (`Quotient.lift`), row 93 (`pencil.kernel.projector`), row 94 (`angular_trace.layout`), row 84 (`build_within_group_system`), row 2 (the `DiffusionMethodSpace` constructor), row 67 (`RieszRaiseOperator`: a RE-TYPING, not a delegate: see 3.6 dead-B). The first red of D6 is a test-only input: no shipped input is rejected by a retired delegate.

### 3.3 The 45 rows kept, with the wave or consumer that keeps each

Gate D1' asserts each is still defined. The derivations rows (A2) are kept because `derivations/` calls them and is closed; the wave rows are kept because a later wave's target is exactly that object; a row kept as a homonym or a gate-pinned survivor says so.

| row | file : symbol | kept for (wave, concept row, or consumer) |
|---|---|---|
| 3 | `geometry/chart.py` `Chart.measure_density` | derivations/ (A2): characteristic/basis.py:285, closure.py:284 |
| 4 | `geometry/chart.py` `Chart.directions_at` | derivations/ (A2): characteristic/reading.py:176,202,276 |
| 5 | `geometry/chart.py` `Chart.line_domain` | derivations/ (A2): characteristic/assembly.py, lines.py, reading.py |
| 7 | `geometry/chart.py` `RadialImage.orbit_coordinate_at` | derivations/ (A2) via Chord/Rays.orbit_coordinate_at |
| 8 | `geometry/chart.py` `AxialImage.orbit_coordinate_at` | as row 7 |
| 9 | `geometry/chord.py` `ConcentricPartition.region_at` | W6, concept row 'Point location, chord' (#596: ConcentricPartition is the one locator) |
| 10 | `geometry/chord.py` `ConcentricPartition.refined_at` | derivations/ (A2): characteristic/reading.py:204,280 |
| 11 | `geometry/chord.py` `ConcentricPartition.coarsening` | derivations/ (A2): characteristic/basis.py:143, reading.py:206, transport.py:679 |
| 12 | `geometry/chord.py` `ConcentricPartition.rays_from` | derivations/ (A2): peierls geometry.py, fn sphere, characteristic |
| 14 | `geometry/chord.py` `Chord.orbit_coordinate_at` | derivations/ (A2) |
| 15 | `geometry/chord.py` `Chord.rays_at` | derivations/ (A2): characteristic/transport.py:453 |
| 16 | `geometry/chord.py` `Chord.rays_at_closest_approach` | derivations/ (A2): characteristic/transport.py:456 |
| 17 | `geometry/chord.py` `Chord.line_integral` | derivations/ (A2): 6 reads |
| 18 | `geometry/chord.py` `Chord.cosine_at` | derivations/ (A2): peierls geometry.py (6) |
| 19 | `geometry/chord.py` `Rays.lengths_to` | derivations/ (A2): peierls geometry.py:910,1181,1263 |
| 20 | `geometry/chord.py` `Rays.to_wall` | derivations/ (A2): 6 reads |
| 21 | `geometry/chord.py` `Rays.orbit_coordinate_at` | derivations/ (A2) |
| 22 | `geometry/chord.py` `Rays.advanced` | derivations/ (A2): characteristic/transport.py:459 |
| 23 | `geometry/coord.py` `CoordSystem.measure_density` | derivations/ (A2) through Chart.measure_density |
| 25 | `geometry/deck.py` `ISOMETRY_CLASSES_BY_DIMENSION` | W6 (declared deferral table, read with DEFERRED_DECKS) |
| 26 | `geometry/deck.py` `DEFERRED_DECKS` | W6: SCOPE-BOUNDARY table named by production's DeckError message (geometry/deck.py:293) |
| 29 | `geometry/structured_geometry.py` `StructuredGeometry.is_hollow` | derivations/ (A2): common/reference_body.py:164,188 |
| 35 | `geometry/transformation.py` `RigidMotion.seated_at` | W6, concept row 'Symmetry group, translations included' (a subgroup of E(d)) |
| 38 | `numerics/basis/base.py` `Basis.synthesize` | W8a, concept row 'Basis' (one table contraction) |
| 47 | `numerics/coupled_system.py` `CoupledSpace.system_slices` | W1, concept row 'Direct sum of spaces' (4 offset loops become one) |
| 48 | `numerics/coupled_system.py` `CoupledSpace.kernel_basis` | W1/W2: the metric-kernel family (concept 'Quotient by a metric kernel', A3) |
| 53 | `numerics/field.py` `Field.cone_violations` | concept row 'Positive cone' (one, landed); W7 'Wall law on its cone' |
| 57 | `numerics/measure.py` `DiscreteMeasure.singular_set_under` | W8a, concept row 'Discrete measure' (the invariance half of the weld) |
| 58 | `numerics/measure.py` `DiscreteMeasure.symmetry_groups` | W8a, as row 57 |
| 59 | `numerics/measure.py` `BundleMeasure` | W8a, concept row 'Discrete measure' ('widened to differing supports') |
| 62 | `numerics/metric.py` `HilbertMetric.kernel_basis` | as row 48 |
| 63 | `numerics/operator.py` `D2` | the TypeVar of LinearOperator.__matmul__'s quoted annotation (numerics/operator.py:929) |
| 71 | `numerics/posing.py` `ALPHA_MAP` | W5, concept row 'Question bound to a system' (spectral map derived from the parameter) |
| 75 | `numerics/space.py` `FunctionSpace.kernel_basis` | as row 48 |
| 78 | `numerics/spaces/full_field_space.py` `FullFieldSpace.kernel_basis` | as row 48 |
| 83 | `sn/angular/closure.py` `MorelMontryAngularSweep.angular_endpoint_defect_per_level` | W7 (#653's radial cases un-defer part of the curvilinear build) / curvilinear SN cluster |
| 100 | `sn/sweep/pairing.py` `pair_diffusion_limit_consistent` | W4, concept rows 'Spatial scheme' and 'Angular closure' |
| 104 | `transport/boundary_conditions/realization.py` `RealizerTable.catalogue` | door: the method tables' admission matrix, asserted through it by 13 test files (5 methods); nothing in W0 opens the realizer-table cycle (3.4) |
| 120 | `transport/operators/multiplication_operator.py` `MultiplicationOperator.from_mesh` | a documented survivor pinned by a design gate (transport/fields/test_meshless_construction.py:245) |
| 127 | `transport/residuals/angular_residual.py` `AngularResidual.relative_to` | W5, concept row 'Iteration record' (one estimator, one norm) |
| 128 | `transport/residuals/scalar_residual.py` `ScalarResidual` | W3, concept row 'Flux and source members' (the role forest) |
| 134 | `numerics/quadrature/registry.py` `_gl1d_invert` | W8a: a value of the required field `QuadratureSpec.degree_of_exactness_for` of the kept catalogue (G1-a); ruled out of W0 by 3.0 item 2 |
| 135 | `numerics/quadrature/registry.py` `_lebedev_invert` | as row 134 |
| 136 | `numerics/quadrature/registry.py` `_ls_sn_invert` | as row 134 |
| 137 | `numerics/quadrature/registry.py` `_product_invert` | as row 134 |


**Assignment of the 45 keep rows to rows modules** (G8; by the file the sub-unit edits; each `KEPT` entry is `(file, qualname, kept_for)` with the reason of the table above): `rows_w0_dead_a.py` holds rows 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 29, 35 (the `geometry/` rows), 57, 58, 59 (`numerics/measure.py`) and 134, 135, 136, 137 (`numerics/quadrature/registry.py`, a file dead-A edits); `rows_w0_dead_b.py` holds rows 38, 47, 48, 53, 62, 71, 75; `rows_w0_roles.py` holds row 63 (`numerics/operator.py`, the file `roles` edits first); `rows_w0_dead_c.py` holds rows 83 and 100; `rows_w0_dead_d.py` holds rows 78 (`numerics/spaces/full_field_space.py`, a `numerics/spaces` file dead-D edits), 104, 120, 127 and 128. That is 30 + 7 + 1 + 2 + 5 = 45 rows. The four rows deferred by 3.1 (52, 99, 106, 116) are ADDITIONAL `KEPT` entries of dead-B, dead-C, dead-D and dead-D (reason `deferred, private state`). `rows_wave2.py` holds the wave-2 survivors. After all four sub-units and `roles`, `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_retirements.py` asserts 45 + 4 = 49 `KEPT` rows defined (its instrument for done-when 1 of 3.7).

### 3.4 The three factory cycles (Table 3 of the dead census): which back reference is removed in W0, and which waits

None is removed in W0. The census (`census_dead_roles_cycles.md` Table 3) classes them; this specification rules each out by the structural reason in the table below, and section 8 decision 6 records the ruling. The orchestrator's decision D3 (amendment of 2026-11-20) is that all three WAIT: the realizer table for W7, `Chord.partition` for W6, `DiscreteMeasurePartition.measure` for W8a (210 test reads consume it). No text of this specification opens a factory cycle in W0.

| cycle | edge needed? (census) | disposition | wave | structural reason |
|---|---|---|---|---|
| realizer table and its realizers (`transport/boundary_conditions/realization.py`: `RealizerTable`, `ResponseRealizer`, `SourceRealizer`, `DeckRealizer`) | NOT needed: no realizer reads its table after class creation; registration happens in `__init_subclass__(table=...)` at 22 production class statements (5 SN, 6 diffusion, 4 CP, 4 MC, 3 MoC) and 32 test class statements | WAITS | W7 (the boundary needs typed faces and the realizer table is the boundary-law ontology, architecture E) | opening it re-orders three method modules (`cp/solver.py`, `mc/solver.py`, `moc/geometry.py`: the table is created first and filled by later class statements) and deletes the open/sealed lifecycle (`SealedTableError`, `TableNotSealedError`); the CP, MC and MoC modules are not sharpened on their own terms in W0 (direction of development), and the lifecycle is the object W7 redesigns; doing it in W0 would edit the boundary package twice (here and again in W7). `[M]` blast: 22 production, 32 test class statements, 20 `RealizerTable(` constructions in 5 test files |
| `ConcentricPartition`, `Chord`, `Rays` (`geometry/chord.py`) | the FIELDS are needed (`Chord.partition` read on 9 lines, `Rays.chord` on 11); the FACTORY METHODS close the cycle | WAITS (declared inside one concept) | W6 (`ConcentricPartition` is the one locator and chord, #596) | opening it means moving the constructors to the products (`Chord.through(partition, line)`), which trades the spelling the plan words on the partition (B5: "the chord on the existing concentric partition"; unit `insulation` spells its chord call as `partition.chord(line)`): the same choice as the algebra's closure; callers of the four factories: production 2 lines (both in `chord.py`), `derivations/` 28 lines, tests 103 lines |
| `DiscreteMeasure` and its partition (`numerics/measure.py`: `DiscreteMeasurePartition.measure`) | NOT needed by production: 0 reads of `.measure` in production, 1466 reads of `label` and `indices` (a runtime spy over 98 tests `[M]`); tests read `.measure` 210 times | WAITS | W8a (concept row "Discrete measure") | removing the `measure` field needs a restriction verb on `DiscreteMeasure` taking the index partition, and the 210 test reads of `.measure` (`[M]` the census, carried: the ten sites named earlier were a sample, not the population) re-point to it (behavioural: rewire); that verb is a design of the measure's invariance half, which W8a unifies, and dead-A already holds `numerics/measure.py` for row 60 (`equispaced`; rows 57 to 59 are KEPT) |

Consequences stated so that 3.3 and 3.4 agree: row 104 `RealizerTable.catalogue` stays (3.0 item 4); `derivations: 0` files edited holds for all five units; `cp/solver.py` is in `insulation`'s set for the three imports only.

### 3.5 Retirement-audit items for the unit (the shared discharge; the sub-unit sections add the row-specific ones)

A.1 and A.2 per row: the appendix columns "production references" and "test files naming it" are the graph and text results; the caveat is `[M]` 8 of the 81 suspect (file, method) pairs were resolved by `nexus callers`, the other 73 cannot be (instance calls through variables and 25 property or field rows are not call edges); the nexus graph is stale (`nexus staleness`: 80 of 147 doc pages stale at the measurement), so an added caller is not seen: the count of callers is a lower bound. A.3: not applicable (no guarded type's constructor is changed). B.4: `nexus dead-references` and `.claude` text grep by tense; B.6: the string residues above; B.8: the paragraph that justified a retired thing (the tense table's `PRESENT` rows are these); B.9: the `:by:` targets (`eq/ref` rows) and the labelled equations `sn-angular-endpoint-defect-eq` (row 83, kept), `bundle-measure-disintegration` (row 59, kept) are NOT touched; C.12/C.13: the per-row fates and the markers listed (rows 6, 36, 44, 45, 72, 73, 81, 85, 86, 91, 93 carry markers); D.14 to D.18 as D6 and D.17 for no re-baselined pin (the unit re-baselines none); E.19: the tense tables; E.20: a rewire onto `@` or `+` of an operator (row 67 `InverseMetricOperator` to `RieszRaiseOperator`, dead-B): the surrogate must honour the contract, and here the contract differs in its ends (3.6 dead-B); E.21: rows retired onto a shared helper (none); E.22: the module `orpheus/sn/loss_representation/assembly.py` is a mass delete: grep the NAME `assembly` and the import strings (`git grep -n "loss_representation.assembly\|loss_representation import assembly" -- orpheus tests tools docs`); G.24 and G.25: one commit per sub-unit, leaves first; the move writes nothing into a moved tree.

### 3.6 Sub-unit notes (specific actions beyond the appendix tables)

**dead-A** (rows 6, 13, 24, 27, 28, 30, 31, 32, 33, 34, 36, 39, 40, 41, 51, 54, 55, 56, 60, 61, 72, 73, 132, 133; rows 134 to 137 are kept, 3.0 item 2).
- The selector closure, by decision G9: the RULED rows retire: `select_quadrature` (73), `SelectionLog` (132, with `SelectionLog.summary`, row 72, and the members that live inside the class, among them `SelectionLog.angular_symmetry`), `QuadratureSelectionError` (133) and `QuadratureSpec.structural_flags` (the ruled closure member of row 73, 3.0 item 3); the nested `_selection_log`, the local `angular_symmetry` and `valid_flags` are inside `select_quadrature`'s body and go with it. NOT touched: the four inverters `_gl1d_invert`, `_lebedev_invert`, `_ls_sn_invert`, `_product_invert` (rows 134 to 137, with their nested helpers `n_min`, `_even_at_least`, `ceiling`, `n_mu_min`, `n_phi_min`), the field `QuadratureSpec.degree_of_exactness_for` that holds them, and `_lebedev_orders` (all kept catalogue data, 3.0 item 2), and `DiscreteMeasure.degree_of_exactness` (KEPT: the measure's defining invariant, named in 11 test files and `discrete_measures.rst`). KEPT and untouched: `AngularSymmetry.admits_symmetry` (G1-b, in the closure of row 73 and named by the ERR-081 legs), `QuadratureSpec.axis_aligned`, `half_range_clean`, `quadrature_registry` and `_registry`. The class docstring of `QuadratureSpec` (the sentence "`structural_flags` reads them off the rule" and the ERR-118 `[M]` paragraph) is a present-tense description of the retired method: rewrite it to the catalogue fields only, keeping the 2026-10-19 measurement as past tense. The `__all__` strings and imports: `orpheus/numerics/quadrature/__init__.py` (the import of `select_quadrature`, `SelectionLog`, `QuadratureSelectionError`, and their `__all__` entries) and `registry.py`'s own `__all__` (`"QuadratureSelectionError"`, `"SelectionLog"`, `"select_quadrature"`). The module docstring of `registry.py` (it describes the selector at length: the `:func:`select_quadrature`` roles at its top, "The :func:`select_quadrature` function returns both ...", "opt-in convenience") is a present-tense description of the retired function: rewrite the module docstring to describe the catalogue only.
- Test files with selector references: `tests/gates/numerics/test_registry.py` (36 occurrences of the three names), `test_angular_symmetry.py` (3), `test_invariance.py` (2), `test_slab_orbit_space.py` (2). Rule: a test whose subject is the selector (stage 0 to 3 behaviour, `SelectionLog`, the rejection reasons) is deleted with its subject (its law, the admission of a rule by a geometry's symmetry, lives in the `admits_domain`/`admits_symmetry`/`domain_refusal` legs that already exist beside it); a test that obtained a rule through the selector is rewired to the factory (`Quadrature.<factory>(...)`) with a literal parameter.
- The geometry rows' markers: rows 6 (`tangencies`: three `verifies('geometry-directions-at')` tests keep their marker, the helper import changes) and 36 (`permutes`: `test_o4_a_duplicated_point_is_refused_and_guard_two_is_the_one_that_fires` keeps `catches('ERR-073')` and still reaches the guard).
- `derivations/` consumers of dead-A rows: row 30 `domain_extent_cm` (3 test files; the 4 tests in `tests/gates/derivations/test_sood2003_to_geometry.py` are behavioural: REWIRE to `breakpoints[-1] - breakpoints[0]`, never delete); `laguerre` in `derivations/` is `np.polynomial.laguerre` (a homonym); no `orpheus/derivations/` file imports a dead-A symbol `[M]` (the review: row 55's derivations hits are homonyms).
- Gate leg for the rows 134 to 137 (kept): the KEPT entries of D-A.7 assert the four helpers still defined, and after the unit `git grep -c "degree_of_exactness_for" -- orpheus/numerics/quadrature/registry.py` still prints the field's declaration and its four assignments (the field has no production reader once the selector is gone; it is data for W8a). The tests that call the field outside the selector's tests are unchanged (their lines are not edited).

**dead-B** (rows 0, 1, 2, 42, 43, 44, 45, 46, 49, 50, 67, 68, 69, 70, 74, 76, 77, 130, 131; row 52 is deferred by D2).
- Eigenvalue cluster: retire `rayleigh_quotient_iteration` (row 50, with its two inner helpers `border` and `step_norm`) and move `direct_eigenvalue`, `dominant_eigenpair`, `_sign_normalised` to `tests/gates/numerics/_dense_eigen_oracle.py`. Delete `TestRayleighQuotientIteration` (6 tests) and the `getattr` binding line in the same commit (the review's finding 1). The three `monkeypatch.setattr(_eig, "dominant_eigenpair", ...)` / string spies (`test_eigenvalue.py` two places, `test_homogeneous.py` one): read each: a spy that substitutes the function to count or fake a call that production no longer makes is deleted with the production path (a counting spy is a call site written as a string, `retirement-audit` members list); one that substitutes the oracle for a test-side call is re-pointed at the oracle module. The three string fixtures in `tests/gates/test_layer_imports.py` (the reference-insulation refused legs importing `dominant_eigenpair`): replace the imported name with a live name of the same module (`power_iteration` from `orpheus.numerics.eigenvalue`); the legs test the import SHAPE, not the name. This is an edit of `test_layer_imports.py` lines OTHER than roles' and insulation's (section 7).
- Row 67 `InverseMetricOperator`: retire, and its rewire is a RE-TYPING, not a delegate swap (decision D5). `[M]` read at `orpheus/numerics/operator.py` (structure check): `InverseMetricOperator` has `domain` and `codomain` both returning its space `V` (an endomorphism `V -> V`) with `apply` and `apply_transpose` both returning `space.apply_inverse_metric(x)`; `RieszRaiseOperator` (`LinearOperator[Domain, Domain]`, in the same file) has domain `V.dual()` and codomain `V` and the same array action. The two are one array action under two typings; the retired one is the mistyped arrow, the pseudo-inverse of the metric typed as an endomorphism of the primal space. The laws of `test_inverse_metric_operator.py` are therefore restated, not copied, for the surviving type: (i) `domain == space.dual()` and `codomain is space` (the old `test_the_operator_is_bound_to_the_space_on_both_ends`, which asserted `operator.codomain is space` on BOTH ends, is INVERTED to state the two ends, and is not kept as written); (ii) symmetry: `apply_transpose(f)` is bit-equal to `apply(f)`, written as `<apply(f), g> = <f, apply(g)>` over `V*` under the canonical pairing; (iii) `G G^+ x = x` off the null space of the metric `G` and exactly 0 on it, over the parametrised spaces including `product(4, 4)`. If a rewired assertion has no meaning on the surviving type it is DELETED with that reason in the test's docstring, not bent; `test_the_pseudo_inverse_masks_the_metrics_null_space` is compared with `test_riesz_legs.py::test_round_trip_is_the_range_projector_on_a_singular_metric` (`[R]`: the bodies were not compared by the specification) and deleted if the second states the same law. D.14 mutation, in-process: `monkeypatch` `RieszRaiseOperator.apply` AND `RieszRaiseOperator.apply_transpose` each to return `-x`; each rewired test must go red. E.20 applies. This file (`operator.py`) is edited by unit `roles` first (section 7).
- Rows 68, 69, 70 (`OperatorPencil.is_regular`, `rhs_rank`, `apply_at`): two SN docs front-matter strings in `docs/theory/methods/sn/index.rst` (the `pencil:` and `splitting:` entries, lines 41 and 46 at base) make PRESENT claims about `is_regular`, `rhs_rank` and `law_residual` (the last is dead-C row 98): rewrite them.
- Rows 76, 77 (`angular_flux_space`, `scalar_flux_space`), 74 (`project_live`): per appendix; row 52 (`Field.where_largest`) is deferred (D2: `_index_tuples`).
- The eigenvalue docs: `docs/api/numerics.rst`, `docs/theory/foundations/infinite_medium.rst` (line naming `rayleigh_quotient_iteration`), `operator_inverse_family.rst`, `diffusion_1d.rst`, `sn/history.rst`.

**dead-C** (rows 37, 80, 81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 141, 142, 143; row 99 is deferred by D2 and row 126 has its home in dead-D).
- `orpheus/sn/loss_representation/assembly.py` is deleted whole (rows 85, 86, 141, 142, 143: `[M]` `ast` over the module's top level: exactly the five functions of those rows and the module's `__all__` assignment, so the move is a whole-module move; one `fuller-view oracle` module); its tests keep four `verifies` markers (`loss-rep-walk-order-rows` among them) on their functions, which move to the new test-side module with the functions; its only production reference is inside the module itself (`[M]` the review's production-reference column: `assembly.py` line 126, a string in the module); `docs/theory/foundations/operator_inverse_family.rst`, `docs/theory/methods/sn/history.rst`, `docs/theory/methods/sn/loss_representation.rst` name it.
- Row 91 `predicted_kernel_dimension` (move: counting-law oracle): five `verifies` markers stay in `tests/gates/sn/operators/test_loss_kernel_gauge.py`, `test_loss_nullspace_reflective_box.py`, `test_predicted_loss_kernel.py`; the function moves to the one module `tests/gates/sn/operators/_kernel_dimension.py` (3 consumers, 3.1 rule 1) as a free function `predicted_kernel_dimension(problem)` importing the private `_anova_dimension` from production (unruled closure member, 3.1 rule 3); `orpheus/sn/operators/loss_kernel_gauge.py` is also edited by unit `posing` (u1.2) in a different function (disjoint; section 7).
- Rows 93 and 94 (`SNProblem.loss_kernel_gauge`, `SNProblem.boundary_face_layout`): rewire 6 and 3 test files; row 93's markers (`test_loss_nullspace_reflective_box.py::test_the_T_plus_R_split_is_exact_per_quadrature`, `test_every_entry_gauges_its_trace.py::test_every_MIRROR_EVEN_functional_is_blind_to_the_gauge`) stay.
- Row 84 `build_coupled_system` (rewire to `build_within_group_system`; production references `orpheus/sn/__init__.py` import and `__all__` string); `orpheus/sn/coupled_system.py` is edited by `roles` first.
- Row 96 `SNSolver.compute_group_absorption_rate`: also edit `tests/gates/cross_method/test_one_parse_on_production.py` `_site_calls` table (the two dotted strings), and rewire the sn test that calls it.
- Rows 80 to 82: `_MMHalfGrid.downstream_per_ordinate`, `.upstream`, `compute_psi_half_per_level` move to `tests/gates/sn/sweep/curvilinear/_psi_half.py` with 8 `verifies` markers staying on `test_diamond.py` and `test_dd_recurrence.py`; the homonym `.upstream` in `derivations/` and `tools/cone_store/placement.py` (parameters and loop variables) is NOT touched.
- The kept neighbours: row 83 `MorelMontryAngularSweep.angular_endpoint_defect_per_level` (W7), row 100 `pair_diffusion_limit_consistent` (W4).

**dead-D** (rows 79, 101, 102, 103, 105, 107, 108, 109, 110, 111, 112, 113, 114, 115, 117, 118, 119, 121, 122, 123, 124, 125, 126, 129, 138, 139, 140; rows 106 and 116 are deferred by D2; row 126 `AngularResidual.balance_map` lives in `transport/residuals/`, so by the package cut it belongs here, and it was listed in dead-C and absent from D-D.1 before: one home now, section 8.3).
- `HarmonicMomentFlux.truncate` (row 111) and `truncated_head` (row 140): ONE commit (the review's closure; `truncate` is `truncated_head`'s only production caller; `numerics/spaces/__init__.py` re-exports it by import and `__all__`; the tests: delete `numerics/test_heads_on_the_axis.py` H7 and the rows of `test_moment_head_axis_built_premise.py` that test it, and the own-SUT tests in 6 files).
- The `TransferOperator` family: rows 123, 124, 125 (`foldable_part`, `residual_part`, `is_foldable_into_sigma_r`), 139 `_sibling`, 138 `foldable_sigma` (re-point `tests/gates/sn/architecture/_config.py` and edit the allow-list of `tests/gates/sn/acceleration/test_dsa_rate.py`; the kept `MaterialXSField.foldable_sigma` stays and a `_HOMONYMS`-style live assertion in the rows module asserts it still exists at `orpheus/transport/mesh/material_xs_field.py`), 121/122 (`LegendreMomentTransfer.on_basis`, `on_frame`: test helpers; production mints through `TransferOperator._moment_transfer`), 119 `IsotropicTransfer.dense_per_material`. `transfer.py` and `isotropic_transfer.py` are edited by `roles` first.
- Row 105 (`FaceField.face_views`): `orpheus/transport/fields/_bases.py` is also edited by `defects` (#709, `MomentField`, disjoint classes): dead-D first. Row 106 (`BoundaryField.from_face_arrays`) is deferred (D2: `cls._face_space_of`).
- Row 115 `FissionKernel.dyad` (rewire `:1055` to `np.outer`; delete `:479-480` in `test_kernels.py`; 3 calls by receiver `[M]`), row 79 `ScalarTraceSpace.slot_shape` (rewire one behavioural test to `layout.faces[face].shape`), rows 101 to 103, 107 to 110, 112 to 114, 117 and 118 (free functions in `test_kernels.py`), 126, 129 per appendix (row 116 is deferred, D2; row 126 moved here from dead-C).
- Row 120 `MultiplicationOperator.from_mesh` is KEPT: `tests/gates/transport/fields/test_meshless_construction.py` fails if it is absent ("the documented survivors still exist"), so the census's move would redden it.

### 3.7 Done when (the `dead` unit; each predicate is one command)

1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_retirements.py` exits 0 and reports the definition, reference and keep legs for every integrated sub-unit (once all four and `roles` are in: the 92 dead rows plus the 5 role rows, 97 `Retired`/`Moved` rows, plus the one ruled closure member `QuadratureSpec.structural_flags`, 98 `Retired`/`Moved` entries, are defined 0 times; the 45 keep rows (the 41 `keep` verdicts and rows 134 to 137) plus the 4 deferred rows, 49 `KEPT` rows, are defined 49 times).
2. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_error_catalogue_reconciles.py tests/gates/test_vv_harness_audit.py tests/gates/test_docstring_xrefs.py tests/gates/test_layer_imports.py tests/gates/test_test_layers.py` exits 0.
3. `.venv/bin/nexus dead-references` (after `sphinx-build -E -W docs docs/_build/html` exits 0) lists none of the retired names (the 98 entries of item 1) as a dead reference introduced by the unit (the pre-edit set diff is empty or smaller).
4. The reference leg of `tests/gates/test_retirements.py` (item 1) reports 0 references to any retired definition by import resolution, with the positive and negative controls of 1.5 in its body.
5. `git grep -n -e rayleigh_quotient_iteration -e select_quadrature -- orpheus tests tools examples ':(exclude)tests/gates/_retired' ':(exclude)tests/gates/test_error_catalogue_reconciles.py'` prints 0 lines (the two exclusions are the rows modules, which name the symbols by design, and the catalogue gate's own synthetic control; control: the same command without the exclusions at `b68f8b8c` prints more than 0). `degree_of_exactness_for` is not in this predicate: the field and its four values stay (3.0 item 2).
6. For each sub-unit, the test files of its appendix pass one at a time (`.venv/bin/python -O -m pytest -p no:cacheprovider <file>`), serially.

## 4. Unit `roles`: the block and system role tags retire, with no replacement

**Outcome.** The tag that no production branch reads (`BlockRole`, `SystemRole`, the metaclass `_BlockRoleMeta`, the three marker classes `BulkOperator`, `FullOperator`, `BoundaryOperator`, the stamp `stamp_boundary_role`, the attributes `block_role` and `system_role`, the two joins `_join_block_roles` and `_join_system_roles`, and `FusedAtom._adopt_expansion_roles`) leaves production, with the tests that only read it, and NOTHING replaces it in W0. The structural reason is not the census of readers but the law the tag does not have: a tag on an operator of the algebra is useful only if it composes, and this one has no composition law. `FULL` "has off-diagonal coupling" (`orpheus/numerics/operator.py`, the `BlockRole` docstring), but the join `BULK ⊔ BOUNDARY = FULL` is returned for a block-diagonal sum `C + B = diag(A_bb, A_ss)`, which has no off-diagonal block, so the join reads a different object (the set of spaces touched) from the one documented (the support in the blocks); the adjoint keeps the tag only because `FULL` forgets the direction of the coupling (`bs` against `sb`); and under composition `{bs}∘{sb} = {bb}` while `{bs}∘{bs} = 0`, with every factor tagged `FULL`. A tag that is a homomorphism for neither `+` (inconsistent with its own docstring) nor `∘` (no law) cannot be derived, so it is retired, and the ruling of 2026-10-20 ("a role is a fact about an atom", carried by `test_c42` of `tests/gates/numerics/test_linear_combination.py`) is REPEALED by this unit for that reason; the commit message says so and gives it; the orchestrator, not the builder, writes the entry for the repeal in the rulings ledger of `.claude/plans/ontology_chart.md` at close-out, with this structural reason (the plan is not in a unit's edit set, so the builder does not edit it). What a later wave derives is the support of an operator over the direct sum `V_bulk ⊕ V_boundary`, an operator datum that W2 re-spells with the operator members; it is not this tag under another name. Posing item u1.6 (the comment at `numerics/operator.py` "the single fact `OperatorSum.apply` dispatches on") is inside this unit and is rewritten with the class (4.1).

**Why nothing replaces it** `[M]` (`census_dead_roles_cycles.md` Table 2; `m_roles1.py`): no production branch reads a role. The only read that decides anything is `_BlockRoleMeta.__instancecheck__` (`numerics/operator.py`), reachable only through `isinstance(x, BulkOperator | FullOperator | BoundaryOperator)`, and production contains 0 such calls (`ast`, also `issubclass` and `match`; positive control: 55 in `tests/`). Attribute reads `.block_role`/`.system_role` (`ast.Attribute`, `Load`): 0. The 10 `getattr(x, "block_role", None)` reads are 8 propagations (`AdjointOperator.__init__` two, `_zero_of` two, `OperatorSum.__init__` two, `FusedAtom._adopt_expansion_roles` two), 1 delegation (`RealizedBoundary.block_role`, `sn/boundary/realized_operators.py`) and the 1 `__instancecheck__`. Implicit readers probed on one fixture (`2.0 * IdentityOperator()` built twice, one copy stamped with both tags): `==`, `hash` and `_algebra_key` read the pair as equal (control: `2*I` against `3*I` is unequal on all three); `content_digest` refuses an operator (`ContentlessError`), so no memo key reads it; `pickle` carries the instance attribute (propagation); `values/content.py` reads `vars(value)` only for subclasses of an array leaf kind. One composite fixture; the dataclass operators stamp a class attribute that is not a field (12 class-level writes, read by `ast`, not probed).

### 4.1 File set (by tree) and the exact edits

Population `[M]` (`scratch/ontology/w0/spec/amend/roles_count.py`, AST over `Name`, `Attribute`, `alias`, `keyword`, `arg`, `def`/`class` names and string constants equal to a name, over `orpheus`, `tests`, `tools`, `examples`, `student_resources`; the 12 names `BlockRole SystemRole _BlockRoleMeta BulkOperator FullOperator BoundaryOperator block_role system_role _join_block_roles _join_system_roles stamp_boundary_role _adopt_expansion_roles`): 392 hits in 36 files, of which production 17 files and tests 19 files (117 hits in the 17 production files and 275 in the 19 test files, `[M]` the structure check's probe over the same 12 names; the earlier "105 in 14 files and 268 in 17" counted only the first 10 names); `derivations/`, `tools/`, `examples/` and `student_resources/` 0; 11 definitions, all of the retired family; `isinstance`/`issubclass` on a marker: production 0, tests 55; positive control: `orpheus/numerics/operator.py:471` `ClassDef BoundaryOperator`. Every done-when of this unit uses these 12 names and these denominators.

Hard, production (17 files; the lists are the appendix R-1; the sites, by file, from the census Table 2.1, symbol decides):
- `orpheus/numerics/operator.py`: the 8 definitions (`BlockRole`, `SystemRole`, `_BlockRoleMeta`, `BulkOperator`, `FullOperator`, `BoundaryOperator`, `_join_block_roles`, `_join_system_roles`); the `__all__` strings naming them; the two `None` defaults on `LinearOperator` (`block_role`, `system_role` class attributes); the propagations `AdjointOperator.__init__` (the `getattr(x, "block_role", None)` read and the assignment), `_zero_of`, `OperatorSum.__init__`, `FusedAtom._adopt_expansion_roles` (the whole method); the comment block `# Block-role classification (Issue #208)` (u1.6): delete ONLY its tag sentences (the lines from "`BlockRole` classifies a leaf by WHICH blocks" to the end of the block) and KEEP the true mathematics, the statement that on `V = V_bulk ⊕ V_boundary` a linear operator is by the biproduct theorem a 2x2 block matrix `A = [A_bb A_bs; A_sb A_ss]` with `A_bb : bulk -> bulk`, `A_bs : boundary -> bulk`, `A_sb : bulk -> boundary`, `A_ss : boundary -> boundary`, retitled `# The block structure of an operator on a direct sum` (the same statement is the displayed matrix of `docs/theory/foundations/coupled_block_operator.rst`, label `coupled-block-2x2`, and of `boundary_conditions.rst` near the line holding `\begin{bmatrix} A_{bb} & A_{bs} \\ A_{sb} & A_{ss} \end{bmatrix}`, `[M]` read); the `BoundaryOperator` docstring is deleted with the class, and its measured lesson (an affine boundary law breaks the Arnoldi relation GMRES's residual depends on: `|B(0)| = q`, `|B(2x) - 2B(x)| = q`, a raised `ConvergenceClaimError`) is ALREADY carried by `docs/theory/foundations/boundary_conditions.rst` (`[M]` read: the table row `‖B(0)‖∞` at about lines 4170 to 4188, the Arnoldi relation paragraph with `ConvergenceClaimError` at about lines 4250 to 4259, and the `P3` entry near line 7081, all under the anchor `bc-affine-source-channel` at about line 3564), so nothing is moved; the `D2` TypeVar is NOT touched (row 63, kept). 72 of the 145 lines naming the tags in production are in this file.
- `orpheus/numerics/coupled_system.py` (import and the `block_role`/`system_role` class attribute on the coupled operator), `orpheus/sn/coupled_system.py` (the `SystemRole` import and the instance write at the quoted line `op.system_role = ...`), `orpheus/diffusion/operators.py` (import; two class-attribute writes), `orpheus/sn/operators/boundary.py` (two imports; five class-attribute writes), `orpheus/sn/operators/radial_characteristic.py` (import; four writes), `orpheus/sn/operators/streaming.py` (import; one write; the call of `_adopt_expansion_roles`), `orpheus/sn/operators/scheduled_invertible.py` (the call of `_adopt_expansion_roles`), `orpheus/sn/boundary/realized_operators.py` (the `RealizedBoundary.block_role` delegating property; the two production calls of `stamp_boundary_role`, in `_attenuated_kernel_operator` and `ZeroMap._derive_inner`: each reduces to its argument), `orpheus/diffusion/boundary_realizer.py` (the call in `_albedo_operator`, which reduces to its argument), `orpheus/transport/boundary_conditions/realization.py` (the definition of `stamp_boundary_role`, its import of `BlockRole`, the instance write, the `__all__` string), `orpheus/transport/boundary_conditions/__init__.py` (the import of `stamp_boundary_role` and its `__all__` string), `orpheus/transport/operators/angular_lift.py`, `isotropic_transfer.py` (import; two writes), `lift.py`, `multiplication_operator.py` (import; one write each), `transfer.py` (the import of `BlockRole`, never used).
- 13 prose files (a docstring or comment names a tag; rewrite or delete, read for tense): `orpheus/diffusion/augmented_mesh.py`, `orpheus/numerics/registry.py` (names `_BlockRoleMeta` as an example in a comment), `orpheus/numerics/spaces/angular_trace_space.py`, `orpheus/sn/splitting.py`, `orpheus/transport/__init__.py`, `orpheus/transport/fields/angular_boundary_source_sink.py`, and 7 test files (below).
Hard, tests (20 files): `tests/gates/diffusion/test_augmented_mesh.py`, `test_boundary_realizer.py`, `test_diffusion_realizer_table.py`, `test_operators.py`; `tests/gates/numerics/_assembly_functor.py`, `test_coupled_operator.py`, `test_linear_combination.py`, `test_operator_protocols.py`; `tests/gates/sn/architecture/test_realized_boundary.py`; `tests/gates/sn/operators/test_capability_survival.py`, `test_linear_combination_sn.py`, `test_operator_block_role.py`, `test_psi_half_coupling.py`, `test_sn_boundary_operator.py`; `tests/gates/sn/solve/test_si_single_primitive_contract.py`; `tests/gates/sn/verification/analytical/test_mms_declared_inflow.py`; `tests/gates/test_layer_imports.py` (the three string fixtures below; integrates FIRST among the units that edit it); `tests/gates/test_wave2_retirement.py`; `tests/gates/transport/test_boundary_conditions_package.py`, `test_bulk_lift.py`. Prose-only tests (docstring or comment, 7): `tests/gates/numerics/test_angular_trace_space.py`, `test_face_streaming_normal.py`; `tests/gates/sn/operators/test_bc_extraction_matvec.py`, `test_declared_law_is_linear.py` (module docstring, two places; the six `catches("ERR-075")` functions in the tests of this family do NOT read the tag: `[M]` 0 role names in the body of any of them), `test_fission_operator.py`; `tests/gates/sn/regression/README.md`; `tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py`.
Hard, docs (8 pages): `docs/api/numerics.rst` (the heading naming the three markers as Protocols), `docs/development/evidence/lessons.md` (one line), `docs/theory/foundations/boundary_conditions.rst` (about 21 lines on the tag or the stamp, 2 on a retired homonym class `BoundaryOperator`, which stay), `coupled_block_operator.rst` (the `SystemRole` lattice section and 5 mentions), `operator_algebra.rst` (two lines), `operator_tensor_network.rst` (two lines, the markers), `docs/theory/methods/diffusion_1d.rst` (two lines), `docs/theory/verification/error_catalog.rst` (ERR-075's three lines are HISTORY and stay; two homonym lines stay). `.claude/agents`, `.claude/rules`, `CLAUDE.md`: 0 hits `[M]` (positive control: the same pattern finds `.claude/skills/vv-principles/error_index.md:126`, generated, a copy of ERR-075's title: unchanged, history). `.claude/agent-memory/` (91 lines in 34 files) and `.claude/plans/` (42 lines in 15 live plans) are NOT in this unit's edit set. The `BoundaryOperator` homonym warning: `BoundaryOperator` is also the name of a class retired at Wave 7 and a suffix of the live `SNBoundaryOperator`, `DiffusionBoundaryOperator`, `SNMaskedBoundaryOperator`, `RadialCharacteristicBoundaryOperator`; no rewrite touches those.
Tense table: Appendix R-3 lists, by file, the lines naming each name and their tense (`PRESENT` rewrite; `role`: drop; `past`: stay).
Further edits that the first draft did not list (amendment of 2026-11-20): (1) `tests/gates/numerics/_assembly_functor.py` line 37, `self._adopt_expansion_roles()` in a test double: delete the call (every consumer of the functor raises `AttributeError` at construction otherwise); (2) the production docstring lines that name the retired names and are not in R-3: `orpheus/numerics/operator.py` (the `:attr:` role to `block_role` and `system_role` at about line 1977, the `:meth:` role to `_adopt_expansion_roles` at about line 2428, a line of the `BoundaryOperator` docstring at about line 498 that goes with the class), `orpheus/transport/boundary_conditions/realization.py` about line 182, `orpheus/sn/boundary/realized_operators.py` about line 214, and the test docstrings `tests/gates/sn/architecture/test_realized_boundary.py` about line 358 and `tests/gates/sn/operators/test_declared_law_is_linear.py` about line 19 (`[M]` each line read and found to contain a name of the 12); (3) the docs lines of R-3b; (4) comments: every comment that names one of the 12 names is read for tense (no census of comments was taken: that is the builder's re-read, with `git grep -n -P` over `orpheus tests` for the 12 names, `[M]` examples met: the `GAP-1` comment of the deleted protocols file, `orpheus/sn/splitting.py` about line 347, `tests/gates/sn/operators/test_fission_operator.py` about line 326).
New: the re-home of `test_prescribed_inflow_realizes_the_same_object_vacuum_does` (see 4.3).

The marker-class rows 64, 65, 66, `_BlockRoleMeta` (145) and `stamp_boundary_role` (144) are `dead/reviewed.json` rows (Appendix R-2) and live in the shared gate's `rows_roles.py`.

### 4.2 The gates

**R1. The retirement rows of the 12 names** (the shared engine of 1.5, `tests/gates/test_retirements.py`, created BEFORE any dead sub-unit by this unit, with the rows module `tests/gates/_retired/rows_w0_roles.py` and the port `rows_wave2.py`). `Retired` rows (definition site and import resolution): in `orpheus/numerics/operator.py` `BlockRole`, `SystemRole`, `_BlockRoleMeta`, `BulkOperator`, `FullOperator`, `BoundaryOperator`, `_join_block_roles`, `_join_system_roles`, `FusedAtom._adopt_expansion_roles`, `LinearOperator.block_role`, `LinearOperator.system_role`; in `orpheus/transport/boundary_conditions/realization.py` `stamp_boundary_role`; in `orpheus/sn/boundary/realized_operators.py` `RealizedBoundary.block_role` (the delegating property: the one `def block_role` of the tree besides the class attribute). `RetiredName` rows (bare name; each with its measured homonym count 0 in the row's comment, `[M]` `scratch/ontology/w0/spec/amend/roles_count.py`): `BlockRole`, `SystemRole`, `block_role`, `system_role`, `_join_block_roles`, `_join_system_roles`, `_adopt_expansion_roles`, `stamp_boundary_role` (the name written in about 20 files, `block_role` 91 hits in 25 files). `BoundaryOperator` is NOT a `RetiredName`: it is a suffix of live class names and a retired Wave-7 class named in docs; it is covered by its `Retired` row, whose definition-site and import-resolution binding cannot be reddened by `SNBoundaryOperator`. The population is every `.py` file of `orpheus`, `derivations`, `tests`, `tools`, `examples` and `student_resources` outside `tests/gates/_retired/`; the gate does not scan `.claude/` or docs (the tense table R-3 covers them). **First red** `[M]`: the 12 names have 117 hits in 17 production files and 275 hits in 19 test files at `b68f8b8c` (the definitions at `numerics/operator.py` among them; probe: `roles_count.py` and the structure check's own AST probe, which agree on the total 392). Negative leg: `isinstance(x, BulkOperator)` in a synthetic source is found. **Mutation** (G6, the root-path mechanism): copy `orpheus/numerics/operator.py` to a temporary root at the same relative path, re-add `block_role = None` to the copy of `LinearOperator`, call the walk on that root, and require red for the `LinearOperator.block_role` row and the `block_role` name row.

**R2. The tests that read the tag** (55 test functions and doubles in 17 files; every assertion of the 55 was read, `roles_tests_read.txt`): delete 42 (41 test functions and the double `_Tagged.__init__`), rewire 8, move 2 (the two assertions of `test_operator_protocols.py` that survive its deletion, below), edit 2 doubles, keep and re-home 1. The full per-function fates are Appendix R-4. Summary: delete `TestBlockRoleEnum`, `TestSystemRoleEnum`, `TestIsinstanceMarkers`, `TestGenericOperatorsAreUnclassified` (apart from two assertions that move), the 11 tag-only tests of `sn/operators/test_operator_block_role.py`, `TestSystemRoleLattice` (four tests) and `TestSystemRoleStamp`, `test_c31` and `test_c42` of `numerics/test_linear_combination.py` WITH the join (their mutation `ROLE_HEAD` has no target left; the `rests_on` line at `test_linear_combination.py` that names `test_c31` goes with its decorated test: of 971 `rests_on` lines in `tests`, 1 names a test this unit deletes), the four tag-only tests that walk every wall (`diffusion/test_augmented_mesh.py`, `diffusion/test_boundary_realizer.py`, `sn/operators/test_operator_block_role.py` two): `[M]` reviewed by coverage (`dynamic_context=test_function`, `python -O -m pytest -p no:cacheprovider` per file, serially; 6 files: 17, 31, 22, 52, 50 and 70 passed): 0 of 4 is the sole exerciser of any production line (the unique lines of `test_augmented_mesh.py:257` (13 lines, `operator.py` role joins, `__add__`, combination flattening, `OperatorSum`) are covered by `test_boundary_realizer.py::TestComposition::test_a_mixed_wall_realizes_to_the_sum_of_the_scalars_action`; the 106 unique lines of `test_operator_block_role.py:179` are covered by `test_sn_boundary_realizer.py`, `test_sn_boundary_operator.py`, `test_sn_realizer_table.py`; the other two 1 and 0 lines; the claim is line coverage over 6 files, not path combinations); rewire the 10 behavioural tests by dropping the tag assertions only (one surviving test, `test_realized_boundary.py:390`, `rests_on` `test_the_algebra_object_delegates`, so that one is rewired, never deleted); edit the two doubles `Mat.__init__` (`numerics/test_linear_combination.py`, `sn/operators/test_linear_combination_sn.py`) to drop the `block_role` parameter; keep and re-home `test_prescribed_inflow_realizes_the_same_object_vacuum_does` (behavioural, names the tag in prose only) into a surviving file of the same directory (`tests/gates/sn/operators/test_sn_boundary_operator.py`), because `test_operator_block_role.py` is otherwise deleted. `[M]` No function among the 55 carries `catches` or `verifies` (AST over decorators and `pytestmark`); the module-level `pytestmark` of the emptied files are tier marks only (`foundation`, `l0`), so no claim marker migrates (C.13). The six `catches("ERR-075")` functions stay green by design and are not edited.

**`test_operator_protocols.py` and the re-home** (decision R2). Audit item D.16(b), measured: the two assertions of `tests/gates/numerics/test_operator_protocols.py` that survive once the tag is dropped are `isinstance(IdentityOperator() + DiagonalOperator(np.ones(3), shape=(3,), spans=(0,)), OperatorSum)` and `isinstance(tensor_a + tensor_b, OperatorSum) and len(.terms) == 2` for two `TensorProductOperator((d0, d1))`. The candidate duplicates do not state them: `tests/gates/numerics/test_operator.py:415-437` asserts `isinstance(.., OperatorSum)` for sums and scalings of MATRIX operands only, and `tests/gates/numerics/test_operator_capability_predicates.py:277` builds `IdentityOperator() + DiagonalOperator(_C, shape=(3,), spans=(0,))` as a fixture row of capability predicates and asserts no `isinstance` or term count (`[M]` read at both sites); no site asserts the two-term `OperatorSum` of two tensor products. They are therefore NOT duplicates: the file `test_operator_protocols.py` is deleted whole, and the two assertions move to `tests/gates/numerics/test_linear_combination.py` under names that state what they assert, `test_a_sum_of_the_identity_and_a_diagonal_keeps_both_terms_as_an_operator_sum` and `test_a_sum_of_two_tensor_products_is_an_operator_sum_of_two_terms`; its module docstring, `_Tagged` double and `GAP-1` comment go with the file (they speak only of roles). The bulk-lift test `tests/gates/transport/test_bulk_lift.py::TestLaws::test_capability_axes_and_role` is renamed `test_capability_axes` (the tag assertions dropped). The re-homed `test_prescribed_inflow_realizes_the_same_object_vacuum_does` needs its two helpers, which live in the deleted file: `_FACE = "xmin"` and `_boundary_under_test(declared)` (`test_operator_block_role.py` lines 127 to 150 at base), whose imports are `WallCondition`, `TranslationDeck`, `Mesher`, `StructuredGeometry`, `CellsByCount`, `SNProblem`, `Quadrature` (imported inside the function) and `placeholder_materials`; both move with the test into `tests/gates/sn/operators/test_sn_boundary_operator.py` (the builder checks that file's imports first). Done-when 3 therefore lists 19 test files, not 20: `test_operator_protocols.py` and `test_operator_block_role.py` no longer exist.

**R3. The three silent fixtures** (names inside strings that stay green when the names are gone; `[M]` `m_roles2.py`: `_check_source` only parses the string, and with `BlockRole` replaced by `NoSuchNameAtAll` the line-644 leg returns the same violation `geometry → numerics`):
- `tests/gates/test_layer_imports.py`: the refused leg `("geometry/deck.py", "from orpheus.numerics.operator import BlockRole\n", "geometry → numerics")` (in the list of refused legs, the quoted `BlockRole` leg): replace the imported name with a live one, `LinearOperator`, keeping the module-level shape; the two admitted legs `("sn/x.py", "...realization import stamp_boundary_role")` and `("diffusion/x.py", "...realization import stamp_boundary_role")` (the lines after the `Layer DAG WP-E1` comment in `_ADMITTED_LEGS`): replace the imported name with a live name of `orpheus.transport.boundary_conditions.realization` (`RealizerTable`). The legs test the import SHAPE and the layer edge, not the name; no gate except the string leg of 1.5 can tell: that leg parses each fixture string as Python and walks it with the one finder, so after the replacement it finds a LIVE name in each of the three fixtures and, before it, it reports the three at the base hash (its first red: the lines at about 644, 667 and 668 of `tests/gates/test_layer_imports.py`).
- `tests/gates/test_wave2_retirement.py`: this unit ports the file's name walks into the one engine (1.5): `RETIRING_NAMES` become `RetiredName` rows and `_HOMONYMS` become `KEPT` entries in `tests/gates/_retired/rows_wave2.py`, and the file's own walks are deleted. The `_HOMONYMS` entry `"stamp_boundary_role": "orpheus/transport/boundary_conditions/realization.py"` and its sentence in the module docstring are therefore not carried into `KEPT` (the stamp retires; `rows_w0_roles.py` carries it); without that omission the live-survivor assertion would redden (`name in spelled` fails: the stamp no longer exists), which was the by-design red of the previous plan.
- `tests/gates/transport/test_boundary_conditions_package.py`: `test_stamp_boundary_role_is_copied_in` is deleted (R2).

### 4.3 What stays bit-identical, laws, raw inputs

No new function or value is created, so no new law is introduced. The unit must change no number: `==`, `hash` and `_algebra_key` of an operator do not read the tag (probed), so operator identity and the memo keys are unchanged; the whole of `tests/gates/numerics/test_linear_combination.py`, `tests/gates/sn/operators` and `tests/gates/diffusion` run green before and after (a gate green before and after is not evidence of anything by X1: the evidence here is the probe above, which names the readers, and R1's red).

### 4.4 Retirement-audit items for the unit

A.1/A.2 as above; B.4: `.venv/bin/nexus dead-references` after a docs rebuild, and a `.claude` text census sorted by tense (the generated index copy of the ERR-075 title and agent memory are history or out of scope); B.6: the three fixtures (R3) and the `__all__` strings; B.7: the docs pages speak of "block role" and "system role" as concepts: read `boundary_conditions.rst` and `coupled_block_operator.rst` for paraphrases (`role` in a column heading, "the tag", "classification", "stamp"), `git grep -n -i -E "block role|system role|BOUNDARY stamp|role lattice" -- docs`; B.8: the `SystemRole` lattice section of `coupled_block_operator.rst` is a paragraph that JUSTIFIED the tag: re-verify its sentences against the code that remains (the lattice's join is retired); C.12/C.13: R2; D.14: no test compares the tag with itself through a wrapper (all 42 deletes read the tag, a tag-only assertion); D.15/D.16: the six ERR-075 catchers do not read the tag; E.19: the tense table (Appendix R-3: PRESENT rewrite, role drop, past stay); E.20: none; E.22 (a mass delete is a retirement: grep the module or script name): the re-homed test `test_prescribed_inflow_realizes_the_same_object_vacuum_does` is named by PATH in six places that no instrument reads: `docs/theory/foundations/boundary_conditions.rst` at about lines 3693 (a vv-status rationale), 4544 and 4624 (a swap-mutation command that discriminates the test), `tests/gates/sn/operators/test_declared_law_is_linear.py` about line 108, `tests/gates/sn/operators/test_capability_survival.py` about line 123 and the deleted protocols file about line 14 (`[M]` `git grep -F test_operator_block_role.py`); each is re-pointed to `test_sn_boundary_operator.py`, and the path inside the command at about line 4624 is re-pointed to `test_sn_boundary_operator.py` by the builder, who also records the command in the sentence around it as DATED history: measured at the old path `test_operator_block_role.py` before the re-home, not re-run (a builder's worktree lacks the `wf25_bcB_mutplug` plugin; the orchestrator does not re-run it either); G.24, the retirement order, as an 8-level dependency order inside ONE commit: (1) the external writers, readers and imports (the 12 class-level writes, the instance write in `orpheus/sn/coupled_system.py`, the `RealizedBoundary.block_role` property at `realized_operators.py`, the imports, the test-double call of F5); (2) the stamp and its re-export (its three callers first, then `stamp_boundary_role` and its `__all__` string); (3) `FusedAtom._adopt_expansion_roles` and the propagations (`AdjointOperator.__init__`, `_zero_of`, `OperatorSum.__init__`) with their two callers; (4) the joins `_join_block_roles` and `_join_system_roles`; (5) the `None` defaults on `LinearOperator`; (6) the three marker classes; (7) `_BlockRoleMeta`; (8) the enums `BlockRole` and `SystemRole` and the `__all__` strings; G.25: none.

### 4.5 Done when

1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_retirements.py tests/gates/test_wave2_retirement.py tests/gates/test_layer_imports.py tests/gates/test_test_layers.py tests/gates/test_error_catalogue_reconciles.py` exits 0 (`test_wave2_retirement.py` keeps the tests that are not name walks, 1.5).
2. The retirement engine's walk over `orpheus`, `derivations`, `tests`, `tools`, `examples` and `student_resources` (outside `tests/gates/_retired/`) reports 0 hits for each of the 12 names; positive control, the same walk at `b68f8b8c` (`git show b68f8b8c:<file>` into a temporary root, or a `git worktree` at that hash) reports 117 hits in 17 production files and 275 hits in 19 test files over the same 12 names.
3. The 19 test files of 4.1 that still exist (`test_operator_protocols.py` and `test_operator_block_role.py` are deleted; `test_linear_combination.py` and `test_sn_boundary_operator.py` gain the moved tests) pass one at a time (serially), the 7 prose-only test files unchanged in behaviour.
4. `sphinx-build -E -W docs docs/_build/html` exits 0; `.venv/bin/nexus dead-references` lists no row for the 12 names beyond the past-tense history of ERR-075.
5. `git grep -n -w -e BlockRole -e SystemRole -e _BlockRoleMeta -e BulkOperator -e FullOperator -e stamp_boundary_role -e block_role -e system_role -e _join_block_roles -e _join_system_roles -e _adopt_expansion_roles -- orpheus tests tools examples student_resources ':(exclude)tests/gates/_retired'` prints 0 lines except the lines of the `BoundaryOperator` homonym class names (`SNBoundaryOperator` and the like are different words and do not match `-w`; the live `BoundaryOperator` word may appear only in a retired-class history line, each read); control: the same command at `b68f8b8c` prints more than 0 (117 + 275 hits over the same names, 4.1). `BoundaryOperator` is checked by the engine's `Retired` row, not by this grep, because the word is a documented homonym.

## 5. Unit `defects`: three refusals, each tagged with the wave that retires it, and the sub-unit `angular-mass` (5.7)

**Outcome.** #709, #710 and #711 are closed by hand with the integrating hash as REFUSALS (binding decision 5): a `raise` that survives `python -O`, tagged `ELEGANCE-DEBT[guard]` with the issue and the wave whose design makes the mistake unspellable; each lands with a gate holding the constructible input it rejects. The root fixes are later waves' (#709 in W3, where a member follows its space; #710 in W1, where a factor's metric has one home derived from its generator; #711 in W8a, where a measure's identity is its content); the wave's done-when is the refusal's retirement. #708 is NOT a refusal and is not in 5.1 to 5.6: its measured verdict (the section's divisor, `math.fsum`, is the right one on the 29 of 120 rules where the two spellings disagree) made it a root fix, sub-unit `angular-mass`, 5.7. #392 is closed as already fixed (5.5). #713 is corrected as an issue and not fixed (thermal hydraulics is out of scope, B13).

### 5.1 File set (by tree)

Hard: `orpheus/transport/fields/_bases.py` (#709: `MomentField.__post_init__`), `orpheus/numerics/axis.py` (#710: the end of `Axis.__post_init__`, and the new property `FromQuadrature.measure`), `orpheus/numerics/generating_measure.py` (#711: the intern table and `GeneratingMeasure.__post_init__`), `tests/gates/test_elegance_debt_is_tagged.py` (the wave leg, 5.3; read-only in the file-set census, edited here by decision 9), and the one existing test the #709 refusal reddens, `tests/gates/transport/fields/test_role_partners.py` (`TestIntoRole::test_moment_family_fields_ride_across`, rewired IN PLACE, 5.2).
New: `tests/gates/transport/fields/test_moment_field_agrees_with_its_space.py` (#709: the space is minted by `tests.gates.transport._phase_space_cases.moment_space`, which a transport-layer gate may import, so the gate lives at the layer of the invariant and builds no `SNProblem`), `tests/gates/numerics/test_axis_weights_agree_with_origin.py` (#710), `tests/gates/numerics/test_generating_measure_name_is_its_content.py` (#711).
Prose-only docs candidates for the three tags (pages naming `exact_recurrence`, `MomentField` or `spatial_moments`; 8 of 153 pages; read for tense, edit only a present-tense claim made false): `docs/theory/conventions/indexing_and_layout.rst`, `docs/theory/foundations/discrete_measures.rst`, `field_algebra.rst`, `frame.rst`, `operator_algebra.rst`, `spaces.rst`, `docs/theory/methods/sn/cartesian_multid.rst`, `docs/theory/methods/sn/history.rst`. Three docstrings that the refusals PROMOTE (`retirement-audit` D.15 and B.8: a duplicated field is now held equal by a refusal, and the prose says it is only checked "separately" or "kept for ergonomics"): the `GeneratingMeasure.name` docstring ("That the two paths agree is a separate, testable claim"), the comment on `MomentField.L` ("kept as a top-level field for ergonomic hot-path read access"), and the comment on the `Axis.weights` field; each is reworded to say that construction now refuses a disagreement and names the tag. The `ELEGANCE-DEBT[guard]` idiom exists in code: 9 production files carry the token today (`git grep -l "ELEGANCE-DEBT\[guard\]" -- orpheus`), among them `transport/fields/_bases.py` (#662 in `spatial_moments_per_axis_of`) and `numerics/measure.py`.
Order: after dead-A (`generating_measure.py`: rows 54, 55 move `HERMITE` and `laguerre` out; the #711 predicate is specified over the factories that remain; and after the precondition unit `quadrature-ties`, which also edits that file, section 7.2) and dead-D (`_bases.py`: row 105; row 106 is deferred).

### 5.2 Each refusal: its predicate, laws, message fragment, and the input that is red today

All three inputs are accepted today, with and without `-O` `[M]` (`m_def1.py`, run twice, `sys.flags.optimize` 0 and 1, identical output). Each tag is a comment immediately above the `raise`, in the form `# ELEGANCE-DEBT[guard] #<issue> (<wave>): <the invariant, one sentence>. Retires when <wave> <structural change that makes the mismatch unspellable>.` so that the token line carries the issue and the wave and a `retires` sentence follows on the same or the next lines (the existing ledger gate requires `#\d+` on the token line and a `retires|retired|until|when` sentence within three lines).

**#709, `MomentField` (W3).** `class MomentField(BulkField)` stores `L: int` and `spatial_moments: int = 1` beside the `space` that already fixes both (`transport/fields/_bases.py`, the class that stores the `L` truncation order). No `__post_init__` on the class today (the base `Field.__post_init__` checks `values.shape == space.shape`). Add one, calling `super().__post_init__()` first. The predicates read the quantity the layout NAMES, never a shape count that merely coincides with it today (`shape[0] == L + 1` holds only because both degree layouts lead with `L + 1`; a packed layout would break it):
- Predicate A: `self.head_layout.L == self.L` (`MomentField.head_layout`, the class's one door to the head's `DegreeLayout = HarmonicLayout | LegendreLayout`, `orpheus/numerics/coefficient_layout.py`, both arms carry `L`); else `ValueError` containing `disagrees with its space`, naming the class, the field `L` and `head_layout.L`.
- Predicate B: `self.spatial_moments_per_axis == self.spatial_moments` (the inherited instance property `spatial_moments_per_axis`, `_bases.py`, which reads the space; the static `spatial_moments_per_axis_of(space)` it calls is not re-spelled); else `ValueError`, same fragment.
- A space that is not a `TensorProductSpace` makes `self.head` / `self.head_layout` raise `TypeError` ("a moment field's space is a tensor product <angular head> ⊗ cells"); that is the right refusal for such a space and is kept (one existing test builds one, below).
Raw inputs: `L` an `int` (>= 0), `spatial_moments` an `int` (>= 1), `space` a `FunctionSpace`. Law: for every constructed moment field, `L` and `spatial_moments` are the ones its space states (so `dataclasses.replace(f, L=...)` re-runs the check: `[M]` `replace` re-runs `__post_init__`; unpickling does not run it). The law the constructor shares with its factories (`from_problem_and_L`, `zeros_for_problem_and_L`, every `replace(..., space=...)` site): each derives `L` from the space it mints, so each passes; the done-when runs the SN solves to say so (5.6 item 4). `[M]` (`scratch/ontology/w0/spec/amend/p709.py`, `-O`): on a slab `MaterialMesh` the layout reads `head_layout.L == 1` and `spatial_moments_per_axis == 1` for a field built with `L=1`, and the same fields built with `L=5` and with `spatial_moments=3` construct today.
First red `[M]`: with `mesh = MaterialMesh(Mesher(StructuredGeometry.slab((0.0, 1.0), (0,), left=0, right=0)).partition(CellsByCount.uniform_width(4)).mesh, placeholder_materials())` (`placeholder_materials` is `tests.gates.transport._materials.placeholder_materials`) and `sp = moment_space(mesh, posing_quadrature(mesh), 1)` (both from `tests.gates.transport._phase_space_cases`): `HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=5)` constructs (`head_layout.L` is 1, the field `L` is 5), and `HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=1, spatial_moments=3)` constructs (the space says 1). Admitted twin the gate keeps: `HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=1)` and `HarmonicMomentSourceSink` likewise.
By-design red `[M]` (census spy over 28 of 941 test files: 30 287 moment-field constructions, 0 `L` mismatches, 0 `spatial_moments` mismatches, 0 production callers reddened): exactly 1 test, `tests/gates/transport/fields/test_role_partners.py::TestIntoRole::test_moment_family_fields_ride_across`, builds `FunctionSpace(name="moments", shape=(2, 3, 2, 4))` and a `HarmonicMomentFlux` on it (the lines `m = HarmonicMomentFlux(values=np.zeros((2, 3, 2, 4)), space=space, L=1)`); `.head` raises `TypeError` there. REWIRE IN PLACE (behavioural: it asserts that `L` and `spatial_moments` ride across `into_role`): build the space with `moment_space(mesh, posing_quadrature(mesh), L)` of `tests.gates.transport._phase_space_cases` (`[M]` the transport layer already builds a legitimate moment space this way: `tests/gates/transport/fields/test_flux_vector_algebra.py` line 125). Direct constructors of the two leaves: production 8 sites in 4 files, tests 12 sites in 4 files (`git grep -n -P '\b(HarmonicMomentFlux|HarmonicMomentSourceSink)\('`); all 4 test files were in the spy. Tag text: `# ELEGANCE-DEBT[guard] #709 (W3): ... Retires when W3 derives L and spatial_moments from the space (the machinery that exists today: MomentField.head_layout and spatial_moments_per_axis; the call-site set W3 re-spells is the 8 production and 12 test constructors of the two leaves).`
Cost `[R]`: `head_layout` runs `moment_head_axis` per construction; not measured.

**#710, `Axis` weights against the origin measure (W1).** `Axis.weights` is an init field stored beside the origin that generated it, while the axis's `nodes` are a derived copy read off the origin: the asymmetry (a derived node set, a passed weight set) IS the defect and names what W1 changes (one factor's weights derived from its generator). `__post_init__` checks shape and finiteness only. At the END of `__post_init__` (after the canonicalisation that stores `None` for an all-ones weight vector and a read-only array otherwise), for the two arms that carry a measure only. The `FromQuadrature` arm gains a property `measure` returning `self.quadrature.measure` (its docstring already says that `FromQuadrature(q)` and `FromMeasure(q.measure)` mint one axis), so that one `match` decides the arm once:

```python
match self.origin:
    case FromMeasure(measure=measure) | FromQuadrature(measure=measure):
        carried = np.ones(shape) if self.weights is None else self.weights
        wanted = np.reshape(np.asarray(measure.weights, dtype=float), shape)
        # ELEGANCE-DEBT[guard] #710 (W1): the weights must be the origin measure's. Retires when W1 derives the factor's weights from its generator (the machinery that exists today: Axis.nodes, read off the origin; the call-site set W1 re-spells is the 5 production and 6 test direct constructors of Axis(..., origin=From*)).
        if not np.array_equal(carried, wanted):
            raise ValueError(f"axis '{self.label}': the weights disagree with the measure its origin records ...")
```
The tag sits immediately above the `raise`'s `if` (one placement for all three refusals: the token line is within three lines of the `raise`). Law: for `FromMeasure(m)` and `FromQuadrature(Q)`, the canonical weights (`None` read as ones) are bit-equal to `m.weights` / `Q.measure.weights`. Equivalent statement: an axis minted by `Axis.of_measure` or `Axis.of_quadrature` is the only construction that passes, except a measure whose weights are all ones (then `None` or an all-ones vector pass). Two invariants the law leaves unstated: the axis's `nodes` come from the origin and its `weights` from the caller, so the predicate holds at construction only; and `FromQuadrature` holds a LIVE reference (`Quadrature` is a mutable dataclass), so a later mutation of `Q` is not seen. Raw inputs: `weights`, `None` or an array of the axis shape; `origin`, a `FromMeasure`/`FromQuadrature` whose measure has `weights` of `n_points` entries. The `FromBasis` arm is NOT covered (`[M]` 1 of 669 constructions differs on purpose: the zero-mass refusal test `tests/gates/numerics/test_axis_marginal_lift.py`, so extending the predicate there reddens 1 test); `FromFrame` (525 constructions) has no measure to compare and `frame_faces.py` re-dresses with `weights=None` on `FromFrame`.
First red `[M]`: `m = gauss_legendre_on_polar_orbit(4, "x")` (it returns a `DiscreteMeasure`, `rules_1d.py`; no `getattr` fallback is needed) then `Axis("angular", (m.n_points,), weights=np.full(m.n_points, 7.0), origin=FromMeasure(m))` constructs (weights[0] = 7.0 against the origin's 0.3478548451...); `Q = Quadrature.gauss_legendre(n_ordinates=4)`, `Axis("angular", (Q.n_ordinates,), weights=np.full(Q.n_ordinates, 7.0), origin=FromQuadrature(Q))` constructs and is unequal to `Axis.of_quadrature(Q)`. By the predicate chosen here (decision 10), `Axis("angular", (m.n_points,), origin=FromMeasure(m))` (weights `None`) is also REFUSED when `m`'s weights are not all ones: it silently dresses a non-counting measure as the counting one, which is the defect. Admitted legs include a RANK-2 axis, `Axis("angular", (2, 2), weights=np.reshape(m.weights, (2, 2)), origin=FromMeasure(m))` (`[M]` it constructs today and `np.array_equal(np.reshape(w, (2, 2)), w)` is False, so it is the leg that reddens a comparison that drops the reshape: with rank 1 every other admitted input is identical with or without it). Spy `[M]`: 1212 of 1212 constructions on the two arms equal the origin measure bit for bit (636 + 576), so 0 tests and 0 production callers redden in 28 files. Direct `Axis(... origin=From*)` constructors: production 5 (3 in `axis.py`, 2 in `frame_faces.py`), tests 6 in 3 files, all in the spy.

**#711, `GeneratingMeasure` identity (W8a).** The class (`@dataclass(frozen=True)`, `name`, `support`, `exact_recurrence` with `compare=False`) compares by `name` alone; the docstring states that as the design; nothing checks the implementation against the name. Add a frozen record `_MeasureContent(support, alpha, beta)` (`support` a `Manifold`, `alpha` and `beta` the `bytes` of the first 8 float64 recurrence coefficients) and a module-level intern table `_NAME_CONTENT: dict[str, _MeasureContent]` above the class (the record names the three fields; the earlier nested-pair annotation disagreed with the flat content), and a `__post_init__`:
- content of a measure: `_MeasureContent(self.support, alpha_8.tobytes(), beta_8.tobytes())` with `alpha_8, beta_8 = self.recurrence(8)` (`K = 8`: the existing `_SYMMETRY_PROBE_ORDER` is 8, "a short prefix decides it"; the table records float64 bytes of the first 8 recurrence coefficients at `_VIEW_DIGITS` = 40, so two bit-equal recurrences are equal);
- if `name` is absent: record it; if present with a different content: `ValueError` containing `already names a different measure`, naming the name.
Law: within a process, one name is one content (`support` and the first 8 recurrence coefficients, bit for bit): name implies content. The converse FAILS and is stated, not hidden: content does not imply name, so two names for one content stay unequal. `[M]` (structure check) `LEGENDRE.on(-1.0, 1.0)` is named `legendre_on[-1.0,1.0]`, has the support and the bit-equal `recurrence(8)` of `LEGENDRE`, and `== LEGENDRE` is False; equality is strictly finer than content, and W8a is where identity becomes the content. Limit stated, not hidden: two measures that agree on the first 8 coefficients and differ later are not told apart. Raw inputs: `name` (any `str`), `support` (a `Manifold`), `exact_recurrence` (a callable `(n, digits) -> (alpha, beta)`, each a list of `n` `mpmath.mpf`). The three factories name by content: `_jacobi_name` returns the classical name for the three classical parameter pairs and `f"jacobi(a={a}, b={b})"` otherwise; `_affine_remap` appends the interval; `laguerre` its parameter, so one name cannot carry two contents through a factory `[R]` (read, not exhaustively tested). The admitted legs include the three classical pairs `[M]` (structure check, `p711.py`, `-O`): `jacobi(0, 0)` against `LEGENDRE`, `jacobi(-0.5, -0.5)` against `CHEBYSHEV_T` and `jacobi(0.5, 0.5)` against `CHEBYSHEV_U` agree on the name, the support and the bytes of `recurrence(8)` (3 of 3). The tag names both interning stores of the module: `Retires when W8a makes a measure's identity its content, which retires this process-global table and the lru_cache stores of _jacobi, _laguerre and _affine_remap (the machinery that exists today: GeneratingMeasure.recurrence; the call-site set W8a re-spells is the 7 direct constructors in generating_measure.py plus the moved HERMITE and laguerre copies)`.
First red `[M]`: `dataclasses.replace(CHEBYSHEV_T, name=LEGENDRE.name) == LEGENDRE` is `True`, the hashes are equal, and the `gauss(3)` weights differ. Spy: 16 constructions in 28 files, 0 violations, 0 tests, 0 production callers (a thin sample). Direct constructors `[M]` (`git grep -E 'GeneratingMeasure\('` outside the module): 0 in tests, tools, examples or production (the one earlier "1 in tests" is a docstring line of `test_legendre_basis.py`); 7 in `generating_measure.py`; the only test-side constructors will be the `HERMITE` and `laguerre` copies that dead-A moves. Cost, a design fact: one mpmath recurrence evaluation (8 coefficients) per construction at the module-level and factory constructions (`[M]` structure check: `recurrence(8)` is total on 34 sampled `jacobi`/`laguerre` parameters at 0.15 ms each, so the eager evaluation refuses no factory input), and a process-global table; both are stated in the tag sentence, neither is a reason. The `HERMITE` and `laguerre` factories leave production in dead-A; the moved test-side copies construct through the same constructor and pass through the refusal.

### 5.3 The gates

For each refusal one gate file (5.1) holding: (1) the refused input of 5.2 run IN A SUBPROCESS under `python -O` (`subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, cwd=<repository root>)`: the child imports `tests.gates...` helpers, which it finds only through its working directory), asserting a non-zero exit and the message fragment on stderr, so that "survives `-O`" is a gate and not prose; (2) the admitted legs of 5.2 in process (for #710 including the rank-2 axis; for #711 including the three classical Jacobi pairs). The earlier AST leg (an `ast.Raise` under an `ast.If`, no `ast.Assert` as the only guard) is DROPPED: the predicate "the only guard of that predicate" has no AST definition (`Axis.__post_init__` already holds 5 `raise` under `if`, so the leg was green on #710's input state), and leg (1) measures survival under `-O` directly (two instruments for one claim). **First red**: each refused input constructs today (`[M]` 5.2; a test-only input in the sense that no shipped construction is rejected: the spies found 0 of 30 287, 1212 and 16). **Mutations**, each with its mechanism: for leg (1), copy `orpheus` to a temporary directory, replace the `raise ValueError(...)` by `assert False, ...` in the copy of the file, set `PYTHONPATH` of the leg-(1) subprocess to that directory and print the child's `module.__file__` first (the `code-search` rule: the child must import the copy); under `-O` the refused input constructs again and leg (1) reddens; for leg (2), an in-process `monkeypatch` of the `__post_init__` of the class with a variant that compares `carried` with `wanted` without the reshape (#710): the rank-2 admitted leg reddens (the rank-1 legs do not, which is why the rank-2 leg exists).

**The wave leg** (`tests/gates/test_elegance_debt_is_tagged.py`; decision 9). The existing gate does not read the wave (`[M]`: it checks `#\d+` on the token line and a `retires|retired|until|when` sentence within three lines; 9 files carry the token today). Add `_WAVE_OF: dict[int, str] = {709: "W3", 710: "W1", 711: "W8a"}` and a test that: for each `(issue, wave)`, at least one token line in `orpheus/` matches `#<issue>(?!\d)` (so `#7090` is not `#709`) and carries the whole word `<wave>` on the SAME line (`re.search(rf"\b{wave}\b", line)`), and no token line matching `#<issue>(?!\d)` omits it. The tag sits immediately above the refusal (one placement for the three). Its positive control (inside the body): a planted token line `# ELEGANCE-DEBT[guard] #709: x. Retires when y.` is a violation, the line with `(W3)` is clean. **First red**: 0 of 3 tag lines exist today. **Mutation** (in-process, root-path walk): delete the `(W3)` from the planted line in a temporary copy of `_bases.py` and call the leg's function on that root: red. **What the leg does and does not do (stated, not hidden):** it makes a tag and its map entry agree, and nothing more; it does not force the retirement. If W1 lands and keeps both the map entry and the refusal, the leg stays green. The ratchet is the wave's own done-when predicate, `710 not in _WAVE_OF`, written into W1's plan; the wave that retires a refusal deletes its map entry and the refusal in the same diff (stated in the leg's docstring). Prose is not enforcement (X3): without this leg the wave is a comment.

### 5.4 Retirement-audit items

A.3: direct constructors of the three guarded types (above: `Axis(...)` 5 production and 6 test sites; the leaves 8 and 12; `GeneratingMeasure(` 7 and 1). The first one to redden in a wider run than the 28-file spy decides: the builder runs the full `tests/gates/numerics` and `tests/gates/transport/fields` directories serially once (not the suite). B.11: the message fragments `disagrees with its space`, `disagree with the measure its origin records`, `already names a different measure`: `git grep -n -F` over `tests/` for each before landing (0 expected). D.18: each refusal is the enforcement that keeps two copies equal (the field and its space): it lands with its witness in the same commit (5.3). C.12: the one rewired test. D.15 and B.8: the three docstrings of 5.1 that the refusals promote; E.19: the 8 docs pages by word match.

### 5.5 Issues

#392 is already fixed `[M]` (`m_def2.py`): `hasattr(SNProblem, "face_areas")` and `hasattr(SNProblem, "delta_A")` are False; control `hasattr(SNProblem, "boundary_face_layout")` True; AST over 260 production files: 0 `DeprecationWarning` name or attribute nodes (the 9 `warn` calls are other warnings; positive control for the scan: 3 test files name `DeprecationWarning`); the retirement comment is in `sn/problem.py` ("P4.1c", 2026-08-27). Close #392 by hand with the verifying hash; no gate is owed (a gate over an absent attribute is green before and after, X1). #713: corrected as an issue with the measurement: `thermal_hydraulics/solver_dae.py` is not a duplicate (14 definitions in each file, 8 identical, 6 differ with 88 changed lines in the diff; `examples/thermal_hydraulics/demo_thermal_hydraulics_dae.py` imports it; `docs/theory/thermal_hydraulics.rst` documents it); the twin is `THParams` and `THResult` only; not touched. #708: 5.7 (a root fix, closed with #288). The three issues' closing comments cite the refusal's file, the tag and the wave.

### 5.6 Done when

1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/transport/fields/test_moment_field_agrees_with_its_space.py tests/gates/numerics/test_axis_weights_agree_with_origin.py tests/gates/numerics/test_generating_measure_name_is_its_content.py tests/gates/test_elegance_debt_is_tagged.py tests/gates/transport/fields/test_role_partners.py` exits 0.
2. A probe file run with `.venv/bin/python -O`: the three inputs of 5.2 each raise `ValueError` (the `TypeError` case of #709 is a separate leg), with a positive control for each (the admitted twin constructs).
3. `git grep -n "ELEGANCE-DEBT\[guard\] #709\|ELEGANCE-DEBT\[guard\] #710\|ELEGANCE-DEBT\[guard\] #711" -- orpheus` prints exactly 3 lines, each carrying its wave (positive control, the same alternation with `#662` in place of the three issues prints the existing `_bases.py` tag line: 1 line, `[M]` structure check).
4. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics tests/gates/transport/fields tests/gates/sn` (directories, serial) exits 0 (the construction population the 28-file spy did not reach: #709's moment fields are built by every SN solve, and the spy reached 6 SN solve files).

### 5.7 Sub-unit `angular-mass`: #708, one definition of the angular mass, and the mixed source sum retired

**Outcome.** #708 is closed by hand with the integrating hash as a ROOT FIX (the orchestrator's ruling of 2026-11-20; section 8.4 items 8.28 and 8.29 and section 8.0 items 35 and 36), not as a refusal, and it settles #288 with it. Four changes, one law. (1) The mass of the angular axis, the divisor of the producer-side `1/W`, has ONE definition: the correctly rounded sum (`math.fsum`) of the stored weights, owned by the angular section `AxisSectionOperator.mass` (`orpheus/numerics/axis_marginal.py`), the object that diffusion-synthetic acceleration (DSA) already reads (`orpheus/sn/acceleration/dsa.py` lines 590 to 591; the brief's path `orpheus/transport/dsa.py` does not exist). (2) `AngularLift.total_weight` reads that mass instead of computing `float(np.asarray(self.frame.measure.weights).sum())`, so every reader of it (`_combine` before its rewrite, the transpose, the transfer's two redistributions) is on the one divisor. (3) `AngularLift._combine` injects the isotropic source through the section and adds within one class, so no `ScalarSourceSink + AngularSourceSink` sum remains anywhere in production. (4) The arm that made that sum, `AngularSourceSink.__add__` and `__radd__` with `AngularSourceSink.CONTAINS`, is deleted, and then `BulkField.CONTAINS` and `BulkField.__add__` (whose only job was the deferral to that arm) are deleted after a census shows 0 users (5.7.2). The mixed sum is then refused by `Field._check_partner`, the existing same-class gate, which is a `raise` and survives `python -O` `[M]` (5.7.2): no new guard, no `ELEGANCE-DEBT` tag, no wave. The fifth production spelling of the same quantity, `RadialCharacteristicReconstruction._ell0_scale`, is re-pointed too (orchestrator's answer, same date).

**Why it is one unit.** The defect is a second spelling of one quantity (the mass) welded to a dunder that adds a scalar rate into an angular rate with no division. Fixing the divisor without the arm leaves the arm, and the arm without the divisor leaves two divisors inside one emitted source (the isotropic part divided by one number, the anisotropic part by another). The root wave in the chart's B3 ("W3, members follow their space") is superseded: the injection is the typed section today, and the sum it replaces is retired now.

#### 5.7.1 File set (by tree)

Production, hard: `orpheus/transport/operators/angular_lift.py` (`total_weight`, `_combine`, one import, two docstrings), `orpheus/transport/fields/angular_source_sink.py` (the arm, `__radd__`, `CONTAINS`, three imports that become unused, the module docstring, `from_isotropic`'s docstring), `orpheus/transport/fields/scalar_source_sink.py` (module docstring and one comment block, no code), `orpheus/transport/fields/_bases.py` (`BulkField.CONTAINS` and `BulkField.__add__`), `orpheus/sn/operators/radial_characteristic.py` (`_ell0_scale`, one line, and `import math` if absent), `orpheus/transport/fields/__init__.py` (docstring only, lines 47 to 53, 5.7.8). Read and NOT edited: `orpheus/transport/operators/transfer.py` (its two readers, lines 690 and 704, read `self.total_weight` and are fixed by the edit of `total_weight`; it carries no sentence naming `np.sum` `[M]` `git grep -n -E 'np\.sum|\.sum\(' -- orpheus/transport/operators/transfer.py` prints 0 lines), `orpheus/sn/acceleration/dsa.py` (already reads the section's mass).
Tests, hard: `tests/gates/sn/primitives/test_typed_source_sinks.py` (8 tests, 5.7.6), `tests/gates/sn/operators/test_n2n_operator.py` (one test deleted, 5.7.6). New: `tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py` (G1), `tests/gates/sn/operators/test_radial_ell0_scale_is_the_correctly_rounded_mass.py` (G2), `tests/gates/transport/fields/test_scalar_source_does_not_add_to_an_angular_one.py` (G3), `tests/gates/_retired/rows_w0_angular_mass.py` (the rows of the shared retirement engine, 1.5; created by this sub-unit, picked up by the engine's `rows_*.py` glob).
Layers: G1 imports `orpheus.data`, `orpheus.geometry`, `orpheus.mesh`, `orpheus.numerics`, `orpheus.transport` only (the same set `tests/gates/transport/test_angular_lift.py` imports); G2 imports `orpheus.sn` and lives in `tests/gates/sn/operators/`; G3 imports `orpheus.transport` and `tests.gates.transport._phase_space_cases`. Run `tests/gates/test_layer_imports.py` and `tests/gates/test_test_layers.py` after creating them.
Docs, with the measured line of each present-tense claim the change makes false, at base `b68f8b8c` (the symbol decides; 5.7.8 gives the edit of each): `docs/theory/foundations/operator_algebra.rst` lines 5736 to 5742, 5793 to 5851 (the whole bullet), 5880 to 5890, 5116 to 5122 and 5142; `docs/theory/foundations/spaces.rst` lines 5961 to 5971, 6465 to 6472, 7043 to 7080; `docs/theory/foundations/field_algebra.rst` lines 1136 to 1146; `docs/theory/conventions/indexing_and_layout.rst` lines 901, 1948, 1961, 1966 to 1995; `docs/architecture/layering.rst` lines 163 to 166; `docs/architecture/conceptual_view.rst` lines 685 to 687; `docs/theory/methods/sn/adjoint.rst` lines 1509 to 1512; `docs/theory/methods/sn/cartesian_multid.rst` line 3144; `docs/theory/methods/sn/slab_multigroup.rst` lines 95 to 125, 432, 729 to 760 (read for tense; edit only a present-tense falsehood).

#### 5.7.2 The premise, measured (every fact a builder needs, with its command)

1. **The mixed sum constructs today and divides by nothing** `[M]` (probe `scratch/ontology/w0/spec/amend/p708/p_refused_input.py`, `.venv/bin/python -O` and without `-O`, identical): with `sn = MaterialMesh(Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(CellsByCount.uniform_width(20)).mesh, {0: get_mixture("A", "2g")})`, `interior = angular_full_field_space(sn, posing_quadrature(sn, 4), posing_scheme(sn)).interior_space`, `scalar_space = AxisRetractionOperator.of(interior, "angular").codomain`, `iso = ScalarSourceSink(values=np.ones(scalar_space.shape), space=scalar_space)`, `aniso = AngularSourceSink.zeros(interior)`: both `iso + aniso` and `aniso + iso` return an `AngularSourceSink` of shape `(4, 2, 20)`. (`angular_full_field_space`, `posing_quadrature`, `posing_scheme` are in `tests/gates/transport/_phase_space_cases.py`.)
2. **The two spellings of the mass differ on 29 of 120 rules** `[M]` (`scratch/ontology/w0/spec/amend/p708/p_rules.py`, `-O`, over the 120-rule population `_RULES` of `tests/gates/numerics/test_axis_marginal_lift.py`, line 398): `float(np.sum(w)) != math.fsum(w)` on exactly GL7, 11, 12, 15, 19, 23, 27, 29, 30, 31, 45, 50, 62; level-symmetric 10, 12, 16; `folded(2,12)`, `folded(4,12)`, `folded(4,16)`, `folded(6,12)`; `product(4,8)`, `product(8,12)`; Lebedev orders 3, 15, 17, 25, 27, 31, 41. On GL7 and GL12, `np.sum` gives `1.9999999999999998` and `fsum` gives `2.0`. The section's mass is already pinned to `math.fsum` on all 120 by `test_g5_the_nodal_mass_is_correctly_rounded` (that gate's arm `mass_numpy_sum` is the 29-red arm). The exact sum of the stored weights was computed with `fractions.Fraction` by the investigation memo (`scratch/ontology/w0/i708/memo.md`, probes `p1_divisors.py`, `p2_definition.py`): `fsum` is the correctly rounded value on 120 of 120 and `np.sum` is never closer to it than `fsum`; the difference is at most 4.7e-16 relative (Lebedev). I re-ran the 29-rule list and the two GL values; I did not re-run the `Fraction` leg.
3. **`np.sum` is storage-order dependent, `fsum` is not** `[M]` (`p_perm.py`): over the 120 rules, reversing the stored weights changes `np.sum` on 13 and sorting them by value changes it on 39; `math.fsum` is unchanged on 120 of 120 under both. (The memo's "83 of 120 over 20 random permutations" and my re-run's 85 depend on the random draw; publish the fixed-permutation counts.) On GL7, `float(w.sum())` is `1.9999999999999998` and `float(np.sort(w).sum())` is `2.0`: a mass that depends on the enumeration is not the measure's mass.
4. **Today `total_weight` is the `np.sum` spelling on a rule where it differs** `[M]` (`p_lift12.py`, `-O`; fixture of 5.7.5 G1): on GL7 and GL12, for each of `FissionOperator`, `ScatteringOperator` and `N2NOperator` built on `slab` GL n, `lift.total_weight` is `1.9999999999999998` while `AxisSectionOperator.of(lift._codomain_interior, "angular").mass` is `2.0`; on GL8 and GL4 both are `2.0`. On the same weights in sorted order (`p_perm_lift.py`: a `Quadrature` whose measure is `dataclasses.replace(q.measure, nodes=nodes[order], weights=weights[order])`, `order = np.argsort(q.weights, kind="stable")`) `total_weight` on GL7 is `2.0`: so the lift's divisor today depends on the enumeration of the rule.
5. **The census of divisor spellings, by `ast`** (`p708_census.py`, the probe of 5.7.9 item 2, `.venv/bin/python -I`): over `orpheus/` it prints 76 rows, 54 under `orpheus/derivations/` and 22 elsewhere. Of the 22: the production divisors of the angular mass are (a) `AngularLift.total_weight`, definition `angular_lift.py` line 323, read at `angular_lift.py` lines 440 (`_combine`) and 466 (`_interior_transpose`) and `transfer.py` lines 690 (`_redistribute_ordinates`) and 704 (`_redistribute_moments`): 1 definition, 4 reads; (b) `radial_characteristic.py` line 1142, `1.0 / float(quadrature.weights.sum())`, read at lines 1211 and 1267; (c) `dsa.py` lines 590 to 591, `self._sum_w = self._section.mass`, read at 692, 714, 717, already the section's mass. So 5 spellings read an `np.sum` divisor before this sub-unit (4 reads of (a), 1 of (b)), and 1 already read the section's; after it 6 of 6 do. The other 16 rows are not divisors of the angular mass: `dsa.py` 304 and 307 and `rules_sphere.py` 459 and 460 are `np.isclose` tolerance guards; `kernel_projection_gauge.py` 313, 317, 367 are Gram norms and a cosine on a coefficient space; `realized_operators.py` 522, `np.sum(in_weights)`, is the mass of the inflow HALF-trace measure (a different measure, not touched, section 8.1); `network.py` 218 and the three `mc/solver.py` rows are other quantities. The 54 `derivations/` rows are references independent by design and closed (binding decision 1 of 1.2): the MMS references spell `float(weights.sum())` and differ from the section by at most 4.7e-16, far below any manufactured-solution tolerance `[R]` (memo; not re-measured). Positive controls of the probe: it must print the three known sites `angular_lift.py` (a `SUM`), `radial_characteristic.py` (a `SUM`) and `dsa.py` (a `READ`), asserted inside the probe.
6. **No mass verb exists on the quadrature or its measure** `[M]` (`git grep -n -E 'def (mass|total_mass|total_weight|weight_sum)\b' -- orpheus` finds `AxisSectionOperator.mass` (`axis_marginal.py` line 308) and `GeneratingMeasure.mass` (`generating_measure.py` line 346, the zeroth moment of a weight function, a different object) and nothing on `Quadrature` or `DiscreteMeasure`; `DiscreteMeasure.integrate` is the only integral, at `measure.py` line 537). So the arc-family scale of `radial_characteristic.py` has no section to read (it is built from a `Quadrature` and a `FullFieldSpace` whose interior is the radial-characteristic split space, not an axis-built angular space) and takes `math.fsum` of the stored weights, which is the section's own definition applied to the same array (the section computes `math.fsum(covector * unit)` with `unit = 1`, `axis_marginal.py` line 479, so the two are bit-equal; G2 pins the equality).
7. **What the mixed sum's callers are** `[M]` (the investigation spy over 28 of 941 test files: 83 000 mixed sums, 82 990 from `angular_lift.py:440`, the other 10 attributed by that census to test code, not individually located; and, after the rewrite of `_combine` of 5.7.3, two in-process variants over the same 21 test files, `python -O`: `scratch/ontology/w0/spec/amend/p708/fix708.py`, which kept the arm and refused with its own `TypeError` and recorded every call into it: 0 calls, 343 passed; and `fix708b.py`, the RULED form, which deletes the five names in process (`AngularSourceSink.__add__`, `__radd__`, `CONTAINS`, `BulkField.__add__`, `CONTAINS`) and applies E1 and E5: 343 passed, 0 failed, serial, 11 minutes). The static census of `CONTAINS` (probe `p708_contains.py`, `ast` over `orpheus`, `tests`, `tools`, `examples`, `student_resources`; 5.7.9 item 3): 3 nodes: `_bases.py` line 373 (the declaration), `_bases.py` line 383 (its only read, in `BulkField.__add__`), `angular_source_sink.py` line 114 (the only non-empty declaration); 0 in `tests/`, `tools/`, `examples/`, `orpheus/derivations/`; positive control: the three known nodes. The only `__add__` definitions in the field family are `numerics/field.py` line 271 (the base), `_bases.py` line 375 (`BulkField`, a pass-through once `CONTAINS` is empty) and `angular_source_sink.py` lines 116 and 149; `BulkField.__add__` falls back to `super().__add__(other)`, which is `Field.__add__`.
8. **After the deletion the mixed sum is refused by an existing `raise`, under `-O`** `[M]` (the same probe, with `AngularSourceSink.CONTAINS = ()` and `__add__`, `__radd__` deleted in process, `python -O`): `iso + aniso` raises `TypeError: ScalarSourceSink [1/cm³/s] arithmetic requires a same-class partner; got AngularSourceSink [1/cm³/s/sr]` and `aniso + iso` raises the mirror `TypeError`. The raise is `Field._check_partner`, `orpheus/numerics/field.py` line 256 (`if not self._is_same_class(other): raise TypeError(...)`, not an `assert`). The fragment `same-class partner` is already pinned by 8 existing tests (`tests/gates/numerics/test_field.py` lines 185, 194, 212, 284; `test_coupled_operator.py` line 392; `test_typed_source_sinks.py` lines 515, 526, 538, the `HarmonicMomentSourceSink` cross-class rows): nothing re-words it.
9. **The rewritten `_combine` is bit-identical to the old one where the divisors agree, and slower by a constant** `[M]` (`p_time.py`, GL8 slab, `ScatteringOperator` `S_L0`, 20 cells, 2 groups): `np.array_equal(old, new)` is True; 7.6 microseconds per call old, 12.2 new (minimum of 5 repeats of 2000 calls). The float operations are the same: `phi / mass` then a broadcast, then one addition of the same two arrays (the sum of two IEEE doubles is commutative).
10. **The traced memo does not serve a stale divisor** `[R]` (read: `orpheus/numerics/traced_memo.py` keys a stored reference on the digest of the normalised source of every first-party def the generating run executed, so an edit to `total_weight`, `_combine` or `_ell0_scale` makes the entries whose run executed them stale and regenerates them; an entry whose run never executed them is unchanged and correct). Not measured; the builder's serial run (5.7.7) is the measurement.

#### 5.7.3 The exact edits

**E1. `orpheus/transport/operators/angular_lift.py`.**
- Import (line 114, `from orpheus.numerics.axis_marginal import AxisRetractionOperator`): add `AxisSectionOperator`.
- `AngularLift.total_weight` (the property at line 319; body line 323) becomes, and is NOT cached (a `functools.cached_property` would fail G1's mass-doubling leg, which reads the section on each call):
```python
    @property
    def total_weight(self) -> float:
        r""":math:`W = \langle m, u\rangle`, the mass of the emitted interior's angular axis (the producer-side :math:`1/W`): the angular section's mass,
        the correctly rounded sum of the stored weights, the one definition every reader of :math:`1/W` divides by (#708)."""
        return AxisSectionOperator.of(self._codomain_interior, "angular").mass
```
(`AxisSectionOperator.of` is memoised per (space, label) by `_lift`, so the read is a lookup. The old docstring's `W = \int_{S^2} d\Omega` was false: the stored weights sum to 2 on a slab rule and to 4 pi only on a sphere rule.)
- `AngularLift._combine` (def at line 426; body at lines 437 to 441) becomes:
```python
        section = AxisSectionOperator.of(self._codomain_interior, "angular")
        injected = AngularSourceSink(values=section.apply(iso.values), space=self._codomain_interior)
        return injected + (aniso if aniso is not None else AngularSourceSink.zeros(self._codomain_interior))
```
and its docstring says: the producer-side combine injects the isotropic source through the angular section `E` (divide by the mass, then scatter across the ordinates) and adds the anisotropic part within one class; the `cast` and the sentence about "the containment dunder's cross-class arm" and "the #288 principled LSP exception" are deleted (`cast` stays imported if used elsewhere in the file; the builder checks). `[M]` this body and the new `total_weight` were run in process (`fix708.py`): 343 tests pass (5.7.7).
- The module docstring lines 15 to 19 (`R_0 E M_0 / W`) are unchanged and true.

**E2. `orpheus/transport/fields/angular_source_sink.py`.** Delete the comment block and declaration `CONTAINS` (lines 108 to 114), the method `__add__` (lines 116 to 147) and `__radd__` (lines 149 to 153), and the section comment "Algebra extensions (over Field)" if nothing remains under it. The imports `replace` (line 62), `BulkField` (line 69) and `ScalarSourceSink` (line 71) become unused: delete each (`ScalarSourceSink` stays named in docstring roles; a docstring role resolves by the Sphinx cross-reference, not by a Python import, but `tests/gates/test_docstring_xrefs.py` resolves roles by importing the target: use the fully qualified role `:class:`~orpheus.transport.fields.scalar_source_sink.ScalarSourceSink``). Rewrite the module docstring (lines 21 to 55, the paragraph "The cross-class arithmetic ... is PRESERVED", the code block "Q_total = aniso + iso", the paragraph "Two named factories complement the dunder", the sentence at lines 54 to 55): the class's algebra is Field's, same-class and closed; a scalar source enters an angular one through `from_isotropic` (the section) and then a within-class add; `ScalarSourceSink + AngularSourceSink` is refused by `Field._check_partner`; the "refined Issue #207 principle" and the cross-class dunder are history: say so in the past tense with the date (2026-11-20) and #708. In `from_isotropic`'s docstring delete the paragraph that begins "It differs from the containment dunder" (line 176 to its end) and the "See also" entry about the dunder (line 202 to its end); keep the sentence that the divisor is the correctly rounded mass.

**E3. `orpheus/transport/fields/scalar_source_sink.py`.** Docstring only. Rewrite the paragraph at lines 19 to 23 and the code block at lines 36 to 50 (they present the cross-class sum as the design) and the comment block at lines 122 to 126 and 134 to 136, as in E2; line 112 ("The cross-class `iso + aniso` dunder of the pre-D-F class is ...") is read, and edited only if it states the dunder in the present tense.

**E4. `orpheus/transport/fields/_bases.py`.** Delete the comment block and `BulkField.CONTAINS` (lines 366 to 373) and `BulkField.__add__` (lines 375 to 387), after the census of 5.7.2 item 7 reads 3 nodes on the tree it is run on. `ClassVar` stays imported if used elsewhere in the file (the builder checks).

**E5. `orpheus/sn/operators/radial_characteristic.py`** (`RadialCharacteristicReconstruction.__init__`, lines 1141 to 1143): `1.0 / float(quadrature.weights.sum()) if arc_family else None` becomes `1.0 / math.fsum(np.asarray(quadrature.weights, dtype=float).ravel()) if arc_family else None`, with `import math` added if the module lacks it. The comment above it (lines 1135 to 1140) says the weight is `1/Σw` of the constant function: add that the sum is correctly rounded, the angular section's definition applied to a quadrature (no section is bound here).

**E6. The rows module `tests/gates/_retired/rows_w0_angular_mass.py`** (data of the shared retirement engine; created by `roles`, 1.5; this sub-unit follows it): five `Retired` rows, each defined today `[M]` (`git grep -n` at the lines above): `Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.__add__')`, `... 'AngularSourceSink.__radd__')`, `... 'AngularSourceSink.CONTAINS')`, `Retired('orpheus/transport/fields/_bases.py', 'BulkField.CONTAINS')`, `... 'BulkField.__add__')`; and the module's `KEPT` is empty (`KEPT = ()`, the engine accepts it). Mutation: as 1.5 (copy `_bases.py` to a temporary root, re-add `CONTAINS = ()` to `BulkField`, call the walk on that root, require red for that row).

#### 5.7.4 The law, whole

Let the angular axis carry the discrete measure `m = sum_n w_n delta_{Omega_n}` of the STORED weights `w_n` (nodal basis, unit `u = 1`). The retraction `R` is `phi -> sum_n w_n psi_n`, the section `E` is `phi -> (u / M) phi` broadcast over the ordinates, and `R` after `E` is the identity IF AND ONLY IF `M = sum_n w_n` exactly: in exact arithmetic `R E phi = phi * (sum_n w_n) / M`. So `M` is the mass of the discrete measure of the stored weights, not the continuum mass (2 or 4 pi; the stored weights differ from it by up to 3.1e-16 relative, `[M]` memo). The best float for that real number is its correctly rounded value, `math.fsum`, which does not depend on the order of the nodes; `float(np.sum(w))` is a pairwise sum in storage order and is not that value on 29 of 120 rules (`[M]` 5.7.2 items 2 and 3). Floating point `R after E` is the identity to the reduction depth and no better (`spaces.rst` lines 6111 to 6117: at most `N/3` ULP on the 120 rules), which is why no gate here compares a round trip bit for bit. The law of the unit: **every divisor of the angular mass in production is `math.fsum` of the stored weights, and the lift reads it from the section**: `lift.total_weight == AxisSectionOperator.of(lift._codomain_interior, "angular").mass == math.fsum(stored weights)` bit for bit, invariant under a permutation of the stored weights. Raw inputs of the laws: the stored weights (a 1-D `float64` array of `N` entries, `N >= 2`); the rules named in 5.7.2 item 2; a permutation of the nodes and weights together. The second law, the retirement: `ScalarSourceSink` and `AngularSourceSink` have different units (`SCALAR_RATE_UNITS` against `ANGULAR_RATE_UNITS`, differing by `1/sr`), so their sum is a dimensional error until the scalar is carried across the angular axis by the section, which divides by the mass: `E` is the only map from the first to the second.

#### 5.7.5 The gates

Each states its input, its first red in today's tree, and its mutation with the mechanism (1.3). What was run, exactly (drafts under `scratch/ontology/w0/spec/amend/p708/`, which a builder does not have, so each gate is described whole here): G1 is `gate1.py`, all three tests, red on today's tree (5 failed) and green under the ruled form (`fix708b.py`, 5 passed), with mutation (a) run (5 failed) and mutation (b) run (red on that reader alone); G2 is `gate2.py`, red on today's tree on the five differing rules (5 failed, 2 passed: the control and the sphere leg pass before and after, by design) and green under `fix708b.py` (7 passed); G3's refusal leg was run on a COPY of the tree with the real text edits of E2 and E4 applied (both orders raise `TypeError ... same-class partner` under `python -O`) and, on that copy with the `raise` of `Field._check_partner` replaced by `assert False, (`, the leg's mutation (5.7.5 G3), and `gate3.py` was run against today's tree (the two refusal legs red: both sums construct; the admitted legs and the within-class leg green). The gate files themselves do not exist in the tree; they are specified here, and a builder runs them.

**G1. `tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py`** (a foundation gate; `pytestmark = pytest.mark.foundation`). Fixture: `_mesh()` as in 5.7.2 item 1; `_posed(n, permuted)` builds `quadrature = posing_quadrature(mesh, n)` (optionally `_permuted(quadrature)`, a `Quadrature(measure=dataclasses.replace(q.measure, nodes=np.asarray(q.nodes)[order], weights=np.asarray(q.weights)[order]))` with `order = np.argsort(q.weights, kind="stable")`), `space = angular_full_field_space(mesh, quadrature, posing_scheme(mesh))`, and the lifts `FissionOperator.from_solver_data(mat_xs=mesh.mat_xs, space=space)` (key `F`) and `ScatteringOperator.from_solver_data(mat_xs=mesh.mat_xs, scattering_order=1, space=space)` (key `S_L1`). It does NOT write into `tests/gates/transport/test_angular_lift.py`'s module cache `_FFS` (a second module's private state). Three tests:
- `test_every_lift_divides_by_the_correctly_rounded_mass[n]`, `n` in (7, 12): first asserts the fixture (`float(w.sum()) != math.fsum(w)`, so the leg cannot pass vacuously), then for each lift `lift.total_weight == AxisSectionOperator.of(lift._codomain_interior, "angular").mass == math.fsum(w)`. **First red** `[M]`: today `1.9999999999999998 == 2.0` fails on `F` for n = 7 and n = 12.
- `test_the_divisor_does_not_depend_on_the_enumeration_order`: `{key: lift.total_weight}` on GL7 equals the dict on the same rule with the permuted enumeration (`_posed(7, permuted=True)`). **First red** `[M]`: today `{'F': 1.9999999999999998, 'S_L1': 1.9999999999999998} != {'F': 2.0, 'S_L1': 2.0}`.
- `test_every_reader_reads_the_sections_mass_object[n]`, `n` in (7, 12): for each lift, take `section = AxisSectionOperator.of(lift._codomain_interior, "angular")`, record the output of each reader on a fixed operand, then inside `monkeypatch.context()` set `section._mass = 2.0 * section.mass`, record again, and assert `np.testing.assert_array_equal(after, before / 2.0)` for each reader (division by `2M` is division by `M` then by 2, exactly in binary floating point; a reader that sums the weights itself does not move). The readers, 6 legs on 2 lifts: on `F` (a `FissionOperator`): `_combine(iso, None).values` (with `iso = ScalarSourceSink(values=standard_normal(lift._scalar_interior_space.shape), space=lift._scalar_interior_space)`) and `apply_transpose(cotangent).interior.values` (`cotangent` a `FullField(interior=AngularSourceSink(...), boundary=AngularBoundaryFlux.zeros(space.trace_space))`); on `S_L1` the same two plus `_redistribute_ordinates(flux).values` (`flux` the interior of a `FullField` of `AngularFlux`) and, on `lift.on_moment_domain()`, `_redistribute_moments(moments).values` with `moments = HarmonicMomentFlux(values=standard_normal(moment_space.shape), space=moment_space, L=1)` where `moment_space = on_moments._domain_interior` and `L = 1` is the fixture's scattering order, followed by `assert moments.head_layout.L == 1` (so the leg holds under #709's `MomentField` refusal, which derives `L` from the space). Each leg first asserts activation (`np.abs(before).max() > 0`). **First red** `[M]` (`p_which.py`): today all 6 legs read `False` (none divides by the section's mass), on GL12. DSA is NOT a leg: it stores `self._sum_w = self._section.mass` at construction (`dsa.py` line 591) and already reads the section; a gate would be green before and after (X1), and `orpheus.sn` is not importable from a transport-layer test. **Mutations**: (a) in process, `monkeypatch` `AngularLift.total_weight` with the old body (`float(np.asarray(self.frame.measure.weights).sum())`) while `_combine` stays on the section: legs 1, 2 and 3 redden, and the reader test reddens on `_interior_transpose` (measured with `fix708.py` and `FIX708_TW=old`: 5 failed); (b) monkeypatch one reader, for example `TransferOperator._redistribute_ordinates`, with a copy of its body that divides by `float(np.asarray(self.frame.measure.weights).sum())`: the reader test reddens on that reader alone.

**G2. `tests/gates/sn/operators/test_radial_ell0_scale_is_the_correctly_rounded_mass.py`.** Input: a cylinder `SNProblem` as `tests/gates/sn/operators/test_radial_characteristic_cylinder_arm.py` builds one (`StructuredGeometry.cylinder((0.0, 2.0), (0,), outer=0)`, `CellsByCount.uniform_width(10)`, `Quadrature.folded_product(a, b)`, a one-group mixture, `conditions_for(mesh, {0: WallCondition.specular(1.0)})`), then `RadialCharacteristicReconstruction(sn.radial_characteristic_field_space, coord=sn.coord, quadrature=sn.quad)`; the test asserts `op._ell0_scale == 1.0 / math.fsum(np.asarray(sn.quad.weights, float))` (a private attribute, the one stored constant the operator folds with; `[M]` no test reads `_ell0_scale` today: `git grep -n '_ell0_scale' -- tests` prints 0 lines). Rules: `(2, 12)`, `(4, 12)`, `(4, 16)`, `(6, 12)` (four of the 29: `np.sum` and `fsum` differ), `(2, 6)` (outside the 120-rule population, differs, and used by `test_quotient_adjoint_sn.py` and `test_dual_metric_sn.py`) and the control `(4, 8)` (they agree). **First red** `[M]` (`p_radial3.py`): today `op._ell0_scale == 1/np.sum` is True and `== 1/fsum` is False on the five (`p_radial3.py` for `(2, 6)` and `(2, 12)`; the other three by the same probe in its first form), both True on `(4, 8)`; for `(4, 12)` the stored scale is `0.07957747154594767` and `1/fsum` is `0.07957747154594769`. **Mutation**: monkeypatch `RadialCharacteristicReconstruction.__init__` with a wrapper that assigns `self._ell0_scale = 1.0 / float(quadrature.weights.sum())` after the original body: the five legs redden and the control does not (which is why the control exists). A leg on a sphere mesh asserts `op._ell0_scale is None` (the non-arc arm is unchanged).

**G3. `tests/gates/transport/fields/test_scalar_source_does_not_add_to_an_angular_one.py`.** The input of 5.7.2 item 1. (1) The refused input, run IN A SUBPROCESS under `python -O` (`subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, cwd=<repository root>)`; the child imports `tests.gates.transport._phase_space_cases`, which it finds through its working directory): `code` builds `iso` and `aniso` and evaluates `iso + aniso`; assert a non-zero exit and `same-class partner` on stderr; a second child does the same for `aniso + iso`. (2) The admitted legs, in process, built on the transport layer's own space (`from_isotropic` needs an SN carrier and lives in `tests/gates/sn`, 5.7.6): `injected = AngularSourceSink(values=AxisSectionOperator.of(interior, "angular").apply(iso_values), space=interior)`, then `injected + AngularSourceSink(values=aniso_values, space=interior)` is an `AngularSourceSink` whose values equal an independent reference, `np.broadcast_to((iso_values / math.fsum(weights))[None], interior.shape) + aniso_values`, bit for bit (the typed injection is the section followed by a within-class add; the reference uses `math.fsum`, the platform-free definition, not another reduction of the weights; `weights` is `posing_quadrature(mesh, 4).weights`); a within-class `AngularSourceSink + AngularSourceSink` is unchanged. (3) A retirement leg: `AngularSourceSink` and `BulkField` carry no attribute `CONTAINS` and `AngularSourceSink.__radd__` is `Field`'s (`'__radd__' not in vars(AngularSourceSink)`), asserted by the engine rows of E6, not duplicated here. **First red** `[M]` (`p_refused_input.py`): today both orders construct and return an `AngularSourceSink`, with and without `-O`; this is the constructible input the refusal rejects. **Mutations**: (i) the `-O` leg: copy `orpheus` to a temporary directory, replace the `raise TypeError(` that follows `if not self._is_same_class(other):` in `Field._check_partner` by `assert False, (` in the copy of `numerics/field.py`, and run the leg-(1) child with `sys.path.insert(0, <copy>)` as its first statement (a `python -c` child puts its working directory at `sys.path[0]`, ahead of `PYTHONPATH`, so `PYTHONPATH` alone imports the main tree: `code-search`), printing the child's `orpheus.__file__` first and asserting it is under the copy; `[M]` under `-O` the `assert` is stripped, the sum then raises `ValueError: ... arithmetic requires equal space` from the NEXT check of `_check_partner` (a second refusal, not a construction), its stderr has no `same-class partner`, and the fragment assertion of leg (1) reddens (this is why the leg asserts the fragment and not only the exit status); (ii) restore the shipped arm in a copy of `angular_source_sink.py` (`CONTAINS = (ScalarSourceSink,)` and the shipped `__add__`): `[M]` the shipped code constructs both sums (5.7.2 item 1), so leg (1) reddens, and the engine rows of E6 redden.

**The error-catalogue and ledger legs.** No `ELEGANCE-DEBT` tag is written and `tests/gates/test_elegance_debt_is_tagged.py` is not edited by this sub-unit (no guard lands: the refusal is `Field`'s existing gate). No ERR entry is owed: no shipped number moves beyond its solve tolerance (the memo's verdict, 5.7.7).

#### 5.7.6 The 8 tests of `test_typed_source_sinks.py`, and the one in `test_n2n_operator.py`

`[M]` Run in process, `-O`, with the E1 rewrite in place and the mixed sum refused, in both variants of 5.7.2 item 7 (`fix708.py` with its own `TypeError`, and `fix708b.py` with the five names deleted: the same 8 tests fail, the same 28 pass): exactly 8 of the 36 tests of `tests/gates/sn/primitives/test_typed_source_sinks.py` fail with `TypeError` and 28 pass: 5 in `TestCrossClassDunder` (the first five below), 2 in `TestBroadcastInjection` and 1 in `Test2DTypedSources`. Fate of each (the file's own fixtures `_slab_mesh`, `_stretched_mesh`, `_2d_mesh` stay; the weights sum in each is `math.fsum(m.quad.weights)`, equal to the `np.sum` value on GL4 and on `level_symmetric(4)` `[M]` `p_rewire.py`, `p_rewire3.py`):
1. `TestCrossClassDunder::test_iso_plus_aniso_returns_per_ordinate`: REWIRED. It pins the un-normalised broadcast that the retired arm performed. It becomes: `combined = AngularSourceSink.from_isotropic(iso_values, m) + aniso`; `isinstance(combined, AngularSourceSink)`; `combined.space is m.angular_bulk_space`; values `np.broadcast_to((iso_values / math.fsum(m.quad.weights))[None], aniso_values.shape) + aniso_values`, `assert_array_equal`. `[M]` bit-equal in process.
2. `::test_aniso_plus_iso_commutative`: DELETED. Its subject is the commutativity of the retired arm; the refusal in both orders is G3 leg (1), which also holds the `-O` leg.
3. `::test_iso_plus_aniso_coherence_is_the_marginal`: REWIRED, same intent (coherence is the space relation). Equal leg: `AngularSourceSink.from_isotropic(np.ones(...), m1) + AngularSourceSink(values=..., space=m2.angular_bulk_space)` with `m2` a twin carrier is legal `[M]` (`p_rewire3.py`); unequal leg: the same sum with `aniso3` on `_stretched_mesh()` raises `ValueError` matching `requires equal space` `[M]` (message `AngularSourceSink arithmetic requires equal space; ...`, raised by `Field._check_partner` after the class check).
4. `::test_pattern_7_normalised_iso_plus_aniso`: REWIRED. `(iso / sum_w) + aniso` is refused; the surviving spelling is the typed injection of 1. The reference is built independently of the section: `np.broadcast_to((iso_values / math.fsum(m.quad.weights))[None], ...).copy()` plus `aniso_values`. The docstring's "mirrors the ScatteringOperator call site at scattering.py:942" is stale (that line no longer exists) and is dropped.
5. `::test_per_ordinate_from_isotropic_alternative_factory`: REWIRED. It compared the factory with the dunder form `(iso / W) + zeros`; it becomes a comparison of `from_isotropic(iso.values, m).values` with the independent `np.broadcast_to((iso_values / math.fsum(m.quad.weights))[None], shape)`, `assert_array_equal` (the section against a platform-free reference, not against another reduction of the same weights: `vv-testing`, X4).
6. `TestBroadcastInjection::test_dunder_broadcast_shape_and_values`: DELETED. It pins the un-normalised broadcast `iso + zeros`, which is the defect; the shape and per-ordinate-slice facts it asserted are covered by 5 (equal values across ordinates) and 1.
7. `::test_dunder_broadcast_owns_data`: REWIRED. It asserted that the dunder's result is writable. It becomes `(AngularSourceSink.from_isotropic(...) + AngularSourceSink.zeros(m.angular_bulk_space)).values.flags.writeable` `[M]` True in process (the `from_isotropic` values alone are a broadcast result with `flags.owndata` False and `writeable` True; the sum allocates).
8. `Test2DTypedSources::test_cross_type_add_2d`: REWIRED to the 2-D mesh: `from_isotropic(iv, m) + aniso` against `np.broadcast_to((iv / math.fsum(m.quad.weights))[None], aniso.values.shape) + aniso.values` `[M]` `array_equal` True (`level_symmetric(4)`, 24 ordinates, `fsum == np.sum == 12.566370614359172`).
Also in that file: the section comment "Cross-class arithmetic via canonical subspace containment (refined Issue #207 principle)" (lines 196 to 202) and the class docstring of `TestBroadcastInjection` are rewritten to say what the class now tests; the module docstring is read for tense. `[M]` `grep -n -E 'verifies\(|catches\(' tests/gates/sn/primitives/test_typed_source_sinks.py tests/gates/sn/operators/test_n2n_operator.py` prints 0 lines (both files carry only a module-level `pytestmark = pytest.mark.foundation` in `test_typed_source_sinks.py`), so the three deleted tests migrate no marker (C.13). The docstring of `TestHarmonicMomentSourceSink::test_cross_class_with_scalar_source_rejected` (lines 529 to 532, present tense: "unlike the Scalar↔Angular pair, the moment leaf carries no subspace-containment dunder") and the module docstring of the file (lines 1 to 10: "the load-bearing cross-type" dunder algebra, "the bit-identity of the new dunder path") are rewritten; the assertion of the former, a `TypeError` matching `same-class partner`, is untouched. `docs/theory/conventions/indexing_and_layout.rst` line 1961 says "(37 cases)": the file has 36 tests `[M]` (`.venv/bin/python -O -m pytest -p no:cacheprovider --collect-only -q tests/gates/sn/primitives/test_typed_source_sinks.py` before the edit; `docs/theory/verification/matrix.rst` line 677 prints 36), and after the edit 34 (2 deleted): the builder writes the count the collection prints.
`tests/gates/sn/operators/test_n2n_operator.py::test_total_weight_is_the_measure_mass` (line 419): DELETED. It compares `n2n.total_weight` with `float(sn.quad.weights.sum())` under `np.testing.assert_allclose` (default `rtol=1e-7`) on GL4, where `np.sum` and `fsum` agree (`2.0`): green before and after the change, so it pins nothing about the divisor (X1), and its subject is G1's first leg, which has a first red. `[M]` `git grep -n 'test_total_weight_is_the_measure_mass' -- . ':!scratch'` is this one definition. The file is also edited by dead-D (row 122's consumer; 7.2): the two edits are disjoint tests. `tests/gates/transport/test_angular_lift.py` line 224 asserts `F.total_weight == S.total_weight == 2.0` on GL8, where both spellings are `2.0`: it stays green and unchanged.

#### 5.7.7 What stays bit-identical, what moves, and the builder's serial run

- **Bit-identical, by construction `[R]` and measured on GL8** (5.7.2 item 9): every number on a rule where `float(np.sum(w)) == math.fsum(w)` (91 of the 120 rules `[M]`: 120 minus the 29 of 5.7.2 item 2): `total_weight`, `_combine`, the transpose and the transfers return the same bits. The section, DSA, `from_isotropic` and the retraction are unchanged.
- **Moves, at the 1e-16 level, on every rule where the two sums differ** (the 29 of the 120-rule population of 5.7.2 item 2, and others outside it: `folded_product(2, 6)` `[M]`, the rule `test_quotient_adjoint_sn.py` and `test_dual_metric_sn.py` build): the divisor moves by at most 4.7e-16 relative (about one ULP of the mass, of either sign: `np.sum` is below `fsum` on GL7 and GL12 and above it on `folded(2,12)` and `folded(4,12)`, `[M]` from the scales of `p_radial.py`). Every emitted source, transpose and transfer built on such a rule moves by that factor. The memo's measurement, with the swap of `total_weight` to `fsum` (state A) applied in process over 356 test files of `tests/gates/sn`, `tests/gates/transport` and `tests/gates/cross_method` (6 493 tests, one file per process, serial): 59 tests in 19 files read a divisor that differs and 0 of 6 493 reddened; one fixed-source and one eigenvalue solve on GL12 moved `k` by at most 1.8e-16 relative and the flux by at most 5.6e-15, below the solve tolerance 1e-13 (`scratch/ontology/w0/i708/memo.md`, sections 3 and 4; NOT re-run by this specification: I did not reproduce the count 59, the 356-file run or the solve table). Measured here, with E1 and the refusal in process (`fix708.py`, the variant that kept the arm; the ruled form's run is the next sentence): 343 tests in the 19 reached files plus `tests/gates/transport/test_angular_lift.py` and `tests/gates/sn/operators/test_n2n_operator.py`, serial, `python -O`: 343 passed, 0 failed (about 13 minutes); and a second in-process run of 5 of those same files (115 tests) with an equality witness on every `total_weight` read (`eqcheck708.py`): 23 516 of 23 516 reads have `math.fsum(self.frame.measure.weights) == AxisSectionOperator.of(...).mass`, over lifts of 4, 6, 8, 12, 16, 32, 64 and 110 ordinates (the fix reads the section, the old body read the frame's measure: they are one array on every lift the suite builds). The ruled form, `fix708b.py`: 343 passed, 0 failed over the same 21 files (`python -O`, 671 s), with the deletion of the five names, E1 and E5 all in process; plus 255 passed over `tests/gates/transport/fields`, `tests/gates/numerics/test_field.py` and the G1 draft, and `test_typed_source_sinks.py` with its 8 expected failures and 28 passes. The G3 refusal leg, a subprocess, cannot see an in-process patch, so it was run on a copy of the tree with the real text edits of E2 and E4 (5.7.5).
- **Pinned comparisons that are bit-exact and therefore must be run by the builder.** A comparison pinned with `assert_array_equal` or a `nulp` band of a source, transpose or transfer on a rule of the 29 is the only kind that can redden; none did in the 343 tests above and none in the memo's 6 493, the memo found 10 tests whose tolerance is tight enough to see a 1e-9 perturbation but admits the 4e-16 one (the three `TestAnisoMomentSourcePath::*bit_identical_to_pre_t3_snapshot` of `tests/gates/sn/operators/test_scattering_operator.py`, `test_apply_isotropic_flux_p0_only`, `test_typed_apply_returns_per_ordinate_already_normalised`, `test_fission_adjoint.py::TestCompositeTransposeArm` [2g] and [4g], and three `test_cylinder_three_way_standoff[6-2-*]`), so those are the ones to read first if the run reddens. The regression snapshot `tests/gates/sn/regression/test_dd_regression.py` reached 0 tests under the swap (none of its cases uses a rule of the 29).
- **Radial scale (E5): which pinned values under `tests/gates/sn` can move** ``[M]` (`scratch/ontology/w0/spec/amend/p708/radial708.py`, a pytest plugin that re-assigns `_ell0_scale = 1.0 / math.fsum(weights)` after the constructor and records the tests where it differs from the shipped value; one file per process, serial, `python -O`, `RADIAL708=fsum`): over 22 files of `tests/gates/sn` that build a cylinder with `Quadrature.folded_product` (chosen by reading `git grep -n 'folded_product(' -- tests/gates/sn`: the files whose rule can be one of the differing folded rules, `(2,12)`, `(4,12)`, `(4,16)`, `(6,12)` of the population and also `(2,6)`, which is outside the 120-rule population and has `np.sum != fsum` `[M]` (`(4,4)`, `(4,6)`, `(2,4)`, `(8,8)`, `(8,16)` agree), added with `test_quotient_adjoint_sn.py`, `test_g_adjoint_reciprocity.py`, `test_dual_metric_sn.py`; the 19 first are `test_unified_matvec_cylinder.py`, `test_pure_L_sigma_free.py`, `test_streaming_operator_decomposition.py`, `test_bc_extraction_matvec.py`, `test_apply_matvec_cylinder_invariants.py`, `test_angular_endpoint_defect.py`, `test_psi_half_positivity.py`, `test_b1pp_verification.py`, `test_affine_carve_baseline.py`, the radial-characteristic mesh and operator files, `tests/gates/sn/regression/test_dd_regression.py`): 22 of 22 files exit 0 with the swap (427 tests); 188 arc-family reconstructions were constructed and 21 tests in 5 files built one on a rule whose `1/np.sum` differs from `1/fsum` (`test_angular_endpoint_defect.py` 2, `test_apply_matvec_cylinder_invariants.py` 3, `test_psi_half_positivity.py` 2, `test_quotient_adjoint_sn.py` 10, `test_dual_metric_sn.py` 4; the scales differ by one or two ULP, about 1.4e-17 to 2.8e-17 absolute at 0.0796); 0 of the 427 tests turned red. The population is the 22 files, not `tests/gates/sn` whole: an attempt to run the directory whole with the plugin did not finish in an hour and was stopped, so a cylinder test that reaches the scale through a rule built in a fixture file I did not list is not covered, and the builder's serial run of `tests/gates/sn` (below) is the measurement for the rest. The snapshot `test_dd_regression.py` reached no differing scale. The rules `(4, 8)` and `(4, 6)` that most cylinder fixtures use agree under both spellings `[M]`.`
- **The builder's serial run**, `.venv/bin/python -O -m pytest -p no:cacheprovider <dir>` one directory at a time, exit 0 each: `tests/gates/transport`, `tests/gates/sn`, `tests/gates/cross_method` (the three directories the memo swapped), then `tests/gates/homogeneous`, `tests/gates/diffusion`, `tests/gates/derivations`, `tests/gates/cp` and the root files `tests/gates/test_*.py`, and `tests/gates/numerics` (the section's own gates and `test_field.py`, which pin `same-class partner`). The list is derived by an `ast` import census over `tests/gates` (`scratch/ontology/w0/spec/amend/p708_dirs.py`; positive control: `transport` must count) of the files that import `orpheus.transport.operators`, `orpheus.transport.fields` or `orpheus.sn` directly, with its denominators: `sn` 267 of 312 files, `transport` 33 of 90, `cross_method` 14 of 25, `homogeneous` 6 of 15, root 6 of 27, `diffusion` 4 of 18, `derivations` 3 of 128, `cp` 1 of 21, `numerics` 0 of 144 (the section's gates are named for another reason), 0 in the other directories; a file that reaches a lift through a helper module it imports is not counted (`[R]`), which is why the three large directories are run whole. `nexus impact AngularLift.total_weight` printed "No upstream dependents found" (the graph does not link the property's readers), so the graph did not decide it. The three new gates and the rewired files first. Run the three new gates and the rewired files first, then `tests/gates/test_layer_imports.py`, `tests/gates/test_test_layers.py`, `tests/gates/test_docstring_xrefs.py`, `tests/gates/test_retirements.py`.

#### 5.7.8 Docs and docstrings made false, with the edit of each (retirement audit E.19 and D.15)

Sort by tense: a present-tense claim that the change makes false is rewritten; a past-tense history line stays and gains a dated note.
- `operator_algebra.rst` lines 5736 to 5742 (present): "The only leaf on the grid that overrides an additive dunder is AngularSourceSink ... it accepts a ScalarSourceSink partner": rewrite to the new fact `[M]`: no leaf on the grid overrides `__add__` or `__sub__` (`git grep -n -E 'def __(add|radd|sub)__' -- orpheus/transport/fields` prints 0 lines after E2 and E4), and #288 is closed by the retirement.
- `operator_algebra.rst` lines 5793 to 5851 (present, then history): the bullet "The `iso + aniso -> AngularSourceSink` source injection is a hand-rolled Representation traversal inside a dunder, and it is endorsed ... so a future reader does not 'tidy it away'". **This is a REPEAL of a recorded ruling** (the endorsement "at creation", the refined #207 principle). Rewrite the bullet in the past tense: the dunder was retired on 2026-11-20 by #708 for the structural reason that it was the pullback `pi^*` spelled outside its type with no division by the mass, it summed two dimensionally different rates, and its only production caller, `_combine`, already divided by the mass first; the injection is the section `E`, `from_isotropic` its producer entry, and `AngularLift._combine` the one combine; keep the ⛔ record of `assemble_per_ordinate_isotropic` (lines 5821 to 5832, past tense, true) and the ⚠ paragraph (lines 5834 to 5851) as dated history. Section 8.4 item 8.29 is the orchestrator-side record of the repeal (as 8.27 is for the role ruling).
- `operator_algebra.rst` lines 5880 to 5890, the table cell (a) of "the five obstructions" (present): its "Two facts make Role's interface role-dependent. First, the Source role owns a cross-representation `__add__`". After E2 that first fact is false; the second fact (role identity is units identity, read by `Field._check_partner` as `type(self) is type(other)`, erased under a phantom parameter) is untouched and still carries the conclusion that Role must be a class. Edit: state that the first fact held until 2026-11-20 and no longer does, that the obstruction stands on the second alone, and read the cell's closing paragraph ("Until 2026-08-19 this cell argued...") for the number of facts it names; the same note in `field_algebra.rst` lines 1136 to 1146 ("AngularSourceSink **does** define its own `__add__`", "the axis that 'changes the arithmetic interface' is now the *source* axis"): rewrite to say it did, until 2026-11-20, and that no leaf now does, so no axis changes the arithmetic interface and the conclusion rests on units identity.
- `operator_algebra.rst` lines 5116 to 5122 and 5142, and `cartesian_multid.rst` line 3144, `adjoint.rst` lines 1509 to 1512 (present): "`W` is read off the retained faces' frame measure ... not off a stored weight vector", "`total_weight` derives `W` from the retained frame's measure": rewrite to "`W` is the angular section's mass on the emitted interior, the correctly rounded sum of the stored weights (`AxisSectionOperator.mass`)".
- `spaces.rst` lines 5961 to 5971 and 6465 to 6472 (present): the `implements::` bullet whose third site is "a hand-written twin on the angular axis ... `AngularSourceSink.__add__` with an isotropic partner adds `other.values[None]`", with the directive `.. implements:: spaces-collapse-adjoint-is-pullback :by: py:method:...AngularSourceSink.__add__`, and the warning "a hand broadcast survives": the target no longer exists, so the directive is a dead reference that `nexus dead-references` reports and `sphinx-build -W -n` may too: delete the directive, change "Implemented by 3 sites" to the true count and record the twin as retired with the date. `spaces.rst` lines 7043 to 7080 (the divisor table): add the row `AngularLift.total_weight` (the lift's divisor until #708; `w.sum()`; differs on 29 of 120) to the list of spellings and a dated sentence that every divisor now reads the section.
- `layering.rst` lines 163 to 166 (present): the bullet "An embedding written in the smaller class ... (`AngularSourceSink` answers `iso + aniso`; the smaller class defers through `BulkField.CONTAINS`)": the example is retired; rewrite the bullet as a past-tense record of one of the "four shapes" the cycle repairs took, saying the mechanism was retired on 2026-11-20 (#708) and naming no `CONTAINS`. `conceptual_view.rst` lines 685 to 687 (present, "adds a hand-written broadcast, pi^* outside its type"): remove `AngularSourceSink.__add__` from the list of hand-rolled sites and say it retired with #708; the debt list is recomputed.
- `indexing_and_layout.rst`: line 901 ("the ONE shared `(iso / W) + aniso` combine", present) becomes "the ONE shared combine: the isotropic part through the angular section, plus the anisotropic part"; line 1944 opens a `.. todo::` item whose bullet at line 1948 is listing what the archivist should write ("The math of why iso + aniso = per-ordinate", "How the typed dunder dissolves ...") and lines 1966 to 1995 are the section "The load-bearing dunder", present tense: rewrite the section as the history of a retired dunder (what it replaced, why it was retired), and drop the two todo bullets that name it; line 1961 count (5.7.6).
- `orpheus/transport/fields/__init__.py` lines 47 to 53 (present): "Cross-class same-units operations (e.g. combining an isotropic source/sink with a per-ordinate one) go through the canonical subspace-containment injection or a NAMED composition factory, not an ad-hoc cross-class dunder": rewrite to "through the angular section (`AngularSourceSink.from_isotropic`) or a named composition factory". `operator_algebra.rst` lines 5015 to 5019 name `TestHarmonicMomentSourceSink` of the test file and stay true; `tests/gates/transport/fields/test_flux_vector_algebra.py` line 312 names the contrast row of `test_typed_source_sinks` as history and stays.
- Docstrings promoted or made false (D.15, B.8): `AngularLift.total_weight` (E1); `AngularLift._combine` (E1); `AngularSourceSink` module, `__add__` (deleted), `from_isotropic` (E2); `ScalarSourceSink` module and comment (E3); `BulkField` has a docstring paragraph at lines 355 to 362 ("Abstract ...") read for a mention; `docs/architecture/layering.rst` as above.
- **Not edited** (checked, `[M]` by the census of 5.7.2): the 54 reference rows under `orpheus/derivations/`; `orpheus/sn/acceleration/dsa.py` (its comment at lines 583 to 590 already says the divisor is the section's mass); `docs/theory/foundations/spaces.rst` lines 6099 to 6110 and 6885 to 6895 (they state the section's divisor is `math.fsum`, true before and after).
- Issues: #708 and #288 are closed by hand with the integrating hash (the closing comment of #288 says the dunder it asked to decide on is retired, and why); `#556` and `#558` are named in #708's body and are not closed.

#### 5.7.9 Retirement-audit items and done-when

A.1 `.venv/bin/nexus callers AngularSourceSink.__add__` after a docs build (the graph is rebuilt by it) must print no caller; A.2 the text surfaces of 5.7.8; A.3 direct constructors: none new; B.4 `nexus dead-references` over `AngularSourceSink.__add__`, `AngularSourceSink.__radd__`, `BulkField.CONTAINS`, `BulkField.__add__` (the `implements::` target of `spaces.rst`); B.6 strings: `git grep -n -E "['\"](CONTAINS|__radd__)['\"]" -- orpheus tests tools` (expected 0; control: a string constant you know, `"__add__"` in `tests/gates/numerics/test_field.py`); B.11 the fragment `same-class partner` is untouched and pinned (5.7.2 item 8); C.12 the 8 and 1 tests of 5.7.6; C.13 `git grep -n -E 'verifies\(|catches\(' -- tests/gates/sn/primitives/test_typed_source_sinks.py tests/gates/sn/operators/test_n2n_operator.py` prints 0 lines (positive control: the same command on `tests/gates/numerics/test_axis_marginal_lift.py` prints the `verifies` markers it carries, or any known marked file); D.14 G1's mutations; D.16 the gate that compared `total_weight` with `np.sum` is deleted with its reason; D.18 the mechanism that kept `total_weight` and `AxisSectionOperator.mass` equal was nothing (they were two computations that agreed except on 29 rules); the witness is G1, written in the same commit that retires the second computation.

Done when (each one command):
1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/transport/test_angular_lift_divides_by_the_sections_mass.py tests/gates/sn/operators/test_radial_ell0_scale_is_the_correctly_rounded_mass.py tests/gates/transport/fields/test_scalar_source_does_not_add_to_an_angular_one.py tests/gates/sn/primitives/test_typed_source_sinks.py tests/gates/sn/operators/test_n2n_operator.py tests/gates/test_retirements.py` exits 0.
2. The census probe, a file written under a quoted heredoc (`<<'PY'`) that walks `orpheus/` (excluding `orpheus/derivations/`, stated) with `ast` and prints every call to `sum`, `np.sum`, `numpy.sum`, `.sum()` or `math.fsum` whose argument or receiver names something containing `weight` or equal to `w` or `wts`, and every load of `total_weight`, `_sum_w`, `sum_w` or `w_norm`: `[M]` today it prints 22 rows outside `derivations/` (5.7.2 item 5); after the sub-unit the rows that are divisors of the angular mass are `angular_lift.py` (the property no longer sums and `_combine` no longer reads `total_weight`: 0 `SUM` rows and 1 `READ` row, the transpose), `radial_characteristic.py` (one `SUM` row, now `math.fsum`), `dsa.py` (unchanged: 3 `READ` rows of `_sum_w` and the two tolerance guards), `transfer.py` (2 `READ` rows); the `np.sum` rows that remain outside `derivations/` are the tolerance guards, the Gram norms, `realized_operators.py` 522 (8.1) and the `mc` rows. Positive control: the probe asserts it finds, on the untouched `dsa.py`, the `READ` of `_sum_w`.
3. The `CONTAINS` census, the same shape: `ast` over `orpheus`, `tests`, `tools`, `examples`, `student_resources` printing every `Name`, `Attribute` and string constant `CONTAINS`; after E2 and E4 it prints 0 rows; positive control: the same probe on a clean checkout of the base commit `b68f8b8c` (a separate worktree, never a stash) prints 3 rows.
4. `git grep -n -E 'def __(add|radd)__' -- orpheus/transport/fields` prints 0 lines; `.venv/bin/python -O -c "from orpheus.transport.fields import AngularSourceSink, ScalarSourceSink"` plus G3 leg (1) shows the refusal.
5. The serial run of 5.7.7 exits 0 per directory; `sphinx-build -E -W docs docs/_build/html` is clean and `.venv/bin/nexus dead-references` lists none of the four names of A.1.

## 6. Unit `posing`: items u1.2, u1.3 and u0.3 only

Binding decision 6: of posing units 0 and 1, W0 keeps u1.2 (the evidence readers become exhaustive), u1.3 (the complex gauge reading) and u0.3 (the `generator_as` residue). u1.6 (the block-role comment) is inside `roles` (section 4). Everything else is a later wave's and is not touched: u0 items 1, 2, 5, 6, 7, 9 (renames in W4/W5), 8 and 10 (W2), 4 (the withdrawal family: 51 of 94 lines are inside `derivations/`, closed by A2, and a rename of a type `derivations/` owns is not "mechanical re-pointing"), u1 items 1 (W4), 4, 5 and 7 (W5), 8 (W1 and W2).

### 6.1 File set

Hard: `orpheus/numerics/convergence.py` (u1.2, the `balance` clause of `warn_if_unconverged`), `orpheus/sn/operators/loss_kernel_gauge.py` (u1.2, `warn_if_gauge_freedom`), `orpheus/numerics/gauge.py` (u1.3, `ScaleGauge.displacement`), `tests/gates/numerics/test_gauge.py` (u1.3), `docs/theory/foundations/spaces.rst` and `docs/theory/methods/sn/history.rst` (u0.3, one present-tense line each, 6.4).
New: `tests/gates/numerics/test_convergence_evidence_reader_is_exhaustive.py` (u1.2, the `convergence.py` reader), `tests/gates/sn/operators/test_gauge_evidence_reader_is_exhaustive.py` (u1.2, the `loss_kernel_gauge.py` reader: a numerics-directory test may not import `orpheus.sn`, so the second reader's gate lives under `tests/gates/sn/`; the measurement stage listed one new file, the layer gate needs two: decision 13).
Read-only consumers (listed so that nobody edits them): `orpheus/numerics/evidence.py` (the type), `orpheus/sn/solver.py` and `orpheus/homogeneous/solver.py` (`ScaleGauge.displacement` is reached in production through `ScaleGauge.apply` at `sn/solver.py` lines 2471 and 2744; `orpheus/homogeneous/solver.py` line 563 is a third; the `.displacement(` calls at `sn/solver.py:811` and `diffusion/solver.py:641` are `KernelProjectionGauge.displacement`, a homonym), `tests/gates/numerics/test_axis_identity.py` (a live negative pin: `assert not hasattr(Axis, "generator_as")`), the seven test files that name `ScaleGauge` (6.3), `tests/gates/numerics/test_matrix_carrier_membership.py` (it names the homonym `KernelProjectionGauge.direct_sum`, not `ScaleGauge`: not a consumer of u1.3), `tests/gates/numerics/test_axis_generator.py` (a docstring that names `generator_as` in the past tense: "it replaced G5_6b's `generator_as` call-site census"; stays). Not edited by this unit although the measurement stage listed it: `docs/theory/foundations/operator_algebra.rst` (the one line naming `generator_as`, "S7a: `Axis.generator` and its narrowing accessor `generator_as` became `Axis.origin`", is history; it stays), so `posing` no longer shares that page with `dead` and `roles`.

### 6.2 u1.2: the evidence readers

`Evidence = Measured | Asserted | NotApplicable | NotYet` (`orpheus/numerics/evidence.py`). A reader that tests `isinstance(x, Measured)` and treats every other value alike is silent on a fifth arm; the union is closed, so a reader is an exhaustive `match` that ends in `assert_never`, which is the idiom the tree declares for closed sums (`orpheus/transport/boundary_conditions/network.py` states the law: "every place that reads a chain does an exhaustive `match` ending in `assert_never` ... a third chain is an edit of the alias, after which every such match reddens"; `[M]` `git grep -n assert_never -- orpheus` finds it in 20 production files, `cp/solver.py:40, 590, 599` among them). The two `match` statements of `orpheus/reference/verification.py` are NOT the pattern to copy: there the `case _` arm catches the legitimate arms `NotApplicable` and `NotYet` and raises `Unestablished`, a domain verdict, not a refusal of a non-arm. The static law is the point: `[M]` (the elegance check's pyright run of 2026-11-20 on a probe) the four-arm reader ending `case _: assert_never(e)` gives 0 errors, the same reader over `Evidence | Fifth` gives `"Fifth" is not assignable to type "Never"`, and `reader(None)` is refused at the call site; `orpheus/numerics` and `orpheus/sn` are inside pyright's focus set (`pyproject.toml` `ignore` omits them) and `tests/gates/test_pyright_ratchet.py` reads it. `typing.assert_never` raises `AssertionError`, which survives `-O`. `[M]` `m_pos1.py`: at both sites every non-`Measured` value falls through silently: for `Asserted`, `NotApplicable`, `NotYet`, a foreign class `NewKind()` and `None`, site 1 produces the same message as `NotApplicable` (no clause, no raise) and site 2 emits 0 warnings against 1 for `Measured(0.5)` (the gauge verdict stubbed to "present"). Three of the five are LEGITIMATE silent arms (`Asserted`, `NotApplicable`, `NotYet`: the clause says nothing rather than "unavailable"); the red inputs are `NewKind()` and `None`, outside the union.

Site 1, `warn_if_unconverged` (signature default `balance_defect: Evidence = NotApplicable("no balance defect is measured at this entry")`; the clause today is `"" if not isinstance(balance_defect, Measured) else f"The returned iterate leaves a per-group balance defect of ‖R_g‖/‖Q_g‖ = ..."`). Replace by:

```python
match balance_defect:
    case Measured(value=value):
        balance = f"The returned iterate leaves a per-group balance defect of ‖R_g‖/‖Q_g‖ = {value:.3e} — a DIAGNOSTIC ... "
    case Asserted() | NotApplicable() | NotYet():
        balance = ""        # absent says nothing: an empty clause cannot be misread as a measurement
    case _:
        assert_never(balance_defect)
```
(keep the whole existing sentence text and the comment block above it; only the dispatch changes; the `value` name replaces `balance_defect.value`; import `Asserted`, `NotApplicable`, `NotYet` and `assert_never` if absent). Site 2, `warn_if_gauge_freedom(problem, correction, ...)` in `orpheus/sn/operators/loss_kernel_gauge.py`: the tail reads `if not verdict.present or not isinstance(correction, Measured): return` then `if correction.value <= _GAUGE_AUDIBLE_FLOOR: return`, and the `if verdict.undetermined:` branch warns and returns BEFORE `correction` is read. Restructure so that the `match` is the FIRST statement of the function body (before `verdict = gauge_freedom(problem)`):

```python
match correction:
    case Measured(value=value):
        reading: float | None = value
    case Asserted() | NotApplicable() | NotYet():
        reading = None
    case _:
        assert_never(correction)
verdict = gauge_freedom(problem)
... the `if verdict.undetermined:` branch, unchanged ...
if not verdict.present or reading is None:
    return
if reading <= _GAUGE_AUDIBLE_FLOOR:
    return
```
and the warning's `{correction.value:.2%}` reads `{reading:.2%}`. Why first: `[M]` (structure check) the undetermined branch returns before `correction` is read, so under any order that puts the `match` after it, `NewKind()` with `undetermined=True` warns and returns with no raise, and every leg whose stub fixes `undetermined=False` stays green; the first statement refuses a foreign kind whichever the verdict, and the four legitimate arms behave as before (a legitimate value produces the same message or the same silence). Behaviour changed: a fifth arm and `None` raise `AssertionError`. `[M]` `git grep -n "balance_defect=\|warn_if_gauge_freedom(" -- orpheus` lists 7 `balance_defect=` sites in 4 files (`cp/solver.py`, `diffusion/solver.py`, `moc/solver.py` pass a literal `NotYet(485, ...)`; `sn/solver.py` five sites pass `exit_report.balance` or `solution.exit_report.balance`) and 5 `warn_if_gauge_freedom` call sites in `sn/solver.py`, each passing `exit_report.gauge` or `solution.exit_report.gauge`; `ExitReport.balance` and `.gauge` are required fields typed `Evidence` (`orpheus/numerics/outcome.py`, the class `ExitReport`), so no call site passes `None` by type, and `cp` and `moc` are pyright-ignored. **No shipped input is rejected by this item: every first red below is a test-only input** (`NewKind` is a test class, `None` is unspellable at every checked caller); the shipped content of the item is the static law (a fifth arm reddens pyright at both readers).

Gates (new files above). Each: the stalled-record fixture of `m_pos1.py` for site 1 (`IterationRecord(label="inner(probe)", criteria=(StoppingCriterion(name="residual", trajectory=tuple(1e-3 for _ in range(40)), tolerance=1e-10),), budget=IterationBudget(40, "max_inner"), iterations_run=40)` and `warn_if_unconverged(stalled, where="probe", balance_defect=ev)` under `warnings.catch_warnings(record=True)`; `NewKind` is a test-only class; the legitimate arms must produce the same message as `NotApplicable("x")`; `Measured(1e-3)` must contain `balance defect`; `NewKind()` and `None` must raise `AssertionError`); for site 2, `K.gauge_freedom` monkeypatched to `lambda p: types.SimpleNamespace(undetermined=False, present=True, because="stub")` and ALSO to a stub with `undetermined=True`, `problem = types.SimpleNamespace(ndim=2)`, `K.warn_if_gauge_freedom(problem, ev, where="probe")`: with `undetermined=False`, `Measured(0.5)` emits 1 warning, `Asserted(1e-9, "t")`, `NotApplicable("x")`, `NotYet(1, "r")` emit 0, `NewKind()` and `None` raise; with `undetermined=True`, `NewKind()` and `None` raise as well (the leg the earlier order fails). Each gate also asserts `set(ARM_FIXTURES) == set(typing.get_args(Evidence))`, where `ARM_FIXTURES` maps each arm to one instance, so that a fifth arm added to the union reddens the gate until a fixture is added. **First red** (a test-only input, stated above): `NewKind()` and `None` are silent today at both sites, and `NewKind()` with `undetermined=True` is silent at site 2. **Mutations**: restore the `isinstance(..., Measured)` form in a copy of the module (a temporary root, the module loaded with `importlib.util.spec_from_file_location`, or an in-process `monkeypatch` of the function with the old body): `NewKind()` is silent again, red; move the site-2 `match` after the `undetermined` branch in the copy: the `undetermined=True` leg reddens; and the static control: add a fifth arm to a temporary copy of `evidence.py` (`class Fifth` joined to the union) and run `.venv/bin/pyright` on the two edited modules with that copy first on the import path: it reports an error at BOTH readers (`"Fifth" is not assignable to type "Never"`), which the unedited tree does not. The refusals are `assert_never` calls (they survive `-O`; the gates run under `-O`). Whether this is "guard" debt: no. A closed union read exhaustively is the law of the type, checked statically, not a guard against a mistake the types could make unspellable (they already make it unspellable); no `ELEGANCE-DEBT` tag. Callers of the two functions and every file that names `Measured` (57 lines in 16 production files; only these two branch with a silent else, `[M]` `git grep -P` for `isinstance(… Measured|Asserted|NotApplicable|NotYet)`, `type(…) is <arm>` and `case <arm>` over `orpheus/`: 3 `isinstance` hits, the third being `Measured.__truediv__` (`evidence.py`; `Measured` defines no `__eq__`); 6 `case` lines, all in `reference/verification.py`) are not edited. The reader MODULES are redesigned later (`convergence.py` by W5's "Iteration record", `loss_kernel_gauge.py` by W2's "Nullspace"); this unit's edit is local.

### 6.3 u1.3: the complex gauge reading

`ScaleGauge.displacement` (`orpheus/numerics/gauge.py`) is `n = float(self.functional(state))` then the zero test and `return self.target / n`. `[M]` (`m_pos2.py`, `-O`): `ScaleGauge(functional=lambda s: np.complex128(1+2j), target=1.0).displacement(None)` returns `1.0` under a `ComplexWarning` (the true value is `0.2-0.4j`: the imaginary part is discarded); a Python `complex` and a 0-d complex array already raise `TypeError` (by accident, inside `float`); `np.complex128(1+0j)` returns `1.0` and also warns `ComplexWarning`. The repository does not turn the warning into an error (`pyproject.toml` `[tool.pytest.ini_options]` has no `filterwarnings`, read lines 99 to 120).
Law (stated whole): the group of the scale section is the multiplicative group of NONZERO REALS under multiplication (`ScaleGauge`'s own module states an `R_+` torsor with a `Z/2` leg, `orpheus/numerics/gauge.py`, the module docstring and the class docstring): the section `t / n(psi)` is defined where the functional reads a nonzero real number, and is undefined at 0, which is the zero refusal that already exists. A reading is a real number; a reading of any complex dtype (`np.complex128`, a Python `complex`, a complex 0-d array), whatever its imaginary part, lies outside the domain of the section and is refused, `ValueError`, fragment `complex reading`, naming the functional. Nothing is widened: no complex carrier is admitted, and the real part is never taken. Then the existing laws: zero reading refused (`ValueError`, the existing message naming the functional), negative reading flips the sign. Replace the first line by:

```python
reading = self.functional(state)
if np.iscomplexobj(reading):
    # ELEGANCE-DEBT[guard] #<issue> (W5): the scale functional is typed Callable[[Any], float] but a producer may return a complex carrier. Retires when the scale functional is typed real by construction (W5): the producer returns a checked real, and the machinery that exists today is the real-codomain functional of the three production constructors.
    raise ValueError(f"ScaleGauge: the functional {self.functional!r} read the complex value {reading!r}; the scale quotient is over the nonzero reals, so a complex reading has no section (complex reading)")
n = float(reading)
```
(the final text of the message is the builder's; the fragment `complex reading` and the naming of the functional are the contract). The refusal is tagged `ELEGANCE-DEBT[guard]`: it is a typed refusal on an alien carrier, and the machinery that would make it unnecessary (a functional typed real by construction) is in scope for W5; the issue is filed with the unit, before its commit (8.1), and the token line carries its number. The refusal precedes the `float` call, so no `ComplexWarning` is emitted for any complex dtype. `[M]` The producers: the three production constructors are `homogeneous/solver.py:562` (`ScaleGauge(production_rate.evaluate, 100.0)`), `sn/solver.py:2461` (`ScaleGauge(_StateProductionRate(solver), 1.0)`) and `sn/solver.py:2741` (`ScaleGauge(ke.compute_production_rate, 1.0)`); `_StateProductionRate.__call__` already returns `float(self.solver.compute_production_rate(phi))` (`sn/solver.py` about line 714), the same lossy `float` one hop earlier, so for that functional the new refusal cannot fire. **No shipped input is rejected by this item: its first red is a test-only input** (`np.complex128(1+2j)` from a test functional), and no in-tree producer of a complex reading is named. This widens nothing: every complex dtype is refused, one branch, and the only change for a carrier that was already refused is the exception type (a Python `complex` or a 0-d complex array, previously `TypeError` from inside `float`, now `ValueError`). Raw input: the functional's return value, any scalar-like. Gate (`tests/gates/numerics/test_gauge.py`): under `warnings.simplefilter("error")` and `python -O` (a subprocess leg with `cwd=<repository root>` for the `-O` claim): `np.complex128(1+2j)` raises `ValueError` with the fragment; `np.complex128(2+0j)` raises `ValueError` with the fragment (today it warns and returns 1.0; no `ComplexWarning` escapes); a real reading is unchanged (`np.float64(2.0)` with target 1.0 returns `0.5`); `1+2j` (Python complex) raises `ValueError` (was `TypeError`: the existing test, if any pins `TypeError` for it, is edited; `git grep -n "TypeError" -- tests/gates/numerics/test_gauge.py` printed nothing at base, so no existing pin changes). **First red** `[M]`: the `np.complex128(1+2j)` call returns `1.0` today (and `np.complex128(2+0j)` returns `1.0` after a `ComplexWarning`, which `warnings.simplefilter("error")` turns into a `ComplexWarning` exception and not the `ValueError` the gate requires). **Mutation**: monkeypatch `ScaleGauge.displacement` in-process with the old body `n = float(self.functional(state))`: red. Callers, `[M]` `git grep -n -P 'ScaleGauge|\.displacement\(' -- orpheus` and `git grep -l -w ScaleGauge -- tests`: production reaches `ScaleGauge.displacement` through `ScaleGauge.apply` at `orpheus/sn/solver.py` lines 2471 and 2744 (`gauge.apply(...)`) and `orpheus/homogeneous/solver.py` line 563 (`gauge.apply(phi_member)`); the `.displacement(` calls at `sn/solver.py:811` and `diffusion/solver.py:641` belong to `KernelProjectionGauge.displacement` (`kernel_projection_gauge.py`), a homonym, and `test_matrix_carrier_membership.py` calls `KernelProjectionGauge.direct_sum`, also the homonym. The 7 test files naming `ScaleGauge`: `tests/gates/numerics/test_gauge.py`, `tests/gates/numerics/test_outcome.py`, `tests/gates/homogeneous/test_eigen_verdict.py`, `tests/gates/homogeneous/test_homogeneous_outcome.py`, `tests/gates/homogeneous/test_kinf_exact_reference.py`, `tests/gates/sn/primitives/test_solution.py`, `tests/gates/sn/architecture/test_step3_solution_anchors.py`. Whether any caller passes a complex state is UNMEASURED `[R]`; none is expected to change, and a hit in the builder's run of those seven files and of `tests/gates/sn/solve` and `tests/gates/numerics` (serial, 6.6 item 5) is the first red.

### 6.4 u0.3: the `generator_as` residue

No gate is red today `[M]` (`m_pos3.py`, Python `re` word-bounded over 6 text suffixes of the tracked tree plus `.claude/skills` and `.claude/rules`): 6 lines in 5 files in the corpus `docs/`, `tests/`, `orpheus/`, `tools/`, `examples/`, `.claude/skills`, `.claude/rules` (suffixes `.py`, `.rst`, `.md`): `docs/theory/foundations/operator_algebra.rst` (1), `docs/theory/foundations/spaces.rst` (2), `docs/theory/methods/sn/history.rst` (1), `tests/gates/numerics/test_axis_generator.py` (1, a docstring), `tests/gates/numerics/test_axis_identity.py` (1, the live negative pin); `hasattr(Axis, "generator_as")` is False (control `hasattr(Axis, "of_measure")` True); AST uses in production: 0. A production-side gate would be green before and after (X1); a docs gate cannot tell tense; and `test_axis_identity.py` must keep naming the string. So this item's evidence is a reading by tense (E.19) plus `dead_references`, and no gate is written. The whole tracked tree holds 29 lines in 16 files (`[M]` `git grep -n -w generator_as`, counted by path: `.claude/plans` 17 lines in 7 files, `.claude/agent-memory` 2 in 2, tracked `scratch/` 4 in 2, `docs/` 4 in 3, `tests/` 2 in 2); only the `docs/` and `tests/` lines are in this unit's reading, and `.claude/plans`, `.claude/agent-memory` and `scratch/` are not in a unit's edit set. The reading, at the base hash, is by the rule the pages themselves declare (`history.rst`: "Every row below is a **dated** milestone and keeps the spelling that was current on its date"; `spaces.rst:8013` sits under the heading "Development history", about line 7675):
- `operator_algebra.rst` ("S7a: ``Axis.generator`` and its narrowing accessor ``generator_as`` became :attr:`Axis.origin`"): past, history, STAYS (not edited).
- `spaces.rst`, the paragraph beginning "The design record of the provenance slot before S7 (the ``generator`` field, its narrowing accessor ``generator_as``, ...": past, a pointer to an archived text, STAYS.
- `spaces.rst`, the line "The recovery goes via one accessor, ``Axis.generator_as`` — a typed narrow whose refusal names both the axis and the asking consumer, and which is load-bearing rather than defensive because ...": a dated history row in the historic present. It is NOT rewritten (rewriting would make a dated row assert a spelling that did not exist on its date). It gets the page's DATED RIDER, in the form `history.rst` already uses at about line 2045 ("⛔ Both halves of that last sentence are dated: ..."): "⛔ dated: S7a (`47e7bfaa`, `refactor(numerics): S7a, an axis is identified by what fixes its space`) retired ``generator_as`` into :attr:`Axis.origin`; the recovery above now reads ``case FromQuadrature(quadrature)`` on the origin." (`[M]` `git log -3 47e7bfaa` finds that commit.)
- `docs/theory/methods/sn/history.rst`, the line "``HarmonicFrame.for_space(angular_space, L)`` is the one blessed chain (space → angular axis → ``generator_as(Quadrature)`` → ``angular_frame(L)`` → ``from_galerkin``)": the same: NOT rewritten, gets the same dated rider with the live step (the builder reads `HarmonicFrame.for_space` first to spell it).
- `tests/gates/numerics/test_axis_generator.py`, the docstring "(it replaced G5_6b's ``generator_as`` call-site census)": past, STAYS. `test_axis_identity.py`: the pin stays.
Two riders are added and 0 lines are rewritten. Done-when predicate (restated with the corpus inlined, no scratch script): `git grep -n -w generator_as -- docs tests orpheus tools examples .claude/skills .claude/rules` prints 6 lines in 5 files at the base hash (the positive control; the 4 `docs/` lines and the 2 `tests/` lines above) and 8 lines in 5 files after the unit (the two riders each name the retired accessor on their own line, 2 more `docs/` lines), each of the 8 a past-tense line, a dated rider or the negative pin. `.venv/bin/nexus dead-references` lists no new row.

### 6.5 Retirement-audit items for the unit

u1.2/u1.3 retire nothing. E.19 applies to u0.3 (the tense reading above); B.4: `dead-references`; B.8: the sentence "load-bearing rather than defensive" justified the accessor and stays as a dated row; D.18: u1.2's refusal is the enforcement that keeps the union closed: its witness is the new gate and the static pyright law, in the same commit.

### 6.6 Done when

1. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics/test_convergence_evidence_reader_is_exhaustive.py tests/gates/sn/operators/test_gauge_evidence_reader_is_exhaustive.py tests/gates/numerics/test_gauge.py tests/gates/numerics/test_outcome.py tests/gates/homogeneous/test_eigen_verdict.py tests/gates/homogeneous/test_homogeneous_outcome.py tests/gates/homogeneous/test_kinf_exact_reference.py tests/gates/sn/primitives/test_solution.py tests/gates/sn/architecture/test_step3_solution_anchors.py` exits 0 (the seven `ScaleGauge` test files, `test_matrix_carrier_membership.py` is not among them: it names the homonym).
2. The gate files of item 1 contain the probes of u1.2 and u1.3: `NewKind()` and `None` raise `AssertionError` at both readers; `np.complex128(1+2j)` raises `ValueError`; `np.complex128(2+0j)` raises `ValueError` under `warnings.simplefilter("error")` (no `ComplexWarning` first); these are the legs of the gates, and the command of item 1 runs them.
3. `git grep -n -w generator_as -- docs tests orpheus tools examples .claude/skills .claude/rules` prints 8 lines in 5 files (6 at the base hash, the positive control), each a past-tense line, a dated rider or the negative pin (6.4).
4. `sphinx-build -E -W docs docs/_build/html` exits 0 and `.venv/bin/nexus dead-references` lists no new row for `generator_as`.
5. `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/sn/solve` and `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics` (directories, serial) exit 0: the run that 6.2 and 6.3 owe for whether a run builds an `ExitReport` carrying `None` or a functional reading a complex value.
6. `.venv/bin/pyright orpheus/numerics/convergence.py orpheus/sn/operators/loss_kernel_gauge.py` reports 0 errors, and `.venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/test_pyright_ratchet.py` exits 0.

## 7. The dependency table and the integration order

`[M]` (`measure_filesets.md`, five edit sets by AST token index over 1453 tracked `*.py` plus word-bounded text over docs): the five units' hard sets hold 232 distinct paths; 20 occur in more than one unit: 18 in two and 2 in three (`tests/gates/test_layer_imports.py`: dead-B, roles, insulation; `docs/theory/foundations/operator_algebra.rst`: dead, roles, and, in the measurement stage's reading, posing). 11 of the 20 are code (7 production, 4 test) and 9 are docs pages. Pairs: dead and roles 14, dead and posing 4, dead and insulation 2, dead and defects 2, roles and insulation 1, roles and posing 1 (the last removed by decision 4 of 1.2, u1.6 inside roles; and `operator_algebra.rst` leaves posing by 6.1, so the triple is now two). The units are NOT file-disjoint. The overlaps are resolved by one integration order (integrations are serial by B9), and builds may run in parallel in worktrees cut from `refactor/ontology`, each rebased onto the tip before its integration. Added by the amendment of 2026-11-20: a precondition unit `quadrature-ties` (2.0) precedes `insulation` and also edits `orpheus/numerics/generating_measure.py`, so that file is held by three units (the precondition, dead-A, `defects`); and `tests/gates/test_layer_imports.py` is edited in the order **roles, insulation, dead-B** (decision G10), after which the later `memo` and `registry` units of part B edit it (rows for `orpheus.memo`, the layer gate extended to the product order) and `tests/gates/numerics/test_concept_order.py` (the tier absorption); those later edits come after ALL of part A and are another specification's. Added by the third pass (2026-11-20): the sub-unit `angular-mass` of `defects` (5.7, the root fix of #708) has its own file set, measured by `scratch/ontology/w0/spec/amend/p708_overlap2.py` (exact-path matches of the file set of 5.7.1 against every other unit's sections and appendices, read line by line, never by basename): 6 code files shared with other units (`angular_lift.py` and `radial_characteristic.py` with `roles`; `fields/_bases.py` with dead-D and `defects`; `fields/angular_source_sink.py`, `test_typed_source_sinks.py` and `test_n2n_operator.py` with dead-D) and 8 docs pages shared (7.2); the new files (three gates and one rows module) are shared with none. The counts of this paragraph (232 paths, 20 in more than one unit) are the five units' and are not recomputed. Added by the fourth pass (2026-11-20): the precondition `quadrature-ties` (2.0) is the sixth unit and integrates first; its file set (2.0.2) is `orpheus/numerics/generating_measure.py`, `orpheus/numerics/recurrence.py`, `tests/gates/numerics/test_generating_measure.py`, `tests/gates/numerics/test_concept_order.py`, `docs/theory/foundations/discrete_measures.rst`, `docs/theory/verification/error_catalog.rst` and two new files (`orpheus/numerics/rational_polynomial.py`, `tests/gates/numerics/test_rational_polynomial.py`); computed by `scratch/ontology/w0/spec/amend/ties/overlap.py` (exact path matches by section), the files it shares are `generating_measure.py` (dead-A rows 54 and 55, `defects` #711), `test_generating_measure.py` (dead-A), `test_concept_order.py` (`insulation`), `discrete_measures.rst` (dead-A, and `defects` by prose) and `error_catalog.rst` (dead-A, `roles`, and `insulation` by read-for-tense); `recurrence.py` and the two new files are in no other unit. The counts of this paragraph's first sentences (232 paths, 20 in more than one unit) do not include the sixth unit and are not recomputed.

### 7.1 Integration order

| step | unit | why it sits here (what it must follow) |
|---|---|---|
| 0 | `quadrature-ties` (precondition, 2.0) | integrates BEFORE every other unit. Edits `orpheus/numerics/generating_measure.py`, `recurrence.py`, `tests/gates/numerics/test_generating_measure.py`, `test_concept_order.py`, `docs/theory/foundations/discrete_measures.rst`, `docs/theory/verification/error_catalog.rst` and adds `rational_polynomial.py` and its gate. `roles` shares none of these files; `insulation` (its tier entry in `test_concept_order.py`, and its G3 odd-n leg, which needs the fixed rule), dead-A (`generating_measure.py` rows 54 and 55, `test_generating_measure.py`, the two docs pages) and `defects` (`generating_measure.py`, #711) rebase onto it |
| 1 | `roles` | edits `numerics/operator.py`, `sn/coupled_system.py`, `transport/operators/transfer.py`, `isotropic_transfer.py`, `tests/gates/test_layer_imports.py` (its three string fixtures), `tests/gates/test_wave2_retirement.py` (the name walks move into the engine), `tests/gates/diffusion/test_boundary_realizer.py`, `tests/gates/sn/operators/test_capability_survival.py`, `test_psi_half_coupling.py`, and creates the ONE shared retirement engine `tests/gates/test_retirements.py` with `tests/gates/_retired/` (types, finder, `rows_wave2.py`, `rows_w0_roles.py`) that the dead sub-units extend |
| 2 | `insulation` | after the precondition unit (the swap to `LEGENDRE.on(a, b).gauss(n)` is total only once it has landed); edits `tests/gates/test_layer_imports.py` (the production-insulation function and legs, after `roles` and before dead-B), `tests/gates/numerics/test_concept_order.py` (the tier entry `"bickley"` at the head of the primitives tuple, one line away from the entry `"rational_polynomial"` that `quadrature-ties` added before `"recurrence"`: disjoint lines, `quadrature-ties` first), `docs/architecture/layering.rst`, `docs/theory/methods/collision_probability.rst`; the pages of its 2.1 prose list, among them `docs/theory/verification/error_catalog.rst`, are read for tense and edited only where a present-tense claim is false (that file is also edited by dead-A and `roles` at other lines); `docs/api/numerics.rst` is shared with dead |
| 3 | `dead-A` | `numerics/quadrature/registry.py`, `measure.py`, `generating_measure.py` (rows 54, 55), `test_registry.py`, the error catalogue gate; after the precondition unit (same file) and before `defects` (the #711 predicate covers the factories that remain) |
| 4 | `dead-B` | follows `roles` (`numerics/operator.py`, row 67 `InverseMetricOperator`) and `insulation` (`test_layer_imports.py`: its three string fixtures in the reference-insulation legs; the file order is roles, insulation, dead-B, decision G10) |
| 5 | `dead-C` | follows `roles` (`sn/coupled_system.py`, row 84; `test_psi_half_coupling.py`, rows 81 and 84); before `posing` (`sn/operators/loss_kernel_gauge.py`) |
| 6 | `dead-D` | follows `roles` (`transfer.py` seven rows, `isotropic_transfer.py` row 119, `test_capability_survival.py` row 121); before `defects` (`transport/fields/_bases.py`, row 105; row 106 is deferred) |
| 7 | `defects` | after dead-A and dead-D (the two files above) |
| 7b | `defects` sub-unit `angular-mass` (5.7) | after `roles` (the shared retirement engine `tests/gates/test_retirements.py` and the rows-module glob, which E6 of 5.7.3 uses; `angular_lift.py`, `radial_characteristic.py`), after dead-D (`fields/_bases.py` row 105, `fields/angular_source_sink.py` row 107, `test_typed_source_sinks.py` rows 107 and 112, `test_n2n_operator.py` row 122) and after `defects` proper (`fields/_bases.py`, class `MomentField`); before `posing` or independent of it (their docs hunks are line-disjoint) |
| 8 | `posing` | after dead-C; `loss_kernel_gauge.py` (row 91 moves `predicted_kernel_dimension` out, u1.2 edits `warn_if_gauge_freedom`: disjoint functions) and the docs pages `spaces.rst`, `history.rst` |

If the harness requires file-disjoint builds (not the case under B9 as ruled), move rows 67, 84, 119, 121, 122, 123, 124, 125, 138, 139 into `roles` (their files are roles files already), rows 54, 55, 105, 106 into `defects`, and row 91 into `posing`; the conflicts left are then the test files and eight docs pages.

### 7.2 The 14 code conflicts, each with the hunks and the resolution (and, after the table, the six of sub-unit `angular-mass`)

| file | unit A: its edit | unit B: its edit | resolution |
|---|---|---|---|
| `orpheus/numerics/operator.py` | roles: the marker classes, the metaclass, the joins, `_adopt_expansion_roles`, the comment of u1.6 | dead-B: row 67 `InverseMetricOperator` | disjoint definitions; roles first, dead-B rebases |
| `orpheus/sn/coupled_system.py` | roles: the `SystemRole` import and the write | dead-C: row 84 `build_coupled_system` | disjoint; roles first |
| `orpheus/transport/operators/transfer.py` | roles: the unused `BlockRole` import | dead-D: rows 121, 122, 123, 124, 125, 138, 139 | the import against seven method bodies; roles first |
| `orpheus/transport/operators/isotropic_transfer.py` | roles: import, two writes | dead-D: row 119 `dense_per_material` | disjoint; roles first |
| `orpheus/transport/fields/_bases.py` | defects: `MomentField.__post_init__` | dead-D: rows 105, 106 | disjoint classes; dead-D first |
| `orpheus/numerics/generating_measure.py` | defects: the intern table and `__post_init__` | dead-A: rows 54, 55 (`HERMITE`, `laguerre` leave) | not disjoint in intent (the #711 predicate must hold for the factories that remain): dead-A first; a correctness dependency, not a textual one |
| `orpheus/numerics/generating_measure.py` (third unit) | `quadrature-ties` (2.0): the exact-tie decision (E3: the call sites of the generic builders at the quoted `exact_recurrence(` lines, `_affine_remap`'s closure, `_rule_at`'s two rounding lines, `_certified_float`, a new `# The exact tie` block, the docstrings of `gauss`, `on` and the module, one message tail in `_correctly_rounded_rule`) | dead-A (rows 54, 55: the `HERMITE` constant and the `laguerre` and `_laguerre` factories leave) and `defects` (#711: the intern table and `__post_init__`, after the field `exact_recurrence`) | `quadrature-ties` first; the hunks are disjoint (the `HERMITE` block and the `laguerre` factories are at the end of the file; `defects` edits the class body, five lines from `quadrature-ties`'s edit of the docstring of the `exact_recurrence` attribute, a textual rebase); dead-A's moved `HERMITE` and `laguerre` fixtures pass `hermite_recurrence` and `laguerre_recurrence(a)` to `GeneratingMeasure(...)` unchanged (the builder names and their position are the same; only the second argument of the builder, which the fixtures never call, changed), so, `[R]` (the text of dead-A's fixtures, D-A.7, was not read for this), dead-A re-spells 0 lines and runs its rows against the new signature; `defects`'s content record reads `self.recurrence(8)`, the float view, which the unit keeps value for value; the order of the three is `quadrature-ties`, dead-A, `defects` |
| `tests/gates/numerics/test_generating_measure.py` | `quadrature-ties`: the promoted gate (G1 to G6), `catches` marker, one test edited (the stand-in lambda of `test_a_rule_that_rounds_differently_at_the_two_precisions_is_refused`) | dead-A: rows 54, 55 (the `HERMITE` and `laguerre` consumers in this file are rewired to the moved fixtures) | disjoint tests (the new tests are appended in the block "the correctly rounded construction refuses what it cannot certify"); `quadrature-ties` first |
| `tests/gates/numerics/test_concept_order.py` | `quadrature-ties`: the tier entry `"rational_polynomial"` | `insulation`: the tier entry `"bickley"` | different lines of one tuple; `quadrature-ties` first |
| `orpheus/sn/operators/loss_kernel_gauge.py` | posing: `warn_if_gauge_freedom` | dead-C: row 91 `predicted_kernel_dimension` | disjoint functions; dead-C first |
| `tests/gates/test_layer_imports.py` (three units) | roles: three string fixtures (a refused leg naming `BlockRole`, two admitted legs naming `stamp_boundary_role`) | insulation: the own-table gate and legs; dead-B: three string fixtures importing `dominant_eigenpair` in the reference-insulation refused legs | disjoint regions of one file; one integrates after the other in the order roles, insulation, dead-B (decision G10); the later `memo` and `registry` units of part B edit the file after all of part A |
| `tests/gates/diffusion/test_boundary_realizer.py` | roles: the `isinstance` test | dead-B: row 2 `minimal` rewire | disjoint tests; roles first |
| `tests/gates/sn/operators/test_capability_survival.py` | roles: five `isinstance` checks | dead-D: row 121 `on_basis` | disjoint; roles first |
| `tests/gates/sn/operators/test_psi_half_coupling.py` | roles | dead-C: rows 81, 84 | disjoint; roles first |

Sub-unit `angular-mass` (5.7), the six code files it shares, each read at the quoted section (`[M]` `p708_overlap2.py` plus the line read; every pair is disjoint regions of one file, so each resolution is the order, never a hand merge):

| file | other unit: its edit | `angular-mass`: its edit | resolution |
|---|---|---|---|
| `orpheus/transport/operators/angular_lift.py` | roles: the unused `BlockRole` import (line 89) | `total_weight`, `_combine`, the `axis_marginal` import line (line 114) | disjoint lines; roles first |
| `orpheus/sn/operators/radial_characteristic.py` | roles: the `SystemRole` import and the class attribute `system_role` of `RadialCharacteristicReconstruction` (about line 1087) and three more writes | `_ell0_scale` in `RadialCharacteristicReconstruction.__init__` (about line 1141) and `import math` | the same class, hunks 50 lines apart; roles first |
| `orpheus/transport/fields/_bases.py` | dead-D: row 105 `FaceField.face_views`; `defects`: `MomentField.__post_init__` | `BulkField.CONTAINS` and `BulkField.__add__` (lines 366 to 387) | three different classes; dead-D, `defects`, then `angular-mass`; its census of `CONTAINS` is re-run on the rebased tree |
| `orpheus/transport/fields/angular_source_sink.py` | dead-D: row 107 `AngularSourceSink.at_ordinate` | the arm, `__radd__`, `CONTAINS`, three imports, docstrings | disjoint methods; dead-D first |
| `tests/gates/sn/primitives/test_typed_source_sinks.py` | dead-D: row 107 deletes the smoke test `test_at_ordinate_selector`; row 112 makes a docstring true (line 442) | 5.7.6: 8 tests rewired or deleted | disjoint tests; dead-D first |
| `tests/gates/sn/operators/test_n2n_operator.py` | dead-D: row 122's consumer rewire | 5.7.6: deletes `test_total_weight_is_the_measure_mass` | disjoint tests; dead-D first |

Its 8 docs pages shared with other units, each edited at its own lines (5.7.8): `operator_algebra.rst` (roles, dead, `defects` by prose, and the five dead sub-units' tense lines), `spaces.rst` (`defects` by prose, posing, the dead sub-units), `field_algebra.rst` (`defects` by prose, dead-B, dead-C, dead-D), `indexing_and_layout.rst` (`defects` by prose, dead-A, dead-C, dead-D), `layering.rst` (insulation, dead-B), `conceptual_view.rst` (insulation, dead-A, dead-B), `adjoint.rst` (dead-A, dead-D), `cartesian_multid.rst` (`defects` by prose, dead-C, dead-D); `slab_multigroup.rst` is read for tense and is in no other unit's set.

The 9 docs conflicts are line-disjoint pages (`docs/architecture/layering.rst`: a dead mention against insulation's gate paragraph; `operator_algebra.rst`: roles lines at the quoted `BlockRole` mentions, dead via its rows; `docs/api/numerics.rst`: roles' heading, insulation's new entry, dead's rows; `error_catalog.rst`: dead-A's ERR-118 and ERR-081 and roles' ERR-075 lines; and others), each edited at its own lines. Hard-by-prose warnings (not conflicts): dead against insulation 8 docs pages, dead against defects 8 pages, roles against dead 3 (`coupled_block_operator.rst`, `operator_tensor_network.rst`, `transport/operators/angular_lift.py`), roles against insulation 1 (`error_catalog.rst`). Hard against suspect: roles against dead 1 (`tests/gates/diffusion/test_diffusion_realizer_table.py`, a `minimal` homonym candidate).

### 7.3 What the orchestrator does at integration, never a unit builder

- Regenerate the generated files once, after the last unit of a batch: `.venv/bin/python -m tools.verification.generate_error_index` (`.claude/skills/vv-principles/error_index.md`) and `.venv/bin/python -m tools.verification.generate_matrix` (`docs/theory/verification/matrix.rst`); a unit that deletes or adds a marked test changes what they print (`quadrature-ties`: the new entry ERR-127 (2.0.10), so arm 3 of `test_error_catalogue_reconciles.py` is red between the unit's commit and the regeneration; dead-A: ERR-081, ERR-118; roles: 41 deleted test functions, the double `_Tagged.__init__` and the deleted files `test_operator_protocols.py` and `test_operator_block_role.py`; the defects gates), so a unit builder editing them would conflict with every other unit.
- Read CI after each push (`gh run list --branch refactor/ontology --limit 3`); a unit is integrated by `git merge --ff-only` after a rebase; issues #709, #710, #711, #392 are closed by hand with the integrating hash; #641 is closed by hand when production's imports of `orpheus.derivations` are 0 and the edge gate is green (decision I5; its item 2, an independent verification chain, is the filed Bickley issue of 8.1), and its closing comment gives the hash and says what the unit did and did not do; #713 gets its correction comment.
- The full bar after the last unit: the suite (about 75 minutes, once), pyright, `sphinx-build -E -W`, `nexus dead-references`, the layer gate, the new gates.
- Sub-unit `angular-mass` (5.7): #708 and #288 are closed by hand with the integrating hash (#288, "the #207 cross-class source-sink dunder is statically unspellable", is settled by retiring the dunder; its closing comment says so); the file `.claude/plans/ontology_chart.md` records that **#708's root wave in the chart's B3 is superseded** (B3 gave it to W3; the root fix lands in W0: `all.json` named W2, "the axis section", and the section was the right object); the generated files are regenerated once at integration (this sub-unit deletes 2 test functions of `test_typed_source_sinks.py` and 1 of `test_n2n_operator.py`, and adds the gates of 5.7.5); the new issue of 8.1 (the half-trace mass) is filed by the orchestrator.
- Step 6 of the wave: the harness texts that name a retired object are edited through `docs/development/` sources and `python -m tools.harness`.

## 8. Decisions for the orchestrator

Each is a choice the inputs did not settle or a disagreement with an input; the structural reason is beside it; none is acted on beyond the specification. Numbering is stable; other sections cite it.

1. **The three re-baselined carry rows live in the test, not in `carry_head.npz`** (2.6.1). The fixture is a capture of a named commit (`meta/commit`, asserted by `test_the_fixture_holds_the_cp_rows`) read by 14 test files; an edit would make it lie about its commit. A 3-row override in `test_cp_boundary_carry.py` carries the reason and touches 0 binary files.
2. **The rounding-certificate defect behind `LEGENDRE.on(a, b).gauss(n)` is fixed at its root BEFORE the rule swap** (2.0; the orchestrator's decision I1, ruled again 2026-11-20 as 8.42 to 8.46). The precondition unit `quadrature-ties` integrates before `insulation`; the reasons, the measurements and the refused alternatives are in 8.44.
3. **The Ki3 copy rebuilds the interpolant from mpmath** (2.6.4), it does not store 64 coefficients: a stored table is an untraceable constant and a second source of truth in a unit whose purpose is a twin held equal by a witness. The witness compares the coefficient arrays bit for bit.
4. **ERR-118 "retired with its code" needs a state the catalogue gate does not have** (3.1, gate D3). Chosen: a `**Retired:** with its code: <names>` paragraph, a `_retired_ids()` that proves each named `file:qualname` pair is not defined in its file (the one finder of 1.5) and carries no marker, arm 2 reading it, and a new test with its own positive control; it changes a CI-run gate file (`test_error_catalogue_reconciles.py`) and joins dead-A's set. Alternative, not chosen: keep ERR-118 caught by the rewired witness `test_lebedev_order_13_is_not_positive_weighted` with a `catches("ERR-118")` marker and re-state the entry as a fixed defect whose subject is the built rule; that leaves no gate change but contradicts the binding decision's wording ("records ERR-118 as retired with its code").
5. **Rows 134 to 137 and the field `QuadratureSpec.degree_of_exactness_for`**: ruled 2026-11-20 at 8.23 (b): the four helpers and the field stay as catalogue data, kept for W8a, and the rows leave W0 (3.0 item 2, 3.3). The draft's first answer, to retire the field with its four values, is withdrawn: the field is an unruled closure member.
6. **None of the three factory cycles of Table 3 is opened in W0** (3.4): the realizer table waits for W7, the chord cycle for W6 (declared inside one concept), the measure partition for W8a. The inputs do not say whether Table 3 is in W0 and the measurement stage's question went unanswered; ruling any cycle in adds a full unit section, `derivations` files (5 for the chord) and files in three other units. Row 104 `RealizerTable.catalogue` is kept and its reason is amended to agree.
7. **`QuadratureSpec.structural_flags`**: ruled 2026-11-20 at 8.24 (a): it retires with `select_quadrature` and ERR-118, as the review's note and decision D4 read (3.0 item 3, 3.1 gate D3). The two fields only it read (`axis_aligned`, `half_range_clean`) stay as catalogue data, flagged for the registry unit (3.0 item 3).
8. **Move destinations and the layer rule** (3.1): a moved symbol lives in the gate directory named for the package the symbol left; ONE module per cluster when two or more test files consume it, inline when one does (decision D1); a moved method becomes a free function of its former receiver (decision D2); `tests/conftest.py` and shared existing helper modules are excluded. Reason: the consumer's (file, package) edge set is unchanged, so `OPEN_TEST_LAYER_VIOLATIONS` is unchanged. A "lowest consumer layer" rule, the previous draft's, was refuted on two members: `ENERGY` (a numerics object with a geometry consumer, `geometry/test_manifold.py`, already a listed open violation) and `assemble_ordinate_blocks` (imported by a derivations test that is already a listed open violation for `orpheus.sn`). `[M]` read at `tests/gates/test_test_layers.py` (the rows for `derivations/test_dsa_production_tie.py`: `orpheus.sn`, `orpheus.transport`; for `geometry/test_manifold.py`: `orpheus.numerics`). The scan was not run on the destinations by the measurement stage; each builder runs it (3.1).
9. **The wave tag of a refusal is enforced by a new leg in `tests/gates/test_elegance_debt_is_tagged.py`** (5.3): the file the measurement stage listed read-only joins defects' set. Without it the wave is prose (X3). The retiring wave deletes the map entry and the refusal in one diff.
10. **#710's predicate refuses `weights=None` against a non-counting measure** (5.2), contradicting the gates memo's "admitted twin that the #710 refusal must keep". Reason: `None` reads as ones, so an axis built that way with `FromMeasure(m)` and non-unit weights silently dresses the measure as the counting measure, which is the defect itself; and the census spied the predicate "`None` reads as ones" (1212 of 1212 constructions equal; 0 tests, 0 production callers reddened in 28 of 941 files). The memo's statement that the twin "constructs today" is true and is not an argument for keeping it.
11. **#709's predicate keeps the `TypeError` from `.head` for a non-tensor-product space and reddens one existing test**, which is rewired IN PLACE with `tests.gates.transport._phase_space_cases.moment_space` (5.2); the #709 gate lives at `tests/gates/transport/fields/`, the layer of the invariant, and builds no `SNProblem`.
12. **#711 is an intern table keyed by name with a content of `(support, first 8 recurrence coefficients as bytes)`** (5.2): the census's "name-to-content" predicate made precise; K = 8 is the module's own `_SYMMETRY_PROBE_ORDER`. Limit: measures agreeing on 8 coefficients are not told apart. A process-global table is a design fact.
13. **u1.2 needs two gate files, not one** (6.1): the second reader (`loss_kernel_gauge.py`, L3) cannot be imported by a numerics-directory test. u1.2's readers end in `assert_never` (a closed union read exhaustively is a law of the type, checked statically by pyright, not a guard), so they carry no `ELEGANCE-DEBT` tag; u1.3's refusal is tagged `ELEGANCE-DEBT[guard]`, retiring when the scale functional is typed real by construction (W5).
14. **Integration order** (7.1): the precondition `quadrature-ties`, roles, insulation, dead-A, dead-B, dead-C, dead-D, defects, posing (`tests/gates/test_layer_imports.py` in the order roles, insulation, dead-B; `generating_measure.py` in the order `quadrature-ties`, dead-A, defects; `test_concept_order.py` in the order `quadrature-ties`, insulation); builds in parallel, rebases serial. The measurement stage's alternative (moving 15 rows between units for file-disjoint builds) is recorded in 7.1 and not taken.
15. **The ghost row `pn`** (`tests/gates/test_layer_imports.py`, the strings `"pn"` in `L3_PACKAGES` and `FORBIDDEN_EDGES`; #411 gap 1; no `orpheus/pn/`) is in no dead row and belongs to the `registry` unit (the layer gate extended to the product order), not to `insulation`; the own-table gate's completeness leg reads the real directory names and is unaffected.
16. **One retirement engine, created by `roles`** (1.5; decisions G1 to G5, G8): `tests/gates/test_retirements.py` with rows modules `tests/gates/_retired/rows_w0_<unit>.py`, three row types, definition-site and import-resolution binding, one finder, a string leg; `tests/gates/test_wave2_retirement.py`'s name walks are ported into it. No `test_w0_retirements.py` and no `_w0_retired/`.
17. **Generated files and harness texts are regenerated once by the orchestrator at integration** (7.3).
18. **u0.3 shrinks to two edited lines** (6.4): of the 6 lines naming `generator_as`, 4 are past tense or the live negative pin; `posing` leaves `operator_algebra.rst`.
19. **For `dead`, the docs authority is the tense table, not the hard-docs list**: the measurement stage's hard docs counts (41 pages for dead, against 99 distinct pages in the census) differ in method; the tense tables are the line-level instrument and their union with the listed pages is the set a builder edits.
20. **#641 is closed by hand** when production's imports of `orpheus.derivations` are 0 and the edge gate is green (decision I5, superseding the draft, which kept it open); its item 2 (the independent verification chain of the kernel: production holds a twin of the reference's interpolant) becomes the filed Bickley issue. The closing comment says so. An accurate Bickley function moves the cylinder k by 1.2e-7 vacuum and 1.9e-7 white `[M]` and makes the cylinder band of `test_cp_first_collision.py` tightenable, a promotion under D.15.
21. **Posing counts**: the measurement stage listed 8 hard files for posing; this specification has 6 (the two `operator_algebra.rst` and `test_axis_generator.py` lines are history, 6.4).

### 8.0 Choices taken where the amendment decisions left a fork (one home per fact; each is the orchestrator's to overturn)

(Numbers 22 to 29 are the items of 8.4, written 8.22 to 8.29, all ruled on 2026-11-20.)

28. **The composite-rule helper's tag is spelled `ELEGANCE-DEBT[guard]`** (2.2 item 2): it is the only spelling of the token that the ledger gate reads (`[M]` 38 occurrences, all `[guard]`); the helper is machinery debt, not a guard, and the grammar has no other token. Decision I2's "guard to machinery" style has no gate that reads it.
29. **Row 126 (`AngularResidual.balance_map`) has one home, dead-D** (3.6): it lives in `transport/residuals/`, so the package cut puts it there; the draft listed it in dead-C's rows and in dead-D's bullet and absent from D-D.1. The alternative was to keep it in dead-C and say why a transport row sits in `sn`.
30. **Four moves are deferred by decision D2** (3.1 rule 2): rows 52, 99, 106, 116, each with the private attribute that no public accessor exposes (8.3). With rows 134 to 137 kept (8.23), W0 counts are 92 rows leaving production (49 retire, 43 move).
31. **`FromQuadrature` gains a `measure` property** (5.2, #710; check_elegance_defects finding 3): one `match` over the two arms replaces a repeated `isinstance`.
32. **The #711 record `_MeasureContent`** (5.2): a frozen record of the content, replacing the nested-pair annotation that disagreed with the flat content.
33. **The u0.3 riders are dated notes appended to two history rows** (6.4, decision P1): the rows are not rewritten.
34. **B5 of the chart ("the rule in the quadrature package") is superseded by binding decision 2** (the composite rule stays in `cp`, tagged): the orchestrator records it in `.claude/plans/ontology_chart.md`.

35. **The fifth spelling is `math.fsum` of the stored weights, not a section** (5.7.2 item 6, E5): `RadialCharacteristicReconstruction` holds a `Quadrature` and a split space, no axis-built angular space, and no mass verb exists on `Quadrature` or `DiscreteMeasure` (`[M]` `git grep` of 5.7.2 item 6), so the section's own definition is applied to the same array and G2 pins the equality. Not taken: minting `DiscreteMeasure.mass` as the one home that the section, `total_weight` and the arc scale would all read. That is the more principled shape (the mass of a discrete measure is a property of the measure, and the section's `math.fsum(covector * unit)` is its modal-axis generalisation); it is a new verb on a core type and is the orchestrator's to rule on; the sub-unit does not block on it.
36. **#708's root wave in the chart's B3 is superseded** (5.7, 7.3): B3 gave the #708 fix to W3 ("a scalar and an angular source live in different spaces"; `all.json` named W2, by the node `restrict`, the axis section); the orchestrator's ruling of 2026-11-20 puts the root fix in W0 and retires the mixed sum, so no wave carries a #708 refusal. The orchestrator records it in `.claude/plans/ontology_chart.md`, as for item 34.

37. **The generic form is an ARITHMETIC OBJECT passed to the builders** (2.0.4; a fork inside ruling 8.43). The second argument of a builder changes from the digit count to an `Arithmetic` (`HighPrecision(digits)` or `RATIONAL`); the builders are decorated by the public `over_arithmetic` (public because `generating_measure.py` also uses it for the remap). Not taken: a type parameter that is the number class (a class cannot spell `pi`, `gamma` or `beta`); a sibling function per arithmetic (a second closed form).
38. **An irrational mass is `None` in `beta[0]`, not an exception and not a second field** (2.0.4). The node decision never reads `beta[0]`; the weight decision reads `None` as "unavailable" and the value is refused as today. Not taken: raising `Irrational` out of the builder (the node decision needs the other coefficients of the same call, and a builder that raises gives it nothing); a parallel list of masses.
39. **The polynomial arithmetic over Q is its own module, `orpheus/numerics/rational_polynomial.py`** (2.0.5 E1), not private helpers of `generating_measure.py`: a polynomial over Q and its gcd are a concept that the file of the Gauss rule does not own, the file is already held by two other units, and the module is a leaf of the primitives tier. Cost: one tier entry in `tests/gates/numerics/test_concept_order.py` (a line `insulation` also edits; the order resolves it, 7.2). Not taken: putting the helpers in `generating_measure.py` (about 80 lines welded to the file), or `sympy` (`git grep -l -E "^\s*(import|from) sympy" -- orpheus/numerics` prints nothing, and its polynomial gcd on degree-63 rationals was not measured).
40. **The masses of Jacobi and Laguerre are rational at non-negative integer parameters** (2.0.4): a CHOICE, not a consequence. `RATIONAL.gamma`, `.beta` and `.power` return the factorial forms for positive integers and raise `Irrational` otherwise. Reason: `jacobi(0, 0) == LEGENDRE` by name (`[M]` `probes_ab_final.py`: True), so without it the two equal measures would round a width tie differently (the weight decision available for one and not for the other); measured with it, `jacobi(0.0, 0.0).on(0.2615032494292951, 0.7646335475579331).gauss(1)` equals the Legendre rule on the same panel. Cost: five lines in `_Rational`. Not taken: leaving every Jacobi and Laguerre mass `None` (two equal measures, two behaviours).
41. **The ERR entry's stub lands with the unit** (2.0.10; a reading of the gate, inside ruling 8.46). `test_every_marker_names_a_declared_entry` is red for a `catches` marker without a declared entry, so the stub (title, status, module, one caught-by line) is in the unit's commit; the next free id at the base hash is ERR-127 and arm 4 decides at integration; arm 3 (the generated index) waits for the orchestrator's regeneration (7.3).

### 8.1 Issues to file or amend as part of W0 (Cardinal Rule 4)

- New: an accurate production Bickley evaluator with an independent verification chain (binding decision 2 and decision I5; `module:cp` and `module:numerics`; relates #641 item 2, #597, #348).
- New, against W8a: `DiscreteMeasure.__add__` refuses to sum the panels of one interval (`direct sum requires equal spaces, got '[0,0.2]' and '[0.2,0.4]'`), so a composite rule is hand-concatenated in `cp` (decision I2; `module:numerics`, `level:L1`, `type:enhancement`). Filed BEFORE `insulation`'s commit: its number is on the helper's `ELEGANCE-DEBT[guard]` line.
- New: `CPParams.n_ki_table` and `CPParams.ki_max` have 0 production readers and `test_verification.py:1100` pins that they do nothing (decision I7: unread rows, not retired in W0; `module:cp`).
- New: `ConcentricPartition.chord` slots are decoded by position at the consumer (`slot_length[..., -n:]` in `CPMesh`); select the outbound slots through `slot_region` and the sense, or add an outbound accessor on `Chord` (`module:geometry`, a later wave).
- New: the f-string head of a dynamic import is a blind spot of `_dynamic_imports` and of the insulation gates that call it (`module:gates`).
- New, filed with unit `posing` before its commit: the `ScaleGauge` complex-reading refusal retires when the scale functional is typed real by construction (W5; `module:numerics`), the number on its `ELEGANCE-DEBT[guard]` line.
- Conditional: if `tests/gates/cp` takes materially longer than at base after `insulation` (`[R]`, unmeasured here), file the first-build cost of the panel rules as a follow-up; no cache is added.
- Amend: #713 (correct the premise with the 88-line, 6-of-14 measurement); #641 (closed by hand, decision I5, with the unit's hash and the filed Bickley issue); #411 (gap 1 to the registry unit).
- New, from the census of 5.7.2 item 5: `orpheus/sn/boundary/realized_operators.py` line 522, `IsotropicEmissionOperator(float(np.sum(in_weights)), ...)`, divides by the mass of the inflow HALF-trace measure with the same storage-order-dependent `np.sum` that #708 retires on the full sphere; it is another measure's mass, not re-pointed in W0 (`module:sn`, `level:L3`, `type:bug`; relates #708); and the possible `DiscreteMeasure.mass` verb of item 35 (`module:numerics`, `level:L1`, `type:enhancement`).
- New, from `quadrature-ties` (2.0.12; `module:numerics`, `level:L1`, `type:enhancement`): (a) an UNEQUAL single-boundary value is refused although a higher working precision would round it; escalate the precision up to a cap that still raises; (b) the weight decision's worst case (a straddle with a trivial gcd, 5.8 s at `n` = 63 `[M]` the memo, not re-run) can skip the gcd when the node is the exact rational centre, by computing that weight in `Fraction`. Neither blocks the unit (no shipped rule straddles one boundary at the shipped radius with an unequal value).
- Close by hand with the integrating hash: #709, #710, #711, #392, #641, and (sub-unit `angular-mass`) #708 and #288.

### 8.2 What this draft did not do

- No tracked file other than this one was written and no pytest run was made. The amendment ran AST and in-process probes, read-only on the tree, under `scratch/ontology/w0/spec/amend/` (`moves.py`, `priv.py`, `roles_count.py`, `verify_roles.py`, `p709.py`, `hollow_cp.py`, `gen_rows.py`); every count in the amended text is either one of those, a measurement of an earlier stage with its source, or marked `[R]`. The unit-level gates and mutations above are specified, not run.
- The 73 suspect method-name pairs, the 377 role tense lines, the full-suite red set, and the hollow-body chord are unmeasured (1.6).
- `retirement-audit` and `coding-elegance` skills: the retirement-audit skill was read as a file and its 25 items transcribed (1.4); the `coding-elegance` skill was not read in full by this stage. Every edit that creates code (the composite-rule helper, `ki3`, the three refusals, the two evidence readers) is for the reviewer to read against it; none introduces a new abstraction beyond a function and a cached constant.

### 8.3 Unread, for a later reviewed read (nothing in this list is touched in W0)

- **Rows deferred by decision D2** (the row stays in production; its body reads private state that no public accessor exposes): row 52 `Field.where_largest` (`Field._index_tuples`, also read by the kept `cone_violations`); row 99 `CollisionCache.reset_build_count` (the private class variable `CollisionCache._build_count`); row 106 `BoundaryField.from_face_arrays` (the private hook `FaceField._face_space_of`); row 116 `MaterialXSField.n2n_matrix` (`MaterialXSField._n2n_dense` and `MaterialXSField._build_dense_caches`). Probe limit: `self` and `cls` receivers only (3.1 rule 2).
- **Closure members that `dead/reviewed.json` does not rule** (decision G9; not retired and not moved; each is an orphan or a dependency of a moved body, and a later reviewed read decides): `_MAX_ELEMENT_ORDER` (row 34's closure), `Descent.max_abs_difference` (row 41), `_laguerre` (row 55; the moved `laguerre` imports it), `Mixture.in_scattering_xs` (row 1), `_anova_dimension` (row 91; the moved `predicted_kernel_dimension` imports it), `BoundaryField.from_face_arrays_on` (row 106), `TransferKernel.emission_matrix` (row 119). KEPT: `DiscreteMeasure.degree_of_exactness` and `AngularSymmetry.admits_symmetry` (row 73's closure), and the four helpers of rows 134 to 137 with the field `QuadratureSpec.degree_of_exactness_for` and `_lebedev_orders` (kept as catalogue data for W8a by 8.23; no longer closure members of a retiring row). `QuadratureSpec.structural_flags` is not in this list: it retires with row 73 (8.24).
- **Suspect pairs and tense tables (decision D6):** the 73 suspect (file, method-name) pairs in 66 test files (homonym or shared test file) and the tense tables D-A.6 to D-D.6 are not read by this specification; the sub-unit file-disjointness of 3.1 is conditional on them, and the serial order dead-A to dead-D makes a shared test file a textual rebase.
- **Unread census rows**: 113 classed by rule and 96 rejected, outside `reviewed.json` (binding decision 3).
- **Row 126's home** is a choice (decision 29): dead-D, by the package cut.

### 8.4 Ruled 2026-11-20 (formerly open; the orchestrator ruled each, the unit text states the ruled option, and the line gives the rejected option's structural reason)

8.22. **I4: the hollow-body gate cannot capture a k.** Ruled (a): the function-level hollow leg of G4 (chord arrays bit-equal for breakpoints `[0.1, 0.3, 0.5]`, `[M]` by the structure check) discharges decision I4; no `edges[0] != 0` refusal is added; a hollow `CPMesh` cannot be built today (#513), so no k exists to capture (2.2 item 1 and 5, 2.4 G4). Rejected (b), deferring I4 to a CP wave that lifts #513: until then no k can exist, so the deferral gates nothing, while the function-level leg already pins the same chord arrays today.
8.23. **Rows 134 to 137 against the field `QuadratureSpec.degree_of_exactness_for` (decision G9).** Ruled (b): rows 134 to 137 leave W0; the four helpers and the field stay as catalogue data, kept for W8a, and move to the KEPT table (3.0 item 2, 3.3, 3.6, D-A.1, D-A.7). W0 counts: 92 rows leave production (49 retire, 43 move), 45 rows kept, 49 `KEPT` entries with the 4 deferred. Rejected (a), retiring the four helpers: they are the values of a required field of a kept catalogue, so they are data and not dead, and retiring them forces literal parameter dicts, a constant derivable from a typed object, written in five places, and takes the degree-to-parameter map's only test.
8.24. **`QuadratureSpec.structural_flags` (decision G9 against decision D4).** Ruled (a): the flag retires with `select_quadrature` and ERR-118; the `**Retired:**` paragraph names both; the ERR-118 witness is rewired to the built rule, `assert not np.all(Quadrature.lebedev(order=13).weights > 0.0)`, with no `catches` marker (3.0 item 3, 3.1 gate D3, D-A.7). Rejected (b), keeping the flag: with the selector gone the method has no consumer and its subject entry is retired, so it would be a surviving accessor owned by nothing, and "retired with its code" would be true of the selector alone.
8.25. **u1.3, a complex reading with a zero imaginary part.** Ruled (b): every complex dtype is refused, nothing is widened (6.3, 6.6 item 2). Rejected (a), admitting a zero-imaginary reading as its real part: it adds a branch and widens a Python `complex` from `TypeError` to a reading for a consumer that does not exist in the tree (the three production functionals return real floats).
8.26. **The swap-mutation command in `boundary_conditions.rst` (about line 4624).** Ruled (a): the builder re-points the path and the page records the command as dated history, measured at the old path (4.4 E.22). Rejected (b), re-running the swap mutation after the re-home: a builder's worktree lacks the scratch plugin, and a re-run couples the page to a measurement taken outside the tree; a command measured at a path that existed then is correct dated history.
8.27. **Record the repeal of the ruling of 2026-10-20** ("a role is a fact about an atom", carried by `test_c42`). Ruled (a): the orchestrator writes the entry in the rulings ledger of `.claude/plans/ontology_chart.md` at close-out with the structural reason of section 4's Outcome; the unit text says so and edits no plan (4 Outcome). Rejected (b), the commit message alone: a commit message is not read at the start of a session, so the repeal would be invisible to the later session that meets the ruling in the ledger.

8.28. **#708's form: the root fix, with the mixed sum retired by deletion.** Ruled (the alternative): delete the `ScalarSourceSink` arm, `__radd__` and `AngularSourceSink.CONTAINS`, and `BulkField.CONTAINS` and `__add__` after a census with a positive control shows 0 remaining users (5.7.2 item 7); `Field._check_partner` refuses the sum under `python -O` (5.7.2 item 8); no tag, no `_WAVE_OF[708]`, `_bases.py` in the file set (7.2). Rejected, an explicit `TypeError` in the arm tagged `ELEGANCE-DEBT[guard] #708 (W3)`: it keeps a guard whose retiring wave is the one the root fix supersedes, adds a `_WAVE_OF` entry, and refuses what an existing gate already refuses.
8.29. **The endorsement of the cross-class dunder and the refined #207 principle are repealed** (5.7.8). Ruled: `operator_algebra.rst` lines 5793 to 5851 recorded the dunder as "endorsed at creation" and not to be tidied away, and `layering.rst` lines 163 to 166 gave it as a shape of cycle repair; the structural reason for the repeal: the dunder was the pullback `pi^*` with no division by the mass, it added a rate in `1/(cm^3 s)` to one in `1/(cm^3 s sr)`, its only production caller divided first and never needed it, and the section is the typed injection. The argument (a) of the five obstructions loses its first fact and stands on the second (the units identity read by `Field._check_partner`). The orchestrator writes the entry in the rulings ledger of `.claude/plans/ontology_chart.md` at close-out, as 8.27 does. Rejected, keeping the endorsement beside the retirement: the page would say in the present tense a thing the tree no longer does, which a later session reads as a live design.

8.42. **An exact value on a float64 rounding tie is decided in rational arithmetic and rounded half to even** (the orchestrator's ruling 1, 2026-11-20; 2.0.3). When the error interval holds exactly one float64 boundary, the code decides whether the exact value EQUALS it (a node: the orthogonal polynomial vanishes there; a weight: the gcd decision), and an equal value is `float(T)`, IEEE 754 `roundTiesToEven`, which is what "correctly rounded to nearest" means; anything else is refused exactly as today. The theory page `docs/theory/foundations/discrete_measures.rst` states the tie rule (2.0.9). Rejected, structural reasons: (b) the remap computing the midpoint more carefully: the midpoint is already exact at 40 and 60 digits, and no interval certificate can round an exact tie whatever its input; (c) a declared refusal (`SCOPE-BOUNDARY[guard]`): the refusal exists so that the code never guesses, an exact tie involves no guess, and it would refuse 384 of 1088 odd-`n` pairs of the measured panels; a catalogue of values known exact: not total (the `n = 3` weight `5s/9` is not dyadic in the reference and lands on a tie).
8.43. **The closed-form recurrences are generic over the number type, with no new stored field** (ruling 2; 2.0.4). One definition per family, evaluated in `HighPrecision(digits)` (today's values: 4368 of 4368 entries equal) or in `RATIONAL` (`fractions.Fraction`); the affine remap maps the rational recurrence exactly. Rejected, structural reasons: a second set of `Fraction` builders (a second closed form per family: the exact and the float value would agree tautologically and could drift apart, X4); a stored `rational_recurrence` field beside `exact_recurrence` (a derived value kept beside its definition, which is what the ruling forbids); deciding in `mpf` at a higher precision (an exact tie is undecidable at any finite precision). The form is the builder's reading of the ruling (8.0 items 37, 38 and 40).
8.44. **The unit `quadrature-ties` integrates BEFORE `insulation`, and first of all** (ruling 3; 2.0, 7.1, 7.2). The order of the units that hold `generating_measure.py` is `quadrature-ties`, dead-A, `defects`. A known defect is fixed before work is built on it: `[M]` 384 of 1088 odd-`n` (panel, n) pairs of 34 panels raise `ArithmeticError` and 0 of 1088 even-`n` pairs (parity is not the law: the width of a panel can be a tie), and (structure check) the carried cylinder and sphere SOLVE today at `n_quad_y` = 3 and 63 while under the per-panel spelling 6 of 12 and 4 of 12 panels raise. So there is no parity refusal, no strict expected failure (the specification never had one) and no deferred issue. Rejected, structural reasons: a parity guard (parity is not the law, and the guard over-refuses and under-refuses); mapping the one `[-1, 1]` rule to each panel (not correctly rounded per panel: 16 nodes and 248 weights differ from the reference at n = 64 and from today's numbers).
8.45. **The unit's gate is the investigation's 12 rows promoted, plus the witnesses it lacked** (ruling 4; 2.0.7). The 12 rows (10 tie rows, one untied control, one anti-guesser row) go into `tests/gates/numerics/test_generating_measure.py`; the existing test whose stand-in for `_certified_float` lacks the new keyword is edited as the investigation said (`**_`); the added witnesses are G4 (the unequal arm of the weight decision, and the decision's domain), G5 (the families other than Legendre), G6 (the generic builders), G7 (the polynomial module). Every mutation of the investigation is listed with its mechanism (2.0.7).
8.46. **An ERR entry is owed** (ruling 5; 2.0.10): a certificate correct for approximate values applied to an exact value on a tie. The catalogue file is `docs/theory/verification/error_catalog.rst`, the anchor form `.. error-entry:: ERR-NNN` (next free id ERR-127 at the base hash), the `catches` marker is on G1; the entry's text is the archivist's at the wave's documentation pass; the unit adds the stub because the catalogue gate's arm 1 requires the marker's target to exist (8.0 item 41).

### 8.5 Open for the orchestrator

8.47. **Does the odd-`n` text of `insulation` stay, and in what form?** The memo's NEEDS list and the brief of this pass say the strict expected failure at `n` = 3, law (f) of 2.3 and the decision about odd `n` (2.6.2) retire with the fix. The specification at `b68f8b8c` has no strict xfail (2.8 item 4: "no xfail is involved"; 2.6.2: "no strict xfail"), so there was none to retire. What this text did: law (f) is a one-sentence pointer to the law of `quadrature-ties` (2.0.3), 2.6.2 is a pointer to 8.44, and the wording that said the odd-`n` leg is red until the precondition lands is corrected. What it kept: G3's leg (`n_per_panel` = 3 and 63 on the carried edges, agreement with the reference rule to `rtol=1e-10`) as a plain test of the helper. Options: (a) keep the leg: it is the consumer-level witness that the swap does not turn a solving input into an error (6 of 12 cylinder panels and 4 of 12 sphere panels are tied), and it cites 2.0.3; (b) delete it, because G1 to G5 of 2.0.7 already pin the totality at the rule level and a leg that cannot redden today after the precondition lands is a test that never was red in the builder's tree (X1: it is red only on the base tree, which a builder never has); (c) also delete the pointer (f). Taken: (a); the orchestrator was asked by message (the premise mismatch) and had not answered when this text was written. **Ruled 2026-11-20 by the orchestrator: option (a).** G3's odd-`n` leg stays as a plain test of the helper; law (f) and 2.6.2 stay as pointers.

## Appendices (inlined from the reviewed inputs, because a builder's worktree has no `scratch/`)
Notation: `file:line kind in scope` in a production-reference column is `[M]` from the review's AST probe (`evidence.py`, 1453 tracked `*.py`; a `import in` is an import statement, `str in` a string constant equal to the name such as an `__all__` member, `name in`/`attr in` a use inside the named scope). Line numbers are as of `d2008d0e` and hold only until an earlier unit integrates; the symbol decides. The `closure` column lists production definitions referenced in the retired body whose every production reference lies in that body (a member nested inside the retired body goes with it; a member that is a separate definition retires or moves only if it is itself a row of `reviewed.json`, decision G9: 3.1 rule 3, section 8.3). Each `Docs` list is the measurement stage's hard docs set; the tense table beside it is the line-level authority and its union with that list is the set to edit. A tense table line is `PRESENT` (a present-tense claim made false: rewrite), `role` (a `:class:`/`:func:`/`:meth:`/`:attr:` role or backticked reference that is dead: drop or re-point to the moved home), `eq/ref` (a `:by:` target or label: re-point), `past` (history: stays), `collision` (the word is not the symbol: leave). `docstring` marks a line inside a docstring of a tracked `.py` file; `text` a line of a `.md`/`.rst` file.

**Corrections that override the tables below** (amendment of 2026-11-20): (1) rows 52, 99, 106 and 116 are NOT moved in W0 (decision D2, 3.1 rule 2; their table rows stay for reference); (2) row 126 `AngularResidual.balance_map` belongs to dead-D, although its table row sits in D-C.1; (3) a `home` or an `inline` in a destination column is superseded by the concrete home of the rows module of D-x.7; (4) "closure" members that are not rows are not touched (G9), except `QuadratureSpec.structural_flags`, ruled with row 73 (3.0 item 3); (5) rows 134 to 137 are KEPT (3.0 item 2): their table rows stay for reference with the verdict `keep (8.23)`, and the counts of this specification exclude them.

### Appendix D-A: dead-A: measure, quadrature selector, geometry, descent

#### D-A.1 Rows

| row | file : symbol | verdict | destination or rewire | test fate (per the review) | markers | production references (file:line kind in scope) | closure | test files naming it |
|---|---|---|---|---|---|---|---|---|
| 6 | `geometry/chart.py` `DirectionDomain.tangencies` | move | tests/gates/geometry helper (closed-form oracle) | 3 tests carry verifies('geometry-directions-at'): marker stays, helper import changes | test_the_tangencies_are_the_closed_form_break_set, test_no_tangency_on_the_slab_or_at_a_stratum, test_the_tangencies_are_where_the_kernels_crossings_change | - | - | derivations/test_characteristic_reading.py, geometry/test_chart_directions.py |
| 13 | `geometry/chord.py` `Chord.slot_start` | retire | rewire to `crossings.parameter[..., :-1]` | rewire geometry/test_chord.py::test_a_subnormal_projected_speed_lifts_to_infinite_lengths_never_nan (behavioural: keep the never-NaN law) | - | - | - | geometry/test_chord.py |
| 24 | `geometry/coord.py` `compute_volumes_1d` | retire | rewire to `coord.measure(edges)` | rewire 4 files; mesh/test_partition.py TestTheOneMeasure lists the name as a string (edit the list); geometry/test_coord.py::test_compute_volumes_1d_is_a_delegate deletes | - | geometry/__init__.py:50 import in ; geometry/__init__.py:84 str in  | - | geometry/test_coord.py, geometry/test_geometry.py, mesh/test_mesher.py, mesh/test_partition.py |
| 27 | `geometry/manifold.py` `UNIT_INTERVAL` | move | a shared test fixture constant (tests/gates conftest or _fixtures) | import moves | - | geometry/manifold.py:131 str in  | - | geometry/test_manifold.py, numerics/test_diagonal_operator.py, numerics/test_generating_measure.py, numerics/test_invariance_composition.py, numerics/test_measure.py, numerics/test_slab_orbit_space.py |
| 28 | `geometry/manifold.py` `barycentre` | retire | rewire to `Quotient.lift` | rewire geometry/test_manifold.py, numerics/test_registry.py imports | - | geometry/manifold.py:138 str in  | - | geometry/test_manifold.py, numerics/test_invariance_composition.py, numerics/test_registry.py, sn/sweep/test_angular_cell_partition.py |
| 30 | `geometry/structured_geometry.py` `StructuredGeometry.domain_extent_cm` | retire | rewire to `breakpoints[-1] - breakpoints[0]` | REWIRE, not delete: derivations/test_sood2003_to_geometry.py (4 tests) are behavioural (the Sood-to-geometry width); geometry/test_structured_geometry.py 1 test rewires | - | - | - | derivations/test_galerkin_spectral_basis_space.py, derivations/test_sood2003_to_geometry.py, geometry/test_structured_geometry.py |
| 31 | `geometry/symmetry.py` `SubgroupOfO3.identity_component` | retire | - | delete TestR1TheIdentityComponentIsTheConnectedHalf (API smoke over the roster) | - | - | - | geometry/test_symmetry.py |
| 32 | `geometry/symmetry.py` `_group_elements` | move | test helper `_realize(tag).elements` | import moves | - | - | - | geometry/test_manifold.py, geometry/test_symmetry.py |
| 33 | `geometry/symmetry.py` `OrbitCertificate.stabilizer_order` | move | test helper (orbit-stabilizer oracle) | - | - | - | - | numerics/test_invariance_verdicts.py, numerics/test_measure.py |
| 34 | `geometry/transformation.py` `RigidMotion.element_order` | move | test helper (group-law oracle) | - | - | - | geometry/transformation.py:197 _MAX_ELEMENT_ORDER | geometry/test_transformation.py |
| 36 | `geometry/transformation.py` `RigidMotion.permutes` | move | test helper over permutation_between (live) | test_transformation.py::test_o4 carries catches('ERR-073'): the marker stays on the test, which still reaches permutation_between's guard | test_o4_a_duplicated_point_is_refused_and_guard_two_is_the_one_that_fires | - | - | geometry/test_deck.py, geometry/test_symmetry.py, geometry/test_transformation.py, sn/operators/test_deck_kernel.py |
| 39 | `numerics/basis/descent.py` `Descent.discriminator` | retire | - | delete the one numerics/test_descent.py reference | - | - | - | data/test_energy_grid.py, data/test_mixture_condense.py, derivations/test_capability_matrices.py, numerics/test_descent.py, numerics/test_registry.py, sn/architecture/test_stage_separation.py, sn/architecture/test_step3_solution_anchors.py, sn/operators/test_fission_adjoint.py, sn/verification/analytical/test_angular_diffusion_limit_consistency.py, sn/verification/mms/test_mms_ld_2d.py, transport/_phase_space_cases.py, transport/fields/test_facefield_hierarchy.py, transport/test_phase_spaces.py |
| 40 | `numerics/basis/descent.py` `Descent.descending_real_slots` | retire | - | none | - | - | - | - |
| 41 | `numerics/basis/descent.py` `Descent.is_isomorphism` | retire | rewire the bit-agreement law into the test | rewire test_the_two_realizations_agree_at_the_BIT (behavioural) | - | - | numerics/basis/descent.py:221 Descent.max_abs_difference | numerics/test_descent.py |
| 51 | `numerics/elementary_operators.py` `DiagonalOperator.from_measure` | retire | - | delete TestDiagonalFromMeasure (2 API-smoke tests) | - | - | - | numerics/test_diagonal_operator.py, transport/frames/test_harmonic_frame.py |
| 54 | `numerics/generating_measure.py` `HERMITE` | move | test fixture | - | - | numerics/generating_measure.py:843 str in  | - | numerics/test_exactness.py, numerics/test_gauss_rules_correctly_rounded.py, numerics/test_generating_measure.py |
| 55 | `numerics/generating_measure.py` `laguerre` | move | test fixture (with HERMITE) | - | - | numerics/generating_measure.py:766 str in _laguerre; numerics/generating_measure.py:847 str in  | numerics/generating_measure.py:751 _laguerre | numerics/test_gauss_rules_correctly_rounded.py, numerics/test_generating_measure.py |
| 56 | `numerics/index_set.py` `ENERGY` | move | test fixture | - | - | numerics/index_set.py:34 str in  | - | geometry/test_manifold.py, numerics/test_basis_domain.py, numerics/test_descent.py, numerics/test_invariance_composition.py |
| 60 | `numerics/measure.py` `equispaced` | move | test fixture factory | - | - | numerics/measure.py:1914 str in  | - | numerics/test_dense_metric.py, numerics/test_invariance_composition.py, numerics/test_measure.py, numerics/test_rules_circle.py, numerics/test_rules_product.py, sn/sweep/test_angular_cell_partition.py, transport/frames/test_binding_tightness.py |
| 61 | `numerics/mesh_free_function.py` `Symbolic.from_srepr` | move | test helper (or the constructor) | - | - | - | - | numerics/test_content_identity_mesh_free.py, numerics/test_mesh_free_function.py |
| 72 | `numerics/quadrature/registry.py` `SelectionLog.summary` | retire | with select_quadrature (G1-c) | see row 73 | test_a_kernel_mode_carries_metric_but_no_NET_CURRENT | - | - | cp/test_cp_realizer_table.py, cross_method/test_admission_matrix_s4.py, cross_method/test_wave2_s4_tables.py, diffusion/_wave2_helpers.py, diffusion/test_diffusion_realizer_table.py, mc/test_mc_boundary_table.py, numerics/test_registry.py, sn/operators/test_loss_nullspace_reflective_box.py, test_docstring_xrefs.py, test_withdrawal.py, tools/test_cone_store_placement.py, tools/test_cone_store_placement_reads.py, transport/test_admission_matrix.py |
| 73 | `numerics/quadrature/registry.py` `select_quadrature` | retire | G1-c | ERR-081: re-point the stage-0 legs at AngularSymmetry.admits_symmetry (kept, G1-b). ERR-118: its subject QuadratureSpec.structural_flags retires with the selector (8.24); the witness test_lebedev_order_13_is_not_positive_weighted is rewired to the built rule, with no catches marker, and the ERR is recorded as retired with its code | TestStageZeroIsTheDescentArrowPlusTheUnspentSymmetry.test_the_rejection_message_names_the_one_failing_clause, test_a_rule_selected_for_positive_weights_has_positive_weights | numerics/quadrature/__init__.py:79 import in ; numerics/quadrature/__init__.py:119 str in ; numerics/quadrature/registry.py:1561 str in  | numerics/quadrature/registry.py:1500 select_quadrature._selection_log; numerics/quadrature/registry.py:1110 AngularSymmetry.admits_symmetry; numerics/quadrature/registry.py:1235 SelectionLog.angular_symmetry, numerics/quadrature/registry.py:1351 select_quadrature.angular_symmetry; numerics/quadrature/registry.py:549 QuadratureSpec.structural_flags; numerics/quadrature/registry.py:1360 select_quadrature.valid_flags | numerics/test_angular_symmetry.py, numerics/test_invariance.py, numerics/test_registry.py, numerics/test_slab_orbit_space.py |
| 132 | `numerics/quadrature/registry.py` `SelectionLog` | retire | G1-c | - | - | numerics/quadrature/__init__.py:79 import in ; numerics/quadrature/__init__.py:111 str in ; numerics/quadrature/registry.py:1273 name in QuadratureSelectionError.__init__; numerics/quadrature/registry.py:1511 name in select_quadrature._selection_log; numerics/quadrature/registry.py:1503 name in select_quadrature._selection_log; numerics/quadrature/registry.py:1284 name in select_quadrature; numerics/quadrature/registry.py:1559 str in  | - | numerics/test_registry.py |
| 133 | `numerics/quadrature/registry.py` `QuadratureSelectionError` | retire | G1-c | see row 73 (ERR-118 catcher catches this error type) | test_a_rule_selected_for_positive_weights_has_positive_weights | numerics/quadrature/__init__.py:79 import in ; numerics/quadrature/__init__.py:109 str in ; numerics/quadrature/registry.py:1523 name in select_quadrature; numerics/quadrature/registry.py:1557 str in  | - | numerics/test_registry.py |
| 134 | `numerics/quadrature/registry.py` `_gl1d_invert` | keep (8.23) | kept as catalogue data for W8a (3.0 item 2) | none: nothing is edited | - | numerics/quadrature/registry.py:821 name in  | numerics/quadrature/registry.py:600 _gl1d_invert.n_min | numerics/test_registry.py |
| 135 | `numerics/quadrature/registry.py` `_lebedev_invert` | keep (8.23) | kept as catalogue data for W8a (3.0 item 2) | none: nothing is edited | - | numerics/quadrature/registry.py:828 name in  | numerics/quadrature/registry.py:385 _lebedev_orders | numerics/test_registry.py |
| 136 | `numerics/quadrature/registry.py` `_ls_sn_invert` | keep (8.23) | kept as catalogue data for W8a (3.0 item 2) | none: nothing is edited | - | numerics/quadrature/registry.py:835 name in  | numerics/quadrature/registry.py:711 _ls_sn_invert._even_at_least; numerics/quadrature/registry.py:722 _ls_sn_invert.ceiling; numerics/measure.py:815 DiscreteMeasure.degree_of_exactness | numerics/test_registry.py |
| 137 | `numerics/quadrature/registry.py` `_product_invert` | keep (8.23) | kept as catalogue data for W8a (3.0 item 2) | none: nothing is edited | - | numerics/quadrature/registry.py:842 name in  | numerics/quadrature/registry.py:767 _product_invert.n_mu_min; numerics/quadrature/registry.py:768 _product_invert.n_phi_min | numerics/test_registry.py |

#### D-A.2 Hard production files

- `orpheus/geometry/__init__.py`
- `orpheus/geometry/chart.py`
- `orpheus/geometry/chord.py`
- `orpheus/geometry/coord.py`
- `orpheus/geometry/manifold.py`
- `orpheus/geometry/structured_geometry.py`
- `orpheus/geometry/symmetry.py`
- `orpheus/geometry/transformation.py`
- `orpheus/numerics/basis/descent.py`
- `orpheus/numerics/elementary_operators.py`
- `orpheus/numerics/generating_measure.py`
- `orpheus/numerics/index_set.py`
- `orpheus/numerics/measure.py`
- `orpheus/numerics/mesh_free_function.py`
- `orpheus/numerics/quadrature/__init__.py`
- `orpheus/numerics/quadrature/registry.py`

#### D-A.3 Hard test files

- `tests/gates/derivations/test_galerkin_spectral_basis_space.py`
- `tests/gates/derivations/test_sood2003_to_geometry.py`
- `tests/gates/geometry/test_chart_directions.py`
- `tests/gates/geometry/test_chord.py`
- `tests/gates/geometry/test_coord.py`
- `tests/gates/geometry/test_deck.py`
- `tests/gates/geometry/test_geometry.py`
- `tests/gates/geometry/test_manifold.py`
- `tests/gates/geometry/test_structured_geometry.py`
- `tests/gates/geometry/test_symmetry.py`
- `tests/gates/geometry/test_transformation.py`
- `tests/gates/mesh/test_mesher.py`
- `tests/gates/mesh/test_partition.py`
- `tests/gates/numerics/test_angular_symmetry.py`
- `tests/gates/numerics/test_basis_domain.py`
- `tests/gates/numerics/test_content_identity_mesh_free.py`
- `tests/gates/numerics/test_dense_metric.py`
- `tests/gates/numerics/test_descent.py`
- `tests/gates/numerics/test_diagonal_operator.py`
- `tests/gates/numerics/test_exactness.py`
- `tests/gates/numerics/test_gauss_rules_correctly_rounded.py`
- `tests/gates/numerics/test_generating_measure.py`
- `tests/gates/numerics/test_invariance.py`
- `tests/gates/numerics/test_invariance_composition.py`
- `tests/gates/numerics/test_measure.py`
- `tests/gates/numerics/test_mesh_free_function.py`
- `tests/gates/numerics/test_registry.py`
- `tests/gates/numerics/test_rules_circle.py`
- `tests/gates/numerics/test_rules_product.py`
- `tests/gates/numerics/test_slab_orbit_space.py`
- `tests/gates/sn/operators/test_deck_kernel.py`
- `tests/gates/sn/sweep/test_angular_cell_partition.py`
- `tests/gates/transport/frames/test_binding_tightness.py`

#### D-A.4 Hard docs files

- `docs/api/discrete_ordinates.rst`
- `docs/api/geometry.rst`
- `docs/api/method_of_characteristics.rst`
- `docs/api/numerics.rst`
- `docs/api/transport.rst`
- `docs/architecture/conceptual_view.rst`
- `docs/development/evidence/lessons.md`
- `docs/theory/conventions/indexing_and_layout.rst`
- `docs/theory/foundations/chart_and_chord.rst`
- `docs/theory/foundations/discrete_measures.rst`
- `docs/theory/foundations/manifolds.rst`
- `docs/theory/foundations/operator_algebra.rst`
- `docs/theory/foundations/spaces.rst`
- `docs/theory/foundations/spherical_harmonics.rst`
- `docs/theory/foundations/structured_geometry.rst`
- `docs/theory/methods/collision_probability.rst`
- `docs/theory/methods/sn/adjoint.rst`
- `docs/theory/methods/sn/angular_quadrature.rst`
- `docs/theory/methods/sn/curvilinear_one_group.rst`
- `docs/theory/methods/sn/index.rst`
- `docs/theory/references/characteristic.rst`
- `docs/theory/references/sood_registry.rst`
- `docs/theory/verification/error_catalog.rst`
- `docs/theory/verification/harness.rst`

#### D-A.5 Suspect test files (opened only when the retired owner class is named; homonyms left alone)

- `tests/gates/cp/test_cp_realizer_table.py`
- `tests/gates/cross_method/test_admission_matrix_s4.py`
- `tests/gates/cross_method/test_wave2_s4_tables.py`
- `tests/gates/data/test_energy_grid.py`
- `tests/gates/data/test_mixture_condense.py`
- `tests/gates/derivations/test_capability_matrices.py`
- `tests/gates/derivations/test_characteristic_reading.py`
- `tests/gates/diffusion/_wave2_helpers.py`
- `tests/gates/diffusion/test_diffusion_realizer_table.py`
- `tests/gates/mc/test_mc_boundary_table.py`
- `tests/gates/numerics/test_invariance_verdicts.py`
- `tests/gates/sn/architecture/test_stage_separation.py`
- `tests/gates/sn/architecture/test_step3_solution_anchors.py`
- `tests/gates/sn/operators/test_fission_adjoint.py`
- `tests/gates/sn/operators/test_loss_nullspace_reflective_box.py`
- `tests/gates/sn/verification/analytical/test_angular_diffusion_limit_consistency.py`
- `tests/gates/sn/verification/mms/test_mms_ld_2d.py`
- `tests/gates/test_docstring_xrefs.py`
- `tests/gates/test_withdrawal.py`
- `tests/gates/tools/test_cone_store_placement.py`
- `tests/gates/tools/test_cone_store_placement_reads.py`
- `tests/gates/transport/_phase_space_cases.py`
- `tests/gates/transport/fields/test_facefield_hierarchy.py`
- `tests/gates/transport/test_admission_matrix.py`
- `tests/gates/transport/test_phase_spaces.py`

#### D-A.6 Tense table (every prose line naming one of the sub-unit's retired names)

- `docs/api/discrete_ordinates.rst` (text): role 88
- `docs/api/geometry.rst` (text): PRESENT 111
- `docs/api/numerics.rst` (text): PRESENT 155,718
- `docs/development/evidence/lessons.md` (text): PRESENT 1075
- `docs/theory/foundations/chart_and_chord.rst` (text): PRESENT 3541; role 863,2683,2684
- `docs/theory/foundations/discrete_measures.rst` (text): PRESENT 779,799,827,3003,3022,3025,3031,3041,3059,3100,3102,3119; role 1897,2779,3058,3093,3095,3126,3127,3351; past 3454; collision 770,1013,3245
- `docs/theory/foundations/manifolds.rst` (text): PRESENT 4133,5279,5428,5910,5916,5922,5936,6132,7366,8354,8779,10641; role 2382,2383,2901,5411,5883,5903,8188,8268,9760,10338,11051; past 29,2450,10951,10974; collision 8580
- `docs/theory/foundations/operator_algebra.rst` (text): PRESENT 3164
- `docs/theory/foundations/spaces.rst` (text): PRESENT 3091
- `docs/theory/foundations/spherical_harmonics.rst` (text): PRESENT 971; role 770
- `docs/theory/foundations/structured_geometry.rst` (text): PRESENT 358,910,1654,3559,3583,3595,3645; role 152,300; past 3952
- `docs/theory/methods/collision_probability.rst` (text): role 116
- `docs/theory/methods/sn/angular_quadrature.rst` (text): role 263
- `docs/theory/references/characteristic.rst` (text): past 8658
- `docs/theory/references/sood_registry.rst` (text): PRESENT 444,554
- `docs/theory/verification/error_catalog.rst` (text): role 6224,6225,7081,7101; past 741,791,829,7092,7492,11623,11634
- `docs/theory/verification/harness.rst` (text): PRESENT 361
- `orpheus/derivations/continuous/fn_method/moment_space.py` (docstring): role 40
- `orpheus/derivations/continuous/galerkin_spectral/basis_space.py` (docstring): collision 682,723
- `orpheus/derivations/continuous/sood_registry/case.py` (docstring): role 197
- `orpheus/geometry/manifold.py` (docstring): PRESENT 323,1657; role 987,1142,1176,1470
- `orpheus/geometry/symmetry.py` (docstring): role 95,1754
- `orpheus/geometry/transformation.py` (docstring): role 69,217,296
- `orpheus/numerics/elementary_operators.py` (docstring): PRESENT 658
- `orpheus/numerics/frame_faces.py` (docstring): PRESENT 92
- `orpheus/numerics/mesh_free_function.py` (docstring): role 157
- `orpheus/numerics/quadrature/directional.py` (docstring): collision 168
- `orpheus/numerics/quadrature/registry.py` (docstring): PRESENT 404,539,1331,1338,1504; role 81,164,165,177,394,457,464,491,759,789,1199,1269,1340; collision 740
- `orpheus/numerics/quadrature/rules_circle.py` (docstring): PRESENT 37; collision 29
- `orpheus/numerics/quadrature/rules_product.py` (docstring): role 81; collision 2,568
- `orpheus/sn/angular/closure.py` (docstring): collision 1132
- `tests/gates/derivations/test_sood2003_to_geometry.py` (docstring): PRESENT 19,20,205
- `tests/gates/geometry/_measure.py` (docstring): PRESENT 15
- `tests/gates/geometry/test_chart_directions.py` (docstring): PRESENT 463
- `tests/gates/geometry/test_coord.py` (docstring): PRESENT 7,59
- `tests/gates/geometry/test_geometry.py` (docstring): role 4
- `tests/gates/geometry/test_manifold.py` (docstring): PRESENT 1771,1881,1882,2981,2984; role 864; past 3171
- `tests/gates/geometry/test_roots_of_unity.py` (docstring): collision 325
- `tests/gates/geometry/test_symmetry.py` (docstring): PRESENT 2089,2307
- `tests/gates/geometry/test_transformation.py` (docstring): PRESENT 1122,1416,1497,1502,1719; past 48
- `tests/gates/mesh/test_partition.py` (docstring): PRESENT 131,145
- `tests/gates/numerics/test_diagonal_operator.py` (docstring): role 25
- `tests/gates/numerics/test_gauss_rules_correctly_rounded.py` (docstring): collision 22
- `tests/gates/numerics/test_measure.py` (docstring): PRESENT 538
- `tests/gates/numerics/test_mesh_free_function.py` (docstring): PRESENT 261
- `tests/gates/numerics/test_registry.py` (docstring): PRESENT 235,486,1568,1570,1738; role 520,663,916,949,1064,1065,1180; past 236,441
- `tests/gates/numerics/test_rules_circle.py` (docstring): PRESENT 332
- `tests/gates/numerics/test_rules_product.py` (docstring): PRESENT 224; collision 3
- `tests/gates/numerics/test_symmetry_exactness.py` (docstring): role 80
- `tests/gates/sn/architecture/test_sn_boundary_binding.py` (docstring): collision 58
- `tests/gates/sn/operators/test_deck_kernel.py` (docstring): role 18,203,324
- `tests/gates/sn/sweep/core/test_sweep_regression.py` (docstring): collision 270
- `tests/gates/transport/_wave2_s4_helpers.py` (docstring): role 16
- `tests/gates/transport/test_admission_matrix.py` (docstring): role 12
### Appendix D-B: dead-B: numerics rest and eigenvalue

#### D-B.1 Rows

| row | file : symbol | verdict | destination or rewire | test fate (per the review) | markers | production references (file:line kind in scope) | closure | test files naming it |
|---|---|---|---|---|---|---|---|---|
| 0 | `data/emission_spectrum.py` `EmissionSpectrum.is_emitting` | retire | - | delete TestIsEmitting (API smoke, data/test_emission_spectrum.py) | - | - | - | data/test_emission_spectrum.py |
| 1 | `data/macro_xs/mixture.py` `Mixture.out_scattering_xs` | retire | - | none | - | - | data/macro_xs/mixture.py:228 Mixture.in_scattering_xs | - |
| 2 | `diffusion/method_space.py` `DiffusionMethodSpace.minimal` | retire | rewire to DiffusionMethodSpace() | rewire 9 test files to the constructor | - | - | - | diffusion/_wave2_helpers.py, diffusion/test_boundary_realizer.py, numerics/test_covector.py, numerics/test_full_field_space.py, numerics/test_registry.py, numerics/test_space_membership.py, sn/primitives/test_method_space.py, transport/test_material_field.py, transport/test_material_mesh.py |
| 42 | `numerics/basis/overlap_basis.py` `OverlapBasis.dominant_column` | move | test helper `overlap_table.argmax(axis=1)` | - | - | - | - | data/test_energy_grid.py |
| 43 | `numerics/basis/overlap_basis.py` `OverlapBasis.fractional_columns` | move | test helper | - | - | - | numerics/basis/overlap_basis.py:181 OverlapBasis.fractional_columns.frac | data/test_energy_grid.py, sn/test_condensation.py |
| 44 | `numerics/basis/spherical_harmonic_basis.py` `SphericalHarmonicBasis.evaluate_from_components` | retire | rewire to `evaluate(np.column_stack(...))` | REWIRE 4 tests (addition theorem, Lebedev orthogonality, on-axis, forged 1-D ordinates); the ERR-080 catcher keeps its marker and still reaches evaluate's refusal | TestEvaluateRefusesDirectionsOffTheSphere.test_the_forged_one_d_ordinates_are_refused | numerics/basis/spherical_harmonic_basis.py:510 attr in MirrorEvenSphericalHarmonicBasis.evaluate_from_components | - | numerics/test_spherical_harmonic_basis.py |
| 45 | `numerics/basis/spherical_harmonic_basis.py` `MirrorEvenSphericalHarmonicBasis.evaluate_from_components` | retire | as row 44 | as row 44 | TestEvaluateRefusesDirectionsOffTheSphere.test_the_forged_one_d_ordinates_are_refused | - | - | numerics/test_spherical_harmonic_basis.py |
| 46 | `numerics/coefficient_layout.py` `TensorLayout.multi_index` | move | test helper (np.unravel_index) | - | - | - | - | numerics/test_interval_legendre_basis.py |
| 52 | `numerics/field.py` `Field.where_largest` | move | test helper (diagnostic) | - | - | - | - | numerics/test_field.py, sn/solve/test_si_diagnostic_trajectory.py |
| 67 | `numerics/operator.py` `InverseMetricOperator` | retire | re-point its law tests at RieszRaiseOperator (numerics/operator.py:1438) | rewire test_inverse_metric_operator.py's laws (inverts off the null space, self-adjoint) onto RieszRaiseOperator if they are not already there; else delete | - | numerics/operator.py:155 str in  | - | numerics/test_inverse_metric_operator.py |
| 68 | `numerics/pencil.py` `OperatorPencil.is_regular` | retire | - | none | - | - | - | - |
| 69 | `numerics/pencil.py` `OperatorPencil.rhs_rank` | retire | - | delete TestLawTheRank | - | - | - | numerics/test_pencil.py |
| 70 | `numerics/pencil.py` `OperatorPencil.apply_at` | retire | - | none | - | - | - | - |
| 74 | `numerics/space.py` `FunctionSpace.project_live` | move | test helper | - | - | - | - | numerics/test_quotient_metric.py, sn/operators/test_quotient_adjoint_sn.py |
| 76 | `numerics/space.py` `angular_flux_space` | retire | - | delete 3 API-smoke tests in numerics/test_space.py | - | numerics/space.py:105 str in  | - | numerics/test_space.py |
| 77 | `numerics/space.py` `scalar_flux_space` | retire | - | delete 1 test | - | numerics/space.py:107 str in  | - | numerics/test_space.py |
| 49 | `numerics/eigenvalue.py` `direct_eigenvalue` | move | test-side dense oracle | also move the string spies: test_eigenvalue.py:561,589 (monkeypatch.setattr 'dominant_eigenpair'), test_homogeneous.py:269; layer-gate fixtures test_layer_imports.py:845,847,861 | - | - | numerics/eigenvalue.py:519 dominant_eigenpair | diffusion/test_solver.py, homogeneous/test_homogeneous.py, numerics/test_eigenvalue.py |
| 50 | `numerics/eigenvalue.py` `rayleigh_quotient_iteration` | retire | - | DELETE TestRayleighQuotientIteration (6 tests, test_eigenvalue.py:773-872) in the same unit, or they skip silently (vv-principles mode 8(6)) | - | - | numerics/eigenvalue.py:808 rayleigh_quotient_iteration.border; numerics/eigenvalue.py:811 rayleigh_quotient_iteration.step_norm | numerics/test_eigenvalue.py |
| 130 | `numerics/eigenvalue.py` `dominant_eigenpair` | move | with direct_eigenvalue | - | - | numerics/eigenvalue.py:698 name in direct_eigenvalue | numerics/eigenvalue.py:577 dominant_eigenpair.dominant | homogeneous/test_homogeneous.py, numerics/test_eigenvalue.py |
| 131 | `numerics/eigenvalue.py` `_sign_normalised` | move | with direct_eigenvalue | - | - | numerics/eigenvalue.py:587 name in dominant_eigenpair; numerics/eigenvalue.py:839 name in rayleigh_quotient_iteration | - | - |

#### D-B.2 Hard production files

- `orpheus/data/emission_spectrum.py`
- `orpheus/data/macro_xs/mixture.py`
- `orpheus/diffusion/method_space.py`
- `orpheus/numerics/basis/overlap_basis.py`
- `orpheus/numerics/basis/spherical_harmonic_basis.py`
- `orpheus/numerics/coefficient_layout.py`
- `orpheus/numerics/eigenvalue.py`
- `orpheus/numerics/field.py`
- `orpheus/numerics/operator.py`
- `orpheus/numerics/pencil.py`
- `orpheus/numerics/space.py`

#### D-B.3 Hard test files

- `tests/gates/data/test_emission_spectrum.py`
- `tests/gates/data/test_energy_grid.py`
- `tests/gates/diffusion/_wave2_helpers.py`
- `tests/gates/diffusion/test_boundary_realizer.py`
- `tests/gates/diffusion/test_solver.py`
- `tests/gates/homogeneous/test_homogeneous.py`
- `tests/gates/numerics/test_eigenvalue.py`
- `tests/gates/numerics/test_field.py`
- `tests/gates/numerics/test_interval_legendre_basis.py`
- `tests/gates/numerics/test_inverse_metric_operator.py`
- `tests/gates/numerics/test_pencil.py`
- `tests/gates/numerics/test_quotient_metric.py`
- `tests/gates/numerics/test_space.py`
- `tests/gates/numerics/test_spherical_harmonic_basis.py`
- `tests/gates/sn/operators/test_quotient_adjoint_sn.py`
- `tests/gates/sn/solve/test_si_diagnostic_trajectory.py`
- `tests/gates/test_layer_imports.py`

#### D-B.4 Hard docs files

- `docs/api/numerics.rst`
- `docs/architecture/conceptual_view.rst`
- `docs/architecture/layering.rst`
- `docs/development/evidence/plan-authoring.md`
- `docs/theory/foundations/boundary_conditions.rst`
- `docs/theory/foundations/field_algebra.rst`
- `docs/theory/foundations/frame.rst`
- `docs/theory/foundations/infinite_medium.rst`
- `docs/theory/foundations/operator_algebra.rst`
- `docs/theory/foundations/operator_inverse_family.rst`
- `docs/theory/foundations/spaces.rst`
- `docs/theory/foundations/spherical_harmonics.rst`
- `docs/theory/methods/diffusion_1d.rst`
- `docs/theory/methods/sn/history.rst`
- `docs/theory/verification/reference_solutions.rst`

#### D-B.5 Suspect test files (opened only when the retired owner class is named; homonyms left alone)

- `tests/gates/sn/test_condensation.py`

#### D-B.6 Tense table (every prose line naming one of the sub-unit's retired names)

- `.claude/skills/vv-principles/scripts/issue226_spectral_invisibility.md` (text): PRESENT 14
- `docs/api/numerics.rst` (text): role 60,65,69,444,446
- `docs/architecture/conceptual_view.rst` (text): role 596
- `docs/architecture/layering.rst` (text): role 726
- `docs/development/evidence/plan-authoring.md` (text): past 550
- `docs/theory/foundations/boundary_conditions.rst` (text): PRESENT 2837,7193; past 7189,7529,8705
- `docs/theory/foundations/cross_section_data.rst` (text): PRESENT 2024
- `docs/theory/foundations/field_algebra.rst` (text): role 1296; past 1661,1707,1812
- `docs/theory/foundations/frame.rst` (text): role 3422,3423,3433,3434,3668,3669,3704; past 5286,5305
- `docs/theory/foundations/infinite_medium.rst` (text): PRESENT 2027,2110,2298,2333,2796,2802; role 72,2029,2220,2320,2323,2742,2749,2753,2790,2809,3517; eq/ref 1252,2032,2688; past 2329
- `docs/theory/foundations/operator_algebra.rst` (text): PRESENT 7503,7986,8357; role 7460,7467; past 11385
- `docs/theory/foundations/operator_inverse_family.rst` (text): PRESENT 1685,1717; role 1701,1702
- `docs/theory/foundations/spaces.rst` (text): PRESENT 3726; role 1406,1407,4367; eq/ref 3743
- `docs/theory/foundations/spherical_harmonics.rst` (text): role 1118
- `docs/theory/methods/diffusion_1d.rst` (text): PRESENT 1185; role 1110; past 2108
- `docs/theory/methods/sn/history.rst` (text): role 2606,3271,3274; past 640,3277
- `docs/theory/methods/sn/index.rst` (text): PRESENT 41
- `docs/theory/verification/reference_solutions.rst` (text): PRESENT 857; role 427,676
- `orpheus/data/emission_spectrum.py` (docstring): PRESENT 15,61
- `orpheus/data/energy_grid.py` (docstring): role 58
- `orpheus/diffusion/method_space.py` (docstring): role 72
- `orpheus/diffusion/solver.py` (docstring): role 114
- `orpheus/numerics/basis/spherical_harmonic_basis.py` (docstring): role 533
- `orpheus/numerics/eigenvalue.py` (docstring): PRESENT 603,666,787; role 513,514,536,659,677,714
- `orpheus/numerics/field.py` (docstring): PRESENT 388
- `orpheus/numerics/frame_faces.py` (docstring): past 218
- `orpheus/numerics/gauge.py` (docstring): PRESENT 95
- `orpheus/numerics/matrix_inverse_operator.py` (docstring): role 55
- `orpheus/numerics/operator.py` (docstring): PRESENT 3209,3231; role 2529
- `tests/gates/data/test_emission_spectrum.py` (docstring): PRESENT 7,179
- `tests/gates/data/test_energy_grid.py` (docstring): PRESENT 15,17,218,221,244,502,539,619
- `tests/gates/diffusion/test_boundary_realizer.py` (docstring): past 18
- `tests/gates/diffusion/test_solver.py` (docstring): role 8
- `tests/gates/homogeneous/test_homogeneous.py` (docstring): PRESENT 242,302,307,308
- `tests/gates/homogeneous/test_homogeneous_problem.py` (docstring): PRESENT 185
- `tests/gates/numerics/test_eigenvalue.py` (docstring): PRESENT 1,20,31,55,543,572,573,646,787,857; role 4,18,435,438,469; past 54,544
- `tests/gates/numerics/test_field.py` (docstring): PRESENT 348
- `tests/gates/numerics/test_inverse_metric_operator.py` (docstring): PRESENT 1
- `tests/gates/numerics/test_quotient_metric.py` (docstring): PRESENT 275; role 5
- `tests/gates/numerics/test_spherical_harmonic_basis.py` (docstring): role 3
- `tests/gates/sn/_test_helpers.py` (docstring): role 889
- `tests/gates/sn/primitives/test_method_space.py` (docstring): PRESENT 6
- `tests/gates/sn/solve/test_si_convergence_diagnostics.py` (docstring): role 141
- `tests/gates/sn/solve/test_si_diagnostic_trajectory.py` (docstring): PRESENT 362; role 27,291; past 11,86
- `tests/gates/sn/test_condensation.py` (docstring): PRESENT 354
### Appendix D-C: dead-C: sn

#### D-C.1 Rows

| row | file : symbol | verdict | destination or rewire | test fate (per the review) | markers | production references (file:line kind in scope) | closure | test files naming it |
|---|---|---|---|---|---|---|---|---|
| 37 | `mesh/axis.py` `n_unknowns_flat` | move | test helper (independent count of the flat size) | - | - | mesh/__init__.py:34 import in ; mesh/__init__.py:92 str in ; sn/problem.py:31 import in  | mesh/axis.py:515 n_unknowns_flat.n_face | sn/primitives/test_axis_primitive.py |
| 80 | `sn/angular/closure.py` `_MMHalfGrid.downstream_per_ordinate` | move | with compute_psi_half_per_level | - | - | - | - | sn/sweep/curvilinear/test_compute_psi_half_per_level.py |
| 81 | `sn/angular/closure.py` `_MMHalfGrid.upstream` | move | with compute_psi_half_per_level | - | TestBitIdenticalSlab.test_slab_first_cell_bit_identical, TestBitIdenticalSlab.test_slab_interior_cell_bit_identical, TestBitIdenticalSlab.test_slab_negative_ordinate_bit_identical, TestCylindricalDegenerate.test_degenerate_cell_synthetic, TestResidual.test_slab_residual_closed_form, TestResidual.test_curvilinear_residual_matches_cell_balance, TestResidual.test_cylindrical_degenerate_residual_matches_cell_balance, test_dd_per_cell_recurrence_matches_symbolic_derivation | - | - | sn/operators/test_psi_half_coupling.py, sn/sweep/core/test_cell_balance_for_streaming.py, sn/sweep/core/test_diamond.py, sn/sweep/core/test_discretization_scheme_protocol.py, sn/sweep/core/test_sweep_graph.py, sn/sweep/curvilinear/test_compute_psi_half_per_level.py, sn/sweep/slab/test_dd_recurrence.py, transport/spatial/test_ld_ubld_symbolic.py, transport/spatial/test_linear_discontinuous.py |
| 82 | `sn/angular/closure.py` `compute_psi_half_per_level` | move | test helper module | - | - | sn/angular/closure.py:2347 str in  | - | sn/_test_helpers.py, sn/sweep/curvilinear/test_compute_psi_half_per_level.py |
| 84 | `sn/coupled_system.py` `build_coupled_system` | retire | rewire to build_within_group_system | rewire test_psi_half_coupling.py TestCoupledBuilder (behavioural) | - | sn/__init__.py:1 import in ; sn/coupled_system.py:220 str in  | - | sn/operators/test_psi_half_coupling.py |
| 85 | `sn/loss_representation/assembly.py` `ordinate_walk_order` | move | with assemble_ordinate_blocks (assembly.py to the test side) | - | test_g2_walk_order_triangularity_is_exact, test_g2_lapack_forward_substitution_equals_sweep | sn/loss_representation/assembly.py:127 str in  | - | sn/sweep/test_assembly_mode.py |
| 86 | `sn/loss_representation/assembly.py` `assemble_ordinate_blocks` | move | assembly.py to the test side (fuller-view oracle) | 4 verifies markers stay; derivations/test_dsa_production_tie.py imports it | test_g1_assembled_matvec_equals_apply, test_g2_walk_order_triangularity_is_exact, test_g2_lapack_forward_substitution_equals_sweep, test_g1_ld_assembled_matvec_equals_apply | sn/loss_representation/assembly.py:126 str in  | sn/loss_representation/assembly.py:186 _probe_coefficient_blocks; sn/loss_representation/assembly.py:368 assemble_ordinate_blocks._to_global_frame; sn/loss_representation/assembly.py:341 assemble_ordinate_blocks.cell_index; sn/loss_representation/assembly.py:364 assemble_ordinate_blocks.dof_signs | derivations/test_dsa_production_tie.py, sn/sweep/core/test_multi_d_reverse_walk.py, sn/sweep/test_assembly_mode.py |
| 87 | `sn/loss_representation/sweep_graph.py` `SweepDependencyGraph.face_in_x` | move | test helper (face_in[0]) | - | - | - | - | sn/primitives/test_dag_ownership.py, sn/sweep/core/test_sweep_graph.py, sn/sweep/core/test_sweep_graph_nd_admission.py |
| 88 | `sn/loss_representation/sweep_graph.py` `SweepDependencyGraph.face_out_x` | move | as 87 | - | - | - | - | sn/primitives/test_dag_ownership.py, sn/sweep/core/test_sweep_graph.py, sn/sweep/core/test_sweep_graph_nd_admission.py |
| 89 | `sn/loss_representation/sweep_graph.py` `SweepDependencyGraph.face_in_y` | move | as 87 | - | - | - | - | sn/primitives/test_dag_ownership.py, sn/sweep/core/test_sweep_graph.py, sn/sweep/core/test_sweep_graph_nd_admission.py |
| 90 | `sn/loss_representation/sweep_graph.py` `SweepDependencyGraph.face_out_y` | move | as 87 | - | - | - | - | sn/primitives/test_dag_ownership.py, sn/sweep/core/test_sweep_graph.py, sn/sweep/core/test_sweep_graph_nd_admission.py |
| 91 | `sn/operators/loss_kernel_gauge.py` `predicted_kernel_dimension` | move | test-side counting-law oracle | 5 verifies markers stay | test_the_dimension_matches_the_combinatorial_counting_law, test_the_singularity_needs_TWO_closed_axis_pairs, test_the_T_plus_R_split_is_exact_per_quadrature, test_S16a_the_DD_left_basis_spans_the_dense_left_null_space, test_S8_the_union_count_and_spans | sn/operators/loss_kernel_gauge.py:237 str in  | sn/operators/loss_kernel_gauge.py:526 _anova_dimension | sn/operators/test_loss_kernel_gauge.py, sn/operators/test_loss_nullspace_reflective_box.py, sn/operators/test_predicted_loss_kernel.py |
| 92 | `sn/problem.py` `SNProblem.n_unknowns_flat` | move | as row 37 | - | - | mesh/__init__.py:34 import in ; mesh/__init__.py:92 str in ; sn/problem.py:31 import in  | - | sn/primitives/test_axis_primitive.py |
| 93 | `sn/problem.py` `SNProblem.loss_kernel_gauge` | retire | rewire to `problem.pencil.kernel.projector` | rewire 6 test files | test_the_T_plus_R_split_is_exact_per_quadrature, test_every_MIRROR_EVEN_functional_is_blind_to_the_gauge | - | - | sn/architecture/test_step2_terminal_object_anchors.py, sn/mesh/test_reflective_axis_pairs.py, sn/operators/test_loss_kernel_gauge.py, sn/operators/test_loss_nullspace_reflective_box.py, sn/solve/test_every_entry_gauges_its_trace.py, sn/solve/test_reflective_si_iteration_budget.py |
| 94 | `sn/problem.py` `SNProblem.boundary_face_layout` | retire | rewire to `problem.angular_trace.layout` | rewire 3 test files | - | - | - | sn/operators/test_snmesh_realizer_wiring.py, sn/primitives/test_axis_native_construction.py, sn/primitives/test_boundary_face_layout.py |
| 95 | `sn/solver.py` `boundary_vs_interior_split` | move | test helper | - | - | - | - | sn/operators/test_typed_residual_evaluation.py |
| 96 | `sn/solver.py` `SNSolver.compute_group_absorption_rate` | retire | - | cross_method/test_one_parse_on_production.py lists it by string in _site_calls: edit the site table; rewire sn/operators/test_solver_components.py | - | - | - | cross_method/test_one_parse_on_production.py, sn/operators/test_solver_components.py |
| 97 | `sn/splitting.py` `LossTerm.is_gain` | retire | - | none | - | - | - | - |
| 98 | `sn/splitting.py` `Splitting.law_residual` | move | test-side law instrument | - | - | - | sn/splitting.py:409 Splitting.law_residual.loss_image | sn/architecture/test_splitting_value.py, sn/architecture/test_step2_terminal_object_anchors.py |
| 99 | `sn/sweep/cache.py` `CollisionCache.reset_build_count` | move | test side (monkeypatch the counter) | rewire the one test | - | - | - | sn/sweep/core/test_cache.py |
| 126 | `transport/residuals/angular_residual.py` `AngularResidual.balance_map` | retire | - | delete 1 test | - | - | - | sn/operators/test_typed_residual_evaluation.py |
| 141 | `sn/loss_representation/assembly.py` `_ordinate_label` | move | with assembly.py | - | - | sn/loss_representation/assembly.py:174 name in ordinate_walk_order; sn/loss_representation/assembly.py:354 name in assemble_ordinate_blocks | - | - |
| 142 | `sn/loss_representation/assembly.py` `_iter_cells_in_walk_order` | move | with assembly.py | - | - | sn/loss_representation/assembly.py:180 name in ordinate_walk_order; sn/loss_representation/assembly.py:403 name in assemble_ordinate_blocks | - | - |
| 143 | `sn/loss_representation/assembly.py` `_probe_coefficient_blocks` | move | with assembly.py | - | - | sn/loss_representation/assembly.py:349 name in assemble_ordinate_blocks | sn/loss_representation/assembly.py:238 _probe_coefficient_blocks.as_cell_columns; sn/loss_representation/assembly.py:243 _probe_coefficient_blocks.as_face_columns; sn/loss_representation/assembly.py:217 _probe_coefficient_blocks.cell_field; sn/loss_representation/assembly.py:214 _probe_coefficient_blocks.cell_tail; sn/loss_representation/assembly.py:226 _probe_coefficient_blocks.face_fields; sn/loss_representation/assembly.py:215 _probe_coefficient_blocks.face_tail; sn/loss_representation/assembly.py:213 _probe_coefficient_blocks.zeros_q | - |

#### D-C.2 Hard production files

- `orpheus/mesh/__init__.py`
- `orpheus/mesh/axis.py`
- `orpheus/sn/__init__.py`
- `orpheus/sn/angular/closure.py`
- `orpheus/sn/coupled_system.py`
- `orpheus/sn/loss_representation/assembly.py`
- `orpheus/sn/loss_representation/sweep_graph.py`
- `orpheus/sn/operators/loss_kernel_gauge.py`
- `orpheus/sn/problem.py`
- `orpheus/sn/solver.py`
- `orpheus/sn/splitting.py`
- `orpheus/sn/sweep/cache.py`
- `orpheus/transport/residuals/angular_residual.py`

#### D-C.3 Hard test files

- `tests/gates/cross_method/test_one_parse_on_production.py`
- `tests/gates/derivations/test_dsa_production_tie.py`
- `tests/gates/sn/_test_helpers.py`
- `tests/gates/sn/architecture/test_splitting_value.py`
- `tests/gates/sn/architecture/test_step2_terminal_object_anchors.py`
- `tests/gates/sn/mesh/test_reflective_axis_pairs.py`
- `tests/gates/sn/operators/test_loss_kernel_gauge.py`
- `tests/gates/sn/operators/test_loss_nullspace_reflective_box.py`
- `tests/gates/sn/operators/test_predicted_loss_kernel.py`
- `tests/gates/sn/operators/test_psi_half_coupling.py`
- `tests/gates/sn/operators/test_snmesh_realizer_wiring.py`
- `tests/gates/sn/operators/test_solver_components.py`
- `tests/gates/sn/operators/test_typed_residual_evaluation.py`
- `tests/gates/sn/primitives/test_axis_native_construction.py`
- `tests/gates/sn/primitives/test_axis_primitive.py`
- `tests/gates/sn/primitives/test_boundary_face_layout.py`
- `tests/gates/sn/primitives/test_dag_ownership.py`
- `tests/gates/sn/sweep/core/test_cache.py`
- `tests/gates/sn/sweep/core/test_multi_d_reverse_walk.py`
- `tests/gates/sn/sweep/core/test_sweep_graph.py`
- `tests/gates/sn/sweep/core/test_sweep_graph_nd_admission.py`
- `tests/gates/sn/sweep/curvilinear/test_compute_psi_half_per_level.py`
- `tests/gates/sn/sweep/test_assembly_mode.py`

#### D-C.4 Hard docs files

- `docs/api/discrete_ordinates.rst`
- `docs/api/mesh.rst`
- `docs/theory/conventions/indexing_and_layout.rst`
- `docs/theory/foundations/boundary_conditions.rst`
- `docs/theory/foundations/field_algebra.rst`
- `docs/theory/foundations/operator_algebra.rst`
- `docs/theory/foundations/operator_inverse_family.rst`
- `docs/theory/foundations/spaces.rst`
- `docs/theory/methods/sn/cartesian_multid.rst`
- `docs/theory/methods/sn/curvilinear_numerics.rst`
- `docs/theory/methods/sn/curvilinear_one_group.rst`
- `docs/theory/methods/sn/history.rst`
- `docs/theory/methods/sn/index.rst`
- `docs/theory/methods/sn/loss_representation.rst`
- `docs/theory/methods/sn/solver.rst`
- `docs/theory/verification/error_catalog.rst`
- `docs/theory/verification/sn.rst`

#### D-C.5 Suspect test files (opened only when the retired owner class is named; homonyms left alone)

- `tests/gates/sn/solve/test_every_entry_gauges_its_trace.py`
- `tests/gates/sn/solve/test_reflective_si_iteration_budget.py`
- `tests/gates/sn/sweep/core/test_cell_balance_for_streaming.py`
- `tests/gates/sn/sweep/core/test_diamond.py`
- `tests/gates/sn/sweep/core/test_discretization_scheme_protocol.py`
- `tests/gates/sn/sweep/slab/test_dd_recurrence.py`
- `tests/gates/transport/spatial/test_ld_ubld_symbolic.py`
- `tests/gates/transport/spatial/test_linear_discontinuous.py`

#### D-C.6 Tense table (every prose line naming one of the sub-unit's retired names)

- `docs/api/discrete_ordinates.rst` (text): PRESENT 199,243; role 209,212,214,225
- `docs/api/mesh.rst` (text): role 367
- `docs/development/evidence/coding-elegance.md` (text): PRESENT 128,530
- `docs/theory/conventions/indexing_and_layout.rst` (text): role 1637
- `docs/theory/foundations/boundary_conditions.rst` (text): PRESENT 9508,9843; role 6386,7468,9389,9439,9582,9863; past 9606,9607
- `docs/theory/foundations/field_algebra.rst` (text): PRESENT 1682; role 1400,1409
- `docs/theory/foundations/manifolds.rst` (text): PRESENT 1172,6861; role 1179
- `docs/theory/foundations/operator_algebra.rst` (text): PRESENT 8458,8577,8740,10459; role 6606,8510,8512,8514,8838
- `docs/theory/foundations/operator_inverse_family.rst` (text): PRESENT 2467; role 1989,2166
- `docs/theory/foundations/spaces.rst` (text): past 1135
- `docs/theory/methods/sn/cartesian_multid.rst` (text): PRESENT 513,514,515,6261,6370; role 6205,6355,6664,6666,6763,6773,6778,6902,6904,6906; eq/ref 6117,6126,6129,6246,6255
- `docs/theory/methods/sn/curvilinear_numerics.rst` (text): role 251,1223
- `docs/theory/methods/sn/curvilinear_one_group.rst` (text): PRESENT 3712,3719,4171,4179,4188; role 3775,4148; eq/ref 3732
- `docs/theory/methods/sn/history.rst` (text): role 95,929,3074,3157; past 673
- `docs/theory/methods/sn/index.rst` (text): PRESENT 46; role 1017
- `docs/theory/methods/sn/loss_representation.rst` (text): PRESENT 214; role 64,922,1119,4322; eq/ref 974,1061,1107
- `docs/theory/methods/sn/solver.rst` (text): PRESENT 1240; role 254,888,1010,1164,2183,2184,3965,3981; past 2491
- `docs/theory/verification/cone_store.rst` (text): PRESENT 812
- `docs/theory/verification/error_catalog.rst` (text): role 8432,8481
- `docs/theory/verification/sn.rst` (text): PRESENT 2716,2744,2793; role 1295,2242,2282,2729,3182,3327,3466
- `orpheus/mesh/axis.py` (docstring): role 23
- `orpheus/numerics/basis/orthonormal_span_basis.py` (docstring): PRESENT 33
- `orpheus/numerics/face_layout.py` (docstring): role 44
- `orpheus/numerics/gauge.py` (docstring): role 42
- `orpheus/numerics/index_set.py` (docstring): PRESENT 43
- `orpheus/numerics/kernel_projection_gauge.py` (docstring): role 39
- `orpheus/numerics/spaces/angular_trace_space.py` (docstring): role 438
- `orpheus/sn/angular/closure.py` (docstring): role 722,738,1457,1627
- `orpheus/sn/coupled_system.py` (docstring): PRESENT 1
- `orpheus/sn/loss_representation/__init__.py` (docstring): PRESENT 797; role 739
- `orpheus/sn/loss_representation/assembly.py` (docstring): role 152
- `orpheus/sn/loss_representation/sweep_graph.py` (docstring): PRESENT 404,408
- `orpheus/sn/operators/__init__.py` (docstring): PRESENT 27; role 28
- `orpheus/sn/operators/loss_kernel_gauge.py` (docstring): role 178
- `orpheus/sn/problem.py` (docstring): PRESENT 653; role 174,854,1019,1067,1080,1490
- `orpheus/sn/solver.py` (docstring): PRESENT 361; role 361,1678
- `orpheus/sn/splitting.py` (docstring): role 22
- `orpheus/transport/fields/angular_boundary_flux.py` (docstring): role 37
- `orpheus/transport/spatial/scheme.py` (docstring): role 1202
- `tests/gates/cross_method/test_one_traced_face_list.py` (docstring): PRESENT 4
- `tests/gates/derivations/test_dsa_production_tie.py` (docstring): PRESENT 9
- `tests/gates/numerics/test_angular_face_trace_space.py` (docstring): PRESENT 77
- `tests/gates/sn/_test_helpers.py` (docstring): PRESENT 595; role 576,1190
- `tests/gates/sn/architecture/test_stage_separation.py` (docstring): role 43
- `tests/gates/sn/architecture/test_step2_terminal_object_anchors.py` (docstring): PRESENT 479; role 350
- `tests/gates/sn/operators/test_loss_kernel_gauge.py` (docstring): role 7,22,262
- `tests/gates/sn/operators/test_loss_nullspace_reflective_box.py` (docstring): PRESENT 123; role 4,8,9,124,125,142,196,261,362,436
- `tests/gates/sn/operators/test_predicted_loss_kernel.py` (docstring): role 6,9
- `tests/gates/sn/operators/test_psi_half_coupling.py` (docstring): PRESENT 3115
- `tests/gates/sn/operators/test_snmesh_realizer_wiring.py` (docstring): PRESENT 8,296
- `tests/gates/sn/operators/test_typed_residual_evaluation.py` (docstring): PRESENT 8,132
- `tests/gates/sn/primitives/test_axis_primitive.py` (docstring): PRESENT 440; role 12
- `tests/gates/sn/primitives/test_boundary_face_layout.py` (docstring): role 1
- `tests/gates/sn/primitives/test_face_name_crosswalk.py` (docstring): PRESENT 66; past 7
- `tests/gates/sn/solve/test_boundary_gs_is_a_coherent_splitting.py` (docstring): role 273
- `tests/gates/sn/solve/test_every_entry_gauges_its_trace.py` (docstring): role 45
- `tests/gates/sn/sweep/core/test_multi_d_reverse_walk.py` (docstring): PRESENT 946
- `tests/gates/sn/sweep/core/test_phase_c_gates.py` (docstring): PRESENT 155
- `tests/gates/sn/sweep/core/test_sweep_graph_nd_admission.py` (docstring): PRESENT 507
- `tests/gates/sn/sweep/curvilinear/test_angular_closure.py` (docstring): role 31,260,287,329,404
- `tests/gates/sn/sweep/curvilinear/test_compute_psi_half_per_level.py` (docstring): PRESENT 1,30,152,178,378,392,445; role 11,139
### Appendix D-D: dead-D: transport

#### D-D.1 Rows

| row | file : symbol | verdict | destination or rewire | test fate (per the review) | markers | production references (file:line kind in scope) | closure | test files naming it |
|---|---|---|---|---|---|---|---|---|
| 79 | `numerics/spaces/scalar_trace_space.py` `ScalarTraceSpace.slot_shape` | retire | rewire to `layout.faces[face].shape` | rewire 1 behavioural test (test_slab_dim_faces_and_area_metric) | - | - | - | transport/fields/test_scalar_boundary_flux.py |
| 101 | `transport/boundary_conditions/admission.py` `Approximated.issues` | retire | - | none | - | - | - | tests/_harness/audit.py, cross_method/test_admission_matrix_s4.py, sn/sweep/curvilinear/test_cyl_sweep_regression.py, sn/sweep/curvilinear/test_sph_sweep_regression.py, test_pending_ports.py, test_withdrawal.py, tools/test_harness_generator.py, transport/_wave2_s4_helpers.py, transport/test_admission_matrix.py, transport/test_admission_mechanism.py |
| 102 | `transport/boundary_conditions/admission.py` `FaceSite.hollow_inner` | retire | - | delete test_hollow_inner_is_derived_never_stored | - | - | - | transport/test_admission_mechanism.py |
| 103 | `transport/boundary_conditions/binding.py` `Binding.faces_of` | retire | rewire into the test | rewire test_binding_invariants (behavioural) | - | - | - | transport/test_binding.py |
| 105 | `transport/fields/_bases.py` `FaceField.face_views` | move | test helper | - | - | - | - | transport/fields/test_angular_boundary_flux.py |
| 106 | `transport/fields/_bases.py` `BoundaryField.from_face_arrays` | move | test helper | - | - | - | transport/fields/_bases.py:1260 BoundaryField.from_face_arrays_on | sn/mesh/test_the_mints_are_the_raw_constructors.py |
| 107 | `transport/fields/angular_source_sink.py` `AngularSourceSink.at_ordinate` | retire | - | delete 1 smoke test | - | - | - | sn/primitives/test_typed_source_sinks.py |
| 108 | `transport/fields/harmonic_moment_flux.py` `HarmonicMomentFlux.l_block` | retire | census G5 | delete own-SUT tests | - | - | - | sn/primitives/test_harmonic_moment_flux.py, transport/fields/test_harmonic_moment_flux.py |
| 109 | `transport/fields/harmonic_moment_flux.py` `HarmonicMomentFlux.isotropic_part` | retire | census G5 | delete own-SUT tests | - | - | - | sn/primitives/test_harmonic_moment_flux.py, transport/fields/test_harmonic_moment_flux.py |
| 110 | `transport/fields/harmonic_moment_flux.py` `HarmonicMomentFlux.anisotropic_part` | retire | census G5 | delete own-SUT tests | - | - | - | sn/primitives/test_harmonic_moment_flux.py, transport/fields/test_harmonic_moment_flux.py |
| 111 | `transport/fields/harmonic_moment_flux.py` `HarmonicMomentFlux.truncate` | retire | census G5 | delete own-SUT tests in 6 files | - | - | transport/fields/harmonic_moment_flux.py:320 HarmonicMomentFlux.truncate.new_head; transport/fields/harmonic_moment_flux.py:328 HarmonicMomentFlux.truncate.new_shape; transport/fields/harmonic_moment_flux.py:329 HarmonicMomentFlux.truncate.new_values | numerics/test_moment_head_axis_built_premise.py, sn/mesh/test_harmonic_moment_flux_on_the_problem.py, sn/mesh/test_moment_space_is_read_off_the_frame_on_the_problem.py, sn/mesh/test_spatial_moment_tail_on_the_problem.py, sn/primitives/test_harmonic_moment_flux.py, transport/fields/test_harmonic_moment_flux.py |
| 112 | `transport/fields/harmonic_moment_source_sink.py` `HarmonicMomentSourceSink._CELL_GROUP_NAME` | retire | - | test_typed_source_sinks.py:442 docstring becomes true | - | - | - | - |
| 113 | `transport/fields/scalar_boundary_flux.py` `ScalarBoundaryFlux.p1_boundary_scalar_flux` | move | test helper `2(J+ + J-)` | - | - | - | - | diffusion/test_solver.py, transport/fields/test_scalar_boundary_flux.py |
| 114 | `transport/fields/scalar_flux.py` `ScalarFlux.at_group` | retire | - | none | - | - | - | - |
| 115 | `transport/kernels.py` `FissionKernel.dyad` | retire | - | rewire :1055 to np.outer; delete :479-480 | - | - | - | cross_method/test_one_parse_on_production.py, homogeneous/test_coda_anchors.py, numerics/_energy_pencil.py, numerics/test_outer_dyad.py, sn/operators/test_fission_adjoint.py, sn/operators/test_fission_kernel_crosscheck.py, sn/operators/test_isotropic_fission.py, transport/test_kernels.py |
| 116 | `transport/mesh/material_xs_field.py` `MaterialXSField.n2n_matrix` | move | test helper (dense oracle) | - | TestProducerSideNormalisation.test_typed_apply_returns_per_ordinate_already_normalised | - | - | homogeneous/test_operator_spaces.py, sn/operators/test_scattering_operator.py, sn/operators/test_solver_components.py, transport/test_kernels.py |
| 117 | `transport/mesh/material_xs_field.py` `MaterialXSField.fission_production_per_material` | move | test helper | - | - | - | - | transport/test_kernels.py |
| 118 | `transport/mesh/material_xs_field.py` `MaterialXSField.chi_per_material` | move | test helper | - | - | - | - | transport/test_kernels.py |
| 119 | `transport/operators/isotropic_transfer.py` `IsotropicTransfer.dense_per_material` | move | test helper | - | - | - | transport/kernels.py:243 TransferKernel.emission_matrix | sn/operators/test_isotropic_scattering.py, transport/test_energy_leaves_as_matrix.py, transport/test_kernels.py |
| 121 | `transport/operators/transfer.py` `LegendreMomentTransfer.on_basis` | move | test helper (production mints through TransferOperator._moment_transfer, transfer.py:575) | - | - | - | - | sn/_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py, sn/mesh/test_moment_space_is_read_off_the_frame_on_the_problem.py, sn/operators/test_capability_survival.py, sn/operators/test_legendre_moment_scattering.py, sn/operators/test_scattering_adjoint.py, sn/primitives/test_harmonic_moment_flux.py, transport/test_tier2_equivalence_s_family.py |
| 122 | `transport/operators/transfer.py` `LegendreMomentTransfer.on_frame` | move | as 121 | - | - | - | - | sn/operators/test_frame_conjugate_carve.py, sn/operators/test_n2n_operator.py |
| 123 | `transport/operators/transfer.py` `TransferOperator.foldable_part` | retire | census G7 | delete own-SUT tests | - | - | - | sn/operators/test_scattering_operator.py |
| 124 | `transport/operators/transfer.py` `TransferOperator.residual_part` | retire | census G7 | as 123 | - | - | - | sn/operators/test_scattering_operator.py |
| 125 | `transport/operators/transfer.py` `TransferOperator.is_foldable_into_sigma_r` | retire | census G7 | as 123 | - | - | - | sn/operators/test_scattering_operator.py |
| 129 | `transport/spatial/_ubld.py` `D1ClosedForm.kernel_rhs` | move | test helper | - | - | - | - | transport/spatial/test_ld_ubld_primitive.py |
| 138 | `transport/operators/transfer.py` `TransferOperator.foldable_sigma` | retire | - | re-point tests/gates/sn/architecture/_config.py:318; edit the allow-list of test_dsa_rate.py:1157 | - | sn/acceleration/dsa.py:258 attr in DSALowOrderSystem.from_problem | - | derivations/test_dsa_production_tie.py, sn/acceleration/test_dsa_low_order.py, sn/acceleration/test_dsa_rate.py, sn/architecture/_config.py, sn/operators/test_scattering_operator.py |
| 139 | `transport/operators/transfer.py` `TransferOperator._sibling` | retire | - | - | - | transport/operators/transfer.py:776 attr in TransferOperator.foldable_part; transport/operators/transfer.py:804 attr in TransferOperator.residual_part | - | - |
| 140 | `numerics/spaces/moment_head.py` `truncated_head` | retire | - | delete numerics/test_heads_on_the_axis.py H7 and test_moment_head_axis_built_premise.py rows that test it | - | numerics/spaces/__init__.py:55 import in ; numerics/spaces/__init__.py:71 str in ; numerics/spaces/moment_head.py:43 str in ; transport/fields/harmonic_moment_flux.py:316 import in HarmonicMomentFlux.truncate; transport/fields/harmonic_moment_flux.py:320 name in HarmonicMomentFlux.truncate | - | numerics/test_heads_on_the_axis.py, numerics/test_moment_head_axis_built_premise.py, sn/primitives/test_harmonic_moment_flux.py |

#### D-D.2 Hard production files

- `orpheus/numerics/spaces/__init__.py`
- `orpheus/numerics/spaces/moment_head.py`
- `orpheus/numerics/spaces/scalar_trace_space.py`
- `orpheus/sn/acceleration/dsa.py`
- `orpheus/transport/boundary_conditions/admission.py`
- `orpheus/transport/boundary_conditions/binding.py`
- `orpheus/transport/fields/_bases.py`
- `orpheus/transport/fields/angular_source_sink.py`
- `orpheus/transport/fields/harmonic_moment_flux.py`
- `orpheus/transport/fields/harmonic_moment_source_sink.py`
- `orpheus/transport/fields/scalar_boundary_flux.py`
- `orpheus/transport/fields/scalar_flux.py`
- `orpheus/transport/kernels.py`
- `orpheus/transport/mesh/material_xs_field.py`
- `orpheus/transport/operators/isotropic_transfer.py`
- `orpheus/transport/operators/transfer.py`
- `orpheus/transport/spatial/_ubld.py`

#### D-D.3 Hard test files

- `tests/gates/numerics/test_heads_on_the_axis.py`
- `tests/gates/numerics/test_moment_head_axis_built_premise.py`
- `tests/gates/sn/_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py`
- `tests/gates/sn/mesh/test_harmonic_moment_flux_on_the_problem.py`
- `tests/gates/sn/mesh/test_moment_space_is_read_off_the_frame_on_the_problem.py`
- `tests/gates/sn/mesh/test_spatial_moment_tail_on_the_problem.py`
- `tests/gates/sn/operators/test_capability_survival.py`
- `tests/gates/sn/operators/test_frame_conjugate_carve.py`
- `tests/gates/sn/operators/test_legendre_moment_scattering.py`
- `tests/gates/sn/operators/test_n2n_operator.py`
- `tests/gates/sn/operators/test_scattering_adjoint.py`
- `tests/gates/sn/operators/test_scattering_operator.py`
- `tests/gates/sn/primitives/test_harmonic_moment_flux.py`
- `tests/gates/sn/primitives/test_typed_source_sinks.py`
- `tests/gates/transport/_wave2_s4_helpers.py`
- `tests/gates/transport/fields/test_harmonic_moment_flux.py`
- `tests/gates/transport/fields/test_scalar_boundary_flux.py`
- `tests/gates/transport/test_admission_mechanism.py`
- `tests/gates/transport/test_binding.py`
- `tests/gates/transport/test_kernels.py`
- `tests/gates/transport/test_tier2_equivalence_s_family.py`

#### D-D.4 Hard docs files

- `docs/theory/conventions/indexing_and_layout.rst`
- `docs/theory/foundations/boundary_conditions.rst`
- `docs/theory/foundations/field_algebra.rst`
- `docs/theory/foundations/frame.rst`
- `docs/theory/foundations/infinite_medium.rst`
- `docs/theory/foundations/operator_inverse_family.rst`
- `docs/theory/foundations/spaces.rst`
- `docs/theory/methods/sn/adjoint.rst`
- `docs/theory/methods/sn/cartesian_multid.rst`
- `docs/theory/methods/sn/slab_one_group.rst`
- `docs/theory/verification/error_catalog.rst`

#### D-D.5 Suspect test files (opened only when the retired owner class is named; homonyms left alone)

- `tests/_harness/audit.py`
- `tests/gates/cross_method/test_admission_matrix_s4.py`
- `tests/gates/cross_method/test_one_parse_on_production.py`
- `tests/gates/derivations/test_dsa_production_tie.py`
- `tests/gates/diffusion/test_solver.py`
- `tests/gates/homogeneous/test_coda_anchors.py`
- `tests/gates/homogeneous/test_operator_spaces.py`
- `tests/gates/numerics/_energy_pencil.py`
- `tests/gates/numerics/test_outer_dyad.py`
- `tests/gates/sn/acceleration/test_dsa_low_order.py`
- `tests/gates/sn/acceleration/test_dsa_rate.py`
- `tests/gates/sn/architecture/_config.py`
- `tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py`
- `tests/gates/sn/operators/test_fission_adjoint.py`
- `tests/gates/sn/operators/test_fission_kernel_crosscheck.py`
- `tests/gates/sn/operators/test_isotropic_fission.py`
- `tests/gates/sn/operators/test_isotropic_scattering.py`
- `tests/gates/sn/operators/test_solver_components.py`
- `tests/gates/sn/sweep/curvilinear/test_cyl_sweep_regression.py`
- `tests/gates/sn/sweep/curvilinear/test_sph_sweep_regression.py`
- `tests/gates/test_pending_ports.py`
- `tests/gates/test_withdrawal.py`
- `tests/gates/tools/test_harness_generator.py`
- `tests/gates/transport/fields/test_angular_boundary_flux.py`
- `tests/gates/transport/spatial/test_ld_ubld_primitive.py`
- `tests/gates/transport/test_admission_matrix.py`
- `tests/gates/transport/test_energy_leaves_as_matrix.py`

#### D-D.6 Tense table (every prose line naming one of the sub-unit's retired names)

- `docs/theory/conventions/indexing_and_layout.rst` (text): PRESENT 1565; role 1231,1276,1277
- `docs/theory/foundations/boundary_conditions.rst` (text): PRESENT 2813; role 3862
- `docs/theory/foundations/field_algebra.rst` (text): role 1783,1784
- `docs/theory/foundations/frame.rst` (text): PRESENT 4843,5153; role 4180,4181,4239,4240,4836,4837,4839,4840,4854,6357; past 4176
- `docs/theory/foundations/infinite_medium.rst` (text): role 1956
- `docs/theory/foundations/operator_algebra.rst` (text): PRESENT 6590,6649
- `docs/theory/foundations/operator_inverse_family.rst` (text): role 1741
- `docs/theory/foundations/spaces.rst` (text): PRESENT 5013,5082,5083; role 1867,5003
- `docs/theory/methods/diffusion_1d.rst` (text): PRESENT 364
- `docs/theory/methods/sn/acceleration.rst` (text): PRESENT 1657; role 1621
- `docs/theory/methods/sn/adjoint.rst` (text): role 926,930,969
- `docs/theory/methods/sn/cartesian_multid.rst` (text): PRESENT 3194,3451; role 2026
- `docs/theory/methods/sn/history.rst` (text): past 1130,3398
- `docs/theory/methods/sn/loss_representation.rst` (text): role 4058,4107
- `docs/theory/methods/sn/slab_one_group.rst` (text): role 907,908
- `docs/theory/verification/error_catalog.rst` (text): past 5748,5769,6952; collision 11436
- `docs/theory/verification/reference_cache.rst` (text): collision 1127
- `orpheus/numerics/basis/base.py` (docstring): role 182
- `orpheus/numerics/elementary_operators.py` (docstring): collision 849
- `orpheus/numerics/face_layout.py` (docstring): role 217
- `orpheus/numerics/frame.py` (docstring): role 389
- `orpheus/numerics/spaces/moment_head.py` (docstring): PRESENT 18
- `orpheus/sn/acceleration/dsa.py` (docstring): role 86
- `orpheus/sn/operators/streaming.py` (docstring): PRESENT 132
- `orpheus/transport/boundary_conditions/admission.py` (docstring): PRESENT 12,42; role 49
- `orpheus/transport/fields/_bases.py` (docstring): PRESENT 39; role 1209,1263
- `orpheus/transport/fields/angular_boundary_flux.py` (docstring): role 58,122
- `orpheus/transport/fields/angular_boundary_source_sink.py` (docstring): role 70,188,233
- `orpheus/transport/fields/harmonic_moment_flux.py` (docstring): PRESENT 210,289; role 209
- `orpheus/transport/fields/scalar_boundary_flux.py` (docstring): role 88
- `orpheus/transport/fields/scalar_boundary_source_sink.py` (docstring): role 82
- `orpheus/transport/fields/scalar_flux.py` (docstring): role 126
- `orpheus/transport/kernels.py` (docstring): role 219,250; collision 57,224
- `orpheus/transport/mesh/material_xs_field.py` (docstring): PRESENT 34,35,676; role 12,13,14,15,72,608
- `orpheus/transport/operators/isotropic_transfer.py` (docstring): role 272
- `orpheus/transport/operators/transfer.py` (docstring): PRESENT 788,789; role 246,757,781,808
- `orpheus/transport/residuals/angular_boundary_residual.py` (docstring): role 102
- `orpheus/transport/spatial/_ubld.py` (docstring): PRESENT 433; role 410
- `tests/gates/numerics/test_matrix_inverse_operator.py` (docstring): PRESENT 7
- `tests/gates/numerics/test_moment_head_axis_built_premise.py` (docstring): role 359
- `tests/gates/numerics/test_traced_memo_template.py` (docstring): collision 968,971
- `tests/gates/sn/acceleration/test_dsa_rate.py` (docstring): PRESENT 1164,1165,1187
- `tests/gates/sn/architecture/_config.py` (docstring): PRESENT 314
- `tests/gates/sn/mesh/test_harmonic_moment_flux_on_the_problem.py` (docstring): PRESENT 11,19; role 3,4
- `tests/gates/sn/mesh/test_moment_space_is_read_off_the_frame_on_the_problem.py` (docstring): PRESENT 8
- `tests/gates/sn/operators/test_collision_operator.py` (docstring): past 5
- `tests/gates/sn/operators/test_removal_form_matvec_sweep.py` (docstring): role 495,496
- `tests/gates/sn/operators/test_scattering_operator.py` (docstring): PRESENT 744,819,895,938,1053,1122,1126,1140,1143,1151
- `tests/gates/sn/operators/test_streaming_operator.py` (docstring): PRESENT 4
- `tests/gates/sn/primitives/test_harmonic_moment_flux.py` (docstring): PRESENT 83,302,303,746,747,748
- `tests/gates/sn/primitives/test_typed_source_sinks.py` (docstring): past 442
- `tests/gates/transport/_carrier_helpers.py` (docstring): PRESENT 20
- `tests/gates/transport/fields/test_angular_boundary_flux.py` (docstring): PRESENT 16
- `tests/gates/transport/fields/test_angular_boundary_source_sink_residual.py` (docstring): PRESENT 15
- `tests/gates/transport/fields/test_harmonic_moment_flux.py` (docstring): PRESENT 11,19,126,189; role 3,4
- `tests/gates/transport/frames/test_moment_space_is_read_off_the_frame.py` (docstring): PRESENT 8
- `tests/gates/transport/test_admission_mechanism.py` (docstring): PRESENT 13,738,742
- `tests/gates/transport/test_binding.py` (docstring): PRESENT 227
- `tests/gates/transport/test_energy_leaves_as_matrix.py` (docstring): PRESENT 62
- `tests/gates/transport/test_kernels.py` (docstring): PRESENT 257
### Appendix R: unit roles

#### R-1 The five dead-table rows owned by roles

| row | file : symbol | verdict | destination or rewire | test fate (per the review) | markers | production references (file:line kind in scope) | closure | test files naming it |
|---|---|---|---|---|---|---|---|---|
| 64 | `numerics/operator.py` `BulkOperator` | retire | Table 2 | see Table 2 verdicts | - | numerics/operator.py:133 str in  | - | numerics/test_operator_protocols.py, sn/operators/test_capability_survival.py, sn/operators/test_operator_block_role.py, sn/operators/test_sn_boundary_operator.py |
| 65 | `numerics/operator.py` `FullOperator` | retire | Table 2 | see Table 2 verdicts | - | numerics/operator.py:134 str in  | - | diffusion/test_operators.py, numerics/test_operator_protocols.py, sn/operators/test_capability_survival.py, sn/operators/test_operator_block_role.py, sn/operators/test_sn_boundary_operator.py |
| 66 | `numerics/operator.py` `BoundaryOperator` | retire | Table 2 | see Table 2 verdicts | - | numerics/operator.py:135 str in  | - | diffusion/test_boundary_realizer.py, diffusion/test_operators.py, numerics/test_operator_protocols.py, sn/operators/test_capability_survival.py, sn/operators/test_operator_block_role.py, sn/operators/test_sn_boundary_operator.py |
| 144 | `transport/boundary_conditions/realization.py` `stamp_boundary_role` | retire | Table 2 | test_wave2_retirement.py:109 and test_layer_imports.py:667-668 string fixtures | - | diffusion/boundary_realizer.py:127 import in ; diffusion/boundary_realizer.py:169 name in _albedo_operator; sn/boundary/realized_operators.py:125 import in ; sn/boundary/realized_operators.py:154 name in _attenuated_kernel_operator; sn/boundary/realized_operators.py:384 name in ZeroMap._derive_inner; transport/boundary_conditions/__init__.py:46 import in ; transport/boundary_conditions/__init__.py:107 str in ; transport/boundary_conditions/realization.py:114 str in  | - | sn/operators/test_capability_survival.py, test_wave2_retirement.py, transport/test_boundary_conditions_package.py |
| 145 | `numerics/operator.py` `_BlockRoleMeta` | retire | Table 2 | - | - | numerics/operator.py:459 name in BulkOperator; numerics/operator.py:465 name in FullOperator; numerics/operator.py:471 name in BoundaryOperator | numerics/operator.py:453 _BlockRoleMeta._role, numerics/operator.py:462 BulkOperator._role, numerics/operator.py:468 FullOperator._role, numerics/operator.py:507 BoundaryOperator._role | - |

#### R-2 Hard files (docs, production, tests; the production and test lists are in 4.1 by symbol)

- `docs/api/numerics.rst`
- `docs/development/evidence/lessons.md`
- `docs/theory/foundations/boundary_conditions.rst`
- `docs/theory/foundations/coupled_block_operator.rst`
- `docs/theory/foundations/operator_algebra.rst`
- `docs/theory/foundations/operator_tensor_network.rst`
- `docs/theory/methods/diffusion_1d.rst`
- `docs/theory/verification/error_catalog.rst`
- `orpheus/diffusion/boundary_realizer.py`
- `orpheus/diffusion/operators.py`
- `orpheus/numerics/coupled_system.py`
- `orpheus/numerics/operator.py`
- `orpheus/sn/boundary/realized_operators.py`
- `orpheus/sn/coupled_system.py`
- `orpheus/sn/operators/boundary.py`
- `orpheus/sn/operators/radial_characteristic.py`
- `orpheus/sn/operators/scheduled_invertible.py`
- `orpheus/sn/operators/streaming.py`
- `orpheus/transport/boundary_conditions/__init__.py`
- `orpheus/transport/boundary_conditions/realization.py`
- `orpheus/transport/operators/angular_lift.py`
- `orpheus/transport/operators/isotropic_transfer.py`
- `orpheus/transport/operators/lift.py`
- `orpheus/transport/operators/multiplication_operator.py`
- `orpheus/transport/operators/transfer.py`
- `tests/gates/diffusion/test_augmented_mesh.py`
- `tests/gates/diffusion/test_boundary_realizer.py`
- `tests/gates/diffusion/test_diffusion_realizer_table.py`
- `tests/gates/diffusion/test_operators.py`
- `tests/gates/numerics/_assembly_functor.py`
- `tests/gates/numerics/test_coupled_operator.py`
- `tests/gates/numerics/test_linear_combination.py`
- `tests/gates/numerics/test_operator_protocols.py`
- `tests/gates/sn/architecture/test_realized_boundary.py`
- `tests/gates/sn/operators/test_capability_survival.py`
- `tests/gates/sn/operators/test_linear_combination_sn.py`
- `tests/gates/sn/operators/test_operator_block_role.py`
- `tests/gates/sn/operators/test_psi_half_coupling.py`
- `tests/gates/sn/operators/test_sn_boundary_operator.py`
- `tests/gates/sn/solve/test_si_single_primitive_contract.py`
- `tests/gates/sn/verification/analytical/test_mms_declared_inflow.py`
- `tests/gates/test_layer_imports.py`
- `tests/gates/test_wave2_retirement.py`
- `tests/gates/transport/test_boundary_conditions_package.py`
- `tests/gates/transport/test_bulk_lift.py`

#### R-2b Prose-only files

- `orpheus/diffusion/augmented_mesh.py`
- `orpheus/numerics/registry.py`
- `orpheus/numerics/spaces/angular_trace_space.py`
- `orpheus/sn/splitting.py`
- `orpheus/transport/__init__.py`
- `orpheus/transport/fields/angular_boundary_source_sink.py`
- `tests/gates/numerics/test_angular_trace_space.py`
- `tests/gates/numerics/test_face_streaming_normal.py`
- `tests/gates/sn/operators/test_bc_extraction_matvec.py`
- `tests/gates/sn/operators/test_declared_law_is_linear.py`
- `tests/gates/sn/operators/test_fission_operator.py`
- `tests/gates/sn/regression/README.md`
- `tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py`

#### R-2c New

- `tests/gates/sn/operators/test_sn_boundary_operator.py` (existing file; gains `test_prescribed_inflow_realizes_the_same_object_vacuum_does`, `_FACE` and `_boundary_under_test`)
- `tests/gates/_retired/__init__.py`, `tests/gates/_retired/_finder.py`, `tests/gates/_retired/rows_w0_roles.py`, `tests/gates/_retired/rows_wave2.py`, `tests/gates/test_retirements.py` (1.5)
- `tests/gates/numerics/test_linear_combination.py` gains the two moved assertions of `test_operator_protocols.py` (R2 paragraph)

#### R-3 Tense table (the 12 names; BoundaryOperator homonyms, retired-homonym lines and ERR-075 history stay)

- `.claude/skills/vv-principles/error_index.md` (text): past 126
- `docs/api/numerics.rst` (text): PRESENT 346
- `docs/theory/foundations/boundary_conditions.rst` (text): PRESENT 4536; role 1514,4125,4541,4761,4913,4922,5284,7072; past 3436,4166,4188,7615,8055,8271,10525
- `docs/theory/foundations/coupled_block_operator.rst` (text): PRESENT 839,1224; role 117,118,142,821,824,833
- `docs/theory/foundations/operator_algebra.rst` (text): PRESENT 2954
- `docs/theory/foundations/operator_tensor_network.rst` (text): PRESENT 137,1125
- `docs/theory/methods/diffusion_1d.rst` (text): PRESENT 373,416
- `docs/theory/verification/error_catalog.rst` (text): past 1378,1386,6340,6369
- `orpheus/diffusion/augmented_mesh.py` (docstring): role 308
- `orpheus/diffusion/boundary_realizer.py` (docstring): role 24,162
- `orpheus/numerics/operator.py` (docstring): PRESENT 436,441; role 227,262,263,296,297,438,439,460,466,472,481,482
- `orpheus/numerics/registry.py` (docstring): PRESENT 52,56; role 55
- `orpheus/numerics/spaces/angular_trace_space.py` (docstring): PRESENT 50
- `orpheus/sn/boundary/realized_operators.py` (docstring): role 80,146
- `orpheus/sn/coupled_system.py` (docstring): PRESENT 38,614
- `orpheus/sn/operators/boundary.py` (docstring): role 42,217
- `orpheus/transport/__init__.py` (docstring): PRESENT 58,59
- `orpheus/transport/boundary_conditions/realization.py` (docstring): PRESENT 69; role 159
- `orpheus/transport/fields/angular_boundary_source_sink.py` (docstring): PRESENT 103
- `tests/gates/diffusion/test_boundary_realizer.py` (docstring): PRESENT 15,140
- `tests/gates/diffusion/test_diffusion_realizer_table.py` (docstring): PRESENT 165,296
- `tests/gates/numerics/test_coupled_operator.py` (docstring): PRESENT 36
- `tests/gates/numerics/test_face_streaming_normal.py` (docstring): PRESENT 7
- `tests/gates/numerics/test_operator_protocols.py` (docstring): role 4,6,7,8,19; past 20
- `tests/gates/sn/architecture/test_realized_boundary.py` (docstring): PRESENT 360
- `tests/gates/sn/operators/test_bc_extraction_matvec.py` (docstring): PRESENT 35
- `tests/gates/sn/operators/test_capability_survival.py` (docstring): PRESENT 148; role 129
- `tests/gates/sn/operators/test_declared_law_is_linear.py` (docstring): PRESENT 14
- `tests/gates/sn/operators/test_operator_block_role.py` (docstring): PRESENT 21,23,191; role 6,8,10,29,176
- `tests/gates/sn/operators/test_psi_half_coupling.py` (docstring): PRESENT 3015,3135,3462; role 3018
- `tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py` (docstring): role 9
- `tests/gates/test_wave2_retirement.py` (docstring): PRESENT 31
- `tests/gates/transport/test_boundary_conditions_package.py` (docstring): PRESENT 21, 202 (line 21 is the bullet of the deleted `test_stamp_boundary_role_is_copied_in` and reads "lives here": present tense, not past)

#### R-3b Lines the tense table missed (`[M]` by reading each at the quoted line, all at `b68f8b8c`; the symbol decides)

Docs (13 lines, 12 to rewrite and 1 homonym): `docs/theory/foundations/boundary_conditions.rst` about lines 4140, 4146, 4164 (a table cell), 4187, 4539, 5286, 5291 (all name `block_role` in PRESENT tense or as a cell); `coupled_block_operator.rst` about lines 427 ("`system_role` tags are runtime") and 837; `operator_algebra.rst` about line 2955 (`head.system_role`); `docs/theory/verification/error_catalog.rst` about line 6381 (ERR-075 history: "no `block_role` test anywhere", a past-tense entry: STAYS, read); `docs/development/evidence/lessons.md` about line 825 (a lessons line: edited at its source only if the page is generated, else in place); `tests/gates/sn/regression/README.md` about line 167 (`BoundaryOperator` here is the retired Wave-7 homonym: STAYS). Docstring lines of 4.1 (2) are the production and test docstring part of the same census gap (7 lines). Comments are not censused (4.1 (4)).


#### R-4 The 55 role-tag test functions, with the review's fate (delete 42 = 41 test functions and the double `_Tagged.__init__`, rewire 8, move 2, edit 2 doubles, keep and re-home 1)

| test (under `tests/gates/`) : line at base | scope | role assertions / all | fate |
|---|---|---|---|
| `diffusion/test_augmented_mesh.py:257` | `TestTheBoundaryStamp.test_every_realized_face_operator_is_a_boundary_leaf_mixed_wall_included` | 1 / 1 | delete |
| `diffusion/test_boundary_realizer.py:139` | `TestAlbedoTable.test_every_realized_operator_is_boundary_block` | 2 / 2 | delete |
| `diffusion/test_diffusion_realizer_table.py:159` | `test_diffusion_vacuum_is_tier_1` | 3 / 10 | rewire: drop the tag assertions |
| `diffusion/test_diffusion_realizer_table.py:291` | `test_diffusion_translation_is_topological` | 1 / 8 | rewire: drop the tag assertions |
| `diffusion/test_operators.py:503` | `TestFamilyLaws.test_loss_block_role_joins_to_full` | 2 / 2 | delete |
| `diffusion/test_operators.py:523` | `TestLeakageOperator.test_block_role_full` | 1 / 1 | delete |
| `diffusion/test_operators.py:622` | `TestDiffusionBoundaryOperator.test_block_role_boundary` | 2 / 2 | delete |
| `numerics/test_coupled_operator.py:935` | `TestSystemRoleStamp.test_coupled_operator_is_stamped_coupled` | 2 / 2 | delete |
| `numerics/test_linear_combination.py:127` | `Mat.__init__` | 0 / 0 | edit: drop the `block_role` parameter of the double |
| `numerics/test_linear_combination.py:811` | `test_c31_roles_fold_over_the_terms_whatever_the_order` | 3 / 3 | delete with the join (its mutation ROLE_HEAD has no target left) |
| `numerics/test_linear_combination.py:872` | `test_c42_a_role_lives_on_the_leaves_and_survives_every_rebuild` | 2 / 2 | delete with the join; the ruling's repeal belongs in the commit |
| `numerics/test_operator_protocols.py:52` | `_Tagged.__init__` | 0 / 0 | delete with its class (it exists to carry a tag) |
| `numerics/test_operator_protocols.py:57` | `TestBlockRoleEnum.test_enum_has_the_three_roles` | 1 / 1 | delete |
| `numerics/test_operator_protocols.py:60` | `TestBlockRoleEnum.test_roles_are_distinct` | 1 / 1 | delete |
| `numerics/test_operator_protocols.py:69` | `TestSystemRoleEnum.test_enum_has_the_three_roles` | 1 / 1 | delete |
| `numerics/test_operator_protocols.py:72` | `TestSystemRoleEnum.test_roles_are_distinct` | 1 / 1 | delete |
| `numerics/test_operator_protocols.py:77` | `TestIsinstanceMarkers.test_bulk_role_matches_bulk_marker_only` | 3 / 3 | delete |
| `numerics/test_operator_protocols.py:83` | `TestIsinstanceMarkers.test_full_role_matches_full_marker_only` | 3 / 3 | delete |
| `numerics/test_operator_protocols.py:89` | `TestIsinstanceMarkers.test_boundary_role_matches_boundary_marker_only` | 3 / 3 | delete |
| `numerics/test_operator_protocols.py:95` | `TestIsinstanceMarkers.test_none_role_matches_no_marker` | 3 / 3 | delete |
| `numerics/test_operator_protocols.py:101` | `TestIsinstanceMarkers.test_object_without_block_role_matches_no_marker` | 3 / 3 | delete |
| `numerics/test_operator_protocols.py:107` | `TestIsinstanceMarkers.test_isinstance_is_value_based_not_inheritance` | 2 / 2 | delete |
| `numerics/test_operator_protocols.py:114` | `TestIsinstanceMarkers.test_classification_is_a_partition` | 0 / 1 | delete |
| `numerics/test_operator_protocols.py:128` | `TestGenericOperatorsAreUnclassified.test_identity_operator_is_unclassified` | 4 / 4 | delete |
| `numerics/test_operator_protocols.py:134` | `TestGenericOperatorsAreUnclassified.test_zero_operator_is_unclassified` | 1 / 1 | delete |
| `numerics/test_operator_protocols.py:137` | `TestGenericOperatorsAreUnclassified.test_operator_sum_is_unclassified_at_O1` | 4 / 5 | MOVE to `numerics/test_linear_combination.py` as `test_a_sum_of_the_identity_and_a_diagonal_keeps_both_terms_as_an_operator_sum`, tag assertions dropped (the file is deleted whole) |
| `numerics/test_operator_protocols.py:157` | `TestGenericOperatorsAreUnclassified.test_tensor_product_operator_is_unclassified_by_default` | 4 / 4 | delete |
| `numerics/test_operator_protocols.py:166` | `TestGenericOperatorsAreUnclassified.test_sum_of_tensor_products_is_unclassified_by_default` | 4 / 5 | MOVE to `numerics/test_linear_combination.py` as `test_a_sum_of_two_tensor_products_is_an_operator_sum_of_two_terms`, tag assertions dropped (the file is deleted whole) |
| `sn/architecture/test_realized_boundary.py:357` | `test_the_algebra_object_delegates` | 1 / 8 | rewire: drop the tag assertions |
| `sn/operators/test_capability_survival.py:137` | `TestTheBlockRoleIsAnInstanceStampNotAClassFact.test_an_unstamped_instance_is_unclassified` | 4 / 4 | delete |
| `sn/operators/test_capability_survival.py:144` | `TestTheBlockRoleIsAnInstanceStampNotAClassFact.test_stamping_one_instance_leaves_the_class_alone` | 4 / 4 | delete |
| `sn/operators/test_linear_combination_sn.py:55` | `Mat.__init__` | 0 / 0 | edit: drop the `block_role` parameter of the double |
| `sn/operators/test_operator_block_role.py:85` | `TestBulkLeaves.test_collision_is_bulk` | 3 / 3 | delete |
| `sn/operators/test_operator_block_role.py:93` | `TestBulkLeaves.test_scattering_is_bulk` | 2 / 2 | delete |
| `sn/operators/test_operator_block_role.py:97` | `TestBulkLeaves.test_fission_is_bulk` | 2 / 2 | delete |
| `sn/operators/test_operator_block_role.py:103` | `TestFullLeaves.test_streaming_is_full` | 3 / 3 | delete |
| `sn/operators/test_operator_block_role.py:111` | `TestFullLeaves.test_invertible_L_plus_C_is_full` | 3 / 3 | delete |
| `sn/operators/test_operator_block_role.py:179` | `TestBoundaryLeaves.test_realized_bc_advertises_boundary_operator` | 4 / 4 | delete |
| `sn/operators/test_operator_block_role.py:187` | `TestBoundaryLeaves.test_prescribed_inflow_realizes_the_same_object_vacuum_does` | 0 / 4 | keep, re-home (its file is otherwise deleted) |
| `sn/operators/test_operator_block_role.py:218` | `TestBoundaryLeaves.test_mesh_bc_forwards_boundary_role` | 3 / 3 | delete |
| `sn/operators/test_operator_block_role.py:259` | `TestComposerRoleDerivation.test_scaled_preserves_role` | 3 / 3 | delete |
| `sn/operators/test_operator_block_role.py:266` | `TestComposerRoleDerivation.test_sum_join_is_union_of_blocks` | 5 / 5 | delete |
| `sn/operators/test_operator_block_role.py:281` | `TestComposerRoleDerivation.test_adjoint_preserves_role` | 6 / 6 | delete |
| `sn/operators/test_operator_block_role.py:293` | `TestComposerRoleDerivation.test_unclassified_operand_propagates_none` | 2 / 2 | delete |
| `sn/operators/test_psi_half_coupling.py:3027` | `TestSystemRoleLattice.test_join_is_the_two_system_union` | 10 / 10 | delete |
| `sn/operators/test_psi_half_coupling.py:3043` | `TestSystemRoleLattice.test_psi_half_blocks_carry_their_system_role` | 4 / 4 | delete |
| `sn/operators/test_psi_half_coupling.py:3054` | `TestSystemRoleLattice.test_model_generic_operators_stay_unclassified` | 4 / 4 | delete |
| `sn/operators/test_psi_half_coupling.py:3063` | `TestSystemRoleLattice.test_role_propagates_through_the_composers` | 3 / 3 | delete |
| `sn/operators/test_psi_half_coupling.py:3129` | `TestCoupledBuilder.test_p1_alignment_by_construction_and_runtime_apply` | 2 / 12 | rewire: drop the tag assertions |
| `sn/operators/test_psi_half_coupling.py:3458` | `TestWithinGroupSystem.test_g_d1_2_record_container_and_identity_pins` | 1 / 20 | rewire: drop the tag assertions |
| `sn/operators/test_sn_boundary_operator.py:201` | `TestContract.test_block_role_is_boundary_and_exclusive` | 4 / 4 | delete |
| `sn/solve/test_si_single_primitive_contract.py:104` | `test_fixed_source_si_and_eigenvalue_inner_share_one_primitive` | 3 / 26 | rewire: drop the tag assertions |
| `sn/verification/analytical/test_mms_declared_inflow.py:497` | `test_B_is_the_zero_morphism_when_every_face_is_prescribed` | 2 / 3 | rewire: drop the 2 tag lines (:523, :527) |
| `transport/test_boundary_conditions_package.py:201` | `test_stamp_boundary_role_is_copied_in` | 3 / 3 | delete |
| `transport/test_bulk_lift.py:163` | `TestLaws.test_capability_axes_and_role` | 1 / 4 | rewire: drop the tag assertions; rename `test_capability_axes` |

#### D-A.7 The rows module `tests/gates/_retired/rows_w0_dead_a.py` (data of the shared gate, 1.5)

Row types are `Retired(file, qualname)`, `Moved(file, qualname, home)` and `RetiredName(name)` of `tests/gates/_retired/__init__.py`; every `Moved` row carries its concrete home and, in its comment, the free function it becomes (decision D2: the object is the first argument and every call site is rewritten); `KEPT` rows are asserted still defined. The sub-unit's rows are listed in the order of its appendix table (`#### D-A.1`).

```python
from tests.gates._retired import Moved, Retired, RetiredName

ROWS = (
    Moved('orpheus/geometry/chart.py', 'DirectionDomain.tangencies', 'tests/gates/geometry/_geometry_oracles.py'),  # row 6; free function tangencies(direction_domain, level)
    Retired('orpheus/geometry/chord.py', 'Chord.slot_start'),
    Retired('orpheus/geometry/coord.py', 'compute_volumes_1d'),
    Moved('orpheus/geometry/manifold.py', 'UNIT_INTERVAL', 'tests/gates/geometry/_manifold_fixtures.py'),  # row 27; a value
    Retired('orpheus/geometry/manifold.py', 'barycentre'),
    Retired('orpheus/geometry/structured_geometry.py', 'StructuredGeometry.domain_extent_cm'),
    Retired('orpheus/geometry/symmetry.py', 'SubgroupOfO3.identity_component'),
    Moved('orpheus/geometry/symmetry.py', '_group_elements', 'tests/gates/geometry/_geometry_oracles.py'),  # row 32; free function _group_elements(tag)
    Moved('orpheus/geometry/symmetry.py', 'OrbitCertificate.stabilizer_order', 'tests/gates/geometry/_geometry_oracles.py'),  # row 33; free function stabilizer_order(orbit_certificate)  # was a property
    Moved('orpheus/geometry/transformation.py', 'RigidMotion.element_order', 'tests/gates/geometry/test_transformation.py'),  # row 34; free function element_order(rigid_motion, *, atol)
    Moved('orpheus/geometry/transformation.py', 'RigidMotion.permutes', 'tests/gates/geometry/_geometry_oracles.py'),  # row 36; free function permutes(rigid_motion, points, *, atol)
    Retired('orpheus/numerics/basis/descent.py', 'Descent.discriminator'),
    Retired('orpheus/numerics/basis/descent.py', 'Descent.descending_real_slots'),
    Retired('orpheus/numerics/basis/descent.py', 'Descent.is_isomorphism'),
    Retired('orpheus/numerics/elementary_operators.py', 'DiagonalOperator.from_measure'),
    Moved('orpheus/numerics/generating_measure.py', 'HERMITE', 'tests/gates/numerics/_measure_fixtures.py'),  # row 54; a value
    Moved('orpheus/numerics/generating_measure.py', 'laguerre', 'tests/gates/numerics/_measure_fixtures.py'),  # row 55; free function laguerre(a)
    Moved('orpheus/numerics/index_set.py', 'ENERGY', 'tests/gates/numerics/_measure_fixtures.py'),  # row 56; a value
    Moved('orpheus/numerics/measure.py', 'equispaced', 'tests/gates/numerics/_measure_fixtures.py'),  # row 60; free function equispaced(a, b, n)
    Moved('orpheus/numerics/mesh_free_function.py', 'Symbolic.from_srepr', 'tests/gates/numerics/_measure_fixtures.py'),  # row 61; free function from_srepr(symbolic, srepr, sympy_version)  # was a classmethod: first argument the class
    Retired('orpheus/numerics/quadrature/registry.py', 'SelectionLog.summary'),
    Retired('orpheus/numerics/quadrature/registry.py', 'select_quadrature'),
    Retired('orpheus/numerics/quadrature/registry.py', 'SelectionLog'),
    Retired('orpheus/numerics/quadrature/registry.py', 'QuadratureSelectionError'),
    Retired('orpheus/numerics/quadrature/registry.py', 'QuadratureSpec.structural_flags'),  # not a row: the ruled closure member of row 73 (3.0 item 3, 8.24)
)

KEPT = (  # (file, qualname, kept_for); asserted still defined
    ('orpheus/geometry/chart.py', 'Chart.measure_density', 'derivations/ (A2): characteristic/basis.py:285, closure.py:284'),  # row 3
    ('orpheus/geometry/chart.py', 'Chart.directions_at', 'derivations/ (A2): characteristic/reading.py:176,202,276'),  # row 4
    ('orpheus/geometry/chart.py', 'Chart.line_domain', 'derivations/ (A2): characteristic/assembly.py, lines.py, reading.py'),  # row 5
    ('orpheus/geometry/chart.py', 'RadialImage.orbit_coordinate_at', 'derivations/ (A2) via Chord/Rays.orbit_coordinate_at'),  # row 7
    ('orpheus/geometry/chart.py', 'AxialImage.orbit_coordinate_at', 'as row 7'),  # row 8
    ('orpheus/geometry/chord.py', 'ConcentricPartition.region_at', "W6, concept row 'Point location, chord' (#596: ConcentricPartition is the one locator)"),  # row 9
    ('orpheus/geometry/chord.py', 'ConcentricPartition.refined_at', 'derivations/ (A2): characteristic/reading.py:204,280'),  # row 10
    ('orpheus/geometry/chord.py', 'ConcentricPartition.coarsening', 'derivations/ (A2): characteristic/basis.py:143, reading.py:206, transport.py:679'),  # row 11
    ('orpheus/geometry/chord.py', 'ConcentricPartition.rays_from', 'derivations/ (A2): peierls geometry.py, fn sphere, characteristic'),  # row 12
    ('orpheus/geometry/chord.py', 'Chord.orbit_coordinate_at', 'derivations/ (A2)'),  # row 14
    ('orpheus/geometry/chord.py', 'Chord.rays_at', 'derivations/ (A2): characteristic/transport.py:453'),  # row 15
    ('orpheus/geometry/chord.py', 'Chord.rays_at_closest_approach', 'derivations/ (A2): characteristic/transport.py:456'),  # row 16
    ('orpheus/geometry/chord.py', 'Chord.line_integral', 'derivations/ (A2): 6 reads'),  # row 17
    ('orpheus/geometry/chord.py', 'Chord.cosine_at', 'derivations/ (A2): peierls geometry.py (6)'),  # row 18
    ('orpheus/geometry/chord.py', 'Rays.lengths_to', 'derivations/ (A2): peierls geometry.py:910,1181,1263'),  # row 19
    ('orpheus/geometry/chord.py', 'Rays.to_wall', 'derivations/ (A2): 6 reads'),  # row 20
    ('orpheus/geometry/chord.py', 'Rays.orbit_coordinate_at', 'derivations/ (A2)'),  # row 21
    ('orpheus/geometry/chord.py', 'Rays.advanced', 'derivations/ (A2): characteristic/transport.py:459'),  # row 22
    ('orpheus/geometry/coord.py', 'CoordSystem.measure_density', 'derivations/ (A2) through Chart.measure_density'),  # row 23
    ('orpheus/geometry/deck.py', 'ISOMETRY_CLASSES_BY_DIMENSION', 'W6 (declared deferral table, read with DEFERRED_DECKS)'),  # row 25
    ('orpheus/geometry/deck.py', 'DEFERRED_DECKS', "W6: SCOPE-BOUNDARY table named by production's DeckError message (geometry/deck.py:293)"),  # row 26
    ('orpheus/geometry/structured_geometry.py', 'StructuredGeometry.is_hollow', 'derivations/ (A2): common/reference_body.py:164,188'),  # row 29
    ('orpheus/geometry/transformation.py', 'RigidMotion.seated_at', "W6, concept row 'Symmetry group, translations included' (a subgroup of E(d))"),  # row 35
    ('orpheus/numerics/measure.py', 'DiscreteMeasure.singular_set_under', "W8a, concept row 'Discrete measure' (the invariance half of the weld)"),  # row 57
    ('orpheus/numerics/measure.py', 'DiscreteMeasure.symmetry_groups', 'W8a, as row 57'),  # row 58
    ('orpheus/numerics/measure.py', 'BundleMeasure', "W8a, concept row 'Discrete measure' ('widened to differing supports')"),  # row 59
    ('orpheus/numerics/quadrature/registry.py', '_gl1d_invert', 'W8a: a value of the required field degree_of_exactness_for of the kept catalogue (3.0 item 2)'),  # row 134
    ('orpheus/numerics/quadrature/registry.py', '_lebedev_invert', 'as row 134'),  # row 135
    ('orpheus/numerics/quadrature/registry.py', '_ls_sn_invert', 'as row 134'),  # row 136
    ('orpheus/numerics/quadrature/registry.py', '_product_invert', 'as row 134'),  # row 137
)
```

#### D-B.7 The rows module `tests/gates/_retired/rows_w0_dead_b.py` (data of the shared gate, 1.5)

Row types are `Retired(file, qualname)`, `Moved(file, qualname, home)` and `RetiredName(name)` of `tests/gates/_retired/__init__.py`; every `Moved` row carries its concrete home and, in its comment, the free function it becomes (decision D2: the object is the first argument and every call site is rewritten); `KEPT` rows are asserted still defined. The sub-unit's rows are listed in the order of its appendix table (`#### D-B.1`).

```python
from tests.gates._retired import Moved, Retired, RetiredName

ROWS = (
    Retired('orpheus/data/emission_spectrum.py', 'EmissionSpectrum.is_emitting'),
    Retired('orpheus/data/macro_xs/mixture.py', 'Mixture.out_scattering_xs'),
    Retired('orpheus/diffusion/method_space.py', 'DiffusionMethodSpace.minimal'),
    Moved('orpheus/numerics/basis/overlap_basis.py', 'OverlapBasis.dominant_column', 'tests/gates/data/test_energy_grid.py'),  # row 42; free function dominant_column(overlap_basis)  # was a property
    Moved('orpheus/numerics/basis/overlap_basis.py', 'OverlapBasis.fractional_columns', 'tests/gates/numerics/_overlap_basis_oracles.py'),  # row 43; free function fractional_columns(overlap_basis)  # was a property
    Retired('orpheus/numerics/basis/spherical_harmonic_basis.py', 'SphericalHarmonicBasis.evaluate_from_components'),
    Retired('orpheus/numerics/basis/spherical_harmonic_basis.py', 'MirrorEvenSphericalHarmonicBasis.evaluate_from_components'),
    Moved('orpheus/numerics/coefficient_layout.py', 'TensorLayout.multi_index', 'tests/gates/numerics/test_interval_legendre_basis.py'),  # row 46; free function multi_index(tensor_layout, flat)
    Moved('orpheus/numerics/eigenvalue.py', 'direct_eigenvalue', 'tests/gates/numerics/_dense_eigen_oracle.py'),  # row 49; free function direct_eigenvalue(A, F, *, imag_tol)
    Retired('orpheus/numerics/eigenvalue.py', 'rayleigh_quotient_iteration'),
    Retired('orpheus/numerics/operator.py', 'InverseMetricOperator'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.is_regular'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.rhs_rank'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.apply_at'),
    Moved('orpheus/numerics/space.py', 'FunctionSpace.project_live', 'tests/gates/numerics/_space_oracles.py'),  # row 74; free function project_live(function_space, x)
    Retired('orpheus/numerics/space.py', 'angular_flux_space'),
    Retired('orpheus/numerics/space.py', 'scalar_flux_space'),
    Moved('orpheus/numerics/eigenvalue.py', 'dominant_eigenpair', 'tests/gates/numerics/_dense_eigen_oracle.py'),  # row 130; free function dominant_eigenpair(M, *, imag_tol)
    Moved('orpheus/numerics/eigenvalue.py', '_sign_normalised', 'tests/gates/numerics/_dense_eigen_oracle.py'),  # row 131; free function _sign_normalised(v)
)

KEPT = (  # (file, qualname, kept_for); asserted still defined
    ('orpheus/numerics/basis/base.py', 'Basis.synthesize', "W8a, concept row 'Basis' (one table contraction)"),  # row 38
    ('orpheus/numerics/coupled_system.py', 'CoupledSpace.system_slices', "W1, concept row 'Direct sum of spaces' (4 offset loops become one)"),  # row 47
    ('orpheus/numerics/coupled_system.py', 'CoupledSpace.kernel_basis', "W1/W2: the metric-kernel family (concept 'Quotient by a metric kernel', A3)"),  # row 48
    ('orpheus/numerics/field.py', 'Field.cone_violations', "concept row 'Positive cone' (one, landed); W7 'Wall law on its cone'"),  # row 53
    ('orpheus/numerics/metric.py', 'HilbertMetric.kernel_basis', 'as row 48'),  # row 62
    ('orpheus/numerics/posing.py', 'ALPHA_MAP', "W5, concept row 'Question bound to a system' (spectral map derived from the parameter)"),  # row 71
    ('orpheus/numerics/space.py', 'FunctionSpace.kernel_basis', 'as row 48'),  # row 75
    ('orpheus/numerics/field.py', 'Field.where_largest', 'deferred, private state (decision D2; section 8.3)'),  # row 52
)
```

#### D-C.7 The rows module `tests/gates/_retired/rows_w0_dead_c.py` (data of the shared gate, 1.5)

Row types are `Retired(file, qualname)`, `Moved(file, qualname, home)` and `RetiredName(name)` of `tests/gates/_retired/__init__.py`; every `Moved` row carries its concrete home and, in its comment, the free function it becomes (decision D2: the object is the first argument and every call site is rewritten); `KEPT` rows are asserted still defined. The sub-unit's rows are listed in the order of its appendix table (`#### D-C.1`).

```python
from tests.gates._retired import Moved, Retired, RetiredName

ROWS = (
    Moved('orpheus/mesh/axis.py', 'n_unknowns_flat', 'tests/gates/sn/primitives/test_axis_primitive.py'),  # row 37; free function n_unknowns_flat(axes, quad, ng)
    Moved('orpheus/sn/angular/closure.py', '_MMHalfGrid.downstream_per_ordinate', 'tests/gates/sn/sweep/curvilinear/_psi_half.py'),  # row 80; free function downstream_per_ordinate(m_m_half_grid)  # was a property
    Moved('orpheus/sn/angular/closure.py', '_MMHalfGrid.upstream', 'tests/gates/sn/sweep/curvilinear/_psi_half.py'),  # row 81; free function upstream(m_m_half_grid, m)
    Moved('orpheus/sn/angular/closure.py', 'compute_psi_half_per_level', 'tests/gates/sn/sweep/curvilinear/_psi_half.py'),  # row 82; free function compute_psi_half_per_level(psi_level, tau_level, *, psi_half_seed)
    Retired('orpheus/sn/coupled_system.py', 'build_coupled_system'),
    Moved('orpheus/sn/loss_representation/assembly.py', 'ordinate_walk_order', 'tests/gates/sn/sweep/_ordinate_assembly.py'),  # row 85; free function ordinate_walk_order(problem, ordinate)
    Moved('orpheus/sn/loss_representation/assembly.py', 'assemble_ordinate_blocks', 'tests/gates/sn/sweep/_ordinate_assembly.py'),  # row 86; free function assemble_ordinate_blocks(problem, ordinate, *, include_collision)
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_in_x', 'tests/gates/sn/sweep/_sweep_graph_faces.py'),  # row 87; free function face_in_x(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_out_x', 'tests/gates/sn/sweep/_sweep_graph_faces.py'),  # row 88; free function face_out_x(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_in_y', 'tests/gates/sn/sweep/_sweep_graph_faces.py'),  # row 89; free function face_in_y(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_out_y', 'tests/gates/sn/sweep/_sweep_graph_faces.py'),  # row 90; free function face_out_y(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/operators/loss_kernel_gauge.py', 'predicted_kernel_dimension', 'tests/gates/sn/operators/_kernel_dimension.py'),  # row 91; free function predicted_kernel_dimension(problem)
    Moved('orpheus/sn/problem.py', 'SNProblem.n_unknowns_flat', 'tests/gates/sn/primitives/test_axis_primitive.py'),  # row 92; free function n_unknowns_flat(s_n_problem)  # was a property
    Retired('orpheus/sn/problem.py', 'SNProblem.loss_kernel_gauge'),
    Retired('orpheus/sn/problem.py', 'SNProblem.boundary_face_layout'),
    Moved('orpheus/sn/solver.py', 'boundary_vs_interior_split', 'tests/gates/sn/operators/test_typed_residual_evaluation.py'),  # row 95; free function boundary_vs_interior_split(residual)
    Retired('orpheus/sn/solver.py', 'SNSolver.compute_group_absorption_rate'),
    Retired('orpheus/sn/splitting.py', 'LossTerm.is_gain'),
    Moved('orpheus/sn/splitting.py', 'Splitting.law_residual', 'tests/gates/sn/architecture/_splitting_oracles.py'),  # row 98; free function law_residual(splitting, state)
    Moved('orpheus/sn/loss_representation/assembly.py', '_ordinate_label', 'tests/gates/sn/sweep/_ordinate_assembly.py'),  # row 141; free function _ordinate_label(problem, ordinate)
    Moved('orpheus/sn/loss_representation/assembly.py', '_iter_cells_in_walk_order', 'tests/gates/sn/sweep/_ordinate_assembly.py'),  # row 142; free function _iter_cells_in_walk_order(graph)
    Moved('orpheus/sn/loss_representation/assembly.py', '_probe_coefficient_blocks', 'tests/gates/sn/sweep/_ordinate_assembly.py'),  # row 143; free function _probe_coefficient_blocks(scheme, s_axes_k, reaction_flat, *, cm, fm, n_cells)
)

KEPT = (  # (file, qualname, kept_for); asserted still defined
    ('orpheus/sn/angular/closure.py', 'MorelMontryAngularSweep.angular_endpoint_defect_per_level', "W7 (#653's radial cases un-defer part of the curvilinear build) / curvilinear SN cluster"),  # row 83
    ('orpheus/sn/sweep/pairing.py', 'pair_diffusion_limit_consistent', "W4, concept rows 'Spatial scheme' and 'Angular closure'"),  # row 100
    ('orpheus/sn/sweep/cache.py', 'CollisionCache.reset_build_count', 'deferred, private state (decision D2; section 8.3)'),  # row 99
)
```

#### D-D.7 The rows module `tests/gates/_retired/rows_w0_dead_d.py` (data of the shared gate, 1.5)

Row types are `Retired(file, qualname)`, `Moved(file, qualname, home)` and `RetiredName(name)` of `tests/gates/_retired/__init__.py`; every `Moved` row carries its concrete home and, in its comment, the free function it becomes (decision D2: the object is the first argument and every call site is rewritten); `KEPT` rows are asserted still defined. The sub-unit's rows are listed in the order of its appendix table (`#### D-D.1`).

```python
from tests.gates._retired import Moved, Retired, RetiredName

ROWS = (
    Retired('orpheus/numerics/spaces/scalar_trace_space.py', 'ScalarTraceSpace.slot_shape'),
    Retired('orpheus/transport/boundary_conditions/admission.py', 'Approximated.issues'),
    Retired('orpheus/transport/boundary_conditions/admission.py', 'FaceSite.hollow_inner'),
    Retired('orpheus/transport/boundary_conditions/binding.py', 'Binding.faces_of'),
    Moved('orpheus/transport/fields/_bases.py', 'FaceField.face_views', 'tests/gates/transport/fields/test_angular_boundary_flux.py'),  # row 105; free function face_views(face_field)  # was a property
    Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.at_ordinate'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.l_block'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.isotropic_part'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.anisotropic_part'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.truncate'),
    Retired('orpheus/transport/fields/harmonic_moment_source_sink.py', 'HarmonicMomentSourceSink._CELL_GROUP_NAME'),
    Moved('orpheus/transport/fields/scalar_boundary_flux.py', 'ScalarBoundaryFlux.p1_boundary_scalar_flux', 'tests/gates/transport/fields/_scalar_boundary_oracles.py'),  # row 113; free function p1_boundary_scalar_flux(scalar_boundary_flux, face)
    Retired('orpheus/transport/fields/scalar_flux.py', 'ScalarFlux.at_group'),
    Retired('orpheus/transport/kernels.py', 'FissionKernel.dyad'),
    Moved('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.fission_production_per_material', 'tests/gates/transport/test_kernels.py'),  # row 117; free function fission_production_per_material(material_x_s_field, material_id)
    Moved('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.chi_per_material', 'tests/gates/transport/test_kernels.py'),  # row 118; free function chi_per_material(material_x_s_field, material_id)
    Moved('orpheus/transport/operators/isotropic_transfer.py', 'IsotropicTransfer.dense_per_material', 'tests/gates/transport/_transfer_oracles.py'),  # row 119; free function dense_per_material(isotropic_transfer)
    Moved('orpheus/transport/operators/transfer.py', 'LegendreMomentTransfer.on_basis', 'tests/gates/transport/_transfer_oracles.py'),  # row 121; free function on_basis(legendre_moment_transfer, transfer, basis, *, skip_l0)  # was a classmethod: first argument the class
    Moved('orpheus/transport/operators/transfer.py', 'LegendreMomentTransfer.on_frame', 'tests/gates/transport/_transfer_oracles.py'),  # row 122; free function on_frame(legendre_moment_transfer, transfer, frame, *, skip_l0)  # was a classmethod: first argument the class
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.foldable_part'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.residual_part'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.is_foldable_into_sigma_r'),
    Retired('orpheus/transport/residuals/angular_residual.py', 'AngularResidual.balance_map'),
    Moved('orpheus/transport/spatial/_ubld.py', 'D1ClosedForm.kernel_rhs', 'tests/gates/transport/spatial/test_ld_ubld_primitive.py'),  # row 129; free function kernel_rhs(d1_closed_form, Q_cells, psi_in)
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.foldable_sigma'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator._sibling'),
    Retired('orpheus/numerics/spaces/moment_head.py', 'truncated_head'),
)

KEPT = (  # (file, qualname, kept_for); asserted still defined
    ('orpheus/numerics/spaces/full_field_space.py', 'FullFieldSpace.kernel_basis', 'as row 48'),  # row 78
    ('orpheus/transport/boundary_conditions/realization.py', 'RealizerTable.catalogue', "door: the method tables' admission matrix, asserted through it by 13 test files (5 methods); nothing in W0 opens the realizer-table cycle (3.4)"),  # row 104
    ('orpheus/transport/operators/multiplication_operator.py', 'MultiplicationOperator.from_mesh', 'a documented survivor pinned by a design gate (transport/fields/test_meshless_construction.py:245)'),  # row 120
    ('orpheus/transport/residuals/angular_residual.py', 'AngularResidual.relative_to', "W5, concept row 'Iteration record' (one estimator, one norm)"),  # row 127
    ('orpheus/transport/residuals/scalar_residual.py', 'ScalarResidual', "W3, concept row 'Flux and source members' (the role forest)"),  # row 128
    ('orpheus/transport/fields/_bases.py', 'BoundaryField.from_face_arrays', 'deferred, private state (decision D2; section 8.3)'),  # row 106
    ('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.n2n_matrix', 'deferred, private state (decision D2; section 8.3)'),  # row 116
)
```


## 9 Amendment record (2026-11-20)

One line per finding of the ten check files `scratch/ontology/w0/spec/check_{elegance,structure}_{insulation,dead,roles,defects,posing}.md` (a builder does not read them; this record is for the orchestrator). The ids are the files' own: the numbered findings (`1`..`9`/`10`, `F1`..`F11`, `S1`..`S9`) and the items they list as decisions or a structural fact (`D1`.., `Decision 1`/`2`, `Structural fact`). The decisions cited (G1 to G10, I1 to I7, D1 to D6, R1 to R4, F1 to F4, P1 to P3) are those of the orchestrator's amendment decisions of 2026-11-20, inlined where they act.

- `check_elegance_insulation.md` 1: resolved at 2.4 G1 (decision I3: `L0_PACKAGES`, `ORPHEUS_TOP_LEVEL - L0_PACKAGES`).
- `check_elegance_insulation.md` 2: resolved at 2.4 G1 (decision I3: `_imported_targets` + `_dynamic_imports`, called).
- `check_elegance_insulation.md` 3: resolved at 2.2 item 2 and 8.1 (decision I2: private helper tagged, W8a issue filed).
- `check_elegance_insulation.md` 4: resolved at 2.2 item 2 (decision I2: no `_gl_panel` cache).
- `check_elegance_insulation.md` 5: resolved at 2.2 item 1, 2.3 (a), 2.7 D.18 (partition is the one edge parse).
- `check_elegance_insulation.md` 6: resolved at 2.2 item 1 and 8.1 (follow-up issue for the slot decode).
- `check_elegance_insulation.md` 7: resolved at 2.2 item 1 and 5, 2.4 G4 (decision I4: no hollow `CPMesh` can be built, so no k exists to capture; the function-level leg discharges it, ruled 8.22 (a)).
- `check_elegance_insulation.md` 8: resolved at 2.8 item 5 (decision I6).
- `check_elegance_insulation.md` 9: resolved at 2.4 G6 (decision I6: override-differs leg, docstring and message).
- `check_elegance_insulation.md` 10: resolved at 2.3, 2.4 G5, 7 intro and 7.2 (decision G10), 7.1 row 2.
- `check_structure_insulation.md` F1: decided otherwise by decision I1 (2.0, 2.3 (f), 2.6.2: root fix first).
- `check_structure_insulation.md` F2: resolved at 2.2 item 2, 8.1 (decision I2).
- `check_structure_insulation.md` F3: resolved at 2.2 item 2 (decision I2).
- `check_structure_insulation.md` F4: resolved at 2.3 (b), 2.4 G5 (c), 2.5 (decision I5).
- `check_structure_insulation.md` F5: resolved at 2.4 G6 (decision I6).
- `check_structure_insulation.md` F6: resolved at 2.4 G1 (decision I3: leg (a) dropped).
- `check_structure_insulation.md` F7: resolved at 2.4 G2 (decision I3: `spec_from_file_location`).
- `check_structure_insulation.md` F8: resolved at 2.4 G1 (decision I3).
- `check_structure_insulation.md` F9: resolved at 2.3 (d).
- `check_structure_insulation.md` F10: resolved at 2.2 item 1 and 2.4 G4 (function-level hollow leg; the captured-k half is discharged by it, ruled 8.22 (a)).
- `check_structure_insulation.md` F11: resolved at 7 intro, 7.2 and 2.4 G1 blind spot (decision G10).
- `check_structure_insulation.md` D1: decided otherwise by decision I1.
- `check_structure_insulation.md` D2: resolved at 8.1 (decision I2).
- `check_structure_insulation.md` D3: resolved at 2.5 and 8.1 (decision I7).
- `check_structure_insulation.md` D4: resolved at 8 item 34 (decision I7).
- `check_elegance_dead.md` F1: resolved at 1.5 (decision G3: definition-site and import-resolution binding).
- `check_elegance_dead.md` F2: resolved at 1.5 (decision G8), D-A.7 to D-D.7, 3.3.
- `check_elegance_dead.md` F3: resolved at 3.1 move rows, rule 1 (decision D1).
- `check_elegance_dead.md` F4: resolved at 3.1 move rows, rule 2 (decision D2).
- `check_elegance_dead.md` F5: resolved at 3.2 D4 and 3.7 item 5 (decision G7).
- `check_elegance_dead.md` F6: resolved at 1.5 and 3.1 gate D3 (decision G4).
- `check_elegance_dead.md` F7: resolved at 3.3 row 104.
- `check_elegance_dead.md` F8: decided otherwise by decision G9 (field is an unruled closure member), and the fork ruled 8.23 (b): the four helpers and the field stay (3.0 item 2).
- `check_elegance_dead.md` F9: resolved at 3.2 D4.
- `check_structure_dead.md` F1: resolved at 1.5 (decision G3).
- `check_structure_dead.md` F2: resolved at 3.2 D4 and 3.7 item 5 (decision G7).
- `check_structure_dead.md` F3: resolved at 3.1 move rows, rule 1 (decision D1).
- `check_structure_dead.md` F4: resolved at 3.1 rule 3 and 8.3 (decision G9).
- `check_structure_dead.md` F5: resolved at 3.6 dead-B row 67 (decision D5).
- `check_structure_dead.md` F6: resolved at 3.1 gate D3 (decision G4).
- `check_structure_dead.md` F7: resolved at 1.5 (KEPT mutation) and 3.2 D6.
- `check_structure_dead.md` F8: resolved at 3.4, 3.6 and 8 item 29 (one home for row 126).
- `check_structure_dead.md` Decision 1: resolved at 3.1 rule 3 and 8.3 (decision G9).
- `check_structure_dead.md` Decision 2: resolved at 3.1 rule 1 (decision D1).
- `check_elegance_roles.md` 1: resolved at 1.5 (decision G1).
- `check_elegance_roles.md` 2: resolved at 1.5 (decision G2).
- `check_elegance_roles.md` 3: resolved at 4.1 (decision R1: biproduct comment kept, Arnoldi lesson located in `boundary_conditions.rst`).
- `check_elegance_roles.md` 4: resolved at 1.5 string leg and 4.2 R3 (decision G5).
- `check_elegance_roles.md` 5: resolved at 4.2 R2 paragraph (decision R2: measured, not duplicates, assertions move).
- `check_elegance_roles.md` 6: resolved at 4.1 population, 4.2 R1, 4.5 item 2 (decision G7).
- `check_elegance_roles.md` 7: resolved at 1.5 and 4.2 R1 (decision G6).
- `check_elegance_roles.md` 8: resolved at 4 Outcome (ruled 8.27 (a): the orchestrator writes the ledger entry).
- `check_elegance_roles.md` 9: resolved at 4.2 R2 summary and 7.3.
- `check_structure_roles.md` F1: resolved at 4.5 item 5 (decision G7).
- `check_structure_roles.md` F2: resolved at 4.1, 4.2 R1, 4.5 item 2 (decision R4).
- `check_structure_roles.md` F3: resolved at 4.1 further edits and R-3b (decision R4).
- `check_structure_roles.md` F4: resolved at 4.4 E.22, 4.2 R2 paragraph, R-2c (decision R4); the path in the mutation command at about line 4624 is re-pointed with a dated-history record, ruled 8.26 (a).
- `check_structure_roles.md` F5: resolved at 4.1 further edits (decision R4).
- `check_structure_roles.md` F6: resolved at 4.4 G.24, the 8-level order (decision R4).
- `check_structure_roles.md` F7: resolved at 4.2 R1 rows and 1.5 node kinds.
- `check_structure_roles.md` F8: resolved at 4.2 R2 paragraph (decision R2).
- `check_structure_roles.md` F9: resolved at 4 Outcome (decision R3).
- `check_structure_roles.md` F10: resolved at 1.5 string leg and population (decision G5).
- `check_structure_roles.md` Structural fact: resolved at 4 Outcome (decision R3).
- `check_elegance_defects.md` 1: resolved at 5.2 #709 predicates (decision F1).
- `check_elegance_defects.md` 2: resolved at 5.1 and 5.2 #709 (decision F1).
- `check_elegance_defects.md` 3: resolved at 5.2 #710 (decision F2: `FromQuadrature.measure`, one `match`).
- `check_elegance_defects.md` 4: resolved at 5.3 mutations (decision G6).
- `check_elegance_defects.md` 5: resolved at 5.3 (decision F4: leg (3) dropped).
- `check_elegance_defects.md` 6: resolved at 5.2 #711 (decision F3: record, Chebyshev legs, counts).
- `check_elegance_defects.md` 7: resolved at 5.1 and 5.4 (D.15 docstrings).
- `check_elegance_defects.md` 8: resolved at 5.3 wave leg (decision F4).
- `check_elegance_defects.md` 9: resolved at 5.2 #711 tag naming both stores (decision F4).
- `check_elegance_defects.md` D1: resolved at 5.2 tag sentences (decision F4: existing machinery and call-site set).
- `check_structure_defects.md` S1: resolved at 5.2 #709 law and 5.6 item 4 (decision F4).
- `check_structure_defects.md` S2: resolved at 5.3 (decision F4: leg (3) dropped).
- `check_structure_defects.md` S3: resolved at 5.2 #710 (decision F2: rank-2 admitted leg).
- `check_structure_defects.md` S4: resolved at 5.2 #709 (decision F1).
- `check_structure_defects.md` S5: resolved at 5.2 #711 (decision F3: failing converse stated).
- `check_structure_defects.md` S6: resolved at 5.2 #710 law (decision F2).
- `check_structure_defects.md` S7: resolved at 5.3 wave leg (decision F4: `(?!\d)`, limits stated).
- `check_structure_defects.md` S8: resolved at 5.2 #710 first red.
- `check_structure_defects.md` S9: resolved at 5.3 (decision F4: `cwd`).
- `check_elegance_posing.md` 1: resolved at 6.4 (decision P1: riders, restated done-when).
- `check_elegance_posing.md` 2: resolved at 6.2 (decision P2: `assert_never`).
- `check_elegance_posing.md` 3: resolved at 6.3 (decision P3: tag, producers named).
- `check_elegance_posing.md` 4: resolved at 6.3 and 6.6 item 2 (ruled 8.25 (b): every complex dtype is refused).
- `check_elegance_posing.md` 5: resolved at 6.6 items 4 to 6.
- `check_structure_posing.md` F1: resolved at 6.1 and 6.3 (decision P3).
- `check_structure_posing.md` F2: resolved at 6.2 (decision P2: `match` first, `undetermined=True` leg).
- `check_structure_posing.md` F3: resolved at 6.2 (decision P2: `assert_never`, arm-set leg).
- `check_structure_posing.md` F4: resolved at 6.2 and 6.3 (decision P3: test-only first reds stated, tag).
- `check_structure_posing.md` F5: resolved at 6.4 (decision P1: inlined corpus, 29 lines in 16 files classified).
- `check_structure_posing.md` F6: resolved at 6.3 (decision P3: nonzero reals under multiplication).
- `check_structure_posing.md` F7: resolved at 6.6 item 5 (decision P3).

Counts. Findings in the check files, by the ids above: `elegance_insulation` 10; `structure_insulation` 15; `elegance_dead` 9; `structure_dead` 10; `elegance_roles` 9; `structure_roles` 11; `elegance_defects` 10; `structure_defects` 9; `elegance_posing` 5; `structure_posing` 7; total 95 (87 numbered findings plus 8 decision or structural-fact items: 4 in `check_structure_insulation.md`, 2 in `check_structure_dead.md`, 1 in `check_elegance_defects.md`, 1 in `check_structure_roles.md`). Lines in this record: 95. Handled 95 of 95: 92 resolved in the text (3 of them after the orchestrator's rulings of 2026-11-20, section 8.4), 3 decided otherwise by an orchestrator decision; 0 open.

### Second pass (2026-11-20)

The orchestrator's six rulings of section 8.4, applied (backup of the previous text: `scratch/ontology/w0/spec/amend/spec_backup_pass2.md`; counts by `scratch/ontology/w0/spec/amend/count_rows.py`).

- 8.22 (a): edited 2.2 item 1 (the "Domain, measured" paragraph), 2.2 item 5 (the hollow body), 2.4 G4, section 8.4, and the section 9 lines of `check_elegance_insulation.md` 7 and `check_structure_insulation.md` F10.
- 8.23 (b): edited 1.2 item 3, 1.5 (first red and G8), 3 Outcome, 3.0 item 2, 3.1 (table row for dead-A, the counts sentence), 3.2, 3.3 (heading, four new rows, rows-module assignment, `KEPT` count), 3.6 dead-A (row list, selector closure, gate leg), 3.7 items 1, 3 and 5, section 8 items 5 and 30 and the 8.3 closure list, section 8.4, Appendix D-A.1 (rows 134 to 137 verdict `keep (8.23)`), the Appendix D-A.1 corrections block, D-A.7 (four `Retired` rows out, four `KEPT` rows in), the section 9 line of `check_elegance_dead.md` F8.
- 8.24 (a): edited 3.0 item 3, 3.1 rule 3 (the one exception), gate D3 (the `**Retired:**` paragraph, the witness, the catalogue prose), 3.6 dead-A (the `QuadratureSpec` docstring), 3.7 item 1, section 8 item 7, 8.3, section 8.4, Appendix D-A.1 row 73, the Appendix corrections block, D-A.7 (one added `Retired` entry, not a row).
- 8.25 (b): edited 6.3 (the law, the code block, the widening sentence, the gate legs, the first red), 6.6 item 2, section 8.4, the section 9 line of `check_elegance_posing.md` 4.
- 8.26 (a): edited 4.4 (E.22, the swap-mutation command), section 8.4, the section 9 line of `check_structure_roles.md` F4.
- 8.27 (a): edited 4 Outcome (the ledger entry is the orchestrator's, at close-out; no plan is edited by the unit), section 8.4, the section 9 line of `check_elegance_roles.md` 8.

### Third pass (2026-11-20): #708

The orchestrator's ruling that #708 is a ROOT FIX in W0, and its two answers (the fifth spelling is in; the form is the alternative: delete the arm, `__radd__` and `CONTAINS`, no tag, no `_WAVE_OF[708]`). Backup of the previous text: `scratch/ontology/w0/spec/amend/spec_backup_pass3.md`; probes and drafts under `scratch/ontology/w0/spec/amend/p708/` (`p_rules.py`, `p_perm.py`, `p_lift12.py`, `p_perm_lift.py`, `p_refused_input.py`, `p_radial3.py`, `p_rewire*.py`, `p_time.py`, `fix708.py`, `eqcheck708.py`, `radial708.py`, `gate1.py`) and `scratch/ontology/w0/spec/amend/p708_census.py`, `p708_contains.py`, `p708_overlap2.py`.

- Section 1: 1.1 (the scope names the sub-unit), 1.1b (the sub-unit's brief), 1.2 item 5 (the sentence "#708 is NOT in this specification" is replaced), 1.6 (what this pass did not measure).
- Section 5: the heading, the Outcome and 5.5 (the sentences saying #708 waits), and the NEW sub-section 5.7 with 5.7.1 to 5.7.9 (file set, premise measured, exact edits E1 to E6, law, gates G1 to G3 with their mutations, the 8 tests and the one in `test_n2n_operator.py`, bit-identity and the serial run, docs and docstrings with lines and tense, retirement audit and done-when). Sections 5.1 to 5.6 are unchanged and keep their numbers.
- Section 7: the intro paragraph, 7.1 (new row `7b`; rows 0 to 8 keep their numbers), 7.2 (heading, a table of the 6 code conflicts of the sub-unit and the list of its 8 docs pages), 7.3 (closing #708 and #288, the supersession of B3's root wave, regeneration).
- Section 8: 8.0 items 35 and 36 (the form of the fifth spelling and the `DiscreteMeasure.mass` question; B3's root wave superseded), 8.4 items 8.28 and 8.29 (the form of the fix; the repeal of the endorsement), 8.1 (the half-trace mass issue, the close list).
- Section 9: this block.


### Fourth pass (2026-11-20): quadrature-ties

The orchestrator's five rulings on the exact-tie unit (recorded as the rulings 8.42 to 8.46 of 8.4, with the builder's forks as 8.0 items 37 to 41) and the fill of the placeholder 2.0. Backup of the previous text: `scratch/ontology/w0/spec/amend/spec_backup_pass4.md`; probes, prototypes and the draft gates under `scratch/ontology/w0/spec/amend/ties/` (`census.py`, `overlap.py`, `generic_recurrence.py`, `compare_generic.py`, `time_builders.py`, `rational_polynomial.py`, `final_decision.py`, `fix_generic.py`, `plugin_final.py`, `probes_ab.py`, `probes_ab_final.py`, `probes_c.py`, `probes_final.py`, `probes_m5.py`, `section_2_0.md`, `apply_pass4.py`), and the memo copies re-run from there (`test_exact_tie_gate.py`, `plugin_mut.py`, `p1_population.py`).

- Section 1: 1.1 (the sixth unit), 1.1b (its brief), 1.6 (what this pass did not measure).
- Section 2: 2.0 written whole (2.0.1 to 2.0.12: the premise measured, the file set, the law with its converse and its domain, the generic form with the six builders and the table of what each admits in rationals, the exact edits E1 to E3 with the code, layers, the gates G1 to G7 with their mutations M1 to M9, bit-identity and the one test edited, docs made false, the ERR entry, done-when, unmeasured); the key `quadrature-odd-n` renamed `quadrature-ties` in sections 1, 2, 5, 7 and 8; 2.3 (f) cut to a one-sentence pointer, 2.6.2 to a pointer (the leg of G3 stays, 8.47), 2.4 G3's first red; the sentences that said the odd-`n` leg is red until the precondition lands corrected.
- Section 7: the intro paragraph, 7.1 steps 0 and 2, 7.2 (heading 12 to 14 conflicts; the `generating_measure.py` row rewritten; new rows for `test_generating_measure.py` and `test_concept_order.py`), 7.3 (ERR-127 in the regeneration).
- Section 8: 8.0 items 37 to 41 (the builder's forks), 8.4 rulings 8.42 to 8.46 (the orchestrator's five rulings, with the rejected options' structural reasons), item 2 and item 14 of the numbered list (pointers), 8.1 (one new issue entry: the escalation and the gcd skip), the new 8.5 "Open for the orchestrator" with 8.47; 1.2 item 9.
- Premise mismatches found and handled: (1) the memo's list and the brief name a strict xfail at `n` = 3 that the specification never had (8.47); (2) the memo's prototype converts the node `mpf` to a `Fraction` through `mpmath.mpf(x)`, which re-rounds to the AMBIENT precision outside a `workdps` block and makes the 1e-20 bracket fail to isolate the node (measured by this pass); production calls it inside the working precision, and the text specifies the exact conversion `x.man_exp` and a gate row that calls the decision outside any precision context; (3) the memo's 12 rows leave the weight decision's "gcd non-trivial, node not a root" arm without a witness (measured: the mutant `return len(g) > 1` keeps 12 of 12 green) and the anti-guesser row stops on a NODE (the weight decision is invoked 0 times): G4 is the added witness; (4) the memo's prototype raises an unlisted refusal (`does not isolate the node`): given a fragment, a Raises line and a witness; (5) the memo's gate is 8 Legendre rows only, so the generic claim had no first red: G5 adds 6 rows of the other families (6 of 6 raise today).
- Re-run by this pass and reproduced: the 12-row gate (today 10 failed, 2 passed; with the fix 12 passed; mutations 10, 7, 3, 1 failed) and the population probe (384 of 1088 odd-`n` pairs, 0 of 1088 even, 0 mismatches of the predicate over 2176 pairs, smallest reproducer `(0.7937005259840998, 0.8735804647362988)` at `n` = 1). Not reproduced: nothing differed; the figures of the memo that were not re-run are listed in 2.0.12.
- Searched for and not found: the memo's "`CPParams.n_quad_y` comment" (`git grep -n -i -E "odd|parity|rounding boundary|tie" -- orpheus/cp/solver.py` prints 3 lines, none about `n_quad_y`; the specification's lines naming `n_quad_y` (2.2, 2.3, 2.4 G4) carry no such comment), and the memo's F1 and D1 (they name findings of a check file that is not part of this specification; the only `F1`/`D1` in the specification are the retirement engine's decision labels).
