# CS3-R completeness census — torsor-era survivors after the cone carve

**Author:** `qa` (independent second filter). **Date:** 2026-08-19.
**Tree state:** `HEAD = f43758d8` (branch `main`); CS3 steps 1–3 =
`c3e66b18` / `993fa280` / `5efd2178`, all ancestors of HEAD.
**Scope:** `orpheus/`, `tests/`, `docs/`, `.claude/{skills,agents,rules}/`,
plus a LIGHT present-tense-only pass over `scratch/cs3_verification_plan.md`
and `scratch/omr_v2_grounding/*.md`.

> ⛔ **INSTRUMENT FAILURE FOUND AND FIXED BEFORE ANY RESULT WAS TRUSTED —
> read this first.** My first two sweeps reported **all-clean** across every
> tree. Both were VOID. Cause: `grep -rn PAT $TREES` where
> `TREES="orpheus tests docs …"`; **zsh does not word-split an unquoted
> `$VAR`**, so the whole string went in as ONE nonexistent path, and the
> `2>/dev/null` I had appended swallowed the `No such file or directory`.
> The failure is silent AND flattering — a dropped tree looks exactly like a
> clean tree. Caught by a positive control (`grep -rl flux` per tree →
> 195/281/72/12/8/2 files), which is now run before every sweep.
> A second hazard, ruled out separately: the interactive `grep` is a shell
> FUNCTION wrapping `ugrep --ignore-files`, i.e. it honours `.gitignore` and
> would hide untracked/ignored files. Every sweep below therefore uses
> **`command grep -rnE`** (real BSD grep, `.gitignore`-blind), cross-checked
> against `git grep`.

---

## (a) Findings table

Verdict key — **MUST-FIX** = present-tense-false / imperative / live Sphinx role
at a deleted target. **HISTORY-OK** = dated past tense or ⛔-marked, keeps the
name as inert literal text. **AMBIGUOUS** = judgement call, reason given.

### A. `orpheus/` — production

| file:line | spelling | fragment | verdict |
|---|---|---|---|
| `orpheus/transport/fields/_bases.py:1160` | `RadialCharacteristicInteriorDisplacement` | *"The concrete role leaves **are** ``RadialCharacteristicInteriorFlux`` (state) **and ``RadialCharacteristicInteriorDisplacement`` (the iterate increment)**. Abstract — instantiate a role leaf."* | **MUST-FIX** — class docstring names a DELETED class as a live role leaf, 26 lines below this file's own corrected note at `:1134` |
| `orpheus/transport/fields/_bases.py:1220` | `RadialCharacteristicBoundaryDisplacement` | *"The concrete role leaves **are** ``RadialCharacteristicBoundaryFlux`` (state) and ``RadialCharacteristicBoundaryDisplacement`` (the iterate increment)."* | **MUST-FIX** — twin of the above |
| `orpheus/transport/radial_characteristic_field.py:28` | `⊖`, displacement | *"ONE composite class carries the flux state …, the source emission …, and **the displacement (minted per block by ``⊖``)**. Role identity lives on the MEMBERS"* | **MUST-FIX** — module docstring, present tense; `⊖` is a retired operator and the displacement role is gone |
| `orpheus/transport/radial_characteristic_field.py:101` | `⊖`, displacement | *"(the marched cells — flux state, source emission, or, **after ``composite ⊖ composite``, an interior displacement**)"* | **MUST-FIX** — class docstring; **self-contradicting**: 5 lines lower the SAME docstring says *"plain V arithmetic per block since campaign 1 CS3"* |
| `orpheus/transport/radial_characteristic_field.py:114` | `⊖`, displacement | *"# flux, source (operator emissions), and **displacement (from ⊖)** composites are all admitted"* | **MUST-FIX** — inline comment stating the guard's rationale in retired terms |
| `orpheus/transport/radial_characteristic_field.py:121`, `:125` | displacement (runtime string) | `f"…must be a RadialCharacteristicInteriorField (its **flux / source / displacement** leaf); got …"` | **MUST-FIX (low sev.)** — user-facing `TypeError` text teaching a 3-role taxonomy that is now 2. ⚠ Retirement rule: grep the shortest distinctive fragment before rewording — see §(c) check T4, **no test pins it** |
| `orpheus/transport/__init__.py:31` | `orpheus.transport.displacements` | *"(``orpheus.transport.displacements`` … — **RETIRED at campaign 1 CS3, 2026-08-19**: flux lives in V …)"* | **HISTORY-OK** — dated, and correctly demoted from `:mod:` to a double-backtick literal |
| `orpheus/transport/fields/_bases.py:18`, `:1134` · `residuals/__init__.py:25` · `full_field.py:111` · `numerics/vector.py:59` · `numerics/coupled_system.py:110,242` · `sn/solver.py:2961` · `fields/_coefficient_role.py:9-83` · `fields/cross_section_field.py:26,97` · `fields/harmonic_moment_flux.py:26` · `fields/radial_characteristic_{interior,boundary}_flux.py:17/18` · `source_sinks/harmonic_moment_source_sink.py:36,126` | `FluxRole`, torsor, Displacement | *"(A fourth role, Displacement, **existed until** campaign 1 CS3 …)"*, *"the ``FluxRole`` affine gate **is retired**"*, *"**historically** the complement of the flux torsor"* | **HISTORY-OK** ×15 — every one dated and past-tensed |
| `orpheus/derivations/discrete/sn/dsa.py:963` · `orpheus/sn/acceleration/dsa.py:404,407` | `displacements` | *"displacements :math:`(d_0, d_1)` to the row sources"*, *"The raw scalar-moment sweep displacements"* | **HISTORY-OK / not-a-hit** — DSA's own physics vocabulary for the sweep delta `(d_0,d_1)`; unrelated to the retired type. See §(b) triage note |
| `orpheus/numerics/measure.py:1276-1333` · `orpheus/data/**` | base point | fiber-bundle base point; sigma-zero base points | **not-a-hit** — two unrelated physics/maths senses |

### B. `tests/`

| file:line | spelling | fragment | verdict |
|---|---|---|---|
| `tests/sn/operators/test_operators_apply_typed.py:311`, `:318` | affine flux STATES, displacements | docstring: *"``A.apply(α·d₁+β·d₂) == …`` **on flux displacements** … since ``α·ψ + β·ψ'`` with ``α+β ≠ 1`` **is illegal** on affine flux STATES (no origin) **but legal on displacements**."* | **MUST-FIX** |
| `tests/sn/operators/test_operators_apply_typed.py:332-336` | `⊖`, torsor form | *"# #208: a general α·ψ₁+β·ψ₂ (α+β≠1) **is illegal** on affine flux STATES; verify linearity via the affine-supported ops … **affine additivity in torsor form** A(ψ₁ + λ(ψ₂⊖ψ₁))"* | **MUST-FIX** + torsor-SHAPED (see §(c) T1) |
| `tests/sn/operators/test_collision_operator.py:208-212` | same | *"# #208: a general α·ψ₁+β·ψ₂ (α+β≠1) **is illegal** on affine flux STATES (no origin); … **affine additivity in torsor form** ψ₁ + λ(ψ₂⊖ψ₁)"* | **MUST-FIX** + torsor-SHAPED |
| `tests/sn/operators/test_scattering_operator.py:460-467` | same | docstring: *"A general ``α·ψ₁ + β·ψ₂`` with ``α+β ≠ 1`` **is illegal** on affine flux STATES (no origin), so linearity is verified by … **affine additivity in torsor form**"* | **MUST-FIX** + torsor-SHAPED |
| `tests/sn/operators/test_streaming_operator.py:300-304` | same | *"# #208: a general α·ψ₁+β·ψ₂ (α+β≠1) **is illegal** on affine flux STATES; … **torsor form** ψ₁ + λ(ψ₂⊖ψ₁)"* | **MUST-FIX** + torsor-SHAPED |
| `tests/sn/sweep/core/test_phase_c_gates.py:258-263` | same | docstring: *"#208 affine reframe — a general ``α·ψ + β·φ`` with ``α+β ≠ 1`` **is illegal** on affine flux STATES (no origin), so linearity is verified by … **affine additivity in torsor form**"* | **MUST-FIX** + torsor-SHAPED |
| `tests/sn/mesh/test_radial_characteristic_split_leaves.py:220-224` | `⊖`, affine points, displacement | class docstring: *"Unlike the flux leaves (**affine points, ``⊖`` mints displacements**), source sums are CLOSED … The Field class-identity gate is what keeps source, flux, and **displacement** arithmetic from mixing"* | **MUST-FIX** — 205 lines below this file's OWN corrected module docstring (`:14-17`) |
| `tests/sn/primitives/test_typed_source_sinks.py:376-378` | `FluxRole`, displacement | *"in deliberate contrast to its **FluxRole sibling** HarmonicMomentFlux, whose ``flux + flux`` **is forbidden** (affine, no origin) and ``flux − flux`` **mints a displacement**"* | **MUST-FIX** — this is the stated RATIONALE for the whole test block; production `harmonic_moment_flux.py:26` says *"it **was** a ``(FluxRole, MomentField)`` carrier"* |
| `tests/transport/fields/test_angular_boundary_flux.py:162` | `AngularBoundaryDisplacement` | `d1 = bf1 - base  # AngularBoundaryDisplacement (V is a vector space)` | **MUST-FIX** — the parenthetical is the NEW truth, the type name beside it is the OLD one; `bf1 - base` returns `AngularBoundaryFlux` |
| `tests/transport/fields/test_angular_boundary_flux.py:346` | `AngularBoundaryDisplacement` | `out = bf1 - bf2  # AngularBoundaryDisplacement — shares the trace layout` | **MUST-FIX** — same |
| `tests/sn/mesh/test_radial_characteristic_split_leaves.py:173` | torsor | `pytest.fail(f"**torsor action** returned {type(recovered).__name__}")` — 2 lines under `# plain V arithmetic` | **MUST-FIX (low sev.)** — vestigial name in a live failure message |
| `tests/transport/test_radial_characteristic_field.py:188` | torsor | `pytest.fail(f"**torsor** returned {type(recovered).__name__}")` — 2 lines under `# plain V arithmetic, blockwise` | **MUST-FIX (low sev.)** — same |
| `tests/transport/fields/test_scalar_boundary_flux.py:199` · `test_angular_boundary_flux.py:102` · `test_radial_characteristic_field.py:171` · `test_radial_characteristic_split_leaves.py:15-17,138` · `test_flux_vector_algebra.py:20-21,177,229,247` · `test_si_diagnostic_trajectory.py:5-15,156` · `test_coefficient_fields.py:13,32,81,98,124-136` · `test_declared_law_is_linear.py:30,46-55,554` · `test_harmonic_moment_flux.py:377` · `test_si_convergence_diagnostics.py:10` | all | *"the #208 affine discipline **minted**"*, *"**Until** CS3 this minted"*, *"**retired** with its type at step 3"*, *"**unlike** FluxRole"* | **HISTORY-OK** ×~28 — dated / past tense throughout |
| `tests/sn/solve/test_affine_carve_bit_identity.py:3,78,82` | `FluxDisplacement`, `MomentDisplacement`, `AngularDisplacement` | *"The affine typing (Piece 1: ``FluxDisplacement`` + the ``flux+flux`` gate), …"* | **AMBIGUOUS** — see §(b); module is a *bit-identity record of the #208 carve*, so the names are the subject. Needs a one-line dated banner to be unambiguously HISTORY-OK |
| `tests/sn/sweep/core/test_ordinate_scan.py:379` | `affine_combination` | `def test_ordinate_scan_affine_combination` | **not-a-hit** — the scan recurrence `a·x+b`, a different "affine" |

### C. `docs/`

| file:line | spelling | fragment | verdict |
|---|---|---|---|
| `docs/theory/foundations/field_algebra.rst` (the whole §"The overturned design", ~`:1200-1530`, plus `:47,49,82-83,132,208,240,331,396,413-415,454,469,484,759,849-863,1032,1210-1325`) | every retired spelling | *"Flux states **were** typed as an **affine space**"*, *"A fourth role, Displacement, **existed from** 2026-06-08 **until** …"* | **HISTORY-OK** — this page IS the adjudication record; uniformly past-tensed |
| `docs/theory/foundations/operator_algebra.rst:3429,3557,3859-3921,4079,4735,5658,5874-5886` | `Displacement`, torsor, affine space | *"⛔ **A fourth column, Displacement, stood here from 2026-06-08 to 2026-08-19**"*, *"⛔ until 2026-08-19 this paragraph read …"* | **HISTORY-OK** — every site ⛔-marked with dates |
| `docs/theory/methods/sn/acceleration.rst:894-899` | torsor action, displacement | *"⛔ **Until 2026-08-19 this read** *'the update is the torsor action …'*; the affine ontology was overturned at campaign-1 CS3"* | **HISTORY-OK** |
| `docs/theory/methods/sn/curvilinear_one_group.rst:4166-4173` | `RadialCharacteristic*Displacement`, `⊖` | *"⛔ **displacement** — ``RadialCharacteristicInteriorDisplacement`` ⊕ … / the affine displacement between two ψ½ states (minted per block by ⊖) / the composite torsor algebra (2.5d d1) — **retired with the row**"* | **HISTORY-OK** — the ⛔ opens the row and *"retired with the row"* closes it |
| `docs/theory/methods/sn/history.rst:1187` · `docs/theory/foundations/boundary_conditions.rst:3957` · `docs/theory/verification/sn.rst:2733` | `FluxDisplacement`, `AngularDisplacement`, `AngularBoundaryDisplacement` | *"**Wave-O** affine-typed operator algebra"*, *"``AngularBoundaryDisplacement`` **existed** too"* | **HISTORY-OK** ×3 — changelog / past tense |
| `docs/theory/fuel_behaviour.rst:303` | `Displacement-Based Constraint` | closed-gap BC4 section heading | **not-a-hit** — mechanical displacement, a different physics |

### D. `.claude/skills`, `.claude/agents`, `.claude/rules`

| file:line | spelling | fragment | verdict |
|---|---|---|---|
| `.claude/agents/explorer/AGENT.md:151-155` | `FluxDisplacement`, Displacement | *"The role grid **is** {Flux, Source/Sink, Residual, **Displacement**} × {Angular, Scalar, Boundary} — … the SI iterate-delta **is** a ``FluxDisplacement`` (affine difference space; ``flux+flux`` **is** a TypeError, ``flux−flux→displacement`` **is** legal)."* | **MUST-FIX — highest severity of the whole census.** Present tense, triple-false, in the standing brief of the project's designated exploration delegate. AGENT.md loads FRESH per dispatch, so every future `explorer` reports the retired taxonomy as current |
| `.claude/agents/cross-domain-attacker/AGENT.md:297-301` | torsor displacement type | *"3. **An iterative increment … typed as the STATE type x**, not a distinct displacement type. TELL: admits illegal ``state+state`` … **FIX: a difference-space / torsor displacement type.**"* | **MUST-FIX** — an **imperative instruction to re-mint the overturned design**. Its own source skill (`cross-domain-frames/reference.md:189-198`, Shape 3) already carries the ⛔ re-pose; this is the un-migrated twin |
| `.claude/agents/elegance-enforcer/AGENT.md:185` | displacement mixin | *"**6. The residual/displacement ASYMMETRY is principled, not duplication.** … **displacements get a shared mixin** because their distinguishing behaviour is their METHODS (contraction diagnostics). … **NO twin path.**"* | **MUST-FIX** — a standing instruction to the elegance GATEKEEPER not to flag a structure that no longer exists (mixin retired; the diagnostics moved to `IterationRecord`) |
| `.claude/skills/coding-elegance/SKILL.md:390` (the ⚠ clause) | meta-claim | *"⚠ … ``cross-domain-frames/reference.md`` (the A.1 row, and the §192/§201 fix-suggestions) **still cites** ``FluxDisplacement`` … and ``numerical-bug-signatures/SKILL.md`` §479/§488 **still credits** the retired type … Both are **stale as of 2026-08-19** and are owed the same dated correction"* | **MUST-FIX (derivative staleness)** — `[M]` BOTH named targets have SINCE been corrected: `reference.md:25` and `:189-198` carry dated ⛔ blocks, and `numerical-bug-signatures:489-491` reads *"born on #208's **retired** ``FluxDisplacement``"*. The correction landed; its own TODO note did not retire, so it now tells readers to distrust two files that are correct |
| `.claude/skills/cross-domain-frames/reference.md:25` (A.1 row), `:189-198` (Shape 3) | torsor, `FluxDisplacement` | *"⛔ The frame is sound; its ORPHEUS worked example **was OVERTURNED** (campaign 1 CS3, 2026-08-19)"*, *"⛔ **Re-posed at** campaign 1 CS3"* | **HISTORY-OK** — correctly and fully migrated |
| `.claude/skills/numerical-bug-signatures/SKILL.md:489-491` | `FluxDisplacement` | *"``IterationRecord.true_error_estimate()`` … derived on the record since campaign 1 CS3; born on #208's **retired** ``FluxDisplacement``"* | **HISTORY-OK** |
| `.claude/skills/coding-elegance/SKILL.md:258`, `:388-390` (body) | `FluxDisplacement`, torsor | *"⛔ **The worked example here was REFUTED on 2026-08-19**"*, *"⛔ **This item read the opposite until 2026-08-19**"* | **HISTORY-OK** |
| `.claude/rules/coding-standards.md:43` | torsor | *"a torsor ``A×V→A`` that must FORBID ``A×A``, versus a vector ``V×V→V``"* | **AMBIGUOUS → keep.** A *generic* encoding-lattice statement (phantom-type-parameter rule); names no ORPHEUS flux type and asserts nothing about flux. Correct as written. Worth ONE clause noting the ORPHEUS instance was overturned, but not a falsehood |

---

## (a-bis) ⚠ RE-BASELINE — the tree moved MID-CENSUS

⚠ **The table above was measured at `HEAD = f43758d8`. It is kept per
`coding-standards` (past-tense history stays), but it is NOT the current
state.** At 18:27:14, while I was reading contexts, the main agent landed
`755f99b5` *"refactor(transport): CS3-R sweep 1 — torsor vocabulary the carve's
own sweeps missed"* (16 files), and at 18:28:36 five further files were dirty in
the shared working tree. Caught by a `sed` reading *"increments"* where my own
grep 4 minutes earlier had read *"displacements"* at the same `file:line`
(`orpheus/sn/acceleration/dsa.py:404`) — i.e. by a discrepancy, not by a check
I had planned. **`[M]` re-measured at `HEAD = 755f99b5`, working tree at
18:28:51.** Every row below was re-run at that instant.

### REMEDIATED concurrently — no longer findings

`orpheus/transport/radial_characteristic_field.py` (all 5 sites incl. the two
runtime strings) · `orpheus/{sn/acceleration,derivations/discrete/sn}/dsa.py`
(*displacements* → *increments*) · `orpheus/numerics/{coupled_system,iteration,
moment_layout}.py` · `orpheus/sn/{solver,operators/boundary,
operators/radial_characteristic}.py` · `orpheus/transport/spatial/scheme.py` ·
`orpheus/transport/source_sinks/radial_characteristic_*.py` ·
`tests/numerics/test_si_diagnostic_trajectory.py` (the *"**Today** the SI loop
mints…"* paragraph → *"**At capture time** (``000cf144``, pre-carve) the SI loop
**minted**…"*, and the two live `:meth:`…Displacement.X`` **roles** demoted to
double-backtick literals — both of my findings there) ·
`tests/sn/mesh/test_radial_characteristic_split_leaves.py:220` ·
`tests/sn/solve/test_affine_carve_bit_identity.py` · and 4 of the 5
*"is illegal on affine flux STATES"* sites (`test_collision_operator.py`,
`test_streaming_operator.py`, `test_scattering_operator.py`,
`test_phase_c_gates.py` — the last two re-spelled to the correct past tense
*"the pre-CS3 blend spelling … could **not**"*).

### ⭐ SURVIVORS — the 10 MUST-FIX rows still standing at `755f99b5` / 18:28:51

| # | file:line | spelling | fragment | why it survives |
|---|---|---|---|---|
| **S1** | `orpheus/transport/fields/_bases.py:1160` | `RadialCharacteristicInteriorDisplacement` | *"The concrete role leaves **are** ``RadialCharacteristicInteriorFlux`` (state) and **``RadialCharacteristicInteriorDisplacement``** (the iterate increment)."* | **PRODUCTION** docstring naming a DELETED class as a live role leaf. The file's own `:18` and `:1134` record the retirement — a **self-contradicting file** (vv #21 aggravator: it can now be cited for either) |
| **S2** | `orpheus/transport/fields/_bases.py:1220` | `RadialCharacteristicBoundaryDisplacement` | twin of S1 | same |
| **S3** | `tests/sn/operators/test_operators_apply_typed.py:311`, `:318` | affine flux STATES, displacements | docstring: *"… **on flux displacements**. … since ``α·ψ + β·ψ'`` with ``α+β ≠ 1`` **is illegal** on affine flux STATES (no origin) **but legal on displacements**."* | ⭐ **the sharpest row in the census.** The body comment 20 lines BELOW was corrected in the concurrent sweep (*"Linearity, stated directly (campaign 1 CS3 — flux lives in V)"*) while the **docstring above it still asserts the opposite** — and the docstring is what a reader, a `--collect-only -v` listing, and a summariser see FIRST |
| **S4** | `tests/sn/primitives/test_typed_source_sinks.py:376-378` | `FluxRole`, displacement | *"in deliberate contrast to its **FluxRole sibling** HarmonicMomentFlux, whose ``flux + flux`` **is forbidden** (affine, no origin) and ``flux − flux`` **mints a displacement**"* | This is the stated RATIONALE for the whole `HarmonicMomentSourceSink` test block — the "deliberate contrast" it rests on no longer exists (production `harmonic_moment_flux.py:26`: *"it **was** a ``(FluxRole, MomentField)`` carrier"*). The block's premise, not a stray comment |
| **S5** | `tests/transport/fields/test_angular_boundary_flux.py:162` | `AngularBoundaryDisplacement` | `d1 = bf1 - base  # AngularBoundaryDisplacement (V is a vector space)` | The parenthetical is the NEW truth; the type name beside it is the OLD one — one comment, both ontologies. `bf1 - base` returns `AngularBoundaryFlux` |
| **S6** | `tests/transport/fields/test_angular_boundary_flux.py:346` | `AngularBoundaryDisplacement` | `out = bf1 - bf2  # AngularBoundaryDisplacement — shares the trace layout` | same |
| **S7** | `tests/transport/test_radial_characteristic_field.py:144-146` | displacement, dead test name | tombstone *"Successors carrying the surviving claims: … role identity (**flux / source / displacement** members) → … ``TestCompositeAlgebra::**test_subtraction_mints_a_displacement_composite_per_block**``"* | **A DEAD REFERENCE INSIDE A RETIREMENT AUDIT'S OWN SUCCESSOR LIST.** `[M]` that test name exists nowhere in `orpheus/ tests/ docs/ .claude/` except this comment (and a stale `docs/_build` graph artifact) — CS3 renamed it `test_subtraction_is_same_typed_per_block`. The artefact whose job is keeping coverage traceable points at nothing |
| **S8** | `tests/transport/test_radial_characteristic_field.py:188` | torsor | `pytest.fail(f"**torsor** returned {type(recovered).__name__}")` — 2 lines under `# plain V arithmetic, blockwise` | vestigial name in a LIVE failure message (low sev.). Sibling at `tests/sn/mesh/test_radial_characteristic_split_leaves.py:173` (*"torsor action returned"*) also survives |
| **S9** | `.claude/agents/explorer/AGENT.md:151-155` | `FluxDisplacement`, Displacement | *"The role grid **is** {Flux, Source/Sink, Residual, **Displacement**} × {Angular, Scalar, Boundary} — … the SI iterate-delta **is** a ``FluxDisplacement`` (affine difference space; ``flux+flux`` **is** a TypeError, ``flux−flux→displacement`` **is** legal)."* | ⭐⭐ **HIGHEST SEVERITY OF THE CENSUS.** Present tense, **triple**-false (grid is 3 roles; the type is deleted; `flux+flux` is now the headline LEGAL op). It sits in the standing brief of the project's designated exploration delegate (`.claude/rules/nexus-tools.md`), and **AGENT.md loads FRESH per dispatch** — so every future `explorer` will report the overturned taxonomy to a main agent as current fact |
| **S10** | `.claude/agents/cross-domain-attacker/AGENT.md:297-301` | torsor displacement type | *"3. An iterative increment … typed as the STATE type x, not a distinct displacement type. TELL: admits illegal ``state+state`` … **FIX: a difference-space / torsor displacement type.**"* | **An imperative instruction to RE-MINT the overturned design.** Its own source skill already carries the ⛔ re-pose (`cross-domain-frames/reference.md:189-198`, Shape 3: *"⛔ Re-posed at campaign 1 CS3 … the home is the ITERATION RECORD, never a state-space type"*) — this is the **un-migrated twin** of a corrected single source (Pattern-2) |
| **S11** | `.claude/agents/elegance-enforcer/AGENT.md:185` | displacement mixin | *"**6. The residual/displacement ASYMMETRY is principled, not duplication.** … **displacements get a shared mixin** because their distinguishing behaviour is their METHODS (contraction diagnostics). … **NO twin path.**"* | A standing instruction to the elegance GATEKEEPER **not to flag** a structure that no longer exists — the mixin is deleted and the contraction diagnostics moved to `IterationRecord`. A ruling that now defends nothing |
| **S12** | `.claude/skills/coding-elegance/SKILL.md:390` (the ⚠ clause) | meta-claim | *"⚠ … ``cross-domain-frames/reference.md`` … **still cites** ``FluxDisplacement`` … and ``numerical-bug-signatures/SKILL.md`` §479/§488 **still credits** the retired type … **Both are stale as of 2026-08-19**"* | **DERIVATIVE staleness — the correction landed, its own TODO did not retire.** `[M]` both named targets are now CORRECT: `reference.md:25` and `:189-198` carry dated ⛔ blocks; `numerical-bug-signatures:489-491` reads *"born on #208's **retired** ``FluxDisplacement``"*. The clause now tells readers to distrust two files that are right |

*(S1–S12 = **12 rows**, of which S8's sibling makes 13 sites; the headline count
below counts the 12 rows.)*

### AMBIGUOUS — argued, not fixed

| file:line | why ambiguous | recommendation |
|---|---|---|
| `.claude/rules/coding-standards.md:43` | *"a torsor ``A×V→A`` that must FORBID ``A×A``, versus a vector ``V×V→V``"* — a **generic** encoding-lattice statement (the phantom-type-parameter rule). It names no ORPHEUS type and asserts nothing about flux; as written it is **correct** | KEEP. Optionally one clause noting the ORPHEUS flux instance was overturned, so a reader does not reach for the torsor arm by analogy |
| `docs/theory/foundations/field_algebra.rst:1093` label `affine-typed-residual-eq` | The plan (`scratch/cs3_verification_plan.md:314`) listed it among *"four labelled equations [that] are APIs"* to retire. `[M]` the other three (`affine-torsor-algebra`, `affine-contraction-ratio`, `affine-true-error`) are **gone with 0 citers**; this one SURVIVES — and the page carries an explicit ruling 10 lines above it: *"The name is stale; the claim is not"*, with the dangling-`:eq:` risk argued | KEEP — a documented deliberate decision, not an oversight. The plan's four-way grouping was over-broad |

---

## (b) Denominator block — what was actually searched

**Instrument.** `command grep -rnE --binary-files=without-match
--exclude-dir=__pycache__ --exclude-dir=_build` (real BSD grep — `.gitignore`-BLIND,
so untracked files are in scope), over an explicit **zsh array** `TREES=(…)`
expanded as `"${TREES[@]}"`. Harness: `/tmp/cs3r_sweep.sh`. **stderr never
suppressed.** Cross-checked with `git grep` (tracked-only).

**Positive control, run before every sweep** (`grep -rl flux` per tree):

| tree | tracked files | files on disk | control: files containing `flux` |
|---|---:|---:|---:|
| `orpheus/` | 349 | 362 | 195 |
| `tests/` | 607 | 607 | 301 |
| `docs/` | 101 | 111 | 72 |
| `.claude/skills/` | 46 | 46 | 12 |
| `.claude/agents/` | 9 | 9 | 8 |
| `.claude/rules/` | 6 | 6 | 2 |

*(on-disk > tracked only in `orpheus/data/micro_xs/*.h5` (12 binaries, excluded
by `--binary-files`), `docs/_generated/*.rst` (9, Sphinx build products) and two
`.DS_Store`. `[M]` `git status --untracked-files=all` over the six trees returns
**no untracked source file**, so the `.gitignore`-blind and `git grep` instruments
see the same corpus — confirmed numerically: on the union pattern
`[A-Za-z]*Displacement|FluxRole|⊖`, `git grep` = **75** and `command grep` =
**75**, exact.)*

**Every pattern actually run** (ERE, all six trees, all at both baselines):

```
transport\.displacements            transport/displacements            displacements
_flux_role                          fields\._flux_role
FluxRole                            [A-Za-z]*Displacement              \bDisplacement\b
FluxDisplacement                    AngularDisplacement                AngularBoundaryDisplacement
MomentDisplacement                  ScalarDisplacement                 ScalarBoundaryDisplacement
RadialCharacteristicInteriorDisplacement   RadialCharacteristicBoundaryDisplacement
[Tt]orsor                           ⊖|\\ominus                         ⊕ (via the ⊖ sweep's context)
affine_combination                  _check_torsor_partner              _mint_displacement
_DISPLACEMENT_CLS                   _flux_displacement_leaf            sibling_of
last_displacement                   _carrier_rep
contraction_ratio[^s]|contraction_ratio$        where_largest          true_error_estimate
[Bb]ase[ _-]?point|basepoint|BASE[ _-]?POINT
no (natural )?origin|natural origin
affine space|affine point|affine subspace|lives in .?A.?[^a-z]|points of an affine
\bdisp\b|\bdisp_|_disp\b|disp_leaf|prev_disp
:(class|meth|mod|attr|func|obj|data|exc):`[^`]*([Dd]isplacement|FluxRole|affine_combination|_flux_role|sibling_of|_check_torsor_partner)
Σλ|sum.?λ|lambda.*== *1\.0|partition of unity|partition-of-unity
mint[a-z]* a[n]? [A-Za-z]*[Dd]isplacement|mint[a-z]* the displacement|displacement mint
affine-torsor-algebra|affine-contraction-ratio|affine-true-error|affine-typed-residual-eq
concrete role leaves|role leaves are
```

⚠ **Filter validated, not merely written down** (the 2026-08-19 plan-authoring
row: a pattern filter fails silently and *flatteringly*). Each of the three
symbol SPELLINGS was proved against a known member before use: the ASCII
identifier (`_flux_displacement_leaf` → hit in
`test_si_diagnostic_trajectory.py`), the Unicode prose form (`⊖` → 6 hits, so
BSD grep's UTF-8 handling is confirmed live, not assumed), and the LaTeX macro
(`\ominus` → hit in `field_algebra.rst:1215`). The retired **NUMBERS**/labels
were grepped separately (P9), because a stale figure outlives every spelling of
its symbol.

**Trees I could NOT cover / deliberately did not:**

- `.claude/plans/` and the rest of `scratch/` — **excluded by brief.**
- ⚠ **`.claude/agent-memory/` — OUT of brief, and it is the largest uncovered
  surface with the SAME hazard class as findings S9–S11.** `[M]` `git grep` on
  the union pattern over all of `.claude/` returns **257** lines vs **75** in the
  three in-scope subtrees; the 182-line remainder is ~20 durable agent-memory
  files, including `explorer/flux_torsor_vs_cone_inventory.md`,
  `elegance-enforcer/issue_208_affine_flux_algebra.md`,
  `cross-domain-attacker/issue_208_delta_psi_affine_frames.md` and
  `numerics-investigator/issue_208_flux_displacement_residual_typing_debug_value.md`.
  Agent memory is read at the START of a dispatch, exactly like AGENT.md — so if
  S9–S11 are worth fixing, this surface is owed the same pass. **I did not audit
  it; I am naming it so its absence is not read as a clean bill.**
- `docs/_build/**` — excluded as build products, but noted: the stale
  `docs/_build/html_p3base/_nexus/graph.json` still carries deleted node ids
  (this is how S7's dead test name was confirmed to have once existed). A Sphinx
  rebuild clears it; nothing to fix.
- `scratch/cs3_verification_plan.md` + `scratch/omr_v2_grounding/*.md` —
  **light present-tense-only pass, as briefed.** They read as pre-carve plans
  with explicit `today | after` columns and are records in the
  `plan-authoring` §3 sense. **No MUST-FIX.** One structural note: the plan's
  §3.6 claim *"four labelled equations are APIs"* is now **3-of-4 retired,
  1 deliberately kept** — see the AMBIGUOUS table.

---

## (c) Torsor-SHAPED patterns — checks run, INCLUDING the clean ones

A clean check is first-class output. Six patterns were looked for; three were
clean, two were real, one is a documented deliberate keep.

### T1 ⭐ REAL, and it was a genuine GATE DOWNGRADE — not just stale prose

**Shape:** an operator-linearity gate that verifies additivity through the
*affine detour* `op(ψ₁ + λ(ψ₂−ψ₁)) = (1−λ)op(ψ₁) + λop(ψ₂)` instead of directly,
because the retired type forbade `ψ+ψ`. Found at **5 sites**
(`test_collision_operator`, `test_streaming_operator`, `test_scattering_operator`,
`test_phase_c_gates`, `test_operators_apply_typed`).

⭐ **This is a `vv-principles` Mode-12 invariance-functional blindness, and I
measured it rather than quoting it.** Affine maps *preserve affine
combinations*, so an affine regression `A(x) = Lx + q` is annihilated by the
detour form exactly — no tolerance, refinement, or fixture change could ever
expose it. `[M]` pure-numpy probe (`n = 6`, `L` random, `λ = 0.7`), no ORPHEUS
import, no file touched:

| form | `q = 0` (honest linear) | `q ≠ 0` (**affine bug**) |
|---|---:|---:|
| RETIRED detour `A(ψ₁+λ(ψ₂−ψ₁))` | `4.440892e-16` | **`4.440892e-16`** ← bit-identical to the control: **exactly blind** |
| CS3 direct `A(ψ₁+ψ₂)` | `8.881784e-16` | **`1.288361e+00`** ← reds |

This independently corroborates the carve's own claim in
`tests/sn/operators/test_declared_law_is_linear.py:53-58` (*"the retired
base-point spelling could NOT catch the affine regression … the direct spelling
is strictly sharper"*), and its measured battery (the `affine` mutation reddened
**19 of 69** rows and *"neither base-point-independence row reddened"*).

✅ **STATUS: repaired in all 5, concurrently.** I verified the fix touched the
**assertion**, not merely the comment — `git show f43758d8:<f>` vs working tree:

```
-  lhs = op.apply(psi1 + lam * (psi2 - psi1))   # (1−λ)ψ₁ + λψ₂, a flux
+  lhs = op.apply(psi1 + psi2)
```

confirmed in all five files. ⚠ **The residue is S3**: `test_operators_apply_typed.py`
got the assertion upgrade and the body-comment rewrite, but its **docstring**
(`:311`, `:318`) still describes the *old, blind* contract — so the file now
MISDESCRIBES a correct gate in the direction that invites someone to "restore"
the weaker form.

### T2 ✅ CLEAN — production difference-then-route

Checked every production site that forms an iterate difference:
`orpheus/numerics/iteration.py:772` (`psi + self.corrector.apply(psi - psi_prev)`)
and `:779` (`_principal_bulk_leaf(psi - psi_prev)`). Both are plain `V`
arithmetic; the leaf finder is the renamed `_principal_bulk_leaf`
(`iteration.py:412`), and `_flux_displacement_leaf` is **gone from the tree**.
No difference is routed through a distinct "not-a-member" type anywhere in
`orpheus/`.

### T3 ✅ CLEAN — the `Σλᵢ = 1` partition-of-unity ceremony

`affine_combination` and every Σλ=1 guard on the flux path are gone. All 40+
surviving "partition of unity" hits belong to the **energy-condensation overlap
table** (`numerics/basis/overlap_basis.py`, `data/energy_grid.py`,
`foundations/frame.rst`) — an unrelated concept with its own live `:label:`
`energy-condensation-partition-of-unity`. **Not a hit; correctly untouched.**

### T4 ✅ CLEAN — the ambiguity my brief warned about is intact

`affine-bc-form` (the affine BOUNDARY law `γ₋ψ = R G γ₊ψ + q`) is **LIVE** with
`:eq:`/`:ref:` citers across `boundary_conditions.rst`, `methods/sn/*.rst`,
`curvilinear_numerics.rst`, `operator_algebra.rst` and `field_algebra.rst:79`.
Untouched by the carve and by any of my findings. The three genuinely retired
labels (`affine-torsor-algebra`, `affine-contraction-ratio`, `affine-true-error`)
are **gone with 0 dangling citers.**

### T5 ✅ CLEAN — retired-message test pins

Before recommending any reword of the retired-taxonomy runtime strings I grepped
the shortest distinctive fragments (`coding-standards`: tests pin SUBSTRINGS, so
grepping your own longer wording is strictly less sensitive). `[M]` **no test
pins them**: the only near hits are
`tests/transport/test_timed_full_field.py:117,144`
(`match="boundary must be a BoundaryField"`), which cannot match
`"…must be a RadialCharacteristicBoundaryField"`. Rewording is safe.
*(Moot for the RC composite — those strings were fixed concurrently — but the
check is what licensed the S1/S2 recommendation.)*

### T6 ✅ CLEAN — live Sphinx roles at a deleted target

`P8` (all 8 Python-domain roles × the retired vocabulary) found exactly **2**
hits, both `:meth:`…Displacement.{true_error_estimate,where_largest}`` in
`tests/numerics/test_si_diagnostic_trajectory.py:9-10`. **Both were demoted to
double-backtick literals in the concurrent sweep.** `[M]` re-run at
`755f99b5`: **0 live roles at a deleted target, tree-wide.**
⚠ Worth recording *why* this class is dangerous and why grep is the only gate:
those roles live under `tests/`, which no `automodule` renders, so **no Sphinx
severity — including `-n` — could ever have warned**.

### T7 ⚠ NOT a torsor shape, but found while checking one — a dead successor pointer

S7: `tests/transport/test_radial_characteristic_field.py:146` names
`test_subtraction_mints_a_displacement_composite_per_block` as the surviving
successor for the "role identity" claim. `[M]` that name exists **nowhere** in
`orpheus/ tests/ docs/ .claude/` outside the comment itself. Reported above.

---

## Summary

**12 MUST-FIX rows survive at `HEAD = 755f99b5`, working tree 18:28:51**
(S1–S12; 13 sites counting S8's sibling). 2 are PRODUCTION docstrings
(`_bases.py`), 6 are in `tests/`, 4 are in `.claude/{agents,skills}`. **The
highest-severity are S9–S11** — three `AGENT.md` files that teach or *instruct*
the overturned design in the present tense to every future dispatch — and
**S3**, a docstring that survived the concurrent repair of its own function body.
Two AMBIGUOUS rows are argued as KEEP. Six torsor-shape checks: **T1 real
(repaired, with the blindness measured), T2/T3/T4/T5/T6 clean.** The one
uncovered surface with the same hazard class as S9–S11 is
`.claude/agent-memory/` (**182 lines, ~20 files**) — out of brief, named so its
absence is not read as clean.

---

### ✅ Final verification stamp

`[M]` **18:34:55, `HEAD = a740d7ba`** (the tree advanced again mid-write —
`f43758d8` → `755f99b5` → `a740d7ba`). All 12 survivor rows re-checked by
predicate at that instant: **10 of 10 predicates still SURVIVE** (S1+S2, S3, S4,
S5+S6, S7, S8, S9, S10, S11, S12). Nothing in this memo's survivor list has been
remediated as of the stamp. Re-run the one-liner below before acting — this tree
is under concurrent edit and every claim here is a point-in-time snapshot:

```zsh
cd /Users/rodrigo/git/nuclear/ORPHEUS
command grep -nE 'RadialCharacteristic(Interior|Boundary)Displacement' orpheus/transport/fields/_bases.py
command grep -n  'illegal on affine flux STATES'  tests/sn/operators/test_operators_apply_typed.py
command grep -n  'FluxRole sibling'               tests/sn/primitives/test_typed_source_sinks.py
command grep -n  '# AngularBoundaryDisplacement'  tests/transport/fields/test_angular_boundary_flux.py
command grep -nE 'test_subtraction_mints_a_displacement|torsor returned' tests/transport/test_radial_characteristic_field.py
command grep -n  'FluxDisplacement'               .claude/agents/explorer/AGENT.md
command grep -n  'torsor displacement type'       .claude/agents/cross-domain-attacker/AGENT.md
command grep -n  'displacements get a shared mixin' .claude/agents/elegance-enforcer/AGENT.md
command grep -n  'still cites .FluxDisplacement'  .claude/skills/coding-elegance/SKILL.md
```

⚠ **Re-stamped at return time: 18:38:43, `HEAD = d79adb27`** (tree advanced
again: `a740d7ba → 36bfeeb2 → d79adb27`). **9 of 10 predicates survive** —
**S11** (`elegance-enforcer/AGENT.md:185`, the displacement-mixin ruling) was
remediated in the interval. Surviving MUST-FIX count: **11 rows / 12 sites**
(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S12). Re-run the block above before
acting.
