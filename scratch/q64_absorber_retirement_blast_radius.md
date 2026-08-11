# Q5.6.4 — the `[½, 1]` absorber retirement: blast-radius audit

**Scope.** Retirement of the cylinder `[½, 1]` absorber in
`orpheus/sn/sweep/pole_angular_closure.py::morel_montry_tau_per_level`
(lines 920–931), which makes that function a pure pass-through to
`morel_montry_tau_raw_per_level`.

**Audit basis.** `.claude/rules/coding-standards.md` three-search retirement
audit (graph callers + text grep + concept grep) + the five lettered questions.

**Provenance.**
- `[M]` branch `refactor/operator-strategy-layers`, HEAD `d68af438`.
- `[M]` Nexus graph built from commit `14d0280b` (dirty), which **IS an
  ancestor of HEAD** (`git merge-base --is-ancestor` → true). The 4 commits
  since the build touch **none** of `orpheus/sn/sweep/`,
  `orpheus/sn/mesh/`, or the cylindrical tests — only
  `orpheus/{cp,diffusion,moc,numerics}` + `orpheus/sn/solver.py` + plan/doc
  files. **⟹ the graph is current for this audit's subject matter.**

_(report written incrementally; sections appended as each search completes)_

---
## Search (1) — Graph callers (Nexus)

`[M]` `mcp__nexus__callers(..., transitive=True)` on both symbols.

### `morel_montry_tau_per_level` — 6 callers, **exactly 1 production**

| caller | file:line | class |
|---|---|---|
| `MorelMontryAngularSweep.__init__` | `orpheus/sn/sweep/pole_angular_closure.py:1116` (call at `:1143`) | **PRODUCTION — the only one** |
| `test_sphere_tau_matches_independent_reference` | `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:57` | test |
| `test_cyl_tau_clamp_is_the_only_difference_from_reference` | `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:85` | test |
| `test_identity_closure_tau_is_neutral_one` | `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:141` | test |
| `test_alpha_and_tau_are_bit_identical_across_tie_breaks` | `tests/sn/verification/mms/test_mms_ordering_blindness.py:226` | test |
| `test_the_degenerate_eta_pair_splits_tau_into_one_and_one_half` | `tests/sn/verification/mms/test_mms_ordering_blindness.py:249` | test |

### `morel_montry_tau_raw_per_level` — 11 callers (transitive), 1 production

Direct: `morel_montry_tau_per_level` (`pole_angular_closure.py:920`) + 5 tests
(`tests/sn/sweep/test_march_start_structure.py:143`,
`tests/sn/sweep/test_tau_arc_wellposedness.py:83/113/165`,
`tests/sn/verification/mms/test_mms_ordering_blindness.py:249`). The other 5 are the
transitive closure through `morel_montry_tau_per_level`.

### ⚠ The graph's known blind spots — cross-checked

The call graph misses **property-reached leaves**, **class-name bypass consumers** and
**direct constructors**. Cross-checks run:

- `[M]` **Property-reached:** `MorelMontryAngularSweep.__init__` stores
  `self._tau_per_level` (`pole_angular_closure.py:1153`, graph degree 1). Every downstream
  τ read is an ATTRIBUTE read on the closure, invisible as a `calls` edge to the producer.
  So "1 production caller" understates the true consumer set: it is 1 production
  **producer-call**, feeding an attribute consumed all over the sweep.
- `[M]` **Reimplementation bypass (the graph can NEVER see this):**
  `tests/sn/sweep/core/_c_surrogate.py:137` recomputes the clamp independently as
  `np.clip(tau_raw[...], 0.5, 1.0)` over the *independent* reference producer
  `orpheus.derivations.discrete.sn.contamination.morel_montry_weights`. Zero graph edges to
  either retired-code symbol, and it is the **highest-value hit in this whole audit** — see
  §B and the MUST-FIX list.
- `[M]` **Docs:** 0 graph edges; 8 `:func:` refs found only by grep (§2).

---

## Search (2) — Text grep of both symbols

Command: `grep -rn '<name>' --include=*.py --include=*.rst --include=*.md .`
(excluding `.venv`, `_build`, `.git`). Classification below.

### Production calls / definitions (`orpheus/`)

| file:line | kind |
|---|---|
| `orpheus/sn/sweep/pole_angular_closure.py:760` | **definition** `morel_montry_tau_raw_per_level` |
| `orpheus/sn/sweep/pole_angular_closure.py:863` | **definition** `morel_montry_tau_per_level` |
| `orpheus/sn/sweep/pole_angular_closure.py:920` | production call (the pass-through) |
| `orpheus/sn/sweep/pole_angular_closure.py:1143` | production call (closure ctor) |
| `orpheus/sn/sweep/pole_angular_closure.py:1633-1634` | **both in `__all__`** (see §E) |
| `orpheus/sn/sweep/pole_angular_closure.py:185, 618, 769, 782, 833, 857, 876, 904` | docstring/comment references |
| `orpheus/geometry/reduced_operator.py:511, 755, 864` | comment references (see §3) |

### Test calls

| file:line | kind |
|---|---|
| `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:45, 66, 99, 172` | test call (import + 3 calls) |
| `tests/sn/verification/mms/test_mms_ordering_blindness.py:228, 240, 261, 264-265` | test call |
| `tests/sn/sweep/test_tau_arc_wellposedness.py:62, 101, 136, 179` | test call (raw only) |
| `tests/sn/sweep/test_march_start_structure.py:38, 158` | test call (raw only) |
| `tests/sn/sweep/core/test_cell_visit_c_stamp.py:47` | prose mention (names it as the path NOT used) |
| `tests/sn/sweep/curvilinear/test_alpha_closed_form.py:231` | prose mention (raw's `eta_edge`) |

### Doc cross-references (Sphinx `:func:` — the SILENT class)

| file:line | target |
|---|---|
| `docs/theory/methods/sn/curvilinear_one_group.rst:517, 538, 799, 1050` | `morel_montry_tau_per_level` |
| `docs/theory/methods/sn/curvilinear_one_group.rst:3890` | `morel_montry_tau_raw_per_level` |
| `docs/theory/foundations/structured_geometry.rst:286, 288` | both (adjacent lines) |
| `docs/theory/verification/sn.rst:1197` | `morel_montry_tau_per_level` |
| `docs/theory/methods/sn/history.rst:983` | `morel_montry_tau_per_level` (past-tense history) |

⚠ **If the collapse DELETES `morel_montry_tau_per_level`**, all 7 present-tense `:func:`
refs above become dangling — and per `coding-standards.md` they render as **plain text with
no `-W` warning** (`orpheus/sn/sweep/` has no `automodule`; verify before relying on `-n`).

### Prose mentions outside code/docs (plans, memory, scratch — informational)

`.claude/plans/quadrature_machinery_campaign.md:1323, 1374, 2672`;
`.claude/plans/phase2_code_prose/{classify_P2-F.md:61-62,89, P2-F_curvilinear.md:80,95}`;
`.claude/plans/{facefield_codim1_design.md:295,481, stencil_assembly_dsa_roadmap.md:885,1939}`;
`.claude/agent-memory/**` (5 files); `scratch/{issue326_*,q5_6_3_flip_call_sites.md}`.
Not blocking (archaeology), but `q5_6_3_flip_call_sites.md:341-383` is the design note that
scoped this very retirement.

**0 hits** in `derivations/` (top-level) and `orpheus/derivations/` for either symbol.

---
## Search (3) — the CONCEPT grep (paraphrase, not symbol)

The absorber is documented by paraphrase far more than by symbol. Six probe families run.
Classification: **(a)** present-tense claim that becomes FALSE ⟹ MUST-FIX, **(b)** past-tense
history that correctly STAYS, **(c)** unrelated.

### 3.1 `max(0.5` / `min(1.0` / `np.clip(..., 0.5, 1.0)` — the literal expression

| file:line | text | class |
|---|---|---|
| `orpheus/sn/sweep/pole_angular_closure.py:929` | the absorber itself | **the retirement** |
| `orpheus/sn/sweep/pole_angular_closure.py:902` | "T27 adjudicated the fused `max(0.5, min(1.0, ·))` as TWO objects" | **(b)** history — but the sentence *continues* into a present-tense "and this absorption, which exists for…" ⟹ split needed |
| `orpheus/sn/sweep/pole_angular_closure.py:924` | "bit-identical to the pre-split inline `max(0.5, min(1.0, tau_raw))`" | deleted with the code |
| **`tests/sn/sweep/core/_c_surrogate.py:137`** | `tau = float(np.clip(tau_raw[...][...], 0.5, 1.0))` | ⛔ **(a) EXECUTABLE — the surrogate ORACLE re-implements the absorber. This is the #1 MUST-FIX.** |
| `tests/sn/sweep/core/_c_surrogate.py:133-135` | "the production cylinder τ is clamped to `[½, 1]`" | **(a)** |
| `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:98` | `tau_clamped_ref = [np.clip(t, 0.5, 1.0) …]` | ⛔ **(a) EXECUTABLE — gate goes RED** |
| `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py:93` | `tau_clamped[n] = max(0.5, min(1.0, raw))` | **(c)** — the **W1 SPHERE** clamp study (sphere was unclamped 2026-06-13). Verify, don't assume: this file is about the sphere, not the cylinder absorber |
| `tests/sn/sweep/test_tau_arc_wellposedness.py:3` | "T27 adjudicated the cylinder's fused `max(0.5, min(1.0, τ_raw))` as TWO objects" | **(b)** history |
| `docs/theory/methods/sn/curvilinear_numerics.rst:2441` | "The τ-clamp (`τ → max(0.5,min(1.0,τ_raw))`) **breaks** the *exact* linear-in-μ threading wherever it is active" | **(a)** present tense |
| `docs/theory/foundations/structured_geometry.rst:349` | "The fused cylinder expression `max(0.5, min(1.0, τ_raw))` **welds**…" | **(a)** present tense |

### 3.2 the interval spelled `[½, 1]` / `[1/2, 1]` / `[\tfrac12, 1]` / `[0.5, 1]`

`orpheus/`: `pole_angular_closure.py:739, 768, 832, 878, 900, 909, 923, 1067`;
`geometry/reduced_operator.py:115, 126, 135, 171, 864`.
`tests/`: `test_curvilinear_aniso_convergence.py:120`; `test_tau_arc_wellposedness.py:15`;
`regression/_generate_snapshots.py:461`; `core/_c_surrogate.py:134`;
`core/test_ordinate_scan_reset.py:264`; `verification/analytical/test_phase_c_crosscheck.py:211`;
`curvilinear/test_compute_psi_half_per_level.py:258`; `verification/mms/test_mms_ordering_blindness.py:254`.
`docs/`: `curvilinear_one_group.rst:69, 74, 498, 518, 549, 1155, 4016, 4084`;
`structured_geometry.rst:220, 328, 335, 351, 352`; `verification/sn.rst:1157, 1158, 1313`;
`foundations/discretization.rst:32` **(c)** — that one is the **LD slope limiter** `w = 1/(1+k) ∈ [½,1]`, a different object entirely.

⚠ `pole_angular_closure.py:1067` — *"reduce to pure DD at τ = 1/2 (the M-M clamp **keeps**
τ ∈ [1/2, 1])"* — **(a)**, and it is in the `MorelMontryAngularSweep` CLASS docstring, ~150
lines from the retired lines. Exactly the "pre-existing line outside the diff" class.
Post-retirement `τ ∈ [1/5, 4/5]` on the folded family, so **the parenthetical inverts**.

⚠ `tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py:258` — *"For any τ ∈ [1/2, 1]
the M-M weighted DD recurrence…"*: a **property-test domain statement**. Check whether the
generated τ range is a `hypothesis`/parametrize domain (⟹ widen to `(0, 1]`) or only prose.

### 3.3 `absorber` (the word) — heavily overloaded

`[M]` ~150 hits, of which the **vast majority are the physics sense** ("pure absorber",
"cavity absorber", "thick absorber") — **(c)**. The τ-sense hits are exactly:
`orpheus/sn/sweep/pole_angular_closure.py:739, 768, 770, 832, 878, 901, 914, 923`;
`tests/sn/verification/mms/test_mms_ordering_blindness.py:254`;
`tests/sn/verification/analytical/test_phase_c_crosscheck.py:211-212`;
`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py:120, 124`;
`tests/sn/sweep/test_tau_arc_wellposedness.py:10, 118, 175`;
`tests/sn/sweep/core/test_ordinate_scan_reset.py:55, 264, 285, 316`;
`tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py:23`;
`tests/sn/regression/_generate_snapshots.py:462-463`;
`docs/theory/methods/sn/curvilinear_one_group.rst:74, 498, 506, 4084`;
`docs/theory/foundations/structured_geometry.rst:390`.

⭐ **A contradiction ALREADY IN THE TREE, in one file.** `tests/sn/sweep/test_tau_arc_wellposedness.py`
says at `:15` *"the `[½, 1]` ABSORPTION — **stays** until the fold wiring (Q5.6)"* (present
tense, becomes false at retirement ⟹ **(a)**) while at `:118` and `:175` it says *"the
**retired** absorber"* (past tense — false **today**, becomes true at retirement). One file,
two tenses, one object. The retirement fixes `:118`/`:175` and breaks `:15`.

### 3.4 `clamp` / `clamped` / `unclamped`

`[M]` ~120 hits. Two large **(c)** families to NOT touch: (i) the **W1 SPHERE unclamp**
narrative (`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py:303-430`,
`docs/theory/verification/sn.rst:1150-1215`) — landed 2026-06-13 at `b2d8a6d`, a different
retirement; (ii) unrelated clamps (`numerics/iteration.py:901`, `data/energy_grid.py:279-289`,
`solver.py:1202`).

**(a)** hits in this family: `orpheus/geometry/reduced_operator.py:126-140` (the
`**CYLINDER** retains the clamp … Removing the clamp here needs a 2-D (η,φ) angular closure
(out of scope), tracked by Issue #229` bullet — **doubly stale**: the fold already removed
the reason); `reduced_operator.py:631` ("the angular closure owns τ (with the cylinder
clamp)"); `reduced_operator.py:864` ("where the cylinder [½, 1] clamp now lives");
`reduced_operator.py:704` ("Bailey-Morel-Chang (2010) Eq. 5 (Morel--Montry τ clamp)" — also
a wrong Eq. number vs. the Eq. 43 used everywhere else);
`orpheus/numerics/spaces/radial_characteristic_space.py:170` ("#229 (the cylinder τ clamp
fact)"); `orpheus/transport/spatial/diamond.py:42` ("justification for the M-M
weighted-diamond τ clamp"); `tests/sn/operators/test_radial_characteristic_metric.py:50`
(same "#229 … cylinder τ clamp fact" phrasing);
`tests/sn/sweep/core/test_cell_visit_c_stamp.py:48, 98`;
`tests/sn/sweep/core/test_diamond.py:584, 615`; `tests/sn/sweep/core/test_ordinate_scan.py:610`.

### 3.5 `tau_raw` / `τ_raw`

`[M]` The notation survives the retirement **as a concept** (it is BMC Eq. 43's fractional
position). But once nothing is clamped, "raw" has no antonym: `docs/` carries 16 `τ_raw`
spellings (`structured_geometry.rst` 5, `curvilinear_one_group.rst` 10,
`angular_quadrature.rst` 1) whose "raw" qualifier becomes vestigial. See §E on the naming.

### 3.6 `0 hits` results — explicitly recorded

- `[M]` `[½, 1]` / `[1/2, 1]` in `derivations/` (top level): **0 hits**.
- `[M]` `morel_montry_tau*` in `derivations/` and `orpheus/derivations/`: **0 hits**.
- `[M]` `morel_montry_tau*` in any `__init__.py` anywhere in `orpheus/`: **0 hits** (§E).
- `[M]` `automodule`/any autodoc directive for `orpheus.sn.sweep.pole_angular_closure` in
  `docs/api/`: **0 hits** ⟹ **the module is NOT rendered**, so Sphinx `-n` would NOT catch a
  dangling `:func:` to it. Grep is the only gate here (`coding-standards.md`, 2026-08-03).

---
## A. Is `assert_carrying_quadrature` really on the cylindrical `SNMesh` path?

**YES — unconditionally, on every constructor.** `[M]`

- **Definition:** `orpheus/sn/sweep/pole_angular_closure.py:672`.
- **The single call site:** `orpheus/sn/mesh/augmented_mesh.py:338`
  (import at `:63`), inside `case CoordSystem.CYLINDRICAL:` of the `match self.coord:`
  block, placed **deliberately AFTER** `cylindrical_streaming(mesh, quadrature)` at `:326`
  so the older structure-less guard keeps ownership of slab/sphere cubatures
  (`augmented_mesh.py:327-337` comment; message-disjointness is asserted in the
  `assert_carrying_quadrature` docstring at `:696-699`).

### The construction path — all three entries converge on ONE body

| entry | file:line | reaches `_init_core`? |
|---|---|---|
| `SNMesh.__init__` (legacy `Mesh1D`/`Mesh2D` surface) | `augmented_mesh.py:202` → `:218` | **yes** |
| `SNMesh.from_axes` (axis-native) | `augmented_mesh.py:634` → `:698` `cls.__new__` → `:699` | **yes** |
| `SNMesh.from_material_mesh` (homogenization re-solve) | `augmented_mesh.py:711` → `:750` `cls.__new__` → `:751` | **yes** |

`_init_core` is defined at `augmented_mesh.py:228`; the guard at `:338` lives inside it.
The two `cls.__new__(cls)` calls look like `__init__`-bypasses but immediately call
`_init_core`, so they are **not** guard-bypasses.

`self.coord` is derived by `MaterialMesh._init_data` from the axes — a caller cannot
declare `CYLINDRICAL` and route around the `case`.

### Is there ANY bypass? — searched, and the answer is "not through the public surface"

`[M]` greps returning **0 production hits**: `SNMesh.__new__` outside the two classmethods;
`class …(SNMesh)` subclasses; post-construction `sn_mesh.quad = …` reassignment;
`dataclasses.replace(...)` on an `SNMesh`. `MorelMontryAngularSweep(...)` is only ever
constructed with an already-built `SNMesh` (`tests/sn/sweep/curvilinear/test_pole_angular_closure.py:74, 84`
— both spherical); the unbound `MorelMontryAngularSweep()` legacy mode was **retired at C5,
2026-07-03**. The `pole_angular_closure=` override (`:225, :706, :758`) passes a CLASS and
does not touch admission.

The only routes that skip the guard are **deliberate test instrumentation**:
`scratch/mutation_battery_admission.py:69, 119, 137, 169` (`mock.patch.object`) and
`tests/sn/mesh/test_cylindrical_quadrature_admission.py:303` (calls the guard directly with
`object()`). Neither is a production path.

---

## B. Does anything call the CLAMPED producer with a NON-admitted cylindrical quadrature?

**YES — 2 tests, and one of them asserts the absorbed values numerically.** `[M]`

| call site | quadrature | admitted? | verdict |
|---|---|---|---|
| `orpheus/sn/sweep/pole_angular_closure.py:1143` (production) | `sn_mesh.quad` | **always admitted** (§A) | safe |
| `tests/.../test_mms_ordering_blindness.py:265` | `Quadrature.product(n_mu=2, n_phi=8)` | **NO** (degenerate levels) | ⛔ **RED** |
| `tests/.../test_mms_ordering_blindness.py:240` | `Quadrature.product(n_mu=4, n_phi=8)` | **NO** | value changes, test likely **stays green** |
| `tests/.../test_tau_producer_equivalence.py:99` | `folded_product(4, n_phi)` | yes | ⛔ **RED** (§below) |
| `tests/.../test_tau_arc_wellposedness.py:101/136/179` | hand-built `SimpleNamespace` + folded | n/a — calls **raw** only | unaffected |
| `tests/sn/sweep/test_march_start_structure.py:158` | hand-built `quadlike` | n/a — **raw** only | unaffected |

### `[M]` measured τ movement (`.venv/bin/python`, this tree, HEAD `d68af438`)

```
folded_product(4,8)  : 4/4 levels MOVE.  raw [0.2195 0.4142 0.5858 0.7805]
                                          cl  [0.5    0.5    0.5858 0.7805]   max|Δτ| = 2.804e-01
folded_product(4,16) : 4/4 levels MOVE.  4 of 8 entries per level are absorbed to 0.5
                                                                                max|Δτ| = 2.953e-01
folded_product(2,4)  : 2/2 levels MOVE.  raw [0.2929 0.7071] → cl [0.5 0.7071] max|Δτ| = 2.071e-01
product(2,8) [refused at admission]: raw [0,1,0,1,…] → cl [0.5,1,0.5,1,…]      max|Δτ| = 5.000e-01
```

⟹ the absorber is **NOT** a rarely-firing safety net on the admitted family: it moves
**every level of every production folded rule**, by up to **0.30 absolute in τ**. Every
cylindrical numeric baseline moves.

### The three tests that go RED by construction

1. ⛔ **`tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:85`
   `test_cyl_tau_clamp_is_the_only_difference_from_reference`** — its whole subject is the
   absorber. It asserts `closure τ == np.clip(ref τ_raw, .5, 1)` (`:114`) and, where the
   clamp bites, `closure τ != ref τ_raw` (`:123`). Post-retirement both invert. Its
   guard-the-guard at `:103` (`any(t.min() < 0.5)`) still passes (0.2195 < 0.5), so it fails
   *loudly on the value assertions*, not vacuously.
   **Migration (not deletion):** rewrite to the SPHERE leg's shape —
   `np.testing.assert_array_equal(closure_τ, morel_montry_weights(quad, "cylindrical"))`,
   0-ULP. That stays a genuine vv-L11 gate (two independent code paths to BMC Eq. 43); the
   NEGATIVE control dies with the transform it pinned, and the module docstring's bullet at
   `:22-24` must go with it.
2. ⛔ **`tests/sn/verification/mms/test_mms_ordering_blindness.py:249`
   `test_the_degenerate_eta_pair_splits_tau_into_one_and_one_half`** — `:268` asserts
   `clamped == [0.5, 1, 0.5, 1, …]`. Post-retirement `clamped is raw == [0, 1, 0, 1, …]`,
   which the `:267` line already asserts. **The test's own name and thesis dissolve**: with
   no absorber there is no "split into one and one half", only the double-cover fingerprint.
   Rename/re-scope to the `[0,1,0,1,…]` fingerprint (which is the durable content — it is
   the guard's documented blind spot) and drop the second assertion + `:254`'s prose.
3. ⛔ **`tests/sn/sweep/core/_c_surrogate.py:137`** — the **independent surrogate oracle**
   for `c_in`/`c_out` replicates the absorber. It is NOT a caller of either retired symbol
   (so **invisible to Nexus AND to a symbol grep**), it is reached by
   `test_cell_visit_c_stamp.py`, `test_diamond.py:584/615`, `test_ordinate_scan.py:610`.
   Leaving the `np.clip` in place turns the oracle into a **wrong** reference and reddens
   every cylinder arm of those three gates; removing it keeps them honest. Prose at `:134-135`
   goes with it.

`[M]` **0 hits** in `derivations/` or `orpheus/derivations/` — no diagnostic calls either
producer. `orpheus/derivations/discrete/sn/contamination.py::morel_montry_weights` (`:156`)
is the independent reference and its cylindrical branch is **already unclamped** (`:32-46`),
so the retirement makes production and reference agree exactly — a *strengthening*, not a
loss, provided item 1 above is migrated rather than deleted.

---

## C. Snapshot / regression baselines that MOVE

`[M]` Every one of these was **pre-declared** by the campaign as moving at 6.4 — the tree
already carries the notice.

| baseline | why it depends on τ | evidence |
|---|---|---|
| `tests/sn/regression/snapshots/cyl_1g_homogeneous_folded_4x8_dd_n20.npz` | DD cylinder solve on `folded_product(4,8)`; τ enters `c_in`/`c_out` on every visit | generator `tests/sn/regression/_generate_snapshots.py:466`; the block comment at `:461-465` says verbatim *"these baselines move again at Q5.6.4 (the absorber retirement…)"* |
| `…/cyl_1g_homogeneous_folded_2x4_dd_n20.npz` | same, `folded_product(2,4)` — `[M]` 2/2 levels move | `_generate_snapshots.py:471` |
| `…/cyl_2g_3reg_folded_4x8_dd_n40.npz` | same, 2G 3-region | `_generate_snapshots.py:476`; also in `_TRUNCATED_INNER_CASES` (`tests/sn/regression/test_dd_regression.py:100`) |
| `…/walk_matvec_cyl_2g.npz` | the `(L+C)` **matvec + adjoint** walk baseline; fixture `_make_cyl` = `folded_product(4,8)` | `tests/sn/regression/_generate_walk_baselines.py:143`; fixture `tests/sn/operators/test_g_adjoint_reciprocity.py:167-186` (docstring `:174` names the file) |

**Consumers that will red until re-capture:** `tests/sn/regression/test_dd_regression.py`
(3 cyl cases), `tests/sn/regression/test_walk_matvec_baselines.py` (the `cyl_2g` row).

**A second, frozen-copy consumer** — easy to miss because it does not read the `.npz`:
`tests/sn/verification/analytical/test_phase_c_crosscheck.py:213-215` hardcodes the three
cylinder `k_eff` values (`1.5`, `1.4999999999999996`, `1.2302082296342958`) into
`_SNAPSHOT_KEFFS`, explicitly *"decoupled from snapshot-file existence"* (`:205-207`). Its
comment at `:211-212` already says they *"move again at Q5.6.4's absorber retirement"*. **The
re-capture MUST update this dict too** — regenerating the `.npz` alone leaves a stale
authoritative-looking literal.

**NOT affected** (verified by geometry, not assumed): all `slab_*`, `sphere_*`, `2d_*`,
`walk_matvec_{slab,sphere,cart2d}` baselines, and `tests/geometry/snapshots/*` (BC
equivalence, Cartesian/Lebedev). Cartesian uses `IdentityAngularClosure` (τ ≡ 1, no
absorber); the sphere branch returns `raw` today and after (`pole_angular_closure.py:921-922`).

**How I know without running them:** the absorber is applied **only** on the
`CoordSystem.CYLINDRICAL` branch (`pole_angular_closure.py:921-931`), τ is consumed by the
closure's `c_in`/`c_out` and the ψ½ recurrence (`_psi_half_grid_single_level:973-977`,
divides by `tau_m`), and the four baselines above are the only cylindrical artefacts in
`tests/**/snapshots/`. Two independent in-tree notices (`_generate_snapshots.py:461-465`,
`test_phase_c_crosscheck.py:211-212`) corroborate.

### Non-snapshot numeric gates that also re-baseline

- `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py:120-125` — the live
  azimuthal-floor ladder `3.538e-3 → 6.782e-4` is stamped *"the [½,1] absorber still live"*
  and explicitly hands the re-measurement to 6.4. **This docstring is the ladder's single
  source** — do not copy its numbers elsewhere (`plan-authoring` §9).
- `tests/sn/sweep/core/test_ordinate_scan_reset.py:283-317`
  `test_pole_resonance_unreachable_on_admitted_family` — a **deliberate 6.4 tripwire**: it
  asserts the ERR-054 pole resonance stays unreachable when `c_out` moves. Green ⟹ nothing
  owed; RED ⟹ 6.4 owes a live reproducer (its own failure message says so). Companion notice
  at `tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py:20-27`.

---
## D. Sphinx pages + sections that need editing

`[M]` Five theory pages. `docs/theory/verification/matrix.rst` is **auto-generated**
(`:1-11` "Do not edit by hand") and needs no manual change.

### D.1 `docs/theory/foundations/structured_geometry.rst` — the equation-of-record page

| line | text | action |
|---|---|---|
| `:220-221` | "unclamped raw weight on the sphere; `[1/2, 1]`-clamped on the cylinder — see `:eq:`morel-montry-clamp`" | **(a)** rewrite |
| `:252` | "the production code applies a **geometry-dependent** clamp on top of it" | **(a)** rewrite |
| **`:250-265`** | **the `:label: morel-montry-clamp` equation** — a two-branch `cases` block whose cylinder arm is `clip(τ_raw, ½, 1)` | ⛔ **the single highest-value edit.** After 6.4 both geometries are `τ = τ_raw`, i.e. there is no `cases` and no split. **Keep the `:label:`** (see D.5) and rewrite the body; the name `morel-montry-clamp` becomes a misnomer — renaming it is a second, `-W`-gated change (`:eq:` refs DO warn) |
| `:266-275` | the `vv-status: morel-montry-clamp documented` rationale — cites `test_tau_producer_equivalence.py` as the verifiable content and names "the production clamp policy" | **(a)** — and the gate it names is one of the RED-by-construction tests (§B.1) |
| `:283-290` | the `_tau-ownership-note`: names `morel_montry_tau_raw_per_level` "(raw)" **and** `morel_montry_tau_per_level` "(the clamp policy applied)" | **(a)** — and §E decides whether both names survive |
| `:334-346` | "**CYLINDER** retains the clamp … *pending retirement* … the clip survives as pure behavioural debt (it still alters every folded level's τ)" | **(a)** → past tense. `[M]` my measurement (§B) **confirms** the "alters every folded level" claim: 4/4 |
| `:347-352` | "**The clamp is TWO objects, and the fold retires one** … The fused cylinder expression `max(0.5, min(1.0, τ_raw))` **welds**…" | **(a)** present tense → past |
| `:386-399` | "…so it **RETIRES** in Q5.6's absorber step (6.4).  (An earlier version of this passage said…)" | **(a)** → "RETIRED at `<hash>`". ⭐ Note this paragraph **already** carries a correction-in-place of its own earlier claim — extend, don't overwrite (`plan-authoring` §3) |
| `:437` | `spherical_streaming` "…and Morel--Montry closure `:eq:`morel-montry-clamp`" | **(a)** if the label is renamed |

### D.2 `docs/theory/methods/sn/curvilinear_one_group.rst` — the biggest surface

| line | text | action |
|---|---|---|
| `:69-77` | Key-Facts bullet: "**Sphere unclamped**; **cylinder clipped** to `[½, 1]` … the `[½, 1]` absorption **retires with** the σ_y fold (Q5.6)" | **(a)** — a **Key Facts** entry, i.e. the first thing a future session reads |
| `:491-514` | "**Cylinder** retains the clip `τ_n → clip(τ_n, 0.5, 1.0)` … so it retires in Q5.6's absorber step (6.4) … the clip is pure behavioural debt" | **(a)** → past tense |
| `:516-519` | "The clamp lives **inside the angular closure**, in `morel_montry_tau_per_level` (sphere unclamped, cylinder clipped to `[½, 1]`)" | **(a)** + a `:func:` ref (§E) |
| `:538-556` | `§sn-tau-closure-owned` (a `.. todo::` block): "Step A is BIT-IDENTICAL … (sphere unclamped, cylinder clamped to `[½, 1]`), so the geometry factory **still bakes** an IDENTICAL τ … the gate pins … to (a) the geometry-factory value (0-ULP) AND (b)…" | ⚠ **ALREADY STALE, pre-existing** — Step C deleted the geometry producer and the `*_equals_geometry_factory_0ulp` legs (`test_tau_producer_equivalence.py:8-15` says so). Fix while you are here |
| `:799` | "…already computes τ from it in `morel_montry_tau_per_level`" | `:func:` ref (§E) |
| `:1044-1062` | "the closure's `morel_montry_tau_per_level` is a line-for-line replica of the geometry factory's τ arithmetic … pinned … to the geometry-factory value (0-ULP) *and* to the structurally-independent reference" | ⚠ same **pre-existing** staleness (leg (a) is gone) |
| `:1095-1104` | the production-stamp catcher paragraph: "re-pointed onto `morel_montry_weights` (**with the cylinder clamp replicated**)" | **(a)** — this is the doc-of-record for `_c_surrogate.py:137` (§B.3) |
| `:1150-1160` | "…with the cylinder clamp `clip(·, ½, 1)` **replicated in the test surrogate**" | **(a)** same |
| `:4010-4016` | "the cylinder's edge-inclusion is a property of the *circle*, **not** the `[½, 1]` Morel–Montry clamp (the clamp **is** a separate recurrence-weight stabiliser…)" | **(a)** — the parenthetical's present tense |
| `:4082-4085` | "…and the τ = 0 division block the `[½, 1]` absorption **existed to hide**" | **(b)** past tense — **STAYS**. Good example of the discriminator |
| `:3890` | `:func:` ref to `morel_montry_tau_raw_per_level` in the R12a trichotomy | (§E) |

### D.3 `docs/theory/methods/sn/curvilinear_numerics.rst`

- `:2438-2446` `.. note::` — "The τ-clamp (`τ → max(0.5,min(1.0,τ_raw))`, Bailey–Morel–Chang)
  **breaks** the *exact* linear-in-μ threading wherever it is active. This is part of the
  residual anisotropic angular floor quantified at Issue #229" — **(a)**, and it is a causal
  claim the retirement **resolves**: the note's own mechanism is what 6.4 removes, so the
  #229 floor prediction (`test_curvilinear_aniso_convergence.py:120-125`) is testable here.
- `:960-970` (the Q5.6.3 closure paragraph, naming `assert_carrying_quadrature`) — **(b)**,
  correct as-is; a good template for the tense the 6.4 edits should adopt.

### D.4 `docs/theory/verification/sn.rst`

- `:1204-1212` — "**Geometry split.** W1 removed the clamp for the **sphere only**. The
  **cylinder keeps** it … the cylinder's real fix is a 2-D (η,φ) closure, **not** unclamping.
  See `:eq:`morel-montry-clamp`… for the equation-of-record carrying both branches" —
  **(a)**, and the *"not unclamping"* clause is the strongest present-tense-false statement
  in the corpus: 6.4 does exactly the thing this sentence says is not the fix. It needs a
  correction-in-place explaining what changed (the σ_y fold removed the `τ_raw = 0`
  structural singularity that made unclamping unsafe), not a deletion.
- `:1197` — the `:func:` ref (§E). `:1313` — "No partition … gives `τ_raw ∈ [½, 1]` with
  bounded edges" — **(c)**, a statement about *partitions*, still true.
- `:1150-1203` (the whole W1 block) — **(b)** sphere history, STAYS.

### D.5 `docs/theory/methods/sn/angular_quadrature.rst` + `history.rst`

- `angular_quadrature.rst:366-378` — **(b)/(c)**: describes what the fold buys
  (`τ_raw ∈ [⅕, ⅘]`, reversal identity, `assert_carrying_quadrature`). Already correct;
  after 6.4 it becomes the *operative* description rather than a property of the raw value.
- `history.rst:983` — **(b)** past-tense milestone naming `morel_montry_tau_per_level`;
  STAYS (but see §E if the symbol is renamed). `history.rst:34-45` — the Q5.6.3 admission
  row; the 6.4 retirement owes a **new row** here (Cardinal Rule 3).

### D.6 Reference-integrity note

`[M]` `:eq:`morel-montry-clamp`` is referenced from **`orpheus/geometry/reduced_operator.py:117`**
(a docstring, in a module with **no** `automodule`) and from `structured_geometry.rst:221, 437`.
`:eq:` refs in `.rst` **do** warn under `-W`; the one inside the un-rendered docstring does
**not**. If the label is renamed, grep — do not trust the build.

---

## E. Public-API surface and the two-function collapse

### E.1 Re-exports

`[M]` **Both names are in `__all__`** of the defining module:
`orpheus/sn/sweep/pole_angular_closure.py:1633` (`morel_montry_tau_per_level`) and
`:1634` (`morel_montry_tau_raw_per_level`).

`[M]` **0 hits** for `morel_montry` in any `__init__.py` under `orpheus/`.
`orpheus/sn/sweep/__init__.py:28` re-exports *from* `pole_angular_closure`, but **not**
these two (it takes the closure classes). So the public surface is exactly
`orpheus.sn.sweep.pole_angular_closure.<name>` — module-qualified, one hop, no package-level
alias to chase.

### E.2 Is the collapse clean? — **YES, behaviourally.** The only open question is the NAME.

`[M]` After the retirement the two functions are **extensionally identical**: same signature,
same return, same exceptions — including the Cartesian `ValueError`, which is *already*
raised by the raw producer (`:856-860`), because `morel_montry_tau_per_level` calls it on
line 920 **before** its own `coord` test. So
`tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:171-172`
(`pytest.raises(ValueError)` on Cartesian) passes under either collapse. Keeping both is a
textbook twin path (Cardinal Rule 2) — the wrapper must go.

**Measured reference sets** (`orpheus/` + `tests/` + `docs/`, excluding `.claude/`, `scratch/`):

| symbol | `orpheus/` | `tests/` | `docs/` | total |
|---|---|---|---|---|
| `morel_montry_tau_per_level` | 10 | 10 | 7 | **27** |
| `morel_montry_tau_raw_per_level` | 7 | 9 | 2 | **18** |

**Option 1 — keep the `_raw_` name, delete the wrapper.** Touches the 27 refs to the deleted
name (7 of them `:func:` in Sphinx, which fail SILENTLY — §3.6). Leaves `raw` as a qualifier
with nothing to be raw *relative to* — a naming-signal regression (`feedback_high_signal_names`).

**Option 2 — RENAME `morel_montry_tau_raw_per_level` → `morel_montry_tau_per_level`, delete
the wrapper.** ⭐ **Recommended.** Touches the 18 refs to `_raw_` (only **2** in `docs/`:
`curvilinear_one_group.rst:3890`, `structured_geometry.rst:286`) and leaves all 7 doc
`:func:` refs to `morel_montry_tau_per_level` **valid by construction** — precisely the class
that Sphinx cannot warn about. Fewer edits, and the surviving name is the honest one: there
is now exactly ONE τ.
- Cost: the `τ_raw` **notation** in prose (`[M]` 16 spellings: `structured_geometry.rst` 5,
  `curvilinear_one_group.rst` 10, `angular_quadrature.rst` 1) keeps a vestigial qualifier.
  That is a doc-notation decision, separable from the code change — but note the R12a
  trichotomy theorem (`curvilinear_one_group.rst:3880-3900`,
  `tests/sn/sweep/test_march_start_structure.py:143`) is *stated* in `τ_raw` and stays
  mathematically correct either way.
- ⚠ Whichever option: `_assert_tau_within_unit_interval` (`:728`) and the raw producer's
  `Raises` section (`:788-807`) currently justify themselves as "promoted from the `[½, 1]`
  absorber's silent absorption" — **(b)** history, keep; but `:768-770` ("Split out of
  `morel_montry_tau_per_level` because the raw value carries structure the absorber
  destroys") is **(a)**: after the collapse there is no sibling to be split from.

---
## Gaps / things I expected and did NOT find

- **No `@pytest.mark.catches(...)` or `@pytest.mark.verifies(...)` on either RED-bound test.**
  `[M]` `test_tau_producer_equivalence.py` carries only `@pytest.mark.foundation` per test;
  `test_mms_ordering_blindness.py` has file-level `pytestmark = pytest.mark.foundation`. So
  the marker-migration obligation of `coding-standards.md` is **empty here** — nothing in
  `tests/_harness/` or the error catalog names them. (The one `verifies` in the neighbourhood
  is `test_tau_arc_wellposedness.py:163 → "morel-montry-folded-arc"`, on a test that
  **survives**.)
- **No `xfail(strict=True)` placeholder** anticipating the retirement anywhere — the
  campaign's usual "the xfail set IS the todo list" technique was not used for 6.4.
- **No issue-number anchor in the code for 6.4 itself.** `#229` is cited at
  `pole_angular_closure.py:909`, `reduced_operator.py:140`, `radial_characteristic_space.py:170`,
  `verification/sn.rst:1318` as the tracker — but **#229 is CLOSED** (a measurement record,
  not a work item; `plan-authoring` §9). Any new prose written at 6.4 should not re-cite it
  as "tracked by".
- **`docs/theory/methods/sn/history.rst` has no 6.4 row yet** — the retirement owes one
  (Cardinal Rule 3), next to the Q5.6.3 admission row at `:34-45`.
- **Two PRE-EXISTING stale doc claims** surfaced by this audit that are *not* caused by the
  retirement but sit inside its edit surface (fix while in the file):
  `curvilinear_one_group.rst:548-556` and `:1050-1058` both still say the producer-equivalence
  gate pins τ "to (a) the geometry-factory value (0-ULP)" — Step C deleted that leg, and the
  test module's own docstring (`test_tau_producer_equivalence.py:8-15`) records the removal.

---

## MUST-FIX list — present-tense claims that become FALSE, ranked

Rank = (will it break/mislead a run?) × (how authoritative does it look?).

| # | file:line | claim | why it ranks here |
|---|---|---|---|
| **1** | `tests/sn/sweep/core/_c_surrogate.py:137` (+ prose `:134-135`) | `np.clip(tau_raw, 0.5, 1.0)` — the absorber **executably re-implemented** in the independent oracle | **Code, not prose.** Invisible to Nexus (no edge) AND to a symbol grep (names neither retired function). Left in place it makes the oracle *wrong* and reddens the cylinder arms of `test_cell_visit_c_stamp.py`, `test_diamond.py:584/615`, `test_ordinate_scan.py:610` — which will read as "the carve broke the sweep", not "the oracle is stale" |
| **2** | `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py:22-24, 85-132` | the whole `test_cyl_tau_clamp_is_the_only_difference_from_reference` gate + its module-docstring bullet | RED by construction. **Migrate to the 0-ULP shape (§B.1), do not delete** — it is the cylinder's only vv-L11 producer gate |
| **3** | `tests/sn/verification/mms/test_mms_ordering_blindness.py:254, 268` | `assert clamped == [0.5, 1, …]` + "the `[1/2, 1]` absorber … turns that into `{1, 1/2}`" | RED by construction; the test's **name and thesis** dissolve, so a silent value edit would leave a misnamed gate |
| **4** | `docs/theory/foundations/structured_geometry.rst:250-265` | the **`:label: morel-montry-clamp` equation-of-record**, cylinder arm `clip(τ_raw, ½, 1)` | It is the *equation of record*, cited from `reduced_operator.py:117` and two other doc sites. A false equation-of-record is the worst class of doc bug in this project |
| **5** | `docs/theory/verification/sn.rst:1204-1212` | "the **cylinder keeps** it … the cylinder's real fix is a 2-D (η,φ) closure, **not** unclamping" | Directly contradicts what 6.4 does. Needs correction-in-place naming what changed (the σ_y fold removed the `τ_raw = 0` singularity), not deletion |
| **6** | `docs/theory/methods/sn/curvilinear_one_group.rst:69-77` | **Key Facts** bullet: "cylinder clipped to `[½, 1]` … retires with the σ_y fold (Q5.6)" | Key Facts is the first thing a future session reads (CLAUDE.md Cardinal Rule 3) |
| **7** | `orpheus/sn/sweep/pole_angular_closure.py:1067` | "(the M-M clamp **keeps** τ ∈ [1/2, 1])" in the `MorelMontryAngularSweep` class docstring | ~150 lines from the diff; inverts (`[M]` τ ∈ [1/5, 4/5] on the fold). Classic pre-existing-line-outside-the-diff |
| **8** | `orpheus/geometry/reduced_operator.py:126-141` | "**CYLINDER** retains the clamp … Removing the clamp here needs a 2-D (η,φ) angular closure (**out of scope**), tracked by Issue #229" | Module-header prose in a **non-`automodule`** file ⟹ no build gate at any severity. Also cites a CLOSED issue as a tracker |
| **9** | `orpheus/geometry/reduced_operator.py:631, 704, 864` | "the angular closure owns τ (with the cylinder clamp)"; "BMC (2010) **Eq. 5** (Morel–Montry τ clamp)"; "where the cylinder `[½, 1]` clamp **now lives**" | Same file; `:704` also carries a wrong equation number (Eq. 5 vs Eq. 43 everywhere else) |
| **10** | `docs/theory/methods/sn/curvilinear_one_group.rst:491-519, 1095-1104, 1150-1160, 4010-4016` | "Cylinder **retains** the clip"; "with the cylinder clamp **replicated**" (×2 — the doc-of-record for MUST-FIX #1); "the clamp **is** a separate recurrence-weight stabiliser" | Bulk prose; `:1095-1104` and `:1150-1160` must move **with** #1 or the surrogate's stale replication stays documented as correct |
| **11** | `docs/theory/methods/sn/curvilinear_numerics.rst:2438-2446` | ".. note:: The τ-clamp … **breaks** the exact linear-in-μ threading wherever it is active" | A causal claim the retirement *resolves*; it names the #229 floor, whose re-measurement is 6.4's acceptance test |
| **12** | `orpheus/sn/sweep/pole_angular_closure.py:768-770, 878, 900-918` | "Split out of `morel_montry_tau_per_level` because the raw value carries structure the absorber destroys"; "applies the cylinder's `[½, 1]` absorber on top"; the whole `Notes` paragraph ending "It RETIRES with the fold wiring (Q5.6); it must survive until then" | Deleted/rewritten with the code. `:739` ("Promoted from the `[½, 1]` absorber's silent absorption") is **(b)** history and **STAYS** |
| **13** | `tests/sn/sweep/test_tau_arc_wellposedness.py:15-19` | "the `[½, 1]` ABSORPTION — **stays** until the fold wiring (Q5.6), because production NODE_ALIGNED cylinders still start on an edge node" | Becomes false; and the SAME file already says "the **retired** absorber" at `:118` and `:175`, so the file is self-contradictory today and self-consistent after |
| **14** | `orpheus/numerics/spaces/radial_characteristic_space.py:170` · `orpheus/transport/spatial/diamond.py:42` · `tests/sn/operators/test_radial_characteristic_metric.py:50` | "#229 (the cylinder τ **clamp fact**)" ×2; "Auxiliary justification for the M-M weighted-diamond τ **clamp**" | Low blast radius, far from the diff, found only by the CONCEPT grep — exactly the class the symbol grep misses |
| **15** | `tests/sn/sweep/core/test_cell_visit_c_stamp.py:48, 98` · `test_diamond.py:584, 615` · `test_ordinate_scan.py:610` | "(clamped for cylinder)" / "(clamped cylinder τ)" in docstrings describing the surrogate | Travel with MUST-FIX #1 |
| **16** | `tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py:258` | "For any τ ∈ [1/2, 1] the M-M weighted DD recurrence…" | **Check whether this is an executable domain** (parametrize/hypothesis bounds) or prose. If executable, widen to `(0, 1]` — the recurrence divides by τ |

### Explicitly NOT to touch (past-tense history — the discriminator working)

`pole_angular_closure.py:739` ("Promoted from the `[½, 1]` absorber's **silent absorption**,
which is **how** T22's mis-ordering laundered…"); `:902` (the T27 adjudication, as history);
`curvilinear_one_group.rst:4082-4085` ("the τ = 0 division block the absorption **existed to
hide**"); `history.rst:983`; the entire W1 **sphere** unclamp record
(`verification/sn.rst:1150-1203`, `test_curvilinear_aniso_convergence.py:303-430`,
`test_w1_clamp_silent_on_flat.py`); `verification/sn.rst:1313` (a statement about
*partitions*); `foundations/discretization.rst:32` (the LD limiter `w ∈ [½,1]` — a different
object).

### Re-capture checklist (the artefacts, in order)

1. `tests/sn/regression/snapshots/cyl_1g_homogeneous_folded_{4x8,2x4}_dd_n20.npz`,
   `cyl_2g_3reg_folded_4x8_dd_n40.npz` — via `tests/sn/regression/_generate_snapshots.py`.
2. `tests/sn/regression/snapshots/walk_matvec_cyl_2g.npz` — via `_generate_walk_baselines.py`.
3. `tests/sn/verification/analytical/test_phase_c_crosscheck.py:213-215` — the **hardcoded**
   `_SNAPSHOT_KEFFS` cylinder literals (does **not** update itself with the `.npz`).
4. `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py:120-125` — re-measure the
   azimuthal-floor ladder under the retired absorber (the 6.4 acceptance test; the docstring
   is the live single source for those numbers).
5. Confirm `tests/sn/sweep/core/test_ordinate_scan_reset.py::test_pole_resonance_unreachable_on_admitted_family`
   is still green (its own message says a RED means 6.4 owes a live ERR-054 reproducer).

---

## Addendum — MUST-FIX #16 upgraded after checking whether it is executable

`[M]` It **is** executable, and it is a **silent coverage gap**, not a red:

```
tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py:255
    @pytest.mark.parametrize("tau_value", [0.5, 0.6, 0.75, 0.9, 1.0])
    def test_recurrence_at_arbitrary_tau(...)
        """For any τ ∈ [1/2, 1] the M-M weighted DD recurrence holds …"""
```

The L0 kernel gate for the ψ½ recurrence sweeps τ over **exactly the absorber's box**. After
6.4 the production τ on the admitted folded family is `[M] [0.2047, 0.7953]` (§B) — so
**every value this L0 gate tests is at or above the production minimum, and the four
production values below ½ are never exercised.** It stays GREEN while silently ceasing to
cover the operating range. Widen the parametrization to include the fold's floor (e.g.
`0.2, 0.2929, 0.4142, …`) and re-word the docstring's domain to `(0, 1]`, in the same commit
as the retirement — otherwise the retirement quietly narrows L0 coverage of its own kernel.

`[M]` Confirmed clean, no action: `orpheus/sn/sweep/__init__.py:28-33` re-exports only
`IdentityAngularClosure`, `MorelMontryAngularSweep`, `PoleAngularClosureBase`,
`default_angular_closure_class` — **neither τ producer**. The public surface is exactly the
two `__all__` entries at `pole_angular_closure.py:1633-1634` (§E.1).
