# Archivist — Memory Index

Slim index. Behavioral lessons live in `lessons.md` (read FIRST each
dispatch). The mechanical build-gating / cross-ref / venv-worktree /
close-out-arc procedure lives in `AGENT.md` ("Build-Gating & Cross-Ref
Reality", "Close-Out Narrative Arc"). The V&V vocabulary lives in the
`vv-principles` / `algebra-of-record` skills. This index holds only
(1) the lessons pointer, (2) git-true active/doc-debt state, (3) durable
doc-architecture reference. Campaign play-by-play is retired — its
behavioral lesson is in `lessons.md`; its landed milestones are in the
SN theory page's "Development history" section.

## 1. Lessons — a HOT digest over a COLD archive

Same hot/cold split as this index: read the digest always, page the archive on
demand. **Never re-summarize a lesson here or in the digest — each layer points
down, it does not copy up.** Counts are deliberately NOT quoted (a frozen number
rots; `grep -c '^## L-0'` answers it).

- [lessons.md](lessons.md) — **HOT digest, read FIRST every dispatch** (one `Read`
  fits it). Every lesson as one imperative rule + its failure→correction core, in
  9 themes: (1) verify against the LIVE tree · (2) the build is blind, grep is the
  gate · (3) a `:label:` is a V&V edge · (4) retirement & staleness · (5) page
  surgery · (6) doc SHAPE per event class · (7) V&V vocabulary curation · (8)
  code-prose rebalance · (9) gates & tooling. Each entry carries a `→ L-0NN`
  pointer into the archive.
- [lessons_archive.md](lessons_archive.md) — **COLD, load on demand** (~250 KB —
  never open whole). One `## L-0NN` section per lesson: war stories, evidence
  tables, `file:line` detail. Open ONLY the section a pointer sent you to; new
  lessons are appended HERE first, then distilled into the digest.

## 2. Active / doc-debt state — git-true

**No owed Sphinx pass on merged work.** Landed SN milestones live in the "Development
history" changelog at `docs/theory/methods/sn/history.rst` (⚠ NOT the pre-split
`discrete_ordinates.rst`, and orphaned July HTML survives in `_build`, so a stale-ref
grep must `test -f` the SOURCE). Active track = **#231** (§3); main agent commits, I
stage + gate.

**In flight (ONE line each; the evidence is in the lesson, the state is in git):**

- **THE NAMING LAW: an orbit space is named by its STABILISER (#432)** — 2026-09-02, same
  branch, #429 tracker 1.9. Docs **UNCOMMITTED** (mine): 9 hand-edited pages **+1144/−355**
  — `manifolds.rst` (a new `-` chapter *An orbit space is named by its STABILISER*: the
  naming law, `orbit_stabiliser`, the construction invariant, the 5-element symbolic table,
  the three near neighbours, the exact invariance criterion, the before/after readings, the
  24-of-24 stage-0 null; the lattice section rewritten with the COMPUTED axial relations and
  11 edges; the compatibility law re-run; the walk's SHRINKING report; a #432 dev-history
  row; the machine header), `frame.rst` (row 7 ✅ + a NEGATIVE row + the ⛔-tombstoned
  warning; the 5/2 → **6/3** `TruncatedBasis` census), `spherical_harmonics.rst`,
  `discrete_measures.rst`, `spaces.rst`, `angular_quadrature.rst`, `cartesian_multid.rst`,
  `operator_algebra.rst`, `error_catalog.rst` (ERR-080 + ERR-072), + regenerated
  `matrix.rst` (sentinels 582 → **584** = my two new `documented` eq-labels, exactly as
  predicted). `-E -W` **0 → 0**, EXIT=0. ⚠ **ONE dead xref BY INSTRUCTION**
  (`SubgroupOfO3.orbit_stabiliser`, the accepted-but-unlanded accessor) — patched gate 1
  dead/2 sites, nexus `dead_references` 1 of 53, the SAME finding from both. ⚠ **REPORTED,
  code-only** (2: a `test_e4` docstring naming `SO2('x')` over a body that builds `O2("x")`;
  a `symmetry.py` comment calling D_∞h "the group a cylinder actually carries" where the
  registry row spends `Trivial`). ⚠ Mid-task DESIGN DELTA replaced the design I had already
  measured; two brief numbers refuted; 8 present-tense-false claims found outside the brief.
  → [[lessons-L86]]

- **THE FIX: a 1-D rule's frame binds the basis its ORBIT SPACE admits — ERR-080 CLOSED** —
  2026-09-02, same branch, #429's fused commit (0.1b + 0.6 + 2.2 + 3.4 + 3.4b). Docs
  **UNCOMMITTED** (mine): 12 hand-edited pages +1756/−240 — `manifolds.rst` (new `=` chapter
  *What descends*: the isotypic probe, the `Descent`, G0; 4 Key Facts bullets re-tensed, status
  YAML, seams row, changelog), `error_catalog.rst` (ERR-080 **FIXED** banner + before/after
  table + every "still OPEN" clause retensed), `frame.rst` (new `frame-g0-descent-arrow` §,
  7-row admit/refuse table, MomentHead; the flagship DENSE witness retracted), `spaces.rst`
  (`spaces-moment-head` §; the frame-square 3-way table re-measured), `spherical_harmonics.rst`
  (new `=` chapter *the 1-D family*), `angular_quadrature.rst`, + `adjoint`/`operator_algebra`/
  `acceleration`/`indexing_and_layout`/`slab_multigroup`/`cartesian_multid`, and regenerated
  `matrix.rst` (sentinels 579 → **582** = my three new `documented` eq-labels). `-E -W` **0 → 0**,
  EXIT=0; patched-xref `DEAD TARGETS: 0` (positive control read 2); nexus `dead_references`
  0/52. ⚠ **REPORTED, code-only** (4 items, headed by a docstring number off ~60× because it
  names the TABLE's observable while the sentence is about the FLUX). ⚠ Two brief numbers
  refuted, one corpus-wide shape contract found stale at 9 sites, 8 self-inflicted
  nested-markup leaks caught and fixed. → [[lessons-L85]]

- **The angular moment space is READ off the frame, never minted from `L`** — 2026-09-02,
  same branch, tracker **2.5** (`TruncatedBasis` Protocol; both `HarmonicFrame` doors; seven
  re-mint sites; `truncated` on the head). Docs **UNCOMMITTED** (mine): `frame.rst` +427 (a new
  `-` section + 6 `~` subsections + 1 Key Facts bullet; ONE new `documented` eq-label
  `moment-space-read-off-the-frame`), `operator_algebra.rst` +24 (a tombstoned `implements`
  body), `manifolds.rst` +23 (status YAML + the lower-bound remedy), `spaces.rst` +11,
  `cartesian_multid.rst` +16, `error_catalog.rst` ERR-080 +49 (**stays OPEN**), regenerated
  `matrix.rst` (sentinels 578 → **579**, exactly as predicted; the +37 tests are the CODE
  side's). `-E -W` **0 → 0**, EXIT=0; xref gate `DEAD TARGETS: 0`; nexus `dead_references`
  0/52. ⚠ **REPORTED, code-only** (4 items, headed by the `harmonic_moment_source_sink`
  docstring twin the step's own sibling correction skipped). ⚠ Three brief numbers did not
  reproduce (its census command, "12 of 12", "96–161 %") and its two candidate host pages were
  both wrong — detail in the lesson. → [[lessons-L84]]

- **A quotient carries TWO coordinate systems — the chart and the section** — 2026-08-31/09-01,
  branch `fix/angular-phantom-support` (#429). Docs **UNCOMMITTED** (mine; main agent commits):
  `manifolds.rst` +1152/−85 + the regenerated `matrix.rst`; code side landed `b55bba56`
  **mid-session**. All gates green and every generated artefact moved exactly as predicted.
  ⚠ **REPORTED, code-only** (3 items: `__all__`, a symbol collision, a hand-typed fixture) —
  detail, the four findings and the numbers are in the lesson. → [[lessons-L80]]

- **A catalogue entry carries its own ARROW and the measure that arrow pushes forward** — 2026-09-02,
  same branch, tracker **3.1** (`Quotient.orbit_coordinates` + the derived `quotient_map`;
  `Quotient.reference`; the registry twin `AngularSymmetry.reference` collapsed). Docs
  **UNCOMMITTED** (mine): `manifolds.rst` +1226, `discrete_measures.rst` +58,
  `error_catalog.rst` ERR-080 +55 (**stays OPEN**), regenerated `matrix.rst` (sentinels
  577 → **578**; the test-count moves are the CODE side's). One new `documented` eq-label
  `manifold-quotient-pushforward`. `-E -W` **0 → 0**, EXIT=0; xref gate 0 dead (stock AND
  patched, split control run); nexus `dead_references` 0/52. ⚠ **REPORTED, code-only** (2
  items: a future-tense docstring at `manifold.py:968`; `reference.support` is the CHART's
  space and nothing gates the pair). ⚠ Three brief numbers corrected (a "pickle round-trip
  equal" that is 1 of 7; `.reference` reads 9 → 10; a `match Quotient` site attributed to
  `barycentre`, which uses `isinstance`) and two widened 5 → 7. → [[lessons-L83]]

- **A map carries its two point sets, so a codomain cannot be forged** — 2026-09-02, same branch,
  tracker **2.3** (`ManifoldMap`; `archimedes`, the orbit retraction, `barycentre`). Docs
  **UNCOMMITTED** (mine): `manifolds.rst` +883/−42 (a new `=` section + 6 `-` subsections, one new
  `documented` eq-label `manifold-map-functoriality`), `discrete_measures.rst` +46, `error_catalog.rst`
  ERR-080 +68 (**stays OPEN**), `spherical_harmonics.rst` +32, regenerated `matrix.rst`
  (10595 → **10616**; sentinels 576 → **577**; `numerics/test_manifold` 56 → **70**,
  `test_rules_product` 38 → **45**). `-E -W` **0 → 0**, EXIT=0. ⚠ **REPORTED, code-only** (4 items:
  a docstring promising a zero-fallback where the body raises, a gate docstring saying "two"
  strict-xfails where three ship, a `Ball` overstatement, an over-indented comment block).
  ⚠ Two brief claims REFUTED and one instrument found broken — the detail is in the lesson.
  → [[lessons-L82]]

- **A basis declares its symmetry by naming what it EATS** — 2026-09-01, same branch, tracker
  **2.1b** (`Basis.invariance_group`, DERIVED + `@final`). Docs **UNCOMMITTED** (mine):
  `manifolds.rst` +591 (a new `=` section + 4 `-` subsections), `discrete_measures.rst` +90 (a
  new HAS/SPENT subsection — the brief's target section did **not exist**), `spaces.rst` +18,
  `error_catalog.rst` +55 (ERR-080 **stays OPEN**), regenerated `matrix.rst` (10584 → **10595**,
  `numerics/test_basis_domain` 13 → **24**, sentinels **576 → 576**). `-E -W` **0 → 0**, EXIT=0.
  ⚠ **REPORTED, code-only**: `orpheus/numerics/measure.py:417` cites `:meth:`reorder``, which
  does not exist. ⚠ Four present-tense-false corpus claims found and tombstoned in place, none of
  them in the brief. → [[lessons-L81]]

- **ERR-026 history block: 29 roles → 13, 15 dead → 0** — branch
  `docs/err026-history-is-not-a-crossref`, 2026-08-18. ⛔ still open, re-confirmed unlanded
  2026-08-24: the `head_role` one-liner (blindness is ROLE-scoped, not `.rst`-scoped) and
  **40 of 100** stale raw file paths in the catalogue. → [[lessons-L62]]

**⏹ MERGED — collapse to one line each; the durable record is the lesson + the tree.**
⚠ This list has frozen "in flight" on landed work **SIX** times; reconcile with `git` first —
`git merge-base --is-ancestor <hash> HEAD`, and `git branch --list <branch>` (a vanished branch
means merged). The 2026-09-01 pass found FOUR stale at once.

- 2026-08-31: **manifolds page** the point-set layer gets a page (`fba4205a`) [[lessons-L79]] ·
  **CS4c step 4** the fission channel becomes two bindings (`4e46dbb9`) [[lessons-L77]] ·
  **(n,2n)** physics-vs-model isotropy (`6906f2a2`) [[lessons-L78]].
- 2026-08-30/29: **P7** the metric becomes an OBJECT (`2ef04dbb`) [[lessons-L76]] ·
  **P4-remainder** the producer binds the axis (`cd176f69`) [[lessons-L75]].
- 2026-08-29: **CS5** axis-generator doctrine (`4e7b8977`/`b0bfc06c`/`cb3cd15b` + docs
  `f8c69117`; all 3 reported code gaps repaired at `cb3cd15b`) [[lessons-L74]] ·
  **§5b P0** carrying-prose sweep (`8d093334`) [[lessons-L73]] · **P4.9b** operator poses
  (`9c3eb60a`) [[lessons-L72]] · **P4.9a** closure owns its march (`ca852c44`/`7a0f434c`)
  [[lessons-L71]].
- 2026-08-28 and earlier: α-dome citation retraction + `alpha-dome-recursion` rename
  [[lessons-L70]] · τ-arity + LD-curvilinear-Padé [[lessons-L69]] · campaign-1 history rows
  `68d265ef` [[lessons-L68]] · Campaign 1 `55bb47b9` [[lessons-L68]]…[[lessons-L63]] ·
  MD→corpus catalogue port `a79f57aa` [[lessons-L61]] · nexus #82 declarations `144cdf51`
  [[lessons-L60]] [[lessons-L59]] · nexus graph-path retirement [[lessons-L58]] · #344
  loss-kernel-gauge [[lessons-L57]] · Q5.6.4 τ/partition carve [[lessons-L54]] · Boundary
  B3.0–B3.2 [[lessons-L42]] · DSA #2 close-out [[lessons-L39]].
- Their CODE-side reports are GitHub's to track, not this file's; the corpus-wide
  RST-nested-markup finding lives on **#379**.

## 3. Durable reference (reusable doc-architecture)

Each entry is a ONE-LINE pointer; the full recipe lives in the linked `feedback_*.md`.

- **Landed-milestone record:** `docs/theory/methods/sn/history.rst`. POINT here instead of
  re-listing campaigns. (This line named the pre-split `discrete_ordinates.rst` until 2026-08-18,
  contradicting §2 four lines up — an index can go stale against ITSELF.)
- **Ontology-overturn rewrite** (a page whose THESIS was refuted): the recipe is
  [[lessons-L63]] — argument-unit not symbol-unit, the 4-way eq-label fate rubric, the
  unlabelled-history-equation trick, the two-sided illegal-states rule. Instance:
  `field_algebra.rst` affine → cone (CS3).
- [canonical-convention-page](feedback_canonical_convention_page.md) — 13-section anatomy
  for a multi-PR migration's canonical theory page + keep/flip rubric (`index_convention.rst`).
- [canonical-axis-convention SSOT section](feedback_canonical_axis_convention_ssot_section.md)
  — SSOT section for an axis-flip enforced at a data-ingest boundary (`cross_section_data.rst`).
- [double-category architecture insight](feedback_double_category_architecture_insight.md) —
  documenting a categorical framing of a SHIPPED type system; impossibility as an
  obstruction table. Instance: (Rep×Role) carrier grid (#268/#261).
- [orientation-axis two-frames doc](feedback_orientation_axis_two_frames_doc.md) — 2×2-face
  operator unification with ORIENTATION as the coherence axis (#280 P2.5e).
- [carrier-grid-typed-seam-layering](feedback_carrier_grid_typed_seam_layering.md) — NxM typed
  grid + seam one layer up; completing one path silently stales a sibling claim (Frame P4).
- [capstone-architecture-page](feedback_capstone_architecture_page.md) — a NEW page for the
  LAYER above per-method pages (cross-ref, don't duplicate). `loss_representations.rst`.
- [capstone-completion-status-reaudit](feedback_capstone_completion_status_reaudit.md) — the
  COMPLETION phase: re-audit ship-state claims; document an unbuilt sibling as a SEAM (P7).
- [capstone-root-cause-ruling](feedback_capstone_root_cause_ruling.md) — retrofitting the
  structural WHY (a theorem) behind a split the docs only ASSERTED (#268 `frame.rst`).
- [operator-classes→frame-faces re-homing](feedback_operator_classes_to_frame_faces_rehoming.md)
  — sweep when standalone operator classes retire into two FACES of one frame (#268 P1).
- [operator-reification/retype doc pattern](feedback_operator_reification_retype_doc_pattern.md)
  — reifying a duck-typed operator; block coisometry `= 4π·I`, never `= I` (#226 step 2).
- [named-family-member theory section](feedback_named_family_member_theory_section.md) — a NEW
  § for a named member of an invariant-keyed operator family (#226 step 4 GreenOperator).
- [step-5b first-consumer close-the-loop](feedback_issue_138_step5b_first_consumer_closeloop.md)
  — wiring the FIRST consumer of a verified-but-unwired type; → vv Mode 12 (#226 step 5b).
- [consumption-mode + capability-axis](feedback_consumption_mode_and_capability_axis.md) — a
  NEW consumption mode on an operator algebra (solve/apply/ASSEMBLE) (#272/#284/#282).
- [algebra-of-record stub→narrative](feedback_stub_to_rich_narrative_expansion.md) —
  SymPy-module-as-canonical-source; stub/expand separation (also lessons L5).
- [solver-replacement campaign close-out](feedback_solver_replacement_campaign_closeout.md) —
  a legacy island solver replaced by the operator-algebra family; LIVE/MOOT split (#290 P8).
- [type-confinement docstring sync](feedback_type_confinement_docstring_sync.md) — code-final
  sync when a carve confines a subtype to one role (P4.5 W-C).
- [Petrov-Galerkin homogenization reframe](feedback_petrov_galerkin_homogenization_reframe.md)
  — THE LIVE recipe: flux-weighting is a TEST weight, not a measure (#268 P3). Supersedes
  [Galerkin-natural-metric](feedback_galerkin_natural_metric_reframe.md) (why-it-was-tried only).
- [domain-op + L2-promotion + asymmetry-law](feedback_domain_op_l2_promotion_asymmetry_law.md)
  — section shape for a domain OPERATION born from an L2 promotion (#267).
- [orbit-space terminology sweep](feedback_orbit_space_terminology.md) — add-aside-then-bridge-then-sweep for a precise math term.
- [auto-generated tables](feedback_autogen_tables.md) — registry-as-SSOT: metadata fn +
  generator + `builder-inited` hook (also lessons L8).
- [audit-then-edit partitions](feedback_audit_partition.md) — the KEEP/RELOCATE/TRIM/REMOVE
  partition table for a read-only doc-cleanup audit.
- [cross-solver unified-law doc architecture](feedback_cross_solver_unified_law_doc_architecture.md)
  — ONE law spanning N solver families: canonical derivation + short sibling spellings (#259/#291).
- **Doc-architecture redesign (#231, OPEN):** the standing target for any "modernize a theory
  page" task — template, machine header, prose rebalancing, V&V slices, bibtex (spec in the
  issue). Phase 1a–1c DONE; **Phase 2 code-prose rebalancing ACTIVE** — P2-A/B/C/D/G done
  (maps in `.claude/plans/phase2_code_prose/`). Five file-classes calibrated:
  teaching-operator = aggressive TWIN-cut; machinery/driver/mesh = small, COMMENTS dominate;
  ABC = leanest; contract-heavy-operator = small. Main agent commits; I stage + gate.
