# SN monolith split — the Phase C execution catalog

> **⚠ RE-SCOPED (user, 2026-07-16).** The SN book is being re-architected around the
> **invariant-first broadening progression** — see **`sn_book_architecture.md`** (CONFIRMED),
> which supersedes the *topic-based* Tier-2/Tier-3 plan below. What **remains valid**: the four
> **Tier-1 shared chapters** already extracted (`angular_quadrature`, `boundary_conditions`,
> `verification`, `history` — orthogonal to the progression) and the **per-chapter recipe +
> lessons** (L-024 / END-boundary / dual-namespace / the `-E` build gate). What is **superseded**:
> the topic chapters (`algebra`/`scattering`/`discrete_balance`/`spatial_closures`/`sweep_*`/
> `curvilinear` as separate topic pages) — that CORE content now gathers into the progression
> chapters (Cartesian slab 1g → energy → multi-D; then curvilinear) via a reorganization +
> connective rewrite, and the build order starts with **`foundations/discretization.rst`**.



**Subordinate to** `documentation_corpus_architecture.md` §5 (the SN mini-book chapter
table), §7 Phase C, §7.2 (Phase B landed), §8 (rulings). This is the **working driver**:
the span→chapter map + per-chapter status. Read the corpus plan for the *why*; read this for
the *what-next*.

**Target:** `docs/theory/methods/sn/index.rst` (20,571 ln at Phase C start) → ~15 chapter
files + a thin router `index.rst`. `loss_representation.rst` is already split (Phase 1a).

---

## Rulings baked in (user, 2026-07-16 — Phase C kickoff)

1. **No filename prefix.** Order lives in the toctree; reorder never renames a file (labels
   are path-immune, but raw filenames are not). Names are greppable: `algebra.rst`,
   `curvilinear.rst`, `angular_quadrature.rst`.
2. **Defer the two cross-part sections.** Phase C stays *inside* `methods/sn/`. "The Transport
   Equation" (→ `foundations/transport_equation.rst`) and "Verification" (→ the `verification/`
   part, Task #10) become **temporary SN chapters** now and relocate later via cheap
   path-immune `git mv`. Smallest blast radius, fully reversible.
3. **Phase C commits are pure structural moves.** Extract the span *verbatim* into the chapter
   file, wire the toctree, fix only L35 prose breaks. The §6 chapter-template polish (synopsis
   retrieval anchors, diagnostic tables, code maps) is **Phase D/G**, not here — this keeps
   each commit a clean cut+paste with a reviewable diff (the Phase B pattern).

---

## The build gate — `sphinx-build -E -W` (NOT plain `-W`)

**Moving a label between files produces a SPURIOUS `duplicate label` warning under an
incremental build.** Sphinx's saved environment (`docs/_build/doctrees`) still has the label
registered to the source page; the incremental re-read adds the destination's copy before
fully clearing the stale source registration → transient "duplicate". A plain `-W` rerun then
flips to exit 0 once the env settles — so **neither incremental result is trustworthy on a
label move.** Every Phase C chapter moves labels ⇒ the authoritative gate is:

```
.venv/bin/sphinx-build -E -W --keep-going -b html docs docs/_build/html
```

`-E` forces a full environment rebuild (no saved cache), so the label registry is consistent.
(Do NOT `rm -rf docs/_build` — it trips the host permission gate; `-E` achieves the same
without deletion.) Confirmed on ch3: `--keep-going` showed the phantom duplicate, `grep -c`
proved the label single-homed, `-E -W` came back clean. [→ fold into lessons as the Phase-C
build-gate note.]

---

## The per-chapter recipe (proven on ch3)

For each chapter, in order:

1. **Re-derive the span.** Line numbers shift with every cut — grep the **section title** (or
   its anchor label, below) to find the current start; find the next same-or-higher `=`-header
   for the end. (`rst_skeleton.py <file>` in the job tmp re-dumps the tree.) **END-boundary
   rule (learned on verification):** if the NEXT section's `.. _anchor:` sits ABOVE its header
   (the common style here), the cut END is `anchor_line − 1`, NOT `header_line − 1` — that
   anchor belongs to the STAYING section; stopping at `header − 1` drags it into the new file.
2. **L35 three-way grep** (the load-bearing move check — `-W` gates the LINK, never the PROSE):
   - inbound refs to every label in the span (`grep -rn '<label>' docs/ orpheus/`) — path-immune,
     they survive, but note directional/qualifier prose *in the referring page*;
   - the span's own **directional prose** (`above|below|earlier|later|preceding|this section`)
     that becomes false once it is its own page;
   - **bystander page-qualifiers** ("`:ref:`X` in :doc:`methods/sn/index`") anywhere tree-wide
     that named the old home.
3. **L34 file-level sweep** (splits create NEW files ⇒ raw-path hazard returns): grep the new
   filename's **last segment** and any raw `methods/sn/…` path strings in prose/docstrings/code.
4. **Cut:** `split_chapter.py SRC START END DEST` (job tmp; byte-exact, prints the new seam).
5. **Wire the toctree** in `index.rst` — insert the new entry in **§5 reading order**.
6. **Fix** any L35 prose breaks + copy any document-local substitution (`|times|`-class) the
   span used but whose `.. |x| ::` definition stays behind.
7. **Gate:** `-E -W` clean (above). Confirm every moved label single-homed **on its anchor
   DEFINITION**: `grep -c '^\.\. _<label>:' <source>` = 0 **(NOT raw `grep -c '<label>'`)** +
   tree-wide `grep -rn '^\.\. _<label>:'` = 1. A router page legitimately KEEPS bare
   `:ref:` back-refs to labels it exported — those are path-immune and STAY; recut ONLY on a
   surviving `.. _label:` DEFINITION or a phantom `-E duplicate label`. (Learned on ch12 BC,
   whose anchor is back-referenced from `index.rst` itself — raw count 2, anchor-def 0.)
   **Dual-namespace (learned on verification):** a name can be BOTH a `.. _X:` section anchor
   (std domain, `:ref:`) AND a `:label: X` equation label (math domain, `:eq:`) — different
   Sphinx domains, so they coexist with no dup-label warning. If a split moves only ONE, single-
   home on the CORRECT namespace: a moving eq-label ⇒ `grep -c ':label: X'` (=0 source, =1 tree);
   a moving section anchor ⇒ `grep -c '^\.\. _X:'`. The same-named member that STAYS reads as
   count 1 in the other namespace and is NOT a recut trigger. The `-E` build is the collision
   oracle (a clean baseline proves the two domains don't collide).
8. **Archivist reports** the diff + grep results + build; **main agent reviews diff-vs-catalog,
   re-runs `-E -W`, commits** with the correct trailer (archivist does NOT commit — Phase B model).

**Delegation (RULED by user, 2026-07-16):** ch3 done by main agent (pattern-prover). Tier-1/2
clean extractions → **archivist**, one chapter per dispatch (resume the same instance across
chapters to keep context; main agent reviews diff-vs-catalog + commits with the correct
trailer — archivist does NOT commit). Tier-3 (dedup-heavy: `algebra`/`iteration`/`solver`)
stays **main-agent + elegance-enforcer** review (it collapses duplicated content, not a move).
**Content ruling: PURE structural moves** — verbatim cut, NO synopsis/router prose added
during the split (that is Phase D). Each chapter = its own reviewable pure-move commit.

---

## The chapter catalog (keyed on stable section title + anchor; line #s advisory only)

`§5#` = chapter number in the corpus plan §5 table. **Tier** = execution order (1 = clean
leaf first, 3 = dedup-heavy last).

| §5# | Chapter file | Source `=`-section(s) — grep this title | Anchor | Tier | Status |
|---|---|---|---|---|---|
| 3 | `angular_quadrature.rst` | "Angular Quadratures" | `quadrature-types` | 1 | ✅ DONE |
| 12 | `boundary_conditions.rst` | "Boundary Conditions" | `boundary-conditions` | 1 | ✅ DONE |
| — | ~~`transport_equation.rst`~~ | "The Transport Equation" — **STAYS in index as SN intro** (user directive; NOT a chapter, NOT → foundations, §4 superseded) | (none) | — | N/A |
| — | `verification.rst` | "Verification" + "Numerical Sensitivities" (temp → V&V part) | (none above h1) | 1 | ✅ DONE |
| — | `history.rst` | "Development history" | `sn-development-history` | 1 | ✅ DONE |
| 5 | `discrete_balance.rst` | "The Discrete Balance Equation" (balance/α/geom-factor/streaming-eq/Morel-Montry dip — up to the WDD weights) | (none) | 2 | pending |
| 6 | `spatial_closures.rst` | rest of "The Discrete Balance Equation" (WDD & Morel-Montry weights, closure substitution) + "Cell update strategies" | `cell-update-strategies` | 2 | pending |
| 9 | `sweep_1d.rst` | "Sweep Algorithm" intro + Cartesian-1D cumprod + affine outflow | `sweep-algorithm` | 2 | pending |
| 7 | `linear_discontinuous.rst` | "Sweep Algorithm" → UBLD tensor-product bilinear cell system | — | 2 | pending |
| 10 | `sweep_multid.rst` | "Sweep Algorithm" → Cartesian-2D wavefront + octant DAG | — | 2 | pending |
| 11 | `curvilinear.rst` | "Sweep Algorithm" → curvilinear sequential sweep + unified dispatch + Carlson D/F + ERR-058 + route-(a); PLUS "The curvilinear anisotropic-MMS floor" | `sn-curvilinear-aniso-norm-reconciliation` | 2 | pending |
| 4 | ~~`scattering_and_moments.rst`~~ | **SUPERSEDED by progression** — "Scattering" + "Scattering and fission as LinearOperators" gathered by `slab_multigroup` @ `b7166ed6`; "Angular windowing" STAYS for `cartesian_multid` | `sn-angular-windowing` (still in index) | — | absorbed |
| 8 | (append to `loss_representation.rst`) | "InvertibleOperator: the streaming-collision operator algebra" | `sn-streaming-operator` | 3 | pending |
| 13 | ~~`iteration.rst`~~ | **SUPERSEDED by progression** — Krylov H1 → `slab_one_group` @ `2c60d6a5`; KEigenvalue/choosing/operand/FEAST/cross-solver → `slab_multigroup` @ `b7166ed6`; the boundary G-S schedule H2 STAYS under `sn-iteration-primitives` for `cartesian_multid` | `sn-iteration-primitives` (still in index) | — | absorbed |
| 14 | `solver.rst` | "SNSolver as an operator-algebra coordinator" + "Consuming the frame in SN" | `sn-solver-operator-algebra-coordinator`, `sn-consuming-the-frame` | 3 | pending |
| 2 | `algebra.rst` ★★ | THE SPINE — synthesized: dedup the `(A−S−F)ψ=q` spine (restated 3× — remaining sites: the SNSolver-coordinator + Angular-windowing H1s). ~~collapse `choosing-inverse-realisation` vs the inverse-family labels~~ (DONE @ `b7166ed6` — compressed into `slab_multigroup`, foundations owns the catalog); ~~`cross-solver-eigenvalue-consumers` vs `matrix-inverse-consumers`~~ (RULED no-collision — different concepts) | (new) | 3 | partially absorbed |
| 15 | `adjoint.rst` | "The scattering adjoint, free from the harmonic frame" (subsection of InvertibleOperator) + new S†/daggered-posing/φ* | (new) | 3 | pending |
| 1 | `placement.rst` | NEW writing (why SN vs CP/MoC/P_N/diffusion, M7) | (new) | — | Phase D/defer |

**Stays in `index.rst` (the router):** title/intro (`theory-discrete-ordinates`), "Synopsis"
(`sn-synopsis`), "Architecture", the "Gotchas" *index* (`sn-gotchas` — content distributes
inline to chapters, index keeps the list per §5), "References", the `.. _sn-chapters:` toctree.

**The Block-1 dedup** (the reason Block 1 folded into Phase C, §8 Phase-B ruling) happens in
the **Tier-3** chapters `algebra`/`iteration`/`solver` — that is where the thrice-restated
spine and the duplicated inverse-realisation labels collapse to one home. Not a move; a merge.

---

## User content directives (2026-07-16) — the method-chapter content contract

Governs how the **algebra / discrete_balance / loss_representation / iteration** chapters are
built (Tier-3 + the two Tier-2 discretization chapters). NOT a Phase-C structural concern —
the content mostly already exists in the monolith; this is the ORGANIZING contract and a
gap-list for Phase D/G. (User, answering the transport_equation question — a directive, not a
placement.)

1. **General → method-specific derivation is per-method, shown in full.** The abstract
   transport equation lives at the ROUTING level (the `methods/index` router + the Phase-H
   root); each method chapter (here S_N) then SHOWS how that general equation becomes THAT
   method's equation, with detailed discretization steps. ⟹ SN's "Transport Equation" section
   (the per-geometry continuous forms) is **SN-book intro and STAYS in the SN book** (feeds the
   `algebra` chapter); it is **NOT minted as a chapter and NOT extracted to
   `foundations/transport_equation.rst`** — that §4 plan is **superseded** for the SN geometry
   forms (the shared root is the more abstract method-invariant operator, Phase H).
2. **Show the loss-representation MATRIX SHAPE.** Each method states the shape of its loss
   representation (SN's `(L+C)` is lower-triangular; others dense; etc.). **The matrix
   representation is the MINIMAL / baseline representation.**
3. **Every non-matrix loss representation must show its strategy to AVOID materializing the
   matrix, and LINK the shape to the strategy** — lower-triangular ⟹ a sweep; a Krylov
   subspace projection; etc.
4. **If a strategy's representation is reached by an operator-algebra operation that POSES the
   problem into the strategy-enabling shape, the operations FROM the general representation TO
   that strategic loss representation must be SHOWN** (the algebra that transforms `A` into the
   sweep-able / Krylov-able form). This is the core job of the `loss_representation` +
   `iteration` + `algebra` chapters. [[feedback-build-the-machinery-operator-algebra]]

**Verification scope (RULED):** `verification.rst` = "Verification" + "Numerical
Sensitivities"; "Consuming the frame in SN" → the `solver` chapter.

---

## Status log

- **ch3 `angular_quadrature.rst`** — ✅ **@ `5a61fdb5`.** Cut 366–500 (135 ln), 1 label
  `quadrature-types`, 1 inbound ref (MoC, path-immune ✓), 0 directional prose, 0 subs.
  L35+L34 clean; `-E -W` 0 warnings. index.rst 20,571 → 20,436. Clean −135/+1. Main-agent
  pattern-prover.
- **ch12 `boundary_conditions.rst`** — ✅ (archivist-executed, first delegated chapter). Cut
  12497–12904 (408 ln), 3 anchors (`boundary-conditions`, `bc-tensor-decompositions`,
  `bc-sn-resolution-table`) + 3 eq-labels. **1 genuine L35 falsehood found + fixed**: a
  bystander `:doc:` page-qualifier in `foundations/boundary_conditions.rst` naming the old
  home → repointed to the new chapter. `-E -W` 0 warnings. index.rst 20,437 → 20,030.
  Sharpening L-024: single-home on the anchor DEFINITION, not raw mentions (BC's anchor is
  back-referenced from `index.rst` — those bare `:ref:`s stay, path-immune). Delegated model
  validated. NB: `methods/sn/boundary_conditions.rst` deliberately parallels
  `foundations/boundary_conditions.rst` (SN realization ↔ general theory) — disambiguate greps
  by full path.
- **`history.rst`** — ✅ @ `80228085` (archivist-executed). Cut 18894–19912 (1,019 ln), 1
  anchor `sn-development-history`, 0 eq-labels, **0 inbound refs, 0 prose fixes** — the
  cleanest cut. Confirmed L-024's `:doc:`-qualifier discriminator: the 9 `:doc:` qualifiers
  in the span all name THIRD (foundations) pages, not the old home, so none went false.
  `-E -W` 0 warnings. index.rst 20,030 → 19,012. Toctree now: angular_quadrature,
  loss_representation, boundary_conditions, history.

- **`verification.rst`** — ✅ (archivist-executed). Cut 16077–18701 (2,625 ln = "Verification"
  + "Numerical Sensitivities"; "Consuming the frame" correctly STAYED for `solver`). 41 labels
  (10 section anchors + 31 eq-labels), **0 prose fixes** (all inbound are bare cross-doc
  `:ref:`/`:eq:` or global label-name strings in the V&V registry/matrix/tests). Two new
  subtleties, both handled + folded into recipe steps 1 & 7: the **END-boundary rule**
  (stopped at 18701 so the next section's `sn-consuming-the-frame` anchor stayed) and the
  **dual-namespace** pair (`sn-mms-{spherical,cylindrical}-aniso-spatial-convergence` = section
  anchor in the sweep area (stays) + eq-label in Verification (moves) — single-homed per
  domain). `-E -W` 0 warnings. index.rst 19,012 → 16,388.

**✅ TIER-1 COMPLETE — 4 chapters (135 + 408 + 1,019 + 2,625 = 4,187 ln), index.rst 20,571 →
16,388 (−20.3%).** Delegated archivist model validated across four scales incl. a 41-label
2,625-ln cut; zero build warnings throughout; recipe hardened with L-024 (+ END-boundary +
dual-namespace). transport_equation ruled to STAY as intro. ~~NEXT = Tier-2~~ **SUPERSEDED by
the progression (`sn_book_architecture.md`)** — the three old Tier-2 boundary decisions carry
over as progression-time decisions: (a → the multi-D chapter), (b) the ~1,356-ln 2-D LD stress
MMS → V&V-vs-`cartesian_multid` home, (c) the sweep_multid/curvilinear interleave (the
"unified sweep dispatch" section).

---

## Phase C progression — slab_one_group disassembly (RATIFIED, user 2026-07-20)

**Rulings:** (1) **Main agent authors** (pattern-prover, mirroring ch3→Tier-1; delegation
revisited after this chapter proves the pattern); elegance-enforcer reviews the draft; user
steers at outline (done) + commit. (2) **Krylov H1 gathers into §7** (the Cartesian core — the
§3.2 step-4 "alternative strategy encodings" completes in the base chapter; its curvilinear
paragraphs stay behind for Part B). (3) One commit per chapter, gathered spans CUT in the same
commit (no twin homes, atomically); interleaved campaign material stays until its chapter.
(4) Names confirmed by use: `slab_one_group.rst`, then `slab_multigroup.rst`,
`cartesian_multid.rst`, `curvilinear_one_group.rst`, `curvilinear_multigroup.rst`.

**The gather map** (monolith line #s advisory, re-derive by title per recipe step 1):

| § | slab_one_group section | Gathers (span) | Cross-links (no re-derivation) |
|---|---|---|---|
| 1 | Synopsis + Key Facts | (new) | — |
| 2 | The posing (slab eqn + invariant) | "The Transport Equation" intro + "Cartesian 1D (Slab)" (269–287) | discretization §1–2; path_integral skeleton |
| 3 | The discrete balance → per-cell update | "Cartesian 1D Balance Equation" (376–405, eq `dd-cartesian-1d`) + slab parts of "The DD recurrence"/"Diamond Difference" (2205–2350) | discretization §3–4 (Step/DD/LD live THERE) |
| 4 | The within-group operator (A = L+C−S; WHY invertible = lower-triangular; matrix = rawest) | "InvertibleOperator" Cartesian core (12967–13255; the historical two-closures + curvilinear H2s STAY) | loss_representation; discretization §7 |
| 5 | The sweep = forward substitution | "Cartesian 1D: Cumprod" (2388–2444) + "Generic affine outflow" (2445–2559) | methods/sn/boundary_conditions |
| 6 | P₀ scattering (S as projection) | "P₀ Isotropic Scattering" (12538–12553) | operator_algebra intrinsic types |
| 7 | Source iteration + Krylov alternative | SI/splitting/ρ=c slab core (14564–14751 + the c-mode-fold impossibility 15006–15080) + Krylov H1 Cartesian core (12416–12496) | iteration leftovers stay for later chapters |
| 8 | Verification hooks + gotchas + code map | (new, thin) | verification.rst (analytical-order ≠ MMS ruling) |

**Leave-behind (stays in monolith):** everything curvilinear (incl. "Substituting the WDD
Closure" 1907–1941 = the α/ΔA/τ update, eq `dd-solve`), multi-D, P_N/multigroup scattering,
the τ-campaign H3s (838–1495), LD-moment S3/D5b blocks, Jacobi-vs-G-S boundary recovery +
reified splitting + ERR-056 (14752–15005, 15081+), solver-coordinator, frame, Gotchas index.

### ✅ slab_one_group.rst — DONE @ `2c60d6a5` (2026-07-20)

Authored ~735 ln; index.rst 16,388 → 15,603 (−828/+43). 14 labels traveled single-homed;
9 labeled eqs character-identical (elegance-enforcer verified); matrix.rst unchanged.
**Adjudications beyond the map** (each verified against code): (1) the DiamondDifference
class-doc (2205–2350) STAYED intact — the chapter cross-links its slab branch (gutting the
3-branch dispatch doc would break D5); (2) the monolith's blend/Péclet passage RETIRED as
superseded — its thin-limit claim (k→∞, w→0) contradicts `_ubld.d1_closed_form`
(k∈[0,1), w∈[½,1]); spectrum single-sourced at `discretization-ld-blend`; (3) the Krylov
H1's "Consistency" DD-vs-FD story + "Round 2 deviation" RETIRED (FD family deleted;
matvec ≡ sweep one walk, #206 Phase C) — chapter states the current truth; (4) "Why not
extract apply_transpose analytically?" RETIRED (pro-dense-probe argument in present tense;
the analytic reverse walk shipped, dense path retired in D-K); (5) "Two Inner Solvers"
stale two-closure claims rewritten to one-walk truth + its broken `Krylov inner solver`_
title-ref repointed; (6) Pautz2002 citation repointed at the staying ERR-056 rule; (7) new
H2 "The boundary Gauss-Seidel schedule (multi-D)" parents the staying G-S H3s; (8) the
(A−S−F)/A=L+C local convention harmonized to the honest algebra in ALL chapter + inserted
text (monolith's remaining old-convention sections harmonize when their chapters come).
L35 breaks found+fixed: 3 (a "below", the title-ref, the Pautz "cited at" semantics).

## Phase C progression — slab_multigroup disassembly (brief queued + executed 2026-07-20)

**Broadens ENERGY** (arch §3.4 row 2): group-to-group + P_N anisotropy (§3.1 `S=R∘Λ∘M`) ·
fission/χ · k-eigenvalue · power iteration. Same ratified per-chapter model: **disassembly →
user steers → main-agent authors → elegance-enforcer review (its L-013 recipe) → ONE atomic
commit (chapter + cuts) → catalog checkpoint.** Recipe = "The per-chapter recipe" above
(re-derive spans by TITLE — the monolith is now 15,603 ln, all pre-cut line #s stale).

**Gather candidates** (grep these stable titles in `methods/sn/index.rst`; adjudicate each
during the disassembly):

- "Multi-Group Extension" (H2 under The Transport Equation) — the multigroup continuous form.
- "Scattering" H1: "Matrix Convention" (SigS[g_from,g_to], the transpose convention,
  `scattering-matrix-convention` ref), "P_N Anisotropic Scattering" (anchor `pn-scattering`),
  "(n,2n) Reactions", "Normalization Chain".
- "Scattering and fission as LinearOperators" H1 (anchor `sn-scattering-fission-operators`):
  apply-only rationale, Pℓ Galerkin projection on Y_ℓ^m, (n,2n) doubling.
- Iteration-primitives rump: "KEigenvalue: outer power iteration", "Choosing the A⁻¹
  realisation (the inverse-operator family)", "Operand requirements", "Forward hook: FEAST",
  "Cross-solver consumers of power_iteration".
- The chapter also OWNS the eigenvalue posing `(L+C−S−B)ψ = (1/k)Fψ` + `K = A⁻¹F`
  (cross-link `foundations/operator_algebra` + `operator_inverse_family`; #226 taxonomy).

**Expected adjudications** (the slab_one_group pattern predicts these):

1. **Twin risk — "Choosing the A⁻¹ realisation" vs `foundations/operator_inverse_family.rst`**
   (a reframe deep-dive): the monolith section predates the reframe — expect overlap; the
   inverse-FAMILY taxonomy lives on foundations, the SN-specific consumption in the chapter.
   Same disposition class as the Péclet retirement (keystone wins; chapter cross-links).
2. **Harmonization**: the KEigenvalue/choosing sections use the OLD local `A = L+C(−S)`
   convention — chapter text harmonizes to the honest algebra (A = L+C−S−B; sweep=(L+C)⁻¹;
   `K = A⁻¹F` per the taxonomy); retained monolith sections keep theirs until their chapter.
3. **Old Tier-3 note**: `cross-solver-eigenvalue-consumers` vs `matrix-inverse-consumers`
   label near-duplication (see the superseded ch2 row above) — check while gathering the
   cross-solver-consumers section.
4. **Angular windowing H1 is NOT multigroup** — it is the 2-D SI-iterate machinery; leave
   for `cartesian_multid`. The scattering-adjoint H2 stays for the adjoint chapter.
5. **P₀ back-link**: slab_one_group §"P₀ scattering" carries the forward note "the
   coefficient becomes the in-scatter sum … next chapter" — slab_multigroup must land that
   promise (and the `pn-scattering` anchor travels to it).

**Gate**: `-E -W --keep-going` (the `-E` phantom-duplicate note above). L35 three-way grep +
L34 filename sweep per recipe; single-home every traveling label on its DEFINITION.

### ✅ slab_multigroup.rst — DONE @ `b7166ed6` (2026-07-20)

Authored 934 ln; index.rst 15,603 → 14,780 (−842/+19). 10 traveled labels single-homed (5
anchors `pn-scattering` / `n2n-reactions` / `sn-scattering-fission-operators` /
`choosing-inverse-realisation` / `cross-solver-eigenvalue-consumers` + 5 eq-labels
`multigroup` / `pn-scatter` / `flux-moments` / `real-spherical-harmonics` /
`addition-theorem`) + 2 new (`sn-slab-multigroup`, `sn-mg-eigenvalue-posing`); 5 labeled
eqs char-identical (enforcer-verified). Lands both slab_one_group forward promises
(in-scatter sum; the eigenvalue posing `(L+C−S−B)ψ=(1/k)Fψ`, `K=A⁻¹F`).

**User rulings (AskUserQuestion, all recommendations accepted):** (1) "Choosing the A⁻¹
realisation" COMPRESSED — foundations/operator_inverse_family owns the family catalog, the
chapter keeps SN consumption (MRO shadow → SweepOperator; ERR-058 reframe note;
rate-not-correctness); (2) cross-solver-consumers TRAVELED to the chapter (Phase-D
foundations promotion stays open); (3) "Backward references — Wave E status" traveled as
the operators' development-history note.

**Adjudications:** (n,2n) twin sections MERGED (physics/storage + algebra-slot rationale →
one section); k∞ block compressed to `mg-eigenvalue-problem` + `scattering-matrix-convention`
links (foundations/infinite_medium owns the derivation); informal iteration-refs list
RETIRED → formal `[TrefethenBau1997]` / `[Polizzi2009]` entries minted in the monolith
bibliography; label near-dup `cross-solver-eigenvalue-consumers` vs
`matrix-inverse-consumers` = NO-COLLISION (power_iteration's five families vs
MatrixInverseOperator's consumer ruling); BiCGSTAB sentence generalized (path retired Wave
E R2); KEigenvalue harmonized WITHOUT falsifying the class contract (constructor operand
``A`` = the invertible loss; B is method-layer — new explanatory paragraph; SN implements
the Protocol directly).

**Defects found + fixed:** pre-existing FALSE page-qualifier (`sn-angular-windowing`
claimed foundations/operator_algebra; actually defined in the monolith) — qualifier
dropped; 2 scattering.py docstring label-homes repointed (`theory-discrete-ordinates` →
`sn-slab-multigroup`); 2 windowing page-qualifiers repointed to the chapter; **enforcer
MUST-FIX** — a stale `M ≈ (A−S−F)⁻¹` display survived the harmonization in the
KrylovAcceleration bullet. Harmonized to the class's variadic algebra-of-record
`M ≈ (A−Σᵢgᵢ)⁻¹` (for SN, `(L+C−S−B)⁻¹`); the enforcer's proposed `(A−S)⁻¹` was
**overridden** — it drops B, contradicting iteration.py's own "the matvec IS the honest
(L+C−S−B)·ψ". The same stale pair fixed at iteration.py:747/758 (they contradicted the
adjacent algebra-of-record).

**Deferred (harmonize when their chapters come):** the iteration.py MODULE-docstring
posings (ln 1/7/32/861) + operator.py:9 keep the dated `(A−S−F)ψ=q_ext` local posing —
internally consistent, no adjacent contradiction; they fall to the algebra/solver-chapter
passes (+ #231 Phase 2 code-prose rebalancing).

## ⏭ ACTIVE = `cartesian_multid.rst` (space broadening) — disassembly RATIFIED (user, 2026-07-21)

**The gather windows** (line #s advisory @ 14,780-ln monolith; re-derive by title):
W1 "Cartesian 2D" transport H2 (276–293, eq `transport-cartesian-2d`) · W2 "Cartesian 2D
Balance Equation" H2 (344–398, `balance-cartesian-2d` + eq `dd-cartesian-2d`) · W3a the
UBLD core (2329–3812: Branch-1/2, S2 wiring, two-moment axes, Σ_s⊗I lift,
SpatialMomentSpace, moment matvec+scan — 11 anchors, 25 eq-labels) · W3b the 2-D LD
stress MMS D5b-S4 (3813–5168, 5 anchors + 6 eq-labels `ld-cartesian-2d*`) · W4
anti-diagonal wavefront (5169–5299, `sweep-wavefront`) · W5 octant DAG (5300–5917,
`sweep-dispatch-relayering` +2) · W7 Angular windowing H1 (12407–13543, 15
`sn-angular-windowing*` anchors + 5 eq-labels) · W8 boundary G-S (13546–13943,
`si-gauss-seidel-reification`; the `sn-iteration-primitives` H1 DISSOLVES — anchor has
ZERO inbound).

**User rulings (2026-07-21):**
1. **W3b → the SN section inside verification** (= `methods/sn/verification.rst`, the
   designated temp→V&V-part home) per the analytical-order≠MMS ruling. **NEW V&V-part
   architecture (user, verbatim intent):** the verification ENTRY page explains the
   correctness/science/logic of verification (what the `vv-principles` skill +
   `reference.md` encode); a SUMMARY page compiles all verification RESULTS across ALL
   transport methods (per-method tables, key metrics); per-method SECTIONS hold the case
   descriptions (SN's gets W3b). Full build-out = the V&V-consolidation task (#10 in the
   session tracker) — record there; this chapter only lands W3b in the SN section.
2. **W6 "Unified sweep dispatch" STAYS** (it carries a SUPERSEDED banner — Wave-D
   archaeology preserved-as-origin; live dispatch = `LossRepresentation` polymorphism in
   `loss_representation.rst`). The chapter authors its dispatch section FRESH against the
   live story, cross-linking loss_representation.
3. **W7 + W8 both travel** (2-D/multi-D iteration machinery; iteration-primitives H1
   dissolves cleanly).
4. **STAGED atomic commits** (deviation from one-commit-per-chapter, for cause at 3,900
   source ln): Stage 1 = posing+balance+sweep [W1/W2/W4/W5] + chapter skeleton; Stage 2 =
   UBLD [W3a → chapter] + W3b → sn/verification.rst; Stage 3 = windowing+G-S [W7/W8].
   Each stage: cut+land atomically (no twin homes at any point), `-E -W` gated,
   enforcer-reviewed (L-013), committed, one-line catalog status appended (compaction
   safety). Full checkpoint after Stage 3.

**Stays (adjudicated):** "Cell update strategies" H1 (cross-geometry strategy contract;
chapter cross-links CellVisit/DAG H2) · the τ-campaign + space⊗angle separability capstone
(684–1843, #236 characterization) · the Transport-Equation / Balance / Sweep-Algorithm H1
headers+intros (parent staying curvilinear H2s; pointers extend) · W6 (ruling 2).

**L35 blast radius (identified pre-cut):** `wavefront_cochain.rst:67,585` page-qualifiers
naming the old home of `sweep-dispatch-relayering` (Stage 1) ·
`tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py:47` docstring +
`orpheus/derivations/continuous/mms/sn.py:1402` comment naming `ld-cartesian-2d`'s home
(Stage 2) · rump pointers per stage. Gate `-E -W --keep-going` throughout.

**Stage 1 ✅ @ `e0a4e8b3` (2026-07-21)** — W1/W2/W4/W5 landed: chapter `cartesian_multid.rst`
962 ln (posing · 2-D balance · anti-diagonal wavefront · octant DAG · fresh thin
representation-layer section); index 14,780 → 13,964. Byte-exact SPLICE method (script
extracts spans; only transforms = underline promotion `-`→`=`/`~`→`-` + one in-span L35
"above" fix + the **Option-Y respelling** of the §15A.2 caption — sweep = `(L+C)⁻¹`, never
"A⁻¹ = sweep"; `A` bound locally for the verbatim eq). 11 labels single-homed. Enforcer:
COMMIT-READY zero MUST-FIX (its dual-`A` CONCERN = the Option-Y fix, applied; Fork-B2
provenance de-twinned to the capstone). L35 ×4 external/in-span; loss_representation intro's
pre-existing half-falsehood fixed. `-E -W` 0 warnings.
**⏭ NEXT = Stage 2:** W3a (UBLD core, grep "Multi-dimensional LD: the tensor-product
bilinear" — anchor `ld-ubld-multidim`) → chapter §LD (insert BEFORE the wavefront section?
NO — reading order ruled: balance → wavefront → LD-as-higher-order-closure AFTER the
dispatch section, before "What broadens next"); W3b (grep "The 2-D Cartesian LD stress MMS",
anchor `ld-cartesian-2d`) → `methods/sn/verification.rst` per ruling 1; Sweep-Algorithm H1
pointer revised ("multi-dimensional LD" clause drops); L35: the test docstring + mms/sn.py
comment sites; expect stale-vs-code adjudications in the S2/S3 "owed" notes — verify against
code at reading time.
