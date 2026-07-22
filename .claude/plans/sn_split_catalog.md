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
| 5 | ~~`discrete_balance.rst`~~ | **SUPERSEDED by progression** — the Cartesian balances → `slab_one_group`/`cartesian_multid`; the curvilinear balance H2s gathered by `curvilinear_one_group` @ `6b8b0ff0` (the Discrete-Balance H1 stays as a 3-chapter router) | (none) | — | absorbed |
| 6 | ~~`spatial_closures.rst`~~ | **SUPERSEDED by progression** — generic closures live in `foundations/discretization`; WDD/Morel-Montry + closure substitution (+ the #236 τ/product capstones per steer Q3) gathered by `curvilinear_one_group` @ `6b8b0ff0`; "Cell update strategies" H1 STAYS (cross-geometry contract, Stage-1 adjudication) | `cell-update-strategies` | — | absorbed |
| 9 | ~~`sweep_1d.rst`~~ | **SUPERSEDED by progression** — Cartesian-1D cumprod + affine outflow gathered by `slab_one_group` @ `2c60d6a5`; the "Sweep Algorithm" H1 intro stays as the router pointer | `sweep-algorithm` (still in index) | — | absorbed |
| 7 | ~~`linear_discontinuous.rst`~~ | **SUPERSEDED by progression** — the UBLD core → `cartesian_multid` @ `eb96f424`; the 2-D LD stress MMS → `verification.rst` (same commit) | — | — | absorbed |
| 10 | ~~`sweep_multid.rst`~~ | **SUPERSEDED by progression** — wavefront + octant DAG → `cartesian_multid` @ `e0a4e8b3` | — | — | absorbed |
| 11 | ~~`curvilinear.rst`~~ | **SUPERSEDED by progression — DONE**: sequential sweep + pole closure + sweep-frame matvec + route-(a) → `curvilinear_one_group` @ `6b5848cd`; ERR-026-WaveE + Phase A + Carlson D/F + ERR-058 + #196 → `curvilinear_numerics` @ `550cee49`; the floor → `verification.rst` @ `50dcdfe1`; the W6 dispatch head STAYS (ruling 2) | `sn-curvilinear-aniso-norm-reconciliation` (→ verification) | — | absorbed |
| 4 | ~~`scattering_and_moments.rst`~~ | **SUPERSEDED by progression** — "Scattering" + "Scattering and fission as LinearOperators" gathered by `slab_multigroup` @ `b7166ed6`; "Angular windowing" → `cartesian_multid` @ `1722d438` | — | — | absorbed |
| 8 | ~~(append to `loss_representation.rst`)~~ | **DONE @ `0b96db7f` (2026-07-21, user-ratified SPLIT-BY-CONTENT — the whole-append superseded)**: the two-closure chronicle + Wave-E record → loss_representation "History and rationale" FIRST H2 (Wave-E its H3); the scattering adjoint → the new `adjoint.rst`; the H1 rump DISSOLVED (router intro's facts live in slab_one_group; refs split — L&M §10 → adjoint, A-L merged, T&B new row). The LAST dual-A local-convention zone harmonized at the move; the chronicle's falsified ERR-026 narrative tombstone-extended (ERR-058/#196 anchors) | — | ✅ DONE |
| 13 | ~~`iteration.rst`~~ | **SUPERSEDED by progression** — Krylov H1 → `slab_one_group` @ `2c60d6a5`; KEigenvalue/choosing/operand/FEAST/cross-solver → `slab_multigroup` @ `b7166ed6`; the boundary G-S schedule → `cartesian_multid` @ `1722d438` (the `sn-iteration-primitives` H1 DISSOLVED; its Wave-E-R1 fact → slab_multigroup) | — | — | absorbed |
| 14 | `solver.rst` | **DONE @ `48750bc4`** — both H1s traveled identity-promoted, harmonized-at-move (Option-Y + current-truth transforms); companion numerics respell @ `19e8b688` | `sn-solver-operator-algebra-coordinator`, `sn-consuming-the-frame` | 3 | ✅ DONE |
| 2 | ~~`algebra.rst`~~ | **CLOSED AS ABSORBED (user-ratified 2026-07-21).** The architecture (sn_book_architecture §4) supersedes the topic chapter; the content contract is realized by DISTRIBUTION — the posing spine single-homed at slab_multigroup `sn-mg-eigenvalue-posing` (`(L+C−S−B)ψ=(1/k)Fψ` · `A_loss=λM` → foundations `eigenvalue-posing` · `K=A⁻¹F` · the resolvent-reuse story); the §3.2 chain (invariant → why-invertible=lower-triangular → matrix=rawest → strategy-encoding) verbatim in slab_one_group "The within-group operator"; the strategy catalog = loss_representation (triangularity theorem, four representations, consumption modes); the shared §2.2 types = foundations/operator_algebra "The intrinsic operator types" (#226, the corrected S=kernel form). The residual restatements (solver framing H2, windowing Key Facts) are harmonized local-consumption one-liners, NOT re-derivations (the third site — iteration-primitives — dissolved at cartesian S3). ~~`choosing-inverse-realisation`~~ DONE @ `b7166ed6`; ~~label near-dup~~ RULED no-collision | (none) | — | ✅ absorbed |
| 15 | `adjoint.rst` | **CHAPTER MINTED @ `0b96db7f`** (extraction half pulled forward at ch8): anchor `sn-adjoint`, thin frame (3-layer scope + 4 Key Facts) + the traveled `sn-scattering-adjoint` section promoted to its spine. **Authoring half PENDING** — daggered posing / ψ* / φ*-consumers land with #276 A4/A5 (+ #51/#281); the traveled body's `(L+C−S)ᵀ` shorthand gets the B-note treatment then (the fresh frame already carries it) | `sn-adjoint`, `sn-scattering-adjoint` | 3 | extraction ✅ / authoring pending |
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

## ✅ `cartesian_multid.rst` (space broadening) — COMPLETE (3 staged commits, 2026-07-21)

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
**Stage 2 ✅ @ `eb96f424` (2026-07-21)** — W3a (UBLD core, 2261–3744 @ the 13,964-ln
monolith) → the chapter between "Choosing the schedule" and "What broadens next" (964 →
2,469 ln); W3b (the stress MMS + Leg A #247 + Leg B #251 + S9 #257, 3745–5100) →
`verification.rst` at the END of the MMS case run, before `sn-composite-fixed-source`
(2,625 → 3,982 ln); index → 11,125. Byte-exact splice (promotion count-asserted 30+24
headers incl. W3b's 4th level `'`→`^`; enforcer certified 31/31 labeled eqs
char-identical, W3b zero drift, hierarchy zero level-jumps). The S3-A "PARTIAL"/"What is
still owed" chronicle got two present-state annotations (every owed item verified LANDED
in code: `spatial_moments=` threaded from field properties; `_slope_fold`/`schur_xV`/
`scan_slope_face_source`/`scan_reconstruct`/`octant_moment_frame_signs`/
`_moment_broadcast_sigma`/`moment_scan_closure` all live) — the owed-lists preserved as
campaign-time boundary records. L35: Sweep-Algorithm H1 pointer extended (chapter +
verification); 11 stale qualifiers repointed (sn.py ×5 + test docstrings ×5, mostly
pre-existing Tier-1 debt caught by the enforcer's blast-radius sweep, + the Stage-1-missed
`loss_representation.rst` `dd-cartesian-2d` page qualifier). `-E -W` 0 warnings ×2; zero
staying-monolith back-refs; all cross-span refs path-immune.

**Deferred → the algebra/harmonization pass:** the dual-`A` overload in the traveled UBLD
text — `ld-ubld-cell-system` (per-cell dense `A = G + F_out + Σ_t M`) vs
`ld-ubld-unified-moment-residual` (operator `A = L+C−S`) — both byte-exact-traveled, both
inline-defined at point of use (enforcer CONCERN, recurring from Stage-1 §15A.2). Remedy
in ONE pass across the UBLD chapter: rename the per-cell matrix `A_cell` (or add a
reconciling note), reserving `A` for the honest operator.

*Method (proven Stages 1–3):* author the chapter/verification insertions as a skeleton with
`<<<SPLICE:Wx>>>` markers, then a python script: extract spans byte-exact from index.rst by
line range (assert content anchors), PROMOTE underlines one level (`-`→`=`, `~`→`-`; only
pure runs below nonblank titles — never transitions after blanks), apply declared in-span
fixes by asserted string-replace, splice markers, cut monolith descending, do pointer
replacements with unique-landmark assertions. Zero retyping ⟹ char-identical math by
construction (the enforcer certifies trivially). Stage-1 reference script:
`stage1_cartesian_multid.py` pattern (job tmp, ephemeral — the pattern above is the record).

**Stage 3 ✅ @ `1722d438` (2026-07-21)** — W7 (Angular windowing H1, 8751–9887 @ the
11,125-ln monolith; **IDENTITY promotion** — the chapter's flat multi-`=` top level absorbs
H1/H2/H3 verbatim, census 1/8/8) + W8 (boundary G-S H2, 9906–10287; promoted `-`→`=`/`~`→`-`/
`^`→`~`, census 1/4/4) → the chapter after UBLD, before "What broadens next". The
`sn-iteration-primitives` H1 wrapper (9888–9905) DISSOLVED — zero inbound; its one
info-bearing fact (Wave E Round 1 #163, the iteration-primitives lift into
`orpheus.numerics.iteration`) preserved as the leading paragraph of slab_multigroup's
"Development history" (completing R1→R2→R3). 7 declared in-span fixes, every one a
PRE-EXISTING cross-file falsehood ("above/below" fossils for content that moved to
`foundations/wavefront_cochain` and `slab_one_group` in earlier extractions; two same-file
qualifiers on `si-gauss-seidel-reification`). One Key Facts bullet (windowing + G-S);
external repoints: `slab_one_group` ×1 + `operator_inverse_family` ×2. 21 labels
single-homed (the 4 `.. vv-status:` markers traveled with their eq-labels). Enforcer:
COMMIT-READY zero MUST-FIX, zero undeclared drift — 5/5 labeled eqs + all 10 math blocks
char-identical; the windowing-retyped → boundary-G-S seam judged BETTER than the wrapper
provided. `-E -W` 0 warnings first pass. index 11,125 → 9,588.

### ✅ CHAPTER COMPLETE — `cartesian_multid.rst`, 3,998 ln

Reading order: posing → 2-D balance → anti-diagonal wavefront → octant DAG → choosing the
schedule (representation layer) → UBLD tensor-product closure → angular windowing → boundary
Gauss-Seidel → what broadens next (curvature = Part B). Across the 3 stages: index.rst
14,780 → 9,588 (−5,192); 78 labels single-homed; every stage enforcer-certified byte-exact
and `-E -W`-gated. Carried deferral: the **dual-`A` one-pass harmonization** (above).

**⏭ IN FLIGHT = `curvilinear_one_group.rst` (Part B) — steered + Stage 1 landed (2026-07-21).**

**The four steer rulings (user, 2026-07-21 — all recommendations ratified):**

1. **Topology: TWO chapters, route-(a) in base.** `curvilinear_one_group` = posing →
   balance → WDD/M-M closure → sequential sweep → pole closure (Phase B) → sweep-frame
   matvec (Phase C) → route-(a) production ψ½. `curvilinear_numerics` (name provisional) =
   the #168 seed-strategy campaign record: ERR-026 Wave-E status, retired Phase A, Carlson
   Phase D, Phase F, ERR-058. Arch §4 open question 2 RESOLVED. (Finding recorded: route-(a)
   is present-state production math — the arch gathers-list "pole/ψ½ start · radial
   characteristics" content — despite the provisional arch row parking it in numerics; its
   banner's "saga documented above" pointers become cross-chapter :ref:s.)
2. **The aniso-MMS floor H1 → `verification.rst`** (three-layer V&V ruling, W3b precedent);
   the same ruling governs Phase C's four nested MMS-gate subsections at authoring time.
3. **The #236 capstones travel with WDD** into `curvilinear_one_group` (S3-A precedent). ✅
   done at Stage 1.
4. **Staging: 4 staged commits** — S1 balance+closure core; S2 sweep+pole+matvec+route-(a);
   S3 the campaign chapter; S4 floor→verification + closeout.
   (Also found at disassembly: ZERO genuine multigroup content in the curvilinear zones —
   `curvilinear_multigroup` receives nothing from this split; authored fresh later.)

**✅ Stage 1 @ `6b8b0ff0` — balance+closure core (monolith 329–1854 → chapter, 9,588 →
8,065).** Skeleton fresh-authored (anchor `sn-curvilinear-one-group`, "angle enters the
walk", §3.2 chain, 6 Key Facts). Promotion {`-`→`=`, `~`→`-`, `^`→`~`} census {8,6,12}; 29
math blocks / 13 eq labels / 7 anchors char-identical + single-homed (enforcer-verified
mechanically: 31 hunks = 26 promotions + 5 declared fixes, zero undeclared drift; zero
dual-`A` concern). Declared fixes F1–F5 all pre-existing falsehoods (retired
`_sweep_1d_*` symbols; stale `derivations.sn_contamination` path; retired
`SNMesh.alpha_half`/`redist_dAw` accessor family → `mesh.reduced.*`). Monolith: Discrete-
Balance H1 → 3-chapter router; M2 cell-update fossil; M3 toctree. External: E1
operator_tensor_network Step-2–3 qualifier; E2 balance.py source-of-truth; E3
test_ld_slope_frame (cartesian-stage miss); + 2 enforcer-caught test stragglers
(test_streaming_equilibrium_curvilinear L0-gate pointer MUST-FIX;
test_space_angle_separability mm-weights half). Enforcer lesson for S2–S4: partition the
tree-wide `methods/sn/index` grep hits by the per-label HEAD-vs-worktree oracle BEFORE the
commit — the test-tree stragglers of THIS stage's move are the recurring miss class.

**✅ Stage 2 @ `6b5848cd` — the sweep core (monolith 8,065 → 5,560; chapter → 4,160).**
Three spans: W2a sequential sweep (738–773, {`-`→`=`}), W2bc Phase B pole closure + Phase C
sweep-frame matvec contiguous (1060–2479, {`~`→`=`, `^`→`-`, `'`→`~`} {2,20,2}), W2d
route-(a) (5274–6328, {`-`→`=`, `~`→`-`}). 19 anchors + 11 eqs single-homed; parity 48=29+19.
**Adjudication (enforcer concurred): Phase C's MMS-gate subsections = tombstoned chronicle —
travel WITH Phase C, NOT verification.rst** (the W3b precedent is for LIVING cases; the floor
goes to verification.rst at S4). F1 retired `_sweep_1d_*` → the live walk; F2/F3 saga
"above/below" cross-file; F4/F5 bc-extraction home = boundary_conditions.rst:352 (NOT
operator_algebra — moved by the corpus campaign; pre-existing). 4 Key Facts bullets added.
Enforcer MUST-FIX ×5: stale `:doc:index` page-qualifiers on moved anchors in
boundary_conditions.rst (655/738/743/3804) + coupled_block_operator.rst (65) — the
foundation-docs inbound-ref class the code/test sweep misses; ADD `docs/theory/foundations/`
to the per-label oracle sweep at S3/S4. Two staying ERR-058 "Phase A–F above" claims made
honest (B/C now in the chapter).

**✅ Stage 3 @ `550cee49` — `curvilinear_numerics.rst` MINTED (2,981 ln; monolith 5,560 →
2,627).** Two spans: the W6-nested prologue (888–1028, {`~`→`=`} {2}) + the contiguous
D/F/ERR-058/#196 block (1029–3822, {`-`→`=`, `~`→`-`, `^`→`~`} {3,40,9}). 18 anchors + 12
eqs single-homed; parity 18/18. FB1 dead `:func:` on retired `_sweep_1d_spherical` →
literal. Router #2; toctree. Externals: bc.rst ×2 Phase-D/F qualifiers → numerics; the
kinf-recovery compound split applied to BOTH test twins (krylov-precond + invertible-operator
— the second enforcer-caught, SHOULD-FIX). Enforcer: COMMIT-READY, ZERO MUST-FIX; W6 rump =
head + 2 H3s, no dangling refs; all Stage-2 findings verified resolved.

**✅ Stage 4 @ `50dcdfe1` — the floor → `verification.rst` (monolith 2,627 → 1,949;
verification 3,982 → 4,662).** First DEMOTION splice: {`=`→`-`, `-`→`~`, `~`→`^`} {1,8,6};
spliced before `ld-cartesian-2d`; 5 anchors + 4 eqs single-homed; F1/F2 "above" fossils →
`:doc:`curvilinear_numerics``. Blurb rewritten honest; 7 external repoints (+ the enforcer
MUST-FIX `structured_geometry.rst:288` — the foundation-doc twin class again). ⚠ the S4
script aborted on the row-8 stale anchor claim (now corrected in the table) — the doc writes
had landed; externals applied by follow-up script, no double-application
(enforcer-triple-checked).

### ✅ PART B COMPLETE — the curvilinear campaign (4 staged commits, 2026-07-21)

`curvilinear_one_group.rst` (4,160 ln: posing → balance → WDD/M-M closure + #236 capstones →
sequential sweep → pole closure B → sweep-frame matvec C → route-(a) ψ½) + NEW
`curvilinear_numerics.rst` (2,981 ln: the #168→#282 seed-strategy campaign record) + the
floor in `verification.rst`. Monolith 9,588 → 1,949 (Synopsis · Architecture · Transport Eq ·
Balance router · Cell-update H1 · Sweep/dispatch rump · InvertibleOperator · SNSolver ·
frame · Gotchas · References · toctree). Commits: S1 `6b8b0ff0` · S2 `6b5848cd` · S3
`550cee49` · S4 `50dcdfe1` (+ checkpoints/lessons between). Enforcer: char-identity verified
mechanically at every stage; every MUST-FIX closed within its stage.

**✅ The dual-`A` one-pass harmonization @ `3872cf9d` (user-ordered, 2026-07-21) — the
standing CONCERN DISCHARGED.** Option-Y realized across the campaign pages: the per-cell
UBLD matrix renamed `A_cell` (eqs + the record docstring in `ld_ubld.py`; the ``A == A``
pin literals KEPT matching the code dict key, bridged); the dated `A = (L+C−S)` binding
de-bound to the explicit equivalence (a type-improvement, enforcer-certified); the W7 +
`_maybe_window` bindings sharpened to "the swept composite — the inner kernel of the honest
(L+C−S−B)⁻¹, never the full solve" (local `A` kept: the swept object ranges over {L+C, M});
Phase-C apply↔sweep respelled `L+C` end-to-end (sentence, gates, intro, Key Facts); the
route-(a) banner carries the section binding (augmented loss `L+C`, System A ⊕ B).
Neighbor sweep clean (slab/loss_representation/foundations all honest-`A`). REMAINING
local-convention zone: the monolith InvertibleOperator H1 (`A = Ω·∇+Σ_t`, "A⁻¹ is the WDD
asymmetric sweep") — harmonize AT its ch8 move to loss_representation. Arch §4 open question 2 RESOLVED
(numerics minted); question 3 (multi-D curvilinear) remains out of scope — nothing in the
monolith. `curvilinear_multigroup` receives NOTHING from the split (zero multigroup content
found) — authored fresh when Part B's energy chapter comes.

**✅ `curvilinear_multigroup` — DONE @ `dc9d9a56` (fresh-authored, 2026-07-21).** Steer
ratified **(a) thin standalone chapter** ("Curvilinear multigroup: the group axis rides the
walk", anchor `sn-curvilinear-multigroup`, ~340 ln); toctree slot between
`curvilinear_one_group` and `curvilinear_numerics`. **The queued brief's parenthetical was
WRONG and the explorer grounding corrected it before authoring**: Morel–Montry τ/c are
GROUP-BLIND angular weights from (μ,w) alone (the #236 Step C ruling) — NOT per-group; the
per-group data is the WDD denominator/attenuation (Σ_t,g·V, the sole group dependence of
the cell update — `cell_balance_for_streaming`) and ψ½. The chapter states the honest
THREE-TIER claim (structure: no ng axis / data: group-diagonal lanes / coupling: S+F only)
plus the monolithic-iteration fact (NO group loop; full multigroup S inside the
within-group system, fission external; NOT group-by-group Gauss-Seidel — the two spellings
of "per group" kept apart). ZERO new eq labels (posing = :eq:`multigroup` ∘
:eq:`balance-general` by reference; 2 unlabeled displays). Enforcer (FRESH instance):
SHIP-WITH-FIXES → all 4 landed in-commit (byte-for-byte overclaim → "same operators,
consumed by composition"; per-geometry resolution-ledger precision; assembly-loop count;
source-wrap) → COMMIT-READY; "perfect Option-Y honest algebra"; every load-bearing claim
verified EXACT vs the live tree. Sibling edits: index router clause + toctree;
curvilinear_one_group:44 live :doc: link; slab_multigroup closing bullets linked. Gate
`-E -W` GREEN. **#302 filed**: the family-wide silent py-domain xref condition (65+
unresolved refs — sweep/spatial/coupled-system modules have no autodoc target; measured
counts on the issue). Agent-memory round @ `1f2f3c7d`. **Part B is now COMPLETE including
its energy rung; the broadening progression (arch §3.4) is fully realized.**

**✅ ch8 — the InvertibleOperator H1 dissolution @ `0b96db7f` (2026-07-21).** Steer:
SPLIT-BY-CONTENT (the row's whole-append superseded by disassembly) + MINT thin
`adjoint.rst` now (ch15's extraction half — the label moves ONCE). The chronicle + Wave-E
→ loss_representation's "History and rationale" (first H2 + its H3), shelved beside the
one-walk theorems they pre-date; Option-Y harmonization of the LAST local-convention zone
(`A = Ω·∇+Σ_t` → `L+C`); the tombstone note extended with the falsification record
(ERR-058/#195 seed-family re-adjudication + #196 SI≡Krylov, both anchors cited).
Current-truth transforms: Wave O/O.2b LANDED the analytic reverse-direction adjoint matvec
— the dense-transpose-fallback claims and the "Forward" bullet updated; 2 fossil pointers
to the retired FD-family docs dropped. **The l1_analytical fossil family** (6 stale path
strings across 4 pages, found by the L34 sweep) fixed in-pass — aniso tripwires →
`verification/mms/test_curvilinear_aniso_convergence.py`, kinf → `verification/analytical/`
(2 foundations pages), the flat-flux-identity mention → retirement note (file retired with
Legacy/BFF, `2f24ae4a`); exactly ONE historical residual by design. index 1,952 → 1,736;
toctree + router updated. Enforcer (resumed instance): **COMMIT-READY, zero MUST-FIX** —
"every moved span char-diffs to exactly the declared transforms"; the one deferrable
CONCERN (B elided in the daggered-posing shorthand) closed in-commit in the fresh frame.
Gate `-E -W` GREEN. Agent-memory @ `b55257a5`.

**✅ ch14 `solver.rst` — DONE @ `48750bc4` + `19e8b688` (2026-07-21).** Steer: all three
recommendations ratified, Q1 UPGRADED by the user to "also respell the code posings" (its
own commit). The two H1s traveled IDENTITY-promoted (flat multi-`=` top level; fresh
opener + Key Facts frame); index 1,736 → 1,149 — **the predicted rump exactly** (Synopsis ·
Architecture · Transport Eq · Cell-update H1 · Sweep rump · Gotchas index · References ·
toctree, all ruled STAYS). Declared transforms: **Option-Y** (the `(A,S,F)` triple →
`(L+C, S, F)`; `(A−S)ψ=q_ext` → `(L+C−S−B)ψ=q_ext`; "A⁻¹ is the WDD asymmetric sweep" →
the resolvent `(L+C)⁻¹`; the H2 title → "The (L + C − S − B)·ψ = (1/k)·F·ψ framing at the
solver level"); **current-truth** (stale ERR-026 + "Reflective-BC equation map only
(Round 3…)" dropped — live solver.py confirms retired; "Two Inner Solvers" honest: T⁻¹ →
`(L+C)⁻¹`, both paths iterate the typed composite, the Krylov system = the honest
`(L+C−S−B)` with gains in the operator, the unverifiable "~0.97 @ 421 groups" figure
dropped for ρ=c → slab_one_group, ERR-053 restart×n_dof memory note); **move-created
fossils** ("the streaming-collision history section above" → `loss-rep-history`;
"elsewhere on this page" → verification); T18 the Gotchas "scale bridge above" →
`sn-keff-estimator`; E1–E3 the MoC/CP/MC `sn-keff-estimator` page-qualifiers → sn/solver;
E4/E5 the pre-existing transport-cartesian test-docstring fossils → slab_one_group. 5
labels single-homed (vv-status markers traveled). **The numerics respell** (`19e8b688`):
the iteration.py/operator.py module-head posings → the variadic `(A − Σᵢgᵢ)`
algebra-of-record (A = the resolvent OPERAND = the literal constructor parameter,
unrenamed; SN binding `A = L+C` ⊕ gains `(S, B)`; fission never a gain in the eigenvalue
posing); the KEigenvalue-layer `(A,S,F)`/`A_loss = A−S` sites deliberately KEPT (B is
method-layer — that class's own harmonized contract) — **the slab_multigroup dated-posing
deferral is CLOSED**. Enforcer (FRESH instance, persists this session): **COMMIT-READY,
zero MUST-FIX** — "exactly 11 non-equal hunks, each a declared transform or fresh-frame
insertion"; the ~380-ln keff record + the whole frame H1 byte-identical; the cross-layer
dual-`A` ruled a NON-finding (the mapping is spelled in-code both directions —
iteration.py's `A = L+C` head binding + its `A_loss = L+C−S−B` comment). Gate `-E -W`
GREEN. Agent-memory @ `8e41e0a1`.

**Phase C's structural-move program is COMPLETE.** Every movable chapter is out of the
monolith; the index is the ruled rump. What remains in the campaign is SYNTHESIS and
polish, not moves.

**✅ ch2 DISPOSITION (user-ratified 2026-07-21): CLOSED AS ABSORBED — no `algebra.rst`.**
The disassembly showed the confirmed architecture's §4 supersession holds AND the four
content directives are realized by distribution (the realization map is in the table row
above). With this, **every row of the chapter catalog is terminal** (✅/absorbed/N/A)
except ch15's #276-blocked authoring half and ch1 `placement.rst` (Phase-D scoped) —
**Phase C is CLOSED end-to-end.**

**✅ Phase D DONE + Phase I LANDED (both 2026-07-21, one session).** Phase I: both surveys
delivered (`.claude/plans/phase_i_survey_{larsen_morel_2010,adams_larsen_2002}.md`, §3.2
rows filled @ `be765fce` — the niche thesis SURVIVED sharpened: verification ZERO in both
reviews; L&M spatial-half-only (BMC absent from its bibliography); A&L confirms no
discrete curvilinear Fourier exists; the A&L p. 140 GMRES-outside-DSA quote licenses
#200/#2; the #2 DSA reading map's four primaries all local). Phase D (ratified order):
(1) `methods/index` §2.3 rebuild @ `8c4fa608` (thin child of the root; the diffusion
moment-reduction gap → **#303**); (2) the SN router @ `366e0bec` — the four P5 tracks +
the 11-row `sn-symptom-table` (every row enforcer-verified to ANSWER its symptom);
(3) foundations/theory-root/verification **NO CHANGE — already at the P2 bar** (enforcer
adversarial CONCUR: candidate additions would TWIN existing routing; verification's
three-layer future = #10); (4) the corpus root's orient-and-delegate paragraph @
`60790f4b`. `placement.rst` deferred to Phase-H adjacency (user-ruled).

**✅ Phase E DONE (2026-07-21, one session) — the conventions part is WHOLE:**
`conventions/{index, notation, normalization, indexing_and_layout,
cross_section_conventions}`. Steer rulings: crosswalk source set = the SURVEYED canon
(Hébert · B&G · Stacey + the two Phase-I reviews), with the user's correction —
**`scratch/literature/Nuclear Computational Science - A Century in Review.pdf` IS
local** (the 2010 Azmy/Sartori volume: the surveyed L&M-2010 chapter + a Lewis
chapter); **INTERNAL CONSISTENCY = the prime directive** ("one concept, one spelling"
across code+docs; a deliberate two-binding case is documented at both ends — the
doctrine now opens `notation.rst`); the α-crosswalk homed at `normalization.rst`
(curvilinear cross-links it, machinery stays); 3 ruled commits landed; main agent
authored the two NEW pages, archivist executed the split (L-013; git rename-detect
97%), the fresh enforcer reviewed everything (READY×4; ONE MUST-FIX total: the ERR-039
catcher citation copied the catalog's own stale Round-1 test-ref — the catalog was the
armed trap, defused @ `d8e0f6bf`).

Commits: `f297d584` fix (infinite_medium stored-SigS triangularity lower→UPPER — its
own equation + the gendf canonical-order test both say upper; machine-header
`HarmonicMomentField`→`Flux`) · `8ef8bbd3` **notation** (symbol table with the
enforcing-test column · the 8-row literature crosswalk incl. the operator-letters row
(A&L A≡I−L⁻¹S ≠ ours), the same-symbol collisions row (L&M-2010 β/α, A&L angle-μ,
mfp-λ), and the ORPHEUS-internal dual-A row as a declared bridge · the 3-question
import boundary · defines `[LarsenMorel2010]`) · `608e9c9b` **normalization** (the W
law + the ERR-025 anatomy · the (2ℓ+1)/W prefactor LEDGER — each constant ONE home;
Hébert's 4π-vs-2 = ONE object (2ℓ+1)/W at two weight sums, code-true at
scattering.py:1695 · the 8-row α-recursion crosswalk = the M5 deliverable, merged from
the catalog four-ways + the researcher's six-row table · enforcer N1 taken: the
part-index M1 bullets thinned to a router pitch, equation numbers single-homed) ·
`d8e0f6bf` chore(vv) catalog test-ref sync · `359d4b4d` the split (all 10 labels stay
on `indexing_and_layout` — 9 external referencers resolve untouched; the one new label
`theory-conventions-cross-sections`) · `93138008` XS-page truth sync (NuSigF→`SigP`;
the retired per-attribute SNSolver surface → the `mat_xs` `MaterialXSField` single
source, PR-TYPED-1/2) · `75f68397` solver.py stale-comment retire · `2abb0b7f`
archivist lesson. New labels: `theory-conventions-notation`,
`notation-{symbol-table,crosswalk,import-boundary}`,
`theory-conventions-normalization`,
`normalization-{weight-sum,prefactor,alpha-crosswalk}`, `sn-alpha-normalization`,
`theory-conventions-cross-sections`.

**✅ Phase F DONE (2026-07-21, one session) — the archaeology sweep, P7.**
~90 retitles across 4 per-part batches, each enforcer-certified READY:
**batch 1 methods/sn @ `41a93d37`** (34 retitles; `curvilinear_numerics` ruled
KEEP-all — a campaign-record chapter BY CHARTER; + 2 stale-deferral content
fixes: T3's removed-xfail terminal note, loss_representation's "waits on #282"
retired; + the test banner @ `b0f7b25a`) · **batch 2 foundations @ `2220fe05`**
(25; the six BC "Step N" pipeline stages ruled false positives; bc:1446's
forward-looking "What remains for Wave O" flipped to the verified discharge
record) · **batch 3 references @ `3fd96012`** (28 + the continuous-µ-retreat
name harmonization; peierls/trajectory/fn ruled DESIGN pages, their
history/close-out records KEPT; the trajectory phase↔geometry rosetta kept as
the information-loss guard; the three "Numerical verification status"
headings de-collided by variant) · **batch 4 tail @ `689b9377`** (3).
Kickoff rulings (user): per-part 4 commits · draft-and-execute (map in each
commit message; enforcer verifies; veto-after). The ratified discriminator
held across ~90 retitles + ~30 KEEPs with zero reversals (enforcer closeout).
Issue-state ground truth pulled via gh (⚠ #9 and #229 were CLOSED though
memory said open — trust git/gh). **The residue a heading-scoped sweep cannot
reach = #304** (P10 label re-namespacing · code-side campaign comments ·
2 prose two-spellings, low/optional).

**⏳ IN FLIGHT = Phase G (backfill) — steer RULED 2026-07-22.** Order =
**G1 equation labels → G2 `:cite:` → G3 `:term:`**; **V&V slices FOLD INTO
task #10** (the three-layer part design dictates the slice format — no
provisional per-page slices). `:cite:` confirmed ruled (Option 1 per corpus
plan §9; bibtex WIRED but INERT — `sphinxcontrib.bibtex` + `refs.bib` with 12
orphan entries, zero `:cite:` uses, NO `.. bibliography::` directive; theory
carries 59 bracket-style defs, ALL UNIQUE keys — #231's "duplicated across
pages" claim is stale, the campaign consolidation dedup'd them). G1 scope:
~450 unlabeled displays (1127 math / 677 labeled), MINUS
thermal_hydraulics(25)/fuel_behaviour(16)/reactor_kinetics(4) — EXCLUDED
under the no-investment rule (subsystem fate undecided) → ~405 to adjudicate
in 3 parallel archivist batches: A = methods/sn (~92) · B =
foundations+methods-non-sn+conventions+root-verification (~106) · C =
references (~194, derivation-heavy — label the skeleton, not every vertebra).

**✅ G1 DONE 2026-07-22** — 135 labels / 405 adjudicated (A 38/96 @
`91477611` · B 66/111 @ `85f2222b` · C 31/198 @ `72146f0a`), ~270
deliberately bare (restatements / intermediate algebra / worked numbers).
Enforcer certify: all three batches CERTIFIED, ZERO must-fix (matrix
orphan-count 108→243 = independent collision-free proof). SC dispositions:
SC-1 multiplication-operator term collision → **ledgered as notation
crosswalk row 9** (domain-inherent: reactor-physics K=A⁻¹F vs
functional-analysis M[f]; label namespaces already disjoint — NO rename);
SC-3 bte/transport-equation + SC-4 pincell/pin-cell → **filed on #304**;
SC-2 balance-prefix (anchor-mirroring, deliberate) + SC-5 accepted as-is.
**vv-status tagging of the new labels DEFERRED to task #10** (per-equation
V&V classification = the slice format; ~43 of B's 66 sit on vv-status
pages, all documentation/definitional class; A flagged
sn-mms-nonvacuum-qext-mg → inherit sibling's `verified` if SymPy-covered;
`resolvent-object-gate` = natural verifies() target for the Mode-12 test).
Cross-doc prose-ref polish candidates (NOT done, out of labels scope):
cartesian_multid ~L3513 → :eq:`pi-r-equals-4pi-i`; sn/verification ~L4559
→ kinf-1g.

**✅ G2 DONE 2026-07-22** — the corpus is on sphinxcontrib-bibtex
end-to-end. Infra @ `30a073fe`: refs.bib 12→59 entries (keys verbatim =
drop-in swap; 33 CrossRef DOIs; **11 field-level errors in the old
definition texts corrected INTO the bib** — the bib outranks any old page
text; ledger = `.claude/agent-memory/literature-researcher/
refs_bib_g2_corrections.md`); conf.py `keylabel` pybtex style →
bibliography labels ARE the BibTeX keys, inline rendering
character-identical to the old brackets. Matrix regen @ `56f0705b` (G1
fallout: orphans 108→243; classification deferred to task #10). The
atomic swap @ `263a7a4b` (49 files, +277/−736): 233 uses → `:cite:` (213
rst + 20 docstring), 78 defs deleted (59 rst + 19 docstring — none of the
docstring modules is automodule'd, so that front is source-consistency),
`references/bibliography.rst` = the single citation home (toctree-wired);
13 emptied References headings removed grep-verified, 1 mixed kept
(singular_eigenfunction); RST footnotes + the two `[A]_` math
pseudo-sites preserved. Consolidations (one concept, one spelling):
PS1982→PomraningSiewert1982 · Sood1999→SoodLA13511_1999 ·
+MetcalfZweifel1968 (NSE 33(3):318–326 = Part II of the 1968 pair,
CrossRef-resolved; companion noted). Citation-hygiene residue for a
future editorial pass (flagged in the ledger, non-blocking):
McCormickKuscer1965 (the 1966 companion is cross-ref'd in the note),
ENDF102 (split keys if decay-data provenance wanted), MATPRO2003 (thin),
Bondarenko1964 (full author list available).

**⏳ G3 IN FLIGHT — `:term:` wiring.** 15 glossary terms, ZERO inbound
today. Policy: first significant mention per page, wrap-don't-reword,
inflections via ``:term:`x <term>``, no headings/math/literals,
sense-checked per site (transport sweep ≠ parameter sweep);
TH/fuel/kinetics excluded (no-investment rule); missing-term candidates
REPORTED, not minted. Then: gate → commit → Phase-G record to GitHub
#231 + tracker close. The CLI task tracker mirrors the queue (task #7 =
Phase G in_progress; #8 = Phase H unblocked; #10 = V&V part + slices +
vv-status tags; #15 = #304 residue; #16 = ch15).

**Then:** H (the root page — Phase-I input IN HAND; **UNBLOCKED**: the
#298/#299 fixes are in-branch @ `639cec9e`/`4425516c`, marked DONE @
`68c39c28`, the issues auto-close at push — the earlier "still gated" note
was the frozen claim; §3.6's 12 MUST-NOT claims still apply; placement.rst
rides its adjacency) · #304 (the Phase-F residue: P10 labels · code comments
· prose two-spellings) · ch15 authoring half (#276 A4/A5-blocked) · #231
Phase 2 (code-prose rebalancing) · task #10 (three-layer V&V part).

*Standing facts:* branch `docs/sn-doc-architecture` UNPUSHED; citations + eq-labels are
project-global; recipe = "The per-chapter recipe" section above (L35 three-way grep + L34
filename sweep per chapter — Stages 2/3 and ch8 found most breaks are PRE-existing fossils
from earlier extractions, so grep every window's directional prose and path strings against
the CURRENT homes of what they cite); the elegance-enforcer review instance persists
within a session — resume it; re-brief fresh only after compaction.
