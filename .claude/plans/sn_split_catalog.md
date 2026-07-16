# SN monolith split — the Phase C execution catalog

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
   for the end. (`rst_skeleton.py <file>` in the job tmp re-dumps the tree.)
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
| — | `transport_equation.rst` | "The Transport Equation" (temp SN ch → foundations) | (none) | 1 | pending |
| — | `verification.rst` | "Verification" + "Numerical Sensitivities" (temp → V&V part) | `sn-keff-estimator` | 1 | pending |
| — | `history.rst` | "Development history" | `sn-development-history` | 1 | ✅ DONE |
| 5 | `discrete_balance.rst` | "The Discrete Balance Equation" (balance/α/geom-factor/streaming-eq/Morel-Montry dip — up to the WDD weights) | (none) | 2 | pending |
| 6 | `spatial_closures.rst` | rest of "The Discrete Balance Equation" (WDD & Morel-Montry weights, closure substitution) + "Cell update strategies" | `cell-update-strategies` | 2 | pending |
| 9 | `sweep_1d.rst` | "Sweep Algorithm" intro + Cartesian-1D cumprod + affine outflow | `sweep-algorithm` | 2 | pending |
| 7 | `linear_discontinuous.rst` | "Sweep Algorithm" → UBLD tensor-product bilinear cell system | — | 2 | pending |
| 10 | `sweep_multid.rst` | "Sweep Algorithm" → Cartesian-2D wavefront + octant DAG | — | 2 | pending |
| 11 | `curvilinear.rst` | "Sweep Algorithm" → curvilinear sequential sweep + unified dispatch + Carlson D/F + ERR-058 + route-(a); PLUS "The curvilinear anisotropic-MMS floor" | `sn-curvilinear-aniso-norm-reconciliation` | 2 | pending |
| 4 | `scattering_and_moments.rst` | "Scattering" + "Scattering and fission as LinearOperators" + "Angular windowing" | `sn-scattering-fission-operators`, `sn-angular-windowing` | 3 | pending |
| 8 | (append to `loss_representation.rst`) | "InvertibleOperator: the streaming-collision operator algebra" | `sn-streaming-operator` | 3 | pending |
| 13 | `iteration.rst` | "Krylov inner solver" + "Iteration primitives (operator algebra)" | `sn-iteration-primitives` | 3 | pending |
| 14 | `solver.rst` | "SNSolver as an operator-algebra coordinator" + "Consuming the frame in SN" | `sn-solver-operator-algebra-coordinator`, `sn-consuming-the-frame` | 3 | pending |
| 2 | `algebra.rst` ★★ | THE SPINE — synthesized: dedup the `(A−S−F)ψ=q` spine (restated 3× at ~13330/13524/13953), collapse `choosing-inverse-realisation` vs `inverse-application-driver`/`green-operator`/`matrix-inverse-operator`, `cross-solver-eigenvalue-consumers` vs `matrix-inverse-consumers` | (new) | 3 | pending |
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

**Progress: 3 chapters (135 + 408 + 1,019 = 1,562 ln), index.rst 20,571 → 19,012 (−7.6%).**
Delegated model validated across three scales, zero warnings throughout.
