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
7. **Gate:** `-E -W` clean (above). Confirm every moved label single-homed (`grep -c` in src = 0).
8. **Archivist reports** the diff + grep results + build; **main agent reviews diff-vs-catalog,
   re-runs `-E -W`, commits** with the correct trailer (archivist does NOT commit — Phase B model).

**Delegation:** ch3 done by main agent (pattern-prover). Tier-1/2 clean extractions →
**archivist**, one chapter per dispatch, brief carries the exact span + this recipe. Tier-3
(dedup-heavy) stays main-agent + elegance review (it collapses duplicated content, not a move).

---

## The chapter catalog (keyed on stable section title + anchor; line #s advisory only)

`§5#` = chapter number in the corpus plan §5 table. **Tier** = execution order (1 = clean
leaf first, 3 = dedup-heavy last).

| §5# | Chapter file | Source `=`-section(s) — grep this title | Anchor | Tier | Status |
|---|---|---|---|---|---|
| 3 | `angular_quadrature.rst` | "Angular Quadratures" | `quadrature-types` | 1 | ✅ DONE |
| 12 | `boundary_conditions.rst` | "Boundary Conditions" | `boundary-conditions` | 1 | pending |
| — | `transport_equation.rst` | "The Transport Equation" (temp SN ch → foundations) | (none) | 1 | pending |
| — | `verification.rst` | "Verification" + "Numerical Sensitivities" (temp → V&V part) | `sn-keff-estimator` | 1 | pending |
| — | `history.rst` | "Development history" | `sn-development-history` | 1 | pending |
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

## Status log

- **ch3 `angular_quadrature.rst`** — ✅ **MERGED-to-branch @ `5a61fdb5`.** Cut 366–500
  (135 ln), 1 label `quadrature-types`, 1 inbound ref (MoC, path-immune ✓), 0 directional
  prose, 0 substitutions. L35+L34 clean; `-E -W` build succeeded (0 warnings). index.rst
  20,571 → 20,436. Diff was a clean −135/+1 (section out, toctree line in). Pattern proven.
