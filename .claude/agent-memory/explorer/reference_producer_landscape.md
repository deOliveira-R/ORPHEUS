---
name: reference-producer-landscape
description: Shape of reference-solution producers under orpheus/derivations (#405 step 2, 2026-09-24): registries, types, duplicated problem spellings, census probes
metadata:
  type: project
---

Tree at 901f64ca (2026-09-24). Re-verify before repeating any count.

- The registry (`reference_values`: `_CASES` 47, `_CONTINUOUS` 16 eager + 13 lazy Peierls) mediates ~1 in 9 reference acquisitions from tests; the heavy generators (Peierls, trajectory resolvent, F_N, Case, Galerkin) are called inline (~759 AST call sites). A cache keyed on the registry alone misses most of the cost.
- Problem description has no single type: geometry kind spelled 7 ways, XS 6 ways, two unrelated `Provenance` classes (continuous_reference vs sood_registry). `ProblemSpec` Literals are violated at construction (`sphere`, `cylinder`, `pin-cell-2d`, `zero_flux`, `white_rank2`).
- The Q-R1 shape already exists one layer up: `MomentSpace`/`Billiard`/`Spectrum`/`BasisSpace` take `(StructuredGeometry, dict[int, Mixture], knobs)`. MMS cases carry the METHOD head (`quadrature`), inverting Q-R1.
- Census traps met: zsh does not word-split `$F` (use an array); `git grep -E` ignores `\b` silently (0 hits, rc 0) — use ugrep or an AST pass; AST attribute-read counts are receiver-unresolved (homonyms `psi`, `tolerance`).

Detail of the run lived in the session scratchpad `refgen_producers.md` (not durable). Related: [[census-predicates-bound-reference-and-activation-traceback]].
