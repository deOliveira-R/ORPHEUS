# Issue #123 — Single-quadrature signal is not closure quality

**The bug class (not a single bug).** Multiple proposed rank-N
closure "improvements" (Direction-C, Direction-Q, Direction-N)
appeared to beat the F.4 baseline at the optically-thick regime
(σ_t·R ≥ 10) on a single high-grade quadrature
(RICH = 4×8×64). When re-run on a structurally different
quadrature (RICH+panels = 5×8×64), the *signed* error of the
closure flipped — and the signed error of the F.4 baseline
itself flipped at a different quadrature-refinement point. The
original "win" was an artefact of the two errors crossing zero
on opposite sides of the chosen quadrature snapshot.

**Evidence that existed before each false win was caught.** A
single-quadrature scan over the canonical 6-point grid showed
each candidate closure beating F.4 in absolute error at every
point. Convergence under mode count looked monotone. The closure
operators were derived from clean physical reasoning (Marshak
DP_N, PCA basis reduction, normalised modal projection). The
pattern was: physically motivated candidate + numerical evidence
on one high-grade quadrature + monotone convergence trace =
"structural improvement."

**Why that evidence didn't catch it.** A single quadrature gives
quadrature noise that is correlated between closure and reference
at the magnitude where rank-N improvements live (≤ 0.01 % k_eff
at the canonical grid). The two errors share an integration grid,
so quadrature artefacts cancel partially in a way that depends
on quadrature choice; switching the quadrature exposes which
sign was real and which was noise. The "L17/L19 quadrature
crossing" pathology is the named failure: closure error and F.4
error cross zero at different quadrature-refinement points, so
any single snapshot picks one sign for each and reports a
structural pattern that does not survive refinement.

**What evidence class would have caught it.** Structural
improvements to closures (vs quadrature noise) are
distinguishable **ONLY** by signed-error stability across
structurally different quadratures. **NEVER** ship a rank-N
closure claim on a single quadrature snapshot — the claim is
statistically degenerate at σ_t·R ≥ 10. **Instead**, require
≥ 2 quadratures (e.g. Gauss-Legendre vs tanh-sinh, or RICH vs
RICH+panels) that show the *same signed error* — not just the
same magnitude — for both the candidate closure AND the F.4
reference. The five L19 conditions in
`assert_rank_n_structural_win` operationalise this gate: strict
beat at every quadrature, no closure sign flip across
quadratures, monotone non-increasing magnitude, AND the F.4
reference itself sign-stable across the same quadratures
(otherwise the reference is in its own quadrature-noise regime
and the comparison is unverifiable). See vv-principles
reference.md §4.4 (semi-analytical correctness ladder, where
quadrature-stability is part of "integrator correctness").

**References.** Issue #123 (open at devcontainer 120 s/run
budget); Issue #121 (closed, Direction-C falsified by this
gate); Issue #122 (closed, Direction-Q falsified by this gate);
numerics-investigator agent memory L3 (lessons.md) and
`direction_n_quadrature_baseline.md`;
`tests/cp/test_peierls_rank_n_protocol.py::assert_rank_n_structural_win`.
