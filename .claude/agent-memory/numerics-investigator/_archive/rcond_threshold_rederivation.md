# Archive — the _DENSE_METRIC_RCOND threshold re-derivation (ERR-080 / #429)

Moved COLD from `lessons.md` 2026-09-21 (verbatim, lines 1168-1228).
Digest successor: L27. Record: `scratch/rcond_rederivation.md` + `scratch/probe_rcond_{01..17}_*.py` (UNTRACKED).

## L27: A "noise mode" reading of a small singular value is a HYPOTHESIS — kill it with the closed form and the RAW TABLE; and a numerical THRESHOLD has two edges, one of which its justifying instrument is structurally blind to

`_DENSE_METRIC_RCOND` re-derivation (2026-08-31, ERR-080/#429). A pinned
`pinv` cutoff was justified by *"one `~1e-16` **noise** mode … `1e-12` sits ~4
orders above the noise floor"*. Every quoted number reproduced; the reading was
inverted. Six transferable points, in order of leverage.

1. ⭐⭐ **Discriminate "round-off residue of an EXACT dependency" from "small
   real mode" by SOLVING for the dependency, then applying the candidate null
   vector to the RAW TABLE — not to the Gram.** Forming `G = AᵀWA` can itself
   manufacture rank loss, so `‖Gv‖ ≈ 0` is weak; `‖Av‖ ≈ 0` says the *columns*
   are dependent as functions on the node set and settles it. `[M]` here the
   closed form fell out in two lines (`Y_2^{+2} = (√3/2)(1−μ²)`,
   `1−μ² = (2/3)(Y_00−Y_20)` ⟹ `v = (1,0,−1,0,−√3)/√5`), matched the SVD null
   vector to `2.2e-16`, and gave `‖Av‖∞ = 1.4e-16`. ⟹ **a "noise floor" is a
   claim about a number that should have a NAME; if it has one, it is not
   noise.**
2. ⭐ **A threshold has TWO edges and the instrument that justifies one is
   usually MONOTONE in it, hence blind to the other** (vv #24(d), zero-set).
   Parseval here is `ratio = 1 − (cᵀG_dropped c)/(cᵀGc)`: *lowering* rcond only
   shrinks `G_dropped`, so a flat `1.000000000` across `[1e-15, 1e-2]` is a
   true reading that carries **zero** information about the lower half. Before
   citing a flat scan as justification, ask which direction the statistic can
   even move in.
3. ⭐ **Then go looking for the guard that already forecloses the other side —
   it is usually there and undocumented.** `[M]` every rcond below `8.696754e-17`
   is REFUSED at construction by the *pair-consistency* guard
   (`_DENSE_METRIC_PENROSE_RTOL`, `|G⁺GG⁺−G⁺| = 7.9e+16` at `1e-18`), so the
   corrupt band is **unreachable, not merely distant**, and `‖G⁺‖₂` is
   bit-constant across the whole legal window. Bisect the boundary; do not
   reason about it. (The same guard is *blind* upward — `G⁺GG⁺ = G⁺` holds for
   a truncated pinv too.)
4. ⚠ **Read which DECOMPOSITION the library actually cuts on, and re-measure
   the residue several ways.** `pinv(hermitian=True)` cuts on `eigh`; the
   shipped comment quoted `svd`. `[M]` the SAME matrix's largest residue reads
   `8.70e-17` (eigh) / `9.71e-18` (eigvalsh, svd-hermitian) / `2.27e-17` (svd)
   — **9.0× spread** — and the bisected refusal boundary equals the *eigh*
   figure to 7 s.f. A round-off number has no stable value; that, not distance,
   is the real argument for a wide margin.
5. ⭐⭐ **A threshold validated on ONE flagship fixture is a claim about that
   fixture's spectral GAP. Census the shipped grid.** `[M]` at the shipped
   rcond: **31 of 105** slab `(order, L)` rows affected (20 raise, 11 breach the
   Parseval gate's own `rtol`), min affected `L = 3`, including the DEFAULT
   `gauss_legendre(16).angular_frame(4)`; **0 of 196** 3-D rows. And the
   mechanism was the opposite of the expected one: **the LIVE modes descend to
   meet the pin**, they are not met by a rising kernel — a 1-D chart mints
   `~L²/2` phantom columns against a fixed node count, the odd-`m` ones carry a
   non-polynomial `√(1−μ²)`, and the spectrum fills in continuously. Where there
   is no gap, **no cutoff is right** and the repair is upstream.
6. ⚠ **Method, twice-bitten: my own instruments died silently and a positive
   control caught both.** (a) A docutils RST check reported "no warnings" on
   deliberately broken input — an unflushed `warning_stream`. (b) A spliced-file
   import raised inside `dataclasses` because the module was not registered in
   `sys.modules` — **my harness, not the candidate**; reporting it as a defect
   would have impeached a correct block (plan-authoring §4's VERIFY clause).
   Also useful: sizing the blast radius of a doc claim needs the *rendering*
   question first — `[M]` `metric.py` has **no `automodule`**, so `-W` cannot
   gate any of this prose at any severity.

Full record + 17 probes: `scratch/rcond_rederivation.md`,
`scratch/probe_rcond_{01..17}_*.py`, block `scratch/rcond_comment_block.txt`.
