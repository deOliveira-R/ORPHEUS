# Retired digest entries (cold; created by the 2026-11-16 distillation)

Verbatim text of digest entries and sub-points that left `lessons.md`, keyed by their stable ID.
Each is carried by the clause the digest's stub names; this file keeps the founding numbers.

## Retired 2026-09-22 (text from `git show 54809381^:.claude/agent-memory/numerics-investigator/lessons.md`)

### M3, M4 (old lines 31–42)

- **M3 — "Converged but wrong" is the SOLVER far more often than the operator.** An unread
  convergence flag reproduces the diverges-with-refinement fingerprint exactly (L10, L10b); a
  hardcoded `restart` truncates the subspace and a grown Krylov composite must re-derive it
  (L15b); the increment is ρ-blind (L11); the "hang" may be the fixture (L9); the iteration
  budget is dimension-dependent (L10c). **Read `converged`/`n_inner` against `max_*` BEFORE any
  tolerance sweep.**
- **M4 — A degenerate fixture annihilates a whole bug class; name the invariant your fixture
  shares with the class, then step outside it.** Flat ψ kills redistribution and routing (L2,
  L6); a uniform MATERIAL field kills an O(1) source/sink swap (L8·1); no trailing moment axis
  kills an axis-targeting bug (L13); a same-direction reduction kills a per-ordinate frame bug
  (L7·1); an even-`n` ladder kills a parity-split mode (L23·8); a zero-valued block kills the
  block's own rows (L17·2).

### M6, M7, M8 (old lines 50–70)

- **M6 — When two of my own measurements contradict, one instrument is wrong FOR THE QUESTION;
  find which directions it probes that production never supplies.** `‖MM⁻¹−I‖` full-space
  against a forward-substitution SUBSPACE inverse (L23·10); a nulp count against the relative
  error (L22·1); `‖Gv‖` against `‖Av‖` for a column dependency (L27·1); a signed projection
  against an approximate weight, worse than no weight (L20·2); which decomposition the library
  actually cuts on (L27·4).
- **M7 — A published claim carries the scope of the analysis that produced it. Evaluate the
  quantity, never the heuristic.** Read the paper's own stated approximation level before
  assuming a code bug (L5); a rule of thumb carries the ORDER it was derived at and its safe
  DIRECTION can invert (L25·5, L21·6); a signature is a consequence of a derivation's scope, not
  evidence about the joint problem (`vv-principles` Mode 8(3)); a literature diagnostic
  transfers only in its LEVEL-LOCAL form (L21·5); a docstring naming a theorem is a checkable
  claim (L19·2); "lumping" names ≥3 disagreeing operations (L25·6); a wrong attribution picks
  the wrong fix even when both fixes help (L4).
- **M8 — The reference is the first suspect, and this project's own corpus can be the
  contamination.** The theory page and the docstring documented the same wrong index (L8); a
  single quadrature can fake a structural win, so demand the same SIGNED error across ≥2 (L3),
  and a cross-family study carries a contamination guard (L16·3); a frozen reference goes stale
  by MAGNITUDE CLASS and the failure message cannot say which side is wrong (L22); gate every
  hand-built fixture on the mixture's own consistency identity before believing either side
  (L20·5).

### L2 (old lines 80–85)

## L2: Curvilinear redistribution is the prime suspect

DUPLICATE: `numerical-bug-signatures` Signature 1 (symptom, mechanism, per-ordinate probe,
catching test, why it hides) and `AGENT.md` Step 5. Diverging-with-refinement keff in
curvilinear geometry means the balance equation is wrong; the α-recursion and `ΔA/w` are where
it lives.

### L8 (full; its last two bullets return to the digest 2026-11-16) (old lines 164–179)

## L8: The project's own theory page can be the contaminated reference

ERR-063 (fission χ read at the SINK node where the SOURCE node belongs). The theory page's
equation AND the solver's docstring math documented the SAME wrong index, so confirming the code
against its own doc would have confirmed the bug. `vv-principles` #6 names reference
contamination and all its tells are *solvers*; this is the sharpening: **our corpus is a
reference too, and a doc plus its own docstring are ONE source, not two.** (`vv-principles`
#6's check since 2026-09-21.)

- Two independent grounds settled it: Hébert 2009 Eq. 3.57/3.58, and a sibling solver of a
  DIFFERENT kernel family that already source-indexed χ.
- **A masked O(1) bug needs a spatially VARYING field to surface** (M4 on a material field); the
  discriminating probe MUTATES a causally-irrelevant region's value and asserts invariance.
- **Grep every sibling ASSEMBLY module and classify each same-bug /
  different-structure-correct — the correct sibling is often the structural ground you need.**
  (A no-re-baseline claim needs a direct old-vs-new value comparison: `vv-principles` #12.)

### L9 (old lines 181–189)

## L9: A reported "hang" may be FIXTURE cost — bound the solver directly FIRST

#212: an "SN heterogeneous-keff non-convergence hang" was a fixture eagerly walking the whole
reference registry into O(minutes) of adaptive-mpmath solves; the solver converged in ~0.3 s.
**A timeout is not a convergence failure until the solver is bounded independently of its data
fixtures.** On a timeout in a test consuming a registry-walk or auto-discovery fixture, FIRST
bound the solver (tiny `max_outer`/`max_inner`, or bypass the fixture by calling the producing
module's builder); if the solver is fast, reproduce the slow path on a fast proxy and bisect
THAT, not the test. Record: `sn_keff_hang_was_eager_registry.md`.

### L12 (old lines 239–250)

## L12: Measure a localized defect's ORDER at the cell, not in the norm

DUPLICATE of the first two probes: `vv-principles` § Bit-identity carries "NEVER call an
OFFLINE-isolated error 'the floor' until it survives (1) an end-to-end swap, (2) a term-silent
control, (3) AMPLIFY", with the ERR-061 tell. The third probe is not there and is the one I keep
re-needing: ⭐ **the first-cell / boundary-row ORDER is the mechanism oracle; the global L2 is
not.** If the boundary cell is already O(h²) with the term off, the term cannot repair an O(h)
deficiency that is not there — and a volume-weighted L2 DILUTES a localized defect (`√V ~ h^1.5`
at the pole), so the order must be measured AT the cell (`vv-principles` #5's check since 2026-09-21).
Records:
`curvilinear_tau_clamp_vs_pole_floor.md`,
`issue_257_s9_ld_boundary_slope_optical_verdict.md`.

### L13 (old lines 252–256)

## L13: "Breaks only for LD / a spectator trailing axis" = a greedy `(Ellipsis, *idx)` index

→ now `numerical-bug-signatures` Signature 11 and ERR-087 (2026-09-22): the mechanism, the
discriminator, the probe, the blind classes and the catching test are the skill's, which this
agent preloads; this line is the pointer. Founding case #276 A2, `0b3275d`.

## Retired or stubbed 2026-11-16 (text from the digest at HEAD `0d68fc28`)

Also carried only here: L1's founding cost (six wrong hypotheses spent on cylindrical DD divergence by guessing before isolating).

### Digest header (with the f36572c8 recovery recipe) (lines 1–11)

# Numerics Investigator — Lessons

Behavioral digest: **what mistake was made, and what correction changed how I work.** Every
measured war story is COLD — six files under `_archive/` (moved verbatim 2026-09-21) plus the
records each entry names. Nothing here restates a rule or a preloaded skill; where a correction
has been uplifted, the entry cites the clause and stops.

**The L-numbers are STABLE identifiers, never renumbered** — `numerical-bug-signatures`
SKILL.md cites "lessons L10/L9" and "lesson L11" by number, and
the #282 probe `diag_282_sphere_repose_convergence.py:23` cited "L14/L15" (the probe is retired at `f36572c8` (R19); recover with `git show f36572c8^:<old path>`). A retired
entry keeps its number as a stub naming where its correction lives.

### Spine M1, M2, M5 (lines 13–37)

## The spine — three meta-lessons

Read these three every dispatch. Open an `L` entry when its instance is in front of you; open
`_archive/` only when a number needs checking.

- **M1 — Ask what KIND of question this is before measuring; the answer is usually a THEOREM, a
  SPECTRUM or a CLOSED FORM, not a re-run.** A rate question is a spectrum question (L19·1);
  "does this contract widen?" is *what commutes* (L25·1); a counting or kernel law is
  combinatorial, so derive it (L24·1); an ownership question is two discriminators (L26·1);
  "can statistic X gate contract Y?" is one number, the TRANSFER GAIN (L20·1); a labeling
  degeneracy is the operator's own symmetry group (L18·4). **Build the difference operator, not
  the operator** (L26·1). → now in `AGENT.md`, "Before the cascade" (2026-09-21).
- **M2 — A claim names a LIMIT; the separating axis is the OTHER discretisation going to zero,
  and your fixture is probably inside the claim's exact regime.** Angular consistency separates
  by `h→0`, not by the parameter it is named after (L21·1); a SEED is an angular closure, so
  sweep `N` at fixed fine mesh (L15a); every curvilinear MMS ansatz is ≤ linear-in-μ = the
  seed's exact regime, so the whole ladder is seed-blind (L14, L16·4). Enumerating the field
  family a closure is EXACT on, and gating outside it, is `vv-principles` Mode 7.
- **M5 — My own instrument has been wrong FIRST, repeatedly, always in the flattering direction.
  Budget every probe for its control.** The varying knob was unreachable (L18); a grep of pytest
  output read zero failures through ANSI codes and `python -c` imported the main tree instead of
  the worktree (L22·5); a probe "checked linearity" on a matrix that was linear by construction
  (L23·9); a docutils check reported no warnings on deliberately broken input (L27·6); a
  monkeypatch onto a `property` was a silent no-op (L26·8). This is `instrument-doctrine` X1
  applied to *analytic* instruments, and **each structural zero needs its OWN control.**

### Retired-IDs line of 2026-09-22 (lines 41–41)

_Retired 2026-09-22 by the agent-definitions audit, each restated by a clause or now carried by the definition: M3, M4, M6, M7, M8, L2, L8, L9, L12, L13. Numbers are not reused._

### L1 (lines 43–47)

## L1: Run the diagnostic cascade in order — no skipping

DUPLICATE: `AGENT.md` § Diagnostic Cascade ("Execute in order … Do NOT skip steps") and
`probe-cascade` § Anti-patterns ("Don't skip probes"). Founding cost kept as the teeth: **six
wrong hypotheses** were spent on cylindrical DD divergence by guessing before isolating.

### L5 table (lines 89–94)

| Observation | Likely cause |
| --- | --- |
| uniform offset across all parameters | code bug (constant, sign, factor) |
| scales with a PHYSICAL parameter (1/d, 1/τ) | paper floor (omitted higher order) |
| scales with a NUMERICAL parameter (n, dps) | quadrature or precision bound |
| insensitive to every numerical knob, exact at one structural limit | fixed approximation level |

### L6 (lines 102–108)

## L6: A curvilinear matvec is verified only against a NON-FLAT per-ordinate reference

DUPLICATE: `AGENT.md` Step 5 carries this as a standing rule and `vv-principles` Mode 7 the
enumerate-the-exact-regime check. Forensics kept: TWO cylinder matvec bugs (ERR-049 bool-mask
scatter ordering; the decoder "analytical extension" O(h) twin path) hid for months behind
flat-ψ-only coverage — `cyl_matvec_twin_path_signatures.md` (#197/#206). The architectural cure
is `coding-elegance` Pattern 2: ONE `SNCellOperator` consumed by sweep and matvec.

### L7·2, L7·3 (lines 119–124)

- **When every component is individually correct but the fixed point is wrong, build a
  structurally-independent from-scratch kernel** — reproducing the wrong value bit-for-bit
  localises the bug to the SHARED math.
- **When the literature says the scheme IS consistent and your faithful implementation is not,
  the bug is a CONVENTION between two correct pieces.** Record:
  `issue_240_d5b_s3_diffusion_limit.md`.

### L10b headline and first bullet (lines 134–143)

**L10b — the flag need not be DISCARDED: recorded-but-unread gives the identical false
fingerprint, and the tell is TOLERANCE-INSENSITIVITY.** A red gate was bit-identical across
`inner_tol ∈ {1e-9 … 1e-15}`, which I read as a real discretization floor; the residual at the
cap was ABOVE every tested tol, so all four runs hit the same `max_inner`.

- ⭐ **A tolerance sweep discriminates only if at least one tol is LOOSER than the residual the
  capped run reaches**; otherwise it is a vacuous knob (M5 in convergence clothing). **Read
  `history.converged` / `n_inner` against `max_*` FIRST** — a plateau at exactly `max_iter − 1`
  is the whole diagnosis (`numerical-bug-signatures` Signature 8 carries this since
  2026-09-21).

### L11 (lines 164–172)

## L11: For a ρ-honest stop or diagnostic, measure the residual `r = Aψ − q`, not `‖Δψ‖`

DUPLICATE: `AGENT.md` Step 3 carries the standing rule and `numerical-bug-signatures` Signature 9
the full symptom, mechanism and discriminator. **STALE API, landing named:** `FluxDisplacement`
was RETIRED 2026-08-19 (`coding-elegance` #18 records the reversal); the diagnostics now live on
the iteration record — `IterationRecord.increment_norms` / `.contraction_ratios` /
`.true_error_estimate()` in `orpheus/numerics/convergence.py`, with `AngularResidual` in
`orpheus/transport/residuals/angular_residual.py` carrying the per-ordinate balance map.
Catalogue: `issue_208_flux_displacement_residual_typing_debug_value.md`.

### L15 header and (a) (lines 190–200)

## L15: A SEED/closure re-pose is ruled PRINCIPLED-vs-REGRESSION by sweeping ANGULAR N, not `h`; and a grown Krylov composite must resize `restart`

`[LANDED a29ab2d]` (#282 route (a)). Cold:
`_archive/krylov_composite_restart_and_stale_references.md`.

- **(a) The N-sweep is the discriminator; the h-sweep at fixed N MISLEADS.** A seed IS an angular
  closure, so two treatments give different `keff(h→0)` at fixed angular order and the SAME
  `(h→0, N→∞)` answer — a gap that never closes under `h` reads FALSELY as a regression.
- **Honest nuance:** at the test's low N the OLD seed was *closer* to the N-converged truth. The
  justification was STRUCTURAL (an honest direct inverse, cold residual 5e5 → 1e-11), not
  accuracy. **Do not dress a structural win as an accuracy win.**

### L16·2 (lines 215–217)

2. **GATE the driver bit-faithful to production** on NON-FLAT vacuum spheres (flat or
   homogeneous-reflective is degenerate, M4) BEFORE swapping anything — else a driver bug
   masquerades as the effect under study.

### L17 dense-probe recipe (lines 240–244)

- **The dense-probe recipe**: assemble `A`, `apply_transpose` and the production metrics as dense
  matrices by unit-vector probing `to_flat`, and check `T == Aᵀ` against **numpy's** transpose (a
  structurally-independent ground, not the operator's own machinery). `T == Aᵀ` ⟹ gauge-free up to
  SPD and reciprocity cannot choose it (the general theorem is L26·4); if `T` carries a weight, the
  metric is PINNED.

### L18 (full, before the 2026-11-16 split) (lines 246–267)

## L18: Adjudicate a LABELING/ORDERING degeneracy with the operator's own SYMMETRY GROUP

#326 `[CLOSED; remediation LANDED dde93b64 — a level's order is the FIBER's.]` Cold:
`_archive/curvilinear_seed_metric_and_ordering.md`.

- **An MMS whose ansatz AND source depend only on the INVARIANTS of the degenerate class is
  EXACTLY blind** — declare the ansatz's invariants BEFORE trusting it as an adjudicator; a
  companion ansatz outside the symmetric sector SEES the defect and still cannot ADJUDICATE.
- **A within-tie permutation cannot move a cumulative-sum coefficient**, so the whole effect is a
  LABELING and the defect magnitude is ordering-invariant.
- **The leak path is any coupling that maps ACROSS the degenerate class, and it need not COMMUTE
  with the tie-break — grep every cross-class index map before calling a relabeling inert.**
- ⭐ **The adjudicating criterion is a SYMMETRY the continuous AND semi-discrete problems both
  have**: no reference solver, no MMS, structurally independent by construction. Verdict here —
  **no ordering is correct, the CLOSURE is broken.**
- **The constructive exit is to fold to the FUNDAMENTAL DOMAIN**: a degeneracy in a sort key
  usually means the discretization carries a redundant symmetric copy, and on the independent half
  the key is monotone and every competing criterion coincides.
- Two promotion-time corrections: **`xfail` swallows FIXTURE SETUP ERRORS too**, so pair the xfail
  row with a **reddenable un-xfailed SIBLING** on the same fixtures (sharpening `vv-principles`
  Mode 8(4)); and a diagnostic that RE-IMPLEMENTS the production kernel must be rewired to CALL it
  at promotion, mutating through the test's OWN import binding (Mode 11).

### L19 header and ·1 (lines 271–281)

## L19: An iterative-solver RATE question is a SPECTRUM question

#341 `[CLOSED; docstring repairs LANDED adc887d6; the octant-order lever is #343, OPEN]`. Cold:
`_archive/rate_spectrum_certificate_and_angular_axis.md`; record `issue_341_boundary_gs_rate.md`.

1. ⭐⭐ **Build the iteration matrix; do not re-time.** Any `x ← A_inv.apply(Σ gᵢ.apply(x))` driver
   IS a linear operator: wrap it over the composite's `to_flat`/`from_flat` as a
   `scipy.sparse.linalg.LinearOperator` and `eigs(which="LM")`. A few hundred sweeps buys the whole
   spectrum, and it is **immune to the stopping test**, so no ρ-blind stop (L11) or `max_inner`
   truncation (L10) can contaminate it. Controls (`instrument-doctrine` X1): reproduce the ρ FITTED
   from a real residual history, plus `G(2x) = 2G(x)` and `G(0) = 0`.

### L19 verdict discipline (lines 300–302)

**Verdict discipline that generalises: a production default must never branch on a variable you
have only CORRELATED.** `ndim` was falsified on both sides by direct measurement, and thickness,
mesh, aspect ratio, quadrature order and `c` all move the sign at fixed `ndim`.

### L20·1 (lines 309–316)

1. ⭐⭐ **The one number that decides any "can statistic X gate contract Y?" question is the TRANSFER
   GAIN `|Δy| / X`, measured across configurations — compute it FIRST, before any threshold hunt.**
   (The `instrument-doctrine` skill's X1 carries it since 2026-09-21.)
   A threshold on `X` bounds `|Δy|` only through that gain, so an unbounded gain means no constant
   exists and the tuning is void (`[M]` gain spread 1.16e+05×; a zero-false-alarm threshold missed
   15 of 16 corrupting cases). This is `vv-principles` Mode 12 read in the MIRROR — the gate is
   SIGHTED on a class the CONTRACT is blind to. **Cure: project onto the functional the contract
   reads.**

### L20·5 (lines 326–332)

5. ⭐ **Gate every verification fixture on the mixture's own consistency identity
   `σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` — an inconsistent mixture makes two legitimate references
   DISAGREE with no bug in either.** `[M]` the brief's "benign" reference was 30 % off because the
   fixture wrote `sig_s` as `[to, from]` while the constructor reads `[from, to]`, so the two solvers
   each honestly reported a different balance. **Never trust a brief's reference value on a
   hand-built mixture until the consistency identity is printed.** (`vv-testing` carries it since
   2026-09-21.)

### L21 header and ·1 (lines 336–344)

#319 / #235, 251 solves all `converged`. `[Phase-0 gates LANDED a3121cfe; #319 and #235 remain OPEN
for later phases.]` Cold: `_archive/rate_spectrum_certificate_and_angular_axis.md`; record
`issue_319_flux_dip_discriminator.md`.

1. ⭐⭐ **When a scheme claims consistency in the limit of ANOTHER discretisation, the axis that
   separates it from a rival is the OTHER mesh going to zero — because the claim is exact only
   there.** `[M]` sweeping optical thickness at fixed cells-per-mfp separated the two τ schemes **not
   at all** (fitted decay rate 0.000 for BOTH); refining `h` at fixed physics separated them without
   bound (3.2× → 204×). **Ask which limit the claim is exact in before choosing the sweep axis.**

### L22·1 (lines 379–383)

1. ⭐⭐ **A nulp count is uninterpretable in BOTH directions — report `max|a−b| / max|b|` FIRST.**
   (`vv-principles` bit-identity since 2026-09-21, with point 3.) The
   received warning is "huge nulp near zero is nothing"; the dual bites as hard. `[M]` `1.04e+15`
   nulp meant the values differ by **8 %** (216/216 elements), while an `array_equal → False` on two
   arrays that PRINT identically was **1 ULP**. One measurement re-sorted the whole investigation.

### L22·3 (lines 386–391)

3. ⭐⭐ **A rule-tier TOLERANCE pin can never warn about a consumer-tier BIT-IDENTITY pin — the gap is
   structural, not an oversight.** `gauss_legendre` is gated against numpy's `leggauss` at `< 8 ulp`
   (correct: neither construction is "the" answer), so a 3-ULP change is INSIDE that contract and
   OUTSIDE every downstream `array_equal`/`sha256` consumer simultaneously. **The only instrument
   that closes it is a byte-level fingerprint beside the tolerance pin**, whose sole job is "the
   bytes moved — re-baseline the consumers".

### L22·5 (lines 396–401)

5. ⚠ **Two probe-harness self-inflicted failures, both flattering** (M5): `grep -cE "^FAILED"` read
   `0` at nine commits including two already measured RED, because pytest's ANSI codes precede
   `FAILED` (use `--color=no`); and `python -c` prepends **CWD** to `sys.path[0]`, AHEAD of
   `PYTHONPATH`, so a worktree probe silently imported the MAIN tree and printed HEAD's values for
   every commit (run probe **script files** outside the repo, and make every probe print
   `module.__file__`; `code-search` carries both since 2026-09-21).

### L23 header and ·2 (lines 405–417)

## L23: A discrete operator's SINGULARITY is a two-question object — refuse the either/or

#344 `[CLOSED; the closed-form kernel and the exit gauge LANDED f934ff57/b51bc802; the
characterization gates LANDED 1a2be025.]` Cold, with every measurement:
`_archive/issue_344_singularity_kernel_and_gauge.md`.

1. ⭐⭐ **When a claim is "the null space IS class X", the decisive number is `dim ker` MINUS `|X|` —
   and you must MEASURE `|X|`, not assume it exists.** One line (`min|Ω·n|`) refuted the brief's
   framing before any solve, because a level-symmetric rule cannot produce `Ω·n = 0` at all; and
   `dim ker = |X| + R` then held exactly on 9 rows, which turned an either/or into a DECOMPOSITION.
2. ⭐⭐ **A dense SVD through the PRODUCTION builders is cheap and settles rank questions that ARPACK
   only BOUNDS** (the prior record's "3 at d=2, ≥6 at d=3" was an `eigs(k=12)` lower bound; the truth
   was 12 and 138). **Report the singular-value GAP** so the rank threshold is visibly not arbitrary.

### L23·9 (lines 446–453)

9. ⭐⭐ **To tell a GAUGE FREEDOM from an INCOHERENT solver, REMOVE THE KERNEL and re-run — and remove
   it ≥3 structurally-different ways**, because "boundary moves, bulk does not" is produced by BOTH.
   ⚠ **The obvious control can be a NON-control**: an even-`n_x` box, where the parity finding says
   the mode is absent, still has `dim ker A = 12` — what is absent is the kernel's EXCITATION by that
   source. **Assert `dim ker == 0` INSIDE the control; never infer it from a deviation being zero.**
   (My own probe was wrong first, twice, both times toward the alarming verdict — M5.) That a
   splitting's fixed point is a MANIFOLD when `A` is singular, so a schedule-invariance gate
   legitimately reds with no bug present, is now `vv-principles` Mode 9's false-RED premise.

### L23·10 (lines 454–459)

10. ⭐⭐ **`‖MM⁻¹−I‖` over the FULL space is the wrong instrument for an iteration's inverse — probe
    the RHS SUBSPACE the driver actually supplies.** `[M]` boundary-G-S reads 3.3e-01, which looks
    like incoherence and contradicts a measured `‖Aψ*−q‖/‖q‖ = 8e-14`; the defect sits exactly in the
    (inflow-row, outflow-column) block, where the driver supplies exactly zero content. **A reified
    forward-substitution "inverse" is a SUBSPACE inverse by construction** — fine for SI, a live
    hazard for a Krylov PRECONDITIONER (M6).

### L24·1 (lines 466–471)

1. ⭐⭐ **A counting law that does NOT depend on a parameter the operator plainly contains is telling
   you the governing equation is COMBINATORIAL — go derive it, do not fit it.** Setting the degenerate
   branch (`ψ_c = 0`) turns DD's closure into an involution, so every face field is a sawtooth;
   substitute into the balance and **every** cross-section, mesh width, weight and area cancels,
   leaving *a sum of functions, each blind to one coordinate and one sign, vanishing identically*.
   Both fitted laws then drop out as theorems. **The parameter-independence WAS the derivation hint.**

### L25·1 (lines 502–508)

1. ⭐⭐ **"Does accessor `f(x)` have to become `f(x, y)`?" is answered by asking what the defining
   conditions COMMUTE with — not by an expansion.** A **scalar convex combination commutes with every
   linear map**, so a quantity defined by membership of a convex set plus exactness of a scalar blend
   has both conditions as the SAME scalar statement in every component of every linear representation,
   and the widened form is an **overdetermined system whose every row returns the same value**. Three
   cheap hypotheses to check: the scalar is independent of the widened index; the projection is
   linear; the set is convex. Minutes, against a multi-day expansion.

### L26·1 (lines 542–545)

1. ⛔ **"Does X influence the adjoint?" is the WRONG question when the adjoint is a reverse-mode VJP** —
   `apply_transpose ≡ apply.T`, so the answer is "everything" and it carries no information. The
   informative question is the **symmetry character of the INCREMENT `∂A/∂k`**, measured by perturbing
   one entry and taking the dense difference. ⟹ **build the DIFFERENCE operator, not the operator.**

### L27·2, L27·3 (lines 603–611)

2. ⭐ **A threshold has TWO edges, and the instrument that justifies one is usually MONOTONE in it, hence
   blind to the other** (`vv-principles` #24(d), zero-set). `[M]` the Parseval statistic read a flat
   `1.000000000` across `[1e-15, 1e-2]` — a TRUE reading carrying **zero** information about the lower
   half. **Before citing a flat scan as justification, ask which direction the statistic can even move
   in.** (`vv-principles` #24(f) since 2026-09-21.)
3. ⭐ **Then go looking for the guard that already forecloses the other side — it is usually there and
   undocumented.** `[M]` every rcond below 8.696754e-17 is REFUSED at construction by a
   pair-consistency guard, so the corrupt band is **unreachable, not merely distant**. **Bisect the
   boundary; do not reason about it.**

### L27·6 (lines 623–626)

6. ⚠ **A doc claim's blast radius needs the RENDERING question first**: `[M]` the module carries no
   `automodule`, so `-W` cannot gate any of that prose at any severity. (My two instrument deaths that
   session — an unflushed docutils `warning_stream`, and a spliced-file import failing inside
   `dataclasses` because the module was not in `sys.modules` — are M5.)

### L28 header and ·2 (lines 628–639)

## L28: Platform drift is proven by INVERTING the frozen snapshots, then sorting gates by whether ULP jitter flips them

Not covered by `numerical-bug-signatures` Signature 10 (stale snapshot) or `vv-principles`
anti-pattern 31 (single-draw bit-exact): those name the classes, not the instrument for a
platform change with no old environment to run. `[M]` 2026-10-01, macOS 27.0.1 Accelerate.
1. **The old primitive output is recoverable from the snapshots it fed.** A coordinate ULP search
   over the few free numbers (GL-4: x1, x2, w1, w2) against the frozen walk-matvec `.npy` found a
   unique rule reproducing 496 of 496 elements; substituting it alone turned 20 of 26 reds green.
   That is a bit-exact positive control, which "it is probably LAPACK" never is.
2. **Jitter a correct input by ±3 ULP over 10 seeds** and count reds per gate: a gate red in some
   seeds and green in others is a noise-verdict gate (two FP association orders, `(a·y)/a == y`,
   a tolerance inside its own noise band), and re-capturing does not cure it.

