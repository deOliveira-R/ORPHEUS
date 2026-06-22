---
name: scan-march-verification
description: S5 #222 scan-march verification plan. Reframes the d-D DD sweep as schedule×backend; row-march+x-scan reuses ordinate_scan per line, unifying CumprodScan+the 2-D window. Wavefront-oracle-pinned principled-equivalence at nulp (NOT bit-identity across schedules); conditioning gate (DEMONSTRATED denormal-cumprod NaN production gap); Mode-9 FP-invariance off the degenerate box; ERR-056 shed-order; ERR-054 pole-reset survives.
metadata:
  type: project
---

S5 verification plan for the **scan-march unification** (issue #222) — the
headline of the SweepStrategy carve (`.claude/plans/sn_sweep_strategy.md` §S5).
Extends the S0 parent memo [[sweep-strategy-carve-verification]] (reuses its
L2/L3/L5/L6/L7/L8/L9 anchor vocabulary). Worktree `worktree-sn-nd-layout`,
HEAD as of 2026-06-10. The carve is an operator-algebra carve crossing a
subsystem boundary (the 2-D sweep schedule ↔ the 1-D `ordinate_scan` backend),
so this plan is the proactive `test-architect` gate before the main agent
writes the scan-march turn-by-turn.

## What the carve IS (the behaviour change the gates target)

The d-D DD spatial sweep is forward-substitution on a lower-triangular operator.
Today's production schedule is the **anti-diagonal wavefront**
(`MovingFrontierWindow` at d=2, `CumprodScan` at d=1). #222 reframes the within-
line DD recurrence as a **first-order linear scan** — the SAME `ordinate_scan`
the 1-D path already calls — so the d-D sweep becomes **scan along the sweep
axis, marched over the transverse axes**: `scan(x)` (d=1) → `scan(x)∘march(y)`
(d=2) → `scan(x)∘march(y,z)` (d=3). This unifies `CumprodScan` + the 2-D window
into ONE scan-march primitive. `FullFieldWavefront` stays the unconditionally-
stable ORACLE the scan-march is pinned against (the S1–S4 oracle-plus-fast-path
pattern). Two ORTHOGONAL axes: **schedule** (anti-diagonal / row-march) ×
**backend** (forward-sub / closed-form-scan / division-free pair-monoid — the
latter two already conditioning-dispatched inside `ordinate_scan`).

## ⭐⭐ CURRENT-BEHAVIOUR SURPRISE (DEMONSTRATED — load-bearing, changes the plan)

**#222's "latent gap" is a LIVE production correctness bug, reproduced.** The
`ordinate_scan` conditioning guard (`scan.py:160`) is
`if np.any(cumprod_a[-1] == 0.0): return _pair_monoid_scan(...)`. It catches an
EXACT zero, but a **denormal-but-nonzero** `cumprod` (where `b / cumprod`
overflows to `inf` and `cumprod · inf = NaN`) slips past the guard and
`ordinate_scan` returns NaN. The scan-march reuses `ordinate_scan` per line, so
it INHERITS this gap. Literal probe stdout (`-W ignore`,
`PYTHONPATH=<worktree>`):

```
=== DENORMAL-CUMPROD GAP (issue #222) — production ordinate_scan ===
chain: nx=312, a=0.1, b=1.0; cumprod[-1]=1.000e-312
cumprod[-1] == 0.0 (the guard test): False  -> guard does NOT fire
cumprod[-1] is denormal (0<|x|<tiny=2.23e-308): True
ordinate_scan all-finite: False  (NaN count=4)
explicit-loop reference all-finite: True, ref[-1]=1.111111
=> PRODUCTION BUG: guard catches EXACT zero, MISSES denormal -> NaN leaks

=== CONTRAST: exact-zero underflow (guard DOES fire, pair-monoid exact) ===
nx=800,a=0.2: cumprod[-1]=0.000e+00 exact_zero=True -> guard fires
ordinate_scan finite=True relerr=1.78e-16
```

(My first probe used `a=0.1, nx=800` which underflows to EXACT zero — guard
fires, no bug. The bug needs `cumprod[-1]` to land IN the denormal band
`(5e-324, 2.23e-308)` without hitting zero: `a=0.1, nx≈306–320`.) Whether the
SN sweep ever drives a single line into this band in production is a separate
question (the well-resolved DD regime keeps `|a|<1` but typically not to
`nx~312` of near-constant attenuation) — but it is reachable in principle, and
the scan-march MULTIPLIES the number of `ordinate_scan` calls (one per transverse
line × octant × group × iteration), widening the exposure. **Decision point for
the main agent (flag, do not silently fix):** either (i) widen the guard to
`np.any(~np.isfinite(b / cumprod_a))` / a denormal-floor test and route to the
pair-monoid, with a regression test pinning the fix; or (ii) explicitly scope it
out with a documented argument that the SN line-recurrence cannot reach the
denormal band, pinned by an `xfail(reason=...)` so the gap is recorded not
forgotten. **NEVER** ship the scan-march leaving this undocumented — a silent
NaN leak is the worst failure class (vv §H3: finite-looking telescoping hides
it). Gate G3.b below pins this regime regardless of which decision wins.

## ⭐ Architecture decision points (ENUMERATE consequences; main agent + user decide)

The plan does NOT pre-decide these; each fork has distinct verification
consequences listed.

**(a) New strategy `ScanMarch` vs `CumprodScan` BECOMES its d=1 instance.**
- *Fork A1 — new `ScanMarch` strategy, `CumprodScan` retained.* Consequence:
  the d=1 equivalence gate (G2 at d=1) must now pin THREE strategies pairwise at
  d=1 (`CumprodScan ≡ ScanMarch ≡ FullFieldWavefront`), or assert `ScanMarch`'s
  d=1 specialization IS `CumprodScan` (delegation, bit-identical). Concept count
  does NOT drop — Cardinal-Rule-2 payoff unrealized. REJECT unless there's a
  measured reason to keep two.
- *Fork A2 — `CumprodScan` retired INTO `ScanMarch`(d=1).* The #222 unification
  payoff. Consequence: `CumprodScan` deletion = TEST MIGRATION
  ([[feedback_retirement_means_test_migration]]) — every `CumprodScan(mesh)`
  reference in `test_wavefront_cumprod_equivalence.py`,
  `test_ordinate_scan_reset.py` solver-leg, the registry `SWEEP_STRATEGIES`, and
  `default_for` rewires to `ScanMarch`. The d=1 equivalence gate becomes
  `ScanMarch(d=1) ≡ FullFieldWavefront(d=1)` at nulp — same shape as today's
  `CumprodScan ≡ spine`. **Audit deliverable:** zero dangling `CumprodScan`
  refs after the carve (grep gate). This is the recommended fork; the verification
  cost is the migration, not new gates.

**(b) `ScanMarch` REPLACES `MovingFrontierWindow` as the d=2 production default
vs coexists.** This is a SEPARATE gated decision — do NOT assume the carve flips
the production default just because the prototype is 1.75× faster single-octant.
- *Fork B1 — coexist, `ScanMarch` opt-in.* `default_for` UNCHANGED at d=2
  (`MovingFrontierWindow` stays first in `SWEEP_STRATEGIES`); `ScanMarch` is
  selectable by tests + a future frontend. Consequence: the bit-identity anchors
  (G1: `test_affine_carve_bit_identity`, A2D-1 source-hash, d=2 window≡full) stay
  green UNTOUCHED — production path byte-for-byte unchanged. The scan-march is
  verified ONLY against the oracle (G2) + conditioning (G3) + Mode-9 (G4); the
  end-to-end solver gates (G6) run with `ScanMarch` forced via the selectable
  API, NOT as the default. LOWEST RISK. Recommended for the FIRST landing.
- *Fork B2 — `ScanMarch` becomes the d=2 default.* `default_for` reorders
  `SWEEP_STRATEGIES` so `ScanMarch` precedes `MovingFrontierWindow` at d=2.
  Consequence: `test_affine_carve_bit_identity` (G1) WILL shift (row-march vs
  anti-diagonal differ at FP-association ~ULP) → the sha256 golden must
  REGENERATE (principled-equivalence, NOT bit-identity — this is the explicit
  "do NOT demand bit-identity across schedules" boundary). The DD regression
  snapshots (G5) shift and regenerate. The A2D-1 source-hash regenerates IF the
  matvec body moves. This fork is a SECOND, separately-gated step — land B1
  first (oracle-pinned, default unchanged), flip to B2 only after G2/G3/G4 are
  green AND a measured end-to-end speedup justifies disturbing the golden.

**(c) Where the per-line `a_attenuation` coefficients get cached (#222 Q2
placement table).** `D = σ_t + Σ sₐ` and the closure multiplier `a` are
flux-independent + iteration-invariant. The 1-D scan already caches them in the
two-stratum cache (`CollisionCache.from_geometry`, `sweep_cache.py:438-449`,
`a_attenuation` on `_coll_cache`); the 2-D wavefront recomputes per cell per
level per sweep. The #222 placement table puts `sₐ` and the cached `a_attenuation`
on the MESH (memoised once per `(mesh, σ_t)`), the DD closure FORMULA on DD, the
per-iteration `b` on the SWEEP. Consequence for verification: a cache placed on
the mesh introduces a STATE object whose correctness needs a FOUNDATION pin —
`@pytest.mark.foundation` test that the cached `a_attenuation` (2-D, per-line)
equals the value the wavefront's inline `cell_kernel_batch` computes for the same
`(geom, σ_t)`, term-for-term (Pattern-2 single-source-of-truth: the cache and
the inline kernel must compute the SAME `a`). If the cache is per-line-recomputed
rather than mesh-memoised, this pin is cheaper but still required. **FLAG:** the
`(d−1)`-slab peak-memory win (#222 scope bullet) requires the coefficients be
computed PER LINE, not as a full-grid `D`/`a` array — Gate G7.b pins that the
peak memory stays `O((d−1)-slab)`, not `O(full grid)`, so a "materialize full-grid
D for convenience" regression is caught.

---

## THE GATE TABLE (headline)

| Gate | Level | Claim layer | Pillar / ground | What it asserts | Config (ndim, groups) |
| ---- | ----- | ----------- | --------------- | --------------- | --------------------- |
| **G1** anchor-preservation | L2/L3 | (all) | bit-identity inheritance | legacy production path byte-identical (Fork B1) | existing |
| **G2** scan-march ≡ oracle | foundation | flux-shape | `FullFieldWavefront` (→k_inf transitively) | `ScanMarch.sweep ≡ FFW.sweep` & `.residual ≡ .residual` at nulp | d=1 + d=2, 2G het aniso |
| **G3** conditioning/robustness | foundation/l0 | (numeric) | explicit-loop serial fold | scan-march finite+correct in underflow; denormal gap; `a=0` pole reset | scan-level + solver |
| **G4** Mode-9 FP-invariance | l1 | eigenvalue + flux-shape | k_inf=1.875 + SI≡Krylov | converged ψ*/k schedule-invariant off the degenerate box | d=2, ≥2G het aniso vacuum |
| **G5** DD regression snapshot | l1 | flux-shape | independently-corroborated snapshot | snapshots shift ~ULP, regenerate (Fork B2 only) | d=2 2G aniso het |
| **G6** end-to-end solver | l1 | eigenvalue | k_inf + SI≡Krylov | SI≡Krylov≡k_inf, keff_2d, fixed-source green w/ scan-march | d=2, ≥2G het |
| **G7** synthetic d=3 + memory | foundation | (structural) | window≡full synthetic-shape idiom | scan-march `scan(x)∘march(y,z)` admits d=3; peak mem = (d−1)-slab | synthetic d=3 (no quad) |

---

## 1. Anchor-preservation set (pins legacy behaviour; MUST stay green through the carve)

Classify each as **bit-identity** (free verification by inheritance), **principled-
equivalence** (nulp), or **value-ground** (structurally-independent reference).
Under **Fork B1 (default unchanged)** ALL of these stay green UNTOUCHED — that is
the whole point of landing the scan-march as opt-in first.

| Anchor | File::test | Class | WRAP-not-relocate? | Note |
| ------ | ---------- | ----- | ------------------ | ---- |
| **L2 PRIMARY** | `tests/sn/solve/test_affine_carve_bit_identity.py::test_converged_flux_bit_identical_after_affine_carve` (3 cases: `si_2d_p1_aniso_het`/`krylov_2d_p1_aniso_het`/`si_slab_2g_het`) | bit-identity (sha256 golden of converged ψ/φ) | n/a | Fork B1: stays byte-identical (default path untouched). Fork B2: REGENERATES (cross-schedule FP-assoc). `-O`-safe (`raise AssertionError`). |
| **L3 A2D-1** | `tests/sn/operators/test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin::test_apply_2d_cartesian_source_hash_unchanged` | bit-identity (sha256 of `inspect.getsource(_apply_2d_cartesian)`) | **YES — WRAP not relocate** | Pins LITERAL source text of the matvec body. If the scan-march matvec is a NEW method (`_apply_2d_cartesian_scanmarch`), `_apply_2d_cartesian` is UNTOUCHED → hash stays free-green. If the scan-march REPLACES the body of `_apply_2d_cartesian` (Fork B2 matvec), regenerate `EXPECTED_SHA256` + add a history line. |
| **L4 DD-regression non-square** | `tests/sn/regression/test_dd_regression.py` (the `2d_*_8x4_het` non-square cases) | principled-equiv (`assert_regression`, SAFETY×conv_tol≈1e-11; `2d_2g_p1_aniso_dd_8x4_het_si` pre-drifts ~6920 ULP) | n/a | The non-square 8×4 mesh catches x↔y swap. See G5. |
| **L5/L6 window≡full** | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py::test_sweep_window_equals_full_field_end_to_end` + `::test_matvec_window_equals_full_field_end_to_end`; `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py::test_solve_window_equals_full_field` + `::test_residual_window_equals_full_field` (d=1/d=2/d=3) | bit-identity (`np.testing.assert_array_equal`) | n/a | These pin window≡full WITHIN the wavefront schedule — UNTOUCHED by the scan-march (a different schedule). They are NOT the scan-march gate; G2 is. |
| **L7 φ=Q/Σ_t d=2** | `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py::test_2d_octant_sweep_closed_form_anchor` (`np.linalg.solve` 2×2) | value-ground (closed-form converged VALUE) | n/a | d=2 structurally-independent converged-value ground. STAYS (vv §1.5: ULP-distance necessary-never-sufficient). |
| **L8 k_inf=1.875 d=1** | `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py::test_cumprod_path_hits_analytical_kinf` (transfer-matrix `kinf_and_spectrum_homogeneous`) | value-ground (closed-form eigenvalue) | n/a | d=1 structurally-independent ground; the scan-march inherits it transitively via G2. |
| **L9 hand-loop x↔y moat** | `tests/sn/sweep_graph` `TestApplyMatchesLegacyInlined::test_per_cell_loop_equivalence` non-square (if present in worktree) + the non-square shapes in `SHAPES` of `test_sweep_graph_window_equivalence.py` ((12,7),(5,9)) | bit-identity | n/a | The x↔y-swap moat (Mode-2 detection). The scan-march's own moat is the d=2 non-square G2 config. |

**Anchor-set acceptance:** under `python -O`, the core anchor batch
(`test_wavefront_cumprod_equivalence` + `test_sweep_graph_window_equivalence` +
`test_affine_carve_bit_identity` + `test_ordinate_scan_reset`) is GREEN at the
current HEAD — verified this session (`PYTHONPATH=<worktree>`):

```
31 passed, 1 warning in 46.39s
```

(The 1 warning is pytest's `-O` advisory that bare `assert`s are ignored —
Mode-8 relevant, see G-discipline below. All selected items green.)

---

## 2. The scan-march ≡ oracle gate (G2) — the PRIMARY correctness gate

`FullFieldWavefront` IS the unconditionally-stable reference (the same DD physics
via a different, always-valid schedule, pinned transitively to k_inf=1.875 via L8
and to φ=Q/Σ_t via L7). **Do NOT propose a new oracle** — vv structural-
independence is satisfied: the scan-march and the wavefront are two valid
topological linearizations of the SAME lower-triangular solve.

**Test file:** extend `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py`
(or a sibling `test_scan_march_equivalence.py` if Fork A1 keeps both) — REUSE its
`_NULP_BOUND`, `_slab_sn_mesh`, `_seeded_inflow` helpers.

**G2.a — sweep equivalence, d=1.** `ScanMarch(mesh).sweep(Q, sig_t, bf) ≡
FullFieldWavefront(mesh).sweep(Q, sig_t, bf)` via
`np.testing.assert_array_almost_equal_nulp(ang_s, ang_m, nulp=_NULP_BOUND)` on
BOTH the angular and scalar outputs. At d=1, `ScanMarch` degenerates to the
`scan(x)` form (= `CumprodScan`'s body, Fork A2). Config: the existing
`_slab_sn_mesh(12, bc=..., ng_key="2g")` — mixture A 2g, GL-S8, heterogeneous
`sig_t` (`rng.uniform(0.3, 3.0, (ng, nx))`), non-uniform source, non-zero seeded
inflow per face. Parametrize `bc ∈ {"vacuum", "reflective"}`.

**G2.b — sweep equivalence, d=2.** Same assertion at d=2:
`ScanMarch(sm2).sweep ≡ FullFieldWavefront(sm2).sweep` at nulp. Config: a 2-D
Cartesian mesh that BREAKS degeneracy (point 4) — fuel|moderator split, ≥2G,
P1-anisotropic, level-symmetric S4 (cylinder/2D need level-structured quad — GL
raises; see S0 memo). Use the `_build_2d` shape from
`test_affine_carve_bit_identity.py` (8×4 fuel|mod, LS-4, 2G) so the config is
already a known non-flat het-aniso case. **MUST be NON-SQUARE (8×4 not 8×8)** —
the non-square mesh is the x↔y-swap moat (Mode-2): if the scan-march swaps the
scan axis with the march axis, the square mesh hides it, the non-square mesh
catches it.

**G2.c — matvec (residual) equivalence (L21: sweep & matvec are ONE operator).**
`ScanMarch.residual(operator, psi) ≡ FullFieldWavefront.residual(operator, psi)`
at nulp, d=1 and d=2. The matvec twin is the apply-direction scan-march; L21
demands sweep≡matvec discipline, so the residual must be pinned identically to
the sweep. This is the gate that prevents the scan-march matvec from drifting
from its sweep (the Phase-F twin-path failure class — `coding-elegance` Pattern 2).

**Nulp bound justification (decided).** Use the existing `_NULP_BOUND = 128` from
`test_wavefront_cumprod_equivalence.py`. Rationale (the file's own derivation,
vv §"Bit-identity vs principled-equivalence" criterion 3): a single-sweep FP
drift is bounded by `(reduction depth) × ULP`; the scan chain's reduction depth
is `nx ≤ 16` plus a ×8 safety factor (2-group coupling + source/BC affine terms)
→ 128. The scan-march at d=2 marches over `ny` lines but each line is an
independent `nx`-chain scan, so the per-element reduction depth is still `~nx`,
NOT `nx·ny` — the bound does NOT need to grow with the transverse extent.
**Decision: keep `nulp=128`.** If d=2 G2.b trips at 128 with a clean
FP-association story, investigate FIRST, then raise to `256` WITH a documented
`(reduction depth)×ULP` recomputation, NEVER silently. **Do NOT** raise the bound
to paper over a real divergence (anti-pattern `coding-elegance` #17).

**Anti-recommendation honored:** G2 is principled-equivalence (nulp), NOT bit-
identity. The scan-march and the wavefront differ at FP-association BY
CONSTRUCTION (row-major processes `(i-1,j)`,`(i,j-1)` before `(i,j)`; anti-
diagonal processes the same dependencies in a different order). Demanding
`np.array_equal` across schedules would be the WRONG gate.

---

## 3. Conditioning / robustness gate (G3)

The closed-form scan backend underflows for optically-thick / long chains
(`cumprod→0→NaN`, #222 Q3: nx=800 → NaN). The division-free pair-monoid backend
is exact and finite where the closed-form fails. The reference is the **explicit
serial-fold loop** (`_explicit_scan_loop` in `test_ordinate_scan_reset.py`) — a
structurally-independent ground (different algorithm, no division), already in
the worktree.

**G3.a — exact-zero underflow (EXISTING, inherited).**
`tests/sn/spatial/test_ordinate_scan_reset.py::TestOrdinateScanReset::test_ordinate_scan_multiple_and_consecutive_resets`
already pins `ordinate_scan` finite + correct at `a=0` (the pair-monoid path)
against the explicit loop. The scan-march reuses `ordinate_scan`, so this anchor
covers the per-line exact-zero case for free. ADD a scan-march-level wrapper that
drives a thick/long transverse line into the exact-zero regime and asserts the
2-D scan-march output is finite + matches the wavefront oracle (G2 config with an
optically-thick `sig_t` that forces `cumprod→0` on at least one line). Verified
contrast stdout (exact-zero → guard fires → pair-monoid exact):

```
nx=800,a=0.2: cumprod[-1]=0.000e+00 exact_zero=True -> guard fires
ordinate_scan finite=True relerr=1.78e-16
```

**G3.b — denormal-cumprod gap (NEW — the DEMONSTRATED production bug, §"surprise").**
A test that drives a single line into the denormal-but-nonzero band
(`cumprod[-1] ∈ (5e-324, 2.23e-308)`, `b/cumprod → inf`) and asserts finiteness.
Construction: `a = full(312, 0.1)`, `b = full(312, 1.0)`, `psi_0 = 1.0` →
`cumprod[-1] = 1e-312` (denormal, `!= 0.0`), guard MISSES, `ordinate_scan`
returns NaN. **This test FAILS against current production** — it is the
regression catcher for the fix (Decision (c) fork (i)), or an `xfail(strict=True,
reason="#222 denormal-cumprod guard gap — scoped out: SN line recurrence cannot
reach the denormal band, see <argument>")` if Decision fork (ii) is chosen.
Either way the gap is PINNED, not forgotten. Reuse `_explicit_scan_loop` as the
finite reference. Tag `@pytest.mark.catches("ERR-NNN")` once the bug gets a
catalog entry (see Self-Improvement below). Demonstrated stdout is in the
§"surprise" code fence above.

**G3.c — `a=0` cylindrical pole reset survives (EXISTING, ERR-054/#209).**
`test_ordinate_scan_reset.py::TestSICylinderResonance::test_si_agrees_with_krylov_at_resonance`
(`@catches("ERR-054")`) pins that the cyl pole-cell `a=0` resonance stays finite
(SI≡Krylov at the resonance). #222 confirmed the `a=0` pole reset SURVIVES the
"no internal boundary" geometry change (it is a geometric zero-area face at r=0,
not a BC). The scan-march is Cartesian-only at d≥2 (the curvilinear pole is a
d=1 `CumprodScan`/`ScanMarch(d=1)` concern), so this anchor stays on the d=1 path
— UNTOUCHED by the d=2 scan-march, but it MUST stay green (the d=1 scan-march, if
Fork A2 retires `CumprodScan`, inherits the cyl path → re-run this with the
`ScanMarch`-backed solver).

---

## 4. vv-principles Mode 9 — FP-invariance of a schedule change (G4)

The row-march schedule is an iteration-INVARIANT change: it MUST NOT move the
converged ψ*/k (only the per-sweep FP-association and the speed). Mode 9 demands
the FP-invariance be verified on a config that BREAKS the degenerate coincidence
— NOT the isotropic-reflective box (which makes the wrong formulation
accidentally exact).

**G4.a — converged value ≡ wavefront-schedule fixed point.** Drive a full
`solve_sn` / `solve_sn_fixed_source` with the scan-march as the inner sweep, and
assert the converged ψ*/k equals the `MovingFrontierWindow`-schedule fixed point
to SOLVER tolerance (NOT nulp — the converged values agree to `SAFETY×conv_tol`,
the schedules differ only in transient FP-association which the outer iteration
washes out). Config (BREAKS degeneracy, per the anti-recommendation): **2-D
Cartesian, ≥2G, fuel|moderator HETEROGENEOUS, P1-ANISOTROPIC, VACUUM (streaming)
on at least one axis pair** — exactly the `si_2d_p1_aniso_het` config from
`test_affine_carve_bit_identity` (8×4, vacuum-x / reflective-y, LS-4, 2G, P1).
The anisotropic + heterogeneous + vacuum drives the flux genuinely non-flat
(redistribution out of cancellation), so a schedule that secretly changes the
operator (not just the FP-order) is caught. Assert `abs(k_scanmarch - k_window)
< 1e-7` AND `assert_allclose(phi_scanmarch_normalized, phi_window_normalized,
rtol=1e-6, atol=1e-8)` (flux SHAPE, mean-normalized — eigenvectors are scale-
free), mirroring `test_si_krylov_heterogeneous_2g_nonflat_flux`. INCLUDE the
non-flat degenerate-gate guard `prof.max()/prof.min() > 1.2` so the test cannot
pass vacuously on an accidentally-flat flux.

**G4.b — ERR-056 shed-order concern (the scan-march's NEW reflective-BC risk).**
The scan-march changes HOW reflective-BC outflow is shed: x-outflow = the LAST
scan value on each line, y-outflow = the LAST row's ψy (vs the wavefront's
per-octant `capture_x`/`capture_y` shed at `sweep.py:956-957`). ERR-056 is the
documented failure class where the OUTFLOW SHED ORDER on a shared face is wrong
for a diagonal/non-axis-aligned cubature (each shared face touched by >1 octant).
**Design a reflective-BC FP-invariance test on a config where shed-order matters:**
a 2-D Cartesian mesh with REFLECTIVE BCs on BOTH axis pairs (so every octant's
outflow on every face is reflected back as another octant's inflow → the shed
order is load-bearing) and a DIAGONAL cubature (`level_symmetric` — O_h, shared
faces; NOT an axis-aligned `product` where each face is one octant). Assert the
scan-march converged ψ* ≡ the wavefront converged ψ* to `SAFETY×conv_tol`.
**LIMITATION to STATE in the test docstring:** in pure 2-D Cartesian with
`level_symmetric(4)`, the in-plane octants are the 4 sign-quadrants and the
shared-face coupling is the reflective inflow=outflow at each axis edge; a
genuinely diagonal cubature with off-axis ordinates sharing a face is more
naturally a 3-D / curvilinear stressor. If the 2-D Cartesian config cannot
realize a true shared-face shed-order ambiguity, STATE that the ERR-056 stressor
is only fully exercised at d=3 (deferred, no 3-D quadrature) and that G4.b at
d=2 pins the reflective shed-order to the extent 2-D Cartesian admits it — do NOT
claim coverage the config cannot deliver.

**Anti-recommendation honored:** G4 runs on anisotropic + heterogeneous +
vacuum/streaming, NEVER the isotropic-reflective box. The eigenvalue leg is
verified via the k_inf analytical anchor + SI≡Krylov, NOT via MMS (MMS is
source-driven, proves zero eigenvalue information).

---

## 5. DD regression snapshot policy (G5)

The scan-march differs from the wavefront left-fold at FP-association (~ULP),
exactly as the 1-D `CumprodScan` already does vs the spine. This is principled-
equivalence, NOT bit-identity — the explicit "do NOT demand bit-identity across
schedules" boundary.

- **Fork B1 (default unchanged):** NO snapshot shifts. `test_dd_regression.py`
  stays green untouched (production path is still the wavefront). The scan-march
  is verified by G2/G3/G4, not by the snapshots.
- **Fork B2 (scan-march becomes the d=2 default):** the d=2 snapshots
  (`2d_1g_dd_15x15`, `2d_2g_LS4_dd_8x4_het_si`, `2d_2g_p1_aniso_dd_8x4_het_si`)
  SHIFT by ~ULP (the row-march FP-association) and REGENERATE via
  `python -m tests.sn.regression._generate_snapshots`. The bound is
  `assert_regression(kind="iterative", conv_tol=...)` → `SAFETY × conv_tol`
  (≈1e-11); the regenerated snapshots stay within it. **Regeneration discipline:**
  before regenerating, the snapshot's converged value MUST be independently
  corroborated (the generator's own narrative requires it: 2G→k_inf=1.875,
  het k_eff against a cross-method, P1 flux against the closed-form flat anchor)
  — so the regeneration is NOT a tautology (re-pinning a possibly-wrong value).
  The `2d_2g_p1_aniso_dd_8x4_het_si` snapshot already pre-drifts ~6920 ULP /
  ~9.8e-13 (Phase-5b/5c inheritance) — well within the ≈1e-11 gate; the scan-
  march adds another schedule-FP-assoc delta of the same magnitude, still within
  tol. **STATE in the commit:** which snapshots regenerated, the measured ULP
  shift, and that the shift is `≤ iteration_count × ULP` (FP-non-associativity,
  dimensionally explainable — vv criterion 3). DO NOT escalate the gate to
  `-W error::DriftWarning` strict for these d=2 cases after a B2 flip (they
  legitimately drift); keep the 1-D / vacuum cases strict where they DON'T drift.

---

## 6. End-to-end solver gates (G6)

These run with the scan-march as the inner sweep — as the DEFAULT (Fork B2) or as
the OPT-IN selectable strategy (Fork B1). **The plan does NOT decide which; the
main agent + user resolve it (decision point (b)).** Under Fork B1 the gates run
with `ScanMarch` forced via the selectable API; under Fork B2 they run as the
default. Either way they MUST stay green.

- **G6.a — SI ≡ Krylov ≡ k_inf (≥2G, the anti-degeneracy gate).**
  `tests/sn/eigenvalue/test_keff_2d.py::test_si_krylov_heterogeneous_2g_nonflat_flux`
  (8×4 fuel|mod, LS-4, 2G, non-flat guarded) — the scan-march inner must reproduce
  SI≡Krylov flux SHAPE agreement + k_inf. **NEVER a 1-group test** (L2: k is
  flux-shape-independent at 1G — the 1-group-degeneracy bar). The k_inf anchor
  (`test_default_entry_hits_kinf`, `test_2g_eigenvector`) supplies the
  structurally-independent eigenvalue ground.
- **G6.b — eigenvalue (`test_keff_2d`).** The full `test_keff_2d` class
  (`test_homogeneous_exact` ≥2G, `test_2g_eigenvector`, `test_default_entry_hits_kinf`,
  `test_si_2d_keff_converges_under_refinement`,
  `test_eigenvalue_jacobi_gauss_seidel_equivalence`) green with the scan-march
  inner.
- **G6.c — fixed-source.** `solve_sn_fixed_source` cases (the
  `test_affine_carve_bit_identity` configs, `test_2d_anisotropic_windowing`,
  `test_fixed_source_2d_equivalence`) green; under Fork B2 the affine-carve golden
  REGENERATES (see G1/G5), under Fork B1 it stays bit-identical.
- **Deselect** `tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`
  (#212 `continuous_get` hang). Held reds untouched by the d=2 scan-march: #206
  cyl-matvec, #195 MMS@160.

---

## 7. Synthetic d=3 march admission (G7)

Prove `scan(x)∘march(y,z)` is genuinely d-general WITHOUT a 3-D quadrature — the
B7-style synthetic-shape idiom already used by
`test_sweep_graph_window_equivalence.py` (`SHAPES` includes `(3,2,3)`,`(4,3,2)`
d=3 synthetic) and `test_sweep_graph_nd_admission`. Separate from the (deferred,
no-3-D-mesh) PERF question.

**G7.a — structure admission.** On a synthetic d=3 shape (arbitrary `N_oct`,
streaming, source — NO real quadrature; the `_inputs(shape, N_oct, ng, seed)`
idiom from `test_sweep_graph_window_equivalence.py`), assert the scan-march
`scan(x)∘march(y,z)` walk produces a finite result that equals the
`FullFieldWavefront` (d-generic spine) `apply` on the same synthetic inputs, via
`np.testing.assert_array_almost_equal_nulp` (cross-schedule → nulp, NOT
`array_equal`). This proves the recursive scan-march construction is d-general:
the `(y,z)`-plane gets its own march schedule. **The d=3 CONTIGUITY / SPEED is
OUT of the correctness gate** — it is the one measured-cost question (#222 scope:
the `(d−1)`-slab memory win holds, the d=3 simplex-surface contiguity/speedup is
profiling, deferred).

**G7.b — peak-memory `(d−1)`-slab invariant.** Mirror
`test_sweep_graph_window_equivalence.py::test_window_backing_is_the_d_minus_1_slab`:
assert the scan-march's per-line coefficient working set is `O((d−1)-slab)`
(`∏_{a<d−1} n_a`), INDEPENDENT of the determined (scan) axis extent. This catches
the "materialize full-grid `D`/`a` for convenience" regression (#222 scope bullet:
compute `D`/`a` PER LINE). If the coefficients are mesh-cached (decision (c)), the
cache is per-`(mesh,σ_t)` and shared across iterations — the per-SWEEP working set
is still the `(d−1)`-slab `b` plus the rolling frontier; pin THAT.

---

## Discipline / G-discipline (Mode 8 — `-O` strip)

EVERY always-on gate (G2/G3/G4/G7 foundation+l0/l1) MUST use
`np.testing.assert_*` / `pytest.fail` / `raise AssertionError`, NEVER a bare
`assert` — the canonical ORPHEUS invocation is `python -O`, which strips bare
`assert` to a NO-OP (Mode 8). The existing anchor files already follow this
(`test_wavefront_cumprod_equivalence.py`, `test_sweep_graph_window_equivalence.py`
docstrings call it out explicitly). ⚠ The S0 memo flagged a LIVE Mode-8 hazard in
`test_unified_sweep_dispatch.py` (bare-assert) — that was MIGRATED in S1; do NOT
re-introduce bare `assert` in any scan-march gate. The `-O` advisory warning in
the anchor-batch run above ("assertions ignored") is pytest noting that SOME
test module somewhere uses bare assert — confirm the scan-march gates are NOT
among them (grep `^\s*assert ` in the new test files; expect zero).

## Self-Improvement triggers (fire BEFORE delivering)

- **NEW failure mode → skill update.** The denormal-cumprod NaN gap (G3.b) is a
  variant of `numerical-bug-signatures` Signature 1 / vv Mode 3 (missing-factor /
  guard-incompleteness) but with a NEW mechanism: a conditioning GUARD that tests
  the wrong predicate (`== 0.0` instead of `~isfinite(b/cumprod)`), so a denormal
  slips through. This is NOT yet in the `vv-principles` failure-mode table nor
  `numerical-bug-signatures`. **Recommended:** when the bug is fixed (decision (c)
  fork (i)), open ERR-NNN in `error_catalog.md` (failure mode: guard-predicate-
  incompleteness; how it hid: the exact-zero pole-reset test passed, the denormal
  band was never probed; catching test: G3.b; lesson: "a conditioning guard must
  test the failure CONDITION (`~isfinite(b/cumprod)`), not a PROXY (`cumprod==0`)
  — the proxy misses the denormal band"). If scoped out (fork ii), it stays a
  skill-table risk row, NOT an ERR (no production bug shipped). NOT appended to
  the skill yet — the disposition (fix vs scope-out) is the main agent's decision
  (c); the ERR/skill update fires at fix-time per "log every caught bug".

## Pre-reads cross-check (file:line, verified this session)

- `orpheus/sn/spatial/scan.py:160` — the `np.any(cumprod_a[-1] == 0.0)` guard (the
  gap's home).
- `orpheus/sn/sweep_cache.py:438-449` — `a_attenuation` two-stratum cache (#222 Q2).
- `orpheus/sn/sweep.py:956-957` — `boundary_flux.face_view(x_out_face)[oct_idx] =
  capture_x` (the wavefront per-octant outflow shed; the ERR-056-relevant seam the
  scan-march replaces with last-scan-value / last-row-ψy).
- `orpheus/sn/sweep_graph.py:849` (`apply_windowed`) / `:929` (`residual_windowed`)
  — the anti-diagonal wavefront walk the scan-march is an alternative schedule to.
- `orpheus/sn/sweep_strategy.py:488` — `SWEEP_STRATEGIES` registry / `default_for`
  (where a `ScanMarch` slots in; selection priority decides Fork B1 vs B2).
- `default_for` confirmed this session: d=1 → `CumprodScan` (`{CumprodScan,
  FullFieldWavefront}` applicable); d=2 → `MovingFrontierWindow` (`{MovingFrontierWindow,
  FullFieldWavefront}`); cyl-1D → `CumprodScan` only.
