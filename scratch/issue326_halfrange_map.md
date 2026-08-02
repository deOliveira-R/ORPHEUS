# Issue #326 — what a HALF-RANGE azimuthal level would touch

**Read-only investigation.** Branch `refactor/operator-strategy-layers` @ `b0a003b4`.

⚠ **The tree MOVED during this investigation (L-007/L-012).** At open, no test file
was dirty; at close, `tests/sn/_test_helpers.py` is **`+118` uncommitted** and three
`#326` test modules are **UNTRACKED and in-flight** — the main session is building
the #326 gate scaffolding concurrently. Every `orpheus/` line reference below was
**re-verified by grep at final read** and is stable; the `tests/` references are
tagged in-flight where they are.

* **No production file was modified** at any point (`git diff --name-only` over
  tracked files: `.claude/*` + `tests/sn/_test_helpers.py` only).
* In-flight, untracked at close: `tests/sn/sweep/curvilinear/test_alpha_closed_form.py`,
  `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py`,
  `tests/sn/verification/mms/test_mms_ordering_blindness.py`,
  `orpheus/numerics/roots_of_unity.py` + `tests/numerics/test_roots_of_unity.py` (#325).

Claims marked `[M]` are **measured** this session with `.venv/bin/python`.
Claims marked `[D]` are **doc/docstring assertions** — reported as CLAIMS and
separately adjudicated against the implementation (`vv-principles` L33).

---

## HEADLINE — the hypothesis SURVIVES, with exactly one genuine break

**The redundancy IS the bug.** Everything in the box holds, and the "ordering
question dissolves" claim is now `[M]`-confirmed on the actual τ machinery. The
sweep/α/τ half of the hypothesis has **no** counter-example.

**The one thing that legitimately needs the full circle is the MOMENT INTEGRAL,
and only its ξ-ODD slots.** `[M]` A naive "halve the nodes, double the weights"
fold reproduces every ξ-EVEN spherical-harmonic moment to 5e-16 and turns the
ξ-ODD moment `φ_1^{ξ}` from its structural `−1.3e-16` into **`+2.94`** — garbage.
Cylindrical P1 scattering is live and tested
(`tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:219-264` sums
**all three** `Y_1^m` slots), so this is a real consumer, not a hypothetical.

The correct statement is not "double the weights" but: **the fold is a restriction
to the ξ-even sector, so the harmonic TRIAL SPACE must be restricted to its ξ-even
subspace too.** The ξ-odd harmonics are not "computed wrongly" on the folded set —
they are *out of the space*. That converts the break from a defect into a typed
design obligation (and ORPHEUS already has the vocabulary: a folded analysis face
with an unfolded reconstruction face is exactly `PetrovGalerkinFrame`; the
symmetrized basis restores `GalerkinFrame`).

**Second-order finding that decides the candidate home:** there are TWO distinct
"half range" constructions and they are NOT interchangeable —

* **(A) fold the existing nodes** onto `ω ∈ [0, π]` (endpoints `ω = 0, π`
  included, trapezoid weights `[w, 2w, …, 2w, w]`) — the construction the tree's
  OWN test already uses (`tests/sn/sweep/curvilinear/test_alpha_closed_form.py:255-268`);
* **(B) Hébert's midpoint half-range**, `2ℱ(p)` nodes strictly inside `(0, π)`.

`[M]` **Both** reproduce the α closed form with the same constant `2κ`
(2.052344 at `n_phi = 8`). What separates them is the **R12a seed predicate**:

| | `τ_raw` on level 0, `n_mu=4, n_phi=8` | `0 < τ_raw[0] < 1` (R12a) |
|---|---|---|
| CURRENT (full circle) | `[0, 1, 0, 1, 0, 1, 0, 1]` | **False** — non-carrying |
| **(A)** fold-existing | `[0, 0.2929, 0.5, 0.7071, 1]` | **False** — non-carrying |
| **(B)** Hébert midpoint | `[0.2195, 0.4142, 0.5858, 0.7805]` | **True** — CARRYING |

⚠ **(B) flips every cylindrical level onto the `radial_characteristic_*` ψ½
route**, which today is exercised by the sphere-GL rule ONLY, and simultaneously
retires the `#280 2.5b` direct-seed diagonal fold for the cylinder. That is a
large, live structural consequence — the single most important thing to decide
before choosing a home. **(A) leaves the seed trichotomy exactly where it is.**

---

## 0. Verification of the brief's premise (Operating Principle 5)

The brief's premise is **current**. The sharpening: part of the work is already
SCAFFOLDED — but **in the WORKING TREE, not at HEAD** (⚠ see the header; these
three modules are untracked and were being written while I audited).

* `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py` *(untracked,
  in-flight)* — carries **three `xfail(strict=True)` rows** documenting the live
  ψ-not-even-in-ξ defect, and its docstring at `:66` already states *"The
  constructive exit is the half-range level (Hebert §3.9.3), under which only one
  member of each pair exists and the symmetry holds by construction."*
* `tests/sn/sweep/curvilinear/test_alpha_closed_form.py:255` *(untracked,
  in-flight)* — a `@pytest.mark.foundation` section **§C "on the HALF range the
  two criteria COINCIDE"**, whose `_half_range_level` helper is construction **(A)**.
* `tests/sn/_test_helpers.py:978` *(uncommitted `+118`)* —
  `PRODUCT_LEVEL_ORDERINGS` + `_product_rule_with_ordering(tie_break,
  exact_nodes=)`, the shared ordering-swap harness.

So the deliverable is **not** "discover the fix" — the fix is named and its α
criterion is already pinned green *in the working tree*. The deliverable is the
**blast radius + home selection**, which is what follows. **Nothing about the fix
is in the theory corpus or on `main` yet** (§6 D4).

---

## 1. The seam inventory — MEASURE vs LEVEL STRUCTURE

### 1a. The split is REAL. Verified, not assumed.

`product_mu_phi` (`orpheus/numerics/quadrature/rules_product.py:56-155`) returns
a 2-tuple and the two halves have **disjoint** consumer sets:

```python
# orpheus/numerics/quadrature/rules_product.py:142-155
    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=SPACE_SPHERE,
        invariance_group=SubgroupOfO3.SO2,
        degree_of_exactness=min(2 * n_mu - 1, n_phi - 1),
    )
    structure = LevelStructure(
        n_levels=n_mu,
        level_indices=level_indices,
        level_mu=mu_gl,
    )
    return measure, structure
```

**Graph evidence** (`mcp__nexus__callers` on `product_mu_phi`, transitive): the
only production caller is `Quadrature.product`
(`orpheus/numerics/quadrature/directional.py:596`). Everything else is a test or
`orpheus/numerics/symmetry.py` (`_octahedral_ops`, `_check_invariance_3d` — these
call it to *verify group invariance*, i.e. they consume the MEASURE).
`scratch/b3_4a_mutations.py` is a scratch mutation harness.

**`LevelStructure` consumers** — `mcp__nexus__context` shows only
`Quadrature.level_structure` + the two rule factories + `WindowedSweep.p`; the
authoritative enumeration is the field grep (L-009: Nexus does not model
dataclass fields, so grep is the primary evidence here). Every production reader
of `level_indices` / `n_levels` / `level_mu`:

| Consumer | What it reads it FOR |
|---|---|
| `orpheus/numerics/quadrature/directional.py:489-523` | the three pass-through properties `n_levels` / `level_indices` / `level_mu` |
| `orpheus/geometry/reduced_operator.py:208`, `:576-582`, `:779`, `:791`, `:807` | **α-recursion**, per-level `ΔA/w`, per-level `−sinθ` start, and the `(mu_level_idx, direction_idx) → global n` decode |
| `orpheus/sn/sweep/pole_angular_closure.py:175`, `:262-269`, `:595`, `:893-934`, `:1019`, `:1070`, `:1184` | **τ_raw / τ**, `c_in`/`c_out`, the per-level ψ½ recurrence, per-level scatter/gather |
| `orpheus/sn/mesh/augmented_mesh.py:1261`, `:1317`, `:1468` | `_representative_ordinate` + the `mu_level_idx → global ordinate` decode used by the walk |
| `orpheus/sn/sweep/cache.py:247-252` | the cached per-level ordinate arrays for the hot path |
| `orpheus/sn/loss_representation/__init__.py:2827`, `:2902`, `:3518`, `:4089-4091`, `:4155`, `:4609-4611`, `:4621` | the curvilinear scan's level loop + the `#280 2.5b` seed fold |
| `orpheus/sn/operators/radial_characteristic.py:780`, `:789`, `:852`, `:859` | the ψ½ radial-characteristic march (`A_BB.solve`) |
| `orpheus/transport/radial_characteristic_field.py:289-315` | the per-level Legendre **moment** fold + the most-inward-ordinate corner datum |
| `orpheus/derivations/discrete/sn/contamination.py:135`, `:189` | a diagnostic |

**Nothing outside the cylindrical/spherical curvilinear path reads
`LevelStructure` at all.** That is the verification the brief asked for.

### 1b. The MEASURE is genuinely NOT cylindrical-specific — halving it is refuted

Three independent pieces of evidence that `product_mu_phi`'s `DiscreteMeasure`
must keep its full-sphere contract:

1. **`Quadrature.product` reaches 2-D Cartesian.** `[M]` grep finds the pairing in
   `tests/sn/sweep/core/test_sweep_regression.py:280` (`Mesh2D` + product), plus
   13 more files that use both `Quadrature.product(` and a `Mesh2D`/`_d3_axes`
   construction (`tests/sn/primitives/test_axis_native_construction.py`,
   `tests/sn/operators/test_sn_boundary_operator.py`,
   `tests/sn/regression/_generate_snapshots.py`,
   `orpheus/derivations/continuous/mms/sn.py`, `orpheus/numerics/operator.py`, …).
   A 2-D Cartesian sweep needs all four `(±η, ±ξ)` sign combinations —
   halving the measure deletes half its octants.
2. **The registry types it as a sphere rule with a declared exactness contract**
   (`orpheus/numerics/quadrature/registry.py:471-482`): `invariance_group=SO2`,
   `degree_of_exactness = min(2n_mu−1, n_phi−1)`, and it is one of only four specs
   `select_quadrature` can return for ANY geometry whose group is a subgroup of
   `SO2` — i.e. `slab`, `sphere`, **and** `cylinder` (`GEOMETRY_GROUPS`, `:492-497`).
3. **`orpheus/numerics/symmetry.py:651` `_check_invariance_3d` and `:802`
   `_octahedral_ops` call `product_mu_phi` directly** to test group invariance. A
   ξ≥0 half-set is invariant under **neither** `SO2` nor `C_{n_phi}` — the declared
   `invariance_group` would become a lie.

`[M]` And the moment probe below is the fourth: the half-set with doubled weights
is not a valid `S²` cubature (it integrates ξ-odd polynomials to nonsense).

**Verdict: the brief's split is CORRECT.** The fold belongs in the level
structure + sweep, not in the `DiscreteMeasure`.

---

## 2. The measurements

### 2a. The three candidate level constructions, level 0, `n_mu=4, n_phi=8`

`[M]` (probe: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe2.py`, `probe4.py`)

```
CURRENT full-circle left-endpoint   (n_mu=2 shown for the phi trace)
 phi/pi  : 1.00 1.25 0.75 1.50 0.50 1.75 0.25 0.00      <- the two halves INTERLEAVED
 eta     : -0.8165 -0.5774 -0.5774 -1.5e-16 5.0e-17 0.5774 0.5774 0.8165
 tau_raw : 0  1  0  1  0  1  0  1                        <- the alternation
 alpha   : 0, 0.6413, 1.0947, 1.5482, 1.5482, 1.5482, 1.0947, 0.6413, 2.2e-16
                                    ^^^^^^^^^^^^^^^^ a PLATEAU: the tied pair advances alpha by 0

OPTION A  fold-existing-nodes (omega in [0,pi], endpoints incl., trapezoid w)
 phi/pi  : 1.00 0.75 0.50 0.25 0.00                      <- strictly monotone
 eta     : -0.8165 -0.5774 5.0e-17 0.5774 0.8165
 w       : 0.7854 1.5708 1.5708 1.5708 0.7854   (sum 6.2832 — PRESERVED)
 tau_raw : 0  0.2929  0.5  0.7071  1
 alpha   : 0, 0.6413, 1.5482, 1.5482, 0.6413, -3.3e-16   (dome, closes)

OPTION B  Hebert midpoint half-range (2F(p) nodes strictly inside (0,pi))
 omega/pi: 0.875 0.625 0.375 0.125
 eta     : -0.7543 -0.3125 0.3125 0.7543
 w       : 1.5708 x4                            (sum 6.2832 — PRESERVED)
 tau_raw : 0.2195 0.4142 0.5858 0.7805
 alpha   : 0, 1.1849, 1.6757, 1.1849, -2.2e-16  (dome, closes)
```

Every claim in the hypothesis box is confirmed on both A and B:

* **η injective within a level** — `[M]` `min|Δη| = 0.239` (A) / `0.442` (B) at
  `n_phi=8`, versus `0` (a bit-exact tie) today. The ordering question dissolves:
  ascending-η **is** descending-ω, uniquely.
* **`eta_edge` no longer collapses onto nodes** — `[M]` `τ_raw` leaves `{0,1}` for
  every interior ordinate; the STEP/DIAMOND alternation is gone.
* **ψ even in ξ by construction** — only one member of each mirror pair exists, so
  the asymmetry is *unspellable*, not merely tested-for. (Pattern 4.)
* **~2× fewer angular unknowns** — `[M]` `n_phi/2 + 1` per level (A) or `n_phi/2`
  (B), versus `n_phi`. At `n_phi = 8` that is 5 or 4 versus 8.

### 2b. α: BOTH constructions reproduce the closed form; α does NOT discriminate

`[M]` With the half-angle boundaries taken at **ω-midpoints** (the tree's own
convention, `test_alpha_closed_form.py:317`), `α / (w_gl · ξ(ω_half))` is constant
to `<1e-12` and equals `2κ`:

```
n_phi   8 -> 2.052344    (kappa = d_omega / (2 sin(d_omega/2)) = 1.026172)
n_phi  16 -> 2.012909
n_phi  32 -> 2.003216
```

identical for A and B (they share `Δω = 2π/n_phi`). **α is therefore NOT the
discriminator between the two candidate half-range rules** — a useful negative
result, because the issue's α evidence is sometimes read as selecting Hébert's
exact node placement. It selects the half RANGE, not the node OFFSET.

### 2c. The edge convention: η-midpoints are RIGHT; BMC cumulative-weight edges are not

`[M]` A side question that had to be closed before recommending anything. BMC
(2010) Eq. 52 defines the level's cell edges by the cumulative normalized weight
(`μ_{1/2} = −√(1−ξ²)` marching up by `W_m`). ORPHEUS's **sphere** branch uses
exactly that (`pole_angular_closure.py:576-579`, `mu_edge[n+1] = mu_edge[n] + w[n]`)
but its **cylinder** branch uses η-**midpoints** (`:599-603`). Applying the BMC
convention to the cosine-spaced azimuthal grid gives `τ_raw` **outside `[0,1]`**:

```
n_phi=16, level 0, OPTION A, w-sum edges:
  tau_raw = [0, -0.1955, -0.3284, -0.0307, 0.5, 1.0307, 1.3284, 1.1955, 1]
```

because a cosine grid is *not* weight-equispaced in η. So the cylinder's
η-midpoint convention is the correct one and is **not** a second latent bug —
this closes that thread. (It also means the `[½, 1]` clamp is doing real work on
BOTH half-range candidates; see §4.)

### 2d. The one genuine break — ξ-ODD spherical-harmonic moments

`[M]` (probe: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe5.py`) On
`Quadrature.product(n_mu=4, n_phi=8)`, with a synthetic **ξ-even** ψ and the fold
`{ξ ≥ 0}` + doubled interior weights (weight sum preserved to 3e-15):

```
 (l,m) parity under xi -> -xi:  l=0 s0 EVEN | l=1 s0 EVEN | l=1 s1 EVEN | l=1 s2 ODD

 moment  l m |    full-sphere       folded(2w)        rel.diff
        0 0 |  6.512041e+00   6.512041e+00    1.4e-16
        1 0 |  2.345248e-01   2.345248e-01    3.6e-16
        1 1 |  7.238071e-01   7.238071e-01    4.6e-16
        1 2 | -1.313633e-16   2.943243e+00    2.2e+16   <-- the ONLY break
```

Every ξ-even slot is exact; the ξ-odd slot is destroyed. **This is where the
hypothesis needs an amendment, not a retraction:** on the folded set the ξ-odd
harmonics have no representation. `Σ_{half} 2w Y_odd ψ` is not "an inaccurate
moment", it is the pairing of ψ with a functional that is identically zero on the
space the fold lives in. The fold must carry that fact structurally.

**Who consumes a ξ-odd cylindrical moment TODAY (so this is not hypothetical):**

* `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:240-256`
  assembles `q_n^{P1} = (1/W)·3·Σ_s1·Σ_{m=0,1,2} Y_1^m(Ω_n) φ_1^m` in CYLINDRICAL
  geometry — the `m_slot = 2` term is the ξ-odd one. The production path it pins
  is `S` at `scattering_order=1` from `build_within_group_system`.
* `orpheus/transport/radial_characteristic_field.py:296-304` takes **Legendre
  moments over each level's η-nodes** — `legvander(mu_p, M_p−1)` with `M_p` the
  level size. Folding a level changes `M_p`, hence the moment ORDER this fold
  reconstructs from. (Not a parity break — an order/rank change. Flagged in §3.)

---

## 3. The consumer blast radius

The single most important structural distinction — **it decides most of §2d**:

> **(a1) fold the SWEEP only** — `LevelStructure` carries the half level, the
> `_OneDimScanWalk` marches `M_p/2` ordinates, and the mirror partners are
> **filled by symmetry** into the full `(N, ng, nx)` buffer before the walk
> returns. The STATE stays `N`-dimensional.
>
> **(a2) fold the STATE** — the angular DOF set itself becomes the half set.
> The "≈2× fewer angular unknowns" bullet of the hypothesis is realized only
> here.

`[M]` Under **(a1)** the §2d moment break **does not occur at all**: the moment
integral is taken on the reconstructed full-`N` ψ, which is now *exactly* ξ-even,
so `φ_1^{ξ}` evaluates to its structural zero by construction rather than by
cancellation. Under **(a2)** the break is live and must be paid for.

### 3a. Consumers that change under (a1) — the sweep-only fold

| Site | What changes |
|---|---|
| `orpheus/numerics/quadrature/rules_product.py:137-140` | the `argsort` becomes a *selection* + monotone order (`ξ ≥ 0`); the sort key becomes injective |
| `orpheus/numerics/quadrature/rules_sphere.py:210-215` | same for `level_symmetric` (a 4-to-1 fold — see §5) |
| `orpheus/numerics/quadrature/rules_sphere.py:111-136` `LevelStructure` | needs the **folded weight** (`w`/`2w`) and the **mirror partner map** as new fields — `level_indices` alone can no longer reconstruct the level's quadrature |
| `orpheus/geometry/reduced_operator.py:779-786` (α), `:791-795` (`ΔA/w`), `:807-809` (`−sinθ`) | unchanged in FORM; the arrays shrink to `M_p/2`. `[M]` The redistribution term `(ΔA/w)(α_out ψ_out − α_in ψ_in)` is **invariant** under a uniform level-weight rescale (α scales with `w`, `ΔA/w` scales as `1/w`), so the doubled weights are self-consistent here |
| `orpheus/sn/sweep/pole_angular_closure.py:588-611` (`τ_raw`) | array shrinks; `[M]` values leave `{0,1}` for every interior ordinate |
| `pole_angular_closure.py:664-675` (the `[½,1]` clamp) | **still active** — `[M]` at `n_phi=16` option A gives `τ_raw = [0, .26, .40, .46, .50, .54, .60, .74, 1]`, so the whole `η < 0` half is still clamped. The clamp's stated rationale ("`τ_raw = 0` bit-exactly ⟹ divide-by-zero") only covers the ONE endpoint. **Re-examining the clamp is in scope of this fix, not out of it.** |
| `pole_angular_closure.py:1012-1029` `_edge_seed_stencil` | the `abs(mu[cand] − mu[m0]) > 1e-14` tie-skip loop becomes dead code (no ties) — `order[1]` always |
| `orpheus/sn/loss_representation/__init__.py:4089-4091`, `:4609-4611` | the level loop shrinks; **plus the NEW mirror-fill** after the walk |
| `orpheus/sn/loss_representation/__init__.py:2827-2845` (`_dag_legs`), `:2902-2915` (`_degenerate_positions`) | the degenerate `\|η\|<1e-15` set halves (one `φ=π/2` node per level instead of two), and vanishes entirely when `n_phi/2` is odd |
| `orpheus/sn/mesh/augmented_mesh.py:1297-1330` `_representative_ordinate`, `:1468` | **no change needed** — all addressing is `(mu_level, local_idx)` and parametric in `M_p` |
| `orpheus/sn/sweep/cache.py:245-255` | `level_visits_iter` covers only the kept half ⟹ `chain_idx[global_n]` is populated for half the ordinates. **This is where the mirror-fill obligation surfaces concretely.** |
| `orpheus/transport/radial_characteristic_field.py:296-304` | `legvander(mu_p, ords.size − 1)` — the Legendre fold order is `M_p − 1`; halving `M_p` **halves the moment order** this reconstructs from. Not a parity bug, a *rank* change. Sphere-only today (§5), but it must be re-derived if a folded cylinder ever carries |
| `orpheus/transport/radial_characteristic_field.py:313-315` | `ords[argmin(mu[ords])]` — collapses to `ords[0]` by the sort contract, and the tie that made it a second free tie-break is **gone** (issue AC #3 satisfied for free) |
| `orpheus/sn/loss_representation/__init__.py:3294`, `:3470`, `:4189`, `:4593` (`reflection_index("x")`) | `[M]` **unchanged and now unambiguous.** The x-reflection closes the `ξ ≥ 0` half for every `n_phi ∈ {4,8,16,32}`. The flagged `rx` vs axis-crossing `(−η,−ξ)` discrepancy — `[M]` 56/64 ordinates differ at `product(4,16)`, the 8 agreeing ones being exactly the `ξ=0` self-mirror pair per level — **dissolves**, because on the folded set `(−η,−ξ) ≡ (−η,+ξ)` by the ξ-evenness that now holds by construction |
| **NEW** — the mirror-fill / lift | a `ξ`-mirror scatter `ψ[ry[n]] ← ψ[n]`. Architecturally this is the **reconstruction face** of a fold/lift pair on the ξ-even sector; ORPHEUS already has the vocabulary (`FrameBase`: analysis = fold, reconstruction = lift) |
| **NEW** — a well-posedness GUARD | the fold is valid **iff the source and the boundary data are ξ-even**. Isotropic sources ✓; `P_ℓ` scattering sources ✓ (they are built from moments of an even ψ); `[D]` both cylindrical MMS ansatzes are functions of `η` and `ξ²` (proved by `tests/sn/verification/mms/test_mms_ordering_blindness.py`). A user-supplied per-ordinate `AngularSourceSink` or a prescribed-inflow trace CAN be ξ-odd — that case must be **refused loudly**, not silently symmetrized |

### 3b. Consumers that change ADDITIONALLY under (a2) — the state fold

Everything in 3a, plus:

| Site | What changes |
|---|---|
| `orpheus/numerics/quadrature/directional.py:401-458` `spherical_harmonics` / `angular_frame` | the analysis face must restrict to the **ξ-even subspace** of the SH basis. Symmetrized basis ⟹ still `GalerkinFrame`; plain-Y reconstruction + folded analysis ⟹ `PetrovGalerkinFrame`. `[M]` Without this, `φ_1^{ξ}` goes `−1.3e-16 → +2.94` |
| `orpheus/transport/operators/scattering.py` (the `S` in the `(L+C), S, B` triple) | consumes `angular_frame`; inherits the restriction |
| `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:240-256` | its hand-reference sums **all three** `Y_1^m` slots; the ξ-odd one must be shown zero |
| `orpheus/numerics/quadrature/directional.py:464-487` `octants` | `[M]` a ξ≥0 set populates 4 of 8 sign-octants. Consumers: `loss_representation/assembly.py:137`, `sweep_schedule.py:105`, `:125` |
| `orpheus/numerics/spaces/angular_trace_space.py:338` `partial_current_metric` (`\|Ω·n\|·w`) | shrinks; `net_current` (`orpheus/sn/solver.py:1410-1423` leakage) stays correct with doubled weights because `Ω·n = η` is ξ-even |
| `FullFieldSpace` / `to_flat` / `n_dof` / Krylov `restart` | the composite resizes — the **ERR-053 family** hazard the user memory flags for exactly this shape of carve |
| every cylindrical regression snapshot | the ordinate axis length changes |

### 3c. What does NOT change under either variant

* The `DiscreteMeasure` — nodes, weights, `invariance_group=SO2`,
  `degree_of_exactness`. (§1b.)
* `orpheus/numerics/symmetry.py:651,802` group-invariance checks.
* Every 2-D / 3-D Cartesian consumer of `Quadrature.product`.
* `Quadrature.reflection_index` — the table is built from the full node set and
  the x-reflection restricts to the folded half `[M]`.
* The slab and sphere paths (`level_structure is None` / single level).

---

## 4. The three candidate homes

### (a) Fold inside `LevelStructure` — **RECOMMENDED**

*The measure is untouched; the level structure carries the half range; the sweep
marches it; the mirror partners are filled by symmetry.*

**Why it is right in the operator-algebra sense.** 1-D cylindrical geometry has
a `Z₂` deck transformation (the `r–z` mirror) acting on `S²`. The physical phase
space is the **quotient** `S²/Z₂`; the full `S²` is its double cover. The
`DiscreteMeasure` is honestly a cubature on the COVER (and must stay so — it is
consumed by geometries whose group is larger). `LevelStructure` is honestly a
description of the QUOTIENT — it exists only for the one geometry that has this
symmetry. **Putting the fold in `LevelStructure` puts the quotient where the
quotient already lives.** No twin path: there remains exactly one product rule,
one level structure, one α producer, one τ producer.

* **Breaks:** the sweep no longer covers all `N` ordinates (§3a `cache.py` row) —
  a mirror-fill/lift must be introduced. `LevelStructure` gains two fields
  (folded weights, mirror map) and stops being a pure *index* structure.
* **Cost:** contained. `[M]` `τ_raw[0] = 0` is preserved under construction (A),
  so the **R12a trichotomy, the `#280 2.5b` direct-seed fold, and the whole
  `radial_characteristic_*` sphere-only story survive untouched.**
* **Twin path:** none.
* **Sub-decision (a1 vs a2)** — fold the sweep only, or the state too. Recommend
  **(a1) first**: it is the correctness fix, it makes ψ-even-in-ξ unspellable, and
  it costs no moment/frame/Krylov surgery. (a2) is a *performance* follow-up whose
  price is the ξ-even trial-space restriction and an `n_dof` resize.
* **Node placement:** use construction **(A)** (fold the existing nodes, endpoints
  included, trapezoid weights). `[M]` It reproduces the α closed form exactly, it
  is the construction the tree's own `foundation` test already pins
  (`test_alpha_closed_form.py:255-268`), and it keeps `τ_raw[0] = 0`.

### (b) A distinct cylindrical quadrature family alongside the sphere rule

*A `cylindrical_half_range(n_mu, n_phi)` rule producing Hébert's `2ℱ(p)`
midpoint nodes on `(0, π)`.*

* **Breaks:** `[M]` **every level becomes R12a-CARRYING** (`τ_raw[0] = 0.2195` at
  `n_phi=8`). That routes every cylindrical solve through
  `RadialCharacteristicOperator.solve` / the ψ½ block — machinery that is
  sphere-only in practice today — and **retires the `#280 2.5b` direct-seed fold
  for the cylinder**. It also changes the sphere measure's nodes (the mirrored
  half-range is the *midpoint* circle rule, offset by `Δφ/2`), so every snapshot
  moves, unless the new family is kept entirely separate from `product`.
* **Cost:** high, and it is *structural* cost, not just churn.
* **Twin path:** **YES, and it is the disqualifying objection.** Against
  `coding-elegance` Pattern 6 (build primitives, not products): a second family
  would be *the same GL×equispaced product rule* with a different azimuthal
  offset and a different level partition. It is not a new primitive — it is the
  same primitive with a different **quotient**. Two rules that must stay
  numerically consistent (the cylinder answer must not depend on which one you
  built) is precisely the twin-source-of-truth smell.
* **When (b) WOULD be right:** if the cylinder genuinely wanted an *open* rule on
  the quotient interval (the sphere's Gauss–Legendre situation). `[M]` It does
  not — the quotient interval's endpoints `ω = 0, π` are the `ξ = 0` directions
  where the redistributed flux `ξψ` **vanishes identically**, so standing on them
  costs nothing and buys `α_{1/2} = α_{M+1/2} = 0` exactly. Closed (endpoint-
  inclusive) is the *natural* rule here, and that is construction (A).

### (c) Change `product_mu_phi` itself to half-range

* **Breaks:** refuted in §1b — four independent ways (2-D Cartesian consumers,
  the registry's `SO2`/degree contract, the two `symmetry.py` invariance checks,
  and the `[M]` ξ-odd cubature failure). The returned object would no longer be a
  sphere cubature while still being typed and registered as one.
* **Cost:** unbounded — it silently changes a shared primitive for consumers that
  have nothing to do with cylindrical geometry.
* **Twin path:** it *creates* one — either a second full-sphere rule has to be
  minted for the non-cylindrical consumers, or `SubgroupOfO3.SO2` becomes a lie.
* **Verdict: refuted.**

### Recommendation, in one line

**(a1) with node construction (A)**: fold in `LevelStructure`, march the half
level, lift by the ξ-mirror, guard ξ-even data, keep the measure a full-sphere
cubature. Then, as a separate follow-up, consider (a2) — and only then does the
ξ-even trial-space restriction have to be built.

---

## 5. The `level_symmetric` question

### 5a. Does it reach the cylindrical sweep in production?

**Two answers, and the distinction matters.**

* **The SELECTOR refuses it.** `[M]`
  `select_quadrature(geometry="cylinder", target_degree=…)` rejects
  `LevelSymmetricSN` with, verbatim from the returned `SelectionLog`:

  ```
  rejected=[('LebedevSphere',    "G mismatch: geometry SO2 is not a subgroup of rule's invariance group Oh"),
            ('LevelSymmetricSN', "G mismatch: geometry SO2 is not a subgroup of rule's invariance group Oh")]
  ```

  (`GEOMETRY_GROUPS["cylinder"] = SO2`; `registry.py:492-497`.) So no
  registry-driven path ever hands LS to a cylinder.
* **But the selector has NO production caller.** `[M]` grep: every
  `select_quadrature` call site outside the module itself is in
  `tests/numerics/test_registry.py`. The registry docstring calls it *"opt-in
  convenience"* (`registry.py:68`). Production constructs
  `SNMesh(mesh, quad, materials)` directly, and
  `cylindrical_streaming` (`orpheus/geometry/reduced_operator.py:763-770`)
  accepts **any** quadrature carrying a `level_structure` — LS included.
* **And LS-on-cylinder is an explicitly SUPPORTED, TESTED configuration.** The
  `#280 2.5b` fold comment (`loss_representation/__init__.py:4142-4150`) reasons
  about it by name (*"Level-symmetric rules have c_in = 0 (dead seed) so the fold
  is an exact no-op there"*), the R12a doc table gives it its own row, and it is
  pinned in `tests/sn/sweep/test_cyl_direct_seed_fold.py`,
  `tests/sn/sweep/curvilinear/test_cyl_sweep_regression.py`,
  `tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py`, and
  `tests/sn/mesh/test_radial_characteristic_*.py`.

**Verdict:** `product` is the *intended* cylindrical rule; LS is reachable,
supported and pinned. A fix that handles only `product` leaves a live worse-case
path (`[M]` LS measures `3.08e-1` on the mirror criterion vs `1.19e-1` for
`product(4,8)`).

### 5b. Does the same fold apply? Yes — but it is 4-to-1, not 2-to-1

`[M]` An LS "level" is keyed on `|μ_z|`
(`rules_sphere.py:213`, `np.abs(np.abs(mu_z_arr) - level_mu_vals[p]) < tol`), so
it contains **both** signs of `μ_z` **and** both signs of `ξ`:

```
LS4 level 0: M=16, 4 distinct eta          LS6 level 0: M=24, 6 distinct eta
LS4 level 1: M= 8, 2 distinct eta          LS6 level 1: M=16, 4 distinct eta
                                            LS6 level 2: M= 8, 2 distinct eta
```

Both redundancies are **physical** for a 1-D infinite cylinder: the `r–z` mirror
gives `ψ` even in `ξ`, and the `z → −z` mirror gives `ψ` even in `μ_z`. Neither
`μ_z` nor `ξ` enters the 1-D cylindrical sweep except through the level LABEL
(`sin_theta` from `mu_z[level_idx[0]]`, which is sign-blind) and the per-ordinate
source. So the 4-to-1 fold `{ξ ≥ 0, μ_z ≥ 0}` with weight ×4 is legitimate.

`[M]` Measured consequence (probe `probe7.py`), with the fold applied and the
level weight sum preserved to 1e-15:

```
        tau_raw NOW (LS4 level 0)   : [1, .5, .5, 0]  x4   (the 4-fold fingerprint)
        tau_raw FOLDED (LS4 level 0): [0.320715, 0.333333, 0.666667, 0.679285]
        tau_raw FOLDED (LS4 level 1): [0.292893, 0.707107]
        tau_raw FOLDED (LS6 level 0): [0.274209, 0.377875, 0.402130, 0.597870, 0.622125, 0.725791]
        R12a `0 < tau_raw[0] < 1` after the fold:  True for EVERY LS level
```

⚠ **The `[1, ½, ½, 0]` repeating pattern the R12a doc table calls a property of
*"level-symmetric rules"* is the 4-fold-covering fingerprint, exactly as
`[0,1,0,1,…]` is the 2-fold one for `product`.** Under the fold the LS row of the
trichotomy **ceases to exist**: LS levels become CARRYING.

**So LS needs its OWN treatment, and it is the harder one.** Unlike `product`
construction (A), there is no LS analogue of "the endpoint `ω = π` is already a
node": `[M]` LS's most-negative folded `η` is strictly inside `(−sinθ, …)`, so
`τ_raw[0] ∈ (0,1)` unavoidably. Folding LS therefore **does** activate the ψ½
`radial_characteristic_*` route for the cylinder — the same consequence that
disqualifies candidate (b) for `product`. Options, in preference order:

1. **Fold `product` (A) now; leave LS unfolded and mark its cylindrical use as
   the known-defective path** — with the `xfail(strict=True)` row at
   `test_azimuthal_mirror_symmetry.py::test_level_symmetric_cylinder_is_even_in_xi`
   kept RED as the standing admission. Cheapest, honest, no half-fix.
2. **Refuse LS on a cylinder at construction** — make `cylindrical_streaming`
   reject an `O_h` quadrature, aligning the direct-construction path with what
   `select_quadrature` already decides. This is a *retirement*, and it has a real
   test-migration bill (§5a lists six modules).
3. **Fold LS 4-to-1 and accept carrying LS levels** — correct, and it exercises
   the ψ½ machinery on a second geometry (arguably an architectural win, since
   that machinery currently has exactly one live instance). Largest blast radius.

This is a genuine user call, not an explorer call (lessons L-004): (1) and (2)
turn on whether LS-on-cylinder is a capability ORPHEUS wants to keep.

---

## 6. The theory corpus — claims vs code (⚠ L33)

### D1 (HEADLINE). `§sn-direct-seed-circle-vs-interval` states an architectural principle whose load-bearing half is FALSIFIED

`docs/theory/methods/sn/curvilinear_one_group.rst:3867-3963`. `[D]` verbatim:

> **Cylinder — the redistribution axis is a circle.** At a fixed polar cosine
> `μ_z = cos θ`, the cylinder redistributes across the **azimuthal angle** `φ`,
> which lives on a **circle** `[0, 2π)` — a *periodic* domain. […] Equispaced
> sampling of a smooth periodic function is the trapezoidal rule, which on a
> circle is **spectrally accurate** […] — there is no accuracy penalty for the
> choice. And crucially, for **even** `n_φ` the grid hits `φ = π` exactly […] the
> starting-edge ordinate `η_0 = −sin θ` therefore lands *exactly on a quadrature
> node* — `τ_raw,0 = 0` — and the seed is a **bulk ordinate for free, at no
> accuracy cost**.

and the closing principle:

> **a periodic redistribution axis gives edge-inclusion for free; an interval
> axis makes you pay for it with a separate seed.**

**Adjudication.** The *α-recursion's* axis is not the circle. `α_{m+1/2} =
α_{m−1/2} − w_m η_m` with `α_{1/2} = 0` is a **cumulative integral along a
contiguous, monotonically-advancing** azimuthal path — a march on an INTERVAL,
with a Dirichlet closure at both ends. `[M]` The as-built order visits

```
φ/π :  1.00  1.25  0.75  1.50  0.50  1.75  0.25  0.00
```

— the two half-circles **interleaved**, monotone on neither branch nor the whole
circle. `[M]` The α "dome" carries a **plateau** (`1.5482, 1.5482, 1.5482` at
`n_mu=2, n_phi=8`): a tied mirror pair advances α by exactly zero, because its
two members' angular cells have zero width. Periodicity is not what gives the
cylinder its `τ_raw,0 = 0`; what gives it that is that **`φ = π` is an endpoint
of the QUOTIENT interval `S¹/Z₂ = [0, π]`, and the endpoint-inclusive grid stands
on it.** The doc's own escape hatch gives the game away — *"It is contingent on
**even** `n_φ` — an odd azimuthal count would miss `φ = π`"* — a genuinely
periodic axis has no distinguished point to miss.

The section's CONCLUSION ("the cylinder does not pay for a seed") survives under
construction (A) and is **falsified** under (B). Its REASON is wrong either way.
The correct statement: *the cylinder's redistribution axis is the interval
`ω ∈ [0, π]` — the quotient of the direction circle by the `ξ`-mirror — and its
endpoints are the `ξ = 0` directions where the redistributed flux `ξψ` vanishes,
so an endpoint-inclusive rule closes the α dome for free.*

The `§sn-direct-seed-lobatto-study` that follows is framed entirely by this
dichotomy and inherits the correction.

### D2. The R12a trichotomy table presents two covering-artefacts as rule properties

`curvilinear_one_group.rst:3826-3856`. `[D]` The `τ_raw = 0` row is attributed to
*"cylinder **product** rules"* and the `τ_raw = 1` row to *"cylinder
**level-symmetric** rules — duplicate-η nodes collapse the midpoint edge onto
η_0"*. The second row **names the degeneracy explicitly** and then treats it as a
property of the RULE rather than of the covering. `[M]` Under the fold, both rows
change (§2a, §5b) and the trichotomy reduces to a dichotomy. The predicate itself
(`τ_raw,0 ∈ (0,1)`) is sound and unaffected — it is the TABLE's causal story that
is contingent.

### D3. The ordering convention is stated in the corpus, and its citation is unsupported in BOTH spellings

* `[D]` `curvilinear_one_group.rst:229-232` — *"On level `p`, the ordinates are
  **sorted by increasing η**"*. **TRUE of the code** (`rules_product.py:139`).
* `[D]` `angular_quadrature.rst:46` and `:65-68` — same claim, twice, with
  `angular_quadrature.rst:66-68` attributing it to `:cite:`BaileyMorelChang2010``.
* `[D]` `rules_product.py:40-44` + `:66-68` attributes the same claim to
  **Bailey, Adams, Yang & Zika (2009) Eq. 50** — the paper
  `curvilinear_one_group.rst:1830-1850` explicitly **retracts** as *"the **wrong
  Bailey paper** — a piecewise-linear FE diffusion paper unrelated to S_N"*.
  That retraction's Phase-B fix list names `reduced_operator`,
  `loss_representation`, `transport.spatial.diamond`, `pole_angular_closure` —
  **`rules_product.py` is not on it** (the module postdates the fix). So a
  retracted citation is live at HEAD in the very function this issue is about.
* And BMC (2010) **Eq. 50** is the two-sided closure `α_{1/2,n} = α_{M+1/2,n} =
  0`; it states no ordering convention. **Neither spelling of the citation
  supports the claim.**

### D4. NO page states the level covers the azimuthal range TWICE — and no page states the mirror theorem

* `[M]` grep of `docs/theory/methods/sn/*.rst` for *"even in"*, *"mirror
  symmetry"*, *"xi-mirror"*, *"azimuthal mirror"* returns **zero hits**. The
  `ψ(η, ξ) = ψ(η, −ξ)` theorem — the criterion the whole issue turns on — is
  documented **only in a test docstring**
  (`tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py:1-80`).
* `[M]` grep for `326` in `docs/theory/methods/sn/*.rst` returns zero hits. The
  defect has no corpus presence at all.
* The nearest thing to a range statement is D1's *"lives on a circle `[0, 2π)`"*,
  which asserts the full circle is **correct**. So the corpus does not merely
  omit the double covering — it **endorses** it.

### D5. Two claims that CHECK OUT

* `[D]` `curvilinear_one_group.rst:240-241` *"Each level's α values form an
  independent dome from `η = −sinθ` to `η = +sinθ`"* and `:245`
  *"`α_{n+1/2} ≥ 0`"* — `[M]` both TRUE today and under both fold constructions
  (the dome closes to `≤ 3.3e-16`). This is the telescoping-invariance of
  lessons L-014: these hold under **every** permutation, which is exactly why
  they never caught the defect.
* `[D]` `angular_quadrature.rst:41-46` — *"the Level-Symmetric quadrature groups
  both `+μ_z` and `−μ_z` hemispheres on the same level (grouped by `|μ_z|`)"* —
  `[M]` TRUE (`rules_sphere.py:213`).

---

## 7. Where the hypothesis breaks — the refutation attempts, scored

| Attempted refutation | Outcome |
|---|---|
| "Halving the measure breaks the sphere-cubature contract" | **SUSTAINED** — but it refutes candidate **(c)**, not the hypothesis. §1b: four independent grounds. |
| "The α recursion needs the full circle to close its dome" | **REFUTED** `[M]` — the dome closes on both half-range constructions to ≤3.3e-16, and only there does α equal its closed form. |
| "The redistribution term is sensitive to the doubled weights" | **REFUTED** `[M]` — `(ΔA/w)(α_out ψ_out − α_in ψ_in)` is invariant under a uniform level-weight rescale. |
| "The `[½,1]` clamp exists because of the tie and will misbehave" | **PARTIAL** `[M]` — the clamp survives and is still active on the whole `η<0` half of a folded level. It is not broken by the fold, but its justification (a single `τ_raw = 0`) no longer matches its reach. **Re-derive it in the same change.** |
| "`reflection_index` / specular BCs break" | **REFUTED** `[M]` — the x-reflection closes the `ξ ≥ 0` half for every `n_phi`; `reflection_index("y")` has no production consumer. |
| "The pole map at `loss_representation:4189` breaks" | **REFUTED, and it is FIXED by the fold** `[M]` — the `(−η,+ξ)` vs `(−η,−ξ)` ambiguity (56/64 ordinates differ at `product(4,16)`) is exactly the ξ-mirror, which becomes the identity on the quotient. |
| "Currents / leakage break" | **REFUTED** — `net_current` reads the `\|Ω·n\|·w = \|η\|·w` metric, ξ-even. |
| "Moment / SH integration breaks" | **SUSTAINED for (a2), VOID for (a1)** `[M]` — the ξ-odd slot `φ_1^{ξ}` goes `−1.3e-16 → +2.94`. Void under (a1) because moments are taken on the lifted full-`N` ψ. |
| "A ξ-ODD source makes the fold invalid" | **SUSTAINED — and this is the real precondition.** Nothing in the tree forbids a user handing a ξ-odd `AngularSourceSink` or prescribed inflow to a cylinder. The folded path must **refuse** it. |
| "`level_symmetric` needs different treatment" | **SUSTAINED** — §5b: 4-to-1, and folding it flips every LS level to R12a-carrying. |

**Net:** the hypothesis is not refuted. Its "2× fewer unknowns" bullet is the one
that carries a real price (`[M]` and the factor is `1.6×` at the production
`n_phi = 8`, `1.88×` at 32 — asymptotically 2, never exactly 2 under
construction (A)). Everything else is a clean win.

---

## 8. Gaps — what I could NOT settle

1. **Whether the `[½,1]` clamp should survive the fold.** `[M]` shows it stays
   active on the whole `η < 0` half. Deciding it needs the flat-flux /
   asymptotic-diffusion-limit argument from BMC §Eq. 43, which is a
   `numerics-investigator` / `literature-researcher` question, not an
   exploration one.
2. **Whether the lift belongs on `LevelStructure`, on the walk, or as a typed
   `Frame`.** All three are spellable; the choice interacts with the (a1)/(a2)
   decision and with `#325`. I recommend the typed-fold/lift pair on principle
   (it makes the quotient a first-class object rather than an implicit index
   trick), but I did not measure a cost for it.
3. **`radial_characteristic_field.py:296-304`'s Legendre fold order.** It uses
   `legvander(mu_p, M_p − 1)`. Halving `M_p` halves the reconstructed order. This
   path is sphere-only today (`n_levels = 1`, `M_p = N`), so the fold does not
   reach it under (a1)+(A) — but it DOES under (b) and under a folded LS. Not
   analysed further.
4. **Whether an odd `n_phi/2` (e.g. `n_phi = 6, 10`) has a hidden case.** `[M]`
   the degenerate `η = 0` node disappears then. `n_phi` odd is already rejected at
   `SNMesh` construction (ERR-042) but `n_phi ≡ 2 (mod 4)` is buildable and I did
   not sweep it.

---

## Re-run of the L-007 / L-012 close check

`git status --short` + `git diff --name-only` re-run at close. **The tree moved**:
`tests/sn/_test_helpers.py` went clean → `+118` and three `#326` test modules
appeared untracked (header). **No production file was modified** at any point, and
none was modified by this investigation (read-only, as briefed).

Every `orpheus/` reference in this document was re-verified by grep at final read:

```
rules_product.py:139         order = np.argsort(mu_x[level_arr])
rules_sphere.py:214          order = np.argsort(mu_x[idx])
reduced_operator.py:785      alpha[m + 1] = alpha[m] - w[m] * eta[m]
pole_angular_closure.py:602  eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])
augmented_mesh.py:815        if eps < float(tau_level[0]) < 1.0 - eps      (R12a)
```

Probes are at `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/probe{1,2,4,5,6,7,8}.py`.
