# A3 `solve_transpose` verification spec — the reverse-DAG inner adjoint solve (#276)

Verification gate chain for **A3** of the adjoint SN transport campaign
(`feature/sn-adjoint-transport`). A3 adds `solve_transpose` — the reverse-DAG
forward-substitution inner adjoint solve on the SN loss `(L+C)` (the analytic
inverse of `(L+C)†`, NOT Krylov-on-transpose, NOT μ-reversal) — and wires
`_AdjointOperator` to advertise `CAP_SOLVE` so `(L+C).H.solve(b)` routes to
`(L+C).solve_transpose(b)`.

Author: test-architect. Design-time only; no production code. Pre-reads cited
inline by `file:line` (trust the code — the campaign-plan STATUS section is
stale, see §0.4).

Canonical gate: `.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE
--timeout=300 -p no:xdist -p no:cacheprovider`.

---

## 0. Structural facts (verified) + what A3 actually changes

### 0.1 The forward chain A3 twins (the single source of truth)

| Piece | Where | What it is |
|---|---|---|
| `(L+C)` = `InvertibleOperator` | `orpheus/sn/operators/streaming.py:509` | the sweep-invertible composite returned by `L + C`; caps `{CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE}` (`:665`) |
| forward matvec `apply` | `streaming.py:712` → `loss_representation.loss_action` | `(L+C)ψ` |
| forward solve `solve` | `streaming.py:761` → `_solve_timed_full_field` (`:862`) → `loss_representation.sweep` (`:1001`) | `(L+C)⁻¹ rhs` — **one** forward-substitution pass (`sweep` at `loss_representation/__init__.py:2160`) |
| reverse matvec `apply_transpose` | `streaming.py:740` → `loss_representation.loss_action_transpose` (`__init__.py:2602`) | `(L+C)ᵀφ` — the reverse-walk (`for i in reversed(cells)` at `:2725`; the curvilinear seed adjoint via the **mirror** ordinate permutation `quad.reflection_index("x")` at `:2658`/`:2758`) |
| **reverse solve `solve_transpose`** | **DOES NOT EXIST** (verified: `grep -rn "solve_transpose\|sweep_transpose" orpheus/` → empty) | **the A3 gap** — `(L+C)⁻ᵀ rhs`, the reverse-DAG forward-substitution; the SOLVE twin of `loss_action_transpose` |

The transpose `(L+C)ᵀ` is upper-triangular (the forward `(L+C)` is lower-tri in
cell-visit DAG order), so forward-substitution on `(L+C)ᵀ` = back-substitution
walking the DAG in **reverse** sweep order — exactly the walk `loss_action_transpose`
already does for the matvec. A3 adds the SOLVE counterpart on the SAME
`LossRepresentation` instance (L21 "matvec ≡ sweep ≡ adjoint-sweep — actions of ONE
operator").

### 0.2 The `.H` wrapper A3 wires (`_AdjointOperator`, `operator.py:615-702`)

- `_AdjointOperator.apply(y) = G_V⁻¹ · inner.apply_transpose(G_W · y)` — the
  Hilbert G-adjoint `A* = G⁻¹AᵀG` (`:668-695`; metric delegated to
  `FunctionSpace.apply_metric` / `apply_inverse_metric`).
- Capability set: advertises `CAP_APPLY` iff inner has `CAP_APPLY_TRANSPOSE`
  (`:646-647`). **It does NOT advertise `CAP_SOLVE`** — the deferral comment is at
  `operator.py:648-651` ("Solve generally does NOT propagate … would need
  `A.H.solve = (A.solve).H`"). The brief cited "~650-651"; the live comment spans
  `648-651` (drift, as anticipated). **There is no `solve` method on
  `_AdjointOperator`** (verified: methods are `__init__`, `domain`, `codomain`,
  `apply`, `apply_transpose` only).

A3 closes this deferral: add `_AdjointOperator.solve` that mirrors `apply` with the
metric roles swapped, and advertise `CAP_SOLVE` **iff inner exposes the reverse
solve**. The math (derive once, gate it):

```
A.H.solve(b)  solves  A* x = b,  A* = G_V⁻¹ Aᵀ G_W
            ⟹ x = (A*)⁻¹ b = G_W⁻¹ A⁻ᵀ G_V b
            ⟹ A.H.solve(b) = inner_codomain.apply_inverse_metric(
                                  inner.solve_transpose(
                                      inner_domain.apply_metric(b)))
```
For the SN `(L+C)` (domain = codomain = `FullFieldSpace`, so `G_V = G_W = G`) this is
`(L+C).H.solve(b) = G⁻¹ · solve_transpose(G · b)` — the exact mirror of
`apply`'s `G⁻¹ · apply_transpose(G · b)` with metric/inverse-metric swapped.

### 0.3 The solve is EXACT (direct), NOT Krylov — the tolerance determinant

The brief asked "iterative-or-exact — determine which." **Determined: direct.**

- **Slab (Cartesian):** `(L+C)` is genuinely lower-triangular; `solve` is **one**
  forward-substitution pass (`sweep` → `_run`, no outer loop). So `solve_transpose`
  is one reverse-substitution pass — **EXACT** up to FP reduction-depth. Tolerance
  is FP-only (`rtol≈1e-11`), **NOT `inner_tol`-scaled.**
- **Curvilinear (sphere/cyl):** the Morel–Montry **Carlson coupled-pole seed reads
  the per-ordinate iterate** (`loss_representation/__init__.py:1005` "Carlson seed
  reads the per-ordinate iterate; lesson L21"), so a cold-start (`initial_guess=None`)
  single call is approximate; the seed is converged by **iterate-threading** (the
  forward analogue loops `LC.solve(rhs, initial_guess=psi)` —
  `test_invertible_operator.py:957-968`, break floor `1e-14`). This is **direct
  machinery iterate-threaded for the seed**, NOT a Krylov `inner_tol`. Converged
  tolerance ≈ `1e-12`; gate at `rtol≈1e-9` (conservative seed-convergence floor).

So: NO `inner_tol` anywhere in A3's intrinsic gate. The reference values are exact
(`np.linalg.solve` on the dense matrix); the only slack is FP / seed-convergence.

### 0.4 Discrepancy flags (trust the code, per `.claude/rules/process-discipline.md`)

1. **`adjoint_sn_transport.md` STATUS (lines 14-62, dated 2026-06-28) is STALE.**
   It frames A2/S† as "forward-migration reverted, Option-2 pending, `full_scatter_kernel`
   UNCOMMITTED." Git refutes: **S† LANDED at `15185e5`** ("SN adjoint S† via
   full_scatter_kernel — closes #118 (#276 P3)") and an entire K_iso refactor
   (`a8bb027`…`9cf70d9`, P1–P5) merged AFTER. Live now:
   `ScatteringOperator.capabilities` carries `CAP_APPLY_TRANSPOSE`
   (`scattering.py:298,454,526`); `S.apply_transpose = (1/W)·full_scatter_kernelᵀ`
   (`scattering.py:1673`); gated by `tests/sn/operators/test_scattering_adjoint.py`.
   `FissionOperator` likewise (F† at `bbd4479`). **⟹ Deliverable 2 (full-loss
   `L+C−S−B` reciprocity) IS reachable now.**
2. **`p6_adjoint_verification_spec.md` §0 ("the adjoint algebra is already fully
   plumbed; P6 adds exactly two leaf transposes") under-counts.** That sentence is
   about the **apply/matvec** adjoint (`apply_transpose`), which IS fully plumbed.
   It does NOT cover the **solve** adjoint — the campaign per-leaf table
   (`adjoint_sn_transport.md:89`) correctly names `(L+C)†⁻¹` "the ONE genuine gap."
   The A0 spec has no intrinsic `solve_transpose` gate; **that is this document.**
3. Code matches the A3 description in every other respect. No code-level
   discrepancy found.

---

## 1. Claim-layer + pillar gate table (vv §1.5, MANDATORY)

| Gate | Claim layer | Pillar | Reference (structurally independent of the SUT) | Why it proves the layer |
|------|-------------|--------|--------------------------------------------------|-------------------------|
| **G1** round-trip `solve_transposeᵀ∘apply_transpose = I` | structural (operator-inverse identity) | closed-form (algebraic) | the identity `M⁻ᵀMᵀ = I` itself; `apply_transpose` independently pinned by `test_g_adjoint_reciprocity` + `test_scattering_adjoint` dense oracle | inverse∘operator = identity is the defining inverse law |
| **G2** dense `(L+C)ᵀ⁻¹` oracle | flux-shape (the solved field) | closed-form | `np.linalg.solve(M.T, b)`, **M = matrix of the FORWARD `apply`** (built by applying `(L+C).apply` to flat basis vectors) — shares NO code with the reverse-walk | M⁻ᵀb is the exact transpose-inverse; M built from forward apply only ⟹ independent of `solve_transpose` AND `apply_transpose` |
| **G3** full-loss G-reciprocity `(L+C−S−B)` | model (adjoint operator) | closed-form (reciprocity identity) | the independent `_g_inner` (`test_g_adjoint_reciprocity.py:154`); both halves independent | `⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G` is the defining G-adjoint law |
| **G4** Mode-11 routing sentinel | structural (call-graph) | — (in-process wrap counter) | monkeypatch counter on `InvertibleOperator.solve_transpose` / `.solve` | proves `.H.solve` executes the reverse solve, not the forward |
| **G5** `.H.solve` value round-trip | flux-shape (metric-conjugated) | closed-form | `(L+C).H.apply((L+C).H.solve(b)) = b` | `A*(A*)⁻¹ = I`; routing `.H.solve` to forward `solve` breaks it O(1) for non-symmetric A |
| **G6** capability flip + forward canaries | structural + regression | — | the capability-set laws + the named 0-ULP canaries | a NET-NEW `CAP_SOLVE` on `(L+C).H` only; forward `solve` untouched |

**Pillar discipline:** every row is closed-form / algebraic-identity / call-graph.
**No MMS** (correct — MMS proves neither operator-inverse identities nor
eigenvalues; vv §pillars). **No eigenvalue claim in A3** (the adjoint-keff gate is
A4/P1.3 downstream; A3 stops at the inner-solve correctness). **Structural
independence:** G2's reference is built ENTIRELY from the forward `apply` — it never
calls `apply_transpose` or `solve_transpose`, so it cannot share a reverse-walk bug
(the strongest available ground; vv L11/§structural-independence). G1's
`apply_transpose` is itself independently pinned, so G1+G2 together close the
"both reverse-walk methods copied the same bug" hole that G1 alone could not.

---

## 2. Config-blindness audit (vv §0.6) — the reverse-DAG-specific traps

The convenient config nulls the EXACT term the reverse-DAG solve is most likely to
get wrong. Each trap, with the activating config:

- **Symmetric `(L+C)` hides forward-DAG-vs-reverse-DAG.** If `(L+C)ᵀ = (L+C)`, then
  `solve_transpose = solve` and the "walked the forward DAG instead of reverse"
  mutation is invisible. `(L+C)` is non-self-adjoint in general (streaming is
  non-symmetric), but to give the mutation O(1) teeth use **heterogeneous
  (cell-varying σ_t) AND ≥2G per-group-varying σ_t** so the operator is strongly
  non-normal. Reuse the `_make_slab(ng=2)` per-group-varying `sig_t`
  (`test_g_adjoint_reciprocity.py:88-90`); make it 2-region heterogeneous.
- **Slab NULLS the μ-reversal (the curvilinear mirror).** The "+μ instead of −μ"
  bug lives ONLY in the curvilinear seed adjoint (`mirror[gd]`, the
  `quad.reflection_index("x")` permutation at `__init__.py:2658,2758`), gated by
  `if curvature != "cartesian"` (`:2751`). On a slab that branch never runs — slab
  is the **degenerate curvilinear case** (lessons L1). **⟹ the SPHERE leg is
  MANDATORY** for the μ-reversal mutation; a slab-only `solve_transpose` gate proves
  NOTHING about the curvilinear adjoint.
- **Flat / zero boundary NULLS the boundary swap.** The reverse sweep reads OUTFLOW
  cotangents and writes INFLOW cotangents (`__init__.py:2677-2695`). It is exercised
  only when `b.boundary ≠ 0`. **⟹ use `_random_composite`** (random bulk AND random
  boundary, `test_g_adjoint_reciprocity.py:173`) for the round-trip / oracle input —
  never a bulk-only or flat field.
- **Symmetric `SigS` hides the S† group-transpose (G3 only).** S† transposes the
  group-transfer `g↔g'`; symmetric `SigS` ⟹ S†=S, invisible (ERR-002 family, vv
  anti-#2/#6). **G3's config MUST use ASYMMETRIC `SigS`** (strong downscatter,
  ~zero upscatter). The bare `(L+C)` (G1/G2) has NO group coupling (collision is a
  per-group diagonal), so group asymmetry there means only per-group-varying σ_t;
  the group-TRANSFER transpose lives in S and is G3's job.
- **1G is degenerate** (vv §1-group): SigS, σ scalars ⟹ S†=S and the group axis
  collapses. **≥2G everywhere.**

**Minimum activating configs (the binding contract):**

- **G1/G2 (intrinsic `solve_transpose`):** VACUUM BC (so `(L+C)` is cleanly
  single-pass invertible — see §0.3), **heterogeneous 2-region σ_t, ≥2G
  per-group-varying**, on **BOTH a slab AND a sphere**. Slab = the tight exact leg
  (catches forward-DAG); sphere = the μ-reversal leg (catches the mirror).
  `_random_composite` input (non-zero boundary).
- **G3 (full-loss reciprocity):** a real solver mixture with **asymmetric `SigS`,
  ≥2G** (2G AND 4G), slab + sphere; `_random_composite` ψ,φ.

> Use **VACUUM**, not the reflective `_make_slab`/`_make_sphere` defaults
> (`test_g_adjoint_reciprocity.py:84,100`), for G1/G2: vacuum makes the bare
> `(L+C)` a clean lower-triangular system so the round-trip and dense oracle are
> exact single-pass (slab) / seed-only-iterated (sphere), with no reflective
> fixed-point iteration confound. Reflective BC belongs to the separate `−B` leaf
> (G3, where the metric reciprocity already handles it).

---

## 3. Deliverable 1 — the intrinsic `solve_transpose` gate (the keystone)

**New file:** `tests/sn/operators/test_loss_transpose_solve.py`. Reuse the mesh
builders + `_random_composite` from `test_g_adjoint_reciprocity.py` (lift them to a
shared helper, or import — do NOT re-derive). `@pytest.mark.foundation` (algebraic
identity, no theory `:label:`). **vv Mode-8:** `np.testing.assert_*` / `pytest.fail`
only (fire under `-O`) — NO bare `assert`.

### G1 — round-trip: `solve_transpose` exactly inverts `apply_transpose`

- **Claim layer / pillar:** structural / closed-form (operator-inverse identity).
- **Config:** vacuum slab (ng=2, het, per-group-varying σ_t) AND vacuum sphere
  (same). `_random_composite` φ (non-flat bulk + boundary).
- **Gate:** both directions —
  - `solve_transpose(apply_transpose(φ)) ≈ φ` (always well-posed: `apply_transpose(φ)`
    is in range by construction).
  - `apply_transpose(solve_transpose(b)) ≈ b` for `b = _random_composite` (well-posed
    for vacuum, where `(L+C)ᵀ` is full-rank).
  - Slab: single call. Sphere: thread the iterate (`initial_guess=`) to convergence,
    mirroring `test_invertible_operator.py:957-968` (the Carlson-seed loop, break
    floor `1e-14`).
- **Tolerance:** slab `rtol=1e-11, atol=1e-12` (exact single-pass); sphere
  `rtol=1e-9` (seed-convergence floor). **NOT `inner_tol`-scaled.**
- **Mutations (each MUST redden under `-O`; monkeypatch in-process — NEVER `git
  checkout` an uncommitted file):**
  1. **Forward-DAG order** — make `solve_transpose` walk the cells in sweep order
     instead of `reversed(...)`. ⟹ it computes `(L+C)⁻¹` not `(L+C)⁻ᵀ` ⟹
     `apply_transpose(solve_transpose(φ)) = (L+C)ᵀ(L+C)⁻¹φ ≠ φ` for the non-symmetric
     (het, 2G-asym) operator → RED. *(Invisible on a symmetric operator — hence the
     het+2G config.)*
  2. **+μ instead of −μ (forget the mirror)** — `mirror[gd] → gd` in the curvilinear
     seed adjoint. ⟹ RED on the **sphere** leg only (slab nulls it). This is the
     config-blindness keystone: slab stays GREEN under this mutation, sphere REDS.
- **Why G1 is not vacuous:** `apply_transpose` is independently correct
  (`test_g_adjoint_reciprocity` pins `(L+C−B)ᵀ`); so `solve_transpose ∘ apply_transpose
  = I` ⟹ `solve_transpose` is the genuine inverse. But G1 alone could pass if
  `solve_transpose` were built by copying `apply_transpose`'s walk with a SHARED bug
  the round-trip can't see → **G2 closes that hole.**

### G2 — dense `(L+C)ᵀ⁻¹` oracle (the structural keystone)

- **Claim layer / pillar:** flux-shape / closed-form.
- **Config:** a **tiny** vacuum mesh (slab `nx=4, ng=2` GL-S4 → ~32 bulk + trace
  DOFs; sphere likewise). het + per-group-varying σ_t.
- **Reference construction (forward-`apply`-only — the structural independence):**
  ```python
  def _dense_forward_matrix(LC, template):
      """M[:, j] = (L+C).apply(e_j) on the flat composite — FORWARD apply only."""
      n = template.to_flat().size           # full_field.py:337 (bulk.ravel ⊕ boundary)
      M = np.empty((n, n))
      for j in range(n):
          e = np.zeros(n); e[j] = 1.0
          M[:, j] = LC.apply(FullField.from_flat(e, template)).to_flat()  # :355
      return M
  ```
  - **Sub-check (the ERR-061 frame-consistency backstop):**
    `matrix(apply_transpose) ≈ M.T` to `rtol=1e-12` — confirms `apply_transpose` is the
    Euclidean transpose of `apply`, the oracle's premise. (Build `matrix(apply_transpose)`
    the same way with `LC.apply_transpose`.)
  - **The oracle:** `x_ref = np.linalg.solve(M.T, b.to_flat())`.
- **Gate:** `solve_transpose(b).to_flat() ≈ x_ref` for `b = _random_composite`.
  Slab: single call. Sphere: iterate-threaded.
- **Tolerance:** slab `rtol=1e-10, atol=1e-12` (sweep-solve vs LU-solve of the same
  system; agreement ~κ(M)·ULP, tiny well-conditioned mesh); sphere `rtol=1e-9`.
- **Mutations:** the SAME two as G1 → RED (forward-DAG: `solve_transpose(b)` matches
  `np.linalg.solve(M, b)` not `np.linalg.solve(M.T, b)`, RED for non-symmetric;
  μ-reversal: sphere RED). PLUS a third unique to G2 confirming structural
  independence: a wrong reverse-walk that happens to satisfy `solve_transpose ∘
  apply_transpose = I` (both miswalk) STILL REDS against `M.T` (M never touched the
  reverse walk).
- **Why G2 is the keystone:** the reference is built from `apply` alone. It is the
  ONLY gate that can catch a bug shared between `apply_transpose` and `solve_transpose`.

---

## 4. Deliverable 2 — confirm/refine P1.1 (full-loss `L+C−S−B` reciprocity)

**Extend** `tests/sn/operators/test_g_adjoint_reciprocity.py` (do NOT fork — it owns
the independent `_g_inner` machinery and the L11 wrong-metric control). The existing
`_loss_operator` (`:185`) builds `L + C − B`; P1.1 extends it to the full
`A_loss = L + C − S − B` now that S† is live (`15185e5`).

### What changes (exactly)

- Add `_loss_operator_with_scattering(solver) → (L+C) − S − B`. The existing
  builders use `placeholder_materials` (no scattering) — S needs a **real** mixture
  with **asymmetric `SigS`, ≥2G**. Reuse the `_mix(p0, p1)` asymmetric-block pattern
  from `test_scattering_adjoint.py:61-69` (or `make_mixture` with asymmetric `sig_s`),
  build an `SNSolver`, take `S = solver.scattering_op`. **S lifts to the FullField
  composite** (`S.apply(FullField) → FullField`, the `@overload` at
  `scattering.py:1659-1665`), so `(L+C) − S − B` is a valid `FullFieldSpace`
  `OperatorSum`; `loss_op.H` propagates through `OperatorSum.adjoint` to each leaf's
  `.H` (S now advertises `CAP_APPLY_TRANSPOSE`).
- **Gate:** `⟨A_loss ψ, φ⟩_G = ⟨ψ, A_loss.H φ⟩_G` with the independent `_g_inner`
  (anti-R1: NOT the production metric), random non-flat `_random_composite` ψ,φ.
- **Per-group, NOT weight-summed (vv L27):** the existing test sums `_g_inner` over
  all (N, ng, spatial). A scattering group-transfer transpose bug lives in the GROUP
  axis and could be masked in the weight sum. Add a per-group breakdown: run the
  reciprocity with **φ a one-hot in each group g** (`φ_g` supported only in group g),
  asserting `⟨A_loss ψ, φ_g⟩_G = ⟨ψ, A_loss.H φ_g⟩_G` for EACH g — a g'→g transpose
  swap then shows in the g-restricted residual. (The leaf-level per-group/per-ordinate
  S† reciprocity is already covered by `test_scattering_adjoint.py::test_S_euclidean_reciprocity`;
  G3's new claim is that S COMPOSES into the full-loss **G-adjoint** via `OperatorSum.H`
  WITH the metric.)
- **Claim layer / pillar:** model / closed-form (reciprocity identity).
- **Tolerance:** `rel < 1e-12` (algebraic; matches the existing gate `:222`).
- **Config:** Mixture asymmetric `SigS` 2G AND 4G; slab + **sphere** (the sphere leg
  proves the curvilinear `Lᵀ` composes with S† under the metric).
- **Mutations:**
  1. **S → S in the adjoint** (drop the group-transpose): monkeypatch
     `ScatteringOperator.apply_transpose → ScatteringOperator.apply`. ⟹ the composite
     reciprocity breaks O(1) with asymmetric SigS → RED. *(Invisible with symmetric
     SigS — hence the config.)*
  2. **`−S` vs `+S` sign error** in the loss posing → RED.
- **Why it is not redundant with G1/G2:** G1/G2 gate the bare `(L+C)` reverse SOLVE;
  G3 gates the full-loss reverse APPLY/adjoint (the `−S` leaf + metric composition).
  Different operators, different claim. G3 is also the gate that proves S† is
  metric-correct inside the full loss (the leaf S† test is Euclidean only).

> **Scope note:** G3 verifies the adjoint **matvec** of the full loss. It does NOT
> need `solve_transpose` (the full `A_loss = L+C−S−B` is NOT sweep-invertible — S
> couples groups; only the bare `(L+C)` has `CAP_SOLVE`). G3 is in A3's scope because
> the brief asks to confirm P1.1 "now that S† exists"; it is logically a matvec gate
> that A3 unblocks, feeding A4's `solve_transpose`-based daggered eigenvalue posing.

---

## 5. Deliverable 3 — `.H.solve` routes to `solve_transpose` (Mode-11 + value)

`(L+C).H.solve` is the daggered inner solve the A4 eigenvalue driver consumes. Two
gates — the Mode-11 split (value catches it on a non-symmetric operator; the
sentinel catches it structurally regardless). Both in `test_loss_transpose_solve.py`.

### G4 — Mode-11 routing sentinel (the new private reader must EXECUTE)

- **Claim layer / pillar:** structural / in-process wrap counter (vv Mode-11
  sharpening — wrap the internal call, strictly stronger than a file-write probe).
- **Gate:** monkeypatch `InvertibleOperator.solve_transpose` to increment a counter
  then call the real method; monkeypatch `InvertibleOperator.solve` (the FORWARD)
  to increment a SEPARATE counter. Call `(L+C).H.solve(b)`. Assert:
  - `solve_transpose` counter `> 0` (the reverse solve IS executed), AND
  - forward `solve` counter `== 0` (the forward solve is NOT touched by `.H.solve`).
- **Why not vacuous:** a bit-identity-preserving regression that wired `.H.solve`
  to the forward `solve` would give a DIFFERENT value than G5 expects — but on a
  near-symmetric fixture G5's value gap could be small; G4 fires regardless of how
  close `(L+C)ᵀ ≈ (L+C)`. Pair G4 (structural) + G5 (value) — the exact split
  Mode-11 exists for.
- **Mutation:** wire `.H.solve` to call `inner.solve` (forward) → G4's
  `solve_transpose counter > 0` assertion REDS (counter stays 0). `-O`-safe
  (counter assertion is `np.testing`/`pytest.fail`, a function call).

### G5 — `.H.solve` value round-trip (the metric-conjugated inverse)

- **Claim layer / pillar:** flux-shape / closed-form.
- **Gate:** `(L+C).H.apply((L+C).H.solve(b)) ≈ b` for `b = _random_composite`
  (`A*(A*)⁻¹ = I`). Slab exact (`rtol=1e-11`), sphere iterate-threaded (`rtol=1e-9`).
- **Mutation (the routing bug at the VALUE level):** `.H.solve` routes to forward
  `inner.solve` instead of `inner.solve_transpose` ⟹ `.H.solve(b) = G⁻¹(L+C)⁻¹G b`
  not `G⁻¹(L+C)⁻ᵀG b` ⟹ `(L+C).H.apply(that) = G⁻¹(L+C)ᵀ(L+C)⁻¹G b ≠ b` for
  non-symmetric `(L+C)` → RED. (Complements G4: G5 needs the het+2G asymmetry; G4
  fires on any config.)
- **Metric-direction mutation (a bug A3's metric wiring can introduce):** swap
  `apply_metric` ↔ `apply_inverse_metric` in `_AdjointOperator.solve` (or apply the
  metric on the wrong side). For the symmetric SN metric `G_V=G_W=G` this still
  reds via the round-trip when the metric is non-trivial (the trace `|Ω·n|·w_n`
  weighting varies across ordinates) — but the SHARPER catcher is the G3 full-loss
  reciprocity machinery's L11 wrong-metric control. Note in the spec that the
  metric direction is gated transitively by G5 + the existing
  `test_wrong_trace_metric_breaks_reciprocity` (`test_g_adjoint_reciprocity.py:268`)
  pattern; if the implementer wants a direct catcher, add an `_AdjointOperator.solve`
  metric-reciprocity probe `⟨(L+C).H.solve(b), c⟩_G` consistency.

---

## 6. Deliverable 4 — forward-safety: the capability flip + the 0-ULP canaries

A3 ADDS `solve_transpose` + `_AdjointOperator.solve`; it does **not** touch the
forward sweep `solve`, `apply`, or `apply_transpose`. The forward path must stay
bit-identical.

### G6a — the capability flip (a NET-NEW `CAP_SOLVE` on `(L+C).H` ONLY)

A3 changes `_AdjointOperator.__init__` to conditionally add `CAP_SOLVE`. This touches
EVERY `.H` in the codebase, so gate the propagation precisely:

- **NEW positive pin:** `(L+C).H` advertises `CAP_SOLVE` (where `(L+C)` is an
  `InvertibleOperator`). Add to `tests/sn/operators/test_capability_survival.py`
  (alongside `test_l_plus_c_is_invertible_with_solve` `:127`).
- **NEW negative pins (forward-safety — the change must NOT over-propagate):**
  `S.H`, `F.H`, and a bare `L.H` (whose inner has NO `solve_transpose`) MUST still
  lack `CAP_SOLVE`. The wiring criterion must be "inner exposes the reverse solve"
  (a distinct capability/method check), NOT merely "inner has `CAP_SOLVE`" — confirm
  the implementer keys it on `solve_transpose` so only `InvertibleOperator.H` flips.
- **Audit (no existing pin to flip):** verified — no test asserts
  `_AdjointOperator`/`.H` lacks `CAP_SOLVE` (the `CAP_SOLVE not in` pins in
  `tests/numerics/test_operator.py:211,228` are `ZeroOperator` / apply-only leaves;
  `test_capability_survival.py:137` is the `(L+C−B)` **OperatorSum** which genuinely
  stays solve-less — UNRELATED to the adjoint, must STAY). So G6a is purely additive.

### G6b — the named forward canaries (MUST stay green / 0-ULP)

| Canary | What it pins | Why A3 must not move it |
|---|---|---|
| `tests/sn/operators/test_invertible_operator.py` | forward `LC.solve(LC.apply(ψ))=ψ`, `Q/Σ_t` recovery (slab/sphere/cyl), `ERR-049` | forward `solve` UNCHANGED |
| `tests/sn/spatial/test_sweep_cache.py` | the sweep dual-view bit-id (`rtol=1e-13`) | the sweep kernel is untouched |
| `tests/sn/operators/test_g_adjoint_reciprocity.py` | `(L+C−B)` G-adjoint reciprocity + L11 control | `apply_transpose` UNCHANGED; G3 ADDS the `−S` leg here, the existing rows stay green |
| `tests/sn/operators/test_scattering_adjoint.py`, `test_fission_adjoint.py` | S†/F† Euclidean reciprocity, capability flips | A3 does not touch S/F |
| `tests/sn/operators/test_capability_survival.py` | the existing `CAP_SOLVE` laws (incl. the `(L+C−B)` solve-less pin `:137`) | only G6a's additive pins change |
| the forward fixed-source + keff suites (`tests/sn/solve/`, `test_kinf_homogeneous`) | the production SI/Krylov/eigenvalue path | the inner forward solve is the SoT, untouched |

- **Bit-identity discipline:** A3 introduces NEW reverse-solve values (re-baseline,
  not bit-id — campaign §"Re-baseline, NOT bit-identity"). But the FORWARD canaries
  above must stay **0-ULP/green**; their RED is the signal that A3 leaked into the
  forward path (a scope violation). Run the canonical `-O -m "not slow"` gate over
  `tests/sn/operators/ tests/sn/spatial/ tests/sn/solve/` post-A3 and confirm green.

---

## 7. Failure-mode table rows (self-improvement — file BEFORE delivering)

New rows for the `vv-principles` failure-mode / `numerical-bug-signatures` table
(the reverse-DAG-solve modes not yet catalogued; promote to ERR-NNN if a real bug is
caught in implementation):

| Failure mode | Test-design row that catches it |
|--------------|---------------------------------|
| **Reverse-DAG solve walks the FORWARD DAG order** (computes `(L+C)⁻¹` not `(L+C)⁻ᵀ`) | G2 dense `M.T` oracle + G1 round-trip; config MUST be non-symmetric (het + ≥2G per-group-varying σ_t) — invisible on a symmetric operator |
| **μ-reversal dropped in the reverse solve** (forget the curvilinear `mirror` permutation) | G1/G2 on the **SPHERE** leg (slab nulls it — the degenerate curvilinear case) |
| **`.H.solve` routes to the forward `solve`** (no reverse) | G4 Mode-11 sentinel (structural, any config) + G5 value round-trip (non-symmetric) |
| **`_AdjointOperator` over-propagates `CAP_SOLVE`** (gives `S.H`/`F.H`/`L.H` a phantom solve) | G6a negative pins (key the flip on `solve_transpose`, not `CAP_SOLVE`) |
| **Metric applied on the wrong side in `.H.solve`** (`apply_metric`↔`apply_inverse_metric` swap) | G5 round-trip (non-trivial trace `|Ω·n|·w_n`) + the L11 wrong-metric control pattern |

---

## 8. Gate dependency DAG + the done criterion

```
G1 round-trip ─┐
               ├─→ (solve_transpose intrinsically correct)
G2 dense M.T ──┘            │
                            ├─→ G4 Mode-11 sentinel ─┐
G3 full-loss reciprocity ───┤                        ├─→ (.H.solve correct ⟹ A4 unblocked)
  (S† composes, metric)     └─→ G5 .H.solve value ───┘
G6a capability flip + G6b forward canaries  (run throughout — forward stays green)
```

**A3 is DONE only when:** G1, G2, G3, G4, G5, G6a all GREEN; every named mutation
reddens its named gate under the canonical `-O` invocation (mutate in-process via
monkeypatch — NEVER `git checkout` an uncommitted file, `.claude/rules/process-discipline.md`);
G6b forward canaries 0-ULP/green; and the SPHERE leg is present in G1/G2/G3 (a
slab-only A3 is REJECTED — it is blind to the μ-reversal, the single hardest term in
the reverse-DAG curvilinear solve).

**Per the standing discipline (AGENT.md §0.5):** a green gate that cannot red is
worse than no gate. Before crediting each gate, name the mutation that reddens it and
confirm it fires under `-O`. The keystone is **G2** (forward-`apply`-only reference —
the only structurally-independent ground for the reverse solve); G1 is its cheap
consistency companion; the SPHERE leg is the config that makes the μ-reversal term
un-hideable.
```
