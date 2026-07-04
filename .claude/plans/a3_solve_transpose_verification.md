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

---

# EXTENSION (2026-07-04, Phase 2.5 P0-B, test-architect) — the #280 forward-refactoring carve, post-assembly world

§§1–8 above remain the base. This extension supersedes §0.1/§0.2 (paths + operator
names) and adds the reverse-SCAN, assembled-oracle, 2.5a-relocation, and
scaffold-order content. The base G1–G6 gate MATH survives; only the operator
vocabulary and the reverse-loop→reverse-scan assumption change.

## §9. Load-bearing reconciliation — the base spec predates two merges

**§9.1 Operator vocabulary is RETIRED (the #226 taxonomy merged).** The base spec's
`InvertibleOperator.solve_transpose` + `CAP_SOLVE_TRANSPOSE` + `_AdjointOperator.solve`
+ `CAP_SOLVE` do NOT exist in the current tree. Map every gate onto the LIVE surface:

| Base-spec name (retired) | CURRENT surface (verified 2026-07-04) | file |
|---|---|---|
| `InvertibleOperator.solve_transpose` (the reverse-scan primitive) | **`SweepOperator.apply_transpose`** — `SweepOperator = (L+C).inverse()`; the reverse-scan IS the inverse operator's transpose matvec `(L+C)⁻ᵀb` | `orpheus/sn/operators/sweep_operator.py` |
| `CAP_SOLVE_TRANSPOSE` string tag | **`SweepOperator.is_adjointable` → True** (today inherits base `False`; the eager `.H` gate raises `MissingAdjoint` until it flips) | `orpheus/numerics/operator.py:538` |
| `_AdjointOperator.solve` + `CAP_SOLVE` | **NOTHING NEW** — `A.inverse().H.apply(b) = G⁺·SweepOperator.apply_transpose(G·b)` falls out of the EXISTING `_AdjointOperator.apply` (already routes `inner.apply_transpose`, `operator.py:1127`) FOR FREE once `SweepOperator.apply_transpose` lands | — |
| the `.H.solve` coherence | **`_AdjointOperator.inverse()` + `is_invertible`** closing the swap law `A.H.inverse() ≡ A.inverse().H` (predicate honesty per the #280 comment: `inner.is_invertible and inner.inverse().is_adjointable`) | the #280 issue comment |

`InverseWrapMixin`'s docstring already names this division ("the adjoint-inverse is
the #280 family, deferred … a reverse-DAG `sweep_transpose` for the direct sweep").
**All gates below are written against the live surface. Do NOT resurrect string tags.**

**§9.2 Current paths (Phase-2a relocation — base-spec line numbers are dead).**
- Production `orpheus/sn/loss_representation/__init__.py`: `loss_action_transpose`
  (1-D adj matvec) `:2777`; `_apply_walk` (1-D fwd matvec) `:2426`; `_run` (1-D solve
  SCAN) `:2994`; `reversed(cells)` `:2903`; `mirror = quad.reflection_index("x")`
  `:2833`; `angular_adjoint` route `:2943`. Multi-D Cartesian adjoint still raises
  (`:2822`).
- Kernel `orpheus/transport/spatial/`: `cell_balance_for_streaming`
  `cell_balance.py:120`; `_cartesian_streaming_diagonal` `diamond.py:363`;
  `residual_kernel_batch` `diamond.py:494` (DD) / `linear_discontinuous.py:646` (LD).
- The 2.5 layout estate `orpheus/sn/spatial/` = {pole_angular_closure,
  psi_half_angle_seed, sweep_cache, scan, pairing} (R9: name kept this phase; 2.5
  owns the layout decision).

## §10. 2.5a apply-loop unification — bit-identity of a RELOCATION (both orientations)

2.5a folds `{_apply_walk (fwd matvec), loss_action_transpose (adj matvec)}` into ONE
orientation-parametrized per-cell loop frame, bit-identical for BOTH orientations.

**§10.1 Refreshed §6A face→canary table (current paths, verified to exist):**

| Face | Canary (current path) | Tol |
|---|---|---|
| `sweep` 1-D solve | `tests/sn/sweep/core/test_sweep_cache.py` (moved from `tests/sn/spatial/`) + slab snapshots | rtol=1e-13 / 1e-12 |
| `loss_action` 1-D fwd matvec (`_apply_walk`) | `tests/sn/operators/test_removal_form_matvec_sweep.py::test_invertible_apply_is_M_of_C_sigma_bit_identical` (slab/sphere/cyl/cart2d 2G) + `test_invertible_operator.py` Q/Σ_t + `tests/sn/sweep/core/test_cell_kernel_batch.py` | **array_equal (0-ULP)** |
| `loss_action_transpose` 1-D adj matvec | `test_removal_form_matvec_sweep.py::test_invertible_apply_transpose_is_M_transpose_of_C_sigma_bit_identical` (slab/sphere/cyl 2G, NOT cart2d) + `tests/sn/sweep/core/test_phase_c_gates.py::test_apply_apply_transpose_reciprocity_under_sweep_frame` (sphere/cyl) + `test_g_adjoint_reciprocity.py` (slab/sphere/cyl 1g+2g) | array_equal / rel<1e-12 |
| one-walk structure (multi-D only) | `tests/sn/operators/test_one_octant_walk.py` — `_interior_walk` spy + AST tripwire `{is_solve,is_apply,is_matvec}` | — |

**§10.2 Sufficiency verdict — the `array_equal` canaries are SELF-REFERENTIAL →
necessary but NOT sufficient for a relocation.** The removal-form gates compare
`op.apply[_transpose]` against the SAME `loss_action[_transpose]` on a fresh
`StreamingOperator(σ_r)` instance (`test_removal_form_matvec_sweep.py:274-286,
314-320`) — an override-owns-it vs leaf-sum-leak discriminator, NOT a value oracle.
A relocation that moves BOTH the SUT path AND the reference path leaves them green
even if values shifted (the twin/Mode-11 hazard). **The sufficient addition (S0.5
pre-carve): a FROZEN pre-carve baseline** via
`tests/sn/regression/_regression_assert.py::assert_regression`
(`--capture-baseline`, `kind="direct"` → nulp(reduction_depth)) of the fwd+adj 1-D
matvec VALUES (slab/sphere/cyl, 2G, `_random_composite` input) — structurally
independent because captured BEFORE the relocation. The self-referential canaries
stay as the override-not-leak guard.

**§10.3 The new shared 1-D loop frame needs its OWN structural pin (none exists
today).** `_OneDimScanWalk` has no `test_one_octant_walk` analog (its three bodies
are separate verbatim relocations). 2.5a mints
`tests/sn/sweep/core/test_one_dim_loop_walk.py`:
1. **Runtime spy (Mode-11-sharpened WRAP):** monkeypatch counter on the ONE shared
   loop method; call the forward matvec AND the adjoint matvec; assert the counter
   fires for BOTH orientations. `np.testing`/`pytest.fail`, never bare `assert`.
2. **AST tripwire — the ORIENTATION sibling of the `is_solve` rule:** `ast.parse`
   the shared frame's source; fail if `is_adjoint`/`is_forward`/`is_transpose`/
   `is_reverse` appear as a real identifier (Name/Attribute/arg, not docstring) —
   demand an orientation-carrying OBJECT (the `_ApplyOperands`/`_SolveOperands`/
   `_SweepEmit` discipline's sibling). Verified: none of those identifiers exist in
   the module today — the tripwire starts clean.

**Config (§0.6 discipline):** het 2-region σ_t + ≥2G per-group-varying
(non-symmetric operator); `_random_composite` NON-ZERO boundary (else the boundary
in↔out swap is nulled); slab + **sphere** + cyl (slab NULLS the μ-reversal mirror).

## §11. G1–G6 refresh for the reverse-SCAN design

The base G1 (round-trip) and G2 (dense `Mᵀ` oracle) survive **unchanged in math** —
reframe `solve_transpose(b)` → `A.inverse().apply_transpose(b)` =
`SweepOperator.apply_transpose(b)`. Three additions:

**§11.1 (a) The assembled-`Mᵀ` oracle — a NEW structurally-independent reference
(Cartesian ONLY).** Phase 2 gives `M` = the per-(ordinate,group) sparse block of
`L+C`, walk-order LOWER-triangular:

```python
permuted = assemble_ordinate_blocks(sn_mesh, n)[g].as_matrix()[np.ix_(order, order)]
x_ref = solve_triangular(permuted.T, b[order], lower=False)  # LAPACK dtrtrs back-sub
# un-permute: SweepOperator.apply_transpose(b)[order] ≡ x_ref  (rtol 1e-11, expect ~6e-16)
```

(`SparseAssembledOperator.apply_transpose` — exact CSR `M.T @ x` — is the matvec
analog usable as a slab canary for `loss_action_transpose` itself.)

**L11 classification (honest):** the assembled blocks are extracted from
`residual_kernel_batch` UNIT PROBES — independent of the REVERSE-walk code, but
sharing the forward KERNEL.
- **Catches (that G2's dense-apply oracle does not):** a wrong TRANSPOSED-SCAN
  COEFFICIENT (a'/b' of the reversed affine recurrence — LAPACK `dtrtrs` shares NO
  code with the ORPHEUS scan) + the exact triangularity certificate + the
  object-level "reverse-substitution IS `sweep_transpose`" discharge (the transpose
  analog of the #284 discharge).
- **Misses (that G2's dense-apply oracle catches):** everything CURVILINEAR —
  assembly refuses curvilinear, so **the sphere μ-reversal/mirror bug stays on
  G2's dense forward-apply keystone.**
- **Gate row:** YES — instantiated on the **DD SLAB** this phase (2.5's transpose
  scope is 1-D DD: LD-slab and multi-D are typed defers, so the LD/2-D oracle rows
  activate only if those deferrals ever close). Fixtures: reuse
  `test_assembly_mode`'s het + non-uniform-h + 2G + vacuum meshes.

**§11.2 (b) Reverse-SCAN failure modes the base (reverse-LOOP) spec did not cover:**

| NEW failure mode | Catcher |
|---|---|
| Wrong transposed scan coefficients (a'/b' of the reversed affine recurrence — scan ≠ loop derivation) | G2 dense-`Mᵀ` + **assembled-`Mᵀ`** (LAPACK trtrs), non-symmetric het+≥2G |
| Two-denom seam: the reverse-scan rides ÷V `residual_kernel_batch` (apply form) instead of ×V `affine_scan_coefficients` (the `_run` form) — same slab denom two ways (#242) | G1 round-trip (denom mismatch breaks `apply_transpose∘sweep_transpose = I`) + G2 |
| Curvilinear per-μ-level ordinate-coupled scan REVERSAL + `angular_adjoint` threading + Carlson mirror-seed cross-sweep coupling | G1/G2 **SPHERE** leg (slab NULLS — degenerate curvilinear) — unchanged keystone |

**§11.3 (c) Mode-12 audit — invariance groups + joint blindness.** The reverse-solve
gates measure the OBJECT (the full solved field), not a scalar functional — preserve
this:
- G1's invariance group: any two-sided inverse pair — a bug corrupting BOTH
  `apply_transpose` and `sweep_transpose` identically composes to I and is invisible
  (the base G1-alone hole; G2 closes it).
- G2 / assembled-`Mᵀ` invariance group: the map `M ↦ Mᵀ`; catches
  forward-DAG-instead-of-reverse, blind to a shared forward-kernel bug.
- **Jointly blind (G1+G2+assembled):** a bug in the SHARED forward kernel
  (`residual_kernel_batch` / `cell_balance_for_streaming` denom) moving `apply`,
  `apply_transpose`, AND `sweep_transpose` together. Closed by the FORWARD ground
  (removal-form array_equal + `test_invertible_operator` Q/Σ_t +
  `test_g_adjoint_reciprocity`) + the one-source teeth
  (`test_assembly_mode::test_teeth_shared_kernel_sign_flip_moves_all_three_modes`).
  Mode-11-adjacent twin-path coverage, not Mode-12; the reverse solve inherits
  forward correctness exactly as assembly (the 3rd consumption mode) did.
- **Hard Mode-12 prohibition:** NEVER credit `sweep_transpose` on `k*` or any
  norm/reaction-rate functional — `eig(Kᵀ) = eig(K)` by construction (the A4
  daggered eigenvalue carries ZERO adjoint mutation coverage on the shared
  spectrum). Every reverse-solve gate compares full fields.

## §12. #282-fix gate plan (CONDITIONAL — pending ruling R10)

**Applies only IF route (a)'s directness lands inside 2.5.**
- The existing characterization
  (`test_assembly_mode::test_282_characterization_spherical_seed_is_a_back_edge`)
  asserts the back edge EXISTS; when the fix makes the spherical walk genuinely
  triangular it REDS with its actionable rewrite message (the loud-flip contract —
  NOT xfail). That flip IS the fix's acceptance evidence.
- **Successor = a spherical G2 gate** (triangularity + LAPACK ≡ sweep on the
  AUGMENTED system). Structure: the μ=−1 zero-weight Carlson starting-direction
  seed rows (Hébert 3.432–3.435) sit FIRST in each μ-level's block; per-ordinate
  cell blocks march after in μ-increasing order; the current back edge (ordinate →
  seed at `mirror`) becomes a within-level FORWARD edge (seed → first-ordinate
  inflow) — `M_aug` block-lower-triangular in [seed-rows-per-level,
  ordinate-blocks] order.
- **Teeth (each must red under `-O`, monkeypatch):** (i) swap the augmented seed's
  coupling direction (feed the LAST ordinate instead of the starting direction) →
  the `triu(P·M_aug·Pᵀ,1)==0` leg reds; (ii) a sign flip in the Hébert 3.432–3.435
  starting-direction solve → the LAPACK ≡ sweep leg reds.
- This fix is what would unblock a CURVILINEAR assembled-`Mᵀ` oracle (§11.1),
  retiring the sphere's dependence on the dense-apply keystone alone.

## §13. Scaffold order (single ordered gate chain)

**G3 lands PRE-carve** — the full-loss `(L+C−S−B)` reciprocity (base Deliverable 2)
hardens the adjoint-matvec surface 2.5a rebuilds; S† is live (`15185e5`), unblocked
today. Canary, not proof of the new primitive.

| Step | Gate | Role |
|---|---|---|
| **S0 pre-carve** | **G3** full-loss `(L+C−S−B)` G-reciprocity — extend `test_g_adjoint_reciprocity.py` `_loss_operator` from `L+C−B` (asymmetric SigS ≥2G real mixture; per-group one-hot φ, vv L27; S lifts to FullField via `OperatorSum.H`). Slab + sphere. PLUS: the LD-slab transpose loud guard (FLAG-2 fix-now) + the `dag_walk_cell_indices ≡ dag_walk` pin. | pre-carve canary |
| **S0.5 pre-carve** | **Frozen baseline** `assert_regression --capture-baseline` of fwd `_apply_walk` + adj `loss_action_transpose` values (slab/sphere/cyl 2G, `_random_composite`) | pre-carve anchor (the relocation's reference) |
| **2.5a** | frozen-baseline array_equal (both orientations) + removal-form array_equal (override-not-leak) + NEW 1-D loop spy + AST tripwire (§10.3). Forward canaries stay green. | relocation proof (bit-id) |
| **2.5b** | **G1** round-trip + **G2** dense-`Mᵀ` (slab + **SPHERE** keystone) + **assembled-`Mᵀ`** (DD slab, §11.1) + §11.2 mutations RED | new-primitive proof (re-baseline, NOT bit-id) |
| **2.5c** | the #280-comment gates in `tests/sn/operators/test_inverse_adjoint_coherence.py`: Gate 1 forward-matvec G-reciprocity pin (`⟨A.apply(ψ), x⟩_G = ⟨ψ, b⟩_G` for `x = A.H.inverse().apply(b)` — never calls the transpose path); Gate 2 swap-law value `A.H.inverse().apply(b) ≡ A.inverse().H.apply(b)` (rtol 1e-12); Gate 3 round-trip `A.H.apply(A.H.inverse().apply(b)) ≈ b`. PLUS Mode-11 wrap sentinel (reverse scan EXECUTES; forward `solve` counter == 0) + mutations M-ADJ-swap (drop the `.H` in `_AdjointOperator.inverse()` → Gates 1+2 red) and M-ADJ-metric (skip `G⁺`/`G` → Gate 1 reds on sphere/cyl ONLY, green slab — the `.H`≠Euclidean discriminator) + predicate flips (`SweepOperator.is_adjointable` True; `InverseOperator`/`GreenOperator` STAY False) + `assert_type(A.H.inverse(), LinearOperator[D,C])` static pin (extend `_composition_algebra_return_type_static_pins`) | wiring proof (structural + value) |

## §14. Failure-mode table rows (append to base §7)

| Failure mode | Test-design row that catches it |
|--------------|---------------------------------|
| Reverse-SCAN wrong transposed scan coefficients (a'/b'; scan ≠ loop) | G2 dense-`Mᵀ` + assembled-`Mᵀ` (LAPACK trtrs), non-symmetric het+≥2G |
| Two-denom seam (reverse-scan rides ÷V apply form vs ×V `_run` form, #242) | G1 round-trip + G2 |
| Relocation moves `op.apply` AND its self-referential reference together | FROZEN pre-carve baseline (§10.2) — NOT the array_equal self-reference |
| Orientation encoded as a boolean flag | AST tripwire `is_adjoint`/`is_forward`/`is_transpose`/`is_reverse` + orientation OBJECT (§10.3) |
| Assembled-`Mᵀ` / G1 / G2 jointly blind to a shared forward-kernel bug | forward-correctness ground + the one-source teeth (§11.3) |
| `sweep_transpose` credited on `k*`/norm | Mode-12: object-level full-field gates ONLY (§11.3) |

## §15. Done criterion (extends base §8)

Phase 2.5's gate chain is DONE only when: G1, G2, assembled-`Mᵀ` (DD slab), G3, and
the 2.5c wiring gates are all GREEN; every §11.2/§14 mutation reddens its named gate
under `-O` (monkeypatch, never `git checkout` an uncommitted file); the frozen
baseline holds array_equal for both orientations (2.5a bit-id); the 1-D loop spy +
AST tripwire are green-on-clean; **the SPHERE leg is present in G1/G2/G3**
(slab-only REJECTED — blind to the μ-reversal). The keystone remains **G2's
dense-apply `Mᵀ` oracle on the SPHERE** (the only structurally-independent ground
where assembly cannot reach); the assembled-`Mᵀ` oracle is the Cartesian L2
cross-check that additionally catches a wrong transposed-scan coefficient.
