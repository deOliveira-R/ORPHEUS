---
name: issue-196-phase-g-step2-5-dd-polymorphism
description: Investigation into whether DiamondDifference's three geometry branches (slab, curvilinear, cylindrical-degenerate) can collapse into ONE polymorphic body. Identifies the minimal geometric-data superset, where the closure splits cleanly from the cell-balance algebra, and the one residual case that genuinely resists unification.
metadata:
  type: project
---

# DD polymorphism — characterising the three branches (Issue #196 Phase G Step 2.5)

User directive: "DD should be DD — geometry becomes DATA carried by `visit`, not a branch inside `DiamondDifference.update`."

This memo investigates whether that is achievable, what the unified body looks like, and what the obstructions are.

Source files inspected (all of them, in full):
- `orpheus/sn/spatial/cell_update.py` (Protocol + dataclasses)
- `orpheus/sn/spatial/diamond.py` (the three `_update_*` static methods, the three `_residual_*` methods, the 2-D `update_batch`)
- `orpheus/sn/spatial/cell_balance.py` (`cell_balance_terms`, `cell_balance_terms_degenerate`, `CellBalanceTerms` dataclass)
- `orpheus/geometry/reduced_operator.py` (the `StreamingTerms` factory — confirms which fields are `None` vs populated per geometry)

---

## 1. The 3-column-by-6-row branch comparison table

Notation: `st = visit.streaming_terms`; `μ ≡ st.abs_mu`; `V ≡ st.volume`; `Σt ≡ total_xs`; `q ≡ source` (already weight-normalised); `ψˢ_in ≡ upstream_state.spatial_upstream`; `ψᵃ_in ≡ upstream_state.angular_upstream`.

| Row | Slab (`_update_slab`) | Curvilinear (`_update_curvilinear`) | Cylindrical-degenerate (`_update_cylindrical_degenerate`) |
|---|---|---|---|
| **1. Inputs consumed** | `st.abs_mu`, `st.chord_length`, `Σt`, `q`, `ψˢ_in`. `ψᵃ_in` UNUSED. `visit.face_area_downstream` UNUSED. | `st.abs_mu`, `st.face_area_inner`, `st.face_area_outer`, `st.delta_A_over_w`, `st.alpha_in`, `st.alpha_out`, `st.tau_mm`, `st.volume`, `visit.face_area_downstream`, `Σt`, `q`, `ψˢ_in`, `ψᵃ_in`. | `st.delta_A_over_w`, `st.alpha_in`, `st.alpha_out`, `st.tau_mm`, `st.volume`, `Σt`, `q`, `ψᵃ_in`. `μ`, face areas, `visit.face_area_downstream`, `ψˢ_in` UNUSED. |
| **2. Geometric quantities** | `chord_length` (=`V` for unit cross-section), `μ`. Angular-redistribution coefficients (`α_in`, `α_out`, `ΔA/w`, `τ_mm`) ABSENT (=`None`). Face areas ABSENT. | All curvature fields populated: `A_inner`, `A_outer`, `A_total = A_inner + A_outer`, `A_downstream` (sweep-resolved), `ΔA/w`, `α_in`, `α_out`, `τ_mm`, `V`. | All curvature fields populated EXCEPT effectively `A_downstream` is unused (its multiplier `μ → 0`). |
| **3. Denominator (symbolic)** | `2μ + chord·Σt`  *(NB: no `V` factor — chord plays both roles since `V = chord` for slab; see Row 6.)* | `2μ·A_down + (ΔA/w)·c_out + Σt·V`, where `c_out = α_out / τ_mm`. | `(ΔA/w)·c_out + Σt·V`. Identical to curvilinear with the `2μ·A_down` term zeroed. |
| **4. Numerator (symbolic)** | Two-step: `ψ_out = a·ψˢ_in + 2·q/denom` with `a = (2μ − chord·Σt)/denom`, then `ψ_avg = ½(ψˢ_in + ψ_out)`. NOT in `(q + numer_upstream)/denom` shape — different operation order to preserve bit-identity. | `(q + numer_upstream) / denom` where `numer_upstream = μ·A_total·ψˢ_in + (ΔA/w)·c_in·ψᵃ_in`, `c_in = (1−τ)/τ·α_out + α_in`. | `(q + numer_upstream) / denom` where `numer_upstream = (ΔA/w)·c_in·ψᵃ_in`. Identical to curvilinear with the `μ·A_total·ψˢ_in` term zeroed. |
| **5. Spatial / angular closure** | Spatial: `ψ_out = a·ψˢ_in + 2q/denom` computed BEFORE `ψ_avg`. Angular: NONE → `outgoing_angular_state = None`. | Spatial: `ψˢ_out = 2·ψ_avg − ψˢ_in` (WDD). Angular: `ψᵃ_out = (ψ_avg − (1−τ)·ψᵃ_in)/τ` (Morel-Montry). | Spatial: NONE → `outgoing_spatial_flux = None`. Angular: same Morel-Montry as curvilinear. |
| **6. Degenerate-axis special handling** | N/A (this IS the slab branch). | N/A (this IS the non-degenerate curvilinear branch). | Triggered by `μ < 1e-15`. Drops `2μ·A_down` from denom, drops `μ·A_total·ψˢ_in` from numer, drops the spatial-closure step entirely. Otherwise identical to curvilinear branch. The cell has NO radial face flow on this ordinate. |

Two observations cascade from this table:
- **The slab branch is algebraically equivalent to a curvilinear branch with `A_in = A_out = 1`, `ΔA/w = 0`, no Morel-Montry closure, but is implemented with a DIFFERENT operation order** (Row 4 / Row 5). The slab path computes `ψ_out` first via the explicit recurrence, then averages; curvilinear computes `ψ_avg` first via the balance solve, then closes. They give the same answer in exact arithmetic but NOT bit-identically — different float intermediates. The module docstring's bit-identity contract is the reason this difference exists.
- **The cylindrical-degenerate branch is a strict sub-case of the curvilinear branch** (Row 3 / Row 4): it is the curvilinear formula evaluated at `μ = 0`, with the spatial-closure step skipped. There is no genuine algebraic break — only the spatial-closure output disappears because the upstream flux gives no information about a flux that never flowed.

---

## 2. Answers to the four design questions

### Q1. Can the three branches collapse into ONE expression?

**Algebraically YES — but with two caveats.**

The slab denominator `2μ + chord·Σt` is `2μ·1 + Σt·chord` — identical in shape to the curvilinear `2μ·A_down + (ΔA/w)·c_out + Σt·V` if you set `A_down = 1` (a "slab face area" of unity), `ΔA/w = 0` (no angular redistribution surface), and `V = chord`. The slab numerator `2q + (2μ − chord·Σt)·ψˢ_in = 2q + (2μ − Σt·V)·ψˢ_in`, divided by `denom`, equals the slab `ψ_out`; averaging with `ψˢ_in` gives `ψ_avg = (q + μ·ψˢ_in)/denom·... ` — and after the algebraic dust settles, this is bit-equivalent to the curvilinear shape `(q + μ·A_total·ψˢ_in)/denom` with `A_total = A_in + A_out = 1 + 1 = 2`. (Because `A_total·μ·ψˢ_in = 2μ·ψˢ_in` and the slab implicit `A_in = A_out = 1` gives a chord of 1 on each side.)

**The two caveats:**

1. **Bit-identity vs. algebraic-identity.** Today's slab branch builds `ψ_out` from `a·ψˢ_in + 2q/denom` and then averages — this preserves bit-equality with the legacy cumprod sweep (sweep.py:117-123, :208-222). The unified `(q + numer_upstream)/denom` shape is algebraically the same answer but different float intermediates → 1-ULP regression on `np.array_equal` tests. **The bit-identity hand-calc tests will fail** unless the slab tests are relaxed to `np.allclose(rtol=1e-13)`. This is an explicit project rule documented in `diamond.py` lines 173-194 ("The bit-identical contract — non-negotiable"). User decision needed: is the unification worth re-baselining the slab bit-identity tests against the new operation order?

2. **The cylindrical-degenerate case requires CONDITIONAL CLOSURE.** The cell-balance algebra is one formula across all three branches (with the appropriate `None`-defaulted fields participating as zeros). The SPATIAL closure (`ψ_out = 2·ψ_avg − ψˢ_in`) is not meaningful in the degenerate axis case — there is no downstream face to write into. Inside the cell-balance algebra this is naturally handled: `2μ·A_down = 0`, so `denom` and `numer_upstream` are well-defined. But the closure step still wants to produce a `ψˢ_out` that the caller mustn't use. The clean answer (see Q4) is to make the closure return `None` for `ψˢ_out` when `2μ·A_down = 0` (equivalently, the `LinearOperator` term that this cell's outflow contributes to the next cell is the zero operator).

**Minimum geometric data structure to subsume the three branches.** Two clean options:

**Option A** — extend `StreamingTerms` so the slab factory populates the curvature fields with neutral values:
```
slab → face_area_inner = 1.0,
       face_area_outer = 1.0,
       delta_A_over_w  = 0.0,
       alpha_in        = 0.0,
       alpha_out       = 0.0,
       tau_mm          = 1.0,   # so c_out = α_out/τ = 0, c_in = α_in = 0
       volume          = chord  (already true)
```
With this, `cell_balance_terms` works on slab too — no special-case helper. The branch on `alpha_in is None` disappears. The cylindrical-degenerate branch becomes the same formula with `abs_mu < 1e-15`.

**Option B** — keep the three explicit branches but at a HIGHER layer (the `LinearOperator` matvec for the cell that this strategy participates in), not inside `DiamondDifference.update`. The strategy stays geometry-agnostic; geometry data flows through `visit.streaming_terms` only.

**Option A is the cleaner answer.** It is consistent with the project rule "geometry as data, not as control-flow" (CLAUDE.md Cardinal Rule 2). The current `alpha_in is None` test is exactly the "stringly-typed dispatch" anti-pattern the `coding-elegance` skill flags — geometry is hiding inside the type-system in a way that forces a runtime branch.

### Q2. What does the closure look like across the three?

The spatial closure `ψˢ_out = 2·ψ_avg − ψˢ_in` is **the same formula in slab and non-degenerate curvilinear**. (In the slab branch today it is hidden — the recurrence's `ψ_out = a·ψˢ_in + 2q/denom` is the WDD closure with `denom`-elimination of `ψ_avg` for bit-equality with cumprod. In exact arithmetic it equals `2·ψ_avg − ψˢ_in` after the average step.)

The angular closure `ψᵃ_out = (ψ_avg − (1−τ)·ψᵃ_in)/τ` (Morel-Montry) is **identical in non-degenerate curvilinear AND cylindrical-degenerate**, and **absent from slab** (no angular redistribution).

**Can both be expressed as one generalised "downstream = 2·avg − upstream"?**

Spatial closure: yes literally. `downstream_spatial = 2·avg − upstream_spatial`. Slab fits when you accept the bit-identity caveat above.

Angular closure: not literally — the Morel-Montry weight `τ` ≠ ½ in general, so the closure is `downstream_angular = (avg − (1−τ)·upstream_angular)/τ`. This reduces to the symmetric `2·avg − upstream` only when `τ = ½` (which is the Diamond Difference angular limit, away from the M-M clamp). The general statement is **"downstream is a τ-weighted antilinear combination of avg and upstream"** — and the spatial closure is the `τ = ½` instance.

A clean unified Protocol-level pattern is therefore:

```
closure(direction)  →  weight τ_direction ∈ (0, 1]
ψ_out_direction     = (ψ_avg − (1 − τ_direction) · ψ_in_direction) / τ_direction
```

with `τ_spatial = ½` (always, for DD) and `τ_angular = st.tau_mm` (the M-M clamp). Slab has only the spatial closure; curvilinear has both; cylindrical-degenerate has only the angular one (the spatial direction has zero weight on this ordinate — `A_down = 0` ⇒ closure produces `None`).

**Both spatial and angular closures are instances of one Protocol** — the "WDD-style face-balance closure": given an average and an upstream, produce a downstream via a τ-weighted antilinear combination. The TWO directions (spatial, angular) are just two consumers of the same algebra, fed by two upstream streams.

### Q3. What changes to `CellVisit` would be needed?

Today `CellVisit` carries:
- `cell_idx: int`
- `streaming_terms: StreamingTerms` (with optional curvature fields)
- `face_area_downstream: float | None` (sweep-resolved)

For DD's body to never inspect `visit.streaming_terms.alpha_in` or `visit.streaming_terms.abs_mu` for branching, the strategy needs to read **uniformly-populated geometric data**. The minimum changes:

1. **Slab `StreamingTerms` factory populates neutral curvature values.** As in Q1 Option A. `alpha_in/out = 0`, `delta_A_over_w = 0`, `tau_mm = 1`, face areas = 1, `volume = chord` (already true). With this, all `assert st.field is not None` in `cell_balance_terms` succeed for slab too.

2. **`face_area_downstream` is always populated.** Slab gets `face_area_downstream = 1.0`; cylindrical-degenerate gets `face_area_downstream = 0.0` (NOT `None`). This is the data that distinguishes "I have a spatial face that wants closing" (slab, non-degenerate curvilinear) from "I don't" (cyl-degenerate). The strategy reads ONE number; geometry chose the value.

3. **(Optional) An explicit `has_angular_redistribution: bool` flag on `StreamingTerms` or `CellVisit`.** Or equivalently: the angular closure's `τ` carries the discrimination — `τ = 1` ⇒ angular closure is the identity (`ψᵃ_out = ψ_avg`), and the caller knows not to write it back if the geometry doesn't track angular state. Today the flag is `upstream_state.angular_upstream is None`. Either is fine; making it data on `StreamingTerms` (e.g. `requires_angular_state: bool`) follows the same "geometry as data" principle.

4. **A `closures: tuple[Closure, ...]` field** (optional, deeper refactor) — `Closure` being a small dataclass `(direction: str, tau: float, has_output: bool)`. Slab carries `(closures = (Closure("spatial", 0.5, True),))`; curvilinear carries both; cyl-degenerate carries only the angular one (`has_output=True` because angular state propagates) AND a `Closure("spatial", 0.5, False)` with `has_output=False` (spatial closure exists in form but produces no output because `A_down = 0`).

### Q4. Is there a clean separation between cell-balance algebra and closure?

**YES — and this is the load-bearing observation for the four-operator algebra.**

The cell-balance algebra (denominator, numerator, solve for `ψ_avg`) is **one formula** across slab / curvilinear / cyl-degenerate once geometry is uniform data. It is the **diagonal action** of (L + C) on `ψ_avg`:
```
(L + C)_cell · ψ_avg = q + numer_upstream
```
with `(L + C)_cell` being `denom` and `numer_upstream` being the contribution from neighbouring cells / angular states (the off-diagonal action L applied to the upstream values).

The closure (spatial `ψˢ_out = 2·ψ_avg − ψˢ_in`, angular `ψᵃ_out = (ψ_avg − (1−τ)·ψᵃ_in)/τ`) is **the propagation step** — it encodes how `ψ_avg` plus upstream face/half-angle values produce the downstream face/half-angle values that the NEXT cell consumes as its upstream. This IS what L's matvec encodes: given an angular flux field, compute the streaming residual face-by-face, which is exactly the antilinear combination above.

**So the four-operator decomposition naturally lifts**:
- L (streaming + closure propagation) → the closure formula. ONE formula, parameterised by `τ` per direction (spatial / angular).
- C (collision) → the `Σt · V · ψ_avg` term in `denom`.
- S (scattering) and F (fission) → not visible at the cell-update layer; they enter via the cell `source` term `q`.

The current `_update_curvilinear` and `_update_cylindrical_degenerate` differ ONLY in the closure side (cyl-deg drops the spatial closure output). The `_update_slab` differs from both in operation order (bit-identity contract) AND in the absence of angular closure. After Q1 Option A unification, **the cell-balance algebra is one function call and the closure is a separate function call (or two — one per direction)**. That is exactly the (L + C) / closure separation the operator algebra wants.

---

## 3. Sketch of unified `DiamondDifference.update` body (≤ 30 lines)

```python
def update(self, visit, total_xs, source, upstream_state):
    st = visit.streaming_terms

    # ── Cell-balance solve: ONE formula, all geometries ─────────
    terms = cell_balance_terms_unified(
        st, visit.face_area_downstream, total_xs, upstream_state,
    )
    psi_avg = (source + terms.numer_upstream) / terms.denom

    # ── Spatial closure (WDD): outputs None if no downstream face
    psi_spat_out = None
    if visit.face_area_downstream > 0.0:
        psi_spat_out = 2.0 * psi_avg - upstream_state.spatial_upstream

    # ── Angular closure (Morel-Montry): outputs None if no
    # angular-state tracking on this geometry.  Detected by
    # `tau_mm == 1.0` (slab neutral) OR by upstream_state.angular_upstream
    # being None.  Either signal works; pick one as canonical.
    psi_angle_out = None
    if upstream_state.angular_upstream is not None:
        tau = st.tau_mm   # populated for all geometries; slab = 1.0
        psi_angle_out = (
            psi_avg - (1.0 - tau) * upstream_state.angular_upstream
        ) / tau

    return CellResult(
        cell_average_flux=psi_avg,
        outgoing_spatial_flux=psi_spat_out,
        outgoing_angular_state=psi_angle_out,
    )
```

No `if alpha_in is None` branch. No `if abs_mu < threshold` branch. The two `if`s remaining (`face_area_downstream > 0`, `angular_upstream is not None`) test "is this direction physically present on this geometry" — they are not geometry dispatch, they are "does this output exist".

---

## 4. Sketch of `cell_balance_terms_unified` (≤ 20 lines)

```python
def cell_balance_terms_unified(st, A_down, total_xs, upstream):
    # All curvature fields populated for slab too (neutral values).
    # Slab factory MUST set: A_in=A_out=1, dA_w=0, alpha_in/out=0,
    # tau_mm=1, V=chord.

    A_total = st.face_area_inner + st.face_area_outer
    c_out   = st.alpha_out / st.tau_mm                    # 0 for slab
    c_in    = (1 - st.tau_mm) / st.tau_mm * st.alpha_out  \
              + st.alpha_in                               # 0 for slab

    denom = (
        2.0 * st.abs_mu * A_down                     # slab: 2μ·1
        + st.delta_A_over_w * c_out                  # 0 for slab
        + total_xs * st.volume                       # slab: Σt·chord
    )

    psi_ang = upstream.angular_upstream
    ang_contrib = 0.0 if psi_ang is None \
                  else st.delta_A_over_w * c_in * psi_ang
    numer_upstream = (
        st.abs_mu * A_total * upstream.spatial_upstream  # slab: 2μ·ψˢ_in
        + ang_contrib                                    # 0 for slab
    )

    return CellBalanceTerms(denom=denom,
                            numer_upstream=numer_upstream,
                            c_in=c_in, c_out=c_out)
```

Single function, ≤ 20 lines, subsumes today's `cell_balance_terms` AND `cell_balance_terms_degenerate` AND the implicit slab algebra.

---

## 5. Pitfalls — known cases where naive unification fails

1. **Slab bit-identity regression.** The hand-calc tests at `tests/sn/spatial/test_diamond.py` assert `np.array_equal` against the cumprod operation order `ψ_out = a·ψˢ_in + 2q/denom; ψ_avg = ½(ψˢ_in + ψ_out)`. The unified body computes `ψ_avg = (q + 2μ·ψˢ_in)/denom; ψ_out = 2·ψ_avg − ψˢ_in`. Same answer in exact arithmetic, different IEEE-754 ULP. Need to either (a) accept and re-baseline tests to `np.allclose(rtol=1e-13)`, or (b) keep the slab op-order as a special case INSIDE `cell_balance_terms_unified` (defeats the purpose). The architecturally clean choice is (a) — the bit-identity contract was for the migration phase, not a forever rule.

2. **Cylindrical-degenerate axis: `2μ·A_down = 0` propagates a zero spatial output naturally**, but the caller must check `psi_spat_out is None` (today: check `outgoing_spatial_flux is None`). The unification preserves the existing API; no caller change needed if the closure step writes `None` when `face_area_downstream == 0.0`. Suggestion: define the threshold once (replace the magic `1e-15` literal); `face_area_downstream == 0.0` is the geometric truth, and the threshold is the floating-point realisation.

3. **`τ_mm = 1` for slab is a synthetic choice.** Real M-M `τ` is in `(½, 1]` after the clamp. Slab geometry has no M-M weight at all — it has no half-angles. Setting `τ_mm = 1` is "the closure is the identity" → safe (it would yield `ψᵃ_out = ψ_avg / 1 = ψ_avg` IF the angular closure were ever evaluated, which it isn't because `upstream.angular_upstream is None` for slab). But it COULD bite a future audit reader who sees `tau_mm = 1.0` on slab and tries to interpret it physically. Document it: a slab `StreamingTerms` carries `tau_mm = 1.0` as a NEUTRAL element of the M-M closure, NOT a physical M-M weight.

4. **Near-zero `μ` at curvilinear poles (not the degenerate axis).** In curvilinear-non-degenerate, `2μ·A_down` can be small but non-zero near the pole; cancellation against `Σt·V` is normal. The unified formula handles this — no special case. The genuine degenerate-axis case (`μ < 1e-15`) is geometric, not a numerical-cancellation artefact: the level's axial cosine `μ_z → 1` zeroes the radial cosine. The `face_area_downstream == 0` (geometric) test is exact; the `abs_mu < 1e-15` (numerical) test is the floating-point realisation of the same fact.

5. **The `face_area_downstream is None` ↔ `0.0` API change.** Today: `None` signals "no spatial face flow". Proposed: `0.0`. Callers in `diamond.py` line 368-371 assert `face_area_downstream is not None` for the non-degenerate curvilinear branch — those asserts need to become `> 0.0`. Single grep + edit; not load-bearing.

6. **The angular closure's `τ = 1.0` synthetic value WOULD divide by 1 cleanly**, but if a future factory sets `τ < 1e-15` (which the M-M clamp at `reduced_operator.py:593` prevents — `tau = max(0.5, min(1.0, tau_raw))`), the closure divides by near-zero. Keep the M-M clamp; the slab synthetic `τ = 1.0` is well within it.

---

## Recommendation summary

The unification is **architecturally clean and supports the four-operator algebra**. The investigation finds NO genuine algebraic obstruction — all three branches collapse into one formula. The two practical caveats are:

- **Slab bit-identity tests** need re-baselining from `np.array_equal` to `np.allclose(rtol=1e-13)`. This is the right trade: the bit-identity contract was migration-phase scaffolding (per the diamond.py module docstring), and Phase G is the migration's natural endpoint.
- **`StreamingTerms` slab factory needs neutral curvature values populated** (Option A in Q1). This is a 6-line edit at `reduced_operator.py:416-421` adding `face_area_inner=1.0, face_area_outer=1.0, delta_A_over_w=0.0, alpha_in=0.0, alpha_out=0.0, tau_mm=1.0`.

After those two edits, `DiamondDifference.update` is the ≤ 30-line body above — geometry-polymorphic by data, not by control flow. The (L + C) / closure separation falls out for free, which is exactly what the four-operator algebra needs.
