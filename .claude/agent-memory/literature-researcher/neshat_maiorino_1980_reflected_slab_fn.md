---
name: Neshat-Maiorino 1980 — F_N method for reflected slab
description: Reflected-slab F_N closure (3 boundary integrals → (3N+2) linear system). Trivial extension of bare-slab F_N. Unlocks ≈12 Sood reflected-slab cases. Companion to sood_fn_method_full_extraction.md.
type: reference
---

# Neshat & Maiorino (1980) Ann. Nucl. Energy 7, 79-81 — extraction memo

PDF: `/workspaces/ORPHEUS/scratch/literature/Neshat-Maiorino(1980)The FN method for solving the critial problem for a slab with a finite reflector.pdf`

3-page Pergamon "Technical Note". Self-contained. Builds directly on Grandjean-Siewert (1979) NSE 69, 161 (bare-slab F_N). Authors: K. Neshat (NCSU) + J.R. Maiorino (visiting from IEA São Paulo).

## 1. Problem geometry & notation (verbatim)

Symmetric slab, 1G isotropic, two regions:
- **Core**: `-τ ≤ x ≤ τ`, `c₁ > 1`.
- **Reflector**: `τ < |x| ≤ b`, `c₂ < 1`. Reflector half-thickness `Δ = b − τ` is GIVEN; critical core half-thickness `τ` is the UNKNOWN.

BCs (Eq. 2): `ψ₁(τ,µ) = ψ₂(τ,µ)` interface continuity; `ψ₂(b,−µ) = 0`, `ψ₂(−b,µ) = 0` (vacuum at outer surfaces). Symmetry exploited (`ψ₁(−τ,µ)` ↔ `ψ₁(τ,−µ)`).

Eigenvalues: `ν₀` for core (`c₁`), `η₀` for reflector (`c₂`). `P₁ = ν₀ ∪ (0,1)`, `P₂ = η₀ ∪ (0,1)`.

## 2. Method specification — verbatim core machinery

**F_N polynomial ansatz** (Eqs. 9a-c — three unknown angular fluxes on the two boundaries):
```
ψ₁(τ, +µ)  = Σ_α a_α µ^α       µ > 0    (core, outgoing toward reflector)
ψ₁(τ, −µ)  = Σ_α b_α µ^α       µ > 0    (core, incoming from reflector)
ψ₂(b, +µ)  = ψ₂(−b, −µ)
            = Σ_α e_α µ^α      µ > 0    (reflector, outgoing through free surface)
```
Three coefficient arrays of length (N+1) each → **3(N+1) unknowns**.

**The three projection equations** (Eqs. 10, 11, 12), obtained by integrating the eigenfunction-weighted angular flux against full-range Case orthogonality (Eqs. 5a, 5b applied at x=±τ and x=±b):

```
Σ a_α A_α(ξ̂) − Σ b_α B_α^(2)(ξ̂) =  e^(−Δ/ξ̂) Σ e_α A_α(ξ̂)         ξ̂ ∈ P₂   (10)

Σ b_α A_α(ξ̂) − Σ a_α B_α^(2)(ξ̂) = −e^(+Δ/ξ̂) Σ e_α B_α^(2)(ξ̂)     ξ̂ ∈ P₂   (11)

Σ b_α A_α(ξ)  − Σ a_α B_α^(1)(ξ)
   = e^(−2τ/ξ) [ Σ a_α A_α(ξ) − Σ b_α B_α^(1)(ξ) ]                  ξ ∈ (0,1) (12)
```

**Critical condition** (Eq. 15) — root-find on τ:
```
e^(−2τ/ν₀) Σ [b_α B_α^(1)(ν₀) − a_α A_α(ν₀)]
        =     Σ [a_α B_α^(1)(ν₀) − b_α A_α(ν₀)]
```

**X-function building blocks** (Eqs. 13-14 — IDENTICAL recursive moment integrals to bare-slab Grandjean-Siewert):
```
A_0(ξ)       = 1 − ξ log(1 + 1/ξ)
A_α(ξ)       = −ξ A_{α−1}(ξ) + 1/(α+1),                 α ≥ 1

B_0^(i)(ξ)   = 2/c_i − 1 − ξ log(1 + 1/ξ),              i = 1,2
B_α^(i)(ξ)   = ξ B_{α−1}^(i)(ξ) − 1/(α+1),              α ≥ 1, i = 1,2
```

That is the WHOLE method. Compare bare-slab F_N: same `A_α`, same `B_α^(i)` (parametrized by `c_i` per region), no new transcendentals.

## 3. Solution algorithm

1. Normalize `a_0 = 1` (overall scale).
2. Solve `F_0` initial guess (Eqs. 16-17 — closed form in `ν₀`, `η₀`, `Δ`) → `τ⁽⁰⁾`.
3. **Each iteration**: pick collocation points
   - `ξ̂ ∈ P₂`: `ξ̂_0 = η₀`, `ξ̂_1 = 0`, `ξ̂_2 = 1`, remaining (N−2) equally spaced in [0,1]
   - `ξ ∈ (0,1)`: `ξ_0 = 0`, `ξ_1 = 1`, remaining equally spaced in [0,1]
4. With current `τ`, build (3N+2)-equation linear system from (10)+(11)+(12) collocated on those points.
5. Solve for `{a_α, b_α, e_α}` (linear with `a_0=1` fixed).
6. Substitute into Eq. 15 → new `τ`.
7. Repeat. Authors report **3-4 iterations** to convergence for all 8 test cases on IBM 370/165, ≤3 s per case for N ≤ 8.

## 4. Convergence (Table 2)

`F_3` already gives 3-4 sig figs for all 8 cases. `F_5` matches Burkart "exact" to 4 digits in 7/8 cases; `F_7` matches to all reported digits (4-5 sig figs) for ALL 8 cases. **N=5-7 sufficient for engineering accuracy; N=8 saturates printed precision.**

Cases studied (Table 1 — `c₁`, `c₂`, `Δ` in mfp):
```
1: c₁=1.01, c₂=0.09, Δ=0.5    → τ_c = 8.3107
2: c₁=1.01, c₂=0.90, Δ=1.0    → τ_c = 7.6778
3: c₁=1.30, c₂=0.09, Δ=0.5    → τ_c = 0.9246
4: c₁=1.30, c₂=0.90, Δ=1.0    → τ_c = 0.6027
5: c₁=1.50, c₂=0.09, Δ=0.5    → τ_c = 0.5943
6: c₁=1.50, c₂=0.90, Δ=1.0    → τ_c = 0.3597
7: c₁=1.91, c₂=0.09, Δ=0.5    → τ_c = 0.3346
8: c₁=1.91, c₂=0.90, Δ=1.0    → τ_c = 0.1893
```

## 5. Mapping to Sood/LA-13511 reflected-slab cases

The `*-N(Δ)-1-0-SL` family in Sood. Composition pairs (`c₁`, `c₂`) inferred from Sood Tables 2-7 cross-section data:

- **`Ua-H2O(N)-1-0-SL`** family: U-235 core (`c₁=1.30`) + H₂O reflector (`c₂≈0.90`) → matches NM Case 4 structure. Multiple Δ values (1, 5, ∞ mfp).
- **`UD2O-H2O(N)-1-0-SL`** family: UD₂O core (`c₁≈1.02`) + H₂O reflector → matches NM Cases 1-2 structure.
- **`Pua-H2O(N)-1-0-SL`** family: Pu-239a core (`c₁=1.50`) + H₂O reflector → matches NM Case 6 structure.
- Various `*-Fe-Na(N)-1-0-SL` and other reflector compositions (`c₂` derived from XS table).

Estimated Sood cases unlocked: **≈10-14 reflected-slab problems** (the entire `*-X(N)-1-0-SL` row of Table 2's catalogue). Won't unlock the **infinite-reflector** ones (`*-IN-1-0-SL`); those are Siewert-Burkart 1975 NSE 58, 253 (Ref. 1 here — uses Chandrasekhar invariance, not F_N, and is structurally distinct). For finite-Δ vacuum-outside cases NM-1980 is THE truth-set source.

## 6. Implementation feasibility — TRIVIAL extension

Per ORPHEUS `fn_method/core/` machinery:

| Required primitive       | Already in `core/`?                                        |
| ------------------------ | ---------------------------------------------------------- |
| `A_α(ξ)` recursion       | YES — `core/moments.py` (eq. 14a-b, identical)             |
| `B_α^(i)(ξ)` recursion   | YES — same file, just parametrize `c_i` per region         |
| Case eigenvalue `ν₀(c)`  | YES — `core/dispersion.py` (call once for `c₁`, once `c₂`) |
| F_N matrix assembly      | YES — `core/fn_matrix.py` (just larger block structure)    |
| Collocation-point picker | LIKELY — extend to 2 separate point sets `P₂`, `(0,1)`     |

**Net incremental effort**: assemble a (3N+2) block matrix from 3 sub-blocks of NM Eqs. 10/11/12. No new special functions. No new transcendentals. No two-region X-function — this is the great architectural win: the X-function is **medium-property local**, so each region uses its own `c_i` independently and the coupling lives entirely in the boundary projection equations.

**Verdict: TRIVIAL extension of `fn_method/slab/one_group.py`.** Estimate <1 day for a competent implementer, including L1 cross-check vs. NM Table 2 + a couple of Sood cases.

**No SymPy walls beyond bare-slab.** The exponentials `exp(±Δ/ξ̂)`, `exp(−2τ/ξ̂)` are clean. The only mpmath fallback already present: transcendental `ν₀(c)` root-find (same as bare slab). System dimension 3(N+1) ≈ 24 at N=7 → numpy/scipy solve handles trivially; no conditioning concern at the Δ values Sood actually tabulates (Δ ≤ 5 mfp typically; thicker reflectors → infinite-reflector limit, which uses different machinery anyway).

**Numerical pathology**: only one — at very large Δ, `exp(+Δ/ξ̂)` in Eq. 11 grows. NM use Δ ≤ 1 mfp in their cases without trouble; for Sood's Δ=5 mfp cases use `c₂≈0.9` → `1/η₀ ≈ 0.5` → `exp(2.5) ≈ 12`, still fine. If Δ→∞ comparison with Siewert-Burkart 1975 is desired, that's the dispatch on a separate paper.

## 7. Architectural seam — recommendation

**Place at `orpheus/derivations/fn_method/slab/reflected.py`** as sibling to `one_group.py`. Reasons:

1. Shares the entire `core/` machinery (moments, dispersion, X-function, F_N matrix builder). 90% code reuse.
2. The two-region structure is a slab-specific concern — the boundary topology (interface at τ, vacuum at b) is what changes, not the angular machinery.
3. A future `fn_method/sphere/reflected.py` (Sjostrand 1986, Sahni-Sjostrand) will mirror this exactly: same `core/`, different boundary projections. The pattern generalizes by geometry, not by "multi-region" abstraction.
4. A `multi_region/` package would prematurely invent an abstraction layer over only 2 examples (slab + sphere reflected); YAGNI applies.

**Public API suggestion**: `solve_critical_reflected_slab(c1, c2, Delta, N=7) -> tau`, mirroring `solve_critical_slab(c, N=7)`. Internally builds a `ReflectedSlabFNMatrix(c1, c2, Delta, N)` from `BareSlabFNMatrix`-compatible primitives.

## 8. Errata check (mental SymPy re-derivation)

- Eqs. 13a-14b (A_α, B_α^(i) recursions) — **clean** ; matches Grandjean-Siewert 1979 verbatim. Already used in our bare-slab solver.
- Eq. 17 (F₀ initial guess `b_0`) — sign and exponential structure consistent (DR-checked: `(1−e^(−2Δ/η₀))/[B_0^(2)(η₀)² − A_0(η₀)² e^(−2Δ/η₀)]` is the standard 2-stream reflector albedo form).
- Eq. 16 (`τ⁽⁰⁾` from `F₀`) — checks out as `−(ν₀/2) log[…]` with the asymmetry between numerator & denominator producing positive `τ` for `c₁>1` (numerator < denominator since `b_0 > 0` and the `−B_0^(1)` term pulls down). Sign convention consistent.
- **No errata flagged**. The paper is short, mature (Pergamon 1980, peer-reviewed), and the only structural risk is a typo in the exponential signs of Eqs. 10/11. Cross-check against the Eq. 15 critical condition: the SAME exponential pattern `exp(−2τ/ξ)` appears, with the same sign convention (decay from the inside, growth from the outside). Self-consistent.

## 9. Supplementary references possibly needed

For the implementation itself: **NONE**. The paper is fully self-contained ON TOP of Grandjean-Siewert 1979 (bare-slab F_N — already in ORPHEUS `core/`) and Case-Zweifel 1967 Ch. 2-3 (eigenfunctions — already known).

If Sood truth-value cross-check requires more than just NM Table 2:
- **Burkart 1975 PhD thesis / Burkart 1976 Trans. ANS 24, 190** (NM Refs. 9-10) — gives the "Exact" column in NM Table 2. NOT needed for the F_N implementation; NM's own F_7 column matches Exact to all printed digits.
- **Siewert & Burkart 1975 NSE 58, 253** (NM Ref. 1) — INFINITE reflector case. Required ONLY if Sood `*-IN-1-0-SL` (infinite-reflector) cases are targeted; uses different (Chandrasekhar invariance) machinery and is OUT OF SCOPE for an F_N extension.
- **Neshat-Siewert-Ishiguro 1977 NSE 62, 330** (NM Ref. 4) — earlier P_N + invariance solution to the same problem. Useful for cross-validation only.

## Action items for user

1. Drop **NOTHING ELSE** for the basic reflected-slab F_N implementation. NM-1980 + the existing `fn_method/core/` is sufficient.
2. If/when reflected-sphere F_N is targeted next: acquire **Sjostrand 1986** (and possibly Sahni-Sjostrand) — the structural analog of NM-1980 but for sphere geometry. Already noted in `sood_fn_method_full_extraction.md` §7.
3. If/when `*-IN-1-0-SL` (infinite-reflector slab) cases are targeted: acquire **Siewert-Burkart 1975 NSE 58, 253**. Different method (invariance principles) — separate solver path.

## Cross-references

- `sood_fn_method_full_extraction.md` — parent context; problem catalogue + 5-case complexity ramp.
- `kaper_lindeman_leaf_1974_fn_method.md` — bare-slab F_N acquisition wall (still open).
- This memo extends the `fn_method/slab/` story from BARE-only to BARE+REFLECTED.
