---
name: sigma-y-orbit-shape
description: The derived SHAPE of the S^2/<sigma_y> orbit space (the shipped cylindrical fold) — closed 2-disk, empty syzygy ideal, det P = 4p3, a CIRCLE-valued singular stratum — and the measured realization option matrix. Read before adding the sigma_y entry to _ORBIT_CATALOGUE.
metadata:
  type: project
---

# `S² / ⟨σ_y⟩` — the derived shape, and the realization fork

**Derived 2026-08-31** by the Procesi–Schwarz procedure
(`.claude/plans/symmetry_quotient_plan.md` §I.1; Procesi & Schwarz, *Invent.
Math.* **81** (1985) 539–554). Full memo + verbatim SymPy output:
`scratch/sigma_y_orbit_derivation.md`; re-runnable probe (57 checks, 9/9
mutations red): `scratch/probe_sigma_y_orbit.py`.
⚠ `scratch/` is UNTRACKED — this note carries the result so it survives a
`git clean`.

**Why:** `Quadrature.folded_product` ends in
`full.quotient(SubgroupOfO3.Mirror("y"))` (`directional.py:845`), so this is the
orbit space every cylindrical SN solve already uses. `Mirror("y").name` is
`'sigma_y'`, so the catalogue key would be `(Sphere, "sigma_y")`.

## The answer

| procedure output | value |
|---|---|
| minimal invariants | `p₁ = x`, `p₂ = z`, `p₃ = y²` (degrees 1,1,2). Complete + free by Molien: `M(t) = 1/((1−t)²(1−t²))` **equals** the Hilbert series of `ℝ[x,z,y²]` |
| syzygy ideal `I` | **`(0)` — empty.** Predicted by Chevalley–Shephard–Todd (`σ_y` is a *reflection*, `det = −1`), confirmed by lex elimination and by Jacobian rank 3 |
| `P = ⟨∇pᵢ,∇pⱼ⟩` | `diag(1, 1, 4p₃)`; eigenvalues `{1 (×2), 4p₃}` — the whole content is the single inequality `p₃ ≥ 0` |
| `det P` | `4p₃`; on `S²`, `4(1 − p₁² − p₂²)` |
| `ℝ³/H` | the closed half-space `{p₃ ≥ 0}` |
| **orbit space** | **the closed unit 2-disk** `{(η, μ) : η² + μ² ≤ 1}`, with `ξ² = p₃ = 1−η²−μ²` eliminated by the sphere relation |
| `dim` | **2 → 2.** `H` is finite, so the dimension does NOT drop — the root of every asymmetry with `S²/SO(2)` |
| **singular stratum** | the **boundary CIRCLE** `η² + μ² = 1`, pre-image `S² ∩ {ξ = 0} = Fix(σ_y)` — the purely-meridional directions |

⛔ **`singular_stratum: tuple[float, ...]` CANNOT hold this.** `S²/SO(2)`'s
stratum is two points; this one is uncountable. The spelling that survives both
is the one the shipped gate already uses: **the stratum is the vanishing locus
of `det_gram`, i.e. the realization's topological boundary**
(`∂[-1,1] = {−1,1}`, `∂D² = S¹`).

## The fork the entry cannot decide alone

`Quotient.contains` delegates to `realization.contains` (`manifold.py:509`) and
`_ambient` reads the realization (`:656`). `[M]` against the shipped
`folded_product(4,8).measure.nodes` `(16,3)`, its `μ_y<0` twins, and the ERR-080
forgery `column_stack([leggauss(8)[0], 0, 0])`:

| realization | shipped | twins | ERR-080 |
|---|---|---|---|
| `SPHERE` | ADMITTED | ⛔ **ADMITTED** | refused |
| the 2-disk (chart codomain — the DOCUMENTED meaning, `manifold.py:476`) | refused (shape) | refused (shape) | refused (shape) — but ⛔ **ADMITTED once charted** |
| a closed hemisphere `μ_y ≥ 0` (a fundamental domain) | ADMITTED | **refused** | **refused** |

* `realization = SPHERE` is **topologically false** (`D² ≇ S²`; `χ = 1` vs `2`)
  and its `contains` is bit-for-bit `SPHERE.contains` on every input — zero
  information about the quotient.
* The disk is what the field's docstring *requires*, and it is **Mode-12 blind
  to ERR-080**: the chart drops `μ_y`, which is exactly what the forgery
  falsifies. It also needs the chart, which `manifolds.rst:1167` records as
  ⛔ *not a slot*.
* The hemisphere is the only option that refuses both wrong inputs. Costs: it
  widens `realization` from *chart codomain* to *fundamental domain*; it is
  canonical only because `σ_y` is a **reflection** (a closed `Cn` sector is not
  a strict fundamental domain — its two edges are identified); and the
  predicate **must be `≥`, not `>`** — see below.

⭐ **`DiscreteMeasure.quotient` (`measure.py:1016-1023`) has already chosen: it
does `nodes[representative]`, a selection of parent nodes with no chart.** Every
measure the tree derives that way is a section in the base's ambient.
`gauss_legendre_on_mu` never goes through it, which is why `S²/SO(2)` ships
chart coordinates and looks unambiguous — and why `[M]` the two support tags
disagree in kind: `'[-1,1]'` (the realization) vs `'S^2/sigma_y'` (the
quotient). That split is live and independent of this entry.

## Two measured facts a gate for this entry needs

1. **The shipped rule never populates the stratum.** `[M]` `μ_y ∈ [+0.1945,
   +0.8688]` on all 16 nodes — the even-`n_φ` staggering makes the fold *free*.
   So any predicate difference confined to `μ_y = 0` is invisible to the
   quadrature fixture.
2. **The closure's edge data sits exactly ON the stratum.** `[M]` nodes at
   `ω/π ∈ {.125,.375,.625,.875}` ⟹ cell edges at `{0,.25,.5,.75,1}`, so the
   outer edges are `ω = 0, π` where `ξ = 0`. There `α[0] = α[M] = 0`, and
   `mu_start_per_level = −√(1−μ_p²)` gives seed points `(η, 0, μ)` on `S²` to
   `0.0` with `μ_y == 0.0` exactly. ⟹ **a strict `μ_y > 0` predicate would
   refuse a direction production marches from.** Gate the stratum on the
   closure, not on the nodes.

⭐ This is the cylindrical analogue of the documented spherical statement that
the `α`-dome closes at `μ = ±1` = the `SO(2)` stratum: **interior carries the
nodes, boundary carries the degenerate closure data** — an orbifold,
discretised. ⚠ `[R]` only: whether the redistribution term **is** the quotient's
connection is unproved, the cylindrical twin of #429 Phase 1.3.

Behavioural lesson distilled from this: [[lessons]] **L-016**.
