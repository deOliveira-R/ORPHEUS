---
name: cs4a-kernel-binding-representation-frames
description: Kernel×Space→Operator IS a *-representation of the SO(3)-commutant algebra ⊕_ℓ M_ng; the four laws and which hypothesis each needs (tightness governs MULTIPLICATIVITY, never the adjoint). CS4a-R attack, 2026-08-21.
metadata:
  type: project
---

# "Kernel × Space → bound operator" is a **\*-representation**, and the laws sort by hypothesis

Attack surface: ORPHEUS campaign 1 CS4a (`orpheus/transport/kernels.py`,
`orpheus/numerics/{space,axis,frame,operator}.py`, `orpheus/homogeneous/solver.py`),
reviewed 2026-08-21 on branch `feature/cs1-energy-space`.

**Why:** this is the durable math, reusable on every future binding question
(CS2's angular frame mint, CS4c's dispatch collapse, Campaign 2's pencil).
**How to apply:** fire it whenever a design says "representation-free data +
a constructor that binds it to a space".

## 1. The exact identification (a theorem, not an analogy)

SO(3) acts on `L²(S²) ⊗ ℝ^ng` (trivially on the energy factor). `L²(S²) = ⊕_ℓ H_ℓ`
with each `H_ℓ` irreducible and multiplicity 1, so on the tensor product `H_ℓ` has
multiplicity `ng`. By Schur the **commutant** is

    𝒜_L  =  ⊕_{ℓ≤L}  M_ng(ℝ)

— which is `ScatteringKernel.moments: tuple[(ng,ng) array, ...]` **verbatim**. So
"representation-free" has an exact meaning: *an element of the commutant algebra,
presented by its isotypic blocks*. The ℓ-index indexes **eigenSPACES (irreps)**,
which are canonical; no basis of `H_ℓ` is chosen. The claim is honest — in ANGLE.

Corollaries that drop out for free:
- 𝒜_L is an **algebra**: `(K∘K')_ℓ = M_ℓM'_ℓ`, `+` componentwise, `†` = per-ℓ
  transpose (energy metric is counting **by the EnergyAxis theorem**).
- `truncated(L')` is the quotient `𝒜_L ↠ 𝒜_{L'}` — a surjective \*-homomorphism
  (projection onto a direct summand), not "the one representation-free morphism"
  by fiat.
- The binding `bind_V = FrameBase.conjugate` (`R∘A∘M`) is a candidate
  **\*-representation** `ρ_V : 𝒜_L → End(V)`. That is the categorical shape.

## 2. ⭐ The four laws and WHICH hypothesis each needs (the reusable table)

Faces: analysis `M : V→M̂`, reconstruction `R : M̂→V`, metrics `G_V`, `G_M̂`.
`bind(K) = R K M`.

| law | necessary AND sufficient condition | name |
|---|---|---|
| `bind(aK₁+bK₂) = a·bind(K₁)+b·bind(K₂)` | — | free (linearity) |
| **`bind(K)† = bind(K†)`** | **`M† ∝ R`** | **GALERKIN**, *not* tightness |
| **`bind(K₁∘K₂) = bind(K₁)∘bind(K₂)`** | **`M R = I_{M̂}`** | **TIGHTNESS / Parseval** |
| `bind(1) = I_V` | `R M = I_V` — impossible when `dim M̂ < dim V` | it is a PROJECTOR |

Proof of row 2 (do not re-derive): need `R† X M† = R X M ∀X`; put `X=|a⟩⟨b|`,
then `a=b` ⟹ `R†a ∥ Ra ∀a` ⟹ `M† = λR`; the λ then cancels on both sides.

⚠ **The trap this table exists to catch:** a charter that writes
`bind(K)† = bind(K†) ⟺ the frame is tight`. Both halves are wrong under the
GALERKIN convention (`M = R†`, which is what `FrameBase.analysis` is): ⟸ is
vacuous (the LHS holds unconditionally), ⟹ is false (a non-tight Galerkin frame
satisfies the LHS). The statement becomes TRUE only if the binding's analysis is
the **canonical dual** `project = G⁻¹∘analysis`: then `M† = R G⁻¹ ∝ R ⟺ G ∝ I ⟺
tight`. ⟹ **an equivalence about a binding is well-posed only after the binding
declares WHICH of the two analysis verbs it uses.** ORPHEUS ships both
(`FrameBase.analysis` vs `FrameBase.project`, `frame.py:195/310`) and
`FrameBase.gram` (`:254`) already computes the cross Gram `MR` — so the tightness
test is one line against machinery that exists.

## 3. The Σw normalisation is a VALUE contract, not a law hypothesis

`E = (1/W)·R†` (the `from_isotropic` 1/W). Then `E† = (1/W)R` and the two 1/W's
**cancel exactly** in `bind(K)† = bind(K†)`. So "Σw = 4π owned once by the
binding" protects the operator's VALUE and has nothing to do with the adjoint law.
The F4-addendum's R/E claims are otherwise all CORRECT: `R† = ` constant
embedding, `Rᵀ = ` w-weighted embedding, differing by G-conjugation. (Verified by
hand; `⟨Rψ,c⟩ = Σw_nψ_n c` ⟹ `(R†c)_n = c`.)

## 4. The ℓ=0 blindness (why the witness cannot adjudicate)

At `L = 0`, `MR = Σ_n w_n/W = 1` for **any** rule and `M† ∝ R` trivially, so both
laws hold for every quadrature — the ℓ=0 binding (fission, isotropic scattering,
the whole CS4a surface) cannot separate Galerkin from tight. Any binding gate
whose only witness is ℓ=0 lands **designed-green** (plan-authoring §6c). The
discriminator needs `L ≥ 1` and two rules, one 2L-exact and one not.

## 5. Frame DIMENSION is the space's, not the kernel's

`Σ_sℓ` is the eigenvalue on the **whole** degree-ℓ eigenspace, whose dimension is
`2ℓ+1` in 3-D, `1` in slab (m=0 only), and the multiplicity of the trivial
`G_geom`-irrep in `H_ℓ` in curvilinear 1-D (the symmetry-quotient result). So a
frame mint sizing itself from `kernel.order` alone is **slab-correct and
multi-D-wrong**. The tree already knows this at `material_xs_field.py:993`
(`n_m = 2*l+1`); the KERNEL does not, and must not.

## 6. Angle-free but energy-BOUND — the asymmetry the docstring does not state

`moments[ℓ]` is an `(ng,ng)` matrix **in a specific group structure**. The kernel
carries `ng` but no `EnergyAxis`, and `assert_energy_extent_conforms`
(`_energy_conformity.py:55`) compares `axis.shape` — a **cardinality** — where the
campaign's own chartered identity doctrine says **content** decides (an
`EnergyAxis` from grid A ≠ one from grid B at equal ng, `axis.py:360`). ⟹ a 2g
kernel built on grid A binds without complaint to a space carrying grid B.
Fix: kernels carry their `EnergyAxis` (same hoisted `from_materials` rule) and the
guard compares axis IDENTITY.

## 7. The isotropic bindings are `Π₀ ⊗ K`, not `K ⊗ I` ⟹ CS2 makes `_agreed_space` wrong

`bind(K)` on the SN composite is NOT a tensor lifting of an energy operator; it is
`E∘K∘R` on the angular leg — the ℓ=0 projector — tensored with `I` elsewhere. That
is a genuine **per-leg** binding, and `operator.py:_agreed_space` (`:338-398`)
states in its own docstring that agreement is the right law only while every bound
factor describes the WHOLE space, and names "an energy-dependent group kernel bound
on its own axis" as the arrival that breaks it. CS2's frame-at-binding IS that
arrival. The campaign plan cites `_agreed_space` only as a guard that flips from
skipped to active (§8 blast radius), never as a LAW with an expiry.

## Related
[[harmonic-frame-ownership-funk-hecke]] (WHO owns the frame — this file is WHAT
the binding must satisfy) · [[unified-frame-api-design]] (the three verbs) ·
[[projection-discipline-hierarchy-frames]] (`project` vs `analysis`) ·
[[xs-coarsening-collapse-marginalize-vs-average]] (the naturality-along-collapse
condition) · [[quadrature-symmetry-quotient-frames]] (the per-ℓ frame dimension in
curvilinear geometry) · [[fission-rank1-normal-form-dead-functional]] (why MPO
does not fire on F).
