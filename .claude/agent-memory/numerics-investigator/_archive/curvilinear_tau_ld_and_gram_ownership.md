# Archive — curvilinear tau x LD arity, and the adjoint/Gram ownership audit

Moved COLD from `lessons.md` 2026-09-21 (verbatim, lines 1015-1165).
Digest successor: L25 / L26. Records: `scratch/tau_under_ld_dip_analysis.md`, `scratch/adjoint_gram_ownership_audit.md`
(both UNTRACKED in the working tree — this archive copy is the durable one).

## L25: An ARITY / "does this contract widen?" question is a THEOREM question — ask what COMMUTES; and a positivity question is a DIFFERENT question the asymptotic expansion is structurally blind to

Curvilinear LD × Morel–Montry τ (2026-08-25/26). Record:
`scratch/tau_under_ld_dip_analysis.md`; probes `scratch/probe_tau_ld_0{1..6}_*.py`
(**106 SymPy checks, 0 failures**). Seven transferable points.

1. ⭐⭐ **"Does accessor `f(x)` have to become `f(x, y)`?" is answered by asking what
   the defining conditions COMMUTE with — not by an expansion.** A **scalar convex
   combination commutes with every linear map**, so if a quantity is defined by (i)
   membership of a convex set and (ii) exactness of a scalar blend, then *both*
   conditions are the SAME scalar statement in every component of every linear
   representation, and the widened form is an **overdetermined system whose every
   row returns the same value**. `[M]` solved independently at `n_mom = 4`: all
   rows give `τ = (μ_n−e_−)/(e_+−e_−)`. Hypotheses to CHECK, all three cheap: the
   scalar is independent of the widened index; the projection is linear; the set is
   convex. Cost: minutes, against a multi-day expansion.
2. ⭐⭐ **A SIGNATURE is not an invariance proof, and a spatial-argument-free
   producer is the classic tell.** "β cannot accept a mesh, therefore β is
   mesh-independent" is `vv-principles` Mode 8's SIGNATURE-tautological class: the
   signature is a consequence of the DERIVATION'S SCOPE (M-M 1984 and BMC 2010 both
   hold space *continuous* and say so), not evidence about the joint problem. The
   conclusion was true; the argument carried zero information. **Always ask what
   analysis produced the object, not what its parameters are.**
3. ⭐⭐ **The asymptotic expansion and the POSITIVE CONE answer different questions,
   and the expansion is structurally blind to the cone.** A sign-alternating
   cell-to-cell mode is EXCLUDED BY THE ANSATZ, so no order of the expansion sees
   it — which is why Palmer–Adams carry "limits to a *stable* diffusion equation"
   as a SEPARATE acceptance criterion. `[M]` transmission sign ladder, all derived
   in one probe: **DD flips at `τ_opt=2` (Padé(1,1)), bare LD at `τ_opt=3`
   (Padé(1,2)), lumped LD NEVER (Padé(0,2), discriminant `1−2<0`)** — and
   `a_LD·τ_opt → −2`, i.e. sign-flipped and decaying only as `1/τ_opt`. ⟹ when a
   scheme is "verified in the diffusion limit", ask which of the two it was.
4. ⭐ **Where a redistribution/coupling operator is a TENSOR PRODUCT `R ⊗ A(θ)`
   with disjoint index sets, EVERY functional of it factors** and the θ-scalar is
   untouched by any change to `R`. `[M]` verified at `n_mom = 1,2,3,4`; `[M]` the
   scalar's free symbols contained no spatial symbol at all — **grep the free
   symbols, it is a one-line proof.** Corollary that pays: the leading-order
   discrete diffusion equation then carries no θ, so a SPATIAL defect cannot be
   repaired by θ and vice versa — **two orthogonal failure modes**, and conflating
   them is the whole confusion.
5. ⭐⭐ **When one knob is exact-by-construction, price the OTHER knob's error in
   the first knob's units — the ratio is usually the finding.** `[M]` with τ at
   Morel–Montry (β ≡ 0), a starting-cosine error of only **1.6 % at S4, 0.05 % at
   S32** reproduces the ENTIRE contamination τ exists to remove, because
   `β_eff(μ_s) = S + (μ_s+1)·S_e` is EXACTLY affine and `|S_e/S_diamond|` GROWS
   with N (24× → 333×). ⟹ the celebrated knob was never the risk; its neighbour
   was. Look for the affine law before assuming two knobs are independent.
6. ⚠ **A published RULE OF THUMB carries the order it was derived at.** M-M's
   summary "the dip is eliminated with any spatial scheme as long as the starting
   flux is not seriously UNDERestimated" is `dβ/dμ_s = S_e > 0`, and `[M]` **`S_e`
   flips sign between N=2 (their own test case, `+9.1e-1`) and N≥4 (`−1.1e-1 …
   −1.7e-5`)** — the safe DIRECTION inverts. Honest companion: `|S_e|` falls 5
   orders over the same range, so the direction inverts while the stakes collapse.
   Same family as L21-6 (a low-order-consistency benefit decaying with angular
   order). **Evaluate the quantity, never the heuristic.**
7. ⚠ **"Lumping" names ≥3 different operations and they disagree.** `[M]`
   Legendre-diagonal lumping destroys `R₀₁` and breaks the per-moment-row
   flat-flux identity; **nodal ROW-SUM lumping (which is what Palmer–Adams's FL
   actually is) PRESERVES `R₀₁` exactly** — only `R₁₁` moves, `ΔA/3 → ΔA` — so the
   identity survives untouched. And row-sum lumping the GRADIENT gives transmission
   **identically 0** (a degenerate scheme), which is why MWS lump everything except
   the gradient. ⟹ a ⛔ banner saying "X may not be lumped" must name the BASIS and
   the MATRIX; mine (inherited) condemned the wrong operation. The real freedom:
   the infinite-medium identity pins only `rowsum(L)`, leaving one parameter per
   row, and a MONOTONE (`A⁻¹ ≥ 0`) member exists inside it at the cost of dropping
   to first order — **the accuracy/positivity trade is a choice of that parameter,
   not a property of "lumping".**

## L26: An OWNERSHIP question ("who owns this quantity?") is answered by TWO measurable discriminators — does it move the METRIC, and is `∂A/∂k` DIAGONAL — never by what the quantity is called

Written 2026-08-26, the adjoint/Gram ownership audit of the 1-D curvilinear SN
streaming path (memo `scratch/adjoint_gram_ownership_audit.md`, probes
`scratch/probe_adjoint_gram_0{1..6}_*.py`).

1. ⛔ **"Does X influence the adjoint?" is the WRONG question when the adjoint is
   a reverse-mode VJP** — `[M]` `apply_transpose ≡ apply.T` to `3.05e-16`, so the
   answer is "everything" and it carries no information. The informative question
   is the **symmetry character of the INCREMENT `∂A/∂k`**, measured by perturbing
   one entry and taking the dense difference. ⟹ **build the difference operator,
   not the operator.**
2. ⭐⭐ **Calibrate the symmetry ratio before reading it.**
   `‖A+Aᵀ‖_F/‖A‖_F = √(2 + 2·tr(A²)/‖A‖²_F)`: **0 = skew, 2 = symmetric,
   √2 = `tr(A²)=0` (no entry pairs with its transpose partner, i.e.
   triangular-like)**. `[M]` controls: skew `0.000000`, symmetric `2.000000`,
   strictly-lower-triangular `1.414214`. Without the calibration a √2 reading
   looks like "very asymmetric" when it actually means *triangular*, which is a
   structural statement about a MARCH, not a magnitude.
3. ⭐⭐ **Two derived constants from the SAME two numbers on ADJACENT lines can sit
   on opposite sides of the self-adjointness split.** `[M]` `c_out = α_out/τ`
   contributes an **exactly diagonal** block (`offdiag = 0.000e+00`, sphere AND
   cylinder — self-adjoint in EVERY metric); `c_in = (1−τ)/τ·α_out + α_in`
   contributes an increment with `tr = 0` and **0 of its 8 nonzeros pairing with
   their transpose**. ⟹ a design that bags them as "the closure constants" fuses
   a reaction-like scalar with a transport-like coupling. **Diagonal ⟹ order-free
   ⟹ may live anywhere; non-diagonal ⟹ welded to the traversal that reads it.**
4. ⭐⭐ **A `⟨Aψ,φ⟩_G = ⟨ψ,A†φ⟩_G` reciprocity gate CANNOT adjudicate the choice of
   `G`.** With `A† ≡ G⁻¹AᵀG` it is an identity for EVERY invertible `G`: `[M]`
   `1.4e-16` under Euclidean `G′=1`, random `G′~U(1,2)` and adversarial
   `G′=(V·w)³`; the mismatch control (A.H on `G`, pairing on `G′`) reads
   `8.22e+00`. So the gate + its wrong-metric control prove **consistency and
   loadedness, never choice** — Mode 12 with the whole invertible-diagonal group
   as the stabiliser. ⟹ **the adjoint does not pin the metric; the physical
   FUNCTIONALS do** (scalar flux `Σw_nψ_n`, reaction rate `∫Σφ dV`), which is
   exactly why the metric belongs to the SPACE and the operator is its consumer.
5. ⭐ **One symbol, four roles, four owners — and the separating experiment is a
   GLOBAL RESCALE through the production constructors, not instance surgery.**
   `[M]` rebuilding at `w → 3.7w` (cylinder/slab, the charts that admit it) scales
   `G` **exactly** by `3.7` and leaves `L` **bit-identical** (`2.06e-16` / `0.0`),
   because `α ∝ w` ⟹ the redistribution sees only the ratio `α/w`. The sphere
   **REFUSES** the same rescale (`τ_raw ∉ [0,1]`) because its angular-cell
   partition is the *cumulative weight* and needs `Σw = 2`. And a SINGLE-ordinate
   `w` perturbation is refused on any curvilinear mesh by the pole mirror's
   weight-preservation contract. ⟹ metric / scale-free ratio / absolute mesh width
   / admission precondition. **Ask which one a consumer means before moving it.**
6. ⭐⭐ **A "Gram" in an OPERATOR is usually the mass matrix under a DIFFERENT
   measure, and naming the measure settles the ownership.** `[M]` (SymPy, exact,
   sphere+cylinder) the curvilinear redistribution Gram is
   `R_kj = ∫ b_k b_j (∇·ê_r) dV` with `∇·ê_r = (d−1)/r`, so `R₀₀ = ΔA` is the
   **divergence theorem**, not a per-chart normalization (the "factor-of-2-absorbed
   α normalization" is an artefact of writing the measure as `r^{d−2}`). It is SPD
   ⟹ a genuine inner product, but `W = M⁻¹R` is **not** `λI` ⟹ not `M` rescaled;
   `W` is the Galerkin compression of the Radon–Nikodym derivative. ⟹ **measure
   from the CHART, basis from the SCHEME ⟹ the home is the (chart × scheme) pair**,
   and it is a *coefficient*, never a metric — `[M]` perturbing it leaves `G`
   bit-identical. ⚠ And when the two axes differ (`n_mom ≠ n_thread`, the
   ONETRAN/Hill column) the object is **rectangular** ⟹ a PAIRING, not a Gram: the
   word over-claims off the diagonal of its own design space.
7. ⛔ **Do NOT justify a transport metric as "the one that makes streaming
   skew-adjoint".** `[M]` `‖ĜL+(ĜL)ᵀ‖/‖ĜL‖ → 1.4172 ≈ √2`, flat over `nx = 4…64`,
   **same number on the slab**: a face-ELIMINATED marching operator is triangular,
   not skew. Skew-adjointness is a property of the (cell ⊕ face) saddle system;
   DD's `ψ_out = 2ψ̄ − ψ_in` substitutes the interior faces out and destroys it.
8. ⚠ **Two structural ZEROES that a naive inventory would call "not used".** `[M]`
   `μ_start` is inert on every shipped 1-D curvilinear fixture (its only consumer
   runs on NON-carrying levels, and every level carries: `{0}` sphere,
   `{0,1,2,3}` `folded_product(4,8)`) — control `4.04e-01` when a level is forced
   non-carrying. And `c_in[m=0]` is inert because it multiplies the ψ½ SEED, which
   the ray-decoupled `(A,A)` block feeds with **zeros** — control: the neighbouring
   index `m=1` moves at `3.08e-04`. ⟹ **each zero needs its OWN control, and the
   cheapest one is the adjacent array entry.** (The instrument itself bit first:
   `sn.__dict__["areas"] = …` is a silent no-op because `areas` is a `property` —
   write the backing field `sn._areas`.)
9. ⭐ **The cylinder is blind to the Gram question — a THIRD member of the family.**
   `[M]` `R/ΔA` is bit-exactly `diag(1, 1/3)` on the cylinder (= the shipped
   moment-axis metric, so it reads as "not its own object") and has off-diagonal
   `h/(6r_c)` on the sphere. Sharper still: **the cylinder's MASS Gram `M/V` IS the
   sphere's REDISTRIBUTION Gram `R/ΔA`, exactly** — the two objects an ownership
   argument is trying to separate are the same matrix one geometry over. Joins the
   known cylinder blindnesses (the slope row reading `0 == 0`; σ_y-folded `β`).
   ⟹ **any curvilinear Gram/measure claim must be witnessed on the SPHERE.**

