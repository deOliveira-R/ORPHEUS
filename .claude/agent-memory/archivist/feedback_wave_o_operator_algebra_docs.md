---
name: wave-o-operator-algebra-docs
description: Consolidated playbook for documenting the Wave-O SN operator-algebra refactor family (Issue #208/#201) into operator_algebra.rst — sibling-operator extraction, operator-output role-typing, face-locus cochain typing, fold→variadic retirement, composite-source API, affine field algebra, layered eigenvalue + G-adjoint, and landed-carve guard-lifts. Each Wave-O step lands a peer ====-section in the bc-extraction family with a shared section template; this note carries the per-step distinctive kernels.
metadata:
  type: feedback
---

# Wave-O operator-algebra docs (Issue #208/#201, branch refactor/field-role-typing)

The whole Wave-O SN role-typing refactor is documented in ONE canonical
home: `docs/theory/operator_algebra.rst`, as a growing family of peer
`====` sections rooted at `.. _bc-extraction:`. The per-solver page
`discrete_ordinates.rst` gets only Key-Facts bullets + retraction notes
that cross-ref the rich content (NEVER duplicate the derivation across two
pages). Universal build-gate / cross-ref / venv discipline now lives in
`AGENT.md` ("Build-Gating & Cross-Ref Reality") — this note carries only
the doc-craft that is specific to these refactors.

## The shared section template (every Wave-O step uses a variant)

A peer `====` section appended to the bc-extraction family, made of:
1. **Lead** — commit chain + Issue + date + a one-liner placing the step
   relative to its predecessor ("the X-locus sibling of …", "completes
   the Y half; the X half landed in `<commit>`", "ORTHOGONAL to …").
2. **Key Facts admonition** (`:class: tip`/`important`) — 5–7 bullets, the
   discoverable summary; ALWAYS includes one honest-scope bullet.
3. **Labeled equations** — each `.. math:: :label:` + a `.. vv-status:
   <label> documented` comment (structural framing, not solver claims).
4. **A list-table** that is the canonical at-a-glance index (role grid /
   block-occupancy / posing table / two-loci / retyped-sites).
5. **The load-bearing rationale** — the WHY, usually a rejected-alternative
   or design-correction catalog (this is the Cardinal-Rule-3 payload).
6. **Numerical-evidence list-table** — bit-identity / principled-equivalence
   grounds, with the closed-form anchor + the "MMS does NOT prove the
   eigenvalue" caveat where an eigenvalue is in play.
7. **Honest-scope `.. warning::`** — what the step does NOT do, framed as
   the attacker's ABSENCE, not a hedge.
8. **Cross-references** — to predecessor sections, the eigenvalue posing
   section, and the generic field-algebra principle in api/numerics.rst.

## Breadcrumb-resolution discipline (the recurring half of every step)

Each Wave-O step DISCHARGES a "minted, no consumer until O.X" / "deferred
to O.4b" / "what remains for O.2" forward-pointer that an EARLIER step
left. Half the doc job is writing the new section; the other half is
RECONCILING the now-stale "future"/"will get for free" framing in the
consumer sections. The idiom: flip the deferral cell/paragraph to
"landed in <step> (:ref:`new-anchor`)", and where the predecessor's prose
was written from ITS vantage, add a forward-pointer clause ("That consumer
has since landed: …") rather than a wholesale rewrite — preserve the
predecessor-era narrative, remove only the now-false present-tense claim.
Grep the breadcrumb token across the WHOLE `docs/` tree, not just the
named pages (the multi-page O.4b-deferral staleness was a 3-PAGE flip:
`discrete_ordinates` + `operator_algebra` + `boundary_conditions`; the
surgical plan named only 2).

## Per-step distinctive kernels (the load-bearing content per step)

### O.4a.2 — sibling-operator extraction (−B BC-extraction)
Promotes a buried sub-action (the `inflow = bc.apply(outflow)` "keystone"
inside the sweep) to a first-class sibling operator on a direct-sum state
`V = V_bulk ⊕ V_inflow ⊕ V_outflow`. Section gets: block matrix with
under/over-braces by block role; the deleted keystone as before/after
code (call out what SURVIVES and WHY — the curvilinear pole seed survives
because it is the r=0 regularity condition, NOT a BC); a **design-correction
catalog** (3 traps: keep-outflow-defect-not-raw / B-projects-to-inflow-row /
seed-from-rhs.boundary). New leaf `SNBoundaryOperator` automodule'd (its
docstring has no `:label:` → safe). Retraction `.. note::` atop the
superseded per-solver Phase-C prose (preserve the trace-contract insight —
the 2-D path still uses it). The carve DELETED `_transport_operator_matvec_unified`
→ grep+repoint 6 dead `:func:` refs.

### B.5.2 — operator-output role-typing (.apply → SourceSink)
TYPE-ONLY retype completing the role grid (`.apply`→SourceSink /
`.solve`→Flux / `from_balance`→Residual). Load-bearing content = the
**two-hat tension + rejected-alternative**: B's output would wear residual
(Krylov matvec) OR source (SI rhs); the class-gate throws on
`Residual + SourceSink` the moment SI rhs adds them. Choosing SourceSink
for ALL outputs dissolves it — B wears ONE hat, BOTH sums close (one
`\underbrace{}` eq showing both closed sums). Quote the USER's verbatim
governing principle ("A residual only arises after we compare an operator
output against something else and get a defect"). The flat-round-trip
escape hatch (`to_flat` ravels type-agnostically; `from_flat` rebuilds off
the FLUX template) is WHY the Krylov path stays Flux. Cite the throwaway
decision instrument `diag_b52_*.py` (READ it; OPT-BSS closes both sums,
OPT-BR throws). The generic class-gate already lives in api/numerics.rst
`_field-algebra:` (`r = Aψ − q`) — cross-ref, don't re-derive.

### #205 Phase 5 — face-locus cochain typing (WavefrontFlux)
Types the interior face state the sweep propagates as `WavefrontFlux`
(primal 1-cochain C¹_int). Biproduct `C¹ = C¹_int ⊕ C¹_∂` mirrors the
cell+trace `V_bulk ⊕ V_boundary` ONE LOCUS DOWN — the KEY insight: the
boundary summand is the SAME at both loci (it IS the domain-edge faces =
the cell-state's trace); only the interior summand differs; boundary
persists, interior is ephemeral. Name the `ι_*`/`ι*` trace operators +
the 2 biproduct laws (`ι*∘ι_*=id` = "absorption=identity", now PROVABLE;
`π_int∘ι_∂=0`). Flux-only-single-role rationale: the interior face flux is
the off-diagonal of the per-octant lower-triangular L_oct, re-formed each
sweep ⟹ transient ⟹ no source/residual role (illegal-states-unrepresentable:
no InteriorFaceResidual). The attacker's field+views-NOT-per-face rejection
(3 grounds). This is the step that owned the 3-PAGE O.4b-deferral staleness
flip.

### O.2a — fold→variadic retirement (honest L+C−S−B)
Retires the transitional `(L, S, F)` triple with B folded into the S slot;
drivers become variadic `Driver(L_resolvent, *gains)`. Name the deleted
symbol (`_scattering_with_boundary_op`) as a PLAIN double-backtick literal
(NOT `:meth:` — it's gone). WHY-variadic = the fixed triple encoded a false
posing-layer distinction at the iteration layer. WHY-B-separate (numbered):
(1) the `|Ω·n|·w` adjoint metric lives on B's TRACE domain not bulk; (2) B
can't join the L+C preconditioner because `OperatorSum` drops `CAP_SOLVE`.
A fold→variadic reassociation `L−(S+B) → (L−S)−B` is a SECOND
principled-equivalence instance (document ALONGSIDE the extraction's drift,
don't overwrite) — cite the 3 criteria + the commit's exact ULP numbers.

### Composite-source API (solve_sn_fixed_source(composite))
NEW ergonomic public form that REUSES an existing carrier (`TimedFullField`,
q = q_bulk ⊕ q_∂) and retires a previously-documented operator-triple
bypass. Load-bearing pedagogy = "we already have the right object/concept —
just better ergonomics to GENERATE it" (Cardinal Rule 2): enumerate what
PRE-DATED the work (carrier, leaf types, affine law, q.boundary slot) and
name the TWO things that were missing (the `prescribed_inflow` generator +
the public entry point). Flip the 4.6-MMS "bypass" narrative to past tense
(retirement = test migration — the new code is what gets tested); KEEP the
source→flux type-bridge concept, reframed as "now INTERNAL to the solve".
The `prescribed_inflow` (SHIPPED known-arrays) vs `from_spec` (DEFERRED lazy
recipe) distinction is unify-after-two, a distinct path not a duplicate.

### O.2 — affine field algebra (FluxDisplacement + typed residual)
The state ⊖ state torsor (affine space A over difference space V) + the
typed `from_balance` residual (box 7). The torsor 4-line `:label:` + the
why-affine rationale (pre-carve `flux+flux` was a literal affine-axiom
violation = the SAME sin #201 fixes for the residual). The load-bearing
NUMERICAL content = the true-error diagnostic: ‖Δψ‖ understates the true
error by 1/(1−ρ) (~100× at c=0.99) — the c→1 false-convergence fix; the
contraction-ratio + true-error `:label:`s. Document the torsor round-trip
subtlety: `ψ₁⊕(ψ₂⊖ψ₁)=ψ₂` is exact in real arithmetic but NOT bit-exact in
IEEE-754 (a+(b−a)≠b; test asserts nulp=8; the individual mint + add ARE
bit-exact) — frame as bit-identity-vs-principled-equivalence at the algebra
level. The ρ≈c contraction test is `l1` closed-form; 1-group is ACCEPTABLE
there because ρ=c is a flux-shape-INDEPENDENT *rate* claim (state this
verbatim so a future reader doesn't flag a Cardinal-6 violation). Pull
numerical bounds from the TEST FILES, not the brief's memo estimates.

### #205 Phase 5a/5c — moment-reduction of the SI iterate
ORTHOGONAL perf carve (not a typing): the SI fixed point lives in moment
space (S consumes ψ only through Mψ), so hold the iterate as
`HarmonicMomentField`, N→(L+1)(2L+1). 5c moves the projection INTO the
anti-diagonal sweep walk (`moment_buf[…] += einsum("nlm,ngd,n->lmgd", …)`),
3.06× peak win. The load-bearing V&V split: the per-STEP source is
BIT-IDENTICAL (M built from S's own quadrature ⇒ 0 ULP, pinned by a
de-risk probe with an INDEPENDENT Bell&Glasstone hand reconstruction at
Q2b); the ONLY non-bit-id change is the SI convergence test moving to
moment-L2 = MORE principled, not a regression. The Y₀⁰=1 convention (ℓ=0
moment IS the scalar flux, read with no rescale). Geometry restriction
(curvilinear NO — the M-M Carlson per-ordinate seed; Krylov NO — iterates
full bulk). The fuller-view post-projection oracle is RETAINED
(aggressive-retirement exception) because the ℓ=0 scalar cross-check is
BLIND to ℓ≥1 (Mode-5/Mode-2 catcher) — `.. warning::` it is procedurally-
not-structurally independent (shared cell kernel), so pair it with the
SI≡Krylov anchor.

### Layered eigenvalue architecture + composite G-adjoint (R5)
These two are COUPLED — write them together. The `eigenvalue-posing`
section (4-layer: leaves/posing/resolvent/algorithm; standard form
Aψ=λMψ; resolvent A_loss⁻¹M; K-row LIVE, α/adjoint/transient/full-spectrum
as documented FUTURE SEAMS with zero scaffolding) leaves THREE abstractions
that R5's composite G-adjoint makes concrete: the "adjoint row (future
seam)", the "metric lives at the leaf" subsection, the Layer-1 G-metric
row. R5 lands `op.H = G⁻¹AᵀG` over a `FullFieldSpace` (bulk⊕boundary) →
the job is half NEW-section (the A†=G⁻¹AᵀG derivation from reciprocity +
per-block ∫dV∫dΩ / ∫_Γ|Ω·n| measure derivations + the pseudo-inverse
exactness proof for the singular trace) and half RECONCILE the eigenvalue
page's stale "future"/"will get for free" framing to "already built and
verified by R5". The C/S/F=None capability-lattice is the load-bearing
subtle point (OperatorSum.domain first-non-None propagates L's metric →
C is G-conjugated for free; C carrying it too would be a double-apply bug;
(L+C-S-F-B).H raises = illegal-states-unrepresentable). The dense-oracle
(`diag_p42_adjoint_oracle.py`, explicit dense Aᵀ + explicit diagonal G,
knows NOTHING of FullFieldSpace.apply_metric) is the structurally-independent
V&V ground — RUN it. power_iteration is CANONICAL not deprecated (general
admits monolithic CP/diffusion matrices with no triple).

### Landed-carve guard-lift (2-D SI Phase A)
Documenting a landed carve that un-gates a stale `NotImplementedError`:
the theory page states the CURRENT invariant + production decision +
architecture WHY, NOT the "we found a stale guard" story (L10, user
explicitly steers this). The only nod to history = a one-line retired-legacy
note so a future grep for the old token (`B1''`) lands somewhere. The stale
claim often lives in the xfail REASON STRING (baked into graph.json only) or
the solver docstring, NOT the source RST — confirm each grep hit's TOPIC
before assuming it's the stale claim; the real gap is often an ABSENCE (the
invariant simply not yet stated), not a wrong sentence to flip.

## Alias-retirement two-sense collision (O.4a.1, a sibling cleanup)
When a retired token (`BoundaryOperator` geometry alias) COLLIDES with a
live forward-looking symbol of the same name (the `numerics.operator`
block-role marker), the sweep is DISAMBIGUATION not rename. Triage by
module: `geometry.boundary.*` ref → canonicalize; `numerics.operator`
block-role ref → LEAVE. Four intentional-remaining-match classes
(rename-table left column + flip surrounding prose / pre-refactor history /
MARKER refs / colliding teaching-example names). Flip the class-docstring
"preserved as deprecated alias" → "retired in O.4a.1, sole importable name".
Collapse stale-comment rationales whose entire premise was the now-gone
alias-equality. Deleted `:class:` xref targets do NOT warn without `-n` —
fix on correctness grounds via grep gate, not the build.

## Source-of-truth discipline (all Wave-O steps)
NO SymPy (architecture, not constructive math). The narrative descends from
the rich CODE docstrings (module/class/property) + the throwaway diagnostic
instruments (`derivations/diagnostics/diag_*.py` — READ and RUN them to
confirm cited numbers) + the commit bodies (which carry the exact ULP
numbers + rationale — QUOTE them) + the cross-domain-attacker frame memos
(where the hardest "why" lives, e.g. the linear-μ-fully-activates argument).
These ARE the algebra-of-record here. Quality scores ran 5/5 across these
tasks except derivation-source (4, no SymPy by nature) and the API-ergonomics
steps (derivation-depth 4, no new symbolic derivation).

Related: [[feedback_nonvacuum_mms_stub_expansion]] (the 4.6 MMS stub this
family's composite-source API drives), [[feedback_bc_trace_law_wave_12]]
(the predecessor 12-wave BC refactor), [[feedback_post_wave_cleanup_docs]]
(the post-wave cleanup close-out variant), [[feedback_canonical_convention_page]]
(the #196 storage-layout flip that PR-INDEX-5 `(N,ng,nx,ny)` rides on).
