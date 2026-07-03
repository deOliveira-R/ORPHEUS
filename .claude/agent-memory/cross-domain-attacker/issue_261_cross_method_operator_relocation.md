---
name: issue-261-cross-method-operator-relocation
description: #261 — relocate C/F/S cross-method operators sn/→transport/. Verdicts on the 4 structural hypotheses (C's sn_mesh asymmetry, the identity-vs-shape guard, coefficient-vs-mat_xs, shared-base) + the __add__ one-sided-dispatch resolution.
metadata:
  type: project
---

# #261 cross-method operator relocation (sn/ → transport/) — structural verdicts

Branch `refactor/operator-inverse-algebra`, main-agent-direct carve. Challenged a
code-verified 4-hypothesis brief on C (`CollisionOperator`), F (`FissionOperator`),
S (`ScatteringOperator`). Branch-verified reads (L-005), NOT Nexus.

**Why:** feeds the #261 relocation of the three cross-method reaction operators
out of `sn/` (which only `transport`-imports them) into `transport/`.
**How to apply:** these are the structural verdicts the carve referenced; reconcile
file:line against the live worktree before acting (process-discipline: merge-status
in memory goes stale).

## The shared morphism (the backbone for all three)
`flux-STATE → rate-density SOURCE/SINK` = the LINEAR PART of an affine bundle map
`E_flux → E_source`, fiber gain = the cross-section (1/cm). This is the
`LinearOperator[Flux, SourceSink]` instantiation (Codomain≠Domain), documented in
`operator.py:370-388`. All three ALSO satisfy `IntegralKernelOperator`
(F via kernel+production_rate, S via kernel). This is a PROTOCOL the three share,
NOT a base class — see H4.

## H1 — C's `sn_mesh` is a SESSION-BORN ASYMMETRY → CONFIRMED
The base `MultiplicationOperator` is ALREADY mesh-free (`multiplication_operator.py:76-80,118`;
`apply` reads `mesh = psi.bulk.mesh` at `:206`). C re-bolted `sn_mesh`
(`operator.py:600`) for exactly three uses, none a genuine apply-time SN-geometry tie:
(a) the legacy ndarray→CrossSectionField wrap (`operator.py:585-594`) — production
BYPASSES it (`solver.py:223/936/1013` all pass `mat_xs.total_cross_section_field`, a
typed `CrossSectionField`); (b) the W-D `domain`/`codomain` = `sn_mesh.full_field_space`
(`operator.py:366`) — but the guard reads only `(name,shape)`, which F/S supply via a
DIRECTLY-HELD `full_field_space` field WITHOUT a mesh (`fission.py:187`/`scattering.py:437`);
(c) the identity invariant (= the artifact, H2). C predates W-D; the W-D comment on S
already tags `full_field_space` "D5 cross-method-honest handle, relocation-ready for
#261" (`scattering.py:436`) — C just never got the same treatment.
**FIX:** C holds `full_field_space` directly (Smell #16 shape-1 collapse — two storage
conventions for ONE quantity, the composite space); the legacy wrap → a
`from_solver_data(*, coefficient, full_field_space=None)` + `from_mesh(sn_mesh, sigma)`
back-compat classmethod, symmetric with F/S.

## H2 — `streaming.sn_mesh is diagonal.sn_mesh` redundant with the W-D guard → REFUTED (remedy), and the instinct is also imprecise
The W-D `OperatorSum` guard checks `(name,shape)` (`operator.py:743/750`;
`FullFieldSpace.__eq__` is `(name,shape)`, block spaces `compare=False`,
`full_field_space.py:70-82,152`). So the guard accepts two DISTINCT-but-same-dimension
meshes — it does NOT subsume object-identity. "Downgrade `is` to `==` because the guard
enforces it" is FALSE: the guard enforces only SHAPE-equality.
The THREE invariant strengths are NESTED:
`object-identity ⊋ geometric-consistency ⊋ shape-equality`.
The ESSENTIAL requirement is the MIDDLE: `InvertibleOperator.solve` threads `self.sigma`
(from `diagonal`, `operator.py:842`) into a sweep whose geometry/scheme/volumes/streaming
ALL come from `self.mesh` (= `streaming.sn_mesh`, `loss_representation.py:1690-1723`, every
`self.mesh.*`). Two meshes of equal shape but different VOLUMES/curvature → shape-guard
passes, σ-field fuses with wrong geometry → silent wrong answer. The current `is`-check is
a blunt proxy for geometric consistency; shape-equality is too weak.
**FIX:** when C drops `sn_mesh` (H1), the `is`-invariant cannot be expressed at all — REPLACE
(not downgrade) it. The natural replacement: `C.coefficient.mesh is L.sn_mesh` — the
`CrossSectionField` ALREADY carries `.mesh` (`cross_section_field.py:85`), so provenance is
recoverable without C holding a redundant `sn_mesh`. Discriminating test: L on mesh_A, C on a
coefficient `from_mesh(mesh_B)` with same shape / different volumes → assert raise OR (if no
raise) `not allclose` vs all-mesh_A reference. A guard-only-checks-shape test PASSES the
mismatch (proves the guard insufficient); the `is`-check REDs it.

## H3 — C holds `CrossSectionField` (not `mat_xs`) GENUINELY → CONFIRMED
C is a general multiplier `M[f]` whose `f` is a single `(ng,*spatial)` field that may be
σ_t OR the algebraically-derived σ_r = σ_t − Σ_s0. `CrossSectionField` carries the
cone/vector-space algebra (signed differences are valid coefficients,
`cross_section_field.py:38-42`) precisely for that derived-σ_r case. F/S read `mat_xs`
because their action is irreducibly PER-MATERIAL (rank-1 χ⊗νΣf; per-ℓ Σ_{s,ℓ}). So the
coefficient-vs-mat_xs split is principled. CAVEAT: σ_r TODAY is `dict[int,np.ndarray]`
(`mat_xs.foldable_sigma()`, `scattering.py:1198`) and NO production site assembles a σ_r
`CollisionOperator` — the #200 fold is documented (`operator.py:702`) but UNWIRED. So
"C serves σ_r" is a LATENT (designed-for) consumer, not live. When #200 lands, type the σ_r
path as a `CrossSectionField` to keep C's coefficient typed.

## H4 — shared `ReactionOperator` base for the three? → DO NOT MINT (L-004 type-theatrics)
The shared morphism is real (flux→source affine-bundle linear part) but the three `apply`
bodies share ZERO code: C = diagonal `M[σ]` (base supplies apply/solve/apply_T once, in
`transport/`); F = rank-1 `χ⊗νΣf` (`TensorProductOperator`); S = non-local `R∘Λ∘M` + the
SO(3) eigenbasis frame (W-E Funk–Hecke: the angular frame is scattering-OWNED, so
frame/quadrature/scattering_order is irreducibly S's). A contract with no shared body is a
PROTOCOL, not a base — and it ALREADY EXISTS: `IntegralKernelOperator` (all three satisfy) +
`LinearOperator[Flux,SourceSink]` (the typed morphism, half-built = P4.5 W-F). An intermediate
ABC hosts no method → ceremony with no referent. F = scattering's rank-1 degenerate (memo
`flux_to_sourcesink_operator_contract_frames.md`). Keep C on `MultiplicationOperator`, F/S on
`LinearOperator` direct. UNIFY via the shared `[Flux,SourceSink]` typing + the existing Protocol.

## The `__add__` layering after relocation → one-sided dispatch on StreamingOperator is FORCED
`transport/ → sn/` is an illegal upward import, so a relocated C cannot `isinstance` on
`StreamingOperator`. Resolution (asymmetric, already half-correct):
- `StreamingOperator.__add__` (stays `sn/`) KEEPS the dispatch (`sn/` importing C from
  `transport/` is legal downward) — already `operator.py:483`.
- `CollisionOperator.__add__` (moves `transport/`) DROPS its `isinstance(StreamingOperator)`
  arm (`operator.py:653`) → falls through to generic `OperatorSum`.
- `C + L` ordering: recovered by `StreamingOperator.__radd__`. BUT the base
  `LinearOperator.__radd__` (`operator.py:405`) returns a plain `OperatorSum(C,L)`, NOT an
  `InvertibleOperator` → `StreamingOperator` MUST OVERRIDE `__radd__` to mirror `__add__`
  (`InvertibleOperator(self, other)` when other is a CollisionOperator). That override is the
  ONE new method the carve requires.
**Invariant after carve:** ALL method-specific algebra dispatch lives on the `sn/`-side
`StreamingOperator` (`__add__`+`__radd__`); the relocated C/F/S carry only generic algebra.
Same arrow-direction the resolvent backbone predicts (L-007): the SN-specific `(L+C)⁻¹≈sweep`
identity is StreamingOperator/InvertibleOperator's to own; C is a generic multiplier that
merely participates.

## Minimal symmetric constructor shape for the three (the deliverable)
Hold typed data + OPTIONAL `full_field_space: FullFieldSpace | None`; read the mesh off the
CARRIER at apply time. F/S already conform; C is the outlier (re-bolted a mesh the
`transport/` base deliberately shed). The carve = make C conform to the pattern
`MultiplicationOperator` ALREADY established and F/S ALREADY follow — NOT design a new base.
Uniform constructor family (greppable): `C/F/S.from_solver_data(*, <data>, full_field_space=None)`.
