---
name: issue-257-s6-integral-kernel-category
description: #257 S6 IntegralKernelOperator §5.6 Kernel refinement Protocol + fission production_rate + scattering kernel — PASS-WITH-NITS verdict
metadata:
  type: project
---

# #257 S6 — IntegralKernelOperator (§5.6 Kernel suffix) — PASS-WITH-NITS

`feature/field-typed-operator-algebra`, HEAD `7b463bb`, UNCOMMITTED, ADDITIVE + BIT-IDENTICAL
(matvec arms untouched). NEW `orpheus/transport/integral_kernel_operator.py` (L2) +
`fission.production_rate` property + `scattering.kernel` property + transport export + docs §5.6 STUB.
Gates verified GREEN during review: 3 S6 test files (20✓: scattering/fission crosscheck +
`test_integral_kernel_category.py`), fission+scattering matvec-arm regression (85✓ →
additive/bit-id confirmed). Closes the §5.6 suffix-law partition: Operator (S3b local
MultiplicationOperator) / **Kernel (S6)** / Functional (S5 disjoint sibling).

**Why:** the campaign's §5.6 middle term. `IntegralKernelOperator(LinearOperator[V], Protocol[V])`
is a REFINEMENT of LinearOperator (still apply+caps, ADDS `kernel` property) — asymmetric to the
S5 Functional which is DISJOINT (speaks `evaluate` not `apply`). Discriminator = LOCALITY (Frame 3):
local/diagonal MultiplicationOperator has NO kernel→NOT a Kernel; nonlocal integral action
(groups for fission, ordinates for Pℓ scattering) HAS a kernel. Two named kernels: fission =
rank-1 `χ⊗νΣf` TensorProductOperator (Wave T); scattering = `R∘Λ∘M` OperatorProduct (S6).

**How to apply:** 4 decisions ALL PASS, 2 do-now doc NITs.

- **D1 `cast(LinearOperator,...)` ×3 in scattering.kernel (`:633-651`) — PASS.** Right tool, masks
  nothing. VERIFIED to root: `ProjectionOperator`/`ReconstructionOperator` (`projection.py:108,208`)
  + `LegendreMomentScattering` (`scattering.py:210`) subclass the UNPARAMETRISED `LinearOperatorMixin`
  (`operator.py:352` `Generic[V]`) with concrete `apply(x: np.ndarray,/)` → pyright genuinely can't
  unify `OperatorProduct[V]`'s V. Precedent `orpheus.data.micro_xs` is a PACKAGE (not .py) — `cast(`
  real at `gendf.py:283,390`/`hdf5_io.py:72-81`, same shape. `cast`>`# type:ignore` (typed assertion
  not suppression, keeps return type honest). The wrong-factor-order risk is RUNTIME, pinned 0-ULP by
  the crosscheck, NOT maskable by cast. Parametrising the mixins = #226 (correctly out of scope).
- **D2 scattering.kernel PARALLEL to `_aniso_source_from_moment_values` (Cardinal-2) — PASS, do NOT
  route arms through kernel.** The note-#2 verified-identical-twin shape, ACCEPTABLE here AND routing
  would be a NET NEGATIVE: (1) `/W` lives at the apply boundary by design (L18/Pattern-7;
  `build_aniso_source:937` divides by sum_w AFTER the map; `kernel.apply` returns pre-`/W`) — routing
  pushes `/W` into the kernel (breaks §5.6 bare-integral reading) or keeps a `/sum_w` wrap (saves
  nothing, adds per-call property-build indirection on the hot matvec); (2) kernel is a FRESH
  3-operator alloc per `@property` access (read-through) → routing imposes that on hot path; (3)
  collapse trigger NOT met (exactly 2 assemblies, leaves are single source, 0-ULP gate reddens on
  drift). Twin is single-sourced AT THE LEAVES (same `M`/`Λ skip_l0`/`R`); kernel = 2nd ASSEMBLY of
  those leaves. `OperatorProduct(R,OperatorProduct(Λ,M)).apply(ψ)=R(Λ(M·ψ))` (`operator.py:826`
  `a.apply(b.apply(x))`) ≡ `build_aniso_source` `M.apply(ψ)` then
  `_aniso_source_from_moment_values=R.apply(Λ.apply(.))` — byte-id pre-`/W`. STANDING: confirm the
  semantic-parallel; the bug habitat (future metric-weighted projection lands in arm not kernel) is
  pinned closed by the 0-ULP equivalence gate.
- **D3 `skip_l0=True` partial kernel — PASS honest.** "kernel is the ℓ≥1 part" honest per locality
  criterion. ℓ=0 P0 in-scatter + (n,2n) documented as LOCAL/separate components of full apply.
  `scattering_order==0→ValueError` (`:653-659`) makes "no aniso kernel" an honest raise (Pattern-4),
  P0-only op correctly NOT an IntegralKernelOperator.
- **D4 Protocol refinement form — PASS elegant+correct.** `LinearOperator` IS `Protocol[V]`
  (`operator.py:260`) so re-declaring `Protocol[V]` in bases is REQUIRED (Protocol-subclassing-Protocol
  idiom, else becomes concrete ABC). Invariant V verbatim correct (`apply(x:V)->V` endo per #256-S3,
  invariant by dual use). `@runtime_checkable` live-verified. `op.kernel`/`op.production_rate` read as
  the §5.6 suffix law. Master standard met.

**Verification SHAPE (note-#3 textbook):** fission crosscheck B.1 = CORRECTNESS hand-derived double-loop
`hand_derived_fission_emission` (NO shared numpy primitive) + role-swap foil; B.2 = `χ·production_rate
≡F.apply` 0-ULP de-risk (Mode-11, reads off live op). Scattering crosscheck = EQUIVALENCE-only de-risk
(`S.kernel.apply≡_aniso_source_from_moment_values(M·ψ)` 0-ULP, Mode-11) + non-degeneracy guard (ℓ≥1
moments non-zero) — physics L1 backing is the EXISTING aniso MMS gate
`test_curvilinear_aniso_scattering_p1.py` (NOT this file). Equivalence-only legit IFF structural physics
anchor exists — it does.

**NITS (do-now, doc-accuracy — the campaign's densest recurring tell):**
- ⭐ N1 (`fission.py:319-326`) production_rate docstring ASSERTS cross-module byte-identity as a
  structural invariant ("same numpy primitive, same axis, same keepdims, 0 ULP" between
  `ProductionRateFunctional.evaluate` in transport/ and `RankOneOperator.apply` in numerics/). 0-ULP is
  a property of THE GATE (B.2) not THE PROSE — a future `RankOneOperator.apply` einsum rewrite stays
  principled-equiv but makes the prose a lie. FIX = demote to test-pinned current fact ("agrees
  byte-for-byte AS PINNED BY test_fission_kernel_crosscheck.py B.2; both currently run...") — bring
  fission to scattering's standard (scattering `:606-609` already cites the gate correctly). The
  recurring docstring-asserts-false-coupling tell ([[issue-236-phase2-stepc-tau-retirement]],
  [[issue-257-s3a-diagonal-engine]]).
- ⭐ N2 (`scattering.py` `_aniso_source_from_moment_values` docstring ~`:492`) — the production twin has
  NO reciprocal pointer back to `kernel`. Per note-#1 standing remedy (reciprocal twin-cross-ref +
  tracked trigger, NOT unification). FIX = add "the §5.6 `kernel` property composes these same R∘Λ
  leaves with M; pinned bit-id by test_scattering_kernel_crosscheck.py; edits here must keep that gate
  green."
- N3 (opt, docs todo) — `SumOfTensorProductsOperator` un-orphaning + carrier-agnostic core relocation
  deferred "per S6 plan" (plan-file pointer = session-scoped) → confirm/file an issue (Cardinal-4,
  anti-#11). Not blocking (stub honest about being a stub).

**STANDING TELLS reinforced:** (a) a `@property` that FRESHLY allocates operators per access is
read-through-correct (frozen-lens, can't go stale) but must NOT be routed into a hot matvec — distinguish
the semantic-reading property from the realization arm; (b) equivalence-only de-risk gate is legit IFF a
structurally-independent physics anchor exists elsewhere (here the aniso MMS gate) — verify it before
accepting a 0-ULP-only crosscheck; (c) cited `cast` precedent may be a PACKAGE not a .py — `find -path`
not `Read` the asserted file; (d) bring sibling docstrings to the same gate-citation standard (fission
lagged scattering on the byte-id-as-asserted-invariant tell).
