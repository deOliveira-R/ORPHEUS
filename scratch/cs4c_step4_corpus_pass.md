# CS4c step-4 corpus pass — findings memo (archivist, 2026-08-31)

Branch `refactor/cs4c-step4-fission-binding`, HEAD `fadad026`.

## Build baseline (re-measured this session, MANDATORY per AGENT.md)

`.venv/bin/python -m sphinx -E -W -b html docs docs/_build/html`
→ **EXIT=0**, `WARNING:|ERROR:|CRITICAL:|SyntaxWarning` count = **0**.
`directives: wrote 412 edges from pending registry`;
`directives: 18 statement(s) declared to have no implementation`.
Acceptance gate: unchanged set (0/0), EXIT=0.

## Census denominator

107 `.rst` files under `docs/` (excluding `_build`), scanned in Python
(`scratch/_cs4c4_census.py`), positive control `.. math::` = 1327/70 files.

| spelling | occ | files |
|---|---:|---:|
| `FissionOperator` | 50 | 19 |
| `IsotropicFission` | **0** | 0 |
| `FissionMaterialField` | **0** | 0 |
| `full_fission_kernel` | **0** | 0 |
| `foldable_sig_s` | **0** | 0 |
| `multiply.outer` | **0** | 0 |
| `FissionKernel` | 3 | 2 |
| `N2NOperator` | 9 | 7 |
| `N2NMomentOperator` | 6 | 2 |
| `production_rate` | 14 | 5 |
| `L+C-S-B` (no N2N) | 118 | 32 |

⭐ The brief's named suspects `foldable_sig_s` / `chi` / `sig_p` read-through
are **0 sites in `docs/`** — those were code-side docstring retirements only.
The real corpus work is the two `.. implements::` declarations whose
transcribed bodies moved, and the F/S/N2N unification the step created.

## Live-code ground truth (verified this session, `hasattr` + `dataclasses.fields`)

- `FissionOperator` fields = `(domain, codomain, energy, frame)`.
  **GONE**: `.chi`, `.sig_p`, `.mat_xs`. `kernel`/`production_rate` survive as
  DELEGATIONS to `energy`.
- `IsotropicFission` fields = `(domain, codomain, fission)`. `IsotropicFission.chi`
  does **not** exist — χ is `self.fission.gather_chi(bulk.shape[1:])`, a GATHER.
- `MaterialXSField.foldable_sig_s` **GONE**; `fission_production_field` survives.
- `N2NOperator` fields = `(domain, codomain, energy, frame)`; `weights` **GONE**.


## Decisive measurement — the "same Python classes" thesis

AST census of production instantiation sites (`orpheus/**/*.py`, `ast.Call`
whose func is a Name or `Cls.classmethod`):

| class | sites | packages |
|---|---:|---|
| `FissionOperator` | **1** | `sn` only (`sn/solver.py:2910`, the eigen-M posing) |
| `IsotropicFission` | **4** | `diffusion`, `homogeneous`, `sn`, `transport` |
| `MultiplicationOperator` | 3 | `diffusion`, `homogeneous`, `sn` |
| `N2NOperator` | 2 | `sn` |
| `IsotropicScattering` / `IsotropicN2N` | 3 / 3 | `diffusion`, `homogeneous`, `transport` |

⟹ THREE sites claim `FissionOperator` is "the *same Python class* instantiated
by SN, diffusion and the infinite-medium solver" — `path_integral.rst:63`,
`path_integral.rst:361`, `foundations/index.rst:12` (plus the
`path_integral.rst:38` machine header `all_three_consumers`). All are now
present-tense-FALSE for F. The thesis SURVIVES and sharpens: the shared class
is `IsotropicFission` (the energy binding); SN adds an angular conjugation.
`MultiplicationOperator`'s half of the same sentence is still true (3/3).

## Fix list (by page)


### Files changed (19 authored + 1 generated)

| file | rationale |
|---|---|
| `docs/theory/foundations/operator_algebra.rst` | `fission-as-dyad`: the two `.. implements::` transcribed bodies were code that no longer exists (`self.chi` is a RETIRED property; `self.mat_xs` gone). Migrated the arithmetic-home declaration to `IsotropicFission.kernel`/`.production_rate`, kept the two `FissionOperator` delegations declared, corrected the χ-column claim ("an array on the operator" → a gather), Source-map note re-pointed, and a dated note on `emission-kernels-btd` (the lens became shared faces). |
| `docs/theory/methods/sn/adjoint.rst` | NEW H2 `sn-fission-binding-adjoint` (2 new eq-labels) — the S/F/N₂ₙ one-shape unification, the factor-reversal transpose, the energy-vs-angular distinction; corrected the N₂ₙ reciprocity closer (evaluated on the frame form, not a w-embedding); measured bit-identity vs ULP note; `F.H` bullet re-derived. |
| `docs/theory/foundations/infinite_medium.rst` | ⛔ **published code block was FALSE** (`FissionOperator.from_solver_data` → live `IsotropicFission.from_material_xs`); two prose sites; `mg-balance` + `two-group-F` declarations re-pointed. |
| `docs/theory/methods/sn/slab_multigroup.rst` | §sn-scattering-fission-operators spelled `A = L+C−S−B` twice while the same page's §480 states `−N₂ₙ` (self-contradicting); `FissionOperator.apply` scalar contract; `multigroup` set 11→12; Key Facts + normalization-chain pointers. |
| `docs/theory/foundations/path_integral.rst` | The "same Python classes" thesis, machine header + 2 prose sites, with the AST census as evidence. |
| `docs/theory/foundations/index.rst` | Same thesis at the foundations root. |
| `docs/theory/foundations/operator_adjoint.rst` | "all five operators carry the same `full_field_space`" — precise note: true of the composite posings, false of the k-outer's binding. |
| `docs/theory/methods/sn/solver.rst` | `SNSolver.fission_op` was `FissionOperator` — live: `IsotropicFission` on the bulk space. |
| `docs/theory/methods/sn/index.rst` | Machine header `composites.A` + the chapter-level `(n,2n)` caveat note (the 37-site pedagogical spelling declared, not silently swept). |
| `docs/theory/conventions/indexing_and_layout.rst` | Layout table split into the two bindings (the scalar row was attributed to the class that now refuses scalars); the F row's dead `to_scalar`. |
| `docs/theory/foundations/frame.rst` | `energy-condensation-chi-collapse` gains the G-F1 gate cross-ref + a NEW `.. important::` stating the χ↔νΣf COUPLED law (1 new eq-label) with the three measured controls and the activation precondition. |
| `docs/theory/methods/sn/history.rst` | The step-4 milestone changelog entry. |
| `docs/theory/foundations/operator_tensor_network.rst`, `cross_section_data.rst`, `boundary_conditions.rst`, `diffusion_1d.rst`, `api/{numerics,transport,homogeneous}.rst` | Kernel-home + binding re-points; N2N added to the boundary-zero table. |
| `docs/theory/verification/matrix.rst` | GENERATED — regenerated by the build. |

## Acceptance evidence

| gate | pre | post |
|---|---|---|
| `sphinx -E -W` EXIT | 0 | **0** |
| `WARNING:`/`ERROR:`/`CRITICAL:`/`SyntaxWarning` | 0 | **0** |
| `directives: wrote N edges` | 412 | **415** (= exactly my +3) |
| `no-implementation` statements | 18 | **18** |
| `check_docstring_xrefs.py` (stock) | 0 dead | **0 dead / 13751 decidable** |
| patched gate (L-062 `head_role` fix) + positive control | control planted → 2 dead found | **0 dead** with control removed |
| `nexus dead_references` | 0/52 | **0 dead / 52 checked** |
| matrix documented-sentinel labels | 571 | **574** (= my +3) |
| `tests/test_vv_harness_audit.py` + `test_pending_ports.py` | — | **16 passed, 5 xfailed** |

## Refuted candidates (checked, left alone, with the reason)

- `docs/theory/verification/error_catalog.rst` (4 fission hits) — all inside
  dated, past-tense `.. error-entry::` bodies (ERR bodies describe a bug *as
  it was*). The retirement rule keeps past-tense history.
- `docs/theory/foundations/coupled_block_operator.rst` — its
  `RadialCharacteristicEmission` / `coupled-ba-emission` prose is about the
  **scattering** K_iso emission, not fission; already correct post-step-3.
- `adjoint.rst`'s "the kernel-generic `RadialCharacteristicEmission` carrying
  the fission kernel" — verified against `sn/solver.py:2946`:
  `RadialCharacteristicEmission(F.kernel, …)` still, and `F.kernel` delegates,
  so the sentence is exactly true. No edit.
- `frame.rst` §`sn-homogenization-chi-collapse` — the SPATIAL χ collapse
  (production-weighted convex average). A different operation from G-F1's
  energy condensation; already correct and already cross-linked.
- `docs/theory/foundations/operator_algebra.rst:379/1538/2192`,
  `indexing_and_layout.rst:57`, `boundary_conditions.rst:2955` — name
  `FissionOperator` as an operator LEAF / an angular-flux-shaped `apply`.
  The angular binding still exists and still has those properties.
- `indexing_and_layout.rst:446/454` — `PR-INDEX-*` changelog rows, past tense.
- `sig_p` in docs (6 hits, 3 files) — all `Mixture.sig_p` (material data), not
  the retired `FissionOperator.sig_p` property. `FissionOperator.chi`/`.sig_p`
  had **0 doc sites**; the retirement was code-side only.

## REPORTED — code-side, not edited by me

1. **`orpheus/transport/operators/n2n.py` module docstring understates its
   own result.** It says the pre-harmonization spelling vs the product
   reversal is *"a pure IEEE-754 order change, principled-equivalent, gated
   at tolerance."* `[M]` 2026-08-31, 1000 draws × `gauss_legendre` n = 2/4/6/
   8/16: `np.array_equal` **1000/1000**, `max|Δ| = 0`. At ℓ = 0 the two outer
   factors degenerate (`R₀ᵀ` = plain ordinate sum, `M₀ᵀ` = per-ordinate ×wₙ),
   so the chain performs the same operations in the same order — there is no
   re-association to make. Not false, but weaker than the truth, and a reader
   could relax a gate on the strength of it. `fission.py`'s equivalent
   sentence IS correct (`[M]` ≤5 ULP, 0/200 bit-equal on 3 rules) — the two
   channels differ in kind and the docstrings currently read alike.
2. **`docs/theory/conventions/indexing_and_layout.rst:812`** claims
   `to_scalar()` / `broadcast_to_ordinates()` "live in the
   `orpheus/sn/typed_fields.py` module (planned)". `[M]` `to_scalar` = 0 hits
   in `orpheus/`; the live spellings are `AngularFlux.integrate_angular()` and
   `assemble_per_ordinate_isotropic`. Self-marked "planned", so not
   present-tense-false about existence — but it belongs to the typed-fields
   thread, not this step. DEFERRED, not fixed.
3. **37 SN-chapter sites still spell `A = L + C − S − B`** without N₂ₙ
   (census: `docs/theory/methods/sn/*.rst`, regex `L\s*\+\s*C\s*-\s*S\s*-\s*B`).
   This is **step-3** residue, not step-4's. I did NOT sweep it: on most of
   those pages Σ₂ₙ ≡ 0 by fixture and the extra term costs pedagogy, and a
   38-site algebra rewrite is a numerics adjudication that must not ride
   inside a fission docs pass. Instead the simplification is now **declared**
   at the SN chapter root (`sn/index.rst` machine header + a `.. note::`) and
   at the multigroup Key Facts, and the two sites that were genuinely
   describing the *shipped composite* were corrected. A dedicated sweep is a
   candidate follow-up issue.
