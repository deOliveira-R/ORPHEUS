---
name: issue-257-s3b-multiplication-operator-closeout
description: #257 S3b — CollisionOperator → thin subclass of the new transport MultiplicationOperator (the §5.7 promotion C = M[σ_t]); the honest CAP_SOLVE-iff-min|σ|>0 spectrum gate; 3 prod sites rewired to the typed CrossSectionField accessor; principled-equiv (0 ULP) + multiplier-algebra law-suite
metadata:
  type: project
---

# #257 S3b — MultiplicationOperator promotion (C = M[σ_t]) CLOSEOUT

`feature/field-typed-operator-algebra` (off `main` @ `05fa1ef`), built on
HEAD `c1da42d` (S3a). NOT committed (main agent reviews [elegance + qa] +
commits). The transport half of the §5.7 operator-as-promotion fold.

## What was built (4 deliverables, all green)

1. **NEW `orpheus/transport/multiplication_operator.py`** — `MultiplicationOperator(LinearOperatorMixin["TimedFullField"])`,
   FIRST operator in `transport/`. The §5.7 promotion `M[f]` of a
   `CrossSectionField`. Stores ONLY `coefficient: CrossSectionField`
   (MESH-FREE — `mesh = psi.bulk.mesh` at apply time). `@dataclass(eq=False)`;
   `capabilities` computed in `__post_init__` by BUILDING the S3a
   `DiagonalOperator(coefficient.values, broadcast_axes=(0,))` engine and
   INHERITING its capability set (single-sources the spectrum law — the
   transport operator and the numerics engine agree by construction, NOT
   a copied predicate). `_engine()` helper (Pattern 2) is the ONE site the
   broadcast lives. `apply` → bulk `AngularSourceSink.from_mesh` (a
   SOURCE; collision turns flux into a rate) + `BoundarySourceSink.zeros_on`;
   `solve` → `AngularFlux.from_mesh` + `BoundaryFlux.zeros_on`;
   `apply_transpose=apply` (self-adjoint, real coeff); `block_role=BlockRole.BULK`.
   Generic `__add__` inherited from the mixin (does NOT reference
   `InvertibleOperator`/`StreamingOperator` — those are L3; transport L2
   cannot import sn). pyright 0/0.

2. **`CollisionOperator` → thin SUBCLASS** (`orpheus/sn/operator.py:510`):
   `class CollisionOperator(MultiplicationOperator)`. NOT a dataclass —
   explicit `__init__(sn_mesh, sigma)` with the back-compat conversion:
   ndarray → `CrossSectionField.from_mesh(np.asarray(sigma), sn_mesh)`
   (the SAME factory the S2 `total_cross_section_field` accessor uses, so
   the two paths produce structurally identical coefficients);
   CrossSectionField → passed straight to base. `self.sn_mesh = sn_mesh`
   (the SAME object — the `InvertibleOperator` `streaming.sn_mesh is
   diagonal.sn_mesh` identity check at `:797` depends on it). `.sigma` is
   a `@property` returning `coefficient.values` (Pattern 2, no duplicate
   storage — `InvertibleOperator` reads `self.diagonal.sigma` at `:854`
   and validates `np.all(sigma > 0)` at `:803`). Override `__add__` keeps
   the `if isinstance(other, StreamingOperator): return
   InvertibleOperator(other, self)` dispatch, else `super().__add__` (now
   the multiplier/mixin generic). apply/solve/apply_transpose INHERITED
   (the duplicated bodies at the old `:592/:620/:651` DELETED — single-source
   the multiply).

3. **3 prod sites rewired** to `total_cross_section_field` (the S2 typed
   accessor): `solver.py:217` (`make_within_group_operators`), `:925`
   (`__init__` `L = Streaming + Collision`), `:1000` (`rebind_sigma_t`).
   The 87 test callers stay on the ndarray back-compat path.

4. **NEW `tests/transport/test_multiplication_operator.py`** (11 tests,
   all `@foundation`, all -O-firing): broadcast oracle (engine ≡ legacy
   `σ[None]·ψ` at VALUES level, `assert_array_equal` PRIMARY, 0 ULP — same
   op) on the discriminating **2-D nx=5≠ny=3, ng=2** carrier; the
   multiplier-algebra law-suite (`M_1=I`, `M_0=ZeroOperator` codomain-aware,
   linearity `M_{af+bg}=aM_f+bM_g` on ≥2G ASYMMETRIC HET, self-adjoint
   `M.H=M`, spectrum `CAP_SOLVE iff min|f|>0`, homomorphism `M_f∘M_g=M_{f·g}`
   at the VALUES level via the engine on the raw `σ·σ'` product array —
   units cm⁻², NOT a CrossSectionField); L0 collision-only streaming-
   equilibrium (`ψ=Q/σ_t`).

## The CAP_SOLVE behavioral-change audit (the brief's ⚠)

GREP'd every production + test path constructing a CollisionOperator and
calling `.solve` / relying on `CAP_SOLVE`. **NOTHING breaks:**
- All 3 PROD sites use `total_cross_section` (σ_t > 0, bounded away from 0).
- The σ_r `CollisionOperator(σ_r)` + `.solve` is ONLY in an `operator.py:720`
  DOCSTRING — NO production caller.
- `InvertibleOperator.__init__` has its OWN stricter construction-time
  `np.all(diagonal.sigma > 0)` check (`:803`) — UNCHANGED, consistent with
  the new gate.
- The ONE test that builds a CollisionOperator with a ZERO-entry σ
  (`test_collision_operator.py::TestSigmaLayout::test_localised_sigma_localised_output`,
  `sigma[:, ix_target]=1.0` else 0) calls ONLY `.apply` — never `.solve`,
  never checks `CAP_SOLVE`. Construction still succeeds (the gate only
  revokes CAP_SOLVE, never blocks construction); apply still works. GREEN.
- `TestCapabilities::test_all_three_capabilities` uses `_sigma_total`
  (0.3 + 0.5·rand, strictly positive) → all 3 caps → GREEN.

The behavioral change is purely ADDITIVE-honesty: σ=0 now revokes
CAP_SOLVE (was silent IEEE NaN). Pattern 4.

## Verification (principled-equivalence gate, user 2026-06-19)

- **Broadcast**: 0 ULP — the generalized engine `expand_dims(coeff,(0,))·x`
  on σ_t IS `sigma[None]·x` (reduction_depth=1, pure broadcast-multiply, no
  sum). `assert_array_equal` holds (the oracle test). Not force-fit — it is
  literally the same op.
- **Structurally-independent refs (kept green)**: `test_kinf_homogeneous.py`
  (≥2G analytical k_∞=νΣf/Σa, eigenvalue pillar, routes σ_t through the
  promoted C in the loss operator) + `test_invertible_operator.py::...test_si_carve_recovers_analytical_kinf`
  + `..._recovers_q_over_sigma_composite` (streaming-equilibrium L+C
  resolvent). All 35+ green.
- **pyright**: 2295 errors / 19 warnings — IDENTICAL to baseline (ZERO
  net-new, ZERO new `# type: ignore`). Both new files 0/0 standalone.
- **Regression subset** (`tests/sn/operators spatial sweep/core solve
  numerics transport`, -O, deselect het keff): 7 failed / 1968 passed.
  The 7 = EXACTLY the documented baseline reds (#250 SPHERE ×5 + #232
  mu_y ×2). +11 = the new tests. NO non-baseline regression.
- **Sphinx**: clean (only pre-existing `mesh.py` paramref error + unrelated
  test SyntaxWarnings). NEW `:label: multiplication-operator-promotion` +
  `:eq: multiplication-operator-action` + `:mod:` cross-refs + archivist
  `.. todo::` in `docs/theory/operator_algebra.rst` — present in built HTML.

## Mutation-verification (the tests are CONSTRAINED, not just green)

- WRONG `broadcast_axes=(1,)` on the nx≠ny carrier → ValueError (RAISES) —
  the broadcast oracle catches an axis-ordering bug. A square mesh would
  silently agree → why nx≠ny+ng=2 is the discriminator (variable-swap mode #2).
- LEGACY always-on gate (monkeypatch `__post_init__` → always CAP_SOLVE) →
  the σ=0 coefficient WOULD advertise CAP_SOLVE → `test_spectrum_cap_solve_iff_min_abs_positive`
  REDS. The honest gate is genuinely constrained.

## NO algebra-of-record manifest (per the plan + prior-closeout precedent)

This is an operator-algebra CARVE (a coefficient→operator promotion + a
numerics-primitive delegation), NOT a new L0/L1 reference solver. There is
NO Branch-1 SymPy / L1 cross-check / capability-matrix owed — the
verification is the multiplier-algebra law-suite (intrinsic properties) +
the EXISTING structurally-independent kinf / streaming-equilibrium refs.
Same posture as #257 S3a, #247, #251 (numerics-primitive widening). NO new
ERR (no real bug caught — the σ=0 gate is a PROACTIVE Pattern-4 hardening,
not a found defect; next free ERR-063 still reserved per #251 note).
archivist DISPATCH emitted for the rich operator_algebra.rst narrative
(`followup: false`).

## LESSON (coding-elegance + the subclass-promotion pattern)

Promoting an existing named leaf onto a new general base via a THIN
subclass: (1) inherit the action (apply/solve/apply_transpose) from the
base — DELETE the duplicated bodies, the multiply lives once; (2) the
subclass adds ONLY the name + the back-compat constructor + the
domain-specific algebra dispatch (`__add__`→InvertibleOperator); (3) keep
the legacy attribute surface alive as PROPERTIES off the base's field
(`.sigma` → `coefficient.values`) so downstream readers (`InvertibleOperator.diagonal.sigma`)
are byte-unchanged with NO duplicate storage (Pattern 2). (4) The honest
capability gate is single-sourced by BUILDING the engine and inheriting its
capability set in `__post_init__`, NOT by re-deriving `min|f|>0` — the
transport operator and the numerics engine then agree by construction.
(5) A non-dataclass subclass of a dataclass base sets the base's field via
`super().__init__(coefficient=...)` after the ndarray↔field conversion —
the base's `__post_init__` then runs the spectrum gate. The
numerics↛transport layer boundary FORBIDS folding the two classes (the
`L+C→InvertibleOperator` L3 dispatch MUST stay on the sn subclass; the
generic `__add__` stays on the L2 base).
