# D-I.2 — test_scattering_operator.py bare-ndarray → typed AngularFlux migration

`refactor/sn-operator-algebra` 2026-05-29. **TESTS-ONLY contribution**;
production operator surgery deferred to the main agent per the brief's
explicit boundary.

## What migrated (count + sample tests)

15 call-site migrations in `tests/sn/test_scattering_operator.py`:

- **`op.apply(psi)` sites — 11 sites** across 5 test classes:
  - `TestProtocolCompliance::test_apply_accepts_psi_shape` — psi now
    `AngularFlux.from_mesh`; assertion on `out.values.shape`.
  - `TestApplySemantics::test_apply_isotropic_flux_p0_only` — typed psi
    + expected ndarray remains, compares `actual.values` to expected.
  - `TestApplySemantics::test_apply_zero_psi_returns_zero` — typed
    zeros AngularFlux → expects `out.values == zeros_like`.
  - `TestApplySemantics::test_apply_linearity` — `alpha * psi1 + beta
    * psi2` via Field dunders (typed AngularFlux); RHS also rebuilt
    via PerOrdinateSource dunders.
  - `TestProducerSideNormalisation::test_typed_apply_returns_per_ordinate_already_normalised`
    — typed `c` constant flux; assertion on `actual.values`.
  - `TestAlgebraicIdentity::_check_identity` helper rewritten +
    5 dependent tests (4 random/uniform/Pℓ≥1/n2n/cross-group identity
    + 1 residual-zero — synthetic standalone solver).
- **`op.build_aniso_source(psi)` sites — 4 sites** in
  `TestAnisotropicScatteringExtraction`:
  - `test_returns_none_for_p0` — typed AngularFlux still returns None.
  - `test_isotropic_flux_zero_aniso_source` — typed psi; assertion on
    `Q_aniso.values`.
  - `test_delegator_matches_operator` — KEY case from the brief.
    Delegator gets bare ndarray (preserves solver.py:1203 contract);
    operator gets typed AngularFlux. Comparison: `out_via_delegator`
    (ndarray) vs `out_via_operator.values` (PerOrdinateSource→ndarray).
  - `test_returns_none_for_no_angular_flux` — `None` short-circuit; no
    migration needed.

Module-level imports added: `AngularFlux`, `PerOrdinateSource` from
`orpheus.transport.fields.angular_flux` / `orpheus.transport.sources`.

## What didn't need migrating (already typed)

- **`tests/sn/test_harmonic_moment_field.py`** — verified at line 454:
  `psi = AngularFlux.from_mesh(psi_values, sn_mesh)` already typed.
  31/31 tests PASS. NO edits made to this file.
- **`TestCompositeInvariants` (4 tests in test_scattering_operator.py)** —
  consume `state = sn_mesh.zeros_timed_full_field()` (TimedFullField).
  Already on the typed dispatch path.
- **`TestBitIdenticalExtractionP0` (5 tests)** — exercise
  `add_iso_source` / `add_n2n_source` helpers that accept bare ndarray
  by design (per-cell P0 + (n,2n) accumulator API). NOT in scope.
- **`TestP0AlgebraicIdentities`** — same; uses `add_iso_source` /
  `add_n2n_source` directly.
- **`TestFoldablePart`, `TestResidualPart`, `TestFoldableSigma`,
  `TestPurity`, `TestIsFoldableIntoSigmaR`** — structural predicates
  on the operator data; no flux argument involved.

## What broke + status

**Surfaced an unstated brief precondition**: the brief states
"AngularFlux arm — find it (likely registered for AngularFlux typed
input; returns PerOrdinateSource)". The AngularFlux `apply.register`
arm DOES NOT EXIST in `orpheus/sn/scattering.py` at HEAD. It was
deleted in commit `1aa6efe` (D-H.2-C3, "retire legacy operator
branches"). Currently registered arms: `TimedFullField`, `ScalarFlux`,
`np.ndarray` (the one being retired in D-I.2). The base `apply`
raises `TypeError("unsupported input type")`.

**Consequence**: 11 of the 91 tests now FAIL after migration with:

```
TypeError: ScatteringOperator.apply: unsupported input type AngularFlux;
expected AngularFlux, ScalarFlux, or numpy.ndarray. Dispatch table is
registered via @singledispatchmethod.
```

This is precisely the contract for the main agent's surgery. After the
main agent ADDs the typed `AngularFlux` arm (returning
`PerOrdinateSource`) and RETIREs the `np.ndarray` arm per the D-I.2
plan, all 11 will turn green by construction.

The `build_aniso_source(AngularFlux)` migrations (4 tests) PASS
unchanged because that method already has a typed `is_typed` branch
returning `PerOrdinateSource` (scattering.py:626).

**No fixes applied by me** — production code is out of scope per
brief's "Do NOT touch orpheus/sn/scattering.py — main agent owns the
operator surgery."

## Test count audit

| State          | Total | PASS | FAIL | Notes                                |
|----------------|-------|------|------|--------------------------------------|
| Pre-migration  | 91    | 91   | 0    | (60 + 31) baseline                   |
| Post-migration | 91    | 80   | 11   | Same cardinality; no xfail; no skip  |

After main-agent surgery: expected 91 PASS / 0 FAIL.

## Verbatim verification (per L12 paste-back)

```
$ /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python -O -m pytest \
      tests/sn/test_scattering_operator.py tests/sn/test_harmonic_moment_field.py \
      --tb=no -v

(tail of output)
======================== short test summary info ========================
FAILED tests/sn/test_scattering_operator.py::TestProtocolCompliance::test_apply_accepts_psi_shape
FAILED tests/sn/test_scattering_operator.py::TestApplySemantics::test_apply_isotropic_flux_p0_only
FAILED tests/sn/test_scattering_operator.py::TestApplySemantics::test_apply_zero_psi_returns_zero
FAILED tests/sn/test_scattering_operator.py::TestApplySemantics::test_apply_linearity
FAILED tests/sn/test_scattering_operator.py::TestProducerSideNormalisation::test_typed_apply_returns_per_ordinate_already_normalised
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_identity_p0_only_random_psi
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_identity_p0_only_uniform_psi
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_identity_with_pl_ge_1
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_identity_with_nonzero_n2n
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_identity_multigroup_cross_group_plus_diagonal
FAILED tests/sn/test_scattering_operator.py::TestAlgebraicIdentity::test_residual_zero_when_p0_diagonal_only_no_n2n
================== 11 failed, 80 passed, 3 warnings in 0.69s ===================
```

All 11 failures share the identical root cause: `TypeError:
ScatteringOperator.apply: unsupported input type AngularFlux`.

## Audit: bare-ndarray sites in the migrated file

```
$ grep -n "op.apply\|op.build_aniso_source\|scattering_op.apply\|
  scattering_op.build_aniso_source" tests/sn/test_scattering_operator.py
```

All remaining sites consume either:
- A typed `AngularFlux` constructed via `AngularFlux.from_mesh`
- A typed `TimedFullField` from `sn_mesh.zeros_timed_full_field()`
- `None` (short-circuit before dispatch)
- Bare ndarray INTO the `_build_aniso_scattering` solver delegator only
  (preserves the bare-ndarray external contract per brief)

Zero bare-ndarray `op.apply(...)` or `op.build_aniso_source(...)`
remain. Migration is complete on the test-side.

## Lessons applied

- **L12 (closeout-time substitution)** — verbatim pytest stdout pasted;
  test cardinality (91 = 60 + 31) confirmed against pre-migration
  baseline; no xfail introduced; no tests deleted.
- **L13 (name existing types to extend)** — consumed `AngularFlux` +
  `PerOrdinateSource` from their canonical modules
  (`orpheus.transport.fields.angular_flux`,
  `orpheus.transport.sources`) without reinvention. Used `Field`
  dunders (`+`, scalar `*`) for the linearity test instead of
  pre-combining ndarrays.
- **vv-principles "Bit-identity vs principled-equivalence"** —
  preserved rtol contracts (1e-13 for apply/aniso, 1e-14 for
  algebraic-identity FP-non-associativity gate, 1e-15 atol for
  residual-zero). NO tolerance relaxation.
- **vv-principles anti-pattern flag** — the brief assumed an
  unverified precondition ("AngularFlux arm — find it"); the audit
  surfaced its absence. Logged here so the main agent can scope the
  surgery correctly: D-I.2 must ADD the AngularFlux arm AND retire
  the ndarray arm in the same atomic commit.

## Risk callouts surfaced

1. **`test_delegator_matches_operator` (line 285)** — chose the
   "migrate, don't delete" branch per brief recommendation. Verified
   the delegator still returns ndarray (existing
   `_build_aniso_scattering` at solver.py contract); operator returns
   `PerOrdinateSource`. The bit-identity gate is now between
   `delegator_output` (ndarray) and `operator_output.values` (ndarray
   unwrap). If the main agent inlines the delegator at solver.py:1203
   later, this test becomes stale; for now it pins the bit-identical
   math.

2. **`AngularFlux + AngularFlux` arithmetic** — verified via
   `orpheus/numerics/field.py:201,205,212,215` (`__add__`, `__sub__`,
   `__mul__`, `__rmul__` on the base Field). The linearity test's
   `alpha * psi1 + beta * psi2` is valid AngularFlux arithmetic via
   inherited dunders. Used the same pattern for the RHS via
   PerOrdinateSource dunders.

3. **`PerOrdinateSource` consumption** — `out.values` for ndarray
   comparisons throughout. Output type is `PerOrdinateSource` (NOT
   AngularFlux) per the brief's "returns PerOrdinateSource" contract.

## Next steps (for the main agent)

1. Add `@apply.register` arm in `orpheus/sn/scattering.py` for
   `AngularFlux` input → `PerOrdinateSource` output (the body
   mirrors the existing `TimedFullField` arm minus the bulk/boundary
   composite wrapping).
2. Retire the `np.ndarray` `@apply.register` arm at
   scattering.py:976-1006 per D-I.2 brief.
3. Retire the `is_typed` ndarray fallthrough at scattering.py:618-671
   inside `build_aniso_source`.
4. Migrate `_build_aniso_scattering` delegator helper at
   `orpheus/sn/solver.py:_build_aniso_scattering` + caller at
   solver.py:1203 to wrap `angular` as `AngularFlux` internally.
5. Re-run `tests/sn/test_scattering_operator.py` + `tests/sn/test_harmonic_moment_field.py`
   — expected outcome: 91 PASS / 0 FAIL.
