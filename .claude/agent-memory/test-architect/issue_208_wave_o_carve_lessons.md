---
name: issue-208-wave-o-carve-lessons
description: Durable verification patterns from the #208/#201 Wave-O operator-typing + affine-flux-algebra carve (LANDED, tests live). Snapshot-on-random-ψ + two-anchor discipline, defect-vs-raw zero-input precision boundary, sha256-golden bit-id gate, affine-torsor algebra gate, residual-before-adjoint sequencing.
metadata:
  type: feedback
---

The #208/#201 Wave-O carve (operator block-role typing → SOURCE-side
dimensional retyping → typed `from_balance` residual → affine
`FluxDisplacement` torsor + flux-add gate) is LANDED on `main`
(`8c2f355`→`63719a2`; tests live: `tests/sn/operators/test_bc_extraction_matvec.py`,
`tests/sn/operators/test_typed_residual_evaluation.py`,
`tests/transport/residuals/test_typed_residuals.py`,
`tests/numerics/test_affine_flux_algebra.py`,
`tests/sn/solve/test_flux_displacement_diagnostics.py`,
`tests/sn/solve/test_affine_carve_bit_identity.py`,
`tests/sn/operators/test_g_adjoint_reciprocity.py`). The test files are
the WHAT. This note keeps the WHY that the test files do not encode.

**1. Snapshot-on-RANDOM-ψ for the bulk + a SEPARATE structurally-independent
VALUE anchor — never one alone.** The curvilinear M-M Carlson bulk matvec
has no cheap closed form, and flat ψ NULLS redistribution (vv §H2), so the
ONLY ψ that stresses the operator is a fixed-seed RANDOM ψ. Gate the bulk
with a committed `.npy` snapshot of `(L+C).apply(random_ψ)` (proves "didn't
move" = bit-identity inheritance) AND a `Q/Σ_t` flat-flux closed-form value
(proves "is correct"). vv §bit-identity criterion 2: ULP-distance is
necessary-NOT-sufficient — the snapshot alone cannot tell you the pre-carve
value was right. This two-anchor pattern is the template for any pure-refactor
matvec carve where the natural input is random.

**2. The `--capture-baseline` snapshot mechanism + its gotcha.** A
`--capture-baseline` flag makes the test WRITE the `.npy` then `pytest.skip`;
flag absent → READ + `assert_array_equal`. **GOTCHA: `pytest_addoption` ONLY
fires in a ROOT conftest (`tests/conftest.py`), never in a test module** —
move the addoption hook to the root before capturing. Capture ONCE at the
pre-carve HEAD, `git add` the `.npy` alongside the first carve commit.

**3. The defect-vs-raw ZERO-INPUT precision boundary (BC-extraction carves).**
When a carve changes a boundary slot from a DEFECT (`outflow − face_value`)
to RAW outflow: for the canonical ZERO-input vacuum test
(`BoundaryFlux.zeros_on`), `face_value = 0` so defect ≡ raw → boundary slot
is bit-identical FOR ZERO INPUT. The two only diverge on NON-zero input — so
the discriminator test MUST feed a non-zero outflow (`TestVacuumBoundaryDefectVsRaw`,
xfail-strict pre-carve, flips at the boundary-emission sub-commit). A
zero-input-only test is a latent dud for this carve class.

**4. The affine torsor algebra gate (FluxDisplacement, #208 box 7+).** The
typed algebra: `flux⊖flux→FluxDisplacement` (affine difference space V);
`flux⊕disp→flux` (torsor action); `flux⊕flux→TypeError` (RAISES — names
`affine_combination` instead); `disp⊕disp→disp`; `scalar·disp→disp`;
`affine_combination(Σλ=1)→flux` (SI relaxation IS this; Σλ≠1 RAISES — pick
ω=0.7/0.3 that sum to EXACTLY 1.0 in float64). Gate per L11 PAIRED: flux+flux
RAISES AND flux+disp WORKS through the same function; affine_combination Σλ=1
WORKS AND Σλ≠1 RAISES. The torsor round-trip `ψ₁⊕(ψ₂⊖ψ₁)==ψ₂` is legitimately
`array_equal`-EXACT (d.values IS the raw single subtraction, not a float
identity claim) — and bit-EXACT here CATCHES a re-projecting mint (a
`FluxDisplacement` that round-trips through a projection drifts off bit-id).
**R-A LOAD-BEARING:** the production SI iterate is a composite `TimedFullField`
(bulk=AngularFlux|HarmonicMomentField, boundary=BoundaryFlux), NOT a bare
AngularFlux — the retyped `__sub__` must compose UP through the composite, so
add a composite-level torsor test (the 5th surface beyond the 4 leaves). The
stopping norm STAYS the flat `_l2_norm(disp.to_flat())`, NOT the metric
`Field.l2` (switching re-interprets conv_tol AND breaks bit-id).

**5. sha256-golden bit-id gate when the existing snapshot gate is too loose.**
The DD regression snapshots gate at `SAFETY×conv_tol≈1e-11`, and the het-aniso
case ALREADY pre-drifts ~6920 ULP from prior carves — INSUFFICIENT to prove
ZERO numerical change for a pure type-relabel. The correct gate is a NEW
capture-before/compare-after sha256-golden (`sha256(psi.to_flat().tobytes())`
+ `repr(keff)==`), generated at the pre-carve HEAD, covering ≥2G+het SI-slab
AND a Krylov-2D-2G-het-aniso case (so the HarmonicMomentField-bulk windowed
iterate is covered). Ship the golden FIRST, with the first carve piece.

**6. Sequencing: residual must exist BEFORE adjoint; build the adjoint ONCE on
the FINAL shape.** The 2-D boundary-residual ADD (O.4b) precedes the adjoint
(O.2) precisely because the residual must be a typed quantity before `A†` can
be a typed defect — do not build the adjoint against an intermediate shape and
re-do it. The −B sign is NOT matvec-observable (it composes at the OperatorSum
residual level) → the reflective sign-flip manifests as MMS divergence at the
SOLVER level, not the matvec level; the matvec-level catch is structural-only
(per-face B-action equivalence pins B's action, not its sign in the residual).

The g-adjoint metric note: the existing Phase-C reciprocity gate was
metric-BLIND on TWO axes (`np.sum(Lpsi*phi.bulk.values)` reads ONLY the bulk,
plain-Euclidean — drops the volume measure V AND the trace block). A
sign-flipped B passed it. The fix rewrites the body to read BOTH blocks via
`op.codomain.inner_product`, with a MANDATORY L11 wrong-metric control
(drop-cosine / drop-weight must break reciprocity by rel>1e-3). White-BC is
the discriminating case (self-adjoint ONLY under the cosine metric; reflective
is a Euclidean involution already). This candidate "metric-blind reciprocity"
mode is NOT skill-appended (abstract until a real B sign-flip escapes a
bulk-Euclidean gate → then it earns an ERR-NNN).

See [[regression-tolerance-design]] (SAFETY×conv_tol vs nulp vs sha256-golden
ladder), [[phase4-46-nonvacuum-mms-ansatz]] (the prescribed-inflow MMS that
4.6 built on the FINAL forward shape). The b5 "dimensional sin" carve
(`IterationResidual.from_balance` + boundary-output retype) was the precursor
that evolved INTO this typed-residual landing — superseded, no separate note.
