---
name: s5-2-scanmarch-end-to-end-gates
description: S5.2 #222 ScanMarch end-to-end FP-invariance/eigenvalue gates — PASS clean; the context-manager-not-fixture forcing pattern + the independent-witness config-duplication exception to Pattern 2.
metadata:
  type: project
---

`tests/sn/solve/test_scan_march_end_to_end.py` (S5.2, #222) — **PASS clean, no conditions.**
The deferred end-to-end gates that drive FULL production solvers with ScanMarch FORCED, after
the sweep/matvec-layer-only G2 nulp oracle. Four tests (G4.a aniso+het+vacuum, G4.b all-reflective
ERR-056 shed pin, G6 k_inf anchor + SI≡Krylov het). Two reusable rulings:

**1. The forcing MUST be a context manager, not a pytest fixture — and this is an ELEGANCE win, not
a mechanical choice.** Each FP-invariance pair needs the REFERENCE leg UNFORCED; a test-scoped
fixture forces BOTH legs and collapses the comparison to scan-march≡scan-march (a vacuous tautology
that PASSES while testing nothing). The `with`-block boundary is the type-level guarantee the two
legs use different schedules. This exact self-comparison bug was caught at design time. When
reviewing any "force representation X and compare to the default" gate, the discriminator is: does
the reference leg run OUTSIDE the forcing scope? A fixture cannot give you that; a context manager can.

**2. The one-patch-forces-everything claim is ENFORCED by S6.5's single-instance collapse — verified
at source, not assumed.** Post-S6.5 `patch.object(lr, "default_for", forced)` reaches all three doors
because each resolves the module attr at call time: (a) `transport_sweep` calls module-local
`default_for` (`loss_representation.py:1613`); (b) the operator cached_property does
`from .loss_representation import default_for` INSIDE the method body (`operator.py:1547`) so it
re-reads the patched attr per construction; (c) the G-S resolvent + Krylov matvec reach it via
`self.streaming.loss_representation` (`operator.py:1948`) — the SAME streaming-leaf cached_property,
no independent `default_for` call. Latent hazard (cached_property escaping the patch if the operator
were built before the `with`) does NOT fire because every test builds the whole solve inside the
`with`. See [[s6-5-one-representation-instance]] — this gate is the consumer-side proof that the
construction-inventory reduction actually single-points the forcing.

**3. NON-VACUITY spy is immune to window interception — the structural check that makes it sound.**
`ScanMarch._sweep_interior` (`loss_representation.py:1232`) is defined ON `ScanMarch`, DISTINCT from
`MovingFrontierWindow._sweep_interior` (`:844`, via `_DAGWavefront`). Patching `ScanMarch._sweep_interior`
therefore CANNOT intercept a window run; if forcing silently fell back, the counter stays 0 and
`_require_ran` fails via `pytest.fail`. When reviewing a "force-and-count" non-vacuity tripwire, ALWAYS
check the spied method is not shared with the fallback path's class — a shared base method makes the
counter fire on the wrong path and produces a false-green.

**4. RULING — independent-witness config duplication is a PRINCIPLED EXCEPTION to Pattern 2, not a
violation.** `_build_2d_aniso_het_vacuum` is a byte-for-byte mirror (NOT import) of the affine golden's
private `_build_2d` (`test_affine_carve_bit_identity.py:84`). Demanding a shared builder is the WRONG
move here: the golden's builder is module-private and sha256 hash-pinned; coupling the FP-invariance
gate to it via import would let a future config edit SILENTLY MOVE BOTH the bit-identity reference and
the FP config simultaneously — masking the golden's drift instead of catching it. The duplication is a
deliberate independent witness; the docstring cross-reference ("mirrors the affine golden") is the
correct discipline. GENERAL RULE: when two tests pin the SAME physical config and one of them is a
FROZEN/hash-pinned reference, independent transcription is correct — sharing converts an observable
divergence into a silent co-movement. Pattern-2 "extract duplication" assumes the copies should track
together; a frozen-reference witness is the case where they must NOT.

**5. The G4.b all-reflective config-fix is the anti-pattern-#17 textbook response.** Uniform source
under all-reflective measured max/min=1.067 → tripped `_assert_nonflat` (`:178`, the non-flat guard).
The fix CONFINED the source to the fuel half (`q_ext[:,:,:4,:]`) to drive a real x-gradient — did NOT
loosen the `ratio<=1.2` threshold. Guard worked, config changed, WHY documented (`:144-150`). The
honest d=2 LIMITATION (full diagonal shared-face ERR-056 stressor needs d=3, no 3-D quad yet) is stated
not papered over. `scalar_flux.values` is `(ng, nx, ny)` so `phi[0].mean(axis=1)` correctly measures the
x-profile gradient — the axis the heterogeneity and half-domain source vary along.

Solver-tol-not-nulp comparison (`rtol=1e-6` vs `inner_tol=1e-12`, ~6 orders headroom) is the correct
schedules-differ-at-FP-association boundary. Bare asserts for numerical claims (rewritten, live under -O
per L26); `pytest.fail` for tripwires.
