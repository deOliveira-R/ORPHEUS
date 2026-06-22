---
name: wavefront-flux-carve-lessons
description: Durable verification patterns from the WavefrontFlux storage-A (type-only cochain) + storage-B (rolling moving-frontier window) SN sweep/matvec carves (LANDED, tests live). Inherited-octant-snapshot discipline, the peak-pin Leg-1-vs-Leg-2 discriminator, nulp-einsum-vs-array_equal-typed-wrapper rule, no-perf-marker L16 diagnostic gate.
metadata:
  type: feedback
---

The WavefrontFlux carves are LANDED on `main`. Storage-A (interior
face-cochain typed as `WavefrontFlux` on `InteriorFaceSpace`) and
storage-B (rolling 2-diagonal `_MovingFrontier` window,
`O(N·ng·nx·ny)`→`O(N·ng·(nx+ny))`) both shipped; the full-field
sweep+matvec ORACLES were RECOVERED (aggressive-retirement EXCEPTION —
the fuller view pins the windowed path). Tests live at
`tests/transport/fields/test_wavefront_flux.py`,
`tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py`,
`tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py`. This note
keeps the WHY.

**1. Storage-B BIT-IDENTITY is INHERITED from the octant snapshot — and the
snapshot's own design makes it perfect for the carve.** The
`test_2d_octant_sweep_equivalence` snapshot pins `angular_flux` +
`scalar_flux` + the FOUR boundary faces (`face_xmin/xmax/ymin/ymax`) and
EXPLICITLY states interior edges are ephemeral ("neither persisted nor
compared"). Storage-B recycles interior slots → interior edges change
representation but the converged outputs + boundary trace are bit-id-preserved
→ the inherited snapshot already gates exactly the right surface. **Tolerance
discipline:** the inherited snapshot keeps `nulp=64` (its drift budget is about
the legacy einsum reduction-tree, ORTHOGONAL to the backing swap); the NEW
de-risk + peak-correctness leg demand `np.array_equal` rolling-vs-full
SAME-RUN (a slot copy IS bit-exact). General rule: **inherited einsum-bearing
snapshots stay nulp(reduction); NEW typed-wrapper pins that wrap the SAME
buffer demand array_equal.**

**2. The peak-memory PIN: Leg-1 (buffer-size ratio) is the sharp gate, Leg-2
(tracemalloc) is diagnostic.** To PROVE a memory win, the load-bearing gate is
Leg-1: the direct buffer `.values.size` ratio `size(48cells)/size(16cells)<4.0`
(rolling ≈ 3.0 linear vs full ≈ 9.0 quadratic) — DETERMINISTIC, PRINCIPLED,
SHARPEST, a foundation gate. Leg-2 end-to-end tracemalloc (`P48 < 0.6×full`) is
LOOSER and carries an HONEST caveat: the `angular_flux (N,ng,nx,ny)` OUTPUT
floor stays full/out-of-scope, so the end-to-end peak drops by the interior
chunk ONLY. Use Leg-1 as the unambiguous CI gate; Leg-2 is diagnostic-if-flaky.
The microbench discriminator that distinguishes storage-A (representation-only,
NO real memory win) from storage-B (real O(nx+ny) memory): tracemalloc peak
ratios across nx=ny∈{16,32,48} — full_interior is clearly quadratic
(~3.79×/2.21× step), rolling-2diag is linear.

**3. SCOPE CORRECTION that recurs: the 2-D MATVEC binds the interior, the 1-D
matvec does NOT.** Storage-A's plan WRONGLY claimed a "matvec mirror" — the
1-D/2-D distinction is load-bearing. The 1-D 2-D-Cartesian matvec
(`_apply_2d_cartesian`) WALKS the interior via `wavefront.face(0)/face(1)`
(NOT edge_view-only — it uses `edge_view` for OUTPUT but the full interior
`.copy()` for the WALK), so storage-B MUST cover the matvec interior backing
too (shared `WavefrontFlux.face()` → one backing swap covers both sweep and
matvec by coding-elegance Pattern 2). The 1-D matvec is a parallel-prefix SCAN
with no interior face buffer → carve is sweep-only there. ALWAYS read HEAD
source to confirm which paths actually bind the interior before scoping —
grep for the buffer accessor, not the plan's line numbers.

**4. L16 perf gate WITHOUT a perf marker → a diagnostic script, not a CI gate.**
There is no `perf` pytest marker in this project (`pyproject` confirmed). So
the L16 perf gate is a `derivations/diagnostics/diag_*_microbench.py` run
per-commit + at close-out, with a STOP circuit-breaker at >2× wall-clock and a
per-line median <5% regression threshold (a zero-copy view ⇒ <1%; a per-cell
Python fold = 10-20× = instant trip). The frontier index map MUST be built at
MESH time, not SWEEP time (per-sweep Python index construction is the trip).
The peak-memory correctness leg is the foundation test (Leg-1 above); the
wall-clock is the diagnostic.

**5. The highest-risk NEW logic in storage-B is the boundary SHED/CAPTURE.**
When the rolling frontier passes a domain edge, its outflow must be CAPTURED
into the BoundaryFlux at the EXACT level the frontier passes the boundary. A
recycled-overwrite that corrupts the boundary trace is SILENT — the interior
looks fine, but the reflective eigenvalue drifts under iteration (Mode-3
fingerprint). The de-risk STOP-gate for it: rolling≡full `array_equal` on
angular + scalar + all 4 faces across vacuum + reflective + 2-sweeps +
multi-octant + HETEROGENEOUS-2-zone (H2), BEFORE wiring.

See [[snapshot-migration-when-production-goes-bare]] (the shared-driver SoT
recipe the octant snapshot uses), [[regression-tolerance-design]] (array_equal vs nulp),
[[phase5a-5c-angular-windowing-lessons]] (the orthogonal angular-windowing
carve that drops the OUTPUT, vs storage-B which keeps outputs bit-id).
