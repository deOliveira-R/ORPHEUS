# Issue #236 Phase 2 (Step A) — τ producer relocation onto the angular closure

**Branch context:** `feature/sn-spatial-angular-product` (the #236 campaign).
**Status:** DONE, NOT committed (the dispatching agent commits after review).
**Date:** 2026-06-17. HOST env, canonical `.venv/bin/python -O -m pytest`.

## One-sentence summary

`MorelMontryAngularSweep` now PRODUCES the Morel–Montry weight τ from the
quadrature `(μ, w, levels)` it already binds (an angular-scheme property),
instead of reading `reduced.tau_mm` / `reduced.tau_mm_per_level` back from the
streaming-geometry factory. BIT-IDENTICAL (0-ULP replica of the factory
arithmetic); the geometry factory still bakes an identical τ for the sweep path
(P1/P2/P3) — the parallel-run-and-compare de-risking move. Step B retires the
geometry-side producer.

## Deliverables (manifest)

1. **τ producer** — NEW module-level function
   `morel_montry_tau_per_level(quad, coord)` in
   `orpheus/sn/spatial/pole_angular_closure.py` (inserted right before the
   `MorelMontryAngularSweep` class, ~λ514). Replicates the factory line-for-line:
   - SPHERE: `reduced_operator.py:681-688` — `mu_edge[0]=-1.0`, weight-sum
     accumulation `mu_edge[n+1]=mu_edge[n]+w[n]`, `tau=(mu-mu_edge)/dmu if
     abs(dmu)>1e-15 else 0.5`. **UNCLAMPED** (W1).
   - CYLINDER: `reduced_operator.py:798-815` — per μ-level, η-midpoint edges
     with ±sinθ endpoints, `tau_raw=...`, then `max(0.5, min(1.0, tau_raw))`.
     **CLAMP RETAINED** (structural τ_raw=0 ÷0 block at the most-inward ordinate).
   - Cartesian → raises `ValueError` (Identity supplies neutral τ by TYPE).
   STRUCTURAL INDEPENDENCE: it's the closure's OWN code — does NOT call
   `contamination.morel_montry_weights` (that is the Leg-1 verification
   reference; using it in production = reference contamination = tautology).
2. **`__init__` flip** — `pole_angular_closure.py` ~λ775: the sphere/cylinder
   arms now set `self._tau_per_level = tau_per_level` (from the producer)
   instead of `(reduced.tau_mm,)` / `tuple(reduced.tau_mm_per_level)`.
   Everything downstream of `_tau_per_level` (the c_in/c_out precompute, the
   recurrence kernels) is UNCHANGED — it now consumes closure-owned τ.
3. **`IdentityAngularClosure`** — NO change. It already exposes
   `self._tau_per_level = (np.ones(self._N),)` (neutral τ=1.0) produced by its
   TYPE, no geometry branch. Confirmed by read (`pole_angular_closure.py:1289`).
4. **Leg-1 gate** — NEW `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`
   (9 tests, all `@pytest.mark.foundation`):
   - sphere closure-τ `np.array_equal` `spherical_streaming(...).tau_mm` [N=8,16]
     + the BOUND `MorelMontryAngularSweep._tau_per_level` equals it (wiring proof).
   - sphere closure-τ `np.testing.assert_array_equal`
     `morel_montry_weights(quad,"spherical")` [N=8,16] (independent ref, 0-ULP).
   - cylinder closure-τ `np.array_equal` `cylindrical_streaming(...).tau_mm_per_level`
     [n_phi=8,16] per level + bound closure (wiring proof).
   - cylinder NEGATIVE control: closure-τ == `clip(reference τ_raw, ½, 1)` AND
     differs from RAW reference where the clamp bites; guard-the-guard asserts
     `any(t.min()<0.5)` so the clamp is exercised.
   - identity: τ=1.0 + producer RAISES on Cartesian (TYPE is the dispatch).
   Mode-8 SAFE: all VALUE checks via `np.array_equal`/`np.testing`; the two
   guard predicates use bare assert (pytest rewrites those in test modules).
   MUTATION-VERIFIED: dropping the cylinder clamp → 4 reds (both the 0-ULP
   factory-equivalence AND the negative control go red), restored to green.
5. **Sphinx stub** — `docs/theory/discrete_ordinates.rst` ~λ956, label
   `.. _sn-tau-closure-owned:`, subsection "τ is an angular-scheme property",
   `.. todo:: Archivist expansion needed.` with `:func:` cross-refs to the
   producer + both factories + the contamination reference + the gate path.
   Build exit 0, NO new warnings (the inconsistent-title ERROR I first
   introduced — `^^^` skipped level 3→5 — was fixed to `~~~`, the file's
   level-3 char). Pre-existing `paramref` ERROR on `mesh.py` is unrelated.
6. **This closeout memo.**
7. **archivist DISPATCH_REQUEST** emitted (followup:false).

## OUT OF SCOPE (Step B/C — NOT touched)

- `tau_mm`/`tau_mm_per_level` NOT removed from `ReducedStreamingOperator`.
- Geometry factory τ production UNCHANGED (still feeds the sweep path P1/P2/P3).
- `cell_balance.py`, `diamond.py`, `sweep_cache.py` NOT touched.
- c_in/c_out four-site duplication NOT collapsed.

## VERIFICATION (paste-back, L12)

NEW gate:
```
tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py ... 9 passed, 1 warning in 0.23s
```
Closure/pole + W1 baseline (test-architect baseline was "66 passed"; suite has
grown to 80 — all green, no regressions):
```
80 passed, 1 warning in 74.44s (0:01:14)
```
Matvec bit-identity twin (P0 now uses closure-τ — the 0-drift proof). 6 passing
matvec tests stayed GREEN; the 27 xfails are PRE-EXISTING route-arounds
(cylinder matvec/sweep WDD divergence #206 + ERR-058), NOT introduced by Step A;
no xpass:
```
6 passed, 27 xfailed, 1 warning in 652.05s (0:10:52)
```
Adjoint reciprocity + curvilinear MMS structural anchors:
```
14 passed, 1 warning in 0.21s
```
Sphinx: `build succeeded.` (exit 0, only pre-existing unrelated test-file
SyntaxWarnings). Stub label/producer/test all render in HTML.

## BIT-IDENTITY VERDICT

The matvec twin stayed BIT-IDENTICAL (0 drift). Closure-τ == geometry-τ
bit-for-bit (Leg-1 `np.array_equal`, 0-ULP), so the P0 matvec contribution that
now reads closure-owned τ computes the same value as HEAD. The 6 matvec tests
that exercise the curvilinear angular branch did not move.

## LESSON (proposed for the carve crosswalk / coding-elegance Pattern 7)

A producer relocation that must stay 0-ULP is a Pattern-7 "normalise at the
definition site" move done by REPLICATION, not by sharing: the closure replicates
the factory arithmetic EXACTLY (accumulation order, the clamp, the ½ fallback)
rather than importing the factory's array, because the GOAL is to make the
closure the eventual sole owner (Step B retires the factory). The teeth of the
de-risk are the Leg-1 0-ULP `np.array_equal` against the factory AND the
structurally-independent `contamination` reference — but the cylinder needs the
NEGATIVE control (closure ≠ raw reference where the clamp bites) because the
clamp is a real producer-side transform that the unclamped reference cannot pin;
without it, a closure that silently dropped the clamp would still pass the
"equals reference" leg on the sphere and would only be caught by the
factory-equivalence leg. The mutation check proved BOTH the factory-equivalence
and the negative control catch a dropped clamp — ship them together.
