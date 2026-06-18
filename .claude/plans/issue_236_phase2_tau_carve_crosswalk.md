# Issue #236 Phase 2 — τ carve: convention crosswalk + 4-leg verification spec

**Status:** SPEC/DESIGN (test-architect, 2026-06-17). No production code, no carve.
**Branch context:** `feature/sn-spatial-angular-product` (the #236 campaign).
**Verified live @ HEAD** (this session); all baselines below are from `.venv/bin/python -O`.

This file is the L17 convention crosswalk + the 4-leg verification spec for the
relocation of the Morel–Montry angular weight τ (`tau_mm`) and its derived
coefficient contribution (`c_in`/`c_out`) from the streaming-geometry factory
onto the `MorelMontryAngularSweep` / `IdentityAngularClosure` angular closure.

---

## 0. The carve in one sentence

τ — `(μ_n − μ_{n−1/2})/(μ_{n+1/2} − μ_{n−1/2})`, BMC 2010 Eq. 43, the UNIQUE
weight exact for a flux linear in μ — is an ANGULAR-scheme property (a function
of `(μ, w)` only) but is today produced inside the streaming-GEOMETRY factory
(`reduced_operator.py`) and baked into FOUR independent rebuild sites of the
derived `c_in`/`c_out`. The carve makes the angular closure PRODUCE τ from the
`(μ, w)` grid it already binds, and consolidates the `c_in`/`c_out` algebra.

---

## 1. THE KEY FINDING — τ is consumed via a FOUR-site c_in/c_out duplication

The task brief named two `c_in`/`c_out` rebuild sites. Live grep found FOUR:

| # | Site | file:line | Path | Owner-candidate |
|---|------|-----------|------|-----------------|
| P0 | closure precompute | `pole_angular_closure.py:692-693` | matvec (`cell_contribution`) | THE intended owner |
| P1 | scalar solve | `cell_balance.py:313-314` | DD `update` (SI sweep) | inline twin |
| P2 | residual | `diamond.py:306-307` | DD `residual` (matvec @ n_mask=1) | inline twin |
| P3 | CollisionCache | `sweep_cache.py:309-310` | CumprodScan precompute | inline twin |

All four compute, byte-for-byte:
```
c_out = alpha_out / tau
c_in  = (1 - tau) / tau * alpha_out + alpha_in
```
P0 reads τ from `reduced.tau_mm` TODAY (line 660/669); P1/P2 read `st.tau_mm`
from the `StreamingTerms` packet; P3 reads `st0.tau_mm` from the first visit.

`diamond.py:295-303` ALREADY documents P2's inline rebuild as a "Phase 6 cleanup
target" — the plan of record is for the sweep route to call
`closure.cell_contribution(...)` like the matvec does (Pattern 2). The τ carve
is the natural moment to discharge that debt.

---

## 2. The τ producer / consumer / reference map (verified line numbers)

PRODUCERS (`orpheus/geometry/reduced_operator.py`):
- `:681-688` SPHERE — weight-sum edges from −1.0, `tau_raw`, **NO clamp** (W1
  declamp; `:669-680` explains). `dmu>1e-15 else 0.5` degeneracy fallback never
  fires on GL.
- `:798-815` CYLINDER — per-level, η-midpoint edges with ±sinθ endpoints,
  `tau = max(0.5, min(1.0, tau_raw))` — **clamp RETAINED** (structural τ_raw=0
  ÷0 block at the most-inward ordinate; `:783-791` + #229).
- `:495` SLAB — synthetic neutral `tau_mm=1.0` (CARTESIAN branch of the
  `streaming_terms()` accessor; no stored array).

PACKING accessor: `reduced_operator.py:486-558` packs per-cell
`StreamingTerms.tau_mm` (sphere `:519`, cyl `:554`, slab `:495`).

CONSUMERS of per-cell `StreamingTerms.tau_mm`:
- `cell_balance.py:306` → c_out/c_in rebuilt `:313-314`, denom `:319` (scalar solve)
- `diamond.py:230` (residual path), `diamond.py:305` (update path)
- `sweep_cache.py:298` → c_out/c_in `:309-310` (CumprodScan precompute)

CONSUMER that today RECEIVES the τ ARRAY (the closure that SHOULD produce it):
- `pole_angular_closure.py:660` `self._tau_per_level = (reduced.tau_mm,)` (sphere)
- `pole_angular_closure.py:669` `= tuple(reduced.tau_mm_per_level)` (cylinder)
- The closure binds `quad = sn_mesh.quad` (`:648`) → it ALREADY HAS `(mu_x,
  weights, level_indices)`, i.e. everything to PRODUCE τ.

STRUCTURALLY-INDEPENDENT τ REFERENCE (independent code path, same BMC weight):
- `contamination.py:156-201` `morel_montry_weights(quad, geometry)` — uses
  `_cell_edge_cosines` (`:43-92`). **UNCLAMPED** for both geometries.
  - SPHERE: bit-identical to producer (both unclamped, both weight-sum edges).
  - CYLINDER: equals producer's PRE-CLAMP `tau_raw` (both η-midpoint edges via
    the SAME `quad.eta == quad.mu_x` column-0 data) — producer then clamps, so
    the reference is valid for `tau_raw`, NOT the shipped clamped value.
  - **No test currently consumes `morel_montry_weights` as a τ oracle** — only
    `contamination_beta` is tested (`test_quadrature.py:375-388`). The
    producer-equivalence test (Leg 1) is the FIRST τ-oracle consumer.

---

## 3. THE CONVENTION CROSSWALK (deliverable A) — see the closeout message.

## 4. THE 4-LEG VERIFICATION SPEC (deliverable B) — see the closeout message.

---

## 5. NEW test skeletons (designed; NOT production)

### 5.1 Producer-equivalence (Leg 1) — NEW file
`tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`

```python
r"""Issue #236 Phase 2 — τ producer-equivalence gate.

The MorelMontryAngularSweep closure must PRODUCE τ that equals the value
the geometry factory bakes today, EXACTLY (0-ULP). Structurally-independent
reference leg: contamination.morel_montry_weights (a different code path to
the SAME BMC-2010-Eq.43 weight).

Mode-8: every assertion is np.testing / pytest.fail (fires under -O).
"""
from __future__ import annotations
import numpy as np
import pytest
from orpheus.numerics.quadrature import Quadrature
from orpheus.derivations.discrete.sn.contamination import morel_montry_weights
# from orpheus.sn.spatial.pole_angular_closure import MorelMontryAngularSweep
# from orpheus.geometry.reduced_operator import spherical_streaming, cylindrical_streaming
# (import the post-carve closure-τ accessor when it exists; see NOTE below.)


@pytest.mark.foundation
@pytest.mark.parametrize("N", [8, 16])  # >=2 quadrature orders
def test_sphere_closure_tau_equals_geometry_factory_0ulp(N):
    """Closure-produced τ == geometry-factory τ, bit-for-bit, sphere."""
    quad = Quadrature.gauss_legendre(N)
    # tau_geom  = <spherical_streaming(...).tau_mm>           # producer today
    # tau_close = <closure.tau_for_level(0)>                  # post-carve owner
    # assert np.array_equal(tau_close, tau_geom)              # 0-ULP, same math
    ...


@pytest.mark.foundation
@pytest.mark.parametrize("N", [8, 16])
def test_sphere_tau_matches_independent_reference(N):
    """Closure-produced τ == contamination.morel_montry_weights (sphere).

    Structurally-independent: contamination.py is a different code path to the
    SAME unclamped BMC weight. Sphere is UNCLAMPED on both sides → 0-ULP.
    """
    quad = Quadrature.gauss_legendre(N)
    tau_ref = morel_montry_weights(quad, "spherical")
    # tau_close = <closure.tau_for_level(0)>
    # np.testing.assert_array_equal(tau_close, tau_ref)   # sphere: exact
    ...


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16])
def test_cyl_closure_tau_equals_geometry_factory_0ulp(n_phi):
    """Closure-produced τ == geometry-factory τ, bit-for-bit, cylinder.

    Cylinder is CLAMPED [½,1] on both producer and (post-carve) closure →
    0-ULP equality is still required (same math, new owner). The clamp is a
    producer-side transform the closure MUST replicate.
    """
    quad = Quadrature.product(n_mu=4, n_phi=n_phi)
    # for p: assert np.array_equal(tau_close[p], tau_geom[p])
    ...


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16])
def test_cyl_tau_clamp_is_the_only_difference_from_reference(n_phi):
    """Cylinder: closure τ == clamp(contamination τ_raw).

    NEGATIVE-CONTROL companion (vv anti-pattern #11): the reference is the
    UNCLAMPED τ_raw; asserting the closure τ equals the RAW reference must
    FAIL where the clamp bites (τ_raw < 0.5 at the most-inward ordinate), and
    equals clamp(τ_raw) everywhere. This pins the clamp as a real transform,
    not an accident.
    """
    quad = Quadrature.product(n_mu=4, n_phi=n_phi)
    tau_raw_ref = morel_montry_weights(quad, "cylindrical")  # list[(M,)], unclamped
    tau_clamped_ref = [np.clip(t, 0.5, 1.0) for t in tau_raw_ref]
    # for p: assert np.array_equal(tau_close[p], tau_clamped_ref[p])
    # at least one level must have min(tau_raw) < 0.5 so the clamp is exercised:
    assert any(t.min() < 0.5 for t in tau_raw_ref), (
        "fixture too weak — clamp never bites; pick a quadrature where "
        "the most-inward ordinate has tau_raw < 0.5"
    )


@pytest.mark.foundation
def test_identity_closure_tau_is_neutral_one():
    """Cartesian: IdentityAngularClosure supplies τ = 1.0 (no redistribution).

    POSITIVE control: the neutral element is produced WITHOUT a geometry branch
    (the closure type IS the dispatch). τ=1.0 ⇒ c_out = α_out/1 = 0 (α=0 slab),
    c_in = 0 — the M-M contribution vanishes identically.
    """
    # tau = <IdentityAngularClosure(...).tau_for_level(0)>
    # np.testing.assert_array_equal(tau, np.ones_like(tau))
    ...
```

NOTE on the closure-τ accessor name: the post-carve API surface that exposes the
closure-produced τ is a DESIGN decision for the implementer (e.g. a
`tau_per_level` property, or `tau_for_level(p)`). The test asserts the VALUE
contract; the accessor name is a fill-in.

### 5.2 c_in/c_out single-source invariant (Leg 1, WIDE-scope only) — NEW
`tests/sn/sweep/core/test_mm_coefficient_single_source.py`

```python
r"""Issue #236 Phase 2 (WIDE scope) — c_in/c_out single-source invariant.

Only meaningful if the carve consolidates c_in/c_out onto ONE owner. The
invariant: every consumer's (c_in, c_out) equals the closure's, bit-for-bit.
A tripwire against the FOUR-site duplication re-opening.
"""
@pytest.mark.foundation
@pytest.mark.parametrize("coord", ["sphere", "cyl"])
def test_all_consumers_read_closure_c_coefficients(coord):
    # closure_c_out, closure_c_in = <closure._c_out_per_level, ._c_in_per_level>
    # for each former rebuild site now reading the closure:
    #   assert np.array_equal(site_c_out, closure_c_out)
    #   assert np.array_equal(site_c_in,  closure_c_in)
    ...
```

---

## 6. Live baselines captured this session (paste-back, L12)

See the closeout message §"Live baselines" for the verbatim stdout fences.

---

## 7. STEP A DONE + Step B/C dependency audit (explorer, 2026-06-17)

**Step A committed** `cdb6488` (feat) + `1a084fe` (chore): the closure PRODUCES τ
(`morel_montry_tau_per_level`, `pole_angular_closure.py:514-600`) and uses it for its own
matvec contribution (P0). Bit-identical; geometry factory still bakes an identical τ for the
sweep path (parallel run, Leg-1-gated). elegance PASS-WITH-NITS + qa SUPPORTED (mutation-proven).

### The explorer audit reshaped Step B/C — TWO findings beyond the crosswalk:

**(α) A FIFTH τ consumer the 4-c-site framing missed.** `DD.update:230`
(`tau = visit.streaming_terms.tau_mm`) reads τ for the ANGULAR WDD recurrence
`(ψ_avg − (1−τ)ψ_in)/τ` — NOT a c-rebuild. So retiring the geometry-side τ producer is
GATED on re-sourcing this consumer's τ from the closure too. The c-fold collapse alone does
NOT remove the geometry τ producer.

**(β) The closure-access seam.** P3 (`sweep_cache.GeometryCoefficients.from_mesh_and_quad`)
already holds `sn_mesh` → reaches `sn_mesh.pole_angular_closure` with ZERO plumbing (the free
seam). P1 (`cell_balance_terms`) + P2 (`DD.residual`) run through the deliberately-stateless
`DD` scheme (no closure ref reachable). The clean threading: **`dag_walk` attaches the
closure's angular contribution (+ closure-τ) onto each `CellVisit`** (the mesh-aware producer
that already carries `streaming_terms`); P1/P2 + DD.update:230 read it off the visit. Keeps
`DD` mesh-free; matches the documented `diamond.py:289-303` TODO. The angular recurrence STAYS
in DD.update — it just reads closure-sourced τ via the visit (no need to carve the recurrence
itself).

**α is OUT of scope** (geometry: the face-area dome recursion `α_{n+1/2}=α_{n-1/2}−w_n μ_n`).
The closure reads α from `reduced` (as today) and produces c from closure-τ + geometry-α —
0-ULP equal to the consumers' inline c. **L20 bonus:** the sole production readers of
`StreamingTerms.alpha_in/out` are the 3 c-sites → after the c-folds retire, those fields orphan
→ retireable in the Step-C tail (they survive `DD.update:230`, which reads `tau_mm` not α).

**Bit-identity (explorer Q2):** P1/P2 per-cell scalar; P0 closure vectorized `(M_p,)` — same
elementwise op order → per-element 0-ULP. P3's per-level→`(N,)` global gather is a pure
permutation (`level_ordinates`) → value-identical. So the c-fold collapse is ~BIT-IDENTICAL,
not principled-equivalence (the test-architect's re-association worry was the denom ASSEMBLY,
which each consumer keeps).

### Refined sub-step plan (full relocation via the CellVisit seam):

- **B1** — P3 (`sweep_cache.py:298,309-310`): source c from `sn_mesh.pole_angular_closure` via a
  NEW PUBLIC accessor (do NOT reach into `_c_in_per_level`). Gather per-level→`(N,)`. Bit-id.
  Gate: `test_cache_populator_matches_cell_balance_terms` 0-ULP + regression snapshots.
- **B2** — `dag_walk`/`CellVisit`: attach the closure's per-(cell,ordinate) angular contribution
  onto the visit; rewire P2 (`diamond.py:306-307`) to read it (delete its inline c-rebuild).
- **B3** — rewire P1 (`cell_balance.py:313-314`) to read the visit's contribution; rewire
  DD.update:230's angular recurrence to read closure-τ off the visit (the 5th consumer).
- **C** — retire the geometry-side τ producers (`reduced_operator.py:681-688` sphere,
  `:798-815` cyl, `:495` slab synthetic) + the `ReducedStreamingOperator.tau_mm`/
  `tau_mm_per_level` fields + `StreamingTerms.tau_mm` + the orphaned `StreamingTerms.alpha_in/out`.
  L20 surfaces enumerated by the explorer (production/test/orphan); migrate the τ-bit-identical
  geometry tests (`test_reduced_operator.py::test_tau_mm_*`) onto the closure producer.

Each sub-step: implement → elegance + qa → commit. Structural anchors (matvec twin, per-ordinate
MMS, adjoint reciprocity) are the correctness floor; DriftWarning-strict (from `sweep/core`+
`solve`, NOT the non-escalating regression dir) is the bit-id gate. The pre-existing
`cyl_1g_homogeneous_product_dd_n20` 3.9e-11 drift is NOT a fresh anchor.

---

## 8. STATUS — B1 + B2 DONE + committed (2026-06-18)

- **B1 DONE** (`b7fed4d` feat / `af1d074` chore): P3 c-fold via the public `c_{in,out}_per_ordinate`
  accessor. Bit-identical 0-ULP. elegance PASS-WITH-NITS + qa SUPPORTED.
- **B2 DONE** (this commit): `CellVisit.c_in/c_out` (defaulted 0.0) + the single production stamp
  `SNMesh._make_cell_visit` (all dag_walk yield sites funnel through it) + `DD.residual` reads
  `visit.c` + the **#226 binding fix** (closure consumers retyped Protocol→ABC `PoleAngularClosureBase`;
  `precompute_psi_state`/`cell_contribution`/`angular_adjoint` declared `@abstractmethod` on the ABC —
  completing the strategy contract, matching `DiscretizationSchemeBase`). Finishing fixes: the
  per-ordinate gather is now **cached** O(1) at closure construction (was O(N²·nx) per-visit re-gather —
  1.46× sweep speedup, value-identical) + a committed production-stamp catcher
  `tests/sn/sweep/core/test_cell_visit_c_stamp.py` (closes the qa vv-**Mode-11** gap — walks real
  visits, asserts `visit.c == inline rebuild` 0-ULP, mutation-verified) + the test surrogate dedup'd to
  `tests/sn/sweep/core/_c_surrogate.py`. Bit-identical 0-ULP; elegance PASS-WITH-NITS + qa SUPPORTED.

### ⚠ KEY FINDING for B3 (qa, verified): `DD.residual` has NO live production callers.
The explorer mis-classified P2 — the live matvec residual uses `cell_contribution` +
`cell_balance_for_streaming` (+ the batched `residual_kernel_batch`); the per-cell `DD.residual` (the
B2 target) is a real but **currently-unused** Krylov/GMRES capability (its only `scheme.py:557` "caller"
is a docstring). So B2 was a **zero-live-risk foundation step** (hence its trivially-green matvec twin).
**The `CellVisit.c` mechanism becomes LIVE in B3** when P1 (`cell_balance_terms` via `DD.update`, which
IS the live sweep path) reads `visit.c`. B3 = `DD.update` reads `visit.c_in/c_out` (already stamped) and
passes them to `cell_balance_terms` instead of `cell_balance_terms` rebuilding from `st.tau_mm`/`alpha`;
PLUS the 5th τ consumer (`DD.update:230` angular recurrence → closure-τ off the visit, which needs a
`CellVisit.tau` field stamped the same way, enabling Step C to drop `StreamingTerms.tau_mm`).

### Owed follow-ups (track, not blockers):
- **Orphaned `PoleAngularClosure` Protocol** (elegance NIT-3 / closes the B1 ARCH-OPP): after B2 retyped
  every consumer to the ABC, the `@runtime_checkable` Protocol has zero production type/isinstance uses
  AND has already diverged (lacks the 3 new strategy methods). Retire it or reconcile with the ABC —
  Step-C tail. FILE/track as `module:sn type:improvement`.
- **Assembly collapse** (elegance NIT-1): `DD.residual`'s scalar `(ΔA/w)·c` assembly twins
  `cell_contribution`'s vectorized one (byte-pinned by the equivalence tests). Collapse when B3/Phase-6
  routes `DD.residual` through `closure.cell_contribution` (the `diamond.py:289-303` TODO).
