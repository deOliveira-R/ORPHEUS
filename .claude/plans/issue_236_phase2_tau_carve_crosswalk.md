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
