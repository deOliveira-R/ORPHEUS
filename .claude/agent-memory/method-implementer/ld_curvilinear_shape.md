---
name: ld-curvilinear-shape
description: Branch-1 rulings for the curvilinear linear-discontinuous (LD) angular-closure member — the one-measure-down Gram, the componentwise tau elimination, and the two seed variants. Read before designing the AngularClosure base class or resuming #158's curvilinear arm.
metadata:
  type: project
---

# The curvilinear LD angular-closure member — derived SHAPE

**Fact.** A Branch-1 (closed-form SymPy) derivation of the 1-D spherical and
cylindrical LD moment balances was completed 2026-08-25 for the
`AngularClosure` base-class design round. Deliverables:

- memo: `scratch/ld_curvilinear_shape_derivation.md`
- probes: `scratch/probe_ld_sphere_moments.py` (`[M]` 9 derivations, 56 checks,
  0 fail, 12 s) · `scratch/probe_ld_cylinder_moments.py` (`[M]` 4 / 20 / 0, 4 s)

⚠ These live in `scratch/` (untracked). If they are gone, the rulings below
stand on their own and the probes are ~20 min to rebuild from this file.

**Why:** the LD member's shape was needed to triangulate the angular-closure
family contract against the shipped DD member. Explicitly NOT the #158
implementation; no production code was written.

**How to apply — the rulings, in dependence order:**

1. ⭐ **The redistribution weight is a GRAM, not a scalar.** The angular term's
   `1/r` eats one power of the volume measure, so
   `R_kj = Delta A * <b_k,b_j>_{d-2}/<b_0,b_0>_{d-2}` (`d = 3` sphere, `2`
   cylinder). `R_00 == Delta A` exactly, so `redist_dAw` is its `(0,0)` corner.
   This one spelling also absorbs the geometry-dependent factor-of-2 in the
   ORPHEUS alpha normalization (sphere needs it, cylinder does not).
2. ⭐⭐ **On the SPHERE `R` is NON-diagonal (`R_01/R_00 = h/(6 r_c)`, `= 1/3` at
   the pole cell) and the off-diagonal is LOAD-BEARING** — per-ordinate
   flat-flux consistency on the SLOPE row cancels only because of it. Lumping
   `R` breaks the canonical L0 identity. On the CYLINDER `R` is diagonal
   (one measure down is flat) and lumping is a no-op — so a per-moment-row
   flat-flux gate must run on the **sphere** or it is vacuous.
3. **`tau` / `alpha` / `c_in` / `c_out` are member-INDEPENDENT.** The WDD
   elimination applies componentwise; the Jacobian of the redistribution w.r.t.
   the moment vector is exactly `c_out * R / w_n` — ONE scalar times the
   geometry matrix. No LD arm is needed in `morel_montry_tau_per_level`,
   `alpha_dome`, or the `c` algebra.
4. **The affine face recurrence survives with a SCALAR `a`** (a 1-D face is a
   point, so the spatial inflow channel is rank 1). The triple widens
   `(a, 1/denom, w)` -> `(a, A^{-1}, w_vec)`; `d1_closed_form` does NOT carry
   over (it assumes `M` diagonal, `G_11 = 0`, `R = 0` — all three false here).
5. ⛔ **Two inequivalent seed (`mu = -1`) variants**, and TWO stacked
   blindnesses. The canonical flat-flux gate cannot see the choice at either
   moment count (DD gap `propto Sigma_t psi_in - q0`, zero at the fixed point;
   both LD variants flat-flux exact). And the LD pole-cell gap is
   `(-h qhat/30, 4 h qhat/15) + O(h^2)` — `propto` the **slope source moment**,
   which ORPHEUS zeroes for the external source (#247), so a `qhat = 0` fixture
   is blind as well. The discriminating fixture needs a non-zero seed slope
   source at `r = 0`.
   The tree ships the slab-march family for DD (`carlson_inward_sweep_from_source`
   `= (h q + 2 psi_in)/(2 + Sigma_t h)`, reproduced exactly). Recommendation is
   convention-continuity, NOT correctness.
6. ⛔ **Highest-risk convention row, with NO shipped catcher:** the half-angle
   thread and the seed must be stored in the **GLOBAL** (r-increasing) frame.
   The thread crosses `mu = 0`, where the sweep direction reverses and the slope
   moment flips sign; DD never meets this because the average moment is
   sign-invariant. This is the #240 D5b-S3 bug class at a new seam.
7. **DD is NOT LD's `p = 0` case** — LD's 1-moment truncation is STEP. A base
   class that assumes otherwise is wrong from line one.
8. **OPEN (deliberately):** whether the flux-dip-eliminating `tau` is
   spatial-scheme independent. A mechanism exists that BMC's purely-angular
   `beta` cannot see (`R_01`-mediated average<->slope coupling, `O(h/r)`,
   peaking at the pole). Settle by redoing the BMC first-order expansion with
   both moments carried; corroborate with a `tau`-sweep on a fixture OUTSIDE
   `span{1, mu}` (the closure is exact there — vv-principles #24(d)).

Related: [[lessons-L014]] (the Gram-corner heuristic), [[lessons-L015]] (derive
the gate's blindness first), [[lessons-L013]] (the discrete curvilinear operator
is not additively separable).
