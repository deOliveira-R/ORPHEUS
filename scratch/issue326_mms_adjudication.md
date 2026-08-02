# Issue #326 — can the curvilinear MMS / L1 suite adjudicate the per-level ordinate ordering?

Status: IN PROGRESS (written incrementally; a mid-run kill preserves what is here).
Branch `refactor/operator-strategy-layers`. Host `.venv`, `python -O -m pytest`, SERIAL.

## The four candidate orderings

| tag | ordering |
|---|---|
| L | legacy — `np.argsort(mu_x)` on trig-evaluated nodes (noise-broken ties) |
| X | `np.lexsort((mu_y, mu_x))` — key made injective by ξ |
| W | azimuthal angle ω increasing — construction order, NO sort |
| S | `np.argsort(mu_x, kind="stable")` |

## Early structural reading (BEFORE any run) — the blindness hypothesis

`build_cylindrical_mms_case` (isotropic) and `build_cylindrical_anisotropic_mms_case`
both use `Quadrature.product(n_mu=4, n_phi=8)`.

The anisotropic ansatz is (`orpheus/derivations/continuous/mms/sn.py:3778`)

    psi_n(r) = (A(r) + B(r)*eta_n) / W,     eta_n = mu_x

and the manufactured source (`:3795`) is

    Q_n(r) = eta_n A' + eta_n^2 B' + xi_n^2 B/r + (SigT-SigS) A + SigT eta_n B

**Every term is a function of `eta_n` and `xi_n**2` only.** The azimuthal mirror
pair `(phi, 2*pi-phi)` shares `eta` (that IS the degeneracy of the sort key) and has
`xi -> -xi`, so it shares `xi**2` too. Therefore:

- `psi_ref` is IDENTICAL on the mirror pair,
- `Q` is IDENTICAL on the mirror pair,
- `w` is IDENTICAL on the mirror pair (product rule: `w = w_gl(mu)*2pi/n_phi`,
  constant within a level).

The mirror pair is therefore INDISTINGUISHABLE to this MMS in every input the sweep
reads. Swapping them is a pure relabeling.

Hypothesis (to be measured, not asserted): the curvilinear MMS is EXACTLY blind to the
ordering degeneracy — vv-principles Mode 7 (the ansatz nulls the very term) rather than
Mode 10 (sub-floor).

---

## [M] MEASURED — the alpha / tau structure per ordering

`orpheus/geometry/reduced_operator.py:778-786` (the cylindrical alpha recursion):

```python
for level_idx in angular_measure.level_indices:
    eta = angular_measure.mu_x[level_idx]
    w   = angular_measure.weights[level_idx]
    alpha = np.zeros(M + 1)
    for m in range(M):
        alpha[m + 1] = alpha[m] - w[m] * eta[m]
```

`orpheus/sn/sweep/pole_angular_closure.py:593-611` (the cylindrical M-M tau):

```python
eta_edge[0] = -sin_theta
for m in range(M - 1):
    eta_edge[m + 1] = 0.5 * (eta[m] + eta[m + 1])
eta_edge[M] = sin_theta
tau[m] = (eta[m] - eta_edge[m]) / deta   if abs(deta) > 1e-15 else 0.5
```

then clamped to `[1/2, 1]` (`:669-675`).

Measured for `product(n_mu=2, n_phi=8)`, level 0, `sin_theta = 0.816497`:

```
--- tag=L  n_phi=8 level 0 (sin_theta=0.816497) ---
   eta      = [-8.164966e-01 -5.773503e-01 -5.773503e-01 -1.499880e-16  4.999600e-17  5.773503e-01  5.773503e-01  8.164966e-01]
   xi       = [ 9.999199e-17 -5.773503e-01  5.773503e-01 -8.164966e-01  8.164966e-01 -5.773503e-01  5.773503e-01  0.000000e+00]
   alpha    = [0.000000e+00 6.412749e-01 1.094725e+00 1.548175e+00 1.548175e+00 1.548175e+00 1.094725e+00 6.412749e-01 2.220446e-16]
   tau_raw  = [0.000000e+00 1.000000e+00 0.000000e+00 1.000000e+00 3.463824e-16 1.000000e+00 9.284885e-16 1.000000e+00]
   tau_clamp= [0.5 1.  0.5 1.  0.5 1.  0.5 1. ]
   alpha min=0.000000e+00  alpha[M]=2.220e-16

--- tag=X  n_phi=8 level 0 (sin_theta=0.816497) ---
   eta      = [-0.816497 -0.57735  -0.57735   0.        0.        0.57735   0.57735   0.816497]
   xi       = [ 0.       -0.57735   0.57735  -0.816497  0.816497 -0.57735   0.57735   0.      ]
   alpha    = [0.       0.641275 1.094725 1.548175 1.548175 1.548175 1.094725 0.641275 0.      ]
   tau_raw  = [0. 1. 0. 1. 0. 1. 0. 1.]
   tau_clamp= [0.5 1.  0.5 1.  0.5 1.  0.5 1. ]
   alpha min=0.000000e+00  alpha[M]=0.000e+00

--- tag=W  n_phi=8 level 0 (sin_theta=0.816497) ---
   eta      = [ 0.816497  0.57735   0.       -0.57735  -0.816497 -0.57735   0.        0.57735 ]
   xi       = [ 0.        0.57735   0.816497  0.57735   0.       -0.57735  -0.816497 -0.57735 ]
   alpha    = [ 0.       -0.641275 -1.094725 -1.094725 -0.641275  0.        0.45345   0.45345   0.      ]
   tau_raw  = [1.079009 0.292893 0.5      0.707107 0.5      0.292893 0.5      0.546918]
   tau_clamp= [1.       0.5      0.5      0.707107 0.5      0.5      0.5      0.546918]
   alpha min=-1.094725e+00  alpha[M]=0.000e+00

--- tag=S  n_phi=8 level 0 (sin_theta=0.816497) ---
   eta      = [-0.816497 -0.57735  -0.57735   0.        0.        0.57735   0.57735   0.816497]
   xi       = [ 0.        0.57735  -0.57735   0.816497 -0.816497  0.57735  -0.57735   0.      ]
   alpha    = [0.       0.641275 1.094725 1.548175 1.548175 1.548175 1.094725 0.641275 0.      ]
   tau_raw  = [0. 1. 0. 1. 0. 1. 0. 1.]
   tau_clamp= [0.5 1.  0.5 1.  0.5 1.  0.5 1. ]
   alpha min=0.000000e+00  alpha[M]=0.000e+00
```

Two facts follow, and they are the whole story:

1. **L, X, S give BIT-IDENTICAL `alpha` and BIT-IDENTICAL `tau`.** `alpha` is a partial
   sum of `w_m * eta_m` and `w` is constant within a product-rule level, so permuting two
   ordinates of EQUAL `eta` cannot move any partial sum. `tau` likewise depends only on
   the SORTED `eta` sequence. The tie-break changes only WHICH GLOBAL ORDINATE INDEX sits
   at each position — i.e. which member of the mirror pair receives `tau = 1` and which
   receives the clamped `tau = 1/2`.
2. **W is not a tie-break variant at all.** It is a different closure: `alpha` goes
   NEGATIVE (min −1.0947, the alpha-dome positivity invariant is broken) and `tau_raw[0]
   = 1.079 > 1` (outside `[0,1]`, i.e. `eta_edge` is non-monotone and the "angular cells"
   have signed widths).

## [M] MEASURED — the MMS ladders under each ordering

`orpheus/derivations/continuous/mms/sn.py` isotropic cylindrical case, the exact ladder
of `tests/sn/verification/mms/test_mms_curvilinear.py::test_sn_cylindrical_mms_converges_second_order`:

```
==============================================================================
CYL isotropic MMS (test_sn_cylindrical_mms_converges_second_order)   n_cells=[20, 40, 80, 160]
==============================================================================
[L] legacy: argsort(mu_x) on trig-evaluated nodes (noise-broken ties)
     level0 = [4, 5, 3, 6, 2, 7, 1, 0]
     errors = [2.160863081427e-03 5.389276744519e-04 1.346506256545e-04
 3.365754646092e-05]
     orders = [2.0034 2.0009 2.0002]   (4.4s)
[X] lexsort((mu_y, mu_x)) on exact nodes — injective key
     level0 = [4, 5, 3, 6, 2, 7, 1, 0]
     errors = [2.160863081428e-03 5.389276744521e-04 1.346506256551e-04
 3.365754646309e-05]
     orders = [2.0034 2.0009 2.0002]   (2.9s)
[W] azimuthal omega increasing — construction order, no sort
     level0 = [0, 1, 2, 3, 4, 5, 6, 7]
     errors = [nan nan nan nan]
     orders = [nan nan nan]   (35.0s)
[S] argsort(mu_x, kind='stable') on exact nodes
     level0 = [4, 3, 5, 2, 6, 1, 7, 0]
     errors = [2.160863081428e-03 5.389276744522e-04 1.346506256550e-04
 3.365754646321e-05]
     orders = [2.0034 2.0009 2.0002]   (2.7s)

  pairwise max |rel diff| (finite legs only):
    L vs X: 6.447141e-11   bit-identical=False
    L vs S: 6.785334e-11   bit-identical=False
    X vs S: 3.381936e-12   bit-identical=False
```

The ANISOTROPIC cylindrical MMS — the ansatz that deliberately ACTIVATES the azimuthal
redistribution term, and whose #229 "floor" IS the M-M eta-thread interpolation floor:

```
==============================================================================
CYL anisotropic MMS (test_sn_cylindrical_aniso_mms_converges_second_order)   n_cells=[10, 20, 40, 80]
==============================================================================
[L] legacy: argsort(mu_x) on trig-evaluated nodes (noise-broken ties)
     level0 = [4, 5, 3, 6, 2, 7, 1, 0]
     errors = [0.022096038325 0.019504300884 0.019116118385 0.019037806484]
     orders = [0.18   0.029  0.0059]   (1.4s)
[X] lexsort((mu_y, mu_x)) on exact nodes — injective key
     level0 = [4, 5, 3, 6, 2, 7, 1, 0]
     errors = [0.022096038325 0.019504300884 0.019116118385 0.019037806484]
     orders = [0.18   0.029  0.0059]   (1.6s)
[W] azimuthal omega increasing — construction order, no sort
     level0 = [0, 1, 2, 3, 4, 5, 6, 7]
     errors = [nan nan nan nan]
     orders = [nan nan nan]   (23.2s)
[S] argsort(mu_x, kind='stable') on exact nodes
     level0 = [4, 3, 5, 2, 6, 1, 7, 0]
     errors = [0.022096038325 0.019504300884 0.019116118385 0.019037806484]
     orders = [0.18   0.029  0.0059]   (1.8s)

  pairwise max |rel diff| (finite legs only):
    L vs X: 2.935039e-14   bit-identical=False
    L vs S: 2.879390e-14   bit-identical=False
    X vs S: 8.635918e-15   bit-identical=False
```

## [M] (1) The existing curvilinear MMS / L1 suite, run under each ordering

`python -O -m pytest -p i326_plugin <files> -v -p no:randomly` with `product_mu_phi`
monkeypatched at both import sites (`rules_product`, `directional`) for the session. Node
VALUES held bit-identical across X / W / S (all built from `roots_of_unity`); L is the
untouched production baseline.

```
mms/test_mms_curvilinear.py::test_sn_spherical_mms_converges_second_order                                                     L=PASSED  X=PASSED  S=PASSED
mms/test_mms_curvilinear.py::test_sn_cylindrical_mms_converges_second_order                                                   L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_cyl_aniso_floor_scales_with_quadrature                                        L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_sn_spherical_angular_convergence_at_fixed_mesh                                L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_sn_spherical_aniso_mms_converges_second_order                                 L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_sn_cylindrical_aniso_mms_converges_second_order                               L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_S32_clean_o_h2_full_ladder                                    L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_S16_coarse_rate_cleaner_unclamped                             L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_floor_scales_with_quadrature                                  L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_operator_admits_mms.py::test_operator_admits_isotropic_mms_per_ordinate[sphere]                          L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_operator_admits_mms.py::test_operator_admits_isotropic_mms_per_ordinate[cylinder]                        L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_operator_admits_anisotropic_mms.py::test_sphere_lc_matvec_nonflat_per_ordinate_decays_under_refinement   L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_operator_admits_anisotropic_mms.py::test_sphere_lc_matvec_redistribution_term_is_active                  L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_scattering_p1.py::test_spherical_p1_source_matches_hand_reference                                  L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_aniso_scattering_p1.py::test_cylindrical_p1_source_matches_hand_reference                                L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_pole_cell_characterization.py::test_sphere_global_L2_second_order_dual_reference                         L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_pole_cell_characterization.py::test_cylinder_global_L2_second_order                                      L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_pole_cell_characterization.py::test_sphere_pole_cell_first_order_and_Linf_dominant                       L=PASSED  X=PASSED  S=PASSED
mms/test_curvilinear_pole_cell_characterization.py::test_cylinder_pole_first_order_vs_volume_average_masked_by_midpoint       L=PASSED  X=PASSED  S=PASSED
analytical/test_cp_standoff_curvilinear.py::test_cross_check_with_cp_1g                                                       L=PASSED  X=PASSED  S=PASSED
analytical/test_cp_standoff_curvilinear.py::test_heterogeneous_sn_vs_cp_cross_check                                           L=PASSED  X=PASSED  S=PASSED
analytical/test_cp_standoff_curvilinear.py::test_cross_check_with_cp_1g_sphere                                                L=PASSED  X=PASSED  S=PASSED

L: 22 passed in 50.55s      X: 22 passed in 44.53s      S: 22 passed in 52.80s
```

**W (azimuthal ordering)** — the run was killed at 90 % (it is slow because every
cylindrical leg diverges to NaN before failing), so there is no pytest summary line. But it
executed 20 of the 22 rows and the partition is unambiguous, matching the NaN measured
directly: **every CYLINDRICAL leg FAILED, every SPHERICAL leg PASSED** (the sphere does not
consume the product rule's `LevelStructure`):

```
mms/test_mms_curvilinear.py::test_sn_spherical_mms_converges_second_order                    PASSED
mms/test_mms_curvilinear.py::test_sn_cylindrical_mms_converges_second_order                  FAILED
mms/test_curvilinear_aniso_convergence.py::test_cyl_aniso_floor_scales_with_quadrature       FAILED
mms/test_curvilinear_aniso_convergence.py::test_sn_spherical_angular_convergence_at_fixed_mesh    PASSED
mms/test_curvilinear_aniso_convergence.py::test_sn_spherical_aniso_mms_converges_second_order     PASSED
mms/test_curvilinear_aniso_convergence.py::test_sn_cylindrical_aniso_mms_converges_second_order   FAILED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_S32_clean_o_h2_full_ladder        PASSED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_S16_coarse_rate_cleaner_unclamped PASSED
mms/test_curvilinear_aniso_convergence.py::test_w1_aniso_sphere_floor_scales_with_quadrature      PASSED
mms/test_curvilinear_operator_admits_mms.py::test_operator_admits_isotropic_mms_per_ordinate[sphere]    PASSED
mms/test_curvilinear_operator_admits_mms.py::test_operator_admits_isotropic_mms_per_ordinate[cylinder]  PASSED   <- MATVEC, not the sweep
mms/test_curvilinear_operator_admits_anisotropic_mms.py::...  PASSED (both, sphere)
mms/test_curvilinear_aniso_scattering_p1.py::...              PASSED (both — source-only hand references)
mms/test_curvilinear_pole_cell_characterization.py::test_sphere_global_L2_second_order_dual_reference   PASSED [ 72%]
mms/test_curvilinear_pole_cell_characterization.py::test_cylinder_global_L2_second_order                FAILED [ 77%]
mms/test_curvilinear_pole_cell_characterization.py::test_sphere_pole_cell_first_order_and_Linf_dominant PASSED [ 81%]
mms/test_curvilinear_pole_cell_characterization.py::test_cylinder_pole_first_order_vs_volume_average_masked_by_midpoint  FAILED [ 86%]
analytical/test_cp_standoff_curvilinear.py::test_cross_check_with_cp_1g                                 FAILED [ 90%]
(killed here — the two remaining rows are the het-cylinder and the sphere CP standoff)
```

Score over the 20 executed rows: **8 FAILED, all cylindrical; 12 PASSED, all spherical** —
with the single instructive exception noted in the table,
`test_operator_admits_isotropic_mms_per_ordinate[cylinder]`, which passes because it
exercises the **matvec**, not the sweep, and so never runs the diverging M-M recurrence.


**22/22 identical, test by test.** The convergence-order numbers behind the two cylindrical
rate gates are given verbatim in the ladder blocks above (orders `[2.0034 2.0009 2.0002]`
for all three orderings on the isotropic ladder; the anisotropic cylinder gate asserts a
FLOOR BAND, not a rate, and reads `[0.0221, 0.0195, 0.0191, 0.0190]` under all three).

`tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py` was EXCLUDED — its
`test_cylinder_l1_sweep_vs_trajectory_resolvent` leg did not finish in >40 min per ordering
(the adaptive-mpmath trajectory-resolvent reference; my lesson L9's fixture-cost class).
It is an L1 standoff against a semi-analytical reference; given the measured 3e-12 / 9e-15
invariance of the underlying solves it would not discriminate either.

### VERDICT on (1) and (2): the MMS CANNOT adjudicate the tie-break

X and S order the same level DIFFERENTLY (`[4,5,3,6,2,7,1,0]` vs `[4,3,5,2,6,1,7,0]`)
and produce the same L2 error ladder to **3.4e-12 (isotropic)** and **8.6e-15
(anisotropic)** — pure FP-reduction noise, orders of magnitude below the 1e-3..1e-2 error
being measured, and NOT converging apart under refinement. There is no ordering under
which the MMS error is strictly smaller or the observed order strictly closer to
theoretical. **The MMS is blind, and it is Mode-7 blind (exact, by ansatz design), not
Mode-10 blind (sub-floor).**

The reason is structural, and it is a two-line argument:

* Every quantity in both cylindrical MMS ansatzes is a function of `eta_n` and `xi_n**2`
  only (`sn.py:3778` `psi_n = A + B*eta_n`; `:3820-3824` `Q_n = eta A' + eta**2 B' +
  xi**2 B/r + (SigT-SigS)A + SigT eta B`). The mirror pair `(phi, 2pi-phi)` shares `eta`
  (that IS the sort-key degeneracy) and shares `xi**2`. So the pair is INDISTINGUISHABLE
  in `psi_ref`, in `Q`, and in `w`.
* `alpha` and `tau` are bit-identical across all ascending-eta orderings (measured above).

So swapping a mirror pair is a **pure relabeling of two ordinates that carry identical
data**: the SET of per-ordinate solutions is unchanged and the scalar flux
`phi = sum_n w_n psi_n` is invariant. This is vv-principles Mode 12 in its exact form —
the measured functional's invariance group (relabeling within an `(eta, Q, w)`-degenerate
class) CONTAINS the mutation class. No tolerance tightening, no mesh refinement, and no
quadrature order can ever expose it through this ansatz.

**Which terms the MMS ansatz nulls** (the Mode-7 declaration the issue asks for): it does
not null the redistribution term (that is genuinely activated by the `B(r)*eta` ansatz and
IS what produces the ~1.9e-2 #229 floor) — it nulls the **xi-ODD sector**. Both ansatzes
are even in `xi`, which is exactly the symmetry that makes the two members of a
degenerate `eta` pair interchangeable.

**W is not a tie-break candidate**: it produces NaN on every cylindrical MMS leg, because
the entire M-M closure downstream (`eta_edge[0] = -sin_theta`, `eta_edge[M] = +sin_theta`,
`mu_start_per_level = -sin_theta` at `reduced_operator.py:808`) is built on the premise
that the level is a MONOTONE ASCENDING march in `eta` spanning `[-sin_theta, +sin_theta]`.
W violates that premise, so `eta_edge` is non-monotone, angular-cell widths are signed,
`tau_raw` leaves `[0,1]` and `alpha` changes sign. W is a different closure, not a
different ordering, and it cannot be dropped in.

---

## [M] (3) The alpha closed form — the claim is CORRECT, with two qualifications

Diagnostic: `derivations/diagnostics/diag_326_alpha_closed_form.py` (35 tests green).

`w_m` is constant within a product-rule level (`w = w_gl(mu_p) * d_omega`,
`rules_product.py:134`), so

    alpha_k = -w_gl * d_omega * sin(theta) * sum_{j<k} cos(omega_j)

and the Dirichlet-kernel identity closes it EXACTLY (not asymptotically):

    alpha_k = -w_gl * kappa * [ xi(omega_{k-1/2}) - xi(omega_{-1/2}) ],
    kappa   = d_omega / (2 sin(d_omega/2)) = 1 + d_omega^2/24 + ...

with `omega_{k-1/2} = (k - 1/2) d_omega` the midpoint between nodes `k-1` and `k`.

MEASURED, `product(n_mu=4, n_phi=*)`, level 0:

```
n_phi=8  order=omega rel_err_vs_closed_form=2.841e-16  alpha_min=-0.237100
n_phi=8  order=eta   rel_err_vs_closed_form=2.414e+00  alpha_min=+0.000000
n_phi=16 order=omega rel_err_vs_closed_form=3.979e-16  alpha_min=-0.209284
n_phi=16 order=eta   rel_err_vs_closed_form=2.414e+00  alpha_min=+0.000000
kappa: {4: 1.1107207, 8: 1.0261722, 16: 1.0064545, 32: 1.0016082, 64: 1.0004017}
```

Two qualifications on "alpha should equal -xi":

* the polar weight `w_gl(mu_p)` scales it (alpha is not `-xi`, it is `-w_gl * xi`), and
* `kappa` is exactly 1 only in the limit — at the production `n_phi = 8` the correction
  is **2.6 %**, so a gate asserting a bare `alpha == -xi` would assert something false.

**But the closed form holds ONLY under the omega ordering** — under the production
ascending-eta ordering it is wrong by a factor 2.414, CONSTANT in `n_phi` (structural,
not a quadrature error). So the closed form does NOT adjudicate in favour of any
production-compatible candidate: it selects W, and W is the one candidate that cannot
run.

### The resolution — on the HALF range the two criteria COINCIDE

`psi` is even in `xi` (see below), so only `omega` in `[0, pi]` is independent. There
`eta = sin(theta) cos(omega)` is strictly MONOTONE, so ascending-eta and descending-omega
are the SAME ordering, there are no ties, and both criteria hold at once.

MEASURED (fold the production level onto `omega in [0, pi]`, weight `2w` for the folded
interior partners, `w` at the two self-paired endpoints `omega = 0, pi`; total weight
preserved exactly):

```
===== n_phi=8 (half-range: 5 ordinates/level) =====
  omega order  = [3.14159265 2.35619449 1.57079633 0.78539816 0.        ]
  eta          = [-5.08374127e-01 -3.59474792e-01  3.11289374e-17  3.59474792e-01  5.08374127e-01]
  min |gap|    = 1.488993e-01   (ties would be 0)
  alpha        = [ 0.00000000e+00  1.38890128e-01  3.35310430e-01  3.35310430e-01  1.38890128e-01 -5.55111512e-17]
  alpha min    = -5.551e-17   alpha[M] = -5.551e-17
  tau_raw      = [0.         0.29289322 0.5        0.70710678 1.        ]
  alpha / xi at interior half-angle boundaries = [0.71391791 0.71391791 0.71391791 0.71391791]
    (w_gl=0.347855;  ratio/w_gl = [2.05234431 2.05234431 2.05234431 2.05234431])

===== n_phi=16 =====
  min |gap|    = 3.869768e-02
  alpha min    = +0.000e+00   alpha[M] = +1.388e-17
  ratio/w_gl   = 2.01290909 (constant across all 8 interior boundaries)

===== n_phi=32 =====
  min |gap|    = 9.768266e-03
  alpha min    = -6.245e-17   alpha[M] = -6.245e-17
  ratio/w_gl   = 2.00321638 (constant across all 16 interior boundaries)
```

So on the half range: **no ties**, `alpha` is a non-negative dome closing to zero, AND
`alpha_{m+1/2} = 2 w_gl kappa_half * xi(omega_{m+1/2})` with **no additive constant**
(the march starts at `omega = pi`, where `xi = 0`) and the factor 2 is the folded weight,
converging as `2.0523 -> 2.0129 -> 2.0032`. The user's closed form is exactly right there.

**The whole #326 degeneracy is an artefact of carrying the redundant half of the circle.**

---

## (4) Does the level span the full circle? YES — and the two half-circles INTERLEAVE

`rules_product.py:115` samples `phi = np.linspace(0, 2*pi, n_phi, endpoint=False)` — the
FULL circle, `n_phi` ordinates per `mu`-level. Sorting that by `eta` interleaves the two
half-circles. Measured level 0 for `n_phi = 8` (production, `eta`-ascending):

```
position : 0      1       2       3       4       5       6       7
omega    : pi     5pi/4   3pi/4   3pi/2   pi/2    7pi/4   pi/4    0
eta      : -.8165 -.5774  -.5774  0       0       +.5774  +.5774  +.8165
xi       :  0     -.5774  +.5774  -.8165  +.8165  -.5774  +.5774   0
```

Consecutive positions alternate between the lower (`xi<0`) and upper (`xi>0`) half-circle.

**Is the interleaving deliberate?** No — nothing in the tree states it, and the theory
page's own framing contradicts the implementation. `docs/theory/methods/sn/curvilinear_one_group.rst`
§`sn-direct-seed-circle-vs-interval` argues at length that

> "the cylinder redistributes across the **azimuthal angle** `varphi`, which lives on a
> **circle** `[0, 2pi)` — a *periodic* domain"

but the implemented closure is an INTERVAL march in `eta`: `eta_edge[0] = -sin_theta`,
`eta_edge[M] = +sin_theta`, `alpha_{1/2} = alpha_{M+1/2} = 0`, `mu_start_per_level =
-sin_theta`. A genuinely periodic axis has no edge and no starting condition. The interval
treatment is only legitimate BECAUSE `psi` is even in `xi`, which splits the circle at the
two `xi = 0` points (`omega = 0, pi`) into two mirror-image intervals — the very symmetry
that makes the mirror pair redundant. The code carries both copies and fuses them into one
`eta`-march; the theory page never discusses the 2-to-1 degeneracy this creates.

`level_symmetric` is worse: the explorer measured its cylinder levels are **4-to-1** over
`eta` (LS4 level 0: 16 ordinates, 4 distinct `eta`), so any ordering rule must answer the
4-fold case, not just the mirror pair.

---

## [M] THE CRITERION THAT DOES BITE — the xi-mirror symmetry, and it says NO ordering is right

Diagnostic: `derivations/diagnostics/diag_326_azimuthal_mirror_symmetry.py`
(6 passed, 3 xfail-strict; every xfail verified under `--runxfail` to red for ITS
documented reason).

1-D cylindrical geometry is invariant under the mirror through the plane spanned by
`z_hat` and `r_hat`, which maps `(eta, xi, mu_z) -> (eta, -xi, mu_z)`. For a `xi`-even
source the exact flux satisfies `psi(r, eta, xi, mu_z) == psi(r, eta, -xi, mu_z)`. The
product rule is closed under that mirror with equal weights (`reflection_index("y")` is an
involution — pinned as a control leg), so the SEMI-DISCRETE problem inherits it EXACTLY.
This is a closed-form statement about the operator's symmetry group: no reference solver,
no manufactured source, no second implementation. Structurally independent by construction.

MEASURED — cylindrical fixed source, **isotropic** source, reflective/vacuum,
`max_n |psi_n - psi_{mirror_y(n)}| / max|psi|`:

```
product(4,8),  nx=20, homogeneous    ->  1.190877e-01   (30.1 % relative to the local psi)
product(4,8),  nx=20, heterogeneous  ->  5.144789e-02
level_symmetric(4), nx=20            ->  3.083566e-01
```

Worst offender located (product(4,8), nx=20, homogeneous):

```
argmax defect at (n,g,cell)= (10, 0, 19)  psi= 0.07340638582  partner psi= 0.05644379524  |d|= 0.01696259
  eta,xi,mu_z of n: 5.758e-17, 0.9404322889, -0.3399810436     cell center r = 1.95
per-cell max local-relative defect:
  cell  0 r=0.050  1.8080e-03      cell 12 r=1.250  6.6757e-02
  cell  4 r=0.450  1.5365e-02      cell 16 r=1.650  1.4455e-01
  cell  8 r=0.850  3.3468e-02      cell 19 r=1.950  3.0052e-01
```

It is an ANGULAR defect, not a spatial one:

```
mirror defect vs nx: {10: 0.11266, 20: 0.11909, 40: 0.12524, 80: 0.12868}   (it GROWS)
defect by n_phi (nx=20, n_mu=4): {4: 1.341e-1, 8: 1.191e-1, 16: 9.665e-2, 32: 6.119e-2, 64: 3.373e-2}
defect by n_mu  (nx=20, n_phi=8): {2: 1.251e-1, 4: 1.191e-1, 8: 1.205e-1, 16: 1.208e-1}   (FLAT)
```

Falls with the AZIMUTHAL order, flat in the POLAR order — the same axis
`test_cyl_aniso_floor_scales_with_quadrature` already identified for the #229 aniso floor.
**The #229 "angular floor" and this symmetry violation are the same defect.**

### The mechanism, and why no ordering fixes it

`morel_montry_tau_raw_per_level` (`pole_angular_closure.py:599-611`) builds `eta_edge` as
midpoints of CONSECUTIVE ordinates. Two mirror partners share `eta` bit-exactly, so their
shared edge collapses ONTO the node — a zero-width angular cell. Measured, level 0 of
`product(2,8)`:

```
tau_raw     = [0. 1. 0. 1. 0. 1. 0. 1.]
tau_clamped = [0.5 1.  0.5 1.  0.5 1.  0.5 1. ]
```

The first partner gets `tau = 1`, the second `tau_raw = 0` which the structural `[1/2, 1]`
clamp (`:667-675`) lifts to `tau = 1/2`. **Two ordinates the geometry says are identical
receive different angular weights, and the sort's tie-break decides which is which.**

Therefore:

* Any ordering ASCENDING in `eta` puts the partners adjacent ⟹ forces the `{1, 1/2}`
  split. L, X and S all give **bit-identical** `alpha` and `tau`, and produce the SAME
  mirror defect `1.190877e-01` — only the LABEL moves.
* Any ordering NOT ascending in `eta` breaks `eta_edge` monotonicity and the alpha dome
  (W → NaN).

⟹ **There is no correct ordering to choose. The defect is the closure.**

---

## [M] What an angularly-non-trivial companion ansatz buys — and what it does not

Diagnostic: `derivations/diagnostics/diag_326_mms_ordering_blindness.py` (10 green).

A `xi`-ODD companion `psi_n = (A(r) + B(r) xi_n)/W` (same `A`, `B` shapes; manufactured
source `Q_n = [eta A' + eta xi B' - eta xi B/r + (SigT-SigS)A + SigT xi B]/W`) leaves the
symmetric sector — its mirror-pair source difference is 3.7e-1. MEASURED, `product(4,8)`:

```
lexsort=[0.07270592 0.07903657 0.0805952  0.08098395] orders=[-0.12044749 -0.02817356 -0.0069421 ]
stable =[0.08769726 0.08273873 0.08151811 0.08121451] orders=[ 0.08396894  0.02144221  0.00538303]
coarse gap=2.0619e-01  fine gap=2.8470e-03
```

So the companion **does see** the tie-break (20.6 % apart at nx=10) — unlike the
production ansatzes, which are blind at 1e-14. But it does **not adjudicate**: both
orderings have spatial order ~0 and converge to the SAME angular floor from opposite
sides (0.28 % apart at nx=80). "Which is smaller" is a coarse-mesh transient, not a limit.

An angular refinement of the same case shows the floor is real and azimuthal:
`n_phi 8/16/32/64 at nx=80 -> 8.098e-2, 4.070e-2, 2.410e-2, 1.165e-2`.

---

## [M] #325 <-> #326 COUPLING — today the tie-break is not even reachable

A control leg in the blindness diagnostic caught a vacuity in my own first probe and turned
it into a result. With production's trig-evaluated nodes
(`np.cos(np.linspace(0, 2*pi, n_phi, endpoint=False))`) the mirror pair's `eta` differ by
~1 ULP:

```
production eta gaps that should be exact ties: [1.11022302e-16 1.24515749e-16 1.11022302e-16]
```

so the "tie" is resolved by the ROUNDING NOISE before any tie-break rule can act, and
lexsort, stable and quicksort all agree — measured identical level orderings. The
tie-break becomes a free variable ONLY once the nodes are algebraically exact
(`roots_of_unity`, i.e. #325). That is why #326 blocks #325 and not the reverse, and it
means every measurement of "ordering alone" must be made on exact nodes.

---

## [M] CORRECTION to the issue's stated mechanism — it is the POLE MAP, not heterogeneity

The issue's table is reproduced EXACTLY (same script structure, exact node values, legacy
key vs exact key, `n_phi=32`, `nx=10`, reflective/reflective, isotropic source):

```
Reproducing the issue table: legacy(B) vs default-argsort(C) ordering, EXACT node values
config                              ORDERING-alone rel (psi)        (phi)
--------------------------------------------------------------------------
1G homogeneous   src=isotropic                    3.7495e-14   8.5487e-15
1G heterogeneous src=isotropic                    2.6073e-02   6.5610e-04
2G homogeneous   src=isotropic                    1.5220e-15   5.8135e-16
2G heterogeneous src=isotropic                    6.3870e-03   1.4108e-04
1G homogeneous   src=random                       1.1359e-01   8.8313e-03
1G heterogeneous src=random                       6.3316e-02   4.1224e-03
2G homogeneous   src=random                       1.3232e-02   7.0746e-04
2G heterogeneous src=random                       1.0078e-02   6.1539e-04
```

(The `psi` column matches the issue's numbers to all printed digits — `2.607e-02`,
`6.387e-03`, `3.75e-14`, `1.522e-15`. The issue's table is the ANGULAR flux; the harness
at `m16_e2e_het.py` walks `interior -> .values`, which returns `AngularFlux.values`.)

**But heterogeneity is NOT the mechanism.** Two facts:

* the effect is present on a HOMOGENEOUS problem the moment the source is `xi`-odd
  (`1G homogeneous src=random -> 1.14e-1`);
* it is absent on a HETEROGENEOUS problem when the data is `xi`-even AND the pole
  coupling is not in play.

The homogeneous+isotropic leg reads ~1e-14 because that problem's exact solution is FLAT
(reflective/reflective, uniform source) — the standard H2 degeneracy, not evidence about
heterogeneity.

The real mechanism, measured:

```
sigma (tie-break permutation) moves 16 of 64 ordinates
does sigma COMMUTE with reflection_index("x") (the pole-continuity map)? False
  non-commuting indices: [2, 6, 7, 9, 10, 14, 18, 22, 23, 25, 26, 30, 34, 38, 39, 41, 42, 46, ...]
   n=2: rx[sigma[n]]=18  vs  sigma[rx[n]]=14   (eta,xi)_n=(+0.7543,+0.3125)
   n=6: rx[sigma[n]]=22  vs  sigma[rx[n]]=10   (eta,xi)_n=(+0.3125,+0.7543)
```

and the direct relabeling test on the 1G heterogeneous isotropic case:

```
RELABELING TEST  psi_C(n) ?= psi_B(sigma(n)):
   max|psi_C - psi_B[sigma]| / max|psi| = 7.576e-03      <-- NOT a relabeling
   max|psi_C - psi_B|        / max|psi| = 2.607e-02
SCALAR FLUX  max|phi_B-phi_C|/max|phi| = 6.561e-04
```

`loss_representation/__init__.py:4189` seeds each outward ordinate's `r=0` inflow from
`pole_outflow[reflection_index("x")[n]]`. A tie-break swap acts inside ONE `eta` class;
the pole map sends the ordinate to the `-eta` class, where nothing forces the tie-break to
have made the same choice. So the tie-break leaks through the pole coupling even when
every ordinate it permutes carries identical data. Without that coupling the swap would be
an exact relabeling and `phi` would be BIT-invariant (which is what the pure-mirror
diagnostic measures: `phi` moves 2.7e-16 on the vacuum-BC isotropic case).

### And the pole map itself carries the same unenforced symmetry

MEASURED, `product(2, 8)`:

```
reflection_index("x")   :  (eta, xi) -> (-eta, +xi)      # what the code uses at :4189
omega -> omega + pi     :  (eta, xi) -> (-eta, -xi)      # the straight-line axis crossing
rx == compose(ry, rot)?  True                            # they differ by EXACTLY the xi mirror
```

A ray crossing the axis keeps its LAB direction while the local `(r_hat, phi_hat)` frame
rotates by `pi`, so both components flip. The shipped map is correct up to the `xi`
symmetry — the same symmetry the solution measurably does not have. I am flagging this as
a structural observation, not a confirmed bug: it needs its own adjudication (it is pinned
by `test_cylinder_pole_map_and_axis_crossing_differ_by_exactly_the_xi_mirror`).

---

## VERDICT

1. **The curvilinear MMS / L1 suite cannot adjudicate the ordering.** X, S and L give the
   same error ladder to 3e-12 (isotropic) / 9e-15 (anisotropic), the same convergence
   orders, and the same magnitude bands. The blindness is EXACT (Mode 7 by ansatz design),
   not sub-floor. Do NOT trust it as the adjudicator.
2. **There is no ordering under which the MMS error is smaller or the order closer to
   theory.** Nor is there under the `xi`-odd companion, which sees the tie-break (20.6 %
   at nx=10) but converges to the SAME angular floor from both sides (0.28 % at nx=80,
   spatial order ~0 for both).
3. **The `alpha` closed form is CORRECT** — `alpha_k = -w_gl kappa [xi_{k-1/2} -
   xi_{-1/2}]`, exact to 3e-16 — but ONLY under the azimuthal ordering, which cannot run
   (alpha changes sign, `tau_raw` leaves `[0,1]`, the solve NaNs). Under the production
   ascending-`eta` ordering the closed form is off by 2.414x, flat in `n_phi`.
4. **The level spans the FULL circle and the two half-circles interleave.** Not
   deliberate: the theory page argues the axis is a periodic circle while the implemented
   closure is an interval march in `eta` bracketed by `+-sin(theta)`.
5. **The criterion that DOES bite is the `xi`-mirror symmetry**, and it says NO candidate
   ordering is correct: the defect magnitude (`1.19e-1`) is IDENTICAL across all
   ascending-`eta` orderings; only the label moves. The defect is the M-M 1-D `eta`-march
   on a level with duplicate `eta` — i.e. #229 — plus the `[1/2, 1]` clamp that turns the
   resulting zero-width angular cell into an arbitrary `{1, 1/2}` split.
6. **On the HALF range every criterion coincides and the ordering is unique.**
   `omega in [0, pi]` is the independent half; `eta` is strictly monotone there (no ties);
   `alpha` is simultaneously a non-negative dome AND exactly `2 w_gl kappa xi` at the
   half-angle boundaries; and the `xi`-mirror symmetry holds by construction because only
   one member of each pair exists. That is the constructive exit.

### Recommendation

Make the key injective (`lexsort`) as the issue's ruling of record says — it is necessary,
it removes the noise-dependence, and it is free. But do NOT record it as "adjudicated
correct": nothing adjudicates it, because there is no correct choice at the level of the
ordering. The adjudicable question is the CLOSURE, and the `xi`-mirror gate is the
instrument for it.

### Convergence with the parallel literature investigation (two independent grounds)

`scratch/issue326_alpha_theory.md` (a parallel, literature-grounded, read-only dispatch)
reached the same place from a structurally independent direction, which is what makes this
a closed case rather than two probes agreeing (`vv-principles` §1):

* its **F1** finds Hebert (2009) Eq. (3.398)/(3.399) *defines* `alpha_{p,q+-1/2} = W_p *
  eta_tangential(omega_{q+-1/2})` — the closed form is the literature's own definition, at
  REAL half-angle boundaries in `omega`. My ground is numerical and independent: the
  Dirichlet-kernel identity, verified to 2.8e-16.
* its **F2/§3** finds Hebert §3.9.3 and BMC Eq. (52) put the level on a HALF range
  `0 <= omega <= pi` with doubled weights, and tabulates three self-consistent designs of
  which ORPHEUS implements none. My ground is the numerical fold: on the half range there
  are no ties, `alpha` is simultaneously a dome and exactly `2 w_gl kappa xi`, and the
  ordering is unique.
* its **§3 design table** and my measurement agree that (B) full-circle-`phi`-monotone has
  a sign-changing `alpha` and `c_out < 0` — measured here as NaN on every cylindrical leg.
* it independently reaches "lexsort makes the order deterministic but not correct", which
  is the same verdict as §VERDICT.1 above.

Where this note adds: the MEASURED blindness of the MMS suite (its §4 predicted it, this
measures it), the xi-mirror symmetry as a reference-free adjudicator with a magnitude
(`1.19e-1`, ordering-INVARIANT), the xi-odd companion's see-but-cannot-adjudicate result,
the half-range numerical fold, the #325 reachability coupling, and the pole-map leak path
+ non-commutation.

### Diagnostics — PROMOTED to the permanent suite 2026-08-01

The `derivations/diagnostics/diag_326_*.py` originals are DELETED (retirement is part of
promotion). Their successors:

* `tests/sn/sweep/curvilinear/test_alpha_closed_form.py` (35: 20 L1 + 15 foundation) — the
  `alpha` closed form, the ordering premise, the half-range coincidence. Now drives the
  PRODUCTION `cylindrical_streaming` recursion (the diagnostic re-implemented it locally —
  a Mode-11 hole closed at promotion). `verifies("alpha-cylindrical", "alpha-recursion")`.
* `tests/sn/sweep/curvilinear/test_azimuthal_mirror_symmetry.py` (9 + 3 xfail-strict;
  5 L1 + 7 foundation) — the `xi`-mirror invariant, the pole-map identity, the
  non-commutation. NO `verifies` edge: the property is currently BROKEN.
* `tests/sn/verification/mms/test_mms_ordering_blindness.py` (10, all foundation) — the
  Mode-7 blindness declaration made executable + the `xi`-odd companion.

The ordering swap is now ONE shared helper,
`tests/sn/_test_helpers.product_level_ordering` (a twin across three modules would have
been a convention that could drift between gates that compare against each other).

