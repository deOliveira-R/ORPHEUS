# The octant face-transmission map — an algebra of record for diamond, step and linear discontinuous

Opened 2026-09-23. Status: ⏹ LANDED 2026-09-23 at `8ff942f9` (the carve); the candidates under "Round 2 and close-out" remain open and unruled. The ontology of the LD comparison is still being searched (`plan-authoring` §0), so this is a living plan, and the design discussion with the user comes before code.

## The instruction (the user, verbatim)

*"The face-transmission proposal you made is accepted, and we shall turn it into an excellent algebra of record improving and extending it (potentially showing how Step and LD do the job comparatively)."* (2026-09-22, ruling R5 of `.claude/plans/vv_suite_layout.md`.) Then, 2026-09-23: *"we will do that after you sort what you need to do with the performance test and prepare to compact context before the algebra of record development of the case."*

## What exists today `[M]` 2026-09-23

- `derivations/sn_dd_face_transmission.py`, the last file under the repository-root `derivations/`. For one Cartesian cell in d dimensions, the diamond closure `psi_out = 2 psi_cell - psi_in` per axis, with `psi_cell` eliminated through the cell balance, gives the d x d map from an octant's d inflow faces to its d outflow faces: `Sigma_DD = (2/D) 1 w^T - I`, with `w_a = 2 |mu_a| A_a` and `D = Sigma_t V + sum_b w_b`. It proves symbolically, per d:
  - Claim 1: `spec(Sigma_DD) = {1 - 2 Sigma_t V / D} u {-1}^(d-1)`;
  - Claim 2 (the control): step closure (`psi_out = psi_cell`) gives `Sigma_step = (1/D) 1 w^T`, with `spec = {(D - Sigma_t V)/D} u {0}^(d-1)`.
  
  It also carries a characteristic-polynomial cross-check (d <= 4) and a numeric spot-check.
- What rests on it:
  - the (d-1)-fold `-1` is the undamped face sawtooth (w-weighted zero-average faces, invisible to collision) that costs 1631 sweeps on the d = 3 all-reflective box, and it is why the default iteration budget targets rho = 0.986 (`orpheus/numerics/convergence.py:289`);
  - negative eigenvalues void Varga's regular-splitting hypothesis, so boundary Gauss-Seidel is not guaranteed faster than Jacobi (#341; the tree said "regular splitting" in 11 places until 2026-08-09);
  - step sends the mode to 0, so the `-1` belongs to diamond's `-psi_in` term, not to transport.
- Cited by: `docs/theory/methods/sn/cartesian_multid.rst` (the label `dd-face-transmission-spectrum` near line 4521, the section `sn-boundary-gs-not-regular` near 4443, and lines 75, 4795, 4801); `docs/theory/foundations/discretization.rst:1270`; `orpheus/numerics/convergence.py:289`; the docstring of `tests/gates/sn/solve/test_reflective_si_iteration_budget.py`.
- **Its defects:**
  - Nothing runs it: the label's V&V status is `documented` (`.. vv-status: dd-face-transmission-spectrum documented`), with no `verifies` edge.
  - Its checks are bare `assert`s, which `python -O` strips.
  - Both its proof and its numeric check build the matrix by hand. Nothing compares against the map the production sweep applies: a transcription gap its own docstring names.
- **The import gate:** `tests/gates/derivations/test_diagnostics_resolve_their_imports.py` scans `derivations/` and names this file as the known member of that root. Moving the file retires that root: drop `DERIVATIONS_ROOT` and its control member in the same commit.
- **The sibling that sets the pattern:** `orpheus/derivations/discrete/sn/sweep_acyclicity.py` ("the algebra of record for SN sweep acyclicity"), gated by `tests/gates/sn/sweep/test_sweep_acyclicity.py`.
- **Production closures:** `DiamondDifference` (`orpheus/transport/spatial/diamond.py:165`) and `LinearDiscontinuous` (`orpheus/transport/spatial/linear_discontinuous.py:245`). There is no step closure in production: step exists only as this derivation's control. LD in d >= 2 is the bilinear (UBLD) closure, with 2^d moments per cell (Branch 1: `orpheus/derivations/discrete/sn/ld_ubld.py`). It rides the DAG wavefront, not the scan-march (#38, #240 D5b).

## Rulings

- **R5** (the user, 2026-09-22): move it to `orpheus/derivations/discrete/sn/` as an algebra of record, improved and extended, possibly comparing step and LD with diamond. The proposal it accepted:
  1. the move, `assert`s becoming raises, and the module returning the matrices and spectra it proves;
  2. a test at `l0` with `verifies("dd-face-transmission-spectrum")`: the symbolic claims for d = 2 to 4, plus a leg that extracts the face-to-face map from the production sweep (unit inflow on each face of one cell, one octant) and asserts it equals `Sigma_DD`; a diamond-to-step mutation must redden it;
  3. repoint the page and `convergence.py` citations, and retire the repository-root `derivations/`.
- Workflow `[R]`: an algebra-of-record carve. Load the `algebra-of-record` skill. The main agent writes the derivation with the user steering; the test-architect designs the gates, with `rests_on` edges; qa and the elegance-enforcer review; the archivist writes the page.

## Open questions, for the discussion at resume

1. **RULED 2026-09-23 (the user chose the recommended option): the full face-moment map of each closure, compared by spectrum.** Diamond and step act on the degree-0 (face-average) part of LD's face space. The comparison reads which eigenvalues have modulus 1, and the spectral radius on the rest. The original question follows. **LD's face map is not d x d.** A bilinear cell's face carries 2^(d-1) moments (its average and its transverse slopes), so the octant map acts on `d * 2^(d-1)` face moments, not d face averages. What is the right object to compare across closures: the full face-moment map; its restriction to face averages (which needs a closure of the transverse moments); or its spectrum on the zero-average subspace that carries diamond's sawtooth? `[HYPOTHESIS]` The comparative question the page needs answered is "does closure X carry an undamped mode at zero leakage?", which is a statement about the spectral radius and the unit-modulus eigenvalues of each map.
2. **Step's algebra of record: RULED** (the user, 2026-09-23, verbatim): *"Step is not *currently* a differenciation scheme, but implementing its algebra of record is a good way to start."* So step is a first-class member of the comparison, with its own algebra of record, and not only a control. Building it as a production closure is not ruled; its algebra of record is the start.
3. **RULED 2026-09-23 (the user accepted the split: "Yes. It's an excellent suggestion! Go ahead").** The derivation carries each scheme's whole cell response: transmission (inflow face moments to outflow), escape (cell source to outflow), inflow to cell moments, and source to cell moments. The transmission T is the block the spectral comparison reads. The production gate extracts T test-side through the existing sweep. Spelling T in production, `T = gamma_out ∘ (trace of the full field) ∘ (L+C)^-1 ∘ (extension of a trace into the full field) ∘ iota_in`, is a separate, later decision, tied to its first consumer: Krylov on the boundary trace for reflective problems, `(I - B T) g = rhs`. Census `[M]` (the explorer, 2026-09-23):
   - `(L+C).inverse()` is `SweepOperator` (`orpheus/sn/operators/streaming.py:774`), and acts on `FullField` (bulk ⊕ trace);
   - the per-face `TraceRestrictionOperator` exists (`orpheus/numerics/operator.py:2999`, bound by `AngularTraceSpace.inflow_restriction` / `outflow_restriction`);
   - two glue maps are missing: `FullField` to its trace member (and its transpose), and a whole-boundary restriction on `AngularBoundaryFlux`.
   
   The original question follows. The original question follows. **The production-extraction leg:**
   - which operator surface exposes the single-cell, single-octant inflow-to-outflow map (the sweep through `StreamingCollisionOperator.solve` on a one-cell all-reflective mesh, or the boundary operator `B` restricted to one octant);
   - whether LD's map can be extracted the same way;
   - the premise to measure first: that a one-cell mesh is admissible.
4. **RULED 2026-09-23 (the user, paraphrased): attempt LD fully symbolically first, and see how far it goes before any numerical fallback; if it is heavy but feasible, cache the result.** The original question follows. **What the extension claims.** For LD: its spectrum in closed form (SymPy), or only numerically for chosen d? The UBLD cell system is `2^d x 2^d` per cell, so a closed form may be heavy for d = 3.

## Measurements at resume `[M]` 2026-09-23

Probe: `scratch`-pad script, one cell, d = 2, `hx = 0.7`, `hy = 1.3`, level-symmetric S4, the first ordinate with `mu_x > 0` and `mu_y > 0`, mixture A 1g, all-vacuum faces, zero volume source. For each inflow face moment, a unit value was set on that moment and one sweep ran (`default_for(problem).sweep`), and the outflow trace written back into the boundary flux was read.

- **A one-cell mesh is admissible** for both closures. Production picks `ScanMarch` for diamond and `MovingFrontierWindow` for LD. The LD face slot carries 2 moments per ordinate at d = 2 (`2^(d-1)`), so its octant map is 4 x 4.
- **Diamond: the extracted map equals `Sigma_DD` exactly** (max abs difference 0.0 at `Sigma_t = 0.9`). Question 3's extraction surface exists already, with no new production API.
- **Spectra of the extracted maps:**

  | `Sigma_t` | diamond | LD |
  |---|---|---|
  | 1e-10 | {1, -1} | {1, 0, -0.469 ± 0.149i} |
  | 0.9 | {0.444, -1} | {0.461, 0, -0.346 ± 0.087i} |
  | 50 | {-0.911, -1} | {0, -0.040, -0.019, -0.025} |

  At d = 2, in this one configuration, LD has no eigenvalue at -1. Its only unit-modulus eigenvalue is the conserved 1 at zero absorption, like step. It has an exact 0 in all three regimes. In the thick cell, diamond's lead eigenvalue also goes to -1, while LD's whole spectrum goes to 0. These are single-configuration readings, not claims yet.

## Symbolic results `[M]` 2026-09-23

Scratchpad probe `ld_facemap_sym2.py`, not yet in the tree. It builds the LD octant map `T = R A^-1 E` from `orpheus/derivations/discrete/sn/ld_ubld.py`'s `assemble_ubld`, with:
- `E`: the inflow lift per axis, `B(-1) = [1, -1]` on the active axis and the 1-D mass on the transverse axes, times `|mu_a|`;
- `R`: the outflow trace, `B(+1) = [1, 1]` on the active axis and the identity on the transverse moments.

The map is homogeneous of degree 0 in `(h_a, mu_a, Sigma_t)`, so it depends only on the optical thicknesses `tau_a = Sigma_t h_a / |mu_a|`. `theta = 1/3`.

- **d = 1:** `T = (6 - 2 tau) / (6 + 4 tau + tau^2)`, the [1/2] Padé approximant of `e^-tau`. For comparison, diamond's `(2 - tau)/(2 + tau)` is [1/1], and step's `1/(1 + tau)` is [0/1]: the three closures are the [0/1], [1/1] and [1/2] Padé approximants of exact transmission `[R]` (the Padé identification is by inspection; not yet asserted).
- **d = 2** (built and factored in 0.1 s; at the probe configuration it matches the production extraction to every printed digit): `charpoly = lambda * (linear factor) * (quadratic factor)`.
  - `lambda = 0`, exactly, for every `tau`.
  - The conserved mode: `lambda_c = 2 (tau1 + tau2)(3 tau1 + 3 tau2 - tau1 tau2) / P1`, with `P1 = tau1^2 tau2^2 + 4 tau1^2 tau2 + 6 tau1^2 + 4 tau1 tau2^2 + 12 tau1 tau2 + 6 tau2^2`. Then `1 - lambda_c = tau1 tau2 (tau1 tau2 + 6 tau1 + 6 tau2) / P1 > 0`, and `1 + lambda_c` has numerator `tau1^2 tau2^2 + 2 tau1^2 tau2 + 2 tau1 tau2^2 + 12 (tau1 + tau2)^2 > 0` (hand-derived; to be asserted). So `|lambda_c| < 1` for `tau > 0`, and `lambda_c -> 1` as `tau -> 0`.
  - The quadratic `a lambda^2 + b lambda + c`, with `a = P1 - 8 tau1 tau2`, `c = 4 tau1 tau2`, `Q(1) = tau1 tau2 (tau1 + 6)(tau2 + 6)`, `Q(-1) = tau1^2 tau2^2 + 2 tau1^2 tau2 + 2 tau1 tau2^2 + 12 (tau1 - tau2)^2 + 4 tau1 tau2`, and `a - c > 0`. The Jury (Schur-Cohn) conditions `|c| < a`, `Q(1) > 0` and `Q(-1) > 0` all hold for every `tau1, tau2 > 0`, so both roots lie strictly inside the unit disk. As `tau -> 0`, `c/a <= 1/4`, so the pair's modulus stays at or below 1/2 `[R]`. As `tau -> inf`, every non-conserved root goes to 0.
  - **Verdict at d = 2 `[M]` symbolic:** LD has no undamped mode. Its only unit-modulus eigenvalue is the conserved 1 in the limit `tau -> 0`, the same as step. Diamond's `-1` is absent. In the thick limit, diamond's lead eigenvalue goes to `-1` as well, while LD's spectrum goes to 0.
- **d = 3, fully symbolic `[M]`** (scratchpad `ld_sylvester.py` and `ld_jury.py`; 22 s). Computing the 12 x 12 rational characteristic polynomial directly did not finish in over 10 minutes. The route that works is Sylvester's determinant identity: `T = R A^-1 E` with `R` of size 12 x 8 and `E` of size 8 x 12, so `det(lambda I_12 - T) = lambda^(12-8) det(lambda A - E R) / det A`, which needs only an 8 x 8 polynomial determinant. The result is `charpoly = lambda^5 * (conserved linear factor) * (three quadratics, one per pair of axes)`, with total degree 12.
  - Every factor's coefficients have no negative term, except `Q(-1)`, which carries two `-20` terms. At d = 3, pair (1, 2): `Q(-1) = [terms with positive coefficients] + 12 (x - y)^2 + 4 x y`, with `x = tau1 tau2` and `y = tau3 (tau1 + tau2)`. The other two pairs are the same identity with the axes permuted. At d = 2 it is `x = tau1`, `y = tau2`. This is a sum-of-squares certificate, so the Jury conditions hold identically, and a SymPy `expand(lhs - rhs) == 0` makes it a proof.
  - **Verdict, d = 1, 2, 3, fully symbolic:** every LD root other than the conserved one lies strictly inside the unit disk for all `tau > 0`. The conserved root lies in (-1, 1) and goes to 1 as `tau -> 0`. LD has no undamped mode.
- **The structural reason `[R]`, now the organising statement of the comparison:** step and LD are upwind closures without feedthrough. Their outflow depends on the inflow only through the cell state, so `T = R A^-1 E` factors through the cell. Its rank is at most the number of cell moments (1 for step, `2^d` for LD), which forces at least `d - 1` zero eigenvalues for step and at least `2^(d-1)(d - 2)` for LD. Diamond has a direct feedthrough, `psi_out = 2 psi_c - psi_in`, so `T_DD = (rank-one map through the cell) - I`, and the `-I` is the `-1`. The page's line that "DD's second-order accuracy and its undamped face sawtooth are one property seen twice" becomes a theorem about feedthrough.
- **A defect in the existing module:** `derivations/sn_dd_face_transmission.py` builds its step control as `Sigma_step = (1/D) 1 w^T`, with diamond's `w_a = 2 |mu_a| A_a` and `D`. Step's balance, `sum_a |mu_a| A_a (psi_c - psi_in_a) + Sigma_t V psi_c = 0`, gives `w'_a = |mu_a| A_a` and `D' = Sigma_t V + sum w'` instead. At d = 1, the module's matrix is `2/(2 + tau)`, but step is `1/(1 + tau)`. So "Claim 2" proves the spectrum of a matrix that is not step, and its lead eigenvalue is wrong by the factor 2 in `w`. The qualitative claim (`{0}^(d-1)`, no `-1`) survives. `cartesian_multid.rst` near line 4574 states it correctly (`w'`, `D'`). The move fixes the module, and a test must pin step's d = 1 value `1/(1 + tau)` against an independent derivation of the balance.

## Build state — 2026-09-23, branch `feature/face-transmission-aor` (uncommitted)

- **Correction to "What exists":** production already spells the face transmission. `_face_transmission_matrix` (`orpheus/transport/spatial/scheme.py:498`) drives each scheme's `cell_kernel_batch` one unit inflow at a time, and `DiscretizationSchemeBase.face_transmission_spectrum` classifies the closure from two probe cells. Its consumer is the loss-kernel gauge (#344); its gates are `tests/gates/transport/spatial/test_face_transmission_damping.py`. The census at the plan's opening missed it, and I told the user T had no production consumer, which was wrong. Two consequences:
  - The Branch-2 counterpart of this algebra of record is `cell_kernel_batch`.
  - The theorem now proves what the two-probe survey only samples: the verdict is cell-independent. Production LD at d = 3 reads `UNDETERMINED`, because its numpy inflow handles only axis ∈ {0, d-1}; the theorem says DAMPED.
  
  Naming defect: `_face_transmission_matrix`'s parameters `w` and `sigma_t_volume` are passed as `s_axes` (`g_a = |mu_a|/Delta_a`) and `reaction_xs` (`Sigma_t`). Fix after the gates land.
- **Branch 1 written:** `orpheus/derivations/discrete/sn/face_transmission.py`. The whole module proves in 19 s `[M]` (`python -O -m orpheus.derivations.discrete.sn.face_transmission`). It holds:
  - `CellResponse` (`A c = E psi_in + S q`, `psi_out = R c + f psi_in`), with `characteristic_factors` computed by Sylvester's identity and a direct `Poly.factor_list`;
  - `transmission_spectrum`, which runs the Jury conditions with an orthant-positivity certificate;
  - `derive_conserved_mode`, `derive_ld_factorisation`, `derive_pade_ladder`, `derive_page_closed_form` and `derive_comparison`.
  
  Variables: `g_a = |mu_a|/Delta_a`, with `Sigma_t = 1`. The conserved eigenvalue tends to 1 as `G -> oo`, the optically thin cell. (An earlier note in this plan wrote `tau -> 0`, which is the same limit.)
- **`ld_ubld.py`:** gained `inflow_lift_axis` (every axis, interior ones included) and `outflow_trace_axis`, and `assemble_inflow_axis` and `downstream_face_trace` now route through them. 35 of 35 of `test_ld_ubld_symbolic.py`, `test_ld_ubld_primitive.py` and `test_face_transmission_damping.py` pass `[M]`.
- **Retired:** `derivations/sn_dd_face_transmission.py` (`git rm`), and the repository-root `derivations/`. The import gate dropped `DERIVATIONS_ROOT` and its control member (12 of 12 pass). `convergence.py:289` is repointed.
- **Pending:**
  - the test-architect's gates (dispatched);
  - qa and the elegance-enforcer;
  - the archivist: the page (`cartesian_multid.rst` near 4443–4600, and lines 75, 4795, 4801; `discretization.rst:1270`) and the vv-status line;
  - the `_face_transmission_matrix` naming fix;
  - `dead_references`, then commit.

## Review round 1 and the restructure — 2026-09-23

- **The reviews:**
  - qa found the mathematics sound: no false certificate, and 2700 random points in g from 1e-6 to 1e6 put no root on the unit circle.
  - The elegance-enforcer found one real defect: the tests re-composed `R A^-1 E + f I`, so the module's own transmission was unpinned for LD at d ≥ 2. Negating it turned 0 of 83 tests red.
  - Both found that nothing checks the certifiers can say False: replacing either with "always True" turned 0 tests red.
  - The elegance-enforcer's insight: f is the transmission in the thick limit, and the 1-D Padé approximant's value at infinity. So diamond's −1 is the trapezoidal rule's undamped stiff mode.
  - Reports: scratchpad `qa_ft_review.md`, `elegance_ft_review.md`.
- **Ruling R6 (the user, 2026-09-23):** the thick-limit and stability claim is added: "Yes, add it".
- **Ruling R7 (the user, 2026-09-23, on naming the inflow term):** "no particular preference between making the derivation match production or production match the derivation. The important thing is that they match, but we should choose the best formulation." I compared three formulations:
  1. the current Kronecker factors;
  2. a trace form, measured at d = 1, 2, 3: `F_out = Σ|μ|γ₊ᵀM⊥γ₊`, `E = |μ|γ₋ᵀM⊥`, `R = γ₊`;
  3. a Petrov–Galerkin family: (trial, test, face space, inflow imposition), with diamond = trial {1, ξ_a}, test {1}, strong inflow; measured to reproduce `Σ_DD` at d = 1, 2, 3.
  
  The user chose 3: "Yes. Let's try that and see how it does in practice." And also: "don't eliminate the unused escape eagerly. Let's see if we could leverage that, in this or another method (and also the other blocks, inflow to cell and source to cell)."
- **The restructure `[M]`:** `face_transmission.py` now builds every closure by integrating its weak form on the reference cell. Nothing is typed, and nothing is imported from `ld_ubld` for the construction.
  - `derive_realization` proves LD's integrated weak form equals `ld_ubld`'s Kronecker UBLD exactly (A, E, S, R, and D = 0) at d = 1, 2, 3. That is two independent Branch-1 routes to one object.
  - Diamond's minimal realization carries the feedthrough D = −I, and step's and LD's carry D = 0.
  - New claims: `derive_stability_function` (Padé type, A-stability, R(∞) = f, T → fI as g → 0), `derive_particle_balance` and `derive_flat_flux_preservation`. The last two read all four blocks: transmission, escape, inflow to cell and source to cell.
  - The whole module proves in 22 s.
  - Other leverage for the blocks `[R]`: they are the response matrix of the interface-current and response-matrix methods (transmission and escape probabilities), for one cell and one ordinate. This is a cross-method link to CP and MoC.
- **Production direction, not in this commit:**
  - #503 (filed 2026-09-23): production's numpy LD inflow cannot drive the interior axis at d = 3. The trace form fixes it by construction.
  - Whether production's closures become one Petrov–Galerkin scheme object is a separate ruling, related to #158.
- **Pending:**
  - the test-architect is adapting the gates (resumed);
  - then qa and the elegance-enforcer re-review;
  - the archivist writes the page, the ERR entry for the predecessor's step defect, and the vv-status change;
  - rebuild the docs, then `dead_references` and `staleness`, then commit.

## Round 2 and close-out — 2026-09-23

- **Round-2 reviews, all fixed inline:**
  - `sp.solve` could silently drop a constraint. It is replaced by the named Schur complement `_static_condensation`, which refuses a singular `C_c`, so the feedthrough is now the identity `D = R_c C_c^-1`.
  - The closures are the data table `SCHEMES`, and `PetrovGalerkinScheme` validates its multi-indices.
  - A non-scalar `D` is refused as a scope edge (`NotImplementedError`, naming the general determinant form that is not built).
  - A-stability is proved by `_is_a_stable`: a Hurwitz denominator and the E-polynomial `E(y) >= 0`. It is its own predicate because the test-architect measured that, inside `derive_stability_function`, the Hurwitz leg could not be reached: every Padé denominator of `e^-tau` has positive coefficients.
  - `inflow_lift_axis` is renamed `inflow_scatter_axis`.
- **Gates `[M]`:** `test_face_transmission_symbolic.py` (107 tests) and `test_face_transmission_xverif.py` (73) pass 187 with 2 strict xfails (#503), including the damping file. The test-architect's battery ran 25 arms; the spec is scratchpad `ft_gate_spec.md`. All 18 `rests_on` ids resolve.
- **The archivist:**
  - the page section `sn-face-transmission` in `cartesian_multid.rst`;
  - ERR-088, the predecessor's step weights, with its catcher `test_step_1d_is_one_over_one_plus_tau_not_the_predecessors_two_over_two_plus_tau`;
  - the `documented` sentinel removed (the audit shows the label with 33 tests);
  - four `refs.bib` entries;
  - `dead_references` 0, with a positive control.
- **Open, for after the merge:** ERR-088's Status line and a `history.rst` row both need the merge hash. The archivist's ready-to-paste row: `2026-09-23 — face transmission becomes an algebra of record: step/diamond/LD as one Petrov–Galerkin family; diamond's −1 is the strong-imposition feedthrough; ERR-088 — <hash> — #341, #503`.
- **Candidates recorded, not ruled:**
  - **Lumped LD** (`a_LLD = 1/(1 + tau + tau^2/2)`, the [0/2] Padé; `discretization-transmission-pade` states it) as a fourth member of the family: a weak Q1 closure with a lumped (quadrature) mass. That would need a quadrature choice in `PetrovGalerkinScheme`. Until it exists, the label `discretization-transmission-pade` is only three quarters verified, so no `verifies` marker is added.
  - **The whole upwind DG(p) ladder** as one claim parametrised by p. At d = 1 it is the [p/p+1] Padé (Lesaint–Raviart Theorem 2). The elegance-enforcer's suggestion.
  - **Production direction:** a production closure as the value (trial, test, face, imposition), which would add step to production and read the spectrum as `R A^-1 E + D`, the gauge's (#344) property. Related to #158; needs the user's ruling. #503, production LD's interior axis at d = 3, is filed.
- **Found in passing:** three `:math:` roles nested inside bold in `discretization.rst` rendered raw (lines ~647, ~812, ~878), fixed. The class is tracked in #422, #424 and #379.

## ⏹ Landed — 2026-09-23

`8ff942f9`: the algebra of record, the gates, the page section `sn-face-transmission`, and ERR-088. The follow-up commit writes the hash into ERR-088's Status line and the SN history row. The open candidates (lumped LD as a fourth member, the DG(p) ladder as one claim, production closures as Petrov–Galerkin values, related to #158) are recorded above and wait for the user's ruling. #503 tracks production LD's interior axis at d = 3.
