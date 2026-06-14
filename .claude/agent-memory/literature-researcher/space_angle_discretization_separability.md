---
name: space-angle-discretization-separability
description: Literature verdict on whether SN spatial differencing (DD/WDD/SC/LD/DGFEM) and angular differencing (quadrature + curvilinear redistribution closure) are independently selectable + combinable. Verdict = TRUE-WITH-CAVEATS. The diffusion-limit consistency literature is explicitly split into a SPATIAL paper (Larsen-Morel-Miller 1987 JCP 69) and an ANGULAR paper (Bailey-Morel-Chang 2010 NSE 165) — proof the two axes carry SEPARATE conditions. Curvilinear couples them at the IMPLEMENTATION (Hébert reuses ONE diamond weight for both face closures) but the math is still two axes.
metadata:
  type: reference
---

# Space ⊗ angle SN discretization separability — literature verdict

**Cite when**: any question about whether the spatial scheme (diamond /
weighted-diamond / step / step-characteristic / linear-discontinuous /
DGFEM) and the angular scheme (quadrature set + the curvilinear
angular-redistribution closure) are independent axes; ORPHEUS issue
#158 (spatial cell-update DD/LD/SC/Step) vs #6 (LD *angular* finite
elements); whether a tensor-product space×angle scheme registry is the
right architecture.

## VERDICT: TRUE-WITH-CAVEATS

Space and angle ARE two mathematically distinct axes. They are
independently selectable AND combinable, with TWO caveats:

1. **Curvilinear IMPLEMENTATION coupling** — in curved geometry the
   angular-redistribution term (1−µ²)/r ∂_µψ is evaluated DURING the
   spatial sweep, and the canonical textbook (Hébert) reuses ONE diamond
   weight for BOTH the spatial face closure AND the angular edge closure.
   This is a *coding* binding, not a *mathematical* necessity — B-M-C
   prove the angular weight can be replaced independently.
2. **Diffusion-limit consistency is a JOINT requirement, but
   FACTORIZES** — the limit needs (spatial condition) AND (angular
   condition), each provable on its own axis. A bad pairing can still
   break the limit even if each scheme is individually fine, because BOTH
   conditions must hold simultaneously. Independence of *selection* ≠
   independence of *consequence*.

## Q1 — Cartesian/slab: space ⊥ angle is TEXTBOOK

Discrete ordinates = angular COLLOCATION (pick N ordinates + weights);
the spatial scheme is applied per-ordinate along the sweep. In
Cartesian/slab the (1−µ²)/r redistribution term is ABSENT, so the two
axes are fully orthogonal — you can freely run S8-DD or S16-LD. This is
the universal textbook structure (Lewis & Miller §3-4; Carlson-Lathrop
1968). NOT controversial.

## Q2 — Affinity map (scheme → axis → strength/weakness)

SPATIAL axis (acts on the streaming/face closure, per-ordinate):
- Diamond-difference (DD): O(h²), but oscillatory / can go NEGATIVE.
- Step / step-characteristic (SC): POSITIVE (Denovo chose SC for
  "linear AND maintains positivity regardless of discretization params,
  provided source positive" — Evans et al. 2010), but O(h).
- Weighted-diamond (WDD/θ-WDD): tunable between DD and step.
- Linear-discontinuous FE (LD) / trilinear-discontinuous (TLD) /
  DGFEM: O(h²) + robust; the diffusion-limit-friendly family.
- Exponential / nonlinear characteristic: positivity + accuracy.

ANGULAR axis:
- Quadrature SET: level-symmetric SN vs Gauss-product vs Galerkin vs
  Lebedev vs quadruple-range product (the latter has "fewer ray
  effects than level-symmetric" — Denovo). This is pure collocation
  point/weight selection.
- Curvilinear angular-redistribution CLOSURE (the second angular
  sub-axis, only in curved geometry): Carlson starting-direction /
  α-recursion, diamond-in-angle vs Morel-Montry weighted-diamond-in-
  angle. THIS is what Bailey-Morel-Chang 2010 and Lathrop 2000 study.
- P_N / spherical-harmonics: a different angular *representation*
  (Galerkin in angle), an alternative to collocation entirely.

⚠ LD APPEARS ON BOTH AXES — this is the ORPHEUS #158-vs-#6 ambiguity
and it is LEGITIMATE, not an error. "Linear-discontinuous" names a
finite-element trial space; you can apply it in SPACE (LD spatial cell
update, #158) OR in ANGLE (LD angular finite elements, #6). They are
DIFFERENT schemes that share a name. Lathrop 2000 explicitly lists
"linear-discontinuous" as one of his five *angular* differencing
schemes; LMM 1987 and every spatial-scheme paper list LD as a *spatial*
scheme. The literature uses the same word for both axes. ORPHEUS #158
and #6 are CONSISTENT with the literature precisely because the
literature also double-uses the term. Flag to the user: the registry
must DISAMBIGUATE by axis (e.g. `SpatialScheme.LINEAR_DISCONTINUOUS`
vs `AngularScheme.LINEAR_DISCONTINUOUS`), never a single `LD` enum.

## Q3 — Curvilinear coupling: the crux. CONFIRMED by PDF read.

`scratch/literature/Hebert(2009)Chapter3.pdf` §3.9.4 (sphere, pp.141-143)
+ §3.9.3 (cylinder, pp.139-140). VERIFIED:

- **Eq. (3.431)** [sphere] diamond closure writes ONE relation that is
  SIMULTANEOUSLY the spatial AND angular diamond:
  `φ_{n,i} = ½(φ_{n,i-1/2}+φ_{n,i+1/2}) = ½(φ_{n-1/2,i}+φ_{n+1/2,i})`.
  The SAME weight ½ closes the radial face pair AND the angular edge
  pair. Cylinder analog = Eq. (3.406).
- **Eq. (3.437)** [sphere] the auxiliary DD relations:
  spatial `φ_{n,i-1/2}=2φ_{n,i}−φ_{n,i+1/2}` AND
  angular `φ_{n+1/2,i}=2φ_{n,i}−φ_{n-1/2,i}` — same factor-2 (=1/τ, τ=½).
- Eq. (3.428) cell-balance carries the redistribution divisor
  `ΔS_i/(2𝒲_n)` — the α-cascade (Eqs. 3.423-3.424) lives INSIDE the
  spatial sweep; the angular edge flux φ_{n±1/2,i} and the spatial face
  flux φ_{n,i±1/2} are solved TOGETHER per cell.

So in the canonical textbook scheme the angular differencing IS bound
to the spatial differencing — by sharing the single diamond weight τ=½.
This is the "single weighted-diamond weight reused for BOTH closures"
the brief describes. It is real.

BUT the binding is UNDONE in the literature:
- **Bailey-Morel-Chang 2010** generalize Eq. (3.431)'s angular ½ to a
  general angular weight τ_m (their Eq. 15: ψ_m = τ_m ψ_{m+1/2} +
  (1−τ_m)ψ_{m-1/2}; τ_m=½ → diamond, τ_m=1 → step) while leaving the
  spatial closure at diamond. THIS IS THE PROOF the angular weight is
  an independent knob. M-M choose τ_m (their Eq. 42) so the FIRST-order
  diffusion limit is preserved — a purely ANGULAR optimization.
- **Lathrop 2000** NSE 134(3):239-264 (DOI 10.13182/NSE00-A2114,
  PAYWALLED, not local) is THE paper that separates the axes for
  verification: he derives 5 ANGULAR schemes "WITHOUT spatial
  differencing" (solves the resulting ODEs with an ODE solver) to
  ISOLATE angular error. Direct proof the angular axis can be studied
  with the spatial axis removed entirely. See
  [[sphere_sn_spatial_order_at_origin]] for the OSTI abstract.

## Q4 — Interference / diffusion limit: the decisive pair. SMOKING GUN.

The diffusion-limit consistency literature is LITERALLY split into a
spatial paper and an angular paper — the strongest possible evidence
the conditions are separable per-axis:

- **SPATIAL condition**: Larsen, Morel & Miller, "Asymptotic solutions
  of numerical transport problems in optically thick, diffusive
  regimes," JCP **69(2):283-324 (1987)**, DOI
  **10.1016/0021-9991(87)90170-7** (⚠ NOT ...171-5; OSTI biblio
  6638316). Part II: Larsen & Morel, JCP 83(2):212-236 (1989), DOI
  10.1016/0021-9991(89)90229-5 (OA, 287 cites). These analyze the
  SPATIAL differencing scheme's diffusion limit (space scaled so cells
  are NOT optically thin). Result: schemes whose spatial limit IS a
  valid diffusion discretization (LD, etc.) are "substantially more
  accurate" than those without (DD without fixup). NOT local — paywalled
  ScienceDirect; OA copies on academia.edu.
- **ANGULAR condition**: Bailey, Morel & Chang, NSE **165(2):149-169
  (2010)**, DOI 10.13182/NSE08-66. LOCAL
  (`scratch/literature/Bailey-Morel-Chang(2010)...pdf`). They analyze
  the SN equations "discretized ONLY in angle" (p.150 step 1 — space
  kept CONTINUOUS), proving the angular axis carries its OWN diffusion
  condition. p.151 RIGHT COLUMN verbatim is the separability statement:
  *"Larsen, Morel, and Miller have shown that if a discretized form ...
  limits to an accurate discretized diffusion equation for the
  leading-order scalar flux ... However, to achieve completely
  consistent behavior ... an accurate diffusion discretization for the
  first-order scalar flux as well. This higher level of accuracy ... has
  generally been neglected in developing SN SPATIAL discretization
  schemes, but we show here that retaining full first-order consistency
  can be important for ANGULAR discretizations."* → spatial half
  (historically studied by LMM) + angular half (the B-M-C contribution),
  TWO SEPARATE conditions on TWO SEPARATE axes.

EXACT condition classification (from B-M-C, VERIFIED by PDF read):
- LEADING-order diffusion limit (ε⁰): preserved by ANY weighted-diamond
  angular scheme (step, diamond, M-M) — Eqs. (23)-(25), (32). "with ANY
  choice of weighting factors" (p.153). So leading-order is robust to
  the angular choice.
- FIRST-order diffusion limit (ε¹): the β-contamination term lives in
  the SECOND-order current Eq. (40); β = the angular sum
  Σ_m µ_m[α_{m+1/2}µ_{m+1/2}−α_{m-1/2}µ_{m-1/2}] (Eq. 41). β=0 ONLY for
  M-M weights (Eq. 42). β≠0 for step/diamond → first-order error → the
  flux DIP. β is a PURELY ANGULAR functional (only α's, µ's, weights).
- The flux DIP (the negative-slope-at-origin pathology) is therefore an
  ANGULAR artifact in curvilinear geometry, NOT a spatial one. B-M-C
  spatially-converged all tests to ISOLATE it (Fig.12 caption). This
  CORRECTS the naive intuition that "negative flux = spatial DD
  problem": the DD NEGATIVE-FLUX pathology in Cartesian IS spatial, but
  the curvilinear ORIGIN dip is angular.

JOINT-or-separable answer: The diffusion limit requires BOTH a valid
spatial condition (LMM) AND a valid angular condition (B-M-C). They
FACTORIZE — each provable on its own axis — but BOTH must hold. So a bad
pairing (e.g. good angular M-M + a spatial scheme that does NOT limit to
diffusion) STILL breaks the limit. Independence of selection is real;
"harvest the benefits of both" requires choosing a pair where BOTH
conditions hold. This is the precise sense in which they can interfere.

## Q5 — Documented mix-and-match codes

**Denovo** (Evans, Stafford, Slaybaugh, Clarno, NUCLEAR TECHNOLOGY
**171(2):171-200 (2010)**, DOI 10.13182/NT171-171, OSTI biblio 984374)
is the canonical selectable-matrix code. Documented options:
- SPATIAL (selectable): diamond-difference, θ-weighted diamond-
  difference, linear-discontinuous FE, trilinear-discontinuous FE,
  step-characteristic.
- ANGULAR (selectable, orthogonal): level-symmetric, Galerkin,
  quadruple-range product quadratures.
- Documented caveat pairing: step-characteristic chosen for positivity;
  quadruple-range product chosen for fewer ray effects than
  level-symmetric. The two menus are independent → space×angle is a
  genuine selectable matrix in a production code.
PARTISN/DANTSYS lineage (Lathrop-Carlson DTF/TWOTRAN) similarly carries
diamond/weighted-diamond spatial × multiple quadrature angular — EXPORT-
CONTROL: do NOT name equations/manual; cite only the open Denovo paper.

## Structural recommendation for ORPHEUS

The literature picture IS a tensor product: **space ⊗ angle**. A
`SpatialScheme` registry (DD/WDD/SC/LD/TLD/DGFEM) and an `AngularScheme`
registry (quadrature × curvilinear-redistribution-closure) that combine
into a `(SpatialScheme, AngularScheme)` pair is the architecture the
literature supports. Caveats to bake into the type system:
1. Disambiguate LD-spatial from LD-angular (never one `LD` enum).
2. The curvilinear redistribution closure is a SECOND angular sub-axis
   (Carlson-DD vs Morel-Montry-WDD), distinct from the quadrature set.
3. A `diffusion_limit_consistent: bool` capability flag should be a
   property of the PAIR, not either scheme alone (factorizes but joint).
4. Cartesian: the redistribution sub-axis collapses (no (1−µ²)/r term),
   so the pair is fully free.

## Cross-refs

- [[sphere_sn_spatial_order_at_origin]] — the three-orders trap (spatial
  h / angular N / diffusion ε); Lathrop 2000 + Wu 1999 details.
- [[spherical_sn_central_cell_spatial_order]] — B-M-C is ANGULAR not
  spatial; Hébert 3.9.4 + Stacey 9.9 spatial closures.
- [[phase_d_carlson_coupled_pole]] — the µ=−1 starting-direction sweep
  (the angular-closure seed), Hébert Eqs. (3.432)-(3.435).
- [[sphere_sn_pole_closure_canonical]] — Hébert §3.9.4 is the canonical
  curvilinear angular pole stencil.
