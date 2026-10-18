# The characteristic references — re-architecting the Variant-α family before patching it

**Status:** LIVING PLAN, ontology being searched (`plan-authoring` §0). Opened 2026-10-05. Implementation waits until the user calls the plan polished.

**Parent:** #405, the slow reference files (`.claude/plans/test_runtime_405.md`, step 3, the chord-oracle hoist). The branch `refactor/chord-oracle-axial-lift` holds commit `a336bde4`: the hoist, bit-identical and reviewed. It is a stepping stone and is NOT merged.

## The ruling that opened this plan

The user ruled on 2026-10-05, answering the two elegance findings on the hoist:
- On where Σ_t lives: "this machinery is old. check if there is a better way to architect the general idea before patching old code."
- On the chord holding a closure: "again... so be tied to legacy code and patching old code. instead, check if there is a better way of architecting this and similar tests."

Read as: do not patch the chord oracles further. First find the right objects for the computation they perform, and then for the tests that verify it. Then decide what the existing code becomes.

## What the machinery computes (the reading to be confirmed or refuted)

The Variant-α references (the package `orpheus/derivations/continuous/trajectory_resolvent/`) solve the k-eigenvalue transport problem with specular or partially reflecting boundaries by the method of characteristics in continuous form:
1. **The angular flux at a point and direction is a backward line integral.** ψ(r, Ω) is the emission density q, attenuated by exp(−τ), integrated along the straight ray from r in direction −Ω to the boundary (the first leg, F), plus the boundary's incoming flux attenuated over the whole leg. This is the attenuated line integral, a Volterra operator along the ray. The project already names the MoC ray as a Volterra operator (CLAUDE.md, "The direction of development").
2. **The boundary closure is a billiard.** A specular or partially reflecting surface returns the ray along a reflected straight line. For a convex body of revolution the reflected ray keeps its impact parameter b, so every later bounce traverses the same chord (the bounce period, B). The infinite sum over bounces is a geometric series, the resolvent (1 − α e^{−τ_period})^{-1}: rank 1 for one surface, rank 2 for two (a shell, a slab).
3. **Every geometry is one line in a plane.** Along a straight ray at impact parameter b, the squared distance to the axis or centre is r² = b² + (s − s₀)². That holds for the sphere (the line in its plane through the centre), the cylinder's in-plane projection, the annulus and the hollow sphere; the slab is the degenerate case. Region boundaries are circles, crossed where b² + (s − s₀)² = R_k².
4. **The symmetry decides the angular fibre.** The sphere's ray depends on (r, μ) only. The cylinder's depends on (r, μ_axial, φ_az), and the axial cosine enters only through the lift λ = 1/√(1−μ²), which multiplies every length and optical depth (`[M]` 2026-10-05, the hoist: bit-identical). The 3D ray space fibres over the in-plane chords.

`[HYPOTHESIS]` The right objects are therefore: a straight line (impact parameter and closest-approach position) crossing a concentric region partition; a medium (Σ_t and the emission density per region) read along it; the attenuated line integral over it, with the lift as a scalar on the chord's measure; and the billiard closure, keyed by the geometry's symmetry group. Seven oracle classes and four segment helpers would be instances of one construction.

## Facts already measured (2026-10-05, the hoist's reviews)

- `CylinderChordOracle` ≡ `MultiRegionCylinderChordOracle` with one region, bit for bit (retired in `a336bde4`).
- At HEAD, `_chord_segments_oracle` (sphere) ≡ `_cylinder_chord_segments_2d` and `_region_at_radius_oracle` ≡ `_region_at_radius_cyl`, identical up to renaming. The sphere and cylinder first-leg helpers agree under μ = cos φ: 0 of 2000 samples differ. The multi-region sphere, routed through `InPlaneChord` with the lift set to 1, gives 0 of 48 entries differing (the elegance-enforcer, `scratch/reference_architecture/p3/oracle_hoist/elegance_report.md`).
- Σ_t is spelled twice on the cylinder oracle, once ignored; the annulus honours the per-call one (elegance F1).
- The `test_*_oracle_bit_equal` rows compare oracles with legacy facades that only construct the oracle: tautologies (X4). `_apply_operator_annulus` is such a facade.
- Every solver gate uses a symmetric Gauss-Legendre axial rule, which hides any mutation that reverses the lift across the cosines (test-architect, `scratch/reference_architecture/p3/oracle_hoist/gates/README.md`).

**Correction (2026-10-05, the inventory):** the sphere and cylinder first-leg helpers agree in structure, not bit for bit: 503 of 2000 samples differ bitwise, 0 differ in region sequence, the largest gap is 4.4e-15 cm. The "0 of 2000" above was a comparison at a tolerance.

## The inventory (explorer, 2026-10-05, at `a336bde4`; full text `scratch/characteristic_architecture/inventory.md`)

The Nexus graph was built at `70da0d4d`, before the hoist; every graph finding was re-checked by AST or grep at `a336bde4`.

1. **No production code uses the family.** 0 of 1156 modules in `orpheus/` outside the package import it (AST; positive control `_aba_reference.py:118`). Tests: 55 import statements in 38 files, plus 4 reached without an import statement (importlib strings in `_trajectory_resolvent_api.py`, dotted names in `_traced_memo_api.py:150-152`, a subprocess string in `test_ambient_state.py:144`, `test_content_identity.py:64` importing the reference suite's `ROSTER`). It imports nothing from `peierls_nystrom/` or `common/kernels`. It is pure reference code (Branch 1).
2. **One computation, many copies.**
   - The region-index, chord-at-impact-parameter and sphere first-leg helpers exist in 3, 3 and 2 copies, identical up to renaming. The copies in `greens_function.py` (`_region_at_radius`, `_trajectory_segments`, `_chord_segments`) are dead.
   - `SphereChordOracle` equals the one-region multi-region sphere oracle to 2.3e-16 (control: a 1e-9 change in Σ_t moves 96 of 96 entries): the sphere's version of the retired cylinder twin.
   - `compute_resolvent_T_rank2` has 0 callers; `apply_variant_alpha_closure_rank2` re-spells its determinant inline.
   - `_scalar_flux_from_psi` has 4 copies; the nested fission-rate functions form two exact pairs; the three `_apply_operator_*` wrappers are clones.
3. **`geometry_kind` is an 8-case string tag**, branched at 3 sites: `_dispatch_critical`, `_dispatch_fixed_source`, `reference._RAYS`. `closure_rank` is 2 only for `slab_asymmetric` (`billiard.py:299`), although the hollow sphere and the annulus also close at rank 2.
4. **The `ChordOracle` Protocol** has no `at=` parameter, and its only users are `isinstance` checks in tests. The multi-region oracles ignore `sigma_t`; every multi-region caller passes `0.0`.
5. **Tests: 370 functions in 41 files, 298 about the family** (claim kinds classified from names and docstrings, spot-read, not all read):
   - 87 SymPy identities on `origins/`, which never call the numeric code;
   - 17 wrapper or `Billiard` delegation tautologies;
   - 6 multi-group at G=1 against single-group (one oracle under two drivers);
   - 27 "a closed body gives k = k∞", which holds for any geometry;
   - 13 use the family as the reference for SN gates (10 cross-checks, 3 records);
   - about 11 compare against an independent value: three closed forms through production primitives, Sood, WM-72 (2), PS-1982, Garcia, the case-truth rows, and the first-leg line integral;
   - 46 are marked `slow`;
   - 3 multi-region cylinder K-uniform-vs-MG gates now compare one oracle class with itself (`retirement-audit` D.14/D.16).
6. **Docs:** `docs/theory/references/trajectory_resolvent.rst`, 6660 lines, has 67 labelled equations, 36 with a `verifies` marker. The core ray equations are among the 31 without one: `peierls-greens-trajectory-integral`, `-bounce-period-integral`, `-L0`, `-Lp`, `-mu-surf`, and the sphere multi-region segment laws. The word "Volterra" does not appear on the page. A dead path, `test_trajectory_resolvent_garcia2021.py`, is cited at `trajectory_resolvent.rst:6592` and `pn_method/README.md:27`.

## Existing machinery elsewhere (explorer, 2026-10-05; full text `scratch/characteristic_architecture/existing_machinery.md`)

1. **The optical depth of a straight line through concentric shells is written six times** and all six agree: `[M]` 400 random cases, seed 20261005, largest relative difference 3.6e-14 (the MoC tracer); control: a 10 % change in one radius moves tau by 3.6e-2. The six: `chord_half_lengths` (`derivations/common/kernels`), Peierls `_chord_tau_mu_sphere`, Peierls `CurvilinearGeometry.optical_depth_along_ray`, MoC `_trace_single_ray`, Variant-alpha `_cylinder_chord_segments_2d`, and the first-leg pair. The cylinder's axial lift is written three ways (Variant-alpha's lambda, `seg.length/sin_p` in `moc/core.py:216`, the Bickley-Naylor Ki_n kernels). The billiard closure is written five ways.
2. **No shared home exists** in L1, L2 or the input layer. The orbit-space machinery covers directions only (all 6 keys of `_ORBIT_CATALOGUE` are `(Sphere, .)`); `conceptual_view.rst` has no row for a ray, chord or line. The nearest candidate, `CurvilinearGeometry` (`peierls_nystrom/geometry.py:248`, L0), is string-tagged (about 66 branch sites, #419), observer-parametrised, has no density read, no lift and no closure, and its family is withdrawn (#506). `chord_half_lengths` is a precursor (lengths at b only); `chord_quadrature` already integrates over impact parameters with breakpoints at the shell radii.
3. **The production ray operator is already placed:** `.claude/plans/posing_sequence.md:1428-1431`, the user's direction of 2026-09-28: the free flight is MoC's Volterra operator, designed once in the MoC campaign with an apply form and a sample form; Monte Carlo inherits it (#534).
4. **The reference/production line.** Production CP (`cp/solver.py:46`) already imports `chord_half_lengths` from L0, and three references do too; its docstring cites an L0 test that does not exist (`test_kernels.py` has 0 hits for "chord"). Variant-alpha imports no other ray code. Shareable without an X4 shared upstream: only the pure geometry (where a line meets the circles, the lengths at b, the measure on lines), and only with its own verification against closed forms (sum of lengths = 2 sqrt(R^2 - b^2), the disc area and sphere volume as integrals over b, Cauchy's mean chord 4V/S). Owned per branch: the attenuated integral, the first-leg/period split, the billiard closure.
5. Outside the question: "lift" already names `Quotient.lift` and `BulkLift`, so the axial factor needs another name; #436 proposes one face-pairing datum for the specular boundary, MoC track links and MC walls.

## Candidate architecture (the main agent, 2026-10-05; for discussion, nothing ruled)

**C1. A shared geometric primitive, verified on its own.** A straight line in a plane, at impact parameter b from the centre, crossing a concentric partition (radii R_1 < ... < R_n): its ordered region sequence with the segment lengths, and the measure on lines (the quadrature over b with breakpoints at the radii). Pure geometry, method-agnostic, so its home is a layer below every method (candidate: `geometry/`, the input layer of shapes; or `numerics/`). Gates: closed forms only (C1 is the trusted line of `algebra-of-record`). Consumers: CP's `chord_half_lengths`, the Peierls walkers, Variant-alpha, later the MoC campaign.

**C2. Branch-owned, built on C1.** In the reference (Variant-alpha), as single implementations each:
- the attenuated line integral of a medium along a chord (the Volterra operator in reference form);
- the axial factor 1/sin(theta) of a cylinder-family ray (renamed off "lift");
- the billiard closure, its rank given by the number of reflecting surfaces the chord meets (rank 1: the solid sphere and cylinder; rank 2: the shells and the slab), so the rank is derived, not tagged.

**C3. The geometry's symmetry picks the angular fibre**, replacing the 8-case `geometry_kind` tag: the sphere's ray depends on (r, mu) through the plane of the ray; the cylinder's on (r, mu_axial, phi_az) through the in-plane chord and the axial factor; the slab is the degenerate line. One oracle, typed by the geometry, instead of 7 classes.

**C4. Tests re-architected with the code.** C1 gets closed-form gates; the attenuated integral gets one value gate against mpmath on a manufactured density; the closure gets its resolvent identity; the solver keeps the independent-value rows (Sood, WM-72, PS-1982, Garcia, case-truth). The delegation tautologies, the twin-class comparisons and the per-geometry copies of one property retire with the code they compared.

## The geometry census (explorer, 2026-10-05, at `a336bde4`; full text `scratch/characteristic_architecture/geometry_census.md`)

Population: the 385 files of `git ls-files 'orpheus/**/*.py'`; B1 is reference code under `derivations/`, B2 production.
1. **The transformations exist; the values they act on do not.** `RigidMotion` (`geometry/transformation.py:361`, 45 tests) with `on_points` (:616) and `on_directions` (:651), and the Householder `reflection` (:809, 5 production callers). `on_directions` has 0 production calls (grep; 20 test lines as control). 0 classes are named Point, Direction, Ray, Line, Surface or Plane.
2. **Surfaces have no form.** 1-D breakpoints (`StructuredGeometry`), tensor `Mesh2D`, circles rasterised by cell centre (`pwr_pin_2d`), MoC's box plus circles, an MC membership protocol, two string-tagged L0 geometries. No quadric or implicit form. CSG appears only as a promise in two docstrings (`structured_geometry.py:18`, `mc/solver.py:57`); the user ruled on 2026-09-24 that CSG is not meshed naively and unstructured meshes come from external meshers.
3. **Ray-surface intersection is almost all reference code.** Outside the SymPy origins: 13 `sqrt(sq - sq)` sites and 14 geometric discriminants (AST). Production has two: MoC `:77` (circle), `:101` (box). Probe F1 (seed 20261005, 1000 draws): Variant-alpha, Peierls `rho_max` and the MoC root agree to 7e-14; control R -> 1.01R moves the answer >= 5.4e-3. The three share no import: independent agreement.
4. **Normals and reflection.** In the spatial layer a normal is never a vector, only an `(axis, sign)` pair. No point-to-surface distance function exists.
5. **Mesh metrics.** One measure definition (`coord.py:191`) plus about 14 inline re-spellings; face areas spelled 4 ways; the "centroid" is the coordinate midpoint (`structured.py:208`); centroid-to-face is h/2 (`diffusion/operators.py:193`). Probe F2: Peierls `shell_volume_integral` matches `CoordSystem.measure` to 2.5e-13. The curvilinear volume centroid: 0 sites. Non-orthogonality: 0 hits (control `orthogon`: 35 files). On every tensor mesh the tree builds these metrics are trivial; they become real with unstructured meshes (#322, #335, #539).
6. **Measures on lines:** three realisations of one measure (`chord_quadrature`, the CP y-quadrature, MoC track spacing and weights). Cauchy's 4V/S is stated in the docs (`collision_probability.rst:391`) and computed nowhere.
7. **Point location has 8 spellings and two boundary conventions.** Probe F3: at r = r_k exactly, Peierls `which_annulus` returns k+1 and five others return k; `MCMesh` is outer-biased by reading (not probed).

## The candidate kernel (the main agent, 2026-10-05; for the user's ruling)

`[HYPOTHESIS]` A concentric geometry is a family of level sets of one function, its radial coordinate rho: rho(x) = |x - c| for the sphere, the distance to the axis |x_perp| for the cylinder, the signed distance x.n for the slab. The project already carries this as `CoordSystem` (CARTESIAN, CYLINDRICAL, SPHERICAL), with its measure, and `StructuredGeometry` is its breakpoints. So:
- **Values:** a point (affine), a unit direction, a line p + t*Omega; the point/direction split as types, acted on by the existing `RigidMotion`.
- **One crossing law:** rho(p + t*Omega)^2 = R_k^2 is a quadratic in t for all three coordinate systems (linear for the slab). For the cylinder its leading coefficient is 1 - mu_axial^2, so the axial factor 1/sin(theta) falls out of the law instead of being a separate lift. The impact parameter is the minimum of rho along the line.
- **A line through a StructuredGeometry** gives the ordered region crossings: the chord, with one boundary-point convention.
- **The measure on lines** (over b, with breakpoints at the radii), gated by Cauchy 4V/S and the volume as an integral over b.
Then the Variant-alpha and Peierls references pose their rays on the same `StructuredGeometry` the SN problem uses, and the 8-case `geometry_kind` tag dissolves into `CoordSystem` plus the boundary laws. Out of scope until a consumer exists: general CSG surfaces and BVH (user ruling 2026-09-24), the finite-volume metrics of unstructured meshes (#322, #335, #539).

## The boundary-point question (the main agent, 2026-10-05; arguments for the user's ruling)

"Which region owns r = r_k" is three questions, by what the algorithm is locating.

**Q-a. A segment of a ray (the open interval between two consecutive crossings).** The ray algorithms (the oracles, MoC, CP chords) need the region of a segment, never of a point: a segment's interior lies in exactly one region, so no convention is involved. `[M]` by reading: the oracles locate the segment's MIDPOINT (`_cylinder_trajectory_segments_2d`, `_chord_segments_oracle`), which is interior unless the segment has zero length. Desirable: the kernel makes the crossing sequence the primary output, with the region attached to each segment by the crossing order (entering a circle from outside moves one region in), so a ray never locates a point at all. A computed crossing radius is r_k to within rounding, on either side by noise; any ray algorithm whose result depends on locating that point is ill-conditioned whatever the convention.

**Q-b. A point with a direction (a particle on a surface after a distance-to-boundary step: MC, a tracker).** Physically the particle is in the region its direction points into: the sign of Omega . grad(rho) decides (outward: region k+1; inward: region k). That is deterministic and needs no epsilon nudge (the classic nudge is a known source of lost particles in Monte Carlo codes). The grazing direction (Omega . grad(rho) = 0) is decided at second order: on a circle or sphere, rho has its MINIMUM along the line at the tangency, so the line lies outside r_k on both sides: region k+1. So the honest locator takes (point, direction), and its answer is derived, not conventional.

**Q-c. A bare point (evaluating a field, binning a tally, the region of a quadrature node).**
- For a field continuous across the interface (the scalar flux; the angular flux along a fixed direction), every owner gives the same value: the convention is invisible.
- For a field discontinuous there (the cross sections, the emission density), the value at r_k is genuinely two-valued, and any convention picks one side arbitrarily. Undesirable: a pointwise comparison between a reference and a solver that picked different owners fails spuriously, which is what 8 spellings with two conventions make possible today (probe F3).
- The better pattern exists already: the producer that places a point records its region (`region_at_node`, returned by the composite per-region quadrature) instead of a consumer re-locating it (`coding-elegance` Pattern 7, normalise at the definition site). A bare-point locator is then a rare path.
- Where one is still needed, the ends decide between the two conventions. Inner-owns, region k = (r_{k-1}, r_k] with region 0 = [0, r_0], partitions the closed domain [0, R] exactly and puts the outer surface r = R in the last region, where the boundary condition acts; 5 of the 6 probed spellings use it. Outer-owns, [r_{k-1}, r_k), leaves the surface itself outside every region and needs a rule for it.
- The fully honest alternative is a typed answer: interior(k) or on_interface(k, k+1), so that a consumer of a discontinuous field must choose a side explicitly.

**Recommendation `[R]`:** the kernel answers Q-a and Q-b by construction (segments carry their region from the crossing order; a point on a surface is located with its direction), so neither ever meets the convention. For Q-c, the producers carry the region of the points they place; the one bare-point locator uses inner-owns, closed at [0, R], as the tree's single convention, and the 8 spellings retire onto it.

## The seed, first-pass design (the main agent, 2026-10-05; under W5 review, nothing built)

Every means below is `[HYPOTHESIS]` until the design review and the user rule on it.

**Home:** a new module in `orpheus/geometry/` beside `transformation.py` (`RigidMotion`, "the Euclidean group E(d)"), numpy-only, no transport vocabulary, no measures of neutrons. `[M]` name check 2026-10-05: 0 classes named Point, Direction, Line, Ray or Chord in `orpheus/`, `tests/`, `.claude/`, `docs/` (excluding `_build`); "Crossing" has 3 code and 6 prose hits, to be read before adopting it.

**Values.**
- `Point`: a point of the affine space R^d. `Direction`: a unit vector. Point - Point is a displacement; Point + t * Direction is a Point; Point + Point is not spelled. `RigidMotion.on_points` / `on_directions` become the action on these types.
- `Line`: a point and a direction, the set {p + t*Omega}. A ray is the line restricted to t >= 0 (or t <= 0 for a backward characteristic); whether that is a type or a parameter is open.

**The concentric geometry.** `CoordSystem` already names the three radial coordinates; it gains rho(x): |x - c| (SPHERICAL), the distance to the axis (CYLINDRICAL), the signed distance n.x (CARTESIAN). A `StructuredGeometry` is its breakpoints r_0 < ... < r_n (hollow when r_0 > 0) with the materials and boundary laws.

**One crossing law.** rho(p + t*Omega)^2 = r_k^2 is a quadratic in t for the sphere and cylinder (its leading coefficient is |Omega_perp|^2: 1 for the sphere, 1 - mu_axial^2 for the cylinder) and linear for the slab. The impact parameter is the minimum of rho along the line. A line through a `StructuredGeometry` gives its chord: the ordered crossings t_0 < t_1 < ..., each open segment carrying its region by the crossing order (ruled 2026-10-05).

**Point location (ruled 2026-10-05).** `locate(point, direction)`: the region the direction points into, by the sign of Omega . grad(rho), grazing to the outer side. `locate(point)`: inner-owns, closed on [0, R], the tree's single convention.

**The measure on lines.** The invariant measure on lines in the plane (impact parameter b times the line's angle); for a concentric geometry the quadrature over b with breakpoints at the radii.

**Gates (closed forms only; the kernel is then a trusted line under `algebra-of-record`, so references and production may share it without an X4 upstream):**
- the chord length through a solid disc or sphere: 2 sqrt(R^2 - b^2), and through shells the difference of two;
- the cylinder: the 3-D chord length is the in-plane length over |Omega_perp|;
- invariance: the chord is unchanged when the line and the geometry move together by any `RigidMotion`, and unchanged when only the line moves by an element of the geometry's symmetry group;
- the measure: the area of a disc and the volume of a sphere as integrals of chord length over the measure on lines; Cauchy's mean chord 4V/S;
- point location: the three rulings, each with the boundary, centre, surface and grazing cases.

**What it replaces, later phases (sizing only, not yet planned):** the 13 hand-written chord square roots and 14 discriminants outside the SymPy origins; the 8 point-location spellings; the 7 oracle classes and the `geometry_kind` tag, by a reference posed on `StructuredGeometry`; the reference test roster (C4).

## W5 review of the first-pass seed (2026-10-05)

Reports: `scratch/characteristic_architecture/w5_elegance.md` (elegance-enforcer), `w5_cross_domain.md` with probe `w5/probe_cross_domain.py` (cross-domain-attacker, seed 20261005), `seed_verification_spec.md` with probes `seed_spec_probes/` (test-architect; constants measured on a 20-line prototype, to be re-measured on the real kernel).

**Refuted in the first-pass design:**
- **"One crossing law rho^2 = r_k^2"** `[REFUTED 2026-10-05]`: with the slab's signed rho, squaring finds spurious roots, 3 of 6 on breakpoints (-1, 0.5, 2) (`[M]` both reviewers). "The impact parameter is the minimum of rho" has no answer on the slab.
- **"The invariant measure on lines in the plane, b times the angle"** `[REFUTED 2026-10-05]`: it drops the factor b on the sphere. Cauchy's mean chord tells them apart (`[M]` cross-domain: sphere b db gives 4/3, the planar measure pi/2; cylinder sin^2 gives 2, sin gives 2.467).
- **"Grazing goes to the outer side, derived not conventional"** `[REFUTED 2026-10-05]` except on the sphere: on a slab with Omega.n = 0, or a cylinder with Omega along the axis, the line lies IN the surface for its whole length (`[M]` cross-domain, every term of g(t) - pi_k vanishes).
- **"`locate(point, direction)`"** `[REFUTED 2026-10-05]` as a spelling: a computed point is never exactly on a surface, so the question belongs to the producer of the crossing (elegance F2).
- **Solving each 3-D line's own quadratic** `[REFUTED 2026-10-05]`: it re-solves the in-plane chord for every axial cosine, discarding the structure `a336bde4` exploits, and divides by zero when the line is parallel to the axis.
- **"Closed on [0, R]"** on a hollow geometry puts the cavity in region 0 (elegance F3): the domain is [r_0, R].

**Established (the facts the revision rests on):**
- The radial coordinate is the chart of the orbit space of the geometry's symmetry group G_c; its stabiliser reproduces `GEOMETRY_ANGULAR_SYMMETRY` on 4 of 4 rows (`[M]` 2026-10-01, `scratch/boundary_ontology/attack_math.md`).
- The polynomial invariant pi = X^T Q0 X (|x|^2, |x_perp|^2, n.x: degree 2, 2, 1) makes the concentric geometry a family of quadrics Q0 - pi_k E; placing it in space is the congruence H^{-T} Q H^{-1} (`[M]` max |dt| 4e-12 over 600 lines; the mutation H^T Q H differs on 522 of 600).
- The cylinder's leading coefficient equals 1 - Omega_z^2 to 2e-16 (`[M]`): the axial factor is the obliquity of the line to the orbit space.
- The invariant measure on lines is dA_perp dOmega; per chart: sphere b db, cylinder db sin^2(theta) d(theta), slab |mu| d(mu); equivalently the next-lower `CoordSystem`'s measure over b (`[M]` elegance: sphere volume via the cylindrical measure off 9e-8, control 0.78).
- A direction is a point of the existing `numerics.manifold.Sphere`, where the quadrature ordinates live.
- `StructuredGeometry` has 0 of 4 fields that place it in space (`[M]` by reading), so no invariance gate can be written on today's data.
- Conditioning (`[M]` test-architect, prototype): `sqrt(R*R - b*b)` errs by up to 2.8e14 ulp near tangency, `sqrt((R - b)*(R + b))` 1.16 ulp; |Omega_perp|^2 as `1 - Omega_z^2` errs 8.9e-5 at |Omega_perp| = 1e-6, as `Omega_x^2 + Omega_y^2` 1e-17; b from `|p|^2 - (p.Omega)^2` errs 7e-6 with the base point 1e6 away; a thin-shell segment length as a difference of crossings loses 2e-8 relative.
- Every concentric chord is symmetric about its closest approach, so full-line assertions cannot see a reversed orientation; the total length cannot see a misplaced interior crossing (test-architect).

## The seed, revised (the main agent, 2026-10-05; for the user's ruling)

1. **Crossings are solved in the orbit space, once.** A line is projected into the orbit space of G_c (the sphere: the line's own plane; the cylinder: the plane perpendicular to the axis; the slab: the normal axis). The chord is computed there once, in the conditioned forms above, with lengths measured from the closest approach; the 3-D lengths are the orbit-space lengths times the obliquity 1/|P Omega| (1 on the sphere, 1/sin(theta) on the cylinder, 1/|mu| on the slab). This is the hoist's structure, named.
2. **Values, batch-shaped:** points and displacements as arrays; directions as points of `numerics.manifold.Sphere`; a line as (direction, moment) Plücker coordinates, so two base points on one line give one line; the impact parameter is |moment|; a ray is a line with a start parameter, not a type.
3. **The crossing value.** The step that hits a surface produces a `Crossing` (parameter, breakpoint index, sense: inward or outward); the region entered is read from it. Segments carry their region by the crossing order. The bare locator is `region_containing(rho)` composed with the chart, inner-owns, on [r_0, R].
4. **Typed regions:** interior(k), the cavity (rho < r_0), the exterior (rho > R), and on_interface(k, k+1) for a line lying in an interface (P Omega = 0 at a breakpoint).
5. **Placement:** a concentric geometry is posed in space by a `RigidMotion` from its canonical frame (the frame `AngularChart`'s columns already declare); the identity by default. The centre and axis are not Enum fields.
6. **The measure on lines** per `CoordSystem` as above, integrated by the existing `chord_quadrature` (with the endpoint substitution: `[M]` plain Gauss-Legendre misses by 4.1e-5 at 16 nodes, 2.7e-16 with it). Gates: per-chart Cauchy mean chord, per-region Dirac-Cauchy (checks the region rule: `[M]` 3e-4 against a 15 % mutation), the volume as an integral over b.
7. **Decided by principle, not for ruling:** segment lengths in the conditioned form (rule: the formula with no cancellation); a non-unit direction is refused, not renormalised (parse at the boundary); the cross-check against today's spellings is a tracked L4 harness with an AST precondition, deleted in each migration commit.

## Implementation constraints found during design

- **The import edge** `geometry` -> `numerics.manifold` is allowed by `tests/gates/test_layer_imports.py` (geometry sits above L1), but `numerics.symmetry` imports `geometry.transformation`, a known package cycle (the file's comment at line 208). The kernel imports submodules (`from orpheus.numerics.manifold import Sphere`), never `from orpheus.numerics import ...`; check by inject-and-run (`plan-authoring` §6d).

## P0 API sketch (the main agent, 2026-10-05; checkpoint with the user before code)

**`Chart`** (`orpheus/geometry/chart.py`), the 1-D charts only (ruled). For E's pair (coordinate system, kept coordinates), the kept coordinates are implied by the coordinate system in 1-D ({x}, {r}, {r}); E adds the field with the (r, z) and 2-D charts.
- `orbit_coordinate(points)`: c(x) = x_0 (slab), sqrt(x_0^2 + x_1^2) (cylinder; the axis is e_z, the chart's column 2), |x| (sphere), in the canonical frame. Vectorised over (..., 3).
- `contains(motion: RigidMotion) -> bool`: G_c = {g in E(3) : c o g = c}, decided in closed form: slab, Q e_x = e_x and t_x = 0; cylinder, Q e_z = +-e_z and the in-plane part of t is 0; sphere, t = 0. (E's attack decided it on 400 random points; the closed form is exact.)
- `singular_strata`: the sphere's centre (isotropy O(3)), the cylinder's axis (D_inf_h); the slab has none. Each a value: the orbit-coordinate value 0.0 and its isotropy `SubgroupOfO3`.
- `measure(edges)`: delegates to `CoordSystem.measure`, the one definition.
- `obliquity(directions)`: 1/|P Omega| (1 on the sphere; 1/sqrt(Omega_x^2 + Omega_y^2) on the cylinder; 1/|Omega_x| on the slab), with |P Omega|^2 summed from components, never 1 - Omega_z^2.
- `line_measure(...)`: the invariant measure on lines pushed to the chart (sphere 2 pi b db; cylinder db per unit height with the direction measure; slab |mu| d mu), integrated by the existing `chord_quadrature`.

**`Line`** (`orpheus/geometry/line.py`): a batch of oriented lines in Plucker coordinates, `direction` (..., 3), unit, refused otherwise; `moment` (..., 3) = p x Omega. `Line.through(points, directions)`; `foot` (the point closest to the origin, Omega x m); `at(t)`; `moved_by(motion)`. Two base points on one line give one line.

**`ConcentricPartition`** (`orpheus/geometry/chord.py`): a `Chart`, the breakpoints r_0 < ... < r_n, and its `pose: RigidMotion` (the identity by default; ruled). Built from a `StructuredGeometry` by `ConcentricPartition.of(geometry, pose=...)`; it reads only the coordinate system and the breakpoints.
- `chord(line) -> Chord`: the line is moved into the canonical frame by the pose's inverse and projected into the orbit space.
- `region_containing(rho)`: inner-owns on [r_0, r_n]; the cavity and the exterior as typed values.

**`Chord`**, fixed shape, batch-friendly:
- curvilinear charts: `impact_parameter` b (...,), `half_chords` h_k (..., n+1) in the conditioned form (absent, NaN-free masked, where b >= r_k, so a tangency is not a crossing), `obliquity`, and the 3-D parameter of the closest approach along the line. Crossings are at s = +-h_k * obliquity about it; each `Crossing` (parameter, breakpoint index, sense) carries the region it enters, from the crossing order. The orbit-space length of region k on one side is (r_{k+1}^2 - r_k^2)/(h_{k+1} + h_k) (no cancellation), or h_{k+1} for the innermost region crossed; the cavity of a hollow body is the segment |s| < h_0.
- the slab: the orbit space is the x axis and the projection is monotone, so crossings are at s_k = (r_k - x_0)/Omega_x, segment lengths (r_{k+1} - r_k)/|Omega_x|; Omega_x = 0 gives the degenerate line (one segment, or on_interface when x_0 is a breakpoint).
- `half_line(start_parameter)`: the crossings beyond a start, for a backward characteristic (the oracle's first leg) or a forward one.
- `on_interface(k, k+1)`: P Omega = 0 with the line on a breakpoint (ruled).

**Gates:** the test-architect's spec (`scratch/characteristic_architecture/seed_verification_spec.md`), with its directional-locator rows (L2, L3, L5, L6) re-posed onto `Crossing` and the typed regions, and its NEEDS answered by the rulings above.

## P0 landed on its branch (2026-10-05)

**Commits on `feature/geometric-kernel-seed`:** `fe0ca696` (the kernel `orpheus/geometry/{chart,line,chord}.py` and 114 gates in `tests/gates/geometry/test_{line,chart,chord,line_measure,kernel_corroboration}.py`), `82ae7013` (the theory page `docs/theory/foundations/chart_and_chord.rst`, API section, concept-table rows, regenerated matrix). Merge status: read git, not this line.

**Evidence** (`[M]` 2026-10-05):
- five gate files 114 passed under `-O`, also with `-W error::RuntimeWarning`; `tests/gates/geometry` 1297 passed; layer, docstring-xref and geometry gates 1790 passed; V&V harness audit 17 passed;
- mutation battery 45 arms, each red on its target row (`scratch/characteristic_architecture/seed_gates/README.md`, v2 section); the subnormal gate's two rows red under the old forms (the main agent's in-process mutations);
- Sphinx `-E -W --keep-going` 0 warnings; `dead_references` 0 of 66;
- pyright 0 errors on `orpheus/geometry/`.

**Review:** qa and elegance-enforcer, two rounds each, every finding resolved (`seed_qa.md`, `seed_elegance.md`). Defects found and fixed during the carve, each now gated: a NaN direction admitted; a NaN coordinate given a region index; an inf direction warning before refusal; posed chord parameters measured on the canonical line; negative exterior codes indexing a material table (the user's ruling: out-of-range codes); 0 * inf in untraversed slots and -inf * 0 in the impact parameter at subnormal |P Omega| (the radial image now carried in orbit-space units).

**Design changes during the carve, ruled or by principle:** the directional locator retired (a `Crossings` entry carries the region entered); `Chart` derived from (kept columns, group) (ruled); `Chart.image` -> `RadialImage`/`AxialImage`, so `Chord` has no optional fields; slots derived from the crossings, `Crossings.region_entered` the one home of the crossing-order rule; the measure's quadrature over b is the consumer's (geometry may not import `derivations/`, where `chord_quadrature` lives) `[REFUTED 2026-10-05]` the sketch's "integrated by the existing chord_quadrature"; `beam_density` over b >= 0; `parallel` decided exactly, no band.

**Deferred, each with its trigger:**
- the beam density's sphere-area constants collapse onto the measure's table: with E's (r, z) chart (comment on #551);
- caching `Chord.slot_length`/`interface` (each read recomputes `region_containing`): when P1 measures a cost;
- the half-chord product over/underflows for radii beyond 1e+-154 (qa: not needed; a constructor refusal would be the type-level fix): no consumer near those scales;
- the theory page's Cauchy and Santalo citations without equation numbers: #579 (W7);
- the migration of the other spellings (13 chord square roots, 14 discriminants, 8 locators, 3 line measures): #578, which is phase P2.

**Next:** P1's design exchange (the Variant-alpha references posed on `StructuredGeometry` through the kernel; one oracle; the closure rank derived; the attenuated integral as the reference's single implementation; the test roster C4), with the user, before any code. `a336bde4` (the hoist) stays on its branch as the speed target.

## Phases (`[HYPOTHESIS]` 2026-10-05, the main agent; the user rules when the plan is polished)

- **P0, the seed (ruled design).** The kernel in `orpheus/geometry/`: the values, the posed concentric geometry, the orbit-space chord with the obliquity, `Crossing`, typed regions, `region_containing`, the measure on lines per chart. Gates: `scratch/characteristic_architecture/seed_verification_spec.md`, re-measured on the real kernel. Theory: a new section on the geometric kernel (the archivist), and a row in `docs/architecture/conceptual_view.rst`. Done when: the spec's gates are green, each with its first red recorded, and the L4 harness shows the kernel agreeing with today's independent spellings.
- **P1, the reference family on the seed (design still open: C2-C4).** One characteristic reference posed on `StructuredGeometry` (its rays from the kernel, its closure rank derived from the reflecting surfaces a chord meets, the attenuated integral as its own single implementation), replacing the 7 oracle classes and the `geometry_kind` tag; the test roster re-architected (C4). Needs its own design exchange before it is built.
- **P2, the other spellings onto the kernel.** The 13 chord square roots and 14 discriminants, the 8 point locators, the 3 measures on lines; reference code first; production CP and MoC in their own campaigns per the sharpening order (a harmonisation onto shared machinery is admitted).
- **P3, #405 resumes:** the slow reference files, re-timed against `a336bde4`'s measured speed.

## Questions to settle before any design (the survey, dispatched 2026-10-05)

1. The family's full inventory: modules, classes, free functions, facades, result types (`billiard.py`, `power_iteration.py`, `reference.py`, the `greens_function_*` solvers), and which consumers outside the package use what (SN verification references, the traced memo, cross-method tests).
2. What ray and characteristic machinery already exists elsewhere, and in what vocabulary: MoC rays (`orpheus/moc/`), the collision-probability chord geometry (`peierls_nystrom`, `derivations/common/kernels/chord_half_lengths`), the symmetry and orbit-space machinery in `numerics/`, any Volterra operator. Is there already a home for "a ray through a concentric partition"?
3. The tests: what each test file of the family verifies (value against an independent reference, an identity, a facade-delegation tautology), and which are worth keeping under a new architecture.

## Rulings ledger

- **2026-10-05, the user, on C1:** "Is this something like a ray tracing? This seems like part of something that calculates lines, intersections and other things like this. Distance between line to shell, might be similar to distance between centroid and face, and between 2 centroids and the deviation from normal. Seems like the seed of something bigger." Read as: C1 is not a chord helper but the seed of a geometric kernel (Euclidean geometry of points, directions, lines and surfaces), with consumers well beyond the Variant-alpha family. The scope widens before any placement is ruled.
- **2026-10-05, the user, on the kernel seed:** "Yes, this seed": points, directions and lines as values; the concentric geometry as the level sets of `CoordSystem`'s radial coordinate; one quadratic crossing law; the measure on lines; in `geometry/` beside `RigidMotion`, with closed-form gates. General CSG and the unstructured finite-volume metrics are named next members, built when a consumer arrives.
- **2026-10-05, the user, on the boundary-point owner:** asked for the arguments from algorithm behaviour before choosing (see "The boundary-point question").
- **2026-10-05, the user, on the boundary-point question:** "This sounds good enough as a seed. I accept the recommendations." Ruled: segments carry their region from the crossing order; a point on a surface is located with its direction (Omega . grad(rho), grazing to the outer side); producers record the region of the points they place; the one bare-point locator is inner-owns, closed on [0, R], and the 8 spellings retire onto it.
- **2026-10-05, the user, on the revised seed:** "Accept": crossings solved once in the orbit space, the obliquity 1/|P Omega| as the factor the three charts share; directions as points of `numerics.manifold.Sphere`; lines as (direction, moment) values; a `Crossing` produced by the step that hits a surface. This supersedes the first-pass design's crossing law, measure and directional locator (see "W5 review", refuted list).
- **2026-10-05, the user, on placement:** "In the seed": a concentric geometry is posed by a `RigidMotion` from its canonical frame, the identity by default.
- **2026-10-05, the user, on a line lying in an interface:** "Typed on_interface": on_interface(k, k+1), never folded into inner-owns.
- **2026-10-05, the user, on the seed and Architecture E's `Chart`:** "P0 mints Chart (1-D charts)". E's `Chart` (#551; `scratch/boundary_ontology/orbifold_architecture_round2.md`: the pair (coordinate system, kept coordinates), realizing G_c, deriving the orbit map and its measure, the singular strata) is pulled forward for cartesian {x}, cylindrical {r}, spherical {r} only; the seed's crossings, obliquity and measure on lines are its verbs. E later adds (r, z), the deck group and the 26-site switch.
- **2026-10-05, the user, on the exterior codes** (both P0 reviewers: -1/-2 are valid numpy indices, `[M]` a hollow sphere's optical depth 7.0 instead of 5.0): "Out-of-range codes n, n+1": regions 0..n-1, the inner exterior n, the outer exterior n+1; the interface in its own field.
- **2026-10-05, the user, on the chart's form:** "Derive it now": `Chart` = (kept columns, linear group), its verbs derived; the five matches on the coordinate system go.
- **2026-10-06, the user, on P1's first-pass design (D1-D5):** all four recommendations accepted. Q1 "Derived rank": the rank is the number of distinct walls a transit touches (0, 1, 2), with one bounce-sum formula for every rank. Q2 "Dense pencil, LAPACK": the operator assembled once as quadrature weights, k from `scipy.linalg.eig`, no power iteration. Q3 "Eigen + fixed source": the fixed-source question is the same matrices with one linear solve; the MR fixed-source solver retires onto it. Q4 "Accept the roster": the test-architect writes the verification spec before any code.
- **2026-10-06, the user, on references and production numerics:** "numerics linear algebra is made to be highly versatile for production. in principle, it would be interesting is there was a small subset that is appropriate for the methods used for reference, but that insulates production from references. Or else, development in production will cause churn in references." Then ruled "Yes, gate-enforced" and "Precursor before P1". References import only the INTERFACE vocabulary from `numerics` (the question, the observables, the answer function, content identity, the traced memo). Their MATHEMATICS lives in a small reference kernel under `orpheus/derivations/common/`, built on numpy, scipy and mpmath with gates of its own. The layer gate (`tests/gates/test_layer_imports.py`) enforces it with an allowlist. A primitive that exists once per branch (the Perron-Frobenius refusal, for one) is a deliberate duplicate across the branch line, recorded by name in `conceptual_view.rst`. `[M]` 2026-10-06, an AST census of `orpheus/derivations/` (control: `traced_memo` found): 17 `numerics` names in 7 files; 4 are mathematics (`eigenvalue.dominant_eigenpair`, `quadrature.Quadrature`, `moment_layout`, `roots_of_unity`, one file each), the rest interface.
- **2026-10-06, the user, the principle behind the insulation:** "the reference methods (like Fn or trajectory resolvent) are highly closed, purpose built and limited scope. They are fundamentally different than production, which is versatile, generalist and large in scope. So they ask for different things from their machinery and the churn should be concentrated on production, whereas once references reach a good architecture, churn should be extremely limited." Consequences the main agent draws `[R]`: (1) a reference depends on slow-moving upstreams (numpy, scipy, mpmath, its own kernel), never on versatile production machinery; (2) the interface vocabulary it shares with its consumers is itself a churn channel, so the allowlist is kept minimal and a change to it is reviewed as a contract change; (3) the P0 geometric kernel (`orpheus/geometry/{chart,line,chord}.py`), which both branches import, must be held to the reference standard of stability: small, closed, changed only by ruling. The principle gets a durable home at P0.5, in `docs/architecture/layering.rst` beside the layer table.
- **2026-10-06, the user, on the W5 and spec questions (first batch):** boundary laws "Extend now (Krein form)": the closure is the boundary resolvent P = P_0 + E (I - T)^{-1} X on the boundary trace space; specular reduces to D2, white (`IsotropicReturn`) is one W x W block per group, periodic through the laws' deck maps; the rank is the size of T per line, derived from the walls after the deck map identifies them (supersedes D1's "rank = distinct walls" where they differ, the periodic slab). Assembly: "Measure accuracy to see if they differ" (points-and-directions collocation versus Galerkin over lines; measurement dispatched, `scratch/characteristic_architecture/p1_assembly/`). Peierls-Nystrom: "Independent; decide #506 later". Eigen scope: "Add modes and adjoint": higher eigenpairs and the adjoint (importance) question are in P1.
- **2026-10-06, the user, second batch:** a line lying in an interface: "Refuse" (the kernel's on_interface, named in the refusal). The basis for q: "Decide from the measurement" (`p1_assembly/`). The 26 duplicated SymPy rows: "Retire the duplicates" (one kept per identity, its markers migrated). Voids: the user asked "Would this be solved structurally if the law was less specialized than I requested (so instead of deriving from bounces)?"; answer `[R]`: no for specular walls (the general resolvent has the same singular I - T on a line whose period is all void with alpha = 1; the singularity is the problem's), yes for white walls (isotropic return couples every line to the material, so nothing is trapped); the structural fix is to spell the closure as the LEAST NON-NEGATIVE solution, the Neumann series sum_n T^n s, equal to the closed form where it converges and exactly 0 on a trapped line (s = 0 there), the alpha -> 1- limit. Ruled "Least solution": Sigma_t = 0 regions admitted; only an external source on a lossless trapped line (infinite least solution) is refused.
- **Consequences for P0.5 of the eigen scope ruling** `[R]`: the reference kernel gains a full-spectrum primitive beside the dominant eigenpair (higher modes may be complex or sign-changing, so no Perron-Frobenius refusal on them), and one adjoint spelling, W^{-1} P^T W with W the quadrature-weight metric.
- **2026-10-06, the user, on the assembly, after the measurement** (`scratch/characteristic_architecture/p1_assembly/report.md`: Galerkin over lines 28x-175x smaller k error than collocation at equal n, order ~5.5 vs ~3.6, a 1G lower bound; per-region polynomial panels 2x-5x better than the cubic spline and 3x-35x cheaper; pointwise flux better under Galerkin only with graded panels): "Lines + graded panels". The pencil (k, higher modes, adjoint) is the Galerkin assembly over the measure on lines; q lives on per-region polynomial panels graded toward walls and interfaces (#566); a reading at a point is the per-point transport of the converged emission (one transport integral, two test measures). Gated: the symmetry defect of the transport block as the under-integration alarm; the 1G Rayleigh-Ritz lower bound as a theorem row. Supersedes D3's collocation. "Add it": one Sood upscatter fixture (URRb or URRc) joins the independent-value rows.
- **2026-10-06, the user, on the insulation gate's scope:** "Continuous references only". The allowlist binds the closed references, `derivations/continuous/` and `derivations/common/`. Exempt by name, with the reason in the gate: `derivations/discrete/` (algebras of record whose subject is a production discretization; `discrete/sn/balance.py` rides production's `roots_of_unity` on purpose) and the MMS harnesses that pose a production method (`continuous/mms/sn.py`, `continuous/mms/moc.py`). Found while reading the import sites `[M]`: `OperatorPencil`, `EigenPosing` and `SourcePosing` are built on `numerics.operator.LinearOperator` (the production operator algebra), so they are NOT interface; the reference spells its own dense pencil (answers W5 elegance E1 the other way).
- **2026-10-06, the user, on the P0.5 elegance review:** the source solve's refusal of a reducible gain: "Build it now". `DensePencil.least_solution` solves on the unknowns the source reaches (`DensePencil.reach`, the downstream closure of the source's support under the pencil's couplings, the Frobenius normal form's reached classes), zero outside, the radius check on the reached block only; a zero source reaches nothing. The insulation gate's width: "Add mesh and methods": closed references never import `mesh`, `transport` or a method package; `sood_registry/builders.py` joins the named exemptions (it builds production CP problems). Also from that review, by principle: `least_solution` became a method of the pencil (one QZ); `FundamentalMode` checks its own invariant; `Fundamental` and `Mode` renamed `FundamentalMode` and `SpectralMode` (homonyms of `numerics.question` names a reference may import); the production twin's over-claiming docstring filed as #580.
- **2026-10-06, the user, on P1's API sketch:** all four recommendations accepted. White block: "One resolvent, two parts" (K_g = K_line + U (I - T_w)^{-1} A U^T, the diffuse walls as a symmetric finite-rank update of the line-diagonal block). Directions at a point: "Chart verb" (`Chart.directions_at(point)`, the fundamental domain of S^2 / Stab(x) with its density; the stabiliser is D_1h, `SubgroupOfO3.Dnh(1)`, already in the lattice; the S^2/D_1h catalogue entry is NOT built: ruled 2026-10-06 after the measurement that the catalogue's barycentre lift is not a right inverse for D_1h's non-linear chart, filed as #581). Transits: "In the kernel" (`Chord.transits`, gated in the kernel). Adjoint eigenpair: "Derivation API only" (no new observable; minting one is owed to #529; `Response` answered through the door). The rest of the sketch (package `characteristic/`, flux coefficients as the unknown, the migration order) stands as written.
- **2026-10-06, the user, on P1 step (b)'s first rung:** Q1, the tag arm: the user refused hoisting production's tag parse ("mixing tags from reference with production would cause reference churn during production development. but we do need a coherent way to parse tags for reference"); ruled "the registry is accepted": the reference keeps its OWN tag registry (vacuum, reflective, partial, white, periodic) mapping each kind, with the wall's context (axis x, outward sign -1 at breakpoint 0 and +1 at n), to a typed law, and `Walls.of` reads every wall, tag or law, through the factor table only; the registry and production's `_law_from_tag` are a declared cross-branch duplicate (conceptual view); a gate pins Walls(tag) == Walls(its law) per kind. Q2 "By its factors": `PrescribedInflow` with `NoSource` is a vacuum wall, a non-zero source refused. Q3 "Accept as sketched": the period chained through the walls' partners, rank derived; `inflow` on optical depths, the least solution.
- **2026-10-06, the user, on P1 step (b)'s second rung:** all four recommendations accepted. Q1: the ladder is re-cut; rung 2 is the panel basis and the transport along one line (B, A, the Volterra triangle, psi on a line), and rung 3 is the line rule, the assembly and `WallCoupling`, all integrals over the measure on lines. Q2: the pieces come from the kernel's chord through the panel ends posed as a refined `ConcentricPartition`, with `Walls.on(partition)` re-keying the walls. Q3: B_k lives on a new `TraversalRule` in `transport.py`. Q4: the basis is discontinuous nodal panels graded toward walls and interfaces (not a singular stratum), and its volume density is derived inside the basis from `measure_constant` and `measure_coordinate`, gated against `Chart.measure`.
- **2026-10-06, the user, on P1 step (b)'s third rung (after the re-measurement):** all four as recommended. Q1: an even basis at a singular stratum (Lagrange in c² on the panel touching the centre or axis; Schwarz's theorem), not a graded b rule. Q2: the line domain is a kernel verb on `Chart` (the orbit space of lines with its density, the lines' analogue of `directions_at`). Q3: C11 is re-posed as closed-body conservation K·1 = W·1/Σ_t, symmetry kept only as a declared-blind foundation row. Q4: the volume density moves to the kernel, the white block carries D = diag(A_w/4), and S and F move to rung 4.
- **2026-10-07, the user, on P1 step (b)'s fourth rung:** all four as recommended. Q1: the emission support is per group and exact, the regions where that group's emission can be non-zero (a scattering or (n,2n) transfer into it, or fission with chi > 0), widened by the posed source; column sets differ by group. Q2: the cross sections are a frozen `RegionCrossSections` read from one `Mixture` per region through `group_emission`, split out of `_infinite_medium_matrices`, refusing anisotropy there. Q3: rung 4 answers in nodal coefficients; the question values and the projection of a mesh-free source go to rung 5. Q4: the system is `GalerkinSystem`. Later the same day, after the flux-form adjoint was refuted (measured): the pencil's unknown is the emission density q on the supports, not the flux (reverses the P1 sketch's item 6).
- **2026-10-06, the user, on P1 step (b)'s third rung API sketch:** both as recommended. `LineRule.transport` is per group (each chunk's geometry rebuilt per group; a cache only if measured worthwhile). A wall with both a specular and a diffuse part stays refused, as the SN realizer refuses it.
- **2026-10-06, the user, on the cylinder's line coordinate (after the rung-3 spec):** the polar angle θ = arccos μ_z, in the kernel's `LineDomain`, where the density is analytic; `directions_at` keeps μ_z.
- **2026-10-06, the user, after rung 3's elegance review:** `compute_areas_1d` (the twin of `measure_density`) retires separately, #584; the white block's near-void conditioning is fixed in rung 3 (the loss-formed I − αT with a balance row).
- **2026-10-06, the user, after rung 3's qa:** every line-measure grading is derived from the group's optical scale (grazing to the thinnest absorbing panel; impact in the chord half-length, hp toward both singularities and exponential at the rim); `LineRule` becomes per group. The regimes are not refused.
- **2026-10-07, the user, on the cylinder's cost:** rung 3 lands with the slow cylinder rows cut to one fixture per law at 8 points; #586 (a non-tensor (b, θ) rule) is the next step, before rung 4.
- **2026-10-07, the user, on P1 step (b)'s fifth rung API sketch:** all four as recommended. Q1: split, 5a (the resolution, the projection, the door, `Eigenvalue` and `FluxIntegral`; `PointValue` refused naming 5b) then 5b (the reading at a point). Q2: a `Response` is answered as the group-transposed forward problem with the detector as its source. Q3: the eigen flux in the gauge ⟨νΣ_f, φ⟩ = 100. Q4: `Nearest(tau)` served, tau read in the k chart, #529 named as the owner of parameter charts.
- **2026-10-07, the user, after rung 5a's reviews (three rulings):** (1) a `Nearest` answer reads `Eigenvalue` only; the flux of a higher mode is refused as a SCOPE-BOUNDARY, since its net fission production can vanish (on a closed homogeneous body every higher mode is biorthogonal to the flat adjoint). (2) The Q3 gauge is corrected: its premise was wrong (100 is the homogeneous solver's production DENSITY); the eigen flux of a finite body has total fission production 1 over the body, the SN solver's `ScaleGauge(production rate, 1.0)` and the trajectory resolvent's. (3) A `Response` answer is the adjoint scalar flux Rψ† (the vocabulary's E^{-†}R), so the transposed problem's source is the retraction of the detector's lift: R R†Σ_d = 4πΣ_d for a table, R f for a symbolic function, the 4π derived from `angular_measure`.
- **2026-10-08, the user, on the gauge (after the archivist found that SN's production counts the (n,2n) emission):** "In principle, both should know what 'production' means (because we need to declare it) so that we can have n2n or not and be consistent." Ruled from the two options that followed: the gauge is DECLARED on the eigen question, `Eigen(parameter, ..., gauge=<CellCoefficient>)`, the same key type as its parameter, resolved by the specification; `None` takes the declared default `EIGEN_GAUGE = CellCoefficient.every(FISSION_EMISSION, N2N_EMISSION)` (production: what fission and (n,2n) emit, SN's functional). What each channel emits per unit flux is defined once, `Channel.emission` (ν per fission, 2 per (n,2n) reaction). The references read it in this rung; SN and the homogeneous solver move onto it in #517.
- **2026-10-08, the user, on rung 5b's API sketch:** all four as recommended. `PointValue` reads the transported emission 𝒦q at the point (iterated Galerkin), not the basis value; the point's directions are the line rule's own lines read at the point's parameters, from the same constructors with the point's c as one more end, and no new `Resolution` field; ψ(x, Ω) lands in 5b, exposed for gates with no observable; `FluxIntegral` keeps the pairing.
- **2026-10-08, the user, after rung 5b's reviews (four rulings):** the general grading law at every impact-panel top lands in 5b; the half-chord is carried in 5b (#590 done in this rung); its kernel shape is an image that carries an exact level (`RadialImage.level`, `level_half_chord`, default today's arithmetic, passed by `ConcentricPartition.chord(line, level=...)`), not a field on `Line`; the cylinder's extra cost is accepted until #587 grades per polar angle.
- **2026-10-05, the user, on sequencing:** design first, then rebuild; `a336bde4` stays unmerged on its branch as a measured speed target.

## The widened question (2026-10-05)

`[HYPOTHESIS]` (the main agent) One structure underlies:
- the distance along a ray p + t*Omega to a surface f(x) = 0, a root of f(p + t*Omega) (a quadratic for the quadrics reactor geometry uses: planes, cylinders, spheres); the chord through concentric shells is this over nested circles;
- the distance from a point to a surface, and its normal grad f;
- specular reflection of a direction about that normal (a Householder map: the boundary laws' reflection, #436's face pairing);
- mesh metrics: cell centroids (first moment over volume), the centroid-to-face distance, the centroid-to-centroid vector and its deviation from the face normal (finite-volume non-orthogonality);
- the point/direction distinction, an affine torsor (`coding-standards`, "Type vs property").
Candidate consumers: MoC tracks, Monte Carlo's distance to boundary, CP chords, the Variant-alpha and Peierls references, the finite-volume metrics of diffusion and SN, the reflective boundary laws. Unmeasured; the census below measures it.

## Log

- 2026-10-05: scope widened by the user's ruling on C1; census of geometric computation dispatched.

- 2026-10-05: plan opened on the user's ruling above.

## ⏸ COMPACTION POINT — 2026-10-06, P0 merged, P1's design exchange next

**State:** `main` is `3db74d59`, pushed. It holds P0: the kernel (`fe0ca696`) and the theory page (`82ae7013`, CI `gates` run 37394764313 green), plus this plan and agent memory (`3db74d59`, `[skip ci]`, pushed separately). The not-slow suite at `82ae7013`: 15 389 passed, 0 failed (`[M]` 2026-10-06, detached worktree). The branch `refactor/chord-oracle-axial-lift` (`a336bde4`, the bit-identical hoist of the old oracles) is parked unmerged as P1's speed target; nothing else is open. The tree is clean apart from `scratch/`.

**Read in order:**
1. this file's "Rulings ledger" (every ruling, dated);
2. "The seed, revised", "P0 API sketch" and "P0 landed on its branch" (what exists, what was refuted, the deferrals);
3. the inventory section and "Candidate architecture" C2-C4 (the reference family P1 rebuilds; `scratch/characteristic_architecture/inventory.md` holds the full inventory and the test classification);
4. the kernel itself: `orpheus/geometry/chart.py`, `chord.py`, `line.py`, and `docs/theory/foundations/chart_and_chord.rst`.

**The next step, P1 (design exchange with the user before any code; the ontology of the reference is still being searched):** the Variant-alpha references posed on `StructuredGeometry` through the kernel. Open questions to bring to the user:
- one reference for all geometries, its rays from `ConcentricPartition.chord` (`lengths_beyond` for the backward first leg; the full chord for the bounce period), replacing the 7 oracle classes and the 8-case `geometry_kind` tag;
- the closure rank derived from the reflecting surfaces a chord meets (rank 1: solid sphere and cylinder; rank 2: shells and the slab), where `billiard.py:299` tags rank 2 for the slab only;
- the attenuated line integral (the reference form of the Volterra operator) as the reference's single implementation, reading the density per region from the chord's slots (exterior codes out of range: a cavity gets a void entry on purpose), with the axial factor through `projected_speed`;
- the test roster (C4): keep the independent-value rows (Sood, WM-72, PS-1982, Garcia, case-truth, the first-leg line integral), retire the 17 delegation tautologies and the twin-class comparisons, re-pose "closed body gives k = k_inf" (holds for any geometry) and the 87 SymPy identities by what they verify;
- performance: the hoist `a336bde4` measured the cylinder A|B|A brute call 1.3 s -> 0.073 s; the rebuild must meet it (the kernel's batch shape is built for it); the slow reference files (#405) are re-timed at P3.

**Lessons from this stretch:**
- A refusal written "any departure > tol" admits NaN; write "all within tol" (four NaN/overflow defects in one kernel, each found by a reviewer's hostile input, not by the closed-form gates).
- Negative sentinel codes are valid numpy indices: an exterior code must be out of range so indexing a per-region table raises.
- A value solved in a canonical frame must report parameters on the caller's object; the identity pose hides the defect.
- At a subnormal scale, hold quantities as finite numerators over one division, so overflow gives a signed infinity and never inf - inf or 0 * inf.

## P1, first-pass design (the main agent, 2026-10-06; for the user's rulings, nothing built)

Every item is `[HYPOTHESIS]` until ruled. It answers the five open questions of the compaction point above.

**D1. The billiard lives in the orbit space, and its rank is derived.** Specular reflection on a surface of the chart's symmetry group maps a line to a congruent line: on the sphere and the cylinder it keeps the impact parameter b and |P Omega| (the cylinder also keeps Omega_z); on the slab it reverses Omega_x. In the orbit space this is a one-dimensional motion of the orbit coordinate c(t) between walls (the geometry's boundary points, each with its law) and, on the radial charts, one smooth turning point at c = b. Define a *transit* as the maximal run of the line inside the domain [r_0, r_n], read off the chord's crossings: it starts at one wall and ends at a wall. Reflection at the end wall continues with the transit reversed. So the period of the unfolded backward characteristic is {transit, reversed transit}, and the closure's rank is the number of distinct walls the transit touches:
- rank 0: a line that reaches no wall (the kernel's `parallel` line, a cylinder ray along the axis or a slab ray with Omega_x = 0); psi is the integral over an infinite path;
- rank 1: both ends on one wall (a solid sphere or cylinder; a shell's ray with b > r_0, which misses the cavity);
- rank 2: the ends on two walls (the slab; a shell's ray with b < r_0).
`[M]` by reading: `billiard.py:299` tags rank 2 for the asymmetric slab only, while the shells branch on b inside their oracles. Under D1 the shells, the slab and the solid bodies are one construction, and `reference_body`'s four shapes are not needed by this family.

**D2. One closure formula for every rank: the geometric sum over the period.** Backward from (x, Omega): the first leg F (from x to the first wall), then the period's transits j = 1..m, each with its source integral B_j (the attenuated integral to the transit's end), its optical depth tau_j and the albedo alpha_j of the wall it ends on. Then psi = F + e^{-tau_F} (sum_{k=1..m} prod_{j<k}(alpha_j e^{-tau_j}) alpha_k B_k) / (1 - prod_j alpha_j e^{-tau_j}). The rank-1 and rank-2 resolvents of `variant_alpha_core.py` are its m = 1 and m = 2 cases; m = 0 is F alone. This is the method of images written for a periodic unfolded path.

**D3. The reference's operator is a discrete measure per evaluation point (the characteristic quadrature).** For a fixed geometry, laws and Sigma_t, psi(x, Omega) is a linear functional of the emission density q: quadrature points t_p on the unfolded path's slots, with weights carrying the attenuation and the closure, give psi = sum_p W_p q(c_p, region_p), the region known from the slot (never located, ruling of 2026-10-05). Integrating over Omega gives the scalar flux at x as a discrete measure on the orbit space. Two consequences:
- the geometry is computed once and reused for every source and every iteration (`[M]` by reading, the inventory's data flow: today every power step rebuilds every chord, per group); the hoist's speed came from the same principle at a smaller scale, one in-plane chord shared across the axial cosines (`[R]` that D3 meets or beats it; measured at build time);
- with q represented on a regionwise basis (today's regionwise cubic spline on composite Gauss nodes), the operator is a dense matrix P_g per group, and the k question is the small dense pencil phi = P (S phi + F phi / k), solved directly by `scipy.linalg.eig` (trusted upstream, below the library line of `algebra-of-record`): no iteration, no tolerance, no `converged` flag.
A reading at any point builds that point's measure and applies it to the converged q: one transport, as today.

**D4. One reference, posed from the specification.** Input: the `StructuredGeometry` (through `ConcentricPartition.of`) with its boundary laws (through `specular_albedos`), the `Materials`, the question, and the angular rule over the directions at a point (reduced by the point's stabiliser, split at the tangency directions b = r_k). The 7 oracle classes, the 15 `solve_greens_function_*` entry points, the 8-case `geometry_kind`, `Billiard`'s routing and `reference._RAYS` retire onto it. The cylinder is not a case: the 3-D line through the partition carries the axial factor through `projected_speed`.

**D5. The test roster (C4).** Keep: the independent-value rows (Sood, WM-72, PS-1982, Garcia, the case-truth rows, the three production-primitive closed forms), the F_N cross rows, the first-leg line-integral row, the k_inf rows as a conservation floor, the reading laws, the refinement ladders, the SN gates that use the family as their reference. New: D1's rank against hand-counted walls; D2 against the explicit sum of the unfolded path truncated far out (an independent spelling); the characteristic quadrature against mpmath on a manufactured q. Retire with the code: the 17 delegation tautologies, the 6 one-oracle-two-drivers rows, the twin-class comparisons. The 87 SymPy identities stay where they verify a law the new code realises, and each is re-pointed at the new code or marked as verifying the theory page only.

**Open for the user:** D1-D2 (derived rank, one formula); D3 (the pencil solved directly versus keeping a power iteration, and whether the reference may use the project's own `power_iteration`, which would share the solver axis with the SN gates it judges: X4); D4's scope (fixed-source questions too, as the MR fixed-source solver does today); D5.

**Ruled 2026-10-06** (see the ledger): D1-D5 as written. Next: a W5 review of D1-D5 (cross-domain-attacker, elegance-enforcer) and the test-architect's verification spec, in parallel; then the API sketch, a checkpoint with the user, then code (P1 is a surgical carve: the main agent writes it).

**P0.5, the reference numerics kernel (ruled 2026-10-06, precursor to P1).** One small carve: the kernel in `orpheus/derivations/common/` (dense dominant eigenpair with the real-eigenpair, single-signed refusal; dense fixed-source solve refusing a spectral radius >= 1; the composite Gauss-Legendre rule already in `quadrature.py`), each with its gates and first red; the 4 mathematical imports moved onto it (or onto its interface, after reading each); the allowlist in the layer gate, its first red today's 4 imports; the conceptual-view row naming the deliberate cross-branch duplicates. P1's D3 then reads "the reference kernel's dense eigenpair", not `numerics.eigenvalue`, and the W5 elegance finding E1 is answered for the question types only (`OperatorPencil`, `EigenPosing`, `SourcePosing` count as interface if they compute nothing: to be read at P0.5).

## P0.5 API sketch (the main agent, 2026-10-06; checkpoint with the user before code)

**Module** `orpheus/derivations/common/dense_pencil.py`, imports numpy and scipy only.

- `DensePencil(loss, production)`: two square float arrays in WEAK form, i.e. the Galerkin bilinear forms a(u_i, u_j) of the loss L and the production F on one nodal basis. The k question is `production c = k loss c`. Refuses non-square, mismatched or non-finite input. Weak form is the representation convention, fixed here once: the adjoint is then the transpose, so no metric is needed for it.
  - `fundamental() -> Eigenpair`: the dominant eigenpair by `scipy.linalg.eig`. It refuses when the dominant eigenvalue is not real, not simple or not positive, or when its eigenvector is not single-signed. The refusal is the Perron-Frobenius / Krein-Rutman contract; nodal coefficients are values, so single-signed coefficients mean a single-signed function.
  - `spectrum() -> tuple[Eigenpair, ...]`: every eigenpair, by decreasing |k|; complex and sign-changing modes allowed (no refusal).
  - `adjoint() -> DensePencil`: `DensePencil(loss.T, production.T)`.
  - `source_solution(source) -> ndarray`: the fixed-source answer `(loss - production) c = source` as the LEAST non-negative solution. It exists iff the fundamental k < 1 (subcritical); otherwise it refuses with the measured k. A pure scattering problem is the same with `production` the scattering form.
- `Eigenpair(k, vector)`: frozen; `vector` normalised by one stated rule (the sign rule `sum >= 0`, scale left to the caller, as the production twin does).
- The production twin `numerics.eigenvalue.dominant_eigenpair` stays where it is: the duplicate across the branch line is deliberate (ruling 2026-10-06), recorded by name in `conceptual_view.rst`.

**Moves.** `derivations/common/eigenvalue.py:226` (the k_inf adjoint spectrum) onto `DensePencil` (the 0-D problem is the 1-node case). The other three mathematical imports are exempt by the scope ruling.

**Gate** in `tests/gates/test_layer_imports.py`: `REFERENCE_INTERFACE_NUMERICS = {traced_memo, content, mesh_free_function, observable, question}` (the measured interface set; a change to it is a contract change); a modules-under-`derivations/continuous` and `derivations/common` rule allowing only those `numerics` modules; the named exemptions. First red: today's `common/eigenvalue.py:226`. Positive control: the census finds the 4 known mathematical imports.

**Docs** (the archivist): the principle and the gate on `docs/architecture/layering.rst`; the conceptual-view row (the reference kernel, its deliberate twins); the kernel's derivation (Perron-Frobenius and Krein-Rutman for the refusal, the least solution and its existence iff subcritical, weak form and the transposed adjoint) on the reference-solutions theory page.

**Gates** (the test-architect's spec, P0.5 section): K1/K1r, K2/K2r, K3, K4, plus the spectrum and the adjoint's biorthogonality (forward and adjoint modes F-orthogonal).

## P0.5 landed on its branch (2026-10-06)

**Commits on `feature/reference-dense-pencil`:** `56328fbd` (the kernel `orpheus/derivations/common/dense_pencil.py`, the moves, the insulation gate, `tests/gates/derivations/test_reference_kernel.py`) and `a5ba5c81` (`layering.rst`, `conceptual_view.rst`, the kernel's theory in `reference_solutions.rst`, the API entry, the regenerated matrix). Merge status: read git, not this line.

**What it is, as built:** `DensePencil(loss, production)` in weak form; `fundamental()` -> `FundamentalMode` (refuses complex, non-positive, not strictly dominant beyond `DOMINANCE_GAP`, or not single-signed beyond the pencil's own `sign_band()`, the eigenvector perturbation bound times `SIGN_SAFETY`); `spectrum()` -> `SpectralMode`s; `adjoint()` the transposed pair; `reach(source)` and `least_solution(source)` (the Neumann sum on the reached block, refusing a reached radius not below `1 - SUBCRITICAL_MARGIN`). Moved onto it: the k_inf forward and adjoint spectra, `kinf_from_cp`, the two diffusion buckling k sites. The gate: closed references import from `numerics` only `question`, `observable`, `mesh_free_function`, `content`, `traced_memo`, and never `mesh`, `transport` or a method package; exempt by name: `mms/sn.py`, `mms/moc.py`, `sood_registry/builders.py`; `derivations/discrete/` out of scope.

**Evidence** (`[M]` 2026-10-06): kernel gates 33 passed under `-O`, each with its first red (test-architect battery `scratch/characteristic_architecture/p05_gates/README.md`; the main agent's in-process reds for the margin, the reach and the sign band); `test_layer_imports.py` 470 passed; the 38 k_inf consumer files 733 passed, 0 failed; qa's old-versus-new extraction wrapper over 902 tests: 252 calls, 0 refusals, largest k difference 2.2e-15; Sphinx `-E -W` 0 warnings; `dead_references` 0 of 66; pyright 0 errors on the changed modules.

**Design changes during the carve, by review or ruling:** `least_solution` became a method (one QZ); the exact reducible case built (the user: "Build it now"); the gate widened to mesh and methods (the user); `Fundamental`/`Mode` renamed `FundamentalMode`/`SpectralMode` (homonyms of `numerics.question`); the subcritical margin (a bare comparison with 1 admitted 473 of 4000 unit-radius draws); the sign band derived from conditioning (a fixed 1e-10 refused 9 of 20 pencils at cond 1e8); complex input refused; the gate reads imported names (`from orpheus import numerics` had passed).

**Filed:** #580 (the production twin's docstring over-claims its contract).

**Deferred, each with its trigger:** a nodal basis is a precondition of the single-sign test (elegance C3): P1's panel basis must be nodal, or the test reads the mode's values through the basis's evaluation map; the 7 trajectory-resolvent k_inf-guess sites retire with P1's rebuild; Krein-Rutman and Perron-Frobenius theorem numbers are not on disk (a W7 task).

**Next:** P1's API sketch, a checkpoint with the user, then code (surgical carve).

## ⏸ COMPACTION POINT — 2026-10-06 (second), P0.5 merged, P1's API sketch next

**State:** `main` is `d00a441b`, pushed. It holds P0 (the geometric kernel, `82ae7013`), P0.5 (the reference numerics kernel `56328fbd`, its docs `a5ba5c81`, CI `gates` run 37416616383 green), and the plan and memory commits. Full not-slow suite at `56328fbd`: 15 446 passed; the 4 failures are the `test_write_guards` worktree artefact. No branch is open apart from the parked hoist `refactor/chord-oracle-axial-lift` (`a336bde4`, P1's speed target: the cylinder A|B|A brute call 1.3 s -> 0.073 s). The tree is clean apart from `scratch/`.

**Read in order:**
1. the "Rulings ledger" entries dated 2026-10-06 (P1's design D1-D5 and every later ruling: the Krein boundary resolvent, Galerkin over lines with graded panels, the least solution, modes and adjoint in scope, voids admitted, a line in an interface refused, the 26 duplicated SymPy rows retired, references insulated from production, the gate's scope and width, the exact reducible least solution);
2. "P1, first-pass design" (D1-D5; D1's rank and D3's collocation are superseded by the rulings);
3. "P0.5 landed" (what the kernel is, its deferrals);
4. the verification spec `scratch/characteristic_architecture/p1_verification_spec.md` (49 gates, the 298-row roster KEEP 122 / RE-POSE 99 / RETIRE 77, the SN consumers' carry-over rule, the performance protocol, X4); the assembly measurement `scratch/characteristic_architecture/p1_assembly/report.md` and its prototype `proto.py`, `cyl.py` (a working Galerkin-over-lines assembly on the sphere, slab and cylinder, numpy and scipy only); the W5 reports `p1_w5_attacker.md`, `p1_w5_elegance.md`;
5. the kernels: `orpheus/geometry/{chart,line,chord}.py` and `orpheus/derivations/common/dense_pencil.py`, with `docs/theory/foundations/chart_and_chord.rst` and the "reference kernel" section of `docs/theory/verification/reference_solutions.rst`.

**The next step, P1's API sketch (a checkpoint with the user before code; P1 is a surgical carve, the main agent writes it).** The objects the rulings imply, each to be named and placed in the sketch:
- the transit: a maximal in-domain run of a line, read from `ConcentricPartition.chord`'s crossings (walls = crossings at the boundary points, keyed by breakpoint index, never tuple position; the albedo paired with the wall the BACKWARD path reflects at, `[M]` the elegance probe: the other pairing moves psi by 0.34);
- the closure as the boundary resolvent P = P_0 + E (I - T)^{-1} X: line-diagonal for specular and partial specular (the period sum over {transit, reversed transit}, the least solution, exactly 0 on a lossless trapped line); one W x W block per group for a white wall (the wall's current couples every line); periodic through the laws' deck maps (rank 1, no reversal, on the periodic slab);
- the basis: per-region nodal (Lagrange) panels graded toward walls and interfaces (#566), nodal because `DensePencil.fundamental`'s single-sign test reads coefficients as values;
- the Galerkin assembly over the measure on lines (the chart's `beam_density`; b split at the radii with the cosine endpoint map; about n + 8 inner points per piece; exponential arc-length panels on thick cylinder chords), its symmetry defect gated as the under-integration alarm, the 1G Rayleigh-Ritz lower bound as a theorem row;
- the questions on `DensePencil`: the k pencil (W - P S, P F), higher modes, the adjoint pencil, the fixed source as `least_solution`;
- the reading at a point: the per-point transport of the converged emission (directions split at tangencies with the cosine map and grazing grading near a surface);
- one reference posed from the `GeometrySpecification`, replacing the 7 oracle classes, the 15 `solve_greens_function_*` entries, the 8-case `geometry_kind`, `Billiard`'s routing and `reference._RAYS`; its home under `orpheus/derivations/continuous/` (a new package; name to be proposed), inside the insulation gate.
Open for the sketch to settle: where the white wall's W x W block lives relative to the line-diagonal specular closure (one boundary-resolvent object with two realizations, or the laws' own responses); the cylinder's direction measure at a point (S^2 / Stab(x), absent from the orbit catalogue, elegance E6); the migration order (build beside the old family, corroborate with the 7 temporary L4 rows, re-point the 13 SN rows, then retire).

**Lessons from this stretch:**
- A refusal comparing a computed quantity with an exact threshold is undecidable at the threshold: give it a measured margin (473 of 4000 unit-radius draws were admitted by a bare comparison with 1).
- A rounding band is a property of the problem, not a constant: derive it from the perturbation bound (a fixed 1e-10 refused 9 of 20 pencils at cond 1e8).
- An import gate that reads only module paths misses `from package import submodule`: read the imported names.
- Two types in different layers that share a name collide in the first module that needs both: grep the interface vocabulary before naming (`Fundamental`, `Mode`).
- Measure a design choice the user asks to see measured before arguing it: the assembly question turned on a 28x-175x accuracy difference no review predicted.

## P1 API sketch (the main agent, 2026-10-06; checkpoint with the user before code)

Every item is a proposal until ruled. It names and places the objects the rulings imply (the compaction point's list), settles the three open questions with a recommendation each, and corrects one statement of the verification spec. Read beside `p1_assembly/proto.py`, which realises most of it as scratch code.

**Package.** `orpheus/derivations/continuous/characteristic/`, inside the insulation gate: imports numpy, scipy, the P0 kernel (`orpheus/geometry/{chart,line,chord}.py`), the boundary-law declarations (`orpheus/geometry/boundary`), `derivations/common/{dense_pencil,quadrature,eigenvalue}` and the interface vocabulary only. The name says what the method integrates along (the characteristic lines); `trajectory_resolvent/` keeps its name until it retires, so the two coexist during the migration without a homonym.

**1. Transits (pure geometry, in the kernel).** `Chord.transits -> Transits`: per line, the maximal runs of slots with an interior region code (< n), each with the breakpoint index of the wall at either end (`[M]` elegance E5's probe: no branch on b against r_0, no inner/outer tag). A wall is a crossing at a boundary point, keyed by its breakpoint index, never by tuple position. Rank-0 lines (the kernel's `parallel`) have no transit. A line lying in an interface is refused by the reference, naming `on_interface` (ruled). This adds one verb to the P0 kernel, which the stability principle says changes only by ruling: this sketch asks for that ruling.

**2. Walls (Branch 1, `walls.py`).** `Wall(breakpoint, specular, diffuse, deck)`: the specular and diffuse return amplitudes and the deck (`Mirror` or `Wrap(partner)`), read once from `geometry.boundaries` zipped with `geometry.boundary_points` by matching the law's own factors (`law.geometry_map`, `law.response_kernel`: `ScalarResponse` under a mirror deck is a specular return of that amplitude; `SpecularReemission` is specular; `LambertianReemission` is diffuse; `PairedDeck` is a wrap). It replaces `specular_albedos` for this family and serves 5 of the 6 transport laws (the attacker's census, `ZeroFluxBoundary` being diffusion's); the sixth, `PrescribedInflow`, needs a boundary-source question the vocabulary does not have (refused, SCOPE-BOUNDARY). A `LawSum` of a specular and a diffuse law reads as both amplitudes on one wall; whether P1 serves it or refuses it is a gate decision at build (B-rows exist for each part alone).

**3. The boundary resolvent (Branch 1, `closure.py`)**, P = P_0 + E (I - T)^{-1} X, realised in two parts:
- `LinePeriod`: per line, the unfolded backward path as the ordered directed transits of one period, each paired with the wall the BACKWARD path reflects at, its specular amplitude and its deck (a mirror reverses the next transit, a wrap continues it in the same direction: the periodic slab is rank 1 with no reversal). The per-direction functional is D2's sum over that period, written as the least solution: the Neumann series, equal to the closed form where it converges and exactly 0 on a trapped line with no source; a lossless trapped line carrying a source is the only refusal. Specular, partial specular, vacuum (amplitude 0, bitwise the vacuum operator, B4a) and periodic are all this part.
- `WallCoupling`: the diffuse part, one W x W block per group (W the walls with a diffuse amplitude), built from the same lines with every diffuse wall treated as absorbing in the line part: the escape functional U (emission to outgoing partial current at each diffuse wall) and the wall-to-wall transmission T_w. Recommendation for the open question "where the white block lives": the transport block is K_g = K_line + U (I - T_w)^{-1} A U^T, a symmetric finite-rank update of the line-diagonal block (A the diffuse amplitudes; reciprocity makes the re-entry functional U^T), so one object assembles both parts and C11's symmetry alarm covers the white law too. The alternative, the laws' own realised responses, is the production route and is refused by the insulation principle.

**4. The basis (Branch 1, `basis.py`).** `PanelBasis(geometry, degree, grading)`: per region, polynomial panels of degree p with nodal Lagrange functions through each panel's Gauss-Legendre points; panel ends graded geometrically toward every wall and interface (ratio and layers in the resolution); no node on an interface. It is nodal, so `DensePencil.fundamental`'s single-sign test reads coefficients as values. It records each node's region and panel, and gives the mass matrix W in the chart's measure (`Chart.measure`'s density, exact by Gauss on each panel). The quadrature is `composite_gauss_legendre` (`derivations/common/quadrature.py:332`); the trajectory family's `_composite_per_region_gl` and `_gauss_legendre_pieces` retire onto it with the family.

**5. The Galerkin assembly over lines (Branch 1, `assembly.py`).** `LineRule(partition, basis, resolution)`: the measure on lines through `Chart.beam_density` (2 pi b db on the sphere, 2|P Omega| per unit height on the cylinder with the polar angle, |Omega_x| on the slab), split at every radius and every panel end, each piece under the cosine endpoint map (`chord_quadrature`'s substitution); along each line, the double arc-length integral of u_i G_L u_j per slot and panel, with about p + 8 points per piece and exponential arc-length panels on optically thick slots. `transport_block(g) -> ndarray` (K_g; the sketch's name, renamed: the landed function, 2026-10-15, is `group_transport(basis, walls, sigma_t, resolution, support) -> GroupTransport`, step 2c); the scattering and fission blocks per node from `_infinite_medium_matrices` per region (so (n,2n) enters as 2 Sigma_2 in S and its refusal retires, elegance E8). Gates: C11's symmetry defect as the under-integration alarm; D11's 1G Rayleigh-Ritz bound.

**6. The questions, on `DensePencil`.** The unknown is the FLUX coefficients phi on every panel; K's columns span only the emission support (regions with scattering, fission or an external source), so a void region never gets an emission column and a trapped line with no source never enters K. Correction to the verification spec §0, which writes the pencil "on the emission coefficients c": with phi the unknown,
- k: `DensePencil(W - K S, K F)`: `.fundamental()` for `Eigen(..., mode=Fundamental())`; `.spectrum()` and the mode nearest tau for `Nearest(tau)`;
- the adjoint: `.adjoint()`, the transposed pair, is the weak-form adjoint in the same nodal coordinates because W and every K_g are symmetric (reciprocity) and S, F act per node;
- fixed source: `DensePencil(W, K (S + F)).least_solution(K q)`, the composition the production half spells `SourcePosing(pencil.at(1), q)`; supercritical refused by the kernel;
- response: `DensePencil(W, K (S + F)).adjoint().least_solution(W r)` for a detector r, the response being y^T K q.
A FluxIntegral of a weight in the basis space is a pairing with the converged answer, w^T K (emission), with no point loop.

**7. The reading at a point (Branch 1, `reading.py`).** `point_flux(x, emission)`: the per-point transport of the converged emission through the same `LinePeriod` and `WallCoupling` (one transport, two test measures, C6), over a direction rule at x split where the kernel's crossing set changes (the tangencies), under the cosine map, with grazing grading near a surface. `psi(x, Omega)` is public (ruled). Open question, the cylinder's direction measure at a point (elegance E6): recommendation, the chart gains `directions_at(point)`, the fundamental domain of S^2 / Stab(x) with its density (sphere: the cosine to the radius on [-1, 1]; cylinder: the axial cosine on [0, 1] times the in-plane angle on [0, pi], the two-mirror stabiliser; slab: the cosine on [-1, 1]), a kernel verb derived from the chart's group, with the two-mirror subgroup added to the `numerics.manifold` catalogue on the production side; the reference builds its own split rule over that domain and never imports numerics. Alternative: the reference states the three domains itself (a per-chart table in Branch 1).

**8. The reference (Branch 1, `reference.py`).** `CharacteristicDerivation(ContentIdentity)`: fields `specification` (a `GeometrySpecification`) and `resolution` (a frozen `Resolution`: degree, grading, line-rule and arc-length points per piece, direction points per piece; one value for the solve and the reading, which dissolves #516's two answers). Construction is the one door: it refuses an infinite medium, anisotropic scattering, an `Eigen` off the physical point or along another parameter than the fission emission, `PrescribedInflow`, a line in an interface; each refusal is a SCOPE-BOUNDARY with its trigger. The solve is a `cached_property`; `evaluate(observable)` is the traced memo, answering `Eigenvalue`, `PointValue`, `FluxIntegral` (and `Ratio` by composition) for the eigen, fixed-source and response questions, uncertified until P4 (#566). `characteristic_reference(specification, resolution) -> ReferenceSolution`, the same consumer surface as `trajectory_resolvent_reference`. Open question, the adjoint eigenpair: the question vocabulary says it belongs to the eigen answer and no observable reads it today; recommendation, P1 exposes it on the derivation's answer for its gates (D13, D14) and mints no observable (that is a contract change to the interface vocabulary, owed to the posing sequence #529).

**9. Migration order.** (a) The kernel verb (`Chord.transits`, `Chart.directions_at` if ruled), with gates. (b) The package, built beside the old family, rung by rung on the spec's ladder: walls and closure (A, B), basis and assembly (C), questions (D, E), reading and door (C1-C2, F). (c) The corroboration file (spec §7, G): the old family against the new on the SN fixtures, temporary L4 rows. (d) The 13 SN rows re-pointed; the two reference keys of `CYLINDER_3REG_RECORD` re-baselined, the three SN keys required unmoved. (e) The old family retired in one commit with its dependency-audit table (`retirement-audit` G.24): the 7 oracle classes, `Billiard`, the 15 `solve_greens_function_*`, `variant_alpha_core`, the family's `power_iteration`, `reference._RAYS`, the 7 k_inf-guess sites, the 77 RETIRE rows; the 26 duplicated SymPy rows retire with their markers migrated; the G rows are deleted in the same commit. (f) The theory page rewritten (the archivist). Each step a checkpoint with the user; P1 is a surgical carve, the main agent writes it.

**Performance.** Acceptance per spec §8: the hoist's cylinder `ABA` brute call (`a336bde4`, 0.073 s recorded, re-measured min of 15 interleaved repeats) is the target for one point reading; the assembly is paid once per specification and resolution (the traced memo).

**Ruled 2026-10-06** (see the ledger): the sketch as written, with the four recommendations. Next: migration step (a), the two kernel verbs with their gates, then the package rung by rung.

## P1 step (a) landed (2026-10-06)

**Commits:** `2b2d7703` (the two kernel verbs and their gates) and `cafdb163` (theory sections, labels `geometry-transits` and `geometry-directions-at`, markers moved to l0). Merge status: read git, not this line.

**What it is, as built:** `Chord.transits` -> `Transits` (at most two runs per line; walls as breakpoint indices; absent transits carry the empty slot range at the chord's end and the wall code n + 1). `Chart.directions_at(c)` -> `DirectionDomain(chart, orbit_coordinate)`, its shape a `DirectionShape` (WHOLE, COSINE, AXIAL_COSINE, ANGLE_AXIAL), deciding stabiliser, axes, bounds, representatives and tangencies by an exhaustive match; the impact parameter is `Chart.image`'s; tangencies scale-free; out-of-box coordinates refused. The shape table is a SCOPE-BOUNDARY retiring onto the S^2/D_1h catalogue entry.

**Rulings during the step:** the S^2/D_1h catalogue entry is not built (the user, after the measurement that the catalogue's barycentre lift is not a right inverse on D_1h's non-linear chart): #581. tangencies(level) includes the grazing value at level = c.

**Evidence** (`[M]` 2026-10-06): gates 168 rows (89 transits, 79 directions), each first red run (test-architect battery 18 arms, `scratch/characteristic_architecture/p1_step_a/battery/`; qa's 3 added rows, `p1_step_a/main/qa_rows_reds.py`); geometry and layer gates 1937 passed; Sphinx `-E -W` clean; `dead_references` 0 of 66; pyright 0 errors.

**Filed:** #581 (catalogue lift on non-linear charts), #582 (the kernel's half-chords at extreme radii).

**Measured, for the reading at a point:** within about sqrt(eps) of grazing at the point's own level no double-precision b resolves the side of the surface (r - b ~ r delta^2 / 2); the lost chord length is about 2 r delta / |P Omega|. The reading's grazing grading must not rely on resolving that band.

**Follow-ups from review, not done (nits):** the absent-wall code reuses the value of `no_interface` (a homonym); `_MAX_TRANSITS = 2` rests on the chord's slot layout, not derived from it.

**Next:** migration step (b), the package `orpheus/derivations/continuous/characteristic/`, rung by rung: walls and the closure (spec rows A and B) first.

## ⏸ COMPACTION POINT — 2026-10-06 (third), P1 step (a) merged, step (b) next

**State:** `main` is `c5a629bc`, pushed, CI `gates` run green on it (the run before, on `cafdb163`, was red on my own SCOPE-BOUNDARY tag spanning more than three lines; fixed in `c5a629bc`). It holds P0, P0.5, the ruled P1 API sketch (`847e8bd9`) and P1 step (a): `Chord.transits` and `Chart.directions_at` (`2b2d7703`), their theory (`cafdb163`). No branch is open apart from the parked hoist `refactor/chord-oracle-axial-lift` (`a336bde4`, the speed target). The tree is clean apart from `scratch/`. Open issues from this stretch: #580, #581, #582.

**Read in order:**
1. the ledger entry "2026-10-06, the user, on P1's API sketch" and the section "P1 API sketch" (items 1-9: transits, walls, the boundary resolvent in two parts, the panel basis, the Galerkin assembly over lines, the questions on `DensePencil` with the flux coefficients as the unknown, the reading, the reference, the migration order);
2. "P1 step (a) landed" (what exists; the grazing-band fact the reading must respect);
3. the verification spec `scratch/characteristic_architecture/p1_verification_spec.md`, rows A and B (walls, transits, the closure) for step (b)'s first rung; the albedo-pairing rule in its §0;
4. the prototype `scratch/characteristic_architecture/p1_assembly/proto.py` (`chord_psi`, `chord_galerkin`, the slab's closure) and `cyl.py`: the working Galerkin-over-lines assembly the package re-spells on the kernel;
5. the kernels: `orpheus/geometry/{chart,chord,line}.py`, `orpheus/derivations/common/dense_pencil.py`, the boundary-law factors `orpheus/geometry/boundary/_factors.py` (`SelfPairedDeck`, `PairedDeck`, `ScalarResponse`, `SpecularReemission`, `LambertianReemission`).

**The next step, P1 step (b), first rung (an API checkpoint with the user before code; surgical carve, the main agent writes it):** the package `orpheus/derivations/continuous/characteristic/`, inside the insulation gate. First rung: `walls.py` (a `Wall` per boundary point, keyed by breakpoint index, read from the law's own factors: specular amplitude, diffuse amplitude, deck mirror or wrap; `PrescribedInflow` refused as a SCOPE-BOUNDARY; a `LawSum` of specular and diffuse decided at build) and the line part of `closure.py` (`LinePeriod`: the backward period's directed transits from `Chord.transits`, each with the wall the BACKWARD path reflects at; the least solution, exactly 0 on a trapped line with no source; the periodic deck continues without reversal). Gates: the spec's A1-A6 and B1-B10 (the test-architect re-specifies them onto the built names first, as in step (a)). Then rungs: `WallCoupling` (white), `basis.py`, `assembly.py`, the questions, `reading.py`, `reference.py`, corroboration, re-pointing, retirement.

**Lessons from this stretch:**
- A ledger gate's format is a contract: run the CI gate set (`tests/gates/test_elegance_debt_is_tagged.py` and its siblings, the list in `.github/workflows`) locally before pushing any SCOPE-BOUNDARY or ELEGANCE-DEBT tag.
- A shape decided by string labels at several sites fails silently where one site lacks a case (`tangencies` returned `[]` on an unknown name while `direction` raised): an Enum with an exhaustive match at each site.
- A second computation of a quantity the kernel already computes (the impact parameter) disagrees in the last ulp on a third of samples: delegate to the kernel.
- `git commit` in the same Bash command as `git switch -c` is refused by the hook, which reads the whole command before the branch exists: create the branch in its own command.
- A product `(r - l)(r + l)` under a square root under- or overflows at extreme scales: form it from the ratio (`#582` holds the kernel's own instance).

## P1 step (b), first rung: API sketch (the main agent, 2026-10-06; checkpoint with the user before code)

Every item is a proposal until it is ruled. The rung covers `walls.py` and the line part of `closure.py`, in the new package `orpheus/derivations/continuous/characteristic/`. It answers the sketch's items 2 and 3 (the line half) and the one open point there: a `LawSum` on one wall.

**Facts measured by reading (2026-10-06, at `3de4ea44`).**
- **A `LawSum` cannot be declared on a geometry.** `parse_boundary_law` (`structured_geometry.py:547`) admits only a `BC` tag or a `BoundaryTraceLaw`. `LawSum` and `LawScaled` are `ContentIdentity` nodes, not laws (`_composition.py:143, 205`). So every wall has exactly one return shape today, and the sketch's "decided at build" question dissolves.
- **What a tag means is spelled at three sites:**
  - `BC.to_alpha` (`_tag.py:97`);
  - `specular_albedo` (`reference_body.py`);
  - `_law_from_tag` (`transport/method.py:268`). Only this one builds typed laws, using the face's axis and outward sign. It sits at L2, which the reference may not import.
  The `"partial"` kind (partial specular) is known only to the first two.
- **Every law reports a `source`.** It is `NoSource` by default (`_base.py:315`). Only `PrescribedInflow` overrides it.

**1. `Wall` and `Walls` (`walls.py`).**
- `Wall(breakpoint, specular, diffuse, partner)` is frozen:
  - `breakpoint` is 0 or n;
  - `specular` is the amplitude returned along the reflected (or, for a wrap, the translated) line;
  - `diffuse` is the amplitude returned isotropically;
  - `partner` is the breakpoint at which the returned path re-enters. It equals `breakpoint` for a mirror and for diffuse re-emission, and is the opposite wall for a wrap.
- Two amplitude fields are kept, rather than a sum type, because a mixed wall is physically legitimate and the resolvent serves it unchanged (the line part reads `specular`, `WallCoupling` reads `diffuse`). It is undeclarable today, not illegal.
- `Walls.of(geometry)` zips `geometry.boundaries` with the breakpoint indices of `geometry.boundary_points`. It reads each typed law through its two factors only:

| response kernel | under a `SelfPairedDeck` mirror | under a `SelfPairedDeck` identity | under a `PairedDeck` wrap |
|---|---|---|---|
| `ScalarResponse(a)` | specular a | a = 0 is vacuum; a > 0 is refused (the shape is unstated, as SN refuses it) | specular a, partner the opposite wall |
| `SpecularReemission(a)` | specular a | specular a (its own motion is the mirror) | refused |
| `LambertianReemission(a)` | diffuse a | diffuse a | refused |

- **Refusals, each a SCOPE-BOUNDARY:**
  - a law whose `source` is not `NoSource`;
  - a wrap whose partner does not wrap back (the deck is an involution);
  - a wrap on a radial chart (a translation is not in the sphere's or the cylinder's group).
- **Arrays.** `Walls` exposes `specular`, `diffuse` and `partner` as arrays indexed by breakpoint. Interior breakpoints, and a solid body's r_0 = 0, are absent. Indexing a non-wall raises, as `Transits`' absent code n + 1 already does.

**2. `LinePeriod` (`closure.py`, the line part).**
- **Built by `LinePeriod.of(chord, walls)`, batched like the chord.**
- **Its traversals.** The period is the ordered sequence of directed traversals of the unfolded line. A traversal is one of the line's transits read forward or reversed. The traversal after one that exits at wall w is the one entering at `walls.partner[w]`, with a forward transit preferred over a reversed one; the two candidates that both enter there trace the same orbit-space path. The sequence closes when it returns to its first traversal.
- **Its rank is derived, never tagged.** The rule gives m = 1 on a solid body, on a shell's ray that misses the cavity and on the periodic slab; m = 2 on a shell's ray through the cavity and on the mirrored slab; m = 0 on a parallel line. There is no branch on the chart or on b against r_0.
- **Fields**, each of shape (..., 2): `transit` (the index into `chord.transits`), `reversed`, `amplitude` (the exit wall's specular amplitude), `present`.
- **`optical_depth(sigma_t)`**: the optical depth of each traversal, given Sigma_t per region, as the sum of the slot lengths times Sigma_t over its slots.
- **`inflow(optical_depth, outflow)`**: the least solution of the cycle in_{k+1} = a_k (e^{-tau_k} in_k + B_k), returning the inflow at each traversal's entry.
  - B_k is the source integral attenuated to the traversal's exit. It is computed by the assembly, with trailing axes for the basis.
  - The input is tau, not the transmission e^{-tau}, so that 1 - P (P the cycle product) is formed as `-expm1(sum(log a_k) - sum tau_k)` without cancellation on a nearly lossless line.
  - Where P = 1 (every amplitude 1, every tau 0: a lossless trapped line), the answer is exactly 0 when the outflow is zero. Otherwise it is refused, which is the only refusal.
  - With amplitude 0 at every wall the result is exactly 0, so the closure adds nothing to the vacuum operator (B4a, bitwise).
- **The albedo pairing (spec §0) is the cycle's indexing**: the inflow to traversal k + 1 carries the amplitude of the wall where traversal k exits, which is the wall where the backward path from k + 1 reflects.

**3. Gates at this rung.** The test-architect re-specifies onto these names the rows that need no basis:
- A1 and A6: the rank, read as the period's length, and the walls;
- A3: the reflection invariants;
- B1: `inflow` against the unfolded wall-by-wall sum on abstract data;
- B3: the pairing zero, at the closure level;
- the least-solution legs of B5b at the closure level;
- B4a's closure leg.
A2's tau half lands here; its B half waits for the assembly. The rows that read psi on real geometry (A2's B half, A4, B2, B4b/c, B5a, B5c, B6, B8, B9, B10) land with the assembly and the reading.

**Open for the user:**
- **Q1, the tag arm.** The reference must read `BC` tags as well as typed laws.
  - (a) Recommended. Move the parse body of `_law_from_tag` (everything after the method's admission check) down to `orpheus/geometry/boundary` as one function from (tag, axis, outward sign) to a typed law, with `"partial"` added as `AlbedoBoundary(a, SpecularReturn(axis))`. `transport.method` calls it, and `walls.py` reads factors only. One parse of a tag's meaning instead of three. `to_alpha` and `specular_albedo` retire with the old family.
  - (b) `walls.py` reads tag kinds itself, as `specular_albedo` does. That is a fourth spelling of the tags' meaning.
  - Sizing of (a): one production commit before the rung (the move, its gates, transport delegating).
- **Q2, `PrescribedInflow`.**
  - (a) Recommended. Read through its factors like every law, so a `PrescribedInflow` with `NoSource` is a vacuum wall, and only a non-zero source is refused.
  - (b) Refuse the class outright, as the sketch wrote.

**Ruled 2026-10-06** (see the ledger): Q1 the reference's own tag registry, tag to typed law, walls read from factors only (not the hoist); Q2 by its factors; Q3 the period and closure as sketched. Next: the test-architect re-specifies the rung's gates onto these names; then the code on `feature/characteristic-walls-closure`.


## P1 step (b), first rung landed (2026-10-06)

**Commits:** `72199f9d` (walls, line closure, 175 gate rows) and `d3db5932` (theory page `docs/theory/references/characteristic.rst`; labels `characteristic-transit-rank` and `characteristic-closure`; 127 rows moved to l0, 48 stay foundation). Read merge status from git, not from this line.

**What it is, as built:**
- **`Walls`** has the fields `walls`, `n_regions` and `chart`, and checks its own invariants: distinct walls at {0, n}, and a wrap only on the slab and only paired. `Wall` checks that its amplitudes lie in [0, 1].
- **The factor table.** A wall is read by these pairs and no others:
  - the identity deck with `ScalarResponse(0)` is vacuum;
  - a mirror or a wrap with `ScalarResponse(1)`;
  - the identity deck with `SpecularReemission(a)` or `LambertianReemission(a)`;
  - every other pair is refused, an attenuated quotient deck included.
- **`TAG_REGISTRY`** is public. Each tag kind takes exactly its declared parameters, so `BC("white", {"albedo": a})` is refused; partial white is spelled `WhiteBoundary(albedo=a)`.
- **`LinePeriod`** has the fields `chord`, `candidate` (0..3), `amplitude` and `present`. `transit`, `reversed`, `rank`, `entry_wall` and `exit_wall` are derived properties.
- **`optical_depth`** refuses a `sigma_t` that is not of shape (n,). It multiplies only on crossed slots with positive cross section, so no inf·0 is ever formed.
- **`inflow`** is one `np.roll` expression for every rank, with 1 - Pi formed by `expm1`; `TrappedSource` is raised for a source on a lossless trapped line.

**Review decisions taken by the main agent (not user rulings; the user may overrule):**
- the elegance review's C1: a quotient deck is admitted only with amplitude 1, because production retired the attenuated mirror on 2026-10-01;
- qa's tag-parameter finding: the reference refuses undeclared parameters. Production silently drops them, filed as #583.

**Evidence** (`[M]` 2026-10-06):
- the rung's gates, the geometry gates, the layer and ledger gates: 2119 passed, 5 skipped;
- CI's remaining gate set: 97 passed;
- the 34-arm battery reddened on every non-null arm (`scratch/characteristic_architecture/p1_step_b1/battery/summary.txt`);
- Sphinx `-E -W` clean, `dead_references` 0 of 66, pyright 0 errors.

**Filed:** #583 (production's tag parse drops undeclared parameters). The kernel's axial overflow warning on subnormal Omega_x is a comment on #582.

**Owed:**
- the 12 optical-depth rows stay `foundation` until a label for the traversal's optical depth exists;
- spec rows B2 (psi against explicitly reflected lines), A2's B half, A4, B4b/c, B5a, B5c, B6, B8, B9 and B10 land with the source integrals, the assembly and the reading.

**Next:** the second rung, the source integrals B_k per traversal on the panel basis (`basis.py`) and `WallCoupling` (the white block). This is an API checkpoint with the user before code.

## ⏸ COMPACTION POINT — 2026-10-06 (fourth), P1 step (b) rung 1 merged, rung 2 next

**State.** `main` is `884ee4a8` and is pushed. CI `gates` is green on it. It holds P0, P0.5, P1 step (a) and P1 step (b) rung 1:
- `72199f9d`: `walls.py`, `closure.py` and 175 gate rows;
- `d3db5932`: the theory page `docs/theory/references/characteristic.rst`;
- `884ee4a8`: the plan and agent memory.
No branch is open apart from the parked hoist `a336bde4` (the speed target) and this compaction branch. The tree is clean apart from `scratch/`. Open issues from this stretch: #580, #581, #582 (with a comment on the axial overflow warning), #583.

**Read in order:**
1. The section "P1 API sketch", items 3-6. Item 3 is the boundary resolvent's two parts, item 4 the panel basis, item 5 the Galerkin assembly over lines, item 6 the questions on `DensePencil`.
2. "P1 step (b), first rung: API sketch" and "P1 step (b), first rung landed": what exists and its names.
3. The theory page `docs/theory/references/characteristic.rst`, sections "What this rung does not compute" and "Gotchas".
4. The verification spec `scratch/characteristic_architecture/p1_verification_spec.md`. Read §0 and rows C, plus B2, B8, B9 and A2's B half.
5. The prototype `scratch/characteristic_architecture/p1_assembly/proto.py`: `chord_psi` (the per-piece emission integral `E`, the attenuation chain `I`), `chord_galerkin`, `slab_mu`/`_finish` (the intra-piece triangle), and `exp_graded`/`graded_rule`. Also `cyl.py`.
6. The reference kernel: `orpheus/derivations/common/{quadrature,dense_pencil}.py`.

**The next step: P1 step (b), rung 2.** It is an API checkpoint with the user before any code, a surgical carve written by the main agent.
- **`basis.py`:** `PanelBasis(geometry, degree, grading)`.
  - Per-region polynomial panels, with nodal Lagrange functions at each panel's Gauss-Legendre points.
  - Panels graded toward walls and interfaces, with no node on an interface.
  - It records each node's region and panel, and gives the mass matrix W in the chart's measure.
- **The source integrals B_k per traversal**, on the basis. They feed `LinePeriod.inflow`.
- **`WallCoupling`:** the white W x W block per group, treating every diffuse wall as absorbing in the line part.
  - U is the escape functional; T_w is the wall-to-wall transmission.
  - K_g = K_line + U (I - T_w)^{-1} A U^T, as ruled.
- **Then** `assembly.py`: the Galerkin double integral over lines, then the questions, the reading, the reference, corroboration, re-pointing and retirement.
- **Open for the API sketch:** where B_k is computed (on `LinePeriod`, or in the assembly) and its arc-length quadrature (graded exponential panels on optically thick slots, as in the prototype). Also: the label for the traversal's optical depth (12 rows wait on it).

**Lessons from this stretch:**
- **Masking a product after forming it still raises.** `np.where(mask, a*b, 0)` evaluates inf·0 and warns, and under `-W error` it raises. Multiply with `where=` and `out=` so the product is never formed.
- **Pad nothing that the kernel made out of range.** Padding a per-region table to cover the exterior codes made a wrong-length table silently valid (qa). Refuse the shape and mask the codes instead.
- **Refuse what the parse does not declare.** A tag parameter the parse does not read was dropped silently in production (#583).
- **Invariants belong in `__post_init__`, not only in the factory.** A directly built `Walls` spelled the states `Walls.of` refused (elegance S1).
- **`np.flip` and `np.roll` agree on an axis of length 2.** A mutation swapping one for the other there is equivalent, so it is not a battery arm.
- **A test can pick its expectation from the code under test.** One row chose its "absorbing traversal" from the code's own amplitudes and stayed green under the pairing mutation. Key expectations to independent data, such as the exit wall.

## P1 step (b), second rung: API sketch (the main agent, 2026-10-06; checkpoint with the user before code)

Every item is a proposal until it is ruled. It covers the panel basis and the transport along one line on that basis. It also proposes re-cutting the ladder.

**Facts measured by reading (2026-10-06, at `5f1e9d25`).**
- **The white block needs the line measure.** Its escape functional U (emission to the outgoing partial current at each diffuse wall) and its wall-to-wall transmission T_w are both integrals over the measure on lines (`Chart.beam_density`). The Galerkin assembly needs that same measure, and so does nothing else in this rung.
- **A panel end is a level set like a breakpoint.** `ConcentricPartition(chart, breakpoints)` admits any strictly increasing breakpoints. The kernel's chord through a partition whose breakpoints are the basis's panel ends therefore yields every piece of every line: one slot per panel crossed, its length cancellation-free, its panel read from `slot_region`, and the closest approach already a slot end. The prototype computes these pieces with its own crossing arithmetic (`sphere_pieces`, `slab_pieces`).
- **The volume density is not exposed.** The kernel has `Chart.measure(edges)` (cell measures, `measure_constant · (T(r_{j+1}) − T(r_j))` with T(r) = r^d) but no density in c.
- **The prototype's per-line Galerkin block has three parts.** For one traversal, with the inflow at its entry written in:
  - V, the Volterra triangle: the double integral of u_i(s) e^{−τ(s′→s)} u_j(s′) over s′ < s;
  - A, the attenuated entry test: the integral of u_i(s) e^{−τ(entry→s)};
  - B, the outflow: the integral of u_j(s′) e^{−τ(s′→exit)}.
  The line's block is V + A ⊗ in, where in is `LinePeriod.inflow` of the B's. A of a traversal equals B of the reversed traversal, and V of the reversed traversal equals Vᵀ (reciprocity along one line).

**0. Re-cutting the ladder (Q1).** The ruled ladder (P1 sketch item 9b) is "basis and assembly (C)". My compaction note put `WallCoupling` in this rung; that was my grouping, not a ruling. Proposed:
- **Rung 2 (this sketch):** the basis and the transport along one line: B, A, V, and ψ along a line. Spec rows A2 (the B half), B5c, C2c, C3, C10 and C12, plus the reciprocity identities above.
- **Rung 3:** the line rule (the measure on lines), the assembly K_line, and `WallCoupling` with it, since all three are integrals over that measure. Spec rows C6 to C9, C11, B4b, B5a and B8's operator half.

**1. `PanelBasis` (`basis.py`).**
- **Panels.** Each region [r_k, r_{k+1}] is cut into panels graded geometrically toward each end that is a wall or an interface: interior panel ends at the depths w·ratio^j, j = 1..`layers`, from each graded end, with w the half width of the region when both its ends are graded and the whole width when one is, the middle left whole (so the two panels at the end have equal width at ratio 1/2: the law places depths, not widths; corrected to the coded law 2026-10-06). An end that is a singular stratum (a solid body's centre or axis) is not graded, because the flux is even and smooth there.
- **Functions.** Nodal Lagrange functions of degree p through each panel's p + 1 Gauss-Legendre points, from `composite_gauss_legendre` (the P0.5 kernel). They are discontinuous across panel ends, and no node lies on a panel end. Because the basis is nodal, `DensePencil.fundamental`'s single-sign test reads coefficients as values.
- **Fields and verbs:**
  - `partition`: the panel ends as a `ConcentricPartition`, a refinement of the geometry's;
  - `region` and `panel` of each node, and `region_of_panel`;
  - `values(c, panel) -> (..., p+1)`: the panel's functions at orbit coordinates;
  - `mass`: W = ∫ u_i u_j dV, block-diagonal by panel, by Gauss with p + 2 points on each panel. The density is `measure_constant · d · c^{d−1}`, derived from the one definition of the measure. A gate checks that the density integrates to `Chart.measure` of each panel.
- **Parameters** (degree, layers, ratio) live in the `Resolution` value. The basis is built from the geometry and that value alone, never from reading points (C6b).

**2. The pieces come from the kernel (Q2).** Every line is chorded through `basis.partition`, not through the geometry's partition. `LinePeriod.of(panel_chord, walls)` then runs on that chord: the walls are the same two ends, re-keyed from breakpoint n to the panel count P. The proposed spelling is `Walls.on(partition)`, which re-checks the invariants. Σ_t per panel is `sigma_t[basis.region_of_panel]`. One chord serves the closure, τ and the source integrals, so no second crossing computation exists.

**3. `TraversalRule` (`transport.py`, new) (Q3).** It is the quadrature along each traversal of a `LinePeriod`, built by `TraversalRule.of(period, basis, sigma_per_panel, resolution)`.
- **Nodes.** Each traversed slot gets Gauss nodes in arc length:
  - a slot ending at the closest approach is mapped by c = b + (c_end − b) u², so its integrand is polynomial in u (the prototype's turning map);
  - a slot of optical length above 2 is cut into exponential panels toward its end (the prototype's `exp_graded`: widths 1/Σ, doubling up to 50 mean free paths).
- **Storage.** Arrays (..., 2, Q), one row per traversal of the period. Q is a fixed maximum over the batch, and an unused node carries weight 0 from a zero-length piece; no out-of-range code is padded. Each node carries its panel and the optical depth from the traversal's entry.
- **Verbs**, each linear in the basis:
  - `outflow() -> (..., 2, N)`: B_k, the source integral attenuated to the traversal's exit. It feeds `LinePeriod.inflow`.
  - `entry_response() -> (..., 2, N)`: A_k, the test u_i attenuated from the entry. Gate: it equals the reversed traversal's `outflow`.
  - `volterra(line_weight) -> (N, N)`: Σ_L w_L Σ_k V_k, the triangle accumulated over the batch with the caller's line weights, never stored per line. The inner integral uses its own Gauss sub-rule from the piece start to each outer node, as the prototype's `chord_psi` does.
  - `angular_flux(t, inflow) -> (..., N)`: ψ at parameters on the line, for the reading (rung 5). It is listed so that the reading reuses the same rule, which C6 requires.
- **Void slots** (Σ = 0) integrate as q·ℓ with no division (B5c); the exponential cut applies only where Σ > 0, using `where=` multiplication as in rung 1.

**4. The inflow argument for the white block (rung 3, named now).** `LinePeriod.inflow(optical_depth, outflow, arriving=0)` will gain the inflow arriving at each traversal's entry from outside the line's own cycle: a diffuse wall's re-emission now, a boundary source later. With an absorbing specular amplitude the outflow term vanishes, but the arriving term does not, so T_w for a pure white wall is spellable. It lands with `WallCoupling`.

**5. The owed label.** The theory page gains the traversal integrals τ_k and B_k as one equation, labelled `characteristic-traversal-integrals`. The 12 optical-depth rows move from `foundation` to `l0` under it, with A2's B half.

**Open for the user:**
- **Q1, the ladder.** (a) Recommended: re-cut as in item 0, with `WallCoupling` in rung 3 beside the line rule. (b) Keep `WallCoupling` in rung 2 and build the line rule with it, so that rung 2 is the whole of C.
- **Q2, where the pieces come from.** (a) Recommended: the panel ends as a refined partition, chorded by the kernel, with `Walls.on(partition)`. (b) Chord the geometry's partition and split each slot at the panel-end crossings inside the reference, which is a second crossing computation beside the kernel's.
- **Q3, where B_k lives.** (a) Recommended: the new `TraversalRule`, because B needs the basis, Σ_t and a quadrature, which `LinePeriod` (geometry and walls) does not hold. (b) As methods on `LinePeriod`.
- **Q4, the basis.** (a) Recommended: discontinuous nodal panels graded toward walls and interfaces and not toward a singular stratum, with the density derived inside the basis from `measure_constant` and `measure_coordinate`. (b) The same, with a kernel verb `Chart.volume_density(c)`, which changes the stable kernel and would need its own ruling.

**Ruled 2026-10-06** (see the ledger): Q1 to Q4 as recommended. Rung 2 is the panel basis and the transport along one line. Rung 3 is the line rule, the assembly and `WallCoupling`. The pieces come from the kernel's chord through the panel partition, with `Walls.on(partition)`. B_k lives on the new `TraversalRule`. The density is derived inside the basis. Next: the test-architect re-specifies the rung's gates onto these names, and the main agent writes the code on `feature/characteristic-basis-transport`.

## P1 step (b), second rung landed (2026-10-06)

**Commits:** `ff979520` (basis, line transport, 328 gate rows) and `fc8ddbd4` (theory sections, label `characteristic-traversal-integrals`, ERR-099 and ERR-100). Read the merge status from git, not from this line.

**What it is, as built:**
- **`PanelBasis(regions, partition, degree)`, built by `of(regions, degree, layers, ratio)`.**
  - Panel ends are at w·ρ^j from each graded end, where w is the half width when both ends are graded.
  - A singular stratum is not graded.
  - `region_of_panel` is derived, so "a panel in the wrong region" cannot be spelled.
  - The mass density is `measure_constant · d · c^(d−1)`.
- **`TraversalRule.of(lines, basis, walls, sigma_t, points, inner_points)`** (the user's ruling after the elegance review's C2). It chords the lines through `basis.partition` itself, and `__post_init__` requires that identity. Its verbs are `outflow`, `entry_response`, `volterra(line_weight)` (over each line's transits read forward) and `angular_flux(t, inflow)`.
- **The internals:**
  - each line keeps only its live pieces, in chord order;
  - one graded body, `_attenuated(slot, start, stop)`, computes every carried integral;
  - B and A are the exit and entry transmissions times the piece integrals;
  - the carried integral is a scan.
- **`LinePeriod.forward_traversal`; `Walls.on(partition)`.**

**Decisions taken by the main agent (not user rulings; the user may overrule):**
- The u-map on turning slots (c = b + (c_far − b)u²) is retired. With grading in place it changed no gated reading (the test-architect's F3). It also leaves a branch point in the Jacobian.
- The grading toward branch points halves rather than quarters. Quartering drifted to 1.5e-14 at 12 points.
- Grading covers every radial slot's end nearer the closest approach, down to c_near/|PΩ| (qa's finding: ERR-099).
- Exit-wall reads are located by the crossing parameter.

**Defects found in review and fixed in the session:**
- ERR-099, the branch-point grading gap: 8.5e-7 on a hollow sphere with cavity radius 0.01.
- ERR-100, the ungraded carry out of a thick slot's middle piece: ψ off by 0.54 at 1000 mfp, by 5.2e-6 at 200 mfp.
- Exit-wall reads refused by rounding: 64 of 119.

**Evidence** (`[M]` 2026-10-06):
- the four characteristic gate files: 503 of 503 under `-O` (189 l0, 314 foundation);
- battery: 30 of 30 arms red (`scratch/characteristic_architecture/p1_step_b2/battery/summary.txt`);
- the ledger, layer and geometry gates: 1946 passed;
- the CI gate set: 1420 passed;
- Sphinx `-E -W` clean; `dead_references` 0 of 66 (a planted control read 1 of 67); pyright 0 errors.

**Deferred, from review (not done):**
- elegance N2: three geometric gradings (`_graded_ends`, the exponential `_toward`, `_branch_edges`), to collapse when rung 5's grazing grading arrives;
- N5: a batched Gauss rule over arrays of intervals is missing from the reference kernel;
- the reframing of B, A, V and ψ as one semiseparable Volterra operator per line (the sweep its matvec, a reversed traversal its transpose);
- `MeasureCoordinate.derivative`, to make when the kernel next opens.

**Next:** rung 3, as an API checkpoint with the user before code:
- the line rule (the measure on lines through `Chart.beam_density`, with the cosine endpoint map);
- the assembly K_line = Σ_L w_L (V_L + A_L ⊗ in_L) over forward traversals;
- `WallCoupling` (U, T_w, and the `arriving` argument of `LinePeriod.inflow`).

## P1 step (b), third rung: the premises re-measured (the main agent, 2026-10-06; before the API sketch)

Rung 3 is the line rule (the measure on lines), the assembly K_line, and `WallCoupling`. Before sketching it, the premises of P1 sketch items 3 and 5, and spec rows C6-C9, C11, B4b, B5a and B8, were measured against the code that now exists.

Every probe is under `scratch/characteristic_architecture/p1_step_b3/`. Each assembles K = Σ_L w_L (V_L + Σ_{forward k} A_k ⊗ in_k) from `TraversalRule`. All runs use `.venv/bin/python -O` at `e96cde98`.

**1. The normalisation holds as derived.** The line weight is (beam density × direction weight) / 4π, folded over the chart's symmetry:
- sphere: one direction, weight 2πb db;
- cylinder: μ_z ∈ [0, 1] doubled, weight |PΩ| db dμ_z;
- slab: μ ∈ [−1, 1], weight |μ|/2 dμ.

With these weights K is the scalar-flux operator for an isotropic emission density. On a closed homogeneous body, K·1 = W·1/Σ_t holds to:
- sphere: 2.5e-13 (`ladder_probe.py`);
- slab: 4.9e-15 (`slab_mu.py`);
- cylinder: 2.6e-10 at 16 impact-parameter points (`cyl.py`).

**2. C11, the symmetry alarm, is designed-green (a stabiliser finding, X1).** K is symmetric to at most 2.5e-16 at every resolution tried, under-integrated ones included:
- the sphere at 8 b-points, where K·1 misses by 5.2e-5;
- the slab with a plain μ rule, which misses by 1.6e-3;
- the cylinder at every rule.

The reason is that each line's quadrature is symmetric under reversing the line. The inbound and outbound halves of a radial chord mirror each other, and the slab's μ rule is symmetric. So every line's block is symmetric whatever the accuracy. The prototype's 1.4e-4 defect came from its own mismatch between inner and outer rules, which the new code does not have.

The closed-body identity K·1 = W·1/Σ_t does redden under under-integration, with the misses measured above. It is the alarm C11 claimed to be.

**3. The singularities in the line measure.**
- **Sphere and cylinder, impact parameter b.**
  - The square-root ends at each radius are absorbed by `chord_quadrature`, the visibility-cone substitution applied per panel. Away from the centre, 16 points give 2e-11 on each panel.
  - The panel touching the centre converges only algebraically (`ladder3.py`, per b-panel against 128 points): 3.5e-8, 1.6e-10, 7.0e-13 at 8, 16, 32 points.
  - The cause is that the basis is polynomial in c, not c², on the panel touching the stratum. The line integral of an odd power of c is an Abel transform, ∫_b c^{2m+1} c dc / √(c² − b²), and it carries a b^{2m+2} log b term.
  - The true flux is a smooth invariant function. By Schwarz's theorem (a smooth O(d)-invariant function is a smooth function of |x|²), it is a smooth function of c² near the stratum, so the odd modes are what the basis adds, not the physics.
  - Grading the b rule geometrically toward 0 takes the closed sphere from 5.3e-8 to 4.4e-11 at 16 points (the remainder is the outer panel). At 8 points both stall near 5e-5, which is the outer panel's convergence and not the centre.
- **Slab, direction cosine μ.** Plain Gauss on [0, 1] reaches 1.6e-3, 1.2e-5, 6.9e-9 at 8, 16, 32 points. Grading geometrically toward μ = 0 (12 halvings) reaches 5.8e-14 at 8 points, with or without reflection (`slab_mu.py`). The prototype's `mu_rule` grading is confirmed.
- **Cylinder, axial cosine μ_z.** No grading is needed. Plain Gauss at 8 points gives the same 2.6e-10 as 16 points or 12 graded layers; the floor is the b rule. The graded version costs 10× more and, run in one batch, exhausted memory (exit 137).

**4. Cost and memory.** One sphere block at N = 32 with 288 lines takes 0.45 s. One axial cosine of the cylinder with about 300 b-lines takes about 0.4 s. A single batch of about 30,000 lines exhausted memory, so the assembly must take its lines in chunks.

The pieces depend on Σ_t, so each group needs its own `TraversalRule`. The geometry (the chord and the period) can be shared across groups through the dataclass constructor, which checks that the partitions are identical.

**5. The white block's formula, corrected.** The measure on lines through a surface element is |Ω·n| dA dΩ, so the escape functional U[i, w] is Σ_L w_L Σ over forward traversals exiting at w of (e^{−τ_k} in_k + B_k)_i, with the same line weights as the assembly. An isotropic inflow of total current J at wall w has ψ = J/(π A_w). Reciprocity then makes its volume response (4/A_w)·U[:, w].

The block is therefore K = K_line + U D^{−1} α (I − T_w α)^{−1} Uᵀ, where:
- D = diag(A_w/4);
- α holds the diffuse amplitudes;
- T_w[w′, w] is the current out at w′ per unit isotropic current in at w, through the line part. It needs the `arriving` argument of `LinePeriod.inflow`. Reciprocity makes D^{−1}T_w symmetric.

The ruled spelling "U (I − T_w)^{−1} A Uᵀ" omits D. A_w is the area of the level set c = r_w, which is the volume density at r_w (4πR², 2πR per unit height, 1 per unit area): `PanelBasis.volume_density` already computes it, so a second consumer of the density now exists.

**6. The cross-section blocks.** `_infinite_medium_matrices` (`derivations/common/eigenvalue.py:41`) returns the 0-D pair (A = diag Σ_t − (Σ_s + 2Σ_2)ᵀ, F). The pencil needs the emission matrix (Σ_s + 2Σ_2)ᵀ on its own, because K already carries Σ_t. Reading it from A would mean subtracting diag Σ_t back out. These blocks are consumed only by the pencil.

**Ruled 2026-10-06** (see the ledger), all four as recommended:
- **Q1, the centre.** The basis is even at a singular stratum. On the panel touching the centre or the axis, the nodal Lagrange functions are polynomials in c² through the same Gauss points. This is a change to the merged `PanelBasis`, and its mass rule gains points (the degree in c doubles).
- **Q2, the line domain.** It is a kernel verb on `Chart`: the orbit space of lines under the chart's group, with its density, gated against `beam_density`. It is the lines' analogue of `directions_at`. The reference builds its quadrature over it:
  - the visibility cone at each panel end in b;
  - geometric grading toward μ = 0 on the slab;
  - plain Gauss in μ_z on the cylinder.
- **Q3, the alarm.** C11 is re-posed as closed-body conservation, K·1 = W·1/Σ_t, on closed homogeneous bodies (specular and white) on every chart. Symmetry stays only as a foundation row, declared blind.
- **Q4, the scope.**
  - The volume density moves to the kernel: a derivative of the measure (`MeasureCoordinate.derivative`, or a `Chart` verb). The basis and the wall area both read it.
  - The white block is K = K_line + U D^{−1} α (I − T_w α)^{−1} Uᵀ with D = diag(A_w/4).
  - S and F, an emission-matrix helper split out of `_infinite_medium_matrices`, move to rung 4 with the pencil.

## ⏸ COMPACTION POINT — 2026-10-06 (fifth), rung 2 merged, rung 3 ruled, API sketch next

**State.** `main` is `e96cde98` plus this plan commit. CI `gates` is green on `e96cde98`. Rungs 1 and 2 of P1 step (b) are merged; rung 3's premises are re-measured and its four questions ruled. Open issues from this stretch: #580 to #583.

**Read in order:**
1. "P1 step (b), second rung landed": what exists and its names.
2. "P1 step (b), third rung: the premises re-measured", with its Ruled line.
3. The P1 API sketch, items 3 and 5: the boundary resolvent's two parts, and the assembly.
4. The theory page `docs/theory/references/characteristic.rst`: the basis and transport sections, and "What the package does not compute".
5. The code: `orpheus/derivations/continuous/characteristic/{basis,transport,closure,walls}.py`, `orpheus/geometry/{chart,coord}.py`.
6. The probes `scratch/characteristic_architecture/p1_step_b3/*.py`. They are the measured templates for the assembly loop, the line weights and the conservation check.

**The next step: rung 3's API sketch,** a checkpoint with the user before code. It must spell:
- **the kernel verbs, with gates:**
  - the measure density, as `MeasureCoordinate.derivative` or a `Chart` verb;
  - `Chart.line_domain()`, with its fundamental domain and density per chart;
- **the even basis at a singular stratum** in `PanelBasis`;
- **`LinePeriod.inflow(..., arriving=...)`;**
- **`LineRule`** in `assembly.py`: a quadrature over the line domain, taken in chunks of lines, with one `TraversalRule` per group sharing the chunk's period;
- **the assembly** K_line, per group;
- **`WallCoupling`:** U, T_w, D, and the corrected formula;
- **the gates:**
  - closed-body conservation (the re-posed C11);
  - C6 to C9;
  - B4b, B5a, and B8's operator half;
  - the self-convergence ladders.

The test-architect then re-specifies the rows, and the main agent writes the code.

**Lessons from this stretch:**
- **The battery must sweep resolutions below the working point.** At 16 points the gates could not see a mechanism that changed the error only at 12 (ratio 4 against halving) or at 8.
- **A symmetry row on a mirror-symmetric quadrature cannot fail.** A designed-green stabiliser (X1) went unnoticed through a whole spec until K was assembled.
- **A basis polynomial in c at a singular stratum spans non-smooth functions.** Schwarz's theorem says the smooth invariant ones are functions of c². The odd modes surfaced as a b² log b term in every line integral.
- **One wide batch can exhaust memory.** Chunk the lines.
- **QA found three real defects that the gates passed:**
  - the branch-point gap: no slot started at a small radius;
  - the thick-slot carry: no slot was thicker than 128 mean free paths;
  - exit-wall rounding: no read at the crossing parameter.
  Each was a missing input region, not a weak tolerance.

## P1 step (b), third rung: API sketch (the main agent, 2026-10-06; checkpoint with the user before code)

Every item is a proposal until ruled. It spells what the fifth compaction point listed, on the names that exist at `e1faa622`. Rung 3 delivers the per-group transport block K_g (line part plus diffuse part) and nothing above it: no pencil, no point reading.

**1. The measure density (kernel, `orpheus/geometry/coord.py` and `chart.py`).**
- `MeasureCoordinate.derivative(r)`: T'(r) = p r^(p−1), the derivative of the one definition T(r) = r^p.
- `CoordSystem.measure_density(r)`: c T'(r), beside `CoordSystem.measure(edges)` = c (T(b) − T(a)), which it differentiates; `Chart.measure_density` delegates, as `Chart.measure` does.
- One function, two readings. It is the volume density in the orbit coordinate, and, by the coarea formula with |∇c| = 1 on all three charts, the area of the level set c = r: 4πr² (sphere), 2πr per unit height (cylinder), 1 per unit area (slab). The basis's mass matrix and the white wall's area both read it, and `PanelBasis.volume_density` retires onto it.
- Gates: the closed forms on each chart; Gauss on [a, b] of the density equals `measure([a, b])` to a few ulp (the derivative and the definition cannot drift).

**2. The line domain (kernel, `Chart.line_domain() -> LineDomain`).** The orbit space of oriented lines under the chart's group G_c, the lines' counterpart of `directions_at`.
- A `LineShape` enum with an exhaustive match at each site (the step (a) lesson):
  - sphere, `IMPACT`: b ∈ [0, ∞); the reversed line is in the same orbit;
  - cylinder, `IMPACT_AXIAL`: b ∈ [0, ∞) × μ_z ∈ [0, 1]; the z-reflection and the C2 axes of D∞h fold μ_z and the orientation;
  - slab, `COSINE`: μ = Ω_x ∈ [−1, 1]; the slab's group fixes the kept space, so μ and −μ are two orbits.
- `density(coordinates)`: the invariant measure dA⊥ dΩ on oriented lines in those coordinates, per unit of the discarded measure (per unit height, per unit area): `beam_density(b, Ω)` times the direction measure the quotient folds (4π on the sphere; 2π × 2 on the cylinder; 2π on the slab). It is not normalised by 4π; that belongs to the consumer (item 5).
- `lines(coordinates) -> Line`: representative lines, as `DirectionDomain.direction` gives representative directions; their impact parameter is `Chart.image`'s, bit for bit.
- The domain is unbounded in b. Which lines meet the body is the partition's question, so the consumer truncates at the outer radius.
- Gates: the density equals `beam_density` times the folded direction measure; and an independent row, Cauchy's formula: the integral of the chord length over the domain is 4π times the body's measure (4πV, per unit height or area), on every chart, which ties the density to the chord with no shared formula.

**3. The even basis at a singular stratum (`PanelBasis`).**
- `PanelBasis.even` (cached, ``(P,)``): true for the panel whose lower end is a singular stratum (`chart.singular_strata`), i.e. the centre panel of a solid sphere or cylinder.
- On that panel the nodes are the same Gauss-Legendre points c_m, and the functions are the Lagrange polynomials in s = c² through s_m = c_m²: they span 1, c², ..., c^(2p), every function of c² of degree p and no odd mode. `values` evaluates the product form in the panel's coordinate (c or c²); everything downstream (`TraversalRule`, the mass) reads `values` and needs no change.
- The mass rule takes 2p + 2 Gauss points on every panel (exact for the even panel's degree 4p + d − 1; one rule, no branch).
- Gates: reproduction of even polynomials of degree 2p on that panel to a few ulp; a fit residual on c (the odd mode is absent); the mass against mpmath; the centre b-panel of the closed sphere converges geometrically where it converged algebraically (`[M]` today 3.5e-8, 1.6e-10, 7.0e-13 at 8, 16, 32 points).
- Re-baseline owed: `tests/gates/derivations/_characteristic_mp.py` builds its own Lagrange functions in c; it gains an even Lagrange in c² for that panel, written independently (it does not call the basis).

**4. `LinePeriod.inflow(optical_depth, outflow, arriving=None)`.**
- `arriving`, ``(..., 2, *rest)``, is a flux injected at each traversal's entry from outside the line part (the diffuse re-entry). It is not multiplied by the wall's amplitude.
- One expression still serves every rank: the flux entering traversal k is a_(k−1) B_(k−1) + s_k, and the inflow is (entering_k + g_(k−1) entering_(k−1)) / (1 − Π), with g the traversal's gain a e^(−τ). `arriving=None` is bitwise today's result.
- The `TrappedSource` refusal covers an arriving flux on a lossless trapped line.
- Gates: the rank-1 and rank-2 closed forms with an arriving flux; `arriving=0` bitwise equal to `None`.

**5. `LineRule` (`assembly.py`): the quadrature over the line domain, in chunks.**
- `LineRule.of(basis, walls, points, grazing_layers, chunk)`: coordinates and weights on `chart.line_domain()`, weight = quadrature weight × density / 4π (K is then the scalar-flux operator of an isotropic emission, the normalisation `[M]` in the premises).
  - b (sphere, cylinder): `chord_quadrature` over the panel ends inside [0, r_n], the visibility-cone substitution on each piece (the square-root end at every radius and panel end); the even basis makes the centre piece smooth in b.
  - μ (slab): geometric grading toward μ = 0, `grazing_layers` halvings on each sign (`[M]` 12 halvings: 5.8e-14 at 8 points).
  - μ_z (cylinder): plain Gauss on [0, 1] (`[M]` no grading needed).
- `chunks()`: the lines in batches of at most `chunk`, each with its `LinePeriod` (one chord per chunk through `basis.partition`).
- `transport(sigma_t, points, inner_points) -> GroupTransport` for one group's `sigma_t` (n,): per chunk, one `TraversalRule(period, basis, ...)` through the dataclass constructor, accumulating:
  - K_line = Σ_L w_L (V_L + Σ_(forward k) A_k ⊗ in_k);
  - the escape U[i, w] = Σ_L w_L Σ_(forward traversals exiting at w) (e^(−τ_k) in_k + B_k)_i, at every diffuse wall;
  - the wall transmission T_w[w′, w]: the current out at w′ per unit isotropic current in at w, through the line part (`arriving` = 1/(π A_w) on the traversals entering at w).
- Per group, not all groups in one pass: the algebra is per group (ruled sketch item 5, `transport_block(g)`, renamed `group_transport` at step 2c). The chunk's chord and period are rebuilt per group; the cost is measured at build, and a cache is added only if the geometry is a measured fraction of the time.

**6. `WallCoupling` (`closure.py`, the diffuse part of the boundary resolvent).**
- `WallCoupling(escape U, transmission T_w, area A_w, diffuse α)`: pure linear algebra; `update` = U D^(−1) α (I − T_w α)^(−1) Uᵀ with D = diag(A_w/4), A_w = `chart.measure_density(r_w)`.
- `GroupTransport.block` = K_line + coupling update; with no diffuse wall the update is exactly zero and the block is K_line, bitwise.
- Refused: I − T_w α singular (every wall diffuse with α = 1 and Σ_t = 0 throughout), the white analogue of `TrappedSource`.

**7. The gates for the rung** (the test-architect re-specifies them onto these names, each with its first red):
- **Closed-body conservation, the re-posed C11:** K·1 = W·1/Σ_t on closed homogeneous bodies, every chart, for the specular mirror, the white wall at α = 1, the periodic slab, the hollow sphere with an inner mirror and with an inner white wall. It reddens if D is dropped from the white update, and under every under-integration measured in the premises.
- **Symmetry, a foundation row declared blind:** K = Kᵀ to 1e-14, stating that each line's rule is reversal-symmetric so the row cannot see under-integration.
- **C9:** the escape probability of a homogeneous body from the vacuum block, against Hebert (sphere), 2E_3 (slab), Bickley (cylinder).
- **The wall transmission against closed forms:** T_w = P_ss on the vacuum-interior sphere and cylinder, 2E_3(τ) face to face on the slab; independent of C9 because it goes through the surface lines and not the volume.
- **C7, operator half:** splitting a region into 2, 3 and 5 equal-material regions leaves the region-to-region transfer rᵀ K q unchanged (piecewise-constant q).
- **B4b and B5a, operator level:** a transparent cavity equals an inner mirror, and a void outer layer is invisible to the specular closure, on the blocks restricted to the material nodes.
- **B8, operator half:** the white wall at α = 1 conserves; at α = 0.5 the white and the specular blocks differ by more than 100 times the ladder step; at 2G, a block built from one group's escape for every group reds.
- **The self-convergence ladders,** in b, μ, μ_z and the inner points, sweeping resolutions below the working point (the rung-2 lesson).
- **Chunk invariance:** the block is independent of the chunk size to 1e-15.
- **Moved to rung 5:** C6 and C8 read the flux at a point, which rung 5 builds; at the operator level C8's reciprocity is the symmetry row, which is blind. D11 (the Rayleigh-Ritz bound) needs the pencil: rung 4.

**A wall with both a specular and a diffuse part (a `LawSum` of the two laws) stays refused.** `[M]` 2026-10-06: the SN realizer refuses a `LawSum` tree too (`orpheus/sn/boundary/realizer.py:819`, "a LawSum/LawScaled tree, which composes laws but is not one, falls through to the loud dispatch-failure raise"), so production cannot pose the case and the reference has no consumer for it (verification conforms to production). `Wall` already carries both amplitudes, so serving it later is one arm of the reader plus its gates.

**Ruled 2026-10-06** (see the ledger), both as recommended: `transport` runs per group, rebuilding each chunk's geometry, with a cache only if the geometry is a measured fraction of the time; the `LawSum` refusal stands. The sketch stands as written otherwise.

**Next:** the test-architect re-specifies the rung's gates (item 7) onto these names, each with its first red; then the main agent writes the code in the order kernel verbs (items 1, 2), the even basis (3), `inflow` (4), `LineRule` and `WallCoupling` (5, 6).

## P1 step (b), third rung: the verification spec and what it refuted (2026-10-06)

The test-architect wrote `scratch/characteristic_architecture/p1_step_b3/spec.md`: 35 new rows (13 foundation, 6 l0, 11 l1, 5 l2) and 8 re-posed functions. The code did not exist yet, so it measured every row on a prototype of the sketch's arithmetic (`p1_step_b3/ta/proto.py`); every number is re-measured on the built code.

**Premises refuted (`[M]` the test-architect's probes `ta/m1`-`m11`, and `ta/m12_theta.py`):**
- **The cylinder's μ_z needs no grading: refuted.** The premise was measured on a closed mirror cylinder, where each line conserves on its own, so no direction rule can show there. On a white wall the escape probability C9 misses Bickley's closed form by 5.4e-4 with plain Gauss in μ_z at 8 points (τ = 0.5). The cause is |PΩ| = √(1 − μ_z²), a square-root endpoint at μ_z = 1. Gauss in the polar angle θ = arccos μ_z, with no grading, gives 2.0e-6, 2.3e-9 and 1.7e-12 at 8, 16 and 32 points in 2 to 5 s; grading toward μ_z = 1 reached 7e-11 at 433 s per block.
- **Conservation cannot see the slab's μ rule** under any law: 33 of 33 runs held to 2.5e-15, plain 4-point rules included. The slab's μ rule is gated by C9 and the wall transmission instead.
- **The symmetry row is blind to the line-measure rules (b, μ, μ_z), not to everything.** It sees a mismatch between the outer and the inner arc-length rules: 1.5e-6 at 4 points.
- **Chunk invariance to 1e-15: refuted.** The slab moves by 2.9e-15, 12.9 eps; the row uses 64 eps.
- **Gauss of the density equals `measure` to a few ulp** holds only with the conditioning factor: 2568 ulp raw on the cylinder, at most 1.1 ulp scaled by 1/(1 + κ).
- **B5a at α = 1 cannot be reached on the full block.** The void layer's basis functions are sources on lossless trapped lines, so `inflow` refuses them.
- **A chunk sized by line count does not bound memory.** 224 cylinder lines at μ_z = 0.999 took 8.2 GB.

**Confirmed:**
- dropping D from the white update reds conservation by 6.6;
- diag(A)⁻¹ T_w is symmetric to 1.6 eps;
- a wall carrying both a specular and a diffuse part conserves through the unchanged formula. It stays refused by ruling.

**The even basis's blast radius:** 29 of today's 503 characteristic rows go red under it. Nine are the turning-grading witnesses on solid bodies. They lose their subject: on the even panel the integrand is polynomial in arc length, so no branch point is left to grade toward. One of the nine is the ERR-099 catcher.

**Decisions:**
- **The user (2026-10-06):** the line domain's cylinder coordinates are (b, θ), with θ ∈ [0, π/2] the polar angle. In those coordinates the density 2 sin²θ · 2π is analytic. `directions_at` keeps μ_z, because a point's direction measure is uniform in μ_z and has no endpoint singularity.
- **The main agent (not a user ruling), on B5a at α = 1:** the block is assembled on the emission support only. Its columns cover the regions with Σ_t > 0, and `transport` takes the support as a region mask. This is the P1 sketch's item 6 (a void region never gets an emission column), brought forward. The full block on a void layer under a mirror stays a refusal row.
- **The main agent, on the turning-grading witnesses:** they move to hollow bodies with a small cavity, where a slot still starts at a small radius. ERR-099 is re-dropped after the move.
- **The main agent, on chunking:** a chunk is bounded by an estimate of its pieces (slots times the arc-length grading the optical length implies), not by its line count. The polar-angle rule removes the extreme lines the 8.2 GB run used, and the budget is measured at build.

## P1 step (b), third rung: built, the elegance review, and two more rulings (2026-10-06)

The production code is written on branch `feature/characteristic-rung3` (main agent): the kernel's `measure_density` and `LineDomain`, the even basis, `inflow(..., arriving=)`, `LineRule` with `GroupTransport`, and `WallCoupling`. `[M]` on the built code (`.venv/bin/python -O`, 16 points, 12/12 inner points):
- **closed-body conservation:** sphere 4.1e-13 (mirror and white), hollow sphere 1.7e-13, cylinder 3.2e-11, slab 8.9e-15 (white and periodic);
- **escape probability C9, wall transmission against P_ss, and reciprocity R = U D⁻¹ (white walls):** sphere and slab ≤ 1.1e-15 at τ 0.5 and 2; cylinder 2.3e-9 and 2.3e-13 at 16 and 32 points; reciprocity ≤ 6.6e-16.

**Decisions taken by the main agent (not user rulings):**
- **The block's columns are the emission support** (regions with Σ_t > 0 by default). Its rows are every panel, so R is computed directly and reciprocity is gated, not assumed (the elegance review withdrew its objection: U over every row would need void sources on trapped lines).
- **The piece budget is 4096 slots per traversal rule.** A slot costs up to about 0.5 MB in the Volterra block's inner rule: a single slab rule of 66 560 slots peaked at 36 GB, its chunks under 5000 at 2.25 GB.
- **From the elegance review, applied inline:**
  - `Wall` refuses specular + diffuse > 1, and a wall that is both (the user's LawSum ruling, now at the type);
  - the slab's grazing rule calls the basis's `graded_ends`;
  - `DirectionDomain` and `LineDomain` share one bounds table and one validator;
  - one stacked `inflow` call (emission columns, then a unit current at each diffuse wall) gives K_line, R, U and T, with D named where it is used;
  - `LineDomain`'s fold table carries a SCOPE-BOUNDARY tag;
  - `GroupTransport` checks its shapes.

**Rulings by the user, 2026-10-06:**
- **`compute_areas_1d` is a twin of `measure_density`:** filed as #584 and retired separately, because it moves bits on the sphere in production SN.
- **The near-void conditioning of the white block is fixed in rung 3.** I − αT was formed by subtraction, so rounding grew as 1/absorption: conservation off by 1.7e-7 at Σ_t = 1e-9 and 1.3e-4 at 1e-12; on a body absorbing nothing the computed T missed 1 by 1.1e-16 and the block came out at 2e16 rather than refusing. The fix:
  - `WallCoupling` carries the loss ℓ, the absorbed and leaked fractions of each injected current (by `expm1`, no subtraction);
  - the diagonal of I − αT is (1 − α) + α(ℓ + the off-diagonal column sum), and the last row is replaced by the balance row (1 − α) + αℓ;
  - the lossless case is exactly singular and refused in `update`; the input predicate first placed in `LineRule.transport` retired.

  Measured: conservation is flat from Σ_t = 1 to 1e-12 (white sphere 5.4e-13; hollow sphere white/white and mirror/white 1.8e-14; slab white/white 1.5e-14). With the loss-formed diagonal alone the two-wall cases still reached 7.7e-6 and 2.2e-5 at 1e-12.

**Open for rung 4 (the elegance review's question):** each group's default support is Σ_t > 0, so groups can have different column sets; the pencil needs one shared emission support (the union over groups, or the regions with any scattering, fission or source).

## P1 step (b), third rung: qa's line-measure defects, and the rules derived from the optical scale (2026-10-06)

**qa found three defects no gate saw** (`scratch/characteristic_architecture/p1_step_b3/qa/`). Each was a line-measure rule with a fixed resolution that did not follow the optical scale:
- **The slab's grazing cosine.** It used a fixed 12 halvings. Near void the collision probability was off by 4e-2 to 5e-1; on thin graded panels a block entry was off by 21%.
- **The cylinder's polar angle.** Plain Gauss: 3.2e-4 at τ = 0.01.
- **The impact parameter** (`chord_quadrature`). It was ungraded toward a thick rim: the sphere's P_ss was off by 4.6e-3 at τ = 100 and 7.3e-1 at τ = 1000. A small first radius left conservation off by 1.1e-8 at 16 points.

Closed-body conservation could not see any of them, because it holds line by line. The closed-form gates covered τ from 0.5 to 8 only.

**Ruled by the user:** derive every grading from the group's optical scale, by the law the traversal rule follows along a line. The alternative, refusing the regimes, was declined. `LineRule` is now per group: `LineRule.of(basis, walls, sigma_t, points, chunk=512, budget=256)` and `.transport(points, inner_points, support)`. `grazing_layers` is gone.

**As built (main agent):**
- **Grazing (the slab's μ on both signs, the cylinder's θ).** The coordinate is halved toward 0 down to the thinnest absorbing panel's optical width. Margins of 0 to 20 extra halvings changed nothing, to 5e-16 at τ from 0.5 to 2.3e-12, so there is no margin (`margin/study.py`).
- **Impact.** Each panel [r_k, r_{k+1}] is integrated in the chord half-length y = √(r_{k+1}² − b²), with three gradings:
  - hp toward y = 0, down to the next radius's branch point, at imaginary distance √(r_{k+2}² − r_{k+1}²);
  - hp toward the lower end, down to the image of b = 0;
  - the exponential ends of the traversal rule (2^j mean free paths of the thickest panel from here outward).

  The panel touching b = 0 takes plain Gauss on its lower half (smooth in b²). A hollow body gets a void impact panel [0, r_0].
- **What qa's small-radius defect actually was.** The cause was not b = 0. The wide middle panels saw the next radius's branch point just past their upper end, at a real distance of 0.07 but much farther in y. Gauss in b on [0.52, 0.88] left 1e-6 at 8 points.
- **The piece budget is 256 slots.** It is both faster and smaller: a two-region white cylinder at 8 points (8512 lines) took 126 s and 7.9 GB at 4096, and 32 s and 0.6 GB at 256.

**Measured after** (`-O`, references in mpmath at 60 digits with the subtractions done there):
- the sphere's P_ss: 2.1e-13 at τ = 100 and 6.4e-11 at τ = 1000 (16 points);
- a 1e-3 cavity: conservation 8.8e-11 at 8 points and 2.9e-13 at 16;
- the cylinder against Bickley: 1.4e-12 and 1.9e-15 at τ = 0.5 (8 and 16 points), 2.9e-14 at τ = 0.01 (8 points);
- the near-void slab: ≤ 5e-16 at Σ_t down to 1e-12;
- conservation at 8 points: sphere 5.6e-11, cylinder 2.1e-10, slab 7.3e-15; the near-void sweep is flat.

**A lesson from the measurement itself.** A reference evaluated in mpmath and subtracted in floating point (1 − P_esc near void) reported a defect of 7.6e-7 that the code did not have. Evaluated at the default 15 digits, it reported 100%. The subtraction belongs in mpmath.

**Not fixed, and recorded:**
- qa's F4: at Σ_t = 1e-315 the block overflows. That is the flux itself, about 1/Σ_t, exceeding double precision, not a rule defect.
- A slab placed at x ≈ 1e6 loses digits to absolute positions: conservation 5.6e-10.

**Re-review after the fix** (2026-10-06):
- **qa confirmed F1 to F3b fixed on its own probes.**
- **A thick slab's transmission was under-resolved (the test-architect).** At 8 points it was off by 4.1e-11 at τ = 8 and 1.8e-4 at τ = 30. In s = 1/v the attenuation e^{−τs} is a layer at the normal direction. The margin study read only τ ≤ 2.3, so it was blind to this.
- **qa's N1:** stopping the grazing halving at τ_min left a block entry off by 2.0e-6.
- **qa's N2:** a cancellation in hi − span refused a cavity of 1e-12.

**Main agent's fixes:**
- **The direction rules grade toward the normal direction too.** They use the exponential ends at 2^k over the body's normal optical depth, on s ∈ [1, ∞).
- **The grazing halving reaches τ_min / 64,** with 64 = `transport.VANISHING_DEPTH`, where e^{−64} vanishes in double precision. It is the traversal rule's own constant.
- **The two cancellations are removed:** the distance is lo²/(hi + span), and b = √(lo² + (span − y)(span + y)).

The slab now matches its closed forms to rounding at 8 points from τ = 2.3e-6 to 100: τ = 8 is 1.2e-15 and τ = 30 is 2.2e-15.

**Cost:**
- **The deeper grazing grading makes chunking matter.** The lines are ordered by projected speed and the budget is 1024. At budget 256 the τ = 0.01 cylinder split into single-line rules and did not finish in 25 minutes; now it takes 57 s and 1.3 GB.
- **A thick cylinder (τ = 30) costs about 1e6 pieces,** 1117 s under load. Its 18 432 lines are the tensor product of the b and θ rules. Filed as #586 (performance). The thick-cylinder rows run at 16 points.

**The last review round and the cost ruling** (2026-10-07):
- **The elegance re-review.** The hp law ("halve until each piece is no wider than its distance to the singularity") was spelled three ways, and the traversal rule's branch grading stopped one halving short. All three gradings (geometric, exponential, hp) moved into `characteristic/grading.py`, which every caller imports. The along-line grading now meets the stated law.
- **R per measured current, reverted.** Dividing R by its own injection tally made the wall area cancel out of the block; the gate on the one density caught it. R stays per nominal unit current (D = A_w/4), and the asymmetry is stated.
- **The slow tier could not run as written.** The three-region test cylinder (radius 2) takes 723 s per block at 8 points, for 38 400 lines (240 × 160). It conserves to 5.4e-13. At 16 points it took over 15 minutes. The slow tier's 18 cylinder rows would have cost hours.

  **Ruled by the user:** rung 3 lands with the slow cylinder rows cut to one fixture per law at 8 points; #586 is the next step, before rung 4. The re-pointing step (d) will need 2-group cylinder references, about 25 minutes each at today's cost.
- **A19, the hp toward the next radius, was not blind.** The test-architect's fixture did not reach it. A hollow sphere (0.2, 0.3, 1.0) graded 10 layers conserves to 9.2e-11 with the grading and 1.7e-3 without, because a wide panel precedes a thin one. A row was added at that input.

## P1 step (b), third rung landed (2026-10-07)

**Commits:** `71a207fa` (code and gates) and `4c996ca1` (theory, labels, ERR-101 to ERR-103, markers). Read the merge status from git, not from this line.

**What exists:**
- **Kernel:**
  - `MeasureCoordinate.derivative`;
  - `CoordSystem` / `Chart.measure_density`;
  - `Chart.line_domain()` returning a `LineDomain` (`LineShape`: IMPACT, IMPACT_POLAR, COSINE), sharing `_AXIS_BOUNDS` and `_in_box` with `DirectionDomain`.
- **Package `orpheus/derivations/continuous/characteristic/`:**
  - `PanelBasis.even`;
  - `LinePeriod.inflow(..., arriving=)`;
  - `WallCoupling(response, escape, transmission, loss, diffuse)` with `returning` and `update`;
  - `LineRule.of(basis, walls, sigma_t, points, chunk=512, budget=1024)` and `.transport(points, inner_points, support)`, which returns a `GroupTransport(line, coupling, support)` with `.block`;
  - `grading.py` (`graded_ends`, `exponential_ends`, `halvings`, `VANISHING_DEPTH`);
  - `Wall` refuses returning more than it receives, and a wall that is both specular and diffuse.

**Evidence** (`[M]` 2026-10-07):
- **Routine gates:** the seven characteristic and geometry gate files pass, 755 of 755 non-slow, in 7.5 minutes.
- **Slow tier:** about 62 minutes.
- **CI gate set:** 1422 passed locally. The catalogue reconciliation and the V&V audit pass.
- **Battery:** 54 arms; 52 redden; 2 are declared blind (N1, masked by the balance row; T2, whose catcher is in rung 2). The table is `scratch/characteristic_architecture/p1_step_b3/gates/battery/battery_table.md`.
- **Build and graph:** Sphinx `-E -W` is clean; `dead_references` reads 0 of 66 (a planted control read 1 of 67); pyright reports 0 errors.

**Filed:** #584 (`compute_areas_1d`, the area twin), #585 (a slab far from the origin), #586 (the cylinder's tensor-product cost).

## ⏸ COMPACTION POINT — 2026-10-07, rung 3 merged, #586 next

**State.** Rungs 1 to 3 of P1 step (b) are merged; read the state from git.

**Read in order:**
1. "third rung landed" (above): what exists and its names.
2. The rung-3 sections from "the premises re-measured" on: the rulings, the refutations and the measurements.
3. `docs/theory/references/characteristic.rst`: the line rule, the assembly, the white walls' coupling.
4. Issue #586, and `orpheus/derivations/continuous/characteristic/assembly.py` (`LineRule.of`, `_impact_rule`, `_grazing_ends`).
5. The cost probes, `scratch/characteristic_architecture/p1_step_b3/cost/`. Run each with `.venv/bin/python -O` from the repository root:
   - `mr3cyl.py <points>`: the three-region cylinder's line count, time and conservation; the 723 s baseline at 8 points;
   - `cylext.py <budget>`: the time and peak memory for a budget;
   - `cyl30.py`: the line and piece counts at τ = 30;
   - `cylmem.py <budget>`: memory against the budget;
   - `thick.py`, `c9.py`: the closed-form checks (escape, transmission), with the references in mpmath at 60 digits;
   - `cons.py`, `band.py`: closed-body conservation, and the near-void sweep.

**The next step: #586** (the user's ruling of 2026-10-07: before rung 4). The cylinder's rule is the tensor product of the impact rule (240 nodes) and the polar rule (160), each graded for the whole body. A three-region cylinder of radius 2 takes 723 s per group block at 8 points (38 400 lines) and conserves to 5.4e-13. Candidates, none measured:
- a non-product (b, θ) rule, grading θ per impact panel by that panel's own optical scale;
- fewer points on the deeply graded pieces;
- vectorising the traversal rule over pieces of like cost.

The acceptance in #586: a τ = 30 two-group block in under a minute at the gates' resolution, with the closed-form rows unchanged. Then rung 4: S and F, the pencil (one shared emission support across groups is an open question), the questions.

**Lessons from this stretch:**
- **A conservation identity that holds line by line cannot see the line measure.** Every line-rule defect qa found was invisible to K·1 = W·1/Σ_t; closed forms at the extreme optical scales saw them.
- **A study that reads totals cannot rule on entries.** The margin study read collision probabilities, and missed both a 2e-6 entry error and the thick-slab layer at the normal direction.
- **A reference subtracted in floating point after the mpmath call reports defects the code does not have.** Subtract in mpmath.
- **A fixed piece budget meets a cost distribution that the grading changes.** Order the lines by cost and re-measure the budget whenever a grading changes.
- **A battery restarted by every production change multiplies its runtime.** Freeze production before the battery, and scope each arm to the rows it can redden.


## #586 landed — the cylinder's cost (2026-10-07)

**The premise was refuted** for the question "what makes the cylinder slow". `[M]` A profile at 8 points, on a sixteenth of the lines, put 78-82 % of the time in `PanelBasis.values`. Of those calls, 98 % came from the Volterra triangle's inner rule (`TraversalRule._attenuated`). That rule was padded to the batch's thickest stretch: 4.2× the live evaluations on the three-region cylinder and 6.1× at τ = 30. The fact this establishes: the line count of the tensor (b, θ) rule was never the leading cost.

**The user's ruling** (2026-10-07): fix the evaluation first and re-measure. The line rule was to be touched only if the τ = 30 two-group block was still over a minute. It was not.

**What landed** (`refactor/characteristic-line-cost`):
- `95d1a511`:
  - `_Packing` packs each line's live entries; `_pieces` and `_attenuated` share it.
  - `_LagrangeTables` are built once.
  - The rule is unchanged: the blocks are bit-identical, or differ by ≤ 2.7e-17 relative from summation order (qa).
- `1d94790e`:
  - The slow tier is restored per the user's ruling: the thick cylinder legs at τ = 30 and 100 at 16 points, which catch ERR-101, and the three-region cylinder's second group at 8 points. The tier is now 20 rows in 36 min 27 s.
  - Fixed `_flat`, which failed on empty reads.
- `3a852eb0`: the theory page, section `characteristic-cylinder-cost`.

**Times `[M]`** (one process per run, budget 1024, which is still the fastest):

| cylinder | before | after |
|---|---|---|
| τ = 30, per group, 8 points | about 215 s on an idle machine | 27.8 s |
| τ = 30, per group, 16 points | — | 102 s |
| τ = 0.01 | 57 s | 15.6 s |
| three regions, 8 points | 723 s | 175 s |
| three regions, 16 points | — | 792 s |

**Filed:** #587, the non-tensor rule. It is now only the lever for the three-region cylinder at 16 points.

**Lesson:** profile before choosing among an issue's candidates. All three of #586's candidates assumed that the line count was the cost; one cProfile run read the cost directly.

## ⏸ COMPACTION POINT — 2026-10-07 (second), #586 merged, rung 4 next

**State.** Rungs 1 to 3 of P1 step (b) and #586 are merged; read the state from git (`main` at the plan commit after `a1066516`; CI passed on `a1066516`).

**Next: rung 4.** S and F, the pencil, and the questions. Nothing about it is designed yet. It is a W3 surgical carve: an API sketch goes to the user (`AskUserQuestion`) before any code.

**Read in order:**
1. "P1 API sketch" (the P1 sketch's items 5 to 9): what the pencil and the questions are for.
2. "P1 step (b), third rung: the premises re-measured": its bullet "S and F, an emission-matrix helper split out of `_infinite_medium_matrices`, move to rung 4 with the pencil".
3. "P1 step (b), third rung: built, the elegance review, and two more rulings": the open question "**Open for rung 4**". Each group's default emission support is Σ_t > 0, so groups can have different column sets. The pencil needs one shared support: the union over groups, or the regions with any scattering, fission or source. This is the first question for the user.
4. The rulings ledger: Q4 of the third rung (S and F move to rung 4; the block carries Σ_t, so the pencil needs the emission on its own, not read back out of the 0-D loss matrix).
5. "P1 step (b), third rung: the verification spec and what it refuted": D11 (the 1-group Rayleigh–Ritz lower bound) needs the pencil and lands in rung 4.
6. `docs/theory/references/characteristic.rst`, the section "What the package does not compute", for the rung-4 items, and the section on the reference kernel, `verification-reference-kernel`: the dense pencil, `EigenPosing` and `SourcePosing`.
7. Code: `orpheus/derivations/continuous/characteristic/assembly.py` (`LineRule`, `GroupTransport`), and `_infinite_medium_matrices` (`orpheus/derivations/common/eigenvalue.py:41`).

**Cost to plan with** `[M]` (8 points, budget 1024):
- the three-region cylinder takes 175 s per group, so a two-group block takes about 6 minutes;
- the one-region cylinder at τ = 30 takes 27.8 s per group;
- a sphere takes under 1 s.

The characteristic slow tier is 20 rows in 36 min 27 s, and the non-slow selection 722 rows in 5 min. Run gates serially.

**#586's probes** are in `scratch/characteristic_architecture/p586/cost/`; run each with `.venv/bin/python -O` from the repository root:
- `prof586*.py`: cProfile on a sixteenth of the lines;
- `base586.py <out.npz>` with `cmp586.py <a.npz> <b.npz>`: blocks on five fixtures, for a before/after bit-identity check;
- `budget586p.py <tau|mr3> <budget> <points>`: full-rule time and peak memory;
- `pad586.py`: the inner rule's padding.

qa's and the elegance review's notes for #586 are in `scratch/characteristic_architecture/p586/`.

**Open issues:** #584 (the twin of `compute_areas_1d`), #585 (a slab far from the origin), #587 (the non-tensor rule).

**Lesson from #586:** profile before choosing among an issue's candidates. The issue's three candidates all presumed one cost, and the measurement named another.

## P1 step (b), fourth rung: API sketch (the main agent, 2026-10-07; checkpoint with the user before code)

Every item is a proposal until ruled. Rung 4 turns the per-group transport blocks of rung 3 into the multigroup Galerkin system and answers its questions in coefficient form. It builds nothing that reads a point (rung 5) and does not map the interface vocabulary's question values onto the system (the door, P1 sketch item 8), unless Q3 rules otherwise.

**The weak form.** With the flux φ = Σ_j φ_j u_j on every panel and every group, and the emission density q_g = Σ_g' S_g←g' φ_g' + (1/k) Σ_g' F_g←g' φ_g' + q_ext,g, the Galerkin equation per group is W φ_g = K_g q_g, with W the mass matrix and K_g the rung-3 block (rows on every panel, columns on group g's emission support). The cross sections are constant on a region and no panel crosses a region, so the emission of a basis function is a basis function times a constant: S and F act per node, exactly. Stacked group-major (the unknown is φ of shape (G, N), flattened):
- the mass W_G = I_G ⊗ W, (GN, GN);
- the transport K = blockdiag(K_g), (GN, ΣM_g);
- the scattering S and fission F per node, (ΣM_g, GN): row (g, i) for i in group g's support reads S[r(i), g, :] and F[r(i), g, :] at node i of every group.

**Why the 1-group bound (D11) holds on this form** `[R]` (the derivation goes to the theory page). In 1 group put M = σ_s + νσ_f/k, constant per region. On the emission support (M > 0) the equation W φ = K M φ, written for q = M φ, is W M⁻¹ q = K_ss q, because M commutes with W there (both are constant per panel). That is the Rayleigh–Ritz Galerkin form of the self-adjoint positive operator K_op in the weight M⁻¹, so the Galerkin λ(k) is below the true λ(k) for every k, and the Galerkin k is below the true k (λ decreases as k grows). Rows off the support are slaved (φ_v = W_vv⁻¹ K_vs q) and do not enter the bound. A column whose emission row is identically zero is multiplied by zero in K S, so it does not change the pencil.

**1. The emission and fission matrices (`derivations/common/eigenvalue.py`).** `group_emission(sig_s, nu_sig_f, chi, sig_2=None) -> GroupEmission(scattering, fission)`, split out of `_infinite_medium_matrices`: scattering = (Σ_s + 2Σ_2)ᵀ (to ← from), fission = χ ⊗ νΣ_f. `_infinite_medium_matrices` becomes (diag Σ_t − scattering, fission), bitwise. The (n,2n) multiplicity literal moves with it, so `tests/gates/transport/test_n2n_multiplicity_census.py`'s `_REFERENCE_LITERALS` row re-points from `_infinite_medium_matrices` to `group_emission`.

**2. The cross sections per region (`characteristic/cross_sections.py`).** `RegionCrossSections(total (n, G), scattering (n, G, G), fission (n, G, G))`, frozen, read-only arrays, shapes checked. `of(mixtures)` reads one `Mixture` per region: `SigT`, `SigS[0]`, `Sig2[0]`, `SigP` (νΣ_f), `chi`, through `group_emission`; it refuses a mixture with a non-zero higher Legendre order of `SigS` or `Sig2`, naming the region and the order (anisotropic scattering is out of P1's scope, SCOPE-BOUNDARY). `emission_support() -> (G, n) bool`: region r emits in group g iff scattering[r, g, :] or fission[r, g, :] is non-zero (Q1).

**3. The system (`characteristic/system.py`).** `GalerkinSystem(basis, cross_sections, groups: tuple[GroupTransport, ...])`, built by `of(basis, walls, cross_sections, points, inner_points, support=None)`: one `LineRule.of(..., sigma_t=total[:, g], points)` and one `.transport(points, inner_points, support[g])` per group (each group graded by its own optical scale, as ruled in rung 3). Its matrices `mass`, `transport`, `scattering`, `fission` (above), and two pencils, both on `DensePencil`:
- `pencil`: `DensePencil(W_G − K S, K F)`, the k question: `.fundamental()`, `.spectrum()` (the higher modes; the mode nearest τ is read from it), `.adjoint()` (the transposed forms; W and each K_g restricted to the support are symmetric by reciprocity, so the transpose is the weak-form adjoint);
- `emission_pencil`: `DensePencil(W_G, K (S + F))`, the source questions: `fixed_source(q) = emission_pencil.least_solution(K q)` and `response(r) = emission_pencil.adjoint().least_solution(W_G r)`, with q and r nodal coefficients of shape (G, N). Every secondary emission is in the gain, so `least_solution`'s subcriticality check covers a closed body made supercritical by (n,2n) with no fission. The k pencil's own `least_solution` would not: its gain is the fission alone, so it would solve that case directly and return a negative flux (the attacker's F5 trap).
- A source with a non-zero coefficient outside the emission support is refused, naming the group and the regions; the caller widens the support (`support=` at `of`).

**4. Cost** `[M]` (#586's table, 8 points): a two-group three-region cylinder is two blocks of 175 s; a sphere is under a second per group. The dense eigenproblem is (GN)², a few hundred unknowns at the gates' resolution, negligible beside the blocks.

**5. The gates (the test-architect re-specifies them onto these names, each with its first red):** D0 (S and F per region against the `Mixture` tables by hand), D11 (the 1-group bound, on a nested degree ladder p = 1 to 4 on fixed panels, which nests the spaces including the even panel), D5 and D5b (a closed body reads k_inf and a flat eigenvector with the 0-D group ratio), D2' (Sood's 2-group critical sizes), D13 (the adjoint), D14 (biorthogonality and response reciprocity), D7 (the pencil residual), E2 (a closed body with a uniform source), E3 (the fixed source F φ_k / k returns φ_k), E6 (the supercritical refusal), and the support rows of Q1.

**Questions for the user:**
- **Q1, the emission support** (the open question of rung 3). (a) Recommended: per group, the regions where that group's emission can be non-zero (a scattering or (n,2n) transfer into g, or fission with χ_g > 0), from `emission_support()`, widened by the posed source; the column sets differ by group. (b) One shared set, the union of (a) over groups. (c) Keep each group's Σ_t > 0. Against (c) `[R]`: a region transparent in group g that scatters or fissions into g has no column, so its emission is silently dropped. Against (b): it assembles columns whose emission is identically zero, and on a region void in g under a mirror such a column is a source on a lossless trapped line, refused by `TrappedSource` although the problem is well posed. The block-level default `support=None` (Σ_t > 0) stays for rung 3's block gates, where no emission exists; the system always passes its support.
- **Q2, the cross-section value.** (a) Recommended: `RegionCrossSections` as item 2, read from mixtures with the anisotropy refusal at that boundary. (b) The system takes the per-region `Mixture`s directly and reads them where it needs them.
- **Q3, the scope.** (a) Recommended: rung 4 answers in coefficients (a source and a detector are nodal coefficients); projecting a mesh-free source onto the basis and mapping `Eigen`, `FixedSource` and `Response` onto the system belong to rung 5 with the door and the point reading. (b) Rung 4 also maps the question values, with the Galerkin projection of a mesh-free source.
- **Q4, the name of the system.** (a) `GalerkinSystem`; (b) `TransportSystem`; (c) `MultigroupSystem`.

**Ruled 2026-10-07** (see the ledger), all four as recommended: Q1 the exact per-group emission support, widened by the posed source; Q2 `RegionCrossSections`; Q3 coefficients only, the question values and the projection of a mesh-free source go to rung 5; Q4 the name `GalerkinSystem`. Next: the test-architect re-specifies the rung's gates (item 5) onto these names, each with its first red, while the main agent writes the code on `feature/characteristic-rung4` in the order items 1, 2, 3.

## P1 step (b), fourth rung: the flux-form adjoint refuted (2026-10-07)

**What was measured** `[M]` (`scratch/characteristic_architecture/p1_step_b4/main/smoke.py`, `qform.py`, `-O`; a 2-group, two-region mirror sphere with upscatter, degree 3, 8 line points, 12/12 along each line):
- **The forward k pencil on the flux, `(W_G − K S, K F)`, is right.** Its k equals k_inf to 1.0e-14, its eigenvector is flat to 2.2e-10 with the 0-D group ratio, and the 1-group bound holds on a vacuum Sood sphere (Ua-1-0-SP at its critical radius): k − 1 = −2.5e-5, −2.1e-6, −4.9e-7, −1.4e-7 at degree 1 to 4, negative and increasing.
- **Its transpose is not the adjoint flux. Refuted:** item 3's claim that "the transposed pencil's eigenvectors are the adjoint flux's coefficients". The transposed eigenvector is flat, but its group ratio is 2.1466, against 1.0733 for the adjoint flux A^(−T) νΣ_f: a factor Σ_t,2/Σ_t,1 = 2.
- **The reason.** The flux-form equation is φ = 𝒦 E φ, with 𝒦 the transport operator and E the emission. Its adjoint is ψ = E* 𝒦 ψ. That ψ is not the adjoint flux φ†; it is E* φ† (in an infinite medium, Σ_t φ†). The adjoint flux is the adjoint of the emission-form equation q = E 𝒦 q, which is φ† = 𝒦 E* φ†.
- **The emission form gives φ† by transposition.** The unknown is the emission density q on the supports. The Galerkin form is W_s q = (S + F/k) K q, with W_s = blockdiag(W restricted to each group's support). Measured on the same sphere:
  - its k equals the flux form's to 4.4e-16;
  - the flux it gives, φ = W_G⁻¹ K q, equals the flux form's mode to 3.2e-15;
  - its transposed eigenvector is flat with the ratio 1.07328385899813, against A^(−T) νΣ_f = 1.07328385899814.
- **The emission form serves the other questions too.**
  - **Source:** W_s q = (S + F) K q + W_s q_ext, then φ = W_G⁻¹ K q.
  - **Detector:** the transpose with K^T r gives φ†, the importance of a source, and the reading is ⟨r, φ(q_ext)⟩ = φ†^T W_s q_ext.
  - **Where each answer lives:** φ† is defined exactly on the supports, which is where a source can be posed (the ruled Q1). The flux is defined everywhere.
  - **Size:** the pencil is ΣM_g, not GN.

This reverses the P1 sketch's item 6 (ruled 2026-10-06: "the unknown is the FLUX coefficients phi"), which had corrected the spec's original "emission coefficients c". The user rules it.

**Ruled 2026-10-07 by the user:** the unknown is the emission density q on each group's support (the emission form). The flux is read as W_G⁻¹ K q; the transposed pencil gives the adjoint flux φ† on the supports. This supersedes the P1 sketch's item 6 on the unknown.

## P1 step (b), fourth rung landed (2026-10-07)

**Commits:** `42222c47` (code and gates) and `d2950bde` (theory). CI `gates` passed on `d2950bde`. Read the merge status from git, not from this line.

**What exists:**
- `group_emission` and `GroupEmission` in `derivations/common/eigenvalue.py`. `_infinite_medium_matrices` reads them and is bitwise unchanged.
- `RegionCrossSections` (`characteristic/cross_sections.py`): `of(mixtures)`, which refuses anisotropy, and `emission_support()`, returning `(n, G)`.
- `TransportResolution(line_points, points, inner_points)` in `assembly.py`.
- `EmissionSpace(supports, region)`, with `restriction`, `restrict` and `split`.
- `GalerkinSystem(basis, walls, cross_sections, resolution, source_regions)`, in `characteristic/system.py`:
  - its blocks (`groups`) are derived from its fields;
  - its matrices: `mass`, `emission_mass`, `transport`, `scattering`, `fission`;
  - its two pencils: `pencil`, `source_pencil`;
  - its questions: `flux`, `fixed_source`, `response`.

**Decisions taken in the build (not user rulings):**
- **From the elegance review:**
  - The blocks are derived, not stored. Measured: blocks handed in from another mixture gave k = 0.108657 against 0.107926, with no error.
  - The emission layout is one object, `EmissionSpace`.
  - The masks are region-major.
  - `emission_pencil` is renamed `source_pencil`.
- **From qa:** the refusal names the regions again, and `extend` was removed because nothing used it.
- **From the test-architect:**
  - Detector reciprocity (D14(ii)) holds for any transport matrix, so it is documented as a transposition check only. D13(iv), the adjoint against the forward solve of the group-transposed problem, carries the adjoint's physics.
  - An under-integrated line rule cannot redden D11, because it moves k without breaking the bound. D11 reds on a non-reciprocal transport matrix instead.
  - D11b runs on ungraded panels. On graded panels the bound's margin falls inside the truth's resolution.

**Evidence** `[M]` 2026-10-07:
- **Gates:** `test_characteristic_system.py` has 80 rows: 78 routine (about 78 s) and 2 slow cylinder rows.
- **Battery:** 23 arms plus controls (`scratch/characteristic_architecture/p1_step_b4/battery/verdicts.md`). The no-arm control reds 0 of 78; the positive control reds 51 of 78. Every row reds under at least one arm.
- **Touched trees:** 3157 passed. The CI set passed locally (627) and on CI.
- **Build and graph:** Sphinx `-W` is clean, and `dead_references` reads 0 of 66.

**Filed:** #588, Sood's UAL-2-0 critical sizes, which are off beyond their printed digits. The D2' rows exclude those two cases.

## ⏸ COMPACTION POINT — 2026-10-07 (third), rung 4 merged, rung 5 next

**State.** Rungs 1 to 4 of P1 step (b) are merged; read the state from git.

**Next: rung 5,** a W3 surgical carve: an API sketch goes to the user before any code. It covers three things:
- **The door:** P1 sketch item 8, `CharacteristicDerivation` posed from a `GeometrySpecification` with a frozen `Resolution`. That `Resolution` now contains the basis's degree and grading, plus a `TransportResolution`. The anisotropy refusal stays at `RegionCrossSections.of`, which the door calls, so there is one door.
- **The mapping of the question values onto `GalerkinSystem`:**
  - `Eigen` with `Fundamental` or `Nearest`;
  - `FixedSource`, with a Galerkin projection of the mesh-free source and its source regions;
  - `Response`.
- **The reading at a point:** P1 sketch item 7, with the spec's C6 and C8.

**Read in order:**
1. The P1 sketch, items 7 and 8.
2. The fourth rung's sections above: its sketch, the refuted flux-form adjoint, and "landed".
3. `docs/theory/references/characteristic.rst`, sections `characteristic-galerkin-system` and "What the package does not compute".
4. `orpheus/numerics/question.py`, the question values.
5. `orpheus/derivations/continuous/characteristic/system.py`.

**Cost to plan with** `[M]` (8 points): a sphere system takes seconds; a two-group three-region cylinder takes about 6 minutes.

**Open issues from this campaign:** #584, #585, #587 and #588. #588 is Sood's UAL-2-0 truths; the D2' rows return when it is settled.

**Probes from rung 4,** under `scratch/characteristic_architecture/p1_step_b4/`. Run each with `.venv/bin/python -O` from the repository root:
- `main/smoke2.py`: the 1G bound ladder, the closed 2G k, flux and adjoint, E2, and reciprocity, on the final API;
- `main/qform.py`: the flux form against the emission form. It uses the pre-refactor API and is kept as the record of the measurement;
- `spec.md`: the rung's verification spec;
- `ta/`: the test-architect's probes;
- `battery/`: the battery plugin and `run.sh`, with `verdicts.md`;
- `qa/`: qa's probes, including the transposed-problem adjoint route;
- `elegance.md`: the elegance review.

**Lessons from this stretch:**
- **Transpose the equation whose adjoint you want.** The weak-form transpose of φ = 𝒦Eφ is E*φ†, not φ†. Smoke-test the adjoint against the 0-D importance A⁻ᵀνΣ_f before building on it.
- **A stored derived field is a door for an inconsistent value.** Blocks handed to the system as a field gave a wrong k in silence. Derive what the posing determines.
- **A reciprocity identity of the form ⟨r, L⁻¹q⟩ = ⟨L⁻ᵀr, q⟩ is true for any L.** The non-tautological check of an adjoint is an independently posed forward problem: the group-transposed one.
- **Do not refactor production while reviewers and a battery are reading it.** The elegance fixes landed under qa and the battery, and both had to re-run. Collect the reviews first, then change the code once.

## P1 step (b), fifth rung: API sketch (the main agent, 2026-10-07; checkpoint with the user before code)

Every item is a proposal until ruled. The surface it plugs into was mapped by the explorer at `4e375fe6`: `scratch/characteristic_architecture/p1_step_b5/surface.md`, every claim with its `file:line`. The facts the sketch rests on:
- **The specification holds ONE question.** `GeometrySpecification(materials, geometry, question)` (`orpheus/specification/specification.py:126`) keeps only the materials the geometry uses and canonicalises the question's keys.
- **A derivation answers that question for every observable.** `Derivation` is a Protocol with one method, `evaluate(observable) -> Evaluation` (`orpheus/reference/solution.py:63`). `evaluate` receives `Eigenvalue`, `FluxIntegral(weight)` or `PointValue(position, group)`; `ReferenceSolution` splits a `Ratio` before `evaluate` sees it, and `admit_observable` has already checked the observable against the problem.
- **The two existing derivations** (`ExactInfiniteMediumDerivation`, `TrajectoryResolventDerivation`) check the question at construction, dispatch on the observable alone, and serve one question only: `Eigen(CellCoefficient.every(Channel.FISSION_EMISSION).resolve(materials))` at the physical point with the fundamental mode. No reference yet serves `Nearest`, an offset point, `FixedSource` or `Response` (`[M]` the explorer: 0 hits outside the specification).
- **The traced memo** rebuilds the receiver from its constructor fields in a fresh interpreter, so the receiver and every field are `ContentIdentity`. No characteristic class is one today.
- **A source and a detector are mesh-free functions:** a `RegionwiseConstant` table `(regions, groups)` on the angle-integrated space, or a `Symbolic` expression per group in `(r, mu, phi)`, a density over dΩ. A source enters phase space through the angular section (its rate kept); a detector through the retraction's adjoint (`orpheus/numerics/mesh_free_function.py:13-30`).
- **Nothing projects a mesh-free function onto a `PanelBasis`,** and nothing reads the transported flux at a point: `TraversalRule.angular_flux`, `LinePeriod.inflow` and `Chart.directions_at` exist, and only tests call them; the diffuse walls' return is assembled only inside `LineRule.transport`.
- **The eigen gauge differs between the two existing references.** The exact homogeneous one fixes ⟨νΣ_f, φ⟩ = 100 (production's gauge, `exact_homogeneous.py:53`); the trajectory resolvent divides by its last fission rate.
- **The chart of a parameter is not minted** (`specification.py:51`, owed to #529). Both existing references return k for `Eigenvalue()`, so k is the chart every consumer reads today.

**1. The resolution (`characteristic/resolution.py`, or in `reference.py`).** `Resolution(degree, layers, ratio, transport: TransportResolution)`, a frozen `ContentIdentity`; `TransportResolution` becomes one too. One value for the solve and, in the reading's sub-rung, the direction rule at a point (a field added then, which changes the schema tag as it should). It dissolves #516's two answers: one resolution, one transport.

**2. The projection (`basis.py`).** `PanelBasis.project(table) -> (G, N)`: the L2 projection c = W⁻¹⟨u_i, f_g⟩ of a function given by its values on each panel's Gauss points, panel by panel (W is block-diagonal by panel). Exact for any per-region polynomial of degree ≤ p, so exact for a `RegionwiseConstant`. The mesh-free function is read into that table at the door, where its role is known:
- **a source** keeps its rate: a `RegionwiseConstant` as given; an isotropic `Symbolic` q(r) is a density over dΩ, so its rate is ∫ q dΩ = 4π q `[R]` (the module docstring's section; the 4π is the angular measure, not a convention). A `Symbolic` source depending on μ or φ is refused: an anisotropic source needs its own first-flight transport (SCOPE-BOUNDARY, the same machinery as the anisotropic emission);
- **a weight or a detector** is the scalar flux's weight, as given.

Because W is block-diagonal by region and every support is a union of regions, restricting the projection to the support and applying W_s equals restricting the load ⟨u_i, f⟩: the system's coefficient API (rung 4) stays unchanged. The source's regions (the support widening, rung 4's Q1) are the regions where a group's projected coefficients are non-zero.

**3. The door (`characteristic/reference.py`).** `CharacteristicDerivation(specification, resolution)`, a frozen `ContentIdentity`. Construction is the one door, each refusal a SCOPE-BOUNDARY or a ValueError naming what it refused:
- an infinite medium (no geometry);
- a mixture with anisotropic scattering or (n,2n) emission, through `RegionCrossSections.of` (rung 4: one door for that refusal);
- `PrescribedInflow` and a wall both specular and diffuse, through `Walls.of` (rung 1);
- an `Eigen` along a parameter other than the fission emission, or at an offset point;
- an anisotropic `Symbolic` source.
The cross sections are `RegionCrossSections.of([materials[m] for m in geometry.mat_ids])`, one region per interval, as `ConcentricPartition.of(geometry)` reads them. The basis is `PanelBasis.of(partition, degree, layers, ratio)`.

The question is resolved once, at construction, into its answer (a `cached_property`, solved on the first `evaluate`); `evaluate` then dispatches on the observable alone, as the two existing derivations do:
- **`Eigen`:** the system with no source regions; `Fundamental` reads `pencil.fundamental()`; `Nearest(tau)` the eigenvalue of `pencil.spectrum()` nearest tau in the k chart (Q4). The flux is `system.flux(q)` in the gauge of Q3.
- **`FixedSource(source)`:** the system posed with the source's regions; the emission is `source_pencil.least_solution(...)` (`fixed_source`), which refuses a supercritical body (E6).
- **`Response(detector)`:** Q2.

The observables, on every answer:
- `Eigenvalue()`: k (eigen answers only; `admit_observable` already refuses it elsewhere);
- `FluxIntegral(w)`: ⟨w, φ⟩ = ⟨P w, φ_h⟩, with P the projection of item 2 and φ_h the Galerkin flux coefficients, read as the pairing `(W c_w)ᵀ φ_h`, no point loop. Since W φ_h = K q, φ_h reproduces every moment of the transported flux 𝒦q against the basis, so the reading is exact for a weight in the basis space (every `RegionwiseConstant`) and carries the weight's projection error otherwise, which falls with the degree `[R]`. The alternative, a volume rule over point readings, waits for the reading;
- `PointValue(position, group)`: the reading at a point (Q1).

`characteristic_reference(specification, resolution) -> ReferenceSolution`, uncertified (`Uncertified`, the certificate waits for P4, #566), the same consumer surface as `trajectory_resolvent_reference`. The `evaluate` is a `@traced_memo`.

**4. The reading at a point** (P1 sketch item 7, ruled 2026-10-06). The scalar flux at x is the transport of the converged emission over `Chart.directions_at(x)`: per direction, the line through x, its chord and period, the attenuated emission up to x (`TraversalRule.angular_flux`), the specular return through `LinePeriod.inflow`, and the diffuse walls' re-emission, which no point reading assembles yet (the white walls' partial currents from `WallCoupling`, transported from the wall to x). The direction rule at x is split where the line's crossing set changes (the tangencies to every interface) and graded near grazing. Gates: C6 (one transport, two test measures: the volume integral of the reading against u_i equals the Galerkin entry), C8 (reciprocity through readings), C7's reading leg, E1 (Garcia Case 1 per point), E4.

**Questions for the user:**
- **Q1, the scope of this rung.** (a) Recommended: split it. 5a is the resolution, the projection, the door and the observables that need no point (`Eigenvalue`, `FluxIntegral`), with `PointValue` refused as a SCOPE-BOUNDARY naming 5b; 5b is the reading at a point, with C6, C8 and E1. The reading is the largest new machinery of P1 (the diffuse return at a point has no assembly today) and earns its own sketch and review. (b) One rung, all of it.
- **Q2, the response's answer.** (a) Recommended: the group-transposed forward problem. For isotropic emission the transport 𝒦 is self-adjoint (reciprocity, the property C11 gates), so the adjoint flux of (Σ_t, S, F) with detector r is the forward flux of (Σ_t, Sᵀ, Fᵀ) with source r: `RegionCrossSections.transposed()` (the adjoint cross sections, `[to, from]` swapped), and the detector posed as that problem's source. Every observable then reads the importance with the forward machinery unchanged: the flux integral everywhere, and in 5b the point reading. Its supports are the adjoint emission supports (where something scatters or fissions out of a group), which is where the importance has an emission. D13(iv), which the test-architect built as this route against `system.response`, then compares the door with the transposition: still two independent routes. (b) `system.response` on the forward system (rung 4's transposition): φ† is known only on the forward emission supports, so a flux integral or a point reading of the importance outside them needs a further transport of the adjoint emission r + (S + F)ᵀφ†, whose support is not the forward one.
- **Q3, the gauge of an eigen flux** (a `FluxIntegral` or a `PointValue` of an eigen answer; a `Ratio` does not depend on it). (a) Recommended: ⟨νΣ_f, φ⟩ = 100 over the body, production's gauge, which the exact homogeneous reference already uses: a closed homogeneous body then reads that reference's flux, a gate across references. (b) Unit total fission production. (c) The trajectory resolvent's (its last fission rate), which no other reference shares. The trajectory resolvent's gauge would remain different until it retires.
- **Q4, `Nearest(tau)`.** (a) Recommended: served, with tau read in the k chart, the chart `Eigenvalue()` is read in by both existing references; the chart is stated once in the door, with #529 named as the owner of parameter charts. (b) Refused as a SCOPE-BOUNDARY until #529 mints the chart of a parameter.

**Not in this rung:** the eigen adjoint as an observable (no observable reads it; a contract change owed to #529); the migration of the SN consumers (P1 step (c) onwards).

**Ruled 2026-10-07** (see the ledger), all four as recommended. Next: 5a on `feature/characteristic-rung5a`; the test-architect specifies its gates while the main agent writes items 1 to 3.

## P1 step (b), rung 5a: the reviews and the fixes (2026-10-07)

**Reviews** (`scratch/characteristic_architecture/p1_step_b5/`): `elegance.md` (three violations, four concerns), `qa.md` (two high, two medium, three low; probes `qa/p1`-`p8`). Both found the higher-mode gauge independently, and so did the test-architect (`ta/probe_n4.py`): a production gauge divides by rounding on a mode whose net production is zero.

**Fixed in one pass, after every review had returned:**
- `RegionCrossSections` stores `spectrum` (χ) and `production` (νΣ_f) and derives `fission`; `transposed()` exchanges them, so the adjoint's production is the forward spectrum (elegance V3, qa F7: `production` had been recovered from the outer product and was wrong on the transposed set).
- A `RegionwiseConstant` is read onto the nodes exactly (`PanelBasis.on_nodes`); `Resolution` refuses `source_points < 2(p + 1)` and a non-integer count (elegance V2, qa F4, F5). `panel_mass` and `project` share one panel rule (elegance C7).
- The answer is typed per question (`_FundamentalAnswer`, `_ModeAnswer`, `_SourceAnswer`); `evaluate` matches on the observable and the answer (elegance C4). The role is an arrow (`_section`, `_pullback`, or none for a weight), not a flag (elegance C5).
- `Nearest` excludes the null-fission cluster: a mode is kept when ‖F K v‖ exceeds the rank tolerance of F K (qa F1: `Nearest(0)` had served a rounding eigenvalue of 1e-20).
- `PRODUCTION_GAUGE` withdrawn with the gauge correction.

**Measured after the fixes** `[M]` (`main/smoke3.py`, `-O`): a closed A sphere's mean flux equals the exact flux /(100 V) to 5e-15; reciprocity ⟨Σ_d, φ(Q)⟩ = ⟨Q, Rψ†⟩/4π to 2e-16; a symbolic detector equals its table to 2e-16; `Nearest(0)` on the vacuum A|B sphere answers 2.3e-4.

**Not fixed here:**
- `GalerkinSystem.source_pencil` spelled beside `pencil` rather than as two splittings of one operator: filed as #589.
- The symbolic evaluator duplicated with the trajectory resolvent's `_radial_weight`: not filed, because the twin retires with that family in P1 step (e).
- Production's `ScaleGauge` refuses only an exactly zero functional (qa F1, the same defect class): not filed. Production gauges only fundamental modes, whose production is positive (Perron–Frobenius), so no realizable state reaches it today (X1); it becomes a defect when production gauges a higher mode.

## P1 step (b), rung 5a landed (2026-10-08)

**Commits:** `fe977a90` (code and gates) and `af286805` (theory). Read the merge status from git, not from this line.

**What exists:**
- `CharacteristicDerivation(specification, resolution)` and `characteristic_reference` (`characteristic/reference.py`). The door refuses at construction. Its answers are typed per question: `_FundamentalAnswer`, `_ModeAnswer` (eigenvalue only) and `_SourceAnswer`. `evaluate` is a traced memo that answers `Eigenvalue` and `FluxIntegral`; `PointValue` is refused, naming 5b.
- `Resolution(degree, layers, ratio, transport, source_points)`, with `source_points ≥ 2(p + 1)`.
- `PanelBasis.on_nodes`, `PanelBasis.project` and `_panel_rule` (shared with `panel_mass`).
- `RegionCrossSections(total, scattering, spectrum, production)`, with `fission` derived and `transposed()` exchanging the spectrum and the production.
- **The declared gauge.** `Eigen.gauge` is resolved by `_canonical_gauge` (`orpheus/specification/specification.py`). The default is `EIGEN_GAUGE = every(FISSION_EMISSION, N2N_EMISSION)`; a problem that produces nothing keeps None. `production_emission` (`derivations/common/eigenvalue.py`) is the references' channel emission, with its own (n,2n) literal registered in the census. The exact infinite medium reads its flux at the declared production density 100; the trajectory resolvent accepts the fission key only.

**Evidence** `[M]` 2026-10-08:
- **Gates:** `test_characteristic_reference.py` has 95 rows. The two batteries (39 arms and 9 gauge arms, `p1_step_b5/ta/battery*`) redden every arm.
- **Test runs:** the CI set passes 1425. The touched trees (derivations, numerics, specification, reference, homogeneous, data, and the SN question rows), excluding slow, give 8303 passed, 246 skipped, 2 xfailed.
- **Build and graph:** Sphinx `-W` is clean, and `dead_references` reads 0 of 66.

**Open, filed:** #589 (the two splittings of one operator). SN and the homogeneous solver read the declared gauge in #517: until then a test judging the homogeneous solver declares the fission gauge.

## ⏸ COMPACTION POINT — 2026-10-08, rung 5a merged, 5b (the reading at a point) next

**State.** Rungs 1 to 4 and 5a of P1 step (b) are merged; read the state from git.

**Next: rung 5b, the reading at a point,** P1 sketch item 7 (ruled 2026-10-06). It is a W3 carve: an API sketch goes to the user before any code. What it must compose, per the explorer's map (`p1_step_b5/surface.md` §6):
- `Chart.directions_at(x)`, whose result has a density and the methods `direction()` and `impact_parameter()`;
- the line through x, with its chord and its `LinePeriod`;
- `TraversalRule.angular_flux(t, inflow)`, the emission attenuated up to x, with `LinePeriod.inflow` for the specular return;
- the diffuse walls' re-emission at x, which nothing assembles today: the white walls' partial currents (`WallCoupling`) transported from each wall to x.

The direction rule at x is split at the tangencies to every interface and graded near grazing. `PointValue` then replaces `_refuse_point_reading`, and the `FluxIntegral` of a symbolic weight could alternatively be read by a volume rule over readings.

Gates owed:
- C6: one transport, two test measures. The volume integral of the reading against u_i equals the Galerkin entry, and one constructor builds both direction rules.
- C8: reciprocity through readings.
- C7's reading leg.
- E1: Garcia Case 1 per point.
- E4: a pure absorber against C1's closed forms.

**Read in order:**
1. The P1 sketch, item 7.
2. This rung's sketch, its reviews, and "landed".
3. `characteristic.rst`, sections `characteristic-angular-flux`, `characteristic-wall-coupling` and `characteristic-door`.
4. `orpheus/geometry/chart.py`, `directions_at`.
5. `characteristic/closure.py` and `transport.py`.
6. The hoist's performance target in spec §8.

**Costs to plan with** `[M]` 2026-10-08: the door's gate file runs in about 60 s, outside slow (95 rows; one slow cylinder row of 27 s). The touched trees, outside slow, take 29 min serially, and the CI set about 1 min.

**Rung 5a's working files,** under `scratch/characteristic_architecture/p1_step_b5/`. Run each probe with `PYTHONPATH=. .venv/bin/python -O` from the repository root:
- `surface.md`: the explorer's map of the interface surface (the specification, the observables, `Derivation`, the traced memo, the mesh-free functions, the pieces of the point reading in §6, the insulation gate);
- `spec.md`: the rung's verification spec, with its measured bands and the battery verdicts (including the gauge-ruling section);
- `qa.md` with `qa/p1`–`p8`, and `elegance.md` with `elegance_probe*.py`: the reviews;
- `main/smoke3.py`: the door on the final API (the closed-body gauge, `Nearest`, reciprocity with its 4π, the symbolic detector, the resolution refusals);
- `ta/battery/` (39 arms) and `ta/battery_gauge/` (9 arms): each has a `run.sh` and its verdicts; `ta/probe_n4.py` is the higher-mode gauge measurement;
- `main/touched.sh`: the pre-merge run (the CI set plus the touched trees).

**Lessons from this stretch:**
- **Measure a ruling's premise before relaying it.** "Production's gauge is 100" named the homogeneous solver's density. "SN's gauge" counts the (n,2n) emission. Each was corrected only by a reviewer reading the code.
- **A gauge functional can vanish on a mode the question can return.** Evaluate it on every mode the question can return, not on the fundamental alone. Three agents found this independently.
- **A shared definition between a reference and production is an X4 exposure.** Share the declaration (which channels count). Keep each channel's physics on each side, as the (n,2n) census requires.

## P1 step (b), rung 5b: API sketch (the main agent, 2026-10-08; checkpoint with the user before code)

Every item is a proposal until ruled. The facts it rests on, read at `a5492017`:
- **The pieces exist and nothing composes them.** `TraversalRule.angular_flux(t, inflow)` (`transport.py:557`) already takes several parameters per line, `(..., q)`, and `LinePeriod.inflow(depth, outflow, arriving)` (`closure.py:201`) already takes an injected flux per traversal entry. The diffuse walls' injection (a unit current entering wall w is carried as `1/D_w`, with `D_w = A_w/4`) and the returned currents `j = α(I − Tα)⁻¹Uᵀq` are computed only inside `LineRule.transport` and `WallCoupling.update`, never exposed. `Chart.directions_at` has 0 production consumers.
- **The gates owed are more than the compaction point listed.** The spec's per-point rows C1 (the singular points), C2 (general points), C2b (near a wall) and C4 (convergence of the reading) have no test yet (`[M]` grep of `tests/gates/derivations/test_characteristic_*.py`: the C-numbered rows there are the door's own content rows, not the spec's). The P1 sketch's migration step (b) puts "reading and door (C1-C2, F)" in this rung. So 5b owes C1, C2, C2b, C4, C6, C7's reading leg, C8, E1 and E4.
- **The Galerkin flux is the L2 projection of the transported flux.** `W φ_h = K q` says `(W φ_h)_i = ∫ u_i 𝒦q`, so `φ_h = P(𝒦q)`. Two consequences: the point value of `φ_h` carries the basis's projection error, and the flux integral of a weight `w` by the pairing errs by `∫(w − Pw)(𝒦q − P𝒦q)`, a product of two projection errors.

**1. What `PointValue` reads: the transported emission (iterated Galerkin).** `φ_g(x) = (𝒦_g q_g)(x)`, with `q` the converged emission the pencil or the source pencil returns (the source included, the gauge's scale applied). Recommended: it is the "one transport, two test measures" of the P1 sketch (item 7, ruled), it is what C1 and C2's 1e-12 bars measure, and C6 then gates `∫ u_i φ(x) dV = (K q)_i`. Alternative: the Galerkin flux `φ_h(x)` read from the basis, which costs nothing but carries the projection error, so C1 and C2 could not be met at the 1e-12 bar without a basis that resolves the flux to that bar.

**2. The direction rule at a point: the line rule's own lines, read at the point (one constructor, two test measures).** A line through the point `x` at orbit coordinate `c` with impact parameter `b ≤ c` is congruent to the line domain's canonical line of the same coordinates; `x` lies on it at the two signed positions `±√(c² − b²)` from the closest approach (`RadialImage.parameter_at`). So:
- **The lines are the line rule's.** The point rule's coordinates come from the same constructors as `LineRule.of` (`_impact_rule` for `b`, `_grazing_ends` for the slab's cosine and the cylinder's polar angle), on the partition with the point's `c` inserted as one more end, keeping `b ≤ c`. The grading the point needs is then the grading the line rule already does: the tangencies `b = r_k` are panel ends under the visibility substitution; the grazing direction at the point is the top panel's end `b = c`, where the substitution variable `y = √(c² − b²)` is `c·|μ|` exactly; and a point near a wall is the impact rule's existing grading toward `y = 0` at the next radius out (C2b's boundary layer, regime 3, needs no new law). On the slab the point splits its panel, so the thinnest optical width that grades the cosine toward 0 includes the point's distance to each face.
- **The weights are the point's measure, `dΩ/4π`, in line coordinates:** the sphere `½|dμ|` per branch with `|dμ| = b db / (c√(c² − b²))`; the cylinder `dα dw/π` per branch with `dα = db/√(c² − b²)` and `dw = sin θ dθ`; the slab `½ dμ`. They differ from the line rule's weights, which is the point of "two test measures".
- **The stratum is the one branch.** At the sphere's centre and the cylinder's axis (`DirectionDomain.shape` WHOLE and AXIAL_COSINE) every direction has `b = 0`: the rule is the line `b = 0` read at its closest approach, a `match` on the closed shape enum as `LineRule.of` matches the line shape.
- **One count.** The point rule uses the transport resolution's `line_points`; `Resolution` gains no field.
Recommended over the alternative the P1 sketch named: a separate rule over `DirectionDomain`'s box coordinates, split at `DirectionDomain.tangencies` and substituted there, with its own `direction_points`. That rule would re-derive the tangency, grazing and near-wall gradings in a second variable, which is the two-angular-rules hazard (regime 10, #516) moved one level down.

**3. The point's measure on the emission, and the diffuse walls exposed.** In `reading.py`:
- `PointRule(rule: LineRule, parameter: (L, q))`, built by `PointRule.of(basis, walls, sigma_t, c, points)`: the lines, their weights `dΩ/4π`, and the point's parameters on each line (`q = 2` on the radial charts, 1 on the slab).
- `PointRule.row(transport: GroupTransport) -> (M,)`: the point's linear functional on the group's emission coefficients (design D3's "discrete measure per evaluation point"), so `φ_g(x) = row · q_g`. Per chunk: `TraversalRule.of(lines)`, the specular inflow `period.inflow(depth, outflow, arriving)` with `arriving` the diffuse walls' injection, then `angular_flux(parameter, inflow)` weighted and summed. The diffuse part enters as the returned currents: `row = row_line + r(x) · α(I − Tα)⁻¹Uᵀ`, with `r(x)` the point's reading of a unit current entering each diffuse wall.
- Refactor in the same commit, so that the block and the reading share one definition (X4, Pattern 2): `WallCoupling` gains the walls' breakpoints and their injection `1/D_w`, and a property `currents = α(I − Tα)⁻¹Uᵀ` `(W, M)`, the balance-row solve moved there; `update` becomes `response @ currents`, bitwise or within the solve's rounding (to be measured; a re-baseline is reported, never silent).
- `GalerkinSystem.point_flux(position, emission) -> (G,)`, each group's row applied to its emission coefficients.

**4. The angular flux at a point (public, ruled 2026-10-06).** `GalerkinSystem.angular_flux(position, directions, emission) -> (..., G)`, for directions on `Chart.directions_at(position)`'s box: the line through the point in each direction (`Line.through`), one `TraversalRule`, read at parameter 0. It is the one consumer of `DirectionDomain` here, and its gate is the point reading against a plain fine rule over the box of `angular_flux` (two line constructions, one transport). The derivation exposes it for its gates and mints no observable (the same footing as the adjoint eigenpair; an angular observable is a vocabulary change owed to #529).

**5. The door.** The answers carry the emission beside the flux, scaled together: `_FundamentalAnswer(k, flux, emission)`, `_SourceAnswer(flux, emission)`. `evaluate(PointValue(position, group))` returns `Uncertified(system.point_flux(position, emission)[group])`; on a `_ModeAnswer` it refuses as the flux integral does (`_refuse_mode_flux`). `_refuse_point_reading` retires. A `Ratio` of point values composes as today. `admit_observable` already confines the position to `[r_0, r_n]` and the group to the problem's.

**6. The flux integral stays the pairing.** Its error for a weight outside the basis space is the product of two projection errors (above), so a volume rule over point readings would buy nothing measurable at the cost of a point loop. The alternative is recorded: the integral read as a volume rule over readings, exact up to that rule.

**Risks, named before the build:**
- **A point on a wall.** `angular_flux` refuses a parameter outside its transit by a comparison with the chord's crossing parameter (`transport.py:577-581`). A point exactly on the outer wall, read inward, sits at the transit's entry, computed by a different formula (regime 1; the old family missed by 6.4e-8 to 9.6e-7 there). The build measures whether the comparison needs the crossing's own parameter in place of `parameter_at`'s.
- **Performance** `[R]`: the reading builds no Volterra block (the assembly's dominant cost) but transports per basis function, `L × J × N` per chunk, chunked by the line rule's budget. The spec §8 target is 0.073 s per point (the hoisted brute call, not a min-of-repeats figure). If the protocol misses it, the lever is to contract the emission before the scan (one function transported, not `N`), at the strategy nearest the hot loop (algebra eager, performance lazy).

**Gates** (the test-architect's spec, from the P1 spec §5 rows): C1, C2, C2b, C4, C6 (with its second clause now structural: one constructor), C7's reading leg, C8, E1 (Garcia Case 1 per point, today's tolerances), E4 (a pure absorber against C1's closed forms), plus the angular-flux row of item 4 and the `currents` refactor's row (the block unchanged). Each with its first red.

**Ruled 2026-10-08** (see the ledger): the sketch as written, the four recommendations.

## P1 step (b), rung 5b: built, and two grading defects of rung 3's line rule found by the reading (2026-10-08)

**Built on `feature/characteristic-point-reading`** (uncommitted until the gates and reviews land):
- `reading.py`: `PointRule` (the line domain's lines through the point, weights dΩ/4π, read at `point_parameters`), `stacked_flux`, `on_emission`, `angular_flux`.
- `GalerkinSystem.point_flux`, `angular_flux`, `source_emission`. The door reads `PointValue` from the answer's emission; `CharacteristicDerivation.angular_flux`; `_refuse_point_reading` retired.
- **Shared with the block, one definition each:**
  - `closure.DiffuseWalls`, the diffuse walls' keying and their injection 1/D;
  - `WallCoupling.currents`, with `update = response @ currents`;
  - `assembly.StackedSources`, the stacked emission and wall sources;
  - `traversal_rules`, the chunking;
  - `OpticalScale`, `impact_panels`, `polar_rule` and `cosine_rule`, the line rule's coordinate rules;
  - `TraversalRule.at`, a `LineReading` whose inflow term pairs with any source set.

**Measured** `[M]` 2026-10-08:
- **Smoke** (`p1_step_b5b/main/smoke.py`, `smoke_door.py`):
  - The point measure sums to 1 within 5e-14 on every chart, at the strata and on walls and interfaces.
  - Closed homogeneous bodies (mirror and white) read the infinite-medium flux within 5e-14 at every point.
  - The vacuum absorber sphere matches its closed form within 5e-14 from the centre to the wall.
  - The door's eigen point readings on a closed body match the infinite medium within 1e-13.
- **The refactor** changed no pinned result: 476 passed (assembly, closure, transport, system, outside slow).

**Defect 1: the rim of a partial mirror.** Found by the test-architect from the reading on the wall; it reaches the block too. On a line grazing the outer wall, the closure 1/(1 − a e^{−2Σy/|PΩ|}) has a pole a distance −|PΩ| ln a/(2Σ) from y = 0. The outermost impact panel had no grading toward it.
- Fix: `rim_distance`, which plays the role of the outermost panel's next radius out.
- Sphere wall reading at 12 line points: a = 0.99 went from 9.3e-5 to 1.2e-13. Block 1ᵀK1 at a = 0.99 went from 5.3e-10 to 3.8e-16.
- A hollow body's inner wall has no such pole (the test-architect, measured with the law removed): a line grazing r₀ from outside misses the cavity.

**Defect 2: grading below the impact parameter's resolution.** Found by the test-architect on the cylinder, where the rim law is scaled by the smallest polar speed (about 1e-7). Nodes landed 3 ulp below R; b rounded to R, so the lines made no crossing. 4 of 32 640 lines through the wall point were affected. The block silently dropped them; the reading refused them.
- Fix: every y-grading distance is floored at 2√(2 r ε(r)) over the first Gauss node's fraction.

**Performance:** a point on a cylinder takes 1.5 to 8 s (one and three regions, 8 line points), against the spec §8 target of 0.073 s. Raising the chunk budget gave 1.6 to 1.8x. The remaining lever is to contract the emission before the scan. Not yet measured by the §8 protocol.

**Gate cost:** the 7 characteristic files, outside slow, pass 891 of 892 in 7 min 17 s. The red is the 5a row pinning the point refusal, which the test-architect re-poses.

## P1 step (b), rung 5b: the reviews, two rulings, and the sketch of the widened rung (2026-10-08)

**Correction to the section above.** The rim law did not fix the cylinder's vacuum wall. The smallest polar speed sampled, about 1e-7, drives the law's distance below the resolution floor, so the floor grades the cylinder whatever the law says (the elegance review's probe 3: the rule with the computed rim equals the rule with rim = 0 on 2 of 2 cylinder fixtures). The pole law acts on the sphere: a = 0.99 went from 9.3e-5 to 1.2e-13.

**qa's findings** (`p1_step_b5b/qa.md`), with the point measure itself clean term by term:
- **F1.** The partial mirror's pole sits at every interior tangency whose outer shells are thin or void. A void outer region at a = 0.99 misses by 5.1e-4 at the interface and by 2.6e-6 in the block's 1ᵀK1, at 8 line points.
- **F2.** The cylinder's polar-speed layer exists at interior radii too: an interface point misses by 1.0e-6 at 8 line points.
- **F3.** A line stored by its rounded impact parameter b loses the half-chord y near grazing. The wall error grows as about 3e-15/(1 − a), whatever the resolution: 6.6e-12 at a = 0.999. The floor patches only the extreme case, b rounding onto R.
- **F4.** No catcher outside slow for the layer law.
- **F5.** Inconsistent refusals: outside the body; the slab's grazing ψ returns 0 silently at |μ| ≤ 1e-100 and is off by 1.6e-3 at 1e-15; `c = 1e-160` gives a NaN weight.

**The elegance review's findings** (`p1_step_b5b/elegance.md`):
- **V1.** `PointRule` duplicates `LineRule`, and the guard has already drifted.
- **V2.** `point_parameters` re-spells the chord's crossing expression.
- **C1.** The walls' fold is written twice.
- **C2.** The shared line rules sit in the wrong module: a private import; `_FULL_SOLID_ANGLE` twice.
- **C4.** The floor and the clamp owe debt tags; filed as #590, since #582 is a different defect.
- **C5.** Bare int tuples are passed where `TransportResolution` exists.
- **Nits:** the `diffuse` alias, `DiffuseWalls.of`.

**Rulings, the user, 2026-10-08:**
1. The general grading law lands in 5b.
2. The half-chord is carried in 5b: #590 is done in this rung, not separately.

### The sketch (for the user's ruling of its shape; nothing built)

**K. The kernel: an image carries an exact level** (P0 kernel, `orpheus/geometry/chart.py` and `chord.py`).
- `RadialImage` gains `level` (a radius r*) and `level_half_chord` (y* = √(r*² − b²), exact), both `(...,)`. `Chart.image` sets r* = b and y* = 0, which is today's arithmetic bit for bit.
- One method forms every half-chord: `RadialImage.half_chord_at(r) = √((r − r*)(r + r*) + y*²)`. A radius is crossed where that square is positive.
- One method forms every parameter at a level, on both image classes: `parameters_at(r, side)`, which is ±`half_chord_at` through `parameter_at` on the radial image, and (r − c_foot)/ċ on the axial one.
- `_radial_chord` and `_axial_chord` read these methods, and so does the reading (this is V2: one definition).
- A caller that knows the exact pair passes it: `ConcentricPartition.chord(line, level=None)` with `level = (r*, y*)` per line, attached to the image after the canonical move. Orbit quantities are pose-invariant.
- `Chart.image` is unchanged for every other caller. Retired: the y-floor and `_B_MARGIN` in `_impact_rule`, and the clamp in `point_parameters`.
- The line rule passes each node's (top, half_chord) from `ImpactRule`. A line through a point in a given direction passes (c, c|Ω_x|/|PΩ|), exact from the direction.

**G. The grading law at every impact-panel top** (replaces `rim_distance`). In panel k, toward y = 0, the distance is

d_k = min(s, s·(−ln a) + τ_out(r_{k+1})) / (2Σ_k),

where:
- s is the slowest projected speed the rule samples (1 on the sphere);
- τ_out is the in-plane optical depth of the shells above r_{k+1}, at b = r_{k+1};
- a is the outer wall's specular amplitude. The −ln a term is present for 0 < a < 1 and absent otherwise: no specular return has no pole, and at a = 1 the closure is regular.
- A void panel (Σ_k = 0) has d = ∞.

The outermost panel (τ_out = 0) is today's rim law, and the next radius out still bounds d. `[R]` cost: on the cylinder every panel top is graded to s ≈ 1e-7 times its own mean free path, about 25 more pieces per panel. The cylinder gets slower until #587 (the non-tensor (b, θ) rule) grades per polar angle.

**R. The review fixes, in the same pass:**
- **V1.** One value for a weighted set of lines: `LineSet(basis, walls, sigma_t, coordinates, levels, weights)`. It holds the guard, the sort and `chunks(resolution)`. Its two constructors are the line domain's measure (`LineSet.through_body`) and a point's (`LineSet.through_point(c)`, which keeps c and derives the side). The block and the point reading are two functionals on it. `cosine_rule` returns the mirrored rule, and one function builds the impact-by-polar product.
- **C1.** `WallCoupling.on_emission` is the one fold; `StackedSources` owns the (M + W) layout.
- **C2.** A module `lines.py`, beneath `assembly.py` and `reading.py`, holds the coordinate rules, `LineSet`, `StackedSources` and the one `_FULL_SOLID_ANGLE`.
- **C5.** `TransportResolution` is passed, not tuples.
- **Nits:** the alias retires; `DiffuseWalls.of`.
- **F5.** One refusal message for a point outside the body. The slab's grazing ψ below its parameter's resolution is refused as the sphere's is (#590's family). `c` below the smallest normal float is refused.

**Gates owed** (the test-architect, resumed by name):
- F1's void and thin outer regions, reading and block;
- F2's cylinder interface point;
- F3's wall at a = 0.999 and 0.99999, now inside the 1e-12 bar;
- a fast catcher for the layer;
- the kernel's level rows: `half_chord_at` exact at the level; `Chart.image`'s default bit for bit today's; the crossings at a level equal to the reading's parameters.

**Ruled 2026-10-08:** K as "the image carries a level"; G with the cylinder's cost accepted until #587; R as written.

## P1 step (b), rung 5b landed (2026-10-08)

**Landed:** code `b76b9a9d` (CI red on its own, see the lessons), then docs `bb84bed7` (CI green; the archivist's sections; ERR-104 to ERR-106 with their `catches` markers). Closes #590. Filed #591, the reading's cost. #585 now owns the slab's grazing refusal.

**What landed:**
- `PointValue` reads 𝒦q at the point, and ψ(x, Ω) is exposed.
- `Lines` is the one weighted set of lines. `LineRule` (the block) and `PointRule` (the reading) are its two roles.
- The grading law at every impact-panel top (`lines.tangency_distances`): the next radius, the layer, and the pole.
- The kernel's exact level (`RadialImage.level`, `chart.half_chord`, `parameters_at`, `chord(line, level)`). `at_level`'s agreement check is measured at worst 8 ulp of r*² over 1 192 112 lines (the archivist), with 64 allowed.
- `WallCoupling.currents` and `on_emission` are the one fold of the walls; `update` is retired.

**Evidence** `[M]` 2026-10-08:
- **Gates:** the reading file and `test_chord_level.py` hold 187 rows; outside slow they pass with the door row (173), and slow 15 of 15. The battery has 26 arms, each reddening its targets.
- **Reviews:** qa, two rounds; elegance, three rounds.
- **Pre-merge:** the CI set passes 1427; the touched trees outside slow pass 9978 (251 skipped, 2 xfailed). Sphinx `-W` clean. `dead_references` 0 of 66.

**Cost, the open item (#591).** On the final tree a cylinder point takes 4.9 s (one region, 19 584 lines) and 42.8 s (three regions, vacuum, 44 992 lines), on a loaded host. That is 68 to 590 times the spec §8 target. The grading law at the tensor rule's slowest polar speed multiplied the line count. #587, the non-tensor (b, θ) rule, is the larger lever; contracting the emission before the scan is the second.

**Lessons from this rung:**
- **The reading found four defects in rung 3's line rule:** the rim pole, the interior pole, the layer, and the lost half-chord. Every one was invisible to the block's own gates at the working resolution, or sat inside their bands. A second test measure over the same transport is a strong instrument.
- **A fix whose law is clamped by a floor can be inert.** On the cylinder the rim law did nothing; the floor did the work. Ablate the law by replacing it with 0 and with ∞ before crediting it (the elegance review's probe 3; its lesson L-041).
- **A floor that patches a representation's resolution is a symptom.** The fix was to carry the exact datum, the half-chord, which retired the floor, its margin and a clamp together.
- **Re-measure a cost after every law that changes the rule's size.** The figure quoted in the issue was stale within the rung (the archivist caught it).
- **A code commit that retires a docs target reds CI on its own.** `b76b9a9d` retired `WallCoupling.update` while an `implements` declaration still named it, and the docs commit re-pointed it. CI was red at `b76b9a9d` (the one Sphinx warning) and green at `bb84bed7`. check: before pushing a code commit that retires or renames a symbol, build Sphinx on that commit's own tree, or move the re-pointing docs line into the code commit.

## ⏸ COMPACTION POINT — 2026-10-08, rung 5b merged; P1 step (b) complete

**State.** P1 step (b), rungs 1 to 5b, is merged; read it from git. The reference reads every observable of its three questions.

**Next, per the P1 sketch's migration order (item 9):**
- step (c): the corroboration file (spec §7, G), the old family against the new on the SN fixtures, temporary L4 rows;
- then (d), the SN rows re-pointed and the record keys re-baselined;
- (e), the old family retired with its dependency-audit table;
- (f), the theory page.

**Before (c), measure where the cost lands.** Spec §7's corroboration rows read:
- G1 and G3, the `ABA` sphere and cylinder k: eigenvalues;
- G2, a shape metric: cell averages through the flux-integral pairing, which reads no point;
- G4, the partial reflectors on a slab and a sphere;
- G5, Garcia Case 1: points, on a sphere (cheap: 56 lines a point).

So #591's point cost reaches (c) only through G5, on a sphere. What can reach it is the cylinder's BLOCK: the grading law (`lines.tangency_distances`) grades every impact-panel top at the tensor rule's slowest polar speed in the line rule too, so the block's line count grew with the reading's. That growth was 3.6 to 5.3 times on the point rule (the archivist); on the block it is not yet measured `[R]`.

The first measurement of the next session:
- the `ABA` cylinder solve (the G3 fixture, `scratch/reference_architecture/p3/oracle_hoist/bitexact.py`'s geometry) on the final tree, wall time and line count, against rung 4's cylinder solve;
- then a ruling on whether #587, the non-tensor (b, θ) rule, precedes (c).

**Read in order, for step (c):**
1. The P1 sketch, item 9 (the migration order).
2. `scratch/characteristic_architecture/p1_verification_spec.md` §7 (the 13 SN rows, the corroboration file's G rows, the AST and spy preconditions, the positive controls) and §9 (X4: the new reference's independence from the SN solver it judges).
3. `characteristic.rst`, sections `characteristic-grading-law`, `characteristic-reading` and `characteristic-door`.
4. P0's L4 corroboration pattern, the model for the file's delete-this-row rows (`retirement-audit` D.14).
5. #587 and #591.

**Costs to plan with** `[M]` 2026-10-08:
- **Gate files:** `test_characteristic_reading.py` with `test_chord_level.py` and the door row, outside slow: 173 rows in about 155 s. Their 15 slow rows take about 26 min.
- **Characteristic and geometry together**, outside slow: 2562 rows in about 10 min.
- **Pre-merge:** the touched trees outside slow (derivations, geometry, numerics, specification, reference, homogeneous, data, and the SN question rows) take 30 min, beyond the 30-minute default background limit, so run them with the 2-hour limit. The CI set takes about 1 min.
- **Sphinx `-W`:** about 30 min.
- **One cylinder point:** 5 to 43 s (#591).

**Open issues from this campaign:** #584, #585, #587, #588, #589, #591.

**Working files:** `scratch/characteristic_architecture/p1_step_b5b/`. It holds:
- `spec.md` and the gates' design;
- `qa.md` with `qa/` and `qa/r2/`;
- `elegance.md` with its probes;
- `ta/battery/` (`run.sh <arm>`, `verdicts.md`);
- the measurement probes under `main/`: `smoke*.py`, `kernel_bitwise.py`, `level_tolerance.py`, `budget_sweep.py`, `touched.sh`;
- `archivist_point_cost*.py`.

**Two process checks this rung paid for:**
- **Build Sphinx on a code commit's own tree before pushing it,** when it retires a docs target (`b76b9a9d` was red on its own).
- **Strip `catches` markers from a code commit whose ERR entries land in the docs commit.** Copy the files aside, remove the marker lines, commit, restore, and check with `diff -q`. Markers and entries then land together, and the reconciliation test stays green at each commit.

## #587 before step (c): the measurement and the ruling (2026-10-08)

**The measurement** `[M]` 2026-10-08, `scratch/characteristic_architecture/p1_step_c/aba_cyl_cost.py` (logs `cost_5a.log`, `cost_head.log`). The `ABA` fixture (`aba_specification`, two groups, reflective) at resolution (3, 2, 0.4, (8, 12, 12), 8), the k eigenvalue through `CharacteristicDerivation`, one process each, the traced memo bypassed. "Rung 5a" is the tree at `fe977a90`, before the grading law; "current" is `11b263a5`.

| body | rung 5a | current |
|---|---|---|
| sphere | 0.7 s; 240 and 256 lines per group block | 0.7 s; 240 and 280 lines |
| cylinder | 384 s; 38 400 and 38 912 lines; k = 1.231729452890078 | 2296 s; 308 480 and 297 920 lines; k = 1.231729452890264 |

The cylinder's lines grew 8.0 times and its time 6.0 times, and k moved by 1.5e-13 relative. On this reflective body (amplitude 1, so no pole) the extra lines come from the layer term s/(2Σ_k), graded at the tensor rule's slowest polar speed at every polar angle. Rung 5a's k lies 5.7e-4 relative from the old family's 1.231036749830859, inside G3's 1.9e-3.

**The ruling** (the user, 2026-10-08): #587 precedes step (c). It is a W3 carve, with an API sketch before any code.

**The probe of the per-angle rule** `[M]` 2026-10-08 (`p1_step_c/iterated_probe.py`, `iterated.log`; `per_theta_counts.py`). The probe patches `LineRule.of` in memory so that each polar node θ_j gets its own impact rule, graded at sin θ_j. On the same fixture it gives:
- 102 320 and 106 784 lines;
- 1275 s;
- k = 1.231729452890254, 8e-15 from today.

The cost by polar angle, for group 0's 160 polar nodes: 51 sit at sin θ < 0.01 and carry 57 536 of the 102 320 lines. The polar rule's own grazing grading, down to τ_min/64 of the whole body, puts them there. That is #587's other half (grade θ per impact panel), which is not measured.

**The sketch, ruled 2026-10-08 ("Approve per-angle").**
1. `tangency_distances` splits in two:
   - `Tangencies.of(chart, ends, sigma, amplitude)` holds, per panel top, the next radius's half-chord, τ_out and the pole's −ln a, from one chord of the tangent lines;
   - `Tangencies.distances(speed)` gives the distance per panel at that speed.
2. `impact_by_polar` retires. `impact_per_polar(impact_at, polar, impact_weights, polar_weights)` builds the iterated rule.
3. `LineRule.of` and `PointRule.of` (with `_impact_below`) take it, and `slowest` retires. The sphere and the slab are unchanged bit for bit.
4. Docs:
   - in `characteristic.rst` (the grading law, the cost sections, the list of what is not built) and in ERR-105, the law's s becomes each line's own speed;
   - the test-architect writes the gates (the structural per-angle gate, red on the tensor rule);
   - C9, T_w, conservation and the ERR-104/105 catchers stay unchanged.

After it lands, the polar half is measured by a probe and comes back as its own question.

Branch `refactor/cylinder-iterated-rule`.

## #587's first half landed (2026-10-09)

`[LANDED 200b4233]` (code) and `1004a666` (docs). CI `gates` is red on `200b4233` alone and green on `1004a666`. The cylinder's line rule is iterated: each polar node's impact rule is graded at its own sin θ. The rule has two parts:
- `lines.ImpactPanels`, with its methods `of`, `distances(speed)` and `rule(points, speed, below)`;
- `lines.impact_per_polar(impact_at, polar)`, which returns plain db·dθ weights; each role then applies its own density.

The ABA cylinder holds 102 320 and 106 784 lines per block, against 308 480 under the tensor rule, for the same k to 8e-15. The sphere and the slab are bit-identical, 86 of 86 rules.

**Evidence:**
- qa: no reading moved by more than 1e-15. Grading at speed 1 moves the interface reading by 1.3e-7. (`p1_step_c/qa.md`)
- The test-architect's battery (`p1_step_c/ta/battery/`, ported to the final names): the five per-angle gate rows are the only catchers of the "slowest" mutation. The ERR-104 and ERR-105 arms redden their catchers.
- Characteristic and geometry outside slow: 2572 passed.
- The CI set plus the traced memo: 1540 passed.
- Sphinx `-W` clean; `dead_references` 0 of 66.

The archivist re-measured the point reading on a loaded host (load 6 to 10):

| cylinder | lines before | lines now | time before | time now |
|---|---|---|---|---|
| one region | 19 584 | 8544 | 4.94 s | 3.08 s |
| three regions | 44 992 | 15 704 | 42.8 s | 23.4 s |

**Found in review and fixed before the commit:**
- The slab's `LineRule.of` raised, because the impact panels were built before the match on the chart (the test-architect).
- The elegance review's three concerns:
  - the panels are now held with their features, so the builder is not twinned;
  - the role's measure is applied after the iterated rule, so it is no longer welded into it;
  - the `zip` now has `strict=True`.
- qa's F2: a speed of 0 gave NaN; it now gives 0.

**A process check this carve paid for:** when a code commit adds or removes gate rows, the regenerated `docs/theory/verification/matrix.rst` belongs in the code commit, not the docs commit. Otherwise CI's "the build changed no tracked file" step reddens the code commit on its own. This is the second rung in a row with a red code commit. Both times the red was a generated or derived docs artefact riding in the docs commit.

**Open in #587:**
- the polar rule graded per impact panel. 51 of the 160 polar nodes sit at sin θ < 0.01 and carry 57 536 of the 102 320 lines. Not measured.
- #587's done-when, the three-region cylinder block at 16 points under about 3 minutes, is not measured.

**Next:** measure the polar half with a probe and bring it back as its own question (the ruling of 2026-10-08). Then P1 step (c). The ABA cylinder k solve now costs about 21 minutes, from the single probe run.

## The polar half refuted; step (c) next (2026-10-09)

The polar half of #587 (grade the polar angle per impact panel) is refuted for the question "does it cut the ABA cylinder's lines". `[M]` 2026-10-09, `p1_step_c/polar_half_scales.py`:
- The panel basis's thinnest panel is the outer-wall panel [1.96, 2.0], 0.02 and 0.04 optical in groups 0 and 1.
- Every line crosses it, so each impact panel's polar rule keeps the body's 160 or 152 nodes.
- The fact it establishes: the polar grazing grading is set by the outer wall's panel.

The remaining idea, a polar weight absorbing sin²θ, is unmeasured and stays in #587.

**Ruled (the user, 2026-10-09):** proceed to P1 step (c). Accept about 21 minutes per ABA cylinder solve, from the single probe run of 1275 s.

## Step (c): the context and four rulings (2026-10-09)

**Context.** The explorer's map is `p1_step_c/explorer_c.md`, with its probes in `p1_step_c/explorer/`.

P0's pattern is `tests/gates/geometry/test_kernel_corroboration.py`. It has three parts:
- an AST precondition, `_assert_independent`;
- a runtime spy, the `kernel_calls` fixture with `_without_kernel`;
- two positive controls.

It carries no markers.

Independence today: 0 of 20 old modules import the new package, and 0 of 12 new modules import the old.

The old door serves G1–G3 only, through `aba_specification`. G4 and G5 call the old array solvers directly.

Preliminary readings, not gates:

| row | measured | spec tolerance | status |
|---|---|---|---|
| G1 | 8.33e-5 | 3.22e-4 | within |
| G2 | 4.05e-5 of M | 1.37e-3 | within |
| G3 | 5.63e-4 | 1.9e-3 | within |
| G4 slab | 3.28e-5 | 2.5e-5 | red |
| G4 sphere | 4.13e-5 | 1e-4 | within |
| G5 interior | ≤ 1.07e-3 | 1.43e-3 | within |
| G5 at r = R | 2.26e-2 | 1.65e-2 | red |

Both reds are the old reference's error:
- The slab ladder converges at second order. Its limit matches the new value to 4e-7, and its finest rung errs by about 3.2e-5 by Richardson.
- At the surface the old value is 2.3e-2 from Garcia, while the new is 1.2e-5 from Garcia.

**Rulings (the user, 2026-10-09):**
1. **G4 slab and G5 surface:** bounded by the old reference's MEASURED error, rounded up to one significant figure (`ceil_one_significant_figure`).
   - The slab bound is Richardson on its ladder: 3.2e-5, giving 4e-5.
   - The surface bound is the old value's distance from Garcia: 2.3e-2, giving 3e-2.
   - Step (d) must re-point the SN Garcia rows, whose tolerance rests on the understated step.
2. **G4 posed with `isotropic_mixture("A")`,** its arrays asserted `array_equal` to the old solver's inputs. Step (d) checks that the SN partial-reflector rows run at scattering order 0.
3. **`slow` marks:** G1–G4 are slow; G5 runs in the default tier.
4. **The file's shape, approved:**
   - no level marker;
   - an AST precondition over the old family's 20 modules;
   - a runtime spy on the new entry points;
   - the old side under `traced_memo.bypass()`;
   - two positive controls, one of them a memoised stand-in that calls the new entry point;
   - P0's `_assert_independent` and spy moved to a shared `tests/gates/_corroboration.py`, used by both files;
   - rows named as declarative sentences.

**Noted for later steps:**
- Spec §9's claim that the reference reads albedos through `specular_albedos` is stale.
- G2's tolerance input (`_trajectory_resolvent_ladders.py:46-49`) was measured on the retired nodal reading.
- Step (e)'s blast radius includes P0's corroboration rows that import `chord_oracle` (their lines 188, 221, 253).

## Step (c) landed (2026-10-09)

`[LANDED 17708ac8]`. CI `gates` is green on the commit by itself: the regenerated `matrix.rst` rode in it, as the process check of #587's landing requires.

**What landed:**
- `tests/gates/derivations/test_characteristic_reference_corroboration.py`: 10 rows, with no level marker. G1–G4 are slow; G5 and the controls run in the default tier.
- The shared `tests/gates/_corroboration.py`. P0's `test_kernel_corroboration.py` now uses it.

**The readings** `[M]` 2026-10-09, the whole file plus P0's: 17 passed in 2667 s.

| row | reading | tolerance |
|---|---|---|
| G1, ABA sphere k | 8.3e-5 | 3.2e-4 |
| G2, ABA sphere shape | 4.0e-5 | 1.4e-3 |
| G3, ABA cylinder k | 5.6e-4 | 1.9e-3 |
| G4 slab | 3.28e-5 | 4e-5 (Richardson at the observed order 1.9465 gives 3.35e-5) |
| G4 sphere | 4.1e-5 | 1e-4 |
| G5 interior | ≤ 1.07e-3 | 1.43e-3 |
| G5 at r = R | 2.26e-2 | 3e-2 (the old value's distance from Garcia) |

So statement (i) of spec §7 holds: the new reference sits inside the old one's own error at every G row. Statement (ii), D9's error estimate of the new reference at each SN row's fixture against the old estimate, belongs to step (d).

**Independence:**
- Static legs: the old family's 20 modules, its three test helpers, and every member of their import closure (153 first-party modules) name nothing of the new package.
- Runtime leg: the old side executes no code object of the new package. The leg counts through `sys.monitoring` under `traced_memo.bypass()`.
- The legs' declared blind spots are in the helper's docstring.

**Review:**
- qa (`p1_step_c/qa_c/`): no old-side code runs outside `bypass()`, and the ceilings redden. Four findings, all fixed before the commit.
- The test-architect's reds: `p1_step_c/ta_c/`.

**Next: step (d).** Re-point the 13 SN rows. Re-baseline `CYLINDER_3REG_RECORD`'s `k_ref` and `phase_c_k_gap`; its three SN keys must not move. Check statement (ii) per row.

Carried into (d):
- The SN Garcia rows' tolerance rests on the understated step (1.65e-2 at the surface, against a measured old error of 2.3e-2).
- Check that the SN partial-reflector rows run at scattering order 0.
- Spec §9's `specular_albedos` claim is stale.
- G2's tolerance input was measured on the retired nodal reading.

Into (e): P0's corroboration rows that import `chord_oracle` are in the blast radius.

## Step (d): the context and three rulings (2026-10-09)

**Context** (`p1_step_d/explorer_d.md`):
- **13 SN consumers** of the old reference family, 17 test ids. 11 of them go through `_aba_reference.aba_reference` / `verify_cylinder_k` (the late import at `_aba_reference.py:118`). The two partial-reflector rows and the parametrized edge row call the old solvers directly.
- **A census** over 1132 tracked files found no 14th SN consumer.
- **No SN Garcia rows exist.** The step (c) carry-over about them is void.
- **Scattering order 0 confirmed.** A spy on the partial-reflector rows' own calls saw `scattering_order=0`.
- **D9 was never built,** and no function returns the new reference's error.
- **The record:** `assert_record` uses a 2e-5 band and has no fingerprint. Under the new reference, `k_ref` moves 5.6e-4, and `phase_c_k_gap` goes from 5.97e-4 to 3.46e-5.

**Hazards:**
- `test_crosscheck_harness.py:103` checks the awaiting rows by name and by "#516" in their reason.
- 8 `rests_on` ids point at old-family tests.
- The labelled equation `sn-curvilinear-trajectory-resolvent-crosscheck` names the old reference and τ_k = 4e-3 in its body.
- The ABA sphere k row tightens from 4e-3 to 4e-5 and reads 1.8e-5.

**Rulings (the user, 2026-10-09):**
1. **Statement (ii)'s error comes from a ladder module,** `tests/gates/derivations/_characteristic_ladders.py`. Per SN fixture it holds the new reference's step to a rung above its working point, as `[M]` literals with a re-measure CLI, consumed by `tolerance_for`. The cylinder's rung above the door default is measured once, offline.
2. **The partial-reflector rows keep their own sum rule,** recomputed from the new reference's error. It is tighter and never looser: the sphere row's band becomes 7e-4 against a reading of 6.54e-4, and its docstring records the margin.
3. **Step (d)'s shape, approved:**
   - `aba_reference` and `verify_cylinder_k` move to `characteristic_reference`: the sphere and slab at p = 5, the cylinder at the door default;
   - the partial-reflector and edge rows move to the new door;
   - `k_ref` and `phase_c_k_gap` are re-baselined, and the three SN keys must not move;
   - the 7 cylinder xfails stay strict (#566);
   - no row is renamed and the #516 reason is not reworded;
   - the 8 `rests_on` ids are re-pointed to the new reference's gates;
   - the docs commit rewrites the labelled equation;
   - the regenerated matrix goes in the code commit.

## Step (d) landed (2026-10-10)

`[LANDED d9425977]` (code; CI green on its own) and `94bfbf51` (docs; CI green).

**The re-point.** The 13 SN reference rows (17 ids) read the characteristic reference. The sphere and slab run at p = 5, the cylinder at `Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)`.

**The new estimates** live in `tests/gates/derivations/_characteristic_ladders.py`, `[M]` tables with an argparse CLI:
- on geometric ladders, the step divided by (1 − r), with r the largest ratio measured on any ladder (0.2957);
- stride 2 for the ABA sphere shape;
- the cylinder's rung (4, 3) with the transport rule held (3399 s), plus the transport step.

**Statement (ii) holds on every row,** and every tolerance tightened or held. The largest tightening: the ABA sphere k, 4e-3 → 4e-5.

**The record:** `k_ref` = 1.2317294528902538, `phase_c_k_gap` = 3.4606e-5. The three SN keys did not move.

**The restructure:**
- `_ladder_rules.py` holds the pure tolerance rules.
- `_trajectory_resolvent_aba.py` holds the old A|B|A spelling and the ERR-094 old arrays, so step (c)'s rows keep comparing the same old spellings.
- `_trajectory_resolvent_ladders.py` now holds only the old reference's tables, so step (e) deletes it whole.

**qa (`p1_step_d/qa_d/`):**
- F1: both sides read the boundary law's response factor, `SpecularReturn.kernel`. It is part of the posed problem, and its external pin is `test_characteristic_walls.py`. Spec §9 is corrected.
- F2: ERR-094's corner moves SN's k toward the reference at this fixture, so the band cannot see it. The leakage re-drop reddens both rows.

**Still `[R]`:**
- the claim that the ABA sphere k row is now a catcher of `tau := 0.7`;
- the cylinder shape estimate's transport term, the k step standing in for a shape step nobody measured.

**Next: step (e).** Retire the old family in one commit, with its dependency-audit table (`retirement-audit` G.24). Its blast radius includes:
- the old family's 20 modules and its own tests (the 77 RETIRE rows of the spec's roster);
- step (c)'s corroboration file, which is deleted;
- the `chord_oracle` imports in `test_kernel_corroboration.py` (lines 156, 193 and 233);
- the Garcia tables imported by `test_characteristic_reading.py:60`;
- the traced-memo client tests;
- the import inside a subprocess string at `test_ambient_state.py:141`;
- the `cross_method` tests;
- `_trajectory_resolvent_ladders.py` and `_trajectory_resolvent_aba.py`.

The explorer's census is in `p1_step_d/explorer_d.md` §2 and in its memory file `characteristic_reference_sn_consumers.md`.

## Step (e): the audit and four rulings (2026-10-10)

**The audit table** is `.claude/plans/characteristic_reference_architecture_dependency_audit.md` (the explorer; probes in `p1_step_e/explorer/`).

**Its main finding:** step (e) cannot be one deletion commit.
- **Free to go:** 12 numeric modules (8 054 lines), with 0 production callers outside the family. No docs directive turns red on its own.
- **Coverage that would be lost:**
  - 123 of the spec's KEEP or RE-POSE rows sit in the 22 test files deleted whole. Their successors (D8, D9, C9, C12, C6b, and the moved `verifies` markers) were never built.
  - ERR-034, ERR-035 and ERR-091 would lose their only catchers.
  - 11 `peierls-greens-*` labels would lose their verifier.
- **What must survive:** the SymPy `origins/` (8 modules). 18 of its 23 `derive_*` functions are kept. 5 retire with the 26 duplicate rows; no KEEP row reaches them.

**Rulings (the user, 2026-10-10):**
1. **Migrate, then delete.**
   - (e1): the successors built on the new reference while the old family still exists, each shown to pass beside its predecessor. The KEEP rows are re-pointed; D8, D9, C9, C12 and C6b are built; the three ERR entries get catchers; the `verifies` markers move.
   - (e2): the deletion commit, driven by the audit table.
2. **`origins/` moves to `characteristic/origins/`.**
3. **Both outside consumers are re-pointed:**
   - `cross_method`'s adapters move to the characteristic reference;
   - the traced-memo rows M4.1b–M4.8 are re-posed on the new reference's memoised `evaluate`.
4. **The small losses:**
   - the backward-segments check is replaced, before the delete, by a hand-computed multi-region segment table (closed-form chord lengths);
   - uncalled helpers (`import_closure`, `assert_closure_independent`) are deleted in (e2);
   - the `_ladder_rules` functions stay while `_characteristic_ladders` uses them.

**Order:**
- (e1a) the `origins/` move: main agent, a production move.
- (e1b) the successors: the test-architect.
- (e2) the delete.

## Step (e1b) landed; #592 ruled before (e2) (2026-10-10)

**(e1b) [LANDED 44303919]**, with the error-index regeneration [LANDED 2f56e69b]. CI was red on 44303919 alone, because the generated `.claude/skills/vv-principles/error_index.md` drifted when the new `catches` markers raised four entries' catcher counts. The lesson: a code commit that changes gate rows carries every file a Sphinx build regenerates, the matrix AND the error index. Run the build before the commit and read `git status`.
- Every one of the 123 rows has a successor. The record, with per-row outcomes and the qa round, is `scratch/characteristic_architecture/p1_step_e/ta_e1b/README.md`; §9 covers the qa round.
- The default tier keeps cylinder eigenvalue coverage. The closed cylinder and the closed annulus read k_inf at rung 2, and Ua-1-0-CY reads k = 1 within 3e-5. Both rows go red by value under a line speed scaled by 1.0001.
- The old default-tier rows whose successors are slow only are listed in §9, each with its reason: 3, 13–15, 22, 27, 65, 113.

**Ruling (the user, 2026-10-10): #592 lands before (e2).** The characteristic solve becomes a traced-memo child entry, so the child-entry feature never has 0 production clients, and the reference meets the reference-cache plan's unit of caching before the family that met it is deleted.

**The API (the user approved it over a module-level `_solve(specification, resolution)`):**
- `@traced_memo def solve(self) -> _Answer` on `CharacteristicDerivation`; it carries the old `answer` body.
- `answer` stays a cached_property, `return self.solve()`.
- `evaluate` and `angular_flux` are unchanged.
- The key is the derivation's content.

**Smoke [M] 2026-10-10** (the tiny sphere, isolated cache root):
- one ChildPin in the reading's manifest;
- read-only arrays;
- the bypass reading is bit-equal.

The gates come from the test-architect (`p1_step_e/ta_592/`). The reviews are qa (`qa_592/`) and the elegance-enforcer (`ee_592/`).

**The reviews of #592 (2026-10-10):**
- **The elegance-enforcer:** no violation.
  - The `solve`/`answer` pair is two concerns, not a weld.
  - Its nit ("once per instance") is applied.
  - Its concern, the theory page crediting `answer` (an `implements` target), went to the archivist.
- **qa:**
  - The payload is exact for all 4 question kinds.
  - The keys are sound, and the refusal order holds.
  - The `answer`-skips-`solve` mutant is blind to the tracked gates (8 of 8 green) and caught by the new M4.3c, M4.11, M4.12 and M4.13 (4 of 4 red).
  - 4 reads of `.answer` in `test_characteristic_reference.py` now reach the real store; they are wrapped in `bypass()`.
- **Both reviews:** a bypassed answer's arrays are writable, while a served one's are read-only. This is filed as #593 (a fix at the definition site, `bypass()` round-tripping the payload).

**A pre-existing defect, fixed first (the user's ruling, 2026-10-10: sort).**
- The traced memo's exact key wrote a frozenset in iteration order. That order follows the caller's hash seed, while generations run under `PYTHONHASHSEED=0`.
- So one call keyed differently across processes `[M]` (qa: 5 unseeded processes split 3 and 2 on an eigen gauge with an (n,2n) material).
- `_Exactly.elements` now sorts by encoding, and the generator contract gains the clause.
- It lands as its own fix commit before #592. Its gate (keys under several seeds) is the test-architect's (`ta_592/README.md` §2).

## #592 landed; (e2) next (2026-10-10)

**Landed:**
- The set-key fix: [LANDED 45405f86].
- #592: [LANDED 23e7711e] (`Closes #592`).
- The docs: [LANDED 875e1054].

CI was green on the code push (run 38085759201). The seed gate's real-case leg reads the reading key (`evaluate`), not `solve`'s, so the fix commit is green without #592. The solve's parent-against-reading split is gated by M4.11b, in the #592 commit.

**Open, filed:** #593. A bypassed answer's arrays are writable, while a served answer's are read-only; the fix is at the definition site.

**Next: (e2), the deletion,** in the order of the audit table (`characteristic_reference_architecture_dependency_audit.md`): leaves first, one commit (`retirement-audit` G.24). It also carries:
- the archivist's catalogue text for ERR-034, ERR-035, ERR-090 and ERR-091, and a possible ERR entry for the old method-of-images k defect (5.3e-4, which both sides shared);
- `dead_references` after the delete;
- the uncalled helpers (`import_closure`, `assert_closure_independent`).

The child-memo feature keeps a production client: the characteristic solve.

## Step (e2) landed: the trajectory-resolvent family is retired (2026-10-10)

**[LANDED 5aa8cb88]**, pushed with CI watched; `Closes #572`.

**What went:**
- the package, 12 modules;
- the 5 duplicate `derive_*` (18 kept);
- 25 test files.

**What changed beside the deletion:**
- 41 test files edited.
- The three `peierls_nystrom` primitive rows were MOVED to `test_peierls_specular_primitives.py`. Their successor map had re-posed them onto characteristic rows that test another object (retirement-audit C.12).
- 7 SN ids and the label `sn-curvilinear-characteristic-reference-crosscheck` were renamed.
- `[M]` `-m "not slow"`: 16 613 passed, 0 failed. `-W` is clean.

**Findings:**
- `test_docstring_xrefs` reds on dead `:mod:` roles even on a page awaiting its rewrite. The 14 roles on `trajectory_resolvent.rst` became literals in this commit.
- qa found two present-tense citations of deleted rows on `readings_and_certificates.rst`; both are fixed.
- qa's design finding, `HollowBody` and `LayeredBody` left with no production reader, is filed as #594: it may be a second classification beside `ConcentricPartition`.
- `dead_references` finds 55 sites, all on `trajectory_resolvent.rst`.

**Next: step (f), the theory page.** The archivist rewrites `docs/theory/references/trajectory_resolvent.rst`, 6 438 lines:
- its 55 dead references;
- the labels that the moved rows still verify (`peierls-greens-slab-T`, `peierls-greens-cylinder-T`);
- what becomes history, against what moves into `characteristic.rst`.

That is P1's last step.

## P1 complete: step (f) landed (2026-10-10)

**[LANDED 54f6c6e7].** `references/characteristic_origins.rst` replaces `trajectory_resolvent.rst`, by the user's ruling: split.
- The page holds 37 live labels:
  - 17 proved symbolically, in `characteristic/origins/`;
  - 17 verified numerically, by the characteristic tests;
  - 3 both.

  An AST census of the markers equals the page's label set.
- The history section records the retired family and its lessons.
- The 31 unverified labels are dropped.
- 43 inbound references are re-pointed.
- `dead_references`: 55 → 0.

**qa's claim check found:**
- **The `peierls-greens-slab-T` hypothesis** (a source symmetric about the mid-plane). It was missing from the label, and its moved marker sat on a row with a non-symmetric source. The fix:
  - the hypothesis is stated;
  - a closed-form row, `test_a_symmetric_slabs_surface_inflow_is_the_one_transit_closed_form`, goes red under an asymmetric source and under ERR-035's heuristic.
- **The T₀₀ rows credited the closure labels.** They now verify their own identity labels, with a new `peierls-greens-cylinder-V-alpha-2`.
- **`derive_alpha_zero_kernel_reduction_slab` wrote ERR-035's refuted closure.** It now writes the one-transit form, and its limit holds either way.

**P1, the record:** steps (a) to (f) are all landed, the last being `54f6c6e7`. The characteristic reference is the one reference. The old family is deleted (`5aa8cb88`), and its solve is a memoised child entry (`23e7711e`).

**Open issues from P1, in GitHub:**
- #585: a slab far from the origin loses digits;
- #587: the non-tensor line rule; its first half landed and its polar half was refuted;
- #589: the k pencil and the source pencil as two splittings;
- #591: the point reading's cost;
- #593: bypass versus served arrays;
- #594: `reference_body` against `ConcentricPartition`.

## ⏸ COMPACTION POINT — 2026-10-10, P1 complete

**Next, by the plan's phases (a `[HYPOTHESIS]` until the user rules):**
- **P2 (#578):** the other spellings onto the geometric kernel: 13 chord square roots, 14 discriminants, 8 point locators, 3 line measures; reference code first.
- **P3:** #405 resumes: the slow reference files, re-timed, through `.claude/plans/reference_cache.md`'s P4.

**The user chose P2 (2026-10-10).**

### P2 resume brief (for the session after compaction)

**State.**
- `main` is `6ac11067`; CI is green on `54f6c6e7`.
- P1 is complete.
- No P2 work has started: there is no branch and no uncommitted change.

**What P2 is.** Issue #578 (read it first: `gh issue view 578`). The hand-written chord, point-location and line-measure spellings migrate onto the geometric kernel in `orpheus/geometry/`:
- `Chart`, `Line` and `ConcentricPartition.chord`;
- `region_containing`, the inner-owns locator;
- `Chart.beam_density`.

The kernel's theory page is `docs/theory/foundations/chart_and_chord.rst`. Each migration:
- deletes its row in `tests/gates/geometry/test_kernel_corroboration.py`;
- adopts the inner-owns convention;
- runs `retirement-audit`.

**What #578's census no longer says right.**
- The census was taken at `a336bde4` (2026-10-05): 13 chord square roots, 14 discriminants, 8 locators and 3 line measures. P1 deleted the trajectory-resolvent family, and with it the 7 chord-oracle classes and their spellings, so the counts are stale.
- `[M]` 2026-10-10, the corroboration harness today:
  - one exit-distance row parametrised over `peierls_rho_max` (`peierls_nystrom.geometry.SPHERE_1D`) and `moc_root` (`moc.geometry._ray_circle_intersections`);
  - one locator row over Peierls `SPHERE_1D`, `mc.solver.ConcentricPinCell` and `moc.geometry._identify_region`.
- So what remains sits mostly in:
  - Peierls Nyström (`derivations/continuous/peierls_nystrom/`), WITHDRAWN under #506; its `CurvilinearGeometry.kind` carve is #419;
  - production CP, MoC and MC.

**The open scope question, the FIRST thing to settle with the user, before any plan.**
- The sharpening order (CLAUDE.md, "The direction of development") defers CP, MoC and MC to their own campaigns, and admits only "a harmonisation onto shared machinery".
- #578 says reference code goes first, and that Peierls follows the scope ruling on #506 and #419.
- So P2 needs a ruling on which spellings are in scope now:
  - (a) reference code only (Peierls, withdrawn but still in the tree; and any `derivations/` site outside Peierls);
  - (b) also the production CP, MoC and MC sites, as a harmonisation onto shared machinery;
  - (c) reference code now, and production through issues filed per method campaign.

**The first concrete step,** before asking: a fresh census by the explorer, against `main` `6ac11067`. It needs:
- #578's predicates: a chord square root `sqrt((R*R|R**2) - b...)`, and a geometric discriminant outside `orpheus/geometry/` and the SymPy origins;
- point locators, by their boundary convention;
- measures on lines;
- an AST census with positive controls, each site tagged by package and by status (reference or production; withdrawn or live; which method);
- counts of the form `k of N`.

Then put the scope question to the user with the counts, then plan mode (the ontology is known: the kernel exists, and P2 re-spells sites onto it).

**Standing constraints** (unchanged from P1; this session's practice):
- W3 surgical carve: the main agent writes production; the test-architect handles gates; qa and the elegance-enforcer review; the archivist writes docs.
- Stage by path, never `git add -A`. Create the branch in its own command (the git-guard refuses otherwise).
- Merge with `--ff-only`, then delete the branch.
- Commit messages through `-F -` with a quoted heredoc.
- A code commit that changes gate rows also carries every file the Sphinx build regenerates: `docs/theory/verification/matrix.rst`, and `.claude/skills/vv-principles/error_index.md`, whose omission turned CI red once this session. Run the build and read `git status` before committing.
- Before committing a deletion, run `tests/gates/test_docstring_xrefs.py`: it rejects dead `:mod:` roles that Sphinx does not see.
- Run tests serially, with `.venv/bin/python -O -m pytest -p no:cacheprovider`. The full non-slow suite takes about 50 minutes `[M]` (16 613 passed at `5aa8cb88`): run it in the background without `-x`.
- `git grep`'s basic regex has no `\|`: use `-E` or Python for a census, with a positive control (it bit this session).
- `rm -rf` is refused by the permission settings: ask the user to run it.

Read this point first on resume.

## P2 — scope rulings and the plan (2026-10-10)

**The census** (`[M]` the explorer, tree `a7aeebb5`; `scratch/characteristic_architecture/p2/census.md`, AST plus text filters, 7 positive controls). Outside the kernel and the SymPy origins:
- **Chord sites, 28:** 22 reference-live (common 4, characteristic 1, fn_method 1, Peierls 16), 2 with only withdrawn consumers, 1 dead, and 3 production (all MoC).
- **Point locators, 14:** production 6 (MoC 1, MC 3, mesh 1, numerics 1); reference-live 7; withdrawn 1.
- **Line-measure sites, 11:** 5 density spellings and 6 quadrature rules.

**The premise correction.** #506 withdraws only the Peierls *solver* half. `CurvilinearGeometry` and the `compute_P_ss_*`, `compute_T_specular_*`, `P_esc*` and `G_bc*` primitives are reference-live.

**The cost** (`[M]`; `scratch/characteristic_architecture/p2/caller_shapes.md`).
- One kernel chord over one line costs 229 µs. Over a batch it costs 0.33–6 µs per line.
- The Peierls scalar primitives cost 0.3–3.8 µs per call, and are called two loops deep.
- So the Peierls layer migrates batched per observer: about +0.16 s in its tests, against +11.7 s for a scalar re-spelling.

**The user's rulings (2026-10-10, AskUserQuestion):**
1. The reference sites move now. Production gets one issue per method campaign (MoC, MC, CP).
2. Generic interval locators (`IndicatorBasis.evaluate`, `mesh/factories.pwr_pin_2d`) are out of scope.
3. The Peierls ray layer is batched per observer: one batched chord per observer (the plan first said one `ConcentricPartition` per `compute_*` call; `rays_from` mints one per observer, 28 µs `[M]`, immaterial beside the chord, so the wording was amended on 2026-10-11). This dissolves the slab/curvilinear branches of the ray methods into the `Chart`, a slice of #419.
4. `chord_half_lengths` is re-spelled on the kernel even though production CP imports it. This is a harmonisation, and CP's numbers move toward the cancellation-free value.

**The plan** (approved; `~/.claude/plans/zesty-dazzling-hanrahan.md`). Branch `refactor/kernel-p2-reference`. The commits, each green on what it touches:
0. the test-architect's specification;
1. K: `impact_rule` onto `half_chord`. This is 0 ulp `[M]`; the 126 memo entries regenerate once, with identical values;
2. F: the Fn sphere reconstruction, batched;
3. C: `chord_half_lengths`;
4. P: the Peierls ray layer, possibly with a `Chord.lengths_between` verb;
5. D: the densities onto `beam_density`, and the deletion of the dead `_chord_tau_mu_sphere`;
6. the issues (MoC, MC, CP), comments on #419 and #506, and the closing of #578;
7. the docs.

**Kept, with the reason:**
- the mpmath-precision chord bodies, because the kernel is float64;
- the quadrature rules, because the kernel owns the density and not the rule;
- the test oracles, because re-pointing them would be demotion.

**CP defects found in P2 (2026-10-10, W2 numerics-investigator; memo `scratch/characteristic_architecture/p2/cp_pi_investigation/memo.md`):**
1. **The sphere's line density.** CP writes it as 2y where it is 2πy, so the sphere's P_cell is π too small, in production and in `flat_source_cp` alike. `[M]`:
   - two independent grounds (the closed form, and a 3-D point escape integral) agree to 4e-14, and both codes read exactly 1/π of them;
   - under white, k moves by +6.4e-5 and the flux shape by up to 13.6 %.
2. **CP's vacuum k omits leakage** in every geometry, the CP twin of ERR-064. `[M]`: slab 1.23 against a true 0.52; sphere 0.59 against 0.018.

**The user's ruling (2026-10-10).** Both are symptoms of CP lacking the shared machinery: a measure on lines, a boundary trace, an operator pencil.
- **Defect 1:** group D SEEDS the measure. CP and the reference read the line density from `Chart.beam_density`, which fixes the π at its root, in a form the CP campaign will build on.
- **No CP plan is opened now.** CP's development is a known latent move; P2 keeps marching, so that CP's campaign starts on a good verification infrastructure.
- **Defect 2:** filed as an issue, with an ERR entry and a strict-xfail gate, until CP's pencil exists.

## P2 step 7, the docs: written to the working tree (the archivist, 2026-10-11)

Not committed; the main agent commits it with the code. Every commit cell below is *(in development)* until then.

**The error catalogue** (`docs/theory/verification/error_catalog.rst`; it ended at ERR-106):
- **ERR-107**: the CP sphere's density `2y` where the measure is `2πy`, in production and in `flat_source_cp`; the root of #283 and #509. Catchers named: `test_cp_first_collision.py::test_the_first_collision_probability_is_the_closed_form` and the large-radius sphere cells of `test_white_boundary_infinite_medium.py`.
- **ERR-108**: the Peierls cylinder-3D response paired the exit distance along +cos ω with τ along −cos ω. `[M]` against qa's 3-D ray trace: the old pairing 1.2 % to 11.0 % off; the one-ray reading within 2.3e-5 at 16 points. Catcher: `test_peierls_ray.py::test_the_cylinder_3d_responses_are_their_characteristic_integrals`.
- **ERR-109**: the `len(radii) == 1` shortcuts ignored `inner_radius`. `[M]` 13 (function, geometry, mode) rows over 7 functions moved 1.4 % to 22 %. Catcher named: `test_peierls_ray.py::test_a_hollow_single_region_body_has_a_void_cavity_in_each_consumer`.
- **ERR-110**: CP's vacuum k omits leakage (#546). No catcher by design: the strict xfail documents it, and its catcher lands with #546's fix. Once the three markers are placed, the index reads 2 uncaught: ERR-103 (its marker sits on one parametrize row, which the graph does not see) and ERR-110.

**The pages:**
- `chart_and_chord.rst`: a new section, `chart-and-chord-rays`, on the rays from a point and the five verbs, with their derivations, the measured guards they retired, and the AST census of consumers (13 modules, 1 of them production). Also: the Key facts, the half-lines, point location, the deferred census, the gates, three refuted designs, three gotchas, and a history row.
- `collision_probability.rst`: the #509 warning, the sphere's weight, and `second-diff-sph` / `self-sph` with `2πy dy`. A new anchor, `cp-sphere-measure-on-lines`, carries the measurement. The cylinder's factor 2 is now explained as S₀, not as the two halves of the chord.
- `peierls_nystrom.rst`: a new section, `peierls-ray-layer`. `cp-unified-outer-integration` and its measure table now carry `2πy`. The retired names are re-tensed, and the slab refusal now holds at every n.
- `peierls.rst`: the ray walker and its short-circuit.
- `characteristic.rst`: the refusal message, `impact_rule`'s inverse map, and `line_integral`.
- `fn_method.rst`: the label `fn-sphere-psi-characteristic` and its code paragraph. It has no `verifies` marker yet; F1 should carry it.

**The hollow responses (the archivist, 2026-10-11, second dispatch; W2 memo `scratch/characteristic_architecture/p2/hollow_investigation/memo.md`).**
- **ERR-111**: `compute_G_bc_outer` counted the rays that meet the cavity first (G/P = 5.20).
- **ERR-112**: the cylinder's `compute_G_bc_inner` used the far-side `Ki_1/d` form.
- **ERR-113**: `W_oo` was −inf for every Σ_t > 0.
- **ERR-114**: the rank-1 cylinder `compute_G_bc` was the surface-centred `Ki_1/d` form. It is now `4·compute_P_esc`. #112's rank-N half (`compute_G_bc_mode`, n ≥ 1) stays open.
- The catchers are the test-architect's rows; their names are reconciled at commit.
- `peierls_nystrom.rst` changes:
  - the F.4 "does NOT close" claim is refuted for the fixed-source identity; k = k∞ is reasoned (`[R]`), not measured;
  - §22.8 lists two surface-centred consumers, not four;
  - the observer rule includes a breakpoint at the observer;
  - `peierls-specular-T-matrix` lost its `documented` status, because its row now verifies it.
- `[M]` After the last `-W` build, the index read 113 entries; after the hollow-response entries it reads 114 entries, 373 catching tests, 2 uncaught: ERR-103 and ERR-110 (`[M]` the `-W` build of 2026-10-11).

## P2 landed (2026-10-11)

**The commits** on `refactor/kernel-p2-reference`: the code `124edf67`, the docs `ba407665`, and this note with the agent memories. The code commit closes #283 and #509.

**The pre-merge evidence:**
- **The slow characteristic rows.** `[M]` 60 passed, 1154 deselected, in 3 h 31 min (`-m slow` over the 15 characteristic and Fn cross-verification files and `test_composition_cost.py`; log `scratch/characteristic_architecture/p2/characteristic_slow.log`).
- **What those rows check, and what they do not.** They call the reference under `traced_memo.bypass()`, so they compute it in-process on the final tree and compare it with their pins. That is the direct check of the pinned values against the re-spelled kernel. They read and write no cache entry. The earlier claim that they regenerate the cache was wrong.
- **The cache is stale and will rebuild lazily.** `[M]` 103 of 103 `evaluate` entries and 24 of 24 `solve` entries in `.cache/references` read Stale. The `solve` entries name `orpheus/geometry/chart.py` (the `half_chord` signature widening), and the `evaluate` entries name their stale `solve` children. The audit script is `scratch/characteristic_architecture/p2/audit_memo.py`. It must run with `PYTHONPATH` set to the repository root, because a pinned path resolves through `sys.path`: run from the scratchpad, it first read every pinned function as "gone". The cache is gitignored, so CI never reads it. The first local test that reads an entry regenerates it.
- **The full non-slow suite** passed before the last fixes (`scratch/characteristic_architecture/p2/full_nonslow_2.log`). The files that failed in it were re-run alone and passed. The last edit changed only an annotation (`half_chord` accepts scalars). The suite was not re-run after it.
- **Sphinx `-E -W`:** exit 0. The incremental `-W` build after the history hashes: exit 0.
- **`dead_references`:** 0 dead of 65 targets, on a graph rebuilt from this tree. Positive control: the graph holds the new `rays_from` and none of the retired Peierls methods.
- **pyright:** 0 new errors. The 12 errors left in the changed files are on lines this branch did not change.
- **Docstring xrefs and the catalogue gate:** 54 passed.
- **The retirement audit:** `.claude/plans/p2_peierls_dependency_audit.md`. 0 of the 18 non-archaeology hits is a present-tense claim or an imperative.

**The issues:**
- **Filed in P2:** #595 (MoC), #596 (MC) and #597 (CP).
- **Commented with the hashes:** #546 (ERR-110 open), #112 (rank 1 fixed, rank N open), #419 (the ray-layer slice done) and #506 (what stays withdrawn).
- **Re-scoped and closed:** #578.

## ⏸ COMPACTION POINT — 2026-10-11, P2 landed

**Next, by the plan's phases** (a `[HYPOTHESIS]` until the user rules): **P3**, where #405 resumes. The slow reference files are re-timed through `.claude/plans/reference_cache.md`'s P4, against the speed of the parked hoist `a336bde4`.

**On resume:**
- Read `git log --oneline -5` and `gh run list --branch main --limit 3`. Do not trust this note's merge state; check it with `git merge-base --is-ancestor 124edf67 main`.
- The characteristic cache is cold. The first non-bypassed run regenerates its 127 entries, which takes hours. Time P3 against a warm cache, or say that it was cold.

## P3 opened (2026-10-11): #405 resumes — living section

**The ruling.** The user, 2026-10-11: "Once everything looks good, let's start P3." P2 was already on `main` (`124edf67`, `ba407665`, `4ebfa7a4`; CI green at `ba407665`), so P3 starts now.

**What P3 owes, read from the plans** (`[M]` read 2026-10-11):
- this plan's phase line: "the slow reference files, re-timed against `a336bde4`'s measured speed";
- the P1 verification spec §8 (`scratch/characteristic_architecture/p1_verification_spec.md`): the protocol for one point reading, the cylinder `ABA` at r = 1.2, old (the hoisted brute call at `a336bde4`) against new, at matched accuracy, min over ≥ 15 alternating subprocesses. It was never run: rung 5b recorded 1.5 to 8 s per cylinder point against the 0.073 s figure, "not yet measured by the §8 protocol";
- `test_runtime_405.md`, "Step 3 moved": step 3 of #405 is done when this phase re-times the slow reference files;
- `reference_cache.md`'s P3 timing protocol (bypassed, cold, warm twice, per file, serially, in a detached worktree; driver `scratch/reference_architecture/p3/timing/run.sh`).

**The file set changed since that protocol ran.** `[M]` 2026-10-11: of its 11 files, the 5 of the deleted families are gone (`test_trajectory_resolvent_reference.py` and the four `test_peierls_greens_function_*` non-symbolic files); 6 remain (`test_aba_specification.py`, `test_crosscheck_harness.py`, `test_phase_c_crosscheck.py`, `test_l1_standoff_slab_cylinder.py`, `test_unified_matvec_cylinder.py`, `test_reference_body.py`). The memo has one client now, `characteristic/reference.py`. The new set is those 6 plus the characteristic consumers: the 11 `test_characteristic_*` files, the two Fn `*_xverif` files, `test_singular_eigenfunction_cylinder_xverif.py` and `test_partial_reflector_resolvent.py`.

**The steps** (`[HYPOTHESIS]` until the measurements are in; each decision after a measurement goes to the user):
1. **The §8 point-reading protocol** (a numerics-investigator; the host idle otherwise). Output: old and new minimum and median, both errors against the truth, the activation counts.
2. **The file re-timing** (bypassed, cold, warm twice) on the new file set at `main`, after step 1, serially in a detached worktree. Output: the table beside P3's of 2026-10-04, and the hottest files, which order `reference_cache.md`'s P4 ("the families, hottest first").
3. **If step 1 misses the target:** the lever named at rung 5b, contracting the emission before the scan (one function transported instead of N), is a code change to the reading. Its design goes to the user before any code.
4. **Close-out:** the parked hoist `a336bde4` retires (its oracles were deleted in P1, so it cannot merge); the branch deletion goes to the user. #405's own remaining items (R2, R3, #504) stay in `test_runtime_405.md`.

**Step 1 measured: the new reading misses the target by about 19×.** `[M]` 2026-10-11. The artefacts are in `scratch/characteristic_architecture/p3/step1/`; the driver is `driver.py c_2_1_2_4_4 15` and its log is `driver_c_2_1_2_4_4.log`. The reading is φ_g(1.2) of the cylinder ABA, both groups counted as one call. Each side ran in 15 alternating subprocesses; each subprocess made one untimed warm-up call and 5 timed calls. Settings: `OMP_NUM_THREADS=1`, Accelerate BLAS, `-O`. The host was not idle (load average 9.6 at the start, 6.5 at the end).

| side | rule | min of the minima | median of the minima | max φ_g error, vs the panel p = 4 truth | φ₁/φ₀ error |
|---|---|---|---|---|---|
| old, `a336bde4` | 48 μ × 256 azimuths, n_traj 64 | 0.122 s | 0.124 s | 1.19e-3 | −3.1e-5 |
| new, `main` | `Resolution(2, 1, 0.4, TransportResolution(2, 4, 4), 6)` | 2.29 s | 2.30 s | 6.8e-4 | −1.3e-3 |

- **The new rule.** It is the coarsest of the 10 rules probed that matches the old accuracy. The coarsest rule probed (one line point) reads in 0.54 s at an error of 5e-2, so no rule choice meets the target.
- **The old error is mostly one gauge factor.** Both groups are off by +1.2e-3 while the ratio is off by 3.1e-5: the old solve ran at the gate quadrature, where k = 1.1703. Under a ratio criterion the miss widens.
- **The controls.** Activation counts were above 0 on every timed call: 512 segment builds per old call and 29 `_stacked_flux` calls per new call. Every reading was bit-equal to the reading of its own setup.
- **Where a new reading spends its time** (`profile_c_2_1_2_4_4.txt`, 2.49 s in total):
  - 55 %: `StackedSources.of` and `outflow`, the full-chord exit and entry integrals, one per basis column;
  - 25 %: the chunking, the chord and its pieces;
  - 12 %: `PointRule.of`, of which 0.18 s is `leggauss` recomputed 2401 times;
  - 8 %: `rule.at`.
- **Contracting the emission before the scan** `[R]`, from the tottime arithmetic, not prototyped. It saves about 0.5 s and leaves about 15× the old time. Each line carries p + 1 = 3 columns, not N, so the per-basis factor is small.
- **The structural lever** `[M]` (`line_census.log`). The new rule chords about 2280 lines per group, with 214 distinct impact parameters over 38 polar angles, and each line is chorded separately. The old side chords 256 in-plane lines and lifts each over its 48 polar cosines almost for free. The impact rule is graded per polar angle (#587), which prevents one in-plane chord from serving every polar angle.
- **The full solve** (reported, not gated): old 7.01 s min (N = 15), new 11.64 s min (N = 5). The two k accuracies differ: the old k is off by 5.0e-2, the new by 3.8e-4.

Step 3's premise is therefore revised: the named lever (contracting the emission) cannot close the gap alone. The candidates go to the user.

**The user's rulings on step 1** (2026-10-11, asked after the measurement):
- Step 3 is the **axial lift**: each in-plane line is chorded once and lifted over every polar angle, as the hoist did. The impact rule is then graded in-plane, not per polar angle, so #587's grading is re-opened. The design goes to the user before any code.
- Step 2 (the file re-timing) runs **after step 3**, so that it times the reading that ships.

**Step 3 measured: the axial lift meets the target in a scratch prototype.** `[M]` 2026-10-11, numerics-investigator. The files are in `scratch/characteristic_architecture/p3/step3/`; no production file changed. The host was heavily loaded (load average 14 to 126), but the old minimum reproduced step 1's to 0.2 %.
- **Shareable fraction** (`q1_split.py`). Of production's 2.30 s reading, 65 % (1.50 s) does not depend on θ once each in-plane line is chorded once: `basis.values` 26 %, the pieces 20.5 %, `PointRule.of` 7 %, the nodes 7 %, the slots, chord and period 4.5 %. 34 % is per (b, θ): the exponentials 22 %, the carried scan 6 %, the scatter 5 %. Of the 56 piece builds, 27 are discarded when the budget halves a chunk.
- **The lift identity.** A 3-D line at polar angle θ is its in-plane line with every optical depth scaled by κ = 1/sin θ and the reading multiplied by κ. On production's own `TraversalRule` the identity holds to 1.2e-7 relative on the matched rule and 1.1e-12 on the panel-4 rule. It is not bit-equal.
- **The prototype** (`lifted.py`):
  - Each in-plane line is chorded once, and θ is a broadcast axis over the exponentials. The specular closure is a·B/(1 − a·e^{−κT}).
  - Gradings: the impact rule is production's, at the slowest polar node. Its layers are nested, so this is the union, and its cost grows with log(1/sin θ_min). Along the line, each slot is graded dyadically from the slowest node's mean free path, with a floor of 100 that moves the result by under 1e-6.
  - Controls: summed over production's own lines, the kernel reproduces production's reading to 3.7e-11. Activation counters were above 0.
- **Timing, the §8 protocol** (`driver_lifted_matched.log`, 15 alternating subprocesses per side). The lifted reading: min 0.0612 s, median 0.0630 s, errors 7.2e-4 and −6.2e-4, ratio error −1.3e-3. The old reading: min 0.1222 s, median 0.1244 s. With 8 points throughout, the lifted reading reproduces the truth to 2e-16 in 14 s, where production takes 355 s.
- **The analytic θ integral (Bickley–Naylor Ki₁)** is accurate (8e-11 at 16 points) but 4 to 6× slower. It needs a Neumann series for the specular closure, which diverges as the cycle product Π → 1. That is the regime the pole grading guards, so it is refuted as the general form.
- **The Galerkin block** `[R]`: the same lift folds the θ-sum of the Volterra triangle into one kernel per in-plane line. The expected gain is up to the number of polar nodes on the shareable 65 %.
- **Not measured:** the gauge-free ratio criterion. It is limited by the solve rung (it needs panel 2 or panel 3, solves of 375 s or 1259 s), not by the reading.

**The user's rulings on the lift's design** (2026-10-11, asked after the measurement):
- **Scope: both the reading and the block.** One lifted-line object serves the point reading and the Galerkin block, so no per-(b, θ) path survives on the cylinder.
- **The polar integral is a numerical polar rule,** broadcast over the shared chord. The Bickley–Naylor form is refuted for this use: its closure series diverges as Π → 1.
- **The along-line grading floor is derived** from the vanishing depth (`grading.VANISHING_DEPTH`), not tuned, and measured on the closed-form rows.

**Step 3's carve** (`[HYPOTHESIS]` until the test-architect's spec is in; branch `refactor/characteristic-axial-lift`):
- **3a. The specification** (test-architect), written to `scratch/characteristic_architecture/p3/step3/gate_spec.md`. It covers:
  - the lift identity as a gate: a 3-D line equals κ times its in-plane line with Σ replaced by κΣ. Its witness is production's per-(b, θ) path, kept as a test oracle until the carve deletes it;
  - the closed-form rows (C9, T_w) and conservation, unchanged within their bands;
  - the re-baseline list: the reading and block values that move, each with its reason under the three-condition acceptance;
  - a mutation that reddens each gate (X1);
  - the memo regeneration: every cylinder entry goes cold once.
- **3b. The lifted line rule.**
  - The in-plane impact rule is graded at the slowest polar node (`ImpactPanels.rule` at sin θ_min; the layers nest, so this is their union). This re-opens #587, whose per-polar grading retires.
  - Along each line, one set of nodes is graded from the slowest speed, with a floor derived from the vanishing depth.
  - θ is a broadcast axis of the exponentials and of the closure a·B/(1 − a·e^{−κT}).
- **3c. Both consumers onto 3b:** `PointRule` (the reading) and `LineRule` (the block, whose Volterra θ-sum folds into one kernel per in-plane line). The cylinder's `impact_per_polar` retires; the sphere and the slab are unchanged.
- **3d.** Fix the 27 of 56 piece builds that `Lines.chunks` throws away when it halves a chunk.
- **3e. Review and docs.** qa and the elegance-enforcer in parallel; the archivist updates `docs/theory/references/characteristic.rst` ("The cylinder's cost on a thick body") and the development history.
- **3f. Re-time.** The §8 protocol on the landed reading, then step 2.

**Step 3a landed: the gate specification** (`[M]` 2026-10-11, the test-architect; `scratch/characteristic_architecture/p3/step3/gate_spec.md`, probes in `.../step3/ta/`). Three findings change the carve:
- **The prototype's closure is wrong between the corners.** `lifted.py` spells a·B/(1 − e^{−κT}), not a·B/(1 − a·e^{−κT}). It is exact at a = 0 and a = 1, and every step-3 measurement used the A|B|A mirror (a = 1), so the numbers above stand. At a = 0.9 a homogeneous cylinder's wall reading is off by 1.1e-2, and at a = 0.5 by 2.0e-2; more points do not reduce it. This is ERR-035's shape. The carve reuses `LinePeriod.inflow`, which is correct, and never re-spells the closure.
- **No fast row sees a cylinder closure defect today.** Under `-m "not slow"`, 16 of 1141 rows of the `test_characteristic_*` files reach the cylinder's rule, and none checks 0 < a < 1 against a reference. The spec adds fast rows: K1c (one line against `MpLine`) and F1, F1b (wall readings against frozen mpmath values).
- **The derived grading law.** The first graded depth is the slowest polar node's own in-plane mean free path, sin θ_min/Σ. A first depth 4× larger cannot be told apart at 8 points; 16× costs 6.7e-11 on a thick wall, and 100× costs up to 1.3e-6. The reach is K = ⌈log₂(VANISHING_DEPTH · sin θ_max/sin θ_min)⌉ + 1 doublings, which reaches 64 mean free paths of the FASTEST node. Stopping at production's 7 doublings misses every closed form, by 7.8e-9 to 2.7e-3. At a single speed the law reduces to production's 7 doublings exactly.
- **#587's retirement costs no accuracy.** Per polar node, the impact rule graded at the slowest speed was worse than the own-speed rule on 0 of 384 nodes. Grading at speed 1 is off by 2e-5 to 8e-4 per node, while the total moves only by 1e-8, so the new gate P1 reads each node.
- **The existing gates.** Simulating the impact half of the carve on today's tree reddens exactly 5 non-slow rows: LR2 (2) and Q4b (3). Both are inverted by the carve and are re-posed in its commit. The 7 slow closed-form cylinder rows tried stay green. The transport rows on 3-D cylinder lines keep gating `TraversalRule` and `angular_flux`, but no longer the cylinder's reading or block, so their scope sentences owe a re-scope.
- **The memo.** Every cylinder entry goes stale. Every solve entry also regenerates if 3b–3d edit a definition it pins (`LineRule.of`, `Lines.chunks`, `StackedSources.of`, `TraversalRule`).

**The user's rulings on the lift's form** (2026-10-11, asked after the spec):
- **The lift is a speed on the traversal rule.** A cylinder 3-D line is its in-plane line traversed at speed v = ds/dt = 1/sin θ, where s is 3-D arc length and t is in-plane arc length. The rule gains one field, the speeds, which broadcast against the line batch; a leading axis of the speeds is a family over the same lines. Every optical depth is Σ·v·Δt, and every length weight along the line is v·Δt; the speed enters there and nowhere else. The pieces, the nodes, the basis values and the period are built once per in-plane line. The closure (`LinePeriod.inflow`) is reused unchanged. The exponential grading becomes graded over the family's range of optical speeds, by the derived law above. The sphere and the slab run at speed 1, and their values must be bit-identical (a gate). `Lines` holds the in-plane lines and the polar rule; its weights are the tensor product of the impact weights and the polar weights times the 3-D density, so the spec's (b, θ) gates re-pose by an accessor. A separate lifted-line class was rejected: it would re-spell the closure, the pieces and the attenuated integral for one geometry, a twin path, and the closure re-spelling is where the prototype's error lived.
- **`angular_flux` stays on the 3-D chord** at speed 1, as the "reading and block" scope ruled. The per-line transport rows against `MpLine` keep covering a production path, and the 3-D chord against the in-plane chord with a speed is an independent check of the lift (spec gate K3).
- **Noted, not in scope:** the slab's line at cosine μ is its normal line at speed 1/|μ|, the same structure. The ruling keeps the slab unchanged; an issue records the option.

**Step 3b–3d written to the working tree; the reading misses the target** (`[M]` 2026-10-11; probes in `scratch/characteristic_architecture/p3/step3/impl/` and `.../baseline/`; the pre-carve tree is the worktree `scratch/characteristic_architecture/p3/step3/wt_head` at `4ebfa7a4`).
- **What is written.** `TraversalRule.speed`, a family axis through the pieces, the attenuated integral, the scan and every functional. `exponential_ends` takes the family's steepest and gentlest attenuation. `Lines.speeds`, with weights `(V, L)`. `axial_lift` replaces `impact_per_polar`. `LinePeriod.inflow` reads its rank from the optical depth. The chunker plans its halving from per-line piece counts (3d). The block folds the family before the basis contraction, as batched matmuls. The packing's scatter is a segmented `reduceat`.
- **Correctness.**
  - The lift identity on single lines: the in-plane line at speed 1/sin θ against the 3-D line at speed 1 agree to ≤ 4.3e-13 (depth, outflow, entry response) and ≤ 2.5e-11 (Volterra) for sin θ from 1 to 0.01.
  - Sphere and slab: 5 bodies × 4 arrays against `4ebfa7a4`. Every reading row and boundary response is bit-identical. Every Volterra block is within 2.9e-16 of its largest entry, a reordered sum.
  - The cylinder block at 4 points sums to within 1.3e-9 of the old block. Its largest entry difference (6 % on one diagonal entry) is the old rule's own error at 4 points: the old block's entry moves from 0.014650 (4 points) to 0.013807 (8 points), and the new block reads 0.013755 at 4 points.
- **Cost.**
  - The cylinder block (MR3, a = 0.6, 4 points, `profile_block*.log`): old 73.6 s; first lifted version 218 s; after the fold, the matmuls and the segmented scatter, 53.7 s.
  - The point reading at step 1's matched rule (`reading_quick.py c_2_1_2_4_4`): 2.57 s, against production's 2.29 s and the target 0.122 s. Its errors are the prototype's: 7.2e-4 and −5.9e-4.
- **Why the reading does not gain** (`reading_census.py`, `reading_profile.py`).
  - Per group the lifted reading carries 126 in-plane lines × 38 speeds = 4788 line-speed pairs, against 2280 3-D lines before. The impact rule graded at the slowest node is the finer rule.
  - Each in-plane line carries up to 342 pieces, because the family's grading reaches 64 mean free paths of the fastest node from both ends of every slot.
  - Production's traversal rule is two-level: every piece's inner rule is graded toward its own end, and a scan carries every basis function along the line. The Volterra block needs that, because it needs the flux at every node. The point reading needs only three targets (the point on each side, and the exit) and one emission. The prototype integrated straight to each target, on one node set graded toward it, with the emission contracted first. In the reading's profile the attenuated integral takes 1.24 s of 2.75 s and the basis-wide scan 0.62 s.

**The user's ruling on the reading** (2026-10-11, asked after the measurement above): **a target-graded reading.** The traversal rule gains a second evaluation, for the flux at a few targets: integrate from the line's nodes straight to each target (the point, the exit), on nodes graded toward it, as the prototype did. The row stays a functional of the emission coefficients, with no scan over every basis function. The block keeps the two-level path, which it needs for the flux at every node. The pointwise reading `at(t)` and `angular_flux` move onto the target-graded evaluation, so that one quadrature serves point values. Rejected for now: landing the carve with the reading unchanged (P3's target stays unmet), and a separate investigation first.

**Step 3a's gate rows are written** (the test-architect, 2026-10-11): the new file `tests/gates/derivations/test_characteristic_lift.py` (K1a–K1d, K2, K3, K4, G1, G1b, F1/F1b/F4, G2, the slow twin, P1, S1), plus the re-posed LR1, LR2', Q4b', M1, M2, M4 and M4.3c, and the transport file's re-scope. 40 rows pass against the working tree as it stood at their run; G2, P1 and S1 have not run yet. S1 compares SHA-256 hashes against a capture from `4ebfa7a4`, so the reading change will move it (see below).

**The target-graded reading, written** (`[M]` 2026-10-11; probes in `scratch/characteristic_architecture/p3/step3/impl/`).
- **A chord through a refinement.** `TraversalRule.on_chord` accepts a chord through any refinement of the basis's panels with the same ends; each region maps to the panel it lies in. `reading.refined_at(basis, c)` inserts the point's level. The reading's lines and `angular_flux`'s rays are chorded through it, so the point is a crossing on each side. The diffuse walls re-key onto the refinement (`DiffuseWalls.on`).
- **Graded piece sets.** A rule's pieces are three sets (`_PieceSet`, by `_Toward`):
  - toward the exits, for the outflow and the flux at a crossing (`at`);
  - toward the entries, for the entry response;
  - toward both ends, for the Volterra integral alone.

  Every slot is cut at its middle. A one-sided set is graded across the whole slot. Each set is built only when a functional reads it, and a solid body's period reads no traversal backward.
- **One-level point values.** The outflow, the entry response and the flux at a crossing are one-level integrals on each piece's own nodes. The flux at a crossing is a sum over the upstream pieces, with no scan over the basis. The block's scan keeps the two-level integral (`_to_end`), because it reads the flux at piece ends inside a slot.
- **Accuracy** (`exits_accuracy.py`, a three-region cylinder's outflow at 8 points against 24): the in-plane family reads 6.6e-15 or better and the single 3-D line 4.3e-13. Two intermediate gradings were refuted, each for a stated reason:
  - grading the whole slot toward its exit alone, with no middle cut: a thin slot stays one piece, 7.9e-11;
  - grading only the half nearer the exit: the far half stays one piece, 1.2e-13 on the family.
- **Sphere and slab against `4ebfa7a4`** (the capture at 4 line points, 6 along the line): reading rows move by 3.3e-11 to 5.7e-11 of their largest entry, a quadrature change; blocks by at most 3.0e-13; boundary responses by at most 3.4e-12. These are re-baselines under the three-condition acceptance, not bit-identity. The test-architect re-poses S1.
- **Cost** (`reading_quick.py c_2_1_2_4_4`, min of 5 after a warm-up, not yet the §8 protocol): 0.39 s, against 2.29 s before the carve and the target 0.122 s. The errors are unchanged: 7.2e-4, −5.9e-4, ratio −1.3e-3. The steps:

  | form | time |
  |---|---|
  | lifted, two-level | 2.57 s |
  | point values one-level, scan-free | 1.16 s |
  | per-panel gather by matmul | 0.94 s |
  | point rule budget 2^18 | 0.72 s |
  | exit-graded pieces | 0.37 to 0.39 s |

  The cost now scales with the pieces, about 300 per in-plane line graded both ways, not with the nodes: 2 points along the line read in 0.64 s, 4 in 0.72 s, before the piece sets. The prototype's 0.061 s used 2 points along the line and a first graded depth of 100 slowest mean free paths, a floor the test-architect refuted at the closed-form bars.
- **The one-level sets are graded at ratio √2** (`grading.ONE_LEVEL_RATIO`; `exponential_ends` takes a `ratio`, default 2). The fast tier (`impl/runs/fast1.log`: 56 failed, 1145 passed) showed P1 off by up to 9.8e-13 per polar node, against its 1e-13 bar. That is the one-level integral on dyadic pieces: over [0, 64] mean free paths at 8 points, ratio 2 misses by 7.4e-13 (7 pieces) and ratio √2 by 2.2e-16 (13 pieces). The Volterra set keeps ratio 2, because its inner rule grades each piece. After the change the in-plane family's outflow reads ≤ 7.8e-15 at 8 points, and the reading costs 0.63 s, up from 0.39 s.
- The other 55 reds: 46 rows read at a parameter inside a slot, which is now refused (the transport file's rows and the old K1b/K1c/K2, which the test-architect is re-posing); the rest are S1's hashes and two `Lines(...)` calls that the run collected before the test edits.
- **The gap to the target, reasoned** `[R]`. Production's grading law is fixed by `VANISHING_DEPTH` and resolves every attenuation to double precision at any resolution, so a coarse resolution still pays for machine-precision grading. The prototype at the derived floor would cost about 1.8× its 0.061 s (its own exponential count), and about 2× more again at ratio √2, so roughly 0.2 s. Production's remaining factor of about 3 over that is implementation: padding, the both-ends set built to size the batch, and the gathers onto the basis columns.
- **The reading after caching the piece integrals:** 0.51 s (`reading_quick.py c_2_1_2_4_4`), with the errors unchanged. The characteristic fast tier is green (`impl/runs/fast2.log`: 1202 passed, 3 skipped, 16.8 min).

**The user's ruling on the reading's cost** (2026-10-12): **accept and continue.** The resolution-independent grading law stands (the "derive it" ruling). P3's reading target of 0.122 s is recorded as **not met**, for this structural reason: production's grading resolves every attenuation to double precision at any resolution. The polar family stretches the grading from the slowest node's first mean free path to 64 of the fastest node's, and the one-level reads need ratio √2. The hoist reached its matched accuracy with a far coarser grading.

Rejected for now:
- coupling the grading to the resolution, which would change the ruled law;
- measuring the prototype at an identical grading first.

Next: 3e (qa and the elegance-enforcer, then the docs), 3f (the slow files re-timed with the §8 protocol), then step 2.

**Step 3e: the reviews, and the fixes applied** (2026-10-12). qa ran 23 mutation arms (`scratch/characteristic_architecture/p3/step3/qa/`); the elegance review read the code (`.../step3/elegance/`). Fixed in the working tree, by the main agent:
- **The speed is a scalar or one family axis `(F,)` shared by every line** (elegance B1, the blocker). Before, a `(V, 1)` speed was a family and a `(V,)` speed on V lines was per-line, told apart by shape alone, so a dropped `[:, None]` was admitted whenever the speed count equalled the chunk's line count. A speed of two or more axes is now refused.
- **The refinement is the kernel's.** `ConcentricPartition.refined_at(c)` and `ConcentricPartition.coarsening(coarse)` (`orpheus/geometry/chord.py`) replace `reading.refined_at`, the inline refinement check in `TraversalRule.__post_init__`, and three spellings of "the coarse region holding each fine region" (`PanelBasis.region_of_panel`, the transport rule, the reading). The wall re-key law is one function, `walls.end_on`, used by `Walls.on` and `DiffuseWalls.on` (elegance S2).
- **`ImpactPanels.lifted_rule(points, polar)`** holds the "grade at sin θ_min" law that the block's and the reading's cylinder arms had each spelled (elegance S3). `ImpactPanels`' parameter is renamed `projected_speed`: "speed" named both v = 1/sin θ (the rule's) and sin θ (the impact rule's), and a swap type-checks (qa finding 4).
- **A batch is sized from the exits-graded set**, which every consumer builds, so the reading no longer builds the both-ends set just to size its chunks (elegance C1). **The integrand v·w·e^{−Σvd} is one function**, `_attenuated_length` (elegance C3).
- **Prose corrected** (qa finding 6, elegance S1 and C4): the transport module's "one body", `_attenuated`'s interval count, the claim that ratio 2 holds for the two-level rule, the cost proxy that is a no-op on the cylinder, and `ImpactPanels`' "no less accurate on 384 of 384", now stated with its predicate (no worse than 1.5 times the node's own rule plus 1e-14).
- The fast gates the fixes touch are green before the ruling below (`fix/smoke1.log`: 1033 passed, 45 deselected, 11.4 min).
- **Not fixed, for a stated reason:**
  - **C2, `at` taking a crossing (a level and a side) instead of a parameter**, is not done in this step: it re-spells the reading's and the tests' calling convention. Today a parameter that is not a crossing is refused at the door.
  - **qa finding 5, the budget not bounding the cylinder**: a rule of one in-plane line holds its whole family (7616 to 24 752 piece slots against the default budget 1024; 52 to 164 MiB per rule at 8 points), so every chunk is one line. The halving planner itself is exact (0 of 576 rules over budget). The budgets are re-measured in 3f.
  - **The folding of `_to_end` into the triangle's `own` read at `upper`** (an elegance NIT): not done.

**The user's ruling on the grading ratio** (2026-10-12, asked after elegance S1): **derive the ratio from the integrating rule's node count.** This supersedes the same day's "resolution-independent law" on this one point. The measurement that prompted it: at any fixed ratio the error is the node count's. At ratio 2 the block's inner rule missed by 1.0e-6 at 4 inner points, 5.8e-13 at 8 and 3.3e-15 at 12 (`elegance/probe_shape.py`); at ratio √2 the one-level reads missed by 5e-9 at 4 points. Rejected: one constant √2 everywhere (the 4-point error stays 5e-9), and keeping 2 and √2 with corrected docstrings.
- **The law** (`grading.ExponentialGrading.of(points)`).
  - The n-point Gauss-Legendre error on e^{−x} over [d, d+h] is at most c_n h^{2n+1} e^{−d}, where c_n = (n!)^4 / ((2n+1)((2n)!)^3).
  - h_n = (ε/c_n)^{1/(2n+1)} is the widest piece at the attenuated end that n points integrate to double precision.
  - The first piece is min(h_n, 2) mean free paths; a stretch no wider is not cut.
  - The geometric pieces have their largest bound at d = 2n+1, so ρ = 1 + e·h_n/(2n+1).
  - The cap of 2 is the nodes' own rounding: an uncapped first piece of 7.9 to 32.6 mean free paths, at 12 to 24 points, missed by 4 to 40 ε.
- **Measured** (`[M]` `ratio/law_capped.log`, the integral of e^{−x} over [0, 64] in double precision):

  | n | first | ρ | pieces | error |
  |---|---|---|---|---|
  | 4 | 0.194 | 1.059 | 103 | 2.7e-15 |
  | 8 | 2 | 1.419 | 11 | 0 |
  | 12 | 2 | 1.856 | 7 | 0 |
  | 16 | 2 | 2.238 | 6 | 2.2e-16 |

  The fixed laws it replaces: at 4 points, 1.1e-6 (ratio 2) and 5.0e-9 (√2); at 8 points, 7.4e-13 and 2.2e-16, with 7 and 13 pieces.
- **Applied to:** every rule whose Gauss nodes are in the exponent's own variable. That is the three piece sets (at `points`), the block's inner rule (at `inner_points`) and the impact rule in y (at `points`). `grading.DYADIC` (first 1, ratio 2) stays for the polar and cosine grazing ends: their nodes are in θ or μ, where the exponent τ/v is not linear. Retired: `ONE_LEVEL_RATIO`, `_THIN`, `_DOUBLINGS`, `exponential_ends`' `ratio`, and `_Toward`'s ratio.
- **The outflow under the derived law** (`[M]` `fix/outflow_accuracy.log`, a three-region cylinder under a partial mirror, against 24 points):

  | points | in-plane family | the 3-D line |
  |---|---|---|
  | 4 | ≤ 8.2e-15 (7.4e-13 at sin θ = 0.01) | up to 6.6e-8 |
  | 8 | ≤ 7.1e-15 | ≤ 4.3e-13 |
  | 16 | ≤ 3.6e-15 | ≤ 3.0e-14 |

  The 3-D line's residual at 4 points is the closest approach's curvature (the hp grading), not the exponential.
- Pending: the block's cost at 4 and 8 points (`fix/ladder_derived.log`); the test-architect's gates and re-baselines; the archivist's pages. The four `TransportResolution(4, 4, 4)` fixtures now pay 103 pieces per graded end.

**Two further user rulings on the derived law** (2026-10-12):
- **Along lines only.**
  - **What I had done.** I applied the derived law to the impact rule's rim grading in y as well; that was my extension of the ruling, not part of it.
  - **What it cost** (`[M]` the archivist, the three-region cylinder's impact rule; working tree against `4ebfa7a4`):

    | line points | nodes, derived law | nodes, before |
    |---|---|---|
    | 2 | 15 066 | 26 |
    | 3 | 1518 | 39 |
    | 4 | 368 | 52 |
    | 8 | 104 | 104 |

    The block's cost probe at 4 points was killed (exit 137) under that law.
  - **Ruled:** the impact rule returns to `grading.DYADIC`. The line quadrature is the reference's convergence parameter, so its coarse rules stay coarse. The derived law applies along lines: the three piece sets at `points`, and the block's inner rule at `inner_points`.
- **The first piece is capped at 1 mean free path, not 2** (`grading._WIDEST_FIRST`).
  - **Why the cap is measured, not derived.** The Gauss bound sees e^{−x} alone, but the integrand is a basis polynomial times it.
  - **The red that prompted it.** At a cap of 2, the test-architect's P1 rows read 2.5e-13 at the wall and 1.2e-13 at an interface, against their 1e-13 band at 8 points; at a cap of 1 they are within it.
  - **Refuted:** a degree-aware derivation. The safe width depends on how fast the polynomial varies over its panel, which the piece does not know (`ratio/first_poly.log`, eps-level and noisy).
  - **The table** (`ratio/law_cap1.log`): 13 pieces at 8 points, 8 at 12, 7 at 16, all exact; 103 at 4 and 23 at 6.

**Step 3e closed** (2026-10-12): code `c9db13e3`, docs `1ad1f557`.
- **Gates.** The test-architect's gates landed, each with a mutation red (`scratch/characteristic_architecture/p3/step3/gate_spec.md` §9). The characteristic fast tier is green (`ta/confirm2/summary.log`, 11 files). K4 went from 523 s to 10 s and is back in the fast tier.
- **Re-drops.** ERR-105, ERR-106 and ERR-035, re-dropped on the landed law, each red their catchers.
- **Issues still to file:**
  - the slab's own speed option (its line at cosine μ is its normal line at speed 1/|μ|);
  - `at` taking a crossing (a level and a side) instead of a parameter (elegance C2);
  - #587: comment on what the lift retired; it stays open for the polar weight absorbing sin²θ.

⏸ **Next: 3f.** Re-time the slow reference files with the §8 protocol, and re-measure the chunk budgets under the derived law; the reading's 0.122 s target is to be re-measured there. Then step 2.

## P3 step 3f (2026-10-12): the landed reading re-timed, and the product-integration prototype

**Landed:** step 3 merged to `main` at `0a058a2a` (`c9db13e3` code, `1ad1f557` docs). The full `-m "not slow"` suite passed before the merge: 16 882 passed, 0 failed. CI passed on `1ad1f557` (run 38331987610). The head commit's `[skip ci]` skipped the whole push, so `1ad1f557` was pushed to a temporary branch to run CI; that branch is deleted. Lesson: a `[skip ci]` plan commit pushed on top of code commits skips their CI too.

**The §8 protocol on the landed reading** `[M]` 2026-10-12 (`scratch/characteristic_architecture/p3/timing/point/`, `driver_c_2_1_2_4_4.log`; 15 alternating subprocesses per side, 1 warm-up and 5 timed calls each, `OMP_NUM_THREADS=1`, `-O`; load average 14 at the start and 8.5 at the end; old tree a detached worktree at `a336bde4`, new tree a detached worktree at `0a058a2a`):

| | old hoist | landed `main` |
|---|---|---|
| point reading, min (median) of the minima | 0.120 s (0.123 s) | 2.72 s (2.77 s) |
| solve at the matched rule `c_2_1_2_4_4` | 7.0 s | 284.7 s (11.6 s at step 1) |
| max φ_g error, φ₁/φ₀ error, against the panel-4 truth | 1.19e-3, −3.1e-5 | 9.3e-4, −1.2e-3 |

**The cause** (`profile_reading.py`): the matched rule integrates each piece with 4 Gauss points, and the derived law at n = 4 has first depth 0.194 and ratio 1.059. Each slot is cut into 254 to 266 ends, and the time is in the piece integrals (`_PieceSet._integral`, `_attenuated_length`, `at`). Nodes per grading (pieces × n over 64 mean free paths times the speed spread S): n = 4: 412 (S = 1) and 652 (S = 30); n = 8: 104 and 184; n = 12: 96 and 168.

**The diagnosis, offered to the user** `[R]`: the step-3 side effects share one cause. One piece mesh along each line resolves two unrelated scales: the geometry (the basis, and the branch points of r(t) = √(b² + t²)) and the attenuation e^{−vΣt}. Gauss–Legendre integrates polynomials, so the exponential has been resolved by cutting the line, and every fix since #587 was that mesh reacting to a speed, a node count or a spread. The proposed fix is **product integration**: interpolate the smooth factor at Legendre nodes shared across every speed, and integrate the attenuation exactly through the Legendre moments ∫₋₁¹ P_k(t) e^{−xt} dt = 2(−1)^k i_k(x) (i_k the modified spherical Bessel function). The error is then the interpolation error of the smooth factor, uniform in x. Checked `[M]`: the identity holds to about 1e-16 against mpmath for k = 0 to 6 and x from 1e-3 to 40. scipy's `spherical_in` loses accuracy at small x and high k (100 % wrong at k = 6, x = 1e-3; 1e-11 at k = 3), so a stable moment evaluation is part of the design.

**The user's ruling** (2026-10-12): caching makes the runtime bearable, but the question was whether a more fundamental fix exists. "Run the scratch prototype and use 2 cases for it. Pick an average case, and pick the worst offender (which is probably a bothersome cylinder)." Step 2 (the file re-timing) waits for the prototype's outcome.

⏸ **Next:** the prototype's report (numerics-investigator, `scratch/characteristic_architecture/p3/product/`), then a design section for the user's ruling before any production code.

### 3f, the prototype's result (2026-10-13)

The numerics-investigator's prototype is in `scratch/characteristic_architecture/p3/product/`. It ran in-process on `0a058a2a`; no file under `orpheus/`, `tests/` or `docs/` changed. Its report is in this section: the harness refused the agent's `report.md`. I re-read the check logs named below.

**The method** `[M]`.
- **Pieces.** Only the geometric ends: the slot, its middle, and production's hp branch halvings, packed by production's `_Packing`. The three piece sets collapse into one.
- **One level (the reading, the start and end weights).** The smooth factor is interpolated at n Gauss–Legendre nodes. The weights are Σ_k C[j,k] · 2(±1)^k ĩ_k(x), with:
  - x = Σvh/2, the half optical width of the piece;
  - ĩ_k = e^{−x} i_k, the scaled modified spherical Bessel function;
  - C, the Legendre coefficients of the Lagrange basis.
- **Two levels (the block's Volterra triangle).** The triangle reduces to the same moment problem: T_ab(x) = ∫₀² C_ab(τ) e^{−xτ} dτ. C_ab is the correlation polynomial of the two interpolants, of degree 2n − 1. It is precomputed once per n and needs the moments up to k = 2n − 1.
- **Reused unchanged.** The scan (`_carried`), `at`, `outflow`, `entry_response` and `LinePeriod.inflow`. The rule is injected through `lines.TraversalRule`.

**The moment evaluation** `[M]` (`moments.py`, `moments_check*.log`).
- **Below a threshold** X(K) = K(K−1)/ln 10: a positive power series gives the top ratio, and the continued-fraction recurrence carries the ratios downward, normalised by ĩ_0 = −expm1(−2x)/(2x).
- **Above the threshold:** forward recurrence from the closed forms of ĩ_0 and ĩ_1.
- **Why this scheme:** it has no cancellation and no overflow.
- **Accuracy against mpmath:** ≤ 2.13e-15 relative for k = 0 to 31, at x = 0, at 89 points over [1e-8, 1e3], and at 1e4 and 1e5. On 1507 dense points across the threshold edges for K = 12 and 16, ≤ 1.42e-15.
- **Mutations:** both arms red (inf): unscaled moments, and forward recurrence everywhere.
- **Premise refuted** (2026-10-13), for the question "is scipy's `spherical_in` unusable at small x?": on scipy 1.17.1 it is accurate to 1.5e-15 for k ≤ 31 and x ≤ 700. The 3f note above is wrong on this point; its probe misread scipy. The real hazard is overflow past x ≈ 709, which the scaled moment avoids.

**The weights** `[M]` (`weights_check.log`).
- One level, against mpmath: ≤ 7.7e-14 at n = 8, relative to ∫|L_a| e^{−xt}.
- The triangle, against a direct two-dimensional composite rule: ≤ 2.6e-15 at n = 4 and n = 8, for x from 0 to 1000.

**The control (X1)** `[M]` (`control_mutations*.log`). The product weights, applied to production's own graded pieces, reproduce production: 0 relative on the reading and 2.0e-16 on the block. Mutation arms:

| arm | reading | block |
|---|---|---|
| E and S swapped | red, 2.9e-5 | |
| moments k ≥ 1 zeroed | red, 1.4e-5 | |
| unscaled moments | red, NaN | |
| the sign (−1)^k dropped in the start weights | green | red, 5.3e-3 |
| the sign dropped in the triangle | | red, 2.5e-2 |
| the triangle transposed | | green |

The two green arms are designed-green:
- The reading reads no start weight, so the start-weight arm cannot change it.
- Every chord of a concentric body is a palindrome in r, so reversing a traversal leaves the sum unchanged. A transposed triangle is invisible on this geometry. A gate for the triangle's orientation therefore needs a line that is not a palindrome: a partial chord, or a point inside the body.

**Case 1, the average case: the A|B|A cylinder, φ(1.2), rule `c_2_1_2_4_4`** `[M]`.
- **Protocol:** the §8 protocol; load average 7.7 to 14.9.
- **Errors:** "err" is measured against production at 16 points on the same lines, groups and emission. It isolates the integration.
- **Truth:** every row reads 9.27e-4 (φ) and −1.22e-3 (φ₁/φ₀) against the panel-4 truth. That error is the line rule's, not the integration's.

| | reading, s | solve, s | reading err | block err | k against landed |
|---|---|---|---|---|---|
| landed, 4 points | 2.753 | 285.9 | 3e-16 | 5e-16 | 0 |
| production, 8 points | 0.578 | 26.3 | 1e-16 | 7e-16 | 3e-15 |
| product, n = 4 | 0.123 | 1.363 | 1.3e-9 | 1.5e-6 | 1.1e-9 |
| product, n = 8 | 0.146 | 1.704 | 4.7e-15 | 6.3e-13 | 1.0e-14 |
| old hoist `a336bde4` | 0.120 | 7.0 | | | |

**Case 2, the worst offender: the MR3 three-region vacuum cylinder read on its wall, rule (3, 0, 16, 12)** `[M]`.
- **How it was chosen.** A census of 8 cylinder fixtures, each with an mpmath or closed-form reference (`census_production*.log`). This case is the costliest on both the block (555 s, 8.8 M slots, 169 M pieces) and the reading (7.45 s). The runner-up block is the white cylinder at τ = 30, at 407 s.
- **The reference that isolates the integration.** The mpmath reference cannot: the fixture's source is constant per region, so every n reads ≤ 1.6e-15 against it. The isolating reference is production at 16 points on the same lines, compared on every basis coefficient (the "row").

| | block, s | reading, s | block err | row err |
|---|---|---|---|---|
| production, 12 points | 572.8 | 7.538 | 2.6e-16 | 3.3e-16 |
| production, 8 points (one run) | 674 | 11.2 | | |
| product, n = 8 | 11.6 | 2.209 | 2.2e-12 | 3.8e-11 |
| product, n = 12 | 15.51 | 2.583 | 1.0e-15 | 2.3e-15 |

The production rows ran 2 or 3 single calls each, not the full protocol.

**What the result says** `[M]`, read with the logs:
- **Product integration takes the attenuation out of the piece mesh, for the reading and for the block.** At matched accuracy (n = 8 on case 1, n = 12 on case 2) the cost drops:

  | | landed | product | factor |
  |---|---|---|---|
  | case 1 reading | 2.75 s | 0.146 s | 19 times |
  | case 1 solve | 286 s | 1.70 s | 168 times |
  | case 2 block | 573 s | 15.5 s | 37 times |
  | case 2 reading | 7.54 s | 2.58 s | 2.9 times |

  Case 1's reading (0.146 s) is within 22 % of the old hoist's 0.120 s, and its solve is 4 times faster than the hoist's 7.0 s.
- **The win is cost, not accuracy, and it changes what the node count means.** The landed reading at 4 points was exact along lines only because its exponential mesh over-resolved the geometry too: 1 145 832 pieces against 7 704 geometric pieces on case 1's reading. On the geometric pieces alone, n is a real convergence parameter: n = 4 leaves 1.5e-6 on case 1's block.
- **The derived law's cost is not monotone in n.** On case 2, production at 8 points costs more than at 12 (674 s against 573 s), so "just raise the points" is not a fix.

**What retires if it lands** (the agent's list; the audit owed under `retirement-audit` G.24 has not been run):
- `ExponentialGrading.of` (the node-count law) and `_WIDEST_FIRST`;
- `exponential_ends` along lines, `_Toward`, and the three piece sets;
- `_attenuated` (the graded inner rule), `_attenuated_length` and `_PieceSet._integral`;
- `TransportResolution.inner_points`: the triangle no longer has an inner rule of its own. This changes the content identity, so every memo key moves once.

**What stays:**
- `DYADIC`, `VANISHING_DEPTH` and `exponential_ends`, for the rules over lines: the polar and cosine grazing ends, and the impact rule's grading toward the rim. These coordinates are not along a line and the attenuation's exponent is not linear in them.
- The hp branch halvings: product integration's accuracy rests on them.

**Rejected candidates**, each for its question:
- **n = 4 kept**, for matched accuracy: the block error is 1.5e-6.
- **"Use 8 production points"**, for cost: cheaper on case 1 and dearer on case 2, as above.
- **Dropping the slot-middle cut**, for cost: 1.9 times fewer pieces and 300 times worse at n = 8.
- **A nested mpmath quadrature for the n = 8 triangle**, as the triangle's check: too slow; replaced by the direct composite rule.

**Open before production code:**
1. **The ruling on what `points` means** (below).
2. **Not run end to end with the product rule:** a hollow body's reading (the `to_start` path) and the partial mirrors (C2b).
3. **The moment evaluation** is about half of case 2's reading (`profile_case2_read.log`). The chunk budgets need re-measuring under the new rule.
4. **Gates for test-architect:** the moments against mpmath; the weights; the control; and a triangle-orientation gate on a line that is not a palindrome. No ERR entry: no defect was found. The scipy premise was a probe error, not a production defect.

**The design question for the user: what `points` (n) means once the attenuation is exact.** With the exponential exact, n integrates only the smooth factor: the basis polynomial times the geometry, analytic on each hp piece. Two readings:
- **(a) n derived from the hp law** `[R]`. The halvings make each piece no wider than its distance to the nearest branch point. In the piece's coordinate on [−1, 1] the singularity therefore lies at 3 or more, which is a Bernstein ellipse of parameter ρ ≥ 3 + √8 ≈ 5.83. Gauss–Legendre's error then falls like ρ^{−2n}, and double precision needs 2n ln ρ ≥ ln(1/ε), so n ≥ 10.4 for a square-root branch point. This is the same shape as the law it replaces (a law derived from double precision), and it agrees with case 2's measured n = 12.
  - It is conservative for case 1, which reached 4.7e-15 at n = 8.
  - The basis degree adds to the count.
  - The rule would be: n = ⌈ln(1/ε) / (2 ln ρ)⌉ plus the basis degree's share. `points` then leaves the resolution, and the memo keys move once.
- **(b) n stays a ladder axis.** The rung names n, and convergence is shown by the ladder against the independent truths, as now. The cost is that a coarse rung (n = 4) is no longer exact along lines, so the low rungs measure the along-line rule as well as the line rule.

My recommendation is (a), with the derivation measured before it is adopted. The measurement: the row error against the hp ratio and n on both cases, and the predicted n against the measured one. The project's line is to derive a law where one exists; (a) keeps the ladder's rungs meaning the line rule alone; and the derivation names its premise (a square-root branch point at a distance of at least one width), which a gate can check.

⏸ **Next:** the user's ruling on (a) or (b). Then a test-architect gate spec (the four gates above), the main agent's carve (W3, surgical: the reference's integration law), and the retirement audit table at `.claude/plans/p3_product_integration_dependency_audit.md`. Step 2, the file re-timing, follows the carve.

**The user's ruling** (2026-10-13): **(a)**, n derived from the hp halving law. `points` leaves the resolution once the derivation is measured.

**The measurement owed before adoption** (the numerics-investigator, resumed):
- the row error on both cases against n and against the halving tightness (width ≤ distance/c, for c = 1, 2 and 4);
- the predicted n against the measured n at 1e-14;
- whether the basis degree's share is additive.

In parallel, the test-architect writes the gate spec for the moments, the weights, the control and the triangle's orientation on a line that is not a palindrome.

⏸ **Next:** both reports; then the law fixed in the plan, the dependency audit table, and the carve.

### 3f, the law and the gate spec (2026-10-13)

**The law, from the numerics-investigator** `[M]` (`product/contour_bound.py`, `contour_p3c2.log`, `contour_x100.log`, `hp_*.log`):
- **The node count is the classical contour bound for interpolatory product rules** (Davis–Rabinowitz; Sloan–Smith). n is the smallest n with min over ρ′ < ρ⊥(c) of (1/2π)∮|f||K_n(x, z)||dz| ≤ ε∫ω|f|.
- **Its inputs, none fitted:**
  - x, the piece's largest optical half-width over the speed family;
  - ρ⊥(c), from the branch point perpendicular to the piece's end, at z = −1 ± 2ic;
  - f, the basis on the piece.
- **The rate depends on x:** about ρ^−2n when the piece is optically thin and ρ^−n when it is thick. On a synthetic piece the bound's n equals the measured n at x ≤ 10.
- **Refuted:** an algebraic n^(−p/2−1) prefactor in p. The dependence on p is the basis's amplitude on the ellipse; it replaces the fitted 1.5p.
- **Even pieces** (the panels touching a solid body's centre) are exact at 2p + 1 nodes: the basis there is a polynomial in (b² + t²).
- **The tightness:** c = 2 is recommended. The slot-middle cut is rejected: it costs more nodes than tightening c.

**The gate spec, from the test-architect** (`product/gate_spec.md`, drafts in `product/ta/`):
1. **The moments overflow at K ≥ 48** (n ≥ 24): the series overflows for x between 709 and X(K). A defect, to be fixed.
2. **The one-level rule misses its closed form on constants by up to 340 ulp.** Projecting onto exact discrete orthogonality brings it to ≤ 3.5 ulp.
3. **The componentwise leg of T7 reads 3.5e-13 against its 1e-13 band**, with a source that vanishes at the wall. `[HYPOTHESIS]` It is the same defect as the investigator's 3e-14 floor at x = 100: signed product weights on thick pieces, whose rounding scales as ε·Σ|w||f|. The investigator is measuring Σ|w|/∫ω and comparing positive-weight routes. **No tolerance is edited.**
4. **Only 2 of the 81 assembled rows see the triangle's orientation:** the two slab rows of T5. Every line set is closed under reversal. The new gate PT1 tests the orientation at the weight level.
5. **No existing row tells n = 8 from n = 12.** The law is pinned by a self-convergence row, PL1, owed with the law.
6. **The control is a one-time carve check, not a gate** (retirement-audit D.14).
7. **Existing gates:** the delete, move and keep lists are in the spec. `characteristic-exponential-grading` has no `verifies` marker.
8. **The ERR re-drops** (ERR-099, ERR-100, ERR-105, ERR-106, ERR-035): each still has a catcher that reddens under the prototype.

⏸ **Next:** the investigator's report on the 11-point table, the stable positive weights (items 1 to 3) and the law's own cost. Then the choice of weights goes to the user, then the carve.

### 3f, the weights: the measured choice (2026-10-13)

The numerics-investigator's item 2 and the law, verified against its logs by the main agent (scratch `p3/product/`; script `routeB_check.py`).

- **The defect of the interpolatory weights is their sign structure** `[M]` (`abs_weights.log`). Σ|w_j| / ∫ω is 1.00 for x ≤ 10 and grows toward the endpoint Lebesgue constant Σ|L_j(1)| (4.5 to 7.9 at n = 8 to 20) for x ≥ 100: 3.5 at x = 100, n = 8; 6.6 at x = 1000, n = 20. This is a property of the rule, not of how its weights are computed. It is the cause of the 3e-14 floor and, `[R]`, of T7's componentwise 3.5e-13.
- **Route B, the Gauss rule of the attenuation weight** (`product3.py`). Per (piece, speed), the m-point Gauss rule of e^{−x(1∓s)}, from Gautschi's modified Chebyshev algorithm and Golub–Welsch: Legendre modified moments for x < m, truncated-Laguerre moments for x ≥ m (`weights_routes.log`, an overlap measured stable). The weights are positive and rescaled to the exact zeroth moment. The source is evaluated at the rule's own nodes, so nothing is interpolated.
  - Against mpmath `[M]` (`routeB_mp.log`): one level 5.5e-15 converged, 1.8e-14 worst (x = 1000, m = 16); two levels (the triangle, relative to ∫|u||v|) 6.7e-16 converged, 1.9e-15 at m = 24, x = 1000. Route A made stable reads 9.2e-14 on the triangle at n = 24, x = 1000 (`weights_stable_check.log`, which measures route A only).
  - Constants: within 2 ulp at every x ∈ [0, 1e3] and m ∈ [2, 24] `[M]` (`constants_B.log`). The earlier "5 ulp" was never logged and is withdrawn.
  - A source vanishing at the wall, x = 1000: 1.8e-14 with the nodes formed through absolute parameters, ≤ 4.4e-16 with the nodes held as offsets from the attenuated end `[M]` (`componentwise.log`, synthetic piece). The offset form needs a production verb, the chord's point-relative orbit coordinate (`RadialImage.orbit_coordinate_at(distance)`). Not yet measured on the test-architect's T7 wall row.
  - Cost: nodes depend on the speed, so the source is evaluated F·m times per piece instead of n. Shared nodes across speeds are lost.
- **The moment overflow at K ≥ 48 is fixed** in `moments.py` (the two series summed and rescaled together): ≤ 4.0e-15 at K = 48 and 64 over x ∈ [1e-6, 3000] `[M]` (`moments_highK.log`).
- **Rejected:** splitting thick pieces at a depth derived from ε. It reintroduces the attenuation mesh that product integration removes, and route B is stable at every x.
- **The law under route B.** For a positive weight the error is at most 2∫ω·E_{2m−1}(f), with E_k(f) ≤ 2M(ρ′)ρ′^{−k}/(ρ′ − 1). m = ⌈(k + 1)/2⌉ with k the least integer such that min over 1 < ρ′ < ρ_sing of 4M(ρ′)ρ′^{−k}/(ρ′ − 1) ≤ ε. m does not depend on x. ρ_sing is from the piece's branch points; M is a closed-form majorant of max |u| on the ellipse (squared for the triangle); an even panel is entire, m = p + 1 for a row and 2p + 1 for the triangle. Measured / sampled bound / closed form, smallest m at 1e-14 (`routeB_sweeps.log`, `routeB_rows.log`, `predict_sep.log`):

  | | c = 1 | c = 2 | c = 4 |
  |---|---|---|---|
  | case 1, p = 2 | 11 / 14 / 16 | 9 / 11 / 12 | 8 / 9 / 10 |
  | case 2, p = 3 | 12 / 15 / 17 | 10 / 11 / 13 | 7 / 9 / 10 |
  | case 2 rows, p = 1 | 9 / 12 / 12 | 7 / 9 / 9 | 6 / 7 / 7 |
  | case 2 rows, p = 5 | 12 / 14 / 16 | 10 / 11 / 12 | 7 / 9 / 10 |

  0 of 11 cells underpredict. The closed form is never below the sampled bound (0 of 15 872 pieces) and costs about 200 flops per piece, once per piece. A sampled contour is not a safe bound (a coarser sampling read below the fine grid on up to 72 pieces). Per-piece m saves 25 to 31 % of nodes against one global m (c = 2: case 1 row 129 against 180, block 93 against 135; case 2 row 56 against 75, block 67 against 91). Not run: p = 5, c = 1, block.

⏸ **The user's ruling on the weights** (route B against a normwise contract or an end-layer cut under route A).

**The user's ruling** (2026-10-13): **route B**, the Gauss rule of the attenuation weight e^{−xt} per (piece, speed), with the nodes held as offsets from the attenuated end through a new chord verb `orbit_coordinate_at(distance)`, and the node count m from the derived law above. The normwise contract and the end-layer cut are not taken; T7's 1e-13 componentwise band stands.

⏸ **Next:** the test-architect re-points the gate spec at route B (T7 under route B with offset nodes, PL1 for the law with its m − 2 red, the kernel gate for `orbit_coordinate_at`, the weight-level orientation row PT1). Then the G.24 dependency table at `.claude/plans/p3_product_integration_dependency_audit.md`, then the W3 carve.

### 3f, the gate spec under route B (2026-10-13)

The test-architect revised `scratch/characteristic_architecture/p3/product/gate_spec.md` for route B (route-A drafts kept in `ta/superseded_routeA/`). Findings, each `[M]` on the prototype under `-O`:

- **Blocking: the moment-route switch is unstable at m ≥ 16 near x = m.** 1805 eps at m = 16, x = m; for m ≥ 20 no stable branch exists between about 0.8m and 1.5m. The earlier overlap measurement never sampled x = m, and the law asks m = 16–17 at c = 1. Sent to the numerics-investigator: a construction stable for every x and m ≤ 24, no cap on m and no fitted switch.
- **Rows.** PW1–PW3, the projection row and the tensor rows die with interpolation. RB1 (positive weights, distinct interior nodes), RB2 (constants within 2 ulp; red without the rescale), RB3 (exactness to degree 2m − 1 against mpmath; red when the switch moves), and a new PT1 (triangle orientation on non-symmetric cubics, 4 of 4 red when transposed) replace them. PM1–PM3 survive under the Legendre route; PM1b (x = 5e-324) waits on whether that input is realizable.
- **The offset verb already exists.** `Rays.orbit_coordinate_at(distance)` (`orpheus/geometry/chord.py`) forms the coordinate from the point, never through the parameter, on both images. `RadialImage.orbit_coordinate_at(t)` takes a parameter, so the proposed name was a homonym. The parameter round trip costs 8e13 to 4.5e15 eps on the slab face and 0.7 eps on radial bodies, so the slab needs the offset form too. `[R]` The carve reads each piece's nodes through a `Rays` minted at the piece's attenuated end; the open question is how a `Rays` is minted at a crossing (the "at taking a crossing" item already on the issue list).
- **PL1.** The m − 2 arm reddens 1 of 4 fixtures at c = 2, because the closed form exceeds the measured m by 2 to 5. The row needs an arm below that margin or a margin assertion. An even panel needs two counts, p + 1 for the row and 2p + 1 for the triangle; the triangle at p + 1 also reddens 7 existing value rows.
- **T7 holds at 1e-13 under route B with every node placement** (offsets ≤ 4.4e-16, absolute 2.8e-14). So T7 does not catch the offset verb (KG1 does) and no end-to-end row catches a missing rescale (RB2 does). T5, T5c and T11 are green on the installed prototype.
- **ERR re-drops:** ERR-105, ERR-106, ERR-035, ERR-100 (re-spelled) and ERR-099 are still caught. F5's ERR-106 marker stays green, as before.
- **The hp walk never terminates on a cavity of radius 1e-12** (QA4 red). Sent to the investigator.
- Owed: a weight-level verb for the triangle so PT1 tests production; a route-B first red for AS2.

⏸ **Next:** the investigator's stable construction, the hp-walk fix and `law.py`; then PL1 finalised, the dependency table, the carve.

### 3f, the route-B blockers closed (2026-10-13)

The investigator's report, checked against its logs (all in `scratch/characteristic_architecture/p3/product/`):

- **The rule, with no switch.** A discretised Stieltjes procedure on the weight in the offset u = x(1 − s) from the attenuated end (`rule_stieltjes.py`). The measure is discretised by Gauss–Legendre panels in u, truncated at U = ln(1/ε) + 10 ln 2; the truncation is derived (every Legendre polynomial in s is bounded by 1, so the dropped tail is below e^{−U}). For 2x ≥ U the truncated measure no longer depends on x, and its rule is one table per m, generated in mpmath (`gen_const_rules.py`, 2m moments to 1e-71 or better for m ≤ 32): generated, not fitted. `[M]` `rb_dense_offsets.log`: RB1 8448 of 8448 (m = 1 to 32, x from 0 to 1e3 plus 61 points in [m/2, 2m]); RB2 ≤ 1.9 ulp; RB3 ≤ 42.4 eps with nodes held as offsets, against the 64-eps band; 0 failures.
- **RB3 must be written on offsets.** Held in s, the same rule reads up to 93 eps (m = 20, x = 1000; `rb3_compare.log`); the rule evaluated exactly at its own double-rounded nodes already reads 66, so the loss is the representation, not the construction.
- **Rejected, each for a stated reason:** the modified Chebyshev algorithm (each moment basis has its own unstable band, and the two bands leave a gap at x ≈ m for m ≥ 16); the Jacobi matrix built in s (54 to 156 eps: the nodes crowd at s = 1 and the recurrence loses digits there); Newton refinement with Christoffel weights (no gain: the error was not in the eigensolver); a third moment basis (still a switch).
- **QA4, the hp walk on a cavity of radius 1e-12.** Root cause: the walk re-read the distance through the line parameter, which rounds the end of a cavity slot onto the closest approach, so the distance read 0 and the walk stalled. Fix: the distance in closed form, D² = (b/|PΩ|)² + (e + a)², with e from the chord's own crossings (forming e from the rounded b lost up to 7e-10 near tangency, ERR-106's mechanism); lines with b = 0 are entire; the step count is bounded at 1 + ⌈ln(c(L + d₀)/d₀)/ln(1 + 1/c)⌉. `[M]` `qa4_fix.log`: conservation 1.42e-15 at c = 1, 2, 4.
- **The law as a function** (`law.py`, `node_counts(...) -> (m_one, m_two)`): an even panel gets p + 1 and 2p + 1; a piece wider than its distance over c is refused (positive control: 68 refusals under the old walk). `[M]` `law_check.log`: equal to `predict.py` on every non-even piece at 12 of 12 (case, p, c) rows (3940 to 22648 pieces per row).
- **The margin**, closed form minus measured m, 1 to 5 per row; measured before the walk fix and with the old rule, so it is being re-measured on the new pair (`margin_B.log`).

⏸ **Next:** `margin_B.log`; the test-architect's RB3 on offsets and PL1 per row; then the G.24 dependency table and the carve.

### 3f, the rulings on the offset's home (2026-10-13)

- **The user's ruling: advance the rays.** "The point a distance d further along this line" is a geometry verb on `Rays`, not arithmetic in the reference: `Rays.advanced(d)` returns the same lines read from points d further along, re-anchored at the new point's own coordinate (on the cylinder and the sphere, the chord's image re-levelled at the new point's radius and half-chord; on the slab, the start shifted). Rejected: rays at crossings only (an hp-walk end lies inside a slot, on no crossing); the long distance from the observer (the offset's precision is lost to rounding against the end's distance).
- **The user's ruling: a crossing is anchored exactly, a slot end is advanced from its crossing.** Advancing from the observer by the summed slot lengths places a wall only to eps × (distance travelled), which near a slab face at grazing incidence is the failure KG1 measured for the parameter route `[R]`. So two verbs: `Chord.rays_at(crossing)` mints the lines read from one crossing per line, with the start, the breakpoint and (radially) the exact half-chord read from the chord's own crossing, so l = 0 is the wall exactly and KG1 stands as written; `Rays.advanced(d)` serves the segment ends inside a slot, advanced from that slot's crossing.
- The dependency table is `.claude/plans/p3_product_integration_dependency_audit.md` (the explorer, 2026-10-13; graph rebuilt the same day).

### 3f, the rule's home and the table retired (2026-10-13)

- **The user's ruling: the reference package.** The Gauss rule of the attenuation weight is a new module beside `grading.py` in `orpheus/derivations/continuous/characteristic/`, not in `numerics/`; the node-count law goes into `grading.py`.
- **The mpmath table does not ship** `[M]` (the numerics-investigator, `p3/product/table_vs_stieltjes.log`, the dense grid of `rb_dense_offsets.py`, 45 to 48 values of x per m). With the table disabled, the Stieltjes construction alone reads RB3 ≤ 61.4 eps over m = 1 to 40 (worst at m = 39; 0 of 40 rows above the 64-eps band). With the table, ≤ 42.4 eps over m = 1 to 32. Where both exist they differ by noise, not trend (Stieltjes better on 6 values of m, worse on 13, by at most 18 eps). The table is therefore a second definition of one quantity (X4), not a repair. The margin narrows for m ≥ 35 (49.5 to 61.4 eps); the law asks at most 17 at c = 1.
- **The batched build, table-free** `[M]` (`rb_batched_noconst.log`, the same grid, the batched construction the product prototype runs, without Newton refinement): read below in the carve's record.
- **Cost** `[M]` (`rule_cost.log`, min over 15 subprocesses of the min of 20 calls, load 11.7): one scalar build 0.63 to 0.67 ms at m = 16 and 0.92 to 0.97 ms at m = 24, independent of x.
- **PL1's safe leg is measured at the law's m** `[M]` (`margin_B.log`, the m = 15 to 20 extension in `mB_case1.log` and `mB_case2.log`): case 1 p = 2 c = 1 and case 2 p = 3 c = 1 read ≤ 1e-14 on block and row at every m from the measured value through 20 (worst 6.1e-15), which covers the law's 16 and 17.
- **The user's rulings on the carve's constants** (2026-10-13): the hp tightness is **c = 2**, fixed in the law (the investigator's recommendation; the node law asks 9 to 13 per piece there). **`TransportResolution` keeps one field**, `line_points`: `points` and `inner_points` leave it, so every memo key moves once.

### 3f, the carve written (2026-10-13, uncommitted on `refactor/p3-routeB-product`)

- **The rule** is `characteristic/attenuation_rule.py`: `AttenuationRule.of(half_width, count)`, nodes as offsets, zero-weight padding to the largest count, and `.triangle()` → `TriangleRule(lead, lag, weight)`, the weight-level verb PT1 asked for. `[M]` (scratchpad `probe_rule.py`): bit-identical to the prototype's batched build without the table at m = 1, 4, 9, 13, 17 over 12 values of x from 0 to 1000; the triangle on two cubics against mpmath at (x, m) = (0, 4), (7, 6), (1000, 7) ≤ 1.7e-16.
- **The law** is in `grading.py`: `TIGHTNESS = 2`, `ALONG_LINE_ERROR = 1e-14`, `branch_walk` (the closed-form walk), `bernstein_parameter`, `gauss_nodes`. The basis's bound on the ellipse is `PanelBasis.bound`. `ExponentialGrading.of` and `_WIDEST_FIRST` are retired; `ExponentialGrading` keeps only DYADIC's first depth and ratio.
- **A new finding, ruled** (the user, 2026-10-13, "Anchor at the closest approach"). A slot that ends at a line's closest approach has no crossing at its near end. Advancing from its far crossing to points near t* forms r² − b² from r_k² by cancellation, and `RadialImage.at_level` refused it on the cylinder `[M]`. The kernel gained `Chord.rays_at_closest_approach(points)`: level (b, 0), side 0. Every piece end that lies on a crossing is read from that crossing exactly. Every other piece end is advanced from its slot's near end, which is a crossing or the closest approach. The slab's slots are advanced from their start crossing.
- **A second kernel ruling** (the user, 2026-10-13, "Rays over a broadcast batch"). `Chord.rays_at(crossing)` takes `(*batch, *k)`, with the chord broadcast as views (`Chord._per_point`). `[M]` (`probe_k.py`): bit-identical to per-column calls on the sphere, the cylinder and the slab, 50 lines each, with k = 3.
- **Transport**: one piece set. `_PieceSet.node_counts` gives (m_one, m_two). `_Toward`, `_attenuated`, `_attenuated_length`, `_nodes`, `_orbit_coordinate`, `_branch_edges` and `_attenuations` are retired. `TraversalRule` lost `points` and `inner_points`. `Lines.chunks()`, `LineRule.transport(support)`, `PointRule.row(transport)` and `angular_flux(...)` lost them too. `TransportResolution(line_points)`.
- **Smoke check** `[M]` (`probe_t.py`; 40 random lines; Σ = (0.6, 13, 0.45); degree 3; against HEAD's rule at 16 and 16 points):

  | body | outflow | entry response | block | pieces | time |
  |---|---|---|---|---|---|
  | sphere | 7.6e-15 | 6.0e-15 | 1.1e-15 | 1096 vs 2122 | 0.11 s vs 0.11 s |
  | cylinder, speeds 1, 3, 40 | 4.9e-14 | 2.7e-14 | 1.5e-15 | 1118 vs 7066 | 0.27 s vs 0.89 s |
  | hollow cylinder | 6.0e-15 | 1.3e-14 | 2.8e-15 | 1210 vs 2684 | 0.11 s vs 0.21 s |
  | slab | 6.2e-15 | 6.2e-15 | 1.9e-15 | 600 vs 1662 | 0.02 s vs 0.14 s |

  The check does not establish which side carries the cylinder's 4.9e-14. The test-architect's mpmath rows decide it.

⏸ **Next:** the test-architect lands the gates (RB1 to RB3 and PT1 on the module, PL1 on the production counts, KG rows for several points per line and for the closest approach) and migrates the signatures. Then the characteristic and geometry suites; the reviews; the archivist; the memo regeneration; the full suite.

### 3f, the reviews and seven rulings (2026-10-13)

**The reviews.** The elegance-enforcer (`p3/product/review/elegance.md`) found 0 violations and five concerns. qa (`review/qa.md`) found the terms right by reading, and its in-process mutations of each term reddened every arm, from 1.8e-5 to O(1). The test-architect's migrated suite found three production defects. Each item below is fixed in the working tree.

- **An unbatched line crashed** in `rays_at`. Cause: the per-point spread skipped numpy scalars. A second cause sat in `TraversalRule._unflat`: `reshape(*())` fails on an empty batch. 98 of the suite's 116 reds `[M]`.
- **A line that makes no crossing aborted the batch**, where HEAD returned zeros. The padded entries asked the kernel to anchor points that do not exist.
- **The level guard tripped on a point rule's line** re-levelled at an inner crossing. The pair agrees with b only as well as the line's current level does.
- **The cylinder's 4.9e-14 is the orbit coordinate's rounding at the nodes**, amplified on the narrow outer panel [1.425, 1.5], and is not the rule's error `[M]` (qa: mpmath 0.72766774696014214; HEAD +3.8e-14 relative, the carve −5.4e-14; the rule at exact coordinates is within 2.8e-16; one ulp at every node bounds it at 8.4e-14). An mpmath row on narrow panels needs a band of about 1e-13.
- **Pieces stored from the slot's start** round-tripped the walk's near-end distances through L − a: 8.0e-4 relative at b/|PΩ| = 1e-13 `[M]` (elegance). Pieces are now signed distances from the slot's anchor (`_Slots.origin`, `first`, `last`).
- **Smaller items, each fixed.** An `at_upper` flag (now the `_End` enum). |PΩ| spelled twice. `gauss_nodes` outside its domain. ⌈(k+1)/2⌉ spelled in two places (now `exact_nodes`). The small-x branch built in s (now in the offset). The false claim 2^20/20! < ε (it is 4.3e-13). Stale prose.

**The rulings (the user, 2026-10-13):**
1. **Live ends only, through `Chord.take`.** Asked to reason from first principles, I found the padding was the defect: a dummy anchor (a foot, a borrowed crossing) places a point that does not exist. `Chord.take(lines)` gives the chord of a sub-batch. Transport reads only the live piece ends, in two groups: from a crossing, and from the closest approach. This reverses two earlier rulings, the broadcast batch (`_per_point`) and `Rays.where`, which was added and removed the same day. `rays_at` takes one crossing per line again.
2. **Zeros, as HEAD, for a line with no crossing.** It follows from ruling 1 by construction, on every chart.
3. **Counts above 40 are refused.** `attenuation_rule.MOST_NODES = 40` is a scope boundary. Basis degree 16 asks 41 on the triangle; m = 64 read 68.5 eps (qa).
4. **The law's target is the rule's band.** `ALONG_LINE_ERROR = 64 eps`. At 1e-14, PL1's safe leg read 1.09e-14 on the cylinder's speed family, under the rule's own floor (61.6 eps at m = 28).
5. **The re-orthogonalisation is removed.** It had no witness: removing it reddened no gate up to m = 80.
6. **The `rays_at` multi-crossing form and `Rays.where`** are superseded by ruling 1.
7. **The level guard's scale is max(r*², the current level²)**, a fix of the guard's reference scale.

⏸ **Next:** the test-architect re-points KG4 and KG5 onto `take`, adds the rows for an unbatched line, the level guard, `MOST_NODES` and `gauss_nodes`' refusals, re-measures PL1 at 64 eps, re-runs both suites and the ERR arms, and re-measures `ABA_CYLINDER_TRANSPORT_STEP`. Then the elegance-enforcer reviews the delta; the archivist; the memo regeneration; the full suite.

### 3f, the level guard's reference and the second and third review rounds (2026-10-14)

- **Ruling 7 is superseded by its root.** The max-scale fix widened the guard's scale but kept the wrong reference. A probe found the cause: a point-rule line at (0.3, y), re-read through `advanced(0)`, disagrees with its own b by 75 to 100 ulp, because b is itself a rounded quantity. `RadialImage.at_level` now holds a new pair to the line's **current** level, b² = (r_c − y_c)(r_c + y_c), within 64 ulp of max(r*², r_c²). The guard now checks that a step keeps the line on its level. It no longer checks agreement with a b recomputed elsewhere. Both measurements are in the method's docstring.
- **The test-architect's last report** `[M]`:
  - the characteristic fast tier, 1507 passed, 3 skipped, 0 failed;
  - geometry, 1615 passed, 5 skipped;
  - pyright 0.
  
  Every ERR arm reddens as designed. The transpose arm reddens only slab T5 and PT1, and the radial rows are declared blind to it. K7 has two first reds. KG1 to KG6 cover `take`, an unbatched line, the closest approach and the level guard.
- **The archivist's report** `[M]`: `sphinx-build -E -W --keep-going` is clean, and `dead_references` reads 0 of 65. Four new labels: `characteristic-attenuation-rule`, `characteristic-node-law`, `characteristic-branch-walk` and `characteristic-anchored-nodes`. `characteristic-exponential-grading` is retired.
- **Filed:**
  - #600: F5's `catches("ERR-106")` is not reddened by the re-drop.
  - #601: refs.bib entries for Gautschi 2004, Golub–Welsch 1969 and Trefethen's ATAP.

⏸ **Next:** the full suite, including the slow tier and the memo regeneration under the moved keys. Then the commit on `refactor/p3-routeB-product`, the merge and CI. Then step 2, re-timing the slow reference files (#591).

### 3f landed (2026-10-14) [LANDED 25b0c38f]

- **The full suite** `[M]` (`python -O -m pytest -p no:cacheprovider tests/`, serial, slow tier included, the traced memo regenerated under the moved keys): 17214 passed, 316 skipped, 60 xfailed, 0 failed, in 4 h 19 min. `sphinx -W --keep-going` is clean. `dead_references` reads 0 of 65, after a control planted on the characteristic page read 1.
- **The slowest characteristic gates** `[M]` (same run, `--durations=40`): the hollow cylinder's emission balance at 555 s and 525 s; the annulus convergence ladder at 450 s; the F4 frozen route at 409 s; the cylinder's emission balances at 377 to 389 s. Step 2 profiles them before choosing any speed-up.

⏸ **Next:** step 2. Re-time the slow reference files (`timing/run.sh`, #591), profile one slow cylinder file with cProfile, and ingest the profile into Nexus. The user asked whether JAX would help. `[R]`: not before the profile shows Python dispatch dominating, for four reasons:
  - float64 has no Metal backend, and on the CPU JAX calls the same LAPACK as numpy;
  - the pieces have ragged node counts, which defeat `jit`;
  - the mpmath tier is untouchable;
  - JAX would be churn in the closed reference layer.

### Step 2, the re-timing and the profile (2026-10-14)

**The re-timing** `[M]`:
- Run: `scratch/characteristic_architecture/p3/timing/run.sh`, 24 files, serial, `OMP_NUM_THREADS=1`, detached worktree at `941723dc`.
- The four passes:
  - bypass: memo plugin bypassed;
  - cold: `.cache/references` emptied;
  - warm1;
  - warm2: `-m "not slow"`.
- Started 02:15, ended 11:12, 96 rows. Every row exits 0, except 4 warm2 rows that exit 5 because `not slow` deselects every test in those files.
- Table: `timing/results.tsv`. Logs, with `--durations=0`: `timing/logs/`.

| pass | total |
|---|---|
| bypass | 11 601 s (3.22 h) |
| cold | 10 329 s (2.87 h) |
| warm1 | 9 307 s (2.59 h) |
| warm2 | 991 s (0.28 h) |

| file | bypass | cold | warm1 | warm2 | 2026-10-04 bypass |
|---|---|---|---|---|---|
| characteristic_assembly | 3559 | 3520 | 3535 | 108 | (new) |
| characteristic_system | 1114 | 1108 | 1105 | 41 | (new) |
| characteristic_lift | 1091 | 1088 | 1090 | 304 | (new) |
| characteristic_albedos | 1039 | 975 | 988 | 90 | (new) |
| characteristic_reading | 1031 | 1028 | 1032 | 131 | (new) |
| characteristic_convergence | 935 | 937 | 933 | 31 | (new) |
| phase_c_crosscheck | 820 | 1072 | 44 | (all slow) | 1112 |
| l1_standoff_slab_cylinder | 736 | 22 | 23 | (all slow) | 1846 |
| unified_matvec_cylinder | 722 | 4.7 | 4.9 | 1.9 | 848 |
| characteristic_independent_references | 213 | 221 | 213 | 17 | (new) |
| characteristic_transport | 167 | 167 | 166 | 168 | (new) |
| characteristic_reference | 107 | 107 | 107 | 55 | (new) |
| each of the other 12 files | ≤ 27 | ≤ 27 | ≤ 27 | ≤ 6 | |

**What the table says:**
- **The traced memo already finishes the older reference files.** Once the cache is warm (6020 entries), phase_c, l1_standoff and unified_matvec run in under a minute. Their bypass times fell 15 % to 60 % against 2026-10-04, from the P3 route-B carve.
- **The characteristic files cost the same warm as in bypass.** About 2.6 h a pass is computation inside the characteristic reference, which is the system under test of those gates. Memoising it would make the gates test the cache. `reference_cache.md` P4, ordering cache families by the hot files, therefore has no hot file left to order `[R]`.
- **The cost sits in the slow tier.** 21 of the 187 assembly tests carry about 3 450 of its 3 559 s (bypass minus warm2). The slowest single tests are cylinder cases:
  - emission balance: hcylinder 540 and 507 s, cylinder 369 to 378 s;
  - annulus convergence 435 s;
  - F4 frozen route 390 s.

**The profile** `[M]`:
- Command: `python -O -m cProfile` on `test_characteristic_assembly.py::test_a_closed_body_conserves_its_emission[cylinder_white-g0]`, at 941723dc.
- Runtime: 391 s under the profiler, against 369 s unprofiled.
- Artefacts: `timing/profile/assembly_cyl_white_g0.prof`, plus `.main.prof` with the paths remapped to the main tree. Ingested into Nexus as run `char_assembly_cyl_white_g0`: 361 nodes, 367 edges.

| where | cumulative | share | calls |
|---|---|---|---|
| `AttenuationRule.of` (rule builds) | 196 s | 50 % | 5 706 → 28 800 `_built` → 287 441 `_stieltjes` |
| inside: `gauss_legendre`, recomputed per `_stieltjes` | 59 s | 15 % | 289 625 `leggauss` |
| inside: `np.linalg.eigh` (Golub–Welsch) | 58 s | 15 % | 287 441 |
| inside: `_jacobi` (Stieltjes recurrence) | 53 s | 14 % | 287 441 |
| `CharacteristicBasis.values` (Lagrange product) | 116 s | 30 % | 7 609, 76 s self |
| `einsum` | 53 s | 13 % | 20 923 |

**The diagnosis** `[M]` from the counts, `[R]` for the remedy:
- **The rule builds are dispatch-bound.**
  - Each `_stieltjes` call costs about 0.65 ms of fixed overhead: one Legendre eigenproblem, one `eigh` call and a Python loop of m steps over small arrays.
  - The batches are small because `of` is called per line (from `transport._one_level` and the exit integral), and then split by node count and panel count.
  - Gauss–Legendre rules of a given order are rebuilt 289 625 times.
- **`basis.values` is compute-bound.** The Lagrange product form costs O(p²) per point.

**The JAX question, answered by measurement:**
- Python dispatch does dominate the rule builds: 287 441 small calls. That is the one place where `jit` would help.
- The same collapse is available in numpy by building the rules in one call per transport rather than one per line.
- The other 45 % is LAPACK `eigh` and memory-bound array arithmetic, which JAX on the CPU does not change.
- `[R]`: no JAX.

**Remedies, cheapest first** `[R]`, unmeasured:
1. **Memoise the Gauss–Legendre nodes by order**, in `derivations/common/quadrature._leggauss`, returning read-only arrays. Bit-identical. About 55 s of 391.
2. **Build the attenuation rules once per transport, over every line and level**, not per line. Batched `eigh` and the row-wise recurrence give each rule the same arithmetic as now, so the result should be bit-identical; that needs a measurement. It turns 287 441 calls into about one per (node count, panel count) pair. Most of the 196 s.
3. **Write `basis.values` as prefix and suffix products**, O(p) per point and still division-free. Not bit-identical, because the product order changes; the difference is at the ulp level and passes through the memo regeneration. About two thirds of the 116 s.

Estimated together: about 391 s down to 120–150 s per slow cylinder test, roughly 2.5×.

⏸ **Next:** the user's ruling on remedies 1–3, and whether 3's ulp-level move is accepted. Then a W3 carve with a bit-identity probe for 1 and 2, and a re-timing of the five hottest characteristic files.

### Step 2b, the remedies measured and the cache boundaries (2026-10-14)

The user's rulings of 2026-10-14:
- "I accept all 3 and ulp change that comes with the 3rd one."
- Then, mid-carve: "attenuation rules and reference results could be 2 different caches", and "Finding cache opportunities and decoupling the impact is an important aspect of speed up".

**The remedies, measured** `[M]` 2026-10-14. The probes are under `scratch/characteristic_architecture/p3/speed/`. Each replays the 1020 `AttenuationRule.of` calls, 909 024 rules, recorded from `test_the_escape_and_transmission_probabilities_are_the_closed_forms[cylinder-tau0.5]`.

1. **Gauss–Legendre memoised by order** (`_leggauss`, `functools.cache`, read-only arrays). The rule build fell from 13.4 s to 8.5 s, bit-identical on every rule (`variants.py`). Kept.
2. **One rule build per transport.** Not built.
   - `[REFUTED 2026-10-14]` as a major remedy: the ceiling, every recorded rule built in one call, is 5.85 s against 8.5 s (`aggregate.py`). It is bit-identical row by row. Once Gauss–Legendre is cached, the arithmetic dominates rather than the per-call overhead.
   - Padding every panel count into one discretisation was slower (9.3 s) and not bit-identical (23 eps).
   - Dropped.
3. **Prefix/suffix products in `basis.values`.** `[REFUTED 2026-10-14]`: 2.2 times slower at the working degree 3, faster only from degree 8 (`values_bench.py`).
   - The real cost was gathering the Lagrange tables per point by panel kind. Applying each kind's table to its own points is 26 % faster and bit-identical (`values_kind.py`). Kept.
   - Rejected: the modal form, `values = P(x) @ C`, the nodal-from-modal matrix built exactly from the Gauss rule's discrete orthogonality. It is 2.2 times faster, but its node values are off by 1e-15 to 1e-14 of either sign (`values_modal.py`), so the basis is no longer exactly one-hot at its nodes, which the dense pencil's single-sign test relies on. That is a semantic change, not the ulp move that was ruled.
- **Net:** the 27 s test takes 20 s (−25 %), and the results are bit-identical, so the accepted ulp change is not used.

**A correction to step 2.** That step said "memoising it would make the gates test the cache". `[REFUTED 2026-10-14]`. A traced-memo entry is pinned on every def that ran and is re-validated against the current source, so a hit is what the current code computes. A child entry is validated recursively, which is exactly the decoupling the user asked for.
- The constraint is granularity: a miss regenerates in a fresh interpreter, about 0.65 s.
- `AttenuationRule.of` is the wrong unit: about 8 ms per call, about 100 MB of rules per 27 s test, and its keys move with anything upstream.

**The boundary census** `[M]` (explorer, `scratch/characteristic_architecture/p3/speed/memo_boundaries.md`). A dry-run plugin ran over 14 files and 1368 tests; 13 failed, and their counts are lower bounds.
- **Where the time is:** about 7800 s of the 9172 s warm1 total is inside the reference stages, and 1350 s is test-side mpmath.
- **Sharing:** 0 of the top 30 tests share a block with another test, so the gain is from re-runs, not from sharing within one run.
- **Ranked boundaries:**
  1. **C1, one group's transport block:** a module-level function equal to `LineRule.of(...).transport(support)`, consumed by `GalerkinSystem.groups` and the test helpers.
     - It serves 5073 s with no gate change. Its payloads round-trip bit-exact, at 23 to 230 kB each and 11.7 MB in total.
     - Body edits to `system`, `reading`, `reference` and `dense_pencil` leave blocks valid.
     - A size threshold of at least 2000 line-speeds keeps 102 blocks (6101 s) and leaves 647 cheap blocks (186 s) unmemoised.
     - **It needs** `PanelBasis`, `ConcentricPartition` and `RigidMotion` to be keyable. All three are `eq=False` today and refused by the key.
     - The assembly test that monkeypatches `MeasureCoordinate.derivative` moves under `bypass()`.
  2. **C2, cylinder point readings:** 60 reads, 10 to 20 s each `[R]`.
  3. **C3, a bypass scoped per function:** the gates behind `bypass()` hold 2726 s. Needs C1 first.
  4. **C4, a block memo keyed at the door** (the fallback): it fits the key today but misses assembly and system (4345 s).
- **A limit:** the package `__init__` imports every submodule, so every memo pins all 13 characteristic module skeletons. Only edits inside function bodies regenerate one stage at a time; an edit to a signature, a constant or an import invalidates every entry.

⏸ **Next:** the user's rulings on the following.
- **(a)** Retire the step-2 claim, so that gates may read memoised blocks.
- **(b)** Content identity for the three types (value equality and the memo's key), or a `canonical=` recipe form.
- **(c)** C1 first, then C2 and C3.

Remedies 1 and 3 [LANDED 6b9fdd9d]. Then: C1 as a W3 carve with its X1 witnesses, then a re-timing.

### Step 2c, the rulings on the cache boundaries and C1 opened (2026-10-14)

The user's rulings of 2026-10-14 ("yes. I agree with your recommendations"):
- **(a)** Gates may read memoised blocks. The step-2 sentence "memoising it would make the gates test the cache" is retired (step 2b gives the reason: a hit is re-validated against every def that ran).
- **(b)** `PanelBasis`, `ConcentricPartition` and `RigidMotion` get value equality, through `ContentIdentity`, so the memo keys them directly. The `canonical=` recipe form is not taken.
- **(c)** The order is C1 (one group's transport block as a memoised function), then C2 (cylinder point rows) and C3 (a bypass scoped per function), then a re-timing.

C1 is a W3 surgical carve on a branch `perf/characteristic-block-memo`; the main agent writes the code. Before the code:
1. the explorer's census of the three types: every constructor, every `==`/`is`/hash use, every container or cache that keys on one, and every field that would be refused as a content part (`scratch/characteristic_architecture/p3/c1/census.md`);
2. the test-architect's gate spec (`scratch/characteristic_architecture/p3/c1/gate_spec.md`).

⏸ **Next:** the census and the gate spec; then the content-identity commit for the three types; then the block function and its consumers.

### Step 2c, C1 [LANDED 9178ad00] (2026-10-15)

Merged to `main` @ `9178ad00`. The full not-slow suite before the merge `[M]` 2026-10-15: 17 080 passed, 1 failed, 261 skipped, 268 deselected, 53 xfailed, 3 924 s. The failure was `test_closed_references_import_only_the_numerics_interface`: `assembly.py` and `walls.py` now import `orpheus.numerics.scalars`, which was outside `REFERENCE_INTERFACE_NUMERICS`. **Contract change:** `scalars` joined the set (the gate's comment, the layering page's table). It is the definition of a value that the content encoder applies, a leaf that `content` already imports, and it parses without computing anything a reference's value depends on. Copying it into `derivations/` would have been a second definition of a real number. After the change, `test_layer_imports.py` passed 471 tests and Sphinx `-W` exited 0.

**What C1 is.** One group's transport block is a traced memo, `group_transport(basis, walls, sigma_t, resolution: TransportResolution, support)` in `assembly.py`, equal to `LineRule.of(basis, walls, sigma_t, resolution.line_points).transport(support)`; `GalerkinSystem.groups` calls it once per group, so each block is a child entry of the solve, and a point reading pins the same block keys. The test-architect's prototype and this plan's earlier text called it `transport_block`; it was renamed `group_transport` (elegance finding 5: the function returns a `GroupTransport`, and `transport_block(...).block` read twice).

**Content identity.** `RigidMotion`, `ConcentricPartition` and `PanelBasis` are `ContentIdentity` values. `RigidMotion` copies its parts, makes them read-only, refuses a non-finite part, and folds `-0.0` through `numerics.scalars.canonical_reals` (so do the partition's breakpoints); its `_exact_key`, `__eq__` and `__hash__` retired onto the encoder. `PanelBasis.degree` and `TransportResolution.line_points` are normalised by `operator.index`; `Wall` normalises its indices by `operator.index` and its amplitudes by `parse_real` (a late fix after qa's key review).

**The elegance review's changes** (`scratch/characteristic_architecture/p3/c1/elegance_review.md`):
1. The resolution is one argument, `TransportResolution`, not loose `points`, `chunk`, `budget`; its guard `line_points >= 1` is no longer bypassed.
2. `chunk` and `budget` are not memo parameters: the entry pins `LineRule.of`'s defaults (512, 1024) through that def's source. `[M]` qa: editing the default chunk to 256 made an entry `Stale`, and restoring it made it a `Hit` again.
3. The canonical forms `_total_cross_sections` (now `canonical_reals`) and `_region_mask` are shared with `LineRule.of` and `LineRule.transport`: one function at each site.
4. `support` is required on `group_transport` and on `LineRule.transport`: the old default `sigma_t > 0` was a second definition of where emission lives, different from the system's `emission_support() | source_regions`.
5. `canonical_reals` replaced the hand-written `+ 0.0` folds.

**The measurements** `[M]` 2026-10-14 (gate_spec §2.1, §7), on the prototype before the re-pose:
- a miss 2.8 to 3.6 s (5 of 5 cheap blocks), a hit 0.012 to 0.017 s;
- the 431 not-slow rows of the system, reading and assembly files: no memo 242 s, cold 1303 s (383 block entries, 30 MB), warm 102 s; 0 of 431 rows change behaviour.
- qa's key review (`qa_review.md` part two): every argument moves the key (6 of 6), the canonical forms are idempotent and lossless; its two over-keys (`-0.0` in sigma_t, a wall's int amplitudes) were closed by the late fixes. `[M]` the archivist's probe on the landed tree: both spellings now key alike, with a one-ulp sigma_t control keying apart.

**Rulings** (the user, 2026-10-15): land without a size threshold; cut the miss cost first, #602 (about 2.3 s of each cheap miss is the recorder and the manifest).

**The battery consequence.** A block read through the memo is generated from the source on disk, so an in-process mutation (monkeypatch, `-p` rebind) does not reach it. Post-C1 batteries and `catches` re-drops on rows that read blocks mutate textually (a copied tree or a worktree) or enter `bypass()` at `pytest_configure`. Written into `vv-principles` #17 (source `docs/development/skills/vv-principles.md`), the qa and test-architect definitions, and the evidence page `docs/development/evidence/vv-anti-patterns.md` ("2026-10-15 the memo-served value").

**Docs** (the archivist): `docs/theory/references/characteristic.rst`, new section `characteristic-block-memo` ("A transport block is an entry of its own"), the support prose, the door's spy, a Key fact, two Gotchas and two History rows (2b, 2c); `docs/theory/verification/reference_cache.rst`, the clients, the miss price, the bypass's battery paragraph, the declared limit and a History row; `docs/theory/foundations/structured_geometry.rst`, the encoder's population.

⏸ **Next:** C2 (the cylinder point rows, 60 reads at 10 to 20 s each `[R]`) and C3 (a bypass scoped per function; the gates behind `bypass()` hold 2726 s), then the re-timing of the hot characteristic files; #602 before the re-timing.

### Step 2c, C2 and C3 measured (2026-10-15)

Two explorer censuses on `main` @ `da74b1ce`. The files are `scratch/characteristic_architecture/p3/c2/census.md` and `scratch/characteristic_architecture/p3/c3/census.md`. Both ran serially with `-O -p no:cacheprovider`, `OMP_NUM_THREADS=1`, on a warm store, and both report the second of two passes.

**The bypass population** `[M]` (AST over `git ls-files tests orpheus`, 1119 files; both call spellings found as positive controls).
- 50 `bypass()` call sites in 16 files, all under `tests/`. 171 tests enter it.
- **Arm P**, a test that monkeypatches: 20 tests, 25.1 s wall, of which 1.6 s is inside memoised calls.
- **Arm N**, a test that patches nothing: 151 tests, 2226 s wall, of which 2049 s is inside memoised calls. By file: convergence 688 s, albedos 746 s, lift 391 s, independent_references 110 s, reference 71 s.
- Four arm-N rows must keep `bypass()`, because their in-process side is the claim: B3, B6b, M4.12 and M4.5c.
- No code comment, commit message or page states why the readers bypass. `[R]` the reason is the cold cost: a miss regenerates in a fresh interpreter, and C1 measured cold 1303 s against no-memo 242 s on its 431 rows.

**C3, a bypass scoped per function.** `[REFUTED 2026-10-15]` as a speed remedy.
- Its ceiling is 1.4 s per pass by the dynamic predicate, and 0.7 s by the store's manifests.
- The reason: only arm P has a patched def to scope by, and its patched defs (`MeasureCoordinate.derivative`, `group_transport`, `GalerkinSystem.pencil`) are reached by nearly every entry. `group_transport` pins `MeasureCoordinate.derivative` in 278 of 278 entries.
- The plan's 2726 s, at step 2b's boundary ranking item 3, was the time `bypass()` keeps unserved (2049 s today), not what scoping would recover.
- **The fact it establishes, for correctness:** a scoped bypass would serve a stale value for a patch on any of the following, because the recorder never pins them:
  - a module constant or registry entry;
  - a class attribute holding a value;
  - a dataclass-generated method (its code is from `<string>`);
  - a C or third-party function (only the version is pinned);
  - a def read but not called;
  - an attribute added and probed by `hasattr`.
- If C3 is ever built: `bypass(of=(original_def, ...))` refuses anything that is not a first-party Python def, and decides by reach through child entries.

**C2, the cylinder point rows** `[M]`.
- The count was right and the time was not: 60 rows, 239 s per pass in all, median 2.65 s, max 26.6 s. 8 of 60 take 10 s or more.
- The plan's "10 to 20 s each" was a residual (test time minus block time) that also held in-process blocks and test-side mpmath.
- The cost is `PointRule.row` (`reading.py:239`). 86 % of it is `AttenuationRule.of` inside the exit integral, and there is no mpmath.
- 40 of the 60 rows, 177 s, run under an arm-N `bypass()`, together with 453 s of their blocks computed in-process.
- `point_flux` computes every group's row, and a door read of `PointValue(x, g)` keeps one: 11 of 34 door rows are discarded.
- **Candidates:**
  - **C2-0:** a door `PointValue(x, g)` reads only group g's row. No key change; about 16 s per pass `[R]`.
  - **C2-A:** a `point_row` memo reading `group_transport` as a child, keyed on `orbit_coordinate` through `canonical_real`. It saves about 62 s warm, or 239 s once arm N is served. Its hazard: under `bypass()` its generation re-assembles the block.
  - **C2-B:** the same memo keyed on the block's content. It needs `GroupTransport`, `WallCoupling` and `DiffuseWalls` as content values.
- **Rejected,** with the reasons in the census:
  - `PointRule.of`;
  - `point_flux` (the system is `eq=False`);
  - the per-chunk kernel;
  - a batch over positions;
  - the `evaluate` memo, which is already a memo but is bypassed.

**The store** still holds 261 entries under the retired name `transport_block`. Nothing reads them.

**The user's ruling (2026-10-15):**
1. C3 is not built.
2. #602 comes first: cut the miss overhead so that a cold run costs about what no memo costs.
3. Then remove `bypass()` from the 147 incidental arm-N readers. The four rows B3, B6b, M4.12 and M4.5c keep it.
4. Measure cold, warm and no-memo before merging.
5. C2-0 lands alongside, since it is independent.

⏸ **Next:** #602 (the miss cost), then C2-0, then dropping `bypass()` from arm N, then the re-timing.

### #602, the miss cost: rounds one and two [LANDED 1d8e2c99] (2026-10-15 to 16)

Cheap miss (5 sphere blocks × 3 repeats, `-O`, serial; `scratch/characteristic_architecture/p3/m602/`): 2.65–2.88 s before, 1.14–1.35 s after round one, 0.66–0.85 s after round two. Three files (system, reading and assembly, 431 not-slow rows, a fresh worktree each): cold 1303 → 538 → 406 s; no memo 236 s; warm 168 s. The warm figure is above C1's 102 s because the landed assembly gates build their blocks in process.

**What landed:**
- **Round one** (`8e02236f`, `e666e58d`): one parse per file, one version lookup per checkout, the child exits with `os._exit`.
- **Two pre-existing defects fixed:**
  - ERR-115, the stale timestamp `.pyc`: a regenerated entry served the old answer and validated. The fix is a pinned-source loader in the generating process.
  - A docstring-only def crashed `_def_pin`.
- **Round two:**
  - D (`0c569df8`): per-file ownership from `RECORD`. It fixes a third defect, namespace tops pinned to the wrong distribution.
  - B (`155fbdbd`): an on-disk source index. Only a generating process writes one, after qa's case of an index written under another code's name. The memo pin names the code that ran (the elegance review's F4).
- **Gates:** `test_traced_memo_compiled_bytes.py`, `_ownership.py`, `_source_index.py`. Every battery arm reddens its row (batteries 0, D, B, R).
- **Full suite:** 17 136 passed.
- **Filed:** #603, nothing prunes the three stores.

**E, forking each miss from a template:** gated on the fork-safety protocol (`round2_gate_spec.md` §3.3), which is running.
- `[M]` A fork after a BLAS call crashes the child (8 OS threads under Accelerate).
- A template that has only imported was safe in 30 of 30 forks.

⏸ **Next:** E's fork-safety verdict. Then, by the user's ruling, C2-0 and dropping `bypass()` from the 147 arm-N readers, then the re-timing.

### #602, cut E ready to land (2026-10-16)

E, forking each miss from a recording template, is written and reviewed on `perf/traced-memo-fork-template` (uncommitted in the main tree at the time of writing; `git log` decides whether it has landed). The chapter of record is `docs/theory/verification/reference_cache.rst`, section "The fork template (#602, cut E)".
- **Fork-safety protocol PASSED** (`[M]` 2026-10-16, `scratch/characteristic_architecture/p3/m602/ta/fork/full/driver.log`): 9 cells (`OMP_NUM_THREADS` × `VECLIB_MAXIMUM_THREADS`, each unset, 1, 8), 3600 forks per cell, 32 400 of 32 400 answered and equal to a fresh interpreter's digest, no new crash report. Control: a template that called BLAS before forking lost 20 of 20 children to a signal.
- **Scope:** only darwin/arm64/Accelerate forks (`FORK_MEASURED`); every other platform spawns fresh.
- **Review fixes:** every input (file bytes, presence, listing) observed at its first touch (ERR-116, pre-existing, also on the fresh route); the run's own output not observed; a run that changes an input it observed refused as `Unpinnable` (qa's F1, F2; N2's absent-then-created cache is refused); a device never read (N3); the environment pin is what the generation was handed (PEP 538 under `LANG=C`); a distribution pinned by its `RECORD` digest; a lock on the template's pipes; `sys.argv` emptied and stdin `/dev/null` in the child; a fresh template's refusal turns forking off for the process.
- **Cost** (`[M]` 2026-10-16, `e_miss_cost.py`, `-O`, five processes per route): a cheap miss 0.293–0.301 s forked against 0.650–0.663 s fresh; the first miss of a process 0.81–1.01 s (it starts the template).
- **Full suite:** 17 181 passed, 0 failed (`-m "not slow"`, `-O`, 53 min).
- **Open, untraced:** with the pinned-source loader removed (arm P-a), E3's def-body case stays green while arm E-c reddens it, so a second mechanism also refuses that stale template (probably the def's source-text pin); not a defect, not traced.

⏸ **Next (the user's ruling):** land E; then C2-0; then drop `bypass()` from the 147 arm-N readers and measure the three files cold, warm and with no memo; then the re-timing.

### C2-0 landed; the bypass drop measured on the not-slow rows and parked (2026-10-16)

**Landed on `main` @ `bb2835d4`** (full suite `-m "not slow"`, `-O`: 17 182 passed, 0 failed, 55 min):
- E (`98dc5242`, docs `ab8217ee`): the fork template, the first-touch inputs, ERR-116. CI run 38934685877 green.
- C2-0 (`410bfca8`): `GalerkinSystem.group_point_flux(position, emission, g)` transports one group; `point_flux` is the stack of it, and the characteristic reference's point read (`reference.py`) calls the one-group primitive. Gate: `test_c2_0_a_point_value_transports_its_own_group_alone`.
- E7c (`452dedb8`): the killed child ends with `os._exit(70)`, not `os.abort()`, which left a macOS crash report per run.
- ERR-034's first red re-aimed (`bb2835d4`): the documented arm edited `chart.AxialImage.orbit_coordinate_at`, which the slab readings never reach (0 of 4 rows red); the re-drop in `chord.Rays.orbit_coordinate_at` reddens 3 of 4 (both mirror cases and the flux row's 1G case; PU2 stays green).

**The bypass drop (arm N, 27 sites in 10 test files)** is committed on `test/arm-n-bypass-drop` (rebased onto `bb2835d4`), UNMERGED. `[M]` 2026-10-16, `scratch/characteristic_architecture/p3/c3/drop/timing/summary.txt`: the 10 files, not slow, a fresh worktree each, serial, `-O`, `OMP_NUM_THREADS=1`, 375 passed and 3 skipped in every run.

| tree | cold | warm | no memo | store entries after cold |
|---|---|---|---|---|
| before the drop (`410bfca8`) | 518 s | 484 s | 505 s | 68 |
| after the drop | 1474 s | 347 s | 525 s | 1819 |

A cold run pays one process per miss (1751 new entries); a warm run saves 137 s. The 40 slow rows of the same files hold most of the census's 2049 s and were not in this measurement.

**The user's ruling (2026-10-16), "measure slow, then decide":** keep the drop unmerged; time the 40 slow rows cold, warm and with no memo before and after; decide per population with both numbers in hand. The driver is `scratch/characteristic_architecture/p3/c3/drop/timing/run_slow.sh` (summary `summary_slow.txt`; worktrees `wt_slow_before` at `bb2835d4`, `wt_slow_after` on the branch).

**Owed by the drop when it lands:** the 22 `verifies`-only rows that reach a dropped site are re-dropped one by one (an in-process mutation no longer reaches a memo-served value, `vv-principles` #17).

⏸ **Next:** the slow-row verdict to the user; then the drop per population; then the re-timing of the three hot files (#602's last done-when; #602 stays open).

### The drop landed everywhere; the warm slow run profiled; gmpy2; a template refusal made audible (2026-10-17)

**The user's ruling (2026-10-16):** drop `bypass()` from arm N everywhere, and time gmpy2 in parallel.

**The drop is on `main`:**
- `6af57441` and `3d3b0f75`; CI 39005451105 green.
- `3d3b0f75` removes `verifies("characteristic-adjoint")` from A1 and A2. These are the 2 markers of 26 that the re-drop could not witness; the other 24 reddened under a textual mutation on a worktree copy (`scratch/characteristic_architecture/p3/c3/drop/redrop/verifies_verdicts.md`).

**The 40 slow rows of the 10 files** `[M]` 2026-10-16:
- Method: a fresh worktree each, serial, `-O`, `OMP_NUM_THREADS=1`, 40 passed in every run.
- Source: `drop/timing/summary_slow.txt`.

| tree | cold | warm | no memo | store entries |
|---|---|---|---|---|
| before the drop (`bb2835d4`) | 3212 s | 2824 s | 3107 s | 14 |
| after the drop | 3272 s | 905 s | 3208 s | 284 |

A cold run costs +60 s; a warm run saves 1919 s.

**The warm slow run, profiled** `[M]` (`drop/timing/profile/warm_slow.prof`; 942 s of wall time, 1290 s under cProfile):
- 82 % of profiled time is `mpmath.quad`, all of it called from the test oracles `_characteristic_point_mp.phi_cylinder_vacuum` and `phi_cylinder`. These oracles are not memoised, so a warm run still pays for them.
- About 85 s is in-process `point_flux` reached through `_given`.
- About 57 s is the PS1982 reference.
- About 10 s is memo lookups.

Three concepts from the profile, all `[HYPOTHESIS]`, to be tested in the order of the time they hold:
1. Store the mpmath oracles as references, i.e. memoise the oracle. This would turn 82 % of warm time into lookups.
2. Route the in-process `point_flux` reads through the memo.
3. Memoise the PS1982 reference.

**gmpy2 (`3ff7866b`), adopted by the user's ruling:**
- The floats are identical on three routes `[M]`, single runs:

  | route | pure Python | gmpy2 |
  |---|---|---|
  | F1 | 11.1 s | 7.8 s |
  | F4 | 26.7 s | 7.7 s |
  | F2 | 241 s | 30 s |

- sympy's ground types become gmpy too.
- The full not-slow suite on the branch matched `main` exactly: 17 182 passed, 261 skipped, 0 failed, 69 min, cold store.

**A refusal by a fresh template is warned (`8f11ad66`, #602):**
- The drop branch's first full run had 3 extra skips (17 179 / 264). A re-run of the same commit matched `main` (17 182 / 261) and logged no refusal, so the cause was not recovered.
- The one mechanism that produces exactly those skips is a fresh template's refusal. It turned forking off for the process with no output.
- It now raises `TemplateRefused` once, with the template's reasons. S4 asserts it, and the arm `0 and warnings.warn` reddens S4.
- The gate run above printed 0 `TemplateRefused`.

**Merged @ `3ff7866b`** (ff), pushed.

⏸ **Next:**
1. The three-file re-timing (system, reading, assembly; not slow; cold, warm, no memo) at `3ff7866b`. This is #602's last done-when. Driver: `scratch/characteristic_architecture/p3/c3/retime/run.sh`; compare `p3/c1/gate_spec.md` §7 (242 / 1303 / 102 s on 2026-10-14).
2. Record the re-timing on the reference-cache theory page.
3. Concepts 1 to 3 above, oracle first.

## The hit cost, then concept 1: the cylinder mpmath routes stored (2026-10-17)

**The hit cost (#602, closed).** A warm hit had come to cost what no memo did: 222 s warm against 221 s with no memo, at `3ff7866b`. Each validation re-read every pinned source and re-found every pinned distribution. `8fbb1d61` keeps each read for the process, guarded by the file's stat identity (`FileIdentity`, `Witness`, `Kept`), and it is the first stage of ruling 15 in `verification_tree_refinement.md`. `[M]` the three files, single runs: no memo 221 s, cold 820 s, warm 163 s. Merged @ `b2037667`.

**Concept 1, landed.** The two costly cylinder routes, `phi_cylinder_vacuum` and `phi_cylinder` in `tests/gates/derivations/_characteristic_point_mp.py`, are `@traced_memo`, with no `lru_cache` above them. The routes that take milliseconds are not stored.
- **A defect it exposed, fixed** (`7ca91c43`). A memo under a `functools.wraps`-style decorator refused to generate: the child found the outer wrapper by name and called the memo again. The child now resolves the memo through `_resolve_memo`. Gates M3.1b and M3.1c.
- **The independence gate is an AST census.** It covers the oracle and every test module it imports, and the one `orpheus` import it allows is the memo.
  - [REFUTED 2026-10-17] the first version compared the module sets that fresh interpreters hold. It was blind on two counts: the memo alone loads 76 orpheus modules, all inside the reference's closure, and an import inside a function body is never run (qa and elegance-enforcer, `scratch/characteristic_architecture/p3/c3/oracle/{qa,elegance}/review.md`).
- **`[M]` timings,** the slow rows of the lift and reading files (23 rows), single serial runs, a fresh worktree each:

  | tree | cold | warm | entries |
  |---|---|---|---|
  | before (`b2037667`) | 1075 s | 353 s | 17 |
  | after | 1044 s | 127 s | 34 |

- **What remains of the warm 127 s** (`oracle/after_warm.log`): in-process reads through `point_flux`, about 28.9 s and 25.6 s on the top two rows; the lift's polar-column quadrature, 14.3 s.
- **Filed:** #605, #606, #607.
- Merged @ `9f924a2d` (`7ca91c43`, `333268fb`, the docs).

⏸ **Next:**
1. Concept 2: the in-process `point_flux` reads through the memo. First measure which rows hold the remaining warm time and whether their reads are a function of content (the key).
2. Concept 3: PS1982 (about 57 s, scipy, returns a frozen dataclass).
3. Close-out asks, held for the user:
   - the worktrees under `p3/c3/{retime,hitcost,hitcost2,oracle}` and the older ones listed in the 2026-10-16 note;
   - the branch `refactor/chord-oracle-axial-lift`.

## Concept 2, and the two more oracle routes (2026-10-17)

**Two more oracle routes stored** (`76d0409d`, the user's ruling "Yes, store them"). The partial-mirror total `specular_cylinder_total` is a stored route. The lift's inline polar density became the stored column `cylinder_polar_density`, one entry per column: a node costs less than a miss, and a column costs 9 to 27 s. `[M]` the 27 rows that read them: cold 162 s, warm 29 s. Merged @ `712293a9` with the docs.

**Concept 2, landed** (the user's ruling "Door, with batteries in bypass"). Eight slow cylinder point-reading rows now read their flux through the door, `CharacteristicDerivation.evaluate(PointValue)`, via `_given_by_the_door`: C2-cylinder ×3, Q2, C2b's wall leg ×2, Q4 `[wall]`, and lift F4.
- `[M]` the test-architect's spec, `p3/c3/concept2/spec.md`:
  - 13 of 13 group readings are within 1 ulp of the in-process reading;
  - under `bypass()`, every re-posed row reddens under at least one arm;
  - the door witness reddens exactly the 10 door rows.
- `[M]` both files on the main tree's store: 240 passed cold for the door keys, in 1 882 s, and warm in 217 s. The 23 slow rows warm take 25.8 s, each door row 0.03 to 0.13 s.
- C6 stays in process: a basis function is not a `RegionwiseConstant` question.
- Battery plugins:
  - b5b enters `bypass()` at configure;
  - the product plugin already pushed `InProcess`, and its ERR-105 run must drop `-m "not slow"`;
  - 11 of b5b's 26 arms no longer install.
- Filed:
  - #608: C2b's documented first red is blind on today's tree.
  - #609: one posing for `_given` and `_given_by_the_door`.
- Merged @ `893d85b5`.

**Concept 3, landed** (the user's ruling "Memo the solver", over splitting out the optical-radius eigenpair).
- `solve_ps1982_vacuum_sphere` carries `@traced_memo`, and the helper's `lru_cache` is removed.
- The census gate is now `test_each_stored_reference_imports_nothing_but_the_memo`, with the rows `[mpmath-routes]` and `[ps1982]`.
- `[M]` The PS-1982 file runs 87 s cold and 2.3 s warm.
- `[M]` The 40 slow rows of the 10 files take 27.1 s warm, against 82.9 s after concept 2 and 2824 s before the drop. 20.7 s of the 27.1 s is C6, in process.
- [REFUTED 2026-10-17] that a docstring edit regenerates the entries. The claim came from qa's 227 s run and was repeated in `da043d1c`'s message. `[M]` A docstring-only edit to the function left the file warm (3.15 s), and so did a module-docstring edit (2.50 s), because the pins strip docstrings (gate M2.1). The 227 s run's cause was not recovered (`concept3/probe_docedit.log`, `probe_moddoc.log`).
- Merged @ `dcf631dd`.

**Concepts 1 to 3 are closed.** What remains in the warm slow run of those 10 files is C6's in-process `PointRule.row`, about 21 s. The assembly file was never among the 10 (see the re-timing below).

⏸ **Next:**
1. C6 stays in process, by the user's ruling (2026-10-17). It is a self-consistency theorem: the code is compared with itself. Storing its 648 point rows would leave a warm run comparing the live block with a stored copy of the code's own output. Its cost goes to ruling 15's verdict revalidation instead (`verification_tree_refinement.md`, "One change-detection layer").
2. The follow-ups: #608, #609, #605, #606, #607.
3. Close-out asks, held for the user. Both were done 2026-10-17: 18 scratch worktrees removed, the dirty ones' diffs saved beside them as `<wt>.uncommitted.patch` and `<wt>.untracked/`; the branch `refactor/chord-oracle-axial-lift` (`a336bde4`, local only) deleted.

## The 24 files re-timed after concepts 1 to 3 (2026-10-17)

`[M]` `scratch/characteristic_architecture/p3/close/timing/` (`run.sh` copied from the 2026-10-14 driver; `results.tsv`; `logs/` with `--durations=0`). The same 24 files, serial, `-O`, `OMP_NUM_THREADS=1`, detached worktree at `a456a9bf`, 17:21 to 23:07. 96 of 96 rows pass; 3 warm2 rows exit 5 (every test in the file is `slow`).

| pass | 2026-10-14 | now |
|---|---|---|
| bypass | 11 601 s | 8 630 s |
| cold | 10 329 s | 8 571 s |
| warm1 | 9 307 s | 3 073 s |
| warm2 (`not slow`) | 991 s | 487 s |

| file, warm1 | 2026-10-14 | now |
|---|---|---|
| characteristic_assembly | 3 535 | 2 579 |
| characteristic_system | 1 105 | 4 |
| characteristic_lift | 1 090 | 134 |
| characteristic_albedos | 988 | 6 |
| characteristic_reading | 1 032 | 88 |
| characteristic_convergence | 933 | 3 |
| characteristic_transport | 166 | 138 |
| characteristic_independent_references | 213 | 4 |
| characteristic_reference | 107 | 15 |

**What it says.**
- The stored references did what they were for. Outside the assembly file, warm1 fell from 5 772 s to 494 s.
- The bypass and cold falls (about 26 % and 17 %) are `[R]` the reference's own speed-ups since 2026-10-14 (the axial lift, the derived node counts), not measured per cause.
- **81 % of what remains warm (2 493 of 3 073 s) is one file's slow tier.** The assembly file's 21 `slow` rows take 2 493 s warm (warm1 2 579 minus warm2 86). 12 cylinder rows carry 2 345 s:
  - the 6 closed-body conservation rows on the three-region cylinders (AS1, `_CLOSED_SLOW`), 277 to 418 s each;
  - the 5 cylinder escape and transmission closed forms, 30 to 97 s;
  - the equal-materials interface on a cylinder, 109 s.
- These rows build the block in process (`LineRule.of(...).transport`, not the stored `group_transport`): the block is their system under test. Like C6, no stored reference reaches them. The two levers are a cheaper cylinder block (#587, the non-tensor (b, θ) rule) and ruling 15's verdict revalidation.
- The remaining warm cost of lift, transport and reading (88 to 138 s each) is also in-process transport, and most of it is outside `slow`.

**#591** (the point reading's cost): its done-when protocol ran at P3 step 3f (`[M]` route B 0.146 s against the old hoist's 0.120 s, within 22 %; the solve 4 times faster). The second half of the done-when, a ruling on the remaining gap, is the user's. **#587**: open, and now the named lever for the remaining 2 345 s.

⏸ **Next:** the user's rulings on #591 (close at 22 %, or keep for the emission contraction) and on #587 (take the non-tensor rule now, or leave the assembly slow tier to ruling 15). Then the verification-tree campaign's second pass (ruling 14).
