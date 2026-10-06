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
