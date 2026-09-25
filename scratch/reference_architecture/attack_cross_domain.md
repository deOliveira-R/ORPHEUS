# Cross-domain attack on `.claude/plans/reference_cache.md` "The architecture as it stands" (W5, #405 step 2)

cross-domain-attacker, 2026-09-24, HEAD 901f64ca. Nexus MCP tools present but not needed (grep/AST/Python probes sufficed).

## PREMISE (step 0, measured before any frame)

- `[M]` Sood LA-13511 (`LA13511_CASES`, 47 cases, runtime dump): `critical_dimension_mfp` set on **25/47**, `k_eff_or_kinf` 47, `scattering_order=1` on 5, geometry `infinite` 22. So for 25 published problems the printed answer is a DIMENSION: the plan lists `La13511Case.to_geometry()` reading the truth as a defect ("posing reads the answer", plan l.93); it is the native shape of a criticality search.
- `[M]` `eigenvalue_kind` literals in `orpheus/`: `k_eff` 6, `k_inf` 2, `c_critical` 4, `none` 1 (`common/solution_types.py:127` `CriticalSolution` documents also `sigma_a`). The eigen-parameter is not always k.
- `[M]` References to a DISCRETE or APPROXIMATE equation exist and are consumed: `sn_slab_1eg_2rg_S8` is "the exact discrete-S_N eigenvalue at Gauss-Legendre order N, consumed only by SN tests" (`continuous/cases/sn.py:706-716`); `flat_source_cp` supplies 27 `_CASES` (the exact solution of the flat-source CP equation, one cell per region), consumed by `tests/gates/cp/test_verification.py`; the Peierls white-BC references solve a rank-N closure (`peierls_nystrom/geometry.py:4891 reflection_white_rank2`, alias `white_rank2 -> white_f4` at :5480).
- `[M]` Operator references consumed by production gates: `discrete/sn/face_transmission.py` ("algebra of record") by `tests/gates/transport/spatial/test_face_transmission_{symbolic,xverif,damping}.py` and `tests/gates/sn/mesh/test_reflective_axis_pairs.py`.
- `[M]` Tail-decay bound on a physical slab flux shape, f(x) = 2 - E2(x) - E2(a-x), a = 2 (the uncollided-like shape; x log x at both free surfaces), single Chebyshev panel, estimate = max |last 4 coefficients|:

  | n | tail estimate | true sup error | true/estimate |
  |---|---|---|---|
  | 16 | 1.8e-3 | 7.2e-3 | 4 |
  | 32 | 9.6e-5 | 1.8e-3 | 19 |
  | 64 | 5.7e-6 | 4.5e-4 | 78 |
  | 128 | 3.6e-7 | 1.1e-4 | 315 |

  Convergence is O(n^-2), and the estimate UNDER-states the error by a factor growing like n. Control (e^-x): estimate at or above round-off, fine. With geometrically graded panels (grading 0.15, 8 levels per end, 16 points per panel, 288 coefficients) the true error is 1.6e-8 `[M]` (inline probe in this session, reproduced in the NEEDS note).
- `[M]` `Mixture` (`orpheus/data/macro_xs/mixture.py:137-145`) carries `SigS` per Legendre order and no velocity field (grep `velocit|inv_v` in the file: 0). No alpha-eigenvalue reference exists in the tree (the `alpha` hits are Variant alpha).

What it collapses the question to: the design fails on four STRUCTURAL assumptions the tree already violates (question derived from source presence; one solution answers one specification; the field bound read from tail decay; the subject of a reference is a problem), and holds on the physics-data axes (anisotropy, groups, source shape).

## Holes, ranked by rewrite risk

| # | hole | concrete case | additive or rewrite | the change now |
|---|---|---|---|---|
| 1 | The question is derived from source presence (2 values). | 25/47 Sood cases answer a critical dimension; 4 `c_critical` sites; an incident-flux boundary (Milne, albedo, Case half-space) is source-driven with NO volumetric `Source`, so the rule calls it an eigenvalue problem; an adjoint (importance) question has no spelling. | **Rewrite** if left: every specification, registry key and the comparison verb key on the 2-valued derivation, and a critical-size problem is a FAMILY of specifications, not one. | A typed `question` field (a sum type) whose DEFAULT is the ruled derivation: `Eigen(parameter)` with parameter in {k, c, ...}; `FixedSource`; reserve `CriticalParameter(breakpoint, condition)`; `adjoint` as a flag on any member. And the comparison verb must take a declared bound transfer (for a critical-size published solution the SUT poses at r* and checks k = 1 within abs(dk/dr) times the printed delta r; the design's "bound on the same observable" cannot compare a k against a dimension bound). |
| 2 | One solution answers one `Specification`; nothing names the EQUATION solved. | The S8-exact Case reference and any continuous reference (F_N, trajectory resolvent) of the same specification differ by the S8 angular error, which is far larger than either bound, so the ruled corroboration rule ("failure makes both `Invalid`") trips two correct references `[R]`. The 27 flat-source CP cases are `Exact` for the flat-source equation, false as `Exact` for the specification. Rank-N white closures answer an approximation of the white law. | **Rewrite** after P4: the registry and certificate key on the specification alone. | Each solution declares `equation = (spec, approximations)`, an ordered tuple of declared approximations (`AngularQuadrature(set, N)`, `ClosureRank(n)`, `FlatSource(partition)`), empty for the continuous problem. Corroboration only between equal equations; the verb refuses a SUT whose declared approximations do not match, unless an arrow with a certified approximation bound is declared. |
| 3 | The field representation and its bound. | (a) `[M]` table above: a free-surface slab flux converges O(n^-2); the tail estimate is 4 to 315 times too small, so the certificate floor ("bound at most a tenth of the tolerance") passes on a false bound. (b) Curvilinear psi(r, mu) with an inner interface r_i has a singular curve at the grazing ray mu = sqrt(1 - (r_i/r)^2) `[R]`, not at mu = 0; a tensor Chebyshev split at mu = 0 does not resolve it. (c) 2-D/3-D angular flux: a 5-D tensor Chebyshev at 32 per axis is 3.4e7 coefficients, about 270 MB per group per region `[R]`, and the exact 2-D psi is discontinuous on rays from every interface corner. | Additive IF consumers never read the coefficients; **rewrite** if tests call a field accessor. | The protocol is `evaluate(functional) -> value with bound`, the representation private. Piecewise Chebyshev on PANELS (breakpoints a superset of the region boundaries, geometrically graded toward free surfaces and interfaces). The bound admissible only when a fitted tail shows geometric decay; otherwise refuse, never report. For 2-D/3-D, a lazy evaluator whose results are cached per functional key. |
| 4 | The subject of a certificate is always a problem. | The face-transmission algebra (4 gate files), P_ij escape probabilities, K_vol, E_n kernels are references for OPERATORS (a morphism parameterised by optical thickness), with no materials-geometry-question. | Additive if the certificate and the verb are generic in their subject. | `ReferenceCertificate[Subject]`, `Specification` one subject type among others; do not hard-type the registry key to `Specification`. |
| 5 | Symmetry quotients. | A reflective half-slab and the full vacuum slab, or `from_homogeneous(width, reflective or periodic)` and the infinite medium: different content identity, same physics; a corroboration or reuse across them is not spellable. | Additive. | Observables stated in coordinates, so they can be pulled back along a declared covering morphism later; name `from_homogeneous` as the first such morphism, not a special case. |

## STRUCTURAL FEATURES

1. The specification is a point in a parameter space (dimension, c, densities) with a symmetry group descending along the overlay chain.
2. The question is a posing over one operator: an eigenvalue of a pencil, a root in a geometry parameter, a fixed source, or the dual of any of these.
3. Solutions sit at different levels of approximation of one continuous problem (closure rank, quadrature, spatial basis).
4. The answer's field has endpoint log singularities (free surfaces, interfaces) and characteristic discontinuity curves.
5. Some references are morphisms, not objects.
6. Bounds propagate along functionals (field sup norm to any cell average).

## ELEGANCE DETECTOR HITS

- A repeated conditional is a missing type: "source absent means eigenvalue" is a tag computed from another field's presence; `eigenvalue_kind` string literals are the same missing type on the solution side.
- A field that must be read to pose the problem (`to_geometry` reads `truth`) marks a missing parameterised object (a family), not a layering defect.
- The word "applicability" carries what should be a typed relation (which equation a solution answers).

## FRAME CANDIDATES

**Fibration of specifications over a parameter space (A.1), with the pencil-posing reading (D2 of the backbone).**
Trigger: features 1 and 2; 25/47 Sood answers are parameters. Reformulation: a `SpecificationFamily` p -> Specification(p) (a geometry breakpoint or a material scalar left free); an eigen-question is a pencil posing A psi = lambda M psi with M in {F (k), S+F (c), 1/v (alpha)}; a critical search is a SECTION: the p* where lambda(p*) = 1. The published r* with delta r restricts to the fibre at r*, and the SUT's k bound there is abs(dlambda/dp) delta r. Payoff: structure-exposing (posing the problem from the answer is a restriction, not a defect); expressive (k, c, alpha, dimension, adjoint as one question type). First test: pose Sood PUa-1-0-SL at its printed critical half-thickness and require the verb to REFUSE to compare k without a declared sensitivity; an implementation that compares k against the dimension's printed precision passes a wrong design and fails this test.

**Poset of approximations (a category of equations, with an error functor on its arrows).**
Trigger: feature 3. Reformulation: equation = spec plus an ordered chain of approximations; arrows carry an approximation bound when one is certified; corroboration and verification are defined within a node or along a bounded arrow. Payoff: structure-exposing (discrete-exact references become first-class, and are consistent with the ban, which forbids fine-mesh runs, not exact solutions of a discrete equation); structurally simpler ("applicability" prose becomes one typed tuple). First test: register the S8 Case reference and an F_N reference for the same 2-region slab; the corroboration must NOT trip. The ruled design trips it.

**Characteristic (impact-parameter) chart for curvilinear psi (A.1, Noether).**
Trigger: feature 4. The impact parameter p = r sqrt(1 - mu^2) is conserved along a ray, so the singular set is the lines p = r_i. Store psi on (r, p) panels split at each r_i. First test: hollow-sphere psi, panels in (r, mu) against panels in (r, p) at equal coefficient count; the (r, p) sup error must be geometric in n, the (r, mu) one algebraic.

**Adaptive piecewise-polynomial approximation with geometric grading (hp, chebfun-style splitting).**
Trigger: feature 4 and the table. Payoff: algorithmic (1e-4 at 128 coefficients to 1.6e-8 at 288, measured); the bound becomes honest. First test: the table itself; an estimator that passes the single-panel E2 flux fails it.

## CROSS-METHOD POLLINATION

From CP (escape probabilities as operator references) into the SN verification of `face_transmission`: both are references parameterised by optical thickness, not by a problem. One `Subject` type serves both.
From MC practice (tally functionals): an observable is a (response function, region) pair, so an adjoint or response reference is a functional of the forward field, not a new problem.

## UNEXPLORED (checked, not triggered)

- Sheaf of certificates: the bound restricts from the field to every functional by the sup norm, which the design already states; no new object. Question: is the certificate a sheaf?
- Continuous energy and ICSBEP physical description: no analytical reference exists in continuous energy. The ANL-7416 "physical situation" record is additive above the specification. Question: rewrite risk. Constraint: `ExperimentalResult` must attach to that physical record, not to the multigroup `Specification`, or model error is folded into the comparison.
- Anisotropic scattering: `Mixture.SigS` per Legendre order `[M]`; the 5 P1 Sood cases are expressible. Question: does the spec express anisotropy?
- Anisotropic or surface sources: a new `Source` member, or a boundary law on the geometry (ruled). Additive, except for the question derivation, which is hole 1.
- Alpha and time: no reference in the tree; they need 1/v (`Mixture` has none `[M]`) and a question member. Additive once hole 1's sum type exists.
- Tensor networks for the 2-D/3-D payload (tensor-train compression of psi): the discontinuity rays destroy low rank. Refuted for storage. A functional-cache is simpler.

## SELF-CORRECTION

None needed.

NEEDS: user rulings on (1) a stored `question` sum type with the derived default (reopens discussion 5's "never stored"); (2) `equation = (spec, approximations)` on every solution before P2 defines corroboration; (3) the panel representation and the refuse-unless-geometric bound rule before P2 fixes the payload. Re-run of the E2 table: the inline probe is reproduced by f(x) = 2 - expn(2, x) - expn(2, a - x), a = 2, Chebyshev nodes, sup error over 20 001 points.
