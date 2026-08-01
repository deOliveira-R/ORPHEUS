# B3.4a doc repair — working notes

Branch `refactor/operator-strategy-layers`. Landed-in-tree (uncommitted) code:
`orpheus/sn/boundary/angular.py`, `orpheus/sn/boundary/realizer.py`.

## Ground truth read from LIVE code (2026-08-01)

### AngularAverageOperator (angular.py)
- `__init__(cos_w, norm, n_inflow)`. `cos_w` shape `(|Γ₊|,)`, STRICTLY positive
  (`if (cos_w <= 0).any(): raise`). Old guard was `(cos_w < 0).any()`.
- Attributes: `n_outflow` (= cos_w.size, DOMAIN) and `n_inflow` (CODOMAIN).
  `n_ordinates` RETIRED.
- `from_quadrature(quadrature, axis, outward_sign)` — signature unchanged.
  Body now: `face = f"{axis}{'max' if outward_sign==+1 else 'min'}"`,
  `omega_dot_n = build_omega_dot_n(quadrature, (face,))[0]`,
  `outflow = flatnonzero(omega_dot_n > +TANGENTIAL_EPS)`,
  `inflow  = flatnonzero(omega_dot_n < -TANGENTIAL_EPS)`.
  Raises if either is empty. `cos_w = weights[outflow] * omega_dot_n[outflow]`.
  Old private test `(outward_sign * mu_n) > 0.0` RETIRED. Old per-axis
  `mu_x/mu_y/mu_z` switch RETIRED (`axis not in AXIS_NAMES` now).
- `apply(psi)` — requires `psi.shape[0] == n_outflow`; returns
  `(n_inflow,) + psi.shape[1:]`.

### IncomingSourceOperator (angular.py)
- `__init__(source, *, n_inflow)`. `inflow_indices`, `n_ordinates`, and the
  internal `_inflow_mask` all RETIRED. Unmasked fallback branch RETIRED.
- `apply(psi_out)` = `source.evaluate((n_inflow,) + psi_out.shape[1:])`.

### realizer.py
- New `_checked_angular_average(law, quadrature, method_space, gamma_out)`:
  rebuilds `Γ₊` from the LAW's own `axis`/`outward_sign` (as a face name),
  compares `np.array_equal(law_outflow, gamma_out.indices)`, raises
  `BoundaryError(..., law="white")` on mismatch. Comparison is on index SETS
  (not sizes) because `|Γ₊| == |Γ₋|` makes a size compare Mode-12 blind.
- white arm: `gamma_out = _outflow_restriction(method_space, "white")` then
  `_checked_angular_average(...) & IdentityOperator()`, albedo fast path.
- prescribed_inflow arm: calls `_outflow_restriction(..., "prescribed_inflow")`
  FOR ITS GUARD (value discarded), requires `inflow_indices`, returns
  `IncomingSourceOperator(law.source, n_inflow=inflow_indices.size)`.
- ⇒ `SNMethodSpace.minimal(quad)` (quadrature-only) can NO LONGER realize
  white or prescribed_inflow (nor vacuum/reflective, since B3.2).

### NOT narrowed (keep honest)
- albedo → ZeroOperator / IdentityOperator / α·(I & I) — full-face endomorphism.
- periodic → PeriodicWrapOperator & IdentityOperator — full-face endomorphism.
  B3.4b / B3.4c, still strict xfails.

## Findings log

(appended as I go)

## MEASURED this session (probes against live tree)

| probe | result |
|---|---|
| `minimal(quad)` realize white | `BoundaryError` "without outflow_indices" |
| `minimal(quad)` realize prescribed (NoSource) | `BoundaryError` "without outflow_indices" |
| `minimal(quad)` realize prescribed (nonzero q) | `BoundarySourceNotOnIncomingTraceError` (ERR-047, fires first) |
| `minimal(quad)` realize albedo/periodic | STILL realizes (un-narrowed) |
| white: Γ₊ probe (4,3) | → (4,3); full-face (8,3) → **REFUSED** ValueError |
| prescribed: Γ₊ probe | → (4,3); full-face → (4,3) — narrowed but does NOT validate |
| albedo0.5 / periodic: full-face | → (8,3) — still ENDOMORPHISMS |
| vacuum / reflective: full-face | → (4,3) silently — still do NOT validate (B3.2 warning stands) |
| `0.3*spec + 0.7*white` via realize_recursively | WORKS, (4,2), matches pointwise sum |
| `0.3*spec + 0.7*albedo(0.5)` | ALSO "works" (4,2) — silently Γ₊→Γ₊, Mode-12 blind |
| SNMesh.BOUNDARY_OPERATOR_REGISTRY | still `{vacuum, reflective}` — the "registry admits only" claim is TRUE |
| tangential @ xmax | GL(8): 0; product(2,4): 4 (>0 vs >eps differ by **2**); lebedev(17): 12 (>0 vs >eps differ by 0); LS(6): 0 |

## OUT OF SCOPE but falsified — FLAG, do not edit
- `orpheus/geometry/boundary/_base.py:326-364` `assert_source_lives_on_incoming_trace`
  docstring: still says the q ∈ Γ₋ guarantee "rests on the realizer masking the
  evaluation" and "the realizer restricts the evaluation to Γ₋ at construction".
  Mask is gone. Also its claim is a PRESENCE check, not an entry-wise validation.
- `tests/geometry/test_bc_equivalence_snapshot.py:101-108` `_MIXED_LAW_XFAIL`
  + `test_boundary.py` twin: reason text is now false (white IS narrowed). The
  row still xfails but for a DIFFERENT reason (`minimal(quad)` now raises) —
  vv Mode-8 class 4, MISATTRIBUTED strict-xfail. Test migration owns it.

## EDITS APPLIED (all gated by `-E -W` EXIT=0, 0 W/E/C, == baseline)

docs/theory/methods/sn/boundary_conditions.rst  — 4 edits (white row, albedo/periodic
  phase tags, prescribed row, the narrowing note)
docs/theory/foundations/boundary_conditions.rst — 14 edits (incl. NEW section
  `bc-narrowing-b34a`)
docs/theory/foundations/operator_algebra.rst    — 1 edit (ERR-047 presence-check bullet)
orpheus/geometry/boundary/white.py              — 1 edit (Γ₊ sum, realization, example, guard)
orpheus/geometry/boundary/prescribed_inflow.py  — 3 edits (R=G=0 ×2, mask, example)

docs/theory/verification/matrix.rst regenerated by the build: tests 7192→7202
(+10) — that is the CONCURRENT test-migration agent's work, not mine. I added no
eq-labels, so no matrix row changed.

## Verified-but-TRUE (looked false on grep, is correct — corpus-fragility signal)
- foundations:1306 "At the time of the extraction the realized per-face law was a
  full-face operator ... an AngularAverageOperator for white" — PAST TENSE. LEAVE.
- foundations:1619 "The extraction achieved that by projecting ...; since B3.2 by
  typing" — past + present, both correct. LEAVE.
- foundations:"the SN registry admits only {vacuum, reflective}" — TRUE
  (SNMesh.BOUNDARY_OPERATOR_REGISTRY verified).
- foundations transpose warning "input is masked to the forward's codomain Γ₋" —
  "masked" = restricted by γ₋. TRUE.
- methods/sn/history.rst:136 "a nonzero boundary source without an inflow mask is
  uncertifiable while a masked one is delivered exactly on Γ₋" — inside a dated
  2026-07-12 changelog row. ARCHAEOLOGY. LEAVE.
- methods/sn/cartesian_multid.rst:3878 "the whole-trace B masked to a per-face set
  of inflow ordinate rows" — that is SNMaskedBoundaryOperator's SCHEDULE row-split,
  a different mask. TRUE.
- methods/sn/curvilinear_one_group.rst:2753 "varying ``n_ordinates``" — a gate's
  sweep parameter, not the retired attribute. TRUE.
- methods/sn/curvilinear_one_group.rst:2524 B3.2 narrowing bullet — TRUE.
- foundations:4051 "the wrong type tag for 'inflow-mask only'" — inside a
  past-tense "the argument WAS" disposition list. LEAVE.
- operator_algebra.rst:3895 "(n_faces, n_ordinates) boolean" — describes the trace
  space's derived mask, not an operator attribute. TRUE.

## THE BRIEF'S MEASURED-EVIDENCE BLOCK WAS PARTLY WRONG — corrected in the corpus

Brief said: "White is bit-identical on gauss_legendre(8) and product(2,4) (x and y
faces), and differs by 1 ULP (1.1e-16 / 5.6e-17) on lebedev(17) and level_symmetric(6)."

REPRODUCED (legacy body reconstructed from the git diff; O(1) random probe, 6 seeds):

| quad | face | mis-admitted by `>0.0` | old vs new | cause |
|---|---|---|---|---|
| gauss_legendre(8) | xmax | 0 | bit-identical all seeds | — |
| level_symmetric(6) | xmax | 0 | bit-identical all seeds | — |
| level_symmetric(6) | ymax | 0 | 1.11e-16 / 5.55e-17 | reduction order |
| lebedev(17) | xmax,ymax | 0 | 1.11e-16 / 5.55e-17 | reduction order |
| product(2,4) | xmax,ymax | **2** | <=1 ULP (seed-dependent) | **the classifier twin** |

So `product(2,4)` is NOT bit-identical, and — more important — the brief's framing
("bit-identical or 1 ULP from reduction order") CONFLATES two mechanically different
effects. On product(2,4) the difference is the BUG, not FP: the mis-admitted rows have
odn = +5.0e-17, cos_w = 7.85e-17 vs norm 2.5651, so Δnorm == 0.0 EXACTLY and the whole
discrepancy is ψ-weighted in the NUMERATOR -> scales with the flux ratio, unbounded by FP.
Reproduced 6.12e-05 at a 1e12 flux ratio (matches angular.py's cited 6.1e-05 exactly).

=> An O(1)-scaled probe CANNOT separate the twin from reduction-order noise. That is
why it survived. Published as a `.. warning::` in the new section.

Full-inventory audit (gauss_legendre 4/8, product 2x4/3x4/4x8, lebedev 9/17,
level_symmetric 4/6 x 6 faces): mis-admission occurs ONLY for `product`, and there
only on xmax/xmin/ymax (ymin has the same tangential count, 0 mis-admitted).
lebedev has 12 tangential/face and mis-admits 0 anywhere; level_symmetric has 0.

## FLAG — falsified, outside my permitted edit scope
1. `orpheus/sn/boundary/angular.py:82-84` — "That is every production quadrature but
   ``gauss_legendre``" is FALSE: level_symmetric(4/6) carries ZERO tangential ordinates
   on every face, and lebedev(9/17) carries 12/face but mis-admits ZERO. The true
   statement is "only the ``product`` family, and only on xmax/xmin/ymax".
   Same line carries a stray `[M]` editorial marker.
2. `orpheus/geometry/boundary/__init__.py:216-250` — four falsified claims:
   - white: "G = G_diff, the cosine-weighted Lambertian average" (B3.0 moved the
     average to R; this is the exact misassignment B3.0 corrected)
   - prescribed: "R = G = 0" (B3.0)
   - prescribed: "returns source.evaluate(psi_out.shape) masked to the face's inflow
     ordinates when the method space carries them" (B3.4a: mask retired)
   - prescribed: "The unmasked branch is reachable only for q == 0 sources" (branch retired)
   - also reflective/white "SN realises to PermutationOperator / AngularAverageOperator"
     omits the `& IdentityOperator()` TP wrap.
3. `orpheus/geometry/boundary/_base.py:334-357` — assert_source_lives_on_incoming_trace
   docstring: q-on-Γ₋ guarantee "rests on the realizer masking the evaluation"; "the
   realizer restricts the evaluation to Γ₋ at construction". Mask retired at B3.4a.
4. `tests/geometry/test_bc_equivalence_snapshot.py` + `tests/geometry/test_boundary.py`
   `_MIXED_LAW_XFAIL` reason is false (white IS narrowed); white rows still call
   `SNMethodSpace.minimal(quad)` which now RAISES. Test migration owns it.
