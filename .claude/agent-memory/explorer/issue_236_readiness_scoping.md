---
name: issue-236-readiness-scoping
description: Done-vs-remaining map for GitHub #236 (realize the spatial-scheme ⊗ angular-scheme discretization product) at main HEAD 3ac96b4 (2026-06-17), post-#240 D5b/D6. The issue's "Current state" table is STALE — cell_update.py is GONE (renamed scheme.py); LinearDiscontinuous EXISTS (slab-only); spatial registry root is DiscretizationSchemeBase; no AngularScheme/SpatialScheme enum anywhere; no diffusion_limit_consistent property; no space×angle separability probes; SpatialMomentSpace + SphericalHarmonicSpace are the typed-space peers.
metadata:
  type: project
---

# Issue #236 readiness/scoping — done-vs-remaining at HEAD `3ac96b4` (2026-06-17)

Scope: the `(spatial-scheme ⊗ angular-scheme)` discretization product. The
issue body's "Current state" table predates the #240 Step-D / D5a / D5b campaign
(all merged to `main`) and is STALE in its file:line claims. This memo verifies
every claim against source.

## STALENESS of the issue body (verify-first findings)

- `orpheus/sn/spatial/cell_update.py` — **GONE**. Renamed to
  `orpheus/sn/spatial/scheme.py`. The base `CellUpdateBase` is now
  `DiscretizationSchemeBase` (`scheme.py:567`); the Protocol is
  `DiscretizationScheme` (`scheme.py:326`).
- `SNMesh(cell_update=...)` injection kwarg — **GONE**. Now
  `SNMesh(scheme=...)` (`geometry.py:203`), default `DiamondDifference()`
  (`geometry.py:271-273`). Read as `mesh.scheme` (a plain instance attr, NOT a
  `@property`).
- "DD only — Step/LD/EC are docstring examples + NotImplementedError" — **STALE**.
  `LinearDiscontinuous` is a REAL, registered occupant (`linear_discontinuous.py:270`,
  key `"linear_discontinuous"`), N-D via UBLD. So the SPATIAL axis now has TWO
  occupants (DD + LD), not one — but LD is **slab/Cartesian-only** (raises on
  curvilinear; see ST2).
- `_run_1d_sweep` — **GONE as a name**. The 1-D walk is now
  `_OneDimScanWalk` (`loss_representation.py:1985`); the line refs `:2004 :2113
  :2167` in the issue body do not map.
- Issue body cross-ref "task #21 (LD-axis type disambiguation)" — **WRONG**.
  GitHub #21 is "MC: Heterogeneous independent reference". There is NO
  LD-axis-disambiguation issue. The LD double-use is split across #6 (LD
  ANGULAR FE) and #158 (LD SPATIAL cell update); the live typed-space follow-up
  is **#246** ("Replace 4 value-based moment-axis shape probes with one typed
  SpatialMomentSpace predicate").

## The two registries at HEAD (ground truth)

| Axis | Registry root (file:line) | `_registry_base` | Keys / occupants | Injection |
|---|---|---|---|---|
| **Spatial scheme** | `DiscretizationSchemeBase(RegistryMixin, ABC)` `scheme.py:567` | `scheme.py:851` | `"diamond_difference"` (`diamond.py:104`), `"linear_discontinuous"` (`linear_discontinuous.py:270`) | `SNMesh(scheme=...)` `geometry.py:203,271` |
| **Angular closure** | `PoleAngularClosureBase(RegistryMixin, ABC)` `pole_angular_closure.py:295` | `pole_angular_closure.py:329` | `"morel_montry_angular_sweep"` (`pole_angular_closure.py:488`), `"identity_angular_closure"` (`pole_angular_closure.py:1189`) | `SNMesh(pole_angular_closure=...)` `geometry.py:204,459-465` |

Both mount the SAME `RegistryMixin` (`orpheus/numerics/registry.py`), self-register
via `key=`, raise `KeyError` on duplicate. They are SEPARATE registries (distinct
`_registry_base` returns), so a `DiscretizationSchemeBase.create("...")` and a
`PoleAngularClosureBase.create("...")` never collide. Default angular closure is
coord-dispatched: Cartesian→`IdentityAngularClosure`, sphere/cyl→
`MorelMontryAngularSweep` (`default_angular_closure_class`, `pole_angular_closure.py:1324`).

Typed spaces (peers under `orpheus/numerics/spaces/`, both `FunctionSpace`
subclasses composing via `*` into `TensorProductSpace`, `find_factor`-queryable):
- `SpatialMomentSpace(per_axis, ndim)` — `spatial_moment_space.py:151`. `shape ==
  (per_axis**ndim,)`. `per_axis` MIRRORS `DiscretizationSchemeBase.spatial_basis_per_axis`.
- `SphericalHarmonicSpace(L)` — `spherical_harmonic_space.py:105`. `shape == (L+1, 2L+1)`.

---

## SUB-TASK 1 — disambiguate LD-spatial vs LD-angular in the type system

**(a) EXISTS.** Two SEPARATE registries (above). The SPATIAL LD occupant is
`LinearDiscontinuous` keyed `"linear_discontinuous"`; there is NO angular-LD
occupant. So the LD-spatial vs LD-angular collision the issue warns about CANNOT
arise today — they would live in different registries with different `create()`
roots. The trait split is also already encoded: `is_affine_scannable` (single-axis
1-D scannability) vs `transverse_coupling_is_facewise` (cross-axis separability) vs
`spatial_basis_per_axis` (moment count) on `DiscretizationScheme`
(`scheme.py:418-422`).

**(b) MISSING.** There is **NO `SpatialScheme` / `AngularScheme` enum** anywhere
(`grep` across `orpheus/`: only `_NamedSubgroup`, `BlockRole`, `AxisCoord`,
`CoordSystem` enums exist). The disambiguation is by STRING KEY in two registries,
not by typed enum members. There is also no angular-LD occupant at all (the angular
axis has DD-equivalent `Identity` + `MorelMontry`, no FE/LD). The issue's
acceptance phrasing "`SpatialScheme.LINEAR_DISCONTINUOUS` ≠
`AngularScheme.LINEAR_DISCONTINUOUS` as distinct symbols" is literally unmet — no
such enums.

**(c) Symbols.** `DiscretizationSchemeBase` `scheme.py:567`; `LinearDiscontinuous`
`linear_discontinuous.py:270`; `DiamondDifference` `diamond.py:104`;
`PoleAngularClosureBase` `pole_angular_closure.py:295`; `RegistryMixin`
`orpheus/numerics/registry.py:63`.

**(d) Cheapest correct first move.** The registries ALREADY disambiguate by
construction (separate roots, separate keys). #236 ST1 is therefore EFFECTIVELY
DONE for "must not have a single LD enum" — there is no enum at all, and the two
LD meanings cannot collide. The cheapest formalizing move, if a typed-symbol
surface is wanted, is to add a thin enum/StrEnum mirroring each registry's keys
(e.g. `SpatialScheme` over `DiscretizationSchemeBase.registry.keys()`,
`AngularScheme` over `PoleAngularClosureBase.registry.keys()`) — but this is
additive sugar over an already-disambiguated registry, NOT load-bearing. Lower-value
than ST2/ST3. Note: an angular-LD occupant (#6) does not exist yet, so the
"distinct LD symbols" demonstration has only one real LD (spatial).

---

## SUB-TASK 3 — extract the curvilinear angular redistribution into a selectable `AngularScheme` (the architectural heart)

**How `MorelMontryAngularSweep` enters the spatial cell-balance (traced):**

1. **Construction.** `SNMesh.__init__` binds `self.pole_angular_closure` after the
   coord block (`geometry.py:459-465`): user-supplied verbatim, else
   `default_angular_closure_class(self.coord)(self)`. For sphere/cyl that is
   `MorelMontryAngularSweep(self)`.

2. **τ enters the SPATIAL denominator.** In `cell_balance_for_streaming` /
   `cell_balance_terms` (`spatial/cell_balance.py`), the M-M weight `tau = st.tau_mm`
   (`cell_balance.py:306`) builds `c_out = alpha_out/tau` (`:313`) and `c_in =
   (1-tau)/tau·alpha_out + alpha_in` (`:314`), and `denom = 2|μ|·A_down + dA_w·c_out
   + Σ_t·V` (`:319`). So τ (the M-M ANGULAR knob, Bailey-Morel-Chang Eq.43, noted in
   the `:312` comment) sits INSIDE the SPATIAL cell-balance denominator that produces
   ψ̄. This is the shared-code coupling the issue and the coupling-memo describe.

3. **The solve-direction sweep body** (`_OneDimScanWalk.sweep`, the curvilinear
   per-ordinate arm `loss_representation.py:3066-3231`):
   - non-degenerate fast path: `ang_contrib = dA_w·c_in · ψ_{n-1/2}`
     (`:3177-3179`) ADDED into the affine scan source `b` via
     `scheme.source_emission(QV_chain + ang_contrib, inv_denom, w)` (`:3190-3192`);
   - the M-M angular thread `ψ_{n+1/2} = τ⁻¹·ψ̄ − mm_a_in·ψ_{n-1/2}` is INLINED at
     `:3216-3221` (reads `geom.tau_inv`, `geom.mm_a_in_coeff`);
   - the spatial closure `ψ̄ = (1-w)ψ_in + w·ψ_out` rides the base coefficient
     methods `scheme.cell_average` / `scheme.source_emission` (`:3190,3210`) — these
     ARE scheme-polymorphic, BUT only DD/LD implement them and LD raises curvilinear.
   - the Carlson coupled-pole seed rides `closure.psi_half_seed` (`:3114`).
   - the DEGENERATE cyl-axis ordinate is the ONE place `scheme.update(...)` is
     called per-cell (`:3162-3167`) with `angular_upstream=psi_angle[:,i]`.

4. **The apply-direction twin** (`_OneDimScanWalk._apply_walk`,
   `loss_representation.py:2151+`): curvilinear is DD-only via
   `cell_balance_for_streaming` density path (`:2375-2456`), with the M-M
   redistribution riding `pole_angular_closure.cell_contribution` (`:2370,:2448`),
   NEVER re-inlined (ERR-058/#195). The transpose carries the angular SECOND
   triangular factor via `closure.angular_adjoint` (`loss_action_transpose`, `:2576+`).

**Is the angular closure selectable?** PARTIALLY YES — `pole_angular_closure` IS a
selectable injection point with two real occupants + swappable `psi_half_seed`
sub-strategy (`AngularEdgeExtrapolation`/`CarlsonInwardSweep`/`ZeroSeed`). But it is
NOT split into the (quadrature-set × redistribution-closure) product the issue asks
for: the quadrature set is `mesh.quad` (a separate object), and the redistribution
closure is `pole_angular_closure`; there is no `AngularScheme` abstraction PAIRING
them. The M-M τ is also still INLINED into the spatial denom (`cell_balance.py:319`)
and the sweep thread (`:3216`), not factored behind a closure call there.

**(a) EXISTS.** `PoleAngularClosure` Protocol (`pole_angular_closure.py:192`) +
`PoleAngularClosureBase` registry + the M-M occupant + the Carlson seed family. The
closure is consumed in BOTH sweep (`:3075,:3114`) and matvec (`:2370`). `τ` is a
populated `StreamingTerms` field (`st.tau_mm`), set by the reduced-operator factory.

**(b) MISSING.** (i) No `AngularScheme` object pairing quadrature-set ×
redistribution-closure as ONE selectable unit. (ii) The M-M `τ` is woven into the
SPATIAL `cell_balance.py` denom and the inlined sweep thread, NOT routed through the
closure at those sites — so swapping the redistribution closure does not currently
re-route the τ-in-denom (the denom reads `st.tau_mm` off the streaming terms
regardless). (iii) No "Carlson-DD vs Morel-Montry-WDD" angular SUB-axis distinct
from the quadrature axis.

**(c) Seam.** The decoupling seam is `cell_balance_for_streaming`/`cell_balance_terms`
(`spatial/cell_balance.py:300-343`) where τ/c_out/c_in are read off `StreamingTerms`
+ the inlined M-M thread at `loss_representation.py:3216-3221`. An `AngularScheme`
abstraction would have to OWN the τ-weight production (today the reduced-operator
factory populates `st.tau_mm`) and the redistribution-thread emission, so the
spatial cell-balance consumes an angular-closure-supplied weight rather than a
pre-baked `st.tau_mm`. The `PoleAngularClosure.cell_contribution` / `precompute_psi_state`
interface (`pole_angular_closure.py:277`) is the right host for the
redistribution-closure half; the quadrature set (`mesh.quad`) is the other factor.

**(d) Cheapest correct first move.** This is the architectural heart and is the
LARGEST remaining piece. Cheapest correct first move: define an `AngularScheme`
dataclass = `(quadrature_set, redistribution_closure)` pairing that `SNMesh` holds
(replacing the two independent `mesh.quad` + `mesh.pole_angular_closure` slots with
one product object), WITHOUT yet moving the τ-in-denom — i.e. make the pairing
explicit and selectable first (a pure type/wiring move, bit-identical), THEN as a
second step relocate the τ-weight production onto the redistribution closure so the
spatial denom reads an angular-supplied weight. Proactively dispatch
`test-architect` (this crosses the spatial↔angular subsystem boundary — the L17
operator-algebra-carve proactive trigger applies).

---

## SUB-TASK 2 (curvilinear) — route a non-DD spatial scheme through the curvilinear `_OneDimScanWalk`

**Where curvilinear-LD is rejected today:** `LinearDiscontinuous` calls
`_require_slab(upstream_state)` in BOTH `update` (`linear_discontinuous.py:389`) and
`residual` (`:415`). `_require_slab` (`:191-205`) RAISES `NotImplementedError` when
`upstream_state.angular_upstream is not None` — i.e. the presence of the M-M angular
thread IS the geometry gate ("curvilinear LD is unpublished and must be derived,
Issue #158 curvilinear arm").

**Is the exclusion in `supports()`/`Compatibility` or a raise?** BOTH, at different
layers:
- The SCHEME-level exclusion is a **raise** inside `LinearDiscontinuous.update/residual`
  via `_require_slab` (`:191`). It is NOT expressed as a `Compatibility`.
- The SELECTION layer reads ONLY spatial traits: `CumprodScan.supports` →
  `mesh.is_1d and mesh.scheme.is_affine_scannable` (`loss_representation.py:849-853`);
  `ScanMarch.supports` → 1-D arm `mesh.scheme.is_affine_scannable`, 2-D arm
  `mesh.scheme.transverse_coupling_is_facewise` (`:1443-1470`). LD IS
  `is_affine_scannable=True`, so a 1-D curvilinear LD mesh would PASS `supports()`
  (it is `is_1d`) and only fail LATER inside the scheme `.update` raise. So the
  curvilinear exclusion is NOT honestly surfaced at selection time — it is a
  runtime raise downstream of a `supports()` that says "ok".

**Trace of the curvilinear scan path:** `default_for(mesh)` (`:1730`) → `CumprodScan`
for any 1-D mesh → `CumprodScan.sweep` → `_OneDimScanWalk(mesh).sweep` (`:874`). The
curvilinear arm (`:3066-3231`) inlines DD's `ordinate_scan` Blelloch recurrence and
threads M-M; the spatial closure is via `scheme.source_emission`/`scheme.cell_average`
base coefficient methods (so DD's `w=½` and a hypothetical curvilinear-LD `w=1/(1+k)`
COULD slot in) BUT the M-M `ang_contrib`/thread is woven around a SINGLE-moment
affine scan (`:3196 ordinate_scan(a, b, psi_in)` is scalar-face), and LD's
multi-moment slope/face-cochain is not carried in this curvilinear arm.

**What routing a non-DD spatial scheme through curvilinear requires:** (1) DERIVE the
curvilinear-LD closure (mathematically unpublished — issue body + `_require_slab`
docstring + `_apply_walk` docstring `:2178` all say so; this is a CORRECTNESS
prerequisite, not just wiring). (2) Carry the LD spatial-moment axis (`2^d` per cell,
`2^{d-1}` per face) THROUGH the curvilinear per-ordinate scan (today scalar-face only,
`:3196,:3051-3059` "the d=1 face cochain is moment-free for both closures"). (3)
Either lift the curvilinear arm onto `scheme.update`/`_ubld_system` per-cell (like the
degenerate path `:3162`) or extend the affine-scan coefficient model to multi-moment
under curvature. (4) Make `supports()` tell the truth — e.g. a
`scheme.supports_curvilinear` trait so `CumprodScan.supports` rejects curvilinear-LD
at selection rather than at the downstream raise.

**(a) EXISTS.** The curvilinear scan path + the DD-only inlining + the slab-LD raise
+ the spatial-closure-via-base-coefficients seam. **(b) MISSING.** The curvilinear-LD
math; the multi-moment curvilinear face cochain; an honest selection-layer exclusion.
**(c) Symbols.** `_require_slab` `linear_discontinuous.py:191`; `_OneDimScanWalk.sweep`
curvilinear arm `loss_representation.py:3066-3231`; `_apply_walk` curvilinear arm
`:2361-2456`; `CumprodScan.supports` `:849`. **(d) Cheapest correct first move:** add
the honest selection-layer exclusion FIRST (a `supports_curvilinear`/equivalent trait
read by `CumprodScan.supports`/`ScanMarch.supports`) so a curvilinear-LD mesh fails at
construction with a `Compatibility` reason instead of a deep runtime raise — cheap,
correct, and unblocks the honest selection surface; the curvilinear-LD math itself is
a #158-curvilinear / #6 research task, NOT cheap. **BLOCKER ordering:** the issue
correctly flags this as the blocker for curvilinear spatial-scheme work; but note the
ROOT blocker is the unpublished curvilinear-LD closure (correctness), the routing is
downstream of that.

---

## SUB-TASK 4 — the pair-level `diffusion_limit_consistent` query

**(a) EXISTS (the ingredients, scattered).**
- SPATIAL diffusion-limit properties are DOCUMENTED in
  `linear_discontinuous.py:7,21-22,149,305` (LM-1987: Step has NO diffusion limit;
  full LD is diffusion-limit-consistent; the flat-cut LD is "O(h²) but NOT
  diffusion-limit-consistent"). The ERR-061 reframe (slope stored in the global-x
  frame so `φ̂ = Σ_n w_n ψ̂_n` does not cancel forward-vs-backward → the
  diffusion-limit slope source `Σ_s·φ̂` is not under-driven) is in the error catalog
  (`.claude/skills/vv-principles/error_catalog.md:4485`) and pinned by
  `tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit`
  (`@catches("ERR-061")`). The `n_spatial_moments` trait the brief names is realized
  as `spatial_basis_per_axis` (`scheme.py:386-404`) + the `SpatialMomentSpace`
  `per_axis` mirror — the UBLD closure (`orpheus/sn/spatial/_ubld.py`) carries the
  `2^d` Kronecker moment basis with the diffusion-limit-load-bearing cross-moment
  `ψ̂_xy` (`linear_discontinuous.py:305`).
- ANGULAR β-condition (Bailey-Morel-Chang Eq.43) is INLINED as the `tau` weight in
  `cell_balance.py:312-314` (comment cites it). NOT a queryable β-property.

**(b) MISSING.** There is **NO `diffusion_limit_consistent` property** anywhere
(`grep`: zero hits in `orpheus/` and `tests/`). Neither the spatial scheme
(`DiscretizationScheme`) nor the angular closure (`PoleAngularClosure`) carries a
diffusion-limit trait. There is no PAIR-level query, and no β-condition property on
the angular side.

**(c) Where a pair-guard would host it.** The `Compatibility`/`supports()` surface
(`loss_representation.py:198`, `LossRepresentation.supports` `:282`) is the natural
host for a JOINT spatial×angular property query — BUT today every `supports()` reads
ONLY `mesh.scheme` (the spatial traits), never `mesh.pole_angular_closure` (verified:
no `supports`/`Compatibility` reads `pole_angular_closure`). So the pair-guard
surface exists structurally but has never been used for a spatial×angular joint
query. **(d) Cheapest correct first move:** add a `diffusion_limit_consistent: bool`
ClassVar to `DiscretizationScheme` (DD: depends; Step: False per LM-1987; LD: True)
and a `beta_first_order_consistent: bool` to `PoleAngularClosure` (M-M β=0 condition),
then a free function `diffusion_limit_consistent(scheme, closure) -> bool` (the
joint AND), exercised by a test on a known-good (LD-spatial, M-M-angular) pair and a
known-bad pair. Cheap, additive, and matches the issue's factorized
spatial-condition × angular-condition framing.

---

## SUB-TASK 5 — the space×angle separability characterization gate

**(a) EXISTS: nothing.** **(b) MISSING: everything.** There are NO space×angle
separability probes in `tests/` or `derivations/diagnostics/`. The
`$CLAUDE_JOB_DIR/.../tmp/diag_sep_{space_angle,cyl,slab,slab_iso}.py` probes the
coupling-memo references lived in an ephemeral job dir and were NEVER promoted —
`find` for `diag_sep*` returns nothing; the only `derivations/diagnostics/` separability-
adjacent file is the user's untracked `diag_s69_scanmarch_vs_window_bench.py` (a
schedule benchmark, NOT a separability gate). The grep hits for "separab"/"cross-term"
in `tests/` are generic tensor-product / operator tests (`test_tensor_product_operator.py`,
`test_measure.py`), NOT the Cartesian-additive / curvilinear-gating space×angle gate.

**(c) Symbols.** None exist. The empirical evidence (Cartesian cross-term ≈ 0;
sphere `E ≈ max(E_space, E_angle)`, |M|/max ≈ 0.37) is captured only in the
coupling-memo, not in a permanent gate. **(d) Cheapest correct first move:**
reconstruct the four probes (the coupling-memo documents the exact construction:
slab-iso control O(h²) bit-identical across N; cylinder azimuthal-floor saturation;
sphere N-gated cross-term) as a permanent characterization test under
`tests/sn/verification/` analogous to the #233 characterization gate. This is a
clean test-architect + numerics task, independent of ST1-ST4. The probes must be
re-derived from scratch (the originals are gone).

---

## ALSO — acceptance criteria (verbatim) vs HEAD

From `gh issue view 236`:
1. "Two disambiguated scheme registries with ≥2 occupants on at least the angular
   axis (DD + M-M exist) and a documented path to ≥2 on the spatial axis."
   → **MET (exceeded).** Two separate registries; angular has 2 (Identity + M-M);
   spatial now has 2 REAL occupants (DD + LD), not just a documented path — LD is
   slab-only but real. (Note: "DD + M-M" in the criterion conflates the spatial DD
   with the angular M-M; the angular occupants are Identity + M-M.)
2. "The curvilinear sweep dispatches the spatial closure through `CellUpdate` (no
   inlined DD-only closure); a non-DD spatial scheme is constructible for
   sphere/cylinder." → **NOT MET.** Curvilinear still inlines DD's scan + M-M thread
   (`loss_representation.py:3066-3231`, `:2361-2456`); the spatial closure rides base
   coefficient methods but only DD works under curvature; LD raises (`_require_slab`).
   No non-DD curvilinear scheme constructible (the closure is unpublished).
3. "A pair-level `diffusion_limit_consistent` query exists and is exercised by a test
   on a known-good and known-bad pair." → **NOT MET.** No such property exists (zero
   grep hits).
4. "The separability characterization gate is green and permanent." → **NOT MET.**
   No separability probes exist at all.
5. "Sphinx documents the product + the interference constraint." → **PARTIAL.** The
   diffusion-limit narrative is in `discrete_ordinates.rst` (LD theory, D6); the
   space×angle product + interference constraint as a named section is not confirmed
   present — would need an archivist pass to verify the
   `sn-curvilinear-aniso-norm-reconciliation` section extension.

**Met at HEAD:** criterion 1 (exceeded). **Unmet:** 2, 3, 4, partially 5.

## ALSO — sibling architecture homes for the (spatial ⊗ angular) pairing object

- **#205** (field storage × role typing — Angular/Scalar/Track/Region ×
  Flux/Source/Residual): a DIFFERENT tensor product (field STORAGE × ROLE, the
  carrier matrix), not the discretization-scheme product. NOT the natural home for
  the (spatial-scheme ⊗ angular-scheme) pairing — it types FIELDS, #236 types
  SCHEMES. Sibling discipline (same `×`-product rigor), not the same product.
- **#219** (grand-report §7 MethodSpace ABC + MethodSpaceBuilder registry +
  GeometrySpec→SpatialMesh→MethodSpace pipeline): the natural ARCHITECTURAL home.
  §7.3 `MethodSpace(Space, ABC)` is exactly the layer that would carry "this mesh's
  chosen spatial scheme × angular scheme" as a selected, validated product. BUT
  #219 is UNBUILT — `grep` confirms 0 hits for `MethodSpace`/`MethodSpaceBuilder` in
  the live SN spine; the only `MethodSpace` class is `SNMethodSpace`
  (`method_space.py:72`), which is a NARROW boundary-realizer argument (mesh + quad +
  face + trace), NOT the grand-report ABC. So #219 DESIGNS the home but has not built
  it.
- **Typed-space layer** (`orpheus/numerics/spaces/`): `SpatialMomentSpace` and
  `SphericalHarmonicSpace` are PEER `FunctionSpace` subclasses that already compose
  via `*` into `TensorProductSpace` (`find_factor`-queryable). This is the EXISTING,
  BUILT tensor-product home at the SPACE level — the pairing object's space factor
  (`SpatialMomentSpace * SphericalHarmonicSpace`) is constructible TODAY. The SCHEME
  pairing (the strategy/selection object) has no home; the SPACE pairing does.

**Verdict on where the pairing lives:** standalone for now (no built home). The
type-hierarchy destination is #219's `MethodSpace` (carries the selected
scheme×closure), with the space-factor side already realizable via
`SpatialMomentSpace * SphericalHarmonicSpace` in the numerics space layer. Until
#219 lands, a standalone `AngularScheme = (quadrature, redistribution_closure)` +
the existing `mesh.scheme` form the de-facto product on `SNMesh`.

## ALSO — #239 / #237 / #238 as prerequisites?

- **#239** (2-D ScanMarch → lift to the group-3 coefficient model so non-DD affine
  schemes ride the 2-D row-march): **PARALLEL/INDEPENDENT** of the curvilinear
  ANGULAR campaign. It is about CARTESIAN 2-D LD (the row-march), not curvilinear
  angular. It IS a prerequisite for "2-D Cartesian non-DD spatial scheme via
  ScanMarch" (ST2's Cartesian sibling) but NOT for the curvilinear angular work (ST3).
  Today `ScanMarch.supports` correctly EXCLUDES 2-D LD (reads
  `transverse_coupling_is_facewise`, `loss_representation.py:1460-1470`), so 2-D LD
  falls to the wavefront — #239 lifts ScanMarch to also carry it.
- **#237** (degenerate-cylinder pure-azimuthal matvec/sweep branch has no test
  stressor): **PARALLEL.** A test-coverage gap on the `scheme.update` degenerate path
  (`loss_representation.py:3148-3172`). Not a prerequisite, but the degenerate path
  IS the one curvilinear site already on `scheme.update`, so #237 hardens the seam
  ST2/ST3 would extend. Worth bundling.
- **#238** (do M_spatial/M_angular_redist as separately-applicable leaves still earn
  their keep post-#206): **PARALLEL/INDEPENDENT.** An operator-leaf consolidation
  question about `loss_action_decomposed` (`:2100-2149`). Orthogonal to the
  scheme-product realization; it touches the SAME `_apply_walk` decomposition the M-M
  thread rides, so a decision here would simplify ST3's matvec seam, but it is not a
  prerequisite.

**Net:** none of #239/#237/#238 BLOCK the angular campaign (ST3). #239 blocks only
the Cartesian-2-D-LD corner (ST2-Cartesian, already handled honestly by the
supports gate). All three are parallel cleanups that would slightly ease the seam.

## Recommended sequencing (cheapest-correct-first)

1. ST5 separability gate (independent, test-only, re-derive the 4 probes) — cheap,
   pins the empirical claim permanently. test-architect + numerics.
2. ST1 enum sugar (optional) + ST4 `diffusion_limit_consistent` pair property
   (additive ClassVars + a free joint query + a known-good/known-bad test) — cheap.
3. ST2 honest selection-layer exclusion (a `supports_curvilinear` trait so
   curvilinear-LD fails at construction, not a deep raise) — cheap, correct.
4. ST3 the architectural heart: `AngularScheme = (quadrature, redistribution_closure)`
   pairing on `SNMesh`, bit-identical wiring FIRST, then relocate τ-weight production
   onto the closure. PROACTIVE `test-architect` dispatch (L17 carve trigger).
5. (research, NOT cheap) curvilinear-LD closure derivation (#158-curvilinear / #6) —
   the true unblock for "non-DD spatial scheme constructible for sphere/cylinder"
   (acceptance criterion 2).
