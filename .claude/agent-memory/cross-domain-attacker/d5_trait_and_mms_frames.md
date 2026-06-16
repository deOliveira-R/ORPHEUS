---
name: d5-trait-and-mms-frames
description: Two design-time frame attacks for #240 Phase 2 Step D5. (1) The scan-march-admission trait — the precise native property separating DD/Step from LD is TRANSVERSE-COUPLING ORDER (the off-scan-axis coupling enters as a 0th-order face trace for DD/Step vs a 1st-order slope moment for LD); equivalently tensor-product SEPARABILITY of the d-D cell closure. Recommend BOTH a scheme-property trait (transverse_coupling_is_facewise / separable_cell_closure) AND keep is_scan_march_compatible as a derived alias if wanted. (2) The multi-D LD MMS ansatz ψ=[A+μx·B+μy·C]/W DOES activate the bilinear per-axis SPATIAL slope rows (verified against the actual LD Schur kernel) PROVIDED B,C carry genuine spatial variation per axis — confirmed it does not collapse to DD; one strengthening needed (a cos·cos cross term in the SLOPE drivers, not only in A) to break a residual x↔y null in the slope-source.
metadata:
  type: project
---

# #240 Phase 2 Step D5 — trait + multi-D LD MMS frame attacks

Grounded in the LIVE source (branch `feature/sn-space-angle-tier2`):
- `orpheus/sn/loss_representation.py` — `ScanMarch.supports` (`:1226-1233`),
  the inline-DD `_sweep_interior` (`:1270-1353`, the transverse fold is at
  `:1330-1331`).
- `orpheus/sn/spatial/diamond.py` — DD `cell_kernel_batch` (`:293-359`, the
  left-fold over independent per-axis face values `:350-352`).
- `orpheus/sn/spatial/linear_discontinuous.py` — `_LDCellTerms` (`:223-259`,
  slope-row `:254-259`), `_ld_source_moments` (`:200-220`), `_kernel_terms`
  d=1 raise (`:430-436`).
- `orpheus/sn/spatial/scheme.py` — `is_affine_scannable` docstring (`:351-362`).
- test-architect spec `.claude/agent-memory/test-architect/d5_nd_polymorphism_verification.md` §3.

---

## FRAME 1 — the scan-march-admission trait

### Structural feature that separates DD/Step from LD

The discriminator is NOT "spatial basis polynomial degree" directly
(that is a proxy), and NOT "1-D affine-scannability" (LD is `True` there —
the bug). The PRECISE native property is **transverse-coupling order** of
the d-D cell closure:

- **DD / Step (slopeless, P0-cell-average closures):** on a scan over axis
  x, the coupling from a NON-scanned axis y enters the x-recurrence as a
  single **0th-order face trace** `ψ_{y,in}` — a known upstream Dirichlet
  value. CONFIRMED at `loss_representation.py:1330-1331`: the y-term is
  `beta += 2·sy2·psi_y_in / D_row`, i.e. the transverse axis is reduced to
  ONE number per (group, x-cell) and folded into the scan's affine source
  `β`. The d-D DD closure is `ψ̄ = (Q + Σ_a 2g_a·ψ_in_a)/(Σ_t + Σ_a 2g_a)`
  (diamond.py:331, the explicit left-fold `:350-352`): a **sum over axes of
  independent per-axis face contributions**. No axis sees another axis's
  internal structure. This is a **rank-1 / tensor-product separable** cell
  closure — the multi-axis update factors into d independent 1-D affine
  recurrences chained through scalar face traces.

- **LD (DG-P1-upwind, carries a per-axis SLOPE moment):** the coupling from
  a non-scanned axis y is a **1st-order object** — the y-slope moment
  `ψ̂_y` — which CANNOT be reduced to a single face value because the cell's
  y-face flux varies linearly along the in-cell x-coordinate, and the
  x-slope row consumes that variation. The d-D LD per-cell system is a
  `(1+d)×(1+d)` block (one average + one slope per axis) whose Schur
  complement does NOT diagonalize across axes: the x-slope row and the
  y-slope row are coupled through the shared average. This is an
  **irreducibly axis-coupled (rank > 1, non-separable)** cell closure. The
  1-D `_LDCellTerms` Schur elimination (`linear_discontinuous.py:295-338`)
  closes BECAUSE there is only one slope to eliminate; in d≥2 there are d
  slopes and the elimination does not produce a single scalar transverse
  trace.

**The native-frame name (numerical-PDE / DG-FEM):** *the cell closure's
transverse coupling is face-wise (a trace, 0th order) vs slope-wise (a
moment, 1st order)*. Equivalently in tensor-product language: *the d-D cell
update is separable (factors into per-axis 1-D affine recurrences chained
by scalar traces) vs non-separable (an irreducible (1+d)-block solve).*

This is the SAME elegance move as Promoted Smell #16: the trait
`is_affine_scannable` is being asked to answer a question it does not
encode. `is_affine_scannable` = "my 1-D closure is a prefix-scannable affine
recurrence" (a statement about ONE axis). Scan-march compatibility = "my d-D
closure is separable into 1-D scans chained by face traces" (a statement
about CROSS-axis structure). One trait, two distinct questions — the
classic over-loaded-trait conflation. LD makes the two answers DIVERGE
(`True` / `False`), which is exactly why the bug is invisible until LD.

### Elegance assessment of the proposed `is_scan_march_compatible`

The proposed `is_scan_march_compatible` (DD/Step True, LD False) is
CORRECT and closes the bug. But it is named for the **strategy**
(ScanMarch), not for the **scheme property**. That is a frame-leak in the
trait name — the same anti-pattern this agent logged in
[[discretization_scheme_naming_signal]] as "frame-leak parameter naming":
naming a model-agnostic interface element after ONE consumer. A scheme is a
mathematical object; `is_scan_march_compatible` ties its vocabulary to one
sweep-strategy class. The moment a SECOND separability-exploiting strategy
appears (a line-SOR, an ADI/dimensional-split preconditioner, a
sweep-direction-batched KBA tiling, the diffusion-DSA consumer that #240 D4
already anticipates), the scheme would carry a trait named after a strategy
it has never heard of. The trait asserts a relationship the scheme does not
own.

The HIGHER-SIGNAL framing is to name the SCHEME PROPERTY:
- `transverse_coupling_is_facewise: bool` (DD/Step True, LD False), or
- `has_separable_cell_closure` / `cell_closure_is_separable: bool`.

Read in the domain: "my multi-axis cell closure couples non-swept axes
only through a face trace (a 0th-order value), so a scan along one axis can
absorb every other axis into its source." That reads like the
advection-discretization domain the project already established (it is the
multi-D analogue of "slopeless"). `is_scan_march_compatible` then becomes a
DERIVED reading — if you want it for the `supports()` call-site to be
self-documenting, define it as an alias/property, not as the primitive:
`is_scan_march_compatible = transverse_coupling_is_facewise` (the strategy
documents WHICH scheme property it needs, the scheme declares the property
in its own vocabulary).

**Why a basis-degree FIELD is the wrong primitive** (checked, rejected as
the trait): "P0-cell-average vs DG-P1" is the CAUSE, but encoding a
polynomial degree on the scheme over-specifies and under-discriminates. It
over-specifies because Step and DD are both "P0" yet differ on the OTHER
traits (`is_positivity_preserving`); the relevant axis here is purely the
transverse-coupling consequence, not the degree itself. It
under-discriminates because a hypothetical future closure could be
high-degree in the SWEPT direction yet still face-wise in the transverse
(a tensor-product DG-in-x / FV-in-y), which would be scan-march-compatible
despite "degree > 0". The degree is the explanation; the **transverse
coupling order is the operative invariant** the predicate must read. Put
the degree story in the DOCSTRING (it is the WHY), the boolean coupling
property in the trait (it is the WHAT `supports()` needs).

### Recommendation (Frame 1)

1. Add the scheme trait `transverse_coupling_is_facewise: ClassVar[bool]`
   (DD/Step = True, LD = False). Name it for the SCHEME property, in the
   numerical-PDE vocabulary the file already speaks. Acceptable alternative
   spelling: `separable_cell_closure`. Do NOT make the primitive
   `is_scan_march_compatible` (strategy-named frame leak) and do NOT make it
   a basis-degree int (over/under-discriminates).
2. KEEP `is_affine_scannable` as the 1-D arm of `ScanMarch.supports`
   (LD stays True there — its 1-D closure IS scannable). Gate the d≥2 arm
   on the NEW property:
   `(is_1d AND is_affine_scannable) OR (is_cartesian AND ndim==2 AND
   transverse_coupling_is_facewise)`.
3. OPTIONALLY expose `is_scan_march_compatible` as a derived alias on the
   strategy side (ScanMarch documents "I need a face-wise-transverse
   scheme") rather than on the scheme — keeps the `supports()` call-site
   readable without leaking the strategy name into the scheme's contract.

### First test (Frame 1) — discriminates the framing

The decisive experiment is NOT the routing gate alone (that catches the bug
under EITHER trait name); it is a SECOND consumer that reads the property.
Concrete: write a one-line probe that asks each scheme for the property and
assert DD/Step report face-wise/separable while LD reports
slope-wise/non-separable, INDEPENDENT of any strategy. If the only way to
get the answer is `ScanMarch.supports(mesh).ok`, the property is
strategy-entangled (the frame leak is real); if the scheme answers it
standalone, the property is in the right frame. This is the
"second-consumer" discriminator from [[discretization_scheme_naming_signal]]:
the trait should be answerable with no strategy in scope.

### Structural attack on current

`is_affine_scannable` is a per-axis (1-D) statement masquerading as a
cross-axis (d-D) statement. The `ScanMarch.supports` predicate at `:1229`
reads a single-axis property and uses it to license a multi-axis schedule.
For every slopeless scheme the two happen to coincide, so the conflation is
INVISIBLE until a scheme appears whose 1-D and d-D answers differ — LD is
precisely that scheme, and the codebase already comments the gap at
`linear_discontinuous.py:406-408` ("a multi-D LD is bilinear — an
independent slope per axis"). The trait does not encode the bilinearity, so
the router cannot see it. (Promoted Smell #16 shape-3 adjacent: one trait
typed to carry two distinct structural claims; the fix is to mint the
missing distinction as its own named property.)

---

## FRAME 2 — the multi-D LD MMS reference

### Does ψ = [A(x,y) + μ_x·B(x,y) + μ_y·C(x,y)]/W activate the bilinear
### per-axis SPATIAL slope rows?

YES — provided B and C carry genuine per-axis SPATIAL variation. The
activation analysis, against the ACTUAL LD kernel:

The LD discrete moments per cell are (average ψ̄, x-slope ψ̂_x, y-slope ψ̂_y),
each per (group, ordinate). The SPATIAL slope moment along axis a is, to
leading order, the cell-integrated `∂ψ/∂x_a` projected onto the LD slope
basis. For the ansatz:
- `∂ψ/∂x = [∂A/∂x + μ_x·∂B/∂x + μ_y·∂C/∂x]/W`
- `∂ψ/∂y = [∂A/∂y + μ_x·∂B/∂y + μ_y·∂C/∂y]/W`

So the x-slope moment is non-zero AS LONG AS `∂A/∂x ≠ 0` OR `∂B/∂x ≠ 0` OR
`∂C/∂x ≠ 0` somewhere — and crucially it is **per-ordinate** (carries an
explicit `μ_x`, `μ_y` weighting), so it is genuinely angular-resolved, not
an angularly-flat field the DD-like average path also captures. The
distinguishing question (vv Mode 7) is whether the slope is a TRIVIAL
consequence of the average — i.e. could a DD run reproduce the same scalar
flux? No: the scalar flux `φ = ∫ψ dμ = A(x,y)` (the `μ_x·B + μ_y·C` terms
integrate to zero over a symmetric quadrature — CONFIRMED this is exactly
what the test-architect spec §3 states and what the LD source split
`_ld_source_moments` assumes). But the per-ORDINATE field that the LD slope
rows discretize is `[A + μ_x·B + μ_y·C]/W`, which is NOT angularly flat —
the LD x-slope row sees `μ_x·∂B/∂x` (a genuinely μ_x-weighted spatial
slope) that the DD cell-average path CANNOT represent (DD has no slope
moment). So the ansatz exercises the NEW bilinear machinery: it is not
nulled by the closure, and DD≡LD would FAIL on it (which is exactly what
gate D5b.5 tests, the routing-flip catcher).

**Crucial subtlety — it does NOT collapse to DD ONLY IF the slope drivers
themselves vary along their own axis.** The dangerous degeneracy: if B were
constant in x (only varying in y) and C constant in y, then the x-slope
`μ_x·∂B/∂x` would vanish and the x-slope row would be driven ONLY by
`∂A/∂x` — an angularly-flat slope that, while non-zero, is the SAME slope
for every ordinate and therefore the weakest possible stress on the
per-ordinate bilinear coupling. The test-architect §3 ansatz AVOIDS this:
`B_g = b0 + b1·sin(2πx/Lx)·sin(πy/Ly)` has `∂B/∂x ≠ 0` and
`C_g = c0 + c1·sin(πx/Lx)·sin(2πy/Ly)` has `∂C/∂y ≠ 0`. So the proposed
ansatz DOES activate the per-ordinate slope rows. GOOD.

### Frame mismatch / null check (the LM-1989 slope-sign trap)

The LM-1989 trap (LD's slope-row SIGN convention — the documented
correctness trap at `linear_discontinuous.py:401-404` and the slope-row site
`_LDCellTerms.slope` `:254-259`) has TWO halves:
1. the slope-UNKNOWN sign (always exercised when the slope is non-trivially
   solved — covered by the ansatz above), and
2. the slope-SOURCE sign `Q̂` (the `_slope_source = θ·ŝ` term, `:336`) —
   exercised ONLY when a non-zero slope-moment source is supplied.

The frame hazard: **the proposed ansatz can leave the slope-SOURCE sign
half UNTESTED, and worse, can hide a residual x↔y asymmetry in the
slope-source.** Two structural issues:

(a) **Whether Q̂≠0 is reached depends on the D5b moment-source posture.** If
D5b threads only the flat (average-moment) source like 1-D Increment B
(`_kernel_terms` `rhs` at `:442` is flat, comment "flat (Q̂=0)"), the
manufactured slope-moment source is DROPPED and only the slope-UNKNOWN sign
is tested — the slope-SOURCE sign stays untested. The test-architect spec
§3 already flags this ("declare the moment-source posture; the MMS Q̂
follows it"). FRAME FINDING: the ansatz is slope-source-READY (its
analytic `∂²ψ/∂x²` etc. give a non-trivial Q̂), but the verification of the
SECOND half of the LM-1989 trap is GATED on D5b actually consuming a
two-moment source in d≥2. If D5b is flat-source in d≥2, document the scope
honestly: "multi-D slope-UNKNOWN verified; multi-D slope-SOURCE sign
deferred to the moment-source increment" — do NOT claim the multi-D LD MMS
closes the LM-1989 trap.

(b) **A residual x↔y symmetry in the SLOPE drivers can null the
cross-axis discrimination of the slope-source even when Q̂≠0.** The A-term
carries the mixed `cos(2πx/Lx)·cos(3πy/Ly)` cross term (different x/y
harmonics — good, breaks x↔y in the AVERAGE). But B and C as written are
each a SEPARABLE product `sin(k_x·x)·sin(k_y·y)` with NO cross-harmonic
mixing INSIDE the slope drivers. The x-slope row's source sees
`∂(μ_x·B)/∂x` and the y-slope row's source sees `∂(μ_y·C)/∂y`; if B and C
are related by an x↔y reflection (which `b1·sin(2πx)sin(πy)` and
`c1·sin(πx)sin(2πy)` ARE, under x↔y + b1↔c1), then a slope-row SIGN bug that
is symmetric under x↔y (e.g. the SAME wrong sign applied to both slope rows
— the most likely transcription error, since both rows share the
`_LDCellTerms.slope` code path) could produce a flux that still converges at
O(h²) because the two symmetric errors do not break the global symmetry the
MMS measures. The mixed cross term lives only in A, which feeds the
AVERAGE rows, not the SLOPE rows.

### Minimal strengthening (Frame 2)

Add a MIXED cross-harmonic term to the SLOPE drivers B and C (not only to
A), with DIFFERENT structure in B vs C so that no x↔y reflection maps one
to the other:
```
B_g(x,y) = b0 + b1·sin(2πx/Lx)·sin(πy/Ly) + b2·cos(πx/Lx)·cos(2πy/Ly)
C_g(x,y) = c0 + c1·sin(πx/Lx)·sin(2πy/Ly) + c2·sin(3πx/Lx)·cos(πy/Ly)
```
with `b2 ≠ c2` and the harmonic pairs chosen so B and C are NOT related by
x↔y. This guarantees: (i) `∂B/∂x` and `∂C/∂y` carry independent
cross-harmonic content → the x-slope-source and y-slope-source are
genuinely independent and a same-sign slope-row bug breaks the measured
flux (no symmetric cancellation); (ii) the per-ordinate field is fully
x↔y-asymmetric in the SLOPE channel, not just the average channel. KEEP the
existing discipline that is already correct: `a0_g > 0` (non-vanishing at
all four edges — the boundary closure stress, CONFIRMED right per
[[phase4_o2b_4_6_mms_ansatz_frame]]); 2G asymmetric downscatter; NON-SQUARE
domain; mixed scales.

The `a0 > 0` choice is independently confirmed against the prior phase-4
MMS frame: the curvilinear hazard there (B(0)=0 pole-regularity) does NOT
apply to a Cartesian 2-D LD — there is no `1/r` redistribution, so B and C
have NO vanishing-at-centre constraint. Cartesian slope drivers are
unconstrained at the boundary; the only boundary requirement is that the
prescribed vacuum-inflow BC be consistent with the (non-vanishing) interior,
which `a0 > 0` + a non-zero slope at the inflow face already supplies.

### Structural independence of the source (L11) — CONFIRMED

The source is derived by substituting ψ into the CONTINUOUS within-group
transport PDE
`μ_x ∂_x ψ + μ_y ∂_y ψ + Σ_t ψ = (1/W)(Σ_s φ + Q^ext)`
and solving symbolically for `Q^ext_n(x,y)` (SymPy, algebra-of-record
Branch 1). This touches NONE of the LD cell-update code: not
`_LDCellTerms`, not `_schur_terms`, not `_kernel_terms`, not the multi-D
bilinear kernel D5b will write. The manufactured `φ = A(x,y)` (isotropic
moment) and the per-ordinate streaming derivative carry the full ansatz —
all from the PDE, none from the discretization. So a sign bug or
index-swap in the LD slope-row code CANNOT also corrupt the reference (the
defining property of an independent MMS reference, L11). CONFIRMED
structurally independent.

### Elegance assessment of the proposed framing (Frame 2)

The μ-linear ansatz `[A + μ_x·B + μ_y·C]/W` is the RIGHT frame for a
BILINEAR (not biquadratic) closure: it is the native truncated-P1 Legendre
element in angle (matching [[phase4_o2b_4_6_mms_ansatz_frame]]'s finding
that `(A+μB)/W` IS the truncated-P1 element and a linear-in-μ input FULLY
probes a linear angular closure — no P2 enrichment needed for the ANGULAR
axis). The per-axis independence (`μ_x·B` separate from `μ_y·C`) correctly
mirrors the per-axis independence of the LD SPATIAL slopes (the
bilinearity: two slopes that do not collapse). The frame is correct; the
only gap is in the SPATIAL structure of the slope drivers (the
cross-harmonic null above), not in the angular ansatz. One caution carried
from the phase-4 frame: the spatial polynomial degree must be high enough
that the quadrature is exact for the manufactured source — the
cross-harmonic strengthening raises the spatial frequency, so confirm the
angular quadrature order resolves the `μ_x·μ_y`-bilinear moments the
streaming term generates (a Mode-7-adjacent quadrature-exactness check,
not a closure-activation issue).

### First test (Frame 2) — discriminates the strengthened from the weak ansatz

Run the WEAK ansatz (B,C separable, x↔y-related) and the STRENGTHENED
ansatz (cross-harmonic, x↔y-broken) against a DELIBERATELY same-sign-flipped
LD slope row (mutate the `_LDCellTerms.slope` sign on BOTH axes
identically). PREDICTION: the weak ansatz still converges at ~O(h²) (the
symmetric error cancels in the measured flux — a false green); the
strengthened ansatz breaks the order or the value band (the bug is caught).
If both catch it, the strengthening is unnecessary; if only the
strengthened one catches it, the strengthening is load-bearing. This is the
experiment that cannot-fail-is-rejected test: it discriminates the two
ansatz designs by their power against the specific LM-1989 same-sign trap.

---

## CROSS-METHOD POLLINATION (current method: SN spatial closure)

- **From diffusion / DSA (the #240 D4 anticipated second consumer):** the
  `transverse_coupling_is_facewise` property is exactly the predicate a
  dimensional-split (ADI) or line-SOR diffusion preconditioner needs to
  decide whether it can sweep one axis at a time. The trait, named for the
  scheme property not the strategy, is reusable by the DSA consumer with no
  rename — the payoff of NOT calling it `is_scan_march_compatible`. First
  test: when the diffusion-DSA call site lands, confirm it reads the SAME
  trait (no new scheme property minted for the same separability question).
- **From MOC (characteristic separability):** MOC's along-track 1-D
  attenuation is the canonical "scan along the characteristic, couple
  transverse via face traces" — the same separability structure. A
  face-wise-transverse SN scheme is the MOC-track analogue; an LD-style
  in-cell slope is the MOC source-region polynomial that breaks pure track
  separability. Confirms the property is method-spanning (the right
  abstraction level for the trait).

---

## UNEXPLORED

- Category theory / operadic composition — no compositional-structure
  trigger beyond the separability already captured; rank/tensor framing is
  the operative one.
- Homology (B∘B≠0) — refuted for this subsystem in
  [[issue_208_operator_algebra_frames]]; no boundary-of-boundary structure
  in the cell closure.
- Group theory / SO(3) — the angular quadrature symmetry is upstream of the
  spatial closure; not the discriminator here (the trait is a spatial-cell
  property).
- Wiener-Hopf / half-space — wrong solver family (a within-cell algebraic
  closure, not a transport half-space); refuted previously
  ([[trajectory_resolvent_foreign_frames]]).
