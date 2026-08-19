# Part III-S (space layer) grounding — plan vs HEAD of main

Agent A (explorer). Reconciles `.claude/plans/orpheus-operator-machinery-report-v2.md`
§I.8-I.10, corrections #16-#19, Part III-S, the F-rows of the stage-inversion table,
and Phase S S1-S7 against **main @ `7aae9bf1` (2026-08-19)**. Audit epoch: plan "read
from the tree 2026-08-08". All `file:line` are at this HEAD; re-derive via Nexus after
edits. Nexus MCP tools were loaded; primary evidence here is grep/read/git/live-probe
because most claims are literal-text or dataclass-field claims (per nexus-tools routing).

## 0. The headline

**[M] The space layer is structurally FROZEN since the audit.** `git log --oneline
--since=2026-08-08 -- orpheus/numerics/spaces/ orpheus/numerics/space.py
orpheus/numerics/coupled_system.py` returns exactly **2 commits**, both verified
docstring-only in this scope: `7c250811` (2026-08-17, 1-line catalogue-path repoint in
spherical_harmonic_space.py) and `c33178ef` (2026-08-12, 3-line citation fix in
radial_characteristic_space.py; its own commit message: "Production diff is
DOCSTRING-ONLY — independently re-verified by AST comparison"). `space.py`,
`coupled_system.py`, `full_field_space.py`, `angular_trace_space.py`,
`scalar_trace_space.py`, `spatial_moment_space.py`: **0 commits** since 2026-08-08.

Consequence: every III-S structural claim survives or fails on what was true at the
audit epoch — nothing landed in between to overturn it. The reconciliation deltas are
of three kinds: (a) **three space files the audit never read** (all predate it — unread,
not new), which **widen** F1/F3; (b) two of the plan's cited sites are **dead factories**
(the live sites are elsewhere, same doctrine); (c) part of S6's target already existed
**pre-audit** under concrete-operator names. Of the brief's four post-audit campaigns,
only B3.3 (`65efbb52`, 2026-08-14, operator.py) and Q6 (measure/registry side) touch this
slice's neighbourhood; G6.x and B3.2/B3.4 landed **before** 2026-08-08 and are already
baked into the audit's praise list; #344 does not touch these files.

## 1. Claim-by-claim verdict table

| # | Plan claim (III-S / corrections) | Verdict | Evidence |
|---|---|---|---|
| P1 | Γ ladder praise: tangential third class first-class | **CONFIRMED** | `angular_trace_space.py:166-171` (TraceRole; "the three do NOT partition"), `:776-784` ([M] counts 0/8/12 on GL4/product(4,4)/lebedev(17)) |
| P2 | Γ ladder: quotient statement (ker G = tangential excess) written down | **CONFIRMED** | `:786-810` (`.. important::` block: `ker G_Γ(f) = Γ(f) \ (Γ₊⊔Γ₋)`; quotient ≅ Γ₊⊕Γ₋) |
| P3 | partial-current metric `\|Ω·n\|⊙w` installed unconditionally + Moore-Penrose pseudo-inverse | **CONFIRMED** | `:276-312` (`_build_trace_metric_weights`), `:382-402` (`partial_current_metric` raises TypeError on None), `space.py:299-322` (`apply_inverse_metric` masked pseudo-inverse) |
| P4 | `omega_dot_n` single source of truth | **CONFIRMED** | `:224-273` (`build_omega_dot_n`, "SINGLE face-name → signed-projection primitive", public per #52/ERR-041), `:345-348` |
| P5 | γ± bound to typed half-spaces refusing cross-face composition despite equal shapes | **CONFIRMED** | `:494-538` (TraceRestrictionOperator bound `Γ(f)→Γ±(f)`, G6.3 step 1 / #330), `:728-742` (`_face_space_name`: [M] `\|Γ₊(xmin)\|=\|Γ₊(xmax)\|` on 3/3 shipped-rule probes). Live probe at HEAD: `Γ₊(xmin) == Γ₊(xmax)` → **False**, `Γ₊(xmin) == Γ₋(xmax)` → **False**, both shape (4,2) |
| P6 | View-G doctrine consistent | **CONFIRMED** | `space.py:131-143` (View-G paragraph, #205/#207; units on Field) |
| P7 | `find_factor` by type; dual idempotence; ⛔ post-mortems | **CONFIRMED** | `space.py:473-516`; `:597-611` (`DualSpace.of` unwraps); `:223-247` (the ⛔ broadcast note) |
| P8 | CoupledSpace/CoupledOperator praise (offsets-as-structure, member-wise metric dispatch, adjoint free, missing block = structural zero) | **CONFIRMED** | `coupled_system.py:539-559` (`system_slices` "This IS the scoped LocalToGlobalMap"), `:573-609` (member-wise dispatch), `:39,650,742,1167` ("None IS the zero map, structurally") |
| F1a | Identity nominal `(name, shape)`, metric `compare=False` | **CONFIRMED** | `space.py:153-173` (quoted §2 below); replicated by explicit `__eq__`/`__hash__` delegation in **6/6** subclasses probed: TensorProductSpace :467-471, DualSpace :591-595, FullFieldSpace :206-210, AngularTraceSpace (inherits; fields compare=False :378-380), SphericalHarmonicSpace :157-161, CoupledSpace :459-463, plus `_RadialCharacteristicSubSpace` :368-372 |
| F1b | Weak-equality "feature" prose present | **CONFIRMED + WIDER** | plan's two quotes verbatim (`angular_trace_space.py:340-343`, `full_field_space.py:74-82`); **two more sites in files the plan never read**: `spatial_moment_space.py:90-99`, `spherical_harmonic_space.py:129-132` (§3.F1) |
| F1c | Stringly `_face_space_name` discriminators live | **CONFIRMED (minting only)** | def `:728`, use `:602`, ref `:767`. [M] grep `\.name ==\|\.name\.startswith` over `orpheus/**/*.py`: the ONLY production hit is `space.py:170` (the `__eq__` itself) — identity is name-ENCODED, nothing name-PARSES (`:753-755`: "nothing dispatches on it") |
| F2 | No axis primitive; metric binding by broadcast convention; ⛔ note; cross-space axis orders | **CONFIRMED, one tense caveat** | §3.F2. The two-methods-two-axes divergence is **historical** (fixed 2026-08-04, pre-audit; ⛔ note `space.py:223-247` says so; [M] live probe: consistent at HEAD) |
| F3 | Three volume-weight answers in three places | **CONFIRMED + WIDER (4-5 doctrines over 6+ sites)** | §3.F3; no consolidation since audit ([M] §0) |
| F4 | TensorProductSpace materializes dense outer weight | **CONFIRMED, measured** | `space.py:397-408` (`np.multiply.outer`, `np.ones(f.shape)` for Euclidean factors); [M] probe: (40,8,2)⊗SH(L=3) → weight tensor (40,8,2,4,7) = 17,920 elements = exactly the state size. Escape: all-Euclidean → None (`:397-398`) |
| F5 | Missing: WithTrace, cone/basis-kind, parity, GroupAxis, 𝔽 | **CONFIRMED, 5/5 absent from the space layer** | §3.F5 (greps quoted); parity/symmetry data EXISTS one level down on measures/registry (partial sibling) |
| Frag | Three parallel product mechanisms; no consolidation | **CONFIRMED (sites refined)** | §3.Frag |
| Drift | Torsor prose verbatim in coupled_system.py | **CONFIRMED** | `:107-110` and `:237-240`, quoted §3.T. Phase 0.7 not executed |
| #16 | "Two copies of ℝⁿ the same regardless of inner product" docstring doctrine | **CONFIRMED still present** | `space.py:145-150` verbatim; 0.8 not executed |
| #17 | Weak equality as feature → defect | **CONFIRMED live** (see F1b) | — |
| #19 | SNMethodSpace misnamed, kept | **CONFIRMED; 0.9 not applied** | `sn/mesh/method_space.py:1` ("SN method-space carrying mesh + quadrature + trace metadata"), class at `:72`, docstring `:73` repeats it. No "generator record" vocabulary anywhere in the file ([M] grep). NB it never subclasses FunctionSpace — the claim is in the NAME + prose only |
| #20 | Phase S prerequisite of 1.1/1.4/1.6 | not adjudicable from the tree (a sequencing ruling); the premises it rests on (F1-F3 live) are confirmed above | — |

## 2. F1 — identity (the load-bearing quotes + live measurements)

(a) `space.py:153-165` at HEAD:

```python
    name: str
    shape: tuple[int, ...]
    # Metadata, NOT identity (see class docstring): two spaces with the same
    # (name, shape) are equal regardless of the installed metric, and hashing
    # is on (name, shape). ``compare=False`` keeps the weights out of the ...
    inner_product_weights: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )
```

with manual `__eq__` (`:167-170`: `self.name == other.name and self.shape == other.shape`)
and `__hash__` (`:172-173`). Class docstring `:145-150` still states the #16 doctrine
verbatim: *"two copies of ℝⁿ are 'the same' space regardless of which inner product is
installed."*

(b) Weak-equality prose, census over `grep -rn "compare equal" orpheus/numerics/` — 7 hits,
4 of them doctrine sites:

- `angular_trace_space.py:341-343` — "two trace spaces on meshes of the same total
  boundary size compare equal regardless of their face decomposition" (the plan's quote).
- `full_field_space.py:75-78` — "two composites over meshes of the same total dimension
  compare equal, so the `OperatorSum` composition guard accepts the full within-group
  loss `L + C - S - B`" (the plan's "confession" quote, verbatim).
- **NEW to the plan** `spatial_moment_space.py:92-99` — "a `(per_axis=4, ndim=1)` space
  and a `(per_axis=2, ndim=2)` space both have `shape == (4,)` and so compare equal.
  The two never coexist on one mesh … so the size-identity convention is the right one."
- **NEW to the plan** `spherical_harmonic_space.py:129-132` — "equal-L spaces compare
  equal even when their `inner_product_weights` arrays are distinct ndarray objects."
  ⚠ Distinguish: this alias is **benign by construction** — L determines the metric via
  `from_L`, so two equal-L spaces cannot carry different metrics; the SpatialMomentSpace
  one is a genuine cross-factorisation alias.

Live measurements at HEAD (`.venv/bin/python` probe, GL8 vs `Quadrature.product(4,2)`
(N=8), 2-face layout (8,2)/(8,2)):

- `AngularTraceSpace(GL8) == AngularTraceSpace(product(4,2))` → **True**, while their
  partial-current metrics are **not** allclose — the plan's "Γ₊(x_max) on
  gauss_legendre(8) equal to the same on any 8-point rule" cross-QUADRATURE aliasing,
  measured.
- `FullFieldSpace.from_blocks(bulk(8,2,10), t) == from_blocks(bulk(4,4,10), t)` →
  **True** (both flat shape (192,)) — cross-decomposition aliasing, measured.
- `SpatialMomentSpace(4,1) == SpatialMomentSpace(2,2)` → **True** — the unread-file
  alias, measured.
- Cross-face: `Γ₊(xmin) == Γ₊(xmax)` → **False** — the hand-minted name discipline holds
  where it was installed.

(c) See table F1c: minting-only; one `.name ==` site in production (the `__eq__`).

(d) **#369 (OPEN, `module:data level:L0`)** is the same defect class on the measure side:
`ReferenceMeasure.name` documented as "canonical mathematical identity — what equality
compares" (`exactness.py:161-164` per the issue) while Lebesgue-on-[-1,1] has **three
mutually-unequal spellings** (`LEGENDRE` / `LEGENDRE.on(-1,1)` / `UniformMeasure` from
`equispaced`); causes: cross-class dataclass `__eq__` returns NotImplemented before any
field is read, and `.on()` mints a new name instead of canonicalising. Its "Suggested
shape" (partition the classes into disjoint domains + canonicalise `.on()` + single-source
the support string) **is S3's interning-at-generators cure one level down**. Map: F1's
"identity by canonical construction" charter (0.8/S3) should CITE #369 as the second
instance — the campaign that fixes space identity and leaves measure identity nominal has
fixed half the defect class. (The issue itself notes the current fence: a registry gate
closing the vocabulary at `{LEGENDRE, UNIFORM_ON_SPHERE}`.)

## 3. F2 — axis/metric binding

(a) The ⛔ note is present at `space.py:223-247` and is **accurate as history**: "This
routed through numpy's default (TRAILING) broadcast **until 2026-08-04**, while
`_diagonal_apply_metric` used the LEADING convention — the same metric applied along
different axes by two methods of one space… [M] 456 vs 552 on a (3,3) probe." At HEAD
both `_diagonal_inner_product` (`:251`) and `_diagonal_apply_metric` (`:297`) and
`_diagonal_apply_inverse_metric` (`:320`) route through the shared `_broadcast_metric`
(`:267-277`, LEADING: pad trailing 1s). [M] live probe: on the SH production layout
`(L+1, 2L+1, ng, nx)` `inner_product(x,y) == sum(apply_metric(x)*y)` to machine
precision at HEAD. ⟹ **the plan's "already bit" phrasing is right; a Phase-S author must
not read F2 as a live two-axes divergence.** The live defect is what the plan's F2
headline says: *no axis primitive — binding is a broadcast CONVENTION* (leading-axis by
position), which no type checks and every space must obey by discipline.

(b) Axis orders actually in use at HEAD (the census; "per space" = the 6 space files +
the 2 production mint sites):

| space / element | layout | where |
|---|---|---|
| SN bulk (production) | **(N, ng, \*spatial)** — ordinates leading | `sn/mesh/augmented_mesh.py:1085` (`bulk_shape = (N, self.ng, *self.spatial_shape)`); metric stored `(N, 1, *spatial)` `:1077-1084`; LD appends a TRAILING moment axis with the moment-mass metric `:1095-1098` |
| AngularField mints (transport) | **(N, ng, \*spatial)** | `transport/fields/_bases.py:374-381` |
| `angular_flux_space` factory | **(n_cells, n_ordinates, n_groups)** — cells leading | `space.py:663-667` — ⚠ **DEAD in production**: [M] grep over `orpheus/ tests/`: consumers = `tests/numerics/test_space.py` only (`:17,152-166`) |
| ScalarField mints | (ng, \*spatial), Euclidean | `_bases.py:497-503` |
| SphericalHarmonicSpace | space shape **(L+1, 2L+1)**; weighted elements are **(L+1, 2L+1, ng, \*spatial)** (harmonic axes leading, metric broadcast over trailing) | `spherical_harmonic_space.py:109-116`; production-layout statement `space.py:234-235` |
| AngularTraceSpace | flat `(total_size,)`; per-face slots **(N, ng[, n_face_cells])**, ordinate axis 0 | `angular_trace_space.py:291-295` |
| ScalarTraceSpace | flat; slots **(2, ng, \*face_spatial)** — component axis 0 | `scalar_trace_space.py:62-67` |
| SpatialMomentSpace | `(per_axis**ndim,)` — a trailing axis composed onto bulk via `*` | `spatial_moment_space.py:68-71` |
| RadialCharacteristic legs | (ng, nx) cells + (ng,) corner per (level, sign) leg | `radial_characteristic_space.py:325-346` |

⟹ the cross-space inconsistency is **real** (bulk ordinate-leading vs SH harmonic-leading
vs scalar group-leading vs trace flat-with-slot-shapes — every inter-space arrow still
carries an implicit transpose no type sees), **but** the plan's specific pairing
"`angular_flux_space` = (cells, ordinates, groups)" cites the dead factory; the live bulk
convention is uniformly (N, ng, \*spatial) across the two production mint sites. Verdict:
**CONFIRMED with the site CHANGED**.

## 3.F3 — volume-weight ownership at HEAD

The plan says "three answers in three places". [M] the census at HEAD is **4-5 doctrines
across 6+ sites** — wider because two of the sites live in files the audit never read:

1. **Euclidean, volumes absorbed into operators** — `space.py:673-678`
   (`scalar_flux_space` docstring: "cell volumes are absorbed into the operator that
   produced φ"; factory dead, see F2b) ; **live twin**: `_bases.py:497-503` ScalarField
   mints with no weights; `spatial_moment_space.py:101-105` ("the cell-volume /
   mass-matrix weighting lives on the UBLD cell operator, not on the field space").
2. **Space owns G_bulk = V_cell·w_n** — `full_field_space.py:19-21,166-169`;
   built at `augmented_mesh.py:1077-1098` (incl. the LD moment-mass tensor factor,
   #310 C2 ruling 3 — a hand-rolled per-axis metric product at the mint site).
3. **Trace refuses None; explicit metrics always** — `angular_trace_space.py:395-402`
   (`partial_current_metric` raises TypeError on None), `:659-675` (`current_space`
   metric "EXPLICITLY unit, never None — None encodes TWO states").
4. **Space owns G_sd = V_cell (ψ½ state metric)** — `radial_characteristic_space.py:334`
   (interior), `:337` (boundary corner: G = V(r=R)), strictly-positive guard `:246`
   (ERR-067 fix `eb5943d6`, 2026-07-06 — **pre-audit**, and it moved the metric TOWARD
   space ownership). File unread by the plan.
5. **Space owns the face-AREA metric** — `scalar_trace_space.py:35-41` (boundary surface
   measure; "angular weights are NOT re-applied here"). File unread by the plan.

Nothing consolidated post-audit ([M] §0: zero semantic commits). B3.3 touched operator.py
only; ERR-067's fix predates the audit. Verdict: **CONFIRMED and worse than stated** —
the trace-layer doctrine the plan calls correct now has TWO more space-owned metrics
agreeing with it (4, 5), so "resolve F3 toward the space" (S2) has more allies and the
same lone holdout (the scalar/bulk-Euclidean arm, doctrine 1).

## 3.F4 — measured (see table). Constructor quote, `space.py:399-408`:

```python
    result: Optional[NDArray] = None
    for f in factors:
        w = (f.inner_product_weights if f.inner_product_weights is not None
             else np.ones(f.shape))
        w = np.broadcast_to(w, f.shape)
        result = w if result is None else np.multiply.outer(result, w)
```

Note the aggravator S1 should name: a single weighted factor forces `np.ones(f.shape)`
materialization for every Euclidean factor too, then the outer product — the O(state)
tensor exists whenever ANY factor is weighted.

## 3.F5 — missing capabilities, each with its grep (all at HEAD)

| capability | verdict | grep |
|---|---|---|
| WithTrace / Data tag | **ABSENT** | `grep -rn "WithTrace" orpheus/` → 0 hits |
| cone / basis-kind (Nodal\|Modal) | **ABSENT** from numerics | `grep -rn "Nodal\|Modal\|basis_kind" orpheus/numerics/` → 0; `grep -rni cone orpheus/numerics/` → 2 hits, both `symmetry.py` "cone point" (unrelated geometric usage) |
| parity data on angular axes (space layer) | **ABSENT** | `grep -rni parity orpheus/numerics/space.py orpheus/numerics/spaces/*.py` → 0. **Partial sibling one level down**: `invariance_group` tag on measures (`measure.py:72`, `quadrature/__init__.py:6-15`) and `half_range_clean` per registry entry (`quadrature/registry.py:475-505,763-793`) — the raw data a `QuadratureAxis.parity` would forget from. ⚠ these are DECLARED tags (construction-time facts), not computed checks |
| GroupAxis | **ABSENT** — energy is a bare `int` in every shape tuple | `grep -rn "GroupAxis\|QuadratureAxis\|SpatialAxis\|HarmonicAxis" orpheus/` → 0; e.g. `augmented_mesh.py:1085` (`self.ng`), `space.py:665`. `EnergyGroupSpace` exists only as "Future direction… NOT shipped" prose (`space.py:45-46,61-63`; `spaces/__init__.py` "Future:") |
| ground field 𝔽 | **ABSENT** | `grep -rni "complex\|sesquilinear\|ground.field" orpheus/numerics/space.py orpheus/numerics/spaces/*.py` → 0 semantic hits; every dtype in the space files is `float` |

**One F5-adjacent thing the plan does not know**: `FunctionSpace` became
carrier-generic — `Carrier = TypeVar("Carrier", default=Any)` (PEP-696, `space.py:92-100`),
with `FullFieldSpace(FunctionSpace[CompositeField])` and
`CoupledSpace(FunctionSpace["CoupledField"])` and a structured-carrier Protocol seam
(`full_field_space.py:113-148`). Pre-audit, but absent from the plan's F5 inventory; it is
the surface S1's `axes` field would be added beneath, and the existing precedent that
"space surface is polymorphic over the element type" — S4's WithTrace stress-test
(per-face polynomial trace content) has this seam available.

## 3.Frag — composition fragmentation

**CONFIRMED: still exactly three parallel ⊗/mint mechanisms + two ⊕ realizations; zero
consolidation since the audit.** Census at HEAD:

1. **Monolithic shape-tuple mints** (never touch product machinery) — TWO live sub-sites:
   `augmented_mesh.py:1099-1102` (`FunctionSpace(name="sn_bulk", shape=(N, ng, *sp, [m]))`
   — a bare base-class mint, not even the space.py factory) and `_bases.py:374-381 /
   :497-503` (per-`_SPACE_NAME` mints for every Angular/Scalar field leaf). The plan's
   named exemplars `angular_flux_space`/`scalar_flux_space` are **dead** (test-only, F2b).
2. **`TensorProductSpace`** — the moment-carrier path only: [M] consumers =
   `transport/fields/harmonic_moment_flux.py`, `_bases.py:563` (space-building),
   `moment_displacement.py`, `harmonic_moment_source_sink.py`, `spatial_moment_space.py`
   (composes via `*`).
3. **`CoupledSpace`** (`coupled_system.py:398`) — the N-system ⊕, `from_systems`.
   Plus **`FullFieldSpace.from_blocks`** — the 2-block bulk⊕trace ⊕ (S1 already knows:
   "FullFieldSpace and the Γ ladder already are direct sums; CoupledSpace becomes the
   unique ⊕ target"), now with a **second instance family** the plan doesn't mention:
   System B's ψ½ composite instantiates the SAME class under
   `name="radial_characteristic"` (`full_field_space.py:180-189,229-234`).

Sibling smell, carrier side: **#297** (three member-contract spellings
Vector/`_is_ravellable`/SystemField) is the same Cardinal-2 shape one type-level down —
a Phase-S consolidation that unifies spaces but not carriers leaves the fragmentation's
other half standing.

## 3.Γ — the ladder after B3 (the brief's git-log question)

`git log --oneline --since=2026-08-08 -- orpheus/numerics/spaces/ orpheus/numerics/space.py
orpheus/numerics/coupled_system.py` → 2 commits, both docstring-only (§0). What DID move
is the **operator side**: `65efbb52` (B3.3, 2026-08-14) retired `IncomingOrdinateMaskTensor`
from `orpheus/numerics/operator.py` (−133 lines). [M] its commit body: B3.2 (pre-audit)
had narrowed the SN boundary law's domain to Γ₊, so vacuum realizes to the honest zero map
Γ₊→Γ₋ and the mask's rows "left the operator's domain entirely"; zero production
construction sites remained; the mask's own docstring had proved it was `I − P_in` with
range Γ₊⊕Γ_tan (four of eight ordinates tangential on a production cylinder).

⟹ **the praise is still accurate, and B3.3 strengthened it**: the retired mask was
precisely the "not-inflow spelled as outflow" hazard the ladder's tangential third class
exists to forbid (`angular_trace_space.py:168-171`); its removal leaves the typed
`TraceRestrictionOperator` route as the only spelling. No post-audit change weakened any
of the five praised properties (P1-P5 table rows). One pre-audit correction the plan may
not have absorbed: the G6.1-era claim that the response R is an endomorphism of Γ₋ was
⛔-refuted in-file on 2026-08-04 (`:709-718` — "no realized response is an endomorphism
of Γ₋… Nothing here is self-adjoint"); the plan's praise text doesn't repeat the refuted
claim, so no action, but any Phase-S prose citing "R : Γ₋ → Γ₋" should read that block.

## 3.T — torsor prose (correction #15 / drift row / 0.7)

Verbatim at HEAD, two sites in `coupled_system.py` (this file only; the tree-wide torsor
question is sibling D's):

- `:107-110` (module docstring, closing prose before References): *"Role semantics stay
  on the members: the machinery never adds arrays — every `+` it evaluates is a member
  `+`, so the **affine flux torsor**, displacement minting, and units law of the member
  family apply unchanged (the same delegation discipline as `Composite._map_binary`)."*
- `:237-240` (CoupledField class docstring): *"every `+` is a member `+`, so the member
  family's role semantics (**affine flux torsor**, displacement minting, units) apply
  unchanged — this class adds structure, never arithmetic."*

**For sibling D**: `grep -rn torsor orpheus/` finds the vocabulary live across
`transport/__init__.py`, `transport/full_field.py:107-110,352-358`,
`transport/fields/_coefficient_role.py:14,65,81`, `transport/displacements/*` (≥13 files)
— the member-family docs correction #15/0.7 says must be swept. Also for D: **#331**
(OPEN) measures the three leaves of `A=(L+C)−S−B` disagreeing on whether their domain
includes the displacement space (L accepts, S and B TypeError) — the domain-identity gap
adjacent to F1 but owned by the flux-ontology scope.

## 4. Phase S rows — target-defect + partially-landed-elsewhere verdicts

| row | target defect still exists? | machinery partially landed under another name? |
|---|---|---|
| **S1** Axis primitive; delete `_broadcast_metric`; kill F4 densification | **YES** — `_broadcast_metric` live at `space.py:267-277`, called from `:251/:297/:320` and re-relied-on by the ladder (`angular_trace_space.py:588-593,825-829`); F4 measured live | ⊕ target exists (CoupledSpace); carrier-generic surface exists (§3.F5 note). ⚠ S1's "deletes `_broadcast_metric`" must also rewire the two docstrings that TEACH the leading-broadcast convention as the contract (`:588-593`, `:825-829`) |
| **S2** leaves on axes; volumes out of operators; GroupAxis; one canonical axis order | **YES** (F3, F5, F2b all live) | **Substantial partial siblings, none an Axis**: `face_layout.py` (index-set layout + the `\|Ω·n\|w` metric-builder `face_streaming_normal` — the trace tier's "axis bookkeeping"); `moment_layout.py` (moment-axis index conventions: `AVERAGE_MOMENT`, tail policy — physics-free, homed in numerics per #245); `exactness.py` + `generating_measure.py` (2026-08-02, pre-audit but plan-unread) + `quadrature/registry.py` (Q6, post-audit `b02dd536` 2026-08-14) carry per-measure **reference/degree/symmetry** metadata (`ReferenceMeasure`, `invariance_group`, `half_range_clean`) — i.e. the GENERATOR side of §I.8's `QuadratureAxis = forget(Quadrature)` is now much richer than at v1; S2's job on the angular axis narrows to *forgetting* (weights + parity closure) rather than inventing. GroupAxis: nothing anywhere ([M] F5) |
| **S3** structural identity via interning | **YES** — F1 measured live | Precedents named by the plan check out: [M] SNMesh caches the trace space (`augmented_mesh.py:766` `angular_trace`), `RegistryMixin` (`numerics/registry.py`) is registration-at-class-creation. **#369 is the measure-side twin** and its proposed fix is the same interning/canonical-construction move — S3 should cite it (§2d) |
| **S4** WithTrace/cone-rep/parity tags | **YES** — all absent ([M] F5) | Parity raw data on measures/registry (F5 row 3); nothing for WithTrace or cone-rep. The plan's own weakest-confidence flag (per-face polynomial trace content vs a boolean) now has a concrete stress case in the tree: `_RadialCharacteristicSubSpace`'s split interior/corner pair — a "trace" whose content is per-(level,sign) legs, not per-face ordinates |
| **S5** ground field 𝔽 + conjugation | **YES** — absent ([M] F5) | none |
| **S6** retract constructors owned by the axis | **YES** as stated (no `Contract`/`Embed`, no minted-together adjoint pair) | **PARTIALLY-LANDED (pre-audit) as concrete operators**: the Lambertian is factored `IsotropicEmissionOperator @ PartialCurrentOperator` (G6.3 step 3, 2026-08-04 — `sn/boundary/angular.py:31-38` module history, `:65`, `:191`; the welded `AngularAverageOperator` retired `b4f0f5c9`); `PartialCurrentOperator` IS the contraction retract π_{\|Ω·n\|w} with "its transpose is a theorem, not a hand-roll" (`:87`); the named intermediate S(f) exists (`angular_trace_space.py:639-699` `current_space`); `TraceRestrictionOperator` is bound to the ladder's typed spaces (G6.3 step 1). What S6 still adds: the axis-OWNED constructor minting (π, ι) together with ι = π^H-up-to-scale **by construction** (today the pairing is delivered by the algebra's `.H` over the bound spaces, per-operator), fission's angular factor and P₀ projection re-expressed through the same retracts, and "partial-current normalization not hand-written anywhere". ⚠ S6's gate text ("white/albedo/vacuum realized as ι∘α∘π") is CLOSER to done than the plan implies — re-derive the delta from `sn/boundary/realizer.py` before scoping |
| **S7** energy V/V* | **YES** — no GroupAxis to hang it on | none |

**§I.10 sibling** (brief's specific question): `DiscreteMeasure.partition_by`
(`measure.py:1042-1112`, landed `f6e6002b` **2026-05-10** — long pre-audit) is a genuine
landed sibling of the Partition layer: disjoint-by-construction, mass-preserving,
label+indices+restricted-measure entries, canonical octant-partition consumer. It lives
on the MEASURE (the generator), consistent with §I.8's forgetful framing — the per-axis
partitions §I.10 assigns to `QuadratureAxis` would forget from exactly this. The middle
layer (BlockDigraph) and Schedule: no trace found in this slice (not searched exhaustively
— operator-side machinery is sibling B's scope).

## 5. What the plan does not know (post-2026-08-08 or unread, in scope)

1. **Three unread space files, all pre-dating the audit** ([M] `git log --follow` first
   commits: `radial_characteristic_space.py` 2026-07-04, `scalar_trace_space.py`
   2026-07-03, `spatial_moment_space.py` 2026-06-16; [M] `grep radial_characteristic\|
   scalar_trace\|spatial_moment` over the plan → 0 hits). Consequences for the audit's
   claims: F1 gains two more documented weak-equality sites (one measured live), F3 gains
   two more space-owned metric doctrines (G_sd = V_cell with a strictly-positive guard;
   face-AREA), and the fragmentation census gains the `_PART`/`_SPACE_NAME` ClassVar
   mint pattern (`radial_characteristic_space.py:349-354`).
2. **The two space.py factories are dead** — `angular_flux_space` / `scalar_flux_space`
   have zero production consumers ([M] grep; `tests/numerics/test_space.py` only). Live
   mints: `augmented_mesh.py:1099` and `_bases.py:381/:503`. Any Phase-S step that
   "wraps old constructors in anonymous axes" (S1) must target the LIVE mints; rewiring
   the dead factories alone would be a no-op that reads as done.
3. **FullFieldSpace is a family-blind class with a second production instance** —
   System B's ψ½ composite under `name="radial_characteristic"` (B.2b DP1;
   `full_field_space.py:180-189,229-234`). S1's "FullFieldSpace … already [is a] direct
   sum" holds twice over; the ⊕-unification target set is larger than one instance.
4. **Q6 measure-side machinery** (partly post-audit): `ReferenceMeasure`/`ExactnessClaim`
   (`exactness.py`, 2026-08-02), the selector's reference conjunct + frozen registry
   (`626dc855`, `b02dd536`, 2026-08-14), `generating_measure.py`. Overlap: §4 S2 row.
   And it minted the F1-class defect #369.
5. **The LD moment-mass metric product** at the bulk mint (`augmented_mesh.py:1095-1098`,
   #310 C2 ruling 3): `G_bulk = V·w_n ⊗ diag(moment mass)` hand-rolled inline — the
   per-factor metric structure F4/S1 wants, already needed and already spelled once,
   at a site the plan's F-rows don't cite.
6. **Post-audit space-layer delta is nil** (§0) — the plan's audit is NOT stale on its
   own files; it is incomplete on the three unread ones.

## 6. Overlap map to open issues (each read via `gh issue view N` 2026-08-19)

| issue | state | relation to this slice |
|---|---|---|
| **#369** | OPEN | F1's defect class on the measure side; fix shape = S3's interning. Cite from 0.8/S3 charter. (§2d) |
| **#330** | OPEN | The operational sibling of S1-S3: "spaces MANDATORY on every operator; the SPACE owns shape + traversal." Confirms half is built (`_AdjointOperator` reads metrics from spaces — the plan's praise P8 agrees). G6.1's `AngularFaceTraceSpace` docstring cites it (`spaces/__init__.py`, `angular_trace_space.py:520`). Phase S landing order should reconcile against #330's measured state rather than re-derive it |
| **#331** | OPEN | Domain-identity gap (displacement space) across L/S/B — F1-adjacent but flux-ontology-owned. **For sibling D** |
| **#295** | OPEN | `LayoutBearingSpace` Protocol — the "layout: FaceLayout" concept spelled in ≥3 spaces + `# type: ignore` at `_bases.py:871` ([M] still live, with getattr guard `:813`). ⚠ Its named third file `starting_direction_space.py` is STALE — renamed to `radial_characteristic_space.py`, key type now `FaceLayout[tuple[int,int]]` (`:363-365`), and the space was SPLIT in two (interior/boundary), so the issue's site list needs a refresh when picked up. Overlaps S1 (layout = index-set half of the Axis job) |
| **#297** | OPEN | Carrier-contract fragmentation (Vector / `_is_ravellable` / SystemField) — the composition-fragmentation twin one level down; a Phase-S scope decision (in or explicitly out) |

## 7. Refuted / changed premises (first-class output, with the structural reason)

1. **REFUTED (site, not substance): "monolithic shape-tuple factories (the SN bulk
   path)" as pointing at `angular_flux_space`/`scalar_flux_space`.** Reason: both
   factories have zero production consumers ([M] grep, test_space.py only); the SN bulk
   path mints a bare `FunctionSpace("sn_bulk", …)` inline at `augmented_mesh.py:1099`
   and the field layer mints per-`_SPACE_NAME` at `_bases.py:381/:503`. The monolithic-
   mint DEFECT is confirmed; a Phase-S step written against the factory sites would fix
   dead code and leave the live defect standing.
2. **CHANGED (tense): F2's two-methods-two-axes divergence is history, not a live
   defect.** Reason: fixed 2026-08-04 (pre-audit) — the ⛔ note itself says "until
   2026-08-04", and [M] the live probe shows `inner_product` ≡ `apply_metric`-pairing on
   the SH production layout at HEAD. What is live is the CONVENTION (leading-broadcast by
   position, no Axis) — which is the F2 headline, so the plan's row stands once its
   parenthetical is read as the failure-class citation it is.
3. **CHANGED (count): F3's "three answers in three places" under-counts.** Reason: two
   unread files each carry a further space-owned metric doctrine (G_sd = V_cell;
   face-AREA), and the explicit-unit `current_space` adds a deliberate fifth spelling of
   intent. The direction of the correction favours the plan's own S2 resolution (more
   spaces already own their metrics; the Euclidean-absorbed arm is the shrinking
   holdout).
4. **CHANGED (scope): S6 is not greenfield.** Reason: the factored Lambertian chain
   (π = `PartialCurrentOperator`, ι-side = `IsotropicEmissionOperator`, named S(f)
   intermediate) landed 2026-08-04 with per-link honest transposes; what remains of S6 is
   the axis-owned constructor + minted-together adjoint pair + the fission/P₀ re-spelling
   — a smaller step than the row reads. (Plan-authoring §6c: S6's gate "an ERR-039-class
   negative control REDs" should be re-derived against the realizer's CURRENT white arm
   before the step is sized.)
5. **NOT refuted but load-bearing to restate: nothing in the brief's landed-campaign
   list invalidates III-S.** Reason: G6.x and B3.2/B3.4 predate the audit (already in
   its praise); B3.3 and Q6 land outside the space files (operator.py / measure layer);
   [M] the space files themselves took two docstring-only commits since 2026-08-08.

— end (agent A, explorer; probes runnable via `.venv/bin/python`; grep commands quoted
inline; graph/workspace: main checkout, no worktree switch needed).
