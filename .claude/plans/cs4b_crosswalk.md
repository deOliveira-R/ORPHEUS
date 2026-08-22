# CS4b convention crosswalk — fields are space elements

**Status: WRITTEN 2026-08-22, pre-code (L17: the crosswalk IS the architecture;
the code is its transcription).** Sources: `.claude/plans/cs4b_fields_design.md`
(rounds 0–3c, the ruling table) + `scratch/cs4b_verification_plan.md` (S1–S7,
gates, the 22-arm table). Every `[M]` below is theirs — cited, not re-measured;
`[M]` `git log 34df88cb..HEAD -- orpheus/ tests/` is EMPTY at write time, so
the plan's anchors hold at HEAD. Boundaries crossed: numerics (space/axis/
frame/field) ↔ transport fields ↔ sn (mesh/operators/solver) ↔ diffusion ↔
homogeneous (D5 wall only).

Per L18 every bridge lands at the PRODUCER. Each section: who produces, who
consumes, where the bridge lives, and the convention pinned in one line.

---

## B1 — space provenance: per-leaf tag mint → carrier-cached axis product

| subsystem | input | internal | output |
|---|---|---|---|
| today's producer `cls._space_for_mesh` (`_bases.py:361,483`) | mesh | `FunctionSpace(name=_SPACE_NAME, shape=…)` — fresh mint per call, NO metric, name = role | uncached tag space; `space ==` twin-mints, `is` fails |
| new producer: the CARRIER | its own axes | `MaterialMesh.bulk_space` (ships, cached): `of_axes(EnergyAxis.from_materials(reachable), Axis("spatial", shape, weights=V, NODAL))`; NEW `SNMesh.angular_bulk_space` (cached): `of_axes(Axis("angular", (N,), weights=quad.weights, NODAL), *bulk_space.axes)` | ONE `is`-stable instance per carrier; `==` by structural content |
| consumers: field factories → (S4) the primary ctor | the carrier | read, never mint | fields hold the space |

**Bridge (S1 mint / S2 re-point):** the mint moves to the carrier — the
definition site. Conventions pinned:

1. **Axis order is `(angular, energy, spatial[, moment])`**, matching the bulk
   tensor `(N, ng, *spatial[, 2^d])` and matching `full_field_space`'s dense
   interior layout. The scalar bulk is the SAME product minus axis 0 — so
   retract-along-angular is literally "drop axis 0", and the angular mint
   REUSES `bulk_space.axes` verbatim (energy-arm single source: the
   `EnergyAxis.from_materials` rule is spelled once, in `bulk_space`; the
   angular mint must NOT respell it). Battery arm M3 (order swap) is the
   control.
2. **The derived space NAME is never pinned** (verif. plan R4: CS2's axis
   subclasses change the digest). Gates pin per-axis
   `label/shape/kind/weights` + relative `is`/`==` — never the literal.
3. **Role leaves the space.** `_SPACE_NAME` retires; the three role spaces
   (`angular_flux`/`angular_source_sink`/`angular_residual`) collapse to ONE
   cached object per family. Consequence bridged at B3 (class gate becomes the
   sole role enforcement, G2.3).

## B2 — the inner product: Euclidean (metric-less) → physical `G = V·w_n[⊗mass]`

| subsystem | input | internal | output |
|---|---|---|---|
| today: bare tag space | values | `inner_product_weights is None` → flat Euclidean | `AngularFlux.l2` = 4.99286… on the §2 fixture |
| after: axis-built space | values | per-axis weights compose the Gram | `l2` = 2.95933… — the **41 %** move, `[M]` verif. §2 |
| consumer `iteration.py:740` | `leaf.l2` | `increment_norms` → ρ trajectory, `true_error_estimate` | diagnostics move; **SI STOP immune** (`[M]` rides `_l2_norm` on bare ndarrays, `:709`) |
| consumer `angular_residual.py:206` | `self.l2 / q_norm` | `relative_to` | the portable residual criterion moves |

**Bridge (S2):** the metric lives ON the space — no consumer ever re-applies
weights (Pattern 7). The move is owned, not incidental: G-M1 (re-derive the
CS3 trajectory against the ρ≈c Adams–Larsen anchor, never old-vs-new), G-M2
(edit the gate's "CS2 owns" block in place), G-M3 (the in-file control must
still red). DriftWarning wall runs at S2 end (R5); a red re-baselines under
the vv three criteria only. D5 stays 8/8 throughout (the homogeneous wall).

## B3 — partner identity: `mesh is` (provenance) → space CONTENT `==` (cached-`is` fast path)

| subsystem | input | internal | output |
|---|---|---|---|
| today: 22 guard sites | two fields | `a.mesh is b.mesh` | twin mesh REFUSED (provenance gate) |
| after (F2 RATIFIED) | two fields | `a.space == b.space` (axis-built structural equality; cached `is` as fast path) | discrimination by CONTENT — §4's table |

**The discrimination convention, per family (`[M]` verif. §4):** volumes / ng
/ quadrature discriminate; **BC never does** (a law changes neither DOFs nor
Gram — laws are operator data; every space law-blind, correctly); the SCALAR
family is additionally quadrature-blind (ratified permission, G2.2 — S4-φ +
S8-φ becomes legal, docstring carries the DOF+Gram rationale).

**Bridges (S3), one per convention consequence:**

- **All 22 arms re-key in S3** with per-arm witnesses whose fixtures differ in
  something the SPACE can see; the 8 unwitnessed arms (O4 O5 O8 O11 O12 O13
  O17 + R1's `sn/solver.py:338`) get MANUFACTURED witnesses in the same step
  — a rewrite of an unwitnessed guard is indistinguishable from its deletion.
- **Role enforcement moves wholly to the CLASS arm** (`Field._check_partner`'s
  `TypeError`, which today fires before the space arm ever sees role): G2.3
  pins it, M5 proves the coverage is NEW (today masked).
- **Composite identity compares BLOCKS, never the composite `==`** —
  `FullFieldSpace.__eq__` is `(name, shape)` with a literal name, i.e.
  block-blind (R2). `Composite.__post_init__` re-keys per-block
  (`interior.space is/== …` + trace), and the `getattr(x,"mesh",None)`
  tolerance dies in S3 (it becomes a silent no-op at S4 otherwise —
  vv #28's temporal twin).
- **Cross-problem discipline** stays at the iteration layer (CS3 relocation);
  the space gate carries arithmetic well-definedness only.

## B4 — the Σw boundary: one "embed" → the frame square's two named up-arrows

The L17 classic (`/W` normalisation). `[M]` verif. §3/§9: the tree ships TWO
up-maps a factor `Σw` apart, both production, plus one Σw-less broadcast with
zero callers:

| spelling today | what it is | convention |
|---|---|---|
| `AngularField._integrate_angular_values` (`_bases.py:409`) ≡ `HarmonicFrame.analyse`@ℓ=0 | **R** — retract, the measure contraction | down; three bit-equal spellings |
| `HarmonicFrame.reconstruct`@ℓ=0 | **R.H** — the SH synthesis (primal up) | NO `÷Σw`; `R∘R.H = Σw·id` |
| `from_isotropic` = `broadcast(φ/Σw)`; `dsa.py:661` hand `÷Σw` | **E** — the section (dual up) | `÷Σw`; `R∘E = id` BIT-EXACT |
| `ScalarSourceSink.as_per_ordinate` (`:198`) | a Σw-LESS broadcast | ⛔ 0 production callers — must not survive as a third normalisation |

**Bridge (S6, Q5 = complete the frame square):** the frame gains its fourth
face `dual_reconstruction = reconstruction ∘ gram.apply_inverse_metric` — the
`÷Σw` is **G₀⁻¹, DERIVED**, never hand-spelled again. The space verbs are its
consumers: `space.retraction(axis) → R`, `space.embedding(axis) → E` (defined
by `R∘E = id`), two DIFFERENT named arrows so the ERR-051 factor-`Σw` class is
unspellable — G6.4 (`R.H == Σw·E`) is the anti-ERR-051 gate, M29/M30 the
convention-partition arms, G6.8 the discrimination row. Every retract/embed
gate row states its Σw convention IN THE ASSERTION.

## B5 — the moment tail: `SpatialMomentSpace` via `*` (mass-less) → a MODAL axis carrying the scheme's mass

| subsystem | input | internal | output |
|---|---|---|---|
| today's field tail (`_compose_spatial_moments`, `_bases.py:216`) | per_axis (caller-selected, default 1) | `space * SpatialMomentSpace.from_per_axis(…)` — the `*` densifier | `axes = None` (axis structure KILLED), no mass anywhere (field spaces have no metric) |
| today's interior tail (`full_field_space`, `augmented_mesh.py:1096-1101`) | `self.scheme` | dense `G_bulk = V·w_n ⊗ moment_mass_diagonal` | scheme-driven, mass-correct |
| after (G1.4/R9 measured form) | scheme | a **MODAL** axis, `weights = scheme.moment_mass_diagonal(ndim)` in the `of_axes` product | Gram ≡ dense interior at ≤1 ULP (`[M]` R9); `has_coordinate_cone → False` |

**Bridge (S1 mint, S2 re-point):** the moment-axis mint is **scheme-owned**
(the mass IS the scheme's — `θ` enters `moment_mass_diagonal`, so
basis↔mass single-sources at the scheme; `full_field_space`'s interior and the
field tail must consume the SAME axis object). Conventions pinned:

1. Construct-general/select-narrow SURVIVES: the caller still selects the
   width (default 1 = no factor, byte-identical); the scheme is never read by
   default at the field layer.
2. **Cone routing (Q6 RATIFIED):** the MODAL kind flips LD bulk fields'
   `cone_violations` from the `None`-legacy arm to the honest typed REFUSAL
   (`field.py:471-479` ships the branch since CS1). DD/Step (all-NODAL) →
   `True` → answers. Two-sided gate + the `field.py:463-475` /
   `scheme.py:530` re-words ride S6/S7; the exact vertex test is **#400**.
3. **Stored-twin heal is BOTH or NEITHER** (R6): `L` and `spatial_moments`
   both become space reads (`space.find_factor(SphericalHarmonicSpace).L`;
   the tail read off the axis) — a partial heal re-opens the desync the
   guarded cache makes unrepresentable today. #399's three members become
   space-derived (`replace(...)` + space-side truncate), never a re-threaded
   stored width.
4. The HARMONIC family's `SphericalHarmonicSpace * cell_group` stays the `*`
   densifier this phase — its `has_coordinate_cone = None` is an ACCIDENT of
   `*`, not a ruling; docstring says so (CS2 axis-ifies it to MODAL/False).

⚠ **Open corner (S2, steer with the user):** a width>1 mint on a carrier with
NO scheme (bare `MaterialMesh`) has no mass source — production never does
this (`[M]` every width>1 origin site reads `sn_mesh.scheme` or propagates a
field's own width), so the honest options are refuse vs counting-weight;
decide at S2, do not discover. Same corner: a requested width ≠ the scheme's
has no shipped mass — today it silently builds a mass-less space of the
requested width; under the axis form it must refuse.

## B6 — construction: mesh-primary factories → the space-primary ctor

| subsystem | input | internal | output |
|---|---|---|---|
| today | `(values, mesh)` | `from_mesh`/`zeros_on`/`from_ndarray` mint the space, pass mesh through; composite `zeros(interior=…, boundary=…, mesh=…)`; `_from_balance` → `cls.from_mesh` (2 `type: ignore`) | field holds `(values, space, mesh)` |
| after | `(values, space)` | `cls(values=…, space=…)`; `Field.zeros(space)`; `_from_balance` → `cls(values=lhs.values − rhs.values, space=lhs.space)` | field IS `(values, space)`; the two ignores retire (ratchet may DECREASE — itself a re-baseline) |

**Bridge order is the §7 mechanics — bodies → consumers → signatures** (a
step-boundary argument, not effort): S2 flips what factories BUILD (signatures
unchanged, zero call sites move, the B2 move isolated); S3 empties `.mesh` of
consumers; S4 deletes the field from the 2 declaration roots (headline gate:
`CrossSectionField(values=…, space=…)` CONSTRUCTS — `[M]` today `TypeError`);
S5 migrates the 909-site sugar spelling (gate = grep/AST count); S6 re-homes
the math-bearing factories (`from_isotropic` → E's consumer;
`from_mesh_and_L` → space construction; `from_face_arrays` → layout assembly).
`MomentResidual` (Q1): the flip dissolves the 2-arg blocker; record-as-choice
at `_from_balance`, mint only with a consumer.

## B7 — the mesh's remaining consumers: field-attribute read → carrier-held

`[M]` 182 production + 116 test non-`self` `.mesh` reads (verif. §7). After
S4 a field cannot answer "which mesh"; the consumer that needs the MESH (not
the space) reads it off the operator/solver that holds it — the #226 F2
pattern. The space answers every STRUCTURAL question (sizing, measure,
membership); `Composite.mesh` (`full_field.py:279`) is re-derived or retired
in S3/S4, and the AST predicate (not text grep — `sn_mesh.mesh` etc. are
legitimate survivors) is the done-when instrument.

## B8 — sentinels and refusals (S7)

| site | today | after |
|---|---|---|
| `sn/mesh/augmented_mesh.py:322` | bare `assert` — messageless plain; **`-O`: deep `AttributeError`** (the canonical runner never runs it) | a real `raise` modeled on `diffusion/augmented_mesh.py:211-218`; gate runs under `-O` (G7.1) |
| `material_mesh.py:207,492` `mesh is None` | TWO meanings welded (d≥3 axis-native vs meshless quotient) | un-welded — `ndim` DISCRIMINATES them (G7.3); named states, not a sentinel |
| `material_mesh.py:517-523` `.areas` | 3 arms share one message, true for 1 | per-arm messages, shortest-fragment pinned (G7.2) |
| `isotropic_scattering.py:203` bare assert | inert under `-O` | triage: type-narrowing (fine) vs admission (→ raise) — one line, R7 |

## Bridge summary — every bridge at the producer (L18)

| # | convention | bridge lives at | step |
|---|---|---|---|
| B1 | space mint + axis order + name discipline | the CARRIER's cached properties | S1/S2 |
| B2 | the physical metric | the SPACE (per-axis weights) | S2 |
| B3 | identity = content; role = class; composite = per-block | the 22 guard sites, re-keyed with witnesses | S3 |
| B4 | Σw: R vs R.H vs E, `÷Σw` = G₀⁻¹ derived | the FRAME (`dual_reconstruction`) + the space verbs | S6 |
| B5 | moment mass; MODAL kind; cone routing | the SCHEME's axis mint | S1/S2 (routing S6/S7) |
| B6 | space-primary construction | the primary ctor; factories retire/re-home | S2→S5, S6 |
| B7 | mesh access | the operator/solver that holds it | S3 |
| B8 | typed refusals under `-O` | each repair site | S7 |

**Steering points reserved for the user (surgical posture):** the B5 open
corner (no-scheme / wrong-width mass source, at S2); Q7's `_StubMesh` seam
(default: `MaterialXSField` takes the space narrowly, at S2/S3); the S2-end
DriftWarning wall verdict if it reds (R5); G-M1's trajectory re-derivation
numbers (S2).
