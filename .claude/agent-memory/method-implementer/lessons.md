# Method Implementer — Lessons

Each lesson names the rule, skill or definition clause that does not already carry it (the workflows rule, invariant 6). Retired, each restated by a clause re-verified 2026-11-16: L-002 (this agent's definition, Done, `npx pyright <touched files>` against a `total: 0` baseline), L-003 (`retirement-audit` G.24), L-004 (`retirement-audit` A.1–A.3), L-005 (`process-discipline`, "Mutation-testing an uncommitted file"), L-007 (`vv-principles` Mode 10), L-012 (`coding-elegance` anti-pattern 20). L-010, L-011, L-014 and L-015 were retired on 2026-09-22 and are restored below with only the part no carrier holds. Numbers are not reused.


## L-003 (residual, restored 2026-11-16)

- check: after deleting a symbol, grep its module for the short name of every import the deleted body used (`spla.`, `sp.`);
  an import with no remaining use is part of the retirement (`retirement-audit` G.24).

## L-001 -- Two typing traps the skill does not list

`coding-elegance` anti-pattern 19 carries the ignore-free fixes for a reduction typed as `int` (`weights @ stacked`), an invariant typevar on a one-directional Protocol (declare the variance) and a stub obscured by a subclass override (`if TYPE_CHECKING:`). Two more:

- **`isinstance(<concrete>, <Protocol>)` in an intrinsic-property test** trips `reportGeneralTypeIssues` "overlaps unsafely" when the concrete type structurally matches. Fix: route the value through an `object`-typed helper; the runtime check is the deliberate point.
- **An attribute meaningful in only one mode** and indexed in the other: keep it the always-valid type (`axis: int`, default 0, ignored when unused), not `T | None`. A legacy attribute meaningful in one mode only becomes a `@property`, not stored state. This is a typing patch, not a design: an `axis=` on a mathematical object is itself a weld (`coding-elegance`, "A weld is the failure to answer question 1"), so first ask whether the mode should be a type.

The symptom-to-fix map for the skill's three: `reportReturnType` on `sum(<arrays>)` is the contraction; `reportInvalidTypeVarUse` on a one-directional Protocol is the variance (input `contravariant=True`, output `covariant=True`; an operator's `apply(x: V) -> V` is invariant by dual use and emits no warning); `reportRedeclaration` on a base stub is `if TYPE_CHECKING:`.

---

## L-006 -- The seam shapes a crosswalk must have a row for

The directive to write the crosswalk first is this agent's definition (Readiness item 3) and `coding-elegance` "Convention crosswalk". Four seams bit in shipped carves, and none is among the skill's seven axes. The typed-field seam differs from the skill's "packed vs typed layout": that axis is how one quantity is stored, this one is the shape contract of the typed space, which the carve must widen. The first, second and fourth generalise (`candidates.md` of the 2026-11-16 distillation proposed them for the skill); the scan seam is LD-specific:

- **×V / ÷V reductions.** The ×V form carries an extra `θ·V`; the ÷V form does not. Tabulate which quantity carries θ and which carries V, then write the single-source helper. A first-draft `V·d2` (missing θ) was caught by the d=1 link test against an independently assembled dense primitive.
- **Mass-weighted versus raw source.** The dense UBLD right-hand side is `M·S`, while `A = (L+C) − S` subtracts the scattering source raw. The matvec residual must be M-normalised, or the slope rows disagree by `θ^{|i|}`; the average row agrees by coincidence because `M_00 = 1` (a degenerate-regime false green).
- **Scalar scan versus moment scan.** The convex face blend `ψ̄ = (1−w)ψ_in + w·ψ_out` decouples from `ψ_out` only for the flat-source LD. With a slope source the parallel-prefix scan propagates the face chain through a slope-augmented `b` and reconstructs `(ψ̄, ψ̂)` per cell by the Schur complement, not by `cell_average`. The closure's validity regime is a crosswalk row.
- **Raw `.values` versus the typed field space.** A producer-side lift on raw ndarrays is cheap and byte-identical, but the iterate feeding it is a typed `Field` with a rigid `(ng, *spatial)` space. A crosswalk that maps the sweep, scattering and source plumbing but skips the field-space shape gate under-scopes the carve: the einsum accepts the wider axis and the typed accumulator cannot hold it without widening the space contract.

check: the d=1 reduction link test against an independently assembled primitive is the catcher for every row.

---

## L-008 -- A source-text sentinel is not rename-transparent

A `sha256(inspect.getsource(...))` pin reddens on any identifier rename inside the pinned body, even when the floating-point reduction tree is untouched. When a rename crosses such a sentinel, the re-hash is mandatory and cites the bit-identity proof (the strict regression gate green at the same ULP values), not an assertion of it.

---

## L-009 -- A codomain re-point's blast radius is three layers wider than the named leaves

Re-pointing operator matvec leaves from a history-bearing type (`TimedFullField`) to its timeless base (`FullField`) is value-neutral only because the iteration driver re-attaches the timed type through the carrier's `__add__` recombine; verify first that `advance` is dead on the production path. The surface is wider than the leaf set the brief names:

(a) the shared producers both operators delegate to (`loss_action`);
(b) any metric or space helper that hardcodes empty metadata on a rebuilt composite: route it through the polymorphic `_recombine` hook;
(c) every test asserting `isinstance(out, TimedFullField)` or `out.history_depth` (dozens; they break by design and are re-pointed to the timeless contract).

The resolvent `solve` is not a base arrow: it produces the driver's iterate and stays timed (covariant, Liskov-safe). Input dispatch on the carrier type (`@singledispatchmethod`) is orthogonal to the output codomain.

---

## L-010 -- A subclass promotion inherits the action and adds only the name

When a brief says "promote leaf X onto general base B", the subclass body is the name, the back-compatible constructor and the domain-specific algebra dispatch, nothing else: inherit `apply`/`solve`/`apply_transpose` and delete the copied bodies. Keep the legacy attribute surface as properties off the base's field (then: `.sigma` returning `coefficient.values`; `[R]` the example's `.sigma` property was not found in the tree on 2026-11-16, `coefficient.values` is live in `orpheus/sn/operators/streaming.py`), so readers are byte-unchanged with no duplicate storage, and single-source the capability gate by building the engine and inheriting its spectrum, not re-deriving it. Before promoting a σ-dependent action, run L-011's probes.

---

## L-011 -- Prove an operator's affinity in σ before a σ-free or leaf-sum carve

A carve that extracts pure streaming from a matvec monolithic in σ is possible without duplicating the discretisation only if the matvec is affine in σ. Two probes before the carve: `loss_action(0, ψ) == loss_action(σ, ψ) − σ·ψ`, and the candidate σ-free action evaluated at two very different σ returns the same field (the decisive σ-leak test). The σ-free primitive is then `loss_action(0, ψ)` named as its own action, and the production composite stays byte-identical. Today's spellings: `loss_action(self, sigma, psi)` and `streaming_action(self, psi)` in `orpheus/sn/loss_representation/__init__.py` (`[M]` 2026-11-16, signatures read; the probes were not re-run). The teeth for the override-versus-leaf-sum twin are `vv-principles` anti-pattern 17(j). A twin that is value-correct by coincidence shipped no wrong value, so it is a `foundation` gate with `verifies(...)`, not a `catches(ERR)`.

---

## L-013 -- The discrete curvilinear WDD/M-M operator is NOT additively separable

A brief asked for the discrete curvilinear matvec as three additive peer leaves `(L_spatial + L_angular + C).apply(ψ)`, forbidding the subtractive form. The continuous operator splits (`Ω·∇ψ + (1−μ²)/r ∂ψ/∂μ`), but the discrete WDD/M-M closure couples the spatial sweep and the angular redistribution in the cell balance; they do not commute as separate matvec leaves, and the landed form was the coupled (subtractive) operator. How to apply: when a brief mandates an additive-leaf decomposition of a curvilinear discrete operator, push back; do not contort the discretisation to fit a continuous-operator algebra the discrete form does not honour (the generic return is this agent's definition, "when the mathematics refuses a structure the plan prescribes").

`[R]` possibly refuted, owed a re-measure before anyone cites it (#502): the redistribution operator was later measured to factor as `R = R_spatial ⊗ A_angular` (`.claude/plans/streaming_path_says_what_it_is.md`, finding F1), and that plan's option "remove the angular obligation from the spatial protocol" treats the two axes as separable. A tensor-product factor of one term is not by itself additive separability of the cell balance; measure which holds before relying on this lesson.

---

## L-014 -- A single-DOF scheme's scalar geometry constant is the corner of a Gram matrix

When widening a scheme from one degree of freedom per cell to N (DD to LD, a scalar to a moment vector), every shipped scalar geometry constant is the `(0,0)` entry of some Gram matrix. Do not widen by multiplying the scalar into each moment row: derive the weak form and read off which inner product, under which measure, each constant belongs to. The scalar says nothing about the off-diagonal, which can be load-bearing; and the new scheme's `p = 0` truncation is often not the shipped scheme (LD truncates to STEP, not DD). check: before a moment-widening carve, tabulate the measure each term carries in the weighted-residual reduction; the table is the crosswalk. The curvilinear instance, with its numbers: `ld_curvilinear_shape.md`, rulings 1, 2 and 7.

---

## L-015 -- Decide a design fork's visibility to the gate by one substitution

`vv-principles` Mode 12 carries the principle (enumerate the stabiliser at gate-design time, the fixture side included). The recipe it does not state: for a design fork, compute `candidate_A − candidate_B` symbolically and substitute the intended gate's fixture. If it vanishes, the gate cannot rank the fork; say so in the plan and name a fixture that does not annihilate it. Two instances: `ld_curvilinear_shape.md`, rulings 2 and 5.

---

## L-016 -- A quotient has two realizations; decide which one a predicate speaks before cataloguing an entry

An orbit space `M/H` is realized two ways:

* **the invariant chart's codomain**: the image `p(M)` of the orbit map, in invariant coordinates; canonical, and what the Procesi-Schwarz procedure emits;
* **a fundamental domain**: a subset of `M` meeting each orbit once, in the base's ambient coordinates; what a folded quadrature ships.

The group decides which one ships:

| | positive-dimensional `H` (e.g. `SO(2)`) | finite `H` (e.g. a mirror) |
|---|---|---|
| `dim(M/H)` | drops | unchanged |
| chart as data | a real reduction (3 floats to 1) | no reduction (3 to 2, third recoverable) |
| canonical section? | no (any half-meridian) | yes for a reflection (the fixed half-space, closed half strict) |
| what ships | the chart | the section |

Three traps, measured on `S^2/<σ_y>` (`[M]` 2026-08-31, `scratch/sigma_y_orbit_derivation.md`, `scratch/probe_sigma_y_orbit.py`):

1. **A chart is blind to any defect in the coordinate it drops** (`vv-principles` Mode 12, applied to a projection). The `(x,y,z) -> (x,z)` chart admits the ERR-080 forgery `(μ,0,0)`: `μ² + 0 ≤ 1` is a legal disk point. check: intersect the chart's kernel with the threat model before adopting chart coordinates for a membership check.
2. **The first catalogued entry's accidental shape becomes the field's type, and the second entry is where you find out.** `singular_stratum: tuple[float, ...]` fitted `S^2/SO(2)`'s two poles and could not hold `σ_y`'s circle; the general spelling is the vanishing locus of `det_gram`. [REMEDIED: `[M]` 2026-11-16, the field now holds that SymPy expression, `orpheus/geometry/manifold.py:1383`; the repairing commit was not identified.]
3. **The producer may already have decided.** `DiscreteMeasure.quotient` (`orpheus/numerics/measure.py:1318`; `[M]` 2026-11-16 by AST) keeps `nodes[representative]`, a selection of parent nodes, so every measure derived that way is a section in the ambient coordinates. Read the producer before ruling on the type.

The counting rider (not carried for code work: `plan-authoring` is path-scoped to `.claude/plans/**`): a minimality count needs its predicate written down. Counting decomposables with `k >= 1` includes the generators themselves and reported "0 new generators in every degree", a clean table saying the empty set generates the ring; `k >= 2` is the predicate.

check: before adding an orbit-space catalogue entry, state which realization the membership predicate speaks and compare it with the nodes the tree's own producer emits. A disagreement is a data-model ruling for the user, not a choice inside the entry.
