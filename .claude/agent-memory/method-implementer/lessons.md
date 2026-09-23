# Method Implementer — Lessons

Lessons each name the rule, skill or definition clause they are not already covered by (the workflows rule, invariant 6). Retired 2026-09-22 by the agent-definitions audit (`.claude/plans/archive/agent_definitions.md`, §8 step 3): L-002, L-003, L-004, L-005, L-007, L-010, L-011, L-012, L-014 and L-015, each restated by a clause of `retirement-audit`, `vv-principles`, `coding-elegance`, `process-discipline` or this agent's definition, or duplicated in `ld_curvilinear_shape.md`; numbers are not reused.

## L-001 -- Principled typing beats suppression, every time

The recurring temptation under a pyright-net-new-0 gate is `# type: ignore`.
It is almost never the right fix; nearly every pyright red I hit had a
*principled* spelling that ALSO read like the math. The catalogue of
traps and their ignore-free fixes:

- **`return sum(generator-of-arrays)` types as `int`** (`reportReturnType`).
  Fix: spell the reduction as the linear-algebra contraction it IS —
  a convex average is `weights @ stacked` (a matvec), not `sum(...)`,
  not `np.sum` (wrong axis). The stacked matrix is the named
  intermediate (Pattern 3). The math reads `Σ_i w_i χ_i`.
- **A one-directional Protocol (`evaluate(x:V)->R`) reusing an invariant
  `V`** reddens `reportInvalidTypeVarUse`. Fix: declare properly-variant
  typevars — input `contravariant=True`, output `covariant=True`. (An
  operator's `apply(x:V)->V` is invariant-by-dual-use and emits NO
  warning — the asymmetry is real, mirror it.)
- **A base abstract method that subclasses override + the base also has a
  concrete method calling it** → a runtime `raise NotImplementedError`
  stub reddens `reportRedeclaration` ("obscured by subclass"). Fix:
  declare the signature under `if TYPE_CHECKING:`; pyright is satisfied,
  the runtime override stands.
- **`isinstance(<concrete>, <Protocol>)` in an intrinsic-property test**
  trips `reportGeneralTypeIssues` "overlaps unsafely" when the concrete
  structurally matches. Fix: route through an `object`-typed helper (the
  runtime check is the deliberate point).
- **A mode-discriminated attribute** indexed in the OTHER mode: keep it
  the always-valid type (`axis: int`, default 0, ignored when unused),
  NOT `T | None`. A legacy attribute meaningful in only one mode → a
  `@property`, not stored state.

How to apply: when pyright reddens on a NEW symbol, find the principled
spelling FIRST (variance, TYPE_CHECKING, the contraction it actually is).
Reach for `# type: ignore` only after proving no such spelling exists —
which, in practice, has not yet happened. (#257 S5/S6/S8b/S4/S3a/S10b.)

---

## L-006 -- Write the convention crosswalk BEFORE the carve — it is the teeth, not ceremony

→ The GENERALIZED directive ("write the crosswalk before any
convention-crossing carve") is now a standing step in AGENT.md
(procedural workflow step 3b). The four forensic instances below stay
here as the recall catalogue of *which seam shapes* bite.

Every multi-convention carve I shipped had its load-bearing bug exactly
at the convention seam, and the crosswalk row (had I written it first)
would have caught it at write-time. The instances, all the same shape:

- **×V / ÷V / ×V-scan reduction conventions:** the ×V form carries an
  extra `θ·V`; the ÷V form does not. Write the "which quantity carries θ,
  which carries V" table FIRST, then the single-source helper. A first-draft
  `V·d2` (missing θ) was caught by the d=1 link-vs-dense-primitive test —
  the link proof is the teeth.
- **A "unify the matvec" carve** inherits a source-normalization mismatch:
  the dense UBLD RHS is mass-weighted (`M·S⃗`) but `A=(L+C)−S` subtracts the
  scattering source RAW at the OperatorSum. The matvec residual must be
  M-normalized or the slope rows disagree by `θ^{|i|}` (the AVERAGE row
  coincidentally agrees because `M_00=1` — a degenerate-regime false green).
- **A moment-carrying parallel-prefix scan** is NOT a drop-in widening of
  the scalar scan: the convex face-blend `ψ̄=(1−w)ψ_in+w·ψ_out` decouples
  from ψ_out ONLY for the flat-source LD. With a slope source the scan must
  propagate the FACE chain via a slope-augmented `b` and reconstruct `(ψ̄,ψ̂)`
  per cell via the Schur — NOT `cell_average`. Map the closure-validity
  regime as an explicit crosswalk row.
- **The TYPED-FIELD space layer is a crosswalk row of its own.** A producer-side
  lift on RAW `.values` ndarrays (below the field-space gate) is cheap and
  byte-identical, but the iterate that FEEDS it lives ABOVE the gate (typed
  `Field` with a rigid `(ng,*spatial)` space). A crosswalk that maps sweep/
  scattering/source plumbing but skips the field-space (Pattern-4 shape) layer
  under-scopes the carve: the einsum accepts the wider axis, the typed
  accumulator cannot hold it without a space-contract widening.

How to apply: for any carve crossing ×V/÷V, mass-weighted/raw, scalar/moment,
or typed/raw boundaries, write the crosswalk table to the plan FIRST. The
load-bearing rows are source-normalization and the field-space shape gate.
The d=1-reduction link test against an independently-assembled primitive IS
the catcher. (#240 D5b S1/unified-matvec/owed2-scan/S3-A; coding-elegance Pattern 7.)

---

## L-008 -- A source-text sentinel (sha256 of getsource) is bit-identity-load-bearing but NOT rename-transparent

A `sha256(inspect.getsource(...))` pin reddens on ANY identifier rename
inside the pinned body, even when the FP reduction tree is provably untouched.

How to apply: when a rename crosses such a sentinel, the re-hash is
mandatory AND must cite the bit-identity proof (the strict regression gate
stayed green at the same within-tol ULP values), not merely assert it.
(#257/#241 reaction_xs rename.)

---

## L-009 -- A codomain re-point is value-neutral only because the driver re-attaches — and its blast radius is 3 layers wider than the named leaves

Re-pointing operator matvec leaves from a history-bearing comonad to its
timeless base is VALUE-NEUTRAL only because the iteration driver re-attaches
the timed type via the carrier's `__add__` recombine (no `advance` call —
verify `advance` is dead in production first; it is, for steady-state SI/Krylov).
The blast radius is THREE layers wider than the named leaf set:

(a) the shared producers that BOTH operators delegate to (`loss_action`);
(b) any metric/space helper that hardcodes empty metadata on a rebuilt
    composite — route it through the polymorphic `_recombine` hook instead;
(c) EVERY test asserting `isinstance(out, TimedType)` / `out.history_depth`
    (dozens — they break by design; re-point them to the timeless contract).

The resolvent `solve` is NOT a base arrow — it produces the driver's iterate
and STAYS timed (covariant, Liskov-safe). Input-dispatch fibration
(`@singledispatchmethod` on carrier TYPE) is ORTHOGONAL to the OUTPUT codomain.

How to apply: before a codomain re-point, enumerate the shared producers, the
metadata-hardcoding helpers, and the isinstance-assertion tests — those are the
real surface, not the leaf set the brief names. (#257 S8a/S8b.)

---

## L-013 -- The discrete curvilinear WDD/M-M operator is NOT additively separable

A brief asked for the discrete curvilinear matvec to be positively constructed
from three additive peer-leaves `(L_spatial + L_angular + C).apply(ψ)`, forbidding
the subtractive form. This is impossible: the continuous operator splits cleanly
(`Ω·∇ψ + (1−μ²)/r ∂ψ/∂μ`), but the DISCRETE WDD/M-M closure couples the spatial
sweep and the angular redistribution in the cell balance — they do not commute as
separate matvec leaves. The honest discrete form is the coupled (subtractive)
operator.

How to apply: when a brief mandates an additive-leaf decomposition of a
curvilinear discrete operator, push back — the discrete closure is coupled. Do
NOT contort the discretization to fit a continuous-operator algebra that the
discrete form does not honor. (#197 typed-cylinder; landed via the subtractive form.)

---

## L-016 -- A quotient has TWO realizations; decide which one a predicate speaks BEFORE cataloguing an entry

An orbit space `M/H` is realizable two ways, and they are different objects
that a single `realization`-style field cannot serve at once:

* **the invariant CHART's codomain** -- the image `p(M)` of the orbit map,
  living in the invariant coordinates. Canonical (derived from the invariants),
  and the thing the Procesi-Schwarz procedure actually emits.
* **a fundamental DOMAIN** -- a subset of `M` meeting each orbit once, living
  in the base's AMBIENT coordinates. This is what a "folded" quadrature ships.

Which one the tree ships is not a style choice; it is forced by the group:

| | positive-dimensional `H` (e.g. `SO(2)`) | finite `H` (e.g. a mirror) |
|---|---|---|
| `dim(M/H)` | **drops** | **unchanged** |
| chart as data | a real reduction (3 floats -> 1) | no reduction (3 -> 2, third recoverable) |
| canonical section? | **no** (any half-meridian; the choice is arbitrary) | **yes** for a REFLECTION (the mirror fixes the half-space, and being fixed pointwise makes the CLOSED half strict) |
| what ships | the chart | the section |

⭐ **Three traps, all measured on `S^2/<sigma_y>` (2026-08-31, memo
`scratch/sigma_y_orbit_derivation.md`, probe `scratch/probe_sigma_y_orbit.py`):**

1. **A chart is Mode-12 BLIND to any defect in the coordinate it drops.** The
   `(x,y,z) -> (x,z)` chart admits the ERR-080 forgery `(mu,0,0)` after
   projection -- `mu^2 + 0 <= 1` is a *legal* disk point, the orbit of
   `(mu, +-sqrt(1-mu^2), 0)`. Not a tolerance question: the projection's kernel
   *is* the error class. ⟹ at design time, intersect the chart's kernel with
   the threat model before adopting chart coordinates for a membership check.
2. **The FIRST catalogued entry's accidental shape becomes the field's TYPE,
   and the second entry is where you find out.** `singular_stratum:
   tuple[float, ...]` fits `S^2/SO(2)`'s two poles and cannot hold `sigma_y`'s
   stratum, which is a CIRCLE. The general spelling is already present and
   un-twinned: the stratum is the vanishing locus of `det_gram` = the
   realization's topological boundary, in both entries.
3. **The PRODUCER may already have decided, silently.** `DiscreteMeasure.
   quotient` does `nodes[representative]` -- a selection of parent nodes, no
   chart -- so *every* measure the tree derives that way is a section in the
   base's ambient, by construction. Read the producer before ruling on the
   type; it is cheaper than the argument.

⚠ And the counting rider that cost a wrong first answer: a *minimality* count
needs its predicate written down. Counting decomposables with `k >= 1` (which
includes the generators themselves) reported "0 new generators in EVERY degree"
-- i.e. that the empty set generates the invariant ring -- as a clean,
internally consistent table. `k >= 2` is the predicate. Same shape as
`plan-authoring` §2's FILTER clause, on an algebra count.

How to apply: before adding any orbit-space catalogue entry, state which
realization the membership predicate speaks and check it against the nodes the
tree's own producer emits. If the two disagree, that is a data-model ruling for
the user, not a choice to make inside the entry.
