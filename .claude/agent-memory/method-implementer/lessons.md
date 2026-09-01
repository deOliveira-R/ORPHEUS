# Method Implementer — Lessons

Behavioral corrections only: "what mistake did I make, and what did I
learn that improved my behaviour?" The HOW lives in the preloaded
skills (`algebra-of-record`, `vv-principles`, `coding-elegance`,
`numerical-bug-signatures`). Campaign play-by-play lives in the topic
files this index points to — NEVER duplicate it here.

The cross-cutting meta-pattern behind half of these — **not done at
"value right" but at "value right AND type-checker clean by construction
AND retirement left no orphan AND a mutation reddens the gate"** — is
now a standing directive in AGENT.md ("Definition of done"). Each lesson
below is one *forensic instance* of that standard: the specific trap, how
it hid, and the principled fix. The generalized rule loads every session
from AGENT.md; these stay as the recall catalogue.

---

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

## L-002 -- Pyright scope traps: promoting code into a scanned dir imports its debt

Two ways "net-new 0" silently breaks:

- **Promoting a diagnostic from an UNSCANNED dir** (`derivations/diagnostics/`)
  **into a SCANNED test dir** (`tests/`) imports its pyright debt — a
  `**dict`-splat into a typed signature is invisible noise in the
  diagnostic but +N net-new once scanned.
- **The aggregate can read 0 while a LOCAL +2 hides behind −2 coincidental
  resolutions elsewhere.** Measure pyright at the SEAM file, not just the
  repo total: a re-inference that cleans up 2 unrelated downstream errors
  can mask 2 genuine new ones at the composition site.

How to apply: when relocating code across the pyright scan boundary,
run pyright on the TOUCHED FILE before and after, not just the repo
count. (#257 S10a/S4.)

---

## L-003 -- Retiring a symbol orphans its sole-use imports — grep after the delete

Consolidating N inline closures while retiring a flat-only public twin
orphans the twin's import: the twin was the only construction site for
that import in the module once the live boundary moves to the caller.
The canonical instance: deleting the flat scipy `LinearOperator` twin
orphaned `import scipy.sparse.linalg as spla` — nothing left used `spla.`.

How to apply: after deleting a symbol, grep the module for the
short-name of every import the deleted body used (`spla.`, `sp.`, the
helper name). A pyright `reportUnusedImport` will also catch it, but the
grep is the habit that makes it automatic. (#257 S7.)

---

## L-004 -- A retirement's blast radius = graph callers AND text grep AND direct constructors

`mcp__nexus__callers` on a retired internal method finds the *graph*
callers but misses two consumer classes:

- **A bypass-trick consumer** that uses the orphan via its CLASS NAME for
  some side purpose (a σ_t=0 validation-bypass in a different test file).
  `callers` on the method caught only the operator-side caller; a repo-wide
  grep for the `_ClassName` surfaced the bypass consumer.
- **DIRECT constructors of a guarded type.** A guard-at-the-data-source
  change has blast radius = EVERY direct `Mixture(...)`/`make_mixture(...)`
  caller, NOT just the factory path (`get_mixture` callers). The
  factory-path audit is necessary but not sufficient.

How to apply: for ANY retirement or guard-at-source, run BOTH graph
callers AND a text grep for the class/method name AND (for a guarded
type) an audit of direct constructors. Three searches, not one. (#238, #257 S10a; reinforces qa/numerics L20.)

---

## L-005 -- Mutation-testing an uncommitted tracked file: monkeypatch in-process, NEVER git checkout

Reverting a mutation in a TRACKED file via `git checkout -- <file>`
wipes uncommitted on-disk WIP back to HEAD — not back to your pre-mutation
state. I lost three edits this way and had to re-apply them.

How to apply: to mutation-test a file you have uncommitted edits in,
monkeypatch in-process (the cleanest) or copy to a tmp path and mutate
the copy. Reserve `git checkout` for files you have NOT touched. (#257 S6.)

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

## L-007 -- A sub-floor term (Mode-10) is closed by STRUCTURAL teeth, not a tightened value band

When a term is activated-but-unconstrained (its code path runs but its
error is O(h²)-small, below the convergence floor), tightening the
converged-flux value band does NOT pin it — the slope-source error is
sub-floor by construction. The complete closure is TWO O(1) structural teeth:

1. Assert the production producer threads the projected moment through at
   MACHINE PRECISION (the production-change proof — catches a regression to
   zeroing), with a structurally-independent (leggauss-only) reference.
2. Assert a CONSUMED source-row sign flip moves the converged answer ≫ solver
   tol (the consumption proof — catches sign-blindness), paired with the FLAT
   no-op leg (scalar/zeroed input → byte-identical) that pins the asymmetry.

A SHARPER case has NO value-improvement leg at all: a boundary-trace slope is
sub-floor for ANY value claim everywhere (no regime makes it the dominant
forcing — a correctly-consumed slope can even make the converged value slightly
WORSE). Do NOT manufacture a value-improvement leg there — it would falsely
RED a correct term. Structural teeth + no-op control is the whole proof.

How to apply: when a term is localized/higher-order-small, do not chase a
value-band gate. Ship producer-threading-at-machine-precision + consumed-flip-≫-tol
+ a no-op control. (#247 Leg A, #251 Leg B, #257 S9; vv Mode-10.)

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

## L-010 -- A subclass-promotion shares the action, adds only the name

Promoting an existing named leaf onto a new general base via a THIN subclass:
inherit the action (apply/solve/apply_transpose) from the base — DELETE the
duplicated bodies (the multiply lives once, Pattern 2); the subclass adds ONLY
the name + back-compat constructor + domain-specific algebra dispatch; keep the
legacy attribute surface alive as PROPERTIES off the base's field
(`.sigma` → `coefficient.values`) so downstream readers are byte-unchanged with
NO duplicate storage; single-source the capability gate by BUILDING the engine
and inheriting its honest spectrum, not re-deriving it.

How to apply: when a brief says "promote X onto base B," the subclass body
should be name + constructor + dispatch ONLY. If you find yourself copying an
action body, stop — inherit it. (#257 S3b; verify affinity empirically first
per L-006/L-011 when the action is σ-dependent.)

---

## L-011 -- Verify operator affinity EMPIRICALLY before a σ-free / leaf-sum carve

A "make L compute pure streaming σ-free" carve over a monolithic-in-σ matvec
is achievable WITHOUT duplicating the discretization IF the matvec is affine in
σ — but prove the affinity with two probes BEFORE the carve:
`loss_action(0) == loss_action(σ) − σ·ψ` AND `streaming(σ_a) == streaming(σ_b)`
for wildly different σ (the second is the decisive σ-leak test). The σ-free
primitive is then `streaming_action := loss_action(0)`, which single-sources the
walk (Pattern 2), names it (Pattern 3), and leaves `loss_action`'s FP tree —
hence the production composite matvec — BYTE-IDENTICAL.

Corollary for the leaf-sum twin (`InvertibleOperator.apply` override vs the
inherited `OperatorSum` leaf sum): both are value-EQUAL when the matvec is affine
in σ, so the teeth gate MUST be `array_equal` (0 ULP), not allclose — only
bit-identity discriminates leak-vs-override. Prove the teeth bite by DISABLING
the override (rename `apply`→`_DISABLED_apply`) and confirming all teeth fail at
the predicted ≤2 ULP, then restore byte-exact.

How to apply: never assume affinity from the literature — run the two probes.
A value-correct-by-coincidence twin is NOT a `catches(ERR)` (no wrong value
shipped) → `foundation` gate, `verifies(...)` only. (#257 S8b; #240 Step B.)

---

## L-012 -- A bare-ndarray return layout lies in its docstring — read the shape, not the prose

A solver-facing `external_source` returned `(N, ng, nx, ny)` (the solver-consumed
shape) while its docstrings claimed `(N, nx, 1, 1)`. A pre-existing cross-check
sliced the WRONG axes (`[:, :, 0, 0]` vs the correct `[:, 0, :, 0]`) and was
already red at clean HEAD — a latent bug masked as a passing-by-omission test.

How to apply: when consuming a bare-ndarray boundary, confirm the ACTUAL shape
(print `.shape` or read the producing code) — never trust the docstring's axis
order. This is the coding-elegance Pattern-13 hazard (shape-in-docstring); the
durable cure is a typed wrapper, but until then, verify the shape at the seam.
(#240/phase4 non-vacuum MMS.)

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

## L-014 -- A single-DOF scheme's scalar geometry factor is a 1x1 GRAM: find the inner product before widening

When widening a scheme from 1 DOF per cell to N (DD -> LD; scalar -> moment
vector), every shipped SCALAR geometry constant is the `(0,0)` corner of some
Gram matrix. Do NOT widen by "multiplying the scalar into each moment row" --
derive the weak form and read off WHICH inner product each constant belongs to.

The forensic instance (curvilinear LD shape derivation, 2026-08-25). The
curvilinear SN redistribution ships `redist_dAw = Delta A / w_n`, a scalar.
Weighting the conservative form by a per-cell basis against `dV = 4 pi r^2 dr`,
the angular term's `1/r` cancels ONE power of the measure, so the half-angle
fluxes enter through the ONE-MEASURE-DOWN Gram
`R_kj = Delta A * <b_k,b_j>_{d-2} / <b_0,b_0>_{d-2}` (`d = 3` sphere, `2`
cylinder -- this spelling also absorbs the geometry-dependent factor-of-2 in the
ORPHEUS alpha normalization). `R_00 == Delta A` exactly, so the shipped scalar
IS the corner; but `R_01/R_00 = h/(6 r_c) != 0` on the SPHERE, and the
off-diagonal is LOAD-BEARING: per-ordinate flat-flux consistency on the SLOPE
row cancels only because of it. Lumping `R` breaks the canonical L0 identity.

Two corollaries that generalise:

* **The scalar tells you nothing about the off-diagonal.** DD's arithmetic is
  compatible with any `R_01`, so no reading of the shipped code can recover it;
  only the weak form can.
* **Check whether the scheme's `p = 0` truncation IS the shipped scheme.** It
  usually is not: LD's 1-moment truncation is STEP (upwind-discontinuous), NOT
  DD (central/diamond). A base class that assumes the shipped member is the new
  member's degenerate case is wrong from the first line.

How to apply: before any moment-widening carve, write the weighted-residual
reduction and tabulate which MEASURE each term carries. The table IS the
convention crosswalk (L-006), and each row names the Gram the widened member
needs. (LD curvilinear shape probe; `scratch/probe_ld_sphere_moments.py` V1/V2/V4/V7.)

---

## L-015 -- Derive a gate's BLINDNESS from the algebra at design time, before writing it

A design choice whose two candidates AGREE at the fixture's fixed point is
invisible to that fixture -- at every tolerance, mesh and order (vv-principles
Mode 12). This is decidable by algebra BEFORE any gate exists: subtract the two
candidates symbolically and substitute the gate's own fixture.

Two instances from one derivation (2026-08-25), both in the curvilinear family's
canonical L0 gate (`streaming-equilibrium` / per-ordinate flat flux):

* **The starting-direction seed.** Two consistent LD discretizations of the
  `mu = -1` march. `[M]` at BOTH moment counts: the DD reduction differs by
  `propto (Sigma_t psi_in - q0)`, IDENTICALLY zero at the flat-flux fixed point;
  the full LD 2x2 has both variants flat-flux EXACT, hence agreeing there
  exactly. ⭐ And the second blindness stacks: the LD pole-cell gap is
  `(-h qhat/30, 4 h qhat/15) + O(h^2)` -- proportional to the SLOPE SOURCE
  MOMENT alone -- so any fixture with `qhat = 0` (which is every ORPHEUS fixture
  today, the external slope source being zeroed, #247) is blind too. Write the
  denominator: a claim measured on the 1-moment reduction is NOT a claim about
  the N-moment scheme, and here the two have DIFFERENT leading terms.
* **The geometry arm.** A per-moment-row flat-flux gate written on the CYLINDER
  is vacuous: one measure down from `r dr` is FLAT, so the redistribution Gram
  is diagonal AND both slope-row terms vanish identically (`0 = 0`). The same
  gate on the SPHERE is a genuine cancellation of two non-zero terms and DOES
  catch a lumped/mis-signed off-diagonal. Same gate, same code, one arm blind.

How to apply: for every design fork, compute `candidate_A - candidate_B`
symbolically and substitute the intended gate's fixture. If it vanishes, say so
in the plan and name the fixture that does NOT annihilate it -- do not ship the
fork under a gate that cannot rank it. Cheap (one `subs`), and it runs before
any test is written. (LD curvilinear shape probe V4/V6 + cylinder C3;
vv-principles Mode 12 / #24(d).)

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
