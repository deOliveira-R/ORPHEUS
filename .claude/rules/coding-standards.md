# Coding standards — the minimum-quality floor

These are **preemptive minimum standards** every contributor (main agent and every
sub-agent) follows by default, regardless of whether they are chasing elegance. They are
the *floor*. For the *ceiling* — recognizing how bad code manifests and refactoring toward
excellence — code-producing agents load the **`coding-elegance` skill**. (See also Cardinal
Rules 1 "Correctness" and 2 "Architecture" in CLAUDE.md.)

## Clean before extending

Before adding a capability to a class/module — especially on a new design — first run a
cleanup pass on that layer: collapse double paths, move concepts to their native place,
delete dead shims, fix twin sources of truth. The new capability should then land as a
**no-op extension through the single generic body**, not as a new arm grafted onto
structural debt.

- When a plan proposes a capability extension, insert a dedicated **cleanup phase before
  the capability phase**. Order findings into must-precede vs independent-polish vs
  explicit-WAIT; gate each cleanup substep as bit-identical where possible.
- Rationale: extending a layer that carries debt grows a third arm on every seam. (C5,
  2026-06-11: `from_axes` round-tripped axes→legacy-mesh→axes, so 3-D admission would have
  needed a new arm in the converter, constructor, AND trace gate — until the cleanup
  inversion let 3-D flow through the one generic body.)

## Type vs property — before minting a type

A representation earns its own **type** only when it is genuinely a different object, not the
same object wearing a label. The decidable (grep-checkable) criterion: mint a type **iff**
(a) there are **≥2 non-isomorphic bases/realizations** of the concept, AND (b) a
**non-identity morphism** is actually applied to it. Otherwise the concept is a **property**
— a field or flag on an existing type.

- If the only "change of basis" is the identity (one realization, no transform), a separate
  "type" is theatrics: it adds ceremony and a conversion seam without making any illegal
  state unrepresentable. Make it a property.
- This is the *type-minting* corollary of the `coding-elegance` "defer abstraction until ≥2
  instances" rule. Worked: an expanded-order spatial moment with a single basis and an
  identity change-of-basis is a `property`, not a `SpatialOrder` type.

## Retire as you go (aggressive retirement)

Superseded code is **noise that obscures signal** — it makes the codebase harder to read,
breeds "which path is canonical?" confusion, and invites accidental extension of the wrong
path. Retirement is a first-class deliverable, not optional cleanup.

- Every refactor that introduces a more elegant pattern **MUST retire its predecessor**.
- Deprecation aliases / compatibility shims live for **one merge cycle only** — remove them
  on the next pass. Never preserve backward-compat code unless the user explicitly
  authorizes it (e.g. a public API on a versioned release).
- Treat the **retirement audit as its own numbered substep**; a plan should carry a
  "retirement list" enumerating what gets deleted (with `file:line`).
- **Retirement means test migration:** retiring a symbol includes **rewiring its tests to
  the successor**, not deleting them with the symbol. Behavioral test (correctness contract)
  → rewire to the new API; API-smoke test (symbol exists) → delete; characterization test
  (e.g. FD-vs-WDD delta) → keep under `tests/<module>/characterization/`. Pure delete-only
  retirement is incomplete — it loses coverage. Inventory with `grep -rn "<symbol>" tests/`.
- **The retirement audit's blast radius is THREE searches, not one** (4–5 agents converged on
  this independently): (1) **graph callers** (`nexus impact`/`callers`) — necessary but NOT
  sufficient; the call graph misses property-reached leaves (`callers()==0` but live via a
  `cached_property`), class-name *bypass* consumers, and direct constructors of a guarded
  type; (2) **text-grep the symbol across code, tests, AND `docs/`** — an unresolved
  `:func:`/`:class:`/`:ref:` cross-reference renders as plain text with **no `-W` warning**,
  so the Sphinx build gate does NOT catch a retirement's doc blast radius; (3) **direct
  constructors** of any guarded type (a guard-at-source change reaches every `T(...)` caller,
  not just the factory path). Run all three, then retire.

### Exception — keep a relinquished *fuller view* as a verification oracle

Retirement targets a *superseded* predecessor (same job, done worse). It does NOT target a
**fuller view of a concept that an optimization relinquished** (full field → rolling window;
full angular flux → moments; dense operator → factored form). The fuller view is a
**verification pathway** that pins the optimized path's reference — keep it.

- Make the keep-vs-retire decision EXPLICIT; never let a fuller view fall out silently or
  sit half-alive (orphaned-but-undeleted). It is either a wired, exercised oracle or it is
  retired.
- The oracle is NOT production-reachable and MUST share the optimized path's kernel (only
  the representation/storage differs), pinned by a permanent **end-to-end** equivalence test
  (`optimized ≡ fuller-view-oracle`), bit-identical or principled-equivalent per
  `vv-principles`.
- **Discriminator:** keep it only if it is a genuine structural reference. A "fuller view"
  that is the SAME math merely *procedurally rearranged*, AND already verified by a
  structurally-independent oracle elsewhere (MMS, closed form), is genuine redundancy —
  retire it. (`vv-principles` L11: procedural ≠ structural independence.) Worked: the 2-D
  rolling-window sweep kept `_sweep_2d_full_field` as a typed-`WavefrontFlux` oracle; the
  per-ℓ scattering kernel was retired as mere procedural rearrangement.
