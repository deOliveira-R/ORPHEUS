# #236 Phase 2 τ-carve STEP B2 — fold the P2 c_in/c_out rebuild off DD.residual onto the CellVisit; complete the #226 Protocol→ABC matvec binding

**Branch** `feature/sn-spatial-angular-product` (B1 = `b7fed4d`; this work NOT committed/staged).
**Date** 2026-06-17. **BIT-IDENTICAL** (0-ULP, proven 3 ways).

## What B2 was

The SECOND of the four `c_in`/`c_out` inline-rebuild sites of the #236 Phase 2
crosswalk (`.claude/plans/issue_236_phase2_tau_carve_crosswalk.md` §7). The
weighted-diamond closure constants `c_out = α_out/τ`, `c_in = (1−τ)/τ·α_out + α_in`
are an ANGULAR-closure property; the closure `MorelMontryAngularSweep` is the
canonical owner (precomputes them per-μ-level at construction →
`_c_in_per_level`/`_c_out_per_level`; B1 minted the PUBLIC `(N,)` accessors
`c_in_per_ordinate`/`c_out_per_ordinate` on the Protocol + the
`PoleAngularClosureBase` ABC). B1 folded site P3 (the CollisionCache populator,
the free seam holding `sn_mesh`). B2 folds site P2 = `DiamondDifference.residual`
(the matvec twin @ `n_mask=1`).

## The design (settled — implemented)

`DiamondDifference` is deliberately STATELESS (reads only `CellVisit` +
`UpstreamState`, never the mesh/closure) → the closure-c cannot reach `DD.residual`
by coupling DD to the closure object (that breaks the spatial⊗angular separation).
So c travels as DATA on the `CellVisit`:

1. **`CellVisit` gained `c_in: float = 0.0` / `c_out: float = 0.0`** (frozen+slotted
   dataclass, `scheme.py:252-256`). Documented in the Attributes docstring as
   angular-closure-owned (#236 Phase 2 B2), distinct in provenance from the
   geometry-owned `streaming_terms`.

2. **NEW single production stamp `SNMesh._make_cell_visit`** (`geometry.py`, just
   before `_iter_cartesian_visits`). ALL FOUR `dag_walk` yield sites
   (`_iter_cartesian_visits`, `_iter_spherical_visits`, `_iter_cylindrical_visits`
   ×2 = degenerate + normal loops) funnel through it (Pattern 2 — the c-lookup lives
   in ONE place, NOT copied into 4 sites). It reads
   `self.pole_angular_closure.c_in_per_ordinate[global_ordinate]` /
   `.c_out_per_ordinate[global_ordinate]`. The GLOBAL-ordinate index is the load-
   bearing crosswalk row: **slab/sphere → `direction_idx`** (global ordinate IS
   ordinate_idx); **cylinder → `level_indices[mu_level_idx][m]` = `global_n`**
   (already computed at the top of `_iter_cylindrical_visits`) — mirrors EXACTLY
   what `ReducedStreamingOperator.streaming_terms` does to read the per-level α/τ.

3. **`DD.residual` reads `visit.c_in`/`visit.c_out`** (`diamond.py:307-309`) —
   replaced the inline `c_out = st.alpha_out/tau` / `c_in = (1-τ)/τ·α_out + α_in`.
   Dropped the now-dead `tau = st.tau_mm` local + the `st.alpha_*` reads. The
   `(ΔA/w)`-scaled assembly that follows is byte-UNCHANGED — only the SOURCE of c
   moved. `DD.update` (`diamond.py:230`, the 5th τ consumer / B3 scope) UNTOUCHED.

4. **#226 Protocol→ABC binding fix.** The matvec (`loss_representation.py`) reads
   `sn_mesh.pole_angular_closure`. Retyped the 4 consumer annotations from the
   narrow Protocol `PoleAngularClosure` → the ABC `PoleAngularClosureBase`
   (`geometry.py` ctor `__init__`, `_init_core`, the `self.pole_angular_closure`
   attribute, `from_axes` + its docstring), removed the now-unused Protocol import.
   **SURPRISE (the dispatch's mental model was incomplete):** retyping to the ABC
   cleared the `level_indices` errors (3) but NOT `precompute_psi_state` /
   `cell_contribution` / `angular_adjoint` — those were defined ONLY on the
   concrete subclasses, NOT declared on the ABC. The principled fix matching the
   dispatch's OWN cited precedent (`DiscretizationSchemeBase` declares
   `update`/`residual` as `@abstractmethod` so `mesh.scheme` consumers see the full
   contract): **declare the 3 strategy methods `@abstractmethod` on
   `PoleAngularClosureBase`** (`pole_angular_closure.py`, after the abstract
   `__call__`). Both concrete closures already implement all 3 (`__abstractmethods__`
   stays empty → still instantiable). Loose params (`psi_state: Any`,
   `precompute_psi_state` keyword args `np.ndarray` matching M-M the canonical
   implementer) to keep both overrides contravariant-compatible (first draft used
   `object`/`|None=None` → 3 NEW override-variance reds; fixed to `Any` + required-
   keyword → 0 new). `psi_half_seed` (line 3166, M-M-only attribute) LEFT — the
   dispatch HARD CONSTRAINT explicitly names it as `PsiHalfAngleSeed` #226 noise to
   route around.

## Defaulted-vs-required CellVisit-field DECISION + rationale

**DEFAULTED `c_in: float = 0.0` / `c_out: float = 0.0`.** Three reasons:
- **Slab/Cartesian correctness is FREE.** The `IdentityAngularClosure` carries
  neutral zeros (`α=0, τ=1 ⇒ c=0`) → `c_*_per_ordinate` returns `(N,)` zeros, so the
  default 0.0 IS the correct slab value. (Production slab still routes through
  `_make_cell_visit` for uniformity; the default just makes the slab path provably
  right even before the stamp.)
- **~25 direct `CellVisit(...)` test constructions still construct** without c.
- **BUT the round-trip contract forced a test-helper update.** The apply-vs-solve
  round-trip `test_residual_zero_at_solved_cell_avg` compares
  `residual(update(q).cell_avg, q) == 0`. `update` reads c via `cell_balance_terms`
  (rebuilds from `streaming_terms` — B3 scope, UNCHANGED in B2); `residual` now reads
  `visit.c_*`. For the round-trip to hold on CURVILINEAR, the fixture visit's c MUST
  equal what `cell_balance_terms` rebuilds. The fixtures build a `CellVisit` from a
  bare `ReducedStreamingOperator` with NO `SNMesh`/closure in scope → they cannot
  read the closure. SOLUTION: a tiny fixture-local surrogate
  `_c_from_streaming_terms(st)` (test_diamond.py) / `_visit_c(st)`
  (test_cell_balance_for_streaming.py) computing the SAME `α_out/τ` / `(1-τ)/τ·α_out+α_in`
  the closure would, stamped onto the curvilinear/slab fixture visits. This is the
  closure's job re-created as a fixture surrogate (it has the same `streaming_terms`
  data the closure derives c from) and pins the B2 round-trip contract honestly.
  Tests updated: `test_diamond.py` (`_slab/_sphere/_cylinder/_cylinder_degenerate
  _visit_inputs` + the inward/outward sphere residual test) ;
  `test_cell_balance_for_streaming.py` (the residual-consumes-helper + round-trip
  tests). `test_ordinate_scan.py` constructs visits but feeds ONLY `update` (which
  reads from `streaming_terms`, never `visit.c_*`) → default 0.0 fine, NO change.

If the fields had been REQUIRED, every test fixture would have to compute c inline
(the same surrogate, but at ~25+ sites) AND the production stamp would still be the
sole real source — required gives nothing over defaulted-0.0 here, and breaks all
existing direct constructions. Defaulted-0.0 is the Pattern-4-respecting choice: the
neutral element is the safe default, the production stamp + the round-trip fixtures
carry the real value.

## BIT-IDENTITY — proven 3 ways (0-ULP)

1. **In-process `abs(visit.c − inline_c(st))` max-diff EXACTLY 0.000e+00** for BOTH
   c_in AND c_out across sphere (single-level) `(nx,quad)∈{(8,GL8),(12,GL6)}` +
   cylinder (multi-level) `{(8,prod2×4),(12,prod2×8)}`, over EVERY visit of EVERY
   ordinate (64/72/64/192 visits). The closure-sourced `visit.c` (via
   `_make_cell_visit` → `c_*_per_ordinate[global_ordinate]`) equals the old inline
   `(1-τ)/τ·α_out+α_in` / `α_out/τ` byte-for-byte. (The closure precomputes the whole
   per-level array with the SAME elementwise op order; the inline did the single
   element; gather is a pure permutation.) Proves the global-ordinate mapping is
   byte-correct for both single-level and multi-level.
2. **Matvec twin** (DD.residual feeds it): `test_unified_matvec_{sphere,cylinder}.py`
   = **6 passed, 27 xfailed** (the 27 = pre-existing #206/ERR-058; NO new failure,
   NO xpass).
3. **DriftWarning-escalating snapshots**: `tests/sn/sweep/core tests/sn/solve` =
   **520 passed, 1 skipped, 4 xfailed, NO DriftWarning fired**. The known
   `cyl_1g_homogeneous_product_dd_n20` 3.9e-11 did NOT grow.

## VERIFICATION (all PASS, L12 paste-back in the closeout message)

- Gate 1 `tests/sn/sweep/core` → **460P/1skip/4xf** (diamond/cell_balance/conformance,
  incl. the curvilinear round-trip + residual-consumes-helper).
- Gate 2 `test_unified_matvec_{sphere,cylinder}` → **6P/27xf** (no new fail/xpass).
- Gate 3 `sweep/core + solve` (DriftWarning escalate) → **520P/1skip/4xf**, no drift.
- Gate 4 anchors (g-adjoint-reciprocity + curvilinear-MMS + tau-producer-equiv) →
  **23P** (g-adjoint exercises `angular_adjoint` through the ABC).
- Gate 5 `sweep/curvilinear -k "not unified_matvec"` → **168P/33desel** (DD.residual's
  other consumer, the 1-D curvilinear sweep).
- Gate 6 pyright: consumers `geometry.py + loss_representation.py` **60→56** total;
  `precompute_psi_state`/`cell_contribution`/`level_indices` "unknown attr" errors ALL
  GONE; only remaining PoleAngularClosure* attr error = `psi_half_seed:3166` (the
  named #226 noise). Closure module `pole_angular_closure.py` 33→30 (my abstract block
  adds 0 errors).
- Sphinx build exit 0; new label `sn-closure-c-on-cellvisit` + `:mod:`/`:meth:`/`:attr:`
  xrefs + archivist TODO render clean in HTML + Nexus.

## Changed files (B2 only)

- `orpheus/sn/spatial/scheme.py` — CellVisit `c_in`/`c_out` fields + docstring.
- `orpheus/sn/geometry.py` — `_make_cell_visit` helper; 4-site rewire; binding-fix
  retypes (Protocol→ABC) + import swap (`StreamingTerms` in, Protocol out).
- `orpheus/sn/spatial/diamond.py` — DD.residual c-rewire (inline → `visit.c_*`).
- `orpheus/sn/spatial/pole_angular_closure.py` — 3 `@abstractmethod` strategy
  declarations on the ABC (binding-fix completion).
- `tests/sn/sweep/core/test_diamond.py` — `_c_from_streaming_terms` surrogate + 4
  fixture-builder + inward/outward-sphere visit stamps.
- `tests/sn/sweep/core/test_cell_balance_for_streaming.py` — `_visit_c` surrogate +
  2 residual-test visit stamps.
- `docs/theory/discrete_ordinates.rst` — `sn-closure-c-on-cellvisit` B2 stub.

(PRE-EXISTING uncommitted, NOT mine: `.claude/skills/vv-principles/*` +
`.claude/agent-memory/elegance-enforcer/*` = prior #240 D5b ERR-060/061/062 work;
LEFT untouched per L28.)

## OWED (NOT B2 — explicitly out of scope)

- **B3** — rewire P1 (`cell_balance_terms`, the DD.update solve path) to read the
  visit's c; re-source DD.update:230's angular-recurrence τ from the closure (the 5th
  consumer). After B3, the 4 c-sites are all folded.
- **C** — retire the geometry-side τ producers + `StreamingTerms.tau_mm`/`alpha_in`/
  `alpha_out` (orphaned once the c-folds land) + migrate the τ-bit-identical geometry
  tests onto the closure producer.
- archivist DISPATCH emitted (followup:false) for the rich narrative of
  `sn-closure-c-on-cellvisit`.

## LESSON (for algebra-of-record / coding-elegance)

A "retype the consumer to the ABC" binding fix is only complete if the ABC actually
CARRIES the full strategy contract. Here the ABC declared only `__call__`; the matvec-
driven methods (`precompute_psi_state`/`cell_contribution`/`angular_adjoint`) lived on
the concrete subclasses, so the retype cleared only the ABC-resident `level_indices`
and the method errors merely changed their class name (Protocol→ABC). The principled
completion — declaring those methods `@abstractmethod` on the ABC — IS the documented
precedent (`DiscretizationSchemeBase` does exactly this for `update`/`residual`), so
the binding fix is "complete the ABC contract", not "retype and hope". WATCH the
override-contravariance trap when minting the abstract signatures: a base param that
is a SUBTYPE of a concrete override (e.g. base `object`, override `tuple[...]`) is an
incompatible override; use `Any` (bidirectionally compatible) or match the canonical
implementer's required-keyword signature. And the c-fold's airtight bit-id proof is
the in-process `np.array_equal`/`abs(diff)==0` on the stamped `visit.c` vs the inline
rebuild at the EXACT sphere(single-level)+cylinder(multi-level) configs — isolating
the closure-source fold from the pre-existing matvec xfails a whole-suite run conflates.
