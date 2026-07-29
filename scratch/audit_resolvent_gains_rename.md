# L20 three-surface dependency audit — `WithinGroupSystem.resolvent` → `implicit_operator`, `.gains` → `explicit_gains`

**Audited at** `refactor/operator-naming-honesty` @ `8654d348`.

> ⚠ **L-007 — the tree MOVED under this audit.** `orpheus/sn/solver.py` carries **uncommitted
> edits** (46 lines, a *different* item of the same naming-honesty campaign: retiring the
> production-dead `self.L`/`self.S`/`self.F` triple and the `rebind_cross_sections` rebuild).
> Line numbers in `solver.py` shifted **+3 then +6** between my first grep and my last. Every
> `solver.py` number below is from the **final re-grep** (the authoritative pass) — but the
> file is live, so **re-grep before editing**. The in-flight diff does **not** touch
> `resolvent`/`gains`; the two changes are independent.

**Target definition** (`orpheus/sn/coupled_system.py:282-328` — this file is NOT under edit):

```python
@dataclass(frozen=True)
class WithinGroupSystem:
    loss: "CoupledOperator"
    space: "CoupledSpace"
    resolvent: "CoupledOperator | InvertibleOperator"     # ← rename target #1  (holds M)
    gains: "tuple[LinearOperator, ...]"                   # ← rename target #2  (holds N)
```

---

## 0. Headline verdicts (read this first)

1. **The premise is TRUE.** `resolvent` holds `M` — the *un-inverted forward* operator of the
   regular splitting `A = M − N`. The resolvent proper is `M⁻¹`, produced downstream by
   `.inverse()` / `.solve`. Read from the constructor: `resolvent=LC` (seedless,
   `coupled_system.py:486`) and `resolvent = CoupledOperator([[LC, A_AB], [None, march]], …)`
   (carrying, `:531`). The tree's own prose already says so: `solver.py:821` calls it
   *"the un-inverted FORWARD ``M``"*.

2. **The "74 test files" figure is a FALSE ALARM — the real number is 12.**
   `grep -rl resolvent tests/` = **72 files**, but only **12** reference the FIELD
   (`system.resolvent` or the mirror-named `_select_si_resolvent`). The other **60** are the
   `trajectory_resolvent` reference-solver corpus and mathematical prose. **The blast radius
   is 5× smaller than the plan review estimated.**

3. **Full (T) change set: ≈75 lines across 20 files** — 16 in `coupled_system.py` (definition
   + builder + its own docstrings), 20 in `solver.py` (+ 20 more if the mirror-named locals
   `_select_si_resolvent`/`base_resolvent`/`final_resolvent` are renamed too), 30 in 12 test
   files, 9 in 6 `.rst` pages, plus ~30 `.claude/` archaeology lines that should NOT be
   rewritten. **The per-file tables in §1/§C are the authoritative list** — the counts here
   are a size estimate, not the change set.

4. **`resolvent` rename: mechanical. `gains` rename: mechanical BUT with two traps.**
   - **Trap 1 (`gains` substring).** `"gains"` is a substring of **`"against"`**
     (a-**gains**-t). There are **679 occurrences of "against"** in `orpheus/` + `tests/`
     `.py` alone. A naive `sed s/gains/explicit_gains/g` produces `aexplicit_gainst`.
     **Never do a bare-substring replace on `gains`.** Use `\bgains\b` with a
     leading-`.`/`=`/`*` anchor, or edit the ~30 sites by hand.
   - **Trap 2 (positional + arity contract).** See §B: the seedless `(S, B_a)` tuple is
     destructured positionally at 3 production sites and ≥8 test sites, with `B_a` LAST and
     one site using `gains[-1]`. A field rename leaves the contract intact, but every one of
     those lines is touched, so a slip there is silent-wrong, not loud-wrong.

5. **`orpheus/diffusion/solver.py` has an unrelated, HONEST `self.resolvent`** — it holds
   `MatrixInverseOperator(...)`, i.e. genuinely `A⁻¹`. Class (O): **do not rename it.**
   It is, however, the strongest evidence the SN field's name is a real collision: two
   solver-side attributes spelled `resolvent`, one holding `M`, one holding `M⁻¹`.

6. **Five prose sites mix (T) with (C)/(O) meanings inside one docstring or one
   paragraph** — the L25 `replace_all` hazard. Enumerated in §1.3. **`replace_all` on
   `resolvent` is UNSAFE in `coupled_system.py`, `solver.py`, `boundary_conditions.rst`, and
   `operator_adjoint.rst`.**

7. **Name-collision check is CLEAN.** `implicit_operator`, `explicit_gains`, `M_operator`:
   **ZERO** occurrences in code, tests, or docs (the only hits are the plan documents that
   propose them). No identifier named `implicit` exists. §D.

8. **The corpus already carries a DOCUMENTED double-binding of "resolvent"** — and the
   rename resolves it. `docs/theory/conventions/notation.rst` row 8 explicitly names this as
   *"the one place ORPHEUS itself carries two bindings of A, by design"*. §D.2. **That row
   must be updated as part of the rename**, or the rename creates a new drift instead of
   closing one.

---

## 1. Class (T) — THE TARGET. The change set.

### 1.1 `orpheus/sn/coupled_system.py` — definition + builder (NOT under concurrent edit)

| file:line | literal text | class | reasoning |
|---|---|:--:|---|
| `coupled_system.py:214` | ``record's ``resolvent``/``gains``.`` | **T** | `build_coupled_system` docstring naming both fields of the returned record. |
| `:290` | ``grid (the equation), ``resolvent``/``gains`` its regular splitting`` | **T** | `WithinGroupSystem` class docstring, names both fields. |
| `:309` | `resolvent : CoupledOperator | InvertibleOperator` | **T** | numpydoc Parameters entry — **must** match the field name or the docstring silently lies. |
| `:317` | `gains : tuple[LinearOperator, ...]` | **T** | numpydoc Parameters entry. |
| **`:327`** | **`resolvent: "CoupledOperator | InvertibleOperator"`** | **T** | **THE FIELD.** Rename anchor #1. |
| **`:328`** | **`gains: "tuple[LinearOperator, ...]"`** | **T** | **THE FIELD.** Rename anchor #2. |
| `:486` | `resolvent=LC,` | **T** | constructor kwarg, seedless arm. |
| `:487` | `gains=(S, B_a),` | **T** | constructor kwarg, seedless arm. **This line IS the `B_a`-LAST convention** (§B). |
| `:491` | `# grid, the gain grid, AND the resolvent grid (single-sourced objects,` | **T** | comment naming the field/grid triple. |
| `:517` | `# The GAIN grid N = M − A: all POSITIVE (rhs gains); the (A,B) slot is` | **T** | comment; "gains" here is the field's concept. |
| `:531` | `resolvent = CoupledOperator(` | **T** | **local variable** feeding the kwarg — rename for coherence (it is the same object). |
| `:537` | `resolvent=resolvent,` | **T** | constructor kwarg + local, carrying arm. |
| `:538` | `gains=(N,),` | **T** | constructor kwarg, carrying arm (arity **1**). |
| `:79` | ``sweepable part inverted every step, ``N`` the lagged coupling gains.`` | **T-prose** | module docstring; this sentence IS the source of both field names. Rewrite to name the new fields. |
| `:102` | ``all POSITIVE here (gains on the rhs: ``rhs = q + N·ψ``)`` | **T-prose** | module docstring, gain-grid sign narration. |
| `:428` | ``(``[[S+B_a, ∅], [Emission, B_b]]`` — gains on the rhs)`` | **T-prose** | builder docstring sign table. |

**(C)/(O) hits inside the SAME file — do NOT rename:**

| file:line | literal text | class | reasoning |
|---|---|:--:|---|
| `coupled_system.py:84` | ``orpheus.sn.solver._within_group_triple`` / ``_lagged_gains`` pair retired`` | **O** | a RETIRED symbol's tombstone. Historical record; renaming it falsifies history. |
| `:295` | ``the ``_within_group_triple``/``_lagged_gains`` construction twin.`` | **O** | same tombstone. |
| **`:401`** | ``resolvent whose ``solve`` is the WDD sweep.`` | **C** ⚠ | **THE COLLISION.** In `build_within_group_system`'s own docstring, describing `L+C` — here "resolvent" is the CORPUS meaning ("the thing whose `.solve` inverts"). It sits **62 lines above** `resolvent=LC` in the same function. A `replace_all` corrupts it. See §1.3. |
| `:180`, `:306` | `…emitted together`, `The coupled carrier space ``loss`` is typed **against**…` | **O** | **substring false positives**: `"against"` contains `"gains"`. Trap 1. |

### 1.2 `orpheus/sn/solver.py` — the consumers (⚠ UNDER CONCURRENT EDIT; re-grep before editing)

| file:line | literal text | class | reasoning |
|---|---|:--:|---|
| `solver.py:838` | `if isinstance(system.resolvent, CoupledOperator):` | **T** | the carrying/seedless structural dispatch — reads the field as the presence discriminator. |
| `:851` | `system.resolvent.inverse(), *system.gains,` | **T ×2** | coupled SI arm. **Unpacking site** (§A). |
| `:854` | `return si, system.resolvent, system.gains, False` | **T ×2** | coupled SI return. |
| `:856` | `S, B = system.gains` | **T** | **arity-2 positional destructure** of the seedless tuple, immediately followed by an `isinstance(S, ScatteringOperator)` drift guard (`:857-862`). §B. |
| `:864` | `system.resolvent, S, B, sn_mesh, inner_schedule,` | **T** | feeds `_select_si_resolvent`. |
| `:1584` | `coupled = isinstance(system.resolvent, CoupledOperator)` | **T** | eigenvalue SI path presence probe. |
| `:1765` | `coupled = isinstance(system.resolvent, CoupledOperator)` | **T** | eigenvalue Krylov path presence probe. |
| `:1795` | `system.resolvent, *system.gains,` | **T ×2** | eigenvalue Krylov matvec operands. **Unpacking site** (§A). |
| `:2146-2148` | `final_resolvent = build_within_group_system(…).resolvent` | **T** (on `.resolvent`) | the post-convergence full-angular reconstruction. The **local** `final_resolvent` (`:2146, 2156, 2202, 2215, 2222`) mirrors the field name — rename for coherence (5 lines). |
| `:2357` | `gain = system.gains[0]` | **T** | adjoint posing; **arity-generic fold** start. §B. |
| `:2358` | `for extra in system.gains[1:]:` | **T** | the fold tail — the ONE production site that is arity-agnostic. |
| `:2369` | `return system.resolvent, gain, F, full_field_zero` | **T** | adjoint posing return, seedless. |
| `:2409` | `return system.resolvent, gain, F_posed, space.zeros()` | **T** | adjoint posing return, carrying. |
| `:3185` | `coupled = isinstance(system.resolvent, CoupledOperator)` | **T** | fixed-source SI presence probe. |
| `:3411` | `coupled = isinstance(system.resolvent, CoupledOperator)` | **T** | fixed-source Krylov presence probe. |
| `:3417` | `system.resolvent, *system.gains,` | **T ×2** | fixed-source Krylov operands. **Unpacking site** (§A). |
| `:712` | ``r"""Pick the ``(resolvent, gains)`` for the SEEDLESS within-group SI`` | **T-prose** | `_select_si_resolvent`'s docstring — names the record's fields it mirrors. |
| `:798` | ``* **coupled** (the record's ``resolvent`` is the triangular`` | **T-prose** | `_within_group_si` docstring, explicit "the record's". |
| `:474` | ``block matvec) and ``*gains`` is ``N`` (the ``(S, B_a)`` pair seedless;`` | **T-prose** | `_within_group_krylov` docstring documenting the fed-in field. |
| `:1549` | ``# since B.2d the system record's ``gains``): each SI iterate`` | **T-prose** | explicit "the system record's". |

**Mirror-named locals / helpers in `solver.py` — (O) by strict classification, but they
exist ONLY because the field is named `resolvent`. Renaming the field and leaving these
creates a two-spelling seam. Recommend renaming in the same change; flagging separately so
the user can decide:**

| file:line | literal text | class | note |
|---|---|:--:|---|
| `:707, 811, 863, 954` | `_select_si_resolvent` (def + 3 references) | **O-mirror** | module-private helper; its docstring `:712` already spells `(resolvent, gains)`. 20 occurrences tree-wide (incl. 6 test files + 3 `.rst`). |
| `:818, 821, 863, 866, 876, 3174, 3180, 3267, 3276` | `base_resolvent` | **O-mirror** | local + return-tuple name for the same `M`. 10 occurrences. |
| `:2146, 2156, 2202, 2215, 2222` | `final_resolvent` | **O-mirror** | local for the same `M`. 5 occurrences. |
| `:459, 514` | `*gains` (parameter of `_within_group_krylov`) | **O-mirror** | the driver's own variadic parameter, fed by `*system.gains`. |

**Pure (O) in `solver.py` — leave alone:**
`:561` ("the windowed resolvent emits"), `:1949/1955` ("the G-S resolvent absorbs"), `:1266`,
`:1502`, `:1950`, `:2143`, `:2223`, `:2326/2329` (`_adjoint_posing_parts` docstring return
name), `:2994`, `:3127`, `:3145`, `:3175`, `:3251/3257`, `:3263`, `:3274` — mathematical prose
about "the resolvent" / "coupling gains" as concepts, not the field.

### 1.3 ⚠ SAME-SENTENCE / SAME-DOCSTRING COLLISIONS — prose rewrite required, NOT a token swap

These are the L25 sites. In each, the (T) meaning and a (C)/(O) meaning appear close enough
that a mechanical swap produces a sentence that is *wrong*, not merely inconsistent.

| # | file:line | the collision |
|---|---|---|
| **1** | `orpheus/sn/coupled_system.py:401` vs `:486`/`:531`/`:537` | Inside **one function** (`build_within_group_system`), `:401` says *"``L + C`` fuses to the InvertibleOperator — the **resolvent** whose ``solve`` is the WDD sweep"* (corpus meaning: the invert-er), while `:486`/`:537` assign `resolvent=` (the field, holding the forward `M`). **Both meanings, one docstring.** Rewrite `:401` to *"the invertible composite whose ``solve`` is the WDD sweep"*. |
| **2** | `orpheus/sn/coupled_system.py:290-291` | *"``resolvent``/``gains`` its regular splitting (Hackbusch 2016 §11 — the drivers iterate ``ψ ← M⁻¹(q + N·ψ)``)"* — the field names and the **inverse** `M⁻¹` in adjacent clauses. The rename is exactly what makes this sentence readable; write it deliberately. |
| **3** | `orpheus/sn/solver.py:812` + `:821` | ``the INVERSE build (``base_resolvent.inverse()``…)`` … ``* ``base_resolvent`` — the un-inverted FORWARD ``M``.`` The docstring is already *apologising* for the name. Rewrite, don't swap. |
| **4** | `docs/theory/foundations/boundary_conditions.rst:683-701` | Four uses in 19 lines: `:683` the LaTeX macro `A_{\rm resolvent}` (driver-operand notation), `:684` "one invertible resolvent A" (corpus), `:699` "the resolvent **already inverted**", `:701` "keeps the *forward* resolvent for its GMRES matvec". Then `:733-736` names the FIELD. Same section, three meanings. |
| **5** | `docs/theory/foundations/operator_adjoint.rst:504` | ``resolvent`` ``(L+C)`` and lagged ``gains`` ``(S, B_a)`` keep ``S`` /` — both field names in one clause, inside an adjoint-metric discussion that also uses "resolvent" mathematically. |

---

## 2. Class (C) — the CORPUS meaning. **Must be left alone.**

The corpus uses "resolvent" for the genuine `A⁻¹`-family object. These are correct
mathematics and renaming any of them would be a regression.

| locus | count | class | reasoning |
|---|--:|:--:|---|
| `docs/theory/foundations/operator_algebra.rst` | 44 | **C** | The canonical definition. `:3960` (Key Facts, eigenvalue posing): *"Its power-method realization is the dominant eigenpair of the **resolvent** `A_loss⁻¹ M`"*; `:239`, `:3971`, `:4012`, `:4023` (equation label `eigen-resolvent`). **This is the meaning the rename defers to.** |
| `docs/theory/conventions/notation.rst:20, 299` | 2 | **C** | The symbol-table row 8 binding: *"`A` the **invertible resolvent operand**"*. See §D.2 — this row is the documented bridge and **must be revised**, not deleted. |
| `docs/theory/methods/sn/slab_multigroup.rst` | 20 | **C** | Layer-3 resolvent / resolvent-eigenproblem narrative (`:609`, `:634`, `:794` `A_loss⁻¹M`, `:845-875` the triple-vs-monolithic split). |
| `docs/api/numerics.rst:61, 70, 122` | 3 | **C** | `A⁻¹F` resolvent prose for `direct_eigenvalue`/`dominant_eigenpair`. |
| `orpheus/numerics/eigenvalue.py` | 15 | **C** | `dominant_eigenpair` / `direct_eigenvalue` docstrings — the materialized `M = A⁻¹F`. Genuinely the inverse. |
| `docs/theory/foundations/operator_inverse_family.rst` | 13 | **C** | The **duck-typed-resolvent retirement** narrative (`:86-94`, `:181-191`) — historical prose about what the word used to mean. Historical record; do not rewrite. |
| `docs/theory/methods/sn/history.rst` (most) | ~9 of 11 | **C/O** | Development-history changelog. `:75` and `:219-220` are (T) (§C); the rest is history. **Do not rewrite history entries** — see §C note. |

---

## 3. Class (O) — other / unrelated. **Must be left alone.**

| locus | count | class | reasoning |
|---|--:|:--:|---|
| `orpheus/derivations/continuous/trajectory_resolvent/**` | ~150 (package) | **O** | The Peierls Variant-α Green's-function reference-solver **package name**. Completely unrelated concept (a trajectory/chord resolvent). |
| `docs/theory/references/trajectory_resolvent.rst` | **219** | **O** | The theory page for that package — the single largest concentration of the word in the tree, and 0 % of it is the target. |
| `tests/derivations/test_*peierls*`, `test_trajectory_resolvent_*` | ~60 files | **O** | The reference-solver test corpus. **This is the bulk of the "74 test files" scare.** |
| `tests/cross_method/{adapters,cases,protocol,test_*}.py` | ~95 | **O** | `method: str = "trajectory_resolvent"` registry keys + adapter class names. A rename here would break the cross-method registry, which is string-keyed. |
| `orpheus/derivations/continuous/spectral_resolvent` + `docs/theory/references/spectral_resolvent.rst` | 10 | **O** | The *β* meaning (spectral resolvent) — a third distinct concept, documented at `docs/theory/verification/reference_solutions.rst:124-135` as **meanings (α) and (β)**. |
| `orpheus/diffusion/solver.py:43, 206, 250, 293` | 4 | **O** | `DiffusionSolver.resolvent` = `MatrixInverseOperator(...)`. **Genuinely `A⁻¹`.** Honest. Leave. |
| `orpheus/sn/operators/radial_characteristic.py:64, 190, 295, 453, 477, 518` | 6 | **O** | "System B's **resolvent** `A_BB⁻¹` as `solve`" — correct corpus usage on a different operator. |
| `orpheus/numerics/iteration.py`, `numerics/operator.py`, `numerics/flat_operator.py`, `numerics/coupled_system.py:963` | ~24 | **O** | "the INVERTIBLE resolvent operand" module-head prose (the notation.rst row-8 binding). See §D.2 — coherence question, not a rename target. |
| `docs/theory/verification/matrix.rst` | 12 | **O** | Auto-generated V&V table rows: equation labels `eigen-resolvent`, `apply-solve-cell-resolvent`, `resolvent-object-gate`, `billiard-resolvent-neumann`, `resolvent-similarity`, `peierls-greens-*-resolvent`. **These are `vv-status` LABEL tokens — renaming any breaks a `:ref:` and the V&V matrix.** |
| `orpheus/geometry/{mesh,structured_geometry}.py`, `docs/.../structured_geometry.rst` | 3 | **O** | "the trajectory_resolvent `Billiard` class" cross-references. |
| `orpheus/derivations/continuous/{singular_eigenfunction,galerkin_spectral,peierls_nystrom,fn_method}/**` | ~30 | **O** | Cross-references to the trajectory/spectral resolvent references. |

---

## A. The `*system.gains` unpacking / indexing sites

`gains: tuple[LinearOperator, ...]` is consumed in **four distinct shapes**. Enumerated so a
future *shape* change (as opposed to this rename) knows its blast radius.

### A.1 Variadic splat `*system.gains` — 3 production + 3 test sites

| file:line | site | shape dependence |
|---|---|---|
| `orpheus/sn/solver.py:851` | `SourceIteration(system.resolvent.inverse(), *system.gains, …)` | **arity-agnostic** — `SourceIteration` is variadic (`iteration.py:592 *gains`). |
| `orpheus/sn/solver.py:1795` | `_within_group_krylov(system.resolvent, *system.gains, …)` (eigenvalue) | **arity-agnostic** — `_within_group_krylov(LC, *gains, …)` (`solver.py:459`). |
| `orpheus/sn/solver.py:3417` | `_within_group_krylov(system.resolvent, *system.gains, …)` (fixed-source) | **arity-agnostic**. |
| `tests/sn/operators/test_psi_half_coupling.py:3590` | `SourceIteration(system.resolvent.inverse(), *system.gains, …)` | arity-agnostic. |
| `tests/sn/operators/test_psi_half_coupling.py:3657` | `_within_group_krylov(system.resolvent, *system.gains, …)` | arity-agnostic. |
| `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:183` | `_within_group_krylov(system.resolvent, *system.gains, …)` | arity-agnostic. |

### A.2 Positional destructure — **arity-2 AND order dependent**

| file:line | site |
|---|---|
| `orpheus/sn/solver.py:856` | `S, B = system.gains` (+ the `isinstance(S, ScatteringOperator)` guard at `:857-862`) |
| `tests/sn/operators/test_apply_full_field_codomain.py:283` | `LC, (S, B) = system.resolvent, system.gains  # seedless slab record shape` |
| `tests/sn/operators/test_typed_residual_evaluation.py:67, 211` | same idiom, twice |
| `tests/sn/solve/test_2d_anisotropic_windowing.py:283` | same idiom (`# seedless 2-D record shape`) |
| `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py:197` | same idiom (`# seedless record shape`) |
| `tests/sn/operators/test_psi_half_coupling.py:3269` | `S_g, B_g = s_system.gains` (after an explicit `len(...) != 2` check at `:3267`) |

### A.3 Index / slice access

| file:line | site | dependence |
|---|---|---|
| `orpheus/sn/solver.py:2357-2358` | `gain = system.gains[0]` / `for extra in system.gains[1:]` | **arity-GENERIC fold** — the one production site that works for any arity. Reduces `N` to a single summed operator for the adjoint posing. |
| `tests/numerics/test_iteration.py:788` | `S_total = system.gains[0] + system.gains[1]` | **arity-2 hard-coded**. |
| `tests/sn/operators/test_psi_half_coupling.py:308, 341, 561, 2149, 2195, 2377, 2872, 3239, 3479` | `system.gains[0]` (the carrying-arm `N` grid) | **arity-1 (carrying)**. `:2377` chains `.gains[0].blocks[1][0]`. |
| **`tests/sn/operators/test_psi_half_coupling.py:593`** | **`B = slab_system.gains[-1]`** | **NEGATIVE index — depends on `B_a` being LAST.** The literal encoding of the documented convention. |
| `tests/sn/operators/test_psi_half_coupling.py:2163` | `if any(isinstance(g, CoupledOperator) for g in slab_system.gains)` | iterates — arity-agnostic. |

### A.4 `dataclasses.replace` kwarg — 1 site

| file:line | site |
|---|---|
| `tests/sn/operators/test_psi_half_coupling.py:2199` | `return replace(system, gains=(crippled,))` — a monkeypatched `build_within_group_system` returning a crippled `N` grid (the emission-unwired negative control). **This is a `gains=` KWARG use in a test** — the kind of site a `system.gains`-only grep misses. Found via the `gains=` pattern. |

**Verdict for §A:** every site is a *name* dependence only. No site depends on the field
being called `gains`; all 30 rewrite mechanically to `explicit_gains`. **No shape change is
implied by the rename.**

---

## B. The presence-dependent arity and the positional convention

### B.1 The two shapes the builder emits

| mesh | `resolvent` (M) | `gains` (N) | arity |
|---|---|---|--:|
| **carrying** (sphere; `radial_characteristic_field_space is not None`) | `CoupledOperator([[LC, A_AB], [None, march]])` (upper-triangular 2×2) | `(N,)` — **ONE** `CoupledOperator` gain grid `[[S+B_a, ∅], [Emission, B_b]]` | **1** |
| **seedless** (slab / cylinder / multi-D Cartesian) | `LC` — a plain `InvertibleOperator` | `(S, B_a)` — the bare pair, **`B_a` LAST** | **2** |

Source: `coupled_system.py:483-488` (seedless) and `:534-539` (carrying). The documented
convention is in the field's own docstring, `coupled_system.py:320-322`:

> ``(S, B_a)`` tuple seedless (``B_a`` LAST — the boundary-gain convention the G-S schedule
> arm parses).

### B.2 The code that PARSES it positionally

**The claim in the brief is correct, and here it is.** The parse is in
`orpheus/sn/solver.py::_within_group_si`, immediately before the G-S selector call:

```python
# solver.py:855-865  (⚠ line numbers live — file under concurrent edit)
    # Seedless: the record's gains are the (S, B_a) pair — loud on drift.
    S, B = system.gains
    if not isinstance(S, ScatteringOperator):
        raise TypeError(
            f"_within_group_si: the seedless record's first gain must be "
            f"the ScatteringOperator (the builder's (S, B_a) convention); "
            f"got {type(S).__name__}."
        )
    base_resolvent, gains = _select_si_resolvent(
        system.resolvent, S, B, sn_mesh, inner_schedule,
    )
```

Three separate contracts are live in those ten lines:

1. **Arity = 2** — `S, B = system.gains` raises `ValueError` on any other length. The
   carrying arm never reaches here (guarded by `isinstance(system.resolvent, CoupledOperator)`
   at `:838`, which returns early at `:854`).
2. **Position 0 = the scattering gain** — enforced *loudly* by the `isinstance` guard
   (`:857-862`). Drift here is a `TypeError`, not a wrong number.
3. **Position 1 (= LAST) = the boundary gain `B_a`** — **NOT** guarded at this site. It is
   guarded one level down, inside `_select_si_resolvent` (`solver.py:766-772`):

```python
# solver.py:766-772
        if not isinstance(B, SNBoundaryOperator):
            raise TypeError(
                "boundary Gauss-Seidel split requires the plain "
                "SNBoundaryOperator (a seedless multi-D Cartesian mesh); got "
                f"{type(B).__name__} — a seed-carrying composite must not reach "
                "the G-S schedule path (RULING P1: gradings live on B_a)."
            )
        parts = B.split(SweepSchedule.gauss_seidel(sn_mesh))
        return LC - parts.lower, (S, parts.upper)
    return LC, (S, B)
```

Note the **`_select_si_resolvent` return re-emits the same `(S, boundary)` positional
convention** (`:774` and `:775`) — so the convention propagates one hop past the record.

**Important scoping fact:** the `B_a`-LAST guard fires **only on the Gauss-Seidel arm**
(multi-D Cartesian). On the Jacobi arm (`:775 return LC, (S, B)`) position 1 is passed
through untyped. So the convention is *loudly* checked at position 0 always, and at
position 1 only under G-S.

### B.3 Which consumers depend on ORDERING vs ARITY

| consumer | depends on ORDER | depends on ARITY | note |
|---|:--:|:--:|---|
| `solver.py:856` `S, B = system.gains` | ✅ (pos 0 = S, guarded) | ✅ (=2) | production, seedless SI |
| `solver.py:766` `isinstance(B, SNBoundaryOperator)` | ✅ (pos 1 = B_a) | — | production, G-S arm only |
| `solver.py:774/775` `return …, (S, parts.upper)` / `(S, B)` | ✅ (re-emits) | ✅ (=2) | production, propagates the convention |
| `solver.py:2357-2358` `gains[0]` + `gains[1:]` fold | ❌ (sum is commutative) | ❌ **generic** | production, adjoint posing — **the arity-safe site** |
| `solver.py:851 / 1795 / 3417` `*system.gains` | ❌ | ❌ | production, variadic splat |
| `test_psi_half_coupling.py:593` `gains[-1]` | ✅ (**B_a LAST**) | ❌ | test, pins the convention |
| `test_psi_half_coupling.py:3267-3269` | ✅ | ✅ (`len != 2` fail) | test, **the convention's dedicated gate** |
| `test_psi_half_coupling.py:3237-3238` | — | ✅ (`len != 1` fail) | test, the carrying-arity gate |
| `test_iteration.py:788` `gains[0] + gains[1]` | ❌ (sum) | ✅ (=2) | test |
| the 5 `LC, (S, B) = …` test idioms | ✅ | ✅ (=2) | tests |

### B.4 Verdict — **can the rename be purely mechanical?**

**YES for the rename itself.** Nothing in the ordering/arity contract is encoded in the
*name* `gains`; it is encoded in the tuple's construction (`coupled_system.py:487`) and its
positional readers. Renaming the field is a pure token substitution at every one of these
sites.

**Two caveats to hand to the implementer:**

1. The field's **docstring at `:320-322` carries the convention.** If the rename rewrites
   that docstring block, the convention sentence (*"`B_a` LAST — the boundary-gain
   convention the G-S schedule arm parses"*) must survive verbatim. It is the only prose
   statement of a contract that three production lines depend on.
2. If the rename is ever *bundled* with a shape change (e.g. "make `explicit_gains` a
   `CoupledOperator` on both arms, dropping the tuple"), the blast radius is the **12 rows
   in §B.3 marked ORDER or ARITY** — 4 production, 8 test. That is a different, larger
   change. **Keep it out of this rename.**

---

## C. Docs surface — every `.rst` occurrence of the target field

Python-domain roles (`:attr:`/`:class:`/`:meth:`) render silently with **no `-W` warning**
unless the build runs `-n`, so the Sphinx gate will not catch a miss here (`coding-standards`
§retirement). None of the (T) doc sites below use a Python-domain *role* — they are all
double-backtick literals — which means they are **even more invisible**: they are plain
inline-literal text and no build check of any kind touches them.

### C.1 (T) — the field named in `.rst`. **8 sites, 6 files.**

| file:line | literal text | class | note |
|---|---|:--:|---|
| `docs/theory/foundations/coupled_block_operator.rst:521` | ``record's ``resolvent``/``gains``. The **iterative** splitting solve`` | **T ×2** | The canonical theory page for this record. |
| `docs/theory/foundations/coupled_block_operator.rst:638` | ``:func:`~orpheus.sn.coupled_system.build_within_group_system` ``.resolvent``` | **T** | Attribute access spelled in prose next to a live `:func:` role. |
| `docs/theory/foundations/boundary_conditions.rst:733` | ``Cartesian) mesh the record degrades to exactly this triple: its ``resolvent``` | **T** | ⚠ §1.3 collision #4 — same section carries 4 (C) uses. |
| `docs/theory/foundations/boundary_conditions.rst:736` | ``sweep) — and its ``gains`` are :math:`N = (S,\ B_a)`, the two lagged`` | **T** | Also states the ordering convention in prose. |
| `docs/theory/foundations/operator_adjoint.rst:504` | ``` ``resolvent`` ``(L+C)`` and lagged ``gains`` ``(S, B_a)`` keep ``S`` /``` | **T ×2** | ⚠ §1.3 collision #5. |
| `docs/theory/methods/sn/loss_representation.rst:1938` | ``within-group resolvent (``build_within_group_system(...).resolvent`` —`` | **T** | ⚠ mixed: "within-group resolvent" (C) then `.resolvent` (T) in ONE sentence. Prose rewrite. |
| `docs/theory/methods/sn/cartesian_multid.rst:3483` | ``applying its ``.resolvent``) — to return the user-facing`` | **T** | attribute access in a `:func:` cross-reference parenthetical. |
| `docs/theory/methods/sn/history.rst:75` | ``` ``.resolvent`` every driver consumes. **(6b, ``015dcc73``)** the``` | **T-hist** | ⚠ see C.2. |
| `docs/theory/methods/sn/history.rst:219-220` | ``:math:`A = M - N` (``resolvent`` = :math:`M`, the sweepable part / inverted each step; ``gains`` = :math:`N`, the lagged couplings)`` | **T-hist ×2** | ⚠ see C.2. |

### C.2 The history-page judgment call

`docs/theory/methods/sn/history.rst` is a **development-history changelog** — entries are
dated, hash-stamped records of what a past commit did. Two options, and the call is the
user's:

- **(a) Leave the history entries** (`:75`, `:219-220`) spelling `resolvent`/`gains`, because
  that is what the field was called when commit `015dcc73` landed, and **add a new entry** for
  this rename. Consistent with how the tombstones at `coupled_system.py:84`/`:295` preserve
  `_within_group_triple`/`_lagged_gains`.
- **(b) Rewrite them**, because the internal-consistency doctrine
  (`notation.rst:14-24`) demands **one spelling across code, tests, and prose** — a reader
  greping `implicit_operator` would find nothing in the history.

**Recommendation: (a) + a new history entry**, and make the new entry the crosswalk
(*"the field formerly spelled `resolvent`"*), so the grep lands. Option (b) falsifies a
hash-stamped record.

### C.3 `.claude/**` — 30 field-naming lines. **Archaeology; do NOT rewrite.**

| file | lines | class | verdict |
|---|---|:--:|---|
| `.claude/plans/sn_operator_realization_and_posing_plan.md` | `:627, 629, 636-637, 803, 950-951, 983, 1146-1147, 1236` | **T-plan** | **The plan proposing this rename.** Already spells the new names. Update its STATUS at merge, per the campaign-plan protocol. |
| `.claude/plans/sn_operator_realization_plan_REVIEW.md` | `:50, 127, 152-154, 350` | **T-plan** | The review that flagged the "74 test files" figure. **This audit corrects it to 12** — worth a one-line note there. |
| `.claude/plans/coupled_block_operator_campaign.md` | `:1479-1480, 1510, 2396` | **T-hist** | The campaign that *created* the record. Archaeology — leave. |
| `.claude/plans/dsa_landing_zone_recon.md` | `:19, 121, 146, 229, 234, 394` | **T-hist** | Recon snapshot. Leave. |
| `.claude/agent-memory/test-architect/coupled_operator_b2d_driver_eviction_verification.md` | `:3, 32, 35, 119, 249, 252, 446, 457` | **T-hist** | The gate spec that pinned `.resolvent`/`.gains`. Leave; but note `:252` records the `(B_a last)` convention — a second copy of the §B contract. |
| `.claude/agent-memory/test-architect/coupled_operator_step5_solve_verification.md` | `:48, 175, 363` | **T-hist** | Leave. |
| `.claude/agent-memory/test-architect/coupled_operator_4e_unweave_verification.md` | `:93` | **T-hist** | Leave. |
| `.claude/agent-memory/elegance-enforcer/coupled_block_boundary_unweld_rulings.md` | `:479, 512` | **T-hist** | Leave. |
| `.claude/agent-memory/archivist/lessons.md:1195, 1203` | | **T-hist** | An archivist lesson literally teaching *"(`.resolvent=(L+C)`, `.gains=(S,B_a)`)"*. **This one WILL mislead a future archivist.** Recommend a one-line update after the rename lands. |
| `.claude/agent-memory/archivist/task55_transport_sweep_docs_pass.md:11, 87` | | **T-hist** | Leave. |
| `.claude/agent-memory/explorer/sn_multigroup_axis_structure.md:41, 57` | | **O** | Names `WithinGroupSystem`/`build_within_group_system`, not the fields. No change. |
| `.claude/agents/explorer/AGENT.md` | 1 | **O** | Says `_within_group_triple` (already stale — the successor is `build_within_group_system`). Unrelated to this rename; already flagged in explorer memory. |
| `.claude/plans/operator_machinery_taxonomy.md:108, 966, 979, 1253` | | **C** ⚠ | **`A.resolvent(z)` — a RESERVED FACTORY NAME.** The #226 taxonomy ruled that `A.resolvent(z) = (A−zI).inverse()` is the future α-eigenvalue/shift-invert factory, with the first-resolvent-identity `R(z)−R(w) = (z−w)R(z)R(w)` as its family law. **Freeing the word `resolvent` from `WithinGroupSystem` is a PREREQUISITE for that factory** — this rename directly unblocks a ratified design. Worth citing in the commit message. |

---

## D. Name-collision check + naming recommendation

### D.1 Collision check — CLEAN

| candidate token | occurrences in `orpheus/`, `tests/`, `docs/` | verdict |
|---|--:|---|
| `implicit_operator` | **0** | free |
| `explicit_gains` | **0** | free |
| `M_operator` | **0** | free |
| `implicit` (as an identifier) | **0** | free — the ~35 hits are English prose (`implicit-reflective` BC, `implicit-zero` boundary, `implicit chain assumption`). No identifier. |
| `explicit` (as an identifier) | 0 | free |

The only tree-wide hits for `implicit_operator`/`explicit_gains` are the two plan documents
proposing them (`.claude/plans/sn_operator_realization_and_posing_plan.md`,
`…_plan_REVIEW.md`).

⚠ **But see Trap 1 (§0.4):** `explicit_gains` is safe as a *target*; the *source* token
`gains` is a substring of `against` (679 occurrences). Anchor the replacement.

### D.2 The vocabulary the tree ALREADY uses for `M`

Ordered by how load-bearing each spelling is:

| spelling | where | weight |
|---|---|--:|
| **"the sweepable part"** | `coupled_system.py:79`, `:310`, `:524`; `docs/…/coupled_block_operator.rst:455`; `docs/…/history.rst:219` | **the dominant ORPHEUS-native phrase — 5 independent sites** |
| **"regular splitting `A = M − N` (Hackbusch 2016 §11)"** | `coupled_system.py:78`, `:120`; `numerics/coupled_system.py:124`; `solver.py:172`, `:724`; `sweep_schedule.py:173`; `coupled_block_operator.rst:413-414`, `:444`; `boundary_conditions.rst:731-732`; `cartesian_multid.rst:67`, `:3809`; `history.rst:218-219` | **the canonical literature frame — 13 sites** |
| "the un-inverted FORWARD `M`" | `solver.py:821` | 1 — an apology for the current name |
| "the invertible **resolvent operand** `A`" | `notation.rst:299`; `numerics/iteration.py:15`, `:29`, `:34`; `numerics/operator.py:17` | 5 — **the numerics-layer binding** |
| "the implicit member `M`" | `.claude/plans/sn_operator_realization_and_posing_plan.md:629` | 1 — the plan only |
| **"lagged"** (for `N`) | the RETIRED helper `_lagged_gains`; `coupled_system.py:79` "the lagged coupling gains"; `:320` "the lagged couplings"; `solver.py:824`; `boundary_conditions.rst:684, 736`; `operator_adjoint.rst:504`; `operator_algebra.rst:190`, `:4121` | **the dominant native word for `N` — 9+ sites** |

**The documented double-binding.** `docs/theory/conventions/notation.rst:14-24` + row 8
(`:292-310`) states outright:

> where two layers must genuinely differ (the honest solver-level operator versus the
> **numerics-layer resolvent operand** — the last crosswalk row below), the binding is stated
> explicitly at **both** ends, so the difference is a documented bridge and not a drift.
> […] Not a canon row — the one place ORPHEUS itself carries two bindings of `A`, by design.

So `WithinGroupSystem.resolvent` is **not** an accident: it inherited the numerics layer's
"resolvent operand" vocabulary. **The rename must therefore update that notation row**, or
the crosswalk will point at a field that no longer exists. This is the single highest-value
doc edit in the change set, and it is the one a `.rst`-grep for `resolvent` will bury under
44 (C) hits in `operator_algebra.rst`.

### D.3 Is `implicit_operator` the best name? — recommendation

**Recommendation: `implicit_operator` / `explicit_gains` — ACCEPT, with one refinement.**

Reasons for:
- It names the **role in the splitting**, which is exactly what the field IS. `M` is solved
  *implicitly* (inverted) each step; `N` is evaluated *explicitly* from the lagged iterate.
  This is standard IMEX / operator-splitting vocabulary and it is the *reason* the splitting
  exists.
- The pair is **symmetric and greppable** — one word-order pattern `<splitting-role><object>`
  (`feedback_naming_consistency_greppable`: the role/qualifier token is always present, so
  grep finds the whole set). `implicit_operator` + `explicit_gains` reads as a set; so would
  `implicit_operator` + `explicit_operator` if `N` were ever folded to one operator.
- Zero collisions, and it **frees the word `resolvent`** for the ratified
  `A.resolvent(z) = (A−zI).inverse()` factory (§C.3, `.claude/plans/operator_machinery_taxonomy.md`).

Alternatives considered and why they lose:

| alternative | for | against |
|---|---|---|
| `sweepable` / `sweepable_part` | **highest ORPHEUS-native signal** (5 sites already say it); a SN reader knows instantly | **SN-specific.** The record is built on the method-generic `CoupledOperator` machinery, and the same splitting shape is the DSA/diffusion consumer's target. A future diffusion/CP `WithinGroupSystem` does not "sweep". Fails `feedback_extract_model_independent_primitive`. |
| `lagged_gains` (for `N`) | **9+ native sites use "lagged"**; the retired helper was literally `_lagged_gains` | (a) collides with the **tombstones** at `coupled_system.py:84`/`:295` — a grep for `lagged_gains` lands on a retirement note; (b) breaks the `implicit`/`explicit` symmetry; (c) "lagged" describes *when* it is evaluated, `explicit` describes *how the splitting classifies it* — the latter is the structural fact. |
| `forward_operator` | matches `solver.py:821` "un-inverted FORWARD M" | "forward" is heavily overloaded in this tree (forward vs **adjoint** transport, forward substitution). Would collide semantically with the `#276` adjoint work. |
| `M` / `N` | shortest, matches the math exactly | one-letter attribute names; unsearchable; the codebase already rejected this (the retired `self.L`/`self.S`/`self.F` triple is being deleted *right now* in the in-flight diff, partly for being "misnamed"). |

**The one refinement:** keep **"the sweepable part"** as the *SN gloss in the docstring*
(`coupled_system.py:310`), so the highest-signal native phrase survives where a SN reader
reads it, while the *field name* carries the model-independent structural role. That is
`feedback_high_signal_names` applied correctly: the precise term in the object's native
vocabulary, at the layer that owns it.

### D.4 The coherence decision the rename must ALSO make

`orpheus/numerics/iteration.py` gives **both driver primitives their own `.gains`
attribute**:

- `SourceIteration.__init__(self, A_inv, *gains, …)` → `self.gains = gains` (`:592`, `:623`),
  read at `:698`.
- `KrylovAcceleration.__init__(self, A, *gains, …)` → `self.gains = gains` (`:862`, `:882`),
  read at `:955`.

If `WithinGroupSystem.gains` becomes `explicit_gains` and the driver attribute stays `gains`,
the tree has **two spellings for the same object across one hop** (`SourceIteration(…,
*system.explicit_gains)` → `self.gains`). That is precisely what the internal-consistency
doctrine (`notation.rst:14-24`, the ratified prime directive) forbids.

Three coherent options — **user's call**:

- **(i) Rename the record field only.** Cheapest. Accept the seam, and *document* it the way
  notation.rst row 8 documents the `A` double-binding (state the bridge at both ends).
- **(ii) Rename both** — `SourceIteration.gains` / `KrylovAcceleration.gains` →
  `explicit_gains` too, and the variadic parameters `*gains` → `*explicit_gains`. Blast
  radius: `numerics/iteration.py` (~10 lines), `numerics/green_operator.py:274`,
  `solver.py:459/514/874/3274`, plus `tests/numerics/test_iteration.py`,
  `tests/numerics/test_green_operator.py`, `tests/sn/solve/test_si_single_primitive_contract.py`
  (its spy signature `_spy_init(self, L, *gains, **kwargs)`), and 4 `.rst` pages that spell
  `Driver(A, *gains)`. Roughly **doubles** the change.
- **(iii) Rename the record field to `explicit_gains` and leave the *numerics driver* alone
  deliberately**, on the grounds that the driver's `gains` is the **numerics-layer** name
  (the same layer whose `A` is the "resolvent operand") and the record's is the
  **solver-layer** name — i.e. exactly the two-layer bridge notation.rst row 8 already
  sanctions. This is (i) with a principled justification rather than an accepted debt.

I recommend **(iii)**, and updating notation.rst row 8 to state *both* bridges
(`A` ↔ `implicit_operator`, `*gains` ↔ `explicit_gains`) in the same row. That closes the
crosswalk instead of leaving it pointing at a dead name.

---

## E. What this audit did NOT find (gaps)

- **No dynamic access.** Zero `getattr(system, "resolvent")` / `getattr(…, "gains")` /
  `**asdict(system)` / `system.__dict__` sites. Grep for `getattr`+`resolvent` and
  `asdict`/`astuple` on `WithinGroupSystem`: empty. The rename cannot break a dynamic path.
- **No serialization.** `WithinGroupSystem` is never pickled, JSON-dumped, or written to a
  baseline file. No on-disk artifact carries the field names.
- **No public API exposure.** `orpheus/sn/coupled_system.py:155-160 __all__` exports the
  *class* and the *builders*, not the field names. No `orpheus/__init__.py` re-export.
  `docs/api/discrete_ordinates.rst` does not autodoc the attributes individually.
- **Nexus graph note.** `mcp__nexus__context` on
  `py:class:orpheus.sn.coupled_system.WithinGroupSystem` returns the class-level edges
  (6 doc pages, 9 implemented equations in `coupled_block_operator.rst`, 3 referencing
  functions) but **does NOT model dataclass fields as nodes** — only
  `py:attribute:…WithinGroupSystem.loss` exists, with degree 2. **The graph is blind to
  `.resolvent` and `.gains` as attributes.** This is the L-001 "graph callers alone
  under-scopes" pattern in its purest form: the graph surface contributed 0 of the 84 (T)
  sites. Text-grep was the primary evidence here, not a cross-check.
- **Sphinx gate will NOT catch a miss.** All 8 (T) `.rst` sites are double-backtick inline
  literals — not `:attr:`/`:class:` roles, not `:ref:`s. No build check of any kind (not even
  `-n` nitpicky) touches them. The `.rst` edits must be verified by grep, not by the build.

---

## F. Recommended execution order

1. `orpheus/sn/coupled_system.py` — the definition, both constructor call sites, the numpydoc
   entries, and the 3 prose sites. **Hand-edit `:401`** (collision #1). **Preserve the
   `B_a`-LAST sentence at `:320-322` verbatim.**
2. `orpheus/sn/solver.py` — **rebase onto the in-flight `self.L`/`self.S`/`self.F` retirement
   first**, then re-grep, then the 26 sites. Decide on the mirror-named locals
   (`_select_si_resolvent`, `base_resolvent`, `final_resolvent`) in the same pass or
   explicitly defer them.
3. The 12 test files (30 sites). Anchor every `gains` replacement (`\bgains\b`) — **679
   `against` landmines**.
4. The 8 `.rst` sites + **the notation.rst row-8 crosswalk** (§D.2) + the history-page
   decision (§C.2).
5. `.claude/agent-memory/archivist/lessons.md:1195-1203` — the one memory that will actively
   mislead a future agent.
6. Green gate: `tests/sn/operators/test_psi_half_coupling.py::TestWithinGroupSystem` is the
   dedicated record gate (arity/order/identity pins at `:3237-3269`) — it is the fastest
   proof the rename is complete and shape-preserving.
