# #343 + #344 — pickup memo, reconciled 2026-08-14

> **Purpose.** Both issues were filed **2026-08-09** and are about to be
> re-evaluated. This memo is the reconciliation `plan-authoring` §7 requires
> *before* designing to them: what still holds, what may have rotted, and what
> to measure first. It is deliberately short — the issue bodies are excellent
> and are the content; this is only the staleness audit around them.
>
> ⚠ Read this **before** `gh issue view 343 / 344`, not after.

---

## The pointers, existence-checked `[M]` 2026-08-14

Everything both issues name **exists**, and — unusually for `scratch/` — all
four instruments are **tracked**:

| artefact | state |
|---|---|
| `scratch/probe_341_iteration_spectrum.py` | EXISTS, tracked, 7026 B |
| `scratch/probe_341_fold_pattern.py` | EXISTS, tracked, 3715 B |
| `scratch/probe_341_octant_model.py` | EXISTS, tracked, 6970 B |
| `scratch/issue_341_gs_jacobi_mechanism.md` | EXISTS, tracked, 51073 B — the analysis of record |
| `SweepSchedule.gauss_seidel` | `orpheus/sn/loss_representation/sweep_schedule.py:144` (class `:122`) |
| `Quadrature.octants` | `orpheus/numerics/quadrature/directional.py:541` |
| `build_within_group_system` | `orpheus/sn/coupled_system.py:425` |
| the splitting-algebra `.. warning::` | `orpheus/sn/coupled_system.py:112` |

⭐ Six **further** `probe_341_*` files exist that neither issue cites —
`antipodal_prediction`, `mode_shape`, `octant_order`, `predicate_sweep`,
`scattering_ratio`, plus `issue_341_docs_sweep.md`. Look there before rebuilding
any instrument.

⚠ **They survived a mass-delete by one day.** `977615f7` (2026-08-09,
*"retire 14 superseded diagnostics; repoint the 2 that are instruments"*) landed
the same day these issues were filed. The `probe_341_*` set came through intact,
but that is the `coding-standards` hazard in the flesh — if any future sweep
targets `scratch/`, these four are **production infrastructure**, not debris.

## ⛔ `#341` is CLOSED — it is a RECORD, not a work item

`gh issue view 341` reads CLOSED (*"Boundary Gauss-Seidel is ~2x SLOWER than
Jacobi at d=3 — the documented rate gain was measured at d=2 and inverts"*).
Both #343 and #344 cite it as the origin of their measurements. Do **not** read
its closed state as "this line of work is done" — it is `plan-authoring` §9(b):
the issue that *characterised* the mechanism, spun out into #343 (the rate
lever) and #344 (the singularity).

---

## ⚠ THE ONE PREMISE MOST LIKELY TO HAVE ROTTED

**#344's headline conjoins a VALUE claim and a TREE-BEHAVIOUR claim, and only
the second has a clock that has been running.** It reads:

> an **11.26 % trace deviation** that the SI residual stopping test
> **structurally cannot see**.

The `11.26 %` is numerics and survives until the numerics move. *"the stopping
test cannot see it"* is a statement about code, and **the #340 convergence
campaign landed directly on that code after the issue was filed**:

| commit | date | what it changed |
|---|---|---|
| `6cb5e519` | 2026-08-10 | the convergence warning speaks for the level that FAILED |
| `28435e11` | 2026-08-10 | a solve standing on a starved inner is no longer silent |
| `4aff0d4e` | 2026-08-10 | the truncation warning says how much it probably cost — **and added the per-group balance-defect projection** |
| `5c4d15f7` | 2026-08-11 | a starved solve announces itself at EVERY public entry |

`4aff0d4e` is the pointed one: it introduced the **balance-defect** reporting
(`‖R_g‖/‖Q_g‖`), which is a *different functional from the SI residual* and is
exactly the kind of quantity that could now observe a trace deviation the
residual cannot. The wide-gate log at HEAD shows it live, e.g.
`ConvergenceWarning: ... leaves a per-group balance defect of ‖R_g‖/‖Q_g‖ = 2.632e-10`.

⟹ **Re-run `probe_341_iteration_spectrum.py` at HEAD before designing.** Split
the headline into two claims and check only the behaviour half first — it is the
cheap one, and if the balance projection *does* see the deviation, #344's framing
("a symptom the stopping test cannot see") changes, and possibly one of its two
readings is already decided. This is `plan-authoring` §2's shelf-life clause with
its own date collision: **filed 08-09, behaviour-changing commits 08-10 and 08-11.**

## ⚠ #343's second-order claim needs the same treatment

> nothing gates it — a future change to `Quadrature.octants`' enumeration would
> silently move the SI rate by up to 2.5×, **with every existing test still
> green** ... and **no gate reads cost**.

Also a tree-behaviour claim, also filed 08-09. Cheap to re-check: is there now
any gate that reads iteration COUNT (as opposed to the converged value)? The
#340 campaign added `IterationRecord` / `fully_converged` machinery and
`IterationBudget` (`0f5ca91c`), so *something* may now observe cost.
⟹ one grep over `tests/` for assertions on iteration counts settles it.

## What is safe to trust

- **#343's spanning measurement** (`n_GS/n_Jacobi ∈ [0.771, 1.892]`, 25/25
  predicate separation, shipped order 24th of 25). Pure linear algebra on the
  assembled iteration matrix; nothing since touched the assembly.
  ⚠ Its configuration includes **level-symmetric S4**, and #337 (`59bb38a0`,
  2026-08-08) changed the LS seed to moment-matched — but that **predates** the
  measurement, so the numbers are already post-#337. No action.
- **#343's scope note**: the fixed point is splitting-invariant (G-S `2.79e-14`,
  Jacobi `4.70e-14`, Krylov `6.56e-13` agree). Nothing is *incorrect* today —
  this is a rate/ownership issue, not a bug.
- **#344's discriminators**, all four, are still the right experiments. In
  particular the cleanest one is unchanged: **a product quadrature with no
  tangential ordinate** should make the singularity vanish entirely under
  reading (1).

## Suggested first moves (cheap → expensive)

1. `gh issue view 341` + skim `scratch/issue_341_gs_jacobi_mechanism.md` — the
   analysis of record for both.
2. Re-run `probe_341_iteration_spectrum.py` at HEAD. Settles #344's behaviour
   half and re-establishes #343's baseline in one shot.
3. Grep `tests/` for iteration-count assertions → settles #343's "no gate reads
   cost".
4. Only then design.

⚠ **Do not start by extracting null vectors** (#344's first bullet). It is the
right experiment and it is not the cheapest one; step 2 may reframe the question
first.

---

## Related open work, for scope awareness

`#369`–`#372` were filed 2026-08-14 from the Q6 campaign and are unrelated to
these two. `#340`'s remainders (`#343`, `#344`, `#348`, `#351`–`#353`,
`#363`–`#366`) are the same family as these; the plan of record for that campaign
is `.claude/plans/archive/nested_iteration_diagnostics.md`.
