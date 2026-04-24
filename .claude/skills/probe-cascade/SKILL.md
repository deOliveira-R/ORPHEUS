---
name: probe-cascade
description: 'Proactively use when a solver/model gives a WRONG-BUT-PLAUSIBLE answer (e.g., 1% off instead of 100%) and the suspect region has many interacting factors (group count × region count × closure × kernel). Formalises the systematic isolation sequence — drop one complication at a time until the minimal reproducer is identified, then audit that minimal configuration for closed-form structure the current code is missing. Examples: "solver disagrees with reference by 1-10%", "K matrix is wrong but not obviously", "convergence plateaus above machine precision", "matches in simple case, fails in compound case".'
---

# Probe Cascade — systematic isolation for plausible-wrong numerical bugs

Named after the Issue #131 investigation (2026-04-24) where a 1.5 %
`k_eff` gap between two slab drivers was narrowed from "multi-region
multi-group × F.4 closure × adaptive-mpmath" to a single
finite-N-GL-where-closed-form-applies bug in
`compute_P_esc_{outer,inner}` via a 5-probe isolation cascade.

## When to use

- Solver disagrees with a trusted reference by a small-but-nonzero
  amount (not obviously catastrophic — e.g. 1e-3 to 1e-1).
- The suspect code path has **multiple independent interacting
  factors** (group count × region count × closure choice × kernel
  type × boundary choice × geometry variant).
- A subset of factor combinations already pass a tight parity test;
  the specific combination that fails has no direct parity test yet.

If the bug is catastrophically wrong (NaN, negative k_eff, off by
orders of magnitude), use direct traceback / `nexus-debugging`
instead — the probe cascade is for subtle bugs in a mostly-correct
solver.

## Core technique — drop one dimension at a time

Enumerate the factors that distinguish the failing case from the
passing one. Run a probe for each factor in isolation:

| Axis         | Passing case | Failing case | Probe                               |
| ------------ | ------------ | ------------ | ----------------------------------- |
| Group count  | 1G           | 2G           | Run 1G at the failing config        |
| Region count | 1 region     | 2 regions    | Run 1 region at the failing config  |
| Closure      | vacuum       | `white_f4`   | Run vacuum at the failing config    |
| Kernel       | `e^(-τ)`     | `Ki₁`        | Compare at matched inputs           |
| BC           | vacuum       | reflective   | Remove BC, check volume kernel only |

Each probe should **change only one factor** vs the passing case.
The first probe that fails localises the bug to that factor (or
its interaction with the other already-confirmed factors).

### The 5-probe structure (minimum)

1. **Probe A — drop factor 1**: isolate the volume-kernel / scatter
   / fission arithmetic from the highest-level factor (MG, BC,
   closure).
2. **Probe B — drop factor 2**: same reasoning one level down.
3. **Probe C — drop factor 3**: continue until factors are pure
   volume-kernel arithmetic.
4. **Probe D — pin the arithmetic**: at the isolated level, build
   a convergence sweep in the likely-offender parameter (quadrature
   node count N, dps precision, mesh refinement). If the error
   plateaus, there's a structural bug (not quadrature
   underconvergence). If it converges, the integral was
   underconvergent — check for a closed form.
5. **Probe E — validate the fix**: write the candidate fix as a
   pure-function probe, compare vs the reference at the failing
   configuration. Must reduce the gap to the reference's precision
   class (typically 1e-10 to machine epsilon).
6. **(Optional) Probe F — end-to-end confirmation**: re-run the
   top-level failing configuration with the fix integrated. Must
   see the same rel_diff reduction as Probe E.

## Directory and filename convention

All probes live in `derivations/diagnostics/diag_{issueN}_probe_{letter}_{descriptor}.py`:

- `diag_slab_issue131_probe_a_1g_2rg_vacuum.py`
- `diag_slab_issue131_probe_b_2g_2rg_vacuum.py`
- `diag_slab_issue131_probe_d_pesc_quadrature.py`
- `diag_slab_issue131_probe_e_closed_form_fix.py`
- `diag_slab_issue131_probe_f_2eg_2rg_f4.py`

Each probe is a standalone pytest test — can be promoted to
permanent tests per `tests/derivations/_promotion_policy.md`.

## Closed-form detection

**Strong bias: when a finite-N GL integral's error plateaus above
machine precision, look for a closed-form integral identity first,
before blaming quadrature convergence rate or adding nodes.**

Pattern that commonly hides closed forms:

- `∫₀¹ exp(-τ/µ) dµ = E₂(τ)` — **any** angular integral with
  µ-independent τ over a factorisable separable integrand.
- `∫₀¹ µⁿ exp(-τ/µ) dµ = Σ Eₙ₊₂`-ish — from integration by parts.
- `∫₀^π sin(θ) exp(-τ/sin(θ)) dθ = 2 Ki₁(τ)` — Bickley-Naylor.
- Polynomial source × exponential kernel = `E_n` moment series
  (see the Peierls moment-form theory for the slab case).

If the integrand has form `f(µ) · exp(-τ(...)/µ)` with τ
independent of µ, **assume a closed form exists** and search for it.
Finite-N GL is only appropriate when τ depends on µ in a way that
prevents factorisation (typically curvilinear chord geometry).

## Output — what to hand back

After the cascade:

1. **The minimal reproducer** — which probe first failed.
2. **The root cause** — specific file:line of the defective code
   and the mathematical explanation.
3. **The fix** — as a diff (patch or edited file), with a re-run
   confirming the top-level parity test now passes at the
   reference precision.
4. **Promotion recommendation** — which probe tests should move to
   `tests/` as permanent regression gates (follow
   `tests/derivations/_promotion_policy.md`).
5. **Related-pattern audit** — grep for the SAME anti-pattern
   elsewhere. If the bug is "finite-N GL for a closed-form
   integral," there may be twins in sibling functions.

## Anti-patterns

- **Don't just bump N.** If the error plateaus, more quadrature
  nodes won't help — go look for the closed form.
- **Don't skip probes.** The temptation is to jump straight to
  the suspected bug. Running the probes in order prevents
  chasing a wrong hypothesis.
- **Don't run full production quadrature in probes.** Use the
  minimum N / dps that reproduces the gap. Full precision is for
  the final validation run.
- **Don't commit probes before triage.** The `_promotion_policy.md`
  decides whether each probe becomes a permanent test or gets
  deleted (the git history retains the investigation trail).

## Reference implementation

Issue #131 (2026-04-24): slab multi-region 2G 1.5 % gap →
5.4 × 10⁻¹⁶ bit-exact after closed-form `½ E₂(τ_total)` replaced
a finite-N GL branch in `compute_P_esc_{outer,inner}` /
`compute_G_bc_{outer,inner}`. See commit `3b0b2c9` and the Sphinx
section `§theory-peierls-slab-polar-g5-diagnosis` for the full
trail, including why each probe fired or didn't, and the specific
closed-form derivation.
