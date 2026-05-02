---
name: algebra-of-record
description: 'PROACTIVELY load when implementing a new numerical method, building a verification-grade reference solver, or reviewing whether a derivation→code pipeline meets project discipline. Codifies the SymPy-as-canonical-algebra pattern: derivations bifurcate at a specialized form into Branch 1 (continuous reference: SymPy, SymPy+mpmath, or MMS) and Branch 2 (discretized production: numpy/scipy). Branch 2 is verified at L1 by cross-check against Branch 1. Sphinx narrates both, with the SymPy module as the canonical algebra-of-record. Read by sub-agents that build, document, or verify reference solvers (method-implementer, archivist, test-architect) and by the main agent when orchestrating that pipeline. Examples: "Where should the SymPy module go", "When do I bifurcate the derivation", "Branch 1 and Branch 2 disagree at L1 — what now", "Why does my SymPy derivation hang on multi-group", "What goes in the Sphinx stub vs the rich narrative".'
---

# Algebra-of-Record — derivation bifurcation discipline

This skill codifies the project's vision for grounding numerical
methods in symbolic mathematics: **SymPy is the canonical
algebra-of-record from which both the verification reference and the
production code descend, with Sphinx documentation as the narrative
that ties them together.**

The discipline maps onto the V&V framework's three pillars (see
`vv-principles`): every reference solver is built via one of three
**Branch-1** strategies (closed-form SymPy / semi-analytical
SymPy+mpmath / MMS), and the production solver lives in **Branch 2**
(numpy/scipy/numba) where it is verified at L1 by structurally
independent cross-check against Branch 1.

This is preloaded by:

- The **method-implementer** sub-agent (builds Branch 2 + writes
  Branch 1 SymPy module + Sphinx stub).
- The **archivist** sub-agent (turns Sphinx stubs into rich
  narrative, with the SymPy module as canonical reference).
- The **test-architect** sub-agent (designs the L1 cross-check
  between branches).
- The **main agent** (orchestrates the pipeline; understands the
  SymPy chokes that constrain how high-dimensional cases get
  verified).

---

## CRITICAL: One-line summary

A derivation starts with basic premises (transport equation, BC, XS
conventions). SymPy carries the algebra forward through specialization
until it reaches a **bifurcation point** — the form where the math
admits both:

- A continuous evaluation (Branch 1: closed-form SymPy / SymPy+mpmath
  / MMS-derived) usable as a verification reference, AND
- A discretization (Branch 2: numpy/scipy production code) usable
  for actual analysis.

After bifurcation, both branches share the same SymPy ancestor, so
disagreement at L1 cross-check exposes a discretization bug, not a
math bug.

---

## CRITICAL: The bifurcation pattern

```
                    transport equation + BC + XS
                              │
                       (SymPy carries algebra)
                              │
                              ▼
                       SPECIALIZED FORM
                       (the bifurcation point)
                              │
                ┌─────────────┴─────────────┐
                ▼                           ▼
            BRANCH 1                    BRANCH 2
       Continuous reference         Discretized production
            │                           │
   Closed-form SymPy               numpy / scipy
   SymPy + mpmath                  vectorized loops
   MMS-derived source              GL/GJ quadrature
            │                           │
            ▼                           ▼
   tests/derivations/              orpheus/<module>/
   test_<name>_symbolic.py         <name>.py
   (foundation-tagged tests)       (production code)
            │                           │
            └─────────┬─────────────────┘
                      ▼
              L1 cross-check test
              tests/derivations/
              test_<name>_xverif.py
              (compares Branch 2 to Branch 1)
                      │
                      ▼
              Sphinx theory page
              docs/theory/<topic>.rst
              (narrates BOTH branches with
              :mod: cross-refs to SymPy module)
```

The skill enforces:

1. **Where to bifurcate** (judgment call, see "The bifurcation point"
   section below).
2. **Folder layout** (Branch 1 vs Branch 2 vs tests vs Sphinx).
3. **Verification function shape** (`derive_*() -> dict`).
4. **Structural-independence discipline** (above vs below the
   trusted-library line).
5. **The minimal-SymPy + scaling-argument pattern** for problems
   SymPy chokes on at full dimension.
6. **Sphinx stub-vs-narrative separation** (method-implementer
   writes stubs; archivist fills with rich narrative).

---

## Branch 1 — the three states of "continuous reference"

Branch 1 is NOT just "pure SymPy". It has three states corresponding
to the three V&V pillars (`vv-principles` § "The three pillars of
verification"):

### State 1A — Closed-form SymPy

The math closes algebraically. `derive_*()` returns symbolic
identities verifiable by `sp.simplify(lhs - rhs) == 0`.

**Use when**: the operator algebra reduces to expressions in
elementary functions or special functions SymPy fully understands
(`exp`, `cosh`, `E_1`, Bessel of small order, polynomial roots up
to degree 4).

**Examples in this project**:

- V_α2 (Plan 2 B3): `T_00^sphere = P_ss^sphere` algebraic identity.
  Both sides reduce to `(1 - (1+2τ_R) e^{-2τ_R})/(2τ_R²)`.
- V_α1 (Plan 2 B3): closed-sphere bounce-sum self-consistency
  `ψ_surf = q/Σ_t` independent of bounce-period chord length.

**Failure mode**: trying to push closed-form too far. If `simplify()`
takes more than ~30 seconds or returns visibly worse-looking
output, you're in the wrong state — drop to State 1B or 1C.

### State 1B — Semi-analytical (SymPy + mpmath)

The math reduces to a single integral that SymPy cannot evaluate
symbolically, but `mpmath.quad` can evaluate to arbitrary precision.

**Use when**: the operator structure is symbolic but the final
quadrature is numerical.

**Examples in this project**:

- PS-1982 vacuum-sphere reference (`ps1982_reference.py`): SymPy
  derives the kernel `[E_1(|r-x|) - E_1(r+x)]`, mpmath-flavored
  scipy.integrate.quad evaluates the Nyström matrix.
- Sanchez 1986 Eq. (A6) for sphere specular: SymPy verifies the
  integrand structure, numerical quadrature evaluates the µ-integral.

**Two-step correctness ladder** (per `vv-principles`
§ "Semi-analytical correctness ladder"):

1. **Integrator correctness**: assumed for tested upstream
   (`scipy.integrate`, `mpmath.quad`); required to be verified for
   custom integrators.
2. **Reduction correctness**: the math from equation to single
   integral is the load-bearing piece. If the reduction is wrong,
   the integral is exact for the wrong equation.

If the integrator is from a trusted upstream library, **the same
integrator can be shared between Branch 1 and Branch 2** — see
"Structural independence applies above the trusted-library line"
below.

### State 1C — MMS (Method of Manufactured Solutions)

When no analytical/semi-analytical reference exists, manufacture a
trial solution `ψ_chosen(r, µ, ...)`, derive the source `Q^ext`
that makes it true symbolically, and use the resulting (ψ, Q^ext)
pair as the reference.

**Use when**: the operator has no closed-form / semi-analytical
reference (typical for high-dimensional, multi-region, anisotropic
configurations).

**MMS operational rules** (from `vv-principles` § "MMS operational
rules"):

- Trial MUST not vanish under derivatives (avoid polynomials).
- Trial MUST be non-trivial at boundaries (test BC handling).
- Trial MUST stress-test the numerical method, not minimize source
  complexity.
- Manufactured source MUST be structurally independent of the code's
  primitives.

**MMS does NOT prove eigenvalues** — it's a source-driven test of
the operator's spatial / angular discretization, not the eigenvalue
iteration. Eigenvalue claims need State 1A or 1B.

### Choosing the state

```
Q: Does the math close in elementary/special functions? → 1A
   Are simplify() calls finishing fast (< 30 s)? Are results
   recognizable? Is the closed form actually shorter than typing
   the result?
       Yes → 1A.
       No  → drop to 1B.

Q: Does the math reduce to a single integral with bounded
   integrand?
       Yes → 1B (SymPy + mpmath/scipy).

Q: No reference exists (or 1A/1B both fail)?
       → 1C (MMS), accepting that the verification claim is
       flux-shape and convergence-order, NOT eigenvalue.
```

---

## Branch 2 — the production discretization

Branch 2 is `numpy` / `scipy` / `numba` production code in
`orpheus/<module>/<name>.py`. It is the code users actually run for
analysis. It descends from the same SymPy ancestor as Branch 1 but
uses discrete data structures (arrays, sparse matrices, grids).

**Architectural rule**: Branch 2 implements the SAME operator algebra
as Branch 1 — only the evaluation strategy differs. If Branch 2's
operator construction is structurally different from Branch 1's
(e.g., uses a different equation), you've bifurcated at the wrong
point.

**Folder layout**:

```
orpheus/<module>/<name>.py              # production solver
orpheus/derivations/<module>/origins/<topic>/<name>.py
                                        # Branch 1 SymPy module
                                        # OR
orpheus/derivations/<module>/<name>_reference.py
                                        # Branch 1 mpmath ref solver
                                        # (when 1B and the ref is
                                        # a full solver, not just
                                        # symbolic identities)
```

---

## The bifurcation point — the load-bearing judgment call

**Bifurcate too early** (before the math reaches its specialized
form) → Branch 1 and Branch 2 solve different things; L1 cross-checks
fail not because of bugs but because the branches don't match.

**Bifurcate too late** (after Branch 1 has done a lot of
"production-style" work) → Branch 2 becomes a near-clone of Branch 1
with shared failure modes; L1 cross-checks pass but mean nothing.

The right bifurcation point is the **first form where**:

1. The continuous form admits closed-form / semi-analytical
   evaluation (Branch 1 is realizable in State 1A or 1B).
2. The discrete form admits standard quadrature/iteration without
   information loss (Branch 2 is realizable as numpy code).

For transport problems specifically, this is typically:

- **AFTER** specialization to geometry (sphere/cylinder/slab) and
  boundary condition (vacuum/specular/white).
- **AFTER** angular reduction (azimuthal integration for
  axially-symmetric, polar integration for slab).
- **BEFORE** any quadrature-rule choice (Gauss-Legendre vs
  Gauss-Jacobi; per-segment vs single-domain).

The SymPy module ends at the bifurcation point. Branch 1 evaluates
the result symbolically/semi-analytically. Branch 2 picks up where
SymPy left off and discretizes.

---

## Verification function pattern

Branch-1 SymPy modules expose verification functions with this shape:

```python
def derive_<identity_name>() -> dict:
    """V_<n> — <identity description>.

    Proves: <math statement, in plain English>.

    Returns dict with:
      - name: verification name (string)
      - <symbolic expressions used>
      - pass: bool (True iff identity verified)
    """
    # SymPy work
    x, y = sp.symbols("x y", positive=True, real=True)
    lhs = ...
    rhs = ...
    diff = sp.simplify(lhs - rhs)
    pass_v = (diff == 0)
    return {
        "name": "V_<n>: <identity>",
        "lhs": lhs,
        "rhs": rhs,
        "diff": diff,
        "pass": pass_v,
    }
```

Test gate at `tests/derivations/test_<name>_symbolic.py`:

```python
@pytest.mark.foundation
def test_v_<n>_<identity>():
    result = derive_<identity_name>()
    assert result["pass"], f"V_<n> failed: {result}"
```

Multi-identity SymPy modules expose multiple `derive_*` functions;
the test file pins one foundation-tagged test per `derive_*`. The
test file's test count = the V_n claim count.

The `@pytest.mark.foundation` tag indicates "software invariant —
verifies the SymPy ground truth, not an L0/L1/L2 claim about a
solver". The V&V audit harness recognizes the tag and does not look
for matching equation labels.

---

## Structural independence applies ABOVE the trusted-library line

L1 cross-check between Branch 1 and Branch 2 is valid only if the
two branches are **structurally independent** — they cannot share
upstream identities or primitives where a bug could hide identically.

**Above the trusted-library line** (math layer): the reduction, the
discretization choice, the operator construction. Branch 1 and Branch
2 must NOT share work here. Sharing produces ERR-032-style
catastrophes (two derivations agreeing because both inherit the same
wrong identity).

**Below the trusted-library line** (numerics primitives): integrators,
quadrature rules, special functions, linear-algebra primitives.
Sharing IS OK if the shared library has independent validation:

- `scipy.integrate.quad`, `scipy.special.exp1`, `scipy.linalg.eig`,
  `numpy.polynomial.legendre.leggauss` — trusted upstream, validated
  by decades of use across the scientific Python ecosystem.
- `mpmath.expint`, `mpmath.quad` — trusted upstream, arbitrary
  precision.
- Project-internal primitives (`orpheus.derivations.common.*`) —
  trusted only if they have their own L1 verification chain in this
  project.

**Concrete rule**: Branch 1 and Branch 2 may share a `scipy.special.exp1`
call. They MUST NOT share a project-internal `compute_my_kernel()`
function unless that function has its own verification chain.

This logic mirrors the V&V framework's two-step semi-analytical
correctness ladder (`vv-principles`): integrator correctness is the
foundation pillar; the reduction correctness is what the
cross-check measures. If the integrator is verified, sharing it is
fine.

**Anti-pattern (ERR-032)**: two SymPy derivations agreed at 1e-39
because both used `∫E_2 = 1 − E_3` (wrong; correct is `½ − E_3`).
The cross-check was procedurally independent (different code paths)
but not structurally independent (same upstream identity). The bug
was caught by a third semi-analytical reference that didn't
inherit the identity.

---

## Minimal-SymPy + scaling argument

SymPy chokes at high dimensionality (see "SymPy choke modes" below).
For realistic transport problems (multi-group + multi-region +
anisotropic), pure-symbolic verification is infeasible. The
discipline:

1. **Verify symbolically on the simplest non-trivial case.** For
   transport: 1-group, 1-region, isotropic-scattering, simplest BC.
   This case captures the operator algebra structure.
2. **Verify Branch 2 against Branch 1 on that simple case**
   (numerical L1 cross-check).
3. **Argue that the operator algebra extends mechanically** to
   higher dimensions. The argument has the form:
   - "The operator construction in Branch 2 uses array operations that
     are dimensional generalizations of the symbolic operations in
     Branch 1."
   - "The discretization choice (e.g., Gauss-Legendre on the same
     interval) is identical at all dimensions."
   - "The only thing that changes with dimension is array shape and
     the number of group/region indices."
4. **Test higher dimensions numerically** against:
   - MMS sources (State 1C) for spatial/angular accuracy.
   - External benchmarks (Garcia 2021, Sood 2003 LANL benchmarks)
     for eigenvalue / flux-shape.
   - Dimensional reductions (e.g., G=1 special case of G≥2 solver)
     bit-for-bit.

The scaling argument is NOT proof — it's a justified extrapolation.
The L1 evidence at the simple case combined with bit-equal reductions
+ external benchmarks at higher cases provides the verification chain.

**Practical implication**: a method-implementer designing a 4G/4R
multi-region prototype should NOT try to write the SymPy derivation
for the full 4G/4R case. They should write the SymPy derivation for
the 1G/1R case, verify it's algebraically correct, then implement
Branch 2 for general G/R and test the G=1/R=1 reduction against
SymPy + the higher-dimensional cases against external references.

---

## SymPy choke modes — a catalog

Knowing where SymPy fails saves you from spending hours fighting it.

### 1. Expression-tree growth is super-linear

Every operation grows the expression tree. `simplify()`,
`integrate()`, `solve()` are O(complex) in tree size.

**Symptom**: A `simplify()` call that completes in 1s for G=2 takes
30s for G=4 and hangs at G=8.

**Workaround**: Don't simplify globally. Use `expand_log`,
`together`, or targeted `subs` to make local progress. Or work in
the smallest dimension that captures the structure.

### 2. Eigenvalue solving has a hard mathematical ceiling

For G×G matrices with symbolic entries, dominant-eigenvalue requires
solving a degree-G characteristic polynomial:

- G=2: quadratic formula. Clean closed form. ✓
- G=3: Cardano's formula. Ugly nested radicals; `simplify` may not
  fully reduce.
- G=4: Ferrari's formula. Near-unreadable; barely useful.
- **G ≥ 5**: NO closed form by Abel-Ruffini theorem. SymPy returns
  `RootOf(...)` placeholders — algebraic numbers defined only
  implicitly. Useless for downstream symbolic manipulation.

**Workaround**: For G ≥ 3, don't try symbolic eigenvalues. Use a
2-group symbolic case for the operator-algebra verification, then
verify higher-G numerically against `np.linalg.eigvals` of the
transfer matrix (e.g., `kinf_and_spectrum_homogeneous` in
`orpheus.derivations.common.eigenvalue`).

### 3. Multi-region piecewise compounds badly

`Piecewise` works for single-piece definitions. Composing two
piecewise functions (piecewise σ_t along a piecewise trajectory)
produces nested `Piecewise` objects whose simplification rules are
incomplete — `integrate()` often returns unevaluated `Integral(...)`.

**Workaround**: Don't push multi-region work through SymPy. Verify
the single-region case symbolically; argue that multi-region
extends via piecewise τ(µ) + composite quadrature; test numerically
against external reference.

### 4. Anisotropic scattering doubles complexity per Legendre order

Each Legendre moment introduces a new tracked angular flux moment +
a coupling equation. The symbolic system grows quadratically with
maximum order.

**Workaround**: Verify L=0 (isotropic) symbolically; argue the
operator extends mechanically; verify L≥1 against MMS or external
benchmarks numerically.

### 5. Performance kills the dev loop

A `simplify()` call that takes 30 seconds means you can't iterate.
You write a derivation, wait 30s, see it failed, fix, wait 30s,
loop. Productivity collapses.

**Workaround**: Pre-emptively avoid the high-cost operations. Don't
call `simplify()` at the end of every derivation step. Use targeted
operations (`expand`, `factor`, `cancel`, `subs`) and only
`simplify()` at the verification-comparison step.

---

## Sphinx stub vs rich narrative — separation of concerns

The Sphinx theory page is the **narrative** that ties Branch 1,
Branch 2, and the L1 cross-check into a coherent mathematical story.
It has TWO modes:

### Stub mode (method-implementer's deliverable)

Written by the agent that builds the prototype. Contains:

- A `:label:` for each verifiable claim (one per `derive_*()`
  function in the SymPy module).
- A `:mod:` cross-reference to the SymPy module.
- A 1-paragraph TODO marker for each label, of the form:

  ```rst
  .. _peierls-greens-V-alpha-1:

  V_α1 — closed-sphere k_inf identity
  ====================================

  .. todo:: Archivist expansion needed.
     The SymPy derivation lives in
     :mod:`orpheus.derivations.continuous.peierls.origins.specular.greens_function`
     (function ``derive_operator_constant_trial_closed_sphere``).
     Test gate:
     :func:`tests.derivations.test_peierls_greens_function_symbolic.test_v_alpha1_overall_pass`.
     Closeout memo:
     ``.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md``.

     Brief: V_α1 algebraically proves
     :math:`(K \cdot \mathrm{const})(r,\mu) = \omega_0 \cdot
     \mathrm{const}` for closed homogeneous sphere with specular BC.
     Three-step proof: (1) surface fixed-point; (2) total-ψ-constant
     across L_first; (3) operator eigenvalue ω_0.
  ```

The stub is what the method-implementer SHIPS. It is sufficient for
the V&V audit harness to recognize the claim, for Nexus to index
the math node, and for downstream code to `:eq:`-cite the label.

### Rich narrative mode (archivist's deliverable)

Written by the archivist, expanding the stub. Contains:

- Full mathematical derivation in math + prose, walking the reader
  through the algebra.
- All `derive_*()` identities written out with their proofs.
- Numerical evidence (cross-verification matrix; convergence tables).
- Design rationale (why this Branch-1 state, why this discretization).
- What was tried and failed (with literature references).
- Cross-document references to related theory pages.

The rich narrative is what Cardinal Rule 3 (Sphinx is the LLM's
brain) calls for. It is a **maximum-effort task** for the archivist.

### Why the separation matters

Conflating the two creates pressure to write the rich page during
implementation, which is wrong:

- Different cognitive modes (synthesis vs implementation).
- Different agent's expertise (archivist has Sphinx + V&V context;
  method-implementer has algebra + numerics context).
- Different time horizons (stubs are quick; rich pages are slow).

**Discipline**: method-implementer produces stubs; archivist
expands stubs. Both reference the same SymPy module as the
algebra-of-record. The transition from stub to rich narrative happens
via a `DISPATCH_REQUEST` to the archivist (see
`subagent-handoff-protocol`).

---

## Anti-patterns

### For method-implementers

- **NEVER skip the SymPy module** "to save time". The SymPy module
  IS the canonical algebra-of-record. The Sphinx page references it.
  The test gate pins it. The L1 cross-check trusts it. Without it,
  there is no verification chain — only "two pieces of code agree",
  which is L4 not L1.
- **NEVER bifurcate before the specialized form.** Branch 1 and
  Branch 2 must descend from the SAME specialized SymPy ancestor.
  Bifurcating earlier means they solve different things.
- **NEVER push SymPy past its choke modes.** If `simplify()` hangs,
  you're in the wrong state. Drop to State 1B (mpmath) or 1C (MMS).
- **NEVER share project-internal primitives between Branch 1 and
  Branch 2.** Sharing trusted-library primitives is fine; sharing
  in-house code creates ERR-032-style hidden dependencies.
- **NEVER write the rich Sphinx narrative yourself.** Stub it and
  dispatch the archivist.

### For archivists

- **NEVER expand a stub without reading the SymPy module first.**
  The SymPy module is the canonical source. Your narrative narrates
  it; it doesn't compete with it.
- **NEVER edit the SymPy module yourself.** If you find an algebra
  error while expanding the stub, return to the user with a
  `DISPATCH_REQUEST` for the method-implementer (or numerics-investigator
  if no method-implementer is available) to fix the SymPy.

### For test-architects

- **NEVER design an L1 cross-check that uses primitives shared
  between Branch 1 and Branch 2.** Above the trusted-library line,
  structural independence is required.
- **NEVER claim L1 verification on the basis of MMS alone.** MMS is
  a flux-shape / convergence-order pillar. Eigenvalue claims need
  closed-form or semi-analytical references.

### For all agents

- **NEVER let "two SymPy derivations agree" satisfy the
  cross-check.** ERR-032 proves that two derivations sharing an
  upstream identity will agree to machine precision while both being
  wrong. Cross-check requires structural independence, which usually
  means Branch 2 (numpy) cross-checking Branch 1 (SymPy or
  semi-analytical).

---

## Worked example — the Plan 2 Variant α pipeline

The Variant α Green's function work in this project followed this
discipline (with a notable A1+A2 deviation that the discipline would
have caught earlier):

**Branch 1 SymPy** (algebra-of-record):

- Path: `orpheus/derivations/continuous/peierls/origins/specular/greens_function.py`
- States used: 1A (closed-form for V_α1, V_α2, V_α3 algebraic identities).
- Verification functions:
  `derive_operator_constant_trial_closed_sphere()`,
  `derive_T00_equals_P_ss_sphere()`,
  `derive_alpha_zero_kernel_reduction()`.
- Test gate: `tests/derivations/test_peierls_greens_function_symbolic.py`
  with one foundation-tagged test per `derive_*`.

**Branch 1 semi-analytical** (PS-1982 reference, State 1B):

- Path: `orpheus/derivations/continuous/peierls/ps1982_reference.py`
- States used: 1B (SymPy-derived kernel `[E_1(|r-x|) - E_1(r+x)]`,
  `scipy.integrate.quad` evaluates Nyström matrix with QAGS log-
  singularity treatment).
- Used as L1 reference for the vacuum-BC case.

**Branch 2 production**:

- Path: `orpheus/derivations/continuous/peierls/greens_function.py`
- Uses numpy + scipy.special.exp1 + scipy.interpolate.CubicSpline +
  `numpy.polynomial.legendre.leggauss`.
- Trajectory + bounce-period integrals along bouncing characteristics.

**L1 cross-check**:

- Closed-sphere V_α1 numerical (`test_peierls_greens_function_solver.py`):
  Branch 2 reproduces V_α1 algebraic identity to machine precision.
- PS-1982 vacuum (`test_peierls_greens_function_xverif_ps1982.py`):
  Branch 2 (vacuum branch, α=0) agrees with PS-1982 reference solver
  to ≤ 1e-4 relative across 4 configurations.
- Multi-group (`test_peierls_greens_function_mg.py`): Branch 2 (G≥1)
  agrees with `kinf_and_spectrum_homogeneous` transfer-matrix at
  closed sphere.
- Multi-region fixed-source (`test_peierls_greens_function_garcia2021.py`):
  Branch 2 (multi-region) agrees with Garcia 2021 Table 5 ppP_N
  benchmark.

**Sphinx**:

- Stub at `docs/theory/peierls_greens.rst` (extended into rich
  narrative by archivist after each phase).

**Discipline deviation that would have been caught**:

- During A1 (vacuum BC implementation), the trajectory direction
  (forward vs backward) was wrong in Branch 2 but Branch 1 V_α1
  algebraically cancelled the error (`L_first` cancels in the
  closure). The bug was caught only by the PS-1982 cross-check
  (Branch 1 State 1B), not by V_α1 (Branch 1 State 1A).
- This is the canonical case for **why structural-independence
  between Branch 1 states matters**: the closed-form-only V_α1 was
  blind to L_first by construction; the semi-analytical PS-1982
  was sensitive to it. **Both states should be implemented when the
  problem admits both.**

The discipline inherits this lesson: **for problems where 1A and
1B are both realizable, implement both as separate Branch-1 forms.
1A catches algebraic errors; 1B catches discretization errors that
1A's algebraic cancellation might mask.**

---

## Pointers

- **V&V principles**: `.claude/skills/vv-principles/SKILL.md`
  — three pillars, structural independence, two-step semi-analytical
  ladder, MMS operational rules.
- **Sub-agent orchestration**: `.claude/skills/subagent-handoff-protocol/SKILL.md`
  — how method-implementer dispatches archivist for the rich
  narrative, how augmentation lets the main agent harvest related
  context during dispatch.
- **Numerical bug signatures**:
  `.claude/skills/numerical-bug-signatures/SKILL.md` — the
  recognition catalogue for "Branch 1 / Branch 2 disagreement"
  fingerprints (constant-factor convention drift, forward-vs-backward
  trajectory, etc.).
- **ERR catalogue**: `.claude/skills/vv-principles/error_catalog.md`
  — every L0-caught bug, including ERR-032 (the canonical
  shared-upstream-identity instance).
- **Plan 2 / Variant α work** (worked example):
  `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
  — running closeout with the full A1-A3 + Plan-(b) pipeline as
  a case study.
