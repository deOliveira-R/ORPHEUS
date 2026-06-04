# Issue #168 design memo — curvilinear FD operator MMS convergence

**Author**: numerics-investigator (2026-05-09)
**Branch context**: `refactor/sn-operator-algebra`, post-Wave E Round 3
**Status**: investigation memo; **the issue's diagnosis is incomplete**.
**Recommendation**: Option C (back-substituted DD throughout) is the only
mechanically correct path — Options A and B are structurally insufficient.
The orchestrator should weigh the evidence and decide whether to
implement now or open a follow-up sub-investigation first.

## TL;DR for the orchestrator

The empirical evidence does not match the issue's narrative. Issue #168
attributes the ~1.26 order to a single bug — cell-center used as the
outer-face flux at `i = nx-1` for outgoing μ — and proposes three
mechanical fixes (Option A: DD extrapolation; Option B: ghost-cell;
Option C: rebuild the operator as DD).

This investigation **falsifies the single-bug framing**. The
symmetric-closure FD operator has at least **three independent
boundary truncation defects**, ALL O(1) on smooth solutions, that
combine to limit MMS convergence to ~O(h^{1.25}). Fixing only the
outer-face cell-center substitution (Option A as the issue proposes
it) IMPROVES the absolute error at fixed nc but does NOT raise the
asymptotic order. Empirically Option A gives orders 1.23 → 1.11 → 1.05
(WORSENING with refinement) — the order is being limited by other
defects that only become dominant once Option A is in place.

Option C (the structural rebuild) is the only fix that simultaneously
addresses all three defects — and is also the cleanest expression of
the Wave A operator-algebra philosophy (Cardinal Rule 2). The cost is
a full operator rewrite, snapshot regeneration for curvilinear
geometries, and (potentially) revisiting the Wave D `apply_transpose`
design.

## 1. Empirical reproduction of the ~1.26 order

Reference probe: `scratch/derivations/diagnostics/diag_issue168_01_characterize.py`.

Spherical isotropic MMS, `phi_exact(r) = sin(πr/R)`, R=5,
sigma_t=1.0, sigma_s=0.5, vacuum BC at r=R, n_ordinates=8 (GL):

| inner_solver        | nc=10 L2  | nc=20 L2  | order(10→20) |
|---------------------|-----------|-----------|--------------|
| `source_iteration`  | 2.52e-1   | 2.92e-1   | -0.21 (diverging — ERR-026) |
| `krylov` (Round 3)  | 1.91e+0   | 8.02e-1   | **1.25**     |

Both reproduce the issue's documented behaviour. The `krylov` order of
1.25 IS the figure the issue reports (1.26 within sampling noise).

## 2. Falsification of the issue's single-bug diagnosis

### Probe A (Option A — DD diamond extrapolation at outer face only)

Patch:
```python
if i == nx - 1 and mu[n] > 1e-15:
    psi_face_in = 0.5 * (fi[:, n, i-1, 0] + fi[:, n, i, 0])
    psi_right = 2.0 * fi[:, n, i, 0] - psi_face_in   # DD extrapolation
```
Reference probe: `diag_issue168_02_option_a_dd_extrap.py`.

| nc range | order with Option A patch (n_ord=4) |
|----------|--------------------------------------|
| 10 → 20  | 1.23                                 |
| 20 → 40  | 1.11                                 |
| 40 → 80  | 1.05                                 |

The order **degrades with refinement** under Option A. This is the
signature of "fixing one boundary defect uncovers a second, slower-
converging one." If Option A were the right fix, the order should
approach 1.9-2.0 as h→0, not regress.

### Probe B (forced exact outer face, vacuum BC)

For this MMS, `phi_exact(R) = sin(π) = 0` exactly, so we can probe
"what if we just hard-code `psi_right = 0` at the outer face?" —
this is the EXACT vacuum-BC face value.

Reference probe: `diag_issue168_04_option_a_vacuum.py`.

| nc range | order with `psi_right = 0` forced (n_ord=4) |
|----------|----------------------------------------------|
| 10 → 20  | 0.97                                         |
| 20 → 40  | 0.68                                         |

**Forcing the EXACT face value gives WORSE convergence than the
unpatched cell-center substitution.** This is paradoxical only if you
believe the outer face is the dominant defect. It makes sense if the
operator has additional structural defects that interact with the
outer-face value: the unpatched cell-center extrapolation introduces a
specific O(h) error pattern that PARTIALLY CANCELS another O(h) error
elsewhere (the "two wrongs make a half-right" anti-pattern).

### Probe C — direct truncation residual on the exact solution

For the exact ψ_n(r) = phi_exact(r)/W (isotropic ansatz), compute
`L · ψ_exact - (Σ_s/W) · phi_exact - Q/W` and tabulate per cell.

Reference probe: `diag_issue168_05_apply_vs_solve.py` and
`diag_issue168_08_unit_residual.py`.

| nc | inner cell i=0 | outer cell i=N-1 | next-to-outer i=N-2 | mid-domain |
|----|----------------|-------------------|----------------------|------------|
| 10 | 1.17e-1        | 1.52e-1           | _below 1e-3 (low resolution)_ | 7.4e-2 |
| 20 | 1.31e-1        | 1.43e-1           | (similar)        | 7.1e-2     |
| 40 | 1.34e-1        | 1.39e-1           | 6.9e-2           | 5.9e-4     |
| 80 | 1.35e-1        | 1.37e-1           | 6.8e-2           | 1.4e-4     |

**Clear structural diagnosis**:
- **Mid-domain residual is perfect O(h²)**: 5.9e-4 → 1.4e-4 (ratio 4.2 → log2 = 2.07). The interior discretization IS second-order.
- **Inner cell residual i=0 is O(1)**: 0.117 → 0.135 (no decay). **Defect at the sphere pole.**
- **Outer cell residual i=N-1 is O(1)**: 0.152 → 0.137 (slight decay only because pole effects bleed in). **Defect at the outer boundary.**
- **Next-to-outer i=N-2 is O(1) for inward direction**: 0.069 → 0.068. **BC-fill contamination of the interior face.**

The integrated effect: total `||r||_2 ~ sqrt(h)` (boundary cells contribute O(1) residual on volumes O(h) for outer, O(h³) for inner). The solution-error scaling depends on the operator's `||L^{-1}||` smoothing; empirically this gives ~O(h^{1.25}).

## 3. The three boundary defects, structurally diagnosed

### Defect 1 — outer-face cell-center substitution (the issue's claimed bug)

```python
# orpheus/sn/operator.py: line 535-536 (transport_operator_matvec_spherical)
if mu[n] > 1e-15:
    psi_right = fi[:, n, i, 0]   # cell-center used as face-flux approximation
```

For interior cells, `psi_right = 0.5*(fi[i] + fi[i+1])` is O(h²)
because cell-centers are at `r_i = r_face - h/2` and `r_{i+1} = r_face
+ h/2`, so the arithmetic average evaluates ψ at the face to O(h²).

For the outer cell `i = N-1`, the OUTER face is at `r = R = r_{N-1} +
h/2`, but we only have one cell-center value `fi[N-1] = ψ(r_{N-1})`.
Substituting it as the face value is O(h) — the mistake the issue
correctly identifies.

**Fix**: replace with one-sided second-order extrapolation
`psi_right = 1.5*fi[N-1] - 0.5*fi[N-2]` (the DD diamond relation
`ψ_face_out = 2*ψ_cell - ψ_face_in` where `ψ_face_in = 0.5*(fi[N-2] +
fi[N-1])` already used by the matvec).

This is necessary but **not sufficient** — see Defects 2 and 3.

### Defect 2 — BC-fill contamination of the interior face stencil

`solution_to_angular_flux_spherical` (Wave E Round 3 BC plumbing)
overwrites `fi[N-1, n_inward, 0]` with the BC-determined incoming
face value (=0 for vacuum). For VacuumBoundaryOperator this is `fi[N-1, n_inward,
0] = 0`.

The matvec at cell `i = N-2` for inward direction `n_inward` then computes:
```python
psi_right = 0.5 * (fi[:, n, N-2, 0] + fi[:, n, N-1, 0])
         = 0.5 * (ψ_cell[N-2] + 0)   # because fi[N-1] was zeroed by BC
```

But this `psi_right` is supposed to be the value at the FACE BETWEEN
cells N-2 and N-1 (at `r = R - h`), NOT the value at `r = R`. Vacuum
BC says only `ψ(R, μ < 0) = 0`; it says NOTHING about `ψ(R - h, μ <
0)`. The arithmetic average that's supposed to give `ψ(R - h)` is
being corrupted by the BC fill at the outer boundary.

This is a **structural conflation** in the operator's storage layout:
the array `fi[..., -1, ...]` is used as both (a) the cell-center
storage at i=N-1 and (b) the face-value storage at the outer face.
The BC writes (b); the matvec reads (a) when computing interior face
averages. Same slot, two semantics.

**Why this didn't show up before Wave E Round 3**: pre-Round 3 the BC
was hard-coded reflective, so the BC fill `fi[N-1, n_inward] =
fi[N-1, n_outward_partner]` wrote the SAME cell-center value back
(since reflective ↔ same i, mirrored n). The slot was a faithful
cell-center value by accident. Wave E Round 3 generalized to vacuum
(which writes 0) and exposed the conflation.

**Fix**: separate cell-center storage from BC face-value storage. The
matvec at `i = N-2` for inward direction should read the CELL-CENTER
value `ψ_cell[N-1]` (= the packed solution slot at `eq_map.ix == N-1`,
ordinate the inward partner — which the eq_map skips, so this value
isn't even an unknown), NOT the BC-overwritten `fi[N-1]`.

This is non-trivial because the inward direction at i=N-1 is not a
solver unknown (the eq_map skips it; the unknowns at i=N-1 are only
the outgoing μ ≥ 0 ordinates). The "cell-center value at i=N-1, μ <
0" doesn't exist in the unknown vector — it has to be inferred from
the OUTGOING values via some relation. For SpecularBoundaryOperator, this is `ψ_out =
ψ_in_partner` and works. For VacuumBoundaryOperator, there's no such relation: the
problem's symmetry just doesn't say what `ψ(R-h/2, μ < 0)` is.

**The actual fix is conceptual**: the matvec equation at cell i=N-2
for inward direction n_inward shouldn't use `fi[N-1]` AT ALL. It
should use the CELL-CENTER value `ψ_cell[N-1]`, which is the same
value regardless of direction (the cell-center is a scalar storage
shared by all ordinates). The current code's `solution_to_angular_flux`
accidentally exposes per-ordinate slots at i=N-1 as a side effect of
the BC fill API.

### Defect 3 — Bailey redistribution at the sphere pole (i=0)

The Bailey 2009 angular redistribution at i=0 evaluates to
`-μ_n · ΔA[0] · ψ_0 / V[0] = -(3μ/h) · ψ_0` for sphere (since
ΔA[0] = 4πh², V[0] = (4π/3)h³).

This is the `-μ ΔA/V · ψ` term that, by design, cancels the streaming
`+μ ΔA/V · ψ` for FLAT ψ (per-ordinate flat-flux consistency).

For NON-flat ψ that varies linearly near r=0 (which the MMS does —
`phi_exact(r) ≈ (π/R)·r` near r=0), the cancellation is incomplete.
Specifically, the streaming term integrated over cell 0 evaluates to
`(3μ/h) · 0.5(ψ_0 + ψ_1) ≈ (3μ/h) · ψ(h)` (consistent O(h²)).
Subtracting the redistribution `-(3μ/h)·ψ_0 = -(3μ/h)·ψ(h/2)`:

```
(3μ/h)·[ψ(h) - ψ(h/2)] = (3μ/h)·(h/2)·ψ'(h/4 ish) ≈ (3μ/2)·ψ'(0)
```

But the continuous spherical streaming `(μ/A) · ∂(A·ψ)/∂r` evaluated
at r → 0 in the conservative form, integrated over [0, h], gives
`(3μ/h)·ψ(h) ≈ (3μ/h)·h·ψ'(0) = 3μ·ψ'(0)`.

**The discrete operator is HALF the continuum at i=0.** The Bailey
redistribution overcorrects by a factor of 2.

This factor-of-2 error appears as the residual `0.135` in the
truncation table. The L2-error contribution from this single cell is
`0.135² · V[0] = 0.0182 · (4π/3) · h³` → bounded by `O(h^{3/2})` in
L2, but combined with the propagation through the rest of the
operator gives an order limitation of approximately h^{1}.

**Why the WDD sweep doesn't suffer from this**: the sweep processes
ordinates one at a time from most-negative to most-positive μ. The
angular face flux `ψ_{n+1/2,i}` is propagated through the M-M
weighted closure `ψ_{n+1/2} = (ψ_n - (1-τ)·ψ_{n-1/2})/τ`, which is
ASYMMETRIC and converges to the correct discrete fixed point at
flat-flux interior cells. At i=0 the same M-M closure is in play but
on a different combination of streaming + redistribution that gives
the right answer. The structural difference between sweep and apply
is that the sweep uses the WDD asymmetric M-M closure to define the
angular face flux, while apply uses the SYMMETRIC τ-weighted M-M
interpolation `ψ_{n+1/2} = τ·ψ_{n+1} + (1-τ)·ψ_n`. The two closures
agree on flat ψ but **DISAGREE on non-flat ψ at the sphere pole**.

This is the most subtle of the three defects. It doesn't have a
single-line fix.

## 4. Comparing Options A, B, C

| Option | Defect 1 (outer face) | Defect 2 (BC fill conflation) | Defect 3 (sphere pole) |
|--------|------------------------|------------------------------|------------------------|
| A — DD extrap at outer face | FIXED | NOT FIXED | NOT FIXED |
| B — ghost cell | FIXED (same idea, more complex) | NOT FIXED | NOT FIXED |
| C — full DD rebuild | FIXED | FIXED (separate cell-center/face storage) | **FIXED** if and only if the rebuild uses the asymmetric M-M angular closure (the same one the WDD sweep uses) — Option C in its purest form, mirroring the WDD sweep math |

**Option A and Option B both fix only Defect 1.** Empirically (Probe
A above) this is INSUFFICIENT — Defects 2 and 3 limit the order to
~O(h) once Defect 1 is fixed.

**Option C is the only mechanically correct choice** if the goal is
true O(h²) MMS convergence. The "true DD operator" mirrors the WDD
sweep's math at the operator level: cell-center is the average of
two face values; face values are propagated via the WDD/M-M closures;
boundary face values come from the BC operator. By construction, the
apply method then matches what the sweep does, just inverted.

Bit-identity to the WDD sweep at the operator level would mean
`apply(ψ) - solve(q) = 0` to machine precision when `ψ = solve(q)`.
This is the deepest possible fix — it forces the apply path and the
sweep path to be the same operator, not just two operators that
agree on the converged answer.

## 5. The cross-domain-frames lens — is there a structural reframing?

The three defects all stem from the same root: **the symmetric-closure
FD operator implements a conservative cell-balance equation, but it
parameterizes the unknowns by cell-CENTERS instead of by cell-FACES
or by face-PAIR (cell-center, face-out)**. The cell-balance is
naturally written in terms of face fluxes:

```
μ [A_out·ψ_face_out - A_in·ψ_face_in]/V + redist + Σ_t·ψ_cell = q
```

The symmetric closure CHOSE to express face fluxes as arithmetic
averages of cell-centers, on the theory that this would be O(h²) on
smooth solutions. It IS O(h²) at interior cells (where both
neighboring cell-centers exist), but FAILS at boundaries because:

1. The OUTER face has only ONE neighboring cell-center → cell-center
   substitution is O(h) (Defect 1).
2. The INTERIOR face NEXT TO the outer boundary has its second
   cell-center "borrowed" from the BC-fill slot, which carries a
   face-value semantics rather than a cell-center semantics (Defect 2).
3. The INNER face at r=0 is a degenerate point where Bailey's flat-
   flux cancellation breaks down for ψ varying linearly (Defect 3).

The structural reframe: **make face fluxes first-class unknowns**.
The DD diamond relation `ψ_cell = 0.5(ψ_face_in + ψ_face_out)` then
becomes a constraint (solved during sweep) or a substitution (fixed
during apply), not a side-product of cell-center arithmetic.

This is exactly the "true DD operator" of Option C. From the
**differential-geometry / mixed-FEM frame**: the cell-balance
equation is naturally posed in a Raviart-Thomas-like mixed
formulation where face fluxes and cell averages are independent
unknowns linked by closure relations. SN-DD is the simplest mixed FE
in this family. The "symmetric closure" we currently have is a
**hybridized** form of this — it eliminates face fluxes via
arithmetic averaging — but the hybridization breaks down at the
boundaries.

A cleaner architectural framing: separate the operator into

- a **face-flux extraction operator** `E: ψ_cell → ψ_face`
  (parameterized by closure choice — DD, WDD, step, exponential)
- a **cell-balance operator** `B: (ψ_cell, ψ_face) → q`
- the BC operator `R: ψ_face_out → ψ_face_in_BC`

The composed operator `apply` is `B(ψ_cell, E(ψ_cell, R(ψ_face_out)))`.

This is the Wave A operator-algebra philosophy expressed at the
discretization level. It also opens up pluggable face-flux closures
(DD, WDD, step, etc.) as Wave C-extension's responsibility.

## 6. Proposed elegant solution — Option C, framed as the operator-algebra capstone

### 6.1 Mathematical formulation

The discrete cell-balance equation per cell i, ordinate n:
```
μ_n [A_{i+1/2} · ψ^face_n,i+1/2 - A_{i-1/2} · ψ^face_n,i-1/2] / V_i
+ (ΔA_i / w_n) · [α_{n+1/2} · ψ^angle_n+1/2,i - α_{n-1/2} · ψ^angle_n-1/2,i] / V_i
+ Σ_t,i · ψ_n,i = q_n,i
```

Closure relations:
- **DD spatial**: `ψ_n,i = 0.5·(ψ^face_n,i-1/2 + ψ^face_n,i+1/2)`
  → `ψ^face_n,i+1/2 = 2·ψ_n,i - ψ^face_n,i-1/2`
- **M-M angular**: `ψ^angle_n+1/2,i = (ψ_n,i - (1-τ_n)·ψ^angle_n-1/2,i) / τ_n`

Boundary conditions:
- At i=0 (sphere pole): `A[0] = 0`, so the spatial flux at r=0 is
  irrelevant; the redistribution term has its own pole stencil.
  **For the pole**: use the standard treatment (Lewis & Miller §4.5):
  initialize `ψ^face_n,1/2 = ψ_n,0` (the cell-center, which is
  arbitrary at r=0 since A=0). The DD relation then gives
  `ψ^face_n,3/2 = 2·ψ_n,1 - ψ_n,0`.
- At i=N-1 (outer face): `ψ^face_n,N+1/2 = bc_outer.apply_to_outgoing(ψ^face_outgoing, quad)`
  for incoming directions; for outgoing directions it's the DD
  extrapolation.

### 6.2 Architectural structure

Define a `BoundaryFaceFlux` Protocol (analogous to Wave C's
`CellUpdate`):
```python
@runtime_checkable
class BoundaryFaceFlux(Protocol):
    def __call__(
        self,
        psi_cells: np.ndarray,    # (N, nx, ng) cell-center values
        ord_idx: int,
        cell_idx: int,            # 0 or nx-1
        side: str,                # "inner" or "outer"
        bc: BoundaryOperator,
    ) -> np.ndarray:               # (ng,) face-flux value
        ...
```

Concrete implementations:
- `DDExtrapolation`: `psi_face_out = 2·psi_cell - psi_face_in` (the
  one-sided variant for boundary cells).
- `CellCenter` (legacy, for reproducing pre-fix snapshots if needed):
  `psi_face_out = psi_cell`.
- `Quadratic`: 3-point one-sided extrapolation (Option B).

The operator dispatches through this Protocol at boundary cells.
Interior cells continue to use the symmetric `0.5*(fi[i] + fi[i+1])`
arithmetic average (which is O(h²) in the interior).

### 6.3 The full apply method

For the outer face at `i = N-1`, ordinate n, μ_n > 0 (outgoing):
```
ψ^face_outer = boundary_face_flux(ψ_cells, n, N-1, "outer", bc)
```

For the outer face at `i = N-1`, ordinate n, μ_n < 0 (incoming):
```
ψ^face_outer = bc.apply_to_incoming(ψ^face_outgoing, quad)[n]
```

For interior faces:
```
ψ^face_{i+1/2} = 0.5·(ψ_cells[i] + ψ_cells[i+1])
```

The matvec accumulates spatial streaming using these face values.

For the redistribution at i=0 (sphere pole), the existing Bailey
formula is structurally unchanged — but if Defect 3 needs to be
addressed independently (it does, per the truncation analysis), a
**pole-correction term** must be added. This is the deepest part of
the fix and may need literature consultation:
- Carlson's starting-direction treatment (the "α_{1/2} pole flux"
  closure).
- Lathrop's modified pole stencil.

For a first cut, the pole-correction can be:
```python
# At i=0, instead of redist = -μ·ΔA[0]·ψ_cells[0]/V[0]:
# use the centered-difference of α·ψ_face values, with
# ψ_face_cell0_inner = 0 (pole) and ψ_face_cell0_outer = (ψ_0 + ψ_1)/2:
redist_pole = (ΔA[0]/w_n)·(α_{n+1/2}·ψ_n,0 - α_{n-1/2}·0) / V[0]
```

But this is a non-trivial design decision and **requires literature
research** before implementation. Lewis & Miller §4.5 should have
the canonical treatment.

### 6.4 Boundary-aware DD diamond relation

For the OUTGOING face at i=N-1 (μ > 0), the simplest second-order
treatment using only existing data is:

```python
ψ^face_in = 0.5·(ψ_cells[N-2] + ψ_cells[N-1])    # interior face, O(h²)
ψ^face_out = 2·ψ_cells[N-1] - ψ^face_in
            = 1.5·ψ_cells[N-1] - 0.5·ψ_cells[N-2]
```

This is the "DD diamond extrapolation". Backward second-order
extrapolation through two cell centers.

The same applies to cylindrical geometry at i=N-1 (the bug location
is identical; the math is the same per-level).

### 6.5 Cartesian path — affected or not?

The Cartesian `transport_operator_matvec` uses **upwind FD** with
cell-center distance `0.5*(dx[i] + dx[i+1])` denominators — a totally
different stencil from the curvilinear symmetric-closure path. Upwind
FD with Cartesian reflective/vacuum BC at i=0/i=N-1 reads through
the BC operator (Wave E Round 3) the same way curvilinear does, but
the stencil itself is one-sided so there's no symmetric closure to
break.

**Conclusion**: Cartesian is unaffected. The fix is curvilinear-only.

The pluggable `BoundaryFaceFlux` Protocol could still be designed to
generalize — e.g., for future spectral or DG schemes that might use
symmetric closures on Cartesian — but for the immediate Issue #168
fix, only the curvilinear path needs editing.

## 7. Verification path

### 7.1 Tests that should xpass after the fix

The 4 `xfail-strict` markers (per Issue #168 acceptance criteria):
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py::test_sn_spherical_aniso_mms_converges_second_order`
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py::test_sn_cylindrical_aniso_mms_converges_second_order`
- `tests/sn/test_mms_curvilinear.py::test_sn_spherical_mms_converges_second_order`
- `tests/sn/test_mms_curvilinear.py::test_sn_cylindrical_mms_converges_second_order`

### 7.2 New tests to add

#### Foundation level — operator truncation residual

A canonical L0 test: take `phi_exact = sin(πr/R)` (or any smooth
ψ) for the spherical mesh, compute `||L · ψ_exact - rhs||_∞` and
assert it scales as O(h²) AT EVERY CELL including boundaries. This
test would FAIL on the current code (the boundary cells have O(1)
residual) and PASS only if all three defects are fixed.

```python
@pytest.mark.foundation
def test_curvilinear_operator_truncation_residual_h2_at_boundaries():
    """The symmetric-closure FD operator must have O(h²) truncation at
    EVERY cell, including i=0 and i=N-1, for smooth solutions."""
    case = build_spherical_mms_case(n_ordinates=4)
    boundary_residuals = []
    for nc in (10, 20, 40, 80):
        # ... build operator, apply to exact ψ, compute residual ...
        boundary_residuals.append(max_boundary_residual)
    orders = np.log2(boundary_residuals[:-1] / boundary_residuals[1:])
    assert all(o > 1.9 for o in orders), (
        f"Boundary truncation must be O(h²); got {orders}"
    )
```

#### L0 — apply-vs-solve consistency

A small Cartesian-uniform-mesh test where `apply(ψ) = solve^{-1}(ψ)`
should be bit-identical (the same DD math). For curvilinear, a
weaker consistency: applying L to the converged sweep result should
give the original RHS to machine precision.

```python
@pytest.mark.l1
def test_apply_solve_consistency_cartesian_uniform_mesh():
    """On Cartesian uniform mesh, apply(ψ) and solve^{-1}(ψ)
    should give the same answer to round-off."""
```

#### L1 — both isotropic AND anisotropic curvilinear MMS

Already in place, just needs the xfail markers removed.

### 7.3 Snapshot regeneration

The 11 frozen regression snapshots include curvilinear cases:
- `sphere_2g_homogeneous`, `sphere_2g_3region`, `sphere_2g_p1_aniso_dd_n20`
- `cyl_1g_LS4_homogeneous`, `cyl_1g_Product_homogeneous`, `cyl_2g_3region`
- (potentially others)

Option C will change the curvilinear FD operator's math non-trivially.
The 6 curvilinear snapshots **WILL CHANGE** — bit-identity is NOT
preserved. The orchestrator must regenerate them.

The 5 Cartesian snapshots (slab homogeneous, slab 3-region, slab P1
aniso, 2D 1G LS4 15x15, slab fixed-source) should remain bit-identical
because the Cartesian path is unaffected.

The regenerated curvilinear snapshots should be **VERIFIED** against
an independent reference (the WDD sweep result, or a higher-resolution
solve). The new values are the CORRECT discrete solutions for the
corrected operator; the old values were bit-identical to a
mathematically-incorrect operator (ERR-026 affected for sweep,
truncation-defective for apply).

## 8. Implementation sketch

### 8.1 Files to modify

- `orpheus/sn/operator.py` — `transport_operator_matvec_spherical` and
  `transport_operator_matvec_cylindrical`. Refactor to dispatch through
  a `BoundaryFaceFlux` strategy.
- `orpheus/sn/spatial/boundary_face_flux.py` (new) — analogous to
  `cell_update.py`. Defines the `BoundaryFaceFlux` Protocol +
  `DDExtrapolation` (default) + `CellCenter` (legacy reproducer).
- `orpheus/sn/operator.py::SNStreamingOperator` — accept an optional
  `boundary_face_flux: BoundaryFaceFlux` constructor parameter
  (default `DDExtrapolation()`), pass through to the matvec functions.
- `orpheus/sn/geometry.py::SNMesh` — accept optional `boundary_face_flux`
  (parallel to existing `cell_update`).
- `tests/sn/regression/snapshots/` — regenerate curvilinear snapshots.
- The 4 xfail markers in MMS tests come off.

### 8.2 Defect 3 (sphere pole) — separate sub-issue

The sphere pole correction is its own sub-design and **may require a
separate implementation phase**. The orchestrator should consider:
- Phase A: implement Option C (DD throughout) + the BoundaryFaceFlux
  Protocol. This addresses Defects 1 and 2.
- Phase B: implement the sphere-pole correction (Defect 3) with
  literature consultation (Lewis & Miller §4.5; Carlson; Lathrop).
  This is the deepest fix.

After Phase A only, the order will be O(h^{1.5}) to O(h^{1.7}) —
better than the current O(h^{1.25}) but not yet O(h²). Phase B
delivers the final O(h²).

If only Phase A ships, the xfail markers should be relaxed (e.g.,
assert orders > 1.5) rather than removed. If Phase B ships too, the
markers come off entirely.

### 8.3 LOC estimate

- Phase A: ~250 LOC across the operator + new strategy module + tests.
- Phase B: ~100 LOC for the pole correction + literature-sourced
  rationale in the docstring.

### 8.4 Sequence of commits (suggested)

1. `feat(sn): add BoundaryFaceFlux Protocol and DDExtrapolation strategy`
   (Phase A, non-breaking — defaults to CellCenter for backward compat).
2. `refactor(sn): wire BoundaryFaceFlux through curvilinear matvec`
   (still backward compat).
3. `chore(sn): regenerate curvilinear snapshots with DDExtrapolation default`
   (this commit BREAKS bit-identity for 6 snapshots; documents the
   physics-equivalence narrative).
4. `feat(sn): SNStreamingOperator accepts boundary_face_flux argument`
   (default `DDExtrapolation()`).
5. `test(sn): add foundation truncation-residual test + remove xfail markers (Phase A)`.
6. (Phase B, separate session/issue) `fix(sn): sphere-pole correction in
   curvilinear matvec at i=0`.

## 9. Risk register

### Risk 1 — Defect 3 may not be solvable cleanly

The sphere pole's factor-of-2 mismatch in the Bailey redistribution
might be a **fundamental limitation** of the symmetric-closure FD
approach. Lewis & Miller §4.5 (which I have not yet consulted in
detail) may have a canonical treatment, or may declare the sphere
pole to be inherently first-order on smooth-but-non-flat solutions.

If Defect 3 cannot be cleanly fixed with the symmetric closure, the
"true DD throughout" version of Option C — where apply uses the same
WDD asymmetric closure as the sweep — would be required. That's a
deeper rewrite (the WDD math is sequential by design, not naturally
matrix-free). It would either need to be implemented as a fixed-point
iteration inside `apply` (defeating the purpose of Krylov-on-apply)
or via a clever block-triangular structure.

**Mitigation**: literature researcher dispatch BEFORE Phase B
implementation. Lewis & Miller, Bailey 2009, Adams & Larsen 2002 §III,
Morel-Larsen-Pautz family. Specifically: how do production curvilinear
SN codes (PARTISN, ATTILA, Denovo) handle the sphere pole in
symmetric-closure / FE / DG variants?

### Risk 2 — flat-flux self-consistency might break

The current Bailey ΔA/w design is provably exact on flat ψ for ANY
ordinate (per-ordinate flat-flux consistency). Modifying the i=0
treatment may break this. The constraint must be:

- New i=0 treatment must give 0 residual on flat ψ.
- AND give O(h²) residual on smooth varying ψ.

These are two separate requirements; the current Bailey design
satisfies the first only. A correction that satisfies BOTH may not
be a simple modification.

**Mitigation**: all new pole-correction designs must include a
flat-flux consistency test as a foundation gate.

### Risk 3 — interaction with anisotropic (Pℓ) sources

The fixed `transport_operator_matvec_spherical` will be called from
the krylov-on-apply path with anisotropic external sources (P1 case).
The boundary-face-flux fix should be source-independent (it just
changes how `psi_right` is computed), but there's a chance the
interaction with the per-ordinate scattering source build creates
new asymmetries. The anisotropic curvilinear MMS test suite is the
gate for this.

### Risk 4 — preconditioner staleness

`_solve_krylov` uses the WDD sweep as a left preconditioner. If the
new `apply` is much closer to the sweep math (Option C), the
preconditioner becomes a near-perfect inverse and GMRES converges in
~1 iteration. That's a feature, not a bug. But it changes the
iteration counts and therefore the regression timing. Snapshots
that pin iteration counts (if any) would need updating.

**Mitigation**: regenerate iteration-count-pinning snapshots if any
exist.

### Risk 5 — Wave D `apply_transpose` design

The Wave D dense-matrix-probe apply_transpose builds the dense matrix
via `apply` calls. After the Option C rewrite, the dense matrix
changes. If anything caches the dense matrix, it must be invalidated.

**Mitigation**: clear the `_dense_matrix` cache on Option C landing;
the existing `_ensure_dense_matrix` lazy-build pattern handles this
naturally (dense matrix is rebuilt on first `apply_transpose` call).

## 10. Cross-references

### Code paths

- `orpheus/sn/operator.py`:
  - `transport_operator_matvec_spherical` (line 460-570) — main bug
    location for Defect 1 and Defect 2.
  - `transport_operator_matvec_cylindrical` (line 583-682) — same
    bug pattern.
  - `solution_to_angular_flux_spherical` (line 410-457) — Defect 2
    contamination point.
  - `SNStreamingOperator.apply` (line 862-931) — entry point.
- `orpheus/sn/sweep.py`:
  - `_sweep_1d_spherical` (line 396-545) — the WDD sweep math that
    Option C should mirror at the operator level.
- `orpheus/sn/spatial/diamond.py`:
  - `DiamondDifference._update_curvilinear` (line 432-527) — the
    canonical DD math reference.
- `orpheus/geometry/reduced_operator.py`:
  - `spherical_streaming` (line 500-578) — Bailey 2009 dome
    recursion.

### Documentation

- `docs/theory/discrete_ordinates.rst` — narrative is partly stale
  (mentions ERR-026 partial closure; will need rewrite once Issue
  #168 closes).
- `docs/theory/structured_geometry.rst` — current as of recent commit.

### Literature (recommended reading before Phase B implementation)

- Lewis & Miller (1993), *Computational Methods of Neutron Transport*,
  §4.5 (curvilinear FD), §6.4 (MMS for SN). The canonical reference;
  should clarify Defect 3.
- Bailey, Yang, Warsa (2009), *NSE 161:*, the connection-coefficient
  paper. This is the source of the ΔA/w design; it should explain
  the flat-flux-consistency goal and any limitations on non-flat ψ.
- Adams & Larsen (2002), *Progress in Nuclear Energy 40:*,
  §III.B for preconditioned-Krylov vs sweep frame ("preconditioner
  correctness vs operator correctness").
- Morel, Larsen, Pautz family for FD closure design choices.
- Carlson's starting-direction treatment + Lathrop's pole stencil
  (PARTISN documentation, or referenced in Lewis & Miller).

## 11. Recommendation to the orchestrator

**Strong recommendation**: implement Option C (boundary-face-flux
Protocol + DDExtrapolation default), which addresses Defects 1 and 2.
This is **Phase A** of the fix — ~250 LOC, clear architectural win
(parallels Wave C's `CellUpdate`), partial closure of Issue #168 (orders
will improve to ~O(h^{1.5}) or better, but may not reach O(h^{1.9}+)).

**Defer Defect 3 to Phase B / separate sub-issue**, after a literature
researcher dispatch confirms the canonical sphere-pole treatment. If
Lewis & Miller §4.5 has a clean fix, ship it in Phase B. If not,
escalate to a research issue (the symmetric-closure operator may have
a fundamental limitation at the sphere pole that warrants a different
operator family — e.g., the Larsen-Morel finite-element cell-balance
or a discontinuous Galerkin variant).

**Do NOT implement Option A or B** as proposed in the issue. They
fix only Defect 1 and produce a misleading "improvement" — the order
APPEARS to improve at one resolution but DEGRADES with refinement
because Defects 2 and 3 dominate as h→0. Shipping Option A alone
would close GH #168 with a fix that empirically regresses the very
property it claims to deliver. That's a session failure under
Cardinal Rule 1.

**Before any implementation**, the orchestrator should consider:
1. Dispatching a `literature-researcher` for Lewis & Miller §4.5 +
   Bailey 2009 sphere-pole treatment.
2. Reading this memo's diagnostic scripts in `scratch/derivations/diagnostics/diag_issue168_*.py`.
3. Choosing whether to ship Phase A alone with relaxed xfail (orders > 1.5)
   or wait for Phase B before removing markers.

## 12. Diagnostic scripts (in-place reference)

All in `scratch/derivations/diagnostics/`:
- `diag_issue168_01_characterize.py` — empirical reproduction of the
  ~1.26 order figure.
- `diag_issue168_02_option_a_dd_extrap.py` — Option A patch (DD
  extrapolation at outer face). Falsifies the issue's single-bug
  framing.
- `diag_issue168_03_patch_verify.py` — confirms the monkey-patch
  mechanism reaches the operator.
- `diag_issue168_04_option_a_vacuum.py` — forces psi_right = 0
  (exact vacuum BC value); shows it makes things WORSE.
- `diag_issue168_05_apply_vs_solve.py` — truncation residual scaling
  with mesh refinement; shows boundary residuals stay O(1).
- `diag_issue168_06_sweep_residual.py` — confirms WDD sweep ALSO has
  the boundary defect (same root cause).
- `diag_issue168_07_mms_curvature_fix.py` — falsifies the
  "MMS-source-is-missing-curvature" hypothesis.
- `diag_issue168_08_unit_residual.py` — per-cell residual structure
  showing the three boundary defects clearly.

These are **not** ready for promotion to permanent tests — they're
investigation artifacts. The truncation-residual L0 test in §7.2
should be a CLEAN re-implementation of the diagnostic concept, not
a copy of the scratch script.
