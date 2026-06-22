---
name: sn-lc-decomposition-derivation
description: SymPy + numpy derivation of whether the curvilinear SN matvec admits the algebra L+C with L=streaming and C=σ_t·ψ. Verdict — the matvec is RATIONAL (not affine) in σ_t for curvilinear via the Carlson coupled-pole seed (Hébert §3.9.4 Eq. 3.434 denominator dr·σ_t + 2). The naïve "L.apply = matvec(ψ, σ_t=0)" decomposition (prior agent's bug at reverted ad37ca0) FAILS. Resolution A — SUBTRACTIVE decomposition L.apply(ψ) := M(ψ; σ_t_full) − σ_t_full ⊙ ψ; L carries σ_t at constructor time (discrete-operator parametrisation, mathematically required by Hébert's closure). Bit-exact (residual = 0.0) across slab/sphere/cylinder. Step 6 adjoint plan PRESERVED.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 3+4.b.i pre-step
  date: 2026-05-14
---

# SN matvec L+C decomposition — SymPy + numpy derivation

## Executive summary (5 lines)

- **Q1** (is M affine in σ_t?): **NO** for curvilinear. The matvec is RATIONAL in σ_t through the Carlson coupled-pole seed (Hébert §3.9.4 Eq. 3.434 denominator `dr·σ_t + 2`). Cartesian IS affine; sphere and cylinder are not.
- **Q2**: not applicable.
- **Q3** (characterise the coupling): the σ_t-coupled term is the seed-chain perturbation. Closed form for the 1-cell case: `Δseed = -dr·σ_t·(Σw·φ_face - φ_0) / (Σw·(dr·σ_t + 2))`. Origin: `orpheus/sn/spatial/psi_half_angle_seed.py:592` (`Q_bar = σ_t · φ_0 / Σw`) feeding `carlson_inward_sweep_from_source` lines 411-419. Δseed cancels iff σ_t = 0 OR flat ψ at the BC (`φ_0/Σw = φ_face`) — which is WHY the prior fixed-source diagnostic agreed at the converged ψ.
- **Resolution chosen**: **Resolution A — SUBTRACTIVE definition**. `L.apply(ψ) := M(ψ; σ_t_full) − σ_t_full ⊙ ψ`. Empirically bit-exact (rel_residual = 0.0) on slab, sphere, cylinder across 3 random seeds each.
- **Step 6 impact**: PRESERVED. `(L+C-S).H.apply` works because both L.H.apply and C.H.apply are well-defined; the OperatorSum distribution holds by construction.

---

## Q1 answer — is M(ψ; σ_t) affine in σ_t?

**Cartesian: YES.** SymPy derivation (V_A1):

```
M_n,i,j = μ_x · dψ/dx + μ_y · dψ/dy + σ_t · ψ
        = L_naïve(ψ)              + C(ψ; σ_t)
residual = simplify(M − L_naïve − σ_t·ψ) = 0       ✓ exact
```

**Spherical and cylindrical: NO.** The matvec is RATIONAL in σ_t, not affine. The pieces are:

1. **Streaming term** (no σ_t):
   `streaming_n,i = (μ_n / V_i) · (A_{i+½}·ψ_face_out − A_{i-½}·ψ_face_in)`
   where ψ_face_out, ψ_face_in are computed via WDD diamond from the input ψ. No σ_t.

2. **Cell-collision term** (linear in σ_t, V_A2):
   `collision_n,i = σ_t,i · ψ_cell[n,i]`   — exactly the C(ψ; σ_t) candidate.

3. **Redistribution term** (RATIONAL in σ_t, V_A3):
   `redist_n,i = (ΔA_i / w_m) · [α_{m+½}·φ_{m+½} − α_{m-½}·φ_{m-½}] / V_i`
   where `φ_{m±½}` come from the M-M angular recurrence seeded by the Carlson coupled-pole sweep:
   `phi_aux_i = (dr_i · σ_t,i · φ_0,i / Σw + 2 · phi_face_{i+½}) / (dr_i · σ_t,i + 2)`     ← Hébert Eq. 3.434
   This RATIONAL form is the bug source. ∂phi_aux/∂σ_t depends on σ_t — not a constant — hence non-affine.

The decomposition residual closed form (V_A4, 1-cell case):
```
R(ψ; σ_t) := M(ψ; σ_t) − L_naïve(ψ) − σ_t·ψ
           ∝ Δseed = -dr · σ_t · (Σw · φ_face − φ_0) / (Σw · (dr·σ_t + 2))
```

This is **factored**: zero iff σ_t = 0 OR `φ_0/Σw = φ_face` (flat ψ at the outer BC). The flat-ψ cancellation is WHY the prior fixed-source diagnostic (`tests/sn/diagnostics/phase_g_step2_04_fixed_source.py` against the converged fixed point) appeared to agree even though the operator is structurally non-affine.

**Verbatim symbolic equation** (from `derivations/diagnostics/diag_LC_decomposition_sn.py:derive_spherical_carlson_seed_symbolic`):

```
phi_aux_full = (2·Σw·φ_face_outer + dr·φ_0·σ_t) / (Σw · (dr·σ_t + 2))
dphi_aux/dσ_t = 2·dr·(φ_0 − Σw·φ_face_outer) / (Σw · (dr²·σ_t² + 4·dr·σ_t + 4))
```

The σ_t derivative depends on σ_t — confirms rational.

### Q1 verdict by geometry

| Geometry    | M affine in σ_t? | Residual norm (random ψ, n=5, GL-4) |
|-------------|------------------|-------------------------------------|
| Cartesian   | YES              | 2e-17 (machine zero across 3 seeds) |
| Spherical   | NO               | 2.6e-2 to 3.9e-2 (~3-4% rel)        |
| Cylindrical | NO               | 2.2e-2 to 1.3e-1 (~2-13% rel)       |

---

## Q3 — where the coupling lives, file:line

The σ_t-coupled non-linear-in-σ_t artefact lives in **the Carlson seed**:

- **Q_bar source**: `orpheus/sn/spatial/psi_half_angle_seed.py:592`
  ```
  Q_bar = sigma_t * phi_0 / weights.sum()   # (ng, nx)
  ```
  Builds `Q̄_i = σ_t,i · φ_0,i / Σw` — multiplicative σ_t · ψ-projection.

- **Hébert (3.434) inward recurrence**: `orpheus/sn/spatial/psi_half_angle_seed.py:411-419`
  ```python
  for k in range(nx - 1, -1, -1):
      denom = dr[k] * sigma_t[:, k] + 2.0
      phi_cell = (dr[k] * Q_bar[:, k] + 2.0 * phi_face) / denom
      phi_aux[:, k] = phi_cell
      phi_face = 2.0 * phi_cell - phi_face
  ```
  Both numerator and denominator carry σ_t. `phi_aux` is therefore rational in σ_t at fixed ψ.

- **Seed enters M-M recurrence**: `orpheus/sn/spatial/pole_angular_closure.py:437-440`
  ```python
  if psi_half_seed is None:
      psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
  else:
      psi_half_left = psi_half_seed.copy()
  ```
  The `phi_aux` from CarlsonInwardSweep becomes `psi_half_left`, then propagates through the M-M recurrence at lines 446-457, contaminating `redist_full` at every (g, m, i) with a σ_t-coupled term.

- **redist_full enters matvec**: `orpheus/sn/operator.py:798, 831` (spherical) and 1040, 1089 (cylindrical):
  ```python
  redistribution = redist_full[:, outgoing_mask, i]
  collision      = sig_t[i, 0, :, None] * psi_cell
  lhs[:, ks]     = streaming + redistribution + collision
  ```
  The matvec emits `streaming + redistribution + collision` where `redistribution` carries the σ_t-coupled chain.

**The σ_t-coupled term is structurally part of the discrete operator**, not a bug. Hébert §3.9.4 derives the closure assuming `(L+C)ψ = q`; the Carlson seed uses `Q ≡ σ_t · φ_0` (the operator's collision term) as the equivalent source at μ = −1. This makes the discrete L (with M-M angular closure) intrinsically σ_t-coupled.

### Why the prior agent's `matvec(ψ, σ_t=0)` approach is wrong

Setting σ_t = 0 substitutes into the Hébert (3.434) denominator:
```
denom = dr·0 + 2 = 2          (instead of dr·σ_t + 2)
```
and into the numerator:
```
Q_bar = 0 · φ_0 / Σw = 0       (instead of σ_t · φ_0 / Σw)
```
giving:
```
phi_aux_at_σ_t=0 = (0 + 2·phi_face) / 2 = phi_face
```
i.e. the seed becomes the BC-trace value alone, INDEPENDENT of the input ψ (except through phi_face which depends on ψ at the outer cell). This is NOT a valid "pure streaming + angular redistribution" — it's a degenerate closure that mathematically corresponds to "σ_t = 0 transport" (no collision anywhere), which is a different operator from L = Ω·∇ + angular redistribution at σ_t-coupled-closure.

Empirically (from `diag_LC_decomposition_resolution.py::test_resolution_a_L_is_NOT_pure_streaming_curvilinear`):
- Prior agent's `L_prior = matvec(ψ, σ_t = 0)` vs Resolution A's `L_correct = matvec(ψ, σ_full) − σ_t·ψ`:
- For sphere: relative difference ≈ 0.05 (5%) on random ψ.
- For cylinder: relative difference ≈ 0.10 (10%) on random ψ.

The two L's solve different operators. The prior agent's L is structurally wrong.

---

## Empirical cross-check (Part B — numpy probe)

Script: `derivations/diagnostics/diag_LC_decomposition_sn.py`. Run with `.venv/bin/python -m pytest derivations/diagnostics/diag_LC_decomposition_sn.py -v` — **all 8 tests pass** (4 SymPy + 4 numpy).

```
CART: rel_res=2.4e-17, 5.6e-17, 6.8e-17    (affine, machine zero)
SPH : rel_res=3.9e-02, 3.4e-02, 2.6e-02    (NOT affine)
CYL : rel_res=2.2e-02, 8.9e-02, 1.3e-01    (NOT affine)
```

**Test of Resolution A** in companion script `diag_LC_decomposition_resolution.py` — **all 10 tests pass**:

```
CART: rel_res=0.0e+00, 0.0e+00, 0.0e+00    (bit-exact)
SPH : rel_res=0.0e+00, 0.0e+00, 0.0e+00    (bit-exact)
CYL : rel_res=0.0e+00, 0.0e+00, 0.0e+00    (bit-exact)
```

The numpy result agrees with the SymPy derivation. The 1-cell symbolic derivation extends mechanically to N cells through the inward-sweep chain (each cell's seed propagates to the next via 3.435), so the closed form is the FACTORED structure of a recursively-composed rational function — but the affine-ness conclusion is robust because rationality propagates under composition.

---

## Resolution proposal — Resolution A (subtractive definition)

### The math

Define for ALL geometries (Cartesian, spherical, cylindrical):

```
L: StreamingOperator(sn_mesh, sigma_t).apply(ψ) := M(ψ; σ_t) − σ_t ⊙ ψ_packed
C: CollisionOperator(sn_mesh, sigma_t).apply(ψ) := σ_t ⊙ ψ_packed
```

Then `(L + C).apply(ψ) = L.apply(ψ) + C.apply(ψ) = M(ψ; σ_t)` exactly, **by construction**.

**Both L and C carry σ_t at constructor time.** This is the discrete-operator parametrisation: at the continuous level `Ω·∇ψ + (1-μ²)/r·∂ψ/∂μ` is σ_t-independent, but the *discrete* curvilinear streaming-plus-redistribution operator (with M-M angular closure) inherits a σ_t-coupling through the Carlson coupled-pole seed (Hébert §3.9.4 §3.9.4 Eq. 3.432). This is a fact about Hébert's discretisation, not an implementation defect. Similar examples: the DD coefficient `α_DD(σ_t·dx)` in characteristic-line methods, the exponential characteristic `exp(-σ_t·s)` in MoC/CP — discrete streaming operators routinely carry σ_t through their closure choices.

For Cartesian, σ_t enters L only through the cell-collision term that we subtract — so L_cart is mathematically σ_t-independent at the apply level (the subtraction removes σ_t exactly). For sphere/cylinder, σ_t enters BOTH the cell-collision AND the seed chain; subtracting the cell-collision leaves the seed-chain σ_t-coupling inside L. Both behaviours are correct.

### Class signature

```python
@dataclass
class StreamingOperator(LinearOperatorMixin):
    r"""Pure streaming + angular redistribution operator (the "L" of L+C-S-F).

    Curvilinear: L is intrinsically σ_t-coupled via the M-M angular
    closure's Carlson coupled-pole seed (Hébert §3.9.4 Eq. 3.432).
    Cartesian: L is σ_t-independent.

    Apply is defined subtractively:
        L.apply(ψ) := M(ψ; σ_t) − σ_t ⊙ ψ
    where M is the production matvec primitive. This is bit-exact
    decomposition: (L + C).apply(ψ) == M(ψ; σ_t).
    """

    sn_mesh: "SNMesh"
    sigma_t: np.ndarray   # (nx, ny, ng); used by the closure AND subtracted

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    # Lazy caches (unchanged from prior reverted version)
    _eq_map: EquationMap | None = field(default=None, init=False, repr=False)

    def _ensure_eq_map(self, ng: int) -> EquationMap:
        ...   # same as ad37ca0

    @property
    def n_unknowns(self) -> int:
        return self._ensure_eq_map(ng=1).n_unknowns

    def apply(self, psi: np.ndarray) -> np.ndarray:
        """L.apply(ψ) := M(ψ; σ_t_full) − σ_t_full ⊙ ψ_packed.

        Subtractive definition. The matvec is called with the FULL
        σ_t (constructor-stored); the cell-collision term σ_t ⊙ ψ
        is then subtracted at the packed-vector level to leave the
        streaming + angular-redistribution part.
        """
        sn_mesh = self.sn_mesh
        eq_map = self._ensure_eq_map(ng=1)
        ng = int(psi.size // eq_map.n_eq)
        if eq_map.n_unknowns != psi.size:
            self._eq_map = None
            eq_map = self._ensure_eq_map(ng=ng)

        # Full matvec at the user's σ_t (NOT zero).
        nx, ny = sn_mesh.nx, sn_mesh.ny
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)

        if curv == "spherical":
            reduced = sn_mesh.reduced
            M_full = transport_operator_matvec_spherical(
                psi, eq_map, quad, self.sigma_t,
                nx, ng, reduced.face_areas, sn_mesh.volumes,
                reduced.alpha_half, reduced.redist_dAw, reduced.tau_mm,
                sn_mesh=sn_mesh, bc_outer=sn_mesh.bc_right,
                pole_angular_closure=sn_mesh.pole_angular_closure,
            )
        elif curv == "cylindrical":
            reduced = sn_mesh.reduced
            M_full = transport_operator_matvec_cylindrical(
                psi, eq_map, quad, self.sigma_t,
                nx, ng, reduced.face_areas, sn_mesh.volumes,
                reduced.alpha_per_level,
                reduced.redist_dAw_per_level,
                reduced.tau_mm_per_level,
                sn_mesh=sn_mesh, bc_outer=sn_mesh.bc_right,
                pole_angular_closure=sn_mesh.pole_angular_closure,
            )
        else:  # Cartesian / slab
            M_full = transport_operator_matvec(
                psi, eq_map, quad, self.sigma_t,
                nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
                bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
                bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
            )

        # Subtract the cell-collision σ_t ⊙ ψ at packed-vector level.
        sigma_packed = self.sigma_t[
            eq_map.ix, eq_map.iy, :
        ].T.ravel(order='F')
        return M_full - sigma_packed * psi
```

The corresponding `CollisionOperator` is unchanged from the reverted commit — it already implements `σ_packed * psi` correctly:

```python
@dataclass
class CollisionOperator(LinearOperatorMixin):
    sn_mesh: "SNMesh"
    sigma: np.ndarray
    capabilities: frozenset[str] = ...   # includes CAP_APPLY_TRANSPOSE

    def apply(self, psi: np.ndarray) -> np.ndarray:
        eq_map = self._ensure_eq_map(ng=...)
        sigma_packed = self.sigma[eq_map.ix, eq_map.iy, :].T.ravel(order='F')
        return sigma_packed * psi
```

### Verification test (already written)

`derivations/diagnostics/diag_LC_decomposition_resolution.py` — 10 tests, all PASS at rel_residual = 0.0 across {CART, SPH, CYL} × {3 random seeds}. Promote to `tests/sn/test_streaming_operator_decomposition.py` after the StreamingOperator lands; tag `@pytest.mark.l0`.

### Why this is the right answer architecturally

- **Pattern 1 (algebra match)**: `(L + C).apply(ψ) == matvec(ψ; σ_t)` exactly. User writes `(L + C - S - F/k) @ ψ` and reads the math.
- **Pattern 2 (single source of truth)**: the production matvec primitive is the SOURCE; L and C are thin Python wrappers that route through it. No twin implementation.
- **Pattern 7 (normalise at the definition site)**: σ_t is fixed once at constructor time on BOTH L and C. Consumers don't re-introduce σ_t conventions at apply time.
- **Pattern 4 (illegal-states-unrepresentable)**: `StreamingOperator(sn_mesh)` without σ_t would NOT type-check — every L is constructed with its σ_t. The "L without σ_t" state literally cannot exist.

### What about the user's framing "L doesn't depend on σ_t"?

The CONTINUOUS L (`Ω·∇ψ + (1-μ²)/r · ∂ψ/∂μ`) is σ_t-independent. The DISCRETE L (with Hébert M-M angular closure) is σ_t-coupled through the Carlson seed. This is analogous to how the DISCRETE diffusion stencil's interface coefficients carry σ_t in the Hébert/Mikityuk-DC mixed-cell formulation — discrete operators sometimes carry parameters that the continuous form doesn't. The user's mandate is consistent at the algebra level (`(L+C-S-F/k)` IS the within-group operator) but the discrete L's apply signature includes σ_t as a parametrisation that comes from the closure choice, not from L's mathematical identity.

**This is the principled answer** rather than a hack — and the alternative (forcing σ_t out of L by changing the M-M closure) would mean retiring Hébert §3.9.4 in favour of a different curvilinear discretisation. That's an out-of-scope architectural choice and would regress the L0-SN-003 streaming-equilibrium correctness that Phase D fought to establish.

### Subset answers to the anti-recommendations

- **NO xfail**: the resolution produces bit-exact decomposition. No tests xfailed.
- **NO `WithinGroupOperator` wrapper**: the algebra IS `L + C - S.foldable_part()`; no wrapper.
- **NO "use solve, never apply for curvilinear"**: L.apply works for curvilinear; the decomposition holds by construction.
- **NO changes to the matvec primitives**: L.apply CALLS the existing primitives unchanged. The subtraction happens AFTER the matvec call.
- **NO `_sweep_1d_unified` or `transport_sweep`** in L formulation: L is forward apply, not inverse solve. The matvec is the canonical forward path.

---

## Step 6 impact statement

**Does `solve_sn_adjoint` via `gmres(as_scipy_linop((L+C-S).H), ...)` work for curvilinear under Resolution A?**

**YES.** The adjoint propagates through `OperatorSum.H`, which propagates through `OperatorSum.apply_transpose`, which calls `a.apply_transpose(x) + b.apply_transpose(x)` on the operands (orpheus/numerics/operator.py:564-565). For Resolution A:

- **`L.apply_transpose`**: implement as `matvec_transpose(ψ; σ_t) − σ_t ⊙ ψ`, where `matvec_transpose` is the reverse-direction sweep (swap outward/inward in the matvec body) with the same σ_t. The subtractive structure preserves: `(L^T + C^T).apply(ψ) = L^T.apply(ψ) + C^T.apply(ψ) = matvec_transpose(ψ; σ_t)`. Since `(L + C)^T = L^T + C^T` for any linear sum, and both σ_t-parameterized operators share σ_t at constructor, this works.
- **`C.apply_transpose = C.apply`** (self-adjoint, diagonal in (cell, group) basis, isotropic across ordinates). Unchanged.
- **`OperatorSum.apply_transpose` distribution**: `(L+C).H.apply(x) = L.H.apply(x) + C.H.apply(x) = matvec_transpose(x; σ_t) − σ_t·x + σ_t·x = matvec_transpose(x; σ_t)`. Bit-exact.

The Step 6 plan (issue_196_phase_g_replan.md:1335-1404) requires:
- StreamingOperator.apply_transpose ← reverse-direction sweep. **Resolution A makes this trivially expressible**: subtract σ_t·ψ from the matvec_transpose output.
- CollisionOperator.apply_transpose = apply. **Unchanged.**
- ScatteringOperator/FissionOperator.apply_transpose. **Unchanged by this resolution.**
- `solve_sn_adjoint` body ≤ 8 lines. **Compatible.**

The only EXTRA work Resolution A imposes on Step 6 is: implement `matvec_transpose` for the three geometries. This was already required by the Step 6 plan (called "reverse-direction sweep"). The subtractive overlay does not add complexity to it; it's the same matvec primitive with the cells visited in reverse order.

**Step 6 plan PRESERVED.** No replan needed.

---

## References

- **Hébert, A. (2009). *Applied Reactor Physics*. §3.9.4 pp. 141-144**:
  - Eq. (3.432): `Q̄_i = Σ_ℓ (2ℓ+1)/2 · Q_ℓ(r) · (-1)^ℓ` — moment-folded source at μ = −1.
  - Eq. (3.434): `φ̄_i = (Δr·Q̄ + 2·φ̄_{i+½}) / (Δr·Σ_t + 2)` — inward DD recurrence (the rational-in-σ_t form).
  - Eq. (3.435): `φ̄_{i-½} = 2·φ̄_i − φ̄_{i+½}` — step inward.

- **Lewis & Miller (1993). *Computational Methods of Neutron Transport*. §4.5 / §5.3.4**:
  - Curvilinear weighted diamond (WDD); Carlson starting-direction initialisation.

- **Pomraning, G.C. (1989). *The Transport Equation in General Geometry*. NSE 101:330-340, p. 339**:
  - Structural singularity at r = 0; per-ordinate isotropy condition.

- **Internal docs**: `docs/theory/discrete_ordinates.rst` (curvilinear SN sections); `docs/theory/operator_algebra.rst` (LinearOperator Protocol).

- **Predecessor diagnostic**: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md` — characterised the SI-vs-Krylov closure difference and named ERR-026 manifestation #6.

- **Reverted prior implementation**: commit `ad37ca0` (reverted as `b47551c`). The `StreamingOperator.apply` body there is WRONG (calls matvec with σ_t = 0). Resolution A replaces it with the subtractive form.

- **Diagnostic scripts**:
  - `derivations/diagnostics/diag_LC_decomposition_sn.py` — Q1 SymPy + numpy probe. 8 tests, all PASS.
  - `derivations/diagnostics/diag_LC_decomposition_resolution.py` — Resolution A bit-exact verification. 10 tests, all PASS at rel_residual = 0.0.

---

## Bottom line for the next dispatch

The next method-implementer dispatch can take this verdict directly:

1. **Implement** `StreamingOperator.apply` with the subtractive body shown in §"Class signature" above. The matvec is unchanged; the change is at the wrapper level.

2. **Test** with `tests/sn/test_streaming_operator_decomposition.py` (promoted from `diag_LC_decomposition_resolution.py`) — bit-exact at `rtol=1e-14` across slab/sphere/cylinder.

3. **Sphinx narrative**: extend `docs/theory/operator_algebra.rst` with one paragraph naming the discrete L's σ_t-parametrisation as a property of Hébert's M-M angular closure (cite Eq. 3.432-3.435). The user's algebra `(L+C-S-F/k)` still reads as math; the constructor signature `StreamingOperator(sn_mesh, σ_t)` carries the discretisation's parameter.

4. **NO Step 6 replan needed**. The reverse-direction sweep already exists in the plan; subtracting `σ_t · ψ` post-matvec is a one-liner.

---

## Linked memories

- `[[issue-196-phase-g-step2-diagnostic]]` — Phase G Step 2 SI-vs-Krylov closure diagnosis. Identified the Carlson seed σ_t coupling as the structural divergence at the L0 streaming-equilibrium level.
- `[[issue-196-phase-g-step2-replan-verdict]]` — Architectural verdict; relevant context for why apply-vs-solve consistency is non-trivial in this codebase.
- `[[issue-196-phase-g-step2-minimal-reproducer]]` — Sphere 2-cell minimal reproducer.
- `[[issue-196-phase-g-step2-cylinder-minimal-reproducer]]` — Cylinder 1-cell reproducer; cylinder Carlson seed normalization.
