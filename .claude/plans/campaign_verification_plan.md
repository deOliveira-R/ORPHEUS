# Verification plan — operator / splitting / realization separation

**Status**: PRE-carve gate specification. Authored by `test-architect`,
2026-07-28, repo `main` @ `b0a003b4`, host venv `.venv/bin/python`,
canonical invocation `python -O -m pytest`.

**Purpose**: this document shapes the implementation plan. Every claim
below marked **MEASURED** was executed against the tree this session; the
measurement command is given so it can be re-run.

---

## 0. HEADLINE — the proposed acceptance criterion is ALREADY GREEN

> **The main agent's plan is anchored on the acceptance criterion being
> RED today. It is not. It is green, at exactly 0.0, and it cannot be
> made red because the invariance it asserts is a SIGNATURE fact, not a
> behavioural one. The anchor must be replaced before the plan is
> approved.**

Three measurements.

### 0.1 The poser takes no strategy argument (structural — cannot red)

```
build_within_group_system(sn_mesh, mat_xs, *, scattering_op, scattering_order)
```

**MEASURED** — `inspect.signature` returns
`['sn_mesh', 'mat_xs', 'scattering_op', 'scattering_order']`. Neither
`inner_schedule` nor `inner_solver` is reachable. `scattering_op` is a
cache seam and `scattering_order` is physics (Legendre truncation), not
strategy. The posed loss `A` is therefore strategy-invariant **by
construction**; a test asserting it is green-on-contact and no mutation
inside the campaign's scope can redden it.

### 0.2 `A == M − N` already holds, exactly, on both arms

Probe: `scratch/_probe_ac2.py` shape — build the record, apply `A` and
`M − ΣNᵢ` to a fixed-seed random **non-flat** state (bulk **and** trace
randomised), heterogeneous 2G, `level_symmetric` S4, mixed vacuum/reflective.

| arm | schedule | `‖Ax − (M−ΣN)x‖∞ / ‖Ax‖∞` | falsifier |
| --- | --- | --- | --- |
| 2-D Cartesian, seedless | `jacobi` | **0.000e+00** | drop `B` from `N` → 5.152e-02 |
| 2-D Cartesian, seedless | `gauss_seidel` | **0.000e+00** | drop `B` from `N` → 5.094e-02 |
| sphere, carrying (2×2 coupled) | n/a (record's own) | **3.545e-17** | sign-flip `N` → 8.569e-03 |

The falsifiers prove the probe is **non-vacuous** — it is simply green.
Note the sphere arm is already `3.5e-17`, **not** 0: the coupled grid
re-associates the block sums. Any `A == M − N` gate is therefore
`array_equal` on the seedless arm and `nulp` on the coupled arm from day
one — never a uniform bit-identity contract.

### 0.3 What IS red today — the honest red set

| # | Red fact | Evidence | Kind |
| --- | --- | --- | --- |
| **R1** | `FissionOperator.domain is None` on the production `from_solver_data` construction | **MEASURED**: `F.domain is None -> True`; every other leaf reports a concrete space (`FullFieldSpace('full_field', shape=(128,))`) | value-red |
| **R2** | `F.H` builds and silently returns a **bare Euclidean transpose** | **MEASURED**: `F.H` → `_AdjointOperator`, `.domain`/`.codomain` both `None`; `_AdjointOperator.apply` (`numerics/operator.py:1225-1229`) applies the metric only when the inner space is non-`None` | value-red, **undocumented on both facades** |
| **R3** | No pencil type exists | `git grep 'GeneralizedEigenPencil\|class .*Pencil'` → **no hits** | *unspellable* — worse than red |
| **R4** | `A − M` **RAISES** | **MEASURED**: `IncompatibleOperatorComposition: OperatorSum requires equal domains; got CoupledSpace('coupled(full_field)', shape=(128,)) and FullFieldSpace('full_field', shape=(128,))` — `A` and its own `M` live on different carriers | structural-red |
| **R5** | `A == M − N` is asserted **nowhere** | no `__post_init__` on `WithinGroupSystem` (`sn/coupled_system.py:283-330`); no test | ungated-but-true |
| **R6** | The carrier guard is **non-uniform** across leaves | **MEASURED**: `StreamingCollisionOperator.apply` raises a typed `TypeError` with a remediation message (`streaming.py:153`); `SNBoundaryOperator.apply` raises a raw `AttributeError: 'CoupledField' object has no attribute 'interior'` (`boundary.py:343`) | structural-red |
| **R7** | **The splitting is a TWIN.** The record advertises one splitting; the G-S driver silently discards it and re-derives another | **MEASURED**: `_within_group_si(..., inner_schedule='gauss_seidel')` → `base is record.implicit_operator` is **False**; driver runs `ScheduledInvertibleOperator` + `SNMaskedBoundaryOperator`, record still says `StreamingCollisionOperator` + `SNBoundaryOperator`. Under `jacobi` it is `True`. | structural-red |
| **R8** | `ρ` is computed and stored but **nothing reads it to decide** | `SourceIteration.contraction_ratios` populated at `numerics/iteration.py:742`; zero production consumers | ungated |
| **R9** | No wall-clock / allocation gate on any operator-composition path | only `tests/sn/sweep/core/test_cache.py:553` (`elapsed_per_sweep_ms < 2.0`, sweep-cache only) | ungated |

**R7 is the load-bearing red.** It is the exact defect the campaign
exists to fix: stage 2 (posing) and stage 3 (splitting) are welded into
one record (`WithinGroupSystem` carries `loss` AND
`implicit_operator`/`explicit_gains`), so there is no named boundary at
which "strategy may enter" can be asserted — and because there is no
boundary, a second splitting grew beside the first and the first went
silently stale.

---

## 1. The campaign acceptance criterion — executable

Replace *"changing a solve strategy must not touch operator construction"*
(green, unfalsifiable) with a **three-legged** criterion. Legs a and c are
RED today; leg b is the regression floor.

**File**: `tests/sn/architecture/test_stage_separation.py` (new package
`tests/sn/architecture/`; `pytestmark = [pytest.mark.foundation]` —
these are software/architecture invariants with no theory `:label:`,
per lessons L9 / `feedback_vv_tagging`).

### AC-a — the strategy-entry boundary is a SIGNATURE fact (RED today)

```
test_poser_signature_admits_no_strategy_parameter
```

Walk the call chain `mesh + materials → posed pencil` by
`inspect.signature`, and assert **no** function on it accepts any of the
strategy tokens:

```python
_STRATEGY_TOKENS = frozenset({
    "inner_solver", "inner_schedule", "max_iter", "max_inner",
    "tol", "inner_tol", "restart", "corrector", "preconditioner",
    "n_dof", "initial_guess",
})
```

Chain (post-carve): `realize_physics(...) → pose_k(...) / pose_alpha(...)`.
The assertion is `_STRATEGY_TOKENS.isdisjoint(params)` for each.

**Why RED today**: the chain does not exist. The single function that
produces the posed `A` (`build_within_group_system`) *also* produces the
splitting, so the gate cannot be scoped — there is no function whose
codomain is "the pencil, and only the pencil". The gate is unwritable
against `main`, which is the strongest possible statement that the
boundary is missing.

**Mutation (proves teeth after the carve)**: add `inner_schedule: str =
"jacobi"` to the poser's signature → RED with
`"strategy token 'inner_schedule' reachable from the poser"`.

### AC-b — one splitting per posed equation (RED today, R7)

```
test_driver_consumes_the_records_own_splitting[jacobi|gauss_seidel]
test_driver_consumes_the_records_own_splitting[si|krylov]
```

For every (geometry × schedule × inner_solver) row, assert the operators
the driver actually runs against are the **same objects** the partition
record advertises:

```python
assert driver_implicit is record.implicit_operator
assert tuple(map(id, driver_gains)) == tuple(map(id, record.explicit_gains))
```

**MEASURED RED today** on `gauss_seidel` (2-D Cartesian, het 2G): `base is
record.implicit_operator` → `False`. Green on `jacobi`. The `jacobi` row
is the **control leg** that pins the asymmetry (lessons L6) — without it,
a change that broke *both* arms would look like "the gate was always
partly red".

**Mutation**: after the carve, re-introduce a driver-side re-split (call
`B.split(...)` inside the driver) → RED on both rows.

### AC-c — the posed object is byte-identical across strategy (the contract)

```
test_posed_pencil_fingerprint_invariant_across_strategy
```

Fingerprint the **object**, not a functional (Mode 12 — see §2.2):

```python
def _fingerprint(op, space, seed=20260728):
    """SHA-256 of the operator's image of a fixed-seed random probe basis.

    NOT a spectrum, NOT a keff, NOT a balance sum — every one of those has
    an invariance group that contains the mutation classes this campaign
    can introduce (similarity, transpose, per-term cancellation).
    """
    rng = np.random.default_rng(seed)
    cols = []
    for _ in range(N_PROBE):                     # N_PROBE = 8
        x = space.zeros(); _randomize(x, rng)    # bulk AND trace, non-flat
        cols.append(_flat(op.apply(x)))
    return hashlib.sha256(np.stack(cols).tobytes()).hexdigest()
```

Assert one hash across the Cartesian product
`{k, alpha} × {forward, adjoint} × {jacobi, gauss_seidel} × {source_iteration, krylov}`
— **16 rows, one hash per (posing, adjointness) pair, invariant over the
4 strategy combos**.

**Why RED today**: unwritable — no pencil object (R3), and `A` and its
`M` are on incompatible carriers (R4), so there is nothing whose hash
could be taken uniformly.

**Tolerance**: this one **is** byte-identity, and legitimately so — the
strategy legs must not change a single float in the posed object, because
strategy is not supposed to reach it at all. This is the one place in the
campaign where `sha256` is the right contract (contrast §4).

### AC-b′ — the value floor (currently green; keep as regression)

```
test_reconstruction_identity_A_equals_M_minus_N[geom × schedule]
```

Retains the §0.2 probe as a permanent gate, with the honest per-arm
tolerance:

| arm | contract |
| --- | --- |
| seedless (slab / 2-D Cartesian) | `np.testing.assert_array_equal` (0 ULP — **MEASURED** exactly 0.0) |
| coupled / carrying (sphere, cylinder) | `assert_array_almost_equal_nulp(nulp=8)` (**MEASURED** 3.5e-17 ≈ 2 ULP; 8 = block-sum reduction depth with headroom) |

Mutations in §7.

---

## 2. Per-stage invariant gates

Every math-bearing type ships a test of its DEFINING law
(`feedback_test_intrinsic_properties`). One file per stage under
`tests/sn/architecture/`.

### 2.1 Stage 1 — physics realization (monomorphic leaves)

**File**: `tests/sn/architecture/test_monomorphic_leaves.py`

Parametrised over the production leaf set `{L, C, S, F, B, T}` × the
geometry ladder (slab / sphere / cylinder / 2-D Cartesian), ≥2G always.

| gate | assertion | tolerance | today |
| --- | --- | --- | --- |
| **G1.1 declared spaces** | `op.domain is not None and op.codomain is not None`, and both are concrete `FunctionSpace` instances (`type(op).domain` annotation is **not** `Optional`) | exact | **RED** — `F.domain is None` (R1) |
| **G1.2 one arrow** | for each leaf, `op.apply(x)` succeeds for exactly ONE carrier type and raises a **typed** error for every other carrier in the registry | exact | **RED** — `F` stands for 7 arrows, `S` for 6 |
| **G1.3 uniform guard** | the refusal is `TypeError` (or a named `IncompatibleCarrier`) whose message contains the operator name AND the expected space; **never** a raw `AttributeError` | `pytest.raises(TypeError, match=...)` | **RED** — `SNBoundaryOperator` raises `AttributeError` (R6) |
| **G1.4 `.H` metric-correct** | G-metric reciprocity `⟨A x, y⟩_G == ⟨x, A.H y⟩_G` on a **non-trivial** metric | `rtol=1e-12` | see below |
| **G1.5 `.H` on an anonymous instance is UNREPRESENTABLE** | constructing a leaf without a space RAISES at construction; `op.H` can never see a `None` space | `pytest.raises` | **RED** — `F.H` builds today and returns a Euclidean transpose (R2) |
| **G1.6 realizer purity** | `realize(law, mesh) ` twice → identical fingerprint (§1 AC-c helper); and `_STRATEGY_TOKENS.isdisjoint(signature)` | sha256 exact | not-yet-existing |

> ⚠ **CORRECTED 2026-07-29 — the criterion below is a PROXY, and it is wrong
> in both directions.** With `A† = G⁻¹AᵀG`, the "drop the metric" mutation is
> invisible **iff `[G, Aᵀ] = 0`** — the commutator, not mesh uniformity.
> **MEASURED:** (a) a *uniform-h* slab under `gauss_legendre(4)` still REDs
> (1.3e-1 / 4.0e-1 / 2.7e-1) because `G = V_cell·w_n` and the **quadrature
> weights** vary — so "non-uniform h" is not what buys the teeth; the truly
> blind fixture needs `G` **globally constant** (`gauss_legendre(2)`, both
> weights exactly 1, `h = 1/√3` so the bulk constant equals the trace
> constant). (b) Conversely a wildly non-uniform metric is still blind for any
> operator that COMMUTES with it: `C` is diagonal (`G⁻¹CG = C` for every
> diagonal `G`) and `B` is a specular permutation preserving `|Ω·n|·w_n`, so
> for those two leaves **no reachable config exists** — their rows are Mode-10
> *exercised-but-unconstrained* and are closed with a second, metric-agnostic
> mutation (double the adjoint → 0.5 exactly). Note also that a
> `level_symmetric` rule has **constant weights**, which alone leaves `S`/`F`
> blind. Landed: `tests/sn/architecture/test_monomorphic_leaves.py` @ `cd7b17cd`;
> generalised into `vv-principles` Mode 12.

**G1.4 config is load-bearing** (lessons L18, ERR-067). The metric MUST be
non-degenerate and non-uniform, or the reciprocity identity is blind by
construction:

- **non-uniform mesh** `h` (so `V_cell` is not a global constant — a
  constant metric cancels from both sides and the gate degenerates to a
  Euclidean transpose check, which is exactly the bug it must catch);
- **spherical or cylindrical**, so `V_cell` spans an order of magnitude;
- **≥2G with asymmetric `SigS`** (`Σ_s[g,g'] ≠ Σ_s[g',g]`), else a
  group-index transpose is invisible (Mode 12);
- for `B`: the trace metric is the cosine-weighted `|Ω·n|·w`, and the
  probe must carry **non-zero inflow on a reflecting face** — a vacuum-only
  probe nulls the reflection partner rows.

**G1.4 control leg** (the ERR-067 closure trap): the gate must run BOTH a
mutated leg (metric dropped → the identity breaks, RED) AND the unmutated
leg (`< rtol`). Without the control, a still-broken baseline mimics
"caught".

**Explicitly OUT of scope for these gates**: the *spelling* of the
dispatch (`singledispatchmethod` + `@overload` vs `match` vs a registry).
Issue **#261** parks that ruling ("deciding the spelling before the cores
move would be premature"). G1.2 gates the **contract** — one carrier
accepted, all others refused with a typed message — and is satisfied by
any spelling. No gate in this plan reads `singledispatchmethod`.

### 2.2 Stage 2 — problem posing

**File**: `tests/sn/architecture/test_posing.py`

| gate | assertion | tolerance | claim layer / pillar |
| --- | --- | --- | --- |
| **G2.1 role table** | `pose_k(leaves).A` is composed of exactly `{L, C, S, B}` and `.M` of exactly `{F}`; `pose_alpha(leaves).A` of exactly `{L, C, S, F, B}` and `.M` of exactly `{T}` — asserted on the composition **tree** (leaf identity set), not on values | exact set equality | foundation |
| **G2.2 F migrates, it does not toggle** | `(pose_alpha(l).A − pose_k(l).A)` applied to a random state `== −F.apply(x)`; AND `_STRATEGY_TOKENS ∪ {"alpha","mode","is_alpha"}` disjoint from `pose_k`'s signature (no boolean-flag α mode) | `nulp(4)` | foundation |
| **G2.3-k μ→k map, closed form** | 2G/4G homogeneous infinite medium: `pencil_k.solve().physical == k_inf` where `k_inf = λ_max(A⁻¹F)` from `orpheus.derivations` | `rtol=1e-12` | **eigenvalue claim / closed-form pillar** |
| **G2.3-α μ→α map, closed form** | 2G/4G homogeneous infinite medium: `pencil_alpha.solve().physical == α*` where `α*` is the largest root of `det(diag(Σ_t + α/v) − Σ_sᵀ − F) = 0`, computed **independently** (SymPy `nroots` / `scipy.linalg.eig` on the G×G generalized pencil `(diagΣ_t − Σ_sᵀ − F, −diag(1/v))`) | `rtol=1e-11` | **eigenvalue claim / closed-form pillar** |
| **G2.4 adjoint = a daggered row** | `pose_k(l).H.A` is `pose_k(l).A.H` **object-wise** (same leaf identities, `.H`-wrapped); AND biorthogonality `⟨φ*, F φ⟩ ≠ 0` with `⟨φ*_i, F φ_j⟩ ≈ 0` for `i≠j` | `rtol=1e-10` | flux-shape |
| **G2.5 λ has a home** | `F.apply(φ)` is `k`-free (signature has no `keff`); the `1/k` scaling is a method **on the pencil**; `compute_fission_source` no longer divides | structural | foundation |

**G2.3-α is the row this campaign must not skip.** The α posing is the
*reason* F must migrate, and there is currently **zero** α coverage. The
generalized `G×G` eigenproblem is a genuine closed-form pillar (structurally
independent — hand-posed matrices from the `Mixture`, `scipy.linalg.eig`,
no ORPHEUS operator touched). Use `get_mixture("A","2g")` and `("A","4g")`;
2G is the minimum (anti-#3) and 4G exposes the multi-root ordering.

**Mode-12 warning, binding**: `k* == k` is `eig(Mᵀ) = eig(M)` **by
construction** and carries ZERO vector information. It is a *necessary*
row only. The committed catcher for the adjoint posing is G2.4's
**biorthogonality** (a bilinear functional, outside the spectrum's
stabiliser). Do **not** let the plan credit `k* == k` as the adjoint gate.

### 2.3 Stage 3 — strategic partitioning

**File**: `tests/sn/architecture/test_partition.py`

| gate | assertion | tolerance | mutation |
| --- | --- | --- | --- |
| **G3.1 completeness** | `Σᵢ Jᵢ Rᵢ == I` — apply to a fixed-seed random state, compare to the state | `assert_array_equal` (pure gather/scatter, no arithmetic — 0 ULP is the correct contract) | drop one block from the sum → RED O(1) |
| **G3.2 biorthogonality** | `Rᵢ Jⱼ == δᵢⱼ I` for every `(i,j)` — the blocks are a genuine partition, not an overlapping cover | `assert_array_equal` | make two blocks share a DOF → the `i≠j` row becomes non-zero → RED |
| **G3.3 block decomposition** | `A == Σᵢⱼ Jᵢ Aᵢⱼ Rⱼ` | `nulp(reduction_depth)`, `reduction_depth = n_blocks²` | transpose a block index (`Aᵢⱼ → Aⱼᵢ`) → RED **only on an asymmetric operator** — see config note |
| **G3.4 splitting** | `A == M − N` (AC-b′) | per-arm, §1 | §7 mutation register |
| **G3.5 piece-splitting** | for every operator split into pieces (`B → B_lower + B_upper`), `piece_sum == whole` | `assert_array_equal` | drop a face's rows without re-complementing → RED (this gate **already exists** as `test_mutation_partition_break_reddens_split`, `tests/sn/solve/test_gauss_seidel_reification.py:383` — **generalise it, do not duplicate it**) |

**G3.3 config note (Mode 12)**: a block-index transposition on a
**symmetric** operator is invisible. The fixture MUST carry
**asymmetric `SigS`** and **non-uniform `h`** so `A ≠ Aᵀ`. Reuse the
existing lever-pulling fixture in `tests/diffusion/test_operators.py`
(its comment already states "asymmetric Σ_s so a transpose is
observable"), and the `_removal_sigmas(sn, seed=...)` helper in
`tests/sn/operators/test_removal_form_matvec_sweep.py:192` for the SN side.

### 2.4 Stage 4 — lowering + solve

**File**: `tests/sn/architecture/test_lowering.py`

| gate | assertion | tolerance | mutation |
| --- | --- | --- | --- |
| **G4.1 acyclic ⟹ one-pass exact** | if the implicit set's block-dependency graph is a DAG, `M.inverse().apply(M.apply(x)) == x` in **one** application, no iteration | `rtol=1e-13` on a **trace-consistent** state | see below |
| **G4.2 cycle ⟹ refusal** | assigning a coupling that closes a cycle to the implicit set RAISES at construction (`pytest.raises(..., match="cycle")`) | exact | remove the cycle check → the one-pass identity silently degrades to O(1) defect → G4.1 REDs (the paired teeth) |
| **G4.3 driver is strategy-only** | `SourceIteration.__init__` / `KrylovAcceleration.__init__` accept `(M⁻¹, *N, tol, max_iter, …)` and **nothing** carrying physics (no mesh, no `mat_xs`, no `Mixture`) | signature | add `mat_xs` → RED |
| **G4.4 M⁻¹ is the inverse of M** (not of a different operator) | `M.inverse()` and `M` are the two faces of ONE object: `M.inverse().apply(M.apply(x)) == x` AND `M.apply(M.inverse().apply(y)) == y` on the source subspace | `rtol=1e-13` | pair `apply` of `M` with `solve` of `M − B_lower` → O(1) defect (**this is exactly the dissolved `_GaussSeidelResolvent` bug, defect 2.667** — the existing `tests/sn/solve/test_gauss_seidel_reification.py` §13.1 gate; generalise) |

**G4.1 domain subtlety (inherited, do not re-discover)**: the sweep
substrate re-derives outflow-definition rows, so `M⁻¹` realizes the
inverse exactly on the **source subspace** `{y : y.outflow-rows = 0}`,
whose `M`-preimage is the **trace-consistent** states. The round-trip
must start from a trace-consistent state (`x.out = streamed(x.interior)`),
not an arbitrary random one. `tests/sn/solve/test_gauss_seidel_reification.py`
has the `_consistent_state(sn, LC, seed=...)` helper — **reuse it**.

---

## 3. The spectral gate — `ρ(M⁻¹N) < 1`

This is the one constraint structure cannot supply, and the campaign
makes it *more* dangerous, not less: an additive design space means many
more splittings are *expressible*, and expressibility is not stability.

**File**: `tests/sn/architecture/test_splitting_spectrum.py`
(`@pytest.mark.slow` on the `c ≥ 0.99` rows only — `pytest.param(..., marks=slow)`,
not the whole function, per lessons L9).

### 3.1 The measurand — two independent estimates, cross-checked

```python
def spectral_radius(M, N, space, *, seed, n_iter=80, anisotropic_seed=True):
    """ρ(M⁻¹N) by power iteration ON THE OPERATOR.

    Structurally INDEPENDENT of the SI driver: it never calls
    SourceIteration; it iterates x <- M.inverse().apply(N.apply(x)) with
    explicit renormalisation and returns the Rayleigh ratio.
    """
```

Cross-check against the driver's own observation:
`si.contraction_ratios[-1]` (`numerics/iteration.py:742` — **already
computed and stored, currently read by nothing**, R8). Agreement gate
`rtol=5e-2` (power-iteration truncation + finite-iteration noise). Two
structurally different observations of the same number; disagreement means
the driver is not running the splitting the record advertises (an R7-class
regression).

### 3.2 The Mode-9 trap, named exactly

**The fully-reflective isotropic box is not merely weak here — for the
σ_r fold it makes `ρ` identically ZERO.** With `N = −Σ_s0(I − P_iso)`
and an isotropic flux, `P_iso ψ = ψ`, so `N ≡ 0` and the splitting looks
*perfect* at every tolerance, in every refinement. The catalogued instance
(#215) is precisely this: exact on the reflective uniform box, 46–56 %
wrong on vacuum/heterogeneous.

Therefore, **binding**:

- the ρ measurement MUST seed the power iteration with an
  **angularly anisotropic** vector (`anisotropic_seed=True`: project OUT
  the isotropic component of the seed, `x ← x − P_iso x`, before
  iterating). An isotropic seed leaves the iteration inside an invariant
  subspace on which `ρ = 0` — the trap in its exact algebraic form;
- the ADMISSION rows run on **heterogeneous 2G + P1 anisotropic
  scattering + mixed vacuum/reflective BC**. `Mixture(...)` must be
  constructed **directly** — `make_mixture` hardcodes `SigL = 0` and
  defaults `Sig2 = 0` (lessons L1), and an all-isotropic library mixture
  nulls the very `(I − P_iso)` term the gate exists to measure;
- an **axis-aligned `product` quadrature is a second degeneracy** for
  angular-schedule splittings (ERR-056). Where the splitting touches an
  angular schedule, use `level_symmetric` or `lebedev`.

### 3.3 The rows

| row | config | claim | assertion | tol |
| --- | --- | --- | --- | --- |
| **S1 analytic anchor, standard SI** | infinite homogeneous (all-reflective), **1G legitimate here** — see claim-layer note | `ρ = c` | `ρ_measured == c` | `rtol=2e-2` |
| **S2 analytic anchor, σ_r fold** | same, **anisotropic seed** | `ρ = c/(1−c)` — the spatially-flat / angularly-anisotropic mode: `L = 0`, `M = σ_r = σ_t(1−c)`, `N = Σ_s0 = cσ_t` ⟹ `ρ = c/(1−c)` | `ρ_measured == c/(1−c)`, sweep `c ∈ {0.2, 0.4, 0.6, 0.8, 0.9}` | `rtol=3e-2` |
| **S3 threshold** | same | `c = 1/2 ⟺ ρ = 1` | `ρ(c=0.45) < 1 < ρ(c=0.55)` — a **bracketing** assertion, not a point value | strict inequality |
| **S4 isotropic-seed CONTROL** | same, `anisotropic_seed=False` | the trap is real and the gate knows it | `ρ_measured < 1e-12` at `c = 0.9`, with `err_msg` naming vv Mode 9 | exact |
| **S5 production divergence** | **S8 slab, heterogeneous 2G, vacuum, P1, `c = 0.9`** | the σ_r fold diverges in production geometry | `ρ > 1` (measured ≈ 6.91 with leakage; the infinite-medium 9.0 is an upper bound — assert `1 < ρ < c/(1−c)`, both sides) | strict |
| **S6 admission** | every splitting the production factory can emit × {slab, sphere, cyl, 2-D Cartesian} × het 2G × P1 × mixed BC | the factory refuses an unstable splitting | `ρ < 1 − MARGIN` (`MARGIN = 0.05`) OR the factory RAISES `UnstableSplitting` | strict |
| **S7 FP-invariance** (paired) | het 2G, vacuum, P1, `lebedev` | a splitting changes the RATE, never `ψ*` | `ψ*(split) == ψ*(plain Jacobi)` at `SAFETY(50) × inner_tol` | see L7 |

**Claim-layer declaration for S1–S4**: these are **convergence-RATE**
claims, not eigenvalue claims. `ρ` is flux-shape-independent by
construction, so 1G is *legitimate* here and the Cardinal Rule (which
bars 1G **eigenvalue** claims) does not apply. This must be stated in the
docstrings, or a reviewer will correctly flag it as an anti-#3 violation.
S5–S7 are ≥2G because they are production-regime claims.

**S4 is the gate's own teeth.** It asserts the Mode-9 trap is *present*
and that the harness's anisotropic seeding is what defeats it. If S4 ever
goes green with `ρ > 0`, the seeding logic changed and S2's anchor is no
longer measuring what it claims.

### 3.4 The `UnstableSplitting` refusal — the production consequence

The gate is only worth its cost if production *acts* on ρ. Recommend
(and gate) that the partition factory carry an **admission check**:

```
test_factory_refuses_unstable_splitting
```

`pytest.raises(UnstableSplitting, match="rho=")` when handed the σ_r fold
at `c = 0.9`, **paired with the positive control** (the same factory,
`c = 0.2`, constructs successfully). Without the control, a factory that
raises on *everything* passes (the L4 `pytest.raises` false-green).

---

## 4. Regression strategy — anchors, and where bit-identity is the WRONG contract

### 4.1 Load-bearing anchors that MUST stay green AND bit-identical

| path | what it pins | contract |
| --- | --- | --- |
| `tests/sn/operators/test_scattering_kernel_crosscheck.py` | the fused ndarray scattering kernel as a **0-ULP oracle** | `array_equal` — **unchanged**. The fused kernel is NOT to be deleted (explicit campaign constraint); it is the fuller-view oracle per `coding-standards.md` |
| `tests/sn/regression/test_dd_regression.py`, `test_walk_matvec_baselines.py`, `tests/sn/regression/conftest.py` | committed snapshots under `-W error::DriftWarning` | unchanged; **but** see §4.3 — an *iterated* snapshot is the wrong bit-identity gate (L7) |
| `tests/sn/solve/test_affine_carve_bit_identity.py`, `tests/sn/sweep/core/test_affine_carve_baseline.py` | the affine carve's byte-identity | unchanged |
| `tests/sn/operators/test_removal_form_matvec_sweep.py` | the σ_r removal-form value gates across slab/sphere/cyl/2-D, ≥2G, vacuum | unchanged — this is the σ_r *value* anchor the spectral gate's S5 sits beside |
| `tests/sn/solve/test_gauss_seidel_reification.py` §13.1/§13.2/§13.3 | `M`'s two faces; `B == B_lower + B_upper`; G-S ≡ Jacobi fixed point | **generalise** into G3.5 / G4.4 / S7 — do not duplicate, do not retire |
| `tests/sn/verification/analytical/test_si_convergence_rate.py` | SI iteration count vs analytic ρ | unchanged; the spectral gate's driver-side cross-check reuses its measurand |
| `tests/sn/eigenvalue/test_keff_2d.py:595` | jacobi-vs-G-S `k_eff` invariance | unchanged |
| `tests/sn/acceleration/test_dsa_low_order.py` | the DSA low-order build (`A_diff` stencil) | unchanged — the campaign must not perturb #2's landed gates |
| `tests/homogeneous/` `test_kinf_exact` | the closed-form `k_inf = λ_max(A⁻¹F)` | unchanged — **this is the structurally-independent anchor for G2.3-k** |

### 4.2 PRINCIPLED re-baseline (NOT bit-identity) — with the 3 criteria discharged

Per `vv-principles` §bit-identity, each of the following changes the FP
reduction tree deliberately. For each, all three criteria are checked:
**(1)** every new intermediate is a named domain quantity, **(2)** a
structurally-independent reference pins the new value, **(3)** the drift
is dimensionally bounded.

| what changes | why not bit-identical | new contract | independent reference (crit-2) | drift bound (crit-3) |
| --- | --- | --- | --- | --- |
| `A_ij = Rᵢ A Jⱼ` replacing the fused `A_AA = LC − S − B_a` | the block composition sums per-block instead of one flat operator sum | `nulp(n_blocks²)` | closed-form `Q/Σ_t` flat fixed-source **and** `k_inf` (two anchors, L2) | `reduction_depth × ULP` |
| `1/k` migrating from `compute_fission_source` INTO the pencil | `F.apply(φ)/k` vs `(F·(1/k)).apply(φ)` — one rounding moves | `nulp(2)` on the fission source; `rtol=1e-12` on `k_eff` | homogeneous `k_inf` closed form | 1–2 ULP per element |
| the coupled-arm `A == M − N` | **already** 3.5e-17 on `main` (**MEASURED**) — the block grid re-associates | `nulp(8)` from day one | the seedless arm's exact-0 row is the contrast | block-sum depth |
| leaf-arrow collapse (7→1 on `F`, 6→1 on `S`) | the surviving arm may not be the one the current production caller hits | `nulp(4)` on `F.apply` / `S.apply`, **anchored** | the retained fused kernel (`test_scattering_kernel_crosscheck.py`) for `S`; `k_inf` for `F` | reduction depth |

**Anti-pattern to refuse**: do **not** propose `array_equal` across the
whole campaign. Two of the four rows above are already non-zero on `main`.

### 4.3 The bit-identity proof must descend below the iterated snapshots

Committed *iterated* DD snapshots already drift thousands of ULP
cross-run under `-W error::DriftWarning` (lessons L7). For any leg of
this campaign claiming "zero numerical change", the proof is a
**single-step DIRECT** snapshot on a fixed-seed random heterogeneous ≥2G
ψ with non-zero inflow — captured **pre-carve** via the root conftest's
`--capture-baseline` flag. Flat ψ nulls redistribution and is not
admissible as the bit-identity probe.

### 4.4 A `catches(ERR-NNN)` audit is owed

> ⚠ **CORRECTED 2026-07-29 — do NOT file a new ERR for the σ_r fold.** This
> section originally read "the σ_r fold (#215) has **no ERR entry** — the
> catalog ends at ERR-069 … file **ERR-070**". That was already false when
> written: **ERR-070 *is* the σ_r fold**, filed 2026-07-26 by the #2 DSA
> battery's seeded mutation (Phase 3c, measured 43.2 % fixed-point shift),
> and ERR-071 is taken too. Executing the original instruction would have
> minted a duplicate ID for a bug that already had one — the plan asserted a
> catalog state it had not read (lessons L33: a doc is a claim, not evidence).

ERR-070's catchers on record are a **value** gate
(`test_dsa_rate.py::TestSigmaRFoldCaught`) and a **structural fence**
(`TestD10RoutingSentinel`, an AST sweep of the foldable accessors'
consumers). This campaign adds a third, **algebraic** one —
`test_stage_separation.py::test_the_sigma_r_fold_is_a_splitting_only_with_its_anisotropic_remainder`
(landed @ `9a546640`), which reds the moment the splitting is *constructed*,
before any solve runs.

So the standing obligation is **narrower than "file an ERR"**: when S5/S6
first catch a σ_r-class regression, attach `@pytest.mark.catches("ERR-070")`
to S5 **only after** mutation-verifying that S5 *specifically* reddens — a
marker is a coverage claim, and an unverified one is a phantom. Do not
attach it to any gate that merely characterises the degeneracy without
detecting the defect.

---

## 5. The performance gate — MANDATORY, with precedent

### 5.1 The precedent

A refactor moved hoisted vectorised work into a per-cell Python fold and
cost **10–20× on slab / ~6× on the suite** with **every correctness gate
green** (`.claude/agent-memory/test-architect/wavefront_flux_carve_lessons.md:64`;
`method-implementer/issue_196_phase_g_step2_5b_closeout.md:110` records the
30.85 s baseline). `A_ij = Rᵢ A Jⱼ` as a composition is the same shape of
risk: each `Rᵢ`/`Jⱼ` is an opportunity to materialise a slice per block
per cell.

### 5.2 The gate has three legs; the STRUCTURAL leg is the real catcher

**File**: `tests/sn/architecture/test_composition_cost.py`

**P-1 — call-count scaling (the catcher; cannot be flaky, cannot be gamed)**

An in-process wrap (autouse fixture / `monkeypatch`) on the **leaf
kernel** counts entries. Then refine the mesh and assert the count does
**not** scale with `n_cells`:

```python
counts = {}
for nx in (20, 40, 80, 160):
    with count_calls(DiamondDifference, "update_batch") as c:   # the batched leaf
        A.apply(x)                                              # ONE apply
    counts[nx] = c.n
# The batched kernel is called O(n_ordinates * n_groups) times, NOT O(nx).
assert counts[160] == counts[20], (
    f"leaf-kernel call count scaled with the mesh: {counts} — a per-cell "
    f"Python fold has replaced a hoisted vectorised call (L16 regression)"
)
```

This is a **Mode-11-style wrap**: it proves the changed path is executed
AND that its call arity is structural, not per-cell. It is deterministic
— no timing noise, no machine dependence, no CI flake.

**P-2 — throughput floor, normalised against an in-process calibration**

Absolute ms thresholds are machine-dependent (the existing
`test_cache.py:553` `< 2.0 ms` gate is the precedent but also the
cautionary tale). Normalise:

```python
calib = time_of(lambda: np.einsum("ij,j->i", Aref, xref), n_reps)  # same FLOP order
t_apply = time_of(lambda: A.apply(x), n_reps)
assert t_apply / calib < COMPOSITION_OVERHEAD_MAX   # start at 12.0
```

`COMPOSITION_OVERHEAD_MAX` is a **claim**, not a magic number (L7): it is
"the block composition costs at most N× the raw contraction of the same
FLOP count". Pin the pre-carve value in the same run and record it in the
docstring. Use `min` of `n_reps` timings, never the mean.

**P-3 — allocation ceiling**

```python
tracemalloc peak bytes for ONE A.apply(x), divided by n_dof * 8
assert peak_per_dof < ALLOC_FACTOR_MAX      # start at 6.0
```

Catches the "materialise a slice per block per cell" mechanism directly,
and unlike wall-clock it is exactly reproducible.

**P-4 — the pre-carve baseline capture**

Before stage 3 lands, run P-1/P-2/P-3 on `main` and commit the numbers
into the test module as named constants with a provenance comment
(`# measured on main @ <sha>, host .venv Py 3.14`). A perf gate whose
baseline was measured *after* the regression is worthless.

### 5.3 Where it must sit in the phase order

> **Stage 3 (strategic partitioning) MUST NOT MERGE without P-1, P-2,
> P-3 green and P-4's baseline captured on the immediately preceding
> commit.**

Stage 3 is where `A_ij = Rᵢ A Jⱼ` first exists, therefore where the
regression can first occur, therefore where the gate must already be
standing. Stages 1–2 change types and composition, not the inner loop;
stage 4 consumes what stage 3 built. **Capture P-4 at the end of stage 2.**

**Run it SERIAL and UNCONTENDED.** xdist is unstable on this venv, and a
contended wall shows up as a false red on the timing legs (P-2). P-1 and
P-3 are contention-immune — if the timing leg is the only red, re-run it
alone before believing it.

---

## 6. Phase-ordering constraints FROM the verification side

The rule: **no phase may land in a state where its own gate is
unwritable.** Each row names what must EXIST before the implementation
step, so the step is verifiable at merge.

| # | Must exist BEFORE | implementation step | why (what becomes unverifiable otherwise) |
| --- | --- | --- | --- |
| **O-1** | AC-b (R7 twin-splitting gate) **written and RED**, with the `jacobi` control leg | any touch to `WithinGroupSystem` or `_select_si_splitting` | R7 is the campaign's motivating defect. If it is fixed *before* the gate is written, the fix is unprovable and a future re-split cannot be caught. Write the red gate first; the fix flips it green. |
| **O-2** | G1.1/G1.2/G1.3 (leaf monomorphism + uniform guard) | the leaf-arrow collapse (7→1 on `F`, 6→1 on `S`) | The collapse RETIRES arms. Retirement without a per-arm refusal gate loses coverage silently (L4: retiring a guard with no negative test makes the successor's teeth net-new, not migrated). Grep first: does ANY test assert the current 7 arrows? |
| **O-3** | G1.5 (`.H` on an anonymous instance RAISES) + G1.4's control leg | making `domain`/`codomain` non-`Optional` | The `None`-space `.H` degradation (R2) is silent today. If the space becomes mandatory before the gate exists, nothing proves the degradation was ever possible — the ERR entry loses its catcher. |
| **O-4** | G2.3-α's independent α reference (`scipy.linalg.eig` on the hand-posed `G×G` pencil) | the F-migration (`pose_alpha`) | MMS cannot prove an eigenvalue (pillar rule). If F migrates before the α closed form exists, the α posing has **no** reference at all — the first α number the code produces would become its own baseline. Build the reference first; it is ~15 lines. |
| **O-5** | G2.1's leaf-identity-set assertion | the pencil type | The pencil is the first object whose *composition* is the claim. Without the tree-level set gate, a leaf silently dropped from `A` is caught only by a value gate, and a dropped `B` moves the value by ~5 % (**MEASURED**: 5.15e-2) — inside some tolerances. |
| **O-6** | **P-4 baseline captured** (§5.2) | stage 3 (`A_ij = Rᵢ A Jⱼ`) | §5.3. Hard blocker. |
| **O-7** | G3.1/G3.2 (`ΣJᵢRᵢ = I`, `RᵢJⱼ = δᵢⱼ`) | the first block-composed operator | These are the partition's *defining laws*. A block algebra shipped without them is a set of arrays that happen to be indexed. |
| **O-8** | S1/S2/S3/S4 (the ρ harness incl. the **anisotropic seed** and the S4 control) | the first splitting the factory can emit BEYOND `{jacobi, gauss_seidel}` | The additive design space is the campaign's *payoff* and its *risk*. The moment a third splitting is expressible, ρ becomes the only thing standing between "expressible" and "convergent". S4 must be green (ρ≈0 under an isotropic seed) before S2's anchor is trusted. |
| **O-9** | S6 admission + `UnstableSplitting` positive control | wiring the factory to ANY user-selectable splitting | A user-facing knob that can silently diverge is worse than no knob. |
| **O-10** | G4.2 (cycle ⟹ refusal) | the schedule realizer | G4.1's "one-pass exact" claim is only meaningful if the acyclicity precondition is enforced; the two are one gate split in half. |
| **O-11** | AC-c fingerprint helper (`_fingerprint`) | any phase claiming strategy-invariance | It is the campaign's contract. It must be a shared helper (one home), not re-spelled per test — a drifting fingerprint helper is a Pattern-2 twin inside the gate itself. |
| **O-12** | The §4.2 re-baseline table ratified by the user | the first non-bit-identical merge | `feedback_principled_over_bit_identical`: re-baselining is a *decision*, made once, in the plan of record — not a per-test tolerance nudge at merge time. |

**Ordering summary**: `O-1 → (O-2, O-3) → (O-4, O-5) → O-11/O-12 →
O-6/O-7 → O-8/O-9/O-10`. Stages 1 and 2 can proceed in parallel with each
other only after O-2/O-3 land; stage 3 is hard-blocked on O-6.

---

## 7. Mutation register — every gate's teeth, with the expected RED signature

Every mutation is applied **in-process by `monkeypatch`**. NEVER
`git checkout` a file carrying uncommitted edits
(`process-discipline.md`). Every assertion uses `np.testing.assert_*` /
`pytest.fail` / explicit `raise` — **never a bare `assert` in a production
or helper module** (Mode 8: `-O` strips it; bare asserts in *collected*
`tests/` modules are pytest-rewritten and DO fire, so they are admissible
inside test bodies only).

| id | gate | exact mutation | expected RED signature |
| --- | --- | --- | --- |
| **M-1** | AC-a | add `inner_schedule: str = "jacobi"` to the poser signature | `strategy token 'inner_schedule' reachable from the poser` |
| **M-2** | AC-b | re-introduce a driver-side `B.split(...)` re-derivation | `driver_implicit is not record.implicit_operator` on BOTH schedule rows (today: red on `gauss_seidel` only) |
| **M-3** | AC-c | change the `restart` used by the Krylov leg | must **NOT** red (strategy may change the algorithm) — the **negative control** proving AC-c fingerprints the *operator*, not the run |
| **M-4** | AC-c | drop `B` from the posed `A` | fingerprint hashes differ across the k/α rows; `assert_array_*` companion moves 5.15e-02 (**MEASURED**) |
| **M-5** | AC-b′ | set `N = 0` (empty `explicit_gains`) | seedless: `‖Ax−Mx‖/‖Ax‖` jumps to O(1) — the `S` term is the bulk of `A`; expect ≫ 0.5 |
| **M-6** | AC-b′ | **the exact σ_r bug**: replace `N = −Σ_s0(I − P_iso)` with `N = −Σ_s0·I` | `A ≠ M − N` on an **anisotropic** state by O(Σ_s0·‖(I−P_iso)ψ‖); **MUST be ≈ 0 on an isotropic state** — ship BOTH legs, the isotropic leg is the Mode-9 control that proves the anisotropic leg is what caught it |
| **M-7** | AC-b′ | flip the sign of one gain (`N → −N` on one block) | sphere carrying arm: 8.569e-03 (**MEASURED**) |
| **M-8** | G1.1 | restore `domain -> FunctionSpace \| None` on one leaf and build it anonymously | `domain is None` → RED (this is **today's state**; the gate's first run is its own mutation proof) |
| **M-9** | G1.3 | make one leaf's carrier refusal a raw `AttributeError` again | `pytest.raises(TypeError, match=...)` fails with `DID NOT RAISE`/wrong type (this is **today's state** for `SNBoundaryOperator`) |
| **M-10** | G1.4 | drop the metric application in `_AdjointOperator.apply` (return `inner.apply_transpose(y)` bare) | reciprocity residual jumps from `<1e-12` to O(V_cell spread) — **only on the non-uniform-h curvilinear config**; the uniform-h leg stays green, which is the config-blindness proof |
| **M-11** | G2.2 | make α a boolean flag on `pose_k` instead of a separate posing | signature gate REDs on `{"alpha","mode","is_alpha"}` |
| **M-12** | G2.3-α | drop the `1/v` weighting from `M_alpha` | α value diverges from the independent `scipy.linalg.eig` root by O(1) |
| **M-13** | G2.4 | replace the adjoint pencil with `A` un-daggered (`F† = F`, `S† = S`) | `k* == k` **stays green** (Mode 12) — biorthogonality REDs. Ship the `k* == k` row explicitly labelled *necessary-not-sufficient* so a reviewer cannot credit it |
| **M-14** | G3.1 | drop one block from `Σ Jᵢ Rᵢ` | `assert_array_equal` fails on that block's DOF slice (O(1)) |
| **M-15** | G3.2 | make two blocks share one DOF | `Rᵢ Jⱼ` for `i≠j` is non-zero on that DOF |
| **M-16** | G3.3 | transpose a block index `Aᵢⱼ → Aⱼᵢ` | REDs **only** with asymmetric `SigS` + non-uniform `h`; the symmetric leg stays green — ship the symmetric leg as the Mode-12 control |
| **M-17** | G3.5 | drop one face's rows from `B_lower` without re-complementing | `B ≠ B_lower + B_upper` (the existing `test_mutation_partition_break_reddens_split` — reuse verbatim) |
| **M-18** | G4.2 | remove the cycle check | G4.1's one-pass identity degrades to O(1) defect |
| **M-19** | G4.4 | pair `M.apply` with `(M − B_lower).solve` | round-trip defect ≈ 2.667 (the documented `_GaussSeidelResolvent` falsifier) |
| **M-20** | P-1 | replace the batched leaf call with a per-cell loop | `counts[160] != counts[20]` — call count scales 8× with the mesh |
| **M-21** | S2 | seed the power iteration isotropically | `ρ → 0` — this IS S4, shipped as a permanent gate rather than a one-off mutation |
| **M-22** | S6 | hand the factory the σ_r fold at `c = 0.9` | `UnstableSplitting` raised (positive), and the `c = 0.2` control constructs (negative) |

---

## 8. Mandatory configuration — one table, checked per gate

Per `AGENT.md` §0.6, do not rely on memory: every gate below names the row
it uses. The convenient config nulls the exact term the campaign is most
likely to get wrong.

| lever | setting | which blindness it breaks |
| --- | --- | --- |
| groups | **≥2G always** for value/eigenvalue claims (4G for G2.3); 1G permitted ONLY for the S1–S4 **rate** claims, with the claim layer declared | `k = νΣf/Σa` is flux-shape-independent (anti-#3) |
| `SigS` | **asymmetric** (`Σ_s[g,g'] ≠ Σ_s[g',g]`) | a group-index transpose / a block-index transpose is invisible on a symmetric operator (Mode 12) |
| mesh | **non-uniform `h`**, **non-square** `nx ≠ ny` in 2-D | an `x↔y` swap; a constant metric cancelling from the reciprocity identity |
| materials | **heterogeneous** (≥2 regions) | flat flux nulls redistribution and every spatial-distribution bug (H2) |
| scattering | **P1 anisotropic**, `Mixture(...)` constructed **DIRECTLY** | `make_mixture` hardcodes `SigL = 0` and defaults `Sig2 = 0`; an all-isotropic fixture nulls `(I − P_iso)` — the exact term the σ_r gate measures (L1) |
| BC | **mixed** vacuum + reflective; at least one reflecting face carrying non-zero inflow | a vacuum-only probe nulls the reflection-partner trace rows `B` acts on |
| quadrature | `level_symmetric` or `lebedev` for anything touching an angular schedule | an axis-aligned `product` quad gives each face one octant — the ERR-056 degeneracy |
| geometry | slab **and** sphere **and** cylinder **and** 2-D Cartesian; slab-only is never sufficient | slab is the degenerate curvilinear case (angular redistribution is a `ZeroOperator`) |
| probe state | fixed-seed **random, non-flat**, bulk AND trace, non-zero inflow | flat ψ nulls the streaming coupling; a zero trace nulls `B` |
| power-iteration seed | **anisotropic** (isotropic component projected out) | `N ≡ 0` on the isotropic subspace for the σ_r fold — ρ measures exactly 0 (§3.2) |

---

## 9. What this plan deliberately does NOT propose

- **Not** bit-identity as a universal contract. §4.2 lists four rows that
  are deliberately non-bit-identical, each with its three `vv-principles`
  criteria discharged; two of them are already non-zero on `main`.
- **Not** the fully-reflective isotropic box as the spectral fixture. §3.2
  shows it drives `ρ` to *identically zero* for the σ_r fold — the Mode-9
  trap in closed form. It appears only as the S4 **control**.
- **Not** the retirement of Pattern M (`singledispatchmethod` + `@overload`).
  #261 parks that ruling; G1.2 gates the one-arrow **contract**, satisfied
  by any spelling. No gate here reads the dispatch mechanism.
- **Not** the deletion of the fused ndarray scattering kernel. It is the
  0-ULP oracle pinned by
  `tests/sn/operators/test_scattering_kernel_crosscheck.py` and it stays,
  as the crit-2 independent reference for the `S`-arrow collapse (§4.2).
- **Not** a `keff(assembled) == keff(apply)`-style gate anywhere. Every
  scalar-functional gate in this plan is labelled necessary-not-sufficient
  and paired with an object-level catcher (Mode 12).
