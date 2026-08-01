# B3.4a test migration — working notes

Branch `refactor/operator-strategy-layers`. Scope: 4 files, 28 measured
failures (brief said 27). All four now GREEN under `python -O -m pytest`:
`110 passed`.

## Production surface (read, 2026-08-01)

- `AngularAverageOperator(cos_w, norm, n_inflow)` — `cos_w` = restriction to
  Γ₊, strictly positive; `n_outflow` = domain size; `n_inflow` = codomain
  size; `n_ordinates` RETIRED. `apply(psi)` requires
  `psi.shape[0] == n_outflow`, returns `(n_inflow,) + psi.shape[1:]`.
- `from_quadrature(quad, axis, outward_sign)` routes through
  `build_omega_dot_n(quad, (face,))` ⇒ a `z` face on GL now raises
  **`Face 'zmax' requires genuine mu_z cosines`**, not `no outgoing`.
- `IncomingSourceOperator(source, *, n_inflow)` — mask RETIRED.
- `SNBoundaryRealizer`: white + prescribed_inflow now call
  `_outflow_restriction` (so `SNMethodSpace.minimal` refuses them); white
  additionally passes `_checked_angular_average` (ERR-041-pattern orientation
  cross-check).

## Measured quadrature facts driving the design

| quadrature | face | N | \|Γ₊\| | \|Γ₋\| | tangential | retired `>0.0` set | agree? |
|---|---|---|---|---|---|---|---|
| lebedev(17) | xmax | 110 | 49 | 49 | 12 | 49 | yes |
| lebedev(11) | zmax | 50 | 21 | 21 | 8 | 21 | yes |
| lebedev(5) | xmax | 14 | 5 | 5 | 4 | 5 | yes |
| level_symmetric(4) | ymin | 24 | 12 | 12 | 0 | 12 | yes |
| level_symmetric(6) | zmax | 48 | 24 | 24 | 0 | 24 | yes |
| gauss_legendre(8) | xmax | 8 | 4 | 4 | 0 | 4 | yes |
| **product(2,4)** | **xmax** | **8** | **2** | **2** | **4** | **4** | **NO** |
| product(2,4) | ymin | 8 | 2 | 2 | 4 | 2 | yes |

- `product(2,4)` @ `xmax` is the ONLY fixture in the tree where the retired
  `> 0.0` classifier disagrees with `> TANGENTIAL_EPS`. It is now the
  activation-guarded discriminator in `test_domain_is_the_tangential_band_gamma_plus`.
- `Σ_{Γ₊} w|Ω·n| == Σ_{Γ₋} w|Ω·n|` BIT-IDENTICALLY on every row above. That
  face symmetry is what makes the re-posed conservation + reciprocity claims
  exact; it is asserted inline as an activation guard.
- `|Γ₊| == |Γ₋|` on EVERY row ⇒ Mode-12: a shape comparison cannot
  distinguish Γ₋ from Γ₊.

## LATENT BUG the new orientation guard surfaced

1. **`tests/sn/operators/test_operator_block_role.py::TestBoundaryLeaves`** —
   `_boundary_method_space()` builds `face="xmin"` while `_LINEAR_LAWS` carries
   `WhiteBoundary(axis="x", outward_sign=+1)` (= `xmax`). MEASURED raise:
   `WhiteBoundary declares axis='x', outward_sign=+1 — i.e. the face 'xmax',
   whose Γ₊ is [2, 3] — but it is being installed where Γ₊ is [0, 1]
   (face='xmin')`. Pre-B3.4a this realized a Lambertian that averaged over the
   INSTALLATION face's INFLOW hemisphere and reported nothing. OUT OF SCOPE
   (5th file) — reported, not fixed.
2. **`orpheus/transport/method.py::_law_from_tag`** — every parameter-free law
   is built as `law_cls()`, so `WhiteBoundary` always gets the default
   `axis="x", outward_sign=+1`, irrespective of the face being resolved
   (reflective correctly uses `AXIS_NAMES[label.axis_index]`). Not reachable
   TODAY: `SNMesh.BOUNDARY_OPERATOR_REGISTRY == {"vacuum", "reflective"}`, so
   `BC("white")` refuses at parse ("SNMesh does not support boundary
   condition 'white' on face 'xmin'"). Fires the day #189 admits white.

## Mutation matrix (`scratch/b3_4a_mutations.py`, in-process plugin)

CONTROL = `110 passed`. Every variant below is a plausible transcription of a
real B3.4a hazard, installed by monkeypatch (no file touched).

| id | wrong code simulated | tests reddened |
|---|---|---|
| W1 | retired `> 0.0` Γ₊ classifier | 1 — `test_domain_is_the_tangential_band_gamma_plus` |
| W2 | `norm = Σ w` (drops \|Ω·n̂\|) | 5 — 3 value + 2 conservation |
| W3 | codomain sized from Γ₊ | **0 — DESIGNED-GREEN (Mode 12), see below** |
| W4 | `apply` broadcasts over the INPUT rows | 1 — `test_codomain_size_follows_n_inflow_not_the_input` |
| W5 | `apply` drops its domain guard | 2 |
| W6 | `cos_w` guard reverts to `>= 0` | 1 — `test_masked_full_face_cos_w_is_refused` |
| W7 | face-sign-dependent `norm` | 3 — both reciprocity rows + the ymin value row |
| Q1 | `q` sized from the INPUT (ERR-047 shape) | 5 |
| Q2 | `n_inflow` guard dropped | 1 |
| R1 | white realizer skips the orientation check | 2 (positive control stays GREEN) |
| R2 | white realizer skips `_outflow_restriction` | 1 |
| R3 | prescribed-inflow skips `_outflow_restriction` | 1 |

W7 is the discriminator proving the reciprocity gate carries teeth the
`+1`-face conservation gates do not.

Q1 initially reddened only 4 — `test_delivered_q_has_no_row_off_the_incoming_trace`
stayed green because its probe was Γ₊-sized and `|Γ₊| == |Γ₋|`. A FULL-FACE
probe leg was added; Q1 now reds it, so the class's `catches("ERR-047")` claim
is LIVE rather than a phantom.

## RESIDUAL BLIND SPOT (honest, irreducible today)

W3 — `from_quadrature` sizing the codomain from Γ₊ instead of Γ₋ — reddens
NOTHING and cannot: `|Γ₊| == |Γ₋|` on every quadrature × face pair in the
tree, so the error class lies inside both the value and the shape functional's
invariance group. The OPERATOR-level twin (W4) IS caught, via a
hand-constructed unequal-size fixture. Making the `from_quadrature` half
structural is B5's job (typing the kernel as `u ⊗ v` with `u` a vector ON the
Γ₋ space).
