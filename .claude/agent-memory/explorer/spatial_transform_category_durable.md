---
name: spatial-transform-category-durable
description: The SPATIAL layer's geometric transformations are Euclidean-group (O(d) x translation) and permutation/sign-realized, NOT O(3) point-group elements — what is expressible in symmetry.py, what is not, and why the translation is the gap.
metadata:
  type: reference
---

# Spatial geometric transformations — the durable category map

Audited 2026-08-03 on `refactor/operator-strategy-layers` (measured, not
read). Structure survives churn; re-derive `file:line` via Nexus.

## The category verdict (the thing to remember)

The spatial layer's transformations are **NOT `O(3)` point-group elements**
in `orpheus/numerics/symmetry.py`'s sense. They split into four kinds:

1. **Index-lattice order reversals** (`arange(n)[::-1]`, `range(...,-1,-1)`,
   `[::-1]` slices) — the reflection `σ_a` of the *cell* lattice.
2. **A face-role involution on a 2-set** (`face_opposite`, `face_in`/`face_out`,
   `select_outer: bool`, the `0 ↔ n_a` edge-index flip).
3. **Sparse axis-aligned normals as a name↔(axis, sign) bijection** — the
   outward normal is never a vector; `Ω·n̂ = sign·μ_axis`.
4. **Integrated Jacobians of a coordinate chart whose chart map is never
   written down** — there is no `x = r cosθ` in `orpheus/geometry/`.

## The measured gap: the TRANSLATION, not the reflection

`symmetry._orbit_closure` on a **centred** 1-D cell lattice returns exactly
the `arange(n)[::-1]` permutation the sweep DAG writes by hand — the
machinery already produces the spatial answer. What kills it is that
production meshes start at `origin = 0` (`Mesh1D.from_geometry`), so
`Mirror('x').is_invariant(mesh.volume_measure)` is **False** for the
production slab AND the sphere. The spatial group is
**`E(d) = O(d) ⋉ ℝ^d`**; `symmetry.py` realizes only the origin-fixing
`O(3)`. The codebase has **no affine/translation type at all** (`origin:
float`, `- pitch/2` inline, the half-cell's symmetry plane is "wherever
x = 0 happens to be").

Sharper reading for CYL/SPH: `r ∈ [0,R]` is the **fundamental domain** of a
quotient, and a fundamental domain is never G-invariant — `False` is the
CORRECT answer. The spatial layer needs a quotient / fundamental-domain
vocabulary, not an invariance predicate.

## The two genuine group objects in the spatial layer

- **`octant_moment_frame_signs`** (`transport/spatial/_ubld.py`) is the exact
  1-D **character representation of `(Z₂)^d`** on the tensor-Legendre spatial
  moment basis — verified a group homomorphism over all `2^d × 2^d` pairs at
  d=1,2,3. Docs call it "an involution" (true per element, undersells the
  object). It was BORN from a bug (ERR-061, the LD diffusion-limit sign
  cancellation) that the group framing would have made unspellable.
- **`quad.reflection_index("x")`** consumed at the r=0 pole handoff
  (`sn/loss_representation/__init__.py`, 4 sites) realizes the *deck
  transformation* of the radial quotient, `ψ(0,+μ) = ψ(0,−μ)`. Measured to
  BE `Mirror("x")` to 0 ULP on gauss_legendre / product / level_symmetric /
  lebedev: bijection, involution, weight-preserving.

## The one place a transformation is USED rather than open-coded

The adjoint (`_reverse_octant_traversal` + `_CellResidualTranspose`) maps the
octant label `o ↦ −o`, **looks up the mirror octant's graph in the same cached
family**, and runs the UNCHANGED forward walk. Verified: `face_in(−o) ==
face_out(o)` exactly at d=1,2,3; `levels(−o) == reverse(levels(o))` **as SETS**
at all d, **elementwise only at d=1** (immaterial — within-level cells are
independent and gather/scatter share one index array). Contrast: the
sweep-direction reversal is spelled **eleven** times overall, four of them
within ~100 lines of one file, one as a literal `x_reverse: bool` parameter.

## Adjacent durable facts

- Coordinate-system identity has **four** parallel encodings: `CoordSystem`
  enum, `AxisCoord` enum, the bare string `SNMesh.curvature ∈ {None,
  "spherical", "cylindrical"}`, and the loss-rep's re-normalization that adds
  `"cartesian"`.
- **No hand-built 2×2/3×3 transformation matrix exists in the spatial layer.**
  The only literal small matrices are LD element matrices (`_GRAD_1D`,
  `_FOUT_1D`). Every transformation matrix in the repo is in `symmetry.py`.
- **No octant/wedge/sector mesh exists.** The only symmetry-reduced meshes are
  CYL/SPH (radial, quotient named in orbit-space language, group never named)
  and `pwr_slab_half_cell` (quotient realized purely as `origin=0` + a default
  reflective BC). `wigner_seitz_pin_cell` is NOT a symmetry reduction — it is
  a square→equal-area-disc *substitution*.
- False positive to not re-chase: `reflect_scan_coefficients` (`diamond.py`)
  is DD's diamond closure `α = −1`, not a geometric reflection.

Full audit + probe script: `scratch/geom_transforms_audit_spatial.md`,
`scratch/_geom_probe_spatial.py`.
