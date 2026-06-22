---
name: issue-240-d5b-s3-a0-spatial-moment-space
description: #240 D5b-S3-A0 elegance review — mint SpatialMomentSpace (typed spatial-moment peer of SphericalHarmonicSpace) + optional field-space factor (default-OFF) + find_factor (#207). PASS-WITH-NITS; 1 must-fix (lost type on find_factor), 1 do-now twin (re-spelled compose gate).
metadata:
  type: project
---

# #240 D5b-S3-A0 — SpatialMomentSpace + typed-field-space widening (PASS-WITH-NITS)

Branch `feature/sn-space-angle-tier2`. Reviewed pre-commit. A foundational typed-field-space
change minting the SPATIAL within-cell tensor-Legendre DG moment axis ({ψ̄,ψ̂_x,ψ̂_y,ψ̂_xy}) as a
first-class typed space (user chose: typed factor, NOT raw int axis), the ANGULAR sibling of
`HarmonicMomentField`'s `SphericalHarmonicSpace`. Construct-general/select-narrow: capability
default-OFF, no production field carries the axis yet (consumer = S3-A).

## VERDICT PASS-WITH-NITS. 2 approval conditions + 1 issue.

### Q2 MUST-FIX (VIOLATION, the headline) — `find_factor` LOSES the queried type.
`numerics/space.py:396` `def find_factor(self, factor_type: type) -> "FunctionSpace"`. Returns the
BASE, so `find_factor(SphericalHarmonicSpace).L` resolves to `object`/error (pyright
`harmonic_moment_field.py:228`/`space.py:521`). The method's whole reason to exist is the typed
bridge composed-space→factor-metadata (`.L`/`.per_axis`); dropping the type defeats it (anti-#13
one rung up). Fix = PEP-695 generic (project is Python 3.14, native, no TypeVar import):
`def find_factor[T: FunctionSpace](self, factor_type: type[T]) -> T:`. Closes #207 WITH the type
not just the value. Zero call-site churn (the one internal consumer `_spatial_moment_tail_of`
reads `.shape`, on the base). RECURRING TELL: a typed-query helper that returns the base instead
of the generic — always demand `type[T]->T`.

### Q4 DO-NOW (CONCERN twin) — compose gate re-spelled in `from_mesh_and_L`.
`_bases.py:194 BulkField._compose_spatial_moments` is the SSOT "append SpatialMomentSpace iff
tail != ()". Angular/Scalar delegate; `harmonic_moment_field.py:266-271` RE-INLINES it (negated
`!= ()` vs helper's `== (): return` — spelling inversion hides the twinhood). HMF IS-A BulkField;
the helper appends to ANY space (build `sh*cell_group` first, then call the staticmethod on it).
"Uniform across 2 of 3 carriers, hand-edited on the 3rd." Collapse: `space =
cls._compose_spatial_moments(sh_space*cell_group_space, mesh, spatial_moments)`. Habitat = a future
gate change (metric on factor / oracle length-1 axis / the Q3 relocation) lands on the helper, the
HMF twin keeps the old behaviour (Phase-F shape).

### Q3 CONCERN-not-blocker (the layering verdict you wanted) — deferred numerics→sn import is RIGHT, but the policy is mis-homed.
The cycle is created by the PACKAGE `orpheus.sn.__init__`→transport→harmonic_moment_field→back,
NOT by `_ubld` (numpy-only leaf, zero SN imports). Deferring to call-time honors "numerics must
not import the SN package at module-load" — CORRECT, better than duplicating the constant. BUT
`AVERAGE_MOMENT=0` + `face_moment_tail` ("append iff >1") are PURE storage-layout conventions
(no σ/quad/transport) — physics-free is the TELL they are numerics concepts mis-homed in
`sn.spatial._ubld`; the deferred import works BECAUSE they're already physics-free. Fully-elegant
end-state: home the layout policy in `numerics`, `_ubld` imports DOWN, deferred import retired.
NOT now (Pattern 6 — the 2nd numerics consumer, the S3-A iterate seam, hasn't landed). FILE a
`module:sn`/`type:improvement` issue. Latent inversion: a 3rd deferred-importer cements the
up-dep at N sites.

### PASS axes (reinforce):
- Q1 true peer: same frozen-dataclass/`from_*`/size-in-shape/`__post_init__` Pattern-4 check/
  `__eq__`-delegation-over-ndarray-hazard idioms as `SphericalHarmonicSpace`; Euclidean
  `inner_product_weights=None` is the HONEST divergence (cell-mass metric on the UBLD op, #207).
- Q5 MomentDisplacement ripple = the affine-torsor symmetry WORKING (inst-knowledge #6): the
  `FluxRole._mint_displacement` field-copy `TypeError` was the type system ENFORCING that φ⊖φ
  lands in the flux's widened space. "Will every future field ripple?" — YES, that's the
  guardrail; a displacement NOT tracking a new flux field admits the wrong-space difference.
  Angular/Scalar Displacement needed no change (their fluxes carry no `spatial_moments` FIELD,
  read the tail off the shared space) — mechanism precise not blunt. The field-vs-space-factor
  distinction is the discriminator: a new dataclass FIELD on a flux ripples to its Displacement;
  a new space FACTOR does not (read off `find_factor`).
- Q6 einsum lift textbook Pattern-7: `"fc->gc"`⇒`"fc...->gc..."` (×3: :503/:533/:599) = Σ_s⊗I
  spectator broadcast at the PRODUCER (`MaterialXSField`, not each consumer); byte-id at no-axis
  (`...` matches nothing). 3 near-verbatim docstrings = acceptable prose dup (not logic).

Gates spot-checked: 34 new foundation green -O; closeout's 513/1/4 DD bit-id negative control +
mutation-verified [ld] byte-id teeth (auto-read scheme → 2 reds) credible.
