---
name: cylindrical-sn-level-order-sensitivity
description: #326 cylindrical SN — why the per-level ordinate ORDER is load-bearing (the sweep reads only (eta,w); a mirror-pair tie splits tau into {1,0.5}), AND the half-range fix map: fold belongs in LevelStructure not the measure, (A) fold-existing vs (B) Hebert-midpoint split by the R12a predicate, and the xi-odd SH moment is the one real break.
metadata:
  type: reference
---

Durable numerical facts behind GitHub issue #326 (per-level ordinate ordering in
the cylindrical SN sweep). Measured 2026-08-01 on `refactor/operator-strategy-layers`.
The SHAPE below survives churn; re-derive `file:line` via Nexus.

## The sweep's angular alphabet is (eta, w) only

The 1-D cylindrical sweep/matvec never reads `xi = mu_y` (`Quadrature.xi`,
`nodes[:,1]`, the azimuthal cosine `Omega . phi-hat`). Grep-verified: `mu_y`/`xi`
appear NOWHERE in `orpheus/transport/` and only in docstrings in `orpheus/sn/`
and `orpheus/geometry/`. `mu_z` is read at exactly one place — `sin_theta` from
`mu_z[level_idx[0]]`, a LEVEL LABEL. `_streaming_axes` is built over
`range(ndim)`, so a 1-D mesh only ever touches axis 0.

Consequence: two ordinates of a **xi-mirror pair** (equal eta, opposite xi —
`quad.reflection_index("y")` is their map) are numerically INDISTINGUISHABLE to
the sweep. They differ only through (a) the per-ordinate source `Q_n` and
(b) per-ordinate inflow/BC data.

## The tie splits tau into {1, 0.5} — decided by 1 ULP of cos(phi)

Within a level, `eta_edge[m+1] = 0.5*(eta[m]+eta[m+1])`. For a mirror pair at
positions (m, m+1) the midpoint edge COLLAPSES onto the node, so
`tau_raw[m] = 1` and `tau_raw[m+1] ~ 1e-16 -> clamped to 0.5`. Which physical
ordinate lands on which side is chosen by `np.argsort` over eta values that are
**1 ULP apart** (trig round-off in `cos(phi)`) — today deterministic and
kind-independent, but only by accident of that round-off. Under #325's
algebraically-exact nodes the ties become EXACT and the (non-stable) argsort
tie-break becomes the decider.

## Measured blast radius of a reorder

Flipping the within-mirror-pair order (cylinder, `product(n_mu=2, n_phi=8)`):

- `alpha_per_level` — **bit-identical** (partial sum of `w*eta`; w constant
  within a level for `product`, eta equal ⇒ every partial sum unchanged).
- `tau` / `c_in` / `c_out` **per ordinate** — PERMUTE (position-stable, but the
  ordinate→position assignment moves).
- Per-ordinate angular flux — an exact PERMUTATION (`max|perm(base)-flip| = 9e-16`),
  but at a FIXED ordinate index it moves by **12.7%**.
- Scalar flux, **isotropic** source — 1–2 ULP (the permutation is invisible to
  `sum_n w_n psi_n`).
- Scalar flux, **xi-dependent** source (`q_n = 1 + xi_n`) — moves by **~20%**.
- The current ordering BREAKS the xi -> -xi mirror symmetry of the discrete
  solution: mirror-pair angular fluxes differ by 2e-3 to 3.7e-3 relative on an
  isotropic source, where the continuous problem demands equality.

## Neighbouring facts worth keeping

- `level_symmetric` cylinder levels are **4-to-1** over eta (LS-4 level 0: 16
  ordinates, 4 distinct eta), so the ordering question is wider there, and its
  `tau_raw[0] = 1` exactly (vs `0` exactly for `product`).
- Cylinder + `product` carries **NO** psi-half seed block:
  `SNMesh.radial_characteristic_levels == ()` because the R12a predicate is
  `0 < tau_raw[0] < 1` and `tau_raw[0] = 0` bit-exactly. The sphere-GL is the
  only carrying instance. So every `radial_characteristic*` seed path is
  sphere-only in practice.
- The degenerate pure-azimuthal ordinates (`|eta| < 1e-15`, `phi = pi/2, 3pi/2`)
  are on NO walk leg and take a volumetric per-cell path.

## The literature ruling: what alpha IS, and why the level must be a HALF range

(Added 2026-08-01 by a second explorer dispatch — the theory-corpus/literature half
of the same issue. Sources spot-verified against the rendered scans.)

**alpha is the TANGENTIAL direction cosine at a half-angle boundary, times the level
weight.** Hebert (2009) Eq. (3.399) literally *defines* it that way
(`alpha_{p,q+-1/2} = W_p * eta^Hebert_{p,q+-1/2}`, where Hebert's `eta` = ORPHEUS `xi`
= `mu_y`), derived from flat-flux preservation in Eqs. (3.397)-(3.398). In ORPHEUS
spelling `alpha_{m+1/2} = -W_p * xi(phi_{m+1/2})`, exactly, with the seed
`alpha_{1/2}=0` being nothing but "the starting boundary is the plane xi=0"
(`phi = pi`). For the equispaced product rule the Dirichlet kernel makes it EXACT:
`alpha_k = -w_gl * kappa * [xi(omega_{k-1/2}) - xi(omega_{-1/2})]`,
`kappa = d_omega/(2 sin(d_omega/2)) -> 1`; nodes sit at the CENTRE of their angular
cell, so boundaries are at `omega_j +- d_omega/2`.

**The recursion is a cumulative integral in `omega`, so the enumeration must be
monotone in `omega` over a CONTIGUOUS interval.** Sorting by `eta = sin(theta)cos(phi)`
is equivalent only on a HALF circle, where `eta` is monotone. Both primary sources put
a level on a half range: Hebert Sect. 3.9.3 ("two octants ... base points in interval
`0 <= omega <= pi`", x4 in the moment sum Eq. 3.401) and Bailey-Morel-Chang (2010)
Eq. (52) (level weights "normalized on each level to sum to `2*sqrt(1-xi_n^2)`" — i.e.
the radial-cosine range covered EXACTLY ONCE). ORPHEUS's `product_mu_phi` puts the
FULL circle `[0, 2pi)` in one level, so it covers it TWICE and the eta-sort
interleaves the two halves. The `tau_raw in {0,1}` trichotomy the codebase documents
as structural is the direct FINGERPRINT of that doubled covering (zero-width angular
cells). ORPHEUS's alpha read at the mirror-PAIR boundaries is bit-exactly the correct
half-range alpha; the intra-pair entries correspond to no real angular boundary.

**Citation hazard, live at HEAD:** `orpheus/numerics/quadrature/rules_product.py:40-44`
and `:66-68` cite "Bailey, Adams, Yang & Zika (2009) ... Eq. 50" for the eta-ordering
convention. `curvilinear_one_group.rst` sect. "Citation correction" already retracts that
paper (wrong Bailey — a PLD-FE diffusion paper) but its fix list predates this module.
And BMC (2010) Eq. (50) is the two-sided CLOSURE `alpha_{1/2,n} = alpha_{M+1/2,n} = 0`
— it says nothing about ordering.

**Every existing alpha gate is structurally blind to the ordering.** Closure is a
telescoping sum (permutation-invariant); the flat-flux proof uses only the recursion
STEP; `tests/.../test_pole_angular_closure.py::TestAlphaRecursionIdentities` builds its
fixtures from `Quadrature.gauss_legendre` only (never the cylindrical `product` rule);
`test_alpha_per_level_bit_identical` is a frozen-RHS route-equivalence. And the
curvilinear MMS is Mode-7 blind (both ansatzes are functions of `eta` and `xi**2`, i.e.
inside the symmetric sector) — so the `psi(xi) == psi(-xi)` mirror invariant is the
adjudicator, not the MMS.

**Reachability bound:** odd `n_phi` is REJECTED at `SNMesh` construction for the
cylinder (`BoundaryGeometryMapNotMeasurePreservingError`, ERR-042 — an odd equispaced
`phi` grid is not closed under the x-mirror, and `r=0` is intrinsically specular). So
`n_phi` is effectively EVEN-only, every level has an exact 2-fold `xi`-mirror pairing
with two self-mirror `xi=0` members (`phi = 0, pi`), and the tied-`argmin`
starting-direction case is LATENT, not live.

## The HALF-RANGE fix: what it touches (measured 2026-08-01, second dispatch)

The constructive exit for #326 is the Hebert half-range level. Two structural
rulings that survive churn:

**1. The fold belongs in `LevelStructure`, NEVER in the `DiscreteMeasure`.**
`product_mu_phi` returns `(DiscreteMeasure, LevelStructure)` and the two halves
have DISJOINT consumers. The measure is a full-S2 cubature consumed by 2-D
Cartesian meshes, by `symmetry.py`'s group-invariance checks, and by the registry
(which types it `SO2` with a declared `degree_of_exactness`). `LevelStructure` is
read by NOTHING outside the curvilinear path. Geometrically: the measure lives on
the COVER S2, `LevelStructure` describes the QUOTIENT S2/Z2 (the r-z mirror) —
so putting the fold in the level structure puts the quotient where it already lives.

**2. TWO half-range constructions exist and the R12a predicate separates them.**
Both reproduce the alpha closed form with the SAME constant `2*kappa` (2.052344 at
n_phi=8), so ALPHA DOES NOT DISCRIMINATE THEM. What does:

- **(A) fold the existing nodes** onto `omega in [0,pi]`, endpoints INCLUDED,
  trapezoid weights `[w, 2w, ..., 2w, w]` -> `tau_raw = [0, .293, .5, .707, 1]`.
  `tau_raw[0] = 0` is PRESERVED, so R12a still says non-carrying: the seed
  trichotomy, the `#280 2.5b` direct-seed fold and the sphere-only
  `radial_characteristic_*` story all survive untouched.
- **(B) Hebert midpoint half-range** (`2F(p)` nodes strictly inside `(0,pi)`) ->
  `tau_raw = [.2195, .4142, .5858, .7805]`, all strictly interior, so
  **EVERY cylindrical level becomes R12a-CARRYING** — the whole psi-half machinery
  goes live on the cylinder and the 2.5b fold retires. Also moves the nodes
  (mirrored back, it is the MIDPOINT circle rule, offset by d_phi/2).

(A) is what the tree's own `test_alpha_closed_form.py::_half_range_level` uses.

**3. The ONE thing that legitimately needs the full circle: xi-ODD SH moments.**
A naive "halve the nodes, double the weights" reproduces every xi-EVEN moment to
5e-16 and turns `phi_1^xi` from its structural `-1.3e-16` into `+2.94`. Live
consumer: cylindrical P1 scattering (`test_curvilinear_aniso_scattering_p1.py`
sums all three `Y_1^m` slots). Correct framing: the fold RESTRICTS to the xi-even
sector, so the harmonic TRIAL SPACE must restrict to its xi-even subspace — the
odd harmonics are out of the space, not miscomputed.

**4. The (a1)/(a2) split is the highest-leverage question.** Fold the SWEEP only
(lift the mirror partners back into the full (N,ng,nx) buffer; state stays N-dim)
=> the moment break VANISHES (moments are taken on the lifted, exactly-even psi),
and nothing in the frame / octants / trace-space / Krylov surface moves. Fold the
STATE (the 2x-fewer-unknowns win, 1.6x at n_phi=8, 1.88x at 32) => the xi-even
trial-space restriction, `octants` 4-of-8, and an `n_dof`/`restart` resize
(ERR-053 family) all become due. Recommend (a1)+(A) first; (a2) is a perf follow-up.

**5. Other measured facts.** The redistribution term `(dA/w)(a_out psi_out -
a_in psi_in)` is INVARIANT under a uniform level-weight rescale (alpha scales with
w, dA/w as 1/w) — so doubled weights are self-consistent in the sweep.
`reflection_index("x")` CLOSES the `xi>=0` half for every n_phi, and the flagged
pole-map ambiguity `(-eta,+xi)` vs `(-eta,-xi)` (56/64 ordinates differ at
product(4,16); the 8 agreeing are the two xi=0 self-mirror members per level)
DISSOLVES under the fold. `net_current` reads the `|Omega.n|*w = |eta|*w` metric —
xi-even, unaffected. The BMC cumulative-weight edge convention gives `tau_raw`
OUTSIDE [0,1] on a cosine grid, so the cylinder's eta-midpoint edges are correct
and are NOT a second latent bug.

**6. `level_symmetric` on a cylinder.** `select_quadrature("cylinder", ...)`
REJECTS it ("geometry SO2 is not a subgroup of rule's invariance group Oh") — but
the selector has NO production caller (opt-in; only tests call it), and
`cylindrical_streaming` accepts any level-structured quadrature, so LS-on-cylinder
is reachable by direct `SNMesh(...)` and is pinned by ~6 test modules. Its level is
4-to-1 (both signs of mu_z AND xi); both redundancies are physical for an infinite
1-D cylinder. Folding it 4-to-1 makes EVERY LS level R12a-carrying (there is no LS
analogue of "the omega=pi endpoint is already a node"), so LS needs its own ruling:
fold-product-only-and-keep-LS-red, refuse-LS-on-cylinder, or accept carrying LS.

**7. DOC FALSIFICATION (L33).** `curvilinear_one_group.rst`
§`sn-direct-seed-circle-vs-interval` states the principle *"a periodic
redistribution axis gives edge-inclusion for free; an interval axis makes you pay
for it with a separate seed"* and attributes the cylinder's `tau_raw,0 = 0` to
PERIODICITY. FALSE: the alpha recursion is a cumulative integral on the QUOTIENT
INTERVAL `omega in [0,pi]`; the eta-ascending march visits `phi/pi = 1.00, 1.25,
0.75, 1.50, 0.50, 1.75, 0.25, 0.00` (two half-circles interleaved) and the alpha
dome carries a PLATEAU. The free edge-inclusion comes from the grid standing on the
QUOTIENT's endpoint (`omega = pi`, where `xi = 0` so the redistributed flux
vanishes), not from periodicity — the doc's own "contingent on even n_phi" escape
hatch gives it away. The section's CONCLUSION survives under (A), is falsified
under (B); its REASON is wrong either way, and the Gauss-Lobatto study that follows
inherits the correction. Also: the theory corpus has ZERO mention of the
psi-even-in-xi theorem or of #326 — it lives only in a test docstring — and
`rules_product.py:40-44` still carries the RETRACTED Bailey-Adams-Yang-Zika 2009
citation (the Phase-B fix list in `curvilinear_one_group.rst` predates the module).
