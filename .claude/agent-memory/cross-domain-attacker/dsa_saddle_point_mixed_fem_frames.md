---
name: dsa-saddle-point-mixed-fem-frames
description: DSA consistency ↔ mixed-FEM inf-sup / CFD pressure–velocity — the saddle-point frame fires on the diffusion/low-order member ONLY; consistency IS Schur-complement inheritance; no primal checkerboard in transport (born staggered)
metadata:
  type: project
---

# DSA consistency ↔ saddle-point / mixed-FEM / CFD pressure–velocity coupling

Deliverable: `.claude/plans/dsa_saddle_point_frame.md` (cited from
#312 LD-arm, #314 2-D DSA, #294 diffusion, #200 Krylov, + a new
k-skeleton design issue). Grounded on the merged DSA #2 campaign
(`docs/theory/methods/sn/acceleration.rst` R4 ruling), the D2
characterization, and the local literature sidecars. **8 hypotheses:
6 VALIDATED, 1 REFINED (Darcy not Stokes), 1 REFUTED-with-reason (no
primal checkerboard).**

## The one durable structural claim

**DSA consistency and mixed-FEM inf-sup stability are the SAME theorem.**
The transport (φ, J) = (scalar flux, current) pair is a discrete
**saddle point** `[[A,Bᵀ],[B,C]]` of **Darcy type**. "Consistent DSA"
= "the low-order operator is the **Schur complement of a compatible
discrete pairing**" — already stated verbatim in the codebase
(`acceleration.rst`: low-order = "Schur complement of a two-moment
(ℓ≤1) Galerkin triple product R₁·A_high·P₁ on the *assembled* DD
operator"). The attack's value was IMPORTABILITY (turn the stated fact
into a cross-domain-portable theorem + pollination), not discovery.

## WHERE the frame fires (the scoping ruling — ties L-007)

The saddle-point/mixed-FEM frame fires on the **diffusion/low-order
member ONLY**, nowhere else in the transport stack. Reason: the mixed
structure appears only under P1 truncation; the primary transport
operator is either the characteristic sweep `L⁻¹` (triangular, no
saddle) or the Peierls integral `I−PL⁻¹Σs` (compact-perturbation-of-
identity, no saddle). This is the **resolvent-backbone exception**
(L-007) seen from the other side: elliptic ⇒ saddle-point-when-mixed
(diffusion); characteristic ⇒ triangular-sweep (SN/MoC/CP). **Never
attack the sweep with inf-sup theory — it has no saddle to stabilize.**

## Durable rulings (per-hypothesis, condensed)

- **(a) REFINED: Darcy, not Stokes.** (J,J)-block is a MASS matrix
  `D⁻¹` (⇒ current Schur-eliminated in CLOSED FORM = Larsen step 4),
  not a vector Laplacian (Stokes A⁻¹ nonlocal, no closed elimination).
  (φ,φ)-block is removal `Σa=Σt(1−c)≥0` (coercive, →0 as c→1), NOT an
  exact-zero pure constraint. Import the **H(div)/mixed-Poisson**
  branch (RT0-P0 stable), NEVER the Stokes velocity-pressure branch
  (Taylor–Hood is inf-sup-stable AND diverges as an accelerator —
  inf-sup ≠ consistency).
- **(b) VALIDATED: Reed = staggering mismatch, not checkerboard.** Two
  individually-stable operators (DD sweep's edge-φ/cell-J implied
  diffusion vs cell-centered diffusion) whose thick-cell limits
  disagree O(1) (Alcouffe p.348 via acceleration.rst). CFD twin =
  inconsistent pressure-Poisson in a SIMPLE/fractional-step loop. D2
  Part C: cell-centered-as-accelerator ρ→54.7 (diverges σₜh≥2);
  derived edge flat ≤0.181.
- **(c) VALIDATED — deepest node: consistency = Schur inheritance.**
  Schur complement of an inf-sup-stable pairing inherits coercivity
  (β²); reduction inherits the staggering for free ⇒ derived edge op
  needs NO harmonic mean (unknowns straddle homogeneous cells) while
  independent cell op REQUIRES it. Sharp witness: **S2-GL anchor** —
  ℓ≤1 reduction is EXACT ⇒ DSA converges in ONE iteration (ρ≈3e-15).
  Hits all 4 elegance criteria. First-test the frame demands: a
  same-skeleton-wrong-coefficients negative control (route β) to
  isolate coefficient-inheritance from skeleton-home (the campaign's
  derived-vs-landed differs in BOTH).
- **(d) VALIDATED: three-tier ladder** fully-consistent (MAC ↔ derived
  edge) / partially-consistent (Rhie–Chow ↔ M4S) / inconsistent
  (collocated ↔ Reed). The D2 M-matrix **sign-flip** = "the
  accelerator must match the operator being accelerated INCLUDING its
  non-monotone defects; a better-behaved diffusion op is a WORSE
  accelerator." Dose-threshold migration (McCoy–Larsen Table II) =
  "inconsistency has no safe dose" = the stabilization-parameter story.
- **(e) VALIDATED: TWO GENUINELY DISTINCT Krylov-rescue modes.**
  (i) inconsistency-divergence = an "unstable acceleration scheme =
  preconditioner that increases the large eigenvalues" (Century p.81);
  Krylov rescues via spectral-radius-IRRELEVANCE (BAND of bad
  eigenvalues), BOUNDED ("amplification too large ⇒ still fails").
  (ii) consistent-DSA discontinuity-degradation = "only ONE eigenvalue
  near unity ruins acceleration but a single large eigenvalue barely
  affects Krylov" (Century p.83, [124]=WWM NSE 147:218); Krylov
  rescues via OUTLIER-DEFLATION (ISOLATED interface eigenvalues).
  Distinct in: what's wrong (discretization vs model), does consistency
  fix it (yes vs no), how Krylov helps (band-tolerance vs outlier-
  deflation). Discriminator = SHAPE of the near-1 spectrum. **1-D
  cannot show mode (ii)** ⇒ Krylov looks useless in 1-D (D13 "+1/+2")
  and becomes essential in 2-D-with-discontinuities (#314).
- **(f) REFUTED: no primal checkerboard in transport.** Transport is
  BORN STAGGERED (DD = the MAC grid of transport: φ on 0-skeleton, J on
  d-cells, maximally separated); removal Σa≥0 coercivizes the φ-block.
  Thick-cell DD oscillation + edge-op sign-flip are DISPERSION defects
  (non-monotone but NONSINGULAR — oscillatory≠null). LD-on-tets SAPD
  loss is a DEFINITENESS defect (indefinite, needs MINRES/GMRES not
  CG), still not a kernel. LMM(1987)/Larsen–Morel(1989) establish
  diffusion-limit ACCURACY (cell-average vs cell-edge limit = the
  staggering), NOT null-mode absence — do not over-cite them as
  stability results. **The CFD failure-mode parallel BREAKS: roles
  correspond, canonical instabilities do not.**
- **(g) VALIDATED: k-skeleton = the right DEC design axis.** Two DEC
  arrangements coexist (not a twin): **covolume** (0-skeleton φ, dual-0
  J = DSA-derived, forced by the sweep) vs **RT0/mixed-FEM** (d-cell φ,
  (d−1)-face J = diffusion module). 2-D DSA scalar is Alcouffe CORNERS
  (0-skeleton), currents cell-center (Adams–Larsen §IV.D.4 9-point).
  THREE warnings: (1) DSA f₀ EARNS a type by the PERSISTENCE criterion;
  (2) DEFER general `TypedCochainField[k]` — only 1 *persistent* client
  (face-flux is a retired transient, doesn't count); named trigger =
  #314 corner-f₀ or nodal-DG; (3) 1-D degeneracy (0-skeleton≡(d−1)-
  faces) HIDES the axis — state the home on the generic 2-D footing.
- **(h) VALIDATED: the angular axis is the same cochain complex.**
  Half-angle faces = (d−1)-cochains on the μ-mesh; redistribution
  `(1−µ²)/r ∂_µ` = angular conservation law; poles µ=±1 = no-flux BC.
  Bailey–Morel–Chang(2010): the angular differencing has its OWN
  diffusion-limit consistency — Morel–Montry weighted-diamond (the
  codebase's `τ_raw∈(0,1)` in `radial_characteristic_space.py`)
  preserves the "Galerkin ℓ≤1 diffusion approximation" to FIRST order
  (full consistency); step/diamond only to leading order (residual
  flux-dip). "Preserving the Galerkin approximation ≡ preserving the
  asymptotic diffusion limit to first order." **The k-skeleton +
  diffusion-limit-consistency machinery is AXIS-AGNOSTIC (spatial LMM
  ∥ angular BMC).** Curvilinear-DSA's missing stability theory (A&L
  p.79 "none exists") — its ANGULAR half is already implemented.

## The promotable doctrine (3rd sighting — lesson-sharpening candidate)

**Persistence, not cochain-degree, is the type-earning criterion for a
mesh field: persistent-iterate-state ⇒ typed cochain; sweep-transient
⇒ native.** Three independent applications of ONE rule: WavefrontFlux
RETIRED (transient, S6.4(f)); ψ½ TYPED (#282 persistent back-edge);
DSA f₀ TO-BE-TYPED (persistent SI correction). Sharpens L-004 (property-
vs-type) with a persistence axis specific to cochain fields. Promote to
a lesson or `coding-elegance` when a 4th sighting lands.

## Pollination (importable — with triggers)

P1 inf-sup/LBB test = fast NECESSARY screen for future (φ,J) pairs
(#314) — but NOT sufficient (Taylor–Hood counterexample). P2 **NDA**
(Hammer–Morel–Wang 2019) = the NONLINEAR sibling (≅ SIMPLE): enforces
consistency nonlinearly ⇒ converges regardless of staggering, handles
VOIDS where linear D→∞; the [124]-answer that avoids Krylov (#294).
P3 Rhie–Chow SCALING recipe for M4S LD arm (#312) — exactly-scaled
compact correction, NOT a φ-null-mode stabilizer. P4 MINRES + block-
saddle-point preconditioners for symmetric-indefinite LD-on-tets
(#200/#312). P5 covolume/MFD branch (Hyman–Shashkov) NOT RT0 for the
2-D corner op (#314). P6 consistent-projection adjoint div/grad
discipline (R=P*) for 2-D DSA (#314).

## Acquisition

HARD (frame depends on it, NOT local, do NOT fetch online — flag user):
**Warsa–Wareing–Morel NSE 147:218** (the [124] degraded-effectiveness
PRIMARY; local secondaries = Century [124] summary + WWM 147:26 carry
claims not proofs). SOFT (if pollination pursued): Boffi–Brezzi–Fortin
(inf-sup), Elman–Silvester–Wathen (saddle-point Krylov/MINRES),
Rhie–Chow 1983 (M4S scaling), Lipnikov–Manzini–Shashkov (MFD/covolume).

## UNEXPLORED (dead frames + reasons)

homology (no ∂²=0; trace±extension = dagger pair not differential) ·
category/double-cat (biproduct already captures it, L-011) · MPO (fixed
2×2, bond-dim degenerate) · multigrid (NATIVE not foreign — DSA IS the
Galerkin coarse op) · Riccati (static KKT, no temporal recursion) ·
symplectic (stationarity not Hamiltonian flow) · Wiener–Hopf (wrong
family) · Schwarz-DD (orthogonal — fires on parallel-sweep impl not
consistency) · Christoffel (fires on the (1−µ²)/r redistribution's
conservation-law face in 4h, not on DSA consistency).
