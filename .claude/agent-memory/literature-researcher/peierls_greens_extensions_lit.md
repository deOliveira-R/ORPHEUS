---
name: Literature for Variant α Green's function extensions across geometries / topologies / BC / MG / multi-region
description: Reference list for expanding the Plan 2 Variant α Green's function approach across multiple fronts. Organised by extension dimension. Indicates which references are locally available in scratch/literature/ vs still need pulling.
type: reference
---

# Literature for Variant α Green's function extensions

**Predecessors**:

- B1 literature memo: `peierls_greens_function_lit.md`
- B2 decision: `../numerics-investigator/peierls_greens_variant_alpha_decision.md`
- B6 closeout: `../numerics-investigator/peierls_greens_phase1_closeout.md`

**Local literature directory**: `/workspaces/ORPHEUS/scratch/literature/`
(actively expanding — user-curated).

**Marker convention**:

- ✓ **LOCAL** — PDF in `scratch/literature/`
- ◇ **WISH** — needs to be pulled (institutional access or DOI)
- — — book / textbook (not a single PDF; full reference)

## Already locally available (as of 2026-05-01)

The user has populated `scratch/literature/` with 22 PDFs covering
most foundational references. These directly enable several
extensions without further library pulls:

| Reference                                                                | File                                                                                        |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| Sanchez 1986 TTSP 14 (sphere, anisotropic, α/β BC)                       | `Sanchez(1986)Integral form of the equation of transfer...`                                 |
| **Pomraning-Siewert 1982 JQSRT 28 (vacuum sphere, isotropic precursor)** | `Pomraning Siewert (1982) On the integral form of the equation of transfer...`              |
| Sanchez 2002 NSE 140 (lattice periodic-trajectory BC)                    | `Sanchez(2002) Treatment of Boundary Conditions in Trajectory-Based Deterministic...`       |
| Sanchez & McCormick 1982 NSE 80 (review)                                 | `Sanchez(1982) A review of neutron transport approximations`                                |
| Sanchez 1977 (2-D heterogeneous CP)                                      | `Sanchez(1977)Approximate Solutions of the Two-Dimensional Integral Transport Equation...` |
| Carlvik 1967 (finite cylinder + cuboid CP)                               | `Calvik(1967)Collision Probabilities for Finite Cylinders and Cuboids` (sic: "Calvik")      |
| Mitsis 1963 ANL-6787 (sphere Case-Zweifel multi-region)                  | `Mitsis(1963) Transport Solutions to the Monoenergetic Critical Problems`                   |
| Brockmann 1981 NSE (multi-layer Green's functions, anisotropic)          | `Brockmann(1981)Treatment of Anisotropic Scattering...`                                     |
| Roy 1989 NSE (DRAGON foundational, 3-D heterogeneous lattices)           | `Roy(1989)A Transport Method for Treating Three-Dimensional Lattices...`                    |
| Halsall 1980 (CACTUS — UK MoC)                                           | `Halsall - 1980 - CACTUS...`                                                                |
| Askew 1972 (characteristics formulation)                                 | `Askew - 1972 - A Charactheristics Formulation...` (sic: "Charactheristics")                |
| Mazumdar 2015 (modern MoC)                                               | `Mazumdar(2015)Solution of neutron transport equation by method of characteristics`         |
| Morel 1989 (hybrid collocation-Galerkin-S_N)                             | `Morel(1989)A Hybrid Collocation-Galerkin-Sn Method...`                                     |
| Benoist 1981 (Wigner-Seitz diffusion-via-IT)                             | `Benoist(1981)Integral Transport Theory Formalism for Diffusion Coefficient...`             |
| Pomraning 1989 (transport in general geometry)                           | `Pomraning(1989)The Transport Equation in General Geometry`                                 |
| Hébert 2009 Chapter 3                                                    | `Hebert(2009)Chapter3` (CP basics + §3.8.5 white BC)                                        |
| Ligou 1982 Chapter 8                                                     | `Ligou(1982)Chapter8`                                                                       |
| Stacey 2007 Chapter 9                                                    | `Stacey(2007)Chapter9`                                                                      |
| Stammler 1983 Chapter 4 (cylinder CP)                                    | `Stammler(1983)Chapter4`                                                                    |
| Stammler 1983 Chapter 6 (sphere CP)                                      | `Stammler(1983)Chapter6`                                                                    |
| Atkinson — *Numerical Solution of Integral Equations of the Second Kind* | `Atkinson - The numerical Solution of Integral Equations...`                                |

## Extension by extension

### 1. Different geometries

#### Cylinder (1-D and 3-D) — Variant α extension

- ✓ **LOCAL** Carlvik 1967 — finite cylinder + cuboid CP, Bickley-Naylor primitives.
- ✓ **LOCAL** Stammler 1983 Chapter 4 — cylinder CP with reflective BC.
- ✓ **LOCAL** Sanchez 1977 NSE — 2-D heterogeneous integral transport (cylinder lattice in cross-section).
- ◇ **WISH** Knyazev, B. A. & Selivanov, A. N. (2014). *Bickley-Naylor functions and their integrals*. Nauka, Moscow / English translation. Modern 3-D Bickley + Ki integral tables. **Critical for cylinder Variant α** (the 3-D Ki correction Knott already needs in ORPHEUS Issue #132 cylinder Hébert).
- ◇ **WISH** Modak, R. S. & Gupta, A. (2003). "Use of the boundary element method in the integral transport equation". *Annals of Nuclear Energy* 30, 943-959. DOI `10.1016/S0306-4549(03)00009-4`. Cylinder IE with reflective BC, BEM formulation.

#### 2-D Cartesian / Hexagonal / Triangular

- ✓ **LOCAL** Sanchez 1977 NSE — directly the 2-D foundational paper.
- ✓ **LOCAL** Halsall 1980 (CACTUS) — 2-D MoC for complicated geometries; relevant to Variant α generalisation to 2-D bouncing trajectories.
- ✓ **LOCAL** Askew 1972 — characteristics formulation (2-D MoC predecessor to CACTUS).
- ✓ **LOCAL** Roy 1989 NSE (DRAGON) — 3-D heterogeneous lattices via transport.
- ◇ **WISH** Sanchez, R. & Mao, L. (1993). "Linear multiple collision probability methods for 2-D heterogeneous geometries". *Annals of Nuclear Energy* 20, 481. **Modern 2-D-CP — directly extends Sanchez 1977**.
- ◇ **WISH** Suslov, I. R. (1997). "Solution of transport equation in 2- and 3-dimensional irregular geometry by the method of characteristics". *Proc. M&C 97 Saratoga Springs*. Modern arbitrary-2-D MoC.

#### Slab refinements + benchmarks

- — Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*, Chapter 1 §1.6-1.7 + Chapter 5 (textbook).
- ✓ **LOCAL** Pomraning 1989 (general geometry transport) — covers slab integral form as a special case.

#### Sphere (Variant α specific extensions)

- ✓ **LOCAL** Sanchez 1986 TTSP 14 — full α/β BC + linearly anisotropic.
- ✓ **LOCAL** **Pomraning-Siewert 1982 JQSRT 28** — **vacuum sphere isotropic precursor**; structurally-independent reference for Variant α vacuum BC test (closes V_α3 numerical gap).
- ✓ **LOCAL** Mitsis 1963 ANL-6787 — sphere Case-Zweifel singular eigenfunction expansion + multi-region sphere.
- ✓ **LOCAL** Brockmann 1981 NSE — multi-layer Green's functions for anisotropic scattering. **Foundational for multi-region sphere Variant α**.
- ✓ **LOCAL** Garcia, R. D. M. (2017). "A P_N particular solution for the radiative transfer equation in spherical geometry". `Garcia(2017)A PN particular solution...`
- ✓ **LOCAL** Garcia, R. D. M. (2019). "A numerically stable spherical harmonics solution for the neutron transport equation in a spherical shell". `Garcia(2019)A numerically stable spherical harmonics solution...spherical shell.`
- ✓ **LOCAL** Garcia, R. D. M. (2019). "On the P_N method in spherical geometry: A stable solution for the exterior of a sphere". `Garcia(2019)On the PN method in spherical geometry a stable solution for the exterior of a sphere`. (cited as 2018 JCTT 47 in the B1 memo — actual J. Comp. Theor. Transport publication date 2019).
- ✓ **LOCAL** Garcia, R. D. M. (2021). "Accurate spherical harmonics solutions for neutron transport problems in multi-region spherical geometry". `Garcia(2021)Accurate spherical harmonics solutions...multi-region`. **Direct attack on Issue #132 Class B MR — structurally independent reference for plan (b) multi-region cross-check**.
- ◇ **WISH** Williams, M. M. R. (1983). "Transport in a multi-region sphere". *Annals of Nuclear Energy* 10, 235.

### 2. Different topologies

#### Periodic / lattice trajectories

- ✓ **LOCAL** Sanchez 2002 NSE 140 — periodic-trajectory closure for 2-D Cartesian / hexagonal lattices.
- ✓ **LOCAL** Halsall 1980 CACTUS — lattice tracking in 2-D.
- ✓ **LOCAL** Askew 1972 — characteristics-based BC handling.
- ◇ **WISH** Wu, G. J. & Roy, R. (2003). "A new characteristics algorithm for 3-D transport calculations". *Annals of Nuclear Energy* 30, 1-16. DRAGON cyclic-tracking machinery.
- ◇ **WISH** Le Tellier, R. & Hébert, A. (2006). "On the integration scheme along a trajectory for the method of characteristics". *Annals of Nuclear Energy* 33, 1260. Cyclic-tracking convergence.
- ◇ **WISH** Cho, N. Z. & Hong, S. G. (1998). "CRX: a code for rectangular and hexagonal lattices". *Annals of Nuclear Energy* 25, 547. Multi-cell periodic.
- ◇ **WISH** Postma, T. & Hennart, J. (1989). "Numerical solution of the integral neutron transport equation in 1-D periodic lattices". *AnNE* 16, 451. **Periodic CP — relevant for Variant α lattice extension**.

#### Multi-pin / heterogeneous (full assembly)

- ✓ **LOCAL** Roy 1989 NSE — DRAGON foundational, 3-D heterogeneous lattices.
- ✓ **LOCAL** Sanchez 1977 NSE — 2-D heterogeneous CP for thermal lattices.
- ✓ **LOCAL** Benoist 1981 — Wigner-Seitz cell formalism for diffusion via integral transport; useful for multi-pin homogenisation.
- ◇ **WISH** Knott, D. & Edenius, M. (1986). "CASMO theory manual". Studsvik (proprietary; not on OpenAlex). Industry-standard heterogeneous lattice methodology.
- ◇ **WISH** Marleau, G. (2001). "DRAGON theory manual". IGE-236, Institut de Génie Nucléaire. (DRAGON code documentation; freely available via Polytechnique Montréal.)

#### Hollow / annular cells

- ✓ **LOCAL** Stammler 1983 Chapters 4 (cylinder) + 6 (sphere) — hollow-cell formulations.
- ◇ **WISH** Carlvik, T. (1965). "A method for calculating collision probabilities in general cylindrical geometry". *Proc. 3rd Geneva Conf.*. Annular cylinder CP.

### 3. Multi-group

#### Foundational

- — Bell & Glasstone 1970, Chapter 7 (textbook).
- — Lewis, E. E. & Miller, W. F. (1984). *Computational Methods of Neutron Transport*. ANS, Chapters 6-7. **Standard reference**.
- ✓ **LOCAL** Stacey 2007 Chapter 9 — multi-group diffusion (companion to integral transport MG).
- ✓ **LOCAL** Hébert 2009 Chapter 3 (CP basics — full book Chapters 6-9 needed for MG-CP / resonance).

#### MG integral form / resonance

- ◇ **WISH** Hébert 2009 *Applied Reactor Physics* Chapters 6 (CP), 8 (resonance self-shielding), 9 (MoC). Hébert Ch.3 already on hand; the full book extends to MG.
- ✓ **LOCAL** Stammler 1983 Chapters 4 + 6 — MG resonance integration with CP at chapter level.
- ◇ **WISH** Gibson, N. A. & Forget, B. (2014). "Resonance treatment using pin-based pointwise energy slowing-down method". *J. Nucl. Sci. Tech.* 51. Subgroup-style MG for Variant α extension.
- ◇ **WISH** Hébert, A. (2017). "Mathematics of Nyström and Galerkin discretizations of integral transport". Survey article.

#### MG iteration acceleration

- ◇ **WISH** Larsen, E. W. (1986). "Diffusion-synthetic acceleration methods for the discrete-ordinates equations". *Transport Theory and Statistical Physics* 13, 107-126.
- ◇ **WISH** Adams, M. L. & Larsen, E. W. (2002). "Fast iterative methods for discrete-ordinates particle transport calculations". *Progress in Nuclear Energy* 40, 3-159. **Comprehensive MG acceleration review**.

### 4. Multi-region (analytic + numerical)

#### Sphere multi-region (high priority — Issue #132)

- ✓ **LOCAL** Mitsis 1963 — multi-layer sphere transport via Case-Zweifel singular eigenfunctions.
- ✓ **LOCAL** Brockmann 1981 — multi-layer Green's functions with anisotropic scattering. **Direct relevance: extends Variant α to multi-region**.
- ◇ **WISH** Garcia 2021 *J. Comp. Phys.* — modern P_N multi-region sphere with reflection.
- ◇ **WISH** Williams 1983 *AnNE* — transport in multi-region sphere.
- ◇ **WISH** Pomraning, G. C. (1991). *Linear Kinetic Theory and Particle Transport in Stochastic Mixtures*. World Scientific. Stochastic-mixture chord-length distributions transferable to deterministic multi-region.

#### Cylinder multi-region

- ✓ **LOCAL** Sanchez 1977 NSE — 2-D heterogeneous integral transport.
- ◇ **WISH** Hébert, A. (2003). "A consistent technique for the global homogenization of a pressurized water reactor assembly". *Nuclear Science and Engineering* 113, 227. Multi-region CP with white BC.

#### Hébert / Marshak multi-rank

- ✓ **LOCAL** Sanchez & McCormick 1982 NSE 80 — review of all closure approximations.
- ✓ **LOCAL** Hébert 2009 Chapter 3 — Hébert §3.8.5 (rank-1 white BC), §3.8.6+ (rank-N).
- ◇ **WISH** Hébert, A. & Saygin, H. (1992). "Development of DRAGON". *Reactor Physics & Reactor Computations*. Mark/Marshak rank-N closure development.

### 5. Different implementations / numerical methods

#### Integral equations (general theory)

- ✓ **LOCAL** Atkinson 1997 — *The Numerical Solution of Integral Equations of the Second Kind*. Cambridge. **THE reference for Nyström / Galerkin / collocation convergence theory**. Directly applicable to Variant α convergence proofs.
- ◇ **WISH** Hackbusch, W. (1995). *Integral Equations: Theory and Numerical Treatment*. Birkhäuser. Better than Atkinson on hypersingular kernels (relevant to Phase 5 retreat post-mortem).
- ◇ **WISH** Kress, R. (2014). *Linear Integral Equations* (3rd ed.). Springer Applied Math. Sciences 82. Excellent on Fredholm theory + spectral compactness.

#### Quadrature (especially µ → 0 singularity)

- ◇ **WISH** Davis, P. J. & Rabinowitz, P. (1984). *Methods of Numerical Integration* (2nd ed.). Academic Press. **Classical reference** — Gauss-Jacobi, weighted GL, singular-integrand recipes.
- ◇ **WISH** Helsing, J. (2009). "Integral equation methods for elliptic problems with boundary conditions of mixed type". *J. Comp. Phys.* 228. Modern singularity subtraction.
- ◇ **WISH** Bremer, J. (2012). "On the Nyström discretization of integral equations". *J. Sci. Comput.* 51, 51-91. **Modern Nyström for nearly-singular kernels**.

#### Power iteration / eigenvalue acceleration

- ◇ **WISH** Saad, Y. (2011). *Numerical Methods for Large Eigenvalue Problems* (revised ed.). SIAM.
- ◇ **WISH** Lehoucq, R. B., Sorensen, D. C. & Yang, C. (1998). *ARPACK Users Guide*. SIAM. Implicitly-restarted Arnoldi for large sparse eigenvalues.
- ◇ **WISH** Warsa, J. S., Wareing, T. A. & Morel, J. E. (2004). "Krylov iterative methods and the degraded effectiveness of diffusion synthetic acceleration". *Nuclear Science and Engineering* 147, 218-248. **Krylov for transport — directly applicable to Variant α acceleration**.
- ◇ **WISH** Hamilton, S. P. & Evans, T. M. (2015). "Efficient solution of the simplified P_N equations". *J. Comp. Phys.* 284. Modern eigensolvers for transport.

#### Method of Characteristics (MoC) for Variant α generalisation

- ✓ **LOCAL** Halsall 1980, Askew 1972, Mazumdar 2015 — MoC foundational + modern.
- ✓ **LOCAL** Morel 1989 — hybrid collocation-Galerkin-S_N (alternative discretisation idea).
- ◇ **WISH** Boyd, W., Forget, B. & Smith, K. (2014). "OpenMOC: A high-performance, open-source method of characteristics code". *Annals of Nuclear Energy* 68. Open-source production MoC.
- ◇ **WISH** Yamamoto, A. (2007). "Generalized coarse-mesh rebalance method for acceleration of MoC". *Nuclear Science and Engineering* 156. Acceleration techniques.

### 6. Different boundary conditions

#### Vacuum (the V&V plan (a) priority)

- ✓ **LOCAL** **Pomraning-Siewert 1982** — vacuum sphere isotropic precursor. **THE structurally-independent reference for Variant α vacuum BC validation**.
- — Bell & Glasstone 1970, Chapter 1 §1.6-1.7.
- ◇ **WISH** Davison, B. (1957). *Neutron Transport Theory*. Oxford. Classical vacuum-sphere reference.

#### White / isotropic reflection

- ✓ **LOCAL** Hébert 2009 §3.8.5 — Hébert's `(1−P_ss)^{−1}` factor.
- ✓ **LOCAL** Stammler 1983 Chapters 4-5 — production CP with white BC.
- ◇ **WISH** Sanchez 1977 NSE 64 (canonical interface-current). Already locally on hand as the 1977 paper.

#### Specular (perfect mirror) — Plan 2 case

- ✓ **LOCAL** Sanchez 1986 TTSP 14 — general α/β.
- ◇ **WISH** Davison 1957 Chapter 11.

#### Albedo (partial reflection) — α + β general

- ✓ **LOCAL** Sanchez 1986 — full α + β branch.
- ◇ **WISH** Pomraning, G. C. & Foglio, S. (1965). "The albedo problem for monoenergetic neutrons". *J. Math. Phys.* 6, 1099. Half-space albedo theory.
- — Williams, M. M. R. (1971). *Mathematical Methods in Particle Transport Theory*. Butterworths. Albedo / Wiener-Hopf chapters.

#### Diffuse-reflective / Lambertian

- — Modest, M. F. (2013). *Radiative Heat Transfer* (3rd ed.). Academic Press, Chapter 6-7. Lambertian neutron BC algebraically identical.
- — Howell, J. R., Siegel, R. & Mengüç, M. P. (2011). *Thermal Radiation Heat Transfer* (5th ed.). CRC Press.

#### Periodic

- ✓ **LOCAL** Sanchez 2002 NSE 140 — already on hand.

#### Mixed / spatially-varying BC

- ◇ **WISH** Williams, M. M. R. (1978). "Generalized albedo theory for neutron transport". *Atomkernenergie* 32. Spatially-varying BC.

### 7. Companion math infrastructure

#### Bickley-Naylor & E_n functions

- ◇ **WISH** Abramowitz, M. & Stegun, I. (1964). *Handbook of Mathematical Functions*. NBS, §5 (exponential integrals). Freely available online.
- ◇ **WISH** Knyazev & Selivanov 2014 (already wish-listed under cylinder).

#### Chandrasekhar H functions

- ◇ **WISH** Chandrasekhar, S. (1960). *Radiative Transfer*. Dover. THE H-function reference.
- ◇ **WISH** Mendelson, M. R. (1961). "On the H-function tables". *J. Quant. Spec. Rad. Transfer* 1.

#### Spectral / functional analysis

- ◇ **WISH** Case, K. M. & Zweifel, P. F. (1967). *Linear Transport Theory*. Addison-Wesley. Singular eigenfunction reference.
- ◇ **WISH** Kato, T. (1995). *Perturbation Theory for Linear Operators*. Springer. Spectral compactness for power-iteration convergence proofs.

### 8. Verification benchmarks

- ◇ **WISH** Sood, A., Forster, R. A. & Parsons, D. K. (2003). "Analytical benchmark test set for criticality code verification". *Progress in Nuclear Energy* 42, 55-106. **The standard analytic benchmark suite**.
- ◇ **WISH** Smith, M. A. et al. (2003). "C5G7 MOX benchmark". OECD-NEA. **Standard MoC benchmark — useful for Variant α extension to lattices**.
- ◇ **WISH** Takeda, T. & Ikeda, H. (1991). "3-D neutron transport benchmarks". OECD-NEA. Fast-reactor benchmarks.
- ◇ **WISH** Larsen, E. W. (1982). "F2P benchmark problems". Spatially-non-uniform problems.

## Priority pull list (top-15 by extension impact)

If only a small batch can be pulled at once, recommended order:

| # | Reference                           | Unlocks                             |
| - | ----------------------------------- | ----------------------------------- |
| 1 | Sood-Forster-Parsons 2003           | Standard analytic benchmark suite   |
| 2 | ~~Garcia 2021 *J. Comp. Phys.* 433~~  | LOCAL — multi-region sphere ref     |
| 3 | Davis-Rabinowitz 1984               | Singular-quadrature recipes         |
| 4 | Knyazev-Selivanov 2014              | **Cylinder Variant α 3-D Ki**       |
| 5 | Hébert 2009 full book (Ch. 6-9)     | MG resonance + production MG-CP     |
| 6 | Sanchez-Mao 1993 *AnNE* 20          | 2-D heterogeneous CP foundational   |
| 7 | Lewis & Miller 1984 book            | MG transport reference              |
| 8 | Bell & Glasstone 1970 book          | All-purpose foundational textbook   |
| 9 | Williams 1971 book                  | Fredholm theory + albedo            |
| 10| Adams-Larsen 2002 *PNE* 40          | MG iteration acceleration survey    |
| 11| Hackbusch 1995 book                 | Hypersingular IE theory             |
| 12| Warsa 2004 *NSE* 147                | Krylov for transport                |
| 13| Williams 1983 *AnNE*                | Multi-region sphere                 |
| 14| Wu-Roy 2003 *AnNE*                  | DRAGON cyclic tracking              |
| 15| C5G7 benchmark + Takeda             | Lattice / fast-reactor benchmarks   |

## Recommended pull order by current Plan 2 follow-on phase

**Phase A — V&V gap closure on existing prototype** (vacuum BC + ≥2G):

- All needed references already locally available: Pomraning-Siewert 1982 (vacuum sphere), Sanchez 1986 (full kernel), Hébert 2009 Chapter 3, Stammler 1983 Chapter 6, Atkinson 1997. **No new pulls blocking**.

**Phase B — Multi-region sphere (Issue #132)**:

- **All required references LOCAL**: Brockmann 1981, Mitsis 1963 (Variant γ ground), Sanchez 1986 (single-region theory), **Garcia 2017/2019/2019/2021** (modern P_N from particular-solution to multi-region with reflection — external numerical reference suite for Variant α multi-region cross-check).
- Useful future pull: Williams 1983 *AnNE*.

**Phase C — Cylinder Variant α** (deferred per user; this is plan (c)):

- Local sufficient: Carlvik 1967, Stammler 1983 Chapter 4, Sanchez 1977 (2-D context).
- Required pull: Knyazev-Selivanov 2014 (3-D Ki tables).

**Phase D — Multi-group / production extensions**:

- Required pulls: Hébert 2009 full book (Ch. 6-9), Lewis & Miller 1984.

**Phase E — 2-D Cartesian / lattice topology** (long-term):

- Sanchez 2002 (already on hand) is the load-bearing literature.
- Pull Sanchez-Mao 1993, Postma-Hennart 1989, OpenMOC benchmarks for cross-validation.

## File index

- This memo: `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/peierls_greens_extensions_lit.md`
- Predecessor (Plan 2 B1): `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
- Local PDF directory: `/workspaces/ORPHEUS/scratch/literature/`
- Plan 2: `/workspaces/ORPHEUS/.claude/plans/peierls-greens-function-approach.md`
