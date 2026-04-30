# Cross-Domain Frame Detection: Reference Tables

The tables in this reference are the data that the
cross-domain-attacker agent uses to match problem features to
mathematical frames. The agent's output quality is bounded by
these tables. They are the primary artifact to grow and
maintain.

---

## Part A: Mathematical field → structural trigger table

Every entry is a (frame, trigger, lever, payoff example) tuple.
A frame without a defensible trigger and concrete lever is
rejected. Each subsection covers one branch of mathematics.

### A.1 — Geometry and topology

| Frame                                         | Trigger conditions                                                                                           | Specific levers                                                               | Elegance payoff example                                                                                                                                                                                                              |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Differential geometry / exterior calculus** | Curvilinear coordinates; conservation laws; flux as a quantity crossing surfaces; divergence-based operators | Connection coefficients, Lie derivatives, Stokes' theorem, differential forms | Cylindrical α-redistribution IS a connection coefficient. Parallel-transporting the ordinate direction along curved coordinates produces exactly the M-M term. Signs become consequences of metric choice, not derivation accidents. |
| **Topology (point-set + algebraic)**          | Multiple geometry variants with shared boundary behavior; periodic / reflective BCs; domain continuation     | Manifolds with boundary, covering spaces, homotopy, quotient spaces           | Slab / annulus / hollow sphere as one parameterized manifold with boundary (validated in ORPHEUS). Reflective BC = quotient by reflection group = manifold with corners.                                                             |
| **de Rham cohomology / FEEC**                 | Conservation that must be exact in discretization; mixed methods; compatible finite elements                 | Exact sequences, discrete cohomology, Whitney forms                           | Compatible discretizations where div-curl-grad relations hold exactly at the discrete level. Prevents spurious modes. Relevant for genuine 3D coupled problems.                                                                      |
| **Symplectic geometry**                       | Phase-space flows; characteristic curves; Hamiltonian structure                                              | Canonical forms, symplectic integrators, Poisson brackets                     | Transport characteristics are Hamiltonian. Symplectic integrators preserve phase-space volume exactly — matters for long trajectories in MC with fields or coupled problems.                                                         |

### A.2 — Algebra and representation

| Frame                                            | Trigger conditions                                                                                                   | Specific levers                                                    | Elegance payoff example                                                                                                                                                                                                                                             |
| ------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Group theory / representation theory**         | Any angular discretization; crystal-lattice problems; rotational / reflection symmetry; full-core reflection         | SO(3), Oh, irreducible representations, Peter-Weyl, Clebsch-Gordan | PN is literally the decomposition of angular flux in irreps of SO(3). Peter-Weyl gives correctness for free. Product-Gauss quadrature violates SO(3); Lebedev respects Oh. Explains their different behavior under rotation.                                        |
| **Crystallographic groups / Bloch theory**       | Repeated-lattice problems; full-core with periodic assemblies; Brillouin zone integration                            | Bloch waves, band structure, Brillouin zone quadrature             | Lattice transport problems solved in the unit cell with Bloch BCs, then integrated over the Brillouin zone. Exact dimensionality reduction.                                                                                                                         |
| **Clebsch-Gordan / spherical harmonics algebra** | Scattering kernel decomposition; anisotropic scattering                                                              | Addition theorem, rotation matrices, 3j-symbols                    | Legendre expansion of the scattering kernel is a Clebsch-Gordan decomposition. Truncation error has spectral meaning.                                                                                                                                               |
| **Tensor networks / tensor-train**               | High-dimensional problems (space × angle × energy); compositional boundary structure; cross-section data compression | MPS / TT, Tucker, CP decomposition, hierarchical Tucker            | Boundary conditions as tensor networks (validated in ORPHEUS). Also: TT-format angular flux for very high-N quadrature. Active research area (Peng, Dektor, Einkemmer).                                                                                             |
| **Category theory**                              | Method / functor composition across geometries; lifting constructions                                                | Functors, natural transformations, monads                          | LOW-CONFIDENCE for reactor physics. Abstract-nonsense payoff in PDEs is usually small. Possible exception: unified software architecture where the same solver object works across geometries via functorial lifting. Mark low-signal until a concrete win appears. |

### A.3 — Analysis (functional, spectral, harmonic)

| Frame                                                          | Trigger conditions                                                                        | Specific levers                                                                                      | Elegance payoff example                                                                                                                                                 |
| -------------------------------------------------------------- | ----------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Spectral theory (compact + positive operators)**             | Any eigenvalue problem; iterative methods; operator structure questions                   | Krein-Rutman theorem, Fredholm theorem, perturbation theory                                          | Krein-Rutman is why k-eigenvalue exists and is real-positive with positive eigenfunction. Compactness of CP operator → discrete spectrum → convergence theory for free. |
| **Spectral theory of multiplication operators**                | Operator of form `(M_φ ψ)(µ) = m(µ) · ψ(µ)` with `m: domain → ℝ`; the integration variable µ is the operator's spectral axis | Essential range, ess_sup, multiplication-operator norm theorem (`‖M_φ‖ = ess_sup |m|`; spectrum = ess_range(m)), Trefethen-Embree pseudospectra | Spectrum equals `ess_range(m)`. `(I − M_φ)^{-1}` is bounded iff `1 ∉ closure(ess_range(m))`. Caveat: when `ess_range(m) ∋ 1` (essential range contains 1 over a non-zero-measure set), matrix-Galerkin discretisation does NOT converge in operator norm — adding modes makes `λ_max → 1` and the inverse blows up. ORPHEUS precedent: Phase 5 continuous-µ specular kernel `1/(1 − α·exp(−2aµ))` blows up as µ → 0+ when α → 1; sphere is structurally pathological, slab is structurally immune. The frame names the per-geometry behaviour as one theorem. |
| **Harmonic analysis (Fourier, spherical harmonics, wavelets)** | Periodic BCs; angular discretization; multi-resolution spatial; boundary layer resolution | DFT for circulant diagonalization, SH for PN, wavelets for MR spatial, Chebyshev for boundary layers | Periodic BCs → circulant matrix → diagonalized by DFT → O(N log N) solves. Wavelets give multi-resolution spatial discretization that auto-adapts to gradients.         |
| **Functional analysis / Sobolev theory**                       | Trace / BC analysis; regularity of solutions; well-posedness                              | Trace theorems, embedding theorems, variational formulations                                         | Trace theorems tell you which BCs are well-posed before you discretize. H(div) spaces are the natural setting for mixed methods.                                        |
| **Approximation theory**                                       | Quadrature design; basis choice; convergence rate analysis                                | Best approximation, Chebyshev-vs-Legendre, Padé, Bernstein                                           | Wrong basis can cost orders of magnitude. Rational (Padé) approximation of exp(−Σt·s) in CP is often better than polynomial.                                            |
| **Complex analysis**                                           | Spectral projections; resonance treatment; analytic continuation                          | Contour integrals, residues, Padé                                                                    | Contour-integral spectral projections (FEAST, Sakurai-Sugiura) for parallel eigenvalue solvers. Resonance self-shielding uses contour methods.                          |

### A.4 — Probability, stochastics, information

| Frame                                         | Trigger conditions                                                                         | Specific levers                                                                      | Elegance payoff example                                                                                                                                                                                                                                          |
| --------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Feynman-Kac / stochastic representation**   | Any transport problem. CP, MOC, MC are all here.                                           | Path integrals, Markov kernels, forward / backward Kolmogorov                        | MASTER BRIDGE. CP = first-collision kernel. MOC = deterministic path integration. MC = stochastic path sampling. Delta tracking = rejection sampling on the path measure. Everything reactor physics does with transport is a Feynman-Kac discretization choice. |
| **Markov chain theory / ergodic theory**      | MC convergence; source iteration as MC; power iteration                                    | Mixing times, spectral gap, coupling arguments, Doeblin                              | Source iteration convergence rate = spectral gap of the transport operator. MC equilibration = Markov mixing time. Same math.                                                                                                                                    |
| **Measure theory / Radon-Nikodym**            | Importance sampling; variance reduction; quadrature as measure choice                      | Change of measure, RN derivatives, Girsanov                                          | Importance sampling is literally a Radon-Nikodym change of measure. Reveals what "optimal" means and gives the zero-variance limit formally.                                                                                                                     |
| **Information theory**                        | Angular moment closures; effective sample size; convergence diagnostics; quadrature design | Maximum entropy, KL divergence, mutual information, rate-distortion                  | MN closures (Minerbo, Dubroca-Feugeas) are MaxEnt under angular moment constraints. Discretization error can be framed as rate-distortion.                                                                                                                       |
| **Number theory / low-discrepancy sequences** | QMC; ray generation for MOC; integration error estimates                                   | Van der Corput, Sobol, Halton, Niederreiter, lattice rules, Koksma-Hlawka            | QMC replaces O(N^(−1/2)) MC with O(N^(−1) (log N)^d). For smooth integrands the gain is enormous. MOC ray distributions from QMC improve convergence.                                                                                                            |
| **Cryptography / PRNG theory**                | RNG quality for large-N MC; reproducible parallel sampling; counter-based RNGs             | BigCrush / TestU01, counter-based RNG (Random123 = reduced AES), hash-based sampling | For 10¹²-sample runs, Mersenne Twister is marginal. Counter-based RNGs (Philox, Threefry) are lightweight block ciphers — parallel-safe, reproducible, cryptographically strong. Reframes RNG as a cryptographic-adversary problem.                              |

### A.5 — Optimization, dynamics, control

| Frame                                   | Trigger conditions                                                                    | Specific levers                                                       | Elegance payoff example                                                                                                                                                            |
| --------------------------------------- | ------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Variational calculus / optimization** | Eigenvalue problems; fixed-point iterations; sensitivity                              | Rayleigh quotient, Lagrangian duality, KKT, saddle-point formulations | k-eigenvalue = Rayleigh quotient extremum → variational bounds → adjoint sensitivity for free. Source iteration = coordinate descent on a specific objective.                      |
| **Dynamical systems theory**            | Fixed-point iterations; convergence / divergence diagnosis; bifurcation near critical | Contraction mapping, Lyapunov functions, center manifold theorem      | Source iteration near criticality has a center manifold — that's why convergence degrades. Explains DSA necessity.                                                                 |
| **Control theory / adjoint calculus**   | Sensitivity analysis; perturbation of any parameter; uncertainty quantification       | Pontryagin's principle, costate equations, adjoint sensitivity        | Keff sensitivity coefficients are adjoint flux weighted with perturbations — identical to costate in optimal control. Full sensitivity-adjoint framework for free.                 |
| **Krylov subspace theory**              | Linear system solves; eigenvalue iteration; operator-based methods                    | GMRES, Arnoldi, Lanczos, CG, JFNK                                     | Transport with Krylov acceleration (Warsa, Guthrie, Adams) beats classical source iteration by orders of magnitude. JFNK treats transport + TH coupling without explicit Jacobian. |

### A.6 — Asymptotic and multiscale

| Frame                                         | Trigger conditions                                                            | Specific levers                                               | Elegance payoff example                                                                                                                                                          |
| --------------------------------------------- | ----------------------------------------------------------------------------- | ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Asymptotic analysis / matched asymptotics** | Boundary layers; diffusion limit; thick-thin transitions                      | Inner / outer expansions, matching, WKB                       | Diffusion is the asymptotic limit of transport at thick-system small-absorption. The asymptotic derivation tells you when diffusion fails and what the correction is (SP3, etc). |
| **Homogenization theory**                     | Heterogeneous-lattice → assembly-level; scale separation                      | Periodic homogenization, two-scale convergence, G-convergence | Assembly-level cross sections are a homogenization problem. Two-scale convergence gives the correct homogenized coefficients, not just volume-averaged.                          |
| **Multiscale numerical analysis**             | Lattice problems with disparate scales; local heterogeneity in global reactor | MsFEM, HMM, GFEM, numerical homogenization                    | Connects to AMG / multi-grid design.                                                                                                                                             |

### A.7 — Integral equations and special structure

| Frame                                          | Trigger conditions                                         | Specific levers                                               | Elegance payoff example                                                                                                                                                         |
| ---------------------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Fredholm integral equations**                | CP method; escape probability; any integral-form transport | Neumann series, Fredholm alternative, compact operator theory | CP is literally a Fredholm integral equation with compact kernel. All spectral / existence theorems apply. Born series = source iteration.                                      |
| **Hilbert-Schmidt / separable kernels**        | Kernel of form `K(x, y, µ) = a(x, µ) · b(y, µ) · c(µ)` integrated over an auxiliary variable µ | Degenerate-kernel quadrature, finite-rank-in-µ representation, SVD of kernel slice | The rank-1 axis is the integration variable µ, NOT the spatial indices x, y. Rank in µ is 1 at every µ; total rank = number of µ-quadrature points. Mode-projecting on x, y is projecting along the wrong axis. ORPHEUS precedent: Phase 5 continuous-µ multi-bounce (M1 — #142) — `K_bc^mb[i,j] = Σ_q w_q G(r_i, µ_q) F(r_j, µ_q) f_mb(µ_q)` replaces matrix-Galerkin `(I − T·R)^{-1}` and removes the unbounded-inverse pathology. |
| **Hierarchical matrices / low-rank structure** | CP matrix, integral operators with smooth kernels          | H-matrix, HSS, HODLR, skeleton decomposition                  | CP matrices are often low-rank between distant regions. H-matrix compression gives O(N log N) matrix-vector products. Can make CP competitive with transport sweeps at large N. |
| **Graph theory / algebraic graph theory**      | Cell connectivity; sweep ordering; parallel scheduling     | Topological sort, graph Laplacian, expanders, DAG scheduling  | Parallel transport sweep = DAG scheduling. Laplacian spectral gap bounds diffusion convergence. Expander-like mesh reduces communication.                                       |

---

## Part B: Cross-method pollination map (reactor physics internal)

The unifying structure (from A.4) is Feynman-Kac; below is the
practical pollination between methods in ORPHEUS.

| Target method                    | Borrows from                 | Specific contribution                                                                                  |
| -------------------------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------ |
| **CP method**                    | MOC                          | Ray-tracing formulation equivalent to chord; alternative CP computation                                |
| CP method                        | MC / Woodcock delta tracking | Rejection sampling for CP in highly heterogeneous media; avoids ray-tracing through pseudo-material    |
| CP method                        | Fredholm theory              | Compactness → spectral theorems; Neumann series → source iteration; low-rank structure for H-matrix CP |
| CP method                        | MaxEnt                       | Angular closure for CP-like methods with partial angular info                                          |
| **MOC**                          | MC                           | QMC ray generation (Sobol, lattice) for better integration                                             |
| MOC                              | Symplectic                   | Symplectic ray integrators for coupled problems with drift                                             |
| MOC                              | Graph theory                 | Ray-sweep scheduling; cyclic ray generation as Euler path                                              |
| **MC**                           | CP                           | Deterministic variance reduction via CP-estimated importance                                           |
| MC                               | QMC / number theory          | Low-discrepancy sampling; randomized QMC for unbiased estimators                                       |
| MC                               | Cryptography                 | Counter-based RNG for parallel reproducibility                                                         |
| MC                               | Adjoint sensitivity          | Adjoint-weighted tallies, CADIS / FW-CADIS, variance reduction by adjoint importance                   |
| MC                               | Markov chain theory          | Autocorrelation, effective sample size, Gelman-Rubin diagnostics                                       |
| **SN**                           | PN / spherical harmonics     | Angular closure; SP_N methods; filtered SN                                                             |
| SN                               | Differential geometry        | Connection coefficients for curvilinear; coordinate-independent formulation                            |
| SN                               | Group theory                 | Lebedev quadrature (Oh-symmetric); product-Gauss violates SO(3)                                        |
| SN                               | Krylov                       | GMRES-accelerated transport; preconditioned Richardson                                                 |
| SN                               | Homogenization               | Subgroup methods, response-matrix coupling                                                             |
| **Diffusion**                    | Asymptotic analysis          | Correct derivation from transport; when it fails (SP3 correction)                                      |
| Diffusion                        | Graph Laplacian              | Spectral properties, multigrid design                                                                  |
| **Eigenvalue / power iteration** | Optimization                 | Rayleigh quotient → variational acceleration; Chebyshev iteration                                      |
| Eigenvalue                       | Dynamical systems            | Center manifold near criticality explains DSA necessity                                                |
| Eigenvalue                       | Krylov                       | Arnoldi / Lanczos acceleration; IRAM for multiple eigenvalues                                          |
| **Sensitivity / adjoint**        | Control theory               | Full Pontryagin framework; second-order sensitivities                                                  |
| Sensitivity                      | Automatic differentiation    | Operator-level AD avoids hand-derivation                                                               |
| **Coupled neutronics + TH**      | JFNK / nonlinear solvers     | Newton on the coupled residual; operator splitting error analysis                                      |
| Coupled                          | Schur complement             | Reduce coupled system to one field exactly                                                             |
| **Resonance self-shielding**     | Complex analysis             | Contour-integral treatment; R-matrix theory                                                            |
| Resonance                        | Homogenization               | Subgroup method as numerical homogenization                                                            |

---

## Part C: Elegance detector

Heuristics for detecting that the current formulation is NOT
yet elegant and a frame-retrieval pass is warranted. If two or
more smells fire, the cross-domain-attacker should probe.

### Structural smells

1. **Repeated near-identical code across geometry variants.**
   Signals a unified geometric frame exists but hasn't been
   found. (Slab / cylinder / sphere with three duplicated sweep
   kernels → differential geometry or topology probe.)
2. **A dense matrix that's being explicitly formed.** Often a
   tensor-network or low-rank structure is hiding. (BC
   decomposition, CP matrix.)
3. **A "magic" correction term with hand-waved derivation.**
   Often a connection coefficient, Jacobian, or measure-change
   in disguise. (Cylindrical redistribution.)
4. **Symmetry present in the problem but absent in the method.**
   Group-theoretic discretization being missed. (Product-Gauss
   on a sphere.)
5. **Convergence rate suboptimal by a known margin.** Often an
   approximation-theoretic or basis choice issue. (Legendre
   when Chebyshev was right.)
6. **Something described as "iterative" without its fixed-point
   structure analyzed.** Dynamical systems / spectral theory
   will tell you the convergence rate in closed form.
7. **"We picked N because it works."** Often a
   representation-theoretic or information-theoretic optimum
   exists. (Quadrature order, moment closure.)
8. **Nested loops with algebraic structure.** Often a tensor
   contraction that could be reformulated.
9. **Boundary handling as a special case added to interior
   logic.** Trace theorem / cohomological structure being
   missed. (BCs as separate code paths instead of part of the
   operator.)
10. **Explicit bookkeeping of which-path / which-history** in a
    Monte Carlo or iterative method. Often has a
    measure-theoretic / Girsanov reformulation.

### Semantic smells

11. **"We had to add a stabilization term."** Usually indicates
    the natural variational setting is wrong. Often needs
    H(div), H(curl), or a mixed formulation.
12. **"We had to add a small ε to avoid division by zero."**
    Often signals the native formulation doesn't have the
    singularity (projective coordinates, exterior calculus).
13. **"This converges because we checked."** No theorem —
    means the contraction / spectral-gap argument hasn't been
    found.
14. **"We normalize at each step to keep it stable."** Often
    Rayleigh quotient structure is being missed.
15. **Rank-N non-monotone convergence.** When increasing mode
    count produces non-monotone error (rank-(N+1) error >
    rank-N error, or convergence flatlines and then worsens),
    Galerkin projection is being applied without a
    Rayleigh-Ritz / variational principle. Rayleigh-Ritz with
    nested subspaces guarantees monotone convergence from
    above (Courant-Fischer min-max); non-monotone behaviour
    means the Ritz frame was abandoned. Probe A.5 variational
    framing before adding more modes — the bug is not in the
    truncation order, it is in the projection. Trigger: any
    rank-truncation method (rank-N closure, finite-mode
    expansion, truncated Karhunen-Loève) where doubling the
    modes worsens the error. ORPHEUS precedent: rank-N closure
    investigation (Direction-C/Q failures, #121, #122 closed)
    where rank-1 F.4 gave 0.003% but rank-2 Marshak gave 1.36%.

---

## Growth Protocol

Tables grow only from evidence. Speculative entries degrade
the signal of the table and lead to noise-output from the
cross-domain-attacker agent.

**Adding a trigger-table entry (Part A):**

1. The entry must have a confirmed or high-prior
   frame-problem match in the ORPHEUS work (or in well-cited
   reactor physics literature)
2. The trigger must be nameable in one sentence
3. The lever must list at least three specific mathematical
   tools within the frame
4. The payoff example must be concrete, not generic

**Adding a pollination entry (Part B):**

1. The borrowing must be concrete — a specific algorithmic
   technique, not a vague analogy
2. The direction must be named (method X borrows from method Y,
   not "X and Y are related")

**Adding an elegance smell (Part C):**

1. The smell must have a paired example from ORPHEUS or
   similar codebases
2. The example must point to a real reformulation (preferably
   in scripts/)

## Version History

- v0 (initial): tables authored by cross-domain-attacker agent
  design session. All entries either validated in ORPHEUS or
  high-prior based on standard mathematical applicability to
  reactor physics.

Record each subsequent edit below with date, entry added /
modified, and justification (validated payoff or high-prior
rationale).

- 2026-04-30:
  - Added Smell #15 (Rank-N non-monotone convergence) to
    Part C. Justification: validated payoff on rank-N closure
    investigation (#121, #122 closed) where rank-2 Marshak
    was worse than rank-1 F.4.
  - Added A.7 row "Hilbert-Schmidt / separable kernels" with
    explicit rank-axis = integration variable. Justification:
    validated payoff on Phase 5 continuous-µ multi-bounce
    (#142) — production M1 reformulation removes
    matrix-Galerkin unbounded-inverse pathology.
  - Added A.3 row "Spectral theory of multiplication operators"
    with `ess_range(m) ∋ 1` non-convergence caveat.
    Justification: validated diagnostic on Phase 5 continuous-µ
    sphere kernel — names the per-geometry pathology in one
    theorem.
  - New precedent file:
    `scripts/validated_hilbert_schmidt_separable.md` (M1/M2/M4
    verification triad).
