# ORPHEUS Operator Machinery: Theory, Architecture, and Action Plan — v2

**Consolidated deep report.** Supersedes v1 (2026-08-02), `orpheus-canonical-naming.md`, and the actionable items of `transport-graph-information-report.md`. This revision absorbs the space-and-cone conversation of 2026-08-06/08: flux ontology (cone, not torsor), the space taxonomy and its audit against the live `numerics/space.py` / `numerics/spaces/` tree, the retract lattice, the cone-theoretic V&V battery, partition ownership, complexification, and the resulting re-sequencing of the action plan. v1 content is retained and silently amended where sharpened; every overturned claim is recorded in Part II.

Reviewed against the live tree (`orpheus/numerics`, `orpheus/sn`, `orpheus/transport`), `docs/theory/foundations/{operator_algebra,operator_inverse_family,path_integral,frame}.rst`, the `#226` rulings in `.claude/agent-memory/`, and `neutron_transport_grand_report_v3.md`. Space-layer audit read from `numerics/space.py`, `numerics/spaces/{angular_trace_space,spherical_harmonic_space,full_field_space}.py`, `numerics/coupled_system.py`, `numerics/operator.py`, `sn/mesh/method_space.py`.

Date: 2026-08-08

> ⚠⚠ **GROUNDED AGAINST THE TREE 2026-08-19 — read Part VI before designing any step.** Four explorer reconciliations measured every tree-facing claim at HEAD `7aae9bf1` (memos: `scratch/omr_v2_grounding/{A,B,C,D}_*.md` — the evidence base under Part VI's summary). Headlines: Part I's theory stands, but §I.3's tree diagnosis was stale at its own write time (it retained v1's operator audit while #261 had landed the binding machinery six weeks earlier), and §I.7's overturn is a pending RULING against an *enforced runtime architecture*, not prose. Part III-S is current but incomplete (three space files never read). Part V's Phases 1–3 overlap the live, tracked, user-ruled `operator_strategy_realization_campaign.md` (P0 merged, P1 live). ⛔/✅/⚠ banners at the affected claims point into Part VI; per plan-authoring §3 the original text stays.

---

# Part I — The theoretical spine

Eleven results carry everything else. Each is stated with its consequence for ORPHEUS. I.1–I.6 are v1's spine, amended in place; I.7–I.11 are new and constitute the **space substrate** — the layer beneath the pipeline that the kernel-binding refactor (`.on(V)`) silently assumes and that the live tree currently gets wrong (Part III-S).

## I.1 One recursion: it is splittings all the way down

Every deterministic transport solve is the same move applied recursively: **partition the operator, invert the part you can, iterate on the rest.**

$$A = M - N, \qquad x_{k+1} = M^{-1}(Nx_k + q), \qquad \text{cost} = \rho(M^{-1}N)$$

The levels differ only in *which cycles* the split resolves and *what you want from the iteration map* $x \mapsto M^{-1}Nx$:

| Level | $M$ | $N$ | Cycle resolved | Want |
|---|---|---|---|---|
| **Pencil** (eigen) | $A(\sigma) = A - \sigma\Pi$ | $\Pi$ | the fission/production loop | the **number** $\rho = k$ (normalize each step) |
| energy | downscatter block-triangle | upscatter | thermal upscatter loop | $\rho < 1$ (fixed point) |
| boundary trace | $(L{+}C{-}B_{\text{lower}})$ | $B_{\text{upper}}$ | reflection loop | $\rho < 1$ |
| angle | $(L{+}C)$ | $S$ | scattering loop | $\rho < 1$ |
| **causal core** | the whole thing | — | none left | terminal: $\rho = 0$ |

Two precise points:

1. **The pencil and the splitting share the iteration map, not the semantics.** A splitting iterates to a fixed point (Neumann series, needs $\rho < 1$); power iteration normalizes each step and converges to the dominant *direction* (works for any $\rho$, and $\rho$ is the answer). Same map, different use. This is exactly why one engine serves both — and why `Splitting.spectral_radius()` is free (§I.5).

2. **Termination invariant.** The recursion bottoms out at the causal core *because quasinilpotency is the only $\rho = 0$*. This gives a checkable well-formedness condition: **a splitting stack is well-formed iff its innermost $M$ is causal (sweepable) or small enough to invert densely.** If you cannot reach a causal or materializable core, you are not done splitting. This should be an assertion, not a convention.

The variable-agnostic construction underneath all rows: compute the SCC condensation of the block digraph, topologically order it, put the block-lower-triangular part in $M$, the rest in $N$. Jacobi and Gauss–Seidel are the two extremes of one construction (trivial order vs. full topological order). Downscatter is not "lagged" — it is *resolved by ordering*, the same operation as the spatial sweep on a different index set. **Four apparently different techniques (spatial sweep, group sweep, boundary schedule, source iteration) are one algorithm parameterized by index set.** The splitting machinery should be written once, variable-agnostically.

**v2 amendment — the three-layer refinement of "block partition."** v1 said the partition belongs to the space, the splitting to the strategy. That survives, sharpened into three data with two owners (full statement §I.10): the *partition* (unordered decomposition, per-axis, space-owned, stage 2), the *block digraph and its topological order* (computed from operator + partition — not space data, because fast→thermal exists only by the scattering kernel's triangularity), and the *schedule* (which blocks implicit/lagged — splitting-owned, stage 5). Conflating the middle layer into either neighbor is where implementations rot. KBA/wavefront parallelism falls out with no new machinery: it is the topological *level sets* of the spatial-partition block digraph per ordinate — the same construction that gives energy-GS its group order. Domain decomposition is a partition of the spatial axis; the machinery is variable-agnostic at the space level too.

## I.2 The Volterra/Fredholm dichotomy — the classification behind the inverse family

The four `#226` inverse members are not four algorithms; they are the classical integral-equation types, keyed by the structure of the forward:

| Forward structure | Analytic class | Inverse realization | ORPHEUS class |
|---|---|---|---|
| diagonal / multiplier | — | pointwise reciprocal | `InverseOperator` |
| **causal / triangular** | **Volterra** | substitution (the sweep) | `SweepOperator` |
| identity − compact | **Fredholm 2nd kind** | Neumann series, needs $\rho < 1$ | `GreenOperator` |
| finite, materializable | dense | LU | `MatrixInverseOperator` |

**The Fredholm row is exact, not an analogy.** Preconditioning the within-group loss:

$$(I - (L{+}C)^{-1}S)\,\psi = (L{+}C)^{-1}q$$

is identity minus compact; `GreenOperator`'s convergence condition $\rho((L{+}C)^{-1}S) < 1$ *is* the Fredholm alternative. CP's $\phi - cP\phi = PS$ is the same equation type — which is why CP's *solve* is well-posed while its *deconvolution* (source recovery) is ill-posed. Both live on the same row.

**The Volterra row carries theorems, not just a label:**

- $(L{+}C)$ is **not** a Volterra operator — it is a first-order differential operator, the *generator*. Its **inverse** is the Volterra one: along a characteristic, $[(L{+}C)^{-1}q](s) = \int_0^s e^{-\tau(s',s)} q(s')\,ds'$, causal kernel. Forward = generator; inverse = propagator.
- **Quasinilpotent on a bounded domain**: $\rho = 0$ exactly; discretely $\rho = O(\Delta x)$ (triangular ⟹ eigenvalues = diagonal = $1/(|\mu|/\Delta x + \Sigma_t)$). Meanwhile $\|K\| = 1/\Sigma_t$, an $O(1)$ mesh-independent quantity. **Total non-normality: the eigenvalues of the sweep tell you nothing about what the sweep does**, and the gap grows under refinement.
- Consequence: **all of source iteration's spectral radius originates in $S$ and the boundary closure — never in streaming.** Multiplying a quasinilpotent by a bounded operator is what lifts the spectrum off zero. The classical bound $\rho_{\rm SI} \le c$ is exactly the submultiplicative $\|K\|\cdot\|S\| = (1/\Sigma_t)(\Sigma_s)$.
- The quasinilpotency mechanism is *causality kills combinatorics*: $\|K^n\| \le (CL)^n/n!$ because the ordered simplex has volume $s^n/n!$. Cyclic graphs grow closed walks like $\rho^n$; DAGs like $1/n!$.

**v2 sharpening — the sweep is literally an inverse, and its type is a theorem.** The pair $(T, \gamma^-) : W \to L^2 \times \Gamma^-$ — "apply collision-streaming, take incoming trace" — is a **bijection**; that is the well-posedness theorem for transport with incoming data. The sweep $\Lambda : L^2 \times \Gamma^- \to W$ is its inverse, characterized by $T\Lambda(q,g) = q$ and $\gamma^-\Lambda(q,g) = g$. The sweep is not "an operator that happens to solve something"; it is *the* inverse arrow of a specific bijection, and its type declares what it consumes (emission density + incoming trace) and what it produces (a $W$-element — something traceable, §I.8). This is also the deepest form of the no-`TrackSweepOperator` ruling (correction #10): MoC inverts the *same bijection* along a different realization of the spatial axis.

**The transition between the rows is one operation with four consequences.** Angular integration ($\int d\Omega$, or equivalently including $S$) simultaneously creates: the cycles ($G(L_{-m}) = G(L_m)^T$; union of a DAG with its transpose is complete — *isotropy creates the cycles*), the compactness (transverse smoothing), the symmetrizability (CP reciprocity ⟹ $DP$ symmetric ⟹ free adjoint), and the ill-posedness of source recovery. **One move separates $S_N$ from CP.**

The exclusivity underneath: **triangular ↔ symmetric** — an operator that is both is diagonal. $S_N$ streaming is maximally non-normal (DAG, $O(N)$, adjoint = second sweep); CP is $D$-self-adjoint (dense, cyclic, $O(N^2)$, adjoint free). You cannot have both, and the choice of which to sacrifice *is* the method taxonomy.

## I.3 Kernel ≠ operator — the Schwartz binding and its correctness condition

**The diagnosis.** `ScatteringOperator`, `FissionOperator`, `CollisionOperator` have no unique $(domain, codomain)$ pair and dispatch on the argument at apply time. They are **kernels**, not operators. The Schwartz kernel theorem: linear maps correspond bijectively to distributions on the product space — *once the spaces are fixed*. A kernel $K(x,y)$ is data; it becomes an operator only via $(Kf)(x) = \int K(x,y) f(y)\,d\mu(y)$, **and $\mu$ is not part of the kernel**. The current apply-time polymorphism is the missing $\mu$ leaking out at call time.

> A kernel is representation-free. An operator is representation-bound. The binding is a choice, and it is not unique. Polymorphism belongs at construction (`.on(V)`), not at apply.

**Two distinct causes — do not unify the fix:**

| Cause | Leaves | Fix |
|---|---|---|
| measure/representation ambiguity (genuinely mathematical) | $S$, $F$ | bind to the frame |
| carrier-shape ambiguity (purely representational) | $C$, ndarray paths | one carrier type at the boundary |

Three kernel kinds, three binding mechanisms: $C$ is a multiplier (binding is a cast), $S/F$ are integral kernels (binding = frame + quadrature measure), $L$ is differential (**binding = the closure** — the Padé entry).

**v2 amendment — who owns the frame and the measure.** v1 said "`.on(V)` supplies $M, R$ from $V$'s frame" without asking whether $V$ can actually answer. The frame, the quadrature measure, the harmonic mass matrix — all of it lives on the **space's axes** (§I.8), and the live space layer currently cannot carry it correctly: identity excludes the metric, metrics bind to axes by broadcast convention, the bulk volume weight is "absorbed into operators" (audit F1–F3, Part III-S). **The space refactor is therefore a prerequisite of the binding refactor, not a parallel track.** `.on(V)` is only as well-defined as $V$'s axes are. This re-sequences Part V.

> ⛔ **PARTIALLY REFUTED 2026-08-19 (memo C).** For the ANGULAR binding the tree ratified the OPPOSITE ownership before this plan was written: the `frame-eigenbasis-ownership` ruling (`d30d4a68`, 2026-06-25, `docs/theory/foundations/frame.rst`) — the scattering operator OWNS its Funk–Hecke eigenbasis frame, and `kernel = frame.conjugate(Λ)` has been the production R∘Λ∘M binding since `cc6d022e` (2026-06-24), with no Phase S. The prerequisite survives only where binding reads genuine SPACE data (bulk metrics, identity, partitions). Adjudicate ownership (Part VI.5 D2) before acting on this paragraph.

**$R\Lambda M$ is the definition of the binding, not an optimization.**

$$S = \underbrace{R}_{\text{synthesis, from } V}\ \underbrace{\Lambda}_{\text{kernel — no measure}}\ \underbrace{M}_{\text{analysis, carries } \mu}$$

$\Lambda$ is what `ScatteringKernel` holds; $M, R$ are what `.on(V)` supplies from $V$'s axes. $F$ is the rank-1 case $|\chi/4\pi\rangle\langle\nu\Sigma_f|$; anisotropic $S$ is rank $(L{+}1)^2$; the analysis functional is precisely the object that needs the measure — and precisely the one currently chosen at apply time. (Escape hatch: a tabulated non-expandable kernel takes $\Lambda$ = full matrix, $M = R = I$; the design exploits finite rank, never assumes it.)

**The binding is Petrov–Galerkin realization, and the spine theorem is its correctness condition.**

$$\text{bind}(K)^\dagger = \text{bind}(K^\dagger) \iff R = M^* \iff \text{the frame is tight}$$

The $B/A$ tightness ratio is exactly the size of the discrepancy when it fails. **This converts ERR-039 from a bug into an instance of a theorem.** Gate every kernel binding on the tightness check and this bug *class* becomes unrepresentable rather than tested-for.

**v2 addendum — the tightness family.** The same $B/A$-shaped invariant ("a computable number certifying when adjoint = dual") now appears at four independent joints, and should be documented as one family: (i) the kernel-binding gate above; (ii) the axis-level rule **`.T` is about labels and layout, `.H` is about metrics** — they coincide iff the axis metric is identity, i.e. the orthonormal/tight case (§I.8); (iii) spectral-mode extraction — biorthogonality $\langle\varphi^\dagger_i, \Pi\varphi_j\rangle = \delta_{ij}$ realizes the dual basis via adjoint modes, and a wrong weighting is the adjoint-vs-dual failure in spectral clothing (§I.11); (iv) the energy axis — multigroup forward fluxes are covariant (group integrals), adjoint fluxes contravariant (spectrum-weighted averages), so the adjoint of the multigroup problem equals the multigroup discretization of the adjoint iff the collapse weightings are bilinearly consistent: **the classic group-collapse adjoint inconsistency is a tightness failure on the energy axis**, sibling to ERR-039 on the angular axis.

**The `LossRepresentation` split rule** stands as in v1: anything that changes the **answer** belongs to the operator (realization); anything that only changes the **cost** belongs to the strategy (solution). Closure → binding of $L$; traversal → strategy; ERR-026 closes structurally. v2 adds the connection to §I.8: **the closure is also what gives the discrete solution space its trace content** — diamond difference's auxiliary edge unknowns exist precisely so the discrete space can be `WithTrace`. Closure choice and trace representation are one decision made at Realization.

## I.4 The pencil and the resolvent

$$\text{Pencil}(A, \Pi): \quad A\psi = \lambda\,\Pi\,\psi, \qquad R(\sigma) = (A - \sigma\Pi)^{-1} = \texttt{pencil.at(σ).inverse()}$$

**The resolvent is not a new type — and not a new space.** At fixed shift, $R(\sigma)$ is an arrow **codomain → domain of the posed operator** — for the fully posed problem, $(L^2 \oplus \Gamma^-) \to W$, the same type $\Lambda$ generalizes. As a family, the resolvent is an analytic operator-valued function on $\mathbb{C} \setminus \text{spec}(A,\Pi)$; the Pencil reifies the family (recipe = structure, `.at(σ)` = evaluation). Four consumers of one object:

| Consumer | Spelling | Status |
|---|---|---|
| $k$ power iteration | `Pencil(A_loss, F).at(0).inverse() @ F` | shipped, as `A_loss.inverse() @ F` |
| $\alpha$ shift-invert | `Pencil(A_prompt, T).at(σ).inverse()` | grand report §36.2, as separate machinery |
| higher harmonics | Riesz projection $\frac{1}{2\pi i}\oint_\Gamma R(\lambda)\,d\lambda$ | not designed |
| **reactor noise** | $R(i\omega)$ — the core's transfer function | not designed |

**v2 amendment — the noise row and the contour row demand complexification.** $\sigma = i\omega$ (and complex contour points, and complex conjugate-pair harmonics of the nonsymmetric pencil) require the space's **ground field** $\mathbb{F} \in \{\mathbb{R}, \mathbb{C}\}$ to be part of the taxonomy. Complexification $V \mapsto V_\mathbb{C}$ is a functor: same axes, same metric data, pairing extended *sesquilinearly* — which means the complexified `.H` performs a conjugation the current all-real metric path never does. It comes with canonical structure worth typing: the conjugation operator $\mathcal{C}$ (antilinear involution) whose fixed points recover $V$; the physical reality constraint on noise, $\delta\varphi(-\omega) = \mathcal{C}\,\delta\varphi(\omega)$, is only spellable if $\mathcal{C}$ exists as an object, and a noise solve violating it has a bug. What does **not** survive complexification is informative: the cone (no order on $\mathbb{C}$ — so no positivity assertions on noise fields, correctly), and hence all of §I.7's Tier-1/2 battery, which must be gated on $\mathbb{F} = \mathbb{R}$. Items 3.2/4.1/4.2 gain this dependency (Part V).

**The shift-absorption rule** stands: multiplier $\Pi$ → folds into the collision diagonal, $M$ changes values, structure intact, **shift is free** ($\alpha$: $\Pi = M[1/v]$). Integral $\Pi$ → rides as another **gain** in $N$; the sweep is untouched but the inner $\rho$ grows with $\sigma$ ($k$/Wielandt — the classical inner/outer trade, which the rule predicts).

**Existing machinery is already the admissibility test.** `StreamingCollisionOperator.__init__` validates $\sigma > 0$ strictly; under an $\alpha$ shift the guard becomes $\Sigma_t - \sigma/v > 0$, failing at the streaming decay rate — the edge of the $\alpha$ essential spectrum in the infinite/periodic case. The error message needs an $\alpha$-aware branch; shift-posability gets a cheap pre-check. *(Bounded-domain caveat held as in v1.)*

**The $k/\alpha$ asymmetry is structural.** $k$ always exists (Krein–Rutman on a positive quasi-compact operator — the compactness hypothesis is delivered by velocity averaging, D4). The $\alpha$-spectrum has a continuous part; discrete $\alpha$'s may be finitely many **or none**; a naive power iteration can converge to nothing or to the continuum edge. Posing-table precondition. **v2 adds the sibling precondition: IRREDUCIBLE** (§I.7 Tier 0) — primitivity of the discrete fission map is a property of the *posed problem*, checkable as connectivity of the fission-coupling block digraph, and it gates which cone assertions are active. Both preconditions live at Posing.

## I.5 Conditioning: three regimes, and the foresight/hindsight connection

| Regime | A priori fraction | Predictor | Where |
|---|---|---|---|
| sweep numerics | **100%** | Padé table: $a(\tau)$, $\tau = \Sigma_t\Delta s/|\mu|$; negative flux iff $\tau > 2$ (diamond) | mesh-generation time, per cell per ordinate |
| sweep parallelism | **100%** | critical-path depth of the cell DAG: $O(3N)$ for $N^3$ structured | mesh topology only |
| outer iteration | **~80%** | $\rho_{\rm SI} = c\,\arctan(\lambda)/\lambda$ at the fundamental buckling; the missing 20% is heterogeneity | one coarse diffusion eigensolve |

**Source iteration is explicit pseudo-time-stepping in collision number**; DSA is the implicit scheme (one global elliptic solve = infinite information speed = iteration count decoupled from domain size). **Perron–Frobenius proves the strictness of $\rho < c$**: the row sum of $(L{+}C)^{-1}S$ is the expected retained secondaries — $c$ in the interior, strictly less within a mean free path of vacuum. **Leakage is row-sum deficiency where the digraph meets the vacuum sink.**

**Foresight and hindsight are the same computation on different operators**: loss of translation invariance is exactly loss of foresight, so the predicted-vs-measured gap is a heterogeneity metric; hindsight is information-bounded ($k$ matvecs ⟹ $2k$ moments; measuring $\rho$ converges no faster than the iteration itself). **Foresight predicts the rate, never the failure mode** (inconsistent DSA is the canonical counterexample).

**v2 amendment — the cone battery completes this section's diagnostics** (full statement §I.7). Three additions upgrade "estimates" to "theorems with runtime witnesses": (i) **Collatz–Wielandt brackets** — every power iterate yields *rigorous two-sided bounds* on $k$, monotonically shrinking under primitivity; strictly stronger than the $k$-history everyone plots, and a bracket violation is an unambiguous bug, no tolerance judgment. (ii) **Birkhoff contraction** — kernel bounds give a certified dominance-ratio ceiling $d \le \tanh(\Delta/4)$ *before running anything*; complementary to Fourier (Fourier gives the expected rate on model problems, Birkhoff a theorem-backed bound on the actual problem; for loosely coupled cores the bound is valid but slack — its value is as a certified ceiling). (iii) **SI cone-monotonicity** — from a zero initial guess, unaccelerated SI converges *pointwise monotonically from below* (Neumann series of positive terms); DSA legitimately breaks the cone-monotone path, so the test doubles as a certificate that the unaccelerated path is truly unaccelerated.

## I.6 Cycles: two canonical moves, rank decides

Stage 2(c) of the splitting has exactly two moves: **split** (lag the closing edges, iterate at cost $\rho$) and **close** (form the once-around monodromy $\Phi$, invert $(I-\Phi)$ on the trace; Woodbury closes an SCC in $r{+}1$ sweeps plus an $r{\times}r$ solve). **The criterion is $\text{rank}(\Phi)$.**

| Configuration | Acyclic? | $\text{rank}(\Phi)$ | Move |
|---|---|---|---|
| vacuum \| vacuum | yes | — | sweep |
| reflective \| vacuum (either side) | yes | — | sweep, correct order first |
| 2-D quarter symmetry (W+S refl, E+N vac) | yes — 4-orthant DAG | — | **one pass, zero boundary iterations** |
| white / albedo per face | no | **1** | close: 2 sweeps + scalar division |
| **periodic, 1-D** | no | **1 per ordinate** | **close: 2 sweeps + one division per ordinate** |
| periodic in x, 2-D | no | 1 per (ordinate, row) | close row-wise in upwind-y order |
| reflective \| reflective | no | $N_{\rm ord}/2$ | split and iterate |

The unifying stage-2(c) statement: **SCC-decompose; per component compute the closing-edge rank; close the low-rank exactly, split only the rest.** `sweep_acyclicity` already computes the components; the rank computation is one addition, generalizing Issue #300.

**v2 note.** The white/albedo rank-1 row now has its structural explanation from the retract lattice (§I.9): the white BC *is* a rank-one-in-angle composite $\iota \circ \alpha \circ \pi$ through the partial-current retract — the monodromy rank is the retract's codimension-collapse made visible at the boundary.

## I.7 Flux ontology: the positive cone, and the V&V battery it generates *(new)*

> ⚠⚠ **RULING PENDING — SCOPE MEASURED 2026-08-19 (memo D).** This overturn's object is not prose: the torsor is a dunder-ENFORCED runtime algebra — [M] 16 production files / 16 test files / 5 doc pages (`FluxRole` mixin; 7 flux ↔ 7 displacement leaf pairs; 2 production consumers of the Displacement type; `field_algebra.rst`, 602 lines, 4 `:eq:`-cited labels). Crux measurements for the ruling: NO site in production, tests, or examples adds two flux-typed fields (superposition lives entirely at the source level — the incumbent blocks no live use case today; the live demand is operator-domain, #331), while the cone side already ships in fragments (`is_positivity_preserving` with DD honestly `False` + a numerical witness, a ψ≥0 realizer refusal, the coefficient cone-as-predicate battery, ray normalization in `power_iteration`). See Part VI.5 D1 + memo D §4/§9 before executing 0.7 or D7.

**The torsor doctrine is overturned.** The codebase's member-family docs (visible verbatim in `coupled_system.py`: "the affine flux torsor ... appl[ies] unchanged") assert flux is a torsor in affine space, motivated by "differences make sense, sums don't." Two facts kill it: flux has a **canonical zero** (vacuum — the absence of neutrons is not a choice of origin), and **superposition is real** (for fixed cross sections the transport operator is linear: fluxes of independent sources in the same medium add; what doesn't add is fluxes of *different problems*, which is fiber discipline, not affine structure). The correct ontology:

> **Flux lives in the positive cone $K$ of an ordered vector space $V$.** Differences $\varphi_1 - \varphi_2$ leave $K$ but stay in $V$ — they are perturbation/error fields, the same tensor type without the positivity predicate. Positive scaling is meaningful; the $k$-eigenfunction is a **ray** in the projective cone $P(K)$, its magnitude fixed only by power normalization. The cone is load-bearing: Krein–Rutman (existence/uniqueness/simplicity/strict positivity of the fundamental mode) is cone theory, and power-iteration convergence is Birkhoff–Hopf contraction in the Hilbert projective metric *on the cone*.

**Cone membership is a predicate on elements, not a subspace, and cone *representability* is an axis property** (§I.8): nodal axes have coordinate cones (positivity coefficient-wise); modal axes do not (harmonic-moment coefficient positivity is neither necessary nor sufficient for function positivity — and this dichotomy *is* the ray-effect/negative-source dichotomy: positivity native to the quadrature axis, rotational equivariance native to the harmonic axis, no angular basis has both). Contraction against a positive weight never degrades cone representability — the ℓ=0 retract of the moment space regains the coordinate cone.

**Positivity-preservation is a property of arrows worth tracking as a flag** — and it is a property of the *realization*, not just the continuous operator: **diamond difference does not preserve the cone**; step and step-characteristics do. A DD sweep is a map $V \to V$ approximating a map $K \to K$, and the Realization stage should record exactly that.

**The battery.** Organized by what each theorem demands, which determines where the test lives. The epistemic gradient — Tier 1 exact inequalities always valid; Tier 2 exact inequalities under a checkable gate; Tier 3 rate bounds under kernel bounds — should be visible in the test suite's organization.

*Tier 0 — preconditions as tests (gates, live at Posing):*

- **M-MATRIX** — inverse-positivity of the discretized $(\Omega\cdot\nabla + \Sigma_t)$: sign-pattern check on small problems, MONOTONE-SOLVE probes on large. Licenses Tier 1's order statements.
- **IRREDUCIBLE** — primitivity of the discrete fission map = connectivity of the fission-coupling block digraph. Void channels, pure absorbers, decoupled regions break it, and when it breaks, CW-MONOTONE and SIGN-STRUCTURE *legitimately* fail. The gate decides which Tier-2 assertions are active per problem.

*Tier 1 — per-apply invariants (positivity only):*

- **CONE-PRESERVE** — every operator flagged positive maps $K \to K$ on random nonnegative inputs. The most violated invariant in deterministic transport (DD); the flag must be allowed to be honestly `False` per realization.
- **HILBERT-NONEXPANSIVE** — any positive $A$: $d_H(A\varphi, A\psi) \le d_H(\varphi,\psi)$ in the Hilbert projective metric $d_H(\varphi,\psi) = \log[(\sup\varphi/\psi)(\sup\psi/\varphi)]$. Assumption-free Birkhoff; two pointwise ratios; scaling-invariant, hence insensitive to normalization convention.
- **MONOTONE-SOLVE** — inverse-positive discretization: $q_1 \le q_2 \Rightarrow \varphi_1 \le \varphi_2$ pointwise. Exact inequality, no tolerance.

*Tier 2 — per-iteration invariants of power iteration (require IRREDUCIBLE):*

- **COLLATZ-WIELANDT-BRACKET** — $\inf_x (F\varphi)(x)/\varphi(x) \le k \le \sup_x (F\varphi)(x)/\varphi(x)$ at *every* iterate, including iterate zero. Rigorous two-sided bounds; the converged $k$ must lie inside every bracket ever computed.
- **CW-MONOTONE** — lower bound nondecreasing, upper nonincreasing; a single inversion is a hard failure.
- **SIGN-STRUCTURE-OF-ERROR** — converged $\varphi$ strictly positive on the irreducible component; the error field $\delta\varphi^n = \varphi^n - \varphi^\infty$ must change sign (only the fundamental is in $K$; the spectral projector is not a positive operator). A solver whose error field is one-signed is converging to the wrong thing.
- **FORWARD-ADJOINT-K** — $F$ and $F^\dagger$ share $k$; adjoint fundamental strictly positive. Exercises `.T`/`.H` on the operator that matters.

*Tier 3 — a priori rate bounds (require kernel bounds; the "80% a priori" regime):*

- **BIRKHOFF-CONTRACTION-RATE** — finite projective diameter $\Delta$ of $F(K)$ ⟹ dominance ratio $d \le \tanh(\Delta/4)$, certified from cross sections and geometry before any run; measured per-iteration contraction never exceeds it.
- **SI-SPECTRAL-RADIUS / SI-MONOTONE** — unaccelerated SI from zero guess: pointwise monotone from below; $\rho_{\rm SI} = c$ exactly in the infinite homogeneous medium; DSA exempted by design.
- **COMPARISON-MONOTONE-K** — raise $\nu\Sigma_f$ pointwise anywhere ⟹ $k$ cannot decrease. Physics regression with theorem backing.

All of Tier 1–3 is gated on $\mathbb{F} = \mathbb{R}$ (§I.4).

## I.8 The space inventory, the two generating theorems, and the axis taxonomy *(new)*

**Continuous inventory — five kinds of space, two theorems generating the arrows.** Bulk data space $L^2(X)$, $X = D \times S^2 \times E$ (sources, emission densities, collision rates). Bulk **solution** space $W = \{\psi \in L^2 : \Omega\cdot\nabla\psi \in L^2\}$ — the graph space of streaming; flux lives in $W$, the things kernels produce live in $L^2$; the distinction is the single most under-implemented fact in transport codes. Trace spaces $\Gamma^\pm = L^2(\partial D \times S^2_\pm \times E;\ |\Omega\cdot n|\,d\gamma)$ — the weight is *part of the space*; the causal typing $\Gamma^+ \to \Gamma^-$ is forced by well-posedness. Moment spaces (images of the angular factor under harmonic analysis; metric = the moment mass matrix). Duals of all, paired by Green's identity. The generators: the **trace theorem** ($\gamma^\pm : W \to \Gamma^\pm$ bounded — traces exist for solutions, not data; you cannot trace an element of $L^2$) and **Green's identity** (the adjointness defect of streaming *is* a difference of trace pairings — duality lives on the boundary; the signed boundary form is a Krein structure whose canonical decomposition is the hemisphere split, which is what reconciles "$\Gamma^\pm$ are different spaces" with "two half-spaces").

**The fully typed fixed point:** $\psi = \Lambda(S\psi + q + k^{-1}F\psi,\ B\gamma^+\psi)$ with $S,F : L^2 \to L^2$ (factoring through moments / ℓ=0), $B : \Gamma^+ \to \Gamma^-$ (a contraction: strict for vacuum/albedo, *isometry* for reflective/periodic — the dissipativity underwriting well-posedness, runtime-checkable as **BC-CONTRACTION**), $\gamma^+ : W \to \Gamma^+$, the parity isometry $R : \Gamma^\pm \to \Gamma^\mp$ explicit and never implicit, and $W \hookrightarrow L^2$ feeding kernels. Every historical mis-composition is a violation of one of these arrows.

**Axis taxonomy — discretize per tensor factor.** An **Axis** carries four things: the **index set** (structured: cells, ordinates, retained $(\ell,m)$, groups, faces); the **metric** (volumes / quadrature weights / moment mass matrix / $|\Omega\cdot n|$-weighted face measure — the metric defines the pairing, hence `.H`, norms, projections); the **basis kind** Nodal | Modal (determines cone representability: `Coordinate` | `Checkable` | `None`); **symmetry data** (parity action on angular axes — a constructible check; if the ordinate set is not parity-closed, discrete $R$ does not exist without interpolation, and adjoint capability on that axis is *conditionally declared*, loudly). Concrete axes: `SpatialAxis`, `QuadratureAxis`, `HarmonicAxis`, `GroupAxis`. There is **no primitive trace axis** — traces are derived.

**A Space is an ordered tensor product of axes plus tags**, with derived constructors: `Trace±(V)` (defined **only** for `WithTrace` spaces), `HalfRange(A, n)`, `Moments(V)`/`Ordinates(V)` (the frame swap, quality invariant $B/A$), `Dual(V)`, and direct sums (partitions, §I.10). **The `WithTrace` tag is the discrete survival of $W$ vs $L^2$**: in finite dimensions the distinction resurfaces as representation content — a WithTrace space's realization includes face/interface data (DD's edge unknowns exist precisely for this); tracing a cell-average-only field requires a **reconstruction operator**, a modeling choice typed as a named arrow, never a canonical map. `Sψ` is a Data-space element — same axes, trace tag off — which makes "you cannot trace an emission density" a type-level fact.

**Where fields live** (the typing that sorts every "what space is this in" question): flux in $K \subset W$; sources/emissions in Data (cone-flagged); the **residual of the posed system in the codomain** $(L^2 \oplus \Gamma^-)$, two blocks (bulk emission mismatch, inflow mismatch), signed, never cone-typed; the **error in the domain** $W$, signed; $\Lambda$ transports residual to error, which is exactly why small-residual and small-error diverge by the sweep's conditioning. Convergence norms must be the codomain's $G$-weighted norms — a Euclidean `np.linalg.norm` on raw arrays is quietly the wrong space.

**Identity and layout doctrine.** Contraction legality is **space identity, not shape equality** ($\Gamma^+$ vs $\Gamma^-$; $N$ ordinates vs $N$ moments; half-range axes of different faces — all shape-equal, all illegal to compose, all invisible to einsum). A space's identity is the tuple of its axes' identities plus tags; an axis's identity is its structural content by canonical construction (interning at the generators). Metric differences imply space differences. **Layout stays out of identity** — axis order in memory, contiguity, batching are Realization concerns (the JAX trace/retrace boundary); spaces own the data operator construction consumes, never operator implementations. The `.T`/`.H` rule (labels/layout vs. metrics) is member (ii) of the tightness family (§I.3).

**Generators and the forgetful decomposition.** Mesh = index set + measure + topology + embedding. The Space needs exactly the first two — `SpatialAxis = forget(Mesh)`; topology feeds Realizations (sweep DAGs, DSA stencils), embedding feeds kernels (track/chord lengths). Symmetrically `QuadratureAxis = forget(Quadrature)` (weights + parity closure; directions remain generator data). The mesh *induces* the spatial discrete measure and is not reducible to it. The per-ordinate sweep DAG is **not** axis data — derived at Realization.

## I.9 The retract lattice: how spaces relate *(new)*

**No subtype hierarchy of spaces — a poset connected by retraction pairs.** Every canonical relation (outflow half of a face slot, scalar inside angular, ℓ=0 inside moments) has one shape: a projection $\pi : V \to V'$, an embedding $\iota : V' \to V$, with $\pi\circ\iota = \mathrm{id}$ and — under the right metrics — $\iota = \pi^H$ up to a scalar. $V'$ is a **retract** of $V$; $\iota\circ\pi$ is the orthogonal projection onto the embedded copy. Two generator families produce every retract: **restriction retracts** (index subsetting — $\Gamma_+$ inside $\Gamma(f)$; `TraceRestrictionOperator` is this) and **contraction retracts** (collapsing an axis against a positive weight — angular→scalar with $w$, angular→partial-current with $|\Omega\cdot n|w$, moments→ℓ=0).

**Scalar space is not a subspace of angular space** — it is a different space connected by a retraction; what *is* a subspace is $\iota(V_{\rm scal})$, the isotropic functions, and conflating the two is the source of the classical $4\pi$ bookkeeping errors. Keeping them typed distinct makes normalization constants *consequences of adjointness* rather than memorized conventions ("uniform in flux" vs "uniform in current" rebroadcast = which weight the embedding is adjoint to, stated in the type).

**Fission and the white BC are the same diagram**, differing in exactly one datum — the collapse weight:

$$F = \underbrace{\chi\,\iota}_{\text{embed}} \circ \underbrace{\nu\Sigma_f}_{\text{act on retract}} \circ \underbrace{\pi_w}_{\text{collapse}}, \qquad B_{\rm white} = \underbrace{\iota}_{\text{embed}} \circ \underbrace{\alpha}_{\text{albedo}} \circ \underbrace{\pi_{|\Omega\cdot n|w}}_{\text{collapse}}$$

Both are positive, rank-one-in-angle composites; the white family (vacuum $\alpha{=}0$, albedo $\alpha{<}1$, white $\alpha{=}1$) is one arrow family parameterized by an endomorphism of the intermediate scalar space. `Contract(axis, weight)` and its metric-adjoint `Embed` are minted together by one constructor **owned by the axis**; the composites above become one-liners, and the adjoint of every such composite is correct by construction — the ERR-039 mechanism, generalized.

**Cones map forward along positive arrows** (functorially): $\pi(K) \subseteq K'$ (onto, for strictly positive weights), $\iota(K') \subseteq K$; the whole white BC is cone-preserving; a contraction retract *upgrades* cone representability.

**Blocks are retract calculus** (bridge to §I.10): $A_{ij} = \pi_i \circ A \circ \iota_j$ — block extraction is conjugation by the partition's restriction retracts, which makes `CoupledOperator` grids *derivable* from any operator plus a partition rather than only hand-assembled.

**Open design point (flagged, not settled):** whether `Contract` pairs carry their **section identity** — this particular $(\pi,\iota)$ being *the* canonical retraction for that axis-weight pair — so the composition algebra can recognize $\pi\circ\iota = \mathrm{id}$ and $\iota\circ\pi$ = projection symbolically. It would let the dunder consolidation collapse composites, but it is the first step toward a rewrite system: a real scope decision, not an obvious win.

## I.10 Partitions: three layers, two owners *(new; finalizes the v1 ruling)*

| Layer | What it is | Owner | Stage |
|---|---|---|---|
| **Partition** | unordered decomposition; always **per-axis** (energy-GS partitions `GroupAxis`, ordinate-Jacobi partitions `QuadratureAxis`, **domain decomposition partitions `SpatialAxis`**, octant scheduling partitions `QuadratureAxis` by a symmetry-compatible coarsening); induces $V = \oplus_i$ by pulling through the tensor product | Space (axis) | 2 |
| **Block digraph + order** | which $A_{ij}$ are structurally nonzero (declared by the kernel: upscatter band, streaming stencil — not probed), SCC condensation, topological order/levels | computed from (operator, partition) | between 3 and 5 |
| **Schedule** | which blocks implicit, which lagged | Splitting | 5 |

The middle layer is the one v1 blurred: fast→thermal order exists only by the scattering kernel's triangularity; the sweep wavefront only by streaming's causality — the order is *operator-derived*, not space data. **KBA is the topological level sets of the spatial-partition block digraph per ordinate** — the same construction as energy-GS's group order; the machinery is variable-agnostic at the space level too. Reifying the middle layer as a `BlockDigraph` object is recommended (moderate confidence): $\alpha$ and time-dependent posings will recompute condensations under shifted diagonals and want to assert "the shift changed values, not the condensation" — the recipe/instance split, one level down.

## I.11 Quotients, the homogeneous solver, and where the modes live *(new)*

**The homogeneous solver is not measureless — it is maximally quotiented.** Translation invariance quotients the spatial axis to a **one-point axis with unit measure** (the "per unit volume" intensive convention *is* the normalized quotient measure); isotropy quotients angle to the ℓ=0 retract. The infinite-medium spectrum problem lives on `GroupAxis` alone tensored with two trivial axes — and everything transfers with **no special cases**: GroupAxis is nodal, so the coordinate cone, CW brackets, ray normalization, and the IRREDUCIBLE gate (group-graph connectivity) all apply verbatim. The buckling/B₁ family is the *intermediate* quotient — a one-dimensional **Modal** spatial axis parameterized by $B$, on which streaming is multiplication by $iB\mu$ — and Fourier convergence analysis is *the solver diagonalized over this quotient family*, computed fiberwise: the a priori regime of §I.5 is the fiber picture. Consequence for the plan: the `domain=None` escape hatch closes not with an ad-hoc 0-D space but with **this** quotient instance (Part V, S2/1.6).

**Where the eigenmodes live.** Mode 0: a ray in $P(K)$. Every higher mode: ambient $V$, or $V_\mathbb{C}$ when complex (conjugate pairs of the nonsymmetric pencil) — same `FunctionSpace`, different predicate values; the question dissolves once cone membership is an element predicate. Three architecture-relevant facts: biorthogonality $\langle\varphi^\dagger_i, \Pi\varphi_j\rangle = \delta_{ij}$ — mode coefficients are extracted through the **dual basis realized by adjoint modes through the pencil's $\Pi$** (tightness-family member (iii); wrong-weighting extraction of the dominance mode is the standard practitioner error); the spectral projector $P_0 = \varphi_0\langle\varphi^\dagger_0, \Pi\,\cdot\rangle$ is **not a positive operator** despite being built from two cone elements (the mechanism behind SIGN-STRUCTURE); the Hilbert-metric contraction story lives entirely on $P(K)$ and sees higher modes only through the dominance ratio.

---

# Part II — Corrections ledger

Everything falsified or materially sharpened. #1–14 from v1, retained; #15–20 new in v2. Several are load-bearing for the action plan.

| # | Claim | Status | Correction |
|---|---|---|---|
| 1 | "CFL is a support condition, not a capacity condition" | **wrong** | Bandwidth is a hard capacity; aliasing is saturation and drives *instability* nonlinearly (Orszag 2/3). Surviving remnant: CFL constrains the *ratio* $\nu = P_x/P_t$ of sampling densities, not their level. |
| 2 | "Implicit ⟹ $A^{-1}$ dense ⟹ one SCC" | **wrong** | Dense-*triangular* exists. Time-implicit $S_N$ adds $1/v\Delta t$ to the diagonal and stays exactly sweepable. Infinite reach ≠ cycles; cycles need bidirectional reach. |
| 3 | "Graph theory is the combinatorial shadow of Feynman–Kac" | **wrong** (corpus L177) | Under the killing split there is no stochastic process to average over; FK is vacuous there. Correct placement: the graph is the combinatorial content of the A1 axis under the killing split. |
| 4 | "Reflective/periodic BCs close tracks into cycles" | **wrong folklore** | `reflective\|vacuum` is acyclic; only opposing/periodic laws create SCCs. Acyclicity is an operator property; triangularity is an (operator, order) property. |
| 5 | "Elliptic ⟹ information destroyed (inverse ill-posed)" as one table cell | **imprecise** | Three objects conflated: forward solve (bounded), adjoint (bounded), source recovery (the only ill-posed one). |
| 6 | R9: `LossRepresentation` → `CausalFactor` | **retracted** | The object is being deleted, not renamed. Contents split by the answer/cost rule (§I.3). ⚠ 2026-08-19: the deletion never happened (see the III-A banner) — the retraction stands, its premise doesn't; the object lives on as the cost-side strategy Protocol. |
| 7 | Pipeline: Declaration → Splitting → Posing | **wrong order** | The splitting splits the *shifted* operator; posing precedes splitting-instantiation. Finalized in Part IV. |
| 8 | "$k$ pencil: shift not absorbable ⟹ sweep is lost" | **overstated** | The sweep is intact; $\sigma F$ rides as a gain. Multiplier shifts modify $M$ (free); integral shifts modify $N$ ($\rho$ penalty). Wielandt confirms. |
| 9 | "The resolvent" for $A_{\rm loss}^{-1}$ | **misnomer, verified live** | Reserve *resolvent* for the $\lambda$-family; the layer is `green`/`loss_inverse`. |
| 10 | `TrackSweepOperator` as MoC's sibling type | **must be rejected** | Same object (the causal/Volterra inverse — v2: the inverse of the *same bijection* $(T,\gamma^-)$, §I.2), different traversal. A traversal strategy on one type. |
| 11 | `LaplaceResolvent(T, A_prompt)` + `AlphaEigenproblem(T, A_prompt)` | **twin** | One object: `Pencil(A_prompt, T)`; the Laplace resolvent is `.at(s).inverse()`, the $\alpha$'s its poles. |
| 12 | `M` symbol | **four live meanings** | Reserve $M$ = splitting operator; production → $\Pi$; loss matvec → $(L{+}C)$; projection → $P_{\rm har}$. |
| 13 | `IntegralKernelOperator` | contradiction in terms — **⛔ CONFLICTED 2026-08-19** | split into `IntegralKernel` (data) and the bound arrow produced by `.on(V)`. ⛔ memo C: the tree ships the explicit counter-doctrine since 2026-06-26 — a `@runtime_checkable` Protocol ("a Kernel REFINES LinearOperator", gated by `tests/transport/test_integral_kernel_category.py`); operators EXPOSE `.kernel`, nothing binds one. Landing #13 as written reverses a ratified design — adjudicate (Part VI.5 D3), don't execute. |
| 14 | $\rho_{\rm SI} = c\arctan(\lambda)/\lambda$ in 3-D | **unverified** | Blocks 2.4's shipping, not its scaffolding. |
| **15** | **"Flux is a torsor in affine space"** (member-family doctrine; verbatim in `coupled_system.py` closing prose) | **overturned** | Flux has a canonical zero (vacuum) and superposition (linearity in a fixed medium); the no-adding intuition was fiber discipline (different problems), not affine structure. Correct: positive cone $K \subset V$; differences in ambient $V$; eigenfunctions rays in $P(K)$ (§I.7). The torsor prose must be swept from the member-family docs or the codebase teaches a retired ontology. ⚠ 2026-08-19 (memo D): "prose sweep" understates by an order of magnitude — the torsor is enforced arithmetic (`_flux_role.py`; [M] 16 production files); see the §I.7 banner + Part VI.5 D1. |
| **16** | **"Two copies of $\mathbb{R}^n$ are 'the same' space regardless of which inner product is installed"** (`FunctionSpace` docstring doctrine) | **overturned** | The metric defines `.H`; two spaces differing only in metric are spaces where the same symbol denotes different operators, and composing across them must be refused. Identity = axes' structural content + tags; metric differences imply space differences (§I.8, audit F1). |
| **17** | **Weak equality as a feature** ("two trace spaces on meshes of the same total boundary size compare equal"; `FullFieldSpace`'s "compare equal, **so** the OperatorSum guard accepts") | **reclassified as defect** | Equality was blunted to compensate for the lack of canonical construction. The cure is interning (one construction site per mathematical space), then strong structural identity — not weakened identity globally. The stringly name-discrimination (`_face_space_name`) is structural identity's job done by hand. |
| **18** | "The homogeneous solver has no spatial or angular measure" | **reframed** | It is the maximally quotiented member of a symmetry-reduction family: one-point spatial axis with *unit* measure (the intensive convention is the normalized quotient measure), ℓ=0 angular retract. No special cases in the space system or the cone battery (§I.11). |
| **19** | `SNMethodSpace` named as a space | **misnamed, kept** | It is the **generator record** — the bundle (mesh, quadrature, face, trace) axes are forgotten *from* (§I.8). Right job, wrong category claim; rename when touched. |
| **20** | v1 Phase 1 sequencing: kernel binding independent of the space layer | **amended** | `.on(V)` reads frames/measures from $V$'s axes; audit F1–F3 show the live space layer cannot answer correctly. The space refactor (Phase S) is a prerequisite of 1.1/1.4/1.6, not a parallel track (§I.3, Part V). ⛔ PARTIALLY REFUTED 2026-08-19 — see the §I.3 banner (operator-owned frames are ratified + production); survives only for space-data-reading items. Part VI.5 D2. |

---

# Part III — ORPHEUS state audit

## III-A. Operator machinery (v1 audit, retained)

**Already built and correct (do not re-plan):**

- The four-member inverse family with `InverseWrapMixin`, structure-keyed `.inverse()`, seeded-apply as structure (#285 closed), the object-identity involution, and the honest ordering rulings. #226 complete through step 6 and clean.
- `B.split(schedule)` → $(B_{\rm lower}, B_{\rm upper})$ and the reified splitting matrix `ScheduledInvertibleOperator` with machine-precision round-trip. **Half of the splitting stage exists**; $N$ is an unreified `*gains` bag — no object owns the pair whose $\rho$ governs convergence. ⚠ IMPRECISE (memo B, 2026-08-19): `WithinGroupSystem` (`sn/coupled_system.py:326`, pre-plan) already owns the (loss, implicit_operator, explicit_gains) triple as a passive RECORD; the live defect is R7 — `_select_si_splitting` re-derives the G-S pair instead of consuming the record (strict-xfail pinned). And `ScheduledInvertibleOperator` is user-ruled (2026-07-28) to DISAPPEAR with campaign P5 — retire, don't rename; build nothing on it.
- The frame hierarchy `FrameBase → PetrovGalerkinFrame → GalerkinFrame` (#268), `scattering.py` documenting $S = R\circ\Lambda\circ M$ via Funk–Hecke. The binding machinery of §I.3 is substantially present; missing: the kernel/operator split around it and the tightness gate.
- `sweep_acyclicity` with SCC computation and the acyclicity gate matrix (1-D).
- The default splitting recipe, implicitly (left-spine head = $M$).

**Known mis-factorings (slated or diagnosed):**

- `KEigenvalue` fuses Problem and Solver; the split into `CriticalityProblem` + `PowerIteration` is the natural landing site for the Pencil — build the pencil *as* the split.
- `power_iteration`'s `EigenvalueSolver` Protocol speaks fission vocabulary; neutralize at the split.
- `LossRepresentation`: deleted; contents split by the answer/cost rule; ERR-026 closes structurally.
  ⛔ STALE AT WRITE (memo C, 2026-08-19): never deleted — alive as the cost-side strategy Protocol (`sn/loss_representation/__init__.py:251`); the answer/cost split happened AROUND it (closure → `SNMesh.scheme` at mesh construction, traversal → strategy via `default_for(mesh)`); ERR-026 was already CLOSED 2026-06-12 (ERR-058 fix, 12 catcher files).
- ERR-039: an instance of the binding-correctness theorem; land the in-flight fix *with* the tightness gate.
  ⛔ STALE AT WRITE (memos B/C, 2026-08-19): fixed 2026-05-10, re-fixed 2026-05-26, merged — nothing was in flight on 2026-08-08. The class closed by TYPED SEPARATION (four separately-typed operators + the `g_C` space metric + generic `.H`), not a tightness gate; no tightness machinery exists at HEAD ([M] grep → 0 hits). Item 1.4 stands alone now (Part VI.3 R5).
- `domain=None` escape hatch: verified live; closes as the quotient instance (§I.11), forcing consumer = homogeneous solver.

**Guide drift (v3 vs. tree):** `BlockOperator` (v3) vs `CoupledOperator` (code) — v3's name is better; §29.7's twin (correction #11); §9's $R\Lambda M$-as-optimization framing; the `SweepOperator`/`TrackSweepOperator` region needs correction #10 written in.

## III-S. Space layer (new audit; read from the tree 2026-08-08)

**Genuinely right — keep, do not churn:**

- **The Γ ladder** is the best-designed part of the space layer: the tangential third class first-class ($\Gamma(f) \ne \Gamma_+ \oplus \Gamma_-$ in storage) with the quotient statement written down (ker $G$ = the tangential excess, so the direct sum holds exactly in the quotient); the partial-current metric $|\Omega\cdot n|\odot w$ installed unconditionally with a Moore–Penrose pseudo-inverse for its degeneracy; `omega_dot_n` as single source of truth; $\gamma^\pm$ bound to typed half-spaces refusing cross-face composition despite equal shapes.
- **View-G doctrine** (space = geometry; units/role on the Field) — consistent with "spaces own metrics, fields own quantities." Keep.
- `find_factor` by type; dual idempotence; the ⛔ post-mortem annotations.
- **`CoupledSpace`/`CoupledOperator`**: offsets-as-structure, member-wise metric dispatch, adjoint free through the generic `_AdjointOperator` wrapper (Mode-12 closed), typed blocks making wrong pairings unconstructable, missing block = structural zero. Its rejected-alternative note (padded `OperatorSum`) is retroactively justified as refusing partition-without-space-ownership.

**Structural flaws, severity-ordered:**

- **F1 (correctness) — nominal `(name, shape)` identity, everything else `compare=False`.** The base doctrine is correction #16; the codebase patches around it stringly (names encoding face/role by hand-minting discipline) and the discipline has documented leaks: cross-mesh/cross-quadrature aliasing stated as a feature (correction #17); $\Gamma_+(x_{\max})$ on `gauss_legendre(8)` equal to the same on any 8-point rule; `FullFieldSpace` equality weakened expressly so the composition guard accepts spaces built at different sites over the same mesh — the confession that canonical construction (interning) is the missing cure.
- **F2 (correctness; already bit) — no axis primitive; metric-to-axis binding by broadcast convention.** The ⛔ note documents the exact generic failure: `inner_product` (trailing broadcast) vs `apply_metric` (leading broadcast) applying *one space's own metric along two different axes* depending on the method called. Compounding: axis order is inconsistent across spaces (`angular_flux_space` = (cells, ordinates, groups); SH production layout = (L+1, 2L+1, ng, *spatial)) — every inter-space arrow carries an implicit transpose no type sees.
- **F3 (correctness-adjacent) — inconsistent metric doctrine in the bulk.** `scalar_flux_space` Euclidean with volumes "absorbed into the operator"; `FullFieldSpace` declares $G_{\rm bulk} = V_{\rm cell}w_n$; the trace layer refuses `None` metrics and spells `ones` deliberately. Three answers to "who owns the volume weight" in three places; adjoint correctness depends on remembering which operators absorbed volumes (ERR-067's breeding ground). The trace layer's doctrine is the correct one; a stage-inversion instance (measure = stage-2 datum living at stage 3).
- **F4 (structure/performance) — `TensorProductSpace` materializes the outer-product weight tensor densely.** O(state)-sized allocation per space storing a metric separable by construction, discarding exactly the per-factor structure the adjoint machinery wants.
- **F5 (missing capability):** no `WithTrace`/Data distinction; no cone/basis-kind metadata (no structural hook for the Tier-1 battery); no parity data (adjoint-capability-conditional-on-parity-closure undeclarable); no `GroupAxis` at all — energy is an anonymous integer in a shape tuple, so the V/V* collapse-consistency hook (tightness family (iv)) has nowhere to live; no ground-field parameter (blocks the noise/contour consumers).
- **Composition fragmentation:** three parallel product mechanisms — monolithic shape-tuple factories (the SN bulk path; never touches the product machinery), `TensorProductSpace` (moment-carrier path), `CoupledSpace` (the ⊕) — plus `SNMethodSpace` claiming the space category while being the generator record (correction #19).
- **Doctrine drift:** the torsor prose (correction #15) verbatim in `coupled_system.py`'s member-family paragraph.

---

# Part IV — The pipeline, finalized

The inversion of Posing before Instantiation stands (the splitting splits the *shifted* operator). v2 amends stage 2's parenthetical and the Posing preconditions.

```
0  DATA          cross sections, nuclear data, generators   CrossSectionField; Mesh, Quadrature,
                                                            GroupStructure = the generator objects
1  KERNEL        representation-free two-point data         ScatteringKernel, FissionKernel, Λ
2  REALIZATION   kernel + frame + measure → arrow           .on(V) → operator with fixed (domain, codomain)
                 AXES forgotten from generators; SPACES     Axis = forget(generator): (labels, measure,
                 minted as axis products + tags; their      kind, parity). V = ⊗ axes + tags (WithTrace,
                 PARTITIONS minted per-axis                 𝔽). Partitions per-axis → ⊕.
                 L's binding carries the CLOSURE            (answer-changing; ERR-026 closes here; the
                                                            closure is also the trace content — §I.3/I.8)
                 realization records cone-preservation      (DD: False; step/SC: True — §I.7)
3  DECLARATION   the operator sums                          A_loss = L + C − S − B;  A_prompt = A_loss − F
4  POSING        the spectral problem (optional)            Pencil(A, Π) + μ→physical map + Π-absorbability
                 + existence preconditions                  k/α asymmetry; IRREDUCIBLE gate; 𝔽 = ℂ
                                                            declared here for noise/contour consumers
5  STRATEGY      the splitting RECIPE (structural,          SplittingRecipe: implicit set, schedule,
                 shift-invariant; declared once)            block order (consumes the BlockDigraph
                                                            computed from operator + partition — §I.10),
                                                            close-vs-split per SCC
6  INSTANTIATION recipe applied to A − σΠ                   Splitting(M, *gains) with values;
                                                            multiplier shifts → into M; integral → into N
7  SOLUTION      iterate / invert / normalize               SourceIteration, PowerIteration, contour
```

**Dependency direction ≠ pipeline order at one point, deliberately:** `pencil.py` imports `splitting.py`, never the reverse. **Well-formedness invariant (assertable):** every splitting stack terminates in a causal or materializable core.

**The stage-inversion diagnostic** — every architectural defect found by either audit is later-stage machinery bound into an earlier stage, or an earlier-stage datum living downstream:

| Defect | Inversion |
|---|---|
| `LossRepresentation` | stage-6 traversal bound at stage 2 |
| $S/F/C$ apply-time dispatch | stage-2 realization deferred to stage 3 |
| `SweepOperator` named for traversal | stage-6 vocabulary on a stage-2/3 object |
| $A_{\rm loss}^{-1}$ called "the resolvent" | stage-4 vocabulary on a stage-6 object |
| `LaplaceResolvent` ≠ `AlphaEigenproblem` | one stage-4 object minted twice |
| closure differing between `apply` and `solve` | stage-2 choice made twice at stage 6/7 |
| **volumes "absorbed into operators" (F3)** | **stage-2 measure living at stage 3** |
| **metric excluded from space identity (F1)** | **stage-2 structure demoted to runtime convention** |
| **sweep DAG contemplated as axis data** | **stage-6 traversal artifact offered to stage 2 — refused by design (§I.8)** |

---

# Part V — Action plan, by implementation order

Ordered by dependency, not value. Slots into the current branch state (`refactor/sn-operator-algebra`; 9.5/9.6 pending; ERR-039 in flight). Each item names its gate. **v2 re-sequencing:** Phase S (space substrate) is inserted as a prerequisite of the remaining Phase-1 items (correction #20); Phase C (cone battery) runs parallel with an immediately-actionable subset; existing items gain the dependency edits marked **[v2]**.

> ⛔ **BRANCH-STATE ARCHAEOLOGY + A LIVE SIBLING (2026-08-19).** `refactor/sn-operator-algebra` is merged history and ERR-039 closed 2026-05. This Part's Phases 1–3 overlap the tracked, user-ruled `.claude/plans/operator_strategy_realization_campaign.md` (P0 merged; P1 "monomorphic domains" live, measured by 12 strict-xfails; [M] 105 passed / 21 xfailed reproduced 2026-08-19): its P2b/P2c own 1.1's binding-collapse territory, P3/P5 own 2.1/2.2/2.10 (PAIRED-CONSTRUCTION + SHAPE-SYMMETRY user rulings; "one interface, two constructors" rules the #296 fork), P4 owns 3.1/3.2 (`GeneralizedEigenPencil` + the L23a additive constraint), and O-4/O-5 precede any α work. Compose with it — do not mint a twin plan of record. Full ownership map: Part VI.2.

## Phase 0 — Now. Doc-only; prevents rework; no code risk.

| # | Item | Where | Gate |
|---|---|---|---|
| 0.1 | Pre-commit the no-`TrackSweepOperator` ruling (v2: cite the $(T,\gamma^-)$-bijection form, §I.2). | `sweep_operator.py` docstring + cross-domain-attacker memory | blocks the fork |
| 0.2 | Fix the `M` collision ($M$ = splitting operator; $\Pi$; $(L{+}C)$; $P_{\rm har}$). | docstrings, grand report §9, `notation.rst` | grep for bare-`M` in prose |
| 0.3 | Retire "resolvent" as the name of $A_{\rm loss}^{-1}$. | `eigenvalue.py`, `iteration.py`, `operator_inverse_family.rst` | the word appears only λ-parameterized |
| 0.4 | Collapse the `LaplaceResolvent`/`AlphaEigenproblem` twin to one `Pencil(A_prompt, T)`. | grand report §17–18, §29.7, §36, §42 | one object, one name |
| 0.5 | Record the finalized pipeline (Part IV as amended) superseding the five-stage order. | grand report §1/§30; `operator_algebra.rst` | stage-inversion table included, with the three v2 rows |
| 0.6 | $k/\alpha$ asymmetry as posing precondition. | posing table | stated before any α solver lands |
| **0.7** | **Retire the torsor doctrine (correction #15).** One paragraph: cone ontology, canonical zero, superposition-in-a-fixed-medium, differences in ambient $V$, rays in $P(K)$. Sweep the "affine flux torsor" prose from `coupled_system.py` and the member-family docs it cites. | `coupled_system.py` closing prose; member-family (`Composite`) docs; a new `docs/theory` cone paragraph (D7 carries the full treatment) | grep for "torsor" returns only the retirement note **⚠ re-scoped 2026-08-19: [M] 16 production + 16 test files + 5 doc pages ENFORCE the torsor — a type-system carve gated on the Part VI.5 D1 ruling, not a prose sweep** |
| **0.8** | **Record the identity doctrine (corrections #16/#17) as the Phase-S charter**: identity = axes' structural content + tags by canonical construction (interning); metric differences imply space differences; layout excluded; the weak-equality "features" reclassified. | `space.py` module docstring (replacing the same-regardless-of-inner-product paragraph); `operator_algebra.rst` | the old doctrine sentence gone; charter cited by S1–S3 |
| **0.9** | **Rename `SNMethodSpace` → generator-record vocabulary in prose** (correction #19; code rename rides the next touch, not a dedicated PR). | `sn/mesh/method_space.py` docstring | docstring no longer claims spacehood |

## Phase S — The space substrate. Prerequisite of remaining Phase-1 items; interleaves with in-flight ERR-039 work. Semantic fixes and mechanical moves kept in separate PRs throughout.

| # | Item | Depends on | Gate |
|---|---|---|---|
| S1 | **`Axis` primitive + factored metrics beneath the existing surface.** `Axis(labels, size, metric_1d, kind: Nodal\|Modal, parity: ParityMap\|None)`; `FunctionSpace` gains `axes` and derives `shape`; `inner_product`/`apply_metric`/`apply_inverse_metric` become per-axis contractions. **Deletes `_broadcast_metric`** (F2's failure class) and F4's densification (weights never outer-multiplied). Includes the ⊕ constructor from day one — `FullFieldSpace` and the Γ ladder already are direct sums; `CoupledSpace` becomes the unique ⊕ target. Old constructors wrap bare shapes in anonymous axes. **⚠ 2026-08-19 (memo A): the LIVE mints are `augmented_mesh.py:1099` + `_bases.py:381/:503`; the space.py factories (`angular_flux_space`/`scalar_flux_space`) are DEAD (test-only) — target the live sites or the step fixes dead code. Also rewire the two docstrings that TEACH the leading-broadcast convention (`angular_trace_space.py:588-593,825-829`).** | 0.8 | bit-identical results suite-wide; `_broadcast_metric` gone; a memory assertion that no dense outer weight tensor is allocated |
| S2 | **Rebuild the leaves on axes; resolve F3 toward the space.** `SpatialAxis.from_mesh` (owns $V_{\rm cell}$ — volumes move *out* of operators; sequence carefully, it moves weights across the operator/space boundary), `QuadratureAxis.from_quadrature` (owns $w$, parity closure), `HarmonicAxis` (metric per ℓ; retires the padded-rectangle metric), **`GroupAxis.from_library`** (finally a real object; carries the V/V* hook of tightness-family (iv) even if unexercised). One canonical mathematical axis order; layout freed for Realization. `angular_flux_space`/`SphericalHarmonicSpace` become thin factories over axis products — collapsing the three composition mechanisms to one. **1.6 lands here as the quotient instance** (§I.11): the 0-D energy-only space = trivial spatial ⊗ ℓ=0 ⊗ `GroupAxis`; the homogeneous solver is the forcing consumer; the `basis_shape=(ng,1)` double-pass retires. | S1 | volumes read from the space by every adjoint path (ERR-067 class); homogeneous path constructs the quotient space with no explicit `basis_shape`; JAX design constraint recorded: axis data static, values traced |
| S3 | **Flip identity to structural (F1).** Interning at the generators (the mesh already caches its trace space — extend the pattern); structural `__eq__`/`__hash__` on axis tuples + tags; `name` demoted to repr; the stringly `_face_space_name` discriminators become redundant. **Intermediate mode first**: structural equality computed and *logged* on mismatch while nominal equality governs — inventory the aliasing sites before the flip (surgical-edit discipline). | S1, S2 | logged-alias inventory empty, then the flip; cross-mesh / cross-quadrature / cross-face composition REDs; `FullFieldSpace`'s weak-equality rationale deleted |
| S4 | **Tags.** `WithTrace` (γ± defined only on tagged spaces; reconstruction Data→WithTrace as a named arrow — the "boundary flux from cell averages by nearest-cell copy" sin becomes a type error); **cone-rep per axis** (`Coordinate\|Checkable\|None` — the structural hook Phase C's tests key on); parity-conditional adjoint declaration (PARITY-CLOSURE loud at construction). Stress-test the WithTrace boolean against the DFEM plan — per-face polynomial trace content may need richer structure than a flag; this is the audit's weakest-confidence design point. | S2 | γ± on a Data-tagged space raises; a DD realization advertises trace content, an FV cell-average space does not; cone-rep flags drive Phase C test dispatch |
| S5 | **Ground field 𝔽 + conjugation $\mathcal{C}$** (§I.4). Complexification functor: same axes, sesquilinear pairing (`.H` conjugates), $\mathcal{C}$ as typed antilinear involution; cone predicates refuse $\mathbb{F} = \mathbb{C}$. | S1 | complexified `inner_product` conjugates (a deliberately complex test vector REDs the real path); reality-involution assertion available for 4.2 |
| S6 | **Retract constructors owned by the axis.** `Restrict`/`Embed` (index subsetting; subsumes `TraceRestrictionOperator`) and `Contract`/`Embed` (weighted collapse; $\iota = \pi^H$ up to scale, minted together). White BC, fission's angular factor, P₀ projection re-expressed as composites through retracts; normalization constants asserted as adjointness consequences. Section-identity recognition ($\pi\circ\iota = \mathrm{id}$ symbolically) explicitly **deferred** — recorded as the open rewrite-system scope decision (§I.9). **⚠ 2026-08-19 (memo A): not greenfield — the factored Lambertian π/ι chain (`IsotropicEmissionOperator @ PartialCurrentOperator`, named S(f) intermediate, per-link honest transposes) landed 2026-08-04 (G6.3); the delta is the axis-owned minted-together (π, ι) pair + the fission/P₀ re-spelling. Re-derive the gate against `sn/boundary/realizer.py`'s current white arm before sizing (§6c).** | S2 | white/albedo/vacuum realized as $\iota\circ\alpha\circ\pi$ with the adjoint correct by construction (an ERR-039-class negative control on the boundary path REDs without the pairing); partial-current normalization not hand-written anywhere |
| S7 | **Energy-axis V/V*** (tightness family (iv)): mark the identification $V \cong V^*$ on `GroupAxis` as carrying the collapse spectra; the group-collapse adjoint-consistency check becomes expressible. Lowest urgency; genuinely separable. | S2, S5 | a bilinearly-inconsistent collapse REDs the check on a two-group toy |

## Phase C — The cone battery (§I.7). Parallel workstream; the tests are the *consumers* that keep the ontology honest.

| # | Item | Depends on | Gate |
|---|---|---|---|
| C1 | **Immediately actionable subset** (no axis machinery needed; runs on raw arrays now): CONE-PRESERVE probes on sweep/S/F applications; COLLATZ-WIELANDT-BRACKET + CW-MONOTONE wired into the power-iteration history; SIGN-STRUCTURE at convergence; SI-MONOTONE on the unaccelerated path; MONOTONE-SOLVE probes. | — | DD documented non-cone-preserving (the flag honestly `False`) with step/SC `True` **(⚠ 2026-08-19: no Step/SC realization exists — only DD + LD ship, both `False`; a §6c gate-without-witness defect unless C1 lands a Step realization or re-scopes. The flag axis + DD's `False` + DD's negative-flux witness already ship — memo D §5. Part VI.3 R12)**; every regression-suite $k$ inside every bracket ever computed; a deliberate positivity bug (negative scattering moment) REDs CONE-PRESERVE |
| C2 | **Gates as preconditions**: IRREDUCIBLE (fission-coupling block-digraph connectivity, computed at Posing; decides which Tier-2 assertions are active) and M-MATRIX (small-problem sign pattern + large-problem probes). | C1; posing table (0.6) | a voided/decoupled toy problem deactivates Tier 2 with a logged reason instead of failing it |
| C3 | **Structural hooks**: cone-preservation flag on realizations (stage 2, Part IV); cone-rep-driven test dispatch (coordinate-wise on nodal axes, evaluation-based otherwise); $\mathbb{F} = \mathbb{R}$ gating. | S4, S5 | the flag asserted against measured behavior per realization; Tier tests skip with typed reason on $\mathbb{C}$ spaces |
| C4 | **Tier 3**: BIRKHOFF-CONTRACTION-RATE (projective diameter from kernel bounds; certified ceiling beside the Fourier estimate — the two logged together extend 2.4's heterogeneity metric); HILBERT-NONEXPANSIVE as the cheap per-apply form; COMPARISON-MONOTONE-K as physics regression. | C1, 2.3 | measured contraction ≤ $\tanh(\Delta/4)$ suite-wide; a νΣ_f-increase case with decreased $k$ REDs |

## Phase 1 — The kernel → operator refactor. **[v2] Now gated on Phase S for its space-reading items.**

| # | Item | Depends on | Gate |
|---|---|---|---|
| 1.1 | `.on(V)` binding for $S$ and $F$: kernel holds $\Lambda$; `.on(V)` supplies $M, R$ **from $V$'s axes** [v2: was "from $V$'s frame"; the axes are where the frame and measure now live]. **⚠ 2026-08-19: R∘Λ∘M is already production, frame-OWNED (`frame.conjugate`, memo C); the surviving apply-time dispatch is deliberate cross-method carrier arms (in-code keep ruling, `scattering.py:1172`). See the §I.3 banner + Part VI.5 D2 before executing.** | frame hierarchy (exists); **S1–S2** | apply-time `singledispatch` deleted from $S$, $F$ |
| 1.2 | $C$ stays a multiplier; carrier discipline only. **✅ already true at write (#261, 2026-06-26; memo C — zero frame imports).** | — | no frame import in `multiplication_operator.py` |
| 1.3 | $L$'s binding carries the closure; ERR-026 closes structurally. [v2: the closure decision is simultaneously the trace-content decision — coordinate with S4.] **⚠ 2026-08-19: ERR-026 closed 2026-06-12; the closure already binds at MESH construction (`SNMesh.scheme`); the trace-content conjecture is now MEASURED by #344 (ker A is a trace object, bulk share 1.1e-28; LD dim ker = 0). Residue: the ANGULAR closure's three spellings (#359). Re-scope against #344/#359 (memo C §3).** | ~~LossRepresentation dismantling~~ (see III-A banner); **S4** for the trace-content half | `apply`/`solve` closure-identity test replaces the ERR-026 narrative |
| 1.4 | Tightness gate on every binding; land with ERR-039, closing the class. [v2: the gate reads the frame from the space — S1/S2; and it is member (i) of the tightness family, documented as such.] **⚠ 2026-08-19: ERR-039 closed 2026-05 by typed separation — no gate anywhere ([M] 0 hits) and no branch to land with; this is a STANDALONE item whose "reads the frame from the space" clause inherits the D2 ruling.** | 1.1; ~~ERR-039 branch~~; **S1–S2** | a deliberately non-tight frame REDs (negative control) |
| 1.5 | Split `IntegralKernelOperator` → `IntegralKernel` / bound arrow. **⛔ conflicted with the ratified Kernel-refines-LinearOperator Protocol — see correction #13 (Part VI.5 D3).** | 1.1 | no kernel claims operatorhood |
| 1.6 | ~~Mint an ad-hoc 0-D `EnergyGroupSpace`~~ **[v2: superseded — lands as S2's quotient instance.]** | S2 | (see S2) |
| 1.7 | Cache kernels, not bound operators: one cross-section set → one kernel → many arrows. **⚠ 2026-08-19: intent satisfied structurally (every kernel is a read-through view on ONE `MaterialXSField`; S/F cached on the solver; L+C deliberately UNcached with an in-code anti-drift ruling, `solver.py:1323-1332`) — close as done-differently or re-scope (memo C).** | 1.1 | binding twice from one kernel shares $\Lambda$ storage |

## Phase 2 — Splitting reification. Independent of Phase 1; consumes Phase S's partition constructor when available.

| # | Item | Depends on | Gate |
|---|---|---|---|
| 2.1 | Reify `Splitting(M, *gains)` — give $N$ an owner. | — | `SourceIteration` constructed from a `Splitting` |
| 2.2 | Promote the implicit recipe to `SplittingRecipe` (structural, shift-invariant; carries schedule, block order, close-vs-split; JAX trace/retrace boundary — S2's recorded constraint applies here too). | 2.1 | zero-strategy path byte-stable |
| 2.3 | `Splitting.spectral_radius()` on the pair. | 2.1; 3.1-adapter | measured $\rho$ on a known-$c$ slab matches Fourier |
| 2.4 | A priori $\rho_{\rm SI}$ (verify the 3-D form first — correction #14); predicted/measured ratio logged as the heterogeneity metric. [v2: log the Birkhoff ceiling beside it — C4.] | diffusion solver; 2.3; C4 | prediction within tolerance homogeneous; ratio logged suite-wide |
| 2.5 | SCC closing-edge rank in `sweep_acyclicity`; close-vs-split per component. | — | §I.6 rank table reproduced |
| 2.6 | `MonodromyOperator` (promote `TrackMonodromy`): white/albedo and 1-D periodic closed exactly. **⛔ 2026-08-19: `TrackMonodromy` never existed in code (grand-report concept only); the real item is #300 — Woodbury close of the rank-1 boundary SCC, blocked on #299 (Part VI.3 R10).** | 2.5 | 1-D periodic solves with zero outer boundary iterations |
| 2.7 | Mesh diagnostics at setup ($\max\tau$ beside the acyclicity certificate; optional critical-path depth). | — | $\tau>2$ mesh warns |
| 2.8 | 2-D quarter-symmetry acyclicity gate (4-orthant DAG, one pass). | — | solver takes the one-pass path |
| 2.9 | CP conditioning fingerprint (cond($P$) grows, cond($I-cP$) bounded). | CP path | refinement-sweep test with both curves |
| **2.10** | **`partition_along(axis, blocks)` on spaces** → `CoupledSpace` whose members share the untouched axes; **block extraction by retract conjugation** $A_{ij} = \pi_i A\iota_j$ — `CoupledOperator` grids derivable from (operator, partition), not only hand-assembled. Energy-GS, ordinate-Jacobi, and domain decomposition all constructed through this one door. | S1 (⊕), S6 | the CP Jacobi/GS pair re-expressed as two schedules over one `partition_along(GroupAxis)`; a hand-assembled grid and the derived grid agree block-wise |
| **2.11** | **Reify `BlockDigraph`** (the middle layer of §I.10): structural nonzeros declared by the kernel, SCC condensation, topological order/levels. `SplittingRecipe` consumes it; KBA = its level sets on the spatial partition. Assert shift-invariance of the condensation under `.at(σ)` diagonal rebinds. | 2.2, 2.10 | energy order derived (not hand-coded) from the upscatter band on a WIMS-structure toy; the shift-invariance assertion REDs on a contrived structure-changing shift |

## Phase 3 — The pencil. Lands as the `KEigenvalue` split.

| # | Item | Depends on | Gate |
|---|---|---|---|
| 3.1 | Split `KEigenvalue` → `CriticalityProblem` (= `Pencil(A_loss, F)` + $\mu\mapsto k$) + `PowerIteration`; neutralize `EigenvalueSolver` vocabulary. [v2: wire CW brackets (C1) into the neutral solver boundary here — the bracket is production/state vocabulary, not fission vocabulary.] **⛔ means refuted 2026-08-19 — L23a: the pencil is ADDITIVE, `power_iteration`'s late binding stays; the neutral seam now exists as `measure_stopping_criteria`/`PowerIterationOutcome` (#340, post-plan). Part VI.3 R9.** | 2.1–2.2; C1 | one power loop preserved; SN $k$ byte-stable; brackets emitted per iterate |
| 3.2 | `Pencil(A, Π).at(σ)` with the absorbability predicate. **[v2: for $\sigma \in \mathbb{C}$, requires S5; the posing declares 𝔽.]** | 3.1; **S5** for complex shifts | $k$ path = `at(0)` bit-stable; Wielandt shift expressible as a gain |
| 3.3 | α-aware guard branch in `StreamingCollisionOperator.__init__`. | 3.2 | shifted construction beyond $v\Sigma_t$ raises the α-branch text |
| 3.4 | α shift-invert through the same machinery; existence precondition (0.6) enforced; **[v2: IRREDUCIBLE gate (C2) consulted at posing]**. | 3.2, 3.3, 0.6, C2 | rod-drop-style case; $k(\alpha)$ equivalence curve runnable |
| 3.5 | Well-formedness assertion (Part IV invariant). | 2.2, 3.2 | a cycle-complete recipe REDs |

## Phase 4 — Consumers unlocked by the pencil.

| # | Item | Depends on |
|---|---|---|
| 4.1 | Higher harmonics via Riesz projection — many `at(λ)`, one recipe. **[v2: complex contour ⟹ S5; mode extraction through the $\Pi$-weighted dual pairing (§I.11), asserted biorthogonal.]** | 3.2, S5 |
| 4.2 | Noise transfer function $R(i\omega)$; point-kinetics zero-power TF as adjoint-weighted rank-1 reduction (verification target). **[v2: S5 hard dependency; the reality involution $\delta\varphi(-\omega) = \mathcal{C}\delta\varphi(\omega)$ is a named assertion on the solve.]** | 3.2, **S5** |
| 4.3 | `CausalOperator` base + `causal_order` at MoC arrival; `SweepOperator` → `CausalInverseOperator` then and only then. | MoC |
| 4.4 | JAX backend alignment: recipe/instance = trace/retrace; **[v2: axis-data-static/values-traced recorded at S2 is the same constraint one level down]**. | 2.2, S2 |

## Doc workstream — parallel to all phases; each item cites its Part I section.

| # | Item | Target |
|---|---|---|
| D1 | Volterra/Fredholm classification; generator-vs-propagator; $\rho({\rm sweep}) = 0$; [v2: + the $(T,\gamma^-)$-bijection statement of the sweep, §I.2] | `operator_inverse_family.rst` |
| D2 | "From condition to rate": Perron–Frobenius strictness; iteration count = mean collision count; a priori estimator; inconsistent-DSA honesty clause; [v2: + CW brackets and the Birkhoff ceiling as the theorem-backed diagnostics, §I.5/I.7] | `path_integral.rst` §4 |
| D3 | Naming-discipline block for **adjoint** — four objects with coincidence conditions; [v2: + the tightness family (i)–(iv) as one documented pattern, §I.3] | `operator_adjoint.rst` |
| D4 | Velocity averaging as the Krein–Rutman compactness hypothesis (verify transport hypotheses before writing). | eigenvalue theory page |
| D5 | The answer/cost rule and the fiber/transversal reading; ERR-026 as closing witness; [v2: + closure = trace content, §I.3/I.8] | `loss_representation.rst` |
| D6 | DSA as an A2 value between SI and Case — route through `proposal-challenger` first. | `path_integral.rst` axis table |
| **D7** | **Flux ontology chapter**: cone, canonical zero, superposition, rays in $P(K)$, cone-membership-as-predicate, the torsor retirement narrative (correction #15), and the battery's epistemic gradient as the reader's map of what each passing test proves. | new `docs/theory/foundations/cone.rst` |
| **D8** | **The space taxonomy**: continuous inventory, the two generating theorems, $W$ vs $L^2$ and the `WithTrace` residue, axes and the forgetful decomposition (Mesh/Quadrature as generators), identity doctrine, where fields live (flux/source/residual/error). | new `docs/theory/foundations/spaces.rst` |
| **D9** | **The retract lattice**: $(\pi,\iota)$ pairs, fission ≡ white BC, normalization-from-adjointness, cone functoriality, blocks as retract conjugation. | `spaces.rst` §2 or standalone |
| **D10** | **The quotient family**: homogeneous solver as maximal quotient, buckling as the Modal intermediate, Fourier analysis as fiberwise diagonalization (correction #18). | `spaces.rst` §3; cross-cite `path_integral.rst` |
| **D11** | **Modes and the dual pairing**: biorthogonality through $\Pi$, the non-positive spectral projector, mode extraction as tightness-family (iii), complexification and the reality involution. | eigenvalue theory page; `operator_adjoint.rst` cross-cite |

## Research tier — recorded, not scheduled.

- Angular-entropy diagnostic $D_{\rm KL}(\psi/\phi\,\|\,1/4\pi)$ paired with the tightness $B/A$ check (representation error vs. method-class error).
- Transfer entropy as measured edge weights for coupled-multiphysics splitting decisions (activates with TITAN/HEST).
- Max-flow/min-cut on contributon current (speculative; no theorem known).
- **[v2]** Section-identity / rewrite-system for retract composites ($\pi\circ\iota \to \mathrm{id}$ symbolic) — the deferred S6 scope decision; revisit when the dunder consolidation (9.6 successor) has stabilized.
- **[v2]** Richer-than-boolean `WithTrace` for DFEM per-face polynomial trace content — S4's flagged stress-test; promote to a plan item when the DFEM campaign is scheduled.

---

# Confidence ledger

**High — derived in-conversation, or verified against the live tree:**
- All v1 high-confidence entries (splitting recursion + termination; Volterra/Fredholm classification; quasinilpotency mechanics; Fredholm status of the preconditioned loss and CP solve; Schwartz-kernel diagnosis; $R\Lambda M$ as binding; tightness as binding-correctness with ERR-039 as instance; shift-absorption rule; α guard identity; $k$-exists/α-may-not; PF strictness; conditioning three-regime split; 1-D periodic closure; quarter-symmetry DAG; pipeline inversion; import direction pencil → splitting)
- Cone ontology and the torsor overturn (§I.7; canonical zero + superposition are decisive)
- The sweep as inverse of the $(T,\gamma^-)$ bijection; the continuous space inventory and both generating theorems (Dautray–Lions/Agoshkov-standard)
- Collatz–Wielandt brackets and their monotonicity under primitivity; Krein–Rutman sign structure; Hilbert-metric nonexpansiveness; the Birkhoff bound's *form* $d \le \tanh(\Delta/4)$
- DD not cone-preserving; step/SC cone-preserving; cone-preservation as a realization property
- Retract structure ($\iota = \pi^H$ up to scale under the right metrics); fission ≡ white BC factoring; cone functoriality along positive arrows
- Partition three-layer split and its ownership; blocks as retract conjugation; KBA as topological levels
- The quotient reading of the homogeneous solver and buckling family
- Complexification requirements for noise/contour consumers; cone non-survival on $\mathbb{C}$; the reality involution as the physical constraint
- The full space-layer audit F1–F5 (each read from the tree 2026-08-08, with the ⛔ notes and weak-equality docstrings as primary evidence)
- The `.T`-labels / `.H`-metrics rule and its coincidence condition

**Moderate — verify before the dependent item ships:**
- $\rho_{\rm SI} = c\arctan(\lambda)/\lambda$ in 3-D (blocks 2.4's shipping)
- Velocity-averaging hypotheses for transport specifically (blocks D4's text)
- $0.2247c$ attribution detail (D6)
- The guard locating the α continuum on bounded domains (3.3's message must not claim it)
- Point-kinetics TF as exactly the adjoint-weighted rank-1 reduction (4.2's target)
- **`WithTrace` as a boolean being sufficient** — the audit's weakest structural bet; DFEM per-face polynomial content may force richer structure (S4 stress-test)
- **`BlockDigraph` reification** (2.11) over keeping the order inside `SplittingRecipe` — recommended for the shift-invariance assertion, but the seam is defensible until Phase 3
- Group-collapse adjoint inconsistency as *exactly* the bilinear-consistency condition stated (S7's toy gate is the verification)

**Judgement calls, flagged as such:**
- Interning + content-fallback for axis identity, over pure content hashing (serialization and any multiprocess JAX story carry the interning table; revisit at 4.4)
- ⊕ constructor in scope for S1 (widens the first PR; justified because `FullFieldSpace` and the Γ ladder already are direct sums)
- F3 resolved toward the space now rather than documenting the Euclidean-bulk convention (the JAX/autodiff argument decides it)
- Eager-logged-then-flip for S3 rather than a hard flip (fits the surgical-edit discipline; costs one intermediate mode)
- Section-identity recognition deferred to research tier rather than built at S6 (rewrite-system scope)
- Keeping `StreamingCollisionOperator`; deferring the `SweepOperator` rename to the MoC relocation; `SplittingRecipe` as object over an `inverse(strategy=)` seam — all as in v1

---

# Part VI — Grounding against the tree (2026-08-19)

The plan-authoring §7 reconciliation of this report against **HEAD `7aae9bf1` (main)**, run before any campaign work. Method: four parallel explorer dispatches; their memos carry the full `[M]` evidence and are the authority under this summary:

- `scratch/omr_v2_grounding/A_space_layer.md` — Part III-S, §I.8–I.10, Phase S
- `scratch/omr_v2_grounding/B_operator_machinery.md` — Part III-A, Phases 2–3, §I.1/I.4/I.6/I.10
- `scratch/omr_v2_grounding/C_kernel_binding.md` — §I.3, corrections #13/#20, Phase 1
- `scratch/omr_v2_grounding/D_flux_ontology.md` — §I.7/correction #15 blast radius, Phase C substrate

⚠ This plan file AND the four memos are **untracked** ([M] `git ls-files --error-unmatch` fails on this file) — commit them before campaign work starts, or the grounding is one `git clean` from gone.

## VI.1 The shape of the verdict

1. **Part I's theory stands where it is mathematics.** Where it makes TREE claims it splits three ways: §I.3's diagnosis describes the **pre-#261 tree** (the "(v1 audit, retained)" text was never §7-reconciled; the binding machinery landed 2026-06-24..26, six weeks before this plan); §I.7's overturn is a **pending ruling** whose object is enforced runtime architecture (banner at §I.7); §I.8–I.11 remain accurate as target-state design, with F1–F5 all measured live.
2. **Part III splits by half.** III-S is CURRENT — the space layer took [M] exactly 2 docstring-only commits since the audit — but INCOMPLETE: three space files were never read (`radial_characteristic_space.py`, `scalar_trace_space.py`, `spatial_moment_space.py`, all pre-dating the audit), which widens F1 (+2 weak-equality sites) and F3 (4–5 metric doctrines over 6+ sites, the space-owned side GROWING). III-A/§I.3 was STALE AT WRITE (banners in place).
3. **Part V's Phases 1–3 substantially overlap the live operator/splitting/realization campaign** (tracked; P0 merged; P1 live, 12 strict-xfails; user rulings: PAIRED CONSTRUCTION, SHAPE SYMMETRY, L23a, "one interface two constructors", "ScheduledInvertibleOperator disappears with P5"). The genuinely NEW v2 contributions: **Phase S**, **Phase C**, **complexification (S5)**, **the tightness gate (1.4)**, **the posing preconditions (0.6/C2)**, the **doc workstream**, and the theory spine.

## VI.2 Ownership map — who owns each Part-V item at HEAD

| v2 item | owner / state at HEAD | disposition |
|---|---|---|
| 0.1–0.9 | ALL unexecuted. 0.1: `TrackSweepOperator` 0 hits anywhere (fork unhappened, ruling unwritten). 0.3: "resolvent"-for-A⁻¹F live in `eigenvalue.py` module prose + `direct_eigenvalue` + `KEigenvalue` docstrings + `convergence.py` ("LU resolvent") — wider than the item states. 0.4: the twin exists as `AlphaEigenproblem` RESERVED in `docs/architecture/layering.rst:127,155,205` (a third naming surface — see D5). 0.7: re-scoped (banner). 0.8: doctrine sentence verbatim at `space.py:145-150`. 0.9: `method_space.py:1,:73` unchanged. 0.2: the tree adds a 5th M (`MultiplicationOperator`) | archivist campaign, after D1/D2 rulings |
| 1.1 | R∘Λ∘M binding LANDED pre-plan (frame-owned, `frame.conjugate`, `cc6d022e`); sibling campaign P2b/P2c own the realizer/dispatch-collapse | residual = carrier-arm collapse (P2c) + D2 ruling |
| 1.2 | ✅ satisfied at write (#261) | close |
| 1.3 | closure → `SNMesh.scheme` (landed); #344 gauge family MEASURES the closure↔trace-content conjecture; residue = #359 (angular closure's 3 spellings) | re-scope against #344/#359 |
| 1.4 | NOBODY — tightness machinery absent ([M] 0 hits); ERR-039 closed by typed separation instead | live, standalone; the one Phase-1 item that is a genuine gap |
| 1.5 | ⛔ conflicted (ratified Kernel-refines-LinearOperator Protocol + category test) | D3 ruling |
| 1.6 | superseded → S2 quotient instance (unchanged) | rides Phase S |
| 1.7 | intent satisfied structurally; L+C deliberately uncached (in-code anti-drift ruling) | close as done-differently, or re-scope |
| 2.1/2.2 | sibling P3 + P5 + PAIRED-CONSTRUCTION ruling; `WithinGroupSystem` owns the RECORD; live defect = R7 (strict-xfail) | fold into sibling campaign |
| 2.3 | sibling P3.3 ("M⁻¹N formable ⟹ ρ probeable"); ⚠ #341: boundary G-S is NOT a regular splitting — no Varga comparison theorems; a per-face closure-damping `spectral_radius` (post-plan, `5def63b0`) exists and is NOT ρ(M⁻¹N) — do not conflate | fold; carry #341/#343/#373 |
| 2.4 | nobody; correction #14 (3-D arctan form) still unverified and still blocks | live |
| 2.5 | `sweep_acyclicity` = derivations oracle, 0 production consumers (#324); rank computed nowhere; the wiring decision (#320) precedes "one addition" | live; wire-then-extend |
| 2.6 | ⛔ premise misnamed; the item is #300 (blocked on #299) | rename to #300 |
| 2.7–2.9 | unaudited this pass | re-derive at pickup |
| 2.10 | sibling P5.1/P5.2 (rules the #296 fork); `DiscreteMeasure.partition_by` (2026-05-10) is the landed measure-level partition sibling | fold |
| 2.11 | nobody; sibling P5.3 nearest ("Schedule = SCC of the implicit set") but names no BlockDigraph object | design decision inside P5 |
| 3.1 | ⛔ means refuted (L23a); the neutral seam = `measure_stopping_criteria`/`PowerIterationOutcome`/`IterationRecord` (#340, post-plan) | fold into sibling P4 |
| 3.2 | sibling P4 (`GeneralizedEigenPencil`, gated on P1); complex σ gated on S5 | fold |
| 3.3 | guard premise confirmed (`streaming.py:597`); new context: `StreamingCollisionOperator` is now an `OperatorSum` subclass with a mesh-identity invariant — a shift rebind must preserve BOTH guards | fold into P4/O-4 |
| 3.4 | sibling O-4 (independent α reference FIRST, before F migrates; needs a velocity fixture — the data layer carries no `v_g`; ⛔⛔ never resolve α by grepping `alpha`) + O-5 | fold |
| 3.5 | #320's `assert_cycles_are_declared` named-and-unbuilt | live |
| 4.1/4.2 | unowned; hard-gated on S5 (complexification) | future |
| 4.3/4.4 | unchanged (MoC-gated / JAX-gated) | future |
| Phase S | **NEW — the report's largest net-new campaign.** All target defects confirmed live at HEAD (F1 measured by live aliasing probes; F4 measured at state-size; F5 5/5 absent). #330 = the doctrine's issue-home; #369 = the measure-side identity twin (cite from S3); #295/#297 = adjacent scope decisions; `FunctionSpace` is now carrier-generic (PEP-696) — the seam S1 builds beneath | charter it, with memo A's corrected site lists |
| Phase C | **NEW — substrate partially ships**: flag axis + DD `False` + negative-flux witness; ψ≥0 realizer refusal; coefficient cone-as-predicate battery; ray normalization; named-criteria trajectories (`power_iteration`) for CW brackets; default-unaccelerated SI for SI-MONOTONE. Missing: any `True`-flag witness (R12), everything Tier-2/3, cone-rep tags (S4) | C1's value-level subset is available NOW and ontology-neutral; the typing half waits on D1. ⚠ name the battery to avoid `tests/sn/sweep/core/test_phase_c_gates.py` (#168's unrelated "Phase C") |
| D1–D11 | unstarted (no `cone.rst`, no `spaces.rst`; `operator_inverse_family.rst` carries no bijection statement) | archivist campaign after rulings |

## VI.3 Refutations and corrections (the banners' index)

- **R1** §I.3 "S/F/C have no unique (domain, codomain)" — CHANGED: declared **Optional** pair since ≤2026-06-26. Honest residue: (a) Optional-ness (= sibling P1's whole job), (b) multi-arm carrier dispatch, kept deliberately.
- **R2** "CollisionOperator" — did not exist at write (collapsed into `MultiplicationOperator`, #261).
- **R3** correction #20 — refuted for the angular binding (`frame-eigenbasis-ownership`); survives for space-data-reading items only.
- **R4** correction #13 — conflicted with the ratified Protocol counter-doctrine.
- **R5** "ERR-039 in flight; land 1.4 with it" — closed 2026-05-10/26 by typed separation; 1.4 stands alone.
- **R6** "ERR-026 closes structurally (via 1.3)" — already closed 2026-06-12 (ERR-058 fix; 12 catcher files).
- **R7** "LossRepresentation: deleted" — alive as the cost-side strategy Protocol; the answer/cost split happened AROUND it.
- **R8** "no object owns the (M, N) pair" — `WithinGroupSystem` owns the RECORD (pre-plan); the defect is R7-the-xfail (the driver re-derives).
- **R9** 3.1's mechanism — refuted by L23a (pencil additive; late binding stays).
- **R10** 2.6 — `TrackMonodromy` never existed; the item is #300.
- **R11** §I.6 "rank is one addition" — the SCC object is an unwired derivations oracle; the #320 wiring decision precedes.
- **R12** C1's step/SC leg — no witness exists ([M] registry = {diamond_difference, linear_discontinuous}, both `False`; `class Step` is a docstring example). §6c: fuse with a Step realization or re-scope.
- **R13** III-S site corrections — named monolithic factories are DEAD (live mints elsewhere); F2's two-axes divergence is HISTORY (fixed 2026-08-04; the live defect is the untyped convention); F3 under-counts in the plan's own favor.
- **R14** correction #15's scope — a type-system carve pending D1 (banner at §I.7). Measured internal tensions of the incumbent (decision material, not verdicts): scalar scaling deliberately legal; zero fluxes freely constructed; DSA's own "the swept vector IS the displacement from zero"; `affine_combination` has 0 production callers.
- **R15** correction #6's premise ("being deleted") — inherits R7; the retraction stands, its premise doesn't.

## VI.4 What the plan does not know (post-2026-08-08 or unread, load-bearing)

1. **#340 convergence rebuild** (2026-08-09..13): `measure_stopping_criteria` replaced `converged()`; `PowerIterationOutcome` + `IterationRecord` tree + typed `IterationBudget` (ERR-079); SourceIteration's stop is the ρ-honest equation residual; the DSA `corrector` seam. Every Phase-2/3/C design wires here, not into the old `(keff, history, flux)` triple.
2. **#344 gauge family** (2026-08-14/15): ker A in closed form; `LossKernelGauge` — a boundary-trace projector built by `frame.conjugate` (the plan's own binding machinery, now load-bearing); `face_transmission_spectrum` (closure ASKED, never tabulated); gauge wired at every solver exit. Measures §I.3's closing conjecture.
3. **Splitting theory post-plan**: #341 (not a regular splitting), #343 (octant order = unowned 2.5× rate lever), #373 (M_GS⁻¹ is a SUBSPACE inverse, full-space defect 3.25e-01).
4. **The assembly axis** (`SupportsAssembly`/`assemble`, sparse sums — pre-plan, absent from III-A): the third structural axis any monodromy materialization routes through.
5. **Space layer**: three unread files; carrier-generic `FunctionSpace`; the LD moment-mass inline metric product at the bulk mint (`augmented_mesh.py:1095-1098`) — the per-factor metric structure S1 wants, already spelled once at a site the F-rows don't cite (memo A §5).
6. **Q6 measure metadata** (`ReferenceMeasure`/`ExactnessClaim`/`invariance_group`/`half_range_clean`): the generator side of `QuadratureAxis = forget(Quadrature)` is now rich — S2's angular job narrows to FORGETTING. Q6 also minted #369, the F1-class defect on measures.
7. **Sibling-campaign facts**: P2a's relocation already happened (S/F/C live in `orpheus/transport/operators/`); the P1 xfail ledger enumerated in memo B §8.8.

## VI.5 Open decisions (the user's, in dependency order)

- **D1 — FLUX ONTOLOGY** (gates 0.7, D7, Phase C's typing half, #331). Torsor (incumbent: enforced, quiescent since 2026-07-26, blocks no live use case — superposition is source-level everywhere) vs cone (§I.7's argument: canonical zero + superposition-in-a-fixed-medium). The live demand is operator-domain: #331 (L accepts displacements, S/B refuse) — either ruling on #331 partially rules the ontology. Options: keep torsor + add cone PREDICATES (most of Phase C is typing-neutral); full cone overturn (blast radius memo D §9 — includes 4 `:eq:`-cited labels and the coding-elegance skill's anti-pattern #18, which cites FluxDisplacement as a positive precedent); or a scoped hybrid ruled at #331. **✅ RULED 2026-08-19 (user criterion: correctness; adjudicated same day): flux lives in the positive cone K ⊂ V — overturn and implement. Cone membership is an element PREDICATE, not a constructor invariant (DD does not preserve K, so a ψ≥0 type would refuse production output); cone-preservation is a realization flag (shipped); the torsor's iterate-hygiene content RELOCATES to the iteration layer (ρ, ‖Δψ‖/(1−ρ) on `IterationRecord`/`SourceIteration`); the fiber discipline (different problems don't mix) stays as space/mesh identity. Executed as phase CS3 of `space_and_kernel_binding_campaign.md`; #331 resolves to "operators are linear on V".**
- **D2 — FRAME/MEASURE OWNERSHIP** (gates 1.1/1.4's spelling, S2's scope). Operator-owned eigenbasis frames are RATIFIED + production; §I.3-v2/S2 propose axis-owned. Adjudicate what Phase S owns: identity + metrics + partitions (uncontested) vs frames (contested). **✅ RULED 2026-08-19: the frame is MINTED AT BINDING — the kernel carries its eigenbasis STRUCTURE (Funk–Hecke), the space's angular axis carries the MEASURE, and the constructor builds the frame from both. Preserves `frame-eigenbasis-ownership`'s content (the frame is the operator's eigenbasis, not a phase-space possession) while Phase S owns identity + metrics + partitions. Frames are constructed, not stored, on spaces.**
- **D3 — KERNEL DOCTRINE** (gates 1.5): "Kernel REFINES LinearOperator" (ratified) vs correction #13's data/arrow split. **✅ RULED 2026-08-19 (user): correction #13's direction — `ScatteringOperator` becomes `ScatteringKernel` (data) and a CONSTRUCTOR binds Kernel × Frame → the fully-bound operator. The "Kernel REFINES LinearOperator" doctrine is overturned; its Protocol + category test re-scope to the new split (campaign phase CS4).**
- **D4 — PLAN-OF-RECORD COMPOSITION + SEQUENCING**: fold Phases 1–3 into the operator-strategy campaign (keeps its P1→P6 order and rulings; v2 grafts = tightness gate, S5 complexification for complex shifts, preconditions 0.6/C2, BlockDigraph question at P5) and charter Phase S + Phase C as the new campaigns. Sequencing recommendation: **P1 first** (small, unblocked, 12-xfail-measured, surgical-carve posture), then Phase S — correction #20's space-first argument fell with R3, and P1's mandatory arrows give S3's identity flip typed consumers everywhere. User to confirm. **✅ RULED 2026-08-19 (user) — SPACE FIRST, overriding the P1-first recommendation: Campaign 1 = `space_and_kernel_binding_campaign.md` (CS1 Energy space, forced by the homogeneous solver → CS2 Spatial + Angular, SN as their composition → CS3 the cone carve → CS4 kernel→operator binding; the sibling's 12 P1 xfails survive as the ledger). Campaign 2 = the LossRepresentation overturn ("an early decision on partitioning" made LATE) → `GeneralizedEigenPencil` + resolvent machinery with complexification + partitioning — chartered after Campaign 1 lands.**
- **D5 — PENCIL NAMING**: `Pencil` (this report) / `GeneralizedEigenPencil` (sibling P4) / `CriticalityProblem`+`AlphaEigenproblem` (reserved in `docs/architecture/layering.rst`) — one name before P4, per [[feedback-high-signal-names]]. **✅ RULED 2026-08-19 (user, by usage): `GeneralizedEigenPencil`. Reconcile `layering.rst`'s reserved names when Campaign 2 lands it.**
- **D6 — C1's step/SC witness** (R12): land a Step realization with C1, or re-scope the gate to DD/LD + flag-axis assertions. **(Still open — decided inside the Phase-C battery's planning, which now follows Campaign 1's CS3 cone carve.)**
