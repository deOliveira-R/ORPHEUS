# ORPHEUS documentation — corpus architecture proposal

**Status:** PROPOSAL for review (2026-07-14). Supersedes nothing; extends issue #231
(whose page-template design is settled) to the **whole-corpus** layer #231 never addressed.
**Evidence:** full 71-page inventory + reference-graph analysis (explorer, 2026-07-14),
recorded inline. **Companion:** `.claude/plans/sn_doc_architecture_231.md` (the in-flight
per-page campaign, which this re-frames).

---

## 0. The one-sentence thesis

> **The corpus loses a competition it should win.** When an agent needs the DD closure,
> Hébert §3.9 is *cheaper to load* than a 19,408-line page — so the agent greps fragments
> instead of reading the frame. The fix is not "write better"; it is to make the **retrieval
> unit** loadable and to **foreground exactly what literature cannot give**.

The reason to read ORPHEUS before Hébert is that Hébert cannot tell you *this project's
conventions*, *which test pins this equation*, *what breaks when the α-dome goes negative*,
or *why the frame is Petrov-Galerkin and not Galerkin*. Those are **the product**. The
standard derivations are connective tissue. Today the product is buried inside the tissue.

---

## 1. Diagnosis — five defects (only one is "too long")

| # | Defect | Evidence |
|---|---|---|
| **D1** | **Granularity.** The retrieval unit is mis-sized. | 10 MONOLITHs (>2000 ln) hold **62,944 ln = 71.6%** of the corpus in **14%** of pages. Only 31% of pages sit in the 300–2000 band. The corpus is **bimodal**: giant or tiny, thin middle. 19,408 ln ≈ ~250k tokens = unloadable ⇒ unread. |
| **D2** | **Namespace + navigation.** No hierarchy ⇒ no routing. | `docs/theory/` is **100% flat** (39 files, 0 subdirs). `theory/index.rst` sets **`:maxdepth: 1`** on all three toctrees ⇒ the sidebar shows page titles only, **never a section**; **8 of 10 monoliths** sit under it. 36 pages compensate with a page-local `.. contents::`. `theory/index.rst` has **exactly one section heading** (`Infrastructure`) but **three captions**, two of which resolve to *other pages* — it is a 1-section page impersonating 3 branches. |
| **D3** | **Adjacency.** Coupled topics are split; the split is fought by the graph. | Community detection (unseeded) finds **C2 = 14 pages / 46,986 ln = 53% of the corpus** — the SN∪algebra∪frame∪BC∪loss-reps core — **split across two of the three asserted branches**. The #1 mutual pair `discrete_ordinates ↔ operator_algebra` (**96 refs**, 4.4× the runner-up) straddles that cut. |
| **D4** | **Role-mixing + mis-filing (bidirectional).** | SN carries **2,164 ln** of operator-algebra content; `operator_algebra` carries **1,137 ln** of SN-only angular windowing (14 `sn-*` labels) and **1,679 ln** of BC (14 `bc-extraction-*` labels); `collision_probability` carries **2,126 ln** of Peierls (**49% of the page**). **130 section headings are named after a closed campaign.** |
| **D5** | **Purpose.** The differentiators are buried, and drift is unpoliced. | **The governing equation contradicts itself across 3 pages** — the SN machine header stated the retired `L: streaming + boundary` / `(L+C−S−F/k)` while the same page at L14056, `boundary_conditions.rst:112`, `operator_algebra.rst:5028` and `diffusion_1d.rst:67` all say **B is a sibling** and **A = L+C−S−B**. Fixed @ `275a753a`; `index_convention.rst` still propagates the stale form at **5 sites**. |

**D5 is the load-bearing finding.** The monolith did not merely fail to inform — it
*misinformed*, in the machine header, which is an agent's first 50 lines. One equation
stated in 6+ places across 2 pages has no single source of truth; that is an architecture
symptom, not a typo.

### Supporting findings (each a concrete fix)

- **Name collision, not overlap** — `homogeneous.rst` (0-D infinite-medium *reactor*, 1,648 ln)
  vs `frame.rst`'s spatial *homogenization* (fine→coarse XS collapse). **Zero** content
  duplication (`homogeni[sz]`: 5 hits vs 229). The hazard is pure retrieval: "homogenization"
  loads the wrong page. **Fix: disambiguate titles.**
- **Three un-linked V&V homes** — `theory/verification.rst` (660) · `verification/*` (1,940) ·
  `testing/*` (978). They point at each other with **raw path strings, not roles**
  (`theory/verification.rst:436`, `verification/matrix.rst:836,1220`) ⇒ no link, no graph
  edge, **no `-W` warning**. `testing/architecture.rst` (633 ln of load-bearing contract) has
  **0 inbound**. Soft contradiction: `theory/verification.rst:15` states "no cross-verification"
  *without* the L4 qualifier its sibling supplies.
- **`glossary.rst` has 0 inbound and 0 outbound** — the TERM layer seeded in Phase 0 is
  unreachable. Wiring `:term:` from pages is what activates it.
- **`.inc.rst` orphans built twice** — the two capability matrices render standalone *and* are
  `include::`d. One-line fix: add `'**/*.inc.rst'` to `exclude_patterns`.
- **No corpus reading order exists.** The *only* prerequisite declaration in 87,930 lines is
  the `depends_on:` in the SN machine header — inside a directive that isn't registered yet.
- **`ERR-NNN` is NOT archaeology** — it is durable V&V vocabulary wired to
  `@pytest.mark.catches`. (477 of 491 raw `R-[0-9]` hits are `ERR-0` false positives; true
  count 7.) Do not sweep it.

---

## 2. Design principles

- **P1 — The chapter is the retrieval unit.** Target **300–1500 lines** (~5–20k tokens): one
  concept-cluster, loadable whole without regret. **>2000 lines is a defect**, not a style.
- **P2 — Router + leaves.** Every part gets a *cheap* index (~150 ln) that is a genuine map.
  The agent loads the router, then **exactly one** leaf. This is the anti-context-bomb shape.
- **P3 — The tree mirrors the code's dependency layering** (`data → geometry → numerics →
  transport → methods/<m>`). One mental model for code and docs; Nexus edges align; knowing
  where the code lives tells you where the docs live.
- **P4 — Adjacency is measured, not asserted.** Part boundaries must not cut a high-coupling
  cluster (§3b of the inventory). Where they must, the coupling is usually an artifact of
  **duplication** — fix the duplication, not the boundary.
- **P5 — Tracks, not one linear order.** A book has one sequence; we can publish several
  (newcomer · modifying-the-sweep · porting-an-equation · debugging-a-wrong-k). This is a
  strict improvement over a book and it is what makes the corpus agent-navigable.
- **P6 — A chapter is self-contained on ONE concept**: its math **+** its ORPHEUS realization
  **+** its gotchas **+** what verifies it. Do **not** split "the math book" from "the
  implementation book" — that manufactures twin paths and violates P4.
- **P7 — Name the design, not the campaign.** Issue refs **only** where they point at
  **open/deferred** work (a live seam). *Never* as the label of a landed design.
  ✗ "the Issue #282 route (a) design marches ψ½…" → ✓ "ψ½ is marched directly from the true
  q½ source." ✓ keep: "the anisotropic angular floor — deferred to #229."
  **A heading named after a closed campaign is the purest form of noise** — it is the
  retrieval unit's own name. (130 such headings today.)
- **P8 — Lead with what literature cannot give.** Every part index opens on conventions,
  failure modes, and the code/test map. That is the moat and the reason to read this first.
- **P9 — General-vs-realization split.** A shared concept lives **once**, in its foundation
  page; each method chapter carries only *its realization* + a link. **Worked precedent:** the
  `frame.rst` ↔ SN-stub split (`3de597a3`) — it works, and it is the template for un-mixing
  D4.
- **P10 — The `:label:` prefix is a contract, not a fossil.** `sn-*` must live on an SN page.
  Re-namespace as content moves, or the oracle becomes a lie.

---

## 3. What books do — and what they miss

*Grounded by a verified literature survey (2026-07-14). Local holdings: **Hébert Ch.3** and
**Stacey Ch.9** (both read); **B&G 1970** (contents verified via Internet Archive).
**Lewis & Miller could not be verified at all** — see §3.4.*

### 3.1 The recurring skeleton (verified across all three texts)

```
(1) continuous 1-D equation, conservation form
(2) angular collocation: ordinates + weights
(3) [1-D] the S_N ≡ P_{N−1} equivalence
(4) slab: spatial differencing (diamond) + the sweep
(5) curvilinear: redistribution → α-recursion → starting direction → sweep
(6) multidimensional: level-symmetric quadrature + 2-D sweep
(7) acceleration (or a pointer elsewhere)
```
**Always angular before spatial. Slab strictly before curvilinear. Acceleration last.**

### 3.2 The page budget — the niche, quantified

| Text | S_N pp | Cylinder? | Adjoint **of S_N**? | Eigenvalue posed independently? | **Verification?** |
|---|---|---|---|---|---|
| **Bell & Glasstone 1970** | **38** | §5.3e "general" only | **YES (§6.2b)** | **YES (§1.5, Ch.1)** | No |
| **Hébert 2009** | 19 | **YES, 5 pp** (only text) | No | No — filed under **CP §3.8.2** | **0 hits** |
| **Stacey 2007** | 17 | **No** | No (Ch.13) | No — filed under **spherical S_N §9.9** | **0 hits** |
| Duderstadt & Hamilton 1976 | **~0** | — | — | — | — |
| Larsen & Morel 2010 (review) | 84 | *unverified — closed access* | | | |
| Adams & Larsen 2002 (review) | 157 (accel.) | *unverified — closed access* | | | |

> **The niche, in the survey's own words:** *"The 1970 text is still the most complete textbook
> treatment of S_N-the-method, and the two modern textbooks are each ~18 pages that omit
> verification and convergence theory entirely. The real modern content moved into two review
> monographs (241 pages combined) that are both closed-access. **That gap — between an 18-page
> textbook section and a 157-page paywalled review — is the niche your mini-book occupies.**"*

### 3.3 What they MISS — verified

| # | Gap | Verified evidence | ORPHEUS's answer |
|---|---|---|---|
| **M1** | **No conventions crosswalk — and the texts contradict each other AND themselves.** | **Hébert Ch.3 has no notation table**: `notation\|nomenclature\|list of symbols` → 7 hits, **all 7 are the phrase "to simplify the notation."** **Within one chapter**: §3.9.1 1-D Gauss-Legendre ⟹ **Σw = 2**; five pages later (3.363)–(3.364) ⟹ **Σw = 1 over the positive octant** — same symbol `w_n`, two normalizations, no note. The **(2ℓ+1) prefactor** carries `4π` in (3.30) but `2` in (3.425) — same object, tied silently to dimensionality. **The scattering arrow points OPPOSITE ways**: Hébert `Σ_s(E ← E′, Ω ← Ω′)` *destination-first*; B&G `σ(x, E′ → E)` and Stacey `Σ_s(Ω′ → Ω)` *source-first*. | **Part 0 `conventions/` + the per-method machine header + a literature crosswalk at every import site.** ⟹ *"code that runs, converges, and is wrong by a constant — the worst failure class there is."* **ERR-025 was a missing 1/W; ERR-039 was the SH prefactor.** |
| **M2** | **No symptom→cause diagnostics** (two partial exceptions). | Best in corpus: **Hébert (3.388)–(3.389)** — a positivity *predicate* + its three named drivers. Stacey §9.10 is one qualitative paragraph, then defers. Stacey states a real failure mode with **no diagnosis and no fix**: *"The synthetic method may even become unstable with small mesh spacing."* | **Gotchas (consequence→manifests→catcher) + ERR catalog + the `symptom → chapter` diagnostic table in each part index.** |
| **M3** | **The four-way blur is VISIBLE IN THE TOC.** | **Hébert files power iteration as §3.8.2 — a subsection of the *collision-probability* method**; §3.9 (S_N) has **no eigenvalue subsection at all**; power iteration is re-derived a **third** time at §3.11.4 for Monte Carlo. A geometry-free, method-free concept instantiated 3× in 3 silos, factored out **zero** times. **Stacey files "Calculation of Criticality" + "Acceleration of Convergence" under §9.9 — *1-D **spherical** S_N***. | **The operator algebra separates posing / discretization / algorithm / resolvent.** ⟹ **the spine, not an appendix.** |
| **M4** | **Verification is a ZERO, not a thin spot.** | `verification\|benchmark\|manufactured solution` over **all 122 pp of Hébert Ch.3 → 0 hits**; over **all 80 pp of Stacey Ch.9 → 0 hits**. B&G: 2 book-wide. | **The V&V ladder + per-equation slices.** *"S_N is a method where wrong code converges."* **The single largest value-add over the entire canon.** |
| **M4b** | **No convergence theory; no asymptotic diffusion limit.** | `spectral radius\|Fourier analysis\|convergence rate` → **0 hits in both**. `diffusion limit` → **0 hits in both**. Both say SI is "slow" when c→1; neither gives ρ ≈ c or any predictor. The Larsen–Morel–Miller admissibility theory is in **neither textbook**. | **The reader learns the closure but not the criterion that SELECTS it.** Our local `Larsen-Morel-Miller(1987)` (spatial) + `Bailey-Morel-Chang(2010)` (angular) factor it per-axis — a first-class **admissibility** section no textbook has. |
| **M5** | **⚠ CORRECTED — curvilinear is NOT hand-waved. It is *inconsistent across the canon*.** | **The same symbol α satisfies FOUR recursions in three texts, none acknowledging another exists:** Stacey (9.213) `α−μw` · **Hébert (3.424) sphere `α−2wμ` (×2 vs Stacey)** · **Hébert (3.399) cylinder `α+wμ` (sign flip vs its OWN sphere, 4 pp apart)** · B&G (5.21) `α−μw(A_{i+1}−A_i)` (**ΔA folded inside α**). All correct in their own conventions. Sub-gaps: **the cylinder is nearly unserved** (Stacey has none; B&G gestures; Hébert's 5 pp is the only real one); **the pole is asserted, not analyzed** (B&G p.232 confesses *"the central flux is not exactly isotropic"*); **ΔA/w is never named as an object**. | **The α-normalization crosswalk table is a literal deliverable** — the highest-density verified cross-source knowledge in the survey, existing nowhere in print. Plus: name ΔA/w; analyse the pole structurally. |
| **M6** | **Adjoint bolted on late — except once, in 1970.** | Hébert Ch.3: 14 "adjoint" lines, **every one outside §3.9**; the adjoint lives in Ch.5 under *diffusion* GPT. Stacey: **one** hit in 80 pp, a forward reference to Ch.13. **B&G Ch.6 §6.2b "One-Speed, P_L, Diffusion, and S_N Theories" constructs the adjoint OF THE S_N DISCRETIZATION** — the only text asking whether the adjoint commutes with discretization. | **Adjoint as a structural thread.** B&G §6.2b is the citable precedent for our `A.H.inverse() ≡ A.inverse().H` swap law; `frame.rst` §2b is the modern instance (the adjoint *forces* Petrov-Galerkin). |
| **M7** | **Methods are siloed; comparative placement ≈ absent.** | Exhaustive grep of Hébert Ch.3 for comparative language: the **only** hits are §3.11.5 comparing *MC estimators to each other*. Five consecutive silos; **no "when do you reach for S_N instead of CP?" paragraph exists.** Only *equivalences* survive (S_N ≡ P_{N−1}). | **`methods/index` — the comparative map.** Our local `Sanchez(1982) A review of neutron transport approximations` is the missing chapter: organized by *approximation*, not by *method*. |
| **M8** | **The selecting invariant is often omitted** (B&G excepted — see §3.5). | — | **Derivations actually done**, opening on the invariant that selects the scheme. |
| **M9** | **No link to running code.** (Books can't.) | — | **Implementation map + Nexus.** |

**Synthesis:** our differentiators are *precisely* what the canon cannot or does not do.
Foreground them; the textbook content is connective tissue. **That is why an agent reads this
first.**

### 3.4 Gaps in our own evidence (honest scope)

- **Lewis & Miller could not be verified at all** (Wiley/ScienceDirect/HathiTrust/Archive all
  403 or absent). It is the text Stacey *explicitly defers to* for acceleration and the likely
  closest structural precedent. **ACTION: add it to `scratch/literature/`.**
- **Larsen & Morel 2010** (84 pp, the post-1968 S_N review) and **Adams & Larsen 2002** (157 pp)
  are closed-access with no OA mirror; their section structure is the most relevant map in
  existence for exactly the parts textbooks omit. Their own abstract states our thesis:
  *"several books and reviews … have been published, but none of these covers the advanced work
  done during the past 20 years."*
- `Stammler(1983)` and `Ligou(1982)` are present locally but **image-only, no OCR layer**.

### 3.5 Pedagogy to steal (verified)

1. **B&G §1.5 — pose the k/α eigenvalue in Chapter 1, before any discretization exists**
   (incl. §1.5f "Comparison of k and α Eigenvalues"). **Alone in the corpus.** Every method
   chapter then *inherits* the posing instead of re-deriving it. ⟹ our **algebra-as-spine**.
2. **B&G §5.3b–c — the three-beat:** state the invariant that *selects* the scheme → derive →
   read every term back physically. B&G's own justification: *"in the absence of such a
   principle, the possible difference equations are so numerous that a good choice is difficult
   to make except by a process of trial and error."*
3. **B&G proves the closure rather than choosing it, then audits it, then confesses.** Where
   Stacey writes *"By choosing α_{1/2} = 0"*, B&G **derives** it as forced by coarse-cell
   balance consistency (5.22), runs a positivity audit, then admits the residual error. ⟹ our
   "test the intrinsic properties" discipline, as prose.
4. **Hébert §3.9.6 — acceleration IS preconditioning of ONE recursion.** `M = D + W⁻¹E`,
   `A = D − W`, and *"**Setting M = D reduces the preconditioned form to the non-accelerated
   scheme.**"* **This is the operator-algebra framing, in a 2009 textbook** — the only instance
   in the canon. **We can cite, not invent.** (Hébert reaches for it only under acceleration; we
   hoist it to the spine.)
5. **Hébert states the design constraint before choosing the machinery:** *"the S_N method must
   be consistent with diffusion theory, even at low N values"* — so the quadrature reads as a
   *consequence*, not a convention.
6. **Hébert's equivalence-as-oracle** (S_N ≡ P_{N−1}): *"any discrepancy … is only due to the
   Marshak versus Mark treatment of vacuum boundary conditions. Infinite lattice solutions are
   identical."* — **a bit-identity claim with ONE named degree of freedom, sitting unexercised
   in prose.** Promote every such equivalence to a *testable predicate + its pinning test*: the
   cheapest V&V content in the survey.
7. **Hébert's executable pedagogy** — a runnable kernel per method; the matvec specified as a
   homework instruction (*"implementing a function `atv` returning an improved flux array…"*).
8. **Stacey's per-context nomenclature figures** — esp. **Fig. 9.20 "Change in angular
   coordinate μ as the neutron moves"** (makes redistribution obvious) and **Fig. 9.21 "Sweep of
   the space–angle mesh"** (*the sweep DAG, already drawn in 2007*). Complementary to — not a
   substitute for — the global crosswalk.
9. **B&G's title makes a distinction the field LOST:** *"Discrete Ordinates **and** Discrete S_N
   Methods"* — the angular collocation vs Carlson's specific difference scheme. **Naming them
   apart is half the fix for M3.**
10. **Hébert factors geometry out of method** — §3.6 "The first-order streaming operator"
    (Cartesian/cylindrical/spherical) sits **before** P_n, CP, S_N, MoC, so all four share one
    derivation. *"The only successful de-duplication in the corpus."* **This is the user's
    insight one level down: Hébert factors streaming by GEOMETRY; the path-integral root (§3.6)
    factors it by METHOD. Same move, one level up — with a canonical precedent.**

---

## 3.6 THE ROOT — the page no book has

*User insight (2026-07-14), adversarially verified (cross-domain-attacker) and **substantially
corrected**. What survives is stronger than what was proposed. **Read the MUST-NOT list before
writing a word of this page.***

### The insight

> The transport methods are not five different subjects. They are **five discretizations of one
> object** — the sum over neutron histories. The reason the operator algebra is powerful is that
> **C, S, F are the same operators in every method**; what varies is how the *streaming/propagator*
> is represented.

**This is not an analogy — it is a shared-code fact.** `MultiplicationOperator`,
`IsotropicScattering`, `IsotropicN2N`, `FissionOperator` are *literally the same Python objects*
imported by SN, diffusion, and homogeneous from `orpheus/transport/operators/`
(`diffusion/solver.py:230-243`; `homogeneous/solver.py:143-146,192-193`). **Keep this; it is the
frame's best asset.**

**Canonical precedent:** Hébert §3.6 factors the streaming operator out (Cartesian/cylindrical/
spherical) *ahead of* P_n, CP, S_N and MoC — *"the only successful de-duplication in the corpus."*
He factors streaming by **geometry**; the root factors it by **method**. Same move, one level up.

### What the root page MAY claim (verified)

1. Every ORPHEUS method solves for the **first moment of one branching PDMP** (piecewise-
   deterministic Markov process). The **many-to-one lemma / spine decomposition** (Lyons–Pemantle–
   Peres 1995; Hardy–Harris 2009) is the theorem that collapses the branching walk to a **single
   weighted path** — *linearity is what makes that collapse legal.* **Fission does not break the
   path reading; it makes the multiplicative functional exceed 1.**
2. The honest name chain — **Hille–Yosida / Dynkin resolvent** (root) → **Dyson–Phillips**
   (time-dependent) → **Neumann–Peierls collision-order series** (stationary) → **Feynman–Kac for
   a PDMP** (the MC reading). ORPHEUS **already has the right name**: `operator_algebra.rst:928`
   says *"the Peierls collision-number expansion."* Keep it.
3. `L` is **unbounded and rank-deficient**; `S`, `F` are **bounded**. *That* is why the splitting
   is around `A₀ = L+C` and why `(L+C)` is bundled as one invertible unit — **not** the scalar
   analogy `(3+5)⁻¹ ≠ 3⁻¹+5⁻¹` at `operator_algebra.rst:734`, nor the resistor identity at `:906`.
   **The page never says the word "unbounded"; that is the actual reason.**
4. **A method is a choice of GENERATOR SPLITTING** `𝒜 = 𝒜₀ + P` — *this is the real root branch,
   and it is better than the generator/propagator dichotomy:*
   | Splitting | Σ_t sits in | S sits in | Series | Methods |
   |---|---|---|---|---|
   | **killing** | the functional `e^{−∫Σ_t}` | the **source** | Dyson = collision-order = **SI** | SN, MoC, CP, Peierls |
   | **jump** | the jump **rate** | the jump **kernel** | none — the process *is* the answer | MC analog |
   | **majorized jump** | the majorant `Σ_maj` | kernel + virtual self-scatter | none | **ORPHEUS MC (Woodcock)** |
   **The Radon–Nikodym derivative between the last two IS the delta-tracking rejection
   probability** — a real bridge, currently un-named in the codebase.
5. **C, S, F are method-invariant reaction operators AT FIXED MULTIGROUP DATA** (scope condition —
   see MUST-NOT #3), demonstrated by shared code.
6. **B factors as method-invariant LAW ∘ method-specific TRACE** (`transport/method.py:199-237`) —
   and **the trace is set by the ANGULAR representation, not by streaming.**
7. **The corrected taxonomy — THREE ORTHOGONAL AXES:**
   | Axis | Values | SN | MoC | CP | MC | Diffusion/P_N | Case/F_N |
   |---|---|---|---|---|---|---|---|
   | **A1** — how `(L+C)⁻¹` is realized | DAG / track / region-pair / sampled / **not at all** | DAG (Padé exp) | track (exact exp) | Peierls kernel | sampled | **not realized (a limit)** | **not realized** |
   | **A2** — where `S` is resummed | outer Neumann (SI) / direct inverse / **exact spectral** | SI or Krylov | SI | SI (1 Jacobi step) | in the process | direct LU | **`Λ(z)` — closed form** |
   | **A3** — angular representation | ordinates / harmonics / continuous / **Case-ν** | ordinates | ordinates+tracks | integrated out | continuous | harmonics | ν-spectrum |
8. **SN and MoC are the SAME side.** `diamond.py:588` with slab neutrals and `τ = Σ_t V/|μ|`:
   `a = (2−τ)/(2+τ) = (1−τ/2)/(1+τ/2)` = **the [1/1] Padé approximant of `e^{−τ}`**; MoC's
   transmission is `e^{−τ}` exactly. **Closures are rational approximants of one propagator** — so
   "negative flux at τ>2" is *the pole of the [1/1] Padé*, not a DD pathology.
9. **Diffusion/P_N are the ONE genuine exception** — a P1/asymptotic **limit**, not a quadrature of
   the resolvent. *That* is why diffusion's solve is elliptic-self-adjoint while every other
   method is characteristic-triangular. **One principle, one exception.**
10. For a **fixed ordinate on a structured Cartesian mesh with a Cartesian closure**, the
    cell-dependency digraph is acyclic **by a lattice product-order argument** (`sign(Ω_x)·i +
    sign(Ω_y)·j` is a strict potential) — a **MESH theorem, not a characteristic theorem**.
    Certified by `triu == 0` (#284) + LAPACK ≡ sweep at ~6e-16.
11. **Acyclicity is a theorem about a (mesh, closure, boundary) TRIPLE, certified per case** —
    and #282 is the record of a defensible closure that **broke** it (cold residual 5.18e5, NaN).
    This is *stronger* than "the sweep is the graph structure," because it is falsifiable and
    ORPHEUS has already falsified an instance.
12. The extraction criterion is an **SCC decomposition** on the `(face, ordinate)` digraph — not a
    boolean. ORPHEUS's own grand-report already names `SweepStrongComponent` /
    `ReflectiveSweepCycle`.

### What the root page MUST NOT claim (each refuted with evidence)

1. **NOT "Feynman–Kac" as the root label.** For `(L+C)⁻¹` the process is **deterministic**
   straight-line motion; the expectation is over a one-point measure — technically true, vacuous.
2. **NOT that the Neumann reading and the Markov-jump reading are "equivalent."** They are
   **different splittings** (MAY #4), and the deterministic/stochastic split falls on that line.
3. **NOT that the difference is "confined to streaming."** **CE-1:** diffusion's `D = 1/(3Σ_tr)`,
   `Σ_tr = Σ_t − Σ_s1_out` (`mixture.py:166,191`) **relocates the ℓ=1 scattering moment into L** —
   and `Σ_t` then sits in **both** `L` and `C`, whereas SN's `L` is the σ-free leaf
   (`operator_algebra.rst:956-957`). **The L/C/S PARTITION is method-dependent.**
   **CE-scope:** multigroup condensation is a solution-weighted **Petrov-Galerkin** projection
   owned by *no* operator — so a 2-group SN's `S` and a 2-group diffusion's `S` are **different
   operators** if each was condensed with its own flux. ⟹ **claim holds at fixed multigroup data.**
4. **NOT that `A = L+C−S−F−B` is universal.** **CE-3b:** CP folds B into `P_inf` as a rank-1 kernel
   update (`cp/solver.py:395`); **Case/F_N have no C/S/F at all** — the ν-spectrum diagonalizes the
   **full** `(L+C−S)`, so `c` is a parameter of `Λ(z)`, not an operator (CE-2).
5. **NOT "B is tied to streaming"** — B's trace is set by the **angular** representation (CE-3).
6. **NOT the generator/propagator dichotomy.** The sweep **IS** `(L+C)⁻¹` exactly; DD is MoC's
   exponential at [1/1]. **SN is "propagator-side" by its own proposer's definition.**
7. **NOT that Case/F_N are "spectral-in-angle."** They are **spectral-in-the-full-generator** — an
   **A2** value (exact resummation of the collision series), not A3.
8. **NOT that the k-eigenvalue path integral converges.** `ρ(A⁻¹F) = k`; for `k>1` the generation
   series **diverges** — the sum-over-histories reading is FALSE until rescaled by `1/k`, which is
   not known a priori. **The k-posing is exactly where the path reading yields to a spectral
   statement** (Krein–Rutman *is* Perron–Frobenius on the mean-offspring operator; α is the
   **Malthusian parameter**).
9. **NOT `ρ = c` as an identity.** It is a **supremum** (`ρ ≤ c`), attained only in the
   infinite-homogeneous limit; finite media with leakage give `ρ < c` strictly. And **the
   hypothesis can be violated by ORPHEUS's own `S = Σ_s0ᵀ + 2Σ₂ᵀ`**: the bound exceeds 1 **iff
   `Σ_2n_out > Σ_c + Σ_L + Σ_f`** — checkable, not vacuous. **State the hypothesis.**
   ((n,2n) breaks sub-stochasticity **before fission does**.)
10. **NOT "reflective B ⇒ cycle ⇒ extraction forced."** ORPHEUS keeps a **specular reflection
    inside the walk**, certified lower-triangular (the r=0 pole mirror,
    `loss_representation/__init__.py:2840-2846,3654-3660`) — a **forward edge**, because the order
    visits μ<0 first. A **single** reflecting face is acyclic; you need a **closed loop** (both
    faces). ⟹ SCC, not boolean.
11. **NOT a fork of the existing taxonomy.** `reference_solvers.rst:161-194` already ships a
    three-meanings split called *"load-bearing for the entire reference-solver section."*
    **Reconcile or retire it — do not fork.** Likewise `peierls.rst:37-39` explicitly denies the
    premise for its own families: *"The two architectures are **not** different discretisations of
    the same operator."*
12. **NOT cite `eigenvalue.py:31,37-43,111`.** It was stale and it manufactured a wrong taxonomy in
    this very session. **Fixed @ (this branch).**

### Blocking defects to fix before the root page cites them

| # | Defect | Evidence |
|---|---|---|
| **X1** | **`eigenvalue.py` docstring claimed CP/diffusion/homogeneous have "no (L,S,F) factorization"** + "CP BiCGSTAB" + "Diffusion BiCGSTAB". All false. **The direct cause of a wrong taxonomy this session.** | **FIXED this branch** |
| **X2** | `operator_algebra.rst:4851` writes `L_full`'s `(s,s)` block as **0**; `:4880-4882` says the inflow row carries the **identity** — and that identity is exactly what makes the trace rows diagonal and the whole thing triangular. | `:4851` vs `:4880-4882` |
| **X3** | `operator_algebra.rst:742` dismisses Sherman–Morrison–Woodbury blanket ("applies only under low-rank structure, which the dense collision diagonal is not") — **true of C, false of B**. A correct dismissal generalized into a blanket one; it **blocks a real algorithmic win** (below). | `:742` vs `cp/solver.py:395` |

### The algorithmic payoff (a real finding, route to an issue)

`cp/solver.py:389-397` computes `P_inf = P_cell + outer(P_out,P_in)/(1−P_inout)` — **Sherman–
Morrison in longhand, unnamed**: the closed-form resolvent of the boundary cycle. **CP *solves* the
SCC that SN *iterates*.** `B` is low-rank per face (white/albedo: **rank 1**; specular slab: an
ordinate permutation, rank N/2). So `(L+C−B)⁻¹` by Woodbury = **two sweeps + one scalar division,
exact, zero boundary iterations** — the classical **response-matrix method**, falling out of the
SCC frame. Blocked by X3.

### What the root page is, structurally

`foundations/path_integral.rst` (~1200 ln) is **the corpus's first page** and the parent of
`methods/index`. It answers: *what is the one object; what is invariant; what varies; on which
axes; and where does each method land.* Every method book's ch.1 then **derives** its method from
the root instead of asserting it. **The eigenvalue posing lives here too** — B&G §1.5 put k/α in
Chapter 1 *before any discretization*, alone in the canon (§3.5.1); it is a property of the
**operator**, so every method inherits it.

---

## 4. The proposed corpus tree

Mirrors the code layering (P3). `→` marks a move; sizes are current line counts.

```
docs/theory/
  index.rst              THE MAP — layers, reading tracks, "start here". Fix the
                         1-section/3-caption defect; raise :maxdepth: to 2.
  glossary.rst           (wire :term: from pages — it has 0 inbound today)

  conventions/           ── PART 0 · the crosswalk (M1 — the biggest differentiator)
    index.rst
    notation.rst           symbol table + the ORPHEUS↔Hébert↔Lewis&Miller crosswalk   [NEW]
    normalization.rst      weight sums, 4π, the (2ℓ+1) prefactor, the 1/W trap         [NEW]
    indexing_and_layout.rst          ← index_convention.rst (layout half, 1,681)
    cross_section_conventions.rst    ← index_convention.rst (XS half)

  foundations/           ── PART 1 · the shared math (numerics + general transport)
    index.rst
    transport_equation.rst           ← SN L259 (99) + homogeneous's mg-eigenvalue    [EXTRACT]
    cross_section_data.rst           (1,381)
    geometry_and_meshes.rst          ← structured_geometry.rst (462)
    discrete_measures.rst            (821)
    spherical_harmonics.rst          (594)
    frame/                           ← frame.rst (3,706) → index + general + PG +
                                       galerkin + homogenization + condensation (~5 ch)
    operator_algebra/                ← operator_algebra.rst (11,435) − 1,137 (SN windowing
                                       → SN) − 1,679 (BC → BC) + 1,654 (from SN) (~9 ch)
    boundary_conditions/             ← boundary_conditions.rst (3,510) + 1,679 re-homed (~4 ch)
    eigenvalue_and_iteration.rst     THE single home for power iteration / posing
                                       (today: 3 homes — SN, operator_algebra, api/numerics)
    infinite_medium.rst              ← homogeneous.rst (1,648), RETITLED (kills the collision)

  methods/
    index.rst              THE COMPARATIVE MAP — why S_N vs CP vs MoC vs P_N (M7)   [NEW]
    sn/                    ← discrete_ordinates (19,408) + loss_representations (2,877)
                             + windowing (1,137) − algebra (1,654) − history (2,794)
                             − verification (3,904)  ≈ 15,000 → ~14 chapters  (§5)
    collision_probability/ ← collision_probability.rst (4,349) − 2,126 Peierls ≈ 2,200 (~3 ch)
    method_of_characteristics.rst   (1,496)
    monte_carlo.rst                 (1,324)
    diffusion.rst                   ← diffusion_1d.rst (1,450)
    spectral/                       ← pn, spn, bn, galerkin_spectral, spectral_collocation,
                                      spectral_resolvent, galerkin_sn_hybrid (7 stubs)

  references/            ── PART · Branch-1 truth set
    index.rst                        ← reference_solvers.rst (349)
    peierls/                         ← peierls (734) + peierls_nystrom (9,241)
                                       + 2,126 from CP  (~9 ch)
    trajectory_resolvent/            ← trajectory_resolvent (5,898) — re-cut BY GEOMETRY
                                       (sphere/cyl/slab/hollow/annulus), not by Phase (~5 ch)
    case_spectrum/                   ← fn_method (2,454) + singular_eigenfunction (2,076)
                                       (mutual cluster #5)  (~4 ch)
    sood_registry.rst                (919)
    escape_probability.rst           (78)

  verification/          ── PART · consolidate the THREE un-linked homes
    index.rst · ladder.rst (← testing/architecture) · matrix.rst (auto)
    reference_solutions.rst · cross_method.rst
    ← absorbs theory/verification.rst (660) + SN's Verification (2,548) + the
      2-D LD stress MMS (1,356)

  (multiphysics/ — fuel_behaviour, thermal_hydraulics, reactor_kinetics: PARKED,
                   fate undecided per the minimal-production-architecture ruling)
```

**Sizing:** ~63,000 monolith lines ÷ ~1,000/chapter ⇒ **~60–65 new chapter pages**. That is
the real scale of this work, and it is why it must be mechanized (§7).

### Why this survives the coupling graph (P4)

The inventory's #1 objection — "any Part boundary separating SN from the algebra will be
fought by 96 cross-references" — **is answered by the deduplication, not by merging them**.
Those 96 edges are inflated by the 2,164 lines of algebra content *sitting on* the SN page.
Once the general algebra lives once in `foundations/operator_algebra/` and SN carries only
its realization + a link (**P9**, the proven `frame.rst`↔SN-stub pattern), the coupling
collapses to a healthy "SN realizes the algebra" edge. **Same instrument, both directions:**
the SN-only windowing (1,137 ln, 14 `sn-*` labels) moves *the other way*, into SN.

---

## 5. The S_N mini-book

**Sequence rationale — the algebra is the spine (M3).** Textbooks build the method then
reveal the structure. We invert it: state `A = L + C − S − B` up front, then **each chapter
builds one operator**, then pose and solve. For a learner this is a roadmap; for an agent
`algebra.rst` is the single highest-value page in the corpus — the whole structural frame in
~600 lines. It is also the fix for the four-way blur books leave.

| # | Chapter | Content (source spans) | ~ln |
|---|---|---|---|
| — | **`index`** ★ | router: machine header · synopsis · annotated chapter map · reading tracks · **symptom→chapter diagnostic table** · V&V summary | 150 |
| 1 | `placement` | why S_N; the trade-space vs CP/MoC/P_N/diffusion (M7) | NEW |
| 2 | **`algebra`** ★★ | **THE SPINE** — A = L+C−S−B; Aψ = (1/k)Fψ; the five operators; posing vs discretization vs algorithm vs resolvent; "what each chapter builds" ← L135 Architecture + realization stubs | 700 |
| 3 | `angular_quadrature` | ordinates, weights, exactness, the families ← L358 | 600 |
| 4 | `scattering_and_moments` | S; Legendre moments; the SH frame; S=R∘Λ∘M; **angular windowing** ← L13030 + 1,137 re-homed | 1,400 |
| 5 | `discrete_balance` | the control-volume derivation; **ΔA/w**; per-ordinate consistency ← L491 | 500 |
| 6 | `spatial_closures` | DD · WDD · Morel–Montry · step ← L871 (1,161) + L2114 (390) | 1,550 |
| 7 | `linear_discontinuous` | the UBLD tensor-product bilinear cell system ← L2685 − MMS | 1,480 |
| 8 | `loss_representation` | (L+C) lower-triangular; solve/apply/**assemble**; the schedules ← loss_representations.rst | 1,500 |
| 9 | `sweep_1d` | cumprod scan; affine outflow ← L2513, L2570 | 400 |
| 10 | `sweep_multid` | wavefront; octant DAG; unified dispatch ← L5525+L5656+L6310 | 1,600 |
| 11 | `curvilinear` | α-dome; **the pole as a structural singularity**; Morel–Montry; the ψ½ starting-direction direct solve ← L6272 + L10810 + L11865 | 1,800 |
| 12 | `boundary_conditions` | S_N's trace realization ← L12624 (→ links foundations/BC) | 400 |
| 13 | `iteration` | SI (ρ=c); Krylov; the DSA seam (→ links foundations/eigenvalue) | 600 |
| 14 | `solver` | `SNSolver` as coordinator ← L14963 | 500 |
| 15 | `adjoint` | the dual operator; S†; daggered posing; φ* | NEW |
| — | `history` | the changelog **+ the 2,794 evicted lines distilled to rows** ← L18292 | 1,000 |

**Gotchas are inline, not a chapter.** The trap lives *next to* the thing it traps (P6); the
part index carries the **index** of them. That index doubles as the **diagnostic router**
(M2) — "k_eff diverges under refinement → ch. 5/11" — which is the first thing a debugging
agent should load, and the single feature no textbook has.

**Exact chapter boundaries come from the Haiku classification catalog (§7), not from this
table.** The spans above are indicative.

---

## 6. Template revision (falls out of parts-vs-chapters)

#231's 9-section template was designed for *page = subject*. With *page = chapter* it splits:

**Part index (the router)** — machine header · synopsis · annotated chapter map · reading
tracks · **symptom→chapter diagnostic table** · V&V summary · history (collapsed).

**Chapter** — one-paragraph synopsis (its own retrieval anchor) · Formulation · Discretisation
& realization (incl. the code map) · **inline Gotchas** · Verification (what pins this) ·
References.

This is a *refinement* of #231, not a reversal: the settled section semantics survive; only
their **home** changes (subject-level → router; concept-level → chapter).

---

## 7. Migration path

Ordered for reversibility and early payoff. **A and B are independently valuable.**

- **Phase A — Skeleton (mechanical, low risk).** Create the directories + part indexes;
  `git mv` pages; fix every `:doc:`/`:ref:`; raise `:maxdepth:`; fix `theory/index`'s
  1-section defect; add `'**/*.inc.rst'` to `exclude_patterns`; retitle
  `homogeneous`→`infinite_medium`. **No content rewriting.** Payoff: navigation exists.
- **Phase B — De-duplicate by the label oracle (P9/P10).** Move the four mis-filed blocks to
  their labelled homes (SN-algebra→foundations; windowing→SN; BC→BC; Peierls→references);
  re-namespace labels as they land. Payoff: the coupling graph becomes honest; D4 closes.
- **Phase C — Split the monoliths into chapters**, driven by the **Haiku fan-out catalog**
  (chunk → classify KEEP / ARCHAEOLOGY→history / DISTILL→gotcha / RELOCATE → aggregate →
  archivist executes → main-agent reviews diff-vs-catalog). Start with SN (§5). **This is
  Phase 1d of the in-flight campaign, re-scoped: "split into chapters" ⊃ "decompose in place."**
- **Phase D — Write the routers** (indexes, tracks, diagnostic tables). Payoff: P2 realized.
- **Phase E — Conventions part** (M1) — mostly new writing; the biggest differentiator.
- **Phase F — Archaeology sweep (P7)** — 130 campaign-named headings; the 3.7:1 (b):(a) ratio.
- **Phase G — Backfill** — equation labels, `:term:` wiring, `:cite:` migration, V&V slices.

---

## 8. Open decisions (for the user)

1. **Chapter filename prefixes** — `01_algebra.rst` (order visible in the tree, but reorder ⇒
   rename ⇒ ref churn) vs **no prefix** (toctree carries order; rename-free; greppable).
   **Recommend: no prefix.**
2. **Depth of the method sub-books** — S_N gets ~14 chapters. Do CP/MoC/MC get the same
   treatment now, or stay single pages until they grow? **Recommend: SN now; others on demand.**
3. **`conventions/` as Part 0** vs folding it into `foundations/`. **Recommend: Part 0** — it
   is the M1 differentiator and the "read first" surface.
4. **`index_convention.rst`'s 5 stale `(L+C−S−F/k)` sites** — fix now (small, but the page is
   slated for absorption) vs fold into Phase A/B. **Recommend: fix now** (Cardinal Rule 1 —
   it is wrong today, and an agent reads it today).
5. **Scope order** — Phase A (navigation) before Phase C (SN split), or SN-first?
   **Recommend: A → B → C.** A is cheap and makes everything else discoverable.

---

## 9. What this does NOT change

The #231 settled design stands: the section semantics, bibtex-from-Zotero, the glossary, the
equation-label policy, the code-prose rebalancing rubric, stub-and-track for Nexus-blocked
automation. This proposal adds the **layer above the page** that #231 assumed and never
specified.
