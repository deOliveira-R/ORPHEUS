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
  closest structural precedent. ~~**ACTION: add it to `scratch/literature/`.**~~
  **RETIRED 2026-07-15 — no OCR'd PDF exists, and the ask was misdirected anyway: the
  acceleration content we wanted is a SECONDARY summary of `Larsen-Morel-Miller(1987)` +
  `Larsen-Morel(1989)` + `Bailey-Morel-Chang(2010)`, all three of which are ALREADY in
  `scratch/literature/`. See §8.3.**
- ~~**Larsen & Morel 2010** (84 pp, the post-1968 S_N review) and **Adams & Larsen 2002** (157 pp)
  are closed-access with no OA mirror~~ — **BOTH NOW LOCAL AND OCR'd (user, 2026-07-16). This
  gap is CLOSED; the §3.2 table's two "unverified — closed access" rows are now VERIFIABLE and
  should be filled by a survey pass.** Their section structure is the most relevant map in
  existence for exactly the parts textbooks omit, and their own abstract states our thesis:
  *"several books and reviews … have been published, but none of these covers the advanced work
  done during the past 20 years."*
  - `scratch/literature/Nuclear Computational Science - A Century in Review.pdf` — the **whole
    Azmy & Sartori book**, 476 pp, searchable. Ch. 1 = Larsen & Morel, pp. 1–84 (TOC-verified:
    ch. 2 opens at p. 85). **Ch. 2 is "Second-Order Neutron Transport Methods" by E. E. Lewis** —
    i.e. the book we chased for the Lewis & Miller structural precedent contains a Lewis chapter,
    on even-parity/second-order transport.
  - `scratch/literature/Adams-Larsen(2002) Fast iterative methods for discrete-ordinates particle
    transport calculations.pdf`.
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

- **Phase A — Skeleton (mechanical, low risk). ✅ DONE @ `08e58ee6` (2026-07-15).**
  Create the directories + part indexes; `git mv` pages; fix every `:doc:`/`:ref:`; raise
  `:maxdepth:`; fix `theory/index`'s 1-section defect; add `'**/*.inc.rst'` to
  `exclude_patterns`; retitle `homogeneous`→`infinite_medium`. **No content rewriting.**
  Payoff: navigation exists. **See §7.1 for what actually landed** — the tree deviates from
  §4 in three places, each for a reason worth keeping.
- **Phase B — De-duplicate by the label oracle (P9/P10).** **REVISED after mapping — see §8
  "RULED (Phase B)".** THREE clean relocations (Peierls→references; windowing→SN; BC→BC),
  sequence 4→2→3, each an independent `-W`-gated commit with zero renames on the moved
  content. The proposed fourth (SN-algebra→foundations) is NOT a mis-file — it is a
  rename-heavy dedup, **de-scoped to the Phase C SN split**. Payoff: the inflated SN↔algebra
  coupling becomes honest; D4 closes.
- **Phase C — Split the monoliths into chapters**, driven by the **Haiku fan-out catalog**
  (chunk → classify KEEP / ARCHAEOLOGY→history / DISTILL→gotcha / RELOCATE → aggregate →
  archivist executes → main-agent reviews diff-vs-catalog). Start with SN (§5). **This is
  Phase 1d of the in-flight campaign, re-scoped: "split into chapters" ⊃ "decompose in place."**
- **Phase D — Write the routers** (indexes, tracks, diagnostic tables). Payoff: P2 realized.
- **Phase E — Conventions part** (M1) — mostly new writing; the biggest differentiator.
- **Phase F — Archaeology sweep (P7)** — 130 campaign-named headings; the 3.7:1 (b):(a) ratio.
- **Phase G — Backfill** — equation labels, `:term:` wiring, `:cite:` migration, V&V slices.
- **Phase H — THE ROOT PAGE** (`foundations/path_integral.rst`, §3.6) — *added 2026-07-16: §3.6
  specifies this page in 24 numbered claims and §4 calls it "the corpus's first page", but the
  migration path never scheduled it. That was a hole, not a decision.* It is the corpus's
  single most distinctive deliverable — the frame from which every method branches. **Gated:**
  §3.6's own "Blocking defects" list says the root page must not cite X2/X3 until they are
  fixed ⟹ **blocked on #298 + #299**; X1 is fixed. **Read §3.6's 12 MUST-NOT claims before
  writing a word.** Wants Phase I's survey as an input (the taxonomy needs the post-1968 map).
- **Phase I — Literature survey pass** (*added 2026-07-16*) — the two §3.2 rows that read
  *"unverified — closed access"* are now local and OCR'd: **Larsen & Morel 2010** (ch. 1 of the
  Azmy & Sartori book, pp. 1–84) and **Adams & Larsen 2002** (157 pp). Fill the table; harvest
  **M4b** (convergence theory / the asymptotic diffusion limit / a ρ ≈ c predictor — a **zero**
  in both textbooks) and **M2** (symptom→cause diagnostics). Larsen & Morel's section structure
  is also the closest structural precedent that exists for the S_N mini-book (§5) — it is what
  the retired Lewis & Miller ask was actually reaching for. **Independent of B/C/D — dispatchable
  in parallel any time.** The DSA primaries (`Alcouffe(1977)`, `Larsen(1982)` I+II,
  `Morel(1982)`, `Hammer-Morel-Wang(2019)`) are for **#2 / stencil-P3**, not this campaign.

---

## 7.1 Phase A — what actually landed (`08e58ee6`, 2026-07-15)

**The tree as built** (33 moves; 39 pages all accounted for):

```
docs/theory/
  index.rst · glossary.rst · verification.rst          <- root
  thermal_hydraulics.rst · fuel_behaviour.rst · reactor_kinetics.rst
  conventions/   index.rst [NEW] · index_convention.rst
  foundations/   index.rst [NEW] · operator_algebra · frame · boundary_conditions
                 cross_section_data · discrete_measures · spherical_harmonics
                 structured_geometry · infinite_medium
  methods/       index.rst (<- transport_methods) · collision_probability
                 method_of_characteristics · monte_carlo · diffusion_1d
                 sn/  index.rst (<- discrete_ordinates) · loss_representation
  references/    index.rst (<- reference_solvers) · peierls · peierls_nystrom
                 trajectory_resolvent · fn_method · singular_eigenfunction
                 galerkin_spectral · sood_registry · escape_probability
                 + the 7 reserved stubs
```

**User rulings (2026-07-15).**

1. **Sub-book dirs: SN only.** SN is the only sub-book §5 actually specs, and it is Phase C's
   first target — so its dir is earned and Phase C becomes pure *addition* (zero `:doc:`
   churn). Every other monolith stays flat in its part; **its dir appears when its split
   does.**
2. **V&V consolidation: deferred to its own phase.** Not mechanical — `docs/verification/
   matrix.rst` is auto-generated (`generate_matrix.py`), `testing/`+`verification/` are
   top-level *siblings* of `theory/` in `docs/index.rst` (folding them under `theory/`
   demotes V&V — an IA judgment), and the three homes soft-contradict each other.
3. **Multiphysics: left flat at the root, no effort** — they may be extracted from the repo.

**Deviations from §4, each with its reason — do NOT "fix" these back:**

| §4 proposed | What landed | Why |
|---|---|---|
| `diffusion.rst ← diffusion_1d.rst` | **kept `diffusion_1d.rst`** | `DiffusionMesh` **refuses a multi-D mesh at construction** (`augmented_mesh.py:130`). The `_1d` is load-bearing; the rename would have made the corpus misinform — the exact D5 disease. |
| `methods/spectral/` ← 7 stubs | **stubs stay in `references/`** | They live under `orpheus.derivations.continuous` and `reference_solvers.rst` calls them *"Reserved **reference** solvers"*. Moving them to `methods/` is a **re-classification** (is P_N a method or a reference?), not a move — a real open question, and not Phase A's business. |
| `geometry_and_meshes.rst ← structured_geometry.rst` | **kept `structured_geometry.rst`** | No justification beyond aesthetics; renames cost ref churn. `infinite_medium` earned its rename (a documented retrieval collision); this one didn't. |
| `theory/verification/` part | **`theory/verification.rst` at root** | Ruling 2. |

**Durable findings (these outlive Phase A):**

- **`-W` gates more than the rules claimed.** Measured, Sphinx 9.1.0: broken `:doc:` →
  `ref.doc` **warns**; broken `:ref:` → `ref.ref` **warns**. Only the *Python-domain* roles
  (`:func:`/`:class:`/`:mod:`) are silent without `-n`. `.claude/rules/coding-standards.md`
  said `:ref:` was silent — **corrected in this commit**. Consequence: **page moves and label
  retirements ARE gated by the build**; only raw text is not.
- **`:ref:` is path-immune** — all 1076 needed zero work. Labels survive any move. The
  restructure surface is `:doc:` + toctree + `.. include::` + raw text, nothing else.
- **A path built from SEGMENTS is invisible to a path-grep.** `REPO_ROOT / "docs" /
  "theory"` carries no `docs/theory/` substring. This was the single biggest hazard and no
  grep in the audit could see it — the **explorer found it by reading the tool**. When
  auditing a move, grep the **last path segment**, not the joined path.
- **All `:doc:` under `docs/theory/` are now ABSOLUTE** (Pattern 7). A future *source* move
  cannot break them — which is exactly what Phase C does. Keep them absolute.
- **The stale-path hazard is not hypothetical — it had already fired 3×**, incl. on this
  branch: the `galerkin_projection.rst`→`frame.rst` rename (`3de297a3`, two days earlier)
  orphaned 2 raw paths; `peierls_greens.rst`→`trajectory_resolvent.rst` orphaned 16; and
  `axis.py` promised a page (`sn_dim_agnostic.rst`) git shows was **never written**, citing a
  plan that no longer exists. **Every rename owes a raw-path sweep** — the build will not
  do it for you.
- **The fix surface excludes archaeology on purpose**: `.claude/{plans,agent-memory,scratch}`
  and `.claude/worktrees` (a separate checkout) were NOT rewritten. A plan's
  "← discrete_ordinates.rst" names a move's *source*; rewriting it corrupts the record.

**Phase A gate (what "done" meant):** `sphinx -W` exit 0 + clean log · 21 tests green ·
pyright 0 in every touched file (`orpheus/` unchanged at its 1 accepted #288 residual) ·
a **filesystem gate over raw path strings** (57 live / 0 dead), since no build checks those.

**⏭ NEXT = Phase B** (de-duplicate by the `:label:` oracle) — the tree now exists to move
content *into*. Then Phase C (split the monoliths, Haiku-catalog-driven, SN first).

---

## 7.2 Phase B — what actually landed (2026-07-16)

**Three clean relocations** (Block 1 de-scoped to Phase C per §8's Phase-B ruling), each an
independent `-W`-gated commit — byte-for-byte moves, zero renames on the moved content:

| Commit | Block | Move | Lines |
|---|---|---|---|
| `87c741c1` | B4 Peierls | `collision_probability.rst` slab/cyl/sphere → `references/peierls.rst` | 2,126 |
| `81f5dd78` | B2 windowing | `operator_algebra.rst` `sn-angular-windowing-*` → `methods/sn/index.rst` | 1,137 |
| `2a06779b` | B3 BC extraction | `operator_algebra.rst` Span A `bc-extraction-*` → `boundary_conditions.rst` | 1,107 |

(Plan-ruling commit `5d0f73fc`.) **Payoff realized:** D4 (role-mixing) closes for these three;
the inflated SN↔algebra and CP↔Peierls coupling edges are now honest.
`operator_algebra.rst` shed **~2,240 ln** (11,435 → 9,196 — windowing + BC Span A);
`collision_probability.rst` shed its 49 % Peierls half (4,349 → 2,223).

**B3's content-purity execution** (the one non-mechanical block): the general composite adjoint
STAYED on `operator_algebra`; its orphaned anchor was renamed `bc-extraction-g-adjoint` →
`g-adjoint` (build-gated — zero tree-wide after), and a two-way cross-ref bridges the split
narrative.

**Durable findings (outlive Phase B):**

- **The `-W` build gates the LINK, not the PROSE.** `:ref:`/`:eq:` are path-immune, so a move
  breaks no link (green build) while leaving two SILENT-prose falsehoods: **directional words**
  ("`:ref:`X` above" → now cross-page) and **page-qualifiers** ("`:ref:`X` in :doc:`<old-home>`").
  Neither is `-W`-visible. Catch them with a **three-way grep per move** — source-page refs,
  moved-block refs, and **bystander qualifiers**. Block 2 had 2 (found post-hoc); Block 3 had 4,
  **3 of them in `sn/index.rst`, a bystander neither source nor dest.** See [[lessons-L35]].
  **This is the load-bearing procedure for Phase C** (~63k lines moving).
- **Document-local RST substitutions do NOT cross a move.** Block 4's moved content used
  `|times|` (defined in the CP footer, staying) → the build ERRORED until the definition was
  copied to the dest page. Citations resolve cross-document via the global index (`[Kress2014]`
  did); substitutions have no global index absent `rst_prolog`. Distinguish a LaTeX `|m|` inside
  a `.. math::` (not a substitution — Block 2's non-issue) from a docutils `|token|`.
- **The label oracle drives the CLEAN case only.** Blocks 2/3/4 were clean because a *foreign*
  prefix on the host page is an unambiguous mis-file. Block 1 was not (honest `sn-*` labels →
  no oracle signal; a rename-heavy dedup) → de-scoped. The oracle finds mis-files; it does not
  adjudicate general-vs-realization.

**Deferred to Phase C/G (correct-but-unpolished, NOT false — not Phase-B debt):**
- 7 bare cross-page `:ref:`s in `operator_algebra` to the moved windowing labels (resolve fine;
  page-qualifying them is polish) + 4 now-redundant on-page self-`:doc:`s inside the moved
  windowing block (Phase C splits that content anyway).
- `|times|` now defined twice (CP + peierls) — a corpus-wide substitution duplication an
  `rst_prolog` consolidation would collapse (Phase-G single-source cleanup).
- `[Kress2014]` defined-in-CP-used-only-in-peierls (Phase-G locality tidy); the CP footer's
  "only citations unique to this page are defined locally" note is now mildly stale in spirit.

**⏭ NEXT = Phase C** (split the monoliths into chapters, Haiku-catalog-driven, SN first) — which
now ALSO absorbs the de-scoped **Block-1 SN-algebra dedup** (§8 Phase-B ruling). Phase I
(literature survey) remains independently dispatchable in parallel.

---

## 8. Decisions

### RULED (user, 2026-07-14)

- ✅ **Scope order: Phase A (navigation) BEFORE Phase C (the SN split).** A → B → C. A is cheap
  and makes everything else discoverable.
- ✅ **The stale `(L+C−S−F/k)` sites: FIXED NOW** — done @ `0ca0d378` (six sites across
  `index_convention.rst`, `diffusion_1d.rst`, `discrete_ordinates.rst`; B added as the fifth
  leaf with a resolve-verified xref). Two sites deliberately kept: a historical narrative
  (`index_convention.rst:395`, Phase-F material) and the *scattering carrier grid's* different
  "four leaves" (`operator_algebra.rst:2513,2541,2549`).
- ✅ **The review's findings: FILED.** **#298** (X2 — `L_full`'s `(s,s)` block written `0` but
  the table says identity; the identity is what makes the walk triangular) · **#299** (X3 — the
  blanket Sherman-Morrison dismissal: true of `C`, false of `B`) · **#300** (the payoff — close
  the boundary SCC in closed form via Woodbury on the rank-1 `B`; CP already does it, SN
  iterates it; blocked on #299). X1 (`eigenvalue.py`) fixed @ `018ecb7b`.

### RULED (user, 2026-07-15 — Phase A)

- ✅ **Depth of the method sub-books: SN now; others on demand.** Landed — only `methods/sn/`
  is a directory. See §7.1 ruling 1.
- ✅ **`conventions/` as Part 0: YES.** Landed as `docs/theory/conventions/` with a router that
  opens on the M1 argument (the canon contradicts itself on weight sums, the `(2ℓ+1)`
  prefactor, and the scattering arrow — each with an ERR catcher). It carries
  `index_convention.rst` today; `notation.rst` + `normalization.rst` are Phase E writing.
- ✅ **V&V consolidation: its own phase, not Phase A.** See §7.1 ruling 2.
- ✅ **Multiphysics: left flat at the theory root, no effort.** See §7.1 ruling 3.

### RULED (user, 2026-07-16 — Phase B)

Phase B mapped (explorer, verified against the current tree). Two forks ruled; the frozen
block-counts corrected. **Phase B = THREE clean relocations, not four.**

- ✅ **Block 1 (SN-algebra → foundations): DE-SCOPED to the Phase C SN split.** It is NOT a
  mis-file — its content carries honest `sn-*`/`si-*` labels, so the `:label:` oracle gives
  no signal. What is there is *duplication* (the `(A−S−F)ψ=q` spine restated 3× at
  `sn/index.rst` ~13330/13524/13953; `choosing-inverse-realisation` dupes
  `inverse-application-driver`/`green-operator`/`matrix-inverse-operator`;
  `cross-solver-eigenvalue-consumers` overlaps `matrix-inverse-consumers`) — a rename-heavy
  dedup, a *different operation* from the three clean relocations. Its content sits on the SN
  monolith that is **Phase C's first split target**, so its cleanup folds INTO that split
  (operate on the monolith once, not twice). Matches the plan's own "Phase C — de-duplicate"
  language.
- ✅ **Block 3 Span B (`bc-extraction-g-adjoint`): STAYS in `operator_algebra` (content-purity
  over the label oracle).** The section's content is the GENERAL composite adjoint
  `A†=G⁻¹AᵀG` (applies to L,C,S,F,B alike); general adjoint machinery belongs on the general
  page. ⟹ Block 3 becomes the **contiguous** cut `4780–5886` (Span A only, ≈1,107 ln); the
  now-orphaned `bc-extraction-g-adjoint` anchor is **renamed** to a general prefix (its
  narrative having left) with inbound refs fixed (`-W`-gated), and a cross-ref bridges the
  moved narrative ↔ the staying g-adjoint.
- **The three clean blocks (sequence 4 → 2 → 3), each an independent `-W`-gated commit, zero
  renames on the moved content:**
  - **B4 Peierls** — `collision_probability.rst:2205–4330` (contiguous, 2,126 ln = 49% ✓) →
    `references/peierls.rst` after `peierls-bc-foundations` §309. 3 anchors + 20 eq-labels +
    `issue-100-retraction` move unchanged. **Audit `peierls_nystrom.rst`'s ~20 `:doc:`→CP
    prose links** — repoint the Peierls-specific ones; a `:doc:` to a still-valid page is NOT
    `-W`-gated (silent drift).
  - **B2 windowing** — `operator_algebra.rst:7538–8674` (core + `windowing-retyped` tail,
    1,137 ln ✓) → `sn/index.rst` after `sn-scattering-adjoint` §13797. 14
    `sn-angular-windowing-*` + 5 eq-labels + `windowing-retyped` move unchanged.
  - **B3 BC extraction** — `operator_algebra.rst:4780–5886` (Span A, contiguous; the
    `affine-typed*` 5887–6408 AND the g-adjoint 6409–6985 STAY) → `boundary_conditions.rst`
    after `bc-realizer-layer` §307. 12 `bc-extraction-*` anchors + 7 eq-labels move unchanged.
- **Frozen-count corrections:** B3 is **12** `bc-extraction-*` anchors, not 14; B4's "3" was
  the *anchor* count (there are **20** eq-labels in the span); the "41 `sn-*`" figure did not
  reproduce. Authoritative per-family counts live in the move-map (this session).

### STILL OPEN

1. **Chapter filename prefixes** — `01_algebra.rst` (order visible in the tree, but reorder ⇒
   rename ⇒ ref churn) vs **no prefix** (toctree carries order; rename-free; greppable).
   **Recommend: no prefix.** *Not blocking until Phase C mints the first chapters.*
2. ~~**Are the 7 spectral stubs methods or references?**~~ **CLOSED (user, 2026-07-15):
   REFERENCES — the Phase A placement stands; §4's `methods/spectral/` was wrong.**
   **The user's criterion:** *"the folders inside `references/` are supposed to be reference
   analytical or semi-analytical solutions that can be used to verify mathematical correctness
   of problems."* **The evidence:** every one of the 17 pages in `references/` has its code home
   under `orpheus/derivations/continuous/` — `bn_method`, `escape_probability`, `fn_method`,
   `galerkin_sn_hybrid`, `galerkin_spectral`, `peierls_nystrom`, `pn_method`,
   `singular_eigenfunction`, `sood_registry`, `spectral_collocation`, `spectral_resolvent`,
   `spn_method`, `trajectory_resolvent` — and **no production package exists for any of them**
   (`orpheus/` ships `cp`, `sn`, `mc`, `moc`, `diffusion`, `homogeneous`; there is no
   `orpheus/pn/`). P_N here is the *continuous* semi-analytical reference (the Garcia stable
   spherical-harmonics benchmarks — six of those papers sit in `scratch/literature/`), **not** a
   discretized production solver. P3 (mirror the code layering) and the criterion agree.
   *Residual nit for Phase F, not a placement issue:* the `pn_method` / `spn_method` READMEs
   frame their subject as "the historical alternative to discrete ordinates" — method prose on a
   reference page.
3. ~~**Lewis & Miller**~~ **CLOSED (2026-07-16): the ask is RETIRED — do not re-raise it, and
   the literature gap it stood for is now FILLED.** L&M itself has no OCR'd PDF (user searched),
   but the ask was misdirected anyway: §3.4 wanted it for *acceleration*, and L&M is a
   **secondary** summary of primaries we already held (`Larsen-Morel-Miller(1987)` — the
   asymptotic/admissibility theory itself, already cited by M4b; `Larsen-Morel(1989)` II;
   `Bailey-Morel-Chang(2010)`; `Lathrop(2000)`). The corpus's own thesis is to beat textbook
   summaries by going to primaries — asking for the summary contradicted it.
   **What the user supplied instead (2026-07-15/16) is strictly better, and it is ALL LOCAL AND
   OCR'd** in `scratch/literature/` (now **73 files** — `.claude/rules/delegation.md`: check
   this folder FIRST):
   - **The complete DSA primary lineage** — `Alcouffe(1977)` (founds it) · `Larsen(1982)` Parts
     **I + II** (unconditional stability) · `Morel(1982)` (**highly anisotropic scattering** —
     directly relevant to ORPHEUS's `S = Σ_s0ᵀ + 2Σ₂ᵀ`) · `Adams-Larsen(2002)` (the 157-pp
     synthesis) · `Hammer-Morel-Wang(2019)` (nonlinear diffusion acceleration, **voids**).
     ⟹ **#2 / stencil-P3 (DSA) now has its sources. No literature ask remains for it.**
   - **`Nuclear Computational Science - A Century in Review.pdf`** — the whole Azmy & Sartori
     book (it was in the user's personal library all along). **Ch. 1 = Larsen & Morel 2010,
     pp. 1–84**, the post-1968 S_N review. This ALSO serves the structural-precedent purpose the
     L&M ask originally had — better, because it covers exactly the post-1968 material the
     textbooks omit.
   **Canonical citation (verified 3 ways — Crossref, the book's title page, its TOC):**
   Larsen, E. W., and J. E. Morel, "Advances in Discrete-Ordinates Methodology," Ch. 1 in
   *Nuclear Computational Science: A Century in Review*, Y. Azmy and E. Sartori (eds.), Springer
   Netherlands, Dordrecht, 2010, pp. 1–84. DOI `10.1007/978-90-481-3411-3_1`; ISBN
   978-90-481-3410-6 (print) / 978-90-481-3411-3 (online). *(Crossref stamps the chapter
   2009-12-24 online-first; cite 2010.)*
   **⚠ Nothing here is read yet. Acquisition ≠ survey.** Filling the §3.2 table's two rows and
   harvesting M4b/M2 is **Phase I** (§7), not a completed fact — do NOT cite these as surveyed,
   or restate their content, until someone has actually read them. ([[lessons-L33]]: a claim is
   not evidence; that applies to a claim about a PDF we merely possess.)

---

## 9. What this does NOT change

The #231 settled design stands: the section semantics, bibtex-from-Zotero, the glossary, the
equation-label policy, the code-prose rebalancing rubric, stub-and-track for Nexus-blocked
automation. This proposal adds the **layer above the page** that #231 assumed and never
specified.
