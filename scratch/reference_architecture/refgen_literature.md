# Verification-case anatomy: literature memo (#405 step 2)
Zotero not queried (no Zotero tools in this agent's list): no annotations checked.
Local library: none of Oberkampf-Roy, Knupp-Salari, Roache, Ganapol 2008, ANL-7416, ICSBEP present (79 PDFs checked by filename).

## 1. Oberkampf & Trucano (2007), "Design of and comparison with verification and validation benchmarks", NEA/CSNI workshop invited paper (NEA No. 6298), pp. 61-103. [M] read (pdftotext; printed p = PDF p + 60)
Local: scratch/literature/Oberkampf-Trucano(2007)Design of and comparison with verification and validation benchmarks.pdf (downloaded from oecd-nea.org, free). Journal twin: Oberkampf & Trucano, Nucl. Eng. Des. 238(3):716-743 (2008) [R: web search only, not DB-verified].
- p.64: "strong-sense benchmark" (SSB, from Oberkampf-Trucano-Hirsch 2004 AMR 57:345) = (a) purpose, (b) precise definition, (c) stated comparison requirements, (d) acceptance criteria defined; promulgated.
- p.70-71 (§2.1.2): AIAA Guide accuracy hierarchy: (1) analytical, (2) highly accurate ODE numerical, (3) highly accurate PDE numerical. "The modeling assumptions must be the same between the benchmark solution and the code being tested."
- p.79 (§3.1): four construction elements: a) purpose/scope, b) mathematical description, c) accuracy assessment + "pedigree of the evidence", d) documentation.
- p.80 (§3.1.1): purpose = TEXT only (searchable DB), 5 facets: physics class; IC/BC+geometry; related applications; TYPE (1 analytical, 2 manufactured, 3 ODE numerical, 4 PDE numerical); features tested.
  Type 1-2 => observed order of accuracy of candidate computable; type 3-4 => "doubtful", only SRQ accuracy comparison.
- p.81-82 (§3.1.2): math description must NOT include any feature of discretization; a) assumptions b) symbols/units c) continuum equations incl. submodels ("neutron cross-section models") d) IC/BC in continuum form "actually used" e) SRQs (dependent vars, functionals, stated in continuum form) f) uncertain inputs as distributions. "unambiguous, reproducible ... must be ruthlessly pursued".
- p.83-85 (§3.1.3): accuracy assessed PER SRQ, as function of space/time/parameters, "definitive pedigree". Type 1: series truncation error cannot be estimated from one more term; integrals/transcendental roots need numerical error estimate. Type 3: integrator order >= 3-4 orders above candidate's formal order; two integrators. Type 4: iterative criteria, >=3 meshes + Richardson, observed order, two markedly different methods near singularities.
- p.86-87 (§3.1.4): documentation: hardware, OS, compiler+options, ARITHMETIC PRECISION, language, run time, authorship; Type 2 source terms in machine-copyable code form + symbolic package version; Type 3 package+version + its own verification record; Type 4 code+version.
- p.88-89 (§3.2): comparisons NOT stored in the benchmark DB (contra Rizzi-Vos). Fixed tolerance "quite arbitrary"; "To say that the accuracy required depends on the application ... defeats the purpose". Pass/fail differs by type: type 1-2 => observed order matches formal (strong) or >0 (weak); else difference vs resolution plot.
- p.96 (§5): DB needs review/approval, configuration management; no detail on versioning/regeneration.
ACCURACY-GRADE: implicit criterion = the benchmark must be accurate enough to measure the candidate's observed order of accuracy (p.88); type hierarchy is the grade.

## 2. Briggs, Scott & Nouri (2003), "The International Criticality Safety Benchmark Evaluation Project", NSE 145(1):1-10, doi:10.13182/NSE03-14. [M] read rendered pp.2-7 (PDF p = printed p + 1)
Local: scratch/literature/Briggs-Scott-Nouri(2003)...pdf (copied from ~/Downloads/NSE Vol_145(1) zip). Companion Dean (2003) NSE 145:20-38 doi:10.13182/NSE03-16 also copied (not yet read).
- p.3: revisions happen; "a revision history is maintained and published"; citations must include the EDITION.
- p.3, Table II p.5: identifier = (Fissile material)-(Physical form)-(Spectrum)-(NNN), SUB- prefix for subcritical: a structured, faceted, immutable ID.
- p.4: five sections + appendix. 1.0 Detailed Description (experiment as performed, sources: logbooks, memos; enough that derivation of Sec 3 is evident); 1.1-1.4. 2.0 Evaluation of Experimental Data (uncertainties quantified on keff; ACCEPTABILITY decided here; unacceptable data are not carried into 3.0/4.0/App A).
- p.6: 3.0 Benchmark Specifications = data to construct a calculational model: 3.1 model + justified simplifications, 3.2 dimensions, 3.3 material atom densities, 3.4 temperatures, 3.5 experimental and BENCHMARK-MODEL keff with uncertainty (adjusted for simplifications). 4.0 Results of Sample Calculations (codes x cross-section sets, Table III). 5.0 References.
- p.7: App. A Typical Input Listings (code version, library, quadrature/scattering order, convergence, mesh; histories). Sample calcs and inputs are NOT peer-reviewed like the rest and are "not ... a validation of the codes". Three-level peer review (internal / independent / working group) checks spec derivable from description, completeness, format.
- p.7 (§V): acceptance grade: uncertainty in calculated keff > 1% "often judged to be unacceptable" as a benchmark. A graded admission criterion (acceptable vs unacceptable) with a number.

## 3. Ganapol (2008), Analytical Benchmarks for Nuclear Engineering Applications: Case Studies in Neutron Transport Theory, NEA/DB/DOC(2008)1, NEA No. 6292, ISBN 978-92-64-99056-2. [M] read (pdftotext + rendered pp.20-21 spot-checked)
Local: scratch/literature/Ganapol(2008)Analytical Benchmarks for Nuclear Engineering Applications.pdf (free, oecd-nea.org db-doc2008-1.pdf). PDF pages below.
- PDF p.13 (Preface xiii): "benchmark quality ... transport theorists later defined as four- or five-place accuracy".
- PDF p.15 (xv): per-benchmark anatomy: reactor-physics context -> classification (Table A) -> general description -> self-contained derivation of the solution representation -> numerical algorithms -> demonstration/trends -> "BENCHMARK QUALIFICATION" table: convergence of flux vs requested error, "to give the reader confidence that the benchmark is indeed accurate to the digits quoted". FORTRAN source + input for tables distributed (App. D = I/O description).
- PDF p.20 (xx): semi-analytical benchmark = "accurate numerical evaluation (usually to four or five correct digits) of an analytical solution representation"; its error is MONITORED ("numerical accountability") vs unmonitored discretization error of numerical benchmarks.
- PDF p.21 (xxi): hierarchy exact-analytical > near-analytical > semi-analytical > purely-numerical; ultra-fine-mesh SN/Pn can match semi-analytical accuracy (Ref [1] Ganapol M&C 2005).
- PDF p.23 (xxiii) Table A: 8-facet classification C1 field, C2 geometry (+qualifiers I/H/2H/S/HE), C3 scattering anisotropy, C4 energy (OG/MG/C), C5 angular source, C6 spatial/temporal source, C7 continuous/discrete treatment of X,A,T,E, C8 METHOD (NLTI, FN, IT, GFM, ...). i.e. problem facets C1-C6 separated from reference-method facets C7-C8.
- PDF p.69: convergence always judged by an "engineering estimate" = relative error between the last two iterates; failures must be DETECTED AND FLAGGED (HT, HR, HJ error reports).
- PDF p.126 Table 3.1.1: qualification = same observable tabulated at eps=1e-3..1e-6; "last column correct to all digits posted".
ACCURACY GRADE: "benchmark quality" = 4-5 correct digits, demonstrated by a qualification table (requested-error ladder), not by assertion. NB Oberkampf-Trucano p.83 warns change-between-last-terms is not a sufficient error estimate for series: tension between the two sources.

## 4. ANL-7416 Supplement 2, Argonne Code Center: Benchmark Problem Book, ANS Math & Comp Division Computational Benchmark Problems Committee, revised June 1977. [M] rendered pages of an excerpt (BSS-15 only; printed p = PDF p + 584)
Local: scratch/literature/ANL-7416-BSS15(corephysics)...pdf (corephysics.com mirror; scanned, no text layer; OCR failed 429 rate limit). OSTI full copies (purl 4498920 = 1968 original; 5037820 = Supp.2) unreachable from this host (connection reset).
- PDF p.3 TOC: III "Guidelines and Format" A Source Situations, B Format (printed pp.8-9) -- NOT in the excerpt: its normative text is UNREAD.
- THREE-LEVEL record, each with its own ID, submitter, date submitted, date accepted/adopted, reviewers:
  (a) PDF p.6 (p.590) BENCHMARK SOURCE SITUATION, ID 15: Descriptive Title, Suggested Function, Configuration, Details (physical situation + data).
  (b) PDF p.7-8 (pp.591-592) BENCHMARK PROBLEM, ID 15-A1, "Source Situation ID.15": Descriptive Title, REDUCTION FROM SOURCE SITUATION (the mathematical model dN/dt = A N and all data tables), EXPECTED PRIMARY RESULTS (observables: 50-day concentrations; "calculational statistics"), SOLUTIONS AVAILABLE (list of solution IDs 15-A1-1..3 with one-line technique).
  (c) PDF p.14-16 (pp.598-600) BENCHMARK PROBLEM SOLUTION, ID 15-A1-1, "Benchmark Problem ID.15-A1": Technique Used (RK-Gill 4th order, time-step rule), Computer, Program, References, Primary Results (Table C.1 at several step sizes, CPU time).
- So: one physical situation -> many mathematical problems -> many solutions per problem; hierarchical IDs; no per-solution accuracy grade seen in this record (only a step-size ladder).

## 5. Sood, Forster & Parsons (1999) LA-13511 (local sidecar) [M] sidecar only
- sidecar p.11-12, Table 1: faceted problem identifier Material-Groups-Scattering-Geometry (+ reflector form), e.g. U-2-0-SP. Same faceting idea as ICSBEP Table II and Ganapol Table A.

## 6. Software: content-addressed caching
### Mokhov, Mitchell & Peyton Jones (2020), "Build systems a la carte: Theory and practice", J. Funct. Prog. 30:e11, doi:10.1017/S0956796820000088 [M] pdftotext, PDF pages
Local: scratch/literature/Mokhov-Mitchell-PeytonJones(2020)...pdf (author copy, ndmitchell.com).
- p.8 (§2.4): Bazel = content-addressable cache (hash(content)->content) + memo table of executed commands with input/output hashes; predicts result hash from dependency hashes.
- p.23 (§5.2) verifying trace: store hashes of deps; if unchanged, skip. Stores no value.
- p.26 (§5.3) constructive trace: also stores the VALUE, may hold many per key, shareable ("cloud").
- p.27 (§5.4) deep constructive trace (Nix, Buck): key on TERMINAL inputs only; disadvantages: tasks MUST be deterministic, and NO early cutoff.
- p.43 (§8.3) impurity: untracked deps (compiler version), non-determinism, volatility; sandboxing never captures "down to the CPU model and microcode".
- p.48 (§8.8) self-tracking: when tasks are functions in a full language, equality of tasks is undecidable; "the only safe approach is to assume (pessimistically) that any change to the build system potentially changes any build task".
### Dolstra (2006), The Purely Functional Software Deployment Model, PhD thesis, Utrecht [M] pdftotext
Local: scratch/literature/Dolstra(2006)...pdf.
- PDF p.27 (§2.1): store path prefix = "cryptographic hash of all inputs involved in building the component" (input-addressed; extensional model, Ch.5).
- PDF p.106 (Ch.5): content-addressability only for sources in extensional model; Ch.6 intensional model = content-addressed for all objects; §6.4.1 equivalence-class collisions (same inputs, different outputs from non-deterministic builds).
### pytest-regressions docs (overview page) [M] via WebFetch summary only
- reference files in a committed data dir, named by test id; missing file => test FAILS and creates it; --force-regen / --regen-all regenerate explicitly. Key = test name, NOT generator content: no staleness detection. Tolerance options not confirmed.

## 7. Not obtained
- Oberkampf & Roy (2010) V&V in Scientific Computing, CUP: paywalled book, not local.
- Knupp & Salari (2003) Verification of Computer Codes in Computational Science and Engineering, CRC: paywalled, not local.
- Roache (1998) V&V in Computational Science and Engineering (Hermosa); Roache (2002) J. Fluids Eng. 124(1):4-10 "Code verification by MMS" [R: web only]: paywalled.
- Salari & Knupp (2000) SAND2000-1444, doi:10.2172/759450: FREE on OSTI but OSTI resets connections from this host; UNT mirror behind bot challenge.
- ANL-7416 §III "Guidelines and Format" (printed pp.8-9): OSTI unreachable; excerpt lacks it.
- Oberkampf-Trucano-Hirsch (2004) Appl. Mech. Rev. 57:345 (origin of SSB four characteristics): not fetched.
