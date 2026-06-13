---
name: spherical-sn-central-cell-spatial-order
description: Canonical literature on spherical SN central-cell (r=0) spatial accuracy, starting direction, cell-center radius, and MMS reference quantity. Hebert 2009 3.9.4 + Stacey 2007 9.9 are the two complete textbook treatments; Bailey-Morel-Chang 2010 is ANGULAR not spatial.
metadata:
  type: project
---

Diagnosing O(h) spatial defect at r=0 central cell in 1-D spherical SN (WDD/DD + Morel-Montry Carlson redistribution). Findings from local PDFs (scratch/literature/).

**Two complete textbook treatments exist; both use plain arithmetic-midpoint diamond at the central cell. NEITHER prescribes a special O(h2) central-cell closure.**

- **Hebert (2009) Applied Reactor Physics §3.9.4** (book pp.142-144, PDF pages idx 74-77). THE canonical curvilinear SN reference.
  - Eq. (3.428) FV balance; Eq. (3.430) DEFINES mesh-centered value as the VOLUME-AVERAGE phi_{n,i} = (4pi/V_i) integral r^2 phi_n dr — mass-weighted, NOT a point value. SAME for cylinder Eq. (3.405) with rho weight.
  - Eq. (3.431) diamond closure phi_{n,i} = 1/2(phi_{n,i-1/2}+phi_{n,i+1/2}) = 1/2(phi_{n-1/2,i}+phi_{n+1/2,i}) — arithmetic mean of BOTH radial AND angular edges.
  - Eqs. (3.432)-(3.435): starting-direction mu=-1 init sweep. Angular redistribution VANISHES (alpha term gone), reduces to a slab-like difference. r_{1/2}=0 stated explicitly above Eq. (3.427).
  - Eq. (3.424) alpha recursion alpha_{n+1/2}=alpha_{n-1/2}-2 W_n mu_n; alpha_1=0.
  - NO convergence-order statement anywhere in §3.9; "first-order" only ever refers to the first-order FORM of the Boltzmann PDE.

- **Stacey (2007) Nuclear Reactor Physics §9.9** (book pp.359-361, PDF pages idx 55-57). Second complete spherical-SN treatment.
  - Eq. (9.218) starting-direction mu=-1 (slab-like, no redistribution); Eq. (9.219) inward sweep; Eq. (9.220) angular diamond; Eq. (9.221) CENTER handled by SYMMETRY condition psi_{N+1-n}^{1/2}=psi_n^{1/2} (NOT a special spatial closure).
  - Eq. (9.205/9.220) arithmetic-midpoint diamond. Eq. (9.210) slab DD error O((Sigma_t Delta/2|mu|)^2). 2-D Cartesian Eq. (9.231) defines volume-averaged flux.
  - NO statement of order degradation at r=0.

- **Bailey-Morel-Chang (2010) NSE 165(2):149-169.** CONCERNS ANGULAR DIFFERENCING DIFFUSION-LIMIT, NOT SPATIAL. The "flux dip / negative slope at the origin" they analyze (Figs. 4-12) is the beta-contamination term in the SECOND-ORDER CURRENT (angular-redistribution artifact), Eq. (40) spherical / (73) cylindrical, beta_spherical Eq. (after 73). MM-WDD forces beta=0; step/diamond have nonzero beta. They EXPLICITLY spatially-converged all tests ("We ensured all numerical solutions are spatially converged") to ISOLATE angular error. DO NOT conflate their origin flux-dip with the spatial volume-weighting O(h) issue.

**Answers:**
1. WDD/DD at r=0 is NOT given a special O(h2) treatment in literature. The cell-AVERAGE vs midpoint-POINT-value gap is O(h) under r^2 weighting and is a COMPARISON subtlety — the cell-average IS the scheme's native unknown (Hebert 3.430). No source frames it as scheme truncation error.
2. Starting direction = mu=-1 inward sweep where redistribution vanishes (Hebert 3.432-3.435 / Stacey 9.218). Center via symmetry (Stacey 9.221), NO parabolic/L'Hopital reconstruction in either text.
3. Hebert uses VOLUME-AVERAGE (3.430), implicitly the centroid weighting. Neither text prescribes an "effective radius" point-eval to restore O(h2). The volume-centroid r_bar = (3/4)(r4_out-r4_in)/(r3_out-r3_in) is the mathematically-correct centroid but is NOT named in Hebert/Stacey.
4. MMS reference quantity SHOULD be the cell-VOLUME-AVERAGE integral phi r^2 dr / integral r^2 dr (matches Hebert 3.430 unknown), NOT phi(midpoint). Comparing cell-avg vs midpoint point-value injects a spurious O(h) at r=0. This is the most likely root of an observed O(h) MMS defect.
5. NO explicit convergence-order statement for the pole cell in Hebert/Stacey. Tier-2 (OpenAlex/web) found NO paper on a special O(h2) spherical-SN central-cell closure — strong evidence the plain DD + starting-direction IS the accepted standard.

Lewis & Miller (1984) NOT in scratch/literature/. Zotero MCP NOT exposed this session (see [[reference-zotero-flakiness]]).
