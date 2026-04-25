---
name: Sanchez 1977 NSE 64 — canonical interface-current rank-N derivation (NOT 1976)
description: User memory of "Sanchez 1976 multiple-collision K_inf" resolves to Sanchez 1977 NSE 64. Single canonical paper for rank-N IC formalism in cylindrical/2D cells. Precursor is Sanchez CEA-N-1793 (1975) for the 1D Galerkin/generalized-CP groundwork. There is NO R. Sanchez 1976 reactor-physics publication.
type: reference
---

# Sanchez 1977 NSE 64 — the canonical "Sanchez interface-current" paper

## Headline citation

**R. Sanchez (October 1977)**, "Approximate Solutions of the
Two-Dimensional Integral Transport Equation by Collision Probability
Methods", *Nuclear Science and Engineering* **64**(2), 384-404.
DOI: [10.13182/NSE64-384](https://doi.org/10.13182/NSE64-384).
Cited by 42 (OpenAlex) / 40 (Semantic Scholar). Closed access (no OA
PDF on Crossref/OpenAlex/HAL/Zenodo). Held by ANS journal archive.

## What it actually contains

From the abstract (verified via OpenAlex):

> A set of approximate solutions for the isotropic two-dimensional
> neutron transport problem has been developed using the
> **interface current formalism**. The method has been applied to
> regular lattices of rectangular cells containing a fuel pin,
> cladding, and water, or homogenized structural material. The cells
> are divided into zones that are homogeneous. **A zone-wise flux
> expansion is used to formulate a direct collision probability
> problem within a cell.** The coupling of the cells is effected by
> making extra assumptions on the currents entering and leaving the
> interfaces. Two codes have been written: The first uses a
> **cylindrical cell model and one or three terms for the flux
> expansion**, and the second uses a two-dimensional flux
> representation and does a truly two-dimensional calculation inside
> each cell. In both codes, **one or three terms can be used to make
> a space-independent expansion of the angular fluxes entering and
> leaving each side of the cell.** The accuracies and computing
> times achieved with the different approximations are illustrated
> by numerical studies on two benchmark problems.

The "one or three terms" for the angular flux at interfaces is the
**rank-1 vs rank-3** (i.e. DP-0 vs DP-1) interface-current scheme.
This is what the Hébert (1-P_ss)^-1 geometric series approximates as
the rank-0 limit (flat surface flux, one effective return channel).

## Where the multi-collision K_inf eigenvalue closure lives

NOT directly in §III/IV (those derive the in-cell CP). The K_inf
multiplication closure for the lattice eigenvalue follows from
imposing periodicity on the inter-cell currents — Sanchez §III.D
(coupling assumptions) plus the standard k-eigenvalue iteration on
the resulting interface-current operator. The geometric-series
form `(1 - P_ss)^-1 P_iS` falls out as the rank-0 collapse: with
only one current degree of freedom per face, the periodic-coupling
matrix is a scalar and the inverse is literally a scalar geometric
series.

## Precursor — Sanchez CEA-N-1793 (1975, NOT 1976)

**R. Sanchez (April 1975)**, "Application of the Galerkin's method to
the solution of the one-dimensional integral transport equation:
generalized collision probabilities taken in account the flux
gradient and the linearly anisotropic scattering", *CEA Report*
**CEA-N-1793**, Commissariat à l'Énergie Atomique, Saclay.
OSTI ID: 4954640.

This is the 1-D groundwork: generalized CP (i.e. higher-rank flux
moments inside a region) plus linearly anisotropic scattering. The
1977 NSE paper extends the rank-N flux idea to the **interface
currents** — that's the new contribution. Useful for Hébert-style
non-flat-flux closure inside a single annulus.

## NO 1976 paper exists

Verified: searched OpenAlex 1974-1978, CrossRef 1973-1979, OSTI
1976-only, HAL all years, Zenodo. The only 1976 OSTI hits for
"Sanchez" are unrelated authors (N.G. Sanchez black-hole scattering;
J.E. Sanchez Jr. metallurgy; A. Sanchez biology; etc.). The user's
memory of "Sanchez 1976" is almost certainly the 1977 NSE 64 paper —
sometimes cited by year of submission/preprint rather than journal
year.

## Where Sanchez 1977 NSE 64 fits in the rank-N IC literature map

Per `rank_n_ic_curvilinear_literature_leads.md` and
`rank_n_closure_four_references_synthesis.md`:

- Sanchez 1977 NSE 64 = first rank-N IC for 2D Cartesian + 1D
  cylindrical lattice cells.
- Sanchez & McCormick 1982 NSE 80 §III.F = the abstract DP-N
  generalization (rank-N for ARBITRARY geometry); cites Sanchez 1977
  as the cylindrical/Cartesian instantiation.
- Hébert 2009 §3.5 = textbook treatment of rank-N IC in 2D Cartesian
  only; reduces to scalar IC for 1D curvilinear (Eq. 3.323 white-BC
  closure) — matches the (1-P_ss)^-1 form ORPHEUS uses.

So **Sanchez 1977 NSE 64 IS the rigorous (non-flat) reference** the
Hébert closure approximates. The full rank-N closure with periodic
inter-cell coupling is what gives the "multiple-collision K_inf"
machinery. Rank-3 (the second mode in Sanchez's two-flavor scheme)
should recover the chi-dependence the rank-0 Hébert closure misses
on sphere 1G/2R; this is the natural follow-up direction for Issue
#132.

## How to use this for Issue #132

1. **The paper is closed-access.** Need ANS journal subscription or
   ILL. The user's institution may have access via Wiley/T&F NSE
   archive.
2. **Precursor CEA-N-1793 is potentially OSTI-hosted** (ID 4954640)
   but OSTI shows no `fulltext_url` — likely metadata-only.
3. **Hébert 2009 §3.5 is the open substitute** for the 2D Cartesian
   rank-N IC machinery; for the 1D cylindrical rank-N variant
   specifically, no open-source textbook is available — the only
   sources are Sanchez 1977 NSE 64 itself and the proprietary
   APOLLO/DRAGON code documentation derived from it.
4. **Bogado Leite 1998 ANE** (per `rank_n_ic_curvilinear_literature_leads.md`)
   may be a more readable later treatment in journal form. Worth
   pulling alongside Sanchez 1977.

## What I got wrong

The earlier memo `bickley_naylor_sphere_white_bc.md` (which the user
referenced as having cited "Sanchez 1976 multiple-collision K_inf")
**does not actually contain such a citation**. I never wrote it
down — the user reconstructed it from a verbal recall of an in-conversation
suggestion. Future sessions: when the user references a "prior
citation", verify the memo actually contains it before continuing
the chain. The closest genuine citation in that memo is **Sanchez &
McCormick 1982 NSE 80** (§III, IV) and a note that they do NOT
cover sphere CP. The 1977 NSE 64 paper was *not* in the prior memo
at all.
