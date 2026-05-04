---
name: Method-canonical naming evidence chain — singular eigenfunction + trajectory resolvent
description: Verbatim quotes from local literature that pin the canonical method names for ORPHEUS folder taxonomy. Cite when the user asks "what does the literature call this method?"
type: reference
---

# Verbatim evidence — method-canonical names from local PDFs

All sources LOCAL in `/workspaces/ORPHEUS/scratch/literature/` —
no paywall hits required.

## "Singular eigenfunction expansion" (Case 1960 family)

The phrase is established and unambiguous across the field. NEVER
call this `case_method` — the literature itself doesn't.

- **Westfall-Metcalf 1973 (NSE 52, 1)**, p. 1, Introduction:
  *"Since the introduction of the **singular eigenfunction expansion
  technique** by Case [Ref. 1] in 1960, a wide variety of transport
  problems have been treated by this method."*
- **Atalay 1997 (Prog. Nucl. Energy 31, 229)**, abstract:
  *"Case's **singular eigenfunction method** is used to formulate the
  criticality conditions."*
- **Mitsis 1963 (ANL-6787)**, abstract:
  *"...developed by the **method of singular expansion modes**."*
- **Sanchez 1977 NSE 64**: "Case singular eigenfunction expansion".

These four papers cover slab + sphere + cylinder + collision-
probability — all four explicitly use the same method-class name.

## "Trajectory-based deterministic transport" + "scattering resolvent"

Sanchez 2002 NSE 140, 23 — the title itself names the family
"trajectory-based deterministic transport methods". Eq. (15) is the
periodic compound trajectory closure with denominator
`1/(1 − ψ_bd(L))` — algebraically identical to ORPHEUS Variant α
T(µ) = 1/(1 − α e^{−τ}).

Sanchez 1986 (TTSP 15, 333) Eq. (A4) is the resolvent
**`T(µ) = 1/[1 − α e^{−τ(µ)}]`** verbatim — the ORPHEUS code
quotes this paper.

Pomraning-Siewert 1982 (JQSRT 28, 503) Eq. (14) gives the same form
with diffuse-reflection extension: `T(µ) = [1 − α e^{−2Rµ}]^{−1}`.

NEVER call this `peierls_greens_function` for the folder name —
"Peierls" is neither in Sanchez 2002 nor Sanchez 1986 nor PS-1982;
it's ORPHEUS-internal framing. Method-canonical: `trajectory_resolvent/`.

## Author registries are NOT method names

User's rule: author registries (e.g., `sood_registry/`, `atalay1997.py`)
are ALLOWED as ENTRY-POINT NAMES because they identify a specific
*body of cases*, not a *method*. The method of record is what
USES the registry.

## How to apply

When asked "what does the literature call X" for an ORPHEUS package:
1. Check this memo for already-vetted quotes.
2. If not here, search local PDFs first
   (`ls /workspaces/ORPHEUS/scratch/literature/`).
3. Quote VERBATIM with paper + page in any naming proposal.
4. Author names → registry of cases is OK; author names → folder
   name for a *method* is NOT.
