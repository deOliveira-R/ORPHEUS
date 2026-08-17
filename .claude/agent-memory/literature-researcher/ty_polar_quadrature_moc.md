---
name: ty-polar-quadrature-moc
description: Tabuchi-Yamamoto (TY) MOC polar quadrature — verified refs (JNST 44(2):129-136 2007 + Knott-Yamamoto Handbook ch.9 2010), the Leonard-McDaniel ancestor in LOCAL Hebert Ch3, the decoded weight convention, and the four facts that need the paper body
metadata:
  type: reference
---

# TY polar quadrature (MOC) — extraction record

Pulled 2026-08-11 for `orpheus/moc/quadrature.py::_TY_TABLES`. Deliverable:
`scratch/literature_yamamoto_ty_quadrature.md`.

## NOT in the local library. Bronze OA — free to read, we are only bot-blocked.

- **Yamamoto, Tabuchi, Sugimura, Ushio, Mori** (5 authors), "Derivation of
  Optimum Polar Angle Quadrature Set for the Method of Characteristics Based
  on Approximation Error for the Bickley Function", *J. Nucl. Sci. Technol.*
  **44**(2):129-136 (2007). DOI `10.3327/jnst.44.129` **and**
  `10.1080/18811248.2007.9711266` — same paper, the first redirects to the
  second. CrossRef-verified on both; CiNii crid `1361137043839127552`.
- **Knott, D. & Yamamoto, A.** (2010) "Lattice Physics Computations",
  ch. 9 of *Handbook of Nuclear Engineering*, Springer US, pp. **913-1239**,
  DOI `10.1007/978-0-387-98149-9_9`. CrossRef `editor: null` — Cacuci is
  UNVERIFIED. Not a journal article.
- **No earlier Tabuchi paper exists.** The hyphenated name is coined IN the
  2007 abstract ("Tabuchi-Yamamoto or the TY quadrature set"); Tabuchi is its
  own 2nd author. Genuine ancestor = **Leonard & McDaniel** (full citation
  still unresolved — Hébert's Ref. [51], and his Ch. 3 scan omits the
  reference list).

⚠ **Access map (do not re-attempt):** J-STAGE 404s on all `article/jnst/44/*`
paths — JNST 1964-2011 was MIGRATED OFF J-STAGE to Taylor & Francis. T&F is
Cloudflare-403 to curl+WebFetch alike, yet Unpaywall says
`bronze / publisher / publishedVersion` ⇒ **free in a browser, not paywalled**.
OpenAlex `any_repository_has_fulltext: false`, one location — no green copy
exists. J-GLOBAL 503s. HAL/OSTI are the wrong corpora (Japanese AESJ work).
⇒ OpenAlex `is_open_access: true` with `?needAccess=true` in the `oa_url` is
the **bronze-OA signature**: free-to-read but bot-blocked, NOT a false OA flag.

## What the abstract gives (VERIFIED) vs what needs the body (UNKNOWN)

VERIFIED from the publisher abstract (reproduced independently by OpenAlex on
both DOIs *and* the CiNii `description` field): objective is **MINIMAX** ("minimizing
the maximum approximation error"); rationale is the MOC≡CP equivalence for 2-D
isotropic scattering; **L = 1, 2, 3** divisions; validated on **C5G7 + 4-loop
PWR whole core**.

⛔ UNKNOWN — abstract is silent, body not read: the **Bickley ORDER** (Ki₂?
Ki₃? joint?), the **x/τ RANGE of the fit**, any **accuracy figure or bound**,
and the **table number**. Do not ship a bound from secondary sources.

## The LOCAL ancestor + the decoded weight convention

`Hebert(2009)Chapter3.pdf` PDF **pp. 88-90** (printed 154-156): Eqs.
**(3.495)-(3.496)** = the **Leonard-McDaniel** family, target named explicitly
as **Ki₃**, with a Matlab `lmcd(nmu)` table of `(z, W)` for nmu = 2, 3, 4,
parameterized by `z = 1/sin θ`. NOT the TY digits.

`Sanchez(2002)...Boundary Conditions` NSE 140 PDF **p. 19** (printed 41,
§VI.A.1): the family's **L² (least-squares)** functional with `x_max` and `n`
("usually 2 or 3") as fixed inputs, constraints `w>0, Σw=1, θ∈(0,π/2)`.
⇒ the family is pinned by **(norm, order n, x-range, constraints)** — the four
things a re-derivation must fix. Sanchez ≠ TY (L² vs minimax); do not transfer
his `x_max`. He also documents the family's hazard: Bickley-optimal sets
integrate `A₂₀` badly ⇒ conservation risk under anisotropic scattering.

⭐ **Convention decoded via the LEONARD pair (both sides published, so this is
verification not inference).** Hébert `lmcd(3)` vs OpenMOC `LeonardPolarQuad`:
`1/z` = OpenMOC `sin θ` to all 6 digits, and `W·z = W/sin θ` = OpenMOC weight
to all 6 digits with `Σ = 1.000000` exactly. Since `∫₀^{π/2} sinθ dθ = 1`,
**the tabulated weights are weights against the measure `sinθ dθ`, total mass
1** — hence the unit sum, and 0.5 after the codes' `/2`. Equivalently
`Ki₃(x) ≈ Σ ŵᵢ sinθᵢ e^{−x/sinθᵢ} = Σ W_λ e^{−z_λ x}`.
⚠ Proven for Leonard only; that TY's *published* table shares it is inference.

## `_TY_TABLES` provenance — the docstring's "Table 2" is probably mis-attributed

The ORPHEUS table is **bit-identical to OpenMOC `TYPolarQuad`** (`sin θ`
0.798184 / 0.363900,0.899900 / 0.166648,0.537707,0.932954; weights
1.0 / 0.212854,0.787146 / 0.046233,0.283619,0.670148, all halved by `/2.0`),
and OpenMOC's `src/Quadrature.h` header cites **only Knott & Yamamoto (2010)**,
never the 2007 paper. ⇒ split the roles: **2007 = derivation, 2010 =
tabulation**. Classic L-003.

## Re-derivation feasibility: measurable in-tree, no literature dependency

Oracle already ships: `orpheus/derivations/common/kernels.py::ki_n` (mpmath,
50 dps, any `n≥1`, `x≥0`) — `[M]` its `cos`-form equals the literature
`sin`-form to `0.00e+00` rel. `ki_n_float` ≈1e-14 for hot loops. Amos 1983
(LOCAL) is an adequate independent oracle (`n≥0`, 18 digits, tested `0≤x≤50`).
⛔ **Lorensi-Azevedo-Sauter 2025 (LOCAL) covers Ki₁ and Ki₂ ONLY** — cannot
serve as the Ki₃ oracle.

`[M]` The stored table vs `ki_n(3,·)`: max\|r\| = 1.279e-02 (L=1), 4.123e-04
(L=2), 3.361e-05 (L=3); the L=3 residual **equioscillates** at ≈3.35e-05 over
6 alternations spanning `x∈[0,~5]` then decays — the Chebyshev signature of an
absolute-error minimax on Ki₃. **Diagnostic only, NOT a paper finding.** It
makes the re-derivation a one-experiment question rather than a literature
dependency.
