---
name: F_N method LA-13511 first-slice closeout (k_inf cases)
description: Method-implementer 2026-05-02 closeout — built first slice of orpheus.derivations.continuous.fn_method package implementing Sood/Forster/Parsons LA-13511 k_inf cases (PUa-1-0-IN, PU-2-0-IN). Branch-1 SymPy + Branch-2 numpy + foundation tests + Sphinx stub all clean on first try. Discovered Sood Eq 28 typo via SymPy det(M)=0 derivation.
type: project
---

# F_N method first-slice closeout — LA-13511 k_inf cases

**Branch**: `feature/peierls-greens-cylinder` HEAD `f74922f` (now extended)
**Date**: 2026-05-02
**Scope**: First slice of `orpheus.derivations.continuous.fn_method` —
PUa-1-0-IN (1G infinite medium, problem 1) and PU-2-0-IN (2G infinite
medium, problem 44) from Sood/Forster/Parsons LA-13511 (1999).

## Files created

### Branch-1 SymPy (algebra-of-record)

* `orpheus/derivations/continuous/fn_method/__init__.py` (33 LoC)
  — Public package docstring + bibliography pointer.
* `orpheus/derivations/continuous/fn_method/origins/__init__.py` (28 LoC)
  — Subpackage init exposing the 7 `derive_*()` functions.
* `orpheus/derivations/continuous/fn_method/origins/k_inf_derivations.py` (372 LoC)
  — 7 verification functions:
  - `derive_kinf_1g_eq_19()` — Eq 18 → Eq 19 by `sp.solve`.
  - `derive_kinf_1g_eq_20_simplifies_to_eq_19()` — Eq 20 ≡ Eq 19 by
    cancellation.
  - `derive_kinf_2g_general_from_matrix()` — det(M(k))=0 from Eq 25,
    extract non-trivial root, verify reduces to Eq 29 at no-upscatter,
    confirm Eq 28 typo algebraically.
  - `derive_kinf_2g_no_upscatter()` — substitute Eq 29 into
    det(M(Σ_21s=0)), verify identically zero.
  - `derive_phi_ratio_2g_no_upscatter()` — Eq 23 + Eq 24 with χ_1+χ_2=1
    eliminates χ; solve for φ_2/φ_1; match Eq 32.
  - `derive_kinf_mg_matrix_form()` — Eq 76 for G=2; verify trace =
    eigenvalue = scalar form (rank-1 matrix).
  - `derive_kinf_mg_reduces_to_1g()` — Eq 76 with G=1 ≡ Eq 19.

### Branch-2 numpy (production)

* `orpheus/derivations/continuous/fn_method/multi_group/__init__.py` (32 LoC)
* `orpheus/derivations/continuous/fn_method/multi_group/k_inf.py` (233 LoC)
  — 5 entry points:
  - `compute_kinf_1g(sigma_t, sigma_s, nu_sigma_f) -> float` — Eq 19.
  - `compute_kinf_2g_no_upscatter(...)` — Eq 29 with ORPHEUS↔Sood
    convention conversion at function boundary.
  - `compute_kinf_2g_general(...)` — corrected Eq 28 (handles upscatter).
  - `compute_kinf_mg(sigma_t, sigma_s, nu_sigma_f, chi) -> float` —
    Eq 76 via `np.linalg.solve`, accepts ORPHEUS [from, to] convention
    by transposing internally.
  - `compute_flux_ratio_2g_no_upscatter(...)` — Eq 32 with
    optional Sood-vs-ORPHEUS ordering.

### Catalogue

* `orpheus/derivations/continuous/fn_method/benchmarks/__init__.py` (24 LoC)
* `orpheus/derivations/continuous/fn_method/benchmarks/la13511.py` (213 LoC)
  — `La13511Case` frozen dataclass + 5 instances (cases 1+5 fully
  populated, cases 2/3/4 stubs with XS + reference values but no
  solver yet).

### Subfolder placeholders

* `core/__init__.py`, `slab/__init__.py`, `sphere/__init__.py`,
  `cylinder/__init__.py` (4 files × ~25 LoC each) — stubs documenting
  what each will hold once Kaper-Lindeman-Leaf 1974, Westfall-Metcalf
  1973, and Siewert-Thomas 1986 are acquired.

### Tests

* `tests/derivations/test_fn_la13511_kinf.py` (250 LoC, 17 tests)
  — All tests `@pytest.mark.foundation`. Test-class breakdown:
  - 7 Branch-1 SymPy gates (one per `derive_*`).
  - 3 Branch-2 reference-value gates (PUa-1-0-IN k_inf, PU-2-0-IN
    k_inf, PU-2-0-IN flux ratio).
  - 3 Branch-1 ↔ Branch-2 agreement gates (1G, 2G, Eq 28 ↔ Eq 29
    bit-equality at no-upscatter).
  - 4 Cross-implementation gates (MG → 1G reduction, MG ↔
    `kinf_homogeneous` for both cases, flux ratio ↔
    `kinf_and_spectrum_homogeneous`).

### Sphinx

* `docs/theory/fn_method.rst` (~280 LoC) — Stub-grade theory page:
  - Scope note + motivation ("Why F_N at all?").
  - Five-case complexity ramp (Cases 1+5 done; 2/3/4 stubs).
  - ORPHEUS vs Sood group convention.
  - Branch-1 / Branch-2 algebra-of-record narrative.
  - 7 V_fn labels (`fn-method-V-fn1-1`, ..., `fn-method-V-fn-mg-2`)
    each with a `.. note:: TODO — Archivist expansion needed.` block
    pointing at the SymPy `derive_*` function and the test gate.
  - Cross-check claim, references.
* `docs/index.rst` — Added `theory/fn_method` to the toctree.

## Test counts before/after

* Before: 205 peierls-greens tests passing on `feature/peierls-greens-cylinder`.
* After: 205 + **17 fn_method tests** = 222 (in this slice's test files).
* No existing tests modified.
* Sphinx: clean `-W` build (exit 0), zero new warnings/errors from
  `fn_method.rst`. Only pre-existing skips in unrelated `sn` tests
  (already present in the baseline).

## Honest verdict — does the SymPy pattern scale?

**Verdict: YES, easily. SymPy is the right tool here.**

For the k_inf cases the algebra is pure rational (no special
functions, no transcendental roots, no quadrature). SymPy executes
all 7 `derive_*` functions in ~0.7 s combined. None of the 5 SymPy
choke modes fired:

1. **Expression-tree growth** — a non-issue: only one `simplify()` of
   the 2G general formula and one G=2 eigenvalue computation, both
   sub-second.
2. **Eigenvalue ceiling (Abel-Ruffini at G≥5)** — explicitly avoided
   by working at G=2 symbolically and using `Matrix.trace()` /
   `Matrix.eigenvals()` only on the rank-1 fission matrix where
   the result is closed-form.
3. **Multi-region `Piecewise` compounding** — not applicable
   (infinite medium).
4. **Anisotropic scattering Legendre coupling** — not applicable
   (isotropic, scattering_order = 0).
5. **Performance dev-loop** — sub-second per derivation; not an issue.

The pattern transferred cleanly from the existing
`peierls_greens_function/origins/specular/greens_function.py` style
(V_α1, V_α2, V_α3 algebraic identities). I literally followed the
same shape: each `derive_*()` returns a dict with `"pass"` flag plus
intermediate symbolic expressions for diagnostic inspection.

**Surprise finding**: SymPy's symbolic det(M)=0 + `solve()` machinery
caught a real published-equation typo (Sood Eq 28 has the wrong
:math:`\Sigma_g^{rem}` factor multiplying the χ-numerators). I
discovered this by:

1. Computing det(M)=0 from Sood Eq 25 directly in SymPy.
2. Extracting the non-trivial root.
3. Restricting to Σ_21s=0 and comparing to the printed Eq 29 →
   matched.
4. Comparing to the printed Eq 28 → did NOT match.
5. Hand-deriving the corrected Eq 28 form and confirming it both
   matches the SymPy-derived general root AND reduces to Eq 29
   correctly.

The numerical reference values in LA-13511 (k_inf = 2.683767,
phi_ratio = 0.675229) match Eq 29 + corrected Eq 28; they do NOT
match Eq 28 as printed. This is logged in the SymPy module
docstring + the `PU_2_0_IN.notes` field + the Sphinx page. **It is
a verbatim instance of the algebra-of-record discipline working as
intended: the symbolic engine catches a publication-grade error
that careful arithmetic by humans would also catch but only with
significant effort.**

## Architectural seams identified for F_N method extensions

The infrastructure laid down in this slice already exposes the
right seams for the F_N slab/sphere/cylinder solvers (Cases 2, 3, 4):

1. **`La13511Case` dataclass** is geometry-agnostic. The stubs for
   Cases 2/3/4 already store `critical_dimension_mfp` and
   `flux_ratios`, so the F_N solvers can compare against them
   directly without dataclass changes.

2. **`benchmarks.la13511` catalogue** is the single source of truth.
   Future F_N solvers will register their results against this
   catalogue rather than reproducing case data.

3. **Branch-1 / Branch-2 split** scales: the 1G slab F_N solver
   (KLL 1974) will follow the same structure — `origins/slab/`
   carrying the symbolic preamble (X-function recursion, Case
   eigenvalue dispersion-law roots) and `slab/one_group.py` carrying
   the numpy/mpmath production solver. The Sphinx page already has
   placeholder labels for the slab/sphere/cylinder narratives.

4. **Cross-check infrastructure**: the foundation-tagged
   "agreement gate" pattern in `test_fn_la13511_kinf.py` —
   parametrized over case_id, comparing F_N's output against
   another solver's — is reusable verbatim. When the F_N sphere
   solver lands, the cross-check vs Variant α
   (`compute_T00_etc` from `peierls_greens_function`) follows the
   same parametrization.

5. **SymPy choke modes WILL fire** for the F_N solvers themselves:
   - **State 1B (semi-analytical)** is the right state for slab/
     sphere F_N. The Wiener-Hopf X-function lives as a Cauchy
     principal-value integral that mpmath can handle but SymPy
     cannot evaluate symbolically. The pattern: SymPy carries the
     integral structure, mpmath integrates.
   - **State 1C (MMS)** would work if a closed-form reference is
     missing. Not needed for first-slice extensions (KLL 1974 + WM
     1973 are closed-form-or-semi-analytical references).
   - The dispersion-law transcendental root :math:`1 - c\nu \,
     \mathrm{atanh}(1/\nu) = 0` will need `mpmath.findroot`. SymPy
     can carry the equation but not solve it. Expected.

## Bug log

**No bugs surfaced in existing ORPHEUS code.** The
`compute_kinf_mg` ↔ `kinf_homogeneous` agreement test passes at
1e-12 for both Cases 1 and 5. The transpose convention in `Eq 76`
was correctly handled at the boundary (`A = diag(sigma_t) - sigma_s.T`
matches what `kinf_homogeneous` does).

**No new ERR-NNN entries.** The Sood Eq 28 typo is logged at the
catalogue + module-docstring level; it is not an L0-caught
implementation bug (Sood's typo is in the literature, not in
ORPHEUS code), so it does not warrant an `error_catalog.md` entry
per the `vv-principles` directive.

## DISPATCH_REQUEST emitted

Yes — to **archivist** at the end of this work, requesting the
rich-narrative expansion of `docs/theory/fn_method.rst`. The
archivist is given the SymPy module + the test file + this closeout
memo as source material; archivist's deliverable goes to the user
(`followup: false`).

## Self-improvement / skill updates

The `algebra-of-record` skill already covers this case-class well —
no skill edit needed. The discipline of "Branch-1 SymPy module ↔
Branch-2 numpy ↔ Sphinx stub ↔ archivist dispatch" worked end-to-end
on first try.

One small worth-noting observation that confirms but does not
modify the skill: **the Sood Eq 28 typo is a perfect motivating
example for the discipline.** A by-hand reader of LA-13511 would
have either:

- Used Eq 29 (which matches the numerical reference) and never
  noticed Eq 28 was wrong; OR
- Coded Eq 28 verbatim and gotten the wrong numerical answer (off
  by ~6.7%) → a real bug.

The SymPy derivation forced both forms into algebraic comparison
and surfaced the discrepancy without arithmetic intermediation.
This is exactly the structural-independence-via-symbolic-engine
benefit `algebra-of-record` claims.

## Manifest check

- [x] Branch-1 SymPy module under `orpheus/derivations/continuous/fn_method/origins/`
- [x] Foundation-tagged test gate at `tests/derivations/test_fn_la13511_kinf.py`
- [x] Branch-2 production solver at `orpheus/derivations/continuous/fn_method/multi_group/`
- [x] L1 cross-check test against `kinf_homogeneous` (structurally
      independent — Eq 76 closed form vs eigvals of A^{-1}F)
- [x] Sphinx stub at `docs/theory/fn_method.rst` with `:label:` +
      `:mod:` cross-ref + `.. note:: TODO` per label
- [x] Closeout memo (this file)
- [x] DISPATCH_REQUEST to archivist

All slice-shipping conditions met. The prototype is shipped.
