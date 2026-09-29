# The cross-section grid stores reactions, and every total is derived where it is used (posing unit 3)

Status: SCHEDULED, not opened (2026-09-29). Issue #526; charter issue #522. Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; ruled content: "The ontology as it stands", "The input layer". Depends on nothing in the posing work; touches `orpheus/data/` and every consumer of the stored totals.

## Goal, in the domain's terms

A material's cross-section data is a (reaction × role) grid of reactions and their multiplicities. A total is a sum over the reactions a method actually models, so it belongs to the consumer that chose the reaction set, not to the data; the library's physical total (MT1) is a separate object, and the gap between the two is a physics fact a user must hear about. Today the total is stored, derived three times (#445), and compared with itself.

## Scope (ruled 2026-09-28)

1. **The grid stores reactions and multiplicities only**: the stored `SigT` and `SigP` retire; every sum over reactions is a derived view.
2. **The imbalanced producers re-spelled class by class**: placeholder scaffolds declare one named removal reaction (radiative capture); the billiard carrier gets a real fission cell `(Σ_f, ν)`; the Atalay encoding (scattering ratio `c_s > 1`) becomes `Σ_f > 0`, `ν = 1 + (c_s − 1)Σ_t/Σ_f`, preserving `c_s` only with `Σ_f := Σ_a`, `Σ_c := 0` (`[M]` otherwise 1.1 becomes 0.81); a pure scatterer with `c_s > 1` is a scattering cell with multiplicity above 1; `make_mixture` stops taking `sig_t`; the balance guard's own tests become the gates that the imbalanced state cannot be written.
3. **Cells named by the reaction they are**: `SigC` is radiative capture (MT102), `SigL` is (n,α) (MT107); a derived sum takes two words.
4. **Two totals**: `derived_SigT`, derived lazily by the consumer of the reaction set (the collision operator; later Monte Carlo's flight structure); `library_SigT` (MT1) at the nuclide level for the σ0 iteration. MT1 is read from the tape (`[M]` 2026-09-28: `gendf.py:476` sums the components today). ONE comparison method on the library takes a derived total, collapses the library total over the same composition, and warns that non-trivial reactions were omitted, with the size.
5. **`NuclideDensity` stores the number density**, `.mass` derived from the atomic weight ratio (if P1 of `reference_cache.md` has not already landed it).
6. Not in scope: the eager reaction-set object (Monte Carlo's campaign); the angular collision operator (#518).
7. Related issues: #445 (`Mixture.SigT` derived three times: absorbed and closed here), #518.

## Gates (from the charter's "data and the balance" list)

The cell → term `array_equal` gate against today's terms; the comparison method's warning with a positive control that omits a reaction; the k direction's negative control (scaling fission's removal share moves k, `[M]` factor 1.57 on the 2-group fixture of `scratch/posing_sequence/channels_attack_memo.md`); `derived_SigT` against the declared operators; each re-spelled producer's k or flux unchanged (bit-identical where the arithmetic allows, else the three-condition principled-equivalence of `vv-principles`).

## Done-when (hypothesis)

No stored `SigT`/`SigP` field on the grid; MT1 read; the comparison method called by the collision operator's construction; the gates green with their first reds; #445 closed.

## Opening obligations

AST census of every `SigT`/`SigP` read and write (fields, not the token), by package; read `orpheus/data/` docs and the cross-section theory pages; test-architect for the re-baselines.

## Sizing

3–4 sessions `[R]`.
