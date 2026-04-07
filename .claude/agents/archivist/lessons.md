# Archivist Lessons

Read this at the START of every invocation.

---

## L1: Build Sphinx first, always

`-W` (warnings-as-errors) catches structural breaks. But zero warnings
does NOT mean cross-refs resolve — Sphinx silently renders unresolvable
`:func:` as plain text. Check whether `nitpicky = True` is in conf.py.
If not, manually verify cross-references against code.

## L2: Duplicate labels are silent

Sphinx does not warn on duplicate `:label:` across RST files. Always
run `grep :label: docs/**/*.rst | sort | uniq -d` before building.
Page-prefix all labels (e.g., `moc-keff-update`, `cp-keff-update`).

## L3: Dead citations recur

Bibliography entries with no body-text citation are common. Check:
grep `[Name20XX]_` in body vs bibliography. Trivially automatable.

## L4: Private functions → code literals, not `:func:`

`:func:`_private_name`` won't resolve unless `:private-members:` is
on the automodule. Use backtick code literals for private functions.

## L5: Verify test counts with `--collect-only`

RST test-count claims go stale fast. Always run
`.venv/bin/python -m pytest --collect-only -q` before documenting
any test count.

## L6: The quality bottleneck is `VerificationCase.latex`

Generated RST quality depends on the derivation modules' LaTeX
output, not on `generate_rst.py` or the RST itself. To enrich
documentation, enrich the derivation modules.

## L7: Discretization-first for simulation modules

TH and RK pages show the continuous PDE but not the FD stencil.
Any textbook has the PDE — the code's value is the specific discrete
form. For simulation modules, derive the discrete equation first,
use the continuous PDE only as context.

## L8: Three-layer tracing

`derivations/ → generate_rst.py → RST`. Quality issues can originate
at any layer. When reviewing generated docs, trace to the generator
AND the data source. When reviewing hand-written RST, check whether
a derivation script should exist but doesn't.

## L9: Autodoc path is step 0

Before suggesting any `:func:` or `:class:` directives, verify the
target module appears in `docs/api/`. If not, all cross-references
are dead on arrival.

## Quality Scores (reference)

| Page | Score | Key gap |
|------|-------|---------|
| discrete_ordinates.rst | 3.5 | 5 equations without derivation scripts |
| verification.rst | 2.3 | Thin content at every layer |
| thermal_hydraulics.rst | 2.2 | No FD discretization, no derivation scripts |
| reactor_kinetics.rst | 1.7 | Same as TH |
| cross_section_data.rst | 2.5 | No derivation script, broken cross-refs |
| collision_probability.rst | 3.7 | No GS convergence comparison table |
| monte_carlo.rst | 4.0 | Autodoc gap (35 non-resolving directives) |
| method_of_characteristics.rst | 8.3 | Core equations hand-written |
