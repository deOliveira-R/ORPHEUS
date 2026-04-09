# Archivist Lessons

Read at START of every invocation.

## Sphinx build traps

- **L1**: `-W` alone is insufficient. Sphinx silently renders unresolvable `:func:` as plain text. Check if `nitpicky = True` is in conf.py; if not, manually verify cross-refs against code.
- **L2**: Sphinx does NOT warn on duplicate `:label:` across files. Run `grep :label: docs/**/*.rst | sort | uniq -d` before building. Always page-prefix labels (e.g., `moc-keff-update`, `cp-keff-update`).
- **L4**: `:func:`_private_name`` won't resolve unless `:private-members:` is on the automodule. Use backtick code literals for private functions.
- **L9**: Before writing any `:func:` or `:class:` directive, verify the target module appears in `docs/api/`. Missing autodoc = dead cross-refs.

## Content quality rules

- **L3**: Check for dead citations: `[Name20XX]_` in body vs bibliography. Common and easy to miss.
- **L5**: Never hardcode test counts in RST. Run `.venv/bin/python -m pytest --collect-only -q` to get current count.
- **L6**: Quality bottleneck is `VerificationCase.latex` in derivation modules, not the RST or generator. To improve docs, improve derivation source.
- **L7**: For simulation modules (TH, RK), document the **discrete FD stencil first** — the continuous PDE is in every textbook. The code's value is the specific discretization.
- **L8**: Three-layer tracing: `derivations/ -> generate_rst.py -> RST`. Quality issues can originate at any layer. Always trace to source.

## Page quality baseline (last audit)

| Page | Score | Key gap |
|------|-------|---------|
| method_of_characteristics.rst | 8.3 | Core equations hand-written (need derivation scripts) |
| monte_carlo.rst | 4.0 | 35 non-resolving autodoc directives |
| collision_probability.rst | 3.7 | Missing GS convergence comparison table |
| discrete_ordinates.rst | 3.5 | 5 equations without derivation scripts |
| cross_section_data.rst | 2.5 | No derivation script, broken cross-refs |
| verification.rst | 2.3 | Thin content at every layer |
| thermal_hydraulics.rst | 2.2 | No FD discretization, no derivation scripts |
| reactor_kinetics.rst | 1.7 | Same gaps as TH |
