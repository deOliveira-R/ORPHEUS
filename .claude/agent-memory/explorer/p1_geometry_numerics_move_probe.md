---
name: p1-geometry-numerics-move-probe
description: P1 layer-DAG (geometry<->numerics) inject-and-run probe recipe and the hidden members a module move hits (path-segment and in-string imports; the digest and memo-source hazards are retired, see body)
metadata:
  type: project
---

State at main f672a2d5 / a846500f (2026-10-22). Re-verify before repeating any count.

**Probe recipe (worked, ~5 min):** rsync `orpheus tests pyproject.toml` into the scratchpad (rm -rf is denied: mkdir+rsync only), apply an AST import rewriter (resolve relative imports against the ORIGINAL module name, split multi-name `from X import a, b` by destination), then a TEXT pass for the same dotted paths (imports inside subprocess strings and dotted string literals are invisible to AST), then run cold `python -c "import <entry>"` with PYTHONPATH=<copy> from /tmp, the AST package census on the copy, and the targeted gates with cwd=<copy>. Scripts: session scratchpad `transform.py`, `simulate.py`, `plug/modmove.py` (a pytest plugin that rewrites `__module__` of classes to simulate a move for digest pins).

**Hidden members a move hits (none is an import edge):**
- PAST (until layer DAG ruling 13, 2026-10-22): `_schema_tag` encoded `cls.__module__`, so a move re-keyed digests. NOW `orpheus.values.content._type_name` names an orpheus class by qualname only (a gate keeps qualnames unique), a foreign type by `module:qualname`: a move changes no digest.
- PAST (until `b1977d2c`): `traced_memo._MEMO_SOURCES` named `content.py` as a sibling filename. NOW it names it through the imported module (`Path(_content.__file__)`); `content.py` and `scalars.py` live in `orpheus/values/`.
- Path-segment literals in tests (`_NUMERICS / "symmetry.py"`, `"orpheus"/"numerics"/"content.py"`): grep the last segment.
- A move-only cut is not acyclic: the geometry door (`structured_geometry.parse_boundary_law`'s isinstance on the law base) keeps the SCC one package up. A package `__init__` re-exporting a name that moved UP creates a new SCC.

Related: [[geometry-numerics-package-cycle]] (the module-grain finding this extends).
