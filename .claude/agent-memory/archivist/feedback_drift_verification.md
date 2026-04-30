---
name: Verify drift before fixing
description: Always diff the user's quoted "stale text" against the committed file before editing — drift may already be fixed
type: feedback
---

When a task description quotes "current stale text" and asks for a
fix, **always verify the quote against the committed file first**.
The user's snapshot may pre-date a fix that already landed.

**Why:** In the Phase G.5 docs reconcile (2026-04-23), drift item 1
quoted "Slab continues to use the native..." but `peierls_unified.rst`
already said "Slab's shipped peierls_slab_2eg_2rg continuous reference
now routes through the unified...path by default" (committed in
`3b0b2c9`/`529cdbe`). The stale quote was the pre-3b0b2c9 state from
commit `aa6ebf0`. Editing blindly would have re-introduced the stale
wording.

**How to apply:**

1. Read the file at the line range first.
2. If the current text matches the user's "before" quote, proceed.
3. If the current text already matches the desired "after" state,
   `git show <prior-commit>:<path>` to confirm the historical drift,
   then **report "already fixed"** rather than editing.
4. The capability-matrix descriptor (drift item 2) was the genuine
   one — `peierls_cases.py::capability_rows` is the single source
   of truth; `_peierls_capability_matrix.inc.rst` regenerates from
   it on `sphinx-build`. Edit the Python, never the .inc.rst.
