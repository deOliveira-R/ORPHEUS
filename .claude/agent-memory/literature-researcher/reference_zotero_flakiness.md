---
name: Zotero MCP server can be unreachable mid-session
description: Observed symptom set for Zotero MCP outage — list_libraries succeeds but search/get_recent return "No items found" or connection refused. Proceed with Tier 2 + disclosure.
type: reference
---

Zotero MCP server failure mode observed 2026-04-17 (Phase-4.2
cylinder research session):

- `zotero_list_libraries` → succeeds, reports full item count (15572).
- `zotero_search_items(query=anything)` → `"No items found matching
  query: '…'"` — this is a *silent* failure: the tool claims the
  query ran but returned zero hits. A 1-word query like `neutron`
  cannot legitimately return 0 hits in a nuclear-engineering
  library, so treat "0 hits on a common term" as a signal the
  server is broken, not that the library lacks the paper.
- `zotero_get_recent`, `zotero_get_collections`,
  `zotero_advanced_search`, `zotero_search_by_citation_key` → all
  return `Error … [Errno 111] Connection refused`.

**Why:** Likely local Zotero desktop daemon not running or its
local HTTP endpoint (port 23119) not forwarding into the
devcontainer reliably. See `reference_zotero_mcp.md` in the
top-level memory for the setup details.

**How to apply:**
1. First attempt: run 2–3 parallel `zotero_search_items` calls
   with single-token queries (e.g. `Sanchez`, `neutron`).
2. If ALL return zero hits AND another endpoint (get_recent,
   get_collections) returns connection-refused, declare the
   server unavailable for this session.
3. Fall back to Tier 2 (OpenAlex, CrossRef, Semantic Scholar) and
   **flag the Zotero unavailability explicitly at the top of
   the deliverable** so the user knows user-annotations were not
   consulted.
4. Do not waste >1 minute retrying — the failure persists for
   many minutes.

**Sharpened 2026-09-22 — the server-down fallback that still reads the
library:** when port 23119 refuses AND no `mcp__zotero__*` tool is exposed,
the library is still searchable READ-ONLY. The data dir is in
`~/Library/Application Support/Zotero/Profiles/*/prefs.js`
(`extensions.zotero.dataDir` = `/Users/rodrigo/Documents/Biblioteca/#Zotero`).
Copy `zotero.sqlite` to the session scratchpad and open it with
`sqlite3.connect("file:<copy>?mode=ro", uri=True)`; titles are
`items ⋈ itemData ⋈ fields(fieldName='title') ⋈ itemDataValues`, creators
`items ⋈ itemCreators ⋈ creators`. Never open the live file (Zotero locks it)
and never write. Annotations are in the same DB but were not needed that
time.
