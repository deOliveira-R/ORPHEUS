"""Aggregate the fanout's inventory and CHECK it against the graph.

A claimed implementer that does not resolve to a real node would become a
dangling `confidence=1.0` edge — worse than the guess it replaces. This
resolves every one before a single directive is written.
"""
import re, sys, json
from collections import Counter, defaultdict
from pathlib import Path
from sphinxcontrib.nexus.export import load_sqlite

OUT = Path("/private/tmp/claude-501/-Users-rodrigo-git-nuclear-ORPHEUS/ffaa72e7-def2-4759-91e6-ff59adf48822/scratchpad/impl_inventory")
g = load_sqlite(Path("/Users/rodrigo/git/nuclear/ORPHEUS/.nexus/graph.db")).nxgraph

# dotted name -> node id, for every code node the ontology allows
# `[edge.implements]` declares domain = function|method|class. Resolve
# against EVERY node type so an ontology violation is distinguishable
# from a symbol that does not exist — nexus#86: the directive path does
# NOT consult the ontology, so a forbidden target ships silently.
ALLOWED = {"function", "method", "class"}
by_name = defaultdict(list)
for n, a in g.nodes(data=True):
    if a.get("name"):
        by_name[a["name"]].append(n)
expected, meta = set(), {}
for f in sorted(OUT.glob("slice_*.json")):
    for rs in json.loads(f.read_text())["pages"].values():
        for r in rs:
            expected.add(r["label"]); meta[r["label"]] = r

rows, seen = [], set()
for f in sorted(OUT.glob("out_*.md")):
    cur = {}
    for line in f.read_text().splitlines():
        if re.match(r"^#{2,4} ", line):
            if cur.get("label"): rows.append(cur)
            lab = line.lstrip("#").strip().strip("`")
            cur = {"label": lab, "slice": f.stem, "impl": []} if lab in expected else {}
        elif ":" in line and cur:
            k, _, v = line.partition(":")
            k = k.strip().lower(); v = v.strip()
            if k == "verdict": cur["verdict"] = v.split()[0] if v else "?"
            elif k == "implementers":
                cur["impl"] = [x.strip() for x in re.split(r"[,\s]+", v.split("#")[0]) if "." in x]
            elif k in ("confidence", "evidence", "note"): cur[k] = v
    if cur.get("label"): rows.append(cur)

for r in rows: seen.add(r["label"])
print(f"rows parsed: {len(rows)}   distinct labels: {len(seen)}")
print(f"expected:    {len(expected)}   MISSING: {len(expected - seen)}   EXTRA: {len(seen - expected)}")
if expected - seen:
    print("  missing sample:", sorted(expected - seen)[:6])
if seen - expected:
    print("  ⚠ EXTRA (not in any slice — hallucinated or renamed):", sorted(seen - expected)[:6])

print("\nverdicts:", dict(Counter(r.get("verdict", "?") for r in rows)))
print("confidence:", dict(Counter(r.get("confidence", "?") for r in rows)))

print("\n=== implementer resolution (the check that matters) ===")
res = Counter(); unresolved = []
n_impl = 0
for r in rows:
    if not r.get("verdict", "").startswith("DECLARABLE"): continue
    if not r["impl"]:
        res["DECLARABLE with NO implementer listed"] += 1; continue
    for sym in r["impl"]:
        n_impl += 1
        hits = by_name.get(sym, [])
        if not hits:
            res["DOES NOT EXIST"] += 1; unresolved.append((r["label"], sym, "no such node"))
            continue
        types = {g.nodes[h].get("type") for h in hits}
        if types & ALLOWED:
            res["ok" if len(hits) == 1 else "ok (ambiguous, >1 node)"] += 1
        else:
            res["WRONG NODE TYPE (ontology forbids)"] += 1
            unresolved.append((r["label"], sym, f"is a {sorted(types)[0]}"))
print(f"  {n_impl} claimed implementers over "
      f"{sum(1 for r in rows if r.get('verdict','').startswith('DECLARABLE'))} DECLARABLE rows")
for k, v in res.most_common(): print(f"    {k:34s} {v:4d}")
if unresolved:
    print(f"\n  ⛔ UNRESOLVED — these would become dangling edges ({len(unresolved)}):")
    for lab, sym, why in unresolved[:20]:
        print(f"      {lab:38s} -> {sym[:52]:52s} [{why}]")

# the "declare EVERY implementer" contract: how many list >1?
multi = Counter(len(r["impl"]) for r in rows if r.get("verdict","").startswith("DECLARABLE"))
print(f"\n  implementers per DECLARABLE row: {dict(sorted(multi.items()))}")


# ── the ASYMMETRIC half ────────────────────────────────────────────────
# A wrong DECLARABLE is caught mechanically above (the symbol resolves or
# it does not). A wrong NOTHING cannot be: it would suppress the guesses
# AND record that nothing implements the equation, hiding a real coverage
# gap behind a confident declaration. It fails in the FLATTERING
# direction, so it needs its own screen.
print("\n=== NOTHING verdicts that the evidence CONTRADICTS ===")
susp = []
for r in rows:
    v = r.get("verdict", "")
    if not v.startswith("NOTHING"): continue
    m = meta.get(r["label"], {})
    ev, guesses, n = (m.get("evidence_candidates") or [],
                      m.get("guessed_implementers") or [], m.get("n_claims", 0))
    why = []
    if ev: why.append(f"{len(ev)} test-EXECUTED candidates")
    if n >= 10: why.append(f"{n} claiming tests")
    if guesses: why.append(f"{len(guesses)} name-match guesses")
    if why: susp.append((r["label"], v, "; ".join(why), (ev or guesses)[:2]))
susp.sort(key=lambda x: -len(x[2]))
print(f"  {len(susp)} of {sum(1 for r in rows if r.get('verdict','').startswith('NOTHING'))}"
      f" NOTHING rows have contrary evidence:")
for lab, v, why, cand in susp[:20]:
    print(f"    {lab:38s} {v:26s} {why}")
    if cand: print(f"        candidates: {', '.join(c[:46] for c in cand)}")
