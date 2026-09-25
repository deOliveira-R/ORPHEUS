import json, collections, re, sys
c = json.load(open(sys.argv[1]))
P = "tests/gates/cp/test_peierls_rank_n_protocol.py"
# the 12 deselected subprocess cases: static TOUCH (they call _run_f4_subprocess)
for fn, ids in [("test_f4_is_sign_stable_at_its_reference_quadrature", 6), ("test_f4_rich_vs_rich_panels_matches_pinned_baseline", 6)]:
    for i in range(ids): c[f"{P}::{fn}[subproc{i}]"] = "TOUCH"
W = {k: (v == "TOUCH") for k, v in c.items()}
files = collections.defaultdict(list)
for k, w in W.items(): files[k.split("::")[0]].append((k, w))
out = {}
tot_w = tot_k = 0
for f, rows in sorted(files.items()):
    nw = sum(w for _, w in rows); nk = len(rows) - nw; tot_w += nw; tot_k += nk
    # group by (class, function) sans params
    grp = collections.defaultdict(list)
    for k, w in rows:
        parts = k.split("::")[1:]
        parts[-1] = re.sub(r"\[.*\]$", "", parts[-1])
        grp[tuple(parts)].append(w)
    mixed_fn = [g for g, ws in grp.items() if 0 < sum(ws) < len(ws)]
    # classes uniform?
    cls = collections.defaultdict(list)
    for g, ws in grp.items():
        if len(g) == 2: cls[g[0]].extend(ws)
    if nk == 0: placement = "FILE"
    else: placement = "PER-TEST"
    out[f] = dict(withdrawn=nw, kept=nk, placement=placement,
        mixed_functions=["::".join(g) for g in mixed_fn],
        withdrawn_units=[], kept_units=[])
    # units: a whole class if uniform-withdrawn, else functions
    done = set()
    for cname, ws in cls.items():
        if all(ws):
            out[f]["withdrawn_units"].append(f"class {cname} ({len(ws)})"); done.add(cname)
        elif not any(ws):
            out[f]["kept_units"].append(f"class {cname} ({len(ws)})"); done.add(cname)
    for g, ws in grp.items():
        if len(g) == 2 and g[0] in done: continue
        name = "::".join(g)
        if all(ws): out[f]["withdrawn_units"].append(f"{name} ({len(ws)})")
        elif not any(ws): out[f]["kept_units"].append(f"{name} ({len(ws)})")
        else: out[f]["withdrawn_units"].append(f"MIXED {name} ({sum(ws)}/{len(ws)})")
print("withdrawn", tot_w, "kept", tot_k, "files", len(files), "FILE-level", sum(1 for v in out.values() if v['placement']=='FILE'))
json.dump(out, open(sys.argv[2], "w"), indent=1)
for f, v in out.items():
    print(f"\n## {f.replace('tests/gates/','')}  W={v['withdrawn']} K={v['kept']} {v['placement']}")
    if v["placement"] != "FILE":
        print("  W:", "; ".join(v["withdrawn_units"]))
        print("  K:", "; ".join(v["kept_units"]))
