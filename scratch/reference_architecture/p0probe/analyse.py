import json, collections, sys
d = json.load(open(sys.argv[1]))
print("rebinds", d["rebinds"])
T = d["tests"]
cls = {}
for nid, v in T.items():
    if nid.startswith("<"): print("NON-TEST hits:", nid, v["hits"]); continue
    touch = bool(v["hits"]) or v.get("wr")
    oc = v.get("outcome", [])
    failed_other = (not touch) and any(o.endswith(":failed") for o in oc)
    cls[nid] = "TOUCH" if touch else ("KEPT-FAIL" if failed_other else ("KEPT-SKIP" if any(o.endswith(":skipped") for o in oc) and not any(o=="call:passed" for o in oc) else "KEPT"))
c = collections.Counter(cls.values()); print(c, len(cls))
byf = collections.defaultdict(collections.Counter)
for nid, k in cls.items(): byf[nid.split("::")[0]][k] += 1
for f in sorted(byf): print(f"{dict(byf[f])}  {f.replace('tests/gates/','')}")
json.dump(cls, open(sys.argv[2], "w"), indent=0)
for nid,k in sorted(cls.items()):
    if k in ("KEPT-FAIL","KEPT-SKIP"): print(k, nid, T[nid]["outcome"])
