import sys, functools, dataclasses
M = sys.monitoring; TID = 3
M.use_tool_id(TID, "refcache-probe")
rec = None
def on_start(code, off):
    if rec is not None: rec.add((code.co_filename.rsplit("/",1)[-1], code.co_qualname))
    return None               # no DISABLE
M.register_callback(TID, M.events.PY_START, on_start)
M.set_events(TID, M.events.PY_START)
@functools.lru_cache
def kernel(n): return helper(n) * 2
def helper(n): return n + 1
@dataclasses.dataclass(frozen=True)
class Spec: a: float
def gen_A(): return kernel(3) + Spec(1.0).a
def gen_B(): return kernel(3)
rec = set(); gen_A(); tA = rec
rec = set(); gen_B(); tB = rec; rec = None
print("trace A:", sorted(tA)); print("trace B:", sorted(tB))
print("B saw kernel/helper?", any("kernel" in q for _, q in tB), any("helper" in q for _, q in tB))
# DISABLE variant
M.set_events(TID, 0); M.free_tool_id(TID)
M.use_tool_id(TID, "p2"); seen=set()
def on2(code, off):
    seen.add(code.co_qualname); return M.DISABLE
M.register_callback(TID, M.events.PY_START, on2); M.set_events(TID, M.events.PY_START)
def f(): return 1
f(); s1=set(seen); seen.clear(); f(); s2=set(seen)
print("DISABLE: first", "f" in s1, "second", "f" in s2)
M.restart_events(); seen.clear(); f(); print("after restart_events:", "f" in seen)
