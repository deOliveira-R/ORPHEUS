"""Dynamic attribution: which collected tests REACH the withdrawn Peierls
Nystrom solver half. Every withdrawn symbol is rebound, in every
sys.modules binding, to a recorder that raises WithdrawnReached."""
import functools, importlib, json, os, sys
import pytest

class WithdrawnReached(BaseException):
    pass

MODS = {
 "orpheus.derivations.continuous.peierls_nystrom.geometry": [
   "K_vol_element_adaptive","build_volume_kernel_adaptive","build_volume_kernel",

   "build_white_bc_correction","build_closure_operator","_build_closure_operator_rank2_white",
   "_build_closure_operator_rank_n_white","build_white_bc_correction_rank_n",
   "_build_slab_per_face_specular_PG","_build_sphere_specular_mode_PG","_build_cylinder_specular_mode_PG",
   "_build_white_rank1_mark_op","_build_white_f4_op","_build_white_hebert_op","_build_specular_op",
   "_build_specular_multibounce_op","_build_full_K_per_group","solve_peierls_mg","solve_peierls_1g"],
 "orpheus.derivations.continuous.peierls_nystrom.slab": [
   "_build_kernel_matrix","_build_system_matrices","solve_peierls_eigenvalue","_build_peierls_slab_case"],
 "orpheus.derivations.continuous.peierls_nystrom.cylinder": [
   "_build_peierls_cylinder_case","_build_peierls_cylinder_hollow_f4_case"],
 "orpheus.derivations.continuous.peierls_nystrom.sphere": [
   "_build_peierls_sphere_case","_build_peierls_sphere_hollow_f4_case"],
 "orpheus.derivations.continuous.peierls_nystrom.cases": [
   "build_two_surface_case","_build_peierls_slab_case_via_unified","build_one_surface_compact_case",
   "_build","_class_a_cases","continuous_cases"],
}
CLASS_INIT = {"orpheus.derivations.continuous.peierls_nystrom.geometry": ["BoundaryClosureOperator"]}

HITS = {}          # nodeid -> set(symbol)
CURRENT = ["<collection>"]
ORIG = {}          # id(orig) -> (wrapper, name)
REBINDS = [0]

def _wrap(name, fn):
    @functools.wraps(fn)
    def w(*a, **k):
        HITS.setdefault(CURRENT[0], set()).add(name)
        raise WithdrawnReached(name)
    return w

def _sweep():
    n = 0
    for mname, mod in list(sys.modules.items()):
        d = getattr(mod, "__dict__", None)
        if not isinstance(d, dict):
            continue
        for k, v in list(d.items()):
            try:
                hit = ORIG.get(id(v))
            except Exception:
                continue
            if hit is not None and hit[2] is v:
                setattr(mod, k, hit[0]); n += 1
    REBINDS[0] += n
    return n

def _install():
    for mname, names in MODS.items():
        mod = sys.modules.get(mname) or importlib.import_module(mname)
        for nm in names:
            fn = getattr(mod, nm)
            if getattr(fn, "__withdrawn_probe__", False):
                continue
            w = _wrap(f"{mname.rsplit('.',1)[1]}.{nm}", fn); w.__withdrawn_probe__ = True
            ORIG[id(fn)] = (w, nm, fn)
    return _sweep()

def pytest_runtest_setup(item):
    n = _install()
    if n: RELOADS.append((item.nodeid, n))

RELOADS = []

def pytest_configure(config):
    _install_first = True
    for mname, names in MODS.items():
        mod = importlib.import_module(mname)
        for nm in names:
            fn = getattr(mod, nm)
            w = _wrap(f"{mname.rsplit('.',1)[1]}.{nm}", fn); w.__withdrawn_probe__ = True
            ORIG[id(fn)] = (w, nm, fn)
    for mname, names in CLASS_INIT.items():
        mod = importlib.import_module(mname)
        for nm in names:
            cls = getattr(mod, nm); init = cls.__init__
            cls.__init__ = _wrap(f"{nm}.__init__", init); REBINDS[0] += 1
    n = _sweep()
    if n < sum(len(v) for v in MODS.values()):
        raise pytest.UsageError(f"probe did not install: {n} rebinds")

def pytest_collection_finish(session):
    _sweep()

@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    CURRENT[0] = item.nodeid
    yield
    CURRENT[0] = "<between>"

OUT = {}
def pytest_runtest_logreport(report):
    rec = OUT.setdefault(report.nodeid, {"outcome": [], "wr": False})
    rec["outcome"].append(f"{report.when}:{report.outcome}")
    if report.failed and report.longrepr is not None and "WithdrawnReached" in str(report.longrepr):
        rec["wr"] = True

def pytest_unconfigure(config):
    path = os.environ.get("P0PROBE_OUT", "probe_out.json")
    res = {nid: {"hits": sorted(HITS.get(nid, ())), **OUT.get(nid, {})} for nid in set(OUT) | set(HITS)}
    json.dump({"rebinds": REBINDS[0], "reinstalls": RELOADS, "tests": res}, open(path, "w"), indent=0)
    print(f"\n[withdrawn_probe] rebinds={REBINDS[0]} tests_recorded={len(res)} "
          f"touch={sum(1 for v in res.values() if v['hits'] or v.get('wr'))}", file=sys.stderr)
