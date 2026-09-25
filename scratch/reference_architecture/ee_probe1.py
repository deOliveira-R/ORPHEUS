import ast, pathlib, collections
root = pathlib.Path("orpheus")
cache_dec = collections.Counter(); cache_sites = []
mpdps = []; lambdify = []; sympyfunc = []; environ = []; numba = []; modcall = []
first = ("orpheus",)
for p in sorted(root.rglob("*.py")):
    if "__pycache__" in p.parts: continue
    src = p.read_text(); tree = ast.parse(src)
    pkg = p.parts[1] if len(p.parts) > 2 else "(root)"
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for d in node.decorator_list:
                s = ast.unparse(d)
                if any(k in s for k in ("lru_cache", "functools.cache", "cache", "cached_property", "cacheit")) and "cached_property" not in s:
                    cache_dec[pkg] += 1; cache_sites.append(f"{p}:{node.lineno} {node.name} @{s}")
        if isinstance(node, ast.Assign):
            for t in node.targets:
                s = ast.unparse(t)
                if s.endswith(".dps") or s.endswith(".prec"):
                    mpdps.append(f"{p}:{node.lineno} {ast.unparse(node)[:80]}")
        if isinstance(node, ast.Call):
            s = ast.unparse(node.func)
            if s.endswith("lambdify"): lambdify.append(f"{p}:{node.lineno}")
            if s in ("os.environ.get", "os.getenv") or s.endswith("environ.get"): environ.append(f"{p}:{node.lineno} {ast.unparse(node)[:90]}")
        if isinstance(node, ast.Subscript) and ast.unparse(node.value) == "os.environ":
            environ.append(f"{p}:{node.lineno} {ast.unparse(node)[:90]}")
        if isinstance(node, ast.ClassDef):
            for b in node.bases:
                bs = ast.unparse(b)
                if bs in ("Function", "sympy.Function", "sp.Function"): sympyfunc.append(f"{p}:{node.lineno} {node.name}({bs})")
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            m = getattr(node, "module", None) or ""
            names = [a.name for a in node.names]
            if "numba" in m or any("numba" in n for n in names) or "ctypes" in m or any(n in ("ctypes","cffi") for n in names):
                numba.append(f"{p}:{node.lineno}")
    # module-level assignment whose value calls a function (runs at import)
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            calls = [ast.unparse(c.func) for c in ast.walk(node.value) if isinstance(c, ast.Call)]
            calls = [c for c in calls if not c.split(".")[0] in ("np","numpy","sp","sympy","mp","mpmath","math","frozenset","tuple","dict","list","set","dataclass","field","TypeVar","NewType","Path","pathlib","re","logging","object","float","int","str","len","range","Literal","Union","Optional","cast","ParamSpec","namedtuple","collections","typing","enum","functools","os","sys","Enum","getLogger")]
            if calls: modcall.append(f"{p}:{node.lineno} {calls[:3]}")
print("CACHE decorators by pkg", dict(cache_dec)); print("\n".join(cache_sites))
print("\nMP dps/prec assignments", len(mpdps)); print("\n".join(mpdps[:40]))
print("\nlambdify calls", len(lambdify), collections.Counter(s.split("/")[1] for s in lambdify))
print("\nsympy Function subclasses", sympyfunc)
print("\nenviron reads", len(environ)); print("\n".join(environ))
print("\nnumba/ctypes", numba)
print("\nmodule-level assigns calling non-stdlib callables", len(modcall), collections.Counter(s.split("/")[1] for s in modcall))
