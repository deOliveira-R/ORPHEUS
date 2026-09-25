import ast, pathlib, sys
ROOT = pathlib.Path("orpheus")
def modpath(name):
    p = pathlib.Path(*name.split("."))
    for c in (p.with_suffix(".py"), p / "__init__.py"):
        if c.is_file(): return c
    return None
def modname(path):
    parts = list(path.with_suffix("").parts)
    if parts[-1] == "__init__": parts = parts[:-1]
    return ".".join(parts)
def imports(path):
    tree = ast.parse(path.read_text())
    pkg = modname(path) if path.name == "__init__.py" else modname(path).rsplit(".", 1)[0]
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names: out.add(a.name)
        elif isinstance(n, ast.ImportFrom):
            if n.level:
                base = pkg.split(".")
                base = base[: len(base) - (n.level - 1)]
                mod = ".".join(base + ([n.module] if n.module else []))
            else:
                mod = n.module or ""
            out.add(mod)
            for a in n.names: out.add(f"{mod}.{a.name}")
    res = set()
    for m in out:
        if not m.startswith("orpheus"): continue
        # add module and all parent packages (their __init__ executes)
        parts = m.split(".")
        for i in range(1, len(parts) + 1):
            p = modpath(".".join(parts[:i]))
            if p: res.add(p)
    return res
def closure(start):
    seen, todo = set(), [start]
    while todo:
        p = todo.pop()
        if p in seen: continue
        seen.add(p); todo.extend(imports(p) - seen)
    return seen
for m in sys.argv[1:]:
    p = modpath(m); c = closure(p)
    pk = {}
    for f in c: pk[f.parts[1] if len(f.parts) > 2 else "(root)"] = pk.get(f.parts[1] if len(f.parts) > 2 else "(root)", 0) + 1
    print(m, len(c), "files,", sum(len(f.read_text().splitlines()) for f in c), "lines", dict(sorted(pk.items())))
