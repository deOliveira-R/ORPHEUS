"""S5.1 call-site rewriter — composite tier (mesh= -> space=), span-exact.

AST offsets -> character spans, applied in reverse order per file:
  1. {FullField,TimedFullField}.zeros(..., mesh=X) -> space=X.full_field_space
  2. RadialCharacteristicField.from_mesh(X)        -> .flux_zeros(X.radial_characteristic_field_space)
  3. RadialCharacteristicField.source_zeros_on(X)  -> .source_zeros(X.radial_characteristic_field_space)
"""
import ast, pathlib, re, sys

DRY = "--apply" not in sys.argv
SKIP = {"orpheus/transport/radial_characteristic_field.py",
        "orpheus/transport/full_field.py",
        "orpheus/transport/timed_full_field.py",
        "scratch"}

def offs(src):
    table, n = [0], 0
    for line in src.splitlines(keepends=True):
        n += len(line); table.append(n)
    return table

edits, warns = 0, 0
for root in ("orpheus", "tests"):
    for p in sorted(pathlib.Path(root).rglob("*.py")):
        rp = str(p)
        if any(rp.startswith(s) for s in SKIP):
            continue
        src = p.read_text()
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        table = offs(src)
        pos = lambda ln, col: table[ln - 1] + col
        spans = []  # (start, end, replacement, lineno)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            if not isinstance(f, ast.Attribute) or not isinstance(f.value, ast.Name):
                continue
            recv, meth = f.value.id, f.attr
            if meth == "zeros" and recv.endswith("FullField"):
                for kw in node.keywords:
                    if kw.arg == "mesh":
                        vs = pos(kw.value.lineno, kw.value.col_offset)
                        ve = pos(kw.value.end_lineno, kw.value.end_col_offset)
                        m = re.search(r"mesh\s*=\s*$", src[:vs])
                        if not m:
                            print(f"WARN kw-anchor {rp}:{node.lineno}"); warns += 1; continue
                        V = src[vs:ve]
                        spans.append((m.start(), ve, f"space={V}.full_field_space", node.lineno))
            elif recv == "RadialCharacteristicField" and meth in ("from_mesh", "source_zeros_on"):
                if len(node.args) != 1 or node.keywords:
                    print(f"WARN odd-args {rp}:{node.lineno}"); warns += 1; continue
                a = node.args[0]
                A = src[pos(a.lineno, a.col_offset):pos(a.end_lineno, a.end_col_offset)]
                cs = pos(node.lineno, node.col_offset)
                ce = pos(node.end_lineno, node.end_col_offset)
                new_meth = "flux_zeros" if meth == "from_mesh" else "source_zeros"
                spans.append((cs, ce,
                    f"RadialCharacteristicField.{new_meth}({A}.radial_characteristic_field_space)",
                    node.lineno))
        if not spans:
            continue
        out = src
        for s, e, new, ln in sorted(spans, reverse=True):
            print(f"{'DRY ' if DRY else ''}{rp}:{ln} [{src[s:e][:50]!r} -> {new[:70]!r}]")
            out = out[:s] + new + out[e:]
            edits += 1
        if not DRY:
            p.write_text(out)
print(f"\n{edits} edits, {warns} warnings ({'dry run' if DRY else 'APPLIED'})")
