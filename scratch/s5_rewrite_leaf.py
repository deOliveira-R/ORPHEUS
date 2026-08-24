"""S5.2/S5.3 leaf-sugar rewriter — span-exact, family-keyed.

  X.from_mesh(V, M)     -> X(values=V, space=M.<family attr>)
  X.zeros_on(M)         -> X.zeros(M.<family attr>)
  X.from_ndarray(A, M)  -> X(values=A, space=M.<family attr>)

Sites passing spatial_moments=, receivers not in the map (cls, aliases),
or extra kwargs are WARNED and left for hand migration.
Usage: s5_rewrite_leaf.py [--apply] [root ...]  (default: orpheus)
"""
import ast, pathlib, sys

DRY = "--apply" not in sys.argv
roots = [a for a in sys.argv[1:] if not a.startswith("--")] or ["orpheus"]

SPACE = {}
for cls in ("AngularFlux", "AngularSourceSink", "AngularResidual"):
    SPACE[cls] = "angular_bulk_space"
for cls in ("ScalarFlux", "ScalarSourceSink", "ScalarResidual", "CrossSectionField"):
    SPACE[cls] = "bulk_space"
for cls in ("AngularBoundaryFlux", "AngularBoundarySourceSink", "AngularBoundaryResidual"):
    SPACE[cls] = "angular_trace"
for cls in ("ScalarBoundaryFlux", "ScalarBoundarySourceSink"):
    SPACE[cls] = "scalar_trace"
for cls in ("RadialCharacteristicInteriorFlux", "RadialCharacteristicInteriorSourceSink", "RadialCharacteristicInteriorResidual"):
    SPACE[cls] = "radial_characteristic_interior_space"
for cls in ("RadialCharacteristicBoundaryFlux", "RadialCharacteristicBoundarySourceSink", "RadialCharacteristicBoundaryResidual"):
    SPACE[cls] = "radial_characteristic_boundary_space"

SKIP = ("orpheus/transport/fields/_bases.py", "scratch")

def offs(src):
    table, n = [0], 0
    for line in src.splitlines(keepends=True):
        n += len(line); table.append(n)
    return table

edits, warns = 0, 0
for root in roots:
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
        seg = lambda n_: src[pos(n_.lineno, n_.col_offset):pos(n_.end_lineno, n_.end_col_offset)]
        spans = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            if not isinstance(f, ast.Attribute):
                continue
            meth = f.attr
            if meth not in ("from_mesh", "zeros_on", "from_ndarray"):
                continue
            recv = f.value.id if isinstance(f.value, ast.Name) else ast.unparse(f.value)
            if recv not in SPACE:
                print(f"HAND unmapped-recv {rp}:{node.lineno} {recv}.{meth}"); warns += 1
                continue
            if any(kw.arg == "spatial_moments" for kw in node.keywords):
                print(f"HAND spatial_moments {rp}:{node.lineno} {recv}.{meth}"); warns += 1
                continue
            if node.keywords:
                print(f"HAND extra-kwargs {rp}:{node.lineno} {recv}.{meth}"); warns += 1
                continue
            attr = SPACE[recv]
            cs, ce = pos(node.lineno, node.col_offset), pos(node.end_lineno, node.end_col_offset)
            if meth == "zeros_on":
                if len(node.args) != 1:
                    print(f"HAND arity {rp}:{node.lineno}"); warns += 1; continue
                M = seg(node.args[0])
                new = f"{recv}.zeros({M}.{attr})"
            else:
                if len(node.args) != 2:
                    print(f"HAND arity {rp}:{node.lineno}"); warns += 1; continue
                V, M = seg(node.args[0]), seg(node.args[1])
                new = f"{recv}(values={V}, space={M}.{attr})"
            spans.append((cs, ce, new, node.lineno))
        if not spans:
            continue
        out = src
        for s, e, new, ln in sorted(spans, reverse=True):
            print(f"{'DRY ' if DRY else ''}{rp}:{ln} -> {new[:90]}")
            out = out[:s] + new + out[e:]
            edits += 1
        if not DRY:
            p.write_text(out)
print(f"\n{edits} edits, {warns} hand-sites ({'dry run' if DRY else 'APPLIED'})")
