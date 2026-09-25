"""Census of tests/gates consumers of orpheus.derivations (Q1)."""
import ast, re, subprocess, json, sys, collections
ROOT = "/Users/rodrigo/git/nuclear/ORPHEUS"
files = subprocess.run(["git","ls-files","tests/gates"],cwd=ROOT,capture_output=True,text=True).stdout.split()
files = [f for f in files if f.endswith(".py")]

PROD_MODS = re.compile(r"^orpheus\.(sn|cp|moc|mc|diffusion|homogeneous)(\.|$)|^orpheus\.transport\.mesh")
PROD_CALL = re.compile(r"\b(SNProblem(\.from_axes|\.from_material_mesh)?|DiffusionMesh(\.from_material_mesh)?|HomogeneousProblem|CPMesh|MaterialMesh|solve_sn|solve_cp|solve_moc|solve_mc|solve_diffusion|solve_homogeneous|solve_\w+_fixed_source|MoCMesh|MCProblem|MoCProblem)\s*\(")
# reference-data -> production feeds
REF_FEED = re.compile(r"\.(problem\.(materials|geometry_params|geometry_type|boundary_conditions|n_groups)|geom_params|geometry_params|build_mesh|build_materials|build_problem|build_sn_problem|to_geometry|material_mesh)\b|\bbuild_(mesh|materials)\(")
CASE_MATERIALS = re.compile(r"\b(case|ref|reference|cr|c|mms|spec|rc|registry_case|vc|refsol|sol_ref)\.materials\b")
HAND = re.compile(r"\b(Mixture|make_mixture|get_mixture|get_materials|Mesh1D|Mesh2D|mesh1d_from_zones|RegionMesh|StructuredGeometry|Region|from_zones|homogeneous_1d|placeholder_materials|material_xs_from_raw)\s*\(|\bedges\s*=\s*np\.")
def imports(tree):
    mods=[]
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom) and n.module:
            mods.append((n.module,[a.name for a in n.names]))
        elif isinstance(n, ast.Import):
            for a in n.names: mods.append((a.name,[]))
    return mods
rows=[]
for f in files:
    src=open(f"{ROOT}/{f}").read()
    try: tree=ast.parse(src)
    except SyntaxError: continue
    mods=imports(tree)
    der=[(m,n) for m,n in mods if m.startswith("orpheus.derivations") or (m=="orpheus" and "derivations" in n)]
    rv = "reference_values" in src
    if not der and not rv: continue
    dermods=sorted({m for m,_ in der})
    xs_only = all(m=="orpheus.derivations.common.xs_library" for m in dermods) and not rv
    # also 'from orpheus.derivations import get_mixture...' style
    if dermods==["orpheus.derivations"]:
        names={x for m,n in der for x in n}
        if names <= {"get_mixture","get_materials","get_xs","make_mixture","validate_all"}: xs_only=True
    prod_imp = any(PROD_MODS.match(m) for m,_ in mods)
    prod_call = len(PROD_CALL.findall(src))
    feed = len(REF_FEED.findall(src)) + len(CASE_MATERIALS.findall(src))
    hand = len(HAND.findall(src))
    has_prod = prod_imp or prod_call>0
    if not has_prod: cls="c"
    elif feed>0 and hand==0: cls="a"
    elif feed>0: cls="a+b"
    elif hand>0: cls="b"
    else: cls="d"
    rows.append(dict(file=f,cls=cls,xs_only=xs_only,prod_imp=prod_imp,prod_call=prod_call,feed=feed,hand=hand,test=f.split("/")[-1].startswith("test_"),tree=f.split("/")[2],dermods=dermods))
json.dump(rows,open(sys.argv[1],"w"),indent=1)
c=collections.Counter((r["cls"],r["xs_only"]) for r in rows if r["test"])
print("N test files:",sum(r["test"] for r in rows),"helpers:",sum(not r["test"] for r in rows))
for k,v in sorted(c.items()): print(k,v)
