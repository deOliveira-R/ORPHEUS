import ast,re,sys,pathlib,json
SOLVE=re.compile(r'^(build_volume_kernel.*|K_vol_element_adaptive|solve_peierls_.*|solve_peierls_eigenvalue|_build_peierls_.*|_build_full_K_per_group|build_white_bc_correction.*|build_closure_operator|_build_closure_operator.*|BoundaryClosureOperator|continuous_case_builders|continuous_cases|_build|_class_[ab]_cases|build_two_surface_case|build_one_surface_compact_case|PeierlsSolution|_build_kernel_matrix|_build_system_matrices|per_observer_angular_assembly|per_surface_centred_angular_assembly|ClosureRecipe|_build_.*_op|_build_.*_PG|compute_K_bc_specular_continuous_mu_sphere|PeierlsSlabSolution)$')
files=[l.strip() for l in open(sys.argv[1])]
res={}
for f in files:
    src=pathlib.Path(f).read_text(); t=ast.parse(src)
    names=set(); modaliases={}
    for n in ast.walk(t):
        if isinstance(n,ast.ImportFrom) and n.module and 'peierls_nystrom' in n.module:
            for a in n.names:
                if a.name in ('geometry','cases','slab','cylinder','sphere','reference','specular'):
                    modaliases[a.asname or a.name]=a.name
                else: names.add(a.name)
        if isinstance(n,ast.Import):
            for a in n.names:
                if 'peierls_nystrom' in a.name: modaliases[a.asname or a.name.split('.')[-1]]=a.name
    for n in ast.walk(t):
        if isinstance(n,ast.Attribute) and isinstance(n.value,ast.Name) and n.value.id in modaliases:
            names.add(n.attr)
    # subprocess strings
    for n in ast.walk(t):
        if isinstance(n,ast.Constant) and isinstance(n.value,str) and 'peierls_nystrom' in n.value:
            for m in re.findall(r'import\s+\(?([^)\n]+)',n.value): 
                for x in re.split(r'[,\s]+',m): 
                    if x: names.add(x)
    reg=sorted(set(re.findall(r'peierls_(?:cyl1D|sph1D|slab)_[A-Za-z0-9_]+',src)))
    solve=sorted(x for x in names if SOLVE.match(x))
    other=sorted(x for x in names if not SOLVE.match(x))
    cls='NYSTROM-SOLVE' if (solve or reg) else ('PN-PRIMITIVE-ONLY' if names else 'NO-PN-IMPORT')
    res[f]=(cls,solve,other,reg)
for f,(c,s,o,r) in res.items(): print(c,f.replace('tests/gates/',''),'| solve:',s,'| other:',o,'| reg:',r[:3])
