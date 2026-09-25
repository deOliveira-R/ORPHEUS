import ast,re,sys,json,collections
SOLVE=re.compile(r'^(build_volume_kernel.*|K_vol_element_adaptive|solve_peierls_.*|solve_peierls_eigenvalue|_build_peierls_.*|_build_full_K_per_group|build_white_bc_correction.*|build_closure_operator|_build_closure_operator.*|BoundaryClosureOperator|continuous_case_builders|continuous_cases|build_two_surface_case|build_one_surface_compact_case|_build_kernel_matrix|_build_system_matrices|per_observer_angular_assembly|per_surface_centred_angular_assembly|solve_ps1982_vacuum_sphere|_build_.*_op|_build_.*_PG|compute_K_bc_specular_continuous_mu_sphere|continuous_get|_class_[ab]_cases)$')
REG=re.compile(r'peierls_(?:cyl1D|sph1D|slab)_')
out={}
for f in open(sys.argv[1]).read().split():
    src=open(f).read(); t=ast.parse(src)
    defs={}  # name -> node (module funcs, class methods keyed Class.meth and meth)
    for n in t.body:
        if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)): defs[n.name]=n
        if isinstance(n,ast.ClassDef):
            for m in n.body:
                if isinstance(m,(ast.FunctionDef,ast.AsyncFunctionDef)): defs[n.name+'.'+m.name]=m; defs.setdefault(m.name,m)
        if isinstance(n,ast.Assign):
            for tg in n.targets:
                if isinstance(tg,ast.Name): defs[tg.id]=n
    # module alias attr: geometry.solve_peierls_1g -> attribute name
    def direct(node):
        s=set(); hit=False
        for x in ast.walk(node):
            if isinstance(x,ast.Name): s.add(x.id)
            if isinstance(x,ast.Attribute): s.add(x.attr)
            if isinstance(x,ast.Constant) and isinstance(x.value,str) and (REG.search(x.value) or 'solve_peierls' in x.value): hit=True
        if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)):
            for a in node.args.args: s.add(a.arg)
        return s,hit
    memo={}
    def touches(name,stack=()):
        if name in memo: return memo[name]
        if name in stack or name not in defs: return False
        s,hit=direct(defs[name])
        r = hit or any(SOLVE.match(x) for x in s) or any(touches(x,stack+(name,)) for x in s if x in defs and x!=name)
        memo[name]=r; return r
    for n in ast.walk(t):
        if isinstance(n,ast.ClassDef):
            clsdec,_=direct(ast.Module(body=n.decorator_list,type_ignores=[]))
            for m in n.body:
                if isinstance(m,ast.FunctionDef) and m.name.startswith('test'):
                    out[f+'::'+n.name+'::'+m.name]=touches(n.name+'.'+m.name) or any(touches(p.arg) for p in m.args.args)
    for n in t.body:
        if isinstance(n,ast.FunctionDef) and n.name.startswith('test'):
            out[f+'::'+n.name]=touches(n.name) or any(touches(p.arg) for p in n.args.args)
json.dump(out,open(sys.argv[2],'w'))
c=collections.Counter()
for k,v in out.items(): c[(k.split('::')[0],v)]+=1
fs=sorted(set(k[0] for k in c))
for f in fs: print(f"{c[(f,True)]:4d} touch {c[(f,False)]:4d} not  {f.replace('tests/gates/','')}")
print('funcs',len(out),'touch',sum(out.values()))
