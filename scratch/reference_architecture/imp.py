import ast, pathlib, sys, collections
root=pathlib.Path('.')
pat=sys.argv[1]
hits=collections.defaultdict(set)
for p in list(root.glob('orpheus/**/*.py'))+list(root.glob('tests/**/*.py'))+list(root.glob('tools/**/*.py'))+list(root.glob('*.py')):
    if '_build' in p.parts: continue
    try: t=ast.parse(p.read_text())
    except Exception as e: print('PARSEFAIL',p,e); continue
    for n in ast.walk(t):
        if isinstance(n,ast.ImportFrom):
            mod=n.module or ''
            if n.level:  # resolve relative
                pk=list(p.with_suffix('').parts[:-n.level]) if n.level>0 else []
                base='.'.join(p.parent.parts[:len(p.parent.parts)-(n.level-1)])
                mod=base+('.'+mod if mod else '')
            names=[a.name for a in n.names]
            full=[mod]+[mod+'.'+x for x in names]
            if any(pat in f for f in full): hits[str(p)].add((n.lineno,mod,','.join(names)))
        elif isinstance(n,ast.Import):
            for a in n.names:
                if pat in a.name: hits[str(p)].add((n.lineno,a.name,''))
        elif isinstance(n,ast.Constant) and isinstance(n.value,str) and pat in n.value and ('import' in n.value or n.value.startswith('orpheus')):
            hits[str(p)].add((n.lineno,'STRING',n.value[:80].replace('\n',' ')))
for k in sorted(hits):
    print(k); 
    for h in sorted(hits[k]): print('   ',h)
print('files',len(hits))
