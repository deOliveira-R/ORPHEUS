# Issue #226 — pyright zero on data/micro_xs (hdf5_io.py + gendf.py)

## Pyright CLI errorCount — verbatim

BEFORE (combined oracle):
```
errorCount: 51
warningCount: 0
```

AFTER (combined oracle, same command):
```
errorCount: 0
warningCount: 0
```

Command: `npx --no-install pyright --outputjson orpheus/data/micro_xs/hdf5_io.py orpheus/data/micro_xs/gendf.py`

## Import smoke (-O)
```
ok
```
`.venv/bin/python -O -c "import orpheus.data.micro_xs.hdf5_io, orpheus.data.micro_xs.gendf; print('ok')"`

## Tests
- `tests/data/test_mixture.py` + `tests/data/test_cross_section_data.py`: 15 passed (these don't
  directly drive gendf/hdf5_io; `micro_xs` there is a local numpy var, `test_cross_section_data`
  imports `isotope` only).
- Added functional verification (ad-hoc, real U_235.GXS): `convert_gxs` → `save_isotope` →
  `load_isotope_h5` round-trips BIT-IDENTICAL on eg/sigC/sigT/nubar/chi/sig2/sigS. Runtime types
  confirmed: `sig2` and every `sigS[j][k]` are genuine `csr_matrix` (matches the Isotope contract).

## Error patterns & fixes
1. h5py untyped __getitem__ (Group|Dataset|Datatype): `typing.cast(h5py.Group/Dataset, ...)` at each
   site, chosen per the surrounding usage (further-indexed/.keys() → Group; sliced/read → Dataset).
   attrs reads cast to `np.floating` for the `float(...)` call. Runtime-inert (cast emits no code).
2. csr_array vs csr_matrix: `coo_matrix(...).tocsr()` returns a real csr_matrix at runtime (matrix
   lineage, `_csr_container == csr_matrix`), but the UNTYPED `tocsr()` body makes pyright infer
   csr_array. Cast `.tocsr()` results to `csr_matrix` (the Isotope contract, unmodifiable — isotope.py
   out of scope). Three sites: hdf5_io `_load_sparse` return; gendf `sig2`; gendf `_init_scattering`.

## gendf possibly-unbound (L160/164/165/168) — bug analysis
`n_lgn`/`n_sig0` in `_extract_mf6` are bound ONLY inside `if m[i,9]==1:` (the section-header / SEQ==1
record), which ALSO increments `n_temp`. They are read at the `if n_temp == temp_idx:` block. Since
every caller passes `temp_idx = i_temp` from `enumerate(temps, start=1)` (so temp_idx >= 1), a read
requires `n_temp >= 1`, which is only reachable AFTER the header bound both vars. So pyright's flow
gap is real structurally but the bad path is unreachable on well-formed input.

Empirically confirmed across 7 isotopes (H_001, U_235, U_238, O_016, B_010, ZR090, NA023): the first
MF=6 record for EVERY extracted MT has SEQ==1 — the header-first invariant holds. **Could NOT have
produced a wrong result at runtime** on real GENDF data; on a malformed/truncated file it would have
raised UnboundLocalError (a crash, not a silent wrong answer — not catches-worthy).

Fix: pre-initialize `n_lgn = n_sig0 = 0` before the loop with a comment documenting that the defaults
are provably never read (guarded by `n_temp == temp_idx >= 1`). Zero runtime change; makes the flow
statically sound and removes the dependent L160 operator error.

## Edit list (file:line: what + why)
- hdf5_io.py:33 — add `from typing import cast` (enable runtime-inert h5py/scipy casts).
- hdf5_io.py:71 — `cast(h5py.Group, f[...])` (temp group is indexed deeper / .attrs).
- hdf5_io.py:73-74 — `float(cast(np.floating, grp.attrs[...]))` (attrs union → ConvertibleToFloat).
- hdf5_io.py:75-85 — `cast(h5py.Dataset, grp[name])[:]` per dataset read (eg/sig0/sigC/sigL/sigF/sigT/nubar/chi).
- hdf5_io.py:90 — `cast(h5py.Group, grp["sigS"])` (indexed deeper + .keys()).
- hdf5_io.py:119-122 — `_load_sparse`: cast parent[name]→Group, row/col/data→Dataset.
- hdf5_io.py:123 — cast `coo_matrix(...).tocsr()` return → csr_matrix (true at runtime; satisfies annotation).
- gendf.py:20 — add `from typing import cast`.
- gendf.py:142 (pre-loop) — init `n_lgn = n_sig0 = 0` + rationale comment (possibly-unbound fix).
- gendf.py:~279 — annotate local `sig2: csr_matrix`; cast the `.tocsr()` branch → csr_matrix.
- gendf.py:~372 — `_init_scattering`: cast the `.tocsr()` assignment → csr_matrix.

## Scope/constraints honored
Only the two named files edited. No commits, no git checkout/restore/stash/reset. No isinstance,
no type:ignore. All casts truthful (verified runtime types). Only behavioral-adjacent change is the
provably-dead `0` defaults — exact runtime behavior preserved.
