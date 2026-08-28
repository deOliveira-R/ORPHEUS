# P4.9a §6b enumeration — call sites and shapes for the unweld (rows 1-3 + the ruled cache edit)

Measured 2026-08-28 at HEAD `5c4f56d7` (P4.3 `da507e3d` is an ancestor; branch `main`).
Method: whole-file enumerations (never `| head`); membership by **Python AST**
(`.venv/bin/python`, script + raw output in the session scratchpad:
`p49a_ast.py` / `p49a_ast_out.txt` / `p49a_triage.py`); 874 files parsed
(orpheus/ + tests/); 4/4 positive controls green (CellVisit ctor in
augmented_mesh, `visit.tau` in diamond, `mm_a_in_coeff` in cache,
`angular_upstream` in cell_balance). String-form residual
(`grep -rnE "['\"]member['\"]"`) run for every shed member + cache constant.
`mcp__nexus__*` tools were PRESENT in the session; the membership questions
below are AST-derived per the brief's mandate, with grep used only for
literal-text location.

## 0 · HEADLINE — where the tree contradicts the brief

1. **Path correction.** The walk lives at
   `orpheus/sn/loss_representation/__init__.py` (4959 lines). There is NO
   `orpheus/sn/operators/loss_representation/` — `orpheus/sn/operators/`
   holds streaming.py, radial_characteristic.py, etc.

2. **Row 1 is incomplete inside diamond.py itself: `DiamondDifference.residual`
   is a SECOND consumer of the shed members** (`visit.c_out` :289,
   `visit.c_in` :290, `upstream_state.angular_upstream` :294). It ALREADY
   routes through `cell_balance_for_streaming` (:304) — but it assembles
   `angular_denom_term = dA_w·c_out` and `angular_numer_upstream =
   dA_w·c_in·ψ_ang` ITSELF from the protocol members (:288-303). Shedding
   row 3's members forces `residual`'s re-plumb in the same step as `update`'s
   (§6b: one interface, one call-site set). Mitigating fact, `[M]` by AST:
   production constructs `UpstreamState` at exactly ONE site
   (`loss_representation/__init__.py:4285`, feeding `update` :4291) ⟹
   **`scheme.residual` and `scheme.update` have zero OTHER production
   callers** — `residual`'s consumers are tests (test_diamond) + the protocol
   contract.

3. **Killing diamond's block does NOT make the M-M relation single-spelled.**
   Production spellings of the update relation after row 1, all outside
   `transport/` (so the done-when grep stays precise AND stays silent about
   them):
   - owner march: `closure.py:1328-1332` (`(ψ_m − (1−τ)ψ_half)/τ`);
   - `_run` fast path: `__init__.py:4348-4352` — `ψ_out = geom.tau_inv·ψ̄ −
     geom.mm_a_in_coeff·ψ_in` (same relation via cache constants);
   - `_run_transpose` hand-written adjoints: `:4653-4655` (degenerate) and
     `:4696-4698` (fast path);
   - the NUMERATOR contribution `(ΔA/w)·c_in·ψ` — the same formula
     `MorelMontryAngularSweep.cell_contribution` owns — is spelled inline at
     `:4306-4308` (`geom.dA_w·geom.c_in`), `:4665` (`kappa`), `:4728`
     (`ang_coeff`).
   `[M]` cache.py's own comment (:361-364) names the scan fast path **"the
   twin of DiamondDifference.update's (ψ̄ − (1−τ)ψ_in)/τ"** — the tree
   already documents the surviving spelling as a twin. Row 1+2 remove the
   `transport/` copy; the walk-side spellings are the "+ruled" edit's and
   §5b/O-3's territory. Enumerated so the phase design decides with open eyes.

4. **The "+ruled" cache edit overturns a RECORDED in-tree rationale.**
   `cache.py:344-371` defends the current split at length: "the closure
   exposes only the PRIMITIVE τ (Pattern 5 — build the primitive, not the
   product)", L16 perf-hoisting, and a byte-identity provenance citing the
   "Step-A Leg-1 producer-equivalence gate". Handing the cache its march
   constants inverts that ruling — the comment block, its Pattern-5 claim,
   and the byte-identity provenance all need rewriting in the same commit.

5. **`c_out` has a THIRD production surface row 3 does not name: the
   `affine_scan_coefficients` kwarg** — base `scheme.py:1394-1406`
   (`c_out: np.ndarray` param), DD override `diamond.py:572` (assembles
   `dA_w * c_out[:, None]` at :676 — an M-M product computed INSIDE
   transport/spatial), LD override `linear_discontinuous.py:759` (shape
   (N, nx), used only in the :802-807 refusal). Sole production caller:
   `cache.py:532-539` passing `c_out=geom.c_out`. The done-when ("no
   tau/c_in/c_out member on any protocol TYPE") does not cover a method
   PARAMETER; the charter's direction-supplier ruling ("assembled face
   coefficients") points at re-posing it (`c_out=`+`dA_w=` →
   `angular_denom_term=`), but that is a design decision, not enumerated
   fact.

6. **Line drift vs the brief:** diamond's M-M comment is at `diamond.py:219`
   (block :219-230; brief said ~:229). `cell_balance_terms` at
   `cell_balance.py:266` ✓. `cache.py:377` (`mm_a_in_coeff`) ✓. The
   degenerate branch's `:4291`/`:4298` ✓ unchanged. LD refusal raise
   :147-153, guard def :140 (brief said :148-153) ✓.

7. **Test-surface triage:** 3 of the brief's 12 files have **zero**
   protocol-member sites (`test_mms_ld_slab`, `test_ld_ubld_primitive`,
   `test_scheme_reaction_rate_contract` — all scan-API `c_out=` only). One
   file NOT on the list carries a cache-constant read
   (`test_loss_transpose_solve.py:404`, `geom.mm_a_in_coeff`), and
   `test_cache.py:172` reads the cache constants by STRING (a field-name
   loop) — invisible to any attribute grep.

---

## 1 · CellVisit construction sites + who computes the passed values

### Production — exactly ONE constructor
| site | fills | value provenance (one hop up) |
|---|---|---|
| `orpheus/sn/mesh/augmented_mesh.py:1668-1675` (`SNMesh._make_cell_visit`, def :1636) | `c_in`, `c_out`, `tau` (+ cell_idx, streaming_terms, face_area_downstream) | `self.pole_angular_closure.c_in_per_ordinate[g]` :1672, `.c_out_per_ordinate[g]` :1673, `.tau_per_ordinate[g]` :1674 — base accessors `closure.py:463/:473/:483`; MM precomputes per level (`c_in = (1−τ)/τ·α_out + α_in` at `closure.py:1593`) from its α-dome/τ at construction; Identity returns neutral (c=0, τ=1) |

Yield sites of `_make_cell_visit`: `augmented_mesh.py:1700, :1733, :1784, :1808`
(all inside `dag_walk`, def :1374). No other production `CellVisit(...)` exists
(`[M]` AST, 874 files).

Angular-member FILLS that are constructions of the sibling types:
- `loss_representation/__init__.py:4285-4288` — `UpstreamState(spatial_upstream=psi_in, angular_upstream=psi_angle[:, i])`; `psi_in` from inflow / zeros / the coupled-pole mirror (:4258-4274); `psi_angle` is the per-level M-M thread buffer (zeroed :4245, written by :4298 and :4353).
- `diamond.py:232-236` — `CellResult(outgoing_angular_state=psi_angle_out)`, computed :225-230 from `visit.tau` + `upstream_state.angular_upstream`.
- `linear_discontinuous.py:379-383` — `CellResult(outgoing_angular_state=None)` literal.

### Tests — 27 `CellVisit(...)` lines (`[M]` AST census)
| file | ctor lines | pass angular kwargs? | value provenance |
|---|---|---|---|
| test_cell_balance_for_streaming.py | 468, 513 | c_in,c_out,tau both | fixture literals (`_CURVILINEAR_TAU` etc.) |
| test_diamond.py | 221, 284, 343, 840 | none (slab arms) | — |
| test_diamond.py | 436, 564, 689, 767 | c_in,c_out,tau | `ref_c_in`/`ref_c_out` re-derived in-test from α/τ literals |
| test_diamond.py | 932, 979, 1032, 1086, 1507, 1512 | c_in,c_out | fixture literals |
| test_discretization_scheme_protocol.py | 445 (tau), 478, 497, 509 | tau only at :445 | literal |
| test_ordinate_scan.py | 533, 570, 618 | c_in,c_out,tau | fixture literals |
| test_dd_recurrence.py | 102 | none | — |
| test_ld_ubld_primitive.py | 382 | none | — |
| test_ld_ubld_symbolic.py | 154 | none | — |
| test_linear_discontinuous.py | 66 | none | — |
| test_cache.py | (none — consumes mesh-stamped visits via `dag_walk`) | reads at :509/:510 | mesh stamp |

---

## 2 · The degenerate branch, exactly

`_OneDimScanWalk` at `orpheus/sn/loss_representation/__init__.py:2815`; the
branch is in `_run` (def :3864). Verbatim, current line numbers:

```
 4276                     # Degenerate cyl-axis ordinate: slow per-cell path.
 4277                     if geom.is_degenerate[global_n]:
 4278                         ordinate_idx = global_n if is_sphere else m_local
 4279                         visits = list(self.mesh.dag_walk(
 4280                             ordinate_idx=ordinate_idx,
 4281                             mu_level_idx=level,
 4282                         ))
 4283                         for visit in visits:
 4284                             i = visit.cell_idx
 4285                             upstream = UpstreamState(
 4286                                 spatial_upstream=psi_in,
 4287                                 angular_upstream=psi_angle[:, i],
 4288                             )
 4289                             # scheme.update expects per-cell (ng,)
 4290                             # arrays — sig_t / source slice on the cell axis.
 4291                             result = scheme.update(
 4292                                 visit=visit,
 4293                                 total_xs=sig_t_p[:, i],
 4294                                 source=QV_full[:, i],
 4295                                 upstream_state=upstream,
 4296                             )
 4297                             psi = result.cell_average_flux           # (ng,)
 4298                             psi_angle[:, i] = result.outgoing_angular_state
 4299                             angular_flux[global_n, :, i] = psi
 4300                             scalar_flux[:, i] += w_n * psi
 4301                         continue
```

What it does: builds a mesh-stamped visit per cell (`dag_walk` →
`_make_cell_visit`, so `c_in`/`c_out`/`tau` arrive stamped from the closure),
hands the M-M thread value `psi_angle[:, i]` in as `angular_upstream`, and
threads `result.outgoing_angular_state` BACK into `psi_angle[:, i]` (:4298) so
the level's next ordinate reads it (:4287) — i.e. the angular march currently
round-trips through the spatial scheme's CellResult. The related degenerate
`psi_in` slot: `elif geom.is_degenerate[global_n]: psi_in = np.zeros(ng)`
(:4260-4265). ⚠ The TRANSPOSE twin (`_run_transpose`, def :4373) has its own
degenerate branch at `:4652-4683` which does NOT call the scheme at all — it
hand-writes the M-M adjoint from cache constants (`tau_inv` :4654,
`mm_a_in_coeff` :4655, `dA_w·c_in` :4665). Any carve of the forward branch
owes the transpose branch a matching answer (today they are already
asymmetric: forward = scheme.update, transpose = inline cache algebra).

**(a) Where the walk gets its scheme.** `_OneDimScanWalk` holds NO scheme
field — it reads `self.mesh.scheme` at `:3921` (`_run`) and `:4420`
(`_run_transpose`), plus `self.mesh.scheme.source_emission/cell_average` at
:4040/:4057/:4318/:4338 and protocol-trait reads in `_apply_walk`
(:3159/:3217/:3221/:3259). `sweep_graph.py`'s `scheme` fields belong to the
2-D DAG cell-op objects — `_CellSolve` (ABC, field :872; subclasses
`_CellSolveAngular` :950 / `_CellSolveMoment` :967), `_CellResidual` (field
:992), `_CellResidualTranspose` (field :1092) — and are FILLED at 7 sites in
the 2-D `_DAGWavefront` strategies, all spelling `scheme=self.mesh.scheme`:
`loss_representation/__init__.py:1635, :1644, :1723, :1807, :2023, :2094,
:2184`. Both roots are `SNMesh.scheme`, filled at
`augmented_mesh.py:268-270` (default `DiamondDifference()`). The 2-D path
consumes only `cell_kernel_batch`/`residual_kernel_batch`(+`_transpose`) —
no `update`, no angular member.

**(b) Where the VECTORIZED branch gets its closure.**
`SNMesh.pole_angular_closure` is built at `augmented_mesh.py:404-424`:
`closure_cls = self._user_supplied_closure or
default_angular_closure_class(self.coord)`, constructed
`closure_cls(angular, pairing)` at :422 (operands from `self.reduced` on 1-D,
or built from `(quad, coord)` on multi-D Cartesian). The matvec
(`_apply_walk`, def :3072) fetches it at `:3147`, calls
`precompute_psi_state(psi_view, radial_characteristic=…)` at `:3195-3200`,
and consumes the state via `cell_contribution` at `:3285` (leg body) and
`:3359` (its OWN degenerate loop — note the matvec already handles degenerate
ordinates through the closure, not through scheme.update). The SOLVE walk
(`_run`) touches the closure OBJECT only as the admission guard
(`isinstance(closure, MorelMontryAngularSweep)` :4212-4218); its numeric path
rides the CACHE constants that `StreamingCoefficientCache.from_mesh_and_quad`
derives from the closure's accessors (cache.py:372-377). ⟹ the planned
`is`-identity gate is well-posed: after the carve, both branches' closure is
`self.mesh.pole_angular_closure` (today the degenerate branch reaches
diamond's inline copy — no closure object — and the gate cannot pass, the
§6c red).

---

## 3 · `closure.cell_contribution` — signature, returns, ALL call sites

Abstract: `closure.py:540-556` —
`cell_contribution(psi_state: Any, cell_idx: int, level_idx: int,
within_positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]`, returning
`(denom_term (n_mask,), upstream_numer_term (ng, n_mask))` = the closure's
additions `(ΔA/w)·c_out` and `(ΔA/w)·c_in·ψ_{m−1/2}` to the cell balance.
Overrides: `MorelMontryAngularSweep` :1786-1820 (reads `self._dAw_per_level`
/ `_c_in_per_level` / `_c_out_per_level` + `psi_state[p].upstream_per_ordinate`
— no M-M argument passed in); `IdentityAngularClosure` :2113-2127 (returns
`(zeros(n), zeros((1, n)))`).

| call site | which half used | context |
|---|---|---|
| `sn/loss_representation/__init__.py:3285` | both | `_apply_walk` leg body (curvilinear matvec arm) → `cell_balance_for_streaming` :3290 |
| `sn/loss_representation/__init__.py:3359` | both | `_apply_walk` degenerate-ordinate loop → `cell_balance_for_streaming` :3367 (A_downstream=0.0) |
| `sn/loss_representation/__init__.py:3653` | denom only | `loss_action_transpose` leg body (ψ-independent coefficient use) |
| `sn/loss_representation/__init__.py:3732` | denom only | `loss_action_transpose` degenerate loop |
| `sn/operators/radial_characteristic.py:895` | numer only | `RadialCharacteristicSeeding.apply` (class :714, def :831) — seed propagation with zeroed bulk (precompute at :887) |
| tests: `tests/sn/operators/test_psi_half_coupling.py:2590` | spy | `_install_closure_spy(…, "cell_contribution")` call-count gate |
| tests: `tests/sn/operators/test_psi_half_coupling.py:2701-2707` | monkeypatch | sign-flip TOOTH for the adjoint gate |

(Prose-only mentions: cell_balance.py:149/:230/:233; test module docstrings
test_cell_balance_for_streaming:13, test_pole_angular_closure:19/:82,
test_cell_visit_c_stamp:21.)

---

## 4 · Every PRODUCTION read of the shed members (`[M]` AST, 874 files)

| member | file:line | context |
|---|---|---|
| `CellVisit.tau` | `transport/spatial/diamond.py:227` | `update`'s M-M block (the twin) |
| `CellVisit.c_in` | `diamond.py:199` | `update` → `cell_balance_terms(c_in=)` kwarg |
| `CellVisit.c_in` | `diamond.py:290` | `residual` — assembles numer contribution itself |
| `CellVisit.c_out` | `diamond.py:200` | `update` → `cell_balance_terms(c_out=)` kwarg |
| `CellVisit.c_out` | `diamond.py:289` | `residual` — assembles denom contribution itself |
| `UpstreamState.angular_upstream` | `diamond.py:226` | `update` — `is not None` geometry branch |
| `UpstreamState.angular_upstream` | `diamond.py:229` | `update` — the M-M formula |
| `UpstreamState.angular_upstream` | `diamond.py:294` | `residual` — numer assembly (None ⇒ zeros) |
| `UpstreamState.angular_upstream` | `transport/spatial/cell_balance.py:352` | `cell_balance_terms` (dies with row 2) |
| `UpstreamState.angular_upstream` | `transport/spatial/linear_discontinuous.py:148` | `_require_slab` — presence-as-signal refusal (see §5) |
| `CellResult.outgoing_angular_state` | `sn/loss_representation/__init__.py:4298` | `_run` degenerate branch (the §2 quote) |

Fills (constructions), for completeness: `augmented_mesh.py:1672-1674`,
`loss_representation/__init__.py:4285-4288`, `diamond.py:232-236`,
`linear_discontinuous.py:379-383`.

**Explicitly NOT protocol members — the CACHE's own fields with the same
names** (receiver `geom = StreamingCoefficientCache`): `geom.c_in`
`__init__.py:4307, :4665, :4728`; `geom.c_out` `cache.py:537`;
`geom.tau_inv` `:4349, :4654, :4698`; `geom.mm_a_in_coeff` `:4350, :4655,
:4697`. These are §6's surface, untouched by row 3.

String-form residual: `[M]` zero production `getattr`/string reads of any
shed member (`grep -rnE "['\"](c_in|c_out|angular_upstream|outgoing_angular_state)['\"]" orpheus/`
→ 0 production hits; the only `"tau"` strings are sympy `Symbol("tau")` in
`orpheus/derivations/` — different referent).

---

## 5 · LD's `angular_upstream` read — a presence SIGNAL, not a pass-through

`linear_discontinuous.py:140-155`, `_require_slab`:
`if upstream_state.angular_upstream is not None: raise NotImplementedError(…)`.
Called from `update` :375 and `residual` :401. Its own docstring: "the
presence of `angular_upstream` is therefore the geometry gate — identical to
the signal DD uses". On the slab the member is `None` ⟹ the guard passes and
the VALUE is never consumed; `update` returns
`CellResult(outgoing_angular_state=None)` (:379-383) unconditionally. So: not
a pass-through — the member's only LD use is its EXISTENCE as a curvilinear
detector. Consequences for row 3:

- Shedding the member deletes `_require_slab`'s signal. The mid-sweep raise
  becomes unspellable — which is ACCEPTABLE only because the refusal is
  triple-guarded upstream: (i) walk admission
  `loss_representation/__init__.py:240` (`mesh.is_cartesian or
  mesh.scheme.supports_curvilinear`; LD `supports_curvilinear=False` at
  :273), (ii) the scan-path guard `affine_scan_coefficients` :802-807 keyed
  on `dA_w`/`c_out` non-neutrality (`np.any(dA_w != 0.0) or np.any(c_out !=
  0.0)`), which itself survives only as long as that kwarg surface survives
  (§0.5). Design must decide `_require_slab`'s fate (retire vs re-key);
  its two negative-leg tests (`test_linear_discontinuous.py:242`,
  `test_mms_ld_slab`-adjacent none) lose their trigger with the member.

---

## 6 · Cache constants: producers and consumers (scopes the "+ruled" edit)

Fields: `StreamingCoefficientCache` `cache.py:239-242` — `c_in (N,)`,
`c_out (N,)`, `tau_inv (N,)`, `mm_a_in_coeff (N,)`.

Populate path (`from_mesh_and_quad`, def :254): `closure =
sn_mesh.pole_angular_closure` :372; `tau = closure.tau_per_ordinate` :373;
`c_out = closure.c_out_per_ordinate` :374; `c_in = closure.c_in_per_ordinate`
:375; **`tau_inv = 1.0 / tau` :376; `mm_a_in_coeff = (1.0 − tau) / tau`
:377** — the two derivations the ruled edit relocates into the closure;
stamped into `cls(…)` :384-400. `is_degenerate` set at :359-364.

| field | production consumers | tests |
|---|---|---|
| `mm_a_in_coeff` | `__init__.py:4350` (`_run` fast-path M-M update), `:4655` (`_run_transpose` degenerate adjoint), `:4697` (`_run_transpose` fast adjoint) | `test_loss_transpose_solve.py:404`; `test_cache.py:114`, `:172` (string-form field loop) |
| `tau_inv` | `:4349`, `:4654`, `:4698` (same three arms) | `test_cache.py:113`, `:172` |
| `c_in` | `:4307` (`_run` fast-path `ang_contrib = dA_w·c_in·ψ`), `:4665` (`kappa`), `:4728` (`ang_coeff`) | `test_cache.py:111, :412, :172` |
| `c_out` | `cache.py:537` — `CollisionCache.from_geometry` → `scheme.affine_scan_coefficients(c_out=geom.c_out)` :532 (the ONLY production read; DD assembles `dA_w·c_out` at `diamond.py:676`) | `test_cache.py:112, :172` |

Scope of the ruled edit: `:376-377` die and the closure exposes the march
constants (the fields and every consumer above stay put); the rationale block
`:344-371` (Pattern-5/L16 defence + "byte-identical to HEAD" provenance +
Step-A Leg-1 gate citation) must be rewritten — it argues the design being
overturned (§0.4).

---

## 7 · Row 3's TEST surface, triaged by MEANING

Predicate: a **real site** = a shed-member kwarg occurrence on a
CellVisit/UpstreamState/CellResult construction (incl. `angular_upstream=None`
— deleting the field makes the kwarg a `TypeError`), OR an attribute read of a
shed member, OR a `cell_balance_terms(c_in=/c_out=)` call / `CellBalanceTerms`
c-field read (dies with row 2). Counted as kwarg/read OCCURRENCES; distinct
lines in brackets. `[M]` AST census + collation (`p49a_triage.py`).

| file | verdict | real sites [lines] | false positives / notes |
|---|---|---|---|
| tests/sn/sweep/core/test_diamond.py | **the heavy migration** — the DD M-M behaviour suite exercises the closure THROUGH `update`; re-target to the closure object | 54 [ctor: 204,271,330,402,436,532,564,670,689,754,758,767,835,923,932,970,979,1019,1032,1081,1086,1500,1507,1512,1547; reads: 242,299,448,579,711,713,1423,1424,1465,1466,1546; row-2: 1418,1460] | none |
| tests/sn/sweep/core/test_ordinate_scan.py | real + row-2 twin gate (:64 family) | 17 [83,94,99,100,107,533,570,618,711,818] | none |
| tests/sn/sweep/core/test_cell_balance_for_streaming.py | real + the row-2 equivalence gate loses its subject (memo §"twin-equivalence") | 16 [143,147,167,178,215,225,468,475,513,519] | `_scalar_to_vectorized_inputs(tau=)` :185,:232,:487 — test-LOCAL helper kwarg (def :114), not protocol |
| tests/sn/sweep/core/test_discretization_scheme_protocol.py | real — the PROTOCOL contract tests; the angular slots of the contract are the direct subject (stub schemes at :89/:192 are duck-typed consumer classes) | 16 [89,168,179,182,188,190,192,350,354,358,369,381,445,454,465,466] | none |
| tests/transport/spatial/test_linear_discontinuous.py | real — mostly `angular_upstream=None` kwargs + the :242 refusal trigger | 10 [116,131,157,172,193,242,255,320,452,457] | `affine_scan_coefficients(c_out=)` :366 — scan-API (§0.5), not row 3 |
| tests/sn/sweep/core/test_cache.py | real (visit reads + row-2 call) + cache-field surface | 6 [438,500,504,509,510] | `geom.c_*` :111-114,:412 and the string loop :172 are CACHE fields (§6), not protocol |
| tests/sn/sweep/core/test_cell_visit_c_stamp.py | real — the file's SUBJECT is the stamp (mesh copies closure values onto visits); row 3 retires the stamp itself, so the file re-derives as a closure-accessor contract | 5 [108,112,116,233] | none |
| tests/sn/sweep/slab/test_dd_recurrence.py | trivial | 1 [99: `angular_upstream=None`] | CellVisit :102 passes no shed kwarg |
| tests/transport/spatial/test_ld_ubld_symbolic.py | trivial | 1 [178: `angular_upstream=None`] | — |
| tests/sn/verification/mms/test_mms_ld_slab.py | **false positive for row 3** | 0 | :269/:277/:278 build `affine_scan_coefficients` kwargs (dict `c_out=`) — the scan-API seam |
| tests/transport/spatial/test_ld_ubld_primitive.py | **false positive for row 3** | 0 | :446 `affine_scan_coefficients(c_out=)`; CellVisit :382 / UpstreamState :402 pass no shed kwarg |
| tests/transport/spatial/test_scheme_reaction_rate_contract.py | **false positive for row 3** | 0 | :116,:126 `affine_scan_coefficients(c_out=)` only |
| (not in the brief's 12) tests/sn/operators/test_loss_transpose_solve.py | cache-constant surface only | 0 | :404 `geom.mm_a_in_coeff` — §6's edit, not row 3 |

Total real row-2+row-3 occurrences across the surface: **126**. The three
scan-API-only files become real §6b members ONLY if the phase re-poses
`affine_scan_coefficients`'s `c_out=`/`dA_w=` kwargs (§0.5).

Row-2 companion (from the pre-measurement memo, not re-derived here): the
twin-equivalence gates naming `cell_balance_terms`
(test_cell_balance_for_streaming:23, test_cache:474, test_diamond:1400/:1444,
test_ordinate_scan:64) lose their comparison subject — each is a
coding-standards rewire, not a delete.

