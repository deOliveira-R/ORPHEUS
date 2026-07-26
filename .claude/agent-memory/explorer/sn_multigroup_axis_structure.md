---
name: sn-multigroup-axis-structure
description: Durable shape of the SN energy-group axis — three-tier group-blindness, no group loop, "within-group" means fission-external; plus the AGENT.md durable-shape drift (within_group_triple → build_within_group_system)
metadata:
  type: project
---

# SN multigroup axis — durable shape (verified 2026-07-21 @ branch docs/sn-doc-architecture)

The SN solve has **NO energy-group loop** on any iterative path. The group axis
is a broadcast array axis of the state (`AngularFlux` canonical `(N, ng, nx[, ny])`;
Krylov `restart = n_dof = N·ng·nx·ny`). The only `for g in range(ng)` in
production sn/ is the dense **assembly** arm (`sn/loss_representation/assembly.py`),
which materializes one sparse block per group over the SAME walk graph — itself
evidence that the transport factor is group-diagonal.

**Three tiers** (the honest decomposition — sharper than "group-blind vs per-group"):

1. **Group-INDEPENDENT structure** (Stratum 1, no `ng` axis): walk order
   (`SNMesh.dag_walk` → `CellVisit` packets; levels/octants), the α-dome
   recursion `α_{n+1/2} = α_{n−1/2} − w_n μ_n` (quadrature-only,
   `geometry/reduced_operator.py`), Morel–Montry τ_n (a function of (μ,w) only —
   #236 ruling), the derived c_in/c_out constants `(N,)`, ΔA/w. In ORPHEUS
   vocabulary **τ/c are ANGULAR closure weights, NOT optical thicknesses** —
   do not confuse with Σ_t·Δr/|μ|.
2. **Group-DIAGONAL data** (Stratum 2, `(N, ng, nx)`, `sn/sweep/cache.py`):
   WDD denominator `2|μ|A_down + angular_denom_term[n] + Σ_t,g·V`
   (`transport/spatial/cell_balance.py::cell_balance_for_streaming` — the ONLY
   group injection is Σ_t,g·V), attenuation `a`, blend weight, and ψ½
   (`carlson_inward_sweep_from_source(Q̄(ng,nx), σ_t(ng,nx), …)` — cell loop,
   groups broadcast). Per-group data riding the group-blind walk; no g→g' mixing.
3. **Group-COUPLING** only in S and F. S = R∘(Λ+N2N)∘M
   (`transport/operators/scattering.py`; Λ = per-ℓ group-transfer matmul, "the
   only group-asymmetric factor"), applied as the lagged gain (SI) / matvec term
   (Krylov) of the ONE monolithic all-groups within-group system. F = rank-1
   χ⊗νΣ_f (`transport/operators/fission.py`), applied ONLY at the eigenvalue
   outer (`compute_fission_source = F.apply(φ)/k`; `power_iteration` in
   `numerics/eigenvalue`), never inside the within-group OperatorSum.

**Naming trap:** ORPHEUS "within-group system" (`build_within_group_system`,
`WithinGroupSystem` record, `orpheus/sn/coupled_system.py`) means
**fission-external**, NOT per-group — it contains the FULL multigroup S. There
is no group Gauss-Seidel / thermal iteration anywhere; the only G-S is over
sweep octants (`_select_si_resolvent`).

**The #196 gate** (= ERR-026 manifestation #7 catchers,
`tests/sn/eigenvalue/test_keff_curvilinear.py`:
`test_si_krylov_eigenvalue_equivalence_{sphere,cylinder}`) runs **2G
heterogeneous** by design ("not a homogeneous/1G degenerate"). Known real
curvilinear×multigroup interactions: sphere-4g-krylov xfail (unpreconditioned
GMRES, #200) and the curvilinear MG inner-tol amplification memo (header of
`tests/sn/verification/analytical/test_kinf_homogeneous.py`).

**AGENT.md drift note:** the durable-shape section still says the within-group
factory is `_within_group_triple` (solver.py). Since the B.2d coupled-block
campaign the ONE construction site is `build_within_group_system`
(`orpheus/sn/coupled_system.py`) returning the `WithinGroupSystem` record
(loss A, space, resolvent M, gains N; `A = M − N` splitting; 2×2 coupled grid
on seed-carrying meshes, 1×1 seedless). `_within_group_si` / `_within_group_krylov`
(solver.py) consume the record. Update AGENT.md's shape section when next revised.

**Second AGENT.md drift (verified 2026-07-26, DSA recon):** the composite field
ctor keyword is `FullField(interior=…, boundary=…)` — NOT `bulk=` as the
durable-shape section spells it (the ATTRIBUTE reads as `.interior`; prose
"bulk ⊕ boundary" is fine, the kwarg is not). Same recon confirmed the P7b
`full_field_space`-not-on-TransportMethod deferral note is still live
(`transport/method.py`, "declare when the DSA driver #2 arrives") and that
`KrylovAcceleration` HAS a `preconditioner=` seam while the SN production path
passes an explicit identity lambda (#200). Full transient map (line numbers):
`.claude/plans/dsa_landing_zone_recon.md` — re-derive after Phase 3 merges.
