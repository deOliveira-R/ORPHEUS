---
name: phase4-46-nonvacuum-mms-ansatz
description: Durable MMS-ansatz lesson — the non-vacuum prescribed-inflow MMS (4.6, LANDED). Every prior SN MMS ansatz VANISHES at boundaries → the q.boundary≠0 inflow path was UNTESTED; the a0>0 anti-bias ansatz, why a0>0 is load-bearing, and the converged-VALUE-not-just-rate check.
metadata:
  type: feedback
---

The non-vacuum prescribed-inflow MMS (Phase 4 O.2b 4.6) is LANDED on `main`
(`3273010`→`d87f3b0`). Tests live at
`tests/sn/verification/analytical/test_mms_prescribed_inflow.py`,
`tests/sn/verification/analytical/test_prescribed_inflow_consistency.py`,
`tests/derivations/test_sn_mms_nonvacuum_symbolic.py`. This note keeps the
durable ansatz-design WHY.

**THE GAP this MMS closed (the reason the ansatz had to be NON-standard).**
EVERY pre-existing SN MMS ansatz VANISHES at the boundary by construction —
slab/sphere/cyl, isotropic and anisotropic, all of `mms/sn.py`: `A=sin(πr/R)→0`
at the face; even the anisotropic `B=(r/R)(1−r/R)cos→B(R)=0`. A vanishing
boundary value means the incoming partial current γ₋ψ ≡ 0 automatically →
the `q.boundary ≠ 0` prescribed-inflow path was NEVER exercised by any MMS.
Prescribed inflow is a `q.boundary` SOURCE SLOT (`BoundarySourceSink` +
`face_view(face)[inflow_indices]=ψ_chosen`) consumed DIRECTLY by `(L+C).solve`
as the sweep inflow seed — NOT a `from_spec`/realizer bridge.

**THE ANSATZ (the anti-simplification-bias choice, vv Mode 7).**
`ψ_chosen = (A(x) + μ·B(x)) / W` with **A = a0 + a1·sin(k_x x), a0 > 0** and
B = b0·cos(k_x x). The φ_chosen = ∫ψ dμ = A(x) (since Σw_n μ_n = 0). **The
a0 > 0 is THE load-bearing bit:** it makes A(face) > 0 → γ₋ψ ≠ 0 → the inflow
slot is genuinely non-zero. NEVER default to `sin(πx/L)` here — it would
re-create the vacuum-automatic gap the MMS exists to close. This is the
concrete override of the MMS simplification bias: the "simplest trig that
satisfies the BCs" (sin, vanishing) is exactly the one that tests NOTHING about
the inflow path.

**THE LOAD-BEARING ASSERTION: converged-VALUE matches φ_chosen, not just the
O(h²) rate.** A solver that silently DROPS `q.boundary` still converges
cleanly at O(h²) — to the VACUUM flux. The rate is necessary-NOT-sufficient
(anti-pattern #5). ONLY the converged-value-vs-A(x)-with-a0>0 check sees the
dropped inflow. Always pair the rate assertion with a value-to-floor assertion
against the imposed φ_chosen, plus an inflow-honoured spot-check.

**Term activation / nulling (declare BOTH per Mode 7).** The slab ansatz
ACTIVATES: streaming μ(A'+μB'), c<1 scatter, the μ² 2nd-moment, the NON-VACUUM
γ₋ψ (the 4.6 novelty). It NULLS: curvilinear angular redistribution (→ a sphere
companion is MANDATORY, NEVER ship slab-only — sphere activates
`(1−μ²)/r ∂ψ/∂μ` with the μ-linear term giving ∂ψ/∂μ=B≠0) and fission (the
mixtures are non-fissile, no keff — MMS does NOT prove eigenvalues). The sphere
companion (T3) ships `xfail(strict, #195)` + `catches(ERR-026)` because the
curvilinear DD convergence it rides on is the still-OPEN #195 gate
(pre-asymptotic L2 + pole-face closure); it flips green when #195 closes. The
≥2G companion (T2) is MANDATORY per Cardinal-6 and exercises the Sᵀ
group-transfer (ERR-002).

**Structural independence (the chain of trust).** φ_chosen is IMPOSED (no
solve); Q_ext is SymPy-derived from the CONTINUOUS operator (new
`derive_nonvacuum_slab_mms`/`_sphere` + a foundation symbolic gate; Branch-2
numerical lambdified bit-equal ≤1e-13); the solver shares NO project-internal
primitive with the reference (numpy arrays only); the chain terminates in
SymPy + the imposed solution, NOT another orpheus solver.

**The latent-dud guard that bit me before (G-4):** an explicitly-reflective
`Mesh2D(bc_ymin=BC('reflective'))` constructed via a dataclass field, NOT a
string kwarg — a reflective mesh built with a `bc='vacuum'` KWARG silently
IGNORES the kwarg (explicit mesh BC wins) and re-solves the wrong problem. The
consistency probe `diag_p46` ran green pre-carve (1-D SI≡Krylov 7.8e-14;
2-D vacuum-x + reflective-y + inflow-xmin SI-Jacobi≡SI-G-S≡Krylov ≤5.6e-13 —
vv Mode-9 for prescribed-inflow confirmed) and was promoted into the
consistency test, then deleted.

See [[issue-208-wave-o-carve-lessons]] (this MMS built on the FINAL forward
shape Wave-O landed), [[regression-tolerance-design]].
