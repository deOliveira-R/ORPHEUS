# Explorer memory index

The durable SHAPE of the SN operator-algebra subsystem now lives in
`.claude/agents/explorer/AGENT.md` ("SN operator-algebra subsystem — durable
shape"). The point-in-time `file:line` carve/surface maps for the #208 / Wave-O /
Phase 2-5 / WavefrontFlux / g-adjoint work were DELETED during the 2026-06-08
hygiene pass — that work is all merged to origin/main, so their line maps were
stale and re-derivable in seconds via Nexus. Keep notes here durable, not transient.

## Durable physics / convention context

- [HarmonicMomentField UNITS convention](harmonic_moment_field_units_convention.md) — why a stored SH moment carries SCALAR-flux units (no-prefactor SH, Y_0^0=1, weights sum to 4π → sr cancels); R≠M*; ERR-039/ERR-051 history.
- [Phase 5 µ-resolved primitive inventory](phase5_mu_resolved_primitive_inventory.md) — µ-resolved vs µ-integrated primitives in peierls_geometry.py for the continuous-µ specular multibounce closure.

## Committed pre-carve dependency audits (durable rationale; line maps may have drifted)

- [D-I.3 dependency audit](D-I.3_dependency_audit.md) — retiring the bare-ndarray adapter arm of StreamingOperator.apply; production+test caller inventory + readiness verdict.
- [A.5 BoundaryFlux→TraceSpace rehome audit](a5_boundary_flux_trace_space_rehome_audit.md) — re-homing BoundaryFlux onto the unified TraceSpace (FaceLayout field→space); consumer inventory + minimal-churn carve.
- [D-G BoundaryFlux consumer audit](dg_boundary_flux_consumer_audit.md) — every BoundaryFlux consumer (prod+tests) with WRITE/mutability sites called out; source of truth for the BoundaryFlux→pure-Field carve.
- [D-H.1b AngularFlux consumer audit](dh1b_angular_flux_consumer_audit.md) — legacy AngularFlux usage → new L2 AngularFlux + TimedFullField migration mapping; Phase-A overlap-risk ordering.
- [FunctionSpace typed-field audit](function_space_typed_field_audit.md) — FunctionSpace as the L1 base for every typed transport field; usage inventory + dunder-algebra patterns + migration targets.

## Architectural-rationale audits (subsystem shape; the WHY behind the typed algebra — mostly LANDED, kept for rationale)

- [Issue #196 SN operator architecture audit](issue_196_sn_operator_architecture_audit.md) — diagnoses SI-vs-Krylov drift as twin procedural impls of the same primitives; proposes the operator-algebra unification that makes twin-path defects structurally impossible.
- [Phase G typed field contracts](typed_field_contracts_for_phase_g.md) — canonical frozen-dataclass contracts for AngularFlux/ScalarFlux/Solution; the (L+C−S−F)ψ=q typed-leaf vocabulary; relationship to Issue #197.
- [Issue #196 Phase G replan — SN structure](issue_196_phase_g_replan_sn_structure.md) — SNMesh IS the SN phase space (no wrapper); four-operator layer status; structure-carve anti-recommendations.
- [Issue #196 Phase G replan — strategy Protocols](issue_196_phase_g_replan_sn_strategies.md) — CellUpdate / PsiHalfAngleSeed / PoleAngularClosure / BoundaryRealizer Protocols; residual math belongs on CellUpdate.residual.
- [Issue #196 Phase G replan — algebra inventory](issue_196_phase_g_replan_algebra.md) — inventory of the LinearOperator/OperatorSum/Scattering(R·Λ·M)/Fission/SourceIteration/KEigenvalue algebra; most already shipped.
- [Issue #196 Phase G Step 2.5 — DD polymorphism](issue_196_phase_g_step2_5_dd_polymorphism.md) — the three DiamondDifference branches (slab/curvilinear/cyl-degenerate); unified-update + clean (L+C)/closure split.
- [Issue #196 Phase G Step 2.5 — DAG-walk topology](issue_196_phase_g_step2_5_dag_walk_topology.md) — iter_cells_by_direction IS the 1-D dag_walk; CellVisit immutable + UpstreamState accumulator; 2-D wavefront is a DIFFERENT (level-batched) shape by design.
- [Issue #196 Phase G Step 2.5 — further-collapse feasibility](issue_196_phase_g_step2_5_further_collapse.md) — which sweep-body collapses are clean Pattern-2 wins vs principled Pattern-6 deferrals (fold-with-accumulator ≠ fold-with-parallel-reduce).
- [Issue #196 Phase G Step 2.5c — immutability audit](issue_196_phase_g_step2_5c_immutability_audit.md) — layered immutable-quantity enumeration (geometry/quadrature/σ_t/per-*) + the σ_t-immutable SweepCoefficientCache proposal.
