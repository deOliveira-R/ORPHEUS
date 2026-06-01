# SN sentinel (canary) regression harness

**Goal (the theory).** "For any regression, the set changes and we know
what." A regression ≈ a behaviour-changing code edit ≈ a *mutant*; "the
set changes" = a sentinel *kills* it; "we know" = the killed sentinel's
`cap()` label + the capability DAG (`sn_test_taxonomy.md`) localize the
blast radius. No small set catches 100 % — sentinels are a fast,
always-on **tripwire + localizer**, NOT proof. The full per-tier suite
stays the deep net. The **mutation score** quantifies how close to "any
regression" the set gets (the honest sensitivity number).

Built ON TOP of the capability taxonomy: the set covers every
capability NODE; a flip localizes via the DAG.

## The marker

`@pytest.mark.sentinel` (registered in `pyproject.toml`). Each sentinel
ALSO keeps its existing `cap()` / `l0` / `l1` / `catches()` marks — the
sentinel mark is additive. Selection: `pytest -m sentinel`.

## CRITICAL: run sentinels WITHOUT `-O`

The canonical ORPHEUS test invocation is `python -O -m pytest`
(`feedback_default_test_mode_is_optimize`). **`-O` strips bare `assert`
statements** — several sentinels (e.g. `test_spherical_alpha_dome_non_negative`,
`test_homogeneous_exact`) assert via bare `assert`, which become NO-OPs
under `-O`: a tripwire that cannot trip. Therefore the sentinel gate
MUST run **without `-O`**:

```sh
env PYTHONPATH=<worktree> <venv>/bin/python -m pytest tests/sn -m sentinel -q -p no:cacheprovider --timeout=120
```

(`np.testing.assert_*` calls are function calls and DO fire under `-O`;
bare `assert` does NOT. The sentinel set mixes both, so `-O` is unsafe
for it.) Measured wall-clock: **15 sentinels, 4.3 s** (no `-O`), all
green.

## Per-capability sentinel census

One sentinel per capability NODE (cheapest sharply-sensitive test;
slow convergence studies excluded). Corners probed per
`numerical-bug-signatures` + `vv-principles` ("bugs cluster at corners").

| cap node | sentinel | cost | catches / pillar | corner probed |
|---|---|---|---|---|
| primitives | `test_quadrature::test_per_ordinate_flat_flux_consistency[SPH,CYL]` | — | ERR-006/007 | curvilinear redistribution per-ordinate |
| primitives | `test_quadrature::test_scattering_source_magnitude` | — | ERR-002 | asymmetric 2G SigS convention |
| primitives | `test_quadrature::test_spherical_alpha_dome_non_negative` | — | Sig5 (α-dome) | α positivity (NaN precursor) |
| primitives | `test_quadrature::test_lebedev_weights_sum_to_4pi` | — | ERR-004 precondition | quadrature normalization (4π vs 2) |
| operators | `test_scattering_operator::TestBitIdenticalExtractionP0::test_add_iso_source_matches_reference` | 0.46 s | ERR-002-class (apply layer) | 2G scatter source apply vs reference |
| sweep_core | `test_sweep_vs_apply_consistency::test_carlson_seed_helper_is_linear_in_Q_bar` | 0.20 s | ERR-026 | matvec≡sweep curvilinear seed linearity |
| sweep_slab | `test_dd_recurrence::test_dd_per_cell_recurrence_matches_symbolic_derivation` | 0.53 s | ERR-025 | slab DD recurrence vs symbolic |
| sweep_curvilinear | `test_streaming_equilibrium_curvilinear::test_sentinel_sphere_streaming_equilibrium` | 0.37 s | ERR-026/048 | r=0 pole + WDD curvilinear fixed point |
| sweep_cartesian_2d | `test_2d_octant_sweep_equivalence::test_2d_octant_sweep_closed_form_anchor` | 0.81 s | closed-form (L1) | 2-D octant + reflective + 2G het |
| solve | `test_fixed_source_g1::…::test_uniform_source_converges_to_q_over_sigma_t[cylinder]` | 0.27 s | Signature 4 | curvilinear Q/Σ_t fixed-source equilibrium |
| eigenvalue | `test_keff_slab::test_homogeneous_exact[sn_slab_2eg_1rg]` | 1.16 s | ERR-001-class | NON-DEGENERATE solver eigenvalue (2G) |
| verification_analytical | `test_kinf_homogeneous::test_sentinel_kinf_slab_2g_krylov` | 1.66 s | eigenvalue PILLAR (closed-form) | 2G closed-form k_inf cross-check |
| verification_mms | `test_mms::test_sn_mms_manufactured_source_vanishes_at_zero_material` | 0.21 s | — | MMS source structure (slab) |
| verification_mms | `test_mms_aniso::test_sn_p1_aniso_mms_source_degrades_to_p0` | 0.21 s | — | P1 anisotropic angular-coupling corner |

15 node-IDs total (`per_ordinate` parametrizes ×2 coords). Every
capability node has ≥1 sentinel; eigenvalue claims are pillared on a
CLOSED-FORM reference (MMS does not prove eigenvalues per
`vv-principles`).

## Mutation-validation (Phase S2) — the sensitivity number

Tool: cosmic-ray (`tests/_mutation/`, Phase S0 verdict). Per-capability
recipe: mutate ONE tier's module, run the sentinel set, `dump` →
classify killed/survived by `definition_name`.

- `diamond.py` (sweep_core/slab/2d cluster) sentinel score: see
  `tests/_mutation/diamond_sentinels.toml` run. The FULL per-tier suite
  scores 99.7 % on diamond.py (Phase S0); the sentinel subset is a
  fraction of that — the honest sensitivity, reported as a number, not
  a "catches any regression" claim.

## Known sentinel coverage notes

- The slab-DD sentinel asserts only `cell_average_flux`, not
  `outgoing_spatial_flux`. The `psi_spat_out` spatial-closure sign flip
  in `update()` is instead caught by the 2-D closed-form anchor
  sentinel (verified: injecting `2·ψ_avg + ψ_in` → anchor fails with
  NaN). The tripwire trips; the localizer points at the sweep cluster.
- Full heterogeneous keff solves (`test_heterogeneous_absolute_keff`,
  ERR-025) are too slow on this branch (>180 s, does not converge) —
  EXCLUDED from sentinels; noted as a gap. The eigenvalue node is
  covered by the cheaper 2G homogeneous keff + the closed-form k_inf.

## Wiring recommendation (pre-commit / CI)

`pytest tests/sn -m sentinel` (no `-O`) is the seconds-fast always-on
gate. Wire as a pre-commit hook and/or the first CI job — if a sentinel
flips, its `cap()` + the DAG name the blast radius before the deep
per-tier suite runs. The deep net (per-capability tiers) runs after.
