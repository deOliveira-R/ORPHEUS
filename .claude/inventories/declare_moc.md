# MoC theory page — `implements` declaration candidates

Source page: `docs/theory/methods/method_of_characteristics.rst`
Derived from the tree at `HEAD = 9ba35a8f` and re-verified at `HEAD = 58e46c6f`
(2026-08-18) — the tree moved mid-dispatch. `[M]` the diff touches this page by
exactly ONE line (`:779`, an error-catalog path in prose) and touches
`orpheus/moc/` **not at all**, so every `file:line` below is current; all 21
`:label:` positions and all 5 `vv-status` positions are unchanged.
NOT derived from `.claude/inventories/implements_declaration_inventory.md`
(refuted on this population).

## Standing facts measured before any verdict below

- ⚠ **The brief's method step 1 does not apply to this page.** The MoC page
  carries **zero** `.. (vv-status rationale) ...` comments. It carries **5**
  bare status stamps of the form `.. vv-status: <label> documented`
  (`transport-equation`:201, `characteristic-ray-parametrization`:214,
  `bickley-integral`:369, `ray-perpendicular-coordinate`:447,
  `moc-homogeneous-weight-degeneracy`:1051) and **none of the 11 equations in
  this brief has one**. `[M]` `grep -n "vv-status" docs/theory/methods/method_of_characteristics.rst`
  → 5 hits. The rationale form exists on other pages
  (`methods/index.rst`, `methods/collision_probability.rst`), so the absence
  here is a property of this page, not of my search. Every `rationale
  comment` field below therefore reads "none", and the verdicts rest on
  steps 2–4 (equation + prose + symbol) instead.
- ⚠ **Nothing is standing down.** `[M]` all 11 equation nodes currently have
  **0 incoming and 0 outgoing `implements` edges** — the name-token heuristic
  produced no guess for any of them (their labels share no token with any
  symbol name). So an incomplete list here costs completeness of the ledger
  but does not orphan an existing guess.
- The MoC package is 1242 lines across 5 files and was read in full, so the
  `orpheus/moc/` half of every enumeration is exhaustive by construction. The
  non-`moc` half came from targeted greps, named per equation.
- Two implementers recur and are worth stating once:
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — the production transport
    sweep (`orpheus/moc/core.py:121-336`).
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` — an **independent**
    per-characteristic sweep used by the L1 MMS gate
    (`orpheus/derivations/continuous/mms/moc.py:249-368`). It re-implements the
    same segment kernel with a manufactured per-segment source; it is imported
    by `tests/moc/test_mms.py`, which carries 3 of the 11 labels.
  - `orpheus.derivations.discrete.moc.equations.*` — SymPy derivations whose
    module docstring calls them "the derivation source-of-truth for the
    documented equations".

---
## attenuation  (41 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none (see standing facts). The nearest
  authored guidance is the page's own `.. note::` at :266-273 —
  *"For small :math:`\tau < 10^{-10}`, the exponential is computed via Taylor
  expansion :math:`1 - e^{-\tau} \approx \tau(1 - \tau/2)` ... This is
  implemented in the inner loop of :meth:`MOCSolver.solve_fixed_source`."*
  That sentence names the implementer outright.
- **what the equation says**: the exact flat-source solution of the
  characteristic ODE across one segment —
  `ψ_out = ψ_in·e^{-τ} + (Q/Σ_t)(1 - e^{-τ})`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:214-220`
    (forward) and `:238-244` (backward) — computes `one_minus_exp` with the
    documented small-τ Taylor guard and then `psi[g] -= dpsi`, which is
    `ψ_in - (ψ_in - Q/Σ)(1-e^{-τ}) = ψ_in e^{-τ} + (Q/Σ)(1-e^{-τ})` exactly.
    The march is spelled through Δψ, but the marched quantity IS ψ_out.
  - `orpheus.derivations.discrete.moc.equations.derive_bar_psi` —
    `orpheus/derivations/discrete/moc/equations.py:33` — builds the attenuation
    solution as a SymPy expression
    (`psi_in*exp(-Sigma_t*s) + (Q/Sigma_t)*(1 - exp(-Sigma_t*s))`) and then
    integrates and `assert`s against the closed form. This is the symbolic
    statement of record for the same formula.
- **considered and EXCLUDED**:
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` — it computes
    `one_minus_exp` and `dpsi` but **never forms ψ_out**: each segment's
    `psi_in` is re-read from the reference solution
    (`mms/moc.py:310-314, 340-344`), so no attenuation march exists there.
    It implements `delta-psi` and `optical-thickness`, not this.
- ⛔ **CORRECTED IN PLACE.** I first wrote that the only other `exp(-tau)` hits
  were "unrelated Peierls/chord-oracle kernels that attenuate but **carry no
  flat-source term**". That is FALSE. `[M]`
  `grep -rn "/ Sigma_t) \* (1 - sp.exp" orpheus/ --include="*.py"` returns **12
  sites** under
  `orpheus/derivations/continuous/trajectory_resolvent/origins/specular/greens_function*.py`
  (slab, slab-asymmetric, cylinder, annulus, hollow sphere — e.g.
  `greens_function_annulus.py:244`, `greens_function_slab.py:220`) that spell
  exactly `(q/Σ_t)(1 - e^{-τ})`, the flat-source half of this equation. They are
  the periodic-trajectory Green's-function derivations and belong to
  `docs/theory/references/trajectory_resolvent.rst`, which carries its own
  labels; none of them is a MoC segment march (no `ψ_in e^{-τ}` carried across
  a segment sequence). I recommend NOT declaring them here — but the reason is
  a scoping ruling, not their absence.
- **confidence**: high on the two listed implementers (the page names the first
  outright, and the algebra is an identity, not a resemblance). What would change
  the enumeration: a ruling that the flat-source uncollided term counts wherever
  it appears, which would add the 12 trajectory-resolvent sites above.

---

## optical-thickness  (41 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none.
- **what the equation says**: the 3-D optical thickness of a 2-D segment,
  `τ_{g,p} = Σ_{t,i,g}·ℓ/sinθ_p` — i.e. the 2-D chord `ℓ` lifted to 3-D by the
  polar projection `1/sinθ_p`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` —
    `orpheus/moc/core.py:214` and `:238` —
    `tau = self.sig_t[rid, g] * seg.length / sin_p`.
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:324` and `:353` —
    `tau = sig_t * seg.length / sin_p`, character-for-character the same law on
    an independent sweep.
- **considered and EXCLUDED (marginal — your call)**:
  - `orpheus.derivations.discrete.moc.equations.derive_bar_psi` —
    `equations.py:30` writes `tau = Sigma_t * L`. That is an optical thickness,
    but with the 3-D path length `L` already substituted; the derivation has no
    polar angle at all, so it does not carry the `ℓ/sinθ_p` projection that is
    this equation's distinguishing content (and is the factor ERR-019 lost).
    Declare it only if you want "any τ" rather than "the MoC polar-projected τ".
- **confidence**: high for the two listed. What would change it: a decision that
  `τ = Σ_t L` counts, which would add `derive_bar_psi`.

---
## scalar-flux-integral  (41 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The nearest authored guidance is the
  claiming test module's own comment, `tests/moc/test_verification.py:42-47`:
  *"labels directly exercised by the L0 single-track tests (attenuation,
  optical-thickness) and by every MOC solve through the angular integral in the
  flux update (scalar-flux-integral)"* — i.e. the page's own author places this
  label on the flux update, not on a separate integrator.
- **what the equation says**: the scalar flux in an FSR is the segment-averaged
  angular flux integrated over all directions and over the region area, divided
  by the area: `φ_i = (1/A_i)∫_{4π}∫_{A_i} ψ̄ dA dΩ`. It is the *continuous
  definition*; `boyd-eq-45` is its quadrature.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:205`
    (the angular measure `4π·ω_a·ω_p·t_s·sinθ_p`, i.e. `dΩ` × the ray's share of
    `dA`), `:221`/`:245` (accumulation over segments = `∫_{A_i}`), and
    `:258-265` (the `/A_i` and the equilibrium term). These four lines ARE the
    double integral; no other object in the tree evaluates it.
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:305` (identical angular measure),
    `:331`/`:363` (accumulation), `:366-368` (`phi = phi_ref_avg +
    delta_phi/(areas*sig_t)`). Its docstring writes the reconstruction out.
    Partial by design: the equilibrium half is taken analytically, the transport
    half comes from the same integral.
- **considered and EXCLUDED**:
  - `orpheus.derivations.discrete.moc.equations.derive_scalar_flux_weight` —
    its docstring states this exact integral and derives its discretisation, but
    the **body is `print` statements only** (`equations.py:76-106`): it computes
    nothing, so no coverage run can ever corroborate a claim routed through it.
    See the note under `boyd-eq-45` — declare it only if you want the derivation
    of record on the chain.
  - `orpheus.moc.solver.solve_moc` post-processing (`solver.py:181-189`) — that
    is `⟨φ⟩` per *material*, a different functional of an already-computed `φ_i`.
- **confidence**: high. What would change it: a ruling that a continuous
  definition may not be `implements`-linked to its own quadrature, which would
  flip this to `NOTHING:definition` and move all 41 claims onto `boyd-eq-45`.

---

## ray-circle  (22 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The claiming test module names the
  gate instead (`tests/moc/test_ray_tracing.py:22-24`): *"ray-circle — directly
  tested by test_ray_circle_{hit,miss,tangent,chord_length} (4 dedicated L0
  tests for the quadratic formula)"*.
- **what the equation says**: a ray `(x_0,y_0) + s(cosφ, sinφ)` meets a circle
  of radius `R_k` centred at `(c_x,c_y)` where `s² + 2bs + c = 0`, with
  `b = (x_0-c_x)cosφ + (y_0-c_y)sinφ`, `c = (x_0-c_x)² + (y_0-c_y)² - R_k²`,
  discriminant `Δ = b² - c`, roots `s_{1,2} = -b ∓ √Δ`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.geometry._ray_circle_intersections` —
    `orpheus/moc/geometry.py:74-95` — the whole equation, term by term:
    `b = dx*cos_phi + dy*sin_phi` (`:89`), `c = dx*dx + dy*dy - radius*radius`
    (`:90`), `disc = b*b - c` (`:91`), `return [-b - sqrt_disc, -b + sqrt_disc]`
    (`:95`), including the `Δ < 0` miss branch (`:92-93`).
- **considered and EXCLUDED**:
  - `orpheus.moc.geometry._trace_single_ray` (`:192-196`) — it *calls* the above
    and applies the page's "only intersections within the cell are retained"
    filter. It implements the segment-construction prose, not the quadratic.
    Marginal; include only if you want the retention rule on this label.
- ⛔ **CORRECTED IN PLACE (my own first measurement was over-claimed).** I first
  wrote "the enumeration is exhaustive: `grep` returns `moc/geometry.py` and
  nothing else in the package". That grep carried `| grep -v derivations`, so
  its denominator was `orpheus/` **minus** `orpheus/derivations/` — and the
  excluded half is not empty. `[M]` re-run without the filter,
  `grep -rn "sqrt_disc" orpheus/ --include="*.py"` returns **~40 further hits**
  across `derivations/continuous/peierls_nystrom/geometry.py`,
  `derivations/continuous/trajectory_resolvent/{chord_oracle,greens_function}.py`
  and `derivations/common/quadrature_recipes.py`. Those solve the SAME quadratic
  in the polar/chord chart (`ρ = -r cos ω ± √Δ`, i.e. `b = r cos ω`,
  `c = r² - R²`) rather than the Cartesian `(x_0,y_0,c_x,c_y)` chart this
  equation is written in, and they belong to the CP / Peierls / trajectory-
  resolvent pages, which carry their own chord labels. I still recommend
  declaring **only** `_ray_circle_intersections`: declaring the chord solvers
  here would make every Peierls/CP chord test adjudicable against the MoC
  pin-cell tracer (the same wrong-granularity error as declaring
  `roots_of_unity` under `azimuthal-angles`).
- **confidence**: high that `_ray_circle_intersections` is the implementer and
  the only one in the MoC ray tracer. Medium that it is the only one you want
  declared — that turns on whether your ontology admits one label per chart or
  one per geometric fact, which is a ruling, not a measurement.

---

## effective-spacing  (22 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Claiming module's comment
  (`tests/moc/test_ray_tracing.py:19-21`): *"effective-spacing — exercised by
  test_volume_conservation which calls MOCMesh.effective_spacing(a_idx) and
  asserts that length * ts * omega_a sums to the geometric region area"*.
- **what the equation says**: rays are laid at `t_k = t_min + (k+½)t_s^eff`
  where `n_rays = ⌈(t_max - t_min)/t_s⌉` and `t_s^eff = (t_max - t_min)/n_rays`
  — the requested spacing rounded UP in count so the rays exactly tile the
  cell width at each azimuthal angle with no gap.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.geometry.MOCMesh._generate_tracks` —
    `orpheus/moc/geometry.py:333-334` is the equation verbatim
    (`n_rays = max(1, int(np.ceil((t_max - t_min) / ts)))`,
    `effective_ts = (t_max - t_min) / n_rays`), `:338` places the rays at the
    midpoints, and `:366-368` stores the per-angle value the solver later reads.
    This is the only site that computes the quantity.
- **considered and EXCLUDED (marginal — your call)**:
  - `orpheus.moc.geometry.MOCMesh.effective_spacing` (`:444-446`) — a pure
    read of `self._effective_spacing[azi_index]`. It is the *public surface*
    the claiming test calls, so declaring it would make the coverage join
    obvious; but it computes nothing. Prefer `_generate_tracks`, which any test
    constructing a `MOCMesh` also executes.
- ⚠ **Side finding worth a look while you are in here.** The value STORED at
  `:366-368` is `(t_max - t_min) / max(1, len(track_indices))` — divided by the
  number of tracks **kept**, after degenerate rays are dropped at `:349-350` —
  while the value used to PLACE the rays at `:334, :338` is
  `(t_max - t_min) / n_rays`. Two different denominators for one documented
  symbol; the weight in `boyd-eq-45` uses the first, the geometry uses the
  second. `[M]` on the shipped default configuration (standard PWR
  Wigner-Seitz pin cell, `n_azi=16`, `n_polar=3`, `ray_spacing=0.05`,
  reflective) they are **bit-identical on 16 of 16 azimuthal indices**, because
  no ray is dropped (`kept == n_rays` everywhere). So the divergence is latent,
  not active — but the equation as documented has one denominator and the code
  has two, and only one of them is what `t_s^eff` is defined to be.
- **confidence**: high on the implementer. What would change it: nothing about
  the symbol; the side finding is a separate (possible) defect, not a doubt
  about the declaration.

---

## pitch-recovery  (22 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. Claiming module's comment
  (`tests/moc/test_ray_tracing.py:25-26`): *"pitch-recovery — directly tested by
  test_moc_mesh_pitch_recovery (asserts pitch = r_cell * sqrt(pi) on factory
  output)"*.
- **what the equation says**: the INVERSE Wigner-Seitz transform — recover the
  square lattice pitch from the equal-area cylinder radius,
  `pitch = r_cell·√π = mesh.edges[-1]·√π`. (The forward direction,
  `r_cell = pitch/√π`, is the separate label `moc-wigner-seitz`.)
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.geometry.MOCMesh.__init__` — `orpheus/moc/geometry.py:272` —
    `self.pitch: float = mesh.edges[-1] * np.sqrt(np.pi)`. Sole site.
    (If a method-level `__init__` target is awkward for your tooling,
    `orpheus.moc.geometry.MOCMesh` resolves as a `class` node and is the same
    claim one level up.)
- **considered and EXCLUDED**:
  - `orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
    (`:438`, `r_cell = pitch / np.sqrt(np.pi)`), `orpheus.mc.solver` (`:92`) and
    `orpheus.derivations.continuous.mms.moc.build_moc_mesh` (`mms/moc.py:409`,
    `ws_r = P / np.sqrt(np.pi)`) — all three are the FORWARD transform
    (`moc-wigner-seitz`), not the recovery.
  - `orpheus.plotting.plot_moc_rays` — its docstring quotes
    `pitch = R_outer * sqrt(pi)` but the body reads `geom.pitch` (`plotting.py:132`).
    A consumer, not an implementer.
- **confidence**: high. `[M]` `grep -rn "sqrt(np.pi)" orpheus/ --include="*.py"`
  returns 5 code sites; 1 is this, 3 are the forward transform, 1 is an unrelated
  Hermite normalisation in `numerics/generating_measure.py:468`.

---
## azimuthal-angles  (14 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. But the page carries an authored
  `.. note::` at `:334-340` that is the load-bearing guidance: *"**The direction
  components are generated from the group action, not by evaluating trig on the
  angles** (issue #325, migrated 2026-08-08). The midpoint grid is the rational
  family :math:`\varphi_m = 2\pi(2m-1)/(4N_\varphi)` — the upper half of the
  STAGGERED periodic trapezoid on :math:`2N_\varphi` points."*
- **what the equation says**: `N_φ` uniformly spaced MIDPOINT azimuthal angles
  on `[0,π)`, `φ_m = π/(2N_φ) + (m-1)π/N_φ`, with equal weights `1/N_φ`
  summing to 1.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.quadrature.MOCQuadrature.create` —
    `orpheus/moc/quadrature.py:147-150`:
    `odd_k = 2*np.arange(n_azi) + 1`; `phi = np.pi*odd_k/(2.0*n_azi)`
    (= `π(2m-1)/(2N)` = the documented `π/(2N) + (m-1)π/N`, identically);
    `omega_azi = np.full(n_azi, 1.0/n_azi)`; and the primary representation
    `cos_phi, sin_phi = roots_of_unity(odd_k, 4*n_azi)`, the same grid as exact
    points. Sole construction site in the tree.
    (`orpheus.moc.quadrature.MOCQuadrature` resolves as a `class` node and its
    docstring states the grid, if you prefer the class-level target.)
- **considered and EXCLUDED**:
  - `orpheus.numerics.roots_of_unity.roots_of_unity` — the shared exact-point
    primitive. It implements the roots of unity, not this grid; it is consumed
    by `numerics/symmetry.py`, `numerics/quadrature/rules_circle.py` and an SN
    derivation as well, so declaring it here would make every SN/symmetry test
    adjudicable against the MoC azimuthal equation. Wrong granularity.
  - `orpheus.numerics.quadrature.rules_circle.periodic_trapezoid` — the generic
    staggered rule the page's note compares against. `[M]` MoC does **not**
    consume it (`grep -rn "roots_of_unity" orpheus/`: `moc/quadrature.py` calls
    the primitive directly).
- ⚠ **Side finding.** `orpheus.derivations.continuous.mms.moc.mms_sweep`
  re-derives the directions by **evaluating trig on the chart**
  (`mms/moc.py:279-280`, `np.cos(quad.phi[a_idx])` / `np.sin(...)`) instead of
  reading the exact `quad.cos_phi` / `quad.sin_phi` that `MOCMesh._generate_tracks`
  uses (`moc/geometry.py:326-327`, explicitly *"never trig on the chart"*).
  `[M]` the two disagree bitwise on **7 of 16 / 15 of 32 / 25 of 64** components
  at `n_azi = 8/16/32` (max abs difference `3.3e-16`), and the trig chart does
  NOT carry the exact azimuthal-mirror identity `cos_phi[n-1-k] == -cos_phi[k]`
  that #325 established (`True` on the points, `False` on the chart, all three
  orders). So the MMS reference sweeps along directions that are not exactly the
  directions its own tracks were traced with. Looks like a residual #325 site.
- **confidence**: high on the declaration; the side finding is independent of it.

---

## characteristic-ode  (4 claims)

- **verdict**: DECLARABLE (with one judgement call, stated)
- **rationale comment on the page**: none — and note the page DID stamp the two
  neighbouring equations `transport-equation` (`:201`) and
  `characteristic-ray-parametrization` (`:214`) as `documented` while leaving
  this one unstamped, so the author distinguished it from those two.
- **what the equation says**: along a characteristic the transport PDE collapses
  to a first-order linear ODE in arc length,
  `dψ_g/ds + Σ_{t,g}(s)ψ_g(s) = Q_g(s)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:318-322` (forward) and
    `:347-351` (backward). This is the only site in the tree that **assembles
    the operator** `(d/ds + Σ_t)` and applies it: `streaming = cos_phi*∂_xφ +
    sin_phi*∂_yφ` is `Ω̂·∇φ` projected on the ray, and
    `q_seg = (1/4π)(φ + sinθ_p·streaming/Σ_t)` is exactly `Q/Σ_t` for
    `Q = Σ_tψ + dψ/ds` with `ψ = φ/4π` and `d/ds = sinθ_p(cosφ ∂_x + sinφ ∂_y)`.
    (The same lines are the page's `moc-mms-qext`; that label is this ODE
    evaluated on the manufactured `ψ_ref`.)
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:214-220`,
    `:238-244` — **judgement call**: it advances ψ across each segment by the
    ODE's exact solution rather than by forming `dψ/ds`. Declare it if you hold
    that "solves the ODE" implements the ODE; DROP it if you hold that the
    solution is owned by the separate `attenuation` label (which is on the same
    lines). ⚠ Dropping it refutes `tests/moc/test_moc.py`'s claim on this label,
    since that module drives `solve_moc` and never touches `mms_sweep`; keeping
    it makes both claiming modules adjudicable.
- **confidence**: high for `mms_sweep`, medium for `solve_fixed_source` — and
  what would change the second is a ruling from you on the sentence above, not
  more evidence from the tree.

---
## bar-psi  (4 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none. The page's own ERR-019 section
  (`:1036-1039`) supplies the equivalent guidance: the missing `sinθ_p` factor
  is *"the 2-D to 3-D projection factor arising from the segment-averaged
  angular flux formula :eq:`bar-psi`"* — i.e. the author attributes the `sinθ_p`
  in the sweep weight to THIS equation.
- **what the equation says**: the segment-averaged angular flux is the
  asymptotic value plus the boundary-flux correction divided by the optical
  thickness, `ψ̄ = Q_{i,g}/Σ_{t,i,g} + Δψ/τ` — exact under the flat-source
  approximation.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.discrete.moc.equations.derive_bar_psi` —
    `orpheus/derivations/discrete/moc/equations.py:26-55` — the explicit
    statement: integrates the attenuation solution over `[0, L]`, divides by
    `L`, forms `expected = Q/Sigma_t + Delta_psi/tau`, and
    `assert sp.simplify(bar_psi - expected) == 0`. The identity of record.
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:205`
    and `:262-265` — implements it in FUSED form (never materialised as a
    variable): the `sin_p` in `weight = 4π·ω_a·ω_p·ts·sin_p` is precisely the
    `Δψ/τ = Δψ·sinθ/(Σ_t ℓ)` term's residue after the `t_s ℓ_k` area factor,
    and the `4π·Q[i,g]` in the flux update is the `Q/Σ_t` term. A sign or factor
    error in `bar-psi` lands exactly here — which is the ERR-019 history.
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:305` (same fused `sin_p`) and
    `:366-368` (`phi = phi_ref_avg + delta_phi/(areas*sig_t)`, where
    `phi_ref_avg` carries the equilibrium term analytically). Its own docstring
    (`:256-270`) writes the split out.
- **confidence**: high for `derive_bar_psi`; medium-high for the two fused
  sites. What would change the latter: a project convention that `implements`
  requires the quantity to be materialised (it is not, in either sweep) — but
  the ERR-019 record is strong positive evidence that both sweeps are
  accountable to this identity, and declaring only the derivation would refute
  both claiming test modules, neither of which imports it.

---

## boyd-eq-45  (4 claims)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none, but the page's **Implementation**
  paragraph (`:743-757`) is explicit: *"In :meth:`MOCSolver.solve_fixed_source`,
  the accumulator ``delta_phi[i, g]`` collects the weighted :math:`\Delta\psi`
  contributions during the sweep. The weight per segment is:
  ``weight = 4.0 * np.pi * omega_a * omega_p * ts * sin_p``"* — a verbatim quote
  of the implementing line.
- **what the equation says**: the discretised flux update,
  `φ_{i,g} = (1/Σ_t)[4π Q_{i,g} + (1/A_i)Σ_{m,p,k} 4π ω_m^a ω_p t_s sinθ_p Δψ]`
  — Boyd 2014 Eq. 45 in ORPHEUS's weight-normalisation convention.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:205`
    (the weight, character-identical to the page's code block), `:221`/`:245`
    (the `Σ_{m,p,k}` accumulation into `delta_phi`), `:254-265` (the update
    `phi_new = (4π Q + delta_phi/A_i)/Σ_t`). The canonical implementer.
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:305` (identical weight),
    `:331`/`:363` (accumulation), `:366-368` (update, with the `4πQ/Σ_t`
    equilibrium term replaced by the analytic `⟨φ_ref⟩_i`). Implements the
    Δψ half of Eq. 45 exactly; this is what `tests/moc/test_mms.py` executes.
- **considered and EXCLUDED**:
  - `orpheus.derivations.discrete.moc.equations.derive_scalar_flux_weight` —
    `equations.py:58-106` derives this equation step by step and ends
    *"DERIVED (algebraic outline): weight = 4*pi*omega_a*omega_p*t_s*sin(theta_p)"*,
    but the **body is nothing but `print` calls** — no SymPy, no assertion, no
    value. It is the derivation of record, not an implementation, and a claim
    routed through it can never be corroborated by a coverage run (no test
    imports it). Declare it only if your ontology wants derivation-of-record
    edges; if so, it also belongs on `scalar-flux-integral`.
  - `orpheus.derivations.discrete.moc.equations.verify_homogeneous_consistency`
    — it does compute `phi_new = 4πQ/Σ_t` symbolically, but that is the
    `Δψ = 0` degenerate branch, which the page gives its OWN label
    (`moc-homogeneous-weight-degeneracy`, `:1046`). It belongs there, not here.
- **confidence**: high.

---

## delta-psi  (1 claim)

- **verdict**: DECLARABLE
- **rationale comment on the page**: none.
- **what the equation says**: the angular-flux drop across a segment,
  `Δψ_{g,p} = ψ_in - ψ_out = (ψ_in - Q_{i,g}/Σ_{t,i,g})(1 - e^{-τ})` — "the
  central building block of the MOC scalar flux update" (page `:263-265`).
- **implementers** (complete list, each verified to resolve):
  - `orpheus.moc.core.MOCSolver.solve_fixed_source` — `orpheus/moc/core.py:219`
    (forward) and `:243` (backward):
    `dpsi = (psi[g] - q_over_sigt[rid, g]) * one_minus_exp`.
  - `orpheus.derivations.continuous.mms.moc.mms_sweep` —
    `orpheus/derivations/continuous/mms/moc.py:330` and `:362`:
    `dpsi = (psi_in - q_seg) * one_minus_exp`. Independent sweep, same law.
  - `orpheus.derivations.discrete.moc.equations.derive_bar_psi` —
    `orpheus/derivations/discrete/moc/equations.py:40`:
    `Delta_psi = (psi_in - Q/Sigma_t) * (1 - sp.exp(-tau))`, the symbolic
    statement it then substitutes and proves.
- **confidence**: high. The enumeration is exhaustive over the whole package:
  `[M]` `grep -rn "dpsi\b|delta_phi\b" orpheus/ --include="*.py"`, run with NO
  directory filter, returns only `orpheus/moc/`,
  `orpheus/derivations/continuous/mms/moc.py`,
  `orpheus/derivations/discrete/moc/` — plus three lines in `orpheus/cp/solver.py`
  (`:812, :828, :841`) where `delta_phi` is CP's outer convergence metric, an
  unrelated quantity sharing the spelling.

---
## Declaration table (copy-paste)

Every target below was resolved against `.nexus/graph.db` with the brief's
script; all are `function` / `method` node types, none is a module.

| equation | claims | implementers |
|---|---|---|
| `attenuation` | 41 | `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.discrete.moc.equations.derive_bar_psi` |
| `optical-thickness` | 41 | `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.continuous.mms.moc.mms_sweep` |
| `scalar-flux-integral` | 41 | `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.continuous.mms.moc.mms_sweep` |
| `ray-circle` | 22 | `orpheus.moc.geometry._ray_circle_intersections` |
| `effective-spacing` | 22 | `orpheus.moc.geometry.MOCMesh._generate_tracks` |
| `pitch-recovery` | 22 | `orpheus.moc.geometry.MOCMesh.__init__` |
| `azimuthal-angles` | 14 | `orpheus.moc.quadrature.MOCQuadrature.create` |
| `characteristic-ode` | 4 | `orpheus.derivations.continuous.mms.moc.mms_sweep` **+ (your ruling)** `orpheus.moc.core.MOCSolver.solve_fixed_source` |
| `bar-psi` | 4 | `orpheus.derivations.discrete.moc.equations.derive_bar_psi`, `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.continuous.mms.moc.mms_sweep` |
| `boyd-eq-45` | 4 | `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.continuous.mms.moc.mms_sweep` |
| `delta-psi` | 1 | `orpheus.moc.core.MOCSolver.solve_fixed_source`, `orpheus.derivations.continuous.mms.moc.mms_sweep`, `orpheus.derivations.discrete.moc.equations.derive_bar_psi` |

---

## Summary

**All 11 equations are DECLARABLE — none is a `NOTHING:<kind>`.** Every one has
a computational site, and for 10 of the 11 the enumeration is settled by reading
rather than by judgement: the MoC package is 1242 lines across 5 files and was
read in full, and the only MoC kernel outside it is the independent MMS sweep
(`[M]` `grep -rn "dpsi\b|delta_phi\b|one_minus_exp" orpheus/ --include="*.py"`,
run with no directory filter, reaches `orpheus/moc/`,
`derivations/continuous/mms/moc.py`, `derivations/discrete/moc/equations.py`, and
three CP convergence-metric lines that share the spelling `delta_phi` and nothing
else).
The single equation I could **not settle without a ruling from you is
`characteristic-ode`**: `mms_sweep` certainly implements it (it is the only site
in the tree that assembles the operator `d/ds + Σ_t` and applies it, to build the
manufactured source), but whether `MOCSolver.solve_fixed_source` also implements
it turns on a labelling question no measurement decides — the solver advances ψ
by the ODE's *exact solution*, and that solution is the separate `attenuation`
label sitting on the same lines. Keeping the solver makes both claiming test
modules adjudicable; dropping it refutes `tests/moc/test_moc.py` on this label
for a reason that is about taste, not coverage. Two further scoping rulings are
flagged in place rather than decided: whether the 12 flat-source
`(q/Σ_t)(1-e^{-τ})` sites in the trajectory-resolvent Green's-function
derivations belong on `attenuation` (I say no — different theory page, own
labels), and whether the ~40 polar-chart chord-quadratic sites belong on
`ray-circle` (same answer, same reason). Finally, two possible defects surfaced
that are independent of the declarations and are worth their own look:
`MOCMesh` computes `t_s^eff` with **two different denominators** (`n_rays` for
ray placement, tracks-kept for the weight — `[M]` bit-identical on 16 of 16
azimuthal indices at the shipped default, so latent not active), and
`mms_sweep` re-derives its sweep directions by **evaluating trig on the angle
chart** (`mms/moc.py:279-280`) instead of reading the exact `#325` points its own
tracks were traced with (`[M]` bitwise-different on 15 of 32 components at
`n_azi=16`; the trig chart loses the exact azimuthal-mirror identity).
