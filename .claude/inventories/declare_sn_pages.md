# `implements` declaration findings — the 17 SN-page equation labels (138 coverage claims)

⚠ **Brief arithmetic note**: the dispatch said "16 equations / 130 claims"; its own table
lists **17** labels across 9 pages, and their claim counts sum to **138**. All 17 are
covered below. (`dd-slab-scalar` 4 + `dd-curvilinear-scalar` 13 + `dd-cylindrical-degenerate` 1
is the row most easily read as one.)

Derived from the tree at HEAD `a1c90aac` (2026-08-18), branch `main`.
**Does NOT use** `.claude/inventories/implements_declaration_inventory.md` (refuted
on this population).

Method per equation: (1) the page's own `.. (vv-status rationale)` comment, if any —
these are frequently the whole answer; (2) the equation + surrounding prose; (3) find
the symbol in the tree; (4) verify it resolves as a `function`/`method`/`class`/`data`
node in `.nexus/graph.db`.

Every symbol listed below was verified to resolve. Node-type is noted where it is
not `function`.

---

## `blelloch-1990-eq-1-5`  (19 claims) — `docs/theory/conventions/indexing_and_layout.rst:271`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none (the nearest, at `:319`, belongs to
  `sn-cell-flatten-roundtrip`, a different equation). The prose *is* the answer:
  *"In numpy this is three ops: ``cumprod(a) * (psi_0 + cumsum(b / cumprod(a)))``"*
  and *"The closed-form scan that :func:`~orpheus.sn.sweep.scan.ordinate_scan`
  evaluates is the **Blelloch §1.5 first-order linear-recurrence form**"*.
- **what the equation says**: the closed form of the affine recurrence
  `ψ[i+1] = a[i]ψ[i] + b[i]` — `ψ[n] = (∏a)·(ψ₀ + Σ b/∏a)`, i.e. the prefix-product
  factorisation that turns forward substitution into three vectorised numpy ops.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.sweep.scan.ordinate_scan` — `orpheus/sn/sweep/scan.py:82`
    (body at `:175-177`: `cumprod_a = np.cumprod(a, axis=0)` then
    `cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))`) — this IS the
    closed form, literally transcribed; the module docstring names Blelloch 1990
    §1.5 as its source.
- **deliberately EXCLUDED, with the structural reason** (so the next audit does
  not re-derive it):
  - `orpheus.sn.sweep.scan._pair_monoid_scan` (`scan.py:190`) computes the same
    `ψ[1…nx]` but its own docstring says *"rather than the Blelloch factored
    form"* — it is a division-free Hillis–Steele prefix scan of the pair monoid,
    the numerically-robust fallback taken when the closed form goes non-finite
    (ERR-057 / #222). It implements the *recurrence*, not *this equation*.
  - `ordinate_scan_transpose` (`scan.py:239`), `CumprodScan`, `ScanMarch`,
    `_OneDimScanWalk` (`orpheus/sn/loss_representation/__init__.py:1359 / 2203 /
    2815`) are all **consumers** — they call `ordinate_scan`; grep confirms no
    second `np.cumprod`/`np.cumsum` of this form anywhere in `orpheus/`.
- **confidence**: high. The page names the function; the function's body is the
  equation; the single-source claim is grep-verified (only one `cumprod`-based
  affine-scan body in the package).

---

## `hebert-3-432`  (7 claims) — `docs/theory/methods/sn/curvilinear_numerics.rst:265`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none on `hebert-3-432` itself. The nearest
  (`:310`) belongs to `hebert-3-433` and reads: *"Literature-transcribed derivation
  step: the Hébert Eq. (3.433) DD cell-balance for the μ=−1 starting-direction sweep
  — a verbatim reference definition en route to the Carlson seed. Its terminal
  result (the seed producing the correct flat-flux / cold-start solution) is tested
  downstream; the hebert-3-43X family's optional explicit wiring is tracked on Issue
  #194."*  ⭐ **The page's own `.. note::` at `:365-378` is the load-bearing
  authored knowledge**: *"The three labels :eq:`hebert-3-432-source`,
  :eq:`hebert-3-434`, :eq:`hebert-3-435` are also declared in the
  :mod:`~orpheus.sn.sweep.psi_half_angle_seed` module docstring (the canonical
  algebra-of-record) … the Sphinx page is the presentation layer for the equations
  the code module owns as source-of-truth."*  (Note it does NOT list `hebert-3-432`
  itself — that one is the continuous ODE, not a march step.)
- **what the equation says**: the continuous starting-direction transport ODE at
  `μ = −1`, where the angular-redistribution coefficient `(1−μ²)` vanishes so the
  balance decouples from the α-cascade: `−∂φ_{−1/2}/∂r + Σ(r)·φ_{−1/2}(r) =
  Σ_ℓ (2ℓ+1)/2 · Q_ℓ(r) · P_ℓ(−1)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
    — `orpheus/sn/operators/radial_characteristic.py:202` (class) — this IS the
    operator form of the equation: its own docstring, *"the banded radial transport
    operator :math:`\mu\,\partial_r + \sigma_t` acting on the starting-direction flux
    :math:`\psi_{1/2}` at the closed :math:`\mu = \pm 1` rays (Hébert §3.9.4)"*, and
    its `References` block cites *"Eqs. (3.432)–(3.435)"*. `.apply` is the LHS action
    `A_BB ψ½`; `.solve` is the exact direct inverse (the march).
  - `orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_forward_residual`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:245` — *"The forward ``(L+C)`` action
    :math:`A_{BB}\psi_{1/2} = (\mu\,\partial_r + \sigma_t)\psi_{1/2}` … the SINGLE
    SOURCE of the ψ½ forward residual"*, i.e. the discrete realisation of this
    equation's left-hand side, shared by `RadialCharacteristicOperator.apply`.
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:112` — the direct solver of this ODE
    (its own docstring: *"The direct solver of the starting-direction ODE along the
    level's own ray"*). It is listed here **as well as** under `hebert-3-434/435`
    because 3-434/435 are its discrete steps and 3-432 is the equation it solves.
- **adjoint-side siblings** (declare only if the project's convention is that a
  transpose implements its equation — I flag rather than assume):
  `radial_characteristic_forward_residual_transpose` (`:298`),
  `RadialCharacteristicOperator.apply_transpose` / `.solve_transpose`.
- **⚠ boundary note for the declarer**: the equation's *right-hand side*
  (the `Σ_ℓ (2ℓ+1)/2 Q_ℓ P_ℓ(−1)` fold) carries its **own** label
  `hebert-3-432-source` on the same page (`:283`) — **not in this batch**. Its
  implementers are
  `orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
  and `…radial_characteristic.RadialCharacteristicReconstruction` (whose docstring
  says *"The 1-D angular reconstruction ℛ (Hébert Eq. 3.432) SAMPLED at the closed
  rays"* — it cites 3.432 but states 3-432-source's formula). If `hebert-3-432` is
  declared to cover the whole ODE, those two belong under `hebert-3-432-source`,
  not here; do not double-declare.
- **confidence**: high for the three listed. Medium on the LHS/RHS boundary — the
  tree's docstrings cite "Eq. 3.432" for the reconstruction too, so if the intent
  is that `hebert-3-432` covers the full ODE including its source, add
  `RadialCharacteristicReconstruction` and `fold_moments_to_radial_characteristic`.
  What would change it: a ruling on whether `hebert-3-432-source` is a *specialisation*
  (own implementers) or a *sub-expression* (shares 3-432's).

---

## `hebert-3-434`  (4 claims) — `docs/theory/methods/sn/curvilinear_numerics.rst:345`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none directly attached; the page's `.. note::`
  (quoted in full under `hebert-3-432` above) explicitly declares this label to be
  owned by the `orpheus.sn.sweep.psi_half_angle_seed` module docstring as *"the
  canonical algebra-of-record"*, with the test wiring
  (`tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py`) tracked at Issue #194.
- **what the equation says**: the DD cell-average update of the inward
  starting-direction march — `φ̄_i = (Δr_i·Q̄_i + 2·φ̄_{i+1/2}) / (Δr_i·Σ_i + 2)`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:112` — the SOLVE. Body at `:177-178`:
    `denom = dr[k]*sigma_t[:,k] + 2.0` ; `phi_cell = (dr[k]*Q_bar[:,k] + 2.0*phi_face)/denom`.
    Verbatim. **The module docstring re-states this equation under
    `:label: hebert-3-434`** (`psi_half_angle_seed.py:21-28`).
  - `orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_residual_march`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:242` — the RESIDUAL (forward-direction)
    form of the same relation: `mᵢ = σᵢcᵢ + (2/Δrᵢ)(cᵢ − fᵢ)`, which is 3-434
    rearranged to `A·φ − q`. Its docstring names it *"the Hébert (3.434) residual"*.
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:187` — the exact reverse-mode adjoint
    of the 3-434 march (`Q_bar[:,k] = (dr[k]/denom)*c_bar`; same `denom`). Adjoint-side;
    include per the project's transpose convention.
- **also touched but NOT an implementer**:
  `radial_characteristic_forward_residual` / `…_transpose` assemble the two legs and
  the corner rows around `radial_characteristic_residual_march`; the 3-434 arithmetic
  is inside the callee, so declaring the callee is the honest single-source choice.
  (They ARE declared under `hebert-3-432`, whose operator statement they realise.)
- **confidence**: high. The page names the owning module; the module's own docstring
  carries the label; the arithmetic is a literal transcription.

---

## `hebert-3-435`  (4 claims) — `docs/theory/methods/sn/curvilinear_numerics.rst:355`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none directly attached; same page `.. note::` as
  above declares the label owned by the `psi_half_angle_seed` module docstring.
- **what the equation says**: the DD face chain — step inward from the cell average
  to the next face: `φ̄_{i−1/2} = 2·φ̄_i − φ̄_{i+1/2}`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:112`, body at `:180`:
    `phi_face = 2.0 * phi_cell - phi_face`, with the inline comment
    *"# Hébert (3.435): step to the next (downwind) face."*  The module docstring
    re-states the equation under `:label: hebert-3-435` (`:30-36`).
  - `orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_residual_march`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:242`, body: `f = 2.0*cells[:,i] - f`
    ("*with the DD face chain ``f ← 2·cᵢ − f``*").
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
    — `orpheus/transport/spatial/scheme.py:1408` (method) — the **`w`-generic**
    primitive `ψ_out = (ψ̄ − (1−w)ψ_in)/w`, of which 3-435 is the `w = ½` instance
    (*"DD's diamond mean ``2ψ̄ − ψ_in`` is the ``w=½`` case (byte-identical: ``÷½``
    is an exact power-of-2 ``×2``)"*). ⚠ **Generic**: it equally realises the WDD /
    LD closures, so declaring it against 3-435 makes 3-435 share a body with
    `wdd-face` and `dd-slab-scalar`. Correct but wide — flagged so the choice is
    deliberate.
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average_transpose`
    — `orpheus/transport/spatial/scheme.py:1424` (method) — the VJP of the same, and
    `psi_half_angle_seed.py:103-107` explicitly binds it as the ψ½ marches' face-chain
    transpose: *"The Hébert (3.435) face chain ``f ← 2·c − f`` IS the diamond
    (``w = ½``) affine chain, so every reverse march below pulls its chain cotangent
    pair through the ONE w-generic VJP primitive at ``w = _DD_W`` (#311)."*
  - `orpheus.transport.spatial.diamond._DD_W` — `orpheus/transport/spatial/diamond.py:102`
    (**`data` node**, `_DD_W: float = 0.5`) — the single-source constant that *makes*
    the generic chain the Hébert 3-435 chain. Declarable per the widened ontology;
    include only if the project wants the constant itself in the ledger.
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:187` — adjoint of the 3-434/3-435 pair
    (it calls the VJP above at `_DD_W`).
- **confidence**: high on the two `psi_half_angle_seed` marches (literal, labelled in
  the module docstring). Medium on the generic `scheme` primitives — they are
  unambiguously the body the chain runs through, but they are `w`-generic, so
  whether the ledger wants a shared body under three labels is a project ruling.

---

## `mm-weights`  (28 claims) — `docs/theory/methods/sn/curvilinear_one_group.rst:588`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none attached to this label (the file's 11
  rationale comments all sit ≥ 577 lines below, on other equations). The prose IS
  strongly directive: the section defines the sphere and cylinder cell-edge
  conventions, states the admissible range as *"predicate **P3**, enforced since
  Q5.5: the producer RAISES on :math:`\tau \notin [0, 1]`"*, and the tree's own
  producer docstring re-states the equation verbatim.
- **what the equation says**: the Morel–Montry angular closure weight —
  `τ_n = (μ_n − μ_{n−1/2}) / (μ_{n+1/2} − μ_{n−1/2})`, i.e. the barycentric
  coordinate of the ordinate between the two edges of its own angular cell
  (BMC 2010 Eq. 43 = Lathrop 2000 Eq. 23).
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
    — `orpheus/sn/sweep/pole_angular_closure.py:1051` — **the single production
    site**, and the only one. Body at `:1153-1159`:
    `dmu = mu_edge[m+1] - mu_edge[m]` ; `tau[m] = (mu[m] - mu_edge[m]) / dmu`.
    Docstring: *"`τ_m = (μ_m − μ_{m−1/2})/(μ_{m+1/2} − μ_{m−1/2})` (Bailey–Morel–Chang
    2010 Eq. 43 = Lathrop 2000 Eq. 23) — the predicate **P2**"*, and *"⭐ **This is
    the whole τ — there is no 'raw' and no 'clamped'.**"*
    ⭐ **ONE generic body serves both geometries** since Q5.6.4 (2026-08-11); before
    that the sphere and cylinder each re-derived edges inline here. `[M]` I grepped
    `orpheus/` for any second `(mu - edge)/(edge - edge)` construction: none.
    `orpheus/derivations/discrete/sn/angular_differencing.py` had a
    `morel_montry_weights` twin — **RETIRED 2026-08-12** (`:184`, `:284`), body was
    a pass-through to this function.
  - `orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`
    — `orpheus/sn/sweep/pole_angular_closure.py:895` — produces the `μ_{n±1/2}` this
    equation is a ratio of, and the page's two paragraphs immediately BELOW the
    equation ("**Spherical cell edges**" — `μ_{1/2} = −1`, `μ_{n+1/2} = μ_{n−1/2} + w_n`;
    "**Cylindrical cell edges**" — the ω-midpoint
    `η_{m+1/2} = sinθ·cos(½(ω_m + ω_{m+1}))`) are *its* two branches, transcribed at
    `:1001-1005` and `:1030-1042`. ⚠ Judgement call: strictly it implements the
    edge convention, not the τ ratio. Declare it here only if the label is read to
    cover the section it heads; otherwise it is a candidate for its own label.
  - `orpheus.sn.sweep.pole_angular_closure._assert_tau_within_unit_interval`
    — `orpheus/sn/sweep/pole_angular_closure.py:863` — enforces the equation's stated
    admissible range `τ ∈ [0,1]` (predicate **P3**), a real `raise` (not a bare
    `assert`, so it survives `python -O`). ⚠ Same judgement call: it implements the
    *admissibility clause* of the label's prose, not the formula.
- **deliberately EXCLUDED, with the structural reason**:
  - `MorelMontryAngularSweep` (`:~1380`) and `PoleAngularClosureBase.tau_per_ordinate`
    (method node) **own and re-lay-out** τ — the class calls the producer once in
    `__init__` (`:1444`) and the property scatters per-level → per-ordinate. Storage
    and view, no arithmetic of this equation.
  - `orpheus/derivations/discrete/sn/balance.py:252` spells
    `psi_a_out_wdd = (psi_avg - (1 - tau)*psi_angle_in)/tau` — that is `wdd-face`
    (the closure SOLVED for the angular face flux), with `tau` a free SymPy symbol;
    it never computes τ from edges.
  - `orpheus.transport.spatial.diamond.DiamondDifference` reads `visit.tau`
    (`diamond.py:227`) — a consumer of the stamped value.
- **INVERSE-form sites** (they solve the same relation for the *edges* given τ — the
  ν-closure; declare only if the ledger admits inverse realisations):
  - `orpheus.derivations.discrete.sn.angular_differencing.morel_montry_beta`
    — `:334`, body `implied[m+1] = (mu[m] - (1.0 - tau[m])*implied[m]) / tau[m]`.
  - `orpheus.derivations.discrete.sn.angular_differencing.nu_closure_residual`
    — `:428`, same recurrence; the page's own principled condemnation of the retired
    `[½,1]` absorber is measured with it.
- **confidence**: high that `morel_montry_tau_per_level` is the complete production
  answer (single-source is asserted in the module and grep-confirmed; the one
  historical twin is retired with a dated note). Medium on the three "supporting
  clause" symbols (`angular_cell_edges_per_level`, `_assert_tau_within_unit_interval`,
  and the two inverse-form diagnostics) — those turn on how wide the label is read,
  which is a project ruling, not a tree fact. What would change it: a convention on
  whether a label covers only its displayed formula or the section's stated
  conventions/admissibility too.

---

## `dd-null-sawtooth`  (2 claims) — `docs/theory/methods/sn/cartesian_multid.rst:4624`

- **verdict**: DECLARABLE
- **rationale comment on the page** (`:4628-4636`, verbatim):
  > *"(vv-status rationale) Structural identity: the null-space specialisation of
  > the already-verified multi-D DD closure (dd-cartesian-2d) at psi_c = 0 — an
  > algebraic rearrangement, not a new solver claim. Its content is that every
  > null vector has this SHAPE, which is asserted end-to-end by the foundation
  > suite tests/sn/operators/test_loss_kernel_gauge.py: the constructed basis is
  > annihilated by the PRODUCTION matvec
  > (test_EVERY_basis_vector_is_annihilated_by_the_production_matvec) and its
  > dimension equals a dense SVD of the assembled operator
  > (test_the_dimension_matches_a_DENSE_SVD_of_the_assembled_operator)."*
- **what the equation says**: every null vector of the multi-D DD operator, restricted
  to axis `a`'s faces, alternates in the face index and factorises —
  `ψⁿ_a(k, i_⊥) = (−1)^k · φⁿ_a(i_⊥)` — which is why the kernel is a *trace* object.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.operators.loss_kernel_gauge._block_support`
    — `orpheus/sn/operators/loss_kernel_gauge.py:852` — evaluates the `(−1)^k` factor
    at each face of the block: `k = 0 if sign < 0 else spatial[axis]` then
    `1.0 if k % 2 == 0 else -1.0` (`:907-909`), with the inline comment
    *"# k = 0 at the min face, k = n_axis at the max face."*
  - `orpheus.sn.operators.loss_kernel_gauge._FacePlacement`
    — `orpheus/sn/operators/loss_kernel_gauge.py:839` (**class**) — carries the value
    as a documented field: `#: ``(-1)^k`` at this face: ``+1`` at the min face
    (``k = 0``), ``(-1)^{n_a}`` at the max face (``k = n_a``) — the sawtooth
    evaluated at the two ends.` The equation is literally this field's docstring.
  - `orpheus.sn.operators.loss_kernel_gauge._build_block_table`
    — `orpheus/sn/operators/loss_kernel_gauge.py:943` — forms the product the
    equation states, `phi = place.sawtooth_sign * array * checker / (|μ_a| · area)`
    (`:990-991`), i.e. maps each generator's `φⁿ_a(i_⊥)` to the face DOFs through the
    sawtooth. Its own docstring names step 2 as
    *"each generator mapped to the block's DOFs through
    :math:`\varphi = Y \cdot (-1)^{\sum i} / (|\mu_a| A_a)`"*.
- **deliberately EXCLUDED, with the structural reason**:
  - `_transverse_factors` (`:921`) computes `(−1)^{Σ_{c≠a} i_c}` — the **transverse**
    checkerboard, which belongs to `dd-null-balance-combinatorial` (the substitution
    image, `:4673` rationale), NOT to the face-index sawtooth. Different index set.
  - `LossKernelGauge` / `LossKernelBasis` consume the table; no sawtooth arithmetic.
- **confidence**: high. The rationale names the owning suite and the owning object;
  the field docstring is a verbatim restatement; there is exactly one `(−1)^k`
  face-index construction in `orpheus/`.

---

## `dd-null-counting-law`  (4 claims) — `docs/theory/methods/sn/cartesian_multid.rst:4732`

- **verdict**: DECLARABLE
- **rationale comment on the page** (`:4741-4751`, verbatim):
  > *"(vv-status rationale) Structural identity: a combinatorial count derived
  > from dd-null-balance-combinatorial, carrying no solver claim. It is
  > nonetheless doubly gated by the foundation suite
  > tests/sn/operators/test_loss_kernel_gauge.py — **the law is evaluated without
  > building a vector by `predicted_kernel_dimension`** and compared BOTH against the
  > rank the construction's SVD finds
  > (test_the_dimension_matches_the_combinatorial_counting_law) and against a
  > dense SVD of the assembled operator
  > (test_the_dimension_matches_a_DENSE_SVD_of_the_assembled_operator), with the
  > two closed-form specialisations pinned separately
  > (test_the_counting_law_reproduces_the_two_closed_form_specialisations)."*
  ⭐ The rationale **names the implementing function outright** — this is the
  authored-knowledge case the method was designed for.
- **what the equation says**: `dim ker A = n_g · Σ_orbits Σ_{U⊆S, |U|≥2} κ(U)·Π_{c∉U} n_c
  + #{tangential trace DOFs}` — the mesh-independent, cross-section-independent
  combinatorial dimension of the DD loss kernel.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.sn.operators.loss_kernel_gauge.predicted_kernel_dimension`
    — `orpheus/sn/operators/loss_kernel_gauge.py:565` — named by the rationale.
    Body is the law term-for-term: the orbit loop (`for _orbit, active in
    _reflection_orbits(sn_mesh)`), the `|U| ≥ 2` subset loop
    (`for size in range(2, len(active)+1): for subset in combinations(...)`), the
    spectator product `Π_{c∉U} n_c`, and the final `* int(sn_mesh.ng)`.
    Docstring: *"``dim ker A`` from the counting law — **without building any
    vector**."*
  - `orpheus.sn.operators.loss_kernel_gauge._anova_dimension`
    — `orpheus/sn/operators/loss_kernel_gauge.py:533` — the `κ(U)` summand, whose
    own docstring re-displays the κ formula (the page shows κ in an **unlabelled**
    `.. math::` immediately above the counting law, so κ has no label of its own and
    its implementation belongs here). Its doctests pin
    `κ((5,7)) = 1` and `κ((3,4,5)) = 11 = n_a+n_b+n_c−1`, i.e. the page's two
    closed-form specialisations.
- **⚠ honest gap the declarer should know**: `predicted_kernel_dimension` implements
  the **first term only**. Its own docstring: *"Counts the **R** component only …
  the tangential component ``T`` contributes one dimension per tangential trace DOF
  and is not gauged."*  So the `+ #{tangential trace DOFs}` term of the displayed
  equation has **no implementer** in `orpheus/`. That is a genuine (documented)
  partial realisation, not a missed symbol.
- **deliberately EXCLUDED, with the structural reason**:
  - `_reflection_orbits` (`:614`) supplies the orbit partition the outer `Σ_orbits`
    runs over, and `gauge_freedom` (`:352`) supplies the early-out; both are inputs
    to the count, not the count.
  - `LossKernelGauge.dimension` (`:1085`) reports the dimension the SVD **built** —
    it is the counting law's independent CHECK (the whole point of the pairing), so
    declaring it here would collapse the two routes the gates exist to compare.
    ⛔ Do not declare it.
- **confidence**: high. The rationale names the function; the body is the law
  term-for-term; the one caveat (the tangential term) is stated by the code itself.

---

## ⚠ Shared premise for the four `dd-*` cell-update labels (read before declaring any)

`docs/theory/methods/sn/index.rst` presents `dd-slab-scalar`, `dd-curvilinear-scalar`
and `dd-cylindrical-degenerate` as **three branches of
`DiamondDifference.update`**, selected by `alpha_in is None` and `abs_mu < 1e-15`.

⛔ **That prose is STALE against the tree** (`[M]` HEAD `a1c90aac`). `update`
(`orpheus/transport/spatial/diamond.py:177`) now carries the docstring *"One body —
no geometry dispatch"*, and the branch collapse happened at Issue #196 Phase G Step
2.5, whose helper module says so explicitly
(`orpheus/transport/spatial/cell_balance.py:39-63`):

> *"Pre-Step-2.5 this lived in two helpers (`cell_balance_terms` for the
> non-degenerate curvilinear branch and `cell_balance_terms_degenerate` for the
> cylindrical pure-azimuthal case) plus an inlined slab recurrence inside
> `DiamondDifference`. Step 2.5 collapses all three into this single helper: •
> **Slab** … • **Curvilinear non-degenerate** … • **Cylindrical pure-azimuthal
> degenerate** …"*

⟹ **the three labels now share the same bodies**, and are distinguished only by the
*data* fed in (`slab_streaming`'s neutral curvature; `A_downstream = 0.0` for the
degenerate cylinder). Declaring them is therefore correct but *not* one-body-per-label.
The page's own three "mirroring ``_sweep_1d_*`` (the dissolved ``sweep.py``) lines
NNN–NNN" citations are dead addresses — `sweep.py` no longer exists.
**This is a doc-drift finding worth its own issue; it does not block the declarations.**

---

## `dd-slab-scalar`  (4 claims) — `docs/theory/methods/sn/index.rst:599`

- **verdict**: DECLARABLE
- **rationale comment on the page**: **none** — `docs/theory/methods/sn/index.rst`
  contains **zero** `.. (vv-status rationale)` comments (`[M]` `grep -c` = 0).
- **what the equation says**: the slab specialisation of the DD cell update — the
  diamond mean `ψ̄_i = ½(ψ_{i−1/2} + ψ_{i+1/2})` together with the face recurrence
  `ψ_{i+1/2} = [(2|μ|−ΔxΣ_t)/(2|μ|+ΔxΣ_t)]·ψ_{i−1/2} + 2QΔx/W/(2|μ|+ΔxΣ_t)`.
- **implementers — TIER 1 (the arithmetic, on the production path)**:
  - `orpheus.transport.spatial.cell_balance.cell_balance_terms`
    — `orpheus/transport/spatial/cell_balance.py:268` — the balance denominator/numerator,
    with the slab collapse spelled in its own inline comments:
    *"reduces to 2|μ|·1 for slab"*, *"slab: |μ|·2·ψ^s_in = 2|μ|·ψ^s_in"*.
  - `orpheus.transport.spatial.diamond.DiamondDifference.update`
    — `orpheus/transport/spatial/diamond.py:177` — `psi_avg = (source +
    terms.numer_upstream) / terms.denom`, then the `w = ½` outflow closure.
  - `orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`
    — `orpheus/transport/spatial/diamond.py:572` — the **production 1-D slab route**
    (`CumprodScan`/`ScanMarch` consume it): `a_attenuation = 2|μ|A_total/denom − 1`,
    which at slab (`A_total = 2`, `A_down = 1`, `dA_w = 0`, `V = Δx`) is *algebraically
    identical* to the page's `a_i` — verified by hand:
    `2|μ|·2/(2|μ|+Σ_tΔx) − 1 = (2|μ|−Σ_tΔx)/(2|μ|+Σ_tΔx)`. Its docstring states the
    collapse: *"For slab the curvature fields are neutral … the denominator collapses
    to the slab form 2|μ_n| + Σ_t V."*
  - `orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
    — `orpheus/transport/spatial/diamond.py:388` — the **dimension-generic Cartesian**
    kernel (`d = 1, 2, 3`); at `d = 1` its
    `ψ_avg = (Q + Σ_a 2g_a ψ_in_a)/(Σ_t + Σ_a 2g_a)` IS this equation with
    `g = |μ|/Δx`. (It is the SOLVE arm of the sweep-DAG walks.)
- **implementers — TIER 2 (`w`-generic primitives this equation instantiates at `w=½`)**:
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average`
    — `scheme.py:1396` — `ψ̄ = (1−w)ψ_in + w·ψ_out`; *"DD's ½(ψ_in + ψ_out) is the
    ``w=½`` special case"* — literally the equation's first line.
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
    — `scheme.py:1408` — `ψ_out = 2ψ̄ − ψ_in` at `w=½`.
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`
    — `scheme.py:1382` — `b = QV·inverse_denom/w`; *"DD's historical 2·QV·inverse_denom
    is the ``w=½`` special case"* — the equation's `2QΔx/W/(…)` term.
  - `orpheus.transport.spatial.diamond._DD_W` (**`data`**) — `diamond.py:102` — the
    `w = ½` that turns each generic body above into *this* equation.
- **implementers — TIER 3 (what makes the slab collapse exact)**:
  - `orpheus.geometry.reduced_operator.slab_streaming` — the factory that populates
    the neutral curvature (`A_in = A_out = 1`, `ΔA/w = 0`, `α = 0`, `τ = 1`) without
    which the one generic body would not reduce to the slab form.
- **confidence**: high that Tier 1 is complete for the production path. Medium on how
  many tiers the declarer wants — Tier 2 bodies are shared with `wdd-face`,
  `hebert-3-435` and every other `w`-scheme, so declaring them here makes those labels
  siblings on one body (correct, but a policy choice).

---

## `dd-curvilinear-scalar`  (13 claims) — `docs/theory/methods/sn/index.rst:639`

- **verdict**: DECLARABLE
- **rationale comment on the page**: **none** (the page has zero rationale comments).
- **what the equation says**: the curvilinear (sphere / non-degenerate cylinder) DD +
  Morel–Montry cell update —
  `ψ̄ = [QV/W + |μ|(A_in+A_out)ψ^s_in + (ΔA/w)c_in ψ_{n−½}] /
        [2|μ|A_out + (ΔA/w)c_out + Σ_t V]`, with the closures `ψ^s_out = 2ψ̄ − ψ^s_in`
  and `ψ_{n+½} = (ψ̄ − (1−τ)ψ_{n−½})/τ`.
- ⭐ **Same formula as `dd-solve`** (`docs/theory/methods/sn/curvilinear_one_group.rst:1859`)
  — `index.rst` restates it in the strategy chapter. **Check whether `dd-solve` is
  already declared before landing this**, or the two labels will duplicate.
- **implementers — TIER 1 (the arithmetic)**:
  - `orpheus.transport.spatial.cell_balance.cell_balance_terms`
    — `orpheus/transport/spatial/cell_balance.py:268`, body at `:343-355`:
    `denom = 2.0*abs_mu*A_downstream + dA_w*c_out + total_xs*V` ;
    `numer_upstream = abs_mu*A_total*ψ^s_in + dA_w*c_in*ψ_ang`. The module docstring's
    "Mathematical content" section displays exactly this equation.
  - `orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
    — `orpheus/transport/spatial/cell_balance.py:123` — the **vectorised (apply-direction)
    twin** of the same algebra, consumed by the SN matvec; Pattern-2 single source.
  - `orpheus.transport.spatial.diamond.DiamondDifference.update`
    — `orpheus/transport/spatial/diamond.py:177` — solves `ψ̄ = (source + numer)/denom`
    and applies **both** closures (`:216-219` spatial, `:224-228` angular
    `(psi_avg − (1−tau)·angular_upstream)/tau`).
  - `orpheus.transport.spatial.diamond.DiamondDifference.residual`
    — `orpheus/transport/spatial/diamond.py:240` — the apply-direction residual
    `denom·ψ̄ − (source + numer_upstream)` of the same system.
  - `orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`
    — `orpheus/transport/spatial/diamond.py:572` — the **curvilinear production route**
    (its σ_t-epoch `denom` carries the `(ΔA/w)c_out` redistribution term explicitly).
- **implementers — TIER 2 (the algebra-of-record / SymPy source of truth)**:
  - `orpheus.derivations.discrete.sn.balance.derive_wdd_solve`
    — `orpheus/derivations/discrete/sn/balance.py:229` — derives and SymPy-verifies
    this exact solved form. The module header declares itself *"the **source of truth**
    for the balance equations in the SN theory chapters under `docs/theory/methods/sn/`"*
    and lists *"5. WDD substitution → solved form (Eq. dd-solve)"*. Its
    `expected_denom` / `expected_numer` are the equation term-for-term.
- **implementers — TIER 3 (`w`-generic primitives at `w=½`)**: same four as
  `dd-slab-scalar` Tier 2 (`cell_average`, `outgoing_face_from_average`,
  `source_emission`, `_DD_W`).
- **⚠ boundary**: the constants `c_out = α_{n+½}/τ`, `c_in = (1−τ)/τ·α_out + α_in`
  carry their **own** label `dd-mm-closure-constants` (`index.rst:628`, **not in this
  batch**); they are produced by the angular closure
  (`MorelMontryAngularSweep._c_out_per_level` / `_c_in_per_level`) and derived by
  `derive_wdd_solve`. `cell_balance_terms` explicitly *"does not derive them"*
  (`:352-356`). Do not attribute the c-constants here.
- **confidence**: high on Tier 1 and Tier 2. The `dd-solve` overlap is the one thing
  that could change the answer — if `dd-solve` is declared elsewhere, this label should
  either share those implementers or be resolved as a duplicate at the doc level.

---

## `dd-cylindrical-degenerate`  (1 claim) — `docs/theory/methods/sn/index.rst:671`

- **verdict**: DECLARABLE
- **rationale comment on the page**: **none** (the page has zero rationale comments).
- **what the equation says**: the cylindrical pure-azimuthal degenerate specialisation
  (`|η| → 0`, no radial face flow) — `denom = (ΔA/w)c_out + Σ_t V`,
  `numer = QV/W + (ΔA/w)c_in ψ_{n−½}`; i.e. the `2|μ|A_out` and `|μ|(A_in+A_out)ψ^s_in`
  terms drop out and no face flux is written.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.transport.spatial.cell_balance.cell_balance_terms`
    — `orpheus/transport/spatial/cell_balance.py:268` — the degenerate case is realised
    **by data, not by a branch**: `A_downstream = 0.0` kills `2|μ|A_down`, and the
    module docstring names it as one of the three collapsed cases —
    *"**Cylindrical pure-azimuthal degenerate** — `CellVisit.face_area_downstream = 0.0`
    (geometric truth: no radial face on this ordinate). The `2|μ|·A_down` term vanishes
    via `A_down = 0` (**not** via the numerical threshold `|μ| < 1e-15`)."*
    ⭐ Note this **contradicts the page**, which still describes the `abs_mu < 1e-15`
    dispatch. The tree's mechanism is the geometric zero.
  - `orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
    — `orpheus/transport/spatial/cell_balance.py:123` — the vectorised twin; its
    `A_downstream` parameter documentation names the case:
    *"broadcast-scalar 0.0 for the degenerate pure-z / pole cell (no outgoing face)"*.
  - `orpheus.transport.spatial.diamond.DiamondDifference.update`
    — `orpheus/transport/spatial/diamond.py:177` — implements the equation's *second*
    half (the "no face-flux write" signal): `if visit.face_area_downstream > 0.0:` …
    else `psi_spat_out` stays `None`, with the comment *"Outputs ``None`` when there is
    no downstream spatial face on this visit (cylindrical pure-azimuthal degenerate:
    face_area_downstream == 0.0)"*. This is the ONLY body that realises the
    `outgoing_spatial_flux = None` contract the page describes.
- **confidence**: high. Two helpers plus the one gate; all three name the degenerate
  case in their own text. ⚠ The `1e-15` threshold the page cites is **gone** from this
  path — declaring is still right, but the page needs the correction.

---

## `dd-recurrence`  (6 claims) — `docs/theory/methods/sn/slab_one_group.rst:358`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none attached to `dd-recurrence` itself. The two
  neighbouring rationales (`:389`, `:404`, on `sweep-cumprod-factors` and
  `sweep-cumprod-solution`) both point at it as their terminal result:
  > *"Not a standalone solver claim; the terminal result — **the recurrence itself** —
  > is pinned against **the symbolic derivation** by
  > `tests/sn/sweep/slab/test_dd_recurrence.py::test_dd_per_cell_recurrence_matches_symbolic_derivation`."*
  ⭐ *"the symbolic derivation"* is the pointer that resolves this equation.
- **what the equation says**: the DD slab face recurrence `ψ_out = a_i·ψ_in + b_i`,
  obtained by substituting the DD closure `ψ_out = 2ψ_avg − ψ_in` into
  `dd-cartesian-1d`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence`
    — `orpheus/derivations/discrete/sn/balance.py:296` — **the algebra of record**. Its
    docstring is the equation (*"Derive Eq. dd-recurrence: ψ_out = a·ψ_in + b"*), it
    performs the substitution symbolically and asserts
    `a·ψ_in + b − ψ_out == 0`, and the module header declares itself *"the **source of
    truth** for the balance equations … If an equation in the RST cannot be derived
    from this script, it must be added here first"*, listing *"6. Cumprod recurrence
    coefficients (Eq. **dd-recurrence**)"*.
  - `orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`
    — `orpheus/transport/spatial/diamond.py:572` — **the production producer of `a`**
    (and of `1/denom`, from which `b` follows through `source_emission`). Docstring:
    *"The single source of the DD affine-recurrence coefficients (Cardinal Rule 2)."*
    ⚠ Note the page separately labels the coefficients `dd-recurrence-coefficients`
    (`slab_one_group.rst:365`, not in this batch) — if that label is declared, this
    symbol arguably belongs there and `dd-recurrence` keeps only the *statement*
    implementers below.
  - `orpheus.sn.sweep.scan.ordinate_scan`
    — `orpheus/sn/sweep/scan.py:82` — **evaluates** the recurrence: its module docstring
    opens *"the canonical closed-form solution to the affine recurrence
    ψ[i+1] = a[i]ψ[i] + b[i]"*. (Shared with `blelloch-1990-eq-1-5`, which is the
    closed FORM; `dd-recurrence` is the recurrence that form solves.)
  - `orpheus.sn.loss_representation.CumprodScan`
    — `orpheus/sn/loss_representation/__init__.py:1359` (**class**) — the loss
    representation the page names by `:meth:`CumprodScan.sweep``; it assembles `a` and
    `b` per chain and runs the scan. Include if the ledger wants the strategy that
    *is* the recurrence's production realisation; exclude if only arithmetic counts.
  - `orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`
    — `orpheus/transport/spatial/diamond.py:694` — the multi-D-Cartesian sibling
    producer of the same `(a, b)` shape for a scanned row. ⚠ Different geometry family
    (`dd-cartesian-2d`'s scan); listed for completeness, declare only if
    `dd-recurrence` is read as geometry-generic.
- **confidence**: high on `derive_cumprod_recurrence` (named as the source of truth,
  by label) and on `affine_scan_coefficients` / `ordinate_scan`. Medium on
  `CumprodScan` and `cartesian_scan_coefficients` — those are the "which layer counts
  as implementing a recurrence" question, which is a project ruling.

---

## `addition-theorem`  (1 claim) — `docs/theory/methods/sn/slab_multigroup.rst:249`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none (the file's 5 rationale comments are all
  ≥ 272 lines below, on the eigenvalue-posing equations). The page's own prose is
  directive: *"The normalisation is the **"no** `4π/(2ℓ+1)` **prefactor"**
  convention under which the addition theorem reads [eq] which is the identity used
  by Eq. :eq:`pn-scatter` to expand the `P_ℓ` scattering kernel as a finite tensor
  product over `m`"*, and it names the gate
  (`test_spherical_harmonics_addition_theorem_L3`).
- **what the equation says**: `Σ_{m=−ℓ}^{ℓ} Y_ℓ^m(Ω)·Y_ℓ^m(Ω′) = P_ℓ(Ω·Ω′)` — the
  property that the project's un-prefactored real-SH normalisation is *chosen* to
  satisfy, and that licenses the `(2ℓ+1)` reconstruction factor.
- ⛔ **DUPLICATE LABEL — check before landing.** The identical identity carries a
  second label, `real-sh-addition-theorem`, at
  `docs/theory/foundations/spherical_harmonics.rst:46` **and inside the code**, at
  `orpheus/numerics/basis/spherical_harmonic_basis.py:24` (the module docstring
  declares `.. math:: :label: real-sh-addition-theorem`). So the tree already
  *declares* an owner for the identity — under the OTHER label. Declaring
  `addition-theorem` against the same symbols is right, but the duplication itself
  is a doc-hygiene finding worth an issue.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
    — `orpheus/numerics/basis/spherical_harmonic_basis.py:120` (**class**) — *"carries
    the SH convention (the no-prefactor normalisation, the addition-theorem factor
    `2ℓ+1`, the continuous Gram diagonal `4π/(2ℓ+1)`)"*, and its module docstring
    carries the equation. This is the object whose normalisation choice IS the theorem.
  - `orpheus.numerics.basis.spherical_harmonic_basis._evaluate_real_sh`
    — `orpheus/numerics/basis/spherical_harmonic_basis.py:501` — the evaluator that
    applies `sqrt(2(ℓ−m)!/(ℓ+m)!)` with the Condon–Shortley phase removed; **that
    normalisation is exactly what makes `Σ_m Y Y = P_ℓ` true**, so it is where the
    identity is realised numerically.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.evaluate`
    — `:186` — the public tabulation `Y[n, l, l+m]` the claiming gate consumes.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.addition_theorem_factor`
    — `:162` (**method** node, a `cached_property`) — the `(2ℓ+1)` array, named for
    this theorem: *"Used by the addition-theorem reconstruction `R = (2ℓ+1)Y`"*.
    ⚠ Note `callers()` under-reports `cached_property` reads (they are property-reached);
    it is live via `reconstruct` / `reconstruct_transpose`.
  - `orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.reconstruct`
    — `:329` — *"The addition-theorem (canonical-dual) synthesis"*; the
    `einsum("nlm,l,lm...->n...", table, addition_theorem_factor, coefficients)` that
    the theorem licenses.
  - `orpheus.numerics.quadrature.directional.Quadrature.spherical_harmonics`
    — `orpheus/numerics/quadrature/directional.py:477` — the surface the claiming
    gate (`tests/sn/operators/test_solver_components.py:564`) actually calls; it
    delegates to the basis table.
- **candidate, weaker**: `MirrorEvenSphericalHarmonicBasis` (`:397`) subclasses the
  basis and masks the σ_y-odd slots — it inherits the convention rather than
  restating it. Declare only if the ledger wants both members of the family.
- **confidence**: high that the SH basis family is the answer and that nothing else
  in `orpheus/` computes `Σ_m Y_ℓ^m Y_ℓ^m`. Medium on WHICH member — the theorem is a
  property of the normalisation, so `SphericalHarmonicBasis` (+ `_evaluate_real_sh`)
  is the tightest honest pair; the other three are consequences. ⚠ Resolve the
  duplicate-label question first.

---

## `normalization-dd-source-coefficient`  (1 claim) — `docs/theory/conventions/normalization.rst:57`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none — `docs/theory/conventions/normalization.rst`
  contains **zero** `.. (vv-status rationale)` comments (`[M]` `grep -c` = 0). But the
  page's prose and the *claiming test* between them settle it exactly (see below).
- **what the equation says**: the 1-D DD fast path's affine source coefficient,
  `b = 2Δx(Q/W)/(2μ + ΔxΣ_t)` — the ERR-025 record: the `1/W` belongs in the source,
  once, and dropping it hid behind eigenvalue invariance.
- **decisive evidence — the claiming gate**:
  `tests/sn/sweep/slab/test_dd_recurrence.py:37`
  `@pytest.mark.verifies("normalization-dd-source-coefficient")` on
  `test_dd_per_cell_recurrence_matches_symbolic_derivation`, which builds
  `source = Q_val * dx_val / W` (`:96`, comment *"The contract source is
  Q · V · weight_norm = Q · dx / W"*), calls `DiamondDifference.update`, and compares
  against `b_sym.subs(..., S: Q_val/W)` from `derive_cumprod_recurrence`.
- **implementers** (complete list, each verified to resolve):
  - `orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`
    — `orpheus/transport/spatial/scheme.py:1382` — `b = QV · inverse_denom / w`.
    At DD (`w = _DD_W = ½`), slab (`V = Δx`, `denom = 2μ + ΔxΣ_t`) and
    `QV = Δx·(Q/W)` this is **exactly** the equation. Its docstring:
    *"DD's historical `2·QV·inverse_denom` is the `w=½` special case."*
  - `orpheus.transport.operators.scattering.ScatteringOperator._assemble_per_ordinate_source`
    — `orpheus/transport/operators/scattering.py:711` — **the `/W` site the page's
    law points at**: *"the division by W happens **once**, at the reconstruction site
    of the scattering application"*. The method's own docstring says
    *"The **single source of truth for the producer-side** `1/W` **combine**… the `/W`
    convention lives HERE, once"*; body at `:756-757`:
    `sum_w = float(self.quadrature.weights.sum())` ; `return (iso / sum_w) + aniso`.
  - `orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence`
    — `orpheus/derivations/discrete/sn/balance.py:296` — the symbolic
    `b = 2S/(2μ/Δx + Σ_t)` that the gate substitutes `S = Q/W` into; the algebra of
    record.
  - `orpheus.transport.spatial.diamond.DiamondDifference.update`
    — `orpheus/transport/spatial/diamond.py:177` — the SUT the gate exercises; it
    consumes the already-`/W`-normalised `source` and folds it as
    `(source + numer_upstream)/denom` (Step 2.5 replaced the explicit `b` with this
    algebraically-identical, ULP-different form — the gate's own note).
  - `orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`
    — `orpheus/transport/spatial/diamond.py:572` — produces the `inverse_denom`
    (`1/(2μ + ΔxΣ_t)` at slab) that `source_emission` multiplies.
- **confidence**: high. The equation is a two-factor product and both factors have a
  single, self-documented home (`source_emission` for the `2ΔxQ/denom` shape,
  `_assemble_per_ordinate_source` for the `/W`), with the symbolic derivation and the
  SUT both named by the one claiming test.

---

## `multigroup`  (41 claims — the largest) — `docs/theory/methods/sn/slab_multigroup.rst:83`

- **verdict**: DECLARABLE
- **rationale comment on the page**: none on this label (the file's 5 rationale
  comments sit at `:521`, `:623`, `:646`, `:692`, `:707` — all on the eigenvalue-posing
  and power-method equations far below). ⭐ The page's **Key Facts** header (`:36-66`)
  is the authored routing, and it enumerates the terms operator by operator (quoted
  per implementer below).
- **what the equation says**: the governing multigroup SN transport equation —
  `streaming + Σ_{t,g}ψ_g = (1/W)[ Σ_{g'} Σ_s^{g'→g} φ_{g'} + (χ_g/k) Σ_{g'} νΣ_{f,g'} φ_{g'} ]`.
  Its organising fact, stated immediately after: **the LHS is diagonal in the group
  index**; all group coupling sits on the RHS inside the sources.
- **implementers, BY TERM** (complete list, each verified to resolve):

  *LHS — `streaming + Σ_{t,g}ψ_g`:*
  - `orpheus.sn.operators.streaming.StreamingOperator`
    — `orpheus/sn/operators/streaming.py:169` (**class**) — the σ-free `L` (#257 S8b:
    *"the streaming leaf reads no σ"*).
  - `orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
    — `orpheus/transport/operators/multiplication_operator.py:130` (**class**) — `C = M[σ_t]`,
    the collision diagonal built from the typed `total_cross_section_field`.
  - `orpheus.sn.operators.streaming.StreamingCollisionOperator`
    — `orpheus/sn/operators/streaming.py:445` (**class**) — the fused, invertible `L+C`
    whose `.solve` IS the WDD sweep. Key Facts: *"the streaming-collision composite
    `L + C` … is a **direct sum over groups** of the one-group operator … (in the
    implementation: one vectorised sweep with a trailing group axis)"* — i.e. this
    class is where the equation's group-diagonality is realised.

  *RHS scattering — `Σ_{g'} Σ_s^{g'→g} φ_{g'}`:*
  - `orpheus.transport.operators.scattering.LegendreMomentScattering`
    — `orpheus/transport/operators/scattering.py:115` (**class**) — the group-to-group
    moment transfer carrying the page's `SigS[g_from, g_to]` **transpose** convention
    (Key Facts: *"the in-scatter source therefore uses the **transpose**,
    `Q = SigS^T φ`"*).
  - `orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
    — `orpheus/transport/operators/isotropic_scattering.py:229` (**class**) — the `ℓ=0`
    energy transfer (`K_iso`'s `Σ_s0` half).
  - `orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
    — `orpheus/transport/operators/isotropic_scattering.py:345` (**class**) — the `2Σ_2n`
    half of the same isotropic kernel.
  - `orpheus.transport.operators.scattering.ScatteringOperator`
    — `orpheus/transport/operators/scattering.py:356` (**class**) — the full `S`,
    including the equation's `1/W` (Key Facts: the factored `S_aniso = R Λ M`).
  - `orpheus.transport.operators.scattering.ScatteringOperator._assemble_per_ordinate_source`
    — `:711` — the equation's leading `(1/W)`, applied **once**
    (`return (iso / sum_w) + aniso`, `:757`).

  *RHS fission — `(χ_g/k) Σ_{g'} νΣ_{f,g'} φ_{g'}`:*
  - `orpheus.transport.operators.fission.FissionOperator`
    — `orpheus/transport/operators/fission.py:147` (**class**) — `F = χ ⊗ νΣ_f`, exactly
    the Key Facts statement *"Fission is **rank-1 in energy**: `(Fφ)_g = χ_g Σ_{g'}
    νΣ_{f,g'} φ_{g'}` — an outer product of spectrum and rate with no useful inverse;
    the operator is `apply`-only by construction."*

  *The assembly (the equation as one object):*
  - `orpheus.sn.coupled_system.build_within_group_system`
    — `orpheus/sn/coupled_system.py:425` — *"The single source of truth (Cardinal
    Rule 2) for the within-group decomposition every solve consumes"*; it builds
    `L`, `C`, `L+C`, `S`, `B_a` and the loss/gain grids. This is the honest
    "the equation, assembled" node.
  - `orpheus.sn.coupled_system.WithinGroupSystem`
    — `orpheus/sn/coupled_system.py` (**class**) — the typed carrier of that
    decomposition.
- **deliberately EXCLUDED, with the structural reason**:
  - `SNBoundaryOperator` realises `B`, which does **not** appear in this equation
    (the page introduces `B` in the *within-group operator equation* one paragraph
    later, and the boundary algebra is `boundary_conditions.rst`'s).
  - `solve_sn` / `SNSolver` are the *drivers*; the `1/k` self-consistency they run is
    the eigenvalue posing, which carries its own labels
    (`sn-mg-eigenvalue-posing`, `eigenvalue-posing`, and the four rationale-marked
    equations at `:521`–`:707`). Declaring the drivers here would swallow those.
- **confidence**: high on the term→operator mapping (each is quoted from the page's own
  Key Facts and each class's own docstring). Medium on **breadth policy**: `multigroup`
  is a *governing* equation, so a ledger convention that a governing equation is
  implemented by its assembler alone would keep only `build_within_group_system`
  (+ maybe `StreamingCollisionOperator`), while a term-wise convention keeps all
  nine. I list all nine because the contract here is complete enumeration; **prune by
  policy, not by search.** What would change it: a project ruling on governing-equation
  granularity.

---

## `phase-f-q-bar-twin-forms`  (1 claim) — `docs/theory/methods/sn/curvilinear_numerics.rst:1462`

- **verdict**: DECLARABLE (for the surviving twin) — ⚠ **with a doc-drift finding**
- **rationale comment on the page**: none attached to this label; the two rationales
  bracketing it belong to its inputs. The nearer one (`:1449`, on
  `phase-f-source-eq-sigt-phi0`) is directly relevant and reads:
  > *"(vv-status rationale) Derivation step: the SI fixed-point identity
  > (Σ_t·φ₀ = Q_1d) that makes the sweep-path Carlson seed canonically equivalent to
  > the apply-path form. Not a standalone claim; verified transitively by the SI
  > convergence infrastructure — off the fixed point the two forms differ by the SI
  > residual, which vanishes at convergence."*
- **what the equation says**: at `μ = −1`, the cell-averaged starting-direction source
  admits two equal-on-the-fixed-point expressions — `Q̄_i = ½Σ_t(r_i)φ_0(r_i)`
  (**apply path**, builds φ₀ from the input ψ) and `Q̄_i = ½Q_1d(r_i)`
  (**sweep path**, takes the within-group source directly).
- ⛔ **ONE OF THE TWO TWINS NO LONGER EXISTS** (`[M]` HEAD `a1c90aac`). The apply-path
  form `Q̄ = ½Σ_tφ₀` was the `CarlsonInwardSweep` **PROXY source**, retired at
  #282 route (a) / #280 Phase 2.5d. Evidence, three independent:
  1. `orpheus/sn/sweep/psi_half_angle_seed.py:51-55` — *"`CarlsonInwardSweep` — this
     module's march driven by the PROXY source `Q̄ = Σ_tφ_0/Σw`, exact only at the
     flat-flux equilibrium (ERR-058 b — O(1) wrong off equilibrium, floored the
     curvilinear MMS at ~0.04 L2 independent of mesh, Issue #195)"*, followed by
     *"Route (a) (#282, ruling R10) **retired the whole strategy family**"*.
  2. The **claiming gate itself** (`tests/sn/sweep/core/test_phase_c_gates.py:888-895`)
     says: *"the seed-strategy zoo is retired … the old test compared it against the
     retired `CarlsonInwardSweep` STRATEGY wrapper on a flat-ψ probe — a
     self-comparison"* — so the gate no longer compares the twins at all.
  3. `[M]` `grep -rn "CarlsonInwardSweep" orpheus/` → **1 hit, past-tense history**;
     no `0.5 * sig_t * phi` construction anywhere in `orpheus/`.
  ⟹ **the page is present-tense-false about this equation** — worth an issue.
- **implementers of the SURVIVING twin** (each verified to resolve):
  - `orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
    — `orpheus/numerics/spaces/radial_characteristic_space.py:525` — the production
    `Q̄(μ=±1)` construction. Its own docstring: *"At `ℓ = 0` this reduces to
    `Q̄ = ½Q_0` — the same convention `carlson_inward_sweep_from_source`'s `Q_bar`
    parameter documents"*, which is the equation's sweep-path line, generalised to
    all `ℓ` (ruling R14). *"The single source of the q½ source construction."*
  - `orpheus.numerics.spaces.radial_characteristic_space._radial_characteristic_reconstruction_weights`
    — `:490` — the `(2ℓ+1)/2·(±1)^ℓ` weight; at `ℓ=0` it IS the equation's `½`, which
    the page derives as *"the Legendre fold weight `(2ℓ+1)/2` at `ℓ = 0` times
    `P_0(−1) = 1`"*. **The single source of the `½`.**
  - `orpheus.sn.operators.radial_characteristic.RadialCharacteristicReconstruction`
    — `orpheus/sn/operators/radial_characteristic.py:959` (**class**) — the operator
    face of the same fold (`A_BA`'s shared factor).
  - `orpheus.transport.radial_characteristic_field.RadialCharacteristicField.source_from_angular`
    — the data-factory that folds a *per-ordinate* source through the same kernel
    after a Legendre projection (*"a different typed input, not a twin"*).
  - `orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
    — `orpheus/sn/sweep/psi_half_angle_seed.py:112` — the consumer whose `Q_bar`
    parameter documentation carries the convention (*"For an isotropic source this is
    `½ Q̄_iso`"*). Include if the ledger wants the convention's consumer; the
    arithmetic is in the two fold symbols above.
- **NOTHING implements the apply-path twin** — retired, deliberately, because it is
  only exact on the fixed point (exactly what the equation says, and what ERR-058b
  measured as O(1) wrong off it).
- **confidence**: high on the retirement (three independent sources, including the
  claiming gate). High on the surviving fold. ⚠ Medium on whether the label should be
  declared at all in its current wording — a label asserting an equivalence between a
  live form and a retired one is arguably a `NOTHING:canonical-form` for its
  apply-path half. **Recommend: fix the page first, then declare the fold.**

---

## `si-spectral-rate`  (1 claim) — `docs/theory/methods/sn/slab_one_group.rst:699`

- **verdict**: DECLARABLE (both halves, by different symbols) — the brief's guess that
  this "may be a reference value" is **half right**: the left identity is a theorem
  with no routine, the right identity is a live production law.
- **rationale comment on the page**: none on this label. The nearest above (`:678`,
  on `si-convergence-residual`) and below (`:769`, on `si-sigma-r-fold-mismatch`)
  belong to other equations. ⭐ The page's `.. note::` immediately BELOW the equation
  is the authored answer and it names a symbol:
  > *"The `c` in :eq:`si-spectral-rate` is the **within-group scattering ratio**
  > `Σ_{s,0}^{g→g}/Σ_{t,g}` … The :meth:`Mixture.scattering_ratio
  > <orpheus.data.macro_xs.mixture.Mixture.scattering_ratio>` property exposes the
  > slightly larger **Case–Zweifel** secondaries-per-collision parameter
  > `c_g = (Σ_{s,g} + νΣ_{f,g})/Σ_{t,g}` … The L1 rate anchor … pins `n_Jacobi`
  > against `log(tol)/log(c_max)` using the Case–Zweifel form."*
- **what the equation says**: two identities —
  (a) `ρ_J = ρ((L+C)^{-1}(S+B)) = c ≡ max_g Σ_{s,g}/Σ_{t,g}` (the Fourier/mode result),
  (b) `n_Jacobi ≈ log ε / log c` (iterations to a tolerance at that rate).
- **implementers** (complete list, each verified to resolve):
  - `orpheus.data.macro_xs.mixture.Mixture.scattering_ratio`
    — `orpheus/data/macro_xs/mixture.py:255` (**method**, a `property`) — the `c`.
    ⚠ **Honest caveat the page itself states**: it returns the Case–Zweifel
    `(Σ_s + νΣ_f)/Σ_t`, *"slightly larger"* than the equation's `Σ_s/Σ_t`; its own
    docstring says *"For purely scattering media this reduces to `Σ_s/Σ_t ≤ 1`."*
    So it implements the equation's `c` exactly in the non-multiplying case and a
    documented generalisation otherwise. The claiming gate
    (`tests/sn/verification/analytical/test_si_convergence_rate.py:289`) uses it:
    `c_max = float(np.max(np.asarray(m.scattering_ratio)))`.
  - `orpheus.numerics.convergence._budget_from_law`
    — `orpheus/numerics/convergence.py:299` — **the single home of identity (b)**:
    `max(1, math.floor((math.log(tolerance) - log_initial)/log_rate) + 2)`, i.e.
    `N = ⌊(ln tol − a)/b⌋ + 2`; at `a = 0`, `b = ln ρ` this is `log ε / log c`. Its
    docstring: *"The geometric budget law, in ONE place … a second spelling of it is
    a twin waiting to drift."*
  - `orpheus.numerics.convergence.default_iteration_budget`
    — `orpheus/numerics/convergence.py:326` — the **a-priori** evaluation of (b)
    (*"Pass a measured `rho` to size a budget for a KNOWN problem"*).
  - `orpheus.numerics.convergence.StoppingCriterion.projected_iterations`
    — `orpheus/numerics/convergence.py:627` (**method**) — the **a-posteriori**
    evaluation, docstring: *"`n = |ln(tol/r_0)| / |ln ρ|`"* — identity (b) with the
    observed `ρ`.
  - `orpheus.numerics.convergence.StoppingCriterion.rate`
    — `orpheus/numerics/convergence.py:618` (**method**) — *"The observed geometric
    decay rate `ρ` over the tail"*: the measured left-hand `ρ_J`.
- **candidates, weaker (declare only if the ledger wants the constant / the resolver)**:
  - `orpheus.numerics.convergence._SERVED_RATE` (**`data`**, `:296`) — the served
    contraction rate `0.986` the a-priori budget is sized against; it is a *promise*
    about `ρ`, not the equation's `c`.
  - `orpheus.numerics.convergence.resolve_iteration_budget` (`:399`) — routes
    `None → default_iteration_budget`; plumbing.
- **what NOTHING implements — and this is deliberate**: the identity
  `ρ((L+C)^{-1}(S+B)) = c` itself. No symbol assembles the SI iteration matrix and
  takes its spectral radius; it is a Fourier-analysis theorem (Lewis & Miller §4.4)
  and is *verified* by counting iterations, not computed. ⛔ Do **not** attribute
  `FaceTransmissionSpectrum.spectral_radius`
  (`orpheus/transport/spatial/scheme.py:311`) here — that is a **different ρ** (the
  DD face-mode sawtooth, `ρ = 1.0` for diamond), and the page's own
  `si-sigma-r-fold-mismatch` rationale warns about exactly this kind of operator
  confusion.
- **confidence**: high. Two identities, two homes, both named in the tree's own
  docstrings, and the one claiming gate exercises both.

---

## `sn-p1-cylinder-hand-ref`  (1 claim) — `docs/theory/verification/sn.rst:1830`

- **verdict**: **NOTHING:`canonical-form`**
- **rationale comment on the page**: none on this label. Its sibling on the SAME page
  (`:2799`) is the closest authored signal and it names the KIND:
  > *"(vv-status rationale) reference: the two hand-derived face Legendre …"*
  — the page uses a `reference` rationale kind for hand-derived oracles.
- **what the equation says**: the structurally-independent hand reference for the
  cylindrical `ℓ=1` scattering source —
  `q_n^{P_1} = (1/W)·3·Σ_{s1}·Σ_m Y_1^m(Ω_n)·φ_1^m` with
  `φ_1^m = Σ_n w_n Y_1^m(Ω_n) ψ_n`.
- **why NOTHING, with positive evidence (not absence of a quick hit)**:
  1. **The page declares it an oracle by construction**: the section is titled *"L0 —
     the operator-admits trick"* and describes the cylinder row as
     *"explicit `Y_1^m` moment-sum, **independent of the production `R Λ M` einsum**"*.
  2. **The claiming test builds it inline and says why**
     (`tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:218-259`):
     *"The hand-reference assembles the moment sum by an INDEPENDENT explicit loop
     over the `spherical_harmonics` table, **NOT the production frame
     analysis/reconstruction faces / `LegendreMomentScattering` einsums** — so a
     transposed or mis-shaped einsum in the production path is detectable."*
     The formula lives at `:241-256` as a double `for` loop inside the test body.
  3. `[M]` searched `orpheus/` for the construction and for a stored reference: no
     `hand_ref`/`hand-ref` symbol, no P1 moment-sum reference constant, and the only
     `sigma_s1` sites are the DSA derivation (`derivations/discrete/sn/dsa.py`) and the
     anisotropic MMS source (`derivations/continuous/mms/sn.py`) — **neither is this
     formula**.
  4. **Declaring the production path here would be a correctness error**, not just a
     stylistic one: `pn-scatter` + `flux-moments` (which the *same test* also claims)
     ARE the production `R Λ M` labels. Pointing `sn-p1-cylinder-hand-ref` at
     `LegendreMomentScattering` / `SphericalHarmonicBasis.reconstruct` would make the
     gate a self-comparison on paper — the exact demotion `coding-standards` warns
     about.
- **KIND justification**: `canonical-form` — *"a form exhibited to show structure that
  no production path takes"*. Production takes the factored `R Λ M` einsum; this
  explicit per-`m` sum is exhibited precisely **because** production does not take it.
  (It is not `identity`/`law`/`definition`: it is a computable expression, just one
  whose only legitimate home is the oracle side of a gate.)
- **⚠ what would change this verdict**: if the project decides a *test-tree* function
  is a legal `implements` origin, the honest implementer is
  `tests.sn.verification.mms.test_curvilinear_aniso_scattering_p1.test_cylindrical_p1_source_matches_hand_reference`
  itself (a `function` node). I do not recommend it — a test is the *claimer*, and an
  `implements` edge from the same node that carries the `verifies` marker makes the
  claim circular.
- **confidence**: high, on positive evidence (the page's own "independent of the
  production einsum" sentence + the test's own "NOT the production … einsums" comment
  + a negative grep). The sibling `sn-p1-sphere-hand-ref` (`:1818`, not in this batch)
  is the identical case.

---

# Close-out — verification of every claim in this file

## Symbol resolution

**63 dotted names cited; 62 resolve as `function`/`method`/`class`/`data`.** The one
non-resolving name is `orpheus.sn.sweep.psi_half_angle_seed` — the **module**, which
appears only in prose (as the label-owning docstring the page points at) and is
**never listed as an implementer** (verified: 0 implementer bullets name it). This is
the ontology refusal the brief warns about, and it is the reason the `hebert-3-43X`
declarations are pinned to the module's **functions** rather than to the module the
page names.

Node types actually used: `function` 26, `method` 21, `class` 14, `data` 2
(`transport.spatial.diamond._DD_W`, `numerics.convergence._SERVED_RATE`).

## Re-run of every claim whose EMPTINESS is the finding

(A "zero hits" claim is the most drift-fragile kind, so each was re-run at close.)

| # | claim | re-run result |
|---|---|---|
| N1 | one Blelloch closed-form body in `orpheus/` | ⚠ **2 files carry `np.cumprod`** — but only `scan.py:177` has the full `cumprod_a * (psi_0 + cumsum(b/cumprod_a))`. `orpheus/sn/sweep/cache.py:505` computes `cumprod_a = np.cumprod(a_attenuation, axis=2)` ALONE — the cached `C_i` prefix product, i.e. `sweep-cumprod-factors` (a **different, unbatched label**), which `affine_scan_coefficients`' docstring explicitly assigns to `CollisionCache`. `np.cumsum` appears once in the whole SN tree, at `scan.py:177`. ✅ claim holds |
| N2 | one τ-from-edges construction | ✅ the 2 grep hits are 2 **lines of the same body** (`pole_angular_closure.py:1154` and `:1158`) |
| N3 | `CarlsonInwardSweep` (the retired `½Σ_tφ₀` proxy) is gone from `orpheus/` | ✅ **1 hit, past-tense history** (`psi_half_angle_seed.py:51`); no live symbol |
| N4 | no P1 hand-reference in `orpheus/` | ✅ 0 |
| N5 | `index.rst` and `normalization.rst` carry no `vv-status rationale` | ✅ 0 and 0 |
| N6 | one `(−1)^k` face-index sawtooth construction | ⚠ **5 files mention "sawtooth"**; 4 are prose about the DD face mode. Only `loss_kernel_gauge.py` constructs it. ✅ claim holds — but see the addendum below |

### Addendum to `dd-null-sawtooth` (from the N6 re-check)

`orpheus/transport/spatial/scheme.py` carries the **precursor** relation, not the
equation: `_face_transmission_matrix` (`:333`) / `_face_transmission_spectrum` (`:372`)
/ `DiscretizationSchemeBase.face_transmission_spectrum` (`:911`) drive one cell with
unit inflows and read `ρ(Σ)`; the `−1` eigenvector they detect is the *cell-local*
`ψ_out = −ψ_in` at `ψ_c = 0` — the page's **unlabelled** `.. math::` at
`cartesian_multid.rst:~4608`, one step above `dd-null-sawtooth`. They do **not** build
the mesh-wide `(−1)^k·φ(i_⊥)` factorisation.
⟹ **do not declare them against `dd-null-sawtooth`**; they are the natural implementers
of that unlabelled precursor, should it ever be given a label.

## Tree state

HEAD moved during this dispatch: `a1c90aac` → `a0a979de` (`docs(rules)` only).
`git status --short` shows **0** modified files under `orpheus/` or `docs/`, so every
`file:LINE` above is current at close.
