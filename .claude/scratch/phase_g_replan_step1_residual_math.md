# Phase G replan Step 1 — salvaged residual-computation math

**Source**: extracted from `orpheus/sn/spatial/operators.py` (committed in `6eeff94`,
to be reverted in Step 0). Step 1 ports these branches into
`DiamondDifference.residual` as three private branches mirroring the
existing `_update_*` branches.

The math is correct (Pattern 2 — both `update` and `residual` consume the
same `cell_balance_terms` helper). Only the *home* changes — the strategy
gets a new method, not a wrapper class.

---

## Slab branch — closed-form

```python
def _residual_slab(
    st,
    cell_avg: np.ndarray,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
    source: np.ndarray,
) -> np.ndarray:
    """Slab DD per-cell balance residual.

    Balance: (2|μ| + chord·Σ_t)·cell_avg = 2|μ|·ψ_in + source.
    Residual = LHS − RHS = 2|μ|·(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source.

    Note: source is already weight-normalised per the CellUpdate contract
    (Q · Δx · weight_norm for slab).
    """
    assert st.abs_mu is not None
    abs_mu = st.abs_mu
    chord = st.chord_length
    psi_in = upstream_state.spatial_upstream
    return (
        2.0 * abs_mu * (cell_avg - psi_in)
        + chord * total_xs * cell_avg
        - source
    )
```

---

## Curvilinear non-degenerate branch

```python
def _residual_curvilinear(
    st,
    A_downstream: float,
    cell_avg: np.ndarray,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
    source: np.ndarray,
) -> np.ndarray:
    """Curvilinear non-degenerate per-cell balance residual.

    Consumes :func:`cell_balance.cell_balance_terms` so the algebra is
    bit-identical to :meth:`DiamondDifference._update_curvilinear` by
    construction (Pattern 2 — qa CONCERN-2).

    Residual = denom · cell_avg − (source + numer_upstream).
    """
    from .cell_balance import cell_balance_terms

    terms = cell_balance_terms(st, A_downstream, total_xs, upstream_state)
    return terms.denom * cell_avg - (source + terms.numer_upstream)
```

---

## Cylindrical pure-azimuthal degenerate branch

```python
def _residual_cylindrical_degenerate(
    st,
    cell_avg: np.ndarray,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
    source: np.ndarray,
) -> np.ndarray:
    """Cylindrical pure-azimuthal degenerate per-cell balance residual.

    No radial face flow (abs_mu < 1e-15). Consumes
    :func:`cell_balance.cell_balance_terms_degenerate`.

    Residual = denom · cell_avg − (source + numer_upstream).
    """
    from .cell_balance import cell_balance_terms_degenerate

    terms = cell_balance_terms_degenerate(st, total_xs, upstream_state)
    return terms.denom * cell_avg - (source + terms.numer_upstream)
```

---

## Public dispatch (apply analogue — to live on `DiamondDifference.residual`)

```python
def residual(
    self,
    cell_avg: np.ndarray,
    visit: CellVisit,
    total_xs: np.ndarray,
    source: np.ndarray,
    upstream_state: UpstreamState,
) -> np.ndarray:
    """Per-cell operator action residual L_cell · ψ̄ − q.

    At the converged cell average (i.e., when cell_avg ==
    update(visit, total_xs, source, upstream_state).cell_average_flux),
    the residual is zero to FP rounding.

    Three branches mirror the corresponding _update_* branches.
    """
    st = visit.streaming_terms

    if st.alpha_in is None:
        return self._residual_slab(st, cell_avg, total_xs, upstream_state, source)

    assert st.abs_mu is not None
    if st.abs_mu < 1e-15:
        return self._residual_cylindrical_degenerate(
            st, cell_avg, total_xs, upstream_state, source,
        )

    assert visit.face_area_downstream is not None
    return self._residual_curvilinear(
        st, visit.face_area_downstream, cell_avg, total_xs, upstream_state, source,
    )
```

---

## Test invariants worth porting to `tests/sn/spatial/test_diamond.py::TestResidual`

1. **Apply-vs-solve round-trip** (Gate 2 from the reverted Step 1).
   For every `(geometry, n_groups, regions, ordinate_index)`:
   - `result = dd.update(visit, total_xs, source, upstream_state)`
   - `r = dd.residual(result.cell_average_flux, visit, total_xs, source, upstream_state)`
   - `np.testing.assert_allclose(r, 0.0, atol=1e-13)`

   Cover slab + sphere + cylinder + cyl-degenerate. Heterogeneous + multi-group bias.

2. **Linearity in `cell_avg`**:
   For two probe points `a, b`:
   - `r_a = dd.residual(a, ...)`, `r_b = dd.residual(b, ...)`
   - `r_λa = dd.residual(λ*a + (1-λ)*b, ...)`
   - `np.testing.assert_allclose(r_λa, λ*r_a + (1-λ)*r_b, rtol=1e-12)`

3. **Bit-identity vs the existing `update` math** when source = 0 and the
   `cell_avg` probe equals what `update` returns. This is automatic from
   the helper-sharing, but worth pinning.

4. **Affine structure**: residual is affine in `source` (linear in
   `cell_avg`, affine in `source`). Add the `source` shift, observe the
   residual shifts by `−source`.

---

## Source attribution

- `orpheus/sn/spatial/operators.py:189-313` — public `apply` dispatch.
- `orpheus/sn/spatial/operators.py:319-341` — `_apply_curvilinear_residual`.
- `orpheus/sn/spatial/operators.py:344-362` — `_apply_cylindrical_degenerate_residual`.

The slab branch was inlined inside the public `apply`; extracted here as
a private static method for the Step 1 port.

**After Step 0 lands** (the file is reverted), Step 1 reintroduces these
three branches as private methods on `DiamondDifference`, alongside the
existing `_update_slab`/`_update_curvilinear`/`_update_cylindrical_degenerate`
methods.
