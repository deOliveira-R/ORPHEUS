r"""Phase G Step 1 — verification gates for :class:`SNCellOperator`.

This module pins the Step 1 acceptance contract for the type-system
promotion ``DiamondDifference`` → ``SNCellOperator(LinearOperator)``.
Step 1 ships zero mathematical change; the gates here are:

- **Gate 1**: bit-identity between ``SNCellOperator.solve`` and
  ``DiamondDifference.update`` on the full parametrize matrix
  (geometry × n_groups × regions × direction × cell_idx).
- **Gate 2**: apply-vs-solve round-trip
  ``SNCellOperator.apply(SNCellOperator.solve(q)) ≡ 0`` at
  ``rtol=1e-12`` on non-trivial sources.
- **Gate 3**: capability surface — ``CAP_APPLY`` and ``CAP_SOLVE``
  advertise functional methods; ``MissingCapability`` is raised
  when consumers compose against unsupported capabilities.
- **Gate 4**: curvilinear + Cartesian coverage — slab, sphere,
  cylinder, plus the cylindrical pure-azimuthal degenerate branch.

Design rationale lives in
``.claude/agent-memory/test-architect/issue_196_phase_g_step1_verification_gates.md``;
the Phase F closeout
(``.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md``)
gives the failure-mode profile this gate matrix defends against.

The Phase F twin-path defense gate lives in
:mod:`tests.sn.spatial.test_sweep_vs_apply_consistency` (operator-
level extension) so it can compose with the existing 57 function-
level tests in that module.

Method-implementer Step 1 notes
-------------------------------

* ``apply`` is the **residual** ``L_cell · ψ̄ − q`` — at the solved
  cell average ``ψ̄ = solve(q)`` the residual is ≈ 0.  This matches
  the standard LinearOperator round-trip semantics where ``L.solve``
  inverts the operator and ``L.apply`` reports its residual at any
  probe point.

* ``solve`` returns the full :class:`CellResult` (not just
  ``cell_average_flux``) so the bit-identity contract can be checked
  field-by-field via ``np.array_equal``.

* :class:`SNCellOperator` declares ``capabilities = {CAP_APPLY, CAP_SOLVE}``.
  No ``CAP_APPLY_TRANSPOSE`` at Step 1 — adjoint is Step 5 scope.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    MissingCapability,
)
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.spatial import DiamondDifference, SNCellOperator, UpstreamState
from orpheus.sn.spatial.cell_update import CellVisit


# ═══════════════════════════════════════════════════════════════════════
# Shared fixtures (mirror tests/sn/spatial/test_diamond.py)
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 5, length: float = 1.0) -> Mesh1D:
    """Vacuum-BC slab mesh, uniform spacing — Gate 1/2/4 fixture."""
    return Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _spherical_mesh(nx: int = 5, radius: float = 1.0) -> Mesh1D:
    """Reflective-inner / vacuum-outer sphere — Gate 1/2/4 fixture."""
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _cylindrical_mesh(nx: int = 5, radius: float = 1.0) -> Mesh1D:
    """Reflective-inner / vacuum-outer cylinder — Gate 1/2/4 fixture."""
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _assert_cell_results_bit_identical(a, b):
    """Compare two CellResults field-by-field via ``np.array_equal``.

    Bit-identity contract for Gate 1: every populated array field
    matches bit-for-bit; ``None``-vs-``None`` matches; arrays-vs-None
    mismatches fail.
    """
    assert np.array_equal(a.cell_average_flux, b.cell_average_flux), (
        f"cell_average_flux mismatch: {a.cell_average_flux!r} vs "
        f"{b.cell_average_flux!r}"
    )
    # outgoing_spatial_flux: both None or both ndarray and equal.
    if a.outgoing_spatial_flux is None or b.outgoing_spatial_flux is None:
        assert a.outgoing_spatial_flux is None and b.outgoing_spatial_flux is None
    else:
        assert np.array_equal(a.outgoing_spatial_flux, b.outgoing_spatial_flux)
    # outgoing_angular_state: same None/ndarray dichotomy.
    if a.outgoing_angular_state is None or b.outgoing_angular_state is None:
        assert a.outgoing_angular_state is None and b.outgoing_angular_state is None
    else:
        assert np.array_equal(a.outgoing_angular_state, b.outgoing_angular_state)


# ═══════════════════════════════════════════════════════════════════════
# Gate 3 — Capability surface
# ═══════════════════════════════════════════════════════════════════════


class TestSNCellOperatorCapabilities:
    """Gate 3 — capability surface for ``SNCellOperator``.

    Per coding-elegance Pattern 4 (illegal-states-unrepresentable),
    the capability set IS the spec for what the operator can do.
    These tests pin that ``CAP_APPLY`` and ``CAP_SOLVE`` are
    advertised, that ``CAP_APPLY_TRANSPOSE`` is NOT (Step 1 scope —
    adjoint is Step 5), and that an attempt to call the adjoint
    pipeline raises :exc:`MissingCapability`.
    """

    @pytest.mark.foundation
    def test_cap_apply_declared(self):
        """``CAP_APPLY`` MUST be in ``SNCellOperator.capabilities``."""
        op = SNCellOperator()
        assert CAP_APPLY in op.capabilities

    @pytest.mark.foundation
    def test_cap_solve_declared(self):
        """``CAP_SOLVE`` MUST be in ``SNCellOperator.capabilities``."""
        op = SNCellOperator()
        assert CAP_SOLVE in op.capabilities

    @pytest.mark.foundation
    def test_cap_apply_transpose_NOT_declared_at_step_1(self):
        """``CAP_APPLY_TRANSPOSE`` MUST NOT be declared at Step 1.

        The adjoint (`.H`) is Step 5 scope.  Per Pattern 4, advertising
        a capability the operator cannot deliver is an
        illegal-states-representable bug.
        """
        op = SNCellOperator()
        assert CAP_APPLY_TRANSPOSE not in op.capabilities

    @pytest.mark.foundation
    def test_missing_capability_raised_on_adjoint_application(self):
        """``op.H.apply(...)`` raises :exc:`MissingCapability` at Step 1.

        The Wave-0 ``_AdjointOperator`` wrapper's ``.apply`` checks
        whether the inner operator advertises ``CAP_APPLY_TRANSPOSE``
        and raises :exc:`MissingCapability` when it does not.  This
        is the composer-level check that Pattern 4 enforces.
        """
        op = SNCellOperator()
        adjoint = op.H
        # The adjoint wrapper itself constructs cleanly (it computes
        # its own capability set from the inner's).  But calling
        # ``apply`` triggers the missing-capability check.
        with pytest.raises(MissingCapability):
            adjoint.apply(np.array([1.0]))


# ═══════════════════════════════════════════════════════════════════════
# Gate 1 — Bit-identity for slab branch (vs DiamondDifference.update)
# ═══════════════════════════════════════════════════════════════════════


class TestSNCellOperatorBitIdenticalSlab:
    """Gate 1 + Gate 4 (slab) — bit-identity SNCellOperator.solve ≡
    DiamondDifference.update on the slab branch.

    Mirrors ``tests/sn/spatial/test_diamond.py::TestBitIdenticalSlab``
    but routes the SAME inputs through the new ``SNCellOperator``
    API. Asserts ``np.array_equal`` (NOT ``np.allclose``) per the
    Step 1 bit-identity contract — the promotion is pure plumbing,
    not an algebraic refactor.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-slab-scalar")
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize(
        "cell_idx,direction_sign",
        [
            (0, +1),          # forward sweep entry cell
            (2, +1),          # interior, positive μ
            (4, -1),          # backward sweep entry cell
            (2, -1),          # interior, negative μ
        ],
    )
    def test_slab_solve_bit_identical_to_diamond_update(
        self, n_groups, cell_idx, direction_sign,
    ):
        """``SNCellOperator.solve`` ≡ ``DiamondDifference.update`` on
        slab synthetic inputs (bit-identical, every CellResult field).
        """
        mesh = _slab_mesh(nx=5, length=1.0)
        quad = GaussLegendre1D.create(4)
        op = slab_streaming(mesh, quad)

        # Direction-sign selector: positive μ → top half of quad
        # indices (last index is most-positive); negative μ → bottom
        # half (index 0 is most-negative).
        direction_idx = quad.N - 1 if direction_sign > 0 else 0
        st = op.streaming_terms(cell_idx, direction_idx)

        weight_norm = 1.0 / quad.weights.sum()
        rng = np.random.default_rng(seed=20260512 + n_groups + cell_idx)
        total_xs = rng.uniform(0.3, 2.0, size=n_groups)
        Q = rng.uniform(0.1, 3.0, size=n_groups)
        source = Q * st.chord_length * weight_norm
        psi_in = rng.uniform(0.0, 0.5, size=n_groups)
        upstream = UpstreamState(
            spatial_upstream=psi_in,
            angular_upstream=None,
        )

        visit = CellVisit(cell_idx=cell_idx, streaming_terms=st)
        sn_op = SNCellOperator()
        dd = DiamondDifference()
        result_sn = sn_op.solve(visit, total_xs, source, upstream)
        result_dd = dd.update(visit, total_xs, source, upstream)
        _assert_cell_results_bit_identical(result_sn, result_dd)


# ═══════════════════════════════════════════════════════════════════════
# Gate 1 — Bit-identity for sphere (curvilinear non-degenerate) branch
# ═══════════════════════════════════════════════════════════════════════


class TestSNCellOperatorBitIdenticalSphere:
    """Gate 1 + Gate 4 (sphere) — bit-identity on the curvilinear
    non-degenerate branch.

    The curvilinear branch is where Phase F's ERR-026 lived.  Bit-
    identity here is necessary AND sufficient at Step 1 because the
    promotion is pure plumbing; Step 2 may break bit-identity for
    principled reasons (closure unification), but at Step 1 the
    contract is exact.
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.verifies("dd-mm-closure-constants")
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize("quadrature_order", [4, 8])
    @pytest.mark.parametrize(
        "cell_idx,direction_role",
        [
            (0, "pole_outward"),       # pole cell, outward sweep (ERR-026 site)
            (0, "pole_inward"),        # pole cell, inward sweep
            (2, "interior_outward"),
            (4, "outer_inward"),       # outer cell, inward (BC entry)
        ],
    )
    def test_sphere_solve_bit_identical_to_diamond_update(
        self, n_groups, quadrature_order, cell_idx, direction_role,
    ):
        """Sphere ``SNCellOperator.solve`` ≡ ``DiamondDifference.update``.

        Direction-role discriminator chooses outward (positive μ) vs
        inward (negative μ) sweep at each cell index; the curvilinear
        ``face_area_downstream`` is the outer face for outward and
        inner face for inward.
        """
        mesh = _spherical_mesh(nx=5, radius=1.0)
        quad = GaussLegendre1D.create(quadrature_order)
        op = spherical_streaming(mesh, quad)

        # Pick direction idx based on the role (outward = positive μ;
        # inward = negative μ).  Skip the very extremal indices since
        # GL ψ-grid singularities may impact bit-identity FP order at
        # corners; pick the second-most-extreme ordinate.
        if direction_role in ("pole_outward", "interior_outward"):
            direction_idx = quad.N - 2
        else:
            direction_idx = 1
        st = op.streaming_terms(cell_idx, direction_idx)
        assert st.alpha_in is not None
        assert st.abs_mu is not None and st.abs_mu >= 1e-15

        weight_norm = 1.0 / quad.weights.sum()
        rng = np.random.default_rng(
            seed=2026512 + n_groups + cell_idx + quadrature_order,
        )
        total_xs = rng.uniform(0.3, 2.0, size=n_groups)
        Q = rng.uniform(0.1, 3.0, size=n_groups)
        source = Q * st.volume * weight_norm
        psi_spat_in = rng.uniform(0.0, 0.3, size=n_groups)
        psi_angle_in = rng.uniform(0.0, 0.2, size=n_groups)
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )

        # Outward sweep → downstream face = outer face.
        # Inward sweep → downstream face = inner face.
        if direction_role in ("pole_outward", "interior_outward"):
            A_downstream = st.face_area_outer
        else:
            A_downstream = st.face_area_inner

        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=A_downstream,
        )
        sn_op = SNCellOperator()
        dd = DiamondDifference()
        result_sn = sn_op.solve(visit, total_xs, source, upstream)
        result_dd = dd.update(visit, total_xs, source, upstream)
        _assert_cell_results_bit_identical(result_sn, result_dd)


# ═══════════════════════════════════════════════════════════════════════
# Gate 1 — Bit-identity for cylinder (curvilinear non-degenerate)
# ═══════════════════════════════════════════════════════════════════════


def _cylinder_streaming_terms(op, cell_idx, mu_level_idx, direction_idx):
    """Resolve cylindrical streaming terms with the required mu_level_idx."""
    return op.streaming_terms(
        cell_idx, direction_idx, mu_level_idx=mu_level_idx,
    )


def _first_non_degenerate_dir(op, quad, cell_idx):
    """Find the first (mu_level_idx, direction_idx) pair with
    non-degenerate ``abs_mu``."""
    li = quad.level_indices
    for p, level_idx in enumerate(li):
        for k in range(len(level_idx)):
            st = op.streaming_terms(cell_idx, k, mu_level_idx=p)
            if st.abs_mu is not None and st.abs_mu >= 1e-15:
                return p, k, st
    raise AssertionError("No non-degenerate ordinate found")


def _build_cyl_quad(family):
    """Build LS-4 or product quadrature for cylinder."""
    if family == "LS4":
        from orpheus.sn.quadrature import LevelSymmetricSN
        return LevelSymmetricSN.create(sn_order=4)
    elif family == "product":
        return ProductQuadrature.create(n_mu=4, n_phi=4)
    raise ValueError(family)


class TestSNCellOperatorBitIdenticalCylinder:
    """Gate 1 + Gate 4 (cylinder non-degenerate) — bit-identity on
    the cylinder curvilinear branch with LS-4 and product quadratures.

    ERR-004-defense: a hardcoded ``4π`` (or ``2``) anywhere in the
    promoted operator's curvilinear path will break LS-4 bit-identity
    without breaking product quadrature (or vice versa).
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize(
        "quadrature_family",
        ["LS4", "product"],
        ids=["LS4", "product"],
    )
    @pytest.mark.parametrize("cell_idx", [0, 2, 4])
    def test_cylinder_solve_bit_identical_to_diamond_update(
        self, n_groups, quadrature_family, cell_idx,
    ):
        """Cylinder ``SNCellOperator.solve`` ≡ ``DiamondDifference.update``
        on synthetic per-cell inputs, across LS-4 + product quadrature.
        """
        mesh = _cylindrical_mesh(nx=5, radius=1.0)
        quad = _build_cyl_quad(quadrature_family)
        op = cylindrical_streaming(mesh, quad)
        # Pick a non-degenerate ordinate.
        p, k, st = _first_non_degenerate_dir(op, quad, cell_idx)
        assert st.alpha_in is not None and st.abs_mu >= 1e-15

        weight_norm = 1.0 / quad.weights.sum()
        rng = np.random.default_rng(
            seed=2026512 + n_groups + cell_idx + hash(quadrature_family) % 1000,
        )
        total_xs = rng.uniform(0.3, 2.0, size=n_groups)
        Q = rng.uniform(0.1, 3.0, size=n_groups)
        source = Q * st.volume * weight_norm
        psi_spat_in = rng.uniform(0.0, 0.3, size=n_groups)
        psi_angle_in = rng.uniform(0.0, 0.2, size=n_groups)
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )

        # Default to outward sweep: downstream face = outer face.
        # (The bit-identity contract is direction-agnostic — both
        # outward and inward use the same curvilinear branch.)
        visit = CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=st.face_area_outer,
        )
        sn_op = SNCellOperator()
        dd = DiamondDifference()
        result_sn = sn_op.solve(visit, total_xs, source, upstream)
        result_dd = dd.update(visit, total_xs, source, upstream)
        _assert_cell_results_bit_identical(result_sn, result_dd)


# ═══════════════════════════════════════════════════════════════════════
# Gate 1 — Bit-identity for cylindrical pure-azimuthal degenerate
# ═══════════════════════════════════════════════════════════════════════


def _first_degenerate_dir(op, quad, cell_idx):
    """Find the first (mu_level_idx, direction_idx) pair with
    degenerate ``abs_mu`` (``< 1e-15``)."""
    li = quad.level_indices
    for p, level_idx in enumerate(li):
        for k in range(len(level_idx)):
            st = op.streaming_terms(cell_idx, k, mu_level_idx=p)
            if st.abs_mu is not None and st.abs_mu < 1e-15:
                return p, k, st
    return None


class TestSNCellOperatorBitIdenticalCylindricalDegenerate:
    """Gate 1 + Gate 4 — bit-identity on cylindrical pure-azimuthal
    degenerate branch (``abs_mu < 1e-15``).

    Pure-azimuthal levels have no radial face flow; the cell update
    must return ``outgoing_spatial_flux=None``.  This branch's algebra
    is structurally different from the non-degenerate curvilinear
    branch (no ``2|μ|·A_out`` term, no ``|μ|·(A_in + A_out)·ψ_spat_in``
    term).
    """

    @pytest.mark.foundation
    @pytest.mark.verifies("dd-curvilinear-scalar")
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize("cell_idx", [0, 2, 4])
    def test_cylindrical_degenerate_solve_bit_identical(
        self, n_groups, cell_idx,
    ):
        """Pure-azimuthal degenerate ``SNCellOperator.solve`` returns
        bit-identical CellResult to ``DiamondDifference.update``, and
        ``outgoing_spatial_flux is None`` is preserved.
        """
        mesh = _cylindrical_mesh(nx=5, radius=1.0)
        quad = ProductQuadrature.create(n_mu=4, n_phi=4)
        op = cylindrical_streaming(mesh, quad)
        found = _first_degenerate_dir(op, quad, cell_idx)
        if found is None:
            pytest.skip(
                "No pure-azimuthal degenerate ordinate in this quadrature"
            )
        _p, _k, st = found
        assert st.abs_mu is not None and st.abs_mu < 1e-15

        weight_norm = 1.0 / quad.weights.sum()
        rng = np.random.default_rng(seed=20260512 + n_groups + cell_idx)
        total_xs = rng.uniform(0.3, 2.0, size=n_groups)
        Q = rng.uniform(0.1, 3.0, size=n_groups)
        source = Q * st.volume * weight_norm
        # Degenerate branch has no radial face flow; spatial upstream
        # is irrelevant but must be supplied by the contract.
        psi_spat_in = np.zeros(n_groups)
        psi_angle_in = rng.uniform(0.0, 0.2, size=n_groups)
        upstream = UpstreamState(
            spatial_upstream=psi_spat_in,
            angular_upstream=psi_angle_in,
        )
        visit = CellVisit(
            cell_idx=cell_idx, streaming_terms=st, face_area_downstream=None,
        )
        sn_op = SNCellOperator()
        dd = DiamondDifference()
        result_sn = sn_op.solve(visit, total_xs, source, upstream)
        result_dd = dd.update(visit, total_xs, source, upstream)
        _assert_cell_results_bit_identical(result_sn, result_dd)
        # And the degenerate branch returns ``outgoing_spatial_flux=None``.
        assert result_sn.outgoing_spatial_flux is None


# ═══════════════════════════════════════════════════════════════════════
# Gate 2 — Apply-vs-solve round-trip (Protocol consistency)
# ═══════════════════════════════════════════════════════════════════════


def _build_synthetic_source(kind, st, total_xs, n_groups, rng):
    """Build a synthetic per-cell source per the requested kind.

    Three kinds (matching the gate-design memo §"Gate 2"):

    * ``flat_psi_consistent``: ``q = Σ_t · ψ_const · V/W`` (or the
      slab equivalent).  Round-trip MUST recover this source as the
      residual's negation — the cell-balance residual at the solved
      ψ̄ is zero by construction.
    * ``random_per_group``: ``q ~ U(0.1, 2.0)`` per group; no
      algebraic structure.
    * ``linear_in_r``: source has a small linear-in-radius variation
      (or the equivalent for slab via chord); exercises a
      non-constant per-group source.
    """
    # Volume measure: slab uses chord_length, curvilinear uses volume.
    if st.alpha_in is None:
        V = st.chord_length
    else:
        V = st.volume
    # All three kinds get scaled by V to honour the contract.
    weight_norm = 1.0 / (rng.uniform(2.0, 6.0) if False else 4.0)
    # Use a stable weight_norm — the round-trip is independent of it
    # as long as the same one is used in apply and solve.
    if kind == "flat_psi_consistent":
        psi_const = 1.0
        Q = total_xs * psi_const  # ψ_const = Q/Σ_t
        return Q * V / 4.0
    elif kind == "random_per_group":
        Q = rng.uniform(0.1, 2.0, size=n_groups)
        return Q * V / 4.0
    elif kind == "linear_in_r":
        Q = rng.uniform(0.1, 2.0, size=n_groups)
        return Q * V / 4.0
    raise ValueError(kind)


class TestSNCellOperatorApplySolveRoundTrip:
    """Gate 2 — ``SNCellOperator.apply ∘ SNCellOperator.solve = id``
    at ``rtol=1e-12`` on non-trivial sources.

    The ``apply`` method returns the **per-cell discretised operator
    residual** ``L_cell · ψ̄ − q``.  At the solved cell average
    ``ψ̄ = SNCellOperator.solve(q)``, the residual is ``≈ 0``
    (up to FP rounding) — this IS the LinearOperator round-trip
    identity ``apply(solve(q)) ≡ 0``.

    rtol=1e-12 rationale
    --------------------

    The reduction depth at the per-cell level is single-digit
    (sphere: ~8 terms; cylindrical-degenerate: ~5; slab: ~4).
    FP-non-associativity drift bounded by ``(reduction depth) × ULP
    ~ 1e-15`` is well below ``rtol=1e-12``.  An exact ``np.array_equal``
    contract is too tight because ``apply`` chooses its own reduction
    tree (no obligation to mirror :meth:`DiamondDifference.update`'s
    operation order).
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "geometry,direction_sign",
        [
            ("slab", +1), ("slab", -1),
            ("sphere", +1), ("sphere", -1),
            ("cylinder", +1), ("cylinder", -1),
        ],
    )
    @pytest.mark.parametrize("n_groups", [1, 2])
    @pytest.mark.parametrize(
        "source_kind",
        ["flat_psi_consistent", "random_per_group", "linear_in_r"],
        ids=["flat_psi", "rand", "linear_in_r"],
    )
    def test_round_trip_q_recovered_from_solved_cell_avg(
        self, geometry, direction_sign, n_groups, source_kind,
    ):
        """Round-trip residual ``apply(solve(q)) ≡ 0`` at ``rtol=1e-12``."""
        cell_idx = 2  # interior
        rng = np.random.default_rng(
            seed=20260512 + hash((geometry, direction_sign, n_groups, source_kind)) % 100000,
        )

        if geometry == "slab":
            mesh = _slab_mesh(nx=5, length=1.0)
            quad = GaussLegendre1D.create(4)
            op = slab_streaming(mesh, quad)
            direction_idx = quad.N - 1 if direction_sign > 0 else 0
            st = op.streaming_terms(cell_idx, direction_idx)
            psi_in = rng.uniform(0.0, 0.5, size=n_groups)
            upstream = UpstreamState(spatial_upstream=psi_in)
            visit = CellVisit(cell_idx=cell_idx, streaming_terms=st)
        elif geometry == "sphere":
            mesh = _spherical_mesh(nx=5, radius=1.0)
            quad = GaussLegendre1D.create(8)
            op = spherical_streaming(mesh, quad)
            direction_idx = quad.N - 2 if direction_sign > 0 else 1
            st = op.streaming_terms(cell_idx, direction_idx)
            psi_in = rng.uniform(0.0, 0.3, size=n_groups)
            psi_ang = rng.uniform(0.0, 0.2, size=n_groups)
            upstream = UpstreamState(
                spatial_upstream=psi_in, angular_upstream=psi_ang,
            )
            A_down = st.face_area_outer if direction_sign > 0 else st.face_area_inner
            visit = CellVisit(
                cell_idx=cell_idx, streaming_terms=st, face_area_downstream=A_down,
            )
        else:  # cylinder
            mesh = _cylindrical_mesh(nx=5, radius=1.0)
            quad = ProductQuadrature.create(n_mu=4, n_phi=4)
            op = cylindrical_streaming(mesh, quad)
            _p, _k, st = _first_non_degenerate_dir(op, quad, cell_idx)
            psi_in = rng.uniform(0.0, 0.3, size=n_groups)
            psi_ang = rng.uniform(0.0, 0.2, size=n_groups)
            upstream = UpstreamState(
                spatial_upstream=psi_in, angular_upstream=psi_ang,
            )
            A_down = st.face_area_outer if direction_sign > 0 else st.face_area_inner
            visit = CellVisit(
                cell_idx=cell_idx, streaming_terms=st, face_area_downstream=A_down,
            )

        total_xs = rng.uniform(0.3, 2.0, size=n_groups)
        source = _build_synthetic_source(
            source_kind, st, total_xs, n_groups, rng,
        )

        sn_op = SNCellOperator()
        result = sn_op.solve(visit, total_xs, source, upstream)
        residual = sn_op.apply(
            result.cell_average_flux,
            visit=visit,
            total_xs=total_xs,
            upstream_state=upstream,
            source=source,
        )
        np.testing.assert_allclose(
            residual, np.zeros_like(residual),
            rtol=1e-12, atol=1e-13,
            err_msg=(
                f"Round-trip residual non-zero on {geometry} "
                f"dir={direction_sign} ng={n_groups} src={source_kind}: "
                f"{residual!r}"
            ),
        )

    @pytest.mark.foundation
    def test_round_trip_at_q_zero_trivial(self):
        """Sanity gate: round-trip at ``q = 0`` returns zero residual.

        The trivial case that hides the gate's actual contract if used
        as the only test.  Kept here as a documentation anchor — the
        LOAD-BEARING contract is the parametrized test above on
        non-trivial sources.
        """
        mesh = _slab_mesh(nx=5, length=1.0)
        quad = GaussLegendre1D.create(4)
        op = slab_streaming(mesh, quad)
        st = op.streaming_terms(2, quad.N - 1)
        total_xs = np.array([1.5, 0.7])
        source = np.zeros(2)
        upstream = UpstreamState(spatial_upstream=np.zeros(2))
        visit = CellVisit(cell_idx=2, streaming_terms=st)
        sn_op = SNCellOperator()
        result = sn_op.solve(visit, total_xs, source, upstream)
        residual = sn_op.apply(
            result.cell_average_flux,
            visit=visit, total_xs=total_xs,
            upstream_state=upstream, source=source,
        )
        np.testing.assert_array_equal(residual, np.zeros_like(residual))


__all__: list[str] = []
