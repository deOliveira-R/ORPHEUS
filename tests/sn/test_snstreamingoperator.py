"""Foundation tests for :class:`orpheus.sn.operator.SNStreamingOperator`.

Round 3 of Wave D of the SN reshape campaign (Issue #160).  The
operator is the Wave D capstone — a unified
:class:`~orpheus.numerics.operator.LinearOperator` for the streaming-
collision operator :math:`L = \\Omega\\cdot\\nabla + \\Sigma_t` with three
capabilities:

* :meth:`SNStreamingOperator.apply` — matrix-free forward action, the
  symmetric closure carried by the historical
  ``transport_operator_matvec_*`` BiCGSTAB FD operator.
* :meth:`SNStreamingOperator.solve` — :math:`L^{-1}\\,q` via the
  Wave D Round 2 unified sweep.
* :meth:`SNStreamingOperator.apply_transpose` — adjoint action via the
  explicit transpose of the dense matrix assembled by probing
  :meth:`apply`.

The test suite below pins three structural-correctness claims:

1. **Capability advertisement** — the operator exposes
   ``{"apply", "solve", "apply_transpose"}`` and is recognised by the
   :class:`~orpheus.numerics.operator.LinearOperator` Protocol.

2. **Bit-identical extraction** — :meth:`apply` produces ``np.array_equal``-
   identical results to the underlying
   ``transport_operator_matvec_*`` functions; :meth:`solve` produces
   ``np.array_equal``-identical results to a direct
   :func:`transport_sweep` call on the same arguments.  This is the
   load-bearing claim that the extraction is correct (the math is moved
   verbatim, not rewritten).

3. **Reciprocity** — :math:`\\langle L\\,\\psi,\\,\\varphi\\rangle =
   \\langle\\psi,\\, L^*\\,\\varphi\\rangle` to round-off
   (``rtol=1e-12``) on synthetic ``(\\psi, \\varphi)`` pairs across
   slab + spherical + cylindrical geometries.

These are foundation-tagged because they pin **software-contract**
claims about the operator API rather than transport-equation
identities; the L1 transport math is verified transitively via the
existing MMS suite (:file:`tests/sn/l1_analytical/`) and the regression
snapshots (:file:`tests/sn/regression/snapshots/`).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperator,
)
from orpheus.sn import SNStreamingOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map,
    build_equation_map_spherical,
    build_equation_map_cylindrical,
    transport_operator_matvec,
    transport_operator_matvec_spherical,
    transport_operator_matvec_cylindrical,
)
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from orpheus.sn.sweep import transport_sweep

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — small problems sized for the dense-transpose probe
# ═══════════════════════════════════════════════════════════════════════
#
# Sizing rationale.  The dense-matrix probe in
# :meth:`SNStreamingOperator._ensure_dense_matrix` is :math:`O(n^2)` in
# both time and space.  Reciprocity tests need the matrix only to verify
# the operator's transpose; problem sizes here are deliberately small
# (``n_unknowns ~ 30-200``) so the dense build is well under a second.
# Larger problems would not strengthen the reciprocity claim — the
# operator's transpose is exact by construction at any size.

def _slab_mesh(nx: int = 4, length: float = 1.0) -> SNMesh:
    """Slab Mesh1D + GL N=4 quadrature, vacuum BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad)


def _spherical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Spherical Mesh1D + GL N=4 quadrature, reflective inner / vacuum outer."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad)


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Cylindrical Mesh1D + Level-Symmetric SN-4 quadrature."""
    # Level-symmetric needs at least 2 levels; SN-4 gives 12 ordinates.
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad)


def _sig_t_uniform(sn_mesh: SNMesh, ng: int = 2, value: float = 0.5) -> np.ndarray:
    """Total cross-section, uniform across all cells/groups."""
    return value * np.ones((sn_mesh.nx, sn_mesh.ny, ng))


# ═══════════════════════════════════════════════════════════════════════
# 1. Capability advertisement
# ═══════════════════════════════════════════════════════════════════════

def test_capabilities_full():
    """SNStreamingOperator advertises apply, solve, apply_transpose."""
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh)
    op = SNStreamingOperator(sn_mesh, sig_t)

    assert op.capabilities == frozenset(
        {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
    )


def test_satisfies_linear_operator_protocol():
    """SNStreamingOperator instances satisfy the runtime-checkable Protocol."""
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh)
    op = SNStreamingOperator(sn_mesh, sig_t)

    # The Protocol uses runtime_checkable; isinstance is the canonical check.
    assert isinstance(op, LinearOperator)


def test_capabilities_is_frozenset():
    """capabilities is hashable / immutable (frozenset)."""
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh)
    op = SNStreamingOperator(sn_mesh, sig_t)

    assert isinstance(op.capabilities, frozenset)
    # Hashability is the reason a frozenset (rather than a plain set) is used.
    hash(op.capabilities)


# ═══════════════════════════════════════════════════════════════════════
# 2. Bit-identical extraction — apply matches transport_operator_matvec_*
# ═══════════════════════════════════════════════════════════════════════
#
# These tests pin the load-bearing claim that
# :meth:`SNStreamingOperator.apply` is a structural extraction of the
# legacy ``transport_operator_matvec_*`` functions, not a re-derivation.
# Disagreement here would mean the SNStreamingOperator math has drifted
# from the BiCGSTAB-consumed FD operator math — a campaign-blocking
# regression.

def test_apply_slab_bit_identical_to_legacy():
    """SNStreamingOperator.apply == transport_operator_matvec on slab.

    Wave E Round 3 update: the legacy call now threads the mesh BCs
    through to ``transport_operator_matvec`` so that
    :class:`SNStreamingOperator.apply` (which always reads BCs from
    the mesh) and the legacy call match.  Pre-Round 3 the legacy call
    hard-coded specular reflective fills, which silently agreed with
    ``apply`` only because the old ``apply`` also hard-coded them.
    """
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    # Build the same EquationMap the legacy BiCGSTAB path would build.
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, 2
    eq_map = build_equation_map(nx, ny, sn_mesh.quad, ng)

    rng = np.random.default_rng(seed=42)
    psi = rng.standard_normal(eq_map.n_unknowns)

    legacy = transport_operator_matvec(
        psi, eq_map, sn_mesh.quad, sig_t,
        nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )
    via_op = op.apply(psi)

    assert np.array_equal(legacy, via_op), (
        "SNStreamingOperator.apply diverged from transport_operator_matvec "
        "on slab: max diff = {:.3e}".format(np.max(np.abs(legacy - via_op)))
    )


def test_apply_spherical_bit_identical_to_legacy():
    """SNStreamingOperator.apply == transport_operator_matvec_spherical.

    Wave E Round 3 update: legacy call threads ``sn_mesh.bc_right`` so
    the test exercises the BC-faithful Round 3 plumbing.
    """
    sn_mesh = _spherical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)

    rng = np.random.default_rng(seed=43)
    psi = rng.standard_normal(eq_map.n_unknowns)

    reduced = sn_mesh.reduced
    legacy = transport_operator_matvec_spherical(
        psi, eq_map, sn_mesh.quad, sig_t,
        nx, ng,
        reduced.face_areas,
        sn_mesh.volumes,
        reduced.alpha_half,
        reduced.redist_dAw,
        reduced.tau_mm,
        bc_outer=sn_mesh.bc_right,
    )
    via_op = op.apply(psi)

    assert np.array_equal(legacy, via_op), (
        "SNStreamingOperator.apply diverged from "
        "transport_operator_matvec_spherical: max diff = {:.3e}".format(
            np.max(np.abs(legacy - via_op))
        )
    )


def test_apply_cylindrical_bit_identical_to_legacy():
    """SNStreamingOperator.apply == transport_operator_matvec_cylindrical.

    Wave E Round 3 update: legacy call threads ``sn_mesh.bc_right``.
    """
    sn_mesh = _cylindrical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_cylindrical(nx, sn_mesh.quad, ng)

    rng = np.random.default_rng(seed=44)
    psi = rng.standard_normal(eq_map.n_unknowns)

    reduced = sn_mesh.reduced
    legacy = transport_operator_matvec_cylindrical(
        psi, eq_map, sn_mesh.quad, sig_t,
        nx, ng,
        reduced.face_areas,
        sn_mesh.volumes,
        reduced.alpha_per_level,
        reduced.redist_dAw_per_level,
        reduced.tau_mm_per_level,
        bc_outer=sn_mesh.bc_right,
    )
    via_op = op.apply(psi)

    assert np.array_equal(legacy, via_op), (
        "SNStreamingOperator.apply diverged from "
        "transport_operator_matvec_cylindrical: max diff = {:.3e}".format(
            np.max(np.abs(legacy - via_op))
        )
    )


def test_apply_2d_cartesian_bit_identical_to_legacy():
    """SNStreamingOperator.apply == transport_operator_matvec on 2-D Cartesian.

    The 2-D Cartesian path exercises the full
    :func:`build_equation_map` filter (including z-hemisphere and
    y-axis reflective BCs) that 1-D slab does not.

    Wave E Round 3 update: legacy call threads all 4 BCs from the mesh.
    """
    nx, ny = 3, 3
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    # Keep N small so dense-matrix probe stays fast.
    from orpheus.sn.quadrature import LebedevSphere
    quad = LebedevSphere.create(order=5)
    sn_mesh = SNMesh(mesh, quad)
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    ng = 2
    eq_map = build_equation_map(nx, ny, quad, ng)

    rng = np.random.default_rng(seed=45)
    psi = rng.standard_normal(eq_map.n_unknowns)

    legacy = transport_operator_matvec(
        psi, eq_map, quad, sig_t,
        nx, ny, ng, sn_mesh.dx, sn_mesh.dy,
        bc_xmin=sn_mesh.bc_xmin, bc_xmax=sn_mesh.bc_xmax,
        bc_ymin=sn_mesh.bc_ymin, bc_ymax=sn_mesh.bc_ymax,
    )
    via_op = op.apply(psi)

    assert np.array_equal(legacy, via_op), (
        "SNStreamingOperator.apply diverged from "
        "transport_operator_matvec on 2-D: max diff = {:.3e}".format(
            np.max(np.abs(legacy - via_op))
        )
    )


# ═══════════════════════════════════════════════════════════════════════
# 3. Bit-identical solve — solve matches transport_sweep
# ═══════════════════════════════════════════════════════════════════════

def test_solve_slab_bit_identical_to_transport_sweep():
    """SNStreamingOperator.solve produces transport_sweep's exact output."""
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=10)
    Q = rng.standard_normal((sn_mesh.nx, sn_mesh.ny, 2))

    # Use independent psi_bc dicts so neither path mutates the other.
    af_legacy, sf_legacy = transport_sweep(Q, sig_t, sn_mesh, {}, None)
    af_op, sf_op = op.solve(Q, psi_bc={}, Q_aniso=None)

    assert np.array_equal(af_legacy, af_op), (
        "SNStreamingOperator.solve angular-flux drift on slab: "
        "max diff = {:.3e}".format(np.max(np.abs(af_legacy - af_op)))
    )
    assert np.array_equal(sf_legacy, sf_op), (
        "SNStreamingOperator.solve scalar-flux drift on slab"
    )


def test_solve_spherical_bit_identical_to_transport_sweep():
    """SNStreamingOperator.solve == transport_sweep on spherical."""
    sn_mesh = _spherical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=11)
    Q = rng.standard_normal((sn_mesh.nx, sn_mesh.ny, 2))

    af_legacy, sf_legacy = transport_sweep(Q, sig_t, sn_mesh, {}, None)
    af_op, sf_op = op.solve(Q, psi_bc={}, Q_aniso=None)

    assert np.array_equal(af_legacy, af_op)
    assert np.array_equal(sf_legacy, sf_op)


def test_solve_cylindrical_bit_identical_to_transport_sweep():
    """SNStreamingOperator.solve == transport_sweep on cylindrical."""
    sn_mesh = _cylindrical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=12)
    Q = rng.standard_normal((sn_mesh.nx, sn_mesh.ny, 2))

    af_legacy, sf_legacy = transport_sweep(Q, sig_t, sn_mesh, {}, None)
    af_op, sf_op = op.solve(Q, psi_bc={}, Q_aniso=None)

    assert np.array_equal(af_legacy, af_op)
    assert np.array_equal(sf_legacy, sf_op)


# ═══════════════════════════════════════════════════════════════════════
# 4. Reciprocity — <L·ψ, φ> == <ψ, L*·φ> to round-off
# ═══════════════════════════════════════════════════════════════════════
#
# The reciprocity identity ⟨Lψ, φ⟩ = ⟨ψ, L*φ⟩ is the defining property
# of the operator-transpose pairing (Lewis & Miller §10).  In ORPHEUS
# the apply_transpose is built from the explicit transpose of the dense
# matrix assembled by probing apply, so the identity holds **by
# construction**: any failure here would mean either
#
#   (a) :meth:`apply` is non-linear, OR
#   (b) the dense-assembly probe code does not faithfully assemble
#       the matrix of :meth:`apply`.
#
# Both are catastrophic operator-correctness failures.  The reciprocity
# test thus serves as a forward-correctness gate on the apply path.

@pytest.mark.verifies("sn-streaming-reciprocity")
@pytest.mark.parametrize(
    "geometry,seed",
    [
        ("slab", 100),
        ("spherical", 200),
        ("cylindrical", 300),
    ],
)
def test_reciprocity_round_off(geometry: str, seed: int):
    """⟨L·ψ, φ⟩ == ⟨ψ, L*·φ⟩ to round-off across all three geometries."""
    if geometry == "slab":
        sn_mesh = _slab_mesh()
    elif geometry == "spherical":
        sn_mesh = _spherical_mesh()
    elif geometry == "cylindrical":
        sn_mesh = _cylindrical_mesh()
    else:
        raise ValueError(f"unknown geometry {geometry!r}")

    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=seed)
    psi = rng.standard_normal(op.n_unknowns)
    phi = rng.standard_normal(op.n_unknowns)

    inner_forward = float(op.apply(psi) @ phi)
    inner_adjoint = float(psi @ op.apply_transpose(phi))

    # Reciprocity should hold to round-off.  rtol=1e-12 absorbs the
    # float-summation drift; atol=1e-13 catches near-zero pathologies
    # (which the random ψ/φ above don't produce, but defend against).
    assert np.isclose(
        inner_forward, inner_adjoint, rtol=1e-12, atol=1e-13,
    ), (
        f"Reciprocity violated on {geometry}: "
        f"<Lψ,φ> = {inner_forward:.15e}, "
        f"<ψ,L*φ> = {inner_adjoint:.15e}, "
        f"rel diff = {abs(inner_forward - inner_adjoint)/max(abs(inner_forward), abs(inner_adjoint)):.3e}"
    )


@pytest.mark.verifies("sn-streaming-reciprocity")
@pytest.mark.parametrize(
    "geometry,seed",
    [
        ("slab", 1000),
        ("spherical", 2000),
        ("cylindrical", 3000),
    ],
)
def test_reciprocity_multiple_pairs(geometry: str, seed: int):
    """Reciprocity holds for **multiple** independent (ψ, φ) pairs.

    A single-pair test could conceivably pass by accident if the inner
    product happens to land on a degenerate axis.  This test draws 10
    independent pairs and gates them all.
    """
    if geometry == "slab":
        sn_mesh = _slab_mesh()
    elif geometry == "spherical":
        sn_mesh = _spherical_mesh()
    elif geometry == "cylindrical":
        sn_mesh = _cylindrical_mesh()
    else:
        raise ValueError(f"unknown geometry {geometry!r}")

    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=seed)
    for trial in range(10):
        psi = rng.standard_normal(op.n_unknowns)
        phi = rng.standard_normal(op.n_unknowns)

        inner_forward = float(op.apply(psi) @ phi)
        inner_adjoint = float(psi @ op.apply_transpose(phi))

        assert np.isclose(
            inner_forward, inner_adjoint, rtol=1e-12, atol=1e-13,
        ), (
            f"Reciprocity violated on {geometry} trial {trial}: "
            f"<Lψ,φ> = {inner_forward:.15e}, "
            f"<ψ,L*φ> = {inner_adjoint:.15e}"
        )


# ═══════════════════════════════════════════════════════════════════════
# 5. Apply linearity — apply(αψ + βφ) == α·apply(ψ) + β·apply(φ)
# ═══════════════════════════════════════════════════════════════════════
#
# Linearity of :meth:`apply` is a precondition for the dense-matrix
# probe to be valid.  Test it explicitly as a defence-in-depth gate.

@pytest.mark.parametrize(
    "geometry",
    ["slab", "spherical", "cylindrical"],
)
def test_apply_is_linear(geometry: str):
    """apply(αψ + βφ) == α·apply(ψ) + β·apply(φ)."""
    if geometry == "slab":
        sn_mesh = _slab_mesh()
    elif geometry == "spherical":
        sn_mesh = _spherical_mesh()
    elif geometry == "cylindrical":
        sn_mesh = _cylindrical_mesh()
    else:
        raise ValueError(f"unknown geometry {geometry!r}")

    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    rng = np.random.default_rng(seed=4242)
    psi = rng.standard_normal(op.n_unknowns)
    phi = rng.standard_normal(op.n_unknowns)
    alpha, beta = 1.7, -2.3

    lhs = op.apply(alpha * psi + beta * phi)
    rhs = alpha * op.apply(psi) + beta * op.apply(phi)

    np.testing.assert_allclose(lhs, rhs, rtol=1e-13, atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# 6. n_unknowns and EquationMap dispatch
# ═══════════════════════════════════════════════════════════════════════

def test_n_unknowns_matches_eq_map_slab():
    """SNStreamingOperator.n_unknowns == build_equation_map(...).n_unknowns."""
    sn_mesh = _slab_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    expected = build_equation_map(
        sn_mesh.nx, sn_mesh.ny, sn_mesh.quad, 2,
    ).n_unknowns
    assert op.n_unknowns == expected


def test_n_unknowns_matches_eq_map_spherical():
    """SNStreamingOperator.n_unknowns matches the spherical EquationMap."""
    sn_mesh = _spherical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    expected = build_equation_map_spherical(
        sn_mesh.nx, sn_mesh.quad, 2,
    ).n_unknowns
    assert op.n_unknowns == expected


def test_n_unknowns_matches_eq_map_cylindrical():
    """SNStreamingOperator.n_unknowns matches the cylindrical EquationMap."""
    sn_mesh = _cylindrical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)

    expected = build_equation_map_cylindrical(
        sn_mesh.nx, sn_mesh.quad, 2,
    ).n_unknowns
    assert op.n_unknowns == expected


# ═══════════════════════════════════════════════════════════════════════
# Issue #168 Phase A — Defect 1 + Defect 2 fixes
# ═══════════════════════════════════════════════════════════════════════
#
# These tests pin the Phase-A boundary-closure behaviour:
#
# * Defect 1 (outer-face cell-centre truncation): ``apply`` at the
#   outer boundary cell uses :class:`DDExtrapolation` (default)
#   producing :math:`1.5 \\psi_{N-1} - 0.5 \\psi_{N-2}` rather than
#   the cell-centre value :math:`\\psi_{N-1}`.
# * Defect 2 (storage-conflation between cell-centre and BC face
#   value): the new ``solution_to_angular_flux_spherical`` returns
#   ``(fi, boundary_face_flux)`` with ``fi`` as PURE cell-centre
#   storage (no BC face values clobbering ``fi[..., -1, 0]``).  The
#   matvec at i=N-2 reads cell-centre vs cell-centre on the
#   interior-face stencil ``0.5*(fi[..., N-2, 0] + fi[..., N-1, 0])``,
#   uncontaminated by BC face values.


def test_solution_to_angular_flux_spherical_returns_single_array():
    """Issue #168 Phase C: returns a single ``fi`` array.

    Phase A returned ``(fi, boundary_face_flux)``; Phase C retires
    the boundary_face_flux companion array — the BC trace law is
    applied INSIDE the matvec body, consuming WDD-propagated outflow
    face values (not the cell-centre slots). The Phase A two-storage
    design is subsumed by the sweep-frame architecture.
    """
    from orpheus.sn.operator import solution_to_angular_flux_spherical

    sn_mesh = _spherical_mesh()
    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)
    psi = np.ones(eq_map.n_unknowns)

    fi = solution_to_angular_flux_spherical(
        psi, eq_map, sn_mesh.quad, nx, ng,
    )
    assert isinstance(fi, np.ndarray)
    assert fi.shape == (ng, sn_mesh.quad.N, nx, 1)


def test_solution_to_angular_flux_spherical_inward_extension():
    """Phase C inward analytical extension: ``fi[..., n_inward, -1, 0]``
    equals the reflected-partner cell-centre value.

    The eq_map excludes inward-at-outer-boundary slots (BC determines
    them); ``solution_to_angular_flux_spherical`` fills them with the
    reflected partner's cell-centre so the WDD recurrence on flat ψ
    preserves the per-ordinate flat-flux invariant. The BC trace law
    is applied separately at the boundary edge of the WDD chain
    inside the matvec.
    """
    from orpheus.sn.operator import solution_to_angular_flux_spherical

    sn_mesh = _spherical_mesh()
    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)
    rng = np.random.default_rng(seed=911)
    psi = rng.standard_normal(eq_map.n_unknowns)

    fi = solution_to_angular_flux_spherical(
        psi, eq_map, sn_mesh.quad, nx, ng,
    )

    ref_x = sn_mesh.quad.reflection_index("x")
    for n in range(sn_mesh.quad.N):
        if sn_mesh.quad.mu_x[n] < -1e-15:
            np.testing.assert_array_equal(
                fi[:, n, nx - 1, 0],
                fi[:, ref_x[n], nx - 1, 0],
                err_msg=f"Phase C inward extension: cell-centre at "
                        f"i=N-1 for inward ordinate {n} must equal "
                        f"the reflected partner's cell-centre.",
            )


def test_apply_spherical_wdd_recurrence_deterministic():
    """Phase C: the WDD sweep-frame matvec is bit-stable (no per-call
    nondeterminism such as random ordering or floating-point sums in
    different orders).

    This is Gate 1.2 (foundation, plan §5): the sweep-frame matvec
    is structurally bit-deterministic. Pinned via np.array_equal on
    repeated apply calls — replaces the Phase A
    ``test_apply_spherical_outer_face_uses_dd_extrapolation`` which
    pinned the now-retired BoundaryFaceFlux Protocol.
    """
    from orpheus.sn.operator import transport_operator_matvec_spherical

    sn_mesh = _spherical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)
    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)

    rng = np.random.default_rng(seed=611)
    psi = rng.standard_normal(eq_map.n_unknowns)

    via_op_1 = op.apply(psi)
    via_op_2 = op.apply(psi)
    np.testing.assert_array_equal(via_op_1, via_op_2)

    # Direct invocation through ``transport_operator_matvec_spherical``
    # must match the dispatcher-via-op path bit-for-bit.
    reduced = sn_mesh.reduced
    via_direct = transport_operator_matvec_spherical(
        psi, eq_map, sn_mesh.quad, sig_t,
        nx, ng,
        reduced.face_areas, sn_mesh.volumes,
        reduced.alpha_half, reduced.redist_dAw, reduced.tau_mm,
        sn_mesh=sn_mesh,
        bc_outer=sn_mesh.bc_right,
        pole_angular_closure=sn_mesh.pole_angular_closure,
    )
    np.testing.assert_array_equal(via_op_1, via_direct)


def test_apply_spherical_constant_flux_yields_zero_collisionless_reflective():
    """Constant ψ + zero Σ_t + REFLECTIVE BC: per-ordinate flat-flux
    invariant holds under Phase C sweep-frame matvec.

    The Phase C rewrite (sweep-frame WDD propagation + BC trace law
    at the boundary edge) preserves the per-ordinate flat-flux
    invariant on a reflective sphere under the default
    LegacyTauSymmetricInterpolation closure. This is Gate 1.1
    parametrised at sphere-legacy-Σ_t=0 (see test_phase_c_gates.py).

    Under vacuum BC the invariant does NOT hold (the BC physically
    removes the inflow → flat ψ cannot persist) — see
    test_apply_spherical_vacuum_bc_residual_is_consistent for the
    vacuum-BC test.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.quadrature import GaussLegendre1D

    nx = 4
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    sn_mesh = SNMesh(mesh, quad)
    ng = 2
    sig_t = np.zeros((nx, sn_mesh.ny, ng))
    op = SNStreamingOperator(sn_mesh, sig_t)
    eq_map = build_equation_map_spherical(nx, quad, ng)

    psi = np.full(eq_map.n_unknowns, 3.14159)
    result = op.apply(psi)
    np.testing.assert_allclose(
        result, np.zeros_like(result),
        atol=1e-12,
        err_msg="Per-ordinate flat-flux invariant violated on "
                "reflective sphere under Phase C sweep-frame matvec.",
    )


def test_apply_spherical_vacuum_bc_residual_is_consistent():
    """Phase C: under vacuum BC, the apply matvec consumes the
    realised BC operator on the WDD-propagated outflow face values.

    On a flat cell-centre ψ with vacuum BC, the residual is NOT
    zero (the BC physically removes inflow → flat ψ cannot persist
    as a transport solution). The Phase C contract is structural:
    the BC trace law's input is the WDD-propagated outflow face
    value, not a cell-centre approximation. Pinned via determinism
    (apply is bit-stable; replaces the Phase A constant-flux-yields-
    zero invariant which relied on the now-retired Defect-2
    cell-centre / BC-face-value separation).
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.quadrature import GaussLegendre1D

    nx = 4
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    sn_mesh = SNMesh(mesh, quad)
    ng = 2
    sig_t = np.zeros((nx, sn_mesh.ny, ng))
    op = SNStreamingOperator(sn_mesh, sig_t)
    eq_map = build_equation_map_spherical(nx, quad, ng)

    # apply(0) = 0 under any BC kind (linearity).
    zero = np.zeros(eq_map.n_unknowns)
    np.testing.assert_array_equal(op.apply(zero), np.zeros_like(zero))
    # apply(constant) ≠ 0 under vacuum BC (BC removes inflow).
    psi = np.full(eq_map.n_unknowns, 1.0)
    result = op.apply(psi)
    assert np.max(np.abs(result)) > 1e-3, (
        "Vacuum BC + flat cell-centre ψ should produce a nonzero "
        "residual (BC physically removes the inflow → flat ψ "
        "cannot persist); got effectively zero, indicating the BC "
        "trace law is not being applied at the boundary edge."
    )
    # The matvec is deterministic (bit-stable).
    np.testing.assert_array_equal(result, op.apply(psi))


def test_apply_spherical_constant_flux_under_morel_montry_canonical_form():
    """Phase B canonical Hébert §3.9.4 form: per-ordinate flat-flux
    balance is **deliberately broken** in exchange for asymptotic
    :math:`\\mathcal{O}(h^2)` accuracy on angularly-varying ψ.  The
    surviving invariant is **angular-integrated**: :math:`\\sum_n
    w_n \\cdot R_{n,i,g} = 0` by α-telescoping (Hébert Eq. 3.420).
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.quadrature import GaussLegendre1D
    from orpheus.sn.spatial.pole_angular_closure import (
        MorelMontryAngularSweep,
    )

    # Phase C: use a reflective sphere so the BC trace law preserves
    # the flat-ψ inflow trace; the angular-integrated invariant then
    # holds globally (not just at interior cells, as under the
    # pre-Phase-C vacuum-BC fixture).
    nx = 4
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    sn_mesh = SNMesh(mesh, quad, pole_angular_closure=MorelMontryAngularSweep())
    ng = 2
    sig_t = np.zeros((nx, sn_mesh.ny, ng))
    op = SNStreamingOperator(sn_mesh, sig_t)
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)

    psi = np.full(eq_map.n_unknowns, 3.14159)
    result_packed = op.apply(psi)
    result = result_packed.reshape(ng, eq_map.n_eq, order='F')
    weights = sn_mesh.quad.weights

    # ── Per-ordinate balance: NOT preserved (deliberately) ────────
    # On constant ψ + zero Σ_t, the Hébert canonical recurrence makes
    # the half-angle face fluxes oscillate (0, 2c, 0, 2c, ...), so
    # the per-ordinate redistribution does NOT cancel the per-ordinate
    # streaming.  This is the canonical Hébert form's structural
    # signature.
    assert np.max(np.abs(result)) > 1.0, (
        "MorelMontryAngularSweep should produce non-zero per-ordinate "
        "residuals on constant ψ (the canonical Hébert form does NOT "
        "preserve per-ordinate flat-flux balance — see Issue #168 "
        "Phase B closeout)."
    )

    # ── Angular-integrated balance: preserved by α-telescoping ────
    # Σ_n w_n × R_{n,i,g} = 0 at INTERIOR cells (where all ordinates
    # are unknowns).  At i=nx-1 the eq_map skips inward ordinates
    # (BC determines them), so the partial sum doesn't telescope.
    integrated = np.zeros((ng, nx))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        integrated[:, i] += weights[n] * result[:, k]
    np.testing.assert_allclose(
        integrated[:, :-1], np.zeros((ng, nx - 1)),
        atol=1e-12,
        err_msg="Angular-integrated flat-flux invariant violated "
                "at interior cells under MorelMontryAngularSweep.",
    )


def test_apply_cylindrical_wdd_recurrence_deterministic():
    """Cylindrical Phase C sweep-frame matvec is bit-stable.

    Replaces the Phase A
    ``test_apply_cylindrical_outer_face_uses_dd_extrapolation``
    which pinned the now-retired BoundaryFaceFlux Protocol.
    """
    from orpheus.sn.operator import transport_operator_matvec_cylindrical

    sn_mesh = _cylindrical_mesh()
    sig_t = _sig_t_uniform(sn_mesh, ng=2)
    op = SNStreamingOperator(sn_mesh, sig_t)
    nx, ng = sn_mesh.nx, 2
    eq_map = build_equation_map_cylindrical(nx, sn_mesh.quad, ng)

    rng = np.random.default_rng(seed=612)
    psi = rng.standard_normal(eq_map.n_unknowns)

    via_op_1 = op.apply(psi)
    via_op_2 = op.apply(psi)
    np.testing.assert_array_equal(via_op_1, via_op_2)

    reduced = sn_mesh.reduced
    via_direct = transport_operator_matvec_cylindrical(
        psi, eq_map, sn_mesh.quad, sig_t,
        nx, ng,
        reduced.face_areas, sn_mesh.volumes,
        reduced.alpha_per_level, reduced.redist_dAw_per_level,
        reduced.tau_mm_per_level,
        sn_mesh=sn_mesh,
        bc_outer=sn_mesh.bc_right,
        pole_angular_closure=sn_mesh.pole_angular_closure,
    )
    np.testing.assert_array_equal(via_op_1, via_direct)


def test_snmesh_no_longer_accepts_boundary_face_flux_kwarg():
    """Phase C: ``boundary_face_flux`` kwarg + field have been retired.

    Constructing SNMesh with ``boundary_face_flux=...`` MUST raise
    TypeError; the attribute MUST be absent from the resulting
    instance. The Phase A Protocol is subsumed by the sweep-frame
    matvec's WDD propagation + BC trace law architecture.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.quadrature import GaussLegendre1D

    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    sn = SNMesh(mesh, quad)
    assert not hasattr(sn, "boundary_face_flux"), (
        "Phase C retired SNMesh.boundary_face_flux; attribute is "
        "back on the instance."
    )
    # Passing the kwarg raises TypeError (the constructor signature
    # no longer accepts it).
    import pytest as _pytest
    with _pytest.raises(TypeError):
        SNMesh(mesh, quad, boundary_face_flux=object())
