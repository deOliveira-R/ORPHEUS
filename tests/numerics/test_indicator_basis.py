r"""Intrinsic-law tests for :class:`~orpheus.numerics.basis.IndicatorBasis`.

The piecewise-constant (P0) cell-indicator basis is the trial side of the
spatial / energy homogenisation :class:`~orpheus.numerics.frame.Frame`.  These
``foundation`` tests pin its DEFINING laws (independent of any consumer):

* **partition of unity** — every point lands in exactly one cell (one-hot rows),
  and the table is bit-exactly the cell-membership matrix;
* **diagonal Gram** — the indicators have disjoint supports, so the discrete
  Gram against any measure is :math:`\mathrm{diag}(m_R)` (region masses);
* **projector idempotence** — :math:`R \circ G^{-1} \circ M` is a projector, and
  it is the identity on fields already piecewise-constant on the cells;
* **Moore–Penrose null** — an empty / zero-mass region normalises to zero, not
  ``nan`` (the homogenisation zero-flux-region behaviour);
* **n-D layout** — the per-axis ``ravel_multi_index("ij")`` flattening agrees
  with ``DiscreteMeasure``'s ``meshgrid("ij").ravel()`` cell ordering, so a
  future n-D ``homogenize`` inherits the correct layout by construction;
* **role separation (negative)** — a mesh is NOT a :class:`Basis`; it *yields*
  one.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.basis import Basis, IndicatorBasis
from orpheus.numerics.measure import DiscreteMeasure

pytestmark = [pytest.mark.foundation]


# ── fixtures ──────────────────────────────────────────────────────────
_FINE_EDGES = np.linspace(0.0, 4.0, 9)                  # 8 fine cells
_FINE_CENTERS = 0.5 * (_FINE_EDGES[:-1] + _FINE_EDGES[1:])
_FINE_VOLS = np.diff(_FINE_EDGES)
_COARSE_EDGES = np.array([0.0, 2.0, 4.0])               # 2 coarse cells


def _basis(coarse_edges=_COARSE_EDGES):
    return IndicatorBasis(edges_per_axis=(np.asarray(coarse_edges, float),))


def _measure(centers=_FINE_CENTERS, vols=_FINE_VOLS):
    return DiscreteMeasure(nodes=centers, weights=vols, support="spatial_R1")


def _independent_membership(coarse_edges, centers):
    """Structurally-independent one-hot membership (own searchsorted)."""
    n_coarse = coarse_edges.size - 1
    cof = np.clip(
        np.searchsorted(coarse_edges, centers, side="right") - 1, 0, n_coarse - 1
    )
    table = np.zeros((centers.size, n_coarse))
    table[np.arange(centers.size), cof] = 1.0
    return table


# ── partition of unity + membership-onehot ────────────────────────────
def test_partition_of_unity():
    """Every fine point in exactly one coarse cell — rows are 0/1 and sum to 1."""
    table = _basis().evaluate(_FINE_CENTERS)
    np.testing.assert_array_equal(table.sum(axis=1), np.ones(_FINE_CENTERS.size))
    assert set(np.unique(table)) <= {0.0, 1.0}


def test_table_is_membership_onehot_bit_exact():
    """evaluate() equals the independent membership matrix — BIT-exact (0/1 ints)."""
    table = _basis().evaluate(_FINE_CENTERS)
    legacy = _independent_membership(_COARSE_EDGES, _FINE_CENTERS)
    np.testing.assert_array_equal(table, legacy)


def test_evaluate_accepts_1d_and_2d_points():
    """A 1-D node array and its (N,1) column form give the same table."""
    b = _basis()
    np.testing.assert_array_equal(
        b.evaluate(_FINE_CENTERS), b.evaluate(_FINE_CENTERS[:, None]),
    )


# ── diagonal Gram ─────────────────────────────────────────────────────
def test_diagonal_gram_disjoint_supports():
    """mass_matrix is diagonal = region masses (indicators have disjoint support)."""
    b, mu = _basis(), _measure()
    gram = b.mass_matrix(mu)
    np.testing.assert_allclose(gram, np.diag(np.diag(gram)), atol=0.0)
    region_mass = _independent_membership(_COARSE_EDGES, _FINE_CENTERS).T @ _FINE_VOLS
    np.testing.assert_allclose(np.diag(gram), region_mass)


# ── projector idempotence ─────────────────────────────────────────────
def _projector(b, table, weights):
    """Π = R ∘ G⁻¹ ∘ M with the region-mass Gram (the homogenisation projector)."""
    region_mass = b.analyze(np.ones(table.shape[0]), table, weights)
    space = dataclasses.replace(b.space, inner_product_weights=region_mass)

    def project(fine_field):
        coeff = space.apply_inverse_metric(b.analyze(fine_field, table, weights))
        return b.reconstruct(coeff, table)

    return project


def test_projector_is_idempotent():
    """R∘G⁻¹∘M is a projector: applying it twice equals applying it once."""
    b, mu = _basis(), _measure()
    table = b.evaluate(mu.nodes)
    project = _projector(b, table, mu.weights)
    rng_field = np.sin(_FINE_CENTERS) + 2.0           # arbitrary fine field
    once = project(rng_field)
    twice = project(once)
    np.testing.assert_allclose(twice, once, rtol=0, atol=1e-13)


def test_projector_is_identity_on_piecewise_constant():
    """A field already piecewise-constant on the coarse cells is unchanged."""
    b, mu = _basis(), _measure()
    table = b.evaluate(mu.nodes)
    project = _projector(b, table, mu.weights)
    coarse_field = b.reconstruct(np.array([3.0, 7.0]), table)   # constant per cell
    np.testing.assert_allclose(project(coarse_field), coarse_field, atol=1e-13)


def test_reconstruct_is_piecewise_constant():
    """reconstruct broadcasts each cell's coefficient onto exactly its points."""
    b = _basis()
    table = b.evaluate(_FINE_CENTERS)
    recon = b.reconstruct(np.array([3.0, 7.0]), table)
    legacy = _independent_membership(_COARSE_EDGES, _FINE_CENTERS)
    expected = legacy @ np.array([3.0, 7.0])
    np.testing.assert_array_equal(recon, expected)


# ── Moore–Penrose null (empty / zero-mass region) ─────────────────────
def test_empty_region_normalises_to_zero_not_nan():
    """A coarse cell containing no points has zero Gram mass → 0, not nan/inf."""
    # coarse cell 0 = [0, 0.1] contains no fine centre (first centre is 0.25).
    b = _basis(np.array([0.0, 0.1, 4.0]))
    mu = _measure()
    table = b.evaluate(mu.nodes)
    region_mass = b.analyze(np.ones(mu.nodes.size), table, mu.weights)
    assert region_mass[0] == 0.0                       # the empty region
    space = dataclasses.replace(b.space, inner_product_weights=region_mass)
    rate = b.analyze(np.ones(mu.nodes.size), table, mu.weights)
    normalised = space.apply_inverse_metric(rate)
    assert np.isfinite(normalised).all()
    assert normalised[0] == 0.0                        # Moore–Penrose null, not nan


# ── n-D layout pin (the future-n-D-homogenize guard) ──────────────────
def test_indicator_2d_matches_volume_measure_layout():
    """A 2-axis IndicatorBasis buckets 2-D cell-centre points via
    ravel_multi_index("ij"), agreeing with meshgrid("ij").ravel() — so the
    membership of the n-D cell-centre grid is the IDENTITY matrix."""
    ex = np.array([0.0, 1.0, 2.0])                     # 2 cells in x
    ey = np.array([0.0, 1.0, 2.0, 3.0])                # 3 cells in y → 6 cells
    b = IndicatorBasis(edges_per_axis=(ex, ey))
    cx = 0.5 * (ex[:-1] + ex[1:])
    cy = 0.5 * (ey[:-1] + ey[1:])
    grid_x, grid_y = np.meshgrid(cx, cy, indexing="ij")
    points = np.stack([grid_x.ravel(), grid_y.ravel()], axis=-1)   # (6, 2), ij order
    table = b.evaluate(points)
    np.testing.assert_array_equal(table, np.eye(6))


def test_evaluate_rejects_wrong_dimension():
    """A 2-axis basis rejects 1-D points (the d-mismatch guard)."""
    b = IndicatorBasis(edges_per_axis=(np.array([0.0, 1.0]), np.array([0.0, 1.0])))
    with pytest.raises(ValueError, match="shape"):
        b.evaluate(np.array([0.5, 0.5]))               # ambiguous: (N,) for d=2


# ── role separation (negative) ────────────────────────────────────────
def test_mesh_is_not_a_basis_but_yields_one():
    """The mesh YIELDS a basis view; it does NOT inherit Basis (role separation)."""
    mesh = Mesh1D(
        edges=_COARSE_EDGES, mat_ids=np.array([0, 1]),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    assert not isinstance(mesh, Basis)
    assert isinstance(mesh.indicator_basis(), Basis)
    # the yielded basis is built from the mesh's own edges
    np.testing.assert_array_equal(
        mesh.indicator_basis().edges_per_axis[0], _COARSE_EDGES,
    )
