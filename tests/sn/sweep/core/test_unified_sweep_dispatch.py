"""Foundation tests for the sweep-strategy dispatch.

Round 2 of Wave D of the SN reshape campaign (Issue #161); migrated by
the **sweep-strategy carve** (C3.4/C3.5, plan ``sn_sweep_strategy.md``).

Historically ``transport_sweep`` chose the 1-D vs 2-D sweep body with a
scattered ``sn_mesh.reduced is not None`` branch, and these tests spied
on the chosen ``_sweep_*`` function being *called* (monkeypatching the
module-level name).  The carve replaced that branch with a first-class,
selectable :class:`~orpheus.sn.loss_representation.LossRepresentation`:
:func:`~orpheus.sn.loss_representation.default_for` picks the strategy, and
``transport_sweep`` delegates to ``strategy.sweep``.  The spy mechanism no
longer applies (the strategy holds its own reference to the sweep body),
so the dispatch contract is now pinned at its single source of truth:

* the SELECTION — ``default_for(mesh)`` returns the right strategy class:
  ALL 1-D meshes (slab, sphere, cylinder; ``is_1d``) →
  :class:`~orpheus.sn.loss_representation.CumprodScan`; multi-D Cartesian →
  :class:`~orpheus.sn.loss_representation.ScanMarch` (the S6.9 Fork-B2
  production default);
* the ROUTING — ``transport_sweep`` delegates to
  ``default_for(mesh).sweep(...)`` exactly once.

``-O``-safe (vv Mode 8): every gate is a ``np.testing`` / ``pytest.fail``
function call, NOT a bare ``assert`` (which the canonical ``python -O``
invocation would strip in any non-rewritten helper).

These are **software-contract** tests — the L1 transport math is
verified by the regression snapshots at
``tests/sn/regression/snapshots/`` (gate the bit-identical contract)
and the Wave C ``DiamondDifference`` hand-calc tests at
``tests/sn/spatial/test_diamond.py`` (gate the per-cell algebra).

Tagged ``@pytest.mark.foundation`` because:

* The dispatch is purely structural — no transport equation identity is
  being verified here.
* L1 transport accuracy is verified transitively via the existing
  MMS suite at ``tests/sn/l1_analytical/``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Mesh2D,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.spatial.linear_discontinuous import LinearDiscontinuous
from orpheus.sn.loss_representation import transport_sweep
from orpheus.sn.loss_representation import (
    CumprodScan,
    MovingFrontierWindow,
    ScanMarch,
    FullFieldWavefront,
    default_for,
)
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.fields.boundary_flux import BoundaryFlux


# ═══════════════════════════════════════════════════════════════════════
# Mesh fixtures
# ═══════════════════════════════════════════════════════════════════════

def _slab_sn_mesh(nx: int = 8, length: float = 1.0) -> SNMesh:
    """Slab SNMesh with vacuum BCs and Gauss-Legendre 1D quadrature."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials())


def _spherical_sn_mesh(nx: int = 8, radius: float = 1.0) -> SNMesh:
    """Spherical SNMesh with reflective inner / vacuum outer BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials())


def _cylindrical_sn_mesh(nx: int = 8, radius: float = 1.0) -> SNMesh:
    """Cylindrical SNMesh with reflective inner / vacuum outer BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _2d_sn_mesh(nx: int = 4, ny: int = 4) -> SNMesh:
    """2-D Cartesian SNMesh with vacuum BCs and Level-Symmetric S4."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _2d_ld_sn_mesh(nx: int = 4, ny: int = 3) -> SNMesh:
    """2-D Cartesian SNMesh carrying a Linear-Discontinuous scheme.

    NON-SQUARE by default (``nx=4, ny=3``) — an x↔y blindness defence (a
    square box can hide an axis swap).  Level-Symmetric S4 supplies genuine
    ``mu_y`` (#214-safe).  This is the mesh whose multi-D LD closure is
    BILINEAR (a per-axis slope moment) and therefore NOT scan-march
    compatible — the #240 D5-0 negative-routing fixture.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(),
                  scheme=LinearDiscontinuous())


# ═══════════════════════════════════════════════════════════════════════
# TestDispatchSelectsStrategy — the SELECTION half of the contract
# ═══════════════════════════════════════════════════════════════════════

class TestDispatchSelectsStrategy:
    """``default_for(mesh)`` selects the geometry-correct sweep strategy.

    The single source of truth that replaced the scattered
    ``reduced is not None`` branch: ALL 1-D meshes (slab, sphere, cylinder)
    → :class:`CumprodScan`; multi-D Cartesian → :class:`ScanMarch` (the
    S6.9 Fork-B2 production default, 2026-06-11 — measured 0.57–0.84× the
    window's sweep time at identical peak memory; the window and the
    full-field oracle stay explicit-select peers).  ``pytest.fail`` on
    mismatch fires under ``python -O``.
    """

    @pytest.mark.foundation
    def test_slab_selects_cumprod_scan(self):
        """Slab (1-D Cartesian) → CumprodScan."""
        sn_mesh = _slab_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("slab fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"slab → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_spherical_selects_cumprod_scan(self):
        """Spherical (1-D curvilinear) → CumprodScan (same as slab)."""
        sn_mesh = _spherical_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("spherical fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"sphere → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_cylindrical_selects_cumprod_scan(self):
        """Cylindrical (1-D curvilinear) → CumprodScan (same as slab/sphere)."""
        sn_mesh = _cylindrical_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("cylindrical fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"cylinder → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_2d_cartesian_selects_scan_march(self):
        """2-D Cartesian → ScanMarch (the S6.9 Fork-B2 production default).

        HISTORY: MovingFrontierWindow was the default through S6.9; the
        2026-06-11 measurement (sweep 0.57-0.84x, matvec 0.55-0.78x,
        identical peak memory) flipped the default.  The window remains a
        selectable peer (explicit construction), pinned by the window-forced
        end-to-end gates in ``test_scan_march_end_to_end.py``.
        """
        sn_mesh = _2d_sn_mesh()
        if sn_mesh.reduced is not None:
            pytest.fail("2-D fixture unexpectedly has reduced is not None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, ScanMarch):
            pytest.fail(
                f"2-D Cartesian → {type(strategy).__name__}, "
                f"expected ScanMarch"
            )


class TestD3SupportsMatrix:
    """C3.6 honest-supports pins at d=3 (test-architect G-c1..c4).

    The ``supports`` predicates read ``is_1d`` / ``is_cartesian`` /
    ``ndim`` (the dimensional narrowing) and — for the scan reps since
    #206 A1 — ``scheme.is_affine_scannable`` (the scheme must be
    scannable), so the PREDICATE pins stay duck-typed (the fake supplies
    an affine-scannable scheme; these pins document the dimensional
    narrowing at the predicate level). Since
    C5.5 (#225) a 3-axis ``SNMesh`` is constructible via
    ``from_axes``, so the SELECTION pin (G-c3) runs the LIVE
    ``default_for`` on a real mesh — it must route d=3 to the
    d-generic ``FullFieldWavefront`` spine, NOT misroute into
    ``ScanMarch``, whose row-march kernels unpack d=2 (its
    ``supports`` was narrowed to tell that truth in C3.6; widen it
    only WITH the scan(x)∘march(y,z) kernel generalization).
    """

    @staticmethod
    def _fake(ndim: int, *, cartesian: bool = True, facewise: bool = True):
        from types import SimpleNamespace
        return SimpleNamespace(
            is_1d=(ndim == 1), is_cartesian=cartesian, ndim=ndim,
            # #206 A1: the scan reps (CumprodScan / ScanMarch) now ALSO
            # gate on ``scheme.is_affine_scannable`` — the fake
            # supplies an affine-scannable scheme so these pins exercise
            # the DIMENSIONAL narrowing, not the scheme gate.
            #
            # #240 D5-0: the d≥2 ScanMarch arm gates on the DISTINCT
            # ``transverse_coupling_is_facewise`` (cross-axis separability),
            # NOT ``is_affine_scannable`` (single-axis prefix scannability).
            # A DD-shaped scheme reports facewise=True (admitted at d=2); an
            # LD-shaped scheme reports facewise=False (refused at d≥2 but
            # still affine-scannable in 1-D).  Defaults to the DD shape so the
            # dimensional-narrowing pins below exercise the right trait.
            scheme=SimpleNamespace(
                is_affine_scannable=True,
                transverse_coupling_is_facewise=facewise,
            ),
        )

    @pytest.mark.foundation
    def test_scan_march_refuses_d3(self):
        """G-c1: the narrowed predicate refuses a 3-axis Cartesian mesh."""
        if ScanMarch.supports(self._fake(3)).ok:
            pytest.fail(
                "ScanMarch.supports admitted a d=3 Cartesian mesh — its "
                "kernels unpack d=2; this would misroute production at d=3"
            )

    @pytest.mark.foundation
    def test_scan_march_keeps_d2_and_1d_coverage(self):
        """G-c2: the narrowing lost NO existing coverage — 2-D Cartesian,
        1-D Cartesian (slab), and 1-D curvilinear all still admitted."""
        for fake, name in (
            (self._fake(2), "2-D Cartesian"),
            (self._fake(1), "slab"),
            (self._fake(1, cartesian=False), "1-D curvilinear"),
        ):
            if not ScanMarch.supports(fake).ok:
                pytest.fail(f"ScanMarch.supports refused {name}")

    @pytest.mark.foundation
    def test_d3_cartesian_falls_through_to_full_field_spine(self):
        """G-c3 (FLIPPED to a real mesh in C5.5/#225): the LIVE
        ``default_for`` selection — not just the predicate walk — lands
        a real 3-axis Cartesian mesh on FullFieldWavefront, the
        never-stuck any-d spine."""
        import numpy as np
        from orpheus.derivations.common.xs_library import make_mixture
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.transport.mesh.axis import AxisMesh
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.sn.loss_representation import (
            FullFieldWavefront,
            default_for,
        )
        mix = make_mixture(
            sig_t=np.array([1.0]), sig_c=np.array([0.5]),
            sig_f=np.array([0.0]), nu=np.array([0.0]),
            chi=np.zeros(1), sig_s=np.array([[0.5]]),  # non-fissile ⇒ null χ (S10a guard)
        )
        mesh = SNMesh.from_axes(
            (
                AxisMesh(edges=np.linspace(0.0, 1.0, 3)),
                AxisMesh(edges=np.linspace(0.0, 1.0, 4)),
                AxisMesh(edges=np.linspace(0.0, 1.0, 5)),
            ),
            Quadrature.level_symmetric(sn_order=4),
            {0: mix},
        )
        selected = default_for(mesh)
        if type(selected) is not FullFieldWavefront:
            pytest.fail(
                f"d=3 Cartesian default_for → "
                f"{type(selected).__name__}, expected FullFieldWavefront "
                f"(the d-generic spine)"
            )

    @pytest.mark.foundation
    def test_d3_curvilinear_admitted_nowhere(self):
        """A hypothetical d=3 curvilinear mesh is refused by every
        registered representation (default_for would raise, not guess)."""
        fake = self._fake(3, cartesian=False)
        from orpheus.sn.loss_representation import LOSS_REPRESENTATIONS
        admitted = [
            cls for cls in LOSS_REPRESENTATIONS if cls.supports(fake).ok
        ]
        if admitted:
            pytest.fail(
                f"d=3 curvilinear admitted by "
                f"{[c.__name__ for c in admitted]} — expected none"
            )

    # ── #240 D5-0 — routing honesty (the confirmed-LIVE 2-D LD misroute) ──
    #
    # Pre-D5-0, ``ScanMarch.supports(2-D LD).ok == True`` so
    # ``default_for(2-D LD) == ScanMarch``; ScanMarch's row-march interior
    # runs INLINE Diamond-Difference with no scheme dispatch, so a 2-D LD
    # mesh SILENTLY computes DD (dropping LD's bilinear slope).  D5-0 narrows
    # the d≥2 ScanMarch arm to read ``transverse_coupling_is_facewise`` (the
    # cross-axis separability the row-march actually requires), which LD
    # leaves at ``False``.  The misroute becomes a refusal: the mesh falls
    # through to the wavefront, whose LD sweep RAISES the honest d=1-only
    # ``NotImplementedError`` (D5b closes the raise).  These gates pin the
    # contract BOTH ways (vv anti-pattern #11): LD must NOT route to
    # ScanMarch (negative), DD/1-D-LD must still route as before (positive).

    @pytest.mark.foundation
    def test_scan_march_refuses_2d_non_facewise_scheme_fake(self):
        """D5-0 (fake): the narrowed predicate refuses a 2-D mesh whose
        scheme reports ``transverse_coupling_is_facewise=False`` (LD-shaped),
        while still admitting the DD-shaped (facewise) 2-D mesh."""
        ld_like = self._fake(2, facewise=False)
        if ScanMarch.supports(ld_like).ok:
            pytest.fail(
                "ScanMarch.supports admitted a 2-D mesh whose scheme is NOT "
                "transverse-coupling-facewise (LD-shaped) — its row-march "
                "runs inline DD, silently dropping the bilinear slope"
            )
        dd_like = self._fake(2, facewise=True)
        if not ScanMarch.supports(dd_like).ok:
            pytest.fail(
                "ScanMarch.supports refused a 2-D facewise (DD-shaped) "
                "scheme — the narrowing lost the production default"
            )

    @pytest.mark.foundation
    def test_scan_march_refuses_2d_ld_real_mesh(self):
        """D5-0 negative (real mesh): a 2-D Cartesian LD ``SNMesh`` (bilinear
        slope coupling) must NOT select ScanMarch — its row-march kernel runs
        inline DD, silently dropping LD's slope.  This is the CONFIRMED-LIVE
        misroute (``default_for(2-D LD) == ScanMarch`` pre-D5-0)."""
        sn = _2d_ld_sn_mesh()
        if ScanMarch.supports(sn).ok:
            pytest.fail(
                "ScanMarch.supports admitted a 2-D LD mesh — its row-march "
                "kernel runs inline DD (loss_representation interior), "
                "silently dropping LD's bilinear slope; the misroute computes "
                "DD, not LD"
            )

    @pytest.mark.foundation
    def test_2d_ld_default_for_routes_to_wavefront(self):
        """D5-0 selection: ``default_for`` on a 2-D LD mesh lands on a
        WAVEFRONT (MovingFrontierWindow / FullFieldWavefront), NEVER
        ScanMarch.  (D5b-S2 lands the multi-D LD ``cell_kernel_batch`` so the
        wavefront now RUNS the bilinear UBLD; what is forbidden remains the
        silent ScanMarch inline-DD path.)"""
        sn = _2d_ld_sn_mesh()
        selected = default_for(sn)
        if isinstance(selected, ScanMarch):
            pytest.fail(
                "2-D LD routed to ScanMarch (inline DD) — the misroute"
            )
        if not isinstance(selected, (MovingFrontierWindow, FullFieldWavefront)):
            pytest.fail(
                f"2-D LD default_for → {type(selected).__name__}, expected a "
                f"wavefront (MovingFrontierWindow / FullFieldWavefront)"
            )

    @pytest.mark.foundation
    def test_2d_ld_sweep_runs_genuine_ld_not_dd(self):
        """D5b-S2 closes the D5-0 honest-interim raise: a 2-D LD fixed-source
        now RUNS (the multi-D bilinear UBLD kernel) and produces a GENUINELY
        DIFFERENT converged flux from DD on the same mesh/materials/source —
        the catcher that LD is computing LD, not silently collapsing to DD (the
        D5-0 misroute would have given DD≡LD).  Inverts the former
        ``test_2d_ld_sweep_raises_not_silently_dd`` (#240 D5b)."""
        import numpy as np
        from orpheus.derivations.common.xs_library import make_mixture
        from orpheus.sn import solve_sn_fixed_source
        from orpheus.sn.spatial import DiamondDifference

        mix = make_mixture(
            sig_t=np.array([1.0]), sig_c=np.array([0.5]),
            sig_f=np.array([0.0]), nu=np.array([0.0]),
            chi=np.zeros(1), sig_s=np.array([[0.5]]),  # non-fissile ⇒ null χ (S10a guard)
        )
        nx, ny = 4, 3                          # coarse — where O(h²) closures diverge
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 1.0, nx + 1),
            edges_y=np.linspace(0.0, 1.0, ny + 1),
            mat_map=np.zeros((nx, ny), dtype=int),
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        N = quad.weights.size
        Q = np.ones((N, 1, nx, ny))            # (N, ng, nx, ny) per-ordinate
        res_ld = solve_sn_fixed_source(
            {0: mix}, mesh, quad, Q,
            scheme=LinearDiscontinuous(), max_inner=200, inner_tol=1e-10,
        )
        res_dd = solve_sn_fixed_source(
            {0: mix}, mesh, quad, Q,
            scheme=DiamondDifference(), max_inner=200, inner_tol=1e-10,
        )
        phi_ld = res_ld.scalar_flux.values
        phi_dd = res_dd.scalar_flux.values
        if not np.all(np.isfinite(phi_ld)):
            pytest.fail("2-D LD fixed-source produced non-finite flux")
        if np.allclose(phi_dd, phi_ld, rtol=1e-3):
            pytest.fail(
                "2-D LD converged flux ≈ DD — LD is silently computing DD "
                "(the D5-0 misroute regressed, or the bilinear closure was "
                "not exercised)"
            )


# ═══════════════════════════════════════════════════════════════════════
# TestSchemeTraitProbe — strategy-free scheme-property probe (D5-0)
# ═══════════════════════════════════════════════════════════════════════


class TestSchemeTraitProbe:
    """The scheme property is answerable with NO strategy in scope.

    The cross-domain-attacker's second-consumer discriminator (#240 D5-0
    Frame 1): if the only way to learn whether a scheme's multi-D transverse
    coupling is facewise/separable is ``ScanMarch.supports(mesh).ok``, the
    property is strategy-entangled (a frame leak).  These probes read the
    trait DIRECTLY off the scheme class — no ``ScanMarch``, no ``supports``,
    no mesh — proving ``transverse_coupling_is_facewise`` is a genuine scheme
    property the diffusion ADI / line-SOR preconditioner (#240's next
    consumer) can reuse with no rename.
    """

    @pytest.mark.foundation
    def test_dd_reports_facewise_transverse_coupling(self):
        """DD's slopeless cell-average closure couples transverse axes by a
        0th-order face value → facewise/separable → True."""
        if DiamondDifference.transverse_coupling_is_facewise is not True:
            pytest.fail(
                "DiamondDifference.transverse_coupling_is_facewise is "
                f"{DiamondDifference.transverse_coupling_is_facewise!r}, "
                "expected True (DD's transverse term is a 0th-order face value)"
            )

    @pytest.mark.foundation
    def test_ld_reports_non_facewise_transverse_coupling(self):
        """LD's bilinear DG-P1 closure couples transverse axes by a 1st-order
        slope moment → non-separable → False (the default it opts out to)."""
        if LinearDiscontinuous.transverse_coupling_is_facewise is not False:
            pytest.fail(
                "LinearDiscontinuous.transverse_coupling_is_facewise is "
                f"{LinearDiscontinuous.transverse_coupling_is_facewise!r}, "
                "expected False (LD's transverse coupling is a 1st-order "
                "slope moment, not a face value)"
            )

    @pytest.mark.foundation
    def test_facewise_distinct_from_affine_scannable(self):
        """The two traits answer DIFFERENT questions: LD is affine-scannable
        (single-axis, True) yet NOT transverse-coupling-facewise (cross-axis,
        False).  If they coincided on LD the conflation that drove the D5-0
        misroute would still be latent."""
        if LinearDiscontinuous.is_affine_scannable is not True:
            pytest.fail("LD.is_affine_scannable expected True (1-D scannable)")
        if LinearDiscontinuous.transverse_coupling_is_facewise is not False:
            pytest.fail(
                "LD.transverse_coupling_is_facewise expected False — the two "
                "traits MUST diverge on LD (single-axis vs cross-axis claim)"
            )
        # DD coincides on both (it satisfies both claims) — the trait split is
        # only OBSERVABLE on a scheme where they diverge, which is LD.
        if not (DiamondDifference.is_affine_scannable
                and DiamondDifference.transverse_coupling_is_facewise):
            pytest.fail("DD expected True on both traits")


# ═══════════════════════════════════════════════════════════════════════
# TestTransportSweepDelegatesToStrategy — the ROUTING half of the contract
# ═══════════════════════════════════════════════════════════════════════

class TestTransportSweepDelegatesToStrategy:
    """``transport_sweep`` routes through ``default_for(mesh).sweep`` once.

    The ROUTING half of the dispatch contract (the SELECTION half is
    :class:`TestDispatchSelectsStrategy`).  Since S6.4(f) ``transport_sweep``
    lives IN ``loss_representation`` and reads the module-global
    ``default_for`` at call time, so patching
    ``loss_representation.default_for`` is seen — the spy confirms the
    dispatcher *delegates to the selected strategy* rather than re-deciding
    the branch itself.  Geometry-agnostic: the same delegation holds for the
    1-D scan and the 2-D window.  (S6.5 scope note: this is the operator-FREE
    functional entry, which legitimately selects per call; the OPERATOR's
    solve consumes its own ``loss_representation`` instance and is pinned by
    ``test_one_representation_instance.py``.)
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "mesh_factory", [_slab_sn_mesh, _2d_sn_mesh], ids=["slab", "cart2d"],
    )
    def test_delegates_to_selected_strategy(self, monkeypatch, mesh_factory):
        """transport_sweep calls the selected strategy's ``sweep`` exactly once."""
        sn_mesh = mesh_factory()
        import orpheus.sn.loss_representation as loss_representation

        selected = type(default_for(sn_mesh)).__name__
        calls = {"sweep": 0}
        N, ng = sn_mesh.quad.N, sn_mesh.ng
        spatial = sn_mesh.spatial_shape

        class _SpyStrategy:
            def sweep(self, *args, **kwargs):
                calls["sweep"] += 1
                return (np.zeros((N, ng, *spatial)), np.zeros((ng, *spatial)))

        monkeypatch.setattr(
            loss_representation, "default_for", lambda mesh: _SpyStrategy(),
        )

        # Σ_t is (ng, *spatial) at any rank — (ng, nx) for 1-D,
        # (ng, nx, ny) for 2-D Cartesian (C5.2: no phantom ny).
        sig_t = np.ones((ng, *spatial))
        source = AngularSourceSink.zeros_on(sn_mesh)
        transport_sweep(source, sig_t, sn_mesh, BoundaryFlux.zeros_on(sn_mesh))

        if calls["sweep"] != 1:
            pytest.fail(
                f"transport_sweep delegated to strategy.sweep "
                f"{calls['sweep']} times (selected {selected}), expected "
                f"exactly 1"
            )


# ═══════════════════════════════════════════════════════════════════════
# TestDefaultDiscretizationScheme
# ═══════════════════════════════════════════════════════════════════════

class TestDefaultDiscretizationScheme:
    """``SNMesh.scheme`` defaults to :class:`DiamondDifference`.

    The default is what guarantees bit-identity with the pre-Wave-D
    sweep — DD's per-cell math is a bit-identical extraction of the
    inlined sweep's algebra.  Wave C-extension will let users pass
    LD / EC / Step strategies via the constructor argument.
    """

    @pytest.mark.foundation
    def test_default_is_diamond_difference(self):
        """No ``scheme`` argument → defaults to DD."""
        sn_mesh = _slab_sn_mesh()
        if not isinstance(sn_mesh.scheme, DiamondDifference):
            pytest.fail(
                f"default scheme is {type(sn_mesh.scheme).__name__}, "
                f"expected DiamondDifference"
            )

    @pytest.mark.foundation
    def test_explicit_scheme_honored(self):
        """User-passed ``scheme`` is stored on the mesh."""
        custom = DiamondDifference()  # the only strategy that ships in Wave D
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 9),
            mat_ids=np.zeros(8, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(), scheme=custom)
        if sn_mesh.scheme is not custom:
            pytest.fail("explicit scheme was not stored on the mesh")


@pytest.mark.foundation
class TestOneDimScanWalkFrame:
    """#206 Phase B: the 1-D sweep is OWNED by ``_OneDimScanWalk`` — the 1-D
    analogue of ``_OctantWalk`` — shared by ``CumprodScan.sweep`` and the
    ``ScanMarch`` 1-D branch (the former module-level helpers
    ``_sweep_1d_unified`` / ``_run_1d_sweep`` / ``_ensure_*_cache`` were
    relocated into it, bit-identical).
    """

    def test_frame_resolves_as_frozen_mesh_holder(self):
        """[foundation] The frame resolves from its new home and is a frozen
        ``mesh``-only dataclass (mirrors ``_OctantWalk``)."""
        import dataclasses

        from orpheus.sn.loss_representation import _OneDimScanWalk

        if not dataclasses.is_dataclass(_OneDimScanWalk):
            pytest.fail("_OneDimScanWalk must be a dataclass (mirror _OctantWalk)")
        names = [f.name for f in dataclasses.fields(_OneDimScanWalk)]
        if names != ["mesh"]:
            pytest.fail(
                f"_OneDimScanWalk fields must be exactly ['mesh'] (a frozen frame "
                f"holding only the mesh, like _OctantWalk); got {names}"
            )
        if not _OneDimScanWalk.__dataclass_params__.frozen:  # type: ignore[attr-defined]
            pytest.fail("_OneDimScanWalk must be frozen (immutable frame)")

    def test_frame_is_kernel_parameterized_not_boolean(self):
        """[foundation] ``_OneDimScanWalk`` forks at the geometry branch (and,
        in Phase C, a kernel object), NEVER a boolean ``is_solve``/``is_apply``
        flag (coding-elegance Smell #3).  This is the twin-path tripwire that
        guards ``_OctantWalk`` (``test_one_octant_walk``), applied to the 1-D
        frame so the Phase-C matvec attachment cannot degrade into a flag.
        AST identifiers only (docstrings naming the anti-pattern don't trip it);
        ``-O``-safe.
        """
        import ast
        import inspect
        import textwrap

        from orpheus.sn.loss_representation import _OneDimScanWalk

        smells = {"is_solve", "is_apply", "is_matvec"}
        tree = ast.parse(textwrap.dedent(inspect.getsource(_OneDimScanWalk)))
        identifiers = {
            node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
        } | {
            node.attr for node in ast.walk(tree) if isinstance(node, ast.Attribute)
        } | {
            node.arg for node in ast.walk(tree)
            if isinstance(node, (ast.arg, ast.keyword)) and node.arg is not None
        }
        offenders = sorted(identifiers & smells)
        if offenders:
            pytest.fail(
                f"_OneDimScanWalk carries boolean direction flag(s) {offenders} — "
                "the 1-D carve degraded into the boolean-flag anti-pattern "
                "(coding-elegance Smell #3). The cell operation MUST attach as a "
                "kernel OBJECT/callable (Phase C), not a flag."
            )
