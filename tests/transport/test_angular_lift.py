r"""**G5.8** — ``AngularLift``'s own laws (the intrinsic-properties standard; CS4c step 5).

The lift of an ENERGY binding onto the angular composite is a mathematical
object with defining laws, gated here on its own terms — not through the
solvers that consume it:

1. **Linearity**, with an ACTIVATION leg (`lessons L40c`: a zero morphism
   satisfies every linearity row with both sides structurally zero).
2. **The ℓ = 0 conjugation identity** ``lift(E)(ψ) = R₀ E M₀ ψ / W``: the
   reaction-rate fast path the base runs (``∫ψ dΩ``, then ``E``, then the
   producer-side ``/W``) against the frame form the transpose is spelled
   with (``full_*_kernel.apply(ψ.values) / W``) — two DIFFERENT reduction
   trees over the same factors. `[M]` 2026-09-04 (the verification plan,
   200-seed sweep, GL8 slab / mixture A 2g / 20 cells): the pure ℓ = 0
   lifts — ``F``, and ``S`` at ``L = 0`` — agree **bit-for-bit, 200/200**;
   ``S`` at ``L = 1`` does NOT (0/200, max |Δ| 2.2e-16), because the ℓ ≥ 1
   half is summed in a different order on the two routes. So the identity
   is pinned ``array_equal`` on the BASE's law and at the draw-stable
   absolute band ``max |Δ| ≤ 2.3e-16`` on the subclass sum (never a nulp
   band — `plan-authoring` 2026-08-28: a nulp band on near-zero outputs
   pins a seed). The partition falls exactly along the base/subclass line
   — the measurement that R-1's split is the right one.
3. **The transpose** is the conjugated product's reversal ``/W``:
   Euclidean reciprocity ``⟨T ψ, χ⟩ = ⟨ψ, Tᵀ χ⟩`` on raw arrays.
4. **The datum → energy derivation**: the energy binding is DERIVED from
   the datum on the CODOMAIN's scalar sub-space (F-1), and the role's own
   ``isotropic_binding`` is what the transfer lift derives.
5. The base is ABSTRACT, and its selection refuses a third interior
   (the admission legs G5.3d exercises on the moment sibling live in
   ``tests/sn/operators/test_moment_domain_binding.py``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.frames.harmonic_frame import HarmonicFrame
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.angular_lift import AngularLift
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicFission,
    IsotropicScattering,
)
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.source_sinks import AngularSourceSink

pytestmark = pytest.mark.foundation

#: The draw-stable band for the ℓ ≥ 1 sum (`[M]` 200 seeds: max |Δ| 2.2e-16).
_L1_ABS_BAND = 2.3e-16


def _sn() -> SNMesh:
    carrier = MaterialMesh.from_materials({0: get_mixture("A", "2g")})
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 21), mat_ids=np.zeros(20, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=8), carrier.materials)


def _mat_xs(sn: SNMesh):
    return sn.material_xs_field()


def _state(sn: SNMesh, seed: int) -> FullField:
    rng = np.random.default_rng(seed)
    interior = sn.full_field_space.interior_space
    assert interior is not None
    return FullField(
        interior=AngularFlux(values=rng.standard_normal(interior.shape), space=interior),
        boundary=AngularBoundaryFlux.zeros(sn.angular_trace),
    )


def _cotangent(sn: SNMesh, seed: int) -> FullField:
    rng = np.random.default_rng(seed)
    interior = sn.full_field_space.interior_space
    assert interior is not None
    return FullField(
        interior=AngularSourceSink(values=rng.standard_normal(interior.shape), space=interior),
        boundary=AngularBoundaryFlux.zeros(sn.angular_trace).into_role(
            __import__("orpheus.transport.fields._bases", fromlist=["FieldRole"]).FieldRole.SOURCE_SINK,
            np.zeros(sn.angular_trace.shape),
        ),
    )


def _lifts(sn: SNMesh) -> dict[str, AngularLift]:
    mat_xs = _mat_xs(sn)
    space = sn.full_field_space
    return {
        "F": FissionOperator.from_solver_data(mat_xs=mat_xs, space=space),
        "S_L0": ScatteringOperator.from_solver_data(mat_xs=mat_xs, scattering_order=0, space=space),
        "S_L1": ScatteringOperator.from_solver_data(mat_xs=mat_xs, scattering_order=1, space=space),
        "N2N_L1": N2NOperator.from_solver_data(mat_xs=mat_xs, scattering_order=1, space=space),
    }


def _frame_form(lift: AngularLift):
    return lift.full_fission_kernel if isinstance(lift, FissionOperator) else lift.full_transfer_kernel  # type: ignore[attr-defined]


@pytest.mark.parametrize("key", ["F", "S_L0", "S_L1", "N2N_L1"])
def test_activation_and_linearity(key):
    sn = _sn()
    lift = _lifts(sn)[key]
    a, b = _state(sn, 1), _state(sn, 2)
    out_a = lift.apply(a).interior.values
    if key == "N2N_L1":
        # mixture A has no (n,2n) data: the lift is honestly the zero
        # morphism here, so the activation leg is the fact recorded, and
        # the linearity row below is structurally vacuous for this key.
        assert not out_a.any()
        return
    assert np.abs(out_a).max() > 0.0, "activation: the lift must move a random flux"
    lhs = lift.apply(2.0 * a - 3.0 * b).interior.values
    rhs = 2.0 * out_a - 3.0 * lift.apply(b).interior.values
    np.testing.assert_allclose(lhs, rhs, rtol=1e-13, atol=1e-15)


@pytest.mark.parametrize("key", ["F", "S_L0"])
def test_the_l0_conjugation_identity_is_bit_exact_on_the_base(key):
    r"""``lift(ψ) == R₀ E M₀ ψ / W`` — the fast path and the frame form agree
    bit-for-bit on the PURE ℓ = 0 lifts (`[M]` 200/200)."""
    sn = _sn()
    lift = _lifts(sn)[key]
    for seed in range(8):
        psi = _state(sn, seed)
        fast = lift.apply(psi).interior.values
        frame = np.asarray(_frame_form(lift).apply(psi.interior.values)) / lift.total_weight
        np.testing.assert_array_equal(fast, frame)


def test_the_l1_sum_agrees_at_the_draw_stable_absolute_band():
    r"""On the ANISOTROPIC subclass the two routes sum ℓ ≥ 1 in different
    orders: not bit-exact (`[M]` 0/200), pinned at the draw-stable absolute
    band. Both facts are the claim — a bit-exact green here would mean the
    ℓ ≥ 1 body is not running."""
    sn = _sn()
    lift = _lifts(sn)["S_L1"]
    assert not lift.is_isotropic
    worst = 0.0
    any_differs = False
    for seed in range(8):
        psi = _state(sn, seed)
        fast = lift.apply(psi).interior.values
        frame = np.asarray(_frame_form(lift).apply(psi.interior.values)) / lift.total_weight
        delta = np.abs(fast - frame).max()
        worst = max(worst, float(delta))
        any_differs |= not np.array_equal(fast, frame)
    assert worst <= _L1_ABS_BAND, f"max |Δ| = {worst:.3e} exceeds the band"
    assert any_differs, "the ℓ ≥ 1 sum came out bit-exact — is the anisotropic body running?"


@pytest.mark.parametrize("key", ["F", "S_L0", "S_L1"])
def test_transpose_reciprocity_on_raw_arrays(key):
    sn = _sn()
    lift = _lifts(sn)[key]
    psi, chi = _state(sn, 5), _cotangent(sn, 6)
    lhs = float(np.vdot(lift.apply(psi).interior.values, chi.interior.values))
    rhs = float(np.vdot(psi.interior.values, lift.apply_transpose(chi).interior.values))
    assert lhs != 0.0
    assert abs(lhs - rhs) <= 1e-13 * abs(lhs), (lhs, rhs)
    assert not lift.apply_transpose(chi).boundary.values.any()


def test_the_energy_binding_is_derived_from_the_datum_on_the_codomain_scalar_subspace():
    sn = _sn()
    lifts = _lifts(sn)
    interior = sn.full_field_space.interior_space
    assert interior is not None and interior.axes is not None
    scalar = FunctionSpace.of_axes(*interior.axes[1:])
    F = lifts["F"]
    assert isinstance(F.isotropic_energy, IsotropicFission)
    assert F.isotropic_energy.fission is F.fission
    assert F.isotropic_energy.domain == scalar and F.isotropic_energy.codomain == scalar
    assert F.isotropic_energy is F.isotropic_energy  # cached once
    S = lifts["S_L1"]
    assert type(S).isotropic_binding is IsotropicScattering
    assert type(S.isotropic_energy) is IsotropicScattering
    assert S.isotropic_energy.transfer.order == 0  # the P0 head, nothing richer
    assert S.isotropic_energy.domain == scalar
    assert F.frame is S.frame.__class__.for_space(interior, 0)  # the interned hub frame
    assert F.total_weight == S.total_weight == 2.0


def test_the_base_is_abstract_and_refuses_a_third_interior():
    sn = _sn()
    interior = sn.full_field_space.interior_space
    assert interior is not None
    frame = HarmonicFrame.for_space(interior, 0)
    space = sn.full_field_space
    with pytest.raises(TypeError, match="abstract"):
        AngularLift(  # type: ignore[abstract]
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space, codomain=space,
        )
    # a THIRD interior: a composite whose bulk is neither face end
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.transport.material_field import FissionMaterialField

    other = FullFieldSpace.from_blocks(
        FunctionSpace(name="third", shape=tuple(interior.shape)), sn.angular_trace,
    )
    with pytest.raises(TypeError, match="neither end of the analysis face"):
        FissionOperator(
            FissionMaterialField.from_material_xs(_mat_xs(sn)),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=other, codomain=space,
        )


def test_the_lift_population_is_the_two_cores_and_their_roles():
    """G5.4b's companion at the base: the lift's subclass population."""
    def walk(c):
        direct = set(c.__subclasses__())
        return direct.union(*(walk(s) for s in direct))

    names = {c.__name__ for c in walk(AngularLift)}
    assert names == {"TransferOperator", "FissionOperator", "ScatteringOperator", "N2NOperator"}, names
