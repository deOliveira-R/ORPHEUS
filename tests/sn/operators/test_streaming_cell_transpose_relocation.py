r"""#310 C1/C2 — the registered transpose-kernel pair: liveness, teeth, retirement.

C1 relocated the hand-coded DD reverse-walk transpose out of
``_OneDimScanWalk``'s ``visit`` closure into the registered scheme kernel
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.streaming_cell_transpose`
and made ``has_transpose_kernel`` DERIVE from the registration
(``DiscretizationSchemeBase.__init_subclass__`` — #310 ruling 2).  C2 gave
the reverse visit the SAME two-arm structure as the forward ``_apply_walk``:
the CARTESIAN arm rides the scheme-uniform batch VJP
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch_transpose`
(DD + LD, no scheme branch — the exact mirror of the forward's kernel arm),
the CURVILINEAR arm keeps the cell-balance kernel with the Morel–Montry
thread; the trait derives from BOTH registrations (the Cartesian batch VJP
always; the curvilinear kernel iff ``supports_curvilinear``).

A green frozen baseline does NOT prove the routing happened: if a lingering
hand-code twin still ran and a registered kernel were dead, the baseline
would stay green on the old code (vv-principles Mode 11).  This file
carries the liveness proof the spec (§2.2 of
``.claude/plans/archive/residue_verification_spec.md``) demands:

* the **two-arm Mode-11 sentinel** — in-process wraps on BOTH registered
  kernels: the batch VJP fires on the slab adjoint and stays silent on
  curvilinear; the cell-balance kernel fires on sphere/cylinder and stays
  silent on slab; the forward ``apply`` is the ``count == 0`` control for
  both (the wraps measure the executed arm, not import side-effects);
* the **angular-thread asymmetry** — ruling 1 keeps the Morel–Montry
  cotangent OUT of the kernels, so the closure's ``angular_adjoint`` wrap
  fires on sphere/cylinder and stays silent on slab;
* the **value teeth** (M-R1a-VALUE + the C2 batch sibling) — a coefficient
  perturbation inside each registered kernel moves the adjoint result O(1)
  (value-load-bearing, not merely called), with a bit-identical unpatched
  recompute as the determinism control;
* the **retirement pins** — the old hand-code spellings are GONE (exactly
  one ``denom * ob`` survivor: the degenerate-ordinate volumetric branch,
  a face-free relation outside the kernel contracts) and the affine-chain
  VJP pair has ONE spelling (#311).

All wraps are in-process ``monkeypatch`` (never ``git checkout`` — the
mutation-hygiene rule).  Foundation: software-structure invariants, no
equation label to link.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import orpheus.sn.loss_representation as loss_representation_pkg
import orpheus.sn.sweep.psi_half_angle_seed as psi_half_angle_seed_mod
import orpheus.transport.spatial.diamond as diamond_mod
from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase
from orpheus.sn.sweep.pole_angular_closure import MorelMontryAngularSweep

from tests.sn.operators.test_g_adjoint_reciprocity import (
    _BUILDERS,
    _loss_operator,
    _random_composite,
)

pytestmark = pytest.mark.foundation


def _counting_wrap(monkeypatch, owner: type, name: str) -> list[int]:
    r"""Wrap ``owner.<name>`` in-process; return the call-count cell."""
    calls: list[int] = []
    original = getattr(owner, name)

    def wrapper(self, *args, **kwargs):
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(owner, name, wrapper)
    return calls


@pytest.mark.parametrize(
    ("case", "expect_batch", "expect_cell_balance"),
    [
        ("slab_2g", True, False),
        ("sphere_2g", False, True),
        ("cyl_product_2g", False, True),
    ],
)
def test_registered_kernels_fire_on_their_arms(
    case, expect_batch, expect_cell_balance, monkeypatch,
):
    r"""Two-arm Mode-11 sentinel: each registered kernel is on ITS arm only.

    The reverse visit mirrors the forward ``_apply_walk``'s two arms (#310
    C2): the slab adjoint rides the scheme-uniform batch VJP
    (``residual_kernel_batch_transpose``) and never the cell-balance kernel;
    the curvilinear adjoint rides ``streaming_cell_transpose`` and never the
    batch VJP.  The forward ``apply`` is the ``count == 0`` control for
    BOTH.  A routing regression (an arm silently falling back to the other
    kernel, or a dead registered kernel behind a lingering hand-code twin)
    reds here even while the frozen baselines stay green (Mode 11).
    """
    sn, sig_t = _BUILDERS[case]()
    A = _loss_operator(sn, sig_t)
    rng = np.random.default_rng(310)
    phi = _random_composite(sn, rng)

    batch_calls = _counting_wrap(
        monkeypatch, DiamondDifference, "residual_kernel_batch_transpose",
    )
    cb_calls = _counting_wrap(
        monkeypatch, DiamondDifference, "streaming_cell_transpose",
    )

    A.apply(phi)
    assert len(batch_calls) == 0 and len(cb_calls) == 0, (
        f"[{case}] forward apply entered a transpose kernel "
        f"(batch={len(batch_calls)}, cell_balance={len(cb_calls)}) — the "
        "sentinels no longer discriminate the adjoint path"
    )

    A.apply_transpose(phi)
    if expect_batch:
        assert len(batch_calls) > 0, (
            f"[{case}] adjoint never entered residual_kernel_batch_transpose "
            "— the Cartesian batch VJP is DEAD on the executed path (Mode 11)"
        )
    else:
        assert len(batch_calls) == 0, (
            f"[{case}] curvilinear adjoint entered the CARTESIAN batch VJP "
            f"({len(batch_calls)} calls) — the arm asymmetry is broken"
        )
    if expect_cell_balance:
        assert len(cb_calls) > 0, (
            f"[{case}] adjoint never entered streaming_cell_transpose — the "
            "curvilinear cell-balance kernel is DEAD on the executed path "
            "(Mode 11)"
        )
    else:
        assert len(cb_calls) == 0, (
            f"[{case}] slab adjoint entered the CURVILINEAR cell-balance "
            f"kernel ({len(cb_calls)} calls) — the arm asymmetry is broken"
        )


@pytest.mark.parametrize(
    ("case", "expect_angular"),
    [("slab_2g", False), ("sphere_2g", True), ("cyl_product_2g", True)],
)
def test_angular_thread_fires_on_curvilinear_only(case, expect_angular, monkeypatch):
    r"""Ruling-1 liveness: the angular thread stays on the CLOSURE, and runs.

    The spatial-only kernel contract keeps the Morel–Montry cotangent with
    ``MorelMontryAngularSweep.angular_adjoint``.  Its wrap fires on the
    curvilinear rows and stays at zero on slab (whose closure is the
    identity class) — the config asymmetry that proves the angular wrap is
    real rather than vacuously green.
    """
    sn, sig_t = _BUILDERS[case]()
    A = _loss_operator(sn, sig_t)
    phi = _random_composite(sn, np.random.default_rng(311))

    calls = _counting_wrap(monkeypatch, MorelMontryAngularSweep, "angular_adjoint")
    A.apply_transpose(phi)

    if expect_angular:
        assert len(calls) > 0, (
            f"[{case}] curvilinear adjoint never entered "
            "MorelMontryAngularSweep.angular_adjoint — the angular thread "
            "is dead on the executed path"
        )
    else:
        assert len(calls) == 0, (
            f"[{case}] slab adjoint entered the Morel–Montry closure "
            f"({len(calls)} calls) — the asymmetry control is broken"
        )


def test_value_tooth_perturbed_kernel_reds():
    r"""M-R1a-VALUE: the relocated coefficient is value-load-bearing.

    Perturb the diamond-chain ``2.0`` to ``1.9`` inside an in-process
    replacement kernel: the adjoint result must move O(1) (a dead kernel or
    a value-inert call site would leave it unchanged).  The unpatched
    recompute is bit-identical to the baseline — the determinism control
    pinning the diff on the patch, not run-to-run noise.
    """
    sn, sig_t = _BUILDERS["sphere_2g"]()
    A = _loss_operator(sn, sig_t)
    phi = _random_composite(sn, np.random.default_rng(312))

    baseline = A.apply_transpose(phi)

    def perturbed(self, *, res_bar, psi_out_bar, denom, abs_mu_A_total, volume):
        psi_bar_cot = 1.9 * psi_out_bar          # the planted M-R1a-VALUE mutation
        psi_in_bar = -psi_out_bar
        psi_bar_cot += denom * res_bar / volume
        psi_in_bar += -(abs_mu_A_total)[None, :] * res_bar / volume
        return psi_bar_cot, psi_in_bar

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(DiamondDifference, "streaming_cell_transpose", perturbed)
        mutated = A.apply_transpose(phi)

    diff = float(
        np.max(np.abs(mutated.interior.values - baseline.interior.values))
    )
    scale = float(np.max(np.abs(baseline.interior.values)))
    assert diff > 1e-3 * scale, (
        f"perturbing the relocated 2.0 coefficient moved the adjoint by only "
        f"{diff:.3e} (scale {scale:.3e}) — the kernel is not value-load-bearing "
        "on the adjoint path"
    )

    recompute = A.apply_transpose(phi)
    assert np.array_equal(
        recompute.interior.values, baseline.interior.values
    ), "unpatched recompute is not bit-identical — the determinism control failed"


def test_value_tooth_perturbed_batch_kernel_reds():
    r"""The C2 batch sibling of M-R1a-VALUE: the batch VJP is value-load-bearing.

    Perturb the diagonal pullback ``denom·r†`` by 5% inside an in-process
    replacement of ``residual_kernel_batch_transpose``: the slab adjoint
    must move O(1) (a dead kernel or a value-inert call site would leave it
    unchanged).  The unpatched recompute is bit-identical to the baseline —
    the determinism control pinning the diff on the patch.
    """
    sn, sig_t = _BUILDERS["slab_2g"]()
    A = _loss_operator(sn, sig_t)
    phi = _random_composite(sn, np.random.default_rng(313))

    baseline = A.apply_transpose(phi)

    def perturbed(self, *, res_bar, psi_out_bar, s_axes, reaction_xs):
        denom, couplings = self._cartesian_streaming_diagonal(
            reaction_xs, s_axes,
        )
        psi_bar_cot = 0.95 * denom * res_bar     # the planted VALUE mutation
        psi_in_cots = []
        for c_a, out_bar_a in zip(couplings, psi_out_bar):
            avg_cot, in_cot = self.outgoing_face_from_average_transpose(
                out_bar_a, 0.5,
            )
            psi_bar_cot = psi_bar_cot + avg_cot
            psi_in_cots.append(in_cot - c_a * res_bar)
        return psi_bar_cot, tuple(psi_in_cots)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(DiamondDifference, "residual_kernel_batch_transpose", perturbed)
        mutated = A.apply_transpose(phi)

    diff = float(
        np.max(np.abs(mutated.interior.values - baseline.interior.values))
    )
    scale = float(np.max(np.abs(baseline.interior.values)))
    assert diff > 1e-3 * scale, (
        f"perturbing the batch VJP's denom pullback moved the slab adjoint by "
        f"only {diff:.3e} (scale {scale:.3e}) — the batch kernel is not "
        "value-load-bearing on the adjoint path"
    )

    recompute = A.apply_transpose(phi)
    assert np.array_equal(
        recompute.interior.values, baseline.interior.values
    ), "unpatched recompute is not bit-identical — the determinism control failed"


def test_retirement_no_handcoded_twin():
    r"""The old hand-code is GONE from the walk (no Mode-11 twin).

    ``2.0 * f_bar`` (the diamond-chain transpose spelling) must not appear
    in ``loss_representation``; ``denom * ob`` must appear EXACTLY once —
    the degenerate-ordinate volumetric branch (`pure-azimuthal set, zero
    face coupling`), which is a face-free relation deliberately outside the
    spatial kernel contract.  A second occurrence is a re-introduced twin.
    """
    src = Path(loss_representation_pkg.__file__).read_text()
    assert "2.0 * f_bar" not in src, (
        "the hand-coded diamond-chain transpose (2.0 * f_bar) is back in "
        "loss_representation — a Mode-11 twin of the registered kernel"
    )
    n_pullback = src.count("denom * ob")
    assert n_pullback == 1, (
        f"expected exactly one 'denom * ob' survivor (the degenerate-ordinate "
        f"volumetric branch), found {n_pullback} — a leg-walk twin of the "
        "registered kernel has been re-introduced (or the degenerate branch "
        "was refactored; update this pin deliberately)"
    )


def test_trait_derives_from_registration():
    r"""Ruling 2: ``has_transpose_kernel`` IS the registration (the law itself).

    The registration must COVER the scheme's forward span (#310 C2): the
    Cartesian batch VJP (``residual_kernel_batch_transpose``) always, plus
    the curvilinear cell-balance VJP (``streaming_cell_transpose``) iff the
    scheme claims ``supports_curvilinear``.  Declared-True-with-no-kernel
    (the pre-2.5a "predicate lie") stays unrepresentable.  The production
    values (DD ``True`` / LD ``False``) are pinned in
    ``test_ld_adjoint_deferral``; this test pins the derivation LAW on
    throwaway subclasses (no registry key — never registered).
    """

    def _batch_vjp(self, *, res_bar, psi_out_bar, s_axes, reaction_xs):
        return res_bar, psi_out_bar

    def _cell_balance_vjp(
        self, *, res_bar, psi_out_bar, denom, abs_mu_A_total, volume,
    ):
        return psi_out_bar, res_bar

    class _NoKernel(DiscretizationSchemeBase):  # type: ignore[abstract]
        pass

    class _BatchOnly(DiscretizationSchemeBase):  # type: ignore[abstract]
        residual_kernel_batch_transpose = _batch_vjp

    class _CellBalanceOnly(DiscretizationSchemeBase):  # type: ignore[abstract]
        streaming_cell_transpose = _cell_balance_vjp

    class _CurvClaimantBatchOnly(DiscretizationSchemeBase):  # type: ignore[abstract]
        supports_curvilinear = True
        residual_kernel_batch_transpose = _batch_vjp

    class _CurvFull(DiscretizationSchemeBase):  # type: ignore[abstract]
        supports_curvilinear = True
        residual_kernel_batch_transpose = _batch_vjp
        streaming_cell_transpose = _cell_balance_vjp

    assert _NoKernel.has_transpose_kernel is False, (
        "a subclass with no transpose-kernel override derived True"
    )
    assert _BatchOnly.has_transpose_kernel is True, (
        "a slab-only scheme registering the batch VJP derived False — the "
        "batch kernel alone covers a non-curvilinear forward span"
    )
    assert _CellBalanceOnly.has_transpose_kernel is False, (
        "a scheme registering ONLY the curvilinear cell-balance kernel "
        "derived True — the Cartesian batch VJP is always required"
    )
    assert _CurvClaimantBatchOnly.has_transpose_kernel is False, (
        "a curvilinear-claiming scheme without streaming_cell_transpose "
        "derived True — its curvilinear reverse arm would be ungated"
    )
    assert _CurvFull.has_transpose_kernel is True, (
        "a curvilinear scheme registering BOTH kernels derived False"
    )

    class _DeclaredLiar(DiscretizationSchemeBase):  # type: ignore[abstract]
        has_transpose_kernel = True  # the pre-2.5a predicate lie …

    assert _DeclaredLiar.has_transpose_kernel is False, (
        "a bare declaration overrode the derivation — the predicate lie is "
        "representable again"
    )


def test_affine_chain_transpose_single_source():
    r"""#311: the affine face-chain VJP pair has ONE spelling — the primitive.

    The three historically open-coded ``(2·f̄, −f̄)`` transpose sites (DD's
    ``streaming_cell_transpose`` + the two ``psi_half_angle_seed`` march
    reversals) all route through
    ``DiscretizationSchemeBase.outgoing_face_from_average_transpose``, and the
    hand-coded pair spellings are GONE.  A re-introduced open-coded pair is
    the single-source debt #311 retired: the LD kernel's ``w = 1/(1+k) ≠ ½``
    means any convention change fixed only in the primitive would silently
    miss a hardcoded ``w=½`` twin.

    #310 C4 added a THIRD diamond.py consumer — the row-march reverse's
    ``reflect_scan_coefficients_transpose`` extracts its ψ̄-independent
    ``(α, β_pullback)`` as the primitive's unit application — which is
    exactly the discipline this pin enforces (ride the primitive, never
    inline the pair), so the DD count moved 2 → 3.
    """
    dd_src = Path(diamond_mod.__file__).read_text()
    seed_src = Path(psi_half_angle_seed_mod.__file__).read_text()
    assert "2.0 * psi_out_bar" not in dd_src, (
        "the hand-transposed diamond-chain pair (2.0 * psi_out_bar) is back "
        "in diamond.py — an open-coded twin of "
        "outgoing_face_from_average_transpose"
    )
    assert "2.0 * f_bar" not in seed_src, (
        "a hand-transposed chain spelling (2.0 * f_bar) is back in "
        "psi_half_angle_seed.py — an open-coded twin of "
        "outgoing_face_from_average_transpose"
    )
    n_dd = dd_src.count("outgoing_face_from_average_transpose(")
    assert n_dd == 3, (
        f"expected exactly three primitive calls in diamond.py (the "
        f"cell-balance kernel's chain pair + the batch VJP's per-axis face "
        f"pair + the reflect-transpose coefficient extraction, #310 C4), "
        f"found {n_dd}"
    )
    n_seed = seed_src.count("outgoing_face_from_average_transpose(")
    assert n_seed == 3, (
        f"expected exactly three primitive calls in psi_half_angle_seed.py "
        f"(the Carlson reversal + the two residual-transpose legs), found "
        f"{n_seed}"
    )
