r"""#310 C1 — the DD kernel-pair relocation: liveness, teeth, and retirement.

The C1 carve relocated the hand-coded DD reverse-walk transpose out of
``_OneDimScanWalk``'s ``visit`` closure into the registered scheme kernel
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.streaming_cell_transpose`,
and made ``has_transpose_kernel`` DERIVE from that registration
(``DiscretizationSchemeBase.__init_subclass__`` — #310 ruling 2).

A green frozen baseline does NOT prove the relocation happened: if a
lingering hand-code twin still ran and the registered kernel were dead, the
baseline would stay green on the old code (vv-principles Mode 11).  This
file carries the liveness proof the spec (§2.2 of
``.claude/plans/residue_verification_spec.md``) demands:

* the **Mode-11 sentinel** — an in-process wrap on the registered kernel,
  ``count > 0`` on the adjoint path with a forward-run ``count == 0``
  control (the wrap measures the path, not import side-effects);
* the **angular-thread asymmetry** — ruling 1 keeps the Morel–Montry
  cotangent OUT of the kernel, so the closure's ``angular_adjoint`` wrap
  fires on sphere/cylinder and stays silent on slab (the config asymmetry
  proving the separation is real on the executed path);
* the **value tooth** (M-R1a-VALUE) — a ±5% perturbation of the relocated
  ``2.0`` coefficient moves the adjoint result O(1) (the kernel is
  value-load-bearing, not merely called), with a bit-identical unpatched
  recompute as the determinism control;
* the **retirement pin** — the old hand-code spellings are GONE from the
  walk (exactly one ``denom * ob`` survivor: the degenerate-ordinate
  volumetric branch, which is a face-free relation and NOT the kernel's).

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

_ADJOINT_CASES = ("slab_2g", "sphere_2g", "cyl_product_2g")


def _counting_wrap(monkeypatch, owner: type, name: str) -> list[int]:
    r"""Wrap ``owner.<name>`` in-process; return the call-count cell."""
    calls: list[int] = []
    original = getattr(owner, name)

    def wrapper(self, *args, **kwargs):
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(owner, name, wrapper)
    return calls


@pytest.mark.parametrize("case", _ADJOINT_CASES)
def test_relocated_kernel_fires_on_adjoint_path(case, monkeypatch):
    r"""Mode-11 sentinel: the registered DD kernel is ON the adjoint path.

    ``count > 0`` during ``(L+C−B).apply_transpose`` on every 1-D geometry,
    and — the discriminating control — ``count == 0`` during the FORWARD
    ``apply`` on the same operator: the sentinel measures the reverse walk
    specifically.  A relocation that left the kernel dead (hand-code twin
    still running) reds here even while the frozen baselines stay green.
    """
    sn, sig_t = _BUILDERS[case]()
    A = _loss_operator(sn, sig_t)
    rng = np.random.default_rng(310)
    phi = _random_composite(sn, rng)

    calls = _counting_wrap(monkeypatch, DiamondDifference, "streaming_cell_transpose")

    A.apply(phi)
    assert len(calls) == 0, (
        f"[{case}] forward apply entered streaming_cell_transpose "
        f"({len(calls)} calls) — the sentinel no longer discriminates the "
        "adjoint path"
    )

    A.apply_transpose(phi)
    assert len(calls) > 0, (
        f"[{case}] adjoint apply_transpose never entered the registered "
        "streaming_cell_transpose — the relocated kernel is DEAD on the "
        "executed path (Mode 11)"
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

    A subclass overriding ``streaming_cell_transpose`` derives ``True``; one
    that does not derives ``False`` — declared-True-with-no-kernel (the
    pre-2.5a "predicate lie") is unrepresentable.  The production values
    (DD ``True`` / LD ``False``) are pinned in ``test_ld_adjoint_deferral``;
    this test pins the derivation LAW on throwaway subclasses (no registry
    key — never registered).
    """

    class _NoKernel(DiscretizationSchemeBase):  # type: ignore[abstract]
        pass

    class _WithKernel(DiscretizationSchemeBase):  # type: ignore[abstract]
        def streaming_cell_transpose(
            self, *, res_bar, psi_out_bar, denom, abs_mu_A_total, volume,
        ):
            return psi_out_bar, res_bar

    assert _NoKernel.has_transpose_kernel is False, (
        "a subclass with no streaming_cell_transpose override derived True"
    )
    assert _WithKernel.has_transpose_kernel is True, (
        "a subclass registering streaming_cell_transpose derived False"
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
    assert n_dd == 1, (
        f"expected exactly one primitive call in diamond.py (the kernel's "
        f"chain pair), found {n_dd}"
    )
    n_seed = seed_src.count("outgoing_face_from_average_transpose(")
    assert n_seed == 3, (
        f"expected exactly three primitive calls in psi_half_angle_seed.py "
        f"(the Carlson reversal + the two residual-transpose legs), found "
        f"{n_seed}"
    )
