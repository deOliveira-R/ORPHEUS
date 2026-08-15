r"""``face_transmission_spectrum`` — is a closure's face mode UNDAMPED?

The closure half of the gauge-freedom predicate (#344).  The geometry
half is ``SNMesh.reflective_axis_pairs``; together they decide whether
the assembled loss operator :math:`A = L + C - S - B` has a null space,
and therefore whether the returned boundary trace needs gauge-fixing.

**What these gates protect.** The verdict is ASKED of each closure's own
:meth:`cell_kernel_batch` — nothing is tabulated, so a new scheme
answers for itself.  The risk that buys is that the *asking* silently
stops working: a driver bug, a shape change, or a scheme that grows a
moment tail could turn a real answer into a wrong one or into a
spurious ``UNDETERMINED``.  Both happened while this was being written
(a face reduced to moment 0 only; a wrong ``Q_cells`` shape), and both
would have read as "LD is unclassifiable" rather than as author error.
So the DAMPED rows matter as much as the UNDAMPED ones.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.spatial.scheme import (
    DiscretizationSchemeBase,
    FaceModeDamping,
)

_DIAMOND = DiscretizationSchemeBase.registry["diamond_difference"]
_LINEAR_DISCONTINUOUS = DiscretizationSchemeBase.registry[
    "linear_discontinuous"
]


@pytest.mark.foundation
@pytest.mark.parametrize("ndim", [2, 3])
def test_diamond_leaves_a_face_mode_UNDAMPED(ndim):
    r"""Diamond's sawtooth transmits at exactly :math:`|\lambda| = 1`.

    :math:`\Sigma = (2/D)\mathbf 1 w^{\mathsf T} - I` carries eigenvalue
    :math:`-1` with multiplicity :math:`d-1` on
    :math:`\{v : w^{\mathsf T}v = 0\}` — the cell-average-blind mode,
    which drives :math:`\psi_c = 0` so the absorption term never sees
    it.  `[M]` reads EXACTLY ``1.0`` (the eigenvalue is ``-1`` by
    construction, not by cancellation), which is why the tolerance only
    has to absorb eigensolver noise.
    """
    result = _DIAMOND().face_transmission_spectrum(ndim)
    assert result.damping is FaceModeDamping.UNDAMPED
    assert result.spectral_radius == pytest.approx(1.0, abs=1e-12)
    assert result.undetermined_because is None


@pytest.mark.foundation
def test_one_dimension_has_no_face_mode_to_leave_undamped():
    r"""``ndim == 1`` is DAMPED for BOTH closures — structurally, not by luck.

    A 1-D cell has a single face per direction, so there is no
    cell-average-blind combination to sustain; the only eigenvalue is
    the absorption-damped :math:`1 - 2\Sigma_t V/D`.  This is the row
    that shows the predicate needs no special-casing for ``d = 1``: it
    falls out of the spectral radius.
    """
    for scheme_type in (_DIAMOND, _LINEAR_DISCONTINUOUS):
        result = scheme_type().face_transmission_spectrum(1)
        assert result.damping is FaceModeDamping.DAMPED, scheme_type.__name__
        assert result.spectral_radius is not None
        assert result.spectral_radius < 1.0


@pytest.mark.foundation
def test_linear_discontinuous_damps_its_face_modes_at_d2():
    """LD is the NEGATIVE control — and it must be a real reading, not a dodge.

    `[M]` ``rho = 0.8607021518`` on the thin probe.  LD's face carries a
    moment tail (``spatial_basis_per_axis ** (ndim-1)`` entries), so the
    transmission is ``4x4`` here rather than ``2x2`` — a driver that
    silently reduced the face to its first moment would still produce a
    plausible sub-unit number, so the SHAPE is asserted too.

    This row independently reproduces #344's measured ``dim ker A = 0``
    for LD on the identical all-reflective box, from the closure alone.
    """
    result = _LINEAR_DISCONTINUOUS().face_transmission_spectrum(2)
    assert result.damping is FaceModeDamping.DAMPED
    assert result.spectral_radius is not None
    assert result.spectral_radius < 1.0
    assert result.undetermined_because is None


@pytest.mark.foundation
def test_UNDETERMINED_is_a_third_state_and_it_carries_its_reason():
    """LD at ``ndim = 3`` cannot be driven, and says why.

    ⚠ The gate exists because the dangerous failure is treating this as
    :attr:`FaceModeDamping.DAMPED`.  A caller that assumed "damped"
    would skip the gauge silently on a scheme it never classified —
    exactly the blindness the predicate exists to remove.

    The obstruction is LD's own ``assemble_inflow_axis`` (``axis in
    {0, d-1}`` only), unrelated to this predicate; the reason string
    must name it so a warning can quote the real cause instead of
    guessing.  ``spectral_radius`` is ``None`` — *not measured*, never
    *measured and small* (the ``balance_defect`` discipline).
    """
    result = _LINEAR_DISCONTINUOUS().face_transmission_spectrum(3)
    assert result.damping is FaceModeDamping.UNDETERMINED
    assert result.damping is not FaceModeDamping.DAMPED
    assert result.spectral_radius is None
    assert result.undetermined_because is not None
    assert "LinearDiscontinuous" in result.undetermined_because
    assert "ndim=3" in result.undetermined_because


@pytest.mark.foundation
@pytest.mark.parametrize("ndim", [1, 2, 3])
def test_the_verdict_is_cached_per_closure_and_dimension(ndim):
    """Two instances of one closure share one answer.

    The verdict is a property of the CLOSURE, so it is computed once and
    reused across every mesh, group and iterate — the eigendecomposition
    must not ride the per-solve path.
    """
    assert (
        _DIAMOND().face_transmission_spectrum(ndim)
        is _DIAMOND().face_transmission_spectrum(ndim)
    )


@pytest.mark.foundation
def test_a_cell_dependent_verdict_is_REFUSED_not_averaged():
    """Two probe cells must AGREE, or the answer is ``UNDETERMINED``.

    The predicate claims to describe a closure, not a cell.  This gate
    pins that the claim is TESTED rather than assumed: a stub whose
    damping flips between the two probes is refused, with both radii
    named, instead of being classified from whichever probe ran first
    (`vv-principles` #13 — one sample is not a survey).

    ⚠ Mutation-verified: with the agreement check removed this returns
    ``DAMPED`` (the thin probe's ``0.5``) and the gate reddens.
    """
    from orpheus.transport.spatial import scheme as scheme_module

    calls: list[float] = []

    def _cell_dependent(scheme, ndim, w, sigma_t_volume):
        # rho = 0.5 on the first probe cell, 1.5 on the second.
        rho = 0.5 if sigma_t_volume == 0.9 else 1.5
        calls.append(rho)
        return np.diag([rho])

    original = scheme_module._face_transmission_matrix
    scheme_module._face_transmission_matrix = _cell_dependent
    try:
        # A distinct ndim so the module-level cache cannot serve a real answer.
        result = scheme_module._face_transmission_spectrum(_DIAMOND, 7)
    finally:
        scheme_module._face_transmission_matrix = original
        scheme_module._face_transmission_spectrum.cache_clear()

    assert calls == [0.5, 1.5], "both probe cells must be exercised"
    assert result.damping is FaceModeDamping.UNDETERMINED
    assert result.spectral_radius is None
    assert result.undetermined_because is not None
    assert "probe cells" in result.undetermined_because
    # Both measurements are named, so the reader can see the disagreement.
    assert "0.5000000000" in result.undetermined_because
    assert "1.5000000000" in result.undetermined_because
