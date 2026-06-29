r"""Analytical reference solutions for the infinite homogeneous medium.

For an infinite homogeneous medium with reflective boundary
conditions the neutron transport equation is degenerate in space —
the flux is spatially flat, and :math:`k_\infty` is the dominant
eigenvalue of :math:`\mathbf{A}^{-1}\mathbf{F}` where

.. math::

    \mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^{T},
    \qquad
    \mathbf{F} = \chi \otimes (\nu\Sigma_f).

This module publishes **both** tiers of reference for 1-group,
2-group, and 4-group test problems:

1. **Legacy** :class:`~orpheus.derivations.common.verification_case.VerificationCase`
   — scalar :math:`k_\infty` only, consumed by the existing legacy
   tests.
2. **Phase-0 continuous**
   :class:`~orpheus.derivations.ContinuousReferenceSolution` —
   carries :math:`k_\infty` **and** the continuous
   :math:`\phi(x, g)` callable (spatially flat, equal to the
   :math:`\ell^{2}`-normalised dominant eigenvector of
   :math:`\mathbf{A}^{-1}\mathbf{F}`) so consumer tests can verify
   flux shape as well as the eigenvalue.

The homogeneous cases are the **simplest possible verification
reference** — no ansatz, no quadrature, no iteration. They land
first in the Phase-1 retrofit plan because they establish the
retrofit pattern (``derive_*`` returns a legacy case alongside
a continuous reference; ``continuous_cases()`` enumerates the
new-form list) without any new numerical machinery.

See :doc:`/verification/reference_solutions` for the contract and
the verification-campaign migration plan.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture

from ...common.eigenvalue import kinf_and_spectrum_homogeneous, kinf_homogeneous
from ...common.continuous_reference import (
    ContinuousReferenceSolution,
    ProblemSpec,
    Provenance,
)
from ...common.verification_case import VerificationCase


def _make_mixture(
    sig_t: np.ndarray,
    sig_c: np.ndarray,
    sig_f: np.ndarray,
    nu: np.ndarray,
    chi: np.ndarray,
    sig_s: np.ndarray,
    sig_2: np.ndarray | None = None,
) -> Mixture:
    """Build a Mixture from N-group arrays.

    The synthetic XS used by analytical-homogeneous derivations are
    abstract — no physical energy grid exists, so ``eg`` is left as
    ``None`` (Phase E of the architectural reset, 2026-05-04).

    ``sig_2`` is the optional ``(ng, ng)`` (n,2n) transfer matrix
    (``[g_from, g_to]``); it defaults to the zero matrix (no (n,2n)),
    which keeps every existing ``derive_*`` call bit-identical. A
    non-zero ``sig_2`` makes the (n,2n) channel — folded into the loss
    matrix as :math:`2\\Sigma_2^T` — materially affect :math:`k_\\infty`,
    de-vacuuming the n2n verification path. Mirrors the ``sig_2=``
    parameter on the sibling :func:`orpheus.derivations.common.xs_library.make_mixture`.
    """
    ng = len(sig_t)
    return Mixture(
        SigC=sig_c.copy(), SigL=np.zeros(ng),
        SigF=sig_f.copy(), SigP=(nu * sig_f).copy(),
        SigT=sig_t.copy(), SigS=[csr_matrix(sig_s)],
        Sig2=csr_matrix(sig_2) if sig_2 is not None else csr_matrix((ng, ng)),
        chi=chi.copy(),
    )


def _build_continuous_homogeneous(
    name: str,
    description: str,
    materials: dict[int, Mixture],
    k_val: float,
    phi_spectrum: np.ndarray,
    provenance: Provenance,
    equation_labels: tuple[str, ...],
) -> ContinuousReferenceSolution:
    r"""Assemble a :class:`ContinuousReferenceSolution` for a homogeneous case.

    The spatial callable ``phi(x, g)`` returns the constant
    ``phi_spectrum[g]`` broadcast over the shape of ``x``. This is
    the honest-to-mathematics continuous solution for the infinite
    homogeneous medium: the flux is spatially flat, so the
    "reference at any point" is literally the same scalar.

    The problem is registered under ``operator_form="homogeneous"``
    — the one operator tag in the taxonomy that is degenerate in
    space and can be consumed by every solver as a sanity check on
    multigroup matrix algebra (which is exactly what the homogeneous
    tests do).
    """
    ng = len(phi_spectrum)

    def phi(x: np.ndarray, g: int = 0) -> np.ndarray:
        return np.full_like(np.asarray(x, dtype=float), phi_spectrum[g])

    return ContinuousReferenceSolution(
        name=name,
        problem=ProblemSpec(
            materials=materials,
            geometry_type="homogeneous",
            geometry_params={},
            boundary_conditions={},  # infinite medium — no boundary
            external_source=None,
            is_eigenvalue=True,
            n_groups=ng,
        ),
        operator_form="homogeneous",
        phi=phi,
        provenance=provenance,
        k_eff=k_val,
        psi=None,  # angular flux is isotropic by symmetry; not needed
        equation_labels=equation_labels,
        vv_level="L1",
        description=description,
        tolerance="< 1e-12",
    )


# ═══════════════════════════════════════════════════════════════════════
# 1-group: k = nu * Sigma_f / Sigma_a  (fully symbolic)
# ═══════════════════════════════════════════════════════════════════════

def derive_1g() -> VerificationCase:
    r"""1-group infinite medium eigenvalue.

    .. math::
        k_\infty = \frac{\nu \Sigma_f}{\Sigma_a}
    """
    nu_s, Sig_f, Sig_a = sp.symbols(r'\nu \Sigma_f \Sigma_a', positive=True)

    k_expr = nu_s * Sig_f / Sig_a

    # Numeric XS
    xs = dict(sig_t=1.0, sig_c=0.2, sig_f=0.3, nu=2.5, sig_s_diag=0.5)
    sig_a_val = xs["sig_c"] + xs["sig_f"]

    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=np.array([xs["sig_t"]]),
        sig_s=np.array([[xs["sig_s_diag"]]]),
        nu_sig_f=np.array([xs["nu"] * xs["sig_f"]]),
        chi=np.array([1.0]),
    )

    latex = (
        r"For a 1-group infinite homogeneous medium, the neutron balance is:"
        "\n\n"
        r".. math::" "\n"
        r"   \Sigma_a \phi = \frac{1}{k} \nu \Sigma_f \phi"
        "\n\n"
        r"which gives the analytical eigenvalue:"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {sp.latex(k_expr)} = "
        rf"\frac{{{xs['nu']} \times {xs['sig_f']}}}"
        rf"{{{sig_a_val}}} = {k_val:.6f}"
    )

    mix = _make_mixture(
        sig_t=np.array([xs["sig_t"]]),
        sig_c=np.array([xs["sig_c"]]),
        sig_f=np.array([xs["sig_f"]]),
        nu=np.array([xs["nu"]]),
        chi=np.array([1.0]),
        sig_s=np.array([[xs["sig_s_diag"]]]),
    )

    return VerificationCase(
        name="homo_1eg",
        k_inf=k_val,
        method="homo",
        geometry="--",
        n_groups=1,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description="1-group infinite medium: k = nu*Sig_f / Sig_a = 1.5",
        tolerance="< 1e-12",
        vv_level="L1",
        equation_labels=("one-group-kinf", "inf-hom-balance"),
    )


# ═══════════════════════════════════════════════════════════════════════
# 2-group: solve characteristic polynomial of 2x2 inv(A)*F
# ═══════════════════════════════════════════════════════════════════════

def derive_2g() -> VerificationCase:
    r"""2-group infinite medium eigenvalue via characteristic polynomial.

    The eigenvalue problem is :math:`\mathbf{A}\phi = \frac{1}{k}\mathbf{F}\phi`
    where :math:`\mathbf{A} = \text{diag}(\Sigma_t) - \Sigma_s^T` and
    :math:`\mathbf{F} = \chi \otimes (\nu\Sigma_f)`.
    """
    # Numeric XS (same as benchmarks.py)
    sig_t = np.array([0.50, 1.00])
    sig_c = np.array([0.01, 0.02])
    sig_f = np.array([0.01, 0.08])
    nu = np.array([2.50, 2.50])
    chi = np.array([1.00, 0.00])
    sig_s = np.array([
        [0.38, 0.10],
        [0.00, 0.90],
    ])

    # Symbolic derivation
    A_sym = sp.Matrix([
        [sig_t[0] - sig_s[0, 0], -sig_s[1, 0]],
        [-sig_s[0, 1], sig_t[1] - sig_s[1, 1]],
    ])
    F_sym = sp.Matrix([
        [chi[0] * nu[0] * sig_f[0], chi[0] * nu[1] * sig_f[1]],
        [chi[1] * nu[0] * sig_f[0], chi[1] * nu[1] * sig_f[1]],
    ])

    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=sig_t, sig_s=sig_s, nu_sig_f=nu * sig_f, chi=chi,
    )

    latex = (
        r"For 2-group infinite medium with downscatter only:"
        "\n\n"
        r".. math::" "\n"
        rf"   \mathbf{{A}} = {sp.latex(A_sym)}, \quad "
        rf"\mathbf{{F}} = {sp.latex(F_sym)}"
        "\n\n"
        r"The eigenvalue is the dominant root of "
        r":math:`\det(\mathbf{A}^{-1}\mathbf{F} - \lambda \mathbf{I}) = 0`:"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_val:.10f}"
    )

    mix = _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s)
    return VerificationCase(
        name="homo_2eg",
        k_inf=k_val,
        method="homo",
        geometry="--",
        n_groups=2,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description="2-group infinite medium (fast + thermal, downscatter only)",
        tolerance="< 1e-12",
        vv_level="L1",
        equation_labels=(
            "matrix-eigenvalue",
            "removal-matrix",
            "fission-matrix",
            "mg-balance",
            "two-group-A",
            "two-group-F",
            "two-group-Ainv",
            "two-group-M",
            "two-group-charpoly",
            "two-group-roots",
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# 4-group: symbolic matrix structure, numeric eigenvalue
# ═══════════════════════════════════════════════════════════════════════

def derive_4g() -> VerificationCase:
    r"""4-group infinite medium eigenvalue.

    Symbolic 4x4 matrix with numeric XS substituted before solving.
    """
    sig_c = np.array([0.01, 0.02, 0.03, 0.05])
    sig_f = np.array([0.005, 0.01, 0.05, 0.10])
    nu = np.array([2.80, 2.60, 2.50, 2.45])
    chi = np.array([0.60, 0.35, 0.05, 0.00])
    sig_s = np.array([
        [0.28, 0.08, 0.02, 0.005],
        [0.00, 0.40, 0.12, 0.06],
        [0.00, 0.00, 0.55, 0.22],
        [0.00, 0.00, 0.00, 0.90],
    ])
    sig_t = sig_c + sig_f + sig_s.sum(axis=1)

    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=sig_t, sig_s=sig_s, nu_sig_f=nu * sig_f, chi=chi,
    )

    latex = (
        r"For 4-group infinite medium (fast :math:`\to` epithermal "
        r":math:`\to` thermal1 :math:`\to` thermal2), "
        r"with downscatter cascade and fission in all groups:"
        "\n\n"
        r".. math::" "\n"
        rf"   \mathbf{{A}} = \text{{diag}}(\Sigma_t) - \Sigma_s^T \in "
        rf"\mathbb{{R}}^{{4 \times 4}}"
        "\n\n"
        r".. math::" "\n"
        rf"   \mathbf{{F}} = \chi \otimes (\nu\Sigma_f) \in "
        rf"\mathbb{{R}}^{{4 \times 4}}"
        "\n\n"
        r"The dominant eigenvalue of :math:`\mathbf{A}^{-1}\mathbf{F}`:"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_val:.10f}"
    )

    mix = _make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s)
    return VerificationCase(
        name="homo_4eg",
        k_inf=k_val,
        method="homo",
        geometry="--",
        n_groups=4,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description="4-group infinite medium (downscatter cascade)",
        tolerance="< 1e-12",
        vv_level="L1",
        equation_labels=(
            "matrix-eigenvalue",
            "removal-matrix",
            "fission-matrix",
            "mg-balance",
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# 2-group with (n,2n): de-vacuums the n2n-in-A convention
# ═══════════════════════════════════════════════════════════════════════

def derive_2g_n2n() -> VerificationCase:
    r"""2-group infinite medium with a non-trivial (n,2n) channel.

    The (n,2n) reaction is a **loss-side multiplicity-2 transfer**: it
    enters the eigenvalue problem ONLY through the loss matrix
    :math:`\mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^T`,
    NOT through the fission production :math:`\mathbf{F} = \chi \otimes
    (\nu\Sigma_f)` — the two emitted neutrons are redistributed by
    :math:`2\Sigma_2` in A, they are not produced with the fission
    spectrum :math:`\chi`.

    This case de-vacuums that convention.  Every other homogeneous
    fixture carries :math:`\Sigma_2 = 0`, so the n2n term is invisible
    (Mode-10 documented-but-unconstrained).  Here :math:`\Sigma_2` is
    asymmetric and non-zero, moving :math:`k_\infty` from ``1.0377``
    (n2n dropped) to ``1.6532`` (correct) — so dropping :math:`2\Sigma_2^T`
    from A, or (the retired bespoke bug) double-counting it into the
    production numerator, reddens the gate far above the FP floor.
    """
    xs = _derive_2g_n2n_inputs()
    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=xs["sig_t"], sig_s=xs["sig_s"],
        nu_sig_f=xs["nu"] * xs["sig_f"], chi=xs["chi"], sig_2=xs["sig_2"],
    )
    latex = (
        r"For a 2-group infinite medium with a non-trivial (n,2n) channel, "
        r"the (n,2n) reaction folds into the loss matrix as a "
        r"multiplicity-2 transfer (it is NOT a production channel):"
        "\n\n"
        r".. math::" "\n"
        r"   \mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^T, "
        r"\quad \mathbf{F} = \chi \otimes (\nu\Sigma_f)"
        "\n\n"
        r"The eigenvalue is the dominant root of "
        r":math:`\det(\mathbf{A}^{-1}\mathbf{F} - \lambda \mathbf{I}) = 0`:"
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_val:.10f}"
    )
    mix = _make_mixture(
        xs["sig_t"], xs["sig_c"], xs["sig_f"],
        xs["nu"], xs["chi"], xs["sig_s"], sig_2=xs["sig_2"],
    )
    return VerificationCase(
        name="homo_2eg_n2n",
        k_inf=k_val,
        method="homo",
        geometry="--",
        n_groups=2,
        n_regions=1,
        materials={0: mix},
        geom_params={},
        latex=latex,
        description=(
            "2-group infinite medium with non-trivial (n,2n) — "
            "de-vacuums the n2n-in-A convention"
        ),
        tolerance="< 1e-12",
        vv_level="L1",
        equation_labels=(
            "matrix-eigenvalue",
            "removal-matrix",
            "fission-matrix",
            "mg-balance",
        ),
    )


def all_cases() -> list[VerificationCase]:
    """Return all homogeneous legacy :class:`VerificationCase` objects."""
    return [derive_1g(), derive_2g(), derive_4g(), derive_2g_n2n()]


# ═══════════════════════════════════════════════════════════════════════
# Phase-0 continuous references
# ═══════════════════════════════════════════════════════════════════════
#
# Each ``derive_Ng_continuous`` builds the same mathematical problem
# as the legacy ``derive_Ng`` above but returns a
# :class:`ContinuousReferenceSolution` carrying both ``k_eff`` and
# the flat spatial flux ``phi(x, g) = phi_spectrum[g]``. The
# spectrum is the :math:`\ell^{2}`-normalised dominant eigenvector
# of :math:`\mathbf{A}^{-1}\mathbf{F}` computed inside
# :func:`kinf_and_spectrum_homogeneous`.
#
# Intentional duplication with the legacy ``derive_Ng``: Phase 1.1
# keeps the legacy call paths untouched so nothing in
# ``tests/homogeneous/`` has to change. Phase 2 may fold them into
# the continuous derivations once every consumer has migrated.


def _derive_1g_inputs():
    """Shared cross sections for the 1-group homogeneous reference."""
    return dict(sig_t=1.0, sig_c=0.2, sig_f=0.3, nu=2.5, sig_s_diag=0.5)


def _derive_2g_inputs():
    """Shared cross sections for the 2-group homogeneous reference."""
    return dict(
        sig_t=np.array([0.50, 1.00]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.50, 2.50]),
        chi=np.array([1.00, 0.00]),
        sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
    )


def _derive_2g_n2n_inputs():
    """2-group homogeneous reference WITH a non-trivial (n,2n) channel.

    Asymmetric ``sig_2`` (rowsum != colsum) so the (n,2n) convention is
    fully exercised: the loss-matrix term :math:`2\\Sigma_2^T` AND the
    rowsum-vs-colsum distinction both matter.  ``sig_t`` is balanced
    (total = capture + fission + scatter-out + n2n-out, per row).
    """
    sig_c = np.array([0.03, 0.06])
    sig_f = np.array([0.012, 0.10])
    nu = np.array([2.50, 2.45])
    chi = np.array([1.00, 0.00])
    sig_s = np.array([[0.45, 0.10], [0.00, 0.82]])
    sig_2 = np.array([[0.010, 0.020], [0.000, 0.005]])  # asymmetric
    sig_t = sig_c + sig_f + sig_s.sum(axis=1) + sig_2.sum(axis=1)
    return dict(
        sig_t=sig_t, sig_c=sig_c, sig_f=sig_f,
        nu=nu, chi=chi, sig_s=sig_s, sig_2=sig_2,
    )


def _derive_4g_inputs():
    """Shared cross sections for the 4-group homogeneous reference."""
    sig_c = np.array([0.01, 0.02, 0.03, 0.05])
    sig_f = np.array([0.005, 0.01, 0.05, 0.10])
    nu = np.array([2.80, 2.60, 2.50, 2.45])
    chi = np.array([0.60, 0.35, 0.05, 0.00])
    sig_s = np.array([
        [0.28, 0.08, 0.02, 0.005],
        [0.00, 0.40, 0.12, 0.06],
        [0.00, 0.00, 0.55, 0.22],
        [0.00, 0.00, 0.00, 0.90],
    ])
    sig_t = sig_c + sig_f + sig_s.sum(axis=1)
    return dict(
        sig_t=sig_t, sig_c=sig_c, sig_f=sig_f,
        nu=nu, chi=chi, sig_s=sig_s,
    )


def derive_1g_continuous() -> ContinuousReferenceSolution:
    r"""Phase-0 continuous reference for the 1-group homogeneous medium.

    The 1-group case is the single cleanest verification reference
    in the whole campaign: the eigenvalue is
    :math:`k_\infty = \nu\Sigma_f / \Sigma_a` in closed form, the
    flux is spatially flat and — with only one group — identically
    equal to 1 after :math:`\ell^{2}` normalisation. A consumer
    test that does not reproduce this to machine precision has a
    bug in its multigroup balance, not in the discretisation.
    """
    xs = _derive_1g_inputs()
    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=np.array([xs["sig_t"]]),
        sig_s=np.array([[xs["sig_s_diag"]]]),
        nu_sig_f=np.array([xs["nu"] * xs["sig_f"]]),
        chi=np.array([1.0]),
    )
    mix = _make_mixture(
        sig_t=np.array([xs["sig_t"]]),
        sig_c=np.array([xs["sig_c"]]),
        sig_f=np.array([xs["sig_f"]]),
        nu=np.array([xs["nu"]]),
        chi=np.array([1.0]),
        sig_s=np.array([[xs["sig_s_diag"]]]),
    )
    return _build_continuous_homogeneous(
        name="homo_1eg",
        description=(
            "1-group infinite medium — k_inf = nu·Sigma_f / Sigma_a. "
            "Phase-0 continuous reference."
        ),
        materials={0: mix},
        k_val=k_val,
        phi_spectrum=phi_spectrum,
        provenance=Provenance(
            citation="Bell & Glasstone 1970, §1.5 Eq. (1.58)",
            derivation_notes=(
                "One-group infinite homogeneous medium with reflective BC. "
                "The transport equation reduces to the scalar balance "
                "Sigma_a·phi = (1/k)·nu·Sigma_f·phi, giving the closed-form "
                "k_inf = nu·Sigma_f / Sigma_a. With a single group the "
                "spectrum is degenerate (phi = 1 after l2 normalisation)."
            ),
            sympy_expression=r"k_\infty = \nu \Sigma_f / \Sigma_a",
            precision_digits=None,  # closed form
        ),
        equation_labels=("one-group-kinf", "inf-hom-balance"),
    )


def derive_2g_continuous() -> ContinuousReferenceSolution:
    r"""Phase-0 continuous reference for the 2-group homogeneous medium.

    The flux spectrum is the non-trivial dominant eigenvector of
    the 2×2 matrix :math:`\mathbf{A}^{-1}\mathbf{F}`. This is the
    **first** reference where a consumer test has something to
    verify beyond :math:`k_{\text{eff}}` — namely the fast/thermal
    split of the converged multigroup flux.
    """
    xs = _derive_2g_inputs()
    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=xs["sig_t"], sig_s=xs["sig_s"],
        nu_sig_f=xs["nu"] * xs["sig_f"], chi=xs["chi"],
    )
    mix = _make_mixture(
        xs["sig_t"], xs["sig_c"], xs["sig_f"],
        xs["nu"], xs["chi"], xs["sig_s"],
    )
    return _build_continuous_homogeneous(
        name="homo_2eg",
        description=(
            "2-group infinite medium (fast + thermal, downscatter only). "
            "Phase-0 continuous reference."
        ),
        materials={0: mix},
        k_val=k_val,
        phi_spectrum=phi_spectrum,
        provenance=Provenance(
            citation="Bell & Glasstone 1970, §7.4",
            derivation_notes=(
                "Two-group reflective infinite medium. The eigenvalue "
                "problem is A·phi = (1/k)·F·phi with A = diag(Sigma_t) - "
                "Sigma_s^T and F = chi ⊗ (nu·Sigma_f). Solved via "
                "characteristic polynomial (closed form for 2x2); "
                "the l2-normalised dominant eigenvector is the flux "
                "spectrum."
            ),
            sympy_expression=None,
            precision_digits=None,
        ),
        equation_labels=(
            "matrix-eigenvalue", "removal-matrix", "fission-matrix",
            "mg-balance", "two-group-A", "two-group-F", "two-group-Ainv",
            "two-group-M", "two-group-charpoly", "two-group-roots",
        ),
    )


def derive_4g_continuous() -> ContinuousReferenceSolution:
    r"""Phase-0 continuous reference for the 4-group homogeneous medium.

    Downscatter cascade with fission in all groups — the dominant
    eigenvector is the first genuinely non-trivial flux spectrum
    in the Phase-0 registry (1g degenerate, 2g has only two
    components). Used to detect bugs in multigroup chi-weighting
    that single out higher-order cascades.
    """
    xs = _derive_4g_inputs()
    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=xs["sig_t"], sig_s=xs["sig_s"],
        nu_sig_f=xs["nu"] * xs["sig_f"], chi=xs["chi"],
    )
    mix = _make_mixture(
        xs["sig_t"], xs["sig_c"], xs["sig_f"],
        xs["nu"], xs["chi"], xs["sig_s"],
    )
    return _build_continuous_homogeneous(
        name="homo_4eg",
        description=(
            "4-group infinite medium (downscatter cascade, fission "
            "in all groups). Phase-0 continuous reference."
        ),
        materials={0: mix},
        k_val=k_val,
        phi_spectrum=phi_spectrum,
        provenance=Provenance(
            citation="Bell & Glasstone 1970, §7.4 (multigroup form)",
            derivation_notes=(
                "Four-group reflective infinite medium with downscatter "
                "cascade. Dense 4x4 A^{-1}·F solved via numpy.linalg.eig; "
                "the dominant eigenvector gives the non-trivial flux "
                "cascade shape."
            ),
            sympy_expression=None,
            precision_digits=None,
        ),
        equation_labels=(
            "matrix-eigenvalue", "removal-matrix", "fission-matrix",
            "mg-balance",
        ),
    )


def derive_2g_n2n_continuous() -> ContinuousReferenceSolution:
    r"""Phase-0 continuous reference for the 2-group (n,2n) homogeneous medium.

    The de-vacuum companion to :func:`derive_2g_continuous`: the dominant
    eigenvector (flux spectrum) of :math:`\mathbf{A}^{-1}\mathbf{F}` with
    the (n,2n) term :math:`2\Sigma_2^T` active in A.  A transpose flip or
    axis bug in the n2n term moves the eigenVECTOR (not just the
    eigenvalue), so this case catches it on the spectrum.
    """
    xs = _derive_2g_n2n_inputs()
    k_val, phi_spectrum = kinf_and_spectrum_homogeneous(
        sig_t=xs["sig_t"], sig_s=xs["sig_s"],
        nu_sig_f=xs["nu"] * xs["sig_f"], chi=xs["chi"], sig_2=xs["sig_2"],
    )
    mix = _make_mixture(
        xs["sig_t"], xs["sig_c"], xs["sig_f"],
        xs["nu"], xs["chi"], xs["sig_s"], sig_2=xs["sig_2"],
    )
    return _build_continuous_homogeneous(
        name="homo_2eg_n2n",
        description=(
            "2-group infinite medium with non-trivial (n,2n). Phase-0 "
            "continuous reference — de-vacuums the n2n-in-A convention."
        ),
        materials={0: mix},
        k_val=k_val,
        phi_spectrum=phi_spectrum,
        provenance=Provenance(
            citation=(
                "Bell & Glasstone 1970, §7.4 + §1.5 "
                "((n,2n) as a loss-side multiplicity transfer)"
            ),
            derivation_notes=(
                "Two-group reflective infinite medium with a non-trivial "
                "(n,2n) transfer. The (n,2n) reaction is a loss-side "
                "multiplicity-2 channel: A = diag(Sigma_t) - (Sigma_s + "
                "2·Sigma_2)^T, F = chi ⊗ (nu·Sigma_f). Production carries "
                "NO (n,2n) term — the two emitted neutrons are redistributed "
                "by 2·Sigma_2 in A, not produced with the fission spectrum. "
                "Asymmetric Sigma_2 (rowsum != colsum) so a transpose or "
                "axis bug is detectable on the eigenvector."
            ),
            sympy_expression=None,
            precision_digits=None,
        ),
        equation_labels=(
            "matrix-eigenvalue", "removal-matrix", "fission-matrix",
            "mg-balance",
        ),
    )


def continuous_cases() -> list[ContinuousReferenceSolution]:
    """Return all homogeneous Phase-0 continuous reference solutions."""
    return [
        derive_1g_continuous(),
        derive_2g_continuous(),
        derive_4g_continuous(),
        derive_2g_n2n_continuous(),
    ]
