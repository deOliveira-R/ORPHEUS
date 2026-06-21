r"""Self-validating value-object for the fission emission spectrum :math:`\chi`.

The fission emission spectrum :math:`\chi_g` is the probability that a
neutron born in fission appears in energy group :math:`g`. As a probability
distribution it carries the **simplex law**

.. math::

    \sum_g \chi_g = 1, \qquad \chi_g \ge 0 \;\; \forall g,

for a *fissile* material. A *non-fissile* material emits no fission
neutrons, so its spectrum is identically zero (the **null spectrum**).

``EmissionSpectrum`` is a :class:`numpy.ndarray` SUBCLASS that carries
these laws as methods (``assert_normalized`` / ``is_emitting`` /
``assert_null``) — so the invariant reads like a law the data carries,
the way an operator carries its adjoint. Every existing call site that
treats :math:`\chi` as a bare array (``chi[None, :, :]`` per-cell
broadcast, ``chi.sum()``, ``chi.copy()``, einsum-feeding the fission
source) keeps working with zero ripple, because an ``EmissionSpectrum``
*is* an ndarray.

Scope (what this type IS and IS NOT)
------------------------------------
This is a SOURCE-LEVEL validated value-object: the simplex / null law is
enforced ONCE, at the data source (``Mixture`` / ``Isotope``
``__post_init__``), NEVER per-cell. It is deliberately NOT a member of
the flux/operator algebra — it has no ``__add__`` gate, no cross-mesh
partner check, no change-of-basis morphism. Contrast the (#263-deferred)
per-cell ``SpectrumField`` that WOULD live in the field-typed algebra:
the per-cell broadcast of :math:`\chi` legitimately holds :math:`\chi=0`
in non-fissile cells, so the simplex law is meaningless there. The law
belongs to the *material*, not to the *cell*.

Tolerance rationale
-------------------
Real GENDF data (``load_isotope('U_235', 294)``) gives
``chi.sum() == 1.0000000000000002`` — an FP residual of :math:`\approx
2.22\times10^{-16}`. ``assert_normalized`` uses ``atol=1e-12`` (~4
orders of magnitude of headroom over the observed residual) while still
flagging a :math:`10^{-6}`-scale physical normalization error.
"""

from __future__ import annotations

import numpy as np


class EmissionSpectrum(np.ndarray):
    r"""A validated fission emission spectrum :math:`\chi` (probability simplex).

    A thin :class:`numpy.ndarray` subclass that carries the simplex /
    null law of the fission spectrum as inspectable methods. Constructed
    from any array-like; the values are NOT modified (this type wraps the
    values, it does not move them).

    Examples
    --------
    >>> chi = EmissionSpectrum(np.array([0.6, 0.35, 0.05]))
    >>> chi.assert_normalized()        # passes: simplex
    >>> chi.is_emitting
    True
    >>> EmissionSpectrum(np.zeros(3)).assert_null()   # passes: null spectrum
    """

    def __new__(cls, values) -> "EmissionSpectrum":
        # ``view(cls)`` reinterprets the existing buffer as this subclass
        # without copying. ``np.asarray`` is a no-op when ``values`` is
        # already an ndarray, so existing arrays are wrapped in place.
        return np.asarray(values).view(cls)

    def __array_finalize__(self, obj) -> None:
        # Called on every construction route (explicit ``__new__``, view,
        # slice, ufunc output). No extra instance state to propagate —
        # the simplex/null law lives in the validators, not in carried
        # attributes — so this is intentionally a no-op. It exists so that
        # slicing / copying / arithmetic return a usable array.
        if obj is None:
            return

    def assert_normalized(self) -> None:
        r"""Raise unless :math:`\sum_g \chi_g \approx 1` AND :math:`\chi_g \ge 0`.

        The two clauses are checked INDEPENDENTLY so a negative-but-
        sums-to-one spectrum (e.g. ``[1.2, -0.2]``) still raises on the
        non-negativity clause.

        Raises
        ------
        ValueError
            If the spectrum is not a valid probability simplex, with the
            actual sum reported in the message.
        """
        total = self.sum()
        if not np.isclose(total, 1.0, atol=1e-12, rtol=0):
            raise ValueError(
                f"EmissionSpectrum is not normalized: sum(chi) = {total!r} "
                f"(expected 1.0 within atol=1e-12)"
            )
        if not np.all(self >= 0):
            raise ValueError(
                f"EmissionSpectrum has negative entries: chi = "
                f"{np.asarray(self)!r} (a probability simplex requires "
                f"chi_g >= 0 for all g)"
            )

    @property
    def is_emitting(self) -> bool:
        r"""``True`` iff the spectrum is non-null (any :math:`\chi_g > 0`).

        A SPECTRUM-level query about the array's own values — it asks
        whether *this* spectrum carries any probability mass, not whether
        the owning material produces fission neutrons. The two coincide
        for a correctly-constructed spectrum, but they are distinct
        concepts: production is a material predicate (``νΣ_f > 0``), and
        the simplex/null law (enforced by :func:`enforce_emission_spectrum`)
        is what ties the producing material to a non-null spectrum.
        """
        return bool(np.any(self > 0))

    def assert_null(self) -> None:
        r"""Raise unless :math:`\chi_g \equiv 0` for every group.

        The non-fissile contract: a non-fissile material emits no fission
        neutrons, so its spectrum is identically zero. These are
        CONSTRUCTED zeros (not computed), so the check is for EXACT zero —
        there is no physical residual to tolerate.

        Raises
        ------
        ValueError
            If any entry is nonzero, with the offending spectrum reported.
        """
        if not np.all(self == 0):
            raise ValueError(
                f"EmissionSpectrum is not null: chi = {np.asarray(self)!r} "
                f"(a non-fissile material must carry an identically-zero "
                f"fission spectrum)"
            )


def enforce_emission_spectrum(chi, *, is_producing: bool) -> EmissionSpectrum:
    r"""Coerce ``chi`` to :class:`EmissionSpectrum` and enforce its law at the source.

    A *producing* material (one that emits fission neutrons, :math:`\nu\Sigma_f
    > 0`) carries a probability simplex; a non-producing material carries the
    null spectrum. The emission spectrum :math:`\chi` is consumed only as a
    fission *source* — the birth density :math:`\chi_g\,\nu\Sigma_{f,g'}\,
    \varphi_{g'}` — so a valid (normalized) spectrum is required exactly where
    production is nonzero, and the null spectrum exactly where it is zero.

    Single-sources the simplex/null law shared by the ``Mixture`` and
    ``Isotope`` ``__post_init__`` (Pattern 2): each container supplies its own
    ``is_producing`` predicate (``νΣ_f > 0``) and routes its ``chi`` through
    this one law.

    Parameters
    ----------
    chi : array-like
        The fission emission spectrum to validate. Wrapped (not copied) into
        an :class:`EmissionSpectrum`.
    is_producing : bool
        ``True`` iff the owning material emits fission neutrons
        (:math:`\nu\Sigma_f > 0`).

    Returns
    -------
    EmissionSpectrum
        The coerced, validated spectrum (same buffer as ``chi``).

    Raises
    ------
    ValueError
        If a producing material's spectrum is not a probability simplex, or a
        non-producing material's spectrum is not identically zero.
    """
    chi = EmissionSpectrum(chi)
    chi.assert_normalized() if is_producing else chi.assert_null()
    return chi
