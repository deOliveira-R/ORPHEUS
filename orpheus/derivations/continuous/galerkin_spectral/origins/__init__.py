r"""Branch-1 SymPy derivations for the Carlvik-Galerkin method.

The single module
:mod:`~orpheus.derivations.continuous.galerkin_spectral.origins.derivations`
carries the algebra-of-record
for the Galerkin spectral expansion of Carlvik's integral equation
:cite:`DahlSjostrand1979`. Every published result that this method produces
is bound, via the foundation-tagged tests at
``tests/derivations/test_galerkin_spectral_symbolic.py``, to a closed-form
SymPy proof in this module.
"""
from __future__ import annotations
