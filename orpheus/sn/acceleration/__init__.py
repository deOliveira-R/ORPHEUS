r"""SN iteration-acceleration machinery.

The first member is consistent Diffusion Synthetic Acceleration (#2):
the :class:`~orpheus.sn.acceleration.dsa.DSALowOrderSystem` low-order
half and the :class:`~orpheus.sn.acceleration.dsa.DSACorrection`
operator both the SI and Krylov postures consume. The package exists on
the SN side by the R4 ruling — an accelerator's low-order coefficients
are properties of the SN discretization (quadrature moments, scheme
:math:`\rho`), not of standalone diffusion physics.
"""

from orpheus.sn.acceleration.dsa import DSACorrection, DSALowOrderSystem

__all__ = ["DSACorrection", "DSALowOrderSystem"]
