"""Macromolecular Rate Theory (MMRT) temperature response for enzyme kcat.

MMRT (Hobbs et al. 2013, ACS Chem. Biol.) extends transition-state theory with a
non-zero, temperature-independent heat-capacity of activation dCp_ddagger. A
negative dCp gives a single-peaked ln(k) vs T curve, which is exactly the shape
of an enzyme thermal performance curve (rise, optimum, decline).

    ln k(T) = ln(kB*T/h)
              + ( -dH0 - dCp*(T - T0) ) / (R*T)
              + (  dS0 + dCp*(ln T - ln T0) ) / R

with dH0, dS0 the enthalpy/entropy of activation referenced at T0.

Reparameterization (the useful bit for sensitivity analysis)
------------------------------------------------------------
An enzyme-constrained GEM already carries a reference kcat that is valid at some
reference temperature T0 (e.g. the 30 C at which many yeast kcats are reported).
We keep that anchor and describe the *shape* of the curve with two interpretable
knobs the sweeps can perturb:

    Topt  : temperature (K) at which kcat peaks
    dCp   : heat capacity of activation (kJ/mol/K, negative) -> curvature/breadth

Given (kcat_ref, T0, Topt, dCp) we solve closed-form for (dH0, dS0) so that:
    * kcat(T0) == kcat_ref      (base model unchanged at its native temperature)
    * d ln k / dT == 0 at Topt  (curve peaks exactly at Topt)

This means a plain enzyme-constrained model is reproduced exactly when the sweep
sets Topt/dCp to their nominal values, and any perturbation is a controlled
deformation around it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Physical constants (SI-ish, energies in kJ/mol so R in kJ/mol/K)
R = 8.314462618e-3       # kJ / mol / K
KB = 1.380649e-23        # J / K
H = 6.62607015e-34       # J s
_LN_KB_OVER_H = np.log(KB / H)   # ln(kB/h), s^-1 K^-1 prefactor


def _ln_prefactor(T):
    """ln(kB*T/h) transition-state prefactor (T in K)."""
    return _LN_KB_OVER_H + np.log(T)


def solve_dH_dS(kcat_ref, T0, Topt, dCp):
    """Solve activation enthalpy/entropy from the interpretable knobs.

    Parameters
    ----------
    kcat_ref : float   reference turnover (1/s) valid at T0
    T0       : float   reference temperature (K)
    Topt     : float   temperature of peak kcat (K)
    dCp      : float   heat capacity of activation (kJ/mol/K, should be < 0)

    Returns
    -------
    (dH0, dS0) referenced at T0, in kJ/mol and kJ/mol/K.
    """
    # Peak condition, d ln k/dT = 0 at Topt (derivation in module docstring):
    #   1 + (dH0 - dCp*T0)/(R*Topt) + dCp/R = 0
    dH0 = dCp * (T0 - Topt) - R * Topt
    # Value anchor at T0 (the dCp*(T-T0) and dCp*ln terms vanish at T=T0):
    #   ln kcat_ref = ln(kB*T0/h) - dH0/(R*T0) + dS0/R
    dS0 = R * (np.log(kcat_ref) - _ln_prefactor(T0) + dH0 / (R * T0))
    return dH0, dS0


@dataclass
class MMRTParams:
    """MMRT parameters for a single enzyme, anchored on a reference kcat.

    Attributes
    ----------
    kcat_ref : reference turnover (1/s) at T0
    T0       : reference temperature (K)
    Topt     : temperature of peak kcat (K)
    dCp      : heat capacity of activation (kJ/mol/K, negative)
    """
    kcat_ref: float
    T0: float
    Topt: float
    dCp: float

    def dH_dS(self):
        return solve_dH_dS(self.kcat_ref, self.T0, self.Topt, self.dCp)

    def kcat(self, T):
        """kcat (1/s) at temperature T (K); accepts scalars or arrays."""
        T = np.asarray(T, dtype=float)
        dH0, dS0 = self.dH_dS()
        lnk = (
            _ln_prefactor(T)
            + (-dH0 - self.dCp * (T - self.T0)) / (R * T)
            + (dS0 + self.dCp * (np.log(T) - np.log(self.T0))) / R
        )
        return np.exp(lnk)

    def relative_kcat(self, T):
        """kcat(T) / kcat_ref -- the multiplicative factor applied to the
        enzyme-constrained model's native coefficients. Equals 1 at T0."""
        return self.kcat(T) / self.kcat_ref


def relative_kcat_vec(T, T0, Topt, dCp):
    """Vectorised kcat(T)/kcat_ref for arrays of enzyme parameters at scalar T.

    Invariant to kcat_ref (it cancels), so only the shape parameters are needed.
    T0, Topt, dCp are arrays (one per enzyme); T is a scalar temperature (K).

        ln rel = ln(T/T0) + (-dH0 - dCp*(T-T0))/(R*T) + dH0/(R*T0)
                 + dCp*ln(T/T0)/R,     with  dH0 = dCp*(T0-Topt) - R*Topt
    """
    T0 = np.asarray(T0, float)
    Topt = np.asarray(Topt, float)
    dCp = np.asarray(dCp, float)
    dH0 = dCp * (T0 - Topt) - R * Topt
    ln_ratio = np.log(T / T0)
    ln_rel = (ln_ratio
              + (-dH0 - dCp * (T - T0)) / (R * T)
              + dH0 / (R * T0)
              + dCp * ln_ratio / R)
    return np.exp(ln_rel)


def kcat_curve(kcat_ref, T0, Topt, dCp, T):
    """Convenience one-shot evaluation of an MMRT kcat curve."""
    return MMRTParams(kcat_ref, T0, Topt, dCp).kcat(T)


def numeric_Topt(kcat_ref, T0, Topt, dCp, T_lo=273.15, T_hi=353.15, n=4001):
    """Numerically locate the peak of the MMRT curve (K).

    Useful as a sanity check that the analytic anchoring reproduces the
    requested Topt, and for curves whose Topt is pushed outside a physical
    range by extreme parameters.
    """
    grid = np.linspace(T_lo, T_hi, n)
    vals = kcat_curve(kcat_ref, T0, Topt, dCp, grid)
    return float(grid[int(np.argmax(vals))])
