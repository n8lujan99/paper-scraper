#potential.py

from __future__ import annotations
import os
import math
from pathlib import Path
from datetime import datetime
import numpy as np
from numpy.polynomial.legendre import legval
from dataclasses import dataclass
from astropy.coordinates import SkyCoord 
import astropy.units as u
import pandas as pd
#THE ONE TRUE FILE
from libshim import register_halo_grid, register_density_grid, fac1
import tables_state as S
from tables import tables, tables_read, tables_write  # etc.
from potential import potential, force
from io_utils import iterread, mlread, weight_read
from orbits import run_single_orbit

def force

def RneeIRhalo(ir: int) -> float:
    """Log-spaced halo r grid from libdefs.h (1-based ir)."""
    if npot <= 1:
        raise ValueError("npot not initialized")
    t = (rlgpotmax - rlgpotmin) * (ir - 1) / float(npot - 1)
    return 10.0 ** (rlgpotmin + t)

def _P(n: int, x: float) -> float:
    c = np.zeros(n + 1, dtype=float); c[n] = 1.0
    return float(legval(x, c))

def potential(r: float, v: float) -> float:
    # rir is the floating radial index on [1, npot]
    rir = (math.log10(r) - rlgpotmin) * (float(npot - 1)) / (rlgpotmax - rlgpotmin) + 1.0
    irlo = int(rir)           # floor in Fortran via INT()
    irup = int(rir + 1.0)
    rlo = RneeIRhalo(irlo_cl)
    rup = RneeIRhalo(irup_cl)
    VVlo, VVup = 0.0, 0.0
    # Fortran index: nleg = n/2 + 1   → Python zero-based row = n//2
    for n in range(0, nlegup + 1, 2):
        row = n // 2
        Pn  = _P(n, v)

        # Only the table lookups are clamped; the r interpolation isn’t.
        if irlo == 0:
            term_lo = 0.0
        else:
            jlo = max(1, min(irlo, npot)) - 1
            term_lo = Pn * tabv[row, jlo]
        jup = max(1, min(irup, npot)) - 1
        term_up = Pn * tabv[row, jup]

        VVlo += term_lo
        VVup += term_up

    VV = VVlo + (VVup - VVlo) * ((r - rlo) / (rup - rlo)) if rup != rlo else VVlo
    VV = -2.0 * math.pi * VV - hole / totlight / gdennorm / r
    return VV
