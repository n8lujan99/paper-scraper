#core.py
# The Control Hub

from __future__ import annotations
import os
# importing all constants
import math
from pathlib import Path
from datetime import datetime
import numpy as np
from numpy.polynomial.legendre import legval
from dataclasses import dataclass
from astropy.coordinates import SkyCoord 
import astropy.units as u
import pandas as pd
import numpy as np
#THE ONE TRUE FILE
from libshim import register_halo_grid, register_density_grid, fac1
import tables_state as S
from tables import tables, tables_read, tables_write  # etc.
from potential import potential, force
from io_utils import iterread, mlread, weight_read
from orbits import run_single_orbit

'''
ospm Table of Contents
  libshim.py              # small: constants, fac1, grid registration, Pl(), mappers
  tables.py               # tables(), denny(), tables_read/write (uses quadpack & shim)
  potential.py            # potential(), force(), potwrite(), (reads tables_state & shim)
  orbits.py               # derivs(), rk4_step(), step(), energy(), ekin(), run_single_orbit()
  io_utils.py             # iterread, mlread, phasewrite/read, weight_read/write, xlib_read, vaniread
  linear_model.py         # build_cmatrix(), new_cmatrix_update(), summ()
'''

def main():
    # 1) load grids
    rhalo = ...   # load from file
    rgrid = ...   # density r grid
    vgrid = ...   # density |sinθ| grid
    register_halo_grid(rhalo)
    register_density_grid(rgrid, vgrid)

    # 2) sizes
    S.nlegup = 8
    S.nlegdim = S.nlegup//2 + 1
    S.npot = len(rhalo)
    S.Nrdat, S.Nvdat = len(rgrid), len(vgrid)

    # 3) allocate tables
    S.tabv  = np.zeros((S.nlegdim, S.npot))
    S.tabfr = np.zeros_like(S.tabv)
    S.tabfv = np.zeros_like(S.tabv)

    # 4) build or read tables
    # tables()  # or tables_read("tables.out")

    # 5) use them
    v = 0.5
    print("Phi(r=1,v=0.5) =", potential(1.0, v))
    stats = run_single_orbit(0.5, 5.0, 0.2, Nsos2=50, Nskip=10)
    print(stats)

if __name__ == "__main__":
    main()
