#tables.py

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



# sizes (set once by core after reading headers or params)
nlegup = None; npot = None; nlegdim = None
Nrdat  = None; Nvdat = None; Nrbey  = None

# arrays (allocated by tables_read / density loaders / tables())
tabv = None   # (nlegdim, npot)
tabfr = None
tabfv = None
rho = None     # (Nrdat, Nvdat)
rhobey = None  # (Nrbey, Nvdat) used if iextra == -1

def tables_read(path: str | Path):
    """
    Fortran TABLESREAD equivalent:
      do n=0,nlegup,2
         nleg=n/2+1
         do irtab=1,npot
            read -> tabv(nleg, irtab)
            read -> tabfr(nleg, irtab)
            read -> tabfv(nleg, irtab)
    File lines are 'i1 i2 x1' triples; ordering matters: v, fr, fv repeating.
    Returns (tabv, tabfr, tabfv) with shapes (nlegs, npot) using 0-based indexing.
    """
    lines = []
    with open(path, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            # tolerate commas or multiple spaces
            parts = s.replace(",", " ").split()
            if len(parts) < 3:
                continue
            i1, i2 = int(float(parts[0])), int(float(parts[1]))  # accept '1.0' too
            x1 = float(parts[2])
            lines.append((i1, i2, x1))

    if len(lines) % 3 != 0:
        raise ValueError(f"tables.out line count {len(lines)} isn’t a multiple of 3; file ordering may be off")

    max_i1 = max(i for i, _, _ in lines)
    max_i2 = max(j for _, j, _ in lines)
    nlegs = max_i1
    npot  = max_i2

    tabv  = np.zeros((nlegs, npot), dtype=float)
    tabfr = np.zeros((nlegs, npot), dtype=float)
    tabfv = np.zeros((nlegs, npot), dtype=float)

    # cycle v → fr → fv, matching the three READs per irtab in Fortran
    for k, (i1, i2, x1) in enumerate(lines):
        i = i1 - 1
        j = i2 - 1
        sel = k % 3
        if sel == 0:
            tabv[i, j] = x1
        elif sel == 1:
            tabfr[i, j] = x1
        else:
            tabfv[i, j] = x1
    return tabv, tabfr, tabfv


def tables_write(path: str | Path, tabv: np.ndarray, tabfr: np.ndarray, tabfv: np.ndarray) -> None:
    """
    Fortran TABLESWRITE:
      do n=0,nlegup,2
        nleg=n/2+1
        do irtab=1,npot
          write nleg, irtab, tabv(nleg,irtab)
          write nleg, irtab, tabfr(nleg,irtab)
          write nleg, irtab, tabfv(nleg,irtab)
    Python expects tab* arrays with shape (nlegs, npot) using 0-based indexing.
    """
    tabv  = np.asarray(tabv,  dtype=float)
    tabfr = np.asarray(tabfr, dtype=float)
    tabfv = np.asarray(tabfv, dtype=float)
    if tabv.shape != tabfr.shape or tabv.shape != tabfv.shape:
        raise ValueError(f"shape mismatch: {tabv.shape=} {tabfr.shape=} {tabfv.shape=}")
    nlegs, npot = tabv.shape

    path = Path(path)
    with path.open("w") as f:
        for n in range(0, 2*nlegs, 2):             # 0,2,4,... maps to nleg = 1..nlegs
            nleg = n // 2 + 1                       # Fortran 1-based row index
            irow = nleg - 1                         # Python 0-based row
            for irtab in range(1, npot + 1):        # Fortran 1..npot
                j = irtab - 1
                f.write(f"{nleg} {irtab} {tabv[irow, j]:.16e}\n")
                f.write(f"{nleg} {irtab} {tabfr[irow, j]:.16e}\n")
                f.write(f"{nleg} {irtab} {tabfv[irow, j]:.16e}\n")

