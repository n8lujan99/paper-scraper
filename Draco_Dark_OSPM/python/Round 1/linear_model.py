#linear_model.py

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

def build_cmatrix(xlib: np.ndarray, v1lib: np.ndarray) -> np.ndarray:
    """
    xlib:  shape (Nbin, Norbit), float64
    v1lib: shape (Nvel, Nvelb, Norbit), float64   # matches Fortran v1lib(ivel, irc, i)
    returns Cm with shape (Nbin + Nvel*Nvelb, Norbit + Nvel*Nvelb)
    """
    Nbin, Norbit = xlib.shape
    Nvel, Nvelb, _ = v1lib.shape
    Nvtot = Nvel * Nvelb
    lda = Nbin + Nvtot
    lds = Norbit + Nvtot
    Cm = np.zeros((lda, lds), dtype=np.float64)
    Cm[:Nbin, :Norbit] = xlib

    # lower-left block: stack v1lib in (irc, ivel) order
    # Fortran row index: Nbin + (irc-1)*Nvel + ivel  with 1-based ivel
    # NumPy zero-based:  Nbin + irc*Nvel + ivel
    row = Nbin
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            Cm[row, :Norbit] = v1lib[ivel, irc, :]
            row += 1

    # lower-right identity
    Cm[Nbin:, Norbit:] = np.eye(Nvtot, dtype=np.float64)
    return Cm

def new_cmatrix_update(Cm: np.ndarray, xlib: np.ndarray, v1lib: np.ndarray, inew_f: int, Norbit_f: int,
                       Nvel: int, Nvelb: int, Nbin: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Mirrors Fortran newcmatrix:
      - swap orbit column 'inew' with column 'Norbit'
      - write that column into Cm
      - Norbit = Norbit - 1
      - zero right block and set new identity at lower-right for size Nvel*Nvelb
    Returns updated (Cm, xlib, v1lib, Norbit_f_new)
    """
    # 1→0 based
    inew = inew_f - 1
    Norbit = Norbit_f
    last  = Norbit - 1
    Nvtot = Nvel * Nvelb
    lda   = Nbin + Nvtot
    lds   = Norbit + Nvtot

    assert xlib.shape == (Nbin, Norbit)
    assert v1lib.shape == (Nvel, Nvelb, Norbit)
    assert Cm.shape[0] == lda and Cm.shape[1] == lds

    # swap xlib(:, inew) <-> xlib(:, last)
    xlib[:, [inew, last]] = xlib[:, [last, inew]]

    # swap v1lib(:,:, inew) <-> v1lib(:,:, last)
    v1lib[:, :, [inew, last]] = v1lib[:, :, [last, inew]]

    # write upper-left column Cm(1:Nbin, inew) = xlib(:, inew)
    Cm[:Nbin, inew] = xlib[:, inew]

    # write lower-left rows for this column, row order j = Nbin + (irc-1)*Nvel + ivel
    row = Nbin
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            Cm[row, inew] = v1lib[ivel, irc, inew]
            row += 1

    # shrink active orbit count
    Norbit_new = Norbit - 1

    # zero the entire right block for the new layout
    Cm[:, Norbit_new:Norbit_new + Nvtot] = 0.0

    # place identity in lower-right of size Nvtot at the new column start
    # rows Nbin: end, cols Norbit_new: end
    Cm[Nbin:, Norbit_new:Norbit_new + Nvtot] = np.eye(Nvtot, dtype=Cm.dtype)

    return Cm, xlib, v1lib, Norbit_new

def summ(w: np.ndarray,
         xlib: np.ndarray,          # (Nbin, Norbit)
         v1lib: np.ndarray,         # (Nvel, Nvelb, Norbit) → v1lib[ivel, irc, iorb]
         cast_to_float32: bool = True) -> np.ndarray:
    """
    Fortran SUMM: summed(1..Nbin)   = Σ_iorb w(iorb) * xlib(ibin, iorb)
                   summed(Nbin+..)  = Σ_iorb w(iorb) * v1lib(ivel, irc, iorb)
    Storage order for the second block matches: i = Nbin + (irc-1)*Nvel + ivel.
    """
    w = np.asarray(w, dtype=np.float64)
    xlib = np.asarray(xlib, dtype=np.float64)
    v1lib = np.asarray(v1lib, dtype=np.float64)

    Nbin, Norbit = xlib.shape
    Nvel, Nvelb, Norbit_v = v1lib.shape
    if Norbit_v != Norbit or w.shape != (Norbit,):
        raise ValueError(f"shape mismatch: w{w.shape}, xlib{ xlib.shape }, v1lib{ v1lib.shape }")

    x_part = xlib @ w
    v_mat  = np.tensordot(v1lib, w, axes=([2], [0]))        # (Nvel, Nvelb)
    v_flat = v_mat.ravel(order='F')                          # irc outer, ivel inner

    out = np.concatenate([x_part, v_flat])
    return out.astype(np.float32) if cast_to_float32 else out