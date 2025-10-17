#I/O_utils.py

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

# Placeholders

'''
def projmass()
def halodens()
def haloread()
def halotables()
def halowrite()
'''

def iterread(path: Path = Path("iter.params")) -> IterParams:
    """Read iteration parameters from iter.params."""
    with open(path, "r") as f:
        lines = f.readlines()

    # Fortran discards the dummy label and keeps the second token
    vals = [line.split()[1] for line in lines if len(line.split()) >= 2]

    return IterParams(
        Niter=int(vals[0]),
        apfac=float(vals[1]),
        apfacmu=float(vals[2]),
        ifit=int(vals[3]),
        alphainc=float(vals[4]),
        alphastop=float(vals[5]),
        fracrv=float(vals[6]),
    )
#itp = iterread()
#print(itp.Niter, itp.apfac)

def mlread(path: Path = Path("ml.out")) -> MLParams:
    """Read ML parameters from ml.out and apply Fortran's xmu transformation."""
    with open(path, "r") as f:
        lines = [line.strip().split() for line in f if line.strip()]

    # Fortran overwrites ent, chi twice; we keep the last read values.
    xmu = float(lines[0][0])
    alp, ent, chi = map(float, lines[1])
    x1, x2, ent2, chi2 = map(float, lines[2])

    # Apply the Fortran transform
    xmu = 1.0 / math.sqrt(xmu)

    return MLParams(xmu, alp, ent2, chi2, x1, x2)    

#ml = mlread()
#print(ml.xmu, ml.alp, ml.chi)


def phasewrite(wphase: np.ndarray, norbtot: int, path: Path = Path("phase.out")) -> None:
    """Write wphase[0:norbtot] to phase.out exactly once, like Fortran."""
    np.savetxt(path, wphase[:norbtot][None, :])  # single line, space-separateddef phasewrite()

def phaseread(norbit_path: Path = Path("norbit"), phase_path: Path = Path("phase.out"),
              Norbm: int = None, Norbitm: int = None) -> PhaseReadResult:
    """
    Read Norb from 'norbit' and weights from 'phase.out'.
    If Norbm != Norbitm, duplicate each weight as in Fortran to build size Norbit=2*Norb.
    Otherwise keep Norbit=Norb.
    """
    with open(norbit_path, "r") as f:
        Norb = int(f.read().strip().split()[0])

    slush = np.loadtxt(phase_path, ndmin=1)
    # tolerate either a row or column; flatten like Fortran read
    slush = np.ravel(slush)
    slush = slush[:Norb]

    if Norbm is not None and Norbitm is not None and Norbm != Norbitm:
        # Fortran: if(Norbm.ne.Norbitm) then Norbit=2*Norb and duplicate each weight
        Norbit = 2 * Norb
        wphase_model = np.repeat(slush, 2)
    else:
        Norbit = Norb
        wphase_model = slush.copy()
    imin = int(np.argmin(wphase_model)) + 1  # 1-based like the print in Fortran
    imax = int(np.argmax(wphase_model)) + 1
    wmin = float(wphase_model[imin - 1])
    wmax = float(wphase_model[imax - 1])
    return PhaseReadResult(wphase=wphase_model, Norb=Norb, Norbit=Norbit, wmin=wmin, wmax=wmax, imin=imin, imax=imax)


def potwrite(b: float, rpotmin: float, rpotmax: float, fac1: float) -> None:
    ddd = math.log(3.6 / b) / 999.0
    ccc = b / 3.0 / math.exp(ddd)
    rs   = np.array([ccc * math.exp(ddd * irr) for irr in range(1, 1001)], dtype=float)
    vv0  = np.array([potential(r, 0.0) for r in rs], dtype=float)
    vv5  = np.array([potential(r, 0.5) for r in rs], dtype=float)
    vv9  = np.array([potential(r, 0.9) for r in rs], dtype=float)
    np.savetxt("potential.out", np.column_stack([rs, vv0, vv5, vv9]), fmt="%.6e")

    r_fixed = 0.5
    th = (np.arange(0, 900, dtype=float) / 10.0) * (math.pi / 180.0)
    v  = np.sin(th)
    V  = np.array([potential(r_fixed, vi) for vi in v], dtype=float)
    deg = th * 180.0 / math.pi
    with open("potential.out", "ab") as f:
        np.savetxt(f, np.column_stack([deg, -V]), fmt="%.6e")

    ipnts = 200
    xr = 10 ** np.linspace(math.log10(rpotmin), math.log10(rpotmax), ipnts)
    Fr = np.array([force(r, 0.0)[0] for r in xr], dtype=float)
    vc = np.sqrt(np.abs(Fr) * xr) * fac1
    np.savetxt("vcirc.out", np.column_stack([np.arange(1, ipnts + 1), xr, vc]), fmt=["%d", "%.6e", "%.6e"])


##########################################################################################################################

def weight_read(path: str | Path, Norb: int, Norbit: int, Nvel: int, Nvelb: int) -> np.ndarray:
    """
    Fortran weightread():
      ilast = Norbit + Nvel*Nvelb
      try to read (ii, wt(ii)) for i=1..ilast
      if full: w(1:ilast) = wt(1:ilast)
      if EOF:  for i=1..Norb: w(2*i-1)=wt(i)/2; w(2*i)=wt(i)/2
               for i=1..Nvel*Nvelb+1: w(Norbit+i)=wt(Norb+i)
    Returns w with length ilast, 1-based semantics reproduced on 0-based array.
    """
    path = Path(path)
    Nvtot = Nvel * Nvelb
    ilast = Norbit + Nvtot

    # stash by 1-based index like Fortran; slot 0 unused
    wt = np.zeros(ilast + 1, dtype=float)
    have = set()

    with path.open("r") as f:
        for line in f:
            s = line.strip().split()
            if len(s) < 2:
                continue
            ii = int(float(s[0]))  # tolerate "1  ..." or "1.0 ..."
            val = float(s[1])
            if 1 <= ii <= ilast:
                wt[ii] = val
                have.add(ii)

    w = np.zeros(ilast + 1, dtype=float)  # 1-based

    # full file path (no EOF in the Fortran loop): every i=1..ilast present
    if len(have) >= ilast and all((i in have) for i in range(1, ilast + 1)):
        w[1:ilast + 1] = wt[1:ilast + 1]
        return w[1:]  # drop dummy 0th

    # EOF path: library-sized weights were supplied
    # duplicate orbit weights into pairs
    for i in range(1, Norb + 1):
        w[2 * i - 1] = wt[i] / 2.0
        w[2 * i]     = wt[i] / 2.0

    # copy nuisance / velocity-bin terms exactly like the Fortran loop 1..Nvtot+1
    # yes, the +1 is preserved on purpose
    for i in range(1, Nvtot + 2):
        idx_src = Norb + i
        idx_dst = Norbit + i
        if idx_dst <= ilast:
            # if the src index wasn’t present in file, it stays 0.0 as in Fortran’s uninitialized wt default
            w[idx_dst] = wt[idx_src]

    return w[1:]

def weight_write(w, Norbit, Nvel, Nvelb, path=Path("weights.out")):
    ilast = int(Norbit + Nvel * Nvelb)
    w = np.asarray(w, dtype=float)
    if w.shape[0] == ilast + 1:      # 1-based style array: ignore w[0]
        vals = w[1:ilast+1]
    else:                            # 0-based style array
        vals = w[:ilast]
    idx = np.arange(1, ilast + 1, dtype=int)
    np.savetxt(path, np.column_stack([idx, vals]), fmt=["%d", "%.10g"])

def vaniout

def vaniread(path: Path | str,
             Norb: int,
             Norbit: int,
             Nrani: int,
             Ntani: int):
    """
    Read 'vani.out' written by the library stage.

    For each orbit i (1..Norb) the file stores, in order:
      v0liba, vr2liba, vt2liba, vrtliba, vp1liba, vp2liba
    each as Nrani*Ntani numbers, list-directed (arbitrary newlines).

    If Norbit != Norb (the model doubles the library), duplicate each orbit
    into columns (2*i-1, 2*i) with vp1 negated in the second copy.

    Returns: v0a, vr2a, vt2a, vp1a, vp2a with shape (Nrani, Ntani, Norbit).
    """
    path = Path(path)
    raw = np.fromfile(path, sep=' ', dtype=float)
    block = Nrani * Ntani
    per_orbit = 6 * block

    if raw.size % per_orbit != 0:
        raise ValueError(f"vani.out length {raw.size} not divisible by 6*Nrani*Ntani={per_orbit}")

    Norb_infer = raw.size // per_orbit
    if Norb is None:
        Norb = Norb_infer
    elif Norb_infer != Norb:
        # Be strict; Fortran would just read sequentially, so mismatch means the file doesn't match your expectations.
        raise ValueError(f"vani.out has {Norb_infer} orbits, but Norb={Norb}")

    # Allocate outputs with the model's column count
    v0a  = np.zeros((Nrani, Ntani, Norbit), dtype=float)
    vr2a = np.zeros_like(v0a)
    vt2a = np.zeros_like(v0a)
    vp1a = np.zeros_like(v0a)
    vp2a = np.zeros_like(v0a)
    # We read vrtliba but do not return it, matching the Fortran signature.
    # If you want it, add another return and fill like the others.

    for i in range(Norb):
        base = i * per_orbit
        # Pull the 6 planes for this orbit, reshaping in Fortran (column-major) order
        a0 = raw[base + 0*block : base + 1*block].reshape(Nrani, Ntani, order='F')  # v0liba
        a1 = raw[base + 1*block : base + 2*block].reshape(Nrani, Ntani, order='F')  # vr2liba
        a2 = raw[base + 2*block : base + 3*block].reshape(Nrani, Ntani, order='F')  # vt2liba
        a3 = raw[base + 3*block : base + 4*block].reshape(Nrani, Ntani, order='F')  # vrtliba (ignored)
        a4 = raw[base + 4*block : base + 5*block].reshape(Nrani, Ntani, order='F')  # vp1liba
        a5 = raw[base + 5*block : base + 6*block].reshape(Nrani, Ntani, order='F')  # vp2liba

        if Norbit == Norb:
            j = i
            v0a[:, :, j]  = a0
            vr2a[:, :, j] = a1
            vt2a[:, :, j] = a2
            vp1a[:, :, j] = a4
            vp2a[:, :, j] = a5
        else:
            j1, j2 = 2*i, 2*i + 1
            v0a[:, :, j1]  = a0
            v0a[:, :, j2]  = a0
            vr2a[:, :, j1] = a1
            vr2a[:, :, j2] = a1
            vt2a[:, :, j1] = a2
            vt2a[:, :, j2] = a2
            vp1a[:, :, j1] = a4
            vp1a[:, :, j2] = -a4   # sign flip matches the Fortran branch
            vp2a[:, :, j1] = a5
            vp2a[:, :, j2] = a5

    return v0a, vr2a, vt2a, vp1a, vp2a    


def xlib_read(path: Path, Norb: int, Norbit: int, Nbin: int, Nvlib: int, dtype=np.float64):
    """
    Read xlib.out and build smlib with shape (Nbin, Norbit).
    For each orbit i:
      - smlib[0, i] = sum(slush[0:Nvlib])     # smear lowest-r bin over all angles
      - for ibin=2..Nbin: smlib[ibin-1, i] = slush[Nvlib + (ibin-2)]
    If Norbit != Norb, duplicate each column i -> columns (2*i-1, 2*i) like the Fortran branch.
    """
    path = Path(path)
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if len(lines) < Norb:
        raise ValueError(f"xlib_read: expected ≥{Norb} lines, found {len(lines)}")

    # allocate target
    if Norbit != Norb:
        smlib = np.zeros((Nbin, 2 * Norb), dtype=dtype)
    else:
        smlib = np.zeros((Nbin, Norb), dtype=dtype)

    for iorb in range(Norb):                      # 0-based
        # Fortran read(71,*) slush  → space-separated numbers on one line
        sl = np.fromstring(lines[iorb], sep=' ', dtype=float)
        need = Nvlib + (Nbin - 1)                 # exactly what the Fortran indexing consumes
        if sl.size < need:
            raise ValueError(f"xlib_read: orbit {iorb+1} line has {sl.size} values, need at least {need}")

        # smear first coarse bin over theta
        first = sl[:Nvlib].sum()
        # copy remaining bins
        rest = sl[Nvlib : Nvlib + (Nbin - 1)]

        if Norbit != Norb:
            j1, j2 = 2*iorb, 2*iorb + 1          # columns 2*iorb-1 and 2*iorb in 1-based
            smlib[0,  j1] = first
            smlib[0,  j2] = first
            smlib[1:, j1] = rest
            smlib[1:, j2] = rest
        else:
            smlib[0,  iorb] = first
            smlib[1:, iorb] = rest

    return smlib