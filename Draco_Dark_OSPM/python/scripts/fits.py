import os
import logging
import numpy as np
import scipy as sc
import pandas as pd
from scipy.interpolate import interp1d
import astropy as ast
from astropy.io import fits
from astropy.io.votable import parse
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import match_coordinates_sky
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle as pckl
import os
from typing import Tuple, Optional, List
from dataclasses import dataclass
from astropy.units import pi as pi

# Fit stuff
######################################################################################################


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Port of the Fortran imfit-style workflow:
- read data from 'in' and limits from 'limits.imfit'
- read run settings from 'imfit.def' (with robust defaults)
- scale, window, and weight points
- fit polynomial of order iord-1 by weighted least squares 
- compute normalized series xo = y / y_fit, honoring the -666 sentinel
- optional region peak extraction from 'regions.dat' using a Tukey biweight
- optional smoothing of those peaks into a new trend line
- write 'out' with columns: w, xo
"""

# ------------------------------
# parameters mirroring Fortran
# ------------------------------
narrm = 400_000
npm   = 50
nl    = narrm

BIG = 1.0e30
TOL = 1.0e-5
ITERMAX = 100        # retained for parity, but we solve poly in closed form
FAC = 0.05           # scaling used on y = xd(:,1)/FAC
HCUT = 3.0e10        # only affected ymax in legacy; kept for completeness

# ------------------------------
# tiny utilities
# ------------------------------
def safe_read_limits(path='limits.imfit', default_ilo=-10000, default_ihi=10000) -> Tuple[int, int]:
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                parts = f.read().strip().split()
                if len(parts) >= 2:
                    return int(float(parts[0])), int(float(parts[1]))
        except Exception:
            pass
    return default_ilo, default_ihi

def parse_imfit_def(path='imfit.def'):
    # Soft-parse free-form keys; fall back to reasonable defaults
    # Fortran calls:
    #   qi1('Order of poly ', ... ) → iord
    #   qr1('Lower cutoff ', ... ) → cutl
    #   qi1('Number for peak values ', ... ) → nbi
    #   qi1('Linear (0) or Smooth (1) ', ... ) → itype
    defaults = dict(iord=3, cutl=-BIG, nbi=5, itype=0)
    if not os.path.exists(path):
        return defaults
    text = open(path, 'r', encoding='utf-8', errors='ignore').read().lower()
    def grab(keys, cast, fallback):
        for k in keys:
            if k in text:
                # take number following key on the same line
                for line in text.splitlines():
                    if k in line:
                        tokens = line.replace('=',' ').split()
                        # find the token right after the key
                        try:
                            idx = [t.lower() for t in tokens].index(k.split()[0])
                        except ValueError:
                            # brute-force scan for the first numeric in line
                            nums = [t for t in tokens if any(c.isdigit() for c in t)]
                            if nums:
                                try:
                                    return cast(nums[-1])
                                except Exception:
                                    break
                        else:
                            # scan right side for first numeric
                            for t in tokens[idx+1:]:
                                try:
                                    return cast(t)
                                except Exception:
                                    continue
        return fallback
    iord = grab(['order of poly'], int, defaults['iord'])
    cutl = grab(['lower cutoff'], float, defaults['cutl'])
    nbi  = grab(['number for peak values'], int, defaults['nbi'])
    itype= grab(['linear (0) or smooth (1)'], int, defaults['itype'])
    return dict(iord=iord, cutl=cutl, nbi=nbi, itype=itype)

def read_in_file(path='in') -> Tuple[np.ndarray, np.ndarray]:
    # Fortran read: each line has two numbers: w, xd(:,1)
    w = []
    xd1 = []
    with open(path, 'r') as f:
        for line in f:
            line=line.strip()
            if not line: 
                continue
            parts = line.replace(',', ' ').split()
            if len(parts) < 2:
                continue
            x1 = float(parts[0]); x2 = float(parts[1])
            w.append(x1); xd1.append(x2)
            if len(w) >= narrm:
                break
    if not w:
        raise RuntimeError("No data read from 'in'. Expect two columns.")
    return np.asarray(w, dtype=np.float64), np.asarray(xd1, dtype=np.float64)

def weighted_polyfit_xpow(x: np.ndarray, y: np.ndarray, sigma: np.ndarray, iord: int) -> np.ndarray:
    # Solve min || W^(1/2) (V a - y) || with V[:,j] = x**j (j=0..iord-1)
    # Matches Fortran funcs: y = sum_j a(j) * x**(j-1)
    V = np.vander(x, N=iord, increasing=True)     # columns: x^0, x^1, ...
    w = np.where(sigma > 0, 1.0 / sigma, 0.0)
    Wsqrt = np.sqrt(w)
    VW = V * Wsqrt[:, None]
    yW = y * Wsqrt
    a, *_ = np.linalg.lstsq(VW, yW, rcond=None)
    return a  # shape (iord,)

def polyval_xpow(a: np.ndarray, x: np.ndarray) -> np.ndarray:
    # a[0] + a[1] x + a[2] x^2 + ...
    V = np.vander(x, N=a.shape[0], increasing=True)
    return V @ a

def tukey_biweight_location(arr: np.ndarray, c: float = 6.0) -> float:
    # robust center; if degenerate, fall back to mean
    arr = np.asarray(arr, dtype=np.float64)
    if arr.size == 0:
        return np.nan
    M = np.median(arr)
    mad = np.median(np.abs(arr - M)) or 1e-12
    u = (arr - M) / (c * mad)
    mask = np.abs(u) < 1.0
    if not np.any(mask):
        return float(M)
    w = (1 - u[mask]**2)**2
    num = np.sum((arr[mask] - M) * w)
    den = np.sum(w) or 1e-12
    return float(M + num / den)

def linear_piecewise(xr: np.ndarray, yr: np.ndarray, xq: np.ndarray) -> np.ndarray:
    # emulate the Fortran linear segment walk with edge clamping
    if xr.size == 0:
        return np.zeros_like(xq)
    order = np.argsort(xr)
    xr2, yr2 = xr[order], yr[order]
    yq = np.interp(xq, xr2, yr2, left=yr2[0], right=yr2[-1])
    return yq

def smoothing_spline(xr: np.ndarray, yr: np.ndarray, xq: np.ndarray, s: Optional[float]) -> np.ndarray:
    # light dependency path: try scipy if available; else fall back to linear
    try:
        from scipy.interpolate import UnivariateSpline
        # s ~ smoothing factor; if None, let library choose; if 0, interpolate
        k = 3 if xr.size >= 4 else min(yr.size-1, 1)
        spl = UnivariateSpline(xr, yr, s=s if s is not None else None, k=k)
        return spl(xq)
    except Exception:
        return linear_piecewise(xr, yr, xq)

def read_regions(path='regions.dat') -> list[Tuple[int,int]]:
    regs = []
    if os.path.exists(path):
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                p = line.replace(',', ' ').split()
                if len(p) >= 2:
                    regs.append((int(round(float(p[0]))), int(round(float(p[1])))))
    return regs

def read_smflat_def(path='smflat.def') -> Optional[float]:
    if not os.path.exists(path):
        return None
    txt = open(path, 'r', encoding='utf-8', errors='ignore').read().lower()
    # look for a standalone number; Fortran asked: 'Enter smoothing val '
    nums = []
    for tok in txt.replace('=',' ').split():
        try:
            nums.append(float(tok))
        except Exception:
            pass
    return nums[-1] if nums else None

# ------------------------------
# main routine
# ------------------------------
def main():
    # stage 1: ingest
    ilo, ihi = safe_read_limits()
    w, xd1 = read_in_file('in')
    ncol = xd1.size
    x = np.arange(1, ncol+1, dtype=np.float64)  # Fortran x(i)=float(i)
    y = xd1 / FAC
    xo = y.copy()
    sig = np.ones_like(y, dtype=np.float64)

    # windowing and zero downweighting
    ymin, ymax = BIG, -BIG
    ilo = max(1, ilo)
    ihi = min(ncol, ihi)
    n0 = 0
    for i in range(ncol):
        if y[i] <= 0.1:
            n0 += 1
            sig[i] = 1000.0
        if (i+1) > ilo and (i+1) < ihi:
            ymin = min(ymin, y[i])
            ymax = max(ymax, y[i])
        else:
            sig[i] = 1000.0
    if n0 > ncol * 0.75:
        print(f"Too Many Zeros: {n0} {ncol}")
        # still proceed, as Fortran jumps to end but we can continue to write scaled output

    # stage 2: settings
    cfg = parse_imfit_def('imfit.def')
    iord = int(cfg['iord'])
    cutl = float(cfg['cutl'])
    nbi  = int(cfg['nbi'])
    itype= int(cfg['itype'])
    if iord < 1 or iord > npm:
        raise ValueError(f"iord={iord} out of range")

    # stage 3: polynomial fit (weighted LS; exact solution, no LM iterations needed)
    a = weighted_polyfit_xpow(x, y, sig, iord)
    yl = polyval_xpow(a, x)

    # stage 4: normalized residual-like quantity xo = y/yl with sentinels
    ylo = yl.copy()
    for i in range(ncol):
        if int(round(y[i])) == -666:
            xo[i] = -666.0
        else:
            # Fortran computed diff and ymax, but only used ymax as a diagnostic
            xo[i] = y[i] / yl[i] if yl[i] != 0 else 0.0

    # stage 5: optional region peak extraction and re-trend
    regs = read_regions('regions.dat')
    if regs:
        xr_list, yr_list = [], []
        for (i1, i2) in regs:
            i1 = max(1, i1); i2 = min(ncol, i2)
            # collect points in region that were not heavily downweighted
            idx = np.arange(i1-1, i2)
            good = sig[idx] < 1000.0
            xs = (idx[good] + 1).astype(np.float64)
            ys = xo[idx][good]               # note: Fortran used xo(j) here
            if xs.size == 0:
                xr_list.append(0.5*(i1+i2))
                yr_list.append(1.0)
                continue
            # sort by ys ascending, then reverse to descending to emulate sort2 + reverse
            order = np.argsort(ys)
            xs_sorted = xs[order][::-1]
            ys_sorted = ys[order][::-1]
            m = int(min(xs_sorted.size, max(1, nbi)))
            xr_list.append(tukey_biweight_location(xs_sorted[:m]))
            yr_list.append(tukey_biweight_location(ys_sorted[:m]))
        xr = np.asarray(xr_list, dtype=np.float64)
        yr = np.asarray(yr_list, dtype=np.float64)

        # produce a new trend yl either smooth or linear piecewise across xr,yr
        if itype == 1:
            sval = read_smflat_def('smflat.def')
            yl2 = smoothing_spline(xr, yr, x, s=sval)
        else:
            yl2 = linear_piecewise(xr, yr, x)

        # redraw normalization against the new trend and compute xo2 like Fortran
        xo2 = np.zeros_like(x, dtype=np.float64)
        for i in range(ncol):
            if int(round(xo[i])) == -666:
                continue
            if yl2[i] > 0:
                xo[i] = xo[i] / yl2[i]
                xo2[i] = yl2[i] * ylo[i]
            else:
                xo[i] = 0.0
                xo2[i] = 0.0

    # stage 6: write result
    with open('out', 'w') as fout:
        for i in range(ncol):
            fout.write(f"{w[i]:.10g} {xo[i]:.10g}\n")

    # parity print (optional)
    print(f"ncol={ncol} iord={iord} cutl={cutl} nbi={nbi} itype={itype}  ilo={ilo} ihi={ihi}")
    print("Wrote 'out' with columns: w, xo")

if __name__ == '__main__':
    main()



