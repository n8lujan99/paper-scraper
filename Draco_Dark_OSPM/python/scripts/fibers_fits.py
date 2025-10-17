#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from typing import Tuple, Optional, List
from dataclasses import dataclass
import math
import argparse

FILE_DRACO = "Draco_stars.csv"
FLUX = "flux"
MAG  = "phot_g_mean_mag"
RA, RA_error, Dec, Dec_error = "RA_deg", "RA_error", "Dec_deg", "Dec_error"
	
Parallax, Parallax_error, PMRA, PMRA_error = "parallax", "parallax_error", "pmra", "pmra_error"
PMDEC, PMDEC_error = "pmdec", "pmdec_error"
Radial_Velocity, Radial_Velocity_error = "radial_velocity", "Radial_Velocity_error"	
L, B, Phot_g_mean_mag, Phot_bp_mean_mag, Phot_rp_mean_mag, bp_rp, ruwe, distance_gspphot, phot_bp_rp_excess_factor = "L","B","phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp", "ruwe", "distance_gspphot", "phot_bp_rp_excess_factor"



"""
This module includes 3 workflows:

1) Estimate position of gravitational center of stellar data then compare to optical center of Draco:
    - Read in Draco_stars utilizing: load_df()
2) imfit: 
   - Read two-column 'in' and limits from 'limits.imfit'
   - Weighted poly fit (order iord-1), robust peak selection from 'regions.dat'
   - Output 'out' with w, xo

3) fibers:
   - Read multi-column 'in' with fiber meta (xr0, xd0, xf, xw, names, az, ferr)
   - Scan amplitude and grid search center using Gaussian or Moffat PSF over fiber apertures
   - Atmospheric differential refraction correction table via adcor()
   - Outputs: 'out2', 'out', 'fac2.out'
"""

# ============================================================================================
# Shared constants
# ============================================================================================
BIG                    = 1.0e30
TOL                    = 1.0e-5
FAC                    = 0.05
HCUT                   = 3.0e10
PI                     = 3.141593
# ============================================================================================
# File names and paths
# ============================================================================================
FILE_Draco             = "Draco_stars"
FILE_IN                = "in"
FILE_LIMITS            = "limits.imfit"
FILE_IMFIT_DEF         = "imfit.def"
FILE_REGIONS           = "regions.dat"
FILE_SMFLAT            = "smflat.def"
FILE_OUT               = "out"

# ============================================================================================
# Center specific
# ============================================================================================

# ============================================================================================
# Fibers specific
# ============================================================================================
FILE_FWHM_USE          = "fwhm.use"
FILE_FWHM_FIX          = "fwhm.fix"
FILE_NORMEXP_OUT       = "normexp.out"
FILE_OUT2              = "out2"
FILE_FAC2              = "fac2.out"
# ============================================================================================
# ---------- IMFIT WORKFLOW (polynomial normalization) ----------
# ============================================================================================

narrm = 400_000
npm   = 50
nl    = narrm




def safe_read_limits(path='FILE_LIMITS', default_ilo=-10000, default_ihi=10000) -> Tuple[int, int]:
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                parts = f.read().strip().split()
                if len(parts) >= 2:
                    return int(float(parts[0])), int(float(parts[1]))
        except Exception:
            pass
    return default_ilo, default_ihi

def parse_imfit_def(path='FILE_IMFIT_DEF'):
    defaults = dict(iord=3, cutl=-BIG, nbi=5, itype=0)
    if not os.path.exists(path):
        return defaults
    text = open(path, 'r', encoding='utf-8', errors='ignore').read().lower()
    def grab(keys, cast, fallback):
        for k in keys:
            if k in text:
                for line in text.splitlines():
                    if k in line:
                        tokens = line.replace('=',' ').split()
                        try:
                            idx = [t.lower() for t in tokens].index(k.split()[0])
                        except ValueError:
                            nums = [t for t in tokens if any(c.isdigit() for c in t)]
                            if nums:
                                try:
                                    return cast(nums[-1])
                                except Exception:
                                    break
                        else:
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

def read_in_file_simple(path='FILE_IN') -> Tuple[np.ndarray, np.ndarray]:
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
    V = np.vander(x, N=iord, increasing=True)
    w = np.where(sigma > 0, 1.0 / sigma, 0.0)
    Wsqrt = np.sqrt(w)
    VW = V * Wsqrt[:, None]
    yW = y * Wsqrt
    a, *_ = np.linalg.lstsq(VW, yW, rcond=None)
    return a

def polyval_xpow(a: np.ndarray, x: np.ndarray) -> np.ndarray:
    V = np.vander(x, N=a.shape[0], increasing=True)
    return V @ a

def tukey_biweight_location(arr: np.ndarray, c: float = 6.0) -> float:
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
    if xr.size == 0:
        return np.zeros_like(xq)
    order = np.argsort(xr)
    xr2, yr2 = xr[order], yr[order]
    yq = np.interp(xq, xr2, yr2, left=yr2[0], right=yr2[-1])
    return yq

def smoothing_spline(xr: np.ndarray, yr: np.ndarray, xq: np.ndarray, s: Optional[float]) -> np.ndarray:
    try:
        from scipy.interpolate import UnivariateSpline
        k = 3 if xr.size >= 4 else min(yr.size-1, 1)
        spl = UnivariateSpline(xr, yr, s=s if s is not None else None, k=k)
        return spl(xq)
    except Exception:
        return linear_piecewise(xr, yr, xq)

def read_regions(path='FILE_REGIONS') -> List[Tuple[int,int]]:
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

def read_smflat_def(path='FILE_SMFLAT') -> Optional[float]:
    if not os.path.exists(path):
        return None
    txt = open(path, 'r', encoding='utf-8', errors='ignore').read().lower()
    nums = []
    for tok in txt.replace('=',' ').split():
        try:
            nums.append(float(tok))
        except Exception:
            pass
    return nums[-1] if nums else None

def run_imfit():
    ilo, ihi = safe_read_limits()
    w, xd1 = read_in_file_simple('in')
    ncol = xd1.size
    x = np.arange(1, ncol+1, dtype=np.float64)
    y = xd1 / FAC
    xo = y.copy()
    sig = np.ones_like(y, dtype=np.float64)

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

    cfg = parse_imfit_def('FILE_IMFIT')
    iord = int(cfg['iord'])
    cutl = float(cfg['cutl'])
    nbi  = int(cfg['nbi'])
    itype= int(cfg['itype'])
    if iord < 1 or iord > npm:
        raise ValueError(f"iord={iord} out of range")

    a = weighted_polyfit_xpow(x, y, sig, iord)
    yl = polyval_xpow(a, x)

    ylo = yl.copy()
    for i in range(ncol):
        if int(round(y[i])) == -666:
            xo[i] = -666.0
        else:
            xo[i] = y[i] / yl[i] if yl[i] != 0 else 0.0

    regs = read_regions('FILE_REGIONS')
    if regs:
        xr_list, yr_list = [], []
        for (i1, i2) in regs:
            i1 = max(1, i1); i2 = min(ncol, i2)
            idx = np.arange(i1-1, i2)
            good = sig[idx] < 1000.0
            xs = (idx[good] + 1).astype(np.float64)
            ys = xo[idx][good]
            if xs.size == 0:
                xr_list.append(0.5*(i1+i2))
                yr_list.append(1.0)
                continue
            order = np.argsort(ys)
            xs_sorted = xs[order][::-1]
            ys_sorted = ys[order][::-1]
            m = int(min(xs_sorted.size, max(1, nbi)))
            xr_list.append(tukey_biweight_location(xs_sorted[:m]))
            yr_list.append(tukey_biweight_location(ys_sorted[:m]))
        xr = np.asarray(xr_list, dtype=np.float64)
        yr = np.asarray(yr_list, dtype=np.float64)
        if itype == 1:
            sval = read_smflat_def('FILE_SMFLAT')
            yl2 = smoothing_spline(xr, yr, x, s=sval)
        else:
            yl2 = linear_piecewise(xr, yr, x)
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

    with open('out', 'w') as fout:
        for i in range(ncol):
            fout.write(f"{w[i]:.10g} {xo[i]:.10g}\n")

    print(f"ncol={ncol} iord={iord} cutl={cutl} nbi={nbi} itype={itype}  ilo={ilo} ihi={ihi}")
    print("Wrote 'out' with columns: w, xo")

# ============================================================================================
# ---------- FIBER WORKFLOW (Gaussian/Moffat over fiber apertures) ----------
# ============================================================================================

@dataclass
class CSigma:
    rsig: float = 0.0   # Gaussian sigma (arcsec)
    fmof: float = 0.0   # Moffat FWHM proxy scale
    bmof: float = 3.5   # Moffat beta
    imoff: int  = 1     # 0→Gaussian, 1→Moffat

cs = CSigma()

def read_fwhm() -> float:
    rfw = 1.55
    try:
        with open('FILE_FWHM_USE','r') as f:
            rfw = float(f.read().split()[0])
    except Exception:
        pass
    try:
        with open('FILE_FWHM_FIX','r') as f:
            rfw0 = float(f.read().split()[0])
        if rfw0 < 0:
            rfw = -rfw0
    except Exception:
        pass
    return rfw

def read_normexp() -> Tuple[int, List[str], np.ndarray]:
    names = ['exp01','exp02','exp03']
    rel = np.array([1.0,1.0,1.0], dtype=float)
    ne = 3
    try:
        with open('FILE_NORMEXP_OUT','r') as f:
            names, vals = [], []
            for _ in range(3):
                line = f.readline()
                if not line:
                    break
                parts = line.split()
                if len(parts) >= 3:
                    names.append(parts[0])
                    vals.append(float(parts[2]))
            if names:
                ne = len(names)
                rel = np.array(vals, dtype=float)
    except Exception:
        pass
    return ne, names, rel

def read_in_table_fibers(nmax=3000):
    xr0=[]; xd0=[]; xf=[]; xw=[]; az=[]; ferr=[]
    an1=[]; an2=[]; an3=[]; an4=[]
    with open('in','r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            x1,x2,x3,x4 = map(float, parts[:4])
            a5,a6,a7,a8 = parts[4:8]
            x9,x10 = map(float, parts[8:10])
            xr0.append(x1); xd0.append(x2); xf.append(x3); xw.append(x4)
            if x3 == 0.0: xw[-1] = 0.0
            an1.append(a5); an2.append(a6); an3.append(a7); an4.append(a8)
            az.append(x9); ferr.append(x10)
            if len(xr0) >= nmax:
                break
    if not xr0:
        return 0, None
    n = len(xr0)
    out = dict(
        xr0=np.array(xr0,dtype=float),
        xd0=np.array(xd0,dtype=float),
        xf=np.array(xf,dtype=float),
        xw=np.array(xw,dtype=float),
        an=np.stack([np.array(an1, dtype=object), np.array(an2, dtype=object),
                    np.array(an3, dtype=object), np.array(an4, dtype=object)], axis=1),
        az=np.array(az,dtype=float),
        ferr=np.array(ferr,dtype=float),
    )
    return n, out

def getchifib(xrs, xds, amps, n, xr, xd, xf, xw, fe, xrel, iflag, an, ip=0):
    rfib = 0.75
    nstep = 30
    xstep = 2.0*rfib/float(nstep-1)
    deltx = PI*rfib*rfib
    area_gaus  = amps * deltx / (2.0*cs.rsig*cs.rsig*PI) if cs.rsig != 0 else 0.0
    am = 4.0*(2.0**(1.0/cs.bmof)-1.0)*(cs.bmof-1.0)/(PI*cs.fmof*cs.fmof)
    area_moff  = amps * deltx * am

    gausa = np.zeros(n, dtype=float)
    da    = np.zeros(n, dtype=float)
    chi   = 0.0

    if ip == 1:
        print("Ifib     Counts      Fit        Distance        C-F           Chi")

    for i in range(n):
        xs = xr[i] - rfib
        ys = xd[i] - rfib
        gaus = 0.0
        xmoff= 0.0
        nsum = 0
        for ix in range(nstep):
            xp = xs + xstep*ix
            for iy in range(nstep):
                yp = ys + xstep*iy
                dist0 = math.hypot(xp - xr[i], yp - xd[i])
                if dist0 < rfib:
                    dist = math.hypot(xp - xrs, yp - xds)
                    g = dist / cs.rsig if cs.rsig != 0 else 0.0
                    gaus  += math.exp(-0.5*g*g) * area_gaus
                    xmoff += area_moff * (1.0 + 4.0*(2.0**(1.0/cs.bmof)-1.0)*(dist/cs.fmof)**2)**(-cs.bmof)
                    nsum  += 1
        if nsum > 0:
            gaus  /= nsum
            xmoff /= nsum
        fit = xmoff if cs.imoff == 1 else gaus
        gausa[i] = fit
        chi1 = 0.0
        if fe[i] > 0:
            chi1 = xw[i]*(fit - xf[i])**2 / (fe[i]**2)
        da[i] = xf[i] - fit
        if iflag[i] == 0:
            chi += chi1
        if ip == 1:
            rad = math.hypot(xr[i]-xrs, xd[i]-xds)
            print(f"{i+1:3d} {xf[i]:12.2f}  {fit:9.2f}  {rad:9.2f}  {da[i]:12.2f}  {chi1:12.2f}  "
                f"{an[i,0]:>17s} {an[i,1]:>8s} {an[i,2]:>3s} {an[i,3]:>5s} {iflag[i]:1d}")

    # coverage fraction
    sumf1 = 0.0
    sumf2 = 0.0
    nfull = 100
    sigfull = 5.0
    xs = xrs - sigfull*cs.rsig
    xe = xrs + sigfull*cs.rsig
    ys = xds - sigfull*cs.rsig
    ye = xds + sigfull*cs.rsig
    for ix in range(nfull):
        xp = xs + (xe-xs)*ix/float(nfull-1)
        for jx in range(nfull):
            yp = ys + (ye-ys)*jx/float(nfull-1)
            covered = False
            for i in range(n):
                if math.hypot(xp - xr[i], yp - xd[i]) < 0.75:
                    covered = True
                    break
            if covered:
                sumf1 += 1.0
            sumf2 += 1.0
    sumrat = (sumf1/sumf2) if sumf2 > 0 else 0.0
    return chi, da, gausa, sumrat

def adcor(nw, wadc, adc, n, xr, xd, az):
    fadc = np.zeros((n, 10), dtype=float)
    rfib = 0.75
    nstep = 50
    xstep = 2.0*rfib/float(nstep-1)
    deltx = PI*rfib*rfib
    area_gaus  = 1.0 * deltx / (2.0*cs.rsig*cs.rsig*PI) if cs.rsig != 0 else 0.0
    am = 4.0*(2.0**(1.0/cs.bmof)-1.0)*(cs.bmof-1.0)/(PI*cs.fmof*cs.fmof)
    area_moff  = deltx * am
    dtr = 180.0 / PI

    for ia in range(nw):
        for i in range(n):
            xaoff = adc[ia] * math.sin(az[i]/dtr)
            yaoff = adc[ia] * math.cos(az[i]/dtr)
            xs = xr[i] - rfib + xaoff
            ys = xd[i] - rfib + yaoff
            gaus = 0.0
            xmoff= 0.0
            nsum = 0
            for ix in range(nstep):
                xp = xs + xstep*ix
                for iy in range(nstep):
                    yp = ys + xstep*iy
                    if math.hypot(xp - xr[i], yp - xd[i]) < rfib:
                        dist = math.hypot(xp - 0.0, yp - 0.0)
                        g = dist / cs.rsig if cs.rsig != 0 else 0.0
                        gaus  += math.exp(-0.5*g*g) * area_gaus
                        xmoff += area_moff * (1.0 + 4.0*(2.0**(1.0/cs.bmof)-1.0)*(dist/cs.fmof)**2)**(-cs.bmof)
                        nsum  += 1
            if nsum > 0:
                gaus  /= nsum
                xmoff /= nsum
            fadc[i, ia] = xmoff if cs.imoff == 1 else gaus
    return fadc

def run_fibers(xrs0: float, xds0: float):
    cs.imoff = 1
    rfw = read_fwhm()
    print(f"Using FWHM: {rfw}")
    cs.rsig = rfw / 2.35
    cs.fmof = rfw
    cs.bmof = 3.5

    ne, aname, relnorm = read_normexp()
    n, dat = read_in_table_fibers(nmax=3000)
    if n == 0:
        print("No rows in 'in'.")
        return

    xr0 = dat['xr0']; xd0 = dat['xd0']; xf = dat['xf']; xw = dat['xw']
    an  = dat['an'];   az  = dat['az'];  ferr = dat['ferr']

    xrel = np.ones(n, dtype=float)
    for i in range(n):
        tag = str(an[i,3])[:5]
        for j, nm in enumerate(aname[:ne]):
            if tag == nm:
                xrel[i] = relnorm[j]
                break

    iflag = np.zeros(n, dtype=int)
    if np.any(xf == 0.0):
        iflag[xf == 0.0] = 1
        xw = xw.copy()
        xw[xf == 0.0] = 0.0

    xr = (xr0 - xrs0) * 3600.0 * math.cos(xds0/57.3)
    xd = (xd0 - xds0) * 3600.0
    xfmax = float(np.max(xf)) if n > 0 else 0.0

    as_ = 40.0
    ae  = max(2000.0, xfmax*5.0)
    if ae < as_:
        print("Problem with negative flux")
        return
    chimin = 1.0e10
    atb = as_
    for ia in range(1, 101):
        at = as_ + (ae - as_)*(ia-1)/float(100-1)
        chi, _, _, _ = getchifib(0.0, 0.0, at, n, xr, xd, xf, xw, ferr, xrel, iflag.copy(), an, ip=0)
        if chi < chimin:
            chimin = chi
            atb = at
    amps = atb
    print("Begin:")
    chi, da, gausa, sumrat = getchifib(0.0, 0.0, amps, n, xr, xd, xf, xw, ferr, xrel, iflag.copy(), an, ip=1)

    nstepa = 30
    nstepc = 21
    xs, xe = -0.2, 0.2
    ys, ye = -0.2, 0.2
    as2, ae2 = 0.7*amps, 1.4*amps
    chimin = 1.0e10
    xtb = 0.0; ytb = 0.0; atb = amps
    for i in range(nstepc):
        xt = xs + (xe - xs)*i/float(nstepc-1)
        for j in range(nstepc):
            yt = ys + (ye - ys)*j/float(nstepc-1)
            for ia in range(nstepa):
                at = as2 + (ae2 - as2)*ia/float(nstepa-1)
                chi, _, _, _ = getchifib(xt, yt, at, n, xr, xd, xf, xw, ferr, xrel, iflag.copy(), an, ip=0)
                if chi < chimin:
                    chimin = chi; xtb = xt; ytb = yt; atb = at

    print("End:")
    chi, da, gausa, sumrat = getchifib(xtb, ytb, atb, n, xr, xd, xf, xw, ferr, xrel, iflag.copy(), an, ip=1)
    print(xtb, ytb, xrs0, xds0, atb/amps)

    xnew = xrs0 + xtb/3600.0/math.cos(xds0/57.3)
    ynew = xds0 + ytb/3600.0
    print("RAnew, DECnew, Amp, chimin, FW")
    print(xnew, ynew, atb, chimin, rfw)

    sumg = float(np.sum(gausa))
    sumd = float(np.sum(xf))
    sumrat = sumg/atb if atb != 0 else 0.0

    with open('out2','w') as f:
        f.write(f"{xnew:9.5f} {ynew:9.5f} {atb:10.2f} {xtb:8.3f} {ytb:8.3f} {chimin:10.2f} {rfw:6.3f} {sumrat:6.3f}\n")

    wadc = np.array([3500., 4000., 4500., 5000., 5500.], dtype=float)
    adc  = np.array([-0.74, -0.40, 0.0, 0.08, 0.20], dtype=float)
    nw   = 5
    fadc = adcor(nw, wadc, adc, n, xr, xd, az)

    with open('FILE_OUT2','w') as f:
        for i in range(n):
            fadcw = [fadc[i,j]/fadc[i,2] if fadc[i,2] != 0 else 0.0 for j in range(nw)]
            f.write(f"  {xr0[i]:10.5f}  {xd0[i]:10.5f} {xf[i]:10.2f} {xw[i]:5.3f}  "
                    f"{an[i,0]:>17s} {an[i,1]:>8s} {an[i,2]:>3s} {an[i,3]:>5s} "
                    f" {iflag[i]:1d} {gausa[i]:11.2f} "
                    f"{fadcw[0]:5.3f} {fadcw[1]:5.3f} {fadcw[2]:5.3f} {fadcw[3]:5.3f} {fadcw[4]:5.3f}\n")

    with open('FILE_FAC2','w') as f:
        scale = atb/sumd if sumd != 0 else 0.0
        f.write(f"{sumrat} {scale}\n")

    print("Fiber coverage:", sumrat, atb/(sumd if sumd != 0 else 1.0))
# ============================================================================================
# Entry point with mode switch
# ============================================================================================

def main():
    parser = argparse.ArgumentParser(description="imfit + fiber workflows")
    parser.add_argument("--mode", choices=["imfit","fibers"], default="imfit", help="Choose the workflow to run.")
    parser.add_argument("--xrs0", type=float, default=0.0, help="Reference RA (deg) for fibers mode.")
    parser.add_argument("--xds0", type=float, default=0.0, help="Reference DEC (deg) for fibers mode.")
    args = parser.parse_args()

    if args.mode == "imfit":
        run_imfit()
    else:
        run_fibers(args.xrs0, args.xds0)

if __name__ == '__main__':
    main()
