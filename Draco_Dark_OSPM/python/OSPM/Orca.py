# 🧠 MAIN ORCHESTRATOR
# Runs build_library(), librarian(), etc.
# continue with core comps orbit ints and whats left in ospm_mods
###################################################################################
# Imports
###################################################################################
import os
import sys
from OSPM.runtime_state import model_params
from OSPM.OSPM_Binning import *
from OSPM.orbit_readers_writers import *

import numpy as np 
import scipy as sp
import pandas as pd 
from scipy.stats import norm
from scipy import constants as phys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.linalg import lu_factor, lu_solve

###################################################################################
M_sun = 1.98847e30                    # Solar mass [kg]
pc = phys.parsec                      # 3.085677581e16 m
kpc = 1e3 * pc
arcsec = phys.arcsec                  # radians per arcsecond (4.848e-6)
pi = np.pi
G_SI = phys.G                         # 6.67430e-11 m^3 / (kg s^2)

GG = G_SI * (1e-3)**2 * (kpc) / M_sun   # 4.30091e-6 kpc (km/s)^2 / Msun

###################################################################################

###################################################################################

# Space grid setup
def spacel():
    pass

def sos2nuclei(Norbit, Nsos=55, nquadlength=10000, dist_lim=0.1, fracshift=0.001, minshift=1e-7, maxshift=1e5, tooclose=1e-8, do_plot_elz=False):
    def read_sos():
        """Read sos1.out and sos2.out -> rs, vs, ts arrays + bins."""
        mlnorm = float(np.loadtxt("gden.norm"))
        df2 = pd.read_csv("sos2.out", delim_whitespace=True, names=["iorb", "iE", "iLz", "dLz", "dE"])
        df1 = pd.read_csv("sos1.out", delim_whitespace=True, names=["rs_in", "vs_in", "ts_in"])

        fac1 = 1.0  # external scaling not applied here
        rmax = 1.0

        rs = np.full((Norbit, Nsos), np.nan)
        vs = np.full_like(rs, np.nan)
        ts = np.full_like(rs, np.nan)
        E_bin = np.zeros(Norbit, int)
        Lz_bin = np.zeros(Norbit, int)
        nsosorb = np.zeros(Norbit, int)

        for i, (iorb, iE, iLz, dLz, dE) in enumerate(df2.values):
            rs_in, vs_in, ts_in = df1.iloc[i]
            iorb = int(iorb) - 1  # Fortran 1-based -> Python 0-based
            nsosorb[iorb] += 1
            idx = nsosorb[iorb] - 1
            rs[iorb, idx] = np.log10(rs_in * rmax)
            vs[iorb, idx] = vs_in * fac1 * mlnorm
            ts[iorb, idx] = ts_in
            E_bin[iorb] = int(iE)
            Lz_bin[iorb] = int(iLz)

        return rs, vs, ts, nsosorb, E_bin, Lz_bin

    # === helper: group orbits into unique (E,Lz) sequences =====================
    def sosseq(E_bin, Lz_bin):
        data = pd.DataFrame({"E": E_bin, "Lz": Lz_bin})
        grouped = data.drop_duplicates().reset_index(drop=True)
        E_seq = grouped["E"].values
        Lz_seq = grouped["Lz"].values
        nseq = len(E_seq)
        seqoforb = np.zeros(len(E_bin), int)
        for iorb, (E, Lz) in enumerate(zip(E_bin, Lz_bin)):
            seqoforb[iorb] = np.where((E_seq == E) & (Lz_seq == Lz))[0][0]
        return nseq, E_seq, Lz_seq, seqoforb

    # === helper: plot ELz like Fortran plotELz (optional) =====================
    def _plot_ELz():
        """Replicates plotELz: scatter of log10(-E) vs ±Lz from integrals.out."""
        # integrals.out lines: xLz, xE, x3, x4; Fortran used E=log10(-x2), Lz=±x1
        arr = np.loadtxt("integrals.out")
        xLz = arr[:, 0]
        xE  = arr[:, 1]
        E_vals = np.log10(-xE)
        Lz_vals = xLz
        # duplicate for symmetric ±Lz, as in Fortran
        E_plot = np.repeat(E_vals, 2)
        Lz_plot = np.empty(2 * len(Lz_vals))
        Lz_plot[0::2] = -Lz_vals
        Lz_plot[1::2] =  Lz_vals

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(E_plot, Lz_plot, s=10)
        ax.set_xlabel(r'$\log_{10}(-E)$')
        ax.set_ylabel(r'$L_z$')
        ax.set_title('ELz diagram (integrals.out)')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        return fig, ax

    # ===================== main routine =======================================
    rs, vs, ts, nsosorb, E_bin, Lz_bin = read_sos()
    nseq, E_seq, Lz_seq, seqoforb = sosseq(E_bin, Lz_bin)

    nuclei_list, vor_cells_list, vor_sos_list, vor_plot_list = [], [], [], []
    seq_summary = []

    for iseq in range(nseq):
        # collect all SOS points for this sequence
        idxs = np.where(seqoforb == iseq)[0]
        x, y, tarr, orbarr = [], [], [], []
        for iorb in idxs:
            count = int(nsosorb[iorb])
            if count <= 0:
                continue
            x.extend(10 ** rs[iorb, :count])
            y.extend(vs[iorb, :count])
            tarr.extend(ts[iorb, :count])
            orbarr.extend([iorb + 1] * count)

        if not x:
            continue

        x = np.asarray(x); y = np.asarray(y)
        tarr = np.asarray(tarr); orbarr = np.asarray(orbarr)

        # normalize x,y by their maxima
        xmax, ymax = np.max(x), np.max(y)
        x = x / xmax
        y = y / ymax

        # --- build envelope (upper ridge), then mirror across v=0
        i_max = np.argmax(y)
        xt = [x[i_max]]
        yt = [y[i_max] + np.clip(fracshift * y[i_max], minshift, maxshift)]

        # to larger radii
        xcur = x[i_max]
        while True:
            idx = np.where(x > xcur)[0]
            if not len(idx):
                break
            i_next = idx[np.argmax(y[idx])]
            if i_next == i_max:
                break
            xcur = x[i_next]
            xt.append(xcur + np.clip(fracshift * xcur, minshift, maxshift))
            yt.append(y[i_next] + np.clip(fracshift * y[i_next], minshift, maxshift))
            if np.isclose(xcur, x.max()):
                break

        # to smaller radii
        xcur = x[i_max]
        while True:
            idx = np.where(x < xcur)[0]
            if not len(idx):
                break
            i_next = idx[np.argmax(y[idx])]
            if i_next == i_max:
                break
            xcur = x[i_next]
            xt.append(xcur * (1 - fracshift))
            yt.append(y[i_next] * (1 + fracshift))
            if np.isclose(xcur, x.min()):
                break

        xt = np.asarray(xt); yt = np.asarray(yt)
        order = np.argsort(xt)
        xt, yt = xt[order], yt[order]

        # interpolate to densify boundary (Fortran's ysmall logic)
        ysmall = max(min((xt[-1] - xt[0]) / 40.0, 1e-2), 1e-3)
        xti, yti = [xt[0]], [yt[0]]
        for i in range(1, len(xt)):
            dx = xt[i] - xt[i - 1]; dy = yt[i] - yt[i - 1]
            steps = int(max(abs(dx), abs(dy)) / ysmall)
            for j in range(1, steps):
                f = j / steps
                xti.append(xt[i - 1] + f * dx)
                yti.append(yt[i - 1] + f * dy)
            xti.append(xt[i]); yti.append(yt[i])
        xt = np.asarray(xti); yt = np.asarray(yti)

        # mirror points with y<dist_lim
        xenv = np.concatenate([1.001 * x[y < dist_lim], xt])
        yenv = np.concatenate([-y[y < dist_lim],       yt])

        # back to physical scale
        x_all = xenv * xmax
        y_all = yenv * ymax
        orb_all = np.concatenate([np.zeros_like(xenv, int), np.ones_like(xt, int)])
        t_all = np.zeros_like(x_all)

        # drop pairs that are too close (scaled by nquadlength like Fortran)
        valid = np.ones(len(x_all), dtype=bool)
        for i in range(len(x_all)):
            if not valid[i]:
                continue
            dx = np.abs(x_all[i] - x_all[i+1:]) / float(nquadlength)
            dy = np.abs(y_all[i] - y_all[i+1:]) / float(nquadlength)
            bad = (dx < tooclose) & (dy < tooclose)
            valid[i+1:][bad] = False

        x_all, y_all, orb_all, t_all = x_all[valid], y_all[valid], orb_all[valid], t_all[valid]

        # assemble outputs analogous to files
        N = len(x_all)
        xjunk = np.full(N, 0.30)
        indices = np.arange(N)
        nuclei   = np.column_stack((x_all, y_all, xjunk, indices))
        vor_cells= np.tile([1, 1, 1, -1.0], (N, 1))
        vor_sos  = np.column_stack((x_all, y_all, xjunk, orb_all, t_all, np.full(N, xmax), np.full(N, ymax)))
        vor_plot = np.column_stack((x_all, y_all, orb_all, np.zeros(N)))

        nuclei_list.append(nuclei)
        vor_cells_list.append(vor_cells)
        vor_sos_list.append(vor_sos)
        vor_plot_list.append(vor_plot)
        seq_summary.append((iseq + 1, E_seq[iseq], Lz_seq[iseq], N))

    seq_df = pd.DataFrame(seq_summary, columns=["seq", "E_seq", "Lz_seq", "N_points"])

    out = {
        "nseq": nseq,
        "nuclei": nuclei_list,
        "vor_cells": vor_cells_list,
        "vor_sos": vor_sos_list,
        "vor_plot": vor_plot_list,
        "seq_info": seq_df,
    }

    if do_plot_elz:
        fig, ax = _plot_ELz()
        out["elz_fig"] = fig
        out["elz_ax"]  = ax

    return out

def velmtod(veld, velm, ad, adfer, xmu, Nvelb):
    """Finds which velocity bins of data correspond to model velocity bin edges."""
    nad = veld.shape[0]
    Nvel = velm.shape[0]

    sumad = np.zeros((Nvel, Nvelb))
    sadfer = np.zeros((Nvel, Nvelb))

    vbin = velm[1] - velm[0]

    # Loop over apertures
    for irc in range(Nvelb):
        # --- Compute data bin edges
        v1 = np.zeros(nad)
        v2 = np.zeros(nad)

        v1[0] = veld[0, irc] * xmu - (veld[1, irc] - veld[0, irc]) * xmu / 2.0
        v2[0] = veld[0, irc] * xmu + (veld[1, irc] - veld[0, irc]) * xmu / 2.0
        v1[-1] = veld[-1, irc] * xmu - (veld[-1, irc] - veld[-2, irc]) * xmu / 2.0
        v2[-1] = veld[-1, irc] * xmu + (veld[-1, irc] - veld[-2, irc]) * xmu / 2.0

        for k in range(1, nad - 1):
            v1[k] = veld[k, irc] * xmu - (veld[k, irc] - veld[k - 1, irc]) * xmu / 2.0
            v2[k] = veld[k, irc] * xmu + (veld[k + 1, irc] - veld[k, irc]) * xmu / 2.0

        sum1 = np.sum(ad[:, irc])  # total amplitude

        # --- Map data bins into model bins
        for i in range(Nvel):
            vpos1 = velm[i] - vbin / 2.0
            vpos2 = velm[i] + vbin / 2.0
            total = 0.0
            totale = 0.0

            for k in range(nad):
                frac = 0.0
                if vpos1 < v1[k] and vpos2 > v2[k]:
                    frac = 1.0
                elif vpos1 > v1[k] and vpos1 < v2[k]:
                    frac = (v2[k] - vpos1) / (v2[k] - v1[k])
                elif vpos2 > v1[k] and vpos2 < v2[k]:
                    frac = (vpos2 - v1[k]) / (v2[k] - v1[k])

                total += frac * ad[k, irc]
                totale += frac * adfer[k, irc]

            if velm[i] < v1[0] or velm[i] > v2[-1]:
                total = 0.0
                totale = -666.0

            sumad[i, irc] = total
            sadfer[i, irc] = totale

    return sumad, sadfer

###################################################################################
# Assemble Matrices
###################################################################################

def cmatrix (xlib, v1lib, Norbit, Norbitm, Nvel, Nvelb, Nvelbm, Nbin):
    # Assembles the matrix in Cm for the matrix Am in spear
    """Assembles the matrix in Cm for the matrix Am in spear."""
    lds = Norbitm + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    Cm = np.zeros((lda, lds), dtype=float)
    Nvtot = Nvel * Nvelb

    # Left two blocks (upper-left and lower-left)
    for i in range(Norbit):
        # upper-left piece
        for j in range(Nbin):
            Cm[j, i] = float(xlib[j, i])

        # lower-left piece
        for irc in range(Nvelb):
            for ivel in range(Nvel):
                j = Nbin + irc * Nvel + ivel
                Cm[j, i] = float(v1lib[ivel, irc, i])

    # Right two pieces (initialize zeros)
    for i in range(Norbit, Norbit + Nvtot):
        for j in range(Nbin + Nvel * Nvelb):
            Cm[j, i] = 0.0

    # Identity entries in lower-right block
    for i in range(Nvtot):
        Cm[Nbin + i, Norbit + i] = 1.0

    return Cm

def getlowright(Cm, w, sumad, Nvel, ircmax, Norbit, Norbitm, Nvelbm, Nbin): 
    """assembles the lower left portion of the matrix Cm which is used to assemble the matrix Am in SPEAR"""
    lds = Norbitm + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    big = 1.0e20

    Nvtot = Nvel * ircmax
    ilast = Norbit + Nvtot      # Fortran's +1 skipped since Python is 0-indexed
    xmu = w[ilast]

    for irc in range(ircmax):
        for ivel in range(Nvel):
            i = Nbin + irc * Nvel + ivel
            Cm[i, ilast] = float(-sumad[ivel, irc])

    return Cm 

def newcmatrix(Cm, xlib, v1lib, inew, Norbit, Nvel, Nvelb, Nbin, Nvelbm):
    """Updates matrix Cm by swapping in orbit 'inew' and reassembling blocks."""
    Nvtot = Nvel * Nvelb
    ilast = Norbit + Nvtot

    # Temporary storage
    xlibtemp = np.zeros(Nbin, dtype=float)
    v1libtemp = np.zeros((Nvel, Nvelb), dtype=float)

    # Copy the orbit column to temporary arrays
    for i in range(Nbin):
        xlibtemp[i] = xlib[i, inew]

    for i in range(Nvelb):
        for ivel in range(Nvel):
            v1libtemp[ivel, i] = v1lib[ivel, i, inew]

    # Swap selected orbit with the last one (Norbit)
    for i in range(Nbin):
        xlib[i, inew] = xlib[i, Norbit]
        xlib[i, Norbit] = xlibtemp[i]

    for i in range(Nvelb):
        for ivel in range(Nvel):
            v1lib[ivel, i, inew] = v1lib[ivel, i, Norbit]
            v1lib[ivel, i, Norbit] = v1libtemp[ivel, i]

    # Upper-left block (density)
    for j in range(Nbin):
        Cm[j, inew] = xlib[j, inew]

    # Lower-left block (velocity)
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            j = Nbin + irc * Nvel + ivel
            Cm[j, inew] = v1lib[ivel, irc, inew]

    # Decrement Norbit (orbit count reduced by one)
    Norbit -= 1

    # Right-side blocks (zeroed region)
    for i in range(Norbit, Norbit + Nvtot):
        for j in range(Nbin + Nvel * Nvelb):
            Cm[j, i] = 0.0

    # Lower-right identity entries
    for i in range(Nvtot):
        Cm[Nbin + i, Norbit + i] = 1.0

    return Cm, xlib, v1lib, Norbit

def spear(Cm, S, dS, ddS, rcond, negweights, apfac2, Norbit, Nvel, Nvelb, Nbin, Nvelbm, w, xlib, v1lib, SMc, sumad, fracnew, apfac):
    """Compute Lagrange multipliers and Newton–Raphson corrections for orbit weights."""

    # --- Define matrix dimensions (same as Fortran parameter() section)
    lds = Norbit + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    Nvtot = Nvel * Nvelb
    Narr = Nbin + Nvtot      # total number of observable constraints
    ilast = Norbit + Nvtot   # total number of active weight elements

    # 1. Construct scaled version of Cm
    Cma = Cm[:Narr, :ilast] / ddS[:ilast]

    # 2. Compute curvature matrix Am = Cm * Cma^T
    Am = Cm[:Narr, :ilast] @ Cma.T

    # 3. Compute residual vector delY
    #     delY contains differences between observed data and model prediction
    #     First part (0:Nbin): surface density residuals
    delY = np.zeros(Narr)
    for ibin in range(Nbin):
        # sum over all orbits contributing to this spatial bin
        sum_val = np.sum(w[:Norbit] * xlib[ibin, :Norbit])
        delY[ibin] = SMc[ibin] - sum_val

    # 4. Second part (Nbin:): velocity distribution residuals
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            i = Nbin + irc * Nvel + ivel
            # predicted velocity contribution from all orbits
            sum_val = np.sum(w[:Norbit] * v1lib[ivel, irc, :Norbit])
            # plus constant velocity term weights
            sum_val += w[Norbit + irc * Nvel + ivel]
            # observed - model difference, scaled by fracnew
            delY[i] = sumad[ivel, irc] * fracnew - sum_val

    # 5. Add the C*dS/ddS correction term
    delY += Cm[:Narr, :ilast] @ (dS[:ilast] / ddS[:ilast])

    # Solve Am * delY = delY  (Fortran's dsico + dsisl)
    # This finds the Lagrange multipliers λ
    lu, piv = lu_factor(Am)
    delY = lu_solve((lu, piv), delY)

    # Adjust orbit weights (Newton–Raphson correction)

    for j in range(ilast):
        sum_val = np.dot(delY, Cm[:Narr, j])
        dw = (sum_val - dS[j]) / ddS[j]
        wnew = w[j] + apfac * dw
        w[j] = wnew

    return w

###################################################################################
# Fundamentals
###################################################################################
def summ(w, xlib, v1lib, Nbin, Norbit, Nvelb, Nvel):
    summed = np.zeros(Nbin + Nvelb * Nvel, dtype=np.float32)

    for ibin in range(Nbin):
        total = 0.0
        for iorb in range(Norbit):
            total += float(w[iorb] * xlib[ibin, iorb])
        summed[ibin] = total

    for irc in range(Nvelb):
        for ivel in range(Nvel):
            i = Nbin + irc * Nvel + ivel
            total = 0.0
            for iorb in range(Norbit):
                total += float(w[iorb] * v1lib[ivel, irc, iorb])
            summed[i] = total

    return summed

def ekin(r, v, Etot, xLz, potential_func): # Kinetic energy
    VV = potential_func(r, v)
    return Etot - VV - (xLz * xLz) / (r * r * (1.0 - v * v) * 2.0)

def force():
    pass

def forcefix():
    pass


def energy(r, th, vr, vth, xLz, potential_func):
    v = np.sin(th)
    vco = np.cos(th)
    vph = xLz / (r * vco)
    VV = potential_func(r, v)
    return 0.5 * (vr**2 + vth**2 + vph**2) + VV


def entropy():
    """Entropy = profit function: entropy term minus χ² terms."""
    pass

def velocity():
    pass

def phase_volume():
    """Calculates orbit phase-space volume."""
    pass

def Lagrangian():
    pass

def potential(r, v, nlegup, npot, rlgpotmin, rlgpotmax, RneeIRhalo, tabv, hole, totlight, gdennorm, Pl):
    rir = (np.log10(r) - rlgpotmin) * (npot - 1) / (rlgpotmax - rlgpotmin) + 1.0
    irlo = int(rir)
    irup = int(rir + 1.0)

    rlo = RneeIRhalo[irlo]
    rup = RneeIRhalo[irup]
    rdiff = r - rlo

    VVlo = 0.0
    VVup = 0.0

    # Linearly interpolate over Legendre orders 0, 2, ..., nlegup
    for n in range(0, nlegup + 1, 2):
        nleg = n // 2
        term = 0.0 if irlo == 0 else Pl(n, v) * tabv[nleg, irlo]
        VVlo += term
        VVup += Pl(n, v) * tabv[nleg, irup]

    VV = VVlo + (VVup - VVlo) / (rup - rlo) * rdiff
    VV = -2.0 * np.pi * VV - hole / (totlight * gdennorm * r)
    return VV

def seecor():
    pass

def radius():
    pass

def density():
    pass

###################################################################################
# Halo Model
###################################################################################

def halodens():
    pass



###################################################################################
# Orbit Model
###################################################################################

def orbits():
    pass

def area():
    pass

def model():
    """Calculates orbit weights matching 2-D profile while maximizing entropy."""
    pass

def projmass():
    """Projects orbit mass distributions."""
    pass


###################################################################################
# Analysis
###################################################################################

def weights():
    pass

def fitting():
    pass


###################################################################################
# Derivative Function (used in RK4 integration)
###################################################################################

def dervis(pos0, xLz, force_func):
    """Compute RK4 derivatives from Legendre polynomial-expanded potential."""
    r, th, vr, vth = pos0
    v = np.sin(th)
    vph = xLz / (r * np.cos(th))

    frL, fvL = force_func(r, v)

    dposdt = np.zeros(4)
    dposdt[0] = vr
    dposdt[1] = vth / r
    dposdt[2] = (vth**2 + vph**2) / r + frL
    dposdt[3] = -np.tan(th) * vph**2 / r - vth * vr / r + fvL
    return dposdt


###################################################################################
# function to compute the LOSVD from orbit weights
###################################################################################
def vor2vphase(sos_sum, dEdLz, I3=None, tiny=1e-10):
    sos_sum = np.array(sos_sum, dtype=float)
    dEdLz = np.array(dEdLz, dtype=float)

    # Handle bad or zero entries
    valid_mask = (sos_sum > 0) & (dEdLz > 0)
    phasevol = np.zeros_like(sos_sum)

    # Compute phase volume: V_i = 1 / (sos_sum * dEdLz)
    with np.errstate(divide="ignore", invalid="ignore"):
        phasevol[valid_mask] = np.abs(1.0 / (sos_sum[valid_mask] * dEdLz[valid_mask]))

    # Replace invalid entries with small fallback value
    phasevol[~valid_mask] = tiny

    # Normalize to match Fortran normalization block
    inv_sum = np.sum(1.0 / np.maximum(phasevol, tiny))
    phasevol *= inv_sum

    return phasevol, valid_mask

def vlook(v1lib, w, Norbit, Nvel, Nvelb):
    """Equivalent to vlook.f — constructs the LOSVD grid from orbit weights."""
    y1 = np.zeros((Nvel, Nvelb))
    for iorb in range(Norbit):
        y1 += w[iorb] * v1lib[:, :, iorb]
    return y1

# Full diagnostic + projection analysis (vlookjm.f)
def vlookjm(v1lib, w, sumad, sadfer, velm, Nvel, Nvelb, Norbit):
    """Equivalent to vlookjm.f — computes projected profiles, H3/H4, and χ²."""

    # === Nested helpers =======================================================
    def get_fwhm(x, y):
        half_max = np.max(y) / 2.0
        indices = np.where(y >= half_max)[0]
        if len(indices) < 2:
            return np.nan
        return x[indices[-1]] - x[indices[0]]

    def gauss_hermite(v, A, v0, sigma, h3, h4):
        """Simple normalized Gauss–Hermite expansion (up to h4)."""
        u = (v - v0) / sigma
        herm3 = (2 * u**3 - 3 * u) / np.sqrt(6)
        herm4 = (4 * u**4 - 12 * u**2 + 3) / np.sqrt(24)
        base = np.exp(-0.5 * u**2)
        return A * base * (1 + h3 * herm3 + h4 * herm4)

    def fit_gauss_hermite(v, f):
        """Fit LOSVD to a Gauss–Hermite expansion."""
        A0 = np.trapz(f, v)
        v0_guess = np.sum(v * f) / np.sum(f)
        sigma_guess = np.sqrt(np.sum((v - v0_guess) ** 2 * f) / np.sum(f))
        popt, _ = curve_fit(
            gauss_hermite, v, f, p0=[A0, v0_guess, sigma_guess, 0.0, 0.0],
            maxfev=10000
        )
        return popt  # (A, v0, sigma, h3, h4)

    def chi_square(model, data, err):
        valid = err > 0
        return np.sum(((model[valid] - data[valid]) / err[valid]) ** 2)

    # === Main computations ====================================================
    y1 = vlook(v1lib, w, Norbit, Nvel, Nvelb)
    sigm, sigd, velmp, veldp, h3m, h4m, h3d, h4d = (
        np.zeros(Nvelb) for _ in range(8)
    )
    chi_arr = np.zeros(Nvelb)

    for ivbin in range(Nvelb):
        y1p = y1[:, ivbin]
        datap = sumad[:, ivbin]
        errv = sadfer[:, ivbin]

        # Model fit
        A, v0, sigma, h3, h4 = fit_gauss_hermite(velm, y1p)
        velmp[ivbin], sigm[ivbin], h3m[ivbin], h4m[ivbin] = v0, sigma, h3, h4

        # Data fit
        A, v0, sigma, h3, h4 = fit_gauss_hermite(velm, datap)
        veldp[ivbin], sigd[ivbin], h3d[ivbin], h4d[ivbin] = v0, sigma, h3, h4

        chi_arr[ivbin] = chi_square(y1p, datap, errv)

    # === Plot example (projected σ) ===========================================
    plt.figure(figsize=(7, 5))
    plt.plot(sigd, label="Observed σ")
    plt.plot(sigm, label="Model σ")
    plt.xlabel("Velocity bin index")
    plt.ylabel("σ (km/s)")
    plt.legend()
    plt.title("Projected Velocity Dispersion")
    plt.show()

    return {
        "velmp": velmp, "sigd": sigd, "sigm": sigm,
        "veldp": veldp, "h3m": h3m, "h4m": h4m,
        "h3d": h3d, "h4d": h4d, "chi2": chi_arr,
    }

###################################################################################
# Plotting / Visualization
###################################################################################

def plot_all():
    pass

def plot_all2():
    pass

def plot_libsos():
    pass


###################################################################################
# Checkers / Filters
###################################################################################

def drop_orbits(w, wphase, Cm, Norbit, Norbito, Norbitm, newcmatrix_func, tiny=1e-30):
    """Filters out negative weights by setting them to small values."""
    iorbiter = np.arange(Norbitm, dtype=int)
    iorb = 0

    while iorb < Norbit:
        if w[iorb] <= tiny:
            # keep pulling last orbit until a non-tiny one remains
            while w[Norbit - 1] <= tiny and Norbit > 1:
                Norbit -= 1

            # swap with last active orbit
            if iorb < Norbit - 1:
                w[iorb], w[Norbit - 1] = w[Norbit - 1], 0.0
                wphase[iorb], wphase[Norbit - 1] = wphase[Norbit - 1], wphase[iorb]

                # update orbit index mapping
                iorbiter[iorb], iorbiter[Norbit - 1] = Norbit - 1, iorb

                # rebuild matrix column
                newcmatrix_func(iorb)

                print(f"Dropped orbit {iorb+1:5d} {Norbit:5d} {Norbito:5d}")

        iorb += 1

    return w, wphase, Cm, Norbit, iorbiter



###################################################################################
# Main library-building routine
###################################################################################

def build_library(GG, totlight, arcsec, Nvel, xorbitmin, xorbitmax, RneeIR, Nrdat):
    """
    Python equivalent of the Fortran PROGRAM library.
    Reads data, builds potential and orbit libraries, and writes output files.
    """

    with open("progressL.out", "w") as f:
        f.write("\n -> reading in data and parameter files\n")

        # --- Read galaxy parameters 
        galaxy_params = read_galaxy(GG, totlight, arcsec, Nvel, xorbitmin, xorbitmax, RneeIR, Nrdat)
        f.write("   - data acquisition complete\n\n")

        # --- Read dark matter halo parameters ---
        halo_params = read_halo()
        f.write(" -> reading dark matter halo parameters\n")
        f.write("   - done\n\n")

        # Setup spatial grid and quadrupole flag
        f.write(" -> setting volume arrays and quadrupole flag\n")
        spacel()
        f.write("   - done\n\n")

        # --- Compute core quantities (requires dL, ratML, vol2d) ---
        f.write(" -> calculating core quantities\n")
        cden = dL[0, 0] * ratML[0, 0] / vol2d[0]
        core = 4.0 * np.pi * cden * rmin**3 / (3.0 + galaxy_params["cslope"])
        f.write("   - done\n\n")

        # --- Calculate density and luminosity distributions ---
        f.write(" -> calculating density distribution\n")
        density()
        write_density()
        write_totml()
        f.write("   - mass and light computed and written to file\n\n")

        # --- Handle extrapolation beyond library field of view ---
        if "iextra" in globals() and iextra == -1:
            rhobeyond()
            f.write("   - depro used outside FOV of the library\n\n")

        # --- Force table calculations ---
        f.write(" -> calculating force tables:\n")
        tables()
        tables_halo()
        tables_add()
        write_tables()
        f.write("   - force tables computed\n\n")

        # --- Non-spherical corrections ---
        if "iquad" in globals() and iquad != 0:
            forcefix()

        # --- Write potential sampling ---
        write_pot()

        # --- Orbit sampling (integral library) ---
        f.write(" -> sampling orbits: calling librarian\n")
        with open("integrals.out", "w"), open("librarian.out", "w"):
            librarian()
        f.write("   - phase space sampled: orbit library created\n\n")

        # --- Write phase-space densities ---
        f.write(" -> writing out phase-space densities\n")
        write_phase()
        f.write("   - phase-space densities written\n\n")

        # --- Write number of orbits to file ---
        with open("norbit", "w") as nfile:
            nfile.write(f"{Norbtot}\n")

        f.write("All tasks completed successfully.\n")

    # --- Return summary of key quantities ---
    return {
    "galaxy": galaxy_params,
    "halo": halo_params,
    "core_density": cden,
    "core_mass": core
}


###############################################################################
# Orchestration
###############################################################################

def librarian():
    """
    Python equivalent of Fortran subroutine LIBRARIAN.
    Samples the phase space available to stars in the potential,
    integrates orbits, and builds the full orbit library.
    """

    # ============================================================
    # SECTION 1. Initialization
    # ============================================================

    print("Initializing orbit library...")
    Norbtot = 0
    wphase[:] = 0.0                      # Phase-space weight array
    idum = -1
    velfac = np.sqrt(GG * totlight * arcsec / distance / angrad)
    tooclose = 0.01
    efac = efac0                         # Integration step control
    v_losvd = 1.0 / vmult
    vstep_min = v_losvd * abs(vfrac)
    facnorb = 1.0
    n_rs0 = Nl                           # number of radial steps
    Eo = 0.0                             # reference energy baseline

    # --- Create output files ---
    open_files = {
        "xlib": open("xlib.out", "w"),
        "v1lib": open("v1lib.out", "w"),
        "vani": open("vani.out", "w"),
        "sumv2": open("sumv2.out", "w"),
        "d3": open("d3.out", "w"),
        "highervani": open("highervani.out", "w"),
        "integrals": open("integrals.out", "w"),
        "librarian": open("librarian.out", "w"),
        "phase": open("phase.out", "w"),
    }

    # ============================================================
    # SECTION 2. Angular momentum (Lz) and Energy loops
    # ============================================================

    for iLz in range(Lz_min, Lz_max, iLz_step):
        rperi = RneeIRperi(iLz)

        for iE in range(iLz + ideltaE, Nrelz, iE_step):
            # --- Compute pericenter/apocenter and potential ---
            r_apo = RneeIRapo(iE)
            if rperi > 1.0 or r_apo < rmin:
                continue  # skip outside range

            phi_peri = potential(rperi, 0.0)
            phi_apo = potential(r_apo, 0.0)

            if phi_peri >= phi_apo:
                phi_peri, phi_apo = phi_apo, phi_peri

            # --- Energy and angular momentum setup ---
            xLz = np.sqrt(2.0 * abs(phi_apo - phi_peri) /
                         (r_apo**2 - rperi**2)) * r_apo * rperi
            Etot = phi_peri + 0.5 * xLz**2 / rperi**2
            dE = Etot - Eo
            iorblo = Norbtot + 1

            # --- Initialize radial sampling ---
            n_rs = 3 if iE == iLz else n_rs0
            norbs = 0
            rsmin, rsmax, vsmin, vsmax = [], [], [], []

            # ====================================================
            # SECTION 3. Launch and integrate orbits
            # ====================================================

            for i_rs in range(n_rs):
                r = radius(rperi, r_apo, i_rs, n_rs)
                if ekin(r, 0.0) < 0:
                    continue

                vt_up = np.sqrt(2.0 * ekin(r, 0.0))
                vr = 0.001
                vth = np.sqrt(max(2.0 * ekin(r, 0.0) - vr**2, 0))
                th = 0.0

                pos0 = np.array([r, th, vr, vth])

                # --- Compute starting energy ---
                Estart = energy(r, th, vr, vth)

                # --- Orbit integration ---
                orbit_result = integrate_orbit(
                    pos0=pos0,
                    efac=efac,
                    epsilon=0.05,
                    xLz=xLz,
                    Etot=Etot
                )

                # --- Check energy conservation and escape ---
                deltaE = orbit_result["deltaE"]
                if abs(deltaE) > deltaE_lim:
                    continue  # skip unstable orbit
                if orbit_result["escaped"]:
                    continue

                # ====================================================
                # SECTION 4. Record orbit data
                # ====================================================

                Norbtot += 1
                norbs += 1

                write_xlib(orbit_result, open_files["xlib"])
                write_v1lib(orbit_result, open_files["v1lib"])
                write_vani(orbit_result, open_files["vani"])
                write_d3(orbit_result, open_files["d3"])
                write_sumv2(orbit_result, open_files["sumv2"])
                write_highervani(orbit_result, open_files["highervani"])
                write_integrals(orbit_result, open_files["integrals"])
                write_phase(orbit_result, open_files["phase"])

            # ====================================================
            # SECTION 5. Post-sequence phase-space volume
            # ====================================================

            compute_phase_volume(iorblo, Norbtot, dLz, dE)

            for i in range(iorblo, Norbtot + 1):
                wphase[i] = 1.0 / max(wphase[i], 1e-12)

    # ============================================================
    # SECTION 6. Finalization
    # ============================================================

    for f in open_files.values():
        f.close()

    np.savetxt("norbit", [Norbtot])
    print(f"Orbit library built successfully with {Norbtot} orbits.")
    return Norbtot

