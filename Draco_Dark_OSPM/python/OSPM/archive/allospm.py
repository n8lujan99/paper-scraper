###################################################################################
# Imports
###################################################################################

import numpy as np 
import scipy as sp
import pandas as pd 
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.linalg import lu_factor, lu_solve


###################################################################################
# Readers
###################################################################################
def galaxyread(GG, totlight, arcsec, Nvel, xorbitmin, xorbitmax, RneeIR, Nrdat):
    """
    Reads galaxy parameters from 'galaxy.params' and 'gal.dat'.
    Combines functionality of Fortran GALAXYREAD and GALREADM.
    """

    # --- Read both input files ---
    gal_params = pd.read_csv("galaxy.params", delim_whitespace=True, header=None, comment="#")
    gal_dat = pd.read_csv("gal.dat", delim_whitespace=True, header=None, comment="#")

    # --- Required shared parameters ---
    hole = float(gal_dat.iloc[0, 1])
    xinclin = float(gal_dat.iloc[1, 1])

    distance = float(gal_params.iloc[1, 1])
    angrad   = float(gal_params.iloc[2, 1])
    vmin     = float(gal_params.iloc[3, 1])
    vmax     = float(gal_params.iloc[4, 1])
    cslope   = float(gal_params.iloc[5, 1])

    # --- Optional fields (galreadm-only) ---
    xmtol = None
    islitsize = None
    if len(gal_params) > 6:
        try:
            xmtol = float(gal_params.iloc[6, 1])
            islitsize = float(gal_params.iloc[7, 1])
        except Exception:
            pass

    # --- Derived parameters ---
    xinclin_rad = np.deg2rad(xinclin)
    sini = np.sin(xinclin_rad)
    cosi = np.cos(xinclin_rad)
    distance *= 1.0e6  # Convert from Mpc → pc

    fac1 = np.sqrt(GG * totlight * arcsec / (distance * angrad))
    vmult = (float(Nvel - 1) / (vmax - vmin)) * fac1
    vadd = 1.0 - float(Nvel - 1) * vmin / (vmax - vmin)

    perimin = xorbitmin / angrad
    apomax = max(xorbitmax / angrad, 1.2)

    rpotmin = min(perimin / 10.0, RneeIR[0] / 10.0)
    rpotmax = max(apomax * 2.0, RneeIR[Nrdat - 1] * 2.0)
    rlgpotmin = np.log10(rpotmin)
    rlgpotmax = np.log10(rpotmax)

    # --- If xmtol is present, compute xmutarget (galreadm feature) ---
    xmutarget = 1.0 / np.sqrt(xmtol) if xmtol is not None else None

    # --- Package results ---
    results = { "hole": hole, "xinclin": xinclin_rad, "sini": sini, "cosi": cosi, "distance": distance, 
            "angrad": angrad, "vmin": vmin, "vmax": vmax, "cslope": cslope, "fac1": fac1, "vmult": vmult, 
            "vadd": vadd, "perimin": perimin, "apomax": apomax, "rpotmin": rpotmin, "rpotmax": rpotmax, 
            "rlgpotmin": rlgpotmin, "rlgpotmax": rlgpotmax }


    if xmtol is not None:
        results["xmtol"] = xmtol
        results["xmutarget"] = xmutarget
    if islitsize is not None:
        results["islitsize"] = islitsize

    return results

def read_halo(filename_halo="halo.dat", filename_gden="gden.norm"):
    with open(filename_halo, "r") as f:
        ihalo = int(f.readline().split()[0])
        dis = float(f.readline().split()[0])
        v0, rc, qdm = map(float, f.readline().split())
        xmgamma, rsgamma, gamma = map(float, f.readline().split())
        cnfw, rsnfw = map(float, f.readline().split())
    with open(filename_gden, "r") as f:
        gdennorm = float(f.readline().split()[0])
    return { "ihalo": ihalo, "dis": dis, "v0": v0, "rc": rc, "qdm": qdm, "xmgamma": xmgamma, 
            "rsgamma": rsgamma, "gamma": gamma, "cnfw": cnfw, "rsnfw": rsnfw, "gdennorm": gdennorm }


def read_weight(w, Norbit, Norbitm, Nvel, Nvelb, filename="weights.out"):
    ilast = Norbit + Nvel * Nvelb
    wt = np.zeros(Norbitm + Nvel * Nvelb)
    try:
        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                ii, val = int(parts[0]), float(parts[1])
                if 1 <= ii <= len(wt):
                    wt[ii - 1] = val
    except FileNotFoundError:
        raise FileNotFoundError(f"{filename} not found")

    if np.any(wt[:Norbit] == 0):
        for i in range(Norbit):
            w[2*i] = wt[i] / 2.0
            w[2*i + 1] = wt[i] / 2.0
        for i in range(Nvel * Nvelb + 1):
            w[Norbit + i] = wt[Norbit + i]
    else:
        for i in range(ilast):
            w[i] = wt[i]
    return w


def read_phase(wphase, Norbm, Norbitm, filename_phase="phase.out", filename_norbit="norbit"):
    with open(filename_norbit, "r") as f:
        Norb = int(f.readline().split()[0])
    print(f"Norbit = {Norb}")
    if Norb > Norbm:
        print("make Norbm bigger")
    Norbit = 2 * Norb if Norbm != Norbitm else Norb

    with open(filename_phase, "r") as f:
        slush = np.fromstring(f.readline().strip(), sep=" ")

    for i in range(Norb):
        if Norbm != Norbitm:
            wphase[2*i] = slush[i]
            wphase[2*i + 1] = slush[i]
        else:
            wphase[i] = slush[i]

    wmin = np.min(wphase[:Norbit])
    wmax = np.max(wphase[:Norbit])
    imin = np.argmin(wphase[:Norbit])
    imax = np.argmax(wphase[:Norbit])
    # print(f"min and max phase vol : {imin+1} {wmin} {imax+1} {wmax}")
    return { "Norb": Norb, "Norbit": Norbit, "wmin": wmin, "wmax": wmax, "imin": imin+1, "imax": imax+1 }


def read_ML(filename="ml.out"):
    with open(filename, "r") as f:
        xmu = float(f.readline().split()[0])
        alp, ent, chi = map(float, f.readline().split())
        x1, x2, ent, chi = map(float, f.readline().split())
    xmu = 1.0 / np.sqrt(xmu)
    return { "xmu": xmu, "alp": alp, "ent": ent, "chi": chi, "x1": x1, "x2": x2 }


def read_iterr(filename="iter.params"):
    with open(filename, "r") as f:
        _, Niter     = f.readline().split()
        _, apfac     = f.readline().split()
        _, apfacmu   = f.readline().split()
        _, ifit      = f.readline().split()
        _, alphainc  = f.readline().split()
        _, alphastop = f.readline().split()
        _, fracrv    = f.readline().split()

    return { "Niter": int(Niter), "apfac": float(apfac), "apfacmu": float(apfacmu), "ifit": int(ifit), "alphainc": float(alphainc),
            "alphastop": float(alphastop), "fracrv": float(fracrv) }


def read_tables(tabv, tabfr, tabfv, nlegup, npot, filename="tables.out"):
    with open(filename, "r") as f:
        for n in range(0, nlegup + 1, 2):
            nleg = n // 2 + 1
            for irtab in range(npot):
                i1, i2, x1 = map(float, f.readline().split())
                i1, i2 = int(i1), int(i2)
                tabv[i1, i2] = x1

                i1, i2, x1 = map(float, f.readline().split())
                i1, i2 = int(i1), int(i2)
                tabfr[i1, i2] = x1

                i1, i2, x1 = map(float, f.readline().split())
                i1, i2 = int(i1), int(i2)
                tabfv[i1, i2] = x1
    return tabv, tabfr, tabfv


def read_vani(v0a, vr2a, vt2a, vp1a, vp2a, Nrani, Ntani, Norb, Norbit, filename="vani.out"):
    with open(filename, "r") as f:
        for iorb in range(Norb):
            line = f.readline()
            if not line: raise ValueError(f"Missing data for orbit {iorb} in {filename}")
            data = np.fromstring(line.strip(), sep=" "); expected = 6 * Nrani * Ntani
            if data.size < expected: raise ValueError(f"Orbit {iorb}: insufficient data (got {data.size}, expected {expected})")
            stride = Nrani * Ntani
            v0liba = data[0*stride:1*stride].reshape(Nrani, Ntani)
            vr2liba = data[1*stride:2*stride].reshape(Nrani, Ntani)
            vt2liba = data[2*stride:3*stride].reshape(Nrani, Ntani)
            vp1liba = data[4*stride:5*stride].reshape(Nrani, Ntani)
            vp2liba = data[5*stride:6*stride].reshape(Nrani, Ntani)
            for ir in range(Nrani):
                for it in range(Ntani):
                    if Norb != Norbit:
                        v0a[ir,it,2*iorb] = v0liba[ir,it]; v0a[ir,it,2*iorb+1] = v0liba[ir,it]
                        vr2a[ir,it,2*iorb] = vr2liba[ir,it]; vr2a[ir,it,2*iorb+1] = vr2liba[ir,it]
                        vt2a[ir,it,2*iorb] = vt2liba[ir,it]; vt2a[ir,it,2*iorb+1] = vt2liba[ir,it]
                        vp1a[ir,it,2*iorb] =  vp1liba[ir,it]; vp1a[ir,it,2*iorb+1] = -vp1liba[ir,it]
                        vp2a[ir,it,2*iorb] = vp2liba[ir,it]; vp2a[ir,it,2*iorb+1] = vp2liba[ir,it]
                    else:
                        v0a[ir,it,iorb] = v0liba[ir,it]; vr2a[ir,it,iorb] = vr2liba[ir,it]
                        vt2a[ir,it,iorb] = vt2liba[ir,it]; vp1a[ir,it,iorb] = vp1liba[ir,it]; vp2a[ir,it,iorb] = vp2liba[ir,it]
    return v0a, vr2a, vt2a, vp1a, vp2a

def highervaniread(Nrani, Ntani, Norb, Norbit, vra, vta, vrtpa, vrpa, vtpa, filename="highervani.out"):
    """
    Reads the higher-order velocity library (from highervani.out)
    and fills the 3D velocity component arrays.
    """
    with open(filename, "r") as f:
        for iorb in range(Norb):
            # Each line in Fortran read(11,*) vrliba, vtliba, vrtpliba, vrpliba, vtpliba
            # Flatten the read and reshape back into (Nrani, Ntani)
            data = np.fromstring(f.readline().strip(), sep=" ")
            expected_size = 5 * Nrani * Ntani
            if data.size < expected_size:
                raise ValueError(f"Orbit {iorb}: insufficient data in {filename}")

            # Split into separate blocks
            stride = Nrani * Ntani
            vrliba  = data[0*stride : 1*stride].reshape(Nrani, Ntani)
            vtliba  = data[1*stride : 2*stride].reshape(Nrani, Ntani)
            vrtpliba = data[2*stride : 3*stride].reshape(Nrani, Ntani)
            vrpliba  = data[3*stride : 4*stride].reshape(Nrani, Ntani)
            vtpliba  = data[4*stride : 5*stride].reshape(Nrani, Ntani)

            # Fill target arrays
            if Norb != Norbit:
                vra[:, :, 2*iorb]     = vrliba
                vra[:, :, 2*iorb + 1] = vrliba
                vta[:, :, 2*iorb]     = vtliba
                vta[:, :, 2*iorb + 1] = vtliba
                vrtpa[:, :, 2*iorb]     = vrtpliba
                vrtpa[:, :, 2*iorb + 1] = vrtpliba
                vrpa[:, :, 2*iorb]     = vrpliba
                vrpa[:, :, 2*iorb + 1] = vrpliba
                vtpa[:, :, 2*iorb]     = vtpliba
                vtpa[:, :, 2*iorb + 1] = vtpliba
            else:
                vra[:, :, iorb]   = vrliba
                vta[:, :, iorb]   = vtliba
                vrtpa[:, :, iorb] = vrtpliba
                vrpa[:, :, iorb]  = vrpliba
                vtpa[:, :, iorb]  = vtpliba

    return vra, vta, vrtpa, vrpa, vtpa

def read_xlib(Nvlib, Nrlib, Nbin, Norb, Norbit, smlib, filename="xlib.out"):
    """
    Reads the orbit library from 'xlib.out' and smears out the lowest-r bin
    information over all angles.
    """
    # Open and read file 
    with open(filename, "r") as f:
        for iorb in range(Norb):
            # Read one orbit’s data into slush
            line = f.readline()
            if not line:
                raise ValueError(f"Missing data for orbit {iorb} in {filename}")
            slush = np.fromstring(line.strip(), sep=" ")

            # Initialize first bin
            if Norb != Norbit:
                smlib[0, 2 * iorb] = 0.0
                smlib[0, 2 * iorb + 1] = 0.0
            else:
                smlib[0, iorb] = 0.0

            # --- Smear out first coarse bin (sum over velocity bins) ---
            for iv in range(Nvlib):
                if Norb != Norbit:
                    smlib[0, 2 * iorb] += slush[iv]
                    smlib[0, 2 * iorb + 1] += slush[iv]
                else:
                    smlib[0, iorb] += slush[iv]

            # --- Copy remaining bins directly ---
            for ibin in range(1, Nbin):
                src_val = slush[ibin + Nvlib - 1]
                if Norb != Norbit:
                    smlib[ibin, 2 * iorb] = src_val
                    smlib[ibin, 2 * iorb + 1] = src_val
                else:
                    smlib[ibin, iorb] = src_val

    return smlib

def d3read(Nvlib, Nrlib, Nbin, Norb, Norbit, d3lib):
    '''
    Reads the 3D orbit library from 'd3.out' and distributes
    (smears) the first coarse theta-bin over all velocity angles.
    '''
    # --- Load file ---
    with open("d3.out", "r") as f:
        for iorb in range(Norb):
            slush = np.fromstring(f.readline().strip(), sep=" ")
            if slush.size < Nvlib * Nrlib:
                raise ValueError(f"Incomplete line in d3.out for orbit {iorb}")

            # Initialize first bin
            if Norb != Norbit:
                d3lib[0, 2 * iorb] = 0.0
                d3lib[0, 2 * iorb + 1] = 0.0
            else:
                d3lib[0, iorb] = 0.0

            # --- Smear first coarse bin (sum over velocity bins) ---
            for iv in range(Nvlib):
                if Norb != Norbit:
                    d3lib[0, 2 * iorb] += slush[iv]
                    d3lib[0, 2 * iorb + 1] += slush[iv]
                else:
                    d3lib[0, iorb] += slush[iv]

            # --- Copy remaining bins directly ---
            for ibin in range(1, Nbin):
                src_val = slush[ibin + Nvlib - 1]
                if Norb != Norbit:
                    d3lib[ibin, 2 * iorb] = src_val
                    d3lib[ibin, 2 * iorb + 1] = src_val
                else:
                    d3lib[ibin, iorb] = src_val

    return d3lib

def read_vlib(Nvel, Nrlib, Nvlib, Nvelb, Norbit, Nstot, areakin, itobin, isee, iminor, angrad, v1lib, seeb=None):
    """
    Reads velocity library files (v1lib.out, v1libsc*.out) and assembles v1lib array.
    """
    # --- Define possible file list
    file_list = [f"v1libsc{i}.out" for i in range(1, 10)]

    # --- Read seeing data (see.dat)
    if seeb is None:
        see_df = pd.read_csv("see.dat", delim_whitespace=True, header=None, skiprows=1)
        seeb = see_df[1].values / angrad

    # --- Initialize arrays
    v1lib.fill(0.0)
    slush = np.zeros(Nvel * Nrlib * Nvlib)
    slusha = np.zeros((Nvel * Nrlib * Nvlib, Nstot))
    izerofwhm = 0

    # --- Open base velocity library
    base_file = open("v1lib.out", "r")

    # --- Open all seeing-convolved libraries once (like Fortran)
    open_files = []
    for i in range(Nstot):
        if seeb[i] > 0:
            fname = file_list[len(open_files)]
            open_files.append(open(fname, "r"))
        else:
            izerofwhm = 1

    # --- Main orbit loop
    for iorb in range(Norbit):
        # Read one line (one orbit) from base file
        slush_line = base_file.readline()
        if not slush_line:
            break
        slush = np.fromstring(slush_line.strip(), sep=" ")

        # Read corresponding lines from each seeing-convolved file
        for isee_idx, f in enumerate(open_files):
            line = f.readline()
            if not line:
                continue
            slush2 = np.fromstring(line.strip(), sep=" ")
            slusha[:, isee_idx] = slush2

        # --- Map data into v1lib
        for ivbin in range(Nvelb):
            for ir in range(Nrlib):
                for iv in range(Nvlib):
                    fracm = areakin[ivbin, itobin[ir, iv]]
                    for ivel in range(Nvel):
                        ic = ir * Nvlib * Nvel + iv * Nvel + ivel
                        if seeb[isee[ivbin]] == 0:
                            value = slush[ic]
                        else:
                            value = slusha[ic, isee[ivbin] - izerofwhm]
                        v1lib[ivel, ivbin, iorb] += value * fracm

    # --- Close all open files (after loops)
    base_file.close()
    for f in open_files:
        f.close()


def sos2nuclei(Norbit, Nsos=55, nquadlength=10000, dist_lim=0.1, fracshift=0.001, minshift=1e-7, maxshift=1e5, tooclose=1e-8, do_plot_elz=False):
    # === helper: read SOS orbit data (sos1.out/sos2.out) ======================
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


###################################################################################
# Writers
###################################################################################

def write_vel # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_Gal_params # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_orbit_lib # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_weight # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_phase # to Phase.out
    pd.DataFrame().to_csv('Phase.out', sep=' ', header=None, index=None)

def write_M/L # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_inter # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_table # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_compare # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def vor2vphase #CELLS2VPHASE: transforms the output of the voronoi-tesselation into phasevols of the orbits


###################################################################################
# Binning
# Functions returning physical quantities from bin indices correspond to BIN CENTER values.
###################################################################################

def IRneeR():
    """Converts physical r (on [0,1]) into radial bin index ir."""
    pass

def RneeIR():
    """Converts radial bin index ir into physical r."""
    pass

def IVneeV():
    """Converts physical v=sin(th) into angular bin index iv."""
    pass

def VneeIV():
    """Converts angular bin index iv into physical v = sin(th)."""
    pass

def mass_bin():
    """Calculates mass in each bin from surface brightness and flattening."""
    pass

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

def force():
    pass

def energy():
    pass

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

def Potential():
    pass

def seecor():
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
# Librarians
###################################################################################

def Librarian(Orbit):
    """Samples phase space and feeds coordinates to ORBIT."""
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


