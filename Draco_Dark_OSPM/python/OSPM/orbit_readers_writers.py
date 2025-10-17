###################################################################################
# Imports
###################################################################################

import numpy as np 
import pandas as pd 


###################################################################################
# Readers
##############################################################

def read_galaxy(GG, totlight, arcsec, Nvel, xorbitmin, xorbitmax, RneeIR, Nrdat):
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

def write_weight(w, Norbit, Nvel, Nvelb, filename="weights.out"):
    with open(filename, "w") as f:
        for i in range(1, Norbit + Nvel * Nvelb + 1):
            f.write(f"{i:6d}  {w[i-1]:12.6e}\n")
    return


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

def write_phase(wphase, Norbtot, filename="phase.out"):
    with open(filename, "w") as f:
        f.write(" ".join(f"{wphase[i]:.6e}" for i in range(Norbtot)) + "\n")
    return


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

def write_tables(tabv, tabfr, tabfv, nlegup, npot, filename="tables.out"):
    with open(filename, "w") as f:
        for n in range(0, nlegup + 1, 2):
            nleg = n // 2 + 1
            for irtab in range(1, npot + 1):
                f.write(f"{nleg:4d}  {irtab:4d}  {tabv[nleg, irtab]:12.6e}\n")
                f.write(f"{nleg:4d}  {irtab:4d}  {tabfr[nleg, irtab]:12.6e}\n")
                f.write(f"{nleg:4d}  {irtab:4d}  {tabfv[nleg, irtab]:12.6e}\n")
    return

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

def read_highervani(Nrani, Ntani, Norb, Norbit, vra, vta, vrtpa, vrpa, vtpa, filename="highervani.out"):
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

def read_d3(Nvlib, Nrlib, Nbin, Norb, Norbit, d3lib):
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

def write_pot(potential_func, force_func, GG, totlight, arcsec, distance, angrad, b, rpotmin, rpotmax, npot=1000):
    with open("potential.out", "w") as f:
        ddd = np.log(3.6 / b) / 999.0
        ccc = b / (3.0 * np.exp(ddd))
        for irr in range(1, 1001):
            r = ccc * np.exp(ddd * irr)
            VV0 = potential_func(r, 0.0)
            VV5 = potential_func(r, 0.5)
            VV9 = potential_func(r, 0.9)
            f.write(f"{r:12.6e}  {VV0:12.6e}  {VV5:12.6e}  {VV9:12.6e}\n")

        r = 0.5
        for ith in range(900):
            th = float(ith) / 10.0 * np.pi / 180.0
            v = np.sin(th)
            VV = potential_func(r, v)
            f.write(f"{180.0 * th / np.pi:12.6e}  {-VV:12.6e}\n")

    with open("vcirc.out", "w") as f2:
        xfac = np.sqrt(GG * totlight * arcsec / (distance * angrad))
        xmin = np.log10(rpotmin)
        xmax = np.log10(rpotmax)
        ipnts = 200
        for ir in range(1, ipnts + 1):
            xr = 10 ** (xmin + (xmax - xmin) / float(ipnts - 1) * float(ir - 1))
            xrforce, xvforce = force_func(xr, 0.0)
            f2.write(f"{ir:4d}  {xr:12.6e}  {np.sqrt(abs(xrforce) * xr) * xfac:12.6e}\n")

    return

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

def write_compare # to 
    pd.DataFrame().to_csv('', sep=' ', header=None, index=None)

def write_density
    pass

def write_totml
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
# Plotting / Visualization
###################################################################################

def plot_all():
    pass

def plot_all2():
    pass

def plot_libsos():
    pass