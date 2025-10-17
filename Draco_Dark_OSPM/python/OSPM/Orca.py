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

