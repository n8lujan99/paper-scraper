'''
ospm/
  libshim.py              # small: constants, fac1, grid registration, Pl(), mappers
  tables_state.py         # *only* shared arrays/sizes for tables & densities
  quadpack.py             # qags/qagse/qk21/qelg/qpsrt/qgaus (the long bits)
  tables.py               # tables(), denny(), tables_read/write (uses quadpack & shim)
  potential.py            # potential(), force(), potwrite() (reads tables_state & shim)
  orbits.py               # derivs(), rk4_step(), step(), energy(), ekin(), run_single_orbit()
  io_utils.py             # iterread, mlread, phasewrite/read, weight_read/write, xlib_read, vaniread
  linear_model.py         # build_cmatrix(), new_cmatrix_update(), summ()
  core.py                 # the control hub / CLI entrypoint
'''

import os 
import math
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval

from libshim import (
    # physical / unit constants
    GG, arcsec, distance, angrad, totlight, vmult, vadd, fac1, epsilon, efac,
    IRneeR, RneeIR, IVneeV, VneeIV,tables_read, tables_write, xlib_read, 
    weight_read, weight_write, vaniread, phaseread, phasewrite, 
    iterread, mlread, build_cmatrix, summ, potential, energy, derivs, step,
)
# The purpose of this document will be to pull all needed files and libraries for the OSPM analysis.
def main():
    ap = argparse.ArgumentParser(
        description="Control hub for OSPM: load tables/weights/library and do basic checks."
    )
    ap.add_argument("--tables", default="tables.out", type=Path, help="path to tables.out")
    ap.add_argument("--weights", default="weights.out", type=Path, help="path to weights.out")
    ap.add_argument("--xlib", default="xlib.out", type=Path, help="path to xlib.out")
    ap.add_argument("--vani", default="vani.out", type=Path, help="path to vani.out (optional)")

    # geometry of your model/library (provide once here)
    ap.add_argument("--Nbin", required=True, type=int, help="number of spatial bins (model)")
    ap.add_argument("--Nvlib", required=True, type=int, help="Nvlib (coarse theta bins) in xlib.out first block")

    # velocity-grid layout for nuisance/LOS blocks
    ap.add_argument("--Nvel", required=True, type=int, help="number of velocity bins per rc-bin")
    ap.add_argument("--Nvelb", required=True, type=int, help="number of radial cells for the velocity grid")

    # library sizing (pull Norb from phase files if present, else allow override)
    ap.add_argument("--norbit", type=int, help="Norbit (model columns). If omitted, inferred from phase/norbit.")
    ap.add_argument("--norb", type=int, help="Norb (library columns). If omitted, inferred from phase/norbit.")
    ap.add_argument("--phase", default="phase.out", type=Path, help="phase.out for optional Norb inference")
    ap.add_argument("--norbit_file", default="norbit", type=Path, help="file with Norb (one int)")

    ap.add_argument("--print-only", action="store_true", help="load & report shapes; don’t do any math")

    args = ap.parse_args()

    # 1) Try to infer Norb/Norbit from phase/norbit (matches your phaseread translation)
    Norb, Norbit = None, None
    try:
        ph = phaseread(norbit_path=args.norbit_file, phase_path=args.phase)
        Norb, Norbit = ph.Norb, ph.Norbit
    except Exception:
        pass
    if args.norb is not None:   Norb = args.norb
    if args.norbit is not None: Norbit = args.norbit
    if Norb is None or Norbit is None:
        raise SystemExit("Could not infer Norb/Norbit. Pass --norb/--norbit or provide phase.out + norbit.")

    # 2) Read potential/force tables
    tabv, tabfr, tabfv = tables_read(args.tables)

    # 3) Read xlib.out → smlib (Nbin, Norbit) with smearing/duplication rules
    smlib = xlib_read(path=args.xlib, Norb=Norb, Norbit=Norbit, Nbin=args.Nbin, Nvlib=args.Nvlib)

    # 4) Read weights
    w = weight_read(path=args.weights, Norb=Norb, Norbit=Norbit, Nvel=args.Nvel, Nvelb=args.Nvelb)

    # 5) Optional: read vani.out (velocity moments per orbit)
    v_read_ok = False
    v0a = vr2a = vt2a = vp1a = vp2a = None
    if args.vani and Path(args.vani).exists():
        # You know Nrani and Ntani; add flags if you want to pipe those through.
        # For now we try to infer by divisibility check later or skip.
        print("[vani] Provide Nrani/Ntani in libshim or extend this CLI to pass them.")
    else:
        print("[vani] Skipping (file not found).")

    # 6) Report shapes & a quick dot (SUMM) if we have everything
    print("== Constants ==")
    print(f"GG={GG}, arcsec={arcsec}, distance={distance}, angrad={angrad}, totlight={totlight}")
    print(f"fac1={fac1}, epsilon={epsilon}, efac={efac}")
    print("== Sizes ==")
    print(f"Norb={Norb}, Norbit={Norbit}, Nbin={args.Nbin}, Nvlib={args.Nvlib}, Nvel={args.Nvel}, Nvelb={args.Nvelb}")
    print("== Loaded ==")
    print(f"tabv/tabfr/tabfv: {tabv.shape} / {tabfr.shape} / {tabfv.shape}")
    print(f"smlib: {smlib.shape}, weights: {w.shape}")

    if args.print-only:
        return

    # If you have a v1lib array available (Nvel, Nvelb, Norbit), you can compute Cm or summed:
    #   summed = summ(w, xlib=smlib, v1lib=v1lib)
    # For now, just demonstrate the spatial part of SUMM (works without v1lib):
    x_part = smlib @ w[:Norbit]
    print(f"spatial SUMM (first 10): {np.asarray(x_part[:10]).ravel()}")

if __name__ == "__main__":
    main()