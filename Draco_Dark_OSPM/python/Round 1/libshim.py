#libshim.py

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

# =================== Paths ===================
CWD          = Path.cwd()
OUTDIR       = (CWD / "products") if (CWD / "products").exists() else (CWD / "research" / "OSPM" / "products")
STAMP        = "20250820"  # YYYYMMDD for the match you want to open
QA_FILE      = OUTDIR / f"draco_spec_qamatch_{STAMP}.csv"
STARLIST_MU  = OUTDIR / f"draco_ospm_starlist_mu_{STAMP}.csv"

# =================== Column order ===================
QA_COLS_ORDER = ['RA_deg','Dec_deg','HJD', 'Vlos','e_Vlos','s_Vlos', '[Fe/H]','Teff','logg',
    'sep_arcsec','opt_source_id','opt_RA_deg','opt_Dec_deg', 'opt_radial_velocity',
    'opt_radial_velocity_error','opt_phot_g_mean_mag','opt_bp_rp','opt_ruwe','opt_rbin_ell_key',
    'opt_tbin_ell_key','opt_r_pc_ell','opt_rbin_circ_key','opt_r_pc_circ']

STARLIST_COLS = ['opt_source_id','RA_deg','Dec_deg','Vlos','e_Vlos', 'opt_r_pc_ell','opt_rbin_ell_key','opt_tbin_ell_key']

@dataclass(frozen=True)
class Config: GG: float; arcsec: float; distance: float; angrad: float; totlight: float; vmult: float; vadd: float
@dataclass
class State: r: float; th: float; vr: float; vth: float
@dataclass(frozen=True)
class BinnerParams: a: float; b: float; cval: float; rmin: float; irmax: int; ivmax: int; irrat: int; ivrat: int; aelz: float; belz: float; celz: float; rlgpotmin: float; rlgpotmax: float; npot: int
@dataclass
class HaloParams: ihalo: int; dis: float; v0: float; rc: float; qdm: float; xmgamma: float; rsgamma: float; gamma: float; cnfw: float; rsnfw: float; gdennorm: float
@dataclass
class IterParams: Niter: int; apfac: float; apfacmu: float; ifit: int; alphainc: float; alphastop: float; fracrv: float
@dataclass
class MLParams: xmu: float; alp: float; ent: float; chi: float; x1: float; x2: float
@dataclass
class PhaseReadResult: wphase: np.ndarray; Norb: int; Norbit: int; wmin: float; wmax: float; imin: int; imax: int

# ===== hard-wired constants =====
GG        = 4.30091e-6
arcsec    = 1.0 / 206265.0
distance  = 1.0
angrad    = arcsec * math.pi / 180.0
totlight  = 1.0
vmult     = 1.0
vadd      = ???      # from tablesread()
epsilon   = 0.025    # time step factor to mirror oneorbit
efac      = 1.0

# derived conversion used all over Fortran: fac1 = sqrt(GG*totlight*arcsec/distance/angrad)
fac1 = math.sqrt(GG * totlight * arcsec / distance / angrad)

# ===== minimal table container with the four mapper methods =====
@dataclass
class Tables:
    Rgrid: np.ndarray
    Vgrid: np.ndarray

    def IRneeR(self, r: float) -> int:
        i = int(np.searchsorted(self.Rgrid, r, side="right"))
        return max(1, min(i, len(self.Rgrid)))  # 1-based

    def RneeIR(self, i: int) -> float:
        idx = max(1, min(i, len(self.Rgrid))) - 1
        return float(self.Rgrid[idx])

    def IVneeV(self, v: float) -> int:
        a = abs(v)
        i = int(np.searchsorted(self.Vgrid, a, side="right"))
        return max(1, min(i, len(self.Vgrid)))  # unsigned bin; add sign in your caller if needed

    def VneeIV(self, i: int) -> float:
        idx = max(1, min(i, len(self.Vgrid))) - 1
        return float(self.Vgrid[idx])

TAB   = Tables(Rgrid=np.geomspace(1e-3, 1e+1, 256), Vgrid=np.linspace(0.0, 1.0, 256))
IRneeR, RneeIR, IVneeV, VneeIV = TAB.IRneeR, TAB.RneeIR, TAB.IVneeV, TAB.VneeIV
print("Initialized")

##########################################################################################################################
# =================== data I/O (QA/starlist) ===================
pd.options.display.width = 160
pd.options.display.max_columns = 200

# ---- table & density sizes (filled after loading) ----
nlegup = None
npot   = None
nlegdim = None
Nrdat  = None
Nvdat  = None
Nrbey  = None

# ---- arrays (assigned by tables_read / density loaders) ----
tabv = tabfr = tabfv = None     # shapes (nlegdim, npot)
rho = rhobey = None             # shapes (Nrdat, Nvdat) and (Nrbey, Nvdat)

# ---- grids & mapping helpers ----
class _HaloGrid:
    def __init__(self, rhalo):
        self.r = np.asarray(rhalo, float)
        self.log10 = np.log10(self.r)
        self.rlgmin = float(self.log10.min())
        self.rlgmax = float(self.log10.max())
    def RneeIR(self, i: int) -> float:
        i = max(1, min(i, len(self.r))); return float(self.r[i-1])
    def IRneeR(self, r: float) -> int:
        # 1-based index on a log10 grid like the Fortran mapping
        t = (math.log10(r) - self.rlgmin) * (len(self.r)-1) / (self.rlgmax - self.rlgmin) + 1.0
        return max(1, min(int(t), len(self.r)))

class _DensityGrid:
    def __init__(self, rgrid, vgrid):
        self.r = np.asarray(rgrid, float)
        self.v = np.asarray(vgrid, float)
    def RneeIR(self, i: int) -> float:
        i = max(1, min(i, len(self.r))); return float(self.r[i-1])
    def IRneeR(self, r: float) -> int:
        j = int(np.searchsorted(self.r, r, side="right")); return max(1, min(j, len(self.r)))
    def VneeIV(self, i: int) -> float:
        i = max(1, min(i, len(self.v))); return float(self.v[i-1])
    def IVneeV(self, v: float) -> int:
        a = abs(v); j = int(np.searchsorted(self.v, a, side="right")); return max(1, min(j, len(self.v)))

_RHALO: _HaloGrid | None = None
_DEN:   _DensityGrid | None = None

def register_halo_grid(rhalo: np.ndarray):
    global _RHALO, npot
    _RHALO = _HaloGrid(rhalo)
    npot = len(rhalo)

def register_density_grid(rgrid: np.ndarray, vgrid: np.ndarray):
    global _DEN, Nrdat, Nvdat
    _DEN = _DensityGrid(rgrid, vgrid)
    Nrdat, Nvdat = len(rgrid), len(vgrid)

# Fortran-style wrappers the rest of the code expects:
def RneeIRhalo(i: int) -> float:  return _RHALO.RneeIR(i)
def IRneeR(r: float) -> int:      return _DEN.IRneeR(r)
def RneeIR(i: int) -> float:      return _DEN.RneeIR(i)
def IVneeV(v: float) -> int:      return _DEN.IVneeV(v)
def VneeIV(i: int) -> float:      return _DEN.VneeIV(i)

# Legendre polynomial P_n (exactly what potential/force/denny need)
from numpy.polynomial.legendre import legval
def Pl(n: int, v: float) -> float:
    c = np.zeros(n + 1, dtype=float); c[n] = 1.0
    return float(legval(v, c))

# =================== potential/tables (READ THESE BEFORE USING) ===================
# Table/global placeholders — fill via tables_read()
# nlegup: highest even Legendre order; nlegs = nlegup//2 + 1
nlegup: int = 0
npot:   int = 0
rlgpotmin: float = 0.0
rlgpotmax: float = 0.0
hole: float = 0.0
gdennorm: float = 1.0

tabv:  np.ndarray | None = None  # shape (nlegs, npot)
tabfr: np.ndarray | None = None  # shape (nlegs, npot)
tabfv: np.ndarray | None = None  # shape (nlegs, npot)

##########################################################################################################################
# Placeholders

def projmass()
def halodens()
def haloread()
def halotables()
def halowrite()


# ---- Launch condition logic from oneorbit.f, simplified and explicit
def launch_from_peri_apo(xperi_arcsec: float, xapo_arcsec: float, v: float) -> tuple[State, float]:

# ---- Main orbit loop with energy/L2 diagnostics just like the Fortran
def run_single_orbit(xperi_arcsec: float, xapo_arcsec: float, v: float, Nsos2: int, Nskip: int, r_escape: float = 10.0):
    
if __name__ == "__main__":
    stats = run_single_orbit(xperi_arcsec=0.5, xapo_arcsec=5.0, v=0.2, Nsos2=50, Nskip=10)
    print(stats)

##########################################################################################################################
# =================== small utils ===================
def filter_weights(w: np.ndarray, tiny: float = 1e-37) -> np.ndarray:
    """Weed out negative or zero weights, replacing them with TINY."""
    return np.maximum(tiny, w)

__all__ = [
    # constants
    "GG","arcsec","distance","angrad","totlight","vmult","vadd","epsilon","efac","fac1",
    # grids
    "IRneeR","RneeIR","IVneeV","VneeIV","Tables","TAB",
    # types
    "State","IterParams","MLParams","PhaseReadResult",
    # tables/potential
    "tables_read","tables_write","potential","force","energy","make_ekin",
    # library helpers
    "build_cmatrix","new_cmatrix_update","summ","xlib_read","vaniread",
    "weight_read","weight_write","phasewrite","phaseread",
    # misc
    "load_qamatch","load_starlist","filter_weights",
]