#Grav/Optical Center
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import astropy as ap 
from dataclasses import dataclass
from matplotlib.patches import Ellipse
import os

# =======================
# CONSTANTS
# =======================

pi = np.pi
RA0_DEG  = 260.0517 #Draco optical center RA in degrees
DEC0_DEG = 57.9153 #Draco optical center Dec in degrees
DRACO_DIST_PC = 76000.0 # Draco optical center (J2000) and adopted system distance
PAR_SNR_MIN = 5.0            # require p / σp >= 5 to use parallax
CLIP_DISTANCE_PC = (20000, 120000)
arcsec_per_rad = 206265.0
D2R, R2R = pi/180.0, 180.0/pi

# =======================
# FILES
# =======================
# Column names present in Draco_stars.csv
FILE_DRACO          = "Draco_stars.csv"
RA_COL              = "RA_deg"
DEC_COL             = "Dec_deg"
RA_ERR_COL          = "RA_error"              # mas
DEC_ERR_COL         = "Dec_error"             # mas
PAR_COL             = "parallax"              # mas
PAR_ERR_COL         = "parallax_error"        # mas
PMRA_COL            = "pmra"
PMRA_ERR_COL        = "pmra_error"
PMDEC_COL           = "pmdec"
PMDEC_ERR_COL       = "pmdec_error"
RV_COL              = "radial_velocity"
RV_ERR_COL          = "radial_velocity_error"
L_COL               = "l"
B_COL               = "b"
G_MAG_COL           = "phot_g_mean_mag"
BP_MAG_COL          = "phot_bp_mean_mag"
RP_MAG_COL          = "phot_rp_mean_mag"
BP_RP_COL           = "bp_rp"
RUWE_COL            = "ruwe"
DIST_GSPP_COL       = "distance_gspphot"
BP_RP_EXCESS_COL    = "phot_bp_rp_excess_factor"
WEIGHT_MODE         = "equal"
ERR_UNITS           = "mas"

# =======================
# CONFIGURATION
# =======================
USE_GSPPHOT_DISTANCE = False         # use per-star distances
USE_PARALLAX = False                 # To ignore per-star distances set false
GSPP_FRAC_ERR = 0.30                 # assume 30% fractional LOS uncertainty for distance_gspphot
CLIP_DISTANCE_PC = (50000, 110000)   # clip distances to [50, 110] kpc to avoid crazy z
PAR_SNR_MIN = 5.0                    # p/sigma_p threshold if USE_PARALLAX=True
TARGET_TOL_PC = 0.01   
EXPORT_FILTERED_CSV = True
MAX_RADIUS_PC = 650.0
OUTPUT_CSV = f"Draco_stars_filtered_{int(MAX_RADIUS_PC)}pc.csv"

# =======================
# Helpers
# =======================
def load_df(path):
    df = pd.read_csv(path, sep=None, engine="python")
    need = {RA_COL, DEC_COL, RA_ERR_COL, DEC_ERR_COL, DIST_GSPP_COL, PAR_COL, PAR_ERR_COL}
    missing = {RA_COL, DEC_COL, RA_ERR_COL, DEC_ERR_COL} - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    # it's fine if parallax or distance_gspphot are NaN; we handle that
    return df

# ---- weighting ----
def choose_weights(df: pd.DataFrame) -> np.ndarray:
    return np.ones(len(df), float)  # equal weights are appropriate for the GLS solve


def distances_with_sigma(df):
    # safe fallbacks if someone forgets the config block
    use_gspp   = globals().get("USE_GSPPHOT_DISTANCE", True)
    gspp_ferr  = globals().get("GSPP_FRAC_ERR", 0.30)
    clip_range = globals().get("CLIP_DISTANCE_PC", (50000, 110000))
    use_par    = globals().get("USE_PARALLAX", False)
    par_snrmin = globals().get("PAR_SNR_MIN", 5.0)

    n = len(df)
    D  = np.full(n, np.nan, float)
    sD = np.full(n, np.nan, float)

    if use_gspp and DIST_GSPP_COL in df.columns:
        D = np.asarray(df[DIST_GSPP_COL], float)
        if clip_range is not None:
            lo, hi = clip_range
            D = np.where(np.isfinite(D), np.clip(D, lo, hi), D)
        sD = np.where(np.isfinite(D), gspp_ferr * D, np.nan)

    if use_par and PAR_COL in df.columns and PAR_ERR_COL in df.columns:
        p  = np.asarray(df[PAR_COL], float)
        pe = np.asarray(df[PAR_ERR_COL], float)
        snr = np.divide(p, pe, out=np.zeros_like(p), where=np.isfinite(p)&np.isfinite(pe)&(pe>0))
        good = (p > 0) & (snr >= par_snrmin)
        D_par  = np.full_like(p, np.nan, dtype=float)
        sD_par = np.full_like(p, np.nan, dtype=float)
        D_par[good]  = 1000.0 / p[good]
        sD_par[good] = 1000.0 * pe[good] / (p[good]**2)
        use = np.isnan(D) & np.isfinite(D_par)
        if clip_range is not None:
            lo, hi = clip_range
            D_par = np.where(np.isfinite(D_par), np.clip(D_par, lo, hi), D_par)
        D[use]  = D_par[use]
        sD[use] = sD_par[use]

    return D, sD

def build_positions_and_cov(df):
    ra  = np.asarray(df[RA_COL], float)*D2R
    dec = np.asarray(df[DEC_COL], float)*D2R
    ra0, dec0 = RA0_DEG*D2R, DEC0_DEG*D2R

    D, sD = distances_with_sigma(df)
    D_use = np.where(np.isfinite(D), D, DRACO_DIST_PC)
    Dref  = DRACO_DIST_PC

    dra, ddec = ra - ra0, dec - dec0
    x = D_use * dra  * np.cos(dec0)     # pc
    y = D_use * ddec                    # pc
    z = D_use - Dref                    # pc

    # RA/Dec uncertainties are in mas → radians
    s_ra_rad  = np.asarray(df[RA_ERR_COL],  float) / 1000.0 / arcsec_per_rad
    s_dec_rad = np.asarray(df[DEC_ERR_COL], float) / 1000.0 / arcsec_per_rad

    Sxx = (D_use*np.cos(dec0))**2 * (s_ra_rad**2)
    Syy = (D_use**2)            * (s_dec_rad**2)
    Sxy = np.zeros_like(Sxx)                       # no RA–Dec corr column in your header
    Szz = np.where(np.isfinite(sD), sD**2, 1e12)   # huge if LOS uncertainty unknown → uninformative

    return x,y,z,Sxx,Syy,Sxy,Szz,Dref

def add_plane_geometry(df: pd.DataFrame) -> pd.DataFrame:
    ra  = pd.to_numeric(df[RA_COL],  errors="coerce").to_numpy(float)*D2R
    dec = pd.to_numeric(df[DEC_COL], errors="coerce").to_numpy(float)*D2R
    dra  = (ra  - RA0_DEG*D2R)*np.cos(DEC0_DEG*D2R)
    ddec =  dec - DEC0_DEG*D2R
    dx_pc = DRACO_DIST_PC * dra
    dy_pc = DRACO_DIST_PC * ddec
    out = df.copy()
    out["dx_pc"] = dx_pc
    out["dy_pc"] = dy_pc
    out["r_pc"]  = np.hypot(dx_pc, dy_pc)
    return out

def build_plane_cov(df: pd.DataFrame):
    # RA/Dec uncertainties are mas → radians
    s_ra = pd.to_numeric(df[RA_ERR_COL],  errors="coerce").to_numpy(float) / 1000.0 / arcsec_per_rad
    s_de = pd.to_numeric(df[DEC_ERR_COL], errors="coerce").to_numpy(float) / 1000.0 / arcsec_per_rad
    fac_x = (DRACO_DIST_PC*np.cos(DEC0_DEG*D2R))**2
    fac_y = (DRACO_DIST_PC)**2
    Sxx = fac_x * (s_ra**2)
    Syy = fac_y * (s_de**2)
    Sxy = np.zeros_like(Sxx)  # no RA–Dec corr in your table
    return Sxx, Syy, Sxy

# --- 2D GLS on the tangent plane (ignore z completely) ---
def gls_center_2d(x, y, Sxx, Syy, Sxy):
    P = np.zeros((2,2)); b = np.zeros(2)
    for xi, yi, sxx, syy, sxy in zip(x, y, Sxx, Syy, Sxy):
        if not np.isfinite(sxx) or not np.isfinite(syy) or sxx<=0 or syy<=0:
            continue
        C = np.array([[sxx, sxy],[sxy, syy]], float)
        C[0,0] += 1e-12; C[1,1] += 1e-12
        W = np.linalg.inv(C)
        P += W
        b += W @ np.array([xi, yi], float)
    mu = np.linalg.solve(P, b)
    Cov = np.linalg.inv(P)
    return mu, Cov

def pc_to_deg(R_pc, D_pc):  # handy for ADQL cones
    return (R_pc / D_pc) * (180.0/np.pi)

# --- 1σ covariance ellipse helper ---
def plot_cov_ellipse(ax, Cov, mean, nsig=1.0, **kwargs):
    vals, vecs = np.linalg.eigh(Cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    width, height = 2*nsig*np.sqrt(vals)
    angle = np.degrees(np.arctan2(vecs[1,0], vecs[0,0]))
    ax.add_patch(Ellipse(xy=mean, width=width, height=height, angle=angle, fill=False, lw=2, **kwargs))

def main():
    df = pd.read_csv(FILE_DRACO, sep=None, engine="python")
    need = {RA_COL, DEC_COL, RA_ERR_COL, DEC_ERR_COL}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = add_plane_geometry(df)
    keep = np.isfinite(df["r_pc"]) & (df["r_pc"] <= MAX_RADIUS_PC)

    if EXPORT_FILTERED_CSV:
        out_path = f"Draco_stars_filtered_{int(MAX_RADIUS_PC)}pc.csv"
        df.loc[keep].to_csv(out_path, index=False)
        print(f"Wrote {keep.sum()}/{len(df)} rows to {out_path}")

    if not np.any(keep):
        raise RuntimeError("No stars left after radius cut; increase MAX_RADIUS_PC.")

    x = df.loc[keep, "dx_pc"].to_numpy(float)
    y = df.loc[keep, "dy_pc"].to_numpy(float)
    Sxx, Syy, Sxy = build_plane_cov(df.loc[keep])

    (mu2, Cov2) = gls_center_2d(x, y, Sxx, Syy, Sxy)
    xc, yc = float(mu2[0]), float(mu2[1])
    sigx, sigy = np.sqrt(np.diag(Cov2))

    # use your constants as named: R2R and arcsec_per_rad
    ra_c  = RA0_DEG + (xc / DRACO_DIST_PC) * R2R / np.cos(DEC0_DEG*D2R)
    dec_c = DEC0_DEG + (yc / DRACO_DIST_PC) * R2R
    sep_pc = float(np.hypot(xc, yc))
    sep_as = float((sep_pc / DRACO_DIST_PC) * arcsec_per_rad)
    print(f"2D GLS center: Δx={xc:.2f}±{sigx:.2f} pc, Δy={yc:.2f}±{sigy:.2f} pc")
    print(f"Offset = {sep_pc:.2f} pc ({sep_as:.2f}\") → RA={ra_c:.6f}°, Dec={dec_c:.6f}°")

    lim = MAX_RADIUS_PC
    fig, ax = plt.subplots(figsize=(7,6))
    ax.scatter(x, y, s=6, alpha=0.6)
    ax.scatter(0, 0, c='k', s=60, marker='x', label='Optical center')
    ax.scatter(xc, yc, c='r', s=70, marker='^', label='2D GLS center')
    plot_cov_ellipse(ax, Cov2, (xc, yc), nsig=1.0, color='r')
    ax.add_artist(plt.Circle((0,0), MAX_RADIUS_PC, fill=False, ls='--', lw=1.2, color='gray'))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('ΔX [pc]'); ax.set_ylabel('ΔY [pc]')
    ax.legend(); plt.tight_layout(); plt.show()

if __name__ == "__main__":
    main()