# ============================================================================================================
# Orbit Super Position Modeling
# ============================================================================================================
from __future__ import annotations
import os
import json
from dataclasses import dataclass, fields
from typing import Callable, Optional, Sequence, Tuple, Union
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import biweight_location, biweight_scale
from pathlib import Path

# ============================================================================================================
# Data binding 
# ============================================================================================================
DRACO_INPUT_FILE: str = "Draco_stars_filtered_650pc.csv"
MAX_RADIUS_PC: float = 650.0
OUTPUT_CSV: str = f"Draco_stars_filtered_{int(MAX_RADIUS_PC)}pc.csv"
Draco_dat = pd.read_csv(DRACO_INPUT_FILE)

_alias = {
    "RA": "RA_deg", "DEC": "Dec_deg", "RA_ERR": "RA_error", "DEC_ERR": "Dec_error",
    "PAR": "parallax", "PAR_ERR": "parallax_error",
    "PMRA": "pmra", "PMRA_ERR": "pmra_error",
    "PMDEC": "pmdec", "PMDEC_ERR": "pmdec_error",
    "RV": "radial_velocity", "RV_ERR": "radial_velocity_error",
    "L_GAL": "l", "B_GAL": "b","G_MAG": "phot_g_mean_mag", 
    "BP_MAG": "phot_bp_mean_mag", "RP_MAG": "phot_rp_mean_mag",
    "BP_RP": "bp_rp", "RUWE": "ruwe", "DIST_GSPP": "distance_gspphot",
    "BP_RP_EXCESS": "phot_bp_rp_excess_factor",
    "dx_pc": "dx_pc", "dy_pc": "dy_pc", "r_pc": "r_pc",
}

def bind_columns(df: pd.DataFrame, strict: bool = False) -> None:
    g = globals()
    for alias, col in _alias.items():
        if col in df.columns:
            g[alias] = df[col]
        elif strict:
            raise KeyError(f"Missing column for alias {alias!r}: expected {col!r}")
bind_columns(Draco_dat, strict=False)
print("Galaxy Data Imported and Bound.")

# ============================================================================================================
# VARIABLES
# ============================================================================================================
nbins            = 18
theta_min_arcsec = 15.0
theta_max_deg    = 0.40
theta_min_arcsec = 5.0
theta_min_deg    = theta_min_arcsec / 3600.0
oversample       = 5
RR_cap_deg       = 0.50
trim             = 3 # drops the first few small-angle bins
K                = 8
# ============================================================================================================
# BINS
# ============================================================================================================
bins_deg = np.logspace(np.log10(theta_min_deg), np.log10(theta_max_deg), nbins + 1).astype(float)
bins_deg[0]  = theta_min_deg
bins_deg[-1] = theta_max_deg
for i in range(1, bins_deg.size):
    if not bins_deg[i] > bins_deg[i-1]:
        bins_deg[i] = np.nextafter(bins_deg[i-1], np.inf)

print("Initalized")

# ============================================================================================================
# Geometry & pair-counting
# ============================================================================================================
def vectorize(RA: np.ndarray, DEC: np.ndarray) -> np.ndarray:
    """Cartesian unit vectors on unit sphere."""
    RA = np.asarray(RA, dtype=float)
    DEC = np.asarray(DEC, dtype=float)
    if RA.shape != DEC.shape:
        raise ValueError(f"RA/DEC shape mismatch: {RA.shape} vs {DEC.shape}")
    ra = np.deg2rad(RA)
    dec = np.deg2rad(DEC)
    c = np.cos(dec)
    return np.c_[c * np.cos(ra), c * np.sin(ra), np.sin(dec)]

def _ensure_bins(bins_deg) -> np.ndarray:
    bins = np.asarray(bins_deg, dtype=float)
    if bins.ndim != 1 or bins.size < 2:
        raise ValueError("bins_deg must be a 1D array of at least 2 edges.")
    if not np.isfinite(bins).all():
        raise ValueError("bins_deg must be finite.")
    if bins[0] < 0.0 or bins[-1] > 180.0:
        raise ValueError("bins_deg must lie within [0, 180] degrees.")
    if not np.all(np.diff(bins) > 0):
        raise ValueError("bins_deg must be strictly increasing.")
    return bins

def _hist_angles_from_tree(vec_src, vec_cat, tree_cat, bins_deg, weights_left=None, weights_right=None, drop_self=False):
    bins = _ensure_bins(bins_deg)
    max_theta = np.deg2rad(bins[-1])
    r_chord = 2.0 * np.sin(max_theta / 2.0)
    counts = np.zeros(len(bins) - 1, dtype=float)

    wl = np.ones(len(vec_src), dtype=float) if weights_left  is None else np.asarray(weights_left,  dtype=float)
    wr = np.ones(len(vec_cat), dtype=float) if weights_right is None else np.asarray(weights_right, dtype=float)
    if wl.shape[0] != len(vec_src) or wr.shape[0] != len(vec_cat):
        raise ValueError("weights_left/weights_right must match src/cat lengths.")

    for i, v in enumerate(vec_src):
        idx = tree_cat.query_ball_point(v, r=r_chord)
        if not idx:
            continue
        if drop_self:
            idx = [j for j in idx if j != i]
            if not idx:
                continue
        dots = np.clip(vec_cat[idx] @ v, -1.0, 1.0)
        theta_deg = np.degrees(np.arccos(dots))
        bin_idx = np.searchsorted(bins, theta_deg, side="right") - 1
        valid = (bin_idx >= 0) & (bin_idx < len(counts))
        if not np.any(valid):
            continue
        w_pairs = wl[i] * wr[np.asarray(idx)[valid]]
        np.add.at(counts, bin_idx[valid], w_pairs)
    return counts

def mean_radec_wrap(ra_deg: np.ndarray, dec_deg: np.ndarray) -> Tuple[float, float]:
    m = np.isfinite(ra_deg) & np.isfinite(dec_deg)
    ra = np.deg2rad(ra_deg[m]); dec = np.deg2rad(dec_deg[m])
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    X, Y, Z = x.mean(), y.mean(), z.mean()
    r = np.hypot(np.hypot(X, Y), Z)
    if r == 0.0:
        raise RuntimeError("Mean direction undefined.")
    ra0 = (np.rad2deg(np.arctan2(Y, X)) + 360.0) % 360.0
    dec0 = np.rad2deg(np.arcsin(Z / r))
    return ra0, dec0

def _bin_pairs_from_sparse(rows, cols, dists, bin_edges_deg, wl, wr, auto=False):
    # chord distance d on unit sphere -> angle theta = 2*arcsin(d/2)
    theta = 2.0 * np.degrees(np.arcsin(0.5 * np.clip(dists, 0.0, 2.0)))
    if auto:
        keep = rows < cols
        rows, cols, theta = rows[keep], cols[keep], theta[keep]
    idx = np.searchsorted(bin_edges_deg, theta, side="right") - 1
    good = (idx >= 0) & (idx < len(bin_edges_deg) - 1)
    if not np.any(good):
        return np.zeros(len(bin_edges_deg) - 1, dtype=float)
    w_pairs = wl[rows[good]] * wr[cols[good]]
    counts = np.zeros(len(bin_edges_deg) - 1, dtype=float)
    np.add.at(counts, idx[good], w_pairs)
    return counts

def pair_counts(ra_deg=None, dec_deg=None, bins_deg=None, weights: Optional[np.ndarray] = None) -> np.ndarray:
    ra = _as_array(ra_deg, "RA"); dec = _as_array(dec_deg, "DEC")
    m = np.isfinite(ra) & np.isfinite(dec)
    vec = vectorize(ra[m], dec[m])
    wl  = np.ones(vec.shape[0], float) if weights is None else np.asarray(weights, float)[m]
    bins = _ensure_bins(bins_deg)
    r_chord = 2.0 * np.sin(np.deg2rad(bins[-1]) / 2.0)

    tree = cKDTree(vec)
    coo  = tree.sparse_distance_matrix(tree, max_distance=r_chord, output_type="coo_matrix")
    # self/duplicate pairs are removed inside _bin_pairs_from_sparse when auto=True
    return _bin_pairs_from_sparse(coo.row, coo.col, coo.data, bins, wl, wl, auto=True)


def pair_counts_cross(ra1_deg=None, dec1_deg=None, ra2_deg=None, dec2_deg=None, bins_deg=None,
                    weights_left: Optional[np.ndarray] = None, weights_right: Optional[np.ndarray] = None) -> np.ndarray:
    ra1 = _as_array(ra1_deg, "RA"); dec1 = _as_array(dec1_deg, "DEC")
    ra2 = _as_array(ra2_deg);        dec2 = _as_array(dec2_deg)
    m1  = np.isfinite(ra1) & np.isfinite(dec1)
    m2  = np.isfinite(ra2) & np.isfinite(dec2)

    v1 = vectorize(ra1[m1], dec1[m1]); v2 = vectorize(ra2[m2], dec2[m2])
    wL = np.ones(v1.shape[0], float) if weights_left  is None else np.asarray(weights_left,  float)[m1]
    wR = np.ones(v2.shape[0], float) if weights_right is None else np.asarray(weights_right, float)[m2]

    bins = _ensure_bins(bins_deg)
    r_chord = 2.0 * np.sin(np.deg2rad(bins[-1]) / 2.0)

    t1 = cKDTree(v1); t2 = cKDTree(v2)
    coo = t1.sparse_distance_matrix(t2, max_distance=r_chord, output_type="coo_matrix")
    return _bin_pairs_from_sparse(coo.row, coo.col, coo.data, bins, wL, wR, auto=False)
    
if "_require" not in globals():
    def _require(name: str):
        if name not in globals():
            avail = sorted([k for k, v in globals().items() if isinstance(v, (np.ndarray, pd.Series))])
            raise NameError(f"Expected global '{name}' not found. Available bound arrays: {avail}")
        return np.asarray(globals()[name], dtype=float)
    
def _as_array(x, fallback_name: Optional[str] = None):
    if x is None and fallback_name is not None:
        return _require(fallback_name)
    if isinstance(x, str):
        return _require(x)
    return np.asarray(x, dtype=float)

def make_randoms(ra0_deg: float, dec0_deg: float, radius_deg: float, n: int, seed: Optional[int] = None):
    rng = np.random.default_rng(seed)
    ra0 = np.deg2rad(ra0_deg); dec0 = np.deg2rad(dec0_deg)
    cospsi_min = np.cos(np.deg2rad(radius_deg))
    u = rng.uniform(cospsi_min, 1.0, size=n)   # cos(psi) uniform on [cos r, 1]
    psi = np.arccos(u)
    phi = rng.uniform(0.0, 2*np.pi, size=n)
    # local cap vector then rotate to (ra0, dec0)
    x = np.sin(psi) * np.cos(phi)
    y = np.sin(psi) * np.sin(phi)
    z = np.cos(psi)
    # rotate: first to dec0 about y, then to ra0 about z
    sx, cx = np.sin(dec0), np.cos(dec0)
    X =  cx * x + sx * z
    Y =  y
    Z = -sx * x + cx * z
    cRA = np.arctan2(Y, X) + ra0
    cDEC = np.arcsin(Z)
    RA_R  = (np.rad2deg(cRA) + 360.0) % 360.0
    DEC_R = np.rad2deg(cDEC)
    globals()["RA_R"] = RA_R
    globals()["DEC_R"] = DEC_R
    return RA_R, DEC_R

def wtheta_landy_szalay(raD=None, decD=None, raR=None, decR=None, bins_deg=None, wD=None, wR=None, return_counts=False):
    raD = _as_array(raD, "RA");   decD = _as_array(decD, "DEC")
    raR = _as_array(raR, "RA_R"); decR = _as_array(decR, "DEC_R")
    if bins_deg is None:
        bins_deg = np.linspace(0.0, 3.0, 51)
    mD = np.isfinite(raD) & np.isfinite(decD)
    mR = np.isfinite(raR) & np.isfinite(decR)
    ND, NR = int(mD.sum()), int(mR.sum())

    if ND < 2 or NR < 2:
        raise ValueError("Need at least 2 data and 2 random points for LS.")
        
    DD_raw = pair_counts(raD[mD], decD[mD], bins_deg, weights=(wD[mD] if wD is not None else None))
    DR_raw = pair_counts_cross(raD[mD], decD[mD], raR[mR], decR[mR], bins_deg, weights_left=(wD[mD] if wD is not None else None), 
            weights_right=(wR[mR] if wR is not None else None))
    RR_raw = pair_counts(raR[mR], decR[mR], bins_deg, weights=(wR[mR] if wR is not None else None))

    DD = DD_raw / (ND * (ND - 1) / 2.0)
    DR = DR_raw / (ND * NR)
    RR = RR_raw / (NR * (NR - 1) / 2.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        w = (DD - 2.0*DR + RR) / RR

    centers = 0.5 * (bins_deg[1:] + bins_deg[:-1])
    keep = np.isfinite(w) & (RR > 0)

    if not return_counts:
        return centers[keep], w[keep]

    out = {
        "DD_pairs": DD_raw[keep],
        "DR_pairs": DR_raw[keep],
        "RR_pairs": RR_raw[keep],
        "DD": DD[keep],
        "DR": DR[keep],
        "RR": RR[keep],
        "keep": keep,  # so you can align with bins if needed
    }
    return centers[keep], w[keep], out

def dedup_close(ra, dec, min_sep_arcsec=2.0):
    vec  = vectorize(ra, dec)
    tree = cKDTree(vec)
    rmin = 2.0*np.sin(np.deg2rad(min_sep_arcsec/3600.0)/2.0)
    pairs = tree.query_pairs(rmin)
    keep = np.ones(ra.size, bool)
    for i, j in pairs:
        keep[j] = False
    return ra[keep], dec[keep]

print("Functions Iitialized")

# ============================================================================================================
# ANALYSIS
# ============================================================================================================
RA      = np.asarray(RA,  dtype=float)
DEC     = np.asarray(DEC, dtype=float)
mD      = np.isfinite(RA) & np.isfinite(DEC)
raD     = RA[mD]
decD    = DEC[mD]

# de-duplicate close pairs
raD, decD = dedup_close(raD, decD, min_sep_arcsec=2.0)
ND = raD.size
if ND < 2:
    raise RuntimeError("Not enough finite RA/Dec in data for example run.")
# Draco center from data
ra0, dec0 = mean_radec_wrap(raD, decD)

# randoms
NR = max(3000, oversample * ND)
RA_R, DEC_R = make_randoms(ra0, dec0, RR_cap_deg, NR, seed=42)

# Landy–Szalay
centers, w = wtheta_landy_szalay(raD=raD, decD=decD, raR=RA_R, decR=DEC_R, bins_deg=bins_deg)
centers_plot, w_plot = centers[trim:], w[trim:]

def _sector_angles(ra, dec, ra0, dec0):
    # small-angle TAN projection about (ra0, dec0)
    dra  = (ra - ra0) * np.cos(np.deg2rad(dec0))
    ddec = (dec - dec0)
    ang  = (np.degrees(np.arctan2(ddec, dra)) + 360.0) % 360.0
    return ang

def LOC(x, axis=None):
    """Robust central value: biweight location; fallback to nanmedian."""
    try:
        return biweight_location(x, axis=axis)
    except Exception:
        return np.nanmedian(x, axis=axis)

def Scale(x, axis=None):
    """Robust scatter: biweight scale; fallback to 1.4826*MAD."""
    try:
        return biweight_scale(x, axis=axis)
    except Exception:
        med = np.nanmedian(x, axis=axis, keepdims=True)
        mad = np.nanmedian(np.abs(x - med), axis=axis)
        return 1.4826 * mad  # Gaussian MAD scale
    
def biwgt(arr):
    arr = np.asarray(arr, dtype=float)
    return LOC(arr), Scale(arr)

def jackknife_wtheta(raD: np.ndarray, decD: np.ndarray, bins_deg: np.ndarray, centers_eval: np.ndarray,
                oversample: int = 20, RR_cap_deg: float = 0.50, n_jack: int = 20, seed: int = 123,
):
    raD = np.asarray(raD, float)
    decD = np.asarray(decD, float)
    good = np.isfinite(raD) & np.isfinite(decD)
    raD, decD = raD[good], decD[good]
    ND = raD.size
    if ND < 3:
        raise ValueError("Not enough data for jackknife.")
    n_jack = int(min(max(2, n_jack), ND))  # clamp

    ra0, dec0 = mean_radec_wrap(raD, decD)
    NR = max(3000, oversample * ND)
    RA_R, DEC_R = make_randoms(ra0, dec0, RR_cap_deg, NR, seed=seed)

    order = np.argsort(raD)
    splits = np.array_split(order, n_jack)

    T = centers_eval.size
    Wk = np.zeros((n_jack, T), float)

    for k, idx_out in enumerate(splits):
        keep = np.ones(ND, bool); keep[idx_out] = False
        c_k, w_k = wtheta_landy_szalay(
            raD=raD[keep], decD=decD[keep],
            raR=RA_R,     decR=DEC_R,
            bins_deg=bins_deg
        )
        wk_interp = np.interp(centers_eval, c_k, w_k, left=w_k[0], right=w_k[-1]).astype(float)
        Wk[k] = wk_interp

    # jackknife mean & covariance
    w_mean = Wk.mean(axis=0)
    delta  = Wk - w_mean
    Cov    = ((n_jack - 1) / n_jack) * (delta.T @ delta)
    w_err  = np.sqrt(np.diag(Cov))
    return w_mean, w_err, Cov


# evaluate errors on your plotted grid (after any trimming you do)
w_loc, w_err, W = jackknife_wtheta(
    raD=raD, decD=decD,
    bins_deg=bins_deg,
    centers_eval=centers_plot,   # same grid you plot/save
    n_jack=K,
    oversample=oversample,
    RR_cap_deg=RR_cap_deg,
    seed=123
)

# optional: small robust rolling smoother (purely for visual guidance)
def rolling_biweight(y, win=3):
    n = len(y); out = np.full(n, np.nan, float); h = win // 2
    for i in range(n):
        i0 = max(0, i - h); i1 = min(n, i + h + 1)
        out[i] = LOC(y[i0:i1])
    return out

w_smooth = rolling_biweight(w_loc, win=3)


print("Analysis Complete")
print(f"ND={ND}, NR={NR}")
print("Centers [deg]:", centers_plot)
print("w(theta):", w_plot)

# ============================================================================================================
# Export 
# ============================================================================================================
out_dir = "Draco_Data"
os.makedirs(out_dir, exist_ok=True)

# Untrimmed (kept bins)
centers_kept, w_kept, counts = wtheta_landy_szalay(
    raD=raD, decD=decD, raR=RA_R, decR=DEC_R, bins_deg=bins_deg, return_counts=True
)
keep = counts["keep"]
theta_lo = bins_deg[:-1][keep]
theta_hi = bins_deg[1:][keep]

with np.errstate(divide="ignore", invalid="ignore"):
    w_err_pois = (1.0 + w_kept) / np.sqrt(np.maximum(counts["RR_pairs"], 1.0))
    w_err_pois[~np.isfinite(w_err_pois)] = np.nan

wtheta_df = pd.DataFrame({
    "theta_lo_deg":       theta_lo,
    "theta_hi_deg":       theta_hi,
    "theta_center_deg":   centers_kept,   # <-- use kept centers
    "w":                  w_kept,
    "w_err_poissonish":   w_err_pois,
    "DD_norm":            counts["DD"],
    "DR_norm":            counts["DR"],
    "RR_norm":            counts["RR"],
    "DD_pairs":           counts["DD_pairs"],
    "DR_pairs":           counts["DR_pairs"],
    "RR_pairs":           counts["RR_pairs"],
})

out_csv = f"{out_dir}/draco_wtheta_smallscales.csv"
wtheta_df.to_csv(out_csv, index=False)

# Jackknife (trimmed grid you evaluated)
jk_df = pd.DataFrame({
    "theta_center_deg": centers_plot,
    "w_trimmed":        w_plot,
    "w_jackknife":      w_loc,
    "w_jackknife_err":  w_err,
})
out_csv_j = f"{out_dir}/draco_wtheta_smallscales_jk.csv"
jk_df.to_csv(out_csv_j, index=False)

# Meta
meta = {
    "ND": int(ND),
    "NR": int(NR),
    "trim": int(trim),
    "nbins": int(nbins),
    "theta_min_arcsec": float(theta_min_arcsec),
    "theta_max_deg": float(theta_max_deg),
    "oversample": int(oversample),
    "RR_cap_deg": float(RR_cap_deg),
    "bins_deg": [float(x) for x in bins_deg.tolist()],
    "keep_idx": np.nonzero(keep)[0].tolist(),
}
out_meta = f"{out_dir}/draco_wtheta_smallscales_meta.json"
with open(out_meta, "w") as f:
    json.dump(meta, f, indent=2)

w_all = pd.read_csv("Draco_Data/draco_wtheta_smallscales.csv")
w_jk  = pd.read_csv("Draco_Data/draco_wtheta_smallscales_jk.csv")

print(f"Saved: {out_csv}")
print(f"Saved: {out_csv_j}")
print(f"Saved: {out_meta}")
print(f"Untrimmed rows: {len(wtheta_df)} | JK rows: {len(jk_df)}")

with open("Draco_Data/draco_wtheta_smallscales_meta.json") as f:
    meta = json.load(f)
# ============================================================================================================
# PLOTTING
# ============================================================================================================
plt.figure(figsize=(6,4))
# untrimmed with Poisson-ish
plt.errorbar(w_all["theta_center_deg"], w_all["w"], yerr=w_all["w_err_poissonish"],
            fmt="o", lw=1, capsize=2, label="Poisson-ish")
# trimmed jackknife
plt.errorbar(w_jk["theta_center_deg"], w_jk["w_jackknife"], yerr=w_jk["w_jackknife_err"],
            fmt="s-", lw=1, capsize=2, label="Jackknife")
plt.axhline(0, ls="--", lw=0.8)
plt.xscale("log")
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$w(\theta)$")
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================================================================
# N2I
# New to implement
# ============================================================================================================

def get_sky(spectra_stack: np.ndarray, method: str = "biweight") -> Tuple[np.ndarray, np.ndarray]:

    if method == "biweight":
        sky = LOC(spectra_stack, axis=0)
        scat = Scale(spectra_stack, axis=0)
    else:
        sky = np.nanmedian(spectra_stack, axis=0)
        scat = 1.4826 * np.nanmedian(np.abs(spectra_stack - sky[None, :]), axis=0)
    return sky, scat


def get_1sigma(spec_path: str | Path, wave: float, nflimh: int = 3, write_to: str | Path = "get1sig.out") -> float:

    df = pd.read_csv(spec_path, sep=r"\s+", header=None, engine="python")
    n = len(df)
    if n < 10:
        f1sig = 0.0
    else:
        w = df.iloc[:, 0].to_numpy(dtype=float)
        x3 = df.iloc[:, 2].to_numpy(dtype=float)
        x9 = df.iloc[:, 8].to_numpy(dtype=float)
        se = np.where(x9 > 0.0, x3 / x9, 0.0)
        iw0 = int(np.argmin(np.abs(w - float(wave))))
        ilo = max(0, iw0 - int(nflimh))
        ihi = min(n - 1, iw0 + int(nflimh))
        f1sig = float(np.sqrt(np.sum(se[ilo:ihi + 1] ** 2)))
    with open(write_to, "w") as fh:
        fh.write(f"{f1sig}\n")
    return f1sig

def get_angle():
    radtodeg = 57.29578
    xd = 1.7
    yd = xd

    out_rows = []
    path = "coords.in"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing input file: {path}")

    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 18:
                # Need a1 + 17 numbers minimum
                continue

            a1 = parts[0]  # unused but preserved for parity
            try:
                vals = list(map(float, parts[1:18]))
            except ValueError:
                # Skip malformed numeric lines
                continue

            x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17 = vals

            if abs(x2) > xd or abs(x3) > yd:
                continue

            cosd = np.cos(x10 / radtodeg)
            di = np.sqrt(x16 * x16 + x17 * x17)
            rdiff = 3600.0 * cosd * (x14 - x9)
            ddiff = 3600.0 * (x15 - x10)
            ds = np.sqrt(rdiff * rdiff + ddiff * ddiff)

            if di == 0.0 or ds == 0.0:
                # Mirror Fortran's undefined division with a quiet NaN rather than crashing
                ati = np.nan
                ats = np.nan
            else:
                ati = radtodeg * np.arctan2(x17 / di, x16 / di)
                ats = radtodeg * np.arctan2(rdiff / ds, ddiff / ds)

            if np.isfinite(ati) and ati < 0.0:
                ati += 360.0
            if np.isfinite(ats) and ats < 0.0:
                ats += 360.0

            ang = (ats if np.isfinite(ats) else 0.0) + 180.0 + (ati if np.isfinite(ati) else 0.0)
            if ang > 360.0:
                ang -= 360.0
            if ang < 0.0:
                ang += 360.0

            ang0 = 270.0 + x12
            if ang0 > 360.0:
                ang0 -= 360.0
            if ang0 < 0.0:
                ang0 += 360.0

            xdum = abs(ang0 - x8)
            if 80.0 < xdum < 270.0:
                angf = x8 + 180.0
                if angf > 360.0:
                    angf -= 360.0
                if angf < 0.0:
                    angf += 360.0
            else:
                angf = x8

            print(ati, ats, ang0, angf, x8)
            out_rows.append((ati, ats, ang0, angf, x8))

    return pd.DataFrame(out_rows, columns=["ati", "ats", "ang0", "angf", "x8"])

'''
# ==========================================
# Translations
# ==========================================

def get_agn():
    """
    Python equivalent of the provided Fortran routine.
    Reads:
    - 'good.spec' for template spectrum (wg, sg)
    - 'in' for header (x1 x2 x3 i4 rad) then data rows (x, y, yn0)
    Produces:
    - 'out2' with sbest, abest
    - 'out3' with wg, y_interpolated(x->wg), stot, sgp, ap
    - 'out4' with x, (y - ap_interp)/xbnorm over full x
    Also plots the raw spectrum and the fitted components.
    Returns (sbest, abest, diff, dict_of_arrays)
    """
    nmax = 10000
    wmin = 22600.0
    wmax = 22800.0
    w1 = 22500.0
    w2 = 24000.0

    idmax = 3

    # Read first line of 'in' for rad, then all subsequent lines into x,y,yn0
    if not os.path.exists("in"):
        raise FileNotFoundError("Missing input file 'in'")
    x_list, y_list, yn0_list = [], [], []
    with open("in", "r") as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError("'in' is empty")
    hdr = lines[0].split()
    if len(hdr) < 5:
        raise ValueError("First line of 'in' must contain at least 5 values: x1 x2 x3 i4 rad")
    # Preserve names (x1,x2,x3,i4,rad) but only 'rad' is used here
    x1, x2, x3, i4, rad = float(hdr[0]), float(hdr[1]), float(hdr[2]), int(float(hdr[3])), float(hdr[4])
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) < 3:
            continue
        xi, yi, yi0 = float(parts[0]), float(parts[1]), float(parts[2])
        x_list.append(xi); y_list.append(yi); yn0_list.append(yi0)
        if len(x_list) >= nmax:
            break
    if not x_list:
        raise ValueError("No data rows found in 'in' after header")
    x = np.asarray(x_list, dtype=float)
    y = np.asarray(y_list, dtype=float)
    yn0 = np.asarray(yn0_list, dtype=float)

    # Compute plot bounds over 22000..24200 as in Fortran
    xmin = 22000.0
    xmax = 24200.0
    mask_plot = (x > xmin) & (x < xmax)
    if np.any(mask_plot):
        ymin = float(np.min(y[mask_plot]))
        ymax = float(np.max(y[mask_plot]))
    else:
        ymin = 1e10
        ymax = 0.0

    # Read good.spec (wg, sg) and build agna scaled/normalized by agn0 slope
    if not os.path.exists("good.spec"):
        raise FileNotFoundError("Missing input file 'good.spec'")
    wg_list, sg_list = [], []
    with open("good.spec", "r") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            wg_list.append(float(parts[0])); sg_list.append(float(parts[1]))
            if len(wg_list) >= nmax:
                break
    if not wg_list:
        raise ValueError("No data found in 'good.spec'")
    wg = np.asarray(wg_list, dtype=float)
    sg = np.asarray(sg_list, dtype=float)

    wave0 = 22000.0
    agn0 = 1250.0
    wave1 = 24000.0
    agn1 = 1750.0
    if rad > 15.0:
        agn1 = 1500.0
    slope = (agn1 - agn0) / (wave1 - wave0)
    agna = agn0 + (wg - wave0) * slope
    agna = agna / agn0

    # Optional smoothing branch (ism == 0 in original, so off by default)
    ism = 0
    if ism == 1:
        ibin = 5
        ib1 = (ibin - 1) // 2
        nbb = 0
        xn = []
        yn = []
        for j in range(0, len(x), ibin):
            nbb += 1
            istart = max(0, j - ib1)
            iend = istart + ibin - 1
            if iend >= len(x):
                iend = len(x) - 1
                istart = max(0, iend - ibin + 1)
            seg = slice(istart, iend + 1)
            yb = y[seg]
            xb = x[seg]
            y_loc, y_scl = biwgt(yb)
            x_loc, x_scl = biwgt(xb)
            yn.append(y_loc); xn.append(x_loc)
        xn = np.asarray(xn); yn = np.asarray(yn)
        # Rebuild y with smoothed values via linear interpolation
        y_sm = np.array([xlinint(xi, xn, yn) for xi in x], dtype=float)
        y = y_sm
        # Recompute plot bounds
        mask_plot = (x > xmin) & (x < xmax)
        if np.any(mask_plot):
            ymin = float(np.min(y[mask_plot])); ymax = float(np.max(y[mask_plot]))
        else:
            ymin = 1e10; ymax = 0.0

    # Plot original spectrum
    plt.figure()
    plt.plot(x, y)
    plt.xlim(xmin, xmax)
    if ymin < ymax:
        plt.ylim(ymin, ymax)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Input spectrum and fitted components")

    # Grid search with iterative refinement (idmax passes)
    ns = 60
    smin = 30.0
    smax = 80.0
    if rad < 5.0:
        smax = 160.0
    na = 50
    amin = 0.0
    amax = 1500.0

    diff = 1e30
    sbest = np.nan
    abest = np.nan
    idone = 0

    def eval_rms(stry, atry):
        # Accumulate squared residuals over fitting window
        mask = (x > xmin) & (x < xmax)
        if not np.any(mask):
            return np.inf
        # Interpolate sg and agna to x[mask]
        s0 = np.array([xlinint(xi, wg, sg) for xi in x[mask]], dtype=float) * stry
        a0 = np.array([xlinint(xi, wg, agna) for xi in x[mask]], dtype=float) * atry
        res2 = (s0 + a0 - y[mask]) ** 2
        xb0, xs0 = biwgt(res2)
        if not np.isfinite(xb0) or xb0 < 0.0:
            return np.inf
        return float(np.sqrt(xb0))

    for _ in range(idmax):
        diff_pass = diff
        for ia in range(na):
            atry = amin + float(ia) * (amax - amin) / float(max(na - 1, 1))
            for i in range(ns):
                stry = smin + float(i) * (smax - smin) / float(max(ns - 1, 1))
                rms = eval_rms(stry, atry)
                if rms < diff:
                    diff = rms
                    sbest = stry
                    abest = atry
        print(idone, sbest, abest, diff)
        idone += 1
        if idone >= idmax:
            break
        # Refine search window around current best
        ns = 60
        smin = sbest * 0.8
        smax = sbest * 1.2
        na = 50
        amin = abest * 0.8
        amax = abest * 1.2
        if abest < 20.0:
            na = 20
            amin = 0.0
            amax = 20.0
            if abest < 9.0:
                amax = 10.0

    # Build fitted components on wg grid
    sgp = sg * sbest
    ap = agna * abest
    stot = sgp + ap

    # Overlay fitted components
    plt.plot(wg, sgp)
    plt.plot(wg, ap)
    plt.plot(wg, stot)
    plt.tight_layout()

    # Write out2: best-fit scales
    with open("out2", "w") as f2:
        f2.write(f"{sbest} {abest}\n")

    # Write out3: wg, y(x->wg), stot, sgp, ap
    with open("out3", "w") as f3:
        for wi, st, ss, aa in zip(wg, stot, sgp, ap):
            y0 = xlinint(wi, x, y)
            f3.write(f"{wi} {y0} {st} {ss} {aa}\n")

    # Write out4: continuum-subtracted and normalized by biweight in [wmin,wmax]
    # Subtract AGN component interpolated to x
    a_on_x = np.array([xlinint(xi, wg, ap) for xi in x], dtype=float)
    y_sub = y - a_on_x
    sel_norm = (x > wmin) & (x < wmax)
    xin = y_sub[sel_norm]
    if xin.size == 0:
        xbnorm = 1.0
    else:
        xbnorm, xsnorm = biwgt(xin)
        if not np.isfinite(xbnorm) or xbnorm == 0.0:
            xbnorm = 1.0
    with open("out4", "w") as f4:
        for xi, yi in zip(x, y_sub / xbnorm):
            f4.write(f"{xi} {yi}\n")

    return sbest, abest, diff, {"x": x, "y": y, "wg": wg, "sgp": sgp, "ap": ap, "stot": stot}

def ata():
    """
    Python equivalent of getata.f.
    Reads filenames from 'inlist', each containing two columns (x1, x2).
    For each row index i, computes the biweight location (xb) and scale (xs)
    across files of the second column, then writes lines to 'ata.out' as:
    w(i) xb xs ntot
    """
    nmax = 10000

    # Use existing biweight if present; otherwise fall back to robust helpers without redefining globals.
    if "biwgt" in globals():
        _biwgt = biwgt
    else:
        def _biwgt(arr):
            return LOC(arr), Scale(arr)

    if not os.path.exists("inlist"):
        raise FileNotFoundError("Missing input file 'inlist'")

    files = []
    with open("inlist", "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            files.append(s)
            if len(files) >= nmax:
                break
    if not files:
        raise ValueError("'inlist' contained no usable filenames")

    waves_list = []
    wo_list = []
    for file1 in files:
        if not os.path.exists(file1):
            raise FileNotFoundError(f"Listed file not found: {file1}")
        w_i = []
        o_i = []
        with open(file1, "r") as f2:
            for ln in f2:
                t = ln.strip()
                if not t or t.startswith("#"):
                    continue
                parts = t.split()
                if len(parts) < 2:
                    continue
                w_i.append(float(parts[0]))
                o_i.append(float(parts[1]))
                if len(w_i) >= nmax:
                    break
        if not w_i:
            raise ValueError(f"No numeric data in '{file1}'")
        waves_list.append(np.asarray(w_i, dtype=float))
        wo_list.append(np.asarray(o_i, dtype=float))

    ntot = len(wo_list)
    w = waves_list[-1]
    n = w.shape[0]

    for idx, arr in enumerate(wo_list):
        if arr.shape[0] != n:
            raise ValueError(f"Length mismatch: '{files[idx]}' has {arr.shape[0]} rows; expected {n}")

    wo_mat = np.vstack(wo_list)  # shape: (ntot, n)

    xb_list = []
    xs_list = []
    # Create new file; error if it already exists to mirror Fortran 'status=new'
    with open("ata.out", "x") as fout:
        for i in range(n):
            xin_col = wo_mat[:, i]
            xb, xs = _biwgt(xin_col)
            xb_list.append(float(xb))
            xs_list.append(float(xs))
            fout.write(f"{w[i]} {xb} {xs} {ntot}\n")

    return pd.DataFrame({"w": w, "xb": np.array(xb_list), "xs": np.array(xs_list), "ntot": ntot})

# =======================
# getdisp / getdispb / getmatch / getsky
# =======================

def get_dispersion(df: pd.DataFrame, col_name: str) -> float:
    """
    Rough dispersion estimator (Å/pix) or (km/s/pix) from a wavelength/pixel column.
    For a vector 'col_name' that should be monotonic in pixel index, return
    robust mean of first-differences.
    """
    w = np.asarray(df[col_name].values, float)
    dw = np.diff(w)
    return float(LOC(dw))

def get_dispersion_b(df: pd.DataFrame, col_name: str) -> Tuple[float, float]:
    """Same as get_dispersion, but also returns a robust scatter."""
    w = np.asarray(df[col_name].values, float)
    dw = np.diff(w)
    return float(LOC(dw)), float(Scale(dw))

def get_match(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    ra1: str, dec1: str,
    ra2: str, dec2: str,
    max_sep_arcsec: float,
    extra_key_1: Optional[str] = None,
    extra_key_2: Optional[str] = None,
    max_abs_diff_extra: Optional[float] = None,
) -> pd.DataFrame:
    """
    Cross-match df1 to df2 on sky within max_sep_arcsec.
    Optionally require |df1[extra_key_1]-df2[extra_key_2]| <= max_abs_diff_extra.
    Returns a dataframe with indices and separations.
    Mirrors getmatch.f logic, but on the sphere.
    """
    c1 = SkyCoord(df1[ra1].values*u.deg, df1[dec1].values*u.deg)
    c2 = SkyCoord(df2[ra2].values*u.deg, df2[dec2].values*u.deg)
    idx2, sep2d, _ = c1.match_to_catalog_sky(c2)
    ok = sep2d.arcsec <= max_sep_arcsec
    if extra_key_1 and extra_key_2 and max_abs_diff_extra is not None:
        ok &= np.abs(df1[extra_key_1].values - df2.iloc[idx2][extra_key_2].values) <= max_abs_diff_extra
    out = pd.DataFrame({
        "idx1": np.arange(len(df1))[ok],
        "idx2": idx2[ok],
        "sep_arcsec": sep2d.arcsec[ok]
    })
    return out


# =======================
# gethist* wrappers
# =======================
def _hist_angles_from_tree(vec_src, vec_cat, tree_cat, bins_deg, weights_left=None, drop_self=False):
    max_theta = np.deg2rad(bins_deg[-1])
    r_chord = 2.0 * np.sin(max_theta/2.0)
    counts = np.zeros(len(bins_deg)-1, dtype=float)
    for i, v in enumerate(vec_src):
        idx = tree_cat.query_ball_point(v, r=r_chord)
        if not idx:
            continue
        if drop_self:
            idx = [j for j in idx if j != i]
            if not idx:
                continue
        dots = np.clip(vec_cat[idx] @ v, -1.0, 1.0)
        theta_deg = np.degrees(np.arccos(dots))
        h, _ = np.histogram(theta_deg, bins=bins_deg)
        counts += h if weights_left is None else h * float(weights_left[i])
    return counts

def get_hist(df, col_name, bins=50, range=None, density=False):
    x = np.asarray(df[col_name].values, float)
    hist, edges = np.histogram(x, bins=bins, range=range, density=density)
    return hist, edges

get_hist0   = get_hist
get_hist2   = get_hist
get_hist_cx = get_hist

def gethist_hdr2():
    """Placeholder for any HDR2-specific binning you had; keep signature freeform."""
    raise NotImplementedError("Implement HDR2-specific histogram here.")

# =======================
# getfield / getfield3 / getfielda
# =======================

def _select_by_radius(df: pd.DataFrame, ra_col: str, dec_col: str,
                    ra0_deg: float, dec0_deg: float, radius_arcsec: float) -> pd.DataFrame:
    c = SkyCoord(df[ra_col].values*u.deg, df[dec_col].values*u.deg)
    c0 = SkyCoord(ra0_deg*u.deg, dec0_deg*u.deg)
    sep = c.separation(c0).arcsec
    return df[sep <= radius_arcsec].copy()

def get_field(df, ra_col, dec_col, center_ra_deg, center_dec_deg, radius_arcsec):
    """Simple radial selection."""
    return _select_by_radius(df, ra_col, dec_col, center_ra_deg, center_dec_deg, radius_arcsec)

def get_field_3(df, ra_col, dec_col, centers_deg: Sequence[Tuple[float,float]], radius_arcsec):
    """Union of 3 (or N) circular fields."""
    mask = np.zeros(len(df), dtype=bool)
    for (ra0, dec0) in centers_deg:
        mask |= _select_by_radius(df, ra_col, dec_col, ra0, dec0, radius_arcsec).index.isin(df.index)
    return df[mask].copy()

def get_field_c(df, ra_col, dec_col, poly_vertices_deg: Sequence[Tuple[float,float]]):
    """Convex polygon cut (approx; project to tangent plane at polygon center)."""
    verts = np.asarray(poly_vertices_deg)
    ra0, dec0 = verts.mean(axis=0)
    c = SkyCoord(df[ra_col].values*u.deg, df[dec_col].values*u.deg)
    c0 = SkyCoord(ra0*u.deg, dec0*u.deg)
    # TAN projection small-angle approx:
    dx = (c.ra - c0.ra).to(u.deg).value * np.cos(np.deg2rad(dec0))
    dy = (c.dec - c0.dec).to(u.deg).value
    from matplotlib.path import Path
    poly = Path((verts[:,0]-ra0, (verts[:,1]-dec0))).transformed(
        np.array([[np.cos(np.deg2rad(dec0)), 0],[0, 1]])  # same approx
    )
    inside = poly.contains_points(np.vstack([dx, dy]).T)
    return df[inside].copy()

# =======================
# getnorm / getnormexp / getgaus / getgpinfo
# =======================

def getnorm(wave: np.ndarray, flux: np.ndarray, t_curve1: np.ndarray, t_curve2: np.ndarray) -> np.ndarray:
    """
    Equivalent of getnorm.f: normalize 'flux' by two smooth curves (interpolated to wave).
    Returns normalized flux.
    """
    from numpy.interp import interp
    y1 = interp(wave, t_curve1[:,0], t_curve1[:,1])
    y2 = interp(wave, t_curve2[:,0], t_curve2[:,1])
    return flux / (y1*y2)

def getnormexp(exp_name: str, calib_table_path: str) -> Tuple[float, Tuple[float,float,float]]:
    """
    Port of getnormexp.f logic: read per-exposure FWHM and 3-channel normalizations.
    Returns (fwhm, (n1,n2,n3)). You’ll adapt I/O to your actual calibration store.
    """
    # Stub: implement your table lookup here
    raise NotImplementedError

def getgaus(x, amp, mu, sigma):
    return amp * np.exp(-0.5*((x-mu)/sigma)**2) / (np.sqrt(2*np.pi)*sigma)

def getgpinfo():
    """Placeholder for ‘gp info’ getter."""
    raise NotImplementedError

# =======================
# getoff* — offsets / registration
# =======================

@dataclass
class OffsetResult:
    dx_arcsec: float
    dy_arcsec: float
    rot_deg: float
    n_pairs: int
    dx_err_arcsec: float
    dy_err_arcsec: float

def get_off(cat1: pd.DataFrame, cat2: pd.DataFrame,
            ra1: str, dec1: str, ra2: str, dec2: str,
            match_arcsec: float = 2.0) -> OffsetResult:
    """
    Rigid translation/rotation fit between two source lists (like getoff/getoff2/getoff3 lineage).
    """
    m = get_match(cat1, cat2, ra1, dec1, ra2, dec2, max_sep_arcsec=match_arcsec)
    if m.empty:
        return OffsetResult(0,0,0,0,np.inf,np.inf)
    p = cat1.iloc[m.idx1][[ra1,dec1]].to_numpy()
    q = cat2.iloc[m.idx2][[ra2,dec2]].to_numpy()
    # Project to tangent plane at the mean:
    ra0, dec0 = p.mean(axis=0)
    c0 = SkyCoord(ra0*u.deg, dec0*u.deg)
    def proj(arr):
        c = SkyCoord(arr[:,0]*u.deg, arr[:,1]*u.deg)
        dx = (c.ra - c0.ra).to(u.arcsec).value * np.cos(np.deg2rad(dec0))
        dy = (c.dec - c0.dec).to(u.arcsec).value
        return np.c_[dx,dy]
    P, Q = proj(p), proj(q)
    # Solve similarity transform Q ≈ R P + t (rotation + translation, no scale)
    Pm = P - P.mean(0); Qm = Q - Q.mean(0)
    U,S,Vt = np.linalg.svd(Pm.T @ Qm)
    R = U @ Vt
    if np.linalg.det(R) < 0:  # enforce proper rotation
        U[:,-1] *= -1
        R = U @ Vt
    t = Q.mean(0) - (R @ P.mean(0))
    rot = np.degrees(np.arctan2(R[1,0], R[0,0]))
    resid = Q - (P @ R.T + t)
    sx, sy = Scale(resid[:,0]), Scale(resid[:,1])
    return OffsetResult(float(t[0]), float(t[1]), float(rot), len(P), float(sx), float(sy))

# convenience aliases
get_off2 = get_off
get_off3 = get_off

# =======================
# getsdssg / getsignif / getsignifa / getstat2
# =======================

def getsdssg(wave: np.ndarray, flux: np.ndarray, sdss_g_curve: np.ndarray) -> float:
    """
    Synthetic SDSS-g flux: integrate flux * T(λ) dλ / ∫T dλ. 'sdss_g_curve' = Nx2 [λ, T].
    """
    lam = sdss_g_curve[:,0]
    T = sdss_g_curve[:,1]
    T_interp = np.interp(wave, lam, T, left=0, right=0)
    num = np.trapz(flux * T_interp, wave)
    den = np.trapz(T_interp, wave)
    return float(num/den) if den > 0 else np.nan

def getsignif(wave: np.ndarray, flux: np.ndarray, err: np.ndarray,
            sigma_pix: float = 2.0) -> pd.DataFrame:
    """
    Scan for Gaussian line significance across spectrum.
    Returns DataFrame with columns: wave0, amp, sn, cont, chi2_red (mimics getsignif.f outputs).
    """
    # simple matched filter with Gaussian kernel
    x = wave
    dx = np.median(np.diff(x))
    halfwin = int(np.ceil(4*sigma_pix))
    kx = np.arange(-halfwin, halfwin+1)
    g = np.exp(-0.5*(kx/sigma_pix)**2)
    g /= np.sum(g)
    # estimate continuum by rolling median
    from scipy.ndimage import median_filter
    cont = median_filter(flux, size=max(9, 2*halfwin+1))
    resid = flux - cont
    # matched filter amplitude and S/N
    amp = np.convolve(resid, g, mode="same")
    noise = np.sqrt(np.convolve(err**2, g**2, mode="same"))
    sn = np.divide(amp, noise, out=np.zeros_like(amp), where=noise>0)
    # pseudo chi2
    chi2 = ((resid/np.maximum(err, 1e-9))**2)
    chi2 = median_filter(chi2, size=max(9, 2*halfwin+1))
    df = pd.DataFrame({"wave0": x, "amp": amp, "sn": sn, "cont": cont, "chi2_red": chi2})
    return df

def getsignifa(signif_tables: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """
    Aggregate many getsignif outputs to robust per-λ summaries (biweight disp & 95th pct),
    like your getsignifa.f.
    """
    # stack by interpolation to common λ grid:
    lam = signif_tables[0]["wave0"].to_numpy()
    stacks = []
    for df in signif_tables:
        if not np.allclose(df["wave0"].to_numpy(), lam):
            # interpolate sn onto common grid
            sn_i = np.interp(lam, df["wave0"].to_numpy(), df["sn"].to_numpy(), left=0, right=0)
        else:
            sn_i = df["sn"].to_numpy()
        stacks.append(sn_i)
    A = np.vstack(stacks)
    disp = Scale(A, axis=0)
    # 95th percentile per λ
    q95 = np.quantile(A, 0.95, axis=0)
    return pd.DataFrame({"wave0": lam, "sn_disp": disp, "sn_95": q95})

def getstat2(df1: pd.DataFrame, df2: pd.DataFrame, sn_col: str, chi_col: str,
            sn1: float = 3.0, sn2: float = 4.0, chi_cut: float = 1.4) -> dict:
    """
    Port of getstat2.f: compare counts above two S/N thresholds after a χ^2 cut,
    and flag if one file is an outlier.
    """
    sel1 = df1[chi_col] < chi_cut
    sel2 = df2[chi_col] < chi_cut
    n1_sn1 = int((df1[sel1][sn_col] > sn1).sum())
    n1_sn2 = int((df1[sel1][sn_col] > sn2).sum())
    n2_sn1 = int((df2[sel2][sn_col] > sn1).sum())
    n2_sn2 = int((df2[sel2][sn_col] > sn2).sum())
    # expected counts in df1 scaled from df2
    x1 = n2_sn1/len(df2) * len(df1)
    x2 = n2_sn2/len(df2) * len(df1)
    igood = int((n1_sn1 > 2*x1) and (n1_sn1 >= 2))
    return dict(n1_sn1=n1_sn1, n1_sn2=n1_sn2, x1_exp=x1, x2_exp=x2, igood=igood)

# =======================
# getvol* — unify into one
# =======================

def compute_effective_volume(
    L_log10: np.ndarray,
    sn_cut_list: Sequence[float],
    z_grid: np.ndarray,
    vol_of_z: Callable[[np.ndarray], np.ndarray],
    FtoL_of_z: Callable[[np.ndarray], np.ndarray],
    completeness: Callable[[float, float, float], float],
    f50_per_field: Sequence[float],
) -> pd.DataFrame:
    """
    Unified port of getvol, getvol0, getvol2, getvol2_mc, getvol3:
    For each (L, S/N cut), integrate over z to get the relative V_eff given
    instrument completeness, flux limit variations (f50 per field), and the F↔L conversion.

    completeness(f, z, sn_cut) should return ∈[0,1].
    """
    L_log10 = np.atleast_1d(L_log10)
    sn_cut_list = list(sn_cut_list)
    out = []
    for snc in sn_cut_list:
        for L in L_log10:
            Llin = 10**(L-40.0)  # convert to 1e-40 cgs units as per your FtoL
            Fconv = FtoL_of_z(z_grid)  # flux per luminosity at z
            F_needed = Llin / Fconv * 1e17  # to match your units
            # average completeness across fields via f50 scaling:
            c_sum = 0.0
            for f50 in f50_per_field:
                # you can bake the sn-dependent rescale into completeness
                c_sum += completeness(F_needed, z_grid, snc)  # vectorized
            c_avg = c_sum / max(1, len(f50_per_field))
            dV = vol_of_z(z_grid)  # differential or cumulative—adapt to your table
            veff = np.trapz(c_avg * dV, z_grid)
            out.append((L, snc, veff))
    df = pd.DataFrame(out, columns=["log10L", "sn_cut", "Veff"])
    return df

# =======================
# gettheta / gettheta_dr / gettheta_w
# =======================

def pair_counts(ra_deg: np.ndarray, dec_deg: np.ndarray,
                bins_deg: np.ndarray, weights: Optional[np.ndarray]=None) -> np.ndarray:
    """
    Fast auto pair counts DD using BallTree with haversine metric.
    bins_deg are θ edges in degrees.
    """
    coords = np.deg2rad(np.c_[dec_deg, ra_deg])
    tree = BallTree(coords, metric='haversine')
    # haversine expects radians; convert deg bins to radians/arc lengths on unit sphere
    bins_rad = np.deg2rad(bins_deg)
    counts = np.zeros(len(bins_deg)-1, dtype=float)
    # approximate counts via two-point: sum over points the histogram of neighbor distances
    for i in range(coords.shape[0]):
        dists, ind = tree.query_radius(coords[i:i+1], r=bins_rad[-1], return_distance=True, sort_results=True)
        if len(ind[0]) == 0: continue
        # exclude self (distance ~0)
        d = dists[0][1:] if len(dists[0])>0 else np.array([])
        if d.size == 0: continue
        if weights is None:
            w = np.ones_like(d)
        else:
            # weight by weight of the central point (or average—mirror of getwtheta_w)
            w = np.full_like(d, weights[i], dtype=float)
        h, _ = np.histogram(d, bins=bins_rad)
        counts += h * w
    return counts

# =======================
# getxtr / getwavo
# =======================

def getwavo(wave_maps: Sequence[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Biweight mean & scatter wavelength solution across many maps (1032x112).
    """
    cube = np.stack(wave_maps, axis=0)  # (N, 1032, 112)
    loc = LOC(cube, axis=0)
    scat = Scale(cube, axis=0)
    return loc, scat

def detrend_rows_by_colbands(images: Sequence[np.ndarray], bands: int = 15,ref_index: int = 0,
                            smooth_box: int = 9) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Approx of getxtr.f:
    - Use first image as reference baseline.
    - For each image, in ~‘bands’ column chunks, compute residual vs ref,
        average across the chunk for each row, smooth along rows, subtract trend.
    - Return (detrended_images, robust_mean, robust_scatter).
    """
    imgs = np.stack(images, axis=0)  # (N, Y, X) ~ (N,112,1032) or (N,1032,112) depending on your order
    if imgs.shape[1] == 1032 and imgs.shape[2] == 112:
        imgs = np.transpose(imgs, (0,2,1))  # (N,112,1032)
    N, H, W = imgs.shape
    ref = imgs[ref_index]
    band_edges = np.linspace(0, W, bands+1, dtype=int)
    from scipy.ndimage import uniform_filter1d
    detr = imgs.copy()
    for n in range(N):
        resid = imgs[n] - ref
        trend = np.zeros_like(resid)
        for b0, b1 in zip(band_edges[:-1], band_edges[1:]):
            chunk = resid[:, b0:b1]
            row_curve = np.nanmean(chunk, axis=1)  # per-row mean residual
            row_curve = uniform_filter1d(row_curve, size=smooth_box, mode='nearest')
            trend[:, b0:b1] = row_curve[:, None]
        detr[n] = imgs[n] - trend
    mean = LOC(detr, axis=0)
    scat = Scale(detr, axis=0)
    return detr, mean, scat

# =======================
# getwtr — visualization
# =======================

def show_rowcol_residuals(images: Sequence[np.ndarray], row: int = 37, col: int = 1020, ref_index: int = 0):
    """
    Emulates getwtr.f quicklook: plot residuals along a fixed row and column vs a reference frame.
    """
    imgs = np.stack(images, axis=0)
    if imgs.shape[1] == 1032 and imgs.shape[2] == 112:
        imgs = np.transpose(imgs, (0,2,1))  # (N,112,1032)
    ref = imgs[ref_index]
    plt.figure(figsize=(10,6))
    # row
    plt.subplot(2,1,1)
    for i, im in enumerate(imgs):
        plt.plot(im[row,:]-ref[row,:], alpha=0.5, lw=1)
    plt.ylabel(f"Δ (row {row})")
    # col
    plt.subplot(2,1,2)
    for i, im in enumerate(imgs):
        plt.plot(im[:,col]-ref[:,col], alpha=0.5, lw=1)
    plt.ylabel(f"Δ (col {col})"); plt.xlabel("Pixel")
    plt.tight_layout()
    plt.show()
'''