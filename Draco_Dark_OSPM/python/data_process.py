import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# ===================== Constants =====================
RA0_DEG  = 260.0517
DEC0_DEG = 57.9153
DRACO_DIST_PC = 76000.0
D2R = np.pi/180.0

# OSPM working radius
R_MAX_PC = 500.0

# Elliptical-binning (tune as needed)
PA_DEG = 0.0          # position angle, east of north
AXIS_RATIO_Q = 0.70   # b/a, 0<q<=1
NBINS_R = 5           # radial bins for elliptical profiles
NBINS_THETA = 3     # angular sectors per annulus (1 = none)

# ===================== Inputs / Outputs =====================
data_file   = "Draco_stars_filtered_650pc.csv"
stars_out   = "draco_clean_data.csv"
prof_circ_out = "draco_profiles_circular.csv"
prof_ell_out  = "draco_profiles_elliptical.csv"
merged_out    = "draco_with_profiles.csv"

# ===================== Load =====================
df = pd.read_csv(data_file)

# --- Galactic l,b ---
coords = SkyCoord(ra=df["RA_deg"].values*u.deg, dec=df["Dec_deg"].values*u.deg, frame="icrs")
df["l_deg"] = coords.galactic.l.deg
df["b_deg"] = coords.galactic.b.deg

# --- Plane geometry (dx,dy,r) if missing ---
if "r_pc" not in df.columns:
    ra  = np.asarray(df["RA_deg"], float)*D2R
    dec = np.asarray(df["Dec_deg"], float)*D2R
    dra  = (ra  - RA0_DEG*D2R)*np.cos(DEC0_DEG*D2R)
    ddec =  dec - DEC0_DEG*D2R
    df["dx_pc"] = DRACO_DIST_PC * dra
    df["dy_pc"] = DRACO_DIST_PC * ddec
    df["r_pc"]  = np.hypot(df["dx_pc"], df["dy_pc"])

# --- Hard radial ceiling for OSPM ---
df = df[df["r_pc"] <= R_MAX_PC].copy()
print(f"Stars within {R_MAX_PC:.0f} pc: {len(df)}")

# --- Quality / membership cuts ---
m = np.isfinite(df.get("radial_velocity")) & np.isfinite(df.get("radial_velocity_error"))
if "ruwe" in df:        m &= (df["ruwe"] < 1.4)
if "pmra_error" in df:  m &= (df["pmra_error"]  < 0.3)
if "pmdec_error" in df: m &= (df["pmdec_error"] < 0.3)
if "radial_velocity_error" in df: m &= (df["radial_velocity_error"] < 5.0)

data = df.loc[m].copy()
if "r_pc" not in data:
    raise RuntimeError("r_pc missing after plane step.")

# --- LOS velocity relative to systemic ---
v_sys = np.nanmedian(data["radial_velocity"]) if "radial_velocity" in data else 0.0
data["vlos"] = data["radial_velocity"] - v_sys

def robust_sigma(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if x.size < 3: return np.nan
    med = np.median(x); mad = np.median(np.abs(x - med))
    return 1.4826*mad

def sigma_err_chi(s, n):
    if not np.isfinite(s) or n is None or n < 3:
        return np.nan
    return s / np.sqrt(2*(n-1))

# ===================== Elliptical annuli =====================
pa = np.deg2rad(PA_DEG)
c, s = np.cos(pa), np.sin(pa)

x_maj =  c*data["dx_pc"].to_numpy(float) + s*data["dy_pc"].to_numpy(float)
y_min = -s*data["dx_pc"].to_numpy(float) + c*data["dy_pc"].to_numpy(float)

R_ell = np.sqrt(x_maj**2 + (y_min/AXIS_RATIO_Q)**2)
theta = np.arctan2(y_min/AXIS_RATIO_Q, x_maj)   # [-pi, pi)

data["R_ell_pc"] = R_ell
data["theta_gal"] = theta
print("theta range (deg):", np.degrees(theta.min()), np.degrees(theta.max()))

# radial edges: equal-count in R_ell (fallback to linear if needed)
rad_edges = np.quantile(R_ell, np.linspace(0, 1, NBINS_R+1))
rad_edges = np.unique(rad_edges)
if len(rad_edges) - 1 < NBINS_R:
    rad_edges = np.linspace(np.nanmin(R_ell), np.nanmax(R_ell), NBINS_R+1)

r_centers_ell = 0.5*(rad_edges[1:]+rad_edges[:-1])
data["rbin_ell"] = pd.cut(R_ell, bins=rad_edges, include_lowest=True, labels=r_centers_ell)

# -------- Theta binning (wedge indices 0..NBINS_THETA-1) --------
if NBINS_THETA > 1:
    th = data["theta_gal"].to_numpy(float)
    th = np.mod(th + np.pi, 2*np.pi)  # map to [0, 2π)
    dth = 2*np.pi / NBINS_THETA
    sector = np.floor(th / dth).astype(int)
    sector = np.clip(sector, 0, NBINS_THETA-1)   # guard rare 2π edge
    data["tbin_ell"] = sector
    # optional: center angle (deg) for readability
    data["tbin_ell_center_deg"] = (sector + 0.5) * (360.0 / NBINS_THETA)
else:
    data["tbin_ell"] = 0
print("unique tbin_ell:", np.unique(data["tbin_ell"], return_counts=True))
grpcols_ell = ["rbin_ell", "tbin_ell"]
prof_ell = data.groupby(grpcols_ell, observed=False).agg(
    r_pc=("R_ell_pc", "median"),
    N=("vlos", "size"),
    sigma_vlos=("vlos", robust_sigma),
    pmra_med=("pmra", "median"),
    pmdec_med=("pmdec", "median"),
    sigma_pmra=("pmra", robust_sigma),
    sigma_pmdec=("pmdec", robust_sigma)
).reset_index()

prof_ell["sigma_vlos_err"] = [sigma_err_chi(sv, n) for sv, n in zip(prof_ell["sigma_vlos"], prof_ell["N"])]
prof_ell.to_csv(prof_ell_out, index=False)
print(f"Elliptical profiles → {prof_ell_out}")

# ===================== Equal-count circular annuli =====================
N_target = 20
nbins = max(1, int(np.ceil(len(data) / N_target)))
q = np.linspace(0.0, 1.0, nbins+1)
edges = np.quantile(data["r_pc"].to_numpy(), q)
edges = np.unique(edges)
if edges.size < 2:
    edges = np.array([data["r_pc"].min(), data["r_pc"].max()])

# keep the interval as a column to merge on
data["rbin_circ"] = pd.cut(data["r_pc"], bins=edges, include_lowest=True)
data["rbin_ell_id"]  = data["rbin_ell"].cat.codes
data["rbin_circ_id"] = data["rbin_circ"].cat.codes
prof_circ = data.groupby("rbin_circ", observed=False).agg(
    r_pc=("r_pc", "median"),
    N=("vlos", "size"),
    sigma_vlos=("vlos", robust_sigma),
    pmra_med=("pmra", "median"),
    pmdec_med=("pmdec", "median"),
    sigma_pmra=("pmra", robust_sigma),
    sigma_pmdec=("pmdec", robust_sigma)
).reset_index()

prof_circ["sigma_vlos_err"] = [sigma_err_chi(sv, n) for sv, n in zip(prof_circ["sigma_vlos"], prof_circ["N"])]
# prune tiny bins if desired
prof_circ = prof_circ.loc[prof_circ["N"] >= 8].reset_index(drop=True)
prof_circ.to_csv(prof_circ_out, index=False)
print(f"Circular profiles → {prof_circ_out}")

# ===================== Merge profiles back onto stars =====================
stars = data.copy()
stars["rbin_circ_key"] = stars["rbin_circ"].astype(str)

pc = prof_circ.copy()
pc["rbin_circ_key"] = pc["rbin_circ"].astype(str)
pc = pc.drop(columns=["rbin_circ"]).rename(columns={
    "r_pc":"r_pc_circ",
    "N":"N_circ",
    "sigma_vlos":"sigma_vlos_circ",
    "sigma_vlos_err":"sigma_vlos_err_circ",
    "pmra_med":"pmra_med_circ",
    "pmdec_med":"pmdec_med_circ",
    "sigma_pmra":"sigma_pmra_circ",
    "sigma_pmdec":"sigma_pmdec_circ",
})
stars = stars.merge(pc, on="rbin_circ_key", how="left")

# elliptical merge keys + rename columns with _ell suffix
stars["rbin_ell_key"] = stars["rbin_ell"].astype(str)
stars["tbin_ell_key"] = stars["tbin_ell"].astype(str)

pe = prof_ell.copy()
pe["rbin_ell_key"] = pe["rbin_ell"].astype(str)
pe["tbin_ell_key"] = pe["tbin_ell"].astype(str)
pe = pe.drop(columns=["rbin_ell", "tbin_ell"]).rename(columns={
    "r_pc":"r_pc_ell",
    "N":"N_ell",
    "sigma_vlos":"sigma_vlos_ell",
    "sigma_vlos_err":"sigma_vlos_err_ell",
    "pmra_med":"pmra_med_ell",
    "pmdec_med":"pmdec_med_ell",
    "sigma_pmra":"sigma_pmra_ell",
    "sigma_pmdec":"sigma_pmdec_ell",
})
print("pre-merge tbin_ell counts:", np.unique(stars["tbin_ell"], return_counts=True))
print("saving to:", merged_out)
print("post-merge tbin_ell counts:", np.unique(stars["tbin_ell"], return_counts=True))

stars = stars.merge(pe, on=["rbin_ell_key","tbin_ell_key"], how="left")

# ===================== Save =====================
stars.to_csv(merged_out, index=False)
data.to_csv(stars_out, index=False)
assert stars["tbin_ell"].notna().all(), "tbin_ell became NaN before save"
assert (stars["tbin_ell"]>=0).all(), "tbin_ell negative? unexpected"
print(f"Wrote cleaned star-level data → {stars_out}")
print(f"Wrote star table WITH circular+elliptical profile columns → {merged_out}")
print(prof_ell.pivot_table(index="rbin_ell", columns="tbin_ell", values="N", fill_value=0))

