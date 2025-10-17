import numpy as np
import pandas as pd
import astropy.units as u
from pathlib import Path
from datetime import datetime
from astropy.coordinates import SkyCoord
import os, glob

# -------------------- paths / run stamp --------------------
optical_file = "draco_with_profiles.csv"   # Gaia DR3 (2016 epoch)
spec_file    = "draco_spec_all.asc"        # 2016-08-04
outdir = Path("products"); outdir.mkdir(exist_ok=True)
stamp = datetime.now().strftime("%Y%m%d")

print(f"PY: {__file__ if '__file__' in globals() else '<stdin>'}")
print(f"CWD: {os.getcwd()}")

# -------------------- load optical (Gaia) --------------------
optical_cols = [
    'source_id','RA_deg','RA_error','Dec_deg','Dec_error','parallax','parallax_error',
    'pmra','pmra_error','pmdec','pmdec_error','radial_velocity','radial_velocity_error',
    'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','bp_rp','ruwe',
    'distance_gspphot','phot_bp_rp_excess_factor',
    'dx_pc','dy_pc','r_pc','l_deg','b_deg','vlos','R_ell_pc','theta_gal',
    'rbin_ell','tbin_ell','tbin_ell_center_deg','rbin_circ','rbin_ell_id','rbin_circ_id',
    'rbin_circ_key','r_pc_circ','N_circ','sigma_vlos_circ','pmra_med_circ','pmdec_med_circ',
    'sigma_pmra_circ','sigma_pmdec_circ','sigma_vlos_err_circ',
    'rbin_ell_key','tbin_ell_key','r_pc_ell','N_ell','sigma_vlos_ell','pmra_med_ell',
    'pmdec_med_ell','sigma_pmra_ell','sigma_pmdec_ell','sigma_vlos_err_ell'
]

def load_optical(path):
    hdr = pd.read_csv(path, nrows=0)
    use = [c for c in optical_cols if c in hdr.columns]
    opt = pd.read_csv(path, usecols=use, dtype={"source_id":"string"})
    for c in use:
        if c != "source_id":
            opt[c] = pd.to_numeric(opt[c], errors="coerce")
    # enforce unique string IDs (avoid Excel-style rounding)
    opt['source_id'] = opt['source_id'].astype('string')
    opt = opt.drop_duplicates(subset='source_id', keep='first').reset_index(drop=True)
    return opt

opt = load_optical(optical_file)
print(f"[OPT] rows={len(opt)} uniq_ids={opt['source_id'].nunique()} RA≈{opt['RA_deg'].median():.3f} Dec≈{opt['Dec_deg'].median():.3f}")

# -------------------- load spectra (ASCII) --------------------
def load_spec_ascii(path):
    raw = pd.read_csv(path, sep=r"\s+", comment="#", header=None, skiprows=1,
                    engine="python").dropna(how="all")
    ncols = raw.shape[1]
    if ncols not in (24, 25):
        raise ValueError(f"Unexpected column count {ncols}; expected 24 or 25")

    arr = raw.values
    out = pd.DataFrame(index=raw.index)
    out['RAh'] = pd.to_numeric(arr[:,0], errors='coerce')
    out['RAm'] = pd.to_numeric(arr[:,1], errors='coerce')
    out['RAs'] = pd.to_numeric(arr[:,2], errors='coerce')

    if ncols == 24:
        out['DEd_signed'] = pd.to_numeric(arr[:,3], errors='coerce')
        out['DEm']        = pd.to_numeric(arr[:,4], errors='coerce')
        out['DEs']        = pd.to_numeric(arr[:,5], errors='coerce')
        out['HJD']        = pd.to_numeric(arr[:,6], errors='coerce')
        t0 = 7
    else:  # 25 columns: explicit sign column
        sgn = pd.Series(arr[:,3]).map({"+":1.0,"-":-1.0,"0":1.0,"1":-1.0}).fillna(1.0).astype(float).to_numpy()
        deg = pd.to_numeric(arr[:,4], errors='coerce')
        out['DEd_signed'] = sgn * deg
        out['DEm']        = pd.to_numeric(arr[:,5], errors='coerce')
        out['DEs']        = pd.to_numeric(arr[:,6], errors='coerce')
        out['HJD']        = pd.to_numeric(arr[:,7], errors='coerce')
        t0 = 8

    tail = ['S/N','Vlos','e_Vlos','s_Vlos','k_Vlos','Teff','e_Teff','s_Teff','k_Teff',
            'logg','e_logg','s_logg','k_logg','[Fe/H]','e_[Fe/H]','s_[Fe/H]','k_[Fe/H]']
    for j, name in enumerate(tail, start=t0):
        out[name] = pd.to_numeric(arr[:,j], errors='coerce')

    # degrees
    rah = out['RAh'].fillna(0.0); ram = out['RAm'].fillna(0.0); ras = out['RAs'].fillna(0.0)
    out['RA_deg'] = 15.0*(rah + ram/60.0 + ras/3600.0)
    des = out['DEs'].copy()
    bad60 = des >= 60.0
    if bad60.any():
        carry = np.floor(des[bad60]/60.0)
        out.loc[bad60,'DEm'] += carry
        des.loc[bad60] -= 60.0*carry
    ded = out['DEd_signed'].fillna(0.0); dem = out['DEm'].fillna(0.0)
    out['Dec_deg'] = ded + np.sign(np.where(ded==0.0, 1.0, ded))*(dem/60.0 + des.fillna(0.0)/3600.0)

    return out[['RA_deg','Dec_deg','HJD'] + tail]

spec = load_spec_ascii(spec_file)
print(f"[SPEC] rows={len(spec)} RA_med={np.nanmedian(spec['RA_deg']):.3f} Dec_med={np.nanmedian(spec['Dec_deg']):.3f}")

# ------------------------- crossmatch + tiny boresight shift -------------------------
def _match_once(opt_d, sp_d, cut_arcsec):
    # guard against NaNs before SkyCoord
    oo = opt_d[['RA_deg','Dec_deg']].apply(pd.to_numeric, errors='coerce')
    so = sp_d[['RA_deg','Dec_deg']].apply(pd.to_numeric, errors='coerce')
    mo = oo['RA_deg'].notna() & oo['Dec_deg'].notna()
    ms = so['RA_deg'].notna() & so['Dec_deg'].notna()
    opt_ok = opt_d.loc[mo]
    spec_ok = sp_d.loc[ms]
    if len(opt_ok) == 0 or len(spec_ok) == 0:
        return pd.DataFrame(), np.array([])

    c_o = SkyCoord(opt_ok['RA_deg'].to_numpy()*u.deg,  opt_ok['Dec_deg'].to_numpy()*u.deg)
    c_s = SkyCoord(spec_ok['RA_deg'].to_numpy()*u.deg, spec_ok['Dec_deg'].to_numpy()*u.deg)
    idx, sep2d, _ = c_s.match_to_catalog_sky(c_o)
    sep = sep2d.to(u.arcsec).value
    m = sep <= cut_arcsec

    a = spec_ok.loc[m].reset_index(drop=True).copy()
    b = opt_ok.iloc[idx[m]].reset_index(drop=True).copy()
    a['sep_arcsec'] = sep[m]
    return pd.concat([a, b.add_prefix('opt_')], axis=1), sep

def _crossmatch(opt_df, spec_df, cut_arcsec=12.0, seed_cut_arcsec=15.0, k_closest=30, estimate_shift=True):
    seed, _ = _match_once(opt_df, spec_df, seed_cut_arcsec)
    base_tight = int((seed['sep_arcsec'] <= cut_arcsec).sum()) if len(seed) else 0

    if estimate_shift and len(seed):
        g = seed.nsmallest(min(k_closest, len(seed)), 'sep_arcsec')
        if len(g):
            dec_rad = np.deg2rad(g['Dec_deg'].to_numpy())
            dx = (g['opt_RA_deg'].to_numpy()  - g['RA_deg'].to_numpy()) * np.cos(dec_rad) * 3600.0
            dy = (g['opt_Dec_deg'].to_numpy() - g['Dec_deg'].to_numpy()) * 3600.0
            dx_med, dy_med = float(np.median(dx)), float(np.median(dy))

            sp2 = spec_df.copy()
            cosd = np.cos(np.deg2rad(sp2['Dec_deg'].to_numpy()))
            cosd[cosd == 0] = 1.0
            sp2['RA_deg']  = sp2['RA_deg']  + (dx_med/3600.0)/cosd
            sp2['Dec_deg'] = sp2['Dec_deg'] +  (dy_med/3600.0)

            tight, _ = _match_once(opt_df, sp2, cut_arcsec)
            if len(tight) >= base_tight:
                print(f"[XM] applied shift dx,dy = {dx_med:.2f}\", {dy_med:.2f}\"")
                return tight, sp2

    tight, _ = _match_once(opt_df, spec_df, cut_arcsec)
    return tight, spec_df

matched_all, spec_aligned = _crossmatch(opt, spec, cut_arcsec=12.0, estimate_shift=True)
matched_best = (matched_all.sort_values('sep_arcsec', kind='mergesort')
                            .drop_duplicates(subset='opt_source_id', keep='first')
                            .reset_index(drop=True))

# enforce ID as string everywhere post-match
for df in (matched_all, matched_best):
    if df is not None and len(df):
        df['opt_source_id'] = df['opt_source_id'].astype('string')

# ------------------ curate per-exposure errors, prune 5σ Gaia conflicts ------------------
if 'Vlos' not in matched_all.columns or len(matched_all) == 0:
    raise RuntimeError("No matched spectra found (matched_all is empty).")

m = np.isfinite(pd.to_numeric(matched_all['Vlos'], errors='coerce'))

e_v = matched_all['e_Vlos'] if 'e_Vlos' in matched_all.columns else pd.Series(np.nan, index=matched_all.index)
s_v = matched_all['s_Vlos'] if 's_Vlos' in matched_all.columns else pd.Series(np.nan, index=matched_all.index)
e_use = np.where(np.isfinite(e_v), e_v, s_v)
e_use = np.where(np.isfinite(e_use) & (e_use > 0), e_use, 10.0)   # floor
e_use = np.clip(e_use, 3.0, 100.0)                                # tame extremes

cur = matched_all.loc[m].copy()
cur['e_use'] = e_use[m]

if {'opt_radial_velocity','opt_radial_velocity_error'}.issubset(cur.columns):
    dv = cur['Vlos'] - cur['opt_radial_velocity']
    dv_err = np.sqrt(cur['e_use']**2 + cur['opt_radial_velocity_error'].fillna(0.0)**2)
    keep = np.isfinite(dv_err) & (dv_err > 0) & (np.abs(dv/dv_err) <= 5)
    cur = cur.loc[keep].copy()

# ------------------------ collapse to one velocity per Gaia star ------------------------
def _collapse_per_star_using(df, e_col='e_use'):
    rows = []
    for sid, g in df.groupby('opt_source_id', dropna=False):
        v = g['Vlos'].to_numpy(float)
        e = g[e_col].to_numpy(float)
        ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
        if ok.sum() == 0:
            continue
        w = 1.0/(e[ok]**2)
        rows.append({
            'opt_source_id': sid,
            'Vlos': float(np.average(v[ok], weights=w)),
            'e_Vlos': float(1.0/np.sqrt(w.sum())),
            'opt_rbin_ell_key': g.get('opt_rbin_ell_key', pd.Series([np.nan])).iloc[0],
            'opt_tbin_ell_key': g.get('opt_tbin_ell_key', pd.Series([np.nan])).iloc[0],
            'opt_r_pc_ell':     g.get('opt_r_pc_ell',     pd.Series([np.nan])).iloc[0],
            'opt_rbin_circ_key':g.get('opt_rbin_circ_key',pd.Series([np.nan])).iloc[0],
            'opt_r_pc_circ':    g.get('opt_r_pc_circ',    pd.Series([np.nan])).iloc[0],
        })
    return pd.DataFrame(rows)

stars_spec = _collapse_per_star_using(cur, e_col='e_use')
print(f"[collapse] stars_with_spec={len(stars_spec)}")

# ----------------------------- optional area moments -----------------------------
def _bin_moments(df, mode='ell'):
    if mode=='ell' and {'opt_rbin_ell_key','opt_tbin_ell_key','opt_r_pc_ell'}.issubset(df.columns):
        groups, rcol = ['opt_rbin_ell_key','opt_tbin_ell_key'], 'opt_r_pc_ell'
    elif mode=='circ' and {'opt_rbin_circ_key','opt_r_pc_circ'}.issubset(df.columns):
        groups, rcol = ['opt_rbin_circ_key'], 'opt_r_pc_circ'
    else:
        return pd.DataFrame()
    rows = []
    for _, g in df.groupby(groups, dropna=False):
        v = g['Vlos'].to_numpy(float); e = g['e_Vlos'].to_numpy(float)
        ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
        if ok.sum() < 2:
            continue
        w = 1.0/(e[ok]**2); wsum = w.sum()
        mu = np.average(v[ok], weights=w)
        var_obs  = np.average((v[ok]-mu)**2, weights=w)
        var_meas = np.average(e[ok]**2,      weights=w)
        sig2 = max(0.0, var_obs - var_meas)
        sig  = np.sqrt(sig2)
        Neff = (wsum**2)/(w**2).sum()
        sig_err = np.sqrt(sig2/(2.0*max(Neff-1.0,1.0))) if sig2>0 else np.nan
        row = {'bin_r_key': g[groups[0]].iloc[0],
                'R_pc': float(np.nanmedian(g[rcol])),
                'N_spec': int(ok.sum()),
                'vlos_mean': float(mu),
                'sigma_vlos': float(sig),
                'sigma_vlos_err': float(sig_err)}
        if len(groups) == 2:
            row['bin_t_key'] = g[groups[1]].iloc[0]
        rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty: 
        return out
    return out.sort_values([c for c in ['R_pc','bin_t_key'] if c in out.columns]).reset_index(drop=True)

bins_ell  = _bin_moments(stars_spec, mode='ell')
bins_circ = _bin_moments(stars_spec, mode='circ')

# -------------------------- OSPM starlist (μ free offset) --------------------------
opt['source_id'] = opt['source_id'].astype('string')
stars_spec['opt_source_id'] = stars_spec['opt_source_id'].astype('string')
opt_full = opt.rename(columns={'source_id':'opt_source_id'})
starlist = stars_spec.merge(opt_full, on='opt_source_id', how='left')

n0 = len(starlist)
starlist = starlist.dropna(subset=['RA_deg','Dec_deg']).copy()
if len(starlist) < n0:
    print(f"[warn] dropped {n0-len(starlist)} rows without RA/Dec after merge")
starlist['Vlos']   = pd.to_numeric(starlist['Vlos'],  errors='coerce')
starlist['e_Vlos'] = pd.to_numeric(starlist['e_Vlos'], errors='coerce')
starlist = starlist[np.isfinite(starlist['Vlos']) & np.isfinite(starlist['e_Vlos']) & (starlist['e_Vlos'] > 0)]

cols_star = ['opt_source_id','RA_deg','Dec_deg','Vlos','e_Vlos','opt_r_pc_ell','opt_rbin_ell_key','opt_tbin_ell_key']
starlist[cols_star].to_csv(outdir / f"draco_ospm_starlist_mu_{stamp}.csv", index=False)
print(f"[OSPM] starlist → {outdir / f'draco_ospm_starlist_mu_{stamp}.csv'}  rows={len(starlist)}")

# ----------------------------- QA export (rich per-exposure) -----------------------------
def write_csv(df, stem, cols=None):
    if df is None or len(df) == 0:
        print(f"[skip] {stem} (empty)"); return
    use = df if cols is None else df[[c for c in cols if c in df.columns]]
    path = outdir / f"{stem}_{stamp}.csv"
    use.to_csv(path, index=False, float_format="%.6f")
    print(f"[ok] {stem} → {path}  rows={len(use)}  cols={use.shape[1]}")

qa_cols  = [
    'opt_source_id','RA_deg','Dec_deg','HJD','S/N',
    'Vlos','e_Vlos','s_Vlos',
    'Teff','e_Teff','s_Teff','k_Teff',
    'logg','e_logg','s_logg','k_logg',
    '[Fe/H]','e_[Fe/H]','s_[Fe/H]','k_[Fe/H]',
    'sep_arcsec',
    'opt_RA_deg','opt_Dec_deg',
    'opt_radial_velocity','opt_radial_velocity_error',
    'opt_phot_g_mean_mag','opt_bp_rp','opt_ruwe',
    'opt_rbin_ell_key','opt_tbin_ell_key','opt_r_pc_ell',
    'opt_rbin_circ_key','opt_r_pc_circ'
]
write_csv(matched_all[[c for c in qa_cols if c in matched_all.columns]], "draco_spec_qamatch")
write_csv(matched_best[[c for c in qa_cols if c in matched_best.columns]], "draco_spec_qamatch_best")
print("[done] outputs in", outdir.resolve())