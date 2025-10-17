# ============================================================================================================
# Orbit Super Position Modeling -Spectrum handler
# ============================================================================================================
from astropy.io import fits
from astropy.stats import biweight_location, biweight_scale
from astropy.io.votable import parse_single_table
import io, numpy as np
from astropy.io.fits.verify import VerifyError
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from typing import Tuple

# ==================================================================================
# USER SUPPLIES LOCAL FILES HERE
# ==================================================================================
input_files = ["spectrum1.fits", "spectrum2.txt", "spectrum3.xml"]
# ============================================================================================================
# Target spectrum acquisition via SSA
# ===========================================================================================================
POS_RA_DEG  = 260.051
POS_DEC_DEG = 57.915
SEARCH_SIZE_DEG = 0.15
BAND_M = (2.2e-6, 2.4e-6)  # K band
MAX_SPECTRA = 12           # how many to download/stack
OUTDIR      = "Draco_Data"
os.makedirs(OUTDIR, exist_ok=True)

def save_good_spec(waveA, flux, loA=22000.0, hiA=24000.0, outfile="Draco_Data/good.spec"):
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    mfin = np.isfinite(waveA) & np.isfinite(flux)
    waveA, flux = waveA[mfin], flux[mfin]
    sel = (waveA >= loA) & (waveA <= hiA)
    if sel.sum() >= 10:
        waveA, flux = waveA[sel], flux[sel]
    med = np.nanmedian(flux) if np.isfinite(flux).any() else 1.0
    if not np.isfinite(med) or med == 0.0:
        med = 1.0
    fluxn = flux / med
    order = np.argsort(waveA)
    np.savetxt(outfile, np.c_[waveA[order], fluxn[order]], fmt="%.2f %.6f")
    print(f"Saved template: {outfile}  (N={order.size})")

def wave_flux_A(path_or_url):
    """
    Robust spectrum reader for FITS, VOTable XML, or ASCII (whitespace/CSV).
    Returns (wave_Angstrom, flux).
    """
    # ---- 1) Try FITS first
    try:
        with fits.open(path_or_url, memmap=False) as hdul:
            hdu = next((h for h in hdul if getattr(h, "data", None) is not None and hasattr(h, "columns")), None)
            if hdu is not None:
                cols_lower = [c.lower() for c in hdu.columns.names]
                def pick(cands):
                    for key in cands:
                        for i, cn in enumerate(cols_lower):
                            if key in cn:
                                return hdu.columns.names[i]
                    return None
                wcol = pick(["wavelength","wave","lambda","lam","spec_val"])
                fcol = pick(["flux","flam","f_lambda","intensity","spec","spectrum"])
                if wcol is None or fcol is None:
                    raise ValueError("No wavelength/flux columns in FITS table.")
                wave = np.asarray(hdu.data[wcol], float).squeeze()
                flux = np.asarray(hdu.data[fcol], float).squeeze()
                unit = (hdu.columns[wcol].unit or "").lower()
                # convert to Å
                if 'nm' in unit:        wave *= 10.0
                elif 'um' in unit or 'micron' in unit:
                    wave *= 1.0e4
                # else assume Å
                return np.asarray(wave, float), np.asarray(flux, float)
            else:
                # image-like with linear wavelength axis
                hdr  = hdul[0].header
                data = np.asarray(hdul[0].data, float).squeeze()
                crval = hdr.get('CRVAL1'); cdelt = hdr.get('CDELT1')
                if crval is None or cdelt is None:
                    raise ValueError("No CRVAL1/CDELT1 in FITS header.")
                pix  = np.arange(data.size)
                wave = crval + cdelt*pix
                cunit = (hdr.get('CUNIT1') or "").lower()
                if 'nm' in cunit:        wave *= 10.0
                elif 'um' in cunit or 'micron' in cunit:
                    wave *= 1.0e4
                return np.asarray(wave, float), data
    except (OSError, VerifyError, ValueError):  # not a FITS or not parseable as FITS
        pass

    # VOTable
    try:
        with open(path_or_url, "rb") as resp:
            raw = resp.read()
        vot = parse_single_table(io.BytesIO(raw))
        tab = vot.to_table()
        lower = {c.lower(): c for c in tab.colnames}
        def pick(names):
            for n in names:
                if n in lower: return lower[n]
            for c in tab.colnames:
                if any(n in c.lower() for n in names): return c
            return None
        wname = pick(["wavelength","wave","lambda","lam","spec_val"])
        fname = pick(["flux","flam","f_lambda","intensity","spec","spectrum"])
        if wname is None or fname is None:
            raise ValueError("No wavelength/flux columns in VOTable.")
        wave = np.asarray(tab[wname], float).squeeze()
        flux = np.asarray(tab[fname], float).squeeze()
        unit = (tab[wname].unit or "").lower()
        if 'nm' in unit: wave *= 10.0
        elif 'um' in unit or 'micron' in unit: wave *= 1.0e4
        elif unit in ('m','meter','metre'): wave *= 1.0e10
        mx = np.nanmax(wave) if np.isfinite(wave).any() else np.nan
        if np.isfinite(mx):
            if mx <= 1.0: wave *= 1.0e10
            elif mx <= 50.0: wave *= 1.0e4
        return wave, flux
    except Exception:
        pass

    # ASCII
    try:
        data = np.genfromtxt(path_or_url, comments="#", delimiter=None)
        if data.ndim == 1 and data.size >= 2:
            data = data.reshape(-1, data.size)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError("ASCII file does not have >=2 columns.")
        wave, flux = data[:,0], data[:,1]
        mx = np.nanmax(wave) if np.isfinite(wave).any() else np.nan
        if np.isfinite(mx):
            if mx <= 1.0: wave *= 1.0e10
            elif mx <= 50.0: wave *= 1.0e4
        return wave, flux
    except Exception as e:
        raise ValueError(f"Unrecognized spectrum format at {path_or_url}: {e}")

def LOC(x, axis=None):
    x = np.asarray(x, float)
    return np.nanmedian(x, axis=axis) if np.isnan(x).any() else biweight_location(x, axis=axis)

def Scale(x, axis=None):
    x = np.asarray(x, float)
    if np.isnan(x).any():
        med = np.nanmedian(x, axis=axis, keepdims=True)
        mad = np.nanmedian(np.abs(x - med), axis=axis)
        return 1.4826 * mad
    return biweight_scale(x, axis=axis)

def get_sky(spectra_stack: np.ndarray, method: str = "biweight") -> Tuple[np.ndarray, np.ndarray]:
    spectra_stack = np.asarray(spectra_stack, dtype=float)
    if spectra_stack.ndim != 2:
        raise ValueError("spectra_stack must be (N_spectra, N_pixels).")
    if method == "biweight":
        sky  = LOC(spectra_stack, axis=0)
        scat = Scale(spectra_stack, axis=0)
    else:
        sky  = np.nanmedian(spectra_stack, axis=0)
        scat = 1.4826 * np.nanmedian(np.abs(spectra_stack - sky[None,:]), axis=0)
    return sky, scat

 
candidates = []
for fn in input_files:
    try:
        waveA, flux = wave_flux_A(fn)
        m = np.isfinite(waveA) & np.isfinite(flux)
        waveA, flux = waveA[m], flux[m]
        order = np.argsort(waveA)
        waveA, flux = waveA[order], flux[order]
        if waveA.size < 20: continue
        inK = (waveA >= 22000.0) & (waveA <= 24000.0)
        if inK.sum() < 5: continue
        overlap = float(waveA[inK].ptp())
        candidates.append((overlap, waveA.size, fn, waveA, flux))
    except Exception as e:
        print(f"[WARN] Failed on {fn}: {e}")

candidates.sort(key=lambda t: (t[0], t[1]), reverse=True)

grid = np.arange(22000.0, 24000.0 + 0.5, 2.0)
stack = []
index_rows = []
for k,(ov,npts,path,wA,f) in enumerate(candidates,1):
    sel = (wA >= grid.min()) & (wA <= grid.max())
    if sel.sum()<5: continue
    wA,f = wA[sel], f[sel]
    med = np.nanmedian(f) if np.isfinite(f).any() else 1.0
    if not np.isfinite(med) or med==0.0: med=1.0
    f = f/med
    fr = np.interp(grid,wA,f,left=np.nan,right=np.nan)
    stack.append(fr)
    index_rows.append(dict(rank=k, overlap_A=ov, n_points=npts, path=path))

if len(stack)==0:
    raise RuntimeError("After resampling, no spectra remained to stack.")

stack = np.vstack(stack)
idx_df = pd.DataFrame(index_rows)
idx_df.to_csv(os.path.join(OUTDIR,"stack_index.csv"),index=False)

sky, scat = get_sky(stack, method="biweight")
sky_df = pd.DataFrame({"wave_A": grid, "sky": sky, "sigma": scat})
sky_csv = os.path.join(OUTDIR, "sky_kband.csv")
sky_df.to_csv(sky_csv, index=False)

best_wA,best_f = candidates[0][3], candidates[0][4]
sel_best = (best_wA >= 22000.0) & (best_wA <= 24000.0)
if sel_best.any():
    best_wA, best_f = best_wA[sel_best], best_f[sel_best]
save_good_spec(best_wA, best_f, loA=22000.0, hiA=24000.0, outfile=os.path.join(OUTDIR,"good.spec"))
print(f"[OK] Stacked {stack.shape[0]} spectra on {grid.size} pixels.")























'''
# ============================================================================================================
# Pull spectrum
# ============================================================================================================
def _infer_wave_A(w):
    w = np.asarray(w, float)
    if w.size == 0 or not np.isfinite(w).any():
        return w
    m = np.nanmax(w)
    # heuristics: Å (~2.2e4), μm (~2.2), m (~2.2e-6)
    if m > 1e3:      # already Å
        return w
    elif m > 1.0:    # μm
        return w * 1.0e4
    else:            # meters
        return w * 1.0e10

def wave_flux_A(path):
    with fits.open(path, memmap=False) as hdul:
        # find first table HDU
        tab = None
        for h in hdul:
            if hasattr(h, "columns"):
                hdu = h; break
        if hdu is None:
            raise ValueError("No table HDU with columns found.")
        cols = [c.lower() for c in hdu.columns.names]
        def pick(keys):
            for key in keys:
                for i, cn in enumerate(cols):
                    if key in cn:
                        return hdu.columns.names[i]
            return None
        wcol = pick(["wave", "wavelength", "lambda", "lam"])
        fcol = pick(["flux", "flam", "f_lambda", "intensity"])
        if wcol is None or fcol is None:
            raise ValueError(f"No wavelength/flux columns in {hdu.columns.names}")
        wave = np.asarray(hdu.data[wcol], float).squeeze()
        flux = np.asarray(hdu.data[fcol], float).squeeze()
    return _infer_wave_A(wave), flux


# ============================================================================================================
# Functions
# ============================================================================================================
def wave_unit_to_angstrom(wave):
    """Return wavelength array in Angstrom from Angstrom/µm/m by heuristic."""
    w = np.asarray(wave, dtype=float)
    if not np.isfinite(w).any():
        return w
    mx = np.nanmax(w)
    # typical scales: 2.2e4 Å, 2.2 µm, 2.2e-6 m
    if mx > 1e3:        # looks like Angstrom already
        return w
    elif mx > 1.0:      # likely microns
        return w * 1.0e4
    else:               # meters
        return w * 1.0e10

def _extract_wave_flux_from_fits(path):
    """Open a FITS spectrum and return (wave[Å], flux). Tries common column names."""
    with fits.open(path, memmap=False) as hdul:
        # find first table HDU with columns
        hdu = None
        for h in hdul:
            if hasattr(h, "columns"):
                hdu = h; break
        if hdu is None:
            raise ValueError("No table HDU with columns found in FITS.")
        tbl = hdu.data
        # find wavelength column
        colnames = [c.lower() for c in hdu.columns.names]
        def _find(names):
            for want in names:
                for i, cn in enumerate(colnames):
                    if want in cn:
                        return hdu.columns.names[i]
            return None
        wcol = _find(["wave", "wavelength", "lambda", "lam"])
        fcol = _find(["flux", "flam", "f_lambda"])
        if wcol is None or fcol is None:
            raise ValueError(f"Could not find wavelength/flux columns in {hdu.columns.names}")
        wave = np.asarray(tbl[wcol], dtype=float).squeeze()
        flux = np.asarray(tbl[fcol], dtype=float).squeeze()
        waveA = wave_unit_to_angstrom(wave)
        return waveA, flux

# ============================================================================================================
# Data binding 
# ============================================================================================================

def get_1sigma(spec_path: str | Path, wave: float, nflimh: int = 3, write_to: str | Path = "get1sig.out") -> float:
    """
    Estimate 1σ at a target wavelength by quadrature-summing (x3/x9)^2 within a small window.

    Expects a whitespace-delimited file with columns where:
    col0 = wavelength, col2 = x3, col8 = x9
    Uses a window of +/- nflimh rows around the nearest wavelength to 'wave'.

    Writes the scalar result to 'write_to' and returns it.
    """
    df = pd.read_csv(spec_path, sep=r"\s+", header=None, engine="python")
    n = len(df)
    if n < 10:
        f1sig = 0.0
    else:
        w  = df.iloc[:, 0].to_numpy(dtype=float)
        x3 = df.iloc[:, 2].to_numpy(dtype=float)
        x9 = df.iloc[:, 8].to_numpy(dtype=float)

        # safe division: se = x3/x9 where x9>0 else 0
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
'''