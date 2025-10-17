# General Purpose Fits reader
# Gaia FITS for testing
#https://gea.esac.esa.int/archive

# Data request for Draco
#def extract_columns(filename, hdu=2, col_indices=None):
#    with fits.open(filename) as hdul:
#        data = hdul[hdu].data
#        return [data.field(col) for col in col_indices]
######################################################################################################


import os
import logging
import numpy as np
import scipy as sc
import pandas as pd
from scipy.interpolate import interp1d
import astropy as ast
from astropy.io import fits
from astropy.io.votable import parse
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import match_coordinates_sky
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle as pckl

######################################################################################################

GAIA_FILENAME  = "1754598165996O-result.fits" #Only shows 5 stars to match within 5"
GAIA_OUTFILE = "gaia_extracted_sources.txt"  
HDU_INDEX = 1  # Confirmed BinTableHDU

stellar_vel = "stellar_vel.txt" # Stellar Data from J.T. Kleyna 2001
votable = parse("simbad_result")
table = votable.get_first_table().to_table()
simbad_df = table.to_pandas()
# Show available columns
print("simbad_df Shape")
print(simbad_df.head())
print(simbad_df.columns)
simbad_df_clean = simbad_df.dropna(subset=["FLUX_V"])
print(simbad_df_clean.head())
# Save the full DataFrame as a text file
simbad_df_clean.to_csv("simbad_full.csv", index=False)

######################################################################################################
# Starting With GAIA
with fits.open(GAIA_FILENAME) as hdul:
    print(f"Number of HDUs: {len(hdul)}\n")
    for i, hdu in enumerate(hdul):
        print(f"HDU {i}: {type(hdu)}  name={hdu.name}")
        print("\n" + "="*60 + "\n")

columns = ["source_id", "ra", "ra_error", "dec", "dec_error", "parallax", "parallax_error",
    "pmra", "pmra_error", "pmdec", "pmdec_error", "radial_velocity", "radial_velocity_error",
    "l", "b", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "ruwe", "distance_gspphot", "phot_bp_rp_excess_factor", "duplicated_source"]

with fits.open(GAIA_FILENAME) as hdul:
    data = hdul[1].data  # Adjust HDU index if needed

    df = pd.DataFrame({col: np.array(data[col]).byteswap().newbyteorder() if data[col].dtype.byteorder == '>' else data[col] for col in columns})
    # Apply filters
    good = (
        (df["ruwe"] < 1.4) &
        (~df["duplicated_source"]) &
        (df["parallax"] > 0) &
        (df["parallax_error"] > 0) &
        (df["parallax"] / df["parallax_error"] > 5))

    clean_df = df[good].copy()
    print(f"Loaded {len(df)} rows, retained {len(clean_df)} after filtering.")

##################################################################################################


def process_stellar_velocity_data(infile="stellar_vel.txt", all_outfile="orbits.csv", good_outfile="G_orbits.csv"):
    columns = ["ID", "V", "RA_h", "RA_m", "RA_s", "Dec_d", "Dec_m", "Dec_s", "v_helio", "v_err", "RTD"]
    stell_vel_df = pd.read_csv(infile, sep=",", names=columns, header =0)
    stell_vel_df = stell_vel_df.sort_values(by="ID").reset_index(drop=True)
    probable_ids = [3, 13, 18, 30, 64, 75, 80, 140, 153, 164, 203, 236, 254, 255,
                    273, 276, 277, 281, 293, 304, 305, 333, 336, 345, 374, 376, 401]
    orbits = stell_vel_df.copy()
    G_orbits = stell_vel_df[~stell_vel_df["ID"].isin(probable_ids)].copy()
    # Save to CSV
    orbits.to_csv(all_outfile, index=False)
    G_orbits.to_csv(good_outfile, index=False)
    return orbits, G_orbits

orbits, G_orbits = process_stellar_velocity_data()
print("All Orbits:")
print(orbits.head())

# Convert RA/Dec from h/m/s and d/m/s to decimal degrees
def convert_radec(row):
    coord = SkyCoord(ra=row["RA_h"] * u.hour + row["RA_m"] * u.minute + row["RA_s"] * u.second,
                     dec=row["Dec_d"] * u.deg + row["Dec_m"] * u.arcmin + row["Dec_s"] * u.arcsec)
    return coord.ra.deg, coord.dec.deg

orbits["RA_deg"], orbits["Dec_deg"] = zip(*orbits.apply(convert_radec, axis=1))
print("RA/Dec deg columns created:")
print(orbits[["RA_deg", "Dec_deg"]].head())

# Rename Gaia columns to match
clean_df = clean_df.rename(columns={"ra": "RA_deg", "dec": "Dec_deg"})

# Filter Gaia stars to Kleyna region
ra_min, ra_max = orbits["RA_deg"].min() - 0.1, orbits["RA_deg"].max() + 0.1
dec_min, dec_max = orbits["Dec_deg"].min() - 0.1, orbits["Dec_deg"].max() + 0.1
gaia_cut = clean_df[(clean_df["RA_deg"] >= ra_min) & (clean_df["RA_deg"] <= ra_max) & (clean_df["Dec_deg"] >= dec_min) & (clean_df["Dec_deg"] <= dec_max)]

print(f"Filtered Gaia stars from {len(clean_df)} → {len(gaia_cut)}")

# Save filtered Gaia and build coordinates
clean_df.to_csv("gaia_clean.csv", index=False)
coords_gaia = SkyCoord(ra=gaia_cut["RA_deg"].values * u.deg, dec=gaia_cut["Dec_deg"].values * u.deg)
print("Cleaned Gaia data saved to gaia_clean.csv")
print(gaia_cut.head())



# ------------------ Matching ------------------
gaia_cut = gaia_cut.rename(columns={"ra": "RA_deg", "dec": "Dec_deg"})

# Set max separation threshold
max_sep = 10.0 * u.arcsec

# Build SkyCoord objects
coords_kleyna = SkyCoord(ra=orbits["RA_deg"].values * u.deg, dec=orbits["Dec_deg"].values * u.deg)
coords_gaia = SkyCoord(ra=gaia_cut["RA_deg"].values * u.deg, dec=gaia_cut["Dec_deg"].values * u.deg)
coords_simbad = SkyCoord(ra=simbad_df["RA_d"].values * u.deg, dec=simbad_df["DEC_d"].values * u.deg)

# Match to Gaia
idx_gaia, d2d_gaia, _ = coords_kleyna.match_to_catalog_sky(coords_gaia)
match_mask_gaia = d2d_gaia < max_sep
matched_orbits_gaia = orbits[match_mask_gaia].reset_index(drop=True)
matched_gaia = gaia_cut.iloc[idx_gaia[match_mask_gaia]].reset_index(drop=True)
joined_df = matched_orbits_gaia.reset_index(drop=True).add_prefix("kleyna_").join(matched_gaia.reset_index(drop=True).add_prefix("gaia_"))
assert matched_orbits_gaia["ID"].is_unique

# Match remaining to SIMBAD
unmatched_orbits = orbits[~match_mask_gaia].reset_index(drop=True)
coords_unmatched = SkyCoord(ra=unmatched_orbits["RA_deg"].values * u.deg, dec=unmatched_orbits["Dec_deg"].values * u.deg)
idx_simbad, d2d_simbad, _ = coords_unmatched.match_to_catalog_sky(coords_simbad)
match_mask_simbad = d2d_simbad < max_sep
matched_orbits_simbad = unmatched_orbits[match_mask_simbad].reset_index(drop=True)
matched_simbad = simbad_df.iloc[idx_simbad[match_mask_simbad]].reset_index(drop=True)
joined_simbad_df = matched_orbits_simbad.reset_index(drop=True).add_prefix("kleyna_").join(matched_simbad.reset_index(drop=True).add_prefix("simbad_"))
assert matched_orbits_simbad["ID"].is_unique

# Final combined DataFrame
combined_df = pd.concat([joined_df, joined_simbad_df], axis=0).reset_index(drop=True)
combined_df.to_csv("combined_kleyna_matches.csv", index=False)
print(f"Final matched sources: {len(combined_df)} saved to combined_kleyna_matches.csv")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


'''
print("combined_df shape:", combined_df.shape)
print(combined_df)
print("matched_orbits_gaia:", len(matched_orbits_gaia))
print("matched_gaia:", len(matched_gaia))
print("joined_df:", len(joined_df))

print("matched_orbits_simbad:", len(matched_orbits_simbad))
print("matched_simbad:", len(matched_simbad))
print("joined_simbad_df:", len(joined_simbad_df))

# Save to file
combined_df.to_csv("combined_matched_sources.csv", index=False)

# specfic to desi_hetdex
def extract_redshift_match(filename, hdu=2, nrows=None):
    """Extract z_dex, z_desi, S/N, fiber info, and offset."""
    with fits.open(filename) as hdul:
        data = hdul[hdu].data
        n_total = len(data)
        n = n_total if nrows is None else min(nrows, n_total)
        output = []
        for i in range(n):
            try:
                idd     = data.field(8)[i]
                zdex    = data.field(19)[i]
                zdesi   = data.field(114)[i]
                sn      = data.field(14)[i]
                ifib    = data.field(14)[i]  # May conflict with sn
                ifibstat= data.field(28)[i]
                dx      = data.field(87)[i]
                dy      = data.field(89)[i]
                dist    = np.sqrt(dx**2 + dy**2)
                output.append((idd, zdex, zdesi, sn, ifib, ifibstat, dist))
            except Exception:
                continue
    np.savetxt("out.txt", output, fmt='%s')
extract_columns("in.fits", hdu=2, col_indices=[7,8,1], nrows=500)
extract_spectrum_vector("in.fits")
extract_redshift_match("desi-hetdex.fits")

file1 = input("Enter FITS filename: ").strip()
campname = input("Enter amplifier name (2 characters): ").strip()

try:
    with fits.open(file1, mode='update') as hdul:
        hdr = hdul[1].header  # Assuming extension 1
        hdr['AMPNAME'] = (campname, 'Amplifier in use')
        hdul.flush()
except Exception as e:
    print(f"Error opening image: {file1}\n{e}")

def batch_extract_spectra(filename, hdu=2, wave_col=20, flux_col=21, err_col=22, id_col=8):
    """Extracts spectra from a FITS catalog and saves them as individual .spec files."""
    from astropy.io import fits
    import numpy as np

    with fits.open(filename) as hdul:
        data = hdul[hdu].data
        for i in range(len(data)):
            try:
                obj_id = int(data.field(id_col)[i])
                wave   = data.field(wave_col)[i]
                flux   = data.field(flux_col)[i]
                fluxe  = data.field(err_col)[i]
                outf   = f"{obj_id:010d}.spec"
                np.savetxt(outf, np.column_stack((wave, flux, fluxe)))
                print(f"Wrote {outf}")
            except Exception as e:
                print(f"Row {i} failed: {e}")

def extract_df_sim(filename, hdu=2, col1=1, col2=2, col3=3, outfile="out.txt"):
    from astropy.io import fits
    import numpy as np

    with fits.open(filename) as hdul:
        data = hdul[hdu].data
        dx1 = data.field(col1).astype(np.float32)
        dx2 = (data.field(col2) / 1e42).astype(np.float32)
        dx3 = (data.field(col3) / 1e-17).astype(np.float32)
        result = np.column_stack((np.arange(1, len(dx1)+1), dx1, dx2, dx3))
        np.savetxt(outfile, result, fmt="%d %.6e %.6e %.6e")
        print(f"Wrote {len(dx1)} rows to {outfile}")

# Set up logging
logging.basicConfig(filename='pipeline.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def safe_open_fits(file_path, mode='readonly'):
    """Open a FITS file safely with logging."""
    if not os.path.exists(file_path):
        logging.error(f"File not found: {file_path}")
        return None
    try:
        return fits.open(file_path, mode=mode, memmap=True)
    except Exception as e:
        logging.error(f"Failed to open FITS file {file_path}: {e}")
        return None

def validate_array_shape(data, narrm=NARRM):
    """Ensure array dimensions are within allowed bounds."""
    if data.ndim == 1:
        data = np.expand_dims(data, axis=1)
    nrow, ncol = data.shape
    if nrow > narrm or ncol > narrm:
        logging.warning(f"Array shape too large: ({nrow}, {ncol}), max allowed is {narrm}")
        return None
    return data

def check_row_continuity(data, icen=ICEN, filename="unknown"):
    """Inspect central row of data for consistency."""
    try:
        if icen >= data.shape[0]:
            logging.warning(f"{filename}: ICEN {icen} exceeds data bounds {data.shape[0]}")
            return
        x1 = data[icen, 0]
        for j in range(1, min(112, data.shape[1])):
            xdiff = data[icen, j] - x1
            x1 = data[icen, j]
            if x1 > 0. and abs(xdiff) < 5.0:
                logging.info(f"{filename[:25]:<25} col={j:3d} Δ={xdiff:5.2f} val={x1:7.2f}")
    except Exception as e:
        logging.error(f"{filename}: row continuity check failed: {e}")

def process_listin_file(list_file=LIST_FILE, narrm=NARRM):
    """Process FITS files listed in a text file and run quality checks."""
    xd = np.zeros((narrm, narrm))
    if not os.path.exists(list_file):
        logging.error(f"Missing list file: {list_file}")
        return

    with open(list_file, 'r') as f:
        for line in f:
            file1 = line.strip()
            if not file1:
                continue

            hdul = safe_open_fits(file1)
            if hdul is None:
                continue

            try:
                data = hdul[1].data if len(hdul) > 1 else hdul[0].data
                data = validate_array_shape(data, narrm=narrm)
                if data is None:
                    continue

                xd[:data.shape[0], :data.shape[1]] = data
                check_row_continuity(xd, icen=ICEN, filename=file1)

            except Exception as e:
                logging.error(f"{file1}: processing failed: {e}")
            finally:
                hdul.close()

process_listin_file()

def inject_ampname(file1, campname):
    """Inject amplifier name into FITS header (extension 1)."""
    try:
        with fits.open(file1, mode='update') as hdul:
            hdr = hdul[1].header  # or adjust if your header is in another HDU
            hdr['AMPNAME'] = (campname, 'Amplifier in use')
            hdul.flush()
    except Exception as e:
        logging.error(f"Error injecting AMPNAME into {file1}: {e}")


def general_fits_reader(list_file=LIST_FILE, inject_amp=False, amp_filename=None, amp_name=None, run_cosmos_filter=False, cosmos_file=None, run_redshift_match=False, redshift_file=None, run_batch_extract=False, batch_file=None, run_vector_extract=False, vector_file=None, vector_row=2, run_df_extract=False, df_file=None, df_outfile="out.txt", cosmos_outfile="out.txt", redshift_outfile="out.txt", verbose=True):

    if inject_amp and amp_filename and amp_name:
        inject_ampname(amp_filename, amp_name)
        if verbose:
            print(f"Injected amplifier name '{amp_name}' into {amp_filename}")

    if run_cosmos_filter and cosmos_file:
        extract_cosmos_sample(cosmos_file, outfile=cosmos_outfile)
        if verbose:
            print(f"Filtered COSMOS catalog written to {cosmos_outfile}")

    if run_redshift_match and redshift_file:
        extract_redshift_match(redshift_file, outfile=redshift_outfile)
        if verbose:
            print(f"Redshift match data written to {redshift_outfile}")

    if run_batch_extract and batch_file:
        batch_extract_spectra(batch_file)
        if verbose:
            print(f"Batch spectra extracted from {batch_file}")

    if run_vector_extract and vector_file:
        extract_spectrum_vector(vector_file, row=vector_row)
        if verbose:
            print(f"Spectrum vector extracted from row {vector_row} of {vector_file}")

    if run_df_extract and df_file:
        extract_df_sim(df_file, outfile=df_outfile)
        if verbose:
            print(f"DF simulation data written to {df_outfile}")

    # Always process listin
    process_listin_file(list_file=list_file)
    if verbose:
        print(f"Finished processing all files listed in {list_file}")

general_fits_reader(True, "image.fits", "A1", True, "COSMOS2020_CLASSIC_R1_v2.0.fits", True, "desi-hetdex.fits", True, "specs.fits", True, "specs.fits", 3, True, "df_sim.fits", "df_output.txt", "cosmos_filtered.txt", "redshift_match.txt", True)


# Requires a COSMOS 2020 catalog
# https://irsa.ipac.caltech.edu/data/COSMOS
# wget https://irsa.ipac.caltech.edu/data/COSMOS/files/cosmos2020/v2.0/COSMOS2020_CLASSIC_R1_v2.0.fits

def extract_cosmos_sample(filename, zmin=1.9, zmax=3.5, gmin=14, gmax=28):
    """Filter COSMOS FITS file based on redshift and g-magnitude."""
    t = Table.read(filename, hdu=2)
    mask = (t["ZPHOT"] > zmin) & (t["ZPHOT"] < zmax) & (t["MAG_G"] > gmin) & (t["MAG_G"] < gmax)
    filtered = t[mask]
    filtered["ID", "RA", "DEC", "MAG_G", "ZPHOT", "SFR"].write("out.txt", format="ascii")
extract_cosmos_sample("COSMOS2020_CLASSIC_R1_v2.0.fits")
'''