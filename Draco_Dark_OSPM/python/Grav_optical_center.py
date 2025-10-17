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
FILE_DRACO          = "Draco_stars_650pc.csv"
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
# getdisp.f / getdispb.f / getfielda.f
# =======================

# continue at getflux.f
# =======================
# Getfield.f  
# =======================

# =======================
# FUNCTIONS
# =======================

# =======================
# FUNCTIONS
# =======================

# =======================
# FUNCTIONS
# =======================

# =======================
# FUNCTIONS
# =======================