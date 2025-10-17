# stellar velocity Calculator
# Utilizes the stars inital and final postions to calculate the angular velocity
# We then make a simple csv file to store the inital and final positions, tangential, 
# and angular velocity that will be used for OSPM. at the end we will do a plot of the projected data

import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

df = pd.read_csv('Draco_stars_data.csv')

# Function to calculate angular displacement between two positions (RA, Dec)
def angular_displacement(ra1, dec1, ra2, dec2):
    ra1, dec1 = np.deg2rad(ra1), np.deg2rad(dec1)
    ra2, dec2 = np.deg2rad(ra2), np.deg2rad(dec2)
    
    # spherical law of cosines formula
    cos_delta_theta = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    
    # Ensure the value is between -1 and 1 due to floating point precision issues
    cos_delta_theta = np.clip(cos_delta_theta, -1.0, 1.0)
    
    # Calculate angular displacement in radians
    delta_theta = np.arccos(cos_delta_theta)
    
    # Convert to arcseconds
    delta_theta_arcsec = np.rad2deg(delta_theta) * 3600
    return delta_theta_arcsec

# Function to compute the angular velocity (omega)
def angular_velocity(delta_theta, delta_t):
    # Angular velocity = angular displacement / time
    omega = delta_theta / delta_t  # rad per year
    return omega

# Function to compute tangential velocity (v_tangential)
def tangential_velocity(omega, distance):
    # v_tangential = r * omega
    v_tangential = distance * omega  # in km/s
    return v_tangential

# Prepare to store output results
results = []

# Process each star's data
for index, row in df.iterrows():
    # Extract data from the row
    ra_initial = row['kleyna_RA_degs']  # Initial RA from Kleyna data
    dec_initial = row['kleyna_Dec_degs']  # Initial Dec from Kleyna data
    ra_final = row['Gaia_ra']  # Final RA from Gaia data
    dec_final = row['dec']  # Final Dec from Gaia data
    distance = 1 / row['parallax'] * u.pc  # Distance from parallax (in parsecs)
    
    # Convert the distance to km (1 parsec ≈ 3.086 x 10^13 km)
    distance_km = distance.to(u.km).value
    
    # Calculate angular displacement in arcseconds
    delta_theta = angular_displacement(ra_initial, dec_initial, ra_final, dec_final)
    
    # Time interval (10 years)
    delta_t = 10  # in years
    
    # Calculate angular velocity (omega) in arcseconds per year
    omega = angular_velocity(delta_theta, delta_t)
    
    # Calculate tangential velocity in km/s
    v_tangential = tangential_velocity(omega, distance_km)
    
    # Calculate circular velocity (same as tangential velocity in this case)
    v_circular = v_tangential  # assuming circular orbit

    # Store results
    results.append({
        'kleyna_ID': row['kleyna_ID'],
        'RA_initial': ra_initial,
        'Dec_initial': dec_initial,
        'RA_final': ra_final,
        'Dec_final': dec_final,
        'distance_km': distance_km,
        'delta_theta_arcsec': delta_theta,
        'omega_arcsec_per_year': omega,
        'v_tangential_km_s': v_tangential,
        'v_circular_km_s': v_circular
    })

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(results)
results_df.to_csv('star_velocities.csv', index=False)

print("Script completed. Results saved in 'star_velocities.csv'.")
