import pandas as pd 
import numpy as np 
import astropy as ast
import astropy.units as u
import astropy.constants as const
import os 


G = const.G
star_data = pd.read_csv('star_velocities.csv')

# Extract the parallax and parallax error values
parallax = star_data['parallax']
parallax_error = star_data['parallax_error']
parallax_over_error = star_data['parallax_over_error']
star_data['distance_km'] = np.nan


# Iterate over the rows and validate parallax, set invalid values to NaN
for index, row in star_data.iterrows():
    # Check for invalid parallax (negative or less than the parallax error in the next row)
    if row['parallax'] <= 0 or (row['parallax'] < row['parallax_error']):
        # Set the distance and parallax-related columns to NaN for invalid rows
        star_data.at[index, 'distance_km'] = np.nan
        star_data.at[index, 'parallax'] = np.nan
        star_data.at[index, 'parallax_error'] = np.nan
        star_data.at[index, 'parallax_over_error'] = np.nan
    else:
        # Calculate valid distance in km for valid parallax values
        star_data.at[index, 'distance_km'] = 1 / row['parallax'] * 1e3  # Distance in km

# Extract circular velocity (in km/s) - no conversion needed since it's already in km/s
v_circ = star_data['v_circular_km_s']  # 'stellar circular velocity in km/s'


def estimate_mass(v_circular, r):
    # Calculate mass in kg (using gravitational constant G)
    M = (v_circular ** 2 * r) / G  # Mass in kg
    return M

def gravitational_potential(M, r):
    # Calculate the gravitational potential in J/kg
    phi = - (G * M) / r  # Gravitational potential in J/kg
    return phi

results = []

for index, row in star_data.iterrows():
    # Skip stars with invalid distance (NaN)
    if np.isnan(row['distance_km']):
        continue
    
    # Estimate mass using circular velocity and distance for each star
    M = estimate_mass(v_circ[index], row['distance_km'])  # Mass in kg
    
    # Calculate gravitational potential for each star at its distance from the center
    phi = gravitational_potential(M, row['distance_km'])  # Gravitational potential in J/kg
    
    # Store results
    results.append({
        'kleyna_ID': row['kleyna_ID'],
        'v_circular_km_s': v_circ[index],
        'estimated_mass_kg': M,  # Mass in kg
        'gravitational_potential_J_kg': phi  # Gravitational potential in J/kg
    })

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(results)
results_df.to_csv('draco_grav.csv', index=False)

print