from astroquery.simbad import Simbad

# Allow longer queries
Simbad.TIMEOUT = 120

# Add the extra fields you need
Simbad.add_votable_fields("ra(d)", "dec(d)", "pmra", "pmdec", "coo_bibcode")

# Object list
objects = [
    "[BFO2002] Dra 10201",
    "BASV 68",
    "SDSS J172015.85+575502.7",
    "SDSS J172007.45+575432.6",
    "[BS61] K",
    "[KWE2002] 254",
    "[BFO2002] Dra 10312"
]

# Query SIMBAD
result = Simbad.query_objects(objects)

# Show the main ID, coords, proper motions, and the bibcode for the source of coords
print(result["MAIN_ID", "RA_d", "DEC_d", "PMRA", "PMDEC", "COO_BIBCODE"])
