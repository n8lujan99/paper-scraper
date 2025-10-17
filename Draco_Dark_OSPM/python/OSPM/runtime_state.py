# runtime_state.py
###################################################################################
# Imports
###################################################################################
import numpy as np 

###################################################################################
model_params = {
    "hole": None,
    "distance": None,
    "angrad": None,
    "vmin": None,
    "vmax": None,
    "cslope": None,
    "xinclin": None,
    "sini": None,
    "cosi": None,
    "vmult": None,
    "perimin": None,
    "apomax": None,
    "rpotmin": None,
    "rpotmax": None,
    "rlgpotmin": None,
    "rlgpotmax": None,
    "gdennorm": None,
    "halo": {
        "ihalo": None,
        "dis": None,
        "v0": None,
        "rc": None,
        "qdm": None,
        "xmgamma": None,
        "rsgamma": None,
        "gamma": None,
        "cnfw": None,
        "rsnfw": None
    }
}

###################################################################################
# Shared arrays or placeholders (optional)
###################################################################################
wphase = np.array([])    # phase-space weights (used by librarian)
Norbtot = 0              # number of orbits generated