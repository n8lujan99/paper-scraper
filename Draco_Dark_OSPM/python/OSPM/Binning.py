# Binning conversions for OSPM

import numpy as np
import pandas as pd


def IRneeR(r, a, b, cval):
    """Converts physical radius r (on [0,1]) into radial bin index ir."""
    return int(np.log(a * r / b + cval) / a + 0.5)

def RneeIR(ir, a, b, cval):
    """Converts radial bin index ir into physical radius r (bin center)."""
    return (b / a) * (np.exp(a * float(ir)) - cval)

def IVneeV(v, ivmax):
    """Converts physical v = sin(th) into angular bin index iv."""
    return int(np.sign(v) * min(round(float(ivmax) * abs(v) + 0.5), ivmax))

def VneeIV(iv, ivmax):
    """Converts angular bin index iv into physical v = sin(th) (bin center)."""
    return (float(iv - 1) + 0.5) / float(ivmax)

def IRCneeR(r, a, b, cval, irrat):
    """Converts physical radius r into coarse radial bin index irc."""
    return (IRneeR(r, a, b, cval) - 1) / irrat + 1

def IVCneeV(v, ivmax, ivrat):
    """Converts physical v = sin(th) into coarse angular bin index ivc."""
    return (IVneeV(v, ivmax) - 1) / ivrat + 1

def RneeIRC(irc, a, b, cval, irrat):
    """Converts coarse radial bin index irc into physical radius (bin center)."""
    return 0.5 * (RneeIR((irc - 1) * irrat + 1, a, b, cval) +
                  RneeIR((irc - 1) * irrat + irrat, a, b, cval))

def VneeIVC(ivc, ivmax, ivrat):
    """Converts coarse angular bin index ivc into physical v = sin(th) (bin center)."""
    return 0.5 * (VneeIV((ivc - 1) * ivrat + 1, ivmax) +
                  VneeIV((ivc - 1) * ivrat + ivrat, ivmax))

def RneeIRhalo(ir, rlgpotmin, rlgpotmax, npot):
    """Converts halo bin index ir into physical radius r (logarithmic grid)."""
    return 10 ** (rlgpotmin + (rlgpotmax - rlgpotmin) / float(npot - 1) * (ir - 1))

def IRneeRhalo(r, rlgpotmin, rlgpotmax, npot):
    """Converts physical halo radius r into halo bin index ir."""
    return round((np.log10(r) - rlgpotmin) /
                 (rlgpotmax - rlgpotmin) * float(npot - 1) + 1.0)

def RneeIRapo(ir, aelz, belz, celz):
    """Converts apocenter bin index ir into physical apocenter radius r."""
    return belz / aelz * (np.exp(aelz * float(ir)) - celz)

def RneeIRperi(ir, aelz, belz, celz):
    """Converts pericenter bin index ir into physical pericenter radius r."""
    return belz / aelz * (np.exp(aelz * float(ir)) - celz)

def IRneeRapo(r, aelz, belz, celz):
    """Converts physical apocenter radius r into apocenter bin index ir."""
    return int(np.log(aelz * r / belz + celz) / aelz + 0.5)

def IRneeRperi(r, aelz, belz, celz):
    """Converts physical pericenter radius r into pericenter bin index ir."""
    return int(np.log(aelz * r / belz + celz) / aelz + 0.5)
