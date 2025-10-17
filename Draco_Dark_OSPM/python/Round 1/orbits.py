#orbits.py

from __future__ import annotations
import os
import math
from pathlib import Path
from datetime import datetime
import numpy as np
from numpy.polynomial.legendre import legval
from dataclasses import dataclass
from astropy.coordinates import SkyCoord 
import astropy.units as u
import pandas as pd
#THE ONE TRUE FILE
from libshim import register_halo_grid, register_density_grid, fac1
import tables_state as S
from tables import tables, tables_read, tables_write  # etc.
from potential import potential, force
from io_utils import iterread, mlread, weight_read
from orbits import run_single_orbit


# --- energy.f translation ---
def energy(r: float, th: float, vr: float, vth: float, Lz: float) -> float:
    v = math.sin(th)           # Fortran v = sin(th)
    vco = math.cos(th)         # vco = cos(th)
    vph = Lz / (r * vco)       # vphi = xLz/r/vco
    VV  = potential(r, v)      # call potential(r, v, VV)
    return 0.5 * (vr*vr + vth*vth + vph*vph) + VV

def make_ekin(Etot: float, Lz: float):
    def _ekin(r: float, v: float) -> float:
        vv = max(-1.0 + 1e-12, min(1.0 - 1e-12, v))
        VV = potential(r, vv)
        return Etot - VV - (Lz*Lz) / (2.0 * r*r * (1.0 - vv*vv))
    return _ekin

# =================== integrator core (derivs uses force) ===================
def derivs(state: State, Lz: float) -> State:
    r, th, vr, vth = state.r, state.th, state.vr, state.vth
    v   = math.sin(th)                 # Fortran v = sin(th)
    cth = math.cos(th)
    vph = Lz / (r * cth)               # vphi = xLz/r/cos(th)
    frL, fvL = force(r, v)             # matches: call force(r,v,frL,fvL)
    drdt   = vr
    dthdt  = vth / r
    dvrdt  = (vth*vth + vph*vph)/r + frL
    dvthdt = -(math.tan(th) * vph*vph)/r - (vth*vr)/r + fvL
    return State(drdt, dthdt, dvrdt, dvthdt)

def rk4_step(state: State, dt: float, Lz: float) -> State:
    def add(a: State, b: State, f: float) -> State:
        return State(a.r + f*b.r, a.th + f*b.th, a.vr + f*b.vr, a.vth + f*b.vth)
    k1 = derivs(state, Lz)
    k2 = derivs(add(state, k1, dt*0.5), Lz)
    k3 = derivs(add(state, k2, dt*0.5), Lz)
    k4 = derivs(add(state, k3, dt), Lz)
    return State(
        state.r   + dt*(k1.r   + 2*k2.r   + 2*k3.r   + k4.r  )/6.0,
        state.th  + dt*(k1.th  + 2*k2.th  + 2*k3.th  + k4.th )/6.0,
        state.vr  + dt*(k1.vr  + 2*k2.vr  + 2*k3.vr  + k4.vr )/6.0,
        state.vth + dt*(k1.vth + 2*k2.vth + 2*k3.vth + k4.vth)/6.0,
    )
def step(state: State, Lz: float, epsilon: float):
    dt = epsilon
    return rk4_step(state, dt, Lz), dt