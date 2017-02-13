# Foregrounds
import numpy as np
import planck
from astropy.constants import h, k_B, c
from astropy import units as u
from IPython.core.debugger import Tracer; debug_here=Tracer()

def synch(f, b=45, Tgal=10.12, beta=-2.55):
    """Get Arcade2 best fit galactic component from Kogut et al, 2011. If f is
    not an astropy units quantity it is assumed in Hz. Galactic latitude b in
    degrees."""
    
    if type(f) is not u.quantity.Quantity:
        f=f*u.Hz
        
    f0=0.31*u.GHz

    Ta = Tgal*u.K*(f/f0)**beta / np.sin(np.deg2rad(np.abs(b))) # Antenna temperature
    Ta = Ta.decompose()
    
    # Convert to inensity
    I=planck.Ta2I(f,Ta)

    return I
    
def dust(f, tau_f0=14.5e-7, T=20.5, beta=1.59):
    """Planck "southern cap" sky model dust."""
    
    if type(f) is not u.quantity.Quantity:
        f=f*u.Hz

    f0=353*u.GHz
    I = tau_f0*planck.planck(f,T*u.K)*((f/f0).decompose().value)**beta

    return I.to(u.W/u.m**2/u.Hz/u.steradian)
 
