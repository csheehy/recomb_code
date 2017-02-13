# Signal to noise calculations
import numpy as np
from astropy.constants import h, k_B, c
from astropy import units as u
from IPython.core.debugger import Tracer; debug_here=Tracer()

def radiometer_equation(Tsys,dt,df):
    """Return sigma_T for input Tsys, integration time dt, and bandwidth
    df. Tsys, dt, df assumed in K, s, and Hz if not astropy unit quantities.
    """

    if type(Tsys) is not u.quantity.Quantity:
        Tsys=Tsys*u.K

    if type(dt) is not u.quantity.Quantity:
        dt=dt*u.s

    if type(df) is not u.quantity.Quantity:
        df=df*u.Hz

    sigma_T = Tsys / np.sqrt(dt * df)

    return sigma_T.to(Tsys.unit)

def Nf2Tn(Nf):
    """Noise factor to noise temperature."""
    
    Tref = 290
    Tn = Tref * (10.0**(Nf/10.0)-1)

    return Tn

def Tn2Nf(Tn):
    """Noise temp to noise factor"""
    
    Tref = 290.0
    Nf = 10.0*np.log10((Tn/Tref)+1)
    
    return Nf


