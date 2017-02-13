# Get atmospheric spectra per speices and plot
import spec
from matplotlib.pyplot import *
from numpy import *
import planck
import foregrounds
from astropy import units as u
from astropy import constants as const
from scipy.stats import chi2
from scipy import optimize as opt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from IPython.core.debugger import Tracer; debug_here=Tracer()

fmin = 1150.0
fmax = 2300.0
porder = 4

######################
# Get recombination lines

xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');

# Interpolate to higher resolution, only needed for high R spectrometer
x = linspace(0.8*fmin, 1.2*fmax, 1000000)
#nbin = 100.
#B = ((fmax-fmin)/nbin)*u.GHz
#be = np.arange(fmin, fmax+1e-2, B.value)
#x = be[0:-1]+B.value/2
y = interp(x, xx[:,0], xx[:,1])

r=spec.spec(x,y); 


#########
# Get photon background

Icmb=planck.planck(r.f*u.GHz,2.725*u.K)
cmb=spec.spec(r.f,Icmb.value)
synch45=foregrounds.synch(r.f*u.GHz,b=45)
synch90=foregrounds.synch(r.f*u.GHz,b=90)
dust=foregrounds.dust(r.f*u.GHz)

# Add in a second dust component that is unknown
dust = dust + foregrounds.dust(r.f*u.GHz, T=18.0, tau_f0 = 1e-7)

Itot = spec.spec(r.f, (Icmb + synch45 + dust).value)



##########################
# Straw man experimental configuration to compute errorbars

# Sounding rocket
#t = 8*u.min
#D = 15*u.imperial.inch
#fwhm = 50*u.deg
#nbin = 5

# Space mission
t = 1*u.year
D = 1*u.m
fwhm = 2*u.deg
nbin = 100.0

B = ((fmax-fmin)/nbin)*u.GHz

# Get bins
be = np.arange(fmin, fmax+1e-2, B.value)
bc = be[0:-1]+B.value/2

sigI = zeros(bc.size)
Ibin = zeros(bc.size)
Itotbin = zeros(bc.size)



###############
# Plot

for k,val in enumerate(bc):

    # Frequency
    nu = val*u.GHz

    binind = ((Itot.f > be[k]) & (Itot.f < be[k+1]))

    # Mean background
    I = mean(Itot.I[binind])*u.W/(u.m**2)/u.Hz/u.steradian

    # Number of photons detected
    N = (I * B * t * pi*(D/2)**2 * pi*(fwhm/2)**2 / (const.h * nu) ).decompose()

    # Error on measured I
    sigI[k] = (I / sqrt(N)).value

    # binned recomb signal
    Ibin[k] = mean(r.I[binind])

    # binned background
    Itotbin[k] = I.value

doplot = False
if doplot:
    clf()
    plot(r.f, r.I)
    xlim(fmin,fmax)
    ylim(0,1e-26)
    
    errorbar(bc, Ibin, yerr=sigI, xerr=B.value/2,fmt='.k')
    
    gca().yaxis.set_major_locator(MaxNLocator(4));show()
    grid()
    
    ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)
    xlabel('frequency (GHz)',fontsize=14)
    

A = pi*(D/2)**2
Omega = pi*(fwhm/2)**2
etendu = (A*Omega).to(u.cm**2 * u.steradian)

rthroat = sqrt(etendu / (pi**2 * u.steradian))

snr = sum(Ibin)/sqrt(sum(sigI**2))

fsky = (pi*(fwhm/2)**2).to(u.steradian)/(4*pi*u.steradian)



################
# Fit spectrum
def Ifunc(be,tau_d, T_d, beta_d, T_s, beta_s, T_cmb):

    df = be[1]-be[0]

    # Generate spectrum at higher resolution
    ff = linspace(0.8*1150e9, 1.2*2300e9, 100000)
    
    Id0 = foregrounds.dust(ff, tau_f0=tau_d, T=T_d, beta=beta_d)
    Is0 = foregrounds.synch(ff, b=45, Tgal=T_s, beta=beta_s)
    Icmb0 = planck.planck(ff, T_cmb)
    Itot0 = (Id0 + Is0 + Icmb0).value
    
    # Average in bins
    Itot = zeros(be.size-1)

    for k in range(be.size-1):
        binind = (ff>be[k]) & (ff<be[k+1])
        Itot[k] = mean(Itot0[binind])
    
    return Itot

def Ifunc_chi2(labels, Itot, be):

    pred = Ifunc(be, *labels)
    chi2 = np.sum((pred - Itot)**2)

    return sqrt(chi2)
        
p  = [14.5e-7, 20.5, 1.59, 10.12, -2.55, 2.725]

sg = [14.6e-7, 20.4, 1.57, 10.1, -2.53, 2.727]

lb = [14.4e-7, 20.3, 1.55, 10.08, -2.57, 2.720]
ub = [14.7e-7, 20.7, 1.61, 10.15, -2.50, 2.730]

ranges = [ [lb[k],ub[k]] for k in range(len(lb)) ]

# Data to fit
y = Itotbin

# Brute force
#labels = opt.brute(Ifunc_chi2, ranges,
#                   args=(Itotbin + Ibin, be),
#                   Ns = 6,
#                   finish = None)

# Cleverer non-linear leastsq
labels, covs = opt.curve_fit(Ifunc, be*1e9, y,
                             p0 = sg,
#                             #bounds=(lb,ub),
                             method='lm',
                             maxfev = 15000)

#res = opt.minimize(Ifunc_chi2, sg,
#                   args = (y, be),
#                   method = 'TNC',
#                   bounds = ranges,
#                   options = {'maxiter':5000})


#args = (y, be)
#res = opt.basinhopping(Ifunc_chi2, sg, T = 0.5, stepsize=0.0001,
#                       minimizer_kwargs = {'args':args})
                                                                
