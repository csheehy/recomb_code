# Get atmospheric spectra per speices and plot
import am_model as am
import spec
from matplotlib.pyplot import *
import planck
import foregrounds
import snr
from astropy import units as u
from scipy.stats import chi2
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from IPython.core.debugger import Tracer; debug_here=Tracer()


getspec=True

if getspec:
    # Get atmospheric profile where I've taken the standard South Pole winter
    # profile and broken it down into constituent species (i.e. "dry_air" -> "o2" +
    # "o2_air" + "co2" etc. This is done by hand but could be automated.
    p=am.readamcfile('../am/south_pole/SPole_winter.amc');p.pwv=500
    #p=am.readamcfile('../am/chajnantor/Chajnantor.amc');p.pwv=500
    #p=am.readamcfile('../am/generic/generic_mid.amc');p.pwv=1
    p.splitlayers();
    p.Nscale={};p.addallnscales()
    
    # Am configuration file
    x=am.am(prof=p,fmin=0,fmax=3000,df=20)
    
    # Loop over species
    sp=p.Nscale.keys()
    
    #########
    # Run am
    s=[]
    for k in sp:
        p.setnscales(0.0);
        
        if k.startswith('h2o'):
            p.Nscale[k]=p.pwv
        else:
            p.Nscale[k]=1
            
        x.callam()
        x.parseresults()
        s.append(spec.spec(x.f,x.I)) # spectral radiance

# Nscale value
Nsc=.0005

########
# Polyfit
porder=4 # 4th order
fmin=9
fmax=12
for k,val in enumerate(s):
    s[k].yfit=s[k].polyfit(order=porder,fmin=fmin,fmax=fmax)

#######
# Get recombination lines and polyfit
xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance
r.yfit=r.polyfit(order=porder,fmin=fmin,fmax=fmax)

xx=np.loadtxt('recomb_spec/HI.dat');
rH=spec.spec(xx[:,0],xx[:,1]); # spectral radiance
rH.yfit=rH.polyfit(order=porder,fmin=fmin,fmax=fmax)


# CMB temperature and foregrounds
I=planck.planck(s[0].f*u.GHz,2.725*u.K)
cmb=spec.spec(s[0].f,I.value)

Itot=np.zeros(np.shape(s[0].I))
for k,val in enumerate(s):
    Itot+=val.I

synch45=foregrounds.synch(cmb.f*u.GHz,45)
synch90=foregrounds.synch(cmb.f*u.GHz,90)
dust=foregrounds.dust(cmb.f*u.GHz)

# Add snr calc
Nmodes=2
dt=90*u.day
df=0.01*u.GHz

# Antenna temperature (in intensity units)
Iback=Itot+(synch90.to(u.W/u.Hz/u.m**2/u.sr)).value+cmb.I

# Noise factor varies linearly with frequency
# c.f. documentation for ATF-35143 amplifier
#Nfy=[0,.05] # Cryogenic
Nfy=[0,.2] # Off shelf -40 C
fx=[0,4] # in GHz

# Get bins
be=np.arange(fmin,fmax+1e-2,df.value)
bc=be[0:-1]+df.value/2

fig1=figure(1,figsize=(12,8))
clf()

ytot=np.zeros(np.shape(s[0].f))
for k,val in enumerate(s):
    if k>4:
        lt='--'
    else:
        lt=':'
    plotind=(val.f>fmin) & (val.f<fmax)
    yy=Nsc*(val.I-val.yfit)
    plot(val.f[plotind],yy[plotind],linewidth=2,label=sp[k],linestyle=lt)
    ytot+=yy

plotind=(r.f>fmin) & (r.f<fmax)
plot(r.f[plotind],r.I[plotind]-r.yfit[plotind],'r',label='recomb lines')
ytot+=np.interp(s[0].f,r.f,r.I-r.yfit)
#plot(s[0].f[plotind],ytot[plotind],'k',linewidth=3,label='total')


xlim(fmin,fmax)
legend(prop={'size':12},ncol=2,loc='upper left')
grid()
title(u'baseline subtracted spectra, atm. components x {:0.1e} (2 modes, -40\N{DEGREE SIGN} C amplifiers)'.format(Nsc),fontsize=12)
xlabel('frequency (GHz)',fontsize=14)
ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)



# Points with error bars
Nf=np.interp(bc,fx,Nfy)
Tn=snr.Nf2Tn(Nf)
In=planck.Ta2I(bc*u.GHz,Tn)

# Interpolate background intensity
Iback2=np.interp(bc,s[0].f,Iback)*u.W/u.m**2/u.Hz/u.sr

Isys=Iback2+In

dI=snr.radiometer_equation(Isys,dt,df)/np.sqrt(Nmodes)

dIback=snr.radiometer_equation(Iback2,dt,df)/np.sqrt(Nmodes)

# Average spectrum in bins
y=np.zeros(np.shape(bc))
    
for k,val in enumerate(bc):
    avgind=(s[0].f>be[k]) & (s[0].f<be[k+1])
    y[k]=np.nanmean(ytot[avgind])
    
errorbar(bc,y,yerr=dI.value,xerr=df.value/2,fmt='.k')
errorbar(bc,y,yerr=dIback.value,fmt='.k')

dosave=True
if dosave:
    fig1.savefig('atm_monitor.pdf',format='pdf')

