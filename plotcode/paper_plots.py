# Get atmospheric spectra per speices and plot
import am_model as am
import spec
from matplotlib.pyplot import *
from numpy import *
import planck
import foregrounds
import snr
import os
from astropy import units as u
from scipy.stats import chi2
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from IPython.core.debugger import Tracer; debug_here=Tracer()

amdir=os.getenv('HOME')+'/am/'

def paper_plots(fig):

    if fig==1:
        plot_total_foregrounds()

    if fig==2:
        plot_atm_byspecies()

    if fig==3:
        plot_sum_species_err()
        
    return

def plot_total_foregrounds():
    # Plot total foregrounds

    # South Pole, .5 mm PWV
    p=am.readamcfile(amdir+'south_pole/SPole_winter.amc')
    p.Nscale['h2o']=500; p.T0=0
    x1=am.am(prof=p,fmin=0,fmax=1000,df=100)
    x1.callam()
    x2=am.am(prof=p,fmin=1000,fmax=10000,df=1000)
    x2.callam()
    sp=spec.spec(hstack((x1.f,x2.f)), hstack((x1.I,x2.I)))

    # Mauna Kea, 1.86 mm PWV
    p=am.readamcfile(amdir+'mk_no_upperstratosphere.amc')
    p.T0=0
    x1=am.am(prof=p,fmin=0,fmax=1000,df=100)
    x1.callam()
    x2=am.am(prof=p,fmin=1000,fmax=10000,df=1000)
    x2.callam()
    mk=spec.spec(hstack((x1.f,x2.f)), hstack((x1.I,x2.I)))

    # Balloon
    p=am.readamcfile(amdir+'o3_only_stratosphere.amc')
    p.T0=0
    x1=am.am(prof=p,fmin=0,fmax=1000,df=100)
    x1.callam()
    x2=am.am(prof=p,fmin=1000,fmax=10000,df=1000)
    x2.callam()
    ba=spec.spec(hstack((x1.f,x2.f)), hstack((x1.I,x2.I)))

    # Get recombination lines
    xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
    r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance

    # CMB temperature
    I=planck.planck(sp.f*u.GHz,2.725*u.K)
    cmb=spec.spec(sp.f,I.value)

    # Plot total background
    synch15=(foregrounds.synch(sp.f*u.GHz,15)).value
    synch90=(foregrounds.synch(sp.f*u.GHz,90)).value

    # Dust frm Planck 2013 XI, Table 3 and Fig 14
    dust_upper=(foregrounds.dust(sp.f*u.GHz, tau_f0=18.5e-7)).value # Eyeballed off Fig 14,
    dust_lower=(foregrounds.dust(sp.f*u.GHz, tau_f0=6.4e-7)).value # Lowest 1% of sky

    close(1)
    fig=figure(1,figsize=(7,6))    
    loglog(ba.f,ba.I,'-',color=(.9,.9,1),label='Balloon',linewidth=2)
    loglog(mk.f,mk.I,'-',color=(.6,.6,1),label='Mauna Kea, pwv=1.9 mm',linewidth=2)
    loglog(sp.f,sp.I,'-',color=(0,0,1),label='South Pole Winter, pwv=0.5 mm',linewidth=2)
    loglog(cmb.f,cmb.I,'g',label='CMB',linewidth=2)
    fill_between(cmb.f,synch90,synch15,color='c',label=u'galactic synch+ff',linewidth=2)
    fill_between(cmb.f,dust_lower,dust_upper,color='#FF8000',label='galactic dust (fsky=0.01-0.50)',linewidth=2)
    loglog(r.f,r.I,'r',label='recombination lines',linewidth=2)
    xlim(.1,10000);ylim(1e-29,1e-11)
    xlabel('frequency (GHz)');
    ylabel(r'Intensity (W/m$^2$/Hz/sr)')
    legend(prop={'size':10},loc='upper left')
    grid()

    fig.savefig('foregrounds.pdf',format='pdf')



def plot_atm_byspecies():

    p=am.readamcfile(amdir+'o3_thermo.amc');
    p.splitlayers();
    p.Nscale={};
    p.addallnscales()
    p.T0=0

    # Am configuration file
    x=am.am(prof=p,fmin=0,fmax=20,df=0.1)

    # Loop over species
    sp=p.Nscale.keys()

    #########
    # Run am
    s=[]
    for k in sp:
        p.setnscales(0.0);
        p.Nscale[k]=1
        x.callam()
        s.append(spec.spec(x.f,x.I)) # spectral radiance

    # Get recombination lines and polyfit
    xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
    r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance

    ##################################################
    # Plot
    close(2)
    fig = figure(2,figsize=(6,6))

    for k,val in enumerate(s):
        plotind=(val.f>fmin) & (val.f<fmax)
        loglog(val.f[plotind],val.I[plotind],label=sp[k],linewidth=1);

    plotind=(r.f>fmin) & (r.f<fmax)
    loglog(r.f[plotind],r.I[plotind],'--',label='recomb line',linewidth=2,color='r');

    xlabel('frequency (GHz)')
    ylabel('Intensity (W/m$^2$/Hz/sr)')
    grid('on')
    grid(b=True, which='minor')
    xlim(0,10)

    # Order legend from highest to lowest
    ind = where(abs(s[0].f - 2) == min(abs(s[0].f - 2)))[0][0]
    y = array([val.I[ind] for k,val in enumerate(s)])
    ind = flipud(argsort(y))

    ax = gca()
    handles, labels = ax.get_legend_handles_labels()

    h = [handles[val] for k,val in enumerate(ind)]
    h.append(handles[-1])
    l = [labels[val] for k,val in enumerate(ind)]
    l.append(labels[-1])

    legend(h,l,loc='lower right',ncol=2)

    fig.savefig('atm_by_species.pdf',format='pdf')


def plot_sum_species_err():
    # Plot the fractional error of spectrum from sum of species vs. total
    # computation.

    p=am.readamcfile(amdir+'o3_thermo.amc')
    p.T0=0
    x=am.am(prof=p,fmin=0,fmax=20,df=20)
    x.callam()
    stot = spec.spec(x.f,x.I)

    p=am.readamcfile(amdir+'o3_thermo.amc')
    p.T0=0
    p.splitlayers();
    p.Nscale={};
    p.addallnscales()

    # Am configuration file
    x=am.am(prof=p,fmin=0,fmax=20,df=20)

    # Loop over species
    sp=p.Nscale.keys()

    #########
    # Run am
    s=[]
    for k in sp:
        p.setnscales(0.0)
        p.Nscale[k]=1
        x.callam()
        s.append(spec.spec(x.f,x.I)) # spectral radiance

    ssum = spec.spec(x.f, zeros(x.f.shape))
    for k,val in enumerate(s):
        ssum.I = ssum.I + val.I

    # Plot
    close(3)
    fig = figure(3,figsize=(6,6))

    plot(ssum.f, (ssum.I-stot.I)/stot.I, linewidth=2)
    xlabel('frequency (GHz)')
    ylabel(r'$(\Sigma I_{i} - I)/I$')
    grid('on')
    tight_layout()

    fig.savefig('sum_species_err.pdf',format='pdf')


def plot_2to4():

    # Get MK atmosphere
    p=am.readamcfile(amdir+'o3_thermo.amc')
    p.T0=0
    p.splitlayers();
    p.Nscale={};
    p.addallnscales()
    x=am.am(prof=p,fmin=2,fmax=4, df=0.01)
    #########
    # Run am
    s=[]
    for k in sp:
        p.setnscales(0.0)
        p.Nscale[k]=1
        x.callam()
        s.append(spec.spec(x.f,x.I)) # spectral radiance

    mk = spec.spec(x.f, zeros(x.f.shape))
    for k,val in enumerate(s):
        mk.I = mk.I + val.I

    
    # Get recombination lines
    xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
    r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance

    # Polyfit
    porder = 4
    mk.polyfit(order=porder, fmin=2, fmax=4)
    r.polyfit(order=porder, fmin=2, fmax=4)
