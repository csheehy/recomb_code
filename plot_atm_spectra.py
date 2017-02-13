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

# Get atmospheric profile where I've taken the standard South Pole winter
# profile and broken it down into constituent species (i.e. "dry_air" -> "o2" +
# "o2_air" + "co2" etc. This is done by hand but could be automated.
amdir='/Users/csheehy/am/'
p=am.readamcfile(amdir+'south_pole/SPole_winter.amc');p.pwv=500
#p=am.readamcfile(amdir+'chajnantor/Chajnantor.amc');p.pwv=500
#p=am.readamcfile(amdir+'generic/generic_mid.amc');p.pwv=1
p.splitlayers();
p.Nscale={};p.addallnscales()
p.T0=0

onlylines=True
if onlylines:
    rmsp=['o2air','o2_uncoupled','o2_coupled','h2o_lines','h2o_continuum']
    for k,val in enumerate(rmsp):
        p.Nscale.pop(val)
        for j,val2 in enumerate(p.layers):
            c=val2.col
            if c.has_key(val):
                c.pop(val)

        
# Am configuration file
x=am.am(prof=p,fmin=0,fmax=3000,df=20)

# Loop over species
sp=p.Nscale.keys()

# This is zero everywhere
#sp.remove('h2o_optical_refractivity')

# Nscale value
Nsc=.01

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
    s.append(spec.spec(x.f,x.I)) # spectral radiance

########
# Polyfit
porder=4 # 4th order
fmin=1.5
fmax=2.4
if fmin==np.round(fmin):
    xt=np.arange(fmin,fmax+1)
else:
    xt=np.linspace(fmin,fmax,2)
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


# CMB temperature
I=planck.planck(s[0].f*u.GHz,2.725*u.K)
cmb=spec.spec(s[0].f,I.value)




##################################################
# Plot
#close(1)

fig1=figure(1,figsize=(12,8))
clf()

for k,val in enumerate(s):
    subplot(3,4,k+1)
    plotind=(val.f>fmin) & (val.f<fmax)
    h1,=plot(val.f[plotind],val.I[plotind],label=sp[k],linewidth=2);
    xlim(fmin,fmax)
    #gca().set_xticks(xt)
    h2,=plot(val.f,val.yfit,'--',linewidth=2)
    if sp[k]=='ch4':
        ylim(0,3e-28)
    gca().yaxis.set_major_locator(MaxNLocator(3));show()
    gca().set_xticks(xt)
    #legend([h1,h2],[sp[k],'best fit p{0}'.format(porder)],prop={'size':10},loc="upper left")
    legend([h1],[sp[k]],prop={'size':10},loc="upper left")
    grid('on')

subplot(3,4,12)
plotind=(r.f>fmin) & (r.f<fmax)
h1,=plot(r.f[plotind],r.I[plotind],label='recomb line',linewidth=2,color='r');
xlim(fmin,fmax)
#gca().set_xticks(xt)
h2,=plot(r.f[plotind],r.yfit[plotind],'--g',linewidth=2);
#legend([h1,h2],['recomb lines','best fit p{0}'.format(porder)],prop={'size':10},loc="upper left")
gca().yaxis.set_major_locator(MaxNLocator(3));show()
gca().set_xticks(xt)
legend([h1],['recomb lines'],prop={'size':10},loc="upper left")
grid('on')
figtext(0.5,0.04,'frequency (GHz)',{'ha':'center'},fontsize=18);
figtext(0.05,0.5,r'Intensity (W/m$^2$/Hz/sr)',{'va':'center'},rotation=90,fontsize=18);



#close(2)
fig2=figure(2,figsize=(12,8))
clf()

for k,val in enumerate(s):
    subplot(3,4,k+1)
    plotind=(val.f>fmin) & (val.f<fmax)
    h1,=plot(val.f[plotind],val.I[plotind]-val.yfit[plotind],
             label='{0} - best fit'.format(sp[k]),linewidth=2);
    xlim(fmin,fmax)
    gca().set_xticks(xt)
    gca().yaxis.set_major_locator(MaxNLocator(3));show()
    legend(prop={'size':10},loc="upper left")
    grid('on')

subplot(3,4,12)
plotind=(r.f>fmin) & (r.f<fmax)
h1,=plot(r.f[plotind],r.I[plotind]-r.yfit[plotind],
         label='recomb line - best fit',linewidth=2,color='r');

xlim(fmin,fmax)
gca().set_xticks(xt)
gca().yaxis.set_major_locator(MaxNLocator(3));show()
legend(prop={'size':10},loc="upper left")
grid('on')
figtext(0.5,0.04,'frequency (GHz)',{'ha':'center'},fontsize=18);
figtext(0.05,0.5,r'Intensity (W/m$^2$/Hz/sr)',{'va':'center'},rotation=90,fontsize=18);



fig3=figure(3,figsize=(12,8))
clf()

for k,val in enumerate(s):
    if k>4:
        lt='--'
    else:
        lt=':'
    plotind=(val.f>fmin) & (val.f<fmax)
    plot(val.f[plotind],Nsc*(val.I[plotind]-val.yfit[plotind]),linewidth=2,label=sp[k],linestyle=lt)
plotind=(r.f>fmin) & (r.f<fmax)
plot(r.f[plotind],r.I[plotind]-r.yfit[plotind],'r',label='recomb lines',linewidth=3)
xlim(fmin,fmax)
legend(prop={'size':12},ncol=2,loc='upper left')
grid()
title('baseline subtracted spectra, atm. components x {:0.1e}'.format(Nsc),fontsize=12)
xlabel('frequency (GHz)',fontsize=14)
ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)


# Plot total background
fig4=figure(4,figsize=(7,6))
clf()

Itot=s[0].I
for k,val in enumerate(s[1:]):
    Itot+=val.I

synch45=foregrounds.synch(cmb.f*u.GHz,45)
synch90=foregrounds.synch(cmb.f*u.GHz,90)
dust=foregrounds.dust(cmb.f*u.GHz)

loglog(s[0].f,Itot,label='South Pole Winter, pwv=0.5 mm',linewidth=2)
loglog(cmb.f,cmb.I,'g',label='CMB',linewidth=2)
loglog(cmb.f,synch45,'c--',label=u'synch+ff (b=45\N{DEGREE SIGN})',linewidth=2)
loglog(cmb.f,synch90,'c',label=u'synch+ff (b=90\N{DEGREE SIGN})',linewidth=2)
loglog(cmb.f,dust,c='#FF8000',label='thermal dust (southern cap)',linewidth=2)
loglog(r.f,r.I,'r',label='recombination lines',linewidth=2)
xlim(.1,3000);ylim(1e-29,1e-11)
xlabel('frequency (GHz)');
ylabel(r'Intensity (W/m$^2$/Hz/sr)')
legend(prop={'size':10},loc='upper left')
grid()


##########################
# Add snr calc
dt=1.0*u.year
nbin=20.0
df=((fmax-fmin)/nbin)*u.GHz

# Antenna temperature (in intensity units)
Iback=Itot+(synch90.to(u.W/u.Hz/u.m**2/u.sr)).value+cmb.I

# Noise factor varies linearly with frequency
# c.f. documentation for ATF-35143 amplifier
#Nfy=[0,.05] # Cryogenic
#Nfy=[0,.2] # Off shelf -40 C
#Nfy_arr=[ [0,0.2],[0,0.05],[0,0.05] ]
Nmodes_arr=[128,1000,10000]
Tarr=np.array([ [1,1],[4,4],[4,4] ])
Nfy_arr=snr.Tn2Nf(Tarr)
fx=[0,4] # in GHz

# Get bins
be=np.arange(fmin,fmax+1e-2,df.value)
bc=be[0:-1]+df.value/2

fig5=figure(5,figsize=(10,7))
clf()


for kk,val in enumerate(Nmodes_arr):


    
    #ax=subplot(1,3,kk+1)
    figure()

    Nmodes=Nmodes_arr[kk]
    Nfy=Nfy_arr[kk]

    Nf=np.interp(bc,fx,Nfy)
    Tn=snr.Nf2Tn(Nf)
    In=planck.Ta2I(bc*u.GHz,Tn)

    # Interpolate background intensity
    Iback2=np.interp(bc,s[0].f,Iback)*u.W/u.m**2/u.Hz/u.sr

    Isys=Iback2+In

    dI=snr.radiometer_equation(Isys,dt,df)/np.sqrt(Nmodes)

    dIback=snr.radiometer_equation(Iback2,dt,df)/np.sqrt(Nmodes)

    plotind=(rH.f>fmin) & (rH.f<fmax)
    plot(rH.f[plotind],rH.I[plotind]-rH.yfit[plotind],'c',label='H recomb',linewidth=3)
    plotind=(r.f>fmin) & (r.f<fmax)
    plot(r.f[plotind],r.I[plotind]-r.yfit[plotind],'r',label='H+He recomb',linewidth=3)
    xlim(fmin,fmax)

    if kk==0:
        legend(prop={'size':12},ncol=2,loc='upper left')
    
    # Average spectrum in bins
    y=np.zeros(np.shape(bc))
    
    for k,val in enumerate(bc):
        avgind=(r.f>be[k]) & (r.f<be[k+1])
        y[k]=np.nanmean(r.I[avgind]-r.yfit[avgind])
        
    errorbar(bc,y,yerr=dI.value,xerr=df.value/2,fmt='.k')
    errorbar(bc,y,yerr=dIback.value,fmt='.k')

    Tsys=planck.I2Ta(bc*1e9,In).value
    tit='{:d} modes, Tsys = {:0.1f} K, {:0.1f} integration'.format(Nmodes_arr[kk],np.mean(Tsys),dt)
    title(tit)
    
    gca().yaxis.set_major_locator(MaxNLocator(4));show()
    grid()


    ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)
    xlabel('frequency (GHz)',fontsize=14)


    
    # Quick monte carlo
    if kk==0:
        #yl=ylim()
        yl=(-6e-29,6e-29)
        err=dI.value
        w=y/err**2 # weighted average using s/var weighting
        #w=y
        
        N=50000
        s_nosig=np.zeros(N)
        s_withsig=np.sum(y*w)/np.sum(w)
        
        for k in range(N):
            yy=np.random.standard_normal(bc.shape)*err
            s_nosig[k]=np.sum(yy*w)/np.sum(w)
            
            pp=np.count_nonzero(s_nosig<s_withsig)/float(N)
            sigma=np.sqrt(chi2.ppf(pp,1))
            
        at=AnchoredText(r'{:02.1f}$\sigma$ detection'.format(sigma),
                        prop=dict(size=10), loc=3)
        gca().add_artist(at)
        draw()

    ylim(yl)




# Print figures
dosave=False

if dosave:
    fig1.savefig('atm_total.pdf',format='pdf')
    fig2.savefig('atm_psub.pdf',format='pdf')
    fig3.savefig('atm_psub_scaled.pdf',format='pdf')
    fig4.savefig('foregrounds.pdf',format='pdf')
    fig5.savefig('snr.pdf',format='pdf')





