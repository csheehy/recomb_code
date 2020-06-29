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

########
# Polyfit
porder=4 # 4th order
fmin=2
fmax=4
xt=[2,3,4]

#######
# Get recombination lines and polyfit
xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance

figure(1)
clf()
r.yfit=r.polyfit(order=porder,fmin=fmin,fmax=fmax)
plotind=(r.f>fmin) & (r.f<fmax)
plot(r.f[plotind],r.I[plotind],'r',label='recombination lines',linewidth=3)
#plot(r.f[plotind],r.yfit[plotind],'r--',label='best fit 4th order polynomial ',linewidth=3)
ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)
xlabel('frequency (GHz)',fontsize=14)
legend(prop={'size':12},ncol=2,loc='upper left')
grid('on')

figure(2)
clf()
r.yfit=r.polyfit(order=porder,fmin=fmin,fmax=fmax)
plotind=(r.f>fmin) & (r.f<fmax)
plot(r.f[plotind],r.I[plotind],'r',label='recombination lines',linewidth=3)
plot(r.f[plotind],r.yfit[plotind],'r--',label='best fit 4th order polynomial ',linewidth=3)
ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)
xlabel('frequency (GHz)',fontsize=14)
legend(prop={'size':12},ncol=2,loc='upper left')
grid('on')


figure(3)
clf()
r.yfit=r.polyfit(order=porder,fmin=fmin,fmax=fmax)
plotind=(r.f>fmin) & (r.f<fmax)
plot(r.f[plotind],r.I[plotind]-r.yfit[plotind],'r',label='recombination lines',linewidth=3)
ylabel(r'Intensity (W/m$^2$/Hz/sr)',fontsize=14)
xlabel('frequency (GHz)',fontsize=14)
legend(prop={'size':12},ncol=2,loc='upper left')
ylim(-6e-29,6e-29)
ax=gca();
ax.set_yticks(np.array([-4,-2,0,2,4])*1e-29)
grid('on')
