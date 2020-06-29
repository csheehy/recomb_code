import numpy as np
import spec
from matplotlib.pyplot import *

#xx=np.loadtxt('../am/south_pole/SPole_winter.amo')
#r=spec.spec(xx[:,0],xx[:,2])

fmin=200
fmax=400
porder=4

xx=np.loadtxt('../am/south_pole/SPole_winter_byspecies.amo')
a=spec.spec(xx[:,0],xx[:,3])
a.fit=a.polyfit(porder,fmin=fmin,fmax=fmax)

xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
fac=(1/100.0)**2 * 1e9
r=spec.spec(xx[:,0],xx[:,1]*fac);
r.fit=r.polyfit(porder,fmin=fmin,fmax=fmax)


clf();
subplot(2,1,1)
semilogx(a.f,a.fit,'--b');
semilogx(r.f,r.fit,'--r');
xlim(fmin,fmax);
ax=gca();ax.set_yscale('log');ax.autoscale(False)
semilogx(a.f,a.I,'b');
semilogx(r.f,r.I,'r')



grid();
grid(which='minor');

#ylim(1e-23,1e-5);



subplot(2,1,2)
semilogx(a.f,a.I-a.fit,'b');
semilogx(r.f,r.I-r.fit,'r')

xlim(fmin,fmax);#ylim(-5e-22,5e-22)
grid();
grid(which='minor');



