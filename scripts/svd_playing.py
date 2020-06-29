import numpy as np
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

datafile = np.load('recomb_plus_atm.npz')
v = datafile['v']
el = datafile['el']
f = datafile['f']
mtype = datafile['mtype']

y1=np.tile(el,(512,1))
y2=np.tile(el**2/100.0,(512,1))

U1,s1,V1=np.linalg.svd(y1,full_matrices=True);
U2,s2,V2=np.linalg.svd(y2,full_matrices=True);
U12,s12,V12=np.linalg.svd(y1+y2,full_matrices=True);

# Now make quadratic component a function of frequency
c=np.tile(np.transpose(np.array(f/4,ndmin=2)),(1,1000))
U12p,s12p,V12p=np.linalg.svd(y1+c*y2,full_matrices=True);

###############################
close(1);figure(1,figsize=(8,8))
subplot(2,1,1)
l1,=plot(y1[0,:],'b',label=r'$\propto$ am');plot(y1[100,:],'.b');
l2,=plot(y2[0,:],'r',label=r'$\propto$ am^2');plot(y2[100,:],'.r')
l3,=plot((y1+y2)[0,:],'g',label=r'$\propto$ am+am$^2$');plot((y1+y2)[100,:],'.g')
l4,=plot((y1+c*y2)[0,:],'c',label=r'$\propto$ am+$f$(freq)*am$^2$');l4p,=plot((y1+c*y2)[100,:],'.c')
xlim(0,200)
leg=legend(handles=[l1,l2,l3,l4],loc='upper right')
gca().add_artist(leg)
legend([l4,l4p],['at freq=2 GHz','at freq=2.4 GHz'],loc='lower right')
xlabel('time');ylabel('Intensity')

subplot(2,1,2)
semilogy(s1[0:9],'b');
semilogy(s2[0:9],'r');
semilogy(s12[0:9],'g');
semilogy(s12p[0:9],'c');
xlabel('SVD component');ylabel('value')

# Atmosphere

datafile_tau = np.load('recomb_plus_atm_tau.npz')
tau0 = datafile_tau['v']

sp='o3'
am=1/np.sin(el*np.pi/180)
#find=0; # 2 GHz
find=256; # 3 Ghz
yy=np.squeeze(v[find,:,np.nonzero(mtype==sp)[0]]);
tau02=np.squeeze(tau0[find,:,np.nonzero(mtype==sp)[0]]);

ind=np.argsort(am)
x=am[ind]
y=yy[ind]
tau=tau02[ind]

I0=np.interp(1.0,x,y)

# Linear fit
plin=np.polyfit(x,y,1)
yfitlin=plin[0]*x+plin[1]

plintau=np.polyfit(x,tau,1)
yfitlintau=plintau[0]*x+plintau[1]


def Ifunc(x, Iz, tauz):
    return Iz*(1-np.exp(-x*tauz))/(1-np.exp(-tauz))
p,pcov=curve_fit(Ifunc,x,y,p0=(I0,tau[0]))
yfit=Ifunc(x,p[0],p[1])

close(2);figure(2,figsize=(14,8))

subplot(2,2,1)
plot(x,tau)
ylabel(r'$\tau$')
title('{0} at 3 GHz'.format(sp))

subplot(2,2,3)
r=tau-yfitlintau;
plot(x,r,'.',label='linear fit resids, rms={:0.2g}'.format(np.std(r)));
legend(loc='lower right')
xlabel('airmass')
ylabel('residuals')


subplot(2,2,2)
plot(x,y)
ylabel('Intensity')
title('{0} at 3 GHz'.format(sp))

ax=subplot(2,2,4)
rlin=y-yfitlin
rexp=y-yfit
plot(x,rlin,'.b',label=r'linear fit resids, rms={:0.2g}'.format(np.std(rlin)))
plot(x,rexp,'.g',label=r'$f=I_z(1-e^{{-\tau_z am}})/(1-e^{{-\tau_z}})$ fit resids, rms={:0.2g}'.format(np.std(rexp)))
legend(loc='lower right')
xlabel('airmass')
ylabel('residuals')
text(0.95,.9,r'best fit $\tau_z=${:0.4g}'.format(p[1]),transform=ax.transAxes,ha='right')
