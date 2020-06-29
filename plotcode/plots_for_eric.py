import numpy as np
from matplotlib.pyplot import *

# Load instance of time-ordered-data object
d=np.load('recomb_plus_atm.npz')
mtype=d['mtype']
v=d['v']
t=d['t']
f=d['f']
el=d['el']


# plot elevation vs. time
figure(1)
clf();plot(t,el);xlabel('t (sec)');ylabel('el (deg)');
show()

# plot a few model components, (2 atmospheric species, a 2.7 K planck
# function, and the recombination lines)
figure(2)
clf()
ext=[t[0],t[-1],f[-1],f[0]]
k=[0,5,12,6]

for l,kk in enumerate(k):
    subplot(np.size(k),1,l+1)
    imshow(v[:,:,kk],extent=ext,aspect='auto');ylabel('f (GHz)');
    title((mtype[kk]) + ' (W/m^2/Hz/sr)')
    if l<3:
        gca().set_xticks([])
    else:
        xlabel('time (s)')
    colorbar()
show()

# sum atmospheric components.  Why oh why are dictionaries unordered?! This
# deserves some dedicated subfunctions to manipulate.
atm=np.sum(v[:,:,[0,1,2,3,4,5,7,8,9,10,11]],axis=2)
figure(3)
clf()
imshow(atm,extent=ext,aspect='auto');title('sum of atmsopheric components');
colorbar()
show()


