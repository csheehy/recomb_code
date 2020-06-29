import sim
import am_model as am
import numpy as np
from matplotlib.pyplot import *


def doplot(prof,sp):

    mod=am.am(prof=prof,df=20,fmin=2,fmax=4)
    zaarr=np.linspace(0,45,100)
    mod.callam()
    mod.parseresults()
    
    Iarr=np.zeros((mod.f.size,zaarr.size))
    tauarr=np.zeros((mod.f.size,zaarr.size))
    
    for k,val in enumerate(zaarr):
        
        mod.za=val
        mod.callam()
        mod.parseresults()
        Iarr[:,k]=mod.I
        tauarr[:,k]=mod.tau
        
    # Secant model
    amm=1/np.cos(zaarr*np.pi/180.)
    
    # Fit line
    find=50
    
    x=tauarr[find,:]
    y=Iarr[find,:]
    
    #pp=np.polyfit(x,y,2)
    #fit=pp[0]*x**2+pp[1]*x+pp[2]
    pp=np.polyfit(x,y,1)
    fit=pp[0]*x+pp[1]
    
    
    clf()
    subplot(2,2,1)
    plot(x,y,'.',label='am output')
    plot(x,fit,label='linear fit')
    xlabel('tau')
    ylabel('Intensity');
    fstr=np.str(mod.f[find])
    title(sp + ', ' + fstr + ' GHz')
    legend(loc='lower right')
    
    subplot(2,2,3)
    plot(x,y-fit,'.')
    ylabel('residuals')
    xlabel('tau')
    
    
    
    x=amm
    y=tauarr[find,:]
    
    pp=np.polyfit(x,y,1)
    fit=pp[0]*x+pp[1]
    
    subplot(2,2,2)
    plot(x,y,'.',label='am output')
    plot(x,fit,label='linear fit')
    ylabel('tau');
    xlabel('airmass (sec model)')
    legend(loc='lower right')
    
    subplot(2,2,4)
    plot(x,y-fit,'.')
    ylabel('residuals')
    xlabel('airmass (sec model)');



###########################
m=sim.skymodel()
    

# o2
close(1);figure(1,figsize=(12,7))
clf()
sp='o2_coupled'
prof=m.model[sp]
doplot(prof,'o2_coupled')
savefig('fig1.png')

# h2o
close(2);figure(2,figsize=(12,7))
clf()
sp='h2o_lines'
prof=m.model[sp]
prof.Nscale['h2o_lines']=500
doplot(prof,'h2o_lines at .5 mm PWV')
savefig('fig2.png')

# h2o matched to o2 intensity
close(3);figure(3,figsize=(12,7))
clf()
prof.Nscale['h2o_lines']=500e3
doplot(prof,'h2o_lines at 500 mm PWV')
savefig('fig3.png')

    
    





