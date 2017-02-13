# Simulation of radiometer data
import numpy as np
import am_model as am
import planck
import spec
import foregrounds
import tod
import copy
import pickle
from scipy.interpolate import interp1d
from matplotlib.pyplot import *
from IPython.core.debugger import Tracer; debug_here=Tracer()

def runsim():

    tt=traj(tint=100);
    p=tt.gettraj()
    d=tod.d(p=p)

    m=skymodel()
    d=m.gensig(d)
    
    #d.save('recomb_plus_atm.npz')
    
    return d

def za2am(za, type=None):
    """Compute airmass for given input tau. Zenith angle in degrees.
    type = 'secant' (default), 'Young'"""

    if type is None:
        type='secant'

    if type == 'secant':
        am=1/np.cos(za*np.pi/180)

    if type == 'Young':
        print 'Not coded yet'

    return am

def Iz2Iam(Iz,tauz,am):
    """Compute intensity as a function of airmass given zenith intensity and
    zenith optical depth."""
    I = Iz*(1-np.exp(-am*tauz))/(1-np.exp(-tauz))
    return I

class traj:

    def __init__(self,dt=None,tint=None,v=None,elrange=None,acc=None):
        """Pointing trajectory of telescope
        dt is sample rate in seconds
        tint = integration time in seconds
        v = velocity in constant velocity section in deg/s
        elrange = elevation range of constant velocity scan
        acc = max acceleration in deg/s^2"""

        # Sample rate in seconds
        if dt is None:
            self.dt=0.1
        else:
            self.dt=dt

        # Integration time in s
        if tint is None:
            self.tint = 100.0
        else:
            self.tint= tint

        if v is None:
            self.v = 10.0
        else:
            self.v = v
        
        # Constant velocity elevation range
        if elrange is None:
            self.elrange = np.array([50.0,85.0])
        else:
            self.elrange = elrange

        # Max acceleration in deg/s^2
        if acc is None:
            self.acc = 15
        else:
            self.acc = acc


    def gettraj(self):
        """Calculate trajectory given pointing parameters, return pointing and
        time classes. Use cosine acceleration to constant velocity scans."""
        
        # Time array
        t=np.arange(0,self.tint,self.dt)

        # Turnaround trajectory
        xsin = np.linspace(0,np.pi,1000)
        cc = self.acc/self.v
        tsin = xsin/cc
        A=self.v/cc
        ysin = A*np.sin(xsin)

        # Total time = upward const + turn + downward const + turn
        tconst = np.diff(self.elrange)/self.v        
        tturn = tsin[-1]-tsin[0]
        ttot = 2*(tconst+tturn)
        
        ################################
        # Calculate trajectory
        tmod = np.mod(t,ttot)
        el = np.zeros_like(t)
        el.fill(np.nan)

        ind=(tmod>=0) & (tmod<tconst) # Upward constant velocty
        el[ind]=self.elrange[0]+tmod[ind]*self.v

        ind=((tmod>=tconst+tturn) & (tmod<2*tconst+tturn)).nonzero()[0] # Downward constant velocity
        if ind.any():
            el[ind]=self.elrange[1]-(tmod[ind]-tmod[ind[0]])*self.v

        ind0=np.isnan(el) # Turnarounds
        ind=(ind0 & (tmod<=0.5*ttot)).nonzero()[0]
        if ind.any():
            el[ind]=np.interp(tmod[ind],tsin+tmod[ind[0]],ysin+self.elrange[1])

        ind=np.isnan(el).nonzero()[0]
        if ind.any():
            el[ind]=np.interp(tmod[ind],tsin+tmod[ind[0]],-ysin+self.elrange[0])

        point=tod.point(0,el,t)
        
        return point

    

class skymodel:
    
    def __init__(self,model=None):
        """Generate sky model. Models must output I(f;d)"""
        
        if model is None:

            self.model={}

            # Atmospheric model
            p=am.readamcfile('../am/south_pole/SPole_winter.amc')

            # Total atm signal
            #p.Nscale['h2o']=500 # 500 um pwv
            #self.model['atm']=p

            # Per species
            p.splitlayers();
            p.Nscale={};
            p.addallnscales()
        
            # Loop over species
            sp=p.Nscale.keys()
    
            for k in sp:
                p.setnscales(0.0);
                
                if k.startswith('h2o'):
                    p.Nscale[k]=500.0 # 500 um PWV
                else:
                    p.Nscale[k]=1
                
                self.model[k]=copy.deepcopy(p)

            # Recomb signal
            xx=np.loadtxt('recomb_spec/HI.HeI.HeII.dat');
            r=spec.spec(xx[:,0],xx[:,1]); # spectral radiance
            self.model['sig']=r

            # CMB
            self.model['planck']=2.725

            # Synchrotron
            self.model['sync']='LFI_SkyMap_030_1024_R1.10_nominal.fits'

        else:
            self.model=model


    def gensig(self,d):
        """Return simulated spectra"""

        nf=d.data.getnf()
        nt=d.point.getnt()

        for k,mtype in enumerate(self.model):

            print mtype

            if isinstance(self.model[mtype],am.profile):
                
                v=np.zeros((nf,nt))

                p=self.model[mtype]
                df=d.data.getdf()*1000.0 # GHz -> MHz
                fmin=d.data.f[0]-df/1000
                fmax=d.data.f[-1]
                za=np.abs(90-d.point.el)
                m=am.am(prof=p,df=df,fmin=fmin,fmax=fmax,za=0)
                
                ## Get spectra for a number of zenith angles
                #zaarr=np.linspace(np.min(za),np.max(za),500)
                #vmod=np.zeros((nf,np.size(zaarr)))
                #for l,z in enumerate(zaarr):
                #    m.za=z
                #    m.callam()
                #    m.parseresults()
                #    vmod[:,l]=np.interp(d.data.f,m.f,m.I)
                ## Interpolate to simulated data
                #for l,f in enumerate(d.data.f):
                #    interpf=interp1d(zaarr,vmod[l,:],kind='cubic')
                #    v[l,:]=interpf(za)

                # I have verified that am scales tau exactly linearly with
                # airmass and assumes a secant airmass model.
                # Get spectrum for zenith and scale to airmass.
                m.callam()
                m.parseresults()
                airmass=za2am(za,'secant')
                for l,aml in enumerate(airmass):
                    v[:,l]=Iz2Iam(m.I,m.tau,aml)

                d.data.append(v,mtype)

            if mtype == 'sig':
                # Recombination line signal
                m=self.model[mtype]
                y=np.interp(d.data.f,m.f,m.I)
                y=np.reshape(y,[y.size,1])
                v=np.tile(y,[1,nt])
                
                d.data.append(v,mtype)

            if mtype == 'planck':
                # Thermal CMB
                I=(planck.planck(d.data.f*1e9, 2.725)).value
                v=np.tile(I.reshape(I.shape+(1,)),(1,nt))
                d.data.append(v,mtype)

            if hasattr(self.model[mtype],'endswith'):
                if self.model[mtype].endswith('.fits'):
                    # Input sky map, assume in galactic coords
                    fname=self.model[mtype]
                    print(fname)
                    
                
        # Record model
        d.m=self.model

        return d
            
