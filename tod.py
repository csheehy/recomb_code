# Time ordered data
import numpy as np
from IPython.core.debugger import Tracer; debug_here=Tracer()

class point:

    def __init__(self,az=None,el=None,t=None):
        """Telescope pointing: 
        azimuth and elevation in degrees
        t in seconds"""

        self.az=az
        self.el=el        
        self.t=t

        # This is only until I can figure out how to compute az,el -> ra,dec
        self.ra = az + np.mod(t, 240*3600 - 4.0*60)*(360.0/(24*3600))
        self.dec = -el

    def getdt(self):
        """Time resolution"""
        return self.t[1]-self.t[0]

    def getnt(self):
        """Number of time elements"""
        return np.size(self.t)
        
class data:

    def __init__(self,v=None,f=None,mtype=None):
        """Data. f in GHz"""

        if v is None:
            self.v=[]
        else:
            self.v=v

        if f is None:
            self.f=np.linspace(2,4,512)
        else:
            self.f=f

        if mtype is None:
            self.mtype=[]
        else:
            self.mtype=mtype
        

    def getdf(self):
        """Frequency bin width in GHz"""
        df=self.f[1]-self.f[0]
        return df

    def getnf(self):
        """Number of frequency elements"""
        return np.size(self.f)

    def append(self,v,mtype):
        """Append data array"""
        if np.size(self.v)==0:
            self.v=v
        else:
            self.v=np.dstack((self.v,v))
            
        self.mtype.append(mtype)
        

class d:

    def __init__(self,p=None,d=None):
        """Full TOD"""
        
        if p is None:
            self.point=point()
        else:
            self.point=p

        if d is None:
            self.data=data()
        else:
            self.data=d

    def save(self,fname):
        v=self.data.v
        f=self.data.f
        el=self.point.el
        t=self.point.t
        mtype=np.array(self.data.mtype)
        
        np.savez(fname,v=v,f=f,el=el,t=t,mtype=mtype)
    

