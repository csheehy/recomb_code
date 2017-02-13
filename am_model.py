# AM atmospheric model class
import numpy as np
from astropy import units as u
from subprocess import Popen, PIPE, STDOUT
from StringIO import StringIO
from IPython.core.debugger import Tracer; debug_here=Tracer()

def readamcfile(fname):
    """Read .amc file into a profile"""
    
    # Skip lines starting with these strings
    skiptypes=['#','!','?','\n','f','output','tol','za','PTmode']
    
    # Get new profile
    p=profile()

    f=open(fname,'r')
    l='a'
    kk=0

    while l != '':
        kk+=1
        l=f.readline()

        for k in skiptypes:
            if l.startswith(k):
                # Skip line
                continue

        if l.startswith('Nscale'):
            x=l.split()
            x[1]=x[1]+' '
            p.addnscale(''.join(x[1:-1]),1.0)
        elif l.startswith('T0'):
            x=l.split()
            p.T0=float(x[1])
        elif l.startswith('layer'):
            # New layer
            ly=layer()

            x=l.split()
            if len(x)==2:
                ly.ltype=x[1]
            
            while bool(l.strip()):
                l=f.readline()
                if l.startswith('Pbase'):
                    ly.Pbase=float(l.split()[1])
                elif l.startswith('Tbase'):
                    ly.Tbase=float(l.split()[1])
                elif l.startswith('lineshape'):
                    ly.lineshape=l.split()[1]
                elif l.startswith('column'):
                    x=l.split()
                    if len(x) == 3:
                        ly.addcolumn(x[1],x[2],'')
                    else:
                        ly.addcolumn(x[1],x[2],float(x[3]))

            # Add layer to profile
            p.addlayer(ly)

    # Close file
    f.close()

    # Return profile
    return p


class layer:

    def __init__(self,Pbase=[],Tbase=[],ltype='',lineshape=''):
        """Layer class. Pbase in mbar, Tbase in K"""
        self.Pbase = Pbase
        self.Tbase = Tbase
        self.lineshape = lineshape
        self.ltype = ltype
        self.col={}

    def addcolumn(self,species,type,val):
        """Add column to layer.
        e.g. layer.addcolumn('ch4','hydrostatic',1e-4)"""
        self.col[species]=[type,val]

    def writecolumn(self,species):
        """Return column string for named species"""
        str='column {0} {1} {2}\n'
        return str.format(species,self.col[species][0],self.col[species][1])

    def writelayer(self):
        """Return layer string for .amc file"""
        str='layer %s\n'          % self.ltype
        str+='Pbase %0.6f mbar\n' % self.Pbase
        str+='Tbase %0.3f K\n'    % self.Tbase
        if bool(self.lineshape):
            str+='lineshape %s\n' % self.lineshape
        for  k in self.col.keys():
            str+=self.writecolumn(k)
        return str

    def splitlayer(self):
        """Split layer into its constituents. For now...
        h2o -> h2o_lines, h2o_continuum
        dry_air -> ch4, co, co2, n2o, o2_coupled, o2_uncoupled, n2air, o2air
        """

        for k in self.col.keys():
            if k=='dry_air':
                # Dry air
                c=self.col.pop(k)
                addc=['ch4', 'co', 'co2', 'n2o', 'o2_coupled', 'o2_uncoupled',
                      'n2air', 'o2air']
                for j in addc:
                    self.addcolumn(j,c[0],c[1])

            if k=='h2o':
                # Dry air
                c=self.col.pop(k)
                addc=['h2o_lines','h2o_continuum']
                for j in addc:
                    self.addcolumn(j,c[0],c[1])


class profile:

    def __init__(self,layers=None,Nscale=None,T0=None):
        """
        Conglomoration of layers to form model atmospheric
        profile. Establishes empty profile if not called with list of
        layers. Frequency resolution, df, given in MHz (default=20). Nscale is
        dictionary of species and Nscale parameters,
        i.e. Nscale=Nscale={'ch4':1,'h20':1e-2}. T0 is background temperature in
        K (default 0). fmin21,fmax in GHz (default 0 and 300 GHz). za in degrees
        (default 0). 
        """
        if layers is None:
            self.layers=[]

        if Nscale is None:
            self.Nscale={}
        
        if T0 is None:
            self.T0=0.0

        # In case only one layer was handed in, make sure it's a list
        if not isinstance(self.layers,list):
            self.layers=[self.layers]

    def clear(self):
        self.layers=[]
        self.Nscale=[]
        self.T0=0.0

    def getuniquecols(self):
        """Return unique column handles in the profile"""
        c=[]
        for k in self.layers:
            c+=k.col.keys()
        # Return unique values
        return list(set(c))

    def addnscale(self,species,nscale):
        """Add (or replace if existing) species Nscale parameter in
        profile.nscale dictionary"""
        self.Nscale[species]=nscale

    def addallnscales(self):
        """Add an Nscale=1 parameter for each column present in the profile."""
        c=self.getuniquecols()
        for k in c:
            self.addnscale(k,1)

    def setnscales(self,val):
        """Set all Nscales to specified value"""
        for k in self.Nscale.keys():
            self.Nscale[k]=val
            
    def splitlayers(self):
        """Split all layers into individual constituents"""
        for k in self.layers:
            k.splitlayer()

    def addlayer(self,layer):
        """Add layer to profile. 
        e.g. profile.addlayer(layer)"""
        self.layers.append(layer)
    
    def writelayers(self):
        """Write all layers to a string"""
        str='';
        for k,val in enumerate(self.layers):
            str+=self.layers[k].writelayer()
            str+='\n'
        return str

    def writeheader(self):
        str='?\n'
        str+='? usage: am file.amc  fmin[GHz]  fmax[GHz]  zenith_angle[deg]\n'
        str+='?\n'
        return str

    def writeio(self):
        str='f {fmin} GHz {fmax} GHz {df} MHz\n'
        str+='output f GHz  tau  Tb K I watt*cm-2*GHz-1*sr-1\n'
        str+='za {za} deg\n'
        str+='tol {tol}\n'
        str+='PTmode Pbase Tbase\n'
        return str
    
    def writenscale(self):
        """Write Nscale string."""
        str=''
        for k,val in enumerate(self.Nscale):
            str+=('Nscale {0} {1}\n').format(val,self.Nscale[val])
        return str

    def writeT0(self):
        """Write background temperature line"""
        str=('T0 {0} K\n').format(self.T0)
        return str
    
    def writeprof(self):
        """Write entire .amc profile file"""
        str=''

        str+=self.writeheader()
        str+='\n'

        str+=self.writeio()
        str+='\n'

        str+=self.writenscale()
        str+='\n'

        str+=self.writeT0()
        str+='\n'
        
        str+=self.writelayers()
        str+='\n'

        return str
    


class am:

    def __init__(self,prof=[],df=20,tol=.0001,fmin=0,fmax=300,za=0):
        """Class for calling am from python. Frequency resolution, df, given in
        MHz (default=20). fmin,fmax in GHz (default 0 and 300 GHz). za in degrees
        (default 0)."""
        self.prof=prof
        self.tol=tol
        self.df=df
        self.fmin=fmin
        self.fmax=fmax
        self.za=za

    def writeamc(self):
        return self.prof.writeprof().format(za=self.za,fmin=self.fmin,
                                            fmax=self.fmax,tol=self.tol,df=self.df) 

    def callam(self):
        """Call and store standard output as am.results"""
        p = Popen(['am', '-'], stdout=PIPE, stdin=PIPE, stderr=PIPE) 
        self.results, self.info = p.communicate(input=self.writeamc())
        self._parseresults()
        
    def _parseresults(self):
        """Parse contents of am.results and store f, tau, Tb (in K)
        and I (W/m^2/Hz/sr)."""
        out=np.genfromtxt(StringIO(self.results));

        self.f=out[:,0] # frequency
        self.tau=out[:,1] # optical depth
        self.Tb=out[:,2] # brightness temperature
        self.I=(out[:,3]*1e-5) # specific intensity

    def writetofile(self,fname):
        f=open(fname,'w')
        f.write(self.results)
        f.close


