# Spectrum class
import numpy as np

class spec:

    def __init__(self,f,I):
        self.f = f
        self.I = I

    def dI(self,n):
        """Compute numerical derivatives of spectrum assuming a uniform
        frequency spacing. Currently not very useful."""
        nshift=20;
        df=self.f[2]-self.f[1];
        y=self.I;
        # Take numerical derivative with a circular shift
        for k in range(0,n):
            y=(y-np.roll(y,nshift));
        # Set the wrapped around points to NaN
        y[0:n*nshift]=np.NaN;
        return y/np.power(df*nshift,n)

    def polyfit(self,order=None,fmin=None,fmax=None):
        """
        Fit a polynomial to log(I) vs. log(f) in the range [fmin,fmax] (default 
        full range). Return best fit over full range.
        
        fmin,fmax = range in which to perform fit (default self.f.min() and
                    self.f.max()) 
        order = polynomial order (default = 4)
        """

        # Default frequency range over which to fit
        if fmin is None:
            fmin=self.f.min()
        if fmax is None:
            fmax=self.f.max()
        
        # Default poly order
        if order is None:
            order=4

        # Do fit
        fitind=(self.f<fmax) & (self.f>fmin)
        x=np.log(self.f[fitind])
        y=np.log(self.I[fitind])
        p=np.polyfit(x,y,order)

        # Return best fit over range and NaN elsewhere
        yfit=np.empty(self.I.shape);yfit.fill(np.nan)
        yfit[fitind]=np.exp(np.polyval(p,x))

        return yfit

