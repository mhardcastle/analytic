# Set up a distribution of sources for simulation

from astropy.table import Table
import numpy as np
from scipy.optimize import brentq
import sys

def rpow(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

def schechter(m,alpha):
    # m is scaled mass, i.e. M/M*
    return m**(-alpha)*np.exp(-m)

def makeschechter(mlower=1e13,mupper=1e15,mstar=3e14,alpha=1.5,size=1):

    sschechter=lambda m: schechter(m/mstar,alpha)/schechter(mlower/mstar,alpha)

    # don't allow it to go all the way to infinity
    upperlimit=1-sschechter(mupper)

    r=np.random.uniform(high=upperlimit,size=size)
    masses=[]
    for v in r:
        m=brentq(lambda m: (1-sschechter(m))-v,mlower,mupper)
        masses.append(m)

    return np.array(masses)

count=10000
lookback=1200 #Myr
max_lifetime=1000 #Myr

t=Table()

t['Q']=10**np.random.uniform(34,40,size=count)
t['M500']=makeschechter(size=count)
t['Tstart']=np.random.uniform(0,lookback,count)
t['lifetime']=np.random.uniform(0,max_lifetime,count)
t['z']=float(sys.argv[1])

t.write('source-table.txt',format='ascii',overwrite=True)
