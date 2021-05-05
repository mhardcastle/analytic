# Example dwarf galaxy profile

from constants import *
from solver import Evolve_RG
import numpy as np
import os.path

xi=0.4
q=2.1
env=Evolve_RG('beta',kT=1e4*boltzmann,p0=10**0.75*3.1e-16,rc=0.8*kpc,beta=0.67,xi=xi,q=q)
outname='dwarf.pickle'

tmin=0
tmax=100
tv=np.logspace(-5,np.log10(tmax),100)*Myr
Q=1e35

env.solve(Q,tv)
env.save(outname)

env.findb()
env.findsynch(150e6)
env.findcorrection((150e6,330e6,1.4e9,5e9,8e9,15e9),do_adiabatic=True)
env.save(outname)
