# Plot a couple of pressure profiles

from constants import *
from solver import Evolve_RG
import numpy as np
import os.path

xi=0.4
labels=['beta','universal']
envs=[Evolve_RG('beta',kT=2.27e3*eV,p0=4e-12,rc=30*kpc,beta=0.67,xi=xi),Evolve_RG('universal',M500=1e14,xi=xi)]

tmin=0
tmax=300
tv=np.logspace(-5,np.log10(tmax),100)*Myr
Q=2e39

for env,l in zip(envs,labels):
    inname='example-'+l+'.pickle'
    outname='example-'+l+'-noad.pickle'
    if os.path.isfile(inname):
        env=Evolve_RG.load(inname)
    else:
        env.solve(Q,tv)

    env.findb()
    env.findsynch(2.1,150e6)
    env.findcorrection((150e6,330e6,1.4e9,5e9,8e9,15e9),do_adiabatic=False)
    env.save(outname)
