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
tv=np.logspace(-6,np.log10(tmax),100)*Myr
Q=2e39

for env,l in zip(envs,labels):
    outname='example-'+l+'.pickle'
    if os.path.isfile(outname):
        envs.append(Evolve_RG.load(outname))
    else:
        env.solve(Q,tv)
        env.save(outname)

