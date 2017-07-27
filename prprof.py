# Plot a couple of pressure profiles

import matplotlib.pyplot as plt
from constants import *
from solver import Evolve_RG
import numpy as np

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)


labels=['$\\beta$ model','Universal']
envs=[Evolve_RG('beta',kT=2.27e3*eV,p0=4e-12,rc=30*kpc,beta=0.67),Evolve_RG('universal',M500=1e14)]

dist=np.logspace(-3,3.3,100)
for env,l in zip(envs,labels):
    plt.plot(dist,env.pr(dist*kpc),label=l)
plt.xlabel('Distance (kpc)')
plt.ylabel('Pressure (Pa)')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)
plt.tight_layout()
plt.savefig('prprof.pdf')
