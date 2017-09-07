# Plot a couple of pressure profiles

from constants import *
from solver import Evolve_RG
import numpy as np
import os.path

xi=0.4
q=2.1
labels=['beta','universal']
envs=[Evolve_RG('beta',kT=2.27e3*eV,p0=4e-12,rc=30*kpc,beta=0.67,xi=xi,q=q),Evolve_RG('universal',M500=1e14,xi=xi,q=q)]

tmin=0
tmax=300
tv=np.logspace(-5,np.log10(tmax),300)*Myr
Q=2e39

for env,l in zip(envs,labels):
    outname='example-tstop-'+l+'.pickle'
    if os.path.isfile(outname):
        env=Evolve_RG.load(outname)
    else:
        env.solve(Q,tv,tstop=100*Myr)
#    try:
#        dummy=env.corrs
#    except:
    env.findb()
    env.findsynch(150e6)
    env.findcorrection((150e6,330e6,1.4e9,5e9,8e9,15e9),do_adiabatic=True)
    env.save(outname)
