from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)


tmin=0
tmax=300
tv=np.logspace(-6,np.log10(tmax),100)*Myr

#env=evolve_rg('universal',M500=1e14)
Q=10**38
for beta in [0.55,0.75,0.90]:
    for rc in [20,30,40]:
        rc*=2.1*kpc

        env=Evolve_RG('beta',kT=2e3*eV,p0=1e-11,rc=rc,beta=beta,do_adiabatic=True)
        env.solve(Q,tv)
        plt.plot(tv/Myr,env.R/kpc,label='$\\beta = %.2f$ $r_c=%.1f$ kpc' % (beta,rc/kpc))

plt.legend(loc=0,fontsize='small')
plt.plot(tv/Myr,tv*3e8/kpc,ls='--')
plt.plot(tv/Myr,tv*env.cs/kpc,ls='--')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-3,1e3])
plt.xlabel('Time (Myr)')
plt.ylabel('Lobe length (kpc)')
plt.show()
