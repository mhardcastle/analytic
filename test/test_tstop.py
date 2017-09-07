from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *

axes=[]
naxes=4
for i in range(naxes):
    axes.append(plt.subplot(naxes/2,naxes/2,i+1))
for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')

tmin=0
tmax=300
tv=np.logspace(-6,np.log10(tmax),100)*Myr

#env=evolve_rg('universal',M500=1e14)
Q=10**38.5
for i in range(2):
    env=Evolve_RG('beta',kT=0.86e3*eV,p0=2e-12,rc=22*1.528714*kpc,beta=0.67,do_adiabatic=True,q=2.1)
    env.solve(Q,tv,tstop=Myr*(tmax-(i*200)))
    env.findb()
    env.findsynch(150e6)
    env.findcorrection((150e6,))

    axes[0].plot(tv/Myr,env.R/kpc)
    axes[1].plot(tv/Myr,env.rlobe,label='i = %i' % i)
    axes[2].plot(tv/Myr,env.m1)
    axes[3].plot(env.R/kpc,env.synch*env.corrs[:,0])

    env.save('save_tstop_%i.pickle' % i)
#axis1.set_ylim(R/kpc/5,5*R/kpc)
#axis2.set_ylim(enow/5,enow*20)
for ax in axes[0:3]:
    ax.set_xlim((tmin,tmax))
    ax.set_xlabel('t/Myr')
    ax.legend(loc=3)
axes[2].plot([tmin,tmax],[1,1],ls='--')
axes[3].set_xlabel('R/kpc')
plt.show()
