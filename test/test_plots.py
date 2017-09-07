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

env=Evolve_RG('beta',kT=0.86e3*eV,p0=2e-12,rc=22*1.528714*kpc,beta=0.67,do_adiabatic=False,q=2.1)
#env=evolve_rg('universal',M500=1e14)
for Q in np.logspace(37.5,39.5,5):
    print 'solving for',Q
    env.solve(Q,tv)
    env.findb()
    env.findsynch(150e6)
    env.findcorrection((150e6,))

    axes[0].plot(tv/Myr,env.R/kpc)
    axes[1].plot(tv/Myr,env.rlobe,label='Q = %.2g W' % Q)
    axes[2].plot(tv/Myr,env.m1)
    axes[3].plot(env.R/kpc,env.synch*env.corrs[:,0],label='Q = %.2g W' % Q)

    env.save('save_%.1f.pickle' % np.log10(Q))
#axis1.set_ylim(R/kpc/5,5*R/kpc)
#axis2.set_ylim(enow/5,enow*20)
for ax in axes[0:3]:
    ax.set_xlim((tmin,tmax))
    ax.set_xlabel('t/Myr')
    ax.legend(loc=3)
axes[2].plot([tmin,tmax],[1,1],ls='--')
axes[3].set_xlabel('R/kpc')
plt.show()
