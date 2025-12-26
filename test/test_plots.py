from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from synch_constants import *

axes=[]
naxes=4
ylabs=['Length/kpc','Axial ratio','Mach number','L_150/(W/Hz)']
for i in range(naxes):
    axes.append(plt.subplot(naxes//2,naxes//2,i+1))
for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')

tmin=10**-6
tmax=300
tv=np.logspace(np.log10(tmin),np.log10(tmax),100)*Myr

env=Evolve_RG('beta',kT=0.86e3*eV,p0=2e-12,rc=22*1.528714*kpc,beta=0.67,do_adiabatic=False,q=2.1)
#env=evolve_rg('universal',M500=1e14)
for Q in np.logspace(37.5,39.5,5):
    print('solving for jet power',Q,'W')
    env.solve(Q,tv)
    env.findb()
    env.findsynch(150e6)
    env.findcorrection((150e6,))

    axes[0].plot(tv/Myr,env.R/kpc,label='Q = %.2g W' % Q)
    axes[1].plot(tv/Myr,env.rlobe,label='Q = %.2g W' % Q)
    axes[2].plot(tv/Myr,env.m1,label='Q = %.2g W' % Q)
    axes[3].plot(env.R/kpc,env.synch*env.corrs[:,0],label='Q = %.2g W' % Q)

    env.save('save_%.1f.pickle' % np.log10(Q))
#axis1.set_ylim(R/kpc/5,5*R/kpc)
#axis2.set_ylim(enow/5,enow*20)
for i,ax in enumerate(axes):
    #ax.set_xlim((tmin,tmax))
    if i!=3: ax.set_xlabel('t/Myr')
    else: ax.set_xlabel('R/kpc')
    ax.set_ylabel(ylabs[i])
    ax.legend(loc=0)
axes[2].plot([tmin,tmax],[1,1],ls='--')
plt.show()
