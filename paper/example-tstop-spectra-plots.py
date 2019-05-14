from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

env=Evolve_RG.load('example-tstop-universal.pickle')

plt.figure(figsize=(12,5.5))

plt.subplot(1,2,1)
nfreq=len(env.freqs)
tv=env.tv
lums=np.zeros_like(env.corrs)
plt.plot(tv/Myr,env.synch,ls='--',color='black',label='150 MHz uncorrected')
for i in range(nfreq):
    lums[:,i]=env.synch*(env.freqs[i]/env.nu_ref)**-env.alpha*env.corrs[:,i]
    plt.plot(tv/Myr,lums[:,i],label='%.0f MHz' % (env.freqs[i]/1e6))
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0,fontsize='small')
plt.xlabel('Time (Myr)')
plt.ylabel('Radio luminosity (W Hz$^{-1}$)')
plt.xlim((1,300))
plt.ylim((1e25,3e29))
plt.subplot(1,2,2)
for i in range(1,nfreq):
    plt.plot(tv/Myr,-np.log(lums[:,i]/lums[:,i-1])/np.log(env.freqs[i]/env.freqs[i-1]),label='$\\alpha_{%.0f}^{%.0f}$' % (env.freqs[i]/1e6,env.freqs[i-1]/1e6))
plt.xscale('log')
plt.legend(loc=0,fontsize='small')
plt.xlabel('Time (Myr)')
plt.xlim((1,300))
plt.ylim((0.5,2.5))
plt.ylabel('Spectral index')
plt.savefig('example-tstop-spectra.pdf')
#plt.show()
