from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

plt.figure(figsize=(12,5.5))
done_uncorr=False
for name,label in zip(['example-universal-noad','example-universal','example-universal-highz'],['No adiabatic, $z=0$','Adiabatic, $z=0$','Adiabatic, $z=2$']):

    env=Evolve_RG.load(name+'.pickle')

    plt.subplot(1,2,1)
    nfreq=len(env.freqs)
    tv=env.tv
    lums=np.zeros_like(env.corrs)
    if not done_uncorr:
        plt.plot(tv/Myr,env.synch,ls='--',color='black',label='150 MHz uncorrected')
    for i in range(nfreq):
        lums[:,i]=env.synch*(env.freqs[i]/env.nu_ref)**-env.alpha*env.corrs[:,i]
    ploti=0
    plt.plot(tv/Myr,lums[:,ploti],label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0,fontsize='small')
    plt.xlabel('Time (Myr)')
    plt.ylabel('150-MHz radio luminosity (W Hz$^{-1}$)')
    plt.subplot(1,2,2)
    if not done_uncorr:
        plt.plot(tv/Myr,[0.55]*len(tv),ls='--',color='black',label='150 MHz uncorrected')
        done_uncorr=True
    ploti=1
    plt.plot(tv/Myr,-np.log(lums[:,ploti]/lums[:,ploti-1])/np.log(env.freqs[ploti]/env.freqs[ploti-1]),label=label)
    plt.xscale('log')
    plt.legend(loc=0,fontsize='small')
    plt.xlabel('Time (Myr)')
    plt.ylabel('Spectral index 150-330 MHz')

plt.savefig('example-spectra-vary.pdf')
