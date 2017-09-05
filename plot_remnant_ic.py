from solver import Evolve_RG
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

q=2.1

env=Evolve_RG.load('example-tstop-universal.pickle')

plt.figure(figsize=(12,6))

for i,z in enumerate([0,2]):

    plt.subplot(1,2,i+1)
    plt.xlim(0,300)
    env.findsynch(q,150e6)
    env.findic(q,2.4e17,z=z)
    if z>0:
        env.findcorrection([150e6],z=z,do_adiabatic=True)
    env.ic_findcorrection([2.4e17],z=z,do_adiabatic=True)
    
    plt.plot(env.tv/Myr,env.synch*env.corrs[:,0]/1e10,label='Synchrotron (scaled)')
    plt.plot(env.tv/Myr,env.ic,label='IC uncorrected')
    plt.plot(env.tv/Myr,env.ic*env.ic_corrs[:,0],label='IC corrected')
    if i==0:
        plt.legend(loc=0)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Luminosity (W Hz$^{-1}$)')
    plt.yscale('log')
    plt.ylim((1e15,3.5e19))

plt.tight_layout()
plt.savefig('remnant_ic.pdf')
plt.show()
