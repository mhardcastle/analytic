from solver import Evolve_RG
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

q=2.1

env=Evolve_RG.load('example-tstop-universal.pickle')
env.gmin=10
env.gmax=1e6
env.q=q

plt.figure(figsize=(12,6))

for i,z in enumerate([0,2]):

    plt.subplot(1,2,i+1)
    plt.xlim(0,300)
    env.findsynch(150e6)
    env.findic(2.4e17,z=z)
    if z>0:
        env.findcorrection([150e6],z=z,do_adiabatic=True)
    env.ic_findcorrection([2.4e17],z=z,do_adiabatic=True)
    
    plt.plot(env.tv/Myr,env.corr_synch[:,0]/1e10,label='Synchrotron (scaled)')
    plt.plot(env.tv/Myr,env.ic,label='IC uncorrected')
    plt.plot(env.tv/Myr,env.corr_ic[:,0],label='IC corrected')
    if i==0:
        plt.legend(loc=0)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Luminosity (W Hz$^{-1}$)')
    plt.yscale('log')
    plt.ylim((1e16,3.5e20))

plt.tight_layout()
plt.savefig('remnant_ic.pdf')
plt.show()
