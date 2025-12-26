from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from synch_constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

env=Evolve_RG.load('example-universal.pickle')

env.z=6

env.finds_loss()
env.findic_loss()
env.findbs_loss()

env.findlosscorrection()

plt.plot(env.tv/Myr,env.loss*env.losscorrs,label='Synchrotron')
plt.plot(env.tv/Myr,env.ic_loss*env.losscorrs,label='Inverse-Compton')
plt.plot(env.tv/Myr,env.bs_loss,label='Bremsstrahlung')

plt.plot(env.tv/Myr,len(env.tv)*[env.Q],ls='--',label='Input power')

plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlim((1e-4,300))
plt.ylim((1e31,3e39))
plt.ylabel('Power (W)')
plt.xlabel('Time (Myr)')
plt.savefig('plot_losses_z6.pdf')
#plt.show()
