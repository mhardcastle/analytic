from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
import os.path

plt.figure(figsize=(12,5.5))

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

i=0
for beta in [0.55,0.75,0.90]:
    for rc in [20,30,40]:
        color=['orange','red','firebrick','yellowgreen','lime','darkgreen','dodgerblue','royalblue','darkblue'][i]
        i+=1

        env=Evolve_RG.load('beta-%.2f-%i.pickle' % (beta,rc))
        tv=env.tv
        plt.subplot(1,2,1)
        plt.plot(tv/Myr,env.R/kpc,color=color,ls='--')

        t=np.loadtxt('/home/mjh/PLUTO/results2/bmb-%i-%i-2s' % (int(beta*100),rc),skiprows=1)
        stu=2.9
        slu=2.1
        plt.plot(t[:,0]*stu,(t[:,24]+t[:,11])*0.5*slu,color=color,label='$\\beta = %.2f$ $r_c=%.0f$ kpc' % (env.beta,env.rc/kpc))

        plt.subplot(1,2,2)
        plt.plot(tv/Myr,env.vt/(kpc**3.0),color=color,ls='--')
        plt.plot(t[:,0]*stu,(t[:,15]+t[:,2]+t[:,14]+t[:,1])*0.5*slu**3.0,color=color,label='$\\beta = %.2f$ $r_c=%.0f$ kpc' % (env.beta,env.rc/kpc))


for i in range(2):
    plt.subplot(1,2,i+1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,100])
    plt.xlabel('Time (Myr)')

plt.subplot(1,2,1)
plt.ylim([10,300])
plt.ylabel('Lobe length (kpc)')
plt.legend(loc=0,fontsize='small')
plt.subplot(1,2,2)
plt.ylim([1e4,1e7])
plt.ylabel('Shock volume (kpc$^3$)')
plt.savefig('hk13-compare.pdf')

