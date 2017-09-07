from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
import os.path

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)


tmin=0
tmax=200
tv=np.logspace(-6,np.log10(tmax),200)*Myr

#env=evolve_rg('universal',M500=1e14)
Q=2e38
envs=[]
for beta in [0.55,0.75,0.90]:
    for rc in [20,30,40]:
        outname='beta-%.2f-%i.pickle' % (beta, rc)
        if os.path.isfile(outname):
            envs.append(Evolve_RG.load(outname))
        else:
            rc*=2.1*kpc
            envs.append(Evolve_RG('beta',kT=2e3*eV,p0=5.76e-12,rc=rc,beta=beta,qfactor=0.5*25*730e3,Gamma=(5.0/3.0)))
            envs[-1].solve(Q,tv)
            envs[-1].save(outname)

for i in range(2):
    for env in envs:
        plt.plot(tv/Myr,env.R/kpc,label='$\\beta = %.2f$ $r_c=%.1f$ kpc' % (env.beta,env.rc/kpc))
    if i==0: plt.plot(tv/Myr,tv*3e8/kpc,ls='--',color='cyan')
    plt.plot(tv/Myr,5*(tv/Myr)**0.6,ls='--',color='blue')
    plt.plot(tv/Myr,tv*envs[-1].cs/kpc,ls='--',color='orange')
    plt.xscale('log')
    plt.yscale('log')
    if i==0:
        plt.ylim([1e-3,1e3])
    else:
        plt.xlim([2,300])
        plt.ylim([10,600])
    plt.xlabel('Time (Myr)')
    plt.ylabel('Lobe length (kpc)')
    plt.legend(loc=0,fontsize='small')

    plt.savefig('hk13beta-%i.pdf' % i)
    plt.clf()
    
