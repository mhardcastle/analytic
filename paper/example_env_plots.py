from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

axes=[]
naxes=4
plt.figure(figsize=(12,9))
for i in range(naxes):
    axes.append(plt.subplot(naxes/2,naxes/2,i+1))
for ax in axes:
    ax.set_xscale('log')

axes[0].set_yscale('log')
axes[2].set_yscale('log')

names=['beta','universal']
colours=['blue','orange']
labels=['$\\beta$ model','Universal']
for n,c,l in zip(names,colours,labels):
    outname='example-'+n+'.pickle'
    env=Evolve_RG.load(outname)
    tv=env.tv

    axes[0].plot(tv/Myr,env.R/kpc,color=c,label=l+' $R$')
    rlobep=np.sqrt(env.vl/(4*np.pi*env.R))
    axes[0].plot(tv/Myr,env.Rp/kpc,ls=':',color=c,label=l+' $R_\\perp$')
    axes[0].plot(tv/Myr,rlobep/kpc,ls='--',color=c,label=l+' $R_{\\perp, lobe}$')
    if n=='universal': axes[0].plot(tv/Myr,40*(tv/Myr)**0.6,ls='-.',color='green',label='$R \propto t^{3/5}$')
    axes[1].plot(tv/Myr,env.rlobe,color=c)
    axes[2].plot(tv/Myr,env.m1,color=c)
    axes[2].plot(tv/Myr,env.mp1,color=c,ls=':')
    axial=2*rlobep/env.R
    axes[3].plot(tv/Myr,axial,color=c,label=l)

#axis1.set_ylim(R/kpc/5,5*R/kpc)
#axis2.set_ylim(enow/5,enow*20)
axlabs=['Source size (kpc)','$V_L/V$','Mach number of expansion, ${\cal M}$','$r_{axial}$']
for ax,l in zip(axes,axlabs):
    ax.set_xlim((tv[0]/Myr,tv[-1]/Myr))
    ax.set_xlabel('$t$ (Myr)')
    ax.set_ylabel(l)

for i in [0,3]:
    axes[i].legend(loc=4,fontsize='small')
#axes[2].plot([tv[0]/Myr,tv[-1]/Myr],[1,1],ls='--')
plt.tight_layout()
plt.savefig('example_env.pdf')
