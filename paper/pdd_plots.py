from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from constants import *
import glob
import matplotlib.cm as cm
from matplotlib import rc
import sys

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

from solver import Evolve_RG

if __name__=='__main__':

    plotnumber=int(sys.argv[1])
    t3c=Table.read('3crr.txt',format='ascii')
    #    t3c=t3c[t3c['Redshift']<0.5]
    
    nfigs=4
    for fig in range(1,nfigs):
        plt.figure(fig)
        plt.xscale('log')
        plt.yscale('log')
    plt.figure(4)
    plt.xscale('log')
    tmin=0
    tmax=300
    tv=np.logspace(-6,np.log10(tmax),100)*Myr

    names=['save_*.pickle','highz*.pickle','z6*.pickle']
    g=glob.glob(names[plotnumber])

    for f in g:
        sdict={}
        env=Evolve_RG.load(f)
        sdict['color']=cm.gist_rainbow((np.log10(env.Q)-36)/4.0)
        sdict['alpha']=((np.log10(env.m500)-11.5)/3.5)**2.0
        plt.figure(1)
        plt.plot(tv/Myr,env.R/kpc,**sdict)
        plt.figure(2)
        plt.plot(env.R/kpc,np.sqrt(env.vl/(np.pi*env.R**3.0)),**sdict)
        plt.figure(3)
        plt.plot(env.R/kpc,env.synch*env.corrs[:,0],zorder=2,**sdict)
        plt.figure(4)
        plt.plot(env.R/kpc,0.55+np.log(env.corrs[:,0]/env.corrs[:,1])/np.log(1400.0/150.0),**sdict)
        #axes[3].plot(env.R/kpc,env.synch)

    #axis1.set_ylim(R/kpc/5,5*R/kpc)
    #axis2.set_ylim(enow/5,enow*20)
    xlabels=['t (Myr)','Lobe length (kpc)','Lobe length (kpc)','Lobe length (kpc)']
    ylabels=['Lobe length (kpc)','Axial ratio','150-MHz radio luminosity (W Hz$^{-1}$)','Spectral index 150-1400 MHz']
    
    plt.figure(1)
    plt.xlim((tmin,tmax))
    for fig in range(4):
        plt.figure(fig+1)
        plt.xlabel(xlabels[fig])
        plt.ylabel(ylabels[fig])
       
#        ax.legend(loc=3)
#    for ax in axes[1:4]:
#        ax.set_xlabel('R/kpc')

    print 'Remember, this is lobe length, NOT source total length'

    plt.figure(3)
    plt.xlim((1e-3,2e3))
    plt.scatter(t3c['Size']/2.0,t3c['L_178']*(178.0/150.0)**t3c['alpha']*np.pi*4,alpha=0.5,color='grey',label='3CRR',zorder=3)
    plt.figure(4)
    plt.scatter(t3c['Size']/2.0,t3c['alpha'],alpha=0.5,color='grey',label='3CRR',zorder=3)
    plt.xlim((1e-3,2e3))
    plt.ylim((0.55,1.2))
    plt.figure(2)
    plt.xlim((1e-3,2e3))
    figs=['lt.pdf','lax.pdf','lllr.pdf','spix.pdf']
    prefixes=['','highz-','z6-']

    for fig in range(4):
        filename=prefixes[plotnumber]+figs[fig]
        plt.figure(fig+1)
        plt.savefig(filename)
    #plt.show()
