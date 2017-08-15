import numpy as np
import matplotlib.pyplot as plt
from constants import *
import glob
import matplotlib.cm as cm

from solver import Evolve_RG

if __name__=='__main__':
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

    g=glob.glob('save*.pickle')

    for f in g:
        sdict={}
        env=Evolve_RG.load(f)
        sdict['color']=cm.gist_rainbow((np.log10(env.Q)-36)/4.0)
        plt.figure(1)
        plt.plot(tv/Myr,env.R/kpc,**sdict)
        plt.figure(2)
        plt.plot(env.R/kpc,np.sqrt(env.vl/(2.0*np.pi*env.R**3.0/3.0)),**sdict)
        plt.figure(3)
        plt.plot(env.R/kpc,env.synch*env.corrs[:,0],**sdict)
        plt.figure(4)
        plt.plot(env.R/kpc,0.55+np.log(env.corrs[:,0]/env.corrs[:,1])/np.log(1400.0/150.0),**sdict)
        #axes[3].plot(env.R/kpc,env.synch)

    #axis1.set_ylim(R/kpc/5,5*R/kpc)
    #axis2.set_ylim(enow/5,enow*20)
    xlabels=['t/Myr','Lobe length/kpc','Lobe length/kpc','Lobe length/kpc']
    ylabels=['Lobe length/kpc','Axial ratio','150-MHz radio luminosity, W/Hz','Spectral index 150-1400 MHz']
    
    plt.figure(1)
    plt.xlim((tmin,tmax))
    for fig in range(4):
        plt.figure(fig+1)
        plt.xlabel(xlabels[fig])
        plt.ylabel(ylabels[fig])
       
#        ax.legend(loc=3)
#    for ax in axes[1:4]:
#        ax.set_xlabel('R/kpc')

    print 'Remember, this is lobe luminosity and lobe length, NOT source properties'
    figs=['talk-lt.png','talk-lax.png','talk-lllr.png','talk-spix.png']
    for fig in range(4):
        plt.figure(fig+1)
        plt.savefig(figs[fig])
    plt.show()
