# Update the distribution table with the values from simulations

from astropy.table import Table
import numpy as np
from solver import Evolve_RG
import os
from tqdm import tqdm

t=Table.read('source-table.txt',format='ascii')

l150=[]
l150_rest=[] # to use in calc of flux -- luminosity at rest-frame freq of observation
d=[]
alpha=[]
live=[]
remnant=[]
u=np.random.uniform(size=len(t))
theta=np.arccos(1-u)

for i,r in enumerate(tqdm(t)):
    inname='run-%i.pickle' % i
    if os.path.isfile(inname):
        env=Evolve_RG.load(inname)
        synch=env.synch[-1]*env.corrs[-1,0]
        live.append(synch>0)
        l150.append(synch)
        l150_rest.append(env.synch[-1]*env.corrs[-1,1]*(1+env.z)**-0.55)
        d.append(env.R[-1])
        if synch>0:
            alpha.append(0.55+np.log(env.corrs[-1,1]/env.corrs[-1,2])/np.log(1400.0/150.0))
        else:
            alpha.append(np.nan)
        remnant.append((r['lifetime']+r['Tstart'])<1000.0)
    else:
        live.append(False)
        l150.append(np.nan)
        l150_rest.append(np.nan)
        d.append(np.nan)
        alpha.append(np.nan)
        remnant.append(False)
t['l150']=l150
t['l150_rest']=l150_rest
t['D']=d
t['alpha']=alpha
t['live']=live
t['remnant']=remnant
t['theta']=theta
t.write('source-table.fits',overwrite=True)
