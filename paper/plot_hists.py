# Plot observational histograms

import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table
rc('font',**{'family':'serif','serif':['Times'],'size':18})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')
t2=Table.read('/home/mjh/distribution2/obs-table.fits')
for t in [t1,t2]:
    t['alpha']=np.where(np.isinf(t['alpha']),10,t['alpha'])

plt.figure(figsize=(18,6))
plt.subplot(1,3,1)
_,bins,_=plt.hist(np.log10(np.sin(t1['theta'])*t1['D']/kpc),alpha=0.5,normed=1,bins=30,color='red',label='Sample i')
plt.hist(np.log10(np.sin(t2['theta'])*t2['D']/kpc),alpha=0.5,normed=1,bins=bins,color='blue',label='Sample ii')
plt.xlabel('$\log_{10}(D_{\\rm projected}/{\\rm kpc})$')
plt.ylabel('Number of objects (normalized)')
plt.legend(loc=0)

plt.subplot(1,3,2)
_,bins,_=plt.hist(np.log10(t1['l150']),alpha=0.5,normed=1,bins=30,color='red')
plt.hist(np.log10(t2['l150']),alpha=0.5,normed=1,bins=bins,color='blue')
plt.xlabel('$\log_{10}(L_{150}/{\\rm W\ Hz^{-1}})$')

plt.subplot(1,3,3)
_,bins,_=plt.hist(t1['alpha'],alpha=0.5,normed=1,bins=30,range=(0.5,2.0),color='red')
plt.hist(t2['alpha'],alpha=0.5,normed=1,bins=bins,color='blue')
plt.xlabel('$\\alpha$')

plt.tight_layout()
plt.savefig('sample-hists.pdf')
