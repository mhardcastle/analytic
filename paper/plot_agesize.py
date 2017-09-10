import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')
t1['Dproj']=t1['D']*np.sin(t1['theta'])/kpc
t1['age']=1000-t1['Tstart']
t1r=t1[t1['remnant']]
t1nr=t1[~t1['remnant']]

#plt.scatter(t1nr['Q'],t1nr['l150'],c=np.log10(t1nr['M500']),alpha=0.6)
#plt.scatter(t1r['Q'],t1r['l150'],c=np.log10(t1r['M500']),alpha=0.6,marker='x')

plt.scatter(t1nr['age'],t1nr['Dproj'],c=np.log10(t1nr['Q']),alpha=0.6,label='Active')
plt.scatter(t1r['age'],t1r['Dproj'],c=np.log10(t1r['Q']),alpha=0.6,marker='x',label='Remnant')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Age (Myr)')
plt.ylabel('Projected size (kpc)')

cb=plt.colorbar()
cb.set_label('$\log_{10}(Q/{\\rm W})$')
plt.legend(loc=0)
plt.tight_layout()
plt.savefig('agesize.pdf')
#plt.show()
