import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')

t1r=t1[t1['remnant']]
t1nr=t1[~t1['remnant']]


#plt.scatter(t1nr['Q'],t1nr['l150'],c=np.log10(t1nr['M500']),alpha=0.6)
#plt.scatter(t1r['Q'],t1r['l150'],c=np.log10(t1r['M500']),alpha=0.6,marker='x')
plt.scatter(t1nr['Q'],t1nr['l150'],c=t1nr['z'],alpha=0.6)
plt.scatter(t1r['Q'],t1r['l150'],c=t1r['z'],alpha=0.6,marker='x')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$Q$')
plt.ylabel('$L_{150}$ (W Hz$^{-1}$)')

# Willott
f=15
q=np.logspace(36,40,100)
l150=4*np.pi*1e28*(q/(3e38*f**1.5))**(7.0/6.0)
plt.plot(q,l150,color='red')

l150=5e27*(q/1e38)
plt.plot(q,l150,color='blue')
l150=5e25*(q/1e38)
plt.plot(q,l150,color='blue')


plt.colorbar()
plt.show()
