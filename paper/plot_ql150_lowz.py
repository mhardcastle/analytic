import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table,vstack
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')

t1=t1[t1['z']<0.5]

t1r=t1[t1['remnant']]
t1nr=t1[~t1['remnant']]

plt.figure(figsize=(6,7))

plt.scatter(t1nr['Q'],t1nr['l150'],c=np.log10(t1nr['M500']),alpha=0.6,label='Active')
plt.scatter(t1r['Q'],t1r['l150'],c=np.log10(t1r['M500']),alpha=0.6,marker='x',label='Remnant')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$Q$ (W)')
plt.ylabel('$L_{150}$ (W Hz$^{-1}$)')

# Willott
#f=15
q=np.logspace(36,40,100)
kwargs={'label':'Willott+ ($f=5$)'}
for f in [5]:
    l150=4*np.pi*1e28*(q/(3e38*f**1.5))**(7.0/6.0)
    plt.plot(q,l150,color='red',**kwargs)
    kwargs={}

# Ineson
l150=4*np.pi*1e28*(q/5e39)**(1.0/0.89)
plt.plot(q,l150,color='orange',label='Ineson+ 17')

kwargs={'label':'Best linear fit'}
for norm in [3e27]:
    l150=norm*(q/1e38)
    plt.plot(q,l150,color='blue',**kwargs)
    kwargs={}

plt.ylim((3e23,3e30))

cb=plt.colorbar(orientation='horizontal')
cb.set_label('$M_{500}$')
plt.legend(loc=0)
plt.tight_layout()
#plt.show()
plt.savefig('plot_ql150_lowz.pdf')
