import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table,vstack
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')

t1=t1[t1['z']<0.5]

t1['LX']=1.74e44*(t1['M500']/2e14)**1.96

t1r=t1[t1['remnant']]
t1nr=t1[~t1['remnant']]

plt.scatter(t1nr['l150'],t1nr['LX'],c=np.log10(t1nr['Q']),alpha=0.6,label='Active')
plt.scatter(t1r['l150'],t1r['LX'],c=np.log10(t1r['Q']),alpha=0.6,marker='x',label='Remnant')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$L_{150}$ (W Hz$^{-1}$)')
plt.ylabel('2-10 keV $L_{\\rm X}$ (erg s$^{-1}$)')

plt.colorbar(label='$\\log_{10}(Q/{\\rm W})$')
plt.tight_layout()
plt.savefig('plot_150m500.pdf')
#plt.show()

