import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table
rc('font',**{'family':'serif','serif':['Times'],'size':12})
rc('text', usetex=True)

t1=Table.read('/home/mjh/distribution/obs-table.fits')
t2=t1[t1['l150']>3e25]
plt.figure(figsize=(12,6))

for i,t in enumerate([t1,t2]):
    plt.subplot(1,2,1+i)
    _,bins,_=plt.hist(t['z'],alpha=0.5,bins=40,color='red',label='Full sample')
    plt.hist(t[t['remnant']]['z'],alpha=0.5,bins=bins,color='gray',label='Remnants')
    plt.xlabel('$z$')
    if i==0:
        plt.ylabel('Number of objects')
        plt.legend(loc=0)

#plt.show()
plt.savefig('remnant-fraction.pdf')
