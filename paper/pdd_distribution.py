import matplotlib.pyplot as plt
from astropy.table import Table
from constants import *

t=Table.read('obs-table.fits')
t=t[t['live']]
plt.xscale('log')
plt.yscale('log')

count=0
colors=['red' if r['remnant'] else 'blue' for r in t]
plt.scatter(t['D']/kpc,t['l150'],c=colors,alpha=0.5)

plt.show()
