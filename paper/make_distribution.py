# Set up a distribution of sources for simulation

from astropy.table import Table
import numpy as np

count=10000
lookback=1000 #Myr
max_lifetime=500 #Myr

t=Table()

t['z']=np.random.uniform(0,4,count)
t['Q']=10**np.random.uniform(36,40,count)
t['M500']=10**np.random.uniform(13,15,count)
t['Tstart']=np.random.uniform(0,lookback,count)
t['lifetime']=10**np.random.uniform(-3,np.log10(max_lifetime),count)

t.write('source-table.txt',format='ascii')
