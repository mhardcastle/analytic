#!/usr/bin/python
from __future__ import division
import numpy as np
from astropy.table import Table
from constants import *

ft1=Table.read('/home/mjh/distribution/source-table.fits')
ft2=Table.read('/home/mjh/distribution2/source-table.fits')

t1=Table.read('/home/mjh/distribution/obs-table.fits')
t2=Table.read('/home/mjh/distribution2/obs-table.fits')

for t in [t1,t2]:
    t['psize']=np.sin(t['theta'])*t['D']

print 'Sample (i) live',np.sum(ft1['live']),'observed',len(t1)
print 'Sample (ii) live',np.sum(ft2['live']),'observed',len(t2)

print 'Sample (i) median projected size',np.median(t1['psize'])/kpc
print 'Sample (ii) median projected size',np.median(t2['psize'])/kpc

print 'Sample (i) median lifetime',np.median(ft1['lifetime'])
print 'Sample (ii) median lifetime',np.median(ft2['lifetime'])

print 'Sample (i) remnant fraction',np.sum(t1['remnant'])/len(t1)
print 'Sample (ii) remnant fraction',np.sum(t2['remnant'])/len(t2)

s1z=t1[t1['z']<1.0]
s1z1=s1z[s1z['l150']>3e25]
print 'Sample (i) remnant fraction with z<1',np.sum(s1z['remnant'])/len(s1z)
print 'Sample (i) remnant fraction with z<1, high L',np.sum(s1z1['remnant'])/len(s1z1)

s1z=t1[(t1['z']>=1.0) & (t1['z']<2.0)]
s1z1=s1z[s1z['l150']>3e25]
print 'Sample (i) remnant fraction with 1<z<2',np.sum(s1z['remnant'])/len(s1z)
print 'Sample (i) remnant fraction with 1<z<2, high L',np.sum(s1z1['remnant'])/len(s1z1)

remnant=t1[t1['remnant']]
print 'Remnants ultra-steep',np.sum(remnant['alpha']>1.5)/len(remnant)
