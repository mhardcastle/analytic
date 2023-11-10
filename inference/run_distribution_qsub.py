#!/usr/bin/python

from solver import Evolve_RG
from astropy.table import Table
import sys
from constants import *
import os

t=Table.read('source-table.txt',format='ascii')
lookback=1200

n=int(sys.argv[1])
r=t[n]
print r
env=Evolve_RG('universal',M500=r['M500'],xi=0.40,q=2.1)
tlimit=lookback-r['Tstart']
tv=np.logspace(-6,np.log10(tlimit),100)*Myr
outfile='run-%i.pickle' % n
if not os.path.isfile(outfile):
    env.solve(r['Q'],tv,tstop=r['lifetime']*Myr)
    env.findb()
    env.findsynch(150e6)
    env.findcorrection((150e6,150e6*(1+r['z']),1400e6*(1+r['z'])),do_adiabatic=True,z=r['z'],timerange=(99,))
    env.save('run-%i.pickle' % n)
