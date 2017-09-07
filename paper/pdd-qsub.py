import numpy as np
import pickle
import os

from constants import *
from solver import Evolve_RG

if __name__=='__main__':
    import sys
    number=int(sys.argv[1])-1
    qnum=number//10
    mnum=number%10
    print qnum, mnum

    tmin=0
    tmax=500
    tv=np.logspace(-6,np.log10(tmax),100)*Myr
    Qs=np.logspace(36,40,13)
    masses=np.logspace(13,15,10)
    Q=Qs[qnum]
    m=masses[mnum]
    print 'Doing the run for Q=',Q,'M=',m
    env=Evolve_RG('universal',M500=m,xi=0.40,q=2.1)
    outfile='save_%.1f_%.1f.pickle' % (np.log10(Q),np.log10(m))
    if not os.path.isfile(outfile):
        print 'solving for',Q,m
        env.solve(Q,tv)
        env.findb()
        env.findsynch(150e6)
        env.findcorrection((150e6,1400e6),do_adiabatic=True)
        env.save(outfile)
