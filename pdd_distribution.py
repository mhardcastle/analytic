import glob
from solver import Evolve_RG
import matplotlib.pyplot as plt
from constants import *

g=glob.glob('run-*.pickle')

plt.xscale('log')
plt.yscale('log')

count=0
for f in g:
    env=Evolve_RG.load(f)
    synch=env.synch[-1]*env.corrs[:,0][-1]
    if synch>0:
        plt.scatter(env.R[-1]/kpc,synch)
        count+=1

print 'Count is',count
plt.show()
