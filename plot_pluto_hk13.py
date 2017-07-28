from solver import Evolve_RG
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
import os.path

rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

beta=0.90
rc=40

env=Evolve_RG.load('beta-%.2f-%i.pickle' % (beta,rc))
tv=env.tv

plt.plot(tv/Myr,env.R/kpc,label='$\\beta = %.2f$ $r_c=%.1f$ kpc' % (env.beta,env.rc/kpc))
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.5,1e3])
plt.xlim([0.1,320])
plt.xlabel('Time (Myr)')
plt.ylabel('Lobe length (kpc)')

t=np.loadtxt('/home/mjh/PLUTO/results2/bmb-%i-%i-2s' % (int(beta*100),rc),skiprows=1)
stu=2.9
slu=2.1
plt.plot(t[:,0]*stu,t[:,11]*slu)
plt.plot(t[:,0]*stu,t[:,24]*slu)


plt.show()
