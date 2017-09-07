import numpy as np
import matplotlib.pyplot as plt
from constants import *
from matplotlib import rc
from astropy.table import Table
rc('font',**{'family':'serif','serif':['Times'],'size':14})
rc('text', usetex=True)

def l150(q,m,m_exponent,l,l_exponent,z,z_exponent):
    return 6e26*(q/1e38)*(m/1e14)**m_exponent*(l/3e22)**l_exponent*(1+z)**z_exponent

t1=Table.read('/home/mjh/distribution/obs-table.fits')
t1=t1[t1['D']<2000*kpc]
t1r=t1[t1['remnant']]
t1nr=t1[~t1['remnant']]

'''
expo=0.33

best=1e100

for i in range(0,100):
    for j in range(0,100):
        lexp=0.01*j
        expo=0.01*i
        calcl=l150(t1nr['Q'],t1nr['M500'],expo,t1nr['D'],-lexp)
        diff=t1nr['l150']-calcl
        ratio=t1nr['l150']/calcl
        sqres=diff**2.0
        if np.sum(sqres)<best:
            best=np.sum(sqres)
            bvals=(lexp,expo)
        print i,lexp,expo,np.mean(diff),np.mean(ratio),np.sum(sqres)
'''

bvals=(0.6,0.5)
lexp,expo=bvals
print 'Best values are',bvals
'''
for i in range(0,100):
    zexp=0.02*i
    calcl=l150(t1nr['Q'],t1nr['M500'],expo,t1nr['D'],-lexp,t1nr['z'],-zexp)
    ratio=np.mean(t1nr['l150']/calcl)
    diff=t1nr['l150']-ratio*calcl
    sqres=diff**2.0
    print i,lexp,expo,zexp,np.mean(diff),np.mean(ratio),np.sum(sqres)

stop
'''
zexp=1.5
lexp,expo=bvals
calcl=l150(t1nr['Q'],t1nr['M500'],expo,t1nr['D'],-lexp,t1nr['z'],-zexp)
l150r=t1nr['l150']/calcl
plt.scatter(t1nr['Q'],l150r,c=t1nr['z'],alpha=0.6)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$Q$')
plt.ylabel('$L_{150}$ (W Hz$^{-1}$)')

q=np.logspace(36,40,100)
plt.plot(q,[10]*len(q),color='blue')
plt.plot(q,[0.1]*len(q),color='blue')


plt.colorbar()
plt.show()
