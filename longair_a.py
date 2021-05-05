from __future__ import division
from scipy.special import gamma
import numpy as np

def a(p):
    a=np.sqrt(np.pi)/2.0
    a*=gamma(p/4 + 19/12)
    a*=gamma(p/4 - 1/12)
    a*=gamma(p/4 + 5/4)
    a/=(p+1)
    a/=gamma(p/4 + 7/4)
    return a

if __name__=='__main__':
    for p in [1,1.5,2,2.5,3,3.5,4,4.5,5]:
        print("%3.1f %5.3f" % (p, a(p)))
        
    
