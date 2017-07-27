from scipy.integrate import odeint,romberg
from scipy.special import kn
import numpy as np
import pickle
import synch

from constants import *

# imported from agecorrection
def intb(st,end,time,bv,bcmb):
    '''
    Find the effective B-field experienced by the electrons
    We integrate by:
    * linearly interpolating for the first part
    * assuming 'end' is an element of 'time' and just summing over all but the first part
    The code works with B^2 (i.e. bv is in T^2, though bcmb is in T)
    '''

    s_index=np.argmax(time>=st)
    t0=time[s_index-1]
    t1=time[s_index]
    b0=bv[s_index-1]
    b1=bv[s_index]
    b=b0*(t1-st)+(0.5*(b1-b0)*(t1*t1-st*st)-(b1-b0)*t0*(t1-st))/(t1-t0)
#    print 'b (start) is',np.sqrt(b/(t1-st))
    eindex=np.argmax(time==end)
#    print s_index,eindex
    for i in range(s_index,eindex):
        b+=0.5*(bv[i+1]+bv[i])*(time[i+1]-time[i])
    b/=(end-st)
    b+=bcmb**2.0
    return np.sqrt(b)

def agecorr_findcorrection(now,freq,time,bm,bcmb,intervals=25,volumes=None,verbose=False,do_adiabatic=False):
    '''
    Find the correction factor to apply to the synchrotron luminosity.
    Optionally takes account of adiabatic expansion

    volumes gives the array of lobe volumes matching time
    '''
    
    bsm=bm**2.0
    endtime=time[now]
    if verbose: print 'endtime is',endtime,'Myr'

    emiss=np.zeros_like(freq)
    uae=np.zeros_like(freq)

    for i in range(intervals+1):

        starttime=endtime*float(i)/float(intervals)
        if verbose: print 'interval',i
        if verbose: print 'starttime is',starttime,'Myr'

        if do_adiabatic:
            vstart=np.interp(starttime,time,volumes)
        
        if (i==intervals):
            b=np.sqrt(bm[now]**2.0+bcmb**2.0)
        else:
            b=intb(starttime,endtime,time,bsm,bcmb)
        if verbose: print 'effective ageing field is',b
        for j,f in enumerate(freq):
            if freq is None:
                emiss[j]=1
                uae[j]=1
            else:
                age=(endtime-starttime)*86400.0*365.0*1.0e6
                if do_adiabatic:
                    age*=(volumes[now]/vstart)**(1.0/3.0)
                synch.setage(age,b)
                emiss[j]+=synch.emiss(1.0,bm[now],f)
                synch.setage(1.0,b)
                uae[j]+=synch.emiss(1.0,bm[now],f)

    return emiss/uae

class Evolve_RG(object):

    def px(self,x):
        return self.P0/((self.c500*x)**self.gamma)/(1+(self.c500*x)**self.alpha)**((self.beta-self.gamma)/self.alpha)

    def app(self,x):
        return 0.10-(self.alpha_p+0.10)*(x/0.5)**3.0/(1.+(x/0.5)**3.0)
    
    def upp(self,r):
        x=r/self.r500
        x=np.where(x<1e-15,1e-15,x)
        # normalization is 1.65e-3 keV/cm^3 but we need Pa
        return 1.65e-3*1e6*1e3*eV*(self.m500/3e14)**(0.66667+self.alpha_p+self.app(x))*self.px(x)

    def nupp(self,r):
        return self.upp(r)/self.kt
    
    def betam(self,r):
        return (1.0+(r/self.rc)**2.0)**(-1.5*self.beta)

    def nbetam(self,r):
        return self.n0*self.betam(r)
    def pbetam(self,r):
        return self.p0*self.betam(r)
    
    def inner_int(self,r,z,R,Rp):
        return 2.0*np.pi*r*self.nr(np.sqrt(z**2.0+r**2.0))

    def outer_int(self,z,R,Rp):
        limit=Rp*np.sqrt((1.0-z**2.0/R**2.0))
        if np.isnan(limit):
            print 'limit is nan:',z,R,Rp
            stop
        return romberg(self.inner_int,0,limit,args=(z,R,Rp),divmax=20)

    def intn(self,R,Rp):
        if R==0 or Rp==0:
            return 0
        result=romberg(self.outer_int,0,R,args=(R,Rp))
        return result

    def vtot(self,R,Rp):
        return (4.0/3.0)*np.pi*R*Rp**2.0

    def vlobe(self,R,Rp,t,N=None):
        if N is None:
            N=self.intn(R,Rp)
        return self.vtot(R,Rp)*(self.xi*self.Q*t/((2.0-self.xi)*self.Q*t+3*N*self.kt))

    def solve_mach(self,p1,p0):
        return self.cs*np.sqrt((1.0/(2.0*gamma))*((gamma+1)*(p1/p0)-(1-gamma)))

    def dL_dt(self,L,t,vl=None):
        R=L[0]
        Rp=L[1]
        if vl is None:
            vl=self.vlobe(R,Rp,t)
        if t<=self.tstop:
            internal=(self.Q*t)/(6*vl)
            ram=(self.Q*R)/(c*vl)
        else:
            internal=(self.Q*tstop)/(6*vl)
            ram=0
        result=np.array([
            self.solve_mach(ram+internal,self.pr(R)),
            self.solve_mach(internal,self.pr(Rp))
            #        solve_hybrid((Q*R)/(c*vl) + (Q*t)/(6*vl),pbetam(R),R,Rp,t),
            #        solve_hybrid((Q*t)/(6*vl),pbetam(Rp),R,Rp,t)
            ])
        result=np.where(result>c,[c,c],result)
        result=np.where(np.isnan(result),[0,0],result)
        print 'Called:',t,R,Rp,vl,result/np.sqrt(5.0*self.kt/(3.0*0.6*1.6e-27))
        return result

    def solve(self,Q,tv,tstop=None):
        self.Q=Q
        self.tv=tv
        if tstop is None:
            self.tstop=tv[-1]*2
        else:
            self.tstop=tstop
        self.results=odeint(self.dL_dt,[c*tv[0],c*tv[0]],tv)
        self.R=self.results[:,0]
        self.Rp=self.results[:,1]
        self.rlobe=[]
        self.m1=[]
        self.vl=[]
        for i in range(len(self.R)):
            vl=self.vlobe(self.R[i],self.Rp[i],tv[i])
            self.vl.append(vl)
            self.rlobe.append(vl/self.vtot(self.R[i],self.Rp[i]))
            self.m1.append(self.dL_dt([self.R[i],self.Rp[i]],tv[i],vl)[0]/self.cs)
        self.vl=np.array(self.vl)
        self.rlobe=np.array(self.rlobe)
        self.m1=np.array(self.m1)
            
    def findsynch(self,q,nu):
        self.q=q
        self.emax=1e8*m_e*c**2.0
        self.emin=10*m_e*c**2.0
        if q==2.0:
            self.I=np.log(self.emax/self.emin)
        else:
            self.I=(1.0/(2.0-q))*(self.emax**(2.0-q)-self.emin**(2.0-q))
        try:
            B=self.B
        except:
            self.findb()
            B=self.B
        # eqs from Longair
        self.synch=2.344e-25*0.4*(self.xi*self.Q*self.tv)*self.B**((q+1.0)/2.0)*(1.253e37/nu)**((q-1)/2.0)/self.I
        print self.synch
            
    def findb(self,eta=0.1):
        B=[]
        for i in range(len(self.tv)):
            vl=self.vl[i]
            E=self.Q*self.tv[i]
            U=E/vl
            B.append(np.sqrt(2*mu0*U*eta/(1+eta)))
        self.B=np.array(B)

    def findcorrection(self,freqs,z=0):
        # adapted from agecorrection.py code
        synch.setspectrum(500,1e6,self.q)

        redshift=z
        cmbtemp=2.73
        boltzmann=1.380658e-23
        planck=6.6260755e-34
        v_c=c

        bcmb=np.sqrt(8.0*np.pi**5.0*((1.0+redshift)*cmbtemp*boltzmann)**4.0/(15*planck**3.0*v_c**3.0)*(2.0*mu0))
        print 'CMB energy density in B-field terms is %g T' % bcmb

        corrs=[]
        for i in range(len(self.tv)):
            corrs.append(agecorr_findcorrection(i,freqs,self.tv/Myr,self.B,bcmb,volumes=self.vl,verbose=False,do_adiabatic=self.do_adiabatic))
        self.corrs=np.array(corrs)
            
    def setfunctions(self):
        if self.env_type=='beta':
            self.nr=self.nbetam
            self.pr=self.pbetam
        elif self.env_type=='universal':
            self.pr=self.upp
            self.nr=self.nupp
        else:
            raise Exception('env_type specified is not recognised')
        
    def __init__(self, env_type, **kwargs):
        # initialize the evolution with an environment specified by env_type
        self.env_type=env_type
        if env_type=='beta':
            self.kt=kwargs['kT']
            self.rc=kwargs['rc']
            self.beta=kwargs['beta']
            self.p0=kwargs['p0']
            self.n0=self.p0/self.kt # particles / m^3
            self.rho0=m0*self.n0 # kg/m^3
        elif env_type=='universal':
            # Universal cluster pressure profile from Arnaud et al
            mass0=3.84e14 # normalization
            self.alpha_p=0.12
            self.c500=1.177
            self.gamma=0.3081
            self.alpha=1.0510
            self.beta=5.4905
            self.P0=8.403
            self.m500=kwargs['M500']
            self.kt=5.0*(self.m500/mass0)**(1.0/1.71)
            print 'Temperature is',self.kt,'keV'
            self.kt*=1e3*eV
            self.r500=1104*kpc*(self.m500/mass0)**0.3333333
        else:
            raise Exception('env_type specified is not recognised')
        self.setfunctions()
        self.cs=np.sqrt(5.0*self.kt/(3.0*m0))
        try:
            self.xi=kwargs['xi']
            print 'Xi set to',self.xi
        except:
            self.xi=0.5

        try:
            self.do_adiabatic=kwargs['do_adiabatic']
        except:
            self.do_adiabatic=False

            
    def __getstate__(self):
        dict=self.__dict__.copy()
        del(dict['nr'])
        del(dict['pr'])
        return dict

    def __setstate__(self,d):
        self.__dict__=d
        self.setfunctions()
        
    def save(self,filename):
        f = file(filename, 'wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    @staticmethod
    def load(filename):
        with file(filename, 'rb') as f:
            return pickle.load(f)
