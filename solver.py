from __future__ import print_function
from scipy.integrate import odeint,romberg
from scipy.special import kn
import numpy as np
import pickle
import synch

from synch_constants import *
from longair_a import a as longair_a

def findbcmb(redshift):
    cmbtemp=2.73
    v_c=c
    return np.sqrt(8.0*np.pi**5.0*((1.0+redshift)*cmbtemp*boltzmann)**4.0/(15*planck**3.0*v_c**3.0)*(2.0*mu0))

# imported from agecorrection
def intb(st,end,time,bv,bcmb):
    '''
    Find the effective B-field experienced by the electrons between
    time 'st' and time 'end' given an array of times 'time' and an
    array of B-field values 'bv'.

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

def loss_findcorrection(now,time,bm,bcmb,intervals=25,volumes=None,verbose=False,do_adiabatic=False,tstop=None):
    
    bsm=bm**2.0
    endtime=time[now]
    if verbose: print('endtime is',endtime,'Myr')

    loss=0
    ual=0
    for i in range(intervals+1):

        starttime=endtime*float(i)/float(intervals)
        if tstop is not None and starttime>tstop:
            starttime=tstop
        if verbose: print('interval',i)
        if verbose: print('starttime is',starttime,'Myr')

        if do_adiabatic:
            vstart=np.interp(starttime,time,volumes)
        
        if (i==intervals):
            b=np.sqrt(bm[now]**2.0+bcmb**2.0)
        else:
            b=intb(starttime,endtime,time,bsm,bcmb)
        if verbose: print('effective ageing field is',b)
        age=(endtime-starttime)*86400.0*365.0*1.0e6
        if do_adiabatic:
            age*=(volumes[now]/vstart)**(1.0/3.0)
            # trapezium rule factors
        if i==0 or i==intervals:
            factor=0.5
        else:
            factor=1.0
        synch.setage(age,b)
        loss+=factor*synch.loss(1.0,bm[now])
        synch.setage(0,b)
        ual+=factor*synch.loss(1.0,bm[now])

    return loss/ual

def agecorr_findcorrection(now,freq,time,bm,bcmb,intervals=25,volumes=None,verbose=False,do_adiabatic=False,tstop=None):
    '''
    Find the correction factor to apply to the synchrotron luminosity.
    Optionally takes account of adiabatic expansion

    'now' is the time index at which we want to compute the correction
    'freq' is an array or list of frequency values in Hz -- we loop over these
    'time' is an array of times (in Myr) at which we have magnetic field strengths calculated
    'bm' are the calculated magnetic field strengths in T
    'intervals' is the number of intervals to break each time range into
    'volumes' gives the array of lobe volumes matching time, only used if 'do_adiabatic' is True
    'tstop' is the stop time in Myr

    This code is called with various values of 'now' and assumes by
    default that the electron population is composed of equal amounts
    of particles injected between t=0 and 'now' -- i.e. 'starttime' varies
    from 0 (particles injected at switchon) to 'endtime' (particles injected
    right now).

    To deal with the situation where there is a tstop we simply prevent starttime from getting any larger than tstop. So no particles younger than this are injected.

    '''
    
    bsm=bm**2.0
    endtime=time[now]
    if verbose: print('endtime is',endtime,'Myr')

    emiss=np.zeros_like(freq)
    uae=np.zeros_like(freq)
    for i in range(intervals+1):

        starttime=endtime*float(i)/float(intervals)
        if tstop is not None and starttime>tstop:
            starttime=tstop
        if verbose: print('interval',i)
        if verbose: print('starttime is',starttime,'Myr')

        if do_adiabatic:
            vstart=np.interp(starttime,time,volumes)
        
        if (i==intervals):
            b=np.sqrt(bm[now]**2.0+bcmb**2.0)
        else:
            b=intb(starttime,endtime,time,bsm,bcmb)
        if verbose: print('effective ageing field is',b)
        for j,f in enumerate(freq):
            if f is None:
                emiss[j]=1
                uae[j]=1
            else:
                age=(endtime-starttime)*86400.0*365.0*1.0e6
                if do_adiabatic:
                    age*=(volumes[now]/vstart)**(1.0/3.0)
                # trapezium rule factors
                if i==0 or i==intervals:
                    factor=0.5
                else:
                    factor=1.0
                synch.setage(age,b)
                emiss[j]+=factor*synch.emiss(1.0,bm[now],f)
                synch.setage(1.0,b)
                uae[j]+=factor*synch.emiss(1.0,bm[now],f)

    return emiss/uae

def ic_agecorr_findcorrection(now,freq,z,time,bm,bcmb,intervals=40,volumes=None,verbose=False,do_adiabatic=False,tstop=None):

    bsm=bm**2.0
    endtime=time[now]
    if verbose: print('endtime is',endtime,'Myr')

    emiss=np.zeros_like(freq)
    uae=np.zeros_like(freq)
    zeroage=np.zeros_like(freq)
    synch.setage(1.0,1e-9)
    for j,f in enumerate(freq):
        zeroage[j]=synch.cmb_ic_emiss(1.0,f,z)

    for i in range(intervals+1):

        starttime=endtime*float(i)/float(intervals)
        if tstop is not None and starttime>tstop:
            starttime=tstop
        if verbose: print('interval',i)
        if verbose: print('starttime is',starttime,'Myr')

        if do_adiabatic:
            vstart=np.interp(starttime,time,volumes)
        
        if (i==intervals):
            b=np.sqrt(bm[now]**2.0+bcmb**2.0)
        else:
            b=intb(starttime,endtime,time,bsm,bcmb)
        if verbose: print('effective ageing field is',b)
        for j,f in enumerate(freq):
            if f is None:
                emiss[j]=1
                uae[j]=1
            else:
                age=(endtime-starttime)*86400.0*365.0*1.0e6
                if do_adiabatic:
                    age*=(volumes[now]/vstart)**(1.0/3.0)
                # trapezium rule factors
                if i==0 or i==intervals:
                    factor=0.5
                else:
                    factor=1.0
                synch.setage(age,b)
                emiss[j]+=factor*synch.cmb_ic_emiss(1.0,f,z)
                uae[j]+=factor*zeroage[j]

    return emiss/uae

def vcomb(v,L):
    R=L[0]
    Rp=L[1]
    return np.sqrt(((Rp/R)*v[0]**2+v[1]**2)/(1+(Rp/R)))

class Evolve_RG(object):
    '''
    Class which sets up and runs the radio galaxy evolution model.
    Methods generally update attributes, see method docstrings for details.
    Public methods include:
    solve             -- solve for the dynamics
    findb             -- find magnetic field strength
    findsynch         -- find synchrotron emission (uncorrected)
    findic            -- find inverse-Compton emission (uncorrected)
    findcorrection    -- find loss corrections for synchrotron
    ic_findcorrection -- find loss corrections for inverse-Compton
    finds_loss        -- find total synchrotron loss rate (uncorrected)
    findic_loss       -- find total IC loss rate (uncorrected)
    findlosscorrection-- find correction factor for losses in synch/IC
    find_shellt       -- find physical conditions in the shocked shell
    findbs_loss       -- find thermal bremsstrahlung losses from shocked shell
    save              -- save the current state
    load              -- load a previously saved state (unbound method)
    '''
    def _px(self,x):
        return self.P0/((self.c500*x)**self.gamma)/(1+(self.c500*x)**self.alpha_u)**((self.beta-self.gamma)/self.alpha_u)

    def _app(self,x):
        return 0.10-(self.alpha_p+0.10)*(x/0.5)**3.0/(1.+(x/0.5)**3.0)
    
    def _upp(self,r):
        x=r/self.r500
        x=np.where(x<1e-15,1e-15,x)
        # normalization is 1.65e-3 keV/cm^3 but we need Pa
        return 1.65e-3*1e6*1e3*eV*(self.m500/3e14)**(0.66667+self.alpha_p+self._app(x))*self._px(x)

    def _nupp(self,r):
        return self._upp(r)/self.kt

    def _upp_floor(self,r):
        p=self._upp(r)
        if self.additive_floor:
            return p+self.pfloor
        return max(p,self.pfloor)

    def _nupp_floor(self,r):
        n=self._nupp(r)
        if self.additive_floor:
            return n+self.nfloor
        return max(n,self.nfloor)
    
    def _betam(self,r):
        return (1.0+(r/self.rc)**2.0)**(-1.5*self.beta)

    def _nbetam(self,r):
        return self.n0*self._betam(r)
    def _pbetam(self,r):
        return self.p0*self._betam(r)
    
    def _inner_int(self,r,z,R,Rp):
        return 2.0*np.pi*r*self.nr(np.sqrt(z**2.0+r**2.0))

    def _outer_int(self,z,R,Rp):
        limit=Rp*np.sqrt((1.0-z**2.0/R**2.0))
        if np.isnan(limit):
            print('limit is nan:',z,R,Rp)
            raise RuntimeError('Integration called with impossible bounds')
        return romberg(self._inner_int,0,limit,args=(z,R,Rp),divmax=20)

    def intn(self,R,Rp):
        ''' Compute number of particles in the lobe given R, Rp '''
        if R<1 or Rp<1:
            return 0
        # factor 2 here since we integrate only one lobe
        result=2.0*romberg(self._outer_int,0,R,args=(R,Rp))
        return result

    def vtot(self,R,Rp):
        '''
        Compute total volume of the shock given R, Rp
        '''
        return (4.0/3.0)*np.pi*R*Rp**2.0

    def vlobe(self,R,Rp,v,t):
        '''
        Compute volume of the lobe given R, Rp, v, t
        '''
        if (R,Rp) in self.ndict:
            N=self.ndict[(R,Rp)]
        else:
            N=self.intn(R,Rp)
            self.ndict[(R,Rp)]=N
        if t>self.tstop:
            time=self.tstop
        else:
            time=t

        # change in pressure_issue
        r=((self.Gamma_j-1)*self.xi*self.Q*time/
           ((self.xi*self.Gamma_j + (1-self.xi)*self.Gamma_s - 1)*self.Q*time +
            N*self.kt - (self.Gamma_s - 1)*N*m0*((v-self.cs)**2.0)/2.0))
        #if r<0:
        #    print 'Warning: vlobe is 0 at t=%g!' % t
        #    r=1e-3 # avoid nans
        #if r>1: r=1
        return self.vtot(R,Rp)*r

    def _solve_rel(self,X):
    # solve the equation Gamma v^2 = X for v
        return np.sqrt(2*np.sqrt(X**4.0+4*c**4.0*X**2.0)-2*X**2.0)/(2.0*c)

    def _solve_rpb(self,pint,pext,dext):
        return self._solve_rel((pint-pext)/dext)

    def _rhp(self,p1,p0):
        return np.sqrt((1.0/(2.0*self.Gamma_s))*((self.Gamma_s+1)*(p1/p0)-(1-self.Gamma_s)))
    
    def _solve_mach(self,p1,p0,do_raise=False):
        # Solve for growth of shocked shell unless pressure is too
        # low. If it is, allow expansion at sound speed to stabilize
        # situation.
        if p1<p0:
            if do_raise:
                raise RuntimeError('Internal pressure has fallen below external pressure')
            #else:
            #    if self.verbose: print('Warning: internal pressure %g has fallen below external pressure %g' % (p1,p0))
            return self.cs
        else:
            return self.cs*self._rhp(p1,p0)

    def dL_dt_vest(self,L,t,v_est):
        R=L[0]
        Rp=L[1]

        vl=self.vlobe(R,Rp,v_est,t)
        self.tempvl=vl
        if t<=self.tstop:
            internal=self.prfactor*(self.xi*self.Q*t)/vl
            ram=self.epsilon*(self.Q*R)/(2*self.qfactor*vl)
        else:
            internal=self.prfactor*(self.xi*self.Q*self.tstop)/vl
            ram=0
        #if self.verbose:
        #    print('dL_dt_Pressures:',ram,internal,self.pr(R),self.pr(Rp))
        try:
            result=np.array([
                self._solve_mach(ram+internal,self.pr(R),do_raise=False),
                self._solve_mach(internal,self.pr(Rp))
                ])
        except RuntimeError:
            print(R,Rp,t,v_est,vl,internal,ram,self.pr(R),self.pr(Rp))
            raise
        result=np.where(result>c,[c,c],result)
        result=np.where(np.isnan(result),[0,0],result)

        return result

    def iter_dLdt(self,L,t,iterlimit=100):
        '''
        Attempt to find a self-consistent velocity solution
        dLdt takes an estimated speed for ke and returns the velocity vector
        iterative solving works by bisection, i.e. we are trying to solve for the difference between the input and output velocity being zero
        '''
        v0=np.array([self.cs,self.cs])
        r0=self.dL_dt_vest(L,t,vcomb(v0,L))
        v2=r0
        r2=self.dL_dt_vest(L,t,vcomb(v2,L))
        v1=(v0+v2)/2.0
        r1=self.dL_dt_vest(L,t,vcomb(v1,L))
        d0=vcomb(v0,L)-vcomb(r0,L)
        d1=vcomb(v1,L)-vcomb(r1,L)
        d2=vcomb(v2,L)-vcomb(r2,L)

        iter=0
        while iter<iterlimit:
            #if self.verbose:
            #    print iter,v1,r1,d1
            if np.sign(d1)==np.sign(d0):
                # new midpoint between 1 and 2
                v0=v1
                r0=r1
                d0=d1
            else:
                v2=v1
                r2=r1
                d2=d1
            v1=(v0+v2)/2.0
            v1=np.where(v1<self.cs,self.cs,v1)
            #print v0,v1,v2
            d1_old=d1
            r1=self.dL_dt_vest(L,t,vcomb(v1,L))
            d1=vcomb(v1,L)-vcomb(r1,L)
            if iter>10 and (np.sum(np.abs(d1_old-d1))/np.sum(v1))<1e-8:
                break

            iter+=1
        if iter==iterlimit:
            print('dLdt: ',t,L,r1,iter)
            raise RuntimeError('Convergence failed')

        if self.verbose:
            print('dLdt: returning:',t,L,r1,iter)
        return r1

    
    def iter_dLdt_old(self,L,t,verbose=False,iterlimit=100):
        '''
        Attempt to find a self-consistent velocity solution
        dLdt takes an estimated speed for ke and returns the velocity vector
        iterative solving works by bisection, i.e. we are trying to minimize
        the difference between the input and the output velocity
        '''
        v0=np.array([0,0])
        r0=self.dL_dt_vest(L,t,vcomb(v0))
        v2=r0
        r2=self.dL_dt_vest(L,t,vcomb(v2))
        v1=(v0+v2)/2.0
        r1=self.dL_dt_vest(L,t,vcomb(v1))
        d0=abs(vcomb(v0)-vcomb(r0))
        d1=abs(vcomb(v1)-vcomb(r1))
        d2=abs(vcomb(v2)-vcomb(r2))

        iter=0
        while iter<iterlimit:
            if verbose: print(iter,v1,r1,d1)
            r1_old=r1
            if d0>d2:
                # new midpoint between 1 and 2
                v0=v1
                r0=r1
                d0=d1
            else:
                v2=v1
                r2=r1
                d2=d1
            v1=(v0+v2)/2.0
            r1=self.dL_dt_vest(L,t,vcomb(v1))
            if (np.sum(np.abs(r1_old-r1))/np.sum(r1))<1e-6:
                break
            d1=abs(vcomb(v1)-vcomb(r1))
            
            iter+=1
        if iter==iterlimit:
            print('dLdt: ',t,L,r1,iter)
            raise RuntimeError('Convergence failed')
        if verbose:
            print('dLdt: returning:',t,L,r1,iter)
        return r1

    def solve(self,Q,tv,tstop=None):
        '''
        Solve for the dynamics of the source.
        Parameters are:
        Q     -- the jet power (W)
        tv    -- a sorted numpy array of times to solve for (s)
        tstop -- a time when the jet switches off. If unset, defaults larger than max(tv)
        Updated attributes:
        tv    -- the time values
        R     -- lobe length (m)
        Rp    -- lobe width (m)
        m1    -- forward Mach number
        mp1   -- transverse Mach number
        vl    -- lobe volume
        vt    -- total volume
        '''
        self.Q=Q
        self.tv=tv
        if tstop is None:
            self.tstop=tv[-1]*2
        else:
            self.tstop=tstop
        if self.verbose:
            print('tstop is',self.tstop)
        self.ndict={}
        self.results=odeint(self.iter_dLdt,[c*tv[0],c*tv[0]],tv)
        self.R=self.results[:,0]
        self.Rp=self.results[:,1]
        # now redo to find speeds etc
        self.rlobe=[]
        self.m1=[]
        self.mp1=[]
        self.vl=[]
        self.vt=[]
        for i in range(len(self.R)):
            vt=self.vtot(self.R[i],self.Rp[i])
            self.vt.append(vt)
            speeds=self.iter_dLdt([self.R[i],self.Rp[i]],tv[i])
            vl=self.tempvl # stored here by dLdt
            self.vl.append(vl)
            self.rlobe.append(vl/vt)
            self.m1.append(speeds[0]/self.cs)
            self.mp1.append(speeds[1]/self.cs)
        self.vl=np.array(self.vl)
        self.vt=np.array(self.vt)
        self.rlobe=np.array(self.rlobe)
        self.m1=np.array(self.m1)
        self.mp1=np.array(self.mp1)

    def find_shellt(self):
        '''Find the shell temperature and other useful physical parameters as
        a function of time, assuming solve has been run. Values before
        and around the transition from relativistic to
        non-relativistic values should not be relied on.
        Updated attributes:
        ns -- an array matching tv of the total particle number N in the shell
        es -- an array matching tv of the total thermal energy in the shell (W)
        ts -- an array matching tv of the temperature in the shell (K)

        '''
        ns=[]
        es=[]
        ts=[]
        speeds=[]
        times=np.where(self.tv<self.tstop,self.tv,self.tstop)
        for i in range(len(self.R)):
            # compute the thermal energy in the shocked shell
            N=self.ndict[(self.R[i],self.Rp[i])]
            # hack for the KE to not count shells expanding at sound speed
            speed=self.cs*(np.array([self.m1[i],self.mp1[i]])-1)
            E = (1-self.xi)*self.Q*times[i] + (1.0/(self.Gamma_s-1))*N*self.kt - 0.5*N*m0*vcomb(speed,[self.R[i],self.Rp[i]])**2.0
            T = (self.Gamma_s-1)*(E/N)/boltzmann
            ns.append(N)
            es.append(E)
            ts.append(T)
        self.ns=np.array(ns)
        self.es=np.array(es)
        self.ts=np.array(ts)
            
    def _init_synch(self):
        ''' Set up the synchrotron parameters '''
        q=self.q
        self.alpha=0.5*(q-1)
        self.emax=self.gmax*m_e*c**2.0
        self.emin=self.gmin*m_e*c**2.0

        if q==2.0:
            self.I=np.log(self.emax/self.emin)
        else:
            self.I=(1.0/(2.0-q))*(self.emax**(2.0-q)-self.emin**(2.0-q))

        if q==3.0:
            self.Iloss=np.log(self.emax/self.emin)
        else:
            self.Iloss=(1.0/(3.0-q))*(self.emax**(3.0-q)-self.emin**(3.0-q))

    def finds_loss(self):
        '''
        Compute the total synchrotron loss rate, neglecting the effect of
        losses that have already happened
        Updated attributes:
        loss -- The synchrotron power (W)

        '''
        
        if self.Gamma_j!=4.0/3.0:
            raise NotImplementedError('Jet fluid adiabatic index is not 4/3: findsynch assumes a relativistic fluid')

        self._init_synch()
        
        try:
            B=self.B
        except AttributeError:
            self.findb()
            B=self.B

        times=np.where(self.tv<self.tstop,self.tv,self.tstop)

        self.loss=self.Iloss*(self.xi*self.Q*times)*(B**2.0/(2*mu0))*4*sigma_T/(3*m_e**2.0*c**3.0*(1+self.zeta+self.kappa)*self.I)

    def findic_loss(self,z=None):
        '''
        Compute the total IC loss rate, neglecting the effect of
        losses that have already happened
        Updated attributes:
        ic_loss -- The IC power (W)

        '''

        if self.Gamma_j!=4.0/3.0:
            raise NotImplementedError('Jet fluid adiabatic index is not 4/3: findsynch assumes a relativistic fluid')

        self._init_synch()

        if z is None:
            z=self.z

        bcmb=findbcmb(z)
        times=np.where(self.tv<self.tstop,self.tv,self.tstop)

        self.ic_loss=self.Iloss*(self.xi*self.Q*times)*(bcmb**2.0/(2*mu0))*4*sigma_T/(3*m_e**2.0*c**3.0*(1+self.zeta+self.kappa)*self.I)

    def findbs_loss(self):
        '''
        Compute the total bremsstrahlung loss rate, neglecting the effect of
        losses that have already happened
        Updated attributes:
        bs_loss -- The thermal bremsstrahlung power (W)

        '''
        try:
            self.ts
        except AttributeError:
            self.find_shellt()
        vs=self.vt-self.vl
        P=1.2*1.44e-40*(1.18/2.18**2.0)*self.ns**2.0*self.ts**(0.5)/vs
        self.bs_loss=P

            
    def findsynch(self,nu):
        '''
        Compute synchrotron emission.
        Parameters:
        nu -- a reference frequency (Hz)
        Updated attributes:
        nu_ref -- the reference frequency
        synch  -- the uncorrected synchrotron luminosity (W/Hz)
        '''
        if self.Gamma_j!=4.0/3.0:
            raise NotImplementedError('Jet fluid adiabatic index is not 4/3: findsynch assumes a relativistic fluid')
        self.nu_ref=nu
        self._init_synch()
        
        try:
            B=self.B
        except AttributeError:
            self.findb()
            B=self.B

        times=np.where(self.tv<self.tstop,self.tv,self.tstop)

        # Longair 2010 eq. 8.130
        self.synch=2.344e-25*longair_a(self.q)*(self.xi*self.Q*times)*self.B**((self.q+1.0)/2.0)*(1.253e37/nu)**((self.q-1)/2.0)/((1+self.zeta+self.kappa)*self.I)

    def findic(self,nu,z=None):
        '''
        Compute inverse-Compton emission
        Parameters:
        nu -- a reference frequency (Hz)
        z  -- the redshift (if unset, use default set on initialization)
        Updated attributes:
        nu_ic_ref -- the reference frequency
        ic        -- the uncorrected inverse-Compton luminosity (W/Hz)
        '''
        if z is None:
            z=self.z
        else:
            self.z=z
            print('findic over-riding previously set z to %f' % z)
        self.nu_ic_ref=nu
        self._init_synch()
        
        times=np.where(self.tv<self.tstop,self.tv,self.tstop)

        self.ic=(self.xi*self.Q*times)/((1+self.zeta+self.kappa)*self.I)
        synch.setspectrum(self.gmin,self.gmax,self.q)
        self.ic*=synch.cmb_ic_emiss(1,nu,z)

    def findb(self):
        '''
        Compute magnetic field strength in the lobes
        Updated attributes:
        B -- The magnetic field strength as a function of time
        '''
        B=[]
        for i in range(len(self.tv)):
            vl=self.vl[i]
            if self.tv[i]>self.tstop:
                t=self.tstop
            else:
                t=self.tv[i]
            E=self.xi*self.Q*t
            U=E/vl
            B.append(np.sqrt(2*mu0*U*self.zeta/(1+self.zeta+self.kappa)))
        self.B=np.array(B)

    def findlosscorrection(self,z=None,do_adiabatic=None,timerange=None):
        '''
        Find corrections to the simple synchrotron loss formula.
        Parameters:
        z            -- the redshift. If undefined defaults to the value already chosen
        do_adiabatic -- boolean specifying whether adiabatic corrections should be done
        timerange    -- a list of the _indices_ of the time values to find the correction for. If unset, do all of them.
        Updated attributes:
        losscorrs    -- the correction factors per frequency (times, freqs)
        '''

        try:
            self.loss
        except AttributeError:
            print('Warning: finds_loss not previously run, running it now')
            self.finds_loss()

        if z is None:
            z=self.z
        else:
            self.z=z
            print('findcorrection over-riding previously set z to %f' % z)

        # adapted from agecorrection.py code
        if timerange is None:
            timerange=range(len(self.tv))
        synch.setspectrum(self.gmin,self.gmax,self.q)
        if do_adiabatic is not None:
            print('Over-riding do_adiabatic setting to',do_adiabatic)
            self.do_adiabatic=do_adiabatic

        bcmb=findbcmb(z)
        if self.verbose:
            print('CMB energy density in B-field terms is %g T' % bcmb)

        corrs=np.ones_like(self.tv)*np.nan
        if self.verbose:
            print('Finding correction factors:')
        for i in timerange:
            if self.verbose: print(self.tv[i])
            corrs[i]=(loss_findcorrection(i,self.tv/Myr,self.B,bcmb,volumes=self.vl,verbose=False,do_adiabatic=self.do_adiabatic,tstop=self.tstop/Myr))
        self.losscorrs=corrs

    def findcorrection(self,freqs,z=None,do_adiabatic=None,timerange=None,intervals=25):
        '''
        Find corrections to the simple synchrotron formula.
        Parameters:
        freqs        -- a list of frequencies at which the correction should be calculated
        z            -- the redshift. If undefined defaults to the value already chosen
        do_adiabatic -- boolean specifying whether adiabatic corrections should be done
        timerange    -- a list of the _indices_ of the time values to find the correction for. If unset, do all of them.
        Updated attributes:
        freqs        -- the frequencies at which corrections were found
        corrs        -- the correction factors per frequency (times, freqs)
        corr_synch   -- the corrected synchrotron luminosity density (W/Hz)
        '''

        try:
            self.synch
        except AttributeError:
            print('Warning: findsynch not previously run, running it now')
            self.findsynch(nu=freqs[0])

        if z is None:
            z=self.z
        else:
            self.z=z
            print('findcorrection over-riding previously set z to %f' % z)

        # adapted from agecorrection.py code
        if timerange is None:
            timerange=range(len(self.tv))
        synch.setspectrum(self.gmin,self.gmax,self.q)
        if do_adiabatic is not None:
            print('Over-riding do_adiabatic setting to',do_adiabatic)
            self.do_adiabatic=do_adiabatic
        self.freqs=freqs

        bcmb=findbcmb(z)
        if self.verbose:
            print('CMB energy density in B-field terms is %g T' % bcmb)

        corrs=np.ones((len(self.tv),len(freqs)))*np.nan
        if self.verbose:
            print('Finding correction factors:')
        for i in timerange:
            if self.verbose: print(self.tv[i])
            corrs[i]=(agecorr_findcorrection(i,freqs,self.tv/Myr,self.B,bcmb,volumes=self.vl,verbose=False,do_adiabatic=self.do_adiabatic,tstop=self.tstop/Myr,intervals=intervals))
        self.corrs=corrs
        cs=np.zeros_like(corrs)
        for i,f in enumerate(self.freqs):
            cs[:,i]=(self.synch*self.corrs[:,i]*(f/self.nu_ref)**-self.alpha)
        self.corr_synch=cs
                      
    def ic_findcorrection(self,freqs,z=None,do_adiabatic=None,timerange=None):
        '''
        Find corrections to the simple inverse-Compton formula.
        Parameters:
        freqs        -- a list of frequencies at which the correction should be calculated
        z            -- the redshift. If undefined defaults to the value already chosen
        do_adiabatic -- boolean specifying whether adiabatic corrections should be done
        timerange    -- a list of the _indices_ of the time values to find the correction for. If unset, do all of them.
        Updated attributes:
        ic_freqs     -- the frequencies at which corrections were found
        ic_corrs     -- the inverse-Compton correction factors per frequency (times, freqs)
        corr_ic      -- the corrected inverse-Compton luminosity density (W/Hz)
        '''

        try:
            self.ic
        except AttributeError:
            print('Warning: findic not previously run, running it now')
            self.findic(nu=freqs[0])


        if z is None:
            z=self.z
        else:
            self.z=z
            print('ic_findcorrection over-riding previously set z to %f' % z)
        if timerange is None:
            timerange=range(len(self.tv))
        synch.setspectrum(self.gmin,self.gmax,self.q)
        if do_adiabatic is not None:
            print('Over-riding do_adiabatic setting to',do_adiabatic)
            self.do_adiabatic=do_adiabatic
        self.ic_freqs=freqs
        bcmb=findbcmb(z)
        if self.verbose:
            print('CMB energy density in B-field terms is %g T' % bcmb)

        corrs=np.ones((len(self.tv),len(freqs)))*np.nan
        if self.verbose:
            print('Finding correction factors:')
        for i in timerange:
            if self.verbose: print(self.tv[i])
            corrs[i]=(ic_agecorr_findcorrection(i,freqs,z,self.tv/Myr,self.B,bcmb,volumes=self.vl,verbose=False,do_adiabatic=self.do_adiabatic,tstop=self.tstop/Myr))
        self.ic_corrs=corrs
        cs=np.zeros_like(corrs)
        for i,f in enumerate(self.ic_freqs):
            cs[:,i]=(self.ic*self.ic_corrs[:,i]*(f/self.nu_ic_ref)**-self.alpha)
        self.corr_ic=cs
           
    def _setfunctions(self):
        if self.env_type=='beta':
            self.nr=self._nbetam
            self.pr=self._pbetam
        elif self.env_type=='universal':
            self.pr=self._upp
            self.nr=self._nupp
        elif self.env_type=='universal_floor':
            self.pr=self._upp_floor
            self.nr=self._nupp_floor
        else:
            raise Exception('env_type specified is not recognised')
        
    def __init__(self, env_type, **kwargs):
        '''
        Keyword arguments:
        env_type     -- What type of environment is required. Currently either 'beta' or 'universal' is allowed.
        The following keywords have sensible defaults (see the paper):
        xi           -- The fraction of energy left in the lobes
        zeta         -- Ratio between field and electron energy density
        kappa        -- Ratio between non-radiating and electron energy density
        epsilon      -- Geometrical factor for momentum flux
        qfactor      -- Energy-momentum conversion factor (m/s)
        Gamma_j        -- Adiabatic index of the lobes and jets. Defaults to 4/3
        Gamma_s        -- Adiabatic index of the shocks. Defaults to 5/3
        do_adiabatic -- Boolean controlling whether adiabatic corrections are applied
        gmin         -- Minimum Lorentz factor of the lobe electrons
        gmax         -- Maximum Lorentz factor of the lobe electrons
        q            -- Injection energy index
        z            -- Source redshift
        verbose      -- Be verbose. Set True to print diagnostics
        Depending on the choice of environment you will need additional keywords.
        For a beta model:
        kT           -- Boltzmann's constant x temperature (J)
        rc           -- Core radius (m)
        beta         -- Beta parameter
        p0           -- Central pressure (Pa)
        For a universal model:
        M500         -- The mass of the environment (solar masses)
        For a universal model with density floor:
        floor        -- minimum number density (particles/m^3)
        
        '''
        keywords = (('xi','Energy fraction in lobes',0.5),
                    ('zeta', 'Magnetic field/electron energy ratio', 0.1),
                    ('kappa', 'Non-radiating particle/electron energy ratio', 0),
                    ('epsilon', 'Geometrical factor for momentum flux', 4),
                    ('qfactor', 'Energy/momentum conversion factor', c),
                    ('Gamma_j', 'Adiabatic index of jet/lobe fluid', 4.0/3.0),
                    ('Gamma_s', 'Adiabatic index of shocked material', 5.0/3.0),
                    ('do_adiabatic', 'Perform the adiabatic corrections in synchrotron and inverse-Compton calculations', True),
                    ('gmin', 'Minimum Lorentz factor of injected electrons', 10),
                    ('gmax', 'Maximum Lorentz factor of injected electrons', 1e6),
                    ('q', 'Power-law index for injected electrons', 2.0),
                    ('z', 'Source redshift', 0),
                    ('verbose', 'Print lots of stuff out', False))
        
        # initialize the evolution with an environment specified by env_type
        self.env_type=env_type
        if env_type=='beta':
            self.kt=kwargs['kT']
            self.rc=kwargs['rc']
            self.beta=kwargs['beta']
            self.p0=kwargs['p0']
            self.n0=self.p0/self.kt # particles / m^3
        elif 'universal' in env_type:
            # Universal cluster pressure profile from Arnaud et al
            mass0=3.84e14 # normalization
            self.alpha_p=0.12
            self.c500=1.177
            self.gamma=0.3081
            self.alpha_u=1.0510
            self.beta=5.4905
            self.P0=8.403
            self.m500=kwargs['M500']
            self.kt=5.0*(self.m500/mass0)**(1.0/1.71)
            if 'verbose' in kwargs and kwargs['verbose']:
                print('Temperature is',self.kt,'keV')
            self.kt*=1e3*eV
            self.r500=1104*kpc*(self.m500/mass0)**0.3333333
            if 'floor' in env_type:
                self.nfloor=kwargs['floor']
                self.pfloor=self.nfloor*self.kt
                if 'additive_floor' in kwargs:
                    self.additive_floor=kwargs['additive_floor']
                    if 'verbose' in kwargs and kwargs['verbose']:
                        print('Setting additive floor to',self.additive_floor)

        self._setfunctions() # raises exception if env_type is not known.

        for k,desc,default in keywords:
            try:
                value=kwargs[k]
            except KeyError:
                value=default
            self.__dict__[k]=value
            if 'verbose' in kwargs and kwargs['verbose']:
                print('%s (%s) is' % (k,desc),value)

        self.cs=np.sqrt(self.Gamma_s*self.kt/m0)
        self.prfactor=self.Gamma_j-1.0
            
    def __getstate__(self):
        # When being pickled, we need to remove the references to
        # functions as they are not picklable. We will restore on unpickle.
        dict=self.__dict__.copy()
        del(dict['nr'])
        del(dict['pr'])
        return dict

    def __setstate__(self,d):
        self.__dict__=d
        self._setfunctions()
        
    def save(self,filename):
        '''
        Save the current state of the object.
        Parameters:
        filename -- a filename to save to
        '''
        with open(filename, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filename):
        '''
        Load a previously saved Evolve_RG object.
        Parameters:
        filename -- the name of a file to load
        '''
        with open(filename, 'rb') as f:
            return pickle.load(f,encoding='latin1')
