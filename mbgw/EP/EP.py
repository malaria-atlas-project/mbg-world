# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from copy import copy
import pymc as pm
import numpy as np
from IPython.Debugger import Pdb
import scipy
from scipy import optimize, integrate

__all__ = ['EP', 'estimate_envelopes']

def compose(*fns):
    def out(res):
        for fn in fns[::-1]:
            res = fn(res)
        return res
    return out

def d(x):
    return x[2]-x[1]

def estimate_envelopes(f, max_guess, env_guess, threshold, ftol=.001):
    m = optimize.fmin(lambda x: -f(x), x0=max_guess, ftol=ftol, full_output=False, disp=0)
    fm = -f(m)
    
    map_to_upper = lambda x: np.exp(x)+m
    map_to_lower = lambda x: m-np.exp(x)
    obj_fn = lambda x: np.abs(x - (-fm-threshold))
    
    lo = optimize.fmin(compose(obj_fn, f, map_to_lower), x0=np.log(env_guess), full_output=False, disp=0)
    hi = optimize.fmin(compose(obj_fn, f, map_to_upper), x0=np.log(env_guess), full_output=False, disp=0)
    
    # from IPython.Debugger import Pdb
    # Pdb(color_scheme='Linux').set_trace()   
    
    return map_to_lower(lo), map_to_upper(hi)

class EP(pm.Sampler):

    def __init__(self, M_pri, C_pri, lp, nug, mu_guess=None, V_guess=None):
        # Number of input points
        self.Nx = len(M_pri)
        # Approximate likelihoods will be of the form s N(mu; theta + epsilon, V)
        # where epsilon ~ N(0,nug)
        if mu_guess is None:
            self.mu = np.zeros(self.Nx,dtype=float)
        else:
            self.mu = mu_guess
        if V_guess is None:
            self.V = np.ones(self.Nx,dtype=float)*1e6
        else:
            self.V = V_guess
        # log(expected likelihood).
        self.p = np.ones(self.Nx, dtype=float)
        # Log-probability functions
        self.lp = lp
        # 'Nugget' for epsilon
        self.nug = np.resize(nug, self.Nx)
        # Prior mean and covariance of theta
        self.M_pri = M_pri
        self.C_pri = C_pri
                      
    def compute_expectation(self, i, N):
        """
        Use rejection sampling to compute moments
        """
        # TODO: This is not working well enough. Spend some cycles here to get as good
        # of an estimate as you can, otherwise your accuracy is severely limited.
        m = self.M[i]
        pri_v = self.C[i,i]
        nug_v = self.nug[i]
        v = pri_v + nug_v
        
        pri_fn = lambda x: pm.normal_like(x, m, 1./v)
        like_fn = lambda x: self.lp[i](np.atleast_1d(x)).squeeze()
        post_fn = lambda x: pri_fn(x) + like_fn(x)
        
        lo, hi = estimate_envelopes(post_fn, m, np.sqrt(v), 13.)
        x = np.linspace(lo, hi, N)
        post_vec = np.exp([post_fn(xi) for xi in x])

        p = integrate.simps(post_vec, dx=d(x))
        
        # The expectations of these functions give the first two moments of x without the nugget.
        # Just passing in lambda x:x and lambda x:x**2 would give the first two moments of x
        # _with_ the nugget.
        nonnug_m1 = lambda x: (m*nug_v + x*pri_v)/v
        nonnug_m2 = lambda x: (nug_v * pri_v)/v + nonnug_m1(x)**2
        funs = [nonnug_m1, nonnug_m2]
        
        # Return E_pri [like_fn(x)] and the posterior expectations of funs(x).        
        moments = []
        for f in funs:
            moments.append(integrate.simps(post_vec * f(x), dx=d(x)) / p)
        
        # from IPython.Debugger import Pdb
        # Pdb(color_scheme='Linux').set_trace()   
        
        return (p,) + tuple(moments)
                        
    def observe(self, i, unobserve=False):
        """
        Imposes observation
            mu[i] ~ N(theta[i], self.V[i] + nugget)
        on current distribution of theta, or 'unimposes' that
        observation if unobserve=True.
        """
        # Checked: unobserve really does invert observe
        if unobserve:   
            if self.V[i]==np.Inf:
                return
            eff_Cobs = (self.C[i,i] - self.nug[i] - self.V[i])
        else:
            eff_Cobs = (self.C[i,i] + self.nug[i] + self.V[i])
        if eff_Cobs==0:
            # Pdb(color_scheme='Linux').set_trace()
            raise RuntimeError, 'Observed covariance is zero.'
        offdiag = np.asarray(self.C[i,:]).ravel()
        new_M = self.M + offdiag / eff_Cobs* (self.mu[i] - self.M[i])
        new_C = self.C - np.outer(offdiag, offdiag) / eff_Cobs
        if np.any(np.isinf(new_M)) or np.any(np.isinf(new_M)):
            # Pdb(color_scheme='Linux').set_trace()  
            raise RuntimeError, 'Infinite covariance or mean.'
        if np.any(np.diag(new_C)<0):
            # Pdb(color_scheme='Linux').set_trace()    
            raise RuntimeError, 'Negative diagonal elements of covariance.'
        self.M = new_M
        self.C = new_C
        self.M.flags['WRITEABLE'] = False
        self.C.flags['WRITEABLE'] = False            
        
    def update_item(self, i, N):
        """
        Update approximate likelihood parameters mu[i], V[i], s[i] according to 
        EP algorithm
        """
        
        # print self.M
        # print self.C
        # print i
        # print
        
        if self.nug[i] == np.Inf:
            return
            
        self.observe(i, unobserve=True)
        self.p[i], m1, m2 = self.compute_expectation(i, N)
        
        if not np.isinf(self.p[i]):

            V_post = m2 - m1*m1
            mu_post = m1
        
            V_pri = self.C[i,i]
            mu_pri = self.M[i]
            
            V_like = V_pri*V_post / (V_pri-V_post)
            mu_like = (mu_post*(V_like+V_pri) - mu_pri*V_like)/V_pri
            self.V[i] = V_like - self.nug[i]            
            self.mu[i] = mu_like

            self.observe(i)
            
            if np.isinf(self.V[i]):
                # Pdb(color_scheme='Linux').set_trace()   
                raise RuntimeError, 'Infinite likelihood varaince.'
            if np.any(np.isnan(self.V)) | np.any(np.isnan(self.mu)):
                # Pdb(color_scheme='Linux').set_trace()
                raise RuntimeError, 'Likelihood variance or mean is nan.'
        else:
            print 'Zero probability'
            self.observe(i)
            self.p[i] = -np.Inf
        
        # print '\t\t',i, self.V[i], self.mu[i]
        
        # if self.C[i,i]<0 or self.V[i] <= 0:
        #     Pdb(color_scheme='Linux').set_trace()
        
    def update_sweep(self, N):
        for i in xrange(self.Nx):
            self.update_item(i, N)
    
    def fit(self, N, tol=.01):
        """
        Stores approximate likelihood parameters mu, V
        and returns approximation of log p(D).
        """
        
        # Make initial observations (usually trivial)
        self.C = self.C_pri.copy()
        self.M = self.M_pri.copy()
        
        for i in xrange(self.Nx):
            self.observe(i)
            if np.any(np.isinf(self.C)) or np.any(np.isinf(self.M)):
                # Pdb(color_scheme='Linux').set_trace()  
                raise RuntimeError, 'C or M is infinite.'
            if np.any(np.diag(self.C)<0):
                # Pdb(color_scheme='Linux').set_trace()   
                raise RuntimeError, 'C has negative diagonal.' 
        
        # To help convergence, use same iid normal samples
        # in all iterations.
        # self.normal_samps = np.random.normal(size=N_samps)
        # self.nug_samps = np.random.normal(size=N_samps)
        
        # Record last values of mu and V to assess convergence.
        last_mu = np.copy(self.mu)
        last_V = np.copy(self.V)
        
        sw = 0
        dmu = np.Inf
        dV = np.Inf
        # print self.mu, self.V
        while max(dmu, dV) > tol:
            sw += 1
            # Update
            self.update_sweep(N)
            # Assess convergence and possibly terminate.
            dmu = np.max(np.abs((self.mu - last_mu)/self.mu))
            if np.isnan(dmu):
                dmu = np.max(np.abs((self.mu - last_mu)*10))
            dV = np.max(np.abs((self.V - last_V)/self.V))

            if sw > 10000 and sw % 100 == 0:
                print 'Too many iterations, randomizing.'
                # Pdb(color_scheme='Linux').set_trace()                   
                print 'mu: %s'%self.mu
                print 'V: %s'%self.V
                self.V = last_V + np.random.random(size=len(self.V)) * dV * self.V
                self.mu = last_mu + np.random.random(size=len(self.mu)) * dmu * self.mu
                
            last_V[:] = self.V[:]
            last_mu[:] = self.mu[:]
        
            # print self.mu, self.V
            
        V = np.array((self.V+self.nug))
        V_ind = V + np.diag(self.C_pri)
        C_joint = self.C_pri+np.diag(V)
        if np.all(V_ind>0):
            joint_term = pm.mv_normal_cov_like(self.mu, self.M_pri, C_joint)
            ind_term = pm.normal_like(self.mu, self.M_pri, 1/V_ind)
            log_ratio = joint_term-ind_term
        # Protect against negative 'variances' (which are acceptable)
        else:
            V = V.astype('complex')
            V_ind = V_ind.astype('complex')
            C_joint = C_joint.astype('complex')

            dev=self.mu-self.M_pri            
            ind_term = -.5*np.sum((np.log(2.*np.pi*V_ind) + dev**2/V_ind))

            val,vec = np.linalg.eig(C_joint)            
            sq = np.asarray(np.dot(dev,vec/np.sqrt(val))).ravel()
            joint_term = -.5*(np.sum((np.log(2.*np.pi*val))) + np.dot(sq, sq))

            log_ratio = np.real(joint_term-ind_term)
        
        if np.sum(self.p) + log_ratio > 10000 or np.any(np.isnan(log_ratio)):
            from IPython.Debugger import Pdb
            Pdb(color_scheme='Linux').set_trace()   
        
        return np.sum(self.p) + log_ratio