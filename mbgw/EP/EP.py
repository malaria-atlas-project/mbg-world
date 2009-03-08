# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from copy import copy
import pymc as pm
import numpy as np
from IPython.Debugger import Pdb

__all__ = ['EP']

class EP(pm.Sampler):

    def __init__(self, M_pri, C_pri, lp, nug):
        # Number of input points
        self.Nx = len(M_pri)
        # Approximate likelihoods will be of the form s N(mu; theta + epsilon, V)
        # where epsilon ~ N(0,nug)
        self.mu = np.zeros(self.Nx,dtype=float)
        self.V = np.ones(self.Nx,dtype=float)*1e6
        # log(expected likelihood).
        self.p = np.ones(self.Nx, dtype=float)
        # Log-probability functions
        self.lp = lp
        # 'Nugget' for epsilon
        self.nug = np.resize(nug, self.Nx)
        # Prior mean and covariance of theta
        self.M_pri = M_pri
        self.C_pri = C_pri
                      
    def compute_expectation(self, i, funs=[lambda x:x, lambda x:x**2]):
        """
        Use rejection sampling to compute moments
        """
        # Pdb(color_scheme='Linux').set_trace()   
        thetas = self.M[i] + np.sqrt(self.C[i,i]) * self.normal_samps
        nug_thetas = thetas + self.nug_samps * np.sqrt(self.nug[i])
        ts = self.lp[i](nug_thetas)
        # ts = np.array([self.lp[i](theta) for theta in nug_thetas]).ravel()
        p_tot = pm.flib.logsum(ts[np.where(1-np.isinf(ts))])
        w = np.exp(ts - p_tot)
        moments = [np.sum(fun(thetas)*w) for fun in funs]
        if np.any(np.isinf(moments)):
            # Pdb(color_scheme='Linux').set_trace()   
            raise RuntimeError, 'Infinite moments.'
        if np.any(np.isnan(moments)):
            # Pdb(color_scheme='Linux').set_trace()   
            raise RuntimeError, 'Some moments are nan.'
        if np.abs(moments[1]-moments[0]**2) < 1e-6:
            # Prior hardly supports posterior, return p_tot=0 and don't update moments.
            p_tot = -np.Inf
            moments = [self.M[i], self.C[i,i]]
        return tuple([p_tot]+moments)
                        
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
        
    def update_item(self, i):
        """
        Update approximate likelihood parameters mu[i], V[i], s[i] according to 
        EP algorithm
        """
        if self.nug[i] == np.Inf:
            return
            
        self.observe(i, unobserve=True)
        self.p[i], m1, m2 = self.compute_expectation(i)
        
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
        
    def update_sweep(self):
        for i in xrange(self.Nx):
            self.update_item(i)
    
    def fit(self, N_samps, tol=.01):
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
        self.normal_samps = np.random.normal(size=N_samps)
        self.nug_samps = np.random.normal(size=N_samps)
        
        # Record last values of mu and V to assess convergence.
        last_mu = np.copy(self.mu)
        last_V = np.copy(self.V)
        
        sw = 0
        dmu = np.Inf
        dV = np.Inf
        while max(dmu, dV) > tol:
            sw += 1
            # Update
            self.update_sweep()
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