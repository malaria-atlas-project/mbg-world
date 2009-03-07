# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import os, sys
import pymc as pm
import numpy as np
import tables as tb
from pylab import *
from EP import EP
from scipy import interpolate as interp
import mbgw
from mbgw import correction_factors
from mbgw.correction_factors import two_ten_factors, known_age_corr_likelihoods_f, stochastic_known_age_corr_likelihoods, age_corr_factors, S_trace, known_age_corr_factors
from mbgw.agepr import a
from IPython.Debugger import Pdb
import time
import matplotlib
matplotlib.interactive(True)

__all__ = ['visualize', 'ages_and_data', 'observed_gp_params', 'pred_samps']

# Use f or f + epsilon as information?
obs_on_f = False

def norm_dens(x,mu,V):
    out = np.exp(-(x-mu)**2/2./V)
    out /= out.max()
    return out

def visualize(M_pri, C_pri, E, lps, nug):
    import matplotlib
    matplotlib.interactive(True)
    k=2
    x = linspace(E.M[k] - 4.*np.sqrt(E.C[k,k]), E.M[k] + 4.*np.sqrt(E.C[k,k]), 100)
    
    pri_real = norm_dens(x, M_pri[k], C_pri[k,k])
    norms = np.random.normal(size=10000)*np.sqrt(nug[k])
    def this_lp(x, k=k, norms=norms):
        return np.array([pm.flib.logsum(lps[k](xi + norms)) - np.log((len(norms))) for xi in x])
    like_real = this_lp(x)
    like_real = np.exp(like_real - like_real.max())
    like_approx = norm_dens(x, E.mu[k], E.V[k] + nug[k])
    post_real = like_real * pri_real

    smoove_mat = np.asarray(pm.gp.cov_funs.gaussian.euclidean(x,x,amp=1,scale=.1))
    smoove_mat /= np.sum(smoove_mat, axis=0)
    like_real = np.dot(smoove_mat, like_real)
    post_real = np.dot(smoove_mat, post_real)
    post_real /= post_real.max()
    like_real /= like_real.max()       

    post_approx = norm_dens(x, E.M[k], E.C[k,k])
    
    # figure(1, figsize=(9,6))
    clf()
    plot(x, pri_real, 'g:', linewidth=2, label='Prior')
    plot(x, like_real, 'b-.', linewidth=2, label='Likelihood')
    plot(x, like_approx, 'r-.', linewidth=2, label='Approx. likelihood')
    plot(x, post_real, 'b-', linewidth=2, label='Posterior')
    plot(x, post_approx, 'r-', linewidth=2, label='Approx. posterior')
    legend(loc=0).legendPatch.set_alpha(0.)
    xlabel(r'$f(x)$')
    axis('tight')

    # Pdb(color_scheme='Linux').set_trace()   
    # savefig('figs/post%i.pdf'%k, transparent=True)    


def ages_and_data(N_exam, N_age_samps, f_samp, i, correction_factor_array, age_lims):
    """
    Called by pred_samps. Simulates data given f.
    """
    
    N_samp = len(f_samp)
    
    # Get a sample for the age distribution in the 2-10 age class
    age_distribution = {}
    for l in set(age_lims):
        age_distribution[l] = S_trace[np.random.randint(S_trace.shape[0]),0,l[0]:l[1]+1]
        age_distribution[l] /= np.sum(age_distribution[l])
    
    # Draw age for each individual, draw an age-correction profile for each location,
    # compute probability of positive for each individual, see how many individuals are
    # positive.
    A = []
    pos = []
    for s in xrange(N_samp):
        A.append(np.array(pm.rcategorical(age_distribution[age_lims[s]], size=N_exam[s]),dtype=int) + age_lims[s][0])
        P_samp = pm.invlogit(f_samp[s].ravel())*correction_factor_array[:,np.random.randint(N_age_samps)][A[-1]]
        pos.append(pm.rbernoulli(P_samp))
    
    return A, pos, age_distribution

def find_closest(x, y, N_nearest):
    big_D = np.empty((x.shape[0], y.shape[0]))
    pm.gp.geo_rad(big_D, x[:,:2], y[:,:2])
    closest = []
    for r in big_D:
        closest = np.hstack((closest, argsort(r)[:N_nearest]))
    closest = np.array(list(set(closest)),dtype=int)
    return closest
    
def observed_gp_params(combo_mesh, tracefile, trace_thin, trace_burn, N_nearest):
    """
    Called by pred_samps.  Generates thinned means and covariances, observed at relevant
    data locations, from an MCMC trace.
    """
    
    C = tracefile.root.chain0.group0.C
    M = tracefile.root.chain0.group0.M
    trace = tracefile.root.chain0.PyMCsamples.cols
    logp_mesh = tracefile.root.metadata.logp_mesh[:]
    data_mesh = tracefile.root.metadata.data_mesh[:]
    
    # Find the nearest N_combo data locations to each new sample location and uniquify.
    if obs_on_f:
        closest = find_closest(combo_mesh, logp_mesh, N_nearest)
    else:
        closest = find_closest(combo_mesh, data_mesh, N_nearest)
    
    # Pull means and covariances out of the trace and observe them on either f or f + epsilon,
    # and return them together with their corresponding nuggets.
    Ms = []
    Cs = []
    Vs = []
    for i in xrange(trace_burn,len(trace),trace_thin):
        # Mean and covariance conditional on existing data
        # TODO: COVARIATES!!!!
        M_here = M[i]        
        C_here = C[i]

        if obs_on_f:
            pm.gp.observe(M_here, C_here, obs_mesh=logp_mesh[closest], obs_vals=trace.f[i][closest])                        
        else:
            pm.gp.observe(M_here, C_here, obs_mesh=data_mesh[closest], obs_vals=trace.eps_p_f[i][closest], obs_V = trace.V[i])                                
            
        Ms.append(M_here)
        Cs.append(C_here)
        Vs.append(trace.V[i])
        
    return Ms, Cs, Vs

# EP_MAP.pred_samps(pred_mesh, samp_mesh, din.n, tracefile, trace_thin, N_outer, N_inner, N_nearest, din.lo_age, din.up_age, correction_factor_arrays)    
def pred_samps(pred_mesh, samp_mesh, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array):
    
    N_age_samps = correction_factor_array.shape[1]
    # print 1, time.time() - t
    obs_on_f = True
    
    combo_mesh = np.vstack((samp_mesh, pred_mesh))
    N_pred = pred_mesh.shape[0]
    N_samp = samp_mesh.shape[0]
    
    model_posteriors = []
    positive_observations = []
    likelihood_means = []
    likelihood_variances = []
    
    # Get mean, covariance and nugget from the MCMC trace.
    # Mean and covariance will have been observed at the most relevant
    # (nearest to prediction and/or sampling locations) data locations.
    # t = time.time()
    Ms, Cs, Vs = observed_gp_params(combo_mesh, tracefile, trace_thin, trace_burn, N_nearest)
    # print 2, time.time() - t
    
    N_param_vals = min(N_param_vals, len(Vs))
    N_per_param = min(N_per_param, len(Vs))
    
    ind_outer = np.array(np.linspace(0,len(Vs)-1,N_param_vals),dtype=int)
    ind_inner = np.array(np.linspace(0,len(Vs)-1,N_per_param),dtype=int)
    
    # Parameter values
    for ii in xrange(N_param_vals):
        print 'Parameter sample %i of %i'%(ii, N_param_vals)
        t = time.time()
        i = ind_outer[ii]        

        model_posteriors.append([])
        likelihood_means.append([])
        likelihood_variances.append([])

        M_pri = Ms[i](samp_mesh)
        C_pri = Cs[i](samp_mesh, samp_mesh)

        # Draw P' from prior.
        f_samp = pm.rmv_normal_cov(M_pri, C_pri + eye(N_samp)*Vs[i])

        # Get ages, number positive, and normalized age distribution for prediction
        ages, positives, age_distribution = ages_and_data(N_exam, N_age_samps, f_samp, i, correction_factor_array, age_lims)
        positive_observations.append(positives)
        
        # Make log-likelihood functions

        # FIXME: Uncomment!
        # marginal_log_likelihoods = known_age_corr_likelihoods_f(positives, ages, correction_factor_array, linspace(-10,10,100), 0)
        marginal_log_likelihoods = [lambda x: np.zeros(np.shape(x))]*len(positives)
        # from IPython.Debugger import Pdb
        # Pdb(color_scheme='Linux').set_trace()


        # Update posterior given simulated data
        this_nug = np.empty(N_samp)
        for jj in xrange(N_per_param):
            j = ind_inner[jj]
            
            M_pri = Ms[j](samp_mesh)
            C_pri = Cs[j](samp_mesh, samp_mesh)
            
            # Fit for a posterior of f + epsilon
            # E = EP(M_pri, C_pri, marginal_log_likelihoods, nug=np.hstack((np.zeros(N_samp), np.ones(N_pred)*np.Inf)))    
            this_nug.fill(Vs[j])
            E = EP(M_pri, C_pri, marginal_log_likelihoods, nug=this_nug)      
            try:          
                this_P = E.fit(10000)
                model_posteriors[-1].append(this_P)
            except:
                cls, inst, tb = sys.exc_info()
                print 'Error: %s'%inst.message
                model_posteriors[-1].append(-np.inf)
            likelihood_means[-1].append(E.mu)
            likelihood_variances[-1].append(E.V)
            # if jj%100==0:
            #     print '\t Sample %i of %i'%(jj,N_per_param)
            #     visualize(M_pri, C_pri, E, marginal_log_likelihoods, this_nug)            
            
            # predictive_samples[-1].append(draw_prediction(M_pred, sig_pred, N_pred, age_distribution, ages, N_age_samps))
        print time.time() - t 
    return ind_outer, ind_inner, Ms, Cs, Vs, np.asarray(likelihood_means), np.asarray(likelihood_variances), np.asarray(model_posteriors)