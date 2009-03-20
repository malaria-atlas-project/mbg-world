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
    """The normal density, without proportionality constant."""
    out = np.exp(-(x-mu)**2/2./V)
    out /= out.max()
    return out

def visualize(M_pri, C_pri, E, lps, nug):
    """
    Shows the EP algorithm's approximate posterior superimposed on the true posterior.
    Requires a fitted EP object as an input.
    """
    import matplotlib
    matplotlib.interactive(True)
    k=0
    
    x = linspace(E.M_pri[k] - 4.*np.sqrt(E.C_pri[k,k]), E.M_pri[k] + 4.*np.sqrt(E.C_pri[k,k]), 500)
    
    pri_real = norm_dens(x, M_pri[k], C_pri[k,k])
    norms = np.random.normal(size=10000)*np.sqrt(nug[k])
    
    def this_lp(x, k=k, norms=norms):
        return np.array([pm.flib.logsum(lps[k](xi + norms)) - np.log((len(norms))) for xi in x])
        
    like_real = this_lp(x)
    where_real_notnan = np.where(1-np.isnan(like_real))
    x_realplot = x[where_real_notnan]
    like_real = like_real[where_real_notnan]
    like_real = np.exp(like_real - like_real.max())
    like_approx = norm_dens(x, E.mu[k], E.V[k] + nug[k])
    post_real = like_real * pri_real[where_real_notnan]
    
    smoove_mat = np.asarray(pm.gp.cov_funs.gaussian.euclidean(x_realplot,x_realplot,amp=1,scale=.1))
    smoove_mat /= np.sum(smoove_mat, axis=0)
    # like_real = np.dot(smoove_mat, like_real)
    # post_real = np.dot(smoove_mat, post_real)
    post_real /= post_real.max() 
    
    post_approx = norm_dens(x, E.M[k], E.C[k,k])
    post_approx2 = pri_real * like_approx
    post_approx2 /= post_approx2.sum()
    
    post_approx /= post_approx.sum()
    post_real /= post_real.sum()
    like_real *= post_real.max()/like_real.max()      
    pri_real *= post_real.max()
    like_approx *= post_approx.max() / like_approx.max()
    
    # figure(1, figsize=(9,6))
    clf()
    plot(x, pri_real, 'g:', linewidth=2, label='Prior')
    plot(x_realplot, like_real, 'b-.', linewidth=2, label='Likelihood')
    plot(x, like_approx, 'r-.', linewidth=2, label='Approx. likelihood')
    plot(x_realplot, post_real, 'b-', linewidth=2, label='Posterior')
    plot(x, post_approx, 'r-', linewidth=2, label='Approx. posterior')
    plot(x, post_approx2, 'g-', linewidth=2, label='Approx. posterior meth 2')
    legend(loc=0).legendPatch.set_alpha(0.)
    xlabel(r'$f(x)$')
    axis('tight')
    m1r = sum(x[where_real_notnan]*post_real)/sum(post_real)
    m1em2 = sum(x*post_approx2)/sum(post_approx2)
    m1em = sum(x*post_approx)/sum(post_approx)
    m2r = sum(x[where_real_notnan]**2*post_real)/sum(post_real)
    m2em = sum(x**2*post_approx)/sum(post_approx)
    m2em2 = sum(x**2*post_approx2)/sum(post_approx2)
    print 'Posterior means: real: %s, EM: %s EM2: %s' % (m1r, m1em, m1em2)
    print 'Posterior variances: real: %s, EM: %s EM2: %s' % (m2r-m1r**2, m2em-m1em**2, m2em2-m1em2**2)
    # title('Prior variance %f, likelihood variance %f'%(C_pri[k,k], E.V[k]+nug[k]))
    
    # Pdb(color_scheme='Linux').set_trace()   
    # savefig('figs/post%i.pdf'%k, transparent=True)    


def ages_and_data(N_exam, f_samp, correction_factor_array, age_lims):
    """Called by pred_samps. Simulates ages of survey participants and data given f."""
    
    N_samp = len(f_samp)
    N_age_samps = correction_factor_array.shape[1]
    
    # Get samples for the age distribution at the observation points.
    age_distribution = []
    for i in xrange(N_samp):
        l = age_lims[i]
        age_distribution.append(S_trace[np.random.randint(S_trace.shape[0]),0,l[0]:l[1]+1])
        age_distribution[-1] /= np.sum(age_distribution[-1])
    
    # Draw age for each individual, draw an age-correction profile for each location,
    # compute probability of positive for each individual, see how many individuals are
    # positive.
    A = []
    pos = []
    for s in xrange(N_samp):
        A.append(np.array(pm.rcategorical(age_distribution[s], size=N_exam[s]),dtype=int) + age_lims[s][0])
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

def simulate_data(M_pri, C_pri, N_samp, V, N_exam, N_age_samps, correction_factor_array, age_lims):
    """Called by pred_samps in the outer loop to simulate data."""
    # Draw P' from prior.
    f_samp = pm.rmv_normal_cov(M_pri, C_pri + eye(N_samp)*V)

    # Get ages, number positive, and normalized age distribution for prediction
    ages, positives, age_distribution = ages_and_data(N_exam, f_samp, correction_factor_array, age_lims)
    
    sig = sqrt(diag(C_pri))
    lo = M_pri - sig*5
    hi = M_pri + sig*5

    # Make log-likelihood functions
    marginal_log_likelihoods = known_age_corr_likelihoods_f(positives, ages, correction_factor_array, linspace(lo.min(),hi.max(),500), 0)
    
    return marginal_log_likelihoods, positives

def fit_EP_to_sim_data(M_pri, C_pri, marginal_log_likelihoods, this_nug, debug=False):
    """
    Called by pred_samps in the inner loop to update the posterior given 
    simulated data.
    """
    E = EP(M_pri, C_pri, marginal_log_likelihoods, nug=this_nug)      
    try:          
        this_P = E.fit(10000, .0001)
    except:
        cls, inst, tb = sys.exc_info()
        print 'Error: %s'%inst.message
        this_P = -np.inf
    if debug:
        visualize(M_pri, C_pri, E, marginal_log_likelihoods, this_nug)
    return this_P, E.mu, E.V

def multi_append(tup, *lists):
    """Utility function for appending a new value to multiple lists."""
    [l.append(t) for l,t in zip(lists,tup)]

def pred_samps(pred_mesh, samp_mesh, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array, debug=False):
    """
    This function is meant to be called directly by the frontend. It returns 
    the prior means and covariances, as well as the products of the EP algorithm, 
    for a single simulated dataset.
    """
    N_age_samps = correction_factor_array.shape[1]
    obs_on_f = True
    N_pred, N_samp = pred_mesh.shape[0], samp_mesh.shape[0]
    model_posteriors, positive_observations, likelihood_means, likelihood_variances = [], [], [], []
    
    # Get mean, covariance and nugget from the MCMC trace.
    # Mean and covariance will have been observed at the most relevant
    # (nearest to prediction and/or sampling locations) data locations.
    Ms, Cs, Vs = observed_gp_params(np.vstack((samp_mesh, pred_mesh)), tracefile, trace_thin, trace_burn, N_nearest)

    N_param_vals, N_per_param = min(N_param_vals, len(Vs)), min(N_per_param, len(Vs))    
    ind_outer, ind_inner = np.array(np.linspace(0,len(Vs)-1,N_param_vals),dtype=int), np.array(np.linspace(0,len(Vs)-1,N_per_param),dtype=int)
    
    # Parameter values
    for ii in xrange(N_param_vals):
        print 'Parameter sample %i of %i'%(ii, N_param_vals)
        t = time.time()
        i = ind_outer[ii]        
        mps, lms, lvs = [], [], []
        
        marginal_log_likelihoods, pos = simulate_data(Ms[i](samp_mesh), Cs[i](samp_mesh, samp_mesh), N_samp, Vs[i], N_exam, N_age_samps, correction_factor_array, age_lims)

        positive_observations.append(pos)
        this_nug = np.empty(N_samp)
        for jj in xrange(N_per_param):
            # Update posterior given simulated data
            j = ind_inner[jj]
            M_pri, C_pri = Ms[j](samp_mesh), Cs[j](samp_mesh, samp_mesh)
            this_nug.fill(Vs[j])
            
            # Fit for a posterior of f + epsilon
            res = fit_EP_to_sim_data(M_pri, C_pri, marginal_log_likelihoods, this_nug, debug)
            multi_append(res, mps, lms, lvs)
        
        multi_append((mps, lms, lvs), model_posteriors, likelihood_means, likelihood_variances)    
        print time.time() - t             

    return ind_outer, ind_inner, Ms, Cs, Vs, np.asarray(likelihood_means), np.asarray(likelihood_variances), np.asarray(model_posteriors)