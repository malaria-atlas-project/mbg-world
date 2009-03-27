# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as pl
import numpy as np
import os
import pickle
import EP_MAP
from IPython.Debugger import Pdb
import pymc as pm
from copy import copy, deepcopy
from tables import openFile, FloatAtom
import mbgw
from mbgw.correction_factors import two_ten_factors, known_age_corr_likelihoods_f, stochastic_known_age_corr_likelihoods, age_corr_factors, S_trace, known_age_corr_factors
from mbgw.age_pr_datasets import a
from mbgw.master_grid import *

__all__ = ['frontend', 'backend', 'update_posterior', 'scratch_cleanup']

def frontend(fun):
	def new_fun(*args):
		return pickle.loads(fun(pickle.dumps(args)))
	return new_fun

def backend(fun):
	def new_fun(*args):
	    if len(args)==0:
	        return pickle.dumps(fun())
	    elif len(args)==1:
    		return pickle.dumps(fun(*pickle.loads(args[0])))
	return new_fun

def age_to_bin(age, a=a):
    """Finds the bin in a to which 'age' belongs."""
    return np.where(a>=age)[0][0]
    
def regularize_ages(lo_age, up_age):
    """
    Takes actual, year ages and maps them into the bins used in the age-correction model.
    This is needed because the age-correction model divides ages into 5-year bins after 15.
    """
    lo_age = np.array(map(age_to_bin, lo_age))
    up_age = np.array(map(age_to_bin, up_age))
    if np.any(up_age<lo_age):
        raise ValueError
    lo_age[np.where(up_age == lo_age)] -= 1
    return lo_age, up_age

def PR_samps(mesh, Ms, Cs, Vs, ind, facs):
    """
    Converts a mean function, covariance function, nugget and array of correction factors
    to a sample for the average of parasite rate over a given spatiotemporal mesh.
    """
    nm = mesh.shape[0]        
    samps = np.empty((len(ind), nm))
    for i in ind:
        C = Cs[i](mesh, mesh)
        C[::nm+1] += Vs[i]
        samps[i,:] = pm.invlogit(pm.mv_normal_cov(Ms[i](mesh), C).ravel()) * facs[A[i]]

    return np.mean(samps,axis=1)
    
def make_pt_fig(cur_val, samps, output_fname, output_path, outfigs_transparent=False, hist_color=(0,.1,.5), line_color='r-'):
    """Creates a png file from a point, writes it to disk and returns the path."""
    output_fname += '.png'
    pl.figure()
    h,b,p=pl.hist(samps,10,normed=True,facecolor=hist_color,histtype='stepfilled')
    pl.xlabel(r'$x$')
    pl.ylabel(r'$p(x)$')
    pl.plot([cur_val, cur_val],[0,h.max()],line_color,linewidth=2,label='Current value')
    pl.title('Utility function: standard deviation')
    l=pl.legend(loc=0)
    l.legendPatch.set_alpha(0)
    pl.savefig(output_path+'/'+output_fname, transparent=outfigs_transparent)
    return '/'.join([os.getcwd(), output_path, output_fname])

@backend
def scratch_cleanup():
    """Cleans out the image cache."""
    for file in os.listdir('web_scratch'):
        os.remove('web_scratch/'+file)

def make_EP_inputs(d):
    """
    Converts an input record array containing lat, lon, year and month columns to
    a spatiotemporal mesh that the gp package can handle.
    """
    samp_mesh = np.vstack((d.lon*deg_to_rad, d.lat*deg_to_rad, d.year + d.month/12. - 2009)).T
    return samp_mesh

def add_times(pred_mesh, nmonths):
    """Expands the prediction mesh to the requested number of months."""
    new_pred_mesh = np.empty((0, 3))
    for i in xrange(pred_mesh.shape[0]):
        p = np.repeat(np.atleast_2d(pred_mesh[i]), nmonths[i], axis=0)
        p[:,2] += np.arange(nmonths[i])/12.
        new_pred_mesh = np.vstack((new_pred_mesh, p))
    return new_pred_mesh
    
def one_point_mean(res, nmonths, i):
    """
    Takes the slice of 'res' that corresponds to a temporal slice at one spatial point,
    and computes the mean, and returns it.
    """
    start = np.sum(nmonths[:i])
    stop = np.sum(nmonths[:i+1])
    return np.mean(res[start:stop])
    
def robust_observe_and_eval(lm, M, C, samp_mesh, nug, pred_mesh):
    """
    Evaluates and observes a mean and covariance, and returns their evaluation on a
    prediction mesh. Works even for negative 'observation variance'.
    """
    nsm = samp_mesh.shape[0]
    if lm is not None:
        # Observe according to EP likelihood
        try:
            pm.gp.observe(M,C,samp_mesh,lm,nug)
            C_mesh = C(pred_mesh, pred_mesh)
            M_mesh = M(pred_mesh)
        except np.linalg.LinAlgError:
            # If there's a problem with the Cholesky-based algorithm (and there easily might be
            # because lv may be large and negative) then do it the boneheaded way.
            C_mesh = C(pred_mesh, pred_mesh)
            M_mesh = M(pred_mesh)
            C_samp = C(samp_mesh, samp_mesh).copy()
            M_samp = M(samp_mesh)

            np.ravel(C_samp)[::nsm+1] += nug # Add lv+V to diagonal of C_samp in place
            C_samp_I = C_samp.I

            C_off = C(samp_mesh, pred_mesh)
            C_mesh -= C_off.T * C_samp_I * C_off
            M_mesh += np.ravel(np.dot(C_off.T * C_samp_I , (lm - M_samp)))

    else:
        C_mesh = C(pred_mesh, pred_mesh)
        M_mesh = M(pred_mesh)
    return M_mesh, C_mesh

def simulate_on_pred(nsamp, Vp, M_mesh, nmonths, fac_array, lo_age, up_age):
    """Simulates data on the prediction mesh over the requested time period."""

    outp = (np.random.normal(size=(nsamp,len(Vp)))*np.sqrt(Vp) + M_mesh)
    outp = pm.flib.invlogit(outp.ravel()).reshape((nsamp,len(Vp))).T

    # For all prediction points,
    for i in xrange(len(nmonths)):
        these_facs = np.empty((nmonths[i],nsamp))
        # For each month at the prediction point, draw a new set of nsamp age-correction factors.
        for month in xrange(nmonths[i]):
            this_fac_array = fac_array[lo_age[i]:up_age[i], np.random.randint(fac_array.shape[1])]
            these_facs[month,:] = this_fac_array[np.random.randint(this_fac_array.shape[0], size=nsamp)]
        outp[np.sum(nmonths[:i]):np.sum(nmonths[:i+1]),:] *= these_facs

    if np.any(np.isnan(outp)):
        raise ValueError, 'NaN in results!'
    return outp

def make_meanifying_matrix(nmonths, Vp):
    """Returns a matrix that takes the required temporal means."""
    meanifying_mat = np.zeros((len(nmonths), len(Vp)))
    for i in xrange(len(nmonths)):
        start = np.sum(nmonths[:i])
        stop = np.sum(nmonths[:i+1])
        meanifying_mat[i,start:stop] = 1./(stop-start)
    return meanifying_mat

def make_justpix_samples(samp_mesh,pred_mesh,M,C,V,fac_array,lm,lv,nmonths,lo_age,up_age,nsamp=10000):
    """
    Converts the outputs of the EP algorithm into predictive samples at the output
    locations averaged over the requested time period.
    """
    npr = pred_mesh.shape[0]
    nsm = samp_mesh.shape[0]
    fac_array = fac_array
    
    if lv is not None:
        obs_nug = lv + V
    else:
        obs_nug = None
    M_mesh, C_mesh = robust_observe_and_eval(lm, M, C, samp_mesh, obs_nug, pred_mesh)
        
    # Sample at prediction points: do jointly
    Vp = np.diag(C_mesh) + V
    outp = simulate_on_pred(nsamp, Vp, M_mesh, nmonths, fac_array, lo_age, up_age)
    
    # Do requested temporal mean to each sample, location and return
    meanifying_mat = make_meanifying_matrix(nmonths, Vp)
    results = np.dot(meanifying_mat, outp)
    
    return results.T

def dtrm_disc_sample(p):
    """
    Like a categorical sample, generated using cdf inversion, where the input
    uniform samples are forced to be actually uniformly spaced over the interval.
    """
    n = len(p)
    pn = np.floor(np.cumsum(p)*n)
    samp = np.empty(n, dtype=int)
    samp[0] = pn[0]
    samp[1:] = np.diff(pn)
    return samp
        
def resample_with_mp(mp):
    """Takes log model posteriors and returns a 'sample' from them."""
    mp = mp - pm.flib.logsum(mp)
    mp = np.exp(mp)

    # Resample the slices of the posterior conditional on simulated dataset i.
    # indices=pm.rcategorical(mp, size=N_inner)
    index_samp = dtrm_disc_sample(mp)
    where_pos = np.where(index_samp>0)
    unique_indices = where_pos[0]
    n_copies = index_samp[where_pos]
    return unique_indices, n_copies
    # return output_info        

def ra_to_mesh(din, dout):
    """Converts the input and output record arrays into meshes."""
    samp_mesh = make_EP_inputs(din)
    pred_mesh = add_times(make_EP_inputs(dout), dout.nmonths)
    return pred_mesh, samp_mesh

def ra_to_age_info(din, dout):
    """Converts the input and output record arrays into age information."""
    lo_age_in, up_age_in = regularize_ages(din.lo_age, din.up_age) 
    lo_age_out, up_age_out = regularize_ages(dout.lo_age, dout.up_age)     
    age_lims = zip(lo_age_in, up_age_in)
    # FIXME: 1000 is a magic number
    correction_factor_array = known_age_corr_factors(np.arange(0,27), 1000)
    return age_lims, correction_factor_array, lo_age_out, up_age_out

def result_containers(utilities, N_outer, N_output, N_utilities):
    """Creates containers to hold the updated posterior."""
    # To hold samples from the predictive distribution of utility
    pred_samps = dict(zip([utility.__name__ for utility in utilities], 
                    [np.empty((N_outer, N_output)) for i in xrange(N_utilities)]))
    # To hold the current utility
    cur_vals =  dict(zip([utility.__name__ for utility in utilities], 
                    [np.empty(N_output) for i in xrange(N_utilities)]))
    # To hold current samples at prediction points.
    cur_samps = np.empty((0, N_output))
    return cur_samps, pred_samps, cur_vals

def create_output_info(N_output, utilities, cur_vals, pred_samps, outfile_path=None):
    """Creates the output_info object to return to the """
    if outfile_path is None:
        outfile_path = mbgw.__path__[0]+'/../testmbgw'
    output_info = []
    for i in xrange(N_output):
        output_info.append({})
        for utility in utilities:
            pl.close('all')
            output_info[i][utility.__name__]=make_pt_fig(cur_vals[utility.__name__][i], pred_samps[utility.__name__][:,i], str(id(output_info[i]))+'_'+utility.__name__, outfile_path+'/figs', outfigs_transparent=True, hist_color='.8', line_color='r-')
    return output_info

def dicts2rec(dicts):
    names = dicts[0].keys()
    arrs =  [np.array([d[key] for d in dicts]) for key in names]
    return np.rec.fromarrays(arrs, names=','.join(names))

# @backend
def update_posterior(input_pts, output_pts, tracefile, trace_thin, trace_burn, N_outer, N_inner, N_nearest, utilities=[np.std], outfile_path=None, nsamp_per_val=1000):
    """
    Inputs and outputs will be pickled as tuples.

    As input, expects two lists:
        input_pts: List of dictionaries of form:
            {'lon': float, 'lat': float, 'month': integer, 'year': integer, 'lo_age': integer, 'up_age': integer, 'n': integer}
        output_pts: List of dictionaries of form:
            {'lon': float, 'lat': float, 'month': integer, 'year': integer, 'nmonths': integer, 'lo_age': integer, 'up_age': integer}
        and one float, 'resolution', which gives the size of a pixel in decimal degrees.

        Defaults for both are:
            - lo_age: 2
            - up_age: 10
            - month: 6 (Jan=1, Dec=12)
            - year: 2009
            - lat, lon, n: no default, required from user.

        Longitudes and latitudes should be in decimal degrees.

    Returns:
        path: string. Path to new image patch.
        llc: (float, float) tuple. Lower left corner of image patch.
        urc: (float, float) tuple. Upper right corner of image patch.
        output_info: list of same length as output pts. Each element is of the form:
            (dict, string), where the dict is tabular information for display and string 
            is the path to a png image.
    """
    
    N_input, N_output = len(input_pts), len(output_pts)
    # Convert dictionaries to record arrays for easier handling
    din = dicts2rec(input_pts)
    dout = dicts2rec(output_pts)

    pred_mesh, samp_mesh = ra_to_mesh(din, dout)

    # Correction factors for each set of age limits.
    age_lims, correction_factor_array, lo_age_out, up_age_out = ra_to_age_info(din, dout)

    # Find posteriors with EP algorithm
    ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors = \
        EP_MAP.pred_samps(pred_mesh, samp_mesh, din.n, tracefile, trace_thin, trace_burn, N_outer, N_inner, N_nearest, age_lims, correction_factor_array)

    # Adjust for failures
    N_outer, N_inner, N_utilities = len(ind_outer), len(ind_inner), len(utilities)
    cur_samps, pred_samps, cur_vals = result_containers(utilities, N_outer, N_output, N_utilities)

    for i in xrange(N_outer):
        
        # Sample from predictive distribution at output points.
        ii = ind_outer[i]
        M, C, V = deepcopy(Ms[ii]), deepcopy(Cs[ii]), Vs[ii]
        new_samps = make_justpix_samples(samp_mesh, pred_mesh, M, C, V, correction_factor_array, None, None, dout.nmonths, lo_age_out, up_age_out, nsamp=nsamp_per_val)
        cur_samps = np.vstack((cur_samps,new_samps))
        
        unique_indices, n_copies = resample_with_mp(model_posteriors[i])
        print '%i indices of %i used.'%(len(unique_indices),N_inner)

        # Sample from predictive distribution at output points conditionally on simulated dataset i.        
        these_samps = np.empty((0,N_output))
        for ui in xrange(len(unique_indices)):
            j = unique_indices[ui]
            
            # Draw samples conditional on simulated dataset i and posterior slice j.
            lm, lv, jj = likelihood_means[i][j], likelihood_variances[i][j], ind_inner[j]
            M, C, V = deepcopy(Ms[jj]), deepcopy(Cs[jj]), Vs[jj]
            # Draw enough values to fill in all the slots that are set to j
            # in the importance resample.
            new_samps = make_justpix_samples(samp_mesh,pred_mesh,M,C,V,correction_factor_array,lm,lv,dout.nmonths,lo_age_out,up_age_out,nsamp=nsamp_per_val*n_copies[ui])
            these_samps = np.vstack((these_samps,new_samps))
        
        # Reduce predictive samples, conditional on simulated dataset i, with the utility functions.    
        for utility in utilities:
            pred_samps[utility.__name__][i,:] = np.apply_along_axis(utility, 0, these_samps)
    
    # Reduce predictive samples, unconditionally on any simulated data, with the utility function.
    for utility in utilities:
        cur_vals[utility.__name__] = np.apply_along_axis(utility, 0, cur_samps)
    
    output_info = create_output_info(N_output, utilities, cur_vals, pred_samps, outfile_path)

    # from IPython.Debugger import Pdb
    # Pdb(color_scheme='Linux').set_trace()