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
from copy import copy
from tables import openFile, FloatAtom
from mbgw.correction_factors import two_ten_factors, known_age_corr_likelihoods_f, stochastic_known_age_corr_likelihoods, age_corr_factors, S_trace, known_age_corr_factors
from mbgw.agepr import a

# __all__ = ['frontend', 'backend', 'PR_samps', 'make_pt_fig', 'make_img_patch', 'scratch_cleanup', 'make_EP_inputs', 'make_pred_meshes', 'make_samples', 'update_posterior']

# - Screen on big mac. Figure out why EP algorithm doesn't produce 
#     reasonable predictive distributions for std even when log-likelihood 
#     function always returns 0 (in this case, you should recover the 
#     current standard deviation exactly.) Debugging techniques:
#    - Rerun visualize in EP_map with the log-likelihood functions returning zero. 
#         You should get the posteriors very close to the prior.
#    - in update_posterior, set model posteriors to a constant. Then regardless of 
#         how badly the model posteriors in the EP are, if the mean and variance are 
#         correct, you should recover the correct predictive standard deviation.


rad_to_km = 6378.1/np.pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./np.pi
deg_to_rad = 1./rad_to_deg

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
    return np.where(a>=age)[0][0]
    
def regularize_ages(lo_age, up_age):
    lo_age= np.array(map(age_to_bin, lo_age))
    up_age = np.array(map(age_to_bin, up_age))
    if np.any(up_age<lo_age):
        raise ValueError
    lo_age[np.where(up_age == lo_age)] -= 1
    return lo_age, up_age

def PR_samps(mesh, Ms, Cs, Vs, ind, facs):
    nm = mesh.shape[0]        
    samps = np.empty((len(ind), nm))
    for i in ind:
        C = Cs[i](mesh, mesh)
        C[::nm+1] += Vs[i]
        samps[i,:] = pm.invlogit(pm.mv_normal_cov(Ms[i](mesh), C).ravel()) * facs[A[i]]

    return np.mean(samps,axis=1)
    
def make_pt_fig(pt, cur_val, samps, output_fname, output_path, outfigs_transparent=False, hist_color=(0,.1,.5), line_color='r-'):
    """
    Creates a png file from a point, writes it to disk and returns the path.
    """
    output_fname += '.png'
    # pl.close('all')
    pl.figure()
    h,b,p=pl.hist(samps,10,normed=True,facecolor=hist_color,histtype='stepfilled')
    pl.xlabel(r'$x$')
    pl.ylabel(r'$p(x)$')
    special_x = np.random.normal()
    print 'Current value: ',cur_val
    pl.plot([cur_val, cur_val],[0,h.max()],line_color,linewidth=2,label='Current value')
    pl.title('Utility function: standard deviation')
    l=pl.legend(loc=0)
    # pl.axis('tight')
    l.legendPatch.set_alpha(0)
    pl.savefig(output_path+'/'+output_fname, transparent=outfigs_transparent)
    return '/'.join([os.getcwd(), output_path, output_fname])

@backend
def scratch_cleanup():
    for file in os.listdir('web_scratch'):
        os.remove('web_scratch/'+file)

def make_EP_inputs(din):
    samp_mesh = np.vstack((din.lon, din.lat, din.year + (din.month-1)/12. - 2009)).T
    return samp_mesh

def add_times(pred_mesh, nmonths):
    """
    Expands the prediction mesh to the requested number of months.
    """
    new_pred_mesh = np.empty((0, 3))
    for i in xrange(pred_mesh.shape[0]):
        p = np.repeat(np.atleast_2d(pred_mesh[i]), nmonths[i], axis=0)
        p[:,2] += np.arange(nmonths[i])
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
        

def make_justpix_samples(samp_mesh,pred_mesh,M,C,V,fac_array,lm,lv,nmonths,lo_age,up_age,nsamp=1000):
    
    npr = pred_mesh.shape[0]
    fac_array = fac_array
    
    if lm is not None:
        # Observe according to EP likelihood
        try:
            pm.gp.observe(M,C,samp_mesh,lm,lv+V)
        except np.linalg.LinAlgError:
            C_old = C
            C = pm.gp.NearlyFullRankCovariance(C_old.eval_fun, **C_old.params)
            C.observe(C_old.obs_mesh, C_old.obs_V)
            pm.gp.observe(M,C,samp_mesh,lm,lv+V)            

    # Sample at prediction points: do jointly
    Vp = np.diag(C(pred_mesh,pred_mesh))
    Vp += V
    outp = (np.random.normal(size=(nsamp,len(Vp)))*np.sqrt(Vp) + M(pred_mesh))
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
    
    # Compress to annual means and return
    results = []
    for i in xrange(len(nmonths)):
        results.append(one_point_mean(outp, nmonths, i))
        
    return results


# @backend
def update_posterior(input_pts, output_pts, tracefile, trace_thin, trace_burn, N_outer, N_inner, N_nearest, utilities=[np.std]):
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
    
    N_input = len(input_pts)
    N_output = len(output_pts)
        
    # Convert dictionaries to record arrays for easier handling
    din = np.rec.fromrecords([input_pt.values() for input_pt in input_pts], names=input_pts[0].keys())
    dout = np.rec.fromrecords([output_pt.values() for output_pt in output_pts], names=output_pts[0].keys())

    for pt in (din,dout):
        for attr in ('lon','lat'):
            pt[attr] *= deg_to_rad
    
    samp_mesh = make_EP_inputs(din)
    pred_mesh = add_times(make_EP_inputs(dout), dout.nmonths)
        
    # Correction factors for each set of age limits.
    lo_age_in, up_age_in = regularize_ages(din.lo_age, din.up_age) 
    lo_age_out, up_age_out = regularize_ages(dout.lo_age, dout.up_age)     
    
    age_lims = zip(lo_age_in, up_age_in)
    correction_factor_array = known_age_corr_factors(np.arange(0,27), 1000)

    # Find posteriors with EP algorithm
    ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors = \
        EP_MAP.pred_samps(pred_mesh, samp_mesh, din.n, tracefile, trace_thin, trace_burn, N_outer, N_inner, N_nearest, age_lims, correction_factor_array)
    
    # Adjust for failures
    N_outer = len(ind_outer)
    N_inner = len(ind_inner)
    
    # Generate current and predictive utilities
    N_utilities = len(utilities)
    # To hold samples from the predictive distribution of utility
    pred_samps = dict(zip([utility.__name__ for utility in utilities], 
                    [np.empty((N_outer, N_output)) for i in xrange(N_utilities)]))
    # To hold the current utility
    cur_vals =  dict(zip([utility.__name__ for utility in utilities], 
                    [np.empty(N_output) for i in xrange(N_utilities)]))
    # To hold current samples at prediction points.
    cur_samps = np.empty((0, N_output))
    
    # Pdb(color_scheme='Linux').set_trace()   
    for i in xrange(N_outer):
        
        # print 'Point predictions: outer %i'%i
        
        # Sample from predictive distribution at output points.
        ii = ind_outer[i]
        M, C, V = copy(Ms[ii]), copy(Cs[ii]), Vs[ii]
        cur_samps = np.vstack((cur_samps,make_justpix_samples(samp_mesh, pred_mesh, M, C, V, correction_factor_array, None, None, dout.nmonths, lo_age_out, up_age_out, nsamp=10)))
        
        mp = model_posteriors[i]
        mp -= pm.flib.logsum(mp)
        mp = np.exp(mp)
        
        # Resample the slices of the posterior conditional on simulated dataset i.
        indices=pm.rcategorical(mp, size=N_inner)
        unique_indices = set(indices)
        print '%i indices of %i used.'%(len(unique_indices),len(indices))

        # Sample from predictive distribution at output points conditionally on simulated dataset i.        
        these_samps = np.empty((0,N_output))        
        for j in unique_indices:
            
            # Draw samples conditional on simulated dataset i and posterior slice j.
            lm, lv = likelihood_means[i][j], likelihood_variances[i][j]
            jj = ind_inner[j]
            M, C, V = copy(Ms[jj]), copy(Cs[jj]), Vs[jj]
            # Draw enough values to fill in all the slots that are set to j
            # in the importance resample.
            n_copies = 1 + np.sum(indices==jj)
            these_samps = np.vstack((these_samps,make_justpix_samples(samp_mesh,pred_mesh,M,C,V,correction_factor_array,lm,lv,dout.nmonths,lo_age_out,up_age_out,nsamp=10*n_copies)))
        
        # Reduce predictive samples, conditional on simulated dataset i, with the utility functions.    
        for utility in utilities:
            pred_samps[utility.__name__][i,:] = np.apply_along_axis(utility, 0, these_samps)
    
    # Reduce predictive samples, unconditionally on any simulated data, with the utility function.
    for utility in utilities:
        cur_vals[utility.__name__] = np.apply_along_axis(utility, 0, cur_samps)
    
    output_info = []
    for i in xrange(N_output):
        output_info.append({})
        for utility in utilities:
            pl.close('all')
            output_info[i][utility.__name__]=make_pt_fig(pt, cur_vals[utility.__name__][i], pred_samps[utility.__name__][:,i],str(id(output_info[i]))+'_'+utility.__name__, 'figs', outfigs_transparent=True, hist_color='.8', line_color='r-')

    from IPython.Debugger import Pdb
    Pdb(color_scheme='Linux').set_trace()   
    # return output_info
    
    # Create output_info
    # dum_objs = []
    # output_info = []
    # pl.figure()
    # for i in xrange(N_output):
    #     dum_objs.append([])
    #     pt = dout[i]
    #     pt_info = {'lon': pt.lon, 'lat': pt.lat, 'year': pt.year, 'lower age': pt.lo_age, 'upper age': pt.up_age, 'random number': str(np.random.random())}
    #     this_out_tup = (pt_info, make_pt_fig(pt, cur_utilities[i], pred_utility_samps[:,i], str(id(dum_objs[-1]))))
    #     output_info.append(this_out_tup)
    # 
    # return (path,(img_lon.min(), img_lat.min()), (img_lon.max(), img_lat.max()), output_info)