# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import matplotlib
matplotlib.interactive(False)
import matplotlib.pyplot as pl
import numpy as np
import os
import pickle
import EP_MAP
from IPython.Debugger import Pdb
import pymc as pm
from copy import copy
from tables import openFile, FloatAtom

__all__ = ['frontend', 'backend', 'PR_samps', 'make_pt_fig', 'make_img_patch', 'scratch_cleanup', 'make_EP_inputs', 'make_pred_meshes', 'make_samples', 'update_posterior']


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

def PR_samps(mesh, Ms, Cs, Vs, ind, facs):
    nm = mesh.shape[0]        
    samps = np.empty((len(ind), nm))
    for i in ind:
        C = Cs[i](mesh, mesh)
        C[::nm+1] += Vs[i]
        samps[i,:] = pm.invlogit(pm.mv_normal_cov(Ms[i](mesh), C)) * facs[A[i]]

    return np.mean(samps,axis=1)
    
def make_pt_fig(pt, cur_val, samps, output_fname, output_path, outfigs_transparent=False, hist_color=(0,.1,.5), line_color='r-'):
    """
    Creates a png file from a point, writes it to disk and returns the path.
    """
    output_fname += '.png'
    pl.close('all')
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
    return '/'.join(os.getcwd(), output_path, output_fname)

def make_img_patch(lon, lat, exp_surf, path, alpha=.8, cmap=matplotlib.cm.hot):
    pl.close('all')
    pl.figure()

    pl.imshow(exp_surf, extent=[lon.min(), lon.max(), lat.min(), lat.max()], origin='lower', alpha=alpha, cmap=cmap)
    
    output_fname = str(id(exp_surf))+'.png'
    pl.axis('off')
    pl.savefig(path+output_fname, transparent=True)
    return os.getcwd() + path + output_fname    

@backend
def scratch_cleanup():
    for file in os.listdir('web_scratch'):
        os.remove('web_scratch/'+file)

def make_EP_inputs(din):
    samp_mesh = np.vstack((din.lon, din.lat, din.year + (din.month-1)/12. - 2009)).T
    return samp_mesh

def make_pred_meshes(dout, img_yr):
    
    pred_mesh = np.repeat(np.vstack((dout.lon, dout.lat, dout.year-2009)),12,1).T
    pred_mesh[:,2] += np.tile(np.arange(12)/12.,len(dout))

    # TODO: return pred_mesh, img_mesh
    main_hdf = openFile(main_hdf_file)
    lon = main_hdf.root.long[:]*np.pi/180.
    lat = main_hdf.root.lat[:]*np.pi/180.
    main_hdf.close()
    
    return pred_mesh,lon,lat

# FIXME: This should use the on-disk realizations Just need to add nugget and compute utility. MUCH easier!    
def make_samples(samp_mesh,pred_mesh,img_mesh,img_slices,img_yr,img_samps,img_i,rt_now,M,C,V,r,lm=None,lv=None,nsamp=1):
    # TODO: take pred_mesh, img_mesh, return needful
    npr = pred_mesh.shape[0]
    
    if lm is not None:
        try:
            pm.gp.observe(M,C,samp_mesh,lm,lv+V)
        # In case of error:
        except np.linalg.LinAlgError:
            C_old = C
            C = pm.gp.NearlyFullRankCovariance(C_old.eval_fun, **C.params)
            C.observe(C_old.obs_mesh, C_old.obs_V)
            pm.gp.observe(M,C,samp_mesh,lm,lv+V)            
    
    # Sample at prediction points: do jointly
    sigp = C(pred_mesh,pred_mesh)
    sigp[::npr+1] += V
    sigp = np.linalg.cholesky(sigp)
    outp = pm.rmv_normal_chol(M(pred_mesh), sigp, size=nsamp)
    outp = (pm.flib.invlogit(outp) * r[np.random.randint(len(r), size=npr*nsamp)]).reshape((nsamp,npr))
    outp = np.mean(outp.reshape((nsamp,-1,12)),axis=2)
    
    # Samples at image points: do pointwise. Eventually do with a block-circulant call.
    nlon, nlat = len(img_mesh[0]), len(img_mesh[1])
    pix_mesh = np.empty((12,3))
    pix_mesh[:,2] = np.arange(12)/12. + img_yr - 2009
    for i in range(nlon)[img_slices[0]]:
        pix_mesh[:,0] = img_mesh[0][i]
        for j in range(nlat)[img_slices[1]]:
            pix_mesh[:,1] = img_mesh[1][j]
            sigp = C(pix_mesh, pix_mesh)
            sigp[::13] += V
            sigp = np.linalg.cholesky(sigp)
            these_samps = pm.rmv_normal_chol(M(pix_mesh), sigp, size=nsamp)
            these_samps = (pm.flib.invlogit(these_samps) * r[np.random.randint(len(r), size=nsamp*12)]).reshape((nsamp,12))
            these_samps = np.mean(these_samps, axis=1)
            img_samps[img_i,i,j] = these_samps[0]
            for k in xrange(nsamp-1):
                img_samps[rt_now[k],i,j] = these_samps[k+1]
                
    # Pdb(color_scheme='Linux').set_trace()   
    return outp

# TODO: Factor some of this shit out into functions!!
# TODO: Get it working from a terminal before you get Will involved!
# TODO: Don't even attempt to make image patch initially.    
@backend
def update_posterior(input_pts, output_pts, img_yr, utility=np.std):
    """
    Inputs and outputs will be pickled as tuples.

    As input, expects two lists:
        input_pts: List of dictionaries of form:
            {'lon': float, 'lat': float, 'month': integer, 'year': integer, 'lo_age': integer, 'up_age': integer, 'n': integer}
        output_pts: List of dictionaries of form:
            {'lon': float, 'lat': float, 'year': integer, 'lo_age': integer, 'up_age': integer}
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

    # FIXME: Check the units... they might be stored correctly already.
    # for pt in (din,dout):
    #     for attr in ('lon','lat'):
    #         pt[attr] *= np.pi/180.
    
    samp_mesh = make_EP_inputs(din)
    ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors = \
        EP_MAP.pred_samps(np.empty((0,3)), samp_mesh, N_exam)
        
    N_inner = len(ind_inner)
    N_outer = len(ind_outer)
    
    # Figure out output bounding box and number of pixels
    pred_mesh,  img_lon, img_lat = make_pred_meshes(dout, img_yr)
    nlon, nlat = len(img_lon), len(img_lat)

    # Age-correction factors
    age_distribution = EP_MAP.S_trace[np.random.randint(EP_MAP.S_trace.shape[0]),0,2:11]
    age_distribution /= np.sum(age_distribution)
    mean_facs=np.dot(age_distribution, EP_MAP.correction_factor_array)
    
    # Generate current and predictive utilities
    cur_samps = np.empty((N_outer, N_output))
    pred_utility_samps = np.empty((N_outer, N_output))
    pred_samps = np.empty((N_inner, N_output))
    
    # Current and predictive expected utilities
    img_scratch = openFile('web_scratch/'+str(id(img_lon))+'.hdf5','w')
    img_exp_utility=img_scratch.createCArray('/','img_exp_utility',FloatAtom(),(nlon,nlat))
    img_samps=img_scratch.createCArray('/','img_samps',FloatAtom(), (N_inner,nlon,nlat))

    # Initialize scratch space
    img_exp_utility[:] = 0

    # Pdb(color_scheme='Linux').set_trace()   
    for i in xrange(N_outer):
        print 'Point predictions: outer %i'%i
        ii = ind_outer[i]
        mp = model_posteriors[i]
        mp -= pm.flib.logsum(mp)
        mp = np.exp(mp)
        # importance-resample ind_outer and uniquify    
        indices=pm.rcategorical(mp, size=N_inner)
        nr, rf, rt, nu, xu, ui = pm.gp.linalg_utils.remove_duplicates(indices)
        
        for j in xrange(1,len(ui)):
            if ui[j]==0:
                break
        rf = rf[:j]
        rt = rt[:j]
        ui = ui[:j]
        
        # FIXME: This should target the appropriate year.
        cur_utilities = np.empty(N_output)
        for j in xrange(N_output):
            this_lon = np.argmin(np.abs(img_lon)-dout.lon[j])
            this_lat = np.argmin(np.abs(img_lat)-dout.lat[j])
            cur_utilities[i] = cur_img.root.data[this_lat, this_lon]
        
        # FIXME: M, C should be observed to make a correct image patch, not just to make
        # correct values at the output points.
        # Current samples
        M, C, V = Ms[ii], Cs[ii], Vs[ii]
        
        jc=0
        # Loop over indices that are represented in the importance resample.
        for j in ui:
            jc+=1
            if jc%100==0:
                print '\t inner %i of %i'%(jc,len(ui))
            
            # Load up mean and covariance parameters
            jj = ind_inner[indices[j]]
            M, C, V = copy(Ms[jj]), copy(Cs[jj]), Vs[jj]
            
            # Draw enough values to fill in all the slots that are set to j
            # in the importance resample.
            n_copies = 1 + np.sum(rf==j)
            rt_now = rt[np.where(rf==j)]
            # TODO: Split into two functions!
            these_samps = make_samples(samp_mesh,pred_mesh,img_mesh,img_slices,img_yr,img_samps,j,rt_now,M,C,V,mean_facs,lm=likelihood_means[i][indices[j]],lv=likelihood_variances[i][indices[j]], nsamp=n_copies)
            pred_samps[j,:] = these_samps[0,:]
            
            for k in xrange(len(rt_now)):
                pred_samps[rt_now[k]] = these_samps[k+1,:]
        
        # FIXME: Why the fuck is utility getting bound to 2009?
        utility = np.std                    
        pred_utility_samps[i,:] = np.apply_along_axis(utility, 0, pred_samps)
        # FIXME: Do this slice-by-slice, will be too big in the full image.
        img_exp_utility += np.apply_along_axis(utility, 0, img_samps)

    img_exp_utility /= float(N_outer)
    
    # Create new patch figure
    path = make_img_patch(img_lon, img_lat, img_exp_utility)

    # Create output_info
    # TODO: MAKE THIS REAL OUTPUT INFO
    dum_objs = []
    output_info = []
    pl.figure()
    for i in xrange(N_output):
        dum_objs.append([])
        pt = dout[i]
        pt_info = {'lon': pt.lon, 'lat': pt.lat, 'year': pt.year, 'lower age': pt.lo_age, 'upper age': pt.up_age, 'random number': str(np.random.random())}
        this_out_tup = (pt_info, make_pt_fig(pt, cur_utilities[i], pred_utility_samps[:,i], str(id(dum_objs[-1]))))
        output_info.append(this_out_tup)

    return (path,(img_lon.min(), img_lat.min()), (img_lon.max(), img_lat.max()), output_info)