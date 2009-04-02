# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

#Change
import pymc as pm
import numpy as np
import os
from copy import copy
from correction_factors import age_corr_likelihoods
from scipy import interpolate as interp
from st_cov_fun import *
import time
import auxiliary_data
# import MAPData
import gc
# from get_covariates import extract_environment_to_hdf5
from tables import ObjectAtom

__all__ = ['st_mean_comp', 'MyGibbsStepper', 'create_model']

def nearest_interp(lon_from, lat_from, data, lon_to, lat_to):
    out = np.empty(len(lon_to))
    for i in xrange(len(lon_to)):
        ilon = np.argmin(np.abs(lon_from - lon_to[i]))
        ilat = np.argmin(np.abs(lat_from - lat_to[i]))
        out[i] =  data[ilat, ilon]
    return out

def st_mean_comp(x, m_const, t_coef):
    lon = x[:,0]
    lat = x[:,1]
    t = x[:,2]
    return m_const + t_coef * t

class MyGibbsStepper(pm.StepMethod):
    """
    A special Gibbs stepper that updates tau and f together give eps_p_f
    and the mean and covariance parameters. Why did I update tau and f together?
    I don't remember. Probably grasping at straws trying to improve mixing
    at some point. This doesn't save any Cholesky decompositions or covariance 
    evaluations.
    """
    def __init__(self, f, tau, V, C_eval, M_eval, logp_mesh, eps_p_f, ti, incomp_jump=False, jump_tau = True):
        self.f = f
        self.tau = tau
        self.V = V
        self.C_eval = C_eval
        self.M_eval = M_eval
        self.logp_mesh = logp_mesh
        self.eps_p_f = eps_p_f
        self.incomp_jump = incomp_jump
        self.ti = ti
        
        self.scratch1 = np.asmatrix(np.empty(self.C_eval.value.shape, order='F'))
        self.scratch2 = np.asmatrix(np.empty(self.C_eval.value.shape, order='F'))
        self.scratch3 = np.empty(self.M_eval.value.shape)      

        self.jump_tau = jump_tau
        if self.jump_tau:
            pm.StepMethod.__init__(self,[f,tau])
        else:
            pm.StepMethod.__init__(self,[f])
            
        self.time = 0.
        self.incomp_time = 0.
    
    def step(self):
        
        fc = pm.gp.fast_matrix_copy
        
        t1 = time.time()
    
        eps_p_f = self.eps_p_f.value
        f = self.f.value
        for i in xrange(len(self.scratch3)):
            self.scratch3[i] = np.sum(eps_p_f[self.ti[i]] - f[i])
    
        if self.jump_tau:
            # tau's full conditional is Gamma.
            amp_alph = eps_p_f.shape[0]/2.+1.
            tau_alph = amp_alph + self.tau.parents['alpha']
            tau_vec = self.scratch3 - self.f.value
            tau_bet = np.dot(tau_vec, tau_vec)/2. + self.tau.parents['beta']
            self.tau.value = pm.rgamma(tau_alph,tau_bet)

        # Compute Cholesky factor of covariance of eps_p_f, C(x,x) + V
        in_chol = fc(self.C_eval.value, self.scratch1)
        for i in xrange(self.logp_mesh.shape[0]):
            in_chol[i,i] += self.V.value / len(self.ti[i])
        info = pm.gp.linalg_utils.dpotrf_wrap(in_chol)
        if info > 0:
            raise np.linalg.LinAlgError

        # Compute covariance of f conditional on eps_p_f.
        offdiag = fc(self.C_eval.value, self.scratch2)
        offdiag = pm.gp.trisolve(in_chol, offdiag, uplo='U', transa='T', inplace=True)

        C_step = offdiag.T * offdiag
        C_step *= -1
        C_step += self.C_eval.value
        
        # Compute mean of f conditional on eps_p_f.
        for i in xrange(len(self.scratch3)):
            self.scratch3[i] = np.mean(eps_p_f[self.ti[i]])
        m_step = self.M_eval.value + np.dot(offdiag.T, pm.gp.trisolve(in_chol,(self.scratch3 - self.M_eval.value),uplo='U',transa='T')).view(np.ndarray).ravel()
        
        sig_step = C_step
        info = pm.gp.linalg_utils.dpotrf_wrap(C_step.T)
        if info > 0:
            print 'WARNING: Cholesky decomposition of full conditional covariance not positive definite.'
            return
            raise np.linalg.LinAlgError
        
        # Update value of f.
        self.f.value = m_step+np.dot(sig_step,np.random.normal(size=sig_step.shape[1])).view(np.ndarray).ravel()
        
        t2 = time.time()
        self.time += t2-t1

region_trans = {'Africa':'AF','Americas':'AM','Asia':'AS'}

def create_model(region_name, all_pts, name, scale_params, amp_params, cpus,
    with_stukel, spatial, chunk, covariate_names, disttol, ttol, 
    AM_delay = 50000, AM_interval=100, AM_sd=.1, crashed_db = None):
    
    # ======================================
    # = Make sure it's safe to make output =
    # ======================================
    
    if not spatial:
        name += '_nonspatial'
    if with_stukel:
        name += '_stukel'
    for cname in covariate_names:
        name += '_%s'%cname
    
    if name + '.hdf5' in os.listdir('.'):
        print
        print """=============
= ATTENTION =
============="""
        print

        OK=False
        while not OK:
            y=raw_input('Database %s already exists.\nDo you want to delete it? Error will be raised otherwise.\n>> ' % (name+'.hdf5'))
            if y.lower() == 'yes':
                print 'OK, moving to trash.'
                os.system('mv %s ~/.Trash'%(name+'.hdf5'))
                OK=True
            elif y.lower() == 'no':
                raise RuntimeError, 'But dash it all! I mean to say, what?'
            else:
                y=raw_input('Please type yes or no.\n>> ')
            
        
    norun_name = '_'.join(name.split('_')[:2])
    
    C_time = [0.]
    f_time = [0.]
    M_time = [0.]
    
    # =============================
    # = Preprocess data, uniquify =
    # =============================
    
    # Convert latitude and longitude from degrees to radians.
    lon = all_pts.LONG*np.pi/180.
    lat = all_pts.LAT*np.pi/180.

    # Convert time to end year - 2009 (no sense forcing mu to adjust by too much).
    # t = all_pts.YEAR_START-2009. + all_pts.MONTH_STAR / 12.
    t = all_pts.TIME - 2009
    
    # Make lon, lat, t triples.
    data_mesh = np.vstack((lon, lat, t)).T   
    
    disttol=disttol/6378.
    ttol=ttol/12.

    # Find near spatiotemporal duplicates.
    if spatial:
        ui = []
        ri = []
        fi = []
        ti = []
        dx = np.empty(1)
        for i in xrange(data_mesh.shape[0]):
            match=False
            for j in xrange(len(ui)):
                pm.gp.geo_rad(dx, data_mesh[i,:2].reshape((1,2)), data_mesh[ui[j],:2].reshape((1,2)))
                dt = abs(t[ui[j]]-t[i])
                
                if dx[0]<disttol and dt<ttol:
                    match=True
                    fi.append(j)
                    ti[j].append(i)
                    ri.append(i)
                    break

            if not match:
                fi.append(len(ui))            
                ui.append(i)
                ti.append([i])
        ui=np.array(ui)
        ti = [np.array(tii) for tii in ti]
        fi = np.array(fi)
        ri = np.array(ri)        
        logp_mesh = data_mesh[ui,:]
        if len(ri)>0:
            repeat_mesh = data_mesh[ri,:]
        else:
            repeat_mesh = np.array([])
    else:
        ui = np.arange(len(t))
        ti = [np.array([uii]) for uii in ui]
        fi = ui
        ri = np.array([])
        logp_mesh = data_mesh
        repeat_mesh = np.array([])
        
    # =====================
    # = Create PyMC model =
    # =====================
    
    init_OK = False
    while not init_OK:
    
        # Flat prior on m_const (mu). 
        m_const = pm.Uninformative('m_const', value=-3.)
        if with_stukel:
            m_const.value = -1.1

        # Flat prior on coefficient of time (k).
        t_coef = pm.Uninformative('t_coef',value=.1)
        if with_stukel:
            t_coef.value = -.4

        # Inverse-gamma prior on nugget variance V.
        tau = pm.Gamma('tau', value=2., alpha=.001, beta=.001/.25)
        V = pm.Lambda('V', lambda tau=tau:1./tau)
    
        vars_to_writeout = ['V', 'm_const', 't_coef']
                
        # Pull out covariate information.
        # The values of covariate_dict are (Stochastic, interpolated covariate) tuples.
        # Interpolation is done to the data mesh.
        covariate_dict = {}
        for cname in covariate_names:
            # hf = openFile(mbgw.__path__[0] + '/auxiliary_data/' + cname + '.hdf5')
            if cname == 'periurb':
                this_interp_covariate = all_pts.URB_CLS==2
                if np.sum(all_pts.URB_CLS==3) < 10:
                    print 'Warning: Very few urban points, using same coefficient for urban and periurban'
                    this_interp_covariate += all_pts.URB_CLS==3
            elif cname == 'urb':
                if np.sum(all_pts.URB_CLS==3) >= 10:
                    this_interp_covariate = all_pts.URB_CLS==3
                else:
                    this_interp_covariate = None
            else:
                this_cov = getattr(auxiliary_data, cname)
                this_interp_covariate = nearest_interp(this_cov.long[:], this_cov.lat[:], this_cov.data, data_mesh[:,0], data_mesh[:,1])            
            if this_interp_covariate is not None:
                this_coef = pm.Uninformative(cname + '_coef', value=0.)
                covariate_dict[cname] = (this_coef, this_interp_covariate)
        
        # Lock down parameters of Stukel's link function to obtain standard logit.
        # These can be freed by removing 'observed' flags, but mixing gets much worse.
        if with_stukel:
            a1 = pm.Uninformative('a1',.5)
            a2 = pm.Uninformative('a2',.8)
        else:
            a1 = pm.Uninformative('a1',0,observed=True)
            a2 = pm.Uninformative('a2',0,observed=True)        

        transformed_spatial_vars=[V]
        if spatial:
            # Make it easier for inc (psi) to jump across 0: let nonmod_inc roam freely over the reals,
            # and mod it by pi to get the 'inc' parameter.
            nonmod_inc = pm.Uninformative('nonmod_inc', value=.5)
            inc = pm.Lambda('inc', lambda nonmod_inc = nonmod_inc: nonmod_inc % np.pi)

            # Use a uniform prior on sqrt ecc (sqrt ???). Using a uniform prior on ecc itself put too little
            # probability mass on appreciable levels of anisotropy.
            sqrt_ecc = pm.Uniform('sqrt_ecc', value=.1, lower=0., upper=1.)
            ecc = pm.Lambda('ecc', lambda s=sqrt_ecc: s**2)

            # Subjective skew-normal prior on amp (the partial sill, tau) in log-space.
            # Parameters are passed in in manual_MCMC_supervisor.
            log_amp = pm.SkewNormal('log_amp',**amp_params)
            amp = pm.Lambda('amp', lambda log_amp = log_amp: np.exp(log_amp))

            # Subjective skew-normal prior on scale (the range, phi_x) in log-space.
            log_scale = pm.SkewNormal('log_scale',**scale_params)
            scale = pm.Lambda('scale', lambda log_scale = log_scale: np.exp(log_scale))

            # Exponential prior on the temporal scale/range, phi_t. Standard one-over-x
            # doesn't work bc data aren't strong enough to prevent collapse to zero.
            scale_t = pm.Exponential('scale_t', .1)

            # Uniform prior on limiting correlation far in the future or past.
            t_lim_corr = pm.Uniform('t_lim_corr',0,1,value=.8)

            # # Uniform prior on sinusoidal fraction in temporal variogram
            sin_frac = pm.Uniform('sin_frac',0,1)
            
            vars_to_writeout.extend(['inc','ecc','amp','scale','scale_t','t_lim_corr','sin_frac'])
            transformed_spatial_vars.extend([inc,ecc,amp,scale])
    
        # Collect stochastic variables with observed=False for the adaptive Metropolis stepper.
        trial_stochs = [v[0] for v in covariate_dict.itervalues()] + [ m_const, tau, a1, a2, t_coef]        
        if spatial:
            trial_stochs = trial_stochs + [nonmod_inc, sqrt_ecc, log_amp, log_scale, scale_t, t_lim_corr, sin_frac]
        nondata_stochs = []
        for stoch in trial_stochs:
            if not stoch.observed:
                nondata_stochs.append(stoch)
                
        # Collect variables to write out

        # The mean of the field
        @pm.deterministic
        def M(m=m_const, tc=t_coef):
            return pm.gp.Mean(st_mean_comp, m_const = m, t_coef = tc)
        
        # The mean, evaluated  at the observation points, plus the covariates    
        @pm.deterministic(trace=False)
        def M_eval(M=M, lpm=logp_mesh, cv=covariate_dict):
            out = M(lpm)
            for c in cv.itervalues():
                out += c[0]*c[1][ui]
            return out

        # Create covariance and MV-normal F if model is spatial.   
        if spatial:
            try:
                # A constraint on the space-time covariance parameters that ensures temporal correlations are 
                # always between -1 and 1.
                @pm.potential
                def st_constraint(sd=.5, sf=sin_frac, tlc=t_lim_corr):    
                    if -sd >= 1./(-sf*(1-tlc)+tlc):
                        return -np.Inf
                    else:
                        return 0.

                # A Deterministic valued as a Covariance object. Uses covariance my_st, defined above. 
                @pm.deterministic
                def C(amp=amp,scale=scale,inc=inc,ecc=ecc,scale_t=scale_t, t_lim_corr=t_lim_corr, sin_frac=sin_frac):
                    return pm.gp.FullRankCovariance(my_st, amp=amp, scale=scale, inc=inc, ecc=ecc,st=scale_t, sd=.5,
                                                    tlc=t_lim_corr, sf = sin_frac, n_threads=cpus)

                # The evaluation of the Covariance object.
                @pm.deterministic(trace=False)
                def C_eval(C=C):
                    return C(logp_mesh, logp_mesh)

                # The field evaluated at the uniquified data locations
                f = pm.MvNormalCov('f',M_eval,C_eval,value=M_eval.value)

                # The field evaluated at all the data locations
                @pm.deterministic(trace=False)
                def f_eval(f=f):
                    return f[fi]

                init_OK = True
            except pm.ZeroProbability, msg:
                print 'Trying again: %s'%msg
                init_OK = False
                gc.collect()

        # if not spatial
        else:
            C=None

            # The field is just the mean, there's no spatially-structured component.
            @pm.deterministic
            def f(M=M_eval):
                return M[fi]
            f_eval = f
            
            init_OK=True            
    
    # ===========================
    # = Create likelihood layer =
    # ===========================
        
    eps_p_f_list = []
    N_pos_list = []
    
    # Obtain the spline representation of the log of the Monte Carlo-integrated 
    # likelihood function at each datapoint. The nodes are at .01,.02,...,.98,.99 .
    junk, splreps = age_corr_likelihoods(all_pts, 10000, np.arange(.01,1.,.01), norun_name)
    for i in xrange(len(splreps)):
        splreps[i] = list(splreps[i])

    # Don't worry, these are just reasonable initial values...
    val_now = pm.logit(np.array(all_pts.PF+1,dtype=float)/(all_pts.EXAMINED+2))
    if with_stukel:
        val_now = pm.stukel_logit(np.array(all_pts.PF+1,dtype=float)/(all_pts.EXAMINED+2), a1.value, a2.value)
    
    if data_mesh.shape[0] % chunk == 0:
        additional_index = 0
    else:
        additional_index = 1
    
    for i in xrange(0,data_mesh.shape[0] / chunk + additional_index):
        
        this_slice = slice(chunk*i, min((i+1)*chunk, data_mesh.shape[0]))

        # epsilon plus f, given f.
        @pm.stochastic(trace=False, dtype=np.float)
        def eps_p_f_now(value=val_now[this_slice], f=f_eval, V=V, this_slice = this_slice):
            return pm.normal_like(value, f[this_slice], 1./V)
        eps_p_f_now.__name__ = "eps_p_f%i"%i
        eps_p_f_list.append(eps_p_f_now)
        
        # The number positive: the data. Uses the spline interpolations of the likelihood
        # functions to compute them.
        try:
            @pm.data
            @pm.stochastic(dtype=np.int)
            def N_pos_now(value = pm.utils.round_array(all_pts.PF[this_slice]), splrep = splreps[this_slice], eps_p_f = eps_p_f_now, a1=a1, a2=a2):
                p_now = pm.flib.stukel_invlogit(eps_p_f, a1, a2)
                out = 0.
                for i in xrange(len(value)):
                    out += interp.splev(p_now[i], splrep[i])
                return out
        except ValueError:
            raise ValueError, 'Log-likelihood is nan at chunk %i'%i

    # Combine the eps_p_f values. This is stupid, I should have just used a Container.
    # I guess this makes it easier to keep traces.
    @pm.deterministic
    def eps_p_f(eps_p_f_list=eps_p_f_list):
        out = np.zeros(data_mesh.shape[0])
        for i in xrange(len(eps_p_f_list)):
            out[chunk*i:min((i+1)*chunk, data_mesh.shape[0])] = eps_p_f_list[i]
        return out

    # ==============================
    # = Create algorithmic objects =
    # ==============================

    # The MCMC object. Put the traces under maximal compression. NB can thin Asia by 50
    # if needed to reduce trace size on disk.
    if crashed_db is None:
        S = pm.MCMC([nondata_stochs,f,eps_p_f_list,tau,eps_p_f,C,M,transformed_spatial_vars],name=name,
                        db='hdf5', complevel=9, tune_interval=100)
    
        # Write metadata into hdf archive
        l=locals()
        hf = S.db._h5file
        hf.createGroup('/','metadata')
        weird_attrs = ['ti','vars_to_writeout','covariate_names','scale_params','amp_params']
        for meta_attr in ['logp_mesh', 'data_mesh','repeat_mesh','ui','fi','with_stukel','spatial','chunk','disttol','ttol','cpus']:
            this_attr = np.atleast_1d(l[meta_attr])
            if len(this_attr) > 0:
                hf.createArray(hf.root.metadata,meta_attr,this_attr)
            else:
                weird_attrs.append(meta_attr)
        for weird_attr in weird_attrs:
            vla=hf.createVLArray(hf.root.metadata,weird_attr,ObjectAtom())
            vla.append(l[weird_attr])
        for cname in covariate_names:
            if covariate_dict.has_key(cname):
                hf.createArray(hf.root.metadata,cname,covariate_dict[cname][1])    
    else:
        print 'Sampler crashed, attempting recovery from database + ' ,crashed_db
        S = pm.MCMC([nondata_stochs,f,eps_p_f_list,tau,eps_p_f,C,M,transformed_spatial_vars],name=name,
                        db=crashed_db, complevel=9, tune_interval=100)

    if spatial:
        S.use_step_method(MyGibbsStepper, f, tau, V, C_eval, M_eval, logp_mesh, eps_p_f, ti, jump_tau=False)
    S.use_step_method(pm.AdaptiveMetropolis, nondata_stochs, delay=AM_delay, interval=AM_interval)
    S.step_method_dict[m_const][0].proposal_sd *= AM_sd
        
    # This is completely unnecessary.
    for i in xrange(len(eps_p_f_list)):
        S.use_step_method(pm.Metropolis, eps_p_f_list[i])
        S.step_method_dict[eps_p_f_list[i]][0].proposal_sd *= .5
        
    S.assign_step_methods()
    
    out = locals()
    for v in covariate_dict.iteritems():
        out[v[0]] = v[1][0]
    return out

def write_out(fname, burn):
    from tables import openFile
    hf = openFile(fname)
    name = fname.split('/')[-1].replace('.hdf5','')
    hr = hf.root.chain0.PyMCsamples.cols

    data_dir_name = name+'_data'
    os.mkdir(data_dir_name)
    os.chdir(data_dir_name)

    trace_len = len(hr.V)
    trace_slice = np.linspace(burn, trace_len-1, 500).astype(int)

    for stoch_name in hr._v_colnames:
        if not stoch_name == 'f':
            dfile = file(stoch_name+'.txt','w')
            try:
                getattr(hr,stoch_name)[:][trace_slice].tofile(dfile,sep=', ')
            except:
                pass
            dfile.write('\n')
            dfile.close()

    dfile = file('f.txt','w')

    for i in trace_slice:
        hr.f[i].tofile(dfile,sep=', ')
        dfile.write('\n')
    dfile.close()
    name = fname

    os.system('cp ../r.txt .')

    readme_str = """
    Data for %s.

    scale.txt: Trace of the 'scale' parameter. Multiply distance matrix by 1/scale before
        overwriting with covariance matrix.
    amp.txt: Trace of the 'amp' parameter. Multiply covariance matrix by amp**2. 
    inc.txt: Trace of the 'inc' anisotropy parameter.
    ecc.txt: Trace of the 'ecc' anisotropy parameter.
    scale_t.txt: Trace of the 'scale_t' parameter. Like scale, but in the time direction.
    t_lim_corr.txt: Limiting correlation in the time direction.
    urb_coef.txt, periurb_coef.txt: Urban and periurban mean boosts
    sin_frac.txt: Relative size of the sinusoidal component in the temporal variogram.

    V.txt: Trace of the nugget variance.

    r.txt: Conversion factor from P' to PfPR210.

    m_const.txt: Trace of prior mean of transformed PfPR.
    t_coef.txt: Trace of the 't_coef' parameter. Coefficient of time in the mean function.

    f.txt: Trace of transformed PfPR, evaluated at the datapoints.
        Each row is one sample from the joint distribution.
        Each column is the trace for f's evaluation at a single datapoint.

    """ % (name)

    readme = file('README.txt','w')
    readme.write(readme_str)

    os.chdir('..')

    os.system('tar -cf %s_data.tar.gz %s --gzip' % (fname[:-5],data_dir_name))
    os.system('chmod ugo+rwx %s_data.tar.gz' % fname[:-5])
    # os.system('rm -r %s' % data_dir_name)
    