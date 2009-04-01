from pylab import *
from numpy.testing import *
import mbgw
from mbgw import EP
import nose,  warnings
from numpy import *
import pymc as pm
from mbgw.EP import ages_and_data
import tables as tb


spatial = True
with_urban=True
with_stukel=False
chunk=2
covariates=False

disttol=10/6378.
ttol=1/12.

N_param_vals = 5
N_per_param = 5
N_nearest = 40
N_age_samps = 1000
lo_age = 2
up_age = 11

trace_thin = 10
trace_burn = 1000

# Store an array of age-correction factors for future use
try:
    hf.close()
except:
    pass

root = mbgw.__path__[0]
hf = tb.openFile(root+'/../datafiles/good-traces/QRYPFPR101108Ken_KenyaThin_Run_11.2.2009_urb_periurb.hdf5')
tracefile = hf

rad_to_km = 6378.1/pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

n_digits = 4
iter = 10000
tol = .0001
    
class test_EP_MAP(TestCase):
# class test_EP_MAP(object):
    
    # N = random.randint(2,20)
    N = 4
    
    f_samp = random.normal(size=N)
    al = random.randint(1,12,size=N)
    al = vstack((al, al+random.randint(30, size=N))).T
    age_lims = [tuple(ali) for ali in al]
    correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 10000)
    M_pri = random.normal(size=N)
    sig_pri = random.normal(size=(N, N))
    C_pri = dot(sig_pri.T, sig_pri)
        
        
    def test_ages_and_data(self):
        "Makes sure that the ages drawn by ages_and_data match the age distributions."
        N_exam = ones(self.N) * 10000
        A, pos, age_distribution = ages_and_data(N_exam, self.f_samp, self.correction_factor_array, self.age_lims)
        empirical_age_distributions = []
        
        for j in xrange(self.N):
            A_ind = A[j]-self.age_lims[j][0]
            empirical_age_distributions.append(array([sum(A_ind==i) for i in xrange(len(age_distribution[j]))]) / 10000.)
            bin_sds = sqrt(age_distribution[j] * (1-age_distribution[j]) / 10000.)
            assert(all(abs(age_distribution[j]-empirical_age_distributions[j]) < 4*bin_sds))
    
    def test_likelihoods(self):
        "Makes sure that the sparsified likelihood function procedure is reasonable."
        import time
        
        nug = random.normal()**2 * .3
        fs = array([0])
        A, pos, age_distribution = ages_and_data([300], fs, self.correction_factor_array, [self.age_lims[0]])

        t1 = time.time()
        lp_small = mbgw.correction_factors.known_age_corr_likelihoods_f(pos, A, self.correction_factor_array, linspace(-5,5,500), nug, 's')[0]
        # print 'Small', time.time() - t1
        t2 = time.time()
        lp_large = mbgw.correction_factors.known_age_corr_likelihoods_f(pos, A, self.correction_factor_array, linspace(-5,5,500), nug, 'l')[0]
        # print 'Large', time.time() - t2
        
        x = linspace(-5, 5, 100)
        small = lp_small(x)
        large = lp_large(x)
        
        assert_almost_equal(small, large, 8)
    
    
    def test_fit(self):
        "Compares EP results with MCMC results with a real age-corrected likelihood."
        nug = random.normal()**2 * .3
        N_exam = ones(self.N) * 1000
        lps, pos = EP.EP_MAP.simulate_data(self.M_pri, self.C_pri, self.N, nug, N_exam, self.correction_factor_array.shape[1], self.correction_factor_array, self.age_lims)

        # Do EP algorithm
        E = EP.EP(self.M_pri, self.C_pri, lps, nug*ones(self.N))
        E.fit(iter, tol=tol)
        # from IPython.Debugger import Pdb
        # Pdb(color_scheme='Linux').set_trace()   
        
        x = pm.MvNormalCov('x', self.M_pri, self.C_pri)
        eps = pm.Normal('eps', x, 1./nug)
        @pm.potential
        def y(eps=eps, lps=lps):
            return sum([lps[i](eps[i]) for i in xrange(len(eps))])
        
        M = pm.MCMC([x,eps,y])
        M.use_step_method(pm.AdaptiveMetropolis, [eps, x])
        M.isample(100000)
        
        post_V = var(M.trace('x')[20000:], axis=0)
        post_M = mean(M.trace('x')[20000:], axis=0)
        
        pm.Matplot.plot(M)
        
        assert_almost_equal(post_M, E.M, 1)        
        assert_almost_equal(post_V, diag(E.C), 1)
                
    def test_pred_samps(self):
        "A dry run in Kenya with only one sample point. This test should not work with N>1."
        
        N = 1
        
        lat_pred = np.atleast_1d(pm.runiform(-5., 5., size=N) * deg_to_rad)
        # lat_pred = array([8.89, 9.5, 1.17, 1.39])
        lon_pred = np.atleast_1d(pm.runiform(33., 40., size=N) * deg_to_rad)
        # lon_pred = array([-1.54, .08, 39.44, 38.12])
        t_pred = np.atleast_1d(array([2007]*N)-2009)

        pred_mesh = vstack((lon_pred, lat_pred, t_pred)).T
        age_lims = [(lo_age, up_age)]*len(lon_pred)

        N_exam = ones(len(lat_pred))*1000
                
        input_pts = [{'lon': lon_pred[i], 'lat': lat_pred[i], 'month': 1, 'year': 2009, 'lo_age': 2, 'up_age': 10, 'n': N_exam[i]}\
                        for i in range(len(lat_pred))]
        output_pts =  [{'lon': lon_pred[i], 'lat': lat_pred[i], 'year': 2009, 'month': 1, 'lo_age': 2, 'up_age': 10, 'nmonths': 2} for i in range(len(lat_pred))]

        correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 1000)

        ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors =\
            mbgw.EP.pred_samps(pred_mesh*deg_to_rad, pred_mesh*deg_to_rad, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array, debug=True)


if __name__ == '__main__':
    # tester = test_EP_MAP()
    # tester.test_pred_samps()
    # tester.test_likelihoods()
    # tester.test_fit()
    # tester.check_ages_and_data()
    # test_EP_MAP().test_low_V()
    # warnings.simplefilter('ignore',  FutureWarning)
    nose.runmodule()


