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

hf = tb.openFile('../datafiles/good-traces/QRYPFPR101108Ken_KenyaThin_Run_11.2.2009_urb_periurb.hdf5')
tracefile = hf

rad_to_km = 6378.1/pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

n_digits = 4
iter = 10000
tol = .0001

# TODO: With large N, correction_factors.known_age_corr_likelihoods_f is by far the bottleneck. You can optimize it there.
    
# class test_mbgw(TestCase):
class test_EP_MAP(object):
    
    # N = random.randint(2,20)
    N = 4
    
    f_samp = random.normal(size=N)
    al = random.randint(1,12,size=N)
    al = vstack((al, al+random.randint(30, size=N))).T
    age_lims = [tuple(ali) for ali in al]
    correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 1000)
    M_pri = random.normal(size=N)
    sig_pri = random.normal(size=(N, N))
    C_pri = dot(sig_pri.T, sig_pri)
        
    def test_ages_and_data(self):
        N_exam = ones(self.N) * 10000
        A, pos, age_distribution = ages_and_data(N_exam, self.f_samp, self.correction_factor_array, self.age_lims)
        empirical_age_distributions = []
        
        for j in xrange(self.N):
            A_ind = A[j]-self.age_lims[j][0]
            empirical_age_distributions.append(array([sum(A_ind==i) for i in xrange(len(age_distribution[j]))]) / 10000.)
            bin_sds = sqrt(age_distribution[j] * (1-age_distribution[j]) / 10000.)
            assert(all(abs(age_distribution[j]-empirical_age_distributions[j]) < 4*bin_sds))
    
    def test_fit(self):
        nug = random.normal()**2 * .3
        N_exam = ones(self.N) * 100
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
        M.isample(20000)
        
        post_V = var(M.trace('x')[5000:], axis=0)
        post_M = mean(M.trace('x')[5000:], axis=0)
        
        pm.Matplot.plot(M)
        
        assert_almost_equal(post_M, E.M, 1)        
        assert_almost_equal(post_V, diag(E.C), 1)
        
        
    def test_pred_samps(self):
        lat_pred = pm.runiform(-5., 5., size=4)
        # lat_pred = array([8.89, 9.5, 1.17, 1.39])
        lon_pred = pm.runform(35., 40., size=4)
        # lon_pred = array([-1.54, .08, 39.44, 38.12])
        t_pred = array([2007]*4)-2009

        pred_mesh = vstack((lon_pred, lat_pred, t_pred)).T
        age_lims = [(lo_age, up_age)]*len(lon_pred)

        N_exam = ones(len(lat_pred))*100
                
        input_pts = [{'lon': lon_pred[i], 'lat': lat_pred[i], 'month': 1, 'year': 2009, 'lo_age': 2, 'up_age': 10, 'n': N_exam[i]}\
                        for i in range(len(lat_pred))]
        output_pts =  [{'lon': lon_pred[i], 'lat': lat_pred[i], 'year': 2009, 'month': 1, 'lo_age': 2, 'up_age': 10, 'nmonths': 2} for i in range(len(lat_pred))]

        correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 1000)

        ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors =\
            mbgw.EP.pred_samps(pred_mesh*deg_to_rad, pred_mesh*deg_to_rad, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array, debug=True)
        # from IPython.Debugger import Pdb
        # Pdb(color_scheme='Linux').set_trace()   



if __name__ == '__main__':
    # tester = test_EP_MAP()
    # tester.test_fit()
    # tester.check_ages_and_data()
    # test_EP_MAP().test_low_V()
    # warnings.simplefilter('ignore',  FutureWarning)
    nose.runmodule()


