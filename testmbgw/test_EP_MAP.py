from numpy.testing import *
import mbgw
from mbgw import EP
import nose,  warnings
from numpy import *
from pymc import normal_like
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

lat_pred = array([8.89, 9.5, 1.17, 1.39])
lon_pred = array([-1.54, .08, 39.44, 38.12])
t_pred = array([2007]*4)-2009

pred_mesh = vstack((lon_pred, lat_pred, t_pred)).T
age_lims = [(lo_age, up_age)]*len(lon_pred)


    
# class test_mbgw(TestCase):
class test_EP_MAP(object):
    
    N = random.randint(2,20)
    f_samp = random.normal(size=N)
    al = random.randint(0,5,size=N)
    al = vstack((al, al+random.randint(30, size=N))).T
    age_lims = [tuple(ali) for ali in al]
    correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 1000)
        
    def test_ages_and_data(self):
        N_exam = ones(self.N) * 10000
        A, pos, age_distribution = ages_and_data(N_exam, self.f_samp, self.correction_factor_array, self.age_lims)
        empirical_age_distributions = []
        for j in xrange(self.N):
            A_ind = A[j]-self.age_lims[j][0]
            empirical_age_distributions.append(array([sum(A_ind==i) for i in xrange(len(age_distribution[j]))]) / 10000.)
            bin_sds = sqrt(age_distribution[j] * (1-age_distribution[j]) / 10000.)
            assert(all(abs(age_distribution[j]-empirical_age_distributions[j]) < 4*bin_sds))
    
    def test_small_sample(self):

        N_exam = array([1,1,1,1])*1
        input_pts = [{'lon': lon_pred[i], 'lat': lat_pred[i], 'month': 1, 'year': 2009, 'lo_age': 2, 'up_age': 10, 'n': N_exam[i]}\
                        for i in range(4)]
        output_pts =  [{'lon': lon_pred[i], 'lat': lat_pred[i], 'year': 2009, 'month': 1, 'lo_age': 2, 'up_age': 10, 'nmonths': 2} for i in range(4)]

        correction_factor_array = mbgw.correction_factors.known_age_corr_factors(arange(0,27), 1000)

        ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors =\
            mbgw.EP.pred_samps(pred_mesh*deg_to_rad, pred_mesh*deg_to_rad, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array, debug=False)



if __name__ == '__main__':
    tester = test_EP_MAP()
    # tester.test_small_sample()
    # tester.check_ages_and_data()
    # test_EP_MAP().test_low_V()
    # warnings.simplefilter('ignore',  FutureWarning)
    # nose.runmodule()


# FIXME: Profiler results are pretty stupid:    
#    ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#   3035518   50.997    0.000  205.972    0.000 correction_factors.py:190(this_fun)
#   9106690   34.756    0.000   38.387    0.000 shape_base.py:242(atleast_1d)
#   6071036   31.072    0.000   61.004    0.000 warnings.py:24(warn)
#   6072224   27.415    0.000   27.415    0.000 {method 'any' of 'numpy.ndarray' objects}
#   3035518   21.052    0.000  268.549    0.000 EP.py:82(<lambda>)
#   3035518   16.332    0.000   77.336    0.000 {scipy.interpolate._fitpack._spl_}
#   3035518   16.196    0.000  104.158    0.000 fitpack.py:434(splev)
#   3035543   11.089    0.000   11.089    0.000 distributions.py:1890(normal_like)
#   6071036    9.440    0.000   15.041    0.000 warnings.py:64(warn_explicit)
#   3035518    9.274    0.000  230.591    0.000 EP.py:81(<lambda>)
#   3035518    5.817    0.000   16.906    0.000 EP.py:80(<lambda>)
#   6072339    5.263    0.000   32.680    0.000 fromnumeric.py:1271(any)
#  12182762    5.148    0.000    5.148    0.000 {isinstance}
#       300    4.366    0.015  275.844    0.919 EP.py:69(compute_expectation)
#  12142141    4.213    0.000    4.213    0.000 {method 'get' of 'dict' objects}
#   6071038    3.392    0.000    3.392    0.000 {method 'lower' of 'str' objects}
#   3035553    3.103    0.000    3.103    0.000 {method 'reshape' of 'numpy.ndarray' objects}
#   6072854    2.928    0.000    2.928    0.000 {issubclass}
#   9108689    2.398    0.000    2.398    0.000 {method 'append' of 'list' objects}
#   6071051    1.994    0.000    1.994    0.000 {method 'endswith' of 'str' objects}
#   6071036    1.879    0.000    1.879    0.000 {sys._getframe}
#         5    1.714    0.343    3.081    0.616 correction_factors.py:132(known_age_corr_likelihoods_f)
# 12152260/12152259    1.677    0.000    1.677    0.000 {len}
#       398    1.419    0.004    1.427    0.004 st_cov_fun.py:26(my_st)
#   3035946    1.337    0.000    1.337    0.000 fromnumeric.py:1006(shape)
#        20    1.168    0.058    1.280    0.064 shape_base.py:11(apply_along_axis)
#       900    1.004    0.001    5.060    0.006 optimize.py:100(fmin)
#   6071036    0.957    0.000    0.957    0.000 {method 'setdefault' of 'dict' objects}
#   3035618    0.849    0.000    0.849    0.000 {method 'squeeze' of 'numpy.ndarray' objects}
#      1800    0.567    0.000    0.606    0.000 quadrature.py:149(_basic_simps)
#       308    0.369    0.001    1.850    0.006 FullRankCovariance.py:64(cholesky)
#       924    0.252    0.000    0.260    0.000 {method '_fillCol' of 'tableExtension.Row' objects}
#       600    0.225    0.000    0.225    0.000 EP.py:96(<lambda>)
#     26924    0.165    0.000    0.165    0.000 EP.py:31(<lambda>)