from numpy.testing import *
import mbgw
from mbgw import EP
import nose,  warnings
from numpy import *
from pymc import normal_like

# ages_and_data(N_exam, N_age_samps, f_samp, correction_factor_array, age_lims)
# find_closest(x, y, N_nearest)
# observed_gp_params(combo_mesh, tracefile, trace_thin, trace_burn, N_nearest)
# simulate_data(M_pri, C_pri, N_samp, V, N_exam, N_age_samps, correction_factor_array, age_lims)
# pred_samps(pred_mesh, samp_mesh, N_exam, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, age_lims, correction_factor_array)
    
# class test_mbgw(TestCase):
class test_EP_MAP(object):
    pass
        
if __name__ == '__main__':
    # test_EP_MAP().test_low_V()
    # warnings.simplefilter('ignore',  FutureWarning)
    nose.runmodule()

