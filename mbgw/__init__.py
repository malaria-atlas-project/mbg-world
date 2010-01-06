# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

#from testmbgw import test
try:
    from testmbgw import test
except ImportError:
    pass

from model import *
import pymc as pm
import numpy as np
import os
from copy import copy
from correction_factors import age_corr_likelihoods, age_corr_factors, two_ten_factors
from scipy import interpolate as interp
from st_cov_fun import *
import time
import auxiliary_data
# import MAPData
import gc
# from get_covariates import extract_environment_to_hdf5
from tables import ObjectAtom
from generic_mbg import *


f_labels = ['eps_p_f']
fs_have_nugget = {'eps_p_f': True}
nugget_labels = {'eps_p_f': 'V'}
M_labels = {'eps_p_f': 'M'}
C_labels = {'eps_p_f': 'C'}
x_labels = {'eps_p_f': 'data_mesh'}
diags_safe = {'eps_p_f': True}

# Extra stuff for predictive ops.
n_facs = 1000
# postproc = invlogit
def map_postproc(eps_p_f, two_ten_facs=two_ten_factors(n_facs)):
    return invlogit(eps_p_f) * two_ten_facs[np.random.randint(n_facs)]
    
metadata_keys = ['ti','fi','ui','with_stukel','chunk','disttol','ttol']

def mcmc_init(M):
    M.use_step_method(FieldStepper, M.f, M.V, M.C_eval, M.M_eval, M.logp_mesh, M.eps_p_f, M.ti)

non_cov_columns = {'lo_age': 'int', 'up_age': 'int', 'pos': 'float', 'neg': 'float'}
