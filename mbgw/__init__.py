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


# import agepr
# import st_cov_fun
# import google_earth
# import os
# import povray
# import joint_simulation
# import final_model
# import correction_factors
# import EP

# err_modules = []
# for module in ['agepr', 'st_cov_fun', 'google_earth', 'os', 'povray', 'joint_simulation', 'final_model', 'correction_factors', 'EP', 'master_grid', 'hdf5_utils']:
#     try:
#         exec('import %s'%module)
#     except:
#         err_modules.append(module)
#         
# if len(err_modules) > 0:
#     print 'Failed to import the following modules: %s'%(', '.join(err_modules))

# import get_covariates
# from writeout_from_hdf5 import *