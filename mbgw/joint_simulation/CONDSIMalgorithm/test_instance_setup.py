print 'Starting test script.py...'

# import libraries
import numpy as np
import pymc as pm
import mbgw
import time
from mbgw.joint_simulation import *
import tables as tb
import mbgw.master_grid as mg
import os,sys
from map_utils import S3
from CONDSIM_params import *

cd ..


import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
from fast_krige import *
import st_cov_fun
from parse_and_check import *
import time
import copy as cp
from mbgw.master_grid import *
from mbgw import auxiliary_data
import os
curpath = os.getcwd()
import mbgw
mbgw_root = __root__ = mbgw.__path__[0]
r_path = mbgw_root+'/joint_simulation/CONDSIMalgorithm'
os.chdir(r_path)
from examineRealization import *
from rpy import r
r.source("CONDSIMpreloop.R")
r.source("CONDSIMmonthloop.R")
r.source('MVRNORM.R')
mvrnormPY=r['MVRNORM']

print 'DONE test script.py'
