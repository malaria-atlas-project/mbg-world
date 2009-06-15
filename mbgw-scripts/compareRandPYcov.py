Rel=0
filename = "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5"

import tables as tb
import numpy as np
from rpy import *

# import R function
r.source('temptestcov.R')
temptestcovPY=r['testRcov']

hf = tb.openFile(filename)    
hr = hf.root

xd=np.array([0.1,0.2,0.3,0.4,0.5])#,0.6,0.7,0.8,0.9,1.0])
yd=np.array([0.1,0.2,0.3,0.4,0.5])#,0.6,0.7,0.8,0.9,1.0])
td=np.array([0.0,0.0,0.0,0.0,0.0])#,0.0,0.0,0.0,0.0,0.0])

Scale=hr.PyMCsamples.col("scale")[Rel]
amp=hr.PyMCsamples.col("amp")[Rel]
inc=hr.PyMCsamples.col("inc")[Rel]
ecc=hr.PyMCsamples.col("ecc")[Rel]
t_lim_corr=hr.PyMCsamples.col("t_lim_corr")[Rel]
scale_t=hr.PyMCsamples.col("scale_t")[Rel]
sin_frac=hr.PyMCsamples.col("sin_frac")[Rel]

CfromR=temptestcovPY(xd,yd,td,Scale,amp,inc,ecc,t_lim_corr,scale_t,sin_frac)

# obtain and overlay theoretical covariance function from input MCMC paramater values
C = hr.group0.C[Rel]
locArray = np.vstack((xd,yd,td)).T
CfromPY = C(locArray, locArray)
CfromPY = np.asarray(CfromPY).squeeze()


print "CfromR"
print CfromR

print "CfromPY"
print CfromPY

