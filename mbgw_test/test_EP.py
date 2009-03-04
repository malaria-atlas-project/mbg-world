import mbgw
from mbgw.EP import *
from mbgw.correction_factors import *
import numpy as np
import tables as tb

def entropy(samps):
    raise NotImplemented

def variance(samps):
    return -np.var(samps)

spatial = True
with_urban=True
with_stukel=False
chunk=2
covariates=False

disttol=10/6378.
ttol=1/12.

N_param_vals = 200
N_per_param = 500
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

# hf = tb.openFile('../Kenya_25_realizations.hdf5')
hf = tb.openFile('../datafiles/good-traces/QRYPFPR101108Ken_KenyaThin_Run_11.2.2009_urb_periurb.hdf5')
tracefile = hf

rad_to_km = 6378.1/np.pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./np.pi
deg_to_rad = 1./rad_to_deg

# # From EP_MAP
# lat_pred = np.array([8.89, 9.5, 1.17, 1.39])
# lon_pred = np.array([-1.54, .08, 39.44, 38.12])
# t_pred = np.array([2007]*4)-2009
# 
# x_pred = np.vstack((lon_pred, lat_pred))*np.pi/180.
# pred_mesh = np.vstack((x_pred, t_pred)).T
# 
# N_exam = np.array([1000,3000,2000,4000])
# utility=variance
# 
# # pos = post_samps(pred_mesh)
# io, ii, ms, cs, lm, lv, p = pred_samps(pred_mesh, pred_mesh, N_exam)

# From frontend_interface

lat_pred = np.array([8.89, 9.5, 1.17, 1.39])
lon_pred = np.array([-1.54, .08, 39.44, 38.12])
t_pred = np.array([2007]*4)-2009

pred_mesh = np.vstack((lon_pred, lat_pred, t_pred)).T

# N_exam = np.array([1000,3000,2000,4000])    
N_exam = np.array([1,1,1,1]) 
input_pts = [{'lon': lon_pred[i], 'lat': lat_pred[i], 'month': 1, 'year': 2009, 'lo_age': 2, 'up_age': 10, 'n': N_exam[i]}\
                for i in range(4)]
output_pts =  [{'lon': lon_pred[i], 'lat': lat_pred[i], 'year': 2009, 'month': 1, 'lo_age': 2, 'up_age': 10, 'nmonths': 2} for i in range(4)]

client_side = frontend(update_posterior)
client_cleanup = frontend(scratch_cleanup)
correction_factor_array = known_age_corr_factors(np.arange(lo_age, up_age), N_age_samps)
# path, llc, urc, output_info = client_side(input_pts, output_pts, 2., 2009)
# ind_outer, ind_inner, Ms, Cs, Vs, likelihood_means, likelihood_variances, model_posteriors =\
#     pred_samps(pred_mesh*deg_to_rad, pred_mesh*deg_to_rad, N_exam, tracefile, trace_thin, N_param_vals, N_per_param, N_nearest, correction_factor_array)

update_posterior(input_pts, output_pts, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, utilities=[np.std])
# update_posterior(input_pts, output_pts, utilities=[np.std, np.var])