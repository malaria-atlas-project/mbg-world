# Master grid:
from numpy import pi

missing_val = -9999

ncols=         8664
nrows=         3384
xllc = -180
yllc = -57
cellsize = 0.04166665

AM_lims = {'topRow': 1890,
'bottomRow': 2140,
'leftCol': 2038,
'rightCol': 3401}

AF_lims = {'topRow': 1393,
'bottomRow': 2745,
'leftCol': 3879,
'rightCol': 5628}

# Two test squares (S1 fits indside S2) in S malawi - includes a high and low focus and an urban area
S1_lims = {'topRow': 2387,
'bottomRow': 2412,
'leftCol': 5147,
'rightCol': 5164}

S2_lims = {'topRow': 2373,
'bottomRow': 2426,
'leftCol': 5120,
'rightCol': 5197}

# TEMP - kenya lims as Af for testing
#AF_lims = {'topRow' : 1901,
#'bottomRow' : 2132,
#'leftCol' : 5135,
#'rightCol' : 5325}

AS_lims = {'topRow': 1018,
'bottomRow': 2512,
'leftCol': 5241,
'rightCol': 8423}

# Kenya limits.
KE_lims = {'topRow' : 1901,
'bottomRow' : 2132,
'leftCol' : 5135,
'rightCol' : 5325}

# very small test case
VS_lims = {'topRow' : 1901,
'bottomRow' : 1908,
'leftCol' : 5135,
'rightCol' : 5141}


rad_to_km = 6378.1
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

##

