# Master grid:
from numpy import pi

missing_val = -9999

ncols=         8664
nrows=         3384
xllc = -180
yllc = -57
cellsize = 0.04166665

AM_lims = {'topRow': 1525,
'bottomRow': 2483,
'leftCol': 2125,
'rightCol': 3296}

AMS1_lims = {'topRow': 2299,
'bottomRow': 2325,
'leftCol': 2784,
'rightCol': 2801}

AF_lims = {'topRow': 1412,
'bottomRow': 2726,
'leftCol': 3886,
'rightCol': 5603}

# Two test squares (S1 fits indside S2) in S malawi - includes a high and low focus and an urban area
S1_lims = {'topRow': 2387,
'bottomRow': 2412,
'leftCol': 5147,
'rightCol': 5164}

S2_lims = {'topRow': 2373,
'bottomRow': 2426,
'leftCol': 5120,
'rightCol': 5197}

#AS_lims = {'topRow': 1018,
#'bottomRow': 2512,
#'leftCol': 5241,
#'rightCol': 8423}


AS1_lims = {'topRow': 1109,
'bottomRow': 1910,
'leftCol': 5575,
'rightCol': 7366}

AS2_lims = {'topRow': 1838,
'bottomRow': 2522,
'leftCol': 6592,
'rightCol': 8419}

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

