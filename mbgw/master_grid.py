# Master grid:
from numpy import pi

missing_val = -9999

ncols=         8664
nrows=         3384
xllc = -180
yllc = -57
cellsize = 0.04166665

AM_lims = {'topRow': 1494,
'bottomRow': 2602,
'leftCol': 2038,
'rightCol': 3401}

AF_lims = {'topRow': 1393,
'bottomRow': 2745,
'leftCol': 3879,
'rightCol': 5628}

AS_lims = {'topRow': 1018,
'bottomRow': 2512,
'leftCol': 5241,
'rightCol': 8423}

# FIXME: Replace with real Kenya limits.
KE_lims = {'topRow' : 789,
'bottomRow' : 1021,
'leftCol' : 2999,
'rightCol' : 3203}

rad_to_km = 6378.1/pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

##

