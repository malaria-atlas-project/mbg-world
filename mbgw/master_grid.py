# Master grid:
from numpy import pi

missing_val = -99999

ncols=         6281
nrows=         1605
xllc = -91.45003542
yllc = -29.23335764
cellsize = 0.04166665

AM_lims = {'topRow': 421,
'bottomRow': 1555,
'leftCol': 1,
'rightCol': 1307}

AF_lims = {'topRow': 1,
'bottomRow': 1605,
'leftCol': 1623,
'rightCol': 3515}

AS_lims = {'topRow': 1,
'bottomRow': 1391,
'leftCol': 3459,
'rightCol': 6281}

# FIXME: Replace with real Kenya limits.
KE_lims = {'topRow' : 789,
'bottomRow' : 1021,
'leftCol' : 2999,
'rightCol' : 3203}

rad_to_km = 6378.1/pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

