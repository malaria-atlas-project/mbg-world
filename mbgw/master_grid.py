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

#AF_lims = {'topRow': 1393,
#'bottomRow': 2745,
#'leftCol': 3879,
#'rightCol': 5628}

#zimbabwe transect low to high
#AF_lims = {'topRow': 2410,
#'bottomRow': 2421,
#'leftCol': 5065,
#'rightCol': 5183}

#zimabwe sqaure - all low - qith some 2001 very low surveys in
#AF_lims = {'topRow': 2489,
#'bottomRow': 2498,
#'leftCol': 5079,
#'rightCol': 5088}

# another test
AF_lims = {'topRow': 2387,
'bottomRow': 2412,
'leftCol': 5147,
'rightCol': 5164}

#AF_lims = {'topRow': 1536,
#'bottomRow': 1776,
#'leftCol': 5000,
#'rightCol': 5250}

AS_lims = {'topRow': 1018,
'bottomRow': 2512,
'leftCol': 5241,
'rightCol': 8423}

# FIXME: Replace with real Kenya limits.
KE_lims = {'topRow' : 789,
'bottomRow' : 1021,
'leftCol' : 2999,
'rightCol' : 3203}


rad_to_km = 6378.1
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

##

