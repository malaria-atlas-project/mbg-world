# import libraries
from mbgw import master_grid
import tables as tb
from numpy import *

# set some parameters
resRatio=5  # ratio between resolution of grids defines in master_grid, and the ones being subsetted here
datafolder = "/home/pwg/mbg-world/datafiles/auxiliary_data/"



################################################################################################################################################
# using row/col numbers of a subset specified for 5km grids (and stored in master_grid), extract the same subset for corresponding 1km grids

def subset1kmgrids (region,lims,gridname):
    # build input path for the global grid
    input_path = datafolder+gridname+".hdf5"

    # build ouput paths for this region and grid 
    output_path = datafolder+gridname+"_"+region+".hdf5"

    # open link to full sized 1km file
    fullHDF5 = tb.openFile(input_path, mode = "r")    

    # get 1km row/cols for this subset from the pre-defined object in mbgw 'master_grid' (NB - the values therein were entered manually after running the R script DefineSubsetCoordinates.R)
    br=(lims['bottomRow'])*resRatio
    tr=(lims['topRow']-1)*resRatio
    lc=(lims['leftCol']-1)*resRatio
    rc=(lims['rightCol'])*resRatio

    # define subsetted lat and long vector
    long = fullHDF5.root.long[lc:rc:1]
    lat = fullHDF5.root.lat[tr:br:1]

    nrows = len(lat)
    ncols = len(long)

    #print(ncols)
    #print(nrows)

    # define subsetted grid
    gridsubset = fullHDF5.root.data[tr:br:1,lc:rc:1]
    #print(mean(gridsubset))
    #print(shape(gridsubset))
    #print (type(gridsubset))
    ##
    ## now build new hdf5 file for subsetted grid
    ##

    # Initialize hdf5 archive
    outHDF5 = tb.openFile(output_path, mode='w', title=gridname+'_'+region)

    # build grid metadata 
    outHDF5.root._v_attrs.asc_file = 'subsetted from hdf5 file: '+input_path
    outHDF5.root._v_attrs.ncols = ncols
    outHDF5.root._v_attrs.nrows = nrows
    outHDF5.root._v_attrs.missing = fullHDF5.root._v_attrs.missing
    outHDF5.root._v_attrs.minx = long.min()
    outHDF5.root._v_attrs.maxx = long.max()
    outHDF5.root._v_attrs.miny = lat.min()
    outHDF5.root._v_attrs.maxy = lat.max()

    # Add longitude and latitude to archive, uncompressed. 
    outHDF5.createArray('/','long',long)
    outHDF5.createArray('/','lat',lat)

    # Add data to archive, heavily compressed, in a chunk array row-by-row (if a row won't fit in memory, the whole array won't fit on disk).
    outHDF5.createCArray('/', 'data', tb.Float64Atom(), (nrows, ncols), filters = tb.Filters(complevel=9, complib='zlib'),chunkshape = (1,ncols))    

    for i in xrange(nrows):
        outHDF5.root.data[i,:] = gridsubset[i,:]

    # close the files
    fullHDF5.close()
    outHDF5.close()
################################################################################################################################################

#subset1kmgrids(region="AM",lims = master_grid.AM_lims,gridname = "salblim1km-e")
#subset1kmgrids(region="AM",lims = master_grid.AM_lims,gridname = "gr001km")
subset1kmgrids(region="AF",lims = master_grid.AF_lims,gridname = "salblim1km-e")
subset1kmgrids(region="AF",lims = master_grid.AF_lims,gridname = "salb1km-e")
subset1kmgrids(region="AF",lims = master_grid.AF_lims,gridname = "lims1km-e")
subset1kmgrids(region="AF",lims = master_grid.AF_lims,gridname = "gr001km")
#subset1kmgrids(region="AS",lims = master_grid.AS_lims,gridname = "salblim1km-e")
#subset1kmgrids(region="AS",lims = master_grid.AS_lims,gridname = "gr001km")
