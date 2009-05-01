from tables import *
from numpy import *

#fname = 'salb1km-e_y-x+.asc'
# path='/home/pwg/mbg-world/datafiles/auxiliary_data/'
# mapView=True


def asc_to_hdf5(fname, path, mapView=False,setNaN=False):
    """
    Extracts long, lat, data from an ascii-format file.
    """
    f = file(fname,'r')
    
    # Extract metadata from asc file.
    ncols = int(f.readline()[14:])
    nrows = int(f.readline()[14:])
    xllcorner = float(f.readline()[14:])
    yllcorner = float(f.readline()[14:])
    cellsize = float(f.readline()[14:])
    NODATA_value = int(f.readline()[14:])
    
    #print('type of NODATA_value')
    #print(type(NODATA_value))

    #print('NODATA_value')
    #print(NODATA_value)
    
    # Make longitude and latitude vectors.    
    long = xllcorner + arange(ncols) * cellsize
    lat = yllcorner + arange(nrows) * cellsize
    if (mapView == True): lat = lat[::-1]
    
    # Initialize hdf5 archive.
    h5file = openFile(path, mode='w', title=fname[:-4] + ' in hdf5 format')

    # Write hdf5 archive metadata.
    h5file.root._v_attrs.asc_file = path + fname
    h5file.root._v_attrs.ncols = ncols
    h5file.root._v_attrs.nrows = nrows
    h5file.root._v_attrs.missing = NODATA_value 
    h5file.root._v_attrs.minx = long.min()
    h5file.root._v_attrs.maxx = long.max()
    h5file.root._v_attrs.miny = lat.min()
    h5file.root._v_attrs.maxy = lat.max()
    h5file.root._v_attrs.cellsize = cellsize
    if mapView==True:h5file.root._v_attrs.order = 'y-x+'
    if mapView==False:h5file.root._v_attrs.order = 'y+x+'    
    
    # Add longitude and latitude to archive, uncompressed. 
    h5file.createArray('/','long',long)
    h5file.createArray('/','lat',lat)
    
    # Add data to archive, heavily compressed, row-by-row (if a row won't fit in memory, the whole array won't fit on disk).
    h5file.createCArray('/', 'data', Float64Atom(), (nrows, ncols), filters = Filters(complevel=9, complib='zlib'))    
    data = h5file.root.data

    # optionally fill rows upwards or downwards, and optionally replace NODATA values with specifed value    
    for i in xrange(nrows):
        temp = fromstring(f.readline(), dtype=float, sep=' ')
        #print('i= '+str(i)+' shape(temp) = '+str(shape(temp))       )
        if setNaN != False:
            temp[temp == NODATA_value] = setNaN
        if mapView==False:
            data[-i-1,:] = temp
        if mapView==True:
            data[i,:] = temp 
    
    return h5file 

#asc_to_hdf5('gr001km_ken.asc', path='/home/pwg/MBGWorld/datafiles/auxiliary_data/')
#asc_to_hdf5('limbnry1km-e_ken.asc', path='/home/pwg/MBGWorld/datafiles/auxiliary_data/')
#asc_to_hdf5('salb1km-e_ken.asc', path='/home/pwg/MBGWorld/datafiles/auxiliary_data/')
#asc_to_hdf5('salblim1km-e_ken.asc', path='/home/pwg/MBGWorld/datafiles/auxiliary_data/',mapView=True)
#asc_to_hdf5('gr001km_ken.asc', path='/home/pwg/MBGWorld/datafiles/auxiliary_data/',mapView=True,setNaN=NaN)

#asc_to_hdf5('salblim1km-e.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/',mapView=True)
#asc_to_hdf5('lims1km-e.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/',mapView=True)
#asc_to_hdf5('gr001km.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/',mapView=True,setNaN=0)
#asc_to_hdf5('lim5kmbnry-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/',mapView=True)
#asc_to_hdf5('gr005km_AF.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/',mapView=True)

#asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/salb1km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/salb1km-e_y-x+.hdf5',mapView=True)
#asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/salb_st_1km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/salb_st_1km-e_y-x+.hdf5',mapView=True) 
asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/st_mask5km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/st_mask5km-e_y-x+.hdf5',mapView=True)
#asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/lims1km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/lims1km-e_y-x+.hdf5',mapView=True)
asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/urb5km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/urb5km-e_y-x+.hdf5',mapView=True)
asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/periurb5km-e_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/periurb5km-e_y-x+.hdf5',mapView=True)
#asc_to_hdf5('/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_y-x+.asc', path='/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_y-x+.hdf5',mapView=True) 




