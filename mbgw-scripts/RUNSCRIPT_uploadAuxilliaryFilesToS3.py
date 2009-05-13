from boto_PYlib import *

S=S3()

uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/gr071km_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)
uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/lims1km-e_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)
uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/gr075km_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)
uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/salb1km-e2_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)
uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/salblim1km-e_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)
uploadFileToBucket(auxiliary_data,"/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/st_mask5km-e_y-x+_AF.hdf5",overwriteContent=False,makeBucket=False,VERBOSE=True)


