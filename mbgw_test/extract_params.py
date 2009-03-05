
# main input hdf5  file of simulated realisations of f
filename = '/home/pwg/Desktop/test_sim_Af.hdf5'

# location for export of raw extractions (as they come off each distributed instance)
exportPath = '/home/pwg/mbg-world/extraction/DistributedOutput/'

# location for export of combined extractions (after distributed files joined by extract_combineDistribExtractions.py)
exportPathCombined = '/home/pwg/mbg-world/extraction/CombinedOutput/'

# input 1km salb raster of unique spatial IDs
salblim1km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/salblim1km-e_ken.hdf5"

# input 1km raster of population per cell
gr001km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/gr001km_ken.hdf5"

# files containing list of unique salb IDs in input raster and pixels per ID : generated as ouptut from FUNexamineSalb
uniqueSalb_path='/home/pwg/mbg-world/extraction/uniqueSalb.txt'
pixelN_path='/home/pwg/mbg-world/extraction/pixelN.txt'

# class definition dictionaries
breaks_MBGW={"BREAKS":[0.,0.05,0.40,1.1],"BREAKNAMES":["lte05","gt05lte40","gt40lte100"],"NAME":"MBGW"}
breaks_MVI={"BREAKS":[0.,0.05,0.30,0.40,0.60,1.1],"BREAKNAMES":["lte05","gt05lte30","gt30lte40","gt40lte60","gt60lte100"],"NAME":"MVI"}
#breaks_LYSENKO={"BREAKS":[0.,0.1,0.5,0.75,1.1],"BREAKNAMES":["lte10","gt10lte500","gt50lte750","gt75lte1000"],"NAME":"LYSENKO"}
#breaks_one={"BREAKS":[0.,1.1],"BREAKNAMES":["all"],"NAME":"ONE"}
#breaksDict={"MBGW":breaks_MBGW,"MVI":breaks_MVI,"LYSENKO":breaks_LYSENKO,"ONE":breaks_one}
breaksDict={"MBGW":breaks_MBGW, "MVI":breaks_MVI}

# ratio of high to low resolution
HiResLowResRatio=5

# parameter determining number of rows of the coarser grid (5km) to process at a time. Default is 1. (currently unsupported)
rowsInslice5km=5 
