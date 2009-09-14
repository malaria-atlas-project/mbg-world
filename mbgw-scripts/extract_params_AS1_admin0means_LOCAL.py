from numpy import array

# set start and end month slice (python array index of months in block to aggregate over - 0 is earliest month in block)
MonthsSlice = slice(0,12,None)
 
# standard path to utility function directory (used to source generic R functions)
utilFolder = '/home/pwg/map_utils/map_utils/'
 
# set path to file containg keys for amazon S3
#keyPath = '/root/s3code.txt' 
 
# main input hdf5  file of simulated realisations of f
realisations_path = '/home/pwg/mbg-world/realizations/AS1/qrypfpr230708_asia_run_1.9.2008_as1_try2/realizations_mem_100000000_QRYPFPR230708_Asia_Run_1.9.2008_iterations_FILESTARTREL_FILEENDREL.hdf5'

# hdf5 file containng traces for burden function
burdentrace_path='/home/pwg/mbg-world/datafiles/burdentraces/CSE_Asia_and_Americas_scale_0.6_model_exp.hdf5'

# location for export of raw extractions (as they come off each distributed instance)
exportPathDistributed_country = '/home/pwg/mbg-world/extractionOutput/distributedoutput_country_as1/admin0/'
exportPathDistributed_perpixel = '/home/pwg/mbg-world/extractionOutput/distributedoutput_perpixel_as1/'
 
# location for export of combined extractions (after distributed files joined by extract_combineDistribExtractions.py)
exportPathCombined_country = '/home/pwg/mbg-world/extractionOutput/combinedoutput_country_as1/admin0/'
exportPathCombined_perpixel = '/home/pwg/mbg-world/extractionOutput/combinedoutput_perpixel_as1/'
 
# input salb rasters of unique spatial IDs (trimmed to limits)
#salblim1km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/NULL.hdf5"
salblim5km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/salblim5km-e_y-x+_AS1.hdf5"

# input full salb rasters
salb5km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/salb5km-e2_y-x+_AS1.hdf5"
 
# input 1km and 5km raster of population per cell
#grump1km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/gr071km_y-x+_AS1.hdf5" 
grump5km_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/gr075km_y-x+_AS1.hdf5"

# input 1k mand 5km rater showing area of each pixel
#pixarea1km_path = "/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/pixarea1km_y-x+_AS1.hdf5"
pixarea5km_path = "/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/pixarea5km_y-x+_AS1.hdf5"
 
# location of 5km hdf5  limits mask
lim5kmbnry_path="/home/pwg/mbg-world/datafiles/auxiliary_data/GridsForCS/st_mask5km-e_y-x+_AS1.hdf5"
 
# files containing list of unique salb IDs in input raster and pixels per ID : generated as ouptut from examineSalb
examineSalbFolder='examinesalb'
uniqueSalb_path='/home/pwg/mbg-world/extractionOutput/examineSalb/uniquesalb_ad0_as1.txt'
pixelN_path='/home/pwg/mbg-world/extractionOutput/examineSalb/pixelN_ad0_as1.txt'
salbpop_path = '/home/pwg/mbg-world/extractionOutput/examineSalb/salbpop_ad0_as1.txt'
salbarea_path = '/home/pwg/mbg-world/extractionOutput/examineSalb/salbarea_ad0_as1.txt'

# age limits for back-transform
a_lo=2
a_hi=10

# what do we want to extract?

PERPIXEL = False  # do we want to extract maps
PERCOUNTRY = True # do we want to aggregate over spatial units like countries or admin1s?

do_AREALMEANPR = False   # area-weighted mean PR (applies to PERCOUNTRY extractions only)
do_POPMEANPR = True      # population-weighted mean PR (applies to PERCOUNTRY extractions only)
do_AREALMEANRo = False   # area-weighted mean Ro (applies to PERCOUNTRY extractions only)
do_POPMEANRo = False     # population-weighted mean Ro(applies to PERCOUNTRY extractions only)
do_PAR = False           # population at risk, will be summarised by breaksDict classes (applies to both PERCOUNTRY and PERPIXEL extractions)
do_BURDEN = False        # burden, will be summarised by breaksDict classes  (applies to both PERCOUNTRY and PERPIXEL extractions)

do_BurdenMap=False       # map of clinical cases per pixel
do_PRMap=True            # map of PfPR per pixel
do_RoMap=True            # map of PfRc per pixel

# class definition dictionaries
#breaks_MBGW={"BREAKS":[0.,0.05,0.40,1.1],"BREAKNAMES":["lte05","gt05lte40","gt40lte100"],"NAME":"MBGW"}
#breaks_MVI={"BREAKS":[0.,0.05,0.30,0.40,0.60,1.1],"BREAKNAMES":["lte05","gt05lte30","gt30lte40","gt40lte60","gt60lte100"],"NAME":"MVI"}
#breaks_NOORKENYA={"BREAKS":[0.,0.01,0.05,0.4,1.1],"BREAKNAMES":["lte1","gt1lte0","gt5lte40","gt40lte100"],"NAME":"NOORKENYA"}
#breaks_LYSENKO={"BREAKS":[0.,0.1,0.5,0.75,1.1],"BREAKNAMES":["lte10","gt10lte50","gt50lte75","gt75lte100"],"NAME":"LYSENKO"}
breaks_one={"BREAKS":[0.,1.1],"BREAKNAMES":["all"],"NAME":"ONE"}

#breaksDict={"MBGW":breaks_MBGW,"MVI":breaks_MVI,"LYSENKO":breaks_LYSENKO,"ONE":breaks_one}
#breaksDict={"MBGW":breaks_MBGW}
breaksDict={"ONE":breaks_one}

# ratio of high to low resolution
HiResLowResRatio_PERPIXEL=1
HiResLowResRatio_PERCOUNTRY=1

# how many rows of the 5km PR surface will we deal with at a time in extractSummaries_country()
rowsPerChunk=1
 
# parameter determining number of rows of the coarser grid (5km) to process at a time. Default is 1. (currently unsupported)
rowsInslice5km=1 
 
# what summary stats do we want to generate for meanPR, PAR etc in each country? If none, then set to None
summaryStats={'mean':'mean','SD':'SD','quantiles':array([0,0.025,0.25,0.5,0.75,0.975,1])}
#summaryStats=None
 
 
 
 
 
 
